#!/bin/Rscript
#library(rjson)
#library(tidyverse)


## Load json files into R
#dr_res_table <- rjson::fromJSON(file="/home/dfgmrtc/Practice/run_folder/n_animal_dr_res.json")
#dr_res_int_table <- rjson::fromJSON(file="/home/dfgmrtc/Practice/run_folder/n_animal_dr_res_int.json")
#lineage_profile_table <- rjson::fromJSON(file="/home/dfgmrtc/Practice/run_folder/n_animal_lineage.json")

########## PROCESS RESISTANCE DATA ###############


### Extract candidates from dr_res_df and append to empty vectors

fill_res_vectors <- function(dr_res_df){
    ## Instantiate empty vectors for all variable/column names in empty dataframe
    sample_name_vec <- character()
    drug_name_vec <- character()
    conf_grading_vec <- character()
    gene_name_vec <- character()
    variant_vec <- character()
    variant_type_vec <- character()
    variant_pos_vec <- character()
    allele_change_vec <- character()
    for (sample_name in names(dr_res_df)){
        for (drug_name in names(dr_res_df[[sample_name]])){
            for (conf_grading in names(dr_res_df[[sample_name]][[drug_name]])){
                for (gene_name in names(dr_res_df[[sample_name]][[drug_name]][[conf_grading]])){
                    for (variant in names(dr_res_df[[sample_name]][[drug_name]][[conf_grading]][[gene_name]])){
                        sample_name_vec <- c(sample_name_vec, sample_name)
                        drug_name_vec <- c(drug_name_vec, drug_name)
                        conf_grading_vec <- c(conf_grading_vec, conf_grading)
                        gene_name_vec <- c(gene_name_vec, gene_name)
                        variant_vec <- c(variant_vec, variant)
                        variant_type_vec <- c(variant_type_vec, dr_res_df[[sample_name]][[drug_name]][[conf_grading]][[gene_name]][[variant]][1])
                        variant_pos_vec <- c(variant_pos_vec, dr_res_df[[sample_name]][[drug_name]][[conf_grading]][[gene_name]][[variant]][2])
                        allele_change_vec <- c(allele_change_vec, dr_res_df[[sample_name]][[drug_name]][[conf_grading]][[gene_name]][[variant]][3])
                        }
                    } 
                }
            }
        }
    ### Construct a vector of lengths each filled vector elements
    vector_length <- c(length(sample_name_vec), length(gene_name_vec), length(variant_vec), length(variant_type_vec),
                   length(allele_change_vec), length(conf_grading_vec), length(drug_name_vec))
    
    dr_res_result_list <- list(vector_length, sample_name_vec, drug_name_vec, conf_grading_vec, gene_name_vec, variant_vec, variant_type_vec, variant_pos_vec, allele_change_vec)
    return(dr_res_result_list)
    }


### Check that all lengths value in vector length are equal. If true, fill the empty dataframe with columns with corresponding vectors
fill_res_data_frame <- function(dr_res_df){
    
    ## Instantiate an empty dataframe
    empty_res_dataframe <- tibble(sample_name = character(), gene_name = character(), variant = character(),
                         variant_type = character(), variant_pos = character(), allele_change = character(), 
                         conf_grading = character(), associated_drug = character())

    dr_res_result_list <- fill_res_vectors(dr_res_df)
    #print(dr_res_result_list[[vector_length]])
    if (length(unique(dr_res_result_list[[1]])) == 1){
        #print("Correct")
        empty_res_dataframe <- add_row(empty_res_dataframe, .rows = length(dr_res_result_list[[2]]))
        empty_res_dataframe$sample_name <- dr_res_result_list[[2]]
        empty_res_dataframe$gene_name <- dr_res_result_list[[5]] 
        empty_res_dataframe$variant <- dr_res_result_list[[6]]
        empty_res_dataframe$variant_type <- dr_res_result_list[[7]]
        empty_res_dataframe$variant_pos <- dr_res_result_list[[8]]
        empty_res_dataframe$allele_change <- dr_res_result_list[[9]]
        empty_res_dataframe$conf_grading <- dr_res_result_list[[4]]
        empty_res_dataframe$associated_drug <- dr_res_result_list[[3]]
    } else {
        #print("something is wrong somewhere")
    }
    filled_res_data_frame <- empty_res_dataframe
    return(filled_res_data_frame)

}

#print(fill_res_data_frame(dr_res_table))
## PROCESS LINEAGE DATA
### Extract candidates from lineage_profile_df and append to empty vectors

fill_lin_vectors <- function(lineage_profile_df){
    lin_sample_names <- names(lineage_profile_df)
    lin_names_vec <- character()
    reference_vec <- character()
    lineage_vec <- character()
    sub_lineage_vec <- character()
    var_position_vec <- character()
    lin_allele_change_vec <- character()
    for (sample_name in lin_sample_names){
        for (reference in names(lineage_profile_df[[sample_name]])){
            for (lineage in names(lineage_profile_df[[sample_name]][[reference]])){
                for (sub_lineage in names(lineage_profile_df[[sample_name]][[reference]][[lineage]])){
                    for (var_pos in names(lineage_profile_df[[sample_name]][[reference]][[lineage]][[sub_lineage]])){
                        for (allele_change in lineage_profile_df[[sample_name]][[reference]][[lineage]][[sub_lineage]][[var_pos]]){
                            lin_names_vec <- c(lin_names_vec, sample_name)
                            reference_vec <- c(reference_vec, reference)
                            lineage_vec <- c(lineage_vec, lineage)
                            sub_lineage_vec <- c(sub_lineage_vec, sub_lineage)
                            var_position_vec <- c(var_position_vec, var_pos)
                            lin_allele_change_vec <- c(lin_allele_change_vec, allele_change)
                        }
                    }
                }
            }
        }
    }

    lineage_vector_length <- c(length(lin_names_vec), length(reference_vec), length(lineage_vec), length(sub_lineage_vec), 
                           length(var_position_vec),
                           length(lin_allele_change_vec))
    lin_result_list <- list(lineage_vector_length, lin_names_vec, reference_vec, lineage_vec, sub_lineage_vec, var_position_vec, lin_allele_change_vec)
    return(lin_result_list)
}

### Check that all lengths value in vector length are equal. If true, fill the empty dataframe with columns with corresponding vectors
fill_lin_dataframe <- function(lineage_profile_df){
    
    # Instatiate an empty dataframe
    empty_lin_dataframe <- tibble(sample_name = character(), reference = character(), lineage = character(), 
                                   sub_lineage = character(), allele_change = character(), var_position = character())

    lin_result_list <- fill_lin_vectors(lineage_profile_df)
    if (length(unique(lin_result_list[[1]])) == 1){
        empty_lin_dataframe <- add_row(empty_lin_dataframe, .rows = length(lin_result_list[[2]]))
        empty_lin_dataframe$sample_name <- lin_result_list[[2]]
        empty_lin_dataframe$reference <- lin_result_list[[3]]
        empty_lin_dataframe$lineage <- lin_result_list[[4]]
        empty_lin_dataframe$sub_lineage <- lin_result_list[[5]]
        empty_lin_dataframe$allele_change <- lin_result_list[[7]]
        empty_lin_dataframe$var_position <- lin_result_list[[6]]
    }

    filled_lin_dataframe <- empty_lin_dataframe

    return(filled_lin_dataframe)
}