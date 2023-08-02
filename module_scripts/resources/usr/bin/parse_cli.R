#!/bin/Rscript

library(rjson)
library(tidyverse)

source("/home/dfgmrtc/Project/parse_JSON.R")






## Load json files into R
dr_res_table <- rjson::fromJSON(file="/home/dfgmrtc/Practice/run_folder/n_animal_dr_res.json")
dr_res_int_table <- rjson::fromJSON(file="/home/dfgmrtc/Practice/run_folder/n_animal_dr_res_int.json")
lineage_profile_table <- rjson::fromJSON(file="/home/dfgmrtc/Practice/run_folder/n_animal_lineage.json")

#print(dr_res_table)
## Convert dr_res_table into a fully fledged dataframe with each row having a unique observation
dr_res_all_df <- fill_res_data_frame(dr_res_table)
#print(dr_res_all_df)
## Convert dr_res_int_table into a fully fledged dataframe with each row having a unique observation
dr_res_int_df <- fill_res_data_frame(dr_res_int_table)

## Convert lineage_profile_table into a fully fledged dataframe with each row having a unique observation
lin_prof_df <- fill_lin_dataframe(lineage_profile_table)


## Merge dr_res_all_df with dr_res_int_df

dr_res_full_df <- rbind(dr_res_all_df, dr_res_int_df)

# Arrange dr_res_full_df by variable "sample_name"
dr_res_full_df <- arrange(dr_res_full_df, sample_name)


#### Further analysis 
## Convert full dr_res profiles into short report with who_res_category

# Step 1: Filter full res_df into a df with conf grading "Assoc_w_R" and "Assoc_w_R_Interim"
res_only_dr_res_full_df <-  dr_res_full_df %>%
    filter(conf_grading == "Assoc_w_R" | conf_grading == "Assoc_w_R_Interim")

#print(res_only_dr_res_full_df)
# Step 2: Create a named_list with key corresponding to "sample_name" and value(s) corresponding to associated drug(s)
# 2a. Instatiate the named_list
sample_name_ass_drug_named_list <- list()
# 2b. Fill up the named list
for (no in 1:nrow(res_only_dr_res_full_df)){
    key <- res_only_dr_res_full_df$sample_name[no]
    value <- res_only_dr_res_full_df$associated_drug[no]
    if(is.null(sample_name_ass_drug_named_list[[key]])){
        sample_name_ass_drug_named_list[[key]] <- value
        #print(sample_name_ass_drug_named_list)
    } else {
        sample_name_ass_drug_named_list[[key]] <- c(sample_name_ass_drug_named_list[[key]], value)
    }

}

# Step 3: Create a vector of unique sample names using full_res_dataframe
all_sample_names <- character()

for (no in 1:nrow(dr_res_full_df)){
    all_sample_names <- c(all_sample_names, dr_res_full_df$sample_name[no])
}

all_sample_names <- unique(all_sample_names)
print(length(all_sample_names))
# Interstep: Create drug resistance patterns that will be used in the next step
#RR <- "^(?!.*INH.*RIF)" # tweaked negative lookahead
#INH-mono <- ".*INH(?!.*RIF)" # Negative lookahead
#EMB-mono <- ".*EMB(?!.*INH.*RIF)" # Negative lookahead
#PZA-mono[["Pyrazinamide-Monoresistant"]] <- "^(?!.*INH.*RIF).*PZA" # tweaked negative lookahead
#MDR <- ".*INH.*RIF"
#"Polydrug-Resistant" = "^(?!.*INH.*RIF).*EMB.*PZA" # tweaked negative lookahead
#pre-XDRf <- "^(?=.*INH.*RIF)(.*MXF|.*CFZ|.*LEV)" # tweaked positive lookahead
#pre-XDRi <- "^(?=.*INH.*RIF)(.*AMI|.*CAP|.*KAN)" # tweaked positive lookahead
#XDR <- "^(?=(.*CFZ|.*INH.*RIF)|(.*INH.*MXF.*RIF)|(.*INH.*LEV.*RIF))(.*AMI|.*CAP|.*KAN)" # tweaked positive lookahead

named_pattern_list <- list("Rifampicin-Resistant" = "^(?!.*INH.*RIF).*RIF", "Isoniazid-Monoresistant" = "^(?!.*EMB.*PZA.*RIF).*INH", 
"Ethambutol-Monoresistant" = ".*EMB(?!(.*INH.*RIF)|(.*INH)|(.*RIF))", "Pyrazinamide-Monoresistant" = "^(?!(.*INH.*RIF)|(.*INH)|(.*RIF)).*PZA", "Polydrug-Resistant" = "^(?!(.*INH.*RIF)|(.*INH)|(.*RIF)).*EMB.*PZA",
"Multidrug-Resistant " = ".*INH.*RIF", "Pre-Extensively Drug-Resistant" = "^(?=.*INH.*RIF)(.*MXF|.*CFZ|.*LEV)", "Pre-Extensively Drug-Resistant" = "^(?=.*INH.*RIF)(.*AMI|.*CAP|.*KAN)", 
"Extensively Drug-Resistant" = "^(?=(.*CFZ.*INH.*RIF)|(.*INH.*MXF.*RIF)|(.*INH.*LEV.*RIF))(.*AMI|.*CAP|.*KAN)")

unnamed_pattern_list <- list("^(?!.*INH.*RIF).*RIF", "^(?!.*EMB.*PZA.*RIF).*INH", ".*EMB(?!(.*INH.*RIF)|(.*INH)|(.*RIF))", "^(?!(.*INH.*RIF)|(.*INH)|(.*RIF)).*PZA", "^(?!(.*INH.*RIF)|(.*INH)|(.*RIF)).*EMB.*PZA", ".*INH.*RIF", "^(?=.*INH.*RIF)(.*MXF|.*CFZ|.*LEV)", "^(?=.*INH.*RIF)(.*AMI|.*CAP|.*KAN)", "^(?=(.*CFZ.*INH.*RIF)|(.*INH.*MXF.*RIF)|(.*INH.*LEV.*RIF))(.*AMI|.*CAP|.*KAN)")

any_match <- function(pattern, ass_drugs){
    grepl(pattern, ass_drugs, perl=TRUE)
}

# Step 4: Create WHO dr_res category matching function
dr_res_matching_func <- function(unnamed_pattern_list, ass_drugs, named_pattern_list){
    res_def_msg <- character()
    if(all(!unlist(lapply(unnamed_pattern_list, any_match, ass_drugs)))){
        res_def_msg_nill <- "Nill"
        return(res_def_msg_nill)
    } else {
        
        for (res_def in names(named_pattern_list)){
            if (grepl(named_pattern_list[[res_def]], ass_drugs, perl=TRUE)){
            res_def_msg <- c(res_def_msg, res_def)
            } 
        }
        
        return(res_def_msg)
    }

    
}
#ass_drugs <- c("BDQ, CFZ, DLM, ETH, EMB, INH, LZD, PZA, RIF, STM")
#print(dr_res_matching_func(unnamed_pattern_list, ass_drugs, named_pattern_list))

# Step 4: Create a new empty dataframe that contains a short drug resistance summary of all samples, one sample per row 
dr_res_short_df <- tibble(sample_name = character(), ass_res_drugs = character(), who_drug_res_profile = character())
dr_res_short_df <- add_row(dr_res_short_df, .rows = length(all_sample_names))
for (no in 1:length(all_sample_names)){
    if (all_sample_names[no] %in% names(sample_name_ass_drug_named_list)){ # i.e if a sample name is present in the keys of our drug res names list
        dr_res_short_df$sample_name[no] <- all_sample_names[no]
        ass_drugs <- paste(unique(sample_name_ass_drug_named_list[[all_sample_names[no]]]), collapse = ", ")  # Collaps all drugs associated with that sample_name in our drug res named list to a comma separated list
        dr_res_short_df$ass_res_drugs[no] <- sprintf("Resistant to %s", ass_drugs) # Include the ass comma separated list into the row corresponding to the ass_drugs column
        print(ass_drugs)
        dr_res_short_df$who_drug_res_profile[no] <- dr_res_matching_func(unnamed_pattern_list, ass_drugs, named_pattern_list)
        #print(dr_res_short_df)
    } else {
        dr_res_short_df$sample_name[no] <- all_sample_names[no]
        dr_res_short_df$ass_res_drugs[no] <- c("Susceptible to all drugs")
        dr_res_short_df$who_drug_res_profile[no] <- c("Drug-Susceptible")
    }
}

#options(dyplr.print_max = Inf)

#as.data.frame(dr_res_short_df)