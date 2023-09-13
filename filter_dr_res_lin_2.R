#!/bin/Rscript

suppressPackageStartupMessages({
    library(rjson)
    library(tidyverse)
    library(argparser)
})


source("/home/dfgmrtc/Workflows/wf-vdg/bin/filter_dr_res_lin_1.R")


## Load json files into R
parser <- arg_parser("Parse command-line arguments")
parser <- add_argument(parser, "--files", type = "character", nargs = "*", help = "Files not correctly parsed")

arguments <- parse_args(parser)

x <- unlist(strsplit(arguments$files, split = ","))


if (!is.null(arguments$files)){
  for (item in unlist(strsplit(arguments$files, split = ","))){
    #print(item)
    if(item == "dr_res.json"){
      dr_res_table <- rjson::fromJSON(file=item)
    } else if(item == "dr_res_int.json"){
      dr_res_int_table <- rjson::fromJSON(file=item)
    } else if(item == "lineage.json"){
      lineage_profile_table <- rjson::fromJSON(file=item)
    }
  }
}

dr_res_full_df <- merge_dr_files(dr_res_table, dr_res_int_table)
dr_res_short_df <- create_short_dr_profile(dr_res_table, dr_res_int_table)
lineage_full_df <- fill_lin_dataframe(lineage_profile_table)
lineage_short_df <- create_short_lin_profile(lineage_profile_table)
#print(as.data.frame(create_short_dr_profile(dr_res_table, dr_res_int_table)))

write_tsv(dr_res_full_df, "detailed_drug_resistance_profile.tsv")
write_tsv(dr_res_short_df, "short_drug_resistance_profile.tsv")
write_tsv(lineage_full_df, "detailed_lineage_assignments.tsv")

dr_res.json dr_res_int.json lineage.json