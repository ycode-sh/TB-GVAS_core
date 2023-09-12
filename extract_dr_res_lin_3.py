#!/bin/python
import argparse
import re
import json
import os
import csv


from  extract_dr_res_lin_1 import make_per_drug_res_dict, define_lin_dict, process_many_vcf_files
from extract_dr_res_lin_2 import Assign_drug_res, Assign_lineage, Confidence_grading, Sample_class, remove_empty_conf_grad_drug_name, create_table_from_filtered_dict

################# PARSE COMMAND LINE ARGUMENTS

parser = argparse.ArgumentParser()
parser.add_argument("all_files", type = str, nargs = "*")

args = parser.parse_args()

command_line_files = args.all_files

vcf_file_list = []
drug_resistance_file_list = []
lineage_file_list = []
sample_name_list = []
drug_name_list = []

sample_name_list_list = []

for any_file in command_line_files:
    if re.search("sample_[0-9]*[i]?[a-z]?_bt_intersect.vcf", any_file): 
        vcf_file_list.append(any_file)
        file_name = os.path.basename(any_file.split("_bt_intersect.vcf", 1)[0])  
        sample_name_list.append(file_name)
    elif re.search("sorted_[A-Z]*.tsv", any_file): 
        drug_resistance_file_list.append(any_file)
        file_name = any_file.split("sorted_", 1)[1].split(".tsv", 1)[0]
        drug_name_list.append(file_name)
    elif re.search("lineage_snp_[0-9]*.tsv", any_file):
        lineage_file_list.append(any_file)
    else:
        print("No file matches")

############ RUN FUNCTIONS ##########################

my_per_drug_res_dict = make_per_drug_res_dict(drug_resistance_file_list)
general_lin_dict = define_lin_dict(lineage_file_list)
per_sample_drug_dict, per_sample_lineage_dict  = process_many_vcf_files(vcf_file_list)



################ VARIABLES ##################


confidence_grading_list = ["Assoc_w_R", "Assoc_w_R_Interim", "combo", "Not_assoc_w_R_Interim","Not_assoc_w_R", "Uncertain_significance"]


main_lineages_list = ["1", "2", "3", "4", "5", "6", "7", "asian_african_1", "asia_ancestral_2", "BOV_AFRI", "BOV", "asia_ancestral_3", "central_asia", 
                                      "asian_african_2/RD142", "pacific_RD150"]

lineage_reference_list = ["coll_2014", "freschi_2020", "shitikov_2017"]


#################### RUN CLASSES #########################

inst = Confidence_grading(my_per_drug_res_dict, ["AMI"], ["Assoc_w_R"])

drug_resistance_instance = Assign_drug_res(my_per_drug_res_dict, drug_name_list, confidence_grading_list, per_sample_drug_dict, sample_name_list)


lineage_ass_instance = Assign_lineage(general_lin_dict, lineage_reference_list, main_lineages_list, per_sample_lineage_dict, sample_name_list)



print(drug_resistance_instance.long_table_dr_call)

print(drug_resistance_instance.long_table_dr_int_var_call)
print(lineage_ass_instance.lineage_long_call)

with open("n_animal_dr_res.json", "w") as json_file:
    json.dump(drug_resistance_instance.filtered_attr_for_long_dr_call, json_file, indent=4)

with open("n_animal_dr_res_int.json", "w") as json_res_int_file:
    json.dump(drug_resistance_instance.samplename_drug_led_int_var_attr, json_res_int_file, indent= 4)

with open("n_animal_lineage.json", "w") as json_lineage_file:
    json.dump(lineage_ass_instance.filtered_for_long_lin_call, json_lineage_file, indent = 4)


## Write sample name in a list
#for name in sample_name_list:
#    sample_name_list_list.append([name])

#file_path = 'sample_names.csv'
#with open(file_path, 'w', newline = '') as csv_file:
#    writer = csv.writer(csv_file)
#    writer.writerows(sample_name_list_list) """


#print(drug_resistance_instance(drug_names = ["RIF", "EMB", "AMI"], sample_names = ['Benchmark_samples2/trimmed_sample_1', 'Benchmark_samples2/trimmed_sample_6'], confidence_list = ["Assoc_w_R"], return_only_p_c_call = True))
