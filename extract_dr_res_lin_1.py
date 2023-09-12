#!/bin/python

import re
import os

# Step 2: Prepare annotation files
# 2A
# This function is working as expected. It takes a list containing dr_resistance files
def define_drug_res_dict(drug_resistance_file):
    confidence_grading_summary_list = ["Assoc_w_R", "Assoc_w_R_Interim", "combo", "Not_assoc_w_R", "Not_assoc_w_R_Interim", "Uncertain_significance"]
    with open(drug_resistance_file, 'r') as drug_file:
        drug_res_dict = {conf_grade:{} for conf_grade in confidence_grading_summary_list}
        for line in drug_file:
            if line.startswith("#"):
                pass
            line_data = line.rstrip().rsplit("\t")
            gene_string = str(line_data[1]).split("_", 1)[0]
            var_string =  str(line_data[1]).split("_", 1)[1]
            conf_string = str(line_data[3])
            pos_string = str(line_data[2])
            if gene_string not in drug_res_dict[conf_string].keys():
                drug_res_dict[conf_string].setdefault(gene_string,{})
                if var_string not in drug_res_dict[conf_string][gene_string].keys():
                    drug_res_dict[conf_string][gene_string].setdefault(var_string, [])
                    drug_res_dict[conf_string][gene_string][var_string].append(pos_string)
                else:
                    drug_res_dict[conf_string][gene_string][var_string].append(pos_string)
            else:
                drug_res_dict[conf_string][gene_string].setdefault(var_string, [])
                drug_res_dict[conf_string][gene_string][var_string].append(pos_string)
    return drug_res_dict 

# This function is working as expected
def make_per_drug_res_dict(drug_resistance_file_list):
    per_drug_dr_res_dict = {}
    for file in drug_resistance_file_list:
        drug_acronym = file.rsplit("sorted_")[1].rsplit(".tsv")[0]
        drug_res_dict = define_drug_res_dict(file)
        if drug_acronym not in per_drug_dr_res_dict.keys():
            per_drug_dr_res_dict.setdefault(drug_acronym, {})
            per_drug_dr_res_dict[drug_acronym] = drug_res_dict
    return per_drug_dr_res_dict


# 2B: Pepare lineage annotation files

def define_lin_dict(lineage_file_list):
    lineage_file = lineage_file_list[0]
    lineage_reference_list = ["coll_2014", "freschi_2020", "shitikov_2017"]
    lin_dict = {lin_ref:{} for lin_ref in lineage_reference_list}

    with open(lineage_file, 'r') as lin_file:
        for line in lin_file:
            data = line.rstrip().rsplit("\t")
            lineage_main = str(data[2]).split(".", 1)[0]
            if lineage_main not in lin_dict[data[1]].keys():
                lin_dict[data[1]].setdefault(lineage_main, {})
                lin_dict[data[1]][lineage_main].setdefault(data[2], {})
                lin_dict[data[1]][lineage_main][data[2]].setdefault(str(data[0]), str(data[3]))
                
            else:
                lin_dict[data[1]][lineage_main].setdefault(data[2], {})
                lin_dict[data[1]][lineage_main][data[2]].setdefault(str(data[0]), str(data[3]))
                
    return lin_dict


# 3. Parse VCF files

# The "conver_p_str" function will convert protein SNPs from HGVS format to WHO mutation catalogue format while "convert_c_str" 
# will convert cdna SNPs and INDELs from HGVS format to WHO mutation catalogue format

def convert_p_str(string):
    fp = ''
    nfp = ''
    snp = ''
    vn = ''
    string = string.split(".", 1)[1].replace("*", "!")
    amino_acids_dict = {"Ala": "A", "Arg": "R", "Asn" : "N", "Asp" : "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", 
                    "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", 
                    "Tyr": "Y", "Val": "V", "del": "del", "!": "!"}
    pattern = r'^([a-zA-Z]{1,3})(\d+)([a-zA-Z]{1,3})?(!)?$'

    match = re.match(pattern, string)
    if match:
        if match.group(3) != "dup" and match.group(3) != "fs":
            fp = match.group(1)
            vn = match.group(2)
            if match.group(3) is None and match.group(4) == "!":
                sp = match.group(4)
            else:   
                sp = match.group(3)
            
            nfp = amino_acids_dict[fp]
            snp = amino_acids_dict[sp]
            
    
    return "".join([nfp, vn, snp])

def convert_c_str(string, allele_change, variant_eff):
    snp = ""
    sn_indel = ""
    indels = ""
    dup_type = ""
    dup_sn_indel = ""
    spl_str = string.split(".", 1)[1].replace("*", "")
    ref = allele_change.split("/", 1)[0].lower()
    alt = allele_change.split("/", 1)[1].lower()
    if "dup" in spl_str:     
        numb = str(int(re.findall(r'\d+', spl_str)[0]) + 1) # Added 1 bcos I realized snpEff miscalculates var_pos in frameshift mutations

        if int(len(ref)) - int(len(alt)) == 1:
            dup_type = "del"
        elif int(len(ref)) - int(len(alt)) == -1:
            dup_type = "ins"
        dup_sn_indel = "_".join([numb, dup_type, "1", ref, alt])
        
    elif re.search("^[-]?[0-9]*[a-zA-Z][>][a-zA-Z]$", spl_str):  
        #print("Yes o")
        pattern = r'^([-]?[0-9]*)([a-zA-Z])(>)([a-zA-Z])$'
        match = re.match(pattern, spl_str)
        #print(match)
        if match:
            snp =  "".join([match.group(2).lower(), match.group(1), match.group(4).lower()])

    elif re.search("^[-]?[0-9]*(del)?(ins)?[a-zA-Z]$", spl_str):   #1471926_ins_2   103delG
        #print("Is getting serious")
        pattern = r'^([-]?[0-9]*)((del)?(ins)?)([a-zA-Z])$'
        #print(pattern)
        match = re.match(pattern, spl_str)
        if match:
            if variant_eff == "frameshift_variant":
                modified_pos = str(int(match.group(1)) + 1) # Added 1 bcos I realized snpEff miscalculates var_pos in frameshift mutations
                sn_indel = "_".join([modified_pos, match.group(2), "1", ref, alt])
            else:
                sn_indel = "_".join([match.group(1), match.group(2), "1", ref, alt])
            #print(sn_indel)
    else:        # 355_356delGC
        #print("This is going to be serious")
        re.search("^[-]?[0-9]*_[-]?[0-9]*(del)?(ins)?[a-zA-Z]*$", spl_str)
        spl_str1 = spl_str.split("_", 1)[1]   # I have been struggling with which one to take, the first or second?
        spl_str1 = re.findall(r'\d+', spl_str1)[0]
        pattern = r'^([-]?[0-9]*_[-]?[0-9]*)((ins)?(del)?)([a-zA-Z]*)$'
        match = re.match(pattern, spl_str)
        if match:
            if match.group(2) == "del":
                spl_str2 = spl_str.split("l", 1)[1]
            elif match.group(2) == "ins":
                spl_str2 = spl_str.split("s", 1)[1]
            else:
                raise TypeError("Variant %s type not known" %match.group(2))
            nN = len(spl_str2)
            if variant_eff == "frameshift_variant":
                modified_pos = str(int(spl_str1) + 1)   # Because snpEff miscalculates frameshifts mutation positions
                indels = "_".join([modified_pos, match.group(2), str(nN), ref, alt])
            else:
                indels = "_".join([spl_str1, match.group(2), str(nN), ref, alt])

    if snp:
        return snp
    elif sn_indel:
        return sn_indel
    elif dup_sn_indel:
        return dup_sn_indel
    else:
        return indels

def convert_n_string(string):
    string = string.split(".", 1)[1].replace("*", "")
    int_snps = ""
    if re.search("^[-]?[0-9]*[a-zA-Z][>][a-zA-Z]$", string):
        pattern = r'^([-]?[0-9]*)([a-zA-Z])(>)([a-zA-Z])$'
        match = re.match(pattern, string)
        if match:
            int_snps =  "".join([match.group(2).lower()])
            #print(int_snps)
    #print(int_snps)
    return int_snps

def modify_drug_res_dict(dr_res_variants, var, allele_change, Gene_name, variant_eff, var_pos, exp_drugs):
    details_list = []
    details_list.append(variant_eff)
    details_list.append(var_pos)
    details_list.append(allele_change)
    details_list.append(exp_drugs)
    dr_res_variants[Gene_name].setdefault(var, [])
    dr_res_variants[Gene_name][var].append(details_list)
    
    #dr_res_variants[Gene_name][var].append(variant_eff)
    #dr_res_variants[Gene_name][var].append(var_pos)
    #dr_res_variants[Gene_name][var].append(allele_change)
    #dr_res_variants[Gene_name][var].append(exp_drugs)
    return dr_res_variants

def process_any_string(dr_res_variants:dict, string, allele_change, Gene_name, variant_eff, var_pos, exp_drugs):
    function_level_dict = {"frameshift_variant":"run_c_string_func", "missense_variant":"run_p_string_func",
                            "synonymous_variant":"run_p_string_func", "upstream_gene_variant": "run_c_string_func",
                            "downstream_gene_variant":"run_c_string_func",
                            "intergenic_region":"run_n_string_func", "intragenic_variant":"run_n_string_func"}

    if function_level_dict[variant_eff] == "run_p_string_func": 
        var = convert_p_str(string)
        dr_res_variants = modify_drug_res_dict(dr_res_variants, var, allele_change, Gene_name, variant_eff,  var_pos, exp_drugs)
    elif function_level_dict[variant_eff] == "run_c_string_func":
        var = convert_c_str(string, allele_change, variant_eff)
        dr_res_variants = modify_drug_res_dict(dr_res_variants, var, allele_change, Gene_name, variant_eff,  var_pos, exp_drugs)
    elif function_level_dict[variant_eff] == "run_n_string_func":
        var = convert_n_string(string)
        dr_res_variants = modify_drug_res_dict(dr_res_variants, var, allele_change, Gene_name, variant_eff,  var_pos, exp_drugs)
    
    else:
        print("This varient_effect:%s is not captured by my code" %variant_eff)
    

    return dr_res_variants


def minos_vcf(data, lineage_positions_dict, dr_res_variants):
    if data[9].split(":")[0] == "0/0":  # Only screen heterogeneous variants further 
        pass
    elif float(data[9].split(":")[1].replace('.', '0')) < 10:   # Only screen variants with DP > 10
        pass
    else:
        if data[10] == "lineage_snps":
            if "/".join([data[3], data[4]]) == "/".join([data[15], data[16]]):
                lineage_positions_dict.setdefault(data[1], "/".join([data[3], data[4]]))
        elif data[10] == "amr_regions":
            for item in range(len(data[7].split("=")[1].split(",")[:])):
                allele_change = "/".join([data[3], data[4]])
                ANN = data[7].split("=")[1].split(",")[item].split("|")[:]
                if ANN[3] in data[14]:
                    dr_res_variants.setdefault(ANN[3], {})
                    if ANN[1] in ["frameshift_variant", "upstream_gene_variant", "downstream_gene_variant", "intergenic_region", "intragenic_variant"]:
                        if ANN[9].split(".", 1)[0] == "c":
                            dr_res_variants = process_any_string(dr_res_variants, ANN[9], allele_change, ANN[3], ANN[1], data[1], data[16])
                        elif ANN[9].split(".", 1)[0] == "n":
                            dr_res_variants = process_any_string(dr_res_variants, ANN[9], allele_change, ANN[3], ANN[1], data[1], data[16])
                        else:
                            print("%s is not captured by c and n string code" %{ANN[9].split(".", 1)[0]})
                    elif ANN[1] in ["missense_variant", "synonymous_variant"]:
                        if ANN[10].split(".", 1)[0] == "p":
                            dr_res_variants = process_any_string(dr_res_variants, ANN[10], allele_change, ANN[3], ANN[1], data[1], data[16])
                        else:
                            string_id = ANN[10].split(".", 1)[0]
                            print("%s is not captured by my code" %string_id)


def bcftools_vcf(data, lineage_positions_dict, dr_res_variants):
    if data[9].split(":")[0] == "0":
        pass
    elif float(data[7].split("DP=")[1].split(";")[0]) < 10:
        pass
    else:
        if data[10] == "lineage_snps":
            if "/".join([data[3], data[4]]) == "/".join([data[16], data[17]]):
                lineage_positions_dict.setdefault(data[1], "/".join([data[3], data[4]]))
        elif data[10] == "amr_regions":
            for item in range(len(data[7].split("ANN=")[1].split(",")[:])):
                allele_change = "/".join([data[3], data[4]])
                ANN = data[7].split("ANN=")[1].split(",")[item].split("|")[:]
                if ANN[3] in data[14]:
                    dr_res_variants.setdefault(ANN[3], {})
                    if ANN[1] in ["frameshift_variant", "upstream_gene_variant", "downstream_gene_variant", "intergenic_region", "intragenic_variant"]:
                        if ANN[9].split(".", 1)[0] == "c":
                            dr_res_variants = process_any_string(dr_res_variants, ANN[9], allele_change, ANN[3], ANN[1], data[1], data[16])
                        elif ANN[9].split(".", 1)[0] == "n":
                            dr_res_variants = process_any_string(dr_res_variants, ANN[9], allele_change, ANN[3], ANN[1], data[1], data[16])
                        else:
                            print("%s is not captured by c and n string code" %{ANN[9].split(".", 1)[0]})
                    elif ANN[1] in ["missense_variant", "synonymous_variant"]:
                        if ANN[10].split(".", 1)[0] == "p":
                            dr_res_variants = process_any_string(dr_res_variants, ANN[10], allele_change, ANN[3], ANN[1], data[1], data[16])
                        else:
                            string_id = ANN[10].split(".", 1)[0]
                            print("%s is not captured by my code" %string_id)

def proces_a_vcf_file(a_vcf_file, variant_caller = "bcftools"):    
    lineage_positions_dict = {}
    dr_res_variants = {}
    with open(a_vcf_file, 'r') as vcf_file:
        for line in vcf_file:
            if line.startswith("#"):
                continue
            data = line.rstrip().rsplit("\t")

            if variant_caller == "bcftools":
                bcftools_vcf(data, lineage_positions_dict, dr_res_variants)
            elif variant_caller == "minos":
                minos_vcf(data, lineage_positions_dict, dr_res_variants)
                            
    return lineage_positions_dict, dr_res_variants
                            
                        
def process_many_vcf_files(vcf_file_list):
    per_file_dr_res_dict = {}
    per_file_lineage_dict = {}
    for a_vcf_file in vcf_file_list:
        vcf_file_name = os.path.basename(a_vcf_file.split("_bt_intersect.vcf")[0])
        lineage_positions_list, dr_res_variants =  proces_a_vcf_file(a_vcf_file)
        per_file_dr_res_dict.setdefault(vcf_file_name, {})
        per_file_dr_res_dict[vcf_file_name] = dr_res_variants
        per_file_lineage_dict.setdefault(vcf_file_name, {})
        per_file_lineage_dict[vcf_file_name] = lineage_positions_list
    return per_file_dr_res_dict, per_file_lineage_dict

