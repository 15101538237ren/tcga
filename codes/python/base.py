# -*- coding:utf-8 -*-
import os, math, re, json, time, pickle, random
import numpy as np

base_dir = os.path.dirname(os.path.dirname(os.getcwd()))

#first level dir
data_dir = os.path.join(base_dir, "data")
global_files_dir = os.path.join(base_dir, "global_files")
figure_dir = os.path.join(base_dir, "figures")

#second level dir
raw_data_dir = "/Users/Ren/PycharmProjects/tcga_raw_data"

#third level dir
dna_methy_data_dir = os.path.join(raw_data_dir, "dna_methy_data")
rna_data_dir = os.path.join(raw_data_dir, "rna")
snv_data_dir = os.path.join(raw_data_dir, "snv")
cnv_data_dir = os.path.join(raw_data_dir, "cnv")
GRCh38_dir = os.path.join(raw_data_dir, "GRCh38")

intermediate_file_dir = os.path.join(data_dir, "intermediate_file")

#fourth level dir
methy_pkl_dir = os.path.join(intermediate_file_dir, "methy_pkl")
manifest_dir = os.path.join(intermediate_file_dir, "manifest")
metadata_dir = os.path.join(intermediate_file_dir, "metadata")

methy_intermidiate_dir = os.path.join(intermediate_file_dir, "methy_intermidiate")
snv_intermidiate_dir = os.path.join(intermediate_file_dir, "snv_intermidiate")
dirs = [methy_pkl_dir, methy_intermidiate_dir, snv_intermidiate_dir]
for dir_name in dirs:
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

#global vars
tumor_stages = ["normal","i","ia","ib","ii","iia","iib","iic","iii","iiia","iiib","iiic","iv","iva","ivb","ivc","x","not reported"]
tumor_stages_n = ["i","ia","ib","ii","iia","iib","iic","iii","iiia","iiib","iiic","iv","iva","ivb","ivc","x"]
tumor_stage_convert = {"normal":"normal","i":"i","ia":"i","ib":"i","ii":"ii","iia":"ii","iib":"ii","iic":"ii","iii":"iii","iiia":"iii","iiib":"iii","iiic":"iii","iv":"iv","iva":"iv","ivb":"iv","ivc":"iv","x":"x","not reported":"not reported"}
merged_stage = ["normal","i","ii","iii","iv","x","not reported"]
merged_stage_n = ["i","ii","iii","iv","x"]

all_cancer_names = ["BRCA", "COAD", "LIHC", "LUAD", "LUSC","BLCA" ,"ESCA","HNSC" ,"KIRC", "KIRP", "PAAD", "READ", "THCA", "STAD","LGG","OV","GBM","LAML", "PRAD","UCEC","SARC", "UVM","CESC", "DLBC"]
cancer_names = ["BRCA", "COAD", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "THCA"]

genome_gene_path = os.path.join(global_files_dir, "gene_with_protein_product.tsv")

def read_whole_genenames(file_path):
    genes = []
    alias_dict = {}
    now_file = open(file_path,'r')
    lines = now_file.readline().split("\r")
    for line in lines:
        contents = line.split("\t")
        gene_name = contents[0]
        if len(contents) > 1:
            alias_dict[gene_name] = gene_name
            alias = contents[1].split("|")
            if len(alias) and alias[0] != "":
                for alias_name in alias:
                    alias_dict[alias_name] = gene_name
        if len(contents) > 2:
            previous_symbols = contents[2].split("|")
            if len(previous_symbols) and previous_symbols[0]!="":
                for previous_symbol in previous_symbols:
                    alias_dict[previous_symbol] = gene_name
        genes.append(gene_name)
    now_file.close()
    return [genes, alias_dict]

[GENOME, alias_dict] = read_whole_genenames(genome_gene_path)
print "length whole genes %d" % len(GENOME)

# 向target_file_path中写target_list的值, 如果index_included=True,则第一列为自增索引, 第二列为value
def write_tab_seperated_file_for_a_list(target_file_path, target_list, index_included=True, sep="\t",line_end = "\n"):
    with open(target_file_path, "w") as outfile:
        if index_included:
            outfile.write(line_end.join([str(gidx+1) + sep + str(value) for gidx, value in enumerate(target_list)]))
        else:
            outfile.write(line_end.join([str(value) for value in target_list]))
        print "write %s successful" % target_file_path

# 读取tab分隔的文件(input_file_path) 的第target_col_index, 返回该列的所有值到一个list
def read_tab_seperated_file_and_get_target_column(target_col_index, input_file_path, sep="\t",line_end = "\n"):
    ret_value_list = []
    with open(input_file_path, "r") as input_file:
        line = input_file.readline()
        while line:
            line_contents = line.split(sep)
            led = line_contents[target_col_index].strip(line_end)
            ret_value_list.append(led)
            line = input_file.readline()
    return ret_value_list
def read_alias_file_and_output_keys_list(input_fp):
    rtn_keys = []
    with open(input_fp,'r') as input_file:
        line = input_file.readline()
        while line:
            contents = line.split("\t")
            gene_name = contents[0]
            rtn_keys.append(gene_name)
            if len(contents) > 1:
                content1 = contents[1].strip("\n")
                alias = content1.split("|")
                if len(alias) and alias[0] != "":
                    for alias_name in alias:
                        rtn_keys.append(alias_name)
            line = input_file.readline()
    return rtn_keys
#input_onco_fp: filepath of onco_gene_file input_tsg_fp: filepath of tumor_suppressed_gene_file, gene_category(0: other, 1: onco, 2: tsg)
def label_TSG_or_OncoGene(input_onco_fp, input_tsg_fp):
    gene_categorys = [0 for item in GENOME]
    onco_keys = read_alias_file_and_output_keys_list(input_onco_fp)
    tsg_keys = read_alias_file_and_output_keys_list(input_tsg_fp)
    for gidx, item in enumerate(GENOME):
        if item in onco_keys:
            gene_categorys[gidx] = 1
        elif item in tsg_keys:
            gene_categorys[gidx] = 2
    return gene_categorys
origin_onco_fp = os.path.join(global_files_dir, "OncoGene_698.tsv")
origin_tsg_fp = os.path.join(global_files_dir, "TSG_1018.tsv")

gene_categorys = label_TSG_or_OncoGene(origin_onco_fp, origin_tsg_fp)

#生成全局统一的gene_index_file,列分别是:gene_idx, gene_name, is_cgi_contained, is_TF_gene, gene_category(0: other, 1: onco, 2: tsg)
def generate_gene_index(gene_idx_path):



    with open(gene_idx_path,"w") as gene_idx_file:
        gene_idx_file.write("\n".join([str(gidx+1) + "\t" + gene for gidx, gene in enumerate(GENOME)]))

#some global variables
gene_idx_path = os.path.join(global_files_dir, "gene_idx.txt")
if not os.path.exists(gene_idx_path):
    generate_gene_index(gene_idx_path)

if __name__ == '__main__':
    pass