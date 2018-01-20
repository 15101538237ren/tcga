# -*- coding:utf-8 -*-
import os, math, re, json, time, pickle, random
import numpy as np

base_dir = os.path.dirname(os.path.dirname(os.getcwd()))

#first level dir
data_dir = os.path.join(base_dir, "data")
global_files_dir = os.path.join(base_dir, "global_files")
codes_files_dir = os.path.join(base_dir, "codes")
figure_dir = os.path.join(base_dir, "figures")

#second level dir
raw_data_dir ="/disk/tcga_raw_data" # "/Users/Ren/PycharmProjects/tcga_raw_data" "/Volumes/Elements/tcga_raw_data"#

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
methy_matlab_data_dir = os.path.join(intermediate_file_dir, "methy_matlab_data")
snv_intermidiate_dir = os.path.join(intermediate_file_dir, "snv_intermidiate")
rna_intermidiate_dir = os.path.join(intermediate_file_dir, "rna_intermidiate")
methy_entropy_dir = os.path.join(intermediate_file_dir, "methy_entropy")
methy_corr_dir = os.path.join(intermediate_file_dir, "methy_corr")

methy_figure_dir = os.path.join(figure_dir, "methy_scatter")
methy_mean_std_dir = os.path.join(intermediate_file_dir, "methy_mean_std")
methy_manifest_path = os.path.join(global_files_dir, "methy_24_cancer_manifest.tsv")
methy_pvalue_dir = os.path.join(intermediate_file_dir, "methy_pvalue")
common_sample_cnt_dir = os.path.join(intermediate_file_dir, "common_sample_cnt")
common_patient_data_dir = os.path.join(intermediate_file_dir, "common_patients_data")

dirs = [methy_pkl_dir, methy_intermidiate_dir,methy_matlab_data_dir, snv_intermidiate_dir, methy_mean_std_dir, methy_entropy_dir, methy_corr_dir, methy_pvalue_dir, common_sample_cnt_dir,common_patient_data_dir]

for dir_name in dirs:
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

promoter_length = 2000
pvalue_significant_threshold = -10
#global vars
tumor_stages = ["normal","i","ia","ib","ii","iia","iib","iic","iii","iiia","iiib","iiic","iv","iva","ivb","ivc","x","not reported"]
tumor_stage_convert = {"normal":"normal","i":"i","ia":"i","ib":"i","ii":"ii","iia":"ii","iib":"ii","iic":"ii","iii":"iii","iiia":"iii","iiib":"iii","iiic":"iii","iv":"iv","iva":"iv","ivb":"iv","ivc":"iv","x":"x","not reported":"not reported"}
merged_stage = ["normal","i","ii","iii","iv","x","not reported"]

methy_and_rna_merged_stages = ["normal","i","ii","iii","iv"]
methy_and_rna_stages = ["normal","i","ia","ib","ii","iia","iib","iic","iii","iiia","iiib","iiic","iv","iva","ivb","ivc"]

mutation_merged_stage = ["i","ii","iii","iv","not reported"]
mutation_stage = ["i","ia","ib","ii","iia","iib","iic","iii","iiia","iiib","iiic","iv","iva","ivb","ivc","not reported"]

all_cancer_names = ["BRCA", "COAD", "LIHC", "LUAD", "LUSC","BLCA" ,"ESCA","HNSC" ,"KIRC", "KIRP", "PAAD", "READ", "THCA", "STAD","LGG","OV","GBM","LAML", "PRAD","UCEC","SARC", "UVM","CESC", "DLBC"]
cancer_names = ["BRCA","COAD", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "THCA"] #

SNP_Ins_Del_classification = {"SNP":0, "INS": 1, "DEL":2}
mutation_classification = {"Frame_Shift_Ins":1, "In_Frame_Ins":2, "Frame_Shift_Del":3, "In_Frame_Del":4,
                           "Missense_Mutation":5, "Translation_Start_Site":6 , "Splice_Region":7, "Splice_Site":8,
                           "3'UTR":9, "5'UTR":10, "3'Flank":11, "5'Flank":12,
                           "Nonsense_Mutation":13, "Nonstop_Mutation":14,
                           "RNA":15, "Silent":16, "IGR":17, "Intron":18}

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

# 通过manifest文件中的对应关系,将下载的文件名filename和uuid对应起来,方便互相查询(uuid->filename, filename->uuid)
def connect_filename_to_uuid():
    uuid_to_filename = {}
    filename_to_uuid = {}
    now_file = open(methy_manifest_path,'r')

    #pass the header
    now_file.readline()

    str_pattern = r'([^\t]+)\t([^\t]+)'
    cancer_pattern = r'jhu-usc.edu_([^\.]+)*'
    uuid_dict = {cancer_name:[] for cancer_name in cancer_names}
    file_name_dict={cancer_name:[] for cancer_name in cancer_names}

    line = now_file.readline()
    while line:
        match_p = re.search(str_pattern, line)
        if match_p:
            uuid = match_p.group(1)
            file_name = match_p.group(2)

            uuid_to_filename[uuid] = file_name
            filename_to_uuid[file_name] = uuid
            try:
                cancer_name = re.search(cancer_pattern, file_name).group(1)
                if cancer_name in cancer_names:
                    uuid_dict[cancer_name].append(uuid)
                    file_name_dict[cancer_name].append(file_name)
            except AttributeError:
                print file_name
            # print "%s\t%s" % (uuid, file_name)

        line=now_file.readline()
    now_file.close()
    print "connect_filename_to_uuid called"
    return [uuid_to_filename, filename_to_uuid, uuid_dict, file_name_dict]

#returned global vars from connect_filename_to_uuid()
[uuid_to_filename, filename_to_uuid, uuid_dict, file_name_dict] = connect_filename_to_uuid()

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

#label whether a gene in GENOME contains a CGI, 1 contained, 0 non-contained
def label_cgi_genes(CGI_genenames_filepath):
    gene_cgi_labels =[0 for item in GENOME]
    cgi_gene_names =[item.strip("\"") for item in read_tab_seperated_file_and_get_target_column(0, CGI_genenames_filepath)]
    for gidx, gene_name in enumerate(GENOME):
        if gene_name in cgi_gene_names:
            gene_cgi_labels[gidx] = 1
    return gene_cgi_labels
# for an alias file, output all key and its alias as a list
def read_alias_file_and_output_keys_list(input_fp):
    rtn_keys = []
    with open(input_fp,'r') as input_file:
        line = input_file.readline().strip("\n")
        while line:
            contents = line.split("\t")
            gene_name = contents[0]
            rtn_keys.append(gene_name)
            if len(contents) > 1:
                content1 = contents[1]
                alias = content1.split("|")
                if len(alias) and alias[0] != "":
                    for alias_name in alias:
                        rtn_keys.append(alias_name)
            line = input_file.readline().strip("\n")
    return rtn_keys

def match_gene_idx_for_genenames(input_fp, out_gidxs_fp):
    gnames = read_tab_seperated_file_and_get_target_column(0, input_fp,line_end='\r\n')
    gidxs = []
    GENOME_idx_dict = {gname: gidx + 1 for gidx, gname in enumerate(GENOME)}
    for gidx, gname in enumerate(gnames):
        try:
            gsymbol = alias_dict[gname]
            gidx_corresponding = GENOME_idx_dict[gsymbol]
        except KeyError,e:
            gidx_corresponding = -1
        gidxs.append(gidx_corresponding)

    with open(out_gidxs_fp,'w') as out_gidxs_file:
        ltws = []
        for gidx, gname in enumerate(gnames):
            ltws.append('\t'.join([str(gidx + 1), str(gidxs[gidx]), gname]))
        out_gidxs_file.write('\n'.join(ltws))
    print "write %s successful!" % out_gidxs_fp

#input_onco_fp: filepath of onco_gene_file input_tsg_fp: filepath of tumor_suppressed_gene_file, gene_category(0: other, 1: onco, 2: tsg)
def label_TSG_or_OncoGene(input_onco_fp, input_tsg_fp):
    gene_categorys = [0 for item in GENOME]
    onco_keys = read_alias_file_and_output_keys_list(input_onco_fp)
    tsg_keys = read_alias_file_and_output_keys_list(input_tsg_fp)
    for gidx, item in enumerate(GENOME):
        gc_nm = 0
        if item in onco_keys:
            gene_categorys[gidx] = 1
            gc_nm += 1
        if item in tsg_keys:
            gene_categorys[gidx] = 2
            gc_nm += 1
        if gc_nm == 2:
            gene_categorys[gidx] = 3
    return gene_categorys

def generate_gene_position_info(gene_body_fp):
    gene_pos_labels = [[-1, -1, -1, -1] for item in GENOME]
    gene_idxs_dict = {item : gidx for gidx, item in enumerate(GENOME)}
    chr_dict = {str(i): i for i in range(1, 23)}
    chr_dict["X"] = 23
    chr_dict["Y"] = 24

    for line in file(gene_body_fp):
        line_items = line.strip("\n").split("\t")
        gene_name = line_items[0].strip('"')
        if gene_name in GENOME:
            chr_str = line_items[1].strip('"')
            if chr_str in chr_dict.keys():
                chr_no = chr_dict[chr_str]
                strand = 1 if line_items[4].strip('"') == "+" else 0
                start = int(line_items[2])
                end = int(line_items[3])
                gene_pos_labels[gene_idxs_dict[gene_name]][:] =  [chr_no, start, end ,strand]
    return gene_pos_labels

def extract_gene_info_from_v22_gtf(v22_gtf_fp):
    gene_pos_labels = [[-1, -1, -1, -1] for item in GENOME]
    gene_idxs_dict = {item: gidx for gidx, item in enumerate(GENOME)}
    chr_dict = {str(i): i for i in range(1, 23)}
    chr_dict["X"] = 23
    chr_dict["Y"] = 24

    gtf_file = open(v22_gtf_fp, "r")

    for i in range(5):
        gtf_file.readline()
    line = gtf_file.readline()

    while line:
        line_contents = line.split("\t")
        feature = line_contents[2]
        if feature == "gene":
            groups = line_contents[-1].split(";")[0:-1]
            group_dict = {}
            for item in groups:
                [k, v] = item.strip().split(" ")
                group_dict[k] = v.replace('\"', "")
            chr_str = line_contents[0].replace("chr", "")
            if "gene_type" in group_dict.keys() and group_dict["gene_type"] == "protein_coding" and "gene_name" in group_dict.keys() and chr_str in chr_dict.keys():
                    gene_name = group_dict["gene_name"]
                    if gene_name in GENOME:
                        chr_no = chr_dict[line_contents[0].replace("chr", "")]
                        start = int(line_contents[3])
                        end = int(line_contents[4])
                        strand = 1 if line_contents[6] == "+" else 0
                        gene_pos_labels[gene_idxs_dict[gene_name]][:] = [chr_no, start, end, strand]
        line = gtf_file.readline()
    return gene_pos_labels
#label whether a gene of GENOME is a TF gene
def label_TF_genes(input_tf_fp):
    tf_labels = [0 for item in GENOME]
    tf_genes = read_tab_seperated_file_and_get_target_column(0,tf_genes_fp)
    for tidx, gene_name in enumerate(GENOME):
        if gene_name in tf_genes:
            tf_labels[tidx] = 1
    return tf_labels

CGI_genenames_fp = os.path.join(global_files_dir, "gene_names_with_CGI.txt")
gene_cgi_labels = label_cgi_genes(CGI_genenames_fp)

tf_genes_fp = os.path.join(global_files_dir, "TF_gene_list.txt")
tf_gene_labels = label_TF_genes(tf_genes_fp)

origin_onco_fp = os.path.join(global_files_dir, "OncoGene_698.tsv")
origin_tsg_fp = os.path.join(global_files_dir, "TSG_1018.tsv")

gene_categorys = label_TSG_or_OncoGene(origin_onco_fp, origin_tsg_fp)

origin_onco_fp_vogelstein = os.path.join(global_files_dir, "Oncogene_Vogelstein.txt")
origin_tsg_fp_vogelstein = os.path.join(global_files_dir, "TSG_Vogelstein.txt")

gene_categorys_vogelstein = label_TSG_or_OncoGene(origin_onco_fp_vogelstein, origin_tsg_fp_vogelstein)

msk_341_fp = os.path.join(global_files_dir, "msk_impact_341_genes.txt")
msk_341_labels = label_cgi_genes(msk_341_fp)

msk_410_fp = os.path.join(global_files_dir, "msk_impact_410_genes.txt")
msk_410_labels = label_cgi_genes(msk_410_fp)

gene_body_fp = os.path.join(global_files_dir, "human_gene_bodys.tsv")
gene_pos_labels = generate_gene_position_info(gene_body_fp)

# v22_gtf_fp = os.path.join(global_files_dir, "gencode.v22.annotation.gtf")
# v22_gene_pos_labels = extract_gene_info_from_v22_gtf(v22_gtf_fp)

gene_pos_labels_used = gene_pos_labels
#生成全局统一的gene_index_file,列分别是:gene_idx, gene_name, is_cgi_contained, is_TF_gene, gene_category(0: other, 1: onco, 2: tsg 3: onco_and_tsg)
def generate_gene_index(gene_idx_fp, gene_label_fp):
    with open(gene_idx_fp,"w") as gene_idx_file:
        ltws = []
        for gidx, gene in enumerate(GENOME):
            ltw = "\t".join([str(gidx + 1), gene])
            ltws.append(ltw)
        gene_idx_file.write("\n".join(ltws))
    with open(gene_label_fp,"w") as gene_label_file:
        ltws = []
        for gidx, gene in enumerate(GENOME):
            chr_no, start, end, strand = gene_pos_labels_used[gidx]
            ltw = "\t".join([str(gidx + 1), str(gene_cgi_labels[gidx]), str(tf_gene_labels[gidx]), str(gene_categorys[gidx]), str(gene_categorys_vogelstein[gidx]), str(msk_341_labels[gidx]), str(msk_410_labels[gidx]), str(chr_no), str(start), str(end), str(strand)])
            ltws.append(ltw)
        gene_label_file.write("\n".join(ltws))
    print "generate_gene_index successful at %s" % gene_idx_fp

with open(os.path.join(global_files_dir, "mutation_classification.txt"),"w") as mutation_classification_file:
    sorted_dict = sorted(mutation_classification.items(), key=lambda d: d[1])
    for k,v in sorted_dict:
        mutation_classification_file.write(str(v) + "\t" + k + "\n")

#some global variables
gene_idx_path = os.path.join(global_files_dir, "gene_idx.txt")
gene_label_path = os.path.join(global_files_dir, "gene_label.dat")
if __name__ == '__main__':
    generate_gene_index(gene_idx_path, gene_label_path)
    # yc_geneset_fp = os.path.join(global_files_dir,'yucheng_candidate_set.txt')
    # out_gidxs_fp = os.path.join(codes_files_dir,'matlab','yucheng_candidate_set.ind')
    # match_gene_idx_for_genenames(yc_geneset_fp,out_gidxs_fp)

