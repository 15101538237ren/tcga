# -*- coding:utf-8 -*-
import os, json
import numpy as np
from base import *
fpkm_file_end = ".FPKM.txt"

is_merge_stage = True
dname = "merged_stage" if is_merge_stage else "stage"
# 获得每种癌症的htseq-count数据的文件名列表
def obtain_fpkm_filelist():
    filename_dict = {}
    for cancer_name in cancer_names:
        cancer_data_dir = os.path.join(rna_data_dir, cancer_name)

        output_cancer_dir = os.path.join(rna_intermidiate_dir, dname, cancer_name)
        if not os.path.exists(output_cancer_dir):
            os.makedirs(output_cancer_dir)

        fpkm_file_names = []
        filenames = os.listdir(cancer_data_dir)
        for filename in filenames:
            if filename.endswith(fpkm_file_end):
                fpkm_file_names.append(filename)
        filename_dict[cancer_name] = fpkm_file_names

        outfile_path = os.path.join(output_cancer_dir, cancer_name + "_fpkm_filelist.txt")
        write_tab_seperated_file_for_a_list(outfile_path, fpkm_file_names, index_included=True)
    print "obtain fpkm filelist successful!"
    return [filename_dict]

# 从meta_data文件中获取所有fpkm文件对应的癌症主分期列表,并将其输出到文件cancer_name_stages.txt中
def obtain_stage_info_from_metadata():
    #make output_dir
    for cancer_name in cancer_names:
        output_cancer_dir = os.path.join(rna_intermidiate_dir, dname, cancer_name)
        if not os.path.exists(output_cancer_dir):
            os.makedirs(output_cancer_dir)

        input_path = os.path.join(output_cancer_dir, cancer_name + "_fpkm_filelist.txt")
        output_path = os.path.join(output_cancer_dir, cancer_name + "_stages.txt")
        metadata_file_path = os.path.join(metadata_dir, "rna_" + cancer_name + "_metadata.json")
        fpkm_file_names = read_tab_seperated_file_and_get_target_column(1, input_path)
        file_name_to_stage = {}
        json_obj = json.load(open(metadata_file_path,'r'))
        stage_list = []
        for fpkm_file_name in fpkm_file_names:
            fpkm_file_name_extend = fpkm_file_name + ".gz"
            stage_save = "not reported"
            for obj in json_obj:
                if obj["file_name"] != fpkm_file_name_extend:
                    continue
                else:
                    entity_submitter_id = obj["associated_entities"][0]["entity_submitter_id"]
                    sample_type = entity_submitter_id.split("-")[3]
                    tumor_type = int(sample_type[0 : -1])
                    normal = 1 if tumor_type > 9 else 0

                    if normal:
                        stage_save = "normal"
                    else:
                        if "cases" in obj.keys():
                            if len(obj["cases"]):
                                if "diagnoses" in obj["cases"][0].keys():
                                    if len(obj["cases"][0]["diagnoses"]):
                                        stage = obj["cases"][0]["diagnoses"][0]["tumor_stage"]
                                        if stage != "not reported":
                                            if is_merge_stage:
                                                stage_save = tumor_stage_convert[stage.split(" ")[1]]
                                            else:
                                                stage_save = stage.split(" ")[1]
                    break
            file_name_to_stage[fpkm_file_name] = stage_save
            stage_list.append(stage_save)
        write_tab_seperated_file_for_a_list(output_path, stage_list, index_included=True)

# 所有fpkm文件的第一列的Ensembl基因列表顺序和数量均相同,要注意最后5行为统计信息,不是基因名
def generate_ensembl_gene_ids_in_a_fpkm_file():
    first_fpkm_file_path = ""
    for cancer_name in cancer_names:
        cancer_data_dir = os.path.join(rna_data_dir, cancer_name)
        filenames = os.listdir(cancer_data_dir)
        flag = False
        for filename in filenames:
            if filename.endswith(fpkm_file_end):
                flag = True
                first_fpkm_file_path = os.path.join(cancer_data_dir, filename)
                break
        if flag:
            break

    ensembl_gene_ids = read_tab_seperated_file_and_get_target_column(0, first_fpkm_file_path)
    outfile_path = os.path.join(rna_intermidiate_dir, dname, "ensembl_gene_ids.txt")
    write_tab_seperated_file_for_a_list(outfile_path, ensembl_gene_ids, index_included=False)

# 根据基因列表文件和gtf注释文件,将gene_id与gene_name进行match,之后输出到gene_symbols.txt中,该文件的索引与gene_id的索引一一对应
def obtain_gene_symbols_and_emsembl_ids_from_gtf(gene_ids_filepath, gtf_filepath, outfile_path):
    ensembl_gene_ids = read_tab_seperated_file_and_get_target_column(0, gene_ids_filepath)
    ensembl_gene_id_to_gene_symbol_dict = {}

    gtf_file = open(gtf_filepath , "r")
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
                [k , v] = item.strip().split(" ")
                group_dict[k] = v.replace('\"',"")
            if "gene_type" in group_dict.keys() and group_dict["gene_type"] == "protein_coding":
                if ("gene_id" in group_dict.keys()) and ("gene_name" in group_dict.keys()):
                    ensembl_gene_id_to_gene_symbol_dict[group_dict["gene_id"]] = group_dict["gene_name"]
        line = gtf_file.readline()
    gene_symbols = []
    for ensembl_gene_id in ensembl_gene_ids:
        if ensembl_gene_id in ensembl_gene_id_to_gene_symbol_dict.keys():
            gene_symbol = ensembl_gene_id_to_gene_symbol_dict[ensembl_gene_id]
        else:
            gene_symbol = "-"
        gene_symbols.append(gene_symbol)
    print "finished handling gtf file, now writing gene_symbols to %s" % outfile_path
    write_tab_seperated_file_for_a_list(outfile_path, gene_symbols, index_included=True)

def connect_whole_genome_to_ensembl_gene_symbol_index(ensembl_gene_symbol_index_path, out_correspondent_index_path):
    ensembl_gene_symbols = read_tab_seperated_file_and_get_target_column(1, ensembl_gene_symbol_index_path)

    #利用基因名查在全基因组中的index
    whole_gene_index_dict = {gene_name: (gindex + 1) for gindex, gene_name in enumerate(GENOME)}

    #给出全基因组中所有基因对应ensembl_gene_symbols文件中该基因的行号, -1代表没有该基因对应的gene_symbol
    ret_index_dict = {(gindex + 1): -1 for gindex, gene_name in enumerate(GENOME)}

    for gindex2, gene_symbol in enumerate(ensembl_gene_symbols):
        if (gene_symbol != "-") and (gene_symbol in alias_dict.keys()):
            gene_name = alias_dict[gene_symbol]
            if gene_name in whole_gene_index_dict.keys():
                origin_gene_idx_in_whole_genome_file = whole_gene_index_dict[gene_name]
                ret_index_dict[origin_gene_idx_in_whole_genome_file] = gindex2 + 1
    ordered_dict= sorted(ret_index_dict.items(), key=lambda d:d[0])

    with open(out_correspondent_index_path, "w") as out_correspondent_index_file:
        ltw = []
        for (origin_index, gene_symbol_index) in ordered_dict:
            ltw.append("\t".join([str(origin_index), str(gene_symbol_index)]))
        out_correspondent_index_file.write("\n".join(ltw))
    print "connect_whole_genome_to_ensembl_gene_symbol_index successful"

#计算一个fpkm文件中每个基因的tpm值
def compute_tpm_for_a_htseq_count_file(fpkm_file_path):
    fpkm_values = [float(value) for value in read_tab_seperated_file_and_get_target_column(1, fpkm_file_path)]
    # sum_of_fpkm = float(np.array(fpkm_values).sum())
    # tpm_values = [(fpkm_value / sum_of_fpkm) * 1000000.0 for fpkm_value in fpkm_values]
    return fpkm_values#tpm_values

#some global variables

gene_symbols_filepath = os.path.join(rna_intermidiate_dir, dname, "gene_symbols.txt")
emsembl_ids_filepath = os.path.join(rna_intermidiate_dir, dname, "ensembl_gene_ids.txt")
gtf_filepath = os.path.join(GRCh38_dir, "gencode.v22.annotation.gtf")

if not os.path.exists(gene_symbols_filepath):
    obtain_fpkm_filelist()
    obtain_stage_info_from_metadata()
    generate_ensembl_gene_ids_in_a_fpkm_file()
    obtain_gene_symbols_and_emsembl_ids_from_gtf(emsembl_ids_filepath, gtf_filepath, gene_symbols_filepath)

correspondent_index_path = os.path.join(rna_intermidiate_dir, dname, "correspondent_index.txt")
if not os.path.exists(correspondent_index_path):
    connect_whole_genome_to_ensembl_gene_symbol_index(gene_symbols_filepath, correspondent_index_path)


def generate_tpm_table_for_each_cancer_and_each_stage():
    correspondent_indexs = [int(value) for value in read_tab_seperated_file_and_get_target_column(1 , correspondent_index_path)]
    gene_idxs = np.array([item for item in range(len(GENOME) + 1)])
    stage_list = merged_stage if is_merge_stage else tumor_stages
    for cancer_name in cancer_names:
        cancer_data_dir = os.path.join(rna_data_dir, cancer_name)
        output_cancer_dir = os.path.join(rna_intermidiate_dir, dname, cancer_name)
        fpkm_filelist_path = os.path.join(output_cancer_dir, cancer_name + "_fpkm_filelist.txt")
        stages_path = os.path.join(output_cancer_dir, cancer_name + "_stages.txt")

        fpkm_filenames = read_tab_seperated_file_and_get_target_column(1, fpkm_filelist_path)
        fpkm_case_ids = [fpkm_filename[0 : -9] for fpkm_filename in fpkm_filenames]
        stages = read_tab_seperated_file_and_get_target_column(1, stages_path)

        stage_to_its_fpkm = { stage: [] for stage in stage_list}

        for hidx, fpkm_case_id in enumerate(fpkm_case_ids):
            stage = stages[hidx]
            stage_to_its_fpkm[stage].append(fpkm_case_id)

        for stage in stage_list:
            out_stage_tpm_data_path = os.path.join(output_cancer_dir, cancer_name + "_" + stage + "_tpm.dat")

            fpkm_case_id_list = stage_to_its_fpkm[stage]
            write_tab_seperated_file_for_a_list(os.path.join(output_cancer_dir, cancer_name + "_" + stage + "_case_ids.txt"), fpkm_case_id_list, index_included=True)

            print "cancer %s, stage %s, #cases %d" %(cancer_name, stage, len(fpkm_case_id_list))
            tpm_matrix = [gene_idxs]
            for fidx ,fpkm_case_id in enumerate(fpkm_case_id_list):
                fpkm_filepath = os.path.join(cancer_data_dir, fpkm_case_id + fpkm_file_end)
                tpm_values = compute_tpm_for_a_htseq_count_file(fpkm_filepath)
                filtered_tpm_values = [fidx + 1]
                for correspondent_index in correspondent_indexs:
                  if correspondent_index < 0:
                      filtered_tpm_values.append(-1)
                  else:
                      filtered_tpm_values.append(tpm_values[correspondent_index - 1])
                tpm_matrix.append(np.array(filtered_tpm_values))
                print "cancer_name %s\tstage %s\ttpm %d / %d" % ( cancer_name, stage, fidx + 1, len(fpkm_case_id_list))
            tpm_matrix = np.array(tpm_matrix).transpose()
            np.savetxt(out_stage_tpm_data_path, tpm_matrix, delimiter="\t")
            print "save %s stage tpm data successful!" % stage
if __name__ == '__main__':
    generate_tpm_table_for_each_cancer_and_each_stage()
    pass

