# -*- coding:utf-8 -*-
import os, json, requests
import numpy as np
from base import *

#输入TCGA的submitter_id列表, 和查询大小(最大300), 返回一个[dict], dict的key:value分别为submitter_id : stage_name
def query_stage_of_an_submitter_id(submitter_ids, query_size):
    cases_endpt = 'https://api.gdc.cancer.gov/cases'
    filt = {"op":"and","content":[{"op":"in","content":{"field":"submitter_id","value":submitter_ids}}]}
    params = {'filters':json.dumps(filt), 'expand':'diagnoses','fields':['diagnoses.submitter_id','diagnoses.tumor_stage'],'pretty':True,'size':query_size}
    # requests URL-encodes automaticallya
    response = requests.get(cases_endpt, params = params)
    ret_obj = response.json()
    # ret_obj = json.dumps(response.json(), indent=2)
    # print ret_obj

    # print json.dumps(ret_obj["data"]["hits"],indent=2)
    submitter_id_to_stage = {}
    for tidx in range(len(submitter_ids)):
        try:
            if "hits" in ret_obj["data"].keys() and  "diagnoses" in ret_obj["data"]["hits"][tidx].keys() and len(ret_obj["data"]["hits"][tidx]["diagnoses"]):
                stage_name = ret_obj["data"]["hits"][tidx]["diagnoses"][0]["tumor_stage"]
                submitter_id = ret_obj["data"]["hits"][tidx]["diagnoses"][0]["submitter_id"]
                stage = stage_name if stage_name == "not reported" else stage_name.split(" ")[1]
                submitter_id_rev = submitter_id.split("_")[0]
                submitter_id_to_stage[submitter_id_rev] = stage
        except IndexError,e:
            print e
    print submitter_id_to_stage
    return [submitter_id_to_stage]

#从突变文件(.maf)获取所有癌症各自的submitter_id列表,并将列表写入到output_cancer_dir的cancer_name_submitter_ids.txt文件中, 文件含有两列,第一列为index, 第二列为submitter_id
def get_maf_submitter_ids(is_merge_stage):
    submitter_dict = {}
    for cancer_name in cancer_names:
        cancer_dir = os.path.join(snv_data_dir, cancer_name)
        if os.path.exists(cancer_dir):
            output_cancer_dir = os.path.join(snv_intermidiate_dir, dname, cancer_name)
            if not os.path.exists(output_cancer_dir):
                os.makedirs(output_cancer_dir)
            outfile_path = os.path.join(output_cancer_dir, cancer_name + "_submitter_ids.txt")
            submitter_dict[cancer_name] = {}
            file_names = os.listdir(cancer_dir)
            for file_name in file_names:
                if file_name.startswith("TCGA." + cancer_name + ".mutect"):
                    file_path = os.path.join(cancer_dir,file_name)
                    snv_file = open(file_path, "r")
                    print "%s start %s" % (cancer_name, file_path)
                    # pass head 6 lines
                    for i in range(6):
                        snv_file.readline()
                    line = snv_file.readline()
                    while line:
                        line_contents = line.split("\t")
                        bar_code = line_contents[15]
                        submitter_id = "-".join(bar_code.split("-")[0:3])
                        if submitter_id not in submitter_dict[cancer_name].keys():
                            submitter_dict[cancer_name][submitter_id] = 1
                        line = snv_file.readline()
                    print "end %s" % file_path

            write_tab_seperated_file_for_a_list(outfile_path, submitter_dict[cancer_name].keys(), index_included= True)

#1. 读取get_maf_submitter_ids函数产生的*_submitter_ids.txt文件, 从而获得该癌症的所有submitter_id列表;
#2. 通过query_stage_of_an_submitter_id函数, 查询query_size大小的submitter_id列表对应的癌症阶段列表(stages)
#3. *_submitter_ids.txt文件中submitter_id列表对应的癌症阶段写到cancer_name_stages.txt中,第一列为idx,对应submitter_ids.txt的idx, 第二列为癌症阶段名称(stage_name)
def get_submitter_id_stages(is_merge_stage):

    query_size = 300
    for cancer_name in cancer_names:
        cancer_dir = os.path.join(snv_data_dir, cancer_name)
        if os.path.exists(cancer_dir):
            output_cancer_dir = os.path.join(snv_intermidiate_dir, dname, cancer_name)
            input_path = os.path.join(output_cancer_dir, cancer_name + "_submitter_ids.txt")
            output_path = os.path.join(output_cancer_dir, cancer_name + "_stages.txt")
            submitter_id_stage_dict = {}
            if is_merge_stage:
                stage_and_its_submitter_dict = {stage:[] for stage in merged_stage}
            else:
                stage_and_its_submitter_dict = {stage:[] for stage in tumor_stages}
            submitter_ids = read_tab_seperated_file_and_get_target_column(1, input_path)
            print cancer_name
            len_submitters = len(submitter_ids)
            print len_submitters
            for i in range(0, len_submitters, query_size):
                sub_submitter_ids = submitter_ids[i:i+query_size]
                stages_of_subids_dict = query_stage_of_an_submitter_id(sub_submitter_ids, query_size)
                if is_merge_stage:
                    for k,v in stages_of_subids_dict[0].items():
                        submitter_id_stage_dict[k] = tumor_stage_convert[v]
                        stage_and_its_submitter_dict[tumor_stage_convert[v]].append(k)
                else:
                    for k,v in stages_of_subids_dict[0].items():
                        submitter_id_stage_dict[k] = v
                        stage_and_its_submitter_dict[v].append(k)
            out_stage_list = merged_stage_n if is_merge_stage else tumor_stages_n
            for stage in out_stage_list:
                if len(stage_and_its_submitter_dict[stage]):
                    output_stage_path = os.path.join(output_cancer_dir, cancer_name + "_" + stage +"_submitter_ids.txt")
                    write_tab_seperated_file_for_a_list(output_stage_path,stage_and_its_submitter_dict[stage], index_included=True)

            for submitter_id in submitter_ids:
                if submitter_id not in submitter_id_stage_dict.keys():
                    submitter_id_stage_dict[submitter_id] = "not reported"
            stages = [submitter_id_stage_dict[submitter_id] for submitter_id in submitter_ids]
            write_tab_seperated_file_for_a_list(output_path, stages, index_included=True)
#1. 从input_cancer_dir中, 输入cancer_name_submitter_ids.txt和cancer_name_stages.txt文件,获取该癌症所有的submitter_id列表和癌症阶段列表
#2. 建立字典对应关系,方便用submitter_id查对应癌症阶段, 即代码中的submitter_id_to_stage_dict
#3. 返回submitter_id_to_stage_dict, 和所有submitter_id列表(按从cancer_name_submitter_ids.txt中读到的顺序)
def input_submitter_id_and_its_stages(input_cancer_dir, cancer_name):
    submitter_id_to_stage_dict = {}
    submitter_id_input_filepath = os.path.join(input_cancer_dir, cancer_name + "_submitter_ids.txt")
    stages_input_filepath = os.path.join(input_cancer_dir, cancer_name + "_stages.txt")
    submitter_ids = []
    stages =[]
    with open(submitter_id_input_filepath,"r") as submitter_id_input_file:
        line = submitter_id_input_file.readline()
        while line:
            submitter_id = line.split("\t")[1]
            if submitter_id.endswith("\n"):
                submitter_id = submitter_id[0:-1]
            submitter_ids.append(submitter_id)
            line = submitter_id_input_file.readline()
    with open(stages_input_filepath,"r") as stages_input_file:
        line = stages_input_file.readline()
        while line:
            stage = line.split("\t")[1]
            if stage.endswith("\n"):
                stage = stage[0:-1]
            stages.append(stage)
            line = stages_input_file.readline()
    for sidx in range(len(submitter_ids)):
        submitter_id_to_stage_dict[submitter_ids[sidx]] = stages[sidx]
    return [submitter_id_to_stage_dict, submitter_ids]

# 运行生成突变数据的流水线: 首先生成gene_id文件,方便对数据文件的gene_id索引
# 从.maf突变数据中,统计会产生编码区mRNA突变或者翻译后氨基酸改变的基因突变数量
# 输出各个癌症阶段对应的基因突变数据, cancer_name_stage_mutation_count.dat, 第一列为gene_index, 第一行为submitter_id的index(对应见cancer_name_submitter_ids.txt文件), 内容为某基因在某样本中影响翻译的基因突变数量
def dna_mutation_data_transform_pipline(is_merge_stage):
    colum_idxs = [0, 8, 15] #要提取maf文件的列编号
    mutation_expections = ["Missense_Mutation", "Translation_Start_Site", "Splice_Region", "Splice_Site", "In_Frame_Del", "In_Frame_Ins", "Frame_Shift_Del", "Frame_Shift_Ins"]
    stage_list = merged_stage if is_merge_stage else tumor_stages
    for cancer_name in cancer_names:
        cancer_dir = os.path.join(snv_data_dir, dname, cancer_name)
        #只对有突变数据的癌症做分析
        if os.path.exists(cancer_dir):
            data_dict = {}
            output_cancer_dir = os.path.join(snv_intermidiate_dir, cancer_name)
            if os.path.exists(output_cancer_dir):
                [submitter_id_to_stage_dict, submitter_ids] = input_submitter_id_and_its_stages(output_cancer_dir, cancer_name)

                for submitter_id in submitter_ids:
                    data_dict[submitter_id] = {gene: 0 for gene in GENOME}

                file_names = os.listdir(cancer_dir)
                for file_name in file_names:
                    if file_name.startswith("TCGA." + cancer_name + ".mutect"):
                        file_path = os.path.join(cancer_dir,file_name)
                        snv_file = open(file_path, "r")
                        print "start %s" % file_path
                        # pass head 6 lines
                        for i in range(6):
                            snv_file.readline()
                        line = snv_file.readline()
                        while line:
                            line_contents = line.split("\t")
                            [gene_name, variation_classification, bar_code] = [line_contents[i] for i in colum_idxs]
                            submitter_id = "-".join(bar_code.split("-")[0:3])
                            stage = submitter_id_to_stage_dict[submitter_id]

                            if gene_name in GENOME and (variation_classification in mutation_expections):
                                data_dict[submitter_id][gene_name] += 1
                            line = snv_file.readline()
                        print "end %s" % file_path

                for cancer_stage in merged_stage_n:
                    stage_submitter_ids_path = os.path.join(output_cancer_dir, cancer_name + "_" + cancer_stage +"_submitter_ids.txt")
                    stage_submitter_ids = read_tab_seperated_file_and_get_target_column(1, stage_submitter_ids_path)
                    mut_data_filepath = os.path.join(output_cancer_dir, cancer_name + "_" + cancer_stage + "_mutation_data.dat")
                    with open(mut_data_filepath,"w") as data_file:
                        data_str = []
                        header = "\t".join([str(item) for item in range(len(stage_submitter_ids) + 1)])
                        data_str.append(header)
                        for gidx, gene in enumerate(GENOME):
                            arr = [gidx + 1]
                            arr.extend([data_dict[submitter_id][gene] for submitter_id in stage_submitter_ids])
                            data_str.append("\t".join([str(item) for item in arr]))
                        data_file.write("\n".join(data_str))

def calc_mutation_rate():
    temp_stages = merged_stage_n[0 : -1]
    for cancer_stage in temp_stages:
        for cancer_name in cancer_names:
            output_cancer_dir = os.path.join(snv_intermidiate_dir, "merged_stage", cancer_name)
            mutation_data_filepath = os.path.join(output_cancer_dir, cancer_name + "_" + cancer_stage + "_mutation_data.dat")
            out_mutation_rate_filepath = os.path.join(output_cancer_dir, cancer_name + "_" + cancer_stage + "_mutation_rate.txt")
            out_mutation_rate_sorted_filepath = os.path.join(output_cancer_dir, cancer_name + "_" + cancer_stage + "_mutation_rate_sorted.txt")
            ltws = []
            mut_dict = {}
            with open(mutation_data_filepath,"r") as data_file:
                data_file.readline()
                line = data_file.readline()
                while line:
                    line_contents = line.split("\t")
                    gene_idx = int(str(line_contents[0]))
                    gene_name = GENOME[gene_idx - 1]
                    line_data = np.sign([int(item.strip("\n")) for item in line_contents[1: -1]])
                    sum_d = line_data.sum()
                    size_d = line_data.size
                    line_mutation_rate = float(sum_d) / float(size_d)
                    ltw = "\t".join([str(gene_idx), str(round(line_mutation_rate, 4))])
                    mut_dict[gene_idx] = line_mutation_rate
                    ltws.append(ltw)
                    line = data_file.readline()
            with open(out_mutation_rate_filepath,"w") as mutation_rate_file:
                mutation_rate_file.write("\n".join(ltws))
            sorted_dict = sorted(mut_dict.items(), key=lambda d: d[1],reverse=True)
            ltws = []
            with open(out_mutation_rate_sorted_filepath,"w") as mutation_rate_sorted_file:
                for (k, v) in sorted_dict:
                    ltws.append("\t".join([str(k), str(v)]))
                mutation_rate_sorted_file.write("\n".join(ltws))
            print "finish %s" % cancer_name

is_merge_stage = False
dname = "merged_stage" if is_merge_stage else "stage"
for cancer_name in cancer_names:
    output_cancer_dir = os.path.join(snv_intermidiate_dir, dname, cancer_name)
    if not os.path.exists(output_cancer_dir):
        os.makedirs(output_cancer_dir)
    submitter_id_input_filepath = os.path.join(output_cancer_dir, cancer_name + "_submitter_ids.txt")
    first_called = False
    if not os.path.exists(submitter_id_input_filepath):
        get_maf_submitter_ids(is_merge_stage)
        first_called = True
    stages_input_filepath = os.path.join(output_cancer_dir, cancer_name + "_stages.txt")
    if not os.path.exists(stages_input_filepath):
        get_submitter_id_stages(is_merge_stage)
        first_called = True
    if first_called:
        break

if __name__ == '__main__':
    # dna_mutation_data_transform_pipline(is_merge_stage)
    calc_mutation_rate()
    #
    pass