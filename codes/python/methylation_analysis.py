# -*- coding:utf-8 -*-
import os, math, re, json, time, pickle, random
import numpy as np
import matplotlib.pyplot as plt
from base import *
methy_manifest_path = os.path.join(global_files_dir, "methy_24_cancer_manifest.tsv")
methy_metadata_path = os.path.join(global_files_dir, "methy_24_cancer_meta.json")
tumor_suppressed_gene_filepath = os.path.join(global_files_dir, "gene_with_protein_product.tsv")

methy_figure_dir = os.path.join(figure_dir, "methy_scatter")

tumor_stages_xaxis = {}
for idx, item in enumerate(tumor_stages):
    tumor_stages_xaxis[item] = idx + 1

tumor_stages_xaxis2 = {}
for idx, item in enumerate(merged_stage):
    tumor_stages_xaxis2[item] = idx + 1

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

#通过json文件, 获取每个uuid对应的一级癌症阶段merged_tumor_stage和二级癌症阶段tumor_stage
def connect_uuid_to_cancer_stage(cancer_name, uuid_list, methy_metadata_path):
    stage_to_uuids = {stage_name:[] for stage_name in tumor_stages}
    merged_stage_to_uuids={stage_name:[] for stage_name in merged_stage}
    uuid_to_stage = {}

    json_obj = json.load(open(methy_metadata_path,'r'))
    cnt_cases = 0
    for uuid in uuid_list:
        filename = uuid_to_filename[uuid]
        sample_type = filename.split(".")[5].split("-")[3]
        tumor_type = int(sample_type[0 : -1])
        normal = 1 if tumor_type > 9 else 0

        stage_save = "not reported"
        if not normal:
            # if cancer, detail stage classification
            for obj in json_obj:
                if obj["file_id"] != uuid:
                    continue
                else:

                    if "cases" in obj.keys():
                        if len(obj["cases"]):
                            if "diagnoses" in obj["cases"][0].keys():
                                if len(obj["cases"][0]["diagnoses"]):
                                    stage = obj["cases"][0]["diagnoses"][0]["tumor_stage"]
                                    if stage != "not reported":
                                        stage_save = stage.split(" ")[1]
                    break
        else:
            stage_save = "normal"
        cnt_cases += 1
        uuid_to_stage[uuid] = stage_save
        stage_to_uuids[stage_save].append(uuid)
        merged_stage_to_uuids[tumor_stage_convert[stage_save]].append(uuid)
    print "cancer_name %s total cases %d" % (cancer_name, cnt_cases)
    return [uuid_to_stage, stage_to_uuids, merged_stage_to_uuids]

#最重要的一个函数, 通过遍历每个下载的tcga甲基化数据文件,将需要的统计量缓存到pkl文件中, load=True直接加载这些缓存好的文件,load=False, 从头计算并缓存, whole_genes=True代表计算全基因组的统计量, 如果计算部分基因集, 如抑癌基因集,则把它设为False
#目前缓存并输出的统计量: profile[gene][stage_idx] = [all cases' methylation values in stage_idx(int) and gene]
#                   : profile_uuid[stage_name] = [all the uuids in stage_name], 每个病人一个uuid, 每个癌症阶段对应的病人uuid不同, uuid顺序与上面每个gene对应癌症阶段的甲基化水平列表的顺序相同.
# 此处输出和缓存的数据的癌症阶段都是按照二级阶段(更细节的分期),如果需要merged_stage,请用下面的convert_origin_profile_into_merged_profile,输入为此函数输出
def gene_and_cancer_stage_profile_of_dna_methy(cancer_name, data_path, pickle_filepath, uuids, load=False, whole_genes= True):

    if not load:
        profile = {}
        profile_uuid = {}
        for tumor_stage in tumor_stages:
            profile_uuid[tumor_stage] = []

        for gene in GENOME:
            profile[gene] = []
            for ts in tumor_stages:
                profile[gene].append([])

        [uuid_to_stage, _, _] = connect_uuid_to_cancer_stage(cancer_name,uuids, methy_metadata_path)
        tot_timelapse = 0.0
        for uidx, uuid in enumerate(uuids):
            t0 = time.time()
            file_path = uuid_to_filename[uuid]

            now_file = open(data_path + file_path,'r')
            now_file.readline()
            line = now_file.readline()
            temp_gene_methy_dict = {}
            for gene_symbol in GENOME:
                temp_gene_methy_dict[gene_symbol] = []
            while line:
                line_contents = line.split("\t")
                try:
                    gene_symbols = line_contents[5].split(";")
                    positions_to_tss = line_contents[8].split(";")
                    beta_val = -1.0 if line_contents[1] == "NA" else float(line_contents[1])
                    gene_types = line_contents[6].split(";")
                    for idx, gene_symbol in enumerate(gene_symbols):
                        if gene_symbol != "." and (-2000 <= int(positions_to_tss[idx]) <= 0) and beta_val > 0.0:
                            if not whole_genes:
                                if (gene_symbol in GENOME):
                                    temp_gene_methy_dict[gene_symbol].append(beta_val)
                                    #one gene only add once for each cpg
                                    break
                            else:
                                if (gene_types[idx] == "protein_coding"):
                                    try:
                                        temp_gene_methy_dict[alias_dict[gene_symbol]].append(beta_val)
                                    except KeyError, e1:
                                        pass
                                        # print "KeyError : %s" % str(e1)
                                    #one gene only add once for each cpg
                                    break
                except IndexError, e:
                    print "line_contents :",
                    print line_contents
                line = now_file.readline()
            now_file.close()

            for gene_symbol in GENOME:
                mean_methy_level_of_this_case = float(np.array(temp_gene_methy_dict[gene_symbol]).mean())
                profile[gene_symbol][tumor_stages_xaxis[uuid_to_stage[uuid]] - 1].append(mean_methy_level_of_this_case)
            profile_uuid[uuid_to_stage[uuid]].append(uuid)

            t1 = time.time()
            ratio = float(uidx + 1.0) / float(len(uuids))
            percentage = ratio * 100.0
            time_this_turn = t1 - t0
            tot_timelapse = tot_timelapse + time_this_turn
            remain_time = (tot_timelapse / ratio) * (1.0 - ratio)
            print "%d/%d, %.2f %%, Total Time: %.2fs, Time Left: %.2fs" %(uidx + 1.0, len(uuids), percentage, tot_timelapse, remain_time)

        pickle_file = open(pickle_filepath, 'wb')
        pickle.dump(profile,pickle_file, -1)
        pickle.dump(profile_uuid,pickle_file, -1)
        pickle_file.close()
    else:
        pickle_file = open(pickle_filepath,"rb")
        profile = pickle.load(pickle_file)
        profile_uuid = pickle.load(pickle_file)
        pickle_file.close()
        print "load pickle file %s finished" % (pickle_filepath)
    return [profile, profile_uuid]

#只进行DNA甲基化数据缓存
def just_calc_methylation_pickle_pipeline():
    for cancer_name in cancer_names:
        print "now start %s" % cancer_name
        data_path = dna_methy_data_dir + os.sep+ cancer_name + os.sep
        pickle_filepath = methy_pkl_dir + os.sep + cancer_name + ".pkl"
        gene_and_cancer_stage_profile_of_dna_methy(cancer_name,data_path, pickle_filepath, uuid_dict[cancer_name], load=False, whole_genes= True)

#some global variables
gene_idx_path = os.path.join(global_files_dir, "gene_idx.txt")
with open(gene_idx_path,"w") as gene_idx_file:
    gene_idx_file.write("\n".join([str(gidx+1) + "\t" + gene for gidx, gene in enumerate(GENOME)]))

if __name__ == '__main__':
    # just_calc_methylation_pickle_pipeline()
    pass