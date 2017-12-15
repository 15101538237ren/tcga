# -*- coding:utf-8 -*-
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from base import *
methy_manifest_path = os.path.join(global_files_dir, "methy_24_cancer_manifest.tsv")
methy_metadata_path = os.path.join(global_files_dir, "methy_24_cancer_meta.json")
tumor_suppressed_gene_filepath = os.path.join(global_files_dir, "gene_with_protein_product.tsv")

tumor_stages_xaxis = {}
for idx, item in enumerate(tumor_stages):
    tumor_stages_xaxis[item] = idx + 1

tumor_stages_xaxis2 = {}
for idx, item in enumerate(merged_stage):
    tumor_stages_xaxis2[item] = idx + 1
is_merge_stage = True
stage_list = methy_and_rna_merged_stages if is_merge_stage else methy_and_rna_stages
dname = "merged_stage" if is_merge_stage else "stage"
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
        profile_cpg = {}
        profile_uuid = {}
        for tumor_stage in tumor_stages:
            profile_uuid[tumor_stage] = []

        for gene in GENOME:
            profile[gene] = []
            profile_cpg[gene] = []
            for ts in tumor_stages:
                profile[gene].append([])
                profile_cpg[gene].append([])

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
                profile_cpg[gene_symbol][tumor_stages_xaxis[uuid_to_stage[uuid]] - 1].append(temp_gene_methy_dict[gene_symbol])
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
        pickle.dump(profile_cpg,pickle_file, -1)
        pickle_file.close()
    else:
        pickle_file = open(pickle_filepath,"rb")
        profile = pickle.load(pickle_file)
        profile_uuid = pickle.load(pickle_file)
        profile_cpg = pickle.load(pickle_file)
        pickle_file.close()
        print "load pickle file %s finished" % (pickle_filepath)
    return [profile, profile_uuid, profile_cpg]

#将缓存的数据中的二级癌症阶段按照一级阶段进行合并,输出格式与gene_and_cancer_stage_profile_of_dna_methy相同
def convert_origin_profile_into_merged_profile(origin_profile_list):
    [origin_profile, origin_profile_uuid, origin_profile_cpg] = origin_profile_list
    new_profile = {gene:[] for gene in GENOME}
    new_profile_uuid = {stage:[] for stage in merged_stage}
    new_profile_cpg = {gene:[] for gene in GENOME}
    for idx, item1 in enumerate(tumor_stages):
        new_profile_uuid[tumor_stage_convert[item1]].extend(origin_profile_uuid[item1])
    for gene in GENOME:
        for stage in merged_stage:
            new_profile[gene].append([])
            new_profile_cpg[gene].append([])

        for idx, item1 in enumerate(tumor_stages):
            # print "%s\t%s\t%d\t%d" % (gene, item1, len(origin_profile[gene][idx]), len(origin_profile_uuid[item1]))
            if len(origin_profile[gene][idx]) == len(origin_profile_uuid[item1]):
                for idx2, item2 in enumerate(origin_profile[gene][idx]):
                    new_profile[gene][tumor_stages_xaxis2[tumor_stage_convert[item1]] - 1].append(item2)
                    new_profile_cpg[gene][tumor_stages_xaxis2[tumor_stage_convert[item1]] - 1].append(origin_profile_cpg[gene][idx][idx2])

    return [new_profile, new_profile_uuid, new_profile_cpg]

#将原来小类stage数据合并到大类stage
def convert_origin_to_new_profile(origin_list_of_a_gene):
    new_list = [[] for item in merged_stage]
    for idx, item1 in enumerate(tumor_stages):
        for item2 in origin_list_of_a_gene[idx]:
            new_list[tumor_stages_xaxis2[tumor_stage_convert[item1]] - 1].append(item2)
    return new_list

def plot_for_each_gene(cancer_name, gene_name, x, y, box_data, c, xrange, xticks, out_fig_path):
    plt.clf()
    plt.cla()
    fig, ax = plt.subplots()
    plt.xticks(xrange, xticks)

    ax.scatter(x, y, color=c, s=2.0)
    ax.boxplot(box_data,sym='')
    ax.set_xlim([0, len(merged_stage) - 0.5])
    ax.set_ylim([0, 1.0])
    ax.set_title(gene_name + " methylation for different cancer stage")
    plt.savefig(out_fig_path)

def merged_stage_scatter_and_box_plot(cancer_name, profile_arr, overwritten=False):
    profile = profile_arr[0]
    new_profile = {}
    figure_dir = methy_figure_dir + os.sep + cancer_name
    if not os.path.exists(figure_dir):
        os.makedirs(figure_dir)
    for gene_idx, gene in enumerate(GENOME):
        if gene in profile.keys():
            out_fig_path = methy_figure_dir + os.sep + cancer_name + os.sep + gene.lower() + '.png'
            if not overwritten:
                if os.path.exists(out_fig_path):
                    continue
            print "now plot scatter of %s %s" %(cancer_name, gene)
            xs = range(1,len(merged_stage)+1)
            new_profile[gene] = convert_origin_to_new_profile(profile[gene])

            new_x_profile = []
            new_y_profile = []
            for idx, arr in enumerate(new_profile[gene]):
                for item1 in arr:
                    ro = random.random() * 0.4 - 0.2
                    new_x_profile.append(idx + 1 + ro)
                    new_y_profile.append(item1)
                plot_for_each_gene(cancer_name, gene, new_x_profile, new_y_profile, new_profile[gene], "blue", xs, merged_stage,out_fig_path)

#put array data to .csv file, one value per line
def save_data_to_file(arr, path, precision = 4):
    file_out = open(path, "w")
    for item in arr:
        file_out.write(str(round(item, precision)) + "\n")
    file_out.close()

#保存cancer_name癌症,out_stage_list中阶段的DNA甲基化数据
def save_gene_methy_data(cancer_name, profile_list, out_stage_list, out_stage_data = False,out_xy=False, out_all_stage=False, target_gene_name=None):
    profile = profile_list[0]
    target_gene_list =[target_gene_name] if target_gene_name else GENOME
    out_methy_cancer_dir = os.path.join(methy_intermidiate_dir, dname, cancer_name)
    if not os.path.exists(out_methy_cancer_dir):
        os.makedirs(out_methy_cancer_dir)
    for gene in target_gene_list:
        if gene in profile.keys():
            gene_data = profile[gene]
            merged_data = []
            ltws_xy = []
            ltws_y = []
            for idx, stage in enumerate(merged_stage):
                if stage in out_stage_list:
                    if out_xy:
                        methy_cases_vals = gene_data[idx]
                        for item_y in methy_cases_vals:
                            ro = random.random()*0.3 - 0.15
                            x = idx + 1 + ro
                            ltw = str(round(float(x), 4)) + "," + str(round(float(item_y), 6)) + "\n"
                            ltws_xy.append(ltw)
                            tmp_stage = "n" if stage == "normal" else stage
                            ltw2 = str(round(float(item_y), 6)) + "," + tmp_stage.ljust(4) + "\n"
                            ltws_y.append(ltw2)
                    stage_data = gene_data[idx]
                    merged_data.extend(stage_data)
                    if out_stage_data:
                        save_data_to_file(stage_data, out_methy_cancer_dir + os.sep + gene + "_" + merged_stage[idx] + "_" + cancer_name + ".dat")
            if out_xy:
                out_xy_path = out_methy_cancer_dir + os.sep + gene + "_xy_" + cancer_name + ".dat"
                out_y_label_path = out_methy_cancer_dir + os.sep + gene + "_y_label_" + cancer_name + ".dat"
                out_xy_file = open(out_xy_path, "w")
                out_y_label_file = open(out_y_label_path, "w")
                out_xy_file.write("\n".join(ltws_xy))
                out_y_label_file.write("\n".join(ltws_y))
                out_xy_file.close()
                out_y_label_file.close()
            if out_all_stage:
                save_data_to_file(merged_data,  out_methy_cancer_dir + os.sep + gene + "_" + "all_stage" + "_" + cancer_name + ".dat")
    print "save methy data successfully!"

#将某癌症数据写入到tsv文件中
def dump_data_into_dat_according_to_cancer_type_and_stage(cancer_name, uuid_list, outdir, profile_list, is_merge_stage=True):
    [profile, profile_uuid,_] = profile_list

    for stage_idx, stage_name in enumerate(stage_list):
        stage_values_length = len(profile['APC'][stage_idx])
        output_cancer_dir = outdir
        if not os.path.exists(output_cancer_dir):
            os.makedirs(output_cancer_dir)
        if stage_name != "not reported" and stage_name in profile_uuid.keys() and len(profile_uuid[stage_name]):
            out_uuid_id_path = os.path.join(output_cancer_dir, cancer_name + "_" + stage_name + "_uuids.txt")
            write_tab_seperated_file_for_a_list(out_uuid_id_path, profile_uuid[stage_name],index_included=True)
            outfile_path =  os.path.join(output_cancer_dir, cancer_name + "_" + stage_name + "_methy_dat.dat")
            outfile = open(outfile_path, "w")
            out_str = []
            header = "\t".join([ str(item) for item in range(len(profile_uuid[stage_name]) + 1)])
            out_str.append(header)
            for gidx, gene in enumerate(GENOME):
                if gene in profile.keys():
                    methy_vals = [-1 for it in range(stage_values_length)]
                    for pidx, item in enumerate(profile[gene][stage_idx]):
                        item_val = round(item, 4)
                        item_str = str(item_val)
                        if item_str == "nan":
                            break
                        else:
                            methy_vals[pidx] = item_val
                    data_str = str(gidx + 1) + "\t" + "\t".join([str(item) for item in methy_vals])
                    out_str.append(data_str)
            outfile.write("\n".join(out_str))
            outfile.close()
    print "%s dump_data_into_dat_according_to_cancer_type_and_stage" % cancer_name

def print_samplesize_of_each_cancer(sample_count_filepath):
    gene_name = "APC"
    merged_stage_not_report = merged_stage[0:-1]
    ltws = []
    header_arr = ["cancer"]
    merged_stage_not_report.append("total")
    header_arr.extend(merged_stage_not_report)

    ltws.append("\t".join(header_arr))
    for cancer_name in cancer_names:
        data_path = dna_methy_data_dir + os.sep+ cancer_name + os.sep
        pickle_filepath = methy_pkl_dir + os.sep + cancer_name + ".pkl"
        temp_profile_list = gene_and_cancer_stage_profile_of_dna_methy(cancer_name,data_path, pickle_filepath, uuid_dict[cancer_name], load=True, whole_genes= True)
        new_profile_list = convert_origin_profile_into_merged_profile(temp_profile_list)
        profile = new_profile_list[0]
        arr = [cancer_name]
        tot = 0
        for sidx,stage in enumerate(merged_stage_not_report):
            stage_cnt = len(profile[gene_name][sidx])
            tot += stage_cnt
            arr.append(str(stage_cnt))
        arr.append(str(tot))
        ltws.append("\t".join(arr))

    with open(sample_count_filepath,"w") as sample_count_file:
        sample_count_file.write("\n".join(ltws))

def calc_cancer_means_and_stds_for_genome(cancer_name, cancer_profile_arr, stage_list, cancer_mean_std_dir):
    cancer_profile = cancer_profile_arr[0]
    len_stages = len(cancer_profile["APC"]) - 1 #去掉not reported
    stage_names =[stage_list[item] for item in range(len_stages)]

    out_stages_fp = os.path.join(cancer_mean_std_dir, cancer_name + "_stages.txt")
    write_tab_seperated_file_for_a_list(out_stages_fp, stage_names, index_included=True)

    mean_arr = [[] for item in GENOME]
    std_arr = [[] for item in GENOME]
    for gidx, gene in enumerate(GENOME):
        for sidx, stage_name in enumerate(stage_names):
            methy_of_this_gene = cancer_profile[gene][sidx]
            mean_arr[gidx].append(np.array(methy_of_this_gene).mean())
            std_arr[gidx].append(np.array(methy_of_this_gene).std())

    out_mean_fp = os.path.join(cancer_mean_std_dir, cancer_name + "_mean.dat")
    out_std_fp = os.path.join(cancer_mean_std_dir, cancer_name + "_std.dat")

    header = "\t".join([str(item) for item in range(len_stages + 1)])

    with open(out_mean_fp, "w") as data_file:
        ltws = [header]
        for gidx, gene in enumerate(GENOME):
            arr = [str(gidx + 1)]
            arr.extend(["-1" if str(mean_arr[gidx][sidx]) == "nan" else str(mean_arr[gidx][sidx]) for sidx in range(len_stages)])
            ltws.append("\t".join(arr))
        data_file.write("\n".join(ltws))

    with open(out_std_fp, "w") as data_file:
        ltws = [header]
        for gidx, gene in enumerate(GENOME):
            arr = [str(gidx + 1)]
            arr.extend(["-1" if str(std_arr[gidx][sidx]) == "nan" else str(std_arr[gidx][sidx]) for sidx in range(len_stages)])
            ltws.append("\t".join(arr))
        data_file.write("\n".join(ltws))

def dump_stage_std_and_mean_pipline():
    for cancer_name in cancer_names:
        data_path = dna_methy_data_dir + os.sep+ cancer_name + os.sep
        pickle_filepath = methy_pkl_dir + os.sep + cancer_name + ".pkl"
        temp_profile_list = gene_and_cancer_stage_profile_of_dna_methy(cancer_name, data_path, pickle_filepath, uuid_dict[cancer_name], load=True, whole_genes= True)
        new_profile_list = convert_origin_profile_into_merged_profile(temp_profile_list)
        profile_list = new_profile_list if is_merge_stage else temp_profile_list
        cancer_mean_std_dir = os.path.join(methy_mean_std_dir, dname, cancer_name)
        if not os.path.exists(cancer_mean_std_dir):
            os.makedirs(cancer_mean_std_dir)
        calc_cancer_means_and_stds_for_genome(cancer_name, profile_list, stage_list, cancer_mean_std_dir)

#生成dna甲基化的dat文件
def dump_data_into_dat_pipepile():
    for cancer_name in cancer_names:
        print "now start %s" % cancer_name
        data_path = dna_methy_data_dir + os.sep+ cancer_name + os.sep
        pickle_filepath = methy_pkl_dir + os.sep + cancer_name + ".pkl"
        temp_profile_list = gene_and_cancer_stage_profile_of_dna_methy(cancer_name, data_path, pickle_filepath, uuid_dict[cancer_name], load=True, whole_genes= True)
        new_profile_list = convert_origin_profile_into_merged_profile(temp_profile_list)
        profile_list = new_profile_list if is_merge_stage else temp_profile_list

        out_dir = os.path.join(methy_intermidiate_dir, dname, cancer_name)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        dump_data_into_dat_according_to_cancer_type_and_stage(cancer_name, uuid_dict[cancer_name], out_dir, profile_list, is_merge_stage=is_merge_stage)

# 保存甲基化数据的pipline
def save_gene_methy_data_pipeline():
    out_stage_list = ["normal","i","ii","iii","iv"]
    target_gene_name = "EYA4"
    for cancer_name in cancer_names:
        if cancer_name == "COAD":
            print "now start %s" % cancer_name
            data_path = dna_methy_data_dir + os.sep+ cancer_name + os.sep
            pickle_filepath = methy_pkl_dir + os.sep + cancer_name + ".pkl"
            temp_profile_list = gene_and_cancer_stage_profile_of_dna_methy(cancer_name,data_path, pickle_filepath, uuid_dict[cancer_name], load=True, whole_genes= True)
            new_profile_list = convert_origin_profile_into_merged_profile(temp_profile_list)
            save_gene_methy_data(cancer_name, new_profile_list, out_stage_list, out_stage_data=False, out_xy=False, out_all_stage=True, target_gene_name=target_gene_name)

#只进行DNA甲基化数据缓存
def just_calc_methylation_pickle_pipeline():
    for cancer_name in cancer_names:
        print "now start %s" % cancer_name
        data_path = dna_methy_data_dir + os.sep+ cancer_name + os.sep
        pickle_filepath = methy_pkl_dir + os.sep + cancer_name + ".pkl"
        gene_and_cancer_stage_profile_of_dna_methy(cancer_name,data_path, pickle_filepath, uuid_dict[cancer_name], load=False, whole_genes= True)

def dump_entropy_into_dat_according_to_cancer_type_and_stage(cancer_name, in_dir, out_dir, bins):
    out_entropy_dat_fp = os.path.join(out_dir, cancer_name + "_entropy.dat")
    entropy_matrix = [[item for item in range(len(stage_list) + 1)]]
    for gene_idx, gene in enumerate(GENOME):
        default_entropy_values = [gene_idx + 1]
        default_entropy_values.extend([-1 for item in stage_list])
        entropy_matrix.append(default_entropy_values)

    for stage_idx, stage_name in enumerate(stage_list):
        input_dat_fp = os.path.join(in_dir, cancer_name + "_" + stage_name + "_methy_dat.dat")
        methy_dat = pd.read_csv(input_dat_fp, sep='\t', lineterminator='\n', header= 0, index_col=0, dtype=np.float64)
        methy_matrix = methy_dat.values

        stage_sample_counts = len(methy_matrix[0])
        for gene_idx, methy_values_of_gene_i in enumerate(methy_matrix):
            #若beta_value矩阵的gene_idx行中没有-1的元素,才计算该基因的信息熵
            if not any(methy_value < 0 for methy_value in methy_values_of_gene_i):
                hist, bin_edges = np.histogram(methy_values_of_gene_i, bins=bins)
                freqencies = hist / float(stage_sample_counts)
                positive_frequencies = [freq for freq in freqencies if freq > 10e-6]
                entropy_of_this_gene = 0.0
                for fidx, freq in enumerate(positive_frequencies):
                    entropy_of_this_gene += -1.0 * freq * math.log(freq)
                # print entropy_of_this_gene
                entropy_matrix[gene_idx + 1][stage_idx + 1] = entropy_of_this_gene
    np.savetxt(out_entropy_dat_fp, np.array(entropy_matrix), delimiter="\t")
    print "save %s successful!" % out_entropy_dat_fp

#根据各癌症各阶段的.dat文件(基因的DNA甲基化水平矩阵),计算并生成entropy矩阵
def dump_entropy_into_dat_pipeline():
    bins = np.linspace(0.0, 1.0, num=51)
    out_dir = os.path.join(methy_entropy_dir, dname)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    stages_fp = os.path.join(methy_entropy_dir, dname, "stage_idx.txt")
    write_tab_seperated_file_for_a_list(stages_fp,stage_list, index_included=True)
    for cancer_name in cancer_names:
        print "now start %s" % cancer_name
        in_dir = os.path.join(methy_intermidiate_dir, dname, cancer_name)
        dump_entropy_into_dat_according_to_cancer_type_and_stage(cancer_name, in_dir, out_dir, bins)
def calc_methy_correlation(cancer_name, methy_in_dir, gidx_in_dir, out_dir , stage_wanted, middle_name):
    input_dat_fp = os.path.join(methy_in_dir, cancer_name + "_" + stage_wanted + "_methy_dat.dat")
    gidx_fp_path = os.path.join(gidx_in_dir, cancer_name + "_" + stage_wanted + "_" + middle_name + "_gidx.txt")
    gidxs =[int(gidx) for gidx in read_tab_seperated_file_and_get_target_column(1, gidx_fp_path)]
    len_gidx = len(gidxs)
    out_corr_dat_fp = os.path.join(out_dir, cancer_name + "_" + stage_wanted + "_" + middle_name + "_corr.dat")
    methy_dat = pd.read_csv(input_dat_fp, sep='\t', lineterminator='\n', header= 0, index_col=0, dtype=np.float64)
    methy_matrix = methy_dat.values
    corr_matrix = np.ones((len_gidx + 1, len_gidx + 1)) * (-10)
    corr_matrix[0][0] = 0
    for mi in range(len_gidx):
        corr_matrix[mi + 1][0] = mi + 1
        corr_matrix[0][mi + 1] = mi + 1

        for mj in range(len_gidx):
            if mi == mj:
                corr_matrix[mi + 1][mj + 1] = 1.0
            elif mj < mi:
                i_arr = methy_matrix[gidxs[mi] - 1]
                if any(methy_value < 0 for methy_value in i_arr):
                    continue

                j_arr = methy_matrix[gidxs[mj] - 1]
                if (mi != mj) and any(methy_value < 0 for methy_value in j_arr):
                    continue
                corr_val = np.corrcoef(i_arr, j_arr)[0, 1]
                corr_matrix[mi + 1][mj + 1] = corr_val
                corr_matrix[mj + 1][mi + 1] = corr_val
        print mi
    np.savetxt(out_corr_dat_fp, corr_matrix, delimiter="\t")
    print "save %s successful!" % out_corr_dat_fp
def calc_methy_correlation_pipeline():
    cancer_name = "COAD"
    middle_name_list = ["pp","pn"]
    for middle_name in middle_name_list:
        out_dir = os.path.join(methy_corr_dir, dname, cancer_name)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        stage_wanted = "i"
        print "now start %s" % cancer_name
        methy_in_dir = os.path.join(methy_intermidiate_dir, dname, cancer_name)
        gidx_in_dir = out_dir
        calc_methy_correlation(cancer_name, methy_in_dir, gidx_in_dir, out_dir, stage_wanted, middle_name)
sample_count_path = os.path.join(global_files_dir, "sample_count.txt")
if not os.path.exists(sample_count_path):
    print_samplesize_of_each_cancer(sample_count_path)

if __name__ == '__main__':
    # dump_entropy_into_dat_pipeline()
    calc_methy_correlation_pipeline()
    pass