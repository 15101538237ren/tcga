# -*- coding:utf-8 -*-
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import numpy.ma as ma
from base import *
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

#通过json文件, 获取每个uuid对应的一级癌症阶段merged_tumor_stage和二级癌症阶段tumor_stage
def connect_uuid_to_cancer_stage_and_age(cancer_name, uuid_list, methy_metadata_path):
    stage_to_uuids = {stage_name:[] for stage_name in tumor_stages}
    merged_stage_to_uuids={stage_name:[] for stage_name in merged_stage}
    uuid_to_stage = {}
    uuid_to_age = {}

    json_obj = json.load(open(methy_metadata_path,'r'))
    cnt_cases = 0
    for uuid in uuid_list:
        filename = uuid_to_filename[uuid]
        sample_type = filename.split(".")[5].split("-")[3]
        tumor_type = int(sample_type[0 : -1])
        normal = 1 if tumor_type > 9 else 0

        stage_save = "not reported"
        age_save = -1
        # if cancer, detail stage classification
        for obj in json_obj:
            if obj["file_id"] != uuid:
                continue
            else:

                if "cases" in obj.keys():
                    if len(obj["cases"]):
                        if "diagnoses" in obj["cases"][0].keys():
                            if len(obj["cases"][0]["diagnoses"]):
                                tmp_age = obj["cases"][0]["diagnoses"][0]["age_at_diagnosis"]
                                age_save = tmp_age / 365.25 if tmp_age else -1

                                if not normal:
                                    stage = obj["cases"][0]["diagnoses"][0]["tumor_stage"]
                                    if stage != "not reported":
                                        stage_save = stage.split(" ")[1]
                                else:
                                    stage_save = "normal"
                break
        cnt_cases += 1
        uuid_to_stage[uuid] = stage_save
        uuid_to_age[uuid] = age_save
        stage_to_uuids[stage_save].append(uuid)
        merged_stage_to_uuids[tumor_stage_convert[stage_save]].append(uuid)
    print "cancer_name %s total cases %d" % (cancer_name, cnt_cases)
    return [uuid_to_stage, uuid_to_age, stage_to_uuids, merged_stage_to_uuids]

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

        [uuid_to_stage,_ , _, _] = connect_uuid_to_cancer_stage_and_age(cancer_name, uuids, methy_metadata_path)
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
                    beta_val = -1.0 if line_contents[1] == "NA" else float(line_contents[1])
                    gene_types = line_contents[6].split(";")
                    for idx, gene_symbol in enumerate(gene_symbols):
                        if (gene_symbol != ".") and (gene_types[idx] == "protein_coding") and (beta_val > 0.0):
                            if not whole_genes:
                                if (gene_symbol in GENOME):
                                    temp_gene_methy_dict[gene_symbol].append(beta_val)
                                    #one gene only add once for each cpg
                                    break
                            else:
                                try:
                                    gname = alias_dict[gene_symbol]
                                    cpg_start = int(line_contents[3])
                                    g_start = gene_infos[gname]["start"]
                                    g_end = gene_infos[gname]["end"]
                                    g_strand = gene_infos[gname]["strand"]
                                    pttss = cpg_start - g_start if g_strand else g_end - cpg_start
                                    # 启动子
                                    if - promoter_length < pttss < 0:
                                        temp_gene_methy_dict[gname].append(beta_val)
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
def save_gene_methy_data(cancer_name, profile_list, out_stage_list, out_stage_data = False,out_xy=False, out_all_stage=False, target_gene_list=None):
    profile = profile_list[0]
    target_gene_list = target_gene_list if target_gene_list else GENOME
    out_methy_cancer_dir = os.path.join(methy_matlab_data_dir, dname, cancer_name)
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
        stage_name_rep = stage_name.replace(" ", "_")
        output_cancer_dir = outdir
        if not os.path.exists(output_cancer_dir):
            os.makedirs(output_cancer_dir)
        if stage_name in profile_uuid.keys() and len(profile_uuid[stage_name]):
            out_uuid_id_path = os.path.join(output_cancer_dir, cancer_name + "_" + stage_name_rep + "_uuids.txt")
            write_tab_seperated_file_for_a_list(out_uuid_id_path, profile_uuid[stage_name],index_included=True)
            outfile_path =  os.path.join(output_cancer_dir, cancer_name + "_" + stage_name_rep + "_methy_dat.dat")
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

    with open(sample_count_filepath, "w") as sample_count_file:
        sample_count_file.write("\n".join(ltws))

def calc_cancer_means_and_stds_for_genome(cancer_name, cancer_profile_arr, stage_list, cancer_mean_std_dir):
    cancer_profile = cancer_profile_arr[0]
    stage_names = stage_list
    len_stages = len(stage_list)

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
        if cancer_name == "COAD":
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
        if cancer_name == "COAD":
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
    vogelstein_genes = [ gene for gidx, gene in enumerate(GENOME) if gene_categorys_vogelstein[gidx] > 0]
    target_gene_list = vogelstein_genes
    for cancer_name in cancer_names:
        print "now start %s" % cancer_name
        data_path = dna_methy_data_dir + os.sep+ cancer_name + os.sep
        pickle_filepath = methy_pkl_dir + os.sep + cancer_name + ".pkl"
        temp_profile_list = gene_and_cancer_stage_profile_of_dna_methy(cancer_name,data_path, pickle_filepath, uuid_dict[cancer_name], load=True, whole_genes= True)
        new_profile_list = convert_origin_profile_into_merged_profile(temp_profile_list)
        save_gene_methy_data(cancer_name, new_profile_list, out_stage_list, out_stage_data=True, out_xy=True, out_all_stage=True,target_gene_list=target_gene_list)

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

def sort_pscore_pipline():
    mid_names = ["p", "n"]
    for cancer_name in cancer_names:
        for mid_name in mid_names:
            print "now %s %s" %(cancer_name, mid_name)
            pscore_dict = {}
            pscore_fp = os.path.join(methy_pvalue_dir, dname, cancer_name, cancer_name + "_" + mid_name + "_score.dat")
            with open(pscore_fp, "r") as pscore_file:
                line = pscore_file.readline()
                while line:
                    line_contents = line.strip("\n").split("\t")
                    gidx = int(line_contents[0])
                    pscore = float(line_contents[-1])
                    pscore_dict[gidx] = pscore
                    line = pscore_file.readline()

            sorted_pscore = sorted(pscore_dict.items(), key=lambda d: d[1], reverse=True)
            out_sorted_pscore_all_file = open(os.path.join(methy_pvalue_dir, dname, cancer_name, cancer_name + "_" + mid_name + "_score_sorted.txt"),"w")
            out_sorted_pscore_onco_file = open(os.path.join(methy_pvalue_dir, dname, cancer_name, cancer_name + "_" + mid_name + "_score_onco.txt"),"w")
            out_sorted_pscore_tsg_file = open(os.path.join(methy_pvalue_dir, dname, cancer_name, cancer_name + "_" + mid_name + "_score_tsg.txt"),"w")
            out_sorted_pscore_both_file = open(os.path.join(methy_pvalue_dir, dname, cancer_name, cancer_name + "_" + mid_name + "_score_both.txt"),"w")
            out_sorted_pscore_other_file = open(os.path.join(methy_pvalue_dir, dname, cancer_name, cancer_name + "_" + mid_name + "_score_other.txt"),"w")

            for k, v in sorted_pscore:
                ltw = "\t".join([str(k), GENOME[k - 1], str(v)]) + "\n"
                if gene_categorys[k - 1] == 1:
                    out_sorted_pscore_onco_file.write(ltw)
                elif gene_categorys[k - 1] == 2:
                    out_sorted_pscore_tsg_file.write(ltw)
                elif gene_categorys[k - 1] == 3:
                    out_sorted_pscore_both_file.write(ltw)
                else:
                    out_sorted_pscore_other_file.write(ltw)
                out_sorted_pscore_all_file.write(ltw)
            print "write sorted files successful"

            out_sorted_pscore_all_file.close()
            out_sorted_pscore_onco_file.close()
            out_sorted_pscore_tsg_file.close()
            out_sorted_pscore_both_file.close()
            out_sorted_pscore_other_file.close()

def dump_sample_entropy_into_dat_according_to_cancer_type_and_stage(cancer_name, in_dir, out_dir, nbins):
    entropy_array = [[] for item in range(len(stage_list))]

    normal_methy_dat_fp = os.path.join(in_dir, cancer_name + "_normal_methy_dat.dat")
    normal_methy_dat = pd.read_csv(normal_methy_dat_fp, sep='\t', lineterminator='\n', header=0, index_col=0, dtype=np.float64)
    normal_ma_mean = ma.masked_less(normal_methy_dat.values, 0).mean(axis= 1)
    normal_ma_mean[normal_ma_mean.mask] = 0.0
    matlab_plot_entropy_array = []

    for stage_idx, stage_name in enumerate(stage_list):
        #创建初始的病人熵数组
        sample_input_fp = os.path.join(in_dir, cancer_name + "_" + stage_name + "_uuids.txt")
        sample_uuids = read_tab_seperated_file_and_get_target_column(1, sample_input_fp)
        for sample_uuid in sample_uuids:
            entropy_array[stage_idx].append(-1)

        input_dat_fp = os.path.join(in_dir, cancer_name + "_" + stage_name + "_methy_dat.dat")
        methy_dat = pd.read_csv(input_dat_fp, sep='\t', lineterminator='\n', header= 0, index_col=0, dtype=np.float64)

        methy_ma_matrix = ma.masked_less(np.transpose(methy_dat.values), 0)
        for sample_idx, sample_methy_array in enumerate(methy_ma_matrix):
            compressed_methy_arr = (sample_methy_array - normal_ma_mean).compressed() #sample_methy_array.compressed()
            hist, bin_edges = np.histogram(compressed_methy_arr, bins=nbins, range=(compressed_methy_arr.min(), compressed_methy_arr.max()))
            freqencies = hist / float(len(compressed_methy_arr))
            positive_frequencies = [freq for freq in freqencies if freq > 10e-9]
            entropy_of_this_sample = 0.0
            for fidx, freq in enumerate(positive_frequencies):
                entropy_of_this_sample += -1.0 * freq * math.log(freq)
            # print entropy_of_this_gene
            entropy_array[stage_idx][sample_idx] = entropy_of_this_sample
            ro = random.random() * 0.4 - 0.2
            matlab_plot_entropy_array.append("\t".join([str(stage_idx + 1 + ro), str(entropy_of_this_sample)]))
        out_entropy_dat_fp = os.path.join(out_dir, cancer_name + "_" + stage_name + "_sample_entropy.dat")
        write_tab_seperated_file_for_a_list(out_entropy_dat_fp, entropy_array[stage_idx], index_included=True)

    out_matlab_plot_entropy_fp = os.path.join(out_dir, cancer_name + "_matlab_sample_entropy.dat")
    with open(out_matlab_plot_entropy_fp,"w") as matlab_plot_entropy_file:
        for mat_item in matlab_plot_entropy_array:
            matlab_plot_entropy_file.write(mat_item + "\n")
        print "write %s successful!" % out_matlab_plot_entropy_fp
#计算每个病人的熵
def dump_sample_entropy_into_dat_pipeline():
    nbins = 50
    out_dir = os.path.join(methy_sample_entropy_dir, dname)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    stages_fp = os.path.join(methy_sample_entropy_dir, dname, "stage_idx.txt")
    write_tab_seperated_file_for_a_list(stages_fp, stage_list, index_included=True)

    for cancer_name in cancer_names:
        print "now start %s" % cancer_name
        in_dir = os.path.join(methy_intermidiate_dir, dname, cancer_name)
        dump_sample_entropy_into_dat_according_to_cancer_type_and_stage(cancer_name,in_dir, out_dir, nbins)

def generate_patient_age_data_pipline():
    for cancer_name in cancer_names:
        [_, uuid_to_age, _, _] = connect_uuid_to_cancer_stage_and_age(cancer_name, uuid_dict[cancer_name], methy_metadata_path)
        cancer_dir = os.path.join(methy_intermidiate_dir, dname, cancer_name)

        for stage_idx, stage_name in enumerate(stage_list):
            # 创建初始的病人熵数组
            sample_input_fp = os.path.join(cancer_dir, cancer_name + "_" + stage_name + "_uuids.txt")
            sample_uuids = read_tab_seperated_file_and_get_target_column(1, sample_input_fp)
            sample_ages = [uuid_to_age[sample_uuid] for sample_uuid in sample_uuids]

            age_output_fp = os.path.join(cancer_dir, cancer_name + "_" + stage_name + "_ages.txt")
            write_tab_seperated_file_for_a_list(age_output_fp, sample_ages,index_included=True)

sample_count_path = os.path.join(global_files_dir, "sample_count.txt")
if not os.path.exists(sample_count_path):
    pass
    # print_samplesize_of_each_cancer(sample_count_path)

if __name__ == '__main__':
    # just_calc_methylation_pickle_pipeline()
    # dump_data_into_dat_pipepile()
    # save_gene_methy_data_pipeline()
    # dump_entropy_into_dat_pipeline()
    # dump_stage_std_and_mean_pipline()
    # calc_methy_correlation_pipeline()
    # sort_pscore_pipline()
    # dump_sample_entropy_into_dat_pipeline()
    generate_patient_age_data_pipline()
    pass