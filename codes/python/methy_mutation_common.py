# -*- coding:utf-8 -*-
from base import *
import pandas as pd
import numpy as np
import time
is_merge_stage = True
common_stages = mutation_merged_stage[0 : -1] if is_merge_stage else mutation_stage[0 : -1]
dname = "merged_stage" if is_merge_stage else "stage"
methy_metadata_path = os.path.join(global_files_dir, "methy_24_cancer_meta.json")
methy_metadata_json_obj = json.load(open(methy_metadata_path, 'r'))
common_cases_path = os.path.join(global_files_dir, "common_cases_of_methy_and_mutation.txt")
common_patient_file_name = 'common_patients_idxs.txt'
#输入TCGA的uuids列表, 和查询大小(最大300), 返回一个list,为对应uuids的submitter_ids
def query_submitter_id_of_a_uuid(uuids):
    submitter_ids = []
    for uuid in uuids:
        for obj in methy_metadata_json_obj:
            if obj["file_id"] != uuid:
                continue
            else:
                submitter_id_full = obj["associated_entities"][0]["entity_submitter_id"]
                submitter_id = "-".join(submitter_id_full.split("-")[0 : 3])
                submitter_ids.append(submitter_id)
    return submitter_ids

def extract_submitter_ids_from_methylation_uuids_and_mutation_submitter_ids():
    ltws = ["\t".join(['cancer', '#methy (cases)', '#mutation', '#total common', '#i', '#ii', '#iii', '#iv'])]

    for cancer_name in cancer_names:
        methy_tot = 0
        mut_tot = 0
        common_tot = 0
        ltw_of_stage_common = []
        for cancer_stage in common_stages:
            cancer_stage_rep = cancer_stage.replace(" ", "_")
            out_idx_dir = os.path.join(common_patient_data_dir, cancer_name, cancer_stage_rep)
            if not os.path.exists(out_idx_dir):
                os.makedirs(out_idx_dir)
            mutation_submitter_ids_fp = os.path.join(snv_intermidiate_dir, dname, cancer_name, cancer_name + "_" + cancer_stage_rep + "_submitter_ids.txt")
            mutation_submitter_ids = read_tab_seperated_file_and_get_target_column(1, mutation_submitter_ids_fp)
            mutation_idx_dict = {mitem: (midx + 1) for midx, mitem in enumerate(mutation_submitter_ids)}
            mut_tot += len(mutation_submitter_ids)

            methy_uuids_fp = os.path.join(methy_intermidiate_dir, dname, cancer_name, cancer_name + "_" + cancer_stage_rep + "_uuids.txt")
            methy_uuids = read_tab_seperated_file_and_get_target_column(1, methy_uuids_fp)

            methy_submitter_ids = query_submitter_id_of_a_uuid(methy_uuids)
            methy_file_name_dict = {submitter_id : uuid_to_filename[methy_uuids[sidx]]  for sidx, submitter_id in enumerate(methy_submitter_ids)}

            methy_idx_dict = {uuid_to_filename[methy_item]: (methy_idx + 1) for methy_idx, methy_item in enumerate(methy_uuids)}

            methy_tot += len(methy_uuids)
            if len(methy_uuids) == len(methy_submitter_ids):
                common_submitter_ids = list(set(methy_submitter_ids).intersection(set(mutation_submitter_ids)))
                common_tot += len(common_submitter_ids)
                common_methy_file_names = [methy_file_name_dict[item] for item in common_submitter_ids]
                ltw_of_stage_common.append(str(len(common_submitter_ids)))

                mut_idxs = [mutation_idx_dict[csid] for csid in common_submitter_ids]

                write_tab_seperated_file_for_a_list(os.path.join(out_idx_dir, 'common_patients_mutation_idxs.txt'), mut_idxs, index_included=True)
                methy_idxs = [methy_idx_dict[methy_file_name_dict[csid]] for csid in common_submitter_ids]

                write_tab_seperated_file_for_a_list(os.path.join(out_idx_dir, 'common_patients_methy_idxs.txt'), methy_idxs, index_included=True)

                out_idx_fp = os.path.join(out_idx_dir, common_patient_file_name)
                with open(out_idx_fp, "w") as out_idx_file:
                    tltws = []
                    for common_idx in range(len(common_submitter_ids)):
                        cidx = str(common_idx + 1)
                        c_submitter_id =  common_submitter_ids[common_idx]
                        c_methy_filename = common_methy_file_names[common_idx]
                        tltws.append("\t".join([cidx, c_submitter_id, c_methy_filename]))
                    out_idx_file.write('\n'.join(tltws))
        ltw_arr = [cancer_name, str(methy_tot), str(mut_tot), str(common_tot)]
        ltw_arr.extend(ltw_of_stage_common)
        ltw = "\t".join(ltw_arr)
        ltws.append(ltw)

    with open(common_cases_path,"w") as common_cases_file:
        common_cases_file.write("\n".join(ltws))
        print "write %s successful" % common_cases_path
    return 0

def obtain_promoter_and_genebody_methy_status():
    gene_infos = {}
    for gidx, gene_name in enumerate(GENOME):
        chr_no, start, end, strand = gene_pos_labels_used[gidx]
        gene_infos[gene_name] = {'chr': chr_no, 'start': start, 'end': end, 'strand': strand}
    tot_time = 0.0
    for cancer_name in cancer_names:
        for cancer_stage in common_stages:
            cancer_stage_rep = cancer_stage.replace(" ", "_")
            out_idx_dir = os.path.join(common_patient_data_dir, cancer_name, cancer_stage_rep)
            out_idx_fp = os.path.join(out_idx_dir, common_patient_file_name)
            common_submitter_ids = read_tab_seperated_file_and_get_target_column(1, out_idx_fp)
            common_filenames = read_tab_seperated_file_and_get_target_column(2, out_idx_fp)

            for sidx, csid in enumerate(common_submitter_ids):
                print "%s %s %d of %d" % (cancer_name, cancer_stage, sidx + 1, len(common_submitter_ids))
                t0 = time.time()
                out_methy_fp = os.path.join(out_idx_dir, str(sidx + 1) + '_methy.tsv')
                methy_dict = {gene_name: {} for gene_name in GENOME}
                methy_fp = os.path.join(dna_methy_data_dir, cancer_name, common_filenames[sidx])
                with open(methy_fp, "r") as methy_file:
                    #因为计算每个基因在promoter(通过distance_to_tss判断)和gene_body(通过cpg位置判断), CpG位点的甲基化水平
                    methy_file.readline()
                    line = methy_file.readline()
                    while line:
                        try:
                            line_contents = line.split("\t")
                            gene_symbols = line_contents[5].split(";")
                            beta_val = -1.0 if line_contents[1] == "NA" else float(line_contents[1])
                            gene_types = line_contents[6].split(";")
                            for idx, gene_symbol in enumerate(gene_symbols):
                                if gene_symbol != "." and gene_types[idx] == "protein_coding" and beta_val > 0.0:
                                    cpg_start = int(line_contents[3])
                                    g_start = gene_infos[gene_symbol]["start"]
                                    g_end = gene_infos[gene_symbol]["end"]
                                    g_strand = gene_infos[gene_symbol]["strand"]
                                    pttss = cpg_start - g_start if g_strand else g_end - cpg_start
                                    #要么是启动子，要么在gene body
                                    if (- promoter_length < pttss < 0) or (g_start <= cpg_start <= g_end):
                                        methy_dict[gene_symbol][pttss] = ",".join([str(pttss), str(cpg_start), str(beta_val)])
                                        break
                        except KeyError, e1:
                            pass
                        line = methy_file.readline()
                with open(out_methy_fp,"w") as out_methy_file:
                    ltws = []
                    for gidx, gene_name in enumerate(GENOME):
                        if methy_dict[gene_name]:
                            sorted_dict = methy_dict[gene_name]
                            sorted_dict = sorted(sorted_dict.items(), key=lambda d: d[0])
                            ltw_t = []
                            for k, v in sorted_dict:
                                ltw_t.append(v)
                            ltw = str(gidx + 1) + "\t" + ";".join(ltw_t)
                        else:
                            ltw = str(gidx + 1)
                        ltws.append(ltw)
                    out_methy_file.write("\n".join(ltws))
                    print "write %s successful" % out_methy_fp
                t1 = time.time()
                t_used_time = t1 - t0
                tot_time += t_used_time
                remain_time = (tot_time / (sidx + 1.0)) * (len(common_submitter_ids) - sidx - 1.0)
                print "%d of %d, %.2f%%, Total Time: %.2f, Time Left: %.2f" % (sidx + 1, len(common_submitter_ids), (sidx + 1.0)/len(common_submitter_ids), tot_time, remain_time)

def obtain_promoter_and_genebody_mutation_status():
    gene_infos = {}
    for gidx, gene_name in enumerate(GENOME):
        chr_no, start, end, strand = gene_pos_labels_used[gidx]
        gene_infos[gene_name] = {'chr': chr_no, 'start': start, 'end': end, 'strand': strand}

    tot_time = 0.0

    for cid, cancer_name in enumerate(cancer_names):
        print "%s %d of %d" % (cancer_name, cid + 1, len(cancer_names))
        t0 = time.time()

        submitter_id_to_stage = {}
        SNP_INS_DEL_dict = [{}, {}, {}]

        for cancer_stage in common_stages:
            cancer_stage_rep = cancer_stage.replace(" ", "_")
            for iclass in range(len(SNP_Ins_Del_classification.keys())):
                SNP_INS_DEL_dict[iclass][cancer_stage_rep] = {}

            out_idx_dir = os.path.join(common_patient_data_dir, cancer_name, cancer_stage_rep)
            out_idx_fp = os.path.join(out_idx_dir, common_patient_file_name)
            common_submitter_ids = read_tab_seperated_file_and_get_target_column(1, out_idx_fp)

            for common_submitter_id in common_submitter_ids:
                submitter_id_to_stage[common_submitter_id] = cancer_stage_rep

                for iclass in range(len(SNP_Ins_Del_classification.keys())):
                    SNP_INS_DEL_dict[iclass][cancer_stage_rep][common_submitter_id] = { gene_name:{} for gene_name in GENOME}

        cancer_mut_data_dir = os.path.join(snv_data_dir, cancer_name)
        file_names = os.listdir(cancer_mut_data_dir)
        for file_name in file_names:
            if file_name.startswith("TCGA." + cancer_name + ".mutect"):
                file_path = os.path.join(cancer_mut_data_dir, file_name)
                with open(file_path, "r") as snv_file:
                    print "start %s" % file_path
                    # pass head 6 lines
                    for i in range(6):
                        snv_file.readline()
                    line = snv_file.readline()
                    while line:
                        try:
                            line_contents = line.split("\t")
                            bar_code = line_contents[15]
                            submitter_id = "-".join(bar_code.split("-")[0:3])
                            #common submitter id
                            if submitter_id in submitter_id_to_stage.keys():
                                variant_type = line_contents[9]
                                if variant_type in SNP_Ins_Del_classification.keys():
                                    gname = line_contents[0]
                                    g_start = gene_infos[gname]["start"]
                                    g_end = gene_infos[gname]["end"]
                                    g_strand = gene_infos[gname]["strand"]

                                    s_stage = submitter_id_to_stage[submitter_id]
                                    variant_class = line_contents[8]

                                    v_start = int(line_contents[5])
                                    v_end = int(line_contents[6])

                                    if variant_type == "SNP":
                                        pttss = v_start - g_start if g_strand else g_end - v_start
                                        HGVScs = line_contents[34].split(">")
                                        if len(HGVScs) >= 2:
                                            change = HGVScs[0][-1] + "->" + HGVScs[1]
                                        else:
                                            change = '-'
                                        v_context = line_contents[111]
                                        if (- promoter_length < pttss < 0) or (g_start <= v_start <= g_end):
                                            SNP_INS_DEL_dict[SNP_Ins_Del_classification[variant_type]][s_stage][
                                                submitter_id][gname][pttss] = {'rel_start':pttss, 'v_start':v_start,
                                                                               'class_code': mutation_classification[variant_class],
                                                                               'class':variant_class, 'context':v_context, 'change': change}
                                    else:
                                        p_start = g_start - promoter_length if g_strand else g_start
                                        p_end = g_end if g_strand else g_end + promoter_length
                                        reference_allele = line_contents[10] if variant_type == "DEL" else line_contents[12]
                                        if (p_start <= v_start) and (v_end <= p_end):
                                            s_to_tss = v_start - g_start if g_strand else g_end - v_start
                                            e_to_tss = v_end - g_start if g_strand else g_end - v_end
                                            SNP_INS_DEL_dict[SNP_Ins_Del_classification[variant_type]][s_stage][
                                                submitter_id][gname][s_to_tss] = {'rel_start': s_to_tss, 'v_start':v_start,
                                                                                  'rel_end' : e_to_tss, 'v_end' : v_end,
                                                                                   'class_code': mutation_classification[variant_class],
                                                                                  'class': variant_class, 'change': reference_allele}
                        except KeyError,e1:
                            pass
                        line = snv_file.readline()
                    print "end %s" % file_path
        for cancer_stage in common_stages:
            cancer_stage_rep = cancer_stage.replace(" ", "_")
            out_idx_dir = os.path.join(common_patient_data_dir, cancer_name, cancer_stage_rep)
            out_idx_fp = os.path.join(out_idx_dir, common_patient_file_name)
            common_submitter_ids = read_tab_seperated_file_and_get_target_column(1, out_idx_fp)

            for sidx, csid in enumerate(common_submitter_ids):

                for key, iclass in SNP_Ins_Del_classification.items():
                    out_mut_fp = os.path.join(out_idx_dir, str(sidx + 1) + '_' + key + '.tsv')
                    with open(out_mut_fp, "w") as out_mut_file:
                        ltws = []
                        for gidx, gene_name in enumerate(GENOME):
                            if len(SNP_INS_DEL_dict[iclass][cancer_stage_rep][csid][gene_name].keys()) > 0:
                                sorted_dict = SNP_INS_DEL_dict[iclass][cancer_stage_rep][csid][gene_name]
                                sorted_dict = sorted(sorted_dict.items(), key=lambda d: d[0])
                                ltw_t = []
                                if key == "SNP":
                                    for k, v in sorted_dict:
                                        ltw_temp = ",".join([str(v['rel_start']),str(v['v_start']), str(v['class_code']), str(v['class']),v['change'],v['context']])
                                        ltw_t.append(ltw_temp)
                                else:
                                    for k, v in sorted_dict:
                                        ltw_temp = ",".join([str(v['rel_start']),str(v['rel_end']),str(v['v_start']),str(v['v_end']), str(v['class_code']) , str(v['class']), v['change']])
                                        ltw_t.append(ltw_temp)
                                ltw = str(gidx + 1) + "\t" + gene_name +"\t" + ";".join(ltw_t)
                                ltws.append(ltw)
                        out_mut_file.write("\n".join(ltws))
                        print "write %s successful" % out_mut_fp

        t1 = time.time()
        t_used_time = t1 - t0
        tot_time += t_used_time
        remain_time = (tot_time / (cid + 1.0)) * (len(cancer_names) - cid - 1.0)
        print "%d of %d, %.2f, Total Time: %.2f, Time Left: %.2f" % (cid + 1, len(cancer_names), (cid + 1.0) / len(cancer_names), tot_time, remain_time)

def compute_common_mutation_or_methy_variation_samples():
    mutation_stage = "i"
    pvalue_label = "p"
    sig_file_names = ["significant_genes"]#,"significant_genes", "significant_no_genes","significant_methy_genes", "significant_mut_genes"
    gene_classification_dir_name = "gene_classification_mp_0.8_mut_0.1"

    for cancer_name in cancer_names:
        if cancer_name == "COAD":
            mutation_dat_fp = os.path.join(snv_intermidiate_dir, dname, cancer_name, cancer_name + "_" + mutation_stage.replace(" ", "_") + "_mutation_data.dat")
            mutation_matrix = pd.read_csv(mutation_dat_fp, sep='\t', lineterminator='\n', header=0, index_col=0, dtype=np.int32).values

            pvalue_dat_fp = os.path.join(methy_pvalue_dir, dname, cancer_name, cancer_name + "_p" + pvalue_label + "_value.dat")
            pvalue_matrix = pd.read_csv(pvalue_dat_fp, sep='\t', lineterminator='\n', header=0, index_col=0, dtype=np.float64).values

            for sig_file_name in sig_file_names:
                significant_gene_ind_fp = os.path.join(intermediate_file_dir, gene_classification_dir_name , sig_file_name + ".ind")
                sig_gidxs = [int(item) - 1 for item in read_tab_seperated_file_and_get_target_column(1, significant_gene_ind_fp)]
                len_gidx = len(sig_gidxs)

                common_mutation_samples_matrix = np.ones((len_gidx + 1, len_gidx + 1)) * (-1)
                common_methy_variation_samples_matrix = np.ones((len_gidx + 1, len_gidx + 1)) * (-1)
                common_mutation_samples_matrix[0][0] = 0
                common_methy_variation_samples_matrix[0][0] = 0

                for mi, gidx1 in enumerate(sig_gidxs):
                    common_mutation_samples_matrix[mi + 1][0] = mi + 1
                    common_mutation_samples_matrix[0][mi + 1] = mi + 1

                    common_methy_variation_samples_matrix[mi + 1][0] = mi + 1
                    common_methy_variation_samples_matrix[0][mi + 1] = mi + 1

                    for mj, gidx2 in enumerate(sig_gidxs):
                        mut_dat_g1 = np.array(mutation_matrix[gidx1])

                        pval_dat_g1 = np.array(pvalue_matrix[gidx1])
                        if mi == mj:
                            common_mutation_samples_matrix[mi + 1][mj + 1] = (mut_dat_g1 > 0).sum()
                            common_methy_variation_samples_matrix[mi + 1][mj + 1] = (pval_dat_g1 < pvalue_significant_threshold).sum()
                        elif mj < mi:
                            mut_dat_g2 = np.array(mutation_matrix[gidx2])
                            common_mut_arr = np.logical_and(mut_dat_g1 > 0, mut_dat_g2 > 0)
                            common_mutation_sample_count = (common_mut_arr > 0).sum()
                            common_mutation_samples_matrix[mi + 1][mj + 1] = common_mutation_sample_count
                            common_mutation_samples_matrix[mj + 1][mi + 1] = common_mutation_sample_count

                            pval_dat_g2 = np.array(pvalue_matrix[gidx2])
                            common_pval_arr = np.logical_and(pval_dat_g1 < pvalue_significant_threshold, pval_dat_g2 < pvalue_significant_threshold)
                            common_methy_variation_sample_count = (common_pval_arr > 0).sum()
                            common_methy_variation_samples_matrix[mi + 1][mj + 1] = common_methy_variation_sample_count
                            common_methy_variation_samples_matrix[mj + 1][mi + 1] = common_methy_variation_sample_count
                    print "%s %d of %d" % (sig_file_name, mi, len_gidx)
                out_dir = os.path.join(common_sample_cnt_dir, dname, cancer_name, sig_file_name)
                if not os.path.exists(out_dir):
                    os.makedirs(out_dir)
                out_gidxs_fp = os.path.join(out_dir, "sig_gidxs.dat")
                write_tab_seperated_file_for_a_list(out_gidxs_fp, read_tab_seperated_file_and_get_target_column(1, significant_gene_ind_fp), index_included=True)

                out_common_mut_fp = os.path.join(out_dir, "common_mutation_samples_cnt.dat")
                np.savetxt(out_common_mut_fp, common_mutation_samples_matrix, delimiter="\t")
                print "save %s successful!" % out_common_mut_fp

                out_common_pval_fp = os.path.join(out_dir, "common_methy_sig_variation_samples_cnt.dat")
                np.savetxt(out_common_pval_fp, common_methy_variation_samples_matrix, delimiter="\t")
                print "save %s successful!" % out_common_pval_fp

def normal_mean_cpg_methy(target_gene_name, cancer_name):
    gene_infos = {}
    for gidx, gene_name in enumerate(GENOME):
        chr_no, start, end, strand = gene_pos_labels_used[gidx]
        gene_infos[gene_name] = {'chr': chr_no, 'start': start, 'end': end, 'strand': strand}

    cancer_stage_rep = "normal"

    out_idx_dir = os.path.join(common_patient_data_dir, cancer_name, cancer_stage_rep)
    if not os.path.exists(out_idx_dir):
        os.makedirs(out_idx_dir)

    #normal期所有病人文件名列表
    uuid_fp = os.path.join(methy_intermidiate_dir, dname, cancer_name, cancer_name + "_" + cancer_stage_rep + "_uuids.txt")
    uuids = read_tab_seperated_file_and_get_target_column(1, uuid_fp)

    common_filenames = [uuid_to_filename[uuid] for uuid in uuids]

    mean_methy_dict = {}
    for fidx, fname in enumerate(common_filenames):
        print "%s %s %d of %d" % (cancer_name, cancer_stage_rep, fidx + 1, len(common_filenames))
        methy_fp = os.path.join(dna_methy_data_dir, cancer_name, fname)
        with open(methy_fp, "r") as methy_file:
            #因为计算每个基因在promoter(通过distance_to_tss判断)和gene_body(通过cpg位置判断), CpG位点的甲基化水平
            methy_file.readline()
            line = methy_file.readline()
            while line:
                try:
                    line_contents = line.split("\t")
                    gene_symbols = line_contents[5].split(";")
                    beta_val = -1.0 if line_contents[1] == "NA" else float(line_contents[1])
                    gene_types = line_contents[6].split(";")
                    for idx, gene_symbol in enumerate(gene_symbols):
                        if gene_symbol == target_gene_name and gene_types[idx] == "protein_coding" and beta_val > 0.0:
                            cpg_start = int(line_contents[3])
                            g_start = gene_infos[gene_symbol]["start"]
                            g_end = gene_infos[gene_symbol]["end"]
                            g_strand = gene_infos[gene_symbol]["strand"]
                            pttss = cpg_start - g_start if g_strand else g_end - cpg_start
                            #要么是启动子，要么在gene body
                            if (- promoter_length < pttss < 0) or (g_start <= cpg_start <= g_end):
                                mean_methy_dict.setdefault(pttss, [])
                                mean_methy_dict[pttss].append(beta_val)
                                break
                except KeyError, e1:
                    pass
                line = methy_file.readline()
    out_mean_methy_fp = os.path.join(out_idx_dir, 'mean_methy.tsv')
    sorted_mean_methy = sorted(mean_methy_dict.items(), key=lambda d: d[0])
    with open(out_mean_methy_fp,"w") as out_methy_file:
        ltws = []
        for k, v in sorted_mean_methy:
            mean_methy = np.array(v).mean()
            ltws.append("\t".join([str(k), str(mean_methy)]))
        out_methy_file.write("\n".join(ltws))
        print "write %s successful" % out_mean_methy_fp

#计算启动子区的平均甲基化水平,输出成<病人编号,该基因启动子区平均甲基化水平>的文件
def mean_methy_of_promoter(target_gene_name, cancer_name, cancer_stage, xshift = 1):
    gene_infos = {}
    for gidx, gene_name in enumerate(GENOME):
        chr_no, start, end, strand = gene_pos_labels_used[gidx]
        gene_infos[gene_name] = {'chr': chr_no, 'start': start, 'end': end, 'strand': strand}

    cancer_stage_rep = cancer_stage.replace(" ", "_")

    out_idx_dir = os.path.join(common_patient_data_dir, cancer_name, cancer_stage_rep)
    if not os.path.exists(out_idx_dir):
        os.makedirs(out_idx_dir)

    #normal期所有病人文件名列表
    uuid_fp = os.path.join(methy_intermidiate_dir, dname, cancer_name, cancer_name + "_" + cancer_stage_rep + "_uuids.txt")
    uuids = read_tab_seperated_file_and_get_target_column(1, uuid_fp)

    common_filenames = [uuid_to_filename[uuid] for uuid in uuids]

    out_mean_methy_fp = os.path.join(out_idx_dir, 'promoter_mean_methy.tsv')
    with open(out_mean_methy_fp, "w") as out_methy_file:
        for fidx, fname in enumerate(common_filenames):
            promoter_methy = {}
            #fid, 即uuid的编号
            patient_id = fidx + 1
            print "%s %s %d of %d" % (cancer_name, cancer_stage_rep, patient_id, len(common_filenames))
            methy_fp = os.path.join(dna_methy_data_dir, cancer_name, fname)
            with open(methy_fp, "r") as methy_file:
                #因为计算每个基因在promoter(通过distance_to_tss判断)和gene_body(通过cpg位置判断), CpG位点的甲基化水平
                methy_file.readline()
                line = methy_file.readline()
                while line:
                    try:
                        line_contents = line.split("\t")
                        gene_symbols = line_contents[5].split(";")
                        beta_val = -1.0 if line_contents[1] == "NA" else float(line_contents[1])
                        gene_types = line_contents[6].split(";")
                        for idx, gene_symbol in enumerate(gene_symbols):
                            if gene_symbol == target_gene_name and gene_types[idx] == "protein_coding" and beta_val > 0.0:
                                cpg_start = int(line_contents[3])
                                g_start = gene_infos[gene_symbol]["start"]
                                g_end = gene_infos[gene_symbol]["end"]
                                g_strand = gene_infos[gene_symbol]["strand"]
                                pttss = cpg_start - g_start if g_strand else g_end - cpg_start
                                #启动子
                                if - promoter_length < pttss < 0:
                                    promoter_methy[pttss] = beta_val
                                    break
                    except KeyError, e1:
                        pass
                    line = methy_file.readline()
            methy_vals = promoter_methy.values()
            mean_methy_val = float(np.array(methy_vals).mean())
            # out_methy_file.write(str(patient_id) + "\t" + str(mean_methy_val) + "\n")
            ro = random.random() * 0.4 - 0.2 + xshift
            out_methy_file.write(str(ro) + "\t" + str(mean_methy_val) + "\n")
    print "write %s successful" % out_mean_methy_fp
if __name__ == '__main__':
    # extract_submitter_ids_from_methylation_uuids_and_mutation_submitter_ids()
    # obtain_promoter_and_genebody_mutation_status()
    # compute_common_mutation_or_methy_variation_samples()
    # normal_mean_cpg_methy("APC", "COAD")
    mean_methy_of_promoter("APC", "COAD", "normal", xshift= 1)
    mean_methy_of_promoter("APC", "COAD", "i", xshift= 2)