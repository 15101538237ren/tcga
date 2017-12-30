# -*- coding:utf-8 -*-
from base import *
import pandas as pd
import numpy as np
is_merge_stage = True
common_stages = mutation_merged_stage[0 : -1] if is_merge_stage else mutation_stage[0 : -1]
dname = "merged_stage" if is_merge_stage else "stage"
methy_metadata_path = os.path.join(global_files_dir, "methy_24_cancer_meta.json")
methy_metadata_json_obj = json.load(open(methy_metadata_path, 'r'))
common_cases_path = os.path.join(global_files_dir, "common_cases_of_methy_and_mutation.txt")

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
    ret_dict = {}
    for cancer_name in cancer_names:
        ret_dict[cancer_name] = {}
        methy_tot = 0
        mut_tot = 0
        common_tot = 0
        ltw_of_stage_common = []
        for cancer_stage in common_stages:
            ret_dict[cancer_name][cancer_stage] = {}
            cancer_stage_rep = cancer_stage.replace(" ", "_")
            out_idx_dir = os.path.join(common_patient_data_dir, cancer_name, cancer_stage_rep)
            if not os.path.exists(out_idx_dir):
                os.makedirs(out_idx_dir)
            mutation_submitter_ids_fp = os.path.join(snv_intermidiate_dir, dname, cancer_name, cancer_name + "_" + cancer_stage_rep + "_submitter_ids.txt")
            mutation_submitter_ids = read_tab_seperated_file_and_get_target_column(1, mutation_submitter_ids_fp)
            mut_tot += len(mutation_submitter_ids)

            methy_uuids_fp = os.path.join(methy_intermidiate_dir, dname, cancer_name, cancer_name + "_" + cancer_stage_rep + "_uuids.txt")
            methy_uuids = read_tab_seperated_file_and_get_target_column(1, methy_uuids_fp)
            methy_submitter_ids = query_submitter_id_of_a_uuid(methy_uuids)
            methy_file_name_dict = {submitter_id : uuid_to_filename[methy_uuids[sidx]]  for sidx, submitter_id in enumerate(methy_submitter_ids)}

            ret_dict[cancer_name][cancer_stage]['mutation_submitter_ids'] = mutation_submitter_ids
            ret_dict[cancer_name][cancer_stage]['methy_uuids'] = methy_uuids
            ret_dict[cancer_name][cancer_stage]['methy_submitter_ids'] = methy_submitter_ids

            methy_tot += len(methy_uuids)
            if len(methy_uuids) == len(methy_submitter_ids):
                common_submitter_ids = list(set(methy_submitter_ids).intersection(set(mutation_submitter_ids)))
                common_tot += len(common_submitter_ids)
                ret_dict[cancer_name][cancer_stage]['common_submitter_ids'] = common_submitter_ids
                ret_dict[cancer_name][cancer_stage]['common_methy_file_names'] = [methy_file_name_dict[item] for item in common_submitter_ids]

                ltw_of_stage_common.append(str(len(common_submitter_ids)))

                out_idx_fp = os.path.join(out_idx_dir, 'common_patients_idx.txt')
                with open(out_idx_fp, "w") as out_idx_file:
                    tltws = []
                    for common_idx in range(len(common_submitter_ids)):
                        cidx = str(common_idx + 1)
                        c_submitter_id =  common_submitter_ids[common_idx]
                        c_methy_filename = ret_dict[cancer_name][cancer_stage]['common_methy_file_names'][common_idx]
                        tltws.append("\t".join([cidx, c_submitter_id, c_methy_filename]))
                    out_idx_file.write('\n'.join(tltws))
        ltw_arr = [cancer_name, str(methy_tot), str(mut_tot), str(common_tot)]
        ltw_arr.extend(ltw_of_stage_common)
        ltw = "\t".join(ltw_arr)
        ltws.append(ltw)

    with open(common_cases_path,"w") as common_cases_file:
        common_cases_file.write("\n".join(ltws))
        print "write %s successful" % common_cases_path
    return [ret_dict]

def obtain_promoter_and_genebody_methy_status(rtn_list):
    ret_dict = rtn_list[0]

    for cancer_name in cancer_names:
        for cancer_stage in common_stages:
            common_submitter_ids = ret_dict[cancer_name][cancer_stage]['common_submitter_ids']
            common_filenames = ret_dict[cancer_name][cancer_stage]['common_methy_file_names']
            ret_dict[cancer_name][cancer_stage]['methy_status'] ={}
            for gidx, gene_name in enumerate(GENOME):
                chr_no, start, end, strand = gene_pos_labels[gidx]
                ret_dict[cancer_name][cancer_stage]['methy_status'][gene_name] = {'chr': chr_no, 'start':start, 'end':end, 'strand': strand, 'submitter_id_promoter_status': {}, 'submitter_id_genebody_status':{}}

            for sidx, csid in enumerate(common_submitter_ids):
                methy_fp = os.path.join(dna_methy_data_dir, cancer_name,common_filenames[sidx])
                with open(methy_fp, "r") as methy_file:
                    #因为计算每个基因在promoter和gene_body, CpG位点在所有样本的平均甲基化水平
                    methy_file.readline()
                    line = methy_file.readline()
                    while line:
                        line_contents = line.split("\t")
                        try:
                            gene_symbols = line_contents[5].split(";")
                            positions_to_tss = line_contents[8].split(";")
                            beta_val = -1.0 if line_contents[1] == "NA" else float(line_contents[1])
                            gene_types = line_contents[6].split(";")
                            for idx, gene_symbol in enumerate(gene_symbols):
                                # if gene_symbol != "." and gene_types[idx] == "protein_coding" and beta_val > 0.0:
                                pass
                        except KeyError, e1:
                            pass

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
if __name__ == '__main__':
    rtn_list = extract_submitter_ids_from_methylation_uuids_and_mutation_submitter_ids()
    # obtain_promoter_and_genebody_methy_status(rtn_list)
    # compute_common_mutation_or_methy_variation_samples()