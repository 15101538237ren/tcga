# -*- coding:utf-8 -*-
from base import *
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
    ltws = ["\t".join(['cancer', '#methy cases', '#mutation cases', '#common cases'])]
    ret_dict = {}
    for cancer_name in cancer_names:
        ret_dict[cancer_name] = {}
        methy_tot = 0
        mut_tot = 0
        common_tot = 0
        for cancer_stage in common_stages:
            ret_dict[cancer_name][cancer_stage] = {}
            cancer_stage_rep = cancer_stage.replace(" ", "_")
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
            ret_dict[cancer_name][cancer_stage]['methy_file_name_dict'] = methy_file_name_dict

            methy_tot += len(methy_uuids)
            if len(methy_uuids) == len(methy_submitter_ids):
                common_submitter_ids = list(set(methy_submitter_ids).intersection(set(mutation_submitter_ids)))
                common_tot += len(common_submitter_ids)
                ret_dict[cancer_name][cancer_stage]['common_submitter_ids'] = common_submitter_ids
                ret_dict[cancer_name][cancer_stage]['common_methy_file_names'] = [methy_file_name_dict[item] for item in common_submitter_ids]

        ltw = "\t".join([cancer_name, str(methy_tot), str(mut_tot), str(common_tot)])
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
                                if gene_symbol != "." and gene_types[idx] == "protein_coding" and beta_val > 0.0:
                                    
if __name__ == '__main__':
    rtn_list = extract_submitter_ids_from_methylation_uuids_and_mutation_submitter_ids()
    obtain_promoter_and_genebody_methy_status(rtn_list)