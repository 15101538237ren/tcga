# -*- coding:utf-8 -*-
import requests
import pandas as pd
import numpy as np
from base import *
is_merge_stage = True
common_stages = mutation_merged_stage[0 : -1] if is_merge_stage else mutation_stage[0 : -1]
dname = "merged_stage" if is_merge_stage else "stage"
methy_metadata_path = os.path.join(global_files_dir, "methy_24_cancer_meta.json")
methy_metadata_json_obj = json.load(open(methy_metadata_path, 'r'))

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
    query_size = 50

    for cancer_name in cancer_names:
        for cancer_stage in common_stages:
            cancer_stage_rep = cancer_stage.replace(" ", "_")
            mutation_submitter_ids_fp = os.path.join(snv_intermidiate_dir, dname, cancer_name, cancer_name + "_" + cancer_stage_rep + "_submitter_ids.txt")
            mutation_submitter_ids = read_tab_seperated_file_and_get_target_column(1, mutation_submitter_ids_fp)

            methy_uuids_fp = os.path.join(methy_intermidiate_dir, dname, cancer_name, cancer_name + "_" + cancer_stage_rep + "_uuids.txt")
            methy_uuids = read_tab_seperated_file_and_get_target_column(1, methy_uuids_fp)
            methy_submitter_ids = query_submitter_id_of_a_uuid(methy_uuids)
            if len(methy_uuids) == len(methy_submitter_ids):
                common_submitter_ids = list(set(methy_submitter_ids).intersection(set(mutation_submitter_ids)))
                print "%s %s common sid %d" % (cancer_name, cancer_stage, len(common_submitter_ids))
if __name__ == '__main__':
    extract_submitter_ids_from_methylation_uuids_and_mutation_submitter_ids()