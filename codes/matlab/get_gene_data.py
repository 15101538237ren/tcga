#coding: utf-8

import os
import pandas as pd
import numpy as np

def get_sample_num(dir_path, cancer_name, stage_name, w_file=True):
	cancer_dir = os.path.join(dir_path, cancer_name)
	#print(cancer_dir)
	stage_dir = os.path.join(cancer_dir, stage_name)
	#print(stage_dir)
	sample_id_file_name = 'common_patients_idxs.txt'
	sample_id_file_path = os.path.join(stage_dir, sample_id_file_name)
	sample_id_file_path = sample_id_file_path.replace('\\','/')
	#print(sample_id_file_path)
	file = open(sample_id_file_path, 'r')
	content = file.readlines()
	file.close()
	if(w_file):
		sample_num_file_name = 'out_sample_num.txt'
		file = open(sample_num_file_name, 'w')
		file.write(str(len(content))+'\n')
		file.close()
		return
	return len(content)

#获取一种gene的一个样本的甲基化水平数据，并写入文件，方便Matlab直接load
def get_sample_methy(dir_path, cancer_name, stage_name, sample_id, gene_id, w_file=True):
	sample_id = int(sample_id)
	gene_id = int(gene_id)
	cancer_dir = os.path.join(dir_path, cancer_name)
	stage_dir = os.path.join(cancer_dir, stage_name)
	methy_file_name = str(sample_id) + '_methy.tsv'
	methy_file_path = os.path.join(stage_dir, methy_file_name)
	out_file_name = 'out_methy.txt'
	out_list = []

	methy_file_path = methy_file_path.replace('\\','/')
	file = open(methy_file_path, 'r')
	content = file.readlines()
	for line in content:
		line_list = [item.strip() for item in line.split('\t')]
		line_gene_id = int(line_list[0])
		if(len(line_list) > 1):
			gene_info = line_list[1]
			if(line_gene_id == gene_id):
				methy_list = [item.strip().split(',') for item in gene_info.split(';')]
				for item in methy_list:
					out_list.append([item[0], item[2]])
				break
	file.close()

	if(w_file):
		file = open(out_file_name, 'w')
		for item in out_list:
			line_str = item[0] + '\t' + item[1] + '\n'
			file.write(line_str)
		file.close()
		return
	return out_list


#获取一种gene的全部样本的平均甲基化水平
def get_all_sample_ave_methy(dir_path, cancer_name, stage_name, gene_id, w_file=True):
	gene_id = int(gene_id)
	cancer_dir = os.path.join(dir_path, cancer_name)
	stage_dir = os.path.join(cancer_dir, stage_name)
	sample_num = get_sample_num(dir_path, cancer_name, stage_name, w_file=False)

	out_file_name = 'out_ave_methy.txt'
	new_methy_dict = {}
	new_num_dict = {}
	for sample_id in range(1, sample_num + 1):
		#print sample_id
		ret_list = get_sample_methy(dir_path, cancer_name, stage_name, sample_id, gene_id, w_file=False)
		#print  ret_list
		if len(ret_list) > 0 and len(ret_list[0]) > 0:
			for item in ret_list:
				pos = int(item[0])
				beta_value = float(item[1])
				if(pos not in new_methy_dict):
					new_methy_dict[pos] = beta_value
					new_num_dict[pos] = 1
				else:
					new_methy_dict[pos] += beta_value
					new_num_dict[pos] += 1

	for k in new_methy_dict.keys():
		new_methy_dict[k] = new_methy_dict[k] / new_num_dict[k]
	#print new_methy_dict
	#print new_num_dict
	new_list = sorted(new_methy_dict.items())

	if(w_file):
		file = open(out_file_name, 'w')
		for item in new_list:
			line_str = str(item[0]) + '\t' + str(item[1]) + '\n'
			file.write(line_str)
		file.close()
		return
	return new_list

def get_all_sample_methy_mat(dir_path, cancer_name, stage_name, gene_id):
	gene_id = int(gene_id)
	cancer_dir = os.path.join(dir_path, cancer_name)
	stage_dir = os.path.join(cancer_dir, stage_name)
	sample_num = get_sample_num(dir_path, cancer_name, stage_name, w_file=False)

	new_methy_dict = {}
	for sample_id in range(1, sample_num + 1):
		#print sample_id
		ret_list = get_sample_methy(dir_path, cancer_name, stage_name, sample_id, gene_id, w_file=False)
		#print  ret_list
		if len(ret_list) > 0 and len(ret_list[0]) > 0:
			for item in ret_list:
				pos = int(item[0])
				beta_value = float(item[1])
				new_methy_dict[pos] = new_methy_dict.get(pos, [])
				new_methy_dict[pos].append(beta_value)

	#print new_methy_dict
	new_list = sorted(new_methy_dict.items())
	out_list = []
	for item_tuple in new_list:
		item_list = [item_tuple[0]]
		item_list.extend(item_tuple[1])
		out_list.append(item_list)
	return np.mat(out_list)


#将每一种类型的突变写入文件, 方便Matlab直接load
#Return的list每个start,end都保证为int(type_num为str)，在写文件时再转换成str
def get_each_mutation(file_path, mutation_type, gene_id, w_file=True, simple_m_type=True):
	gene_id_idx = 0
	gene_info_idx = 2
	mutation_start_idx = 0
	mutation_end_idx = 1
	mutation_snp_type_num_idx = 2
	mutation_ins_del_type_num_idx = 4
	out_file_name = 'out_' + mutation_type + '.txt'
	out_list = []
	file_path = file_path.replace('\\','/')
	file = open(file_path, 'r')
	content = file.readlines()
	for line in content:
		line_list = [item.strip() for item in line.split('\t')]
		line_gene_id = int(line_list[gene_id_idx])
		gene_info = line_list[gene_info_idx]
		if(line_gene_id == gene_id):
			mutation_list = [item for item in gene_info.split(',')]
			if(mutation_type == 'SNP'):
				x = int(mutation_list[mutation_start_idx])
				type_num = mutation_list[mutation_snp_type_num_idx]
				pos_list = [x]
				if(simple_m_type == False):
					pos_list = [x, x, type_num]
			else:
				x = min(int(mutation_list[mutation_start_idx]), int(mutation_list[mutation_end_idx]))
				y = max(int(mutation_list[mutation_start_idx]), int(mutation_list[mutation_end_idx]))
				pos_list = [num for num in range(x, y + 1)]
				if(simple_m_type == False):
					type_num = mutation_list[mutation_ins_del_type_num_idx]
					#print type_num, type(type_num)
					pos_list= [x, y, type_num]
			out_list.append(pos_list)
			break
	file.close()

	if(w_file):
		file = open(out_file_name, 'w')
		for item in out_list:
			file.write(str(item[0]))
			for i in range(1, len(item)):
				line_str = '\t' + str(item[i])
				file.write(line_str)
			file.write('\n')
		file.close()
		return
	return out_list


#获取一个gene的一个样本中的突变情况，并写入文件, 方便Matlab直接load
def get_sample_mutation(dir_path, cancer_name, stage_name, sample_id, gene_id):
	sample_id = int(sample_id)
	gene_id = int(gene_id)
	cancer_dir = os.path.join(dir_path, cancer_name)
	stage_dir = os.path.join(cancer_dir, stage_name)
	INS_file_name = str(sample_id) + '_INS.tsv'
	DEL_file_name = str(sample_id) + '_DEL.tsv'
	SNP_file_name = str(sample_id) + '_SNP.tsv'

	INS_file_path = os.path.join(stage_dir, INS_file_name)
	DEL_file_path = os.path.join(stage_dir, DEL_file_name)
	SNP_file_path = os.path.join(stage_dir, SNP_file_name)

	get_each_mutation(INS_file_path, 'INS', gene_id)
	get_each_mutation(DEL_file_path, 'DEL', gene_id)
	get_each_mutation(SNP_file_path, 'SNP', gene_id)


#将原始list转换成(pos, mutation_num)的形式
def reorder_list(origin_list):
	new_dict = {}
	if(len(origin_list)==0):
		return origin_list
	for item_list in origin_list:
		for item in item_list:
			if(item not in new_dict):
				new_dict[item] = 1
			else:
				new_dict[item] += 1
	new_list = sorted(new_dict.items())
	return new_list

def write_file_for_reordered_list(file_path, input_list):
	file_path = file_path.replace('\\','/')
	file = open(file_path, 'w')
	for item in input_list:
		mutation_pos = str(item[0])
		file.write(mutation_pos)
		for num in item[1:]:
			#print num, type(num)
			line_str = '\t' + str(num)
			file.write(line_str)
		file.write('\n')
	file.close()

def get_mutation_list(stage_dir, gene_id, sample_num, simple_m_type=True):
	out_ins_list = []
	out_del_list = []
	out_snp_list = []

	for sample_id in range(1, sample_num + 1):
		INS_file_name = str(sample_id) + '_INS.tsv'
		DEL_file_name = str(sample_id) + '_DEL.tsv'
		SNP_file_name = str(sample_id) + '_SNP.tsv'

		INS_file_path = os.path.join(stage_dir, INS_file_name)
		DEL_file_path = os.path.join(stage_dir, DEL_file_name)
		SNP_file_path = os.path.join(stage_dir, SNP_file_name)

		out_ins = get_each_mutation(INS_file_path, 'INS', gene_id, w_file=False, simple_m_type=simple_m_type)
		out_del = get_each_mutation(DEL_file_path, 'DEL', gene_id, w_file=False, simple_m_type=simple_m_type)
		out_snp = get_each_mutation(SNP_file_path, 'SNP', gene_id, w_file=False, simple_m_type=simple_m_type)

		if(len(out_ins) > 0 and len(out_ins[0]) > 0):
			out_ins_list.append(out_ins[0])
		if(len(out_del) > 0 and len(out_del[0]) > 0):
			out_del_list.append(out_del[0])
		if(len(out_snp) > 0 and len(out_snp[0]) > 0):
			out_snp_list.append(out_snp[0])

	return [out_ins_list, out_del_list, out_snp_list]

def get_mutation_list_by_sample(stage_dir, gene_id, sample_num, simple_m_type=True):
	out_ins_list = []
	out_del_list = []
	out_snp_list = []

	for sample_id in range(1, sample_num + 1):
		INS_file_name = str(sample_id) + '_INS.tsv'
		DEL_file_name = str(sample_id) + '_DEL.tsv'
		SNP_file_name = str(sample_id) + '_SNP.tsv'

		INS_file_path = os.path.join(stage_dir, INS_file_name)
		DEL_file_path = os.path.join(stage_dir, DEL_file_name)
		SNP_file_path = os.path.join(stage_dir, SNP_file_name)

		out_ins = get_each_mutation(INS_file_path, 'INS', gene_id, w_file=False, simple_m_type=simple_m_type)
		out_del = get_each_mutation(DEL_file_path, 'DEL', gene_id, w_file=False, simple_m_type=simple_m_type)
		out_snp = get_each_mutation(SNP_file_path, 'SNP', gene_id, w_file=False, simple_m_type=simple_m_type)

		if(len(out_ins) > 0 and len(out_ins[0]) > 0):
			out_ins_list.append(out_ins[0])
		else:
			out_ins_list.append([])
		if(len(out_del) > 0 and len(out_del[0]) > 0):
			out_del_list.append(out_del[0])
		else:
			out_del_list.append([])
		if(len(out_snp) > 0 and len(out_snp[0]) > 0):
			out_snp_list.append(out_snp[0])
		else:
			out_snp_list.append([])
	return [out_ins_list, out_del_list, out_snp_list]


#获取一个gene的全部样本对应的mutation数量
def get_all_sample_mutation(dir_path, cancer_name, stage_name, gene_id, w_file=True, simple_m_type=True, loop_method=get_mutation_list):
	gene_id = int(gene_id)
	cancer_dir = os.path.join(dir_path, cancer_name)
	stage_dir = os.path.join(cancer_dir, stage_name)
	sample_num = get_sample_num(dir_path, cancer_name, stage_name, w_file=False)
	
	[out_ins_list, out_del_list, out_snp_list] = loop_method(stage_dir, gene_id, sample_num, simple_m_type=simple_m_type)
	#print 'out_ins_list:', out_ins_list
	#print 'out_del_list:', out_del_list
	#print 'out_snp_list:', out_snp_list

	if(simple_m_type):
		out_ins_list = reorder_list(out_ins_list)
		out_del_list = reorder_list(out_del_list)
		out_snp_list = reorder_list(out_snp_list)
	if(w_file):
		out_ins_file_name = 'out_INS_num.txt'
		out_del_file_name = 'out_DEL_num.txt'
		out_snp_file_name = 'out_SNP_num.txt'
		write_file_for_reordered_list(out_ins_file_name, out_ins_list)
		write_file_for_reordered_list(out_del_file_name, out_del_list)
		write_file_for_reordered_list(out_snp_file_name, out_snp_list)
		return
	return {'INS': out_ins_list, 'DEL':out_del_list, 'SNP':out_snp_list}

def get_methy_idx_mapping_list(mapping_file_path):
	mapping_list = []
	mapping_file_path = mapping_file_path.replace('\\','/')
	file = open(mapping_file_path, 'r')
	content = file.readlines()
	for line in content:
		line_list = [item.strip() for item in line.split('\t')]
		common_id = line_list[0]
		origin_methy_sample_id = line_list[1]
		mapping_list.append((common_id, origin_methy_sample_id))
	file.close()
	return mapping_list


def get_sample_methy_p_value(common_dir_path, p_value_dir_path, cancer_name, stage_name, gene_id):
	gene_id = int(gene_id)
	cancer_dir = os.path.join(common_dir_path, cancer_name)
	stage_dir = os.path.join(cancer_dir, stage_name)
	methy_mapping_file_name = 'common_patients_methy_idxs.txt'
	methy_mapping_file_path = os.path.join(stage_dir, methy_mapping_file_name)
	methy_idx_mapping_list = get_methy_idx_mapping_list(methy_mapping_file_path)

	p_value_file_name = cancer_name + '_pp_value.dat'
	p_value_file_path = os.path.join(p_value_dir_path, p_value_file_name)
	p_value_list = []

	#print 'methy_idx:', methy_idx_mapping_list

	df = pd.read_csv(p_value_file_path, sep='\t', header=0)

	#print df.head()

	for item in methy_idx_mapping_list:
		common_id = item[0]
		origin_sample_id = item[1]
		p_value = df.loc[gene_id - 1, origin_sample_id]
		p_value_list.append(p_value)
	#print df.head()

	return p_value_list


def get_mutation_classification(dir_path):
	file_name = 'mutation_classification.txt'
	file_path = os.path.join(dir_path, file_name)
	file_path = file_path.replace('\\','/')
	mutation_class_id_idx = 0
	mutation_class_name_idx = 1
	mutation_class = {}

	file_path = file_path.replace('\\','/')
	file = open(file_path, 'r')
	content = file.readlines()
	for line in content:
		line_list = [item.strip() for item in line.split('\t')]
		type_id = line_list[mutation_class_id_idx]
		type_name = line_list[mutation_class_name_idx]
		mutation_class[type_id] = type_name
	file.close()
	return mutation_class

def get_gene_name(dir_path, gene_id):
	file_name = 'gene_idx.txt'
	file_path = os.path.join(dir_path, file_name)
	file_path = file_path.replace('\\','/')
	gene_id_idx = 0
	gene_name_idx = 1
	gene_dict = {}

	file = open(file_path, 'r')
	content = file.readlines()
	for line in content:
		line_list = [item.strip() for item in line.split('\t')]
		gene_idx = line_list[gene_id_idx]
		gene_name = line_list[gene_name_idx]
		gene_dict[gene_idx] = gene_name
	file.close()
	return gene_dict[str(gene_id)]

def get_mutation_info_output(mutation_list, mutation_class_dict):
	if(len(mutation_list) == 0):
		return ''
	start_pos = str(mutation_list[0])
	end_pos = str(mutation_list[1])
	mutation_name = mutation_class_dict[mutation_list[2]]
	out_str = start_pos + ',' + end_pos + ',' + mutation_name + ';'
	return out_str

def write_methy_mutation_file(cancer_name, gene_name, sample_num, p_value_list, mutation_info_dict, mutation_class_dict):

	file_path = cancer_name + '_' + gene_name + '_Sample_Methy_Mutation_Info.tsv'
	file = open(file_path, 'w')

	for i in range(1, sample_num + 1):
		sample_id = str(i)
		pp_value = str(p_value_list[i - 1])
		sample_ins_list = mutation_info_dict['INS'][i - 1]
		sample_del_list = mutation_info_dict['DEL'][i - 1]
		sample_snp_list = mutation_info_dict['SNP'][i - 1]
		line_str = sample_id + '\t' + pp_value
		mutation_str = get_mutation_info_output(sample_ins_list, mutation_class_dict) + \
			get_mutation_info_output(sample_del_list, mutation_class_dict) + \
			get_mutation_info_output(sample_snp_list, mutation_class_dict)
		line_str = line_str + '\t' + mutation_str + '\n'
		file.write(line_str)

	file.close()

def output_methy_mutation_file(global_file_dir, common_dir_path, p_value_dir_path, cancer_name, stage_name, gene_id):
	sample_num = get_sample_num(common_dir_path, cancer_name, stage_name, w_file=False)
	p_value_list = get_sample_methy_p_value(common_dir_path, p_value_dir_path, cancer_name, stage_name, gene_id)
	mutation_info_dict = get_all_sample_mutation(common_dir_path, cancer_name, stage_name, gene_id, w_file=False, \
		simple_m_type=False, loop_method=get_mutation_list_by_sample)
	mutation_class_dict = get_mutation_classification(global_file_dir)
	gene_name = get_gene_name(global_file_dir, gene_id)
	#print mutation_info_dict
	write_methy_mutation_file(cancer_name, gene_name, sample_num, p_value_list, mutation_info_dict, mutation_class_dict)
	
