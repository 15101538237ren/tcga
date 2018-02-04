#coding: utf-8

import sys
sys.path.append('../matlab')
import get_gene_data
reload(get_gene_data)
from sklearn.cluster import Birch    
from sklearn.cluster import KMeans

import numpy as np

#注意参数stage要将"i_th"和"normal_th"当成list或tuple传入
def get_cpi_sample_mat(dir_path, cancer_name, stage_name, gene_id):
	i_th_sample_mat = get_gene_data.get_all_sample_methy_mat(dir_path, cancer_name, stage_name['i'], gene_id)
	normal_sample_mat = get_gene_data.get_all_sample_methy_mat(dir_path, cancer_name, stage_name['normal'], gene_id)
	cpi_sample_mat = np.hstack((i_th_sample_mat, normal_sample_mat))  #按行合并，扩展列数
	# Kmeans聚类  
	clf = KMeans(n_clusters=2)
	s = clf.fit(cpi_sample_mat) #加载数据集合
	numSamples=len(cpi_sample_mat)
	centroids = clf.labels_
	print centroids,type(centroids) #显示中心点
	print clf.inertia_  #显示聚类效果 
