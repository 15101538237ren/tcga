#coding: utf-8

import sys
sys.path.append('../matlab')
import get_gene_data
reload(get_gene_data)
from sklearn.cluster import Birch    
from sklearn.cluster import KMeans

import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

#注意参数stage要将"i_th"和"normal_th"当成list或tuple传入
def get_cpi_sample_mat(dir_path, cancer_name, stage_name, gene_id):
	i_th_sample_mat = get_gene_data.get_all_sample_methy_mat(dir_path, cancer_name, stage_name[0], gene_id)
	normal_sample_mat = get_gene_data.get_all_sample_methy_mat(dir_path, cancer_name, stage_name[1], gene_id)
	cpi_sample_mat = np.hstack((i_th_sample_mat[:, 1:], normal_sample_mat[:, 1:]))  #按行合并，扩展列数
	cpi_sample_mat = cpi_sample_mat.transpose()

	with open('test.csv', 'w') as f:
		for i in range(cpi_sample_mat.shape[0]):
			f.write(str(cpi_sample_mat[i, 0]))
			for j in range(1, cpi_sample_mat.shape[1]):
				f.write(', ' + str(cpi_sample_mat[i, j]))
			f.write('\n')

	print cpi_sample_mat.shape
	# Kmeans聚类  
	clf = KMeans(n_clusters=2)
	s = clf.fit(cpi_sample_mat) #加载数据集合
	#print clf.cluster_centers_
	centroids = clf.labels_
	#print centroids,type(centroids) #显示中心点
	print clf.inertia_  #显示聚类效果

	pca = PCA(n_components=2)
	res = pca.fit(cpi_sample_mat)
	#print res
	print pca.explained_variance_ratio_
	new_cpi_sample_mat =  pca.transform(cpi_sample_mat)
	print new_cpi_sample_mat[:, 0]
	#print 'shape :', new_cpi_sample_mat.shape
	with open('test1.csv', 'w') as f:
		for i in range(new_cpi_sample_mat.shape[0]):
			f.write(str(new_cpi_sample_mat[i, 0]))
			for j in range(1, new_cpi_sample_mat.shape[1]):
				f.write(', ' + str(new_cpi_sample_mat[i, j]))
			f.write('\n')
	fig = plt.figure()
	ax = fig.add_subplot(111)
	scatter1 = ax.scatter(new_cpi_sample_mat[:42, 0], new_cpi_sample_mat[:42, 1], label='Line 1', s=7, c='red')
	scatter2 = ax.scatter(new_cpi_sample_mat[43:, 0], new_cpi_sample_mat[43:, 1], label='Line 2', s=7, c='blue')
	ax.set_xlabel('PCA1')
	ax.set_ylabel('PCA2')
	ax.set_title('COAD APC cluster PCA visualization')
	ax.legend([scatter1, scatter2], ['I-th sample', 'Normal sample'])
	plt.show()
