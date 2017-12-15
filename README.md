1. 熵的数据在/disk/tcga/data/intermediate_file/methy_entropy/merged_stage​
    非法熵用-1表示，癌症阶段列表index在文件夹下的stage_idx.txt

2. log10(p-value)的数据在：/disk/tcga/data/intermediate_file/methy_pvalue/merged_stage/​​​
    非法的log10(p-value)设置为了10（相当于pvalue = 10^10）
    *_pp_value.dat​和_pn_value.dat分别代表p+​value和p-value的数据
    *_p_score.dat​和_n_score.dat分别代表m+score和m-score的数据，该文件最后一列为score值。

3. 基因索引在:global_files/gene_idx.txt
    第1列: 从1开始的基因索引编号
    第2列: HUGO 基因名, 参考链接: https://www.genenames.org/

4. 基因标注文件:global_files/gene_label.dat
    第1列: 基因索引,同gene_idx.txt的索引
    第2列: 当前基因是否在上游2000bp-TSS范围内包含CpG岛, 1包含,0不包含
    第3列: 当前基因是否是转录因子, 参考链接:http://www.tfcheckpoint.org/
    第4列: 基因的抑癌,原癌类别(ref TSG2.0, OnCoGene): 1:仅为原癌基因, 2.仅为抑癌基因 3.既为原癌又为抑癌基因 0:其他基因
    第5列: 基因的抑癌,原癌类别(ref Vogelstein): 1:仅为原癌基因, 2.仅为抑癌基因 3.既为原癌又为抑癌基因 0:其他基因

## Means and Stds of methylation data

3.data/intermediate_file/methy_mean_std/merged_stage/Cancer/

3.1 *_stages.txt

	list the stage index and its names

3.2 *_mean.dat
	
	row idx: gene idx
	
	column idx: stage idx

	data: mean of the beta values for the all cases in a specific stage

3.3 *_std.dat
	
	row idx: gene idx
	
	column idx: stage idx

	data: the standard deviation of the beta values for the all cases in a specific stage

3.4 codes/matlab/plot_methy_mean_and_std.m
	The matlab code to plot the mean and std of methylation data

3.5 figures/methy_mean_and_std/Cancer.eps
	The output figures of the plot_methy_mean_and_std.m program

RNA Data

4. data/intermediate_file/rna_intermidiate/merged_stage/

4.1 ensembl_id_idx.txt
	
	1st column: gene index
	
	2nd column: the corresponding ensembl ids

4.2 Cancer/Cancer_stage_case_ids.txt
	
	the case ids of each stage of a particular cancer

4.3 Cancer/Cancer_stage_fpkm.dat
	
	row idx: gene idx
	
	column idx: case idx

	data: the fpkm data (not comparable between cases) extracted from the original *.FPKM.txt, which is downloaded directly from TCGA

4.3 Cancer/Cancer_stage_tpm.dat
	
	row idx: gene idx
	
	column idx: case idx

	data: the tpm (relative value) data which is more comparable than fpkm data between cases, the total value of each column is 10^6.

