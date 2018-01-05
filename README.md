global_files/

----------------gene_idx.txt
				基因索引 数据

				第1列: 从1开始的基因索引编号
				第2列: HUGO 基因名, 参考链接: https://www.genenames.org/


----------------gene_label.dat
				基因标注 数据

				第1列: 基因索引,同gene_idx.txt的索引
			    第2列: 当前基因是否在上游2000bp-TSS范围内包含CpG岛, 1包含,0不包含
			    第3列: 当前基因是否是转录因子, 参考链接:http://www.tfcheckpoint.org/
			    第4列: 基因的抑癌,原癌类别(ref TSG2.0, OnCoGene): 1:仅为原癌基因, 2.仅为抑癌基因 3.既为原癌又为抑癌基因 0:其他基因
			    第5列: 基因的抑癌,原癌类别(ref Vogelstein): 1:仅为原癌基因, 2.仅为抑癌基因 3.既为原癌又为抑癌基因 0:其他基因
			    第6列: MSK-IMPACT 341个基因的标注, 1表示该基因在这341个基因中, 0表示不在
			    第7列: MSK-IMPACT 410个基因的标注, 1表示该基因在这410个基因中, 0表示不在
			    第8列: 染色体编号, 取值范围1-24。其中23表示X染色体, 24表示Y染色体。非法值为-1
			    第9列: 基因起始位点.(但是不能直接当成TSS。因为如果是反链，需要用第10列的值做TSS), 非法值为-1
			    第10列: 基因结束位点, 非法值为-1
			    第11列: 正反链, 正链 1, 反链0, 非法值为-1

----------------methylation_sample_count.txt
				甲基化样本数 数据

				第1列: 癌症名称
				第2列: normal 样本数
				第3列: i期 样本数
				第4列: ii期 样本数
				第5列: iii期 样本数
				第6列: iv期 样本数
				第7列: x期 样本数
				第8列: normal至x期总样本数

----------------mutation_sample_count.txt
				突变样本数 数据

				第1列: 癌症名称
				第2列: normal 样本数
				第3列: i期 样本数
				第4列: ii期 样本数
				第5列: iii期 样本数
				第6列: iv期 样本数
				第7列: x期 样本数
				第8列: normal至x期总样本数


----------------common_cases_of_methy_and_mutation.txt
				甲基化和突变拥有的共同突变样本数 数据
				
				第1列: 癌症名称
				第2列: 甲基化总样本数
				第3列: 突变总样本数
				第4列: 甲基化-突变共同总样本数
				第5列: i期的甲基化-突变共同样本数
				第6列: ii期的甲基化-突变共同样本数
				第7列: iii期的甲基化-突变共同样本数
				第8列: iv期的甲基化-突变共同样本数


data/intermediate_file/

----------------common_patients_data
				甲基化、突变共同样本的突变和甲基化 数据

				----common_patients_idxs.txt
				共同病人索引(index) 数据

				第1列: patient index (Pidx, 病人编号)
				第2列: TCGA submitter id
				第3列: 对应甲基化数据的文件名

				----common_patients_methy_idxs.txt
				甲基化和突变共同病人的methylation case index (甲基化编号) 数据, 用于从原有甲基化数据中直接提取共同病人的子数据

				第1列: patient index (Pidx, 病人编号)
				第2列: 对应甲基化数据的case编号，即: methy_intermidiate/*/*/*_uuids.txt 中的编号

				----common_patients_mutation_idxs.txt
				甲基化和突变共同病人的mutation case index (突变编号) 数据, 用于从原有突变数据中直接提取共同病人的子数据

				第1列: patient index (Pidx, 病人编号)
				第2列: 对应甲基化数据的case编号，即: methy_intermidiate/*/*/*_uuids.txt 中的编号

				----j_SNP.tsv
				第j个甲基化、突变的共同病人(j属于Pidxs)的(启动子区和gene body区)的点突变数据 (\t分隔符)
				
				第1列: gene index
				第2列: gene name
				第3-n列: 若干点突变分别的详细信息, 每个点突变用 "\t" 分开, 详细信息用 "," 分开。
						点突变详细信息的数据顺序:
							1. 突变位点相对于TSS的位置p: p<0 即在启动子区; p > 0, 在 gene body 区。
							2. 突变在基因组上的位置
							3. 


1. DNA甲基化熵
	位置: /disk/tcga/data/intermediate_file/methy_entropy/merged_stage​
    
    非法熵用-1表示，癌症阶段列表index在文件夹下的stage_idx.txt

2. 甲基化i期相对normal的log10(p-value)
	位置: /disk/tcga/data/intermediate_file/methy_pvalue/merged_stage/​​​
    
    非法的log10(p-value)设置为了10（相当于pvalue = 10^10）
    
    *_pp_value.dat​和_pn_value.dat分别代表p+​value和p-value的数据
    *_p_score.dat​和_n_score.dat分别代表m+score和m-score的数据，该文件最后一列为score值。


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

