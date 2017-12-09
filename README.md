## Gene Index file
1.global_files/gene_idx.txt

The meaning of the columns in gene_idx.txt is as following:

(1) : gene_index, start from 1

(2) : HUGO genename, ref: https://www.genenames.org/

## Gene Label file

2.global_files/gene_label.dat

The meaning of the columns in gene_label.dat is as following:

(1) : gene_index, the same as gene_idx.txt

(2) : whether this gene contains CGI before 2000bp of its TSS

(3) : whether this gene is a DNA-binding RNA polymerase II transcription factor, ref: http://www.tfcheckpoint.org/

(4) : the category of this gene, 0: other, 1: onco gene only, 2: tsg only, 3: both oncogene and tsg

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

