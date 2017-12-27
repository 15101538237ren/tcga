function plot_common_sample_count()
clc;
clear all;
cancers={'COAD'};
fpre = '../../data/intermediate_file/';
base_path = strcat(fpre, 'common_sample_cnt/merged_stage/');
gene_idx_filepath = '../../global_files/gene_idx.txt';
[gidxs, gene_names] = textread(gene_idx_filepath,'%d\t%s');
gene_types = {'significant_genes';'significant_no_genes'};
for i=1 : length(cancers)
   cancer_name = char(cancers(i));
   for j=1 : length(gene_types)
       gene_type = char(gene_types(j));
       dir_path = strcat(base_path,cancer_name,'/',gene_type);
       methy_data = load(strcat(dir_path,'/','common_methy_sig_variation_samples_cnt.dat'));
       mut_data = load(strcat(dir_path,'/','common_mutation_samples_cnt.dat'));
       sig_gidxs = load(strcat(dir_path,'/','sig_gidxs.dat'));
       len_of_gene = length(sig_gidxs(:,1))
       HeatMap(methy_data(2:end,2:end));
       HeatMap(mut_data(2:end,2:end));
       disp('hello');
   end
end
end