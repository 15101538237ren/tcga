function [gene_indexs, gene_names] = match_gene_name()
gene_idx_filepath = '../../global_files/gene_idx.txt';
[gidxs, gene_names] = textread(gene_idx_filepath,'%d\t%s');
target_gene_idxs = load('significant_genes.ind');
target_gene_names = gene_names(target_gene_idxs);
fid = fopen('target_gene_index_and_names.txt','w');
for i=1:length(target_gene_idxs)
    fprintf(fid,'%d\t%s\n',target_gene_idxs(i,1),char(target_gene_names{i,:}));
end
fclose(fid);
end