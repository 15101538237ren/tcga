function [gene_indexs, gene_names] = match_gene_name(target_gene_idxs,fname)
global base_path
gene_idx_filepath = '../../global_files/gene_idx.txt';
[gidxs, gene_names] = textread(gene_idx_filepath,'%d\t%s');
%target_gene_idxs = load('significant_genes.ind');
target_gene_names = gene_names(target_gene_idxs);
fid = fopen(strcat(base_path,fname),'w');
for i=1:length(target_gene_idxs)
    fprintf(fid,'%d\t%d\t%s\n',i,target_gene_idxs(i,1),char(target_gene_names{i,:}));
end
fclose(fid);
end