function plot_mutation_classification_pipline()
target_gene_idx_fp = './target_gene_list.tsv';
[gene_indexs, gene_names] = textread(target_gene_idx_fp,'%d\t%s');
cancer_name = 'COAD';
for idx=1 : length(gene_indexs)
    gene_name = char(gene_names(idx));
    gene_id = gene_indexs(idx);
    get_gene_info(cancer_name, gene_name, gene_id);
    plot_mutation_classification(cancer_name, gene_name, gene_id);
    plot_mutation_classification_single_normal(cancer_name, gene_name, gene_id);
    plot_mutation_classification_single_i(cancer_name, gene_name, gene_id);
end
end