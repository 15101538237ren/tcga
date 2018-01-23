global mp_score_threshold mutation_rate_threshold
mp_score_threshold = 0.8;
mutation_rate_threshold = 0.1;
genes={'BRCA';'COAD';'LIHC';'KIRC';'KIRP';'LUAD';'LUSC';'THCA'};
fpre = '../../data/intermediate_file/';
base_path = strcat(fpre, 'gene_classification_mp_',num2str(mp_score_threshold),'_mut_',num2str(mutation_rate_threshold),'/');
if ~exist(base_path)
        mkdir(base_path);
end
J=[];
Methy_gidxs = [];
Mut_gidxs = [];
for i=1:8
    for j=1:3
        A=load(strcat(base_path, char(genes(i)),'/',char(genes(i)),'_genome_class_',num2str(j),'.dat'));
        J=union(J,A);
        if (j == 1)
            Methy_gidxs = union(Methy_gidxs,A);
        elseif (j == 2)
            Mut_gidxs= union(Mut_gidxs,A);
        elseif (j==3)
            Methy_gidxs = union(Methy_gidxs,A);
            Mut_gidxs= union(Mut_gidxs,A);
        end
    end
end
J0=sort(J);
dlmwrite(strcat(base_path,'genes_sig.ind'),J0);
match_gene_name(J0,strcat(base_path,'significant_genes_all.ind'));
match_gene_name(Methy_gidxs,strcat(base_path,'significant_methy_genes.ind'));
match_gene_name(Mut_gidxs,strcat(base_path,'significant_mut_genes.ind'));
L=load('../../global_files/gene_label.dat');
Onco=1;
Tsg=2;
Both=3;
Onco_Vogel_indexs = find(L(:,5)==Onco);
TSG_Vogel_indexs = find(L(:,5)==Tsg);
match_gene_name(Onco_Vogel_indexs,strcat(base_path,'Onco_TSG_Vogel.ind'));
match_gene_name(TSG_Vogel_indexs,strcat(base_path,'TSG_Vogel.ind'));