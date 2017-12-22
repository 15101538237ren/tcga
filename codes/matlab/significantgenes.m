function significantgenes()
genes={'BRCA';'COAD';'LIHC';'KIRC';'KIRP';'LUAD';'LUSC';'THCA'};
fpre = '../../data/intermediate_file/';
mp_score_threshold = 0.8;
mutation_rate_threshold = 0.2;
base_path = strcat(fpre, 'gene_classification_mp_',num2str(mp_score_threshold),'_mut_',num2str(mutation_rate_threshold),'/');
%base_path = '../../data/intermediate_file/gene_classification/';
J=[];
for i=1:8
    for j=1:3
        A=load(strcat(base_path, char(genes(i)),'/',char(genes(i)),'_genome_class_',num2str(j),'.dat'));
        J=union(J,A);
    end
end
end
J0=sort(J);
dlmwrite(strcat(base_path,'genes_sig.ind'),J0);
