function gene_classification()
clc;
clear all;
close all;
global stage mp_score_threshold mutation_rate_threshold
mp_score_threshold = 0.8;
mutation_rate_threshold = 0.1;

stage = 'i';
cancer_names = {'BRCA';'COAD'; 'LIHC'; 'KIRC'; 'KIRP'; 'LUAD'; 'LUSC'; 'THCA'};
fpre = '../../data/intermediate_file/';
gene_class_path = strcat(fpre, 'gene_classification_mp_',num2str(mp_score_threshold),'_mut_',num2str(mutation_rate_threshold),'/');
class_by_gene_category_path = strcat(fpre, 'class_by_gene_category/');
class_by_cancer_path = strcat(fpre, 'class_by_cancer/');

figure_path = '../../figures/gene_classification/';
paths = {gene_class_path; figure_path;class_by_gene_category_path;class_by_cancer_path};
for i = 1: length(paths)
    if ~exist(char(paths(i)))
        mkdir(char(paths(i)));
    end
end

gene_label = load(strcat('../../global_files/gene_label.dat'));
gene_category_names = {'oncogene_698';'tsg_1018';'oncogene_vogelstein';'tsg_vogelstein';'genome'};

OncoGene_698 = find(gene_label(:,4) == 1 | gene_label(:,4) == 3);
TSG_1018 = find(gene_label(:,4) == 2 | gene_label(:,4) == 3);
Oncogene_Vogelstein = find(gene_label(:,5) == 1 | gene_label(:,5) == 3);
TSG_Vogelstein = find(gene_label(:,5) == 2 | gene_label(:,5) == 3);

len_cancers = length(cancer_names);

genome_class_dat = zeros(len_cancers, 4);
OG_class_dat = zeros(len_cancers, 4);
TSG_class_dat = zeros(len_cancers, 4);
OGV_class_dat = zeros(len_cancers, 4);
TSGV_class_dat = zeros(len_cancers, 4);

class_names = {' none significant'; ' methylation significant';' mutatation significant';' both significant'};

% 1.For each cancer, load M+ score, entropy, mutation rate
for i = 1 : len_cancers
    cancer_name = char(cancer_names(i));
    out_path = strcat(gene_class_path ,cancer_name,'/');
    cancer_paths = {out_path};
    for k = 1: length(cancer_paths)
        if ~exist(char(cancer_paths(k)))
            mkdir(char(cancer_paths(k)));
        end
    end
    entropy = load(strcat(fpre, 'methy_entropy/merged_stage/', cancer_name, '_entropy.dat'));
    entropy_diff = entropy(2:end, 3) - entropy(2:end,2);
    mp_score = load(strcat(fpre, 'methy_pvalue/merged_stage/',cancer_name,'/',cancer_name, '_p_score.dat'));
    mutation_rate  = load(strcat(fpre,'snv_intermidiate/merged_stage/',cancer_name,'/',cancer_name,'_',stage, '_mutation_rate.txt'));
    % class 0
    none_significants = find(~(entropy_diff(:, 1) > 1 & mp_score(:, 4) > mp_score_threshold) & mutation_rate(:, 2) < mutation_rate_threshold);
    % class 1
    methy_significants = find(entropy_diff(:, 1) > 1 & mp_score(:, 4) > mp_score_threshold & mutation_rate(:, 2) < mutation_rate_threshold);
    % class 2
    mutation_significants = find(~(entropy_diff(:, 1) > 1 & mp_score(:, 4) > mp_score_threshold) & mutation_rate(:, 2) > mutation_rate_threshold);
    % class 3
    both_significants = find(entropy_diff(:, 1) > 1 & mp_score(:, 4) > mp_score_threshold & mutation_rate(:, 2) > mutation_rate_threshold);
    length_vector_genome = zeros(1,4);
    for j = 0:3
        fid = fopen(strcat(out_path,cancer_name, '_', 'genome_class_', num2str(j) ,'.dat'),'w');
        switch j
            case 0
                dat = none_significants;
            case 1
                dat = methy_significants;
            case 2
                dat = mutation_significants;
            case 3
                dat = both_significants;
        end
        fprintf(fid,'%d\n', dat);
        fclose(fid);
        length_vector_genome(1,j+1) = length(dat);
        genome_class_dat(i,j+1) = length(dat);
    end
    sample_cnt_matrix = [0:3;length_vector_genome];
    fid = fopen(strcat(out_path,'genome_gcounts_of_class.dat'),'w');
    fprintf(fid,'%d\t%d\n', sample_cnt_matrix);
    fclose(fid);
    
    for k =1: 4
        switch k
            case 1
                genelist = OncoGene_698;
            case 2
                genelist = TSG_1018;
            case 3
                genelist = Oncogene_Vogelstein;
            case 4
                genelist = TSG_Vogelstein;
        end
        length_vector = zeros(1,4);
        for j = 0:3
            fid = fopen(strcat(out_path,cancer_name, '_', char(gene_category_names(k)),'_class_', num2str(j) ,'.dat'),'w');
            switch j
                case 0
                    dat = intersect(none_significants,genelist);
                case 1
                    dat = intersect(methy_significants,genelist);
                case 2
                    dat = intersect(mutation_significants,genelist);
                case 3
                    dat = intersect(both_significants,genelist);
            end
            fprintf(fid,'%d\n', dat);
            fclose(fid);
            len_dat = length(dat);
            length_vector(1,j+1) = len_dat;
            switch k
                case 1
                    OG_class_dat(i,j+1) = len_dat;
                case 2
                    TSG_class_dat(i,j+1) = len_dat;
                case 3
                    OGV_class_dat(i,j+1) = len_dat;
                case 4
                    TSGV_class_dat(i,j+1) = len_dat;
            end
        end
        sample_cnt_matrix = [0:3;length_vector];
        fid = fopen(strcat(out_path,char(gene_category_names(k)),'_gcounts_of_class.dat'),'w');
        fprintf(fid,'%d\t%d\n', sample_cnt_matrix);
        fclose(fid);
    end
end
close all;
end
