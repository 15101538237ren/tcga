function gene_classification()
clc;
clear all;
close all;
global stage

stage = 'i';
cancer_names = {'BRCA';'COAD'; 'LIHC'; 'KIRC'; 'KIRP'; 'LUAD'; 'LUSC'; 'THCA'};
fpre = '../../data/intermediate_file/';
gene_class_path = strcat(fpre, 'gene_classification/');
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

% 1.For each cancer, load M+ score, entropy, mutation rate
for i = 1 : len_cancers
    cancer_name = char(cancer_names(i));
    out_path = strcat(fpre, 'gene_classification/',cancer_name,'/');
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
    none_significants = find(~(entropy_diff(:, 1) > 1 & mp_score(:, 4) > 0.3) & mutation_rate(:, 2) < 0.2);
    % class 1
    methy_significants = find(entropy_diff(:, 1) > 1 & mp_score(:, 4) > 0.3 & mutation_rate(:, 2) < 0.2);
    % class 2
    mutation_significants = find(~(entropy_diff(:, 1) > 1 & mp_score(:, 4) > 0.3) & mutation_rate(:, 2) > 0.2);
    % class 3
    both_significants = find(entropy_diff(:, 1) > 1 & mp_score(:, 4) > 0.3 & mutation_rate(:, 2) > 0.2);
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
            if isempty(dat)
                len_dat = 0;
            else
                len_dat = length(dat);
            end
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

nrow = 2;
ncol = 4;
for k = 1 : 5
    switch k
        case 1
            dat = OG_class_dat;
        case 2
            dat = TSG_class_dat;
        case 3
            dat = OGV_class_dat;
        case 4
            dat = TSGV_class_dat;
        case 5
            dat = genome_class_dat;
    end
    
    fig = figure(k);
    hold on;
    % loop cancer
    for i = 1 : 8
        subplot(nrow, ncol, i);
        dat_val = dat(i,:);
        pie(dat_val);
        colormap([
                1,0,0 %第一个是红的
                0,1,0 %第二个是绿的
                0,0,1 %第三个是蓝的
                0.98,0.5,0.04 %第四个是橙的
                ]);
        title(char(cancer_names(i)));
    end
    fig_path = strcat(figure_path, char(gene_category_names(k)),'_pie.eps');
    disp(fig_path);
    exportfig(fig,fig_path,'color','cmyk','fontmode','fixed','fontsize',10);
end

for k = 1 : 5
    switch k
        case 1
            dat = OG_class_dat;
        case 2
            dat = TSG_class_dat;
        case 3
            dat = OGV_class_dat;
        case 4
            dat = TSGV_class_dat;
        case 5
            dat = genome_class_dat;
    end
    cancer_dat = zeros(len_cancers + 1, 5);
    cancer_dat(1, 2 : 5) = 1 : 4;
    cancer_dat(2:end, 1) = 1 : len_cancers;
    cancer_dat(2:end, 2:end) = dat;
    cancer_dat = cancer_dat';
    out_dat_fp = strcat(class_by_gene_category_path, char(gene_category_names(k)),'.dat');
    fid = fopen(out_dat_fp,'w');
    fprintf(fid,'%d\t%d\t%d\t%d\t%d\n', cancer_dat);
    fclose(fid);
end

for j = 1 : 8
    category_dat = zeros(6, 5);
    category_dat(1, 2 : 5) = 1 : 4;
    category_dat(2:end, 1) = 1 : 5;
    for k = 1 : 5
        switch k
            case 1
                dat = OG_class_dat;
            case 2
                dat = TSG_class_dat;
            case 3
                dat = OGV_class_dat;
            case 4
                dat = TSGV_class_dat;
            case 5
                dat = genome_class_dat;
        end
        category_dat(k + 1, 2 : end) = dat(j,:);
    end
    category_dat = category_dat';
    out_dat_fp = strcat(class_by_cancer_path, char(cancer_names(j)),'.dat');
    fid = fopen(out_dat_fp,'w');
    fprintf(fid,'%d\t%d\t%d\t%d\t%d\n', category_dat);
    fclose(fid);
end
close all;
end