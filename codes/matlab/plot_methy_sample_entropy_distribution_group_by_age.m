function plot_methy_sample_entropy_distribution_group_by_age()
clc;
clear all;
close all;

sample_entropy_dir= '../../data/intermediate_file/methy_sample_entropy/merged_stage/';
age_data_dir= '../../data/intermediate_file/methy_intermidiate/merged_stage/';
cancer_names = {'BRCA'; 'COAD'; 'KIRC'; 'KIRP'; 'LIHC'; 'LUAD'; 'LUSC'; 'THCA'};
stage_names = {'normal';'i';'ii';'iii';'iv'};
ncol = 4;
nrow = 2;
figure_dir = '../../figures/sample_entropy/';
if ~exist(figure_dir)
    mkdir(figure_dir);
end

dx=0.04;
min_entropy = 1.5;
max_entropy = 3.0;
for i = 1 : length(stage_names)
    fig=figure(i);
    clf();
    stage_name = char(stage_names(i));
    
    for j= 1 : length(cancer_names)
        cancer_name = char(cancer_names(j));
        subplot(nrow, ncol, j);
        sample_entropy= load(strcat(sample_entropy_dir, cancer_name, '_', stage_name ,'_sample_entropy.dat'));
        age = load(strcat(age_data_dir, cancer_name, '/', cancer_name, '_', stage_name ,'_ages.txt'));
        
        young = find(age(:, 2) > 0 & age(:, 2) < 40);
        mid = find(age(:, 2) >= 40 & age(:, 2) < 60);
        high = find(age(:, 2) >= 60);
        
        hold on;
        [u, x] = hist(sample_entropy(young, 2), min_entropy: dx: max_entropy);
        plot(x, u/(dx*sum(u)),'-r');
        
        [u, x] = hist(sample_entropy(mid, 2), min_entropy: dx: max_entropy);
        plot(x, u/(dx*sum(u)),'-g');
        
        [u, x] = hist(sample_entropy(high, 2), min_entropy: dx: max_entropy);
        plot(x, u/(dx*sum(u)),'-b');
        if j == 1
           legend({'< 40', '40~60', '> 60'}, 'Box','off');
        end
        title(cancer_name);
        box on;
        axis([min_entropy max_entropy 0.0 20.0]);
        if i > 4
            xlabel('entropy');
        end
        if i==1 || i == 5
            ylabel('cases (%)');
        end
    end
    fig_save_path = strcat(figure_dir, 'sample_entropy_distribution_group_by_age_', stage_name,'.eps');
    exportfig(fig, fig_save_path , 'FontMode', 'fixed', 'color', 'cmyk','width',8, 'height', 4, 'FontSize', 12,'Resolution',300,'LineWidth',0.5);
end
close all;
