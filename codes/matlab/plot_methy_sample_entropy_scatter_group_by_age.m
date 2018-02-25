function plot_methy_sample_entropy_scatter_group_by_age()
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
min_entropy = 1.0;
max_entropy = 4.0;

lengend_names = {'40-', '40-60', '60+'};
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
        x_young = zeros(length(young), 1) + rand(length(young), 1)*0.4 + 0.8;
        x_mid = zeros(length(mid), 1) + rand(length(mid), 1)*0.4 + 1.8;
        x_high = zeros(length(high), 1) + rand(length(high), 1)*0.4 + 2.8;
        
        hold on;
        plot(x_young, sample_entropy(young, 2),'r.');
        plot(x_mid, sample_entropy(mid, 2),'g.');
        plot(x_high, sample_entropy(high, 2),'b.');
        
        if j == 1
           legend(lengend_names, 'Box','off');
        end
        title(cancer_name);
        box on;
        axis([0 3.5 min_entropy max_entropy]);
        set(gca, 'XTick', []);
        if i==1 || i == 5
            ylabel('entropy');
        end
    end
    fig_save_path = strcat(figure_dir, 'sample_entropy_scatter_group_by_age_', stage_name,'.eps');
    exportfig(fig, fig_save_path , 'FontMode', 'fixed', 'color', 'cmyk','width',8, 'height', 4, 'FontSize', 12,'Resolution',300,'LineWidth',0.5);
end
close all;
