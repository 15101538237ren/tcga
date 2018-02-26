function plot_scatter_of_sample_entropy_and_age()
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

fig=figure(1);
clf();
min_entropy = 1.0;
max_entropy = 4.0;

for i= 1 : length(cancer_names)
    cancer_name = char(cancer_names(i));
    subplot(nrow, ncol, i);
    hold on;
    for j = 1 : length(stage_names)
        stage_name = char(stage_names(j));
        sample_entropy= load(strcat(sample_entropy_dir, cancer_name, '_', stage_name ,'_sample_entropy.dat'));
        age = load(strcat(age_data_dir, cancer_name, '/', cancer_name, '_', stage_name ,'_ages.txt'));
        if j == 1
            point_style = 'b.';
        else
            point_style = 'r.';
        end
        plot(age(:, 2), sample_entropy(:, 2), point_style);
    end
    
    title(cancer_name);
    box on;
    axis([0 100 min_entropy max_entropy]);
    set(gca, 'XTick', [20 40 60 80 100]);
    if i > 4
        xlabel('age');
    end
    
    if i==1 || i==5
        ylabel('entropy');
    end
    if i ==1
        legend({'normal','cancer'},'Box','off');
    end
end
fig_save_path = strcat(figure_dir, 'scatter_of_sample_entropy_and_age.eps');
exportfig(fig, fig_save_path , 'FontMode', 'fixed', 'color', 'cmyk','width',8, 'height', 4, 'FontSize', 12,'Resolution',300,'LineWidth',0.5);
close all;