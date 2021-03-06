function plot_methy_sample_entropy_distribution()
clc;
clear all;
close all;

sample_entropy_dir= '../../data/intermediate_file/methy_sample_entropy/merged_stage/';
cancer_names = {'BRCA'; 'COAD'; 'KIRC'; 'KIRP'; 'LIHC'; 'LUAD'; 'LUSC'; 'THCA'};
stage_names = {'n';'i';'ii';'iii';'iv'};
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
    sample_entropy= load(strcat(sample_entropy_dir, cancer_name, '_matlab_sample_entropy.dat'));
    plot(sample_entropy(:, 1), sample_entropy(:, 2),'b.');
    title(cancer_name);
    box on;
    set(gca, 'XTick', [1 2 3 4 5]);
    set(gca,'XTicklabel', stage_names);
    axis([0 5.5 min_entropy max_entropy]);
end
fig_save_path = strcat(figure_dir, 'sample_entropy_scatter.eps');
exportfig(fig, fig_save_path , 'FontMode', 'fixed', 'color', 'cmyk','width',8, 'height', 4, 'FontSize', 12,'Resolution',300,'LineWidth',0.5);
close all;