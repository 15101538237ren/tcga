function plot_methy_sample_entropy_scatter()
clc;
clear all;
close all;

sample_entropy_dir= '../../data/intermediate_file/methy_sample_entropy/merged_stage/';
cancer_names = {'BRCA'; 'COAD'; 'KIRC'; 'KIRP'; 'LIHC'; 'LUAD'; 'LUSC'; 'THCA'};

figure_dir = '../../figures/sample_entropy/';
if ~exist(figure_dir)
    mkdir(figure_dir);
end
for i= 1 : length(cancer_names)
    cancer_name = char(cancer_names(i));
    fig=figure(i);
    clf();
    sample_entropy= load(strcat(sample_entropy_dir, cancer_name, '_matlab_sample_entropy.dat'));
    subplot(2,2,2);
    plot(sample_entropy(:, 1), sample_entropy(:, 2),'b.');
    set(gca, 'XTick', [1 2 3 4 5]);
    set(gca,'XTicklabel',{'n';'i';'ii';'iii';'iv'});
    axis([0 5.5 1.0 4.0]);
    title(['Sample Entropy of ', cancer_name]);
    fig_save_path = strcat(figure_dir, cancer_name, '.eps');
    exportfig(fig, fig_save_path , 'FontMode', 'fixed', 'color', 'cmyk','FontSize', 12,'Resolution',300,'LineWidth',0.5);
end
close all;
   

