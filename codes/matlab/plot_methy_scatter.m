function plot_methy_scatter()
clc;
clear all;
close all;

methy_data_dir= '../../data/intermediate_file/common_patients_data/';
figure_dir = '../../figures/methy_scatter/';
if ~exist(figure_dir)
    mkdir(figure_dir);
end

cancer_name = 'COAD';
gene_name = 'APC';

fig=figure(1);
clf();
normal_dat= load(strcat(methy_data_dir,'/', cancer_name, '/normal','/promoter_mean_methy.tsv'));
i_dat= load(strcat(methy_data_dir,'/', cancer_name, '/i','/promoter_mean_methy.tsv'));

subplot(2,2,2);
plot(normal_dat(:, 1), normal_dat(:, 2),'b.');
hold on;
plot(i_dat(:, 1), i_dat(:, 2),'r.');

set(gca, 'YTick', [0 0.3 0.6 0.9]);
set(gca, 'XTick', [1 2]);
set(gca,'XTicklabel',{'normal';'i'});
axis([0 2.5 0.0 1.0]);
title(cancer_name);
fig_save_path = strcat(figure_dir, gene_name, '.eps');
exportfig(fig, fig_save_path , 'FontMode', 'fixed', 'color', 'cmyk','FontSize', 12,'Resolution',300,'LineWidth',0.5);
close all;
   

