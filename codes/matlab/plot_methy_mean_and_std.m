function plot_methy_mean_and_std()
% This function plot the scatter of the methylation means and its standard variations of different stages for all genes

clc;
clear all;
close all;

% 1. change work directory
mean_std_data_path = '../../data/intermediate_file/methy_mean_std/merged_stage/';
figure_path = '../../figures/methy_mean_and_std/';
cancer_names = {'BRCA'; 'COAD'; 'KIRC'; 'KIRP'; 'LIHC'; 'LUAD'; 'LUSC'; 'THCA'};
size_cancer_names = size(cancer_names);

stage_names = {'normal';'i'}; %;'ii';'iii';'iv'
size_stage = size(stage_names);

stage_markers = {'g.';'r.'};%;'b.';'k.';'c.'

% counter for figures
fig_counter = 0;

for i=1: size_cancer_names(1)
    cancer_name = char(cancer_names(i));
    fig_counter = fig_counter + 1;
    fig = figure(fig_counter);
    clf();
    mean_data_raw = load(fullfile(mean_std_data_path, cancer_name,[cancer_name, '_mean.dat']));
    std_data_raw = load(fullfile(mean_std_data_path, cancer_name,[cancer_name, '_std.dat']));
    mean_data = mean_data_raw(2:end,2:end);
    std_data = std_data_raw(2:end,2:end);
    hold on;
    %iterate stage
    for j= 1: size_stage(1)
        plot(mean_data(:,j),(std_data(:,j).^2)./(mean_data(:,j).^2),char(stage_markers(j))); 
    end
    axis([0 1.0 0.0 5.0]);
    title(cancer_name);
    xlabel('\mu');
    ylabel('\sigma^2/\mu^2');
    fig_save_path = [figure_path, cancer_name,'.eps'];
    exportfig(fig, fig_save_path , 'FontMode', 'fixed', 'color', 'cmyk', 'FontSize', 12,'Resolution',300,'LineWidth',0.5);
end
close all;
end

