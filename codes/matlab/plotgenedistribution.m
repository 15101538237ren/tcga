function plotgenedistribution()
clc;
clear all;
close all;
methy_data_dir = '../../data/intermediate_file/methy_matlab_data/';
fig_dir = '../../figures/methy_scatter/TSG_Vogelstein/';

if ~exist(fig_dir) 
    mkdir(fig_dir)
end 

cancer_names={'BRCA';'COAD';'LIHC';'LUAD';'LUSC'}; %
cancer_markers = {'r-s';'g-p';'b-o';'k-x';'c-*'};
mean_std_names = {'mean';'std'};
stage_names = {'normal';'i';'ii';'iii';'iv'};%;'x'
gene_idx_filepath = '../../global_files/TSG_Vogelstein.txt';
gene_names = textread(gene_idx_filepath,'%s');

rows = 5;
cols = 7;
dx=0.01;
fig_counter = 0;

alw = 0.75;    % AxesLineWidth
fsz = 16;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize

for i=1: length(gene_names)
    gene_name = char(gene_names(i));
    fig_counter = fig_counter + 1;
    fig=figure(fig_counter);
    clf();
    hold on;
    % 第一张图, 每个基因一个图, 每个癌症一个子图, 画散点和箱线图
    for j=1:length(cancer_names)
        cancer_name = char(cancer_names(j));
        data1_path= strcat(methy_data_dir,'merged_stage/',cancer_name,'/', gene_name, '_xy_', cancer_name,'.dat');
        data1 = load(data1_path);
        subplot(rows,cols,cols*(j-1)+1);
        plot(data1(:, 1), data1(:, 2),'b.');
        set(gca, 'YTick', [0 0.3 0.6 0.9]);
        set(gca, 'XTick', [1 2 3 4 5 6]);
        set(gca,'XTicklabel',{'n';'i';'ii';'iii';'iv';'x'});
        axis([0 6.5 0.0 1.0]);
        title(cancer_name);
        set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
    end
    hold on;
    for j=1:length(cancer_names)
        cancer_name = char(cancer_names(j));
        data2_path= strcat(methy_data_dir,'merged_stage/',cancer_name,'/', gene_name, '_all_stage_', cancer_name,'.dat');
        data2 = load(data2_path);
        subplot(rows,cols,cols*(j-1)+2);
        [u,x]=hist(data2, 0:dx:1);
        plot(x,u/(dx*sum(u)));
        box on;
        set(gca, 'XTick', [0 0.3 0.6 0.9]);
        %axis([0 6.5 0.0 1.0]);
        title('Distribution');
        set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
    end
    
    for j=1:length(stage_names)
        hold on;
        stage_name = char(stage_names(j));
        for k=1:length(cancer_names)
            cancer_name = char(cancer_names(k));
            data3_path= strcat(methy_data_dir,'merged_stage/',cancer_name,'/', gene_name, '_', stage_name, '_', cancer_name,'.dat');
            data3 = load(data3_path);
            subplot(rows,cols,cols*(k-1)+2+j);
            hold on;
            [u,x]=hist(data3,0:dx:1);
            bar(x,u/(dx*sum(u)));
            set(gca, 'XTick', [0 0.3 0.6 0.9]);
            xlim([0 1]);
            ylim([0 60]);
            title(stage_name);
            box on;
            set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
        end
    end
     set(fig,'paperunits','centimeter')
     set(fig,'papersize',[80,50])
     set(fig,'paperposition',[0 0 80 50]);
     fig_save_path = strcat(fig_dir, gene_name,'.pdf');
     print(fig,fig_save_path,'-dpdf','-opengl','-r300');
%     exportfig(fig, fig_save_path , 'FontMode', 'fixed', 'color', 'cmyk','width',10, 'height',6,'FontSize', 12,'Resolution',300,'LineWidth',0.5);
    close all;
    %
end

