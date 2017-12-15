function plotdistribution()
clc;
clear all;
close all;

methy_data_dir = '';

cancer_names={'COAD'};
size_cancer_names = size(cancer_names);

stage_names = {'normal';'i';'ii';'iii';'iv'};%;'x'
gene_idx_filepath = '../../global_files/gene_idx.txt';

[gene_indexs, gene_names] = textread(gene_idx_filepath,'%d\t%s');
target_gene_index = 5237;
gene_name = char(gene_names(target_gene_index));

dx=0.01;

fig=figure(1);
clf();

for j=1:size_cancer_names(1)
    cancer_name = char(cancer_names(j));
    data2_path= [methy_data_dir, gene_name, '_all_stage_', cancer_name,'.dat'];
    data2 = load(data2_path);
%     subplot(2,2,1);
    [u,x]=hist(data2, 0:dx:1);
%     plot(x,u/(dx*sum(u)));
%     box on;
%     set(gca, 'XTick', [0 0.3 0.6 0.9]);
    %axis([0 6.5 0.0 1.0]);
%     title(['(a) ',gene_name, ' Distribution']);
    
    subplot(2,2,2);
    histogram(data2,x,'FaceColor','blue');
    box on;
    set(gca, 'XTick', [0 0.3 0.6 0.9]);
    title([gene_name, ' Histogram']);
    fig_save_path = [cancer_name, '_', gene_name,'_histogram.eps'];
    exportfig(fig, fig_save_path , 'FontMode', 'fixed', 'color', 'cmyk','FontSize', 12,'Resolution',300,'LineWidth',0.5);
end
close all;

