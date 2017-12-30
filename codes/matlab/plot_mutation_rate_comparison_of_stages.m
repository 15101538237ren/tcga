function plot_entropy_comparison_of_stages()
clc;
clear all;
close all;
mutation_base_dir = '../../data/intermediate_file/snv_intermidiate/merged_stage/';
fig_base = '../../figures/mutation_rate_variation/';

if ~exist(fig_base)
    mkdir(fig_base);
end

cancers={'COAD'};
gidxs = [831,17125,8406,12171];
stages = {'i';'ii';'iii';'iv'};
gene_idx_filepath = '../../global_files/gene_idx.txt';
[all_gidxs, gene_names] = textread(gene_idx_filepath,'%d\t%s');
c=['y','b','g','r'];
%L=load('../../global_files/gene_label.dat');

for i = 1: length(cancers)
    cancer_name = char(cancers(i));
    M = zeros(length(all_gidxs),length(stages));
    for j = 1 : length(stages)
        stname = char(stages(j));
        mut_dat=load(strcat(mutation_base_dir, cancer_name,'/',cancer_name,'_',stname,'_mutation_rate.txt'));
        M(:, j) = mut_dat(:, 2)';
    end
    fig = figure(i);
    
    subplot(2,2,1);
    hold on;
    plot(M(:,1),M(:,2),'.');
    plot([0 1], [0 1],'color','r');
    xlim([0 1]);
    ylim([0 1]);
    box on;
    xlabel('i','FontSize',12,'FontWeight','bold');
    ylabel('ii','FontSize',12,'FontWeight','bold');
    title('Mutation Rate of i and ii');
    
    subplot(2,2,2);
    hold on;
    plot(M(:,1),M(:,3),'.');
    plot([0 1], [0 1],'color','r');
    xlim([0 1]);
    ylim([0 1]);
    box on;
    xlabel('i','FontSize',12,'FontWeight','bold');
    ylabel('iii','FontSize',12,'FontWeight','bold');
    title('Mutation Rate of i and iii');
    
    subplot(2,2,3);
    hold on;
    plot(M(:,1),M(:,4),'.');
    plot([0 1], [0 1],'color','r');
    xlim([0 1]);
    ylim([0 1]);
    box on;
    xlabel('i','FontSize',12,'FontWeight','bold');
    ylabel('iv','FontSize',12,'FontWeight','bold');
    title('Mutation Rate of i and iv');
    
    subplot(2,2,4);
    hold on;
    A = 1:length(stages);
    for k=1:length(gidxs)
        B = M(k,:);
        plot(A,B,strcat(c(k),'--o'),'markerfacecolor',c(k),'markersize',5)
    end
    legend('APC','TP53','KRAS','PIK3CA');
    title('Mutation Rate Variation of Key Genes');
    xlabel('Stage','FontSize',12,'FontWeight','bold');
    xlim([1 4.2]);
    ylim([0 0.2]);
    box on;
    set(gca,'xtick',[1,2,3,4]);
    set(gca,'ytick',[0.05,0.1,0.15,0.2]);
    set(gca,'xticklabel',{'i','ii','iii','iv'});
    ylabel('Mutation Rate','FontSize',12,'FontWeight','bold');
    exportfig(fig,strcat(fig_base,cancer_name,'.eps'),'color','cmyk');
    close all;
end