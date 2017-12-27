function plot_entropy_comparison_of_stages()
global entropy_base_dir
clc;
clear all;
close all;
entropy_base_dir = '../../data/intermediate_file/methy_entropy/merged_stage/';
figure_base_dir = '../../figures/entropy_comparison/';

if ~exist(figure_base_dir)
    mkdir(figure_base_dir);
end

cancers={'COAD'};
for i = 1: length(cancers)
    cancer_name = char(cancers(i));
    fig = figure(i);
    E = load(strcat(entropy_base_dir,cancer_name,'_entropy.dat'));
    E0=E(2:end,2:end);
    
    subplot(2,2,1);
    hold on;
    plot(E0(:,1),E0(:,2),'.');
    plot([0 4], [0 4],'color','r');
    plot([0 4],[1 5],'color','r','linestyle','--');
    xlim([0 3]);
    ylim([0 3.5]);
    box on;
    xlabel('normal');
    ylabel('i');
    title('Entropy of normal and i');
    
    subplot(2,2,2);
    hold on;
    plot(E0(:,2),E0(:,3),'.');
    plot([0 4], [0 4],'color','r');
    plot([0 4],[1 5],'color','r','linestyle','--');
    xlim([0 3]);
    ylim([0 3.5]);
    box on;
    xlabel('i');
    ylabel('ii');
    title('Entropy of i and ii');
    
    subplot(2,2,3);
    hold on;
    plot(E0(:,2),E0(:,4),'.');
    plot([0 4], [0 4],'color','r');
    plot([0 4],[1 5],'color','r','linestyle','--');
    xlim([0 3]);
    ylim([0 3.5]);
    box on;
    xlabel('i');
    ylabel('iii');
    title('Entropy of i and iii');
    
    subplot(2,2,4);
    hold on;
    plot(E0(:,2),E0(:,5),'.');
    plot([0 4], [0 4],'color','r');
    plot([0 4],[1 5],'color','r','linestyle','--');
    xlim([0 3]);
    ylim([0 3.5]);
    box on;
    xlabel('i');
    ylabel('iv');
    title('Entropy of i and iv');
    exportfig(fig,strcat(figure_base_dir,cancer_name,'.eps'),'color','cmyk');
    close all;
end
end