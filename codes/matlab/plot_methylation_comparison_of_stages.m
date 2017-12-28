function plot_entropy_comparison_of_stages()
global entropy_base_dir
clc;
clear all;
close all;
entropy_base_dir = '../../data/intermediate_file/methy_intermidiate/merged_stage/';
figure_base_dir = '../../figures/methylation_comparison/';

if ~exist(figure_base_dir)
    mkdir(figure_base_dir);
end

cancers={'COAD'};
for i = 1: length(cancers)
    cancer_name = char(cancers(i));
    fig = figure(i);
    Mn = load(strcat(entropy_base_dir,cancer_name,'/', cancer_name,'_normal_methy_dat.dat'));
    M1 = load(strcat(entropy_base_dir,cancer_name,'/', cancer_name,'_i_methy_dat.dat'));
    M2 = load(strcat(entropy_base_dir,cancer_name,'/', cancer_name,'_ii_methy_dat.dat'));
    M3 = load(strcat(entropy_base_dir,cancer_name,'/', cancer_name,'_iii_methy_dat.dat'));
    M4 = load(strcat(entropy_base_dir,cancer_name,'/', cancer_name,'_iv_methy_dat.dat'));
    
    Mn0 = Mn(2:end,2:end);
    M10 = M1(2:end,2:end);
    M20 = M2(2:end,2:end);
    M30 = M3(2:end,2:end);
    M40 = M4(2:end,2:end);
    
    Mn0 = mean(Mn0,2);
    M10 = mean(M10,2);
    M20 = mean(M20,2);
    M30 = mean(M30,2);
    M40 = mean(M40,2);
    
    subplot(2,2,1);
    hold on;
    plot(Mn0(:,1),M10(:,1),'.');
    plot([0 1], [0 1],'color','r');
    xlim([0 1]);
    ylim([0 1]);
    box on;
    xlabel('normal','FontSize',12,'FontWeight','bold');
    ylabel('i','FontSize',12,'FontWeight','bold');
    title('Methy of normal and i');
    
    subplot(2,2,2);
    hold on;
    plot(M10(:,1),M20(:,1),'.');
    plot([0 1], [0 1],'color','r');
    xlim([0 1]);
    ylim([0 1]);
    box on;
    xlabel('i','FontSize',12,'FontWeight','bold');
    ylabel('ii','FontSize',12,'FontWeight','bold');
    title('Methy of i and ii');
    
    subplot(2,2,3);
    hold on;
    plot(M10(:,1),M30(:,1),'.');
    plot([0 1], [0 1],'color','r');
    xlim([0 1]);
    ylim([0 1]);
    box on;
    xlabel('i','FontSize',12,'FontWeight','bold');
    ylabel('iii','FontSize',12,'FontWeight','bold');
    title('Methy of i and iii');
    
    subplot(2,2,4);
    hold on;
    plot(M10(:,1),M40(:,1),'.');
    plot([0 1], [0 1],'color','r');
    xlim([0 1]);
    ylim([0 1]);
    box on;
    xlabel('i','FontSize',12,'FontWeight','bold');
    ylabel('iv','FontSize',12,'FontWeight','bold');
    title('Methy of i and iv');
    exportfig(fig,strcat(figure_base_dir,cancer_name,'.eps'),'color','cmyk');
    close all;
end
end