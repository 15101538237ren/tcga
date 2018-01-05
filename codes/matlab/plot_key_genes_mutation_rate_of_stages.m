function plot_mutation_rate_varation_among_stages()
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
gnames = {'APC';'TP53';'KRAS';'PIK3CA'};
for i=1:length(cancers)
    cancer_name = char(cancers(i));
    mut_rate_matrix = zeros(length(gidxs),length(stages));
    for j = 1 : length(stages)
        stname = char(stages(j));
        mut_dat=load(strcat(mutation_base_dir, cancer_name,'/',cancer_name,'_',stname,'_mutation_rate.txt'));
        mut_target = mut_dat(gidxs, 2);
        mut_rate_matrix(:, j) = mut_target';
    end
    
    fig = figure(i);
    clf();
    
    for k=1:length(gidxs)
        subplot(2,2,k);
        hold on;
        A = 1:length(stages);
        B = mut_rate_matrix(k,:);
        plot(A,B,'k--o','markerfacecolor','k','markersize',5)
        title(char(gnames(k)));
        xlabel('Stage','FontSize',12,'FontWeight','bold');
        xlim([0 4.5]);
        ylim([0 1]);
        box on;
        set(gca,'xtick',[1,2,3,4]);
        set(gca,'ytick',[0.2, 0.4, 0.6, 0.8, 1.0]);
        set(gca,'xticklabel',{'i','ii','iii','iv'});
        if (mod(k,2) == 1)
            ylabel('Mutation Rate','FontSize',12,'FontWeight','bold');
        end
    end
    exportfig(fig,strcat(fig_base, cancer_name,'_keygenes.eps'),'color','cmyk');
    close all;
end
end