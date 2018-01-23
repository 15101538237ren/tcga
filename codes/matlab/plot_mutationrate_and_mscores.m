function plot_mutationrate_and_mscores()
global base_path
clc;
clear all;
close all;
base_path = './';
mutation_base_dir = '../../data/intermediate_file/snv_intermidiate/merged_stage/';
pval_base_dir = '../../data/intermediate_file/methy_pvalue/merged_stage/';

% mutation_base_dir = 'I:/intermediate_file/snv_intermidiate/merged_stage/';
% pval_base_dir = 'I:/intermediate_file/methy_pvalue/merged_stage/';
figure_base_dir = '../../figures/mutation_rate_and_mscore/';
cancer_name='COAD';
figdir = strcat(figure_base_dir, cancer_name,'/');

if ~exist(figdir)
    mkdir(figdir);
end
mp_score_threshold = 0.8;
mutation_rate_threshold = 0.1;
fpre = '../../data/intermediate_file/';
sig_base_dir = strcat(fpre, 'gene_classification_mp_',num2str(mp_score_threshold),'_mut_',num2str(mutation_rate_threshold),'/');
gene_types = {'significant_genes';'significant_no_genes';'TSG_Vogel';'Onco_TSG_Vogel';};
gene_type_names = {'Significant TSG & Onco Genes of Zhao list';'Significant Other Genes';'Vogelstein TSG list';' Vogelstein Onco Genes list'};
fig_names = {'Significant_Genes';'Vogelstein_Genes'};
plot_region = [0.3 0.1 0.4 0.85];
mid = 0.5;
dis_l = 0.02;
dis_r = 0.02;

for i = 1: 2
    fig = figure(i);
    clf();
    for j = 1 : 2
        kk = (i - 1) * 2 + j ;
        if (kk > 2)
           dis_l = 0.025;
           dis_r = 0.02;
        end
        gtname = char(gene_types(kk));
        [order, gidxs, gene_names] = textread(strcat(sig_base_dir, gtname,'.ind'),'%d\t%d\t%s');
        len_gidxs = length(gidxs);
        mp_score = load(strcat(pval_base_dir,cancer_name,'/', cancer_name,'_p_score.dat'));
        mn_score = load(strcat(pval_base_dir,cancer_name,'/', cancer_name,'_n_score.dat'));
        mut_rate = load(strcat(mutation_base_dir,cancer_name,'/', cancer_name,'_i_mutation_rate.txt'));

        mp_score = mp_score(:, 4);
        mn_score = mn_score(:, 4);
        mut_rate = mut_rate(:, 2);
        subplot(1,2,j);
        %axes('position',plot_region);
        hold on;
        plot([mid mid],[0 1],'color','k');
        bellow_y = -0.01;
        bellow_yy = -0.04;
        plot([0 1],[bellow_y bellow_y],'color','k');
        xs = [0, 0.25,0.5,0.75,1];
        xtext = {'1';'0.5';'0';'0.5';'1'};
        fs = 8;
        for xi = 1: length(xs)
          text(xs(xi),bellow_yy,char(xtext(xi)),'FontSize',fs);
        end
        mut_height = 1.0 / (3 * len_gidxs - 1); % the height of each mutation of a gene
        mscore_height = 1.0 / (3 * len_gidxs - 1); % the height of each score

        for k = 1: len_gidxs
            gidx = gidxs(k);
            gene_name = char(gene_names(k));
            % plot mutation rate bars
            mut_rate_of_gene = mut_rate(gidx);
            mu_xl = mid - (mut_rate_of_gene / 2.0);
            mu_xr = mid;
            mu_yu = 1.0 - (3*k - 3) * mut_height;
            mu_yd = 1.0 - (3*k - 1) * mut_height;
            p_mut= fill([mu_xl, mu_xl, mu_xr, mu_xr],[mu_yd, mu_yu, mu_yu, mu_yd],'r','linestyle','none');
            if (mut_rate_of_gene > 0.1 && kk ~= 2)
                hold on;
                text(mu_xl-0.02-dis_l * length(gene_name),(mu_yu + mu_yd)/2.0,gene_name,'FontSize',6,'FontWeight','bold');
            end

            % plot m^+score bars
            mpscore_of_gene = mp_score(gidx);
            if mpscore_of_gene < 0
                mpscore_of_gene = 0;
            end
            mp_xl = mid;
            mp_xr = mid + (mpscore_of_gene / 2.0);
            mp_yu = 1.0 - (3*k - 3) * mut_height;
            mp_yd = 1.0 - (3*k - 2) * mut_height;
            p_mp = fill([mp_xl, mp_xl, mp_xr, mp_xr],[mp_yd, mp_yu, mp_yu, mp_yd],'b','linestyle','none');
            if (mpscore_of_gene > 0.3 && kk ~= 2)
                hold on;
                text(mp_xr + dis_r,(mp_yu + mp_yd)/2.0,gene_name,'FontSize',6,'FontWeight','bold');
            end

            % plot m^-score bars
            mnscore_of_gene = mn_score(gidx);
            if mnscore_of_gene < 0
                mnscore_of_gene = 0;
            end
            mn_xl = mid;
            mn_xr = mid + (mnscore_of_gene / 2.0);
            mn_yu = 1.0 - (3*k - 2) * mut_height;
            mn_yd = 1.0 - (3*k - 1) * mut_height;
            p_mn = fill([mn_xl, mn_xl, mn_xr, mn_xr],[mn_yd, mn_yu, mn_yu, mn_yd],'g','linestyle','none');
            if (mnscore_of_gene > 0.2 && kk ~= 2)
                hold on;
                text(mn_xr + dis_r,(mn_yu + mn_yd)/2.0, gene_name,'FontSize',6,'FontWeight','bold');
            end
        end
        if (i==1 && j == 2)
            legend([p_mut,p_mp,p_mn],'Mutation Rate','M^+Score','M^-Score','Location','northwestoutside');
        end
        xlim([0 1]);
        ylim([bellow_y 1]);
        set(gca,'ytick',[]);
        set(gca,'xtick',[0,0.25,0.5,0.75,1]);
        set(gca,'xticklabel',{'1','0.5','0','0.5','1'});
        box off;
        axis off;
        title(char(gene_type_names(kk)),'FontWeight','normal')
    end
    exportfig(fig,strcat(figdir, char(fig_names(i)),'.eps'),'color','cmyk');
    %print(fig,'a.pdf','-dpdf','-opengl');
    close all;
end
end