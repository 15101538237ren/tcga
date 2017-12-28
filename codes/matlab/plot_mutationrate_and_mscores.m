function plot_mutationrate_and_mscores()
global base_path
clc;
clear all;
close all;
base_path = './';
mutation_base_dir = '../../data/intermediate_file/snv_intermidiate/merged_stage/';
pval_base_dir = '../../data/intermediate_file/methy_pvalue/merged_stage/';
figure_base_dir = '../../figures/mutation_rate_and_mscore/';
L=load('../../global_files/gene_label.dat');
Onco=1;
Tsg=2;
Both=3;
Onco_Vogel_indexs = find(L(:,5)==Onco);
TSG_Vogel_indexs = find(L(:,5)==Tsg);
match_gene_name(Onco_Vogel_indexs,'Onco_TSG_Vogel.ind');
match_gene_name(TSG_Vogel_indexs,'TSG_Vogel.ind');

cancers={'COAD'};
gene_types = {'significant_genes';'TSG_Vogel';'Onco_TSG_Vogel';'significant_no_genes'};
gene_type_names = {'Significant TSG & Onco Genes of Zhao list';'Vogelstein TSG list';' Vogelstein Onco Genes list'; 'Significant Other Genes'};
plot_region = [0.3 0.1 0.4 0.85];
mid = 0.5;
dis_l = 0.035;
dis_r = 0.02;
for i = 1: length(gene_types)
    gtname = char(gene_types(i));
    [order, gidxs, gene_names] = textread(strcat(gtname,'.ind'),'%d\t%d\t%s');
    len_gidxs = length(gidxs);
    for j = 1: length(cancers)
        cancer_name = char(cancers(j));
        figdir = strcat(figure_base_dir, cancer_name,'/');
        if ~exist(figdir)
            mkdir(figdir);
        end
        fig = figure(j);
        mp_score = load(strcat(pval_base_dir,cancer_name,'/', cancer_name,'_p_score.dat'));
        mn_score = load(strcat(pval_base_dir,cancer_name,'/', cancer_name,'_n_score.dat'));
        mut_rate = load(strcat(mutation_base_dir,cancer_name,'/', cancer_name,'_i_mutation_rate.txt'));

        mp_score = mp_score(:, 4);
        mn_score = mn_score(:, 4);
        mut_rate = mut_rate(:, 2);

        axes('position',plot_region);
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
            if (mut_rate_of_gene > 0.1 && i ~= length(gene_types))
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
            if (mpscore_of_gene > 0.3 && i ~= length(gene_types))
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
            if (mnscore_of_gene > 0.2 && i ~= length(gene_types))
                hold on;
                text(mn_xr + dis_r,(mn_yu + mn_yd)/2.0, gene_name,'FontSize',6,'FontWeight','bold');
            end
        end
        legend([p_mut,p_mp,p_mn],'Mutation Rate','M^+Score','M^-Score','Location','bestoutside');
        xlim([0 1]);
        ylim([bellow_y 1]);
        set(gca,'ytick',[]);
        set(gca,'xtick',[0,0.25,0.5,0.75,1]);
        set(gca,'xticklabel',{'1','0.5','0','0.5','1'});
        box off;
        axis off;
        title(char(gene_type_names(i)),'FontWeight','normal')
        exportfig(fig,strcat(figdir, gtname,'.eps'),'color','cmyk');
        close all;
    end
end
end