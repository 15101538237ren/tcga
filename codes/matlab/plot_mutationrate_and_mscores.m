function plot_mutationrate_and_mscores()
clc;
clear all;
close all;
mutation_base_dir = '../../data/intermediate_file/snv_intermidiate/merged_stage/';
pval_base_dir = '../../data/intermediate_file/methy_pvalue/merged_stage/';
figure_base_dir = '../../figures/mutation_rate_and_mscore/';

if ~exist(figure_base_dir)
    mkdir(figure_base_dir);
end

cancers={'COAD'};
gene_types = {'significant_genes';'significant_no_genes'};
plot_region = [0.3 0.1 0.4 0.85];
mid = 0.5;
for i = 1: length(gene_types)
    gtname = char(gene_types(i));
    [order, gidxs, gene_names] = textread(strcat(gtname,'.ind'),'%d\t%d\t%s');
    len_gidxs = length(gidxs);
    for j = 1: length(cancers)
        cancer_name = char(cancers(j));
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
        mut_height = 1.0 / (2 * len_gidxs - 1); % the height of each mutation of a gene
        mscore_height = 1.0 / (3 * len_gidxs - 1); % the height of each score
        % plot mutation rate bars
        for k = 1: len_gidxs
            gidx = gidxs(k);
            mut_rate_of_gene = mut_rate(gidx);
            xl = mid - (mut_rate_of_gene / 2.0);
            xr = mid;
            yu = 1.0 - (2*k - 2) * mut_height;
            yd = 1.0 - (2*k - 1) * mut_height;
            fill([xl, xl, xr, xr],[yd, yu, yu, yd],'r','linestyle','none');
        end
        
        % plot m^+score bars
        for k = 1: len_gidxs
            gidx = gidxs(k);
            mpscore_of_gene = mp_score(gidx);
            xl = mid - (mut_rate_of_gene / 2.0);
            xr = mid;
            yu = 1.0 - (2*k - 2) * mut_height;
            yd = 1.0 - (2*k - 1) * mut_height;
            fill([xl, xl, xr, xr],[yd, yu, yu, yd],'b','linestyle','none');
        end
        
        % plot m^-score bars
        for k = 1: len_gidxs
            gidx = gidxs(k);
            mnscore_of_gene = mn_score(gidx);
            xl = mid - (mut_rate_of_gene / 2.0);
            xr = mid;
            yu = 1.0 - (2*k - 2) * mut_height;
            yd = 1.0 - (2*k - 1) * mut_height;
            fill([xl, xl, xr, xr],[yd, yu, yu, yd],'g','linestyle','none');
        end
        
        xlim([0 1]);
        ylim([0 1]);
        hold on;
    end
end
end