function plot_mutation_classification_single_i()
global sample_base_dir;
global sample_num;
global cancer_name;
global gene_name;
global gene_id;
global tss_len;
global figdir;
global gene_body_len;
global cancer_stage;
global out_methy_file;
global out_ins_file;
global out_del_file;
global out_snp_file;
global normal_mean_methy_file;
%clc;
%clear;
%close all;

sub_plot_num = 4;
fig_num = floor(sample_num/sub_plot_num);


for k=1:fig_num
    fig = figure(k);
    clf();
    subplot(sub_plot_num+1,1,1);
    hold on;
    normal_methy_list = load(normal_mean_methy_file);
    J = find(normal_methy_list(:,1) < 0);
    normal_methy_list(J,1) = normal_methy_list(J,1) / tss_len;
    J = find(normal_methy_list(:,1) >= 0);
    normal_methy_list(J,1) = normal_methy_list(J,1) / gene_body_len;
    bar(normal_methy_list(:,1), normal_methy_list(:,2), 'k', 'LineWidth',1.5);
    title(strcat(gene_name, ' normal average beta-value'));
    xlim([-0.2 1]);
    ylim([0 1]);

    for i=1:sub_plot_num
        sample_id=sub_plot_num*(k-1)+i;
        subplot(sub_plot_num+1,1,i+1);
        hold on;
        title(strcat(num2str(sample_id)));
        xlim([-0.2 1]);
        ylim([0 1]);
        py.get_gene_data.get_sample_methy(sample_base_dir, cancer_name, cancer_stage, sample_id, gene_id, out_methy_file);
        py.get_gene_data.get_sample_mutation(sample_base_dir, cancer_name, cancer_stage, sample_id, gene_id);
        methy_list = load(out_methy_file);
        ins_list = load(out_ins_file);
        del_list = load(out_del_file);
        snp_list = load(out_snp_file);
        
        J = find(methy_list(:,1) < 0);
        methy_list(J,1) = methy_list(J,1) / tss_len;
        J = find(methy_list(:,1) >= 0);
        methy_list(J,1) = methy_list(J,1) / gene_body_len;
        bar(methy_list(:,1), methy_list(:,2), 'k', 'LineWidth',1.5);
        hold on;
        
        ins_list(ins_list <= 0) = ins_list(ins_list <= 0) / tss_len;
        ins_list(ins_list > 0) = ins_list(ins_list > 0) / gene_body_len;
        del_list(del_list <= 0) = del_list(del_list <= 0) / tss_len;
        del_list(del_list > 0) = del_list(del_list > 0) / gene_body_len;
        snp_list(snp_list <= 0) = snp_list(snp_list <= 0) / tss_len;
        snp_list(snp_list > 0) = snp_list(snp_list > 0) / gene_body_len;
        
        y1 = ones(1, size(ins_list, 2)) * 0.05;
        y2 = ones(1, size(del_list, 2)) * 0.2;
        y3 = ones(1, size(snp_list, 2)) * 0.35;
        yy3 = ones(1, size(snp_list, 2)) * 0.4;
        plot(ins_list, y1, '-v', 'MarkerSize',4, 'Markerfacecolor','r','markeredgecolor','r'); 
        plot(del_list, y2, '-^', 'MarkerSize',4, 'Markerfacecolor','b','markeredgecolor','b');
        plot(snp_list, y3, '-o', 'MarkerSize',4, 'Markerfacecolor','m','markeredgecolor','m'); 
        
    end
    print(fig,strcat(figdir, gene_name, '_', num2str(k), '_compare.pdf'),'-dpdf','-opengl');
end
close all;
delete('*.txt');
end