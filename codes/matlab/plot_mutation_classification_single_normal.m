function plot_mutation_classification_single_normal(cancer_name, gene_name, gene_id)
global sample_base_dir;
global normal_sample_num;
global tss_len;
global figdir;
global gene_body_len;
global out_normal_methy_file;
global normal_mean_methy_file;
global normal_stage;
%clc;
%clear;
%close all;

sub_plot_num = 4;
fig_num = floor(normal_sample_num/sub_plot_num);

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
        py.get_gene_data.get_sample_methy(sample_base_dir, cancer_name, normal_stage, sample_id, gene_id, out_normal_methy_file);
        methy_list = load(out_normal_methy_file);
        
        J = find(methy_list(:,1) < 0);
        methy_list(J,1) = methy_list(J,1) / tss_len;
        J = find(methy_list(:,1) >= 0);
        methy_list(J,1) = methy_list(J,1) / gene_body_len;
        bar(methy_list(:,1), methy_list(:,2), 'k', 'LineWidth',1.5);
        
    end
    print(fig,strcat(figdir, gene_name, '_', num2str(k), '_normal.pdf'),'-dpdf','-opengl');
end
close all;
%delete('*.txt');
end