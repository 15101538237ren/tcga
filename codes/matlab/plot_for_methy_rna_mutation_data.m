function plot_scatter_for_three_types_data()
clc;
clear all;
close all;

gene_idx_filepath = '../../global_files/gene_idx.txt';
gene_names = extract_gene_idx_and_genename(gene_idx_filepath);

target_gene_idxs = load('../../global_files/target_gene_idx.dat');
size_target_gene = size(target_gene_idxs);
intermediate_dir = '../../data/intermediate_file/';
is_merge_stage = 'merged_stage/';

data_dir_names = {'methy_intermidiate/';'rna_intermidiate/';'snv_intermidiate/'};

cancer_names = {'BRCA';'COAD'; 'LIHC'};%; ; 'KIRC'; 'KIRP'; 'LUAD'; 'LUSC'; 'THCA'
size_cancer_names = size(cancer_names);

stage_names = {'normal';'i';'ii';'iii';'iv'};
x_tick_max = 5.5;
x_tick_nums = [1 2,3 4 5];
x_tick_names = {'n';'i';'ii';'iii';'iv'};
size_stage = size(stage_names);
stage_markers = {'b.';'r.';'g.';'k.';'c.'};

fig_dir = '../../figures/methy_rna_mutation/';

data_types = {'Methylation';'mRNA';'SNV'};
size_data_types = size(data_types);
file_ends = {'_methy_dat.dat';'_tpm.dat';'_mutation_rate.txt'};

fig_counter = 0;
nrow = 5;
ncol = 3;
for j = 1 : size_target_gene(1)
    fig_counter = fig_counter + 1;
    fig=figure(fig_counter);
    clf();
    gene_name = char(gene_names(target_gene_idxs(j , 1)));
    fig_path = [fig_dir, gene_name,'.eps'];
    hold on;
    for i = 1 : size_cancer_names(1)
        cancer_name = char(cancer_names(i));
        %1. plot methylation beta-value scatter for each stage
        subplot(nrow, ncol, ncol*(i-1) + 1);
        for k = 1 : size_stage(1)
            stage_name = char(stage_names(k));
            methy_data_for_stage = load([intermediate_dir, char(data_dir_names(1)), is_merge_stage, cancer_name , '/', cancer_name, '_', stage_name, char(file_ends(1))]);
            methy_data_of_gene = methy_data_for_stage(target_gene_idxs(j , 1) + 1, :)';
            size_y_data = size(methy_data_of_gene);
            x_data = ones(size_y_data(1), 1) .* k  + (0.4 * rand(size_y_data(1), 1) - 0.2);
            plot(x_data(:, 1), methy_data_of_gene(:, 1), char(stage_markers(k)));
            hold on;
        end
        
        set(gca, 'YTick', [0 0.3 0.6 0.9]);
        set(gca, 'XTick', x_tick_nums); %  
        set(gca,'XTicklabel', x_tick_names); % 
        axis([0 x_tick_max 0.0 1.0]);
        ylabel(cancer_name);
        box on;
        if i == 1
            title('Methylation');
        end
        hold on;
        subplot(nrow, ncol, ncol*(i-1) + 2);
        %2. plot log(mRNA-value) for each stage
        for k = 1 : size_stage(1)
             stage_name = char(stage_names(k));
             rna_data_for_stage = load([intermediate_dir, char(data_dir_names(2)), is_merge_stage, cancer_name , '/', cancer_name, '_', stage_name, char(file_ends(2))]);
             rna_data_of_gene = rna_data_for_stage(target_gene_idxs(j , 1) + 1, :)';
             size_y_data = size(rna_data_of_gene);
             x_data = ones(size_y_data(1), 1) .* k + (0.4 * rand(size_y_data(1), 1) - 0.2);
             plot(x_data(:, 1), log10(rna_data_of_gene(:, 1)), char(stage_markers(k)));
             hold on;
        end
        set(gca, 'YTick', [-2, -1, 0, 1, 2]);
        set(gca, 'XTick', x_tick_nums); %  
        set(gca,'XTicklabel', x_tick_names); % 
        axis([0 x_tick_max -2.5 2.5]);
        box on;
        if i == 1
            title('log_{10}(mRNA)');
        end
        
        hold on;
        subplot(nrow, ncol, ncol*(i-1) + 3);
        %3. plot mutation rate for each stage
        for k = 1 : size_stage(1)
             if k == 1
                 continue;
             end
             stage_name = char(stage_names(k));
             mutation_rate_for_stage = load([intermediate_dir, char(data_dir_names(3)), is_merge_stage, cancer_name , '/', cancer_name, '_', stage_name, char(file_ends(3))]);
             mutation_rate_of_gene = mutation_rate_for_stage(target_gene_idxs(j , 1) + 1, :)';
             size_y_data = size(mutation_rate_of_gene);
             x_data = ones(size_y_data(1), 1) .* k + (0.4 * rand(size_y_data(1), 1) - 0.2);
             plot(x_data(:, 1), mutation_rate_of_gene(:, 1), char(stage_markers(k)),'markersize',20);
             hold on;
        end
        set(gca, 'YTick', [-2, -1, 0, 1, 2]);
        set(gca, 'XTick', x_tick_nums); %  
        set(gca,'XTicklabel',x_tick_names); % 
        axis([0 x_tick_max 0.0 1.0]);
        box on;
        if i == 1
            title('SNV Mutation Rate');
        end
    end
    set(gcf,'Name',gene_name);
    exportfig(fig, fig_path , 'FontMode', 'fixed', 'color', 'cmyk','width',6 ,'height', 9 , 'FontSize', 12,'Resolution',300,'LineWidth',0.5); % , 
    close all;
end
end