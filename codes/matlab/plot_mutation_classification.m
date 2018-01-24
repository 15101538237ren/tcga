global base_path
clc;
clear all;
close all;
base_path = './';

sample_base_dir = 'I:/intermediate_file/common_patients_data/';
pval_base_dir = 'I:/intermediate_file/methy_pvalue/merged_stage/';
figure_base_dir = '../../figures/common_patients_data/';
global_file_dir = '../../global_files/';

gene_idx_file_path = strcat(global_file_dir, 'gene_idx.txt');
gene_label_path = strcat(global_file_dir, 'gene_label.dat');

mutation_classification_file_path = strcat(global_file_dir, 'mutation_classification.txt');
mutation_classification_list = cell(1);
num = 0;
fid = fopen(mutation_classification_file_path, 'r');
while ~feof(fid)                                      % 判断是否为文件末尾               
    tline=fgetl(fid);                                 % 从文件读行
    tmp_line = regexp(tline, '\t', 'split');
    num = num + 1;
    mutation_classification_list{num} = char(tmp_line(1,2));
end
fclose(fid);

cancer_name='COAD';
cancer_stage='i';
gene_id = 831;  % APC gene 831
mutation_type_num = 14;
sample_id_list = zeros(1,100);
sample_num = 0;

cancer_sample_dir = strcat(sample_base_dir, cancer_name, '/', cancer_stage, '/');
figdir = strcat(figure_base_dir, cancer_name,'/');
sample_id_file_name = 'common_patients_idxs.txt';
sample_id_file_path = strcat(cancer_sample_dir, sample_id_file_name);

if ~exist(figdir)
    mkdir(figdir);
end

gene_label_info = load(gene_label_path);
gene_body_absolute_start = gene_label_info(gene_id, 9);
gene_body_absolute_end = gene_label_info(gene_id, 10);
tss_area_pos = -2000;
tss_len = 2000;
global gene_body_len;
gene_body_len = gene_body_absolute_end - gene_body_absolute_start;
subfigure_num = 2;

fid = fopen(sample_id_file_path, 'r');
while ~feof(fid)                                      % 判断是否为文件末尾               
    tline=fgetl(fid);                                 % 从文件读行
    tmp_line = regexp(tline, '\t', 'split');
    sample_num = sample_num + 1;
    sample_id_list(1, sample_num) = str2num(char(tmp_line(1,1)));
end
fclose(fid);

for sampleId = 1:sample_num
%sampleId = 2;
    methy_file_path = strcat(cancer_sample_dir, num2str(sampleId), '_methy.tsv');
    ins_file_path = strcat(cancer_sample_dir, num2str(sampleId), '_INS.tsv');
    del_file_path = strcat(cancer_sample_dir, num2str(sampleId), '_DEL.tsv');
    snp_file_path = strcat(cancer_sample_dir, num2str(sampleId), '_SNP.tsv');
    
    fig = figure(sampleId);
    clf();
    title(num2str(sampleId));
    cur_fig_num = subfigure_num;
    space_height = 0.03;
    left = 0.03;
    height = (1 - (subfigure_num + 1) * space_height) / subfigure_num;
    bottom = cur_fig_num * space_height + (cur_fig_num - 1) * height;
    width = 0.94;
    
    sub_left = -0.5;
    sub_mid = 0;
    sub_right = 1;
    sub_down = 0;
    sub_down_y = -0.4;
    sub_up = 1;
    
    % first figure
    axe_region = [left bottom width height];
    plot_region = [sub_left sub_mid sub_right sub_down sub_down_y sub_up];
    idx = subfigure_num - cur_fig_num + 1;
    plotCoord(idx, axe_region, plot_region);
    
    % second figure
    cur_fig_num = cur_fig_num  - 1;
    bottom = cur_fig_num * space_height + (cur_fig_num - 1) * height;
    axe_region = [left bottom width height];
    idx = subfigure_num - cur_fig_num + 1;
    plotCoord(idx, axe_region, plot_region);
    plot_methy(methy_file_path, gene_id, tss_len, sub_down);
    
%     for typeId = 1:mutation_type_num
%         cur_fig_num = cur_fig_num  - 1;
%         idx = subfigure_num - cur_fig_num + 1;
%         bottom = cur_fig_num * space_height + (cur_fig_num - 1) * height;
%         axe_region = [left bottom width height];
%         plotCoord(idx, axe_region, plot_region);
%         mutation_type = mutation_classification_list{typeId};
%         plot_mutation(ins_file_path, gene_id, typeId, 0);
%         plot_mutation(del_file_path, gene_id, typeId, 1);
%         plot_mutation(snp_file_path, gene_id, typeId, 2);
%     end
    %set(gcf,'visible','off');
%    print(fig,strcat(figdir, num2str(sampleId), '.pdf'),'-dpdf','-opengl');
    close all;
    %break;
end

function [] = plotCoord(fig_num, axe_region, plot_region)
    left = axe_region(1);
    bottom = axe_region(2);
    width = axe_region(3);
    height = axe_region(4);
    
    if(size(plot_region,1) > 0)
        sub_left = plot_region(1);
        sub_mid = plot_region(2);
        sub_right = plot_region(3);
        sub_down = plot_region(4);
        sub_down_y = plot_region(5);
        sub_up = plot_region(6);
    end
    axes('position', [left bottom  width height]);
    hold on;
    plot([sub_left sub_right],[sub_down sub_down],'color','k');
    xs = [sub_left sub_mid sub_right];
    global gene_body_len;
    xtext = {'-2000'; '0(TSS)'; num2str(gene_body_len)};
    fs = 8;
    for xi = 1: length(xs)
      text(xs(xi)-0.04,sub_down_y,char(xtext(xi)),'FontSize',fs);
      %if(fig_num == 1)
      %    plot([xs(xi) xs(xi)],[sub_down 0.3],'color','k');
      %end
    end
    if(fig_num == 1)
        plot([sub_mid sub_mid], [sub_down, 0.5], 'color', 'k');
        annotation('arrow',[left+width / 3 left + width / 3 + 0.15],[bottom + height/2 bottom + height/2]); 
    end
    xlim([sub_left sub_right]);
    ylim([sub_down sub_up]);
    box off;
    axis off;
end

function [] = plot_methy(methy_file_path, gene_id, tss_len, sub_down)
    global gene_body_len;
    fid = fopen(methy_file_path, 'r');
    cnt = 0;
    while(~feof(fid))
        tline = fgetl(fid);
        tmp_line = regexp(tline, '\t', 'split');
        cnt = cnt + 1;
        if(cnt == gene_id)
            gene_methy_list = regexp(tmp_line{1,2}, ';', 'split');
            pos_num = size(gene_methy_list, 2);
            for j = 1:pos_num
                pos_info_list = regexp(gene_methy_list{1,j}, ',', 'split');
                relative_pos = str2num(pos_info_list{1, 1});
                beta_value = str2double(pos_info_list{1, 3});
                if(relative_pos < 0)
                    coor_x = relative_pos / tss_len;
                    plot([coor_x coor_x], [sub_down beta_value], 'color', 'k');
                else
                    coor_x = relative_pos / gene_body_len;
                    plot([coor_x coor_x], [sub_down beta_value], 'color', 'k');
                end
            end
            break
        end
    end
    fclose(fid);
end


function [] = plot_mutation(mutation_file_path, gene_id, mutation_type_num, file_type, sub_down)
    global gene_body_len;
    fid = fopen(mutation_file_path, 'r');
    cnt = 0;
    while(~feof(fid))
        tline = fgetl(fid);
        if(isfloat(tline) && tline == -1)
            break
        end
        tmp_line = regexp(tline, '\t', 'split');
        idx = str2num(tmp_line{1,1});
        if(idx == gene_id)
            %fprintf('');
            gene_mutation_info = regexp(tmp_line{1,3}, ',', 'split');
            if(file_type == 0)  %INS
                start_pos = str2num(gene_mutation_info{1,1});
                end_pos = str2num(gene_mutation_info{1,2});
                tmp_type = str2num(gene_mutation_info{1,5});
                
                if(tmp_type ~= mutation_type_num)
                    break
                end
                fprintf('YES\n');
                
                x = (start_pos :end_pos) / gene_body_len;
                y = ones(1, end_pos - start_pos + 1) * 0.2;
                plot(x, y, '-vr');      
            elseif(file_type == 1)  %DEL
                start_pos = str2num(gene_mutation_info{1,1});
                end_pos = str2num(gene_mutation_info{1,2});
                tmp_type = str2num(gene_mutation_info{1,5});
                if(tmp_type ~= mutation_type_num)
                    break
                end
                x = (start_pos :end_pos) / gene_body_len;
                y = ones(1, end_pos - start_pos + 1) * 0.2;
                plot(x, y, '-^r', 'MarkerSize',10);
            else  %SNP
                pos = str2num(gene_mutation_info{1,1});
                tmp_type = str2num(gene_mutation_info{1,3});
                if(tmp_type ~= mutation_type_num)
                    break
                end
                x = pos / gene_body_len;
                y = 0.2;
                plot(x, y, 'v', 'MarkerSize',10);
            end
            break
        end
    end
    fclose(fid);
end
