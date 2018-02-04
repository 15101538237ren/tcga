function plot_mutation_classification()
global sample_base_dir;
global sample_num;
global cancer_name;
global gene_id;
global gene_name;
global tss_len;
global figdir;
global gene_body_len;
global cancer_stage;
global out_ave_methy_file;
% global out_ins_num_file;
% global out_del_num_file;
% global out_snp_num_file;
global normal_mean_methy_file;

clc;
%clear;
close all;

fig = figure();
clf();


py.get_gene_data.get_all_sample_methy(sample_base_dir, cancer_name, cancer_stage, gene_id, sample_num);
py.get_gene_data.get_all_sample_mutation(sample_base_dir, cancer_name, cancer_stage, gene_id, sample_num);

normal_methy_list = load(normal_mean_methy_file);
J = find(normal_methy_list(:,1) < 0);
normal_methy_list(J,1) = normal_methy_list(J,1) / tss_len;
J = find(normal_methy_list(:,1) >= 0);
normal_methy_list(J,1) = normal_methy_list(J,1) / gene_body_len;

methy_list = load(out_ave_methy_file);
% ins_list = load(out_ins_num_file);
% del_list = load(out_del_num_file);
% snp_list = load(out_snp_num_file);

J = find(methy_list(:,1) < 0);
methy_list(J,1) = methy_list(J,1) / tss_len;
J = find(methy_list(:,1) >= 0);
methy_list(J,1) = methy_list(J,1) / gene_body_len;

subplot(2,1,1);
bar(normal_methy_list(:,1), normal_methy_list(:,2), 'k', 'LineWidth',1.5);
title_str = sprintf('%s %s normal average beta-value', cancer_name, gene_name);
title([title_str]);
xlim([-0.2 1]);
ylim([0 1]);
hold on;


subplot(2,1,2);
bar(methy_list(:,1), methy_list(:,2), 'k', 'LineWidth',1.5);
title_str = sprintf('%s %s i-th average beta-value', cancer_name, gene_name);
title([title_str]);
xlim([-0.2 1]);
ylim([0 1]);
% %two subfigure
% %average methy
% subplot(4,1,1);
% bar(methy_list(:,1), methy_list(:,2), 'k', 'LineWidth',1.5);
% title(['COAD APC i-th average beta-value']);
% xlim([-0.2 1]);
% ylim([0 1]);
% hold on;
% 
% %mutation list, first column: pos, second column: count
% J = find(ins_list(:,1) < 0);
% ins_list(J,1) = ins_list(J,1) / tss_len;
% J = find(ins_list(:,1) >= 0);
% ins_list(J,1) = ins_list(J,1) / gene_body_len;
% J = find(del_list(:,1) < 0);
% del_list(J,1) = del_list(J,1) / tss_len;
% J = find(del_list(:,1) >= 0);
% del_list(J,1) = del_list(J,1) / gene_body_len;
% J = find(snp_list(:,1) < 0);
% snp_list(J,1) = snp_list(J,1) / tss_len;
% J = find(snp_list(:,1) >= 0);
% snp_list(J,1) = snp_list(J,1) / gene_body_len;
% 
% %y1 = ones(1, size(ins_list, 1)) * 0.05;
% %y2 = ones(1, size(del_list, 1)) * 0.2;
% %y3 = ones(1, size(snp_list, 1)) * 0.35;
% 
% pos1 = 0.08;  %start_height_ins/snp
% pos2 = 0.03;  %start_height_del
% ins = 0.03;  %marker_height_interval
% 
% %mutation num
% subplot(4,1,2);
% height = max(ins_list(:,2));
% for i = 1:height
%     J = find(ins_list(:,2) >= i);
%     plot(ins_list(J,1), i*(pos1 + (i - 1)* ins), 'v', 'MarkerSize',4, 'Markerfacecolor','r','markeredgecolor','r');
%     hold on;
% end
% title(['COAD APC INS mutation']);
% xlim([-0.2 1]);
% ylim([0 1]);
% 
% subplot(4,1,3);
% height = max(del_list(:,2));
% for i = 1:height
%     J = find(del_list(:,2) >= i);
%     plot(del_list(J,1), i*(pos2 + (i - 1)* ins), '^', 'MarkerSize',4, 'Markerfacecolor','b','markeredgecolor','b');
%     hold on;
% end
% title(['COAD APC DEL mutation']);
% xlim([-0.2 1]);
% ylim([0 1]);
% 
% subplot(4,1,4);
% height = max(snp_list(:,2));
% for i = 1:height
%     J = find(snp_list(:,2) >= i);
%     plot(snp_list(J,1), i*(pos1 + (i - 1)* ins), 'v', 'MarkerSize',4, 'Markerfacecolor','m','markeredgecolor','m');
%     hold on;
% end
% title(['COAD APC SNP mutation']);
% xlim([-0.2 1]);
% ylim([0 1]);
% 
print(fig,strcat(figdir, 'COAD_',gene_name,'_avg_methy_and_mutation.pdf'),'-dpdf','-opengl');
close all;
delete('*.txt');
end