function plot_mutation_methy_comparasion()
G=load('../../global_files/gene_label.dat');
p_score_path = '../../data/intermediate_file/methy_pvalue/merged_stage/COAD/COAD_p_score.dat';
SP=load(p_score_path);
mutation_rate_path = '../../data/intermediate_file/snv_intermidiate/merged_stage/COAD/COAD_i_mutation_rate.txt';
MR=load(mutation_rate_path);

xticks = 0:0.2:1;
ZSP = SP(:,4);
ZSP(find(ZSP<0)) = 0.0;

fig=figure(1);
[sSP,sidxs] = sort(ZSP,'descend');
subplot(2,2,1);
plot(MR(sidxs,2),'r.');
hold on;
plot(sSP,'b.'); 
legend('Mutation Rate','M^+ Score');
set(gca,'ytick',xticks);
ylim([0 1]);
xlabel('Gene');

figdir = '../../figures/plotmut_methy_comparasion/';

if ~exist(figdir)
    mkdir(figdir);
end

exportfig(fig,strcat(figdir, 'COAD.eps'),'color','cmyk','fontmode','fixed','fontsize',10);
close all;
end