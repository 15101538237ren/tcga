function plotpvalue()
clc;
clear all;
close all;
global An Ai Aii Aiii Aiv PP PN Sp Sn E cancer_name
figure_path = '../../figures/plot_pvalue/';
if ~exist(figure_path) 
    mkdir(figure_path)
end 
cancer_name = 'COAD';
stage_of_pscore_observed = 'i';
out_significant_pscore_gidx_dir = ['../../data/intermediate_file/methy_corr/merged_stage/',cancer_name,'/'];
if ~exist(out_significant_pscore_gidx_dir)
    mkdir(out_significant_pscore_gidx_dir)
end
methy_dat_path = ['../../data/intermediate_file/methy_intermidiate/merged_stage/',cancer_name,'/'];
An=load(strcat(methy_dat_path,cancer_name,'_normal_methy_dat.dat'));
Ai=load(strcat(methy_dat_path,cancer_name,'_i_methy_dat.dat'));
Aii=load(strcat(methy_dat_path,cancer_name,'_ii_methy_dat.dat'));
Aiii=load(strcat(methy_dat_path,cancer_name,'_iii_methy_dat.dat'));
Aiv=load(strcat(methy_dat_path,cancer_name,'_iv_methy_dat.dat'));

pvalue_dat_path = ['../../data/intermediate_file/methy_pvalue/merged_stage/',cancer_name,'/'];
PP=load(strcat(pvalue_dat_path,cancer_name,'_pp_value.dat'));
PN=load(strcat(pvalue_dat_path,cancer_name,'_pn_value.dat'));
Sp=load(strcat(pvalue_dat_path,cancer_name,'_p_score.dat'));
Sn=load(strcat(pvalue_dat_path,cancer_name,'_n_score.dat'));

entropy_path = ['../../data/intermediate_file/methy_entropy/merged_stage/',cancer_name,'_entropy.dat'];
E=load(entropy_path);

m=50; % Totoal sample number;

fig=figure(1);
clf();
gene_name = 'WT1';
WT1=18267;

subplot(2,2,1);
hold on;
dx=0.03;
[u,x]=hist(An(WT1,2:end),0:dx:1);
plot(x,100*u/sum(u));

[u,x]=hist(Ai(WT1,2:end),0:dx:1);
plot(x,100*u/sum(u));

[u,x]=hist(Aii(WT1,2:end),0:dx:1);
plot(x,100*u/sum(u));

[u,x]=hist(Aiii(WT1,2:end),0:dx:1);
plot(x,100*u/sum(u));

[u,x]=hist(Aiv(WT1,2:end),0:dx:1);
plot(x,100*u/sum(u));

legend('normal','i','ii','iii','iv');
box on;
xlabel('beta-value');
ylabel('cases (%)');

xlim([0 0.7]);
ylim([0 80]);

% n=size(An,2)-1;
% r0=0.1*rand(1,n);
% plot(0.2+r0,An(WT1,2:end),'k.');
% n=size(Ai,2)-1;
% r0=0.1*rand(1,n);
% plot(0.6+r0,Ai(WT1,2:end),'k.');
% box on;
% xlim([0 1]);
% xlabel('stage');
% ylabel('Beta-value');
title('(a) Distribuiton of Beta-values of WT1');

subplot(2,2,2);
hold on;
E0=E(2:end,2:end);
plot([0 20000], [1 1],'color','r','linestyle','--');
plot([0 20000], [0 0],'color','r');
plot(sort(E0(:,2)-E0(:,1),'descend'),'k.');
%plot([0 20000], [1 1],'color','r','linestyle,'--');
%plot([0 4],[1 5],'color','r','linestyle','--');
%xlim([0 3]);
ylim([-1 5]);
box on;
xlabel('genes');
ylabel('\Delta Entropy');
title('(b) Entropy increase');

subplot(2,2,3);
hold on;
G=load('../../global_files/gene_label.dat');
J=find(E0(:,2)-E0(:,1)>1 & E0(:,2)>0 & E0(:,1)>0 & G(:,2)>0);
Sp0=Sp(J,:);
Sn0=Sn(J,:);
[B,J]=sort(Sp0(:,4),'descend');
plot(Sp0(J,4),'r.');
plot(Sn0(J,4),'b.');
line([0 3500],[0.3 0.3],'color','k','linestyle','--');
J0=find(Sp0(:,4)>0.3);
sp0_gene_idxs = Sp0(J0,1);

n0=size(J0,1);
idxs1 = 1 : n0;
out_mat1 = [idxs1', sp0_gene_idxs(:, 1)]';
out_sp0_gidx_fp = [out_significant_pscore_gidx_dir, cancer_name,'_',stage_of_pscore_observed, '_pp_gidx.txt'];
fid1 = fopen(out_sp0_gidx_fp,'w');
fprintf(fid1,'%d\t%d\n', out_mat1);
fclose(fid1);

text(200, 0.9,strcat(num2str(n0),' genes'),'color','r','fontweight','bold');
%text(200, 0.83,'with M^+-score > 0.3','color','r','fontweight','bold');

J0=find(Sn0(:,4)>0.3);
n0=size(J0,1);
sn0_gene_idxs = Sn0(J0,1);
idxs2 = 1 : n0;
out_mat2 = [idxs2', sn0_gene_idxs(:, 1)]';
out_sn0_gidx_fp = [out_significant_pscore_gidx_dir, cancer_name,'_',stage_of_pscore_observed, '_pn_gidx.txt'];
fid2 = fopen(out_sn0_gidx_fp,'w');
fprintf(fid2,'%d\t%d\n', out_mat2);
fclose(fid2);

text(1000, 0.9,strcat(num2str(n0),' genes'),'color','b','fontweight','bold');
%text(800, 0.83,'with M^--score > 0.3','color','b','fontweight','bold');
xlim([0 1800]);
ylim([0 1]);
xlabel('genes');
ylabel('M-score');
box on;
title('(c) M-Score');

subplot(2,2,4);
J=find(E0(:,2)-E0(:,1)>1 & E0(:,2)>0 & E0(:,1)>0 & G(:,2)>0);
PP0=PP(J+1,2:end);
PN0=PN(J+1,2:end);
Z=zeros(m,2);
for i=1:m
    J0=find(PP0(:,i)<-10);
    n0=size(J0,1);
    Z(i,1)=n0;
    J0=find(PN0(:,i)<-10);
    n0=size(J0,1);
    Z(i,2)=n0;
end
b=bar(Z,'stacked');
b(2).FaceColor='yellow';
legend('log(P^+) <-10','log(P^-)<-10');
xlabel('Samples');
ylabel('genes');
xlim([0 51]);
ylim([0 2200]);
title('(d) Number of abnormal genes in each sample');

exportfig(fig,strcat(figure_path,gene_name,'.eps'),'color','cmyk','fontmode','fixed','fontsize',10);
close all;

end