function plotmethy()
G=load('../../global_files/gene_label.dat');
p_score_path = '../../data/intermediate_file/methy_pvalue/merged_stage/COAD/COAD_p_score.dat';
n_score_path = '../../data/intermediate_file/methy_pvalue/merged_stage/COAD/COAD_n_score.dat';
SP=load(p_score_path);
SN=load(n_score_path);

xticks = 0:0.2:1;

fig=figure(1);

subplot(2,3,1);
J=find(G(:,2)==1 & SP(:,4)>0 & G(:,4)==1);
plot(sort(SP(J,4),'descend'),'k-','linewidth',2);

K=SP(J,:);
[X,J0]=sort(K(:,4),'descend');
J0(1:5)
K(J0(1:5),1)
K(J0(1:5),4)
set(gca,'ytick',xticks);
xlim([0 160]);
ylim([0 1]);
xlabel('Gene');
ylabel('M^+ Score');
title('(a) Onco genes');

subplot(2,3,2);
J=find(G(:,2)==1 & SP(:,4)>0 & G(:,4)==2);
plot(sort(SP(J,4),'descend'),'k-','linewidth',2);
set(gca,'ytick',xticks);
xlim([0 300]);
ylim([0 1]);
xlabel('Gene');
ylabel('M^+ Score');
title('(b) TSG genes');

subplot(2,3,3);
J=find(G(:,2)==1 & SP(:,4)>0 & G(:,4)==0);
plot(sort(SP(J,4),'descend'),'k-','linewidth',2);
set(gca,'ytick',xticks);
xlim([0 4500]);
ylim([0 1]);
xlabel('Gene');
ylabel('M^+ Score');
title('(c) Other genes');

subplot(2,3,4);
J=find(G(:,2)==1 & SN(:,4)>0 & G(:,4)==1);
plot(sort(SN(J,4),'descend'),'k-','linewidth',2);
set(gca,'ytick',xticks);
xlim([0 70]);
ylim([0 1]);
xlabel('Gene');
ylabel('M^- Score');
title('(d) Onco genes');

subplot(2,3,5);
J=find(G(:,2)==1 & SN(:,4)>0 & G(:,4)==2);
plot(sort(SN(J,4),'descend'),'k-','linewidth',2);
set(gca,'ytick',xticks);
xlim([0 100]);
ylim([0 1]);
xlabel('Gene');
ylabel('M^- Score');
title('(e) TSG genes');

subplot(2,3,6);
J=find(G(:,2)==1 & SP(:,4)>0 & G(:,4)==0);
plot(sort(SN(J,4),'descend'),'k-','linewidth',2);
set(gca,'ytick',xticks);
xlim([0 600]);
ylim([0 1]);
xlabel('Gene');
ylabel('M^- Score');
title('(f) Other genes');

figdir = '../../figures/plotmethy/';

if ~exist(figdir)
    mkdir(figdir);
end

exportfig(fig,strcat(figdir, 'MScore.eps'),'color','cmyk','fontmode','fixed','fontsize',10);
close all;
end