function plotpvalue()
global An Ai Aii Aiii Aiv PP PN Sp Sn E cancer_name
cancer_name = 'COAD';

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
APC=832;

subplot(2,2,1);
hold on;
dx=0.03;
[u,x]=hist(An(APC,2:end),0:dx:1);
plot(x,100*u/sum(u));

[u,x]=hist(Ai(APC,2:end),0:dx:1);
plot(x,100*u/sum(u));

[u,x]=hist(Aii(APC,2:end),0:dx:1);
plot(x,100*u/sum(u));

[u,x]=hist(Aiii(APC,2:end),0:dx:1);
plot(x,100*u/sum(u));

[u,x]=hist(Aiv(APC,2:end),0:dx:1);
plot(x,100*u/sum(u));

legend('normal','i','ii','iii','iv');
box on;
xlabel('beta-value');
ylabel('cases (%)');


xlim([0 0.7]);
ylim([0 80]);
text(-0.03, 95,'(a)','fontsize',12,'fontweight','bold');

% n=size(An,2)-1;
% r0=0.1*rand(1,n);
% plot(0.2+r0,An(APC,2:end),'k.');
% n=size(Ai,2)-1;
% r0=0.1*rand(1,n);
% plot(0.6+r0,Ai(APC,2:end),'k.');
% box on;
% xlim([0 1]);
% xlabel('stage');
% ylabel('Beta-value');
title('Distribuiton of Beta-values of APC');

subplot(2,2,2);
hold on;
E0=E(2:end,2:end);
plot(E0(:,1),E0(:,2),'.');
plot([0 4], [0 4],'color','r');
plot([0 4],[1 5],'color','r','linestyle','--');
xlim([0 3]);
ylim([0 3.5]);
text(0, 4.15,'(b)','fontsize',12,'fontweight','bold');
box on;
xlabel('normal');
ylabel('i');
title('Entropy');

subplot(2,2,3);
hold on;
J=find(E0(:,2)-E0(:,1)>1 & E0(:,2)>0 & E0(:,1)>0);
Sp0=Sp(J,:);
Sn0=Sn(J,:);
[B,J]=sort(Sp0(:,4),'descend');
plot(Sp0(J,4),'r.');
plot(Sn0(J,4),'b.');
line([0 3500],[0.3 0.3],'color','k','linestyle','--');
J0=find(Sp0(:,4)>0.3);
n0=size(J0,1);
text(400, 0.9,strcat(num2str(n0),' genes'),'color','r','fontweight','bold');
text(400, 0.83,'with M^+-score > 0.3','color','r','fontweight','bold');

J0=find(Sn0(:,4)>0.3);
n0=size(J0,1);
text(2000, 0.9,strcat(num2str(n0),' genes'),'color','b','fontweight','bold');
text(2000, 0.83,'with M^--score > 0.3','color','b','fontweight','bold');

ylim([0 1]);
text(-200, 1.2,'(c)','fontsize',12,'fontweight','bold');
xlabel('genes');
ylabel('M-score');
box on;
title('M-Score');

subplot(2,2,4);
J=find(E0(:,2)-E0(:,1)>1 & E0(:,2)>0 & E0(:,1)>0);
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
bar(Z,'stacked');
legend('log(P^+) <-10','log(P^-)<-10');
xlabel('Samples');
ylabel('genes');
xlim([0 51]);
ylim([0 2200]);
title('Number of abnormal genes in each sample');
text(-3, 2650,'(d)','fontsize',12,'fontweight','bold');

exportfig(fig,'APC.eps','color','cmyk');%,'FontSize', 12,'Resolution',300,'LineWidth',0.5
close all;
end