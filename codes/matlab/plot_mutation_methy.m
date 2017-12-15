function mutations()
global cancer_name

mutation_dat_dir = '../../data/intermediate_file/snv_intermidiate/merged_stage/';
cancer_name = 'COAD';
R=load(strcat(mutation_dat_dir, cancer_name,'/', cancer_name,'_i_mutation_rate.txt'));

fig=figure(1);
clf();

subplot(1,4,1:2);
hold on;
Z=zeros(20,3);
pvalue_dir = '../../data/intermediate_file/methy_pvalue/merged_stage/';
P=load(strcat(pvalue_dir, cancer_name,'/', cancer_name,'_p_score.dat'));
N=load(strcat(pvalue_dir, cancer_name,'/', cancer_name,'_n_score.dat'));
[X,J]=sort(R(:,2),'descend');
J0=J(1:20);
for i=1:19
    Z(i,1)=X(i);
    Z(i,2)=P(J0(i),4);
    Z(i,3)=N(J0(i),4);
end
b=bar(Z,'stacked');
b(1).FaceColor='blue';
b(2).FaceColor='yellow';
b(3).FaceColor='red';
xlim([0 20]);
ylim([0 1.5]);
legend('mutation rate','M^+','M^-','Location','northwest','Orientation','horizontal');
box on;
text(0.5,1.05,'TTN');
text(2.5,0.9,'APC');
text(9.5,1.01,'ZFHX4');
text(14.5,1.01,'RYR2');
xlabel('genes');
ylabel('values');
title('(a) Mutation rate and M-Scores \newline of the 20 highest mutant genes'); %0, 1.7, '                    '
hold on;
methy_dat_dir = '../../data/intermediate_file/methy_intermidiate/merged_stage/';
An=load(strcat(methy_dat_dir, cancer_name,'/', cancer_name,'_normal_methy_dat.dat'));
Ai=load(strcat(methy_dat_dir, cancer_name,'/', cancer_name,'_i_methy_dat.dat'));

TTN=17529;
APC=831;
RYR2=14195;
ZFHX4=18516;
J1=[17529,831,18516,14195];
names={'TTN';'APC';'ZFHX4';'RYR2'};
r0=0.1*rand(1,size(An,2)-1);
r1=0.1*rand(1,size(Ai,2)-1);
for i=1:4
    if(i<=2)
        subplot(2,4,2+i);
    else
        subplot(2,4,4+i);
    end
    hold on;
    plot(0.1+r0,An(J1(i)+1,2:end),'k.','markersize',10);
    plot(0.3+r1,Ai(J1(i)+1,2:end),'k.','markersize',10);
    box on;
    if(i==1)
        text(0.03,0.1,char(names(i)));
    else
        text(0.03,0.8,char(names(i)));
    end
    xlim([0 0.5]);
    ylim([0 1]);
    set(gca,'xtick',[0.15,0.35]);
    set(gca,'xticklabel',{''});
    if(i>2)
        set(gca,'xticklabel',{'normal','i'});
    end
    if(i==1 || i==3)
        ylabel('beta-value');
    end
    if(i==1)
        text(0.15, 1.25,'(b) DNA methylation of selected genes','FontWeight','bold');%, 
    end
end

exportfig(fig,'mutation_methy.eps','color','cmyk','width',9,'height',3,'FontSize', 0.9);
close all;
end