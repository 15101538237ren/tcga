function plotclassification()
global L1
genes={'BRCA';'COAD';'LIHC';'KIRC';'KIRP';'LUAD';'LUSC';'THCA'};
fpre = '../../data/intermediate_file/';
mp_score_threshold = 0.8;
mutation_rate_threshold = 0.1;
base_path = strcat(fpre, 'gene_classification_mp_',num2str(mp_score_threshold),'_mut_',num2str(mutation_rate_threshold),'/');
%base_path = '../../data/intermediate_file/gene_classification/';
J0=load(strcat(base_path,'genes_sig.ind'));
L=load('../../global_files/gene_label.dat');
fig=figure(1);
clf();

Onco=1;
Tsg=2;
Both=3;
Zo=find(L(:,4)==Onco);
Zt=find(L(:,4)==Tsg);
Zb=find(L(:,4)==Both);
Z0=find(L(:,4)==0);

Vo=find(L(:,5)==Onco);
Vt=find(L(:,5)==Tsg);
Vb=find(L(:,5)==Both);
V0=find(L(:,5)==0);

H1=find(ismember(J0,Zo)==1);
H2=find(ismember(J0,Zt)==1);
H3=find(ismember(J0,Zb)==1);
H4=find(ismember(J0,Z0)==1);

hy=0.7;
c=['y','b','g','r'];
col=['k','r','g','w'];

H=[H1;H2;H3];
n=size(H,1);
L1=L(J0(H),:);


axes('position',[0.2 0.92 hy 0.02]);
hold on;
box off;
axis off;

fill([0 0 0.05 0.05],[0 1 1 0],'b','linestyle','none');
text(0.07,0.3,'Onco');

fill([0.2 0.2 0.25 0.25],[0 1 1 0],'g','linestyle','none');
text(0.27,0.3,'TSG');

fill([0.45 0.45 0.5 0.5],[0 1 1 0],'r','linestyle','none');
text(0.52,0.3,'Onco & TSG');

fill([0.8 0.8 0.85 0.85],[0 1 1 0],'y','linestyle','none');
text(0.87,0.3,'Other');

xlim([0 1]);
ylim([0 1]);

axes('position',[0.2 0.8 hy 0.1]);
hold on;
axis off;
for k=1:n
    fill([k-1, k-1,k,k],[0, 0.5, 0.5, 0],c(1+L1(k,4)),'linestyle','none');
    fill([k-1, k-1,k,k],[0.5, 1,1,0.5],c(1+L1(k,5)),'linestyle','none');
end

xlim([0 n]);
ylim([0 1]);

text(-floor(0.116*n),0.2,'Zhao list','fontsize',8);
text(-floor(0.116*n),0.6,'Vogel list','fontsize',8);

axes('position',[0.2 0.48 hy 0.3]);
hold on;
axis on;

for i=1:8
    A1=load(strcat(base_path, char(genes(i)),'/',char(genes(i)),'_genome_class_',num2str(1),'.dat')); 
    A2=load(strcat(base_path, char(genes(i)),'/',char(genes(i)),'_genome_class_',num2str(2),'.dat')); 
    A3=load(strcat(base_path, char(genes(i)),'/',char(genes(i)),'_genome_class_',num2str(3),'.dat')); 
    for k=1:n
        if(ismember(L1(k,1),A1)==1)
            fill([k-1, k-1,k,k],[i-1, i,i,i-1],col(1),'linestyle','none');
        elseif(ismember(L1(k,1),A2)==1)
            fill([k-1, k-1,k,k],[i-1, i,i,i-1],col(2),'linestyle','none');
        elseif(ismember(L1(k,1),A3)==1)
            fill([k-1, k-1,k,k],[i-1, i,i,i-1],col(3),'linestyle','none');
        end
    end
    text(-floor(0.097*n),i-0.8,char(genes(i)),'fontsize',8);
end
box on;
xlim([0 n]);
ylim([0 8]);
xlabel('genes');
set(gca,'yticklabel',{''});

axes('position',[0.2 0.4 hy 0.02]);
hold on;
box off;
axis off;

fill([0 0 0.01 0.01],[0 1 1 0],'k','linestyle','none');
text(0.03,0.3,'DNA Methylation');

fill([0.4 0.4 0.41 0.41],[0 1 1 0],'r','linestyle','none');
text(0.43,0.3,'Mutation');

fill([0.7 0.7 0.71 0.71],[0 1 1 0],'g','linestyle','none');
text(0.73,0.3,'Both');

xlim([0 1]);
ylim([0 1]);

exportfig(fig,strcat(base_path,'significantgenes','_mp_',num2str(mp_score_threshold),'_mut_',num2str(mutation_rate_threshold),'.eps'),'color','cmyk','fontmode','scaled','fontsize',1.0);

clf();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H=H4;
n=size(H,1);
L1=L(J0(H),:);


axes('position',[0.2 0.92 hy 0.02]);
hold on;
box off;
axis off;

fill([0 0 0.05 0.05],[0 1 1 0],'b','linestyle','none');
text(0.07,0.3,'Onco');

fill([0.2 0.2 0.25 0.25],[0 1 1 0],'g','linestyle','none');
text(0.27,0.3,'TSG');

fill([0.45 0.45 0.5 0.5],[0 1 1 0],'r','linestyle','none');
text(0.52,0.3,'Onco & TSG');

fill([0.8 0.8 0.85 0.85],[0 1 1 0],'y','linestyle','none');
text(0.87,0.3,'Other');

xlim([0 1]);
ylim([0 1]);

axes('position',[0.2 0.8 hy 0.1]);
hold on;
axis off;
for k=1:n
    fill([k-1, k-1,k,k],[0, 0.5, 0.5, 0],c(1+L1(k,4)),'linestyle','none');
    fill([k-1, k-1,k,k],[0.5, 1,1,0.5],c(1+L1(k,5)),'linestyle','none');
end

xlim([0 n]);
ylim([0 1]);

text(-floor(0.116*n),0.2,'Zhao list','fontsize',8);
text(-floor(0.116*n),0.6,'Vogel list','fontsize',8);

axes('position',[0.2 0.48 hy 0.3]);
hold on;
axis on;

for i=1:8
    A1=load(strcat(base_path, char(genes(i)),'/',char(genes(i)),'_genome_class_',num2str(1),'.dat')); 
    A2=load(strcat(base_path, char(genes(i)),'/',char(genes(i)),'_genome_class_',num2str(2),'.dat')); 
    A3=load(strcat(base_path, char(genes(i)),'/',char(genes(i)),'_genome_class_',num2str(3),'.dat')); 
    for k=1:n
        if(ismember(L1(k,1),A1)==1)
            fill([k-1, k-1,k,k],[i-1, i,i,i-1],col(1),'linestyle','none');
        elseif(ismember(L1(k,1),A2)==1)
            fill([k-1, k-1,k,k],[i-1, i,i,i-1],col(2),'linestyle','none');
        elseif(ismember(L1(k,1),A3)==1)
            fill([k-1, k-1,k,k],[i-1, i,i,i-1],col(3),'linestyle','none');
        end
    end
    text(-floor(0.097*n),i-0.8,char(genes(i)),'fontsize',8);
end
box on;
xlim([0 n]);
ylim([0 8]);
xlabel('genes');
set(gca,'yticklabel',{''});

axes('position',[0.2 0.4 hy 0.02]);
hold on;
box off;
axis off;

fill([0 0 0.01 0.01],[0 1 1 0],'k','linestyle','none');
text(0.03,0.3,'DNA Methylation');

fill([0.4 0.4 0.41 0.41],[0 1 1 0],'r','linestyle','none');
text(0.43,0.3,'Mutation');

fill([0.7 0.7 0.71 0.71],[0 1 1 0],'g','linestyle','none');
text(0.73,0.3,'Both');

xlim([0 1]);
ylim([0 1]);

exportfig(fig,strcat(base_path,'significantgenes-no','_mp_',num2str(mp_score_threshold),'_mut_',num2str(mutation_rate_threshold),'.eps'),'color','cmyk','fontmode','scaled','fontsize',1.0);
close all;
end