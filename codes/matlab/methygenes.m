function methygenes()
type = 3;
%An=load('../COAD/COAD_normal_methy_dat.dat');
%Ai=load('../COAD/COAD_i_methy_dat.dat');

An=load('../../data/intermediate_file/methy_intermidiate/merged_stage/COAD/COAD_normal_methy_dat.dat');
Ai=load('../../data/intermediate_file/methy_intermidiate/merged_stage/COAD/COAD_i_methy_dat.dat');

dx=0.02;

figdir = '../../figures/methygenes/';

if ~exist(figdir)
    mkdir(figdir);
end

if(type == 1)
    %Oncognes
    typeStr = 'Onco';
    geneid=[5729,574,8853,3825];  
    names={'FGF5'; 'ALK'; 'LMO1'; 'CTNND2'};
elseif(type == 2)
    %TSG
    typeStr = 'TSG';
    geneid=[16538,1697,13165,14615];  
    names={'THBD'; 'BTG4'; 'PTPRT'; 'SFRP2'};
elseif(type == 3)
    %Other
    typeStr = 'Other';
    geneid=[3464,4536,15457,15571];
    names={'COL25A1', 'DOK6', 'SORCS1', 'SPATA32'};
else
    %Both
    typeStr = 'Both';
    geneid=[18266,5821,7602,1125];
    names={'WT1','FLT3','IKZF1','ASCL1'};
end
lx=[0.3 0.70 0.3 0.7];
ly=[0.7 0.75 0.25 0.25];
fig=figure(1);
clf();
for i=1:4
    subplot(2,2,i);
    hold on;
    X=An(geneid(i)+1,2:end);
    Y=Ai(geneid(i)+1,2:end);
    [u,x]=hist(X,0:dx:1);
    plot(x,100*u/sum(u),'b-');
    [u,x]=hist(Y,0:dx:1);
    plot(x,100*u/sum(u),'r-');
    xlabel('beta value');
    ylabel('cases (%)');
    xlim([0 1]);
    ylim([0 100]);
    box on;
    title(char(names(i)));
    
    axes('position',[lx(i), ly(i), 0.15, 0.15]);
    hold on;
    n=size(X,2);
    z=0.1*rand(1,n);
    plot(0.1+z,X,'b.');
    n=size(Y,2);
    z=0.1*rand(1,n);
    plot(0.3+z,Y,'r.');
    xlim([0 0.5]);
    ylim([0 1]);
    box on;
    set(gca,'fontsize',6);
    set(gca,'xtick',[0.15,0.35]);
    set(gca,'xticklabel',{'normal','i'});
end
exportfig(fig,strcat(figdir,typeStr,'.eps'),'color','cmyk','fontmode','fixed','fontsize',10);
close all;
end