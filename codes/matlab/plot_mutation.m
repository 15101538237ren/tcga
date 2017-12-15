function plotmutation()
global M0 cancer_name

mutation_dat_dir = '../../data/intermediate_file/snv_intermidiate/merged_stage/';
cancer_name = 'COAD';
M0=load(strcat(mutation_dat_dir, cancer_name,'/', cancer_name,'_i_mutation_data.dat'));
R=load(strcat(mutation_dat_dir, cancer_name,'/', cancer_name,'_i_mutation_rate.txt'));

M=M0(2:end,2:end);
m=size(M,2);
X0=sum(M);
[X,J]=sort(X0);

fig=figure(1);
clf();

stages={'i';'ii';'iii';'iv'};
col={'b','g','r','k'};
subplot(3,1,1);
hold on;
for k=1:1
    S0=load(strcat(mutation_dat_dir, cancer_name,'/', cancer_name,'_',char(stages(k)),'_mutation_data.dat'));
    S=sum(S0(2:end,2:end));
    plot(sort(S),'-o','markerfacecolor',char(col(k)),'markersize',5);
end
box on;
xlabel('Sample');
ylabel('Number of mutant genes');
xlim([0 65]);
title('(a) Numer of mutant genes in each COAD stage i sample.');

subplot(3,1,2);
hold on;
for i=1:m
    j=J(i);  % The j'th sample
    J0=find(M(:,j)>0); % The genes in the j'th sample with mutation
    n=size(J0,1);      % number of genes with mutation.
    x=i-0.8*(rand(n,1)-0.5);
    R0=R(J0,2);
    plot(sort(x),sort(R0),'.');
end
box on;
xlabel('Sample');
ylabel('mutation rates');
xlim([0 65]);
ylim([0 0.7]);
title('(b) Mutation rates of each mutant gene in each COAD stage i sample');

mutation_dat_dir = '../../data/intermediate_file/snv_intermidiate/merged_stage/';
M0=load(strcat(mutation_dat_dir, cancer_name,'/', cancer_name,'_i_mutation_data.dat'));
M=M0(2:end,2:end);
m=size(M,2);
X0=sum(M);
[X,J]=sort(R(:,2),'descend');
J0=J(1:20)
M1=sign(M(J0,:));

subplot(3,1,3);
pcolor(M1);
colormap winter;
xlabel('Samples');
ylabel('Genes');
title('(c) Mutations of genes in each sample');

exportfig(fig,'mutation.eps','color','cmyk');
close all;
end