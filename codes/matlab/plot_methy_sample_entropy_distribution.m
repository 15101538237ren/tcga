function plot_methy_sample_entropy_distribution()
clc;
clear all;
close all;

sample_entropy_dir= '../../data/intermediate_file/methy_sample_entropy/merged_stage/';
cancer_names = {'BRCA'; 'COAD'; 'KIRC'; 'KIRP'; 'LIHC'; 'LUAD'; 'LUSC'; 'THCA'};
stage_names = {'normal';'i';'ii';'iii';'iv'};
ncol = 4;
nrow = 2;
figure_dir = '../../figures/sample_entropy/';
if ~exist(figure_dir)
    mkdir(figure_dir);
end

fig=figure(1);
clf();
dx=0.04;
min_entropy = 1.5;
max_entropy = 3.0;
for i= 1 : length(cancer_names)
    cancer_name = char(cancer_names(i));
    subplot(nrow, ncol, i);
    hold on;
    for j = 1 : length(stage_names)
        sample_entropy= load(strcat(sample_entropy_dir, cancer_name, '_', char(stage_names(j)),'_sample_entropy.dat'));
        [u, x]=hist(sample_entropy(:,2), min_entropy: dx: max_entropy);
        plot(x, u/(dx*sum(u)));
    end
    if i == 2
       legend({'n','i','ii','iii','iv'},'Box','off');
    end
    title(cancer_name);
    box on;
    axis([min_entropy max_entropy 0.0 8.0]);
    if i > 4
        xlabel('entropy');
    end
    if i==1 || i == 5
        ylabel('cases (%)');
    end
end
fig_save_path = strcat(figure_dir, 'sample_entropy_distribution.eps');
exportfig(fig, fig_save_path , 'FontMode', 'fixed', 'color', 'cmyk','width',8, 'height', 4, 'FontSize', 12,'Resolution',300,'LineWidth',0.5);
close all;
