function plot_methy_mean_and_std()
% This function plot the scatter of the methylation means and its standard
% variations of different stages for all genes
%plot_mean_std();
plot_mean_dist();
end

function plot_mean_std()
fig=figure(1);
clf();

% 1. change work directory
mean_std_data_path = '../../data/intermediate_file/methy_mean_std_data/merged_stage/';
figure_path = '../../figures/methy_mean_and_std_figure/';
cancer_names = {'BRCA'; 'COAD'; 'KIRC'; 'KIRP'; 'LIHC'; 'LUAD'; 'LUSC'; 'THCA'};
size_cancer_names = size(cancer_names);

stage_names = {'normal';'i'}; %;'ii';'iii';'iv'
size_stage = size(stage_names);

stage_markers = {'g.';'r.'};%;'b.';'k.';'c.'

% counter for figures
fig_counter = 0;

for i=1: size_cancer_names(1)
    cancer_name = char(cancer_names(i));
    fig_counter = fig_counter + 1;
    mean_data_raw = load(fullfile(mean_std_data_path, cancer_name,[cancer_name, '_mean.dat']));
    std_data_raw = load(fullfile(mean_std_data_path, cancer_name,[cancer_name, '_std.dat']));
    mean_data = mean_data_raw(2:end,2:end);
    std_data = std_data_raw(2:end,2:end);
    %iterate stage
    for j= 1: size_stage(1)
        subplot(4,4,2*(i-1)+j);
%        plot(mean_data(:,j),(std_data(:,j).^2)./(mean_data(:,j).^2),char(stage_markers(j))); 
        plot(mean_data(:,j),std_data(:,j),char(stage_markers(j)),'markersize',3); 
 
        xlabel('\mu');
        ylabel('\sigma');
        axis([0 1.0 0.0 0.3]);
        title(strcat(cancer_name,'-',stage_names(j)));
    end
end
fig_save_path = [figure_path, 'cv2.eps'];
exportfig(fig, fig_save_path , 'FontMode', 'fixed', 'color', 'cmyk', 'FontSize', 12,'Resolution',300,'LineWidth',0.5);
end

function plot_mean_dist()
fig=figure(1);
clf();

% 1. change work directory
mean_std_data_path = 'methy_mean_std_data/merged_stage/';
figure_path = 'methy_mean_and_std_figure/';
cancer_names = {'BRCA'; 'COAD'; 'KIRC'; 'KIRP'; 'LIHC'; 'LUAD'; 'LUSC'; 'THCA'};
size_cancer_names = size(cancer_names);

stage_names = {'normal';'i'}; %;'ii';'iii';'iv'
size_stage = size(stage_names);

stage_markers = {'g-';'r-'};%;'b.';'k.';'c.'

% counter for figures
fig_counter = 0;
dx=0.02;

for i=1: size_cancer_names(1)
    cancer_name = char(cancer_names(i));
    fig_counter = fig_counter + 1;
    mean_data_raw = load(fullfile(mean_std_data_path, cancer_name,[cancer_name, '_mean.dat']));
    std_data_raw = load(fullfile(mean_std_data_path, cancer_name,[cancer_name, '_std.dat']));
    mean_data = mean_data_raw(2:end,2:end);
    %iterate stage
    subplot(4,2,i);
    hold on;
    
    for j= 1: size_stage(1)
        [u,x]=hist(mean_data(:,j),0:dx:1);
        plot(x,u/(dx*sum(u)),char(stage_markers(j)));
    end
    legend('normal','stage i');
    box on;
    xlabel('\beta');
    ylabel('# Samples');
    xlim([0 1]);
    title(cancer_name);

end
fig_save_path = [figure_path, 'mean.eps'];
exportfig(fig, fig_save_path , 'FontMode', 'fixed', 'color', 'cmyk', 'FontSize', 12,'Resolution',300,'LineWidth',0.5);
end

