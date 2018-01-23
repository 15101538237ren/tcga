clc;
clear all;
close all;
figure_path = '../../figures/entropy/';
if ~exist(figure_path) 
    mkdir(figure_path)
end 

entropy_base_dir = '../../data/intermediate_file/methy_entropy/merged_stage/';

gene_label_path = '../../global_files/gene_label.dat';
gene_label = [zeros(1,11); load(gene_label_path)];

gene_category_col_idx = 4; % the column index of the gene category in the gene_label.dat
gene_categories = [[0, 0];[1, 3];[2, 3];[0, 0]]; %Genome doesn't have category, therefore, the 1st element in this array have no use ; Onco category: 1 or 3; TSG: 2 or 3; Other: 0; %;
n_gene_categories = 4;
figure_filename_start = {'Genome';'Onco';'TSG';'Other'};

cancername={'BRCA'; 'COAD'; 'KIRC'; 'KIRP'; 'LIHC'; 'LUAD'; 'LUSC'; 'THCA'};

for gc = 1 : n_gene_categories
    fig=figure(gc);
    clf();
    dx=0.02;
    % filter the rows of gene_category
    if gc ==2 || gc==3
        rows_filtered = [find(gene_label(:, gene_category_col_idx)==gene_categories(gc , 1));find(gene_label(:, gene_category_col_idx)==gene_categories(gc , 2))]; %
    elseif gc == 4
        rows_filtered = find(gene_label(:, gene_category_col_idx)==gene_categories(gc , 1));
    else
        rows_filtered = find(gene_label(:, 1)>=0);
    end
    
    for i=1:8
        A=load(strcat(entropy_base_dir, char(cancername(i)),'_entropy.dat'));
        X=zeros(4,2);
        fl_i = floor(i/5);
        subplot(4,4,4 * fl_i + i);
        hold on;
        for k=1:4
            J=find(A(:,k+1)>0);
            intersection_rows = intersect(J, rows_filtered);
            [u,x]=hist(A(intersection_rows,k+1),0:dx:4);
            X(k,1)=k;
            X(k,2)=mean(A(intersection_rows,k+1));
            plot(x,100*u/sum(u));
        end
        if(i==1)
            h0=legend('normal','i','ii','iii');
            set(h0,'box','off');
        end
        box on;
        xlabel('entropy');
        ylim([0 4]);
        set(gca,'ytick',[1,2,3,4]);
        if (mod(i, 4) == 1)
            ylabel('cases (%)');
        end
        set(gca,'xtick',[1,2,3,4]);
        title(cancername(i));
        
        subplot(4,4, 4 * fl_i + i + 4);
        plot(X(:,1),X(:,2),'k--o','markerfacecolor','k','markersize',5);
        xlim([0 5]);
        ylim([1 2.5]);
        set(gca,'xtick',[1,2,3,4]);
        set(gca,'ytick',[1,1.5,2,2.5]);
        set(gca,'xticklabel',{'normal','i','ii','iii'});
        if (mod(i, 4) == 1)
            ylabel('entropy');
        end
    end
    exportfig(fig,strcat(figure_path, char(figure_filename_start(gc)),'_entropy.eps'),'color','cmyk');
end
close all;