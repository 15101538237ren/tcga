library(gplots)
library(RColorBrewer)
setwd('~/PycharmProjects/tcga')
cancers = c('COAD')
base_path = 'data/intermediate_file/common_sample_cnt/merged_stage';
gene_idx_filepath = 'global_files/gene_idx.txt';
gene_idxs = read.csv(gene_idx_filepath,header = F,na.strings = "-1",sep = '\t',stringsAsFactors = F)
gene_types = c('significant_genes_all') #'significant_genes','significant_no_genes','significant_methy_genes','significant_mut_genes'
my_palette = colorRampPalette(c("green", "yellow", "red"))(n = 300)
figure_dir = 'figures/common_samples'
cexrow = 0.8
ks = 0.8
pointsz = 12
res = 300
wid = 50
key_title = 'Common Samples'
for (i in 1 : length(cancers))
{
  
  for (j in 1 : length(gene_types))
  {
    if (j ==2 || j == 5)
    {
      wid = 50
    }
    else if (j ==3 || j == 4)
    {
      wid = 40
    }
    cancer_fig_fp = paste(figure_dir, cancers[i], gene_types[j], sep = '/')
    if (!file.exists(cancer_fig_fp))
    {
      dir.create(file.path(cancer_fig_fp), recursive=T)
    }

    dir_path = paste(base_path,cancers[i],gene_types[j], sep = '/')
    sig_gene_fp = paste(dir_path, 'sig_gidxs.dat', sep = '/')
    sig_gidxs = read.csv(sig_gene_fp, header = F, sep='\t',stringsAsFactors = F)
    gene_names = as.character(gene_idxs[ sig_gidxs[ , 2], 2])
    
    mut_fp = paste(dir_path, 'common_mutation_samples_cnt.dat', sep = '/')
    mut_data = read.csv(mut_fp, header = T, row.names = 1, sep='\t') 
    rownames(mut_data)  = gene_names
    colnames(mut_data) = gene_names
    mut_fig_fp = paste(cancer_fig_fp,'mut_heatmap.png', sep = '/')
    png(mut_fig_fp,width = wid*res,height = wid*res, res = res,pointsize = pointsz)
    gplots:::heatmap.2(as.matrix(mut_data),Rowv=F,Colv=F,col = my_palette,tracecol=NA,density.info="none",keysize = ks,cexRow=cexrow,cexCol=cexrow)
    dev.off()
    
    methy_fp = paste(dir_path, 'common_methy_sig_variation_samples_cnt.dat', sep = '/')
    methy_data = read.csv(methy_fp, header = T, row.names = 1, sep='\t')
    rownames(methy_data)  = gene_names
    colnames(methy_data) = gene_names
    
    methy_fig_fp = paste(cancer_fig_fp,'methy_heatmap.png', sep = '/')
    png(methy_fig_fp,width = wid*res,height = wid*res, res = res,pointsize = pointsz)
    gplots:::heatmap.2(as.matrix(methy_data),Rowv=F,Colv=F,col = my_palette,tracecol=NA,density.info="none",keysize = ks,cexRow=cexrow,cexCol=cexrow)
    dev.off()
  }
}
 