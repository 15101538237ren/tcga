library(GenomicFeatures)
data_path = "~/PycharmProjects/tcga_raw_data/GRCh38/"
out_path = "~/PycharmProjects/tcga/global_files/"
txdb <- makeTxDbFromGFF(paste(data_path, "Homo_sapiens.GRCh38.90.gtf",sep = ""), format="gtf")
genes <- genes(txdb)
gene_df <- as.data.frame(genes)
ensembl_ids<-gene_df[,6]
library(biomaRt)
mart<-useMart("ensembl")
mart<- useDataset("hsapiens_gene_ensembl", mart)
attributs_df<- as.data.frame(listAttributes(mart)) #list all available attributes

genes_table <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "hgnc_symbol", "description"), values= ensembl_ids, mart= mart)
merged_gene_df<-merge(gene_df, genes_table, by.x="gene_id", by.y="ensembl_gene_id")
df_out <- merged_gene_df[ ,c(7,2,3,4,6)]
write.table(df_out, file= paste(out_path, "human_gene_bodys.tsv",sep = ""), col.names=F, row.names=F, sep="\t")