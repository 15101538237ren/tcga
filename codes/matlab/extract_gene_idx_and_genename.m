function gene_names = extract_gene_idx_and_genename(file_path)
[idxs, gene_names] = textread(file_path,'%d\t%s');
end