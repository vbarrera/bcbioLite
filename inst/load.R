library(bcbioLite)
library(SummarizedExperiment)


se = bcbreader("final_project_path")
se = bcbrun(se)

library(AnnotationHub)
## Load the annotation resource.
ah <- AnnotationHub()
ahDb <- query(ah, pattern =c("Homo sapiens", "EnsDb") )
ahEdb <- ahDb[[length(ahDb)]]
rows = genes(ahEdb) %>% 
    as.data.frame() %>% 
    group_by(gene_id, gene_biotype, chrom=seqnames, start, end, strand) %>% 
    summarise(gene_name = paste(unique(gene_name), collapse = ","),
              entrezid = paste(unique(entrezid), collapse = ","),
              description = paste(unique(description), collapse = ",")) %>% 
    mutate(gene_name = ifelse(gene_name=="", gene_id, gene_name)) %>% 
    as.data.frame()
row.names(rows) = rows[["gene_id"]]
rowData(se) = rows[names(se),]

saveRDS(se, file = "data/se.rds")
