## Enrichment

library(clusterProfiler)
library(enrichplot)


BiocManager::install("org.Mm.eg.db", update = FALSE)
library(org.Mm.eg.db)
AnnotationDbi::keytypes(org.Mm.eg.db)


tum_down  <- subset(limma_de,
                    limma_de$logFC < -1 
                    & limma_de$adj.P.Val <  0.05)
tum_down_genes <- rownames(tum_down)


## no limma genes so using seuart 

tum_down  <- subset(tum_vs_norm,
                    tum_vs_norm$avg_log2FC < -1 
                    & tum_vs_norm$p_val_adj <  0.05)
tum_down_genes <- rownames(tum_down)


?enrichGO


tum_vs_norm_go <- clusterProfiler::enrichGO(tum_down_genes,
                                            "org.Mm.eg.db",
                                            keyType = "SYMBOL",
                                            ont = "BP",
                                            minGSSize = 50)

View(tum_vs_norm_go@result)

# Remove redundant gene sets
enr_go <- clusterProfiler::simplify(tum_vs_norm_go)
View(enr_go@result)

gc()
enrichplot::emapplot(enrichplot::pairwise_termsim(enr_go),
                     showCategory = 30,
                     cex.params = list(category_label = 0.5))


# Error in `geom_edge_link()`:
#   ! Problem while converting geom to grob.
# ℹ Error occurred in the 1st layer.
# Caused by error:
#   ! vector memory limit of 32.0 Gb reached, see mem.maxVSize()
# Run `rlang::last_trace()` to see where the error occurred.


gmt <- msigdbr::msigdbr(species = "Mus musculus")


tum_vs_norm_enrich <- clusterProfiler::enricher(gene = tum_down_genes,
                                                universe = rownames(proB),
                                                pAdjustMethod = "BH",
                                                pvalueCutoff  = 0.05,
                                                qvalueCutoff  = 0.05,
                                                TERM2GENE = gmt[,c("gs_name", "gene_symbol")])

View(tum_vs_norm_enrich@result[tum_vs_norm_enrich@result$p.adjust < 0.05,])


gc()

