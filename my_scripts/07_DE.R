## Differential Gene Expression

## To try all the DE methods(https://www.nature.com/articles/nmeth.4612)

# https://github.com/csoneson/conquer_comparison in scripts


library(Seurat)
library(edgeR) # BiocManager::install("edgeR")
library(limma)
library(dplyr)
library(scuttle)

seu <- readRDS("my_objects/seu_day2-4.rds")



?Seurat::FindAllMarkers
## test.use can use other types of tests

de_genes <- Seurat::FindAllMarkers(seu,  min.pct = 0.25,
                                   only.pos = TRUE)


top_specific_markers <- de_genes %>%
  group_by(cluster) %>%
  top_n(3, avg_log2FC)


dittoSeq::dittoDotPlot(seu,
                       vars = unique(top_specific_markers$gene), 
                       group.by = "RNA_snn_res.0.3")

tcell_genes <- c("Kit", "Ano1")

de_genes[de_genes$gene %in% tcell_genes,] |> knitr::kable()


seu <- Seurat::SetIdent(seu, value = "SingleR_annot")
seu@meta.data


deg_cd8_cd4 <- Seurat::FindMarkers(seu,
                                   ident.1 = "Fibroblasts",
                                   ident.2 = "Endothelial cells",
                                   group.by = seu$SingleR_annot,
                                   test.use = "wilcox")
deg_cd8_cd4 <- subset(deg_cd8_cd4, deg_cd8_cd4$p_val_adj<0.05)

View(deg_cd8_cd4)

## avg_logFC: log fold-chage of the average expression between the two groups. 
# Positive values indicate that the gene is more highly expressed in the first group

deg_cd8_cd4[c("Ptprb", "Cdh5", "Fabp4"),]


Seurat::VlnPlot(seu, 
                features = c("Ptprb", "Cdh5", "Fabp4"),
                idents = c("Fibroblasts", "Endothelial cells"))


# The Wilcoxon test implemented in FindMarkers does not allow you to test for complex design (eg factorial experiments) or to include batch as a covariate. 

###########################################
### Using Limma

## Testing Sal's data here
seu <- readRDS("../../Projects/Sal_Baker/scRNAseq_2024/obj.rds")

Seurat::DimPlot(seu, group.by = "orig.ident")

table(seu@meta.data$group) # If we had different types I will use group
head(seu@meta.data)
tail(seu@meta.data)
# Add a new column 'type' to the metadata based on orig.ident
seu$type <- ifelse(seu$tissue %in% c("Control.Region"), "control", "treatment")
head(seu@meta.data)

Seurat::DimPlot(seu, group.by = "type")

## Testing between two "types" control vs treatment

proB <- seu

#taking the data 
Seurat::DefaultAssay(proB) <- "RNA"
Seurat::Idents(proB) <- proB$orig.ident

## add the patient id also for paired DGE
proB$patient.id<-proB$orig.ident
#proB$patient.id<-sapply(strsplit(proB$patient.id, "-"), '[', 2)

head(proB@meta.data)

## Here we do perform pseudo-bulk:
##first a mandatory column of sample needs to be added to the meta data that is the grouping factor, should be the samples
proB$sample <- factor(proB$orig.ident)

# aggergate the cells per sampple
bulk <- Seurat::AggregateExpression(proB, group.by = "sample",
                                    return.seurat = TRUE,
                                    assay = "RNA")

# create a metadata data frame based on the aggregated cells
meta_data <- unique(proB@meta.data[, c("orig.ident",
                                       "sample", 
                                       "type",
                                       "group"
                                       )])
rownames(meta_data) <- meta_data$orig.ident

print(colnames(bulk))
# "sample-1" "sample-2" "sample-3" "sample-4"
## Need to change names 
# Replace hyphens with underscores in column names
#colnames(bulk) <- gsub("-", "_", colnames(bulk))
#print(colnames(bulk))

bulk@meta.data <- meta_data[colnames(bulk), ]
print(colnames(bulk))
##have a look at the counts
counts <- Seurat::GetAssayData(bulk, layer = "counts") |> as.matrix()

head(counts)

#have a look at the colData of our new object summed, can you see type and 
#patient.id are there
head(bulk@meta.data)


# orig.ident   sample      type
# sample_1   sample_1 sample_1   control
# sample_2   sample_2 sample_2   control
# sample_3   sample_3 sample_3 treatment
# sample_4   sample_4 sample_4 treatment


# orig.ident
# Obstruction.recovery-post.2.weeks.-.control.region.-.rep1         Obstruction.recovery-post.2.weeks.-.control.region.-.rep1
# Obstruction.recovery-post.2.weeks.-.control.region.-.rep2         Obstruction.recovery-post.2.weeks.-.control.region.-.rep2
# Obstruction.recovery-post.2.weeks.-.obstruction.region.-.rep1 Obstruction.recovery-post.2.weeks.-.obstruction.region.-.rep1
# Obstruction.recovery-post.2.weeks.-.obstruction.region.-.rep2 Obstruction.recovery-post.2.weeks.-.obstruction.region.-.rep2
# sample
# Obstruction.recovery-post.2.weeks.-.control.region.-.rep1         Obstruction.recovery-post.2.weeks.-.control.region.-.rep1
# Obstruction.recovery-post.2.weeks.-.control.region.-.rep2         Obstruction.recovery-post.2.weeks.-.control.region.-.rep2
# Obstruction.recovery-post.2.weeks.-.obstruction.region.-.rep1 Obstruction.recovery-post.2.weeks.-.obstruction.region.-.rep1
# Obstruction.recovery-post.2.weeks.-.obstruction.region.-.rep2 Obstruction.recovery-post.2.weeks.-.obstruction.region.-.rep2
# type
# Obstruction.recovery-post.2.weeks.-.control.region.-.rep1       control
# Obstruction.recovery-post.2.weeks.-.control.region.-.rep2       control
# Obstruction.recovery-post.2.weeks.-.obstruction.region.-.rep1 treatment
# Obstruction.recovery-post.2.weeks.-.obstruction.region.-.rep2 treatment


#As in the standard limma analysis generate a DGE object

y <- edgeR::DGEList(counts, samples = bulk@meta.data)

##filter lowly expressed (recommanded for limma)
keep <- edgeR::filterByExpr(y, group = bulk$type)
y <- y[keep,]

##see how many genes were kept 
summary(keep)

# Mode   FALSE    TRUE 
# logical   11556   14334 


# Mode   FALSE    TRUE 
# logical    4312   12658

## Create the design matrix and include the technology as a covariate:
design <- model.matrix(~0 + y$samples$type + y$samples$group) # This is for a factorial design. 
#design <- model.matrix(~0 + y$samples$type)

# Have a look
design

# change column/rownames names to more simple group names: 
colnames(design) <- make.names(c("tissue.Control", "tissue.Obst"
                                 ,
                                 #"patient1", why is this one missing??
                                 "Partial", "Recovery"
                                 ))
rownames(design) <- rownames(y$samples)
design

contrast.mat <- limma::makeContrasts(tissue.Obst - tissue.Control,
                                     levels = design)

dge <- edgeR::calcNormFactors(y)  

#Do limma
vm <- limma::voom(dge, design = design, plot = TRUE)


fit <- limma::lmFit(vm, design = design)
fit.contrasts <- limma::contrasts.fit(fit, contrast.mat)
fit.contrasts <- limma::eBayes(fit.contrasts)

# Show the top differentially expressed genes:
limma::topTable(fit.contrasts, number = 10, sort.by = "P")

limma_de <- limma::topTable(fit.contrasts, number = Inf, sort.by = "P")
length(which(limma_de$adj.P.Val<0.05))

Seurat::VlnPlot(proB, "Elavl4", split.by = "type")
Seurat::VlnPlot(proB, "Chrna7", split.by = "type")


tum_vs_norm <- Seurat::FindMarkers(proB, 
                                   ident.1 = "control", 
                                   ident.2 = "treatment", 
                                   group.by = "type")
tum_vs_norm <- subset(tum_vs_norm, tum_vs_norm$p_val_adj<0.05)

dim(tum_vs_norm) 
# 2932    5


merge_limma_FindMarkers <- merge(tum_vs_norm, limma_de, by="row.names",
                                 all.x=T)

par(mar=c(4,4,4,4))
plot(merge_limma_FindMarkers$avg_log2FC,
     merge_limma_FindMarkers$logFC,
     xlab="log2FC Wilcoxon", ylab="log2FC limma",
     pch=15, cex=0.5)
abline(a=0, b=1, col="red")

