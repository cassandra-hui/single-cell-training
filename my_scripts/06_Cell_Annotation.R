library(celldex)
library(SingleR)

###################################
# You can decide on manual or automatic cell annotations

# SingleR uses references you select
# Made by Seurat with Human and Mouse data (https://azimuth.hubmapconsortium.org/)

# With mannual annotation you can use, Seurat AddModuleScore()
# But if you want to incorporate negative genes you need to use UCell (on Bioconductor, https://github.com/carmonalab/UCell)


###################################

# Set default 
seu <- Seurat::SetIdent(seu, value = seu$RNA_snn_res.0.3)

# Use original coutn data not intergrated data
DefaultAssay(seu) <- "RNA"

Seurat::FeaturePlot(seu, "C3")


## Manually define cell types

tcell_genes <- c("Il7r", "Ltb", "Trac", "Cd3d")
monocyte_genes <- c("Cd14", "Cst3", "Cd68", "Ctss")
Seurat::FeaturePlot(seu, tcell_genes, ncol=2)
Seurat::FeaturePlot(seu, monocyte_genes, ncol=2)

Seurat::VlnPlot(seu,
                features = tcell_genes,
                ncol = 2)

Seurat::VlnPlot(seu,
                features = monocyte_genes,
                ncol = 2)


# Calculate the score
seu <- Seurat::AddModuleScore(seu,
                              features = list(monocyte_genes),
                              name = "monocyte_genes")

# a column was added to seu@meta.data.
seu@meta.data

Seurat::FeaturePlot(seu, "monocyte_genes1")
Seurat::VlnPlot(seu,
                "monocyte_genes1")


### Cell Cylce Scoring
# Extract built in genes
s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes


seu <- Seurat::CellCycleScoring(seu,
                                s.features = s.genes,
                                g2m.features = g2m.genes)

seu@meta.data
## This still worked, but gave soem errors I thought because the gene names are not all uppercase 

Seurat::DimPlot(seu, group.by = "Phase")


## Cell annotation using SingleR

# Need to select a dataset from scRNAseq pacakge in Bioconductor
# Let's try one from 'celldex'
#ref <- celldex::NovershternHematopoieticData()
library(celldex)
ref <- MouseRNAseqData()
class(ref)
table(ref$label.main)


seu_SingleR <- SingleR::SingleR(test = Seurat::GetAssayData(seu),
                                ref = ref,
                                labels = ref$label.main)

head(seu_SingleR)

SingleR::plotScoreHeatmap(seu_SingleR)
SingleR::plotDeltaDistribution(seu_SingleR)


## Remove cell types with less than 10 cells so not to clogg our plots

singleR_labels <- seu_SingleR$labels
t <- table(singleR_labels)
other <- names(t)[t < 10]
singleR_labels[singleR_labels %in% other] <- "none"

## Add to seu object
seu$SingleR_annot <- singleR_labels


dittoSeq::dittoDimPlot(seu, "SingleR_annot", size = 0.7)
dittoSeq::dittoBarPlot(seu, var = "SingleR_annot", group.by = "orig.ident")

dittoSeq::dittoBarPlot(seu, 
                       var = "SingleR_annot", 
                       group.by = "RNA_snn_res.0.3")

saveRDS(seu, "seu_day2-4.rds")

