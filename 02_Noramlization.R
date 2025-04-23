library(Seurat)

# Chcek data before and after Normalization

Seurat::GetAssayData(seu)[15:30,1:30]  

# LOC103027786 . . . . . . . 2 . . . . . . . . . . . . . . . . . . . 1 1 .
# LOC125804717 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# LOC125781472 . . . . . . . . . . . . . . . . . . . . . . . . 1 . . . . .
# LOC111189986 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# LOC111189988 . . . . 1 . . . . . . . . . . . . . . . . . . . . . . . . .
# LOC125804618 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# LOC103047048 . . . . . . . . . . . 1 . . . . . . . . . . . . . . . . . .


seu <- Seurat::NormalizeData(seu,
                             normalization.method = "LogNormalize",
                             scale.factor = 10000)


Seurat::GetAssayData(seu)[15:30,1:30]  

# LOC103027786 . . . . .        . . 2.406118 . . . .        . . . . . . . . . . . . .       .
# LOC125804717 . . . . .        . . .        . . . .        . . . . . . . . . . . . .       .
# LOC125781472 . . . . .        . . .        . . . .        . . . . . . . . . . . . 1.25941 .
# LOC111189986 . . . . .        . . .        . . . .        . . . . . . . . . . . . .       .
# LOC111189988 . . . . 2.624965 . . .        . . . .        . . . . . . . . . . . . .       .
# LOC125804618 . . . . .        . . .        . . . .        . . . . . . . . . . . . .       .
# LOC103047048 . . . . .        . . .        . . . 1.555387 . . . . . . . . . . . . .       .


seu <- Seurat::FindVariableFeatures(seu,
                                    selection.method = "vst",
                                    nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(Seurat::VariableFeatures(seu), 10)
top10

vf_plot <- Seurat::VariableFeaturePlot(seu)
Seurat::LabelPoints(plot = vf_plot,
                    points = top10, repel = TRUE)

# When using repel, set xnudge and ynudge to 0 for optimal results
# Warning message:
#   In scale_x_log10() :
#   log-10 transformation introduced infinite values.

## Increase plotting size

par(mar = c(7, 8, 7, 4)) # doesn't work, save bigger instead
vf_plot <- Seurat::VariableFeaturePlot(seu)
plot<- Seurat::LabelPoints(plot = vf_plot,
                    points = top10, repel = TRUE)
ggsave('plot.pdf', plot, width = 10, height = 8)


## Scaling

# Next, we apply scaling, a linear transformation that is a standard pre-processing step prior to dimensional reduction techniques like PCA. The ScaleData() function
# 
# shifts the expression of each gene, so that the mean expression across cells is 0
# scales the expression of each gene, so that the variance across cells is 1
# This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate. The results of this are stored in seu$RNA@scale.data


seu <- Seurat::ScaleData(seu)

## OR
seu <- Seurat::SCTransform(seu)

## This errored for Misty's data, maybe it is too big. 

# And it will add an extra assay to the object. names(seu@assays) returns:
#   
#   [1] "RNA" "SCT"
# Meaning that a whole new assay was added (including the sparse matrices with counts, normalized data and scaled data).
names(seu@assays)

## Warning
# Running SCTransform will change @active.assay into SCT(in stead of RNA; check it with DefaultAssay(seu)). 
# This assay is used as a default for following function calls. To change the active assay to RNA run:

DefaultAssay(seu) <- "RNA"

# Save
saveRDS(seu, "seu_day1-3.rds")


# And clear environment

rm(list = ls())
gc()
.rs.restartR()

