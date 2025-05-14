#Seurat

library(Seurat)
library(ggplot2)

data_path <- "~/Documents/Projects/Misty_Riddle/scRNAseq_2024"
dataset_loc <- "~/Documents/Projects/Misty_Riddle/scRNAseq_2024"
ids <- c("M", "P", "S", "T")

# Workshop codes but needs barcodes.tsv.gz
########## 

#sample_info <- read.csv("course_data/sample_info_course.csv")

#datadirs <- file.path(data_path, ids,
#                      "outs")
#names(datadirs) <- gsub("_", "-", sample_info$SampleName)
#datadirs <- datadirs[1:3]
#datadirs
#sparse_matrix <- Seurat::Read10X(data.dir = datadirs)
#######

### Metrics
#####
d10x.metrics <- lapply(ids, function(i){
  metrics <- read.csv(file.path(dataset_loc,paste0(i,"/outs"),"metrics_summary.csv"), colClasses = "character")
})


experiment.metrics <- do.call("rbind", d10x.metrics)
rownames(experiment.metrics) <- ids
sequencing.metrics <- data.frame(t(experiment.metrics[,c(1:19)]))
####

# Data
#####
d10x.data <- lapply(ids, function(i){
  d10x <- Read10X_h5(file.path(dataset_loc, i, "/outs","filtered_feature_bc_matrix.h5"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="_")
  d10x
})
names(d10x.data) <- ids

str(d10x.data)

sparse_matrix <- do.call("cbind", d10x.data)
#sparse_matrix[c("slc25a35", "LOC125803521", "rangrf"), 1:50]


seu <- Seurat::CreateSeuratObject(sparse_matrix,
                                  project = "pbmmc",
                                  min.cells = 3,
                                  min.features = 100,
                                  names.field = 2,          # For some reason we need this last two lines for samples to be orig.ident
                                  names.delim = "\\_")

seu

# retesting with Caroline


head(seu@assays$RNA@layers$counts)
View(seu)
head(seu@meta.data)
tail(seu@meta.data)


hist(seu$nCount_RNA, breaks = 400)
hist(seu$nCount_RNA, breaks = 600, xlim = range(1, 12000))  # to see lower end better

hist(seu$nFeature_RNA, breaks = 100)
hist(seu$nFeature_RNA, breaks = 400, xlim = range(1, 3000))  # to see lower end better

Seurat::FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

Seurat::VlnPlot(seu, features = c("nCount_RNA",
                                  "nFeature_RNA"))


# mitochondrial genes
seu <- Seurat::PercentageFeatureSet(seu, 
                                    pattern = "^mt-", 
                                    #pattern = "^MT-", 
                                    col.name = "percent.mito")

# ribosomal genes
seu <- Seurat::PercentageFeatureSet(seu, 
                                    pattern = "^RP[SL]",
                                    col.name = "percent.ribo")

# hemoglobin genes (but not HBP)
seu <- Seurat::PercentageFeatureSet(seu,
                                    pattern = "^HB[^(P)]",
                                    col.name = "percent.globin")


Seurat::VlnPlot(seu, features = c("percent.mito",
                                  "percent.ribo",
                                  "percent.globin"))


Seurat::FeatureScatter(seu, 
                       feature1 = "percent.mito", 
                       feature2 = "nFeature_RNA")

# Used to much memory. 
# Below is to view the highest expression genes to determine if they should be removed

################

library(ggplot2)
library(Matrix)
library(Seurat)

most_expressed_boxplot <- function(object, ngenes = 20){
  
  # matrix of raw counts
  cts <- Seurat::GetAssayData(object, assay = "RNA", layer = "counts")
  
  # get percentage/cell
  cts <- t(cts)/colSums(cts)*100
  medians <- apply(cts, 2, median)
  
  # get top n genes
  most_expressed <- order(medians, decreasing = T)[ngenes:1]
  most_exp_matrix <- as.matrix((cts[,most_expressed]))
  
  # prepare for plotting
  most_exp_df <- stack(as.data.frame(most_exp_matrix))
  colnames(most_exp_df) <- c("perc_total", "gene")
  
  # boxplot with ggplot2
  boxplot <- ggplot(most_exp_df, aes(x=gene, y=perc_total)) +
    geom_boxplot() +
    coord_flip()
  return(boxplot)
}

most_expressed_boxplot(seu, 20)

#####################











