## Trajectory Monocle3

seu <- readRDS("my_objects/seu_day2-4.rds")
library(monocle3)

# create gene metadata data.frame
feature_names <- as.data.frame(rownames(seu))
rownames(feature_names) <- rownames(seu)
colnames(feature_names) <- "gene_short_name"

# initiate monocle object from seurat count table 
seu_monocl <- monocle3::new_cell_data_set(Seurat::GetAssayData(seu,
                                                               layer = "counts"),
                                          cell_metadata = seu@meta.data,
                                          gene_metadata = feature_names)

?preprocess_cds

## Preprocessing by normalizing and scaling the raw counts again ?
seu_monocl <- monocle3::preprocess_cds(seu_monocl)

# Check elbow plot
monocle3::plot_pc_variance_explained(seu_monocl)


# UMAP through monocl
seu_monocl <- monocle3::reduce_dimension(seu_monocl, reduction_method = "UMAP")

monocle3::plot_cells(seu_monocl, 
                     color_cells_by = "RNA_snn_res.0.3", 
                     cell_size = 1, 
                     show_trajectory_graph = FALSE)


monocle3::plot_cells(seu_monocl, genes = "Cd79a", 
                     cell_size = 1,
                     show_trajectory_graph = FALSE,
                     scale_to_range = FALSE)



monocle3::plot_cells(seu_monocl, 
                     color_cells_by = "SingleR_annot", 
                     cell_size = 1, 
                     show_trajectory_graph = FALSE)


## Cluster Cells

seu_monocl <- monocle3::cluster_cells(seu_monocl, resolution=0.00025)
monocle3::plot_cells(seu_monocl, label_cell_groups = F)

seu_monocl <- monocle3::learn_graph(seu_monocl)
monocle3::plot_cells(seu_monocl)

monocle3::plot_cells(seu_monocl, genes = c("Kit", "Ano1"),
                     show_trajectory_graph = FALSE, 
                     cell_size = 0.7, group_label_size = 4)



# Select inital cells
seu_monocl <- monocle3::order_cells(seu_monocl)

monocle3::plot_cells(seu_monocl,
                     color_cells_by = "pseudotime",
                     label_cell_groups=F,
                     label_leaves=F,
                     label_branch_points=FALSE,
                     graph_label_size=1.5, cell_size = 1)



# Select an entire area to extract
seuB <- choose_cells(seu_monocl) ## Interactive



plot_cells(seuB, show_trajectory_graph = FALSE, cell_size = 1
           ,
           color_cells_by = "pseudotime"
           )

# Cells aren't colored in a way that allows them to be grouped.
## not working

# https://github.com/cole-trapnell-lab/monocle3/issues/171


root_group = colnames(seu_monocl)[clusters(seu_monocl) == 6]
cds = order_cells(seu_monocl, root_cells = root_group)
plot_cells(cds, color_cells_by = "pseudotime")


seuB <- cds
pr_test <- graph_test(seuB, 
                      cores=4, 
                      neighbor_graph = "principal_graph")

# order by test statistic
pr_test <- pr_test[order(pr_test$morans_test_statistic, 
                         decreasing = TRUE),]

View(pr_test)


goi <- c("Cd34", "Ms4a1", "Igll1", "Igll5", 
         "Mki67", "Cks2")
plot_cells(seuB, label_cell_groups=FALSE, genes = goi,
           show_trajectory_graph=FALSE, cell_size = 1)


seuB@colData$monocle_cluster <- clusters(seuB)

plot_genes_in_pseudotime(subset(seuB, 
                                rowData(seuB)$gene_short_name %in% goi),
                         min_expr=0.5, color_cells_by = "monocle_cluster")
