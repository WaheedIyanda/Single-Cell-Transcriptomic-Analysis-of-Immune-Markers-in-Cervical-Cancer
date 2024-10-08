# Loading the important Library
library(hdf5r)
library(Seurat)
library(ggplot2)
library(dplyr)

# Loading and processing our Data (scRNA-seq data from cervical cancer tissue from 10x genomic dataset)
cervcal <- Read10X_h5("Data/9k_Cervical_Cancer_scFFPE_count_filtered_feature_bc_matrix.h5")

# Creating a Seurat object
cervical_object <- CreateSeuratObject(counts = cervcal, min.cells=3, min.feature= 200)

# Quality Control
# Add Metadata: percentage of mitochondrial gene
cervical_object[["percent.mt"]] <- PercentageFeatureSet(cervical_object, pattern = "^MT-")

# Violin plot visualization of QC 
VlnPlot(cervical_object, feature = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

# Scatter plot
FeatureScatter(cervical_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method ="lm")

# Filtering cells based on QC metrics
cervical_object <- subset(cervical_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize The Data
cervical_object <- NormalizeData(cervical_object)

# Variable feature identification
cervical_object <- FindVariableFeatures(cervical_object, selection.method = "vst", nfeatures = 2000)

# Viewing top 10 most variable genes
top10 <- head(VariableFeatures(cervical_object), 10)
variable_plot <- VariableFeaturePlot(cervical_object)

# Gene labeling
LabelPoints(plot = variable_plot, points = top10, repel = TRUE)

# Dimensionality reduction
cervical_object <- ScaleData(cervical_object)
cervical_object <- RunPCA(cervical_object, features= VariableFeatures(object = cervical_object))

# Visualization
ElbowPlot(cervical_object)

# Clustering Analysis
# Finding Nearest Neighbors
cervical_object <- FindNeighbors(cervical_object, dims = 1:10)
cervical_object <- FindClusters(cervical_object, resolution = 0.5)

# Visualizing using UMAP
cervical_object <- RunUMAP(cervical_object, dims= 1:10)
DimPlot(cervical_object, reduction = "umap", label = TRUE)

# Identifying marker genes by annotation
marker <- FindAllMarkers(cervical_object, only.pos= TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Display top 10 markers
top_markers <- marker %>% group_by(cluster) %>% top_n(n= 10, wt = avg_log2FC)

# Feature plots for specific markers
FeaturePlot(cervical_object, features = c("EPCAM", "CD3D", "PTPRC", "CD68", "CD163", "MS4A1", "KRT19"))
FeaturePlot(cervical_object, features = "EPCAM")
FeaturePlot(cervical_object, features = "CD3D")
FeaturePlot(cervical_object, features = "PTPRC")
FeaturePlot(cervical_object, features = "CD68")
FeaturePlot(cervical_object, features = "CD163")
FeaturePlot(cervical_object, features = "MS4A1")
FeaturePlot(cervical_object, features = "KRT19")
FeaturePlot(cervical_object, features = "VIM")

# Heatmap visualization
DoHeatmap(cervical_object, features = top_markers$gene, size = 20) + NoLegend()

# Exploring the immunological landscape and tumor microenvironment
immune_markers <- c("CD3D", "CD4", "CD8A", "FOXP3", "CD19", "MS4A1", 
                    "CD68", "CD163", "LYZ", "CD163", "IFI30", "C1QA", "TRBC1", "ITGAX")

marker_names <- c("T-cell receptor CD3 delta", "T-cell co-receptor CD4", 
                  "T-cell co-receptor CD8A", "Forkhead box P3", 
                  "B-lymphocyte marker CD19", "B-lymphocyte marker MS4A1", 
                  "Macrophage marker CD68", "Macrophage marker CD163", 
                  "Lysozyme", "Interferon gamma-inducible protein 30", 
                  "Complement component C1q A chain", 
                  "T-cell receptor beta chain", 
                  "Integrin alpha X")

for (i in seq_along(immune_markers)) {
  plot <- FeaturePlot(cervical_object, features = immune_markers[i], pt.size = 1) +
    ggtitle(marker_names[i]) +  
    theme_minimal() +  
    theme(plot.title = element_text(hjust = 0.5)) 
  print(plot)
}

# Quantitative Analysis (Average expression of selected immune markers)
avg_expression <- AverageExpression(cervical_object, features = immune_markers)

# Convert to a data frame for easier visualization
avg_expression_df <- as.data.frame(avg_expression$RNA)
print(avg_expression_df)

# Printed average expression value organized into New DataFrame
data_matrix <- matrix(c(
  0.3996149, 0.35619481, 0.44643158, 0.447712111, 8.56844783, 0.2830056,
  5.6470090, 0.32634211, 0.40406524, 0.159537399, 0.42612160, 0.3094180,
  12.2354256, 0.84746057, 0.89134201, 0.292503841, 0.68784562, 0.3447863,
  1.4387382, 0.51210539, 0.56019587, 0.075757379, 0.75480796, 0.2001858,
  3.6737802, 0.27935001, 0.30312693, 0.173255185, 0.16804421, 0.1548049,
  3.1475932, 0.26684849, 0.30995636, 0.257535107, 5.50659722, 0.3445585,
  0.2526186, 0.56025744, 0.16882076, 0.615739552, 13.67868361, 0.2251964,
  0.9224197, 0.87096206, 0.72021423, 1.266798108, 29.32077828, 1.2352009,
  0.1746619, 0.07724593, 0.03438756, 0.008303717, 0.07520176, 0.0000000,
  1.2334088, 0.68304987, 0.46638363, 0.628477714, 35.93755228, 0.4507245,
  0.2254856, 0.35030460, 0.69904925, 0.502906109, 3.19412913, 0.1785921,
  1.0683809, 1.11655199, 1.17984637, 1.440126733, 22.37539611, 0.8083363,
  2.4217344, 0.11480146, 0.22098862, 0.058838883, 0.12847649, 0.0244810,
  
  25.83731260, 0.2178339, 0.14001310, 0.32080730, 4.744000, 0.00000000,
  0.10345553, 0.7049577, 0.08847937, 0.07418926, 1.831371, 0.19778481,
  0.28052935, 2.0411570, 0.10517468, 0.81337060, 3.525386, 0.94651792,
  0.02522787, 20.9053813, 0.00000000, 0.36350307, 1.716631, 0.00000000,
  0.01764677, 0.5969689, 0.04725959, 0.09058221, 1.367831, 0.00000000,
  8.13400719, 0.6438944, 0.09072342, 11.22677464, 5.287958, 0.04876821,
  31.16546445, 0.1191825, 0.06272044, 0.52766648, 4.174224, 0.23162027,
  64.60123233, 0.8516187, 0.
  
  ###Assigning row and column names
  rownames(data_matrix) <- c("C1QA", "CD8A", "TRBC1", "MS4A1", "CD3D", "CD4", "CD163", "LYZ", "CD19", "ITGAX", "CD68", "IFI30", "FOXP3")
  colnames(data_matrix) <- c("g0", "g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8", "g9", "g10", "g11")
  
  #### creating the dataFrame
  expression_df <- as.data.frame(data_matrix)
  
  ####Print the table
  print(expression_df)
  
  library(reshape2)
  expression_df$Gene <- rownames(expression_df)
  expression_long <- melt(expression_df, id.vars = "Gene")
  mean_expression <- expression_long %>%
    group_by(Gene, variable) %>%
    summarize(Mean = mean(value), .groups = 'drop')
  
  ggplot(mean_expression, aes(x = Gene, y = Mean, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Average Gene Expression by Group", x = "Gene", y = "Average Expression Value") +
    theme_minimal() +
    scale_fill_brewer(palette = "Paired") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  colnames(expression_long) <- c("Gene", "Group", "Expression")
  dim(expression_long)
  head(expression_long)
  
  ggplot(expression_long, aes(x = Gene, y = Expression)) +
    geom_boxplot() +
    labs(title = "Gene Expression Boxplot", x = "Gene", y = "Expression Value") +
    theme_minimal()
  
  