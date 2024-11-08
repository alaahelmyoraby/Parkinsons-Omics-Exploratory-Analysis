# **Parkinsons-Omics-Exploratory-Analysis**

This repository contains an R script for performing a variety of data analysis tasks on transcriptomic data. The primary goal is to analyze gene expression data, perform exploratory analysis, conduct PCA, and visualize the results. Additionally, heatmaps are generated to visualize the expression of the most variable genes.

## **Files Included**
- `Parkinson_exp.csv`: Gene expression data file (rows represent genes and columns represent samples).
- `Parkinson_phenotable.csv`: Phenotypic data file (contains sample annotations such as disease type).

## **Required Libraries**
To run the analysis, you'll need to install the following R libraries. You can install them via the `install.packages()` function or `Bioconductor` for certain libraries:

```r
install.packages(c("dplyr", "ggplot2", "matrixStats", "ComplexHeatmap", "circlize", "rgl", "ggfortify"))
```

### **Libraries Overview**
- **dplyr**: Data manipulation.
- **ggplot2**: Data visualization.
- **matrixStats**: Statistical calculations like row variance.
- **ComplexHeatmap**: Heatmap generation.
- **circlize**: Color functions for heatmaps.
- **rgl**: 3D visualization (used for 3D PCA).
- **ggfortify**: Simplified plotting of PCA results.

## **Steps for Running the Analysis**

### 1. **Load the Data**
- The expression data (`Parkinson_exp.csv`) and phenotype data (`Parkinson_phenotable.csv`) are read into R. Both files are tab-separated and should be placed in the same directory as the script.

```r
exp_file <- read.csv("Parkinson_exp.csv", header = TRUE, row.names = 1, sep = "\t")
annot_file <- read.csv("Parkinson_phenotable.csv", header = TRUE, sep = "\t")
```

### 2. **Data Quality Check**
- Missing values are checked in both the gene expression and phenotype data.
- Basic exploratory analysis is performed to check the structure of the data (`str()`) and summary statistics (`summary()`).

```r
sum(is.na(exp_file))  # Check missing values in expression data
sum(is.na(annot_file))  # Check missing values in phenotype data
```

### 3. **Annotation Validation**
- Ensure that all sample IDs in the expression data match those in the phenotype data.
  
```r
all(colnames(exp_file) %in% annot_file$sample.id)
```

### 4. **Merge Expression and Annotation Data**
- The expression data is transposed and merged with the phenotype data based on the sample ID.
  
```r
exp_t_df <- as.data.frame(t(exp_file))
exp_t_df$sample.id <- row.names(exp_t_df)
merged <- merge(exp_t_df, annot_file, by = "sample.id")
```

### 5. **Exploratory Data Visualizations**
- **Boxplot**: Visualizes the distribution of expression values across samples.
- **Density Plot**: Shows the density distribution of expression levels across samples.
  
```r
boxplot(as.matrix(exp_file), main = "Boxplot of Expression Data", xlab = "Samples", ylab = "Expression Levels")
ggplot(dat, aes(x = values, fill = ind)) + 
  geom_density(alpha = 0.3) + 
  labs(title = "Density Plot of Expression Data", x = "Expression Level", y = "Density")
```

### 6. **Principal Component Analysis (PCA)**
- PCA is performed on the expression data to reduce dimensionality and visualize the variance across samples.
  
#### **PCA Plot (Manual)**
- In the original script, PCA is performed manually using `prcomp()` and plotted using `ggplot2`.
  
```r
PCA_res <- prcomp(PCA_mat, center = TRUE, scale. = TRUE)
PCA_df <- as.data.frame(PCA_res$x)
PCA_df$sample.type <- merged$sample.type
ggplot(PCA_df, aes(x = PC1, y = PC2, color = sample.type)) +
  geom_point() +
  labs(title = "PCA Plot", x = "PC1", y = "PC2") +
  theme_dark()
```

#### **Enhanced PCA Plot (Using `ggfortify`)**
- The enhanced version uses `ggfortify`'s `autoplot()` function for a more streamlined PCA visualization, which automatically handles color and labeling.

```r
autoplot(PCA_res, data = merged, colour = 'sample.type', frame = TRUE, frame.type = "norm") +
  labs(title = "PCA Plot with ggfortify", x = "PC1", y = "PC2")
```

### 7. **3D PCA Visualization**
- A 3D PCA plot is generated using the `rgl` package to visualize the first three principal components.
  
```r
pca_3d <- prcomp(t(as.matrix(exp_file)), scale. = TRUE)
mycolors <- ifelse(merged$sample.type == "PD", "red", "lightgreen")
plot3d(pca_3d$x[, 1:3], col = mycolors, size = 12, type = "s", main = "3D PCA Plot")
```

### 8. **Identifying Most Variable Genes**
- The most variable genes are identified by calculating the row variance of the expression matrix. A list of the top 10 most variable genes is generated.

```r
genevariance <- rowVars(exp_mat)
sorted_genevariance <- sort(genevariance, decreasing = TRUE)
top_10 <- names(sorted_genevariance)[1:10]
top_10_var <- sorted_genevariance[1:10]
```

### 9. **Heatmap of Most Variable Genes**
- A heatmap is generated for the top 100 most variable genes using the `ComplexHeatmap` package.
  
```r
top_gene_names <- names(sorted_genevariance)[1:100]
top_exp_data <- exp_mat[top_gene_names, ]
heatmap_absexp <- Heatmap(as.matrix(top_exp_data), 
                          name = "Expression", 
                          cluster_rows = TRUE, 
                          cluster_columns = TRUE)
draw(heatmap_absexp, height = unit(15, "cm"))
```

### 10. **Z-Score Normalization and Heatmap**
- Expression data is normalized using Z-scores and a heatmap is generated with sample annotations.
  
```r
scaled_exp_data <- t(scale(t(top_exp_data)))
column_ha <- HeatmapAnnotation(Sample.type = merged$sample.type)
col_fun <- colorRamp2(c(min(scaled_exp_data, na.rm = TRUE), 
                        mean(scaled_exp_data, na.rm = TRUE), 
                        max(scaled_exp_data, na.rm = TRUE)), 
                      c("navy", "white", "pink"))
Heatmap(scaled_exp_data, name = "Z-Score", top_annotation = column_ha, 
        col = col_fun, show_row_names = TRUE, show_column_names = TRUE)
```
