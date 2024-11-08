# Parkinsons-Omics-Exploratory-Analysis

This project involves exploratory data analysis of a Parkinson's disease dataset, specifically focusing on gene expression data. The workflow includes data loading, quality checks, PCA analysis, identifying variable genes, and visualization through heatmaps.

## Requirements

This analysis requires the following R packages:
- `dplyr`
- `ggplot2`
- `matrixStats`
- `ComplexHeatmap`
- `circlize`
- `rgl` (for 3D plotting)
- `ggfortify` (for enhanced PCA plotting)
- `scatterplot3d` (for additional 3D PCA options)

Make sure these libraries are installed in your R environment:
```r
install.packages(c("dplyr", "ggplot2", "matrixStats", "ComplexHeatmap", "circlize", "rgl", "ggfortify", "scatterplot3d"))
```

## Steps in the Analysis

1. **Set Working Directory**  
   Set the working directory to the folder containing the data files. Adjust the path as necessary:
   ```r
   setwd("D:/EG-COMPBIO/MODA/lec1/Assignment")
   ```

2. **Load the Data**  
   Load gene expression data and phenotype data. We use `read.csv` to import the expression and annotation files:
   ```r
   exp_file <- read.csv("Parkinson_exp.csv", header = TRUE, row.names = 1, sep = "\t")
   annot_file <- read.csv("Parkinson_phenotable.csv", header = TRUE, sep = "\t")
   ```

3. **Data Quality Check**  
   Check for missing values in the datasets:
   ```r
   sum(is.na(exp_file))   # Check missing values in expression data
   sum(is.na(annot_file)) # Check missing values in phenotype data
   ```

4. **Exploratory Analysis**  
   Summarize the structure of the datasets:
   ```r
   str(exp_file)
   str(annot_file)
   summary(exp_file)
   summary(annot_file)
   ```

5. **Annotation Validation**  
   Ensure that all samples are correctly annotated by checking if sample IDs in the expression data match those in the metadata file:
   ```r
   all(colnames(exp_file) %in% annot_file$sample.id)
   ```

6. **Check Annotation and Merge Data**  
   This step merges the expression data and annotation files to ensure that all samples have been correctly labeled as "tumor" or "normal." It's essential to confirm that the sample IDs are consistent between the two files.
   ```r
   exp_t_df <- as.data.frame(t(exp_file))
   exp_t_df$sample.id <- row.names(exp_t_df)
   merged <- merge(exp_t_df, annot_file, by = "sample.id")
   table(merged$sample.type)  # Check sample types (tumor or normal)
   ```

7. **Boxplot and Density Plot**  
   Create a boxplot and density plot for the expression data to visualize its distribution:
   ```r
   boxplot(as.matrix(exp_file), main = "Boxplot of Expression Data", xlab = "Samples", ylab = "Expression Levels")
   dat <- stack(as.data.frame(exp_file))
   ggplot(dat, aes(x = values, fill = ind)) + 
     geom_density(alpha = 0.3) + 
     labs(title = "Density Plot of Expression Data", x = "Expression Level", y = "Density")
   ```

8. **PCA Analysis**  
   Perform Principal Component Analysis (PCA) for dimensionality reduction and visualize in 2D and 3D.

   - **2D PCA:**
     ```r
     PCA_mat <- merged[, sapply(merged, is.numeric)]
     PCA_res <- prcomp(PCA_mat, center = TRUE, scale. = TRUE)
     PCA_df <- as.data.frame(PCA_res$x)
     PCA_df$sample.type <- merged$sample.type
     
     ggplot(PCA_df, aes(x = PC1, y = PC2, color = sample.type)) +
       geom_point() +
       labs(title = "PCA Plot", x = "PC1", y = "PC2") +
       theme_dark()
     ```

   - **3D PCA:**
     ```r
     pca_3d <- prcomp(t(as.matrix(exp_file)), scale. = TRUE)
     mycolors <- ifelse(merged$sample.type == "PD", "red", "lightgreen")
     plot3d(pca_3d$x[,1:3], col = mycolors, size = 12, type = "s", main = "3D PCA Plot")
     ```

9. **Identify Most Variable Genes**  
   Identify the top 10 genes with the highest variance to understand which genes show the most variation across samples:
   ```r
   exp_mat <- as.matrix(exp_file)
   genevariance <- rowVars(exp_mat)
   names(genevariance) <- rownames(exp_mat)
   sorted_genevariance <- sort(genevariance, decreasing = TRUE)
   top_10 <- names(sorted_genevariance)[1:10]
   top_10_var <- sorted_genevariance[1:10]
   top_10_gen_val <- data.frame(Genes = top_10, Variance = top_10_var)
   print(top_10_gen_val)
   ```

10. **Heatmap of Top 100 Most Variable Genes**  
   Plot a heatmap for the top 100 most variable genes to observe their expression across samples.

   - **Heatmap with Absolute Expression:**
     ```r
     top_gene_names <- names(sorted_genevariance)[1:100]
     top_exp_data <- exp_mat[top_gene_names, ]
     heatmap_absexp <- Heatmap(as.matrix(top_exp_data), 
                               name = "Expression",
                               cluster_rows = TRUE, 
                               cluster_columns = TRUE,
                               show_row_names = TRUE, 
                               show_column_names = TRUE)
     draw(heatmap_absexp, height = unit(15, "cm"))
     ```

   - **Z-Score Scaled Heatmap with Annotations:**
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

## Summary
This code provides a comprehensive exploratory analysis pipeline for gene expression data. By following the steps above, one can visualize sample distributions, perform dimensionality reduction, and analyze gene variability, helping uncover patterns in complex biological data.

