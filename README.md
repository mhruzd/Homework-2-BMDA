**This is the code that I used for Homework 2 for Biomedical Data Analysis Course. I had to use R version 3.5.1 on the Rutgers Amarel desktop. There are some work arounds, since this version of R was not up to date.**


*Align and filter the data:*

Load counts data from `features_combined.txt`
```{r}
counts <- as.matrix(read.table("features_combined.txt", header = TRUE, row.names = 1, sep = "\t"))
```

Show genes x samples to check to see if it loads proverly 
```{r}
dim(counts)

head(colnames(counts))
```

Load the metadata file (`meta_data.txt`) as the `colData`
```{r}
col_data <- read.table("meta_data.txt", header = TRUE, row.names = 1, sep = "", stringsAsFactors = FALSE)
```

Align metadata with counts 
```{r}
col_data <- col_data[colnames(counts), ]
```

Get a `summarizedexperiment`-like data set saved (with this version of R, we cannot load the `summarizedexperiment` package, so this is a workaround)
```{r}
se_like <- list(
  assays = list(counts = counts),
  colData = col_data
)
```

Remove `tb_hiv_art` samples from `colData`
```{r}
col_data_filtered <- col_data[col_data$Disease != "tb_hiv_art", , drop = FALSE]

counts_filtered <- counts[, rownames(col_data_filtered)]

all(colnames(counts_filtered) == rownames(col_data_filtered))
```

Remove 0 counts across all samples
```{r}
counts_filtered <- counts_filtered[rowSums(counts_filtered) > 0, ]
```


*Generate dimension reduction plots applying PCA*

Convert to CPM manually 
```{r}
lib_sizes <- colSums(counts_filtered)

cpm_matrix <- t(t(counts_filtered) / lib_sizes * 1e6)
```

Log-transform with pseudocount
```{r}
log_cpm <- log2(cpm_matrix + 1)
```

Create another (final) `SummarizedExperiment` with the filtered data and the logCPM
```{r}
se_final <- list(
  assays = list(
    counts = counts_filtered,
    logCPM = log_cpm
  ),
  colData = col_data_filtered
)
```

Transpose the sampels as rows and genes as columns since PCA expects the variables to be presented as columns (otherwise will not run)
```{r}
logCPM_t <- t(se_final$assays$logCPM)
```

Run PCA using `procomp`
```{r}
pca_res <- prcomp(logCPM_t, center = TRUE, scale. = TRUE)

summary(pca_res)
```


*Run MDS from PCA Distances*

Calculate distances between the samples in the PCA space (using the first 10 PCs)
```{r}
pca_dist <- dist(pca_res$x[,1:10])
```

Run classic MDS (since we cannot use `UMAP` with this version of R)
```{r}
mds_res <- cmdscale(pca_dist, k = 2)
```

Plot MDS colored by disease status
```{r}
plot(mds_res[,1], mds_res[,2],
     col = as.factor(se_final$colData$Disease),
     pch = 19,
     xlab = "MDS 1",
     ylab = "MDS 2",
     main = "MDS on PCA distances")

legend("topright", legend = unique(se_final$colData$Disease),
       col = 1:length(unique(se_final$colData$Disease)), pch = 19)
```

*Conduct the same analysis on the log CPM values*

Load library for `limma` and assign vairables to run the analysis 
```{r}
library(limma)

counts <- se_final$assays$counts

meta <- se_final$colData
```

Filter for `TB-HIV` and `HIV-only` disease status
```{r}
keep <- meta$Disease %in% c("tb_hiv", "hiv_only")
counts_sub <- counts[, keep]
meta_sub <- meta[keep, , drop = FALSE]
```

Make the group variable a factor with `hiv_only` as the reference 
```{r}
group <- factor(meta_sub$Disease, levels = c("hiv_only", "tb_hiv"))
```

Check the sample alignment to make sure that everything looks good 
```{r}
stopifnot(all(colnames(counts_sub) == rownames(meta_sub)))
```

Manually calculate the library sises 
```{r}
lib_sizes <- colSums(counts_sub)
```

Crteate a logCPM matrix manually (like before)
```{r}
cpm_matrix <- t(t(counts_sub) / lib_sizes * 1e6)

log_cpm <- log2(cpm_matrix + 1)
```

Use `limma` on the logCPM matrix 
```{r}
design <- model.matrix(~ group)
fit <- lmFit(log_cpm, design)
fit <- eBayes(fit)
```

Get the top 50 genes
```{r}
res <- topTable(fit, coef = 2, number = Inf, sort.by = "P")

top50 <- head(res, 50)

write.csv(top50, file = "top50_limma_DE_genes.csv")
```


*Create a heatmap plot of the results*

Pull the top 50 genes from the logCPM matrix that was just calculated with `limma` in the prior step
```{r}
top_genes <- rownames(top50)

log_cpm <- se_final$assays$logCPM

log_cpm_top50 <- log_cpm[top_genes, ]
```

Get the disease status for the columns
```{r}
group <- se_final$colData$Disease

group <- group[colnames(log_cpm_top50)]
```

Assign colors to the groups for the different dieases statuses
```{r}
group_colors <- as.character(factor(group, 
                         levels = c("hiv_only", "tb_hiv", "tb_hiv_art"),
                         labels = c("blue", "red", "purple")))
```

Get the matrix of colors to plot as a color bar 
```{r}
col_annotation <- matrix(group_colors, nrow = 1)

layout(matrix(c(1, 2), nrow = 2), heights = c(1, 10))

par(mar = c(0, 5, 2, 2))  
image(1:ncol(log_cpm_top50), 1, col_annotation,
      col = c("blue", "red", "purple"),
      axes = FALSE, xlab = "", ylab = "")
title("Disease Status")
```

Plot heatmap using the `heatmap` package
```{r}
par(mar = c(0, 5, 2, 2))  
image(1:ncol(log_cpm_top50), 1, matrix(group_colors, nrow = 1),
      col = c("blue", "red", "purple"),
      axes = FALSE, xlab = "", ylab = "")
title("Disease Status")

par(mar = c(5, 5, 2, 2))
heatmap(as.matrix(log_cpm_top50),
        scale = "row",
        Colv = NA,  # skip clustering samples (optional)
        Rowv = TRUE,
        col = heat.colors(256),
        labCol = FALSE,  # hide sample names to keep plot clean
        margins = c(5, 10))

dev.off()
```
