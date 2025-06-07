**This is the code that I used for Homework 2 for Biomedical Data Analysis Course. I had to use R version 3.5.1 on the Rutgers Amarel desktop. There are some work arounds, since this version of R was not up to date.**

Load counts data from features_combined.txt
```{r}
counts <- as.matrix(read.table("features_combined.txt", header = TRUE, row.names = 1, sep = "\t"))
```

Show genes x samples to check to see if it loads proverly 
```{r}
dim(counts)

head(colnames(counts))
```

Load the metadata file (meta_data.txt) as the colData
```{r}
col_data <- read.table("meta_data.txt", header = TRUE, row.names = 1, sep = "", stringsAsFactors = FALSE)
```

Align metadata with counts 
```{r}
col_data <- col_data[colnames(counts), ]
```

Get a summarizedexperiment-like data set saved (with this version of R, we cannot load the summarizedexperiment package, so this is a workaround)
```{r}
se_like <- list(
  assays = list(counts = counts),
  colData = col_data
)
```

Remove "tb_hiv_art" samples from colData
```{r}
col_data_filtered <- col_data[col_data$Disease != "tb_hiv_art", , drop = FALSE]

counts_filtered <- counts[, rownames(col_data_filtered)]

all(colnames(counts_filtered) == rownames(col_data_filtered))  
```

Remove 0 counts across all samples
```{r}
counts_filtered <- counts_filtered[rowSums(counts_filtered) > 0, ]
```
