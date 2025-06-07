This is the code that I used for Homework 2 for Biomedical Data Analysis Course. I had to use R version 3.5.1 on the Rutgers Amarel desktop. There are some work arounds, since this version of R was not up to date. 

Load counts data from features_combined.txt
```{r}
counts <- as.matrix(read.table("features_combined.txt", header = TRUE, row.names = 1, sep = "\t"))
```
