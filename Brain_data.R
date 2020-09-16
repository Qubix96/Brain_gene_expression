library(tidyverse)
source('bio_function.R')
Sys.setenv(LANG = "en")
set.seed(123)

#Data preparing

brain_expression <- read.csv('Brain_dev_data/expression_matrix.csv', header = F)
brain_rows <- read.csv('Brain_dev_data/rows_metadata.csv') 
brain_cols <- read.csv('Brain_dev_data/columns_metadata.csv')

check <- brain_expression
{

#DF data

brain_expression$V1 <- brain_rows$gene_symbol

brain_expression <- brain_expression %>% 
  filter(V1 %in% c('HTT', 'ATXN7', 'ATXN1', 'ATXN2', 'ATXN3', 'TBP', 'ATN1', 'AR'))


##Matrix data - mean(genes)


brain_expression_matrix_all <- as.matrix(brain_expression[,-1])

rownames(brain_expression_matrix_all) <- brain_expression$V1
colnames(brain_expression_matrix_all) <- brain_cols$structure_name


df <- as.data.frame(brain_expression_matrix_all)

df$V1 <- brain_expression$V1



aggregate_all('df' , '[1:524]', 'V1')


matrix_genes_all <- as.matrix(df[,-1])
rownames(matrix_genes_all) <- df$Group.1


##Matrix data mean(brain_parts)

t <- t(matrix_genes_all)

df <- as.data.frame(t, col.names = F)

df$col <- brain_cols$structure_name

aggregate_all('df', '[1:8]', 'col')

matrix_genes_all <- as.matrix(df[,-1])

rownames(matrix_genes_all) <- df$Group.1

matrix_genes_all <- t(matrix_genes_all)
}

##Analysis

#Heatmaps 

jpeg("Brain_dev_figures/heatmap_all_brainparts.jpeg" , units="in", width=15, height=20, res=300)


heatmap(matrix_genes_all)

dev.off()


jpeg("Brain_dev_figures/pheatmap_all_brainparts.jpeg" , units="in", width=20, height=8, res=300)

pheatmap::pheatmap(matrix_genes_all)

dev.off()


jpeg("Brain_dev_figures/corplot_all_brainparts.jpeg" , units="in", width=20, height=8, res=300)

corrplot::corrplot(matrix_genes_all, is.corr=FALSE)

dev.off()


##############################################################################

#DF data - HTT, ATXN7, ATXN3

{
brain_expression$V1 <- brain_rows$gene_symbol

brain_expression <- brain_expression %>% 
  filter(V1 %in% c('HTT', 'ATXN7', 'ATXN3'))


##Matrix data- mean(HTT, ATXN7, ATXN3)


brain_expression_matrix <- as.matrix(brain_expression[,-1])
  
  rownames(brain_expression_matrix) <- brain_expression$V1
  colnames(brain_expression_matrix) <- brain_cols$structure_name
  
  
  df <- as.data.frame(brain_expression_matrix)
  
  df$V1 <- brain_expression$V1
  
  
  
  aggregate_all('df' , '[1:524]', 'V1')
  
  
  matrix_genes <- as.matrix(df[,-1])
  rownames(matrix_genes) <- df$Group.1
  
  
##Matrix data mean(brain_parts)

t <- t(matrix_genes)

df <- as.data.frame(t, col.names = F)

df$col <- brain_cols$structure_name

aggregate_all('df', '[1:3]', 'col')

matrix_genes <- as.matrix(df[,-1])

rownames(matrix_genes) <- df$Group.1

matrix_genes <- t(matrix_genes)
}

##Analysis

#Heatmaps 

jpeg("Brain_dev_figures/heatmap_all_brainparts_HTT_ATXN7_ATXN3.jpeg" , units="in", width=15, height=20, res=300)


heatmap(matrix_genes)

dev.off()


jpeg("Brain_dev_figures/pheatmap_all_brainparts_HTT_ATXN7_ATXN3.jpeg" , units="in", width=20, height=8, res=300)

pheatmap::pheatmap(matrix_genes)

dev.off()


jpeg("Brain_dev_figures/corplot_all_brainparts_HTT_ATXN7_ATXN3.jpeg" , units="in", width=20, height=8, res=300)

corrplot::corrplot(matrix_genes, is.corr=FALSE)

dev.off()
