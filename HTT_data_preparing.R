library(tidyverse)
source('bio_function.R')
Sys.setenv(LANG = "en")
set.seed(123)


##DATA load and preparing - scRNAseq (previous normalized)

data2 <- read.csv('Data/trimmed_means.csv')
metadata <- read.csv('Data/metadata.csv')


data3 <- data2 %>% 
  select(cluster_label, HTT, ATXN7, ATXN1, ATXN2, ATXN3, TBP, ATN1, AR)

data4 <- data3 %>% 
  left_join(metadata, by = 'cluster_label')

data5 <- data4 %>% 
  select(AR, HTT, ATXN7, ATXN1, ATXN2, ATXN3, TBP, ATN1, cell_type_alias_label) %>% 
  distinct()


matrix1 <- as.matrix(data5[,-9])

rownames(matrix1) <- data5$cell_type_alias_label

#Heatmap generate

heatmap(matrix1)


#Correlation

corelation_genes <- cor(matrix1, method = "spearman", use = "everything")


library(corrplot)

corrplot(corelation_genes, method = 'circle')

matrix_cells <- t(matrix1)


##

library("pheatmap")

pheatmap(matrix_cells , scale = 'column')


###Change data into inh and ex



data6 <- data5

source("bio_function.R")

change_value('data6', 'cell_type_alias_label', 'Exc', 'Exc')
change_value('data6', 'cell_type_alias_label', 'Inh', 'Inh')

matrix_inf_exc <- as.matrix(data6[,-9])

rownames(matrix_inf_exc) <- as.factor(data6$cell_type_alias_label)

data6 <- data6 %>% 
  filter(cell_type_alias_label %in% c('Exc', 'Inh'))

data6$cell_type <- as.numeric(as.factor(data6$cell_type_alias_label))

data6 <- data6 %>% 
  select(-cell_type_alias_label)

cell_corelation_gene <- cor(data6)

corrplot(cell_corelation_gene)


## Ploty    

data7 <- data5

source("bio_function.R")

change_value('data7', 'cell_type_alias_label', 'Exc', 'Exc')
change_value('data7', 'cell_type_alias_label', 'Inh', 'Inh')

matrix_inf_exc <- as.matrix(data6[,-9])

data7 <- data7 %>% 
  filter(cell_type_alias_label %in% c('Exc', 'Inh'))

## HTT

ggplot(data = data7) +
  geom_violin(aes(cell_type_alias_label,
                  HTT, fill = cell_type_alias_label))+
  ylab("Normalized_expression") +
  xlab("Cell type")+
  ggtitle("HTT") +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))

ggplot(data = data7) +
  geom_boxplot(aes(cell_type_alias_label, 
                   HTT, fill = cell_type_alias_label))+
  ylab("Normalized_expression") +
  xlab("Cell type")+
  ggtitle("HTT") +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))

htt_anova <- aov(HTT ~ cell_type_alias_label, data = data7)

summary(htt_anova)

##ATXN7

ggplot(data = data7) +
  geom_violin(aes(cell_type_alias_label,
                  ATXN7, fill = cell_type_alias_label))+
  ylab("Normalized_expression") +
  xlab("Cell type")+
  ggtitle("ATXN7") +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))

ggplot(data = data7) +
  geom_boxplot(aes(cell_type_alias_label, 
                   ATXN7, fill = cell_type_alias_label))+
  ylab("Normalized_expression") +
  xlab("Cell type")+
  ggtitle("ATXN7") +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))

ATXN7_anova <- aov(ATXN7 ~ cell_type_alias_label, data = data7)

summary(ATXN7_anova)

##ATXN3

ggplot(data = data7) +
  geom_violin(aes(cell_type_alias_label,
                ATXN3, fill = cell_type_alias_label))+
  ylab("Normalized_expression") +
  xlab("Cell type")+
  ggtitle("ATXN3") +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))

ggplot(data = data7) +
  geom_boxplot(aes(cell_type_alias_label, 
                   ATXN3, fill = cell_type_alias_label))+
  ylab("Normalized_expression") +
  xlab("Cell type")+
  ggtitle("ATXN3") +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))

ATXN3_anova <- aov(ATXN3 ~ cell_type_alias_label, data = data7)

summary(ATXN3_anova)

##ATN1

ggplot(data = data7) +
  geom_violin(aes(cell_type_alias_label,
                  ATN1, fill = cell_type_alias_label))+
  ylab("Normalized_expression") +
  xlab("Cell type")+
  ggtitle("ATN1") +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))

ggplot(data = data7) +
  geom_boxplot(aes(cell_type_alias_label, 
                   ATN1, fill = cell_type_alias_label))+
  ylab("Normalized_expression") +
  xlab("Cell type")+
  ggtitle("ATN1") +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))

ATN1_anova <- aov(ATN1 ~ cell_type_alias_label, data = data7)

summary(ATN1_anova)

##TBP

ggplot(data = data7) +
  geom_violin(aes(cell_type_alias_label,
                  TBP, fill = cell_type_alias_label))+
  ylab("Normalized_expression") +
  xlab("Cell type")+
  ggtitle("TBP") +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))

ggplot(data = data7) +
  geom_boxplot(aes(cell_type_alias_label, 
                   TBP, fill = cell_type_alias_label))+
  ylab("Normalized_expression") +
  xlab("Cell type")+
  ggtitle("TBP") +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))

TBP_anova <- aov(TBP ~ cell_type_alias_label, data = data7)

summary(TBP_anova)



