library(tidyverse)
source('bio_function.R')
Sys.setenv(LANG = "en")
set.seed(123)


##DATA load and preparing - scRNAseq (previous normalized)

{
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

}

#Heatmap generate

jpeg("Figures/heatmap.jpeg" , units="in", width=15, height=10, res=300)


heatmap(matrix1)

dev.off()

#Correlation

corelation_genes <- cor(matrix1, method = "spearman", use = "everything")


library(corrplot)

jpeg("Figures/corplot.jpeg" , units="in", width=10, height=10, res=300)


corrplot(corelation_genes, method = 'circle')

dev.off()

matrix_cells <- t(matrix1)


##

library("pheatmap")

jpeg("Figures/pheatmap.jpeg" , units="in", width=20, height=10, res=300)

pheatmap(matrix_cells , scale = 'column')

dev.off()

###Change data into inh and ex


{
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

}

cell_corelation_gene <- cor(data6)

jpeg("Figures/corplotwithcelltype.jpeg" , units="in", width=10, height=10, res=300)

corrplot(cell_corelation_gene)

dev.off()

## Ploty    


{
data8 <- data5

data8$cell_type <- data5$cell_type_alias_label

change_value('data8', 'cell_type_alias_label', 'Exc', 'Exc')
change_value('data8', 'cell_type_alias_label', 'Inh', 'Inh')


data8 <- data8 %>% 
  filter(cell_type_alias_label %in% c('Exc', 'Inh'))
}

## HTT
{
ggplot(data = data8) +
  geom_violin(aes(cell_type_alias_label,
                  HTT, fill = cell_type_alias_label))+
  ylab("Normalized_expression") +
  xlab("Cell type")+
  ggtitle("HTT") +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))

ggsave(filename = "Figures/violinHTT.jpeg", dpi = 600)

ggplot(data = data8) +
  geom_boxplot(aes(cell_type_alias_label, 
                   HTT, fill = cell_type_alias_label))+
  ylab("Normalized_expression") +
  xlab("Cell type")+
  ggtitle("HTT") +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))

ggsave(filename = "Figures/boxHTT.jpeg", dpi = 600)

htt_anova <- aov(HTT ~ cell_type_alias_label, data = data8)

summary(htt_anova)

HTT <- data8 %>% 
  select(HTT, cell_type_alias_label) %>% 
  filter(cell_type_alias_label == 'Inh')

summarytools::descr(HTT$HTT)

HTT_upregulated <- data8 %>% 
  select(HTT, cell_type_alias_label, cell_type) %>% 
  filter(cell_type_alias_label == 'Inh', HTT >= 6.28)

##ATXN7

ggplot(data = data8) +
  geom_violin(aes(cell_type_alias_label,
                  ATXN7, fill = cell_type_alias_label))+
  ylab("Normalized_expression") +
  xlab("Cell type")+
  ggtitle("ATXN7") +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))

ggsave(filename = "Figures/violinATXN7.jpeg", dpi = 600)

ggplot(data = data8) +
  geom_boxplot(aes(cell_type_alias_label, 
                   ATXN7, fill = cell_type_alias_label))+
  ylab("Normalized_expression") +
  xlab("Cell type")+
  ggtitle("ATXN7") +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))

ggsave(filename = "Figures/boxATXN7.jpeg", dpi = 600)

ATXN7_anova <- aov(ATXN7 ~ cell_type_alias_label, data = data8)

summary(ATXN7_anova)

ATXN7 <- data8 %>% 
  select(ATXN7, cell_type_alias_label) %>% 
  filter(cell_type_alias_label == 'Inh')

summarytools::descr(ATXN7$ATXN7)

ATXN7_upregulated <- data8 %>% 
  select(ATXN7, cell_type_alias_label, cell_type) %>% 
  filter(cell_type_alias_label == 'Inh', ATXN7 >= 3.97)

##ATXN3

ggplot(data = data8) +
  geom_violin(aes(cell_type_alias_label,
                ATXN3, fill = cell_type_alias_label))+
  ylab("Normalized_expression") +
  xlab("Cell type")+
  ggtitle("ATXN3") +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))

ggsave(filename = "Figures/violinATXN3.jpeg", dpi = 600)

ggplot(data = data8) +
  geom_boxplot(aes(cell_type_alias_label, 
                   ATXN3, fill = cell_type_alias_label))+
  ylab("Normalized_expression") +
  xlab("Cell type")+
  ggtitle("ATXN3") +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))

ggsave(filename = "Figures/boxATXN3.jpeg", dpi = 600)

ATXN3_anova <- aov(ATXN3 ~ cell_type_alias_label, data = data8)

summary(ATXN3_anova)

ATXN3 <- data8 %>% 
  select(ATXN3, cell_type_alias_label) %>% 
  filter(cell_type_alias_label == 'Inh')

summarytools::descr(ATXN3$ATXN3)

ATXN3_upregulated <- data8 %>% 
  select(ATXN3, cell_type_alias_label, cell_type) %>% 
  filter(cell_type_alias_label == 'Inh', ATXN3 >= 1.59)

##ATN1

ggplot(data = data8) +
  geom_violin(aes(cell_type_alias_label,
                  ATN1, fill = cell_type_alias_label))+
  ylab("Normalized_expression") +
  xlab("Cell type")+
  ggtitle("ATN1") +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))

ggsave(filename = "Figures/violinATN1.jpeg", dpi = 600)

ggplot(data = data8) +
  geom_boxplot(aes(cell_type_alias_label, 
                   ATN1, fill = cell_type_alias_label))+
  ylab("Normalized_expression") +
  xlab("Cell type")+
  ggtitle("ATN1") +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))

ggsave(filename = "Figures/boxATN1.jpeg", dpi = 600)

ATN1_anova <- aov(ATN1 ~ cell_type_alias_label, data = data8)

summary(ATN1_anova)

ATN1 <- data8 %>% 
  select(ATN1, cell_type_alias_label) %>% 
  filter(cell_type_alias_label == 'Inh')

summarytools::descr(ATN1$ATN1)

ATN1_upregulated <- data8 %>% 
  select(ATN1, cell_type_alias_label, cell_type) %>% 
  filter(cell_type_alias_label == 'Inh', ATN1 >= 1.62)

##TBP

ggplot(data = data8) +
  geom_violin(aes(cell_type_alias_label,
                  TBP, fill = cell_type_alias_label))+
  ylab("Normalized_expression") +
  xlab("Cell type")+
  ggtitle("TBP") +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))

ggsave(filename = "Figures/violinTBP.jpeg", dpi = 600)

ggplot(data = data8) +
  geom_boxplot(aes(cell_type_alias_label, 
                   TBP, fill = cell_type_alias_label))+
  ylab("Normalized_expression") +
  xlab("Cell type")+
  ggtitle("TBP") +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))

ggsave(filename = "Figures/boxTBP.jpeg", dpi = 600)

TBP_anova <- aov(TBP ~ cell_type_alias_label, data = data8)

summary(TBP_anova)

TBP <- data8 %>% 
  select(TBP, cell_type_alias_label) %>% 
  filter(cell_type_alias_label == 'Inh')

summarytools::descr(TBP$TBP)

TBP_upregulated <- data8 %>% 
  select(TBP, cell_type_alias_label, cell_type) %>% 
  filter(cell_type_alias_label == 'Inh', TBP > 0)

}
##Loop for - choose mutual cell type

for (i in HTT_upregulated$cell_type){
  if (i %in% ATXN3_upregulated$cell_type){
    if (i %in% ATXN7_upregulated$cell_type){
      if (i %in% ATN1_upregulated$cell_type){
        if (i %in% TBP_upregulated$cell_type)
        print(i)
      }
    }
    }
}



cell_type <- data4 %>% 
  filter(cell_type_alias_label %in% c('Inh L3-5 PVALB ISG20', 'Inh L3-6 VIP ZIM2-AS1')) %>% 
  distinct()

data9 <- data5 %>% 
  filter(cell_type_alias_label %in% c('Inh L3-5 PVALB ISG20', 'Inh L3-6 VIP ZIM2-AS1')) %>% 
  distinct()

data9 <- data9 %>% 
  select(HTT, ATXN7, ATXN3, TBP, ATN1, cell_type_alias_label)


data9 <- data9 %>% 
  pivot_longer(cols = c(HTT, ATXN7, ATXN3, TBP, ATN1))


data9$mean <- NaN

data9$mean[data9$name == 'HTT'] <- as.numeric(4.29)
data9$mean[data9$name == 'ATXN3'] <- as.numeric(1.15)
data9$mean[data9$name == 'ATXN7'] <- as.numeric(3.61)
data9$mean[data9$name == 'ATN1'] <- as.numeric(1.31)
data9$mean[data9$name == 'TBP'] <- as.numeric(0.13)


ggplot() +
  geom_col(data9, mapping = aes(x = name, y = value, fill = name))+
  geom_col(data9, mapping = aes(x =name,y = mean, fill = name), color = 'black')+
  facet_grid(~cell_type_alias_label)

ggsave(filename = "Figures/all_genes_mutual_inh_overexpression.jpeg", dpi = 600)

##Select mutual cell type for ATXN2, ATXN7 and HTT

for (i in HTT_upregulated$cell_type){
  if (i %in% ATXN3_upregulated$cell_type){
    if (i %in% ATXN7_upregulated$cell_type){
          print(i)
      }
    }
}


data10 <- data5 %>% 
  filter(cell_type_alias_label %in% c("Inh L3-6 VIP ZIM2-AS1", "Inh L5 PVALB LRIG3",
                                      "Inh L2-5 PVALB HHIPL1", "Inh L5-6 PVALB KCNIP2",
                                      "Inh L5-6 SST DNAJC14", "Inh L5-6 PVALB ZFPM2-AS1",
                                      "Inh L3-5 PVALB ISG20", "Inh L1-6 SST NPY")) %>% 
  distinct()

data10 <- data10 %>% 
  select(HTT, ATXN7, ATXN3, cell_type_alias_label)


data10 <- data10 %>% 
  pivot_longer(cols = c(HTT, ATXN7, ATXN3))


data10$mean <- NaN

data10$mean[data10$name == 'HTT'] <- as.numeric(4.29)
data10$mean[data10$name == 'ATXN3'] <- as.numeric(1.15)
data10$mean[data10$name == 'ATXN7'] <- as.numeric(3.61)
data10$mean[data10$name == 'ATN1'] <- as.numeric(1.31)
data10$mean[data10$name == 'TBP'] <- as.numeric(0.13)


ggplot() +
  geom_col(data10, mapping = aes(x = name, y = value, fill = name))+
  geom_col(data10, mapping = aes(x =name,y = mean, fill = name), color = 'black')+
  facet_wrap(~cell_type_alias_label)


ggsave(filename = "Figures/ATXN3_ATXN7_HTT_mutual_inh_overexpression.jpeg", dpi = 600)


##Corplot heatmap all cell type

matrix_last <- data5 %>% 
  select(ATXN3, ATXN7, HTT)

matrix_last <- as.matrix(matrix_last)

rownames(matrix_last) <- data5$cell_type_alias_label

matrix_last_t <- t(matrix_last)


jpeg("Figures/pheatmap_ATXN3_ATXN7_HTT.jpeg" , units="in", width=50, height=15, res=300)


pheatmap::pheatmap(matrix_last_t)

dev.off()

jpeg("Figures/corplot_all_cell_ATXN3_ATXN7_HTT.jpeg" , units="in", width=50, height=15, res=300)

corrplot(matrix_last_t, is.corr=FALSE)

dev.off()


##Corplot heatmap inh cell type



df_inh_cell <- data8 %>% 
  filter(cell_type_alias_label == 'Inh') %>% 
  select(ATXN3, ATXN7, HTT, cell_type)

matrix_inh_cell <- as.matrix(df_inh_cell[,-4])

rownames(matrix_inh_cell) <- df_inh_cell$cell_type

matrix_inh_cell_t <- t(matrix_inh_cell)


jpeg("Figures/pheatmap_inh_cell_ATXN3_ATXN7_HTT.jpeg" , units="in", width=16, height=8, res=300)


pheatmap::pheatmap(matrix_inh_cell_t)

dev.off()

jpeg("Figures/corplot_inh_cell_ATXN3_ATXN7_HTT.jpeg" , units="in", width=25, height=8, res=300)

corrplot(matrix_inh_cell_t, is.corr=FALSE)

dev.off()
