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
  
  
}



## Selecting Cells Type    


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

HTT <- data8 %>% 
  select(HTT, cell_type_alias_label) %>% 
  filter(cell_type_alias_label == 'Exc')

summarytools::descr(HTT$HTT)

HTT_downregulated <- data8 %>% 
  select(HTT, cell_type_alias_label, cell_type) %>% 
  filter(cell_type_alias_label == 'Exc', HTT <= 6.39)

##ATXN7

ATXN7 <- data8 %>% 
  select(ATXN7, cell_type_alias_label) %>% 
  filter(cell_type_alias_label == 'Exc')

summarytools::descr(ATXN7$ATXN7)

ATXN7_downregulated <- data8 %>% 
  select(ATXN7, cell_type_alias_label, cell_type) %>% 
  filter(cell_type_alias_label == 'Exc', ATXN7 <= 6.33)

##ATXN3


ATXN3 <- data8 %>% 
  select(ATXN3, cell_type_alias_label) %>% 
  filter(cell_type_alias_label == 'Exc')

summarytools::descr(ATXN3$ATXN3)

ATXN3_downregulated <- data8 %>% 
  select(ATXN3, cell_type_alias_label, cell_type) %>% 
  filter(cell_type_alias_label == 'Exc', ATXN3 <= 2.30)

##ATN1

ATN1 <- data8 %>% 
  select(ATN1, cell_type_alias_label) %>% 
  filter(cell_type_alias_label == 'Exc')

summarytools::descr(ATN1$ATN1)

ATN1_downregulated <- data8 %>% 
  select(ATN1, cell_type_alias_label, cell_type) %>% 
  filter(cell_type_alias_label == 'Exc', ATN1 <= 1.42)

##TBP

TBP <- data8 %>% 
  select(TBP, cell_type_alias_label) %>% 
  filter(cell_type_alias_label == 'Exc')

summarytools::descr(TBP$TBP)

TBP_downregulated <- data8 %>% 
  select(TBP, cell_type_alias_label, cell_type) %>% 
  filter(cell_type_alias_label == 'Exc', TBP <= 0.33)

}
##Loop for - choose mutual cell type

for (i in HTT_downregulated$cell_type){
  if (i %in% ATXN3_downregulated$cell_type){
    if (i %in% ATXN7_downregulated$cell_type){
      if (i %in% ATN1_downregulated$cell_type){
        if (i %in% TBP_downregulated$cell_type)
          print(i)
      }
    }
  }
}




##Select mutual cell type for ATXN2, ATXN7 and HTT

for (i in HTT_downregulated$cell_type){
  if (i %in% ATXN3_downregulated$cell_type){
    if (i %in% ATXN7_downregulated$cell_type){
      print(i)
    }
  }
}


data10 <- data5 %>% 
  filter(cell_type_alias_label %in% c("Exc L2-3 RORB PTPN3")) 

data10 <- data10 %>% 
  select(HTT, ATXN7, ATXN3, cell_type_alias_label)


data10 <- data10 %>% 
  pivot_longer(cols = c(HTT, ATXN7, ATXN3))


data10$mean <- NaN

data10$mean[data10$name == 'HTT'] <- as.numeric(6.51)
data10$mean[data10$name == 'ATXN3'] <- as.numeric(2.83)
data10$mean[data10$name == 'ATXN7'] <- as.numeric(6.30)
data10$mean[data10$name == 'ATN1'] <- as.numeric(1.67)
data10$mean[data10$name == 'TBP'] <- as.numeric(0.59)


ggplot() +
  geom_col(data10, mapping = aes(x =name,y = mean),fill = 'white', color = 'red')+
  geom_col(data10, mapping = aes(x = name, y = value, fill = name))+
  facet_wrap(~cell_type_alias_label)


ggsave(filename = "Figures/exc_mutual_cell.jpeg", dpi = 600)



##Corplot heatmap exc cell type



df_inh_cell <- data8 %>% 
  filter(cell_type_alias_label == 'Exc') %>% 
  select(ATXN3, ATXN7, HTT, cell_type)

matrix_inh_cell <- as.matrix(df_inh_cell[,-4])

rownames(matrix_inh_cell) <- df_inh_cell$cell_type

matrix_inh_cell_t <- t(matrix_inh_cell)


jpeg("Figures/pheatmap_exc_cell.jpeg" , units="in", width=16, height=8, res=300)


pheatmap::pheatmap(matrix_inh_cell_t)

dev.off()

jpeg("Figures/corplot_exc_cell.jpeg" , units="in", width=25, height=8, res=300)

corrplot(matrix_inh_cell_t, is.corr=FALSE)

dev.off()
