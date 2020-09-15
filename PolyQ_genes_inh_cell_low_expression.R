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

## Select low expression cells   


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
  filter(cell_type_alias_label == 'Inh')

summarytools::descr(HTT$HTT)

HTT_upregulated <- data8 %>% 
  select(HTT, cell_type_alias_label, cell_type) %>% 
  filter(cell_type_alias_label == 'Inh', HTT <= 3.12)

##ATXN7

ATXN7 <- data8 %>% 
  select(ATXN7, cell_type_alias_label) %>% 
  filter(cell_type_alias_label == 'Inh')

summarytools::descr(ATXN7$ATXN7)

ATXN7_upregulated <- data8 %>% 
  select(ATXN7, cell_type_alias_label, cell_type) %>% 
  filter(cell_type_alias_label == 'Inh', ATXN7 <= 2.56)

##ATXN3

ATXN3 <- data8 %>% 
  select(ATXN3, cell_type_alias_label) %>% 
  filter(cell_type_alias_label == 'Inh')

summarytools::descr(ATXN3$ATXN3)

ATXN3_upregulated <- data8 %>% 
  select(ATXN3, cell_type_alias_label, cell_type) %>% 
  filter(cell_type_alias_label == 'Inh', ATXN3 <= 0.63)

##ATN1

ATN1 <- data8 %>% 
  select(ATN1, cell_type_alias_label) %>% 
  filter(cell_type_alias_label == 'Inh')

summarytools::descr(ATN1$ATN1)

ATN1_upregulated <- data8 %>% 
  select(ATN1, cell_type_alias_label, cell_type) %>% 
  filter(cell_type_alias_label == 'Inh', ATN1 <= 0.90)

##TBP

TBP <- data8 %>% 
  select(TBP, cell_type_alias_label) %>% 
  filter(cell_type_alias_label == 'Inh')

summarytools::descr(TBP$TBP)

TBP_upregulated <- data8 %>% 
  select(TBP, cell_type_alias_label, cell_type) %>% 
  filter(cell_type_alias_label == 'Inh', TBP <= 0)
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
  filter(cell_type_alias_label %in% 
           c("Inh L1 PAX6 CHRFAM7A", "Inh L1-3 SST FAM20A",
             "Inh L1-3 VIP HSPB6", "Inh L3-5 VIP TAC3", 
             "Inh L1-6 LAMP5 NES")) %>% 
  distinct()

data9 <- data5 %>% 
  filter(cell_type_alias_label %in% c("Inh L1 PAX6 CHRFAM7A", "Inh L1-3 SST FAM20A",
                                       "Inh L1-3 VIP HSPB6", "Inh L3-5 VIP TAC3", 
                                       "Inh L1-6 LAMP5 NES")) %>% 
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
  geom_col(data9, mapping = aes(x =name,y = mean), fill = 'white', color = 'black')+
  geom_col(data9, mapping = aes(x = name, y = value, fill = name))+
  facet_wrap(~cell_type_alias_label)

ggsave(filename = "Figures/all_genes_inh_low_expressed.jpeg", dpi = 600)

##Select mutual cell type for ATXN2, ATXN7 and HTT

for (i in HTT_upregulated$cell_type){
  if (i %in% ATXN3_upregulated$cell_type){
    if (i %in% ATXN7_upregulated$cell_type){
      print(i)
    }
  }
}


data10 <- data5 %>% 
  filter(cell_type_alias_label %in% c("Inh L1 PAX6 CHRFAM7A", "Inh L1-3 SST FAM20A",
                                      "Inh L1-3 VIP HSPB6", "Inh L3-5 VIP TAC3", 
                                      "Inh L1-6 LAMP5 NES", "Inh L1 SST P4HA3",
                                      "Inh L1-5 VIP SMOC1", "Inh L2-5 VIP BSPRY")) %>% 
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
  geom_col(data10, mapping = aes(x =name,y = mean), fill = 'white', color = 'black')+
  geom_col(data10, mapping = aes(x = name, y = value, fill = name))+
  facet_wrap(~cell_type_alias_label)



ggsave(filename = "Figures/ATXN3_ATXN7_HTT_mutual_inh_low_expression.jpeg", dpi = 600)

