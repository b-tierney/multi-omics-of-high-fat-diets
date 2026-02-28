# generate intermediate files for regression etc data packet

library(tidyverse)
library(UpSetR)
library(stringr)
library(cowplot)
library(ggplot2)
library(vegan)
library(tidyr)
library(randomForest)
library(mlbench)
library(caret)
library(e1071)
library(magrittr)
library(ggnewscale)
library(ggalluvial)
library(phylobase)
library(reshape2)
library(irlba)
library(harmony)
library(ComplexUpset)
library(broom)
library(umap)
library(ecodist)
library(readxl)
library(pheatmap)
library(rlang) 
library(circlize)
library(metafor)
library(meta)
library(lme4)
library(ggpubr)
library(lmerTest)
library(broom.mixed)

theme_set(theme_cowplot())

### LOADING

### non-reversal, overall data 
setwd('~/Dropbox (Mason Lab)/mouse_diet_semir/')

# load in data, compute overall diversities, compute overall abundances per diet

data_standard = read.csv('data_packet/kraken_standard/merged_bracken_species.tsv',sep='\t') 

# filter based on number of reads and prevalence
data_standard_num = data_standard %>% select(name,all_of(grep('bracken_num',colnames(.))))%>% melt
counts = data_standard_num %>% ungroup %>% select(name,value) %>% group_by(name) %>% filter(value>50) %>% mutate(count = if_else(value>0,1,0))  %>% dplyr::summarise(s = sum(count))
tokeep = counts %>% filter(s > 10) %>% select(name) %>% unlist %>% unname

data_standard = data_standard %>% select(name,all_of(grep('bracken_frac',colnames(.)))) %>% filter(name %in% tokeep)

colnames(data_standard) = gsub('\\.bracken_frac','',colnames(data_standard))
colnames(data_standard) = gsub('\\.report','',colnames(data_standard))
data_melted = melt(data_standard)
colnames(data_melted) = c('SPECIES','UID','ABUNDANCE')
sample_data = read.csv('metadata/sample_data',sep='\t',header=T)
metadata = read.csv('metadata/metadata_20221701.csv')
metadata = inner_join(sample_data,metadata,by='FileID')  %>% mutate(animal_id = paste(Cohort,Sex,cage,Diet,animal,sep='-')) #%>% mutate(animal_id = paste(cage,animal,sep='-'))

# COHORT -- DIET -- ANIMAL -- REVERSAL

reversals4month = metadata %>% filter(grepl('M4R',timepoint) & !grepl('M9R',FileID)) %>% select(animal_id) %>% unlist %>% unname %>% unique
reversals9month = metadata %>% filter(grepl('M9R',timepoint) & !grepl('M4R',FileID)) %>% select(animal_id) %>% unlist %>% unname %>% unique
metadata = metadata %>% mutate(REVERSAL = if_else(animal_id %in% reversals4month,'REVERSAL_4MONTH','NON_REVERSAL'))
metadata = metadata %>% mutate(REVERSAL = if_else(animal_id %in% reversals9month,'REVERSAL_9MONTH',REVERSAL)) %>% filter(animal_id != 'Cohort2-M-Coconut_3-Coconut Oil-4', animal_id!="Cohort2-M-Keto_4-Ketogenic-7")
timerecode = read.csv('metadata/recode_time.csv')
timerecode$TIME_RECODE_OVERALL = timerecode$timepoint
timerecode$TIME_RECODE_OVERALL = gsub('M4R2','7',timerecode$timepoint)
timerecode$TIME_RECODE_OVERALL = gsub('M9R0','10',timerecode$TIME_RECODE_OVERALL)
timerecode$TIME_RECODE_OVERALL = gsub('M9R2','12',timerecode$TIME_RECODE_OVERALL)
timerecode$TIME_RECODE_OVERALL = gsub('M4R5','10',timerecode$TIME_RECODE_OVERALL)
timerecode$TIME_RECODE_OVERALL = gsub('M4R7','12',timerecode$TIME_RECODE_OVERALL)
timerecode$TIME_RECODE_OVERALL = gsub('D0','0',timerecode$TIME_RECODE_OVERALL)
timerecode$TIME_RECODE_OVERALL = gsub('M0','1',timerecode$TIME_RECODE_OVERALL)
timerecode$TIME_RECODE_OVERALL = gsub('M10','11',timerecode$TIME_RECODE_OVERALL)
timerecode$TIME_RECODE_OVERALL = gsub('M11','12',timerecode$TIME_RECODE_OVERALL)
timerecode$TIME_RECODE_OVERALL = gsub('M1','2',timerecode$TIME_RECODE_OVERALL)
timerecode$TIME_RECODE_OVERALL = gsub('M4','5',timerecode$TIME_RECODE_OVERALL)
timerecode$TIME_RECODE_OVERALL = gsub('M6','7',timerecode$TIME_RECODE_OVERALL)
timerecode$TIME_RECODE_OVERALL = gsub('M9','10',timerecode$TIME_RECODE_OVERALL)
timerecode$TIME_RECODE_OVERALL = as.numeric(timerecode$TIME_RECODE_OVERALL)
timerecode$TIME_ON_DIET = timerecode$timepoint
timerecode$TIME_ON_DIET = gsub('M4R2','5',timerecode$timepoint)
timerecode$TIME_ON_DIET = gsub('M9R0','10',timerecode$TIME_ON_DIET)
timerecode$TIME_ON_DIET = gsub('M9R2','10',timerecode$TIME_ON_DIET)
timerecode$TIME_ON_DIET = gsub('M4R5','5',timerecode$TIME_ON_DIET)
timerecode$TIME_ON_DIET = gsub('M4R7','5',timerecode$TIME_ON_DIET)
timerecode$TIME_ON_DIET = gsub('D0','0',timerecode$TIME_ON_DIET)
timerecode$TIME_ON_DIET = gsub('M0','1',timerecode$TIME_ON_DIET)
timerecode$TIME_ON_DIET = gsub('M10','11',timerecode$TIME_ON_DIET)
timerecode$TIME_ON_DIET = gsub('M11','12',timerecode$TIME_ON_DIET)
timerecode$TIME_ON_DIET = gsub('M1','2',timerecode$TIME_ON_DIET)
timerecode$TIME_ON_DIET = gsub('M4','5',timerecode$TIME_ON_DIET)
timerecode$TIME_ON_DIET = gsub('M6','7',timerecode$TIME_ON_DIET)
timerecode$TIME_ON_DIET = gsub('M9','10',timerecode$TIME_ON_DIET)
timerecode$TIME_ON_DIET = as.numeric(timerecode$TIME_ON_DIET)
timerecode$TIME_OFF_DIET = timerecode$timepoint
timerecode$TIME_OFF_DIET = gsub('M4R2','2',timerecode$timepoint)
timerecode$TIME_OFF_DIET = gsub('M9R0','0',timerecode$TIME_OFF_DIET)
timerecode$TIME_OFF_DIET = gsub('M9R2','2',timerecode$TIME_OFF_DIET)
timerecode$TIME_OFF_DIET = gsub('M4R5','5',timerecode$TIME_OFF_DIET)
timerecode$TIME_OFF_DIET = gsub('M4R7','7',timerecode$TIME_OFF_DIET)
timerecode$TIME_OFF_DIET = gsub('D0','0',timerecode$TIME_OFF_DIET)
timerecode$TIME_OFF_DIET = gsub('M0','0',timerecode$TIME_OFF_DIET)
timerecode$TIME_OFF_DIET = gsub('M10','0',timerecode$TIME_OFF_DIET)
timerecode$TIME_OFF_DIET = gsub('M11','0',timerecode$TIME_OFF_DIET)
timerecode$TIME_OFF_DIET = gsub('M1','0',timerecode$TIME_OFF_DIET)
timerecode$TIME_OFF_DIET = gsub('M4','0',timerecode$TIME_OFF_DIET)
timerecode$TIME_OFF_DIET = gsub('M6','0',timerecode$TIME_OFF_DIET)
timerecode$TIME_OFF_DIET = gsub('M9','0',timerecode$TIME_OFF_DIET)
timerecode$TIME_OFF_DIET = as.numeric(timerecode$TIME_OFF_DIET)
metadata = left_join(metadata,timerecode)

#metadata = metadata %>% mutate(foo=Diet)%>% pivot_wider(names_from = Diet, values_from = Diet, values_fn = function(x) 1, values_fill = 0) %>% rename(NONE="NA",Diet=foo)

#metadata = metadata %>% mutate(PREPOST_REVERSAL_4MONTH = if_else(Diet == 'Control' & REVERSAL != "NON_REVERSAL",0,as.numeric(PREPOST_REVERSAL_4MONTH))) %>% mutate(PREPOST_REVERSAL_9MONTH = if_else(Diet == 'Control' & REVERSAL != "NON_REVERSAL",0,as.numeric(PREPOST_REVERSAL_9MONTH))) 
metadata = metadata %>% mutate(KETO_OTHER_CONTROL = if_else(Diet != 'Ketogenic' & Diet !='Control','Other',Diet))

merged_data = left_join(data_melted,metadata,by=c('UID'='SampleID'))
merged_data = merged_data %>% select(UID,ABUNDANCE,SPECIES,Timepoint,Diet,Sex,Batch,Cohort,cage,animal,weight,REVERSAL,TIME_ON_DIET,TIME_OFF_DIET,KETO_OTHER_CONTROL,animal_id,all_of(colnames(timerecode)[2:ncol(timerecode)]))
colnames(merged_data) = toupper(colnames(merged_data))
merged_data = merged_data %>% filter(!is.na(SPECIES)) %>% filter(!is.na(UID))  %>% mutate(TAXUNIQUE = SPECIES)
merged_data = merged_data %>% mutate(PRE_POST_REVERSAL = if_else(DIET == 'Control','ON DIET',PRE_POST_REVERSAL))

#### ALSO MERGE XTREE DATA IN NOW

xtree = read.csv('~/Dropbox (Mason Lab)/mouse_diet_semir/xtree_merged/GTDB_.1_.05_metagenomics_species_ra.tsv',sep='\t') %>% filter(rowSums(.) != 0) %>% t %>% data.frame(check.names=F) # %>% select(-Unknown)
xtree = xtree[rowSums(xtree) %>% data.frame %>% filter(.!=0) %>% rownames,]
xtree = xtree %>% rownames_to_column('UID') %>% melt
colnames(xtree) = c('UID','SPECIES','ABUNDANCE')
merged_data_xtree = left_join(xtree,metadata,by=c('UID'='SampleID'))  %>% select(-uid,-diet,-batch,-sex,-timepoint)
colnames(merged_data_xtree)  = toupper(colnames(merged_data_xtree))

#### AND NOW DO METAPHLAN
metaph = read.table('~/Dropbox (Mason Lab)/mouse_diet_semir/data_packet/metaphlan2/metaphlan_merged.tsv',sep='\t',header=T) 
colnames(metaph) = gsub('_metaphlan','',colnames(metaph)) 
metaph = melt(metaph)
colnames(metaph) = c('SPECIES','UID','ABUNDANCE')
merged_data_metaph = left_join(metaph,metadata,by=c('UID'='SampleID'))  %>% select(-uid,-diet,-batch,-sex,-timepoint)
colnames(merged_data_metaph)  = toupper(colnames(merged_data_metaph))

### PROCESSING
process_tax_data <- function(merged_data,type){
  ### data processing -- compute alpha and beta diversities for all averaged samples
  setwd('~/Dropbox (Mason Lab)/mouse_diet_semir/')
  
  ## ALPHA
  norm_counts_wide = dcast(data = merged_data %>% dplyr::select(UID,ABUNDANCE,SPECIES), SPECIES ~ UID,value.var = 'ABUNDANCE') %>% column_to_rownames('SPECIES')
  norm_counts_wide[is.na(norm_counts_wide)] = 0
  norm_counts_wide = norm_counts_wide[colSums(norm_counts_wide)!=0] 
  annotationdata = merged_data  %>% select(-ABUNDANCE,-SPECIES) %>% distinct
  
  richness = purrr::map(colnames(norm_counts_wide), function(x) length(which(norm_counts_wide[,x]!=0))) %>% unlist %>% unname
  shannon = purrr::map(colnames(norm_counts_wide), function(x) vegan::diversity(norm_counts_wide[,x],index='shannon')) %>% unlist %>% unname
  simpson = purrr::map(colnames(norm_counts_wide), function(x) vegan::diversity(norm_counts_wide[,x],index='simpson')) %>% unlist %>% unname
  
  divs = as.data.frame(list(richness = richness,simpson = simpson,shannon = shannon))
  rownames(divs) = colnames(norm_counts_wide)  
  divs_mdat = left_join(divs %>% rownames_to_column('UID'),annotationdata)
  colnames(divs_mdat) = toupper(colnames(divs_mdat))
  
  divs_mdat = divs_mdat%>% data.frame %>% filter(!is.na(DIET),!is.na(TIMEPOINT))
  
  divs_mdat$DIET = factor(divs_mdat$DIET,levels=c('Control','Coconut Oil','Fish Oil','Ketogenic','Lard','Milkfat','Olive Oil','Palm Oil'))
  saveRDS(divs_mdat,paste('intermediate_files/',type,'_alpha_divs.rds',sep=''))
  
  ## BETA 
  
  norm_counts_wide_beta = dcast(data = merged_data %>% select(UID,ABUNDANCE,SPECIES), SPECIES ~ UID,value.var = 'ABUNDANCE') %>% column_to_rownames('SPECIES')
  norm_counts_wide_beta[is.na(norm_counts_wide_beta)] = 0
  beta_dist = as.matrix(vegdist(t(norm_counts_wide_beta),index = "bray")) %>% data.frame(check.names=F)
  beta_dist2 = beta_dist %>% rownames_to_column('UID')
  beta_dist_anno = left_join(beta_dist2,annotationdata %>% select(UID,DIET,SEX,BATCH,COHORT,CAGE,ANIMAL,WEIGHT,REVERSAL,ANIMAL_ID,TIMEPOINT) %>% distinct) 
  beta_dist_anno_melt = melt(beta_dist_anno,id.vars = c("UID","DIET", "SEX" ,"BATCH" ,"COHORT" ,"CAGE", "ANIMAL", "WEIGHT", "REVERSAL", "ANIMAL_ID",'TIMEPOINT'))
  beta_dist_anno_melt = left_join(beta_dist_anno_melt,annotationdata %>% select(UID,TIMEPOINT,ANIMAL_ID) %>% distinct %>% dplyr::rename(variable=UID,ANIMAL_ID2=ANIMAL_ID,TIMEPOINT_2=TIMEPOINT)) 
  
  saveRDS(divs_mdat,paste('intermediate_files/',type,'_beta_divs.rds',sep=''))
  
  #timepoints = unique(beta_dist_anno_sub$TIMEPOINT)
  #collist = list(DIET = c("Coconut Oil" = "red", "Control" = "black", "Ketogenic" = "blue", "Lard" = "purple", "Palm Oil" = "green", "Milkfat" = "grey","Olive Oil" = "orange", "Fish Oil" = "brown"))
  
  ### FILTERING
  
  ## only take cohort 1/cohort 2 NON reversal samples and average across all timepoints
  setwd('~/Dropbox (Mason Lab)/mouse_diet_semir/')
  
  merged_data_c12 = merged_data %>% filter(REVERSAL == 'NON_REVERSAL',COHORT =='Cohort1' | COHORT =='Cohort2')
  saveRDS(merged_data_c12,paste('intermediate_files/',type,'_merged_data_c12.rds',sep=''))
  merged_data_c12_rev = merged_data %>% filter(COHORT =='Cohort1' | COHORT =='Cohort2')
  
  saveRDS(merged_data_c12_rev,paste('intermediate_files/',type,'_all_merged_data_reversals.rds',sep=''))
  saveRDS(merged_data,paste('intermediate_files/',type,'_all_merged_data_complete.rds',sep=''))
  merged_data_c12_avg = merged_data_c12  %>% group_by(COHORT,SPECIES,DIET,ANIMAL_ID,SEX) %>% dplyr::summarise(MEAN_WEIGHT = mean(WEIGHT,na.rm=T),SD_WEIGHT = sd(WEIGHT),MEAN_ABUNDANCE = mean(ABUNDANCE,na.rm=T),SD_ABUNDANCE = sd(ABUNDANCE,na.rm=TRUE))
  
  ## do the same thing with beta and alpha diversity
  alpha_c12 = divs_mdat %>% filter(REVERSAL == 'NON_REVERSAL',COHORT =='Cohort1' | COHORT =='Cohort2')
  #alpha_c12 = alpha_c12  %>% group_by(COHORT,DIET,ANIMAL_ID,SEX) %>% summarise(MEAN_WEIGHT = mean(WEIGHT,na.rm=T),SD_WEIGHT = sd(WEIGHT),MEAN_SHANNON = mean(SHANNON,na.rm=T),SD_SHANNON = sd(SHANNON,na.rm=TRUE),MEAN_SIMPSON = mean(SIMPSON,na.rm=T),SD_SIMPSON = sd(SIMPSON,na.rm=TRUE),MEAN_RICHNESS = mean(RICHNESS,na.rm=T),SD_RICHNESS = sd(RICHNESS,na.rm=TRUE))
  
  beta_c12_intraind = beta_dist_anno_melt %>% filter(REVERSAL == 'NON_REVERSAL',COHORT =='Cohort1' | COHORT =='Cohort2') %>% mutate(AID2 = ANIMAL_ID)%>% filter(ANIMAL_ID == ANIMAL_ID2, UID !=variable) #%>% select(-variable,-UID) %>% group_by(DIET,SEX,COHORT,CAGE,ANIMAL,REVERSAL,ANIMAL_ID) 
  
  # save tmp files so it can be loaded later if desired 
  saveRDS(merged_data_c12_avg,paste('intermediate_files/',type,'_average_abundances_cohort12.rds'))
  saveRDS(alpha_c12,paste('intermediate_files/',type,'_alpha_avg_cohort12.rds'))
  saveRDS(beta_c12_intraind,paste('intermediate_files/',type,'_betadiv_med_cohort12.rds'))
  
  saveRDS(merged_data_c12_avg,paste('intermediate_files/',type,'_average_abundances_cohort12.rds',sep=''))
  saveRDS(alpha_c12,paste('intermediate_files/',type,'_alpha_avg_cohort12.rds',sep=''))
  saveRDS(beta_c12_intraind,paste('intermediate_files/',type,'_betadiv_med_cohort12.rds',sep=''))
  
  ### ANALYSIS -- AVERAGE/MEDIAN DIVERSITY
  
  # show average alpha diversity by diet by cohort
  ggplot(data = alpha_c12 %>% select(UID,SIMPSON,SHANNON,RICHNESS,DIET,COHORT)%>% distinct %>% melt, aes(x = DIET, y = value,fill=COHORT)) +geom_boxplot() + facet_wrap(.~variable,scales='free') +theme(axis.text.x = element_text(angle = 60,hjust=1))
  
  # median beta diversity by diet
  beta_c12_intraind$UID = as.character(beta_c12_intraind$UID)
  beta_c12_intraind_arr = beta_c12_intraind %>% select(ANIMAL_ID,UID,value,DIET,COHORT)  %>% arrange(COHORT,DIET,value) %>% mutate(val = seq(1,nrow(beta_c12_intraind)))
  
  beta_c12_intraind_arr$ANIMAL_ID = fct_reorder(beta_c12_intraind_arr$ANIMAL_ID,beta_c12_intraind_arr$val)
  
  ggplot(data =beta_c12_intraind_arr, aes(x = DIET, y = value, group= ANIMAL_ID)) +geom_boxplot() + xlab('') + facet_wrap(COHORT~.,nrow=2)+theme(legend.position="none")
  
  # median abundance by diet
  specieslist=unique(merged_data_c12_avg$SPECIES)
  
  merged_data_c12_avg$DIET = factor(merged_data_c12_avg$DIET,levels=c('Control','Coconut Oil','Fish Oil','Ketogenic','Lard','Milkfat','Olive Oil','Palm Oil'))
  
  saveRDS(merged_data_c12_avg,paste('intermediate_files/',type,'_merged_data_c12_avg.rds',sep=''))
  
  merged_data_c1_avg_wide = merged_data_c12_avg %>% filter(COHORT=='Cohort1') %>% ungroup %>% select(ANIMAL_ID,SPECIES,DIET,SEX,MEAN_ABUNDANCE,MEAN_WEIGHT) %>% dcast(ANIMAL_ID +MEAN_WEIGHT+ DIET + SEX ~ SPECIES, value.var = 'MEAN_ABUNDANCE') %>% select(-ANIMAL_ID)
  merged_data_c1_avg_wide[is.na(merged_data_c1_avg_wide)]=0
  saveRDS(merged_data_c1_avg_wide,paste('intermediate_files/',type,'_merged_data_c1_avg_wide,rds',sep=''))
}

process_tax_data(merged_data,'kraken2')
process_tax_data(merged_data_xtree,'xtree')
process_tax_data(merged_data_metaph,'metaphlan4')

# humann3 pathway output
h3p = read.csv('data_packet/processed_humann_abundance_tables/humann_pathabundance.tsv',sep='\t')
colnames(h3p)[1]='PATHWAY'
data_melted = melt(h3p)
data_melted$variable = gsub('_all_reads_Abundance','',data_melted$variable)
colnames(data_melted) = c('PATHWAY','UID','ABUNDANCE')
sample_data = read.csv('metadata/sample_data',sep='\t',header=T)
metadata = read.csv('metadata/metadata_20221701.csv')
metadata = inner_join(sample_data,metadata,by='FileID') %>% mutate(animal_id = paste(cage,animal,sep='-'))
reversals4month = metadata %>% filter(grepl('M4R',timepoint)) %>% select(animal_id) %>% unlist %>% unname %>% unique
reversals9month = metadata %>% filter(grepl('M9R',timepoint)) %>% select(animal_id) %>% unlist %>% unname %>% unique
metadata = metadata %>% mutate(REVERSAL = if_else(animal_id %in% reversals4month,'REVERSAL_4MONTH','NON_REVERSAL'))
metadata = metadata %>% mutate(REVERSAL = if_else(animal_id %in% reversals9month,'REVERSAL_9MONTH',REVERSAL))
timerecode = read.csv('metadata/recode_time.csv')
timerecode$TIME_RECODE_OVERALL = timerecode$timepoint
timerecode$TIME_RECODE_OVERALL = gsub('M4R2','7',timerecode$timepoint)
timerecode$TIME_RECODE_OVERALL = gsub('M9R0','10',timerecode$TIME_RECODE_OVERALL)
timerecode$TIME_RECODE_OVERALL = gsub('M9R2','12',timerecode$TIME_RECODE_OVERALL)
timerecode$TIME_RECODE_OVERALL = gsub('M4R5','10',timerecode$TIME_RECODE_OVERALL)
timerecode$TIME_RECODE_OVERALL = gsub('M4R7','12',timerecode$TIME_RECODE_OVERALL)
timerecode$TIME_RECODE_OVERALL = gsub('D0','0',timerecode$TIME_RECODE_OVERALL)
timerecode$TIME_RECODE_OVERALL = gsub('M0','1',timerecode$TIME_RECODE_OVERALL)
timerecode$TIME_RECODE_OVERALL = gsub('M10','11',timerecode$TIME_RECODE_OVERALL)
timerecode$TIME_RECODE_OVERALL = gsub('M11','12',timerecode$TIME_RECODE_OVERALL)
timerecode$TIME_RECODE_OVERALL = gsub('M1','2',timerecode$TIME_RECODE_OVERALL)
timerecode$TIME_RECODE_OVERALL = gsub('M4','5',timerecode$TIME_RECODE_OVERALL)
timerecode$TIME_RECODE_OVERALL = gsub('M6','7',timerecode$TIME_RECODE_OVERALL)
timerecode$TIME_RECODE_OVERALL = gsub('M9','10',timerecode$TIME_RECODE_OVERALL)
timerecode$TIME_RECODE_OVERALL = as.numeric(timerecode$TIME_RECODE_OVERALL)
timerecode$TIME_ON_DIET = timerecode$timepoint
timerecode$TIME_ON_DIET = gsub('M4R2','5',timerecode$timepoint)
timerecode$TIME_ON_DIET = gsub('M9R0','10',timerecode$TIME_ON_DIET)
timerecode$TIME_ON_DIET = gsub('M9R2','10',timerecode$TIME_ON_DIET)
timerecode$TIME_ON_DIET = gsub('M4R5','5',timerecode$TIME_ON_DIET)
timerecode$TIME_ON_DIET = gsub('M4R7','5',timerecode$TIME_ON_DIET)
timerecode$TIME_ON_DIET = gsub('D0','0',timerecode$TIME_ON_DIET)
timerecode$TIME_ON_DIET = gsub('M0','1',timerecode$TIME_ON_DIET)
timerecode$TIME_ON_DIET = gsub('M10','11',timerecode$TIME_ON_DIET)
timerecode$TIME_ON_DIET = gsub('M11','12',timerecode$TIME_ON_DIET)
timerecode$TIME_ON_DIET = gsub('M1','2',timerecode$TIME_ON_DIET)
timerecode$TIME_ON_DIET = gsub('M4','5',timerecode$TIME_ON_DIET)
timerecode$TIME_ON_DIET = gsub('M6','7',timerecode$TIME_ON_DIET)
timerecode$TIME_ON_DIET = gsub('M9','10',timerecode$TIME_ON_DIET)
timerecode$TIME_ON_DIET = as.numeric(timerecode$TIME_ON_DIET)
timerecode$TIME_OFF_DIET = timerecode$timepoint
timerecode$TIME_OFF_DIET = gsub('M4R2','2',timerecode$timepoint)
timerecode$TIME_OFF_DIET = gsub('M9R0','0',timerecode$TIME_OFF_DIET)
timerecode$TIME_OFF_DIET = gsub('M9R2','2',timerecode$TIME_OFF_DIET)
timerecode$TIME_OFF_DIET = gsub('M4R5','5',timerecode$TIME_OFF_DIET)
timerecode$TIME_OFF_DIET = gsub('M4R7','7',timerecode$TIME_OFF_DIET)
timerecode$TIME_OFF_DIET = gsub('D0','0',timerecode$TIME_OFF_DIET)
timerecode$TIME_OFF_DIET = gsub('M0','0',timerecode$TIME_OFF_DIET)
timerecode$TIME_OFF_DIET = gsub('M10','0',timerecode$TIME_OFF_DIET)
timerecode$TIME_OFF_DIET = gsub('M11','0',timerecode$TIME_OFF_DIET)
timerecode$TIME_OFF_DIET = gsub('M1','0',timerecode$TIME_OFF_DIET)
timerecode$TIME_OFF_DIET = gsub('M4','0',timerecode$TIME_OFF_DIET)
timerecode$TIME_OFF_DIET = gsub('M6','0',timerecode$TIME_OFF_DIET)
timerecode$TIME_OFF_DIET = gsub('M9','0',timerecode$TIME_OFF_DIET)
timerecode$TIME_OFF_DIET = as.numeric(timerecode$TIME_OFF_DIET)
metadata = left_join(metadata,timerecode)
metadata = metadata %>% mutate(KETO_OTHER_CONTROL = if_else(Diet != 'Ketogenic' & Diet !='Control','Other',Diet))
metadata = left_join(metadata,timerecode)
merged_data = left_join(data_melted,metadata,by=c('UID'='SampleID'))
merged_data = merged_data %>% select(UID,ABUNDANCE,PATHWAY,Timepoint,Diet,Sex,Batch,Cohort,cage,animal,weight,KETO_OTHER_CONTROL,REVERSAL,animal_id,all_of(colnames(timerecode)[2:ncol(timerecode)]))
colnames(merged_data) = toupper(colnames(merged_data))
merged_data = merged_data %>% filter(!is.na(PATHWAY)) %>% filter(!is.na(UID)) 


process_tax_data(merged_data %>% dplyr::rename(SPECIES=PATHWAY),'pathway')

#saveRDS(merged_data,'intermediate_files/humann3_path_abundances_merged.rds')

