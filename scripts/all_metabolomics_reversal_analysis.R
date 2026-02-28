

library(tidyverse)
library(furrr)
library(broom)
library(ComplexUpset)
library(circlize)
library(progressr)
library(glmnet)
library(ggrepel)
library(reshape2)
library(caret)
library(tidygraph)
library(RColorBrewer)
library(igraph)
library(ComplexHeatmap)
library(ggplot2)

# Clean metabolite names - remove [POS]>, [NEG]>, and _C0 suffixes
clean_metabolite_name <- function(name) {
  name %>%
    gsub("^\\[POS\\]>", "", .) %>%
    gsub("^\\[NEG\\]>", "", .) %>%
    gsub("_C0$", "", .) %>%
    trimws()
}


regress <- function(s,data){
  print(s)
  data_sub = data %>% filter(Metabolite == s)
  reg1 = glm(data = data_sub , ValueLog10 ~ SEX + Diet) 
  reg1 = broom::tidy(reg1) %>% mutate(variable = s)
  return(reg1)
}


# filter out the pathways??
get_metabolite_reps <- function(data, correlation_threshold) {
  correlation_matrix <- cor(t(data))
  hc <- hclust(as.dist(1 - correlation_matrix), method = "complete")
  groups <- cutree(hc, h = 1 - correlation_threshold)
  unique_groups <- unique(groups)
  representative_metabolites <- vector("list", length(unique_groups))
  for (i in seq_along(unique_groups)) {
    group_indices <- which(groups == unique_groups[i])
    cluster_metabolites <- rownames(data)[group_indices]
    filtered_metabolites <- cluster_metabolites[!grepl("^X-", cluster_metabolites)]
    if (length(filtered_metabolites) == 0 && length(cluster_metabolites) == 1) {
      filtered_metabolites <- cluster_metabolites
    }
    if (length(filtered_metabolites) > 0) {
      max_prevalence_metabolite <- (rowSums(data[group_indices, ] != 0)) %>% data.frame %>% rownames_to_column('name') %>% filter(!grepl("^X-", name)) %>% arrange(desc(.)) %>% head(1)
      max_prevalence_metabolite = max_prevalence_metabolite[[1]]
      representative_metabolites[[i]] <- max_prevalence_metabolite
    }
    if (length(filtered_metabolites) == 0 & length(cluster_metabolites)>0) {
      max_prevalence_metabolite <- cluster_metabolites[which.max(rowSums(data[group_indices, ] != 0))]
      representative_metabolites[[i]] <- max_prevalence_metabolite
    }
  }
  representative_metabolites = unlist(unname(representative_metabolites))
  return(representative_metabolites)
}


colorsDiet = c(
  Control = "#377eb8",
  Lard = "#e41a1c",
  Milkfat = "#800026",
  MilkFat = "#800026",
  Ketogenic = "#ff7f00",
  `Fish Oil` = "#f781bf",
  `Coconut Oil` = "#a65628",
  `Olive Oil` = "#984ea3",
  `Palm Oil` = "#4daf4a",
  `FishOil` = "#f781bf",
  `CoconutOil` = "#a65628",
  `OliveOil` = "#984ea3",
  `PalmOil` = "#4daf4a"
)
setwd('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/')

##### get sample mdat
merged_data = readRDS('intermediate_files/metaphlan4_all_merged_data_complete.rds') %>% filter(grepl('t__',SPECIES)) %>% mutate(SGB = strsplit(SPECIES,'t__') %>% map_chr(2))
merged_data = merged_data %>% mutate(ANIMAL_ID = paste(CAGE,ANIMAL,sep='-'))
merged_data = merged_data %>% filter(SGB %in% (merged_data %>% dcast(SGB ~ UID,value.var = 'ABUNDANCE') %>% column_to_rownames('SGB') %>% rowSums %>% data.frame %>% filter(.!=0) %>% rownames))
merged_data$DIET = factor(merged_data$DIET,levels=c('Control','Coconut Oil','Fish Oil','Ketogenic','Lard','Milkfat','Olive Oil','Palm Oil'))
merged_data$ANIMAL_ID = as.factor(merged_data$ANIMAL_ID)
merged_data_tpm = merged_data %>% filter(COHORT == 'Cohort2')  #%>% mutate(SampleID =  strsplit(' ',ANIMAL_ID) %>% map_chr(2))
merged_data_tpm = merged_data_tpm %>% select(ANIMAL_ID,ABUNDANCE,WEIGHT,SEX ,TIME_RECODE_OVERALL ,DIET,ANIMAL_ID,SGB,COHORT)
merged_data_tpm_avg = merged_data_tpm  %>% dplyr::group_by(ANIMAL_ID,SGB,SEX,DIET) %>% summarise(WEIGHT = mean(WEIGHT,na.rm=T),MEAN_ABUNDANCE = mean(ABUNDANCE,na.rm=T),SD_ABUNDANCE = sd(ABUNDANCE,na.rm=TRUE))
merged_data_tpm_avg$SampleID = gsub('-','_',merged_data_tpm_avg$ANIMAL_ID)
samplemdat = merged_data_tpm_avg %>% ungroup %>% select(SampleID,SEX,WEIGHT,DIET) %>% unique %>% distinct

##### SERUM

metabolites = read.csv('data_packet/metabolites/PlasmaMBX_2023_B2_NormLog_long.txt',sep='\t') 
metabolites$SampleID = gsub('.rev','',metabolites$SampleID )
metabolites$SampleID = gsub('_F_M6','',metabolites$SampleID )
metabolites$SampleID = gsub('_F_M11','',metabolites$SampleID )
metabolites$SampleID = gsub('_M_M11','',metabolites$SampleID )
metabolites$SampleID = gsub('_M_M6','',metabolites$SampleID )
metabolites$SampleID = gsub('_M8_2_22_2020','',metabolites$SampleID )
metabolites$SampleID = gsub('Coco_','Coconut_',metabolites$SampleID )
metabolites$SampleID = gsub('Fish.Oil','Fish Oil',metabolites$SampleID)
metabolites2 = metabolites %>% filter(Time == 'M6reversal' | Time == 'M11reversal') %>% select(-Time) %>% group_by(Metabolite,SampleID,Diet,Sex) %>% summarise(Value = mean(Value)) %>% ungroup %>% dplyr::rename(ValueLog10 = Value) %>% mutate(rawval = 10^(ValueLog10))  %>% select(Metabolite,SampleID,ValueLog10,Sex,Diet,rawval) %>% dplyr::rename(SEX = Sex)

metabolites2$Diet = factor(metabolites2$Diet,levels=c('Control','CoconutOil','FishOil','Ketogenic','Lard','MilkFat','OliveOil','PalmOil'))

### filter metabolites based on prevalence -- acting under assumption minimum value is the zero val
tokeep = metabolites2 %>% group_by(Metabolite) %>% summarise(minval = min(rawval),count = sum(rawval>minval)) %>% arrange(count) %>% filter(count>10) %>% select(Metabolite) %>% unlist %>% unname

metabolites2 = metabolites2 %>% filter(Metabolite %in% tokeep)

tmp = metabolites2 %>% dcast(Metabolite ~ SampleID,value.var ='ValueLog10') %>% column_to_rownames('Metabolite')
repmets =  get_metabolite_reps(tmp,.75)


date()
output = purrr::map(unique(metabolites2$Metabolite), function(x) regress(x,metabolites2))
output = bind_rows(output) %>% filter(term != '(Intercept)') %>% mutate(BH = p.adjust(p.value,method="BH"))
output$variable = clean_metabolite_name(output$variable)
date()
print('Finished linear modeling...')
write.csv(output,'association_output/metabolic_pathway_regression_output_dietSERUMMET_REVERSALS.csv')

##### FECAL

metabolites = read.csv('data_packet/metabolites/FecalMetabolite_logTransformed_long.txt',sep='\t') %>% filter(Time == 'M4R3') %>% select(-Time) %>% mutate(rawval = 10^(ValueLog10)) %>% mutate()
metabolites$SampleID = gsub('_M4R3_1_24_2020','',metabolites$SampleID )
metabolites$SampleID = gsub('Fish.Oil','Fish Oil',metabolites$SampleID)
metabolites$Diet = factor(metabolites$Diet,levels=c('Control','CoconutOil','FishOil','Ketogenic','Lard','MilkFat','OliveOil','PalmOil'))

### filter metabolites based on prevalence -- acting under assumption minimum value is the zero val
tokeep = metabolites %>% group_by(Metabolite) %>% summarise(minval = min(rawval),count = sum(rawval>minval)) %>% arrange(count) %>% filter(count>10) %>% select(Metabolite) %>% unlist %>% unname

metabolites = metabolites %>% filter(Metabolite %in% tokeep)
metabolites2 = inner_join(metabolites,samplemdat%>% ungroup %>% select(SEX,SampleID) %>% distinct)

# filter out the pathways??
get_metabolite_reps <- function(data, correlation_threshold) {
  correlation_matrix <- cor(t(data))
  hc <- hclust(as.dist(1 - correlation_matrix), method = "complete")
  groups <- cutree(hc, h = 1 - correlation_threshold)
  unique_groups <- unique(groups)
  representative_metabolites <- vector("list", length(unique_groups))
  for (i in seq_along(unique_groups)) {
    group_indices <- which(groups == unique_groups[i])
    cluster_metabolites <- rownames(data)[group_indices]
    filtered_metabolites <- cluster_metabolites[!grepl("^X-", cluster_metabolites)]
    if (length(filtered_metabolites) == 0 && length(cluster_metabolites) == 1) {
      filtered_metabolites <- cluster_metabolites
    }
    if (length(filtered_metabolites) > 0) {
      max_prevalence_metabolite <- (rowSums(data[group_indices, ] != 0)) %>% data.frame %>% rownames_to_column('name') %>% filter(!grepl("^X-", name)) %>% arrange(desc(.)) %>% head(1)
      max_prevalence_metabolite = max_prevalence_metabolite[[1]]
      representative_metabolites[[i]] <- max_prevalence_metabolite
    }
    if (length(filtered_metabolites) == 0 & length(cluster_metabolites)>0) {
      max_prevalence_metabolite <- cluster_metabolites[which.max(rowSums(data[group_indices, ] != 0))]
      representative_metabolites[[i]] <- max_prevalence_metabolite
    }
  }
  representative_metabolites = unlist(unname(representative_metabolites))
  return(representative_metabolites)
}

tmp = metabolites2 %>% dcast(Metabolite ~ SampleID,value.var ='ValueLog10') %>% column_to_rownames('Metabolite')
repmets =  get_metabolite_reps(tmp,.75)


date()
output = purrr::map(unique(metabolites2$Metabolite), function(x) regress(x,metabolites2))
output = bind_rows(output) %>% filter(term != '(Intercept)') %>% mutate(BH = p.adjust(p.value,method="BH"))
output$variable = clean_metabolite_name(output$variable)
date()
print('Finished linear modeling...')
write.csv(output,'association_output/metabolic_pathway_regression_output_diet_REVERSALS.csv')


#### LIPID

metabolites = read.csv('data_packet/metabolites/Diet_lipidomics_batch2.txt',sep='\t') %>% filter(Time == 'M4R4') %>% select(-Time) %>% mutate(rawval = 10^(ValueLog10))  %>% dplyr::rename(Metabolite = Compound) %>% select(Metabolite,SampleID,ValueLog10,Diet,Type,rawval) %>% mutate(SampleID = gsub('\\+','_',strsplit(SampleID,' ') %>% map_chr(1))) %>% group_by(Metabolite,SampleID,Diet,Type) %>% summarise(ValueLog10 = mean(ValueLog10),rawval = mean(rawval))
metabolites$SampleID = gsub('Fish.Oil','Fish Oil',metabolites$SampleID)

metabolites2 = inner_join(metabolites,samplemdat %>% ungroup %>% select(SampleID,SEX) %>% distinct)
metabolites2$Diet = factor(metabolites2$Diet,levels=c('Control','CoconutOil','FishOil','Ketogenic','Lard','MilkFat','OliveOil','PalmOil'))

### filter metabolites based on prevalence -- acting under assumption minimum value is the zero val
tokeep = metabolites2 %>% group_by(Metabolite) %>% summarise(minval = min(rawval),count = sum(rawval>minval)) %>% arrange(count) %>% filter(count>10) %>% select(Metabolite) %>% unlist %>% unname

metabolites2 = metabolites2 %>% filter(Metabolite %in% tokeep)

# filter out the pathways??
get_metabolite_reps <- function(data, correlation_threshold) {
  correlation_matrix <- cor(t(data))
  hc <- hclust(as.dist(1 - correlation_matrix), method = "complete")
  groups <- cutree(hc, h = 1 - correlation_threshold)
  unique_groups <- unique(groups)
  representative_metabolites <- vector("list", length(unique_groups))
  for (i in seq_along(unique_groups)) {
    group_indices <- which(groups == unique_groups[i])
    cluster_metabolites <- rownames(data)[group_indices]
    filtered_metabolites <- cluster_metabolites[!grepl("^X-", cluster_metabolites)]
    if (length(filtered_metabolites) == 0 && length(cluster_metabolites) == 1) {
      filtered_metabolites <- cluster_metabolites
    }
    if (length(filtered_metabolites) > 0) {
      max_prevalence_metabolite <- (rowSums(data[group_indices, ] != 0)) %>% data.frame %>% rownames_to_column('name') %>% filter(!grepl("^X-", name)) %>% arrange(desc(.)) %>% head(1)
      max_prevalence_metabolite = max_prevalence_metabolite[[1]]
      representative_metabolites[[i]] <- max_prevalence_metabolite
    }
    if (length(filtered_metabolites) == 0 & length(cluster_metabolites)>0) {
      max_prevalence_metabolite <- cluster_metabolites[which.max(rowSums(data[group_indices, ] != 0))]
      representative_metabolites[[i]] <- max_prevalence_metabolite
    }
  }
  representative_metabolites = unlist(unname(representative_metabolites))
  return(representative_metabolites)
}

tmp = metabolites2 %>% dcast(Metabolite ~ SampleID,value.var ='ValueLog10') %>% column_to_rownames('Metabolite')
repmets =  get_metabolite_reps(tmp,.75)

date()
output = purrr::map(unique(metabolites2$Metabolite), function(x) regress(x,metabolites2))
output = bind_rows(output) %>% filter(term != '(Intercept)') %>% mutate(BH = p.adjust(p.value,method="BH"))
output$variable = clean_metabolite_name(output$variable)
date()
print('Finished linear modeling...')
write.csv(output,'association_output/metabolic_pathway_regression_output_dietLIPIDS_REVERSALS.csv')

####### INTEGRATE AND COMPARE REVERSAL NON REVERSAL REGRESSION OUTPUT

serum_nr = read.csv('association_output/metabolic_pathway_regression_output_dietSERUMMET.csv') %>% filter(BH<0.05) %>% select(variable,term,estimate) %>% dplyr::rename(estimate_nr = estimate)
fecal_nr = read.csv('association_output/metabolic_pathway_regression_output_diet.csv')  %>% filter(BH<0.05) %>% select(variable,term,estimate) %>% dplyr::rename(estimate_nr = estimate)
lipids_nr = read.csv('association_output/metabolic_pathway_regression_output_dietLIPIDS.csv')  %>% filter(BH<0.05) %>% select(variable,term,estimate) %>% dplyr::rename(estimate_nr = estimate)

serum_r = read.csv('association_output/metabolic_pathway_regression_output_dietSERUMMET_REVERSALS.csv') %>% filter(BH<0.05) %>% select(variable,term,estimate) %>% dplyr::rename(estimate_r = estimate)
fecal_r = read.csv('association_output/metabolic_pathway_regression_output_diet_REVERSALS.csv')  %>% filter(BH<0.05) %>% select(variable,term,estimate) %>% dplyr::rename(estimate_r = estimate)
lipids_r = read.csv('association_output/metabolic_pathway_regression_output_dietLIPIDS_REVERSALS.csv')  %>% filter(BH<0.05) %>% select(variable,term,estimate) %>% dplyr::rename(estimate_r = estimate)

serum = left_join(serum_nr,serum_r) %>% mutate(r = if_else(estimate_nr>0 & estimate_r > 0 ,'Persistently Increased',NA))  %>% mutate(r = if_else( estimate_nr<0 & estimate_r< 0 ,'Persistently Decreased',r)) %>% mutate(dir = if_else(estimate_nr>0 ,'Positive',NA)) %>% mutate(dir = if_else(estimate_nr<0 ,'Negative',dir))%>% mutate(type = 'Serum Metabolites')
fecal = left_join(fecal_nr,fecal_r) %>% mutate(r = if_else(estimate_nr>0 & estimate_r > 0 ,'Persistently Increased',NA))  %>% mutate(r = if_else( estimate_nr<0 & estimate_r< 0 ,'Persistently Decreased',r)) %>% mutate(dir = if_else(estimate_nr>0 ,'Positive',NA)) %>% mutate(dir = if_else(estimate_nr<0 ,'Negative',dir))%>% mutate(type = 'Fecal Metabolites')
lipids = left_join(lipids_nr,lipids_r) %>% mutate(r = if_else(estimate_nr>0 & estimate_r > 0 ,'Persistently Increased',NA))  %>% mutate(r = if_else( estimate_nr<0 & estimate_r< 0 ,'Persistently Decreased',r)) %>% mutate(dir = if_else(estimate_nr>0 ,'Positive',NA)) %>% mutate(dir = if_else(estimate_nr<0 ,'Negative',dir)) %>% mutate(type = 'Lipids')

all = bind_rows(serum,fecal,lipids)
all$term = gsub('Diet','',all$term)
all$term = factor(all$term,levels = c('SEXM','Ketogenic','FishOil','PalmOil','OliveOil','CoconutOil','MilkFat','Lard'))
all_nr = all %>% group_by(type,term,dir) %>% count()
all_r = all %>% group_by(type,term,r) %>% count() %>% filter(!is.na(r))


##### BARPLOTS

ggplot(data = all_nr,aes(x=term,y=n,fill=dir)) +cowplot::theme_cowplot()+ geom_bar(stat='identity') + scale_fill_manual(values = c('darkblue','#BE1E2D')) + ylab('# associated features') + theme(legend.title = element_blank(),axis.title.y = element_blank(),legend.position = 'bottom')+ facet_grid(. ~ type) + theme(axis.text.x = element_text(angle = 60,hjust =1 ),axis.title.x =  element_blank())
ggsave('./plots/metabolite_association_overview.pdf',width=12,height=3)


ggplot(data = all_r,aes(x=term,y=n,fill=r)) +cowplot::theme_cowplot()+ geom_bar(stat='identity') + scale_fill_manual(values = c('darkblue','#BE1E2D')) + ylab('# associated features') + theme(legend.title = element_blank(),axis.title.y = element_blank(),legend.position = 'bottom')+ facet_grid(. ~ type,scales= 'free_x') + theme(axis.text.x = element_text(angle = 60,hjust =1 ),axis.title.x =  element_blank())
ggsave('./plots/metabolite_reversal_persistence.pdf',width=6,height=4)

##### FDR-based persistent metabolites dotplot (no functions)
allsig_fdr = all %>% filter(!is.na(r),term!='SEXM')
ggplot(allsig_fdr, aes(x=fct_reorder(variable,estimate_r), y=estimate_r,
                      color=term, shape=type, size=2)) +
  cowplot::theme_cowplot() + geom_point(alpha=0.85) +
  xlab('Metabolite') + ylab('Beta (reversal model)') +
  ggtitle('Persistently altered metabolites post-reversal') + geom_hline(linetype = 'dashed',yintercept = 0) + 
  theme(legend.position='bottom', legend.title=element_blank(),
        axis.text.x = element_text(angle=60, hjust=1))
ggsave('./plots/metabolite_reversal_dotplot_persistent_fdr.pdf',width=12,height=6)

##### Recompute persistence using nominal p-values (not BH), redo barplot + dotplot
serum_nr_full = read.csv('association_output/metabolic_pathway_regression_output_dietSERUMMET.csv')
serum_r_full  = read.csv('association_output/metabolic_pathway_regression_output_dietSERUMMET_REVERSALS.csv')
fecal_nr_full = read.csv('association_output/metabolic_pathway_regression_output_diet.csv')
fecal_r_full  = read.csv('association_output/metabolic_pathway_regression_output_diet_REVERSALS.csv')
lipids_nr_full = read.csv('association_output/metabolic_pathway_regression_output_dietLIPIDS.csv')
lipids_r_full  = read.csv('association_output/metabolic_pathway_regression_output_dietLIPIDS_REVERSALS.csv')

serum_nr_p = serum_nr_full %>% filter(p.value<0.05,term!='SEXM') %>% dplyr::rename(estimate_nr=estimate,p_nr=p.value)
serum_r_p  = serum_r_full  %>% filter(p.value<0.05,term!='SEXM') %>% dplyr::rename(estimate_r=estimate,p_r=p.value)
fecal_nr_p = fecal_nr_full %>% filter(p.value<0.05,term!='SEXM') %>% dplyr::rename(estimate_nr=estimate,p_nr=p.value)
fecal_r_p  = fecal_r_full  %>% filter(p.value<0.05,term!='SEXM') %>% dplyr::rename(estimate_r=estimate,p_r=p.value)
lipids_nr_p = lipids_nr_full %>% filter(p.value<0.05,term!='SEXM') %>% dplyr::rename(estimate_nr=estimate,p_nr=p.value)
lipids_r_p  = lipids_r_full  %>% filter(p.value<0.05,term!='SEXM') %>% dplyr::rename(estimate_r=estimate,p_r=p.value)

serum_p = left_join(serum_nr_p,serum_r_p,by=c('variable','term')) %>% mutate(r_p = if_else(estimate_nr>0 & estimate_r>0,'Persistently Increased', if_else(estimate_nr<0 & estimate_r<0,'Persistently Decreased',NA_character_)), type='Serum Metabolites')
fecal_p = left_join(fecal_nr_p,fecal_r_p,by=c('variable','term')) %>% mutate(r_p = if_else(estimate_nr>0 & estimate_r>0,'Persistently Increased', if_else(estimate_nr<0 & estimate_r<0,'Persistently Decreased',NA_character_)), type='Fecal Metabolites')
lipids_p = left_join(lipids_nr_p,lipids_r_p,by=c('variable','term')) %>% mutate(r_p = if_else(estimate_nr>0 & estimate_r>0,'Persistently Increased', if_else(estimate_nr<0 & estimate_r<0,'Persistently Decreased',NA_character_)), type='Lipids')

all_p = bind_rows(serum_p,fecal_p,lipids_p) %>% mutate(term=gsub('Diet','',term))

all_r_p = all_p %>% filter(!is.na(r_p)) %>% group_by(type,term,r_p) %>% count()
ggplot(all_r_p,aes(x=term,y=n,fill=r_p)) +
  cowplot::theme_cowplot() + geom_bar(stat='identity') +
  facet_grid(.~type, scales='free_x') +
  ylab('# persistent (p<0.05)') + xlab('Diet') +
  scale_fill_manual(values=c('Persistently Increased'='#BE1E2D','Persistently Decreased'='darkblue')) +
  theme(axis.text.x = element_text(angle=60,hjust=1), legend.title=element_blank(), legend.position='bottom')
ggsave('./plots/metabolite_reversal_persistence_pvalue_barplot.pdf',width=16,height=6)

allsig_p = all_p %>% filter(!is.na(r_p))
ggplot(allsig_p, aes(x=fct_reorder(variable,estimate_r), y=estimate_r,
                     color=term, shape=type, size=2)) +
  cowplot::theme_cowplot() + geom_point(alpha=0.85) +
  geom_hline(linetype = 'dashed', yintercept = 0) +
  xlab('Metabolite') + ylab('Beta (reversal model)') +
  ggtitle('Persistently altered metabolites post-reversal (nominal significance)') +
  theme(legend.position='bottom', legend.title=element_blank(),
        axis.text.x = element_text(angle=60, hjust=1))
ggsave('./plots/metabolite_reversal_dotplot_persistent_pvalue.pdf',width=16,height=6)

##### Microbe–metabolite links (lasso coefficients) for persistent metabolites
# load lasso coefficients (fecal, serum, lipid) and map SGB to GTDB species
gtdbmap = read.delim('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/mpa_vJan21_CHOCOPhlAnSGB_202103_SGB2GTDB.tsv', sep='\t', header=FALSE)
colnames(gtdbmap) = c('SGB','GTDB_TAX')

lasso_files = c(
  fecal = 'association_output/metabolic_metaphlan4_regression_output_lasso.csv',
  serum = 'association_output/metabolic_metaphlan4_regression_output_lassoSERUMMET.csv',
  lipid = 'association_output/metabolic_metaphlan4_regression_output_lassoLIPIDS.csv'
)

read_lasso_long = function(path, source_label){
  df = read.csv(path, check.names = FALSE)
  colnames(df)[1] = 'SGB'
  reshape2::melt(df, id.vars = 'SGB', variable.name = 'Metabolite', value.name = 'Coefficient') %>%
    mutate(lasso_source = source_label)
}

lasso_long = purrr::map2(lasso_files, names(lasso_files), read_lasso_long) %>%
  bind_rows() %>%
  left_join(gtdbmap, by='SGB') %>%
  separate(GTDB_TAX, into=c("Domain","Phylum","Class","Order","Family","Genus","Species"), sep=";", remove=TRUE, extra="merge")

# pairwise associations for significance markers
pairwise_sig = read.csv('association_output/metabolic_metaphlan4_regression_output_pairwise.csv') %>%
  left_join(gtdbmap, by='SGB') %>%
  separate(GTDB_TAX, into=c("Domain","Phylum","Class","Order","Family","Genus","Species"), sep=";", remove=TRUE, extra="merge") %>%
  mutate(BH = p.adjust(p.value, method='BH')) %>%
  filter(!is.na(Species), !grepl(' sp', Species))

speciesdata = lasso_long %>%
  filter(Metabolite %in% allsig_fdr$variable) %>%
  mutate(direction = if_else(Coefficient > 0, "Positive", "Negative"),
         size_val = abs(Coefficient)) %>%
  filter(size_val > 0.1) %>%  # drop near-zero lasso effects
  filter(!is.na(Species), !grepl(' sp', Species))

##### Gridded microbe–metabolite dotplot with fixed metabolite axis (FDR set)
met_order_fdr = allsig_fdr %>% pull(variable) %>% unique()
species_order_fdr = speciesdata %>% arrange(Species) %>% pull(Species) %>% unique()

grid_dot_fdr = speciesdata %>%
  mutate(Metabolite = factor(Metabolite, levels = met_order_fdr),
         Species = factor(Species, levels = species_order_fdr))

p_sig_fdr = pairwise_sig %>%
  filter(Metabolite %in% allsig_fdr$variable,
         Species %in% species_order_fdr) %>%
  mutate(sig_lab = case_when(BH < 0.05 ~ "*",
                             p.value < 0.05 ~ ".",
                             TRUE ~ NA_character_),
         Metabolite = factor(Metabolite, levels = met_order_fdr),
         Species = factor(Species, levels = species_order_fdr))

ggplot(grid_dot_fdr, aes(x=Metabolite, y=Species, size=size_val, shape=direction, color=direction)) +
  cowplot::theme_cowplot() + geom_point(alpha=0.9) +
  scale_shape_manual(values=c('Positive'=15,'Negative'=17), name='Direction') +
  scale_color_manual(values=c('Positive'='#BE1E2D','Negative'='darkblue'), name='Direction') +
  scale_x_discrete(drop=FALSE) +
  scale_y_discrete(drop=FALSE) +
  ylab('Microbe') + xlab('Metabolite (persistent)') +
  ggtitle('Microbe–metabolite associations (persistent metabolites, FDR set)') +
  theme(axis.text.x = element_text(angle=60, hjust=1),
        legend.position='bottom', legend.title=element_blank())
p_micro_dot_fdr <- last_plot()
ggsave('./plots/microbe_metabolite_persistent_metabolites_grid_fixed_fdr.pdf', width=14, height=10)

# stack persistent metabolite plot (p_persist) over the microbe grid, drop x labels from grid
p_persist <- ggplot(allsig_fdr %>% mutate(Metabolite = factor(variable, levels = met_order_fdr)),
                    aes(x=Metabolite, y=estimate_r,
                        color=term, shape=type, size=2)) +
  cowplot::theme_cowplot() + geom_point(alpha=0.85) +
  geom_hline(linetype = 'dashed', yintercept = 0) +
  xlab('Metabolite') + ylab('Beta (reversal model)') +
  ggtitle('Persistently altered metabolites post-reversal') +
  theme(legend.position='bottom', legend.title=element_blank(),
        axis.text.x = element_text(angle=60, hjust=1))

p_micro_dot_fdr_nox = p_micro_dot_fdr + theme(axis.title.x = element_blank(),
                                               axis.text.x = element_blank(),
                                               axis.ticks.x = element_blank())
# shared legend, reverse stack order (microbe grid on top, metabolites on bottom) — FDR set
p_persist_noleg = p_persist + theme(legend.position='none')
p_micro_noleg   = p_micro_dot_fdr_nox + theme(legend.position='none')
leg_dir   = cowplot::get_legend(p_micro_dot_fdr + theme(legend.position='bottom'))
leg_diet  = cowplot::get_legend(p_persist + theme(legend.position='bottom'))
stack_body_fdr  = cowplot::plot_grid(p_micro_noleg, p_persist_noleg, ncol=1, align='v', rel_heights=c(1,1.6))
stack_legends_fdr = cowplot::plot_grid(leg_dir, leg_diet, ncol=1, rel_heights=c(1,1))
stack_micro_fdr = cowplot::plot_grid(stack_body_fdr, stack_legends_fdr, ncol=1, rel_heights=c(10,2))
ggsave('./plots/microbe_metabolite_persistent_metabolites_stacked_grid_fdr.pdf', stack_micro_fdr, width=14, height=12)

# shared legend, reverse stack order — nominal p set (lasso grid)
p_persist_p <- NULL
p2_grid_p <- NULL
met_order_p = allsig_p %>% arrange(estimate_r) %>% pull(variable) %>% unique()

speciesdata_p = lasso_long %>%
  filter(Metabolite %in% allsig_p$variable) %>%
  mutate(direction = if_else(Coefficient > 0, "Positive", "Negative"),
         size_val = abs(Coefficient)) %>%
  filter(size_val > 0.1) %>%  # drop near-zero lasso effects
  filter(!is.na(Species), !grepl(' sp', Species))

species_order_p = speciesdata_p %>% arrange(Species) %>% pull(Species) %>% unique()
grid_p = speciesdata_p %>%
  mutate(Metabolite = factor(Metabolite, levels = met_order_p),
         Species = factor(Species, levels = species_order_p))

p_sig_p = pairwise_sig %>%
  filter(Metabolite %in% allsig_p$variable,
         Species %in% species_order_p) %>%
  mutate(sig_lab = case_when(BH < 0.05 ~ "*",
                             p.value < 0.05 ~ ".",
                             TRUE ~ NA_character_),
         Metabolite = factor(Metabolite, levels = met_order_p),
         Species = factor(Species, levels = species_order_p))
p_persist_p <- ggplot(allsig_p %>% mutate(Metabolite = factor(variable, levels = met_order_p)),
                     aes(x=Metabolite, y=estimate_r,
                        color=term, shape=type, size=2)) +
  cowplot::theme_cowplot() + geom_point(alpha=0.85) +
  geom_hline(linetype = 'dashed', yintercept = 0) +
  xlab('Metabolite') + ylab('Beta (reversal model)') +
  ggtitle('Persistently altered metabolites post-reversal (nominal)') +
  theme(legend.position='bottom', legend.title=element_blank(),
        axis.text.x = element_text(angle=60, hjust=1))

p2_grid_p <- grid_p %>% mutate(Metabolite=factor(Metabolite,levels=met_order_p)) %>%
  ggplot(aes(x=Metabolite, y=Species, size=size_val, shape=direction, color=direction)) +
  cowplot::theme_cowplot() + geom_point(alpha=0.8) +
  scale_shape_manual(values=c('Positive'=15,'Negative'=17), name='Direction') +
  scale_color_manual(values=c('Positive'='#BE1E2D','Negative'='darkblue'), name='Direction') +
  scale_x_discrete(drop=FALSE) +
  ylab('Microbe') + xlab('Metabolite') +
  ggtitle('Microbe–metabolite links (nominal, lasso Coeff)') +
  theme(axis.text.x = element_text(angle=60, hjust=1),
        legend.position='bottom', legend.title=element_blank())

p_persist_p_noleg = p_persist_p + theme(legend.position='none')
p_micro_p_noleg   = p2_grid_p + theme(legend.position='none',
                                      axis.title.x = element_blank(),
                                      axis.text.x = element_blank(),
                                      axis.ticks.x = element_blank())
leg_dir_p   = cowplot::get_legend(p2_grid_p + theme(legend.position='bottom'))
leg_diet_p  = cowplot::get_legend(p_persist_p + theme(legend.position='bottom'))
stack_body_p  = cowplot::plot_grid(p_micro_p_noleg, p_persist_p_noleg, ncol=1, align='v', rel_heights=c(1,1.6))
stack_legends_p = cowplot::plot_grid(leg_dir_p, leg_diet_p, ncol=1, rel_heights=c(1,1))
stack_p = cowplot::plot_grid(stack_body_p, stack_legends_p, ncol=1, rel_heights=c(10,2))
ggsave('./plots/microbe_metabolite_persistent_metabolites_stacked_grid_p.pdf', stack_p, width=14, height=12)







