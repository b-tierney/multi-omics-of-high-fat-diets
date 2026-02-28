# pathway vs metabolite corr


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

# Clean metabolite names - remove [POS]>, [NEG]>, and _C0 suffixes
clean_metabolite_name <- function(name) {
  name %>%
    gsub("^\\[POS\\]>", "", .) %>%
    gsub("^\\[NEG\\]>", "", .) %>%
    gsub("_C0$", "", .) %>%
    trimws()
}
library(tidygraph)
library(RColorBrewer)
library(igraph)
library(ComplexHeatmap)
library(ggplot2)

setwd('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/')

### MERGE IN GTDB SGB MAPPING
gtdbmap = read.delim('mpa_vJan21_CHOCOPhlAnSGB_202103_SGB2GTDB.tsv',sep='\t',header=F)
colnames(gtdbmap) = c('SGB','GTDB_TAX')



regress <- function(s,data){
  print(s)
  data_sub = data %>% filter(Metabolite == s)
  reg1 = glm(data = data_sub , ValueLog10 ~ SEX + Diet) 
  reg1 = broom::tidy(reg1) %>% mutate(variable = s)
  return(reg1)
}

regress_microbe_metabolite <- function(p,m,data_p,data_m){
  data_p_sub = data_p %>% filter(PATHWAY == p) %>% select(MEAN_ABUNDANCE,SampleID)
  data_m_sub = data_m %>% filter(Metabolite == m)
  for_regression = inner_join(data_m_sub,data_p_sub,by='SampleID')
  reg1 = glm(data = for_regression , ValueLog10 ~ SEX + MEAN_ABUNDANCE) 
  reg1 = broom::tidy(reg1) %>% mutate(pathway = p,metabolite = m)
  return(reg1)
}

regress_microbe_metabolite2 <- function(m,data_p,data_m,minval){
  print(m)
  data_p_sub = data_p%>% select(PATHWAY,MEAN_ABUNDANCE,SEX,SampleID) %>% dcast(SampleID + SEX~ PATHWAY,value.var ='MEAN_ABUNDANCE')
  data_m_sub = data_m %>% filter(Metabolite == m)
  for_regression = inner_join(data_m_sub %>% select(-SEX),data_p_sub,by='SampleID') %>% select(-rawval)
  for_regression[is.na(for_regression)] = minval
  X <- for_regression %>% select(-Metabolite,-ValueLog10,-SampleID,-Diet) %>% mutate(SEX = as.numeric(as.factor(SEX))-1)
  y <- for_regression$ValueLog10
  preprocessParams<-preProcess(X, method = c("center", "scale"))
  X <- predict(preprocessParams, X)
  lasso <- 'none'
  try({lasso<-coef(cv.glmnet(as.matrix(X),y)) %>% as.matrix %>% as.data.frame %>% dplyr::rename(!!m := lambda.1se) %>% rownames_to_column('features')},silent=TRUE)
  return(lasso)
}

regress_microbe_metabolite3 <- function(m,data_p,data_m,minval){
  print(m)
  data_p_sub = data_p%>% select(PATHWAY,MEAN_ABUNDANCE,SEX,SampleID) %>% dcast(SampleID + SEX~ PATHWAY,value.var ='MEAN_ABUNDANCE')
  data_m_sub = data_m %>% filter(Metabolite == m)
  for_regression = inner_join(data_m_sub,data_p_sub,by='SampleID')
  # for_regression[is.na(for_regression)] = minval
  output <- purrr::map(colnames(for_regression)[6:ncol(for_regression)], function(x) glm(as.formula(paste("ValueLog10 ~", x)), data = for_regression) %>% tidy() %>% filter(term != '(Intercept)') %>% mutate(SGB = x, Metabolite = m)) %>% bind_rows() %>% select(-term)
  return(output)
}

setwd('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/')


##### RUN WITH LASSO METABOLITES AGAINST METAPHLAN4
merged_data = readRDS('intermediate_files/metaphlan4_merged_data_c12.rds') %>% filter(grepl('t__',SPECIES)) %>% mutate(SGB = strsplit(SPECIES,'t__') %>% map_chr(2))
merged_data = merged_data %>% mutate(ANIMAL_ID = paste(CAGE,ANIMAL,sep='-'))
merged_data = merged_data %>% filter(SGB %in% (merged_data %>% dcast(SGB ~ UID,value.var = 'ABUNDANCE') %>% column_to_rownames('SGB') %>% rowSums %>% data.frame %>% filter(.!=0) %>% rownames))
merged_data$DIET = factor(merged_data$DIET,levels=c('Control','Coconut Oil','Fish Oil','Ketogenic','Lard','Milkfat','Olive Oil','Palm Oil'))
merged_data$ANIMAL_ID = as.factor(merged_data$ANIMAL_ID)
merged_data_tpm = merged_data %>% filter(REVERSAL == 'NON_REVERSAL',COHORT == 'Cohort2') %>% filter(TIMEPOINT == 'M6' | TIMEPOINT == 'M9') #%>% mutate(SampleID =  strsplit(' ',ANIMAL_ID) %>% map_chr(2))
merged_data_tpm = merged_data_tpm %>% select(ANIMAL_ID,ABUNDANCE,WEIGHT,SEX ,TIME_RECODE_OVERALL ,DIET,ANIMAL_ID,SGB,COHORT)
merged_data_tpm_avg = merged_data_tpm  %>% dplyr::group_by(ANIMAL_ID,SGB,SEX,DIET) %>% summarise(WEIGHT = mean(WEIGHT,na.rm=T),MEAN_ABUNDANCE = mean(ABUNDANCE,na.rm=T),SD_ABUNDANCE = sd(ABUNDANCE,na.rm=TRUE))
merged_data_tpm_avg$SampleID = gsub('-','_',merged_data_tpm_avg$ANIMAL_ID)
samplemdat = merged_data_tpm_avg %>% select(SampleID,SEX,WEIGHT,DIET) %>% unique
merged_data_tpm_avg = merged_data_tpm_avg %>% ungroup %>% select(-ANIMAL_ID)
counts = merged_data_tpm_avg %>% ungroup %>% select(SGB,MEAN_ABUNDANCE) %>% group_by(SGB) %>% mutate(count = if_else(MEAN_ABUNDANCE>0,1,0))  %>% dplyr::summarise(s = sum(count))
counts = counts %>% filter(s > 3) %>% select(SGB) %>% unlist %>% unname
merged_data_tpm_avg = merged_data_tpm_avg %>% filter(SGB %in% counts)
minval = log10(min(merged_data_tpm_avg$MEAN_ABUNDANCE[merged_data_tpm_avg$MEAN_ABUNDANCE>0]))
merged_data_tpm_avg$MEAN_ABUNDANCE = log10(merged_data_tpm_avg$MEAN_ABUNDANCE + min(merged_data_tpm_avg$MEAN_ABUNDANCE[merged_data_tpm_avg$MEAN_ABUNDANCE>0]))

### GET REPRESENTATIVE TAXA

get_tax_reps <- function(data,correlation_threshold){
  correlation_matrix <- cor(t(data))
  hc <- hclust(as.dist(1 - correlation_matrix), method = "complete")
  groups <- cutree(hc, h = 1 - correlation_threshold)
  unique_groups <- unique(groups)
  representative_species <- vector("list", length(unique_groups))
  for (i in seq_along(unique_groups)) {
    group_indices <- which(groups == unique_groups[i])
    cluster_species <- rownames(data)[group_indices]
    interpretable_species <- cluster_species[!grepl("sp[0-9]+$", cluster_species)]
    if (length(interpretable_species) == 1) {
      representative_species[[i]] <- interpretable_species
      next
    } 
    if (length(cluster_species)>0) {
      max_prevalence_species <- cluster_species[which.max(rowSums(data[group_indices, ] != min(data)))]
      representative_species[[i]] <- max_prevalence_species    
    }
  }
  representative_species = unlist(unname(representative_species))
  representative_species <- representative_species[!grepl("s__$", representative_species)]
  return((representative_species))
}

tmp = left_join(merged_data_tpm_avg,gtdbmap) %>% dcast(GTDB_TAX ~ SampleID,value.var = 'MEAN_ABUNDANCE') %>% column_to_rownames('GTDB_TAX')
repspec = get_tax_reps(tmp,.75)

tokeep= gtdbmap %>% filter(GTDB_TAX %in% repspec) %>% select(SGB) %>% unique %>% unlist %>% unname
merged_data_tpm_avg = merged_data_tpm_avg %>% filter(SGB %in% tokeep)


### LOAD SERUMMET DATA
metabolites = read.csv('data_packet/metabolites/PlasmaMBX_2023_B2_NormLog_long.txt',sep='\t') 
metabolites$SampleID = gsub('_F_M6','',metabolites$SampleID )
metabolites$SampleID = gsub('_F_M11','',metabolites$SampleID )
metabolites$SampleID = gsub('_M_M11','',metabolites$SampleID )
metabolites$SampleID = gsub('_M_M6','',metabolites$SampleID )
metabolites$SampleID = gsub('_M8_2_22_2020','',metabolites$SampleID )
metabolites$SampleID = gsub('Coco_','Coconut_',metabolites$SampleID )
metabolites$SampleID = gsub('Fish.Oil','Fish Oil',metabolites$SampleID)
metabolites = metabolites %>% filter(Time == 'M6' | Time == 'M11') %>% select(-Time) %>% group_by(Metabolite,SampleID,Diet,Sex) %>% summarise(Value = mean(Value)) %>% ungroup %>% dplyr::rename(ValueLog10 = Value) %>% mutate(rawval = 10^(ValueLog10))  %>% select(Metabolite,SampleID,ValueLog10,Sex,Diet,rawval)

metabolites2 = inner_join(metabolites,samplemdat %>% ungroup %>% select(SampleID,SEX) %>% distinct)
metabolites2$Diet = factor(metabolites2$Diet,levels=c('Control','CoconutOil','FishOil','Ketogenic','Lard','MilkFat','OliveOil','PalmOil'))

##### DBRDA

merged_data_cd_wide = dcast(data = metabolites2 ,Diet +SEX + SampleID ~ Metabolite,value.var = 'ValueLog10')

library(vegan)

temp1 = merged_data_cd_wide  %>% mutate(rn = paste(Diet,SampleID,SEX,sep='---'))%>% column_to_rownames("rn") %>% select(-Diet,-SampleID,-SEX)
temp2 = merged_data_cd_wide %>% mutate(rn = paste(Diet,SampleID,SEX,sep='---'))%>% column_to_rownames("rn") %>% select(Diet)

common_rows <- intersect(rownames(temp1), rownames(temp2))
temp1 <- temp1[common_rows, ]
temp2 <- temp2[common_rows, ,drop=F]

rdaout1 = vegan::dbrda(temp1  ~ .,data = temp2,na.rm = TRUE)

site_scores <- vegan::scores(rdaout1)$sites %>% data.frame %>% rownames_to_column('sample_id') %>% mutate(DIET = strsplit(sample_id,'---') %>% map_chr(1),SampleID = strsplit(sample_id,'---') %>% map_chr(2)) %>% select(-sample_id)
#site_scores = site_scores %>% rownames_to_column('sampleid')
#site_scores = inner_join(site_scores,mdat)

biplot_scores <- vegan::scores(rdaout1)$biplot %>% data.frame
rownames(biplot_scores) = gsub('Diet','',rownames(biplot_scores)) 
biplot_scores = biplot_scores%>% rownames_to_column('var')

centroids = vegan::scores(rdaout1)$centroid %>% data.frame
rownames(centroids) = gsub('Diet','',rownames(centroids)) 
centroids = centroids%>% rownames_to_column('var')

control_centroid1 = centroids %>% filter(var == 'Control') %>% select(dbRDA1) %>% unlist %>% unname
control_centroid2 = centroids %>% filter(var == 'Control') %>% select(dbRDA2) %>% unlist %>% unname

ggplot()+geom_point(data = site_scores %>% data.frame,size=2, aes(x = dbRDA1, y = dbRDA2,color=DIET)) +theme_minimal() + labs(title = "Serum Metabolites", x = "RDA1", y = "RDA2")+ theme(legend.position = 'bottom') + geom_segment(data = biplot_scores %>% data.frame,size=.5, aes(x = 0, y = 0, xend = dbRDA1*1, yend = dbRDA2*1,color = var), arrow = arrow(length = unit(0.1, "inches")))+stat_ellipse(data = site_scores %>% data.frame, aes(x = dbRDA1, y = dbRDA2,color = DIET)) + scale_color_brewer(palette = 'Set1') + scale_color_manual(values = colorsDiet)
ggsave('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/biplot_dbrda_SERUMMET.pdf',width=6,height=6)


### filter metabolites based on prevalence -- acting under assumption minimum value is the zero val
tokeep = metabolites2 %>% group_by(Metabolite) %>% summarise(minval = min(rawval),count = sum(rawval>minval)) %>% arrange(count) %>% filter(count>2) %>% select(Metabolite) %>% unlist %>% unname

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
write.csv(output,'association_output/metabolic_pathway_regression_output_dietSERUMMET.csv')

metabolites2 = inner_join(metabolites,samplemdat %>% ungroup %>% select(SampleID,SEX) %>% distinct) %>% select(-Sex)
metabolites2$Diet = factor(metabolites2$Diet,levels=c('Control','CoconutOil','FishOil','Ketogenic','Lard','MilkFat','OliveOil','PalmOil'))

date()
output = purrr::map(unique(metabolites2$Metabolite), function(x) regress(x,metabolites2))
output = bind_rows(output) %>% filter(term != '(Intercept)') %>% mutate(BH = p.adjust(p.value,method="BH"))
output$variable = clean_metabolite_name(output$variable)
date()
print('Finished linear modeling...')
write.csv(output,'association_output/metabolic_metaphlan4_regression_output_dietSERUMMET.csv')

date()
output = purrr::map(unique(metabolites2$Metabolite), function(x) regress_microbe_metabolite2(x,merged_data_tpm_avg %>% rename(PATHWAY=SGB),metabolites2,minval))
output2 = purrr::map(output, function(x) if(typeof(x)!='character'){return(x)})
output2 = output2[which(!sapply(output2, is.null))]
output2 = output2 %>% reduce(left_join,by='features')
date()
print('Finished linear modeling...')

all_output = output2 %>% filter(features!='(Intercept)') %>% column_to_rownames('features')
write.csv(all_output,'association_output/metabolic_metaphlan4_regression_output_lassoSERUMMET.csv')

date()
output = purrr::map(unique(metabolites2$Metabolite), function(x) regress_microbe_metabolite3(x,merged_data_tpm_avg %>% mutate(PATHWAY = SGB) ,metabolites2 %>% select(-SEX,-rawval),minval)) 
output2 = purrr::map(output, function(x) if(typeof(x)!='character'){return(x)})
output2 = output2[which(!sapply(output2, is.null))]
output2 = output2 %>% bind_rows() %>% mutate(BH = p.adjust(p.value))
date()
print('Finished linear modeling...')

write.csv(output2,'~/HMS Dropbox/Braden Tierney/mouse_diet_semir/association_output/metabolic_metaphlan4_regression_output_pairwiseSERUMMET.csv')

#### PLOT METABOLITE OUTPUT BY DIET VOLCANO
res = read.csv('association_output/metabolic_metaphlan4_regression_output_dietSERUMMET.csv') %>% mutate(direction = if_else(estimate>0,'positive',if_else(estimate<0,'negative',NA)))
labelsub1 = res  %>% filter(BH<0.05) %>% filter(!grepl('^X-',variable)) %>% group_by(term,direction) %>% slice_min(estimate,n=5) %>% select(term,variable) %>% mutate(tolabel=TRUE)
labelsub2 = res  %>% filter(BH<0.05) %>% filter(!grepl('^X-',variable)) %>% group_by(term,direction) %>% slice_max(estimate,n=5) %>% select(term,variable) %>% mutate(tolabel=TRUE)

res = left_join(res,bind_rows(labelsub1,labelsub2)) 
res$tolabel[is.na(res$tolabel)] = FALSE

pval_threshold = 0.05
fc_threshold = 1
fc_label_threshold = 1

res = res %>% mutate(Significance = if_else(estimate > 0 & BH<0.05,'Increased','Decreased'))
res$Significance[res$BH>0.05] = 'Insignificant'

res$term  = gsub('Diet','',res$term)
res$term  = gsub('SEXM','Sex',res$term)

# number of associated metabolites
nummet = res %>% filter(BH<0.05) %>% select(term,direction) %>% table %>% data.frame
ggplot(data = nummet,aes(y=fct_reorder(term,Freq),x=Freq,fill=direction)) +cowplot::theme_cowplot()+ geom_bar(stat='identity') + scale_fill_manual(values = c('darkblue','#BE1E2D')) + xlab('# associated metabolites') + theme(legend.title = element_blank(),axis.title.y = element_blank(),legend.position = 'bottom')
ggsave('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/METABOLITE_diet_metabolite_association_countsSERUMMET.pdf',width=3,height=2.5)

# upsetR
dietsub_ur = res %>% filter(BH<0.05) %>% mutate(val = 1,direction = if_else(estimate>0,'positive',if_else(estimate<0,'negative',NA))) %>% reshape2::dcast(variable + direction  ~ term,value.var = 'val') 
dietsub_ur[is.na(dietsub_ur)] = 0
#pdf('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/METABOLITE_upset_dietSERUMMET.pdf',width=6,height=4)
#print(ComplexUpset::upset(dietsub_ur,sort_sets=FALSE,set_sizes = FALSE,min_size=3,intersect = colnames(dietsub_ur)[3:ncol(dietsub_ur)],base_annotations=list('Intersection size'=intersection_size(counts=T,mapping=aes(fill=direction))+theme(legend.position = 'None')+scale_fill_manual(values = c('darkblue','#BE1E2D'))),width_ratio=0.1))
#dev.off()
#png('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/METABOLITE_upset_diet.png')
#print(ComplexUpset::upset(dietsub_ur,sort_sets=FALSE,set_sizes = FALSE,min_size=3,intersect = colnames(dietsub_ur)[3:ncol(dietsub_ur)],base_annotations=list('Intersection size'=intersection_size(counts=T,mapping=aes(fill=direction))+theme(legend.position = 'None')+scale_fill_manual(values = c('darkblue','#BE1E2D'))),width_ratio=0.1))
#dev.off()

###  heatmap 

foo = res %>% filter(BH<0.05)  %>% dcast(term ~ variable,value.var = 'estimate') %>% column_to_rownames('term')
foo[is.na(foo)] = 0

col_fun = circlize::colorRamp2(c(-4,0,4), rev(c("#d73027",'#000000',"#4575b4")))

#Heatmap(foo,show_column_names = F,col=col_fun)

### heatmap, log2fc

merged_data_cd_sub_wide = metabolites2 %>% group_by(Metabolite,Diet) %>% summarise(mean = mean(rawval)) %>% dcast(Metabolite ~ Diet, value.var = 'mean')

minval=.5

colnames(merged_data_cd_sub_wide) = c('Metabolite','Control','Coconut Oil','Fish Oil','Ketogenic','Lard','Milkfat','Olive Oil','Palm Oil')


merged_data_cd_sub_wide$`Coconut Oil` = log2((merged_data_cd_sub_wide$`Coconut Oil` + minval) /  (merged_data_cd_sub_wide$`Control`+minval))
merged_data_cd_sub_wide$`Ketogenic` = log2((merged_data_cd_sub_wide$`Ketogenic` + minval) /  (merged_data_cd_sub_wide$`Control`+minval))
merged_data_cd_sub_wide$`Olive Oil` = log2((merged_data_cd_sub_wide$`Olive Oil` + minval) /  (merged_data_cd_sub_wide$`Control`+minval))
merged_data_cd_sub_wide$`Palm Oil` = log2((merged_data_cd_sub_wide$`Palm Oil` + minval) /  (merged_data_cd_sub_wide$`Control`+minval))
merged_data_cd_sub_wide$`Fish Oil` = log2((merged_data_cd_sub_wide$`Fish Oil` + minval) /  (merged_data_cd_sub_wide$`Control`+minval))
merged_data_cd_sub_wide$`Lard` = log2((merged_data_cd_sub_wide$`Lard` + minval) /  (merged_data_cd_sub_wide$`Control`+minval))
merged_data_cd_sub_wide$`Milkfat` = log2((merged_data_cd_sub_wide$`Milkfat` + minval) /  (merged_data_cd_sub_wide$`Control`+minval))

merged_data_cd_sub_wide = merged_data_cd_sub_wide %>% select(-Control) %>% filter(Metabolite %in% colnames(foo)) 
merged_data_cd_sub_wide_s = merged_data_cd_sub_wide
#merged_data_cd_sub_wide = merged_data_cd_sub_wide %>% column_to_rownames('Metabolite')

#merged_data_cd_sub_wide$RowNames <- ifelse(grepl("^X-", merged_data_cd_sub_wide$Metabolite), "", merged_data_cd_sub_wide$Metabolite)

#rownames(merged_data_cd_sub_wide) <- NULL

col_fun = circlize::colorRamp2(c(-5,-.0001,0,.0001,5), (c('darkblue','lightblue','white','pink',"#d73027")))

pdf('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/metabolite_heatmap_unlabeled_clusteredSERUMMET.pdf',width=14,height=4)
Heatmap(merged_data_cd_sub_wide_s  %>% column_to_rownames('Metabolite') %>% t,col=col_fun,cluster_columns = T,column_names_rot = 60)
dev.off()

# for metabolites associated with a diet, filter for associated bugs based on lasso
ressig = (res %>% filter(BH<0.05) %>% select(variable) %>% unlist %>% unname %>% unique)
metphaphlanlasso = read.csv('association_output/metabolic_metaphlan4_regression_output_lassoSERUMMET.csv',check.names=F) %>% melt() 
colnames(metphaphlanlasso)[1] = 'SGB'
metphaphlanlasso = metphaphlanlasso %>% filter(variable %in% ressig) %>% reshape2::dcast(SGB ~ variable,value.var = 'value') #%>% column_to_rownames('SGB') 
#metphaphlanlasso = metphaphlanlasso[,colSums(metphaphlanlasso)>0]
metphaphlanlasso = left_join(metphaphlanlasso,gtdbmap) %>% separate(GTDB_TAX, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = TRUE, extra = "merge") %>% select(-Domain,-Phylum,-Class,-Family,-Genus,-Class,-Order,-SGB) %>% filter(!is.na(Species)) %>% reshape2::melt(id.vars = c('Species')) %>% filter(Species != 's__')
#rownames(metphaphlanlasso) =gsub(".*\\|", "", rownames(metphaphlanlasso))
#rownames(metphaphlanlasso) = gsub('s__','',rownames(metphaphlanlasso))

bacmostmet = metphaphlanlasso %>% reshape2::melt() %>% data.frame %>% mutate(direction = if_else(value>0,'positive',if_else(value<0,'negative',NA))) %>% group_by(direction) %>% filter(value!=0) # %>% select(SGB) 
#bacmostmet = left_join(bacmostmet,gtdbmap)  %>% separate(GTDB_TAX, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = TRUE, extra = "merge")  %>% filter(Species!= 's__') 


bacmostmet = bacmostmet %>% select(Species,direction) %>% group_by(Species,direction) %>% table %>% data.frame 

metsub = bacmostmet%>% mutate(Species = gsub('s__','',Species)) %>% dplyr::group_by(direction)%>% dplyr::filter(Freq>0)%>% arrange(desc(Freq)) %>% ungroup %>% dplyr::group_by(Species) %>% mutate(sumval = sum(Freq)) 

tokeep = metsub %>% ungroup %>% select(Species,sumval) %>% distinct %>% arrange(desc(sumval)) %>% head(75)

metsub = metsub %>% filter(Species %in% tokeep$Species)

ggplot(data = metsub %>% filter(Freq>1),aes(x = fct_reorder(Species,sumval),y = Freq+1)) + geom_bar(stat='identity') + cowplot::theme_cowplot() + theme(axis.title.x =  element_blank(),axis.text.x = element_text(angle=60,hjust=1)) + ggtitle('Bacteria with the most fecal metabolomic associations') +scale_fill_brewer(palette = 'Set3') + facet_grid(direction ~ .) + scale_y_log10() + theme(legend.position = 'bottom',legend.title = element_blank())#+scale_fill_manual(values = c('darkblue','#BE1E2D')) 
ggsave('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/bacteria_with_most_metabolites_metaphlan4SERUMMET.pdf',width=16,height=5)

bacregout = read.delim('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/association_output/metabolic_metaphlan4_regression_output_pairwiseSERUMMET.csv',sep=',')  %>% mutate(BH = p.adjust(p.value))
gtdbmap = read.delim('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/mpa_vJan21_CHOCOPhlAnSGB_202103_SGB2GTDB.tsv',sep='\t',header=F)
colnames(gtdbmap) = c('SGB','GTDB_TAX')
bacregout = left_join(bacregout,gtdbmap)  %>% separate(GTDB_TAX, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = TRUE, extra = "merge") %>% mutate(direction = if_else(estimate>0,'positive',if_else(estimate<0,'negative',NA)))

bacmostmet2 = bacregout %>% mutate(BH=p.adjust(p.value))%>% filter(BH<0.1) %>% mutate(direction = if_else(estimate>0,'positive',if_else(estimate<0,'negative',NA))) %>% group_by(direction) %>% filter(estimate!=0)  %>% select(Species) %>% table %>% data.frame 
ggplot(data = bacmostmet2%>% mutate(Species = gsub('s__','',Species)) %>% group_by(direction)%>% filter(Freq>0)%>% arrange(desc(Freq)) %>% ungroup %>% group_by(Species) %>% mutate(sumval = sum(Freq)),aes(x = fct_reorder(Species,sumval),y = Freq,fill=direction)) + geom_bar(stat='identity') + cowplot::theme_cowplot() + theme(axis.title.x =  element_blank(),axis.text.x = element_text(angle=60,hjust=1)) + ggtitle('Bacteria with most metabolite associations') +scale_fill_manual(values = c('#BE1E2D','darkblue'))
ggsave('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/bacteria_with_most_metabolites_metaphlan4_FDRSERUMMET.pdf',width=10,height=4)


#metphaphlanlasso  = metphaphlanlasso %>% select(-all_of(colnames(metphaphlanlasso)[grepl("^X-", colnames(metphaphlanlasso))]))
#pdf('~/HMS Dropbox/Braden TierneySuper Mice Diet Omics Project/Metagenomics/Metabolomic_Associations_Nov23/metabolite_microbe_heatmapSERUMMET.pdf',width=16,height=10)
#Heatmap(metphaphlanlasso)
#dev.off()
#png('~/HMS Dropbox/Braden TierneySuper Mice Diet Omics Project/Metagenomics/Metabolomic_Associations_Nov23/metabolite_microbe_heatmap.png')
#Heatmap(metphaphlanlasso)
#dev.off()

#### COMBINED MICROBE METABOLITE HEATMAP --- LOOKING FOR METABOLITES BOTH ASSOCIATED WITH BACTERIA AND DIET -- THEN PLOT A FEW OF THE MOST INTERESTING IN SCATTERPLOTS

# get metabolites of interest

ofinterest = res %>% filter(term!="term") %>% filter(BH<0.1) %>% filter(!grepl('^X-',variable)) %>% filter(term!='Sex')%>% group_by(term) #%>% slice_min(BH,n=20)

metabolites3 = metabolites2 %>% filter(Metabolite %in% unique(ofinterest$variable)) %>% group_by(Metabolite,Diet) %>% summarise(value = mean(rawval)) %>% dcast(Metabolite ~ Diet,value.var = 'value')

metabolites3$`Coconut Oil` = log2((metabolites3$`CoconutOil`) /  (metabolites3$`Control`))
metabolites3$`Ketogenic` = log2((metabolites3$`Ketogenic`) /  (metabolites3$`Control`))
metabolites3$`Olive Oil` = log2((metabolites3$`OliveOil`) /  (metabolites3$`Control`))
metabolites3$`Palm Oil` = log2((metabolites3$`PalmOil`) /  (metabolites3$`Control`))
metabolites3$`Fish Oil` = log2((metabolites3$`FishOil`) /  (metabolites3$`Control`))
metabolites3$`Lard` = log2((metabolites3$`Lard`) /  (metabolites3$`Control`))
metabolites3$`Milkfat` = log2((metabolites3$`MilkFat`) /  (metabolites3$`Control`))

metabolites3 = metabolites3 %>% select(Metabolite,`Coconut Oil`,Ketogenic,`Olive Oil`,`Fish Oil`,`Lard`,`Milkfat`,`Palm Oil`) %>% column_to_rownames('Metabolite')

sigmat = ofinterest %>% dcast(data = .,variable ~ term, value.var = 'BH') %>% filter(variable %in% rownames(metabolites3)) %>% column_to_rownames('variable')

sigmat3=sigmat
sigmat3[sigmat<=0.05] = "."
sigmat3[sigmat>0.05] = ""
sigmat3[is.na(sigmat)] = ''

sigmat3 = sigmat3[rownames(metabolites3),]
sigmat3 = (sigmat3)

sigmat4 = sigmat3 # save for later

col_fun1 = circlize::colorRamp2(c(-5,.1,0,.1,5), rev(c("#d73027",'orange','#000000',"#4575b4",'darkblue')))
pdf('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/metabolite_diet_associations_hmSERUMMET.pdf',width=5,height=12)
Heatmap((metabolites3),row_dend_side = 'left',row_names_side = 'right',col=col_fun1,column_names_rot = 90, cell_fun = function(j, i, x, y, width, height, fill) { grid.text(sigmat3[i, j], x, y, gp = gpar(fontsize = 16, col = "white")) })
dev.off()

# bacterial heatmap


# load in the diet associated microbes

longtax = read.csv('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/association_output/metaphlan4_merged_cohort_regression_output_LONGTAX_AUG23.csv') %>% filter(grepl('t__',variable))%>% mutate(SGB = strsplit(variable,'t__') %>% map_chr(2))%>% mutate(BY = p.adjust(p.value,method = 'BY'))%>% mutate(BH = p.adjust(p.value,method = 'BH'))
longtax$term = gsub('DIET','',longtax$term)
longtax$variable = as.character(longtax$variable)
longtax = left_join(longtax,gtdbmap)
diet = read.csv('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/association_output/metaphlan4_merged_cohort_regression_output_DIET_AUG23.csv')%>% filter(grepl('t__',variable))%>% mutate(SGB = strsplit(variable,'t__') %>% map_chr(2))%>% mutate(BY = p.adjust(p.value,method = 'BY'))%>% mutate(BH = p.adjust(p.value,method = 'BH'))
diet$term = gsub('DIET','',diet$term)
diet$variable = as.character(diet$variable)
diet = left_join(diet,gtdbmap)

dietsub = diet %>% filter(BH<0.1) %>% filter(grepl('s__', GTDB_TAX) & !grepl('s__$', GTDB_TAX)) %>% mutate(microbe = strsplit(GTDB_TAX, 's__') %>% map_chr(2))%>% mutate(diet = term)
diettimesub = longtax %>% filter(grepl(':',term)) %>% mutate(term = gsub('TIME_RECODE_OVERALL:','Time X ',term)) %>% filter(BH<0.1) %>% filter(grepl('s__', GTDB_TAX) & !grepl('s__$', GTDB_TAX)) %>% mutate(microbe = strsplit(GTDB_TAX, 's__') %>% map_chr(2))  %>% mutate(diet = strsplit(term,'X ') %>% map_chr(2))
bugsofinterest = bind_rows(dietsub %>% select(diet,microbe),diettimesub%>% select(diet,microbe)) %>% distinct %>% rename(DIET = diet,Species = microbe) %>% mutate(DIET = as.factor(DIET))

#### generate giant heatmap
bacregout = read.delim('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/association_output/metabolic_metaphlan4_regression_output_pairwiseSERUMMET.csv',sep=',')  %>% mutate(BH = p.adjust(p.value))
gtdbmap = read.delim('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/mpa_vJan21_CHOCOPhlAnSGB_202103_SGB2GTDB.tsv',sep='\t',header=F)
colnames(gtdbmap) = c('SGB','GTDB_TAX')
bacregout = left_join(bacregout,gtdbmap)  %>% separate(GTDB_TAX, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = TRUE, extra = "merge") %>% mutate(direction = if_else(estimate>0,'positive',if_else(estimate<0,'negative',NA)))

tokeep = unique(bugsofinterest$Species)

bacmetsub2 = bacregout %>% mutate(Species = gsub('s__','',Species)) %>% filter(!is.na(Species),Species!='') %>% filter(Species %in% tokeep)%>% filter(Metabolite %in% rownames(metabolites3)) %>% mutate(BH = p.adjust(p.value))  %>% filter(BH<0.1) 

bacmetsub2 = bacmetsub2 %>% filter(Species %in% tokeep)

bacmetsub3= bacmetsub2 %>% dcast(Metabolite ~ Species,value.var = 'estimate') %>% column_to_rownames('Metabolite')
bacmetsub3[is.na(bacmetsub3)] = 0

sigmat = bacmetsub2 %>% dcast(data = .,Metabolite ~ Species, value.var = 'BH') %>% column_to_rownames('Metabolite')

sigmat3=sigmat
sigmat3[sigmat<=0.05] = "*"
sigmat3[sigmat>0.05] = ""
sigmat3[is.na(sigmat)] = ""

sigmat3 = sigmat3[rownames(bacmetsub3),]

col_fun = circlize::colorRamp2(c(-.5,0,.5), c("#d73027",'#000000',"#4575b4"))

metabolites4 = metabolites3 %>% rownames_to_column('Metabolite') %>% filter(Metabolite %in% rownames(bacmetsub3))
metabolites4 = metabolites4 %>% column_to_rownames('Metabolite')
#bacmetsub3 = bacmetsub3 %>% column_to_rownames('Metabolite')
bacmetsub3 = bacmetsub3[rownames(metabolites4),]
#bacmetsub3[is.na(bacmetsub3)] = 0

a = Heatmap(bacmetsub3,row_dend_side = 'left',row_names_side = 'right',col=col_fun,column_names_rot = 90,row_names_rot = 45, cell_fun = function(j, i, x, y, width, height, fill) { grid.text(sigmat3[i, j], x, y, gp = gpar(fontsize = 8, col = "white")) })
b = Heatmap(metabolites4,col=col_fun1,cell_fun = function(j, i, x, y, width, height, fill) { grid.text(sigmat4[i, j], x, y, gp = gpar(fontsize = 8, col = "white")) })

pdf('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/bacterial_metabonlite_intersectionSERUMMET.pdf',width=9,height=20)
b+a
dev.off()
#%>% group_by(Species,metabolite) %>% summarise(ABUNDANCE = mean(value))  %>% ungroup %>% filter(Species != 's__') %>% mutate(Species = gsub('s__','',Species)) %>% dcast(metabolite ~ Species, value.var = 'ABUNDANCE') 

#%>% data.frame(check.names=F) %>% t %>% data.frame(check.names=F)%>% rownames_to_column('metabolite') %>% data.frame %>% melt

#bacmetsub = left_join(bacmetsub,gtdbmap,by=c('variable'='SGB')) %>% separate(GTDB_TAX, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = TRUE, extra = "merge") %>% group_by(Species,metabolite) %>% summarise(ABUNDANCE = mean(value))  %>% ungroup %>% filter(Species != 's__') %>% mutate(Species = gsub('s__','',Species)) %>% dcast(metabolite ~ Species, value.var = 'ABUNDANCE') 

##### look at just diet associated microbes and their relevant metabolites



# filter for bugs of interest
bacregoutsub = bacregout %>% mutate(Species = gsub('s__','',Species)) %>% filter(Species %in% unique(bugsofinterest$Species)) %>% filter(!grepl('^X-',Metabolite))%>% mutate(BH = p.adjust(p.value)) 

forhm = bacregoutsub %>% filter(BH<.1) %>% dcast(Species ~ Metabolite,value.var = 'estimate') %>% column_to_rownames('Species') 
forhm[is.na(forhm)] = 0

Heatmap(t(forhm))

# try a metabolite heatmap for the lasso %>% filter(Species %in% bugsofinterest$Species)

tmp = metphaphlanlasso %>% mutate(Species = gsub('s__','',Species)) %>% filter(!grepl('^X-',variable)) %>% group_by(variable) %>% mutate(val = if_else(value>0,1,0),total = val) %>% filter(total!=0) %>% dcast(Species ~ variable,value.var = 'value')%>% column_to_rownames("Species") 
tmp[is.na(tmp)] = 0

pdf('~/testSERUMMET.pdf',width=40,height=40)
Heatmap(tmp)
dev.off()

library(ggraph)

##### CEM VOLCANO PLOT

qCutoff = 0.05
estimateCutoff = 0.0 # cutoff for estimate, not log2FC
minSignificance = 10
showTopN = 8
maxOverlaps = 10

dirColors = c("Increased"="#ff964a", "Insignificant"="#999999", "Decreased"="#ae91e3")
dirSize = c("Increased"=3, "Decreased"=3, "Insignificant"=1)
dirShape = c(Increased=21, Decreased=21, Insignificant=16)
labelMetabolites = c() # metabolites to manually label

#res = res[ res$term != "Sex", ]
res$Metabolite = res$variable
res$Significance = "Insignificant"
res$Significance[ which(res$BH < qCutoff & res$estimate > estimateCutoff)] = "Increased"
res$Significance[ which(res$BH < qCutoff & res$estimate < -estimateCutoff)] = "Decreased"


resLabel = data.frame()
for(curComp in unique(res$term))
{	
  curRes = res[ res$term %in% curComp, ]
  curName = curComp
  
  # label only vars that are named metabolites
  resFilter = curRes[ !grepl("^X-", curRes$Metabolite), ]
  
  # Get top and bottom metabolites by log2FC and/or pval
  resUp = resFilter[ which(resFilter$Significance == "Increased"), ]
  resDn = resFilter[ which(resFilter$Significance == "Decreased"), ]
  topMetabolitesUp = min(showTopN, nrow(resUp))
  topMetabolitesDn = min(showTopN, nrow(resDn))
  topMetabolites = resUp[order(resUp$estimate, decreasing=T), ][1:topMetabolitesUp,]
  topMetabolites = rbind(topMetabolites, resUp[order(resUp$BH), ][1:topMetabolitesUp,])
  topMetabolites = rbind(topMetabolites, resDn[order(resDn$estimate, decreasing=F), ][1:topMetabolitesDn,])
  topMetabolites = rbind(topMetabolites, resDn[order(resDn$BH), ][1:topMetabolitesDn,])
  topMetabolites = rbind(topMetabolites, curRes[ curRes$Metabolite %in% labelMetabolites, ])
  topMetabolites = unique(topMetabolites)
  topMetabolites = topMetabolites[ !is.na(topMetabolites$BH), ]
  resLabel = rbind(resLabel, topMetabolites)
}

nudgeMax = 3 #max(-log10(resLabel$padj), na.rm=T) + 1
nudger = function(x, lim) { return( ifelse( (lim-x) > 0, (lim-x), 0) ) }
ggplot(res, aes(x=estimate, y=-log10(BH), fill=Significance)) + geom_hline(yintercept=0) + 
  geom_hline(yintercept=-log10(qCutoff), linetype="dashed") + 
  geom_vline(xintercept=c(-estimateCutoff, estimateCutoff), linetype="dashed") +  
  geom_point(aes(size=Significance, shape=Significance)) + 
  geom_label_repel(data=resLabel, aes(label=Metabolite), nudge_y=nudger(-log10(resLabel$BH), nudgeMax), segment.size = 0.5, segment.alpha=0.6, size=3, box.padding = 0.1, label.padding = 0.1, label.r = 0.05,  max.overlaps = maxOverlaps) + 
  scale_size_manual(values=dirSize) + 
  scale_fill_manual(name="Metabolite Direction", values=dirColors) + 
  scale_shape_manual(values=dirShape) + 
  scale_x_continuous(name="log2 Fold-Change") + 
  scale_y_continuous("-log10( adjusted p-value)", limits=c(0, min(minSignificance, max( -log10(res$BH), na.rm=T))), oob=scales::squish) + 
  facet_wrap(~term, ncol=3) +
  cowplot::theme_cowplot() + theme(legend.position="bottom") +
  ggtitle(paste0("Volcano plot\n-log10(qval) > ", minSignificance, " are squished"))

ggsave("~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/diet_metabolite_volcano_cem_SERUMMET.pdf", width=18, height=18)


