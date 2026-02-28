library(tidyverse)
library(furrr)
library(broom)
library(progressr)
library(ggplot2)
library(dynamicTreeCut)
library(ComplexUpset)
library(ComplexHeatmap)
library(cowplot)
library(reshape2)

minval=0.5

colorsDiet = c(
  Control = "#377eb8",
  Lard = "#e41a1c",
  Milkfat = "#800026",
  Ketogenic = "#ff7f00",
  `Fish Oil` = "#f781bf",
  `Coconut Oil` = "#a65628",
  `Olive Oil` = "#984ea3",
  `Palm Oil` = "#4daf4a"
)


regress <- function(s,data){
  data_sub = data %>% filter(SPECIES == s)
  reg1 = glm(data = data_sub , MEAN_ABUNDANCE ~ DIET*COHORT) 
  reg1 = broom::tidy(reg1) %>% mutate(variable = s)
  return(reg1)
}

setwd('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/')

merged_data = readRDS('intermediate_files/metaphlan4_all_merged_data_complete.rds')
merged_data_cd = merged_data%>% filter(REVERSAL == 'NON_REVERSAL',COHORT =='Cohort1' | COHORT =='Cohort2') %>% filter(grepl('t__',SPECIES)) 
merged_data_cd$SPECIES = as.character(merged_data_cd$SPECIES)
merged_data_cd$COHORT = as.factor(merged_data_cd$COHORT)
merged_data_cd$DIET = factor(merged_data_cd$DIET,levels=c('Control','Coconut Oil','Fish Oil','Ketogenic','Lard','Milkfat','Olive Oil','Palm Oil'))
merged_data_cd$ANIMAL_ID = as.factor(merged_data_cd$ANIMAL_ID)

merged_data_cd = merged_data_cd %>% select(ABUNDANCE,WEIGHT,SEX ,TIME_RECODE_OVERALL ,DIET,ANIMAL_ID,SPECIES,COHORT) %>% mutate(SGB = strsplit(SPECIES,'t__') %>% map_chr(2))

### LOAD DATA
sex = read.csv('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/association_output/metaphlan4_merged_cohort_regression_output_SEX_AUG23.csv') %>% filter(grepl('t__',variable))%>% mutate(SGB = strsplit(variable,'t__') %>% map_chr(2)) %>% mutate(BY = p.adjust(p.value,method = 'BY'))%>% mutate(BH = p.adjust(p.value,method = 'BH'))
sex$term = gsub('SEXM','Sex',sex$term)
sex$term = gsub('DIET:','Diet X ',sex$term)
sex$term = gsub('DIET','',sex$term)
sex$variable = as.character(sex$variable)

metaphlanregs=readRDS('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/association_output/metaphlan4_reversal_associations.rds')
diet = metaphlanregs[[1]] %>% distinct %>% filter(grepl('t__',variable))%>% mutate(SGB = strsplit(variable,'t__') %>% map_chr(2))%>% mutate(BY = p.adjust(p.value,method = 'BY'))%>% mutate(BH = p.adjust(p.value,method = 'BH')) %>% select(-BY,-BONFERRONI)
diet$term = gsub('DIET','',diet$term)
diet$variable = as.character(diet$variable)

longtax = read.csv('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/association_output/metaphlan4_merged_cohort_regression_output_LONGTAX_AUG23.csv') %>% filter(grepl('t__',variable))%>% mutate(SGB = strsplit(variable,'t__') %>% map_chr(2))%>% mutate(BY = p.adjust(p.value,method = 'BY'))%>% mutate(BH = p.adjust(p.value,method = 'BH'))
longtax$term = gsub('DIET','',longtax$term)
longtax$variable = as.character(longtax$variable)

### MERGE IN GTDB SGB MAPPING
gtdbmap = read.delim('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/mpa_vJan21_CHOCOPhlAnSGB_202103_SGB2GTDB.tsv',sep='\t',header=F)
colnames(gtdbmap) = c('SGB','GTDB_TAX')
 
#longtax = left_join(longtax,gtdbmap) %>% mutate(direction = if_else(estimate>0,'Positive','Negative'))
sex = left_join(sex,gtdbmap) %>% mutate(direction = if_else(estimate>0,'Positive','Negative'))
diet = left_join(diet,gtdbmap) %>% mutate(direction = if_else(estimate>0,'Positive','Negative'))

#### JUST COUNTING BUGS THAT WERE ASSOCIATED WITH AT LEAST TWO TIMEPOINTS IN SAME DIRECTION
tokeep = diet %>% mutate(dir = if_else(estimate>0,'Positive','Negative')) %>% filter(BH<0.05) %>% dplyr::group_by(term,dir,SGB) %>% dplyr::count() %>% filter(n>1) %>% select(-n) %>% distinct%>% ungroup %>% select(-dir) %>% group_by(term,SGB) %>% dplyr::count() %>% filter(n==1) %>% select(-n) 
diet = inner_join(diet,tokeep)

longtax = left_join(longtax,gtdbmap) %>% mutate(direction = if_else(estimate>0,'Positive','Negative'))

### SUMMARY BARPLOTS PER REGRESSION -- SEX, TIME, DIET, DIET-TIME

tokeep = diet %>% mutate(dir = if_else(estimate>0,'Positive','Negative')) %>% filter(BH<0.05) %>% dplyr::group_by(term,dir,SGB) %>% dplyr::count() %>% filter(n>1) %>% select(-n) %>% distinct%>% ungroup %>% select(-dir) %>% group_by(term,SGB) %>% dplyr::count() %>% filter(n==1) %>% select(-n) 


diet_counted = inner_join(diet,tokeep) %>% filter(BH<0.05) %>% mutate(dir = if_else(estimate>0,'Positive','Negative')) %>% select(GTDB_TAX,term,dir) %>% distinct%>% group_by(term,dir) %>% dplyr::count() %>% mutate(regtype = 'Diet')
sex_counted = sex %>% filter(term == 'Sex',BH<0.05) %>% mutate(dir = if_else(estimate>0,'Positive','Negative'))  %>% select(term,GTDB_TAX,dir) %>% distinct%>% group_by(term,dir) %>% dplyr::count() %>% mutate(regtype = 'Sex')
sex_diet_counted = sex %>% filter(grepl(':',term),BH<0.05) %>% mutate(dir = if_else(estimate>0,'Positive','Negative'))  %>% select(term,GTDB_TAX,dir) %>% distinct %>%distinct %>% group_by(term,dir) %>% dplyr::count() %>% mutate(regtype = 'Sex X Diet')
time_counted = longtax %>% filter(BH<0.05,term == 'TIME_RECODE_OVERALL') %>% mutate(dir = if_else(estimate>0,'Positive','Negative')) %>% select(term,GTDB_TAX,dir) %>% distinct%>% group_by(term,dir) %>% dplyr::count() %>% mutate(regtype = 'Time')
#diettime_counted = longtax %>% filter(BH<0.1,grepl(':',term)) %>% mutate(dir = if_else(estimate>0,'Positive','Negative')) %>% group_by(term,dir) %>% dplyr::count() %>% mutate(regtype = 'Time X Diet')

summaryplot = bind_rows(diet_counted,sex_counted,time_counted,sex_diet_counted)
summaryplot$term = gsub('Sex','Sex (Male)',summaryplot$term)
summaryplot$term = gsub('TIME_RECODE_OVERALL:','Time X ',summaryplot$term)
summaryplot$term = gsub('Diet:','Diet X ',summaryplot$term)
summaryplot$term = gsub('Sex \\(Male\\):','Sex (Male) X ',summaryplot$term)
summaryplot$term = gsub('TIME_RECODE_OVERALL','Time',summaryplot$term)

summaryplot$term = factor(summaryplot$term,levels = c('Sex (Male)','Time','Coconut Oil','Fish Oil','Ketogenic','Lard','Milkfat','Olive Oil','Palm Oil','Sex (Male) X Coconut Oil','Sex (Male) X Fish Oil','Sex (Male) X Ketogenic','Sex (Male) X Lard','Sex (Male) X Milkfat','Sex (Male) X Olive Oil','Sex (Male) X Palm Oil'))
summaryplot$regtype = factor(summaryplot$regtype,levels = c('Sex','Sex X Diet','Time','Diet','Time X Diet'))


ggplot(summaryplot %>% filter(regtype != 'Sex X Diet'),aes(x = term,y =n,fill = dir)) + geom_bar(stat='identity',position = 'dodge') + facet_grid( . ~ regtype,scales='free_x',space= 'free') + theme_minimal() + theme(axis.text.x = element_text(angle = 60,hjust=1),legend.title=element_blank(),legend.position='none') + ylab('Count') + xlab('') +scale_fill_manual(values = c('darkblue','#BE1E2D'))
ggsave('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/overall_regression_output_AUG23.pdf',width=3,height=3.5)


dietsub = diet %>% filter(BH<0.05) %>% filter(grepl('s__', GTDB_TAX) & !grepl('s__$', GTDB_TAX)) %>% mutate(microbe = strsplit(GTDB_TAX, 's__') %>% map_chr(2))%>% mutate(diet = term) %>% distinct

# upsetR
dietsub_ur = dietsub %>% filter(BH<0.05) %>% mutate(val = 1,direction = if_else(estimate>0,'Positive','Negative')) 
tokeep = diet %>% mutate(dir = if_else(estimate>0,'Positive','Negative')) %>% filter(BH<0.05) %>% dplyr::group_by(term,dir,SGB) %>% dplyr::count() %>% filter(n>1) %>% select(-n) %>% distinct%>% ungroup %>% select(-dir) %>% group_by(term,SGB) %>% dplyr::count() %>% filter(n==1) %>% select(-n) 
dietsub_ur = inner_join(dietsub_ur,tokeep)%>% select(microbe,direction,diet,val) %>% distinct %>% reshape2::dcast(microbe + direction  ~ diet,value.var = 'val') 

### check for consistency in direction across diets
#### THIS ONE -- NEGATIVE IS BLUE

dietsub_ur[is.na(dietsub_ur)] = 0
pdf('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/upset_diet_SEPT23.pdf',width=5,height=4)
#print(ComplexUpset::upset(dietsub_ur,sort_sets=FALSE,set_sizes = FALSE,min_size=2,intersect = colnames(dietsub_ur)[3:ncol(dietsub_ur)],base_annotations=list('Intersection size'=intersection_size(counts=T,mapping=aes(fill=direction))+theme(legend.position = 'None')+scale_fill_manual(values = c('darkblue','#BE1E2D'))),width_ratio=0.1))
ComplexUpset::upset(dietsub_ur %>% group_by(microbe) %>% select(-direction) %>% summarise(across(everything(), sum, na.rm = TRUE)),sort_sets=FALSE,set_sizes = FALSE,min_size=2,intersect = colnames(dietsub_ur)[3:ncol(dietsub_ur)],base_annotations=list('Intersection size'=intersection_size(counts=T)+theme(legend.position = 'None')+scale_fill_manual(values = c('darkblue','#BE1E2D'))),width_ratio=0.1)
dev.off()

saveRDS(dietsub_ur,'~/HMS Dropbox/Braden Tierney/mouse_diet_semir/dietsub_uroverviewdat.rds')

# HEATMAP
bugsofinterest = bind_rows(inner_join(dietsub,tokeep) %>% select(variable,diet,SGB)) %>% distinct %>% rename(DIET = diet,SPECIES = variable)  %>% select(SGB) %>% unique %>% unlist %>% unname
merged_data_cd_sub = merged_data_cd%>% filter(COHORT == 'Cohort2',SGB %in% bugsofinterest) %>% select(-WEIGHT,-COHORT,-SEX) %>% filter(ANIMAL_ID!= 'Cohort2-F-Cont_1-Control-1' & TIME_RECODE_OVERALL!=10) %>% filter(TIME_RECODE_OVERALL>0)#%>% filter(TIME_RECODE_OVERALL == 12)
merged_data_cd_sub = left_join(merged_data_cd_sub,gtdbmap)
merged_data_cd_sub2 = merged_data_cd_sub %>% filter(grepl('s__', GTDB_TAX) & !grepl('s__$', GTDB_TAX)) %>% mutate(SPECIES = strsplit(GTDB_TAX, 's__') %>% map_chr(2)) %>% group_by(SPECIES,ANIMAL_ID,DIET) %>% summarise(ABUNDANCE = mean(ABUNDANCE))  %>% filter(!grepl(' sp',SPECIES))
merged_data_cd_sub_wide = merged_data_cd_sub2 %>% reshape2::dcast(SPECIES + ANIMAL_ID ~ DIET,value.var = 'ABUNDANCE')
merged_data_cd_sub_wide[is.na(merged_data_cd_sub_wide)] = 0

merged_data_cd_sub_wide$`Coconut Oil` = log2((merged_data_cd_sub_wide$`Coconut Oil` + minval) /  (merged_data_cd_sub_wide$`Control`+minval))
merged_data_cd_sub_wide$`Ketogenic` = log2((merged_data_cd_sub_wide$`Ketogenic` + minval) /  (merged_data_cd_sub_wide$`Control`+minval))
merged_data_cd_sub_wide$`Olive Oil` = log2((merged_data_cd_sub_wide$`Olive Oil` + minval) /  (merged_data_cd_sub_wide$`Control`+minval))
merged_data_cd_sub_wide$`Palm Oil` = log2((merged_data_cd_sub_wide$`Palm Oil` + minval) /  (merged_data_cd_sub_wide$`Control`+minval))
merged_data_cd_sub_wide$`Fish Oil` = log2((merged_data_cd_sub_wide$`Fish Oil` + minval) /  (merged_data_cd_sub_wide$`Control`+minval))
merged_data_cd_sub_wide$`Lard` = log2((merged_data_cd_sub_wide$`Lard` + minval) /  (merged_data_cd_sub_wide$`Control`+minval))
merged_data_cd_sub_wide$`Milkfat` = log2((merged_data_cd_sub_wide$`Milkfat` + minval) /  (merged_data_cd_sub_wide$`Control`+minval))

merged_data_cd_sub_wide = merged_data_cd_sub_wide %>% select(-Control,-ANIMAL_ID)%>% melt() %>% group_by(SPECIES,variable) %>% summarise(value = mean(value))

merged_data_cd_sub_wide_av = merged_data_cd_sub_wide %>% reshape2::dcast(SPECIES ~ variable,value.var = 'value') %>% column_to_rownames("SPECIES")

bugsofinterest = bind_rows(inner_join(dietsub,tokeep) %>% select(microbe,diet,BH)) %>% distinct %>% rename(DIET = diet,SPECIES = microbe) %>% filter(!grepl(' sp',SPECIES)) %>% mutate(ind = if_else(BH<=0.05,"*",'.')) %>% group_by(DIET,SPECIES) %>% slice_min(BH,n=1) %>% select(DIET,SPECIES,ind) %>% distinct %>% group_by(DIET,SPECIES) %>% reshape2::dcast(DIET ~ SPECIES,value.var = 'ind') 
#bugsofinterest = bind_rows(dietsub %>% select(microbe,diet,BH),diettimesub%>% select(microbe,diet,BH)) %>% distinct %>% rename(DIET = diet,SPECIES = microbe) %>% filter(!grepl(' sp',SPECIES)) %>% mutate(ind = if_else(BH<=0.05,"*",'.')) %>% group_by(DIET,SPECIES) %>% slice_min(BH,n=1) %>% select(DIET,SPECIES,ind) %>% distinct %>% group_by(DIET,SPECIES) %>% reshape2::dcast(DIET ~ SPECIES,value.var = 'ind') 
#bugsofinterest$ind[bugsofinterest$BH<=0.1] = "."
#bugsofinterest$ind[bugsofinterest$BH<=0.05] = "*"
bugsofinterest[is.na(bugsofinterest!='*')] = ''
bugsofinterest = bugsofinterest %>% column_to_rownames('DIET')

bugsofinterest = bugsofinterest[colnames(merged_data_cd_sub_wide_av),rownames(merged_data_cd_sub_wide_av)]

col_fun = circlize::colorRamp2(c(-.4,-0.001,0,0.001,.4), c("darkblue",'lightblue',"white","pink","#BE1E2D"))

pdf('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/heatmap_overall_diet.pdf',width=7,height=5.5)
Heatmap(t(merged_data_cd_sub_wide_av),column_names_rot = 60,col=col_fun, cell_fun = function(j, i, x, y, width, height, fill) { grid.text(bugsofinterest[i, j], x, y, gp = gpar(fontsize = 14, col = "black"))})
dev.off()

saveRDS(merged_data_cd_sub_wide_av,'~/HMS Dropbox/Braden Tierney/mouse_diet_semir/heatmapoverviewdat.rds')


# testing a second heatmap vis w/ collapsed species
bugsofinterest = bind_rows(inner_join(dietsub,tokeep) %>% select(variable,diet,SGB)) %>% distinct %>% rename(DIET = diet,SPECIES = variable)  %>% select(SGB) %>% unique %>% unlist %>% unname
merged_data_cd_sub = merged_data_cd%>% filter(COHORT == 'Cohort2',SGB %in% bugsofinterest) %>% select(-WEIGHT,-COHORT,-SEX) %>% filter(ANIMAL_ID!= 'Cohort2-F-Cont_1-Control-1' & TIME_RECODE_OVERALL!=10) %>% filter(TIME_RECODE_OVERALL>0)#%>% filter(TIME_RECODE_OVERALL == 12)
merged_data_cd_sub = left_join(merged_data_cd_sub,gtdbmap)
rawdat = merged_data_cd_sub %>% filter(grepl('s__', GTDB_TAX) & !grepl('s__$', GTDB_TAX)) %>% mutate(SPECIES = strsplit(GTDB_TAX, 's__') %>% map_chr(2)) 
counts = rawdat %>% filter(grepl(' sp',SPECIES)) %>% mutate(Genus = strsplit(SPECIES,' ') %>% map_chr(1)) %>% mutate(temp = paste(Genus,'sp.')) %>% select(temp,SPECIES) %>% distinct %>% group_by(temp) %>% mutate(foo=n()) %>% arrange(temp) %>% filter(foo>1) %>% ungroup %>% select(SPECIES) %>% unlist %>% unname %>% unique
sub1 = rawdat %>% filter(SPECIES %in% counts) %>% filter(grepl(' sp',SPECIES)) %>% mutate(Genus = strsplit(SPECIES,' ') %>% map_chr(1)) %>% mutate(SPECIES = paste(Genus,'sp.')) 
sub1 = sub1 %>% group_by(SPECIES,DIET,TIME_RECODE_OVERALL) %>% summarise(ABUNDANCE = mean(ABUNDANCE))  
merged_data_cd_sub2 = bind_rows(rawdat %>% filter(!(SPECIES %in% counts)) %>% group_by(SPECIES,TIME_RECODE_OVERALL,DIET) %>% summarise(ABUNDANCE = mean(ABUNDANCE)) ,sub1) %>% mutate(SPECIES = gsub('s__','',SPECIES)) %>% filter(SPECIES!='') #%>% filter(TIME_RECODE_OVERALL==12)

merged_data_cd_sub_wide = merged_data_cd_sub2 %>% reshape2::dcast(SPECIES + TIME_RECODE_OVERALL ~ DIET,value.var = 'ABUNDANCE')
merged_data_cd_sub_wide[is.na(merged_data_cd_sub_wide)] = 0

merged_data_cd_sub_wide$`Coconut Oil` = log2((merged_data_cd_sub_wide$`Coconut Oil` + minval) /  (merged_data_cd_sub_wide$`Control`+minval))
merged_data_cd_sub_wide$`Ketogenic` = log2((merged_data_cd_sub_wide$`Ketogenic` + minval) /  (merged_data_cd_sub_wide$`Control`+minval))
merged_data_cd_sub_wide$`Olive Oil` = log2((merged_data_cd_sub_wide$`Olive Oil` + minval) /  (merged_data_cd_sub_wide$`Control`+minval))
merged_data_cd_sub_wide$`Palm Oil` = log2((merged_data_cd_sub_wide$`Palm Oil` + minval) /  (merged_data_cd_sub_wide$`Control`+minval))
merged_data_cd_sub_wide$`Fish Oil` = log2((merged_data_cd_sub_wide$`Fish Oil` + minval) /  (merged_data_cd_sub_wide$`Control`+minval))
merged_data_cd_sub_wide$`Lard` = log2((merged_data_cd_sub_wide$`Lard` + minval) /  (merged_data_cd_sub_wide$`Control`+minval))
merged_data_cd_sub_wide$`Milkfat` = log2((merged_data_cd_sub_wide$`Milkfat` + minval) /  (merged_data_cd_sub_wide$`Control`+minval))
merged_data_cd_sub_wide = merged_data_cd_sub_wide %>% select(-Control,-TIME_RECODE_OVERALL)%>% melt() %>% group_by(SPECIES,variable) %>% summarise(value = mean(value))
merged_data_cd_sub_wide_av = merged_data_cd_sub_wide %>% reshape2::dcast(SPECIES ~ variable,value.var = 'value') %>% column_to_rownames("SPECIES")

bugsofinterest = bind_rows(inner_join(dietsub,tokeep) %>% select(microbe,diet,BH)) %>% distinct %>% rename(DIET = diet,SPECIES = microbe) %>% mutate(ind = if_else(BH<0.05,"*",'.')) %>% group_by(DIET,SPECIES) %>% slice_min(BH,n=1) %>% select(DIET,SPECIES,ind) %>% distinct %>% group_by(DIET,SPECIES)# %>% reshape2::dcast(DIET ~ SPECIES,value.var = 'ind') 
counts = bugsofinterest %>% ungroup %>% filter(grepl(' sp',SPECIES)) %>% mutate(Genus = strsplit(SPECIES,' ') %>% map_chr(1)) %>% mutate(temp = paste(Genus,'sp.')) %>% select(temp,SPECIES) %>% distinct %>% group_by(temp) %>% mutate(foo=n()) %>% arrange(temp) %>% filter(foo>1) %>% ungroup %>% select(SPECIES) %>% unlist %>% unname %>% unique

sub1 = bugsofinterest %>% filter(SPECIES %in% counts) %>% mutate(Genus = strsplit(SPECIES,' ') %>% map_chr(1)) %>% mutate(SPECIES = paste(Genus,'sp.')) %>% select(SPECIES,DIET,ind) 
rawdat2 = bind_rows(bugsofinterest %>% filter(!(SPECIES %in% counts)) %>% select(SPECIES,DIET,ind),sub1) %>% mutate(SPECIES = gsub('s__','',SPECIES)) %>% filter(SPECIES!='')  %>% distinct
sigmat = rawdat2 %>% reshape2::dcast(data = .,SPECIES ~ DIET, value.var = 'ind') %>% filter(SPECIES %in% rownames(merged_data_cd_sub_wide_av))%>% column_to_rownames('SPECIES') 
sigmat[is.na(sigmat!='*')] = ''
sigmat = sigmat %>% t 
sigmat = sigmat[colnames(merged_data_cd_sub_wide_av),rownames(merged_data_cd_sub_wide_av)]
#### ADDRESS CHOICE TO SET 0 LOG2FC IN METHODS TO 0 
sigmat[t(merged_data_cd_sub_wide_av)==0]=''

col_fun = circlize::colorRamp2(c(-5,-0.0001,0,0.0001,3), c("darkblue",'lightblue',"white","pink","#BE1E2D"))

pdf('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/heatmap_overall_diet_2.pdf',width=18,height=5.5)
Heatmap(t(merged_data_cd_sub_wide_av),column_names_rot = 60,col=col_fun, cell_fun = function(j, i, x, y, width, height, fill) { grid.text(sigmat[i, j], x, y, gp = gpar(fontsize = 14, col = "black"))})
dev.off()

saveRDS(merged_data_cd_sub_wide_av,'~/HMS Dropbox/Braden Tierney/mouse_diet_semir/heatmapoverviewdat_2.rds')

#### also add in a dbrda -- diet by direction for F1

toremove = merged_data_cd%>% select(-WEIGHT) %>% distinct %>% group_by(DIET,TIME_RECODE_OVERALL,ANIMAL_ID,SGB) %>% dplyr::count() %>% filter(n>2) %>% ungroup %>% select(ANIMAL_ID) %>% distinct %>% unlist %>% unname

merged_data_cd_wide = reshape2::dcast(data = merged_data_cd %>% filter(!(ANIMAL_ID %in% toremove)) %>% distinct%>% filter(COHORT == 'Cohort2') %>% mutate(temp = paste(ANIMAL_ID,WEIGHT,TIME_RECODE_OVERALL)) %>% filter(temp !='Cohort2-F-Cont_1-Control-1 NA 10'),DIET + TIME_RECODE_OVERALL + ANIMAL_ID ~ SGB,value.var = 'ABUNDANCE')

library(vegan)

temp1 = merged_data_cd_wide %>% mutate(rn = paste(DIET,ANIMAL_ID,TIME_RECODE_OVERALL,sep='-----')) %>% column_to_rownames("rn") %>% select(-DIET,-ANIMAL_ID,-TIME_RECODE_OVERALL)
temp1 = log(temp1 + min(temp1[temp1>0]))
temp2 = merged_data_cd_wide %>% mutate(rn = paste(DIET,ANIMAL_ID,TIME_RECODE_OVERALL,sep='-----')) %>% column_to_rownames("rn") %>% select(DIET)

common_rows <- intersect(rownames(temp1), rownames(temp2))
temp1 <- temp1[common_rows, ]
temp2 <- temp2[common_rows, ,drop=F]

rdaout1 = vegan::dbrda(temp1  ~ .,data = temp2)

site_scores <- vegan::scores(rdaout1)$sites %>% data.frame %>% rownames_to_column('sample_id') %>% mutate(DIET = strsplit(sample_id,'-----') %>% map_chr(1),ANIMAL_ID = strsplit(sample_id,'-----') %>% map_chr(2),TIME_RECODE_OVERALL = strsplit(sample_id,'-----') %>% map_chr(3))
#site_scores = site_scores %>% rownames_to_column('sampleid')
#site_scores = inner_join(site_scores,mdat)

biplot_scores <- vegan::scores(rdaout1)$biplot %>% data.frame
rownames(biplot_scores) = gsub('DIET','',rownames(biplot_scores)) 
biplot_scores = biplot_scores%>% rownames_to_column('var')

centroids = vegan::scores(rdaout1)$centroid %>% data.frame
rownames(centroids) = gsub('DIET','',rownames(centroids)) 
centroids = centroids%>% rownames_to_column('var')

control_centroid1 = centroids %>% filter(var == 'Control') %>% select(dbRDA1) %>% unlist %>% unname
control_centroid2 = centroids %>% filter(var == 'Control') %>% select(dbRDA2) %>% unlist %>% unname

ggplot() + geom_line(data = site_scores,aes(x = dbRDA1, y = dbRDA2,group=ANIMAL_ID),color = 'grey',alpha=.4)+geom_point(data = site_scores %>% data.frame,size=2, aes(x = dbRDA1, y = dbRDA2,color=DIET,alpha=as.numeric(TIME_RECODE_OVERALL))) +theme_minimal() + labs(title = "Bacterial db-RDA", x = "RDA1", y = "RDA2")+ theme(legend.position = 'bottom') + geom_segment(data = biplot_scores %>% data.frame,size=.5, aes(x = 0, y = 0, xend = dbRDA1*1, yend = dbRDA2*1,color = var), arrow = arrow(length = unit(0.1, "inches")))+stat_ellipse(data = site_scores %>% data.frame, aes(x = dbRDA1, y = dbRDA2,color = DIET)) + scale_color_brewer(palette = 'Set1') + scale_color_manual(values = colorsDiet)
ggsave('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/biplot_dbrda_mar24_bacterial.pdf',width=6,height=6)
#### big phylogenetic tree figure with metaphlan species? https://github.com/biobakery/MetaPhlAn/blob/master/metaphlan/utils/mpa_vJan21_CHOCOPhlAnSGB_202103.nwk -- replace heatmap, 

dietfams = diet %>% filter(BH<0.05) %>% select(SGB) %>% unique %>% unlist %>% unname
sexfams = sex%>% select(SGB) %>% unique %>% unlist %>% unname

sgbsofint = c(dietfams) %>% unique

library(ape)
library(ggtree)
library(phylobase)
metatree = read.tree('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/mpa_vJan21_CHOCOPhlAnSGB_202103.nwk')

mapfile = gtdbmap  %>% filter(SGB %in% sgbsofint) %>% separate(GTDB_TAX, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = TRUE, extra = "merge")  %>% filter(Family!= 'f__') 

sgbsrawdat = left_join(merged_data_cd,gtdbmap)  %>% filter(SGB %in% sgbsofint)%>% separate(GTDB_TAX, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = TRUE, extra = "merge") %>% group_by(Family,DIET) %>% summarise(ABUNDANCE = mean(ABUNDANCE)) 

mapfile = mapfile %>% group_by(Family)%>% sample_n(1) 

sgbsrawdat = sgbsrawdat %>% filter(Family %in% mapfile$Family)

mapfile$SGB = gsub('SGB','',mapfile$SGB)
metatree = keep.tip(tip = intersect(mapfile$SGB,metatree$tip.label),phy = metatree)
mapfile= mapfile %>% column_to_rownames('SGB')
mapfile = mapfile[metatree$tip.label,,drop=T]
metatree$tip.label = mapfile$Family

sgbsrawdat_wide = sgbsrawdat %>% reshape2::dcast(Family ~ DIET,value.var = 'ABUNDANCE') %>% column_to_rownames('Family')
sgbsrawdat_wide[metatree$tip.label,] = sgbsrawdat_wide

minval = .1#min(sgbsrawdat_wide[sgbsrawdat_wide>0])

sgbsrawdat_wide$`Coconut Oil` = log2((sgbsrawdat_wide$`Coconut Oil` + minval) /  (sgbsrawdat_wide$`Control`+minval))
sgbsrawdat_wide$`Ketogenic` = log2((sgbsrawdat_wide$`Ketogenic` + minval) /  (sgbsrawdat_wide$`Control`+minval))
sgbsrawdat_wide$`Olive Oil` = log2((sgbsrawdat_wide$`Olive Oil` + minval) /  (sgbsrawdat_wide$`Control`+minval))
sgbsrawdat_wide$`Palm Oil` = log2((sgbsrawdat_wide$`Palm Oil` + minval) /  (sgbsrawdat_wide$`Control`+minval))
sgbsrawdat_wide$`Fish Oil` = log2((sgbsrawdat_wide$`Fish Oil` + minval) /  (sgbsrawdat_wide$`Control`+minval))
sgbsrawdat_wide$`Lard` = log2((sgbsrawdat_wide$`Lard` + minval) /  (sgbsrawdat_wide$`Control`+minval))
sgbsrawdat_wide$`Milkfat` = log2((sgbsrawdat_wide$`Milkfat` + minval) /  (sgbsrawdat_wide$`Control`+minval))

sgbsrawdat_wide = sgbsrawdat_wide %>% select(-Control)
rownames(sgbsrawdat_wide) = gsub('f__','',rownames(sgbsrawdat_wide))

col_fun = circlize::colorRamp2(c(-.4,-0.001,0,0.001,.4), c("darkblue",'lightblue',"white","pink","#BE1E2D"))

metatree$tip.label = gsub('f__','',metatree$tip.label)
mapfile$Family = gsub('f__','',mapfile$Family)
mapfile$Phylum = gsub('p__','',mapfile$Phylum)
tree = ggtree(metatree)
tree$data = left_join(tree$data,mapfile %>% select(Family,Phylum),by = c('label' = 'Family'))
tree0 = tree  + geom_tippoint(aes(color = Phylum),size=4,shape='triangle') + geom_tiplab(aes(label=label))
tree1=gheatmap(tree0, sgbsrawdat_wide, offset = .5, width = .3, colnames_angle = 90, hjust = 1) +scale_x_ggtree() +scale_fill_gradientn(colors = c("darkblue",'lightblue',"white","pink","#BE1E2D"),values = scales::rescale(c(-2,-0.0001,0,.0001,2)), limits = c(-2, 2),guide = "colourbar")
tree2 = tree1 +layout_rectangular()+geom_treescale()#+ ggtree::geom_text(aes(label=label2),color='red')
ggsave(plot=tree2,'~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/family_tree_rect_baselinediet.pdf',width=8,height=7)

##### DIVERSITY ANALYSIS

# GET AND PLOT ALPHA DIVERSITIES
divs = merged_data_cd %>% filter(COHORT == 'Cohort2')%>% select(-SPECIES,-SGB,-COHORT) %>% group_by(across(-c(ABUNDANCE))) %>% summarise(shannon = vegan::diversity(ABUNDANCE, index = "shannon"), simpson = vegan::diversity(ABUNDANCE, index = "simpson"), richness = sum(ABUNDANCE > 0))

ggplot(data = divs,aes(x = TIME_RECODE_OVERALL,y = shannon)) + cowplot::theme_cowplot() + facet_grid(. ~ DIET) + geom_line(aes(group= ANIMAL_ID))+ geom_point()

ggplot(data = divs,aes(x = factor(DIET),y = shannon,fill=factor(DIET))) + cowplot::theme_cowplot() + geom_boxplot() + scale_fill_brewer(palette =  'Set1') + ggtitle('Alpha diversity over time') + facet_grid(. ~ factor(TIME_RECODE_OVERALL),scales='free_x')+ theme(axis.text.x = element_text(angle = 60,hjust=1)) + ylim(0,7) + ggpubr::stat_compare_means(comparisons = list(c('Control','Coconut Oil'),c('Control','Fish Oil'),c('Control','Ketogenic'),c('Control','Lard'),c('Control','Milkfat'),c('Control','Olive Oil'),c('Control','Palm Oil')),method = 'wilcox.test')
ggsave('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/dietbaseline_diversity_overtime_shannon.pdf',width = 10,height=6)

ggplot(data = divs,aes(x = factor(TIME_RECODE_OVERALL),y = shannon,fill=factor(DIET))) + cowplot::theme_cowplot() + geom_boxplot() + scale_fill_brewer(palette =  'Set1') + ggtitle('Alpha diversity over time') + facet_grid(. ~ factor(DIET),scales='free_x')+ theme(axis.text.x = element_text(angle = 60,hjust=1)) + ylim(0,7) + ggpubr::stat_compare_means(comparisons = list(c('Control','Coconut Oil'),c('Control','Fish Oil'),c('Control','Ketogenic'),c('Control','Lard'),c('Control','Milkfat'),c('Control','Olive Oil'),c('Control','Palm Oil')),method = 'wilcox.test')
ggsave('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/dietbaseline_diversity_overtime_shannon.pdf',width = 10,height=6)

bar = divs %>% group_by(DIET,TIME_RECODE_OVERALL) %>% summarise(mean = mean(shannon),sd = sd(shannon))
baz = bar %>% filter(DIET == 'Control') %>% summarise(val = mean(mean))
control_mean = baz$val[[1]]

#### THIS ONE
ggplot(data = bar,aes(x = factor(TIME_RECODE_OVERALL),y = mean,color=factor(DIET)))+geom_point() + facet_wrap(factor(DIET) ~ .,scales='free',nrow=2) +geom_errorbar(aes(ymin= mean-sd,ymax = mean+sd),width=.1) + geom_line(aes(group = paste(DIET)))+ cowplot::theme_cowplot()  + ggtitle('Alpha diversity over time')+ theme(axis.text.x = element_text(hjust=1))  + geom_hline(color='gray',linetype = 'dashed',yintercept = control_mean) + ylab('Shannon Diversity') + xlab('Months')+ scale_color_manual(values=colorsDiet) + theme(legend.position = 'bottom')
ggsave('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/dietbaseline_diversity_overtime_shannon_facet_20240418.pdf',width = 10,height=4)



bar = divs %>% group_by(DIET,TIME_RECODE_OVERALL,SEX) %>% summarise(mean = mean(shannon),sd = sd(shannon))
baz = bar %>% filter(DIET == 'Control') %>% summarise(val = mean(mean))
control_mean = baz$val[[1]]

#### THIS ONE
ggplot(data = bar,aes(x = factor(TIME_RECODE_OVERALL),y = mean,color=factor(DIET)))+geom_point() + facet_wrap(factor(DIET) ~ SEX,scales='free',nrow=2) +geom_errorbar(aes(ymin= mean-sd,ymax = mean+sd),width=.1) + geom_line(aes(group = paste(DIET)))+ cowplot::theme_cowplot()  + ggtitle('Alpha diversity over time')+ theme(axis.text.x = element_text(hjust=1))  + geom_hline(color='gray',linetype = 'dashed',yintercept = control_mean) + ylab('Shannon Diversity') + xlab('Months')+ scale_color_manual(values=colorsDiet) + theme(legend.position = 'bottom')
ggsave('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/dietbaseline_diversity_overtime_shannon_facet_20240418.pdf',width = 12,height=6)


#ggplot(data = divs,aes(x = factor(DIET),y = simpson,fill=factor(DIET))) + cowplot::theme_cowplot() + geom_boxplot() + scale_fill_brewer(palette =  'Set1') + ggtitle('Simpson diversity over time') + facet_grid(. ~ factor(TIME_RECODE_OVERALL),scales='free_x')+ theme(axis.text.x = element_text(angle = 60,hjust=1))  + ggpubr::stat_compare_means(comparisons = list(c('Control','Coconut Oil'),c('Control','Fish Oil'),c('Control','Ketogenic'),c('Control','Lard'),c('Control','Milkfat'),c('Control','Olive Oil'),c('Control','Palm Oil')),method = 'wilcox.test') 
#ggsave('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/dietbaseline_diversity_overtime_simpson.pdf',width = 10,height=6)

ggplot(data = divs,aes(x = factor(DIET),y = richness,fill=factor(DIET))) + cowplot::theme_cowplot() + geom_boxplot() + scale_fill_brewer(palette =  'Set1') + ggtitle('Richness over time') + facet_grid(. ~ factor(TIME_RECODE_OVERALL),scales='free_x')+ theme(axis.text.x = element_text(angle = 60,hjust=1))  + ggpubr::stat_compare_means(comparisons = list(c('Control','Coconut Oil'),c('Control','Fish Oil'),c('Control','Ketogenic'),c('Control','Lard'),c('Control','Milkfat'),c('Control','Olive Oil'),c('Control','Palm Oil')),method = 'wilcox.test')
ggsave('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/dietbaseline_diversity_overtime_richness.pdf',width = 10,height=6)


# GET WITHIN DIET AND TIMEPOINT BETA DIVERSITIES
# beta div 
merged1_sub = merged_data_cd %>% filter(COHORT == 'Cohort2')%>% reshape2::dcast(ANIMAL_ID + WEIGHT + SEX +  TIME_RECODE_OVERALL ~ SGB,value.var ='ABUNDANCE') %>% mutate(sid = paste(ANIMAL_ID,WEIGHT,SEX,TIME_RECODE_OVERALL,sep='-----')) %>% select(-ANIMAL_ID,-TIME_RECODE_OVERALL,-SEX,-WEIGHT) 

### REMOVE WEIRD EMPTY SAMPLE
merged1_sub = merged1_sub %>% filter(sid!='Cohort2-F-Keto_2-Ketogenic-6-----NA-----F-----2') %>% column_to_rownames('sid')

bray = as.data.frame(as.matrix(vegan::vegdist(merged1_sub))) %>% rownames_to_column('sample_id')

bray_melted = melt(bray) #%>% filter(sample_id != variable) 

mdattemp = merged_data_cd %>% select(ANIMAL_ID,WEIGHT,SEX,TIME_RECODE_OVERALL,DIET) %>% distinct %>% mutate(sample_id =paste(ANIMAL_ID,WEIGHT,SEX,TIME_RECODE_OVERALL,sep='-----'))
braym = left_join(bray_melted,mdattemp) %>% dplyr::rename(sample_id1 = sample_id,sample_id = variable)
braym = left_join(braym,mdattemp,by= c('sample_id')) %>% dplyr::rename(sample_id2 = sample_id)

######## ACROSS MOUSE BETADIV
# filter for relevant comparisions
braym_sub =braym %>% filter(DIET.x == DIET.y,TIME_RECODE_OVERALL.x == TIME_RECODE_OVERALL.y,ANIMAL_ID.x != ANIMAL_ID.y)

average_values = braym_sub %>% group_by(DIET.x,TIME_RECODE_OVERALL.x) %>% summarize(avg_value = mean(value, na.rm = TRUE),totalcomps = sum(value>0)) 
average_values$xaxis = paste(average_values$DIET.x,average_values$TIME_RECODE_OVERALL.x)

braym_overall = left_join(braym_sub,average_values)
braym_overall$xaxis = paste(braym_overall$DIET.x,braym_overall$TIME_RECODE_OVERALL.x)

# Reorder the levels of site_name_description.x based on avg_value
braym_overall$xaxis <- factor(braym_overall$xaxis,levels = average_values$xaxis[order(-average_values$avg_value)])
average_values$xaxis <- factor(average_values$xaxis,levels = average_values$xaxis[order(-average_values$avg_value)])

# Create the plot
plot1 = ggplot(data = braym_overall, aes(x = xaxis, y = value,fill=avg_value)) + ggbeeswarm::geom_quasirandom(alpha=.3)+ geom_boxplot(alpha=.7) + cowplot::theme_cowplot() + theme(axis.text.x = element_blank(),axis.title.x = element_blank())+ scale_fill_gradient(low = "blue", high = "red") + ylab('Bray-curtis distance')+ theme(axis.title.y = element_text(angle = 0, vjust = 0.5))+ theme(panel.spacing = unit(0, "lines"), plot.margin = unit(c(0,0,0,0), "lines"))
plot2 = ggplot(data = braym_overall %>% select(xaxis,avg_value) %>% distinct %>% reshape2::melt(id.vars = c('xaxis')) %>% mutate(group = 'a'),aes(x = xaxis,group=variable,y= value)) + geom_line(aes(group = variable))  + cowplot::theme_cowplot()+ theme(axis.text.x = element_blank(),axis.title.x = element_blank()) + ylab('Avg. within-site variation')+ theme(axis.title.y = element_text(angle = 0, vjust = 0.5))+ theme(panel.spacing = unit(0, "lines"), plot.margin = unit(c(0,0,0,0), "lines"))

plot3 = ggplot(data = braym_overall %>% select(xaxis,TIME_RECODE_OVERALL.x) %>% distinct,aes(x = xaxis,y=1,alpha = TIME_RECODE_OVERALL.x)) + geom_tile()  + theme(axis.text.x = element_blank()) + cowplot::theme_cowplot() + ylab('Location type') + theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.title.x = element_blank())+ theme(axis.title.y = element_text(angle = 0, vjust = 0.5),legend.position = 'bottom',legend.text.align = 0,legend.title = element_blank())+ theme(panel.spacing = unit(0, "lines"), plot.margin = unit(c(0,0,0,0), "lines")) + scale_fill_brewer(palette = 'Set1')
plot4 = ggplot(data = braym_overall %>% select(xaxis,DIET.x) %>% distinct,aes(x = xaxis,y=1,fill = DIET.x)) + geom_tile()  + theme(axis.text.x = element_blank()) + cowplot::theme_cowplot() + ylab('Location type') + theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.title.x = element_blank())+ theme(axis.title.y = element_text(angle = 0, vjust = 0.5),legend.position = 'bottom',legend.text.align = 0,legend.title = element_blank())+ theme(panel.spacing = unit(0, "lines"), plot.margin = unit(c(0,0,0,0), "lines"))+ scale_fill_brewer(palette = 'Set1')

combined_plot <- plot1 / (plot2 / plot3 / plot4) + patchwork::plot_layout(heights = c(5, 2.3))

ggsave('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/diet_baseline_across_mouse_betadiv.pdf',width=10,height=8)


##### ACROSS MOUSE BETADIV VISUALIZE DIFFERENTLY
bar = braym_overall %>% group_by(DIET.x,TIME_RECODE_OVERALL.x) %>% summarise(mean = mean(value),sd = sd(value))
baz = bar %>% ungroup %>% filter(TIME_RECODE_OVERALL.x == 0) %>% summarise(val = mean(mean))
control_mean =  baz$val[[1]]

#### THIS ONE

ggplot(data = bar,aes(x = factor(TIME_RECODE_OVERALL.x),y = mean,color=factor(DIET.x)))+geom_point() + facet_wrap(factor(DIET.x) ~ .,scales='fixed',nrow=2) +geom_errorbar(aes(ymin= mean-sd,ymax = mean+sd),width=.1) + geom_line(aes(group = paste(DIET.x)))+ cowplot::theme_cowplot()  + ggtitle('Intra-mouse dissimilarity over time')+ theme(axis.text.x = element_text(hjust=1))  + geom_hline(color='gray',linetype = 'dashed',yintercept = control_mean) + ylab('Beta Diversity') + xlab('Months')+ scale_color_manual(values=colorsDiet)+ theme(legend.position = 'bottom')+ ylim(0,1)
ggsave('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/dietbeta_diversity_overtime_shannon_facet_20240418.pdf',width = 10,height=4)


##### ACROSS MOUSE BETADIV VISUALIZE DIFFERENTLY
bar = braym_overall %>% group_by(DIET.x,TIME_RECODE_OVERALL.x,SEX.x) %>% summarise(mean = mean(value),sd = sd(value))
baz = bar %>% ungroup %>% filter(TIME_RECODE_OVERALL.x == 0) %>% summarise(val = mean(mean))
control_mean =  baz$val[[1]]

#### THIS ONE

ggplot(data = bar,aes(x = factor(TIME_RECODE_OVERALL.x),y = mean,color=factor(DIET.x)))+geom_point() + facet_wrap(factor(DIET.x) ~ SEX.x,scales='fixed',nrow=2) +geom_errorbar(aes(ymin= mean-sd,ymax = mean+sd),width=.1) + geom_line(aes(group = paste(DIET.x)))+ cowplot::theme_cowplot()  + ggtitle('Intra-mouse dissimilarity over time')+ theme(axis.text.x = element_text(hjust=1))  + geom_hline(color='gray',linetype = 'dashed',yintercept = control_mean) + ylab('Beta Diversity') + xlab('Months')+ scale_color_manual(values=colorsDiet)+ theme(legend.position = 'bottom') + ylim(0,1)
ggsave('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/dietbeta_diversity_overtime_shannon_facet_20240418.pdf',width = 12,height=6)




ggplot(data = braym_overall, aes(x = factor(TIME_RECODE_OVERALL.x), y = value)) + ggbeeswarm::geom_quasirandom(alpha=.3) + geom_boxplot(alpha=.7,aes(color=avg_value)) + cowplot::theme_cowplot() + facet_grid(. ~ DIET.x) + scale_fill_gradient(low = "blue", high = "red") + ylab('Bray-curtis distance') + xlab('Timepoint') + theme(axis.title.y = element_text(angle = 90, vjust = 0.5))
ggsave('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/diet_baseline_across_mouse_betadiv_facetplot.pdf',width=10,height=4)


######## WITHIN MOUSE BETADIV
# filter for relevant comparisions
braym_sub =braym %>% filter(DIET.x == DIET.y,ANIMAL_ID.x == ANIMAL_ID.y)%>% filter(value!=0)

average_values = braym_sub %>% group_by(ANIMAL_ID.x,DIET.x) %>% summarize(avg_value = mean(value, na.rm = TRUE),totalcomps = sum(value>0)) 
average_values$xaxis = paste(average_values$ANIMAL_ID.x)

braym_overall = left_join(braym_sub,average_values)
braym_overall$xaxis = paste(braym_overall$ANIMAL_ID.x)

# Reorder the levels of site_name_description.x based on avg_value
braym_overall$xaxis <- factor(braym_overall$xaxis,levels = average_values$xaxis[order(-average_values$avg_value)])
average_values$xaxis <- factor(average_values$xaxis,levels = average_values$xaxis[order(-average_values$avg_value)])

# Create the plot
plot1 = ggplot(data = braym_overall, aes(x = xaxis, y = value,alpha=(TIME_RECODE_OVERALL.x/12),fill=avg_value)) + ggbeeswarm::geom_quasirandom()+ geom_boxplot(alpha=.7) + cowplot::theme_cowplot() + theme(axis.text.x = element_blank(),axis.title.x = element_blank())+ scale_fill_gradient(low = "blue", high = "red") + ylab('Bray-curtis distance')+ theme(axis.title.y = element_text(angle = 0, vjust = 0.5))+ theme(panel.spacing = unit(0, "lines"), plot.margin = unit(c(0,0,0,0), "lines"))
plot2 = ggplot(data = braym_overall %>% select(xaxis,avg_value) %>% distinct %>% reshape2::melt(id.vars = c('xaxis')) %>% mutate(group = 'a'),aes(x = xaxis,group=variable,y= value)) + geom_line(aes(group = variable))  + cowplot::theme_cowplot()+ theme(axis.text.x = element_blank(),axis.title.x = element_blank()) + ylab('Avg. within-site variation')+ theme(axis.title.y = element_text(angle = 0, vjust = 0.5))+ theme(panel.spacing = unit(0, "lines"), plot.margin = unit(c(0,0,0,0), "lines"))
#plot3 = ggplot(data = braym_overall %>% select(xaxis,TIME_RECODE_OVERALL.x) %>% distinct,aes(x = xaxis,y=1,alpha = TIME_RECODE_OVERALL.x)) + geom_tile()  + theme(axis.text.x = element_blank()) + cowplot::theme_cowplot() + ylab('Location type') + theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.title.x = element_blank())+ theme(axis.title.y = element_text(angle = 0, vjust = 0.5),legend.position = 'bottom',legend.text.align = 0,legend.title = element_blank())+ theme(panel.spacing = unit(0, "lines"), plot.margin = unit(c(0,0,0,0), "lines")) + scale_fill_brewer(palette = 'Set1')
plot4 = ggplot(data = braym_overall %>% select(xaxis,DIET.x) %>% distinct,aes(x = xaxis,y=1,fill = DIET.x)) + geom_tile()  + theme(axis.text.x = element_blank()) + cowplot::theme_cowplot() + ylab('Location type') + theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.title.x = element_blank())+ theme(axis.title.y = element_text(angle = 0, vjust = 0.5),legend.position = 'bottom',legend.text.align = 0,legend.title = element_blank())+ theme(panel.spacing = unit(0, "lines"), plot.margin = unit(c(0,0,0,0), "lines"))+ scale_fill_brewer(palette = 'Set1')

combined_plot <- plot1 / (plot2 / plot4) + patchwork::plot_layout(heights = c(5, 2.3))

ggsave('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/diet_baseline_within_mouse_betadiv.pdf',width=10,height=8)

### alt vis

#### THIS ONE
average_values = average_values %>% group_by(DIET.x) %>% mutate(med2 = median(avg_value))

ggplot(data = average_values, aes(y = fct_reorder(DIET.x,med2), x = avg_value,color=DIET.x))+ ggbeeswarm::geom_quasirandom(size=2) + geom_boxplot(alpha=.7) + cowplot::theme_cowplot() + scale_fill_gradient(low = "blue", high = "red") + xlab('Bray-curtis distance') + ylab('Diet') + theme(axis.title.y = element_text(angle = 90, vjust = 0.5),legend.position = 'none') + ggtitle('Within-Subject Variation') + scale_color_manual(values = colorsDiet)
ggsave('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/diet_baseline_within_mouse_betadiv_singleplot.pdf',width=4,height=4)

# time plot

ggplot(dietsub%>% filter(timepoint>0) %>% group_by(timepoint, diet) %>% count(direction) %>%  mutate(percentage = n / sum(n) * 100), aes(x = factor(timepoint), y = percentage,color=diet)) +geom_line(aes(linetype = direction, group = interaction(direction, diet))) +scale_color_manual(values = colorsDiet) + cowplot::theme_cowplot()

ggsave('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/timeplot_taxa_ondiet.pdf',width=6,height=4)


# time plot reversals --- 

# need a measure of proportion recovered -- how many taxa came back after reversal

# HEATMAP -- REV
diet = metaphlanregs[[2]] %>% distinct %>% filter(grepl('t__',variable))%>% mutate(SGB = strsplit(variable,'t__') %>% map_chr(2))%>% mutate(BY = p.adjust(p.value,method = 'BY'))%>% mutate(BH = p.adjust(p.value,method = 'BH')) %>% select(-BY,-BONFERRONI)
diet$term = gsub('DIET','',diet$term)
diet$variable = as.character(diet$variable)

dietrevsub = diet %>% filter(BH<0.05) %>% filter(grepl('t__',variable)) %>% mutate(direction = if_else(estimate>0,'Positive','Negative'))
dietrevsub$diet = gsub('DIET','',dietrevsub$term)

tokeep = left_join(dietrevsub %>% mutate(SGB = strsplit(variable,'t__') %>% map_chr(2)),gtdbmap) %>% select(term,SGB) %>% distinct

#tokeep = bind_rows(tokeep,tokeep2) %>% distinct

merged_data = readRDS('intermediate_files/metaphlan4_all_merged_data_complete.rds')
merged_data_cd = merged_data%>% filter(REVERSAL == 'REVERSAL' | COHORT =='Cohort2') %>% filter(grepl('t__',SPECIES)) 
merged_data_cd$SPECIES = as.character(merged_data_cd$SPECIES)
merged_data_cd$COHORT = as.factor(merged_data_cd$COHORT)
merged_data_cd$DIET = factor(merged_data_cd$DIET,levels=c('Control','Coconut Oil','Fish Oil','Ketogenic','Lard','Milkfat','Olive Oil','Palm Oil'))
merged_data_cd$ANIMAL_ID = as.factor(merged_data_cd$ANIMAL_ID)

merged_data_cd = merged_data_cd %>% select(REVERSAL,ABUNDANCE,WEIGHT,SEX ,TIME_RECODE_OVERALL ,DIET,ANIMAL_ID,SPECIES,COHORT) %>% mutate(SGB = strsplit(SPECIES,'t__') %>% map_chr(2))

merged_data_cd_sub = merged_data_cd %>% filter(SGB %in% unique(tokeep$SGB)) %>% select(-WEIGHT,-SEX) %>% filter(ANIMAL_ID!= 'Cohort2-F-Cont_1-Control-1' & TIME_RECODE_OVERALL!=10)# %>% filter(TIME_RECODE_OVERALL>0)#%>% filter(TIME_RECODE_OVERALL == 12)

merged_data_cd_sub2_c2 = merged_data_cd_sub %>% filter(COHORT == 'Cohort2' & REVERSAL == 'NON_REVERSAL'| COHORT == 'Cohort2' & REVERSAL == 'REVERSAL_4MONTH' & TIME_RECODE_OVERALL<=5 | COHORT == 'Cohort2' & REVERSAL == 'REVERSAL_9MONTH' & TIME_RECODE_OVERALL<=9 ) %>% group_by(DIET,TIME_RECODE_OVERALL,COHORT,SPECIES) %>% summarise(ABUNDANCE = mean(ABUNDANCE)) %>% reshape2::dcast(TIME_RECODE_OVERALL + SPECIES ~ DIET ,value.var = 'ABUNDANCE')
merged_data_cd_sub2_r4 = merged_data_cd_sub %>% filter( COHORT == 'Cohort2' & REVERSAL == 'REVERSAL_4MONTH' & TIME_RECODE_OVERALL>5) %>% group_by(DIET,TIME_RECODE_OVERALL,COHORT,SPECIES) %>% summarise(ABUNDANCE = mean(ABUNDANCE)) %>% reshape2::dcast(TIME_RECODE_OVERALL + SPECIES ~ DIET,value.var = 'ABUNDANCE')
merged_data_cd_sub2_r9 = merged_data_cd_sub %>% filter( COHORT == 'Cohort2' & REVERSAL == 'REVERSAL_9MONTH' & TIME_RECODE_OVERALL>9 ) %>% group_by(DIET,TIME_RECODE_OVERALL,COHORT,SPECIES) %>% summarise(ABUNDANCE = mean(ABUNDANCE)) %>% reshape2::dcast(TIME_RECODE_OVERALL + SPECIES ~ DIET,value.var = 'ABUNDANCE')

merged_data_cd_sub2_c2$`Coconut Oil` = log2((merged_data_cd_sub2_c2$`Coconut Oil` + minval) /  (merged_data_cd_sub2_c2$`Control`+minval))
merged_data_cd_sub2_c2$`Ketogenic` = log2((merged_data_cd_sub2_c2$`Ketogenic` + minval) /  (merged_data_cd_sub2_c2$`Control`+minval))
merged_data_cd_sub2_c2$`Olive Oil` = log2((merged_data_cd_sub2_c2$`Olive Oil` + minval) /  (merged_data_cd_sub2_c2$`Control`+minval))
merged_data_cd_sub2_c2$`Palm Oil` = log2((merged_data_cd_sub2_c2$`Palm Oil` + minval) /  (merged_data_cd_sub2_c2$`Control`+minval))
merged_data_cd_sub2_c2$`Fish Oil` = log2((merged_data_cd_sub2_c2$`Fish Oil` + minval) /  (merged_data_cd_sub2_c2$`Control`+minval))
merged_data_cd_sub2_c2$`Lard` = log2((merged_data_cd_sub2_c2$`Lard` + minval) /  (merged_data_cd_sub2_c2$`Control`+minval))
merged_data_cd_sub2_c2$`Milkfat` = log2((merged_data_cd_sub2_c2$`Milkfat` + minval) /  (merged_data_cd_sub2_c2$`Control`+minval))

merged_data_cd_sub2_c2 = merged_data_cd_sub2_c2 %>% select(-Control)%>% reshape2::melt(id.vars = c('SPECIES','TIME_RECODE_OVERALL'))%>% mutate(val = 'ON DIET')

merged_data_cd_sub2_r4$`Coconut Oil` = log2((merged_data_cd_sub2_r4$`Coconut Oil` + minval) /  (merged_data_cd_sub2_r4$`Control`+minval))
merged_data_cd_sub2_r4$`Ketogenic` = log2((merged_data_cd_sub2_r4$`Ketogenic` + minval) /  (merged_data_cd_sub2_r4$`Control`+minval))
merged_data_cd_sub2_r4$`Olive Oil` = log2((merged_data_cd_sub2_r4$`Olive Oil` + minval) /  (merged_data_cd_sub2_r4$`Control`+minval))
merged_data_cd_sub2_r4$`Palm Oil` = log2((merged_data_cd_sub2_r4$`Palm Oil` + minval) /  (merged_data_cd_sub2_r4$`Control`+minval))
merged_data_cd_sub2_r4$`Fish Oil` = log2((merged_data_cd_sub2_r4$`Fish Oil` + minval) /  (merged_data_cd_sub2_r4$`Control`+minval))
merged_data_cd_sub2_r4$`Lard` = log2((merged_data_cd_sub2_r4$`Lard` + minval) /  (merged_data_cd_sub2_r4$`Control`+minval))
merged_data_cd_sub2_r4$`Milkfat` = log2((merged_data_cd_sub2_r4$`Milkfat` + minval) /  (merged_data_cd_sub2_r4$`Control`+minval))

merged_data_cd_sub2_r4 = merged_data_cd_sub2_r4 %>% select(-Control)%>% reshape2::melt(id.vars = c('SPECIES','TIME_RECODE_OVERALL'))  %>% mutate(val = '4MR')

merged_data_cd_sub2_r9$`Coconut Oil` = log2((merged_data_cd_sub2_r9$`Coconut Oil` + minval) /  (merged_data_cd_sub2_r9$`Control`+minval))
merged_data_cd_sub2_r9$`Ketogenic` = log2((merged_data_cd_sub2_r9$`Ketogenic` + minval) /  (merged_data_cd_sub2_r9$`Control`+minval))
merged_data_cd_sub2_r9$`Olive Oil` = log2((merged_data_cd_sub2_r9$`Olive Oil` + minval) /  (merged_data_cd_sub2_r9$`Control`+minval))
merged_data_cd_sub2_r9$`Palm Oil` = log2((merged_data_cd_sub2_r9$`Palm Oil` + minval) /  (merged_data_cd_sub2_r9$`Control`+minval))
merged_data_cd_sub2_r9$`Fish Oil` = log2((merged_data_cd_sub2_r9$`Fish Oil` + minval) /  (merged_data_cd_sub2_r9$`Control`+minval))
merged_data_cd_sub2_r9$`Lard` = log2((merged_data_cd_sub2_r9$`Lard` + minval) /  (merged_data_cd_sub2_r9$`Control`+minval))
merged_data_cd_sub2_r9$`Milkfat` = log2((merged_data_cd_sub2_r9$`Milkfat` + minval) /  (merged_data_cd_sub2_r9$`Control`+minval))

merged_data_cd_sub2_r9 = merged_data_cd_sub2_r9 %>% select(-Control) %>% reshape2::melt(id.vars = c('SPECIES','TIME_RECODE_OVERALL')) %>% mutate(val = '9MR')

revhm = bind_rows(merged_data_cd_sub2_c2,merged_data_cd_sub2_r4,merged_data_cd_sub2_r9)

revhm_w = revhm %>% reshape2::dcast(SPECIES + variable ~ val + TIME_RECODE_OVERALL,value.var = 'value') %>% mutate(rn = paste(SPECIES,variable)) %>% column_to_rownames('rn') %>% select(-SPECIES) %>% select(`ON DIET_0`,`ON DIET_2`,`ON DIET_5`,`ON DIET_7`,`ON DIET_12`,`4MR_7`,`4MR_12`,`9MR_12`)

col_fun = circlize::colorRamp2(c(-5,-0.001,0,0.001,5), c("darkblue",'lightblue',"white","pink","#BE1E2D"))

toremove = rowSums(revhm_w) %>% data.frame %>% filter(.==0) %>% rownames
revhm_w  = revhm_w %>% rownames_to_column('foo') %>% filter(!(foo %in% toremove)) %>% column_to_rownames('foo')

pdf('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/overall_rev.pdf',width=6,height=3)
Heatmap(revhm_w,show_row_names = F,show_column_names=T,cluster_columns=F,cluster_column_slices = F,column_split = factor(c('ON DIET','ON DIET','ON DIET','ON DIET','ON DIET','4 MONTH REVERSAL','4 MONTH REVERSAL','9 MONTH REVERSAL'),levels = c('ON DIET','4 MONTH REVERSAL','9 MONTH REVERSAL')),column_names_rot = 0,col=col_fun,name = 'log2(Fold Change)')
dev.off()

###### SUMMARY STATS

# HEATMAP -- REV
diet = metaphlanregs[[2]] %>% distinct %>% filter(grepl('t__',variable))%>% mutate(SGB = strsplit(variable,'t__') %>% map_chr(2))%>% mutate(BY = p.adjust(p.value,method = 'BY'))%>% mutate(BH = p.adjust(p.value,method = 'BH')) %>% select(-BY,-BONFERRONI)
diet$term = gsub('DIET','',diet$term)
diet$variable = as.character(diet$variable)

dietrevsub = diet %>% filter(BH<0.05) %>% filter(grepl('t__',variable)) %>% mutate(direction = if_else(estimate>0,'Positive','Negative'))
dietrevsub$diet = gsub('DIET','',dietrevsub$term)

tokeep = left_join(dietrevsub %>% mutate(SGB = strsplit(variable,'t__') %>% map_chr(2)),gtdbmap) %>% select(term,SGB) %>% distinct

#tokeep = bind_rows(tokeep,tokeep2) %>% distinct

merged_data = readRDS('intermediate_files/metaphlan4_all_merged_data_complete.rds')
merged_data_cd = merged_data%>% filter(REVERSAL == 'REVERSAL' | COHORT =='Cohort2') %>% filter(grepl('t__',SPECIES)) 
merged_data_cd$SPECIES = as.character(merged_data_cd$SPECIES)
merged_data_cd$COHORT = as.factor(merged_data_cd$COHORT)
merged_data_cd$DIET = factor(merged_data_cd$DIET,levels=c('Control','Coconut Oil','Fish Oil','Ketogenic','Lard','Milkfat','Olive Oil','Palm Oil'))
merged_data_cd$ANIMAL_ID = as.factor(merged_data_cd$ANIMAL_ID)

merged_data_cd = merged_data_cd %>% select(REVERSAL,ABUNDANCE,WEIGHT,SEX ,TIME_RECODE_OVERALL ,DIET,ANIMAL_ID,SPECIES,COHORT) %>% mutate(SGB = strsplit(SPECIES,'t__') %>% map_chr(2))

merged_data_cd_sub = merged_data_cd %>% filter(SGB %in% unique(tokeep$SGB)) %>% select(-WEIGHT,-SEX) %>% filter(ANIMAL_ID!= 'Cohort2-F-Cont_1-Control-1' & TIME_RECODE_OVERALL!=10)# %>% filter(TIME_RECODE_OVERALL>0)#%>% filter(TIME_RECODE_OVERALL == 12)

merged_data_cd_sub2_c2 = merged_data_cd_sub %>% filter(COHORT == 'Cohort2' & REVERSAL == 'NON_REVERSAL'| COHORT == 'Cohort2' & REVERSAL == 'REVERSAL_4MONTH' & TIME_RECODE_OVERALL<=5 | COHORT == 'Cohort2' & REVERSAL == 'REVERSAL_9MONTH' & TIME_RECODE_OVERALL<=9 ) %>% group_by(DIET,TIME_RECODE_OVERALL,COHORT,SPECIES) %>% summarise(ABUNDANCE = mean(ABUNDANCE)) %>% reshape2::dcast(TIME_RECODE_OVERALL + SPECIES ~ DIET ,value.var = 'ABUNDANCE')
merged_data_cd_sub2_r4 = merged_data_cd_sub %>% filter( COHORT == 'Cohort2' & REVERSAL == 'REVERSAL_4MONTH' & TIME_RECODE_OVERALL>5) %>% group_by(DIET,TIME_RECODE_OVERALL,COHORT,SPECIES) %>% summarise(ABUNDANCE = mean(ABUNDANCE)) %>% reshape2::dcast(TIME_RECODE_OVERALL + SPECIES ~ DIET,value.var = 'ABUNDANCE')
merged_data_cd_sub2_r9 = merged_data_cd_sub %>% filter( COHORT == 'Cohort2' & REVERSAL == 'REVERSAL_9MONTH' & TIME_RECODE_OVERALL>9 ) %>% group_by(DIET,TIME_RECODE_OVERALL,COHORT,SPECIES) %>% summarise(ABUNDANCE = mean(ABUNDANCE)) %>% reshape2::dcast(TIME_RECODE_OVERALL + SPECIES ~ DIET,value.var = 'ABUNDANCE')

merged_data_cd_sub2_c2$`Coconut Oil` = log2((merged_data_cd_sub2_c2$`Coconut Oil` + minval) /  (merged_data_cd_sub2_c2$`Control`+minval))
merged_data_cd_sub2_c2$`Ketogenic` = log2((merged_data_cd_sub2_c2$`Ketogenic` + minval) /  (merged_data_cd_sub2_c2$`Control`+minval))
merged_data_cd_sub2_c2$`Olive Oil` = log2((merged_data_cd_sub2_c2$`Olive Oil` + minval) /  (merged_data_cd_sub2_c2$`Control`+minval))
merged_data_cd_sub2_c2$`Palm Oil` = log2((merged_data_cd_sub2_c2$`Palm Oil` + minval) /  (merged_data_cd_sub2_c2$`Control`+minval))
merged_data_cd_sub2_c2$`Fish Oil` = log2((merged_data_cd_sub2_c2$`Fish Oil` + minval) /  (merged_data_cd_sub2_c2$`Control`+minval))
merged_data_cd_sub2_c2$`Lard` = log2((merged_data_cd_sub2_c2$`Lard` + minval) /  (merged_data_cd_sub2_c2$`Control`+minval))
merged_data_cd_sub2_c2$`Milkfat` = log2((merged_data_cd_sub2_c2$`Milkfat` + minval) /  (merged_data_cd_sub2_c2$`Control`+minval))

merged_data_cd_sub2_c2 = merged_data_cd_sub2_c2 %>% select(-Control)%>% reshape2::melt(id.vars = c('SPECIES','TIME_RECODE_OVERALL'))%>% mutate(val = 'ON DIET')

merged_data_cd_sub2_r4$`Coconut Oil` = log2((merged_data_cd_sub2_r4$`Coconut Oil` + minval) /  (merged_data_cd_sub2_r4$`Control`+minval))
merged_data_cd_sub2_r4$`Ketogenic` = log2((merged_data_cd_sub2_r4$`Ketogenic` + minval) /  (merged_data_cd_sub2_r4$`Control`+minval))
merged_data_cd_sub2_r4$`Olive Oil` = log2((merged_data_cd_sub2_r4$`Olive Oil` + minval) /  (merged_data_cd_sub2_r4$`Control`+minval))
merged_data_cd_sub2_r4$`Palm Oil` = log2((merged_data_cd_sub2_r4$`Palm Oil` + minval) /  (merged_data_cd_sub2_r4$`Control`+minval))
merged_data_cd_sub2_r4$`Fish Oil` = log2((merged_data_cd_sub2_r4$`Fish Oil` + minval) /  (merged_data_cd_sub2_r4$`Control`+minval))
merged_data_cd_sub2_r4$`Lard` = log2((merged_data_cd_sub2_r4$`Lard` + minval) /  (merged_data_cd_sub2_r4$`Control`+minval))
merged_data_cd_sub2_r4$`Milkfat` = log2((merged_data_cd_sub2_r4$`Milkfat` + minval) /  (merged_data_cd_sub2_r4$`Control`+minval))

merged_data_cd_sub2_r4 = merged_data_cd_sub2_r4 %>% select(-Control)%>% reshape2::melt(id.vars = c('SPECIES','TIME_RECODE_OVERALL'))  %>% mutate(val = '4MR')

merged_data_cd_sub2_r9$`Coconut Oil` = log2((merged_data_cd_sub2_r9$`Coconut Oil` + minval) /  (merged_data_cd_sub2_r9$`Control`+minval))
merged_data_cd_sub2_r9$`Ketogenic` = log2((merged_data_cd_sub2_r9$`Ketogenic` + minval) /  (merged_data_cd_sub2_r9$`Control`+minval))
merged_data_cd_sub2_r9$`Olive Oil` = log2((merged_data_cd_sub2_r9$`Olive Oil` + minval) /  (merged_data_cd_sub2_r9$`Control`+minval))
merged_data_cd_sub2_r9$`Palm Oil` = log2((merged_data_cd_sub2_r9$`Palm Oil` + minval) /  (merged_data_cd_sub2_r9$`Control`+minval))
merged_data_cd_sub2_r9$`Fish Oil` = log2((merged_data_cd_sub2_r9$`Fish Oil` + minval) /  (merged_data_cd_sub2_r9$`Control`+minval))
merged_data_cd_sub2_r9$`Lard` = log2((merged_data_cd_sub2_r9$`Lard` + minval) /  (merged_data_cd_sub2_r9$`Control`+minval))
merged_data_cd_sub2_r9$`Milkfat` = log2((merged_data_cd_sub2_r9$`Milkfat` + minval) /  (merged_data_cd_sub2_r9$`Control`+minval))

merged_data_cd_sub2_r9 = merged_data_cd_sub2_r9 %>% select(-Control) %>% reshape2::melt(id.vars = c('SPECIES','TIME_RECODE_OVERALL')) %>% mutate(val = '9MR')

revhm = bind_rows(merged_data_cd_sub2_c2,merged_data_cd_sub2_r4,merged_data_cd_sub2_r9)

#%>% filter(`4MR` != 0, `ON DIET` != 0, `9MR` != 0)
overallplot = revhm %>% group_by(val,SPECIES,variable) %>% mutate(value = if_else(abs(value)<0.29729974,0,value)) %>% summarise(mean = mean(value)) %>% reshape2::dcast(SPECIES + variable ~ val,value.var = 'mean') %>% mutate(od_dir = if_else(`ON DIET`<0,-1,1)) %>% mutate(od_dir =  if_else(`ON DIET` == 0,0,od_dir)) %>% mutate(dir4 = if_else(`4MR`<0,-1,1))%>% mutate(dir9 = if_else(`9MR`<0,-1,1)) %>% mutate(dir9 =  if_else(`dir9` == 0,0,dir9)) %>% mutate(dir4 =  if_else(`dir4` == 0,0,dir4)) %>% mutate(`4 Month` = if_else(od_dir == dir4,'Persistent','Reverted'),`9 Month` = if_else(od_dir == dir9,'Persistent','Reverted')) %>% mutate(dir = if_else(`ON DIET`<0,'Depletion','Enrichment'))%>% select(dir,variable,`4 Month`,`9 Month`) %>% dplyr::rename(diet = variable) %>% melt(id.vars = c('diet','dir')) 

overallplot_perc = overallplot %>% group_by(diet,variable,value,dir) %>% count() %>% group_by(diet,variable) %>% mutate(foo = sum(n)) %>% mutate(frac = n/foo) %>% group_by(variable,value,dir) %>% mutate(av = mean(frac))

ggplot(overallplot_perc, aes(x = diet, y = frac * 100, fill = value)) + geom_bar(stat = 'identity', position = 'dodge') + facet_grid(dir ~ variable) + theme_cowplot() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1), legend.position = 'bottom') + scale_fill_manual(values = c('black', 'grey')) + ylab('%') + geom_hline(data = overallplot_perc[overallplot_perc$value == 'Persistent' & overallplot_perc$dir == 'Depletion',], aes(yintercept = av * 100), color = 'black', linetype = "dashed") + geom_hline(data = overallplot_perc[overallplot_perc$value == 'Persistent' & overallplot_perc$dir == 'Enrichment',], aes(yintercept = av * 100), color = 'black', linetype = "dashed") + geom_hline(data = overallplot_perc[overallplot_perc$value == 'Reverted' & overallplot_perc$dir == 'Depletion',], aes(yintercept = av * 100), color = 'grey', linetype = "dashed") + geom_hline(data = overallplot_perc[overallplot_perc$value == 'Reverted' & overallplot_perc$dir == 'Enrichment',], aes(yintercept = av * 100), color = 'grey', linetype = "dashed")
ggsave('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/reversal_summary.pdf',width=7,height=4)

### EXTRA

merged_data_cd_sub2 = merged_data_cd_sub %>% mutate(SPECIES = strsplit(SPECIES, 't__') %>% map_chr(2)) %>% group_by(SPECIES,ANIMAL_ID,DIET) %>% summarise(ABUNDANCE = mean(ABUNDANCE))  %>% filter(!grepl(' sp',SPECIES))
merged_data_cd_sub_wide = merged_data_cd_sub2 %>% reshape2::dcast(SPECIES + ANIMAL_ID ~ DIET,value.var = 'ABUNDANCE')
merged_data_cd_sub_wide[is.na(merged_data_cd_sub_wide)] = 0

merged_data_cd_sub_wide$`Coconut Oil` = log2((merged_data_cd_sub_wide$`Coconut Oil` + minval) /  (merged_data_cd_sub_wide$`Control`+minval))
merged_data_cd_sub_wide$`Ketogenic` = log2((merged_data_cd_sub_wide$`Ketogenic` + minval) /  (merged_data_cd_sub_wide$`Control`+minval))
merged_data_cd_sub_wide$`Olive Oil` = log2((merged_data_cd_sub_wide$`Olive Oil` + minval) /  (merged_data_cd_sub_wide$`Control`+minval))
merged_data_cd_sub_wide$`Palm Oil` = log2((merged_data_cd_sub_wide$`Palm Oil` + minval) /  (merged_data_cd_sub_wide$`Control`+minval))
merged_data_cd_sub_wide$`Fish Oil` = log2((merged_data_cd_sub_wide$`Fish Oil` + minval) /  (merged_data_cd_sub_wide$`Control`+minval))
merged_data_cd_sub_wide$`Lard` = log2((merged_data_cd_sub_wide$`Lard` + minval) /  (merged_data_cd_sub_wide$`Control`+minval))
merged_data_cd_sub_wide$`Milkfat` = log2((merged_data_cd_sub_wide$`Milkfat` + minval) /  (merged_data_cd_sub_wide$`Control`+minval))

merged_data_cd_sub_wide = merged_data_cd_sub_wide %>% select(-Control,-ANIMAL_ID)%>% melt() %>% group_by(SPECIES,variable) %>% summarise(value = mean(value))

merged_data_cd_sub_wide_av = merged_data_cd_sub_wide %>% reshape2::dcast(SPECIES ~ variable,value.var = 'value') %>% column_to_rownames("SPECIES")

bugsofinterest = bind_rows(inner_join(dietsub,tokeep) %>% select(microbe,diet,BH)) %>% distinct %>% rename(DIET = diet,SPECIES = microbe) %>% filter(!grepl(' sp',SPECIES)) %>% mutate(ind = if_else(BH<=0.05,"*",'.')) %>% group_by(DIET,SPECIES) %>% slice_min(BH,n=1) %>% select(DIET,SPECIES,ind) %>% distinct %>% group_by(DIET,SPECIES) %>% reshape2::dcast(DIET ~ SPECIES,value.var = 'ind') 
#bugsofinterest = bind_rows(dietsub %>% select(microbe,diet,BH),diettimesub%>% select(microbe,diet,BH)) %>% distinct %>% rename(DIET = diet,SPECIES = microbe) %>% filter(!grepl(' sp',SPECIES)) %>% mutate(ind = if_else(BH<=0.05,"*",'.')) %>% group_by(DIET,SPECIES) %>% slice_min(BH,n=1) %>% select(DIET,SPECIES,ind) %>% distinct %>% group_by(DIET,SPECIES) %>% reshape2::dcast(DIET ~ SPECIES,value.var = 'ind') 
#bugsofinterest$ind[bugsofinterest$BH<=0.1] = "."
#bugsofinterest$ind[bugsofinterest$BH<=0.05] = "*"
bugsofinterest[is.na(bugsofinterest!='*')] = ''
bugsofinterest = bugsofinterest %>% column_to_rownames('DIET')

bugsofinterest = bugsofinterest[colnames(merged_data_cd_sub_wide_av),rownames(merged_data_cd_sub_wide_av)]

col_fun = circlize::colorRamp2(c(-.4,-0.001,0,0.001,.4), c("darkblue",'lightblue',"white","pink","#BE1E2D"))

pdf('~/HMS Dropbox/Braden Tierney/mouse_diet_semir/plots/heatmap_overall_diet.pdf',width=7,height=5.5)
Heatmap(t(merged_data_cd_sub_wide_av),column_names_rot = 60,col=col_fun, cell_fun = function(j, i, x, y, width, height, fill) { grid.text(bugsofinterest[i, j], x, y, gp = gpar(fontsize = 14, col = "black"))})
dev.off()

saveRDS(merged_data_cd_sub_wide_av,'~/HMS Dropbox/Braden Tierney/mouse_diet_semir/heatmapoverviewdat.rds')


