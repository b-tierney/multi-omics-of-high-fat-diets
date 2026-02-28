# clean vs dirty cohort

library(tidyverse)
library(furrr)
library(broom)
library(progressr)
library(ComplexHeatmap)
library(reshape2)
library(ggrepel)
library(patchwork)

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

setwd('~/Dropbox (Mason Lab)/mouse_diet_semir/') 

merged_data = readRDS('intermediate_files/metaphlan4_all_merged_data_complete.rds')
merged_data_cd = merged_data%>% filter(REVERSAL == 'NON_REVERSAL',COHORT =='Cohort2_Dirty' | COHORT =='Cohort2') %>% filter(grepl('t__',SPECIES))
merged_data_cd$SPECIES = as.character(merged_data_cd$SPECIES)
merged_data_cd$COHORT = as.factor(merged_data_cd$COHORT)


merged_data_CLEAN = merged_data%>% filter(REVERSAL == 'NON_REVERSAL',COHORT =='Cohort2') %>% filter(grepl('t__',SPECIES)) 
merged_data_CLEAN$SPECIES = as.character(merged_data_CLEAN$SPECIES)
merged_data_CLEAN$COHORT = as.factor(merged_data_CLEAN$COHORT)
merged_data_CLEAN$DIET = factor(merged_data_CLEAN$DIET,levels=c('Control','Coconut Oil','Fish Oil','Ketogenic','Lard','Milkfat','Olive Oil','Palm Oil'))
merged_data_CLEAN$ANIMAL_ID = as.factor(merged_data_CLEAN$ANIMAL_ID)

# upset

ups = merged_data %>% filter(REVERSAL == 'NON_REVERSAL')%>% filter(grepl('t__',SPECIES))%>% filter(ABUNDANCE>0) %>% filter(COHORT == 'Cohort2' | COHORT == 'Cohort2_Dirty') %>% mutate(pres = if_else(ABUNDANCE>0,1,0)) %>% select(SPECIES,COHORT,pres) %>% distinct %>% dcast(SPECIES ~ COHORT,value.var = 'pres') %>% column_to_rownames('SPECIES')
ups[is.na(ups)] = 0

pdf('~/Dropbox (Mason Lab)/mouse_diet_semir/plots/upset_cleandirt.pdf',width=6,height=3)
ComplexUpset::upset(ups,intersect = colnames(ups))
dev.off()

#### longitudinal analysis

library(lme4)
library(lmerTest)

library(broom.mixed)

myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}

regress_lmer_cleandirtoverall <- function(s,data){
  data_sub = data %>% filter(SPECIES == s) %>% filter(DIET == 'Control')
  reg1 = lmer(data = data_sub , ABUNDANCE ~ COHORT + (1|ANIMAL_ID)) 
  conv_message=as.character(summary(reg1)[['optinfo']]$conv$lme4$messages)
  warning = as.character(summary(reg1)[['optinfo']]$warnings)
  if(length(conv_message)==0){
    conv_message=''
  }
  if(length(warning)==0){
    warning=''
  }
  reg1 = broom.mixed::tidy(reg1) %>% mutate(variable = s,warn = list(warning), conv_message = list(conv_message))
  return(reg1)
}

regress_lmer_cleandirtdietcont <- function(s,data){
  data_sub = data %>% filter(SPECIES == s,COHORT == 'Cohort2_Dirty')
  data_sub = data_sub %>% select(ABUNDANCE,DIET,ANIMAL_ID) %>% group_by(ANIMAL_ID) %>% mutate(ABUNDANCE = mean(ABUNDANCE))
  reg1 = lm(data = data_sub, (ABUNDANCE) ~ DIET) 
  #reg1 = lm(data = data_sub, ABUNDANCE ~ DIET) 
  # conv_message=as.character(summary(reg1)[['optinfo']]$conv$lme4$messages)
#  warning = as.character(summary(reg1)[['optinfo']]$warnings)
 # if(length(conv_message)==0){
  #  conv_message=''
  #}
  #if(length(warning)==0){
  #  warning=''
  #}
  #reg1 = broom.mixed::tidy(reg1) %>% mutate(variable = s,warn = list(warning), conv_message = list(conv_message))
  reg1 = broom.mixed::tidy(reg1) %>% mutate(variable = s,warn = '', conv_message = '')
  return(reg1)
}


regress_lmer_cleandirtdiet <- function(s,data){
  data_sub = data %>% filter(SPECIES == s)
  reg1 = lmer(data = data_sub , ABUNDANCE ~ DIET*COHORT + (1|ANIMAL_ID)) 
  conv_message=as.character(summary(reg1)[['optinfo']]$conv$lme4$messages)
  warning = as.character(summary(reg1)[['optinfo']]$warnings)
  if(length(conv_message)==0){
    conv_message=''
  }
  if(length(warning)==0){
    warning=''
  }
  reg1 = broom.mixed::tidy(reg1) %>% mutate(variable = s,warn = list(warning), conv_message = list(conv_message))
  return(reg1)
}

merged_data_cd$DIET = factor(merged_data_cd$DIET,levels=c('Control','Coconut Oil','Fish Oil','Ketogenic','Lard','Milkfat','Olive Oil','Palm Oil'))
merged_data_cd$ANIMAL_ID = as.factor(merged_data_cd$ANIMAL_ID)

merged_data_cd = merged_data_cd %>% select(ABUNDANCE,WEIGHT,SEX ,TIME_RECODE_OVERALL ,DIET,ANIMAL_ID,SPECIES,COHORT)

# filter based on prevalence
counts = merged_data_cd %>% ungroup %>% select(SPECIES,ABUNDANCE) %>% dplyr::group_by(SPECIES) %>% mutate(count = if_else(ABUNDANCE>0,1,0))  %>% dplyr::summarise(s = sum(count))
counts = counts %>% filter(s > 3) %>% select(SPECIES) %>% unlist %>% unname
merged_data_cd = merged_data_cd %>% filter(SPECIES %in% counts)

# log transform
merged_data_cd_UNLOGGED = merged_data_cd
merged_data_cd$ABUNDANCE = merged_data_cd$ABUNDANCE + 0.5
merged_data_cd$ABUNDANCE = log(merged_data_cd$ABUNDANCE)
#plan(multisession, workers = 8)


specieslist = merged_data_cd%>% ungroup %>% select(SPECIES) %>% filter(grepl('t__',SPECIES))%>% distinct %>% unlist %>% unname #%>% filter(!grepl('GGB',SPECIES),!grepl('bacterium',SPECIES))  %>% distinct %>% unlist %>% unname
print('Starting mixed modeling...')
date()
output = purrr::map(specieslist, function(x) suppressMessages(regress_lmer_cleandirtoverall(x,merged_data_cd)))
output1 = purrr::map(specieslist, function(x) suppressMessages(regress_lmer_cleandirtdietcont(x,merged_data_cd)))
output2 = purrr::map(specieslist, function(x) suppressMessages(regress_lmer_cleandirtdiet(x,merged_data_cd)))
date()
print('Finished mixed modeling...')

cleandirtyoverall = output %>% bind_rows() %>% filter(term!='(Intercept)') %>% mutate(BH = p.adjust(p.value),BH= p.adjust(p.value,method = 'BH'))
cleandirtyoverall$term = gsub('COHORT','',cleandirtyoverall$term)
cleandirtyoverall$variable = as.character(cleandirtyoverall$variable)

cleandirtydietjustdirty = output1%>% bind_rows() %>% filter(grepl('DIET',term))%>% mutate(BH = p.adjust(p.value),BH= p.adjust(p.value,method = 'BH'))
cleandirtydietjustdirty$term = gsub('DIET','',cleandirtydietjustdirty$term)
cleandirtydietjustdirty$variable = as.character(cleandirtydietjustdirty$variable)

cleandirtydiet = output2%>% bind_rows() %>% mutate(BH = p.adjust(p.value),BH= p.adjust(p.value,method = 'BH'))
cleandirtydiet$term = gsub('DIET','',cleandirtydiet$term)
cleandirtydiet$term = gsub('COHORT','',cleandirtydiet$term)
cleandirtydiet$variable = as.character(cleandirtydiet$variable)

#### MERGE IN GTDB DATA
cleandirtyoverall = cleandirtyoverall %>% mutate(SGB = strsplit(variable,'t__') %>% map_chr(2))
cleandirtydiet = cleandirtydiet %>% mutate(SGB = strsplit(variable,'t__') %>% map_chr(2))
cleandirtydietjustdirty = cleandirtydietjustdirty %>% mutate(SGB = strsplit(variable,'t__') %>% map_chr(2))

gtdbmap = read.delim('~/Dropbox (Mason Lab)/mouse_diet_semir/mpa_vJan21_CHOCOPhlAnSGB_202103_SGB2GTDB.tsv',sep='\t',header=F)
colnames(gtdbmap) = c('SGB','GTDB_TAX')

cleandirtyoverall = left_join(cleandirtyoverall,gtdbmap) %>% separate(GTDB_TAX, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = TRUE, extra = "merge") 
cleandirtydiet = left_join(cleandirtydiet,gtdbmap) %>% separate(GTDB_TAX, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = TRUE, extra = "merge") 
cleandirtydietjustdirty = left_join(cleandirtydietjustdirty,gtdbmap) %>% separate(GTDB_TAX, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = TRUE, extra = "merge") 
 
write.csv(cleandirtyoverall %>% select(-conv_message,-warn),'~/Dropbox (Mason Lab)/mouse_diet_semir/association_output/clean_dirty_overall_mar24.csv')
write.csv(cleandirtydiet %>% select(-conv_message,-warn),'~/Dropbox (Mason Lab)/mouse_diet_semir/association_output/clean_dirty_WITHINDIET_mar24.csv')
write.csv(cleandirtydietjustdirty %>% select(-conv_message,-warn),'~/Dropbox (Mason Lab)/mouse_diet_semir/association_output/clean_dirty_DIET_mar24.csv')


###### CLEAN DIRTY TREE

library(ape)
library(ggtree)
library(phylobase)

dietfams = cleandirtyoverall %>% filter(BH<0.05) %>% select(SGB) %>% unique %>% unlist %>% unname

sgbsofint = c(dietfams) %>% unique

metatree = read.tree('~/Dropbox (Mason Lab)/mouse_diet_semir/mpa_vJan21_CHOCOPhlAnSGB_202103.nwk')

mapfile = gtdbmap  %>% filter(SGB %in% sgbsofint) %>% separate(GTDB_TAX, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = TRUE, extra = "merge")  %>% filter(Family!= 'f__') 

#### reload merged data
merged_data = readRDS('intermediate_files/metaphlan4_all_merged_data_complete.rds')
merged_data_cd = merged_data%>% filter(REVERSAL == 'NON_REVERSAL',COHORT =='Cohort2_Dirty' | COHORT =='Cohort2') %>% filter(grepl('t__',SPECIES))
merged_data_cd$SPECIES = as.character(merged_data_cd$SPECIES)
merged_data_cd$COHORT = as.factor(merged_data_cd$COHORT)

sgbsrawdat = left_join(merged_data_cd %>% mutate(SGB = strsplit(SPECIES,'t__') %>% map_chr(2)),gtdbmap)  %>% filter(SGB %in% sgbsofint)%>% separate(GTDB_TAX, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = TRUE, extra = "merge") %>% group_by(Family,DIET,COHORT) %>% summarise(ABUNDANCE = mean(ABUNDANCE)) 

mapfile = mapfile %>% group_by(Family)%>% sample_n(1) 

sgbsrawdat = sgbsrawdat %>% filter(Family %in% mapfile$Family)

mapfile$SGB = gsub('SGB','',mapfile$SGB)
metatree = keep.tip(tip = intersect(mapfile$SGB,metatree$tip.label),phy = metatree)
mapfile= mapfile %>% column_to_rownames('SGB')
mapfile = mapfile[metatree$tip.label,,drop=T]
metatree$tip.label = mapfile$Family

sgbsrawdat_wide = sgbsrawdat %>% filter(DIET == 'Control') %>% dcast(Family ~ COHORT,value.var = 'ABUNDANCE') %>% column_to_rownames('Family')
sgbsrawdat_wide = sgbsrawdat_wide[metatree$tip.label,]

minval = 0.5#min(sgbsrawdat_wide[sgbsrawdat_wide>0])

sgbsrawdat_wide$`Clean vs. Dirty` = log2((sgbsrawdat_wide$`Cohort2_Dirty` + minval) /  (sgbsrawdat_wide$`Cohort2` + minval))

#sgbsrawdat_wide$`Coconut Oil` = log2((sgbsrawdat_wide$`Coconut Oil`) /  (sgbsrawdat_wide$`Control`))
#sgbsrawdat_wide$`Ketogenic` = log2((sgbsrawdat_wide$`Ketogenic`) /  (sgbsrawdat_wide$`Control`))
#sgbsrawdat_wide$`Olive Oil` = log2((sgbsrawdat_wide$`Olive Oil`) /  (sgbsrawdat_wide$`Control`))
#sgbsrawdat_wide$`Palm Oil` = log2((sgbsrawdat_wide$`Palm Oil`) /  (sgbsrawdat_wide$`Control`))
#sgbsrawdat_wide$`Fish Oil` = log2((sgbsrawdat_wide$`Fish Oil`) /  (sgbsrawdat_wide$`Control`))
#sgbsrawdat_wide$`Lard` = log2((sgbsrawdat_wide$`Lard`) /  (sgbsrawdat_wide$`Control`))
#sgbsrawdat_wide$`Milkfat` = log2((sgbsrawdat_wide$`Milkfat`) /  (sgbsrawdat_wide$`Control`))

sgbsrawdat_wide = sgbsrawdat_wide %>% select(`Clean vs. Dirty`)

tree = ggtree(metatree)
tree$data = left_join(tree$data,mapfile %>% select(Family,Phylum),by = c('label' = 'Family'))
tree0 = tree  + geom_tippoint(aes(color = Phylum),size=4,shape='triangle') + geom_tiplab(aes(label=label))
tree1=gheatmap(tree0, sgbsrawdat_wide, offset = .5, width = .3, colnames_angle = 90, hjust = 1) +scale_x_ggtree() + scale_fill_viridis_c()+scale_fill_gradientn(colors = c("#d73027", "#f46d43", "#fdae61", "gray", "#abd9e9", "#74add1", "#4575b4"),values = scales::rescale(c(-15, -10,-0.000001, 0, 0.000001, 10, 15)), limits = c(-15, 15),guide = "colourbar")
tree2 = tree1 + layout_rectangular() +geom_treescale()#+ ggtree::geom_text(aes(label=label2),color='red')
ggsave(plot=tree2,'~/Dropbox (Mason Lab)/mouse_diet_semir/plots/family_tree_rect_cleandirty.pdf',width=8,height=10)

##### VOLCANO PLOT
cleandirtyoverall2 = cleandirtyoverall %>% mutate(dir = if_else(estimate>0,1,-1)) %>% filter(effect == 'fixed') %>% filter(Species != 's__')
cleandirtyoverall2 = cleandirtyoverall2 %>% mutate(microbe = strsplit(variable,'s__') %>% map_chr(2))
cleandirtyoverall2$abs_estimate <- abs(cleandirtyoverall2$estimate)
top_microbes <- cleandirtyoverall2 %>% group_by(dir) %>% arrange(desc(abs_estimate)) %>% filter( !grepl('CAG', Species), !grepl('UBA', Species), !grepl('NM', Species), !grepl('1XD', Species)) %>% slice_head(n = 10)
top_microbes2 <- cleandirtyoverall2 %>% group_by(dir) %>% arrange(BH) %>% filter( !grepl('CAG', Species), !grepl('UBA', Species), !grepl('NM', Species), !grepl('1XD', Species)) %>% slice_head(n = 10)
top_microbes = bind_rows(top_microbes, top_microbes2) %>% distinct
cleandirtyoverall2$label_microbe <- gsub('s__', '', ifelse(cleandirtyoverall2$Species %in% top_microbes$Species, cleandirtyoverall2$Species, NA))
ggplot(cleandirtyoverall2, aes(x=estimate, y=-log10(BH))) + geom_point(aes(color=BH), alpha=0.6) + geom_hline(yintercept = -log10(0.05), linetype="dashed", color="red") + geom_text_repel(aes(label=label_microbe), size=3, direction='y', max.overlaps = Inf, nudge_x = .4, nudge_y = .4, segment.color = "gray", segment.size = .5) + theme_minimal() + labs(title="Species increased in the dirty cohort", x="Estimate", y="-log10(BH)") + scale_color_gradient(low="blue", high="red") + theme(legend.position = 'bottom')
ggsave('~/Dropbox (Mason Lab)/mouse_diet_semir/plots/cleandirty_species_volcano.pdf',width=6,height=6)

#### SHARED VS UNIQUE RESPONSES

####### HEATMAP 0 == bugs found only in dirty cohort respond across diet

rawdat = left_join(merged_data_cd %>% filter(COHORT == 'Cohort2_Dirty' | COHORT == 'Cohort2') %>% select(SPECIES,DIET,COHORT,ABUNDANCE) %>% mutate(SGB = strsplit(SPECIES,'t__') %>% map_chr(2)) %>% select(-SPECIES),gtdbmap)  %>% separate(GTDB_TAX, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = TRUE, extra = "merge") %>% ungroup
counts = rawdat %>% filter(grepl(' sp',Species)) %>% mutate(Genus = strsplit(Species,' ') %>% map_chr(1)) %>% mutate(temp = paste(Genus,'sp.')) %>% select(temp,Species) %>% distinct %>% group_by(temp) %>% mutate(foo=n()) %>% arrange(temp) %>% filter(foo>1) %>% ungroup %>% select(Species) %>% unlist %>% unname %>% unique
sub1 = rawdat %>% filter(Species %in% counts)  %>% mutate(Genus = strsplit(Species,' ') %>% map_chr(1)) %>% mutate(Species = paste(Genus,'sp.')) %>% select(Species,DIET,COHORT,ABUNDANCE) %>% mutate(presence = if_else(ABUNDANCE>0,1,0)) 
rawdat2 = bind_rows(rawdat %>% filter(!(Species %in% counts)) %>% select(Species,DIET,COHORT,ABUNDANCE) %>% mutate(presence = if_else(ABUNDANCE>0,1,0)) ,sub1) %>% mutate(Species = gsub('s__','',Species)) %>% filter(Species!='') 
hmdat = rawdat2 %>% select(Species,DIET,ABUNDANCE)
presencespecies = rawdat2 %>% filter(ABUNDANCE>0) %>% select(-ABUNDANCE) %>% distinct%>% dcast(Species ~ COHORT,value.var = 'presence') 
presencespecies_2 = rawdat %>% mutate(presence = if_else(ABUNDANCE>0,1,0))  %>% filter(ABUNDANCE>0) %>% select(-ABUNDANCE) %>% distinct%>% dcast(Species ~ COHORT,value.var = 'presence') 
presencespecies[is.na(presencespecies)] = 0
presencespecies_2[is.na(presencespecies_2)] = 0
forlater = presencespecies
presencespecies = presencespecies %>% filter(Cohort2_Dirty > 0 & Cohort2 == 0) 

hmdat2 = inner_join(hmdat,presencespecies %>% select(Species))
hmdat2 = hmdat2  %>% group_by(Species,DIET) %>% summarise(ABUNDANCE = mean(ABUNDANCE,na.rm=T)) %>% dcast(Species ~ DIET,value.var = 'ABUNDANCE') %>% column_to_rownames('Species')
hmdat2[is.na(hmdat2)] = 0

minval = 0.5# min(hmdat2[hmdat2>0])

hmdat2$`Coconut Oil` = log2((hmdat2$`Coconut Oil`+minval) /  (hmdat2$`Control`+minval))
hmdat2$`Ketogenic` = log2((hmdat2$`Ketogenic`+minval) /  (hmdat2$`Control`+minval))
hmdat2$`Olive Oil` = log2((hmdat2$`Olive Oil`+minval) /  (hmdat2$`Control`+minval))
hmdat2$`Palm Oil` = log2((hmdat2$`Palm Oil`+minval) /  (hmdat2$`Control`+minval))
hmdat2$`Fish Oil` = log2((hmdat2$`Fish Oil`+minval) /  (hmdat2$`Control`+minval))
hmdat2$`Lard` = log2((hmdat2$`Lard`+minval) /  (hmdat2$`Control`+minval))
hmdat2$`Milkfat` = log2((hmdat2$`Milkfat`+minval) /  (hmdat2$`Control`+minval))

hmdat2 = hmdat2 %>% select(-Control)

cleandirtydietjustdirty2 = cleandirtydietjustdirty %>% filter(BH<0.05) %>% dplyr::rename(DIET = term)
sub1 = cleandirtydietjustdirty2 %>% filter(Species %in% counts) %>% mutate(Genus = strsplit(Species,' ') %>% map_chr(1)) %>% mutate(Species = paste(Genus,'sp.')) %>% select(Species,DIET,BH) 
rawdat2 = bind_rows(cleandirtydietjustdirty2 %>% filter(!(Species %in% counts)) %>% select(Species,DIET,BH),sub1) %>% mutate(Species = gsub('s__','',Species)) %>% filter(Species!='') 
cleandirtydietjustdirty2 = bind_rows(cleandirtydietjustdirty2 %>% filter(!(Species %in% counts)) %>% select(Species,DIET,BH),sub1) %>% mutate(Species = gsub('s__','',Species)) %>% filter(Species!='') 
cleandirtydietjustdirty2=cleandirtydietjustdirty2 %>% group_by(Species,DIET) %>% slice_min(n =1, BH)

sigmat = cleandirtydietjustdirty2 %>% dcast(data = .,Species ~ DIET, value.var = 'BH') %>% filter(Species %in% rownames(hmdat2))%>% column_to_rownames('Species') 

hmdat2 = hmdat2 %>% rownames_to_column('Species') %>% filter(Species %in% rownames(sigmat)) %>% column_to_rownames('Species')

sigmat3=sigmat
sigmat3[sigmat<=0.05] = "*"
sigmat3[sigmat>0.05] = ""
sigmat3[is.na(sigmat)] = ''

sigmat3 = sigmat3[rownames(hmdat2),]
sigmat3 = t(sigmat3)

col_fun = circlize::colorRamp2(c(-2,-.00001,0,.00001,2), c("darkblue",'lightblue',"white","pink","#BE1E2D"))

pdf('~/Dropbox (Mason Lab)/mouse_diet_semir/plots/dirty_unique_differences.pdf',width=5,height=6)
Heatmap((hmdat2),col=col_fun,column_names_rot = 90, cell_fun = function(j, i, x, y, width, height, fill) { grid.text(t(sigmat3)[i, j], x, y, gp = gpar(fontsize = 6, col = "black")) })
dev.off()

uni = cleandirtydietjustdirty %>% filter(BH<0.05) %>% select(-variable) %>% dplyr::rename(variable = term) %>% select(Species,variable) %>% mutate(value = 1) %>% dplyr::rename(species = Species)

# HEATMAP 3 -- MERGED 

metaphlanregs=readRDS('~/Dropbox (Mason Lab)/mouse_diet_semir/association_output/metaphlan4_reversal_associations.rds')
diet = metaphlanregs[[1]] %>% distinct %>% filter(grepl('t__',variable))%>% mutate(SGB = strsplit(variable,'t__') %>% map_chr(2))%>% mutate(BY = p.adjust(p.value,method = 'BY'))%>% mutate(BH = p.adjust(p.value,method = 'BH')) %>% select(-BY,-BONFERRONI) %>% filter(BH<0.05) %>% filter(timepoint>0)
diet$term = gsub('DIET','',diet$term)
diet$variable = as.character(diet$variable)
tokeep = diet%>% mutate(dir = if_else(estimate>0,'Positive','Negative')) %>% filter(BH<0.05) %>% dplyr::group_by(term,dir,SGB) %>% dplyr::count() %>% filter(n>1) %>% select(-n) %>% distinct%>% ungroup %>% select(-dir) %>% group_by(term,SGB) %>% dplyr::count() %>% filter(n==1) %>% select(-n) 
cohort2ofint = inner_join(diet,tokeep) %>% group_by(SGB,term) %>% dplyr::summarise(mean = mean(estimate))  %>% mutate(dir = if_else(mean>0,'Positive','Negative'))%>% select(SGB,term,dir) %>% dplyr::rename(DIET = term)  %>% distinct

# parse dirty associations
cleandirtydietjustdirtyplot = cleandirtydietjustdirty %>% filter(BH<0.05)%>% mutate(DIET =strsplit(term,':') %>% map_chr(1)) %>% mutate(SGB = strsplit(variable,'t__') %>% map_chr(2))%>% mutate(dir = if_else(estimate>0,'Positive','Negative')) %>% select(SGB,term,dir) %>% dplyr::rename(DIET = term) %>% distinct
#cleandirtydiet1= cleandirtydiet %>% filter(!grepl('Cohort2_Dirty',term)) %>% filter(term!='(Intercept)') %>% filter(effect == 'fixed') %>% filter(BH<0.05)%>% mutate(DIET =strsplit(term,':') %>% map_chr(1)) %>% mutate(SGB = strsplit(variable,'t__') %>% map_chr(2))%>% mutate(dir = if_else(estimate>0,'Positive','Negative')) %>% select(SGB,term,dir) %>% dplyr::rename(DIET = term)

#dirtyofint = bind_rows(cleandirtydiet1,cleandirtydietjustdirtyplot) %>% distinct

# merge to find taxa of interest to plot

mutualofint = bind_rows(cleandirtydietjustdirtyplot %>% mutate(cohort = paste(DIET,'-- H+')),cohort2ofint  %>% mutate(cohort = paste(DIET,'-- OC'))) %>% mutate(sig  = "*")

# NOW LOAD IN THE CLEAN DATA

# get relevant abundances

clean_hmdat = left_join(merged_data_CLEAN   %>% mutate(SGB = strsplit(SPECIES,'t__') %>% map_chr(2)),gtdbmap)  %>% filter(SGB %in% unique(mutualofint$SGB)) %>% separate(GTDB_TAX, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = TRUE, extra = "merge") %>% group_by(Species,DIET) %>% summarise(ABUNDANCE = mean(ABUNDANCE))  %>% ungroup %>% filter(Species != 's__') %>% dcast(Species ~ DIET, value.var = 'ABUNDANCE') 

clean_and_dirt_hmdat = left_join(merged_data_cd %>% mutate(SGB = strsplit(SPECIES,'t__') %>% map_chr(2)),gtdbmap)  %>% filter(SGB %in% unique(mutualofint$SGB)) %>% separate(GTDB_TAX, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = TRUE, extra = "merge") %>% group_by(Species,DIET) %>% summarise(ABUNDANCE = mean(ABUNDANCE))  %>% ungroup %>% filter(Species != 's__') %>% dcast(Species ~ DIET, value.var = 'ABUNDANCE') 


#mutualofint_conserved = mutualofint %>% filter(schema == 'conserved')

mutualofint = left_join(mutualofint,gtdbmap) %>% filter(SGB %in% unique(mutualofint$SGB)) %>% separate(GTDB_TAX, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = TRUE, extra = "merge") %>% select(Species,cohort,DIET,sig) %>% filter(Species != 's__') %>% distinct 
mutualofint_w = mutualofint %>% dcast(Species~ cohort, value.var = 'sig')
mutualofint_w[is.na(mutualofint_w)] = ''

clean_hmdat = clean_hmdat %>% column_to_rownames('Species')
clean_and_dirt_hmdat = clean_and_dirt_hmdat %>% column_to_rownames('Species')
mutualofint_w = mutualofint_w %>% column_to_rownames('Species')

minval = 0.5 #$min(clean_and_dirt_hmdat[clean_and_dirt_hmdat>0])

clean_and_dirt_hmdat$`Coconut Oil` = log2((clean_and_dirt_hmdat$`Coconut Oil`+minval) /  (clean_and_dirt_hmdat$`Control`+minval))
clean_and_dirt_hmdat$`Ketogenic` = log2((clean_and_dirt_hmdat$`Ketogenic`+minval) /  (clean_and_dirt_hmdat$`Control`+minval))
clean_and_dirt_hmdat$`Olive Oil` = log2((clean_and_dirt_hmdat$`Olive Oil`+minval) /  (clean_and_dirt_hmdat$`Control`+minval))
clean_and_dirt_hmdat$`Palm Oil` = log2((clean_and_dirt_hmdat$`Palm Oil`+minval) /  (clean_and_dirt_hmdat$`Control`+minval))
clean_and_dirt_hmdat$`Fish Oil` = log2((clean_and_dirt_hmdat$`Fish Oil`+minval) /  (clean_and_dirt_hmdat$`Control`+minval))
clean_and_dirt_hmdat$`Lard` = log2((clean_and_dirt_hmdat$`Lard`+minval) /  (clean_and_dirt_hmdat$`Control`+minval))
clean_and_dirt_hmdat$`Milkfat` = log2((clean_and_dirt_hmdat$`Milkfat`+minval) /  (clean_and_dirt_hmdat$`Control`+minval))

clean_and_dirt_hmdat = clean_and_dirt_hmdat %>% select(-Control)
colnames(clean_and_dirt_hmdat) = paste(colnames(clean_and_dirt_hmdat),'-- H+')

clean_hmdat$`Coconut Oil` = log2((clean_hmdat$`Coconut Oil`+minval) /  (clean_hmdat$`Control`+minval))
clean_hmdat$`Ketogenic` = log2((clean_hmdat$`Ketogenic`+minval) /  (clean_hmdat$`Control`+minval))
clean_hmdat$`Olive Oil` = log2((clean_hmdat$`Olive Oil`+minval) /  (clean_hmdat$`Control`+minval))
clean_hmdat$`Palm Oil` = log2((clean_hmdat$`Palm Oil`+minval) /  (clean_hmdat$`Control`+minval))
clean_hmdat$`Fish Oil` = log2((clean_hmdat$`Fish Oil`+minval) /  (clean_hmdat$`Control`+minval))
clean_hmdat$`Lard` = log2((clean_hmdat$`Lard`+minval) /  (clean_hmdat$`Control`+minval))
clean_hmdat$`Milkfat` = log2((clean_hmdat$`Milkfat`+minval) /  (clean_hmdat$`Control`+minval))

clean_hmdat = clean_hmdat %>% select(-Control) 
colnames(clean_hmdat) = paste(colnames(clean_hmdat),'-- OC')

temp = presencespecies_2 %>% filter(Cohort2>0,Cohort2_Dirty>0) %>% select(Species) %>% unlist %>% unname %>% unique
clean_and_dirt_hmdat = clean_and_dirt_hmdat %>% data.frame(check.names=F)%>% rownames_to_column('tem2p') %>% filter((tem2p %in% temp)) %>% column_to_rownames('tem2p')
clean_hmdat = clean_hmdat %>% data.frame(check.names=F)%>% rownames_to_column('tem2p') %>% filter((tem2p %in% temp)) %>% column_to_rownames('tem2p')
mutualofint_w = mutualofint_w %>% rownames_to_column('tem2p') %>% filter((tem2p %in% temp)) %>% column_to_rownames('tem2p')

rownames(clean_hmdat) = gsub('s__','',rownames(clean_hmdat))
rownames(clean_and_dirt_hmdat) = gsub('s__','',rownames(clean_and_dirt_hmdat))
rownames(mutualofint_w) = gsub('s__','',rownames(mutualofint_w))

forhm = bind_rows(data.frame(t(clean_hmdat),check.names=F),data.frame(t(clean_and_dirt_hmdat),check.names=F))
mutualofint_w = t(mutualofint_w)

forhm = forhm[order(rownames(forhm)),]
mutualofint_w = mutualofint_w[order(rownames(mutualofint_w)),]


h1 = forhm %>% t %>% data.frame(check.names=F) %>% rownames_to_column('Species') %>% reshape2::melt() %>% mutate(cohort = strsplit(as.character(variable),' -- ') %>% map_chr(2))%>% mutate(diet = strsplit(as.character(variable),' -- ') %>% map_chr(1)) %>% mutate(value = sign(value)) %>% dcast(Species + diet~ cohort, value.var ='value') %>% filter(`H+` == `OC`) %>% select(Species) %>% table() %>% data.frame %>% filter(Freq==7)  %>% select(".") %>% unlist %>% unname

forhm1 = forhm %>% select(all_of(h1))
forhm2 = forhm %>% select(!all_of(h1))

mutualofint_w1 = mutualofint_w  %>% data.frame(check.names=F)%>% dplyr::select(all_of(h1))
mutualofint_w2 = mutualofint_w %>% data.frame(check.names=F) %>% dplyr::select(!all_of(h1))

pdf('~/Dropbox (Mason Lab)/mouse_diet_semir/plots/clean_dirty_direction_comparison_taxa_MATCHES.pdf',width=6,height=10)
Heatmap(t(forhm1),column_split = c('Coconut Oil','Coconut Oil','Fish Oil','Fish Oil','Ketogenic','Ketogenic','Lard','Lard','Milkfat','Milkfat','Olive Oil','Olive Oil','Palm Oil','Palm Oil'),cluster_columns = F,col=col_fun,column_names_rot = 60, cell_fun = function(j, i, x, y, width, height, fill) { grid.text(t(mutualofint_w1)[i, j], x, y, gp = gpar(fontsize = 10, col = "black")) })
dev.off()

pdf('~/Dropbox (Mason Lab)/mouse_diet_semir/plots/clean_dirty_direction_comparison_taxa_DIFF.pdf',width=6,height=10)
Heatmap(t(forhm2),column_split = c('Coconut Oil','Coconut Oil','Fish Oil','Fish Oil','Ketogenic','Ketogenic','Lard','Lard','Milkfat','Milkfat','Olive Oil','Olive Oil','Palm Oil','Palm Oil'),cluster_columns = F,col=col_fun,column_names_rot = 60, cell_fun = function(j, i, x, y, width, height, fill) { grid.text(t(mutualofint_w)[i, j], x, y, gp = gpar(fontsize = 10, col = "black")) })
dev.off()

### counts


b = mutualofint %>% dplyr::rename(temp=cohort) %>% mutate(cohort = strsplit(temp,' -- ') %>% map_chr(2)) %>% mutate(diet = strsplit(temp,' -- ') %>% map_chr(1)) %>% mutate(value = 1) %>% dcast(Species + diet ~ cohort,value.var = 'value') %>% mutate(`H+` = if_else(`H+`>0,1,if_else(`H+`<0,-1,0))) %>% mutate(`OC` = if_else(`OC`>0,1,if_else(`OC`<0,-1,0)))

a1 = b %>% filter(`H+`==`OC`) %>% select(Species,diet) %>% distinct %>% mutate(category = 'Conserved')
b = b %>% filter(is.na(`H+`) | is.na(OC)) %>% select(Species,diet) %>% distinct %>% mutate(category = 'Distinct')

a2 =  mutualofint %>% dplyr::rename(temp=cohort) %>% mutate(cohort = strsplit(temp,' -- ') %>% map_chr(2)) %>% mutate(diet = strsplit(temp,' -- ') %>% map_chr(1)) %>% mutate(value = 1) %>% dcast(Species + diet ~ cohort,value.var = 'value') %>% mutate(`H+` = if_else(`H+`>0,1,if_else(`H+`<0,-1,0))) %>% mutate(`OC` = if_else(`OC`>0,1,if_else(`OC`<0,-1,0)))%>% filter(`H+`==`OC`) %>% select(Species,diet) %>% distinct %>% mutate(category = 'Conserved')

c = bind_rows(a1,a2,b) %>% group_by(category,diet) %>% summarise(count = n())

uni = uni %>% group_by(variable) %>% summarise(count=n()) %>% mutate(category = 'H+ Unique') %>% dplyr::rename(diet = variable)

c = bind_rows(c,uni)

ggplot(c,aes(x=diet,y=count,fill=diet))  + cowplot::theme_cowplot() + geom_bar(stat='identity',position='dodge') + facet_grid(.~category) + theme(axis.text.x  = element_text(angle=60,hjust=1),legend.position = 'none') + scale_fill_manual(values = colorsDiet) + xlab('') + ylab('Count')
ggsave('~/Dropbox (Mason Lab)/mouse_diet_semir/plots/conserved_distinct_counts.pdf',width=7,height=4)

##### DIVERSITY ANALYSIS
# GET AND PLOT ALPHA DIVERSITIES
merged_data = readRDS('intermediate_files/metaphlan4_all_merged_data_complete.rds')
merged_data_cd = merged_data%>% filter(REVERSAL == 'NON_REVERSAL',COHORT =='Cohort2' | COHORT =='Cohort2_Dirty') %>% filter(grepl('t__',SPECIES)) 
merged_data_cd$SPECIES = as.character(merged_data_cd$SPECIES)
merged_data_cd$COHORT = as.factor(merged_data_cd$COHORT)
merged_data_cd$DIET = factor(merged_data_cd$DIET,levels=c('Control','Coconut Oil','Fish Oil','Ketogenic','Lard','Milkfat','Olive Oil','Palm Oil'))
merged_data_cd$ANIMAL_ID = as.factor(merged_data_cd$ANIMAL_ID)

merged_data_cd = merged_data_cd %>% select(ABUNDANCE,WEIGHT,SEX ,TIME_RECODE_OVERALL ,DIET,ANIMAL_ID,SPECIES,COHORT) %>% mutate(SGB = strsplit(SPECIES,'t__') %>% map_chr(2))

divs = merged_data_cd %>% group_by(COHORT) %>% select(-SPECIES,-SGB) %>% group_by(COHORT,across(-c(ABUNDANCE))) %>% summarise(shannon = vegan::diversity(ABUNDANCE, index = "shannon"), simpson = vegan::diversity(ABUNDANCE, index = "simpson"), richness = sum(ABUNDANCE > 0))

bar = divs %>% group_by(DIET,TIME_RECODE_OVERALL,SEX,COHORT)
baz = bar  %>% summarise(mean = mean(shannon),sd = sd(shannon))%>% filter(DIET == 'Control')%>% group_by(COHORT)%>% summarise(val = mean(mean))
#baz = bar%>% ungroup %>% filter(DIET == 'Control') %>% summarise(val = mean(mean))
#control_mean = baz$val[[1]]

#### THIS ONE 
ggplot(data = bar %>% filter(TIME_RECODE_OVERALL == 5 | TIME_RECODE_OVERALL == 7), aes(x = factor(COHORT), y = shannon, shape = SEX, color = factor(DIET))) + ggbeeswarm::geom_quasirandom() + facet_wrap(factor(DIET) ~ ., scales = 'fixed', ncol = 2) + cowplot::theme_cowplot() + geom_boxplot(alpha=0.2,aes(group = COHORT)) + ggtitle('Alpha diversity, clean vs. dirty cohorts') + theme(axis.text.x = element_text(hjust = 1, angle = 60)) + geom_hline(data = baz, aes(linetype = COHORT, yintercept = val), color = 'gray') + ylab('Shannon Diversity') + xlab('Months') + scale_color_manual(values = colorsDiet) + theme(legend.position = 'bottom') + ggpubr::stat_compare_means(aes(group = COHORT), method = "t.test", label = "p.format", label.x.npc = "middle", label.y.npc = "top", hjust = .1)
ggsave('~/Dropbox (Mason Lab)/mouse_diet_semir/plots/dietbaseline_diversity_cleandirt_shannon_facet_20240418.pdf',width = 5,height=8)

# GET WITHIN DIET AND TIMEPOINT BETA DIVERSITIES
# beta div 
merged1_sub = merged_data_cd%>% dcast(ANIMAL_ID + WEIGHT + SEX +  TIME_RECODE_OVERALL ~ SGB,value.var ='ABUNDANCE') %>% mutate(sid = paste(ANIMAL_ID,WEIGHT,SEX,TIME_RECODE_OVERALL,sep='-----')) %>% select(-ANIMAL_ID,-TIME_RECODE_OVERALL,-SEX,-WEIGHT) 

### REMOVE WEIRD EMPTY SAMPLE
merged1_sub = merged1_sub %>% filter(sid!='Cohort2-F-Keto_2-Ketogenic-6-----NA-----F-----2') %>% column_to_rownames('sid')

bray = as.data.frame(as.matrix(vegan::vegdist(merged1_sub))) %>% rownames_to_column('sample_id')

bray_melted = melt(bray) #%>% filter(sample_id != variable) 

mdattemp = merged_data_cd %>% select(ANIMAL_ID,WEIGHT,SEX,TIME_RECODE_OVERALL,DIET,COHORT) %>% distinct %>% mutate(sample_id =paste(ANIMAL_ID,WEIGHT,SEX,TIME_RECODE_OVERALL,sep='-----'))
braym = left_join(bray_melted,mdattemp) %>% dplyr::rename(sample_id1 = sample_id,sample_id = variable)
braym = left_join(braym,mdattemp,by= c('sample_id')) %>% dplyr::rename(sample_id2 = sample_id)

######## ACROSS MOUSE BETADIV
# filter for relevant comparisions
braym_sub = braym %>% filter(DIET.x == DIET.y,TIME_RECODE_OVERALL.x == TIME_RECODE_OVERALL.y,COHORT.x == COHORT.y,ANIMAL_ID.x != ANIMAL_ID.y)

average_values = braym_sub %>% group_by(DIET.x,TIME_RECODE_OVERALL.x,COHORT.x) %>% summarize(avg_value = mean(value, na.rm = TRUE),totalcomps = sum(value>0)) 
average_values$xaxis = paste(average_values$DIET.x,average_values$TIME_RECODE_OVERALL.x,average_values$COHORT.x)

braym_overall = left_join(braym_sub,average_values)
braym_overall$xaxis = paste(braym_overall$DIET.x,braym_overall$TIME_RECODE_OVERALL.x)

braym_overall = braym_overall %>% filter(TIME_RECODE_OVERALL.x == 5 | TIME_RECODE_OVERALL.x == 7)

# Reorder the levels of site_name_description.x based on avg_value
braym_overall$xaxis <- factor(braym_overall$xaxis,levels = average_values$xaxis[order(-average_values$avg_value)])
average_values$xaxis <- factor(average_values$xaxis,levels = average_values$xaxis[order(-average_values$avg_value)])

# Create the plot
#plot1 = ggplot(data = braym_overall, aes(x = xaxis, y = value,fill=avg_value)) + ggbeeswarm::geom_quasirandom(alpha=.3)+ geom_boxplot(alpha=.7) + cowplot::theme_cowplot() + theme(axis.text.x = element_blank(),axis.title.x = element_blank())+ scale_fill_gradient(low = "blue", high = "red") + ylab('Bray-curtis distance')+ theme(axis.title.y = element_text(angle = 0, vjust = 0.5))+ theme(panel.spacing = unit(0, "lines"), plot.margin = unit(c(0,0,0,0), "lines"))
#plot2 = ggplot(data = braym_overall %>% select(xaxis,avg_value) %>% distinct %>% reshape2::melt(id.vars = c('xaxis')) %>% mutate(group = 'a'),aes(x = xaxis,group=variable,y= value)) + geom_line(aes(group = variable))  + cowplot::theme_cowplot()+ theme(axis.text.x = element_blank(),axis.title.x = element_blank()) + ylab('Avg. within-site variation')+ theme(axis.title.y = element_text(angle = 0, vjust = 0.5))+ theme(panel.spacing = unit(0, "lines"), plot.margin = unit(c(0,0,0,0), "lines"))
#plot3 = ggplot(data = braym_overall %>% select(xaxis,TIME_RECODE_OVERALL.x) %>% distinct,aes(x = xaxis,y=1,alpha = TIME_RECODE_OVERALL.x)) + geom_tile()  + theme(axis.text.x = element_blank()) + cowplot::theme_cowplot() + ylab('Location type') + theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.title.x = element_blank())+ theme(axis.title.y = element_text(angle = 0, vjust = 0.5),legend.position = 'bottom',legend.text.align = 0,legend.title = element_blank())+ theme(panel.spacing = unit(0, "lines"), plot.margin = unit(c(0,0,0,0), "lines")) + scale_fill_brewer(palette = 'Set1')
#plot4 = ggplot(data = braym_overall %>% select(xaxis,DIET.x) %>% distinct,aes(x = xaxis,y=1,fill = DIET.x)) + geom_tile()  + theme(axis.text.x = element_blank()) + cowplot::theme_cowplot() + ylab('Location type') + theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.title.x = element_blank())+ theme(axis.title.y = element_text(angle = 0, vjust = 0.5),legend.position = 'bottom',legend.text.align = 0,legend.title = element_blank())+ theme(panel.spacing = unit(0, "lines"), plot.margin = unit(c(0,0,0,0), "lines"))+ scale_fill_brewer(palette = 'Set1')

#combined_plot <- plot1 / (plot2 / plot3 / plot4) + plot_layout(heights = c(5, 2.3))

#ggsave('~/Dropbox (Mason Lab)/mouse_diet_semir/plots/diet_baseline_across_mouse_betadiv_DIRTY.pdf',width=10,height=8)

##### ACROSS MOUSE BETADIV VISUALIZE DIFFERENTLY
#ggplot(data = braym_overall, aes(x = factor(TIME_RECODE_OVERALL.x), y = value, fill=avg_value)) + ggbeeswarm::geom_quasirandom(alpha=.3) + geom_boxplot(alpha=.7) + cowplot::theme_cowplot() + facet_grid(. ~ DIET.x) + scale_fill_gradient(low = "blue", high = "red") + ylab('Bray-curtis distance') + xlab('Timepoint') + theme(axis.title.y = element_text(angle = 90, vjust = 0.5))
#ggsave('~/Dropbox (Mason Lab)/mouse_diet_semir/plots/diet_baseline_across_mouse_betadiv_facetplot_DIRTY.pdf',width=8,height=4)


######## WITHIN MOUSE BETADIV
# filter for relevant comparisions
braym_sub =braym %>% filter(DIET.x == DIET.y,COHORT.x == COHORT.y,ANIMAL_ID.x == ANIMAL_ID.y,COHORT.x == COHORT.y)%>% filter(value!=0)

average_values = braym_sub %>% group_by(ANIMAL_ID.x,DIET.x,COHORT.x) %>% summarize(avg_value = mean(value, na.rm = TRUE),totalcomps = sum(value>0)) 
average_values$xaxis = paste(average_values$ANIMAL_ID.x,average_values$COHORT.x)

braym_overall = left_join(braym_sub,average_values)
braym_overall$xaxis = paste(braym_overall$ANIMAL_ID.x)

# Reorder the levels of site_name_description.x based on avg_value
braym_overall$xaxis <- factor(braym_overall$xaxis,levels = average_values$xaxis[order(-average_values$avg_value)])
average_values$xaxis <- factor(average_values$xaxis,levels = average_values$xaxis[order(-average_values$avg_value)])

# Create the plot
#plot1 = ggplot(data = braym_overall, aes(x = xaxis, y = value,alpha=(TIME_RECODE_OVERALL.x/12),fill=avg_value)) + ggbeeswarm::geom_quasirandom()+ geom_boxplot(alpha=.7) + cowplot::theme_cowplot() + theme(axis.text.x = element_blank(),axis.title.x = element_blank())+ scale_fill_gradient(low = "blue", high = "red") + ylab('Bray-curtis distance')+ theme(axis.title.y = element_text(angle = 0, vjust = 0.5))+ theme(panel.spacing = unit(0, "lines"), plot.margin = unit(c(0,0,0,0), "lines"))
#plot2 = ggplot(data = braym_overall %>% select(xaxis,avg_value) %>% distinct %>% reshape2::melt(id.vars = c('xaxis')) %>% mutate(group = 'a'),aes(x = xaxis,group=variable,y= value)) + geom_line(aes(group = variable))  + cowplot::theme_cowplot()+ theme(axis.text.x = element_blank(),axis.title.x = element_blank()) + ylab('Avg. within-site variation')+ theme(axis.title.y = element_text(angle = 0, vjust = 0.5))+ theme(panel.spacing = unit(0, "lines"), plot.margin = unit(c(0,0,0,0), "lines"))
#plot3 = ggplot(data = braym_overall %>% select(xaxis,TIME_RECODE_OVERALL.x) %>% distinct,aes(x = xaxis,y=1,alpha = TIME_RECODE_OVERALL.x)) + geom_tile()  + theme(axis.text.x = element_blank()) + cowplot::theme_cowplot() + ylab('Location type') + theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.title.x = element_blank())+ theme(axis.title.y = element_text(angle = 0, vjust = 0.5),legend.position = 'bottom',legend.text.align = 0,legend.title = element_blank())+ theme(panel.spacing = unit(0, "lines"), plot.margin = unit(c(0,0,0,0), "lines")) + scale_fill_brewer(palette = 'Set1')
#plot4 = ggplot(data = braym_overall %>% select(xaxis,DIET.x) %>% distinct,aes(x = xaxis,y=1,fill = DIET.x)) + geom_tile()  + theme(axis.text.x = element_blank()) + cowplot::theme_cowplot() + ylab('Location type') + theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.title.x = element_blank())+ theme(axis.title.y = element_text(angle = 0, vjust = 0.5),legend.position = 'bottom',legend.text.align = 0,legend.title = element_blank())+ theme(panel.spacing = unit(0, "lines"), plot.margin = unit(c(0,0,0,0), "lines"))+ scale_fill_brewer(palette = 'Set1')

#combined_plot <- plot1 / (plot2 / plot4) + plot_layout(heights = c(5, 2.3))

#ggsave('~/Dropbox (Mason Lab)/mouse_diet_semir/plots/diet_baseline_within_mouse_betadiv_dirty.pdf',width=10,height=8)

### alt vis
average_values = average_values %>% group_by(DIET.x) %>% mutate(med2 = median(avg_value))

ggplot(data = average_values, aes(y = fct_reorder(DIET.x,med2),color= COHORT.x, x = avg_value))+ ggbeeswarm::geom_quasirandom(size=2) + geom_boxplot(alpha=.7) + cowplot::theme_cowplot() + scale_fill_gradient(low = "blue", high = "red") + xlab('Bray-curtis distance') + ylab('Diet') + theme(axis.title.y = element_text(angle = 90, vjust = 0.5)) + ggtitle('Within-Mouse Variation')
ggsave('~/Dropbox (Mason Lab)/mouse_diet_semir/plots/diet_baseline_within_mouse_betadiv_singleplot_dirty.pdf',width=4,height=4)

### another alt vis

#### THIS ONE --- NEED TO GET NON DIRTY DATA HERE OTO
##### ACROSS MOUSE BETADIV VISUALIZE DIFFERENTLY
bar = braym_overall %>% group_by(DIET.x,SEX.x,TIME_RECODE_OVERALL.x,COHORT.x) #%>% summarise(mean = mean(value),sd = sd(value))
baz = bar %>% summarise(mean = mean(value),sd = sd(value))%>% filter(DIET.x == 'Control')%>% group_by(COHORT.x)%>% summarise(val = mean(mean))

#### THIS ONE
ggplot(data = bar %>% filter(TIME_RECODE_OVERALL.x == 5 | TIME_RECODE_OVERALL.x == 7,TIME_RECODE_OVERALL.y == 5 | TIME_RECODE_OVERALL.y == 7),aes(x = factor(COHORT.x),y = value,color=factor(DIET.x)))+geom_point() + facet_wrap(factor(DIET.x) ~ SEX.x,scales='fixed',nrow=2)  + ggbeeswarm::geom_quasirandom() + geom_boxplot(alpha=0.5)+ ggtitle('Beta diversity, clean vs. dirty cohorts')  + geom_hline(data = baz,aes(linetype = COHORT.x,yintercept = val),color='gray') + ylab('Beta Diversity') + xlab('Months')+ scale_color_manual(values=colorsDiet) + theme(legend.position = 'bottom')+ cowplot::theme_cowplot()+ theme(axis.text.x = element_text(hjust=1,angle=60))


#+geom_errorbar(aes(ymin= mean-sd,ymax = mean+sd),width=.1) + geom_line(aes(linetype = COHORT.x,group= COHORT.x))+ cowplot::theme_cowplot()  + ggtitle('Beta diversity, clean vs. dirty cohorts')+ theme(axis.text.x = element_text(hjust=1))  + geom_hline(data = baz,aes(linetype = COHORT.x,yintercept = val),color='gray') + ylab('Beta Diversity') + xlab('Months')+ scale_color_manual(values=colorsDiet) + theme(legend.position = 'bottom')

ggsave('~/Dropbox (Mason Lab)/mouse_diet_semir/plots/dietbeta_diversity_cleandirt_facet_20240418.pdf',width = 10,height=4)

##### CEM VOLCANO PLOT THIS ONE

res = cleandirtyoverall2

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
res$Species = gsub('s__','',res$Species)
res$Metabolite = res$Species
res$Significance = "Insignificant"
res$Significance[ which(res$BH < qCutoff & res$estimate > estimateCutoff)] = "Increased"
res$Significance[ which(res$BH < qCutoff & res$estimate < -estimateCutoff)] = "Decreased"

top_microbes$Species = gsub('s__','',top_microbes$Species)
top_microbes$Metabolite = top_microbes$Species
top_microbes$Significance = "Insignificant"
top_microbes$Significance[ which(top_microbes$BH < qCutoff & top_microbes$estimate > estimateCutoff)] = "Increased"
top_microbes$Significance[ which(top_microbes$BH < qCutoff & top_microbes$estimate < -estimateCutoff)] = "Decreased"

nudgeMax = 3 #max(-log10(resLabel$padj), na.rm=T) + 1
nudger = function(x, lim) { return( ifelse( (lim-x) > 0, (lim-x), 0) ) }
ggplot(res, aes(x=estimate, y=-log10(BH), fill=Significance)) + geom_hline(yintercept=0) + 
  geom_hline(yintercept=-log10(qCutoff), linetype="dashed") + 
  geom_vline(xintercept=c(-estimateCutoff, estimateCutoff), linetype="dashed") +  
  geom_point(aes(size=Significance, shape=Significance)) + 
  geom_label_repel(data=top_microbes, aes(label=Species), nudge_y=nudger(-log10(top_microbes$BH), nudgeMax), segment.size = 0.5, box.alpha=0.5,segment.alpha=0.6, size=6, box.padding = 0.1, label.padding = 0.1, label.r = 0.05,  max.overlaps = 1000000) + 
  scale_size_manual(values=dirSize) + 
  scale_fill_manual(name="Species Direction", values=dirColors) + 
  scale_shape_manual(values=dirShape) + 
  scale_x_continuous(name="Beta Coefficient") + 
  scale_y_continuous("-log10( adjusted p-value)", oob=scales::squish) + 
  cowplot::theme_cowplot() + theme(legend.position="bottom") +
  ggtitle(paste0("Volcano plot\n-log10(qval) > ", minSignificance, " are squished"))

ggsave("~/Dropbox (Mason Lab)/mouse_diet_semir/plots/cleandirty_volcano_cem.pdf", width=14, height=8)


