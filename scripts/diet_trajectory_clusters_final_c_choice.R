library(ggplot2)
library(data.table)
library(matrixStats)
library(plyr)
library(reshape2)
library(dplyr)
library(TSP)
library(ape)
library(ggtree)
library(phylobase)
library(tidyverse)
library(furrr)
library(broom)
library(progressr)
library(dynamicTreeCut)
library(ComplexUpset)
library(ComplexHeatmap)
library(cowplot)

membcut=0.2

##### FIGURE OUT REMAP OF DEFINING INCLUSION
##### ADD IN % OF SPECIES DETECTED/REVERTED/PERSISTENT TO THIS FIGURE AND F1


plotSpeciesnonrev = function(datLfcTime, speciesList)
{
  speciesNames =gsub('\\|t','',splitGetFromEnd(speciesList, "__", 2))
  speciesNames = gsub("_", " ", speciesNames)
  speciesNames = paste0(speciesNames, collapse=", ")
  
  curData = datLfcTime[ datLfcTime$SPECIES %in% speciesList, ]
  curDataM = reshape2::melt(curData, c("SPECIES", "Rank", "Diet"), variable.name="Type", value.name="log2FC")
  curDataM$group = paste0(curDataM$SPECIES, "_", splitGet(as.character(curDataM$Type), "_", 1))
  curDataM$Linetype = "1solid"
  
  curDataM2 = rbind(curDataM)
  curDataM2$group = paste0(curDataM2$cluster, "_", splitGet(as.character(curDataM2$Type), "_", 1))
  curDataM2 = curDataM2[ order(curDataM2$Type), ]
  curDataM2$Type = gsub("Pre", "", curDataM2$Type)
  
  curDataM2$Type = factor(curDataM2$Type,levels = c('D0','M1','M4','M6','M9','M11'))
  
 print(ggplot(curDataM2, aes(x = as.factor(Type), y = log2FC)) + 
    geom_line(aes(colour = Diet, group = Diet), size=1.1) +
    scale_colour_manual(values=colorsDiet) +
    scale_x_discrete(name="Timepoint") +
    ylab("log2 fold change") +  
    geom_hline(yintercept=0) + 
    ggtitle(speciesNames) +
    theme_cem + theme(legend.position="none"))
  ggsave(paste0('~/Dropbox (Mason Lab)/mouse_diet_semir/plots/',speciesList[[1]],'nonrev.pdf'),width=5,height=5)
}

### MERGE IN GTDB SGB MAPPING
gtdbmap = read.delim('~/Dropbox (Mason Lab)/mouse_diet_semir/mpa_vJan21_CHOCOPhlAnSGB_202103_SGB2GTDB.tsv',sep='\t',header=F)
colnames(gtdbmap) = c('SGB','GTDB_TAX')

# todo here -- genus-level tree, highlight specific clusters of interest, filter for membership in the final plot

metatree = read.tree('~/Dropbox (Mason Lab)/mouse_diet_semir/mpa_vJan21_CHOCOPhlAnSGB_202103.nwk')
mapfile = gtdbmap %>% separate(GTDB_TAX, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = TRUE, extra = "merge")  %>% filter(Family!= 'f__') 
mapfile$SGB = gsub('SGB','',mapfile$SGB)
mapfile$SGB = gsub('_group','',mapfile$SGB)
metatree = keep.tip(tip = intersect(mapfile$SGB,metatree$tip.label),phy = metatree)
mapfile= mapfile %>% column_to_rownames('SGB')
mapfile = mapfile[metatree$tip.label,,drop=T]
metatree$tip.label = mapfile$Species

## Functions
{
  theme_linedraw2 = theme_linedraw() + theme(strip.background=element_rect(fill="grey80", colour="grey50", size=0.2), strip.text.x=element_text(colour="black"), strip.text.y=element_text(colour="black"))
  
  theme_cem = theme_linedraw2 + theme(legend.position="right", legend.title=element_text(size=15), legend.text=element_text(size=14), axis.text.x = element_text(size=12, angle=-90, hjust=0, vjust=0.5), axis.text.y=element_text(size=12), axis.title=element_text(size=15), axis.title.y=element_text(vjust=1), plot.title = element_text(size=18, vjust=1.5), strip.background = element_rect(fill="#EEEEEE"), strip.text = element_text(size = 11), panel.grid.major = element_line(colour = "grey98"), panel.grid.minor = element_blank())
  
  getListIndex = function(x, n)
  {
    if(n > length(x)) return("") else return(x[[n]])
  }
  splitGet=function (strings, sep, n)
  {
    splitVals = strsplit(as.character(strings), sep)
    return ( sapply(splitVals, getListIndex, n) )
  }
  
  splitGetFromEnd=function (strings, sep, n)
  {
    splitVals = strsplit(as.character(strings), sep)
    strLength = sapply(splitVals,length) - n + 1
    return(mapply("[[", splitVals, strLength))
  }
  
  write.tsv = function(object, file, row.names=T, quote=F, col.names=NA)
  {
    if(row.names==F & is.na(col.names)) col.names=T
    write.table(object, file, sep="\t", quote=quote, row.names=row.names, col.names=col.names)
  }
  
  PairwiseCompare = function(matDat, samplesCase, samplesCtrl, colSuffix ="", operation=c("logratio", "subtract", "divide"), summaryFun=c("median", "mean", "sum"))
  {
    allCols = paste0(c(samplesCase, samplesCtrl), colSuffix)
    missingCol = setdiff(allCols, colnames(matDat))
    if(length(missingCol) > 0)
    {
      stop(paste0("Missing columns: ", paste0(missingCol, collapse=", ")))
    }
    matDat = as.data.frame(matDat)
    tempDat = data.frame(row.names=rownames(matDat))
    for(curSample1 in samplesCase)
    {
      for(curSample2 in samplesCtrl)
      {
        curCol1 = paste0(curSample1, colSuffix)
        curCol2 = paste0(curSample2, colSuffix)
        if(operation == "subtract")
        {
          curColDat = matDat[, curCol1] - matDat[, curCol2]
        }else if(operation == "divide")
        {
          curColDat = matDat[, curCol1] / matDat[, curCol2]
        }else if(operation == "logratio")
        {
          curColDat = log2(matDat[, curCol1] / matDat[, curCol2])
        }
        tempDat = cbind(tempDat, curColDat)
      }
    }
    res = NULL
    if(summaryFun == "mean")
    {
      res = rowMeans(as.matrix(tempDat))
    }else if(summaryFun == "median")
    {
      res = rowMedians(as.matrix(tempDat))
    }else if(summaryFun == "mean")
    {
      res = rowSums(as.matrix(tempDat))
    }
    return(res)
  }
  
  #helper function for the within sum of squared error
  sumsqr = function(x, clusters)
  {
    sumsqr = function(x) sum(scale(x, scale = FALSE)^2)
    wss = sapply(split(as.data.frame(x), clusters), sumsqr)
    return(wss)
  }
  
  #get the wss for repeated clustering
  iterate_fcm_WSS = function(df,m)
  {
    totss = numeric()
    for (i in 2:20)
    {
      FCMresults = cmeans(df, centers = i, m = m)
      totss[i] = sum(sumsqr(df, FCMresults$cluster))
    }
    return(totss)
  }
  
  mestimate = function(df)
  {
    N = dim(df)[[1]]
    D = dim(df)[[2]]
    m.sj = 1 + (1418/N + 22.05)*D^(-2) + (12.33/N +0.243)*D^(-0.0406*log(N) - 0.1134)
    return(m.sj)
  }
  
  getTSPOrder = function(dat)
  {
    library(TSP)
    curDist = dist(dat)
    curDist = curDist / max(curDist)
    curTsp = TSP(curDist)
    curTsp = insert_dummy(curTsp, label = "cut")
    #curTour = solve_TSP(curTsp, method="concorde", exe="~/.local/bin/concorde.exe")
    curTour = solve_TSP(curTsp, method="nn")
    curTour = cut_tour(curTour, "cut")
    curOrder = labels(curTour)
    return(curOrder)
  }
  
}

colorsDiet = c(
  Control = "#377eb8",
  Lard = "#e41a1c",
  Milkfat = "#800026",
  Ketogenic = "#ff7f00",
  FishOil = "#f781bf",
  CoconutOil = "#a65628",
  OliveOil = "#984ea3",
  PalmOil = "#4daf4a"
)

colorsDiet2 = c(
  Control = "#377eb8",
  Lard = "#e41a1c",
  Milkfat = "#800026",
  Ketogenic = "#ff7f00",
  `Fish Oil` = "#f781bf",
  `Coconut Oil` = "#a65628",
  `Olive Oil` = "#984ea3",
  `Palm Oil` = "#4daf4a"
)


## Clustering/Plotting Non-reversal of interest
{
  outName = "Trajectory_nonreversal_species_pseudo0.1"
  comparisonOrder = c("D0", "M1", "M4", "M6", "M9", "M11")
  
  setwd('~/Dropbox (Mason Lab)/mouse_diet_semir/cluster_analysis_may24/')
  
  merged_data = readRDS('~/Dropbox (Mason Lab)/mouse_diet_semir/intermediate_files//metaphlan4_all_merged_data_complete.rds')
  merged_data_cd = merged_data%>% filter(REVERSAL == 'NON_REVERSAL', COHORT == 'Cohort2') %>% filter(grepl('t__',SPECIES)) %>% mutate(SGB = strsplit(SPECIES,'t__') %>% map_chr(2)) # %>% filter(SPECIES %in% bugsofinterest)
  
  metaphlanregs=readRDS('~/Dropbox (Mason Lab)/mouse_diet_semir/association_output/metaphlan4_reversal_associations.rds')
  diet = metaphlanregs[[1]] %>% filter(grepl('t__',variable))%>% mutate(SGB = strsplit(variable,'t__') %>% map_chr(2))%>% mutate(BY = p.adjust(p.value,method = 'BY'))%>% mutate(BH = p.adjust(p.value,method = 'BH')) %>% select(-BY,-BONFERRONI)
  tokeep = diet %>% mutate(dir = if_else(estimate>0,'Positive','Negative')) %>% filter(BH<0.05) %>% dplyr::group_by(term,dir,SGB) %>% dplyr::count() %>% filter(n>2) %>% select(-n) %>% distinct%>% ungroup %>% select(-dir) %>% group_by(term,SGB) %>% dplyr::count() %>% filter(n==1) %>% select(-n) 
  diet = inner_join(diet,tokeep)
  diet$term = gsub('DIET','',diet$term)
  diet$variable = as.character(diet$variable)
  tokeep  = diet %>% filter(dietstatus == 'ON DIET',BH<0.05) %>% select(SGB) %>% distinct %>% unlist %>% unname
  merged_data_cd =merged_data_cd %>% filter(SGB %in% tokeep)
  
  merged_data_cd$SPECIES = as.character(merged_data_cd$SPECIES)
  merged_data_cd$COHORT = as.factor(merged_data_cd$COHORT)
  merged_data_cd$DIET = factor(merged_data_cd$DIET, levels = c('Control','Coconut Oil','Fish Oil','Ketogenic','Lard','Milkfat','Olive Oil','Palm Oil'))
  merged_data_cd$ANIMAL_ID = as.factor(merged_data_cd$ANIMAL_ID)
  
  merged_data_cd$DietTime = paste0( gsub(" ", "", as.character(merged_data_cd$DIET)), "_", merged_data_cd$TIMEPOINT)
  merged_data_cd$Rank = splitGet( splitGetFromEnd(merged_data_cd$SPECIES, "\\|", 1), "__", 1)
  #merged_data_cd$Rank = plyr::mapvalues(merged_data_cd$Rank, c("k", "p", "c", "o", "f", "g", "s"), c("1K", "2P", "3C", "4O", "5F", "6G", "8T"))
  merged_data_cd$Rank = plyr::mapvalues(merged_data_cd$Rank, c("t"), c("8T"))
  merged_data_cd$TIMEPOINT = factor(merged_data_cd$TIMEPOINT, levels = c("D0", "M1", "M4", "M6", "M9", "M11"))
  metadata = unique(merged_data_cd[, c("UID", "TIMEPOINT", "DIET", "SEX", "BATCH", "COHORT", "DietTime")])
  metadata = metadata[ order(metadata$DIET, metadata$TIMEPOINT), ]
  
  pseudocount = 0.5
  
  ## Calculate log2FC based on median of log2FC values across all replicate pairs
  datWide = reshape2::dcast(merged_data_cd, SPECIES+Rank~UID, value.var="ABUNDANCE")
  datWide[, 3:ncol(datWide)] = datWide[, 3:ncol(datWide)] + pseudocount
  datLfc = datWide[, c(1,2)]
  for(curDiet in unique(metadata$DIET))
  {
    if(curDiet == "Control") next
    for(curTime in unique(metadata$TIMEPOINT))
    {
      curCol = paste0( gsub(" ", "", as.character(curDiet)), "_", curTime)
      samplesCase = metadata$UID[ metadata$DIET %in% curDiet & metadata$TIMEPOINT %in% curTime]
      samplesCtrl = metadata$UID[ metadata$DIET %in% "Control" & metadata$TIMEPOINT %in% curTime]
      
      curLfc = data.frame(row.names = datWide[, 1])
      curLfc[, curCol] = PairwiseCompare(datWide, samplesCase, samplesCtrl, operation="logratio", summaryFun="median" )
      datLfc = cbind(datLfc, curLfc)
    }
  }
  
  
  datLfcLong = reshape2::melt(datLfc, c("SPECIES", "Rank"), variable.name="DietTime", value.name="Log2FC")
  datLfcLong$Diet = splitGet(datLfcLong$DietTime, "_", 1)
  datLfcLong$Time = splitGet(datLfcLong$DietTime, "_", 2)
  
  datLfcTime = reshape2::dcast(datLfcLong, SPECIES+Rank+Diet~Time, value.var="Log2FC")
  write.tsv(datLfcTime, paste0(outName, "_log2FC_wide.txt"), row.names = F)
  
  datLfcTime = datLfcTime[ datLfcTime$Rank != "1K", ]
  
  ## ?? Filter to species or not ??
  datLfcTime = datLfcTime[ datLfcTime$Rank == "8T", ] 
  
  rownames(datLfcTime) = paste0(datLfcTime$SPECIES, ":", datLfcTime$Diet)
  datLfcTime2 = datLfcTime[, -c(1:3)]
  datLfcTime2 = datLfcTime2[ rowSds(as.matrix(datLfcTime2)) > 0, ] ## remove rows with all zeroes
  #####
  
  
  library(e1071)
  library(tidyr)
  library(dplyr)
  library(pheatmap)
  
  dataList = list(LFC_sig=datLfcTime2)
  kList = c(20)
  distMethod = "euclidean"
  for(k in kList)
  {
    for(i in 1:length(dataList))
    {
      curData = dataList[[i]]
      curName = names(dataList)[i]
      curOut = paste0(outName, "_cmeans_", curName, "_k", k, "_" , distMethod)
      curData = curData[ complete.cases(curData), ]
      
      ## Normalize/standardize data by some metric?
      curData2 = curData
      #curData2 = t(scale(t(curData)))
      #curData = curData - rep.col(rowMeans(curData), ncol(curData))
      
      ## M estimate
      m = mestimate(curData2)
      
      if( ! file.exists(paste0(curOut, ".rds")) )
      {
        fcm_results = cmeans(curData2, centers = k, m = m, iter.max = 2000, dist=)
        center_order = getTSPOrder(fcm_results$centers)
        fcm_results$centers = fcm_results$centers[ center_order, ]
        rownames(fcm_results$centers) = paste0("C_", sprintf("%02d", 1:k))
        fcm_results$membership = fcm_results$membership[, center_order]
        colnames(fcm_results$membership) = paste0("C_", sprintf("%02d", 1:k))
        
        fcm_results$cluster = mapvalues(fcm_results$cluster, center_order, paste0("C_", sprintf("%02d", 1:k)))
        fcm_results$size = fcm_results$size[as.numeric(center_order)]
        saveRDS(fcm_results, paste0(curOut, ".rds"))
      }else
      {	
        fcm_results = readRDS(paste0(curOut, ".rds"))
      }
      
      
      fcm_centroids = fcm_results$centers
      fcm_centroids_df = data.frame(fcm_centroids)
      fcm_centroids_df$cluster = row.names(fcm_centroids_df)
      centroids_long = tidyr::gather(fcm_centroids_df, "Type", "log2FC", 1:ncol(fcm_centroids))
      centroids_long$Type = factor(centroids_long$Type, levels = comparisonOrder)
      
      #start with the input data # tmploc
      fcm_plotting_df = data.frame(curData)
      
      #add genes
      fcm_plotting_df$gene = row.names(fcm_plotting_df)
      
      #bind cluster assinment
      fcm_plotting_df$cluster = fcm_results$cluster
      #fetch the membership for each gene/top scoring cluster
      fcm_plotting_df$membership = sapply(1:length(fcm_plotting_df$cluster), function(row)
      {
        clust = fcm_plotting_df$cluster[row]
        fcm_results$membership[row, clust]
      })
      
      
      cluster_plot_df = fcm_plotting_df %>% group_by(cluster) %>%
        dplyr::select(., 1:ncol(curData), membership, gene, cluster) %>%
        tidyr::gather(., "Type", "log2FC", 1:ncol(curData))
      
      cluster_plot_df = as.data.frame(cluster_plot_df)
      cluster_plot_df$Type = factor(cluster_plot_df$Type, levels = comparisonOrder)
      
      #order the dataframe by score
      cluster_plot_df = cluster_plot_df[order(cluster_plot_df$membership),]
      #set the order by setting the factors using forcats
      cluster_plot_df$gene = forcats::fct_inorder(cluster_plot_df$gene)
      
      cluster_plot_df$group = paste0(cluster_plot_df$gene)
      cluster_plot_df$Type = factor(cluster_plot_df$Type, levels = comparisonOrder)
      cluster_plot_df$Diet = splitGet(cluster_plot_df$gene, ":", 2)
      cluster_plot_df = cluster_plot_df[ order(cluster_plot_df$Type), ]
      
      write.tsv(cluster_plot_df, paste0(curOut, "_long.txt"), row.names=F)
      
      #### filtering for membereship
      cluster_plot_df = cluster_plot_df %>% filter(membership>membcut)
      
      ## cluster centroids as trajectory lines
      core = centroids_long
      core$group = paste0(core$cluster)
      
      #####
     # clusterstokeep =c("C_01", "C_02", "C_03", "C_04", "C_05", "C_06", "C_07", "C_08", "C_09", "C_10", "C_11", "C_12", "C_13", "C_14", "C_15", "C_16", "C_17", "C_18", "C_19", "C_20")
      clusterstokeep = c("C_08","C_09","C_12","C_10","C_13","C_19","C_20","C_16","C_18","C_15","C_17","C_06","C_14","C_05")
      cluster_plot_df = cluster_plot_df %>% filter(cluster %in% clusterstokeep)
      cluster_plot_df$cluster = factor(cluster_plot_df$cluster,levels = clusterstokeep) 
      
      core = core %>% filter(cluster %in% clusterstokeep)
      core$cluster = factor(core$cluster,levels = clusterstokeep) 
      
      ggplot(cluster_plot_df, aes(x = Type, y = log2FC)) + 
        geom_line(aes(colour = membership, group = group)) +
        scale_colour_gradientn(colours=c('yellow1','red2')) +
        geom_line(data=core, aes(Type, log2FC, group = group), color="black",inherit.aes=FALSE, linewidth=0.8) +
        scale_x_discrete(name="Timepoint", limits=comparisonOrder) +
        ylab("log2 fold change vs time-matched control diet") +
        labs(title= paste0("All clusters log2FC vs time-matched Control Diet"), color = "Score") + facet_wrap(~cluster, scales="free",nrow=2) + theme_cem + 
        geom_hline(yintercept=0) + theme(legend.position = 'none')
      ggsave(paste0(curOut, "_clusterProfile_FORFIG.pdf"), width=12, height=4)
    
      
      ### NOW TO GENERATE SOME OTHER INFO
      cluster_plot_df = cluster_plot_df %>% mutate(SGB = strsplit(group,':') %>% map_chr(1) %>% strsplit('t__') %>% map_chr(2))
      cluster_plot_df = left_join(cluster_plot_df,gtdbmap)  %>% separate(GTDB_TAX, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = TRUE, extra = "merge")
      # get total number of lines in each cluster
      clustersizes = cluster_plot_df %>% select(cluster,group) %>% distinct %>% group_by(cluster) %>% select(-group)%>% table %>% data.frame %>% dplyr::rename(cluster = ".",cluster_counts = Freq)
      # find the number of unique species per cluster
      speciescounts = cluster_plot_df %>% select(-Type,-log2FC) %>% distinct %>% select(cluster,Species) %>% table %>% data.frame %>% filter(Freq>0)
      colnames(speciescounts) = c('cluster','Species','Species_Count')
      # fund the number of unique diets per cluster
      dietcounts = cluster_plot_df %>% select(-Type,-log2FC) %>% distinct %>% select(cluster,Diet) %>% table %>% data.frame %>% filter(Freq>0)
      colnames(dietcounts) = c('cluster','Diet','Diet_Count')
      clustersumsspecies = left_join(clustersizes,speciescounts) %>% mutate(species_frac = Species_Count/cluster_counts) %>% filter(Species != 's__')
      clustersumsdiets = left_join(clustersizes,dietcounts) %>% mutate(diet_frac = Diet_Count/cluster_counts)
      # generate tree heatmap thing -- nonreversal
        #wide = clustersumsspecies %>% filter(Species_Count>1) %>% group_by(cluster) %>% slice_max(n=5,species_frac) %>% reshape2::dcast(cluster ~ Species,value.var = 'species_frac') %>% column_to_rownames('cluster')
        wide = clustersumsspecies %>% group_by(cluster) %>% reshape2::dcast(cluster ~ Species,value.var = 'species_frac') %>% column_to_rownames('cluster')
        wide[is.na(wide)]= 0
        metatreesub = keep.tip(tip = colnames(wide),phy = metatree)
        wide = wide %>% t() %>% data.frame(check.names = F)#%>% rownames_to_column('id')
        mapfilesub = mapfile %>% filter(Species %in% rownames(wide)) %>% select(Species,Phylum,Class,Order,Family,Genus) 
        rownames(mapfilesub) =NULL
        mapfilesub = mapfilesub %>% column_to_rownames('Species')
        mapfilesub[metatreesub$tip.label,] = mapfilesub
        wide[metatreesub$tip.label,] = wide
        mapfilesub = mapfilesub %>% rownames_to_column('Species')
        tree = ggtree(metatreesub,layout = 'dendrogram') 
        tree = tree %<+% mapfilesub
        tree = tree + geom_tiplab(align = TRUE,aes(color=Phylum)) + geom_treescale()
        tree1 = gheatmap(tree, wide, offset = .5, width = 1, colnames_angle = 0, hjust = 1)  + scale_fill_gradientn(colors = c('white', "red", "darkblue"), values = scales::rescale(c(0, .5, 1)), limits = c(0, 1), guide = "colourbar")
        ggsave(plot = tree1,path.expand(paste0("~/Dropbox (Mason Lab)/mouse_diet_semir/cluster_analysis_may24/individual_plots/nonreversal/",as.character(k),'/clustertree.pdf')),width=14,height=10)
        # now make the tree at the genus level
        clustersumsspeciesgen = clustersumsspecies %>% mutate(genus = strsplit(as.character(Species),' ') %>% map_chr(1) %>% gsub('s__','',.)) %>% group_by(genus,cluster) %>% mutate(Genus_Count = sum(Species_Count),genus_frac = Genus_Count/cluster_counts) %>% ungroup %>% group_by(cluster,genus) %>% sample_n(1) ##### FIXING
        metatreesub = keep.tip(tip = as.character(clustersumsspeciesgen$Species),phy = metatree)
        wide = clustersumsspeciesgen%>% group_by(cluster) %>% reshape2::dcast(cluster ~ genus,value.var = 'species_frac') %>% column_to_rownames('cluster')
        wide[is.na(wide)]= 0
        metatreesub$tip.label = strsplit(metatreesub$tip.label,' ') %>% map_chr(1) %>% gsub('s__','',.)
        metatreesub = keep.tip(tip = colnames(wide),phy = metatreesub)
        wide = wide %>% t() %>% data.frame(check.names = F)#%>% rownames_to_column('id')
        mapfilesub = mapfile %>% mutate(Genus = gsub('g__','',Genus))%>% filter(Genus %in% rownames(wide)) %>% select(Genus,Phylum,Class,Order,Family) %>% distinct 
        rownames(mapfilesub) =NULL
        mapfilesub = mapfilesub %>% column_to_rownames('Genus')
        mapfilesub[metatreesub$tip.label,] = mapfilesub
        wide[metatreesub$tip.label,] = wide
        mapfilesub = mapfilesub %>% rownames_to_column('Genus')
        tree = ggtree(metatreesub,layout = 'dendrogram') 
        tree = tree %<+% mapfilesub
        tree = tree + geom_tiplab(align = TRUE,aes(color=Phylum)) + geom_treescale()
        tree1 = gheatmap(tree, wide, offset = .5, width = 1, colnames_angle = 0, hjust = 1)  + scale_fill_gradientn(colors = c('white', "red", "darkblue"), values = scales::rescale(c(0, .5, 1)), limits = c(0, 1), guide = "colourbar")
        ggsave(plot = tree1,path.expand(paste0("~/Dropbox (Mason Lab)/mouse_diet_semir/cluster_analysis_may24/individual_plots/nonreversal/",as.character(k),'/clustertree_GENUS.pdf')),width=14,height=10)
    }
    # make an associated stacked barplot
    longplot = ggplot(clustersumsdiets,aes(y = factor(cluster,levels=rev(clusterstokeep)), x = diet_frac*100, fill = Diet)) + geom_bar(stat = 'identity') + theme_minimal() + scale_fill_manual(values=colorsDiet) + theme(legend.position = 'none') + xlab('%') + ylab('')
    ggsave(plot = longplot,path.expand(paste0("~/Dropbox (Mason Lab)/mouse_diet_semir/cluster_analysis_may24/individual_plots/nonreversal/",as.character(k),'/clustercomp.pdf')),width=6,height=6)
  }
}

#plotSpeciesnonrev(datLfcTime, c("k__Bacteria|p__Firmicutes|c__Erysipelotrichia|o__Erysipelotrichales|f__Erysipelotrichaceae|g__Dubosiella|s__Dubosiella_newyorkensis|t__SGB6789"))

#plotSpeciesnonrev(datLfcTime, c("k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Peptostreptococcaceae|g__Clostridioides|s__Clostridioides_difficile|t__SGB6136"))

#plotSpeciesnonrev(datLfcTime, c("k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Rikenellaceae|g__Alistipes|s__Alistipes_finegoldii|t__SGB2301"))

#plotSpeciesnonrev(datLfcTime, c("k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Oscillospiraceae|g__Oscillibacter|s__Oscillibacter_sp_1_3|t__SGB7266"))

#plotSpeciesnonrev(datLfcTime, c("k__Bacteria|p__Verrucomicrobia|c__Verrucomicrobiae|o__Verrucomicrobiales|f__Akkermansiaceae|g__Akkermansia|s__Akkermansia_muciniphila|t__SGB9226"))

#plotSpeciesnonrev(datLfcTime, c("k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae|g__Ligilactobacillus|s__Ligilactobacillus_murinus|t__SGB7077"))

#plotSpeciesnonrev(datLfcTime, c("k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae|g__Lactobacillus|s__Lactobacillus_johnsonii|t__SGB7041"))

#plotSpeciesnonrev(datLfcTime, c("k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Lactococcus|s__Lactococcus_lactis|t__SGB7985"))

#plotSpeciesnonrev(datLfcTime, c("k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Lactococcus|s__Lactococcus_lactis|t__SGB7984"))

#plotSpeciesnonrev(datLfcTime, c("k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae|g__Limosilactobacillus|s__Limosilactobacillus_reuteri|t__SGB7095"))

#plotSpeciesnonrev(datLfcTime, c("k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_thetaiotaomicron|t__SGB1861"))

#plotSpeciesnonrev(datLfcTime, c("k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Bifidobacteriales|f__Bifidobacteriaceae|g__Bifidobacterium|s__Bifidobacterium_pseudolongum|t__SGB17279"))

#################### Clustering/Plotting including reversal ##########################
{
  outName = "Trajectory_reversal_species_pseudo0.1"
  
  comparisonOrder = c("OnDiet_D0", "OnDiet_M1", "OnDiet_M4", "OnDiet_M6", "OnDiet_M9", "OnDiet_M11", "Rev4MR_M4R2", "Rev4MR_M4R5", "Rev4MR_M4R7", "Rev9MR_M9R2")
  plotOrder = c("OnDiet_D0", "OnDiet_M1", "OnDiet_M4", "OnDiet_M6", "OnDiet_M9", "OnDiet_M11",
                "Rev4MR_D0", "Rev4MR_M1", "Rev4MR_M4", "Rev4MR_M4R2", "Rev4MR_M4R5", "Rev4MR_M4R7",
                "Rev9MR_D0", "Rev9MR_M1", "Rev9MR_M4", "Rev9MR_M6", "Rev9MR_M9", "Rev9MR_M9R2")
  
  setwd('~/Dropbox (Mason Lab)/mouse_diet_semir/cluster_analysis_may24/')
  
  merged_data = readRDS('~/Dropbox (Mason Lab)/mouse_diet_semir/intermediate_files/metaphlan4_all_merged_data_complete.rds')
  
  merged_data_cd = merged_data%>% filter(COHORT == 'Cohort2') %>% filter(grepl('t__',SPECIES)) %>% mutate(SGB = strsplit(SPECIES,'t__') %>% map_chr(2))
  #### OF INTEREST REVERSAL
  metaphlanregs=readRDS('~/Dropbox (Mason Lab)/mouse_diet_semir/association_output/metaphlan4_reversal_associations.rds')
  
  ###### NEED TO CHECK THIS
  diet = metaphlanregs[[1]] %>% distinct %>% filter(grepl('t__',variable))%>% mutate(SGB = strsplit(variable,'t__') %>% map_chr(2))%>% mutate(BY = p.adjust(p.value,method = 'BY'))%>% mutate(BH = p.adjust(p.value,method = 'BH')) %>% select(-BY,-BONFERRONI) %>% filter(BH<0.05)

  #  diet = metaphlanregs[[2]] %>% filter(grepl('t__',variable))%>% mutate(SGB = strsplit(variable,'t__') %>% map_chr(2))%>% mutate(BY = p.adjust(p.value,method = 'BY'))%>% mutate(BH = p.adjust(p.value,method = 'BH')) %>% select(-BY,-BONFERRONI)
   # diet = metaphlanregs[[2]] %>% filter(grepl('t__',variable))%>% mutate(SGB = strsplit(variable,'t__') %>% map_chr(2))%>% mutate(BY = p.adjust(p.value,method = 'BY'))%>% mutate(BH = p.adjust(p.value,method = 'BH')) %>% select(-BY,-BONFERRONI)
 # tokeep = diet %>% mutate(dir = if_else(estimate > 0, 'Positive', 'Negative')) %>% filter(BH < 0.05) %>% dplyr::group_by(term, dir, SGB, dietstatus) %>% dplyr::count() %>% filter((dietstatus == 'POST 4 MONTH REVERSAL' & n > 1) | (dietstatus == 'POST 9 MONTH REVERSAL' & n > 0)) %>% select(-n) %>% distinct() %>% ungroup() %>% select(-dir) %>% group_by(term, SGB) %>% dplyr::count() %>% filter(n == 1) %>% select(-n)  %>% distinct
  #diet = inner_join(diet,tokeep)
  # tokeep = diet %>% mutate(dir = if_else(estimate>0,'Positive','Negative')) %>% filter(BH<0.05) %>% dplyr::group_by(term,dir,SGB) %>% dplyr::count() %>% filter(n>1) %>% select(-n) %>% distinct%>% ungroup %>% select(-dir) %>% group_by(term,SGB) %>% dplyr::count() %>% filter(n==1) %>% select(-n) 
  
  diet$term = gsub('DIET','',diet$term)
  diet$variable = as.character(diet$variable)
  tokeep  = diet %>% filter(BH<0.05) %>% select(SGB) %>% distinct %>% unlist %>% unname
  merged_data_cd =merged_data_cd %>% filter(SGB %in% tokeep)
  
  merged_data_cd$SPECIES = as.character(merged_data_cd$SPECIES)
  merged_data_cd$COHORT = as.factor(merged_data_cd$COHORT)
  merged_data_cd$DIET = factor(merged_data_cd$DIET, levels = c('Control','Coconut Oil','Fish Oil','Ketogenic','Lard','Milkfat','Olive Oil','Palm Oil'))
  merged_data_cd$ANIMAL_ID = as.factor(merged_data_cd$ANIMAL_ID)
  
  merged_data_cd$Rank = splitGet( splitGetFromEnd(merged_data_cd$SPECIES, "\\|", 1), "__", 1)
  #merged_data_cd$Rank = plyr::mapvalues(merged_data_cd$Rank, c("k", "p", "c", "o", "f", "g", "s"), c("1K", "2P", "3C", "4O", "5F", "6G", "8T"))
  merged_data_cd$Rank = plyr::mapvalues(merged_data_cd$Rank, c("t"), c("8T"))
  
  ## ??
  merged_data_cd$TIMEPOINT = mapvalues(merged_data_cd$TIMEPOINT, "M9R0", "M9") ## convert M9R0 to non-reversal?
  
  #merged_data_cd$Type[ merged_data_cd$REVERSAL %in% "NON_REVERSAL"] = "OnDiet"
  merged_data_cd$Type = "OnDiet"
  merged_data_cd$Type[ merged_data_cd$PREPOST_REVERSAL_4MONTH %in% 1	] = "Rev4MR"
  merged_data_cd$Type[ merged_data_cd$PREPOST_REVERSAL_9MONTH %in% 1 & merged_data_cd$TIMEPOINT != "M9" ] = "Rev9MR"
  
  merged_data_cd$TimeType = paste0(merged_data_cd$Type, "_", merged_data_cd$TIMEPOINT)
  merged_data_cd$TimeType = factor(merged_data_cd$TimeType, levels = comparisonOrder)
  merged_data_cd$DietTime = paste0( gsub(" ", "", as.character(merged_data_cd$DIET)), "_", merged_data_cd$TimeType)	
  
  metadata = unique(merged_data_cd[, c("UID", "TIMEPOINT", "DIET", "SEX", "BATCH", "COHORT", "DietTime", "Type", "TimeType")])
  metadata = metadata[ order(metadata$DIET, metadata$TimeType), ]
  
  pseudocount = 0.5
  
  ## Calculate log2FC based on median of log2FC values across all replicate pairs
  datWide = reshape2::dcast(merged_data_cd, SPECIES+Rank~UID, value.var="ABUNDANCE")
  datWide[, 3:ncol(datWide)] = datWide[, 3:ncol(datWide)] + pseudocount
  datLfc = datWide[, c(1,2)]
  for(curDiet in unique(metadata$DIET))
  {
    if(curDiet == "Control") next
    for(curTimeType in unique(metadata$TimeType))
    {
      curCol = paste0( gsub(" ", "", as.character(curDiet)), "_", curTimeType)
      
      curTime = splitGet(curTimeType, "_", 2)
      samplesCase = metadata$UID[ metadata$DIET %in% curDiet & metadata$TIMEPOINT %in% curTime]
      
      curTimeCtrl = plyr::mapvalues(curTime, c("M4R2", "M9R0", "M4R5", "M4R7", "M9R2"), c("M6", "M9", "M9", "M11", "M11"), warn_missing=F) ## Convert reversal time to chronological time for controls
      samplesCtrl = metadata$UID[ metadata$DIET %in% "Control" & metadata$TIMEPOINT %in% curTimeCtrl]
      
      curLfc = data.frame(row.names = datWide[, 1])
      curLfc[, curCol] = PairwiseCompare(datWide, samplesCase, samplesCtrl, operation="logratio", summaryFun="median" )
      datLfc = cbind(datLfc, curLfc)
    }
  }
  
  datLfcLong = reshape2::melt(datLfc, c("SPECIES", "Rank"), variable.name="DietTime", value.name="Log2FC")
  datLfcLong$Diet = splitGet(datLfcLong$DietTime, "_", 1)
  datLfcLong$Type = splitGet(datLfcLong$DietTime, "_", 2)
  datLfcLong$Time = splitGet(datLfcLong$DietTime, "_", 3)
  datLfcLong$TimeType = paste0(datLfcLong$Type, "_", datLfcLong$Time)
  
  datLfcTime = reshape2::dcast(datLfcLong, SPECIES+Rank+Diet~TimeType, value.var="Log2FC")
  write.tsv(datLfcTime, paste0(outName, "_log2FC_wide.txt"), row.names = F)
  
  datLfcTime = datLfcTime[ datLfcTime$Rank != "1K", ]
  
  ## ?? Filter to species or not ??
  datLfcTime = datLfcTime[ datLfcTime$Rank == "8T", ] 
  
  rownames(datLfcTime) = paste0(datLfcTime$SPECIES, ":", datLfcTime$Diet)
  datLfcTime2 = datLfcTime[, -c(1:3)]
  datLfcTime2 = datLfcTime2[ rowSds(as.matrix(datLfcTime2)) > 0, ] ## remove rows with all zeroes
  #####
  
  
  library(e1071)
  library(tidyr)
  library(dplyr)
  library(pheatmap)
  
  dataList = list(LFC_sig=datLfcTime2)
  kList = c(20)
  distMethod = "euclidean"
  for(k in kList)
  {
    for(i in 1:length(dataList))
    {
      curData = dataList[[i]]
      curName = names(dataList)[i]
      curOut = paste0(outName, "_cmeans_", curName, "_k", k, "_" , distMethod)
      curData = curData[ complete.cases(curData), ]
      
      ## Normalize/standardize data by some metric?
      curData2 = curData
      #curData2 = t(scale(t(curData)))
      #curData = curData - rep.col(rowMeans(curData), ncol(curData))
      
      ## M estimate
      m = mestimate(curData2)
      
      if( ! file.exists(paste0(curOut, ".rds")) )
      {
        fcm_results = cmeans(curData2, centers = k, m = m, iter.max = 2000, dist=)
        center_order = getTSPOrder(fcm_results$centers)
        fcm_results$centers = fcm_results$centers[ center_order, ]
        rownames(fcm_results$centers) = paste0("C_", sprintf("%02d", 1:k))
        fcm_results$membership = fcm_results$membership[, center_order]
        colnames(fcm_results$membership) = paste0("C_", sprintf("%02d", 1:k))
        
        fcm_results$cluster = mapvalues(fcm_results$cluster, center_order, paste0("C_", sprintf("%02d", 1:k)))
        fcm_results$size = fcm_results$size[as.numeric(center_order)]
        saveRDS(fcm_results, paste0(curOut, ".rds"))
      }else
      {	
        fcm_results = readRDS(paste0(curOut, ".rds"))
      }
      
      
      fcm_centroids = fcm_results$centers
      fcm_centroids_df = data.frame(fcm_centroids)
      fcm_centroids_df$cluster = row.names(fcm_centroids_df)
      centroids_long = tidyr::gather(fcm_centroids_df, "Type", "log2FC", 1:ncol(fcm_centroids))
      centroids_long$Type = factor(centroids_long$Type, levels = plotOrder)
      
      #start with the input data
      fcm_plotting_df = data.frame(curData)
      
      #add genes
      fcm_plotting_df$gene = row.names(fcm_plotting_df)
      
      #bind cluster assinment
      fcm_plotting_df$cluster = fcm_results$cluster
      #fetch the membership for each gene/top scoring cluster
      fcm_plotting_df$membership = sapply(1:length(fcm_plotting_df$cluster), function(row)
      {
        clust = fcm_plotting_df$cluster[row]
        fcm_results$membership[row, clust]
      })
      
      
      cluster_plot_df = fcm_plotting_df %>% group_by(cluster) %>%
        dplyr::select(., 1:ncol(curData), membership, gene, cluster) %>%
        tidyr::gather(., "Type", "log2FC", 1:ncol(curData))
      
      cluster_plot_df = as.data.frame(cluster_plot_df)
      cluster_plot_df$Type = factor(cluster_plot_df$Type, levels = comparisonOrder)
      
      #order the dataframe by score
      cluster_plot_df = cluster_plot_df[order(cluster_plot_df$membership),]
      #set the order by setting the factors using forcats
      cluster_plot_df$gene = forcats::fct_inorder(cluster_plot_df$gene)
      
      cluster_plot_df$group = paste0(cluster_plot_df$gene, "_", splitGet(as.character(cluster_plot_df$Type), "_", 1))
      cluster_plot_df$Type = factor(cluster_plot_df$Type, levels = comparisonOrder)
      cluster_plot_df$Diet = splitGet(cluster_plot_df$gene, ":", 2)
      cluster_plot_df = cluster_plot_df[ order(cluster_plot_df$Type), ]
      
      write.tsv(cluster_plot_df, paste0(curOut, "_long.txt"), row.names=F)
      
      
      ## cluster centroids as trajectory lines
      ## Hacky code to clone pre-reversal timepoint rows for the reversal samples so they can be plotted separately
      core = centroids_long
      core$Linetype = "1solid"
      
      # Clone for M4R
      coreTmp1a = core[ grep("OnDiet", core$Type),]
      coreTmp1a = coreTmp1a[ grep("(D0|M1|M4)", coreTmp1a$Type), ]
      coreTmp1a$Type = gsub("OnDiet", "Rev4MRPre", as.character(coreTmp1a$Type))
      coreTmp1a$Linetype = "2dashed"
      
      coreTmp1b = core[ grep("OnDiet", core$Type),]
      coreTmp1b = coreTmp1b[ grep("(M4)", coreTmp1b$Type), ]
      coreTmp1b$Type = gsub("OnDiet", "Rev4MR", as.character(coreTmp1b$Type))
      coreTmp1b$Linetype = "1solid"
      
      # Clone for M9R
      coreTmp2a = core[ grep("OnDiet", core$Type),]
      coreTmp2a = coreTmp2a[ grep("(D0|M1|M4|M6|M9)", coreTmp2a$Type), ]
      coreTmp2a$Type = gsub("OnDiet", "Rev9MRPre", as.character(coreTmp2a$Type))
      coreTmp2a$Linetype = "2dashed"
      
      coreTmp2b = core[ grep("OnDiet", core$Type),]
      coreTmp2b = coreTmp2b[ grep("(M9)", coreTmp2b$Type), ]
      coreTmp2b$Type = gsub("OnDiet", "Rev9MR", as.character(coreTmp2b$Type))
      coreTmp2b$Linetype = "1solid"
      
      # merge all
      core = rbind(core, coreTmp1a, coreTmp1b, coreTmp2a, coreTmp2b)
      #core$Type = factor(core$Type, levels = plotOrder)
      core$group = paste0(core$cluster, "_", splitGet(as.character(core$Type), "_", 1))
      core = core[ order(core$Type), ]
      core$Type = gsub("Pre", "", core$Type)
      
      ## More hacky code to clone pre-reversal timepoint rows for the reversal samples so they can be plotted separately 
      cluster_plot_df2Tmp1 = cluster_plot_df[ cluster_plot_df$Type == "OnDiet_M4", ]
      cluster_plot_df2Tmp1$Type = "Rev4MR_M4"
      cluster_plot_df2Tmp1$group = gsub("OnDiet", "Rev4MR", cluster_plot_df2Tmp1$group)
      
      cluster_plot_df2Tmp2 = cluster_plot_df[ cluster_plot_df$Type == "OnDiet_M9", ]
      cluster_plot_df2Tmp2$Type = "Rev9MR_M9"
      cluster_plot_df2Tmp2$group = gsub("OnDiet", "Rev9MR", cluster_plot_df2Tmp2$group)
      
      cluster_plot_df2 = rbind(cluster_plot_df, cluster_plot_df2Tmp1, cluster_plot_df2Tmp2 )
      
      ### filter 
      cluster_plot_df2 = cluster_plot_df2 %>% filter(membership>membcut)
      
      #clusterstokeep = c("C_01", "C_02", "C_03", "C_04", "C_05", "C_06", "C_07", "C_08", "C_09", "C_10", "C_11", "C_12", "C_13", "C_14", "C_15", "C_16", "C_17", "C_18", "C_19", "C_20")
      clusterstokeep =c("C_13","C_02","C_15","C_08","C_19","C_12","C_18","C_14","C_10","C_11","C_16","C_03")

      cluster_plot_df2 = cluster_plot_df2 %>% filter(cluster %in% clusterstokeep)
      cluster_plot_df2$cluster = factor(cluster_plot_df2$cluster,levels = clusterstokeep) 
      
      core = core %>% filter(cluster %in% clusterstokeep)
      core$cluster = factor(core$cluster,levels = clusterstokeep) 
      
      ggplot(cluster_plot_df2, aes(x = Type, y = log2FC)) + 
        geom_line(aes(colour = membership, group = group)) +
        scale_colour_gradientn(colours=c('yellow1','red2')) +
        geom_line(data=core, aes(Type, log2FC, group = group, linetype=Linetype), color="black",inherit.aes=FALSE, linewidth=0.8) +
        scale_x_discrete(name="Timepoint", limits=plotOrder) +
        ylab("log2 fold change vs time-matched control diet") +
        labs(title= paste0("All clusters log2FC vs time-matched Control Diet"), color = "Score") + facet_wrap(~cluster, nrow=2,scales="free") + theme_cem + 
        geom_hline(yintercept=0) 
      ggsave(paste0(curOut, "_clusterProfile_FORFIG.pdf"), width=20, height=8)

      ### NOW TO GENERATE SOME OTHER INFO
      cluster_plot_df2 = cluster_plot_df2 %>% mutate(SGB = strsplit(group,':') %>% map_chr(1) %>% strsplit('t__') %>% map_chr(2))
      cluster_plot_df2 = left_join(cluster_plot_df2,gtdbmap)  %>% separate(GTDB_TAX, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = TRUE, extra = "merge")
      # get total number of lines in each cluster
      clustersizes = cluster_plot_df2 %>% select(cluster,group) %>% distinct %>% group_by(cluster) %>% select(-group)%>% table %>% data.frame %>% dplyr::rename(cluster = ".",cluster_counts = Freq)
      # find the number of unique species per cluster
      speciescounts = cluster_plot_df2 %>% select(-Type,-log2FC) %>% distinct %>% select(cluster,Species) %>% table %>% data.frame %>% filter(Freq>0)
      colnames(speciescounts) = c('cluster','Species','Species_Count')
      # fund the number of unique diets per cluster
      dietcounts = cluster_plot_df2 %>% select(-Type,-log2FC) %>% distinct %>% select(cluster,Diet) %>% table %>% data.frame %>% filter(Freq>0)
      colnames(dietcounts) = c('cluster','Diet','Diet_Count')
      clustersumsspecies = left_join(clustersizes,speciescounts) %>% mutate(species_frac = Species_Count/cluster_counts) %>% filter(Species != 's__')
      clustersumsdiets = left_join(clustersizes,dietcounts) %>% mutate(diet_frac = Diet_Count/cluster_counts)
      
      # make an associated stacked barplot
      long = dietcounts %>% dplyr::group_by(cluster) %>% mutate(total = sum(Diet_Count),prop = Diet_Count/total)
      longplot = ggplot(long %>% mutate(count=strsplit(as.character(cluster),'_') %>% map_chr(2) %>% as.numeric),aes(y = fct_reorder(cluster,(desc(count))), x = prop*100, fill = Diet)) + geom_bar(stat = 'identity') + theme_minimal() + scale_fill_manual(values=colorsDiet) + theme(legend.position = 'none') + xlab('%') + ylab('')
      ggsave(plot = longplot,path.expand(paste0("~/Dropbox (Mason Lab)/mouse_diet_semir/cluster_analysis_may24/individual_plots/reversal/",as.character(k),'/clustercomp_',as.character(k),'.pdf')),width=6,height=6)
      # generate tree heatmap thing -- reversal
      #wide = clustersumsspecies %>% filter(Species_Count>1) %>% group_by(cluster) %>% slice_max(n=5,species_frac) %>% reshape2::dcast(cluster ~ Species,value.var = 'species_frac') %>% column_to_rownames('cluster')
      wide = clustersumsspecies %>% group_by(cluster)  %>% reshape2::dcast(cluster ~ Species,value.var = 'species_frac') %>% column_to_rownames('cluster')
      wide[is.na(wide)]= 0
      metatreesub = keep.tip(tip = colnames(wide),phy = metatree)
      wide = wide %>% t() %>% data.frame(check.names = F)#%>% rownames_to_column('id')
      mapfilesub = mapfile %>% filter(Species %in% rownames(wide)) %>% select(Species,Phylum,Class,Order,Family,Genus) 
      rownames(mapfilesub) =NULL
      mapfilesub = mapfilesub %>% column_to_rownames('Species')
      mapfilesub[metatreesub$tip.label,] = mapfilesub
      wide[metatreesub$tip.label,] = wide
      mapfilesub = mapfilesub %>% rownames_to_column('Species')
      tree = ggtree(metatreesub,layout = 'dendrogram') 
      tree = tree %<+% mapfilesub
      tree = tree + geom_tiplab(align = TRUE,aes(color=Phylum)) + geom_treescale()
      tree1 = gheatmap(tree, wide, offset = .5, width = 1, colnames_angle = 0, hjust = 1)  + scale_fill_gradientn(colors = c('white', "red", "darkblue"), values = scales::rescale(c(0, .5, 1)), limits = c(0, 1), guide = "colourbar")
      ggsave(plot = tree1,path.expand(paste0("~/Dropbox (Mason Lab)/mouse_diet_semir/cluster_analysis_may24/individual_plots/reversal/",as.character(k),'/clustertree.pdf')),width=14,height=10)
      
      # now make the tree at the genus level
      clustersumsspeciesgen = clustersumsspecies %>% mutate(genus = strsplit(as.character(Species),' ') %>% map_chr(1) %>% gsub('s__','',.)) %>% group_by(genus,cluster) %>% mutate(Genus_Count = sum(Species_Count),genus_frac = Genus_Count/cluster_counts) %>% ungroup %>% group_by(cluster,genus) %>% sample_n(1) 
      metatreesub = keep.tip(tip = as.character(clustersumsspeciesgen$Species),phy = metatree)
      wide = clustersumsspeciesgen %>% group_by(cluster) %>% reshape2::dcast(cluster ~ genus,value.var = 'species_frac') %>% column_to_rownames('cluster')
      wide[is.na(wide)]= 0
      metatreesub$tip.label = strsplit(metatreesub$tip.label,' ') %>% map_chr(1) %>% gsub('s__','',.)
      metatreesub = keep.tip(tip = colnames(wide),phy = metatreesub)
      wide = wide %>% t() %>% data.frame(check.names = F)#%>% rownames_to_column('id')
      mapfilesub = mapfile %>% mutate(Genus = gsub('g__','',Genus))%>% filter(Genus %in% rownames(wide)) %>% select(Genus,Phylum,Class,Order,Family) %>% distinct 
      rownames(mapfilesub) =NULL
      mapfilesub = mapfilesub %>% column_to_rownames('Genus')
      mapfilesub[metatreesub$tip.label,] = mapfilesub
      wide[metatreesub$tip.label,] = wide
      mapfilesub = mapfilesub %>% rownames_to_column('Genus')
      tree = ggtree(metatreesub,layout = 'dendrogram') 
      tree = tree %<+% mapfilesub
      tree = tree + geom_tiplab(align = TRUE,aes(color=Phylum)) + geom_treescale()
      tree1 = gheatmap(tree, wide, offset = .5, width = 1, colnames_angle = 0, hjust = 1)  + scale_fill_gradientn(colors = c('white', "red", "darkblue"), values = scales::rescale(c(0, .5, 1)), limits = c(0, 1), guide = "colourbar")
      ggsave(plot = tree1,path.expand(paste0("~/Dropbox (Mason Lab)/mouse_diet_semir/cluster_analysis_may24/individual_plots/reversal/",as.character(k),'/clustertree_GENUS.pdf')),width=14,height=10)
      
      # now for the reversals pull out some specific bugs
      diet1 = metaphlanregs[[1]] %>% filter(grepl('t__',variable))%>% mutate(SGB = strsplit(variable,'t__') %>% map_chr(2))%>% mutate(BY = p.adjust(p.value,method = 'BY'))%>% mutate(BH = p.adjust(p.value,method = 'BH')) %>% select(-BY,-BONFERRONI)
      diet2 = metaphlanregs[[2]] %>% filter(grepl('t__',variable))%>% mutate(SGB = strsplit(variable,'t__') %>% map_chr(2))%>% mutate(BY = p.adjust(p.value,method = 'BY'))%>% mutate(BH = p.adjust(p.value,method = 'BH')) %>% select(-BY,-BONFERRONI)
      
      tokeep2 = left_join(bind_rows(diet1,diet2),gtdbmap) %>% separate(GTDB_TAX, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = TRUE, extra = "merge")  %>% filter(BH<0.05) %>% filter(Species %in% c('s__Alistipes finegoldii','s__Limosilactobacillus reuteri','s__Lactobacillus johnsonii','s__Bifidobacterium globosum','s__Akkermansia muciniphila','s__Ligilactobacillus murinus','s__Bacteroides thetaiotaomicron','s__Clostridioides difficile')) %>% select(Species,term) %>% dplyr::rename(Diet = term)%>% distinct
      tokeep2$Diet = gsub('DIET','',tokeep2$Diet)
      cluster_plot_df3 = inner_join(cluster_plot_df2,tokeep2)
      #core2 = core %>% filter(Species %in% tokeep2)
      ggplot(cluster_plot_df3, aes(x = Type, y = log2FC)) +  geom_line(aes(colour = Diet, group = group)) +scale_colour_manual(values=colorsDiet) + scale_x_discrete(name="Timepoint", limits=plotOrder) + ylab("log2 fold change vs time-matched control diet") + facet_wrap(~Species, scales="free") + theme_cem + geom_hline(yintercept=0) 
      longplot = ggplot(clustersumsdiets,aes(y = factor(cluster,levels=rev(clusterstokeep)), x = diet_frac*100, fill = Diet)) + geom_bar(stat = 'identity') + theme_minimal() + scale_fill_manual(values=colorsDiet) + theme(legend.position = 'none') + xlab('%') + ylab('')
      ggsave(plot = longplot,path.expand(paste0("~/Dropbox (Mason Lab)/mouse_diet_semir/cluster_analysis_may24/individual_plots/reversal/",as.character(k),'/clustercomp.pdf')),width=3,height=7)
      }
  }
}

### COUNT PERSISTENT VS REVERTED

diet = metaphlanregs[[1]] %>% filter(grepl('t__',variable))%>% mutate(SGB = strsplit(variable,'t__') %>% map_chr(2))%>% mutate(BY = p.adjust(p.value,method = 'BY'))%>% mutate(BH = p.adjust(p.value,method = 'BH')) %>% select(-BY,-BONFERRONI)
tokeep = diet %>% mutate(dir = if_else(estimate>0,'Positive','Negative')) %>% filter(BH<0.05) %>% dplyr::group_by(term,dir,SGB) %>% dplyr::count() %>% filter(n>1) %>% select(-n) %>% distinct%>% ungroup %>% select(-dir) %>% group_by(term,SGB) %>% dplyr::count() %>% filter(n==1) %>% select(-n) 
diet = inner_join(diet,tokeep) 
diet = diet %>% mutate()

diet2 = metaphlanregs[[2]] %>% filter(grepl('t__',variable))%>% mutate(SGB = strsplit(variable,'t__') %>% map_chr(2))%>% mutate(BY = p.adjust(p.value,method = 'BY'))%>% mutate(BH = p.adjust(p.value,method = 'BH')) %>% select(-BY,-BONFERRONI)
tokeep = diet2 %>% mutate(dir = if_else(estimate > 0, 'Positive', 'Negative')) %>% filter(BH < 0.05) %>% dplyr::group_by(term, dir, SGB, dietstatus) %>% dplyr::count() %>% filter((dietstatus == 'POST 4 MONTH REVERSAL' & n > 1) | (dietstatus == 'POST 9 MONTH REVERSAL' & n > 0)) %>% select(-n) %>% distinct() %>% ungroup() %>% select(-dir) %>% group_by(term, SGB, dietstatus) %>% dplyr::count() %>% filter(n == 1) %>% select(-n)
diet2 = inner_join(diet2,tokeep)

diet_sub = diet %>% filter(BH<0.05) %>% select(term,SGB,estimate) %>% mutate(direction_pre = if_else(estimate > 0, 1, -1)) %>% select(-estimate)%>% ungroup %>% distinct %>% mutate(dietstatus = 'POST 4 MONTH REVERSAL')
diet_sub_tmp = diet %>% filter(BH<0.05) %>% select(term,SGB,estimate) %>% mutate(direction_pre = if_else(estimate > 0, 1, -1)) %>% select(-estimate)%>% ungroup %>% distinct %>% mutate(dietstatus = 'POST 9 MONTH REVERSAL')
diet_sub = bind_rows(diet_sub,diet_sub_tmp)
diet2_sub = diet2 %>% filter(BH<0.05) %>% select(term,SGB,estimate,dietstatus) %>% mutate(direction_post = if_else(estimate > 0, 1, -1)) %>% select(-estimate)%>% ungroup%>% distinct

merged = dplyr::left_join(diet_sub,diet2_sub) %>% distinct %>% dplyr::group_by(SGB,term)%>% mutate(direction_post= if_else(is.na(direction_post),0,direction_post)) %>% mutate(type = if_else(direction_pre == -1 & direction_post==-1,'Persistently Decreased','none')) %>% mutate(type = if_else(direction_pre == 1 & direction_post== 1,'Persistently Increased',type)) %>% mutate(type = if_else(direction_pre != direction_post,'Reverted',type)) %>% ungroup

summarised = merged %>% dplyr::group_by(term,type,dietstatus) %>% dplyr::summarise(n = n())
summarised$term = gsub('DIET','',summarised$term)

all_categories <- c("Persistently Decreased", "Persistently Increased", "Reverted")
all_combinations <- summarised %>% ungroup %>% dplyr::distinct(term, dietstatus) %>% tidyr::crossing(type = all_categories, .name_repair = "unique")

summarised_full <- all_combinations %>% dplyr::left_join(summarised, by = c("term", "type", "dietstatus")) %>% dplyr::mutate(n = ifelse(is.na(n), 0, n))

# by diet
ggplot(summarised_full,aes(y = (n+1),x = term,fill = factor(type))) + geom_bar(stat='identity',position = 'dodge',width=0.9) + facet_wrap(. ~ dietstatus,scales='free_x',nrow=1) + theme_cowplot()+ scale_fill_brewer(palette = 'Set1') + theme(axis.text.x = element_text(angle = 60,hjust = 1)) + scale_y_log10() + xlab('') + ylab('# of Species')
ggsave('~/Dropbox (Mason Lab)/mouse_diet_semir/plots/diet_reversion_summary.pdf',width=8,height=2.5)

# overall
dat = summarised_full %>% dplyr::group_by(dietstatus,type) %>% dplyr::summarise(total = sum(n))
ggplot(dat,aes(y = (total+1),x = type,fill = factor(type))) + geom_bar(stat='identity',position = 'dodge',width=0.9) + theme_cowplot()+ scale_fill_brewer(palette = 'Set1') + theme(axis.text.x = element_text(angle = 60,hjust = 1)) + scale_y_log10() + theme(legend.position = 'none') + xlab('') + ylab('# of Species')+facet_wrap(. ~ dietstatus,scales='free_x',nrow=1)
ggsave('~/Dropbox (Mason Lab)/mouse_diet_semir/plots/overall_reversion_summary.pdf',width=2,height=2.5)

### Get overlap in bugs between categories at 4 and 9 month
diet_sub = diet %>% filter(BH<0.05) %>% select(term,SGB,estimate,variable) %>% mutate(direction_pre = if_else(estimate > 0, 1, -1)) %>% select(-estimate)%>% ungroup %>% distinct %>% mutate(dietstatus = 'POST 4 MONTH REVERSAL')
diet_sub_tmp = diet %>% filter(BH<0.05) %>% select(term,SGB,estimate,variable) %>% mutate(direction_pre = if_else(estimate > 0, 1, -1)) %>% select(-estimate)%>% ungroup %>% distinct %>% mutate(dietstatus = 'POST 9 MONTH REVERSAL')
diet_sub = bind_rows(diet_sub,diet_sub_tmp)
diet2_sub = diet2 %>% filter(BH<0.05) %>% select(term,SGB,estimate,dietstatus,variable) %>% mutate(direction_post = if_else(estimate > 0, 1, -1)) %>% select(-estimate)%>% ungroup%>% distinct

merged = dplyr::left_join(diet_sub,diet2_sub) %>% distinct %>% dplyr::group_by(SGB,term)%>% mutate(direction_post= if_else(is.na(direction_post),0,direction_post)) %>% mutate(type = if_else(direction_pre == -1 & direction_post==-1,'Persistently Decreased','none')) %>% mutate(type = if_else(direction_pre == 1 & direction_post== 1,'Persistently Increased',type)) %>% mutate(type = if_else(direction_pre != direction_post,'Reverted',type)) %>% ungroup
merged$term = gsub('DIET','',merged$term)

wide = merged  %>% mutate(val = 1,termterm = paste(term,dietstatus,sep=' -- ')) %>% select(SGB,type,termterm,val) %>% distinct%>% reshape2::dcast(SGB + type ~ termterm,value.var = 'val')
wide[is.na(wide)]=0

#### UPSETs
diet_sets <- list(
  c("Ketogenic -- POST 4 MONTH REVERSAL", "Ketogenic -- POST 9 MONTH REVERSAL"),
  "Ketogenic -- POST 4 MONTH REVERSAL", "Ketogenic -- POST 9 MONTH REVERSAL",
  c("Coconut Oil -- POST 4 MONTH REVERSAL", "Coconut Oil -- POST 9 MONTH REVERSAL"), 
  "Coconut Oil -- POST 4 MONTH REVERSAL", "Coconut Oil -- POST 9 MONTH REVERSAL",
  c("Lard -- POST 4 MONTH REVERSAL", "Lard -- POST 9 MONTH REVERSAL"),
  "Lard -- POST 4 MONTH REVERSAL", "Lard -- POST 9 MONTH REVERSAL",
  c("Fish Oil -- POST 4 MONTH REVERSAL", "Fish Oil -- POST 9 MONTH REVERSAL"),
  "Fish Oil -- POST 4 MONTH REVERSAL", "Fish Oil -- POST 9 MONTH REVERSAL",
  c("Olive Oil -- POST 4 MONTH REVERSAL", "Olive Oil -- POST 9 MONTH REVERSAL"),
  "Olive Oil -- POST 4 MONTH REVERSAL", "Olive Oil -- POST 9 MONTH REVERSAL",
  c("Palm Oil -- POST 4 MONTH REVERSAL", "Palm Oil -- POST 9 MONTH REVERSAL"),
  "Palm Oil -- POST 4 MONTH REVERSAL", "Palm Oil -- POST 9 MONTH REVERSAL",
  c("Milkfat -- POST 4 MONTH REVERSAL", "Milkfat -- POST 9 MONTH REVERSAL"),
  "Milkfat -- POST 4 MONTH REVERSAL", "Milkfat -- POST 9 MONTH REVERSAL"
)


pdf('~/Dropbox (Mason Lab)/mouse_diet_semir/plots/upset_reversal.pdf',width=7,height=7)
ComplexUpset::upset(data = (wide),min_size=2,intersect = colnames(wide)[3:16],intersections = diet_sets,base_annotations=list('Intersection size'=intersection_size(mapping=aes(fill = type))+theme(legend.position = 'None')+scale_fill_brewer(palette = 'Set1')),width_ratio=0.1,sort_intersections=F)
dev.off()


plotSpecies = function(datLfcTime, speciesList)
{
  speciesNames =gsub('\\|t','',splitGetFromEnd(speciesList, "__", 2))
  speciesNames = gsub("_", " ", speciesNames)
  speciesNames = paste0(speciesNames, collapse=", ")
  
  curData = datLfcTime[ datLfcTime$SPECIES %in% speciesList, ]
  curDataM = reshape2::melt(curData, c("SPECIES", "Rank", "Diet"), variable.name="Type", value.name="log2FC")
  curDataM$group = paste0(curDataM$SPECIES, "_", splitGet(as.character(curDataM$Type), "_", 1))
  curDataM$Linetype = "1solid"
  
  # Clone for M4R
  curDataMTmp1a = curDataM[ grep("OnDiet", curDataM$Type),]
  curDataMTmp1a = curDataMTmp1a[ grep("(D0|M1|M4)", curDataMTmp1a$Type), ]
  curDataMTmp1a$Type = gsub("OnDiet", "Rev4MRPre", as.character(curDataMTmp1a$Type))
  curDataMTmp1a$Linetype = "2dashed"
  
  curDataMTmp1b = curDataM[ grep("OnDiet", curDataM$Type),]
  curDataMTmp1b = curDataMTmp1b[ grep("(M4)", curDataMTmp1b$Type), ]
  curDataMTmp1b$Type = gsub("OnDiet", "Rev4MR", as.character(curDataMTmp1b$Type))
  curDataMTmp1b$Linetype = "1solid"
  
  # Clone for M9R
  curDataMTmp2a = curDataM[ grep("OnDiet", curDataM$Type),]
  curDataMTmp2a = curDataMTmp2a[ grep("(D0|M1|M4|M6|M9)", curDataMTmp2a$Type), ]
  curDataMTmp2a$Type = gsub("OnDiet", "Rev9MRPre", as.character(curDataMTmp2a$Type))
  curDataMTmp2a$Linetype = "2dashed"
  
  curDataMTmp2b = curDataM[ grep("OnDiet", curDataM$Type),]
  curDataMTmp2b = curDataMTmp2b[ grep("(M9)", curDataMTmp2b$Type), ]
  curDataMTmp2b$Type = gsub("OnDiet", "Rev9MR", as.character(curDataMTmp2b$Type))
  curDataMTmp2b$Linetype = "1solid"
  
  curDataM2 = rbind(curDataM, curDataMTmp1a, curDataMTmp1b, curDataMTmp2a, curDataMTmp2b)
  curDataM2$group = paste0(curDataM2$cluster, "_", splitGet(as.character(curDataM2$Type), "_", 1))
  curDataM2 = curDataM2[ order(curDataM2$Type), ]
  curDataM2$Type = gsub("Pre", "", curDataM2$Type)
  
  
  out = ggplot(curDataM2, aes(x = Type, y = log2FC)) +
    geom_line(aes(colour = Diet, group = group, linetype=Linetype), size=1.1) +
    scale_colour_manual(values=colorsDiet2) +
    scale_x_discrete(name="Timepoint", limits=plotOrder) +
    ylab("log2 fold change vs time-matched control diet") +
    facet_wrap(~Diet, scales="free",nrow=1) +  
    geom_hline(yintercept=0) + 
    ggtitle(speciesNames) +
    theme_cem + theme(legend.position="none")+ theme(axis.title.x = element_blank(),axis.text.x = element_blank())
  #ggsave(paste0('~/Dropbox (Mason Lab)/mouse_diet_semir/plots/',speciesList[[1]],'.pdf'),width=24,height=4)
  return(out)
}


tokeep3= diet %>% filter(BH<0.05,dietstatus == 'ON DIET') %>% select(SGB,term) %>% mutate(term = gsub('DIET','',term)) %>% dplyr::rename(Diet = term)  %>% select(Diet,SGB) %>% distinct
datLfcTime_2 = datLfcTime %>% mutate(SGB = strsplit(SPECIES,'t__') %>% map_chr(2))
datLfcTime_2$Diet = gsub('CoconutOil','Coconut Oil',datLfcTime_2$Diet)
datLfcTime_2$Diet = gsub('FishOil','Fish Oil',datLfcTime_2$Diet)
datLfcTime_2$Diet = gsub('OliveOil','Olive Oil',datLfcTime_2$Diet)
datLfcTime_2$Diet = gsub('PalmOil','Palm Oil',datLfcTime_2$Diet)
 
datLfcTime_2 = left_join(datLfcTime_2,tokeep3) %>% select(-SGB)

pdf('~/Dropbox (Mason Lab)/mouse_diet_semir/plots/speciesl2fc/afiegold.pdf',width = 24,height=4)
a=plotSpecies(datLfcTime_2, c("k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Rikenellaceae|g__Alistipes|s__Alistipes_finegoldii|t__SGB2301"))
dev.off()

pdf('~/Dropbox (Mason Lab)/mouse_diet_semir/plots/speciesl2fc/akkmun.pdf',width = 24,height=4)
b=plotSpecies(datLfcTime_2,c(c("k__Bacteria|p__Verrucomicrobia|c__Verrucomicrobiae|o__Verrucomicrobiales|f__Akkermansiaceae|g__Akkermansia|s__Akkermansia_muciniphila|t__SGB9226")))
dev.off()

pdf('~/Dropbox (Mason Lab)/mouse_diet_semir/plots/speciesl2fc/ljohns.pdf',width = 24,height=4)
c=plotSpecies(datLfcTime_2,c("k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae|g__Lactobacillus|s__Lactobacillus_johnsonii|t__SGB7041"))
dev.off()

pdf('~/Dropbox (Mason Lab)/mouse_diet_semir/plots/speciesl2fc/bpseudo.pdf',width = 24,height=4)
d=plotSpecies(datLfcTime_2,c("k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Bifidobacteriales|f__Bifidobacteriaceae|g__Bifidobacterium|s__Bifidobacterium_pseudolongum|t__SGB17279"))
dev.off()

pdf('~/Dropbox (Mason Lab)/mouse_diet_semir/plots/speciesl2fc/btheta.pdf',width = 24,height=4)
e=plotSpecies(datLfcTime_2,c("k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_thetaiotaomicron|t__SGB1861"))
dev.off()


library(patchwork)

pdf('~/Dropbox (Mason Lab)/mouse_diet_semir/plots/speciesl2fc/combined_species.pdf',width=24,height=14)
(a / b / c / d / e) & theme(axis.text.x = element_blank(),axis.title.y = element_blank(), strip.background = element_blank(),strip.text = element_blank(),strip.text.x = element_blank(),axis.text.y = element_text(size=20))
dev.off()













