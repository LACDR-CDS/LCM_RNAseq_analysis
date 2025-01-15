

## Data analysis pipeline BioSpyder datasets for LCM samples ##

# 0 - Data assembly, normalization, QCs

# Lukas Wijaya

# Version 5

## Steps:
# Loading of data
# Library size calculation and filtering
# Exluding specific samples for further analysis (NoCells, and others if necessary)
# Normalization of raw counts using DESEq2 package
# Quality control plots (library size and count distribution, replicate correlation, variance across means)
# DEGs determination: Log2FC and p-value calculation using DESEq2 package
# Filtering of DEGs (based on p-value/Log2FC)
# Assessing gene expression for specific genes
# Principal component analysis
# Hierarchical clustering of samples and genes in heatmap
# Saving data

#############################################

rm(list=ls(all=TRUE))


###### INPUT USER ######
NoCells <- c()

MetaName <- "meta_data.txt"      #use template and save as txt file
CountThreshold <- 10                              #Exclude samples lower than this number for %mapped
VariableQCs <- "LOCATION_ID"                                #Variable to check distribution Library size and counts (choose meanID, CELL_ID, TREATMENT, SAMPLE_ID, EXP_ID, TIMEPOINT or REPLICATE)
RemoveSamples <- NoCells                             #Give list of SAMPLE_IDs for samples which you would like to remove for further analysis (for example MAQC, NoCells samples, outliers)

NormMethod <- "CPM"                                    #Fill in desired normalization method prior to Deseq2 DEG calculation: "CPM" or "DESeq2"
Threshold_padj <- 0.1                                     #Filtering of genes with p-adj lower than threshold
Threshold_log2FC <- 0.1                                   #Filtering of genes with absolute log2FC higher than threshold (genes lower than -0.5, and genes higher than 0.5)
Filtering <- "padj&log2FC"                                      #Fill in filtering method ("padj", "log2FC" or "padj&log2FC")
CheckGenes <- NA #Check expression of individual genes among samples (Give gene symbols)
                                                          #You can list as many genes as you like, however graph may become too full
Geneset <- NA                 #Fill in file name of specific gene list of interest (GeneSymbols, txt format). If not used, fill in NA.
                                                          

# Working directories #

CountTableDir <- "H:/Experiment_file_backup/experiment_2020/20200501_LCM_RealExperiment/Data/Count_table" #Folder containing only raw count table files (csv or excel files)
inputDir <-  'H:/Experiment_file_backup/experiment_2020/20200501_LCM_RealExperiment/Data/' #Store metadata file and raw count files here
outputDir_GR <- 'H:/Experiment_file_backup/experiment_2020/20200501_LCM_RealExperiment/Data/Graphs'

scriptDir <- "H:/Experiment_file_backup/experiment_2020/20200501_LCM_RealExperiment/Script"
source(file.path(scriptDir, "A_BioSpyder analysis pipeline_Lib_DnRData.R"))


## Assembly of Data ##

######  Load data & library size filter ###### 
meta <- read.delim(paste0(inputDir, MetaName), stringsAsFactors = FALSE)
rownames(meta) <- meta[,1]
count <- LoadCountData(CountTableDir)
count$probe <- rownames(count)



meta$meanID <- paste0(meta$EXP_ID, "_",
                      meta$LOCATION_ID, "_",
                      meta$TREATMENT, "_conc",
                      meta$CONCENTRATION,"_TP", 
                      meta$TIMEPOINT, "_rep",
                      meta$REPLICATE)

meta$SAMPLENAME <- paste0(meta$EXP_ID, "_",
                          meta$LOCATION_ID, "_",
                          meta$TREATMENT, "_TP", 
                          meta$TIMEPOINT, "_conc",
                          meta$CONCENTRATION, "_rep",
                          meta$REPLICATE)




if(!all(rownames(meta) %in% colnames(count)) |
   !all(colnames(count) %in% rownames(meta))){
  warning("No identical sample names (meta & counts)")
  print(paste0("Unmatched samples in count table: ", paste0(c(colnames(count)[which(!colnames(count) %in% rownames(meta))]), collapse = ", ")))
  print(paste0("Unmatched samples in meta: ", paste0(c(rownames(meta)[which(!rownames(meta) %in% colnames(count))]), collapse = ", ")))
} 

count <- count[, rownames(meta)]

## Library size filter ##

meta$LIB_SIZE <- colSums(count, na.rm = TRUE)
meta$ReadDepth <- meta$LIB_SIZE / nrow(count)
summary(meta$LIB_SIZE)
summary(meta$ReadDepth)
meta_LowLibSize <- meta[which(meta$Percentage_mapped < CountThreshold), ]
CELLID_LowLibSize <- as.data.frame(table(Percentage_mapped $CELL_ID))
meta_Filtered <- meta[which(meta$Percentage_mapped > CountThreshold), ] 
count_Filtered <- count[, rownames(meta_Filtered)]

## Check if CountTable contains NAs ##

count_NA <- count_Filtered[rowSums(is.na(count_Filtered)) > 0, colSums(is.na(count_Filtered)) > 0]
if(nrow(count_NA) == 0){
  print("No NAs in counts table")
} else {
  warning("NAs present in counts table")
  print(table(count_NA)); print(dim(count_NA))
}

count_Filtered <- na.omit(count_Filtered)

## Excluding specific samples for further analysis

count_Filtered <- count_Filtered[, rownames(meta_Filtered)]

## Number of genes with at least one zero count among samples

NrGenes_zero <- sum(apply(count_Filtered, MARGIN = 1, function(x) sum(x == 0)) > 0)
message(paste0(NrGenes_zero, "/", nrow(count_Filtered), " genes have at least one zero count among samples and will not be used for size factor calculation"))

##correlation check rawcount###
require(PerformanceAnalytics)
require(Hmisc)


meta_filetered_whole <- meta_Filtered[meta_Filtered$LOCATION_ID == "Whole",]
meta_filetered_CPT <- meta_Filtered[meta_Filtered$LOCATION_ID == "CPT",]
meta_filetered_PPT <- meta_Filtered[meta_Filtered$LOCATION_ID == "PPT",]
meta_filetered_GM <- meta_Filtered[meta_Filtered$LOCATION_ID == "GM",]


count_Filtered_whole <- count[, rownames(meta_filetered_whole)]

for (i in 1 : length(unique(meta_filetered_whole$CONCENTRATION))){
  CON <- unique(meta_filetered_whole$CONCENTRATION)[i]
  
  meta_filtered_sub <- meta_filetered_whole[meta_filetered_whole$CONCENTRATION == CON, ]
  count_filtered_sub <- count_Filtered_whole[, rownames(meta_filtered_sub)]
  
  res2 <- rcorr(as.matrix(count_filtered_sub))
  rrSquare <- as.data.frame(res2[["r"]])
  
  colBreaks <- seq(-1,1, length.out = 100)
  
  my.colors  <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
  
  pdf(paste0(file = "//VUW/PerSonal$/HomeS/W/wijayalS/DeSktop/ExperimentS/2020/20200501_LCM_RealExperiment/Correlation_beforeMergin/correlation_matrix_replicate_",CON,"_whole_.pdf"), height =15, width = 15)
  print(pheatmap(rrSquare, breaks = colBreaks, color = my.colors))
  dev.off()
}



count_Filtered_GM <- count[, rownames(meta_filetered_GM)]

for (i in 1 : length(unique(meta_filetered_GM$CONCENTRATION))){
  CON <- unique(meta_filetered_whole$CONCENTRATION)[i]
  
  meta_filtered_sub <- meta_filetered_GM[meta_filetered_GM$CONCENTRATION == CON, ]
  count_filtered_sub <- count_Filtered_GM[, rownames(meta_filtered_sub)]
  
  res2 <- rcorr(as.matrix(count_filtered_sub))
  rrSquare <- as.data.frame(res2[["r"]])
  
  colBreaks <- seq(-1,1, length.out = 100)
  
  my.colors  <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
  
  pdf(paste0(file = "//VUW/PerSonal$/HomeS/W/wijayalS/DeSktop/ExperimentS/2020/20200501_LCM_RealExperiment/Correlation_beforeMergin/correlation_matrix_replicate_",CON,"_GM_.pdf"), height =15, width = 15)
  print(pheatmap(rrSquare, breaks = colBreaks, color = my.colors))
  dev.off()
}


count_Filtered_CPT <- count[, rownames(meta_filetered_CPT)]

for (i in 1 : length(unique(meta_filetered_CPT$CONCENTRATION))){
  CON <- unique(meta_filetered_CPT$CONCENTRATION)[i]
  
  meta_filtered_sub <- meta_filetered_CPT[meta_filetered_CPT$CONCENTRATION == CON, ]
  count_filtered_sub <- count_Filtered_CPT[, rownames(meta_filtered_sub)]
  
  res2 <- rcorr(as.matrix(count_filtered_sub))
  rrSquare <- as.data.frame(res2[["r"]])
  
  colBreaks <- seq(-1,1, length.out = 100)
  
  my.colors  <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
  
  pdf(paste0(file = "//VUW/PerSonal$/HomeS/W/wijayalS/DeSktop/ExperimentS/2020/20200501_LCM_RealExperiment/Correlation_beforeMergin/correlation_matrix_replicate_",CON,"_CPT_.pdf"), height =15, width = 15)
  print(pheatmap(rrSquare, breaks = colBreaks, color = my.colors))
  dev.off()
}

count_Filtered_PPT <- count[, rownames(meta_filetered_PPT)]


for (i in 1 : length(unique(meta_filetered_PPT$CONCENTRATION))){
  CON <- unique(meta_filetered_PPT$CONCENTRATION)[i]
  
  meta_filtered_sub <- meta_filetered_PPT[meta_filetered_PPT$CONCENTRATION == CON, ]
  count_filtered_sub <- count_Filtered_PPT[, rownames(meta_filtered_sub)]
  
  res2 <- rcorr(as.matrix(count_filtered_sub))
  rrSquare <- as.data.frame(res2[["r"]])
  
  colBreaks <- seq(-1,1, length.out = 100)
  
  my.colors  <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
  
  pdf(paste0(file = "//VUW/PerSonal$/HomeS/W/wijayalS/DeSktop/ExperimentS/2020/20200501_LCM_RealExperiment/Correlation_beforeMergin/correlation_matrix_replicate_",CON,"_PPT_.pdf"), height =15, width = 15)
  print(pheatmap(rrSquare, breaks = colBreaks, color = my.colors))
  dev.off()
}

###diStribution###

count_Filtered_LCM <- count_Filtered[, !grepl("Slide", colnameS(count_Filtered))]
meta_Filtered_LCM <- meta_Filtered[!grepl("Slide", meta_Filtered$SAMPLE_ID),]
meta_Filtered_clean <- meta_Filtered[meta_Filtered$SAMPLE_ID %in% colnames(Count_Filtered_Clean),]

CPM <- apply(Count_Filtered_Clean, 2, function(x) (x/sum(x))*1000000)

dds <- DESeqDataSetFromMatrix(countData = Count_Filtered_Clean,
                              colData = meta_Filtered_clean,
                              design = ~ meanID)

if (NormMethod == "DESeq2"){ 
  dds <- estimateSizeFactors(dds) 
} else if (NormMethod == "CPM"){
  sizeFactors(dds) <- colSums(Count_Filtered_Clean)/1000000
} else { 
  warning("Incorrect input for 'NormMethod', please double check above!")
}

Norm <- as.data.frame(counts(dds, normalized = TRUE)) #Normalized by Size factorS
log2Norm <- log2(Norm + 1) #Log2 normalization of Size factor corrected countS

# CountS & normalized countS diStribution amoung CELL IDS #
count_long <- melt(as.matrix(Count_Filtered_Clean))
colnames(count_long) <- c("Probe_ID", "SAMPLE_ID", "COUNTS")
count_long <- left_join(count_long, meta_Filtered[,c("SAMPLE_ID", "SAMPLENAME", "CELL_ID", "EXP_ID", "TIMEPOINT", "TREATMENT", "CONCENTRATION", "REPLICATE", "meanID", "LOCATION_ID")])

log2Norm_long <- melt(as.matrix(log2Norm))
colnames(log2Norm_long) <- c("Probe_ID", "SAMPLE_ID", "COUNTS")
log2Norm_long <- left_join(log2Norm_long, meta_Filtered[,c("SAMPLE_ID", "SAMPLENAME", "CELL_ID", "EXP_ID", "TIMEPOINT", "TREATMENT", "CONCENTRATION", "REPLICATE", "meanID", "LOCATION_ID")])

# Counts & normalized counts distribution amoung CELL IDs #
for(i in 1:length(unique(count_long$CONCENTRATION))){
  DOSE <- unique(count_long$CONCENTRATION)[i] 
  count_long_sub <- count_long[count_long$CONCENTRATION == DOSE,]
  
  count_long_sub$VariableQCs = count_long_sub[,which(names(count_long_sub)==VariableQCs)]
  pdf(paste0(outputDir_GR, "Distribution-Counts_aftercorrCleanup_version2_", VariableQCs,"_", DOSE,".pdf"), width = 6, height = 4)
  print(ggplot(count_long_sub, aes(x = reorder(VariableQCs, COUNTS, FUN = median), y = COUNTS+1)) +
          geom_boxplot(size = 0.3, outlier.size = 0.5) +
          scale_y_log10(limits = c(1, max(count_long_sub$COUNTS))) +
          theme_classic() +
          theme(plot.title = element_text(size=14, face="bold", vjust = 2, hjust = 0.5), 
                axis.title.x = element_text(size=12, vjust = 0.25),
                axis.title.y = element_text(size=12, vjust = 1),
                axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
                axis.text.y = element_text(size=12)) +
          ggtitle("Distribution counts") + ylab('counts') + xlab(VariableQCs))
  dev.off()
  
}

for(i in 1:length(unique(log2Norm_long$CONCENTRATION))){
  DOSE <- unique(log2Norm_long$CONCENTRATION)[i] 
  log2Norm_long_sub <- log2Norm_long[log2Norm_long$CONCENTRATION == DOSE,]
  
  
  log2Norm_long_sub$VariableQCs = log2Norm_long_sub[,which(names(log2Norm_long_sub)==VariableQCs)] 
  pdf(paste0(outputDir_GR, "Distribution-CountsNorm_AfterCorrCleanup_version2_", DOSE,".pdf"), width = 6, height = 4)
  
  print(ggplot(log2Norm_long_sub, aes(x = reorder(VariableQCs, COUNTS, FUN = median), y = COUNTS)) +
          geom_boxplot(size = 0.3, outlier.size = 0.5) +
          theme_classic() +
          theme(plot.title = element_text(size=14, face="bold", vjust = 2, hjust = 0.5), 
                axis.title.x = element_text(size=12, vjust = 0.25),
                axis.title.y = element_text(size=12, vjust = 1),
                axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
                axis.text.y = element_text(size=12)) +
          ggtitle("Distribution Normalized counts") + ylab('Normalized counts') + xlab(VariableQCs))
  dev.off()
  
}

##Eliminating the outlierS baSed on the correlation analySiS##

Eliminated_samples <- c("S_2001_1_PPT", "S_2007_3_PPT", "S_2003_3_PPT", "S_4005_1_PPT", "S_4001_3_PPT", "S_5001_4_PPT", "S_2005_4_CPT", "S_2003_4_CPT",
                        "S_4007_3_CPT", "S_5001_3_CPT", "S_5003_1_GM", "S_1003_2_GM","S_1005_2_GM", "S_2005_1_GM", "S_4003_1_GM")

Count_Filtered_Clean <- t(count_Filtered)
Count_Filtered_Clean <- Count_Filtered_Clean[!rownames(Count_Filtered_Clean) %in% Eliminated_samples,]
Count_Filtered_Clean <- t(Count_Filtered_Clean)
Count_Filtered_Clean <- as.data.frame(Count_Filtered_Clean)

###summing the count from technical replicates###
##CELL_ID == Technical_replicate##

count_long <- melt(as.matrix(Count_Filtered_Clean))
colnames(count_long) <- c("Probe_ID", "SAMPLE_ID", "COUNTS")
count_long <- left_join(count_long, meta[,c("SAMPLE_ID", "SAMPLENAME", "CELL_ID", "EXP_ID", "TIMEPOINT", "TREATMENT", "CONCENTRATION", "REPLICATE", "meanID", "LOCATION_ID")])

count_long$sampleID <- count_long$SAMPLE_ID
count_long$sampleID <-  sub("(_[^_]+)_.*", "\\1", count_long$sampleID)
count_long$sampleID <- paste0(count_long$sampleID,"_", count_long$LOCATION_ID, sep = "")
count_long$sampleID <- gsub("Whole", "slide", count_long$sampleID)
count_long$SAMPLE_ID <- NULL
count_long$SAMPLENAME <- NULL
count_long <- aggregate(COUNTS ~ Probe_ID + sampleID + EXP_ID + TIMEPOINT + TREATMENT + CONCENTRATION + REPLICATE + meanID + LOCATION_ID, data = count_long, sum)
count_Filtered <- count_long[,c("Probe_ID", "sampleID", "COUNTS")]
count_Filtered <- spread(count_Filtered,sampleID, COUNTS)
rownames(count_Filtered) <- count_Filtered$Probe_ID
count_Filtered$Probe_ID <- NULL

names(count_long)[2] <- "SAMPLE_ID"
count_long$Probe_ID <- NULL
count_long$COUNTS <- NULL
meta_Filtered <- count_long[!duplicated(count_long$SAMPLE_ID),]
rownames(meta_Filtered) <- meta_Filtered$SAMPLE_ID

idx <- sapply(rownames(meta_Filtered), function(x) {
  which(colnames(count_Filtered) == x)
})
count_Filtered <- count_Filtered[,idx]

##transforming the data##

count_Filtered <- count_Filtered+1


###base Mean calculation##
CPM <- apply(count_Filtered, 2, function(x) (x/sum(x))*1000000)

dds <- DESeqDataSetFromMatrix(countData = count_Filtered,
                              colData = meta_Filtered,
                              design = ~ meanID)
#par(mar=c(8,5,2,2))
#boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

if (NormMethod == "DESeq2"){ 
  dds <- estimateSizeFactors(dds) 
} else if (NormMethod == "CPM"){
  sizeFactors(dds) <- colSums(count_Filtered)/1000000
} else { 
  warning("Incorrect input for 'NormMethod', please double check above!")
}

baseMean <-  as.data.frame(rowMeans(counts(dds, normalized=TRUE)))
summary(baseMean)
baseMean$Probe_ID <- rownames(baseMean)
Low_Expr_Gene <- baseMean[baseMean$`rowMeans(counts(dds, normalized = TRUE))` < 0.355,]

count_Filtered$probe <- rownames(count_Filtered)
count_Filtered <- count_Filtered[!count_Filtered$probe %in% Low_Expr_Gene$Probe_ID,]
count_Filtered$probe <- NULL


###Normalized count##
conc <- colnames(count_Filtered)
conc <- gsub("[0-9]", "", sapply(strsplit(conc, '_', perl = T), '[', 3))
conc <- unique(conc)

paths <- conc

calc_norm_gf <- function(paths){
  toMatch <- c(paths)
  count_sub_com <- count_Filtered[, grepl(paste0(toMatch), names(count_Filtered))]
  
  metafiltered_sub <- meta_Filtered[meta_Filtered$SAMPLE_ID %in% names(count_sub_com),]
  
  dds <- DESeqDataSetFromMatrix(countData = count_sub_com,
                                colData = metafiltered_sub,
                                design = ~ meanID)
  
  if (NormMethod == "DESeq2"){ 
    dds <- estimateSizeFactors(dds) 
  } else if (NormMethod == "CPM"){
    sizeFactors(dds) <- colSums(count_sub_com)/1000000
  } else { 
    warning("Incorrect input for 'NormMethod', please double check above!")
  }
  
  
  dds_compound <- DESeq(dds)
}

dds <- lapply(paths, calc_norm_gf)


Norm_1 <- as.data.frame(counts(dds[[1]], normalized = TRUE)) #Normalized by size factors
Norm_2 <- as.data.frame(counts(dds[[2]], normalized = TRUE)) #Normalized by size factors
Norm_3 <- as.data.frame(counts(dds[[3]], normalized = TRUE)) #Normalized by size factors
Norm_4 <- as.data.frame(counts(dds[[4]], normalized = TRUE)) #Normalized by size factors
Norm <- cbind(Norm_1, Norm_2, Norm_3, Norm_4)
log2Norm <- log2(Norm + 1) #Log2 normalization of size factor corrected counts

###retrieving the log2FC###
meta_filetered_whole <- meta_Filtered[meta_Filtered$LOCATION_ID == "Whole",]
ControlSample_whole <- "EUT_076_LCM_Whole_CIS_conc0_TP24"

meta_filetered_CPT <- meta_Filtered[meta_Filtered$LOCATION_ID == "CPT",]
ControlSample_CPT <- "EUT_076_LCM_CPT_CIS_conc0_TP24"

meta_filetered_PPT <- meta_Filtered[meta_Filtered$LOCATION_ID == "PPT",]
ControlSample_PPT <- "EUT_076_LCM_PPT_CIS_conc0_TP24"

meta_filetered_GM <- meta_Filtered[meta_Filtered$LOCATION_ID == "GM",]
ControlSample_GM <- "EUT_076_LCM_GM_CIS_conc0_TP24"


samples_whole <- unique(meta_filetered_whole$meanID[which(!meta_filetered_whole$meanID %in% ControlSample_whole)])
samples_CPT <- unique(meta_filetered_CPT$meanID[which(!meta_filetered_CPT$meanID %in% ControlSample_CPT)])
samples_PPT <- unique(meta_filetered_PPT$meanID[which(!meta_filetered_PPT$meanID %in% ControlSample_PPT)])
samples_GM <- unique(meta_filetered_GM$meanID[which(!meta_filetered_GM$meanID %in% ControlSample_GM)])

data_list <- c("samples_whole","samples_CPT","samples_PPT","samples_GM")

paths <- data_list

calc_summaries <- function( paths ){
  
  if (paths == "samples_whole") {
    df_DEGs <- c() 
    for(y in 1:length(samples_whole)){
      sample <- samples_whole[y]
      print(paste(y, sample, ControlSample_whole, sep = ", "))
      rslt_tmp <- as.data.frame(results(dds[[4]], contrast = c("meanID", sample, ControlSample_whole)))
      rslt_tmp$meanID <- sample
      rslt_tmp$ref_meanID <- ControlSample_whole
      rslt_tmp$Probe_ID <- rownames(rslt_tmp)
      rownames(rslt_tmp) <- c()
      df_DEGs <- rbind(df_DEGs, rslt_tmp)
    }
    return(df_DEGs)
  }else if (paths == "samples_CPT"){
    df_DEGs <- c()
    for(y in 1:length(samples_CPT)){
      sample <- samples_CPT[y]
      print(paste(y, sample, ControlSample_CPT, sep = ", "))
      rslt_tmp <- as.data.frame(results(dds[[1]], contrast = c("meanID", sample, ControlSample_CPT)))
      rslt_tmp$meanID <- sample
      rslt_tmp$ref_meanID <- ControlSample_CPT
      rslt_tmp$Probe_ID <- rownames(rslt_tmp)
      rownames(rslt_tmp) <- c()
      df_DEGs <- rbind(df_DEGs, rslt_tmp)
    }
    return(df_DEGs)
  } else if (paths == "samples_PPT"){
    df_DEGs <- c()  
    for(y in 1:length(samples_PPT)){
      sample <- samples_PPT[y]
      print(paste(y, sample, ControlSample_PPT, sep = ", "))
      rslt_tmp <- as.data.frame(results(dds[[3]], contrast = c("meanID", sample, ControlSample_PPT)))
      rslt_tmp$meanID <- sample
      rslt_tmp$ref_meanID <- ControlSample_PPT
      rslt_tmp$Probe_ID <- rownames(rslt_tmp)
      rownames(rslt_tmp) <- c()
      df_DEGs <- rbind(df_DEGs, rslt_tmp)
    }
    return(df_DEGs) 
  } else if (paths == "samples_GM"){
    df_DEGs <- c()
    for(y in 1:length(samples_GM)){
      sample <- samples_GM[y]
      print(paste(y, sample, ControlSample_GM, sep = ", "))
      rslt_tmp <- as.data.frame(results(dds[[2]], contrast = c("meanID", sample, ControlSample_GM)))
      rslt_tmp$meanID <- sample
      rslt_tmp$ref_meanID <- ControlSample_GM
      rslt_tmp$Probe_ID <- rownames(rslt_tmp)
      rownames(rslt_tmp) <- c()
      df_DEGs <- rbind(df_DEGs, rslt_tmp)
    }
    return(df_DEGs)
  } else { stop ("data matrix is not found")}
}

summary_lists <- lapply(paths, calc_summaries)
df_DEGs <- do.call(rbind,summary_lists)


df_DEGs <- unique(left_join(df_DEGs, meta_Filtered[,c("EXP_ID","TIMEPOINT","TREATMENT","CONCENTRATION","meanID","LOCATION_ID" , "REPLICATE")], by = "meanID"))
df_DEGs <- separate(df_DEGs, c("Probe_ID"), into = c("GeneSymbol", "ProbeNr"), sep = "_", remove = FALSE)

count_long <- melt(as.matrix(count_Filtered))
colnames(count_long) <- c("Probe_ID", "SAMPLE_ID", "COUNTS")
count_long <- left_join(count_long, meta_Filtered[,c("SAMPLE_ID", "EXP_ID", "TIMEPOINT", "TREATMENT", "CONCENTRATION", "REPLICATE", "meanID", "LOCATION_ID")])

Norm_long <- melt(as.matrix(Norm))
colnames(Norm_long) <- c("Probe_ID", "SAMPLE_ID", "COUNTS")
Norm_long <- left_join(Norm_long, meta_Filtered[,c("SAMPLE_ID", "EXP_ID", "TIMEPOINT", "TREATMENT", "CONCENTRATION", "REPLICATE", "meanID", "LOCATION_ID")])

log2Norm_long <- melt(as.matrix(log2Norm))
colnames(log2Norm_long) <- c("Probe_ID", "SAMPLE_ID", "COUNTS")
log2Norm_long <- left_join(log2Norm_long, meta_Filtered[,c("SAMPLE_ID", "EXP_ID", "TIMEPOINT", "TREATMENT", "CONCENTRATION", "REPLICATE", "meanID", "LOCATION_ID")])

meta_Filtered$TIMEPOINT <- as.factor(meta_Filtered$TIMEPOINT)
meta_Filtered$REPLICATE <- as.factor(meta_Filtered$REPLICATE)

# Library size distribution #
for(i in 1:length(unique(meta$CONCENTRATION))){
 DOSE <- unique(meta$CONCENTRATION)[i] 
meta_sub <- meta[meta$CONCENTRATION == DOSE,]
 
meta_sub$VariableQCs = meta_sub[,which(names(meta_sub)==VariableQCs)]
pdf(paste0(outputDir_GR, "LibSize_", VariableQCs,"_", DOSE,".pdf"), width = 6, height = 4)
print(ggplot(meta_sub, aes(x = reorder(VariableQCs, LIB_SIZE, FUN = median), y = LIB_SIZE)) +
        geom_boxplot(size = 0.3, outlier.size = 0.5) +
        scale_y_log10(limits = c(1, max(meta_sub$LIB_SIZE))) +
        theme_classic() +
        theme(plot.title = element_text(size=14, face="bold", vjust = 2, hjust = 0.5), 
              axis.title.x = element_text(size=12, vjust = 0.25),
              axis.title.y = element_text(size=12, vjust = 1),
              axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
              axis.text.y = element_text(size=12)) +
        ggtitle("Library size distribution") + ylab('Library size') + xlab(VariableQCs))
dev.off()

}

# Counts & normalized counts distribution for each concentration #

for(i in 1:length(unique(count_long$CONCENTRATION))){
  DOSE <- unique(count_long$CONCENTRATION)[i] 
  count_long_sub <- count_long[count_long$CONCENTRATION == DOSE,]
  
count_long_sub$VariableQCs = count_long_sub[,which(names(count_long_sub)==VariableQCs)]
pdf(paste0(outputDir_GR, "Distribution-Counts_CleanVersion6_TechCombine", VariableQCs,"_", DOSE,".pdf"), width = 6, height = 4)
print(ggplot(count_long_sub, aes(x = reorder(VariableQCs, COUNTS, FUN = median), y = COUNTS+1)) +
        geom_boxplot(size = 0.3, outlier.size = 0.5) +
        scale_y_log10(limits = c(1, max(count_long_sub$COUNTS))) +
        theme_classic() +
        theme(plot.title = element_text(size=14, face="bold", vjust = 2, hjust = 0.5), 
              axis.title.x = element_text(size=12, vjust = 0.25),
              axis.title.y = element_text(size=12, vjust = 1),
              axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
              axis.text.y = element_text(size=12)) +
        ggtitle("Distribution counts") + ylab('counts') + xlab(VariableQCs))
dev.off()

}

for(i in 1:length(unique(log2Norm_long$CONCENTRATION))){
  DOSE <- unique(log2Norm_long$CONCENTRATION)[i] 
  log2Norm_long_sub <- log2Norm_long[log2Norm_long$CONCENTRATION == DOSE,]
  
  
log2Norm_long_sub$VariableQCs = log2Norm_long_sub[,which(names(log2Norm_long_sub)==VariableQCs)] 
pdf(paste0(outputDir_GR, "Distribution-CountsNorm_CleanVersion6_TechCombine_", DOSE,".pdf"), width = 6, height = 4)


print(ggplot(log2Norm_long_sub, aes(x = reorder(VariableQCs, COUNTS, FUN = median), y = COUNTS)) +
        geom_boxplot(size = 0.3, outlier.size = 0.5) +
        theme_classic() +
        theme(plot.title = element_text(size=14, face="bold", vjust = 2, hjust = 0.5), 
              axis.title.x = element_text(size=12, vjust = 0.25),
              axis.title.y = element_text(size=12, vjust = 1),
              axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
              axis.text.y = element_text(size=12)) +
        ggtitle("Distribution Normalized counts") + ylab('Normalized counts') + xlab(VariableQCs))
dev.off()

}


# Replicate correlation to meanID#
for(i in 1:length(unique(meta_Filtered$CONCENTRATION))){
  DOSE <- unique(meta_Filtered$CONCENTRATION)[i] 

  meta_sub <- meta_Filtered[meta_Filtered$CONCENTRATION == DOSE,]

    
repCors <- lapply(unique(meta_sub[, "meanID"]), function(expID) {
#expID <- "EUT_076_LCM_Whole_CIS_conc0_TP24"
  replicateIDs <- rownames(meta_sub)[which(meta_sub[, "meanID"] == expID)]
  if(length(replicateIDs)>1){
    experimentMean <- rowMeans(log2Norm[, replicateIDs])
    apply(log2Norm[, replicateIDs], 2, function(.expID) cor(.expID, experimentMean))
  }
})

repCors <- lapply(repCors,function(y){cbind("pearsonR" = y, "SAMPLE_ID" = names(y))})
repCors <- as.data.frame(do.call(rbind, repCors))
repCors$pearsonR <- as.numeric(as.character(repCors$pearsonR))
repCors$clr <- repCors$pearsonR > 0.8
repCors <- repCors[which(!is.na(repCors$pearsonR)),]
repCors <- left_join(repCors, meta_sub)

repCors$VariableQCs = repCors[,which(names(repCors)==VariableQCs)]
pdf(paste0(outputDir_GR, "RepCor_PearsonR_CorrClean_ver6", VariableQCs,"_",DOSE ,".pdf"), width = 10, height = 4)
print(ggplot(repCors, aes(x = reorder(VariableQCs, pearsonR, FUN = min), y = pearsonR)) +
  geom_hline(yintercept = 0.8, colour = 'grey75', size = 0.75) +
  geom_point(aes(color = clr)) +
  theme_classic() +
  theme(plot.title = element_text(size=14, face="bold", vjust = 2, hjust = 0.5), 
        axis.title.x = element_text(size=12, vjust = 0.25),
        axis.title.y = element_text(size=12, vjust = 1),
        axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=12))+
  ggtitle("Replicate correlation") + ylab('PearsonR') + xlab(VariableQCs))
dev.off()

}

repCors_lowPR <- filter(repCors, clr == FALSE)


##Number of differential expressed genes plotting##

df_nrDEGs <- ddply(df_DEGs, .(meanID, ref_meanID, CELL_ID, TREATMENT, TIMEPOINT, EXP_ID), summarize, 
                   NrDEGs = sum(padj < Threshold_padj , na.rm = TRUE),
                   NrDEGsFC = sum(padj < Threshold_padj & abs(log2FoldChange) > Threshold_log2FC , na.rm = TRUE))



df_DEGs_FCmatrix <- dcast(Count_batch2, Probe_ID ~ meanID, value.var = "log2FoldChange")
rownames(df_DEGs_FCmatrix) <- df_DEGs_FCmatrix[,1]
df_DEGs_FCmatrix <- df_DEGs_FCmatrix[,-1]

meta_Filtered_meanID <- unique(meta_Filtered[, c("CELL_ID", "TREATMENT", "CONCENTRATION", "TIMEPOINT", "EXP_ID", "meanID")])
rownames(meta_Filtered_meanID) <- meta_Filtered_meanID$meanID

## Separate text data frames for each sample for Ingenuity Pathway Analysis ##
## Distribution baseMean ##

for(i in 1:length(unique(df_DEGs_filtered$CONCENTRATION))){
  DOSE <- unique(df_DEGs_filtered$CONCENTRATION)[i] 
  count_long_sub <- df_DEGs_filtered[df_DEGs_filtered$CONCENTRATION == DOSE,]

  count_long_sub$VariableQCs = count_long_sub[,which(names(count_long_sub)==VariableQCs)]
pdf(paste0(outputDir_GR, "Histogram_baseMean_", ControlSample,"_",DOSE, ".pdf"), width = 15, height = 3)
print(ggplot(count_long_sub, aes(x = baseMean)) +
  geom_histogram(binwidth = 0.1) +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0), trans = "log10") +
  theme(plot.title = element_text(size=14, face="bold", vjust = 2, hjust = 0.5), 
        axis.title.x = element_text(size=12, vjust = 0.25),
        axis.title.y = element_text(size=12, vjust = 1),
        axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=8)) +
  ggtitle("Histogram baseMean") + ylab('Count') + xlab('Log10 baseMean') +
  facet_grid(.~ VariableQCs))
dev.off()

}

## Distribution of log2FC and p-values ##
for(i in 1:length(unique(df_DEGs_filtered$CONCENTRATION))){
  DOSE <- unique(df_DEGs_filtered$CONCENTRATION)[i] 
  count_long_sub <- df_DEGs_filtered[df_DEGs_filtered$CONCENTRATION == DOSE,]

pdf(paste0(outputDir_GR, "Log2FC_vs_baseMean_CorrClean_ver6", DOSE, ".pdf"), width = 15, height = 3)
print(ggplot(count_long_sub, aes(x = baseMean, y  = log2FoldChange)) +
  geom_point(aes(color = padj < Threshold_padj), size = 0.2, alpha = 0.3) +
  geom_hline(aes(yintercept = Threshold_log2FC), color = "dodgerblue", size = 0.6, alpha = 0.5) +
  geom_hline(aes(yintercept = -Threshold_log2FC), color = "dodgerblue", size = 0.6, alpha = 0.5) +
  scale_color_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "red", `NA` = "darkgrey")) +
  scale_x_continuous(expand = c(0, 0), trans = "log10") +
  scale_y_continuous(expand = c(0, 0), limits = c(-20, 20)) +
  theme_bw() +
  theme(plot.title = element_text(size=14, face="bold", vjust = 2, hjust = 0.5), 
        axis.title.x = element_text(size=12, vjust = 0.25),
        axis.title.y = element_text(size=12, vjust = 1),
        axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=8)) +
  ggtitle("Log2FC vs baseMean") + ylab('Log2FC (vs control)') + xlab('Log10 baseMean') +
  facet_grid(.~ COMPOUNDS))
dev.off()

}

for(i in 1:length(unique(df_DEGs_filtered$CONCENTRATION))){
  DOSE <- unique(df_DEGs_filtered$CONCENTRATION)[i] 
  count_long_sub <- df_DEGs_filtered[df_DEGs_filtered$CONCENTRATION == DOSE,]
  
pdf(paste0(outputDir_GR, "padj_vs_log2FC_CorrClean_ver6", DOSE, ".pdf"), width = 15, height = 3)
print(ggplot(count_long_sub, aes(x = log2FoldChange, y  = -log10(padj))) +
  geom_point(aes(color = padj < Threshold_padj), size = 0.2, alpha = 0.3) +
  geom_vline(aes(xintercept = Threshold_log2FC), color = "dodgerblue", size = 0.6, alpha = 0.5) +
  geom_vline(aes(xintercept = -Threshold_log2FC), color = "dodgerblue", size = 0.6, alpha = 0.5) +
  scale_color_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "red", `NA` = "darkgrey")) +
  scale_x_continuous(expand = c(0, 0), limits = c(-10, 10)) +
  scale_y_continuous(expand = c(0, 0), limits =c(0, 15)) +
  theme_bw() +
  theme(plot.title = element_text(size=14, face="bold", vjust = 2, hjust = 0.5), 
        axis.title.x = element_text(size=12, vjust = 0.25),
        axis.title.y = element_text(size=12, vjust = 1),
        axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=8)) +
  ggtitle("Adjusted p-value vs Log2FC") + ylab('-log10 adjusted p-value') + xlab('log2FC') +
  facet_grid(.~ LOCATION_ID))
dev.off()
}


###PCA plots##

  
df_log2Norm_pca <- PCADF(log2Norm)
PCAplots(df_log2Norm_pca, meta_Filtered, colourVar = "Cluster", Clusts = FALSE, FileID = "Log2-AllSamples-AllGenes_th_new", width = 10, height = 8, sizeP = 0.2, Dataclass = "log2", a = 0.7)

pcaCoV(df_log2Norm_pca, meta_Filtered, filename = "Log2-AllSamples-AllGenes")


## Saving data ##	

setwd(inputDir)

write.table(meta_Filtered, file = 'meta_Filtered_.txt', sep = '\t')
write.table(count, file = 'count.txt', sep = '\t')
write.table(count_Filtered, file = 'count_Filtered.txt', sep = '\t')

write.table(Norm, file = 'count_Filtered_CPMnorm.txt', sep = '\t')
write.table(log2Norm, file = 'count_Filtered_log2CPMnorm.txt', sep = '\t')

write.table(df_DEGs, file = 'DEGsfiltered_vs_correctControl_c.txt', sep = "\t")

