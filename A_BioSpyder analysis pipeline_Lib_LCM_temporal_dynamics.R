

## Data analysis pipeline BioSpyder datasets ##

# A - Library

# Marije Niemeijer & Giulia Callegaro
# 11/07/2018

# Version 1

###############################################

###### Settings ######

options(stringsAsFactors = FALSE)
options(scipen=999)

## Libraries ##

library("gridExtra")
library("stringr")
library("ggplot2")
library("pheatmap")
library("reshape2")
library("RColorBrewer")
library("plyr")
library("dplyr")
library("tidyr")
library("colorspace")
library("scales")
library("data.table")
library("DESeq2") # bioconductor 
library("compare")
library("readxl")
library("PoiClaClu")
library("hexbin")
library("ggalt")
library("vsn")# bioconductor
library("org.Hs.eg.db")# bioconductor
library("AnnotationDbi")# bioconductor

#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")


###### Functions ######

# Function: Loading of csv or xlsx raw count tables supplied by BioSpyder
LoadCountData <- function(FilePath) # Directory of count table files (allows to load multiple count tables at once)
  {
  setwd(FilePath)
  CountData <- list()
  
  for(i in 1:length(dir(FilePath))){
    x <- dir(FilePath)[i]
    if(grepl("xlsx", x)){
      tmp  <- read_excel(x, sheet = 1)
    } else {
      tmp  <- read.csv(x)
    }
    
    CountData[i] <- list(as.data.frame(tmp))
  }  
  
  probeNames <- CountData[[1]][[1]]
  
  for(i in 1:length(CountData)) {
    rownames(CountData[[i]]) <- probeNames
    if (!isTRUE(all.equal(rownames(CountData[[i]]),CountData[[i]][[1]]))) { 
      mismatches <- paste(which(rownames(CountData[[i]]) != CountData[[i]][[1]]), collapse = ",")
      print(paste0("ProbeIDs do not match between count tables: ", i)); print(mismatches)
    }
    
    CountData[[i]] <- CountData[[i]][-1]
  }
  CountData <- do.call(cbind, CountData)
  CountData
}

# Function: Subset count table based on gene list
PathwayProbesCounts <- function(Genelist = "OX_genes.txt", #Variable name or name of txt file of list of genes (containing either probes or genesymbols)
                                CountTable = df,           #Variable name of count table
                                GeneID = "Probes",         #Type of input gene list, fill in either "Probes" or "GeneSymbols"
                                ProbeCol = "Probe_ID",     #Name of column in count table of probes or fill in "rownames" when count table is a matrix
                                match = "Entrez",          #Match genes based on Entrez ID or identical probes, fill in "Entrez" or "Probe"
                                Dir = GeneListDir          #When using txt file as input for list of genes, fill in directory of txt file 
                                ){ 
  
  if(grepl(".txt$", Genelist[1])){
    Genes <- read.table(paste0(Dir, Genelist), sep = "\t", header = TRUE)	
    GenesOfInterest <- Genes[,1]
  } else {
    GenesOfInterest <- Genelist
  }
  
  if(match == "Entrez"){
    if(GeneID == "Probes"){
      GenesOfInterest <- sub("_.*", "", GenesOfInterest)
    }
    GenesOfInterest <- unname(mapIds(org.Hs.eg.db, keys= GenesOfInterest, column="ENTREZID", keytype="ALIAS", multiVals="first"))
    
    if (ProbeCol == "rownames"){
      CountTableProbes <- sub("_.*", "", rownames(CountTable))
    } else {
      CountTableProbes <- sub("_.*", "", CountTable[,ProbeCol])
    }
    CountTableProbes <- unname(mapIds(org.Hs.eg.db, keys= CountTableProbes, column="ENTREZID", keytype="ALIAS", multiVals="first"))
  } 
  
  if(match == "Probe"){
    if(ProbeCol == "rownames"){
      CountTableProbes <- rownames(CountTable)
    } else {
      CountTableProbes <- CountTable[,ProbeCol]
    }
  }
  
  df_sub <- CountTable[which(CountTableProbes %in% GenesOfInterest),]
  return(df_sub)
}

# Function: PCA using prcomp function
PCADF <- function(DF = log2Norm) #Variable name of normalized and transformed count table for PCA
  {
  pca <- prcomp(t(DF))
  return(pca)
}

# Function: Plotting PCA results (PC1-5) 
PCAplots <- function(pca = df_pca,            #Variable of PCA results derived from PCADF function
                     metaDF = meta_Filtered,  #Variable name of meta data table with samples as rownames (identical as column names in count table)
                     colourVar = "TREATMENT", #Fill in variable (meta data column name) for which you want unique colors in PCA plot
                     Clusts = FALSE,          #When TRUE, clusters will be plotted (based on k-means clustering)
                     NrClusts = 3,            #When Clusts = TRUE, specify number of groups to cluster samples
                     FileID = "Allsamples",   #Unique name to save plot
                     width = 10,              #Width of plot
                     height = 8,              #Height of plot
                     sizeP = 0.2,             #size of points in PCA plot
                     Dataclass = "log2",      #Type of data supplied to PCA either log2 normalized counts or log2FC, fill in "log2" or "log2FC"
                     a = 0.7){                #Transparency of dots in plot
  
  if(Dataclass == "log2"){
    Vari <- c("SAMPLE_ID",  "TREATMENT", "REPLICATE", "EXP_ID","TIMEPOINT", "LOCATION_ID", "CONCENTRATION")
    ID <- "SAMPLE_ID"
  }
  
  if(Dataclass == "log2FC"){
    Vari <- c("SAMPLE_ID",  "TREATMENT", "REPLICATE", "EXP_ID","TIMEPOINT", "LOCATION_ID", "CONCENTRATION", "CELL_ID", "percetangeMapped_cat")
    ID <- "meanID"
  }
  
  pca_eigs <- pca$sdev^2
  pca_proportion <- pca_eigs / sum(pca_eigs)
  pca_proportion <- data.frame("variable" = c("PC1", "PC2", "PC3", "PC4", "PC5"), "proportion" = pca_proportion[1:5])
  pca_clust <- melt(kmeans(pca$x[,1:5], centers = NrClusts)$cluster)
  pca_clust[,ID] <- rownames(pca_clust)
  colnames(pca_clust)[1] <- "clust"
  
  pca_df <- data.frame(pca$x[which(grepl(paste(metaDF[,ID], collapse = "|"), rownames(pca$x))), 1:5])
  pca_df <- data.frame(metaDF[rownames(pca_df),], pca_df)
  pca_Long <- melt(pca_df, id.vars = Vari, measure.vars = grep("PC", colnames(pca_df), value = TRUE))
  pca_Long <- left_join(pca_Long, pca_proportion)
  pca_Long$PlotLabels <- paste0(pca_Long$variable, " (Variance: ", 100*round(pca_Long$proportion, digits = 3), "%)")
  pca_Long <- left_join(pca_Long, pca_clust)
  
  pca_GR <- do.call(rbind, lapply(1:5, function(i) {
    pca_Long$x <-   pca_Long$value[pca_Long$variable == paste0('PC', i)]
    pca_Long$xPC <- factor(paste0('PC', i), ordered = TRUE)
    pca_Long$yPC <- factor(pca_Long$variable, ordered = TRUE)
    pca_Long
  })) 
  
  levels(pca_GR$xPC) <- unique(pca_GR$PlotLabels) 
  
  pca_GR$colourVar <- pca_GR[,which(names(pca_GR)==colourVar)]
  
  if(Clusts == TRUE){
    cl <- "Clusts"
    pdf(paste0(outputDir_GR, "PCA_", FileID, "_", colourVar, "_", cl,"_",".pdf"), width = width, height = height)
    print(ggplot(pca_GR, aes(x = x, y = value)) +
            geom_encircle(aes(group = factor(clust), fill = factor(clust)), s_shape = 0.7, expand = 0, color = "#666666", alpha = a) +    
            geom_point(aes(colour = colourVar, shape = LOCATION_ID, size = CONCENTRATION), alpha = a) +
            facet_grid(yPC ~ xPC, scales = 'free') +
            theme_bw() + scale_size_continuous(range = c(0.3, 1))+
            theme(panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank()) +
            labs(x = "", 
                 y = "", 
                 title = paste0("PCA_", FileID, "_", colourVar, "_", cl),
                 colour = paste0(colourVar)))
    dev.off()

  } else { 
    cl <- "NoClusts"
    pdf(paste0(outputDir_GR, "PCA_", FileID, "_", colourVar, "_", cl,"_",".pdf"), width = width, height = height)
    print(ggplot(pca_GR, aes(x = x, y = value)) +
            geom_point(aes(colour = colourVar,shape = LOCATION_ID, size = CONCENTRATION), alpha = a) +
            facet_grid(yPC ~ xPC, scales = 'free') +
            theme_bw() + scale_size_continuous(range = c(0.3, 1))+
            theme(panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank()) +
            labs(x = "", 
                 y = "", 
                 title = paste0("PCA_", FileID, "_", colourVar, "_", cl),
                 colour = paste0(colourVar)))
    dev.off()
  }
}

# Function: Returns dataframe of top genes contributing to PC1-5 (based on prcomp-rotation)
PCAcontrVar <- function(pcaDF = df_pca, #Variable of PCA results derived from PCADF function
                        NrGenes = 25    #Number of top genes in output data frame
                        ){
  tmp <- data.frame(PC1Genes = names(pcaDF$rotation[rev(order(abs(pcaDF$rotation[,1])))[1:NrGenes],1]),
                    PC1rotation = unname(pcaDF$rotation[rev(order(abs(pcaDF$rotation[,1])))[1:NrGenes],1]),
                    PC2Genes = names(pcaDF$rotation[rev(order(abs(pcaDF$rotation[,2])))[1:NrGenes],2]),
                    PC2rotation = unname(pcaDF$rotation[rev(order(abs(pcaDF$rotation[,2])))[1:NrGenes],2]),
                    PC3Genes = names(pcaDF$rotation[rev(order(abs(pcaDF$rotation[,3])))[1:NrGenes],3]),
                    PC3rotation = unname(pcaDF$rotation[rev(order(abs(pcaDF$rotation[,3])))[1:NrGenes],3]),
                    PC4Genes = names(pcaDF$rotation[rev(order(abs(pcaDF$rotation[,4])))[1:NrGenes],4]),
                    PC4rotation = unname(pcaDF$rotation[rev(order(abs(pcaDF$rotation[,4])))[1:NrGenes],4]),
                    PC5Genes = names(pcaDF$rotation[rev(order(abs(pcaDF$rotation[,5])))[1:NrGenes],5]),
                    PC5rotation = unname(pcaDF$rotation[rev(order(abs(pcaDF$rotation[,5])))[1:NrGenes],5]))
  return(tmp)
}

# Function: Heatmap of contribution of variables to PC1-5 (based on Pearson correlation)
pcaCoV <- function(pcaDF = df_pca,         #Variable of PCA results derived from PCADF function
                   metaDF = meta_Filtered, #Variable name of meta data table with samples as rownames (identical as column names in count table)
                   filename = "Allsamples" #Unique name to save heatmap
                   ){
  
  pcaDF <- data.frame(metaDF[rownames(pcaDF$x),], pcaDF$x[, 1:5])
  pcaDF[, which(grepl("character|logical|factor", sapply(pcaDF, class)))] <- apply(pcaDF[,which(grepl("character|logical|factor", sapply(pcaDF, class)))], 2, function(x) as.numeric(as.factor(x)))
  pearsonMatrix <- na.omit(cor(pcaDF, method = 'pearson')[grep('PC', colnames(pcaDF), value = TRUE, invert = TRUE), 
                                                          grep('PC', colnames(pcaDF), value = TRUE)])
  
  pheatmap(pearsonMatrix, main = paste0('Pearson correlation -', filename),
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           breaks = c(seq((-1 * 1), 0, length.out = ceiling(100/2) + 1), 
                      seq(1/100, 1, length.out = floor(100/2))),
           cellwidth = 25, 
           cellheight = 15,
           width = 7,
           file = paste0(outputDir_GR, "Heatmap_PCACoV_", filename, ".pdf"))
}

# Function: Hierarchical clustering of samples and genes based on expression using Pearson correlation 
Heatmap <- function(Data = log2Norm,         #Variable name of normalized-transformed count table or Log2FCs (as matrix: with probes as rownames and samples as column names)
                    Meta = meta_Filtered,    #Variable name of meta data table with samples as rownames (identical as column names in count table)
                    FC = FALSE,              #When TRUE, count table are log2FCs. Otherwise fill in FALSE
                    heightGene = 0.2,        #Height of squares in heatmap
                    widthGenes = 1,          #Width of squares in heatmap
                    FileName = "Allsamples", #Unique name to save heatmap
                    GeneName = FALSE         #If TRUE, gene names are shown in heatmap
                    ) {
  
  Meta <- Meta[colnames(Data),]
  
  if(FC == FALSE){
    breaksList <- seq(min(Data, na.rm = TRUE)*1.1, max(Data, na.rm = TRUE)*1.1, by = 0.01)
    meta_Filtered_ann <- Meta[, c("EXP_ID", "CELL_ID", "REPLICATE", "TREATMENT", "CONCENTRATION", "TIMEPOINT")]
    hm_color <- rev(heat_hcl(length(breaksList), c = c(80, 30), l = c(30, 90), power = c(1/5, 1.3)))
  } else {
    meta_Filtered_ann <- Meta[, c("EXP_ID", "CELL_ID", "TREATMENT", "CONCENTRATION", "TIMEPOINT")]
    Data[is.na(Data)] <- 0
    breaksList <- c(seq(min(Data), -0.5, length.out=50), seq(-0.5, 0.5,length.out=50)[-1],
                    seq(0.5, max(Data),length.out=50)[-1])
    hm_color <- rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(150))
  }
  
  meta_Filtered_ann[] <- lapply(meta_Filtered_ann, factor)
  df_hm <- Data[,rownames(meta_Filtered_ann)]
  
  df_hm <- df_hm[apply(df_hm, MARGIN = 1, FUN = function(x) sd(x, na.rm = TRUE) != 0),]
  df_hm <- df_hm[, apply(df_hm, MARGIN = 2, FUN = function(x) sd(x, na.rm = TRUE) != 0)]
  df_hm <- df_hm[, apply(df_hm, MARGIN = 2, FUN = function(x) all(is.finite(x)))]
  
  if(nrow(df_hm) > 5000){
    stop("Too many genes (rows > 4000) in heatmap, subset your data")
  }
  
  pheatmap(as.matrix(df_hm),
           color = hm_color,
           breaks = breaksList,
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           clustering_method = "complete", 
           annotation_col = meta_Filtered_ann,
           display_numbers = FALSE, 
           fontsize_row = 6,
           fontsize_col = 6,
           fontsize_number = 3, 
           border_color = NA, 
           cluster_cols = TRUE,
           show_colnames = TRUE,
           show_rownames = GeneName,
           cellwidth = widthGenes,
           cellheight = heightGene,
           treeheight_row = 25,
           treeheight_col = 30,
           main = FileName,
           filename = paste0(outputDir_GR, paste("Heatmap_", FileName, ".pdf")))
}


