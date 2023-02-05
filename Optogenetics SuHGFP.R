# clear the environment: 
rm(list = ls())
# clear the console: ctrl+L
# clear all the plots: 
dev.off()

#### Packages for this script ####
install.packages("ggplot2")
install.packages("cowplot")
install.packages("dplyr")
install.packages("stringr")
install.packages("ggsignif")
install.packages("ggpubr")
install.packages("gridExtra")
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyverse)
library(stringr)
library(ggsignif)
library(ggpubr)
library(gridExtra)
#### end ####

#### Data analysis functions ####
## Calculate the standard error for each row of a data frame 
SEMrows <- function(df){
  dft <- as.data.frame(t(df))
  rm(df)
  dft <- dft[-1,]
  dft <- na.omit(dft)
  output <- matrix(ncol = length(colnames(dft)), nrow = 1)
  for (n in 1:length(colnames(dft))) {
    dftn <- dft[[n]]
    output[,n] <- sqrt(var(dftn)/length(dftn))
    rm(dftn)
  }
  output <- as.data.frame(t(output))
  return(output$V1)
}

MeanDF <- function(DataFrame){
  MeanList <- list()
  for (l in 1:length(DataFrame)) {
    df <- DataFrame[[l]]
    dfmeanse <- data.frame(Distance=df[["Distance (um)"]], Mean=rowMeans(df[,-1], na.rm = TRUE))
    dfmeanse$SEM <- SEMrows(df)
    Nnumber <- ncol(df)-1
    dfmeanse$Genotype <- paste(names(DataFrame)[l]," (N = ", Nnumber, ")", sep = "")
    MeanList[[l]] <- dfmeanse
    names(MeanList)[l] <- names(DataFrame[l])
    rm(dfmeanse,Nnumber)
  }
  Output <- dplyr::bind_rows(MeanList, .id = 'Name')
  rm(MeanList)
  return(Output)
}

NonListMeanDF <- function(df, Genotype){
  dfmeanse <- data.frame(Distance=df[["Distance (um)"]], Mean=rowMeans(df[,-1], na.rm = TRUE))
  dfmeanse$SEM <- SEMrows(df[,-1])
  Nnumber <- ncol(df)-1
  dfmeanse$Genotype <- paste(Genotype," (N = ", Nnumber, ")", sep = "")
  rm(Nnumber)
  return(dfmeanse)
}

# Create a list of dataframes, named for each truncation, where the data contained
# is the wild type data collected alongside that truncation only.
SepDates <- function(DataList,GenotypeToSep){
  Dates <- list()
  for (g in 1:length(DataList)) {
    df <- DataList[[g]]
    ColNames <- colnames(df)
    ColNames <- ColNames[-1]
    for (c in 1:length(ColNames)) {
      SplitName <- str_split(ColNames[c]," ")
      SplitName <- SplitName[[1]]
      ColNames[c] <- SplitName[1]
    }
    dfDates <- unique(ColNames)
    Dates[[g]] <- dfDates
    names(Dates)[g] <- names(DataList[g])
    rm(dfDates, ColNames, SplitName)
  }
  
  data <- DataList[[GenotypeToSep]]
  dataD <- as.data.frame(data[,1])
  names(dataD)[names(dataD) == "data[, 1]"] <- names(data)[1]
  
  SubDates <- Dates
  SubDates[[GenotypeToSep]] <- NULL
  
  Output <- list()
  for (c in 1:length(SubDates)) {
    GenotypeDates <- SubDates[[c]]
    dataG <- dataD
    for (d in 1:length(GenotypeDates)) {
      df <- data[,grepl(GenotypeDates[d], names(data))]
      dataG <- cbind(dataG, df)
      rm(df)
    }
    Output[[c]] <- dataG
    names(Output)[c] <- names(SubDates[c])
    rm(dataG,GenotypeDates)
  }
  rm(data, dataD, SubDates, Dates)
  return(Output)
}

# Create a list of dataframes, named for each truncation, where each dataframe
# has the data of the truncation and wild type data from the corresponding days
# organised into the means and SEM for all the bands imaged.
WTvsMutByDate <- function(SeparatedWTdata, AllData){
  OutputList <- list()
  for (l in 1:length(SeparatedWTdata)) {
    Genotype <- names(SeparatedWTdata)[l]
    WT <- NonListMeanDF(SeparatedWTdata[[l]],"Wild Type")
    Mutant <- NonListMeanDF(AllData[[Genotype]],Genotype)
    df <- rbind(WT,Mutant)
    OutputList[[l]] <- df
    names(OutputList)[l] <- names(SeparatedWTdata[l])
    rm(df, WT, Mutant)
  }
  return(OutputList)
}

# This function will take the RelInt values and separate them into a list of
# dateframes for each truncation. The dataframes are named for their mutant
# and contain RelInt values for the truncation and wild type data taken alongside.
RelIntSepDates <- function(RelIntDataFrame, GenotypeToSep) {
  Mutants <- dplyr::filter(RelIntDataFrame, Genotype != GenotypeToSep)
  WildTypes <- dplyr::filter(RelIntDataFrame, Genotype == GenotypeToSep)
  OutputList <- list()
  for (t in 1:length(unique(Mutants$Genotype))) {
    MutantData <- dplyr::filter(RelIntDataFrame, Genotype == unique(Mutants$Genotype)[t])
    StrSplit <- str_split(MutantData$Nucleus, " ")
    Dates <- c()
    for (s in 1:length(StrSplit)) {
      SplitName <- StrSplit[[s]]
      Date <- SplitName[1]
      Dates <- append(Dates,Date)
      rm(Date, SplitName)
    }
    UniqueDates <- unique(Dates)
    rm(Dates)
    df <- WildTypes[grep(UniqueDates[1], WildTypes$Nucleus),]
    for (d in 2:length(UniqueDates)) {
      Nextdf <- WildTypes[grep(UniqueDates[d], WildTypes$Nucleus),]
      df <- rbind(df,Nextdf)
    }
    df <- rbind(df,MutantData)
    OutputList[[t]] <- df
    names(OutputList)[t] <- unique(Mutants$Genotype)[t]
    rm(df)
  }
  return(OutputList)
}
#### end ####

#### Plotting functions ####
## Define a basic theme for the plots
theme_jmt <- function(){theme(
  text = element_text(family = "sans", size = 25), # set default font to be Arial
  plot.title = element_text(hjust=0.5, face = 'plain'), # align title to centre and make bold
  panel.grid.major = element_line(colour = "grey80", size = 0.1), # strip major gridlines
  panel.grid.minor = element_blank(), # strip minor gridlines
  panel.background = element_blank(), # strip background
  axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5),  # tilt the x axis
  axis.line = element_line(colour = "black", size = 1.5), # add axis lines in black
  axis.ticks = element_line(colour = "black", size = 1.5),
  axis.ticks.length=unit(2,"mm")) # add tick marks in black
}

## Where required, labels are the Title, X axis and Y axis labels defined in a character string (in order)
## e.g. Labels <- c(title = "Normalised to Background value", x = "Distance (um)", y = "Normalised fluorescence")
## these are then read into the function.

## Define a function that plots a dataframe of distance, mean, sem, grouping by the IncubationTime
## as a fourth column
band_plot <- function(Data, Labels){
  ggplot(Data, aes(Distance,Mean,group = Genotype,col=Genotype,shape=Genotype,ymin=Mean-SEM,ymax=Mean+SEM)) + theme_jmt() +
    geom_ribbon(linetype = "dashed", fill = "grey90",alpha = 0.5, size = 0.25) + geom_line(size = 0.5) + geom_point() +
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]]) +
    scale_color_brewer(palette = "Set2")
}

## Define a function to plot the relative intensities of the bands
RelInt_boxviolin <- function(Data, Labels){
  ggplot(Data, aes(Genotype,RelativeBandIntensity,col=Genotype,shape=Genotype)) + theme_jmt() +
    geom_violin() + geom_boxplot(width=0.1, show.legend = FALSE,colour = "Grey60",alpha = 0.5, outlier.shape = NA) + geom_jitter() +
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]]) + 
    scale_color_brewer(palette = "Set2")
}

RelInt_crossbar <- function(Data, Labels){
  ggplot(Data,aes(Genotype,RelativeBandIntensity,col=Genotype,shape=Genotype)) + theme_jmt() +
    geom_jitter(position = position_jitterdodge()) +
    stat_summary(fun="mean", geom="crossbar", width=0.5) + theme(legend.key=element_blank())  + 
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]]) + 
    scale_color_brewer(palette = "Set2")
}

#### end ####

# Set location of data and where plots will be saved to:
mainDir <- "/Users/jonathantownson/Documents/PhD/Images/Analysed image data/Optogenetics data/newBLITz in light incubator/SuHGFP with locus tag"
setwd(mainDir)
# Load the background values
DF <- read.csv("/Users/jonathantownson/Documents/PhD/Images/Analysed image data/Optogenetics data/newBLITz in light incubator/SuHGFP with locus tag/SuHGFP background values.csv", head = TRUE, sep = ",")
MamBackgroundValues <- data.frame(Nucleus=DF[,1], MamBG=rowMeans(DF[,-1]))
rm(DF)
# Create the path where the data is stored
DataPath <- "/Users/jonathantownson/Documents/PhD/Images/Analysed image data/Optogenetics data/newBLITz in light incubator/SuHGFP with locus tag/SuH"

#### Read in and normalise the data ####
# Create the folder names for the data
FolderNames <- list.dirs(DataPath, full.names = FALSE, recursive = FALSE)
# Read in all the csv files for each condition
MamRawData <- list()
for (f in 1:length(FolderNames)){
  Path <- paste(DataPath,FolderNames[f], sep ="/")
  Files <- list.files(Path, pattern = "*.csv")
  FilePaths <- list()
  for (i in 1:length(Files)){
    FilePaths[[i]] <- paste(Path,Files[i], sep ="/")
  } 
  rm(Path)
  rm(Files)
  df <- read.csv(FilePaths[[1]])
  names(df)[names(df) == 'Y'] <- str_sub(basename(FilePaths[[1]]), end=-9)
  for (i in 2:length(FilePaths)){
    Nextdf <- read.csv(FilePaths[[i]])
    df <- merge(df,Nextdf, by="X", all = TRUE)
    names(df)[names(df) == 'Y'] <- str_sub(basename(FilePaths[[i]]), end=-9)
    rm(Nextdf)
  }
  names(df)[names(df) == 'X'] <- "Distance (um)"
  MamRawData[[f]] <- df
  rm(df)
  rm(FilePaths)
  names(MamRawData)[f] <- FolderNames[f] 
}

rm(FolderNames)

MamminsixSub <- list()
for (s in 1:length(MamRawData)) {
  df <- MamRawData[[s]]
  ColNames <- colnames(df)
  for (t in 2:length(ColNames)) {
    dfOrdered <- df[order(df[ColNames[[t]]])[1:6],]
    NormValue <- mean(dfOrdered[[ColNames[[t]]]])
    rm(dfOrdered)
    df[ColNames[[t]]] <- df[ColNames[[t]]] - NormValue
    rm(NormValue)
  }
  MamminsixSub[[s]] <- df
  names(MamminsixSub)[s] <- names(MamRawData[s])
  rm(df,ColNames)
}

MamBgNorm <- list()
for (s in 1:length(MamRawData)) {
  df <- MamRawData[[s]]
  ColNames <- colnames(df)
  for (t in 2:length(ColNames)) {
    NormValue <- MamBackgroundValues$MamBG[which(grepl(ColNames[[t]], MamBackgroundValues$Nucleus))]
    df[ColNames[[t]]] <- df[ColNames[[t]]] / NormValue
    rm(NormValue)
  }
  MamBgNorm[[s]] <- df
  names(MamBgNorm)[s] <- names(MamRawData[s])
  rm(df,ColNames)
}

MeanMamRawData <- MeanDF(MamRawData)
MeanMamBgNorm <- MeanDF(MamBgNorm)
MeanMamMinSixSub <- MeanDF(MamminsixSub)
#### end ####

#### Calculate fold increase of band over background ####

RelIntMam <- NULL
for (s in 1:length(MamRawData)) {
  df <- MamRawData[[s]]
  df$`Distance (um)` <- NULL
  Band <- colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2))))
  dfOrdered <- NULL
  for (c in 1:length(colnames(df))) {
    Column <- as.data.frame(df[,c])
    MinSix <- Column[order(Column$`df[, c]`),][1:6]
    dfOrdered <- cbind(dfOrdered,MinSix)
  }
  Background <- colMeans(dfOrdered)
  RelInt <- as.data.frame(Band/Background)
  RelInt$Genotype <- names(MamRawData[s])
  RelInt$Nucleus <- rownames(RelInt)
  RelInt <- rename(RelInt, RelativeBandIntensity = `Band/Background`)
  RelIntMam <- rbind(RelIntMam,RelInt)
  rm(RelInt, Band, Background, df)
}
## Here we will use the mean of the 5 middle values and divide by the mean of the 5 values at
## each edge of the band (10 values in total).
RelIntMamEdge <- NULL
for (s in 1:length(MamRawData)) {
  df <- MamRawData[[s]]
  df$`Distance (um)` <- NULL
  Band <- colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2))))
  Background <- colMeans(rbind(slice_head(df, n=5),slice_tail(df,n=5)))
  RelInt <- as.data.frame(Band/Background)
  RelInt$Genotype <- names(MamRawData[s])
  RelInt$Nucleus <- rownames(RelInt)
  RelInt <- rename(RelInt, RelativeBandIntensity = `Band/Background`)
  RelIntMamEdge <- rbind(RelIntMamEdge,RelInt)
  rm(RelInt, Band, Background, df)
}

RelIntMamBg <- NULL
for (s in 1:length(MamRawData)) {
  df <- MamRawData[[s]]
  df$`Distance (um)` <- NULL
  ColNames <- colnames(df)
  BandValues <- as.data.frame(colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2)))))
  for (t in 1:length(ColNames)) {
    Band <- BandValues[t,]
    Background <- MamBackgroundValues$MamBG[which(grepl(ColNames[[t]], MamBackgroundValues$Nucleus))]
    RelInt <- as.data.frame(Band/Background)
    RelInt$Genotype <- names(MamRawData[s])
    RelInt$Nucleus <- ColNames[[t]]
    RelInt <- rename(RelInt, RelativeBandIntensity = `Band/Background`)
    rm(Band, Background, BandSubBg)
    RelIntMamBg <- rbind(RelIntMamBg,RelInt)
  }
  rm(RelInt, BandValues, df, ColNames)
}

RelIntMamBg0 <- NULL
for (s in 1:length(MamRawData)) {
  df <- MamRawData[[s]]
  df$`Distance (um)` <- NULL
  ColNames <- colnames(df)
  BandValues <- as.data.frame(colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2)))))
  for (t in 1:length(ColNames)) {
    Band <- BandValues[t,]
    Background <- MamBackgroundValues$MamBG[which(grepl(ColNames[[t]], MamBackgroundValues$Nucleus))]
    BandSubBg <- Band - Background
    RelInt <- as.data.frame(BandSubBg/Band)
    RelInt$Genotype <- names(MamRawData[s])
    RelInt$Nucleus <- ColNames[[t]]
    RelInt <- rename(RelInt, RelativeBandIntensity = `BandSubBg/Band`)
    rm(Band, Background, BandSubBg)
    RelIntMamBg0 <- rbind(RelIntMamBg0,RelInt)
  }
  rm(RelInt, BandValues, df, ColNames)
}
RelIntMamMinSixSub <- NULL
for (s in 1:length(MamRawData)) {
  df <- MamRawData[[s]]
  df$`Distance (um)` <- NULL
  Band <- colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2))))
  dfOrdered <- NULL
  for (c in 1:length(colnames(df))) {
    Column <- as.data.frame(df[,c])
    MinSix <- Column[order(Column$`df[, c]`),][1:6]
    dfOrdered <- cbind(dfOrdered,MinSix)
  }
  Background <- colMeans(dfOrdered)
  RelInt <- as.data.frame(Band-Background)
  RelInt$Genotype <- names(MamRawData[s])
  RelInt$Nucleus <- rownames(RelInt)
  RelInt <- rename(RelInt, RelativeBandIntensity = `Band - Background`)
  RelIntMamMinSixSub <- rbind(RelIntMamMinSixSub,RelInt)
  rm(RelInt, Band, Background, df,dfOrdered,c,Column,MinSix)
}
RelIntChangeMam <- NULL
for (s in 1:length(MamRawData)) {
  df <- MamRawData[[s]]
  df$`Distance (um)` <- NULL
  Band <- colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2))))
  dfOrdered <- NULL
  for (c in 1:length(colnames(df))) {
    Column <- as.data.frame(df[,c])
    MinSix <- Column[order(Column$`df[, c]`),][1:6]
    dfOrdered <- cbind(dfOrdered,MinSix)
  }
  Background <- colMeans(dfOrdered)
  RelInt <- as.data.frame((Band-Background)/Background)
  RelInt$Genotype <- names(MamRawData[s])
  RelInt$Nucleus <- rownames(RelInt)
  RelInt <- rename(RelInt, RelativeBandIntensity = `(Band - Background)/Background`)
  RelIntChangeMam <- rbind(RelIntChangeMam,RelInt)
  rm(RelInt, Band, Background, df, c, MinSix, dfOrdered, Column)
}
RelIntChangeMamBg <- NULL
for (s in 1:length(MamRawData)) {
  df <- MamRawData[[s]]
  df$`Distance (um)` <- NULL
  ColNames <- colnames(df)
  BandValues <- as.data.frame(colMeans(slice(df,((ceiling(nrow(df)/2)-2):(ceiling(nrow(df)/2)+2)))))
  for (t in 1:length(ColNames)) {
    Band <- BandValues[t,]
    Background <- MamBackgroundValues$MamBG[which(grepl(ColNames[[t]], MamBackgroundValues$Nucleus))]
    RelInt <- as.data.frame((Band-Background)/Background)
    RelInt$Genotype <- names(MamRawData[s])
    RelInt$Nucleus <- ColNames[[t]]
    RelInt <- rename(RelInt, RelativeBandIntensity = `(Band - Background)/Background`)
    rm(Band, Background, BandSubBg)
    RelIntChangeMamBg <- rbind(RelIntChangeMamBg,RelInt)
  }
  rm(RelInt, BandValues, df, ColNames)
}


#### end ####

#### Saved plots ####
setwd("/Users/jonathantownson/Documents/PhD/CSL localisation Paper/Figures for paper/Files/Fig 4")
## Create and save band plots of the bands compared to their own wt data 

df <- MeanMamRawData %>% dplyr::filter(Name == "Dark Control" | Name == "2 hours")
XAxisOrder <- c("Dark Control (N = 24)", "2 hours (N = 33)")
df <- df %>% mutate(Genotype = fct_relevel(Genotype, XAxisOrder))

Labs <- c(title = expression(paste("eGFP::CSL at ",italic("E(spl)-C"))), x = "Distance (μm)", y = "eGFP::CSL \nFluorescence (AU)")
band_plot(df, Labs) + theme(legend.position = "bottom", legend.title = element_blank()) + 
  guides(colour = guide_legend(nrow=2, byrow = TRUE)) + 
  scale_y_continuous(breaks = seq(0, 800, by = 200), limits = c(0,800))
ggsave("RawData Band plot of 2 hours and dark control CSL.jpg", device = "jpeg", dpi = "retina",
         width = 15, height = 15, units = "cm")

df <- MeanMamMinSixSub %>% dplyr::filter(Name == "Dark Control" | Name == "2 hours")
XAxisOrder <- c("Dark Control (N = 24)", "2 hours (N = 33)")
df <- df %>% mutate(Genotype = fct_relevel(Genotype, XAxisOrder))

Labs <- c(title = " ", x = "Distance (μm)", y = "EGFP::CSL Fluorescence (AU)")
band_plot(df, Labs) + theme(legend.position = "bottom", legend.title = element_blank()) + 
  guides(colour = guide_legend(nrow=1, byrow = TRUE)) +
  scale_y_continuous(breaks = seq(0, 800, by = 200), limits = c(0,600))
ggsave("MinSix Band plot of 2 hours and dark control CSL.jpg", device = "jpeg", dpi = "retina",
       width = 12, height = 12, units = "cm")

#### end ####
