# clear the environment: 
rm(list = ls())
# clear the console: ctrl+L
# clear all the plots: 
dev.off()

#### Packages for this script ####
install.packages("ggplot2")
install.packages("dplyr")
install.packages("stringr")
install.packages("forcats")
library(ggplot2)
library(dplyr)
library(stringr)
library(forcats)
library(ggsignif)
#### end ####

# Set location of data and where plots will be saved to:
setwd("/Users/jonathantownson/Documents/PhD/Images/Analysed image data/Optogenetics data/newBLITz in light incubator/Nuclear levels csv files")
# Create the path where the data is stored
DataPath <- "/Users/jonathantownson/Documents/PhD/Images/Analysed image data/Optogenetics data/newBLITz in light incubator/Nuclear levels csv files"

#### Read in and normalise the data ####
# Create the folder names for the data
FolderNames <- list.dirs(DataPath, full.names = FALSE, recursive = FALSE)
# Read in all the csv files for each condition
SuHRawData <- list()
for (f in 1:length(FolderNames)){
  Path <- paste(DataPath,FolderNames[f], sep ="/")
  Files <- list.files(Path, pattern = "*.csv")
  FilePaths <- list()
  for (i in 1:length(Files)){
    if (grepl("SuH",Files[i])) {
      FilePaths[[i]] <- paste(Path,Files[i], sep ="/")
    }
  } 
  FilePaths[sapply(FilePaths, is.null)] <- NULL
  rm(Path)
  rm(Files)
  Data <- list()
  for (i in 1:length(FilePaths)){
    df <- read.csv(FilePaths[[i]])
    dfOdd <- df[1:nrow(df)%%2!=0,]
    for (c in 2:length(colnames(dfOdd))) {
      names(dfOdd)[names(dfOdd) == colnames(dfOdd)[c]] <- paste("Cell", str_sub(colnames(dfOdd)[c]))
    }
    rm(c)
    dfOdd$Cell <- ceiling(dfOdd$X/2)
    dfOdd$X <- NULL
    dfEven <- df[1:nrow(df)%%2==0,]
    for (c in 2:length(colnames(dfEven))) {
      names(dfEven)[names(dfEven) == colnames(dfEven)[c]] <- paste("Nucleus", str_sub(colnames(dfEven)[c]))
    }
    dfEven$Cell <- ceiling(dfEven$X/2)
    dfEven$X <- NULL
    df <- merge(dfOdd, dfEven, by = "Cell")
    rm(dfOdd,dfEven,c)
    df$'Cytoplasm Area' <- df$`Cell Area` - df$`Nucleus Area`
    df$'Cytoplasm IntDen' <- df$`Cell IntDen` - df$`Nucleus IntDen`
    df$'Cytoplasm Mean' <- df$`Cytoplasm IntDen` / df$`Cytoplasm Area`
    df$NCratio <- df$`Nucleus Mean` / df$`Cytoplasm Mean`
    df$Gland <- str_sub(basename(FilePaths[[i]]), end=-5)
    df$IncubationTime <- FolderNames[f]
    Data[[i]] <- df
    rm(df)
  }
  Data <- do.call(rbind.data.frame, Data)
  SuHRawData[[f]] <- Data
  rm(Data,i)
  rm(FilePaths)
  names(SuHRawData)[f] <- FolderNames[f] 
}
rm(f)

NICDRawData <- list()
for (f in 1:length(FolderNames)){
  Path <- paste(DataPath,FolderNames[f], sep ="/")
  Files <- list.files(Path, pattern = "*.csv")
  FilePaths <- list()
  for (i in 1:length(Files)){
    if (grepl("NICD",Files[i])) {
      FilePaths[[i]] <- paste(Path,Files[i], sep ="/")
    }
  } 
  FilePaths[sapply(FilePaths, is.null)] <- NULL
  rm(Path)
  rm(Files)
  Data <- list()
  for (i in 1:length(FilePaths)){
    df <- read.csv(FilePaths[[i]])
    dfOdd <- df[1:nrow(df)%%2!=0,]
    for (c in 2:length(colnames(dfOdd))) {
      names(dfOdd)[names(dfOdd) == colnames(dfOdd)[c]] <- paste("Cell", str_sub(colnames(dfOdd)[c]))
    }
    rm(c)
    dfOdd$Cell <- ceiling(dfOdd$X/2)
    dfOdd$X <- NULL
    dfEven <- df[1:nrow(df)%%2==0,]
    for (c in 2:length(colnames(dfEven))) {
      names(dfEven)[names(dfEven) == colnames(dfEven)[c]] <- paste("Nucleus", str_sub(colnames(dfEven)[c]))
    }
    dfEven$Cell <- ceiling(dfEven$X/2)
    dfEven$X <- NULL
    df <- merge(dfOdd, dfEven, by = "Cell")
    rm(dfOdd,dfEven,c)
    df$'Cytoplasm Area' <- df$`Cell Area` - df$`Nucleus Area`
    df$'Cytoplasm IntDen' <- df$`Cell IntDen` - df$`Nucleus IntDen`
    df$'Cytoplasm Mean' <- df$`Cytoplasm IntDen` / df$`Cytoplasm Area`
    df$NCratio <- df$`Nucleus Mean` / df$`Cytoplasm Mean`
    df$Gland <- str_sub(basename(FilePaths[[i]]), end=-5)
    df$IncubationTime <- FolderNames[f]
    Data[[i]] <- df
    rm(df)
  }
  Data <- do.call(rbind.data.frame, Data)
  NICDRawData[[f]] <- Data
  rm(Data, i)
  rm(FilePaths)
  names(NICDRawData)[f] <- FolderNames[f] 
}

XAxisOrder <- str_sort(FolderNames, numeric = TRUE)

rm(FolderNames,f)


NICDdata <- do.call(rbind.data.frame, NICDRawData)
NICDdata <- dplyr::select(NICDdata,c("Cell","Nucleus Area","Nucleus Mean","Cytoplasm Area","Cytoplasm Mean","NCratio","Gland","IncubationTime"))
NICDdata <- NICDdata %>% mutate(IncubationTime = fct_relevel(IncubationTime, XAxisOrder)) 

SuHdata <- do.call(rbind.data.frame, SuHRawData)
SuHdata <- dplyr::select(SuHdata,c("Cell","Nucleus Area","Nucleus Mean","Cytoplasm Area","Cytoplasm Mean","NCratio","Gland","IncubationTime"))
SuHdata <- SuHdata %>% mutate(IncubationTime = fct_relevel(IncubationTime, XAxisOrder)) 

#### end ####

#### Statistics functions ####

## The first thing we look for is normality of the data. Note that anything with
## a sample size greater than 30 can be considered normal because of the central
## limit theorem.

## We can do density plots of the band intensities to make sure the data is
## normally distributed

# This function will create a JPEG file of data density on y axis, it will need
# a list of data, the name of the variable density we are looking at, something
# to separate the plots by (e.g. Genotype), and the name of the JPEG file
NormDis_plots <- function(DataList, VariableDens, SepBy, MainTitle){
  Plots <- list()
  for (p in 1:length(DataList)) {
    df <- DataList[[p]]
    Plot <- ggdensity(df, x = VariableDens, facet.by = SepBy)
    Plots[[p]] <- Plot
    rm(Plot)
  }
  ggsave(paste(MainTitle,".jpg"), arrangeGrob(grobs = Plots, ncol = 2),
         device = "jpeg", dpi = "retina", width = 24,
         height = 9*ceiling(length(Plots)/2), units = "cm")
  rm(Plots)
}

## A quantile-quantile plot can also be used to check for normal distribution
## the data should align with the 45 degree line if it is normal.

# This function will create a JPEG file of QQ plots with sample on y axis, it will need
# a list of data, the name of the variable we are looking at, something
# to separate the plots by (e.g. Genotype), and the name of the JPEG file
QQ_plots <- function(DataList, VariableDens, SepBy, MainTitle){
  Plots <- list()
  for (p in 1:length(DataList)) {
    df <- DataList[[p]]
    Plot <- ggqqplot(df, x = VariableDens, facet.by = SepBy)
    Plots[[p]] <- Plot
    rm(Plot)
  }
  ggsave(paste(MainTitle,".jpg"), arrangeGrob(grobs = Plots, ncol = 2),
         device = "jpeg", dpi = "retina", width = 24,
         height = 9*ceiling(length(Plots)/2), units = "cm")
  rm(Plots)
}

## Perform the Shapiro-Wilk test of normality on samples. In this test the null
## hypothesis is that the data is normal and therefore if the p value <0.05 the 
## data should be considered to NOT be normal. The test works best on small sample
## sizes and should be considered alongside the graphs.

# This function will create a csv file with some summary statistics for a list
# including the Shapiro value. It also creates a list of these summaries for
# viewing in R.

VariableList <- c("Nucleus Area","Nucleus Mean","Cytoplasm Area","Cytoplasm Mean","NCratio")
DataStatSumTestListFun <- function(DataList,Filename){
  OutputList <- list()
  for (p in 1:length(DataList)) {
    df <- na.omit(DataList[[p]])
    dfStatistics <- data.frame()
    for (g in 1:length(unique(df$IncubationTime))) {
      for (v in 1:length(VariableList)) {
        Testdf <- dplyr::filter(df,IncubationTime == unique(df$IncubationTime)[g])
        MeanValue <- mean(Testdf[[VariableList[v]]])
        MedianValue <- median(Testdf[[VariableList[v]]])
        VarValue <- var(Testdf[[VariableList[v]]])
        StandardDeviation <- sd(Testdf[[VariableList[v]]])
        StandardError <- StandardDeviation/sqrt(length(Testdf[[VariableList[v]]]))
        result <- shapiro.test(Testdf[[VariableList[v]]])
        SWpvalue <- result[["p.value"]]
        StatsResults <- c(unique(df$IncubationTime)[g],VariableList[v],MeanValue,MedianValue,VarValue,StandardDeviation,StandardError,SWpvalue)
        dfStatistics <- rbind(dfStatistics,StatsResults)
        rm(Testdf,result,MeanValue,MedianValue,VarValue,StandardDeviation,StandardError,SWpvalue,StatsResults)
      }
    }
    names(dfStatistics) <- c("IncubationTime","Variable","Mean","Median","Variance","Standard Deviation","Standard Error","Shapiro-Wilk p-value")
    OutputList[[p]] <- dfStatistics
    names(OutputList)[p] <- names(DataList[p])
    rm(df,dfStatistics)
  }
  SaveResult <- as.data.frame(OutputList)
  SaveResult <- t(SaveResult)
  write.csv(SaveResult, file = paste(Filename,".csv"))
  return(OutputList)
}

## Finally we want to generate a set of p-values for the data which can be checked for
## significance and compared to the summaries generated above to know which is most
## relevant. For normal data with equal variance it is the t test, unequal variance use
## the Welch's t test. If the data is nonparametric (a.k.a. not normal), then it is
## the Mann-Whitney test.

Stats_tests <- function(DataList, Filename){
  dfStatistics <- data.frame()
  for (p in 1:(length(DataList))) {
    df1 <- na.omit(DataList[[p]])
    for (q in 1:(length(DataList))) { 
      df2 <- na.omit(DataList[[q]])
      Comparison <- paste(unique(df1$IncubationTime), "vs", unique(df2$IncubationTime))
      StudentT <- t.test(df1$NCratio, df2$NCratio, var.equal = TRUE)
      StudentT <- StudentT[["p.value"]]
      WelchT <- t.test(df1$NCratio, df2$NCratio, var.equal = FALSE)
      WelchT <- WelchT[["p.value"]]
      MannWhitney <- wilcox.test(df1$NCratio,df2$NCratio)
      MannWhitney <- MannWhitney[["p.value"]]
      SummaryVec <- c(Comparison, StudentT, WelchT, MannWhitney)
      dfStatistics <- rbind(dfStatistics, SummaryVec)
      rm(SummaryVec,df2,Comparison, StudentT, WelchT, MannWhitney)
    }
    rm(df1,q)
  }
  names(dfStatistics) <- c("Comparison","StudentT","WelchT","MannWhitney")
  write.csv(dfStatistics, file = paste(Filename,".csv"))
  return(dfStatistics)
}
NuclearMean_Stats_tests <- function(DataList, Filename){
  dfStatistics <- data.frame()
  for (p in 1:(length(DataList))) {
    df1 <- na.omit(DataList[[p]])
    for (q in 1:(length(DataList))) { 
      df2 <- na.omit(DataList[[q]])
      Comparison <- paste(unique(df1$IncubationTime), "vs", unique(df2$IncubationTime))
      StudentT <- t.test(df1$`Nucleus Mean`, df2$`Nucleus Mean`, var.equal = TRUE)
      StudentT <- StudentT[["p.value"]]
      WelchT <- t.test(df1$`Nucleus Mean`, df2$`Nucleus Mean`, var.equal = FALSE)
      WelchT <- WelchT[["p.value"]]
      MannWhitney <- wilcox.test(df1$`Nucleus Mean`,df2$`Nucleus Mean`)
      MannWhitney <- MannWhitney[["p.value"]]
      SummaryVec <- c(Comparison, StudentT, WelchT, MannWhitney)
      dfStatistics <- rbind(dfStatistics, SummaryVec)
      rm(SummaryVec,df2,Comparison, StudentT, WelchT, MannWhitney)
    }
    rm(df1,q)
  }
  names(dfStatistics) <- c("Comparison","StudentT","WelchT","MannWhitney")
  write.csv(dfStatistics, file = paste("NuclearMean",Filename,".csv"))
  return(dfStatistics)
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
  axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),  # tilt the x axis
  axis.line = element_line(colour = "black", size = 1.5), # add axis lines in black
  axis.ticks = element_line(colour = "black", size = 1.5),
  axis.ticks.length=unit(2,"mm")) # add tick marks in black
}

NCratio_crossbar <- function(Data, Labels) {
  ggplot(Data,aes(IncubationTime,NCratio,col=IncubationTime,shape=IncubationTime)) + theme_jmt() +
    geom_jitter(width=0.25) +
    stat_summary(fun="mean", geom="crossbar", width=0.5) + theme(legend.key=element_blank())  + 
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]]) + 
    scale_color_brewer(palette = "Set2")
}

NucleusMean_crossbar <- function(Data, Labels) {
  ggplot(Data,aes(IncubationTime,`Nucleus Mean`,col=IncubationTime,shape=IncubationTime)) + theme_jmt() +
    geom_jitter(width=0.25) +
    stat_summary(fun="mean", geom="crossbar", width=0.5) + theme(legend.key=element_blank())  + 
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]]) + 
    scale_color_brewer(palette = "Set2")
}

CytoplasmMean_crossbar <- function(Data, Labels) {
  ggplot(Data,aes(IncubationTime,`Cytoplasm Mean`,col=IncubationTime,shape=IncubationTime)) + theme_jmt() +
    geom_jitter(width=0.25) +
    stat_summary(fun="mean", geom="crossbar", width=0.5) + theme(legend.key=element_blank())  + 
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]]) + 
    scale_color_brewer(palette = "Set2")
}

#### end ####

DarkControl <- NICDdata %>% dplyr::filter(IncubationTime == "Dark control")
Light2 <- NICDdata %>% dplyr::filter(IncubationTime == "2 hours light")
Light5 <- NICDdata %>% dplyr::filter(IncubationTime == "5 hours light")
Light8 <- NICDdata %>% dplyr::filter(IncubationTime == "8 hours light")
Light12 <- NICDdata %>% dplyr::filter(IncubationTime == "12 hours light")
Light24 <- NICDdata %>% dplyr::filter(IncubationTime == "24 hours light")

df1 <- rbind(DarkControl,Light2)
df1list <- list()
for (g in 1:length(unique(df1$IncubationTime))) {
  df1list[[g]] <- df1[str_detect(df1$IncubationTime, as.character(unique(df1$IncubationTime)[g])),]
  names(df1list)[g] <- as.character(unique(df1$IncubationTime)[g])
}
df1 <- df1 %>% mutate(IncubationTime = fct_relevel(IncubationTime,"Dark control", "2 hours light"))
levels(df1$IncubationTime) <- gsub("2 hours light", "2hr Light", levels(df1$IncubationTime))

setwd("/Users/jonathantownson/Documents/PhD/CSL localisation Paper/Figures for paper/Files/Fig 4")

Labs <- c(title = "OptIC-Notch{ω}", x = "Incubation Time", y = "Nuclear:Cytoplasmic ratio")
NCratio_crossbar(df1, Labs) + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5), axis.title.x=element_blank(), legend.position = "none") +
  geom_signif(comparisons = list(c("Dark control","2hr Light")),y_position = 5.5, 
              map_signif_level=TRUE, test = "t.test", color = "black") + 
  scale_y_continuous(breaks = seq(0, 6, by = 2), limits = c(0,6))
ggsave("NCratio of NICDmCherry.jpg", device = "jpeg", dpi = "retina",
width = 15, height = 15, units = "cm")

NormDis_plots(df1list, "NCratio", "IncubationTime","Normal distribution of NICD NCratio plots for newblitz with cry in dark or 2hr light")
QQ_plots(df1list, "NCratio", "IncubationTime",paste("QQ plots of NICD NCratio plots for newblitz with cry in dark or 2hr light"))
DataStatSumTestListFun(df1list,paste("Summary of NICD NCratio newblitz with cry in dark or 2hr light"))
Stats_tests(df1list,paste("SigStats of NICD NCratio for newblitz with cry in dark or 2hr light"))
rm(df1, df1list)