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
install.packages("ggsignif")
install.packages("ggpubr")
install.packages("gridExtra")
library(ggplot2)
library(dplyr)
library(stringr)
library(forcats)
library(ggsignif)
library(ggpubr)
library(gridExtra)
#### end ####

# Co factor to be analysed
CSLFactor <- "SuH"
CSLLabel <- "Halo::CSL"
NICDFactor <- "NICD"
NICDLabel <- "OptIC-Notch"
# Set location of data and where plots will be saved to:
mainDir <- "/Users/jonathantownson/Documents/PhD/Images/Analysed image data/SuH localisation/Halo gland image quantifications"
setwd(mainDir)

#### Read in and normalise the CSL data ####
# Create the path where the data is stored
DataPath <- paste(mainDir, CSLFactor, sep = "/")
# Create the folder names for the data
FolderNames <- list.dirs(DataPath, full.names = FALSE, recursive = FALSE)
# Read in all the csv files for each condition
SuHRawData <- list()
for (f in 1:length(FolderNames)){
  Path <- paste(DataPath,FolderNames[f], sep ="/")
  Files <- list.files(Path, pattern = "*.csv")
  FilePaths <- list()
  for (i in 1:length(Files)){
    FilePaths[[i]] <- paste(Path,Files[i], sep ="/")
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
    df$`Gland Mean` <- mean(df$`Cell Mean`)
    df$`Gland IntDen` <- sum(df$`Cell RawIntDen`)
    df$`Cytoplasm Area` <- df$`Cell Area` - df$`Nucleus Area`
    df$`Cytoplasm IntDen` <- df$`Cell IntDen` - df$`Nucleus IntDen`
    df$`Cytoplasm Mean` <- df$`Cytoplasm IntDen` / df$`Cytoplasm Area`
    df$`GlandCyt Mean` <- mean(df$`Cytoplasm Mean`)
    df$`GlandCyt IntDen` <- sum(df$`Cytoplasm IntDen`)
    df$NCratio <- df$`Nucleus Mean` / df$`Cytoplasm Mean`
    df$NGratio <- df$`Nucleus Mean` / df$`Gland Mean`
    df$NGCytratio <- df$`Nucleus Mean` / df$`GlandCyt Mean`
    df$NuclearRelativeChange <- (df$`Nucleus Mean`-df$`Cytoplasm Mean`) / df$`Cytoplasm Mean`
    df$NuclearPercent <- (df$`Nucleus Mean`/(df$`Cell Mean`+df$`Nucleus Mean`))*100
    df$NuclearPercentIntDen <- (df$`Nucleus IntDen`/(df$`Cell IntDen`))*100
    df$NuclearPercentGlandIntDen <- (df$`Nucleus IntDen`/(df$`Gland IntDen`))*100
    df$NuclearPercentGlandCytIntDen <- (df$`Nucleus IntDen`/(df$`GlandCyt IntDen`))*100
    df$NucAreaPerc <- (df$`Nucleus Area`/(df$`Cell Area`))*100
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

rm(FolderNames,f)

SuHdata <- do.call(rbind.data.frame, SuHRawData)
SuHdata <- dplyr::select(SuHdata,c("Cell","Gland","IncubationTime","Nucleus Area","Nucleus Mean","Cytoplasm Area",
                                   "Cytoplasm Mean","Gland Mean","GlandCyt Mean","NCratio","NGratio","NGCytratio",
                                   "NuclearRelativeChange","NucAreaPerc","NuclearPercent","NuclearPercentIntDen",
                                   "NuclearPercentGlandIntDen","NuclearPercentGlandCytIntDen"))


#### end ####

#### Read in and normalise the NICD data ####
# Create the path where the data is stored
DataPath <- paste(mainDir, NICDFactor, sep = "/")
# Create the folder names for the data
FolderNames <- list.dirs(DataPath, full.names = FALSE, recursive = FALSE)
# Read in all the csv files for each condition
NICDRawData <- list()
for (f in 1:length(FolderNames)){
  Path <- paste(DataPath,FolderNames[f], sep ="/")
  Files <- list.files(Path, pattern = "*.csv")
  FilePaths <- list()
  for (i in 1:length(Files)){
    FilePaths[[i]] <- paste(Path,Files[i], sep ="/")
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
    df$`Gland Mean` <- mean(df$`Cell Mean`)
    df$`Gland IntDen` <- sum(df$`Cell RawIntDen`)
    df$`Cytoplasm Area` <- df$`Cell Area` - df$`Nucleus Area`
    df$`Cytoplasm IntDen` <- df$`Cell IntDen` - df$`Nucleus IntDen`
    df$`Cytoplasm Mean` <- df$`Cytoplasm IntDen` / df$`Cytoplasm Area`
    df$`GlandCyt Mean` <- mean(df$`Cytoplasm Mean`)
    df$`GlandCyt IntDen` <- sum(df$`Cytoplasm IntDen`)
    df$NCratio <- df$`Nucleus Mean` / df$`Cytoplasm Mean`
    df$NGratio <- df$`Nucleus Mean` / df$`Gland Mean`
    df$NGCytratio <- df$`Nucleus Mean` / df$`GlandCyt Mean`
    df$NuclearRelativeChange <- (df$`Nucleus Mean`-df$`Cytoplasm Mean`) / df$`Cytoplasm Mean`
    df$NuclearPercent <- (df$`Nucleus Mean`/(df$`Cell Mean`+df$`Nucleus Mean`))*100
    df$NuclearPercentIntDen <- (df$`Nucleus IntDen`/(df$`Cell IntDen`))*100
    df$NuclearPercentGlandIntDen <- (df$`Nucleus IntDen`/(df$`Gland IntDen`))*100
    df$NuclearPercentGlandCytIntDen <- (df$`Nucleus IntDen`/(df$`GlandCyt IntDen`))*100
    df$NucAreaPerc <- (df$`Nucleus Area`/(df$`Cell Area`))*100
    df$Gland <- str_sub(basename(FilePaths[[i]]), end=-5)
    df$IncubationTime <- FolderNames[f]
    Data[[i]] <- df
    rm(df)
  }
  Data <- do.call(rbind.data.frame, Data)
  NICDRawData[[f]] <- Data
  rm(Data,i)
  rm(FilePaths)
  names(NICDRawData)[f] <- FolderNames[f] 
}

rm(FolderNames,f)

NICDdata <- do.call(rbind.data.frame, NICDRawData)
NICDdata <- dplyr::select(NICDdata,c("Cell","Gland","IncubationTime","Nucleus Area","Nucleus Mean","Cytoplasm Area",
                                   "Cytoplasm Mean","Gland Mean","GlandCyt Mean","NCratio","NGratio","NGCytratio",
                                   "NuclearRelativeChange","NucAreaPerc","NuclearPercent","NuclearPercentIntDen",
                                   "NuclearPercentGlandIntDen","NuclearPercentGlandCytIntDen"))


#### end ####

#### Calculate relative changes of NCratio in data ####
InitialNCratioValues <- SuHRawData[["010 min"]]$NCratio
InitialNGratioValues <- SuHRawData[["010 min"]]$NGratio

SuHRelChangeData <- list()
for (f in 1:length(SuHRawData)) {
  df <- SuHRawData[[names(SuHRawData[f])]]
  df <- dplyr::select(df,c("Cell","Gland","IncubationTime","Nucleus Mean","Cytoplasm Mean","NCratio","NGratio"))
  df$RelChangeNCratio <- (df$NCratio-InitialNCratioValues)/InitialNCratioValues
  df$RelChangeNGratio <- (df$NGratio-InitialNGratioValues)/InitialNGratioValues
  SuHRelChangeData[[f]] <- df
  rm(df)
  names(SuHRelChangeData)[f] <- names(SuHRawData)[f]
}
SuHRelChange <- do.call(rbind.data.frame, SuHRelChangeData)

InitialNCratioValues <- NICDRawData[["010 min"]]$NCratio
InitialNGratioValues <- NICDRawData[["010 min"]]$NGratio

NICDRelChangeData <- list()
for (f in 2:length(NICDRawData)) {
  df <- NICDRawData[[names(NICDRawData[f])]]
  df <- dplyr::select(df,c("Cell","Gland","IncubationTime","Nucleus Mean","Cytoplasm Mean","NCratio","NGratio"))
  df$RelChangeNCratio <- (df$NCratio-InitialNCratioValues)/InitialNCratioValues
  df$RelChangeNGratio <- (df$NGratio-InitialNGratioValues)/InitialNGratioValues
  NICDRelChangeData[[f]] <- df
  rm(df)
  names(NICDRelChangeData)[f] <- names(NICDRawData)[f]
}
NICDRelChange <- do.call(rbind.data.frame, NICDRelChangeData)

rm(InitialNCratioValues,InitialNGratioValues)

df <- subset(SuHRelChange, IncubationTime!= "000 min")
df1 <- subset(NICDRelChange, IncubationTime!= "000 min")

RelativeChanges <- dplyr::select(df,c("Cell","Gland","IncubationTime","Nucleus Mean","Cytoplasm Mean","NCratio","NGratio"))
RelativeChanges$SuHRelChangeNCratio <- df$RelChangeNCratio
RelativeChanges$SuHRelChangeNGratio <- df$RelChangeNGratio
RelativeChanges$NICDRelChangeNCratio <- df1$RelChangeNCratio
RelativeChanges$NICDRelChangeNGratio <- df1$RelChangeNGratio

rm(df, df1)

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
  ggplot(Data,aes(IncubationTime,NCratio)) + theme_jmt() +
    geom_jitter(width=2.5, col = "darkgray") +
    stat_summary(fun="mean", geom="crossbar", width=5) + theme(legend.key=element_blank())  + 
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]]) + 
    scale_color_brewer(palette = "Set2") + 
    scale_x_continuous(breaks = seq(0, 120, by = 20), limits = c(-5,125))
}
NCratio_crossbar_GlandCol <- function(Data, Labels) {
  ggplot(Data,aes(IncubationTime,NCratio,col=Gland,shape=Gland)) + theme_jmt() +
    geom_jitter(width=2.5) +
    stat_summary(fun="mean", geom="crossbar", width=5) + theme(legend.key=element_blank())  + 
    labs(title = Labels[[1]], x = Labels[[2]], y = Labels[[3]]) + 
    scale_color_brewer(palette = "Set2") + 
    scale_x_continuous(breaks = seq(0, 120, by = 20), limits = c(-5,125))
}
#### end ####

#### Saved plots ####
setwd("/Users/jonathantownson/Documents/PhD/CSL localisation Paper/Figures for paper/Files/Fig 2")
##CSL
df <- SuHdata
df$IncubationTime <- as.numeric(gsub(".*?([0-9]+).*", "\\1", df$IncubationTime))
Labs <- c(title = CSLLabel, x = "Time (mins)", y = "Nuclear:Cytoplasmic ratio")
NCratio_crossbar(df, Labs) + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
ggsave(paste("OldBLITzGFP with CRY activation over 2hr NCratio of",CSLFactor,".jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")
NCratio_crossbar_GlandCol(df, Labs) + theme(axis.text.x = element_text(angle = 0, vjust = 2, hjust=0.5))
ggsave(paste("OldBLITzGFP with CRY activation over 2hr glands coloured NCratio of",CSLFactor,".jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

##NICD
df <- NICDdata[-(1:5), , drop = FALSE]
df$IncubationTime <- as.numeric(gsub(".*?([0-9]+).*", "\\1", df$IncubationTime))
Labs <- c(title = NICDLabel, x = "Time (mins)", y = "Nuclear:Cytoplasmic ratio")
NCratio_crossbar(df, Labs) + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
ggsave(paste("OldBLITzGFP with CRY activation over 2hr NCratio of",NICDFactor,".jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")
NCratio_crossbar_GlandCol(df, Labs) + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
ggsave(paste("OldBLITzGFP with CRY activation over 2hr glands coloured NCratio of",NICDFactor,".jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

rm(df)
## Relative change correlation

sp <- ggscatter(RelativeChanges, x = "NICDRelChangeNCratio", y = "SuHRelChangeNCratio",
      add = "reg.line",  # Add regressin line
      add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
      conf.int = TRUE) + # Add confidence interval
  theme_jmt() + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5), legend.position = "none") + 
  labs(title = "CSL vs OptIC-Notch", x = paste("Relative change in \n",NICDLabel, "N:C ratio"), y = paste("Relative change in \n",CSLLabel, "N:C ratio"))
ggsave(paste("Correlation in relative change of NCratio between",NICDFactor," and",CSLFactor,".jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

# Add correlation coefficient
sp + stat_cor(method = "spearman", label.x = 0.75, label.y = 1.5)
ggsave(paste("Correlation in relative change of NCratio between",NICDFactor," and",CSLFactor,"with Spearman test.jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

cor.test(RelativeChanges$SuHRelChangeNCratio, RelativeChanges$NICDRelChangeNCratio, method = "pearson")
cor.test(RelativeChanges$SuHRelChangeNCratio, RelativeChanges$NICDRelChangeNCratio, method = "spearman")

setwd(mainDir)
#### end ####

#### Statistics Tests and Plots ####
subDir <- paste("Stats")
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))
rm(subDir)

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

NormDis_plots(SuHRawData, "NCratio", "IncubationTime",paste("Normal distribution of NCratio plots for",CSLFactor))
NormDis_plots(NICDRawData, "NCratio", "IncubationTime",paste("Normal distribution of NCratio plots for",NICDFactor))

NormDis_plots(SuHRelChangeData[-1], "RelChangeNCratio", "IncubationTime",paste("Normal distribution of RelChangeNCratio plots for",CSLFactor))
NormDis_plots(NICDRelChangeData[-1], "RelChangeNCratio", "IncubationTime",paste("Normal distribution of RelChangeNCratio plots for",NICDFactor))

ggdensity(RelativeChanges, x = "SuHRelChangeNCratio")
ggsave(paste("Normal distribution of combine RelChangeNCratio plots for",CSLFactor,".jpg"), device = "jpeg", dpi = "retina",
width = 15, height = 15, units = "cm")
ggdensity(RelativeChanges, x = "NICDRelChangeNCratio")
ggsave(paste("Normal distribution of combine RelChangeNCratio plots for",NICDFactor,".jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

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

QQ_plots(SuHRawData, "NCratio", "IncubationTime",paste("QQ Plots of",CSLFactor))
QQ_plots(NICDRawData, "NCratio", "IncubationTime",paste("QQ Plots of",NICDFactor))

QQ_plots(SuHRelChangeData[-1], "RelChangeNCratio", "IncubationTime",paste("QQ Plots of RelChangeNCratio for ",CSLFactor))
QQ_plots(NICDRelChangeData[-1], "RelChangeNCratio", "IncubationTime",paste("QQ Plots of RelChangeNCratio for ",NICDFactor))

ggqqplot(RelativeChanges, x = "SuHRelChangeNCratio")
ggsave(paste("QQ Plots of combined RelChangeNCratio plots for",CSLFactor,".jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")
ggqqplot(RelativeChanges, x = "NICDRelChangeNCratio")
ggsave(paste("QQ Plots of combined RelChangeNCratio plots for",NICDFactor,".jpg"), device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")
## Perform the Shapiro-Wilk test of normality on samples. In this test the null
## hypothesis is that the data is normal and therefore if the p value <0.05 the 
## data should be considered to NOT be normal. The test works best on small sample
## sizes and should be considered alongside the graphs.

# This function will create a csv file with some summary statistics for a list
# including the Shapiro value. It also creates a list of these summaries for
# viewing in R.

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

VariableList <- c("Nucleus Area","Nucleus Mean","Cytoplasm Area","Cytoplasm Mean","NCratio","NGratio","NuclearPercent","NuclearRelativeChange")
DataStatSumTestListFun(SuHRawData,paste("Summary of",CSLFactor, "at each interval"))
DataStatSumTestListFun(NICDRawData,paste("Summary of",NICDFactor, "at each interval"))

VariableList <- c("RelChangeNCratio","RelChangeNGratio")
DataStatSumTestListFun(SuHRelChangeData[-1],paste("Summary of Rel change for",CSLFactor, "at each interval"))
DataStatSumTestListFun(NICDRelChangeData[-1],paste("Summary of Rel change for",NICDFactor, "at each interval"))

VariableList <- c("SuHRelChangeNCratio","SuHRelChangeNGratio","NICDRelChangeNCratio","NICDRelChangeNGratio")
dfStatistics <- data.frame()
for (v in 1:length(VariableList)) {
  Testdf <- RelativeChanges
  MeanValue <- mean(Testdf[[VariableList[v]]])
  MedianValue <- median(Testdf[[VariableList[v]]])
  VarValue <- var(Testdf[[VariableList[v]]])
  StandardDeviation <- sd(Testdf[[VariableList[v]]])
  StandardError <- StandardDeviation/sqrt(length(Testdf[[VariableList[v]]]))
  result <- shapiro.test(Testdf[[VariableList[v]]])
  SWpvalue <- result[["p.value"]]
  StatsResults <- c(VariableList[v],MeanValue,MedianValue,VarValue,StandardDeviation,StandardError,SWpvalue)
  dfStatistics <- rbind(dfStatistics,StatsResults)
  names(dfStatistics) <- c("Variable","Mean","Median","Variance","Standard Deviation","Standard Error","Shapiro-Wilk p-value")
  rm(Testdf,result,MeanValue,MedianValue,VarValue,StandardDeviation,StandardError,SWpvalue,StatsResults)
}
write.csv(dfStatistics, file = "Summary of combined Relative changes.csv")
rm(dfStatistics)
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

SigStatsSuH <- Stats_tests(SuHRawData,paste("SigStats comparing NCratio of CSL"))
SigStatsNICD <- Stats_tests(NICDRawData,paste("SigStats comparing NCratio of NICD"))

## Correlation stats
dfStatistics <- data.frame()
Spearman <- cor.test(RelativeChanges$SuHRelChangeNCratio, RelativeChanges$NICDRelChangeNCratio, method = "spearman")
SpearmanRho <- Spearman[["estimate"]][["rho"]]
SpearmanP <- Spearman[["p.value"]]
Pearson <- cor.test(RelativeChanges$SuHRelChangeNCratio, RelativeChanges$NICDRelChangeNCratio, method = "pearson")
PearsonR <- Pearson[["estimate"]][["cor"]]
PearsonP <- Pearson[["p.value"]]
Pearson95CIlower <- Pearson[["conf.int"]][1]
Pearson95CIupper <- Pearson[["conf.int"]][2]
SummaryVec <- c(SpearmanRho, SpearmanP, PearsonR, PearsonP, Pearson95CIlower, Pearson95CIupper)
dfStatistics <- rbind(dfStatistics, SummaryVec)
names(dfStatistics) <- c("SpearmanRho", "SpearmanP", "PearsonR", "PearsonP", "Pearson95CIlower", "Pearson95CIupper")
write.csv(dfStatistics, "Stats for RelChange between SuH and NICD correlation.csv")

rm(Spearman, SpearmanRho, SpearmanP, Pearson, PearsonR, PearsonP, 
   Pearson95CIlower, Pearson95CIupper, SummaryVec, dfStatistics)

setwd(mainDir)
#### end ####

