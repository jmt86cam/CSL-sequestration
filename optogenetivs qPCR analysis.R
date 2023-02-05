# clear the environment: 
rm(list = ls())
# clear the console: ctrl+L
# clear all the plots: 
dev.off()

#### Packages for this script ####
install.packages("tidyr")
install.packages("ggplot2")
install.packages("cowplot")
install.packages("forcats")
install.packages("dplyr")
library(tidyr)
library(ggplot2)
library(cowplot)
library(forcats)
library(dplyr)
library(ggsignif)
#### end ####

# Set where plots will be saved to:
setwd("/Users/jonathantownson/Documents/PhD/qPCR data/Optogenetic experiments")
# Set the path to the folder where the data csv files are stored:
DataPath <- "/Users/jonathantownson/Documents/PhD/qPCR data/Optogenetic experiments/Cleaned data for R"
# Instruct the script what the control primers being used are
ControlPrimers <- "RP49"
# Instruct the script of any Genotype ordering required e.g. control then experimental condition
XAxisOrder <- c("NoCRY Dark", "NoCRY 24hr Light", "Opto Dark", "Opto 1hr Light", "Opto 2hr Light", "Opto 8hr Light", "Opto 24hr Light")

#### Read in the data ####
## Create the folder names for the data
FolderFiles <- list.files(DataPath, pattern = "*.csv")
## Read in all the files in the folder
Files <- list()
# Comment out one of the below, first for loop only keeps files where StDev <1
# Second for loop keeps all files to be merged into raw data
for (f in 1:length(FolderFiles)){
  df <- read.csv(paste(DataPath,FolderFiles[f], sep ="/"))
  df <- dplyr::filter(df,Cp < 40)
  TechRepNum <- count(df, Name, Primers)
  ## Calculate mean and StDev of technical replicates
  MeanCp <- aggregate(df$Cp,by=list(Name=df$Name,Primers=df$Primers),data=df,FUN = mean)
  colnames(MeanCp)[colnames(MeanCp) == 'x'] <- 'Cp'
  StDevCp <- aggregate(df$Cp,by=list(Name=df$Name,Primers=df$Primers),FUN = sd)
  MeanCp <- merge(MeanCp,StDevCp,by=c("Name","Primers"))
  colnames(MeanCp)[colnames(MeanCp) == 'x'] <- 'StDevCp'
  MeanCp <- merge(MeanCp,TechRepNum,by=c("Name","Primers"))
  colnames(MeanCp)[colnames(MeanCp) == 'n'] <- 'TechRepNum'
  rm(StDevCp)
  rm(TechRepNum)
  rm(df)
  MeanCp$StDevCp[is.na(MeanCp$StDevCp)] <- 0
  ## Select data with StDev<1 and Cp less than 40
  Data <- dplyr::filter(MeanCp,StDevCp < 1)
  Data <- subset(Data, select = -StDevCp)
  Data$qPCRdate <- FolderFiles[[f]]
  Files[[f]] <- Data
  rm(MeanCp)
  rm(Data)
}
for (f in 1:length(FolderFiles)){
  df <- read.csv(paste(DataPath,FolderFiles[f], sep ="/"))
  df <- dplyr::filter(df,Cp < 40)
  TechRepNum <- count(df, Name, Primers)
  ## Calculate mean and StDev of technical replicates
  MeanCp <- aggregate(df$Cp,by=list(Name=df$Name,Primers=df$Primers),data=df,FUN = mean)
  colnames(MeanCp)[colnames(MeanCp) == 'x'] <- 'Cp'
  StDevCp <- aggregate(df$Cp,by=list(Name=df$Name,Primers=df$Primers),FUN = sd)
  MeanCp <- merge(MeanCp,StDevCp,by=c("Name","Primers"))
  colnames(MeanCp)[colnames(MeanCp) == 'x'] <- 'StDevCp'
  MeanCp <- merge(MeanCp,TechRepNum,by=c("Name","Primers"))
  colnames(MeanCp)[colnames(MeanCp) == 'n'] <- 'TechRepNum'
  rm(StDevCp)
  rm(TechRepNum)
  rm(df)
  MeanCp$StDevCp[is.na(MeanCp$StDevCp)] <- 0
  ## Select data with StDev<1 and Cp less than 40
  Data <- subset(MeanCp, select = -StDevCp)
  Data$qPCRdate <- FolderFiles[[f]]
  Files[[f]] <- Data
  rm(MeanCp)
  rm(Data)
} 
RawData <- bind_rows(Files)
## Get list of primers in dataset
PrimerList <- unique(RawData$Primers)
Genotypes <- unique(as.character(RawData$Genotype))
ExperimentPrimers <- setdiff(PrimerList,ControlPrimers)
## Create date and genotype columns from Name
RawData$NameTwo <- RawData$Name
## I name my samples in format DDMMYY_Genotype but others might not, the if statements should detect this accordingly.
if (grepl("_",RawData$NameTwo)) {
  RawData <- separate(RawData, NameTwo, into = c("NameThree", "NameFour") ,sep = "_")
  if (!is.na(as.numeric(RawData$NameThree))){
    RawData$Date <- RawData$NameThree
    RawData$Genotype <- RawData$NameFour
    RawData$NameThree <- NULL
    RawData$NameFour <- NULL
   } else {
    RawData$Date <- RawData$NameFour
    RawData$Genotype <- RawData$NameThree
    RawData$NameThree <- NULL
    RawData$NameFour <- NULL
   }
  NameGenotypeDate <- subset(RawData, select = c(Name, Genotype, Date))
  NameGenotypeDate <- NameGenotypeDate[!duplicated(NameGenotypeDate$Name), ]
} else if (grepl(" ",RawData$NameTwo)) {
  RawData <- separate(RawData, NameTwo, into = c("NameThree", "NameFour") ,sep = " ")
  if (!is.na(as.numeric(RawData$NameThree))){
    RawData$Date <- RawData$NameThree
    RawData$Genotype <- RawData$NameFour
    RawData$NameThree <- NULL
    RawData$NameFour <- NULL
  } else {
    RawData$Date <- RawData$NameFour
    RawData$Genotype <- RawData$NameThree
    RawData$NameThree <- NULL
    RawData$NameFour <- NULL
  }
  NameGenotypeDate <- subset(RawData, select = c(Name, Genotype, Date))
  NameGenotypeDate <- NameGenotypeDate[!duplicated(NameGenotypeDate$Name), ]
} else {
  RawData$Genotype <- RawData$NameTwo
  RawData$NameTwo <- NULL
  NameGenotypeDate <- subset(RawData, select = c(Name, Genotype))
  NameGenotypeDate <- NameGenotypeDate[!duplicated(NameGenotypeDate$Name), ]
}
# this line determines the order factors are plotted on the x axis, any not mentioned are determined automatically
RawData <- RawData %>% mutate(Genotype = fct_relevel(Genotype, XAxisOrder)) 
## Make each primer its own column
DifferentPrimers <- list()
for (i in 1:as.numeric(length(PrimerList))){
  df <- subset(RawData, RawData$Primers == PrimerList[i])
  DifferentPrimers[[i]] <- df
  rm(df)
}
for (i in 1:as.numeric(length(DifferentPrimers))){
  colnames(DifferentPrimers[[i]])[colnames(DifferentPrimers[[i]]) == "Cp"] <- unique(DifferentPrimers[[i]]$Primers)
  DifferentPrimers[[i]] <- subset(DifferentPrimers[[i]], select = -c(Primers,Date,Genotype,TechRepNum))
}
SeparatedPrimers <- DifferentPrimers[[1]]
for (i in 2:as.numeric(length(DifferentPrimers))){
  SeparatedPrimers <- merge(SeparatedPrimers,DifferentPrimers[[i]], by= c("Name","qPCRdate"), all = TRUE)
}
rm(DifferentPrimers)
RawDataSeparatedPrimers <- merge(SeparatedPrimers,NameGenotypeDate,by="Name")
rm(SeparatedPrimers)

## Calculate mean of samples across different qPCR runs, getting rid of any outliers
StDevCp <- aggregate(RawData$Cp,by=list(Name=RawData$Name,Primers=RawData$Primers),FUN = sd)
StDevCp$x[is.na(StDevCp$x)] <- 0
df <- merge(StDevCp,RawData,by= c("Name","Primers"), all = TRUE)
Sb1 <- subset(df,x<1)
Sb2 <- subset(df,x>=1)
Sb2 <- subset(Sb2, TechRepNum == 2)
Sb2$x <- NULL
StDevSb2 <- aggregate(Sb2$Cp,by=list(Name=Sb2$Name,Primers=Sb2$Primers),FUN = sd)
StDevSb2$x[is.na(StDevSb2$x)] <- 0
Sb2 <- merge(StDevSb2,Sb2,by= c("Name","Primers"), all = TRUE)
rm(StDevSb2,df,StDevCp)
Sb2 <- dplyr::filter(Sb2, x < 1)
UsefulData <- bind_rows(Sb1,Sb2)
rm(Sb1,Sb2)
UsefulData <- subset(UsefulData, select = -c(x,TechRepNum,qPCRdate))
UsefulData <- aggregate(UsefulData$Cp,by=list(Name=UsefulData$Name,Primers=UsefulData$Primers),data=UsefulData,FUN = mean)
colnames(UsefulData)[colnames(UsefulData) == 'x'] <- 'MeanCp'
UsefulData <- merge(UsefulData,NameGenotypeDate,by="Name")

## Make each primer its own column
DifferentPrimers <- list()
for (i in 1:as.numeric(length(PrimerList))){
  df <- subset(UsefulData, UsefulData$Primers == PrimerList[i])
  DifferentPrimers[[i]] <- df
  rm(df)
}
for (i in 1:as.numeric(length(DifferentPrimers))){
  colnames(DifferentPrimers[[i]])[colnames(DifferentPrimers[[i]]) == "MeanCp"] <- unique(DifferentPrimers[[i]]$Primers)
  DifferentPrimers[[i]] <- subset(DifferentPrimers[[i]], select = -c(Primers,Date,Genotype))
}
SeparatedPrimers <- DifferentPrimers[[1]]
for (i in 2:as.numeric(length(DifferentPrimers))){
  SeparatedPrimers <- merge(SeparatedPrimers,DifferentPrimers[[i]], by= "Name", all = TRUE)
}
rm(DifferentPrimers)
UsefulDataSeparatedPrimers <- merge(SeparatedPrimers,NameGenotypeDate,by="Name")
rm(SeparatedPrimers)

## Convert Cp values to relative mRNA values
Data <- UsefulData
Data$RelativemRNA <- 2^-Data$MeanCp
Data <- subset(Data, select = -MeanCp)
DifferentPrimers <- list()
for (i in 1:as.numeric(length(PrimerList))){
  df <- subset(Data, Data$Primers == PrimerList[i])
  DifferentPrimers[[i]] <- df
  rm(df)
}
for (i in 1:as.numeric(length(DifferentPrimers))){
  colnames(DifferentPrimers[[i]])[colnames(DifferentPrimers[[i]]) == "RelativemRNA"] <- unique(DifferentPrimers[[i]]$Primers)
  DifferentPrimers[[i]] <- subset(DifferentPrimers[[i]], select = -c(Primers,Date,Genotype))
}
SeparatedPrimers <- DifferentPrimers[[1]]
for (i in 2:as.numeric(length(DifferentPrimers))){
  SeparatedPrimers <- merge(SeparatedPrimers,DifferentPrimers[[i]], by= "Name", all = TRUE)
}
rm(DifferentPrimers)
DataSeparatedPrimers <- merge(SeparatedPrimers,NameGenotypeDate,by="Name")
rm(SeparatedPrimers)

## Normalise data with control target genes
# Create two dataframes named after the control primers
NormalisedData <- list()
for (i in 1:as.numeric(length(ControlPrimers))){
  df <- subset(DataSeparatedPrimers, select = c(Name, Date, Genotype))
  NormalisedData[[i]] <- df
  assign(ControlPrimers[i],NormalisedData[[i]])
  rm(df)
}
rm(NormalisedData)
# To each control primer dataframe add the experiment primer mRNA values and normalise by dividing with the control primer mRNA value
for (i in 1:as.numeric(length(ControlPrimers))){
  df <- get(ControlPrimers[i])
  ControlValues <- DataSeparatedPrimers[[ControlPrimers[i]]]
  for (k in 1:as.numeric(length(ExperimentPrimers))){
    df[[ExperimentPrimers[k]]] <- DataSeparatedPrimers[[ExperimentPrimers[k]]]/ControlValues
  }
  df <- df %>% mutate(Genotype = fct_relevel(Genotype, XAxisOrder)) # this line determines the order factors are plotted on the x axis, any not mentioned are determined automatically
  assign(ControlPrimers[i],df)
  rm(df)
  rm(ControlValues)
}
rm(NameGenotypeDate)
#### end ####

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

VariableList <- c("m3","mBeta")
DataStatSumTestListFun <- function(DataList,Filename){
  OutputList <- list()
  for (p in 1:length(DataList)) {
    df <- na.omit(DataList[[p]])
    dfStatistics <- data.frame()
    for (g in 1:length(unique(df$Genotype))) {
      for (v in 1:length(VariableList)) {
        Testdf <- dplyr::filter(df,Genotype == unique(df$Genotype)[g])
        MeanValue <- mean(Testdf[[VariableList[v]]])
        MedianValue <- median(Testdf[[VariableList[v]]])
        VarValue <- var(Testdf[[VariableList[v]]])
        StandardDeviation <- sd(Testdf[[VariableList[v]]])
        StandardError <- StandardDeviation/sqrt(length(Testdf[[VariableList[v]]]))
        result <- shapiro.test(Testdf[[VariableList[v]]])
        SWpvalue <- result[["p.value"]]
        StatsResults <- c(unique(df$Genotype)[g],VariableList[v],MeanValue,MedianValue,VarValue,StandardDeviation,StandardError,SWpvalue)
        dfStatistics <- rbind(dfStatistics,StatsResults)
        rm(Testdf,result,MeanValue,MedianValue,VarValue,StandardDeviation,StandardError,SWpvalue,StatsResults)
      }
    }
    names(dfStatistics) <- c("Genotype","Variable","Mean","Median","Variance","Standard Deviation","Standard Error","Shapiro-Wilk p-value")
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

Stats_tests <- function(DataList, Filename, Gene){
  dfStatistics <- data.frame()
  for (p in 1:(length(DataList))) {
    df1 <- na.omit(DataList[[p]])
    for (q in 1:(length(DataList))) { 
      df2 <- na.omit(DataList[[q]])
      Comparison <- paste(unique(df1$Genotype), "vs", unique(df2$Genotype))
      StudentT <- t.test(df1[,Gene], df2[,Gene], var.equal = TRUE)
      StudentT <- StudentT[["p.value"]]
      WelchT <- t.test(df1[,Gene], df2[,Gene], var.equal = FALSE)
      WelchT <- WelchT[["p.value"]]
      MannWhitney <- wilcox.test(df1[,Gene],df2[,Gene])
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
#### end ####

#### Bar graphs with/out statistics and saved ####
setwd("/Users/jonathantownson/Documents/PhD/CSL localisation Paper/Figures for paper/Files/Fig 4")

df <- filter(RP49, Genotype == "Opto Dark" | Genotype == "Opto 24hr Light")
levels(df$Genotype) <- gsub("Opto Dark", "Dark", levels(df$Genotype))
levels(df$Genotype) <- gsub("Opto 24hr Light", "Light", levels(df$Genotype))
df1list <- list()
for (g in 1:length(unique(df$Genotype))) {
  df1list[[g]] <- df[str_detect(df$Genotype, as.character(unique(df$Genotype)[g])),]
  names(df1list)[g] <- as.character(unique(df$Genotype)[g])
}
ggplot(df, aes(x=Genotype, y=m3, fill = Genotype)) + theme_jmt() + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) + theme(axis.title.x=element_blank()) +  
  stat_summary(fun.y = mean, geom = "bar",show.legend = FALSE) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.4, colour = "grey50") +
  geom_signif(comparisons = list(c("Dark","Light")), y_position = 0.035,
              map_signif_level=TRUE, test = "t.test", color = "black") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = expression(paste("Expression of ", italic("E(spl)m3"))), y = expression(paste("Relative ", italic("m3"), " mRNA level")))
ggsave("m3 expression of optogenetics.jpg", device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

ggplot(df, aes(x=Genotype, y=mBeta, fill = Genotype)) + theme_jmt() + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) + theme(axis.title.x=element_blank()) + 
  stat_summary(fun.y = mean, geom = "bar",show.legend = FALSE) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.4, colour = "grey50") +
  geom_signif(comparisons = list(c("Dark","Light")), y_position = 0.17,
              map_signif_level=TRUE, test = "t.test", color = "black") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = expression(paste("Expression of ", italic("E(spl)mβ"))), y = expression(paste("Relative ", italic("mβ"), " mRNA level")))
ggsave("mBeta expression of optogenetics.jpg", device = "jpeg", dpi = "retina",
       width = 15, height = 15, units = "cm")

NormDis_plots(df1list, "m3", "Genotype", "Normal distribution of m3")
QQ_plots(df1list, "m3", "Genotype","QQ plots of m3")
DataStatSumTestListFun(df1list,"Summary of m3")
Stats_tests(df1list,"SigStats of m3","m3")

NormDis_plots(df1list, "mBeta", "Genotype", "Normal distribution of mBeta")
QQ_plots(df1list, "mBeta", "Genotype","QQ plots of mBeta")
DataStatSumTestListFun(df1list,"Summary of mBeta")
Stats_tests(df1list,"SigStats of mBeta","mBeta")

#### end ####