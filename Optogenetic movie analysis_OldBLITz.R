# clear the environment: 
rm(list = ls())
# clear the console: ctrl+L
# clear all the plots: 
dev.off()

#### Packages for this script ####
install.packages("ggplot2")
install.packages("dplyr")
install.packages("tibble")
library(ggplot2)
library(cowplot)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(tibble)
#### end ####

#### Data analysis functions ####
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
#### end ####

# Set location of data and where plots will be saved to:
mainDir <- "/Users/jonathantownson/Documents/PhD/Images/Analysed image data/Optogenetics data/R analysis/OldBlitz/mCherry"
plotDir <- "/Users/jonathantownson/Documents/PhD/CSL localisation Paper/Figures for paper/Files/Fig 2"
setwd(mainDir )
# Create the path where the data is stored
DataPath <- paste(mainDir,"Cleaned csv files",sep = "/")
# Provided the Background files folder is the only one in DataPath, create a path to the background files
BackgroundFilesPath <- paste(mainDir,"Background files",sep = "/")
# Provide the length of the experiment in minutes
TimeMin <- 10

#### Read in the raw data ####
RawData <- list()
Files <- list.files(DataPath, pattern = "*.csv")
for (f in 1:length(Files)) {
  # Create Filepath for each gland
  FilePath <- paste(DataPath,Files[f],sep = "/")
  
  # Extract the file conditions in the name of the file. These are added to the
  # data in lines 72-78
  Filename <- str_sub(Files[f], end = -5)
  FileConds <- str_split(Filename, "_")[[1]]
  
  # Read in the data for each gland
  df <- read.csv(FilePath)
  names(df)[names(df) == 'X'] <- 'Frame'
  
  # Read in the background values for the gland
  BGfilename <- paste(Filename,"BG.csv")
  BGvalues <- read.csv(paste(BackgroundFilesPath,BGfilename,sep = "/"))
  BGvalues$X <- NULL
  for (c in 1:length(colnames(BGvalues))) {
    names(BGvalues)[names(BGvalues) == colnames(BGvalues)[c]] <- paste("BG", str_sub(colnames(BGvalues)[c], end = -2))
  }
  rm(BGfilename)
  
  # Work out how many cells were measured in the gland, here we subtract 1 because the frame column the same for all cells
  NumberOfCells <- (length(unique(parse_number(colnames(df)[-1]))))/2
  
  # Make each gland a list of cells within it
  Gland <- list()
  h <- 1
  for (g in 1:NumberOfCells) {
    # Rename each column name with Cell/Nucleus + measurement + cellnumber
    df1 <- df
    df1$Frame <- NULL
    CellValues <- df1[,c(h:(h+5))]
    CellValuesColNames <- colnames(CellValues)
    for (c in 1:length(CellValuesColNames)) {
      names(CellValues)[names(CellValues) == CellValuesColNames[c]] <- paste("Cell", str_extract(CellValuesColNames[c], "^\\D+"))
    }
    
    NucleusValues <- df1[,c((h+6):(h+11))]
    NucleusValuesColNames <- colnames(NucleusValues)
    for (c in 1:length(NucleusValuesColNames)) {
      names(NucleusValues)[names(NucleusValues) == NucleusValuesColNames[c]] <- paste("Nucleus", str_extract(CellValuesColNames[c], "^\\D+"))
    }
    
    # Combine the cell, nucleus and background values into one dataframe for each
    # cell in the gland
    Cell <- cbind(CellValues, NucleusValues, BGvalues)
    rm(CellValues,NucleusValues,CellValuesColNames,NucleusValuesColNames)
    
    # Add the frame number to the dataframe
    Cell$Frame <- df$Frame
    Cell$CellNumber <- as.numeric(g)
    # Add columns for the file conditions extracted earlier
    Cell$Date <- FileConds[1]
    Cell$Construct <- FileConds[2]
    Cell$Gland <- parse_number(FileConds[3])
    Cell$ArgonPercent <- parse_number(FileConds[4])
    Cell$LaserActivation <- FileConds[5]
    Cell$Minutes <- parse_number(FileConds[6])
    
    # Add each cell dataframe to the list of cells in the gland
    Gland[[g]] <- Cell
    rm(Cell)
    h <- h+12
  }
  
  # Add each gland (aka file) to the raw data list
  RawData[[f]] <- Gland
  # Rename the gland the filename
  names(RawData)[f] <- Filename
  
  rm(df, Gland, BGvalues, FileConds, NumberOfCells, Filename, FilePath, h, df1)
}
rm(c,f,g)
#### end ####

#### Perform calculations and keep the useful bits of data ####
UsefulData <- list()
for (g in 1:length(RawData)) {
  Gland <- RawData[[g]]
  UsefulGland <- list()
  for (c in 1:length(Gland)) {
    df <- Gland[[c]]
    # Calculate the signal coming from the cytoplasm (effectively subtracting the nuclear signal from cell signal)
    df$'Cytoplasm Mean' <- (df$`Cell IntDen`-df$`Nucleus IntDen`)/(df$`Cell Area`-df$`Nucleus Area`)
    # Keep only columns that will be useful in the future
    Usefuldf <- select(df,c(Frame:Minutes, 'Nucleus Mean', 'Cytoplasm Mean', 'BG Mean'))
    # Subtract the background and remove BG column
    Usefuldf$`Nucleus Mean` <- Usefuldf$`Nucleus Mean` - Usefuldf$`BG Mean`
    Usefuldf$`Cytoplasm Mean` <- Usefuldf$`Cytoplasm Mean` - Usefuldf$`BG Mean`
    Usefuldf$`BG Mean` <- NULL
    ## Calculate means of tracking nuclear fluorescence
    # with NC ratio
    Usefuldf$NCratio <- Usefuldf$`Nucleus Mean` / Usefuldf$`Cytoplasm Mean`
    Usefuldf$NCratioFoldChange <- Usefuldf$NCratio / Usefuldf$NCratio[1]
    Usefuldf$NCratioSubStart <- Usefuldf$NCratio - Usefuldf$NCratio[1]
    # as nuclear accumulation
    Usefuldf$NuclearAccumulation <- Usefuldf$`Nucleus Mean` / Usefuldf$`Nucleus Mean`[1]
    Usefuldf$NuclearCumulativeIncrease <- Usefuldf$`Nucleus Mean` - Usefuldf$`Nucleus Mean`[1]
    # as nuclear increment
    Usefuldf$NuclearIncrement[1] <- (Usefuldf$`Nucleus Mean`[1] - Usefuldf$`Nucleus Mean`[1])/Usefuldf$`Nucleus Mean`[1]
    for (r in 2:length(Usefuldf$`Nucleus Mean`)) {
      Usefuldf$NuclearIncrement[r] <- (Usefuldf$`Nucleus Mean`[r] - Usefuldf$`Nucleus Mean`[r-1])/Usefuldf$`Nucleus Mean`[1]
    }
    UsefulGland[[c]] <- Usefuldf
    rm(df, Usefuldf, r)
  }
  UsefulData[[g]] <- UsefulGland
  names(UsefulData)[g] <- names(RawData)[g]
  rm(Gland, UsefulGland, c)
}
rm(g)
#### end ####

#### Combine all the data into one data frame ####
DataInFiles <- list()
for (g in 1:length(UsefulData)) {
  File <- UsefulData[[g]]
  File <- do.call(rbind.data.frame,File)
  DataInFiles[[g]] <- File
  names(DataInFiles)[g] <- names(UsefulData[g])
  rm(File)
}
AllData <- do.call(rbind.data.frame,DataInFiles)
rm(g, DataInFiles)

DataByActivation <- list()
for (g in 1:length(unique(AllData$LaserActivation))) {
  df <- dplyr::filter(AllData, LaserActivation == unique(AllData$LaserActivation)[g])
  dfVariables <- select(df,`Nucleus Mean`:NuclearIncrement)
  df <- select(df,c(Frame, CellNumber))
  df$Name <- paste(as.data.frame(str_split_fixed(rownames(df), fixed("."), 2))$V1," Cell",df$CellNumber)
  df <- select(df,-CellNumber)
  Data <- list()
  for (d in 1:length(colnames(dfVariables))) {
    df1 <- cbind(select(df,c(Frame,Name)),select(dfVariables,colnames(dfVariables)[d]))
    df1 <- df1 %>% spread(Name,colnames(dfVariables)[d])
    df2 <- as.data.frame(df1$Frame)
    df1$Frame<-NULL
    df2$Mean <- rowMeans(df1)
    df2$SEM <- SEMrows(df1)
    names(df2)[names(df2) == 'df1$Frame'] <- 'Frame'
    names(df2)[names(df2) == 'Mean'] <- paste("Mean",colnames(dfVariables)[d])
    names(df2)[names(df2) == 'SEM'] <- paste("SEM",colnames(dfVariables)[d])
    Data[[d]] <- df2
    names(Data)[d] <- colnames(dfVariables)[d]
    rm(df1,df2)
  }
  FinalDF <- Data[[1]]
  for (f in 2:length(Data)) {
    Nextdf <- Data[[f]]
    FinalDF <- merge(FinalDF, Nextdf, by = "Frame")
  }
  DataByActivation[[g]] <- FinalDF
  names(DataByActivation)[g] <- unique(AllData$LaserActivation)[g]
  rm(df,dfVariables,Data,FinalDF,d,f,Nextdf)
}
AllSummarisedData <- do.call(rbind.data.frame,DataByActivation)
AllSummarisedData$LaserActivation <- as.data.frame(str_split_fixed(rownames(AllSummarisedData), fixed("."), 2))$V1
rm(g)
#### end ####

#### Calculate the mean and SEM for each activated gland ####
SummarisedGlands <- list()
for (g in 1:length(UsefulData)) {
  Gland <- UsefulData[[g]]
  df <- Gland[[1]]
  SummarisedDF <- select(df,c(Frame,Date:Minutes))
  MeasurementsOnly <- select(df,-c(Date:Minutes))
  for (c in 2:length(Gland)) {
    df <- Gland[[c]]
    df <- select(df,-c(Date:Minutes))
    MeasurementsOnly <- rbind(MeasurementsOnly,df)
  }
  rm(df,c)
  
  for (m in 3:length(colnames(MeasurementsOnly))) {
    df <- select(MeasurementsOnly,c(Frame, CellNumber,colnames(MeasurementsOnly)[m]))
    df <- df %>% spread(CellNumber, colnames(MeasurementsOnly)[m])
    MeanDF <- as.data.frame(df$Frame)
    df$Frame <- NULL
    MeanDF$Mean <- rowMeans(df)
    MeanDF$SEM <- SEMrows(df)
    names(MeanDF)[names(MeanDF) == 'Mean'] <- paste(colnames(MeasurementsOnly)[m] ,"Mean")
    names(MeanDF)[names(MeanDF) == 'SEM'] <- paste(colnames(MeasurementsOnly)[m] ,"SEM")
    MeanDF$`df$Frame` <- NULL
    SummarisedDF <- cbind(SummarisedDF,MeanDF)
  }
  SummarisedGlands[[g]] <- SummarisedDF
  rm(df, SummarisedDF, m, MeanDF, MeasurementsOnly)
  names(SummarisedGlands)[g] <- names(UsefulData)[g]
}
rm(Gland, g)

SummarisedGlands <- do.call(rbind,SummarisedGlands)
SummarisedGlands <- rownames_to_column(SummarisedGlands)
names(SummarisedGlands)[names(SummarisedGlands) == 'rowname'] <- "File"
#### end ####

setwd(plotDir)

#### Comparing laser power ####
df <- filter(AllSummarisedData,LaserActivation %in% c("458nm10","458nm5","458nm20"))
df$'Time (minutes)' <- df$Frame*(TimeMin/(length(df$Frame)/length(unique(df$LaserActivation))))
df$LaserActivation <- gsub("458nm10", "1x", df$LaserActivation)
df$LaserActivation <- gsub("458nm20", "2x", df$LaserActivation)
df$LaserActivation <- gsub("458nm5", "0.5x", df$LaserActivation)
names(df)[names(df) == 'LaserActivation'] <- 'Power'

ggplot(df, aes(`Time (minutes)`,`Mean NuclearAccumulation`,group = Power, col=Power,shape=Power,
                              ymin=`Mean NuclearAccumulation`-`SEM NuclearAccumulation`,ymax=`Mean NuclearAccumulation`+`SEM NuclearAccumulation`)) + theme_jmt() +
  geom_ribbon(linetype = "dashed", fill = "grey90",alpha = 0.5, size = 0.25) + geom_line(size = 0.5) +
  labs(title = "", x = "Time (minutes)", y = "Nuclear Accumulation") +
  scale_color_brewer(palette = "Set2")
ggsave(paste("Nuclear Accumulation of different powers.jpg"), device = "jpeg", 
       width = 12, height = 12, units = "cm", dpi = "retina")

rm(df)
#### end ####

#### Comparing wavelengths ####
df <- filter(AllSummarisedData,LaserActivation %in% c("458nm10","476nm10","488nm10","496nm10","NoActivation"))
df$'Time (minutes)' <- df$Frame*(TimeMin/(length(df$Frame)/length(unique(df$LaserActivation))))
df$LaserActivation <- gsub("458nm10", "458nm", df$LaserActivation)
df$LaserActivation <- gsub("476nm10", "476nm", df$LaserActivation)
df$LaserActivation <- gsub("488nm10", "488nm", df$LaserActivation)
df$LaserActivation <- gsub("496nm10", "496nm", df$LaserActivation)
df$LaserActivation <- gsub("NoActivation", "Dark control", df$LaserActivation)
names(df)[names(df) == 'LaserActivation'] <- 'Wavelength'

ggplot(df, aes(`Time (minutes)`,`Mean NuclearAccumulation`,group = Wavelength, col=Wavelength,shape=Wavelength,
               ymin=`Mean NuclearAccumulation`-`SEM NuclearAccumulation`,ymax=`Mean NuclearAccumulation`+`SEM NuclearAccumulation`)) + theme_jmt() +
  geom_ribbon(linetype = "dashed", fill = "grey90",alpha = 0.5, size = 0.25, show.legend = FALSE) + geom_line(size = 0.5) +
  guides(color=guide_legend(override.aes=list(fill=NA, size = 2.5))) +
  labs(title = "OptIC-Notch activation", x = "Time (mins)", y = "Nuclear Accumulation") +
  scale_color_brewer(palette = "Set2")
ggsave(paste("Nuclear Accumulation of different wavelengths.jpg"), device = "jpeg", 
       width = 16, height = 16, units = "cm", dpi = "retina")

rm(df)
#### end ####
