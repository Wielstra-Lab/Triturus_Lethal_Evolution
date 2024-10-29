args_in <- commandArgs(trailingOnly = TRUE)

library(ggplot2)

target <- args_in[1]
output <- args_in[2]
peakplot <- paste(output,"Peak_Histogram.png", sep = "")
peak100plot <- paste(output,"Peak_Region_Histogram.png", sep = "")

BigDepth <- read.table(target, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

MainFrame <- data.frame(Target=character(),PeakCoverage=integer(),PeakRegion=integer())

StartSeq <- 1
Peak100 <- 0
Peak <- 0
OldName <- BigDepth[1,1]

for(x in 1:nrow(BigDepth)){
  
  if(BigDepth[x,1] == OldName){
    
    if(BigDepth[x,3] > Peak){
      
      Peak <- BigDepth[x,3]
    }
    
    if(x - StartSeq > 98){
    
      CovA <- BigDepth[x,3]
      CovB <- BigDepth[x - 49, 3]
      CovC <- BigDepth[x - 99, 3]
      RegMin <- min(CovA,CovB,CovC)
      if(RegMin > Peak100){Peak100 <- RegMin}
    }
  }
  else{
    
    NewFrame <- data.frame(Target=OldName,PeakCoverage=Peak,PeakRegion=Peak100)
    MainFrame <- rbind(MainFrame,NewFrame)
    
    StartSeq <- x
    OldName <- BigDepth[x,1]
    Peak <- BigDepth[x,3]
    Peak100 <- 0
  }
}
NewFrame <- data.frame(Target=OldName,PeakCoverage=Peak,PeakRegion=Peak100)
MainFrame <- rbind(MainFrame,NewFrame)

write.table(MainFrame, file = output, row.name = FALSE)

MedPeak <- median(MainFrame[,2])
MeanPeak <- mean(MainFrame[,2])
MedPeak100 <- median(MainFrame[,3])
MeanPeak100 <- mean(MainFrame[,3])
NoTar <- nrow(MainFrame)

print(paste("Input file:", target))
print(paste("Number of target sequences with at least on read:", NoTar))
print(paste("Mean peak coverage:", MeanPeak))
print(paste("Median peak coverage:", MedPeak))
print(paste("Mean peak 100 bp region coverage:", MeanPeak100))
print(paste("Median peak 100 bp region  coverage:", MedPeak100))
print("")
print("IGNORE THESE WARNINGS FOR NOW")
print("")

Peak_Lim <- quantile(MainFrame[,2], 0.95) * 1.5
Peak_width <- floor(Peak_Lim / 50) + 1
Peak100_Lim <- quantile(MainFrame[,3], 0.95) * 1.5
Peak100_width <- floor(Peak100_Lim / 50) + 1

PeakHist <- ggplot(data = MainFrame, aes(x = PeakCoverage))
PeakHist + geom_histogram(binwidth = Peak_width, color = "Darkslategrey", fill = "Darkslategrey", alpha = 0.5) + xlim(0,Peak_Lim) + theme_minimal()
ggsave(peakplot)

Peak100Hist <- ggplot(data = MainFrame, aes(x = PeakRegion))
Peak100Hist + geom_histogram(binwidth = Peak100_width, color = "Darkslategrey", fill = "Darkslategrey", alpha = 0.5) + xlim(0,Peak100_Lim) + theme_minimal()
ggsave(peak100plot)


