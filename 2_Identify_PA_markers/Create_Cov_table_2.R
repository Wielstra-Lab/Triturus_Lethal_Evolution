args_in <- commandArgs(trailingOnly = TRUE)

### Create coverage table: takes the .PeakCoverage files created by Peakloop2.R and 
### builds a table of the coverage of every defined target in every define sample
### (coverage is taken as the best covered 100bp region in each target) 

target_list_file <- args_in[1]
sample_list_file <- args_in[2]
sample_dir <- args_in[3]
output_name <- args_in[4]

Master_Table <- read.table(target_list_file)

Sample_list_table <- read.table(sample_list_file)
Sample_List <- Sample_list_table$V1

print(Sample_List)

for (i in 1:length(Sample_List)) {
  
  Sample_Path <- paste(sample_dir,"/",Sample_List[i],".PeakCoverage", sep = "")
  Sample_Code <- substr(Sample_List[i],1,7)
  
  Sample_Table <- read.table(Sample_Path, header = TRUE)
  
  for (irow in 1:nrow(Sample_Table)) {
    Sample_Table[irow,1] <- sub("TRINITY_", "",Sample_Table[irow,1])
  }
  
  Master_Table[Sample_Code] <- NA
  x <- 0
  
  for (j in 1:nrow(Master_Table)) {
    if (Master_Table[j,1] == Sample_Table[(j-x),1]) {
      Master_Table[j,Sample_Code] <- Sample_Table[(j-x),3]
    }
    else {
      Master_Table[j,Sample_Code] <- 0
      x <- x + 1
    }
  }
  
  print(paste("done with: ", Sample_Code, sep = ""))
}

Master_Table_2 <- Master_Table[,-1]
rownames(Master_Table_2) <- Master_Table[,1]

write.table(Master_Table_2, output_name)