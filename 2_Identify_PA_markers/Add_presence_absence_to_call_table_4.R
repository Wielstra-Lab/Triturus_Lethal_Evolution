### Script takes a table of coverage for a linkage map family and converts coverage/presence absence data, into psuedoSNPs in the parent.call format used by LepMAP 3

## This version is modified to use whole sample set for projecting average coverage - hopefully should make it a little more accurate

args_in <- commandArgs(trailingOnly = TRUE)

### First load coverage table and list of Presence/Absence markers

input_file <- args_in[1]
PA_marker_file <- args_in[2]
output_file <- args_in[3]

Coverage_table_full <- read.table(input_file,stringsAsFactors = F)

PA_marker_table <- read.table(PA_marker_file)
PA_marker_List <- PA_marker_table$V1

### created expected coverage table based on marker and sample average coverage

Sample_ave_cov <- as.numeric(colMeans(Coverage_table_full))

Sample_set_ave_cov <- mean(Sample_ave_cov)

Marker_ave_cov <- as.numeric(rowMeans(Coverage_table_full))

Expected_coverage_table_full <- Coverage_table_full

for(i in 1:ncol(Expected_coverage_table_full)){
  Expected_coverage_table_full[,i] <- as.integer(Marker_ave_cov * (Sample_ave_cov[i] / Sample_set_ave_cov))
}

### subset the table for only the markers we're interested in

Coverage_table <- Coverage_table_full[PA_marker_List,]
Expected_coverage_table <- Expected_coverage_table_full[PA_marker_List,]

### subtract 0.1 of expected coverage away from measured coverage, with min of 0 
### (removes effects of mismapped reads, contamination etc.)

Processed_cov_table <- Coverage_table - (Expected_coverage_table * 0.1)

Homo_screen_table <- (Coverage_table - (Expected_coverage_table * 1.2)) / Expected_coverage_table

### use these codes for presence and absence 

Absence_code <- "1.0 0 0 0 0 0 0 0 0 0"
Presence_code <- "0 1.0 0 0 0 0 0 0 0 0"
Homo_code <- "0 0 0 0 1.0 0 0 0 0 0"

#Absence_code <- "0"
#Presence_code <- "1"
#Homo_code <- "2"

### Make a LepMAP format table of two columns (CHR and POS) and the same number of rows as the number of makers

PA_marker_format <- Processed_cov_table[,1:2]
colnames(PA_marker_format) <- c("CHR","POS")

### Combine the format table and the processed coverage table

PA_marker_calls <- cbind(PA_marker_format,Processed_cov_table)

### loop over ever marker in the table

for(i in 1:nrow(PA_marker_calls)){
  
  ### CHR is the marker name and POS is set to 100
  
  PA_marker_calls[i,1] <- row.names(PA_marker_calls)[i]
  PA_marker_calls[i,2] <- 100
  
  ### Loop over every sample, if the sample has 0 coverage for this marker replace with absence code, otherwise, presence code
  for(j in 3:ncol(PA_marker_calls)){
    if(as.numeric(PA_marker_calls[i,j]) <= 0){
      PA_marker_calls[i,j] <- Absence_code
    }
    else if(as.numeric(Homo_screen_table[i,j - 2]) >= 0){
      Homo_prob <- min(1, as.numeric(Homo_screen_table[i,j - 2]))
      Het_prob <- 1 - Homo_prob
      PA_marker_calls[i,j] <- paste("0", Het_prob, "0 0", Homo_prob, "0 0 0 0 0", sep = " ")
      }
    
    else{PA_marker_calls[i,j] <- Presence_code}
  }
}

### All done, output the table without row or col names

write.table(PA_marker_calls,output_file,sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE) 
