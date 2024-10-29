### Script written to find up to two classes of tightly linked presence absence markers based on coverage data from target capture

args_in <- commandArgs(trailingOnly = TRUE)

Coverage_file <- args_in[1]
Output <- args_in[2]

### First load in the table of coverage values per marker and per sample

Coverage_table <- read.table(Coverage_file, stringsAsFactors = F)

### Then calculate the average coverage per sample and per marker (stored as vectors), and the average of the averages

Sample_ave_cov <- as.numeric(colMeans(Coverage_table))

Sample_set_ave_cov <- mean(Sample_ave_cov)

Marker_ave_cov <- as.numeric(rowMeans(Coverage_table))

### For all positions multiply the sample mean (divided by the total average) by the marker mean to calculate the expected coverage

Expected_coverage_table <- Coverage_table

for(i in 1:ncol(Expected_coverage_table)){
  Expected_coverage_table[,i] <- as.integer(Marker_ave_cov * (Sample_ave_cov[i] / Sample_set_ave_cov))
}

### Subtract the measured values from the expected values to create a table of differences
### Subtract 10 to suppress the effect of low coverage in the next step


Difference_raw <- ((Expected_coverage_table - Coverage_table) - 10)

### Create a table of the difference divided by the measured coverage
### Effectively this normalizes for different coverage between markers

Difference_pro <- (Difference_raw / (Coverage_table + 0.1))

### Calculate magnitude of the top 10th percentile of the normalized differences per marker

for(i in 1:nrow(Difference_pro)){
  Difference_pro[i,"Tentiles"] <- quantile(as.numeric(Difference_pro[i,]), probs = 0.9, na.rm = TRUE)
}

### Select for only markers where more that 10% show a significant difference (0.5) between measured and expected coverage

Difference_final <- subset(Difference_pro,Tentiles > 0.5)

### Get the names of these markers and make table of their measured coverage

Candidate_markers <- row.names(Difference_final)

Candidate_coverage <- Coverage_table[Candidate_markers,]

### subtract 20% of the marker average coverage from each value (min = 0)

Candidate_ave_coverage_0.1 <- as.integer(rowMeans(Candidate_coverage) / 5)

for(i in 1:ncol(Candidate_coverage)){
  Candidate_coverage[,i] <- pmax(Candidate_coverage[,i] - Candidate_ave_coverage_0.1, 0)
}

### Count the number of samples with no read for each marker

Zero_scores <- rowSums(Candidate_coverage == 0, na.rm = TRUE)

### Filter for markers which are absent in at least 25% of all samples

Filtered_candidate_coverage <- subset(Candidate_coverage, rowSums(Candidate_coverage == 0, na.rm = TRUE) > ncol(Candidate_coverage) / 4)

### Identify list of type 1 and type 2 markers

Type_1_markers <- c()
Not_type_1_markers <- c()

# Define the first marker as a type 1 marker

Type_1_markers <- c(Type_1_markers, rownames(Filtered_candidate_coverage)[1])

# Loop through every marker after the first

for( marker in 2:nrow(Filtered_candidate_coverage)){
  
  #count the number of times the coverage of a sample in marker N matches the coverage of that sample in the first marker
  
  Matches <- 0
  
  for(samp in 1:ncol(Filtered_candidate_coverage)){
    if(Filtered_candidate_coverage[1,samp] == Filtered_candidate_coverage[marker,samp]){
      Matches <- Matches + 1
    }
  }
  
  # If the number of matches is greater than the number of samples / 4, add it to the list of type 1 markers
  # If not add to the list of Not_class_1_markers
  
  if(Matches > ncol(Candidate_coverage) / 4){
    Type_1_markers <- c(Type_1_markers, rownames(Filtered_candidate_coverage)[marker])
  }
  else{
    Not_type_1_markers <- c(Not_type_1_markers, rownames(Filtered_candidate_coverage)[marker])
  }
}

### Create table of all the non type 1 markers, then test to see if they fit type 2

Not_type_1_table <- Filtered_candidate_coverage[Not_type_1_markers,]

Type_2_markers <- c()
Not_type_2_markers <- c()

Type_2_markers <- c(Type_2_markers, rownames(Not_type_1_table)[1])

for( marker in 2:nrow(Not_type_1_table)){
  
  Matches <- 0
  
  for(samp in 1:ncol(Not_type_1_table)){
    if(Not_type_1_table[1,samp] == Not_type_1_table[marker,samp]){
      Matches <- Matches + 1
    }
  }
  
  if(Matches > ncol(Candidate_coverage) / 4){
    Type_2_markers <- c(Type_2_markers, rownames(Not_type_1_table)[marker])
  }
  else{
    Not_type_2_markers <- c(Not_type_2_markers, rownames(Not_type_1_table)[marker])
  }
}

### Make tables for markers

Type_1_table <- Filtered_candidate_coverage[Type_1_markers,]
Type_2_table <- Filtered_candidate_coverage[Type_2_markers,]
Anomolous_table <- Filtered_candidate_coverage[Not_type_2_markers,]

### Make lists of samples, for all 3 possible genotypes

Type_2_homo_list <- c()
Type_1_homo_list <- c()
Hetero_list <- c()

### Loop of every sample, for each category of marker, if a sample has 0 coverage for more than half the markers, 
### assign sample to opposite category. If sample fits in neither category, assign to heterozygote list.

for(sample in 1:ncol(Filtered_candidate_coverage)){
  
  if( sum(Type_1_table[,sample] == 0) > nrow(Type_1_table) / 2){
    Type_2_homo_list <- c(Type_2_homo_list, colnames(Filtered_candidate_coverage)[sample])
  }
  else if( sum(Type_2_table[,sample] == 0) > nrow(Type_2_table) / 2){
    Type_1_homo_list <- c(Type_1_homo_list, colnames(Filtered_candidate_coverage)[sample])
  }
  else{ Hetero_list <- c(Hetero_list, colnames(Filtered_candidate_coverage)[sample]) }
}

Hetero_samples <- Filtered_candidate_coverage[,Hetero_list]
Type_1_homo_samples <- Filtered_candidate_coverage[,Type_1_homo_list]
Type_2_homo_samples <- Filtered_candidate_coverage[,Type_2_homo_list]

### Make outputs

write(Type_1_homo_list, paste(Output, "Type_1_homo_samples", sep = "_"))
write(Type_2_homo_list, paste(Output, "Type_2_homo_samples", sep = "_"))
write(Hetero_list, paste(Output, "Arrested_hetrozygote_samples", sep = "_"))

write(Type_1_markers, paste(Output, "Type_1_markers", sep = "_"))
write(Type_2_markers, paste(Output, "Type_2_markers", sep = "_"))
write(Not_type_2_markers, paste(Output, "Anomolous_markers", sep = "_"))