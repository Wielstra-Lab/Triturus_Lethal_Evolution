### Script designed for sanity checking presence absence markers

# outputs stats for all markers and a list of markers screened for: absence in at least 20 samples,
# median coverage of at least 5, and mean coverage of at least 5 

args_in <- commandArgs(trailingOnly = TRUE)

### First load coverage table and list of Presence/Absence markers

input_file <- args_in[1]
PA_marker_file <- args_in[2]
output_file <- args_in[3]

output_stats <- paste(output_file, "_stats.txt", sep = "")
output_list <- paste(output_file, "_screened.txt", sep = "")

Coverage_table_full <- read.table(input_file,stringsAsFactors = F)

PA_marker_table <- read.table(PA_marker_file)
PA_marker_List <- PA_marker_table$V1

#Subset to only markers of interest, and make table of means, medians and zeros

Coverage_table <- Coverage_table_full[PA_marker_List,]

Stats_table <- Coverage_table[,FALSE]

Stats_table$mean_cov <- apply(Coverage_table, 1, mean, na.rm=TRUE)
Stats_table$median_cov <- apply(Coverage_table, 1, median, na.rm=TRUE)
Stats_table$zero_cov <- apply(Coverage_table == 0, 1, sum, na.rm=TRUE)

# Screen markers for absence in at least 20 samples

Screened_table <- subset(Stats_table, zero_cov >= 20 & mean_cov >= 5 & median_cov >= 5)
Screened_markers <- rownames(Screened_table)

# report number of markers excluded for various deficiencies 

print(paste("Total number of raw markers:", nrow(Stats_table)))
print(paste("Total number of screened markers:", nrow(Screened_table)))

print(paste("markers present in too many samples: ", sum(Stats_table$zero_cov < 20)))
print(paste("markers with too low mean coverage: ", sum(Stats_table$mean_cov < 5)))
print(paste("markers with too low median coverage: ", sum(Stats_table$median_cov < 5)))


# write the full stats table and list of screened markers

write.table(Stats_table, output_stats, sep = "\t", quote = FALSE)
write(Screened_markers, output_list)