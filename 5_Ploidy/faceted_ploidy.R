library(ggplot2)
library(ggtext)

# Make normalised data

table_SNP_ratios <- function(table_file,gene_type,sample_type) {
  
  SNP_table <- read.table(table_file)
  SNP_table_norm <- SNP_table
  sample_means <- colMeans(SNP_table[,2:11])
  
  for (i in 2:11){
    SNP_table_norm[,i] <- SNP_table[,i] / sample_means[i-1]
  }
  
  SNP_table_norm$average <- rowMeans(SNP_table_norm[,2:11])
  SNP_table_norm[is.na(SNP_table_norm)] <- 0
  SNP_table_norm$gene_class <- gene_type
  SNP_table_norm$sample_class <- sample_type
  return(SNP_table_norm[,c(1,12,13,14)])
}

Master_table <- table_SNP_ratios("A28_ST", "A_linked", "AA")
Master_table <- rbind(Master_table, table_SNP_ratios("A28_hatch", "A_linked", "AB"))
Master_table <- rbind(Master_table, table_SNP_ratios("A28_FT", "A_linked", "BB"))
Master_table <- rbind(Master_table, table_SNP_ratios("B35_ST", "B_linked", "AA"))
Master_table <- rbind(Master_table, table_SNP_ratios("B35_hatch", "B_linked", "AB"))
Master_table <- rbind(Master_table, table_SNP_ratios("B35_FT", "B_linked", "BB"))
Master_table <- rbind(Master_table, table_SNP_ratios("ALL_ST", "All", "AA"))
Master_table <- rbind(Master_table, table_SNP_ratios("ALL_hatch", "All", "AB"))
Master_table <- rbind(Master_table, table_SNP_ratios("ALL_FT", "All", "BB"))

Master_table$gene_class_f <- factor(Master_table$gene_class, levels = c("All","A_linked","B_linked"))
Master_table$sample_class_f <- factor(Master_table$sample_class, levels = c("AB","AA","BB"))

### load nQuire estimates:

nQuire_table <- read.table("Ploidy_ave.txt", header = TRUE)

nQuire_table$ave_cov <- c(84.5, 42.6, 0.1, 58.6, 58.6, 58.6, 0.1, 32.1, 69.7)
nQuire_table$ave_SNPs <- c(123, 166.6, 0.0, 10756, 12484, 10574, 0.5, 168.4, 104.0)

Ploidy_labels <- c("unkown","Diploid","nothing","diploid", "<b>diploid: 0.75,</b> <br> triploid: 0.03 <br> tetraploid: 0.05","diploid", "nothing", "diploid", "unkown")

Annotations <- data.frame(text = Ploidy_labels)

Annotations$gene_class <- c("All","All","All","A_linked","A_linked","A_linked","B_linked","B_linked","B_linked")
Annotations$sample_class <- c("AA","AB","BB","AA","AB","BB","AA","AB","BB")
Annotations$average <- 3
Annotations$V1 <- 20


for(i in 1:9){
  score_table <- nQuire_table[nQuire_table$gene_class == Annotations$gene_class[i] & nQuire_table$sample_class == Annotations$sample_class[i] ,]
  
  if(is.na(score_table$Diploid[1])){
    score_label <- "**insufficient data**"
  }
  
  else if(score_table$Diploid[1] > 2 * score_table$Triploid[1] & score_table$Diploid[1] > 2 * score_table$Tetraploid[1]){
    score_label <- paste("<b>Diploid:",score_table$Diploid[1],"</b> <br>Triploid:", score_table$Triploid[1], "<br> Tetraploid:", score_table$Tetraploid[1])  
  }
  else{
    score_label <- paste("Diploid:",score_table$Diploid[1]," <br>Triploid:", score_table$Triploid[1], "<br> Tetraploid:", score_table$Tetraploid[1])
  }
  Annotations$text[i] <- score_label
  
  Annotations$text_2[i] <- paste("Coverage: ", score_table$ave_cov[1], " <br>SNPs:", score_table$ave_SNPs[1])
}



Annotations$gene_class_f <- factor(Annotations$gene_class, levels = c("All","A_linked","B_linked"))
Annotations$sample_class_f <- factor(Annotations$sample_class, levels = c("AB","AA","BB"))

gene_labels <- c('All' = "All genes", 'A_linked' = "A-linked genes", 'B_linked' = "B-linked genes", "AA" ="AA", "AB" = "AB", "BB" = "BB")

ggplot(Master_table, aes(x = V1, y = average)) + geom_area() + ylim(0,4) +
  facet_grid(sample_class_f ~ gene_class_f, switch = 'y', labeller = as_labeller(gene_labels)) +
  geom_richtext(data = Annotations, label = Annotations$text, family = "serif", hjust = "left", fill = NA, label.color = NA, size = 3) +
  geom_richtext(data = Annotations, label = Annotations$text_2, family = "serif", hjust = "right", x = 80, fill = NA, label.color = NA, size = 3) +
  xlab("% allele") + ylab("proportion of SNPs in genotype") + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.y.left = element_text(angle = 0, size = 10, face = 'bold'),
        strip.text.x = element_text(size = 10, face = 'bold'), 
        strip.background = element_blank())
  

