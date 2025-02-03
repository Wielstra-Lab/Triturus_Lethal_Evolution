### Visuliser ###

args <- commandArgs(trailingOnly = TRUE)

require(grid)
require(gridExtra)
require(ggplot2)

### get arguments ###

if(length(args) == 0) {
  print("no arguments supplied")
  input_file <- args[1]
}

if(length(args) == 1) {
  input_file <- args[1]
  output_file <- "end"
  num_maps <- 1
  make_graph <- FALSE
  vis_generations <- "end"
  
} else if(length(args) == 2) {
    input_file <- args[1]
    output_file <- args[2]
    num_maps <- 1
    make_graph <- FALSE
    vis_generations <- "end"
    
} else if(length(args) == 3) {
    input_file <- args[1]
    output_file <- args[2]
    num_maps <- args[3]
    make_graph <- FALSE
    vis_generations <- "end"
    
} else if(length(args) == 4) {
    input_file <- args[1]
    output_file <- args[2]
    num_maps <- args[3]
    make_graph <- args [4]
    vis_generations <- "end"
    
} else {
    input_file <- args[1]
    output_file <- args[2]
    num_maps <- args[3]
    make_graph <- args[4]
    
    vis_generations <- args[5:length(args)]
    
    if(length(vis_generations) != num_maps & length(vis_generations) > 1) {
      num_maps <- length(vis_generations)
      print("adjusting number of maps to match generations requested")
      }
}

pick_generations <- function(num_maps, final_generation) {
  
  print(paste("final_gen:", final_generation))
  
  print(paste("numbers of maps:", num_maps ))
  
  ## picks generations to display if none specified 
  
  if(as.integer(num_maps) == 1) {return(final_generation)}
  
  root_n <- 3
  final_root <- as.integer(final_generation) ^ (1/ root_n)
  gen_roots <- 0:(as.integer(num_maps) - 1) * (final_root / (as.integer(num_maps) - 1))
  raw_gens <- gen_roots ^ root_n
  refined_gens <- as.integer(signif(raw_gens, digits = 1))
  return(refined_gens)
}

visulise_landscape <- function(gen_num) {
  
  Landscape <- state_list$landscape
  
  pops_total <- as.numeric(state_list$pop_history[gen_num + 1,])
  pops_BL <- as.numeric(state_list$BL_history[gen_num + 1,])
  pops_NN <- as.numeric(state_list$NN_history[gen_num + 1,])
  
  Landscape$colour <- ifelse(pops_total > 0, '#b845ff', 'white')
  Landscape$colour <- ifelse(pops_BL > 0.8 * pops_total, '#a2f76e', Landscape$colour)
  Landscape$colour <- ifelse(pops_NN > 0.8 * pops_total, 'black', Landscape$colour)
  
  pond_circles <- gTree()
  
  for(i in 1:nrow(Landscape)){
    
    pond_dot <- circleGrob(x = Landscape[i,"grid_x"] * 0.95, y = Landscape[i,"grid_y"] * 0.95 + 0.05,
                           r = Landscape[i,"size_rt"] * pond_scale, gp = gpar(fill = Landscape[i,"colour"]))
    pond_circles <- addGrob(pond_circles, pond_dot)
  }
  
  
  #label_border <- rectGrob(x = 0.09, y = 0.1, just = "left", width = 0.24, height = 0.08)
  gen_label <- textGrob(label = paste("generation:", gen_num), x = 0.1, y = 0.05, just = 'left', gp=gpar(size = 3))
  #pond_circles <- addGrob(pond_circles, label_border)
  pond_circles <- addGrob(pond_circles, gen_label)
  
  return(gList(pond_circles))
  
  #grid.newpage()
  #grid.draw(pond_circles)
  #grid.rect(x = 0.09, y = 0.1, just = "left", width = 0.24, height = 0.08)
  #grid.text(label = paste("generation:", gen_num), x = 0.1, y = 0.1, just = 'left')
}

pop_graph <- function(max_gen = 100) {
  
  gen <- 0:max_gen
  sum_pops_total <- rowSums(state_list$pop_history[0:max_gen + 1,])
  sum_pops_BL <- rowSums(state_list$BL_history[0:max_gen + 1,])
  sum_pops_NN <- rowSums(state_list$NN_history[0:max_gen + 1,])
  sum_pops_other <- sum_pops_total - (sum_pops_BL + sum_pops_NN)
  
  # print(sum_pops_total)
  # print(length(sum_pops_other))
  
  pop_table <- data.frame(gen,sum_pops_total,sum_pops_BL,sum_pops_NN,sum_pops_other)
  
  pop_chart <- ggplot(data = pop_table, aes(x = gen, y = sum_pops_total)) + 
    #geom_line(color = "grey", linewidth = 1, linetype = 'dotted') + 
    geom_line(y = pop_table$sum_pops_NN, color = "black", linewidth = 1) +
    geom_line(y = pop_table$sum_pops_BL, color = "#a2f76e", linewidth = 1) +
    geom_line(y = pop_table$sum_pops_other, color = "#b845ff", linewidth = 1) +
    labs(x = "Generation", y = "Population") +
    scale_x_continuous(expand = c(0,0), limits = c(0,final_generation)) +
    scale_y_continuous(expand = c(0.005,0), limits = c(0,max_pop_geno * 1.05)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
          panel.border =  element_rect(colour = "black", fill=NA))
  
  return(gList(ggplotGrob(pop_chart)))
  
  #ggplot(data = pop_sums_table, aes(x = gen, y = sum_pops_total)) + geom_line(x = gen, y = sum_pops_total, color = "red", linewidth = 1)
  
}

load(input_file)
final_generation <- nrow(state_list$pop_history) - 1

if (output_file == "end") {output_file <- paste(input_file, "_gen_", final_generation, ".pdf", sep = "")}
if (vis_generations[1] == "end") {vis_generations <- pick_generations(num_maps, final_generation)}
if (length(vis_generations) == 1) {vis_generations <- pick_generations(num_maps, vis_generations)}

if (make_graph == TRUE) {total_figs <- as.integer(num_maps) + 1} else {total_figs <- as.integer(num_maps)}

pond_scale <- sqrt(5000 / max(state_list$landscape$size))

max_NN <- max(rowSums(state_list$NN_history))
max_BL <- max(rowSums(state_list$BL_history))
max_pop_geno <- max(c(max_NN, max_BL))

print(paste("input_file:", input_file ))
print(paste("output_file:", output_file ))
print(paste("numbers of maps:", num_maps ))
print(paste("draw graph:", make_graph ))
print(paste("vis_generations:", vis_generations ))
print(paste("total_figs:", total_figs ))

print(paste("pond_scale:", pond_scale ))

fig_list <- list()

if(num_maps > 0) {
  for (fig in 1:length(vis_generations)) {
    map_plot <- visulise_landscape(vis_generations[fig])
    fig_list <- c(fig_list, map_plot )
  }
}

if (make_graph == TRUE) {
  pop_plot <- pop_graph(max_gen = final_generation)
  fig_list <- c(fig_list, pop_plot)
}

if(total_figs == 1) {
  pdf_col = 1
  pdf_width = 7
  pdf_height = 7
} else {
  pdf_col = 2
  pdf_width = 8.25 
  pdf_height = (total_figs %/% 2) * (11.75 / 3)
}

pdf(output_file, width = pdf_width, height = pdf_height)
grid.arrange(grobs = fig_list, pdf_col = 2)
dev.off()

