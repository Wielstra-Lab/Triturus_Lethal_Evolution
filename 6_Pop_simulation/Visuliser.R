### Visuliser ###

# Produces a image of the landscape simulate in BL sim

args <- commandArgs(trailingOnly = TRUE)

### get arguments ###

Save_file <- args[1]
Generation <- args[2]
Output_file <- args[3]

require(grid)
require(gridExtra)
require(ggplot2)

visulise_landscape <- function(gen_num) {
  
  Landscape <- state_list$landscape
  
  pops_total <- as.numeric(state_list$pop_history[gen_num,])
  pops_BL <- as.numeric(state_list$BL_history[gen_num,])
  pops_NN <- as.numeric(state_list$NN_history[gen_num,])
  
  Landscape$colour <- ifelse(pops_total > 0, 'purple', 'white')
  Landscape$colour <- ifelse(pops_BL > 0.8 * pops_total, 'red', Landscape$colour)
  Landscape$colour <- ifelse(pops_NN > 0.8 * pops_total, 'blue', Landscape$colour)
  
  pond_circles <- gTree()
  
  for(i in 1:nrow(Landscape)){
    
    pond_dot <- circleGrob(x = Landscape[i,"grid_x"], y = Landscape[i,"grid_y"],
                           r = Landscape[i,"size_rt"], gp = gpar(fill = Landscape[i,"colour"]))
    pond_circles <- addGrob(pond_circles, pond_dot)
  }
  
  gen_label <- textGrob(label = paste("generation:", gen_num), x = 0.1, y = 0.1, just = 'left', gp=gpar(size = 3))
  pond_circles <- addGrob(pond_circles, gen_label)

}

load(Save_file)
grid.draw(visulise_landscape(Generation))
