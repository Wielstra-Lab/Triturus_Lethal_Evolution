### Balanced Lethal System Evolution Simulation ###

# Version 1.0
#
# Script models the evolution of a population (of Triturus newts) with on locus with 3 possible alleles
# and so 9 possible genotypes - where the original state is NN and alleles A and B are derived.
# This is indented to simulate a balanced lethal system so genotypes AA and BB have 0 survival.
#
# Normally a single individual in an NN population will be mutated to AB, however this script can also
# simulate secondary contact between two populations, one with genotype NN fixed and the other fixed AB
#
# The simulation area consists of a landscape of between 50 to 150 sites (ponds) each with its own 
# sub-population, where juvinial individuals can disperse between them. Ponds are located randomly within
# a 5 x 5 km grid, and max dispersal distance is 1 km (but chance of dispersal falls off linearly with
# distance).
#
# Input parameters are the number of generations of each iteration of the simulation, the number of
# iterations, the fitness values for the various genotypes, and whether secondary contact is modeled. 
#
# Each simulation has its own unique name and the random seed used can be set manually if desired. 
# Overall statistics are recorded and the "successful" (where the A and B alleles survive) iterations 
# are saved. This results in a save file that can be visualized with the accompanying script.

args <- commandArgs(trailingOnly = TRUE)

### get arguments ###

Sim_name <- args[1]
number_of_iterations <- args[2]
number_of_generations <- args[3]
genotype_file <- args[4]
input_seed <- args[5]
secondary_contact <- args[6]

### set seed parameters ###

if(input_seed == "random"){
  set.seed(NULL)
  base_seed <- sample(1:10000000,1)
}else{
  base_seed <- input_seed
}

### set simulation area parameters ###

area_x_length <- 5000
area_y_length <- 5000
min_ponds <- 50
max_ponds <- 150

start_from_centre <- FALSE

juvinial_cull_ratio <- 20 # The reciprocal of this is the proportion of hatched eggs which survive year 0

### load or generate genotypes ###

genotypes <- c("AA","AB","AN","BA","BB","BN","NA","NB","NN")

embryo_fitness <- c(0, 100, 75, 100, 0, 75, 75, 75, 100)
juvinial_survival <- c(0, 20, 15, 20, 0, 15, 15, 15, 20)
adult_survival <- c(0, 80, 60, 80, 0, 60, 60, 60, 80)
fecundicity <- c(0, 200, 150, 200, 0, 150, 150, 150, 200)
attraction <- c(0, 100, 75, 100, 0, 75, 75, 75, 100)

if(genotype_file == "default"){
  genotype_table <- data.frame(genotypes, embryo_fitness, juvinial_survival, adult_survival, fecundicity, attraction)
}else{ 
  genotype_table <- read.table(genotype_file, na.strings = "'NA'")
}

### make migrant template and pond class ###

master_migrant_template <- data.frame(sex = character(), age = integer(), genotype = character(), origin = integer(), 
                                      destination = integer(), disperse = logical())

setClass("pond", slots = list(pond_ID = 'numeric', size = 'numeric', routes = 'numeric', col_gen = 'numeric',
                              virgin = 'logical', route_rating = 'numeric',  pond_table = 'data.frame'))

### simulation  sub-functions ###

breed_clutch <- function(father, mother) {
  
  match(mother, genotypes)
  
  clutch_size <- genotype_table$fecundicity[match(mother, genotypes)]
  
  father_genotype_1 <- substring(father, 1, 1)
  father_genotype_2 <- substring(father, 2)
  
  mother_genotype_1 <- substring(mother, 1, 1)
  mother_genotype_2 <- substring(mother, 2)
  
  mix_A <- paste(father_genotype_1, mother_genotype_1, sep = "")
  mix_B <- paste(father_genotype_1, mother_genotype_2, sep = "")
  mix_C <- paste(father_genotype_2, mother_genotype_1, sep = "")
  mix_D <- paste(father_genotype_2, mother_genotype_2, sep = "")
  
  clutch_punnet <- c(mix_A, mix_B, mix_C, mix_D)
  
  punnet_codes <- match(clutch_punnet,genotypes)
  
  clutch_geno_codes <- sample(punnet_codes, clutch_size, replace = TRUE)
  clutch_fitness <- embryo_fitness[clutch_geno_codes]
  
  clutch_table <- data.frame(clutch_geno_codes, clutch_fitness)

}

breed_pond <- function(pond){
  
  newt_table <- pond@pond_table
  
  female_table <- newt_table[which(newt_table$sex == 'female' & newt_table$age > 1),]
  male_table <- newt_table[which(newt_table$sex == 'male' & newt_table$age > 1),]
  
  if(nrow(male_table) == 0 || nrow(female_table) == 0) {return(newt_table)}
  
  num_AB <- sum(newt_table$genotype == 'BA' | newt_table$genotype == 'AB')
  num_NN <- sum(newt_table$genotype == 'NN')
  
  if(num_NN == nrow(newt_table)) {

    num_offspring <- 200 * nrow(female_table)
    
    clutch_geno_codes <- rep(9, num_offspring)
    clutch_fitness <- rep(100, num_offspring)
    
    survivors_table <- data.frame(clutch_geno_codes, clutch_fitness)
  }
  
  else if (num_AB == nrow(newt_table)) {
    
    num_offspring <- 200 * nrow(female_table)
    
    clutch_geno_codes <- sample(c(2,4), num_offspring, replace = TRUE)
    clutch_fitness <- sample(c(0,100), num_offspring, replace = TRUE)
    
    offspring_table <- data.frame(clutch_geno_codes, clutch_fitness)
    survivors_table <- offspring_table[which(offspring_table$clutch_fitness == 100),]
  }
  
  else{

    male_table$attraction <- genotype_table$attraction[match(male_table$genotype,genotypes)]
    
    female_table$mate <- sample(1:nrow(male_table),nrow(female_table),replace = TRUE, prob = male_table$attraction)
    
    offspring_table <- data.frame(geno_codes = integer(), fitness = numeric())
    
    for(i in 1:nrow(female_table)) {
      father_id <- female_table[i,"mate"]
      new_clutch_table <- breed_clutch(female_table[i,"genotype"],male_table[father_id,"genotype"])
      offspring_table <- rbind(offspring_table,new_clutch_table)
    }
    
    offspring_table$luck <- sample(0:99, nrow(offspring_table), replace = TRUE)
    offspring_table$surival <- offspring_table$luck + offspring_table$clutch_fitness
    
    survivors_table <- offspring_table[which(offspring_table$surival >= 100),]
  }
  
  survivors_table$luck <- sample(0:99, nrow(survivors_table), replace = TRUE)
  
  carrying_threshold <- 100 - (pond@size * 100 / nrow(survivors_table))
  
  survivors_table_2 <- survivors_table[which(survivors_table$luck >= carrying_threshold),]
  
  juvinial_num <- as.integer(nrow(survivors_table_2) / juvinial_cull_ratio)
  
  survivors_table_3 <- survivors_table_2[sample(1:nrow(survivors_table_2),juvinial_num, replace = FALSE),]
  
  survivors_table_3$sex <- sample(c('male','female'),nrow(survivors_table_3), replace = TRUE)
  survivors_table_3$age <- 0
  survivors_table_3$genotype <- genotype_table[survivors_table_3$clutch_geno_codes,'genotypes']
  survivors_table_3$origin <- pond@pond_ID

  return(rbind (newt_table, survivors_table_3[,c('sex','age','genotype','origin')]))
}

age_pond <- function(pond){
  
  pond_newts <- pond@pond_table
  
  pond_newts$juvinial_survival <- genotype_table$juvinial_survival[match(pond_newts$genotype,genotypes)]
  pond_newts$adult_survival <- genotype_table$adult_survival[match(pond_newts$genotype,genotypes)]
  
  # survival chances for newts of age 0 are 100% at this stage since 
  # egg to juvinial mortality already accounted for in breed_pond
  
  pond_newts$survival <- ifelse(pond_newts$age > 1, pond_newts$adult_survival, 
                                ifelse(pond_newts$age == 1, pond_newts$juvinial_survival, 100))
  
  pond_newts$luck <- sample(1:100, nrow(pond_newts), replace = TRUE)
  pond_newts <- pond_newts[which(pond_newts$survival >= pond_newts$luck),]
  
  pond_newts$age <- pond_newts$age + 1
  pond_newts <- pond_newts[,c('sex','age','genotype','origin')]
  return(pond_newts)
}

disperse_newts <- function(pond){
  
  pond_newts <- pond@pond_table
  
  if(length(pond@routes) == 1) {return(pond_newts)}
  
  pond_newts$disperse <- sample(c(TRUE,FALSE), nrow(pond_newts), replace = TRUE)
  migrants <- pond_newts[which(pond_newts$disperse == TRUE & pond_newts$age == 1),]
  
  migrants$destination <- sample(pond@routes, nrow(migrants), replace = TRUE, prob = pond@route_rating)
  master_migrant_table <<- rbind(master_migrant_table, migrants)
  
  pond_newts <- pond_newts[which(pond_newts$disperse == FALSE | pond_newts$age != 1),c('sex','age','genotype','origin')]
  return(pond_newts)
}

collect_newts <- function(pond,generation){
  
  incoming_migrants <- master_migrant_table[which(master_migrant_table$destination == pond@pond_ID),]
  
    if(nrow(incoming_migrants) > 0 & pond@virgin == TRUE) {
      pond@virgin <- FALSE
      pond@col_gen <- generation
      }
  
  pond@pond_table <- (rbind(pond@pond_table, incoming_migrants[,c('sex','age','genotype','origin')]))
  return(pond)
}

### landscape generation functions ### 

generate_landscape <- function(num_ponds = 100, area_x = 5000, area_y = 5000){
  
  pond_ID <- c(1:num_ponds)
  
  table_of_ponds <- data.frame(pond_ID)
  table_of_ponds$x_coord <- round(runif(n = num_ponds, min = 0, max = area_x))
  table_of_ponds$y_coord <- round(runif(n = num_ponds, min = 0, max = area_y))
  table_of_ponds$size <- sample(c(500,750,1000,1500,2000,2500,3000,4000,5000), num_ponds, replace = TRUE, prob = c(9:1))
  
  table_of_ponds$grid_x <- (table_of_ponds$x_coord / (area_x * 1.2)) + 0.1
  table_of_ponds$grid_y <- (table_of_ponds$y_coord / (area_x * 1.2)) + 0.1
  table_of_ponds$size_rt <- (sqrt(table_of_ponds$size) / 2000)
  
  table_of_ponds$total_pop <- 0
  table_of_ponds$num_NN <- 0
  table_of_ponds$num_BL <- 0
  table_of_ponds$num_other <- 0
  
  return(table_of_ponds)
}

generate_ponds <- function(table_of_ponds) {
  
  num_ponds <- nrow(table_of_ponds)
  
  all_ponds <- list()
  
  empty_pond <- data.frame(sex = character(), age = numeric(), genotype = character(), origin = numeric())
  
  distance_matrix <- matrix(nrow = num_ponds, ncol = num_ponds)
  table_of_ponds$num_neighbours <- NA
  
  for( i in 1:num_ponds) {
    x_distances <- table_of_ponds$x_coord - table_of_ponds$x_coord[i]
    y_distances <- table_of_ponds$y_coord - table_of_ponds$y_coord[i]
    distances <- sqrt(x_distances^2 + y_distances^2)
    distance_matrix[,i] <- round(distances)
    table_of_ponds$num_neighbours[i] <- sum(distances < 1000) - 1
  }
  
  for( i in 1:num_ponds) {
    pond_neigbours <- which(distance_matrix[,i] < 1000)
    pond_neigbour_distances <- distance_matrix[pond_neigbours,i]
    pond_neighbour_sizes <- table_of_ponds$size[pond_neigbours]
    pond_neighbour_weights <- round((1000 - pond_neigbour_distances) * pond_neighbour_sizes / 1000)
    
    pond_X <- new('pond', pond_ID = i, routes = pond_neigbours, route_rating = pond_neighbour_weights, 
                  col_gen = 0, virgin = TRUE, pond_table = empty_pond, size = table_of_ponds$size[i])
    all_ponds <- c(all_ponds, pond_X)
  }
  
  return(all_ponds)
}

### core simulation function ###

run_simulation <- function(landscape, pond_list, generations = 10, start_gen = 0, save_run = TRUE, 
                            run_num = 0, run_name = "Run", save_directory = "Sim_saves"){
  
  gen_counter <- start_gen
  all_ponds <- pond_list
  BL_exist <- FALSE
  BL_start <- 5
  
  BL_start_gen <- 0
  BL_start_total_pop <- 0
  BL_start_pond <- 0
  BL_start_pond_pop <- 0
  BL_start_colonies <- 0
  Peak_BL <- 0
  
  finish_condition <- 0
  
  for(pond_num in 1:length(all_ponds)){
    if( nrow(all_ponds[[pond_num]]@pond_table) > 0 ){
      
      all_ponds[[pond_num]]@col_gen <- gen_counter
      all_ponds[[pond_num]]@virgin <- FALSE
    }
  }
  
  ### intialise history ###
  
  pond_pop_hist <- data.frame(matrix(nrow = 0, ncol = length(all_ponds)))
  
  pond_BL_hist <- pond_pop_hist
  pond_NN_hist <- pond_pop_hist
  
  pond_pop_vec <- c()
  pond_BL_vec <- c()
  pond_NN_vec <- c()
  
  for(pond_num in 1:length(all_ponds)){
    
    pond_pop <- nrow(all_ponds[[pond_num]]@pond_table)
    pond_BL <- sum(all_ponds[[pond_num]]@pond_table$genotype == 'AB') + sum(all_ponds[[pond_num]]@pond_table$genotype == 'BA')
    pond_NN <- sum(all_ponds[[pond_num]]@pond_table$genotype == 'NN')
    
    pond_pop_vec <- c(pond_pop_vec, pond_pop)
    pond_BL_vec <- c(pond_BL_vec, pond_BL)
    pond_NN_vec <- c(pond_NN_vec, pond_NN)
  }
  
  pond_pop_hist <- rbind(pond_pop_hist, pond_pop_vec)
  pond_BL_hist <- rbind(pond_BL_hist, pond_BL_vec)
  pond_NN_hist <- rbind(pond_NN_hist, pond_NN_vec)
  
  colnames(pond_pop_hist) <- 1:length(all_ponds)
  colnames(pond_BL_hist) <- 1:length(all_ponds)
  colnames(pond_NN_hist) <- 1:length(all_ponds)
  
  if(sum(pond_BL_vec) > 0 ) {BL_exist <- TRUE}
  

  for(i in 1:generations) {
    
    ### Iterate landscape ###
    
    gen_counter <- gen_counter + 1
    
    master_migrant_table <<- master_migrant_template
    
    for(pond_num in 1:length(all_ponds)){
      
      all_ponds[[pond_num]]@pond_table <- breed_pond(all_ponds[[pond_num]])
      all_ponds[[pond_num]]@pond_table <- disperse_newts(all_ponds[[pond_num]])
    }
    
    for(pond_num in 1:length(all_ponds)){
      
      all_ponds[[pond_num]] <- collect_newts(all_ponds[[pond_num]], gen_counter)
      all_ponds[[pond_num]]@pond_table <- age_pond(all_ponds[[pond_num]])
    }
    
    
    ### count stats ###
    
    BL_total <- 0
    NN_total <- 0
    Pop_total <- 0
    Other_total <- 0
    
    pond_pop_vec <- c()
    pond_BL_vec <- c()
    pond_NN_vec <- c()
    
    for(pond_num in 1:length(all_ponds)){
      
      pond_pop <- nrow(all_ponds[[pond_num]]@pond_table)
      pond_BL <- sum(all_ponds[[pond_num]]@pond_table$genotype == 'AB') + sum(all_ponds[[pond_num]]@pond_table$genotype == 'BA')
      pond_NN <- sum(all_ponds[[pond_num]]@pond_table$genotype == 'NN')
      
      pond_pop_vec <- c(pond_pop_vec, pond_pop)
      pond_BL_vec <- c(pond_BL_vec, pond_BL)
      pond_NN_vec <- c(pond_NN_vec, pond_NN)
      
      BL_total <- BL_total + pond_BL
      NN_total <- NN_total + pond_NN 
      Pop_total <- Pop_total + pond_pop
    }
    
    pond_pop_hist <- rbind(pond_pop_hist, pond_pop_vec)
    pond_BL_hist <- rbind(pond_BL_hist, pond_BL_vec)
    pond_NN_hist <- rbind(pond_NN_hist, pond_NN_vec)
    
    num_colonies <- sum(pond_pop_vec > 0)
    
    if(BL_total > Peak_BL){ Peak_BL <- BL_total}
    
    Other_total <- Pop_total - (BL_total + NN_total)
    
    ### Initiate Balanced lethal system ###
    
    if( gen_counter >= BL_start & BL_exist == FALSE) {
      
      virgin_vec <- c()
      col_gen_vec <- c()
      for(pond_num in 1:length(all_ponds)){ 
        col_gen_vec <- c(col_gen_vec, all_ponds[[pond_num]]@col_gen)
      }
      
      if(gen_counter %in% col_gen_vec){
        new_colonies <- which(col_gen_vec == gen_counter)
        BL_start_pond <- sample(new_colonies,1)
        
        if(nrow(all_ponds[[BL_start_pond]]@pond_table) > 0 & col_gen_vec[BL_start_pond] == gen_counter){
          all_ponds[[BL_start_pond ]]@pond_table$genotype[1] <- 'AB'
          print(paste("BL begins at pond:", BL_start_pond, "at gen:", gen_counter))
          
          BL_exist <- TRUE
          BL_start_gen <- gen_counter
          BL_start_total_pop <- Pop_total
          BL_start_pond_pop <- nrow(all_ponds[[BL_start_pond]]@pond_table)
          
          BL_start_colonies <- num_colonies
        }
      }
    }  
    
    #print(paste("Gen:", gen_counter ,"Total pop:", Pop_total, "Total BL, NN, mixed:", BL_total, NN_total, Other_total))
    
    if(Pop_total == 0){ 
      finish_message <- "Total extinction"
      finish_condition <- 1 
    }
    
    if(Pop_total == BL_total){ 
      finish_message <- "Balanced lethal system fixed"
      finish_condition <- 2
    }
    
    if(Pop_total == NN_total & BL_exist == TRUE & gen_counter > BL_start_gen){ 
      finish_message <- paste("Ancestral chromosome fixed at gen:", gen_counter, "Peak_BL:", Peak_BL)
      finish_condition <- 3
    }
    
    if(gen_counter == 25 & BL_exist == FALSE){ 
      finish_message <- paste("Failed to spawn balanced lethal system")
      finish_condition <- 4
    }
    
    if(finish_condition > 0){
      print(finish_message)
      return(c(finish_condition, gen_counter,  BL_start_gen, BL_start_total_pop, BL_start_pond_pop, BL_start_colonies, num_colonies, Peak_BL))
    }
  }
  
  finish_message <- "still going!"
  print(finish_message)
  
  ### save state ###
  
  save_name <- paste(save_directory, "/", run_name, "_run_", run_num, sep = "")
  
  run_info <- list(base_seed, iteration_seed, run_name, run_num, start_from_centre, juvinial_cull_ratio,
                   BL_start_gen, BL_start_total_pop, BL_start_pond, BL_start_pond_pop)
  names(run_info) <- c("base_seed", "iteration_seed", "run_name", "run_num", "start_from_centre", "juvinial_cull_ratio",
                       "BL_start_gen", "BL_start_total_pop", "BL_start_pond", "BL_start_pond_pop")
  
  state_list <<- list(run_info, genotype_table, landscape, pond_pop_hist, pond_BL_hist, pond_NN_hist, all_ponds)
  names(state_list) <<- c("run_info", "genotype_info", "landscape", "pop_history", "BL_history", "NN_history", "pond_data")
  
  if(save_run == TRUE) {save(state_list, file = save_name)}
  
  print(finish_message)
  return(c(finish_condition, gen_counter,  BL_start_gen, BL_start_total_pop, BL_start_pond_pop, BL_start_colonies, num_colonies, Peak_BL))
}

### main program ###

print(paste("Sim name: ", Sim_name))
print(paste("Number of iterations: ", number_of_iterations))
print(paste("Base seed: ", base_seed))
print(paste("start from center:", start_from_centre))
print(paste("genotypes:", genotype_file))
print("genotype_information:")
print(genotype_table)

dir.create(Sim_name)

outcomes_file <- paste(Sim_name, "_outcomes.txt", sep = "")

outcomes <- data.frame(final_state = integer(), final_gen = integer(), BL_start = integer(),
                       start_pop_tot = integer(), start_pond_pop = integer(), colonies_at_BL_start = integer(),
                       colonies_at_end = integer(), BL_peak = integer(), total_ponds = integer())

### simulation loop ###

for(iteration in 1:number_of_iterations){
  
  iteration_seed <- (base_seed + iteration)
  set.seed(iteration_seed)
  
  print(iteration)
  
  iteration_pond_num <- sample(min_ponds:max_ponds, 1)
  iteration_landscape <- generate_landscape(num_ponds = iteration_pond_num, area_x = area_x_length, area_y = area_y_length)
  
  centre_pond <- which.min(sqrt((iteration_landscape$x_coord - area_x_length / 2) ^ 2 + (iteration_landscape$y_coord - area_y_length / 2)) ^ 2)
  westmost_pond <- which.min(iteration_landscape$x_coord)
  eastmost_pond <- which.max(iteration_landscape$x_coord)
  
  if(start_from_centre == TRUE) {start_pond <- centre_pond}
  else {start_pond <- westmost_pond}

  start_sex <- rep(c("male","female"), each=10)
  start_age <- rep(2, each=20)
  start_genotype <- rep("NN", each=20)
  start_origin <- rep(start_pond, each=20)
  newt_table_start <- data.frame(sex=start_sex, age=start_age, genotype=start_genotype, origin=start_origin)
  newt_table_start_AB <- newt_table_start
  newt_table_start_AB$genotype <- "AB"
  
  iteration_pond_list <- generate_ponds(iteration_landscape)
  iteration_pond_list[[start_pond]]@pond_table <- newt_table_start
  iteration_pond_list[[start_pond]]@virgin <- FALSE
  
  if(secondary_contact == TRUE) {
  
    iteration_pond_list[[eastmost_pond]]@pond_table <- newt_table_start_AB
    iteration_pond_list[[eastmost_pond]]@virgin <- FALSE
  }
  
  stats_out <- run_simulation(iteration_landscape, iteration_pond_list, number_of_generations, save_directory = Sim_name, run_name = Sim_name, run_num = iteration)
  outcomes[iteration,] <- c(stats_out, iteration_pond_num)
}

write.table(outcomes, outcomes_file, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE, na = "NA")
