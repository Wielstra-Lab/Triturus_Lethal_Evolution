### Balanced Lethal System Evolution Simulation ###

# Version 3
#
# Changes from version 2: 
#
# Carrying capactiy now reworked to be based on adult population size, adults now have chance to 
# skip breeding if the pond is over-capacity, so that on average the number of breeding adults will not
# exceed the carrying capactity
#
# breed_pond function reworked to remove loops and so improve performance
#
# Added connect_landscape to make sure a graph of the landscape (with ponds as nodes and viable dispersal routes as edges)
# has a single component, so eventually a population can spread from one pond to any other pond - previously the starting population
# would sometimes begin in an isolated pond or cluster of ponds, and be unable to reach the rest of the landscape
#


args <- commandArgs(trailingOnly = TRUE)

### get arguments ###

Sim_name <- args[1]
number_of_iterations <- args[2]
number_of_generations <- args[3]
genotype_file <- args[4]
input_seed <- args[5]
secondary_contact <- args[6]


### set seed parameters ###

# if a seed is specifed in the input, use it, otherwise pick one randomly
# the seed then determines the pseudo-random numbers used in all parts of the simulation

if(input_seed == "random"){
  set.seed(NULL)
  base_seed <- sample(1:10000000,1)
}else{
  base_seed <- as.numeric(input_seed)
}

### set simulation area parameters ###

area_x_length <- 5000
area_y_length <- 5000
min_ponds <- 50
max_ponds <- 150

start_from_centre <- FALSE
optimise_AB <- TRUE # this enables a shortcut in the breed_ponds function to improve performance, but must be disabled if AA or BB genotypes are of non-zero fitness

juvinial_cull_ratio <- 20 # The reciprocal of this is the proportion of hatched eggs which survive year 0
disperse_percentage <- 50

### load or generate genotypes ###

# the default hybrid fitness is 0.75, if a different value is specified then all AN and BN genotype fitness parameters are scaled by that value
# alternatively a file can be specifed to give custom values for all fitness parameters for all genotypes, individualy

hybrid_fitness <- 0.75

if (!is.na(as.numeric(genotype_file))) { 
  hybrid_fitness <- as.numeric(genotype_file)
  print(paste("hybrid fitness scaled to:", hybrid_fitness))
  }

genotypes <- c("AA","AB","AN","BA","BB","BN","NA","NB","NN")

HF <- hybrid_fitness

genotype_fitness <- c(0, 1, HF, 1, 0, HF, HF, HF, 1)

embryo_fitness <- genotype_fitness * 100
juvinial_survival <- genotype_fitness * 20
adult_survival <- genotype_fitness * 80
fecundicity <- genotype_fitness * 200
attraction <- genotype_fitness * 100

genotype_table <- data.frame(genotypes, embryo_fitness, juvinial_survival, adult_survival, fecundicity, attraction)

# load genotype table file if hybrid_fitness is not default or a numeric value

if(genotype_file != "default" && is.na(as.numeric(genotype_file))){
  genotype_table <- read.table(genotype_file, na.strings = "'NA'")
}

# calculate relative hybrid fitness - this assumes NN has the dafult parameters and all hybrid genotypes have the same parameters

female_hybrid_fitness <- fecundicity[3] * embryo_fitness[3] *  juvinial_survival[3] / (1 -  (adult_survival[3] * 0.01))
male_hybrid_fitness <- attraction[3] * embryo_fitness[3] *  juvinial_survival[3] / (1 -  (adult_survival[3] * 0.01))

female_ancestral_fitness <- fecundicity[9] * embryo_fitness[9] *  juvinial_survival[9] / (1 -  (adult_survival[9] * 0.01))
male_ancestral_fitness <- attraction[9] * embryo_fitness[9] *  juvinial_survival[9] / (1 -  (adult_survival[9] * 0.01))

female_hybrid_relative_fitness <- female_hybrid_fitness / female_ancestral_fitness
male_hybrid_relative_fitness <- male_hybrid_fitness / male_ancestral_fitness

### make migrant template and pond class ###

# master_migrant_table will act as a clearning house for dispercing individuals

master_migrant_template <- data.frame(sex = character(), age = integer(), genotype = character(), origin = integer(), 
                                      destination = integer(), disperse = logical())

setClass("pond", slots = list(pond_ID = 'numeric', size = 'numeric', routes = 'numeric', col_gen = 'numeric',
                              virgin = 'logical', route_rating = 'numeric',  pond_table = 'data.frame'))

### simulation  sub-functions ###

breed_clutch <- function(father, mother) {

  # takes maternal and paternal genotypes and then returns a random (Mendillian) vector of offspring genotypes
  
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
  
  return(clutch_geno_codes)
}

breed_pond <- function(pond){

  # Takes the current inhabitants of a pond, simulates breeding and then returns the original inhabitants + offspring
  # first make sure the pond acctualy has at least one adult of each sex, otherwise just return the existing inhabitants
  
  ad_males <- which(pond@pond_table$sex == 'male' & pond@pond_table$age > 1)
  ad_females <- which(pond@pond_table$sex == 'female' & pond@pond_table$age > 1)
  
  if(length(ad_males) == 0 || length(ad_females) == 0) {return(pond@pond_table)}
  
  breeding_threshold <- 100 * (1 - (pond@size / (length(ad_males) + length(ad_males))))
  
  # if the number of adults is greater than pond size, randomly reduce the chance of each adult breeding (so on average Nbreeders = pond size)
  
  breeding_impulse_males <- sample(1:100, length(ad_males), replace = TRUE)
  ad_males <- ad_males[which(breeding_impulse_males > breeding_threshold)]
  
  breeding_impulse_females <- sample(1:100, length(ad_females), replace = TRUE)
  ad_females <- ad_females[which(breeding_impulse_females > breeding_threshold)]
  
  breeders <- c(ad_males, ad_females)
  
  if(length(ad_males) == 0 || length(ad_females) == 0) {return(pond@pond_table)}

  num_AB <- sum(pond@pond_table$genotype[breeders] == 'BA' | pond@pond_table$genotype[breeders] == 'AB')
  num_NN <- sum(pond@pond_table$genotype[breeders] == 'NN')
  
  # if all breeders are genotype NN then all offspring must be too, so no need to bother simulating Mendillian genetics
  
  if(num_NN == length(breeders)) {

    num_offspring <- genotype_table$fecundicity[9] * length(ad_females)
    
    clutch_geno_codes <- rep(9, num_offspring)
  }
  
  # Likewise if all breeders are AB heterozygotes, all surviving offspring will be too (and there will be half as many of them)
  
  else if (num_AB == length(breeders) & optimise_AB == TRUE) {
    
    num_offspring <- (genotype_table$fecundicity[2] + genotype_table$fecundicity[1]) * length(ad_females) / 2
    
    clutch_geno_codes <- sample(c(2,4), num_offspring, replace = TRUE)
  }
  
  # if breeders are of mixed population, then we must simulate Mendillian segregation, and per-genotype embryonics mortality 
  
  else {
    
    male_genotypes <- pond@pond_table$genotype[ad_males]
    female_genotypes <- pond@pond_table$genotype[ad_females]

    male_attraction <- genotype_table$attraction[match(male_genotypes, genotypes)]
    
    female_mates <- sample(length(ad_males), length(ad_females), replace = TRUE, prob = male_attraction)
    
    clutch_geno_codes <- c()
    
    for(i in 1:length(female_mates)) {
      father_id <- female_mates[i]
      clutch_geno_codes <- c(clutch_geno_codes, breed_clutch(female_genotypes[i], male_genotypes[father_id]))
    }
    
    clutch_fitness <- genotype_table$embryo_fitness[clutch_geno_codes]
    
    clutch_luck <- sample(0:99, length(clutch_fitness), replace = TRUE)
    
    clutch_survial <- clutch_fitness + clutch_luck
    
    clutch_geno_codes <- clutch_geno_codes[which(clutch_survial >= 100)]
  }
  
  juvinial_num <- as.integer(length(clutch_geno_codes) / juvinial_cull_ratio) # simulate 95% hatchling mortallity 
  
  survivors <- clutch_geno_codes[sample(clutch_geno_codes, juvinial_num, replace = FALSE)]
  
  survivors_sex <- sample(c('male','female'), length(survivors), replace = TRUE)
  survivors_age <- rep(0, length(survivors))
  survivors_genotype <- genotype_table$genotypes[survivors]
  survivors_origin <- rep(pond@pond_ID, length(survivors))
  
  survivors_table <- data.frame(survivors_sex, survivors_age, survivors_genotype, survivors_origin)
  
  colnames(survivors_table) <- c('sex','age','genotype','origin')

  return(rbind (pond@pond_table, survivors_table))
}

age_pond <- function(pond){

  # take the inhibitants of a pond, randomly cull according to age and genotype and return 
  
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

  # take the inhibitants of a pond, randomly select 50% of 1 year olds to disperse, select their destination 
  # and rbind to the master migrant table, then return remaining inhabitants
  
  pond_newts <- pond@pond_table
  
  if(length(pond@routes) == 1) {return(pond_newts)}
  
  pond_newts$disperse <- sample(c(TRUE,FALSE), nrow(pond_newts), replace = TRUE, prob = c(disperse_percentage, 100 - disperse_percentage))
  migrants <- pond_newts[which(pond_newts$disperse == TRUE & pond_newts$age == 1),]
  
  # select destinations randomly using pre-calculated weights 
  
  migrants$destination <- sample(pond@routes, nrow(migrants), replace = TRUE, prob = pond@route_rating)
  master_migrant_table <<- rbind(master_migrant_table, migrants)
  
  pond_newts <- pond_newts[which(pond_newts$disperse == FALSE | pond_newts$age != 1),c('sex','age','genotype','origin')]
  return(pond_newts)
}

collect_newts <- function(pond,generation){

  # scoop up any incoming migrants from the master migrant table, add them to the ponds inhibatant table, and return
  
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
  
  # randomly generates a landscape, the position and size of ponds 
  # (also calculates graphical parameters for the visuliser script - should probably remove this code to a different function/script)
  
  pond_ID <- c(1:num_ponds)
  
  table_of_ponds <- data.frame(pond_ID)
  table_of_ponds$x_coord <- round(runif(n = num_ponds, min = 0, max = area_x))
  table_of_ponds$y_coord <- round(runif(n = num_ponds, min = 0, max = area_y))
  table_of_ponds$size <- sample(c(5,7,10,15,20,30,40,50,60), num_ponds, replace = TRUE, prob = c(9:1)) # this closely approximates Zipf's law 
  
  table_of_ponds$grid_x <- (table_of_ponds$x_coord / (area_x * 1.2)) + 0.1
  table_of_ponds$grid_y <- (table_of_ponds$y_coord / (area_x * 1.2)) + 0.1
  table_of_ponds$size_rt <- (sqrt(table_of_ponds$size) / 2000)
  
  table_of_ponds$total_pop <- 0
  table_of_ponds$num_NN <- 0
  table_of_ponds$num_BL <- 0
  table_of_ponds$num_other <- 0
  
  return(table_of_ponds)
}

connect_landscape <- function(table_of_ponds, max_distance = 1000) {

  # Returns TRUE if the landscape is fully traversable for a given dispersal radius
  
  num_ponds <- nrow(table_of_ponds)
  
  distance_matrix <- matrix(nrow = num_ponds, ncol = num_ponds)
  
  for( i in 1:num_ponds) {
    x_distances <- table_of_ponds$x_coord - table_of_ponds$x_coord[i]
    y_distances <- table_of_ponds$y_coord - table_of_ponds$y_coord[i]
    distances <- sqrt(x_distances^2 + y_distances^2)
    distance_matrix[,i] <- round(distances)
  }
  
  cluster <- c(1)
  
  cluster_index <- 1
  
  while(cluster_index <= length(cluster)) {
    
    pond_i <- cluster[cluster_index]
    pond_neigbours <- which(distance_matrix[,pond_i] < max_distance)
    cluster <- c( cluster, pond_neigbours[!(pond_neigbours %in% cluster)])
    cluster_index <- cluster_index + 1
  }
  
  #print(paste("number of ponds:", num_ponds))
  #print(paste("cluster size:", length(cluster)))
  
  if(length(cluster) == num_ponds) {
    #print ("fully connected map")
    return(TRUE)
  }
  else{return(FALSE)}
}

generate_ponds <- function(table_of_ponds) {

  # initialises empty populations tables for all ponds, also determines which ponds are neighbours and weights dispersal routes by size and distance
  
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

run_simulation <- function(landscape, pond_list, generations = 10, start_gen = 0, save_run = 2, 
                            run_num = 0, run_name = "Run", save_directory = "Sim_saves", BL_start_pop_thresh = 200){
  
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
  pond_AN_hist <- pond_pop_hist
  pond_BN_hist <- pond_pop_hist
  
  pond_pop_vec <- c()
  pond_BL_vec <- c()
  pond_NN_vec <- c()
  pond_AN_vec <- c()
  pond_BN_vec <- c()
  
  for(pond_num in 1:length(all_ponds)){
    
    pond_pop <- nrow(all_ponds[[pond_num]]@pond_table)
    pond_BL <- sum(all_ponds[[pond_num]]@pond_table$genotype == 'AB') + sum(all_ponds[[pond_num]]@pond_table$genotype == 'BA')
    pond_NN <- sum(all_ponds[[pond_num]]@pond_table$genotype == 'NN')
    pond_AN <- sum(all_ponds[[pond_num]]@pond_table$genotype == 'AN') + sum(all_ponds[[pond_num]]@pond_table$genotype == 'NA')
    pond_BN <- sum(all_ponds[[pond_num]]@pond_table$genotype == 'BN') + sum(all_ponds[[pond_num]]@pond_table$genotype == 'NB')
    
    pond_pop_vec <- c(pond_pop_vec, pond_pop)
    pond_BL_vec <- c(pond_BL_vec, pond_BL)
    pond_NN_vec <- c(pond_NN_vec, pond_NN)
    pond_AN_vec <- c(pond_AN_vec, pond_AN)
    pond_BN_vec <- c(pond_BN_vec, pond_BN)
  }
  
  pond_pop_hist <- rbind(pond_pop_hist, pond_pop_vec)
  pond_BL_hist <- rbind(pond_BL_hist, pond_BL_vec)
  pond_NN_hist <- rbind(pond_NN_hist, pond_NN_vec)
  pond_AN_hist <- rbind(pond_AN_hist, pond_AN_vec)
  pond_BN_hist <- rbind(pond_BN_hist, pond_BN_vec)
  
  colnames(pond_pop_hist) <- 1:length(all_ponds)
  colnames(pond_BL_hist) <- 1:length(all_ponds)
  colnames(pond_NN_hist) <- 1:length(all_ponds)
  colnames(pond_AN_hist) <- 1:length(all_ponds)
  colnames(pond_BN_hist) <- 1:length(all_ponds)
  
  if(sum(pond_BL_vec) > 0 ) {BL_exist <- TRUE}
  
  ### Iterate landscape ###

  for(i in 1:generations) {
    
    gen_counter <- gen_counter + 1
    
    if (gen_counter %% 50 == 0) {print(paste(run_num, "gen", gen_counter))}
    
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
    
    pond_pop_vec <- c()
    pond_BL_vec <- c()
    pond_NN_vec <- c()
    pond_AN_vec <- c()
    pond_BN_vec <- c()
    
    for(pond_num in 1:length(all_ponds)){
      
      pond_pop <- nrow(all_ponds[[pond_num]]@pond_table)
      pond_BL <- sum(all_ponds[[pond_num]]@pond_table$genotype == 'AB') + sum(all_ponds[[pond_num]]@pond_table$genotype == 'BA')
      pond_NN <- sum(all_ponds[[pond_num]]@pond_table$genotype == 'NN')
      pond_AN <- sum(all_ponds[[pond_num]]@pond_table$genotype == 'AN') + sum(all_ponds[[pond_num]]@pond_table$genotype == 'NA')
      pond_BN <- sum(all_ponds[[pond_num]]@pond_table$genotype == 'BN') + sum(all_ponds[[pond_num]]@pond_table$genotype == 'NB')
      
      pond_pop_vec <- c(pond_pop_vec, pond_pop)
      pond_BL_vec <- c(pond_BL_vec, pond_BL)
      pond_NN_vec <- c(pond_NN_vec, pond_NN)
      pond_AN_vec <- c(pond_AN_vec, pond_AN)
      pond_BN_vec <- c(pond_BN_vec, pond_BN)
    }
    
    pond_pop_hist <- rbind(pond_pop_hist, pond_pop_vec)
    pond_BL_hist <- rbind(pond_BL_hist, pond_BL_vec)
    pond_NN_hist <- rbind(pond_NN_hist, pond_NN_vec)
    
    Pop_total <- sum(pond_pop_vec)
    BL_total <- sum(pond_BL_vec)
    NN_total <- sum(pond_NN_vec)
    AN_total <- sum(pond_AN_vec)
    BN_total <- sum(pond_BN_vec)
    
    num_colonies <- sum(pond_pop_vec > 0)
    
    if(BL_total > Peak_BL){ Peak_BL <- BL_total}
    
    Other_total <- Pop_total - (BL_total + NN_total)
    
    fresh_colonies <- which(pond_pop_hist[i,] < 1 & pond_pop_hist[i + 1,] > 0)
    
    ### Initiate Balanced lethal system ###
    
    if( gen_counter >= BL_start & BL_exist == FALSE & length(fresh_colonies) > 0 & Pop_total >= BL_start_pop_thresh ) {
      
      if(length(fresh_colonies) > 1) {BL_start_pond <- sample(fresh_colonies,1)}
      else {BL_start_pond <- fresh_colonies[1]}

      all_ponds[[BL_start_pond ]]@pond_table$genotype[1] <- 'AB'
      
      print(paste("BL begins at pond:", BL_start_pond, "at gen:", gen_counter))
          
      BL_exist <- TRUE
      BL_start_gen <- gen_counter
      BL_start_total_pop <- Pop_total
      BL_start_pond_pop <- nrow(all_ponds[[BL_start_pond]]@pond_table)
      BL_start_colonies <- num_colonies
    }  
    
    #print(paste("Gen:", gen_counter ,"Total pop:", Pop_total, "Total BL, NN, mixed:", BL_total, NN_total, Other_total))
    
    ### check for finishing conditions ###
    
    if(Pop_total == 0){ 
      finish_message <- paste("Total extinction at gen:", gen_counter, "Peak_BL:", Peak_BL)
      finish_condition <- 1 
      if(save_run > 0) {save_run <- TRUE}
      break
    }
    
    if(Pop_total == BL_total && BL_total > 0){ 
      finish_message <- paste("Balanced lethal system fixed at gen:", gen_counter, "Peak_BL:", Peak_BL)
      finish_condition <- 2
      if(save_run > 0) {save_run <- TRUE}
      break
    }
    
    if(AN_total + BL_total + BN_total == 0 && BL_exist == TRUE && gen_counter > BL_start_gen){ 
      finish_message <- paste("Ancestral chromosome fixed at gen:", gen_counter, "Peak_BL:", Peak_BL)
      finish_condition <- 3
      if(save_run > 0 && Peak_BL > 200) {save_run <- TRUE}
      break
    }
    
    if(AN_total + BL_total == 0 && BL_exist == TRUE && gen_counter > BL_start_gen){ 
      finish_message <- paste("A chromosome extinct at gen:", gen_counter, "Peak_BL:", Peak_BL)
      finish_condition <- 3
      if(save_run > 0 && Peak_BL > 200) {save_run <- TRUE}
      break
    }
    
    if(BN_total + BL_total == 0 && BL_exist == TRUE && gen_counter > BL_start_gen){ 
      finish_message <- paste("B chromosome extinct at gen:", gen_counter, "Peak_BL:", Peak_BL)
      finish_condition <- 3
      if(save_run > 0 && Peak_BL > 200) {save_run <- TRUE}
      break
    }
    
    if(gen_counter == 25 && BL_exist == FALSE){ 
      finish_message <- paste("Failed to spawn balanced lethal system")
      finish_condition <- 4
      break
    }
  }
  
  ## After simulation iteration finishes ## 
  
  if(finish_condition == 0) {
    finish_message <- "still going!"
    if(save_run > 0) {save_run <- TRUE}
  }
  
  if(save_run > 1) {save_run <- FALSE}
  
  print(finish_message)
  
  ### save state ###
  
  if(save_run == TRUE) {
    
    print("saving")
  
    save_name <- paste(save_directory, "/", run_name, "_run_", run_num, "_class_", finish_condition, sep = "")
    
    run_info <- list(base_seed, iteration_seed, run_name, run_num, start_from_centre, juvinial_cull_ratio,
                     BL_start_gen, BL_start_total_pop, BL_start_pond, BL_start_pond_pop)
    names(run_info) <- c("base_seed", "iteration_seed", "run_name", "run_num", "start_from_centre", "juvinial_cull_ratio",
                         "BL_start_gen", "BL_start_total_pop", "BL_start_pond", "BL_start_pond_pop")
    
    state_list <<- list(run_info, genotype_table, landscape, pond_pop_hist, pond_BL_hist, pond_NN_hist, all_ponds)
    names(state_list) <<- c("run_info", "genotype_info", "landscape", "pop_history", "BL_history", "NN_history", "pond_data")
    
    save(state_list, file = save_name)
  }
  
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

print(paste("female hybrid fitness:", female_hybrid_relative_fitness))
print(paste("male hybrid fitness:", male_hybrid_relative_fitness))

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
  
  while( connect_landscape(iteration_landscape, 900) == FALSE) {
    iteration_landscape <- generate_landscape(num_ponds = iteration_pond_num, area_x = area_x_length, area_y = area_y_length)
  }
  
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