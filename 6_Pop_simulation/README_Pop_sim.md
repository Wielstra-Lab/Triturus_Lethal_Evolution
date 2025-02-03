### Simulation of Balanced Lethal Evolution ###

BL_sim_3.R is a script designed to model the evolution of a balanced lethal system in a population 
(of Triturus newts) with on locus with 3 possible alleles and so 9 possible genotypes - where the 
original state is NN and alleles A and B are derived. This is indented to simulate a balanced lethal 
system so genotypes AA and BB have 0 survival.

The simulation area consists of a landscape of between 50 to 150 sites (ponds) each with its own 
sub-population, where juvinial individuals can disperse between them. Ponds are located randomly within
a 5 x 5 km grid, and max dispersal distance is 1 km (but chance of dispersal falls off linearly with
distance).

Each simulation starts with 1 pond (on the left edge of the landscape) populated with 20 (10 male, 
10 female) adult individuals of the ancestral genotype (NN), with the other ponds being empty.
Normally a single individual in the NN population will be mutated to AB, however this script can also
simulate secondary contact between two populations, one with genotype NN fixed and the other fixed AB, 
in this case a second pond (on the right edge of the landscape) starts populated with 20 indivduals of 
genotype AB.

Input parameters are the number of generations of each iteration of the simulation, the number of
iterations, the fitness values for the various genotypes, and whether secondary contact is modeled. 

Each simulation has its own unique name and the random seed used can be set manually if desired. 
Overall statistics are recorded and the "successful" (where the A and B alleles survive) iterations 
are saved. This results in a save file that can be visualized with the accompanying script.

To avoid useless simulation time, each simulation is terminated after either the A or B chromosome
has gone extinct. To avoid creating many thousands of save files, simulations are only saved under
interesting conditions (i.e. when the population of AB individuals has reach at least 50 or the maximum 
number of generations has been reached without the balanced lethal system going extinct).

### Arguments ###

BL_sim_3.R accepts 6 arguments, which must be given in order.

1) Simulation name - determines the name of a the save file 
2) Number of itererations - the number of replicates simulated in this run of the simulation
3) Number of generations - the maximum number of generations in each replicate
4) Genotype fitness - determines the fitness parameters of the 'hybrid' (AN and BN) individuals (string "default" sets hybrid fitness to 0.75)
5) Input seed - if an integer is given then that will determine all 'random' value used in the simulation, allowing for exact repeats (string "random' will pick a seed randomly)
6) Secondary contact - if TRUE the a population with the AB genotype fixed is spawned on the right side of the landscape

### Example ###

``` Rscript BL_sim_3.R example_sim 100 1000 default random FALSE ```

Will simulate 100 replicates for 1000 generations each with the default (0.75) hybrid fitness, a random seed, and will not model seconday contact

if this run generated a random seed of 132421 and we wanted to repeat this run exactly, we could run:  

``` Rscript BL_sim_3.R example_sim 100 1000 default 132421 FALSE ```

### Runs used in study ###
 
``` Rscript BL_sim_3.R Sim_3_sec_hyb100_T1 100 1000 1.0 random TRUE ```  
100 replicates of secondary contact for up to 1000 generations with hybrid fitness equal to ancestral (1.0)  
_this run generated random seed: 4037782_


``` Rscript BL_sim_3.R Sim_3_sec_contact_T1 100 1000 default random TRUE ```  
100 replicates of secondary contact for up to 1000 generations with hybrid fitness set to default (0.75)  
_this run generated random seed: 6122263_

``` for rep in {1..5}; do Rscript ../BL_sim_3.R Sim_3_T1_$rep 10000 1000 default random FALSE; done ```  
5 x 10000 replicates of balanced lethal initation for up to 1000 generations with hybrid fitness set to default (0.75)  
_Sim_3_T1_1 generated random seed: 1088689_  
_Sim_3_T1_2 generated random seed: 6036874_  
_Sim_3_T1_3 generated random seed: 4955070_  
_Sim_3_T1_4 generated random seed: 96996_  
_Sim_3_T1_5 generated random seed: 8259011_  


### Outcomes ###

Each run of BL_sim_3.R produces an outcomes.txt file giving details of each replicate in the simulation run, the columns are:  

1) Replicate number
2) Final state (0 = simulation ran until max generation, 1 = Total extinction, 2 = BL system fixed, 3 = BL system extinct, 4 = BL mutation never occured)
3) Final generation of simulation
4) Generation Balanced Lethal mutation occured
5) Total Population at generation of Balanced Lethal mutation 
6) Population in the pond where the Balanced Lethal mutation occured at the generation it occured
7) Number of ponds colonised at the generation of Balanced Lethal mutation 
8) Number of ponds colonised at the end of Simulation
9) Number of ponds colonised at the generation where the number of individuals with the balanced lethal genotype reached its peak
10) Peak number of individuals with the balanced lethal genotype
11) Total number of ponds in the landscape

### Visulisation ###

The save files produced by BL_sim_3.R can be visualized with the script Visulise_2.R. This takes 5 arguments:

1) The input save file path
2) The output (pdf) file path
3) The desired number of maps
4) Whether a total population graph is desired (defaults to FALSE)
5) The final generation to draw a map for (defaults to final generation if none given)

So to make the plots shown in Fig. S6 from the save file 'Sim_3_sec_contact_T1_run_19_class_0' we ran:  

``` Rscript Visulise_2.R Sim_3_sec_contact_T1_run_19_class_0 Fig_S6_plots.pdf 5 TRUE ```

And to make the plots shown in Fig. S7 from the save file 'Sim_3_sec_contact_T1_run_19_class_0' we ran:  

``` Rscript Visulise_2.R Sim_3_T1_1_run_978_class_0 Fig_S7_plots.pdf 5 TRUE ```

### Dependancies ###

BL_sim_3.R uses only base R   
Visulise_2.R uses the ggplot2 grid and gridExtra packages



