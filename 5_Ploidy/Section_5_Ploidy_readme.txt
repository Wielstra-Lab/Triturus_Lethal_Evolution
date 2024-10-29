### Section 5: Ploidy analysis ###

.bam files from 30 F1 Triturus ivanbureschi x Triturus macedonicus samples (10 of each genotype) are segmented into sets of A-linked, B-linked and All genes.  
Each file is then processed in via the following pathway 

### Create the .bin file with nQuire, miniumum coverage (-q) is set to 10 and minimum mapping quality (-c) to 20 

for file in *chr12*.bam; do ~/nQuire/nQuire create -b $file -c 20 -q 10 -o ${file%.bam}.c20q10.X -x; done

### .bin file then denoised

for file in *.c20q10.X.bin; do ~/nQuire/nQuire denoise $file -o ${file%.bin}.denoised; done

### Assess the ploidy level using the histotest function and the denoised input .bin files

for file in *.c20q10.X.bin.denoised; do ~/nQuire/nQuire histotest $file >> bin_histotest_c20q10.txt; done

### Construct histograms for bam files

for file in *.c20q10.X.bin.denoised; do ~/nQuire/nQuire histo $file | cut -f 2 

### paste these together for multi-sample histograms for each gene class and genotype

### Run script faceted_ploidy.R to visualise the histograms

faceted_ploidy.R