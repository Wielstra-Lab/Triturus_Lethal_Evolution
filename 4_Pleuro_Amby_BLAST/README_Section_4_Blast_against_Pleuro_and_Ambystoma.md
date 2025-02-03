## Section 5: BLAST searches ###

#### The markers on the Triturus linkage map are BLASTed against the genomes of Pleurodeles Waltl and Ambystoma Mexicanum

 First BLAST databases are constructed:

``` makeblastdb -dbtype nucl -in aPleWal1.pri.20220803.fasta -blastdb_version 4 ```

``` makeblastdb -dbtype nucl -in GCA_002915635.3_AmbMex60DD_genomic.fna -blastdb_version 4 ```

 Then a pipeline is used to BLAST the .fasta containing the mapped markers against the databases, this includes the script Filter_blast_2.R, used to filter paralogous hits.  
 
 An  overview is as follows:
 
 ``` 
 blastn -query Trit_map_Main.trimmed.sequences.fasta  -db /data1/francejm/Ambystoma/GCA_002915635.3_AmbMex60DD_genomic.fna -outfmt 6 -word_size 15 -evalue 1e-10 -num_threads 4 > Trit_TC_amby_blast_4_out_raw.txt  
 
 sort -k1,1 -k12,12nr Trit_TC_amby_blast_4_out_raw.txt > Trit_TC_amby_blast_4_out_sorted.txt  
 
 Rscript ~/JaFrance/Linkage_Map_Scripts/Filter_blast_2.R Trit_TC_amby_blast_4_out_sorted.txt Trit_TC_amby_blast_4_out_filtered.txt 
 ```   
 
 ##### Slurm sheduler files for each BLAST search:  

_Trit_Amby_blast_slurm_4_

_Trit_Pluro_blast_slurm_1_

