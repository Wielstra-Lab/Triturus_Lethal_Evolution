### Section 5: BLASTing ###

# The markers on the Triturus linkage map are BLASTed against the genomes of Pleurodeles Waltl and Ambystoma Mexicanum

# First BLAST databases are constructed:

makeblastdb -dbtype nucl -in aPleWal1.pri.20220803.fasta -blastdb_version 4

makeblastdb -dbtype nucl -in GCA_002915635.3_AmbMex60DD_genomic.fna -blastdb_version 4

# Then a pipeline is used to BLAST the .fasta containing the mapped markers against the databases, this includes the script Filter_blast_2.R, used to filter paralogous hits

Trit_Amby_blast_slurm_4

Trit_Pluro_blast_slurm_1

