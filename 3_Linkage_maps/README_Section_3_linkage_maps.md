
### Filter Triturus and Lissotriton joint VCFs for use in LepMAP3 linkage analysis ###


##### GATK calls missing data as '0/0' in VCF files, instead of the standard './.', to filter for the missing data we need to remedy this, using BCFtools

```
bcftools +setGT Trit_joint_vcf/Main_output.g.vcf -- -t q -n . -i 'FMT/DP=0' > Trit_joint_vcf/Triturus_joint.g.vcf

bcftools +setGT Liss_joint_vcf/Main_output.g.vcf -- -t q -n . -i 'FMT/DP=0' > Liss_joint_vcf/Lissotriton_joint.g.vcf
```

 Now strictly filter the corrected joint VCF for minimum depth, mapping quality, missing data, minor allele frequency and remove indels and thin to a single SNP per marker

```
vcftools --vcf Triturus_joint.g.vcf --recode --recode-INFO-all --out Triturus_joint.vcf.filtered --maf 0.4 --min-meanDP 10 --minQ 20 --max-missing 0.9 --remove-indels --thin 1000

After filtering, kept 212 out of 212 Individuals
After filtering, kept 4939 out of a possible 126886 Sites

vcftools --vcf Lissotriton_joint.g.vcf --recode --recode-INFO-all --out Lissotriton_joint.vcf.filtered --maf 0.4 --min-meanDP 10 --minQ 20 --max-missing 0.9 --remove-indels --thin 1000

After filtering, kept 205 out of 205 Individuals
After filtering, kept 4540 out of a possible 216530 Sites
```

### Build Lissotriton linkage map ###

Build map using LepMAP3 pipeline, using settings specified in the SLURM job:  

``` Liss_linkage_map_All_slurm ```

 Separate chromosomes: lodLimit = 22 distortionLod = 1  
 Join singles: lodLimit = 15 lodDifference = 5 iterate = 1  
 Order markers: sexAveraged = 1 numMergeIterations = 12 numPolishIterations = 8 minError = 0.01  

Map is then ordered by length (as opposed to number of markers)

###  Build Main Triturus linkage map ###

 Build map using modified LepMAP3 pipeline including presence/absence call data:

``` Triturus_linkage_map_All_slurm ```

Separate chromosomes: lodLimit = 27 distortionLod = 1
Join singles: lodLimit = 20 lodDifference = 5 iterate = 1
Order markers: sexAveraged = 1 numMergeIterations = 12 numPolishIterations = 8 minError = 0.01

Identify which linkage group presence/absence markers are clustered on:

``` grep ^DN Trit_map_all.ordered.mapped ```

_63 presence/absence markers mapped, all on raw linkage group 3_  

Re-order with presence/absence group (raw group 3) prioritized as group 1, and the other groups ordered by length


###  Chromosome 1 Triturus linkage maps ###


 Make list of SNPs and presence/absence markers mapped in group 1:  

```
grep 3$ Trit_map_all.ordered.mapped | grep TRIN | cut -f 1 > Trit_map_Main.group_1_SNPs

grep ^DN Trit_map_all.ordered.mapped | cut -f 1 > Trit_map_Main.group_1_PA
```

Subset input VCF with only group 1 SNPs, and exclude presence/absence markers which failed to map. Copy, compress, and index the Triturus joint vcf so BCFtools can process it, then use BCFtools view, to extract desired SNP data

```
cp Triturus_joint.vcf.filtered.recode.vcf Triturus_joint.vcf.filtered.recode.vcf.copy
bgzip Triturus_joint.vcf.filtered.recode.vcf.copy
tabix Triturus_joint.vcf.filtered.recode.vcf.copy.gz

bcftools view -R Trit_map_Main/Trit_map_Main.group_1_SNPs -O v -o Triturus_joint.vcf.C1_SNPs.vcf Triturus_joint.vcf.filtered.recode.vcf.copy.gz
```

Filter the PA call files for presence/absence markers which mapped 

```
cat Trit_map_Main/Trit_map_Main.group_1_PA | while read line; do grep $line Triturus_PA_screened_call_type_1; done > Triturus_PA_screened_call_type_1_mapped

cat Trit_map_Main/Trit_map_Main.group_1_PA | while read line; do grep $line Triturus_PA_screened_call_type_2; done > Triturus_PA_screened_call_type_2_mapped
```

Create two new linkage maps based only on group 1 SNPs, with type 1 and type 2 presence absence markers included respectively using same settings as for the main map

```
Triturus_linkage_map_type1_slurm

Triturus_linkage_map_type2_slurm
```

