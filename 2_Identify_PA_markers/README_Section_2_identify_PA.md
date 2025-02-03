
## SECTION 2 Identification of presence absence markers in Triturus data set, using per maker coverage data ##


Group all coverage data into a single directory and correct names to match samples ID  

```
cat Triturus_map_samples.txt  while read pool data cov type samp; do cp $poolDedup$data.dedup.bam_PeakCoverage Triturus_map_cov$samp.PeakCoverage; done
```

 Run Create_Cov_table_2.R on coverage data from arrested embryos to make coverage matrix, and again for all Triturus samples  

```
Rscript Create_Cov_table_2.R Target_names_short.txt Arrested_Triturus_samples.txt Triturus_map_cov Triturus_arrested_coverage_table.txt
Rscript Create_Cov_table_2.R Target_names_short.txt All_Triturus_samples_vcf_ordered.txt Triturus_map_cov Triturus_all_vcf_coverage_table.txt
```

Run Select_presence_absence_4.R to create lists of Type 1 and 2 markers, as well identify which embryos are heterozygote or Type 1 or Type 2 homozygote   

```
Rscript Select_presence_absence_1.R Triturus_arrested_coverage_table.txt Triturus_PA/Triturus
```

Run Screen_cov_stats.R on both type 1 and type 2 marker lists, to double check for any markers selected in error in the previous step, combine lists into one  

```
Rscript Screen_cov_stats.R Triturus_all_vcf_coverage_table.txt Triturus_PA/Triturus_Type_1_markers Triturus_PA/Triturus_Type_1_markers

Rscript Screen_cov_stats.R Triturus_all_vcf_coverage_table.txt Triturus_PA/Triturus_Type_2_markers Triturus_PA/Triturus_Type_2_markers

cat Triturus_PA/Triturus_Type_1_markers_screened.txt Triturus_PA/Triturus_Type_2_markers_screened.txt  Triturus_PA/Triturus_PA_markers_screened.txt
```

Run Add_presence_absence_to_call_table_4.R to convert the coverage data A and B linked presence absence markers from of the entire Triturus linkage map family into pseudoSNP calls in a format suitable for LepMAP3  

```
Rscript Add_presence_absence_to_call_table_4.R Triturus_all_vcf_coverage_table.txt Triturus_PA/Triturus_PA_markers_screened.txt Triturus_PA_screened_call
```
