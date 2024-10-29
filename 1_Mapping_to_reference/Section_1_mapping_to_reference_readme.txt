###################################
### SECTION 1: Mapping raw data ###
################################### 

# The reference sequences used can be found in the associated Zenodo repository: 10.5281/zenodo.14008530 

# Raw reads were initially mapped against reference sequence with script Pipeline_1.pl e.g.:

perl Pipeline_1.pl -d "/data1/francejm/SuperPoolA" -b "/home/francejm/BB/bbmap" -c -v 

# This produced a .bam file and a a g.vcf for each sequenced library
# Coverage was assessed with script PeakShell_2 and Peakloop2.R e.g:

cd /data1/francejm/Final_map/SuperPoolA/Dedup
for bam in *bam; do sh PeakShell_2 $bam; done

### g.vcf files (and index files) were moved into a single directory for each mapping family, and names corrected to match sample IDs

cat Triturus_map_samples.txt | while read pool data cov type samp; do cp ~/data1/Final_map/$pool/VCF/$data.raw.g.vcf ~/data1/Final_map/Triturus_map_vcf/$samp.raw.g.vcf; done
cat Triturus_map_samples.txt | while read pool data cov type samp; do cp ~/data1/Final_map/$pool/VCF/$data.raw.g.vcf.idx ~/data1/Final_map/Triturus_map_vcf/$samp.raw.g.vcf.idx; done

cat Lissotriton_map_samples.txt | while read pool data cov type samp; do cp ~/data1/Final_map/$pool/VCF/$data.raw.g.vcf ~/data1/Final_map/Lissotriton_map_vcf/$samp.raw.g.vcf; done
cat Lissotriton_map_samples.txt | while read pool data cov type samp; do cp ~/data1/Final_map/$pool/VCF/$data.raw.g.vcf.idx ~/data1/Final_map/Lissotriton_map_vcf/$samp.raw.g.vcf.idx; done

### Joint genotyping for each mapping family performed with script Run_genomics_DB.pl 

perl Run_genomics_DB.pl -d ~/data1/Final_map/Triturus_map_vcf 
perl Run_genomics_DB.pl -d ~/data1/Final_map/Lissotriton_map_vcf/ 
