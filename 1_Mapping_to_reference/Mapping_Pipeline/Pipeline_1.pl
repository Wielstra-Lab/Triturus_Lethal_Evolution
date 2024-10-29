#!/usr/bin/perl -w                              ###this enables useful warnings
#Written by M.C. de Visser (with input from James France)
#Based on scripts written by E. McCartney-Melstad

# Master Pipeline 1.0 script: 
# Will take a folder with a number of fastq.gz files (compressed, demultiplexed illumina NGS data) from target capture samples.
# Will clean the data, map the data to the specified reference, call variants, and combine resulting vcfs into a single g.vcf for all samples.
# Ideally run on HPC, Leiden University's ALICE cluster used for developement (ALICE slurm script attached)

####################
### Dependencies ###
####################

# Script requires the following to be installed and active:

# Perl
# BBmap 
# Trimmomatic
# BWA
# SAMtools
# VCFtools
# GATK

##########################
### Prepare enviroment ###
##########################

### Load Perl modules ###

use strict;                    # uses ' vars' , ' refs' , and ' subs' , to generate errors in case empty variables, symbolic references or bareword identifiers are used in improper ways
use warnings;                  # if triggered, it will indicate a ' problem '  exists, generates an error
use Parallel::ForkManager;     # used to perform a number of repetitive operations withing a single Perl script
use Cwd;
use Getopt::Std;
use File::Copy;

### Get Options ###

our $opt_c = 0;
# If -c is used, all the intermediate files (and directories) before the final deduplicated bams will be deleted after use

our $opt_v = 0;
# If -v is used, the script will not attempt to generate the final combined g.vcf

our $opt_d = getcwd; 
# The directory with the raw fastaq.gz data 
# If there is no "raw data" directory specified the script will default to the current working directory

our $opt_t = "_R1"; 
# The terminator which marks the end of the 'sample name' in the raw file names
# If no terminator is specified then the script will use characters before the string "_R1" as the sample name

our $opt_a = "/data1/projects/pi-vissermcde/Triturus_reference/TruSeq2-PE_trimmomatic.fa";
# The path to the adapter file used in the trimming step (if none specified the above default is used)

our $opt_r = "/data1/projects/pi-vissermcde/Triturus_reference/triturusTargetsV1Oneliners_7139unique.fasta"; 
# The path to the reference fasta file that the data will be mapped against (if none specified the above default is used)

our $opt_i = "NULL"; 
# The path to the interval file used for making the combined gvcf, if none specified an file will be made from the reference fasta

our $opt_b = "~/bbmap";
# The path to the installation of bbmap (if none specified the above default is used)

our $opt_o = "Main_output.g.vcf";
# The name (and path to) of the final output combined g.vcf generated at the end of the script (if none specified the above default is used)

our $opt_p = 15;
# The number of simultaneous processes to use in the steps managed by parallel fork manager (default is 15)

our $opt_x = "~/data1";
# The directory used for temprorary files used by some programmes (default is your data1 directory)

getopts("cvd:t:a:r:i:b:o:p:x:");

### Set major variables ###

my $adapterFile = $opt_a;
my $RefFile = $opt_r;
my $Ref_intervals = $opt_i;
my $BBdukscript = $opt_b . "/bbduk.sh";
my $BBrepairscript = $opt_b . "/repair.sh";
my $output_vcf = $opt_o;

# Make sure you use the reference file that you are interested in. We use the exon reference with 7139 unique targets:

### Assign and create directories ###

my $raw_dir = $opt_d;
my $temp_dir = $opt_x;

my $unzipped_dir = getcwd . "/Unzip";
my $BBduk_dir = getcwd . "/BBduk";
my $trimmed_dir = getcwd . "/Trimmed";
my $repair_dir = getcwd . "/Repaired";
my $mapped_dir = getcwd . "/Mapped";
my $dedup_dir = getcwd . "/Dedup";
my $vcf_dir = getcwd . "/VCF";

mkdir($unzipped_dir);
mkdir($BBduk_dir);
mkdir($trimmed_dir);
mkdir($repair_dir);
mkdir($mapped_dir);
mkdir($dedup_dir);
mkdir($vcf_dir);

### print settings as specified ### 

print "~~~ Settings as specified: ~~~ \n";
print "raw directory: " . $raw_dir . "\n";
print "sample name terminator string: " . $opt_t . "\n";
print "adapter file path: " . $adapterFile . "\n";
print "reference file path: " . $RefFile . "\n";
print "interval file path: " . $Ref_intervals . "\n";
print "BBduk script path: " . $BBdukscript . "\n";
print "BBrepair script path: " . $BBrepairscript . "\n";
print "output combined g.vcf name: " . $output_vcf . "\n";
print "number of threads: " . $opt_p . "\n";
print "temp directory: " . $temp_dir . "\n";
print "Cull intermediates: " . $opt_c . "\n";
print "terminate before combining vcfs: " . $opt_v . "\n";
print "\n";

#######################
### Specify samples ###
#######################

# This will generate a list of sample names from the file names in the "raw data" directory
# Any file name with the string "_R1" will be translated into a sample name and added to the list
# Sample names will be a shortened version of the file name, either up to the specified terminator character or to the first 20 characters

opendir DIR,$raw_dir;
my @Rawfiles = readdir DIR;
closedir DIR;

my @samples;

foreach my $file (@Rawfiles) {

  my $name;
	
	if ($file =~ /_R1/) {
		
		if ($opt_t eq "NULL" or $file !~ /$opt_t/) { $name = substr($file, 0, 20); } # if terminator is "NULL" or missing from file the sample name is first 20 characters
		else { ($name) = ($file =~ /^(.*)$opt_t/); } # gets the sample name before terminator
		
		push(@samples, $name);
	  }
	} 

### List samples ###

print "Processing " . scalar(@samples) . " samples\n";
foreach my $sample (@samples) {print $sample . "\n";}

######################
### Copy and Unzip ###
######################

# Copy files to the unzipped directory
# this also adjusts the file name to the format Sample_name_R1/2_fastq.gz

print "copying and unzipping raw data \n";
system("DATE=\$(date) \n echo \"Time is now: \$DATE\" ");
 
foreach my $sample (@samples) {
	
  foreach my $file (@Rawfiles) {
		
		if ($file =~ /$sample.*R1/) {
    
      my $old_location = $raw_dir . "/" . $file;
      my $new_location =  $unzipped_dir . "/" . $sample . "_R1.fastq.gz" ;

      print "Copying " . $old_location . " to " . $new_location . "\n";
     
      copy($old_location, $new_location );
      }   
    
		elsif ($file =~ /$sample.*R2/) {
    
      my $old_location = $raw_dir . "/" . $file;
      my $new_location =  $unzipped_dir . "/" . $sample . "_R2.fastq.gz" ;

      print "Copying " . $old_location . " to " . $new_location . "\n";
     
      copy($old_location, $new_location );
      }
    }
  }

# Unzip files in the unzipped directory
# Drops the .gz from the file name

opendir DIR,$unzipped_dir;
my @zipped_files = readdir DIR;
closedir DIR;

my $unzipFM = Parallel::ForkManager->new($opt_p);
foreach my $file (@zipped_files) {
    $unzipFM->start and next;
    print "unzipping " . $file . "\n";
    system("gunzip " . $unzipped_dir . "/" . $file);
    $unzipFM->finish;
}
$unzipFM->wait_all_children();

print "finished unzipping files \n";


####################################
### Trim 151st bp off with BBDuk ###
####################################

# First we want to get rid of the 151st bp that many reads have, which could be prone to error.

my @BBdukCommands;
foreach my $sample (@samples) {
    my $R1 = $unzipped_dir . "/" . $sample . "_R1.fastq";
    my $R2 = $unzipped_dir . "/" . $sample . "_R2.fastq";
    
    my $BBdukBaseNameR1 = $BBduk_dir . "/" . $sample . "_150_R1.fastq";
    my $BBdukBaseNameR2 = $BBduk_dir . "/" . $sample . "_150_R2.fastq";

    push(@BBdukCommands, "$BBdukscript -Xmx4g in=$R1 out=$BBdukBaseNameR1 ftr=149");
    push(@BBdukCommands, "$BBdukscript -Xmx4g in=$R2 out=$BBdukBaseNameR2 ftr=149");
}

print ">>>>Running all BBduk commands\n";

my $BBdukFM = Parallel::ForkManager->new($opt_p);
foreach my $BBdukCommand(@BBdukCommands) {
    $BBdukFM->start and next;
    print "Running the following command: \n$BBdukCommand\n";
    system($BBdukCommand);
    $BBdukFM->finish;
}
$BBdukFM->wait_all_children();

print "\n\n>>>>Finished running all BBduk commands\n\n";

if ($opt_c == 1) {
    print "\n\n>>>> Deleting unzipped directory \n\n";
    system("rm -r $unzipped_dir");
}


################
### TRIMMING ###
################

# Here we trim universal adapters and low quality bases/reads
# We use Trimmomatic to do so (this has been adjusted from the version by E. McCartney-Melstad, which used Skewer

my $Trim_FM = Parallel::ForkManager->new($opt_p);
print "\n\nTrimming for adapter contamination, read length and quality using Trimmomatic!\n\n";


foreach my $sample (@samples) {
    $Trim_FM->start and next;
    
    my $R1 = $BBduk_dir . "/" . $sample . "_150_R1.fastq";
    my $R2 = $BBduk_dir . "/" . $sample . "_150_R2.fastq";
    
    unless (-e $adapterFile) {die "$adapterFile not present!\n";}
    
    my $OP1 = $trimmed_dir . "/" . $sample . "_150_R1_P_trimmomatic.fastq";
    my $OU1 = $trimmed_dir . "/" . $sample . "_150_R1_U_trimmomatic.fastq";
    my $OP2 = $trimmed_dir . "/" . $sample . "_150_R2_P_trimmomatic.fastq";
    my $OU2 = $trimmed_dir . "/" . $sample . "_150_R2_U_trimmomatic.fastq"; 
    
# These are some default/basic settings that can be altered if needed. 
# In the paper by E. McCartney-Melstad et al, a trailing threshold of 15 was used, but we use 5 for both leading and trailing
    system("DATE=\$(date) \n echo \"Begin trimming $sample at \$DATE\" ");
    system("java -jar \$EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE $R1 $R2 $OP1 $OU1 $OP2 $OU2 ILLUMINACLIP:$adapterFile:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:5:20 MINLEN:50 -threads 12");
    system("DATE=\$(date) \n echo \"Done trimming $sample at \$DATE\" ");

    $Trim_FM->finish;
}
$Trim_FM->wait_all_children();
print "\n\n >>>>Finished trimming with Trimmomatic\n\n";

if ($opt_c == 1) {
    print "\n\n>>>> Deleting BBduk directory \n\n";
    system("rm -r $BBduk_dir");
}


###########################
### SORT TRIMMED FASTQs ###
###########################

# This will make sure the R1 and R2 FASTQ files obtained after BBduk + Trimming are actually sorted/in sync again.
# This is necessary because otherwise BWA will not recognize read pairs that belong to each other, and hence mapping may fail


print "\n\n >>>>Sorting trimmed fastqs using BBmap\n\n";

my $BBrepairFM = Parallel::ForkManager->new($opt_p);

foreach my $sample (@samples) {
    $BBrepairFM ->start and next;

    my $read_1_in = $trimmed_dir . "/" . $sample . "_150_R1_P_trimmomatic.fastq";
    my $read_2_in = $trimmed_dir . "/" . $sample . "_150_R2_P_trimmomatic.fastq";
    
    my $read_1_out = $repair_dir . "/" . $sample . "_150_R1_P_trimmomatic_fixed.fastq";
    my $read_2_out = $repair_dir . "/" . $sample . "_150_R2_P_trimmomatic_fixed.fastq";

    system($BBrepairscript . " in1=" . $read_1_in . " in2=" . $read_2_in . " out1=" . $read_1_out . " out2=" . $read_2_out);
    
    $BBrepairFM ->finish;
}
$BBrepairFM ->wait_all_children();

print "\n\n >>>>Finished sorting trimmed fastqs using BBmap\n\n";
system("DATE=\$(date) \n echo \"Done sorting at \$DATE\" ");

if ($opt_c == 1) {
    print "\n\n>>>> Deleting trimming directory \n\n";
    system("rm -r $trimmed_dir"); 
}


###############
### MAPPING ###
###############

### Here we align the reads to a reference / we map the reads and create BAM files, which are compressed SAM files
### Adjust the forkmanager number appropriately.


print "\n\n >>>>Mapping reads using BWA MEM\n\n";

my $bwaFM = Parallel::ForkManager->new($opt_p);
foreach my $sample (@samples) {
    $bwaFM->start and next;
    
    my $R1 = $repair_dir . "/" . $sample . "_150_R1_P_trimmomatic_fixed.fastq";
    my $R2 = $repair_dir . "/" . $sample . "_150_R2_P_trimmomatic_fixed.fastq";
    my $bam = $mapped_dir . "/" . $sample . ".bam";

    system("DATE=\$(date) \n echo \"Mapping $sample at \$DATE\" ");
    system("bwa mem -M -v 1 $RefFile $R1 $R2 | samtools view -bS - > $bam");
    system("DATE=\$(date) \n echo \"Mapped $sample at \$DATE\" ");
    
    $bwaFM->finish; 
}
$bwaFM->wait_all_children();

print "\n\n >>>>Finished mapping reads using BWA mem\n\n";

if ($opt_c == 1) {
    print "\n\n>>>> Deleting repair directory \n\n";
    system("rm -r $repair_dir");
}


########################################
### Add read-groups and de-duplicate ### 
########################################

# Here we'll add RG information using picard AddOrReplaceReadGroups and mark duplicates, adjust the forkmanager number appropriately.
# Be aware that with these settings, Picard does not distinguish between different runs or lanes
# so if you leave these settings for the whole bulk/all your samples, there will be no correction of in-between-run technical biases.
# If you do want to do that, you will need to adjust the RGLB and PU information in more detail and run it for each sample/batch/run separately/specifically.
# Here We use RGLB -lib1, so everything is assumed to be a unique/one library (if one sample occurs twice under a different name, it will not be recognized as the same sample obviously
# also if you combined reads from different runs of one sample into one fastq file, of course Picard is also unable to recognize reads from one or the other run, so be careful with that).

my $addReplaceFM = Parallel::ForkManager->new($opt_p);
print "\n\nAdding read groups and marking duplicates with picard\n\n";

foreach my $sample (@samples) {
    $addReplaceFM->start and next;

    my $input = $mapped_dir . "/" . $sample . ".bam";
    my $output = $mapped_dir . "/" . $sample . ".RG.bam";
    my $SM = $sample;
    my $RGLB = $SM . '-lib1';
    my $RGID = $sample;

    system("DATE=\$(date) \n echo \"Adding read groups to $sample at \$DATE\" ");
    system("gatk AddOrReplaceReadGroups I=$input O=$output RGLB=$RGLB RGPL=ILLUMINA RGSM=$SM RGID=$RGID RGPU=NA SORT_ORDER=coordinate");
    system("DATE=\$(date) \n echo \"Done adding read groups, and begining deduplication for $sample at \$DATE\" ");
    
    my $MDout = $dedup_dir . "/" . $sample . ".dedup.bam";
    my $metrics = $dedup_dir . "/" . $sample . ".metrics";
    
    system("gatk MarkDuplicates I=$output O=$MDout M=$metrics");
    system("samtools index $MDout");
    system("DATE=\$(date) \n echo \"Done deduplicating $sample at \$DATE\" ");
    
    if ($opt_c == 1) {unlink($input, $output);}
    
    $addReplaceFM->finish;
}

$addReplaceFM->wait_all_children();
print "\n\n >>>>Finished adding read groups and marking duplicates with picard\n\n";

if ($opt_c == 1) {
    print "\n\n>>>> Deleting mapping directory \n\n";
    system("rm -r $mapped_dir");
}


#######################
### VARIANT CALLING ### 
#######################

# Here we will do the first round of variant calling, which will generate raw g.vcf files
# Adjust the forkmanager number accordingly.

print "\n\nRunning haplotypeCaller to generate pre-BQSR g.vcfs\n\n";

my $hapCallerFM = Parallel::ForkManager->new($opt_p);
foreach my $sample (@samples) {
    $hapCallerFM->start and next;
    
    my $inputBAM = $dedup_dir . "/" . $sample . ".dedup.bam";
    my $gVCF = $vcf_dir . "/" . $sample . ".raw.g.vcf";

    system("DATE=\$(date) \n echo \"Begining calling haplotypes for $sample at \$DATE\" ");
    system("gatk --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' HaplotypeCaller --reference $RefFile -ERC GVCF -I $inputBAM -O $gVCF --tmp-dir $temp_dir");
    system("DATE=\$(date) \n echo \"Done calling haplotypes for $sample at \$DATE\" "); 
    
    $hapCallerFM->finish;
}

$hapCallerFM->wait_all_children();

print "\n\nFinished running first round of haplotypeCaller to generate pre-BQSR, raw gvcfs\n\n";

if ($opt_v == 1) {
    print "\n\n>>>> Ending script as commanded before joint genotyping \n\n";
    die;
}


########################
### GenomicsDBImport ### 
########################

### Make map file ###

my $map_file = "map_file.txt";

foreach my $sample (@samples) {

    my $gVCF = $vcf_dir . "/" . $sample . ".raw.g.vcf";
    system("echo -e $sample'\t'$gVCF >> $map_file");
}

### Make interval file ###
# If no interval file is specified one will be made from the reference fasta

if ($Ref_intervals eq "NULL") {
    
    print "no interval specified, creating from reference file \n";
    print "old ref file :" . $Ref_intervals . "\n";
    $Ref_intervals = "interval_file.list";
    print "new ref file :" . $Ref_intervals . "\n";
    system("cat $RefFile | grep '>' | sed 's/>//' | sort > $Ref_intervals");
}

### Combine VCFs into database ###

print "\n\nCombining gvcfs with GDBI to generate pre-BQSR ms-gVCF file to use for joint-genotype calling\n\n";
system("DATE=\$(date) \n echo \"Start at \$DATE\" ");

my $database = getcwd . "/gvcf_db";

system("gatk --java-options '-Xms800m -Xmx100g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' GenomicsDBImport --sample-name-map $map_file --max-num-intervals-to-import-in-parallel 20 --genomicsdb-workspace-path $database --intervals  $Ref_intervals --merge-contigs-into-num-partitions 1 --tmp-dir $temp_dir");

print "\n\nFinished combining gvcfs with GDBI to generate ms-gVCF file to use for joint-genotype calling\n\n";

## Here we'll perform joint-genotype calling as to truly 'fill in' the merged gVCF file with the correct genotypes/SNP calls ###

print "\n\nPerforming joint-genotype calling using the created database\n\n";
system("DATE=\$(date) \n echo \"Start at \$DATE\" ");

system("gatk --java-options '-Xms800m -Xmx110g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' GenotypeGVCFs -R $RefFile -V gendb://$database -O $output_vcf --tmp-dir $temp_dir");

print "\n\nFinished joint-genotype calling, created pre-BQSR raw ms-gVCF\n\n";
system("DATE=\$(date) \n echo \"Finished at \$DATE\" ");