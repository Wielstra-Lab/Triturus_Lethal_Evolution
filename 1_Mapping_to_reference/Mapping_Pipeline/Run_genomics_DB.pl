########################
### GenomicsDBImport ### 
########################

# This script is designed for combining multiple VCFs into a single joint VCF, if this hasn't already been done in the main pipeline

### Dependancies ###

# Perl must be installed (obviously)
# GATK must be installed and available as the path "gatk"

### Load Perl modules ###

use strict;                    # uses ' vars' , ' refs' , and ' subs' , to generate errors in case empty variables, symbolic references or bareword identifiers are used in improper ways
use warnings;                  # if triggered, it will indicate a ' problem '  exists, generates an error
use Cwd;
use Getopt::Std;

### Get Options ###

our $opt_d = getcwd; 
# The directory with the VCF data 
# If there is no "raw data" directory specified the script will default the current working directory

our $opt_t = ".raw.g.vcf"; 
# The terminator which marks the end of the 'sample name' in the VCF file names
# If no terminator is specified then the script will use characters before the string ".raw.g.vcf" as the sample name

our $opt_r = "/data1/projects/pi-vissermcde/Triturus_reference/triturusTargetsV1Oneliners_7139unique.fasta"; 
# The path to the reference fasta file that the data will be mapped against (if none specified the above default is used)

our $opt_i = "NULL"; 
# The path to the interval file used for making the combined gvcf, if none specified an file will be made from the reference fasta

our $opt_o = "Main_output.g.vcf";
# The name (and path to) of the final output combined g.vcf generated at the end of the script (if none specified the above default is used)

our $opt_x = "~/data1";
# The directory used for temprorary files used by some programmes (default is your data1 directory)

getopts("d:t:r:i:o:x:");

### Set major variables ###

my $vcf_dir = $opt_d;
my $terminator = $opt_t;
my $RefFile = $opt_r;
my $Ref_intervals = $opt_i;
my $output_vcf = $opt_o;
my $temp_dir = $opt_x;

### print settings as specified ### 

print "~~~ Settings as specified: ~~~ \n";
print "VCF directory: " . $vcf_dir . "\n";
print "sample name terminator string: " . $opt_t . "\n";
print "reference file path: " . $RefFile . "\n";
print "interval file path: " . $Ref_intervals . "\n";
print "output combined g.vcf name: " . $output_vcf . "\n";
print "temp directory: " . $temp_dir . "\n";
print "\n";

### make list of samples ###

opendir DIR,$vcf_dir;
my @vcfFiles = readdir DIR;
closedir DIR;

my @samples;

foreach my $file (@vcfFiles) {

    my $name;
	
    if ($file =~ /$terminator$/) { 
      ($name) = ($file =~ /^(.*)$terminator/);
      print "file name: " . $file . "\n Sample name: " . $name . "\n";
      push(@samples, $name);
    }
}

print "Samples included: \n";

foreach my $sample (@samples) {print $sample . "\n";}


### Make map file ###

print "Making Map File \n";

my $map_file = "map_file.txt";

foreach my $sample (@samples) {

    my $gVCF = $vcf_dir . "/" . $sample . ".raw.g.vcf";
    system("echo -e $sample'\t'$gVCF >> $map_file");
}

### Make interval file ###
# If no interval file is specified one will be made from the reference fasta

print "Making Interval File \n";

if ($Ref_intervals eq "NULL") {
    
    print "no interval specified, creating from reference file \n";
    print "old ref file: " . $Ref_intervals . "\n";
    $Ref_intervals = "interval_file.list";
    print "new ref file: " . $Ref_intervals . "\n";
    system("cat $RefFile | grep '>' | sed 's/>//' | sort > $Ref_intervals");
}

print $Ref_intervals . "\n";

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