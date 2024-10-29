
PIPELINE README:

The Pipeline_1.pl script is designed to process raw target capture data from multiple samples (typically as part
of the same project). The script will clean the data, map the samples against a reference genome, detect varients
from the reference genome and finally build a multisample g.vcf file containing all identified varients.

Raw data should come in the fastq.gz format (which is typically the format they are received in from, for example,
BaseClear). The raw data is also expected to come as two seperate files for each sample: R1 and R2. 

The script is capable of identifying and processing all files in a directory automatically, however to make things 
run smoothly it is best to ONLY have files you want to process with in this data directory (if there are additional
files hanging around the script might try to ingest them!)

Raw files typically come in a format like this:

BW_0406_JF-superpool-A_80025-307_TCGATCTCCATAGAAGCCAT_L001_R1_001_BHMGVYDRXY.filt.fastq.gz
BW_0406_JF-superpool-A_80025-307_TCGATCTCCATAGAAGCCAT_L001_R2_001_BHMGVYDRXY.filt.fastq.gz
BW_0407_JF-superpool-A_80025-308_TCGATCTCCACACTTAGGTA_L001_R1_001_BHMGVYDRXY.filt.fastq.gz
BW_0407_JF-superpool-A_80025-308_TCGATCTCCACACTTAGGTA_L001_R2_001_BHMGVYDRXY.filt.fastq.gz
BW_0408_JF-superpool-A_80025-309_TCGATCTCCAATAGGTAAGG_L001_R1_001_BHMGVYDRXY.filt.fastq.gz
BW_0408_JF-superpool-A_80025-309_TCGATCTCCAATAGGTAAGG_L001_R2_001_BHMGVYDRXY.filt.fastq.gz

It's important that the files have the "_R1" and "_R2" sections in as that's how the script will identify them as the
R1 and R2 files! (Also, having these strings of characters multiple times in the file name might cause issues, so avoid that).


RUNNING THE SCRIPT:

You probably want to run the script through SLURM, a template SLURM command script is included in the WiesltraLab GitHub. Make
sure to change your jobname, email address and set your raw data folder and any terminator you want to use (see below).

The script is included in the WielstraLab GitHub repository in the folder Pipeline. If you have linked this repository to your ALICE
/home/ directory you can run the script using the command:

perl ~/WielstraLab/Pipeline/Pipeline_1.pl

If you don't have the repository linked then run the script from wherever you installed it.

To run the script Perl must be loaded as well as several other modules (listed below). If running with the SLURM template, these will 
be loaded for you, except for BBmap which you will have to install yourself.

The script will use whatever directory it is called from as a working directory, so take that into account (you can set SLURM to change
directories before running the command).


USING A TERMINATOR TO SHORTEN FILE NAMES:

The script can deal with file names like the above just fine, however in this case it would use everything before the "_R1_" as 
the sample name. Which can make for quite long names e.g. 

BW_0406_JF-superpool-A_80025-307_TCGATCTCCATAGAAGCCAT_L001

To make things more readble we could use a terminator option to shorten the sample name. For all of these samples the code is
in the first part of the file name e.g. "BW_0406" and after that all contain the characters "_JF-superpool-A_80025". So if we
set the terminator to be "_JF" then the sample names would be shortened to:

BW_0406
BW_0407
BW_0408

Which is much more readable. Make sure that shortening the sample names in this keeps the samples unique (i.e. we should
not use "4" as a terminator in this example as then all these files will become "BW_0", which won't be very useful!).
Also try to make sure a terminator is in ALL the filenames, the script should be able to adapt if it isn't but it can make 
for some weird sample names.


USING OPTIONS WITH THE SCRIPT:

There are several settings that can be invoked when calling the script. In the SLURM template a few are used and you can
set the raw data directory and terminator in the lines above the main command. For a more complete overview of the script options
see the Mapping Pipeline Settings Guide in the WielstraLab GitHub repository.


DEPENDENCIES:

The script requires the following to be installed and active:

Perl
BBmap 
Trimmomatic
BWA
SAMtools
VCFtools
GATK

