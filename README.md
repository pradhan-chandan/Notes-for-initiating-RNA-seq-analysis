# Notes-for-initiating-RNA-seq-analysis-from-plant-fungal-system
Practical workflow documentation for initiating RNA-Seq data analysis from raw FASTQ files to get hands-on practice



For unzipping zipped or compressed files using terminal
	Type “gunzip [filename.format]” like “gunzip Ecoracana_560.fa.gz”



For fastqc running
	To generate fastqc results for a number of rna_seq files at a time:
Keep all the files in a single folder
And run fastqc with “fastqc *.fq.gz”
.fq.gz signifies the file format which is saved in your folder
The program will run the files based on alphabetical order and will store the output files in the same directory.



For trimming the sequences using trimmomatic
	For trimming your seq_reads to remove the unwanted bases or adaptor contaminations from it, use the following code:

“trimmomatic PE -threads 4 sample_R1.fastq.gz sample_R2.fastq.gz sample_R1_paired.fastq.gz sample_R1_unpaired.fastq.gz sample_R2_paired.fastq.gz sample_R2_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:50”

trimmomatic PE: This specifies that you are running Trimmomatic in paired-end mode. If you are using single-end reads, replace PE with SE.
-threads 4: This specifies the number of threads to use. Adjust this number based on your CPU cores.
input_forward.fq.gz and input_reverse.fq.gz: The input files for forward and reverse reads.
output_forward_paired.fq.gz, output_forward_unpaired.fq.gz, output_reverse_paired.fq.gz, output_reverse_unpaired.fq.gz: The output files for paired and unpaired reads.
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10: Adapter clipping, using the TruSeq3-PE.fa file for the adapter sequences.
LEADING:5: Remove leading low-quality bases.
TRAILING:5: Remove trailing low-quality bases.
SLIDINGWINDOW:4:20: Perform sliding window trimming, cutting once the average quality within the window falls below a threshold.
MINLEN:50: Drop reads that are shorter than 50 bases after trimming.



For performing genome indexing
	To create genome index files of the organisms:
Prepare separate folders for each organism
Download the fasta file ‘.fa’ and annotation files ‘.gtf/.gff’ into the respective folders from “https://plants.ensembl.org/index.html” or  “https://fungi.ensembl.org/index.html” or “https://phytozome-next.jgi.doe.gov/” or NCBI or “https://asia.ensembl.org/index.html” 
Now move to the terminal or command prompt
Download Anaconda latest version (always check what processor you’re using when working with mac, as it creates mess with compatibility and deny installation) to use ‘conda command line’
Create channels, create environment (env)
Activate the env, install the required packages (take codes to install from anaconda website)/ “https://anaconda.org/bioconda/star”
Follow the below steps to continue

Commands to perform genome indexing
	For --genomesDAindexNbases the number is max 14, log2(genome length)/2-1), (14, 10.8(round down)). First do the log2 value to the genome size in bytes, then divide by 2 and then subtract 1 from it. Round down the value and if possible take even one lower value.

“STAR --runThreadN 8 --runMode genomeGenerate --genomeDir [path to store the indexed files] --genomeFastaFiles [path to .fa file of genome] --sjdbGTFfile [path to gtf/gff file of genome] --sjdbGTFtagExonParentTranscript Parent --genomeSAindexNbases 14 --sjdbOverhang 149”

If you’re already in the folder where the genome data files are present then use a dot ‘.’ instead of the entire directory path.

“STAR --runThreadN 8 --runMode genomeGenerate --genomeDir . --genomeFastaFiles ./Setaria_italica.Setaria_italica_v2.0.dna.toplevel.fa --sjdbGTFfile ./Setaria_italica.Setaria_italica_v2.0.59.gff3 --sjdbGTFtagExonParentTranscript Parent --genomeSAindexNbases 13 --sjdbOverhang 149”

Use ——genomeSAindexNbases for Oryza sativa - 13
 Setaria italica - 13
 Eleusine coracana - 14
 Hordeum vulgare - 14
 Magnaporthe oryzae - 11 

In some cases, in the GTF file, the exon regions are not mentioned properly in the annotated files and the exon column is missing in the annotated files. In that case “--sjdbGTFtagExonParentTranscript Parent” command will not work, because this command is designed naturally to catch the exon column in the annotated files. Hence, to perform the indexing with these types of files, you need first check the terminology used instead of exon in the files by opening it. Either it is mentioned to be CDS or transcripts (if only there is no gene transcript column). During this error use the following command: “--sjdbGTFfeatureExon CDS”


And for certain things you’ll need an indexed file of the genome fasta file. Hence, to do that, follow the command “samtools faidx [path to your fasta file]”.



For aligning the rna_seq files to the indexed files using STAR aligner
	For making bam files (alignment files), follow the instructions:
Create a folder to store the BAM files and you can also store the BED files over there as well
Run the command

“STAR --runThreadN 4 --genomeDir [path of genome indexed files] --readFilesCommand gunzip -c --readFilesIn [path to .fq.gz file 1] [path to .fq.gz file 2] --outFileNamePrefix [output file path with file name] --alignIntronMin 20 --alignIntronMax 10000 --outSAMtype BAM SortedByCoordinate”

To align a single fastq file to the genome index of an organism use the below code: “c”
You need to optimise the max intron length during alignment depending on the sample or organism. 20 is the minimum and max should go up to 10000. If you set the max value very less, then you lose some introns. And if you set the value to be much bigger than needed then you’ll find misalignments.
To visualise the aligned file, use IGV software
To get some information about the created .bam file you can use the command “samtools flagstat [bamfilepath] > [path to output file].txt” eg. “samtools flagstat /Users/chandan/Desktop/Analysis_CKPRNA/BAM_BEDfiles/OSC241_10000maxAligned.sortedByCoord.out.bam > /Users/chandan/Desktop/Analysis_CKPRNA/Allignedtext/OSC241_10000maxAligned.txt”
Check properly if you have a python version compatible with openssl version (it’s a tool complementary to flagstat command), otherwise the above command will  not work. If this happens and the env you have created doesn’t contain the required python version (most probably the lower version 3.8), then create another conda env and install the tool with compatible python version from conda forge. “conda create -n sub -c conda-forge python=3.8” and install samtools in this environment to run the above command on flagstat.



For visualising the .BAM file for right alignment
	To check if the aligned files are correct, having proper description of intron lengths, we need to first see the aligned .bam file in IGV software:
Download the IGV software
Create index files of the .bam files first by running this code “samtools index <bamfilepath>” eg. “samtools index /Users/chandan/Desktop/Analysis_CKPRNA/BAM_BEDfiles/F6/CKPRNA1/CKPRNA1_10000maxAligned.sortedByCoord.out.bam”. This will create a small file with .bai file format
Then go to IGV> genomes> load the .fa file of your genome of choice
Go to files> and load the .gff/.gtf file for visualising the gene lengths
Go to files> load the index file of .bam



For separating aligned and unaligned reads from a bam file 
	To separate the aligned and unaligned reads from a .bam file using samtools, use to the below command
 “samtools view -b -f 4 input.bam -o unaligned.bam” or 
“samtools view -b -F 4 input.bam -o aligned.bam”

-f flag with a value of 4. The -f 4 flag selects reads that are not aligned (unmapped) and -F flag with a value of 4. The -F 4 flag filters out reads that are aligned (mapped)


To check if the output bam file contains any read or not use samtools flagstat command or view command.

“samtools flagstat [path to bam file]”

This command will show you the read numbers in terminal

“samtools view [path to bam file] | head”
This will also show the information in the terminal. But if you get a empty return or just the header, then the .bam file is empty



For converting .bam file to .fastq file
	This step can be done using samtools by the following command
“samtools fastq /Users/chandan/Desktop/Analysis_CKPRNA/CKPRNA14_10000max_unalligned.bam > CKPRNA14_10000max_unalligned.fastq”



To read the details of alignment 
	To see the percentage of alignment and the number of reads aligned to the annotated file use the command ‘cat’ eg. “cat [path to Log.final.out file]”
While aligning the fastq files to the genome index it should have created multiple files including .bam file. The file named .final.out contains the information of alignment

For aligning the RNA_seq reads with genome using hisat2
	To use hisat2 to build the aligning file:
Install hisat2 in the conda environment (see conda website for command to install hisat2)
Expand your seq reads (.fq files) bu using “gunzip command” eg. “gunzip /Users/chandan/Desktop/Analysis_CKPRNA/for_hisat2_allign/CKPRNA4-LGB2424-RL1_L1_1.fq.gz”
Build the genome index by the following command. “hisat2-build /Users/chandan/Desktop/Analysis_CKPRNA/for_hisat2_allign/gen2/Magnaporthe_oryzae.MG8.dna.toplevel.fa /Users/chandan/Desktop/Analysis_CKPRNA/for_hisat2_allign/gen2/Mo_index”. The last directory input_index is the name in which you would like to create index files.
To align the sequences with the genome indexed files, use the command “hisat2 -p 32 -x /Users/chandan/Desktop/Analysis_CKPRNA/for_hisat2_allign/gen2/Mo_index -1 /Users/chandan/Desktop/Analysis_CKPRNA/for_hisat2_allign/CKPRNA4-LGB2424-RL1_L1_1.fq /Users/chandan/Desktop/Analysis_CKPRNA/for_hisat2_allign/CKPRNA4-LGB2424-RL1_L1_2.fq --max-intronlen 10000 -S /Users/chandan/Desktop/Analysis_CKPRNA/for_hisat2_allign/CKPRNA4Mo_hisat2_10000max.sam”
To convert and sort the .sam files to .bam files, use the following command “samtools view -bS filename.sam | samtools sort -o filename.sorted.bam”



# Run STAR alignment
STAR --runThreadN ${THREADS} --genomeDir /"${GENOME_DIR}" --readFilesCommand gunzip -c --readFilesIn "path to R1 file" "path to R2 file" --outFileNamePrefix "SAMPLE_OUTPUT_DIR" --alignIntronMin 20 --alignIntronMax 10000 --outSAMtype BAM SortedByCoordinate

    #STAR --runThreadN ${THREADS} \
         #--genomeDir "${GENOME_DIR}" \
         #--readFilesIn "${R1_FILE}" "${R2_FILE}" \
         #--readFilesCommand zcat \
         #--outFileNamePrefix "${SAMPLE_OUTPUT_DIR}/" \
         #--outSAMtype BAM SortedByCoordinate \
         #--quantMode GeneCounts
        
    echo "Alignment completed for ${SAMPLE_NAME}"




For converting .fna files
	While downloading genome files from NCBI, it provides genome sequences in .fna format (.fna format is a rigid format, only able to handle nucleotide sequences. Whereas, .fa file can handle both nucleotides and protein sequences). To convert the .fna file to .fa file, use this command: “mv”.

First cd to the location where you have stored the sequence and then type “mv sequence.fna sequence.fa”

Importantly for some of the sequences, the NCBI website provides two types of files, initiating with GCA_ and GCF_. In these scenarios, always go with GCF_ files. GCA files are those which have been submitted to GeneBank from the submitters as it is. While on the other hand, GCF files are those, which NCBI has cross checked the data and is more accurate than what submitted to GeneBank.



For Intron length distribution curve
	For obtaining a graph to visualise intron length distribution in a genome of an organism, follow the below steps:
Create an environment in conda 
Then install this “conda install -c bioconda misopy”
There are scripts available in both; perl command and python command. Whichever is feasible, run that in the terminal

To run the perl command “perl [path to script] [path to .gff file] > [path to store output file]/file_name.txt”



Learning R programming and studio is important to plot graphs (very important)
Codes to plot a few types of graphs and work with table values in R-studio.



Colours have been picked up by hexadecimal colours (go to the website and copy the #6-digit number and paste it in R-script)



To find the maximum value of a column in a table using R
Command “max(data$MaxValue, na.rm = TRUE)”

To find the position of maximum value of a column in a table using R
Command “which.max(data$MaxValue)”



To run your seqs through SpliSER, follow the commands

For looking at strandedness: craig says if u open the bam in IGV and colour reads by first of pair strand, and then look as a + strand gene (load gff file to see +/-), there are two scenarios:

#open IGV software > open your bam file to see the alignment > right click on the reads > choose “colour alignment by” > choose first-of-pair-strand > check for +ve strands (->) arrows and -ve strands (<-)  arrows.

#if your sequencing is stranded sequenced, then you’ll see both red and blue colour or else a single color.

CASE1. if the reads in the + strand gene are red (+strand) then use regtools -s 2 (RF)
CASE2. if the reads in the + strand gene are blue (-strand) then use regtools -s 1 (FR)

THEN later in spliser you will use the OPPOSITE.

If using CASE1 for beds, use SpliSER -s fr
If using CASE2 for beds, use SpliSER -s rf

If you use the wrong options, then genes will not be assigned and there will be lots of NA gene names. If the data is unstranded (mix of blue and red in IGV) then use -s 0 for bed.




bed file command
“regtools junctions extract -m 20 -M 10000 -o [output path] -s 1 [input path of .bam files]” #the input folder also needs to have the bam indexed file to make .bed files


If running in WSL linux platform then use the code as follows:
regtools junctions extract -m 20 -M 10000 -o /mnt/e/Analysis_CKPRNA/BAM_BEDfiles/BED_files/ECC12R1_EC_10000maxAligned.sortedByCoord.out.bed -s FR /mnt/e/Analysis_CKPRNA/BAM_BEDfiles/bam_stored/ECC12R1_EC_10000maxAligned.sortedByCoord.out.bam



spliser process command
“python [path to SpliSER_v0_1_8.py] process -B [path to bam file] -b [path to bed file] --isStranded -s rf -o [path to output file as SpliSERout.10k] -A [path to gtf/gff file]”



spliser combined
“python [path to SpliSER_v0_1_8.py] combine -S [ path to RNAseqAnalysis/samplefile] --isStranded -s rf -o [path to Combined_Batch1and2_10k_SpliSER.combined.tsv]”



spliser output
“python [path to SpliSER_v0_1_8.py] output -S [path to samplefile] -C [path to Combined_Batch1and2_10k_SpliSER.combined.tsv] -t DiffSpliSER -o [path to Output_Batch1and2_10k_SpliSER]”




For working with subreads for differential gene expression
It uses feature count command to perform the function, and the command is:
“featureCounts -T 8 \
   -a [path to gtf/gff file] \
   -t gene -g ID \
   -o [path to output file and name.txt] \
   -p --countReadPairs\
   [path to bam file]”

If you want to use to do the same function with all the bam files in a directory then mention:
“featureCounts -T 8 -a [path to gtf/gff file] -t gene -g ID -o [path to output file and name.txt] -p --countReadPairs [path to folder of bam files]/*.bam”

If you want to assign featurecounts to count multi-mapping and multi-overlapping reads then use:
-M (for multi-mapping), and 
--fraction (for multi-overlapping)
“featureCounts -T 8 -a [path to gtf/gff file] -t gene -g ID -o [path to output file and name.txt] -p --countReadPairs -M --fraction [path to folder of bam files]/*.bam”

If working in linux platform (or WSL system), then first move to the directory where the code is stored and then run the script like the below example:
(geneanalysis) chandan10@LAPTOP-7J3VU1K9:/mnt/e/Analysis_CKPRNA/codes/ubuntu_codes$ ./subread_count.sh




For checking and separating mismatched reads from from the files of a paired-end sequenced file
In some instances the percentage of unmapped reads (too short, eg, less than 300 in my case as the paired end reads are of 150 base pairs) can go really high because the lines in the corresponding fastq files are not matching and can give this type of output.
Hence, to check if there are any mismatched lines in the fastq files, first install fastq-pair by “conda install bioconda::fastq-pair”

Then run the command to have a look into the corresponding files: “fastq_pair file1.fastq file2.fastq”. This will generate 4 output files, 2 of each having matching reads and 2 files having separated reads from each of the input files that are not matching with the other file. To read the files containing information of lines use the command “wc -l file.single.fq”. If any output number pops up that means that the number of lines in that particular file type doesn’t have a match in the other file.

You can use the files with name paired.fastq as they contain the common reads or matching reads from both the input files.




Creating indexed file for big genome organism (Hordeum vulgare)
When creating indexed file of BAM files of a larger genome organism for example barley, the normal samtools command to create a .bai file will not work and will give an error saying : Region 536989088..536989114 cannot be stored in a bai index. Try using a csi index. Read 'A00599:622:HTNHHDSXC:1:1136:21052:3129' with ref_name='NC_058519.1', ref_length=665585731, flags=163, pos=536989089 cannot be indexed.

Chromosome of contigs up to ~536 million bases can only be handled or stored in .bai file. But if the region is more than the limit, then you have to store the indexed file in the .csi file. To create a .csi file, follow the commands: “samtools index -c [path to .bam file]”. Probably .bed files can be created using .csi files, but IGV software doesn’t work properly with .csi files.

Alternatively, either you can remove the region/contigs/chromosome if possible from the bam files or else you have to convert the BAM files to CRAM files and then make a .bai file out of it.
