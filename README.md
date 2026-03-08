# group1
LIFE750 group 1 poser data
1.	Raw Data

###### using the genomics env #####
###### use comet activate genomics

# Obtaining the data
mkdir 1-Raw
cd 1-Raw
wget -r -np -nH --cut-dirs=3 -A "*.fastq.gz" https://cgr.liv.ac.uk/454/acdarby/LIFE750/

2.	FastQC
#R1 Fastqc
mkdir R1_fastqc
fastqc -t 3 -o R1_fastqc *R1.fastq.gz
#R2 Fastqc
mkdir R2_fastqc
fastqc -t 3 -o R2_fastqc *R2.fastq.gz

3.	MultiQC summary report
#R1 Multiqc report
mkdir R1_fastqc/multiqc
multiqc -o R1_fastqc/multiqc R1_fastqc
#R2 Multiqc report
mkdir R2_fastqc/multiqc
multiqc -o R2_fastqc/multiqc R2_fastqc
#View the reports
firefox R1_fastqc/multiqc/multiqc_report.html \
R2_fastqc/multiqc/multiqc_report.html &

4.	Trimming
# Removing adapters and low quality bases 
# Make new directory for trimmed reads
cd 
mkdir 2-Trimmed
cd 2-Trimmed

#Run trim galore
trim_galore --paired --quality 20 --stringency 4 \
../1-Raw/K1_R1.fastq.gz ../1-Raw/K1_R2.fastq.gz

# Rename the files 
mv K1_R1_val_1.fq.gz K1_R1.fq.gz
mv K1_R2_val_2.fq.gz K1_R2.fq.gz
mv K2_R1_val_1.fq.gz K2_R1.fq.gz
mv K2_R2_val_2.fq.gz K2_R2.fq.gz
mv W1_R1_val_1.fq.gz W1_R1.fq.gz
mv W1_R2_val_2.fq.gz W1_R2.fq.gz

#Further fastqc and multiqc on trimmed reads
mkdir R1_fastqc
fastqc -t 3 -o R1_fastqc *R1.fq.gz
mkdir R1_fastqc/multiqc
multiqc -o R1_fastqc/multiqc R1_fastqc
mkdir R2_fastqc
fastqc -t 3 -o R2_fastqc *R2.fq.gz
mkdir R2_fastqc/multiqc
multiqc -o R2_fastqc/multiqc R2_fastqc

5.	Host removal
# Copy host reference fasta file into the working directory 
mkdir 2-5_Host_removed
cd 2-5_Host_removed
wget -r -np -nH --cut-dirs=3 -A "GRCh38_slice.fasta" https://cgr.liv.ac.uk/454/acdarby/LIFE750/

# Index the reference 
bowtie2-build GRCh38_slice.fasta GRCh38_slice.fasta

# Align to host ref
bowtie2 -x GRCh38_slice.fasta   -1 ../K1_R1.fq.gz -2 ../K1_R2.fq.gz -p 12 2> K1_bowtie2_out.log | samtools view -b -S -h > K1_mapped.bam

# Use samtools to extract reads that did not map to host genome 
samtools fastq -f 4   -1 ../K1_R1.u.fq   -2 ../K1_R2.u.fq   K1_mapped.bam

# Repairing reads – remove unpaired using BBTools
repair.sh \
  in1=../K1_R1.u.fq \
  in2=../K1_R2.u.fq \
  out1=K1_R1.final.fastq \
  out2=K1_R2.final.fastq \
  outs=K1_singletons.fastq

6.	Taxonomic Profiling
# Make new directory for taxonomy 
cd 
mkdir 3-Taxonomy
cd 3-Taxonomy

# Kraken2
# Use the trimmed reads
kraken2 --paired \
  --db  ~/kraken2_db\
  --output K1.kraken \
  --report K1.kreport2 \
  ~/2-Trimmed/K1_R1.fq.gz ~/2-Trimmed/K1_R2.fq.gz
# should produce output file (.kraken) and report file (.kreport2)

7.	Data visualisation
# Krona
ktImportTaxonomy -o kraken2.krona.html *.kreport2
ktImportTaxonomy.sh
firefox kraken2.krona.html &

# Pavian
# In R (type “R” into Ubuntu)
# In R console:
if (!require(remotes)) { install.packages("remotes") }
remotes::install_github("fbreitwieser/pavian")
pavian::runApp(port=5000)
# import .kreport2 files to view

8.	Bracken
bracken -d ~/kraken2_db/\
  -i K1.kreport2 \ # Kraken 2 report file used as input
  -o K1.bracken \ # output bracken file 
  -w K1.breport2 \ # output in kraken format
  -r 100 \ # read length 100bp for classification (based on original library, how many paired end reads)
  -l S \ # taxonomic rank for abundance estimation 
  -t 5 # minimum number of reads needed for classification at specific rank

# Repeat for all other samples
# K2
bracken -d ~/kraken2_db/ \
  -i K2.kreport2 \
  -o K2.bracken \
  -r 100 \
  -l S \
  -t 5

# W1
bracken -d ~/kraken2_db/ \
  -i W1.kreport2 \
  -o W1.bracken \
  -r 100 \
  -l S \
  -t 5

# merging bracken outputs
combine_bracken_outputs.py --files [KW]*.bracken -o all.bracken

# extract abundance columns
seq -s , 4 2 50
bracken_num_columns=$(seq -s , 4 2 50)
echo $bracken_num_columns
cut -f 1,$bracken_num_columns all.bracken > all_num.bracken
less all_num.bracken

9.	Metagenome assembly 

# stitching read pairs using FLASH 
# make new directory 
cd ~
#Make and move into new directory
mkdir 5-Stitched
cd 5-Stitched
#Run flash
flash  -o K1 -z -t 12 -d . \
../2-Trimmed/K1_R1.fq.gz ../2-Trimmed/K1_R2.fq.gz

# FLASH output
ls
less K1.histogram 


# MEGAHIT
# make new directory
cd ..
mkdir 6-Assembly
cd 6-Assembly

# run MEGAHIT using newly stitched read data 
megahit \
-r ../5-Stitched/K1.extendedFrags.fastq.gz \
-1 ../5-Stitched/K1.notCombined_1.fastq.gz \
-2 ../5-Stitched/K1.notCombined_2.fastq.gz \
-o K1 \
-t 12 \
--k-list 29,49,69,89,109,129,149,169,189

# view output
less K1/final.contigs.fa

# QUAST 
# make new directory 
mkdir -p quast/K1

# run QUAST
quast -o quast/K1 K1/final.contigs.fa

# Visualise
firefox quast/K1/report.html

10.	Genome Binning

# MetaBAT2
# make new directory 
mkdir -p ~/7-Binning/K1
#Move into it
cd ~/7-Binning/K1

# index the assembly file prior to alignment 
bwa index ~/6-Assembly/K1/final.contigs.fa

# Alignment of trimmed paired reads 
bwa mem ~/6-Assembly/K1/final.contigs.fa \
~/2-Trimmed/K1_R1.fq.gz ~/2-Trimmed/K1_R2.fq.gz > \
K1.sam

# Sam to sorted bam files 
samtools view -bu K1.sam > K1.bam
# Created sorted bam file
samtools sort K1.bam > K1.sort.bam

# summarise contigs depths 
jgi_summarize_bam_contig_depths --outputDepth K1.depth.txt K1.sort.bam
# view
Less K1.depth.txt

# view summary in R
R
df <- read.table("K1.depth.txt", header=TRUE, check.names=FALSE)
summary(df)

# remove rows with information on contigs shorter than 1500 
df_min1500len <- df[df$contigLen >= 1500,]
summary(df_min1500len)
q()

# running MetaBAT2
#make a diretcory for the bins
mkdir bins
#Run MetaBAT2
metabat2 \
--inFile ~/6-Assembly/K1/final.contigs.fa \
--outFile bins/K1 \
--abdFile K1.depth.txt \
--minContig 1500

# view output
ls bins

# CheckM
# activate checkM  environment
mamba create -n checkm2 -c bioconda -c conda-forge checkm2
mamba activate checkm2
# go into this directory 
cd ~/7-Binning/K1/

# run CheckM 
checkm2 database –download

ls checkm2 predict --threads 20 --input bins -x fa --output-directory checkm_output

# view output
cd ~/7-Binning
mkdir K1_fullset
cd K1_fullset
 https://cgr.liv.ac.uk/454/acdarby/LIFE750/Workshops/7-Binning/K1_fullset

less quality_report.tsv

cat quality_report.tsv | tr "\t" "," > quality_report.csv

# create new quality file 
echo "Quality" > quality_report.csv

# Calculate quality with awk
# extract the completeness and contamination fields/columns
awk -F, '{print $12,$13}' quality_report.csv # enter specific fields 
awk -F, 'NR>1 {print $12,$13}' quality_report.csv
# overall quality calculation and append info to file 
awk -F, 'NR>1 {print $12 - (5 * $13)}' quality_report.csv >> MAGS_quality.csv
less MAGS_quality.csv
#Add quality to the checkm results file 
paste -d "," quality_report.csv MAGS_quality.csv > MAGS_checkm_quality.csv


11.	Assembly functional annotation 
# Bakta
# Make new directory
mkdir ~/8-Annotation
cd ~/8-Annotation

# results
```bash
wget \
  --recursive \
  --no-parent \
  --no-host-directories \
  --cut-dirs=4 \
  https://cgr.liv.ac.uk/454/acdarby/LIFE750/Workshop01/8-Annoataion/

# view summary file
less -S K1.1/K1.1.txt

# view gff file
less K1.1/K1.1.gff3

# search for specific annotations e.g. ATP-binding protein
grep "ATP-binding protein" */*gff3 | less -S

# EC extraction 
Mkdir EC

# create EC annotation files 
ls -1 */*gff3 | sed "s|.*/||" | sed "s/.gff3//" | while read bin
do
cat ${bin}/${bin}.gff3 | grep "EC:" | cut -f 9 | sed "s/ID=//" | \
sed "s/;.*EC:/\t/" | sed "s/,.*//" > EC/${bin}.ec
done

# naming each bin
#List all gff files on one (-1) each
ls -1 */*gff3
#Remove the directory name
#You can use any character as the divider in sed
#Useful when you want to move slashes from a file name
ls -1 */*gff3 | sed "s|.*/||"
#List all the file prefixes (on one line each)
ls -1 */*gff3 | sed "s|.*/||" | sed "s/.gff3//"

# specify input and output files 
ls -1 */*gff3 | sed "s|.*/||" | sed "s/.gff3//" | while read bin
do
echo "${bin}/${bin}.gff3 ../EC/${bin}.ec"
done

cat ${bin}.gff3 | grep "EC:" | cut -f 9 | \
sed "s/ID=//" | sed "s/;.*;Dbxref=/,/" | \
sed "s/,.*EC:/\t/" | sed "s/,.*//" >../EC/${bin}.ec

#Grab every line containing "EC:"
cat K1.1/K1.1.gff3 | grep "EC:" | head

#Cut out the 9th column/field (-f) (i.e. only keep the 9th column)
#This is the attributes field in GFF3
#This contains a plethora of information including the EC annotation if present
#cut uses tabs as the default column/field delimiter
cat K1.1/K1.1.gff3 | grep "EC:" | cut -f 9 | head

#The gff3 attributes field starts with the ID
#We want to keep this but remove the "ID=" part
cat K1.1/K1.1.gff3 | grep "EC:" | cut -f 9 | sed "s/ID=//" | head

#We don't want any of the info between the ID and the EC number
#Therefore we want to remove everything (.*) between
# the first ";" (at the end of the ID info)
# and "EC="
#We'll replace this with a \t to seprarate the ID and EC
# columns with a tab (required by MinPath)
cat K1.1/K1.1.gff3 | grep "EC:" | cut -f 9 |  sed "s/ID=//" | \
sed "s/;.*EC:/\t/" | head

#Finally remove all the info after the EC number
#This info will be after the last ,
cat K1.1/K1.1.gff3 | grep "EC:" | cut -f 9 | sed "s/ID=//" | \
sed "s/;.*EC:/\t/" | sed "s/,.*//" | head

# running minpath
cd ./EC
mkdir ../MetaCyc

# loop through the file suffixes to run Minpath
ls -1 *ec | sed "s/.ec//" | while read bin
do
python /pub14/tea/nsc206/git_installs/Minpath/MinPath.py \
-ec ${bin}.ec \
-report ../MetaCyc/${bin}.minpath
Done

# view output
#Change directory
cd ../MetaCyc
#List contents
Ls
less K1.22.minpath

# KEGG info 
#Change to correct directory
cd ~/8-Annotation
#Make directory for KEGGs
mkdir KEGG
#KEGG extractions
ls -1 */*gff3 | sed "s|.*/||" | sed "s/.gff3//" | while read bin
do
cat ${bin}/${bin}.gff3 | grep "KEGG:" | cut -f 9 | sed "s/ID=//" | \
sed "s/;.*KEGG:/\t/" | sed "s/,.*//" > KEGG/${bin}.kegg
done
#Make directory for KEGG minpath output
mkdir KEGG_minpath
#Change directory to KEGG
cd KEGG
#Run MinPath
ls -1 *kegg | sed "s/.kegg//" | while read bin
do
python /pub14/tea/nsc206/git_installs/Minpath/MinPath.py \
-ko ${bin}.kegg \
-report ../KEGG_minpath/${bin}.minpath
done
<img width="451" height="697" alt="image" src="https://github.com/user-attachments/assets/eb473b42-e3f9-405e-a16b-2f90837f2588" />
