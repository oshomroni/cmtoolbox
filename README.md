# CMToolBox
CMToolBox – container to support chromatin modification data analysis for the data descriptor paper describing ChIP-, MedIP and RNA-seq datasets on brain samples from different tissues and cell-types.

## References:
**Network analysis of epigenomic and transcriptomic tissue- and cell-type-specific data from learning mice**

<Author_list>

Scientific Data, <publication_date>

doi: <DOI_ACCESSION>

<PATH_TO_PUBLICATION>

**Methyl-CpG-binding domain sequencing reveals a prognostic methylation signature in neuroblastoma**

Rashi Halder, Magali Hennion, Ramon O Vidal, Orr Shomroni, Raza-Ur Rahman, Ashish Rajput, Tonatiuh Pena Centeno, Frauke van Bebber, Vincenzo Capece, Julio C Garcia Vizcaino, Anna-Lena Schuetz, Susanne Burkhardt, Eva Benito, Magdalena Navarro Sala, Sanaz Bahari Javan, Christian Haass, Bettina Schmid, Andre Fischer & Stefan Bonn

Nature, 2015

doi:10.1038/nn.4194

http://www.nature.com/neuro/journal/v19/n1/full/nn.4194.html#supplementary-information
## Tools:
- sratoolkit (v2.3.5) – NCBI toolkit to work with files from SRA (e.g. convert SRA to FASTQ)
- fastqc (v0.10.1) – QC on raw sequencing reads (FASTQ files)
- bowtie2 (v2.0.2) – maps reads to the reference genome (mm10 index files included as well) (used for ChIP-seq and MedIP-seq data)
- STAR (v2.3.0e_r291) – maps reads to the reference genome (mm10 index files included as well) (used for RNA-seq data)
- samtools (v0.1.18) – tools to manipulate SAM files (used to sort/index BAM files)
- MEDIPS (v1.16.0) – an R package used for analyzing data derived from methylated DNA immunoprecipitation
- chequeR (v0.99) – an in-house R package used to calculate enrichment statistics for ChIP-seq data

## Usage:
Since this toolbox works with different kinds of data, mainly deep sequencing of chromatin immunoprecipitation (ChIP-seq), methylated DNA immunoprecipitation (MedIP-seq) and RNA (RNA-seq), the steps differ slightly for the different data types. The analyses can be divided into two main parts, being “pre-processing” and “differential expression”. The pre-processing is fairly similar between the different data types, involving alignment, quality control (QC) and visualization of the data.
## Pre-processing
### Steps:
1. Get the SRA-file (SRA repository, NCBI) and convert it to FASTQ file (single-end) (SRAtoolkit)
2. Map reads to the mouse mm10 reference genome
  1. Using 2-mismatch seed alignment (bowtie2) (ChIP- and MedIP-seq)
  2. Using gapped alignment (STAR) (RNA-seq)
3. Sort and index the SAM/BAM files (samtools)
4. Establish the quality of the reads in the aligned BAM files (FastQC)
5. Filter data to keep high-quality uniquely- and multi-mapped reads (awk) (ChIP- and MedIP-seq)
6. Merge replicate BAM files (samtools)
7. Convert BAM files into WIG files (MEDIPS), and converting those into bigwig files (wigToBigWig)
8. Produce enrichment statistics (ChIP- and MedIP-seq only)
  1. Fragment size, normalized strand cross correlation (NSC) and relative strand cross correlation (RSC) (chequeR) (ChIP-seq)
  2. Saturation correlation (MEDIPS) (ChIP- and MedIP-seq)
  3. Pearson correlation (MEDIPS) (ChIP- and MedIP-seq)

### Script:
**1. Get the SRA-file (SRA repository, NCBI) and convert it to FASTQ files (single-end) (SRAtoolkit)**

**get the SRA file from SRA archives through FTP**

H3K4me3 file: CON-01H-CA1-NEU-H3K4ME3-1
```
curl ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX143%2FSRX1430120/SRR2922205/SRR2922205.sra > CON-01H-CA1-NEU-H3K4ME3-1.sra
```
H3K4me3 file: CON-01H-CA1-NEU-H3K4ME3-2
```
curl ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX143%2FSRX1430121/SRR2922206/SRR2922206.sra > CON-01H-CA1-NEU-H3K4ME3-2.sra
```
RNA file:
```
curl ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX143%2FSRX1430252/SRR2922335/SRR2922335.sra > CON-01H-CA1-1.sra
```
output: downloaded SRA file

SRA to gzipped FASTQ file with single-end reads
```
fastq-dump CON-01H-CA1-NEU-H3K4ME3-1.sra CON-01H-CA1-NEU-H3K4ME3-2.sra CON-01H-CA1-1.sra --gzip
```
output: 4 fastq.gz files containing single-end reads, corresponding to ChIP- (Input and H3K4me3), MedIP- and RNA-seq data, respectively
**2. Map reads to the mouse mm10 reference genome**
create directory for bam files
```
mkdir bam
```
  **1. For ChIP- and MedIP-seq data, use bowtie2 with 8 cores (-p 8)**

read mapping using reference genome mm10 (-x bowtie2/mm10; bowtie2 index files available in the repository)
```
bowtie2 -x bowtie2/mm10 -t -U CON-01H-CA1-NEU-H3K4ME3-1.fastq.gz -S bam/CON-01H-CA1-NEU-H3K4ME3-1.sam -p 8
```
  **2. For RNA-seq data, use STAR with 8 cores (--runThreadN 8), 2 maximum mismatches per reads (--outFilterMismatchNmax 2), spliced alignments for un-stranded data with XS strand attribute (--outSAMstrandField intronMotif) and length of the donor/acceptor sequence on each side of the junctions set to read_length-1 (--sjdbOverhang 49)**
  
read mapping using reference genome mm10 (--genomeDir mm10; STAR index files available in the repository)
```
STAR --outFileNamePrefix bam/CON-01H-CA1-1 --genomeDir star/sr --readFilesIn CON-01H-CA1-1.fastq.gz --runThreadN 8 --readFilesCommand zcat --outStd SAM --outFilterMismatchNmax 2 --outSAMstrandField intronMotif --sjdbOverhang 49 > bam/CON-01H-CA1-1.sam
```
output: SAM file

**3. Sort and index the SAM/BAM files (samtools)**

Generate sorted BAM file using samtools
```
samtools view -S -h -F 4 -bt bam/CON-01H-CA1-NEU-H3K4ME3-1.sam | samtools sort - bam/CON-01H-CA1-NEU-H3K4ME3-1
```
generate index for BAM file using samtools
```
samtools index bam/CON-01H-CA1-NEU-H3K4ME3-1.bam
```
output: sorted, indexed BAM file

if you want to, you could remove the SRA, FastQ and SAM files
```
rm *.sra *.fastq.gz bam/*.sam
```
**4. Establish the quality of the reads in the aligned BAM files (FastQC)**

make output directory
```
mkdir fastqc
```
FastQC analysis, results in the created output directory
```
fastqc bam/CON-01H-CA1-NEU-H3K4ME3-1.bam --outdir fastqc/
```
output: in fastqc directory: a HTML file and a ZIP file, containing the QC report and the images

**5. Filter data to keep high-quality uniquely- and multi-mapped reads (awk) (ChIP- and MedIP-seq)**

create new directory for high-quality BAM files
```
mkdir bam-highquality
```
create new sorted, indexed BAM files with high-quality reads (MAPQ values 0, 2, 3 and 4)
```
samtools view -h bam/CON-01H-CA1-NEU-H3K4ME3-1.bam | awk '$5 != 0 && $5 != 2 && $5 != 3 && $5 != 4 {print}' > bam-highquality/CON-01H-CA1-NEU-H3K4ME3-1.sam
samtools view -bS bam-highquality/CON-01H-CA1-NEU-H3K4ME3-1.sam | samtools sort - bam-highquality/CON-01H-CA1-NEU-H3K4ME3-1
samtools index bam-highquality/CON-01H-CA1-NEU-H3K4ME3-1.bam
```
output: sorted, indexed BAM file with high-quality reads

**6. Merge replicate BAM files (samtools)**

create new directory for merged BAM files (for RNA-seq use regular BAM files, for ChIP- and MedIP-seq use high-quality BAM files)
```
mkdir bam-merged-highquality
```
merge replicate BAM files
```
samtools merge bam-highquality/CON-01H-CA1-NEU-H3K4ME3.bam bam-highquality/CON-01H-CA1-NEU-H3K4ME3-1.bam bam-highquality/CON-01H-CA1-NEU-H3K4ME3-2.bam
```
generate index of merged BAM file
```
samtools index bam-highquality/CON-01H-CA1-NEU-H3K4ME3.bam
```
output: merged sorted, indexed BAM file

**7. Convert BAM files into WIG files (MEDIPS, wigToBigWig)**

Create new directory for bigWig files
```
mkdir bw-highquality
```
Create text file containing all BAM files that will be converted into bigWig files
```
echo -e "bam-highquality/CON-01H-CA1-NEU-H3K4ME3-1.bam \n bam-highquality/CON-01H-CA1-NEU-H3K4ME3-2.bam" > inputFiles.txt
```
Run Rscript to convert bam files to wig files
```
Rscript bam2wig.R inputFiles.txt bw-highquality 2
```
Get the wigToBigWig script 
```
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
```
Get chromosome sizes for mouse
```
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes
```
Convert the WIG files to bigWig files
```
./wigToBigWig bw-highquality/CON-01H-CA1-NEU-H3K4ME3-1.50.wig mm10.chrom.sizes bw-highquality/CON-01H-CA1-NEU-H3K4ME3-1.bw
```
**8. Produce enrichment statistics (ChIP- and MedIP-seq only)**

create new directory for enrichment QC
```
mkdir enrich_QC
```
Run fragment size estimation and NSC and RSC calculations
```
Rscript EstimateFragSizesChequeR.R inputFiles.txt mm10_50mer_uniquely_mappable.bw enrich_QC 2
```
calculate saturation correlation for each sample
```
Rscript SaturationCorrelation.R inputFiles.txt enrich_QC 2
```
calculate Pearson correlation between replicates
```
Rscript PairwiseCorrelationCalculation.R inputFiles.txt enrich_QC 2
```
