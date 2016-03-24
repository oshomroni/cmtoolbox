# CMWorkflow
CMWorkflow – this container supports a data descriptor paper explaining the data analysis of chromatin immunoprecipitation (ChIP-seq), methylated DNA immunoprecipitation (MedIP-seq) and RNA (RNA-seq) datasets on brain samples from different tissues and cell-types. The analysis consists of two main steps, involving pre-processing and downstream analyses. The pre-processing involves alignment of reads to the target genome, quality control (QC) and visualization of the alignment data. The second step involves using the aligned files to run downstream analyses, including (1) finding differential histone post-translational modifications (DHPTMs) for ChIP-seq data, (2) predicting cis-regulatory modules (CRMs) using ChIP-seq data, (3) finding differentially methylated regions (DMRs) for MedIP-seq data and (4) finding differentially expressed genes (DEGs) and exons (DEEs) for RNA-seq data.

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
## Tools
### General:
- R scripting front-end version (v3.1.1; 2014-07-10) – executes R scripts that contain functions based in various R packages
- Python (v 2.7.6) – used to run functions written in python

### Pre-processing:
**Note**: unless otherwise stated, the tools below are used for all data types (ChIP-, MedIP- and RNA-seq). For example, sratoolkit is used for all data types, but Bowtie2 is used for ChIP- and MedIP-seq, whereas chequeR is used for ChIP-seq only.
- SRA toolkit (v2.5.0) – manipulates SRA files, specifically to convert them to FASTQ files
- Bowtie2 (v2.0.2) – maps FASTQ reads to reference genome to obtain Sequence Alignment/Map (SAM) files (used for ChIP- and MedIP-seq)
- STAR (v2.3.0) – maps FASTQ reads to reference genome to obtain SAM files  (used for RNA-seq)
- Samtools (v0.1.18) – manipulates SAM files to generate sorted, indexed Binary Alignment/Map (BAM) files
- FastQC (v0.10.1) – performs QC on raw sequencing reads in FASTQ files
- MEDIPS (v1.14.0) – an R package that generates files for data visualization (used for all data types) and QC (used for ChIP- and MedIP-seq)
- chequeR (v0.99.0) – an R package that calculates enrichment statistics (used for ChIP-seq)

### Downstream analyses:
**Note**: the tools below are used for specific downstream analyses, and in brackets it is mentioned which ones.
- MACS2 (v2.0.10.20131028 tag: beta) – calls peaks (used for finding DHPTM)
- DESeq2 (v1.4.5) – an R package that runs differential expression/enrichment analyses (used for finding DHPTMs and DEGs)
-	RFECS (v) – identifies enhancers from chromatin state using random forest algorithm (used for predicting CRMs)
-	RSEG (v0.4.4) – calls broad and dispersed peaks (used for predicting CRMs)
-	MEDIPS (v1.14.0) – an R package that calculates the differential enrichment in methylated regions between different conditions (used for finding DMRs)
-	featureCounts (v1.4.6) – counts the number of reads in all genes (used for finding DEGs)
-	dexseq_count.py (v1.10.8) – a python code within the DEXSeq R package; counts the number of reads in all exons (used for finding DEEs)
-	DEXSeq (v1.10.8) – an R package that runs differential expression analysis for exons (used for finding DEEs)

## Usage:
The following sections contain the specific steps used to (1) generate the pre-processed data and (2) run the downstream analyses.

**Note**: all commands were run on a UNIX machine, using Red Hat 4.1.2-50 operating system.

## Pre-processing
### Steps:
1.	Get the SRA-file (SRA repository, NCBI) and convert it to a single-end FASTQ file (SRA toolkit)
2.	Map FASTQ reads to the mouse mm10 reference genome to obtain a SAM file
 1.	Using 2-mismatch seed alignment (Bowtie2) (ChIP- and MedIP-seq)
 2.	Using gapped alignment (STAR) (RNA-seq)
3.	Generate BAM file from SAM file, then sort and index the BAM file (samtools)
4.	Establish the quality of the reads in the aligned BAM files (FastQC)
5.	Filter data to keep high-quality uniquely- and multi-mapped reads (awk and samtools) (ChIP- and MedIP-seq)
6.	Merge replicate BAM files (samtools)
7.	Convert BAM file into WIG file (MEDIPS), and convert it into a bigWig file (wigToBigWig)
8.	Produce enrichment statistics (ChIP- and MedIP-seq)
 1.	Fragment size, normalized strand cross correlation (NSC) and relative strand cross correlation (RSC) (chequeR) (ChIP-seq)
 2.	Saturation correlation (MEDIPS) (ChIP- and MedIP-seq)
 3. Pearson correlation (MEDIPS) (ChIP- and MedIP-seq)

### Script:
**Note**: before running the script, establish a local directory where you will store your input and output. The directory should contain a sub-directory “suppFiles” for the supplementary material from the Github repository, as well as a sub-directory “inHouseCodes” for the scripts used to run some of the analyses (also in Github repository).

The data we consider in this script is of two types:

1.	Replicate ChIP-seq files, where “CON” refers to context, “01H” refers to the mice being sacrificed 1 hour after training, “CA1” and “NEU” refer to the cells being neuronal from CA1 and “H3K4ME3” refers to the chromatin modification histone 3 lysine tri-methyl
2.  RNA-seq file, where “CON” refers to context, “01H” refers to the mice being sacrificed 1 hour after training and “CA1” refers to the cells being from CA1

**1. Get the SRA-file (SRA repository, NCBI) and convert it to a single-end FASTQ file (SRA toolkit)**

Get the SRA file from SRA archives through FTP

Download SRA file corresponding to sample CON-01H-CA1-NEU-H3K4ME3-1
```
curl ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX143%2FSRX1430120/SRR2922205/SRR2922205.sra > CON-01H-CA1-NEU-H3K4ME3-1.sra
```
Download SRA file corresponding to sample CON-01H-CA1-NEU-H3K4ME3-2
```
curl ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX143%2FSRX1430121/SRR2922206/SRR2922206.sra > CON-01H-CA1-NEU-H3K4ME3-2.sra
```
Download SRA file corresponding to sample CON-01H-CA1-1
```
curl ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX143%2FSRX1430252/SRR2922335/SRR2922335.sra > CON-01H-CA1-1.sra
```
**Output**: downloaded SRA file

Convert SRA to gzipped FASTQ file with single-end reads
```
fastq-dump CON-01H-CA1-NEU-H3K4ME3-1.sra CON-01H-CA1-NEU-H3K4ME3-2.sra CON-01H-CA1-1.sra --gzip
```
**Output**: fastq.gz file containing single-end reads

**2.	Map FASTQ reads to the mouse mm10 reference genome to obtain a SAM file**

create directory for bam files
```
mkdir bam
```
For ChIP- and MedIP-seq data, run Bowtie2

Use Bowtie2 with "-x" for reference genome mm10 (bowtie2 index files available as supplementary files in this repository), "-t" to return runtime of algorithm, "-U" as input FastQ file, "-S" as output SAM file and "-p" as number of threads 

```
bowtie2 -x suppFiles/bowtie2/mm10 -t -U CON-01H-CA1-NEU-H3K4ME3-1.fastq.gz -S bam/CON-01H-CA1-NEU-H3K4ME3-1.sam -p 8
```
For RNA-seq data, run STAR
  
Use STAR with "--outFileNamePrefix" as output prefix of the alignment log files (including directory and file name), "--genomeDir" as reference genome mm10 (STAR index files available as supplementary files in this respository), "--readFilesIn" as input FastQ file, "--runThreadN" as number of threads, "--readFilesCommand zcat" to indicate the input files are zipped, "--outStd SAM", "--outFilterMismatchNmax 2" to allow 2 maximum mismatches per read, "--outSAMstrandField intronMotif" to allow spliced alignments for un-stranded data with XS strand attribute and "--sjdbOverhang 49" as length of the donor/acceptor sequence on each side of the junctions
```
STAR --outFileNamePrefix bam/CON-01H-CA1-1 --genomeDir suppFiles/star/mm10_sr --readFilesIn CON-01H-CA1-1.fastq.gz --runThreadN 8 --readFilesCommand zcat --outStd SAM --outFilterMismatchNmax 2 --outSAMstrandField intronMotif --sjdbOverhang 49 > bam/CON-01H-CA1-1.sam
```
**Output**: SAM file

**3. Generate BAM file from SAM file, then sort and index the BAM file (samtools)**

Generate sorted BAM file without unaligned reads (-F 4) using samtools
```
samtools view -S -h -F 4 -bt bam/CON-01H-CA1-NEU-H3K4ME3-1.sam | samtools sort - bam/CON-01H-CA1-NEU-H3K4ME3-1
```
Generate index for BAM file using samtools
```
samtools index bam/CON-01H-CA1-NEU-H3K4ME3-1.bam
```
**Output**: sorted, indexed BAM file

**4. Establish the quality of the reads in the aligned BAM files (FastQC)**

make output directory
```
mkdir fastqc
```
FastQC analysis, results in the created output directory
```
fastqc bam/CON-01H-CA1-NEU-H3K4ME3-1.bam --outdir fastqc/
```
**Output**: in fastqc directory: a HTML file and a ZIP file, containing the QC report and the images

**5.	Filter data to keep high-quality uniquely- and multi-mapped reads (awk and samtools) (ChIP- and MedIP-seq)**

Create new directory for high-quality BAM files
```
mkdir bam-highquality
```
Create new sorted, indexed BAM files with high-quality reads (reads with MAPQ values 0, 2, 3 and 4 are filtered out)
```
samtools view -h bam/CON-01H-CA1-NEU-H3K4ME3-1.bam | awk '$5 != 0 && $5 != 2 && $5 != 3 && $5 != 4 {print}' > bam-highquality/CON-01H-CA1-NEU-H3K4ME3-1.sam
samtools view -bS bam-highquality/CON-01H-CA1-NEU-H3K4ME3-1.sam | samtools sort - bam-highquality/CON-01H-CA1-NEU-H3K4ME3-1
samtools index bam-highquality/CON-01H-CA1-NEU-H3K4ME3-1.bam
```
**Ouput**: sorted, indexed BAM file with high-quality reads

**6. Merge replicate BAM files (samtools)**

Create new directory for merged BAM files
```
mkdir bam-merged
```
Merge replicate BAM files, where for RNA-seq use regular BAM files and for ChIP- and MedIP-seq use high-quality BAM files. The command is inserted as the path to the merged BAM file, followed by all the replicate files to be merged
```
samtools merge bam-merged/CON-01H-CA1-NEU-H3K4ME3.bam bam-highquality/CON-01H-CA1-NEU-H3K4ME3-1.bam bam-highquality/CON-01H-CA1-NEU-H3K4ME3-2.bam
```
generate index of merged BAM file
```
samtools index bam-merged/CON-01H-CA1-NEU-H3K4ME3.bam
```
**Output**: merged, sorted, indexed BAM file

**7. Convert BAM file into WIG file (MEDIPS), and convert it into a bigWig file (wigToBigWig)**

Create new directory for bigWig files
```
mkdir bw-highquality
```
Create text file containing paths to input BAM files
```
echo -e "bam-highquality/CON-01H-CA1-NEU-H3K4ME3-1.bam\n bam-highquality/CON-01H-CA1-NEU-H3K4ME3-2.bam" > inputFiles.txt
```
Run Rscript to convert BAM files to wig files, using "inputFiles.txt" as a list of BAM files converted to WIG files, "bw-highquality" as the output directory and "2" as the number of cores to run the analysis (when enough cores are available, set number of cores to number of samples to run the analysis as fast as possible)
```
Rscript bam2wig.R inputFiles.txt bw-highquality 2
```
Convert the WIG files to bigWig files, where "CON-01H-CA1-NEU-H3K4ME3-1.50.wig" is the input file, “mm10.chrom.sizes” is a tab-delimited file containing sizes of different mice chromosomes (available in supplementary files in this repository), and "CON-01H-CA1-NEU-H3K4ME3-1.bw" is the output file
```
./suppFiles/wigToBigWig bw-highquality/CON-01H-CA1-NEU-H3K4ME3-1.50.wig suppFiles/mm10.chrom.sizes bw-highquality/CON-01H-CA1-NEU-H3K4ME3-1.bw
```
**Output**: a bigWig file which can be uploaded to our online genome browser (https://memory-epigenome-browser.dzne.de/JBrowse-1.11.4/index.html?data=sample_data/json/mm10)

**8. Produce enrichment statistics (ChIP- and MedIP-seq)**

create new directory for enrichment QC
```
mkdir enrich_QC
```
Run fragment size estimation and NSC and RSC calculations, using "inputFiles.txt" as a list of BAM files tested for their enrichment, "mm10\_50mer\_uniquely\_mappable.bw" as the mappability file used to focus the enrichment analysis on uniquely mapped regions, "enrich\_QC" as the output directory and “2” as the number of cores to run the analysis (when enough cores are available, set number of cores to number of samples to run the analysis as fast as possible)
```
Rscript inHouseCode/EstimateFragSizesChequeR.R inputFiles.txt suppFiles/mm10_50mer_uniquely_mappable.bw enrich_QC 2
```
**Output**: enrichment plots for run samples, as well as a table with RSC and NSC statistics for run samples. These results should correspond to the results presented in Supp. Tab. 3 in the main paper (Halder et. al., 2015)

Calculate saturation correlation for each sample, using "inputFiles.txt" as a list of BAM files tested for their enrichment, "enrich_QC" as the output directory and "2" as the number of cores to run the analysis
```
Rscript inHouseCode/SaturationCorrelation.R inputFiles.txt enrich_QC 2
```
**Output**: saturation correlation plots for run samples, as well as a table with number of reads used to obtain the saturation correlation and the saturation correlation values itself for run samples. These results should correspond to the results presented in Supp. Tab. 3 in the main paper (Halder et. al., 2015)

Calculate Pearson correlation between replicates, using "inputFiles.txt" as a list of BAM files tested for their enrichment, "enrich_QC" as the output directory and "2" as the number of cores to run the analysis
```
Rscript inHouseCode/PairwiseCorrelationCalculation.R inputFiles.txt enrich_QC 2
```
**Output**: a table with the Pearson correlations between all run samples, where the particular interest is in those values between replicate samples. These results should correspond to the results presented in Supp. Tab. 3 in the main paper (Halder et. al., 2015)

## Downstream analyses
### Analyses used:
1.	Finding DHPTMs (ChIP-seq)
2.	Predicting CRMs (ChIP-seq)
3.	Finding DMRs (MedIP-seq)
4.	Finding DEGs (RNA-seq)
5.	Finding DEEs (RNA-seq)

### Finding DHPTMs (ChIP-seq)
#### Steps:
1.	Estimate shift sizes of merged BAM files (chequeR)
2.	Call peaks to identify enriched regions using previously determined fragment sizes (MACS2)
3.	Establish differentially enriched regions (DESeq2)
#### Script:
1.	Estimate shift sizes of merged BAM files (chequeR)

Create file with paths to merged high-quality BAM files that will be used to estimate shift sizes
```
echo -e "BAM-HIGHQUALITY/NAI-01H-CA1-NEU-H3K4ME1.bam\nBAM-HIGHQUALITY/SHC-01H-CA1-NEU-H3K4ME1.bam" > inputFiles.txt
```
Run fragment size estimation on merged BAM file, using "inputFiles.txt" as a list of BAM files the shift sizes are estimated for, "mm10\_50mer\_uniquely\_mappable.bw" as the mappability file used to focus the shift size estimation on uniquely mapped regions, "." as the output directory and "2" as the number of cores to run the analysis
```
Rscript inHouseCodes/EstimateFragSizesChequeR.R inputFiles.txt suppFiles/mm10_50mer_uniquely_mappable.bw . 2
```
2.	Call peaks to identify enriched regions using previously determined fragment sizes (MACS2)
create peak directory
```
mkdir chip_peaks
```
Call narrow peaks for the merged BAM files, using mouse genome (-g mm), without model creation (--nomodel), shift size used equal to the previously found fragment size divided by 2 (--shiftsize fragment_size/2) and minimum false-discovery rate (FDR) cutoff for peak detection at 0.05 (-q 0.05)
```
fragSize1= fragment size for NAI-01H-CA1-NEU-H3K4ME1.bam
fragSize2=fragment size for SHC-01H-CA1-NEU-H3K4ME1.bam 
macs2 callpeak -t BAM-HIGHQUALITY-MERGED/NAI-01H-CA1-NEU-H3K4ME1.bam -f BAM -g mm -n chip_peaks/NAI-01H-CA1-NEU-H3K4ME1.narrow --nomodel --shiftsize fragSize1/2 -q 0.05
macs2 callpeak -t BAM-HIGHQUALITY-MERGED/SHC-01H-CA1-NEU-H3K4ME1.bam -f BAM -g mm -n chip_peaks/SHC-01H-CA1-NEU-H3K4ME1.narrow --nomodel --shiftsize fragSize2/2 -q 0.05
```
3.	Establish differentially enriched regions (DESeq2)
Run DE analysis for geneBody, where the DHPTMs are determined for entire genes from TSS (transcript start site) to TTS (transcript termination site)

Generate text file with high-quality BAM files for gene-body DE analysis
```
echo -e "BAM-HIGHQUALITY/NAI-01H-CA1-NEU-H3K79ME3-1.bam\nBAM-HIGHQUALITY/NAI-01H-CA1-NEU-H3K79ME3-2.bam\nBAM-HIGHQUALITY/SHC-01H-CA1-NEU-H3K79ME3-1.bam\nBAM-HIGHQUALITY/SHC-01H-CA1-NEU-H3K79ME3-2.bam" > geneBody.txt
```
Create output directory
```
mkdir geneBodyDE
```
Run DE analysis for genes, using "geneBody.txt" as a list of BAM files used to run the DE analysis (first two files are control, second two files are treatment; there are always 2 replicates for the ChIP-seq data), "geneBodyDE" as the output directory, "1" as the number of as the number of cores to run the analysis and "0.1" and "0" as the adjusted p-value and log2 fold change thresholds, respectively, used to determine significant genes
```
Rscript inHouseCodes/finalGeneBodyDE.R geneBody.txt geneBodyDE 1 0.1 0
```
Run DE analysis for merged peak files, where the DHPTMs are determined for previously established peak regions

Generate text file with high-quality BAM files for peak DE analysis
```
echo -e "BAM-HIGHQUALITY/NAI-01H-CA1-NEU-H3K4ME1-1.bam\nBAM-HIGHQUALITY/NAI-01H-CA1-NEU-H3K4ME1-2.bam\nBAM-HIGHQUALITY/SHC-01H-CA1-NEU-H3K4ME1-1.bam\nBAM-HIGHQUALITY/SHC-01H-CA1-NEU-H3K4ME1-2.bam" > peaks.txt
```
Create output directory
```
mkdir peaksDE
```
Run DE analysis for peaks, using "peaks.txt" as a list of BAM files used to run the DE analysis, "chip\_peaks" as the directory containing the peak files used for the analysis, "peaksDE" as the output directory, "1" as the number of as the number of cores to run the analysis and "0.1" and "0" as the adjusted p-value and log2 fold change thresholds, respectively, used to determine significant genes
```
Rscript inHouseCodes/finalPeaksDE.R peaks.txt chip_peaks peaksDE 1 0.1 0
```
Run DE analysis for TSS regions, where the DHPTMs are determined for particular regions around the TSS (in our case, we use [-500,+1000])

Generate text file with high-quality BAM files for TSS DE analysis

```
echo -e "BAM-HIGHQUALITY/NAI-01H-CA1-NEU-H3K4ME3-1.bam\nBAM-HIGHQUALITY/NAI-01H-CA1-NEU-H3K4ME3-2.bam\nBAM-HIGHQUALITY/SHC-01H-CA1-NEU-H3K4ME3-1.bam\nBAM-HIGHQUALITY/SHC-01H-CA1-NEU-H3K4ME3-2.bam" > tss.txt
```
Create output directory
```
mkdir tssDE
```
Run DE analysis for TSS regions, using “tss.txt” as a list of BAM files used to run the DE analysis, “tssDE” as the output directory, “1” as the number of as the number of cores to run the analysis, “0.1” and “0” as the adjusted p-value and log2 fold change thresholds, respectively, used to determine significant genes, and “500” and “1000” as the upstream and downstream limits around the TSS
```
Rscript finalTssDE.R tss.txt tssDE 1 0.1 0 500 1000
```
**Output**: raw tables containing complete statistics of DE analyses of any of the regions (gene body, TSS regions or peaks), where one table contains results for all regions, and another for regions with average read count higher than or equal to 2\*region\_width/50 for either control or treatment samples. In addition, the code would output a table showing the overall number of DHPTMs for each analysis under regular (all regions) or filtered (read count >= 2*region_width/50) conditions, equivalent to Supp. Table 8 in the main paper (Halder et. al., 2015)

### Prediction CRMs (ChIP-seq)
#### Steps:
1.	Create positive CRM-predicting set (bedtools)
2.	Create negative CRM-predicting set (bedtools)
3.	Train sets for CRM prediction (RFECS)

### Finding DMRs (MedIP-seq)
#### Steps:
1.	Identify differentially methylated regions (DMRs) (MEDIPS)
Script:
1.	Identify differentially methylated regions (DMRs) (MEDIPS)
create directory for DMR results
```
mkdir dmr
```
Get DMR results using medipDE.R script
```
Rscript medipDE.R DMRs CON-01H-ACC-NEU-MC-1.bam,CON-01H-ACC-NEU-MC-2.bam SHC-01H-ACC-NEU-MC-1.bam,SHC-01H-ACC-NEU-MC-2.bam TRUE 250 0 BSgenome.Mmusculus.UCSC.mm10 700 5 dmrFile BAM-HIGHQUALITY
```
### Finding DEGs (RNA-seq)
#### Steps:
1.	Count number of reads in genes (featureCounts)
2.	Identify differentially expressed genes (DESeq)
Script:
1.	Count number of reads in genes (featureCounts)
Finding DEEs (RNA-seq)
Steps:
1.	Count number of reads in exons (dexseq_count.py)
2.	Identify differentially expressed genes (DEXSeq)
Script:
1.	Count number of reads in exons (dexseq_count.py)
