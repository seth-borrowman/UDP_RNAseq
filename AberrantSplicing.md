# Detecting Aberrant Splicing in RNA data

## Purpose
This document should aid a user in detecing aberrant splicing in RNA-seq data. The user should have followed the previous guide (README.md) to retreive RNAseq data from IIHG and should have some familiarity with Argon HPC.

## Create an index
To begin, we must create an index to align RNAseq data to. If you have used Salmon to align the RNAseq data to a transcriptome, the process is similar, but distinct. The Salmon index cannot be used, and we will create a new index and alignment using Star.
### Download necessary files
You will first need to download a fasta (or fa) and GTF file of the human genome. There are multiple sources for this; I've used the GENCODE primary assembly found [here](https://www.gencodegenes.org/human/).
This can be downloaded directly through the HPC using the `curl` command.
Star requires that the files not be compressed. It's likely they were downloaded in .gz format. This can be decompressed with `gzip -d filename`.
### qlogin
I've found that the memory required for the next operations is too much to be run normally. I've had multiple qsub jobs crash while running, but have found luck with qlogin. Before beginning to run STAR, simply enter the command `qlogin`.
### Star indexing
Documentation for Star can be found [here](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf). The command can either be entered into the command line, or be submitted as a job. I'll provide methods for how I used the command line.
1. Navigate to where the fasta and GTF files are and make a new directory (`mkdir`).
2. Load Star by entering `module load star`.
3. Create the index with the following code (changed as needed):
```
STAR --runMode genomeGenerate \
--runThreadN 8 \
--genomeDir StarIndex \
--genomeFastaFiles GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile gencode.v44.primary_assembly.annotation.gtf \
--sjdbOverhang 99 \
--genomeSAsparseD 2
```
- runThreadN is the number of parallel threads to use.
- genomeDir is whichever directory you want the output files to be saved in.
- genomeFastaFiles is the name/location of the fasta file you've downloaded.
- sjdbGTFfile is the name/location of the GTF file you've downloaded.
- sjdbOverhang is the read length from RNAseq - 1. Assume 99 is correct if you do not know the read length.

You will know the indexing has completed successfully when it outputs "finished successfully". For me, this took ~30 minutes. Depending on the RAM used and number of threads available (runThreadN), it is possible it could take much longer.

## Align your reads to the index
Once you have the index made, you need to align your RNAseq reads to it to create a BAM file. The process to do this is very similar to creating the index.
1. Use `qlogin` to start a session.
2. Navigate to where you want the BAM files to be output.
3. Load Star with `module load star`.
4. Generate BAM files using the following code:
```
STAR --runThreadN 8 \
--genomeDir /Shared/lss_chndr/UDP_Research/RNAseq/Indices/StarIndex \
--readFilesIn UIRDB20230003_R1.fastq.gz UIRDB20230003_R2.fastq.gz \
--readFilesCommand "gzip -d" \
--outSAMattrRGline ID:UIRDB20230003 \
--outSAMtype BAM SortedByCoordinate
```
- genomeDir is the path to the directory where you stored your index files previously.
- readFilesIn Is the paths/names of the fastq files to read in. In this example, we have two paired-end reads.
- readFilesCommand is necessary if the fastq files are compressed. It specifies the command Star should use to decompress the files.
- outSAMattrRGline is the IDs for each the sample.
- outSAMtype specifies the output file type. For future steps, we need a BAM file. I've specified here for the BAM file to be sorted.

## Create BAM index files
Many tools also require an individual BAM index file (.bai) for each BAM file. This can be done very easily in Argon using samtools.
1. `module load samtools`
2. `parallel  samtools index ::: *.bam`
This will create a .bai file for each BAM file in your current folder.
