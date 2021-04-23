# Tillandsia whole genome data: from sequencer to VCF
Updated: 9/3/2021

Required: [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), [deML](https://github.com/grenaud/deML), [samtools](https://github.com/samtools/samtools),[bamtools](https://github.com/pezmaster31/bamtools),[bedtools](https://github.com/arq5x/bedtools2),[trim_galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

Process raw whole-genome-sequencing data to variant call file. Quality control, pre-processing and variant calling.

Overall, the data was sequenced in over 4 different lanes. Here I'm presenting an example of processing a single lane (called Lib3) but the same applies to all lanes. 

## List of metadata files
- P5P7_for_demulti - tab delimited file with sample name and indices used. See description in the [deML format examples](https://github.com/grenaud/deML)
- ind_list.txt - list of sample names

## demultiplex

VBCF gave us the raw library in BAM format so the code fits that. --bamtags tells deML where the index tags are, and this may depend on your data. Output is same file, with sample name indicated in RG field. 

```bash
> deML -i P5P7_for_demulti --bamtags BC,QT,B2,Q2 -o Lib3-demulti_deML.bam -s demult_stats.txt -e demult_unassigned.txt Lib3_raw.bam
```
Seperating the demultiplexed file to many BAM files, converting to FASTQ & generating a report while at it. Previously used samtools for this step but no longer do due to [this behaviour of adding unassigned reads](https://github.com/samtools/samtools/issues/896) 

```bash
while read ind; do
	bamtools filter -tag RG:"$ind" -in Lib3-demulti_deML.bam -out "$ind"_Lib3.bam;
	samtools stats "$ind"_Lib3.bam > "$ind"_raw_stats.txt;
	bedtools bamtofastq -i "$ind"_Lib3.bam -fq "$ind"_Lib3_R1.fq -fq2 "$ind"_Lib3_R2.fq;
done < ind_list.txt
```

## trim & produce fastqc reports 
_
Compare the raw fastqc reports with the trimmed. This script uses 8 cores but use whatever is available to you.
the output from trim_galore is, by default, files named "$ind"\_Lib3_R1_val_1.fq, "$ind"\_Lib3_R2_val_2.fq for paired end data.

```bash
while read ind; do
	fastqc -t 8 -f fastqc "$ind"_Lib3_R*.fq 
	trim_galore --fastqc --cores 8 --retain_unpaired --paired "$ind"_Lib3_R1.fq "$ind"_Lib3_R2.fq;
done < ind_list.txt
```

## align to reference and prepare for variant call

Prepare your reference by indexing with bowtie2. Ours is a *Tillandsia fasciculata* reference in fasta format:

```bash
bowtie2-build tillandsia_fasciculata_assembly.sorted.fasta Tfasc_bowtie2_index
```


This script was adapted to run on as an array the Vienna Scientific Cluster (VSC) and will probably need to be adapted to run on other system.
I love to keep the bowtie output as a log, hence the 2> redirect, which can be removed.

```bash
#!/bin/bash
#
#SBATCH -J bowtie2
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks-per-core=1
#SBATCH -p mem_xxx
#SBATCH --array=74,77-79,80-99

conda activate my-env

i=$SLURM_ARRAY_TASK_ID

bowtie2 --very-sensitive-local -x /gpfs/data/fs71400/yardeni/Tillandsia_ref/Tfasc_bowtie2_index -1 "$i"_Lib3_R1_val_1.fq -2 "$i"B_Lib5_R2_val_2.fq -S "$i"_aligned_Tfas.sam -p 16 2> "$i"B_bowtie2.log;

```



