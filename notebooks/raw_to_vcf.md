# Tillandsia whole genome data: from sequencer to VCF
Updated: 28/06/2021

Required: [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), [deML](https://github.com/grenaud/deML), [samtools](https://github.com/samtools/samtools),[bamtools](https://github.com/pezmaster31/bamtools), [bedtools](https://github.com/arq5x/bedtools2),[trim_galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/), [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [picard-tools](https://broadinstitute.github.io/picard/), [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

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

source /home/fs71400/yardeni/.bashr
conda init bash  
conda activate my-env

i=$SLURM_ARRAY_TASK_ID

#map to reference
bowtie2 --very-sensitive-local -x /gpfs/data/fs71400/yardeni/Tillandsia_ref/Tfasc_bowtie2_index -1 "$i"_Lib3_R1_val_1.fq -2 "$i"B_Lib5_R2_val_2.fq -S "$i"_aligned_Tfas.sam -p 16 2> "$i"B_bowtie2.log;
#collect stats on alignment
samtools stats "$ind"_aligned_Tfasc.sam > "$ind"_aligned_Tfasc_sam_stats.txt;
#convert to bam, keep uniquely mapped reads only, filter by mq>10 & sort
samtools view -h -b -q 10 "$ind"_aligned_Tfasc.sam | samtools sort -o "$ind"_asmq10.bam;
#add read group info. This is specific to my data, other users would want to modify
picard AddOrReplaceReadGroups I="$ind"_asmq10.bam o="$ind"_asmq10rg.bam RGLB=WGD RGPL=illumina RGPU=Lib3 RGSM="$ind" RGID="$ind";
#mark duplicates. notice duplicates here are MARKED NOT REMOVED
picard MarkDuplicates I="$ind"_asmq10rg.bam o="$ind"_asmq10rgd.bam M="$ind"_dup_metrics.txt;
```
## per-sample variant calling using HaploTypeCaller

This script was adapted to run on the Vienna Scientific Cluster (VSC) and in addition designed to run 9 tasks on 9 cores, using an array. So is the first value on the array is 70, this will run on numbers 70-79.

```bash
#!/bin/bash                                                                                                                                                                                                    
#                                                                                                                                                                                               
#SBATCH -J GATK                                                                                  
#SBATCH -N 1               
#SBATCH --ntasks-per-node=9         
#SBATCH --ntasks-per-core=1
#SBATCH -p mem_xxx
#SBATCH --array=70,80,90                                                                                                                                                                                                                                                                                                       
source /home/fs71400/yardeni/.bashr
conda init bash  
conda activate gatk-vcf                                                                                                                                                                                                                                                                                                                                                                                          
#go to wd                                                                                                                                                                                                     
cd /gpfs/data/fs71400/yardeni/WGS/vcf;                                                                                                                                                                             
                                                                                                                                                                                                                   
j=$SLURM_ARRAY_TASK_ID     
((j+=9))                                                                                                                                                                                                           
                                                                                                                                                                                                                   
for ((i=$SLURM_ARRAY_TASK_ID; i<j; i++)                                                                                                                                        
do
gatk --java-options "-Xmx16G" HaplotypeCaller -R /gpfs/data/fs71400/yardeni/Tillandsia_ref/tillandsia_fasciculata_assembly.sorted.fasta -I "$i"B_asmq10rgd.bam -O "$i"B_Tfas.g.vcf -ERC GVCF &
done                
wait                                                                                                                                                                                            
```

## Join all per-sample vcf