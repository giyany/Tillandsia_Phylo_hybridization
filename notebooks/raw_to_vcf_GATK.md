# Tillandsia whole genome data: from sequencer to VCF
Updated: 21/10/2021

Required: [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), [deML](https://github.com/grenaud/deML), [samtools](https://github.com/samtools/samtools),[bamtools](https://github.com/pezmaster31/bamtools), [bedtools](https://github.com/arq5x/bedtools2),[trim_galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/), [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [picard-tools](https://broadinstitute.github.io/picard/), [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller), [bcftools](https://samtools.github.io/bcftools/bcftools.html), [vcftools](http://vcftools.sourceforge.net/), [SnpSift](http://pcingola.github.io/SnpEff/)

#Overview: 

Process raw whole-genome-sequencing data to variant call file. Quality control, pre-processing, variant calling and filtering.

The data was sequenced in over 4 different lanes. Here I'm presenting an example of processing a single lane (called Lib3) but the same applies to all lanes. 

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

## per-sample variant calling using HaploTypeCaller, includes non-variant sites with --output-mode.

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
gatk --java-options "-Xmx16G" HaplotypeCaller -R /gpfs/data/fs71400/yardeni/Tillandsia_ref/tillandsia_fasciculata_assembly.sorted.fasta -I "$i"B_asmq10rgd.bam -O "$i"B_Tfas.g.vcf -ERC GVCF --output-mode EMIT_ALL_CONFIDENT_SITES &
done                
wait                                                                                                                                                                                            
```

## Join all per-sample vcf

Using GenotypeGVCFs to join all per-sample vcf into one file. For the process to finish within a mortal lifetime, I generated a databse using GenomicsDBimport. First, generate a sample map:

```bash

for file in *.g.vcf;
do echo ${file/.g.vcf/} $file | sed 's/ /\t/g' >> sample_map;
done
```

To generate a bed regions file according to coverage using frebayes' cov_to_regions.py.

```bash
#!/bin/bash
#
#SBATCH -J freebayes_cov_to_reg
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks-per-core=1 
#SBATCH -p mem_xxx
#SBATCH --qos p71400_xxx                                                                                                                                                                     
                                                                                                                                                                                                                   
source /home/fs71400yardeni/.bashr
conda init bash   
conda activate freebyes-env
                                                                                                                                                                                                                   
## move to working directory                                                                                                                                                                                       
                                                                                                                                                                                                                   
cd /gpfs/data/fs71400/yardeni/Tillandsia_ref/

bamtools coverage -in ../WGS/mapped_bams/19B_asmq10rgd.bam | coverage_to_regions.py tillandsia_fasciculata_assembly.sorted.fasta.fai 10000 > Tfas.fa.10k.regions.byCoverage 
```

then separated by chr and corrected fields to bed format:

```bash
while read chr; do cat Tfas.fa.10k.regions.byCoverage | grep "$chr:" >> Tfas.fa.10k.regions.byCoverage_"$chr"; done < Tfas_scaffold_names.txt;
for file in Tfas.fa.10k.regions.byCoverage_*; do sed -e 's/\:/\'$'\t/g' -e 's/\-/\'$'\t/g' "$file" > "$file".bed; done
```

Perform joint calling and include non-variant sites. The array contains all the Chr numbers in my reference:

```bash

#!/bin/bash
#
#SBATCH -J GATK4_JointHC
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH -p mem_xxx
#SBATCH --qos=p71400_xxx
#SBATCH --array=1-2321%4

source /home/fs71400/yardeni/.bashrc                                                                                     
conda init bash                                                                                                   
conda activate gatk-vcf

# go to working directory

cd /gpfs/data/fs71400/yardeni/WGS/vcf

# Run GenomicsDBImport

$gatk --java-options "-Xmx32g -Xms4g" GenomicsDBImport \
      --genomicsdb-workspace-path GenomicsDB_${SLURM_ARRAY_TASK_ID} \
      --sample-name-map sample_map \
      --tmp-dir /gpfs/data/fs71400/yardeni/WGS/vcf/tmp/ \
      -L /gpfs/data/fs71400/yardeni/Tillandsia_ref/Tfas_cov_Scaffolds/Tfas.500.regions.byCov.Scaffold_Scaffold_${SLURM_ARRAY_TASK_ID}.bed \
      --reader-threads 8

# For each database directory make a vcf file using GenotypeGVCfs and concatenate all vcf files into one using bcftools

$gatk --java-options "-Xmx8g" GenotypeGVCFs \
      -R /gpfs/data/fs71400/yardeni/Tillandsia_ref/tillandsia_fasciculata_assembly.sorted.fasta  \
      -V gendb://GenomicsDB_${SLURM_ARRAY_TASK_ID} \
      -O PerChr.Variants.AllTillandsias.WGS.${SLURM_ARRAY_TASK_ID}.vcf.gz \
      --tmp-dir /gpfs/data/fs71400/yardeni/WGS/vcf/tmp/

```
concatenate all vcf files into one vcf and generate some stats

```bash
conda activate my-env                                                                                                                                                                                             
                                                                                                                                                                                                                   
bcftools concat $(for file in *.vcf.gz; do echo "$file "; done) > Variants.AllTillandsias.allChr.WGS.raw.full.vcf                                                                                                  
                                                                                                                                                                                                                   
#generate stats                                                                                                                                                                                                    
                                                                                                                                                                                                                   
bcftools stats Variants.AllTillandsias.allChr.WGS.raw.full.vcf > Variants.AllTillandsias.WGS.raw.full.vcf_stats.txt                                                                                                
SnpSift tstv Variants.AllTillandsias.allChr.WGS.raw.full.vcf > Variants.AllTillandsias.WGS.raw.full.snpsift_stats.txt
vcftools 
```


## Filtering (basic)

```bash

#set DP=0 to missing (./.) then use bcftools to allow max 20% missing data, exclude indels and 3bp around indels. The missing data filtering is not final, but only makes the file more managable.

echo "make file manageable, filter around SNPs, missing 20% allowed"
vcftools --gzvcf Variants.AllTillandsias.allChr.WGS.raw.full.vcf.gz --minDP 1 --recode --recode-INFO-all -c | bcftools filter --SnpGap 3 -i 'F_MISSING<0.2' | bcftools view --exclude-types indels -e 'ALT="*"' -Oz > Variants.AllTillandsias.allChr.WGS.3bpindel.0.2missing.vcf.gz
tabix Variants.AllTillandsias.allChr.WGS.3bpindel.0.2missing.vcf.gz;
echo "generate stats"
bcftools stats Variants.AllTillandsias.allChr.WGS.3bpindel.0.2missing.vcf.gz > Variants.AllTillandsias.allChr.WGS.3bpindel.0.2missing.vcf_stats.txt
SnpSift tstv Variants.AllTillandsias.allChr.WGS.3bpindel.0.2missing.vcf.gz > Variants.AllTillandsias.allChr.WGS.3bpindel.0.2missing.vcf.gz_snpsift.txt

#removing known TEs
echo "remove TEs"
bedtools intersect -v -header -a Variants.AllTillandsias.allChr.WGS.3bpindel.0.2missing.vcf.gz -b /gpfs/data/fs71400/yardeni/Tillandsia_ref//tillandsia_fasciculata_assembly.sorted.25_scaffolds.fasta.mod.EDTA.intact.no_unknown.gff3 | bgzip > Variants.AllTillandsias.allChr.WGS.3bpindel.0.2missing.No_TEs_EDTA.vcf.gz
tabix Variants.AllTillandsias.allChr.WGS.3bpindel.0.2missing.No_TEs_EDTA.vcf.gz

conda activate gatk-vcf
echo "filter for quality"
gatk VariantFiltration \
-R /gpfs/data/fs71400/yardeni/Tillandsia_ref/tillandsia_fasciculata_assembly.sorted.fasta \
-V Variants.AllTillandsias.allChr.WGS.3bpindel.0.2missing.No_TEs_EDTA.vcf.gz \
--filter-name "mq" \
--filter-expression "MQ < 20.00" \
--genotype-filter-expression "DP < 4" \
--genotype-filter-name "genotypedepth" \
--filter-name "QD" \
--filter-expression "QD < 4.00" \
--filter-name "FS" \
--filter-expression "FS > 40.00" \
--filter-name "SOR" \
--filter-expression "SOR > 3.00" \
--output Variants.AllTillandsias.allChr.WGS.3bpindel.0.2missing.No_TEs_EDTA.MQ15_DP4.vcf.gz

echo "remove sites that had no genotype calls:"

conda activate gatk-vcf

gatk SelectVariants \
        -V  Variants.AllTillandsias.allChr.WGS.3bpindel.0.2missing.No_TEs_EDTA.MQ15_DP4.vcf.gz \
        --set-filtered-gt-to-nocall \
        --remove-unused-alternates \
        --output  Variants.AllTillandsias.allChr.WGS.3bpindel.0.2missing.No_TEs_EDTA.MQ15_DP4.NoCall.vcf.gz

conda activate my-env

### I want to filter for MAF, but this filter will remove all non-variant sites. In order to keep the nonvariant sites, and since the filter is only relevant for variant sites,
### in the next step I seperate the vcf into variant and non-variant. I then produce stats to see that the numbers make sense, and only then proceed filtering.

bcftools view -C0 Variants.AllTillandsias.allChr.WGS.3bpindel.0.2missing.No_TEs_EDTA.MQ15_DP4.NoCall.vcf.gz -Oz > Variants.AllTillandsias.allChr.WGS.3bpindel.0.2missing.No_TEs_EDTA.MQ15_DP4.NoCall.nonvariant.vcf.gz
tabix Variants.AllTillandsias.allChr.WGS.3bpindel.0.2missing.No_TEs_EDTA.MQ15_DP4.NoCall.nonvariant.vcf.gz

###filter MAF, this file will results in only variants

bcftools view -i 'MAF>0.045' Variants.AllTillandsias.allChr.WGS.3bpindel.0.2missing.No_TEs_EDTA.MQ15_DP4.NoCall.vcf.gz -Oz > Variants.AllTillandsias.allChr.WGS.3bpindel.0.2missing.No_TEs_EDTA.MQ15_DP4.MAF0.045.NoCall.onlyvariants.vcf.gz
tabix Variants.AllTillandsias.allChr.WGS.3bpindel.0.2missing.No_TEs_EDTA.MQ15_DP4.MAF0.045.NoCall.onlyvariants.vcf.gz

###merge the files

bcftools concat -a -D Variants.AllTillandsias.allChr.WGS.3bpindel.0.2missing.No_TEs_EDTA.MQ15_DP4.NoCall.nonvariant.vcf.gz Variants.AllTillandsias.allChr.WGS.3bpindel.0.2missing.No_TEs_EDTA.MQ15_DP4.MAF0.045.NoCall.onlyvariants.vcf.gz -Oz
 > Variants.AllTillandsias.allChr.WGS.3bpindel.0.2missing.No_TEs_EDTA.MQ15_DP4.MAF0.045.NoCall.vcf.gz

echo "last filter for missing and leave only sites that pass"

bcftools filter -i 'F_MISSING<0.2' Variants.AllTillandsias.allChr.WGS.3bpindel.0.2missing.No_TEs_EDTA.MQ15_DP4.MAF0.045.NoCall.vcf.gz | bcftools view -f .,PASS --trim-alt-alleles -Oz > Variants.AllTillandsias.allChr.WGS.3bpindel.0.2missin
g.No_TEs_EDTA.MQ15_DP4.MAF0.045.NoCall.PASS.vcf.gz

bcftools stats Variants.AllTillandsias.allChr.WGS.3bpindel.0.2missing.No_TEs_EDTA.MQ15_DP4.MAF0.045.NoCall.PASS.vcf.gz > Variants.AllTillandsias.allChr.WGS.3bpindel.0.2missing.No_TEs_EDTA.MQ15_DP4.MAF0.045.NoCall.PASS.vcf.gz_stats.txt
SnpSift tstv Variants.AllTillandsias.allChr.WGS.3bpindel.0.2missing.No_TEs_EDTA.MQ15_DP4.MAF0.045.NoCall.PASS.vcf.gz > Variants.AllTillandsias.allChr.WGS.3bpindel.0.2missing.No_TEs_EDTA.MQ15_DP4.MAF0.045.NoCall.PASS.vcf.gz_snpsift.txt

