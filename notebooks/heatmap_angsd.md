# Tillandsia whole genome data: heatmap for major clades
Updated: 1/6/2021

Required: [ANGSD](https://github.com/ANGSD/angsd), [pcangsd](http://www.popgen.dk/software/index.php/PCAngsd)

Create a heatmap of relatedness with angsd. We did this before merging libraries of the same accession which were prepared separately, to ensure they're the same samples (lab switches happen).

## List of metadata files
- bam_list.txt - list of bam files. Create with:
```bash
> ls ./*rg.bam | sed -e 's/\.//' -e 's/\///g' > bam_list.txt
```
- ind_list.txt - list of your samples, same order as in bam_list.

## run ANGSD

infer gentype likelihood. using GATK model (-GL), maximum likelihood approach to choose the major and minor alleles (doMajorMinor), assume known major and minor alleles (-doMaf), minInd adjusted to number of samples.

```bash
> angsd -bam bam_list.txt -GL 2 -doMajorMinor 1 -P 16 -doMaf 1 -SNP_pval 1e-6 -minMapQ 20 -minQ 20 -minInd 10 -minMaf 0.027 -doGlf 2 -out gatk_all_samples  
```

calculate covariate matrix 
```python
> python pcangsd.py -beagle gatk_all_samples.beagle.gz -o gatk_all_samples_pcangsd_out -minMaf 0.027
```
plot heatmap. The clade I work on is clade K, hence the name

```R
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(hrbrthemes)
library(reshape2)

gatk_cladeK<-as.matrix(read.table("gatk_all_samples.cov"))
inds<-read.table("ind_list.txt")
inds<-paste0(inds[,1])

#replace diagonal with NA
cladeK_covar<- gatk_cladeK %>% 
  replace(., col(.) == row(.), NA)
  
#assign names
rownames(cladeK_covar)<-inds
colnames(cladeK_covar)<-inds

#melt for ggplot2's needs. This changes the format of your dataset
melted_cladeK<-melt(cladeK_covar)


ggplot(melted_cladeK, aes(Var1, Var2, fill= value)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu") +
  theme_minimal()
```
