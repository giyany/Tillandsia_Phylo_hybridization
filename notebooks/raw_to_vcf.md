# Tillandsia whole genome data: from sequencer to VCF
Updated: 9/3/2021

Required: [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), [deML](https://github.com/grenaud/deML)

Description and goals

## List of metadata files
- P5P7_for_demulti - tab delimited file with accession and indices used. See description in the [deML format examples](https://github.com/grenaud/deML)

## demultiplex

VBCF gave us 
```bash
> deML -i P5P7_for_demulti --bamtags BC,QT,B2,Q2 -o Lib3-demulti_deML.bam -s demult_stats.txt -e demult_unassigned.txt Lib3_raw_full.bam
```
ddd
