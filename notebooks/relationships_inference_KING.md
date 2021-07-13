# Tillandsia whole genome data: using KING for relationship inference
Updated: 12/07/2021

Required: [plink2](https://www.cog-genomics.org/plink/2.0/)

Looking at relatedness before merging files of samples that we expect are the same accession. Let's estimate kinship coefficients.

## List of metadata files
- Variants.AllTillandsias.allChr.WGS.raw.full.vcf - variant calls, NOT PRUNED. I left indels.

#produce PLINK binary format: binary genotype file, a family file, and a map file.


ink2 --vcf Variants.AllTillandsias.allChr.WGS.raw.full.vcf --make-bed --allow-extra-chr --max-alleles 2 --out Variants.AllTilandsias.raw.full 



plink2 --vcf Variants.AllTillandsias.allChr.WGS.raw.full.vcf --make-bed --allow-extra-chr --max-alleles 2 --out Variants.AllTilandsias.raw.full 




