#!/bin/bash
#
#SBATCH -J vcf_filter
#SBATCH -N 2
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks-per-core=1

source .bashrc 
conda init bash
conda activate KING-env

#go to wd

cd /vcf;

plink2 --vcf Variants.AllTillandsias.allChr.WGS.raw.full.vcf --make-bed --allow-extra-chr --max-alleles 2 --set-all-var-ids @:# --out Variants.AllTillandsias.3bpindel.0.2miss.plink

mv Variants.AllTillandsias.3bpindel.0.2miss.plink.* ../KING_analysis/;
cd ../KING_analysis/
mv Variants.AllTillandsias.3bpindel.0.2miss.plink.bim Variants.AllTillandsias.3bpindel.0.2miss.plink.bim.old;

#Change chr names for king
sed -e 's/Scaffold\_//g' Variants.AllTillandsias.3bpindel.0.2miss.plink.bim.old > Variants.AllTillandsias.3bpindel.0.2miss.plink.bim;

king -b Variants.AllTillandsias.3bpindel.0.2miss.plink.bed --kinship
