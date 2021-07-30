#!/bin/bash
#
#SBATCH -J vcf_filter
#SBATCH -N 2
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks-per-core=1
#SBATCH -p mem_0096
#SBATCH --qos p71400_0096

source /home/fs71400/yardeni/.bashrc 
conda init bash
conda activate KING-env

#go to wd

cd /gpfs/data/fs71400/yardeni/WGS/vcf;

plink2 --vcf Variants.AllTillandsias.allChr.WGS.raw.full.vcf --make-bed --allow-extra-chr --max-alleles 2 --out Variants.AllTilandsias.raw.full
