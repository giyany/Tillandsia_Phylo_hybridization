#!/bin/bash
#
#SBATCH -J iqtree
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -p mem_xxx
#SBATCH --qos p71400_xxx
#SBATCH --array=1174,1383,15,2041,2073,2096,2199,2225,2287,2312,2313,2314,2315,2316,2317,2318,2319,2320,2321,326,346,53,658,819,845

#module load ### 
source /home/.bashrc

i=$SLURM_ARRAY_TASK_ID
conda init bash
conda activate my-env

## move to working directory ##

cd /gpfs/data/fs71400/yardeni/WGS/vcf/per_chr

#convert vcf to phylip file

#python ../../../Programs/vcf2phylip-master/vcf2phylip.py -i Variants.Only.AllTillandsias.Scaffold_"$i".WGS.3bpindel.0.2missing.No_TEs_EDTA.MQ15_DP4.MAF0.045.NoCall.PASS.vcf.gz -f;

#python ../../../Programs/raxml_ascbias-master/ascbias.py -p Variants.Only.AllTillandsias.Scaffold_"$i".WGS.3bpindel.0.2missing.No_TEs_EDTA.MQ15_DP4.MAF0.045.NoCall.PASS.min4.phy -o Variants.Only.AllTillandsias.Scaffold_"$i".WGS.3bpindel.0.2missing.No_TEs_EDTA.MQ15_DP4.MAF0.045.NoCall.PASS.ascbias.min4.phy;

#mv Variants.Only.AllTillandsias.Scaffold_"$i".WGS.3bpindel.0.2missing.No_TEs_EDTA.MQ15_DP4.MAF0.045.NoCall.PASS.ascbias.min4.phy ../../iqtree;

cd ../../iqtree;

conda activate iqtree

#estimate best-fit model, construct tree and assess branch-support with ultra-fast bootstrap estimation. Using a file with all sites.
iqtree -s Variants.Only.AllTillandsias.Scaffold_"$i".WGS.3bpindel.0.2missing.No_TEs_EDTA.MQ15_DP4.MAF0.045.NoCall.PASS.ascbias.min4.phy -m TVMe+ASC+R2 -T AUTO -B 1000 --redo

#construct ML tree and determine best-fit model
#iqtree -s Variants.Only.Tillandsia.65indv.1outgrp.WGS.3bpindel.0.2missing.No_TEs_EDTA.MQ15_DP4.MAF0.045.NoCall.PASS.min4.ascbias.modelfind.phy -T 32 -m MFP
