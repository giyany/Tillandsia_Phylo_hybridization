#!/bin/bash

# create genomic windows from VCF file, construct maximum-likelihood tree for each window. Concatenate the output gene trees and use ASTRAL to infer a species tree
#define
dir=./
bedfile=RefBeet-3.0.bed
vcf=/project/cultibv/yardeni/EasternBeets/call_RefBeet-3.0/EasternBeets.filter.harder.full.vcf.gz
run="wins.10k"
size=10000
minsnp=40
vcf2phylip_path=/home/yardeni/Programs/vcf2phylip
ascbias_path=/home/yardeni/Programs/raxml_ascbias
fastagap=/home/yardeni/Programs/fastagap
phylip2fasta=/home/yardeni/Programs/Sequence-manipulation


#create bed with windows
echo 'generating bed file'
 
#bedtools makewindows -b $bedfile -w $size > "$bedfile"."$run".bed

echo 'creating all windows'

#awk '{print $1":"$2"-"$3}' "$bedfile"."$run".bed | parallel --jobs 5 "bcftools view -r {} $vcf -Oz > {}.window.vcf.gz"

echo 'calculating SNP content... (this can take a while)'

#for file in *.window.vcf.gz; do printf "$file\t" >> all_vcf_win_site_count.txt; zcat "$file" | grep -v '#' | wc -l >> all_vcf_win_site_count.txt; done

echo 'remove windows with less than '$minsnp' sites'

#awk -v min="$minsnp" '$2<min {print $1}' all_vcf_win_site_count.txt > wins_to_remove.txt
#awk -v min="$minsnp" '$2>=min {print $1}' all_vcf_win_site_count.txt > wins_stay.txt

#while read file; do rm $file; done < wins_to_remove.txt;

echo 'convert to phylip'

#awk '{print $1":"$2"-"$3}' "$bedfile"."$run".bed | parallel --jobs 5 "$vcf2phylip_path/vcf2phylip.py -i {}.window.vcf.gz --resolve-IUPAC; python $ascbias_path/ascbias.py -p {}.window.min4.phy -o {}.window.ascbias.phy"

awk '{print $1":"$2"-"$3}' "$bedfile"."$run".bed | parallel --jobs 5 "python $ascbias_path/ascbias.py -p {}.window.min4.phy -o {}.window.ascbias.phy"

for file in *.felsenstein; do rm $file; done
for file in *.stamatakis;do rm $file; done

echo 'find and remove phylip windows with too few variant sites'

for file in *.ascbias.phy; do printf "$file\t" >> all_phy_win_site_count.txt; awk 'NR==1 {print $2}' "$file" >> all_phy_win_site_count.txt; done

awk -v min="$minsnp" '$2<min {print $1}' all_phy_win_site_count.txt > wins_to_remove2.txt;
while read file; do rm $file; done < wins_to_remove2.txt

echo 'find and remove taxa with fragmantary data (gaps) above 20%, 35%, 50%, 65%'

#convert phylip to fasta

ls -l | grep ascbias.phy | awk '{print $9}' | sed 's/\.phy//g' >> phy_wins_stay.txt;
while read file; do perl $phylip2fasta/Phylip2Fasta.pl $file.phy $file.fasta; done < phy_wins_stay.txt

while read file; do perl $fastagap/fastagap.pl -PA=20 --missing="N" $file.fasta > $file.gap20.fasta; done < phy_wins_stay.txt
mkdir gap20/
mv *gap20.fasta gap20/
while read file; do perl $fastagap/fastagap.pl -PA=35 --missing="N" $file.fasta > $file.gap35.fasta; done < phy_wins_stay.txt
mkdir gap35/
mv *gap35.fasta gap35/
while read file; do perl $fastagap/fastagap.pl -PA=50 --missing="N" $file.fasta > $file.gap50.fasta; done < phy_wins_stay.txt
mkdir gap50/
mv *gap50.fasta gap50/
while read file; do perl $fastagap/fastagap.pl -PA=65 --missing="N" $file.fasta > $file.gap65.fasta; done < phy_wins_stay.tx
mkdir gap65/
mv *gap65.fasta gap65/

echo 'inferring trees on all windows (this can take a while)'
echo 'note that logs will be written into *.log file, not stdout'

cd gap20/
ls *.fasta | parallel --jobs 5 "iqtree2 -s {} -B 1000 --quiet"
cd ..
cd gap35/
ls *.fasta | parallel --jobs 5 "iqtree2 -s {} -B 1000 --quiet"
cd ..
cd gap50/
ls *.fasta | parallel --jobs 5 "iqtree2 -s {} -B 1000 --quiet"
cd ..
cd gap65/
ls *.fasta | parallel --jobs 5 "iqtree2 -s {} -B 1000 --quiet"
cd ..

