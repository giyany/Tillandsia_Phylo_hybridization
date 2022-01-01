#take vcf and bed files, create genomic windows in phylip, remove those with <n SNPs, estimate ML tree for each window

#!/bin/bash
#
#SBATCH -J makewins
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks-per-core=1
#SBATCH -p 
#SBATCH --qos

source /home/yardeni/.bashrc 
conda init bash

#directory
#dir=tmp4
#bed file used to define windows. I'm working with 10kb wins today
bed=10kb_wins_all.bed
#file I wanna slice
vcf=Variants.Only.twisst.run1.25scaffold.vcf.gz
#name the run because I like order
run="run4"
#min num of variable sites allowed. Any window with fewer SNPs will be discarded
minsnp=40

#conda activate my-env

#cd "$dir"

#bes is tab-delimited and I replaced the tabs so bcftools can read it and our filenames behave
sed 's/\t/:/' "$bed" | sed 's/\t/-/' | sort > winfile.sorted.txt

##create all windows
echo 'creating all windows'
while read region; do 
	echo "$region"; bcftools view -r "$region" ../Variants.Only.twisst.run1.25scaffold.vcf.gz | bgzip > "$region".window.vcf.gz; 
done < winfile.sorted.txt

#count snp nr.
echo 'counting snp nr'
while read region; do bcftools view "$region".window.vcf.gz| grep -v "#" | wc -l >> win_vcf_snpnr.txt; done < winfile.sorted.txt

#append into a very useful file

paste winfile.sorted.txt win_vcf_snpnr.txt > vcf_site_nr.txt

#remove those with <minsnp
echo 'removing wins with not enough snps'
awk -v min="$minsnp" '$2<min {print $1}' vcf_site_nr.txt > wins_to_remove.txt
awk -v min="$minsnp" '$2>=min {print $1}' vcf_site_nr.txt > wins_stay.txt

while read file; do rm "$file"*; done < wins_to_remove.txt

#convert to phylip and remove invariable sites created after IUPAC resolution, so we can use ascertainment bias correction for tree construction;
echo 'converting to phylip'
while read region; do 
	python /gpfs/data/fs71400/yardeni/Programs/vcf2phylip-master/vcf2phylip.py -i "$region".window.vcf.gz --resolve-IUPAC; 
	python /gpfs/data/fs71400/yardeni/Programs/raxml_ascbias-master/ascbias.py -p "$region".window.min4.phy -o "$region".window.ascbias.phy;
done < winfile.sorted.txt

rm *.felsenstein
rm *.stamatakis
#we're going to count the nr. of sites and again remove those with too few sites.
#grab the number of variable sites for each window from phylip head

for file in *.ascbias.phy; do
	awk 'NR==1 {print $2}' "$file" >> var_site_nr2.txt;
done

#append into a very useful file

paste wins_stay.txt var_site_nr2.txt > win_phylip_snp_nr.txt

##remove wins with <nr or snps

awk -v min="$minsnp" '$2<min' win_phylip_snp_nr.txt > wins_to_remove2.txt
while read file; do rm "$file"*; done < wins_to_remove2.txt
awk -v min="$minsnp" '$2>=min {print $1}' win_phylip_snp_nr.txt > wins_stay.txt

###construct tree and assess branch-support with ultra-fast bootstrap estimation. Don't estimate model, define it with -m. 
#My model is TVMe with 2 rates and ascertainment bias correction 

conda activate iqtree

for file in *ascbias.phy
do
iqtree -s "$file" -m TVMe+ASC+R2 -T 16 -B 1000
done

