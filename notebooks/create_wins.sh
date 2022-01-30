#!/bin/bash
#
#SBATCH -J makewins
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks-per-core=1
#SBATCH -p mem_0384
#SBATCH --qos p71400_0384

source /home/fs71400/yardeni/.bashrc 
conda init bash

#go to wd

#directory
dir=/gpfs/data/fs71400/yardeni/WGS/astral/gen_wins_noSA/
#file used to define windows. I'm working with 100kb wins today
winfile=100kb_wins_all.sorted.txt
#file I wanna slice
vcf=../Only.Variants.onlyMexican.allChr.WGS.3bpindel.0.2missing.No_TEs_EDTA.MQ15_DP4.MAF0.045.NoCall.PASS.vcf.gz
#name the run because I like order
run="wins.100kb.no.SA"
#min num of variable sites allowed. Any window with fewer SNPs will be discarded
minsnp=40

conda activate my-env

cd "$dir"

rm *.vcf.gz
rm *.txt
rm *.phy
### my file was tab-delimited and I replaced the tabs so bcftools can read it
### sed 's/\t/:/' ../tillandsia_fasciculata_25_scaffold.windows.100k.bed | sed 's/\t/-/' | sort > 100kb_wins_all.sorted.txt

## create all windows
echo 'creating all windows'

while read region; do 
	echo "$region"; bcftools view -r "$region" "$vcf" | bgzip > "$region".window.vcf.gz; 
done < ../"$winfile"

mkdir vcf_backup
cp *.vcf.gz vcf_backup;

#### count snp nr.
echo 'counting snp nr'
while read region; do bcftools view "$region".window.vcf.gz| grep -v "#" | wc -l >> all_win_vcf_snpnr.txt; done < ../"$winfile"

#### append into a very useful file

paste ../"$winfile" all_win_vcf_snpnr.txt > all_vcf_files_with_snpnr.txt

####remove those with <minsnp
echo 'removing wins with not enough snps'
awk -v min="$minsnp" '$2<min {print $1}' all_vcf_files_with_snpnr.txt > wins_to_remove.txt
awk -v min="$minsnp" '$2>=min {print $1}' all_vcf_files_with_snpnr.txt > wins_stay.txt

while read file; do rm "$file"*; done < wins_to_remove.txt

### convert to phylip and remove invariable sites created after IUPAC resolution, so we can use ascertainment bias correction for tree construction;
echo 'converting to phylip'
while read region; do 
	python /gpfs/data/fs71400/yardeni/Programs/vcf2phylip-master/vcf2phylip.py -i "$region".window.vcf.gz --resolve-IUPAC; 
	python /gpfs/data/fs71400/yardeni/Programs/raxml_ascbias-master/ascbias.py -p "$region".window.min4.phy -o "$region".window.ascbias.phy;
done < ../"$winfile"

rm *.felsenstein
rm *.stamatakis

### we're going to count the nr. of sites and again remove those with too few sites.
### grab the number of variable sites for each window from phylip head

for file in *.ascbias.phy; do
	awk 'NR==1 {print $2}' "$file" >> all_phy_snp_nr.txt;
done

####append into a very useful file

paste wins_stay.txt all_phy_snp_nr.txt > win_phylip_snp_nr.txt

#####remove wins with <nr or snps

awk -v min="$minsnp" '$2<min {print $1}' win_phylip_snp_nr.txt > wins_to_remove2.txt
while read file; do rm "$file"*; done < wins_to_remove2.txt
awk -v min="$minsnp" '$2>=min {print $1}' win_phylip_snp_nr.txt > wins_stay.txt

conda activate iqtree

for file in *ascbias.phy
do
iqtree -s "$file" -m TVMe+ASC+R2 -T 16 -B 1000
done

#not all the files may be used to create windows, because sometimes there is literally only missing data for one of the samples and iqtree will not estimate a tree then.
#I'll collect the final list of trees:

printf "scaffold\tstart\tend\n" > output."$run".data.tsv
ls *.treefile | cut -d"." -f 1 | sed -e 's/:/\t/' -e 's/-/\t/' >> output."$run".data.tsv

