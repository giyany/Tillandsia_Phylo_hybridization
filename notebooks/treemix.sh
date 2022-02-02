#ftools --gzvcf Treemix.64indv.nomissing.thinned.biallelic.vcf.gz --plink --mac 2 --max-alleles 2 --chrom-map Treemix.64indv.nomissing.thinned.biallelic.chrom-map.txt --out Treemix.64indv.nomissing.thinned.biallelic
his is not edited yet
#requires clustl file: tab-delimited, 1st and 2nd column are the indv name, 3rd column is the species assignment

#file produced: I removed indvs with missing data >10%, leaving me with 64 indvs. only biallelic SNPs, no missing data. Didn't do LD pruning with plink, but thinned the data every 1000snps with vcftools. This leaft 9.8k sites
#following https://github.com/speciationgenomics/scripts/blob/master/vcf2treemix.sh
#create chromosome map for the next stage (vcftools and other apps don't like our chromosome names)
bcftools view -H Treemix.64indv.nomissing.thinned.biallelic.vcf.gz | cut -f 1 | uniq | awk '{print $0"\t"$0}' > Treemix.64indv.nomissing.thinned.biallelic.chrom-map.txt

#get map and ped files for plink using the chromosome map we just created
vcftools --gzvcf Treemix.64indv.nomissing.thinned.biallelic.vcf.gz --plink --mac 2 --max-alleles 2 --chrom-map Treemix.64indv.nomissing.thinned.biallelic.chrom-map.txt --out Treemix.64indv.nomissing.thinned.biallelic

#adjust the map to allow for non-human chromosome names

awk -F"\t" '{
        split($2,chr,":")
	$1="1"
	$2="1:"chr[2]
	        print $0
	}' Treemix.64indv.nomissing.thinned.biallelic.map > better.map

mv better.map Treemix.64indv.nomissing.thinned.biallelic.map

#convert files to stratified frq files
plink --file Treemix.64indv.nomissing.thinned.biallelic --make-bed --out Treemix.64indv.nomissing.thinned.biallelic --allow-no-sex --allow-extra-chr 0

plink -bfile Treemix.64indv.nomissing.thinned.biallelic --freq --missing --within clustfile.nomissing.txt --out Treemix.64indv.nomissing.thinned.biallelic --allow-no-sex --allow-extra-chr 0

gzip Treemix.64indv.nomissing.thinned.biallelic.frq.strat

./plink2treemix.py Treemix.64indv.nomissing.thinned.biallelic.frq.strat.gz Treemix.64indv.nomissing.thinned.biallelic.frq

gunzip Treemix.64indv.nomissing.thinned.biallelic.frq.gz
gunzip Treemix.64indv.nomissing.thinned.biallelic.frq.strat.gz

awk 'BEGIN{print "scaffold_pos\tscaffold\tpos"}{split($2,pos,":");print $2"\t"pos[1]"\t"pos[2]}' Treemix.64indv.nomissing.thinned.biallelic.map > Treemix.64indv.nomissing.thinned.biallelic.position

paste Treemix.64indv.nomissing.thinned.biallelic.position Treemix.64indv.nomissing.thinned.biallelic.frq > Treemix.64indv.nomissing.thinned.biallelic.frequencies

awk '{printf $0
for(i = 4; i <= NF; i++){
	split($i,values,",")
	if((values[1]+values[2])>0) freq=values[1]/(values[1]+values[2])
	else freq=0
	printf freq"\t" 
} printf "\n"}'  Treemix.64indv.nomissing.thinned.biallelic.frequencies >  Treemix.64indv.nomissing.thinned.biallelic.frequencies2

mv Treemix.64indv.nomissing.thinned.biallelic.frequencies2 Treemix.64indv.nomissing.thinned.biallelic.frequencies

awk 'BEGIN{scaffold="";pos=0;newpos=0}
{if($2==scaffold){newpos=pos+$3}else{scaffold=$2;pos=newpos};chpos=pos+$3;print $0,chpos}' \
	Treemix.64indv.nomissing.thinned.biallelic.frequencies > Treemix.64indv.nomissing.thinned.biallelic.frequencies.newpos

gzip Treemix.64indv.nomissing.thinned.biallelic.frq

