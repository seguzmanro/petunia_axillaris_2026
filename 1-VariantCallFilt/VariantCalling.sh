export REF_GENOME="/data/users/guzmans/petunia_axillaris_2026/Peaxi162_genome.fa"
PREFIX_SP="Paxil"
PREFIX_MISS="M095"

mkdir -p 01_GBSX 02_fastq-mcf 03_fastq-stats 04_bwa 05_samtools 06_bamaddrg 07_freebayes

cd 01_GBSX

# Demultiplexing
java -jar GBSX_v1.2.jar --Demultiplexer -f1 ../Plate_5_R1.fastq.gz -i ../barcodes.txt -gzip true

cd ../02_fastq-mcf

# Use fastq-mcf to clean raw reads and not allow per-base quality < 30 nor reads shorter than 50
ls ../01_GBSX/ | grep fastq | sed 's/.fastq.gz//' | awk 'BEGIN{ print "#!/bin/bash\n\n"}{print "fastq-mcf -q 30 -l 50 -o "$1"_Q30L50.fq ../../IlluminaAdapters_V2.fasta ../01_GBSX/"$1".fastq.gz\n"}' > Run_ProcFiles.sh

bash Run_ProcFiles.sh

cd ../03_fastq-stats

# Get basic stats for processed reads from each sample
ls ../02_fastq-mcf | grep _Q30L50.fq | sed 's/_Q30L50.fq//' | awk 'BEGIN{ print "#!/bin/bash\n\n"}{print "fastq-stats ../02_fastq-mcf/"$1"_Q30L50.fq > "$1"_Q30L50.stats.txt"}' > Run_GetQ30L50stats.sh

bash Run_GetQ30L50stats.sh

grep reads *stats.txt | sed 's/.R1_Q30L50.stats.txt:reads//' > reads_count.txt
grep "total bases" *stats.txt | sed 's/.R1_Q30L50.stats.txt:total bases//' > totalbases_count.txt
join reads_count.txt totalbases_count.txt | sed -r 's/ /\t/g' > Q30L50reads_tbases.txt

cd ../04_bwa

bwa index $REF_GENOME

# Map the reads from each sample to the reference genome
ls ../02_fastq-mcf/ | grep _Q30L50.fq | sed 's/_Q30L50.fq//' | awk 'BEGIN{ print "#!/bin/bash\n\n"}{print "bwa mem -t 64 ", ENVIRON["REF_GENOME"], " ../02_fastq-mcf/"$1"_Q30L50.fq | samtools view -F 4 -Sb -o "$1".bam -"}' > Run_Mapping.sh

bash Run_Mapping.sh

cd ../05_samtools

# Sort individual bam files
ls ../04_bwa | grep .bam | sed 's/.bam//' | awk 'BEGIN{ print "#!/bin/bash\n\n"}{print "samtools sort -@ 64 -o "$1".bam ../04_bwa/"$1".bam"}' > Run_SortBam.sh

bash Run_SortBam.sh

cd ../06_bamaddrg

# Use bamaddrg to merge all bam files (I installed bamaddrg through conda and called it here using its full path)
ls ../05_samtools | grep bam | sed -r 's/.bam//' | awk 'BEGIN{ print "#!/bin/bash\n\n"; a="/data/users/guzmans/soft/miniconda3/envs/bamaddrg/bin/bamaddrg"}{a=a" -b ../05_samtools/"$1".bam -s "$1;}END{ print a " > merge_paxil.bam"}' > Run_bamaddrg.sh

bash Run_bamaddrg.sh

cd ../07_freebayes

samtools index ../06_bamaddrg/merge_paxil.bam

# Call Variants using freebayes
freebayes -b ../06_bamaddrg/merge_paxil.bam -f $REF_GENOME -v ${PREFIX_SP}Variants.vcf --min-mapping-quality 20 --min-base-quality 30 -C 10

# Remove variants other than SNPs
awk '{ if ($8 !~ /mnp/) print $0}' ${PREFIX_SP}Variants.vcf | awk '{ if ($8 !~ /complex/) print $0}' | awk '{ if ($8 !~ /ins/) print $0}' | awk '{ if ($8 !~ /del/) print $0 }' | awk '{ if ($8 !~ /snp,/) print $0}' > ${PREFIX_SP}OnlySNPs.vcf
echo $(grep -v '#' ${PREFIX_SP}OnlySNPs.vcf | wc -l)

vcf-stats ${PREFIX_SP}OnlySNPs.vcf -p stats/${PREFIX_SP}OnlySNPs

# Make sure we allow for 5% missing data at most by loci
vcftools --vcf ${PREFIX_SP}OnlySNPs.vcf --recode --recode-INFO-all --max-missing 0.95 --maf 0.01 --out ${PREFIX_SP}_${PREFIX_MISS}_noMono
echo $(grep -v '#' ${PREFIX_SP}_${PREFIX_MISS}_noMono.recode.vcf | wc -l)

vcf-stats ${PREFIX_SP}_${PREFIX_MISS}_noMono.recode.vcf -p stats/${PREFIX_SP}_${PREFIX_MISS}_noMono

# Filter bad individuals by the number total number of SNPs available:
awk '{ if ($1 < 10000) print $2 }' stats/${PREFIX_SP}_${PREFIX_MISS}_noMono.counts | sed 's/samples\///' > bad_indv
vcftools --vcf ${PREFIX_SP}_${PREFIX_MISS}_noMono.recode.vcf --remove bad_indv --recode --recode-INFO-all --out ${PREFIX_SP}_${PREFIX_MISS}_indvFilt

# Compute R^2 for all sites in 50kb windows
vcftools --vcf ${PREFIX_SP}_${PREFIX_MISS}_indvFilt.recode.vcf --ld-window-bp 50000 --geno-r2 --max-alleles 2 --min-alleles 2 --out ${PREFIX_SP}_${PREFIX_MISS}_indvFilt_LD

# Remove sites with LD R^2 > 0.1:
awk 'NR==1 || $NF > 0.1' ${PREFIX_SP}_${PREFIX_MISS}_indvFilt_LD.geno.ld | cut -f1,2 > sites_ld_higher_01
vcftools --vcf ${PREFIX_SP}_${PREFIX_MISS}_indvFilt.recode.vcf --exclude-positions sites_ld_higher_01 --recode --recode-INFO-all --out ${PREFIX_SP}_${PREFIX_MISS}_noLD

# rename sample names in all vcf files to make sure they coincide exactly with the map files throughout analyses:
for vcf in $(ls *vcf); do sed -i 's/.R1//g' $vcf ; done