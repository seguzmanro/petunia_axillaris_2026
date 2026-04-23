mkdir -p Paxil_M095_PutatNeutral
populations -V ../../2-OutlierDetection/Paxil_M095_PutatNeutral.recode.vcf.gz -M ../../Paxil_PopGroupMap.STACKS2.txt -O Paxil_M095_PutatNeutral --threads 16

mkdir -p Paxil_M095_noLD
populations -V ../../1-VariantCallFilt/07_freebayes/Paxil_M095_noLD.recode.vcf.gz -M ../../Paxil_PopGroupMap.STACKS2.txt -O Paxil_M095_noLD --threads 16