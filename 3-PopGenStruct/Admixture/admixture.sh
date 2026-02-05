for k in $(seq 11)
do
admixture Paxil_M095_PutatNeutral.PLINK.bed $k --cv -j16 -B1000 | tee Paxil_M095_PutatNeutral__Admixture_log.${k}.out
done
