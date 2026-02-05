for a in $(ls *txt)
do
PREFIX=$(basename $a .txt)
echo "starting run"
echo $PREFIX
bayescenv ../../2-OutlierDetection/Bayescan/Paxil_M095_noLD_Bayescan -env $a -od ./Results -o ${PREFIX}'_results' -nbp 40 -threads 64
echo "run"
echo $PREFIX
echo "DONE"
echo " "
done
