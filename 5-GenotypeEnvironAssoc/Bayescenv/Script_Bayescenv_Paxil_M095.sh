for a in $(ls *txt)
do
echo "starting run"
echo $a
/data/users/guzmans/bayescenv-1.1/bin/linux64/bayescenv /data/users/guzmans/Bayescenv_paxil_M095_010624/Paxil_M095_noMono_noLD01_BAYESCAN_input -env $a -od /data/users/guzmans/Bayescenv_paxil_M095_010624/Results -o $a'_results' -threads 12
echo "run"
echo $a
echo "finished"
echo " "
done
