for a in $(ls *txt)
do
PREFIX=$(basename $a .txt)
echo "starting run"
echo $PREFIX
/data/users/guzmans/soft/bayescenv-1.1/bin/linux64/bayescenv ../../2-OutlierDetection/Bayescan/Paxil_M095_noLD_Bayescan -env $a -od ./Results -o ${PREFIX}'_results' -nbp 40 -threads 32
echo "run"
echo $PREFIX
echo "DONE"
echo " "
python3 /data/users/guzmans/send_mail.py -s 'Bayescenv RUN Complete' -m $PREFIX 
done
