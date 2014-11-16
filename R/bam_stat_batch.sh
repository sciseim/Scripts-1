for i in $HOME/Data/Ouyang/PEC_CT_CGATGT_L007_R1_all $HOME/Data/Ouyang/PEC_CoCl2_GCCAAT_L007_R1_all
do
    cd $i
    bam_stat.py -i $i/accepted_hits.bam > bam_stat.txt
    cd ../
    read -p "Press any key to continue..." -n 1 -s
done
