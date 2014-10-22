for i in Irma-*
do
    if [ -d "$i" ]
    then 
          echo "./"$i"/"$i"_transcripts.gtf" >> assembly_list.txt
    fi
done
