###-----Add a prefix (folder name) to all files in the folder.
###-----Zhen Li

for i in Irma-*
do
    if [ -d "$i" ]
    then
          cd $i
          for a in *
          do
               name=$i"_"$a
               mv ./$a ./$name 
          done
          cd ../
    fi
done
