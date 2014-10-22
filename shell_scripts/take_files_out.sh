for i in Irma-*
do
    if [ -d "$i" ]
    then
          cd $i
          mv ./$i/* ./
          cd ../
    fi
done
