for i in ./Irma-*
do
    if [ -d "$i" ]
    then 
          name=${i%%\_*}
          name=${name#\./}
          mv ./$i ./$name 
    fi
done
