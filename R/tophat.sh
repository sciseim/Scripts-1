
for i in *.fastq
do
    dir=${i%%.fastq}
	mkdir $dir
    tophat -o $dir -g 1 --library-type fr-firststrand -G ~/genome_and_annotation/gencode.vM3.annotation.gtf -p 30 ~/bin/bowtie2-2.1.0/indexes/mm10 $i 
done
