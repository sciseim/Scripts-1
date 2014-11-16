##----running cufflinks

for i in $HOME/Data/Ouyang/PEC*/accepted_hits.bam
do
    dir=${i%%.fastq}
    cufflinks -o $dir --library-type fr-firststrand -G ~/genome_and_annotation/gencode.vM3.annotation.gtf -p 30 $i 
done