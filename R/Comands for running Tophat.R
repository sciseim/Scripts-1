##-----------------Commands for running Tophat
##-----------------121913

##-----------------combine all file into one
cat *.fastq.gz >> AF3.fastq.gz

##-----------------running tophat
tophat -g 1 -G ~/gencode.vM2.annotation.gtf -p 10 ~/bin/bowtie2-2.1.0/indexes/mm10 AF3.fastq.gz