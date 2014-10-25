###Use this script to parse align_summary.txt file from tophat into a .txt file 
###containing informaiton on read depth.
###Zhen Li

import os

###Setup output file
cwd = os.path.basename(os.getcwd())
outfile = cwd + '_read_depth.txt' #Name of output file
outf = open(outfile, "w")
outf.write('Sample_Id \t Input_reads \t Mapped_reads \t Percentage_mapped \n')

###Iterate through subdirectories (1 level down) to find align_summary.txt file 
rootdir = os.getcwd()
dirs = os.walk(rootdir).next()[1]

for d in dirs:
    infile = rootdir + '/' + d + '/' + d + '_align_summary.txt'
    with open(infile) as inf:
        f = inf.read()
        line_words = f.splitlines()
        word_2 = line_words[1].split() 
        input_reads = word_2[2]         #Parse out input reads number
        word_3 = line_words[2].split()
        mapped_reads = word_3[2]       #Parse out output reads number
        percentage_mapped = float(mapped_reads) / float(input_reads)
        #Write to output file        
        outf.write('%s \t %s \t %s \t %s \n' %(d, 
                                               input_reads, 
                                               mapped_reads, 
                                               percentage_mapped))
