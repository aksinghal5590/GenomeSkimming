# input is a fastq file
from Bio import SeqIO
import os, glob

## read kmerfile into a list
#kmerfile = ["GCG", "AAA"]
kmerfile = []
for record in SeqIO.parse("kmers.fa", "fasta"):
    kmerfile.append(str(record.seq))


## read readfile into reads list
os.chdir("reads_fa")
for filename in glob.glob("indiv_*.fa"):
    print("opening" + filename)

    reads = [] # fastq file

    for record in SeqIO.parse(filename, "fasta"):
        reads.append(str(record.seq))
    
    print("number of reads: " + str(len(reads)))
    #print(reads[0])
    
    #reads = ["AAAGCGGTT", "GCG", "ACTGCG", "GGTAAA"]
    kmerDict = {} # (key, value) --> (kmer, list of indices of reads from reads)
    kmerlen = len(kmerfile[0])
    
    # read kmer file and fill dictionary
    for kmer in kmerfile:
        kmerDict[kmer] = []
        
    print("kmerDictLen: " + str(len(kmerDict.keys())))
    
    for read in range(0, len(reads)):
        if kmerlen > 0:
            for start in range(0, len(reads[read])-kmerlen+1):
                readWindow = reads[read][start:start+kmerlen]
                if kmerDict.get(readWindow) is not None:
                    val = kmerDict.get(readWindow)
                    val.append((read, start))
                    kmerDict[readWindow] = val     
    
                
    f = open(filename+"_kmerDict.txt", "w")
    f.write("\nReads List\n")
    for read in reads:
        f.write(read+"\n")
        
    f.write("\nKmerDictionary (key, value)\n")
    for key in kmerDict.keys():
        finalString = ""
        for val in kmerDict[key]:
            finalString += str(val) + "\t"
        f.write(key+"\t"+finalString+"\n")
                    
    f.close()

os.chdir("..")
