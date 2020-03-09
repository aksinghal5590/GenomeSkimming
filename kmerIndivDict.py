# input is a fastq file
from Bio import SeqIO
import os, glob

## read kmerfile into a list
#kmerfile = ["GCG", "AAA"]
kmerfile = []
for record in SeqIO.parse("kmers.fa", "fasta"):
    kmerfile.append(str(record.seq))
kmerfile.sort()

## read readfile into reads list
os.chdir("reads_fa")

for filename in glob.glob("indiv_*.fa"):
    print("opening" + filename)

    reads = [] # fastq file
    #filename = "indiv_0.fa"
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
    
    # create the kmer --> list of (reads) dictionary
    for read in range(0, len(reads)):
        if kmerlen > 0:
            for start in range(0, len(reads[read])-kmerlen+1):
                readWindow = reads[read][start:start+kmerlen]
                if kmerDict.get(readWindow) is not None:
                    val = kmerDict.get(readWindow)
                    val.append((read, start))
                    kmerDict[readWindow] = val 
    
    print("removing kmer collisions with the exact same set of reads...")
    # create the list(reads-->kmer) dictionary,
    # remove any kmers with collisions from the kmerDict dictionary
    # remove kmers with more than 8x coverage
    kmersToRemove = []
    MAXCOVERAGE = 8
    
    readsListDict = {}
    for kmer in kmerDict.keys():
        if len(kmerDict[kmer]) > MAXCOVERAGE:
            kmersToRemove.append(kmer)
        else:
            readList = tuple(set([read for (read, start) in kmerDict[kmer]]))
            if readsListDict.get(readList) is not None:
                kmersToRemove.append(kmer)
            else:
                readsListDict[readList] = kmer
    
    print("kmersToRemove: " + str(len(kmersToRemove)))
            
    for kmer in kmersToRemove:
        del kmerDict[kmer] # this will throw an error if we try to remove a key that is not in the dictionary
                        
    print("NEW kmerDictLen: " + str(len(kmerDict.keys())))
    
    # smaller code for testing purposes only
    #kmerDict = {"AAA":[(0, 6), (1, 2), (2, 3)]}
    #reads = ["TTTTTCAAA", "TCAAA-CGT", "AGCAAAAGGT"]
    #print(len(kmerDict.keys()))
    
    emptyCtr = 0
    for kmer in kmerDict.keys():
        if len(kmerDict[kmer]) > 0:
            max_bases_before = -1
            max_bases_after = -1
            max_bases_before_READ_ind = -1
            max_bases_after_READ_ind = -1
            readsUnalignedList = [reads[r] for (r, s) in kmerDict[kmer]]
            startPosList = [s for (r, s) in kmerDict[kmer]]
            max_bases_before = max(startPosList)
            for s in range (0, len(startPosList)):
                if startPosList[s] == max_bases_before:
                    max_bases_before_READ_ind = s
                if len(readsUnalignedList[s]) - startPosList[s] > max_bases_after:
                    max_bases_after = len(readsUnalignedList[s]) - startPosList[s]
                    max_bases_after_READ_ind = s
        
            # now combine the two reads
            contig = ""
            #                                                  read[:]
            contig += readsUnalignedList[max_bases_before_READ_ind][0:startPosList[max_bases_before_READ_ind]]
            #                                                 read[:]
            contig += readsUnalignedList[max_bases_after_READ_ind][startPosList[max_bases_after_READ_ind]:]
        
            kmerDict[kmer] = contig
        else:
            emptyCtr += 1
            
    print("emptyCtr: " + str(emptyCtr))
    print("kmerDict.keys() len: "+ str(len(kmerDict.keys())))
    
    # write sorted kmers to file
    f = open(filename+"_kmerDict.txt", "w")
    for kmer in kmerfile:
        toWrite = kmer + "\t"
        if kmerDict.get(kmer) is not None and kmerDict.get(kmer) != []:
            toWrite += kmerDict[kmer] + "\n"
        else:
            toWrite += "\n"
        f.write(toWrite)
    f.close()



#f = open(filename+"_kmerDict.txt", "w")
#f.write("\nReads List\n")
#for read in reads:
#    f.write(read+"\n")
#    
#f.write("\nKmerDictionary (key, value)\n")
#for key in kmerDict.keys():
#    finalString = ""
#    for val in kmerDict[key]:
#        finalString += str(val) + "\t"
#    f.write(key+"\t"+finalString+"\n")
#                
#f.close()

os.chdir("..")
