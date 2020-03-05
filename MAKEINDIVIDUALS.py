
#==============================================================================
# read in the reference genome as a single string
#==============================================================================
import os, glob
import Bio
from Bio import SeqIO

refgenome = ""

filepath = os.getcwd() + "/chr1.fa"
for record in SeqIO.parse(filepath, "fasta"):
    refgenome += record.seq
    
print(type(refgenome))

#==============================================================================
# only need this once: write reference genome as a single text string to a file
#==============================================================================
#f = open("ref_genome_c_elegans.txt", "w+")
#f.write(str(refgenome))
#f.close()

#==============================================================================
# read in the vcf file and create a dictionary of format:
# (key, value) -> (position, [G, C, C, G...] list of length # individuals, 40 in this case)
#==============================================================================

# helper function for when we have two rows in the vcf file at the same position
# this means we have multiple alternative alleles
# so we need to combine the two
def combineRows(ref, oldValues, newValues):
    new = oldValues
    for i in range(0, len(oldValues)):
        if oldValues[i] != newValues[i] and newValues[i] != ref:
            new[i] = newValues[i]
    return new      
#print(combineRows("A", ["A", "A", "G"], ["A", "C", "A"]))

vcfDict = {}
f = open(os.getcwd() + "/I.vcf") # "/smaller.txt") 
for line in f:
    parsed = line.strip().replace('|', '\t').split('\t')
    if len(parsed) == 29+20 and parsed[0] != "#CHROM":
        key = parsed[1]
        ref = parsed[3]
        alt = parsed[4]
        value = []
        for i in range(9,len(parsed)):
            if parsed[i] == "0":
                value.append(ref)
            else:
                value.append(alt)
        if vcfDict.get(key) is not None:
            vcfDict[key] = combineRows(ref, vcfDict.get(key), value)                
        else:
            vcfDict[key] = value
f.close()

print(len(vcfDict.keys()))

#==============================================================================
# for each value in the dictionary, replace those values in the reference genome
# and write to a file
#==============================================================================
# for small-scale testing purposes only, use "smaller.txt" which uses the first three lines of SNP matrix
#refgenome = "oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo"
#indivGenome = list(refgenome)

indivID = 0
#indivGenome = list(refgenome.tostring()) # needs to be of type list to put in mutations
NUMBEROFINDIVIDUALS = 40 # for future use, please change from this to length(value) in vcfDict
while indivID < 40:
    indivGenome = list(str(refgenome)) 
    
    for snp in vcfDict.keys():
        indivGenome[int(snp)-1] = vcfDict[snp][indivID]
    
    f = open("indiv_"+str(indivID)+".txt", "w+")
    f.write("".join(indivGenome))
    f.close()
    
    indivID += 1



