import os
import glob

def create_reads(ref_dir):
    for filename in glob.glob(ref_dir + '/*.fa'):
        outfileprefix = filename.split('/')[1].split('.')[0]
        cmd1 = 'art_illumina -ss HS25 -i ' + filename + ' -o reads/' + outfileprefix + ' -l 150 -f 5 -k 0'
        #os.system(cmd1)
        cmd2 = 'jellyfish count -C -m 20 -s 100M -t 8 reads/' + outfileprefix + '.fq'
        #os.system(cmd2)
        cmd3 = 'jellyfish dump mer_counts.jf -L 2 -o kmers/' + outfileprefix + '.fa'
        #os.system(cmd3)
        cmd4 = 'rm reads/' + outfileprefix + '.aln'
        cmd5 = 'rm reads/' + outfileprefix + '.fq'
        cmd6 = 'rm mer_counts.jf'
        #os.system(cmd4)
        #os.system(cmd5)
        #os.system(cmd6)


create_reads('indivGenotypes')
