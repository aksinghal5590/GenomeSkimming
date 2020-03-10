import os, sys
import glob

def create_reads(ref_dir):
	for filename in glob.glob(ref_dir + '/*.fa'):
		outfileprefix = filename.split('/')[1].split('.')[0]
		cmd1 = 'art_illumina -ss HS25 -i ' + filename + ' -o reads/' + outfileprefix + ' -l 150 -f 5 -k 0'
		os.system(cmd1)
		cmd2 = 'rm reads/' + outfileprefix + '.aln'
		os.system(cmd2)


def exec_jellyfish():
	cmd1 = 'jellyfish count -C -m 25 -s 100M -t 8 reads_fa/*.fa'
	cmd2 = 'jellyfish dump mer_counts.jf -L 20 -o kmers.fa'
	os.system(cmd1)
	os.system(cmd2)


print(int(sys.argv[1]))
if int(sys.argv[1]) == 1:
	create_reads('indivGenotypes')
if int(sys.argv[1]) == 2:
	exec_jellyfish()
