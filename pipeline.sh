gunzip resources/X.vcf.gz
gunzip resources/C_elegans_dna.fa.gz
run getChromosome('X', 'resources/C_elegans_dna.fa')
mkdir chrX
mv resources/chrX.fa chrX/
mv resources/X.vcf chrX/
python3 MAKEINDIVIDUALS.py resources/chrX/chrX.fa resources/chrX/X.vcf
mv indiv_*.txt resources/chrX/indivGenotypes
run create_seq_fasta('indivGenotypes')
rm -rf resources/chrX/indivGenotypes/
mkdir resources/chrX/indivGenotypes
mv resources/chrX/indiv* resources/chrX/indivGenotypes/
cp cr_read.py resources/chrX/
cd resources/chrX
mkdir reads
python3 cr_read.py 1
mkdir reads_fa
run convert_to_fasta()
python3 cr_read.py 2
cp ../../kmerIndivDict.py .
python3 kmerIndivDict.py
mkdir kmer_contigs
mv reads_fa/*.txt kmer_contigs/
run SNPFormation.ipynb