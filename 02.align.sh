# align with hmmer

# multiple alignment using muscle
module load muscle/5.1
muscle -align 40_refs.fa -output 40_refs.afa

# running hmmer
module load hmmer/3.3.2

# build reference
hmmbuild 40_refs.hmm 40_refs.afa

# DNA similarity search
nhmmer --notextw --tblout nhmmer_q0.res.tab 40_refs.hmm all_inserts_q0.fasta > nhmmer_q0.res.txt

# align with blastn
module load blast/2.15.0+
makeblastdb -in all_inserts_q0.fasta -dbtype nucl
blastn -db all_inserts_q0.fasta -query FMR1_BisTemplate.fa -out ref_to_all_inserts_q0.txt -line_length 1000 -max_target_seqs 5000
