#Create pooled spacer sets

### For phage
# 1. Combine C1 files
cat FASTAmaps/B+P_b_+1_C1_mapped_fcl2_C1.fastq FASTAmaps/B+P_e_+1_C1_mapped_fcl2_C1.fastq > Pooled_spacers/C1_phage_pool.fastq

# 2. Combine C2 files
cat FASTAmaps/B+P_b_+1_C2_mapped_fcl2_C2.fastq FASTAmaps/B+P_d_+1_C2_mapped_fcl2_C2.fastq  FASTAmaps/B+P_e_+1_C2_mapped_fcl2_C2.fastq > Pooled_spacers/C2_phage_pool.fastq

# 3. Cluster the combined C1 and C2 files (separately) to get unique spacers from the pools. Threshold 0.8 is still good to avoid 100% overlaps with snps
#cd-hit-est -i Pooled_spacers/C1_phage_pool.fastq -o Pooled_spacers/C1_phage_pool_clustered.fastq -c 0.8 -n 5 -M 16000 -T 8
#cd-hit-est -i Pooled_spacers/C2_phage_pool.fastq -o Pooled_spacers/C2_phage_pool_clustered.fastq -c 0.8 -n 5 -M 16000 -T 8

cp Pooled_spacers/C1_phage_pool.fastq Pooled_spacers/C1_phage_pool_clustered.fastq
cp Pooled_spacers/C2_phage_pool.fastq Pooled_spacers/C2_phage_pool_clustered.fastq

### For self

# 1. Combine C1 files
cat FASTAmaps/B+P_b_+1_C1_mapped_b185_C1.fastq FASTAmaps/B+P_e_+1_C1_mapped_b185_C1.fastq > Pooled_spacers/C1_b185_pool.fastq

# 2. Combine C2 files
cat FASTAmaps/B+P_b_+1_C2_mapped_b185_C2.fastq FASTAmaps/B+P_d_+1_C2_mapped_b185_C2.fastq  FASTAmaps/B+P_e_+1_C2_mapped_b185_C2.fastq > Pooled_spacers/C2_b185_pool.fastq

# 3. Cluster the combined C1 and C2 files (separately) to get unique spacers from the pools. Threshold 0.8 is still good to avoid 100% overlaps with snps
#cd-hit-est -i Pooled_spacers/C1_b185_pool.fastq -o Pooled_spacers/C1_b185_pool_clustered.fastq -c 0.8 -n 5 -M 16000 -T 8
#cd-hit-est -i Pooled_spacers/C2_b185_pool.fastq -o Pooled_spacers/C2_b185_pool_clustered.fastq -c 0.8 -n 5 -M 16000 -T 8

cp Pooled_spacers/C1_b185_pool.fastq Pooled_spacers/C1_b185_pool_clustered.fastq
cp Pooled_spacers/C2_b185_pool.fastq Pooled_spacers/C2_b185_pool_clustered.fastq
