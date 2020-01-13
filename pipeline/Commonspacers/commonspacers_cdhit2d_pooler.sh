#First, pool all spacer within a locus together (i.e concatenate replicates)

#Define file names
c1_pool=Pooled_spacers/C1_phage_pool.fastq
c2_pool=Pooled_spacers/C2_phage_pool.fastq
c1_pool_cluster=Pooled_spacers/C1_phage_pool_clustered.fastq
c2_pool_cluster=Pooled_spacers/C2_phage_pool_clustered.fastq
allPooled=Pooled_spacers/AllPooled_cluster.fastq

# 1. Combine C1 files
cat FASTAmaps/B+P_b_+1_C1_mapped_fcl2_C1.fastq FASTAmaps/B+P_e_+1_C1_mapped_fcl2_C1.fastq > $c1_pool

# 2. Combine C2 files
cat FASTAmaps/B+P_b_+1_C2_mapped_fcl2_C2.fastq FASTAmaps/B+P_d_+1_C2_mapped_fcl2_C2.fastq FASTAmaps/B+P_e_+1_C2_mapped_fcl2_C2.fastq > $c2_pool

# 3. Cluster the combined C1 and C2 files (separately) to get unique spacers from the pools
cd-hit-est -i $c1_pool -o $c1_pool_cluster -c 0.8 -n 5 -M 16000 -T 8
cd-hit-est -i $c2_pool -o $c2_pool_cluster -c 0.8 -n 5 -M 16000 -T 8



c1=$(echo $(cat $c1_pool_cluster|wc -l)/4|bc)
echo "II-C clustered spacer pool size: $c1"

c2=$(echo $(cat $c2_pool_cluster|wc -l)/4|bc)
echo "VI-B clustered spacer pool size: $c2"

cd-hit-est-2d -i $c1_pool_cluster -i2 $c2_pool_cluster -o $allPooled -c 0.9 -n 5 -M 16000 -T 8

commonSpacers=$(grep ">Cluster" $allPooled.clstr -c)
echo "Cluster entries: $commonSpacers"
echo "Total VI-B on FCL2 spacers: $c2"
#uniqueToC2_b=$(grep "@" 2d_b.fastq -c)
uniqueToC2=$(echo $(cat $allPooled|wc -l)/4|bc)
echo "Spacers unique to VI-B: $uniqueToC2"
shared=$(echo "$c2 - $uniqueToC2" | bc)
echo "Spacers shared between loci: $shared"


prop_c2=$(echo "scale=2 ; 100 - $uniqueToC2 / $c2 * 100" | bc) #Ratio of unique VI-B spacers and shared spacers
echo "Proportion of shared VI-B spacers: $prop_c2%"
echo "-------"
