#For simulating spacer distribution on the phage genome

declare -a clusteredList
cycleCounter=0

#Define files
c1_pool_cluster=Pooled_spacers/C1_phage_pool_clustered.fastq #Actual C1 spacers to which compare simulated ones

for i in {1..1000}
do
python spacersimulation_locations/random_spacer_sampling_fcl2_v2.py --genome fcl2.fasta --mode random -c 430 -ps NNNNNTAAA -pp downstream > /dev/null 2>&1
cd-hit-est-2d -i $c1_pool_cluster -i2 randomizedSpacers.fasta -o spacersimulation_locations/clust.clstr -c 0.90 -M 16000 -T 8 > /dev/null 2>&1
#cd-hit-est-2d -i ../SAM/CorrectedSpacers/B+P_b_+1_C1_mapped_fcl2_C1.fasta -i2 randomizedSpacers.fasta -o clust.clstr -c 0.8 -M 16000 -T 8 > /dev/null 2>&1
clusters=$(grep ">" -c spacersimulation_locations/clust.clstr)
#echo $clusters
clusteredList[$cycleCounter]=$clusters

cycleCounter=$((++cycleCounter))
done

printf "${clusteredList[*]}\n" > spacersimulation_locations/random_set.txt

#Same run for PAM version of the analysis
cycleCounter=0

for i in {1..1000}
do
python spacersimulation_locations/random_spacer_sampling_fcl2_v2.py --genome fcl2.fasta --mode PAM -c 430 -ps NNNNNTAAA -pp downstream > /dev/null 2>&1
cd-hit-est-2d -i $c1_pool_cluster -i2 randomizedSpacers.fasta -o spacersimulation_locations/clust.clstr -c 0.90 -M 16000 -T 8 > /dev/null 2>&1
#cd-hit-est-2d -i ../SAM/CorrectedSpacers/B+P_b_+1_C1_mapped_fcl2_C1.fasta -i2 randomizedSpacers.fasta -o clust.clstr -c 0.8 -M 16000 -T 8 > /dev/null 2>&1
clusters=$(grep ">" -c spacersimulation_locations/clust.clstr)
#echo $clusters
clusteredList[$cycleCounter]=$clusters

cycleCounter=$((++cycleCounter))
done

printf "${clusteredList[*]}\n" > spacersimulation_locations/PAM_set.txt

#cd-hit-est-2d -i ../Commonspacers/B+P_b_+1_C1_mapped_fcl2_C1.fastq -i2 100percent_testset.fasta -o clust.clstr -c 1.0 -M 16000 -T 8
