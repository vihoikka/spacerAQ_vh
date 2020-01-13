declare -a reads=("C1_phage_pool_clustered.fastq" #Pooled files
                "C2_phage_pool_clustered.fastq"
                "C1_b185_pool_clustered.fastq" #Pooled files
                "C2_b185_pool_clustered.fastq"
                )

declare -a loci=("C1" #Determine the locus of each file above. Retain same order
                "C2"
                "C1"
                "C2")

declare -a targetgenomes=("fcl2" #Determine the locus of each file above. Retain same order
                "fcl2"
                "b185"
                "b185")

targetgenome1="fcl2"
targetgenome2="b185"
SAMdir="Pooled_spacers/SAM"
FASTAdir="Pooled_spacers/FASTAmaps"
lociCounter=0
genomeCounter=0

#Prepare file for ORF targeting stats
echo "sample,locus,target,f_hits,r_hits,total_hits,mRNA_hits,mRNA:DNA_ratio,intergenics" > Pooled_spacers/ORF_targets/ORF_spacer_stats_pooled.csv


for read in "${reads[@]}"
do
  inputFile=$read
	locus=${loci[${lociCounter}]}
  targetgenome=${targetgenomes[$genomeCounter]}

  echo !!! Running Bowtie2...
  bowtie2 -x $targetgenome -U Pooled_spacers/$inputFile -S ${SAMdir}/${inputFile}_map_${targetgenome}_${locus}.sam -f -q --very-sensitive-local --score-min G,10,8

  #Characterize ORF targeting with ORFmapper for genome 1
	python ../python/ORFs_spacers.py ${SAMdir}/${inputFile}_map_${targetgenome}_${locus}.sam ${targetgenome}_ORFs.csv Pooled_spacers/ORF_targets/${inputFile}_${targetgenome}_ORFs_${locus}.csv

  f_hits=$(awk -F, '{print $2}' Pooled_spacers/ORF_targets/${inputFile}_${targetgenome}_ORFs_${locus}.csv_ORF_spacer_stats.csv | grep -c 'f')

	r_hits=$(awk -F, '{print $2}' Pooled_spacers/ORF_targets/${inputFile}_${targetgenome}_ORFs_${locus}.csv_ORF_spacer_stats.csv | grep -c 'r')

	total_hits=$(grep -c "" Pooled_spacers/ORF_targets/${inputFile}_${targetgenome}_ORFs_${locus}.csv_ORF_spacer_stats.csv)

	mRNA_hits=$(awk -F, '{print $5}' Pooled_spacers/ORF_targets/${inputFile}_${targetgenome}_ORFs_${locus}.csv_ORF_spacer_stats.csv | grep -c '1')

	mRNA_ratio=$((100*$mRNA_hits/$total_hits))

	intergenics=$(awk -F, '{print $6}' Pooled_spacers/ORF_targets/${inputFile}_${targetgenome}_ORFs_${locus}.csv_ORF_spacer_stats.csv | grep -c 'yes')

	echo "${inputFile},${locus},${targetgenome},${f_hits},${r_hits},${total_hits},${mRNA_hits},${mRNA_ratio},${intergenics}" >> Pooled_spacers/ORF_targets/ORF_spacer_stats_pooled.csv

	#Write target results to csv file
#	echo "${read},${locus},${genome1spacers},${genome2spacers},${dmSpacers}" >> targets.csv

  lociCounter=$((lociCounter+1))
  genomeCounter=$((genomeCounter+1))


done
