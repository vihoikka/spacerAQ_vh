#Ville Hoikkala 2020 (GNU GENERAL PUBLIC LICENSE v3, 29 June 2007)

#To make this work, make sure you have:
# - maintained the folder structure
# - python3
# - added data files in folder raw_reads
# - the following dependencies: bowtie2, cd-hit, samtools, weblogo and Trimmomatic (requires trimmomatic-0.36.jar in root folder)
# - Bowtie2 indexes for self and phage genomes

declare -a reads=("B+P_b_+1_C1" #Read files
                "B+P_e_+1_C1"
                "B+P_b_+1_C2"
                "B+P_d_+1_C2"
                "B+P_e_+1_C2"
                )

declare -a loci=("C1" #Determine the locus of each file above. Retain same order
            	"C1"
                "C2"
                "C2"
                "C2"
                )

#Prepare output files
echo "sample,locus,phage,self,other" > targets.csv
echo "sample,locus,target,canonical_PAM,4bp_bridge_PAM,6bp_bridge_PAM,seven_PAM,total" > PAMs/PAMs.csv
echo "sample,locus,target,f_hits,r_hits,total_hits,mRNA_hits,mRNA:DNA_ratio,intergenics" > ORFs/ORF_spacer_stats.csv

inputDir="raw_reads"
TrimmomaticOutputDir="TRIMMOMATIC"
SPACERSdir="SPACERS"
SAMdir="SAM"
FASTAdir="FASTAmaps"
spacersOutputFile=${inputFile}_spacers
clusteredSpacers=${inputFile}_spacers_clustered
targetgenome1="fcl2"
targetgenome2="b185"
genome1ORFs="FCL2_ORFs"
genome2ORFs="b185_ORFs"
#locus="C2"
binsize=3000
numberOfBins=0
binFraction=0.03
PAMsize=15
lociCounter=0

# canPAMrev="......TTTA....."
# fourbpBridgePAMrev=".......TTTA...."
# sixbpBridgePAMrev=".....TTTA......"
# sevenPAMrev="....TTTA......."
#
# canPAM=".....TAAA......"
# fourbpBridgePAM="....TAAA......."
# sixbpBridgePAM="......TAAA....."
# sevenPAM=".......TAAA...."

canPAMrev="......TTTA....."
fourbpBridgePAMrev=".......TTTA...."
sixbpBridgePAMrev=".....TTTA......"
sevenPAMrev=".....TTTA......."

canPAM=".....TAAA......"
fourbpBridgePAM="....TAAA......."
sixbpBridgePAM="......TAAA....."
sevenPAM=".......TAAA...."

for read in "${reads[@]}"
do
	inputFile=$read
	spacersOutputFile=${inputFile}_spacers
	clusteredSpacers=${inputFile}_spacers_clustered
	locus=${loci[${lociCounter}]}

	#Dispose bad quality reads with Trimmomatic
	java -jar ../trimmomatic-0.36.jar SE -phred33 ${inputDir}/$inputFile.fastq ${TrimmomaticOutputDir}/${inputFile}_trimmed.fastq SLIDINGWINDOW:3:21 MINLEN:100 TRAILING:23

	#Extract spacers from the reads with spExt and store as fastq. Discard command removes existing wt-spacers. Nodup prevents from removing duplicates (use dup to remove them)
	python ../python/spExtSets_0319.py ${TrimmomaticOutputDir}/${inputFile}_trimmed.fastq ${SPACERSdir}/${spacersOutputFile}_${locus} $locus nodup discard

	#Remove duplicates, i.e. cluster similar spacers using CD-hit
	cd-hit-est -i $SPACERSdir/${spacersOutputFile}_${locus}.fastq -o $SPACERSdir/${clusteredSpacers}_${locus}.fastq -c 0.8 -n 5 -M 16000 -T 8

  #Also make a fasta file from clustered spacers for later processing with R
  seqtk seq -a $SPACERSdir/${clusteredSpacers}_${locus}.fastq > $SPACERSdir/${clusteredSpacers}_${locus}.fasta

  #cd-hit-dup -i $SPACERSdir/${spacersOutputFile}_${locus}.fastq -o $SPACERSdir/${clusteredSpacers}_${locus}.fastq -c 0.8 -n 5 -M 16000 -T 8

  #uncomment the line below to retain duplicates
	#cp ${SPACERSdir}/${spacersOutputFile}_${locus}.fasta ${SPACERSdir}/${clusteredSpacers}_${locus}.fasta

	#Count the number of spacers
	#uniqueSpacers=$(grep -c "@" ${SPACERSdir}/${clusteredSpacers}_${locus}.fastq)
  #uniqueSpacers=$(grep -c "@" ${SPACERSdir}/${clusteredSpacers}_${locus}.fastq)
  uniqueSpacers=$(echo $(cat ${SPACERSdir}/${clusteredSpacers}_${locus}.fastq|wc -l)/4|bc)


	#Map the spacers onto first target genome and output SAM file
	echo !!! Running Bowtie2...
	bowtie2 -x $targetgenome1 -U ${SPACERSdir}/${clusteredSpacers}_${locus}.fastq -S ${SAMdir}/${inputFile}_map_${targetgenome1}_${locus}.sam -f -q --very-sensitive-local --score-min G,10,8

  echo Counting spacers on first genome...
	#Count number of hits on first genome
	genome1spacers=$(awk '$2==0 || $2==16' ${SAMdir}/${inputFile}_map_${targetgenome1}_${locus}.sam | wc -l)

  echo !!! Extracting non-mapped reads from SAM
	#This extracts all non-mapped reads from the SAM-file
	awk '($2 != 0) && ($2 != 16)' ${SAMdir}/${inputFile}_map_${targetgenome1}_${locus}.sam > ${SAMdir}/${inputFile}_unmapped_${targetgenome1}_${locus}.sam
  samtools bam2fq ${SAMdir}/${inputFile}_unmapped_${targetgenome1}_${locus}.sam > ${FASTAdir}/${inputFile}_unmapped_${targetgenome1}_${locus}.fastq

  #Make a fastq from genome1-mapped spacers
  awk '($2==0 || $2==16 || NR < 5)' ${SAMdir}/${inputFile}_map_${targetgenome1}_${locus}.sam > ${SAMdir}/${inputFile}_mapped_${targetgenome1}_${locus}.sam
  #Transform the genome1 -mapped reads to fastq
  head ${SAMdir}/${inputFile}_mapped_${targetgenome1}_${locus}.sam
	samtools bam2fq ${SAMdir}/${inputFile}_mapped_${targetgenome1}_${locus}.sam > ${FASTAdir}/${inputFile}_mapped_${targetgenome1}_${locus}.fastq


	#Map the unmapped spacers from previous step onto second target genome and output SAM file
	echo Running Bowtie2
	bowtie2 -x $targetgenome2 -U ${FASTAdir}/${inputFile}_unmapped_${targetgenome1}_${locus}.fastq -S ${SAMdir}/${inputFile}_map_${targetgenome2}_${locus}.sam -f --very-sensitive-local -q -N 1

	#Count number of hits on second genome
	genome2spacers=$(awk '$2==0 || $2==16' ${SAMdir}/${inputFile}_map_${targetgenome2}_${locus}.sam | wc -l)

	#This extracts all non-mapped reads from the SAM-file
	awk '($2 != 0) && ($2 != 16)' ${SAMdir}/${inputFile}_map_${targetgenome2}_${locus}.sam > ${SAMdir}/${inputFile}_unmapped_final_${locus}.sam
  #Transform the final unmapped SAM-reads to fasta
	samtools bam2fq ${SAMdir}/${inputFile}_unmapped_final_${locus}.sam > ${FASTAdir}/${inputFile}_unmapped_final_${locus}.fastq

  #Make a fastq from genome1-mapped spacers
  awk '($2==0 || $2==16 || NR < 5)' ${SAMdir}/${inputFile}_map_${targetgenome2}_${locus}.sam > ${SAMdir}/${inputFile}_mapped_${targetgenome2}_${locus}.sam
  #Transform the genome2 -mapped reads to fastq
	samtools bam2fq ${SAMdir}/${inputFile}_mapped_${targetgenome2}_${locus}.sam > ${FASTAdir}/${inputFile}_mapped_${targetgenome2}_${locus}.fastq



	dmSpacers=$(echo $(cat ${FASTAdir}/${inputFile}_unmapped_final_${locus}.fastq|wc -l)/4|bc)

	#Write target results to csv file
	echo "${read},${locus},${genome1spacers},${genome2spacers},${dmSpacers}" >> targets.csv

	#Bin the hits from SAM files with spBinner and store results as .csv
	python ../python/spBinner.py ${SAMdir}/${inputFile}_map_${targetgenome1}_${locus}.sam $binsize 0 percent Binned/${inputFile}_map_${targetgenome1}_binned_${locus} ${targetgenome1}.fasta $numberOfBins $binFraction

	#Bin the hits from SAM files with spBinner and store results as .csv
	python ../python/spBinner.py ${SAMdir}/${inputFile}_map_${targetgenome2}_${locus}.sam $binsize 0 percent Binned/${inputFile}_map_${targetgenome2}_binned_${locus} ${targetgenome2}.fasta $numberOfBins $binFraction

	#locate and extract PAM sequences with pamSeek (phage)
	python ../python/pamSeek.py ${SAMdir}/${inputFile}_map_${targetgenome1}_${locus}.sam fcl2 PAMs/fasta/$inputFile 0 $PAMsize

	#locate and extract PAM sequences with pamSeek (self)
	python ../python/pamSeek.py ${SAMdir}/${inputFile}_map_${targetgenome2}_${locus}.sam b185 PAMs/fasta/$inputFile 0 $PAMsize

	#Fetch PAM data and write to csv file and weblogo PNGs

	if [ $locus == C1 ]
	then

	totalPAMcount=$(grep -c ">" PAMs/fasta/${inputFile}_PAMs_downstream_fcl2.fasta)
	canPAMcount=$(grep -c $canPAM PAMs/fasta/${inputFile}_PAMs_downstream_fcl2.fasta)
	sevenPAMcount=$(grep -c $sevenPAM PAMs/fasta/${inputFile}_PAMs_downstream_fcl2.fasta)
	fourbpbridgePAMcount=$(grep -c $fourbpBridgePAM PAMs/fasta/${inputFile}_PAMs_downstream_fcl2.fasta)
	sixbpBridgePAMcount=$(grep -c $sixbpBridgePAM PAMs/fasta/${inputFile}_PAMs_downstream_fcl2.fasta)

		echo "${inputFile},${locus},${targetgenome1},${canPAMcount},${fourbpbridgePAMcount},${sixbpBridgePAMcount},${sevenPAMcount},${totalPAMcount}" >> PAMs/PAMs.csv

	totalPAMcount=$(grep -c ">" PAMs/fasta/${inputFile}_PAMs_downstream_b185.fasta)
	canPAMcount=$(grep -c $canPAM PAMs/fasta/${inputFile}_PAMs_downstream_b185.fasta)
	sevenPAMcount=$(grep -c $sevenPAM PAMs/fasta/${inputFile}_PAMs_downstream_b185.fasta)
	fourbpbridgePAMcount=$(grep -c $fourbpBridgePAM PAMs/fasta/${inputFile}_PAMs_downstream_b185.fasta)
	sixbpBridgePAMcount=$(grep -c $sixbpBridgePAM PAMs/fasta/${inputFile}_PAMs_downstream_b185.fasta)

	echo "${inputFile},${locus},${targetgenome2},${canPAMcount},${fourbpbridgePAMcount},${sixbpBridgePAMcount},${sevenPAMcount},${totalPAMcount}" >> PAMs/PAMs.csv
	fi


	if [ $locus == C2 ]
	then

	totalPAMcount=$(grep -c ">" PAMs/fasta/${inputFile}_PAMs_upstream_fcl2.fasta)
	canPAMcount=$(grep -c $canPAMrev PAMs/fasta/${inputFile}_PAMs_upstream_fcl2.fasta)
	sevenPAMcount=$(grep -c $sevenPAMrev PAMs/fasta/${inputFile}_PAMs_upstream_fcl2.fasta)
	fourbpbridgePAMcount=$(grep -c $fourbpBridgePAMrev PAMs/fasta/${inputFile}_PAMs_upstream_fcl2.fasta)
	sixbpBridgePAMcount=$(grep -c $sixbpBridgePAMrev PAMs/fasta/${inputFile}_PAMs_upstream_fcl2.fasta)

		echo "${inputFile},${locus},${targetgenome1},${canPAMcount},${fourbpbridgePAMcount},${sixbpBridgePAMcount},${sevenPAMcount},${totalPAMcount}" >> PAMs/PAMs.csv

	totalPAMcount=$(grep -c ">" PAMs/fasta/${inputFile}_PAMs_upstream_b185.fasta)
	canPAMcount=$(grep -c $canPAMrev PAMs/fasta/${inputFile}_PAMs_upstream_b185.fasta)
	sevenPAMcount=$(grep -c $sevenPAMrev PAMs/fasta/${inputFile}_PAMs_upstream_b185.fasta)
	fourbpbridgePAMcount=$(grep -c $fourbpBridgePAMrev PAMs/fasta/${inputFile}_PAMs_upstream_b185.fasta)
	sixbpBridgePAMcount=$(grep -c $sixbpBridgePAMrev PAMs/fasta/${inputFile}_PAMs_upstream_b185.fasta)

		echo "${inputFile},${locus},${targetgenome2},${canPAMcount},${fourbpbridgePAMcount},${sixbpBridgePAMcount},${sevenPAMcount},${totalPAMcount}" >> PAMs/PAMs.csv

	fi



	weblogo --units 'bits' --size large --format png_print --title "${inputFile} PAMs upstream ${targetgenome1}" < PAMs/fasta/${inputFile}_PAMs_upstream_fcl2.fasta > PAMs/${inputFile}_PAMs_upstream_${targetgenome1}.png || exit

	weblogo --units 'bits' --size large --format png_print --title "${inputFile} PAMs downstream ${targetgenome1}" < PAMs/fasta/${inputFile}_PAMs_downstream_fcl2.fasta > PAMs/${inputFile}_PAMs_downstream_${targetgenome1}.png || exit

	weblogo --units 'bits' --size large --format png_print --title "${inputFile} PAMs upstream ${targetgenome2}" < PAMs/fasta/${inputFile}_PAMs_upstream_b185.fasta > PAMs/${inputFile}_PAMs_upstream_${targetgenome2}.png || exit

	weblogo --units 'bits' --size large --format png_print --title "${inputFile} PAMs downstream ${targetgenome2}" < PAMs/fasta/${inputFile}_PAMs_downstream_b185.fasta > PAMs/${inputFile}_PAMs_downstream_${targetgenome2}.png || exit

	#Characterize ORF targeting with ORFmapper for genome 1
	python ../python/ORFs_spacers.py ${SAMdir}/${inputFile}_map_${targetgenome1}_${locus}.sam ${genome1ORFs}.csv ORFs/${inputFile}_ORFs_${genome1ORFs}_${locus}.csv

	#Characterize ORF targeting with ORFmapper for genome 2
	python ../python/ORFs_spacers.py ${SAMdir}/${inputFile}_map_${targetgenome2}_${locus}.sam ${genome2ORFs}.csv ORFs/${inputFile}_ORFs_${genome2ORFs}_${locus}.csv

	#Stats about spacers on ORFS for both target genomes. Mine the spacer orf csv files.

	f_hits=$(awk -F, '{print $2}' ORFs/${inputFile}_ORFs_${genome1ORFs}_${locus}.csv_ORF_spacer_stats.csv | grep -c 'f')

	r_hits=$(awk -F, '{print $2}' ORFs/${inputFile}_ORFs_${genome1ORFs}_${locus}.csv_ORF_spacer_stats.csv | grep -c 'r')

	total_hits=$(grep -c "" ORFs/${inputFile}_ORFs_${genome1ORFs}_${locus}.csv_ORF_spacer_stats.csv)

	mRNA_hits=$(awk -F, '{print $5}' ORFs/${inputFile}_ORFs_${genome1ORFs}_${locus}.csv_ORF_spacer_stats.csv | grep -c '1')

	mRNA_ratio=$((100*$mRNA_hits/$total_hits))

	intergenics=$(awk -F, '{print $6}' ORFs/${inputFile}_ORFs_${genome1ORFs}_${locus}.csv_ORF_spacer_stats.csv | grep -c 'yes')

	echo "${inputFile},${locus},${targetgenome1},${f_hits},${r_hits},${total_hits},${mRNA_hits},${mRNA_ratio},${intergenics}" >> ORFs/ORF_spacer_stats.csv


	#Repeat for other genome


	f_hits=$(awk -F, '{print $2}' ORFs/${inputFile}_ORFs_${genome2ORFs}_${locus}.csv_ORF_spacer_stats.csv | grep -c 'f')

	r_hits=$(awk -F, '{print $2}' ORFs/${inputFile}_ORFs_${genome2ORFs}_${locus}.csv_ORF_spacer_stats.csv | grep -c 'r')

	total_hits=$(grep -c "" ORFs/${inputFile}_ORFs_${genome2ORFs}_${locus}.csv_ORF_spacer_stats.csv)

	mRNA_hits=$(awk -F, '{print $5}' ORFs/${inputFile}_ORFs_${genome2ORFs}_${locus}.csv_ORF_spacer_stats.csv | grep -c '1')

	mRNA_ratio=$((100*$mRNA_hits/$total_hits))

	intergenics=$(awk -F, '{print $6}' ORFs/${inputFile}_ORFs_${genome2ORFs}_${locus}.csv_ORF_spacer_stats.csv | grep -c 'yes')

	echo "${inputFile},${locus},${targetgenome2},${f_hits},${r_hits},${total_hits},${mRNA_hits},${mRNA_ratio},${intergenics}" >> ORFs/ORF_spacer_stats.csv


	g1Per=$(echo "scale=3; $genome1spacers / $uniqueSpacers * 100" | bc)
	g2Per=$(echo "scale=3; $genome2spacers / $uniqueSpacers * 100" | bc)
	dmPer=$(echo "scale=3; $dmSpacers / $uniqueSpacers * 100" | bc)

	echo "Unique spacers: $uniqueSpacers"
	echo "Spacers mapped to ${targetgenome1}:${genome1spacers} (${g1Per}%)"
	echo "Spacers mapped to ${targetgenome2}:${genome2spacers} (${g2Per}%)"
	echo "Spacers mapped to nothing: ${dmSpacers} (${dmPer}%)"

	lociCounter=$((lociCounter+1))
done

#Calculate number of shared spacers
sh Commonspacers/commonspacers_cdhit2d_pooler.sh

#Simulate spacer acquisition on phage genome in random positions or next to PAMs. This takes a long time to run and is therefore commented
#sh spacersimulation_locations/spacersim_v2_131219.sh

#Create spacer pools
sh poolSpacers_171219.sh

#Perform mRNA targeting analysis on pooled spacers
sh pooled_spacers_ORFs.sh
