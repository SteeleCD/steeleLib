#!/bin/bash
# alignONB
# directory info
dataDir="/media/DATA1/ONB_project_exome/chris"
softDir="/media/DATA1/Analysis_tools"
refDir="/media/DATA1/Analysis_tools/GenomeAnalysisTK-3.3/resources/human_g1k_v37"
#refFile="$refDir/human_g1k_v37_decoy.fasta"
refFile="$refDir/human_g1k_v37.fasta"
indelsFile="$softDir/GenomeAnalysisTK-3.3/resources/Mills_and_1000G_gold_standard.indels.b37.vcf"
snpFile="$softDir/GenomeAnalysisTK-3.3/resources/1000G_phase1.snps.high_confidence.b37.vcf"
indelsFile2="$softDir/GenomeAnalysisTK-3.3/resources/1000G_phase1.indels.b37.vcf"
intervalList="$softDir/GenomeAnalysisTK-3.3/resources/Indel.intervals"
captureIntervals="$dataDir/SeqCap_exome.intervals"
# samples
samples=$(ls $dataDir --ignore=*.sh --ignore=*.R --ignore=*.interval --ignore=docker --ignore=sequenzaOut --ignore=test --ignore=wget-log --ignore *.intervals --ignore FQC)
echo $samples

# generate realignment intervals
# java -Xmx1g -jar $softDir/GATK-3.5/GenomeAnalysisTK.jar \
#  -T RealignerTargetCreator \
#  -R $refFile \
#  -o $outDir/$samp/output.intervals \
#tmp=$(echo $samples|awk '{print $18 " " $24}') # subset to two samples
tmp=$(echo $samples|awk '{print $8}')
samples=$tmp
echo $samples 
for samp in $samples
do
	outDir="$dataDir/$samp"
	echo $outDir
	# get read group info
	 declare $( awk -F, 'FNR == 2 {print "id="$3}' $outDir/SampleSheet.csv )
	 declare $( awk -F, 'FNR == 2 {print "sm="$4}' $outDir/SampleSheet.csv )
	 declare $( awk -F, 'FNR == 2 {print "lb="$8}' $outDir/SampleSheet.csv )
	 declare $( awk -F, 'FNR == 2 {print "	pu="$2}' $outDir/SampleSheet.csv )
	# get fastQs
	fastQs=$(ls $outDir/*.fastq.gz)
	echo $fastQs
	set -- $fastQs
	# generate SAM with bwa mem
	bwa mem -aM -R "@RG\tID:"$id"\tSM:"$sm"\tPL:ILLUMINA\tLB:"$lb"\tPU:"$pu $refFile $1 $2 > $outDir/aligned.sam
	# SAM -> BAM
	samtools view -bS  $outDir/aligned.sam > $outDir/aligned.bam
	# sort BAM
	samtools sort -o $outDir/sorted.bam $outDir/aligned.bam 
	# index BAM
	samtools index $outDir/sorted.bam $outDir/sorted.bai
	# mark duplicates with picard
	java -jar $softDir/picard-tools-1.128/picard.jar MarkDuplicates \
	INPUT=$outDir/sorted.bam \
	OUTPUT=$outDir/remDups.bam \
	METRICS_FILE=$outDir/duplicatesMetrics.txt \
	REMOVE_DUPLICATES=true
	# reindex
	samtools index $outDir/remDups.bam $outDir/remDups.bai
	# realign BAM around indels
	java -Xmx4g -Djava.io.tmpdir=/tmp \
  	-jar $softDir/GATK-3.5/GenomeAnalysisTK.jar \
	--filter_bases_not_stored \
  	-I $outDir/remDups.bam \
  	-R $refFile \
  	-T IndelRealigner \
  	-targetIntervals $intervalList \
  	-o $outDir/realigned.bam \
  	-known $indelsFile \
	-known $indelsFile2
	# base recalibration info
	java -jar $softDir/GATK-3.5/GenomeAnalysisTK.jar \
   	-T BaseRecalibrator \
   	-R $refFile \
   	-I $outDir/realigned.bam \
   	-knownSites $indelsFile \
	-knownSites $snpFile \
   	-o $outDir/recal_data.table #\
	#-L $captureIntervals
	# base quality score recalibrator
	java -jar $softDir/GATK-3.5/GenomeAnalysisTK.jar \
   	-T PrintReads \
   	-R $refFile \
   	-I $outDir/realigned.bam \
   	-BQSR $outDir/recal_data.table \
   	-o $outDir/recalibrated.bam
	# reindex
	samtools index $outDir/recalibrated.bam $outDir/recalibrated.bai
	# base quality plots
	java -jar $softDir/GATK-3.5/GenomeAnalysisTK.jar \
	-T BaseRecalibrator \
	-R $refFile \
	-I $outDir/realigned.bam \
	-knownSites $indelsFile \
	-knownSites $snpFile \
	-BQSR $outDir/recal_data.table \
	-o $outDir/recal_data2.table #\
	#-L $captureIntervals
	java -jar $softDir/GATK-3.5/GenomeAnalysisTK.jar \
	-T AnalyzeCovariates \
	-R $refFile \
	-before $outDir/recal_data.table \
	-after $outDir/recal_data2.table \
	-plots $outDir/recal_plots.pdf 
	# remove intermediates
	if [ ! -f $outDir/aligned.sam ]; then
		echo "Failed at fastQ->aligned.sam"
		echo $outDir
		exit
	fi
	if [ ! -f $outDir/aligned.bam ]; then
                echo "Failed at aligned.sam->aligned.bam"
		echo $outDir
		exit
        fi
	rm $outDir/aligned.sam
	if [ ! -f $outDir/sorted.bam ]; then
                echo "Failed at aligned.bam->sorted.bam"
		echo $outDir
                exit
        fi
	rm $outDir/aligned.bam
	if [ ! -f $outDir/sorted.bai ]; then
                echo "Failed at sorted.bam->sorted.bai"
		echo $outDir
                exit
        fi
	if [ ! -f $outDir/remDups.bam ]; then
                echo "Failed at sorted.bam->remDups.bam"
		echo $outDir
                exit
        fi
	rm $outDir/sorted.bam
	rm $outDir/sorted.bai
	if [ ! -f $outDir/remDups.bai ]; then
                echo "Failed at remDups.bam->remDups.bai"
		echo $outDir
                exit
        fi
        if [ ! -f $outDir/realigned.bam ]; then
                echo "Failed at remDups.bam->realigned.bam"
		echo $outDir
                exit
        fi
	if [ ! -f $outDir/realigned.bai ]; then
                echo "Failed at realigned.bam->realigned.bai"
		echo $outDir
                exit
        fi
	rm $outDir/remDups.bam
	rm $outDir/remDups.bai
	rm $outDir/realigned.bam
	rm $outDir/realigned.bai
done




