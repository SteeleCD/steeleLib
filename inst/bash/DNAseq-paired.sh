#!/bin/bash
# align and call variants for paired tumor normal genomes/exomes

# reset POSIX
OPTIND=1

# set variables
dataDir=""
refFile=""
softDir=""
indelsFile=""
indelsFile2=""
snpFile=""
intervalList=""
captureIntervals=""
normString="N"
tumourString="T"
uniqueString1=""
uniqueString2="_"
picardDir="picard-tools-1.128"
GATKdir="GATK-3.5"

# get variables
while getopts "h?vf:" opt; do
    case "$opt" in
    h|\?)
        echo "help"
        exit 0
        ;;
    d)  dataDir=$OPTARG
        ;;
    r)  refFile=$OPTARG
        ;;
    s)  softDir=$OPTARG
        ;;
    i)  indelsFile=$OPTARG
        ;;
    I)  indelsFile2=$OPTARG
        ;;
    S)  snpFile=$OPTARG
        ;;
    l)  intervalList=$OPTARG
        ;;
    c)  captureIntervals=$OPTARG
        ;;
    n)  normString=$OPTARG
        ;;
    t)  tumourString=$OPTARG
        ;;
    u)  uniqueString1=$OPTARG
        ;;
    U)  uniqueString2=$OPTARG
        ;;
    p)  picardDir=$OPTARG
        ;;
    g)  GATKdir=$OPTARG
        ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

# quit if any variables missing 
if [ -z $dataDir ]; then
	echo "No d argument (data directory)"
	exit
fi
if [ -z $refFile ]; then
	echo "No r argument (reference fasta)"
	exit
fi
if [ -z $softDir ]; then
	echo "No s argument (software directory)"
	exit
fi
if [ -z $indelsFile ]; then
	echo "No i argument (indels vcf)"
	exit
fi
if [ -z $indelsFile2 ]; then
	echo "No I argument (indels vcf)"
	exit
fi
if [ -z $snpFile ]; then
	echo "No S argument (SNP vcf)"
	exit
fi
if [ -z $intervalList ]; then
	echo "No l argument (indel interval list)"
	exit
fi

# samples
samples=$(ls $dataDir)
echo $samples
# assume cases are N (normal) and T (tumour)
cases="$normString $tumourString"
# loop over individuals
for ind in $individuals
do
	# loop over cases (N or T)
	for case in $cases
		do
		# ====================================================================
		# get information for sample
		# ====================================================================
		outDir="$dataDir/$ind"
		echo $outDir
		# get read group info
		 declare $( awk -F, 'FNR == 2 {print "id="$3}' $outDir/SampleSheet$case.csv )
		 declare $( awk -F, 'FNR == 2 {print "sm="$4}' $outDir/SampleSheet$case.csv )
		 declare $( awk -F, 'FNR == 2 {print "lb="$8}' $outDir/SampleSheet$case.csv )
		 declare $( awk -F, 'FNR == 2 {print "	pu="$2}' $outDir/SampleSheet$case.csv )
		# get fastQs normal
		fastQs=$(ls $outDir/*$uniqueString1$case$uniqueString2*.fastq.gz)
		echo $fastQs
		set -- $fastQs
		# ====================================================================
		# generate SAM with bwa mem
		# ====================================================================
		bwa mem -aM -R "@RG\tID:"$id"\tSM:"$sm"\tPL:ILLUMINA\tLB:"$lb"\tPU:"$pu \
			 $refFile $1 $2 > $outDir/aligned$case.sam
		# check if complete
		if [ ! -f $outDir/aligned$case.sam ]; then
			echo "Failed at fastQ->aligned.sam"
			echo $outDir
			exit
		fi
		# ====================================================================
		#                    SAM -> BAM
		# ====================================================================
		samtools view -bS  $outDir/aligned$case.sam > $outDir/aligned$case.bam
		# check if complete
		if [ ! -f $outDir/aligned$case.bam ]; then
	                echo "Failed at aligned.sam->aligned.bam"
			echo $outDir
			exit
	        fi
		rm $outDir/aligned.sam
		# ====================================================================
		# sort BAM
		# ====================================================================
		samtools sort -o $outDir/sorted$case.bam $outDir/aligned$case.bam 
		# check if complete
		if [ ! -f $outDir/sorted$case.bam ]; then
	                echo "Failed at aligned.bam->sorted.bam"
			echo $outDir
	                exit
	        fi
		rm $outDir/aligned.bam
		# ====================================================================
		# index BAM
		# ====================================================================
		samtools index $outDir/sorted$case.bam $outDir/sorted$case.bai
		# check if complete
		if [ ! -f $outDir/sorted$case.bai ]; then
	                echo "Failed at sorted.bam->sorted.bai"
			echo $outDir
	                exit
	        fi
		# ====================================================================
		# mark duplicates with picard
		# ====================================================================
		java -jar $softDir/$picardDir/picard.jar MarkDuplicates \
		INPUT=$outDir/sorted$case.bam \
		OUTPUT=$outDir/remDups$case.bam \
		METRICS_FILE=$outDir/duplicatesMetrics$case.txt \
		REMOVE_DUPLICATES=true
		# check if complete
		if [ ! -f $outDir/remDups$case.bam ]; then
	                echo "Failed at sorted.bam->remDups.bam"
			echo $outDir
	                exit
	        fi
		# ====================================================================
		# reindex
		# ====================================================================
		samtools index $outDir/remDups$case.bam $outDir/remDups$case.bai
		# check if complete
		if [ ! -f $outDir/remDups$case.bai ]; then
	                echo "Failed at remDups.bam->remDups.bai"
			echo $outDir
	                exit
	        fi
		rm $outDir/sorted$case.bam
		rm $outDir/sorted$case.bai
		# ====================================================================
		# realign BAM around indels
		# ====================================================================
		java -Xmx4g -Djava.io.tmpdir=/tmp \
	  	-jar $softDir/$GATKdir/GenomeAnalysisTK.jar \
	  	-T IndelRealigner \
		--filter_bases_not_stored \
	  	-I $outDir/remDups$case.bam \
	  	-R $refFile \
	  	-targetIntervals $intervalList \
	  	-o $outDir/realigned$case.bam \
	  	-known $indelsFile \
		-known $indelsFile2
		# check if complete
        	if [ ! -f $outDir/realigned$case.bam ]; then
	                echo "Failed at remDups.bam->realigned.bam"
			echo $outDir
	                exit
	        fi
		# ====================================================================
		# base recalibration info
		# ====================================================================
		java -jar $softDir/$GATKdir/GenomeAnalysisTK.jar \
	   	-T BaseRecalibrator \
	   	-R $refFile \
	   	-I $outDir/realigned$case.bam \
	   	-knownSites $indelsFile \
		-knownSites $snpFile \
	   	-o $outDir/recal_data$case.table #\
		#-L $captureIntervals
		# check if complete
		if [ ! -f $outDir/realigned$case.bai ]; then
	                echo "Failed at realigned.bam->realigned.bai"
			echo $outDir
	                exit
	        fi
		rm $outDir/remDups$case.bam
		rm $outDir/remDups$case.bai
		# ====================================================================
		# base quality score recalibrator
		# ====================================================================
		java -jar $softDir/$GATKdir/GenomeAnalysisTK.jar \
	   	-T PrintReads \
	   	-R $refFile \
	   	-I $outDir/realigned$case.bam \
	   	-BQSR $outDir/recal_data$case.table \
	   	-o $outDir/recalibrated$case.bam
		# check if complete
		if [ ! -f $outDir/realigned$case.bai ]; then
	                echo "Failed at realigned.bam->recalibrated.bam"
			echo $outDir
	                exit
	        fi
		# ====================================================================
		# reindex
		# ====================================================================
		samtools index $outDir/recalibrated$case.bam $outDir/recalibrated$case.bai
		# check if complete
		if [ ! -f $outDir/realigned$case.bai ]; then
	                echo "Failed at recalibrated.bam->recalibrated.bai"
			echo $outDir
	                exit
	        fi
		rm $outDir/realigned$case.bam
		rm $outDir/realigned$case.bai
		# ====================================================================
		# base quality plots
		# ====================================================================
		java -jar $softDir/$GATKdir/GenomeAnalysisTK.jar \
		-T BaseRecalibrator \
		-R $refFile \
		-I $outDir/realigned$case.bam \
		-knownSites $indelsFile \
		-knownSites $snpFile \
		-BQSR $outDir/recal_data$case.table \
		-o $outDir/recal_data2$case.table #\
		#-L $captureIntervals
		java -jar $softDir/$GATKdir/GenomeAnalysisTK.jar \
		-T AnalyzeCovariates \
		-R $refFile \
		-before $outDir/recal_data$case.table \
		-after $outDir/recal_data2$case.table \
		-plots $outDir/recal_plots$case.pdf 
	done
	# ====================================================================
	# run MuTect2
	# ====================================================================
	java -jar $softDir/$GATKdir/GenomeAnalysisTK.jar \
		-T MuTect2 \
		-R $refFile \
		-I:tumor $outDir/realignedT.bam \
		-I:normal $outDir/realignedN.bam \
		-o $outDir/output.vcf #\
		#-L $captureIntervals
done




