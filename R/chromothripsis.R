# function to test whether a single chromosome has exponential segment lengths
segLengthsExponential = function(seg,startCol=3,endCol=4) 
  {
  breakpoints = sapply(1:(nrow(seg)-1),
                       FUN=function(i) (seg[i+1,startCol]+seg[i,endCol])/2)
  seglengths = diff(breakpoints) # segment lengths
  # are segment lengths exponentially distributed?
  test = ks.test(seglengths, "pexp", 1/mean(seglengths)) # p>0.05 indicates that segLengths fit exponential distr
  return(test$p.value) # p>0.05 indicates exponential  (suggesting chromothripsis)
  }

# randomness of DNA fragment joins
# counts of +/+ -/- +/- -/+ should be random (1/4,1/4,1/4,1/4)
randomJoins = function(bedpe,direction1col=9,direction2col=10)
  {
  joins = paste0(bedpe[,direction1col],bedpe[,direction2col])
  counts = table(joins)
  if(length(counts)<4) counts = c(counts,rep(0,4-length(counts)))
  # goodness of fit test to multinomial
  test = chisq.test(counts,p=rep(0.25,4)) # p>0.05 indicates that counts fit multinomial distr
  return(test$p.value) # p>0.05 indicates multinomial  (suggesting chromothripsis)
  }

# randomness of DNA fragment order
# two sides of each breakpoint should be random draws from all breakpoint positions
randomOrder = function(bedpe,chromCol1=1,posCol1=2,chromCol2=4,posCol2=5,nSims=1000)
  {
  # breakpoints
  breakpoints1 = bedpe[,c(chromCol1,posCol1)]
  breakpoints2 = bedpe[,c(chromCol2,posCol2)]
  colnames(breakpoints1)=colnames(breakpoints2)=c("chrom","pos")
  breakpoints=rbind(breakpoints1,breakpoints2)
  # order breakpoints
  breakpoints = breakpoints[order(breakpoints[,1],breakpoints[,2]),]
  breakpoints = unique(breakpoints)
  # indices
  indices = apply(bedpe,MARGIN=1,FUN=function(x)
    {
    x = gsub(" ","",x)
      c(
      which(breakpoints[,1]==unlist(x[chromCol1])&
              as.numeric(breakpoints[,2])==as.numeric(x[posCol1])),
      which(breakpoints[,1]==unlist(x[chromCol2])&
              as.numeric(breakpoints[,2])==as.numeric(x[posCol2]))
      )})
  indicesScore = mean(abs(indices[2,]-indices[1,]))
  # monte carlo simulations
  sims = replicate(nSims,abs(diff(sample(nrow(breakpoints),2))))
  # p value
  sum(sims>indicesScore)/nSims # p>0.05 indicates random draw (suggesting chromothripsis)
  }

# ability to walk chromosome
walkChrom = function()
  {

  }

# combining p-values with fishers method
fishersMethod = function(Ps)
  {
  test=-2*sum(log(Ps))
  pchisq(test,df=2*length(Ps),lower.tail=FALSE)
  }

# run a single test
runSingle = function(bedpe,direction1col=9,direction2col=10,chromCol1=1,posCol1=2,chromCol2=4,posCol2=5,nSims=1000,
		seg,startCol=3,endCol=4)
	{
	dobedpe = nrow(bedpe)>0
	#doseg = nrow(seg)>2
	#if(!doseg) return(NA)
	#P1 = segLengthsExponential(seg,
	#	startCol=startCol,
	#	endCol=endCol)
	if(dobedpe)
		{
		P2 = randomJoins(bedpe,
			direction1col=direction1col,
			direction2col=direction2col)
		P3 = randomOrder(bedpe,
			chromCol1=chromCol1,
			posCol1=posCol1,
			chromCol2=chromCol2,
			posCol2=posCol2,
			nSims=nSims)
		return(c(P2,P3))
		} else {
		return(NA)
		}
	}

# split into windows, then check for chromothripsis
splitWindow = function(bedpe,seg,size=1e7,gap=1e6,chromCol=2,startCol=3,endCol=4,chromCol1=1,posCol1=2,chromCol2=4,posCol2=5,direction1col=9,direction2col=10)
	{
	if(nrow(seg)<3) return(NA)
	# p value for exponential distribution of segments
	P1 = segLengthsExponential(seg,
		startCol=startCol,
		endCol=endCol)
	if(nrow(bedpe)==0) return(P1)
	# get windows of size
	chrom = paste0(unique(seg[,chromCol]))
	chromSize = range(c(bedpe[which(paste0(bedpe[,chromCol1])==chrom),posCol1],
		bedpe[which(bedpe[,chromCol2]==chrom),posCol2],
		seg[,startCol],
		seg[,endCol]))
	split = seq(from=min(chromSize),to=max(chromSize)-size,by=gap)
	if((max(split)+size)<max(chromSize)) 
		{
		split = c(split,split[length(split)]+gap)
		names(split)[length(split)] = split[length(split)]
		}
	# get Granges
	segGrange = as(paste0("chr",
			seg[,chromCol],":",
			seg[,startCol],"-",
			seg[,endCol]),
		"GRanges")
	bedpeGrange1 = as(paste0("chr",
			bedpe[,chromCol1],":",
			bedpe[,posCol1],"-",
			bedpe[,posCol1]),
		"GRanges")
	bedpeGrange2 = as(paste0("chr",
			bedpe[,chromCol2],":",
			bedpe[,posCol2],"-",
			bedpe[,posCol2]),
		"GRanges")
	# check for chromothripsis in each window
	res = sapply(1:length(split),FUN=function(i)
		{
		checkGrange = as(paste0("chr",chrom,":",
					split[i],"-",
					split[i]+size),
				"GRanges")
		segIndex = findOverlaps(segGrange,checkGrange)
		segIndex = segIndex@from
		bedpeIndex1 = findOverlaps(bedpeGrange1,checkGrange)
		bedpeIndex2 = findOverlaps(bedpeGrange2,checkGrange)
		bedpeIndex = unique(bedpeIndex1@from,bedpeIndex2@from)
		P = runSingle(bedpe=bedpe[bedpeIndex,],
			seg=seg[segIndex,],
			startCol=startCol,
			endCol=endCol,
			chromCol1=chromCol1,
			posCol1=posCol1,
			chromCol2=chromCol2,
			posCol2=posCol2,
			direction1col=direction1col,
			direction2col=direction2col)
		if(!any(is.na(P))) 
			{
			return(fishersMethod(c(P1,P)))
			} else {
			return(NA)
			}
		})
	if(all(is.na(res))) return(P1)
	res[which(is.na(res))] = 0
	names(res) = split
	return(res)
	}

# function to get runs and clean output
getRuns = function(chromScores,chrom,samp,size)
	{
	if(length(chromScores)>1)
		{
		chromBool = chromScores>0.05
		if(!any(chromBool)) return(NULL)
		runs = rle(chromBool)
		ends = cumsum(runs$lengths)
		starts = c(1,ends[-length(ends)]+1)
		windowStarts = as.numeric(names(chromBool)[starts[which(runs$values==TRUE)]])
		windowEnds = as.numeric(names(chromBool)[ends[which(runs$values==TRUE)]])
		windowEnds = windowEnds+size
		return(cbind(samp,chrom,windowStarts,windowEnds))
		} else {
		#return(cbind(samp,chrom))
		# only return results if have looked at translocations
		return(NULL)
		}
	}

# read in a file
readFile = function(file,head)
	{
	ending = rev(strsplit(file,split="[.]")[[1]])[1]
	if(ending=="csv") return(read.csv(file,head=head))
	if(ending%in%c("txt","tsv")) return(read.table(file,sep="\t",head=head))
	return(read.table(file,sep="\t",head=head))	
	}

# function to run whole chromothripsis analysis
chromothripsis = function(segFile, # combined seg file
			bedpeFile=NULL, # directory of separate bedpes, or single bedpe file
			size=1e7, # window size
			gap=1e6, # gap between sliding windows
			chromCol=2, # seg chrom col
			startCol=3, # seg start col
			endCol=4, # seg end col
			bedpeChromCol1=1, # bedpe chrom1 col
			bedpePosCol1=2, # bedpe pos1 col
			direction1col=9, # orientation of first partner
			bedpeChromCol2=4, # bedpe chrom2 col
			bedpePosCol2=5, # bedpe pos2col
			direction2col=10, # orientation of second partner
			segSampleCol=1, # seg sample col
			bedpeSampleCol=1, # seg sample col
			doParallel=FALSE, # parallel computation
			nCores = NULL,	# number of cores
			samplesToRun = NULL, # which samples to run
			chromsToRun = NULL, # which chromosomes to run
			sepbedpe=TRUE, # Are bedpes separate
			bedpeHead=FALSE, # does bedpe file have header
			segHead=TRUE, # does seg file have header
			bedpeEnding=".brass.annot.bedpe.gz"
			) 
	{
	if(doParallel&is.null(nCores)) nCores = detectCores()
	# read in seg file
	seg = readFile(segFile,segHead)
	if(!is.null(samplesToRun)) seg = seg[which(seg[,segSampleCol]%in%samplesToRun),]
	samples = unique(seg[,segSampleCol])
	# read in bedpe
	if(!sepbedpe) allbedpe = readFile(bedpeFile,bedpeHead)
	# run analysis per sample per chromosome
	# loop over samples
	Ps = sapply(samples,FUN=function(y)
		{
		print(y)
		# load bedpe for this sample
		if(sepbedpe)
			{
			bedpe = readFile(paste0(bedpeFile,"/",y,bedpeEnding),bedpeHead)
			} else {
			bedpe = allbedpe[which(paste0(allbedpe[,bedpeSampleCol])==paste0(y)),]
			}
		# data munging
		sampleIndex = paste0(seg[,segSampleCol])==paste0(y)
		subSeg = seg[which(sampleIndex),]
		if(is.null(chromsToRun))
			{
			chromosomes = unique(subSeg[,chromCol])
			} else {
			chromosomes = chromsToRun
			} 
		# chromothripsis calculation
		if(doParallel)
			{
			# loop over chromosomes
			res = mclapply(chromosomes,FUN=function(x)
				{
				print(x)
				# keep bedpe rows for this chms
				index1=paste0(bedpe[,bedpeChromCol1])==paste0(x)
				index2=paste0(bedpe[,bedpeChromCol2])==paste0(x)
				indexBedpe = which(index1|index2)
				if(length(indexBedpe)==0) return(NULL)
				# keep seg rows for this chms
				indexSeg = which(paste0(subSeg[,chromCol])==paste0(x))
				# check for chromothripsis
				chromScores = splitWindow(bedpe=bedpe[indexBedpe,],
					seg=subSeg[indexSeg,],
					size=size,
					chromCol=chromCol,
					startCol=startCol,
					endCol=endCol,
					chromCol1=bedpeChromCol1,
					posCol1=bedpePosCol1,
					chromCol2=bedpeChromCol2,
					posCol2=bedpePosCol2,
					direction1col=direction1col,
					direction2col=direction2col)
				# just output regions that are chromothriptic
				getRuns(chromScores,paste0(x),paste0(y),size)},mc.cores=nCores)	
			} else {
			# loop over chromosomes
			res = sapply(chromosomes,FUN=function(x)
				{
				print(x)
				# keep bedpe rows for this chms
				index1=paste0(bedpe[,bedpeChromCol1])==paste0(x)
				index2=paste0(bedpe[,bedpeChromCol2])==paste0(x)
				indexBedpe = which(index1|index2)
				if(length(indexBedpe)==0) return(NULL)
				# keep seg rows for this chms
				indexSeg = which(paste0(subSeg[,chromCol])==paste0(x))
				# check for chromothripsis
				chromScores = splitWindow(bedpe=bedpe[indexBedpe,],
					seg=subSeg[indexSeg,],
					size=size,
					chromCol=chromCol,
					startCol=startCol,
					endCol=endCol,
					chromCol1=bedpeChromCol1,
					posCol1=bedpePosCol1,
					chromCol2=bedpeChromCol2,
					posCol2=bedpePosCol2,
					direction1col=direction1col,
					direction2col=direction2col)
				# just output regions that are chromothriptic
				getRuns(chromScores,paste0(x),paste0(y),size)},simplify=FALSE)
			}
		names(res) = chromosomes
		return(res)
		},simplify=FALSE)
	names(Ps) = samples
	return(Ps)
	}

