# function to test whether a single chromosome has exponential segment lengths
segLengthsExponential = function(seg,startCol=3,endCol=4) 
  {
  breakpoints = sapply(1:(nrow(seg)-1),
                       FUN=function(i) (seg[i+1,startCol]+seg[i,endCol])/2)
  seglengths = diff(breakpoints) # segment lengths
  # are segment lengths exponentially distributed?
  test = ks.test(seglengths, "pexp", 1/mean(seglengths)) # p>0.05 indicates that segLengths fit exponential distr
  return(1-test$p.value) # small p-value indicates exponential
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
  return(1-test$p.value) # small p value indicates multinomial
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
  sum(sims>indicesScore)/nSims
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
splitWindow = function(bedpe,seg,size=50000000,chromCol=2,startCol=3,endCol=4,chromCol1=1,posCol1=2,chromCol2=4,posCol2=5)
	{
	# p value for exponential distribution of segments
	P1 = steeleLib:::segLengthsExponential(seg,
		startCol=startCol,
		endCol=endCol)
	if(nrow(bedpe)==0) return(P1)
	# get windows of size
	chrom = paste0(unique(seg[,chromCol]))
	chromSize = range(c(bedpe[which(paste0(bedpe[,chromCol1])==chrom),posCol1],
		bedpe[which(bedpe[,chromCol2]==chrom),posCol2],
		seg[,startCol],
		seg[,endCol]))
	split = seq(from=chromSize[1],to=chromSize[2],by=size)
	if(max(split)<chromSize[2]) split = c(split,split[length(split)]+size)
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
	res = sapply(1:(length(split)-1),FUN=function(i)
		{
		checkGrange = as(paste0("chr",chrom,":",
					split[i],"-",
					split[i+1]),
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
			posCol2=posCol2)
		if(!is.na(P)) 
			{
			return(fishersMethod(c(P1,P)))
			} else {
			return(NA)
			}
		})
	if(all(is.na(res))) return(P1)
	names(res) = split[-length(split)]
	return(res)
	}

# function to run whole chromothripsis analysis
chromothripsis = function(segFile,bedpeDir,
			size=50000000, # window size
			chromCol=2, # seg chrom col
			startCol=3, # seg start col
			endCol=4, # seg end col
			bedpeChromCol1=1, # bedpe chrom1 col
			bedpePosCol1=2, # bedpe pos1 col
			bedpeChromCol2=4, # bedpe chrom2 col
			bedpePosCol2=5, # bedpe pos2col
			segSampleCol=1, # seg sample col
			doParallel=FALSE,
			nCores = NULL
			) 
	{
	if(doParallel&is.null(nCores)) nCores = detectCores()
	# read in seg file
	seg = read.table(segFile,sep="\t")
	samples = unique(seg[,segSampleCol])
	# run analysis per sample per chromosome
	# loop over samples
	Ps = sapply(samples,FUN=function(y)
		{
		print(y)
		# load bedpe for this sample
		bedpe = read.table(paste0(bedpeDir,"/",y,".brass.annot.bedpe.gz"),sep="\t")
		# data munging
		sampleIndex = paste0(seg[,segSampleCol])==paste0(y)
		subSeg = seg[which(sampleIndex),]
		chromosomes = unique(subSeg[,chromCol]) 
		# loop over chromosomes
		if(doParallel)
			{
			res = mclapply(chromosomes,FUN=function(x)
				{
				print(x)
				# keep bedpe rows for this chms
				index1=paste0(bedpe[,bedpeChromCol1])==paste0(x)
				index2=paste0(bedpe[,bedpeChromCol2])==paste0(x)
				indexBedpe = which(index1|index2)
				# keep seg rows for this chms
				indexSeg = which(paste0(subSeg[,chromCol])==paste0(x))
				# check for chromothripsis
				splitWindow(bedpe=bedpe[indexBedpe,],
					seg=subSeg[indexSeg,],
					size=size,
					chromCol=chromCol,
					startCol=startCol,
					endCol=endCol,
					chromCol1=bedpeChromCol1,
					posCol1=bedpePosCol1,
					chromCol2=bedpeChromCol2,
					posCol2=bedpePosCol2)},mc.cores=nCores)	
			} else {
			res = sapply(chromosomes,FUN=function(x)
				{
				print(x)
				# keep bedpe rows for this chms
				index1=paste0(bedpe[,bedpeChromCol1])==paste0(x)
				index2=paste0(bedpe[,bedpeChromCol2])==paste0(x)
				indexBedpe = which(index1|index2)
				# keep seg rows for this chms
				indexSeg = which(paste0(subSeg[,chromCol])==paste0(x))
				# check for chromothripsis
				splitWindow(bedpe=bedpe[indexBedpe,],
					seg=subSeg[indexSeg,],
					size=size,
					chromCol=chromCol,
					startCol=startCol,
					endCol=endCol,
					chromCol1=bedpeChromCol1,
					posCol1=bedpePosCol1,
					chromCol2=bedpeChromCol2,
					posCol2=bedpePosCol2)},simplify=FALSE)
			}
		names(res) = chromosomes
		return(res)
		},simplify=FALSE)
	names(Ps) = samples
	return(Ps)
	}

