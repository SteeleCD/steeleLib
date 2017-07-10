# function to test whether a single chromosome has exponential segment lengths
segLengthsExponential = function(seg,startCol=3,endCol=4) 
  {
  breakpoints = sapply(1:(nrow(seg)-1),
                       FUN=function(i) (seg[i+1,startCol]+seg[i,endCol])/2)
  seglengths = diff(breakpoints) # segment lengths
  test = ks.test(seglengths, "pexp", 1/mean(seglengths)) # are segment lengths exponentially distributed?
  return(test$p.value) # small p-value indicates exponential
  }

# randomness of DNA fragment joins
# counts of +/+ -/- +/- -/+ should be random (1/4,1/4,1/4,1/4)
randomJoins = function(bedpe,direction1col=9,direction2col=10)
  {
  joins = paste0(bedpe[,direction1col],bedpe[,direction2col])
  counts = table(joins)
  if(length(counts)<4) counts = c(counts,rep(0,4-length(counts)))
  1-dmultinom(counts,prob=rep(0.25,4)) # small densities of dmultinom = not random
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
	P1 = segLengthsExponential(seg,
		startCol=startCol,
		endCol=endCol)
	# get windows of size
	chrom = unique(seg[,chromCol])
	chromSize = range(c(bedpe[which(bedpe[,chromCol1]==chrom),posCol1],
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
	names(res) = split[-length(split)]
	return(res)
	}

# function to run whole chromothripsis analysis
chromothripsis = function(segFile,bedpeDir,
			size=50000000, # window size
			chromCol=2, # seg chrom col
			startCol=3, # seg start col
			endCol=4, # seg end col
			chromCol1=1, # bedpe chrom1 col
			posCol1=2, # bedpe pos1 col
			chromCol2=4, # bedpe chrom2 col
			posCol2=5, # bedpe pos2col
			segSampleCol=1 # seg sample col
			) 
	{
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
		# loop over chromosomes
		res = sapply(unique(bedpe[,c(chromCol1,chromCol2)]),FUN=function(x)
			{
			print(x)
			splitWindow(bedpe=bedpe[which(paste0(bedpe[,bedpeChromCol1])==paste0(x)|
					paste0(bedpe[,bedpeChromCol2])==paste0(x)),],
				seg=seg[which(paste0(seg[,segSampleCol])==paste0(y)&
					paste0(seg[,chromCol])==paste0(x)),],
				size=size,
				chromCol=chromCol,
				startCol=startCol,
				endCol=endCol,
				chromCol1=chromCol1,
				posCol1=posCol1,
				chromCol2=chromCol2,
				posCol2=posCol2)})
		names(res) = samples
		return(res)
		},simplify=FALSE)
	names(Ps) = samples
	return(Ps)
	}

