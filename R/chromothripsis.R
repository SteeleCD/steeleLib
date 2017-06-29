# function to test whether a single chromosome has exponential segment lengths
segLengthsExponential = function(seg,startCol=3,endCol=4) 
  {
  breakpoints = sapply(1:(nrow(seg)-1),
                       FUN=function(i) (seg[i+1,startCol]+seg[i,endCol])/2)
  seglengths = diff(breakpoints) # segment lengths
  test = ks.test(seglengths, "pexp", 1/mean(seglengths)) # are segment lengths exponentially distributed?
  return(test$p.value)
  }

# randomness of DNA fragment joins
# counts of +/+ -/- +/- -/+ should be random (1/4,1/4,1/4,1/4)
randomJoins = function(bedpe,direction1col=9,direction2col=10)
  {
  joins = paste0(bedpe[,direction1col],bedpe[,direction2col])
  counts = table(joins)
  if(length(counts)<4) counts = c(counts,rep(0,4-length(counts)))
  dmultinom(counts,prob=rep(0.25,4))
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

