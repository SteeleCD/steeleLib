# get bump limits using a sliding window and change of gradient
getBumpLimsWindow = function(betas,windowSize=50,negCount=c(-40),gradLimit=10,retLim=TRUE)
{
  # GET DENSITY
  dens = density(betas)
  # get gradients of density
  index1 = 2:length(dens$y)
  gradient = (dens$y[index1]-dens$y[index1-1])/(dens$x[index1]-dens$x[index1-1])
  # get change in gradients of density
  index2 = 2:length(gradient)
  d2 = (gradient[index2]-gradient[index2-1])/(dens$x[index2+1]-dens$x[index2])
  # new method - sliding window of change in gradient
  negIndex = sapply(1:(length(d2)-windowSize),FUN=function(x) sum(sign(d2[x:(x+windowSize)])))
  # check that enough of window is negative
  test = which(negIndex<negCount)
  if(length(test)>0)
  {
    # window limits
    tmp = range(test)
    tmp[2] = tmp[2]+windowSize
    # remove those that have too large a gradient
    if(any(abs(gradient[tmp[1]:tmp[2]])>gradLimit)) tmp=NA
  } else {
    # no bump
    tmp = NA
  }
  # plot
  par(mfrow=c(3,1),mar=c(3,4,0,0))
  # plot density
  plot(dens,main=NA)
  abline(v=dens$x[tmp],lty=2)
  # plot gradient
  plot(dens$x[index1],gradient,col=as.factor(sign(gradient)))
  abline(h=0,lty=2,col="gray")
  abline(v=dens$x[tmp],lty=2)
  # plot change in gradient
  plot(dens$x[index2+1],d2,col=as.factor(sign(d2)))
  abline(h=0,lty=2,col="gray")
  abline(v=dens$x[tmp],lty=2)
  # return bum limits
  if(retLim)
  {
  if(!is.na(tmp))
  {
    out = list(lims=dens$x[tmp],relHeight=max(dens$y[tmp[1]:tmp[2]])/max(dens$y))
  } else {out=NA}
  } else {
    if(!is.na(tmp))
    {
      out = TRUE
    } else {out=FALSE}
  }
  return(out)
}


# get bump limits using change in sign of gradient
getBumpLims = function(betas)
{
  # GET DENSITY
  dens = density(betas)
  index = 2:length(dens$y)
  # get gradients of density
  gradient = (dens$y[index]-dens$y[index-1])/(dens$x[index]-dens$x[index-1])
  # signs of gradient
  signs = sign(gradient)
  # minimums of density
  mins = which(signs[2:length(signs)]-signs[1:(length(signs)-1)]==2)+1
  plot(dens)
  abline(v=dens$x[mins],lty=2)
  # check if have bump
  check = length(mins)>1
  if(check) 
  {
    # get maximums
    maxs = which(signs[2:length(signs)]-signs[1:(length(signs)-1)]==c(-2))+1
    maxValue = dens$y[maxs[2]]
    minValues = dens$x[mins]
    minValue = min(minValues)
    bumpHeight = maxValue-minValue
    relHeight = bumpHeight/max(dens$y)
    return(list(lims=dens$x[mins],height=bumpHeight,relHeight=relHeight))
  } else {
    return(NA)
  }
}
