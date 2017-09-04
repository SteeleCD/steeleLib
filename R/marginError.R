marginError=function(p,n)
	{
	1.96*sqrt((p*(1-p))/n)
	}

#p = (0:1000)/1000
#n = 100:1000
#res = sapply(p,FUN=function(x) sapply(n,FUN=function(y) marginError(x,y)))
#colnames(res) = p
#rownames(res) = n

#colours = colorRampPalette(c("purple","blue","white","orange","red"))(100)
#image(x=n,y=p,res,col=colours,main="Margin of error")
#contour(x=n,y=p,res,add=TRUE)
