# function to fit simple linear model and plot confidence intervals
plotLM = function(x,y,xlim=NULL,ylim=NULL,legloc="topleft",...)
        {
        # linear fit
        fit = lm(y~x)
        # set predict limits
        if(is.null(xlim))
                {
                newx = seq(from=min(x,na.rm=TRUE),to=max(x,na.rm=TRUE),length.out=100)
                } else {
                newx = seq(from=min(xlim,na.rm=TRUE),to=max(xlim,na.rm=TRUE),length.out=100)
                }
        # confidence intervals
        pred.conf = predict(fit, data.frame(x=newx), interval="confidence")
        # prediction intervals
        pred.pred = predict(fit, data.frame(x=newx), interval="prediction") 
        # plot data
        plot(x,y,pch=4,xlim=xlim,ylim=ylim,...)
        # plot fit
        abline(fit,lwd=2)
        # plot x=y
        abline(a=0,b=1,lty=2)
        # plot intervals
        polygon(x=c(newx,rev(newx)),
                y=c(pred.conf[,"lwr"],rev(pred.conf[,"upr"])),
                col=rgb(0.5,0.5,0.5,0.2),
                border=NA)
        polygon(x=c(newx,rev(newx)),
                y=c(pred.pred[,"lwr"],rev(pred.pred[,"upr"])),
                col=rgb(0.5,0.5,0.5,0.2),
                border=NA)
        # plot legend
        legend(legloc,
                legend=c("Observations","Linear fit","x=y","95% CI","95% PI"),
                lty=c(NA,1,2,NA,NA),
                lwd=c(NA,2,1,NA,NA),
                pch=c(4,NA,NA,15,15),
                col=c(1,1,1,rgb(0.5,0.5,0.5,0.4),rgb(0.5,0.5,0.5,0.2)))
        }
