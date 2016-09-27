# slightly altered ChAMP
library(limma)
library(DNAcopy)
library(preprocessCore)
library(impute)
library(marray)
library(sva)
library(IlluminaHumanMethylation450kmanifest)
library(minfi)
library(ChAMPdata)
library(Illumina450ProbeVariants.db)
library(wateRmelon)
library(plyr)
library(GenomicRanges)
library(RefFreeEWAS)
library(isva)
library(qvalue)
library(bumphunter)
library(doParallel)
library(quadprog)

# ==========================================================================
# CHANGEPOINTS
# ==========================================================================
changepoints <- function(genomdat, data.type="logratio", alpha=0.01, weights=
                         NULL, sbdry, sbn, nperm=10000, p.method="hybrid", 
                         min.width=2, kmax=25, nmin=200, trimmed.SD=NULL, 
                         undo.splits="none", undo.prune=0.05, undo.SD=3,
                         verbose=1, ngrid=100, tol=1e-6)
  {
    n <- length(genomdat)
    if (missing(trimmed.SD)) trimmed.SD <- mad(diff(genomdat))/sqrt(2)
#   start with the whole 
    seg.end <- c(0,n)
    k <- length(seg.end)
    change.loc <- NULL
    weighted <- ifelse(is.null(weights), FALSE, TRUE) 
    while (k > 1)
      {
        current.n <- seg.end[k]-seg.end[k-1]
        if (verbose>=3) cat(".... current segment:",seg.end[k-1]+1,"-",seg.end[k],"\n")
        if(current.n >= 2*min.width) {
          current.genomdat <- genomdat[(seg.end[k-1]+1):seg.end[k]]
#   check whether hybrid method needs to be used
          hybrid <- FALSE
          delta <- 0
          if ((p.method=="hybrid") & (nmin < current.n)) {
            hybrid <- TRUE
            delta <- (kmax+1)/current.n
          }
#   call the changepoint routine
          if (weighted) {
#   get the weights for the current set of probes
            current.wts <- weights[(seg.end[k-1]+1):seg.end[k]]
            current.rwts <- sqrt(current.wts)
            current.cwts <- cumsum(current.wts)/sqrt(sum(current.wts))
#   if all values of current.genomdat are the same don't segment
            if (isTRUE(all.equal(diff(range(current.genomdat)), 0))) {
              zzz <- list()
              zzz$ncpt <- 0
            } else {
#   centering the current data will save a lot of computations later
              current.avg <- sum(current.genomdat*current.wts)/sum(current.wts)
              current.genomdat <- current.genomdat - current.avg
#   need total sum of squares too
              current.tss <- sum(current.wts*(current.genomdat^2))
              zzz <- .Fortran("wfindcpt",
                              n=as.integer(current.n),
                              x=as.double(current.genomdat),
                              tss=as.double(current.tss),
                              wts=as.double(current.wts),
                              rwts=as.double(current.rwts),
                              cwts=as.double(current.cwts),
                              px=double(current.n),
                              sx=double(current.n),
                              nperm=as.integer(nperm),
                              cpval=as.double(alpha),
                              ncpt=integer(1),
                              icpt=integer(2),
                              hybrid=as.logical(hybrid),
                              al0=as.integer(min.width),
                              hk=as.integer(kmax),
                              mncwt=double(kmax),
                              delta=as.double(delta),
                              ngrid=as.integer(ngrid),
                              sbn=as.integer(sbn),
                              sbdry=as.integer(sbdry),
                              tol= as.double(tol),
                              PACKAGE="DNAcopy")
            }
          } else { 
#   if all values of current.genomdat are the same don't segment
            if (isTRUE(all.equal(diff(range(current.genomdat)), 0))) {
              zzz <- list()
              zzz$ncpt <- 0
            } else {
#   centering the current data will save a lot of computations later
              current.avg <- mean(current.genomdat)
              current.genomdat <- current.genomdat - current.avg
#   need total sum of squares too
              current.tss <- sum(current.genomdat^2)
              zzz <- .Fortran("fndcpt",
                              n=as.integer(current.n),
                              x=as.double(current.genomdat),
                              tss=as.double(current.tss),
                              px=double(current.n),
                              sx=double(current.n),
                              nperm=as.integer(nperm),
                              cpval=as.double(alpha),
                              ncpt=integer(1),
                              icpt=integer(2),
                              ibin=as.logical(data.type=="binary"),
                              hybrid=as.logical(hybrid),
                              al0=as.integer(min.width),
                              hk=as.integer(kmax),
                              delta=as.double(delta),
                              ngrid=as.integer(ngrid),
                              sbn=as.integer(sbn),
                              sbdry=as.integer(sbdry),
                              tol= as.double(tol),
                              PACKAGE="DNAcopy")
            }
          }
        } else {
          zzz <- list()
          zzz$ncpt <- 0
        }
        if(zzz$ncpt==0) change.loc <- c(change.loc,seg.end[k])
        seg.end <- switch(1+zzz$ncpt,seg.end[-k],
                          c(seg.end[1:(k-1)],seg.end[k-1]+zzz$icpt[1],seg.end[k]),
                          c(seg.end[1:(k-1)],seg.end[k-1]+zzz$icpt,seg.end[k]))
        k <- length(seg.end)
        if(verbose>=3) cat(".... segments to go:",seg.end,"\n")
      }
    seg.ends <- rev(change.loc)
    nseg <- length(seg.ends)
    lseg <- diff(c(0,seg.ends))
    if (nseg > 1) {
        if (undo.splits == "prune") {
            lseg <- DNAcopy:::changepoints.prune(genomdat, lseg, undo.prune)
        }
        if (undo.splits == "sdundo") {
            lseg <- DNAcopy:::changepoints.sdundo(genomdat, lseg, trimmed.SD, undo.SD)
        }
    }
    segmeans <- 0*lseg
	# MINE
    segsds <- 0*lseg
    ll <- uu <- 0
    for(i in 1:length(lseg)) {
      uu <- uu + lseg[i]
      if (weighted) {
        segmeans[i] <- sum(genomdat[(ll+1):uu]*weights[(ll+1):uu])/sum(weights[(ll+1):uu])
	# MINE
	nzero = length(which(weights[(ll+1):uu]!=0))
	segsds[i] <- sum(weights[(ll+1):uu]*(genomdat[(ll+1):uu]-segmeans[i])^2)/(((nzero-1)/nzero)*sum(weights[(ll+1):uu]))
      } else {
        segmeans[i] <- mean(genomdat[(ll+1):uu])
	# MINE
	segsds[i] <- sd(genomdat[(ll+1):uu])
      }
      ll <- uu
    }
    list("lseg" = lseg, "segmeans" = segmeans, "segsds" = segsds)
  }


# ==========================================================================
# SEGMENT
# ==========================================================================

segment <- function(x, weights=NULL, alpha=0.01, nperm=10000, p.method= 
                    c("hybrid","perm"), min.width=2, kmax=25, nmin=200, 
                    eta=0.05, sbdry=NULL, trim = 0.025, undo.splits=
                    c("none","prune", "sdundo"), undo.prune=0.05, undo.SD=3,
                    verbose=1)
  {
    if (!inherits(x, 'CNA')) stop("First arg must be a copy number array object")
    call <- match.call()
    if (min.width < 2 | min.width > 5) stop("minimum segment width should be between 2 and 5")
    if (nmin < 4*kmax) stop("nmin should be >= 4*kmax")
    if (missing(sbdry)) {
      if (nperm==10000 & alpha==0.01 & eta==0.05) {
        if (!exists("default.DNAcopy.bdry")) data(default.DNAcopy.bdry, package="DNAcopy",envir=environment())
        sbdry <- get("default.DNAcopy.bdry", envir=environment())
      } else {
        max.ones <- floor(nperm*alpha) + 1
        sbdry <- getbdry(eta, nperm, max.ones)
      }
    }
    weighted <- ifelse(missing(weights), FALSE, TRUE)
#   rudimentary error checking for weights
    if (weighted) {
      if (length(weights) != nrow(x)) stop("length of weights should be the same as the number of probes")
      if (min(weights) <= 0) stop("all weights should be positive")
    }
    sbn <- length(sbdry)
    nsample <- ncol(x)-2
    sampleid <- colnames(x)[-(1:2)]
    uchrom <- unique(x$chrom)
    data.type <- attr(x, "data.type")
    p.method <- match.arg(p.method)
    undo.splits <- match.arg(undo.splits)
    segres <- list()
    segres$data <- x
    allsegs <- list()
    allsegs$ID <- NULL
    allsegs$chrom <- NULL
    allsegs$loc.start <- NULL
    allsegs$loc.end <- NULL
    allsegs$num.mark <- NULL
    allsegs$seg.mean <- NULL
    # MINE
    allsegs$seg.sd <- NULL
    segRows <- list()
    segRows$startRow <- NULL
    segRows$endRow <- NULL
    for (isamp in 1:nsample) {
      if (verbose>=1) cat(paste("Analyzing:", sampleid[isamp],"\n"))
      genomdati <- x[,isamp+2]
      ina <- which(is.finite(genomdati))
      genomdati <- genomdati[ina]
      trimmed.SD <- sqrt(DNAcopy:::trimmed.variance(genomdati, trim))
      chromi <- x$chrom[ina]
#      maploci <- x$maploc[ina]
      if (weighted) {
        wghts <- weights[ina]
      } else {
        wghts <- NULL
      }
      sample.lsegs <- NULL
      sample.segmeans <- NULL
      # MINE
      sample.segsds <- NULL
      for (ic in uchrom) {
        if (verbose>=2) cat(paste("  current chromosome:", ic, "\n"))
        segci <- changepoints(genomdati[chromi==ic], data.type, alpha, wghts,
                              sbdry, sbn, nperm, p.method, min.width, kmax,
                              nmin, trimmed.SD, undo.splits, undo.prune,
                              undo.SD, verbose)
        sample.lsegs <- c(sample.lsegs, segci$lseg)
        sample.segmeans <- c(sample.segmeans, segci$segmeans)
	# MINE
	sample.segsds <- c(sample.segsds, segci$segsds)
      }
      sample.nseg <- length(sample.lsegs)
      sample.segs.start <- ina[cumsum(c(1,sample.lsegs[-sample.nseg]))]
      sample.segs.end <- ina[cumsum(sample.lsegs)]
      allsegs$ID <- c(allsegs$ID, rep(isamp,sample.nseg))
      allsegs$chrom <- c(allsegs$chrom, x$chrom[sample.segs.end])
      allsegs$loc.start <- c(allsegs$loc.start, x$maploc[sample.segs.start])
      allsegs$loc.end <- c(allsegs$loc.end, x$maploc[sample.segs.end])
      allsegs$num.mark <- c(allsegs$num.mark, sample.lsegs)
      allsegs$seg.mean <- c(allsegs$seg.mean, sample.segmeans)
      # MINE
      allsegs$seg.sd <- c(allsegs$seg.sd, sample.segsds)
      segRows$startRow <- c(segRows$startRow, sample.segs.start)
      segRows$endRow <- c(segRows$endRow, sample.segs.end)
    }
    allsegs$ID <- sampleid[allsegs$ID]
    allsegs$seg.mean <- round(allsegs$seg.mean, 4)
    allsegs <- as.data.frame(allsegs)
    allsegs$ID <- as.character(allsegs$ID)
    segres$output <- allsegs
    segres$segRows <- as.data.frame(segRows)
    segres$call <- call
    if (weighted) segres$weights <- weights
    class(segres) <- "DNAcopy"
    segres
  }


mychamp.CNA <-
function(intensity=myLoad$intensity, pd=myLoad$pd, loadFile=FALSE, batchCorrect=TRUE, file="intensity.txt", resultsDir=paste(getwd(),"resultsChamp",sep="/"), sampleCNA=TRUE,plotSample=TRUE, filterXY=TRUE, groupFreqPlots=TRUE,freqThreshold=0.3, control=TRUE, controlGroup="Control",arraytype="450K",seg.min.width=2,seg.nperm=10000)
{
    if(arraytype=="EPIC"){
        data(probe.features.epic)
    }else{
        data(probe.features)
    }
	#normalize.quantiles<-NULL
    #rm(normalize.quantiles)
	#control.intsqnlog<-NULL
	#CNA<-NA
	#rm(CNA)
	#smooth.CNA<-NA
	#rm(smooth.CNA)
# MINE
	#segment<-NA
	#rm(segment)
	
	message("Run champ.CNA")
	newDir=paste(resultsDir,"CNA",sep="/")
    
	if(!file.exists(resultsDir)){dir.create(resultsDir)}
	if(!file.exists(newDir))
	{
		dir.create(newDir)
	}
		
	if(loadFile)
	{
		ints = read.table(file,row.names =T, sep = "\t")
        if(filterXY)
        {
            autosomes=probe.features[!probe.features$CHR %in% c("X","Y"), ]
            ints=ints[row.names(ints) %in% row.names(autosomes), ]
        }
	}else
	{
		ints=intensity
	}
	
	if(control)
	{
		        	
        if(controlGroup != "champCtl" & !(controlGroup %in% pd$Sample_Group))
        {
        	message("You have chosen ", controlGroup, " as the reference and this does not exist in your sample sheet (column Sample_Group). The analysis will run with ChAMP blood controls.")
        	controlGroup="champCtl"
        }
        
        if(controlGroup == "champCtl")
        {
            #load("champBloodCtls.Rdata")
            data(champBloodCtls)
            ctlIntensity=bloodCtl$intensity
            ctlIntensity=ctlIntensity[which(row.names(ctlIntensity) %in% row.names(ints)),]

	#MINE
	bloodCtl$pd[,3] = "champCtl"
		
            ints=cbind(ints,ctlIntensity)
            pd=rbind(pd,bloodCtl$pd)
            batchCorrect=F
        }
	}
	
	#Extracts names of samples 
	names<-colnames(ints)

	#Quantile normalises intensities	
	intsqn<-preprocessCore::normalize.quantiles(as.matrix(ints))
	colnames(intsqn)<-names


	#Calculates Log2
	intsqnlog<-log2(intsqn)
    	
	if(batchCorrect)
	{
		message("Run batch correction for CNA")
		combat=champ.runCombat(beta.c=intsqnlog,pd=pd,logitTrans=FALSE)
		intsqnlog=combat$beta
	}

	if(control)
	{
        	#separates case from control(reference sample/samples)        	
        	message("champ.CNA is using the samples you have defined as ", controlGroup," as the reference for calculating copy number aberrations.")
        	if(controlGroup=="champCtl")
        	{
        		message("As you are using the ChAMP controls Combat cannot adjust for batch effects. Batch effects may affect your dataset.")
        	}
        
        	
        	controlSamples= pd[which(pd$Sample_Group==controlGroup),]
        	caseSamples= pd[which(pd$Sample_Group != controlGroup),]
		case.intsqnlog<-intsqnlog[,which(colnames(intsqnlog) %in% caseSamples$Sample_Name)]
	        control.intsqnlog<-intsqnlog[,which(colnames(intsqnlog) %in% controlSamples$Sample_Name)]
       	 	control.intsqnlog<-rowMeans(control.intsqnlog)
        
        	intsqnlogratio<-case.intsqnlog
        	for(i in 1:ncol(case.intsqnlog))
        	{
            		intsqnlogratio[,i]<-case.intsqnlog[,i]-control.intsqnlog
        	}
			
	}else
	{
			message("champ.CNA is using an average of all your samples as the reference for calculating copy number aberrations.")
			#Creates alternate reference sample from rowMeans if proper reference /control is not available 
			case.intsqnlog<-intsqnlog[,1:length(names)]
			ref.intsqnlog<-rowMeans(intsqnlog)
   
        	intsqnlogratio<-intsqnlog
        	for(i in 1:ncol(case.intsqnlog))
        	{
            		intsqnlogratio[,i]<-case.intsqnlog[,i]-ref.intsqnlog
        	}
	}

	ints <- data.frame(ints, probe.features$MAPINFO[match(rownames(ints), rownames(probe.features))])
	names(ints)[length(ints)] <- "MAPINFO"
	ints <- data.frame(ints, probe.features$CHR[match(rownames(ints), rownames(probe.features))])
	names(ints)[length(ints)] <- "CHR"

	#Replaces Chr X and Y with 23 and 24
	levels(ints$CHR)[levels(ints$CHR)=='X']='23'
	levels(ints$CHR)[levels(ints$CHR)=='Y']='24'

	#converts CHR factors to numeric 
	CHR<-as.numeric(levels(ints$CHR))[ints$CHR]

	#need to copy in MAPINFO
	ints$MAPINFO<-as.numeric(ints$MAPINFO)
	MAPINFO=probe.features$MAPINFO[match(rownames(ints), rownames(probe.features))]
	MAPINFO<-as.numeric(MAPINFO)

	#Runs CNA and generates individual DNA Copy Number profiles
	if(sampleCNA)
	{
		message("Saving Copy Number information for each Sample")
		for(i in 1:ncol(case.intsqnlog))
		{
			CNA.object <- DNAcopy::CNA(cbind(intsqnlogratio[,i]), CHR, MAPINFO ,data.type = "logratio", sampleid = paste(colnames(case.intsqnlog)[i],"qn"))
			smoothed.CNA.object <- smooth.CNA(CNA.object)
			segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose = 1,alpha=0.001, undo.splits="sdundo", undo.SD=2,min.width=seg.min.width,nperm=seg.nperm)
			if(plotSample)
			{
				imageName<-paste(colnames(case.intsqnlog)[i],"qn.jpg",sep="")
				imageName=paste(newDir,imageName,sep="/")
				jpeg(imageName)
				plot(segment.smoothed.CNA.object, plot.type = "w", ylim=c(-6,6))
				dev.off()
			}
			seg<-segment.smoothed.CNA.object$output
			table_name<-paste(newDir,"/",colnames(case.intsqnlog)[i],"qn.txt",sep="")
			write.table(seg,table_name, sep="\t", col.names=T, row.names=F, quote=FALSE)
		}
	}
	
	##group Frequency plots
	if(groupFreqPlots)
	{
	message("Saving frequency plots for each group")
        
        if(control)
        {
            groups = unique(pd$Sample_Group)
	    groups = groups[!groups==controlGroup]
        }else
        {
            groups = unique(pd$Sample_Group)
        }
	
		for(g in 1:length(groups))
		{
		
			pd_group = pd[which(pd$Sample_Group==groups[g]),]
			data_group=intsqnlogratio[,which(colnames(intsqnlogratio) %in% pd_group$Sample_Name)]
			ints_group=ints[,which(colnames(ints) %in% pd_group$Sample_Name)]
			row.names(ints_group)=row.names(ints)
	
			group.CNA.object <- CNA(data_group, CHR, MAPINFO,data.type = "logratio", sampleid = paste(paste(pd_group$Sample_Name,pd_group$Sample_Group),"qn"))
			group.smoothed.CNA.object <- smooth.CNA(group.CNA.object)
			group.segment.smoothed.CNA.object <- segment(group.smoothed.CNA.object, verbose = 1,alpha=0.001, undo.splits="sdundo", undo.SD=2)
		
			group.freq = glFrequency(group.segment.smoothed.CNA.object,freqThreshold)		
	
			#begin plot
			ints = ints[order(ints$CHR,ints$MAPINFO),]
			labels_chr <- data.matrix(summary(as.factor(ints$CHR)))

			test1<- data.frame(labels_chr,row.names(labels_chr) )
			test <- data.frame(unique(CHR))
			colnames(test) = c("chr")
			colnames(test1) = c("count","chr")
			F1 <- merge(test,test1, by="chr", sort=T)
			for(i in 2:length(row.names(F1))){F1[i,2] = F1[i-1,2] + F1[i,2] ; }

			F1$label <- NULL ; F1[1,3] <- F1[1,2] / 2 ;	
			for (i in 2:length(row.names(F1))){ F1[i,3] <- (F1[i,2]+F1[i-1,2])/2; }
	
			y1=group.freq$gain
			y2=group.freq$loss

			imageName1=paste(groups[g],"_","FreqPlot.pdf",sep="")
			imageName1=paste(newDir,imageName1,sep="/")
			graphTitle = paste("Frequency Plot of ",groups[g]," Samples",sep="")
		
			pdf(imageName1, width = 10.0, height = 9.0)
		
			#plot gain
			plot(y1, type='h',  xaxt="n",  yaxt="n", col="green", main = graphTitle , ylim=range(-1, 1), xlab='Chromosome Number',  ylab=paste('Fraction of Samples with Gain or Loss (n=',dim(data_group)[2],")",sep=""),xaxs = "i", yaxs = "i")

			#plot loss
			points(y2, type='h', col="red")

			#label for chromosomes
			x= c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1)
			y= c(1:length(F1[,2]))
			axis(1, at = c(F1[,2]), labels =FALSE, tick = TRUE, las = 1, col = "black", lty = "dotted", tck = 1 );
			axis(1, at = c(F1[,3]), labels =F1$chr, tick = FALSE );
			axis(2, at = c(x), labels=x, tick = TRUE, las = 1, col = "black", lty = "dotted", tck = 1 );
			dev.off()
		}
	}

}


