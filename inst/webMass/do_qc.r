
    measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	# positive #################################################################
	if(any(measurements[,"Mode"]=="positive") & (length(measurements[,"Mode"])>1)){
		filed<-list.files(file.path(logfile[[1]],"peaklist"))
		iles<-c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
		quant<-matrix(ncol=length(iles),nrow=length(filed),0)
		labelit<-c()
		at<-1
		rownames(quant)<-filed
		colnames(quant)<-as.character(iles)
		colorit<-sample(colors(),length(filed),replace=TRUE)
		for(i in 1:length(filed)){
			if(any(measurements[,"ID"]==filed[i])){
				if(  measurements[measurements[,"ID"]==filed[i],"Mode"]=="positive" ){
					if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
					if(any(objects()=="peaklist")){rm(peaklist)}
					load(file.path(logfile[[1]],"peaklist",filed[i]),envir=as.environment(".GlobalEnv"))
					int<-log10(peaklist[,3])
					int<-c(int-mean(int))
					int<-c(int/var(int))
					quant[at,]<-quantile(int,iles)
					at<-(at+1)
					labelit<-c(labelit,filed[i])
				}
			}
		}
		rm(peaklist,envir=as.environment(".GlobalEnv"))
		# first QC plot
		path=file.path(logfile[[1]],"pics","plotQCa_pos")
			png(filename = path, bg = "white")
			plot.new()
			plot.window(xlim=c(min(iles),max(iles)),ylim=c(min(quant),max(quant)))
			title(xlab="Quantiles",ylab="Norm. log intensity")
			for(n in 1:length(filed)){
				lines(iles,quant[n,],col=colorit[n],lwd=2)
				#lines(iles,quant[n,],col="darkgrey",lwd=2)				
			}
			box();axis(1);axis(2);
			#for(n in 1:length(filed)){
			#	if(maxit1[n]>=cut2|maxit2[n]>=cut1){
			#		lines(iles,quant[n,],col="red",lwd=2)
			#	}
			#}
		dev.off()
		expr1p<-list(src=file.path(logfile[[1]],"pics","plotQCa_pos"))
		output$plotQCa_pos<-renderImage(expr1p, deleteFile = FALSE)	
		# second QC plot	
		maxit1<-rep(0,length(filed))
		maxit2<-rep(0,length(filed))
		for(i in 1:length(labelit)){
		  atmax<-c();
		  for(j in 1:length(iles)){
			atmax<-c(atmax,abs(quant[i,j]-median(quant[-i,j])))
		  }
		  maxit1[i]<-max(atmax);
		  maxit2[i]<-median(atmax);
		}
		cut1<-(mean(maxit2)+6*sd(maxit2))
		cut2<-(mean(maxit1)+6*sd(maxit1))
		path=file.path(logfile[[1]],"pics","plotQCb_pos")
			png(filename = path, bg = "white")
			plot(maxit1,maxit2,pch=19,cex=1,col="white",
				xlab="Maximum quantile deviation",
				ylab="Median quantile deviation",
				xlim=c(0,max(c(cut2,maxit1))),
				ylim=c(0,max(c(cut1,maxit2)))	
			)
			abline(h=cut1,col="darkgrey",lty=2)
			abline(v=cut2,col="darkgrey",lty=2)
			text(maxit1,maxit2,labels=labelit,col=colorit)
			#points(maxit1,maxit2,pch=19,col="grey",cex=0.8)
			#points(maxit1[maxit1>=cut2|maxit2>=cut1],maxit2[maxit1>=cut2|maxit2>=cut1],pch=19,col="red")
		dev.off() 
		expr2p<-list(src=file.path(logfile[[1]],"pics","plotQCb_pos"))
		output$plotQCb_pos<-renderImage(expr2p, deleteFile = FALSE)
		# save results #############################################################
		# (1) quantile data
		qc_pos<-list(0)
		qc_pos[[1]]<-quant;
		qc_pos[[2]]<-maxit1;
		qc_pos[[3]]<-maxit2;
		save(qc_pos,file=file.path(logfile[[1]],"results","qc_pos"))
		# (2) mark measurements
		measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		measurements[,"checked"]<-"TRUE";
		filedcut<-filed[maxit1>cut2 & maxit2>cut1]
		if(length(filedcut)>0){
		  for(i in 1:length(filedcut)){
			measurements[measurements[,"ID"]==filedcut[i],"checked"]<-"FALSE";
			measurements[measurements[,"ID"]==filedcut[i],"include"]<-"FALSE";
		  }
		}
		write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
		measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	}
	# negative #################################################################
	if(any(measurements[,"Mode"]=="negative") & (length(measurements[,"Mode"])>1)){
		filed<-list.files(file.path(logfile[[1]],"peaklist"))
		iles<-c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
		quant<-matrix(ncol=length(iles),nrow=length(filed),0)
		labelit<-c();
		at<-1
		rownames(quant)<-filed
		colnames(quant)<-as.character(iles)
		colorit<-sample(colors(),length(filed),replace=TRUE)
		for(i in 1:length(filed)){
			if(any(measurements[,"ID"]==filed[i])){
				if(  measurements[measurements[,"ID"]==filed[i],"Mode"]=="negative" ){
					if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
					if(any(objects()=="peaklist")){rm(peaklist)}
					load(file.path(logfile[[1]],"peaklist","",filed[i]),envir=as.environment(".GlobalEnv"))
					int<-log10(peaklist[,3])
					int<-c(int-mean(int))
					int<-c(int/var(int))
					quant[at,]<-quantile(int,iles)
					at<-c(at+1)
					labelit<-c(labelit,filed[i])
				}
			}
		}
		rm(peaklist,envir=as.environment(".GlobalEnv"))
		# first QC plot
		path=file.path(logfile[[1]],"pics","plotQCa_neg")
			png(filename = path, bg = "white")
			plot.new()
			plot.window(xlim=c(min(iles),max(iles)),ylim=c(min(quant),max(quant)))
			title(xlab="Quantiles",ylab="Norm. log intensity")
			for(n in 1:length(filed)){
				lines(iles,quant[n,],col=colorit[n],lwd=2)
			}
			box();axis(1);axis(2);
		dev.off()
		expr1n<-list(src=file.path(logfile[[1]],"pics","plotQCa_neg"))
		output$plotQCa_neg<-renderImage(expr1n, deleteFile = FALSE)	
		# second QC plot	
		maxit1<-rep(0,length(filed))
		maxit2<-rep(0,length(filed))
		for(i in 1:length(labelit)){
		  atmax<-c();
		  for(j in 1:length(iles)){
			atmax<-c(atmax,abs(quant[i,j]-median(quant[-i,j])))
		  }
		  maxit1[i]<-max(atmax);
		  maxit2[i]<-median(atmax);
		}
		cut1<-(mean(maxit2)+6*sd(maxit2))
		cut2<-(mean(maxit1)+6*sd(maxit1))
		path=file.path(logfile[[1]],"pics","plotQCb_neg")
			png(filename = path, bg = "white")
			plot(maxit1,maxit2,pch=19,cex=1,xlab="Maximum quantile deviation",ylab="Median quantile deviation",col="white")
			abline(h=cut1,col="darkgrey",lty=2)
			abline(v=cut2,col="darkgrey",lty=2)
			text(maxit1,maxit2,labels=labelit,col=colorit)
		dev.off() 
		expr2n<-list(src=file.path(logfile[[1]],"pics","plotQCb_neg"))
		output$plotQCb_neg<-renderImage(expr2n, deleteFile = FALSE)
		# save results #############################################################
		# (1) quantile data
		qc_neg<-list(0)
		qc_neg[[1]]<-quant;
		qc_neg[[2]]<-maxit1;
		qc_neg[[3]]<-maxit2;
		save(qc_neg,file=file.path(logfile[[1]],"results","qc_neg"))
		# (2) mark measurements
		measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		measurements[,"checked"]<-"TRUE";
		filedcut<-filed[maxit1>cut2 & maxit2>cut1]
		if(length(filedcut)>0){
		  for(i in 1:length(filedcut)){
			measurements[measurements[,"ID"]==filedcut[i],"checked"]<-"FALSE";
			measurements[measurements[,"ID"]==filedcut[i],"include"]<-"FALSE";
		  }
		}
		write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
		measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	}
