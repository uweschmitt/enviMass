#################################################################################
# get LOD- & 90percentile-Intensity for picked peaks ############################
those<-list.files(file.path(logfile$project_folder,"peaklist", fsep = "\\"))
measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");

# MAKE this run file-wise
if(length(those)>0){
	LOD_splined_new<-list()
	if(file.exists(file.path(logfile$project_folder,"results","LOD","LOD_splined"))){ # load existing model
		load(file.path(logfile$project_folder,"results","LOD","LOD_splined"))
	}
	at<-1;
	cat(" | ")
	for(i in 1:length(those)){
		if(!any(measurements[,"ID"]==those[i])){cat("\n orphaned peaklist found.");next;} # not in list of measurements?
		if(measurements[measurements[,"ID"]==those[i],"include"]=="FALSE"){next}
		if(measurements[measurements[,"ID"]==those[i],"LOD"]=="TRUE"){ # ONLY valid for old files! Not for newly loaded ones!
			if(any(names(LOD_splined)==paste("LOD_",those[i],sep=""))){ # copy old model, already done
				copy_this<-which(names(LOD_splined)==paste("LOD_",those[i],sep=""))
				LOD_splined_new[[at]]<-LOD_splined[[copy_this]]
				names(LOD_splined_new)[at]<-paste("LOD_",those[i],sep="")
				at<-(at+1)
				cat("\nCopied LOD model.")
				next;
			}else{
				cat("\nMake new LOD model.")
			}
		} 		
		if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="peaklist")){rm(peaklist)}
		load(file.path(logfile$project_folder,"peaklist",those[i]),envir=as.environment(".GlobalEnv"))
		#peaklist<-peaklist[peaklist[,colnames(peaklist)=="keep"]==1,,drop=FALSE]
		if(length(peaklist[,1])==0){next}
		# LOD ###################################################################
		his<-hist(peaklist[,"RT"],breaks=100,plot=FALSE)
		get_int<-c()
		get_ret<-c()
		get_w<-c()
		for(j in 2:length(his$breaks)){
			ret<-peaklist[(peaklist[,"RT"]>=his$breaks[j-1] & peaklist[,"RT"]<his$breaks[j]),"RT"]
			int<-log10(peaklist[(peaklist[,"RT"]>=his$breaks[j-1] & peaklist[,"RT"]<his$breaks[j]),"int_corr"])
			ret<-ret[order(int,decreasing=FALSE)]
			int<-int[order(int,decreasing=FALSE)]
			getit<-ceiling(length(int)*0.1)
			get_int<-c(get_int,int[getit])
			get_ret<-c(get_ret,ret[getit])
			if(length(ret)>0){get_w<-c(get_w,length(int))}
		}
		model<-smooth.spline(x=get_ret,y=get_int)	
		assign(paste("LOD_",those[i],sep=""),model);rm(model)
		LOD_splined_new[[at]]<-get(paste("LOD_",those[i],sep=""))
		names(LOD_splined_new)[at]<-paste("LOD_",those[i],sep="")
		at<-(at+1)
		if(TRUE){
			png(file=file.path(logfile$project_folder,"results","LOD",paste("plot_LOD_",those[i],".png",sep="")),
				width = 550, height = 250)
				par(mar=c(4.2,4.2,0.2,0.2))
				plot(peaklist[,5],log10(peaklist[,13]),pch=19,cex=0.6,col="darkgrey",xlab="RT",
					ylab=expression(log[10]*paste(" Intensity",sep=" "))
				)
				his<-hist(peaklist[,5],breaks=100,plot=FALSE)
				abline(v=his$breaks,col="grey")
				lines(get_ret,predict(get(paste("LOD_",those[i],sep="")))$y,col="red",lwd=2)
				points(get_ret,get_int,col="black",pch=19,cex=0.5)
				box()
			dev.off();
			cat("!")
		}
		measurements[measurements[,"ID"]==those[i],"LOD"]<-"TRUE";
	}
	LOD_splined<-LOD_splined_new
	save(LOD_splined,file=file.path(logfile$project_folder,"results","LOD","LOD_splined"))
	write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
}
#################################################################################

