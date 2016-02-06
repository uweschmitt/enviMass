####################################################################
# delete old LOD gams & results ####################################
those<-list.files(file.path(logfile$project_folder,"results","LOD"))
if(length(those)>0){
	for(i in 1:length(those)){
		file.remove(file.path(logfile$project_folder,"results","LOD",those[i]))
	}
}

#################################################################################
# get LOD- & 90percentile-Intensity for picked peaks ############################
those<-list.files(file.path(logfile$project_folder,"peaklist", fsep = "\\"))

# MAKE this run file-wise
if(length(those)>0){

	LOD_splined<-list()
	at<-1;
	cat("\n")
	for(i in 1:length(those)){
		if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="peaklist")){rm(peaklist)}
		load(file.path(logfile$project_folder,"peaklist",those[i]),envir=as.environment(".GlobalEnv"))
		peaklist<-peaklist[peaklist[,colnames(peaklist)=="keep"]==1,,drop=FALSE]
		if(length(peaklist[,1])==0){next}
		# LOD ###################################################################
		get_int<-c()
		get_ret<-c()
		get_w<-c()
		for(j in 2:length(his$breaks)){
			ret<-peaklist[(peaklist[,5]>=his$breaks[j-1] & peaklist[,5]<his$breaks[j]),5]
			int<-log10(peaklist[(peaklist[,5]>=his$breaks[j-1] & peaklist[,5]<his$breaks[j]),13])
			ret<-ret[order(int,decreasing=FALSE)]
			int<-int[order(int,decreasing=FALSE)]
			getit<-ceiling(length(int)*0.1)
			get_int<-c(get_int,int[getit])
			get_ret<-c(get_ret,ret[getit])
			if(length(ret)>0){get_w<-c(get_w,length(int))}
		}
		model<-smooth.spline(x=get_ret,y=get_int)	
		assign(paste("LOD_",those[i],sep=""),model);rm(model)
		LOD_splined[[at]]<-get(paste("LOD_",those[i],sep=""))
		names(LOD_splined[at])<-paste("LOD_",those[i],sep="")
		at<-(at+1)
		if(FALSE){
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
			dev.off();cat(".")
		}
	}
	save(LOD_splined,file=file.path(logfile$project_folder,"results","LOD","LOD_splined"))

}
#################################################################################

