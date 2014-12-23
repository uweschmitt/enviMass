#################################################################################
# get LOD- & 90percentile-Intensity for picked peaks ############################
those<-list.files(file.path(logfile$project_folder,"peaklist", fsep = "\\"))
for(i in 1:length(those)){
	if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="peaklist")){rm(peaklist)}
	load(file.path(logfile$project_folder,"peaklist",those[i], fsep = "\\"),envir=as.environment(".GlobalEnv"))
	#plot(peaklist[,5],log10(peaklist[,13]),pch=19,cex=0.6,col="darkgrey",xlab="RT",
	#	ylab=expression(log[10]*paste(" Intensity",sep=" "))
	#)
	his<-hist(peaklist[,5],breaks=100,plot=FALSE)
	#abline(v=his$breaks,col="grey")
	# LOD #######################################################################
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
	#points(get_ret,get_int,col="black",pch=19,cex=0.5)
	model<-gam(get_int~s(get_ret),weights=get_w,sp=.1)
	assign(paste("LOD_",those[i],sep=""),model)
	#lines(get_ret,predict(get(paste("LOD_",those[i],sep=""))),col="red",lwd=3)
	save(
		list=paste("LOD_",those[i],sep="")[[1]],
		file=file.path(logfile$project_folder,"results",paste("LOD_",those[i],sep=""),fsep = "\\"),
		precheck=FALSE
	)
	# 90% ######################################################################
	get_int<-c()
	get_ret<-c()
	get_w<-c()
	for(j in 2:length(his$breaks)){
		ret<-peaklist[(peaklist[,5]>=his$breaks[j-1] & peaklist[,5]<his$breaks[j]),5]
		int<-log10(peaklist[(peaklist[,5]>=his$breaks[j-1] & peaklist[,5]<his$breaks[j]),13])
		ret<-ret[order(int,decreasing=FALSE)]
		int<-int[order(int,decreasing=FALSE)]
		getit<-ceiling(length(int)*0.95)
		get_int<-c(get_int,int[getit])
		get_ret<-c(get_ret,ret[getit])
		if(length(ret)>0){get_w<-c(get_w,length(int))}
	}
	#points(get_ret,get_int,col="black",pch=19,cex=0.5)
	model<-gam(get_int~s(get_ret),weights=get_w,sp=1)
	assign(paste("90perc_",those[i],sep=""),model)
	#lines(get_ret,predict(get(paste("90perc_",those[i],sep=""))),col="darkgreen",lwd=3)
	save(
		list=paste("90perc_",those[i],sep="")[[1]],
		file=file.path(logfile$project_folder,"results",paste("90perc_",those[i],sep=""),fsep = "\\"),
		precheck=FALSE
	)
	#Sys.sleep(2)	
}
#################################################################################
# load all LOD- & 90%-gams ######################################################
those<-list.files(file.path(logfile$project_folder,"results", fsep = "\\"))
for(i in 1:length(those)){
	if(grepl("LOD_",those[i])){
		load(file=file.path(logfile$project_folder,"results",those[i],fsep = "\\"),envir=as.environment(".GlobalEnv"))
		next;
	}
	if(grepl("90perc_",those[i])){
			load(file=file.path(logfile$project_folder,"results",those[i],fsep = "\\"),envir=as.environment(".GlobalEnv"))
	}
}
#################################################################################
	