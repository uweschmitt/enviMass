
measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
measurements<-measurements[measurements[,8]=="TRUE",]
replic<-(measurements$tag3[measurements$tag3!="FALSE"])
replic<-replic[duplicated(replic)]
replic<-unique(replic)
if(length(replic)>0){
	for(i in 1:length(replic)){
		IDs<-measurements$ID[measurements$tag3==replic[i]]
		for(j in 1:length(IDs)){
			if(file.exists(file.path(logfile[[1]],"peaklist",as.character(IDs[j])))){
				load(file=file.path(logfile[[1]],"peaklist",as.character(IDs[j])),envir=as.environment(".GlobalEnv"),verbose=FALSE);
				keep<-rep(1,length(peaklist[,1]))
				peaklist[,colnames(peaklist)=="keep"]<-keep
				save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(IDs[j])))
				rm(peaklist)
			}else{
				cat("Missing peaklist detected in done_replicates. Debug!?")
			}
		}
	}
}

