measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
types<-measurements[,"Type"]
IDs<-list.files(file.path(logfile[[1]],"peaklist"))
if(length(IDs)>0){
	for(i in 1:length(IDs)){
		if(any(measurements[,"ID"]==IDs[i])){
			if(types[measurements[,"ID"]==IDs[i]]=="sample"){
				load(file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])),envir=as.environment(".GlobalEnv"),verbose=FALSE);
				keep_2<-rep(1,length(peaklist[,1])) # 1 == TRUE
				peaklist[,colnames(peaklist)=="keep_2"]<-keep_2
				save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])))
				rm(peaklist)
			}
		}else{
			cat("\n Orphaned peaklist detected - from an older workflow run?")
		}
	}
}
