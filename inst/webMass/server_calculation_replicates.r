
if(  
	#(logfile$workflow[2]=="yes" && logfile$summary[5,2]=="FALSE")  || 
	#(logfile$Tasks_to_redo[3]=="TRUE") 
	FALSE
){

    measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	# define replicate groups by tag3 entries ###############################################################
	if(any(duplicated(measurements$tag3[measurements$tag3!="FALSE"])){
	
		replic<-(measurements$tag3[measurements$tag3!="FALSE"])
		replic<-replic[duplicated(replic)]
		replic<-unique(replic)
		for(i in 1:length(replic)){
			IDs<-measurements$ID[measurements$tag3==replic[i]]
			if(any(duplicated(IDs))){stop("replicates: non-unique IDs found!")} # should not happen anyway
	
			replicated<-list()
			for(j in 1:length(IDs)){
				load(file=file.path(logfile[[1]],"peaklist",as.character(IDs[j])),envir=as.environment(".GlobalEnv"),verbose=FALSE);
				if(ppm){
				
				
				}else{
				
				
				}
				
				
				
				peaklist[,c(12,13,14)]
				
				
				replicated[[j]]<-rep(FALSE,length(peaklist[,1]))
			}
			for(j in 1:length()){
			
	
			}
	
	
	
	
		}
	
	
	
	
	
	}
	





}else{
if(FALSE){
	if(logfile$workflow[2]=="no"){
		logfile$summary[5,2]<<-"FALSE";
		logfile$summary[5,2]<-"FALSE";
		path=file.path(logfile[[1]],"pics","recal_none")
			png(filename = path, bg = "white")
			plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
			dev.off()
		    exprrec<-list(src=path)
			output$recal_pic<-renderImage(exprrec, deleteFile = FALSE);		
			output$peakhist_pic<-renderImage(exprrec, deleteFile = FALSE);
			output$peakmzRT_pic<-renderImage(exprrec, deleteFile = FALSE);	
	}
	logfile$Tasks_to_redo[3]<-"FALSE";
	logfile$Tasks_to_redo[3]<<-"FALSE";
    save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
    summa[5,2]<-"skipped"
	summa[5,2]<<-"skipped"
	output$summa_html<<-renderText(summary_html(summa));
    cat("Mass recalibration skipped \n");
    output$dowhat<<-renderText("Recalibration skipped ... wait")
}
}




