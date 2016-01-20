
if(  
	#(logfile$workflow[2]=="yes" && logfile$summary[5,2]=="FALSE")  || 
	#(logfile$Tasks_to_redo[3]=="TRUE") 
	FALSE
){

    ############################################################################
	# for IS compounds #########################################################
	# positive ionization ######################################################
	if(TRUE){
	
		if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="peaklist")){rm(peaklist)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_pos")){rm(profileList_pos,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="profileList_pos")){rm(profileList_pos)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="pattern_pos_IS")){rm(pattern_pos_IS,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="pattern_pos_IS")){rm(pattern_pos_IS)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="patternRT_pos_IS")){rm(patternRT_pos_IS,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="patternRT_pos_IS")){rm(patternRT_pos_IS)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="patternDelRT_pos_IS")){rm(patternDelRT_pos_IS,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="patternDelRT_pos_IS")){rm(patternDelRT_pos_IS)}
		load(file=file.path(as.character(logfile[[1]]),"results","profileList_pos"),envir=as.environment(".GlobalEnv"));	
		load(file=file.path(logfile[[1]],"results","pattern_pos_IS"),envir=as.environment(".GlobalEnv"));
		load(file=file.path(logfile[[1]],"results","patternRT_pos_IS"),envir=as.environment(".GlobalEnv"));
		load(file=file.path(logfile[[1]],"results","patternDelRT_pos_IS"),envir=as.environment(".GlobalEnv"));
		
		
		peaks<-profileList_pos[[7]];
		peaklist<-peaks[,c(14,16,15)];
				
		# screen centroids
		count_nonmax<-0
		for(i in 1:length(getit_1)){
			count_nonmax<-(count_nonmax+
				length(pattern_pos_IS[[i]][,1])
			)
		}
		other_peak_mass<-rep(0,count_nonmax)
		other_peak_ID<-rep(0,count_nonmax)
		other_peak_number<-rep(0,count_nonmax)
		other_peak_RT<-rep(0,count_nonmax)
		other_peak_dRT<-rep(0,count_nonmax)
		at_ID<-1
		for(i in 1:length(getit_1)){
			if(getit_1[i]!="FALSE"){
				n<-length(pattern_pos_IS[[i]][,1])
				if(n>0){ # more than one centroid per compound?
					other_peak_mass[at_ID:(at_ID+n-1)]<-pattern_pos_IS[[i]][,1]
					other_peak_ID[at_ID:(at_ID+n-1)]<-i
					other_peak_number[at_ID:(at_ID+n-1)]<-(1:n)
					other_peak_RT[at_ID:(at_ID+n-1)]<-patternRT_pos_IS[i]
					other_peak_dRT[at_ID:(at_ID+n-1)]<-patternDelRT_pos_IS[i]
					at_ID<-(at_ID+n)
				}
			}
		}
		getit <- search_peak(
			peaklist, 
			other_peak_mass, 
			dmz=5, 
			ppm=TRUE, 
			RT = other_peak_RT, 
			dRT=other_peak_dRT)	





			
		
	
	
	}
	
	
    ############################################################################
	# then on targets pos, before swithcing to negative list - upload takes too long
	
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
	output$summar<<-renderTable(summa);
    cat("Mass recalibration skipped \n");
    output$dowhat<<-renderText("Recalibration skipped ... wait")
}
}



