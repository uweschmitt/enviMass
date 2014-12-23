if(
	(logfile$workflow[9]=="yes" && logfile$summary[8,2]=="FALSE") || (logfile$Tasks_to_redo[6]=="TRUE") 
){

    measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	if(any(measurements[,4]=="positive")){
		if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="peaklist")){rm(peaklist)}
		profileList_pos<-startprofiles(
							logfile,
							frac=FALSE,
							sets=as.numeric(logfile$parameters$prof_maxfiles),
							progbar=logfile$parameters$progressBar,
							ion_mode="positive",
							until=logfile$parameters$upto_file
						)
		profileList_pos<-agglomer(
							profileList_pos,
							dmass=(as.numeric(logfile$parameters$prof_dmz)+1),
							ppm=as.logical(as.character(logfile$parameters$prof_ppm)),
							dret=(as.numeric(logfile$parameters$prof_drt)+10)
						)
		profileList_pos<-partcluster(
							profileList_pos,
							dmass=as.numeric(logfile$parameters$prof_dmz),
							ppm=as.logical(as.character(logfile$parameters$prof_ppm)),
							dret=as.numeric(logfile$parameters$prof_drt),
							from=FALSE,
							to=FALSE,
							progbar=logfile$parameters$progressBar,
							plotit=FALSE
						)
		profileList_pos<<-profileList_pos
		save(profileList_pos,file=file.path(as.character(logfile[[1]]),"results","profileList_pos"));
	}
	if(any(measurements[,4]=="negative")){
		if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="peaklist")){rm(peaklist)}
		profileList_neg<-startprofiles(
							logfile,
							frac=FALSE,
							sets=as.numeric(logfile$parameters$prof_maxfiles),
							progbar=logfile$parameters$progressBar,
							ion_mode="negative",
							until=logfile$parameters$upto_file
						)
		profileList_neg<-agglomer(
							profileList_neg,
							dmass=(as.numeric(logfile$parameters$prof_dmz)+1),
							ppm=as.logical(as.character(logfile$parameters$prof_ppm)),
							dret=(as.numeric(logfile$parameters$prof_drt)+10)
						)
		profileList_neg<-partcluster(
							profileList_neg,
							dmass=as.numeric(logfile$parameters$prof_dmz),
							ppm=as.logical(as.character(logfile$parameters$prof_ppm)),
							dret=as.numeric(logfile$parameters$prof_drt),
							from=FALSE,
							to=FALSE,
							progbar=logfile$parameters$progressBar,
							plotit=FALSE
						)
		profileList_neg<<-profileList_neg
		save(profileList_neg,file=file.path(as.character(logfile[[1]]),"results","profileList_neg"));
	}
	logfile$summary[8,2]<<-"TRUE";
    logfile$summary[8,2]<-"TRUE";
	logfile$Tasks_to_redo[6]<-"FALSE";
	logfile$Tasks_to_redo[6]<<-"FALSE";
    save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
    summa[8,2]<-"done"
	summa[8,2]<<-"done"
	output$summar<<-renderTable(summa)
    cat("Profile extraction done \n");
    output$dowhat<<-renderText("Profile extraction ... wait")

}else{

	if(logfile$workflow[9]=="no"){
			logfile$summary[8,2]<<-"FALSE";
			logfile$summary[8,2]<-"FALSE";
	}
	logfile$Tasks_to_redo[6]<-"FALSE";
	logfile$Tasks_to_redo[6]<<-"FALSE";	
    save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
    summa[8,2]<-"skipped";
    summa[8,2]<<-"skipped";
	output$summar<<-renderTable(summa)
    cat("Profile extraction skipped \n");
    output$dowhat<<-renderText("Profile extraction skipped ... wait")

}