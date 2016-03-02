
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
		profileList_pos<-enviMass:::in_blind(profileList_pos)
		profileList_pos<<-profileList_pos
		save(profileList_pos,file=file.path(as.character(logfile[[1]]),"results","profileList_pos"));
		profpeaks_pos<-enviMass:::profiletopeak(profileList_pos,progbar=logfile$parameters[21])		
		profpeaks_pos<-profpeaks_pos[order(profpeaks_pos[,13],decreasing=TRUE),]
		profpeaks_pos<<-profpeaks_pos;
		save(profpeaks_pos,file=file.path(as.character(logfile[[1]]),"results","profpeaks_pos"));
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
		profileList_neg<-enviMass:::in_blind(profileList_neg)
		profileList_neg<<-profileList_neg
		save(profileList_neg,file=file.path(as.character(logfile[[1]]),"results","profileList_neg"));
		profpeaks_neg<-enviMass:::profiletopeak(profileList_neg,progbar=logfile$parameters[21])
		profpeaks_neg<-profpeaks_neg[order(profpeaks_neg[,13],decreasing=TRUE),]
		profpeaks_neg<<-profpeaks_neg;
		save(profpeaks_neg,file=file.path(as.character(logfile[[1]]),"results","profpeaks_neg"));
	}

