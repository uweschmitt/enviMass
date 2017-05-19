	
	############################################################################################
	# REMOVE OLD RESULTS #######################################################################
	if(file.exists(file.path(as.character(logfile[[1]]),"results","profileList_pos"))){file.remove(file.path(as.character(logfile[[1]]),"results","profileList_pos"))}
	if(file.exists(file.path(as.character(logfile[[1]]),"results","profileList_pos_copy"))){file.remove(file.path(as.character(logfile[[1]]),"results","profileList_pos_copy"))}
	if(file.exists(file.path(as.character(logfile[[1]]),"results","profileList_neg"))){file.remove(file.path(as.character(logfile[[1]]),"results","profileList_neg"))}
	if(file.exists(file.path(as.character(logfile[[1]]),"results","profileList_neg_copy"))){file.remove(file.path(as.character(logfile[[1]]),"results","profileList_neg_copy"))}
	
	############################################################################################
	# FILTER FILES #############################################################################
    measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	measurements<-measurements[measurements[,"include"]=="TRUE",,drop=FALSE]
	if(logfile$parameters$prof_select=="TRUE"){
		measurements<-measurements[measurements[,names(measurements)=="profiled"]=="TRUE",]	
	}
	############################################################################################

	
	if(any(measurements[,"Mode"]=="positive")){
		if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="peaklist")){rm(peaklist)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_pos")){rm(profileList_pos,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="profileList_pos")){rm(profileList_pos)}		
		if(file.exists(file.path(as.character(logfile[[1]]),"results","profileList_pos"))){file.remove(file.path(as.character(logfile[[1]]),"results","profileList_pos"))}
		if(file.exists(file.path(as.character(logfile[[1]]),"results","profileList_pos_copy"))){file.remove(file.path(as.character(logfile[[1]]),"results","profileList_pos_copy"))}
		profileList_pos<-startprofiles(
							logfile,
							frac=FALSE,
							sets=as.numeric(logfile$parameters$prof_maxfiles),
							progbar=logfile$parameters$progressBar,
							ion_mode="positive",
							until=logfile$parameters$upto_file,
							selective=logfile$parameters$prof_select,
							types=c("sample","blank","spiked")
						)
		if(any(profileList_pos[[2]][,2]==0)){stop("\n issue in do_profiling: zero intensities detected - resolve issue!")}
		profileList_pos<-agglomer(
							profileList_pos,
							dmass=(as.numeric(logfile$parameters$prof_dmz)+1),
							ppm=as.logical(as.character(logfile$parameters$prof_ppm)),
							dret=(as.numeric(logfile$parameters$prof_drt)+10)
						)
		if((logfile$workflow[names(logfile$workflow)=="replicates"]=="no")||(logfile$parameters$replicates_prof=="no")){ 				
			profileList_pos<-partcluster(
								profileList=profileList_pos,
								dmass=as.numeric(logfile$parameters$prof_dmz),
								ppm=as.logical(as.character(logfile$parameters$prof_ppm)),
								dret=as.numeric(logfile$parameters$prof_drt),
								from=FALSE,
								to=FALSE,
								progbar=logfile$parameters$progressBar,
								#progbar=TRUE,
								plot_it=FALSE,
								replicates=FALSE,
								IDs=FALSE,
								with_test=FALSE
							)
		}else{ # run a profiling in the replicate groups first
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			measurements<-measurements[measurements[,"include"]=="TRUE",]
			measurements<-measurements[measurements[,"Mode"]=="positive",]
			replicates<-measurements$tag3
			IDs<-measurements$ID
			profileList_pos<-partcluster(
								profileList=profileList_pos,
								dmass=as.numeric(logfile$parameters$prof_dmz),
								ppm=as.logical(as.character(logfile$parameters$prof_ppm)),
								dret=as.numeric(logfile$parameters$prof_drt),
								from=FALSE,
								to=FALSE,
								progbar=logfile$parameters$progressBar,
								plot_it=FALSE,
								replicates=replicates,
								IDs
							)
		}
		profileList_pos<-enviMass:::in_blind(profileList_pos)
		profileList_pos<<-profileList_pos
		save(profileList_pos,file=file.path(as.character(logfile[[1]]),"results","profileList_pos"),compress=FALSE);
		profileList_pos_copy<-profileList_pos
		save(profileList_pos_copy,file=file.path(as.character(logfile[[1]]),"results","profileList_pos_copy"),compress=FALSE); # used for screening - does not include modifications of downstream compound subtraction		
		profpeaks_pos<-enviMass:::profiletopeak(profileList_pos,progbar=logfile$parameters$progressBar)		
		profpeaks_pos<-profpeaks_pos[order(profpeaks_pos[,13],decreasing=TRUE),];
		profpeaks_pos<<-profpeaks_pos;
		save(profpeaks_pos,file=file.path(as.character(logfile[[1]]),"results","profpeaks_pos"));
		links_peaks_pos<-list(); # each entry with 6 lists itself: targets, IS, EIC_correl, isotop, adducts, homol
		save(links_peaks_pos,file=file.path(as.character(logfile[[1]]),"results","links_peaks_pos"));				
	}
	if(any(measurements[,"Mode"]=="negative")){
		if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="peaklist")){rm(peaklist)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_neg")){rm(profileList_neg,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="profileList_neg")){rm(profileList_neg)}	
		profileList_neg<-startprofiles(
							logfile,
							frac=FALSE,
							sets=as.numeric(logfile$parameters$prof_maxfiles),
							progbar=logfile$parameters$progressBar,
							ion_mode="negative",
							until=logfile$parameters$upto_file,
							selective=logfile$parameters$prof_select,
							types=c("sample","blank","spiked")
						)
		if(any(profileList_neg[[2]][,2]==0)){stop("\n issue in do_profiling: zero intensities detected - resolve issue!")}
		profileList_neg<-agglomer(
							profileList_neg,
							dmass=(as.numeric(logfile$parameters$prof_dmz)+1),
							ppm=as.logical(as.character(logfile$parameters$prof_ppm)),
							dret=(as.numeric(logfile$parameters$prof_drt)+10)
						)
		if((logfile$workflow[names(logfile$workflow)=="replicates"]=="no")||(logfile$parameters$replicates_prof=="no")){ 				
			profileList_neg<-partcluster(
								profileList=profileList_neg,
								dmass=as.numeric(logfile$parameters$prof_dmz),
								ppm=as.logical(as.character(logfile$parameters$prof_ppm)),
								dret=as.numeric(logfile$parameters$prof_drt),
								from=FALSE,
								to=FALSE,
								progbar=logfile$parameters$progressBar,
								plot_it=FALSE,
								replicates=FALSE,
								IDs=FALSE
							)
		}else{ # run a profiling in the replicate groups first
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			measurements<-measurements[measurements[,"include"]=="TRUE",]
			measurements<-measurements[measurements[,"Mode"]=="negative",]
			replicates<-measurements$tag3
			IDs<-measurements$ID
			profileList_neg<-partcluster(
								profileList=profileList_neg,
								dmass=as.numeric(logfile$parameters$prof_dmz),
								ppm=as.logical(as.character(logfile$parameters$prof_ppm)),
								dret=as.numeric(logfile$parameters$prof_drt),
								from=FALSE,
								to=FALSE,
								progbar=logfile$parameters$progressBar,
								plot_it=FALSE,
								replicates=replicates,
								IDs
							)
		}
		profileList_neg<-enviMass:::in_blind(profileList_neg)
		profileList_neg<<-profileList_neg
		save(profileList_neg,file=file.path(as.character(logfile[[1]]),"results","profileList_neg"),compress=FALSE);
		profileList_neg_copy<-profileList_neg
		save(profileList_neg_copy,file=file.path(as.character(logfile[[1]]),"results","profileList_neg_copy"),compress=FALSE); # used for screening - does not include modifications of downstream compound subtraction			
		profpeaks_neg<-enviMass:::profiletopeak(profileList_neg,progbar=logfile$parameters$progressBar)
		profpeaks_neg<-profpeaks_neg[order(profpeaks_neg[,13],decreasing=TRUE),]
		profpeaks_neg<<-profpeaks_neg;
		save(profpeaks_neg,file=file.path(as.character(logfile[[1]]),"results","profpeaks_neg"));
		links_peaks_neg<-list(); # each entry with 6 lists itself: targets, IS, EIC_correl, isotop, adducts, homol
		save(links_peaks_neg,file=file.path(as.character(logfile[[1]]),"results","links_peaks_neg"));		
	}

