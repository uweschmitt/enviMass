#######################################################################################################
# POSITIVE IONIZATION #################################################################################
if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_pos")){rm(profileList_pos,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="profileList_pos")){rm(profileList_pos)}
if(any(objects(envir=as.environment(".GlobalEnv"))=="links_peaks_pos")){rm(links_peaks_pos,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="links_peaks_pos")){rm(links_peaks_pos)}
if(
	file.exists(file.path(as.character(logfile[[1]]),"results","profileList_pos")) 
){

	load(file=file.path(as.character(logfile[[1]]),"results","profileList_pos"),envir=as.environment(".GlobalEnv"));
	if( any(profileList_pos[[2]][,9]==0) ){ # if any blind peaks exist	
		###############################################################################################
		keep_peaks<-rep(TRUE,length(profileList_neg[[2]][,1]))
		keep_peaks[profileList_neg[[2]][,9]==0]<-FALSE
		profileList_pos[[2]]<-profileList_pos[[2]][keep_peaks,]
		profileList_pos[[6]]<-0  # index_agglom
		profileList_pos[[7]]<-0  # index_prof
		profileList_pos[[8]]<-0  # parameters
		profileList_pos[[1]]<-data.frame(TRUE,FALSE,FALSE,FALSE)    # state
		colnames(profileList_pos[[1]])<-c("peaks?","agglom?","profiling","trends?")
		###############################################################################################
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
		profileList_pos<-enviMass:::in_blind(profileList_pos) # mark profile if from a blank/blind file
		profileList_pos<<-profileList_pos
		save(profileList_pos,file=file.path(as.character(logfile[[1]]),"results","profileList_pos"));
		profpeaks_pos<-enviMass:::profiletopeak(profileList_pos,progbar=logfile$parameters[21])		
		profpeaks_pos<-profpeaks_pos[order(profpeaks_pos[,13],decreasing=TRUE),];
		profpeaks_pos<<-profpeaks_pos;
		save(profpeaks_pos,file=file.path(as.character(logfile[[1]]),"results","profpeaks_pos"));
		save(links_peaks_pos,file=file.path(as.character(logfile[[1]]),"results","links_peaks_pos"));	
		cat(paste("\nBlind in profile subtraction, positive: ",round((sum(!keep_peaks)/length(keep_peaks)*100),digits=3),"% of peaks removed",sep=""))
		###############################################################################################
		save(profileList_pos,file=file.path(as.character(logfile[[1]]),"results","profileList_pos"))
	}else{
		cat("\n do_blind_subtr: nothing to subtract, positive ioniz.")
	}
}
#######################################################################################################

#######################################################################################################
# NEGATIVE IONIZATION #################################################################################
if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_neg")){rm(profileList_pos,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="profileList_neg")){rm(profileList_pos)}
if(any(objects(envir=as.environment(".GlobalEnv"))=="links_peaks_neg")){rm(links_peaks_pos,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="links_peaks_neg")){rm(links_peaks_neg)}
if(
	file.exists(file.path(as.character(logfile[[1]]),"results","profileList_neg"))
){

	load(file=file.path(as.character(logfile[[1]]),"results","profileList_neg"),envir=as.environment(".GlobalEnv"));	
	if( any(profileList_neg[[2]][,9]==0) ){ # if any links exist
		###############################################################################################
		keep_peaks<-rep(TRUE,length(profileList_neg[[2]][,1]))
		keep_peaks[profileList_neg[[2]][,9]==0]<-FALSE
		profileList_neg[[2]]<-profileList_neg[[2]][keep_peaks,]
		profileList_neg[[6]]<-0  # index_agglom
		profileList_neg[[7]]<-0  # index_prof
		profileList_neg[[8]]<-0  # parameters
		profileList_neg[[1]]<-data.frame(TRUE,FALSE,FALSE,FALSE)    # state
		colnames(profileList_neg[[1]])<-c("peaks?","agglom?","profiling","trends?")
		###############################################################################################
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
		profileList_neg<-enviMass:::in_blind(profileList_neg) # mark profile if from a blank/blind file
		profileList_neg<<-profileList_neg
		save(profileList_neg,file=file.path(as.character(logfile[[1]]),"results","profileList_neg"));
		profpeaks_neg<-enviMass:::profiletopeak(profileList_neg,progbar=logfile$parameters[21])		
		profpeaks_neg<-profpeaks_neg[order(profpeaks_neg[,13],decreasing=TRUE),];
		profpeaks_neg<<-profpeaks_neg;
		save(profpeaks_neg,file=file.path(as.character(logfile[[1]]),"results","profpeaks_neg"));
		save(links_peaks_neg,file=file.path(as.character(logfile[[1]]),"results","links_peaks_neg"));	
		cat(paste("\nBlind in profile subtraction, negative: ",round((sum(!keep_peaks)/length(keep_peaks)*100),digits=3),"% of peaks removed",sep=""))
		###############################################################################################
		save(profileList_neg,file=file.path(as.character(logfile[[1]]),"results","profileList_neg"))
	}else{
		cat("\n do_blind_subtr: nothing to subtract, negative ioniz.")
	}
}
#######################################################################################################



