#######################################################################################################
# POSITIVE IONIZATION #################################################################################
if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_pos")){rm(profileList_pos,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="profileList_pos")){rm(profileList_pos)}
if(any(objects(envir=as.environment(".GlobalEnv"))=="links_peaks_pos")){rm(links_peaks_pos,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="links_peaks_pos")){rm(links_peaks_pos)}
if(
	file.exists(file.path(as.character(logfile[[1]]),"results","profileList_pos")) &
	file.exists(file.path(as.character(logfile[[1]]),"results","links_peaks_pos"))	
){

	load(file=file.path(as.character(logfile[[1]]),"results","links_peaks_pos"),envir=as.environment(".GlobalEnv"));
	if( length(links_peaks_pos)>0 ){ # if any links exist
		load(file=file.path(as.character(logfile[[1]]),"results","profileList_pos"),envir=as.environment(".GlobalEnv"));	
		###############################################################################################
		keep_peaks<-rep(TRUE,length(profileList_pos[[2]][,1]))
		for(i in 1:length(profileList_pos[[2]][,1])){
			if(profileList_pos[[2]][i,5]!=0){
				if(length(links_peaks_pos[[profileList_pos[[2]][i,5]]][[2]])>0){
					keep_peaks[i]<-FALSE
					links_peaks_pos[[profileList_pos[[2]][i,5]]][[2]]<-list() # = IS entry[[2]] in links_peaks_pos
				}
			}
		}
		profileList_pos[[2]]<<-profileList_pos[[2]][keep_peaks,]
		profileList_pos[[6]]<<-0  # index_agglom
		profileList_pos[[7]]<<-0  # index_prof
		profileList_pos[[8]]<<-0  # parameters
		profileList_pos[[1]]<<-data.frame(TRUE,FALSE,FALSE,FALSE)    # state
		colnames(profileList_pos[[1]])<<-c("peaks?","agglom?","profiling","trends?")
		###############################################################################################
		profileList_pos<<-agglomer(
							profileList_pos,
							dmass=(as.numeric(logfile$parameters$prof_dmz)+1),
							ppm=as.logical(as.character(logfile$parameters$prof_ppm)),
							dret=(as.numeric(logfile$parameters$prof_drt)+10)
						)
		if(logfile$parameters[[91]]=="no"){ 				
			profileList_pos<<-partcluster(
								profileList=profileList_pos,
								dmass=as.numeric(logfile$parameters$prof_dmz),
								ppm=as.logical(as.character(logfile$parameters$prof_ppm)),
								dret=as.numeric(logfile$parameters$prof_drt),
								from=FALSE,
								to=FALSE,
								progbar=logfile$parameters$progressBar,
								plotit=FALSE,
								replicates=FALSE
							)
		}else{ # run a profiling in the replicate groups first
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			measurements<-measurements[measurements[,8]=="TRUE",]
			measurements<-measurements[measurements[,4]=="positive",]
			replicates<-measurements$tag3
			IDs<-measurements$ID
			profileList_pos<<-partcluster(
								profileList=profileList_pos,
								dmass=as.numeric(logfile$parameters$prof_dmz),
								ppm=as.logical(as.character(logfile$parameters$prof_ppm)),
								dret=as.numeric(logfile$parameters$prof_drt),
								from=FALSE,
								to=FALSE,
								progbar=logfile$parameters$progressBar,
								plotit=FALSE,
								replicates=replicates,
								IDs
							)
		}
		profileList_pos<<-enviMass:::in_blind(profileList_pos)
		save(profileList_pos,file=file.path(as.character(logfile[[1]]),"results","profileList_pos"));
		profpeaks_pos<<-enviMass:::profiletopeak(profileList_pos,progbar=logfile$parameters[21])		
		profpeaks_pos<<-profpeaks_pos[order(profpeaks_pos[,13],decreasing=TRUE),];
		save(profpeaks_pos,file=file.path(as.character(logfile[[1]]),"results","profpeaks_pos"));
		save(links_peaks_pos,file=file.path(as.character(logfile[[1]]),"results","links_peaks_pos"));	
		cat(paste("\nIS subtraction: ",round((sum(!keep_peaks)/length(keep_peaks)*100),digits=3),"% of peaks removed\n",sep=""))
		###############################################################################################
		save(profileList_pos,file=file.path(as.character(logfile[[1]]),"results","profileList_pos"))
	}
}
#######################################################################################################

#######################################################################################################
# NEGATIVE IONIZATION #################################################################################
if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_neg")){rm(profileList_neg,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="profileList_neg")){rm(profileList_neg)}
if(any(objects(envir=as.environment(".GlobalEnv"))=="links_peaks_neg")){rm(links_peaks_neg,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="links_peaks_neg")){rm(links_peaks_neg)}
if(
	file.exists(file.path(as.character(logfile[[1]]),"results","profileList_neg")) &
	file.exists(file.path(as.character(logfile[[1]]),"results","links_peaks_neg"))	
){

	load(file=file.path(as.character(logfile[[1]]),"results","links_peaks_neg"),envir=as.environment(".GlobalEnv"));
	if( length(links_peaks_neg)>0 ){ # if any links exist
		load(file=file.path(as.character(logfile[[1]]),"results","profileList_neg"),envir=as.environment(".GlobalEnv"));	
		###############################################################################################
		keep_peaks<-rep(TRUE,length(profileList_neg[[2]][,1]))
		for(i in 1:length(profileList_neg[[2]][,1])){
			if(profileList_neg[[2]][i,5]!=0){
				if(length(links_peaks_neg[[profileList_neg[[2]][i,5]]][[2]])>0){
					keep_peaks[i]<-FALSE
					links_peaks_neg[[profileList_neg[[2]][i,5]]][[2]]<-list() # = IS entry[[2]] in links_peaks_pos
				}
			}
		}
		profileList_neg[[2]]<<-profileList_neg[[2]][keep_peaks,]
		profileList_neg[[6]]<<-0  # index_agglom
		profileList_neg[[7]]<<-0  # index_prof
		profileList_neg[[8]]<<-0  # parameters
		profileList_neg[[1]]<<-data.frame(TRUE,FALSE,FALSE,FALSE)    # state
		colnames(profileList_neg[[1]])<<-c("peaks?","agglom?","profiling","trends?")
		###############################################################################################
		profileList_neg<<-agglomer(
							profileList_neg,
							dmass=(as.numeric(logfile$parameters$prof_dmz)+1),
							ppm=as.logical(as.character(logfile$parameters$prof_ppm)),
							dret=(as.numeric(logfile$parameters$prof_drt)+10)
						)
		if(logfile$parameters[[91]]=="no"){ 				
			profileList_neg<<-partcluster(
								profileList_neg,
								dmass=as.numeric(logfile$parameters$prof_dmz),
								ppm=as.logical(as.character(logfile$parameters$prof_ppm)),
								dret=as.numeric(logfile$parameters$prof_drt),
								from=FALSE,
								to=FALSE,
								progbar=logfile$parameters$progressBar,
								plotit=FALSE
							)
		}else{ # run a profiling in the replicate groups first
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			measurements<-measurements[measurements[,8]=="TRUE",]
			measurements<-measurements[measurements[,4]=="negative",]
			replicates<-measurements$tag3
			IDs<-measurements$ID
			profileList_neg<<-partcluster(
								profileList=profileList_neg,
								dmass=as.numeric(logfile$parameters$prof_dmz),
								ppm=as.logical(as.character(logfile$parameters$prof_ppm)),
								dret=as.numeric(logfile$parameters$prof_drt),
								from=FALSE,
								to=FALSE,
								progbar=logfile$parameters$progressBar,
								plotit=FALSE,
								replicates=replicates,
								IDs
							)
		}
		profileList_neg<<-enviMass:::in_blind(profileList_neg)
		save(profileList_neg,file=file.path(as.character(logfile[[1]]),"results","profileList_neg"));
		profpeaks_neg<<-enviMass:::profiletopeak(profileList_neg,progbar=logfile$parameters[21])		
		profpeaks_neg<<-profpeaks_neg[order(profpeaks_neg[,13],decreasing=TRUE),];
		save(profpeaks_neg,file=file.path(as.character(logfile[[1]]),"results","profpeaks_neg"));
		save(links_peaks_neg,file=file.path(as.character(logfile[[1]]),"results","links_peaks_neg"));	
		cat(paste("\nIS subtraction: ",round((sum(!keep_peaks)/length(keep_peaks)*100),digits=3),"% of peaks removed\n",sep=""))
		###############################################################################################
		save(profileList_neg,file=file.path(as.character(logfile[[1]]),"results","profileList_neg"))
	}
}
#######################################################################################################
