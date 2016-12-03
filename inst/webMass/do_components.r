# Componentization ##################################################################


##################################################################################
# POSITIVE #######################################################################
if( 
	file.exists(file.path(logfile[[1]],"results","profileList_pos")) 
){

	##############################################################################	
	if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_pos")){rm(profileList_pos,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="profileList_pos")){rm(profileList_pos)}	
	if(any(objects(envir=as.environment(".GlobalEnv"))=="links_peaks_pos")){rm(links_peaks_pos,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="links_peaks_pos")){rm(links_peaks_pos)}				
	if(any(objects(envir=as.environment(".GlobalEnv"))=="links_profiles_pos")){rm(links_profiles_pos,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="links_profiles_pos")){rm(links_profiles_pos)}					
	load(file.path(as.character(logfile[[1]]),"results","profileList_pos"),envir=as.environment(".GlobalEnv"));
	# links_peaks_pos<-list(); # each entry with 6 lists itself: targets, IS, EIC_correl, isotop, adducts, homol
	load(file.path(as.character(logfile[[1]]),"results","links_peaks_pos"));	
	links_profiles_pos<-list(); # each entry with 6 lists itself: targets, IS, EIC_correl, isotop, adducts, homol
	measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	peaks<-profileList_pos[["peaks"]][,c("sampleIDs","peakIDs","profileIDs","RT")] # to retrieve relations with, sampleID, peakID, profileID, RT
	ord<-order(peaks[,"sampleIDs"],peaks[,"peakIDs"],peaks[,"profileIDs"],decreasing=FALSE)
	peaks<-peaks[ord,]
	use_entries_peaks<-enviMass:::find_empty(links_peaks_pos)
	use_entries_profiles<-enviMass:::find_empty(links_profiles_pos)
	profileList_pos[["index_prof"]][,"links"]<-0
	##############################################################################	

	##############################################################################	
	# annotate target & IS screening matches stored in links_peaks_pos ###########
	if(	
		(
			(logfile$workflow[names(logfile$workflow)=="subtr"]=="no") |	
			(
				((logfile$workflow[names(logfile$workflow)=="IS_screen"]=="yes") &  (logfile$parameters$subtr_IS!="yes")) |
				((logfile$workflow[names(logfile$workflow)=="target_screen"]=="yes") & (logfile$parameters$subtr_target!="yes"))	
			) 
		) & 
		(length(links_peaks_pos)>0) # anything screened?
	){
		cat("\n  Annotation of screening results to profiles")
		for(i in 1:length(profileList_pos[[7]][,1])){
			if(
				any(profileList_pos[[2]][
					(profileList_pos[[7]][i,"start_ID"]:profileList_pos[[7]][i,"end_ID"]),"links"
				]!=0)			
			){
			
				###################################################################
				# add a new link to the profile ###################################
				if( profileList_pos[[7]][i,"links"]==0 ){ 	# establish a new link ...
					if(length(use_entries_profiles)>0){
						at_entry<-use_entries_profiles[1]
						use_entries<-use_entries_profiles[-1]
					}else{
						at_entry<-(length(links_profiles_pos)+1)
					}
					links_profiles_pos[[at_entry]]<-list()
					links_profiles_pos[[at_entry]][[1]]<-list() # target
					links_profiles_pos[[at_entry]][[2]]<-list()	# IS
					links_profiles_pos[[at_entry]][[3]]<-list()	# EIC_correl
						links_profiles_pos[[at_entry]][[3]]<-matrix(ncol=4,nrow=0)
						colnames(links_profiles_pos[[at_entry]][[3]])<-c("linked profile","link counts","ref","no-link counts")
					links_profiles_pos[[at_entry]][[4]]<-list()	# isotop
						links_profiles_pos[[at_entry]][[4]]<-matrix(ncol=3,nrow=0)
						colnames(links_profiles_pos[[at_entry]][[4]])<-c("linked profile","link counts","ref")					
					links_profiles_pos[[at_entry]][[5]]<-list()	# adducts
						links_profiles_pos[[at_entry]][[5]]<-matrix(ncol=3,nrow=0)
						colnames(links_profiles_pos[[at_entry]][[5]])<-c("linked profile","link counts","ref")					
					links_profiles_pos[[at_entry]][[6]]<-list()	# homol		
						links_profiles_pos[[at_entry]][[6]]<-matrix(ncol=3,nrow=0)
						colnames(links_profiles_pos[[at_entry]][[6]])<-c("linked profile","link counts","ref")
					links_profiles_pos[[at_entry]][[7]]<-list()	# group		
					names(links_profiles_pos[[at_entry]])<-c("targ","IS","EIC","isot","adduc","homol","group")
					names(links_profiles_pos)[at_entry]<-as.character(i)
					profileList_pos[[7]][i,"links"]<-at_entry						
				}else{
					at_entry<-profileList_pos[[7]][i,"links"]
				}				
				###################################################################
				# search peaks & their links to compounds #########################
				for(j in (profileList_pos[[7]][i,"start_ID"]:profileList_pos[[7]][i,"end_ID"])){
					if(profileList_pos[[2]][j,"links"]!=0){
						# add IS link #############################################
						if(	length(links_peaks_pos[[profileList_pos[[2]][j,"links"]]][[2]])>0 ){
							for(k in 1:length(links_peaks_pos[[profileList_pos[[2]][j,"links"]]][[2]])){ # if several compound matches exist for this peak
								if( length(links_profiles_pos[[at_entry]][[2]])==0 ){ # make a new entry for profile link to IS
									links_profiles_pos[[at_entry]][[2]]<-
										data.frame(
											links_peaks_pos[[profileList_pos[[2]][j,"links"]]][[2]][[k]],1,stringsAsFactors = FALSE
										)
									names(links_profiles_pos[[at_entry]][[2]])<-c("Compound","Counts")
								}else{
									at<-which(links_profiles_pos[[at_entry]][[2]][,1]==links_peaks_pos[[profileList_pos[[2]][j,"links"]]][[2]][[k]])
									if(length(at)>0){ 	# increment existing link ...
										links_profiles_pos[[at_entry]][[2]][at,2]<-(links_profiles_pos[[at_entry]][[2]][at,2]+1)
									}else{	# or just add a new one?
										links_profiles_pos[[at_entry]][[2]]<-data.frame(
											c(links_profiles_pos[[at_entry]][[2]][,1], links_peaks_pos[[profileList_pos[[2]][j,"links"]]][[2]][[k]]),
											c(links_profiles_pos[[at_entry]][[2]][,2],1),
											stringsAsFactors = FALSE
										)									
										names(links_profiles_pos[[at_entry]][[2]])<-c("Compound","Counts")
									}
								}
							}	
						}					
						# add target link #########################################
						if(	length(links_peaks_pos[[profileList_pos[[2]][j,"links"]]][[1]])>0 ){
							for(k in 1:length(links_peaks_pos[[profileList_pos[[2]][j,"links"]]][[1]])){ # if several compound matches exist for this peak
								if( length(links_profiles_pos[[at_entry]][[1]])==0 ){ # make a new entry for profile link to IS
									links_profiles_pos[[at_entry]][[1]]<-
										data.frame(
											links_peaks_pos[[profileList_pos[[2]][j,"links"]]][[1]][[k]],1,stringsAsFactors = FALSE
										)
									names(links_profiles_pos[[at_entry]][[1]])<-c("Compound","Counts")
								}else{
									at<-which(links_profiles_pos[[at_entry]][[1]][,1]==links_peaks_pos[[profileList_pos[[2]][j,"links"]]][[1]][[k]])
									if(length(at)>0){ 	# increment existing link ...
										links_profiles_pos[[at_entry]][[1]][at,2]<-(links_profiles_pos[[at_entry]][[1]][at,2]+1)
									}else{	# or just add a new one?
										links_profiles_pos[[at_entry]][[1]]<-data.frame(
											c(links_profiles_pos[[at_entry]][[1]][,1], links_peaks_pos[[profileList_pos[[2]][j,"links"]]][[1]][[k]]),
											c(links_profiles_pos[[at_entry]][[1]][,2],1),
											stringsAsFactors = FALSE
										)									
										names(links_profiles_pos[[at_entry]][[1]])<-c("Compound","Counts")
									}
								}
							}					
							
						}
						###########################################################
					}
				}
				###################################################################
				
			}
		}
	}
	rm(links_peaks_pos)
	##############################################################################	
	
	##############################################################################	
	# for each profile, filter all relations to other profiles ###################
			
	##############################################################################	
	
	##############################################################################
	# save! ######################################################################
	save(profileList_pos,file=file.path(as.character(logfile[[1]]),"results","profileList_pos"));
	save(links_profiles_pos,file=file.path(as.character(logfile[[1]]),"results","links_profiles_pos"));	
	##############################################################################	
	
}






##################################################################################
# Negative #######################################################################
if( 
	file.exists(file.path(logfile[[1]],"results","profileList_neg")) 
){

	##############################################################################	
	if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_neg")){rm(profileList_neg,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="profileList_neg")){rm(profileList_neg)}	
	if(any(objects(envir=as.environment(".GlobalEnv"))=="links_peaks_neg")){rm(links_peaks_neg,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="links_peaks_neg")){rm(links_peaks_neg)}				
	if(any(objects(envir=as.environment(".GlobalEnv"))=="links_profiles_neg")){rm(links_profiles_neg,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="links_profiles_neg")){rm(links_profiles_neg)}					
	load(file.path(as.character(logfile[[1]]),"results","profileList_neg"),envir=as.environment(".GlobalEnv"));
	# links_peaks_neg<-list(); # each entry with 6 lists itself: targets, IS, EIC_correl, isotop, adducts, homol
	load(file.path(as.character(logfile[[1]]),"results","links_peaks_neg"));	
	links_profiles_neg<-list(); # each entry with 6 lists itself: targets, IS, EIC_correl, isotop, adducts, homol
	measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	peaks<-profileList_neg[["peaks"]][,c("sampleIDs","peakIDs","profileIDs","RT")] # to retrieve relations with, sampleID, peakID, profileID, RT
	ord<-order(peaks[,"sampleIDs"],peaks[,"peakIDs"],peaks[,"profileIDs"],decreasing=FALSE)
	peaks<-peaks[ord,]
	use_entries_peaks<-enviMass:::find_empty(links_peaks_neg)
	use_entries_profiles<-enviMass:::find_empty(links_profiles_neg)
	profileList_neg[["index_prof"]][,"links"]<-0
	##############################################################################	

	##############################################################################	
	# annotate target & IS screening matches stored in links_peaks_neg ###########
	if(	
		(
			(logfile$workflow[names(logfile$workflow)=="subtr"]=="no") |	
			(
				((logfile$workflow[names(logfile$workflow)=="IS_screen"]=="yes") &  (logfile$parameters$subtr_IS!="yes")) |
				((logfile$workflow[names(logfile$workflow)=="target_screen"]=="yes") & (logfile$parameters$subtr_target!="yes"))	
			) 
		) & 
		(length(links_peaks_neg)>0) # anything screened?
	){
		cat("\n  Annotation of screening results to profiles")
		for(i in 1:length(profileList_neg[[7]][,1])){
			if(
				any(profileList_neg[[2]][
					(profileList_neg[[7]][i,"start_ID"]:profileList_neg[[7]][i,"end_ID"]),"links"
				]!=0)			
			){
			
				###################################################################
				# add a new link to the profile ###################################
				if( profileList_neg[[7]][i,"links"]==0 ){ 	# establish a new link ...
					if(length(use_entries_profiles)>0){
						at_entry<-use_entries_profiles[1]
						use_entries<-use_entries_profiles[-1]
					}else{
						at_entry<-(length(links_profiles_neg)+1)
					}
					links_profiles_neg[[at_entry]]<-list()
					links_profiles_neg[[at_entry]][[1]]<-list() # target
					links_profiles_neg[[at_entry]][[2]]<-list()	# IS
					links_profiles_neg[[at_entry]][[3]]<-list()	# EIC_correl
						links_profiles_neg[[at_entry]][[3]]<-matrix(ncol=4,nrow=0)
						colnames(links_profiles_neg[[at_entry]][[3]])<-c("linked profile","link counts","ref","no-link counts")
					links_profiles_neg[[at_entry]][[4]]<-list()	# isotop
						links_profiles_neg[[at_entry]][[4]]<-matrix(ncol=3,nrow=0)
						colnames(links_profiles_neg[[at_entry]][[4]])<-c("linked profile","link counts","ref")					
					links_profiles_neg[[at_entry]][[5]]<-list()	# adducts
						links_profiles_neg[[at_entry]][[5]]<-matrix(ncol=3,nrow=0)
						colnames(links_profiles_neg[[at_entry]][[5]])<-c("linked profile","link counts","ref")					
					links_profiles_neg[[at_entry]][[6]]<-list()	# homol		
						links_profiles_neg[[at_entry]][[6]]<-matrix(ncol=3,nrow=0)
						colnames(links_profiles_neg[[at_entry]][[6]])<-c("linked profile","link counts","ref")
					links_profiles_neg[[at_entry]][[7]]<-list()	# group		
					names(links_profiles_neg[[at_entry]])<-c("targ","IS","EIC","isot","adduc","homol","group")
					names(links_profiles_neg)[at_entry]<-as.character(i)
					profileList_neg[[7]][i,"links"]<-at_entry						
				}else{
					at_entry<-profileList_neg[[7]][i,"links"]
				}				
				###################################################################
				# search peaks & their links to compounds #########################
				for(j in (profileList_neg[[7]][i,"start_ID"]:profileList_neg[[7]][i,"end_ID"])){
					if(profileList_neg[[2]][j,"links"]!=0){
						# add IS link #############################################
						if(	length(links_peaks_neg[[profileList_neg[[2]][j,"links"]]][[2]])>0 ){
							for(k in 1:length(links_peaks_neg[[profileList_neg[[2]][j,"links"]]][[2]])){ # if several compound matches exist for this peak
								if( length(links_profiles_neg[[at_entry]][[2]])==0 ){ # make a new entry for profile link to IS
									links_profiles_neg[[at_entry]][[2]]<-
										data.frame(
											links_peaks_neg[[profileList_neg[[2]][j,"links"]]][[2]][[k]],1,stringsAsFactors = FALSE
										)
									names(links_profiles_neg[[at_entry]][[2]])<-c("Compound","Counts")
								}else{
									at<-which(links_profiles_neg[[at_entry]][[2]][,1]==links_peaks_neg[[profileList_neg[[2]][j,"links"]]][[2]][[k]])
									if(length(at)>0){ 	# increment existing link ...
										links_profiles_neg[[at_entry]][[2]][at,2]<-(links_profiles_neg[[at_entry]][[2]][at,2]+1)
									}else{	# or just add a new one?
										links_profiles_neg[[at_entry]][[2]]<-data.frame(
											c(links_profiles_neg[[at_entry]][[2]][,1], links_peaks_neg[[profileList_neg[[2]][j,"links"]]][[2]][[k]]),
											c(links_profiles_neg[[at_entry]][[2]][,2],1),
											stringsAsFactors = FALSE
										)									
										names(links_profiles_neg[[at_entry]][[2]])<-c("Compound","Counts")
									}
								}
							}	
						}					
						# add target link #########################################
						if(	length(links_peaks_neg[[profileList_neg[[2]][j,"links"]]][[1]])>0 ){
							for(k in 1:length(links_peaks_neg[[profileList_neg[[2]][j,"links"]]][[1]])){ # if several compound matches exist for this peak
								if( length(links_profiles_neg[[at_entry]][[1]])==0 ){ # make a new entry for profile link to IS
									links_profiles_neg[[at_entry]][[1]]<-
										data.frame(
											links_peaks_neg[[profileList_neg[[2]][j,"links"]]][[1]][[k]],1,stringsAsFactors = FALSE
										)
									names(links_profiles_neg[[at_entry]][[1]])<-c("Compound","Counts")
								}else{
									at<-which(links_profiles_neg[[at_entry]][[1]][,1]==links_peaks_neg[[profileList_neg[[2]][j,"links"]]][[1]][[k]])
									if(length(at)>0){ 	# increment existing link ...
										links_profiles_neg[[at_entry]][[1]][at,2]<-(links_profiles_neg[[at_entry]][[1]][at,2]+1)
									}else{	# or just add a new one?
										links_profiles_neg[[at_entry]][[1]]<-data.frame(
											c(links_profiles_neg[[at_entry]][[1]][,1], links_peaks_neg[[profileList_neg[[2]][j,"links"]]][[1]][[k]]),
											c(links_profiles_neg[[at_entry]][[1]][,2],1),
											stringsAsFactors = FALSE
										)									
										names(links_profiles_neg[[at_entry]][[1]])<-c("Compound","Counts")
									}
								}
							}					
							
						}
						###########################################################
					}
				}
				###################################################################
				
			}
		}
	}
	rm(links_peaks_neg)
	##############################################################################	
	
	##############################################################################	
	# for each profile, filter all relations to other profiles ###################
			
	##############################################################################	
	
	##############################################################################
	# save! ######################################################################
	save(profileList_neg,file=file.path(as.character(logfile[[1]]),"results","profileList_neg"));
	save(links_profiles_neg,file=file.path(as.character(logfile[[1]]),"results","links_profiles_neg"));	
	##############################################################################	
	
}

