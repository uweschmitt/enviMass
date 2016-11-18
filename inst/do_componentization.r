# Profile componentization #######################################################

# -> link upstream to do_profiling.r and screening and all grouping and trend detection
# -> check: profiling should have no downstream dependency on grouping, yes? (seemingly, not)
# -> check subtract: ... if still runs with subtr.r enabled
# -> make adduct node dependent on EIC_correlation (done)
# -> remove componentization results during file removal (done)
# -> run groupings only if componentization is enabled
# -> EIC impacts blind subtraction? depends matrix ok?


# new parameters:
comp_EIC_freq<-.9
comp_num_cor<-3
comp_prof_cor<-.9
comp_max_delRT<-30
comp_method<-2

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
	##############################################################################	

	##############################################################################	
	# for each profile, filter all relations to other profiles ###################
	peaks<-profileList_pos[["peaks"]][,c("sampleIDs","peakIDs","profileIDs","RT")] # to retrieve relations with, sampleID, peakID, profileID, RT
	ord<-order(peaks[,"sampleIDs"],peaks[,"peakIDs"],peaks[,"profileIDs"],decreasing=FALSE)
	peaks<-peaks[ord,]
	use_entries_peaks<-enviMass:::find_empty(links_peaks_pos)
	use_entries_profiles<-enviMass:::find_empty(links_profiles_pos)
	profileList_pos[["index_prof"]][,"links"]<-0
		
if(FALSE){		
	find_mass<-242.0708
	find_mass<-244.068
	find_mass<-246.0649
	find_RT<-(2.45*60)	
	
	finds<-which(
		(profileList_pos[["peaks"]][,"m/z"]>=(find_mass-.001)) &
		(profileList_pos[["peaks"]][,"m/z"]<=(find_mass+.001)) &
		(profileList_pos[["peaks"]][,"RT"]>=(find_RT-10)) &
		(profileList_pos[["peaks"]][,"RT"]<=(find_RT+10)) 
	)	
		
	profileList_pos[["peaks"]][finds,]	
		
	at_ID<-3078
	profileList_pos[["index_prof"]][at_ID,"profile_ID"]
	at_link<-profileList_pos[["index_prof"]][at_ID,"links"]
	links_profiles_pos[[at_link]]	
}		
		
	# insert EIC links ###########################################################
	forIDs<-profileList_pos[["sampleID"]]
	for_files<-list.files(file.path(logfile[[1]],"results","componentization","EIC_corr"))
	keep<-match(forIDs,for_files) # which files are available?
	TRUE_IDs<-(measurements$ID[measurements$adducts=="TRUE"]) # for files which have run through that step
	keep2<-match(forIDs,for_files)
	forIDs<-forIDs[!is.na(keep) & !is.na(keep2)]
	not_found1<-0
	inserted1<-0
	if(length(forIDs)>0){
		cat("\n Retrieving EIC correlation links ")
		for(i in 1:length(forIDs)){
			cat(".")
			load(file=file.path(logfile[[1]],"results","componentization","EIC_corr",forIDs[i]))
			
# any((EIC_pairs[,1]==10487) & (EIC_pairs[,2]==10948)) # 7
# any((EIC_pairs[,1]==10577) & (EIC_pairs[,2]==11047)) # 5		
# any((EIC_pairs[,1]==10648) & (EIC_pairs[,2]==11090)) # 6
	
	
			EIC_pairs<-EIC_pairs[EIC_pairs[,4]>=logfile$parameters$external$EICor_mincor,,drop=FALSE]		
			if(length(EIC_pairs[,1])==0){next}
			# find profiles for first peak
			get1<-cbind(
				rep(as.numeric(forIDs[i]),length(EIC_pairs[,1])),EIC_pairs[,1] # first column sorted correctly
			)
			found1<-enviMass:::rows_compare(get1,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=FALSE,column_order_b=FALSE,get_index=TRUE)
			# find profiles for second peak
			get2<-cbind(
				rep(as.numeric(forIDs[i]),length(EIC_pairs[,2])),EIC_pairs[,2] # second column requires sorting
			)
			found2<-enviMass:::rows_compare(get2,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=TRUE,column_order_b=FALSE,get_index=TRUE)
			for(j in 1:length(found1)){ # insert links
				# insert PROFILE LINKS #############################################
				if(found1[j]==0){not_found1<-(not_found1+1);next} # e.g., peak blind-removed 
				if(found2[j]==0){not_found1<-(not_found1+1);next}		
				inserted1<-(inserted1+1);				
				# (1) insert link to second profile for the first profile
				prof1<-peaks[found1[j],"profileIDs"][[1]]
				prof2<-peaks[found2[j],"profileIDs"][[1]]
				# enable check ... pairs musst have very similar retention times
				if(FALSE){
					cat("\n");
					cat(peaks[found1[j],"RT"]);cat(" - ")
					cat(peaks[found2[j],"RT"])
				}
				if(profileList_pos[[7]][prof1,"profile_ID"]!=prof1){stop("\nComponentization: debug me, #1!")}				
				if(profileList_pos[[7]][prof1,"links"]==0){ 	# establish a new link ...
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
					names(links_profiles_pos)[at_entry]<-as.character(prof1)
					profileList_pos[[7]][prof1,"links"]<-at_entry						
				}else{
					at_entry<-profileList_pos[[7]][prof1,"links"]
				}
				here<-which(links_profiles_pos[[at_entry]][[3]][,"linked profile"]==prof2)
				if(length(here)==0){			
					links_profiles_pos[[at_entry]][[3]]<-rbind(
						links_profiles_pos[[at_entry]][[3]], c(prof2,1,0,0)
					)
				}else{
					links_profiles_pos[[at_entry]][[3]][here,"link counts"]<-(links_profiles_pos[[at_entry]][[3]][here,"link counts"]+1)
				}
				# (2) insert link to first profile for the second profile
				if(profileList_pos[[7]][prof2,"profile_ID"]!=prof2){stop("\nComponentization: debug me, #1!")}				
				if(profileList_pos[[7]][prof2,"links"]==0){ 	# establish a new link ...
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
					names(links_profiles_pos)[at_entry]<-as.character(prof2)
					profileList_pos[[7]][prof2,"links"]<-at_entry						
				}else{
					at_entry<-profileList_pos[[7]][prof2,"links"]
				}
				here<-which(links_profiles_pos[[at_entry]][[3]][,"linked profile"]==prof1)
				if(length(here)==0){
					links_profiles_pos[[at_entry]][[3]]<-rbind(
						links_profiles_pos[[at_entry]][[3]], c(prof1,1,0,0)
					)
				}else{
					links_profiles_pos[[at_entry]][[3]][here,"link counts"]<-(links_profiles_pos[[at_entry]][[3]][here,"link counts"]+1)
				}				
				# insert PEAK LINKS ##############################################

			}
		}
		cat(" done.")
	}
	# insert EIC non-links #######################################################
	forIDs<-profileList_pos[["sampleID"]]
	for_files<-list.files(file.path(logfile[[1]],"results","componentization","EIC_corr"))
	keep<-match(forIDs,for_files) # which files are available?
	TRUE_IDs<-(measurements$ID[measurements$adducts=="TRUE"]) # for files which have run through that step
	keep2<-match(forIDs,for_files)
	forIDs<-forIDs[!is.na(keep) & !is.na(keep2)]
	not_found2<-0
	inserted2<-0
	if(length(forIDs)>0){
		cat("\n Retrieving EIC correlation non-links ")
		for(i in 1:length(forIDs)){
			cat(".")
			load(file=file.path(logfile[[1]],"results","componentization","EIC_corr",forIDs[i]))
			EIC_pairs<-EIC_pairs[EIC_pairs[,4]<logfile$parameters$external$EICor_mincor,,drop=FALSE]		
			if(length(EIC_pairs[,1])==0){next}
			# find profiles for first peak
			get1<-cbind(
				rep(as.numeric(forIDs[i]),length(EIC_pairs[,1])),EIC_pairs[,1] # first column sorted correctly
			)
			found1<-enviMass:::rows_compare(get1,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=FALSE,column_order_b=FALSE,get_index=TRUE)
			# find profiles for second peak
			get2<-cbind(
				rep(as.numeric(forIDs[i]),length(EIC_pairs[,2])),EIC_pairs[,2] # second column requires sorting
			)
			found2<-enviMass:::rows_compare(get2,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=TRUE,column_order_b=FALSE,get_index=TRUE)
			for(j in 1:length(found1)){ # insert links
				# insert PROFILE LINKS #############################################
				if(found1[j]==0){not_found2<-(not_found2+1);next} # e.g., peak blind-removed 
				if(found2[j]==0){not_found2<-(not_found2+1);next}			
				# (1) insert link to second profile for the first profile
				prof1<-peaks[found1[j],"profileIDs"][[1]]
				prof2<-peaks[found2[j],"profileIDs"][[1]]
				# profiles have to have positive entries from above
				if(profileList_pos[[7]][prof1,"links"]==0) next
				if(profileList_pos[[7]][prof2,"links"]==0) next
				inserted2<-(inserted2+1);
				# enable check ... pairs musst have very similar retention times
				if(FALSE){
					cat("\n");
					cat(peaks[found1[j],"RT"]);cat(" - ")
					cat(peaks[found2[j],"RT"])
				}
				if(profileList_pos[[7]][prof1,"profile_ID"]!=prof1){stop("\nComponentization: debug me, #1!")}				
				at_entry<-profileList_pos[[7]][prof1,"links"]
				if( # no non-link required if no link inserted before
					!any(links_profiles_pos[[at_entry]]$EIC[,1]==prof2) 
				){ next }	
				here<-which(links_profiles_pos[[at_entry]][[3]][,"linked profile"]==prof2)
				if(length(here)==0){			
					links_profiles_pos[[at_entry]][[3]]<-rbind(
						links_profiles_pos[[at_entry]][[3]], c(prof2,0,0,1)
					)
				}else{
					links_profiles_pos[[at_entry]][[3]][here,"ref"]<-(links_profiles_pos[[at_entry]][[3]][here,"ref"]+1)
				}
				# (2) insert link to first profile for the second profile
				if(profileList_pos[[7]][prof2,4]!=prof2){stop("\nComponentization: debug me, #1!")}				
				at_entry<-profileList_pos[[7]][prof2,"links"] # must exist already - see above no non-link "next"
# GAPPED?
				here<-which(links_profiles_pos[[at_entry]][[3]][,"linked profile"]==prof1)
				if(length(here)==0){
					links_profiles_pos[[at_entry]][[3]]<-rbind(
						links_profiles_pos[[at_entry]][[3]],c(prof1,0,0,1)
					)
				}else{
					links_profiles_pos[[at_entry]][[3]][here,"ref"]<-(links_profiles_pos[[at_entry]][[3]][here,"ref"]+1)
				}				
				# insert PEAK LINKS ##############################################

			}
		}
		cat(" done.")
	}
	# insert isotopologue links ##################################################
	forIDs<-profileList_pos[["sampleID"]]
	for_files<-list.files(file.path(logfile[[1]],"results","componentization","isotopologues"))
	keep<-match(forIDs,for_files) # which files are available?
	TRUE_IDs<-(measurements$ID[measurements$adducts=="TRUE"]) # for files which have run through that step
	keep2<-match(forIDs,for_files)
	forIDs<-forIDs[!is.na(keep) & !is.na(keep2)]
	not_found3<-0
	inserted3<-0
	if(length(forIDs)>0){
		cat("\n Retrieving isotopologue links ")
		for(i in 1:length(forIDs)){
			cat(".")
			load(file=file.path(logfile[[1]],"results","componentization","isotopologues",forIDs[i]))
			if(length(Isot_pairs[,1])==0){next}
			# find profiles for first peak
			get1<-cbind(
				rep(as.numeric(forIDs[i]),length(Isot_pairs[,1])),Isot_pairs[,1]
			)
			found1<-enviMass:::rows_compare(get1,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=FALSE,column_order_b=FALSE,get_index=TRUE)
			# find profiles for second peak
			get2<-cbind(
				rep(as.numeric(forIDs[i]),length(Isot_pairs[,2])),Isot_pairs[,2]
			)
			found2<-enviMass:::rows_compare(get2,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=TRUE,column_order_b=FALSE,get_index=TRUE)
			for(j in 1:length(found1)){ # insert links
				# insert PROFILE LINKS #############################################
				if(found1[j]==0){not_found3<-(not_found3+1);next} # e.g., peak blind-removed 
				if(found2[j]==0){not_found3<-(not_found3+1);next}			
				inserted3<-(inserted3+1);				
				# enable check ... pairs musst have very similar retention times
				if(FALSE){
					cat("\n");
					cat(peaks[found1[j],"RT"]);cat(" - ")
					cat(peaks[found2[j],"RT"])				
				}
				# (1) insert link to second profile for the first profile
				prof1<-peaks[found1[j],"profileIDs"][[1]]
				prof2<-peaks[found2[j],"profileIDs"][[1]]				
				if(profileList_pos[[7]][prof1,"profile_ID"]!=prof1){stop("\nComponentization: debug me, #1!")}				
				if(profileList_pos[[7]][prof1,"links"]==0){ 	# establish a new link ...
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
					names(links_profiles_pos)[at_entry]<-as.character(prof1)
					profileList_pos[[7]][prof1,"links"]<-at_entry						
				}else{
					at_entry<-profileList_pos[[7]][prof1,"links"]
				}
				here<-which(links_profiles_pos[[at_entry]][[4]][,"linked profile"]==prof2)
				if(length(here)==0){
					links_profiles_pos[[at_entry]][[4]]<-rbind(
						links_profiles_pos[[at_entry]][[4]], c(prof2,1,0)
					)
				}else{
					links_profiles_pos[[at_entry]][[4]][here,"link counts"]<-(links_profiles_pos[[at_entry]][[4]][here,"link counts"]+1)
				}
				# (2) insert link to first profile for the second profile
				if(profileList_pos[[7]][prof2,"profile_ID"]!=prof2){stop("\nComponentization: debug me, #1!")}				
				if(profileList_pos[[7]][prof2,"links"]==0){ 	# establish a new link ...
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
					names(links_profiles_pos)[at_entry]<-as.character(prof2)
					profileList_pos[[7]][prof2,"links"]<-at_entry						
				}else{
					at_entry<-profileList_pos[[7]][prof2,"links"]
				}
				here<-which(links_profiles_pos[[at_entry]][[4]][,"linked profile"]==prof1)
				if(length(here)==0){
					links_profiles_pos[[at_entry]][[4]]<-rbind(
						links_profiles_pos[[at_entry]][[4]], c(prof1,1,0)
					)
				}else{
					links_profiles_pos[[at_entry]][[4]][here,"link counts"]<-(links_profiles_pos[[at_entry]][[4]][here,"link counts"]+1)
				}				
				# insert PEAK LINKS ##############################################
				
			}
		}
		cat(" done.")
	}
	# insert adduct links ########################################################
	forIDs<-profileList_pos[["sampleID"]]
	for_files<-list.files(file.path(logfile[[1]],"results","componentization","adducts"))
	keep<-match(forIDs,for_files) # which files are available?
	TRUE_IDs<-(measurements$ID[measurements$adducts=="TRUE"]) # for files which have run through that step
	keep2<-match(forIDs,for_files)
	forIDs<-forIDs[!is.na(keep) & !is.na(keep2)]
	not_found4<-0
	inserted4<-0
	if(length(forIDs)>0){
		cat("\n Retrieving adduct links ")
		for(i in 1:length(forIDs)){
			cat(".")
			load(file=file.path(logfile[[1]],"results","componentization","adducts",forIDs[i]))
			if(length(Adduct_pairs[,1])==0){next}
			# find profiles for first peak
			get1<-cbind(
				rep(as.numeric(forIDs[i]),length(Adduct_pairs[,1])),Adduct_pairs[,1]
			)
			found1<-enviMass:::rows_compare(get1,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=FALSE,column_order_b=FALSE,get_index=TRUE)
			# find profiles for second peak
			get2<-cbind(
				rep(as.numeric(forIDs[i]),length(Adduct_pairs[,2])),Adduct_pairs[,2]
			)
			found2<-enviMass:::rows_compare(get2,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=TRUE,column_order_b=FALSE,get_index=TRUE)
			for(j in 1:length(found1)){ # insert links	
				if(found1[j]==0){not_found4<-(not_found4+1);next} # e.g., peak blind-removed 
				if(found2[j]==0){not_found4<-(not_found4+1);next}		
				inserted4<-(inserted4+1);
				# enable check ... pairs musst have very similar retention times
				if(FALSE){
					cat("\n");
					cat(peaks[found1[j],"RT"]);cat(" - ")
					cat(peaks[found2[j],"RT"])				
				}
				# insert PROFILE LINKS #############################################	
				# (1) insert link to second profile for the first profile
				prof1<-peaks[found1[j],"profileIDs"][[1]]
				prof2<-peaks[found2[j],"profileIDs"][[1]]
				if(profileList_pos[[7]][prof1,"profile_ID"]!=prof1){stop("\nComponentization: debug me, #1!")}				
				if(profileList_pos[[7]][prof1,"links"]==0){ 	# establish a new link ...
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
					names(links_profiles_pos)[at_entry]<-as.character(prof1)
					profileList_pos[[7]][prof1,"links"]<-at_entry						
				}else{
					at_entry<-profileList_pos[[7]][prof1,"links"]
				}
				here<-which(links_profiles_pos[[at_entry]][[5]][,"linked profile"]==prof2)
				if(length(here)==0){
					links_profiles_pos[[at_entry]][[5]]<-rbind(
						links_profiles_pos[[at_entry]][[5]], c(prof2,1,0)
					)
				}else{
					links_profiles_pos[[at_entry]][[5]][here,"link counts"]<-(links_profiles_pos[[at_entry]][[5]][here,"link counts"]+1)
				}
				# (2) insert link to first profile for the second profile
				if(profileList_pos[[7]][prof2,"profile_ID"]!=prof2){stop("\nComponentization: debug me, #1!")}				
				if(profileList_pos[[7]][prof2,"links"]==0){ 	# establish a new link ...
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
					names(links_profiles_pos)[at_entry]<-as.character(prof2)
					profileList_pos[[7]][prof2,"links"]<-at_entry						
				}else{
					at_entry<-profileList_pos[[7]][prof2,"links"]
				}
				here<-which(links_profiles_pos[[at_entry]][[5]][,"linked profile"]==prof1)
				if(length(here)==0){
					links_profiles_pos[[at_entry]][[5]]<-rbind(
						links_profiles_pos[[at_entry]][[5]], c(prof1,1,0)
					)
				}else{
					links_profiles_pos[[at_entry]][[5]][here,"link counts"]<-(links_profiles_pos[[at_entry]][[5]][here,"link counts"]+1)
				}				
				# insert PEAK LINKS ##############################################
	
			}
		}
		cat(" done.")
	}	
	# insert homologen links #####################################################
	forIDs<-profileList_pos[["sampleID"]]
	for_files<-list.files(file.path(logfile[[1]],"results","componentization","homologues"))
	keep<-match(forIDs,for_files) # which files are available?
	TRUE_IDs<-(measurements$ID[measurements$homologues=="TRUE"]) # for files which have run through that step
	keep2<-match(forIDs,for_files)
	forIDs<-forIDs[!is.na(keep) & !is.na(keep2)]
	not_found5<-0
	inserted5<-0
	if(length(forIDs)>0){
		cat("\n Retrieving homologue links ")
		for(i in 1:length(forIDs)){
			cat(".")
			load(file=file.path(logfile[[1]],"results","componentization","homologues",forIDs[i]))	
			if(length(Homol_groups[,1])==0){next}	
			keep_group<-rep(TRUE,max(Homol_groups[,3])) # make sure all peaks in a group are still present
			# find profiles for first peak
			get1<-cbind(
				rep(as.numeric(forIDs[i]),length(Homol_groups[,1])),Homol_groups[,1]
			)
			found1<-enviMass:::rows_compare(get1,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=FALSE,column_order_b=FALSE,get_index=TRUE)
			keep_group[Homol_groups[found1==0,3]]<-FALSE
			if(!any(keep_group)){next} # none left with all peaks still embedded into profiles
			# find profiles for second peak
			get2<-cbind(
				rep(as.numeric(forIDs[i]),length(Homol_groups[,2])),Homol_groups[,2]
			)
			found2<-enviMass:::rows_compare(get2,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=TRUE,column_order_b=FALSE,get_index=TRUE)
			if(!any(keep_group)){next} # none left with all peaks still embedded into profiles
			for(j in 1:length(found1)){ # insert links	
				if(found1[j]==0){not_found5<-(not_found5+1);next} # e.g., peak blind-removed 
				if(found2[j]==0){not_found5<-(not_found5+1);next}
				if(Homol_groups[j,3]!=Homol_groups[j,3]){stop("\n Debug componentization for its homologue relations!")}
				if(!keep_group[Homol_groups[j,3]]){next} # not all peaks for this group were embedded into profiles, e.g., peak blind-removed
				inserted5<-(inserted5+1);
				# enable check ... pairs musst have very similar retention times
				if(FALSE){
					cat("\n");
					cat(peaks[found1[j],"RT"]);cat(" - ")
					cat(peaks[found2[j],"RT"])				
				}
				# insert PROFILE LINKS #############################################	
				# (1) insert link to second profile for the first profile
				prof1<-peaks[found1[j],"profileIDs"][[1]]
				prof2<-peaks[found2[j],"profileIDs"][[1]]
				if(profileList_pos[[7]][prof1,"profile_ID"]!=prof1){stop("\nComponentization: debug me, #1!")}				
				if(profileList_pos[[7]][prof1,"links"]==0){ 	# establish a new link ...
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
					names(links_profiles_pos)[at_entry]<-as.character(prof1)
					profileList_pos[[7]][prof1,"links"]<-at_entry						
				}else{
					at_entry<-profileList_pos[[7]][prof1,"links"]
				}
				here<-which(links_profiles_pos[[at_entry]][["homol"]][,"linked profile"]==prof2)
				if(length(here)==0){
					links_profiles_pos[[at_entry]][["homol"]]<-rbind(
						links_profiles_pos[[at_entry]][["homol"]], c(prof2,1,0)
					)
				}else{
					links_profiles_pos[[at_entry]][["homol"]][here,"link counts"]<-(links_profiles_pos[[at_entry]][["homol"]][here,"link counts"]+1)
				}
				# (2) insert link to first profile for the second profile
				if(profileList_pos[[7]][prof2,"profile_ID"]!=prof2){stop("\nComponentization: debug me, #1!")}				
				if(profileList_pos[[7]][prof2,"links"]==0){ 	# establish a new link ...
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
					names(links_profiles_pos)[at_entry]<-as.character(prof2)
					profileList_pos[[7]][prof2,"links"]<-at_entry						
				}else{
					at_entry<-profileList_pos[[7]][prof2,"links"]
				}
				here<-which(links_profiles_pos[[at_entry]][["homol"]][,"linked profile"]==prof1)
				if(length(here)==0){
					links_profiles_pos[[at_entry]][["homol"]]<-rbind(
						links_profiles_pos[[at_entry]][["homol"]], c(prof1,1,0)
					)
				}else{
					links_profiles_pos[[at_entry]][["homol"]][here,"link counts"]<-(links_profiles_pos[[at_entry]][["homol"]][here,"link counts"]+1)
				}				
			}	
		}
	}
	##############################################################################	
	# summary output #############################################################
	got_links<-data.frame(
		c(inserted1,inserted2,inserted3,inserted4,inserted5),
		c(not_found1,not_found2,not_found3,not_found4,not_found5),
		row.names=c("EIC","non-EIC","Isotopologues","Adducts","Homologues")
	)
	names(got_links)<-c("Accepted links","Rejected links")
	cat("\n\n");print(got_links);cat("\n\n");
	##############################################################################	
	# filter links based on co-occurences etc. ###################################
	if(TRUE){
		# (1) check EIC relations: 
		# (1.1) weight detected positive against detected negative (= within RT window) cases
		ii<-0;aa<-0
		for(n in 1:length(links_profiles_pos)){
			if(length(links_profiles_pos[[n]]$EIC[,1])>0){
				for(m in 1:length(links_profiles_pos[[n]]$EIC[,1])){
					########################################################
					if(links_profiles_pos[[n]]$EIC[m,3]==0){
						ii<-(ii+1)
						if(
							(links_profiles_pos[[n]]$EIC[m,2]/(links_profiles_pos[[n]]$EIC[m,2]+links_profiles_pos[[n]]$EIC[m,4]))<comp_EIC_freq
						){
							aa<-(aa+1)
							# mark lack of relation in this profile ...
							links_profiles_pos[[n]]$EIC[m,3]<-(-3)
							# ... and the linked profile
							links_profiles_pos[[	
									profileList_pos[[7]][
										links_profiles_pos[[n]]$EIC[m,1]
									,18]
							]]$EIC[
								links_profiles_pos[[	
									profileList_pos[[7]][
										links_profiles_pos[[n]]$EIC[m,1]
									,18]
								]]$EIC[,1]==as.numeric(names(links_profiles_pos)[n])
							,3]<-(-3)			
						}
					}
					########################################################
				}
			}
		}	
		# (1.2) remaining: weight detected positive against all negative (= within & outside RT window) cases
		bb<-0
		for(n in 1:length(links_profiles_pos)){
			if(length(links_profiles_pos[[n]]$EIC[,"linked profile"])>0){
				for(m in 1:length(links_profiles_pos[[n]]$EIC[,"linked profile"])){
					########################################################
					if(links_profiles_pos[[n]]$EIC[m,"ref"]==0){
						prof1<-as.numeric(names(links_profiles_pos)[n])
						prof2<-links_profiles_pos[[n]]$EIC[m,"linked profile"]
						del1<-profileList_pos[[7]][prof1,3][[1]]
						if( (links_profiles_pos[[n]]$EIC[m,2]+links_profiles_pos[[n]]$EIC[m,"no-link counts"])==del1 ){ # no other peaks left
							next
						}						
						del2<-profileList_pos[[7]][prof2,3][[1]]
						if( (links_profiles_pos[[n]]$EIC[m,2]+links_profiles_pos[[n]]$EIC[m,"no-link counts"])==del2 ){ # no other peaks left
							next
						}	
						these<-profileList_pos[[2]][
							profileList_pos[[7]][prof1,1]:profileList_pos[[7]][prof1,2]
						,6]
						those<-profileList_pos[[2]][
							profileList_pos[[7]][prof2,1]:profileList_pos[[7]][prof2,2]
						,6]
						if(del1<del2){
							num_sam<-sum(!is.na(match(these,those)))
						}else{
							num_sam<-sum(!is.na(match(those,these)))
						}	
						if((links_profiles_pos[[n]]$EIC[m,2]+links_profiles_pos[[n]]$EIC[m,"no-link counts"])==num_sam){ # no other peaks left
							next
						}
						if( (links_profiles_pos[[n]]$EIC[m,2]/num_sam)<comp_EIC_freq ){
							bb<-(bb+1)
							# mark lack of relation in this profile prof1 ...
							links_profiles_pos[[n]]$EIC[m,"ref"]<-(-5)
							# ... and the linked profile prof2
							links_profiles_pos[[	
									profileList_pos[[7]][prof2,18]
							]]$EIC[
								links_profiles_pos[[	
									profileList_pos[[7]][prof2,18]
								]]$EIC[,1]==prof1
							,"ref"]<-(-5)						
						}
					}
					########################################################
				}
			}
		}	
		#(aa+bb)/ii
		# (2) apply results to isotopologues & adducts
		dd<-0;ee<-0
		for(n in 1:length(links_profiles_pos)){
			if(length(links_profiles_pos[[n]]$EIC[,1])>0){	
				if(length(links_profiles_pos[[n]]$adduc[,1])>0){
					links_profiles_pos[[n]]$adduc[ # mark adduct links with bad EIC correlation
						!is.na(match(links_profiles_pos[[n]]$adduc[,1],links_profiles_pos[[n]]$EIC[
							links_profiles_pos[[n]]$EIC[,3]<0
						,1]))
					,3]<-(-1)
					dd<-(dd+ sum(links_profiles_pos[[n]]$adduc[,3]<0) )
					links_profiles_pos[[n]]$EIC[ # mark adducts in EIC links
						!is.na(match(links_profiles_pos[[n]]$EIC[,1],links_profiles_pos[[n]]$adduc[,1]))
					,3]<-(-2)					
				}
				if(length(links_profiles_pos[[n]]$isot[,1])>0){
					links_profiles_pos[[n]]$isot[ # mark isotopologue links with bad EIC correlation
						!is.na(match(links_profiles_pos[[n]]$isot[,1],links_profiles_pos[[n]]$EIC[
							links_profiles_pos[[n]]$EIC[,3]<0
						,1]))
					,3]<-(-1)
					ee<-(ee+ sum(links_profiles_pos[[n]]$isot[,3]<0) )
					links_profiles_pos[[n]]$EIC[ # mark isotopologue in EIC links
						!is.na(match(links_profiles_pos[[n]]$EIC[,1],links_profiles_pos[[n]]$isot[,1]))
					,3]<-(-1)					
					
				}				
			}
		}
		# (3) clean & sort remaining relations 
		for(n in 1:length(links_profiles_pos)){
			if(length(links_profiles_pos[[n]]$EIC[,1])>0){
				links_profiles_pos[[n]]$EIC<-links_profiles_pos[[n]]$EIC[links_profiles_pos[[n]]$EIC[,3]>=0,,drop=FALSE]
			}		
			if(length(links_profiles_pos[[n]]$EIC[,1])>0){
				links_profiles_pos[[n]]$EIC<-links_profiles_pos[[n]]$EIC[order(links_profiles_pos[[n]]$EIC[,1],decreasing=FALSE),,drop=FALSE]
				links_profiles_pos[[n]]$EIC[,3]<-0;
			}	
			if(length(links_profiles_pos[[n]]$adduc[,1])>0){
				links_profiles_pos[[n]]$adduc<-links_profiles_pos[[n]]$adduc[links_profiles_pos[[n]]$adduc[,3]>=0,,drop=FALSE]
			}
			if(length(links_profiles_pos[[n]]$adduc[,1])>0){
				links_profiles_pos[[n]]$adduc<-links_profiles_pos[[n]]$adduc[order(links_profiles_pos[[n]]$adduc[,1],decreasing=FALSE),,drop=FALSE]			
				links_profiles_pos[[n]]$adduc[,3]<-0;
			}
			if(length(links_profiles_pos[[n]]$isot[,1])>0){
				links_profiles_pos[[n]]$isot<-links_profiles_pos[[n]]$isot[links_profiles_pos[[n]]$isot[,3]>=0,,drop=FALSE]			
			}	
			if(length(links_profiles_pos[[n]]$isot[,1])>0){
				links_profiles_pos[[n]]$isot<-links_profiles_pos[[n]]$isot[order(links_profiles_pos[[n]]$isot[,1],decreasing=FALSE),,drop=FALSE]
				links_profiles_pos[[n]]$isot[,3]<-0
			}		
			if(length(links_profiles_pos[[n]]$homol[,1])>1){
				links_profiles_pos[[n]]$homol<-links_profiles_pos[[n]]$homol[order(links_profiles_pos[[n]]$homol[,1],decreasing=FALSE),,drop=FALSE]
				links_profiles_pos[[n]]$homol[,3]<-0
			}							
		}
		# (4) correlate remaining relations ...
		ff<-0
		mark_it<-TRUE # set non-correlatable links to 1000
		plot_it<-FALSE
		for(n in 1:length(links_profiles_pos)){
			if(length(links_profiles_pos[[n]]$EIC[,1])>0){
				ff<-(ff+length(links_profiles_pos[[n]]$EIC[,1]))
				for(m in 1:length(links_profiles_pos[[n]]$EIC[,1])){
					if(links_profiles_pos[[n]]$EIC[m,3]>0){next}; # done before via linked profile		
					prof1<-as.numeric(names(links_profiles_pos)[n])
					prof2<-links_profiles_pos[[n]]$EIC[m,1]					
					if(links_profiles_pos[[n]]$EIC[m,2]<comp_num_cor){
						if(mark_it){
							links_profiles_pos[[n]]$EIC[m,3]<-1000
							# ... and the linked profile prof2
							links_profiles_pos[[	
									profileList_pos[[7]][prof2,18]
							]]$EIC[
								links_profiles_pos[[	
									profileList_pos[[7]][prof2,18]
								]]$EIC[,1]==prof1
							,3]<-1000					
						}
						next
					}
					del1<-profileList_pos[[7]][prof1,3][[1]]					
					del2<-profileList_pos[[7]][prof2,3][[1]]	
					these<-profileList_pos[[2]][profileList_pos[[7]][prof1,1]:profileList_pos[[7]][prof1,2],6]
					those<-profileList_pos[[2]][profileList_pos[[7]][prof2,1]:profileList_pos[[7]][prof2,2],6]
					this<-match(these,those)
					int1<-((profileList_pos[[2]][profileList_pos[[7]][prof1,1]:profileList_pos[[7]][prof1,2],2])[!is.na(this)])
					int2<-((profileList_pos[[2]][profileList_pos[[7]][prof2,1]:profileList_pos[[7]][prof2,2],2])[this[!is.na(this)]])					
					correl<-cor(int1,int2)
					if(is.na(correl)){stop("Correl. issue#1")}
					# mark lack of relation in this profile prof1 ...
					links_profiles_pos[[n]]$EIC[m,3]<-correl
					# ... and the linked profile prof2
					links_profiles_pos[[	
							profileList_pos[[7]][prof2,18]
					]]$EIC[
						links_profiles_pos[[	
							profileList_pos[[7]][prof2,18]
						]]$EIC[,1]==prof1
					,3]<-correl		
					if(plot_it){
						delmass<-round(abs((profileList_pos[[7]][prof1,14][[1]])-(profileList_pos[[7]][prof2,14][[1]])),digits=3)
						plot(
							int1,int2,
							pch=19,cex=.8,
							xlab="",ylab="",main=paste("EIC",round(correl,digits=3),delmass,sep=" / "))
						invisible(readline(prompt="Press [enter] to continue"))
					}
				}
			}	
			if(length(links_profiles_pos[[n]]$adduc[,1])>0){
				ff<-(ff+length(links_profiles_pos[[n]]$adduc[,1]))
				for(m in 1:length(links_profiles_pos[[n]]$adduc[,1])){
					if(links_profiles_pos[[n]]$adduc[m,3]>0){next}; # done before via linked profile
					prof1<-as.numeric(names(links_profiles_pos)[n])
					prof2<-links_profiles_pos[[n]]$adduc[m,1]
					if(links_profiles_pos[[n]]$adduc[m,2]<comp_num_cor){						
						if(mark_it){
							links_profiles_pos[[n]]$adduc[m,3]<-1000
							# ... and the linked profile prof2
							links_profiles_pos[[	
									profileList_pos[[7]][prof2,18]
							]]$adduc[
								links_profiles_pos[[	
									profileList_pos[[7]][prof2,18]
								]]$adduc[,1]==prof1
							,3]<-1000								
						}
						next
					}
					del1<-profileList_pos[[7]][prof1,3][[1]]					
					del2<-profileList_pos[[7]][prof2,3][[1]]	
					these<-profileList_pos[[2]][profileList_pos[[7]][prof1,1]:profileList_pos[[7]][prof1,2],6]
					those<-profileList_pos[[2]][profileList_pos[[7]][prof2,1]:profileList_pos[[7]][prof2,2],6]
					this<-match(these,those)
					int1<-((profileList_pos[[2]][profileList_pos[[7]][prof1,1]:profileList_pos[[7]][prof1,2],2])[!is.na(this)])
					int2<-((profileList_pos[[2]][profileList_pos[[7]][prof2,1]:profileList_pos[[7]][prof2,2],2])[this[!is.na(this)]])					
					correl<-cor(int1,int2)
					if(is.na(correl)){stop("Correl. issue#2")}
					# mark lack of relation in this profile prof1 ...
					links_profiles_pos[[n]]$adduc[m,3]<-correl
					# ... and the linked profile prof2
					links_profiles_pos[[	
							profileList_pos[[7]][prof2,18]
					]]$adduc[
						links_profiles_pos[[	
							profileList_pos[[7]][prof2,18]
						]]$adduc[,1]==prof1
					,3]<-correl		
					if(plot_it){
						delmass<-round(abs((profileList_pos[[7]][prof1,14][[1]])-(profileList_pos[[7]][prof2,14][[1]])),digits=3)
						plot(
							int1,int2,
							pch=19,cex=.8,
							xlab="",ylab="",main=paste("Adduc",round(correl,digits=3),delmass,sep=" / "))
						invisible(readline(prompt="Press [enter] to continue"))
					}					
				}
			}
			if(length(links_profiles_pos[[n]]$isot[,1])>0){
				ff<-(ff+length(links_profiles_pos[[n]]$isot[,1]))			
				for(m in 1:length(links_profiles_pos[[n]]$isot[,1])){
					if(links_profiles_pos[[n]]$isot[m,3]>0){next}; # done before via linked profile
					prof1<-as.numeric(names(links_profiles_pos)[n])
					prof2<-links_profiles_pos[[n]]$isot[m,1]
					if(links_profiles_pos[[n]]$isot[m,2]<comp_num_cor){
						if(mark_it){
								links_profiles_pos[[n]]$isot[m,3]<-1000
								# ... and the linked profile prof2
								links_profiles_pos[[	
										profileList_pos[[7]][prof2,18]
								]]$isot[
									links_profiles_pos[[	
										profileList_pos[[7]][prof2,18]
									]]$isot[,1]==prof1
								,3]<-1000
						}
						next
					}
					del1<-profileList_pos[[7]][prof1,3][[1]]					
					del2<-profileList_pos[[7]][prof2,3][[1]]	
					these<-profileList_pos[[2]][profileList_pos[[7]][prof1,1]:profileList_pos[[7]][prof1,2],6]
					those<-profileList_pos[[2]][profileList_pos[[7]][prof2,1]:profileList_pos[[7]][prof2,2],6]
					this<-match(these,those)
					int1<-((profileList_pos[[2]][profileList_pos[[7]][prof1,1]:profileList_pos[[7]][prof1,2],2])[!is.na(this)])
					int2<-((profileList_pos[[2]][profileList_pos[[7]][prof2,1]:profileList_pos[[7]][prof2,2],2])[this[!is.na(this)]])					
					correl<-cor(int1,int2)
					if(is.na(correl)){stop("Correl. issue#3")}
					# mark lack of relation in this profile prof1 ...
					links_profiles_pos[[n]]$isot[m,3]<-correl
					# ... and the linked profile prof2
					links_profiles_pos[[	
							profileList_pos[[7]][prof2,18]
					]]$isot[
						links_profiles_pos[[	
							profileList_pos[[7]][prof2,18]
						]]$isot[,1]==prof1
					,3]<-correl	
					if(plot_it){
						delmass<-round(abs((profileList_pos[[7]][prof1,14][[1]])-(profileList_pos[[7]][prof2,14][[1]])),digits=3)
						plot(
							int1,int2,
							pch=19,cex=.8,
							xlab="",ylab="",main=paste("Isot",round(correl,digits=3),delmass,sep=" / "))
						invisible(readline(prompt="Press [enter] to continue"))
					}					
				}
			}				
			if(length(links_profiles_pos[[n]]$homol[,1])>0){
				ff<-(ff+length(links_profiles_pos[[n]]$homol[,1]))			
				for(m in 1:length(links_profiles_pos[[n]]$homol[,1])){
					if(links_profiles_pos[[n]]$homol[m,3]>0){next}; # done before via linked profile
					prof1<-as.numeric(names(links_profiles_pos)[n])
					prof2<-links_profiles_pos[[n]]$homol[m,1]
					if(links_profiles_pos[[n]]$homol[m,2]<comp_num_cor){
						if(mark_it){
							links_profiles_pos[[n]]$homol[m,3]<-1000
							# ... and the linked profile prof2
							links_profiles_pos[[	
									profileList_pos[["index_prof"]][prof2,18]
							]]$homol[
								links_profiles_pos[[	
									profileList_pos[["index_prof"]][prof2,18]
								]]$homol[,1]==prof1
							,3]<-1000						
						}
						next;
					}
					del1<-profileList_pos[["index_prof"]][prof1,"number_peaks_total"][[1]]					
					del2<-profileList_pos[["index_prof"]][prof2,"number_peaks_total"][[1]]	
					these<-profileList_pos[["peaks"]][profileList_pos[["index_prof"]][prof1,1]:profileList_pos[["index_prof"]][prof1,2],6]
					those<-profileList_pos[["peaks"]][profileList_pos[["index_prof"]][prof2,1]:profileList_pos[["index_prof"]][prof2,2],6]
					this<-match(these,those)
					if(sum(!is.na(this))<comp_num_cor){ # in case: profiles linked more than once in the same file ...
						if(mark_it){
							links_profiles_pos[[n]]$homol[m,3]<-1000
							# ... and the linked profile prof2
							links_profiles_pos[[	
									profileList_pos[["index_prof"]][prof2,18]
							]]$homol[
								links_profiles_pos[[	
									profileList_pos[["index_prof"]][prof2,18]
								]]$homol[,1]==prof1
							,3]<-1000							
						}
						next;					
					}
					int1<-((profileList_pos[["peaks"]][profileList_pos[["index_prof"]][prof1,1]:profileList_pos[["index_prof"]][prof1,2],2])[!is.na(this)])
					int2<-((profileList_pos[["peaks"]][profileList_pos[["index_prof"]][prof2,1]:profileList_pos[["index_prof"]][prof2,2],2])[this[!is.na(this)]])					
					correl<-cor(int1,int2)
					if(is.na(correl)){stop("Correl. issue#4")}
					# mark lack of relation in this profile prof1 ...
					links_profiles_pos[[n]]$homol[m,3]<-correl
					# ... and the linked profile prof2
					links_profiles_pos[[	
							profileList_pos[["index_prof"]][prof2,18]
					]]$homol[
						links_profiles_pos[[	
							profileList_pos[["index_prof"]][prof2,18]
						]]$homol[,1]==prof1
					,3]<-correl	
					if(plot_it){
						delmass<-round(abs((profileList_pos[["index_prof"]][prof1,14][[1]])-(profileList_pos[["index_prof"]][prof2,14][[1]])),digits=3)
						plot(
							int1,int2,
							pch=19,cex=.8,
							xlab="",ylab="",main=paste("Homologues",round(correl,digits=3),delmass,sep=" / "))
						invisible(readline(prompt="Press [enter] to continue"))
					}					
				}
			}					
		}
	}
	##############################################################################
	#save(links_profiles_pos,file=file.path(as.character(logfile[[1]]),"results","componentization","links_profiles_pos"),compress=FALSE)
	#save(profileList_pos,file=file.path(as.character(logfile[[1]]),"results","profileList_pos"),compress=FALSE);
	#load(file=file.path(as.character(logfile[[1]]),"results","componentization","links_profiles_pos"));
	##############################################################################
	# build components ###########################################################
	plot_it<-FALSE
#############################
	
	verbose<-FALSE # DELETE LATER !!!!!! DEFINED IN SERVER:R

#############################
	with_test<-FALSE
	skip_peaks<-FALSE
	maxint<-rep(0,length(profileList_pos[["index_prof"]][,1]))
	in_group<-rep(0,length(profileList_pos[["index_prof"]][,1]))
	max_group<-rep(0,length(profileList_pos[["index_prof"]][,1]))
	do_group<-1
	for(i in 1:length(profileList_pos[["index_prof"]][,1])){
		maxint[i]<-max(
			profileList_pos[[2]][
				profileList_pos[["index_prof"]][i,1]:profileList_pos[["index_prof"]][i,2]
			,2]
		)
	}
	######################
	along<-order(maxint,decreasing=TRUE)
	for(i in 1:length(along)){
		if(in_group[along[i]]!=0){next}
		max_group[along[i]]<-1
		# (1) get in-/directly related isotopologue peaks ########################	
		prof_isot_IDs<-enviMass:::get_isotopol(
			profileList=profileList_pos,
			prof_ID=along[i],
			links_profiles=links_profiles_pos,
			min_peaks=3,
			skip_peaks=skip_peaks,
			min_cor=.9,
			with_test=with_test,
			only_direct=FALSE,
			del_RT=30,
			omit_profiles=in_group
		)		
		# collect adducts for these isotopologue peaks ###########################
		prof_adduct_IDs<-c()
		for(j in 1:length(prof_isot_IDs)){
			got_adducts<-enviMass:::get_adducts(
				profileList=profileList_pos,
				prof_ID=prof_isot_IDs[j],
				links_profiles=links_profiles_pos,
				min_peaks=3,
				skip_peaks=skip_peaks,
				min_cor=.9,
				with_test=with_test,
				omit_profiles=in_group
			)
			prof_adduct_IDs<-c(prof_adduct_IDs,got_adducts)
		}
		prof_all_IDs<-c(prof_isot_IDs,prof_adduct_IDs)
		if(plot_it){
			enviMass:::plot_components(
				profileList=profileList_pos,
				prof_IDs=prof_all_IDs,
				links_profiles=links_profiles_pos,
				what="relations",
				#what="profiles",
				xlim=FALSE,ylim=FALSE,await_input=TRUE,
				skipit=TRUE,
				min_peaks=3
			)	
		}	
		in_group[prof_all_IDs]<-do_group
		do_group<-(do_group+1)
	}
	cat(paste("\n Reductive fraction of profiles: ",
		as.character(1-round(sum(max_group==1)/length(max_group),digits=2))
		,sep=""))
	##############################################################################
	components<-list(0)
	components[[1]]<-cbind(profileList_pos[["index_prof"]][,"profile_ID"],maxint,in_group,max_group)
	colnames(components[[1]])<-c("Profile ID","max_int","Group number","Maximum in group")
	groups<-vector("list",max(in_group))
	for(i in 1:length(in_group)){ # make main profile come first!
		if(max_group[i]==1){
			groups[[in_group[i]]]<-c(i,groups[[in_group[i]]])
		}else{
			groups[[in_group[i]]]<-c(groups[[in_group[i]]],i)
		}
	}
	components[[2]]<-groups
	names(components)<-c("Profiles","Components")
	save(components,file=file.path(as.character(logfile[[1]]),"results","componentization","components"),compress=FALSE)
	# save peak_links
	rm(measurements)
	##############################################################################
	
}
##################################################################################	
