	########################################################################################################	
	# clean ! ##############################################################################################
	those<-list.files(file.path(logfile$project_folder,"results","screening"))
	if(length(those)>0){
		for(i in 1:length(those)){
			if(grepl("IS",those[i])){ # distinguish wfrom targets
				file.remove(file.path(logfile$project_folder,"results","screening",those[i]))
			}
		}
	}
	########################################################################################################	

	########################################################################################################
	# load available LOD smoothing spline models ###########################################################
	if(any(objects(envir=as.environment(".GlobalEnv"))=="LOD_splined")){rm(LOD_splined,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="LOD_splined")){rm(LOD_splined)}	
	if(file.exists(file=file.path(logfile$project_folder,"results","LOD","LOD_splined"))){
		load(file=file.path(logfile$project_folder,"results","LOD","LOD_splined"));
		do_LOD<-TRUE
	}else{
		do_LOD<-FALSE	
	}		
	########################################################################################################

	########################################################################################################	
	########################################################################################################
	# IS screening on positive ionization ##################################################################
	if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="peaklist")){rm(peaklist)}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_pos")){rm(profileList_pos,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="profileList_pos")){rm(profileList_pos)}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="links_peaks_pos")){rm(links_peaks_pos,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="links_peaks_pos")){rm(links_peaks_pos)}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="pattern_pos_IS")){rm(pattern_pos_IS,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="pattern_pos_IS")){rm(pattern_pos_IS)}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="patternRT_pos_IS")){rm(patternRT_pos_IS,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="patternRT_pos_IS")){rm(patternRT_pos_IS)}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="patternDelRT_pos_IS")){rm(patternDelRT_pos_IS,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="patternDelRT_pos_IS")){rm(patternDelRT_pos_IS)}
	if(
		file.exists(file.path(as.character(logfile[[1]]),"results","profileList_pos")) &
		file.exists(file.path(logfile[[1]],"results","pattern_pos_IS"))
	){

		load(file=file.path(as.character(logfile[[1]]),"results","profileList_pos"),envir=as.environment(".GlobalEnv"));
		load(file=file.path(as.character(logfile[[1]]),"results","links_peaks_pos"),envir=as.environment(".GlobalEnv"));		
		load(file=file.path(logfile[[1]],"results","pattern_pos_IS"),envir=as.environment(".GlobalEnv"));
		pattern<<-pattern_pos_IS;rm(pattern_pos_IS,envir=as.environment(".GlobalEnv"));
		load(file=file.path(logfile[[1]],"results","patternRT_pos_IS"),envir=as.environment(".GlobalEnv"));
		pattern_RT<<-patternRT_pos_IS;rm(patternRT_pos_IS,envir=as.environment(".GlobalEnv"));
		load(file=file.path(logfile[[1]],"results","patternDelRT_pos_IS"),envir=as.environment(".GlobalEnv"));
		pattern_delRT<<-patternDelRT_pos_IS;rm(patternDelRT_pos_IS,envir=as.environment(".GlobalEnv"));
		
		mztol<-as.numeric(logfile$parameters$IS_dmz)				# m/z tolerance ...
		ppm<-as.logical(as.character(logfile$parameters$IS_ppm))	# ... given in pppm?
		cutint<-as.numeric(logfile$parameters$IS_intcut)			# Lower intensity threhold
		int_tol<-as.numeric(logfile$parameters$IS_inttol)			# Intensity tolerance %
		#RT_tol_outside<-as.numeric(logfile$parameters$IS_drt1)		# RT tolerance of peaks in sample relative to their expected RT [s] - incorporated via pattern_delRT during isotopologue pattern generation
		RT_tol_inside<-as.numeric(logfile$parameters$IS_drt2)		# RT tolerance of peaks within an isotope pattern [s]
		cut_score<-as.numeric(logfile$parameters$IS_w1)	
		
		peaks<-profileList_pos[[7]];
		peaklist<-peaks[,c(14,16,15)];
		# screen centroids
		count_nonmax<-0
		for(i in 1:length(pattern)){
			count_nonmax<-(count_nonmax+
				length(pattern[[i]][,1])
			)
		}
		centro_mass<-rep(0,count_nonmax)
		centro_ID<-rep(0,count_nonmax)
		centro_maxpeak<-rep(FALSE,count_nonmax)		
		centro_number<-rep(0,count_nonmax)
		centro_RT<-rep(0,count_nonmax)
		centro_dRT<-rep(0,count_nonmax)
		at_ID<-1
		screen_list<-as.list(rep("FALSE",length(pattern)))
		for(i in 1:length(pattern)){
			n<-length(pattern[[i]][,1])
			centro_mass[at_ID:(at_ID+n-1)]<-pattern[[i]][,1]
			centro_ID[at_ID:(at_ID+n-1)]<-i
			centro_maxpeak[at_ID:(at_ID+n-1)]<-(pattern[[i]][,2]==max(pattern[[i]][,2]))			
			centro_number[at_ID:(at_ID+n-1)]<-(1:n)
			centro_RT[at_ID:(at_ID+n-1)]<-pattern_RT[i]
			centro_dRT[at_ID:(at_ID+n-1)]<-pattern_delRT[i]
			screen_list[[i]]<-as.list(rep("FALSE",n))
			at_ID<-(at_ID+n)
		}
		getit <- search_peak( 
			peaklist, 
			centro_mass, 
			dmz=mztol*4, # precheck for profiles
			ppm=ppm, 
			RT=centro_RT, 
			dRT=centro_dRT
		)	
		if(as.character(logfile$parameters$screen_IS_maxonly)=="TRUE"){ # only retain max_peak-results?
			getit[!centro_maxpeak]<-"FALSE"
		}
		if(FALSE){ # debug - retain results of getit only for a selected compound pattern
			at<-25
			getit[centro_ID!=at]<-"FALSE"
		}
		for(i in 1:length(getit)){ # transfer to a fist list of compoundadduct x centroids
			screen_list[[centro_ID[i]]][[centro_number[i]]]<-getit[i]
		}
		# resort to a full result list: pattern x sample x (centroids,matches) ( = peak index in profileList_pos)
		IS_pos_screen_listed<-list()  # default: no match at all
		for(i in 1:length(screen_list)){ # over compound x adduct = i
			if(any(is.na(screen_list[[i]]==FALSE))){
				IS_pos_screen_listed[[i]]<-list() # m-level		
				for(j in 1:length(screen_list[[i]])){ # over its centroids = j
					if(screen_list[[i]][[j]]!="FALSE"){ 
						profs<-as.numeric(strsplit(screen_list[[i]][[j]]," / ")[[1]])
						for(k in 1:length(profs)){ # over their matched profile peaks = k		
							if(profileList_pos[[7]][profs[k],4]!=profs[k]){cat("\n debug me: profile ID mismatch");stop();} # just a check
							for(m in profileList_pos[[7]][profs[k],1]:profileList_pos[[7]][profs[k],2]){ # over their sample peaks			
								delmass<-abs(profileList_pos[[2]][m,1]-pattern[[i]][j,1])		
								if(!ppm){
									if(delmass>mztol){next}
								}else{
									if((delmass*1E6/pattern[[i]][j,1])>mztol){next}
								}
								if(length(IS_pos_screen_listed[[i]])<profileList_pos[[2]][m,6][[1]] ){
									IS_pos_screen_listed[[i]][[profileList_pos[[2]][m,6][[1]]]]<-matrix(ncol=2,nrow=0)	# sample level
								}
								IS_pos_screen_listed[[i]][[profileList_pos[[2]][m,6][[1]]]]<-rbind(
									IS_pos_screen_listed[[i]][[profileList_pos[[2]][m,6][[1]]]],c(j,m)
								)
							}							
						}
					}
				}
			}else{
				IS_pos_screen_listed[[i]]<-numeric(0)	
			}
		}
		# decompose ###########################################################################		
		if(logfile$parameters$screen_IS_cutit=="TRUE"){
			use_score_cut<-TRUE;
			score_cut<-cut_score
		}else{
			use_score_cut<-FALSE;
			score_cut<-0		
		}
		many<-0
		many_unamb<-0
		res_IS_pos_screen<-list()  # default: no match at all
		if(length(IS_pos_screen_listed)>0){
			for(i in 1:length(IS_pos_screen_listed)){ # i - on compound_adduct
				if(length(IS_pos_screen_listed[[i]])>0){
					res_IS_pos_screen[[i]]<-list()
					for(m in 1:length(IS_pos_screen_listed[[i]])){ # m - sample
						if(length(IS_pos_screen_listed[[i]][[m]])>0){
							if(do_LOD){
								with_model<-which(names(LOD_splined)==paste("LOD_",m,sep=""))
								if(length(with_model)>0){						
									use_cutint<-10^(predict(LOD_splined[[with_model]],pattern_RT[i])$y)
								}else{
									cat("\n Missing LOD model; using default intensity thershold. Debug?")
									use_cutint<-cutint;
								}
							}else{
								use_cutint<-cutint
							}
							combination_matches<-recomb_score(
								cent_peak_mat=IS_pos_screen_listed[[i]][[m]],
								pattern_compound=pattern[[i]],
								profileList=profileList_pos,
								LOD=use_cutint,
								RT_tol_inside=RT_tol_inside,
								int_tol=int_tol,
								use_score_cut=FALSE,
								score_cut=score_cut,
								plotit=FALSE,
								verbose=FALSE
							)
							res_IS_pos_screen[[i]][[m]]<-combination_matches
							if(length(combination_matches)>1){many_unamb<-(many_unamb+1)}
							many<-(many+1)
						}
					}
				}else{
					res_IS_pos_screen[[i]]<-numeric(0)
				}
			}
			names(res_IS_pos_screen)<-names(IS_pos_screen_listed)
		}
		# save list ########################################################################################
		save(res_IS_pos_screen,file=file.path(logfile$project_folder,"results","screening","res_IS_pos_screen"))
		# assemble output table of length(list) ############################################################
		# iterator m is directly equal to the sample ID ####################################################
		if( length(IS_pos_screen_listed)>0 ){
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			measurements<-measurements[measurements[,8]=="TRUE",]
			intstand<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
			results_screen_IS_pos<-enviMass:::get_screening_results(
				screened_listed=res_IS_pos_screen,
				pattern=pattern,
				profileList=profileList_pos,
				measurements_table=measurements,
				compound_table=intstand,
				cut_score=cut_score
			)
			# create links in link_list_pos for peaks in profileList_pos = used for tracing back ###########
			use_entries<-find_empty(links_peaks_pos)
			for(i in 1:length(res_IS_pos_screen)){
				if(length(res_IS_pos_screen[[i]])>0){
					for(m in 1:length(res_IS_pos_screen[[i]])){
						if(length(res_IS_pos_screen[[i]][[m]])>0){
							for(k in 1:length(res_IS_pos_screen[[i]][[m]])){
								local_score<-0
								if(!is.na(res_IS_pos_screen[[i]][[m]][[k]]$score_1)){
									local_score<-(local_score+res_IS_pos_screen[[i]][[m]][[k]]$score_1)
								}
								if( (local_score>=1) || (is.na(res_IS_pos_screen[[i]][[m]][[k]]$score_1)) ){
									if(!is.na(res_IS_pos_screen[[i]][[m]][[k]]$score_2)){
										local_score<-(local_score+res_IS_pos_screen[[i]][[m]][[k]]$score_2)
									}
								}									
								if(local_score>=cut_score){
									for(a in 1:length(res_IS_pos_screen[[i]][[m]][[k]]$Peaks[,2]) ){
										# no entry for this peak in profileList<->links_peaks_pos exists yet
										if(profileList_pos[[2]][res_IS_pos_screen[[i]][[m]][[k]]$Peaks[a,2],5]==0){ 
											if(length(use_entries)>0){
												at_entry<-use_entries
												use_entries<-use_entries[-1]
											}else{
												at_entry<-(length(links_peaks_pos)+1)
											}
											links_peaks_pos[[at_entry]]<-list()
											links_peaks_pos[[at_entry]][[1]]<-list()
											links_peaks_pos[[at_entry]][[2]]<-list()
											links_peaks_pos[[at_entry]][[3]]<-list()
											profileList_pos[[2]][res_IS_pos_screen[[i]][[m]][[k]]$Peaks[a,2],5]<<-at_entry
											links_peaks_pos[[at_entry]][[1]][[1]]<-names(pattern)[i]
										# or expand existing entry
										}else{
											at_entry<-profileList_pos[[2]][res_IS_pos_screen[[i]][[m]][[k]]$Peaks[a,2],5]
											at_list<-(length(links_peaks_pos[[at_entry]][[1]])+1)
											links_peaks_pos[[at_entry]][[1]][[at_list]]<-names(pattern)[i]
										}
									}
								}
							}
						}	
					}
				}
			}
			save(profileList_pos,file=file.path(as.character(logfile[[1]]),"results","profileList_pos"))
			save(links_peaks_pos,file=file.path(as.character(logfile[[1]]),"results","links_peaks_pos"))
			################################################################################################
			save(results_screen_IS_pos,file=file.path(logfile$project_folder,"results","screening","results_screen_IS_pos"))
			rm(measurements,intstand,results_screen_IS_pos);
		}
		####################################################################################################
		rm(getit,IS_pos_screen_listed,res_IS_pos_screen)
		rm(pattern,pattern_RT,pattern_delRT,envir=as.environment(".GlobalEnv"))
		####################################################################################################
		
}		
	########################################################################################################
	########################################################################################################
	
	########################################################################################################	
	########################################################################################################
	# IS screening on negative ionization ##################################################################
	if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="peaklist")){rm(peaklist)}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_neg")){rm(profileList_neg,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="profileList_neg")){rm(profileList_neg)}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="links_peaks_neg")){rm(links_peaks_neg,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="links_peaks_neg")){rm(links_peaks_neg)}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="pattern_neg_IS")){rm(pattern_neg_IS,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="pattern_neg_IS")){rm(pattern_neg_IS)}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="patternRT_neg_IS")){rm(patternRT_neg_IS,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="patternRT_neg_IS")){rm(patternRT_neg_IS)}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="patternDelRT_neg_IS")){rm(patternDelRT_neg_IS,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="patternDelRT_neg_IS")){rm(patternDelRT_neg_IS)}
	
	if(
		file.exists(file.path(as.character(logfile[[1]]),"results","profileList_neg")) &
		file.exists(file.path(logfile[[1]],"results","pattern_neg_IS"))
	){

		load(file=file.path(as.character(logfile[[1]]),"results","profileList_neg"),envir=as.environment(".GlobalEnv"));	
		load(file=file.path(as.character(logfile[[1]]),"results","links_peaks_neg"),envir=as.environment(".GlobalEnv"));
		load(file=file.path(logfile[[1]],"results","pattern_neg_IS"),envir=as.environment(".GlobalEnv"));
		pattern<<-pattern_neg_IS;rm(pattern_neg_IS,envir=as.environment(".GlobalEnv"));
		load(file=file.path(logfile[[1]],"results","patternRT_neg_IS"),envir=as.environment(".GlobalEnv"));
		pattern_RT<<-patternRT_neg_IS;rm(patternRT_neg_IS,envir=as.environment(".GlobalEnv"));
		load(file=file.path(logfile[[1]],"results","patternDelRT_neg_IS"),envir=as.environment(".GlobalEnv"));
		pattern_delRT<<-patternDelRT_neg_IS;rm(patternDelRT_neg_IS,envir=as.environment(".GlobalEnv"));
		
		mztol<-as.numeric(logfile$parameters$IS_dmz)				# m/z tolerance ...
		ppm<-as.logical(as.character(logfile$parameters$IS_ppm))	# ... given in pppm?
		cutint<-as.numeric(logfile$parameters$IS_intcut)			# Lower intensity threhold
		int_tol<-as.numeric(logfile$parameters$IS_inttol)			# Intensity tolerance %
		#RT_tol_outside<-as.numeric(logfile$parameters$IS_drt1)		# RT tolerance of peaks in sample relative to their expected RT [s] - incorporated via pattern_delRT during isotopologue pattern generation
		RT_tol_inside<-as.numeric(logfile$parameters$IS_drt2)		# RT tolerance of peaks within an isotope pattern [s]
		cut_score<-as.numeric(logfile$parameters$IS_w1)	
		
		peaks<-profileList_neg[[7]];
		peaklist<-peaks[,c(14,16,15)];
		# screen centroids
		count_nonmax<-0
		for(i in 1:length(pattern)){
			count_nonmax<-(count_nonmax+
				length(pattern[[i]][,1])
			)
		}
		centro_mass<-rep(0,count_nonmax)
		centro_ID<-rep(0,count_nonmax)
		centro_maxpeak<-rep(FALSE,count_nonmax)	
		centro_number<-rep(0,count_nonmax)
		centro_RT<-rep(0,count_nonmax)
		centro_dRT<-rep(0,count_nonmax)
		at_ID<-1
		screen_list<-as.list(rep("FALSE",length(pattern)))
		for(i in 1:length(pattern)){
			n<-length(pattern[[i]][,1])
			centro_mass[at_ID:(at_ID+n-1)]<-pattern[[i]][,1]
			centro_ID[at_ID:(at_ID+n-1)]<-i
			centro_maxpeak[at_ID:(at_ID+n-1)]<-(pattern[[i]][,2]==max(pattern[[i]][,2]))	
			centro_number[at_ID:(at_ID+n-1)]<-(1:n)
			centro_RT[at_ID:(at_ID+n-1)]<-pattern_RT[i]
			centro_dRT[at_ID:(at_ID+n-1)]<-pattern_delRT[i]
			screen_list[[i]]<-as.list(rep("FALSE",n))
			at_ID<-(at_ID+n)
		}
		getit <- search_peak( 
			peaklist, 
			centro_mass, 
			dmz=mztol*4, # precheck
			ppm=ppm, 
			RT=centro_RT, 
			dRT=centro_dRT
		)	
		if(as.character(logfile$parameters$screen_IS_maxonly)=="TRUE"){ # only retain max_peak-results?
			getit[!centro_maxpeak]<-"FALSE"
		}
		for(i in 1:length(getit)){ # transfer to a fist list of compoundadduct x centroids
			screen_list[[centro_ID[i]]][[centro_number[i]]]<-getit[i]
		}
		# resort to a full result list: pattern x sample x (centroids,matches) ( = peak index in profileList_neg)
		IS_neg_screen_listed<-list()  # default: no match at all
		for(i in 1:length(screen_list)){ # over compound x adduct = i
			if(any(is.na(screen_list[[i]]==FALSE))){
				IS_neg_screen_listed[[i]]<-list() # m-level		
				for(j in 1:length(screen_list[[i]])){ # over its centroids = j
					if(screen_list[[i]][[j]]!="FALSE"){ 
						profs<-as.numeric(strsplit(screen_list[[i]][[j]]," / ")[[1]])
						for(k in 1:length(profs)){ # over their matched profile peaks = k
							if(profileList_neg[[7]][profs[k],4]!=profs[k]){cat("\n debug me: profile ID mismatch");stop();} # just a check
							for(m in profileList_neg[[7]][profs[k],1]:profileList_neg[[7]][profs[k],2]){ # over their sample peaks
								delmass<-abs(profileList_neg[[2]][m,1]-pattern[[i]][j,1])		
								if(!ppm){
									if(delmass>mztol){next}
								}else{
									if(delmass*1E6/pattern[[i]][j,1]>mztol){next}
								}
								if(length(IS_neg_screen_listed[[i]])<profileList_neg[[2]][m,6][[1]] ){
									IS_neg_screen_listed[[i]][[profileList_neg[[2]][m,6][[1]]]]<-matrix(ncol=2,nrow=0)	# sample level
								}
								IS_neg_screen_listed[[i]][[profileList_neg[[2]][m,6][[1]]]]<-rbind(
									IS_neg_screen_listed[[i]][[profileList_neg[[2]][m,6][[1]]]],c(j,m)
								)
							}							
						}
					}
				}
			}else{
				IS_neg_screen_listed[[i]]<-numeric(0)
			}
		}
		# decompose ###########################################################################		
		if(logfile$parameters$screen_IS_cutit=="TRUE"){
			use_score_cut<-TRUE;
			score_cut<-cut_score
		}else{
			use_score_cut<-FALSE;
			score_cut<-0		
		}
		many<-0
		many_unamb<-0
		res_IS_neg_screen<-list()  # default: no match at all
		if(length(IS_neg_screen_listed)>0){
			j<-1
			for(i in j:length(IS_neg_screen_listed)){ # i - on compound_adduct
				if(length(IS_neg_screen_listed[[i]])>0){	
					res_IS_neg_screen[[i]]<-list()
					for(m in 1:length(IS_neg_screen_listed[[i]])){ # m - sample
						if(length(IS_neg_screen_listed[[i]][[m]])>0){
							if(do_LOD){
								with_model<-which(names(LOD_splined)==paste("LOD_",m,sep=""))
								if(length(with_model)>0){						
									use_cutint<-10^(predict(LOD_splined[[with_model]],pattern_RT[i])$y)
								}else{
									cat("\n Missing LOD model; using default intensity thershold. Debug?")
									use_cutint<-cutint;
								}
							}else{
								use_cutint<-cutint
							}
							combination_matches<-recomb_score(
								cent_peak_mat=IS_neg_screen_listed[[i]][[m]],
								pattern_compound=pattern[[i]],
								profileList=profileList_neg,
								LOD=use_cutint,
								RT_tol_inside=RT_tol_inside,
								int_tol=int_tol,
								use_score_cut=use_score_cut,
								score_cut=score_cut,
								plotit=FALSE,
								verbose=FALSE
							)
							res_IS_neg_screen[[i]][[m]]<-combination_matches
							if(length(combination_matches)>1){many_unamb<-(many_unamb+1)}
							many<-(many+1)
						}
					}
				}else{
					res_IS_neg_screen[[i]]<-numeric(0)
				}
			}
			names(res_IS_neg_screen)<-names(IS_neg_screen_listed)
		}
		# save list ########################################################################################
		save(res_IS_neg_screen,file=file.path(logfile$project_folder,"results","screening","res_IS_neg_screen"))
		# assemble output table of length(list) ############################################################
		# iterator m is directly equal to the sample ID ####################################################
		if(length(IS_neg_screen_listed)>0){
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			measurements<-measurements[measurements[,8]=="TRUE",]
			intstand<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
			results_screen_IS_neg<-enviMass:::get_screening_results(
				screened_listed=res_IS_neg_screen,
				pattern=pattern,
				profileList=profileList_neg,
				measurements_table=measurements,
				compound_table=intstand,
				cut_score=cut_score
			)
			# create links in link_list_pos for peaks in profileList_pos = used for tracing back ###########
			use_entries<-find_empty(links_peaks_neg)
			for(i in 1:length(res_IS_neg_screen)){
				if(length(res_IS_neg_screen[[i]])>0){
					for(m in 1:length(res_IS_neg_screen[[i]])){
						if(length(res_IS_neg_screen[[i]][[m]])>0){
							for(k in 1:length(res_IS_neg_screen[[i]][[m]])){
								local_score<-0
								if(!is.na(res_IS_neg_screen[[i]][[m]][[k]]$score_1)){
									local_score<-(local_score+res_IS_neg_screen[[i]][[m]][[k]]$score_1)
								}
								if( (local_score>=1) || (is.na(res_IS_neg_screen[[i]][[m]][[k]]$score_1)) ){
									if(!is.na(res_IS_neg_screen[[i]][[m]][[k]]$score_2)){
										local_score<-(local_score+res_IS_neg_screen[[i]][[m]][[k]]$score_2)
									}
								}									
								if(local_score>=cut_score){
									for(a in 1:length(res_IS_neg_screen[[i]][[m]][[k]]$Peaks[,2]) ){
										# no entry for this peak in profileList<->links_peaks_pos exists yet
										if(profileList_neg[[2]][res_IS_neg_screen[[i]][[m]][[k]]$Peaks[a,2],5]==0){ 
											if(length(use_entries)>0){
												at_entry<-use_entries
												use_entries<-use_entries[-1]
											}else{
												at_entry<-(length(links_peaks_neg)+1)
											}
											links_peaks_neg[[at_entry]]<-list()
											links_peaks_neg[[at_entry]][[1]]<-list()
											links_peaks_neg[[at_entry]][[2]]<-list()
											links_peaks_neg[[at_entry]][[3]]<-list()
											profileList_neg[[2]][res_IS_neg_screen[[i]][[m]][[k]]$Peaks[a,2],5]<<-at_entry
											links_peaks_neg[[at_entry]][[1]][[1]]<-names(pattern)[i]
										# or expand existing entry
										}else{
											at_entry<-profileList_neg[[2]][res_IS_neg_screen[[i]][[m]][[k]]$Peaks[a,2],5]
											at_list<-(length(links_peaks_neg[[at_entry]][[1]])+1)
											links_peaks_neg[[at_entry]][[1]][[at_list]]<-names(pattern)[i]
										}
									}
								}
							}
						}	
					}
				}
			}
			save(profileList_neg,file=file.path(as.character(logfile[[1]]),"results","profileList_neg"))
			save(links_peaks_neg,file=file.path(as.character(logfile[[1]]),"results","links_peaks_neg"))
			################################################################################################
			save(results_screen_IS_neg,file=file.path(logfile$project_folder,"results","screening","results_screen_IS_neg"))
			rm(measurements,intstand,results_screen_IS_neg);
		}
		rm(getit,IS_neg_screen_listed,res_IS_neg_screen);
		####################################################################################################
		rm(pattern,pattern_RT,pattern_delRT,envir=as.environment(".GlobalEnv"))

}		
	########################################################################################################
	########################################################################################################
	
	
	
	