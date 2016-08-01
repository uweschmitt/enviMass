
	######################################################################################################################
	# load available LOD smoothing spline models #########################################################################
	if(any(objects(envir=as.environment(".GlobalEnv"))=="LOD_splined")){rm(LOD_splined,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="LOD_splined")){rm(LOD_splined)}	
	if(file.exists(file=file.path(logfile$project_folder,"results","LOD","LOD_splined"))){
		load(file=file.path(logfile$project_folder,"results","LOD","LOD_splined"));
		do_LOD<-TRUE
	}else{
		do_LOD<-FALSE	
	}		
	######################################################################################################################
    measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	measurements<-measurements[measurements[,8]=="TRUE",,drop=FALSE]
	measurements<-measurements[measurements[,3]=="calibration",,drop=FALSE]	
	######################################################################################################################
	
	
	# POSITIVE ###########################################################################################################
	if(any(measurements[,4]=="positive" & measurements[,3]=="calibration")){

		##################################################################################################################	
		# CLEAN ALL PREVIOUS RESULTS #####################################################################################
		if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="peaklist")){rm(peaklist)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_pos_cal")){rm(profileList_pos_cal,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="profileList_pos_cal")){rm(profileList_pos_cal)}		
		if(any(objects(envir=as.environment(".GlobalEnv"))=="pattern_pos_IS")){rm(pattern_pos_IS,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="pattern_pos_IS")){rm(pattern_pos_IS)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="patternRT_pos_IS")){rm(patternRT_pos_IS,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="patternRT_pos_IS")){rm(patternRT_pos_IS)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="patternDelRT_pos_IS")){rm(patternDelRT_pos_IS,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="patternDelRT_pos_IS")){rm(patternDelRT_pos_IS)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="pattern")){rm(pattern,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="pattern")){rm(pattern)}	
		if(any(objects(envir=as.environment(".GlobalEnv"))=="pattern_RT")){rm(pattern_RT,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="pattern_RT")){rm(pattern_RT)}		
		if(any(objects(envir=as.environment(".GlobalEnv"))=="pattern_delRT")){rm(pattern_delRT,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="pattern_delRT")){rm(pattern_delRT)}		
		if(file.exists(file=file.path(as.character(logfile[[1]]),"quantification","profileList_pos_cal"))){
			file.remove(file=file.path(as.character(logfile[[1]]),"quantification","profileList_pos_cal"))
		}	
		if(file.exists(file=file.path(as.character(logfile[[1]]),"quantification","res_IS_pos_screen_cal"))){
			file.remove(file=file.path(as.character(logfile[[1]]),"quantification","res_IS_pos_screen_cal"))
		}
		if(file.exists(file=file.path(as.character(logfile[[1]]),"quantification","res_target_pos_screen_cal"))){
			file.remove(file=file.path(as.character(logfile[[1]]),"quantification","res_target_pos_screen_cal"))
		}
		if(file.exists(file=file.path(as.character(logfile[[1]]),"quantification","results_screen_IS_pos_cal"))){
			file.remove(file=file.path(as.character(logfile[[1]]),"quantification","results_screen_IS_pos_cal"))
		}
		if(file.exists(file=file.path(as.character(logfile[[1]]),"quantification","results_screen_target_pos_cal"))){
			file.remove(file=file.path(as.character(logfile[[1]]),"quantification","results_screen_target_pos_cal"))
		}
		##################################################################################################################
		# PROFILING ######################################################################################################		
		profileList_pos_cal<-startprofiles(
							logfile,
							frac=FALSE,
							sets=FALSE,
							progbar=logfile$parameters$progressBar,
							ion_mode="positive",
							until=FALSE,
							selective=FALSE,
							types=c("calibration")
						)							
		profileList_pos_cal<-agglomer(
							profileList_pos_cal,
							dmass=(as.numeric(logfile$parameters$prof_dmz)+1),
							ppm=as.logical(as.character(logfile$parameters$prof_ppm)),
							dret=(as.numeric(logfile$parameters$prof_drt)+10)
						)
		profileList_pos_cal<-partcluster(
							profileList=profileList_pos_cal,
							dmass=as.numeric(logfile$parameters$prof_dmz),
							ppm=as.logical(as.character(logfile$parameters$prof_ppm)),
							dret=as.numeric(logfile$parameters$prof_drt),
							from=FALSE,
							to=FALSE,
							progbar=logfile$parameters$progressBar,
							plotit=FALSE,
							replicates=FALSE
						)
		profileList_pos_cal<<-profileList_pos_cal
		save(profileList_pos_cal,file=file.path(as.character(logfile[[1]]),"quantification","profileList_pos_cal"),compress=FALSE);
		##################################################################################################################
		# IS SCREENING ###################################################################################################
		#load(file=file.path(as.character(logfile[[1]]),"quantification","profileList_pos_cal"),envir=as.environment(".GlobalEnv"));
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
		# include: restrict pattern only to the adducts used for screening! - but then they dont show up the screening results ...
		
		peaks<-profileList_pos_cal[[7]];
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
		if( as.character(logfile$parameters$screen_IS_maxonly)=="TRUE" ){ # only retain max_peak-results?
			getit[!centro_maxpeak]<-"FALSE"
		}
		for(i in 1:length(getit)){ # transfer to a fist list of compoundadduct x centroids
			screen_list[[centro_ID[i]]][[centro_number[i]]]<-getit[i]
		}
		# resort to a full result list: pattern x sample x (centroids,matches) ( = peak index in profileList_pos_cal)		
		IS_pos_screen_listed_cal<-list()  # default: no match at all	
		set_ID<-seq(1:length(profileList_pos_cal[[4]]))	
		for(i in 1:length(screen_list)){ # over compound x adduct = i
			if(any(is.na(screen_list[[i]]==FALSE))){
				IS_pos_screen_listed_cal[[i]]<-list() # m-level		
				for(j in 1:length(screen_list[[i]])){ # over its centroids = j
					if(screen_list[[i]][[j]]!="FALSE"){ 
						profs<-as.numeric(strsplit(screen_list[[i]][[j]]," / ")[[1]])
						for(k in 1:length(profs)){ # over their matched profile peaks = k		
							if(profileList_pos_cal[[7]][profs[k],4]!=profs[k]){cat("\n debug me: profile ID mismatch");stop();} # a check
							for(m in profileList_pos_cal[[7]][profs[k],1]:profileList_pos_cal[[7]][profs[k],2]){ # over their sample peaks			
								delmass<-abs(profileList_pos_cal[[2]][m,1]-pattern[[i]][j,1])		
								if(!ppm){
									if(delmass>mztol){next}
								}else{
									if((delmass*1E6/pattern[[i]][j,1])>mztol){next}
								}
								at_ID<-set_ID[profileList_pos_cal[[4]]==as.character(profileList_pos_cal[[2]][m,6])]								
								if(length(IS_pos_screen_listed_cal[[i]])<at_ID){							
									IS_pos_screen_listed_cal[[i]][[at_ID]]<-matrix(ncol=2,nrow=0)	
								}else{
									if(length(IS_pos_screen_listed_cal[[i]][[at_ID]])==0){
										IS_pos_screen_listed_cal[[i]][[at_ID]]<-matrix(ncol=2,nrow=0)	
									}
								}
								IS_pos_screen_listed_cal[[i]][[at_ID]]<-rbind(
									IS_pos_screen_listed_cal[[i]][[at_ID]],c(j,m)
								)					
								colnames(IS_pos_screen_listed_cal[[i]][[at_ID]])<-c(as.character(profileList_pos_cal[[4]][at_ID]),"")								
							}							
						}
					}
				}
			}else{
				IS_pos_screen_listed_cal[[i]]<-numeric(0)	
			}
		}	
		# decompose ###########################################################################		
		if( logfile$parameters$screen_IS_cutit=="TRUE" ){
			use_score_cut<-TRUE;
			score_cut<-cut_score
		}else{
			use_score_cut<-FALSE;
			score_cut<-0		
		}
		many<-0
		many_unamb<-0
		res_IS_pos_screen_cal<-list()  # default: no match at all
		if(length(IS_pos_screen_listed_cal)>0){
			for(i in 1:length(IS_pos_screen_listed_cal)){ # i - on compound_adduct
				if(length(IS_pos_screen_listed_cal[[i]])>0){
					res_IS_pos_screen_cal[[i]]<-list()
					for(m in 1:length(IS_pos_screen_listed_cal[[i]])){ # m (relates to IDs in profileList_pos_cal[[4]])				
						at_ID<-set_ID[profileList_pos_cal[[4]]==colnames(IS_pos_screen_listed_cal[[i]][[m]])[1]]					
						if(length(IS_pos_screen_listed_cal[[i]][[m]])>0){
							if(do_LOD){							
								with_model<-which(names(LOD_splined)==paste("LOD_",colnames(IS_pos_screen_listed_cal[[i]][[m]])[1],sep=""))						
								if(length(with_model)>0){						
									use_cutint<-10^(predict(LOD_splined[[with_model]],pattern_RT[i])$y)
								}else{
									cat("\n Missing LOD model; using default intensity threshold. Debug?")
									use_cutint<-cutint;
								}
							}else{
								use_cutint<-cutint
							}								
							combination_matches<-recomb_score(
								cent_peak_mat=IS_pos_screen_listed_cal[[i]][[m]],
								pattern_compound=pattern[[i]],
								profileList=profileList_pos_cal,
								LOD=use_cutint,
								RT_tol_inside=RT_tol_inside,
								int_tol=int_tol,
								use_score_cut=use_score_cut,
								score_cut=score_cut,
								plotit=FALSE,
								verbose=FALSE,
								RT_seperate=TRUE
							)
							for(k in 1:length(combination_matches)){ # add file ID
								combination_matches[[k]][[10]]<-colnames(IS_pos_screen_listed_cal[[i]][[m]])[1]
								names(combination_matches[[k]])[10]<-"file_ID"
							}
							res_IS_pos_screen_cal[[i]][[at_ID]]<-combination_matches
							if(length(combination_matches)>1){many_unamb<-(many_unamb+1)}
							many<-(many+1)
						}
					}
				}else{
					res_IS_pos_screen_cal[[i]]<-numeric(0)
				}
			}
			names(res_IS_pos_screen_cal)<-names(IS_pos_screen_listed_cal)
		}
		##################################################################################################################
		# save list ######################################################################################################
		names(res_IS_pos_screen_cal)<-names(pattern)
		save(res_IS_pos_screen_cal,file=file.path(logfile$project_folder,"quantification","res_IS_pos_screen_cal"))
		# assemble output table of length(list) ##########################################################################
		# iterator m is directly equal to the sample ID ##################################################################
		if( length(IS_pos_screen_listed_cal)>0 ){
			intstand<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
			results_screen_IS_pos_cal<-enviMass:::get_screening_results(
				screened_listed=res_IS_pos_screen_cal,
				pattern=pattern,
				at_RT=pattern_RT,
				profileList=profileList_pos_cal,
				measurements_table=measurements,
				compound_table=intstand,
				cut_score=cut_score
			)
			save(results_screen_IS_pos_cal,file=file.path(logfile$project_folder,"quantification","results_screen_IS_pos_cal"))
			rm(intstand,results_screen_IS_pos_cal);
		}
		##################################################################################################################
		rm(getit,IS_pos_screen_listed_cal,res_IS_pos_screen_cal)
		rm(pattern,pattern_RT,pattern_delRT,envir=as.environment(".GlobalEnv"))
		##################################################################################################################

		##################################################################################################################
		# TARGET SCREENING ###############################################################################################
		#load(file=file.path(as.character(logfile[[1]]),"quantification","profileList_pos_cal"),envir=as.environment(".GlobalEnv"));
		load(file=file.path(logfile[[1]],"results","pattern_pos_target"),envir=as.environment(".GlobalEnv"));
		pattern<<-pattern_pos_target;rm(pattern_pos_target,envir=as.environment(".GlobalEnv"));
		load(file=file.path(logfile[[1]],"results","patternRT_pos_target"),envir=as.environment(".GlobalEnv"));
		pattern_RT<<-patternRT_pos_target;rm(patternRT_pos_target,envir=as.environment(".GlobalEnv"));
		load(file=file.path(logfile[[1]],"results","patternDelRT_pos_target"),envir=as.environment(".GlobalEnv"));
		pattern_delRT<<-patternDelRT_pos_target;rm(patternDelRT_pos_target,envir=as.environment(".GlobalEnv"));
		targets<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");
		
		mztol<-as.numeric(logfile$parameters$tar_dmz)				# m/z tolerance ...
		ppm<-as.logical(as.character(logfile$parameters$tar_ppm))	# ... given in pppm?
		cutint<-as.numeric(logfile$parameters$tar_intcut)			# Lower intensity threhold
		int_tol<-as.numeric(logfile$parameters$tar_inttol)			# Intensity tolerance %
		#RT_tol_outside<-as.numeric(logfile$parameters$target_drt1)		# RT tolerance of peaks in sample relative to their expected RT [s] - incorporated via pattern_delRT during isotopologue pattern generation
		RT_tol_inside<-as.numeric(logfile$parameters$tar_drt2)		# RT tolerance of peaks within an isotope pattern [s]
		cut_score<-as.numeric(logfile$parameters$tar_w1)	
		
		peaks<-profileList_pos_cal[[7]];
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
			forID<-strsplit(names(pattern)[i],"_")[[1]][1]
			if(targets[targets[,1]==forID,6]!="FALSE"){ # only include compounds used for quantification
				centro_mass[at_ID:(at_ID+n-1)]<-pattern[[i]][,1]
				centro_maxpeak[at_ID:(at_ID+n-1)]<-(pattern[[i]][,2]==max(pattern[[i]][,2]))			
				centro_RT[at_ID:(at_ID+n-1)]<-pattern_RT[i]
				centro_dRT[at_ID:(at_ID+n-1)]<-pattern_delRT[i]
			}
			screen_list[[i]]<-as.list(rep("FALSE",n))
			centro_ID[at_ID:(at_ID+n-1)]<-i
			centro_number[at_ID:(at_ID+n-1)]<-(1:n)
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
		if( as.character(logfile$parameters$screen_target_maxonly)=="TRUE" ){ # only retain max_peak-results?
			getit[!centro_maxpeak]<-"FALSE"
		}
		for(i in 1:length(getit)){ # transfer to a fist list of compoundadduct x centroids
			screen_list[[centro_ID[i]]][[centro_number[i]]]<-getit[i]
		}
		# resort to a full result list: pattern x sample x (centroids,matches) ( = peak index in profileList_pos_cal)		
		target_pos_screen_listed_cal<-list()  # default: no match at all	
		set_ID<-seq(1:length(profileList_pos_cal[[4]]))	
		for(i in 1:length(screen_list)){ # over compound x adduct = i
			if(any(is.na(screen_list[[i]]==FALSE))){
				target_pos_screen_listed_cal[[i]]<-list() # m-level		
				for(j in 1:length(screen_list[[i]])){ # over its centroids = j
					if(screen_list[[i]][[j]]!="FALSE"){ 
						profs<-as.numeric(strsplit(screen_list[[i]][[j]]," / ")[[1]])
						for(k in 1:length(profs)){ # over their matched profile peaks = k		
							if(profileList_pos_cal[[7]][profs[k],4]!=profs[k]){cat("\n debug me: profile ID mismatch");stop();} # a check
							for(m in profileList_pos_cal[[7]][profs[k],1]:profileList_pos_cal[[7]][profs[k],2]){ # over their sample peaks			
								delmass<-abs(profileList_pos_cal[[2]][m,1]-pattern[[i]][j,1])		
								if(!ppm){
									if(delmass>mztol){next}
								}else{
									if((delmass*1E6/pattern[[i]][j,1])>mztol){next}
								}
								at_ID<-set_ID[profileList_pos_cal[[4]]==as.character(profileList_pos_cal[[2]][m,6])]								
								if(length(target_pos_screen_listed_cal[[i]])<at_ID){							
									target_pos_screen_listed_cal[[i]][[at_ID]]<-matrix(ncol=2,nrow=0)	
								}else{
									if(length(target_pos_screen_listed_cal[[i]][[at_ID]])==0){
										target_pos_screen_listed_cal[[i]][[at_ID]]<-matrix(ncol=2,nrow=0)	
									}
								}
								target_pos_screen_listed_cal[[i]][[at_ID]]<-rbind(
									target_pos_screen_listed_cal[[i]][[at_ID]],c(j,m)
								)					
								colnames(target_pos_screen_listed_cal[[i]][[at_ID]])<-c(as.character(profileList_pos_cal[[4]][at_ID]),"")							
							}							
						}
					}
				}
			}else{
				target_pos_screen_listed_cal[[i]]<-numeric(0)	
			}
		}	
		# decompose ###########################################################################		
		if( logfile$parameters$screen_target_cutit=="TRUE" ){
			use_score_cut<-TRUE;
			score_cut<-cut_score
		}else{
			use_score_cut<-FALSE;
			score_cut<-0		
		}
		many<-0
		many_unamb<-0
		res_target_pos_screen_cal<-list()  # default: no match at all
		if(length(target_pos_screen_listed_cal)>0){
			for(i in 1:length(target_pos_screen_listed_cal)){ # i - on compound_adduct
				if(length(target_pos_screen_listed_cal[[i]])>0){
					res_target_pos_screen_cal[[i]]<-list()
					for(m in 1:length(target_pos_screen_listed_cal[[i]])){ # m (relates to IDs in profileList_pos_cal[[4]])					
						at_ID<-set_ID[profileList_pos_cal[[4]]==colnames(target_pos_screen_listed_cal[[i]][[m]])[1]]					
						if(length(target_pos_screen_listed_cal[[i]][[m]])>0){
							if(do_LOD){							
								with_model<-which(names(LOD_splined)==paste("LOD_",colnames(target_pos_screen_listed_cal[[i]][[m]])[1],sep=""))					
								if(length(with_model)>0){						
									use_cutint<-10^(predict(LOD_splined[[with_model]],pattern_RT[i])$y)
								}else{
									cat("\n Missing LOD model; using default intensity threshold. Debug?")
									use_cutint<-cutint;
								}
							}else{
								use_cutint<-cutint
							}							
							combination_matches<-recomb_score(
								cent_peak_mat=target_pos_screen_listed_cal[[i]][[m]],
								pattern_compound=pattern[[i]],
								profileList=profileList_pos_cal,
								LOD=use_cutint,
								RT_tol_inside=RT_tol_inside,
								int_tol=int_tol,
								use_score_cut=use_score_cut,
								score_cut=score_cut,
								plotit=FALSE,
								verbose=FALSE,
								RT_seperate=TRUE
							)
							for(k in 1:length(combination_matches)){ # add file ID
								combination_matches[[k]][[10]]<-colnames(target_pos_screen_listed_cal[[i]][[m]])[1]
								names(combination_matches[[k]])[10]<-"file_ID"
							}
							res_target_pos_screen_cal[[i]][[at_ID]]<-combination_matches
							if(length(combination_matches)>1){many_unamb<-(many_unamb+1)}
							many<-(many+1)
						}
					}
				}else{
					res_target_pos_screen_cal[[i]]<-numeric(0)
				}
			}
			names(res_target_pos_screen_cal)<-names(target_pos_screen_listed_cal)
		}
		##################################################################################################################
		# save list ######################################################################################################
		names(res_target_pos_screen_cal)<-names(pattern)
		save(res_target_pos_screen_cal,file=file.path(logfile$project_folder,"quantification","res_target_pos_screen_cal"))
		# assemble output table of length(list) ##########################################################################
		# iterator m is directly equal to the sample ID ##################################################################
		if( length(target_pos_screen_listed_cal)>0 ){
			targets<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");
			results_screen_target_pos_cal<-enviMass:::get_screening_results(
				screened_listed=res_target_pos_screen_cal,
				pattern=pattern,
				at_RT=pattern_RT,
				profileList=profileList_pos_cal,
				measurements_table=measurements,
				compound_table=targets,
				cut_score=cut_score
			)
			save(results_screen_target_pos_cal,file=file.path(logfile$project_folder,"quantification","results_screen_target_pos_cal"))
			rm(targets,results_screen_target_pos_cal);
		}
		##################################################################################################################
		rm(getit,target_pos_screen_listed_cal,res_target_pos_screen_cal)
		rm(pattern,pattern_RT,pattern_delRT,envir=as.environment(".GlobalEnv"))
		##################################################################################################################

	} # on positive calibration

	
	# NEGATIVE ###########################################################################################################
	if(any(measurements[,4]=="negative" & measurements[,3]=="calibration")){

		##################################################################################################################	
		# CLEAN ALL PREVIOUS RESULTS #####################################################################################
		if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="peaklist")){rm(peaklist)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_neg_cal")){rm(profileList_neg_cal,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="profileList_neg_cal")){rm(profileList_neg_cal)}		
		if(any(objects(envir=as.environment(".GlobalEnv"))=="pattern_neg_IS")){rm(pattern_neg_IS,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="pattern_neg_IS")){rm(pattern_neg_IS)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="patternRT_neg_IS")){rm(patternRT_neg_IS,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="patternRT_neg_IS")){rm(patternRT_neg_IS)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="patternDelRT_neg_IS")){rm(patternDelRT_neg_IS,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="patternDelRT_neg_IS")){rm(patternDelRT_neg_IS)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="pattern")){rm(pattern,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="pattern")){rm(pattern)}	
		if(any(objects(envir=as.environment(".GlobalEnv"))=="pattern_RT")){rm(pattern_RT,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="pattern_RT")){rm(pattern_RT)}		
		if(any(objects(envir=as.environment(".GlobalEnv"))=="pattern_delRT")){rm(pattern_delRT,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="pattern_delRT")){rm(pattern_delRT)}		
		if(file.exists(file=file.path(as.character(logfile[[1]]),"quantification","profileList_neg_cal"))){
			file.remove(file=file.path(as.character(logfile[[1]]),"quantification","profileList_neg_cal"))
		}	
		if(file.exists(file=file.path(as.character(logfile[[1]]),"quantification","res_IS_neg_screen_cal"))){
			file.remove(file=file.path(as.character(logfile[[1]]),"quantification","res_IS_neg_screen_cal"))
		}
		if(file.exists(file=file.path(as.character(logfile[[1]]),"quantification","res_target_neg_screen_cal"))){
			file.remove(file=file.path(as.character(logfile[[1]]),"quantification","res_target_neg_screen_cal"))
		}
		if(file.exists(file=file.path(as.character(logfile[[1]]),"quantification","results_screen_IS_neg_cal"))){
			file.remove(file=file.path(as.character(logfile[[1]]),"quantification","results_screen_IS_neg_cal"))
		}
		if(file.exists(file=file.path(as.character(logfile[[1]]),"quantification","results_screen_target_neg_cal"))){
			file.remove(file=file.path(as.character(logfile[[1]]),"quantification","results_screen_target_neg_cal"))
		}
		##################################################################################################################
		# PROFILING ######################################################################################################
		profileList_neg_cal<-startprofiles(
							logfile,
							frac=FALSE,
							sets=FALSE,
							progbar=logfile$parameters$progressBar,
							ion_mode="negative",
							until=FALSE,
							selective=FALSE,
							types=c("calibration")
						)
		profileList_neg_cal<-agglomer(
							profileList_neg_cal,
							dmass=(as.numeric(logfile$parameters$prof_dmz)+1),
							ppm=as.logical(as.character(logfile$parameters$prof_ppm)),
							dret=(as.numeric(logfile$parameters$prof_drt)+10)
						)
		profileList_neg_cal<-partcluster(
							profileList=profileList_neg_cal,
							dmass=as.numeric(logfile$parameters$prof_dmz),
							ppm=as.logical(as.character(logfile$parameters$prof_ppm)),
							dret=as.numeric(logfile$parameters$prof_drt),
							from=FALSE,
							to=FALSE,
							progbar=logfile$parameters$progressBar,
							plotit=FALSE,
							replicates=FALSE
						)
		profileList_neg_cal<<-profileList_neg_cal
		save(profileList_neg_cal,file=file.path(as.character(logfile[[1]]),"quantification","profileList_neg_cal"),compress=FALSE);
		##################################################################################################################
		# IS SCREENING ###################################################################################################
		#load(file=file.path(as.character(logfile[[1]]),"quantification","profileList_neg_cal"),envir=as.environment(".GlobalEnv"));
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
		
		peaks<-profileList_neg_cal[[7]];
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
		if( as.character(logfile$parameters$screen_IS_maxonly)=="TRUE" ){ # only retain max_peak-results?
			getit[!centro_maxpeak]<-"FALSE"
		}
		for(i in 1:length(getit)){ # transfer to a fist list of compoundadduct x centroids
			screen_list[[centro_ID[i]]][[centro_number[i]]]<-getit[i]
		}
		# resort to a full result list: pattern x sample x (centroids,matches) ( = peak index in profileList_neg_cal)		
		IS_neg_screen_listed_cal<-list()  # default: no match at all	
		set_ID<-seq(1:length(profileList_neg_cal[[4]]))		
		for(i in 1:length(screen_list)){ # over compound x adduct = i
			if(any(is.na(screen_list[[i]]==FALSE))){
				IS_neg_screen_listed_cal[[i]]<-list() # m-level		
				for(j in 1:length(screen_list[[i]])){ # over its centroids = j
					if(screen_list[[i]][[j]]!="FALSE"){ 
						profs<-as.numeric(strsplit(screen_list[[i]][[j]]," / ")[[1]])
						for(k in 1:length(profs)){ # over their matched profile peaks = k		
							if(profileList_neg_cal[[7]][profs[k],4]!=profs[k]){cat("\n debug me: profile ID mismatch");stop();} # a check
							for(m in profileList_neg_cal[[7]][profs[k],1]:profileList_neg_cal[[7]][profs[k],2]){ # over their sample peaks			
								delmass<-abs(profileList_neg_cal[[2]][m,1]-pattern[[i]][j,1])		
								if(!ppm){
									if(delmass>mztol){next}
								}else{
									if((delmass*1E6/pattern[[i]][j,1])>mztol){next}
								}
								at_ID<-set_ID[profileList_neg_cal[[4]]==as.character(profileList_neg_cal[[2]][m,6])]								
								if(length(IS_neg_screen_listed_cal[[i]])<at_ID){							
									IS_neg_screen_listed_cal[[i]][[at_ID]]<-matrix(ncol=2,nrow=0)	
								}else{
									if(length(IS_neg_screen_listed_cal[[i]][[at_ID]])==0){
										IS_neg_screen_listed_cal[[i]][[at_ID]]<-matrix(ncol=2,nrow=0)	
									}
								}
								IS_neg_screen_listed_cal[[i]][[at_ID]]<-rbind(
									IS_neg_screen_listed_cal[[i]][[at_ID]],c(j,m)
								)					
								colnames(IS_neg_screen_listed_cal[[i]][[at_ID]])<-c(as.character(profileList_neg_cal[[4]][at_ID]),"")								
							}							
						}
					}
				}
			}else{
				IS_neg_screen_listed_cal[[i]]<-numeric(0)	
			}
		}	
		# decompose ###########################################################################		
		if( logfile$parameters$screen_IS_cutit=="TRUE" ){
			use_score_cut<-TRUE;
			score_cut<-cut_score
		}else{
			use_score_cut<-FALSE;
			score_cut<-0		
		}
		many<-0
		many_unamb<-0
		res_IS_neg_screen_cal<-list()  # default: no match at all
		if(length(IS_neg_screen_listed_cal)>0){
			for(i in 1:length(IS_neg_screen_listed_cal)){ # i - on compound_adduct
				if(length(IS_neg_screen_listed_cal[[i]])>0){
					res_IS_neg_screen_cal[[i]]<-list()
					for(m in 1:length(IS_neg_screen_listed_cal[[i]])){ # m (relates to IDs in profileList_neg_cal[[4]])					
						at_ID<-set_ID[profileList_neg_cal[[4]]==colnames(IS_neg_screen_listed_cal[[i]][[m]])[1]]					
						if(length(IS_neg_screen_listed_cal[[i]][[m]])>0){
							if(do_LOD){							
								with_model<-which(names(LOD_splined)==paste("LOD_",colnames(IS_neg_screen_listed_cal[[i]][[m]])[1],sep=""))							
								if(length(with_model)>0){						
									use_cutint<-10^(predict(LOD_splined[[with_model]],pattern_RT[i])$y)
								}else{
									cat("\n Missing LOD model; using default intensity threshold. Debug?")
									use_cutint<-cutint;
								}
							}else{
								use_cutint<-cutint
							}							
							combination_matches<-recomb_score(
								cent_peak_mat=IS_neg_screen_listed_cal[[i]][[m]],
								pattern_compound=pattern[[i]],
								profileList=profileList_neg_cal,
								LOD=use_cutint,
								RT_tol_inside=RT_tol_inside,
								int_tol=int_tol,
								use_score_cut=use_score_cut,
								score_cut=score_cut,
								plotit=FALSE,
								verbose=FALSE,
								RT_seperate=TRUE
							)
							for(k in 1:length(combination_matches)){ # add file ID
								combination_matches[[k]][[10]]<-colnames(IS_neg_screen_listed_cal[[i]][[m]])[1]
								names(combination_matches[[k]])[10]<-"file_ID"
							}
							res_IS_neg_screen_cal[[i]][[at_ID]]<-combination_matches
							if(length(combination_matches)>1){many_unamb<-(many_unamb+1)}
							many<-(many+1)
						}
					}
				}else{
					res_IS_neg_screen_cal[[i]]<-numeric(0)
				}
			}
			names(res_IS_neg_screen_cal)<-names(IS_neg_screen_listed_cal)
		}
	
		##################################################################################################################
		# save list ######################################################################################################
		names(res_IS_neg_screen_cal)<-names(pattern)
		save(res_IS_neg_screen_cal,file=file.path(logfile$project_folder,"quantification","res_IS_neg_screen_cal"))
		# assemble output table of length(list) ##########################################################################
		# iterator m is directly equal to the sample ID ##################################################################
		if( length(IS_neg_screen_listed_cal)>0 ){
			intstand<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
			results_screen_IS_neg_cal<-enviMass:::get_screening_results(
				screened_listed=res_IS_neg_screen_cal,
				pattern=pattern,
				at_RT=pattern_RT,
				profileList=profileList_neg_cal,
				measurements_table=measurements,
				compound_table=intstand,
				cut_score=cut_score
			)
			save(results_screen_IS_neg_cal,file=file.path(logfile$project_folder,"quantification","results_screen_IS_neg_cal"))
			rm(intstand,results_screen_IS_neg_cal);
		}
		##################################################################################################################
		rm(getit,IS_neg_screen_listed_cal,res_IS_neg_screen_cal)
		rm(pattern,pattern_RT,pattern_delRT,envir=as.environment(".GlobalEnv"))
		##################################################################################################################

		##################################################################################################################
		# TARGET SCREENING ###############################################################################################
		#load(file=file.path(as.character(logfile[[1]]),"quantification","profileList_neg_cal"),envir=as.environment(".GlobalEnv"));
		load(file=file.path(logfile[[1]],"results","pattern_neg_target"),envir=as.environment(".GlobalEnv"));
		pattern<<-pattern_neg_target;rm(pattern_neg_target,envir=as.environment(".GlobalEnv"));
		load(file=file.path(logfile[[1]],"results","patternRT_neg_target"),envir=as.environment(".GlobalEnv"));
		pattern_RT<<-patternRT_neg_target;rm(patternRT_neg_target,envir=as.environment(".GlobalEnv"));
		load(file=file.path(logfile[[1]],"results","patternDelRT_neg_target"),envir=as.environment(".GlobalEnv"));
		pattern_delRT<<-patternDelRT_neg_target;rm(patternDelRT_neg_target,envir=as.environment(".GlobalEnv"));
		targets<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");
		
		mztol<-as.numeric(logfile$parameters$tar_dmz)				# m/z tolerance ...
		ppm<-as.logical(as.character(logfile$parameters$tar_ppm))	# ... given in pppm?
		cutint<-as.numeric(logfile$parameters$tar_intcut)			# Lower intensity threhold
		int_tol<-as.numeric(logfile$parameters$tar_inttol)			# Intensity tolerance %
		#RT_tol_outside<-as.numeric(logfile$parameters$target_drt1)		# RT tolerance of peaks in sample relative to their expected RT [s] - incorporated via pattern_delRT during isotopologue pattern generation
		RT_tol_inside<-as.numeric(logfile$parameters$tar_drt2)		# RT tolerance of peaks within an isotope pattern [s]
		cut_score<-as.numeric(logfile$parameters$tar_w1)	
		
		peaks<-profileList_neg_cal[[7]];
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
			forID<-strsplit(names(pattern)[i],"_")[[1]][1]
			if(targets[targets[,1]==forID,6]!="FALSE"){ # only include compounds used for quantification
				centro_mass[at_ID:(at_ID+n-1)]<-pattern[[i]][,1]
				centro_maxpeak[at_ID:(at_ID+n-1)]<-(pattern[[i]][,2]==max(pattern[[i]][,2]))			
				centro_RT[at_ID:(at_ID+n-1)]<-pattern_RT[i]
				centro_dRT[at_ID:(at_ID+n-1)]<-pattern_delRT[i]
			}
			screen_list[[i]]<-as.list(rep("FALSE",n))
			centro_ID[at_ID:(at_ID+n-1)]<-i
			centro_number[at_ID:(at_ID+n-1)]<-(1:n)
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
		if( as.character(logfile$parameters$screen_target_maxonly)=="TRUE" ){ # only retain max_peak-results?
			getit[!centro_maxpeak]<-"FALSE"
		}
		for(i in 1:length(getit)){ # transfer to a fist list of compoundadduct x centroids
			screen_list[[centro_ID[i]]][[centro_number[i]]]<-getit[i]
		}
		# resort to a full result list: pattern x sample x (centroids,matches) ( = peak index in profileList_neg_cal)		
		target_neg_screen_listed_cal<-list()  # default: no match at all	
		set_ID<-seq(1:length(profileList_neg_cal[[4]]))	
		for(i in 1:length(screen_list)){ # over compound x adduct = i
			if(any(is.na(screen_list[[i]]==FALSE))){
				target_neg_screen_listed_cal[[i]]<-list() # m-level		
				for(j in 1:length(screen_list[[i]])){ # over its centroids = j
					if(screen_list[[i]][[j]]!="FALSE"){ 
						profs<-as.numeric(strsplit(screen_list[[i]][[j]]," / ")[[1]])
						for(k in 1:length(profs)){ # over their matched profile peaks = k		
							if(profileList_neg_cal[[7]][profs[k],4]!=profs[k]){cat("\n debug me: profile ID mismatch");stop();} # a check
							for(m in profileList_neg_cal[[7]][profs[k],1]:profileList_neg_cal[[7]][profs[k],2]){ # over their sample peaks			
								delmass<-abs(profileList_neg_cal[[2]][m,1]-pattern[[i]][j,1])		
								if(!ppm){
									if(delmass>mztol){next}
								}else{
									if((delmass*1E6/pattern[[i]][j,1])>mztol){next}
								}
								at_ID<-set_ID[profileList_neg_cal[[4]]==as.character(profileList_neg_cal[[2]][m,6])]								
								if(length(target_neg_screen_listed_cal[[i]])<at_ID){							
									target_neg_screen_listed_cal[[i]][[at_ID]]<-matrix(ncol=2,nrow=0)	
								}else{
									if(length(target_neg_screen_listed_cal[[i]][[at_ID]])==0){
										target_neg_screen_listed_cal[[i]][[at_ID]]<-matrix(ncol=2,nrow=0)	
									}
								}
								target_neg_screen_listed_cal[[i]][[at_ID]]<-rbind(
									target_neg_screen_listed_cal[[i]][[at_ID]],c(j,m)
								)					
								colnames(target_neg_screen_listed_cal[[i]][[at_ID]])<-c(as.character(profileList_neg_cal[[4]][at_ID]),"")								
							}							
						}
					}
				}
			}else{
				target_neg_screen_listed_cal[[i]]<-numeric(0)	
			}
		}	
		# decompose ###########################################################################		
		if( logfile$parameters$screen_target_cutit=="TRUE" ){
			use_score_cut<-TRUE;
			score_cut<-cut_score
		}else{
			use_score_cut<-FALSE;
			score_cut<-0		
		}
		many<-0
		many_unamb<-0
		res_target_neg_screen_cal<-list()  # default: no match at all
		if(length(target_neg_screen_listed_cal)>0){
			for(i in 1:length(target_neg_screen_listed_cal)){ # i - on compound_adduct
				if(length(target_neg_screen_listed_cal[[i]])>0){
					res_target_neg_screen_cal[[i]]<-list()
					for(m in 1:length(target_neg_screen_listed_cal[[i]])){ # m (relates to IDs in profileList_neg_cal[[4]])					
						at_ID<-set_ID[profileList_neg_cal[[4]]==colnames(target_neg_screen_listed_cal[[i]][[m]])[1]]					
						if(length(target_neg_screen_listed_cal[[i]][[m]])>0){
							if(do_LOD){								
								with_model<-which(names(LOD_splined)==paste("LOD_",colnames(target_neg_screen_listed_cal[[i]][[m]])[1],sep=""))						
								if(length(with_model)>0){						
									use_cutint<-10^(predict(LOD_splined[[with_model]],pattern_RT[i])$y)
								}else{
									cat("\n Missing LOD model; using default intensity threshold. Debug?")
									use_cutint<-cutint;
								}
							}else{
								use_cutint<-cutint
							}						
							combination_matches<-recomb_score(
								cent_peak_mat=target_neg_screen_listed_cal[[i]][[m]],
								pattern_compound=pattern[[i]],
								profileList=profileList_neg_cal,
								LOD=use_cutint,
								RT_tol_inside=RT_tol_inside,
								int_tol=int_tol,
								use_score_cut=use_score_cut,
								score_cut=score_cut,
								plotit=FALSE,
								verbose=FALSE,
								RT_seperate=TRUE
							)
							for(k in 1:length(combination_matches)){ # add file ID
								combination_matches[[k]][[10]]<-colnames(target_neg_screen_listed_cal[[i]][[m]])[1]
								names(combination_matches[[k]])[10]<-"file_ID"
							}
							res_target_neg_screen_cal[[i]][[at_ID]]<-combination_matches
							if(length(combination_matches)>1){many_unamb<-(many_unamb+1)}
							many<-(many+1)
						}
					}
				}else{
					res_target_neg_screen_cal[[i]]<-numeric(0)
				}
			}
			names(res_target_neg_screen_cal)<-names(target_neg_screen_listed_cal)
		}
		##################################################################################################################
		# save list ######################################################################################################
		names(res_target_neg_screen_cal)<-names(pattern)
		save(res_target_neg_screen_cal,file=file.path(logfile$project_folder,"quantification","res_target_neg_screen_cal"))
		# assemble output table of length(list) ##########################################################################
		# iterator m is directly equal to the sample ID ##################################################################
		if( length(target_neg_screen_listed_cal)>0 ){
			targets<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");
			results_screen_target_neg_cal<-enviMass:::get_screening_results(
				screened_listed=res_target_neg_screen_cal,
				pattern=pattern,
				at_RT=pattern_RT,
				profileList=profileList_neg_cal,
				measurements_table=measurements,
				compound_table=targets,
				cut_score=cut_score
			)
			save(results_screen_target_neg_cal,file=file.path(logfile$project_folder,"quantification","results_screen_target_neg_cal"))
			rm(targets,results_screen_target_neg_cal);
		}
		##################################################################################################################
		rm(getit,target_neg_screen_listed_cal,res_target_neg_screen_cal)
		rm(pattern,pattern_RT,pattern_delRT,envir=as.environment(".GlobalEnv"))
		##################################################################################################################

	} # on negative calibration




