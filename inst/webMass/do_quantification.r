
	
		
	# POSITIVE IONIZATION ##################################################################
	# check if any calibration files are available at all ##################################
	all_files<-list.files(file.path(logfile$project_folder,"quantification"))
	got_models_pos<-FALSE
	got_models_neg<-FALSE	
	if(any(grepl("cal_models_pos_",all_files))){got_models_pos<-TRUE}
	if(any(grepl("cal_models_neg_",all_files))){got_models_neg<-TRUE}
	
	if(
		(file.exists(file=file.path(logfile$project_folder,"results","screening","results_screen_target_pos"))) &	# the summary screening table for targets
		(file.exists(file=file.path(logfile$project_folder,"results","screening","results_screen_IS_pos"))) &		# the summary screening table for IS compounds	
		got_models_pos
	){
if(FALSE){	


		# LOAD DATA ########################################################################
		load(file=file.path(logfile$project_folder,"results","screening","res_target_pos_screen"))		
		load(file=file.path(logfile$project_folder,"results","screening","results_screen_target_pos")) # requires new entries to [[1]] and [[2]]
		load(file=file.path(logfile$project_folder,"results","screening","res_IS_pos_screen"))
		target_table<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");	
		IS_table<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
		measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		# WHICH measurements belongs to which calibration set, if at all? ##################
		measurements<-measurements[measurements[,4]=="positive",,drop=FALSE]
		latest_ID<-get_latestID(measurements)
		cal_files<-measurements[measurements[,3]=="calibration",,drop=FALSE]
		cal_files<-unique(cal_files[,c(20,6,7,22,23),drop=FALSE])
		starttime<-as.difftime(cal_files[,3]);startdate<-as.Date(cal_files[,2]);
		numstart<-(as.numeric(startdate)+as.numeric(starttime/24))		
		endtime<-as.difftime(cal_files[,5]);enddate<-as.Date(cal_files[,4]);
		numend<-(as.numeric(enddate)+as.numeric(endtime/24))		
		use_files<-measurements[measurements[,3]!="calibration",,drop=FALSE]
		use_group<-c()
		for(i in 1:length(use_files[,1])){ # determine which file belongs to which calibration group
			attime<-as.difftime(use_files[i,7]);atdate<-as.Date(use_files[i,6]);
			numuse<-(as.numeric(atdate)+as.numeric(attime/24))		
			if(any((numuse>=numstart) & (numuse<=numend))){
				use_group<-c(use_group,
					cal_files[(numuse>=numstart) & (numuse<=numend),1]
				)
			}else{
				use_group<-c(use_group,"FALSE")
			}
		}		
		# load available calibration models into a single list - check availability ########
		cal_models_pos_used<-list()
		use_group_names<-unique(use_group)
		use_group_names<-use_group_names[use_group_names!="FALSE"]
		if(length(use_group_names)>0){	
			cat("\nLoading calibration models ...")
			for(i in 1:length(use_group_names)){
				if(file.exists(file.path(logfile[[1]],"quantification",paste("cal_models_pos",use_group_names[i],sep="_")))){
					load(file=file.path(logfile[[1]],"quantification",paste("cal_models_pos",use_group_names[i],sep="_")))
					at<-(length(cal_models_pos_used)+1)
					cal_models_pos_used[[at]]<-cal_models_pos[[1]]
					names(cal_models_pos_used)[at]<-names(cal_models_pos)
					rm(cal_models_pos)
				}
			}
			cat(" done.")
		}
		# QUANTIFY #########################################################################
		system.time({
		if(length(cal_models_pos_used)>0){ # no calibration models? 
			res_IS_names<-rep("",length(res_IS_pos_screen))
			res_IS_adduct<-rep("",length(res_IS_pos_screen))
			for(i in 1:length(res_IS_pos_screen)){
				res_IS_names[i]<-strsplit(names(res_IS_pos_screen)[i],"_")[[1]][1]
				res_IS_adduct[i]<-strsplit(names(res_IS_pos_screen)[i],"_")[[1]][2]
			}
			found_which<-list() # save indices to write faster into summary table
			if(length(res_target_pos_screen)>0){
				for(i in 1:length(res_target_pos_screen)){
					if(length(res_target_pos_screen[[i]])>0){ # anything screened?
						at_ID<-strsplit(names(res_target_pos_screen)[i],"_")[[1]][1]
						if(target_table[target_table[,1]==at_ID,6]!="FALSE"){ # target linked to an IS compound?
							at_adduct_target<-strsplit(names(res_target_pos_screen)[i],"_")[[1]][2]
							if(target_table[target_table[,1]==at_ID,20]!=at_adduct_target){next} # relevant quantification adduct?
							at_IS<-target_table[target_table[,1]==at_ID,6]
							at_peak_target<-as.numeric(target_table[target_table[,1]==at_ID,21])
							for(j in 1:length(res_target_pos_screen[[i]])){ # over files
								if(length(res_target_pos_screen[[i]][[j]])>0){ # at this file, target
									at_sample<-res_target_pos_screen[[i]][[j]][[1]]$file_ID
									# check: does a calibration model exist?		
									at_group<-which(names(cal_models_pos_used)==use_group[use_files[,1]==at_sample])
									if(length(at_group)==0){next} # calibration group available?
									at_group_model<-which(names(cal_models_pos_used[[at_group]])==paste(at_IS,at_ID,sep="_"))					
									if(length(at_group_model)==0){next} # model in that group available?
									at_adduct_IS<-IS_table[IS_table[,1]==at_IS,19] 	# get relevant IS adduct
									at_peak_IS<-IS_table[IS_table[,1]==at_IS,20] 	# get relevant IS peak
									get_conc<-c()
									for(k in 1:length(res_target_pos_screen[[i]][[j]])){ # all matches for this file
										use_target_peak<-which(res_target_pos_screen[[i]][[j]][[k]]$Peaks[,1]==at_peak_target)
										if(length(use_target_peak)>0){
											# find IS results for this target
											at_IS_entry<-which((res_IS_names==at_IS)&(res_IS_adduct==at_adduct_IS))
											at_IS_entry_sample<-which(names(res_IS_pos_screen[[at_IS_entry]])==at_sample)
											if(length(at_IS_entry_sample)==0){next}
											# get IS intensity bounds
											low_bound<-as.numeric(IS_table[IS_table[,1]==at_IS,17])
											if(low_bound!=0){low_bound<-(10^low_bound)}
											high_bound<-as.numeric(IS_table[IS_table[,1]==at_IS,18])
											if(high_bound!=0){high_bound<-(10^high_bound)}
											for(w in 1:length(res_IS_pos_screen[[at_IS_entry]][[at_IS_entry_sample]])){
												# check IS peak
												use_IS_peak<-which(res_IS_pos_screen[[at_IS_entry]][[at_IS_entry_sample]][[w]]$Peaks[,1]==at_peak_IS)
												if(length(use_IS_peak)==0){next}
												if(res_IS_pos_screen[[at_IS_entry]][[at_IS_entry_sample]][[w]]$Peaks[use_IS_peak,2]>high_bound){next} # out of intensity bounds?
												if(res_IS_pos_screen[[at_IS_entry]][[at_IS_entry_sample]][[w]]$Peaks[use_IS_peak,2]<low_bound){next}  # out of intensity bounds?
												int_IS<-(res_IS_pos_screen[[at_IS_entry]][[at_IS_entry_sample]][[w]]$Peaks[use_IS_peak,2])
												int_target<-(res_target_pos_screen[[i]][[j]][[k]]$Peaks[use_target_peak,2])
												int_rat<-(int_target/int_IS)
												new_conc<-c()
												if(cal_models_pos_used[[at_group]][[at_group_model]]$call=="resp ~ 0 + lin"){ # linear, 0-intercept
													new_conc<-(
														cal_models_pos_used[[at_group]][[at_group_model]]$coefficients[[1]]*int_rat
													)
												}
												if(cal_models_pos_used[[at_group]][[at_group_model]]$call=="resp ~ lin"){ # linear, with intercept
													new_conc<-(
														cal_models_pos_used[[at_group]][[at_group_model]]$coefficients[[1]]+(cal_models_pos_used[[at_group]][[at_group_model]]$coefficients[[2]]*int_rat)
													)
												}
												if(cal_models_pos_used[[at_group]][[at_group_model]]$call=="resp ~ 0 + lin + quad"){ # quadratic, 0-intercept
													new_conc<-(
														(cal_models_pos_used[[at_group]][[at_group_model]]$coefficients[[1]]*int_rat)+(cal_models_pos_used[[at_group]][[at_group_model]]$coefficients[[2]]*(int_rat^2))
													)
												}
												if(cal_models_pos_used[[at_group]][[at_group_model]]$call=="resp ~ lin + quad"){ # quadratic, 0-intercept
													new_conc<-(
														(cal_models_pos_used[[at_group]][[at_group_model]]$coefficients[[1]])+
														(cal_models_pos_used[[at_group]][[at_group_model]]$coefficients[[2]]*int_rat)+
														(cal_models_pos_used[[at_group]][[at_group_model]]$coefficients[[3]]*(int_rat^2))														
													)
												}
												cat(".")							
												get_conc<-c(get_conc,new_conc)				
											}
											if(length(get_conc)==0){next}
											res_target_pos_screen[[i]][[j]][[k]]$conc<-get_conc
											found_which[[length(found_which)+1]]<-c(i,j,k)
											# i = target/addduct: names(res_target_pos_screen)[i]
											# j = file: res_target_pos_screen[[i]][[j]][[k]]$file_ID
											# k = all matches: res_target_pos_screen[[i]][[j]][[k]]
										}	
									}
								}						
							}
						}
					}		
				}
			}
		}	
		# MAKE ENTRY INTO SUMMARY TABLE ####################################################
		if(length(found_which)>0){
			for(m in 1:length(found_which)){
				at_ID<-strsplit(names(res_target_pos_screen)[found_which[[m]][1]],"_")[[1]][1]
				at_adduct<-strsplit(names(res_target_pos_screen)[found_which[[m]][1]],"_")[[1]][2]
				for(n in 1:length(res_target_pos_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$conc)){ # over the different concentration values found per target
					# on results_screen_target_pos[[1]] - each adduct itemized #############
					at_entry<-which((results_screen_target_pos[[1]][,1]==at_ID) & (results_screen_target_pos[[1]][,3]==at_adduct))
					if( is.na(results_screen_target_pos[[1]][at_entry,11]) ){
						results_screen_target_pos[[1]][at_entry,11]<-round(
							res_target_pos_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$conc[n],
							digits=2)
					}else{
						if(
							results_screen_target_pos[[1]][at_entry,11]<
							res_target_pos_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$conc[n]
						){
							results_screen_target_pos[[1]][at_entry,11]<-
							res_target_pos_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$conc[n]					
						}
					}	
					if(res_target_pos_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$file_ID==latest_ID){ # update latest conc.
						if( is.na(results_screen_target_pos[[1]][at_entry,12]) ){
							results_screen_target_pos[[1]][at_entry,12]<-round(
								res_target_pos_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$conc[n],
								digits=2)
						}else{
							if(
								results_screen_target_pos[[1]][at_entry,12]<
								res_target_pos_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$conc[n]
							){
								results_screen_target_pos[[1]][at_entry,12]<-
								res_target_pos_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$conc[n]					
							}
						}	
					}
					# on results_screen_target_pos[[2]] - adducts per compound summarized #	
					at_entry<-which((results_screen_target_pos[[2]][,1]==at_ID))
					if( is.na(results_screen_target_pos[[2]][at_entry,8]) ){
						results_screen_target_pos[[2]][at_entry,8]<-round(
							res_target_pos_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$conc[n],
							digits=2)
					}else{
						if(
							results_screen_target_pos[[2]][at_entry,8]<
							res_target_pos_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$conc[n]
						){
							results_screen_target_pos[[2]][at_entry,8]<-
							res_target_pos_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$conc[n]					
						}
					}	
					if(res_target_pos_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$file_ID==latest_ID){ # update latest conc.
						if( is.na(results_screen_target_pos[[2]][at_entry,9]) ){
							results_screen_target_pos[[2]][at_entry,9]<-round(
								res_target_pos_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$conc[n],
								digits=2)
						}else{
							if(
								results_screen_target_pos[[2]][at_entry,9]<
								res_target_pos_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$conc[n]
							){
								results_screen_target_pos[[2]][at_entry,9]<-
								res_target_pos_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$conc[n]					
							}
						}	
					}
				}	
			}
		}
		})
		# INSERT & SAVE RESULTS ############################################################
		save(res_target_pos_screen,file=file.path(logfile$project_folder,"results","screening","res_target_pos_screen"))
		save(results_screen_target_pos,file=file.path(logfile$project_folder,"results","screening","results_screen_target_pos"))


	} # if FALSE
	
}	
	
	
	# NEGATIVE IONIZATION ##################################################################
	

