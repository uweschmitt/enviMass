
	# 
	
	# POSITIVE IONIZATION ##################################################################
	if(
		(file.exists(file=file.path(logfile$project_folder,"results","screening","results_screen_target_pos"))) &	# the summary screening table for targets
		(file.exists(file=file.path(logfile$project_folder,"results","screening","results_screen_IS_pos"))) &		# the summary screening table for IS compounds	
		(file.exists(file=file.path(logfile$project_folder,"quantification","cal_models_pos")))		
	){
	
		# LOAD DATA ########################################################################
		load(file=file.path(logfile$project_folder,"results","screening","res_target_pos_screen"))		
		load(file=file.path(logfile$project_folder,"results","screening","results_screen_target_pos")) # requires new entries to [[1]] and [[2]]
		load(file=file.path(logfile$project_folder,"results","screening","res_IS_pos_screen"))
		load(file=file.path(logfile$project_folder,"quantification","cal_models_pos"))
		target_table<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");	
		IS_table<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
		measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		# WHICH measurements belongs to which calibration set, if at all? ##################
		measurements<-measurements[measurements[,4]=="positive",,drop=FALSE]
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
				use_group<-c(use_group,"FALSE"
				)
			}
		}		
		# QUANTIFY #########################################################################
		#system.time({
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
						if(target_table[target_table[,1]==at_ID,20]==at_adduct_target){ # relevant quantification adduct?
							at_IS<-target_table[target_table[,1]==at_ID,6]
							at_peak_target<-as.numeric(target_table[target_table[,1]==at_ID,21])
							for(j in 1:length(res_target_pos_screen[[i]])){ # over files
								if(length(res_target_pos_screen[[i]][[j]])>0){ # at this file, target
									at_sample<-res_target_pos_screen[[i]][[j]][[1]]$file_ID
									# check: does a calibration model exist?		
									at_group<-which(names(cal_models_pos)==use_group[use_files[,1]==at_sample])
									if(length(at_group)==0){next} # calibration group available?
									at_group_model<-which(names(cal_models_pos[[at_group]])==paste(at_IS,at_ID,sep="_"))					
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
											high_bound<-as.numeric(IS_table[IS_table[,1]==at_IS,18])
											for(w in 1:length(res_IS_pos_screen[[at_IS_entry]][[at_IS_entry_sample]])){
												# check IS peak
												use_IS_peak<-which(res_IS_pos_screen[[at_IS_entry]][[at_IS_entry_sample]][[w]]$Peaks[,1]==at_peak_IS)
												if(length(use_IS_peak)==0){next}
												if(res_IS_pos_screen[[at_IS_entry]][[at_IS_entry_sample]][[w]]$Peaks[use_IS_peak,2]>high_bound){next} # out of intensity bounds?
												if(res_IS_pos_screen[[at_IS_entry]][[at_IS_entry_sample]][[w]]$Peaks[use_IS_peak,2]<low_bound){next}  # out of intensity bounds?
												int_IS<-(res_IS_pos_screen[[at_IS_entry]][[at_IS_entry_sample]][[w]]$Peaks[use_IS_peak,2])
												int_target<-(res_target_pos_screen[[i]][[j]][[k]]$Peaks[use_target_peak,2])
												int_rat<-(int_target/int_IS)
												if(length(cal_models_pos[[at_group]][[at_group_model]][[1]])==2){ # for the linear model
													new_conc<-predict(cal_models_pos[[at_group]][[at_group_model]],list(lin=int_rat))
												}																				
												if(length(cal_models_pos[[at_group]][[at_group_model]][[1]])==3){ # for the quadratic model
													int_rat2<-(int_rat^2)
													new_conc<-predict(cal_models_pos[[at_group]][[at_group_model]],list(lin=int_rat,quad=int_rat2))
												}
												cat(".")							
												get_conc<-c(get_conc,new_conc)#stop("HALTED")
											}
											res_target_pos_screen[[i]][[j]][[k]]$conc<-get_conc
											found_which[[length(found_which)+1]]<-c(i,j,k)
										}	
									}
								}						
							}
						}
					}
				}		
			}
		}
		#})	
		# MAKE ENTRY INTO SUMMARY TABLE ####################################################
		if(length(found_which)>0){
		
# BAUSTELLE		
		
		
		
		}
		# INSERT & SAVE RESULTS ############################################################
		save(res_target_pos_screen,file=file.path(logfile$project_folder,"results","screening","res_target_pos_screen"))
		save(results_screen_target_pos,file=file.path(logfile$project_folder,"results","screening","results_screen_target_pos"))
	}
	
	
	
	
	# NEGATIVE IONIZATION ##################################################################
	

