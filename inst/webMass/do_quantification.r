
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
		measurements<-measurements[measurements[,3]=="calibration",,drop=FALSE]
		cal_files<-measurements[measurements[,4]=="positive",,drop=FALSE]
		cal_files<-unique(cal_files[,c(20,22,23),drop=FALSE])
		
		
		
# add to checkproject: (1) dates for calibration files in a group all equal / (2) no overlaps in dates among calibration groups
		
system.time({		
		# QUANTIFY #########################################################################
		if(length(res_target_pos_screen)>0){
			for(i in 1:length(res_target_pos_screen)){
				if(length(res_target_pos_screen[[i]])>0){ # anything screened?
					at_ID<-strsplit(names(res_target_pos_screen)[i],"_")[[1]][1]
					if(target_table[target_table[,1]==at_ID,6]!="FALSE"){ # target linked to an IS compound?
						at_adduct_target<-strsplit(names(res_target_pos_screen)[i],"_")[[1]][2]
						if(target_table[target_table[,1]==at_ID,20]==at_adduct_target){ # relevant quantification adduct?
							at_IS<-target_table[target_table[,1]==at_ID,6]
							at_peak_target<-as.numeric(target_table[target_table[,1]==at_ID,21])
							for(j in 1:length(res_target_pos_screen[[i]])){
								if(length(res_target_pos_screen[[i]][[j]])>0){ # at this sample
									at_sample<-res_target_pos_screen[[i]][[j]][[1]]$file_ID
# check: does a calibration model exist?								
								
									at_adduct_IS<-IS_table[IS_table[,1]==at_IS,19] 	# get relevant IS adduct
									at_peak_IS<-IS_table[IS_table[,1]==at_IS,20] 	# get relevant IS peak
									for(k in 1:length(res_target_pos_screen[[i]][[j]])){ # all matches for this file
										if(any(res_target_pos_screen[[i]][[j]][[k]]$Peaks[,1]==at_peak_target)){
										
											#stop()					
										
										
										
										}	
									}
								}						
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
	}
	
	
	
	
	# NEGATIVE IONIZATION ##################################################################
	

