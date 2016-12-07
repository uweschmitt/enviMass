
	########################################################################################
	# remove SOME old results ##############################################################
	those<-list(0)
	those[[1]]<-file.path(logfile[[1]],"quantification","target_quant_table_pos")		
	those[[2]]<-file.path(logfile[[1]],"quantification","target_quant_table_pos_warn")				
	those[[3]]<-file.path(logfile[[1]],"quantification","target_quant_table_neg")			
	those[[4]]<-file.path(logfile[[1]],"quantification","target_quant_table_neg_warn")	
	for(n in 1:length(those)){
		if(file.exists(those[[n]])){
			file.remove(those[[n]])
		}
	}	
	rm(those)
	########################################################################################
	
	########################################################################################
	# check if any calibration files are available at all ##################################
	all_files<-list.files(file.path(logfile$project_folder,"quantification"))
	got_models_pos<-FALSE
	got_models_neg<-FALSE	
	if(any(grepl("cal_models_pos_",all_files))){got_models_pos<-TRUE}
	if(any(grepl("cal_models_neg_",all_files))){got_models_neg<-TRUE}	
	########################################################################################

	
	# POSITIVE IONIZATION ##################################################################	
	if(
		(file.exists(file=file.path(logfile$project_folder,"results","screening","results_screen_target_pos"))) &	# the summary screening table for targets
		(file.exists(file=file.path(logfile$project_folder,"results","screening","results_screen_IS_pos"))) &		# the summary screening table for IS compounds	
		got_models_pos
	){	

		# LOAD DATA ########################################################################
		load(file=file.path(logfile$project_folder,"results","screening","res_target_pos_screen"))		
		load(file=file.path(logfile$project_folder,"results","screening","results_screen_target_pos")) # requires to fill entries in [[1]] and [[2]]
		load(file=file.path(logfile$project_folder,"results","screening","res_IS_pos_screen"))
		target_table<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");	
		IS_table<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
		measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");		
		# WHICH measurements belongs to which calibration set, if at all? ##################
		measurements<-measurements[measurements[,"Mode"]=="positive",,drop=FALSE]
		latest_ID<-get_latestID(measurements)
		cal_files<-measurements[measurements[,"Type"]=="calibration",,drop=FALSE]
		cal_files<-unique(cal_files[,c(20,6,7,22,23),drop=FALSE])
		starttime<-as.difftime(cal_files[,3]);startdate<-as.Date(cal_files[,2]);
		numstart<-(as.numeric(startdate)+as.numeric(starttime/24))		
		endtime<-as.difftime(cal_files[,5]);enddate<-as.Date(cal_files[,4]);
		numend<-(as.numeric(enddate)+as.numeric(endtime/24))		
		use_files<-measurements[(measurements[,"Type"]!="calibration"),,drop=FALSE] # MUST retain spiked files!
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
					source(file=file.path(logfile[[1]],"quantification",paste("cal_models_pos",use_group_names[i],sep="_")),local=as.environment(".GlobalEnv"))
					at<-(length(cal_models_pos_used)+1)
					cal_models_pos_used[[at]]<-cal_models_pos[[1]]
					names(cal_models_pos_used)[at]<-names(cal_models_pos)
					rm(cal_models_pos,envir=as.environment(".GlobalEnv"))
				}
			}
			cat(" done.\n")
		}
		# REMOVE PREVIOUS RESULTS & RE-INITIATE ############################################
		if(file.exists(file.path(logfile[[1]],"quantification","target_quant_table_pos"))){
			file.remove(file.path(logfile[[1]],"quantification","target_quant_table_pos"))
		}
		those_files<-use_files[use_group!="FALSE",,drop=FALSE] # only files covered by a calibration model are considered!
		if(length(those_files[,1])==0){
			shinyjs:::info("No samples to be quantified (positive mode)? Check if timing of samples are covered by calibration file sets periods; otherwise consider removing the quantification step from your workflow!")
			stop("\nNo samples to be quantified (positive mode)? Check if timing of samples are covered by calibration file sets periods; otherwise consider removing the quantification step from your workflow!")
		}
		# sort those files by decreasing date!
		atdate<-those_files[,6]
		atdate<-as.Date(atdate);
		attime<-those_files[,7]
		attime<-as.difftime(attime);
		ord<-order(as.numeric(atdate),as.numeric(attime),as.numeric(those_files[,1]),decreasing=TRUE);
		those_files<-those_files[ord,,drop=FALSE]	
		if(logfile$parameters$quant_files_included!="FALSE"){ # restrict number of files to inlcude
			if(logfile$parameters$quant_files_included<length(those_files[,1])){
				those_files<-those_files[1:logfile$parameters$quant_files_included,,drop=FALSE]		
			}
		}
		those_targets<-target_table[target_table[,"ID_internal_standard"]!="FALSE",,drop=FALSE]
		those_targets<-those_targets[those_targets[,"ion_mode"]=="positive",,drop=FALSE]
		target_quant_table_pos<-matrix(nrow=(length(those_targets[,1])+5),ncol=(length(those_files[,1])+2),"")
		colnames(target_quant_table_pos)<-c("Target ID","Target name",those_files[,1])
		rownames(target_quant_table_pos)<-c("Name","Type","Date","Time","Custom ID",rep("",length(those_targets[,1])))
		target_quant_table_pos[1,]<-c("","",as.character(those_files[,"Name"]))
		target_quant_table_pos[2,]<-c("","",as.character(those_files[,"Type"]))
		target_quant_table_pos[3,]<-c("","",as.character(those_files[,"Date"]))
		target_quant_table_pos[4,]<-c("","",as.character(those_files[,"Time"]))	
		target_quant_table_pos[5,]<-c("","",as.character(those_files[,"ID_2"]))			
		target_quant_table_pos[,1]<-c("","","","","",those_targets[,1])
		target_quant_table_pos[,2]<-c("","","","","",those_targets[,2])
		target_quant_table_pos_warn<-target_quant_table_pos
		target_quant_table_pos_warn[6:length(target_quant_table_pos_warn[,1]),3:length(target_quant_table_pos_warn[1,])]<-"0"
		# QUANTIFY #########################################################################
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
					at_adduct_target<-strsplit(names(res_target_pos_screen)[i],"_")[[1]][2]
					at_ID<-strsplit(names(res_target_pos_screen)[i],"_")[[1]][1]
					if(target_table[target_table[,"ID"]==at_ID,"Quant_adduct"]!=at_adduct_target){next} # relevant quantification adduct?
					at_1<-which(target_quant_table_pos[,"Target ID"]==at_ID)
					if(length(at_1)==0){next} # in case: target not in table
					if(length(res_target_pos_screen[[i]])>0){ # anything screened in any file?
						if(target_table[target_table[,"ID"]==at_ID,"ID_internal_standard"]!="FALSE"){ # target linked to an IS compound?
							at_IS<-target_table[target_table[,"ID"]==at_ID,"ID_internal_standard"]
							at_peak_target<-as.numeric(target_table[target_table[,"ID"]==at_ID,"Quant_peak"])
							rule_target<-target_table[target_table[,"ID"]==at_ID,"Quant_rule"] # "most intense peak","closest RT","closest m/z"
							rule_IS<-IS_table[IS_table[,"ID"]==at_IS,"Quant_rule"]		
							RT_target<-as.numeric(target_table[target_table[,"ID"]==at_ID,"RT"])
							RT_IS<-as.numeric(IS_table[IS_table[,"ID"]==at_IS,"RT"])
							at_adduct_IS<-IS_table[IS_table[,"ID"]==at_IS,"Quant_adduct"] 	# get relevant IS adduct
							at_peak_IS<-IS_table[IS_table[,"ID"]==at_IS,"Quant_peak"] 	# get relevant IS peak
							at_IS_entry<-which((res_IS_names==at_IS)&(res_IS_adduct==at_adduct_IS))
							low_bound<-as.numeric(IS_table[IS_table[,"ID"]==at_IS,"Lower_intensity_bound"])
							if(low_bound!=0){low_bound<-(10^low_bound)}
							high_bound<-as.numeric(IS_table[IS_table[,"ID"]==at_IS,"Upper_intensity_bound"])
							if(high_bound!=0){high_bound<-(10^high_bound)}								
							for(j in 1:length(res_target_pos_screen[[i]])){ # over files
								at_sample<-res_target_pos_screen[[i]][[j]][[1]]$file_ID
								at_IS_entry_sample<-which(names(res_IS_pos_screen[[at_IS_entry]])==at_sample)
								if(length(at_IS_entry_sample)==0){next} # "No ISTD screening matches for this file."
								at_2<-which(colnames(target_quant_table_pos)==at_sample)
								if(length(at_2)==0){next} # file not in table 
								if(length(res_target_pos_screen[[i]][[j]])>0){ # at this file, target
									# check: does a calibration model exist?		
									at_group<-which(names(cal_models_pos_used)==use_group[use_files[,1]==at_sample])
									if(length(at_group)==0){ next } # calibration group available?						
									at_group_model<-which(names(cal_models_pos_used[[at_group]])==paste("_",at_IS,"_",at_ID,"_",sep=""))					
									if(length(at_group_model)==0){ next }; # model in that group available?										
									# collect intensities for target ########################### 
									int_target<-c()
									int_target_rank<-c()
									at_k<-c()
									for(k in 1:length(res_target_pos_screen[[i]][[j]])){ # all matches for this file
										use_target_peak<-which(res_target_pos_screen[[i]][[j]][[k]]$Peaks[,1]==at_peak_target)
										if(length(use_target_peak)>0){									
											int_target<-c(int_target,
												(res_target_pos_screen[[i]][[j]][[k]]$Intensity[use_target_peak])
											)
											if(rule_target=="most intense peak"){
												int_target_rank<-c(int_target_rank,
													(res_target_pos_screen[[i]][[j]][[k]]$Intensity[use_target_peak])
												)
											}
											if(rule_target=="closest RT"){
												int_target_rank<-c(int_target_rank,
													1/abs(RT_target-(res_target_pos_screen[[i]][[j]][[k]]$RT[use_target_peak]))
												)											
											}											
											if(rule_target=="closest m/z"){
												int_target_rank<-c(int_target_rank,
													1/abs(res_target_pos_screen[[i]][[j]][[k]]$ppm[use_target_peak])
												)											
											}		
											at_k<-c(at_k,k)
										}
									}	
									if(length(int_target)==0){next}
									# collect intensities for ISTD ##############################
									int_IS<-c()
									int_IS_rank<-c()
									for(k in 1:length(res_IS_pos_screen[[at_IS_entry]][[at_IS_entry_sample]]) ){ # all matches for this file
										use_IS_peak<-which(res_IS_pos_screen[[at_IS_entry]][[at_IS_entry_sample]][[k]]$Peaks[,1]==at_peak_IS)
										if(length(use_IS_peak)>0){									
											int_IS<-c(int_IS,
												(res_IS_pos_screen[[at_IS_entry]][[at_IS_entry_sample]][[k]]$Intensity[use_IS_peak])
											)
											if(rule_IS=="most intense peak"){
												int_IS_rank<-c(int_IS_rank,
													(res_IS_pos_screen[[at_IS_entry]][[at_IS_entry_sample]][[k]]$Intensity[use_IS_peak])
												)
											}
											if(rule_IS=="closest RT"){
												int_IS_rank<-c(int_IS_rank,
													1/abs(RT_IS-(res_IS_pos_screen[[at_IS_entry]][[at_IS_entry_sample]][[k]]$RT[use_IS_peak]))
												)											
											}											
											if(rule_IS=="closest m/z"){
												int_IS_rank<-c(int_IS_rank,
													1/abs(res_IS_pos_screen[[at_IS_entry]][[at_IS_entry_sample]][[k]]$ppm[use_IS_peak])
												)											
											}													
										}
									}	
									if(length(int_IS)==0){next}
									# re-rank intensities #########################################
									int_target<-(int_target[order(int_target_rank,decreasing=TRUE)])
									int_IS<-(int_IS[order(int_IS_rank,decreasing=TRUE)])
									# calculate concentrations ####################################
									get_conc<-c()
									for_k<-c()
									for(a in 1:length(int_target)){
										for(b in 1:length(int_IS)){
											for_k<-c(for_k,at_k[a])
											int_rat<-(int_target[a]/int_IS[b])
											if(int_rat<cal_models_pos_used[[at_group]][[at_group_model]]$low_bound){cat(",");next}		
											if(int_rat>cal_models_pos_used[[at_group]][[at_group_model]]$up_bound){cat(",");next}											
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
									}
									# make entries for further below usage in summary table
									if(length(get_conc)==0){next}
									for(a in 1:length(get_conc)){
											res_target_pos_screen[[i]][[j]][[for_k[a]]]$conc<-c(
												res_target_pos_screen[[i]][[j]][[for_k[a]]]$conc,
												get_conc[a]
											)
											found_which[[length(found_which)+1]]<-c(i,j,for_k[a])	
											# i = target/addduct: names(res_target_pos_screen)[i]
											# j = file: res_target_pos_screen[[i]][[j]][[k]]$file_ID
											# k = all matches: res_target_pos_screen[[i]][[j]][[k]]											
									}																	
									if(length(get_conc)==0){next}
									#if(length(get_conc)>1){stop()}
									get_conc<-round(get_conc,digits=1)
									get_conc<-unique(get_conc)							
									at_3<-which(those_targets$ID==at_ID)
									target_quant_table_pos_warn[at_1,at_2]<-"1"
									warn_conc_1<-those_targets$warn_1[at_3]
									if(warn_conc_1!="FALSE"){
										warn_conc_1<-as.numeric(warn_conc_1)
										#if(any(get_conc>=warn_conc_1)){
										if(any(get_conc[1]>=warn_conc_1)){
											target_quant_table_pos_warn[at_1,at_2]<-"2"																				
										}
									}
									warn_conc_2<-those_targets$warn_2[at_3]
									if(warn_conc_2!="FALSE"){
										warn_conc_2<-as.numeric(warn_conc_2)
										#if(any(get_conc>=warn_conc_2)){
										if(any(get_conc[1]>=warn_conc_2)){
											target_quant_table_pos_warn[at_1,at_2]<-"3"
										}
									}														
									get_conc<-as.character(get_conc)
									get_conc<-paste(get_conc,collapse=", ")
									target_quant_table_pos[at_1,at_2]<-get_conc								
								}						
							}						
						}
					}		
				}
			}
		}	
		# sort target_quant_table_pos ######################################################	target_quant_table_pos[1:20,1:10]
		if(length(target_quant_table_pos[,1])>5){
			this<-which(target_quant_table_pos[2,]=="sample")[1]
			if(is.na(this)){this<-3}
			splitted<-strsplit(target_quant_table_pos[-c(1:5),this],",")
			splitted<-lapply(splitted,"as.numeric")
			for(i in 1:length(splitted)){
				if(identical(splitted[[i]],numeric(0))){
					splitted[[i]]<-0
				}else{
					splitted[[i]]<-max(splitted[[i]])
				}
			}
			splitted<-unlist(splitted)
			ord<-order(splitted,decreasing=TRUE)
			target_quant_table_pos[6:length(target_quant_table_pos[,1]),]<-
				(target_quant_table_pos[6:length(target_quant_table_pos[,1]),,drop=FALSE][ord,,drop=FALSE])
			target_quant_table_pos_warn[6:length(target_quant_table_pos_warn[,1]),]<-
				(target_quant_table_pos_warn[6:length(target_quant_table_pos_warn[,1]),,drop=FALSE][ord,,drop=FALSE])		
		}
		# state why missing quantification arose ###########################################
		if(length(target_quant_table_pos[,1])>5){		
			for(i in 6:length(target_quant_table_pos[,1])){ # over targets
				if(any(target_quant_table_pos[i,]=="")){
					# info target:
					at_ID_tar<-(target_quant_table_pos[i,"Target ID"])
					at_adduct_tar<-target_table[target_table[,"ID"]==at_ID_tar,"Quant_adduct"]
					at_peak_tar<-as.numeric(target_table[target_table[,"ID"]==at_ID_tar,"Quant_peak"])
					find_name<-paste(at_ID_tar,"_",at_adduct_tar,"_",sep="")
					at_res_tar<-which(substr(names(res_target_pos_screen),1,nchar(find_name))==find_name)
					# info IS:					
					at_ID_IS<-target_table[target_table[,"ID"]==at_ID_tar,"ID_internal_standard"]					
					at_adduct_IS<-IS_table[IS_table[,"ID"]==at_ID_IS,"Quant_adduct"]
					at_peak_IS<-IS_table[IS_table[,"ID"]==at_ID_IS,"Quant_peak"]
					find_name<-paste(at_ID_IS,"_",at_adduct_IS,"_",sep="")
					at_res_IS<-which(substr(names(res_IS_pos_screen),1,nchar(find_name))==find_name)
					if(length(at_res_tar)==0){stop("Debug in do_quantification.r required!")}
					for(j in 3:length(target_quant_table_pos[1,])){ # over samples
						if(target_quant_table_pos[i,j]!=""){next}
						reason<-"!"
						# calibration models cover file (just a check, should always be the case as such files are excluded further above)
						at_sample<-colnames(target_quant_table_pos)[j]
						at_group<-which(names(cal_models_pos_used)==use_group[use_files[,1]==at_sample])
						if(length(at_group)==0){ 
							reason<-paste(reason,"no quantific.",sep=" / ")
						} 		
						# if above correct, then: calibration models for this target exist at all ... ?									
						if(!any(grepl(paste("_",at_ID_tar,"_",sep=""),names(cal_models_pos_used[[at_group]])))){ 
							reason<-paste(reason,"missing calibration model",sep=" / ")
						# ... or just not the one with the ISTD specified in the target table?
						}else{
							at_group_model<-which(names(cal_models_pos_used[[at_group]])==paste("_",at_ID_IS,"_",at_ID_tar,"_",sep=""))
							if(length(at_group_model)==0){
								reason<-paste(reason,"only model w/ incorrect ISTD",sep=" / ")
							}
						}				
						# target screening results exist, correct peak found?
						at_sam<-which(names(res_target_pos_screen[[at_res_tar]])==at_sample)
						if(length(at_sam)==0){
							reason<-paste(reason,"no target matches",sep=" / ")
						}else{
							if(length(res_target_pos_screen[[at_res_tar]][[at_sam]])==0){
								reason<-paste(reason,"no target matches",sep=" / ")
							}else{
								found_it<-FALSE
								for(k in 1:length(res_target_pos_screen[[at_res_tar]][[at_sam]])){
									if(any(res_target_pos_screen[[at_res_tar]][[at_sam]][[k]]$Peaks[,1]==at_peak_tar)){
										found_it<-TRUE;
										break;
									}
								}
								if(!found_it){
									reason<-paste(reason,"no target match for quant. peak",sep=" / ")
								}
							}
						}
						# IS screening results exist, correct peak found?
						at_sam<-which(names(res_IS_pos_screen[[at_res_IS]])==at_sample)
						if(length(at_sam)==0){
							reason<-paste(reason,"no ISTD matches",sep=" / ")
						}else{
							if(length(res_IS_pos_screen[[at_res_IS]][[at_sam]])==0){
								reason<-paste(reason,"no ISTD matches at all",sep=" / ")
							}else{
								found_it<-FALSE
								for(k in 1:length(res_IS_pos_screen[[at_res_IS]][[at_sam]])){
									if(any(res_IS_pos_screen[[at_res_IS]][[at_sam]][[k]]$Peaks[,1]==at_peak_IS)){
										found_it<-TRUE;
										break;
									}
								}
								if(!found_it){				
									reason<-paste(reason,"no ISTD match for quant. peak",sep=" / ")
								}
							}
						}					
						if(reason=="!"){
							reason<-paste(reason,"ratio possibly out of bounds",sep=" / ")
						}
						target_quant_table_pos[i,j]<-reason
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
		# INSERT & SAVE RESULTS ############################################################
		save(res_target_pos_screen,file=file.path(logfile$project_folder,"results","screening","res_target_pos_screen"))
		save(results_screen_target_pos,file=file.path(logfile$project_folder,"results","screening","results_screen_target_pos"))
		save(target_quant_table_pos,file=file.path(logfile[[1]],"quantification","target_quant_table_pos"))
		save(target_quant_table_pos_warn,file=file.path(logfile[[1]],"quantification","target_quant_table_pos_warn"))
	
}	
	

	
	# NEGATIVE IONIZATION ##################################################################	
	if(
		(file.exists(file=file.path(logfile$project_folder,"results","screening","results_screen_target_neg"))) &	# the summary screening table for targets
		(file.exists(file=file.path(logfile$project_folder,"results","screening","results_screen_IS_neg"))) &		# the summary screening table for IS compounds	
		got_models_neg
	){	

		# LOAD DATA ########################################################################
		load(file=file.path(logfile$project_folder,"results","screening","res_target_neg_screen"))		
		load(file=file.path(logfile$project_folder,"results","screening","results_screen_target_neg")) # requires to fill entries in [[1]] and [[2]]
		load(file=file.path(logfile$project_folder,"results","screening","res_IS_neg_screen"))
		target_table<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");	
		IS_table<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
		measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");		
		# WHICH measurements belongs to which calibration set, if at all? ##################
		measurements<-measurements[measurements[,"Mode"]=="negative",,drop=FALSE]
		latest_ID<-get_latestID(measurements)
		cal_files<-measurements[measurements[,"Type"]=="calibration",,drop=FALSE]
		cal_files<-unique(cal_files[,c(20,6,7,22,23),drop=FALSE])
		starttime<-as.difftime(cal_files[,3]);startdate<-as.Date(cal_files[,2]);
		numstart<-(as.numeric(startdate)+as.numeric(starttime/24))		
		endtime<-as.difftime(cal_files[,5]);enddate<-as.Date(cal_files[,4]);
		numend<-(as.numeric(enddate)+as.numeric(endtime/24))		
		use_files<-measurements[(measurements[,"Type"]!="calibration"),,drop=FALSE] # MUST retain spiked files!
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
		cal_models_neg_used<-list()
		use_group_names<-unique(use_group)
		use_group_names<-use_group_names[use_group_names!="FALSE"]
		if(length(use_group_names)>0){	
			cat("\nLoading calibration models ...")
			for(i in 1:length(use_group_names)){
				if(file.exists(file.path(logfile[[1]],"quantification",paste("cal_models_neg",use_group_names[i],sep="_")))){
					source(file=file.path(logfile[[1]],"quantification",paste("cal_models_neg",use_group_names[i],sep="_")),local=as.environment(".GlobalEnv"))
					at<-(length(cal_models_neg_used)+1)
					cal_models_neg_used[[at]]<-cal_models_neg[[1]]
					names(cal_models_neg_used)[at]<-names(cal_models_neg)
					rm(cal_models_neg,envir=as.environment(".GlobalEnv"))
				}
			}
			cat(" done.\n")
		}
		# REMOVE PREVIOUS RESULTS & RE-INITIATE ############################################
		if(file.exists(file.path(logfile[[1]],"quantification","target_quant_table_neg"))){
			file.remove(file.path(logfile[[1]],"quantification","target_quant_table_neg"))
		}
		those_files<-use_files[use_group!="FALSE",,drop=FALSE] # only files covered by a calibration model are considered!
		if(length(those_files[,1])==0){
			shinyjs:::info("No samples to be quantified (negative mode)? Check if timing of samples are covered by calibration file sets periods; otherwise consider removing the quantification step from your workflow!")
			stop("\nNo samples to be quantified (negative mode)? Check if timing of samples are covered by calibration file sets periods; otherwise consider removing the quantification step from your workflow!")
		}
		# sort those files by decreasing date!
		atdate<-those_files[,6]
		atdate<-as.Date(atdate);
		attime<-those_files[,7]
		attime<-as.difftime(attime);
		ord<-order(as.numeric(atdate),as.numeric(attime),as.numeric(those_files[,1]),decreasing=TRUE);
		those_files<-those_files[ord,,drop=FALSE]	
		if(logfile$parameters$quant_files_included!="FALSE"){ # restrict number of files to include
			if(logfile$parameters$quant_files_included<length(those_files[,1])){
				those_files<-those_files[1:logfile$parameters$quant_files_included,,drop=FALSE]		
			}
		}
		those_targets<-target_table[target_table[,"ID_internal_standard"]!="FALSE",,drop=FALSE]
		those_targets<-those_targets[those_targets[,"ion_mode"]=="negative",,drop=FALSE]
		target_quant_table_neg<-matrix(nrow=(length(those_targets[,1])+5),ncol=(length(those_files[,1])+2),"")
		colnames(target_quant_table_neg)<-c("Target ID","Target name",those_files[,1])
		rownames(target_quant_table_neg)<-c("Name","Type","Date","Time","Custom ID",rep("",length(those_targets[,1])))
		target_quant_table_neg[1,]<-c("","",as.character(those_files[,"Name"]))
		target_quant_table_neg[2,]<-c("","",as.character(those_files[,"Type"]))
		target_quant_table_neg[3,]<-c("","",as.character(those_files[,"Date"]))
		target_quant_table_neg[4,]<-c("","",as.character(those_files[,"Time"]))	
		target_quant_table_neg[5,]<-c("","",as.character(those_files[,"ID_2"]))				
		target_quant_table_neg[,1]<-c("","","","","",those_targets[,1])
		target_quant_table_neg[,2]<-c("","","","","",those_targets[,2])
		target_quant_table_neg_warn<-target_quant_table_neg
		target_quant_table_neg_warn[6:length(target_quant_table_neg_warn[,1]),3:length(target_quant_table_neg_warn[1,])]<-"0"
		# QUANTIFY #########################################################################
		if(length(cal_models_neg_used)>0){ # no calibration models? 
			res_IS_names<-rep("",length(res_IS_neg_screen))
			res_IS_adduct<-rep("",length(res_IS_neg_screen))
			for(i in 1:length(res_IS_neg_screen)){
				res_IS_names[i]<-strsplit(names(res_IS_neg_screen)[i],"_")[[1]][1]
				res_IS_adduct[i]<-strsplit(names(res_IS_neg_screen)[i],"_")[[1]][2]
			}
			found_which<-list() # save indices to write faster into summary table
			if(length(res_target_neg_screen)>0){
				for(i in 1:length(res_target_neg_screen)){
					at_adduct_target<-strsplit(names(res_target_neg_screen)[i],"_")[[1]][2]
					at_ID<-strsplit(names(res_target_neg_screen)[i],"_")[[1]][1]
					if(target_table[target_table[,"ID"]==at_ID,"Quant_adduct"]!=at_adduct_target){next} # relevant quantification adduct?
					at_1<-which(target_quant_table_neg[,"Target ID"]==at_ID)
					if(length(at_1)==0){next} # in case: target not in table
					if(length(res_target_neg_screen[[i]])>0){ # anything screened in any file?
						if(target_table[target_table[,"ID"]==at_ID,"ID_internal_standard"]!="FALSE"){ # target linked to an IS compound?
							at_IS<-target_table[target_table[,"ID"]==at_ID,"ID_internal_standard"]
							at_peak_target<-as.numeric(target_table[target_table[,"ID"]==at_ID,"Quant_peak"])
							rule_target<-target_table[target_table[,"ID"]==at_ID,"Quant_rule"] # "most intense peak","closest RT","closest m/z"
							rule_IS<-IS_table[IS_table[,"ID"]==at_IS,"Quant_rule"]		
							RT_target<-as.numeric(target_table[target_table[,"ID"]==at_ID,"RT"])
							RT_IS<-as.numeric(IS_table[IS_table[,"ID"]==at_IS,"RT"])
							at_adduct_IS<-IS_table[IS_table[,"ID"]==at_IS,"Quant_adduct"] 	# get relevant IS adduct
							at_peak_IS<-IS_table[IS_table[,"ID"]==at_IS,"Quant_peak"] 	# get relevant IS peak
							at_IS_entry<-which((res_IS_names==at_IS)&(res_IS_adduct==at_adduct_IS))
							low_bound<-as.numeric(IS_table[IS_table[,"ID"]==at_IS,"Lower_intensity_bound"])
							if(low_bound!=0){low_bound<-(10^low_bound)}
							high_bound<-as.numeric(IS_table[IS_table[,"ID"]==at_IS,"Upper_intensity_bound"])
							if(high_bound!=0){high_bound<-(10^high_bound)}								
							for(j in 1:length(res_target_neg_screen[[i]])){ # over files
								at_sample<-res_target_neg_screen[[i]][[j]][[1]]$file_ID
								at_IS_entry_sample<-which(names(res_IS_neg_screen[[at_IS_entry]])==at_sample)
								if(length(at_IS_entry_sample)==0){next} # "No ISTD screening matches for this file."
								at_2<-which(colnames(target_quant_table_neg)==at_sample)
								if(length(at_2)==0){next} # file not in table 
								if(length(res_target_neg_screen[[i]][[j]])>0){ # at this file, target
									# check: does a calibration model exist?		
									at_group<-which(names(cal_models_neg_used)==use_group[use_files[,1]==at_sample])
									if(length(at_group)==0){ next } # calibration group available?						
									at_group_model<-which(names(cal_models_neg_used[[at_group]])==paste("_",at_IS,"_",at_ID,"_",sep=""))					
									if(length(at_group_model)==0){ next }; # model in that group available?										
									# collect intensities for target ########################### 
									int_target<-c()
									int_target_rank<-c()
									at_k<-c()
									for(k in 1:length(res_target_neg_screen[[i]][[j]])){ # all matches for this file
										use_target_peak<-which(res_target_neg_screen[[i]][[j]][[k]]$Peaks[,1]==at_peak_target)
										if(length(use_target_peak)>0){									
											int_target<-c(int_target,
												(res_target_neg_screen[[i]][[j]][[k]]$Intensity[use_target_peak])
											)
											if(rule_target=="most intense peak"){
												int_target_rank<-c(int_target_rank,
													(res_target_neg_screen[[i]][[j]][[k]]$Intensity[use_target_peak])
												)
											}
											if(rule_target=="closest RT"){
												int_target_rank<-c(int_target_rank,
													1/abs(RT_target-(res_target_neg_screen[[i]][[j]][[k]]$RT[use_target_peak]))
												)											
											}											
											if(rule_target=="closest m/z"){
												int_target_rank<-c(int_target_rank,
													1/abs(res_target_neg_screen[[i]][[j]][[k]]$ppm[use_target_peak])
												)											
											}		
											at_k<-c(at_k,k)
										}
									}	
									if(length(int_target)==0){next}
									# collect intensities for ISTD ##############################
									int_IS<-c()
									int_IS_rank<-c()
									for(k in 1:length(res_IS_neg_screen[[at_IS_entry]][[at_IS_entry_sample]]) ){ # all matches for this file
										use_IS_peak<-which(res_IS_neg_screen[[at_IS_entry]][[at_IS_entry_sample]][[k]]$Peaks[,1]==at_peak_IS)
										if(length(use_IS_peak)>0){									
											int_IS<-c(int_IS,
												(res_IS_neg_screen[[at_IS_entry]][[at_IS_entry_sample]][[k]]$Intensity[use_IS_peak])
											)
											if(rule_IS=="most intense peak"){
												int_IS_rank<-c(int_IS_rank,
													(res_IS_neg_screen[[at_IS_entry]][[at_IS_entry_sample]][[k]]$Intensity[use_IS_peak])
												)
											}
											if(rule_IS=="closest RT"){
												int_IS_rank<-c(int_IS_rank,
													1/abs(RT_IS-(res_IS_neg_screen[[at_IS_entry]][[at_IS_entry_sample]][[k]]$RT[use_IS_peak]))
												)											
											}											
											if(rule_IS=="closest m/z"){
												int_IS_rank<-c(int_IS_rank,
													1/abs(res_IS_neg_screen[[at_IS_entry]][[at_IS_entry_sample]][[k]]$ppm[use_IS_peak])
												)											
											}													
										}
									}	
									if(length(int_IS)==0){next}
									# re-rank intensities #########################################
									int_target<-(int_target[order(int_target_rank,decreasing=TRUE)])
									int_IS<-(int_IS[order(int_IS_rank,decreasing=TRUE)])
									# calculate concentrations ####################################
									get_conc<-c()
									for_k<-c()
									for(a in 1:length(int_target)){
										for(b in 1:length(int_IS)){
											for_k<-c(for_k,at_k[a])
											int_rat<-(int_target[a]/int_IS[b])
											if(int_rat<cal_models_neg_used[[at_group]][[at_group_model]]$low_bound){cat(",");next}		
											if(int_rat>cal_models_neg_used[[at_group]][[at_group_model]]$up_bound){cat(",");next}	
											if(cal_models_neg_used[[at_group]][[at_group_model]]$call=="resp ~ 0 + lin"){ # linear, 0-intercept
												new_conc<-(
													cal_models_neg_used[[at_group]][[at_group_model]]$coefficients[[1]]*int_rat
												)
											}
											if(cal_models_neg_used[[at_group]][[at_group_model]]$call=="resp ~ lin"){ # linear, with intercept
												new_conc<-(
													cal_models_neg_used[[at_group]][[at_group_model]]$coefficients[[1]]+(cal_models_neg_used[[at_group]][[at_group_model]]$coefficients[[2]]*int_rat)
												)
											}
											if(cal_models_neg_used[[at_group]][[at_group_model]]$call=="resp ~ 0 + lin + quad"){ # quadratic, 0-intercept
												new_conc<-(
													(cal_models_neg_used[[at_group]][[at_group_model]]$coefficients[[1]]*int_rat)+(cal_models_neg_used[[at_group]][[at_group_model]]$coefficients[[2]]*(int_rat^2))
												)
											}
											if(cal_models_neg_used[[at_group]][[at_group_model]]$call=="resp ~ lin + quad"){ # quadratic, 0-intercept
												new_conc<-(
													(cal_models_neg_used[[at_group]][[at_group_model]]$coefficients[[1]])+
													(cal_models_neg_used[[at_group]][[at_group_model]]$coefficients[[2]]*int_rat)+
													(cal_models_neg_used[[at_group]][[at_group_model]]$coefficients[[3]]*(int_rat^2))														
												)
											}
											cat(".")							
											get_conc<-c(get_conc,new_conc)				
										}
									}
									if(length(get_conc)==0){next}
									# make entries for further below uasge in summary table
									for(a in 1:length(get_conc)){
											res_target_neg_screen[[i]][[j]][[for_k[a]]]$conc<-c(
												res_target_neg_screen[[i]][[j]][[for_k[a]]]$conc,
												get_conc[a]
											)
											found_which[[length(found_which)+1]]<-c(i,j,for_k[a])	
											# i = target/addduct: names(res_target_neg_screen)[i]
											# j = file: res_target_neg_screen[[i]][[j]][[k]]$file_ID
											# k = all matches: res_target_neg_screen[[i]][[j]][[k]]											
									}																	
									if(length(get_conc)==0){next}
									#if(length(get_conc)>1){stop()}
									get_conc<-round(get_conc,digits=1)
									get_conc<-unique(get_conc)							
									at_3<-which(those_targets$ID==at_ID)
									target_quant_table_neg_warn[at_1,at_2]<-"1"
									warn_conc_1<-those_targets$warn_1[at_3]
									if(warn_conc_1!="FALSE"){
										warn_conc_1<-as.numeric(warn_conc_1)
										#if(any(get_conc>=warn_conc_1)){
										if(any(get_conc[1]>=warn_conc_1)){
											target_quant_table_neg_warn[at_1,at_2]<-"2"																				
										}
									}
									warn_conc_2<-those_targets$warn_2[at_3]
									if(warn_conc_2!="FALSE"){
										warn_conc_2<-as.numeric(warn_conc_2)
										#if(any(get_conc>=warn_conc_2)){
										if(any(get_conc[1]>=warn_conc_2)){
											target_quant_table_neg_warn[at_1,at_2]<-"3"
										}
									}														
									get_conc<-as.character(get_conc)
									get_conc<-paste(get_conc,collapse=", ")
									target_quant_table_neg[at_1,at_2]<-get_conc								
								}						
							}						
						}
					}		
				}
			}
		}	
		# sort target_quant_table_neg ######################################################
		if(length(target_quant_table_neg[,1])>5){
			this<-which(target_quant_table_neg[2,]=="sample")[1]
			if(is.na(this)){this<-3}
			splitted<-strsplit(target_quant_table_neg[-c(1:5),this],",")
			splitted<-lapply(splitted,"as.numeric")
			for(i in 1:length(splitted)){
				if(identical(splitted[[i]],numeric(0))){
					splitted[[i]]<-0
				}else{
					splitted[[i]]<-max(splitted[[i]])
				}
			}
			splitted<-unlist(splitted)
			ord<-order(splitted,decreasing=TRUE)
			target_quant_table_neg[6:length(target_quant_table_neg[,1]),]<-(target_quant_table_neg[6:length(target_quant_table_neg[,1]),,drop=FALSE][ord,,drop=FALSE])
			target_quant_table_neg_warn[6:length(target_quant_table_neg_warn[,1]),]<-(target_quant_table_neg_warn[6:length(target_quant_table_neg_warn[,1]),,drop=FALSE][ord,,drop=FALSE])
		}
		# state why missing quantification arose ###########################################
		if(length(target_quant_table_neg[,1])>5){		
			for(i in 6:length(target_quant_table_neg[,1])){ # over targets
				if(any(target_quant_table_neg[i,]=="")){
					# info target:
					at_ID_tar<-(target_quant_table_neg[i,"Target ID"])
					at_adduct_tar<-target_table[target_table[,"ID"]==at_ID_tar,"Quant_adduct"]
					at_peak_tar<-as.numeric(target_table[target_table[,"ID"]==at_ID_tar,"Quant_peak"])
					find_name<-paste(at_ID_tar,"_",at_adduct_tar,"_",sep="")
					at_res_tar<-which(substr(names(res_target_neg_screen),1,nchar(find_name))==find_name)
					# info IS:					
					at_ID_IS<-target_table[target_table[,"ID"]==at_ID_tar,"ID_internal_standard"]					
					at_adduct_IS<-IS_table[IS_table[,"ID"]==at_ID_IS,"Quant_adduct"]
					at_peak_IS<-IS_table[IS_table[,"ID"]==at_ID_IS,"Quant_peak"]
					find_name<-paste(at_ID_IS,"_",at_adduct_IS,"_",sep="")
					at_res_IS<-which(substr(names(res_IS_neg_screen),1,nchar(find_name))==find_name)
					if(length(at_res_tar)==0){stop("Debug in do_quantification.r required!")}
					for(j in 3:length(target_quant_table_neg[1,])){ # over samples
						if(target_quant_table_neg[i,j]!=""){next}
						reason<-"!"
						# calibration models cover file (just a check, should always be the case as such files are excluded further above)
						at_sample<-colnames(target_quant_table_neg)[j]
						at_group<-which(names(cal_models_neg_used)==use_group[use_files[,1]==at_sample])
						if(length(at_group)==0){ 
							reason<-paste(reason,"no quantific.",sep=" / ")
						} 		
						# if above correct, then: calibration models for this target exist at all ... ?					
						if(!any(grepl(paste("_",at_ID_tar,"_",sep=""),names(cal_models_neg_used[[at_group]])))){ 
							reason<-paste(reason,"missing calibration model",sep=" / ")
						# ... or just not the one with the ISTD specified in the target table?
						}else{
							at_group_model<-which(names(cal_models_neg_used[[at_group]])==paste("_",at_ID_IS,"_",at_ID_tar,"_",sep=""))
							if(length(at_group_model)==0){
								reason<-paste(reason,"incorrect calibration model",sep=" / ")
							}
						}
						# target screening results exist, correct peak found?
						at_sam<-which(names(res_target_neg_screen[[at_res_tar]])==at_sample)
						if(length(at_sam)==0){
							reason<-paste(reason,"no target matches",sep=" / ")
						}else{
							if(length(res_target_neg_screen[[at_res_tar]][[at_sam]])==0){
								reason<-paste(reason,"no target matches",sep=" / ")
							}else{
								found_it<-FALSE
								for(k in 1:length(res_target_neg_screen[[at_res_tar]][[at_sam]])){
									if(any(res_target_neg_screen[[at_res_tar]][[at_sam]][[k]]$Peaks[,1]==at_peak_tar)){
										found_it<-TRUE;
										break;
									}
								}
								if(!found_it){
									reason<-paste(reason,"no target match for quant. peak",sep=" / ")
								}
							}
						}
						# IS screening results exist, correct peak found?
						at_sam<-which(names(res_IS_neg_screen[[at_res_IS]])==at_sample)
						if(length(at_sam)==0){
							reason<-paste(reason,"no ISTD matches",sep=" / ")
						}else{
							if(length(res_IS_neg_screen[[at_res_IS]][[at_sam]])==0){
								reason<-paste(reason,"no ISTD matches at all",sep=" / ")
							}else{
								found_it<-FALSE
								for(k in 1:length(res_IS_neg_screen[[at_res_IS]][[at_sam]])){
									if(any(res_IS_neg_screen[[at_res_IS]][[at_sam]][[k]]$Peaks[,1]==at_peak_IS)){
										found_it<-TRUE;
										break;
									}
								}
								if(!found_it){				
									reason<-paste(reason,"no ISTD match for quant. peak",sep=" / ")
								}
							}
						}					
						if(reason=="!"){
							reason<-paste(reason,"ratio possibly out of bounds",sep=" / ")
						}
						target_quant_table_neg[i,j]<-reason
					}
				}
			}
		}
		# MAKE ENTRY INTO SUMMARY TABLE ####################################################
		if(length(found_which)>0){
			for(m in 1:length(found_which)){
				at_ID<-strsplit(names(res_target_neg_screen)[found_which[[m]][1]],"_")[[1]][1]
				at_adduct<-strsplit(names(res_target_neg_screen)[found_which[[m]][1]],"_")[[1]][2]
				for(n in 1:length(res_target_neg_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$conc)){ # over the different concentration values found per target
					# on results_screen_target_neg[[1]] - each adduct itemized #############
					at_entry<-which((results_screen_target_neg[[1]][,1]==at_ID) & (results_screen_target_neg[[1]][,3]==at_adduct))
					if( is.na(results_screen_target_neg[[1]][at_entry,11]) ){
						results_screen_target_neg[[1]][at_entry,11]<-round(
							res_target_neg_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$conc[n],
							digits=2)
					}else{
						if(
							results_screen_target_neg[[1]][at_entry,11]<
							res_target_neg_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$conc[n]
						){
							results_screen_target_neg[[1]][at_entry,11]<-
							res_target_neg_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$conc[n]					
						}
					}	
					if(res_target_neg_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$file_ID==latest_ID){ # update latest conc.
						if( is.na(results_screen_target_neg[[1]][at_entry,12]) ){
							results_screen_target_neg[[1]][at_entry,12]<-round(
								res_target_neg_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$conc[n],
								digits=2)
						}else{
							if(
								results_screen_target_neg[[1]][at_entry,12]<
								res_target_neg_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$conc[n]
							){
								results_screen_target_neg[[1]][at_entry,12]<-
								res_target_neg_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$conc[n]					
							}
						}	
					}
					# on results_screen_target_neg[[2]] - adducts per compound summarized #	
					at_entry<-which((results_screen_target_neg[[2]][,1]==at_ID))
					if( is.na(results_screen_target_neg[[2]][at_entry,8]) ){
						results_screen_target_neg[[2]][at_entry,8]<-round(
							res_target_neg_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$conc[n],
							digits=2)
					}else{
						if(
							results_screen_target_neg[[2]][at_entry,8]<
							res_target_neg_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$conc[n]
						){
							results_screen_target_neg[[2]][at_entry,8]<-
							res_target_neg_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$conc[n]					
						}
					}	
					if(res_target_neg_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$file_ID==latest_ID){ # update latest conc.
						if( is.na(results_screen_target_neg[[2]][at_entry,9]) ){
							results_screen_target_neg[[2]][at_entry,9]<-round(
								res_target_neg_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$conc[n],
								digits=2)
						}else{
							if(
								results_screen_target_neg[[2]][at_entry,9]<
								res_target_neg_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$conc[n]
							){
								results_screen_target_neg[[2]][at_entry,9]<-
								res_target_neg_screen [[ found_which[[m]][1] ]] [[ found_which[[m]][2] ]] [[ found_which[[m]][3] ]]$conc[n]					
							}
						}	
					}
				}	
			}
		}
		# INSERT & SAVE RESULTS ############################################################
		save(res_target_neg_screen,file=file.path(logfile$project_folder,"results","screening","res_target_neg_screen"))
		save(results_screen_target_neg,file=file.path(logfile$project_folder,"results","screening","results_screen_target_neg"))
		save(target_quant_table_neg,file=file.path(logfile[[1]],"quantification","target_quant_table_neg"))
		save(target_quant_table_neg_warn,file=file.path(logfile[[1]],"quantification","target_quant_table_neg_warn"))
	
}	
	


