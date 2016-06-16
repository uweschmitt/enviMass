if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_obs_screening.r!")}

###########################################################################################################
# SPECIFY IONIZATION MODE #################################################################################
observe({ 
	input$Ion_mode_Cal 
	init$b
	if(isolate(init$a)=="TRUE"){
		if(isolate(input$Ion_mode_Cal)=="positive"){
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			measurements<-measurements[measurements[,3]=="calibration",,drop=FALSE]
			measurements<-measurements[measurements[,4]=="positive",,drop=FALSE]
			if(length(measurements[,1])>0){
				those<-unique(measurements$tag2)
				if(all(those!="FALSE")){
					those<-c("none",those)
					updateSelectInput(session,"Cal_file_set","Specify calibration file group",choices = those, selected = those[1])
				}else{ # trigger warning
					cat("all calibration groups must have a tag2 other than FALSE!")
				}
			}
		}
		if(isolate(input$Ion_mode_Cal)=="negative"){
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			measurements<-measurements[measurements[,3]=="calibration",,drop=FALSE]
			measurements<-measurements[measurements[,4]=="negative",,drop=FALSE]
			if(length(measurements[,1])>0){
				those<-unique(measurements$tag2)
				if(all(those!="FALSE")){
					those<-c("none",those)
					updateSelectInput(session,"Cal_file_set","Specify calibration file group",choices = those, selected = those[1])
				}else{ # trigger warning
					cat("all calibration groups must have a tag2 other than FALSE!")
				}
			}
		}
	}
})
###########################################################################################################

###########################################################################################################
# SPECIFY CALIBRATION GROUP ###############################################################################
observe({ 
	input$Cal_file_set
	init$b
	if(isolate(init$a)=="TRUE"){
		if(isolate(input$Cal_file_set)!="none"){
			if(
				(isolate(input$Ion_mode_Cal)=="positive")&
				(file.exists(file.path(logfile[[1]],"quantification","results_screen_IS_pos_cal")))&
				(file.exists(file.path(logfile[[1]],"quantification","results_screen_target_pos_cal")))
			){
				load(file=file.path(logfile[[1]],"quantification","results_screen_IS_pos_cal"),envir=as.environment(".GlobalEnv"));	
				load(file=file.path(logfile[[1]],"quantification","profileList_pos_cal"),envir=as.environment(".GlobalEnv"));	
				load(file=file.path(logfile[[1]],"quantification","results_screen_target_pos_cal"),envir=as.environment(".GlobalEnv"));	
				load(file=file.path(logfile[[1]],"quantification","res_IS_pos_screen_cal"),envir=as.environment(".GlobalEnv"));	
				load(file=file.path(logfile[[1]],"quantification","res_target_pos_screen_cal"),envir=as.environment(".GlobalEnv"));	
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				measurements<-measurements[measurements[,3]=="calibration",,drop=FALSE]
				measurements<-measurements[measurements[,4]=="positive",,drop=FALSE]
				measurements<<-measurements
				targets<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");
				intstand<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
				targets<-targets[targets[,8]=="positive",,drop=FALSE]
				intstand<-intstand[intstand[,7]=="positive",,drop=FALSE]
				targets<-targets[targets[,6]!="FALSE",,drop=FALSE]
				targets<<-targets
				intstand<<-intstand
				# Update IS compounds
				IS_names<-unique(results_screen_IS_pos_cal[[1]][,2])
				IS_names<-IS_names[order(IS_names)]
				IS_names<-c("none",IS_names)
				updateSelectInput(session,inputId="Cal_IS_name",label="Name",choices=IS_names,selected = IS_names[1])
				IS_IDs<-unique(results_screen_IS_pos_cal[[1]][,1])		
				IS_IDs<-IS_IDs[order(IS_IDs)]
				IS_IDs<-c("none",IS_IDs)
				updateSelectInput(session,inputId="Cal_IS_ID",label="Internal standard ID",choices=IS_IDs,selected = IS_IDs[1])
				# Update target compounds
				target_names<-unique(results_screen_target_pos_cal[[1]][,2])
				target_names<-target_names[order(target_names)]
				target_names<-c("none",target_names)
				updateSelectInput(session,inputId="Cal_target_name",label="Name",choices=target_names,selected = target_names[1])
				target_IDs<-unique(results_screen_target_pos_cal[[1]][,1])		
				target_IDs<-target_IDs[order(target_IDs)]
				target_IDs<-c("none",target_IDs)
				updateSelectInput(session,inputId="Cal_target_ID",label="Target ID",choices=target_IDs,selected = target_IDs[1])
			}else{ # not available
				if((isolate(input$Ion_mode_Cal)!="negative")){
					updateSelectInput(session,inputId="Cal_IS_name",choices="none",selected = "none")
					updateSelectInput(session,inputId="Cal_IS_ID",choices="none",selected = "none")			
					updateSelectInput(session,inputId="Cal_target_name",choices="none",selected = "none")
					updateSelectInput(session,inputId="Cal_target_ID",choices="none",selected = "none")	
				}
			}
			if(
				(isolate(input$Ion_mode_Cal)=="negative")&
				(file.exists(file.path(logfile[[1]],"quantification","results_screen_IS_neg_cal")))&
				(file.exists(file.path(logfile[[1]],"quantification","results_screen_target_neg_cal")))			
			){
				load(file=file.path(logfile[[1]],"quantification","profileList_neg_cal"),envir=as.environment(".GlobalEnv"));	
				load(file=file.path(logfile[[1]],"quantification","results_screen_target_neg_cal"),envir=as.environment(".GlobalEnv"));	
				load(file=file.path(logfile[[1]],"quantification","results_screen_IS_neg_cal"),envir=as.environment(".GlobalEnv"));	
				load(file=file.path(logfile[[1]],"quantification","res_IS_neg_screen_cal"),envir=as.environment(".GlobalEnv"));	
				load(file=file.path(logfile[[1]],"quantification","res_target_neg_screen_cal"),envir=as.environment(".GlobalEnv"));	
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				measurements<-measurements[measurements[,3]=="calibration",,drop=FALSE]
				measurements<-measurements[measurements[,4]=="negative",,drop=FALSE]
				measurements<<-measurements
				targets<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");
				intstand<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
				targets<-targets[targets[,8]=="negative",,drop=FALSE]
				intstand<-intstand[intstand[,7]=="negative",,drop=FALSE]
				targets<-targets[targets[,6]!="FALSE",,drop=FALSE]
				targets<<-targets
				intstand<<-intstand
				# Update IS compounds
				IS_names<-unique(results_screen_IS_neg_cal[[1]][
					results_screen_IS_neg_cal[[1]][,4]!=0
				,2])
				IS_names<-IS_names[order(IS_names)]
				IS_names<-c("none",IS_names)
				updateSelectInput(session,inputId="Cal_IS_name",label="Name",choices=IS_names,selected = IS_names[1])
				IS_IDs<-unique(results_screen_IS_neg_cal[[1]][,1])		
				IS_IDs<-IS_IDs[order(IS_IDs)]
				IS_IDs<-c("none",IS_IDs)
				updateSelectInput(session,inputId="Cal_IS_ID",label="Internal standard ID",choices=IS_IDs,selected = IS_IDs[1])
				# Update target compounds
				target_names<-unique(results_screen_target_neg_cal[[1]][,2])
				target_names<-target_names[order(target_names)]
				target_names<-c("none",target_names)
				updateSelectInput(session,inputId="Cal_target_name",label="Name",choices=target_names,selected = target_names[1])
				target_IDs<-unique(results_screen_target_neg_cal[[1]][,1])		
				target_IDs<-target_IDs[order(target_IDs)]
				target_IDs<-c("none",target_IDs)
				updateSelectInput(session,inputId="Cal_target_ID",label="Target ID",choices=target_IDs,selected = target_IDs[1])
			}else{ # not available
				if((isolate(input$Ion_mode_Cal)!="positive")){
					updateSelectInput(session,inputId="Cal_IS_name",choices="none",selected = "none")
					updateSelectInput(session,inputId="Cal_IS_ID",choices="none",selected = "none")			
					updateSelectInput(session,inputId="Cal_target_name",choices="none",selected = "none")
					updateSelectInput(session,inputId="Cal_target_ID",choices="none",selected = "none")		
				}
			}
		}
	}
})
###########################################################################################################

###########################################################################################################
# MATCH COMPOUND NAMES AND ID SELECTIONS ##################################################################
observe({ 
	input$Cal_IS_ID
	init$b
	if((isolate(init$a)=="TRUE")){
		if(isolate(input$Cal_IS_ID)!="none"){
				if(isolate(input$Ion_mode_Cal)=="positive"){
					use_this_name<-unique(results_screen_IS_pos_cal[[1]][
						as.character(results_screen_IS_pos_cal[[1]][,1])==isolate(input$Cal_IS_ID)
					,2,drop=FALSE])
					updateSelectInput(session,inputId="Cal_IS_name",selected = as.character(use_this_name))
				}
				if(isolate(input$Ion_mode_Cal)=="negative"){
					use_this_name<-unique(results_screen_IS_neg_cal[[1]][
						as.character(results_screen_IS_neg_cal[[1]][,1])==isolate(input$Cal_IS_ID)
					,2,drop=FALSE])
					updateSelectInput(session,inputId="Cal_IS_name",selected = as.character(use_this_name))
				}
		}else{
			updateSelectInput(session,inputId="Cal_IS_name",selected="none")
		}
	}
})

observe({ 
	input$Cal_IS_name
	init$b
	if((isolate(init$a)=="TRUE")&(isolate(input$Cal_IS_name)!="none")){
			if(isolate(input$Ion_mode_Cal)=="positive"){
				use_this_ID<-unique(results_screen_IS_pos_cal[[1]][
					as.character(results_screen_IS_pos_cal[[1]][,2])==isolate(input$Cal_IS_name)
				,1,drop=FALSE])
				updateSelectInput(session,inputId="Cal_IS_ID",selected = as.character(use_this_ID))
			}
			if(isolate(input$Ion_mode_Cal)=="negative"){
				use_this_ID<-unique(results_screen_IS_neg_cal[[1]][
					as.character(results_screen_IS_neg_cal[[1]][,2])==isolate(input$Cal_IS_name)
				,1,drop=FALSE])
				updateSelectInput(session,inputId="Cal_IS_ID",selected = as.character(use_this_ID))
			}
	}	
})

observe({ 
	input$Cal_target_ID
	init$b
	if((isolate(init$a)=="TRUE")){
		if(isolate(input$Cal_target_ID)!="none"){
				if(isolate(input$Ion_mode_Cal)=="positive"){
					use_this_name<-unique(results_screen_target_pos_cal[[1]][
						as.character(results_screen_target_pos_cal[[1]][,1])==isolate(input$Cal_target_ID)
					,2,drop=FALSE])				
					updateSelectInput(session,inputId="Cal_target_name",selected = as.character(use_this_name))
				}
				if(isolate(input$Ion_mode_Cal)=="negative"){
					load(file=file.path(logfile[[1]],"quantification","results_screen_target_neg_cal"),envir=as.environment(".GlobalEnv"));	
					use_this_name<-unique(results_screen_target_neg_cal[[1]][
						as.character(results_screen_target_neg_cal[[1]][,1])==isolate(input$Cal_target_ID)
					,2,drop=FALSE])
					updateSelectInput(session,inputId="Cal_target_name",selected = as.character(use_this_name))
				}
		}else{	
			updateSelectInput(session,inputId="Cal_target_name",selected="none")
		}
	}	
})

observe({ 
	input$Cal_target_name
	init$b
	if((isolate(init$a)=="TRUE")&(isolate(input$Cal_target_name)!="none")){
			if(isolate(input$Ion_mode_Cal)=="positive"){
				use_this_ID<-unique(results_screen_target_pos_cal[[1]][
					as.character(results_screen_target_pos_cal[[1]][,2])==isolate(input$Cal_target_name)
				,1,drop=FALSE])
				updateSelectInput(session,inputId="Cal_target_ID",selected = as.character(use_this_ID))
			}
			if(isolate(input$Ion_mode_Cal)=="negative"){
				load(file=file.path(logfile[[1]],"quantification","results_screen_IS_neg_cal"),envir=as.environment(".GlobalEnv"));	
				use_this_ID<-unique(results_screen_target_neg_cal[[1]][
					as.character(results_screen_target_neg_cal[[1]][,2])==isolate(input$Cal_target_name)
				,1,drop=FALSE])
				updateSelectInput(session,inputId="Cal_target_ID",selected = as.character(use_this_ID))
			}
	}	
})
###########################################################################################################

###########################################################################################################
# FORWARD BACKWARD ETC ####################################################################################
observe({ 
	init$b
	input$Cal_next
	if((isolate(init$a)=="TRUE")){
		is_at_targetID<-(isolate(input$Cal_target_ID))
		is_at_ISID<-(isolate(input$Cal_IS_ID))
		in_table<-which( ((targets[,1]==is_at_targetID)&(targets[,6]==is_at_ISID)) )
		at_target_ID<-"none"
		at_IS_ID<-"none"
		if( (length(in_table)==0) || (is.na(in_table[1])) ){ # match with is_at_ISID not existing in table; reset to first available entry
			in_table<-which(targets[,1]==is_at_targetID)
			in_table<-in_table[1]
			if( (length(in_table)==0) || (is.na(in_table[1])) ){ # for any invalid entry, start from beginning
				in_table<-1
			}
			at_target_ID<-as.character(targets[in_table,1])
			at_IS_ID<-as.character(targets[in_table,6])	
		}else{ # match existing; get next entry
			if(in_table<length(targets[,1])){
				in_table<-(in_table+1)
				at_target_ID<-as.character(targets[in_table,1])
				at_IS_ID<-as.character(targets[in_table,6])			
			}else{
				at_target_ID<-"none"
				at_IS_ID<-"none"
			}
		}
		updateSelectInput(session, inputId="Cal_target_ID", selected = at_target_ID)
		updateSelectInput(session, inputId="Cal_IS_ID", selected = at_IS_ID)
	}

})

observe({ 
	init$b
	input$Cal_previous
	if((isolate(init$a)=="TRUE")&(isolate(input$Cal_target_name)!="none")){
		is_at_targetID<-(isolate(input$Cal_target_ID))
		is_at_ISID<-(isolate(input$Cal_IS_ID))
		in_table<-which( ((targets[,1]==is_at_targetID)&(targets[,6]==is_at_ISID)) )
		at_target_ID<-"none"
		at_IS_ID<-"none"
		if( (length(in_table)==0) || (is.na(in_table[1])) ){ # match with is_at_ISID not existing in table; reset to first available entry
			in_table<-which(targets[,1]==is_at_targetID)
			in_table<-in_table[1]
			if( (length(in_table)==0) || (is.na(in_table[1])) ){ # for any invalid entry, start from beginning
				in_table<-1
			}
			at_target_ID<-as.character(targets[in_table,1])
			at_IS_ID<-as.character(targets[in_table,6])	
		}else{ # match existing; get next entry
			if(in_table>1){
				in_table<-(in_table-1)
				at_target_ID<-as.character(targets[in_table,1])
				at_IS_ID<-as.character(targets[in_table,6])			
			}else{
				at_target_ID<-"none"
				at_IS_ID<-"none"
			}
		}	
		updateSelectInput(session, inputId="Cal_target_ID", selected = at_target_ID)
		updateSelectInput(session, inputId="Cal_IS_ID", selected = at_IS_ID)	
	}
})

observe({ 
	init$b
	input$Cal_first
	if((isolate(init$a)=="TRUE")){
		at_target_ID<-"none"
		at_IS_ID<-"none"
		if( (length(targets[,1])>0)){ # match with is_at_ISID not existing in table; reset to first available entry
			in_table<-1
			at_target_ID<-as.character(targets[in_table,1])
			at_IS_ID<-as.character(targets[in_table,6])	
		}
		updateSelectInput(session, inputId="Cal_target_ID", selected = at_target_ID)
		updateSelectInput(session, inputId="Cal_IS_ID", selected = at_IS_ID)
	}
})

observe({ 
	init$b
	input$Cal_last
	if((isolate(init$a)=="TRUE")){
		at_target_ID<-"none"
		at_IS_ID<-"none"
		if( (length(targets[,1])>0)){ # match with is_at_ISID not existing in table; reset to first available entry
			in_table<-length(targets[,1])
			at_target_ID<-as.character(targets[in_table,1])
			at_IS_ID<-as.character(targets[in_table,6])	
		}
		updateSelectInput(session, inputId="Cal_target_ID", selected = at_target_ID)
		updateSelectInput(session, inputId="Cal_IS_ID", selected = at_IS_ID)
	}
})

###########################################################################################################


###########################################################################################################
# RETRIEVE SETS & PLOT THEM ###############################################################################
observe({ 
	init$b
	input$Cal_IS_ID
	input$Cal_IS_name
	input$Cal_target_ID
	input$Cal_target_name
	if((isolate(init$a)=="TRUE")){	
			if(isolate(input$Ion_mode_Cal)=="positive"){	
				IS_ID<-isolate(input$Cal_IS_ID)
				target_ID<-isolate(input$Cal_target_ID)
				at_Cal<-isolate(input$Cal_file_set)
				#target_ID<-"4"; IS_ID<-"693"			
				# extract IS peaks ######################################################
				IS_adduct<-intstand[intstand[,1]==IS_ID,19]
				IS_peak<-as.numeric(intstand[intstand[,1]==IS_ID,20])
				at_entry<-FALSE
				for(i in 1:length(names(res_IS_pos_screen_cal))){ # where?
					if(
						(strsplit(names(res_IS_pos_screen_cal)[i],"_")[[1]][1]==IS_ID) &
						(strsplit(names(res_IS_pos_screen_cal)[i],"_")[[1]][2]==IS_adduct)
					){
						at_entry<-i;break;
					}
				}
				IS_in_file<-c()
				IS_with_peak<-c()
				IS_with_score<-c()
				if(length(res_IS_pos_screen_cal[[at_entry]])>0){
					for(j in 1:length(res_IS_pos_screen_cal[[at_entry]])){
						if(length(res_IS_pos_screen_cal[[at_entry]][[j]])>0){						
							if(measurements[measurements[,1]==res_IS_pos_screen_cal[[at_entry]][[j]][[1]]$file_ID,]$tag2==at_Cal){
								for(k in 1:length(res_IS_pos_screen_cal[[at_entry]][[j]])){
									if(any(res_IS_pos_screen_cal[[at_entry]][[j]][[k]]$Peaks[,1]==IS_peak)){
										that<-which(res_IS_pos_screen_cal[[at_entry]][[j]][[k]]$Peaks[,1]==IS_peak)
										IS_in_file<-c(IS_in_file,res_IS_pos_screen_cal[[at_entry]][[j]][[k]]$file_ID)
										IS_with_peak<-c(IS_with_peak,res_IS_pos_screen_cal[[at_entry]][[j]][[k]]$Peaks[that,2])
										IS_with_score<-c(IS_with_score,res_IS_pos_screen_cal[[at_entry]][[j]][[k]]$score_1)	
									}
								}
							}
						}
					}
				}
				# extract target peaks ##################################################
				target_adduct<-targets[targets[,1]==target_ID,20]
				target_peak<-as.numeric(targets[targets[,1]==target_ID,21])
				at_entry<-FALSE
				for(i in 1:length(names(res_target_pos_screen_cal))){ # where?
					if(
						(strsplit(names(res_target_pos_screen_cal)[i],"_")[[1]][1]==target_ID) &
						(strsplit(names(res_target_pos_screen_cal)[i],"_")[[1]][2]==target_adduct)
					){
						at_entry<-i;break;
					}
				}
				target_in_file<-c()
				target_with_peak<-c()
				target_with_score<-c()
				if(length(res_target_pos_screen_cal[[at_entry]])>0){
					for(j in 1:length(res_target_pos_screen_cal[[at_entry]])){
						if(length(res_target_pos_screen_cal[[at_entry]][[j]])>0){						
							if(measurements[measurements[,1]==res_target_pos_screen_cal[[at_entry]][[j]][[1]]$file_ID,]$tag2==at_Cal){
								for(k in 1:length(res_target_pos_screen_cal[[at_entry]][[j]])){
									if(any(res_target_pos_screen_cal[[at_entry]][[j]][[k]]$Peaks[,1]==target_peak)){
										that<-which(res_target_pos_screen_cal[[at_entry]][[j]][[k]]$Peaks[,1]==target_peak)
										target_in_file<-c(target_in_file,res_target_pos_screen_cal[[at_entry]][[j]][[k]]$file_ID)
										target_with_peak<-c(target_with_peak,res_target_pos_screen_cal[[at_entry]][[j]][[k]]$Peaks[that,2])
										target_with_score<-c(target_with_score,res_target_pos_screen_cal[[at_entry]][[j]][[k]]$score_1)
									}
								}
							}
						}
					}
				}
				# derive pairs ##########################################################
				mat_cal<-matrix(nrow=0,ncol=7)
				colnames(mat_cal)<-c("intensity_target","intensity_IS","ratio","concentration","score_target","score_IS","used?")
				if(	(length(target_in_file)>0) & (length(IS_in_file)>0)	){
					for(i in 1:length(target_in_file)){
						if(any(IS_in_file==target_in_file[i])){
							those<-which(IS_in_file==target_in_file[i])
							mat_cal<-rbind(mat_cal,
								cbind(
									rep(profileList_pos_cal[[2]][target_with_peak[i],2],length(those)), # intensity target
									(profileList_pos_cal[[2]][IS_with_peak[those],2]), 					# intensity IS
									((profileList_pos_cal[[2]][IS_with_peak[those],2])/profileList_pos_cal[[2]][target_with_peak[i],2]), # ration
									rep(as.numeric(
										measurements[measurements[,1]==target_in_file[i],]$tag1	
									),length(those)), # concentration
									rep(target_with_score[i],length(those)), # score target
									IS_with_score[those],
									rep(1,length(those)) # used?
								)
							)
						}
					}
				}
				print(mat_cal)
				# filter IS intensities #################################################
				
								
	
			}
			if(isolate(input$Ion_mode_Cal)=="negative"){	
# TO BE COMPLETED	
			}
	}	
})
###########################################################################################################

  
 if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_obs_screening.r!")}