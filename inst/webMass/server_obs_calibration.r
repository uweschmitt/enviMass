if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_obs_screening.r!")}

ranges_cal_plot <- reactiveValues(x = NULL, y = NULL)
dd	<-	reactiveValues() # reactive value ...
dd$d<-matrix(ncol=2,nrow=0,0) 	# ... to indicate & contain model save / removal
dd$entry<-"none_none"			# ... non-redundant changes in the compounds selection
redo_cal<-reactiveValues() 	# reactive value to indicate model save / removal
redo_cal$a<-1

###########################################################################################################
# SPECIFY IONIZATION MODE #################################################################################
observe({ # - A
	input$Ion_mode_Cal 
	if(verbose){cat("\n in A")}
	if(isolate(init$a)=="TRUE"){
		if(logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="calibration"]=="FALSE"){ # e.g., if files were changed / added / ...
			if(verbose){cat("\n in A_1")}
			if(isolate(input$Ion_mode_Cal)=="positive"){
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				measurements<-measurements[measurements[,"Type"]=="calibration",,drop=FALSE]
				measurements<-measurements[measurements[,"Mode"]=="positive",,drop=FALSE]
				measurements<<-measurements
				targets<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");
				intstand<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
				targets<-targets[targets[,"ion_mode"]=="positive",,drop=FALSE]
				intstand<-intstand[intstand[,"ion_mode"]=="positive",,drop=FALSE]
				targets<-targets[targets[,"ID_internal_standard"]!="FALSE",,drop=FALSE] # MUST have an ISTD associated!
				targets<<-targets
				intstand<<-intstand
				if(length(measurements[,"ID"])>0 & length(targets[,"ID"])>0 & length(intstand[,"ID"])>0 ){
					those<-unique(measurements$tag2)
					if(all(those!="FALSE")){
						those<-c("none",those)
						updateSelectInput(session,"Cal_file_set","Specify calibration file group",choices = those, selected = those[1])
					}else{ # trigger warning
						cat("all calibration groups must have a tag2 other than FALSE!")
					}
				}
				if(length(targets[,"ID"])==0 || length(intstand[,"ID"])==0 ){
					shinyjs:::info("No valid targets and/or internal standard compounds found for a quantification!");
				}
			}
			if(isolate(input$Ion_mode_Cal)=="negative"){
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				measurements<-measurements[measurements[,"Type"]=="calibration",,drop=FALSE]
				measurements<-measurements[measurements[,"Mode"]=="negative",,drop=FALSE]
				measurements<<-measurements
				targets<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");
				intstand<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
				targets<-targets[targets[,"ion_mode"]=="negative",,drop=FALSE]
				intstand<-intstand[intstand[,"ion_mode"]=="negative",,drop=FALSE]
				targets<-targets[targets[,"ID_internal_standard"]!="FALSE",,drop=FALSE] # MUST have an ISTD associated!
				targets<<-targets
				intstand<<-intstand
				if(length(measurements[,"ID"])>0 & length(targets[,"ID"])>0 & length(intstand[,"ID"])>0 ){
					those<-unique(measurements$tag2)
					if(all(those!="FALSE")){
						those<-c("none",those)
						updateSelectInput(session,"Cal_file_set","Specify calibration file group",choices = those, selected = those[1])
					}else{ # trigger warning
						cat("all calibration groups must have a tag2 other than FALSE!")
					}
				}
				if(length(targets[,"ID"])==0 || length(intstand[,"ID"])==0 ){
					shinyjs:::info("No valid targets and/or internal standard compounds found for a quantification!");
				}
			}
		}else{
			if(isolate(input$Ion_mode_Cal)!="none"){
				shinyjs:::info("Calibration files have been modified or compounds added. Workflow recalculation including the calibration step (enabled?) required.");
				cat("\n Calibration files have been modified or compounds added. Recalculation required!")
			}
		}	
	}
})
###########################################################################################################

###########################################################################################################
# SPECIFY CALIBRATION GROUP ###############################################################################
observe({ # - B
	input$Cal_file_set
	if(verbose){cat("\n in B")}
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
				# Update IS compounds; only use compounds that were screened & appear in the compound tables
				IS_names<-unique(intstand[,"Name"])
				IS_names<-IS_names[order(IS_names)]
				IS_names<-c("none",IS_names)
				updateSelectInput(session,inputId="Cal_IS_name",label="Internal standard name",choices=IS_names,selected = IS_names[1])
				IS_IDs<-unique(intstand[,"ID"])
				IS_IDs<-IS_IDs[order(IS_IDs)]
				IS_IDs<-c("none",IS_IDs)
				updateSelectInput(session,inputId="Cal_IS_ID",label="Internal standard ID",choices=IS_IDs,selected = IS_IDs[1])
				# Update target compounds; only use compounds that were screened & appear in the compound tables
				target_names<-unique(targets[,"Name"])
				target_names<-target_names[order(target_names)]
				target_names<-c("none",target_names)
				updateSelectInput(session,inputId="Cal_target_name",label="Target name",choices=target_names,selected = target_names[1])	
				target_IDs<-unique(targets[,"ID"])
				target_IDs<-target_IDs[order(target_IDs)]
				target_IDs<-c("none",target_IDs)
				updateSelectInput(session,inputId="Cal_target_ID",label="Target ID",choices=target_IDs,selected = target_IDs[1])
				# clean, load & update available calibration models
				if(any(objects()=="cal_models_pos")){rm(cal_models_pos)}
				if(any(objects(envir=as.environment(".GlobalEnv"))=="cal_models_pos")){rm(cal_models_pos,envir=as.environment(".GlobalEnv"))}	
				if(!file.exists(file.path(logfile[[1]],"quantification",paste("cal_models_pos_",isolate(input$Cal_file_set),sep="")))){
					cal_models_pos<<-list()
					cal_models_pos[[1]]<<-list()
					names(cal_models_pos)<<-isolate(input$Cal_file_set)
					dump("cal_models_pos",file=file.path(logfile[[1]],"quantification",paste("cal_models_pos_",isolate(input$Cal_file_set),sep="")),envir=as.environment(".GlobalEnv"));			
					if(verbose){cat("\n Made new cal_models_pos_...")}
				}else{
					source(file=file.path(logfile[[1]],"quantification",paste("cal_models_pos_",isolate(input$Cal_file_set),sep="")),local=as.environment(".GlobalEnv"));
				}	
				# how many models left to calibrate? concerns only targets with a valid ISTD link
				at_Cal<-isolate(input$Cal_file_set)
				use_cal<-which(names(cal_models_pos)==at_Cal) # well, the first entry ... just in case different calibration groups are merged into a list at some point (= makes saving too slow).
				if(length(use_cal)>0){
					output$number_missing_models<-renderText({
						paste("Number of Targets with missing calibration models: ",
							as.character(sum(is.na(match(
								paste("_",targets[,"ID_internal_standard"],"_",targets[,"ID"],"_",sep=""),
								names(cal_models_pos[[use_cal]])
							))))
						,sep="")
					})	
				}else{
					output$number_missing_models<-renderText({
						paste("Number of Targets with missing calibration models: ",
							as.character(length(targets[,1]))
						,sep="")
					})				
				}
			}else{ # not available
				if((isolate(input$Ion_mode_Cal)!="negative")){
					updateSelectInput(session,inputId="Cal_IS_name",choices="none",selected = "none")
					updateSelectInput(session,inputId="Cal_IS_ID",choices="none",selected = "none")			
					updateSelectInput(session,inputId="Cal_target_name",choices="none",selected = "none")
					updateSelectInput(session,inputId="Cal_target_ID",choices="none",selected = "none")	
					updateSelectInput(session,inputId="Cal_file_set",selected = "none")	
					shinyjs:::info("No screening results for this calibration file set (positive mode) found - have you run the workflow with the calibration step enabled before?")
					output$number_missing_models<-renderText({"Screening results are missing."})
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
				# Update IS compounds; only use compounds that were screened & appear in the compound tables
				IS_names<-unique(intstand[,2])
				IS_names<-IS_names[order(IS_names)]
				IS_names<-c("none",IS_names)
				updateSelectInput(session,inputId="Cal_IS_name",label="Internal standard name",choices=IS_names,selected = IS_names[1])	
				IS_IDs<-unique(intstand[,1])
				IS_IDs<-IS_IDs[order(IS_IDs)]
				IS_IDs<-c("none",IS_IDs)
				updateSelectInput(session,inputId="Cal_IS_ID",label="Internal standard ID",choices=IS_IDs,selected = IS_IDs[1])
				# Update target compounds; only use compounds that were screened & appear in the compound tables
				target_names<-unique(targets[,2])
				target_names<-target_names[order(target_names)]
				target_names<-c("none",target_names)
				updateSelectInput(session,inputId="Cal_target_name",label="Target name",choices=target_names,selected = target_names[1])
				target_IDs<-unique(targets[,1])				
				target_IDs<-target_IDs[order(target_IDs)]
				target_IDs<-c("none",target_IDs)
				updateSelectInput(session,inputId="Cal_target_ID",label="Target ID",choices=target_IDs,selected = target_IDs[1])
				# clean, load & update available calibration models
				if(any(objects()=="cal_models_neg")){rm(cal_models_neg)}
				if(any(objects(envir=as.environment(".GlobalEnv"))=="cal_models_neg")){rm(cal_models_neg,envir=as.environment(".GlobalEnv"))}	
				if(!file.exists(file.path(logfile[[1]],"quantification",paste("cal_models_neg_",isolate(input$Cal_file_set),sep="")))){
					cal_models_neg<<-list()
					cal_models_neg[[1]]<<-list()
					names(cal_models_neg)<<-isolate(input$Cal_file_set)
					dump("cal_models_neg",file=file.path(logfile[[1]],"quantification",paste("cal_models_neg_",isolate(input$Cal_file_set),sep="")),envir=as.environment(".GlobalEnv"));			
					if(verbose){cat("\n Made new cal_models_neg_...")}
				}
				source(file=file.path(logfile[[1]],"quantification",paste("cal_models_neg_",isolate(input$Cal_file_set),sep="")),local=as.environment(".GlobalEnv"));	
				# how many models left to calibrate? concerns only targets with a valid ISTD link
				at_Cal<-isolate(input$Cal_file_set)
				use_cal<-which(names(cal_models_neg)==at_Cal) # well, the first entry ... just in case different calibration groups are merged into a list at some point (= makes saving too slow).
				if(length(use_cal)>0){
					output$number_missing_models<-renderText({
						paste("Number of Targets with missing calibration models: ",
							as.character(sum(is.na(match(
								paste("_",targets[,"ID_internal_standard"],"_",targets[,"ID"],"_",sep=""),
								names(cal_models_neg[[use_cal]])
							))))
						,sep="")
					})	
				}else{
					output$number_missing_models<-renderText({
						paste("Number of Targets with missing calibration models: ",
							as.character(length(targets[,1]))
						,sep="")
					})				
				}
			}else{ # not available
				if((isolate(input$Ion_mode_Cal)!="positive")){
					updateSelectInput(session,inputId="Cal_IS_name",choices="none",selected = "none")
					updateSelectInput(session,inputId="Cal_IS_ID",choices="none",selected = "none")			
					updateSelectInput(session,inputId="Cal_target_name",choices="none",selected = "none")
					updateSelectInput(session,inputId="Cal_target_ID",choices="none",selected = "none")	
					updateSelectInput(session,inputId="Cal_file_set",selected = "none")	
					shinyjs:::info("No screening results for this calibration file set (negative mode) found - have you run the workflow with the calibration step enabled before?")
					output$number_missing_models<-renderText({"Screening results are missing."})					
				}
			}
		}
	}
})
###########################################################################################################

###########################################################################################################
# DELETE MODELS FOR A CALIBRATION GROUP ###################################################################
observe({ # Update target name & IS_ID - C
	input$Cal_file_set_delete
	if((isolate(init$a)=="TRUE")){	
		at_Cal<-isolate(input$Cal_file_set)
		if(isolate(input$Ion_mode_Cal)=="positive"){
			use_cal<-which(names(cal_models_pos)==at_Cal) # well, the first entry ... just in case different calibration groups are merged into a list at some point (= makes saving too slow).
			if(use_cal!="none"){
				cal_models_pos[[use_cal]]<<-list()
				dump("cal_models_pos",file=file.path(logfile[[1]],"quantification",paste("cal_models_pos_",isolate(input$Cal_file_set),sep="")),envir=as.environment(".GlobalEnv"));			
			}
		}
		if(isolate(input$Ion_mode_Cal)=="negative"){
			use_cal<-which(names(cal_models_neg)==at_Cal) # well, the first entry ... just in case different calibration groups are merged into a list at some point (= makes saving too slow).
			if(use_cal!="none"){
				cal_models_neg[[use_cal]]<<-list()
				dump("cal_models_neg",file=file.path(logfile[[1]],"quantification",paste("cal_models_neg_",isolate(input$Cal_file_set),sep="")),envir=as.environment(".GlobalEnv"));			
			}
		}
	}
})
###########################################################################################################

###########################################################################################################
# MATCH COMPOUND NAMES AND ID SELECTIONS ##################################################################
observe({ # Update target name & IS_ID - C
	input$Cal_target_ID
	if(verbose){cat("\n in C")}
	if((isolate(init$a)=="TRUE")){
		if(isolate(input$Cal_target_ID)!="none"){
			use_this_name<-targets[targets[,1]==isolate(input$Cal_target_ID),"Name"]
			updateSelectInput(session,inputId="Cal_target_name",selected = as.character(use_this_name))
			new_IS_ID<-targets[targets[,1]==isolate(input$Cal_target_ID),"ID_internal_standard"]
			if(verbose){cat(" - ");cat(new_IS_ID)}
			if(length(new_IS_ID)>0){ # just in case ...
				if(new_IS_ID!="FALSE"){
					if(new_IS_ID!=isolate(input$Cal_IS_ID)){
						updateSelectInput(session,inputId="Cal_IS_ID",selected = as.character(new_IS_ID))
					}else{ # if IS_ID the same; update on dd$d must be made here, not in observer for input$Cal_IS_ID which would remain unresponise (no value change)
						isolate(dd$entry<-paste("_",input$Cal_IS_ID,"_",input$Cal_target_ID,"_",sep=""))
					}	
				}
			}		
		}else{	
			updateSelectInput(session,inputId="Cal_target_name",selected="none")
		}
	}	
})

observe({ # Update target ID - D
	input$Cal_target_name
	if(verbose){cat("\n in D")}
	if((isolate(init$a)=="TRUE")&(isolate(input$Cal_target_name)!="none")){
		use_this_ID<-use_this_name<-targets[targets[,2]==isolate(input$Cal_target_name),1]
		updateSelectInput(session,inputId="Cal_target_ID",selected = as.character(use_this_ID))
	}	
})

observe({ # Update IS name - E
	input$Cal_IS_ID
	if(verbose){cat("\n in E")}
	if((isolate(init$a)=="TRUE")){
		if(isolate(input$Cal_IS_ID)!="none"){
			use_this_name<-intstand[intstand[,1]==isolate(input$Cal_IS_ID),2]
			updateSelectInput(session,inputId="Cal_IS_name",selected = as.character(use_this_name))
			if(isolate(dd$entry!=paste("_",input$Cal_IS_ID,"_",input$Cal_target_ID,"_",sep=""))){
				isolate(dd$entry<-paste("_",input$Cal_IS_ID,"_",input$Cal_target_ID,"_",sep=""))
			}
		}else{
			updateSelectInput(session,inputId="Cal_IS_name",selected="none")
		}
	}
})

observe({ # Update IS ID - F
	input$Cal_IS_name
	if(verbose){cat("\n in F")}
	if((isolate(init$a)=="TRUE")&(isolate(input$Cal_IS_name)!="none")){
		use_this_ID<-intstand[intstand[,2]==isolate(input$Cal_IS_name),1]
		updateSelectInput(session,inputId="Cal_IS_ID",selected = as.character(use_this_ID))
	}	
})
###########################################################################################################

###########################################################################################################
# FORWARD BACKWARD ETC ####################################################################################
observe({ # - G
	input$Cal_next
	if(verbose){cat("\n in   G")}
	if((isolate(init$a)=="TRUE") & (isolate(input$Cal_file_set)!="none")){
		is_at_targetID<-(isolate(input$Cal_target_ID))
		is_at_ISID<-(isolate(input$Cal_IS_ID))
		in_table<-which( ((targets[,"ID"]==is_at_targetID)&(targets[,"ID_internal_standard"]==is_at_ISID)) )
		at_target_ID<-"none"
		if( (length(in_table)==0) || (is.na(in_table[1])) ){ # match with is_at_ISID not existing in table; reset to first available entry
			if(verbose){cat("\n in G_1")}
			in_table<-which(targets[,1]==is_at_targetID)
			in_table<-in_table[1] # should only be one entry!
			if( (length(in_table)==0) || (is.na(in_table[1])) ){ # for any invalid entry, start from beginning
				if(verbose){cat("\n in G_1_1")}
				in_table<-1
			}
			if(in_table<length(targets[,"ID"])){
				in_table<-(in_table+1)		
			}
			at_target_ID<-as.character(targets[in_table,"ID"])	
		}else{ # match existing; get next entry
			if(verbose){cat("\n in G_2")}
			if(in_table<length(targets[,"ID"])){
				in_table<-(in_table+1)
				at_target_ID<-as.character(targets[in_table,"ID"])			
			}else{
				at_target_ID<-"none"
			}
		}
		updateSelectInput(session, inputId="Cal_target_ID", selected = at_target_ID)
	}
})

observe({ # - H
	input$Cal_previous
	if(verbose){cat("\n in   H")}
	if((isolate(init$a)=="TRUE")&(isolate(input$Cal_target_name)!="none")& (isolate(input$Cal_file_set)!="none")){
		is_at_targetID<-(isolate(input$Cal_target_ID))
		is_at_ISID<-(isolate(input$Cal_IS_ID))
		in_table<-which( ((targets[,1]==is_at_targetID)&(targets[,6]==is_at_ISID)) )
		at_target_ID<-"none"
		if( (length(in_table)==0) || (is.na(in_table[1])) ){ # match with is_at_ISID not existing in table; reset to first available entry
			in_table<-which(targets[,1]==is_at_targetID)
			in_table<-in_table[1] # should only be one entry!
			if( (length(in_table)==0) || (is.na(in_table[1])) ){ # for any invalid entry, start from beginning
				in_table<-1
			}
			if(in_table>1){
				in_table<-(in_table-1)		
			}
			at_target_ID<-as.character(targets[in_table,1])
			at_target_ID<-as.character(targets[in_table,1])	
		}else{ # match existing; get next entry
			if(in_table>1){
				in_table<-(in_table-1)
				at_target_ID<-as.character(targets[in_table,1])			
			}else{
				at_target_ID<-"none"
			}
		}	
		updateSelectInput(session, inputId="Cal_target_ID", selected = at_target_ID)
	}
})

observe({ # - I
	input$Cal_first
	if(verbose){cat("\n in   I")}
	if((isolate(init$a)=="TRUE")& (isolate(input$Cal_file_set)!="none")){
		at_target_ID<-"none"
		if(length(targets[,1])>0){ # match with is_at_ISID not existing in table; reset to first available entry
			in_table<-1
			at_target_ID<-as.character(targets[in_table,1])	
		}
		updateSelectInput(session, inputId="Cal_target_ID", selected = at_target_ID)
	}
})

observe({ # find last entry - J
	input$Cal_last
	if(verbose){cat("\n in   J")}
	if((isolate(init$a)=="TRUE")& (isolate(input$Cal_file_set)!="none")){
		at_target_ID<-"none"
		if(length(targets[,1])>0){ # match with is_at_ISID not existing in table; reset to first available entry
			in_table<-length(targets[,1])
			at_target_ID<-as.character(targets[in_table,1])
		}
		updateSelectInput(session, inputId="Cal_target_ID", selected = at_target_ID)
	}
})

observe({ # find next entry with a missing calibration model
	input$Cal_next_missing
	if(verbose){cat("\n in  Missing")}
	if((isolate(init$a)=="TRUE")& (isolate(input$Cal_file_set)!="none")){
		if(length(targets[,1])>0){ # match with is_at_ISID not existing in table; reset to first available entry
			at_Cal<-isolate(input$Cal_file_set)
			# POSITIVE MODE
			if(isolate(input$Ion_mode_Cal)=="positive"){
				use_cal<-which(names(cal_models_pos)==at_Cal) # well, the first entry ... just in case different calibration groups are merged into a list at some point (= makes saving too slow).
				if(length(names(cal_models_pos[[use_cal]]))>0){
					in_table<-length(targets[,1])
					found_one<-FALSE
					at<-0
					is_at_targetID<-(isolate(input$Cal_target_ID))
					if(is_at_targetID=="none"){
						b<-1
					}else{
						b<-which(targets[,"ID"]==is_at_targetID)
						b<-(b+1)
						if(b>in_table){b<-1} # start from beginning of table
					}
					for(i in b:in_table){
						at<-i
						if(targets[i,"ID_internal_standard"]=="FALSE"){next}
						if(!any(names(cal_models_pos[[use_cal]])==paste("_",targets[i,"ID_internal_standard"],"_",targets[i,"ID"],"_",sep=""))){
							found_one<-TRUE
							break;
						}
					}
					if(found_one){
						at_target_ID<-as.character(targets[i,"ID"])
						updateSelectInput(session, inputId="Cal_target_ID", selected = at_target_ID)
					}else{
						at_target_ID<-"none"
						updateSelectInput(session, inputId="Cal_target_ID", selected = at_target_ID)
						shinyjs:::info("No Target compounds without calibration model left");
					}
				}else{
					at_target_ID<-as.character(targets[1,"ID"])
					updateSelectInput(session, inputId="Cal_target_ID", selected = at_target_ID)			
				}
			}
			# NEGATIVE MODE
			if(isolate(input$Ion_mode_Cal)=="negative"){
				use_cal<-which(names(cal_models_neg)==at_Cal) # well, the first entry ... just in case different calibration groups are merged into a list at some point (= makes saving too slow).
				if(length(names(cal_models_neg[[use_cal]]))>0){
					in_table<-length(targets[,1])
					found_one<-FALSE
					at<-0
					is_at_targetID<-(isolate(input$Cal_target_ID))
					if(is_at_targetID=="none"){
						b<-1
					}else{
						b<-which(targets[,"ID"]==is_at_targetID)
						b<-(b+1)
						if(b>in_table){b<-1} # start from beginning of table
					}
					for(i in b:in_table){
						at<-i
						if(targets[i,"ID_internal_standard"]=="FALSE"){next}
						if(!any(names(cal_models_neg[[use_cal]])==paste("_",targets[i,"ID_internal_standard"],"_",targets[i,"ID"],"_",sep=""))){
							found_one<-TRUE
							break;
						}
					}
					if(found_one){
						at_target_ID<-as.character(targets[i,"ID"])
						updateSelectInput(session, inputId="Cal_target_ID", selected = at_target_ID)
					}else{
						at_target_ID<-"none"
						updateSelectInput(session, inputId="Cal_target_ID", selected = at_target_ID)
						shinyjs:::info("No Target compounds without calibration model left");
					}
				}else{
					at_target_ID<-as.character(targets[1,"ID"])
					updateSelectInput(session, inputId="Cal_target_ID", selected = at_target_ID)			
				}
			}
		}
	}
})

observe({ # find first entry with a missing calibration model
	input$Cal_missing
	if(verbose){cat("\n in  Missing")}
	if((isolate(init$a)=="TRUE")& (isolate(input$Cal_file_set)!="none")){
		at_target_ID<-"none"
		if(length(targets[,1])>0){ # match with is_at_ISID not existing in table; reset to first available entry
			at_Cal<-isolate(input$Cal_file_set)
			# POSITIVE MODE
			if(isolate(input$Ion_mode_Cal)=="positive"){
				use_cal<-which(names(cal_models_pos)==at_Cal) # well, the first entry ... just in case different calibration groups are merged into a list at some point (= makes saving too slow).
				if(length(names(cal_models_pos[[use_cal]]))>0){
					in_table<-length(targets[,1])
					found_one<-FALSE
					at<-0
					for(i in 1:in_table){
						at<-i
						if(targets[i,"ID_internal_standard"]=="FALSE"){next}
						if(!any(names(cal_models_pos[[use_cal]])==paste("_",targets[i,"ID_internal_standard"],"_",targets[i,"ID"],"_",sep=""))){
							found_one<-TRUE
							break;
						}
					}
					if(found_one){
						at_target_ID<-as.character(targets[i,"ID"])
						updateSelectInput(session, inputId="Cal_target_ID", selected = at_target_ID)
					}else{
						at_target_ID<-"none"
						updateSelectInput(session, inputId="Cal_target_ID", selected = at_target_ID)
						shinyjs:::info("No Target compounds without calibration model left");
					}
				}else{
					at_target_ID<-as.character(targets[1,"ID"])
					updateSelectInput(session, inputId="Cal_target_ID", selected = at_target_ID)			
				}
			}
			# NEGATIVE MODE
			if(isolate(input$Ion_mode_Cal)=="negative"){
				use_cal<-which(names(cal_models_neg)==at_Cal) # well, the first entry ... just in case different calibration groups are merged into a list at some point (= makes saving too slow).
				if(length(names(cal_models_neg[[use_cal]]))>0){
					in_table<-length(targets[,1])
					found_one<-FALSE
					at<-0
					for(i in 1:in_table){
						at<-i
						if(targets[i,"ID_internal_standard"]=="FALSE"){next}
						if(!any(names(cal_models_neg[[use_cal]])==paste("_",targets[i,"ID_internal_standard"],"_",targets[i,"ID"],"_",sep=""))){
							found_one<-TRUE
							break;
						}
					}
					if(found_one){
						at_target_ID<-as.character(targets[i,"ID"])
						updateSelectInput(session, inputId="Cal_target_ID", selected = at_target_ID)
					}else{
						at_target_ID<-"none"
						updateSelectInput(session, inputId="Cal_target_ID", selected = at_target_ID)
						shinyjs:::info("No Target compounds without calibration model left");
					}
				}else{
					at_target_ID<-as.character(targets[1,"ID"])
					updateSelectInput(session, inputId="Cal_target_ID", selected = at_target_ID)			
				}
			}
		}
	}
})

###########################################################################################################

###########################################################################################################
# RELOAD DATA (e.g., when IS intensity ranges were changed) ###############################################
observe({ # - Reload
	input$reload_Cal
	if(isolate(init$a)=="TRUE"){
		if(isolate(input$Cal_file_set)!="none"){
			if(verbose){cat("\nReloaded:")}
			if(
				(isolate(input$Ion_mode_Cal)=="positive")
			){
				targets<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");
				intstand<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
				targets<-targets[targets[,8]=="positive",,drop=FALSE]
				intstand<-intstand[intstand[,7]=="positive",,drop=FALSE]
				targets<-targets[targets[,6]!="FALSE",,drop=FALSE]
				targets<<-targets
				intstand<<-intstand
				if(verbose){cat(" positive")}
			}
			if(
				(isolate(input$Ion_mode_Cal)=="negative")		
			){
				targets<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");
				intstand<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
				targets<-targets[targets[,8]=="negative",,drop=FALSE]
				intstand<-intstand[intstand[,7]=="negative",,drop=FALSE]
				targets<-targets[targets[,6]!="FALSE",,drop=FALSE]
				targets<<-targets
				intstand<<-intstand
				if(verbose){cat(" negative")}
			}		
		}
	}
})	
###########################################################################################################

###########################################################################################################
# RETRIEVE SETS ###########################################################################################
observe({ # - K
	dd$entry
	input$reload_Cal
	input$cal_model_bound_low
	input$cal_model_bound_low_value
	input$cal_model_bound_up
	input$cal_model_bound_up_value
	if(verbose){cat("\n in K")}
		if(isolate(init$a)=="TRUE"){
		# anything available from calibration screening? Results still in measurements?
		something<-FALSE
		if(isolate(input$Ion_mode_Cal)=="positive"){ # check if compounds in tables are also found in screening outcomes!
			if(
				any(as.character(results_screen_target_pos_cal[[1]][,1])==isolate(input$Cal_target_ID)) &
				any(as.character(results_screen_IS_pos_cal[[1]][,1])==isolate(input$Cal_IS_ID)) 
			){
				something<-TRUE
			}
		}
		if(isolate(input$Ion_mode_Cal)=="negative"){ # check if compounds in tables are also found in screening outcomes!
			if(
				any(as.character(results_screen_target_neg_cal[[1]][,1])==isolate(input$Cal_target_ID)) &
				any(as.character(results_screen_IS_neg_cal[[1]][,1])==isolate(input$Cal_IS_ID))
			){
				something<-TRUE
			}
		}	
		if(verbose){cat(something)}
		if(
			(isolate(init$a)=="TRUE") & 
			(something) & 
			(isolate(input$Cal_IS_ID)!="none") & 
			(isolate(input$Cal_target_ID)!="none") &  
			(isolate(input$Cal_file_set)!="none")
		){	
			# POSITIVE ##################################################################
			if(isolate(input$Ion_mode_Cal)=="positive"){	
				if(verbose){cat("\n in K_1")}
				IS_ID<-isolate(input$Cal_IS_ID)
				target_ID<-isolate(input$Cal_target_ID)
				at_Cal<-isolate(input$Cal_file_set)
				#target_ID<-"315"; IS_ID<-"74";at_Cal<-"A"			
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
							if(measurements[measurements[,"ID"]==res_IS_pos_screen_cal[[at_entry]][[j]][[1]]$file_ID,]$tag2==at_Cal){
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
				if(verbose){cat("\n in K_2")}
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
				if(verbose){cat("\n in K_2_1")}
				target_in_file<-c()
				target_with_peak<-c()
				target_with_score<-c()
				if(length(res_target_pos_screen_cal[[at_entry]])>0){
					for(j in 1:length(res_target_pos_screen_cal[[at_entry]])){
						if(length(res_target_pos_screen_cal[[at_entry]][[j]])>0){						
							if(measurements[measurements[,"ID"]==res_target_pos_screen_cal[[at_entry]][[j]][[1]]$file_ID,]$tag2==at_Cal){
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
				if(verbose){cat("\n in K_3")}
				# derive pairs ##########################################################
				mat_cal<-matrix(nrow=0,ncol=13)
				colnames(mat_cal)<-c("Pair #","Target intensity","IS intensity","Intensity ratio","Concentration","Target score","IS score","Used?",
					"Target RT [s]","IS RT [s]","Target peak ID","IS peak ID","File ID")
				if(	(length(target_in_file)>0) & (length(IS_in_file)>0)	){
					for(i in 1:length(target_in_file)){
						if(any(IS_in_file==target_in_file[i])){
							those<-which(IS_in_file==target_in_file[i])
							mat_cal<-rbind(mat_cal,
								cbind(
									rep(1,length(those)),
									round(rep(profileList_pos_cal[[2]][target_with_peak[i],2],length(those)),digits=0), # intensity target
									round(profileList_pos_cal[[2]][IS_with_peak[those],2],digits=0), 					# intensity IS
									round((profileList_pos_cal[[2]][target_with_peak[i],2]/(profileList_pos_cal[[2]][IS_with_peak[those],2])),digits=3), # ratio									
									rep(as.numeric(
										measurements[measurements[,"ID"]==target_in_file[i],]$tag1	
									),length(those)), # concentration
									rep(target_with_score[i],length(those)), # score target
									IS_with_score[those],
									rep(1,length(those)), # used?
									round(rep(profileList_pos_cal[[2]][target_with_peak[i],3],length(those)),digits=2), # Target RT
									round(rep(profileList_pos_cal[[2]][IS_with_peak[those],3],length(those)),digits=2),	# IS RT
									round(rep(profileList_pos_cal[[2]][target_with_peak[i],4],length(those)),digits=0),	# Target peak ID
									round(rep(profileList_pos_cal[[2]][IS_with_peak[those],4],length(those)),digits=0),	# IS peak ID
									round(rep(profileList_pos_cal[[2]][target_with_peak[i],6],length(those)),digits=0)	# File ID
								,deparse.level = 0),
							deparse.level = 0)
						}
					}
				}
				rownames(mat_cal)<-NULL
				if(verbose){cat("\n in K_4")}
				# filter ################################################################
				mat_cal<-mat_cal[!duplicated(mat_cal),,drop=FALSE] # same entries? - remove
				mat_cal<-mat_cal[
					order(mat_cal[,5],mat_cal[,11],mat_cal[,12],mat_cal[,7],mat_cal[,6],decreasing=TRUE)
				,,drop=FALSE]
				mat_cal<-mat_cal[!duplicated(mat_cal[,c(11,12),drop=FALSE]),,drop=FALSE] # same peaks in different combinations - remove
				mat_cal[,1]<-(1:length(mat_cal[,1]))
				min_int<-as.numeric(intstand[intstand[,1]==IS_ID,17])
				if(min_int!=0){min_int<-10^min_int}
				max_int<-as.numeric(intstand[intstand[,1]==IS_ID,18])
				max_int<-(max_int^10)
				mat_cal[mat_cal[,3]<min_int,8]<-0
				mat_cal[mat_cal[,3]>max_int,8]<-0
				# adapt point selection to existing model (if any) ######################
				use_cal<-which(names(cal_models_pos)==at_Cal) # well, the first entry ... just in case different calibration groups are merged into a list at some point (= makes saving too slow).
				if(length(names(cal_models_pos[[use_cal]]))>0){				
					use_precision<-isolate(input$use_precision)
					at_model<-which(names(cal_models_pos[[use_cal]])==paste("_",IS_ID,"_",target_ID,"_",sep=""))
					if(length(at_model)>0){
						cal_models_pos[[use_cal]][[at_model]]$data
						for(k in 1:length(mat_cal[,1])){
							if(!any(
								(cal_models_pos[[use_cal]][[at_model]]$data$resp==mat_cal[k,"Concentration"]) &
								(abs(cal_models_pos[[use_cal]][[at_model]]$data$lin-mat_cal[k,"Intensity ratio"])<use_precision)	
							)){
								mat_cal[k,"Used?"]<-0
								if(verbose){cat(".")}
							}else{
								mat_cal[k,"Used?"]<-1
								if(verbose){cat("*")}							
							}
						}
					}
				}
				# adapt point selection also to new ratio bounds ########################
				if(isolate(input$cal_model_bound_low)){	
					mat_cal[
						mat_cal[,"Intensity ratio"]<isolate(input$cal_model_bound_low_value)
					,"Used?"]<-0
				}
				if(isolate(input$cal_model_bound_up)){
					mat_cal[
						mat_cal[,"Intensity ratio"]>isolate(input$cal_model_bound_up_value)
					,"Used?"]<-0					
				}						
				#########################################################################
				mat_cal<<-mat_cal
				isolate(dd$d<-mat_cal)
			}
			# NEGATIVE ##################################################################
			if(isolate(input$Ion_mode_Cal)=="negative"){	
				if(verbose){cat("\n in K_negative_1")}
				IS_ID<-isolate(input$Cal_IS_ID)
				target_ID<-isolate(input$Cal_target_ID)
				at_Cal<-isolate(input$Cal_file_set)
				#target_ID<-"4"; IS_ID<-"693";at_Cal<-"A"			
				if(verbose){cat("\n in K_negative_1_2")}
				# extract IS peaks ######################################################
				IS_adduct<-intstand[intstand[,1]==IS_ID,19]
				IS_peak<-as.numeric(intstand[intstand[,1]==IS_ID,20])
				at_entry<-FALSE
				for(i in 1:length(names(res_IS_neg_screen_cal))){ # where?
					if(
						(strsplit(names(res_IS_neg_screen_cal)[i],"_")[[1]][1]==IS_ID) &
						(strsplit(names(res_IS_neg_screen_cal)[i],"_")[[1]][2]==IS_adduct)
					){
						at_entry<-i;break;
					}
				}
				if(verbose){cat("\n in K_negative_1_3")}
				IS_in_file<-c()
				IS_with_peak<-c()
				IS_with_score<-c()
				if(length(res_IS_neg_screen_cal[[at_entry]])>0){
					for(j in 1:length(res_IS_neg_screen_cal[[at_entry]])){
						if(length(res_IS_neg_screen_cal[[at_entry]][[j]])>0){						
							if(measurements[measurements[,"ID"]==res_IS_neg_screen_cal[[at_entry]][[j]][[1]]$file_ID,]$tag2==at_Cal){
								for(k in 1:length(res_IS_neg_screen_cal[[at_entry]][[j]])){
									if(any(res_IS_neg_screen_cal[[at_entry]][[j]][[k]]$Peaks[,1]==IS_peak)){
										that<-which(res_IS_neg_screen_cal[[at_entry]][[j]][[k]]$Peaks[,1]==IS_peak)
										IS_in_file<-c(IS_in_file,res_IS_neg_screen_cal[[at_entry]][[j]][[k]]$file_ID)
										IS_with_peak<-c(IS_with_peak,res_IS_neg_screen_cal[[at_entry]][[j]][[k]]$Peaks[that,2])
										IS_with_score<-c(IS_with_score,res_IS_neg_screen_cal[[at_entry]][[j]][[k]]$score_1)	
									}
								}
							}
						}
					}
				}
				if(verbose){cat("\n in K_negative_5")}
				# extract target peaks ##################################################
				target_adduct<-targets[targets[,1]==target_ID,20]
				target_peak<-as.numeric(targets[targets[,1]==target_ID,21])
				at_entry<-FALSE
				for(i in 1:length(names(res_target_neg_screen_cal))){ # where?
					if(
						(strsplit(names(res_target_neg_screen_cal)[i],"_")[[1]][1]==target_ID) &
						(strsplit(names(res_target_neg_screen_cal)[i],"_")[[1]][2]==target_adduct)
					){
						at_entry<-i;break;
					}
				}
				if(verbose){cat("\n in K_negative_5_1")}
				target_in_file<-c()
				target_with_peak<-c()
				target_with_score<-c()
				if(length(res_target_neg_screen_cal[[at_entry]])>0){
					for(j in 1:length(res_target_neg_screen_cal[[at_entry]])){
						if(length(res_target_neg_screen_cal[[at_entry]][[j]])>0){						
							if(measurements[measurements[,"ID"]==res_target_neg_screen_cal[[at_entry]][[j]][[1]]$file_ID,]$tag2==at_Cal){
								for(k in 1:length(res_target_neg_screen_cal[[at_entry]][[j]])){
									if(any(res_target_neg_screen_cal[[at_entry]][[j]][[k]]$Peaks[,1]==target_peak)){
										that<-which(res_target_neg_screen_cal[[at_entry]][[j]][[k]]$Peaks[,1]==target_peak)
										target_in_file<-c(target_in_file,res_target_neg_screen_cal[[at_entry]][[j]][[k]]$file_ID)
										target_with_peak<-c(target_with_peak,res_target_neg_screen_cal[[at_entry]][[j]][[k]]$Peaks[that,2])
										target_with_score<-c(target_with_score,res_target_neg_screen_cal[[at_entry]][[j]][[k]]$score_1)
									}
								}
							}
						}
					}
				}
				if(verbose){cat("\n in K_negative_6")}
				# derive pairs ##########################################################
				mat_cal<-matrix(nrow=0,ncol=13)
				colnames(mat_cal)<-c("Pair #","Target intensity","IS intensity","Intensity ratio","Concentration","Target score","IS score","Used?",
				"Target RT [s]","IS RT [s]","Target peak ID","IS peak ID","File ID")
				if(	(length(target_in_file)>0) & (length(IS_in_file)>0)	){
					for(i in 1:length(target_in_file)){
						if(any(IS_in_file==target_in_file[i])){
							those<-which(IS_in_file==target_in_file[i])
							mat_cal<-rbind(mat_cal,
								cbind(
									rep(1,length(those)),
									round(rep(profileList_neg_cal[[2]][target_with_peak[i],2],length(those)),digits=0), # intensity target
									round(profileList_neg_cal[[2]][IS_with_peak[those],2],digits=0), 					# intensity IS
									round((profileList_neg_cal[[2]][target_with_peak[i],2]/(profileList_neg_cal[[2]][IS_with_peak[those],2])),digits=2), # ratio									
									rep(as.numeric(
										measurements[measurements[,"ID"]==target_in_file[i],]$tag1	
									),length(those)), # concentration
									rep(target_with_score[i],length(those)), # score target
									IS_with_score[those],
									rep(1,length(those)), # used?
									round(rep(profileList_neg_cal[[2]][target_with_peak[i],3],length(those)),digits=2), # Target RT
									round(rep(profileList_neg_cal[[2]][IS_with_peak[those],3],length(those)),digits=2),	# IS RT
									round(rep(profileList_neg_cal[[2]][target_with_peak[i],4],length(those)),digits=0),	# Target peak ID
									round(rep(profileList_neg_cal[[2]][IS_with_peak[those],4],length(those)),digits=0),	# IS peak ID
									round(rep(profileList_neg_cal[[2]][target_with_peak[i],6],length(those)),digits=0)	# File ID
								,deparse.level = 0),
							deparse.level = 0)
						}
					}
				}
				rownames(mat_cal)<-NULL
				if(verbose){cat("\n in K_negative_7")}
				# filter ################################################################
				mat_cal<-mat_cal[!duplicated(mat_cal),,drop=FALSE] # same entries? - remove
				mat_cal<-mat_cal[
					order(mat_cal[,5],mat_cal[,11],mat_cal[,12],mat_cal[,7],mat_cal[,6],decreasing=TRUE)
				,]
				mat_cal<-mat_cal[!duplicated(mat_cal[,c(11,12)]),,drop=FALSE] # same peaks in different combinations - remove
				mat_cal[,1]<-(1:length(mat_cal[,1]))
				min_int<-as.numeric(intstand[intstand[,1]==IS_ID,17])
				if(min_int!=0){min_int<-10^min_int}
				max_int<-as.numeric(intstand[intstand[,1]==IS_ID,18])
				max_int<-(max_int^10)
				mat_cal[mat_cal[,3]<min_int,8]<-0
				mat_cal[mat_cal[,3]>max_int,8]<-0
				# adapt point selection to existing model (if any) ######################
				use_cal<-which(names(cal_models_neg)==at_Cal) # well, the first entry ... just in case different calibration groups are merged into a list at some point (= makes saving too slow).
				if(length(names(cal_models_neg[[use_cal]]))>0){				
					use_precision<-isolate(input$use_precision)
					at_model<-which(names(cal_models_neg[[use_cal]])==paste("_",IS_ID,"_",target_ID,"_",sep=""))
					if(length(at_model)>0){
						cal_models_neg[[use_cal]][[at_model]]$data
						for(k in 1:length(mat_cal[,1])){
							if(!any(
								(cal_models_neg[[use_cal]][[at_model]]$data$resp==mat_cal[k,"Concentration"]) &
								(abs(cal_models_neg[[use_cal]][[at_model]]$data$lin-mat_cal[k,"Intensity ratio"])<use_precision)	
							)){
								mat_cal[k,"Used?"]<-0
								if(verbose){cat(".")}
							}
						}
					}
				}
				#########################################################################
				mat_cal<<-mat_cal
				isolate(dd$d<-mat_cal)
			}
		}else{
			isolate(dd$d<-matrix(ncol=2,nrow=0,0))
		}	
	}
})				
###########################################################################################################				

###########################################################################################################
# PLOT SETS, MAKE MODELS & OUTPUT TABLE ###################################################################
observe({ # - L output table
	redo_cal$a
	dd$d
	if(verbose){cat("\n in L")}
	if((isolate(init$a)=="TRUE") & (isolate(input$Cal_IS_ID)!="none")  & (isolate(input$Cal_target_ID)!="none")& (isolate(input$Cal_file_set)!="none")){					
		# generate outputs ######################################################
		if(length(isolate(dd$d[,1]))>0){
			output$cal_table <- DT::renderDataTable(
				datatable(
					isolate(dd$d),selection =c('single'),options = list(lengthMenu = c(25,50,100))
				)%>%
				formatStyle('Used?',backgroundColor = styleInterval(0.5, c('orange', 'lightgreen')))
			)
		}else{
			output$cal_table <- renderDataTable(
				DT::datatable(as.data.frame(cbind("")),selection = 'single',rownames=FALSE,colnames="No target screening results available")
			)
		}	
	}	
})

observe({ # - M plot
	input$cal_model
	input$cal_model_0intercept
	input$cal_model_weight
	input$cal_model_bound_low
	input$cal_model_bound_low_value
	input$cal_model_bound_up
	input$cal_model_bound_up_value
	redo_cal$a
	dd$d
	if(verbose){cat("\n in M: plotting")}
	if((isolate(init$a)=="TRUE") & (isolate(input$Cal_IS_ID)!="none") & (isolate(input$Cal_target_ID)!="none")& (isolate(input$Cal_file_set)!="none")){					
		# generate outputs ######################################################
		if(length(isolate(dd$d[,1]))>1){
			output$cal_plot <- renderPlot({
				plot(isolate(dd$d[,"Concentration"]), isolate(dd$d[,"Intensity ratio"]),
					xlab="Concentration",ylab="Intensity ratio",pch=19,
					xlim=ranges_cal_plot$y,ylim=ranges_cal_plot$x,
					main="Draw rectangles and double-click into them to zoom, double-click again to zoom out.",cex.main=1,col="white")#,yaxs="i",xaxs="i")
				abline(h=0,col="lightgrey")
				abline(v=0,col="lightgrey")				
				points(isolate(dd$d[dd$d[,"Used?"]==1,"Concentration"]), isolate(dd$d[dd$d[,"Used?"]==1,"Intensity ratio"]),col="black",pch=19)
				points(isolate(dd$d[dd$d[,"Used?"]==0,"Concentration"]), isolate(dd$d[dd$d[,"Used?"]==0,"Intensity ratio"]),col="gray",pch=19)		
				# add lower & upper bounds
				if(isolate(input$cal_model_bound_low)){
					abline(h=isolate(input$cal_model_bound_low_value),col="red",lty=2)
				}
				if(isolate(input$cal_model_bound_up)){
					abline(h=isolate(input$cal_model_bound_up_value),col="red",lty=3)
				}					
				if((!is.null(ranges_cal_plot$x))||(!is.null(ranges_cal_plot$y))){
					mtext("Now zoomed in",side=3,col="gray")
				}
				if(sum(isolate(dd$d[,"Used?"]))>=2){ # at least two data points!
					resp<-(isolate(dd$d[dd$d[,"Used?"]==1,"Concentration"]))	# concentration
					lin<-(isolate(dd$d[dd$d[,"Used?"]==1,"Intensity ratio"])) 	# Intensity ratio target/IS
					if(isolate(input$cal_model_weight)){
						use_weights<-(isolate(dd$d[dd$d[,"Used?"]==1,"Target intensity"])) 
						use_weights<-(1/use_weights)
					}else{
						use_weights<-NULL
					}
					if(isolate(input$cal_model)=="linear"){
						if(isolate(input$cal_model_0intercept)){
							cal_model<<-lm(resp~0+lin,weights=use_weights)
							abline(
								a=0,
								b=(1/cal_model$coefficients[1]),
								col="red",lwd=2)
						}else{
							cal_model<<-lm(resp~lin,weights=use_weights)
							abline(
								a=(-1*cal_model$coefficients[1]/cal_model$coefficients[2]),
								b=(1/cal_model$coefficients[2]),
								col="red",lwd=2)
						}
					}else{
						quad<-(lin^2)
						if(isolate(input$cal_model_0intercept)){
							cal_model<<-lm(resp~0+lin+quad,weights=use_weights)					
						}else{
							cal_model<<-lm(resp~lin+quad,weights=use_weights)							
						}
						for_y<-seq(from=0,to=(max(lin)*2),length.out=100)
						for_y2<-(for_y^2)
						for_x<-predict(cal_model,list(lin=for_y,quad=for_y2))
						lines(for_x,for_y,col="red",lwd=2)
					}					
				}else{
					if(any(objects()=="cal_model")){cat("\n cal_model in invalid environment found! DEBUG!");rm(cal_model)}
					if(any(objects(envir=as.environment(".GlobalEnv"))=="cal_model")){rm(cal_model,envir=as.environment(".GlobalEnv"))}	
				}
				# add saved calibration model ############################################################
				IS_ID<-isolate(input$Cal_IS_ID)
				target_ID<-isolate(input$Cal_target_ID)
				at_Cal<-isolate(input$Cal_file_set)
				#IS_ID<-"693";target_ID<-"4";at_Cal<-"new_group"
				use_name<-paste("_",IS_ID,"_",target_ID,"_",sep="")
				if(isolate(input$Ion_mode_Cal)=="positive"){	
					if(verbose){cat("\n in M_positive")}
					use_cal<-which(names(cal_models_pos)==at_Cal) # = 1
					if(length(use_cal)==0){stop("\nNo cal_models_pos list available. DEBUG!")} # should not happen at all!
					if(any(names(cal_models_pos[[use_cal]])==use_name)){
						use_entry<-which((names(cal_models_pos[[use_cal]])==use_name))
						if(cal_models_pos[[use_cal]][[use_entry]]$call=="resp ~ 0 + lin"){ # linear, 0-intercept
							abline(
								a=0,
								b=(1/cal_models_pos[[use_cal]][[use_entry]]$coefficients[1]),
								col="gray",lwd=1,lty=2)
						}
						if(cal_models_pos[[use_cal]][[use_entry]]$call=="resp ~ lin"){ # linear, with intercept
							abline(
								a=(-1*cal_models_pos[[use_cal]][[use_entry]]$coefficients[1]/cal_models_pos[[use_cal]][[use_entry]]$coefficients[2]),
								b=(1/cal_models_pos[[use_cal]][[use_entry]]$coefficients[2]),
								col="gray",lwd=1,lty=2)						
						}
						if(cal_models_pos[[use_cal]][[use_entry]]$call=="resp ~ 0 + lin + quad"){ # quadratic, 0-intercept
							for_y<-seq(from=0,to=(max(isolate(dd$d[,4])*2)),length.out=100)
							for_y2<-(for_y^2)							
							for_x<-(
								(cal_models_pos[[use_cal]][[use_entry]]$coefficients[1]*for_y)+
								(cal_models_pos[[use_cal]][[use_entry]]$coefficients[2]*for_y2)
							)
							lines(for_x,for_y,col="gray",lwd=1,lty=2)						
						}
						if(cal_models_pos[[use_cal]][[use_entry]]$call=="resp ~ lin + quad"){ # quadratic, 0-intercept
							for_y<-seq(from=0,to=(max(isolate(dd$d[,4])*2)),length.out=100)
							for_y2<-(for_y^2)							
							for_x<-(
								(cal_models_pos[[use_cal]][[use_entry]]$coefficients[1]) +
								(cal_models_pos[[use_cal]][[use_entry]]$coefficients[2]*for_y)+
								(cal_models_pos[[use_cal]][[use_entry]]$coefficients[3]*for_y2)
							)
							lines(for_x,for_y,col="gray",lwd=1,lty=2)						
						}
						if(cal_models_pos[[use_cal]][[use_entry]]$low_bound!=(-Inf)){
							abline(h=cal_models_pos[[use_cal]][[use_entry]]$low_bound,col="gray",lty=2)
						}
						if(cal_models_pos[[use_cal]][[use_entry]]$up_bound!=Inf){
							abline(h=cal_models_pos[[use_cal]][[use_entry]]$up_bound,col="gray",lty=3)
						}								
						points(
							cal_models_pos[[use_cal]][[use_entry]]$data[,1],
							cal_models_pos[[use_cal]][[use_entry]]$data[,2],
							pch=0,col="gray",cex=2.5
						)
						mtext("Saved calibration model available",side=3,col="gray",line=-1)
					}else{
						mtext("No saved calibration model available",side=3,col="gray",line=-1)
					}
				}
				if(isolate(input$Ion_mode_Cal)=="negative"){
					if(verbose){cat("\n in M_negative")}
					use_cal<-which(names(cal_models_neg)==at_Cal)				
					if(length(use_cal)==0){stop("\nNo cal_models_neg list available. DEBUG!")} # should not happen at all!
					if(any(names(cal_models_neg[[use_cal]])==use_name)){
						if(verbose){cat("\n in M_negative_inside")}
						use_entry<-which((names(cal_models_neg[[use_cal]])==use_name))
						if(cal_models_neg[[use_cal]][[use_entry]]$call=="resp ~ 0 + lin"){ # linear, 0-intercept
							abline(
								a=0,
								b=(1/cal_models_neg[[use_cal]][[use_entry]]$coefficients[1]),
								col="gray",lwd=1,lty=2)
						}
						if(cal_models_neg[[use_cal]][[use_entry]]$call=="resp ~ lin"){ # linear, with intercept
							abline(
								a=(-1*cal_models_neg[[use_cal]][[use_entry]]$coefficients[1]/cal_models_neg[[use_cal]][[use_entry]]$coefficients[2]),
								b=(1/cal_models_neg[[use_cal]][[use_entry]]$coefficients[2]),
								col="gray",lwd=1,lty=2)						
						}
						if(cal_models_neg[[use_cal]][[use_entry]]$call=="resp ~ 0 + lin + quad"){ # quadratic, 0-intercept
							for_y<-seq(from=0,to=(max(isolate(dd$d[,4])*2)),length.out=100)
							for_y2<-(for_y^2)							
							for_x<-(
								(cal_models_neg[[use_cal]][[use_entry]]$coefficients[1]*for_y)+
								(cal_models_neg[[use_cal]][[use_entry]]$coefficients[2]*for_y2)
							)
							lines(for_x,for_y,col="gray",lwd=1,lty=2)						
						}
						if(cal_models_neg[[use_cal]][[use_entry]]$call=="resp ~ lin + quad"){ # quadratic, 0-intercept
							for_y<-seq(from=0,to=(max(isolate(dd$d[,4])*2)),length.out=100)
							for_y2<-(for_y^2)							
							for_x<-(
								(cal_models_neg[[use_cal]][[use_entry]]$coefficients[1]) +
								(cal_models_neg[[use_cal]][[use_entry]]$coefficients[2]*for_y)+
								(cal_models_neg[[use_cal]][[use_entry]]$coefficients[3]*for_y2)
							)
							lines(for_x,for_y,col="gray",lwd=1,lty=2)						
						}
						if(cal_models_neg[[use_cal]][[use_entry]]$low_bound!=(-Inf)){
							abline(h=cal_models_neg[[use_cal]][[use_entry]]$low_bound,col="gray",lty=2)
						}
						if(cal_models_neg[[use_cal]][[use_entry]]$up_bound!=Inf){
							abline(h=cal_models_neg[[use_cal]][[use_entry]]$up_bound,col="gray",lty=3)
						}					
						points(
							cal_models_neg[[use_cal]][[use_entry]]$data[,1],
							cal_models_neg[[use_cal]][[use_entry]]$data[,2],
							pch=0,col="gray",cex=2.5
						)
						mtext("Saved calibration model available",side=3,col="gray",line=-1)
					}else{
						mtext("No saved calibration model available",side=3,col="gray",line=-1)
					}					
				}					
			})
			output$cal_model_summary<-renderText({		
				if(verbose){cat("\n in M_text_output")}
				if(cal_model$call[[2]]=="resp ~ 0 + lin"){ # linear, 0-intercept
					coefi<-cal_model$coefficient
					printthis<-paste("Concentration = ",as.character(round(coefi[[1]],digits=5)),"*Ratio",sep="")
					if(verbose){cat("\n in M_text_1")}
				}
				if(cal_model$call[[2]]=="resp ~ lin"){ # linear, with intercept
					coefi<-cal_model$coefficient
					printthis<-paste("Concentration = ",as.character(round(coefi[[1]],digits=5))," + ",as.character(round(coefi[[2]],digits=5)),"*Ratio",sep="")
					if(verbose){cat("\n in M_text_2")}
				}
				if(cal_model$call[[2]]=="resp ~ 0 + lin + quad"){ # linear, 0-intercept
					coefi<-cal_model$coefficient
					printthis<-paste("Concentration = ",
						as.character(round(coefi[[1]],digits=5)),"*Ratio + ",
						as.character(round(coefi[[2]],digits=5)),"*Ratio^2",					
					sep="")
				}
				if(cal_model$call[[2]]=="resp ~ lin + quad"){ # linear, 0-intercept
					coefi<-cal_model$coefficient
					printthis<-paste("Concentration = ",
						as.character(round(coefi[[1]],digits=5))," + ",
						as.character(round(coefi[[2]],digits=5)),"*Ratio + ",
						as.character(round(coefi[[3]],digits=5)),"*Ratio^2",					
					sep="")
				}				
				return(printthis)			
			})
		}else{
			output$cal_plot <- renderPlot({
				plot(0.5,0.5,col="white",xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,"Not enough calibration data available")
			})
			output$cal_model_summary<-renderText({"Not enough calibration data available"})
			cal_model<<-NA
		}	
	}	
})

###########################################################################################################

###########################################################################################################
# PLOT ZOOM & TABLE SELECTION #############################################################################

# When a double-click happens, check if there's a brush on the plot.
# If so, zoom to the brush bounds; if not, reset the zoom.
observeEvent(input$cal_plot_dblclick, { # - N
	if(verbose){cat("\n in N")}
    brush <- input$cal_plot_brush
    if (!is.null(brush)) {
		ranges_cal_plot$y <- c(brush$xmin, brush$xmax)
		ranges_cal_plot$x <- c(brush$ymin, brush$ymax)
    } else {
		ranges_cal_plot$x <- NULL
		ranges_cal_plot$y <- NULL
    }
})

observeEvent( input$cal_table_rows_selected,{ # - O
	if(verbose){cat("\n in O")}
	cat(paste("\nModified calibration table in row",as.character(input$cal_table_rows_selected)))
	if((isolate(init$a)=="TRUE")&(all(input$cal_table_rows_selected>0))){
		if(dd$d[input$cal_table_rows_selected,8]==1){
			dd$d[input$cal_table_rows_selected,8]<-0
		}else{
			dd$d[input$cal_table_rows_selected,8]<-1
		}
	}
})
###########################################################################################################

###########################################################################################################
# SAVE & REMOVE MODELS ####################################################################################
observe({ # - P
	input$save_Cal
	if(verbose){cat("\n in P")}
	if((isolate(init$a)=="TRUE")){
		if(sum(isolate(dd$d[,8]))>=2){ # not to save am invalid model ...
			if(isolate(input$Ion_mode_Cal)=="positive"){
				IS_ID<-isolate(input$Cal_IS_ID)
				target_ID<-isolate(input$Cal_target_ID)
				at_Cal<-isolate(input$Cal_file_set)
				#IS_ID<-"693";target_ID<-"4";at_Cal<-"A"
				use_cal<-which(names(cal_models_pos)==at_Cal) # well, the first entry ... just in case different calibration groups are merged into a list at some point (= makes saving too slow).
				use_name<-paste("_",IS_ID,"_",target_ID,"_",sep="")
				if(!any(names(cal_models_pos[[use_cal]])==use_name)){ # model not yet existing
					make_entry<-1
					if(length(cal_models_pos[[use_cal]])>0){
						gapped<-FALSE
						for(make_entry in 1:length(cal_models_pos[[use_cal]])){
							if(length(cal_models_pos[[use_cal]][[make_entry]])==0){
								gapped<-TRUE
								break
							}
						}
						if(!gapped){
							make_entry<-(length(cal_models_pos[[use_cal]])+1)
						}
					}
				}else{
					make_entry<-which((names(cal_models_pos[[use_cal]])==use_name))
				}	
				cal_models_pos[[use_cal]][[make_entry]]<<-list()
				cal_models_pos[[use_cal]][[make_entry]]$call<<-cal_model$call[[2]]
				cal_models_pos[[use_cal]][[make_entry]]$coefficients<<-cal_model$coefficients
				cal_models_pos[[use_cal]][[make_entry]]$data<<-cal_model$model
				if(isolate(input$cal_model_weight)){
					cal_models_pos[[use_cal]][[make_entry]]$weighted<<-TRUE	
				}else{
					cal_models_pos[[use_cal]][[make_entry]]$weighted<<-FALSE
				}
				if(isolate(input$cal_model_bound_low)){				
					cal_models_pos[[use_cal]][[make_entry]]$low_bound<<-isolate(input$cal_model_bound_low_value)			
				}else{
					cal_models_pos[[use_cal]][[make_entry]]$low_bound<<-(-Inf)
				}
				if(isolate(input$cal_model_bound_up)){				
					cal_models_pos[[use_cal]][[make_entry]]$up_bound<<-isolate(input$cal_model_bound_up_value)			
				}else{
					cal_models_pos[[use_cal]][[make_entry]]$up_bound<<-(Inf)
				}						
				names(cal_models_pos[[use_cal]])[make_entry]<<-use_name	
				dump("cal_models_pos",file=file.path(logfile[[1]],"quantification",paste("cal_models_pos_",isolate(input$Cal_file_set),sep="")),envir=as.environment(".GlobalEnv"));					
				cat("\n Calibration model saved")
				# count missing models
				at_Cal<-isolate(input$Cal_file_set)
				use_cal<-which(names(cal_models_pos)==at_Cal) # well, the first entry ... just in case different calibration groups are merged into a list at some point (= makes saving too slow).
				if(length(use_cal)>0){
					output$number_missing_models<-renderText({
						paste("Number of Targets with missing calibration models: ",
							as.character(sum(is.na(match(
								paste("_",targets[,"ID_internal_standard"],"_",targets[,"ID"],"_",sep=""),
								names(cal_models_pos[[use_cal]])
							))))
						,sep="")
					})	
				}else{
					output$number_missing_models<-renderText({
						paste("Number of Targets with missing calibration models: ",
							as.character(length(targets[,1]))
						,sep="")
					})				
				}
			}
			if(isolate(input$Ion_mode_Cal)=="negative"){
				IS_ID<-isolate(input$Cal_IS_ID)
				target_ID<-isolate(input$Cal_target_ID)
				at_Cal<-isolate(input$Cal_file_set)
				use_cal<-which(names(cal_models_neg)==at_Cal)
				use_name<-paste("_",IS_ID,"_",target_ID,"_",sep="")
				if(!any(names(cal_models_neg[[use_cal]])==use_name)){ # model not yet existing
					make_entry<-1
					if(length(cal_models_neg[[use_cal]])>0){
						gapped<-FALSE
						for(make_entry in 1:length(cal_models_neg[[use_cal]])){
							if(length(cal_models_neg[[use_cal]][[make_entry]])==0){
								gapped<-TRUE
								break
							}
						}
						if(!gapped){
							make_entry<-(length(cal_models_neg[[use_cal]])+1)
						}
					}
				}else{
					make_entry<-which((names(cal_models_neg[[use_cal]])==use_name))
				}			
				cal_models_neg[[use_cal]][[make_entry]]<<-list()
				cal_models_neg[[use_cal]][[make_entry]]$call<<-cal_model$call[[2]]
				cal_models_neg[[use_cal]][[make_entry]]$coefficients<<-cal_model$coefficients
				cal_models_neg[[use_cal]][[make_entry]]$data<<-cal_model$model
				if(isolate(input$cal_model_weight)){
					cal_models_neg[[use_cal]][[make_entry]]$weighted<<-TRUE	
				}else{
					cal_models_neg[[use_cal]][[make_entry]]$weighted<<-FALSE
				}
				if(isolate(input$cal_model_bound_low)){				
					cal_models_neg[[use_cal]][[make_entry]]$low_bound<<-isolate(input$cal_model_bound_low_value)			
				}else{
					cal_models_neg[[use_cal]][[make_entry]]$low_bound<<-(-Inf)
				}
				if(isolate(input$cal_model_bound_up)){				
					cal_models_neg[[use_cal]][[make_entry]]$up_bound<<-isolate(input$cal_model_bound_up_value)			
				}else{
					cal_models_neg[[use_cal]][[make_entry]]$up_bound<<-(Inf)
				}				
				names(cal_models_neg[[use_cal]])[make_entry]<<-use_name
				dump("cal_models_neg",file=file.path(logfile[[1]],"quantification",paste("cal_models_neg_",isolate(input$Cal_file_set),sep="")),envir=as.environment(".GlobalEnv"));					
				cat("\n Calibration model saved")
				# count missing models
				at_Cal<-isolate(input$Cal_file_set)
				use_cal<-which(names(cal_models_neg)==at_Cal) # well, the first entry ... just in case different calibration groups are merged into a list at some point (= makes saving too slow).
				if(length(use_cal)>0){
					output$number_missing_models<-renderText({
						paste("Number of Targets with missing calibration models: ",
							as.character(sum(is.na(match(
								paste("_",targets[,"ID_internal_standard"],"_",targets[,"ID"],"_",sep=""),
								names(cal_models_neg[[use_cal]])
							))))
						,sep="")
					})	
				}else{
					output$number_missing_models<-renderText({
						paste("Number of Targets with missing calibration models: ",
							as.character(length(targets[,1]))
						,sep="")
					})				
				}
			}
			enviMass:::workflow_set(down="quantification",check_node=TRUE,single_file=FALSE,except="calibration")	
			isolate(redo_cal$a<-(redo_cal$a+1))
		}	
	}
})

observe({ # - Q
	input$remove_Cal
	if(verbose){cat("\n in Q")}
	if((isolate(init$a)=="TRUE")){
		if(isolate(input$Ion_mode_Cal)=="positive"){
			IS_ID<-isolate(input$Cal_IS_ID)
			target_ID<-isolate(input$Cal_target_ID)
			at_Cal<-isolate(input$Cal_file_set)
			use_cal<-which(names(cal_models_pos)==at_Cal)
			if(length(use_cal)>0){
				use_name<-paste("_",IS_ID,"_",target_ID,"_",sep="")
				if(any(names(cal_models_pos[[use_cal]])==use_name)){ # model not yet existing
					delete_entry<-which((names(cal_models_pos[[use_cal]])==use_name))
					cal_models_pos[[use_cal]][[delete_entry]]<<-list()
					names(cal_models_pos[[use_cal]])[delete_entry]<<-"_"
					if(length(cal_models_pos[[use_cal]][[delete_entry]])!=0){
						stop("\n Calibration model delete fucked up. Debug me")
					}
				}
				dump("cal_models_pos",file=file.path(logfile[[1]],"quantification",paste("cal_models_pos_",isolate(input$Cal_file_set),sep="")),envir=as.environment(".GlobalEnv"));					
				# count missing models
				at_Cal<-isolate(input$Cal_file_set)
				use_cal<-which(names(cal_models_pos)==at_Cal) # well, the first entry ... just in case different calibration groups are merged into a list at some point (= makes saving too slow).
				if(length(use_cal)>0){
					output$number_missing_models<-renderText({
						paste("Number of Targets with missing calibration models: ",
							as.character(sum(is.na(match(
								paste("_",targets[,"ID_internal_standard"],"_",targets[,"ID"],"_",sep=""),
								names(cal_models_pos[[use_cal]])
							))))
						,sep="")
					})	
				}else{
					output$number_missing_models<-renderText({
						paste("Number of Targets with missing calibration models: ",
							as.character(length(targets[,1]))
						,sep="")
					})				
				}
			}else{
				cat("\n Nothing to remove ...")
			}
		}
		if(isolate(input$Ion_mode_Cal)=="negative"){
			IS_ID<-isolate(input$Cal_IS_ID)
			target_ID<-isolate(input$Cal_target_ID)
			at_Cal<-isolate(input$Cal_file_set)
			use_cal<-which(names(cal_models_neg)==at_Cal)
			if(length(use_cal)>0){
				use_name<-paste("_",IS_ID,"_",target_ID,"_",sep="")
				if(any(names(cal_models_neg[[use_cal]])==use_name)){ # model not yet existing
					delete_entry<-which((names(cal_models_neg[[use_cal]])==use_name))
					cal_models_neg[[use_cal]][[delete_entry]]<<-list()
					names(cal_models_neg[[use_cal]])[delete_entry]<<-"_"
					if(length(cal_models_neg[[use_cal]][[delete_entry]])!=0){
						stop("\n Calibration model delete fucked up. Debug me")
					}
				}
				dump("cal_models_neg",file=file.path(logfile[[1]],"quantification",paste("cal_models_neg_",isolate(input$Cal_file_set),sep="")));					
				# count missing models
				at_Cal<-isolate(input$Cal_file_set)
				use_cal<-which(names(cal_models_neg)==at_Cal) # well, the first entry ... just in case different calibration groups are merged into a list at some point (= makes saving too slow).
				if(length(use_cal)>0){
					output$number_missing_models<-renderText({
						paste("Number of Targets with missing calibration models: ",
							as.character(sum(is.na(match(
								paste("_",targets[,"ID_internal_standard"],"_",targets[,"ID"],"_",sep=""),
								names(cal_models_neg[[use_cal]])
							))))
						,sep="")
					})	
				}else{
					output$number_missing_models<-renderText({
						paste("Number of Targets with missing calibration models: ",
							as.character(length(targets[,1]))
						,sep="")
					})				
				}
			}else{
				cat("\n Nothing to remove ...")
			}			
		}
		enviMass:::workflow_set(down="quantification",check_node=TRUE,single_file=FALSE,except="calibration")	
		isolate(redo_cal$a<-(redo_cal$a+1))
	}
})
###########################################################################################################

###########################################################################################################
# SET NEW IS AS DEFAULT FOR QUANTIFICATION ################################################################
observe({
	input$use_Cal
	if(verbose){cat("\n in R")}
	if(
		isolate(input$Cal_IS_ID!="none") & isolate(input$Cal_target_ID!="none") & 
		isolate(input$Cal_IS_name!="none") & isolate(input$Cal_target_name!="none") 					
	){
		change_IS_target_ID<-isolate(input$Cal_target_ID)
		change_IS_IS_ID<-isolate(input$Cal_IS_ID)
		this<-which(targets[,1]==change_IS_target_ID)
		targets[this,6]<-as.character(change_IS_IS_ID)
		targets[this,6]<<-as.character(change_IS_IS_ID)
		write.table(targets,file=file.path(logfile[[1]],"dataframes","targets.txt"),row.names=FALSE,sep="\t",quote=FALSE)
		output$targets<<-DT::renderDataTable(read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character"));      
		cat("\n Changed/set default IS to be used in quantification for the selected target.")
		if(isolate(input$Ion_mode_Cal)=="positive"){		
			# count missing models
			at_Cal<-isolate(input$Cal_file_set)
			use_cal<-which(names(cal_models_pos)==at_Cal) # well, the first entry ... just in case different calibration groups are merged into a list at some point (= makes saving too slow).
			if(length(use_cal)>0){
				output$number_missing_models<-renderText({
					paste("Number of Targets with missing calibration models: ",
						as.character(sum(is.na(match(
							paste("_",targets[,"ID_internal_standard"],"_",targets[,"ID"],"_",sep=""),
							names(cal_models_pos[[use_cal]])
						))))
					,sep="")
				})	
			}else{
				output$number_missing_models<-renderText({
					paste("Number of Targets with missing calibration models: ",
						as.character(length(targets[,1]))
					,sep="")
				})				
			}
		}
		if(isolate(input$Ion_mode_Cal)=="negative"){		
			# count missing models
			at_Cal<-isolate(input$Cal_file_set)
			use_cal<-which(names(cal_models_neg)==at_Cal) # well, the first entry ... just in case different calibration groups are merged into a list at some point (= makes saving too slow).
			if(length(use_cal)>0){
				output$number_missing_models<-renderText({
					paste("Number of Targets with missing calibration models: ",
						as.character(sum(is.na(match(
							paste("_",targets[,"ID_internal_standard"],"_",targets[,"ID"],"_",sep=""),
							names(cal_models_neg[[use_cal]])
						))))
					,sep="")
				})	
			}else{
				output$number_missing_models<-renderText({
					paste("Number of Targets with missing calibration models: ",
						as.character(length(targets[,1]))
					,sep="")
				})				
			}
		}	
	}
})
###########################################################################################################

if(any(ls()=="logfile")){stop("\n illegal logfile detected in server_obs_calibration.r!")}
if(any(ls()=="cal_models_neg")){stop("\n illegal cal_models_neg detected in server_obs_calibration.r!")}
if(any(ls()=="cal_models_pos")){stop("\n illegal cal_models_pos detected in server_obs_calibration.r!")} 
 
 
 
 
 
 
 
 
 
