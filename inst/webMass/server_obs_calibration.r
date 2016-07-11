if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_obs_screening.r!")}
verbose<-TRUE
###########################################################################################################
# SPECIFY IONIZATION MODE #################################################################################
observe({ 
	input$Ion_mode_Cal 
	if(verbose){cat("\n in A")}
	if(isolate(init$a)=="TRUE"){
	if(logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="calibration"]!="TRUE"){ # e.g., if files were changed / added / ...
		if(verbose){cat("\n in A_1")}
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
	}else{
		info("Calibration files have been added or modified or copied. Workflow recalculation including the calibration required.");
		cat("\n Calibration files have been added or modified or copied. Recalculation required!")
	}	
	}
})
###########################################################################################################

###########################################################################################################
# SPECIFY CALIBRATION GROUP ###############################################################################
observe({ 
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
				# load & update available calibration models
				if(!file.exists(file.path(logfile[[1]],"quantification",paste("cal_models_pos_",isolate(input$Cal_file_set),sep="")))){
					cal_models_pos<<-list()
					cal_models_pos[[1]]<<-list()
					names(cal_models_pos)<<-isolate(input$Cal_file_set)
					save(cal_models_pos,file=file.path(logfile[[1]],"quantification",paste("cal_models_pos_",isolate(input$Cal_file_set),sep="")),envir=as.environment(".GlobalEnv"));			
				}
				load(file=file.path(logfile[[1]],"quantification",paste("cal_models_pos_",isolate(input$Cal_file_set),sep="")),envir=as.environment(".GlobalEnv"));	
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
				# load & update available calibration models
				if(!file.exists(file.path(logfile[[1]],"quantification",paste("cal_models_neg_",isolate(input$Cal_file_set),sep="")))){
					cal_models_neg<<-list()
					cal_models_neg[[1]]<<-list()
					names(cal_models_neg)<<-isolate(input$Cal_file_set)
					save(cal_models_neg,file=file.path(logfile[[1]],"quantification",paste("cal_models_neg_",isolate(input$Cal_file_set),sep="")),envir=as.environment(".GlobalEnv"));			
				}
				load(file=file.path(logfile[[1]],"quantification",paste("cal_models_neg_",isolate(input$Cal_file_set),sep="")),envir=as.environment(".GlobalEnv"));	
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
	input$Cal_target_ID
	if(verbose){cat("\n in C")}
	if((isolate(init$a)=="TRUE")){
		if(isolate(input$Cal_target_ID)!="none"){
			use_this_name<-targets[targets[,1]==isolate(input$Cal_target_ID),2]
			updateSelectInput(session,inputId="Cal_target_name",selected = as.character(use_this_name))
			new_IS_ID<-targets[targets[,1]==isolate(input$Cal_target_ID),6]
			if(verbose){cat(" -p- ");cat(new_IS_ID)}
			if(length(new_IS_ID)>0){ # just in case ...
				if(new_IS_ID!="FALSE"){
					updateSelectInput(session,inputId="Cal_IS_ID",selected = as.character(new_IS_ID))
				}
			}
		}else{	
			updateSelectInput(session,inputId="Cal_target_name",selected="none")
		}
	}	
})

observe({ 
	input$Cal_target_name
	if(verbose){cat("\n in D")}
	if((isolate(init$a)=="TRUE")&(isolate(input$Cal_target_name)!="none")){
		use_this_ID<-use_this_name<-targets[targets[,2]==isolate(input$Cal_target_name),1]
		updateSelectInput(session,inputId="Cal_target_ID",selected = as.character(use_this_ID))
	}	
})

observe({ 
	input$Cal_IS_ID
	if(verbose){cat("\n in E")}
	if((isolate(init$a)=="TRUE")){
		if(isolate(input$Cal_IS_ID)!="none"){
			use_this_name<-intstand[intstand[,1]==isolate(input$Cal_IS_ID),2]
			updateSelectInput(session,inputId="Cal_IS_name",selected = as.character(use_this_name))
		}else{
			updateSelectInput(session,inputId="Cal_IS_name",selected="none")
		}
	}
})

observe({ 
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
observe({ 
	input$Cal_next
	if(verbose){cat("\n in G")}
	if((isolate(init$a)=="TRUE") & (isolate(input$Cal_file_set)!="none")){
		is_at_targetID<-(isolate(input$Cal_target_ID))
		is_at_ISID<-(isolate(input$Cal_IS_ID))
		in_table<-which( ((targets[,1]==is_at_targetID)&(targets[,6]==is_at_ISID)) )
		at_target_ID<-"none"
		at_IS_ID<-"none"
		if( (length(in_table)==0) || (is.na(in_table[1])) ){ # match with is_at_ISID not existing in table; reset to first available entry
			if(verbose){cat("\n in G_1")}
			in_table<-which(targets[,1]==is_at_targetID)
			in_table<-in_table[1]
			if( (length(in_table)==0) || (is.na(in_table[1])) ){ # for any invalid entry, start from beginning
				if(verbose){cat("\n in G_1_1")}
				in_table<-1
			}
			at_target_ID<-as.character(targets[in_table,1])
			at_IS_ID<-as.character(targets[in_table,6])	
		}else{ # match existing; get next entry
			if(verbose){cat("\n in G_2")}
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
	input$Cal_previous
	if(verbose){cat("\n in H")}
	if((isolate(init$a)=="TRUE")&(isolate(input$Cal_target_name)!="none")& (isolate(input$Cal_file_set)!="none")){
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
	input$Cal_first
	if(verbose){cat("\n in I")}
	if((isolate(init$a)=="TRUE")& (isolate(input$Cal_file_set)!="none")){
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
	input$Cal_last
	if(verbose){cat("\n in J")}
	if((isolate(init$a)=="TRUE")& (isolate(input$Cal_file_set)!="none")){
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
# RETRIEVE SETS ###########################################################################################
ranges_cal_plot <- reactiveValues(x = NULL, y = NULL)
dd	<-	reactiveValues()
observe({ 
	input$Cal_IS_ID
	input$Cal_IS_name
	input$Cal_target_ID
	input$Cal_target_name
	if(verbose){cat("\n in K")}
	# anything available from calibration screening? Results still in measurements?
	nothing<-TRUE
	if(isolate(input$Ion_mode_Cal)=="positive"){ # check if compounds in tables are also found in screening outcomes!
		if(
			any(as.character(results_screen_target_pos_cal[[1]][,1])==isolate(input$Cal_target_ID)) &
			any(as.character(results_screen_IS_pos_cal[[1]][,1])==isolate(input$Cal_IS_ID)) 
		){
			nothing<-FALSE
		}
	}
	if(isolate(input$Ion_mode_Cal)=="negative"){ # check if compounds in tables are also found in screening outcomes!
		if(
			any(as.character(results_screen_target_neg_cal[[1]][,1])==isolate(input$Cal_target_ID)) &
			any(as.character(results_screen_IS_neg_cal[[1]][,1])==isolate(input$Cal_IS_ID))
		){
			nothing<-FALSE
		}
	}	
	if(verbose){cat(nothing)}
	if(
		(isolate(init$a)=="TRUE") & 
		(!nothing) & 
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
			#target_ID<-"4"; IS_ID<-"693";at_Cal<-"A"			
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
									measurements[measurements[,1]==target_in_file[i],]$tag1	
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
			mat_cal<-mat_cal[!duplicated(mat_cal),,drop=FALSE] # same peaks in different combinations - remove
			mat_cal[,1]<-(1:length(mat_cal[,1]))
			min_int<-as.numeric(intstand[intstand[,1]==IS_ID,17])
			max_int<-as.numeric(intstand[intstand[,1]==IS_ID,18])
			mat_cal[log10(mat_cal[,2])<min_int,8]<-0
			mat_cal[log10(mat_cal[,2])>max_int,8]<-0
			mat_cal<<-mat_cal
			dd$d<-mat_cal
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
						if(measurements[measurements[,1]==res_IS_neg_screen_cal[[at_entry]][[j]][[1]]$file_ID,]$tag2==at_Cal){
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
						if(measurements[measurements[,1]==res_target_neg_screen_cal[[at_entry]][[j]][[1]]$file_ID,]$tag2==at_Cal){
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
									measurements[measurements[,1]==target_in_file[i],]$tag1	
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
			mat_cal<-mat_cal[!duplicated(mat_cal),,drop=FALSE] # same peaks in different combinations - remove
			mat_cal[,1]<-(1:length(mat_cal[,1]))
			min_int<-as.numeric(intstand[intstand[,1]==IS_ID,17])
			max_int<-as.numeric(intstand[intstand[,1]==IS_ID,18])
			mat_cal[mat_cal[,2]<min_int,8]<-0
			mat_cal[mat_cal[,2]>max_int,8]<-0
			mat_cal<<-mat_cal
			dd$d<-mat_cal
		}
	}else{
		dd$d<-matrix(ncol=2,nrow=0,0)
	}	
})				
###########################################################################################################				

###########################################################################################################
# PLOT SETS & OUTPUT TABLE ################################################################################

redo_cal<-reactiveValues() 	# reactive value to indicate model save / removal
redo_cal$a<-1
observe({ 
	input$Cal_IS_ID
	input$Cal_IS_name
	input$Cal_target_ID
	input$Cal_target_name
	redo_cal$a
	if(verbose){cat("\n in L")}
	if((isolate(init$a)=="TRUE") & (isolate(input$Cal_IS_ID)!="none")  & (isolate(input$Cal_target_ID)!="none")& (isolate(input$Cal_file_set)!="none")){					
		# generate outputs ######################################################
		if(length(dd$d[,1])>0){
			output$cal_table <- DT::renderDataTable(
				datatable(
					dd$d,selection =c('single'),options = list(lengthMenu = c(25,50,100))
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

observe({ 
	input$Cal_IS_ID
	input$Cal_IS_name
	input$Cal_target_ID
	input$Cal_target_name
	input$cal_model
	input$cal_model_0intercept
	redo_cal$a
	if(verbose){cat("\n in M")}
	if((isolate(init$a)=="TRUE") & (isolate(input$Cal_IS_ID)!="none") & (isolate(input$Cal_target_ID)!="none")& (isolate(input$Cal_file_set)!="none")){					
		# generate outputs ######################################################
		if(length(dd$d[,1])>1){
			output$cal_plot <- renderPlot({
				plot(dd$d[,5], dd$d[,4],
					xlab="Concentration",ylab="Intensity ratio",pch=19,
					xlim=ranges_cal_plot$y,ylim=ranges_cal_plot$x,
					main="Draw rectangles and double-click into them to zoom, double-click again to zoom out.",cex.main=1,col="white")#,yaxs="i",xaxs="i")
				abline(h=0,col="lightgrey")
				abline(v=0,col="lightgrey")				
				points(dd$d[dd$d[,8]==1,5], dd$d[dd$d[,8]==1,4],col="black",pch=19)
				points(dd$d[dd$d[,8]==0,5], dd$d[dd$d[,8]==0,4],col="gray",pch=19)				
				if((!is.null(ranges_cal_plot$x))||(!is.null(ranges_cal_plot$y))){
					mtext("Now zoomed in",side=3,col="gray")
				}
				if(sum(dd$d[,8])>=2){ # at least two data points!
					lin<-(dd$d[dd$d[,8]==1,5])
					resp<-(dd$d[dd$d[,8]==1,4])
					if(isolate(input$cal_model)=="linear"){
						if(isolate(input$cal_model_0intercept)){
							cal_model<<-lm(resp~0+lin)
						}else{
							cal_model<<-lm(resp~lin)
						}
						abline(cal_model,col="red",lwd=2)
					}else{
						quad<-((dd$d[dd$d[,8]==1,5])^2)
						if(isolate(input$cal_model_0intercept)){
							cal_model<<-lm(resp~0+lin+quad)					
						}else{
							cal_model<<-lm(resp~lin+quad)							
						}
						for_x<-seq(from=0,to=(max(dd$d[,5])*2),length.out=100)
						for_x2<-(for_x^2)
						for_y<-predict(cal_model,list(lin=for_x,quad=for_x2))
						lines(for_x,for_y,col="red",lwd=2)
					}
				}
				# add saved calibration model ############################################################
				IS_ID<-isolate(input$Cal_IS_ID)
				target_ID<-isolate(input$Cal_target_ID)
				at_Cal<-isolate(input$Cal_file_set)
				#IS_ID<-"693";target_ID<-"4";at_Cal<-"new_group"
				use_name<-paste(IS_ID,target_ID,sep="_")
				if(isolate(input$Ion_mode_Cal)=="positive"){	
					if(verbose){cat("\n in M_positive")}
					use_cal<-which(names(cal_models_pos)==at_Cal) # = 1
					if(any(names(cal_models_pos[[use_cal]])==use_name)){
						use_entry<-which((names(cal_models_pos[[use_cal]])==use_name))
						if(cal_models_pos[[use_cal]][[use_entry]]$call=="resp ~ 0 + lin"){ # linear, 0-intercept
							abline(
								a=0,
								b=cal_models_pos[[use_cal]][[use_entry]]$coefficients,
								col="gray",lwd=1,lty=2)
						}
						if(cal_models_pos[[use_cal]][[use_entry]]$call=="resp ~ lin"){ # linear, with intercept
							abline(
								a=cal_models_pos[[use_cal]][[use_entry]]$coefficients[1],
								b=cal_models_pos[[use_cal]][[use_entry]]$coefficients[2],
								col="gray",lwd=1,lty=2)						
						}
						if(cal_models_pos[[use_cal]][[use_entry]]$call=="resp ~ 0 + lin + quad"){ # quadratic, 0-intercept
							for_x<-seq(from=0,to=(max(dd$d[,5]*2)),length.out=100)
							for_x2<-(for_x^2)							
							for_y<-(
								(cal_models_pos[[use_cal]][[use_entry]]$coefficients[1]*for_x)+
								(cal_models_pos[[use_cal]][[use_entry]]$coefficients[2]*for_x2)
							)
							lines(for_x,for_y,col="gray",lwd=1,lty=2)						
						}
						if(cal_models_pos[[use_cal]][[use_entry]]$call=="resp ~ lin + quad"){ # quadratic, 0-intercept
							for_x<-seq(from=0,to=(max(dd$d[,5]*2)),length.out=100)
							for_x2<-(for_x^2)							
							for_y<-(
								(cal_models_pos[[use_cal]][[use_entry]]$coefficients[1]) +
								(cal_models_pos[[use_cal]][[use_entry]]$coefficients[2]*for_x)+
								(cal_models_pos[[use_cal]][[use_entry]]$coefficients[3]*for_x2)
							)
							lines(for_x,for_y,col="gray",lwd=1,lty=2)						
						}
						points(
							cal_models_pos[[use_cal]][[use_entry]]$data[,2],
							cal_models_pos[[use_cal]][[use_entry]]$data[,1],
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
if(length(use_cal)!=0){
					if(any(names(cal_models_neg[[use_cal]])==use_name)){
						if(verbose){cat("\n in M_negative_inside")}
						use_entry<-which((names(cal_models_neg[[use_cal]])==use_name))
						if(cal_models_neg[[use_cal]][[use_entry]]$call=="resp ~ 0 + lin"){ # linear, 0-intercept
							abline(
								a=0,
								b=cal_models_neg[[use_cal]][[use_entry]]$coefficients,
								col="gray",lwd=1,lty=2)
						}
						if(cal_models_neg[[use_cal]][[use_entry]]$call=="resp ~ lin"){ # linear, with intercept
							abline(
								a=cal_models_neg[[use_cal]][[use_entry]]$coefficients[1],
								b=cal_models_neg[[use_cal]][[use_entry]]$coefficients[2],
								col="gray",lwd=1,lty=2)						
						}
						if(cal_models_neg[[use_cal]][[use_entry]]$call=="resp ~ 0 + lin + quad"){ # quadratic, 0-intercept
							for_x<-seq(from=0,to=(max(dd$d[,5]*2)),length.out=100)
							for_x2<-(for_x^2)							
							for_y<-(
								(cal_models_neg[[use_cal]][[use_entry]]$coefficients[1]*for_x)+
								(cal_models_neg[[use_cal]][[use_entry]]$coefficients[2]*for_x2)
							)
							lines(for_x,for_y,col="gray",lwd=1,lty=2)						
						}
						if(cal_models_neg[[use_cal]][[use_entry]]$call=="resp ~ lin + quad"){ # quadratic, 0-intercept
							for_x<-seq(from=0,to=(max(dd$d[,5]*2)),length.out=100)
							for_x2<-(for_x^2)							
							for_y<-(
								(cal_models_neg[[use_cal]][[use_entry]]$coefficients[1]) +
								(cal_models_neg[[use_cal]][[use_entry]]$coefficients[2]*for_x)+
								(cal_models_neg[[use_cal]][[use_entry]]$coefficients[3]*for_x2)
							)
							lines(for_x,for_y,col="gray",lwd=1,lty=2)						
						}
						points(
							cal_models_neg[[use_cal]][[use_entry]]$data[,2],
							cal_models_neg[[use_cal]][[use_entry]]$data[,1],
							pch=0,col="gray",cex=2.5
						)
						mtext("Saved calibration model available",side=3,col="gray",line=-1)
					}else{
						mtext("No saved calibration model available",side=3,col="gray",line=-1)
					}					
				}		
}else{
	mtext("No saved calibration model available",side=3,col="gray",line=-1)
}				
			})
			output$cal_model_summary<-renderText({		
				if(verbose){cat("\n in M_text_output")}
				if(cal_model$call[[2]]=="resp ~ 0 + lin"){ # linear, 0-intercept
					coefi<-cal_model$coefficient
					printthis<-paste("Ratio = ",as.character(round(coefi[[1]],digits=5)),"*Concentration",sep="")
					if(verbose){cat("\n in M_text_1")}
				}
				if(cal_model$call[[2]]=="resp ~ lin"){ # linear, with intercept
					coefi<-cal_model$coefficient
					printthis<-paste("Ratio = ",as.character(round(coefi[[1]],digits=5))," + ",as.character(round(coefi[[2]],digits=5)),"*Concentration",sep="")
					if(verbose){cat("\n in M_text_2")}
				}
				if(cal_model$call[[2]]=="resp ~ 0 + lin + quad"){ # linear, 0-intercept
					coefi<-cal_model$coefficient
					printthis<-paste("Ratio = ",
						as.character(round(coefi[[1]],digits=5)),"*Concentration + ",
						as.character(round(coefi[[2]],digits=5)),"*Concentration^2",					
					sep="")
				}
				if(cal_model$call[[2]]=="resp ~ lin + quad"){ # linear, 0-intercept
					coefi<-cal_model$coefficient
					printthis<-paste("Ratio = ",
						as.character(round(coefi[[1]],digits=5))," + ",
						as.character(round(coefi[[2]],digits=5)),"*Concentration + ",
						as.character(round(coefi[[3]],digits=5)),"*Concentration^2",					
					sep="")
				}				
				return(printthis)			
			})
			}else{
			output$cal_plot <- renderPlot({
				plot(0.5,0.5,col="white",xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,"Not enough calibration data available")
			})
		}	
	}	
})

###########################################################################################################

###########################################################################################################
# PLOT ZOOM & TABLE SELECTION #############################################################################

# When a double-click happens, check if there's a brush on the plot.
# If so, zoom to the brush bounds; if not, reset the zoom.
ranges_cal_plot <- reactiveValues(x = NULL, y = NULL)
observeEvent(input$cal_plot_dblclick, {
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

observeEvent(
	input$cal_table_rows_selected,{
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
observe({ 
	input$save_Cal
	if(verbose){cat("\n in P")}
	if((isolate(init$a)=="TRUE")){
		if(isolate(input$Ion_mode_Cal)=="positive"){
			IS_ID<-isolate(input$Cal_IS_ID)
			target_ID<-isolate(input$Cal_target_ID)
			at_Cal<-isolate(input$Cal_file_set)
			#IS_ID<-"693";target_ID<-"4";at_Cal<-"A"
			use_cal<-which(names(cal_models_pos)==at_Cal) # well, the first entry ... just in case different calibration groups are merged into a list at some point (= makes saving too slow).
			use_name<-paste(IS_ID,target_ID,sep="_")
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
			names(cal_models_pos[[use_cal]])[make_entry]<<-use_name	
			save(cal_models_pos,file=file.path(logfile[[1]],"quantification",paste("cal_models_pos_",isolate(input$Cal_file_set),sep="")),envir=as.environment(".GlobalEnv"));					
			cat("\n Calibration model saved")
		}
		if(isolate(input$Ion_mode_Cal)=="negative"){
			IS_ID<-isolate(input$Cal_IS_ID)
			target_ID<-isolate(input$Cal_target_ID)
			at_Cal<-isolate(input$Cal_file_set)
			use_cal<-which(names(cal_models_neg)==at_Cal)
			use_name<-paste(IS_ID,target_ID,sep="_")
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
			assign("cal_model_saved",cal_model)
			cal_models_neg[[use_cal]][[make_entry]]<<-list()
			cal_models_neg[[use_cal]][[make_entry]]$call<<-cal_model$call[[2]]
			cal_models_neg[[use_cal]][[make_entry]]$coefficients<<-cal_model$coefficients
			cal_models_neg[[use_cal]][[make_entry]]$data<<-cal_model$model
			names(cal_models_neg[[use_cal]])[make_entry]<<-use_name
			save(cal_models_neg,file=file.path(logfile[[1]],"quantification",paste("cal_models_neg_",isolate(input$Cal_file_set),sep="")),envir=as.environment(".GlobalEnv"));					
			cat("\n Calibration model saved")
		}
		enviMass:::workflow_set(down="quantification",check_node=TRUE,single_file=FALSE)	
		isolate(redo_cal$a<-(redo_cal$a+1))
	}
})

observe({ 
	input$remove_Cal
	if(verbose){cat("\n in Q")}
	if((isolate(init$a)=="TRUE")){
		if(isolate(input$Ion_mode_Cal)=="positive"){
			IS_ID<-isolate(input$Cal_IS_ID)
			target_ID<-isolate(input$Cal_target_ID)
			at_Cal<-isolate(input$Cal_file_set)
			use_cal<-which(names(cal_models_pos)==at_Cal)
			use_name<-paste(IS_ID,target_ID,sep="_")
			if(any(names(cal_models_pos[[use_cal]])==use_name)){ # model not yet existing
				delete_entry<-which((names(cal_models_pos[[use_cal]])==use_name))
				cal_models_pos[[use_cal]][[delete_entry]]<<-list()
				names(cal_models_pos[[use_cal]])[delete_entry]<<-"_"
				if(length(cal_models_pos[[use_cal]][[delete_entry]])!=0){
					stop("\n Calibration model delete fucked up. Debug me")
				}
			}
			save(cal_models_pos,file=file.path(logfile[[1]],"quantification",paste("cal_models_pos_",isolate(input$Cal_file_set),sep="")),envir=as.environment(".GlobalEnv"));					
		}
		if(isolate(input$Ion_mode_Cal)=="negative"){
			IS_ID<-isolate(input$Cal_IS_ID)
			target_ID<-isolate(input$Cal_target_ID)
			at_Cal<-isolate(input$Cal_file_set)
			use_cal<-which(names(cal_models_neg)==at_Cal)
			use_name<-paste(IS_ID,target_ID,sep="_")
			if(any(names(cal_models_neg[[use_cal]])==use_name)){ # model not yet existing
				delete_entry<-which((names(cal_models_neg[[use_cal]])==use_name))
				cal_models_neg[[use_cal]][[delete_entry]]<<-list()
				names(cal_models_neg[[use_cal]])[delete_entry]<<-"_"
				if(length(cal_models_neg[[use_cal]][[delete_entry]])!=0){
					stop("\n Calibration model delete fucked up. Debug me")
				}
			}
			save(cal_models_neg,file=file.path(logfile[[1]],"quantification",paste("cal_models_neg_",isolate(input$Cal_file_set),sep="")),envir=as.environment(".GlobalEnv"));					
		}
		enviMass:::workflow_set(down="quantification",check_node=TRUE,single_file=FALSE)	
		isolate(redo_cal$a<-(redo_cal$a+1))
	}
})
###########################################################################################################

if(any(ls()=="logfile")){stop("\n illegal logfile detected in server_obs_calibration.r!")}
if(any(ls()=="cal_models_neg")){stop("\n illegal cal_models_neg detected in server_obs_calibration.r!")}
if(any(ls()=="cal_models_pos")){stop("\n illegal cal_models_pos detected in server_obs_calibration.r!")} 
 
 
 
 
 
 
 
 
 