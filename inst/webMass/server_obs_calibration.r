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
				# Update IS compounds
				load(file=file.path(logfile[[1]],"quantification","results_screen_IS_pos_cal"),envir=as.environment(".GlobalEnv"));	
				IS_names<-unique(results_screen_IS_pos_cal[[1]][
					results_screen_IS_pos_cal[[1]][,4]!=0
				,2])
				IS_names<-IS_names[order(IS_names)]
				IS_names<-c("none",IS_names)
				updateSelectInput(session,inputId="Cal_IS_name",label="Name",choices=IS_names,selected = IS_names[1])
				IS_IDs<-unique(results_screen_IS_pos_cal[[1]][
					results_screen_IS_pos_cal[[1]][,4]!=0
				,1])		
				IS_IDs<-IS_IDs[order(IS_IDs)]
				IS_IDs<-c("none",IS_IDs)
				updateSelectInput(session,inputId="Cal_IS_ID",label="Internal standard ID",choices=IS_IDs,selected = IS_IDs[1])
				# Update target compounds
				load(file=file.path(logfile[[1]],"quantification","results_screen_target_pos_cal"),envir=as.environment(".GlobalEnv"));	
				target_names<-unique(results_screen_target_pos_cal[[1]][
					results_screen_target_pos_cal[[1]][,4]!=0
				,2])
				target_names<-target_names[order(target_names)]
				target_names<-c("none",target_names)
				updateSelectInput(session,inputId="Cal_target_name",label="Name",choices=target_names,selected = target_names[1])
				target_IDs<-unique(results_screen_target_pos_cal[[1]][
					results_screen_target_pos_cal[[1]][,4]!=0
				,1])		
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
				# Update IS compounds
				load(file=file.path(logfile[[1]],"quantification","results_screen_IS_neg_cal"),envir=as.environment(".GlobalEnv"));	
				IS_names<-unique(results_screen_IS_neg_cal[[1]][
					results_screen_IS_neg_cal[[1]][,4]!=0
				,2])
				IS_names<-IS_names[order(IS_names)]
				IS_names<-c("none",IS_names)
				updateSelectInput(session,inputId="Cal_IS_name",label="Name",choices=IS_names,selected = IS_names[1])
				IS_IDs<-unique(results_screen_IS_neg_cal[[1]][
					results_screen_IS_neg_cal[[1]][,4]!=0
				,1])		
				IS_IDs<-IS_IDs[order(IS_IDs)]
				IS_IDs<-c("none",IS_IDs)
				updateSelectInput(session,inputId="Cal_IS_ID",label="Internal standard ID",choices=IS_IDs,selected = IS_IDs[1])
				# Update target compounds
				load(file=file.path(logfile[[1]],"quantification","results_screen_target_neg_cal"),envir=as.environment(".GlobalEnv"));	
				target_names<-unique(results_screen_target_neg_cal[[1]][
					results_screen_target_neg_cal[[1]][,4]!=0
				,2])
				target_names<-target_names[order(target_names)]
				target_names<-c("none",target_names)
				updateSelectInput(session,inputId="Cal_target_name",label="Name",choices=target_names,selected = target_names[1])
				target_IDs<-unique(results_screen_target_neg_cal[[1]][
					results_screen_target_neg_cal[[1]][,4]!=0
				,1])		
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
	if((isolate(init$a)=="TRUE")&(isolate(input$Cal_IS_ID)!="none")){
			if(isolate(input$Ion_mode_Cal)=="positive"){
				load(file=file.path(logfile[[1]],"quantification","results_screen_IS_pos_cal"),envir=as.environment(".GlobalEnv"));	
				use_this_name<-unique(results_screen_IS_pos_cal[[1]][
					results_screen_IS_pos_cal[[1]][,1]==isolate(input$Cal_IS_ID)
				,2,drop=FALSE])
				updateSelectInput(session,inputId="Cal_IS_name",selected = use_this_name)
			}
			if(isolate(input$Ion_mode_Cal)=="negative"){
				load(file=file.path(logfile[[1]],"quantification","results_screen_IS_neg_cal"),envir=as.environment(".GlobalEnv"));	
				use_this_name<-unique(results_screen_IS_neg_cal[[1]][
					results_screen_IS_neg_cal[[1]][,1]==isolate(input$Cal_IS_ID)
				,2,drop=FALSE])
				updateSelectInput(session,inputId="Cal_IS_name",selected = use_this_name)
			}
	}	
})

observe({ 
	input$Cal_IS_name
	init$b
	if((isolate(init$a)=="TRUE")&(isolate(input$Cal_IS_name)!="none")){
			if(isolate(input$Ion_mode_Cal)=="positive"){
				load(file=file.path(logfile[[1]],"quantification","results_screen_IS_pos_cal"),envir=as.environment(".GlobalEnv"));	
				use_this_ID<-unique(results_screen_IS_pos_cal[[1]][
					results_screen_IS_pos_cal[[1]][,2]==isolate(input$Cal_IS_name)
				,1,drop=FALSE])
				updateSelectInput(session,inputId="Cal_IS_ID",selected = use_this_ID)
			}
			if(isolate(input$Ion_mode_Cal)=="negative"){
				load(file=file.path(logfile[[1]],"quantification","results_screen_IS_neg_cal"),envir=as.environment(".GlobalEnv"));	
				use_this_ID<-unique(results_screen_IS_neg_cal[[1]][
					results_screen_IS_neg_cal[[1]][,2]==isolate(input$Cal_IS_name)
				,1,drop=FALSE])
				updateSelectInput(session,inputId="Cal_IS_ID",selected = use_this_ID)
			}
	}	
})

observe({ 
	input$Cal_target_ID
	init$b
	if((isolate(init$a)=="TRUE")&(isolate(input$Cal_target_ID)!="none")){
			if(isolate(input$Ion_mode_Cal)=="positive"){
				load(file=file.path(logfile[[1]],"quantification","results_screen_target_pos_cal"),envir=as.environment(".GlobalEnv"));	
				use_this_name<-unique(results_screen_target_pos_cal[[1]][
					results_screen_target_pos_cal[[1]][,1]==isolate(input$Cal_target_ID)
				,2,drop=FALSE])
				updateSelectInput(session,inputId="Cal_target_name",selected = use_this_name)
			}
			if(isolate(input$Ion_mode_Cal)=="negative"){
				load(file=file.path(logfile[[1]],"quantification","results_screen_target_neg_cal"),envir=as.environment(".GlobalEnv"));	
				use_this_name<-unique(results_screen_target_neg_cal[[1]][
					results_screen_target_neg_cal[[1]][,1]==isolate(input$Cal_target_ID)
				,2,drop=FALSE])
				updateSelectInput(session,inputId="Cal_target_name",selected = use_this_name)
			}
	}	
})

observe({ 
	input$Cal_target_name
	init$b
	if((isolate(init$a)=="TRUE")&(isolate(input$Cal_target_name)!="none")){
			if(isolate(input$Ion_mode_Cal)=="positive"){
				load(file=file.path(logfile[[1]],"quantification","results_screen_target_pos_cal"),envir=as.environment(".GlobalEnv"));	
				use_this_ID<-unique(results_screen_target_pos_cal[[1]][
					results_screen_target_pos_cal[[1]][,2]==isolate(input$Cal_target_name)
				,1,drop=FALSE])
				updateSelectInput(session,inputId="Cal_target_ID",selected = use_this_ID)
			}
			if(isolate(input$Ion_mode_Cal)=="negative"){
				load(file=file.path(logfile[[1]],"quantification","results_screen_IS_neg_cal"),envir=as.environment(".GlobalEnv"));	
				use_this_ID<-unique(results_screen_target_neg_cal[[1]][
					results_screen_target_neg_cal[[1]][,2]==isolate(input$Cal_target_name)
				,1,drop=FALSE])
				updateSelectInput(session,inputId="Cal_target_ID",selected = use_this_ID)
			}
	}	
})
###########################################################################################################

###########################################################################################################
# GET ADDUCTS AND ISOTOPOLOGUE PEAKS ######################################################################
# IS COMPOUNDS ############################################################################################
observe({ 
	init$b
	input$Cal_IS_name
	input$Cal_IS_ID
	if((isolate(init$a)=="TRUE") & (isolate(input$Cal_IS_ID)!="none") & (isolate(input$Cal_IS_name)!="none")){
		if(isolate(input$Ion_mode_Cal)=="positive"){
			# update adduct selection #####################################################################
			load(file=file.path(logfile[[1]],"quantification","results_screen_IS_pos_cal"),envir=as.environment(".GlobalEnv"));	
			at_ID<-isolate(input$Cal_IS_ID)
			use_adducts<-unique(results_screen_IS_pos_cal[[1]][
				results_screen_IS_pos_cal[[1]][,1]==at_ID
			,3,drop=FALSE])
			use_adducts<-unique(c("none",c(use_adducts)$adduct))
			updateSelectInput(session,inputId="Cal_IS_adduct",choices=use_adducts,selected = use_adducts[1])
			# update peak selection #######################################################################
			load(file=file.path(logfile[[1]],"results","pattern_pos_IS"),envir=as.environment(".GlobalEnv"),verbose=TRUE);
			use_peaks<-c()
			named<-strsplit(names(pattern_pos_IS),"_")
			for(i in 1:length(named)){
				if(named[[i]][1]==at_ID){
					use_peaks<-c(use_peaks,
						1:length(pattern_pos_IS[[i]][,1])
					)
				}
			}
			use_peaks<-c("none",as.character(unique(use_peaks)))
			updateSelectInput(session,inputId="Cal_IS_peak",choices=use_peaks,selected=use_peaks[1])
		}
		if(isolate(input$Ion_mode_Cal)=="negative"){
			# update adduct selection #####################################################################
			load(file=file.path(logfile[[1]],"quantification","results_screen_IS_neg_cal"),envir=as.environment(".GlobalEnv"));	
			at_ID<-isolate(input$Cal_IS_ID)
			use_adducts<-unique(results_screen_IS_neg_cal[[1]][
				results_screen_IS_neg_cal[[1]][,1]==at_ID
			,3,drop=FALSE])
			use_adducts<-unique(c("none",c(use_adducts)$adduct))
			updateSelectInput(session,inputId="Cal_IS_adduct",choices=use_adducts,selected = use_adducts[1])
			# update peak selection #######################################################################
			load(file=file.path(logfile[[1]],"results","pattern_neg_IS"),envir=as.environment(".GlobalEnv"),verbose=TRUE);
			use_peaks<-c()
			named<-strsplit(names(pattern_neg_IS),"_")
			for(i in 1:length(named)){
				if(named[[i]][1]==at_ID){
					use_peaks<-c(use_peaks,
						1:length(pattern_neg_IS[[i]][,1])
					)
				}
			}
			use_peaks<-c("none",as.character(unique(use_peaks)))
			updateSelectInput(session,inputId="Cal_IS_peak",choices=use_peaks,selected=use_peaks[1])
		}
	}

})
# TARGET COMPOUNDS ########################################################################################
observe({ 
	init$b
	input$Cal_target_name
	input$Cal_target_ID
	if((isolate(init$a)=="TRUE") & (isolate(input$Cal_target_ID)!="none") & (isolate(input$Cal_target_name)!="none")){
		if(isolate(input$Ion_mode_Cal)=="positive"){
			# update adduct selection #####################################################################
			load(file=file.path(logfile[[1]],"quantification","results_screen_target_pos_cal"),envir=as.environment(".GlobalEnv"));	
			at_ID<-isolate(input$Cal_target_ID)
			use_adducts<-unique(results_screen_target_pos_cal[[1]][
				results_screen_target_pos_cal[[1]][,1]==at_ID
			,3,drop=FALSE])
			use_adducts<-unique(c("none",c(use_adducts)$adduct))
			updateSelectInput(session,inputId="Cal_target_adduct",choices=use_adducts,selected = use_adducts[1])
			# update peak selection #######################################################################
			load(file=file.path(logfile[[1]],"results","pattern_pos_target"),envir=as.environment(".GlobalEnv"),verbose=TRUE);
			use_peaks<-c()
			named<-strsplit(names(pattern_pos_target),"_")
			for(i in 1:length(named)){
				if(named[[i]][1]==at_ID){
					use_peaks<-c(use_peaks,
						1:length(pattern_pos_target[[i]][,1])
					)
				}
			}
			use_peaks<-c("none",as.character(unique(use_peaks)))
			updateSelectInput(session,inputId="Cal_target_peak",choices=use_peaks,selected=use_peaks[1])
		}
		if(isolate(input$Ion_mode_Cal)=="negative"){
			# update adduct selection #####################################################################
			load(file=file.path(logfile[[1]],"quantification","results_screen_target_neg_cal"),envir=as.environment(".GlobalEnv"));	
			at_ID<-isolate(input$Cal_target_ID)
			use_adducts<-unique(results_screen_target_neg_cal[[1]][
				results_screen_target_neg_cal[[1]][,1]==at_ID
			,3,drop=FALSE])
			use_adducts<-unique(c("none",c(use_adducts)$adduct))
			updateSelectInput(session,inputId="Cal_target_adduct",choices=use_adducts,selected = use_adducts[1])
			# update peak selection #######################################################################
			load(file=file.path(logfile[[1]],"results","pattern_neg_target"),envir=as.environment(".GlobalEnv"),verbose=TRUE);
			use_peaks<-c()
			named<-strsplit(names(pattern_neg_target),"_")
			for(i in 1:length(named)){
				if(named[[i]][1]==at_ID){
					use_peaks<-c(use_peaks,
						1:length(pattern_neg_target[[i]][,1])
					)
				}
			}
			use_peaks<-c("none",as.character(unique(use_peaks)))
			updateSelectInput(session,inputId="Cal_target_peak",choices=use_peaks,selected=use_peaks[1])
		}
	}

})
###########################################################################################################

###########################################################################################################
# RETRIEVE SETS & PLOT THEM ###############################################################################
observe({ 
	init$b
	input$Cal_target_name
	input$Cal_target_ID
	input$Cal_target_adduct
	input$Cal_target_peak
	input$Cal_IS_name
	input$Cal_IS_ID
	input$Cal_IS_adduct
	input$Cal_IS_peak
	if(
		(isolate(init$a)=="TRUE") &
		(isolate(input$Cal_target_name)!="none") &
		(isolate(input$Cal_target_ID)!="none") &
		(isolate(input$Cal_target_adduct)!="none") &
		(isolate(input$Cal_target_peak)!="none") &
		(isolate(input$Cal_IS_name)!="none") &
		(isolate(input$Cal_IS_ID)!="none") &
		(isolate(input$Cal_IS_adduct)!="none") &
		(isolate(input$Cal_IS_peak)!="none") &
		(isolate(input$Cal_file_set)!="none")		
	){
		if(isolate(input$Ion_mode_Cal)=="positive"){

			load(file=file.path(logfile[[1]],"quantification","profileList_pos_cal"),envir=as.environment(".GlobalEnv"));	
			load(file=file.path(logfile[[1]],"quantification","res_IS_pos_screen_cal"),envir=as.environment(".GlobalEnv"));	
			load(file=file.path(logfile[[1]],"quantification","res_target_pos_screen_cal"),envir=as.environment(".GlobalEnv"));	
			load(file=file.path(logfile[[1]],"results","pattern_pos_target"),envir=as.environment(".GlobalEnv"),verbose=TRUE);
			load(file=file.path(logfile[[1]],"results","pattern_pos_IS"),envir=as.environment(".GlobalEnv"),verbose=TRUE);

			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			measurements<-measurements[measurements[,3]=="calibration",,drop=FALSE]
			measurements<-measurements[measurements[,4]=="positive",,drop=FALSE]
			measurements<-measurements[measurements$tag2==isolate(input$Cal_file_set),,drop=FALSE]
			#measurements<-measurements[measurements$tag2=="A",,drop=FALSE]
			IDs<-as.numeric(measurements[,1])
			calib_group<-isolate(input$Cal_file_set)
			# extract results for the IS ##############################################
			named<-strsplit(names(pattern_pos_IS),"_")
			at_IS_ID<-isolate(input$Cal_IS_ID)
			at_IS_adduct<-isolate(input$Cal_IS_adduct)
			at_IS_peak<-as.numeric(isolate(input$Cal_IS_peak))
			for(i in 1:length(named)){
				if(named[[i]][1]==at_IS_ID){
					if(named[[i]][2]==at_IS_adduct){
						break; # i is the list entry to look at
					}
				}
			}
			IS_in_file<-c()
			IS_intens<-c()
			if(length(res_IS_pos_screen_cal[[i]])>0){	
				for(j in 1:length(res_IS_pos_screen_cal[[i]])){  # use above i
					if(length(res_IS_pos_screen_cal[[i]][[j]])>0){
						if(any(IDs==j)){ # its a calibration file in the concerned group
							for(k in 1:length(res_IS_pos_screen_cal[[i]][[j]])){
								if(any(res_IS_pos_screen_cal[[i]][[j]][[k]]$Peaks[,1]==at_IS_peak)){
									IS_in_file<-c(IS_in_file,j)
									IS_intens<-c(IS_intens,
										res_IS_pos_screen_cal[[i]][[j]][[k]]$Intensity[
											res_IS_pos_screen_cal[[i]][[j]][[k]]$Peaks[,1]==at_IS_peak
										]
									)
								}
							}
						}
					}
				}
			}	
			# extract results for the target ##############################################
			named<-strsplit(names(pattern_pos_target),"_")
			at_target_ID<-isolate(input$Cal_target_ID)
			at_target_adduct<-isolate(input$Cal_target_adduct)
			at_target_peak<-as.numeric(isolate(input$Cal_target_peak))
			for(i in 1:length(named)){
				if(named[[i]][1]==at_target_ID){
					if(named[[i]][2]==at_target_adduct){
						break; # i is the list entry to look at
					}
				}
			}
			target_in_file<-c()
			target_intens<-c()
			if(length(res_target_pos_screen_cal[[i]])>0){
				for(j in 1:length(res_target_pos_screen_cal[[i]])){ # use above i
					if(length(res_target_pos_screen_cal[[i]][[j]])>0){				
						if(any(IDs==j)){ # its a calibration file in the concerned group					
							for(k in 1:length(res_target_pos_screen_cal[[i]][[j]])){						
								if(any(res_target_pos_screen_cal[[i]][[j]][[k]]$Peaks[,1]==at_target_peak)){							
									target_in_file<-c(target_in_file,j)
									target_intens<-c(target_intens,
										res_target_pos_screen_cal[[i]][[j]][[k]]$Intensity[
											res_target_pos_screen_cal[[i]][[j]][[k]]$Peaks[,1]==at_target_peak
										]
									)
								}
							}
						}
					}
				}
			}
			# match target & IS ##########################################################
			if(TRUE){#debug
				target_in_file<<-target_in_file
				target_intens<<-target_intens
				IS_in_file<<-IS_in_file
				IS_intens<<-IS_intens
				print(target_in_file)
				print(IS_in_file)
			}
			output$target_in_file <- renderText({ paste(target_in_file,collapse=",") })	
			output$IS_in_file <- renderText({ paste(IS_in_file,collapse=",") })
			
			
			
		}









	}

})
###########################################################################################################

  
 if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_obs_screening.r!")}