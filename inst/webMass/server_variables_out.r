# write variables to logfile ###################################################
# observe what has changed and reset Tasks_to_do accordingly ###################
observe({
    input$savepar;
    input$saveflow;
    if(	((exists("logfile")) & (isolate(input$savepar)))||
		((exists("logfile")) & (isolate(input$saveflow)))
	){
		do_debug<-FALSE
		########################################################################
	    if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_variable_out.r!")}
		########################################################################   
		createAlert(session,anchorId = "alert_1", alertId="a1", title = NULL, 
              content="Changes in the workflow settings require a project recalculation to become effective.",
              style = "warning",append=FALSE)
		createAlert(session,anchorId = "alert_2", alertId="a2", title = NULL, 
              content="Changes in the parameter settings require a project recalculation to become effective.",
              style = "warning",append=FALSE)
		if(do_debug){cat("\n at_1")}	  
		######################################################################## 		

		########################################################################     
        # parameter settings ###################################################
        # PW path ##############################################################
		logfile[[4]]<<-as.character(isolate(input$PWpath));
		# peak picking #########################################################
		at1<-logfile$parameters[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)];
		logfile$parameters$peak_MSlevel<<-as.character(1); # MSlevel=1
		logfile$parameters$peak_drtgap<<-as.character(isolate(input$peak_drtgap))
		logfile$parameters$peak_dmzdens<<-as.character(isolate(input$peak_dmzdens))		
		logfile$parameters$peak_minpeak<<-as.character(isolate(input$peak_minpeak))		
		logfile$parameters$peak_drtsmall2<<-as.character(isolate(input$peak_drtsmall2))		
		logfile$parameters$peak_drtfill<<-as.character(isolate(input$peak_drtfill))		
		logfile$parameters$peak_drtdens2<<-as.character(isolate(input$peak_drtdens2))
		logfile$parameters$peak_minint_log10<<-as.character(isolate(input$peak_minint))
		logfile$parameters$peak_SN<<-as.character(isolate(input$peak_SN))
		logfile$parameters$peak_SB<<-as.character(isolate(input$peak_SB))
		logfile$parameters$peak_recurs<<-as.character(isolate(input$peak_recurs))
		logfile$parameters$peak_ended<<-as.character(isolate(input$peak_ended))		
		logfile$parameters$peak_weight<<-as.character(isolate(input$peak_weight))
		logfile$parameters$peak_maxint_log10<<-as.character(isolate(input$peak_maxint))
		at2<-logfile$parameters[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)]
		if(any(is.na(at2))){stop("\nThere was an issue reading out the new settings - maybe comma / dot separation was not fullfilled?")}		
		if(!enviMass:::comp_list(at1,at2)){ 
			enviMass:::workflow_set(down="peakpicking",check_node=TRUE)
		}
		if(do_debug){cat("\n at_2")}
        # progbar? #############################################################
        logfile$parameters$progressBar<<-as.character(isolate(input$progressbar));
        # isotope pattern calculation ##########################################
		# adducts ################################################################		
		at1<-logfile$parameters[22];
        logfile$parameters[[22]]<<-as.character(isolate(input$resolution));
		at2<-logfile$parameters[22];
		at3<-logfile$adducts_pos
		logfile$adducts_pos<<-as.character(isolate(input$adducts_pos))
		at4<-logfile$adducts_pos
		at5<-logfile$adducts_neg
		logfile$adducts_neg<<-as.character(isolate(input$adducts_neg))
		at6<-logfile$adducts_neg
		if(any(is.na(at2))){stop("\nThere was an issue reading out the new settings - maybe comma / dot separation was not fullfilled?")}		
		if( any(is.na(match(at2,at1))) || any(is.na(match(at1,at2))) || any(is.na(match(at3,at4))) || any(is.na(match(at4,at3))) || any(is.na(match(at5,at6))) || any(is.na(match(at6,at5))) ){ 
			enviMass:::workflow_set(down="pattern",check_node=TRUE)
		}
		if(do_debug){cat("\n at_3")}
        # recalibration ######################################################## 
        at1<-logfile$parameters[c(30,31,32,33)];
		logfile$parameters[[30]]<<-as.character(isolate(input$recal_what))
        logfile$parameters[[31]]<<-as.character(isolate(input$recal_dmz))
        logfile$parameters[[32]]<<-as.character(isolate(input$recal_ppm))
        logfile$parameters[[33]]<<-as.character(isolate(input$recal_drt))
		at2<-logfile$parameters[c(30,31,32,33)];
		if(any(is.na(at2))){stop("\nThere was an issue reading out the new settings - maybe comma / dot separation was not fullfilled?")}		
		if(!enviMass:::comp_list(at1,at2)){ 
			enviMass:::workflow_set(down="recal",check_node=TRUE)
		}
		if(do_debug){cat("\n at_4")}
		# replicates ###########################################################
		at1<-logfile$parameters[c(15,16,17,18,19)];
		logfile$parameters[[15]]<<-as.character(isolate(input$replicate_dmz))
		logfile$parameters[[16]]<<-as.character(isolate(input$replicate_ppm))
		#logfile$parameters[[17]]<<-as.character(isolate(input$replicate_recalib))
		logfile$parameters[[18]]<<-as.character(isolate(input$replicate_delRT))	
		logfile$parameters[[19]]<<-as.character(isolate(input$replicate_IS_dInt))	
		at2<-logfile$parameters[c(15,16,17,18,19)];
		if(any(is.na(at2))){stop("\nThere was an issue reading out the new settings - maybe comma / dot separation was not fullfilled?")}		
		if(!enviMass:::comp_list(at1,at2)){ 
			enviMass:::workflow_set(down="replicates",check_node=TRUE)		
		}
		if(do_debug){cat("\n at_5")}
		# profiling ############################################################
        at1<-c(logfile$parameters[c(38,39,40,41,90,91)],logfile$parameters$upto_file);
		logfile$parameters[[38]]<<-as.character(isolate(input$prof_sets))
		logfile$parameters$upto_file<<-as.character(isolate(input$upto_file))
        logfile$parameters[[39]]<<-as.character(isolate(input$prof_dmz))
        logfile$parameters[[40]]<<-as.character(isolate(input$prof_ppm))
        logfile$parameters[[41]]<<-as.character(isolate(input$prof_drt))
        logfile$parameters[[90]]<<-as.character(isolate(input$prof_select))	
        logfile$parameters[[91]]<<-as.character(isolate(input$replicates_prof))		
		at2<-c(logfile$parameters[c(38,39,40,41,90,91)],logfile$parameters$upto_file);
		if(any(is.na(at2))){stop("\nThere was an issue reading out the new settings - maybe comma / dot separation was not fullfilled?")}		
		if(!enviMass:::comp_list(at1,at2)){ 
			enviMass:::workflow_set(down="profiling",check_node=TRUE)		
		}
		if(do_debug){cat("\n at_6")}
		# trend detection ######################################################
		at1<-logfile$parameters[c(29,34,35,36)]
		logfile$parameters[[29]]<<-as.character(isolate(input$notrend))
		logfile$parameters[[34]]<<-as.character(isolate(input$trend_lags))
		logfile$parameters[[35]]<<-as.character(isolate(input$trend_thres))
		logfile$parameters[[36]]<<-as.character(isolate(input$trend_blind))
		at2<-logfile$parameters[c(29,34,35,36)]
		if(any(is.na(at2))){stop("\nThere was an issue reading out the new settings - maybe comma / dot separation was not fullfilled?")}		
		if(!enviMass:::comp_list(at1,at2)){ 
			enviMass:::workflow_set(down="trendblind",check_node=TRUE)
		}		
		if(do_debug){cat("\n at_7")}
		# blind subtraction ####################################################		
		at1<-logfile$parameters[c(37,82,83,84,85,86,87,88,89)]
		logfile$parameters[[37]]<<-as.character(isolate(input$blind_fold))
		logfile$parameters[[82]]<<-as.character(isolate(input$blind_dmz))
		logfile$parameters[[83]]<<-as.character(isolate(input$blind_ppm))
		logfile$parameters[[84]]<<-as.character(isolate(input$blind_drt))	
		logfile$parameters[[85]]<<-as.character(isolate(input$subtract_pos_bydate))		
		logfile$parameters[[86]]<<-as.character(isolate(input$subtract_pos_byfile))		
		logfile$parameters[[87]]<<-as.character(isolate(input$subtract_neg_bydate))			
		logfile$parameters[[88]]<<-as.character(isolate(input$subtract_neg_byfile))	
		logfile$parameters[[89]]<<-as.character(isolate(input$blind_omit))
		at2<-logfile$parameters[c(37,82,83,84,85,86,87,88,89)]
		if(any(is.na(at2))){stop("\nThere was an issue reading out the new settings - maybe comma / dot separation was not fullfilled?")}		
		if(!enviMass:::comp_list(at1,at2)){ # both steps take partly the same parameters! 
			enviMass:::workflow_set(down="blinds",check_node=TRUE)
		}		
		if(!enviMass:::comp_list(at1,at2)){ 
			enviMass:::workflow_set(down="trendblind",check_node=TRUE)
		}		
		if(do_debug){cat("\n at_8a")}
		# subtraction files
		at1<-logfile$Positive_subtraction_files
		logfile$Positive_subtraction_files<<-c(isolate(input$files_pos_select_subtract),"FALSE")
		at2<-logfile$Positive_subtraction_files
		if(any(is.na(at2))){stop("\nThere was an issue reading out the new settings - maybe comma / dot separation was not fullfilled?")}		
		if(!enviMass:::comp_list(at1,at2,as_pairs=FALSE)){ # both steps take partly the same parameters! 
			enviMass:::workflow_set(down="blinds",check_node=TRUE)
		}	
		at1<-logfile$Negative_subtraction_files
		logfile$Negative_subtraction_files<<-c(isolate(input$files_neg_select_subtract),"FALSE")
		at2<-logfile$Negative_subtraction_files
		if(any(is.na(at2))){stop("\nThere was an issue reading out the new settings - maybe comma / dot separation was not fullfilled?")}		
		if(!enviMass:::comp_list(at1,at2,as_pairs=FALSE)){ # both steps take partly the same parameters! 
			enviMass:::workflow_set(down="blinds",check_node=TRUE)
		}		
		if(do_debug){cat("\n at_8b")}
		########################################################################  		
		# IS screening #########################################################
		at1<-logfile$parameters[c(42,43,45:51)]
		logfile$parameters[[42]]<<-as.character(isolate(input$screen_IS_delRT))
		logfile$parameters[[43]]<<-as.character(isolate(input$screen_IS_dRTwithin))
		logfile$parameters[[45]]<<-as.character(isolate(input$screen_IS_dmz))
		logfile$parameters[[46]]<<-as.character(isolate(input$screen_IS_ppm))
		logfile$parameters[[47]]<<-as.character(isolate(input$screen_IS_dInt))
		logfile$parameters[[48]]<<-as.character(isolate(input$screen_IS_Intcut))
		use_value<-as.character(isolate(input$screen_IS_w1))
		if(is.na(use_value)){stop("\n\nInvalid parameter - dot instead of comma-separated numeric input used?")}
		logfile$parameters[[49]]<<-as.character(isolate(input$screen_IS_w1))		
		logfile$parameters[[50]]<<-as.character(isolate(input$screen_IS_cutit))		
		logfile$parameters[[51]]<<-as.character(isolate(input$screen_IS_maxonly))			
		at2<-logfile$parameters[c(42,43,45:51)]
		if(any(is.na(at2))){stop("\nThere was an issue reading out the new settings - maybe comma / dot separation was not fullfilled?")}		
		if(!enviMass:::comp_list(at1,at2)){ 
			enviMass:::workflow_set(down="IS_screen",check_node=TRUE) 
			enviMass:::workflow_set(down="pattern",check_node=TRUE) # because delRT is set during pattern calculation
		}	
		if(do_debug){cat("\n at_9")}
		# target screening #####################################################
		at1<-logfile$parameters[c(55,56,58:64)]
		logfile$parameters[[55]]<<-as.character(isolate(input$screen_target_delRT))
		logfile$parameters[[56]]<<-as.character(isolate(input$screen_target_dRTwithin))
		logfile$parameters[[58]]<<-as.character(isolate(input$screen_target_dmz))
		logfile$parameters[[59]]<<-as.character(isolate(input$screen_target_ppm))
		logfile$parameters[[60]]<<-as.character(isolate(input$screen_target_dInt))
		logfile$parameters[[61]]<<-as.character(isolate(input$screen_target_Intcut))
		use_value<-as.character(isolate(input$screen_target_w1))
		if(is.na(use_value)){stop("\n\nInvalid parameter - dot instead of comma-separated numeric input used?")}		
		logfile$parameters[[62]]<<-as.character(isolate(input$screen_target_w1))	
		logfile$parameters[[63]]<<-as.character(isolate(input$screen_target_cutit))		
		logfile$parameters[[64]]<<-as.character(isolate(input$screen_target_maxonly))			
		at2<-logfile$parameters[c(55,56,58:64)]
		if(any(is.na(at2))){stop("\nThere was an issue reading out the new settings - maybe comma / dot separation was not fullfilled?")}		
		if(!enviMass:::comp_list(at1,at2)){ 
			enviMass:::workflow_set(down="target_screen",check_node=TRUE)
			enviMass:::workflow_set(down="pattern",check_node=TRUE) # because delRT is set during pattern calculation
		}
		if(do_debug){cat("\n at_10")}
		# IS-Normalization #####################################################
		at1<-logfile$parameters[c(70:78)]
		logfile$parameters[[70]]<<-as.character(isolate(input$profnorm_cover_files))
		logfile$parameters[[71]]<<-as.character(isolate(input$profnorm_cover_isccount))
		logfile$parameters[[72]]<<-as.character(isolate(input$profnorm_use_blank))
		logfile$parameters[[73]]<<-as.character(isolate(input$profnorm_use_blank_sample))
		logfile$parameters[[74]]<<-as.character(isolate(input$profnorm_use_blank_samplecount))
		logfile$parameters[[75]]<<-as.character(isolate(input$profnorm_use_nonblank))
		logfile$parameters[[76]]<<-as.character(isolate(input$profnorm_use_nonblank_sample))
		logfile$parameters[[77]]<<-as.character(isolate(input$profnorm_use_nonblank_samplecount))	
		logfile$parameters[[78]]<<-as.character(isolate(input$profnorm_threshold))
		at2<-logfile$parameters[c(70:78)]
		if(any(is.na(at2))){stop("\nThere was an issue reading out the new settings - maybe comma / dot separation was not fullfilled?")}		
		if(!enviMass:::comp_list(at1,at2)){ 
			enviMass:::workflow_set(down="IS_normaliz",check_node=TRUE)		
		}	
		if(do_debug){cat("\n at_11")}
		########################################################################   
		
		##########################################################################
		# workflow settings ######################################################
		# qc #####################################################################
		at1<-logfile$workflow[names(logfile$workflow)=="qc"]; 
		logfile$workflow[names(logfile$workflow)=="qc"]<<-as.character(isolate(input$qc));
		at2<-logfile$workflow[names(logfile$workflow)=="qc"];
		if(any(is.na(at2))){stop("\nThere was an issue reading out the new settings - maybe comma / dot separation was not fullfilled?")}		
		if(at1!=at2){
			enviMass:::workflow_set(down="qc",check_node=FALSE)
		}
		if(do_debug){cat("\n at_12")}
		# recal ##################################################################	
		at1<-logfile$workflow[names(logfile$workflow)=="recal"];
		logfile$workflow[names(logfile$workflow)=="recal"]<<-as.character(isolate(input$recal));
		at2<-logfile$workflow[names(logfile$workflow)=="recal"];
		if(any(is.na(at2))){stop("\nThere was an issue reading out the new settings - maybe comma / dot separation was not fullfilled?")}		
		if(at1!=at2){
			enviMass:::workflow_set(down="recal",check_node=FALSE)
		}
		if(do_debug){cat("\n at_13")}
		# blinds #################################################################	
		at1<-logfile$workflow[names(logfile$workflow)=="blinds"];
		logfile$workflow[names(logfile$workflow)=="blinds"]<<-as.character(isolate(input$blind_filter));
		at2<-logfile$workflow[names(logfile$workflow)=="blinds"];
		if(any(is.na(at2))){stop("\nThere was an issue reading out the new settings - maybe comma / dot separation was not fullfilled?")}		
		if(at1!=at2){
			enviMass:::workflow_set(down="blinds",check_node=FALSE)
		}		
		if(do_debug){cat("\n at_14")}
		# replicates #############################################################	
		at1<-logfile$workflow[names(logfile$workflow)=="replicates"];
		logfile$workflow[names(logfile$workflow)=="replicates"]<<-as.character(isolate(input$replicates));
		at2<-logfile$workflow[names(logfile$workflow)=="replicates"];		
		if(any(is.na(at2))){stop("\nThere was an issue reading out the new settings - maybe comma / dot separation was not fullfilled?")}		
		if(at1!=at2){
			enviMass:::workflow_set(down="replicates",check_node=FALSE)
		}
		if(do_debug){cat("\n at_15")}
		# align ##################################################################
		# norm ###################################################################
		at1<-logfile$workflow[names(logfile$workflow)=="norm"];
		logfile$workflow[names(logfile$workflow)=="norm"]<<-as.character(isolate(input$intnorm));
		at2<-logfile$workflow[names(logfile$workflow)=="norm"];
		if(any(is.na(at2))){stop("\nThere was an issue reading out the new settings - maybe comma / dot separation was not fullfilled?")}		
		if(at1!=at2){
			enviMass:::workflow_set(down="norm",check_node=FALSE)	
		}
		if(do_debug){cat("\n at_16")}
		# profiling ##############################################################
		at1<-logfile$workflow[names(logfile$workflow)=="profiling"]; 
		logfile$workflow[names(logfile$workflow)=="profiling"]<<-as.character(isolate(input$profiled));
		at2<-logfile$workflow[names(logfile$workflow)=="profiling"];
		if(any(is.na(at2))){stop("\nThere was an issue reading out the new settings - maybe comma / dot separation was not fullfilled?")}		
		if(at1!=at2){
			enviMass:::workflow_set(down="profiling",check_node=FALSE)		
		}
		if(do_debug){cat("\n at_17")}
		# trendblind #############################################################
		at1<-logfile$workflow[names(logfile$workflow)=="trendblind"];
		logfile$workflow[names(logfile$workflow)=="trendblind"]<<-as.character(isolate(input$trenddetect));
		at2<-logfile$workflow[names(logfile$workflow)=="trendblind"];
		if(any(is.na(at2))){stop("\nThere was an issue reading out the new settings - maybe comma / dot separation was not fullfilled?")}		
		if(at1!=at2){
			enviMass:::workflow_set(down="trendblind",check_node=FALSE)
		}
		if(do_debug){cat("\n at_18")}
		# IS-normalization #######################################################
		at1<-logfile$workflow[names(logfile$workflow)=="IS_normaliz"];
		logfile$workflow[names(logfile$workflow)=="IS_normaliz"]<<-as.character(isolate(input$profnorm));
		at2<-logfile$workflow[names(logfile$workflow)=="IS_normaliz"];
		if(any(is.na(at2))){stop("\nThere was an issue reading out the new settings - maybe comma / dot separation was not fullfilled?")}		
		if(at1!=at2){
			enviMass:::workflow_set(down="IS_normaliz",check_node=FALSE)
		}		
		if(do_debug){cat("\n at_19")}
		# IS screening ###########################################################
		at1<-logfile$workflow[names(logfile$workflow)=="IS_screen"];
		logfile$workflow[names(logfile$workflow)=="IS_screen"]<<-as.character(isolate(input$screen_IS));
		at2<-logfile$workflow[names(logfile$workflow)=="IS_screen"];
		if(any(is.na(at2))){stop("\nThere was an issue reading out the new settings - maybe comma / dot separation was not fullfilled?")}		
		if(at1!=at2){
			enviMass:::workflow_set(down="IS_screen",check_node=FALSE)		
		}
		if(do_debug){cat("\n at_20")}
		# target screening #######################################################
		at1<-logfile$workflow[names(logfile$workflow)=="target_screen"];
		logfile$workflow[names(logfile$workflow)=="target_screen"]<<-as.character(isolate(input$screen_target));
		at2<-logfile$workflow[names(logfile$workflow)=="target_screen"];
		if(any(is.na(at2))){stop("\nThere was an issue reading out the new settings - maybe comma / dot separation was not fullfilled?")}		
		if(at1!=at2){
			enviMass:::workflow_set(down="target_screen",check_node=FALSE)		
		}
		if(do_debug){cat("\n at_21")}
		# LOD ###################################################################
		at1<-logfile$workflow[names(logfile$workflow)=="LOD"];
		logfile$workflow[names(logfile$workflow)=="LOD"]<<-as.character(isolate(input$LOD_interpol));
		at2<-logfile$workflow[names(logfile$workflow)=="LOD"];
		if(at1!=at2){
			enviMass:::workflow_set(down="LOD",check_node=FALSE)		
		}			
		if(do_debug){cat("\n at_22")}
		# quantification #######################################################
		at1<-logfile$workflow[names(logfile$workflow)=="quantification"];
		logfile$workflow[names(logfile$workflow)=="quantification"]<<-as.character(isolate(input$quantif));
		at2<-logfile$workflow[names(logfile$workflow)=="quantification"];
		if(any(is.na(at2))){stop("\nThere was an issue reading out the new settings - maybe comma / dot separation was not fullfilled?")}		
		if(at1!=at2){
			enviMass:::workflow_set(down="quantification",check_node=FALSE)		
		}		
		if(do_debug){cat("\n at_23")}
		# IS subtraction #########################################################
		at1<-logfile$workflow[names(logfile$workflow)=="IS_subtr"];
		logfile$workflow[names(logfile$workflow)=="IS_subtr"]<<-as.character(isolate(input$subtr_IS));
		at2<-logfile$workflow[names(logfile$workflow)=="IS_subtr"];
		if(any(is.na(at2))){stop("\nThere was an issue reading out the new settings - maybe comma / dot separation was not fullfilled?")}		
		if(at1!=at2){
			enviMass:::workflow_set(down="IS_subtr",check_node=FALSE)		
		}				
		if(do_debug){cat("\n at_24")}
		# target subtraction #####################################################			
		at1<-logfile$workflow[names(logfile$workflow)=="target_subtr"];
		logfile$workflow[names(logfile$workflow)=="target_subtr"]<<-as.character(isolate(input$subtr_target));
		at2<-logfile$workflow[names(logfile$workflow)=="target_subtr"];
		if(any(is.na(at2))){stop("\nThere was an issue reading out the new settings - maybe comma / dot separation was not fullfilled?")}		
		if(at1!=at2){
			enviMass:::workflow_set(down="target_subtr",check_node=FALSE)		
		}	
		if(do_debug){cat("\n at_25")}
		# blind subtraction ######################################################			
		at1<-logfile$workflow[names(logfile$workflow)=="blind_subtr"];
		logfile$workflow[names(logfile$workflow)=="blind_subtr"]<<-as.character(isolate(input$subtr_blind));
		at2<-logfile$workflow[names(logfile$workflow)=="blind_subtr"];
		if(any(is.na(at2))){stop("\nThere was an issue reading out the new settings - maybe comma / dot separation was not fullfilled?")}		
		if(at1!=at2){
			enviMass:::workflow_set(down="blind_subtr",check_node=FALSE)		
		}	
		if(do_debug){cat("\n at_26")}
		##########################################################################	

		##########################################################################
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		cat("settings changed \n");
		output$dowhat<<-renderText("Project settings modified");
		##########################################################################
		if(any(ls()=="logfile")){stop("\n illegal logfile detected #2 in server_variable_out.r!")}
		if(do_debug){cat("\n at_27")}
		
    }

})
################################################################################