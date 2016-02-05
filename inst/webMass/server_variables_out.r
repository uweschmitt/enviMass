# write variables to logfile ###################################################
# observe what has changed and reset Tasks_to_do accordingly ###################
observe({
    input$savepar;
    input$saveflow;
    if(	((exists("logfile")) & (isolate(input$savepar)))||
		((exists("logfile")) & (isolate(input$saveflow)))
	){
		########################################################################
	    if(any(ls()=="logfile")){stop("illegal logfile detected #1 in server_variable_out.r!")}
		########################################################################   
		createAlert(session,anchorId = "alert_1", alertId="a1", title = NULL, 
              content="Changes in the workflow settings require a project recalculation to become effective.",
              style = "warning",append=FALSE)
		createAlert(session,anchorId = "alert_2", alertId="a2", title = NULL, 
              content="Changes in the parameter settings require a project recalculation to become effective.",
              style = "warning",append=FALSE)
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
		if(any(is.na(match(at2,at1)))){ 
			workflow_set(down="peakpick",check_node=TRUE)
		}
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
		if( any(is.na(match(at2,at1))) || any(is.na(match(at1,at2))) || any(is.na(match(at3,at4))) || any(is.na(match(at4,at3))) || any(is.na(match(at5,at6))) || any(is.na(match(at6,at5))) ){ 
			workflow_set(down="pattern",check_node=TRUE)
		}
        # recalibration ######################################################## 
        at1<-logfile$parameters[c(30,31,32,33)];
		logfile$parameters[[30]]<<-as.character(isolate(input$recal_what))
        logfile$parameters[[31]]<<-as.character(isolate(input$recal_dmz))
        logfile$parameters[[32]]<<-as.character(isolate(input$recal_ppm))
        logfile$parameters[[33]]<<-as.character(isolate(input$recal_drt))
		at2<-logfile$parameters[c(30,31,32,33)];
		if(any(is.na(match(at2,at1)))){ 
			workflow_set(down="recal",check_node=TRUE)
		}
		# replicates ###########################################################
		at1<-logfile$parameters[c(15,16,17,18)];
		logfile$parameters[[15]]<<-as.character(isolate(input$replicate_dmz))
		logfile$parameters[[16]]<<-as.character(isolate(input$replicate_ppm))
		#logfile$parameters[[17]]<<-as.character(isolate(input$replicate_recalib))
		logfile$parameters[[18]]<<-as.character(isolate(input$replicate_delRT))		
		at2<-logfile$parameters[c(15,16,17,18)];
		if(any(is.na(match(at2,at1)))){ 
			workflow_set(down="replicates",check_node=TRUE)		
		}
		# profiling ############################################################
        at1<-c(logfile$parameters[c(38,39,40,41)],logfile$parameters$upto_file);
		logfile$parameters[[38]]<<-as.character(isolate(input$prof_sets))
		logfile$parameters$upto_file<<-as.character(isolate(input$upto_file))
        logfile$parameters[[39]]<<-as.character(isolate(input$prof_dmz))
        logfile$parameters[[40]]<<-as.character(isolate(input$prof_ppm))
        logfile$parameters[[41]]<<-as.character(isolate(input$prof_drt))
		at2<-c(logfile$parameters[c(38,39,40,41)],logfile$parameters$upto_file);
		if(any(is.na(match(at2,at1)))){ 
			workflow_set(down="profiling",check_node=TRUE)		
		}
		# trend detection ######################################################
		at1<-logfile$parameters[c(29,34,35)]
		logfile$parameters[[29]]<<-as.character(isolate(input$notrend))
		logfile$parameters[[34]]<<-as.character(isolate(input$trend_lags))
		logfile$parameters[[35]]<<-as.character(isolate(input$trend_thres))
		at2<-logfile$parameters[c(29,34,35)]
		if(any(is.na(match(at2,at1)))){ 
			workflow_set(down="trendblind",check_node=TRUE)
		}		
		# blind subtraction ####################################################		
		at1<-logfile$parameters[c(36,37)]
		logfile$parameters[[36]]<<-as.character(isolate(input$blind_do))
		logfile$parameters[[37]]<<-as.character(isolate(input$blind_fold))
		at2<-logfile$parameters[c(36,37)]
		if(any(is.na(match(at2,at1)))){ 
			workflow_set(down="trendblind",check_node=TRUE)
		}		
		# IS screening #########################################################
		at1<-logfile$parameters[c(42:51)]
		logfile$parameters[[42]]<<-as.character(isolate(input$screen_IS_delRT))
		logfile$parameters[[43]]<<-as.character(isolate(input$screen_IS_dRTwithin))
		logfile$parameters[[44]]<<-as.character(isolate(input$screen_IS_dRTblank))
		logfile$parameters[[45]]<<-as.character(isolate(input$screen_IS_dmz))
		logfile$parameters[[46]]<<-as.character(isolate(input$screen_IS_ppm))
		logfile$parameters[[47]]<<-as.character(isolate(input$screen_IS_dInt))
		logfile$parameters[[48]]<<-as.character(isolate(input$screen_IS_Intcut))
		logfile$parameters[[49]]<<-as.character(isolate(input$screen_IS_w1))
		logfile$parameters[[50]]<<-as.character(isolate(input$screen_IS_w2))
		logfile$parameters[[51]]<<-as.character(isolate(input$screen_IS_w3))		
		at2<-logfile$parameters[c(42:51)]
		if(any(is.na(match(at2,at1)))){ 
			workflow_set(down="IS_screen",check_node=TRUE)
		}
		# target screening #####################################################
		at1<-logfile$parameters[c(55:64)]
		logfile$parameters[[55]]<<-as.character(isolate(input$screen_target_delRT))
		logfile$parameters[[56]]<<-as.character(isolate(input$screen_target_dRTwithin))
		logfile$parameters[[57]]<<-as.character(isolate(input$screen_target_dRTblank))
		logfile$parameters[[58]]<<-as.character(isolate(input$screen_target_dmz))
		logfile$parameters[[59]]<<-as.character(isolate(input$screen_target_ppm))
		logfile$parameters[[60]]<<-as.character(isolate(input$screen_target_dInt))
		logfile$parameters[[61]]<<-as.character(isolate(input$screen_target_Intcut))
		logfile$parameters[[62]]<<-as.character(isolate(input$screen_target_w1))
		logfile$parameters[[63]]<<-as.character(isolate(input$screen_target_w2))
		logfile$parameters[[64]]<<-as.character(isolate(input$screen_target_w3))		
		at2<-logfile$parameters[c(55:64)]
		if(any(is.na(match(at2,at1)))){ 
			workflow_set(down="target_screen",check_node=TRUE)
		}
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
		if(any(is.na(match(at2,at1)))){ 
			workflow_set(down="norm_prof",check_node=TRUE)		
		}	
		########################################################################   
	
		##########################################################################
		# workflow settings ######################################################
		# QC #####################################################################
		at1<-logfile$workflow[names(logfile$workflow)=="qc"]; 
		logfile$workflow[names(logfile$workflow)=="qc"]<<-as.character(isolate(input$qc));
		at2<-logfile$workflow[names(logfile$workflow)=="qc"];
		if(at1!=at2){
			workflow_set(down="QC",check_node=FALSE)
		}
		# recal ##################################################################	
		at1<-logfile$workflow[names(logfile$workflow)=="recal"];
		logfile$workflow[names(logfile$workflow)=="recal"]<<-as.character(isolate(input$recal));
		at2<-logfile$workflow[names(logfile$workflow)=="recal"];
		if(at1!=at2){
			workflow_set(down="recal",check_node=FALSE)
		}
		# replicates #############################################################	
		at1<-logfile$workflow[names(logfile$workflow)=="replicates"];
		logfile$workflow[names(logfile$workflow)=="replicates"]<<-as.character(isolate(input$replicates));
		at2<-logfile$workflow[names(logfile$workflow)=="replicates"];		
		if(at1!=at2){
			workflow_set(down="replicates",check_node=FALSE)
		}
		# align ##################################################################
		# norm ###################################################################
		at1<-logfile$workflow[names(logfile$workflow)=="norm"];
		logfile$workflow[names(logfile$workflow)=="norm"]<<-as.character(isolate(input$intnorm));
		at2<-logfile$workflow[names(logfile$workflow)=="norm"];
		if(at1!=at2){
			workflow_set(down="normalize",check_node=FALSE)	
		}
		# profiling ##############################################################
		at1<-logfile$workflow[names(logfile$workflow)=="profiled"]; 
		logfile$workflow[names(logfile$workflow)=="profiled"]<<-as.character(isolate(input$profiled));
		at2<-logfile$workflow[names(logfile$workflow)=="profiled"];
		if(at1!=at2){
			workflow_set(down="profiling",check_node=FALSE)		
		}
		# trendblind #############################################################
		at1<-logfile$workflow[names(logfile$workflow)=="trenddetect"];
		logfile$workflow[names(logfile$workflow)=="trenddetect"]<<-as.character(isolate(input$trenddetect));
		at2<-logfile$workflow[names(logfile$workflow)=="trenddetect"];
		if(at1!=at2){
			workflow_set(down="trendblind",check_node=FALSE)
		}
		# IS-normalization #######################################################
		at1<-logfile$workflow[names(logfile$workflow)=="profnorm"];
		logfile$workflow[names(logfile$workflow)=="profnorm"]<<-as.character(isolate(input$profnorm));
		at2<-logfile$workflow[names(logfile$workflow)=="profnorm"];
		if(at1!=at2){
			workflow_set(down="norm_prof",check_node=FALSE)
		}		
		# IS screening ###########################################################
		at1<-logfile$workflow[names(logfile$workflow)=="screen_IS"];
		logfile$workflow[names(logfile$workflow)=="screen_IS"]<<-as.character(isolate(input$screen_IS));
		at2<-logfile$workflow[names(logfile$workflow)=="screen_IS"];
		if(at1!=at2){
			workflow_set(down="IS_screen",check_node=FALSE)		
		}
		# target screening #######################################################
		at1<-logfile$workflow[names(logfile$workflow)=="screen_target"];
		logfile$workflow[names(logfile$workflow)=="screen_target"]<<-as.character(isolate(input$screen_target));
		at2<-logfile$workflow[names(logfile$workflow)=="screen_target"];
		if(at1!=at2){
			workflow_set(down="target_screen",check_node=FALSE)		
		}
		##########################################################################	

		##########################################################################
		cat("settings changed \n");
		output$dowhat<<-renderText("Project settings modified");
		##########################################################################
		if(any(ls()=="logfile")){stop("illegal logfile detected #2 in server_variable_out.r!")}
		
		
    }

})
################################################################################