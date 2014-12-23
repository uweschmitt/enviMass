# write variables to logfile ###################################################
# observe what has changed and reset Tasks_to_do accordingly ###################
observe({
    input$savepar;
    input$saveflow;
    if(	((exists("logfile")) & (isolate(input$savepar)))||
		((exists("logfile")) & (isolate(input$saveflow)))
	){
		createAlert(session,inputId = "alert_1", alertId="a1", title = NULL, 
              message="Changes in the workflow settings require a project recalculation to become effective.",
              type = "warning",append=FALSE)
		createAlert(session,inputId = "alert_2", alertId="a2", title = NULL, 
              message="Changes in the parameter settings require a project recalculation to become effective.",
              type = "warning",append=FALSE)
		##########################################################################      
        # parameter settings #####################################################
        # PW path ##############################################################
		logfile[[4]]<-as.character(isolate(input$PWpath));
		# peak picking #########################################################
		at1<-logfile$parameters[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)];
		logfile$parameters$peak_MSlevel<-as.character(1); # MSlevel=1
		logfile$parameters$peak_drtgap<-as.character(isolate(input$peak_drtgap))
		logfile$parameters$peak_dmzdens<-as.character(isolate(input$peak_dmzdens))		
		logfile$parameters$peak_minpeak<-as.character(isolate(input$peak_minpeak))		
		logfile$parameters$peak_drtsmall2<-as.character(isolate(input$peak_drtsmall2))		
		logfile$parameters$peak_drtfill<-as.character(isolate(input$peak_drtfill))		
		logfile$parameters$peak_drtdens2<-as.character(isolate(input$peak_drtdens2))
		logfile$parameters$peak_minint_log10<-as.character(isolate(input$peak_minint))
		logfile$parameters$peak_SN<-as.character(isolate(input$peak_SN))
		logfile$parameters$peak_SB<-as.character(isolate(input$peak_SB))
		logfile$parameters$peak_recurs<-as.character(isolate(input$peak_recurs))
		logfile$parameters$peak_ended<-as.character(isolate(input$peak_ended))		
		logfile$parameters$peak_weight<-as.character(isolate(input$peak_weight))
		logfile$parameters$peak_maxint_log10<-as.character(isolate(input$peak_maxint))
		at2<-logfile$parameters[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)]
		if(any(is.na(match(at2,at1)))){ 
			# must rerun: everything
			# must rerun: peak picking
			logfile$Tasks_to_redo[1]<-TRUE;
			measurements1<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			measurements1[,10]<-FALSE;
			write.csv(measurements1,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
			rm(measurements1);			
			if(logfile$workflow[1]=="yes"){# must rerun: qc
				logfile$Tasks_to_redo[2]<-TRUE;			
			}
			if(logfile$workflow[2]=="yes"){# must rerun: recal
				logfile$Tasks_to_redo[3]<-TRUE;			
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				measurements[,12]<-FALSE;
				write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);				
			}
			if(logfile$workflow[3]=="yes"){# must rerun: align
				logfile$Tasks_to_redo[5]<-TRUE;			
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				measurements[,13]<-FALSE;
				write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
			}
			if(logfile$workflow[4]=="yes"){# must rerun: norm
				logfile$Tasks_to_redo[4]<-TRUE;			
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				measurements[,14]<-FALSE;
				write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
			}
			# skip is_pattern
			# skip target_pattern
			if(logfile$workflow[7]=="TRUE"){# must rerun: component isotopologues
				logfile$Tasks_to_redo[13]<-TRUE;	
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				measurements[,16]<-FALSE;
				write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);				
			}
			if(logfile$workflow[8]=="TRUE"){# must rerun: component adducts
				logfile$Tasks_to_redo[13]<-TRUE;
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				measurements[,16]<-FALSE;
				write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);				
			}
			if(logfile$workflow[9]=="yes"){# must rerun: profiled
				logfile$Tasks_to_redo[6]<-TRUE;				
			}
			if(logfile$workflow[10]=="yes"){# must rerun: trendblind
				logfile$Tasks_to_redo[7]<-TRUE;			
			}
			if(logfile$workflow[11]=="TRUE"){# must rerun: screen_IS_sam  
				logfile$Tasks_to_redo[10]<-TRUE;			
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				measurements[,17]<-FALSE;
				write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);												
			}
			if(logfile$workflow[12]=="TRUE"){# must rerun: screen_target_sam
				logfile$Tasks_to_redo[9]<-TRUE;			
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				measurements[,18]<-FALSE;
				write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);												
			}		
			if(logfile$workflow[13]=="TRUE"){# must rerun: screen_IS_comp
				logfile$Tasks_to_redo[12]<-TRUE;			
			}			
			if(logfile$workflow[14]=="TRUE"){# must rerun: screen_target_comp
				logfile$Tasks_to_redo[11]<-TRUE;		
			}
			if(logfile$workflow[15]=="yes"){# must rerun: profnorm - IS-based normalization
				logfile$Tasks_to_redo[15]<-TRUE;			
			}			
			if(logfile$workflow[16]=="yes"){# must rerun: homol
				logfile$Tasks_to_redo[14]<-TRUE;			
			}			
			if(logfile$workflow[17]=="yes"){# must rerun: massdef
				logfile$Tasks_to_redo[16]<-TRUE;			
			}			
		}
        # progbar? #############################################################
        logfile$parameters$progressBar<-as.character(isolate(input$progressbar));
        # isotope pattern calculation ##########################################
		# adducts ################################################################		
		at1<-logfile$parameters[22];
        logfile$parameters[[22]]<-as.character(isolate(input$resolution));
		at2<-logfile$parameters[22];
		at3<-logfile$adducts_pos
		logfile$adducts_pos<-as.character(isolate(input$adducts_pos))
		at4<-logfile$adducts_pos
		at5<-logfile$adducts_neg
		logfile$adducts_neg<-as.character(isolate(input$adducts_neg))
		at6<-logfile$adducts_neg
		if( any(is.na(match(at2,at1))) || any(is.na(match(at1,at2))) || any(is.na(match(at3,at4))) || any(is.na(match(at4,at3))) || any(is.na(match(at5,at6))) || any(is.na(match(at6,at5))) ){ 
			if(logfile$workflow[2]=="yes"){# must rerun: recal
				logfile$Tasks_to_redo[3]<-TRUE;			
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				measurements[,12]<-FALSE;
				write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);				
			}
			# must not rerun: align
			# must not rerun: norm
			# must rerun: is_pattern / target_pattern
			logfile$Tasks_to_redo[8]<-TRUE;		
			# must not rerun: component isotopologues
			# must not rerun: component adducts
			if(logfile$workflow[9]=="yes"){# must rerun: profiled
				logfile$Tasks_to_redo[6]<-TRUE;				
			}
			if(logfile$workflow[10]=="yes"){# must rerun: trendblind
				logfile$Tasks_to_redo[7]<-TRUE;			
			}
			if(logfile$workflow[11]=="TRUE"){# must rerun: screen_IS_sam  
				logfile$Tasks_to_redo[10]<-TRUE;			
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				measurements[,17]<-FALSE;
				write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);												
			}
			if(logfile$workflow[12]=="TRUE"){# must rerun: screen_target_sam
				logfile$Tasks_to_redo[9]<-TRUE;			
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				measurements[,18]<-FALSE;
				write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);												
			}		
			if(logfile$workflow[13]=="TRUE"){# must rerun: screen_IS_comp
				logfile$Tasks_to_redo[12]<-TRUE;			
			}			
			if(logfile$workflow[14]=="TRUE"){# must rerun: screen_target_comp
				logfile$Tasks_to_redo[11]<-TRUE;			
			}
			if(logfile$workflow[15]=="yes"){# must rerun: profnorm - IS-based normalization
				logfile$Tasks_to_redo[15]<-TRUE;			
			}			
			if(logfile$workflow[16]=="yes"){# must rerun: homol
				logfile$Tasks_to_redo[14]<-TRUE;			
			}			
			if(logfile$workflow[17]=="yes"){# must rerun: massdef
				logfile$Tasks_to_redo[16]<-TRUE;		
			}			
		}
        # recalibration ######################################################## 
        at1<-logfile$parameters[c(30,31,32,33)];
		logfile$parameters[[30]]<-as.character(isolate(input$recal_what))
        logfile$parameters[[31]]<-as.character(isolate(input$recal_dmz))
        logfile$parameters[[32]]<-as.character(isolate(input$recal_ppm))
        logfile$parameters[[33]]<-as.character(isolate(input$recal_drt))
		at2<-logfile$parameters[c(30,31,32,33)];
		if(any(is.na(match(at2,at1)))){ 
			# must rerun: peak picking
			# must not rerun: qc
			if(logfile$workflow[2]=="yes"){# must rerun: recal
				logfile$Tasks_to_redo[3]<-TRUE;			
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				measurements[,12]<-FALSE;
				write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);				
			}
			if(logfile$workflow[3]=="yes"){# must rerun: align
				logfile$Tasks_to_redo[5]<-TRUE;			
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				measurements[,13]<-FALSE;
				write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
			}
			# must not rerun: norm
			# skip is_pattern
			# skip target_pattern
			if(logfile$workflow[7]=="TRUE"){# must rerun: component isotopologues - but NOT on file level (measurements[,16])!
				logfile$Tasks_to_redo[13]<-TRUE;					
			}
			if(logfile$workflow[8]=="TRUE"){# must rerun: component adducts - but NOT on file level (measurements[,16])!
				logfile$Tasks_to_redo[13]<-TRUE;				
			}			
			if(logfile$workflow[9]=="yes"){# must rerun: profiled
				logfile$Tasks_to_redo[6]<-TRUE;			
			}
			if(logfile$workflow[10]=="yes"){# must rerun: trendblind
				logfile$Tasks_to_redo[7]<-TRUE;		
			}
			if(logfile$workflow[11]=="TRUE"){# must rerun: screen_IS_sam  
				logfile$Tasks_to_redo[10]<-TRUE;			
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				measurements[,17]<-FALSE;
				write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);												
			}
			if(logfile$workflow[12]=="TRUE"){# must rerun: screen_target_sam
				logfile$Tasks_to_redo[9]<-TRUE;			
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				measurements[,18]<-FALSE;
				write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);												
			}		
			if(logfile$workflow[13]=="TRUE"){# must rerun: screen_IS_comp
				logfile$Tasks_to_redo[12]<-TRUE;			
			}			
			if(logfile$workflow[14]=="TRUE"){# must rerun: screen_target_comp
				logfile$Tasks_to_redo[11]<-TRUE;			
			}
			if(logfile$workflow[15]=="yes"){# must rerun: profnorm - IS-based normalization
				logfile$Tasks_to_redo[15]<-TRUE;			
			}			
			if(logfile$workflow[16]=="yes"){# must rerun: homol
				logfile$Tasks_to_redo[14]<-TRUE;		
			}			
			if(logfile$workflow[17]=="yes"){# must rerun: massdef
				logfile$Tasks_to_redo[16]<-TRUE;		
			}			
		}
		# profiling ############################################################
        at1<-c(logfile$parameters[c(38,39,40,41)],logfile$parameters$upto_file);
		logfile$parameters[[38]]<-as.character(isolate(input$prof_sets))
		logfile$parameters$upto_file<-as.character(isolate(input$upto_file))
        logfile$parameters[[39]]<-as.character(isolate(input$prof_dmz))
        logfile$parameters[[40]]<-as.character(isolate(input$prof_ppm))
        logfile$parameters[[41]]<-as.character(isolate(input$prof_drt))
		at2<-c(logfile$parameters[c(38,39,40,41)],logfile$parameters$upto_file);
		if(any(is.na(match(at2,at1)))){ 
			# must not rerun: qc
			# must not rerun: recal
			# must not rerun: align
			# must not rerun: norm
			# skip is_pattern
			# skip target_pattern
			if(logfile$workflow[7]=="TRUE"){# must rerun: component isotopologues - but NOT on file level (measurements[,16])!
				logfile$Tasks_to_redo[13]<-TRUE;				
			}
			if(logfile$workflow[8]=="TRUE"){# must rerun: component adducts - but NOT on file level (measurements[,16])!
				logfile$Tasks_to_redo[13]<-TRUE;				
			}			
			if(logfile$workflow[9]=="yes"){# must rerun: profiled
				logfile$Tasks_to_redo[6]<-TRUE;			
			}
			if(logfile$workflow[10]=="yes"){# must rerun: trendblind
				logfile$Tasks_to_redo[7]<-TRUE;			
			}
			# must not rerun: screen_IS_sam  
			# must not rerun: screen_target_sam
			if(logfile$workflow[13]=="TRUE"){# must rerun: screen_IS_comp
				logfile$Tasks_to_redo[12]<-TRUE;			
			}			
			if(logfile$workflow[14]=="TRUE"){# must rerun: screen_target_comp
				logfile$Tasks_to_redo[11]<-TRUE;		
			}
			if(logfile$workflow[15]=="yes"){# must rerun: profnorm - IS-based normalization
				logfile$Tasks_to_redo[15]<-TRUE;		
			}			
			if(logfile$workflow[16]=="yes"){# must rerun: homol
				logfile$Tasks_to_redo[14]<-TRUE;			
			}			
			if(logfile$workflow[17]=="yes"){# must rerun: massdef
				logfile$Tasks_to_redo[16]<-TRUE;		
			}			
		}
		# trend detection ######################################################
		at1<-logfile$parameters[c(34,35)]
		logfile$parameters[[34]]<-as.character(isolate(input$trend_lags))
		logfile$parameters[[35]]<-as.character(isolate(input$trend_thres))
		at2<-logfile$parameters[c(34,35)]
		if(any(is.na(match(at2,at1)))){ 
			# must not rerun: qc
			# must not rerun: recal
			# must not rerun: align
			# must not rerun: norm
			# skip is_pattern
			# skip target_pattern
			if(logfile$workflow[7]=="TRUE"){# must rerun: component isotopologues - but NOT on file level (measurements[,16])!
				logfile$Tasks_to_redo[13]<-TRUE;				
			}
			if(logfile$workflow[8]=="TRUE"){# must rerun: component adducts - but NOT on file level (measurements[,16])!
				logfile$Tasks_to_redo[13]<-TRUE;			
			}			
			# must not rerun: profiled
			if(logfile$workflow[10]=="yes"){# must rerun: trendblind
				logfile$Tasks_to_redo[7]<-TRUE;			
			}
			# must not rerun: screen_IS_sam  
			# must not rerun: screen_target_sam
			if(logfile$workflow[13]=="TRUE"){# must rerun: screen_IS_comp
				logfile$Tasks_to_redo[12]<-TRUE;			
			}			
			if(logfile$workflow[14]=="TRUE"){# must rerun: screen_target_comp
				logfile$Tasks_to_redo[11]<-TRUE;	
			}
			if(logfile$workflow[15]=="yes"){# must rerun: profnorm - IS-based normalization
				logfile$Tasks_to_redo[15]<-TRUE;		
			}			
			if(logfile$workflow[16]=="yes"){# must rerun: homol
				logfile$Tasks_to_redo[14]<-TRUE;		
			}			
			if(logfile$workflow[17]=="yes"){# must rerun: massdef
				logfile$Tasks_to_redo[16]<-TRUE;	
			}			
		}		
		# blind subtraction ####################################################		
		at1<-logfile$parameters[c(36,37)]
		logfile$parameters[[36]]<-as.character(isolate(input$blind_do))
		logfile$parameters[[37]]<-as.character(isolate(input$blind_fold))
		at2<-logfile$parameters[c(36,37)]
		if(any(is.na(match(at2,at1)))){ 
			# must not rerun: qc
			# must not rerun: recal
			# must not rerun: align
			# must not rerun: norm
			# skip is_pattern
			# skip target_pattern
			if(logfile$workflow[7]=="TRUE"){# must rerun: component isotopologues - but NOT on file level (measurements[,16])!
				logfile$Tasks_to_redo[13]<-TRUE;					
			}
			if(logfile$workflow[8]=="TRUE"){# must rerun: component adducts - but NOT on file level (measurements[,16])!
				logfile$Tasks_to_redo[13]<-TRUE;			
			}			
			# must not rerun: profiled
			if(logfile$workflow[10]=="yes"){# must rerun: trendblind
				logfile$Tasks_to_redo[7]<-TRUE;		
			}
			# must not rerun: screen_IS_sam  
			# must not rerun: screen_target_sam
			if(logfile$workflow[13]=="TRUE"){# must rerun: screen_IS_comp
				logfile$Tasks_to_redo[12]<-TRUE;			
			}			
			if(logfile$workflow[14]=="TRUE"){# must rerun: screen_target_comp
				logfile$Tasks_to_redo[11]<-TRUE;			
			}
			if(logfile$workflow[15]=="yes"){# must rerun: profnorm - IS-based normalization
				logfile$Tasks_to_redo[15]<-TRUE;			
			}			
			if(logfile$workflow[16]=="yes"){# must rerun: homol
				logfile$Tasks_to_redo[14]<-TRUE;	
			}			
			if(logfile$workflow[17]=="yes"){# must rerun: massdef
				logfile$Tasks_to_redo[16]<-TRUE;
			}			
		}		
		# IS screening #########################################################
		at1<-logfile$parameters[c(42:51)]
		logfile$parameters[[42]]<-as.character(isolate(input$screen_IS_delRT))
		logfile$parameters[[43]]<-as.character(isolate(input$screen_IS_dRTwithin))
		logfile$parameters[[44]]<-as.character(isolate(input$screen_IS_dRTblank))
		logfile$parameters[[45]]<-as.character(isolate(input$screen_IS_dmz))
		logfile$parameters[[46]]<-as.character(isolate(input$screen_IS_ppm))
		logfile$parameters[[47]]<-as.character(isolate(input$screen_IS_dInt))
		logfile$parameters[[48]]<-as.character(isolate(input$screen_IS_Intcut))
		logfile$parameters[[49]]<-as.character(isolate(input$screen_IS_w1))
		logfile$parameters[[50]]<-as.character(isolate(input$screen_IS_w2))
		logfile$parameters[[51]]<-as.character(isolate(input$screen_IS_w3))		
		at2<-logfile$parameters[c(42:51)]
		if(any(is.na(match(at2,at1)))){ 
			# must not rerun: qc
			# must not rerun: recal
			# must not rerun: align
			# must not rerun: norm
			# skip is_pattern
			# skip target_pattern
			# must not rerun: component isotopologues
			# must not rerun: component adducts
			# must not rerun: profiled
			# must not rerun: trendblind
			if(logfile$workflow[11]=="TRUE"){# must rerun: screen_IS_sam  
				logfile$Tasks_to_redo[10]<-TRUE;			
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				measurements[,17]<-FALSE;
				write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);												
			}
			# must not rerun: screen_target_sam		
			if(logfile$workflow[13]=="TRUE"){# must rerun: screen_IS_comp
				logfile$Tasks_to_redo[12]<-TRUE;			
			}			
			# must not rerun: screen_target_comp
			# must not rerun: profnorm - IS-based normalization
			# must not rerun: homol
			# must not rerun: massdef
		}
		# target screening #####################################################
		at1<-logfile$parameters[c(55:64)]
		logfile$parameters[[55]]<-as.character(isolate(input$screen_target_delRT))
		logfile$parameters[[56]]<-as.character(isolate(input$screen_target_dRTwithin))
		logfile$parameters[[57]]<-as.character(isolate(input$screen_target_dRTblank))
		logfile$parameters[[58]]<-as.character(isolate(input$screen_target_dmz))
		logfile$parameters[[59]]<-as.character(isolate(input$screen_target_ppm))
		logfile$parameters[[60]]<-as.character(isolate(input$screen_target_dInt))
		logfile$parameters[[61]]<-as.character(isolate(input$screen_target_Intcut))
		logfile$parameters[[62]]<-as.character(isolate(input$screen_target_w1))
		logfile$parameters[[63]]<-as.character(isolate(input$screen_target_w2))
		logfile$parameters[[64]]<-as.character(isolate(input$screen_target_w3))		
		at2<-logfile$parameters[c(55:64)]
		if(any(is.na(match(at2,at1)))){ 
			# must not rerun: qc
			# must not rerun: recal
			# must not rerun: align
			# must not rerun: norm
			# skip is_pattern
			# skip target_pattern
			# must not rerun: component isotopologues
			# must not rerun: component adducts
			# must not rerun: profiled
			# must not rerun: trendblind
			# must rerun: screen_IS_sam  
			# must not rerun: screen_target_sam	
			if(logfile$workflow[12]=="TRUE"){# must rerun: screen_target_sam
				logfile$Tasks_to_redo[9]<-TRUE;			
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				measurements[,18]<-FALSE;
				write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);												
			}							
			# must not rerun: screen_IS_comp
			# must  not rerun: screen_target_comp
			if(logfile$workflow[14]=="TRUE"){# must rerun: screen_target_comp
				logfile$Tasks_to_redo[11]<-TRUE;			
			}
			# must not rerun: profnorm - IS-based normalization
			# must not rerun: homol
			# must not rerun: massdef		
		}
		# IS-Normalization #####################################################
		at1<-logfile$parameters[c(70:78)]
		logfile$parameters[[70]]<-as.character(isolate(input$profnorm_cover_files))
		logfile$parameters[[71]]<-as.character(isolate(input$profnorm_cover_isccount))
		logfile$parameters[[72]]<-as.character(isolate(input$profnorm_use_blank))
		logfile$parameters[[73]]<-as.character(isolate(input$profnorm_use_blank_sample))
		logfile$parameters[[74]]<-as.character(isolate(input$profnorm_use_blank_samplecount))
		logfile$parameters[[75]]<-as.character(isolate(input$profnorm_use_nonblank))
		logfile$parameters[[76]]<-as.character(isolate(input$profnorm_use_nonblank_sample))
		logfile$parameters[[77]]<-as.character(isolate(input$profnorm_use_nonblank_samplecount))	
		logfile$parameters[[78]]<-as.character(isolate(input$profnorm_threshold))
		at2<-logfile$parameters[c(70:78)]
		if(any(is.na(match(at2,at1)))){ 
			# must not rerun: qc
			# must not rerun: recal
			# must not rerun: align
			# must not rerun: norm
			# skip is_pattern
			# skip target_pattern
			if(logfile$workflow[7]=="TRUE"){# must rerun: component isotopologues - but NOT on file level (measurements[,16])!
				logfile$Tasks_to_redo[13]<-TRUE;					
			}
			if(logfile$workflow[8]=="TRUE"){# must rerun: component adducts - but NOT on file level (measurements[,16])!
				logfile$Tasks_to_redo[13]<-TRUE;				
			}
			# must not rerun: profiled
			if(logfile$workflow[10]=="yes"){# must rerun: trendblind
				logfile$Tasks_to_redo[7]<-TRUE;			
			}
			# must not rerun: screen_IS_sam  
			# must not rerun: screen_target_sam
			if(logfile$workflow[13]=="TRUE"){# must rerun: screen_IS_comp
				logfile$Tasks_to_redo[12]<-TRUE;			
			}			
			if(logfile$workflow[14]=="TRUE"){# must rerun: screen_target_comp
				logfile$Tasks_to_redo[11]<-TRUE;			
			}
			if(logfile$workflow[15]=="yes"){# must rerun: profnorm - IS-based normalization
				logfile$Tasks_to_redo[15]<-TRUE;			
			}			
			if(logfile$workflow[16]=="yes"){# must rerun: homol
				logfile$Tasks_to_redo[14]<-TRUE;			
			}			
			if(logfile$workflow[17]=="yes"){# must rerun: massdef
				logfile$Tasks_to_redo[16]<-TRUE;			
			}				
		}		
		# Componentization #####################################################
		at1<-"1"
		at2<-"1"
		if(any(is.na(match(at2,at1)))){ 
			# must not rerun: qc
			# must not rerun: recal
			# must not rerun: align
			# must not rerun: norm
			# skip is_pattern
			# skip target_pattern
			if(logfile$workflow[7]=="TRUE"){ # must rerun: component isotopologues
				logfile$Tasks_to_redo[13]<-TRUE;	
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				measurements[,16]<-FALSE;
				write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);				
			}
			if(logfile$workflow[8]=="TRUE"){ # must rerun: component adducts
				logfile$Tasks_to_redo[13]<-TRUE;
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				measurements[,16]<-FALSE;
				write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);				
			}
			# must not rerun: profiled
			# must not rerun: trendblind
			# must not rerun: screen_IS_sam  
			# must not rerun: screen_target_sam
			if(logfile$workflow[13]=="TRUE"){# must rerun: screen_IS_comp
				logfile$Tasks_to_redo[12]<-TRUE;			
			}			
			if(logfile$workflow[14]=="TRUE"){# must rerun: screen_target_comp
				logfile$Tasks_to_redo[11]<-TRUE;			
			}
			# must not rerun: profnorm - IS-based normalization
			if(logfile$workflow[16]=="yes"){# must rerun: homol
				logfile$Tasks_to_redo[14]<-TRUE;		
			}			
			if(logfile$workflow[17]=="yes"){# must rerun: massdef
				logfile$Tasks_to_redo[16]<-TRUE;		
			}			
		}
		# Homologue series detection ###########################################
		at1<-"1"
		at2<-"1"
		if(any(is.na(match(at2,at1)))){ 
			# must not rerun: qc
			# must not rerun: recal
			# must not rerun: align
			# must not rerun: norm
			# skip is_pattern
			# skip target_pattern
			# must not rerun: component isotopologues
			# must not rerun: component adducts
			# must not rerun: profiled
			# must not rerun: trendblind
			# must not rerun: screen_IS_sam  
			# must not rerun: screen_target_sam
			# must not rerun: screen_IS_comp
			# must not rerun: screen_target_comp
			# must not rerun: profnorm - IS-based normalization
			if(logfile$workflow[16]=="yes"){# must rerun: homol
				logfile$Tasks_to_redo[14]<-TRUE;			
			}			
			# must not rerun: massdef
		}
		########################################################################	
			

		##########################################################################
		# workflow settings ######################################################
		# QC #####################################################################
		at1<-logfile$workflow[1]; 
		logfile$workflow[1]<-as.character(isolate(input$qc));
		at2<-logfile$workflow[1];
		if(at1!=at2){
			if(logfile$workflow[1]=="yes"){# must rerun: qc
				logfile$Tasks_to_redo[2]<-TRUE;			
			}else{
				logfile$Tasks_to_redo[2]<-FALSE; # overwrite a remaining task					
			}			
			if(logfile$workflow[2]=="yes"){# must rerun: recal
				logfile$Tasks_to_redo[3]<-TRUE;			
			}
			if(logfile$workflow[3]=="yes"){# must rerun: align
				logfile$Tasks_to_redo[5]<-TRUE;			
			}
			if(logfile$workflow[4]=="yes"){# must rerun: norm
				logfile$Tasks_to_redo[4]<-TRUE;			
			}
			# skip is_pattern
			# skip target_pattern
			if(logfile$workflow[7]=="TRUE"){# must rerun: component isotopologues
				logfile$Tasks_to_redo[13]<-TRUE;
			}
			if(logfile$workflow[8]=="TRUE"){# must rerun: component adducts
				logfile$Tasks_to_redo[13]<-TRUE;
			}
			if(logfile$workflow[9]=="yes"){# must rerun: profiled
				logfile$Tasks_to_redo[6]<-TRUE;		
			}
			if(logfile$workflow[10]=="yes"){# must rerun: trendblind
				logfile$Tasks_to_redo[7]<-TRUE;			
			}
			if(logfile$workflow[11]=="TRUE"){# must rerun: screen_IS_sam  
				logfile$Tasks_to_redo[10]<-TRUE;		
			}
			if(logfile$workflow[12]=="TRUE"){# must rerun: screen_target_sam
				logfile$Tasks_to_redo[9]<-TRUE;		
			}		
			if(logfile$workflow[13]=="TRUE"){# must rerun: screen_IS_comp
				logfile$Tasks_to_redo[12]<-TRUE;			
			}			
			if(logfile$workflow[14]=="TRUE"){# must rerun: screen_target_comp
				logfile$Tasks_to_redo[11]<-TRUE;		
			}
			if(logfile$workflow[15]=="yes"){# must rerun: profnorm - IS-based normalization
				logfile$Tasks_to_redo[15]<-TRUE;		
			}			
			if(logfile$workflow[16]=="yes"){# must rerun: homol
				logfile$Tasks_to_redo[14]<-TRUE;			
			}			
			if(logfile$workflow[17]=="yes"){# must rerun: massdef
				logfile$Tasks_to_redo[16]<-TRUE;			
			}			
		}
		# recal ##################################################################	
		at1<-logfile$workflow[2];
		logfile$workflow[2]<-as.character(isolate(input$recal));
		at2<-logfile$workflow[2];
		if(at1!=at2){
			if(logfile$workflow[2]=="yes"){# must rerun: recal
				logfile$Tasks_to_redo[3]<-TRUE;			
			}else{
				logfile$Tasks_to_redo[3]<-FALSE;				
			}
			if(logfile$workflow[3]=="yes"){# must rerun: align
				logfile$Tasks_to_redo[5]<-TRUE;			
			}
			if(logfile$workflow[4]=="yes"){# must rerun: norm
				logfile$Tasks_to_redo[4]<-TRUE;			
			}
			# skip is_pattern
			# skip target_pattern
			if(logfile$workflow[7]=="TRUE"){# must rerun: component isotopologues
				logfile$Tasks_to_redo[13]<-TRUE;
			}
			if(logfile$workflow[8]=="TRUE"){# must rerun: component adducts
				logfile$Tasks_to_redo[13]<-TRUE;
			}
			if(logfile$workflow[9]=="yes"){# must rerun: profiled
				logfile$Tasks_to_redo[6]<-TRUE;			
			}
			if(logfile$workflow[10]=="yes"){# must rerun: trendblind
				logfile$Tasks_to_redo[7]<-TRUE;		
			}
			if(logfile$workflow[11]=="TRUE"){# must rerun: screen_IS_sam  
				logfile$Tasks_to_redo[10]<-TRUE;		
			}
			if(logfile$workflow[12]=="TRUE"){# must rerun: screen_target_sam
				logfile$Tasks_to_redo[9]<-TRUE;			
			}		
			if(logfile$workflow[13]=="TRUE"){# must rerun: screen_IS_comp
				logfile$Tasks_to_redo[12]<-TRUE;		
			}			
			if(logfile$workflow[14]=="TRUE"){# must rerun: screen_target_comp
				logfile$Tasks_to_redo[11]<-TRUE;	
			}
			if(logfile$workflow[15]=="yes"){# must rerun: profnorm - IS-based normalization
				logfile$Tasks_to_redo[15]<-TRUE;		
			}			
			if(logfile$workflow[16]=="yes"){# must rerun: homol
				logfile$Tasks_to_redo[14]<-TRUE;			
			}			
			if(logfile$workflow[17]=="yes"){# must rerun: massdef
				logfile$Tasks_to_redo[16]<-TRUE;		
			}			
		}
		# recal ##################################################################	
		# align ##################################################################
		#at1<-logfile[[6]][3]; # align ###########################################
		#logfile[[6]][3]<-as.character(isolate(input$RTalign));
		#at2<-logfile[[6]][3];
		#if(at1!=at2){
		#	logfile$Tasks_to_redo[3:7]<-rep(TRUE,length(3:7));
		#}	
		# norm ###################################################################
		at1<-logfile$workflow[4];
		logfile$workflow[4]<-as.character(isolate(input$intnorm));
		at2<-logfile$workflow[4];
		if(at1!=at2){		
			if(logfile$workflow[4]=="yes"){# must rerun: norm
				logfile$Tasks_to_redo[4]<-TRUE;			
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				measurements[,14]<-FALSE;
				write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
			}else{
				logfile$Tasks_to_redo[4]<-FALSE;		
			}
			# skip is_pattern
			# skip target_pattern
			if(logfile$workflow[7]=="TRUE"){# must rerun: component isotopologues
				logfile$Tasks_to_redo[13]<-TRUE;
			}
			if(logfile$workflow[8]=="TRUE"){# must rerun: component adducts
				logfile$Tasks_to_redo[13]<-TRUE;
			}
			if(logfile$workflow[9]=="yes"){# must rerun: profiled
				logfile$Tasks_to_redo[6]<-TRUE;		
			}
			if(logfile$workflow[10]=="yes"){# must rerun: trendblind
				logfile$Tasks_to_redo[7]<-TRUE;		
			}
			if(logfile$workflow[11]=="TRUE"){# must rerun: screen_IS_sam  
				logfile$Tasks_to_redo[10]<-TRUE;		
			}
			if(logfile$workflow[12]=="TRUE"){# must rerun: screen_target_sam
				logfile$Tasks_to_redo[9]<-TRUE;			
			}		
			if(logfile$workflow[13]=="TRUE"){# must rerun: screen_IS_comp
				logfile$Tasks_to_redo[12]<-TRUE;		
			}			
			if(logfile$workflow[14]=="TRUE"){# must rerun: screen_target_comp
				logfile$Tasks_to_redo[11]<-TRUE;		
			}
			if(logfile$workflow[15]=="yes"){# must rerun: profnorm - IS-based normalization
				logfile$Tasks_to_redo[15]<-TRUE;			
			}			
			if(logfile$workflow[16]=="yes"){# must rerun: homol
				logfile$Tasks_to_redo[14]<-TRUE;			
			}			
			if(logfile$workflow[17]=="yes"){# must rerun: massdef
				logfile$Tasks_to_redo[16]<-TRUE;		
			}			
		}
		# profiling ##############################################################
		at1<-logfile$workflow[9]; 
		logfile$workflow[9]<-as.character(isolate(input$profiled));
		at2<-logfile$workflow[9];
		if(at1!=at2){		
			# skip is_pattern
			# skip target_pattern
			if(logfile$workflow[7]=="TRUE"){# must rerun: component isotopologues
				logfile$Tasks_to_redo[13]<-TRUE;
			}
			if(logfile$workflow[8]=="TRUE"){# must rerun: component adducts
				logfile$Tasks_to_redo[13]<-TRUE;
			}
			if(logfile$workflow[9]=="yes"){# must rerun: profiled
				logfile$Tasks_to_redo[6]<-TRUE;		
			}else{
				logfile$Tasks_to_redo[6]<-FALSE;			
			}
			if(logfile$workflow[10]=="yes"){# must rerun: trendblind
				logfile$Tasks_to_redo[7]<-TRUE;			
			}
			if(logfile$workflow[11]=="TRUE"){# must rerun: screen_IS_sam  
				logfile$Tasks_to_redo[10]<-TRUE;		
			}
			if(logfile$workflow[12]=="TRUE"){# must rerun: screen_target_sam
				logfile$Tasks_to_redo[9]<-TRUE;			
			}		
			if(logfile$workflow[13]=="TRUE"){# must rerun: screen_IS_comp
				logfile$Tasks_to_redo[12]<-TRUE;		
			}			
			if(logfile$workflow[14]=="TRUE"){# must rerun: screen_target_comp
				logfile$Tasks_to_redo[11]<-TRUE;			
			}
			if(logfile$workflow[15]=="yes"){# must rerun: profnorm - IS-based normalization
				logfile$Tasks_to_redo[15]<-TRUE;			
			}			
			if(logfile$workflow[16]=="yes"){# must rerun: homol
				logfile$Tasks_to_redo[14]<-TRUE;		
			}			
			if(logfile$workflow[17]=="yes"){# must rerun: massdef
				logfile$Tasks_to_redo[16]<-TRUE;		
			}			
		}
		# trendblind #############################################################
		at1<-logfile$workflow[10];
		logfile$workflow[10]<-as.character(isolate(input$trenddetect));
		at2<-logfile$workflow[10];
		if(at1!=at2){
			cat("\nchanged trend+blind workflow setting!")
			# skip is_pattern
			# skip target_pattern
			if(logfile$workflow[7]=="TRUE"){# must rerun: component isotopologues
				logfile$Tasks_to_redo[13]<-TRUE;	
			}
			if(logfile$workflow[8]=="TRUE"){# must rerun: component adducts
				logfile$Tasks_to_redo[13]<-TRUE;
			}
			if(logfile$workflow[9]=="yes"){# must rerun: profiled
				logfile$Tasks_to_redo[6]<-TRUE;			
			}
			if(logfile$workflow[10]=="yes"){# must rerun: trendblind
				logfile$Tasks_to_redo[7]<-TRUE;		
			}else{
				logfile$Tasks_to_redo[7]<-FALSE;				
			}
			if(logfile$workflow[11]=="TRUE"){# must rerun: screen_IS_sam  
				logfile$Tasks_to_redo[10]<-TRUE;	
			}
			if(logfile$workflow[12]=="TRUE"){# must rerun: screen_target_sam
				logfile$Tasks_to_redo[9]<-TRUE;			
			}		
			if(logfile$workflow[13]=="TRUE"){# must rerun: screen_IS_comp
				logfile$Tasks_to_redo[12]<-TRUE;			
			}			
			if(logfile$workflow[14]=="TRUE"){# must rerun: screen_target_comp
				logfile$Tasks_to_redo[11]<-TRUE;		
			}
			if(logfile$workflow[15]=="yes"){# must rerun: profnorm - IS-based normalization
				logfile$Tasks_to_redo[15]<-TRUE;		
			}			
			if(logfile$workflow[16]=="yes"){# must rerun: homol
				logfile$Tasks_to_redo[14]<-TRUE;		
			}			
			if(logfile$workflow[17]=="yes"){# must rerun: massdef
				logfile$Tasks_to_redo[16]<-TRUE;		
			}			
		}
		# IS-normalization #######################################################
		at1<-logfile$workflow[15];
		logfile$workflow[15]<-as.character(isolate(input$profnorm));
		at2<-logfile$workflow[15];
		if(at1!=at2){
			# skip is_pattern
			# skip target_pattern
			if(logfile$workflow[7]=="TRUE"){# must rerun: component isotopologues
				logfile$Tasks_to_redo[13]<-TRUE;
			}
			if(logfile$workflow[8]=="TRUE"){# must rerun: component adducts
				logfile$Tasks_to_redo[13]<-TRUE;
			}
			if(logfile$workflow[11]=="TRUE"){# must rerun: screen_IS_sam  
				logfile$Tasks_to_redo[10]<-TRUE;		
			}
			if(logfile$workflow[12]=="TRUE"){# must rerun: screen_target_sam
				logfile$Tasks_to_redo[9]<-TRUE;		
			}		
			if(logfile$workflow[13]=="TRUE"){# must rerun: screen_IS_comp
				logfile$Tasks_to_redo[12]<-TRUE;		
			}			
			if(logfile$workflow[14]=="TRUE"){# must rerun: screen_target_comp
				logfile$Tasks_to_redo[11]<-TRUE;		
			}
			if(logfile$workflow[15]=="yes"){# must rerun: profnorm - IS-based normalization
				logfile$Tasks_to_redo[15]<-TRUE;		
			}else{
				logfile$Tasks_to_redo[15]<-FALSE;			
			}			
			if(logfile$workflow[16]=="yes"){# must rerun: homol
				logfile$Tasks_to_redo[14]<-TRUE;			
			}			
			if(logfile$workflow[17]=="yes"){# must rerun: massdef
				logfile$Tasks_to_redo[16]<-TRUE;		
			}			
		}		
		##########################################################################	
		# set initial workflow dependencies ######################################
		if(logfile$workflow[2]=="yes"){ # recalibration depends...
			logfile$workflow[5]<-"yes"
			logfile$workflow[6]<-"yes"
		}	
		if(logfile$workflow[10]=="yes"){ # trendblind depends... 
			logfile$workflow[9]<-"yes"
		}	
		#if(logfile$workflow[13]=="yes"){ # screen_IS_comp  depends... 
		#	logfile$workflow[9]<-"yes"
		#	logfile$workflow[10]<-"yes"
		#	if( (logfile$workflow[7]!="TRUE") && (logfile$workflow[8]!="TRUE") ){ # use both as default 	
		#		logfile$workflow[7]<-"TRUE"
		#		logfile$workflow[8]<-"TRUE"
		#	}
		#}	
		#if(logfile$workflow[14]=="yes"){ # screen_IS_comp  depends...
		#	logfile$workflow[9]<-"yes"
		#	logfile$workflow[10]<-"yes"
		#	if( (logfile$workflow[7]!="TRUE") && (logfile$workflow[8]!="TRUE") ){ # use both as default 	
		#		logfile$workflow[7]<-"TRUE"
		#		logfile$workflow[8]<-"TRUE"
		#	}
		#}	
		#if(logfile$workflow[15]=="yes"){ # profnorm depends...
		#	logfile$workflow[9]<-"yes"
		#	logfile$workflow[10]<-"yes"
		#}	
		#if(logfile$workflow[16]=="yes"){ # homol depends...
		#	logfile$workflow[9]<-"yes"
		#	logfile$workflow[10]<-"yes"
		#}	
		#if(logfile$workflow[17]=="yes"){ # mass defect depends...
		#	logfile$workflow[9]<-"yes"
		#	logfile$workflow[10]<-"yes"
		#}		
		#if(logfile$workflow[7]=="TRUE"){ # Comp_isotop depends...
		#	logfile$workflow[9]<-"yes"
		#	logfile$workflow[10]<-"yes"
		#}	
		#if(logfile$workflow[8]=="TRUE"){ # Comp_add depends...
		#	logfile$workflow[9]<-"yes"
		#	logfile$workflow[10]<-"yes"
		#}	
   	  ##########################################################################
      save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
      logfile<<-logfile;
      cat("settings changed \n");
      output$dowhat<<-renderText("Project settings modified");
    }

})
################################################################################