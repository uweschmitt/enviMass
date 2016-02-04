 observe({ # Set calculation counter "do_flow"
    input$Calc
    if(input$Calc){
      do_flow<<-1
    }
})


#observe({ # Run calculations
maincalc<-reactive({
    input$Calc
    if(input$Calc){
      say<-enviMass:::checkproject(logfile,isotopes,adducts,skipcheck=isolate(input$do_project_check));
      output$dowhat<<-renderText(say)
      if(say=="Project consistent"){
	  
        ########################################################################
        # restart logfile[[3]] & mark data availability ########################        
        if(do_flow==1){
			####################################################################
			# erase all previous results #######################################
		  
			
			####################################################################
			closeAlert(session, alertId="a3")
			logfile$summary[1,2]<<-c(TRUE);
			save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
			summa<<-logfile$summary
			summa[1,2]<<-"ok"
			summa[-1,2]<<-"..."
			output$summa_html<-renderText(summary_html(summa));
        }
        ########################################################################
        # peak picking - always run ############################################
        if(do_flow==2){
			workflow_node(
				"peakpicking","Peak pick?","peakpick","Peak picking",
				path_do="do_peakpick.R",path_undo=FALSE,logfile,output,input
			)  
        }
        ########################################################################
        # quality check ########################################################
        if(do_flow==3){           
          	workflow_node(
				"qc","QC?","QC","Quality control",
				path_do="do_qc.R",path_undo="dont_qc.R",logfile,output,input
			)   
        }
        ########################################################################
        # isotope pattern calculation ##########################################
        if(do_flow==4){
		    workflow_node(
				"pattern","Isotope pattern?","pattern","Compound isotope patterns",
				path_do="do_pattern.R",path_undo=FALSE,logfile,output,input
			)  		            
        }
        ########################################################################
        # m/z recalibration ####################################################
        if(do_flow==5){
			workflow_node(
				"recal","m/z recal.?","recal","Mass recalibration",
				path_do="do_recal.R",path_undo="dont_recal.R",logfile,output,input
			)     
        }
        ########################################################################
        # RT alignment #########################################################
        if(do_flow==6){
		  cat("Alignment skipped \n");
          output$dowhat<<-renderText("Alignment skipped ... wait");
        }
        ########################################################################
		# Replicates ###########################################################
        if(do_flow==7){
			workflow_node(
				"replicates","Replicate filter","replicates","Replicate filter",
				path_do="do_replicates.R",path_undo="dont_replicates.R",logfile,output,input
			)  	
		}			
        ########################################################################
        # intensity normalization ##############################################
        if(do_flow==8){
			workflow_node(
				"norm","Intensity norm.?","normalize","Median intensity normalization",
				path_do="do_normaliz.R",path_undo="dont_normaliz.R",logfile,output,input
			)  
        }
        ########################################################################
        # profiling ############################################################
        if(do_flow==9){
			workflow_node(
				"profiled","Profiled?","profiling","Profile extraction",
				path_do="do_profiling.R",path_undo="dont_profiling.R",logfile,output,input
			)  
        }
        ########################################################################
		# IS screening #########################################################
        if(do_flow==10){
			#source("server_calculation_screening_IS.R", local=TRUE);		
		}		
        ########################################################################
		# target screening #####################################################
        if(do_flow==11){
			#source("server_calculation_screening_target.R", local=TRUE);		
		}			
		########################################################################
		# IS-normalization #####################################################
        if(do_flow==12){
			workflow_node(
				"profnorm","IS norm.?","norm_prof","IS-based intensity normalization",
				path_do="do_IS_normaliz.R",path_undo="dont_IS_normaliz.R",logfile,output,input
			)  
		}
        ########################################################################
        # trend / blind ########################################################
        if(do_flow==13){
			workflow_node(
				"trenddetect","Trend+Blind?","trendblind","Trend detection and blind subtraction",
				path_do="do_trendblind.R",path_undo="dont_trendblind.R",logfile,output,input
			)  	
        }
        ########################################################################
        # 
        if(do_flow==14){
##	
        }
        ########################################################################
        # 
        if(do_flow==15){
##	
        }
        ########################################################################
        # 
        if(do_flow==16){
##	
        }
        ########################################################################
        # 
        if(do_flow==17){
##	
        }
        ########################################################################
		
		
        ########################################################################
        # make function reiterate ##############################################
        if(do_flow==18){
			if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_pos")){rm(profileList_pos,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="profileList_pos")){rm(profileList_pos)}
			if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_neg")){rm(profileList_neg,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="profileList_neg")){rm(profileList_neg)}
			if(any(objects(envir=as.environment(".GlobalEnv"))=="profpeaks_pos")){rm(profpeaks_pos,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="profpeaks_pos")){rm(profpeaks_pos)}
			if(any(objects(envir=as.environment(".GlobalEnv"))=="profpeaks_neg")){rm(profpeaks_neg,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="profpeaks_neg")){rm(profpeaks_neg)}
			logfile$Tasks_to_redo[1:length(logfile$Tasks_to_redo)]<-"FALSE"	
			save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));      
			if(file.exists(file.path(as.character(logfile[[1]]),"results","profileList_pos"))){
				load(file.path(as.character(logfile[[1]]),"results","profileList_pos"),envir=as.environment(".GlobalEnv"),verbose=TRUE)
			}
			if(file.exists(file.path(as.character(logfile[[1]]),"results","profileList_neg"))){
				load(file.path(as.character(logfile[[1]]),"results","profileList_neg"),envir=as.environment(".GlobalEnv"),verbose=TRUE)
			}
			if(file.exists(file.path(as.character(logfile[[1]]),"results","profpeaks_pos"))){
				load(file.path(as.character(logfile[[1]]),"results","profpeaks_pos"),envir=as.environment(".GlobalEnv"),verbose=TRUE)
			}
			if(file.exists(file.path(as.character(logfile[[1]]),"results","profpeaks_neg"))){
				load(file.path(as.character(logfile[[1]]),"results","profpeaks_neg"),envir=as.environment(".GlobalEnv"),verbose=TRUE)
			}
        }
        do_flow<<-(do_flow+1);
		if(do_flow==19){
			output$summa_html<<-renderText(summary_html(logfile$summary));
		}
        if(do_flow<20){
			invalidateLater(500, session=NULL)
			cat("Calculating...");
			return("Calculating...")
		}else{		
			output$summa_html<<-renderText(summary_html(logfile$summary));
			isolate(init$b<-(init$b+1))
			cat("Calculations completed \n")
			return("Calculations completed \n")
		}
        ########################################################################
      }else{
        cat("Project inconsistent\n");
		return("Project inconsistent\n")
      }
    }
})
output$had_calculated<-renderText(paste(maincalc()))  
