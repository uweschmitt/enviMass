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
      say<-enviMass:::checkproject(logfile,isotopes,adducts);
      output$dowhat<<-renderText(say)
      if(say=="Project consistent"){
	  
        ########################################################################
        # restart logfile[[3]] & mark data availability ########################        
        if(do_flow==1){
			######################################################################
			# erase all previous results #########################################
		  
			
			######################################################################
			closeAlert(session, alertId="a3")
			logfile$summary[1,2]<<-c(TRUE);
			save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
			summa<<-logfile$summary
			summa[1,2]<<-"ok"
			summa[-1,2]<<-"..."
			output$summar<<-renderTable(summa[c(1,2,3,4,5,7,8,9,10),]);
        }
        ########################################################################
        # peak picking - always run ############################################
        if(do_flow==2){
          source("server_calculation_peakpick.R",local=TRUE);    
        }
        ########################################################################
        # quality check ########################################################
        if(do_flow==3){           
          source("server_calculation_qc.R",local=TRUE);    
        }
        ########################################################################
        # isotope pattern calculation ##########################################
        if(do_flow==4){
          source("server_calculation_pattern.R",local=TRUE);            
        }
        ########################################################################
        # m/z recalibration ####################################################
        if(do_flow==5){
          source("server_calculation_recal.R", local=TRUE);    
        }
        ########################################################################
        # RT alignment #########################################################
        if(do_flow==6){
		  cat("Alignment skipped \n");
          output$dowhat<<-renderText("Alignment skipped ... wait")
        }
        ########################################################################
        # intensity normalization ##############################################
        if(do_flow==7){
			source("server_calculation_normaliz.R", local=TRUE);
        }
        ########################################################################
        # profiling ############################################################
        if(do_flow==8){
			source("server_calculation_profiling.R", local=TRUE);
        }
		########################################################################
		# IS-normalization #####################################################
        if(do_flow==9){
			source("server_calculation_IS_normaliz.R", local=TRUE);
		}
        ########################################################################
        # trend / blind ########################################################
        if(do_flow==10){
			source("server_calculation_trendblind.R", local=TRUE);		
        }
        ########################################################################
		# IS file screen #######################################################
        if(do_flow==11){
##		
		}
        ########################################################################
		# target file screen ###################################################
        if(do_flow==12){
##		
		}
        ########################################################################
        # componentization #####################################################
        if(do_flow==13){
##	
        }
        ########################################################################
        # IS comp-screen? ######################################################
        if(do_flow==14){
##	
        }
        ########################################################################
        # Target comp-screen? ##################################################
        if(do_flow==15){
##	
        }
        ########################################################################
        # Homologues ###########################################################
        if(do_flow==16){
##	
        }
        ########################################################################
        # Mass defect ##########################################################
        if(do_flow==17){
##	
        }
        ########################################################################
		
		
        ########################################################################
        # make function reiterate ##############################################
        if(do_flow==19){
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
		if(do_flow==18){output$summar<<-renderTable(logfile$summary)}
        if(do_flow<20){
			invalidateLater(500, session=NULL)
			cat("Calculating...");
			return("Calculating...")
		}else{		
			output$summar<<-renderTable(summa[c(1,2,3,4,5,7,8,9,10),]);
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
