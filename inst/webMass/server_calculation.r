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

	  say<-enviMass:::checkproject(isotopes,adducts,skipcheck=isolate(input$do_project_check));
	  output$dowhat<<-renderText(say)
      if(say=="Project consistent"){
	  
		if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_calculation.r!")}
        ########################################################################
        # restart logfile[[3]] & mark data availability ########################        
        if(do_flow==1){
			####################################################################
			# erase all previous results #######################################
		  
			
			# adapt upstream requirements if set to "no" #######################
			must<-logfile[[12]]
			for(i in 1:length(must[1,])){
				for(j in 1:length(must[,1])){	
					if(must[j,i]==1){
						if(logfile$workflow[names(logfile$workflow)==(colnames(must)[i])]=="yes"){
							if(logfile$workflow[names(logfile$workflow)==(rownames(must)[j])]!="yes"){
								#cat("\n",colnames(must)[i]," depends on ",rownames(must)[j])
								logfile$workflow[names(logfile$workflow)==(rownames(must)[j])]<-"yes"
							}
						}
					}
				}
			}
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
				name_workflow="peakpicking",
				name_summary="peakpicking",
				name_redo="peakpicking",
				"Peak picking",
				path_do="do_peakpicking.R",path_undo=FALSE,session,output,input
			)  
        }
        ########################################################################
        # quality check ########################################################
        if(do_flow==3){           
          	workflow_node(
				"qc","qc","qc","Quality control",
				path_do="do_qc.R",path_undo="dont_qc.R",session,output,input
			)   
        }
        ########################################################################
        # isotope pattern calculation ##########################################
        if(do_flow==4){
		    workflow_node(
				"pattern","pattern","pattern","Compound isotope patterns",
				path_do="do_pattern.R",path_undo=FALSE,session,output,input
			)  		            
        }
        ########################################################################
        # m/z recalibration ####################################################
        if(do_flow==5){
			workflow_node(
				"recal","recal","recal","Mass recalibration",
				path_do="do_recal.R",path_undo="dont_recal.R",session,output,input
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
				"blinds","blinds","blinds","Blind filter",
				path_do="do_blind.R",path_undo="dont_blind.R",session,output,input
			)  	
		}			
        ########################################################################
		# Replicates ###########################################################
        if(do_flow==8){
			workflow_node(
				"replicates","replicates","replicates","Replicate filter",
				path_do="do_replicates.R",path_undo="dont_replicates.R",session,output,input
			)  	
		}			
        ########################################################################
        # intensity normalization ##############################################
        if(do_flow==9){
			workflow_node(
				"norm","norm","norm","Median intensity normalization",
				path_do="do_norm.R",path_undo="dont_norm.R",session,output,input
			)  
        }
        ########################################################################
        # profiling ############################################################
        if(do_flow==10){
			workflow_node(
				"profiling","profiling","profiling","Profile extraction",
				path_do="do_profiling.R",path_undo="dont_profiling.R",session,output,input
			)  
        }
        ########################################################################
		# LOD ##################################################################
        if(do_flow==11){
			workflow_node(
				"LOD","LOD","LOD","LOD interpolation",
				path_do="do_LOD.R",path_undo="dont_LOD.R",session,output,input
			)  	
		}				
        ########################################################################
		# IS screening #########################################################
        if(do_flow==12){
			workflow_node(
				"IS_screen","IS_screen","IS_screen","IS screening",
				path_do="do_IS_screen.R",path_undo="dont_IS_screen.R",session,output,input
			)  	
		}		
        ########################################################################
		# target screening #####################################################
        if(do_flow==13){
			workflow_node(
				"target_screen","target_screen","target_screen","Target screening",
				path_do="do_target_screen.R",path_undo="dont_target_screen.R",session,output,input
			)  		
		}			
		########################################################################
		# IS-normalization #####################################################
        if(do_flow==14){
			workflow_node(
				"IS_normaliz","IS_normaliz","IS_normaliz","IS-based intensity normalization",
				path_do="do_IS_normaliz.R",path_undo="dont_IS_normaliz.R",session,output,input
			)  
		}
		########################################################################
		# Quantification #######################################################
        if(do_flow==15){
			workflow_node(
				"quantification","quantification","quantification","Quantification",
				path_do="do_quantification.R",path_undo="dont_quantification.R",session,output,input
			)  	
		}		
        ########################################################################
        # trend / blind ########################################################
        if(do_flow==16){
			workflow_node(
				"trendblind","trendblind","trendblind","Trend detection and blind subtraction",
				path_do="do_trendblind.R",path_undo="dont_trendblind.R",session,output,input
			)  	
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
			logfile$Tasks_to_redo[1:length(logfile$Tasks_to_redo)]<<-"FALSE"	
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
			if(any(ls()=="logfile")){stop("\n illegal logfile detected #2 in server_calculation.r!")}
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
