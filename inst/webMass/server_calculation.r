mainchecked<-reactive({
    input$Check
    if(input$Check){
		say<-enviMass:::check_project(isotopes,adducts,skipcheck=isolate(input$do_project_check),ignorefiles=isolate(input$ignore_large_files),write_tables=FALSE);
		output$dowhat<<-renderText(say)
		if(say=="Project consistent"){
			cat("Project consistent\n");
			return("Project consistent\n");
		}else{
			cat("Project inconsistent\n");
			shinyjs:::info(say);
			return("Project inconsistent\n");		
		}
	}
})
output$had_checked<-renderText(paste(mainchecked()))  

observe({ # Set calculation counter "do_flow"
    input$Calc
    if(input$Calc){
      do_flow<<-0
	  time_start<<-Sys.time()
    }
})

#observe({ # Run calculations
maincalc<-reactive({
    input$Calc
    if(input$Calc){

		if(do_flow==0){	# check only once, initially at do_flow==0! really?
			say<-enviMass:::check_project(isotopes,adducts,skipcheck=isolate(input$do_project_check),ignorefiles=isolate(input$ignore_large_files),write_tables=FALSE); # because of write_tables=TRUE only here, this check must remain here!
			output$dowhat<<-renderText(say)
			updateSelectInput(session,inputId="Ion_mode_Cal",selected = "none")			
			updateSelectInput(session,inputId="Cal_file_set",selected = "none")
			updateSelectInput(session,inputId="Pos_compound_select",selected = "Choose")
			updateSelectInput(session,inputId="Neg_compound_select",selected = "Choose")
		}else{
			say<-"Project consistent"
		}
      if(say=="Project consistent"){
	  
		if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_calculation.r!")}
        ########################################################################
        # restart logfile[[3]] & mark data availability ########################        
        if(do_flow==0){
			####################################################################
			# erase all previous results #######################################

			# adapt upstream requirements if set to "no" #######################
			must<-logfile[[12]]
			for(i in 1:length(must[1,])){
				for(j in 1:length(must[,1])){	
					if(must[j,i]==1){
						if(logfile$workflow[names(logfile$workflow)==(colnames(must)[i])]=="yes"){
							if(logfile$workflow[names(logfile$workflow)==(rownames(must)[j])]!="yes"){
								cat("\n",colnames(must)[i]," depends on ",rownames(must)[j])
								logfile$workflow[names(logfile$workflow)==(rownames(must)[j])]<<-"yes"
							}
						}
					}
				}
			}
			save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
			####################################################################
			closeAlert(session, alertId="a3")
			logfile$summary[1,2]<<-c(TRUE);
			save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
			summa<<-logfile$summary
			summa[,2]<<-"..."
			output$summa_html<-renderText(enviMass:::summary_html(summa));
        }
		########################################################################
		# run calculations #####################################################
        if(do_flow>0 & do_flow<=length(logfile$summary[,1])){	
			at_node<-logfile$summary[do_flow,1]
			at_script_do<-paste("do_",at_node,".R",sep="")
			if(!file.exists(at_script_do)){
				at_script_do<-FALSE
			}
			at_script_dont<-paste("dont_",at_node,".R",sep="")
			if(!file.exists(at_script_dont)){
				at_script_dont<-FALSE
			}
			enviMass:::workflow_node(
				at_node,at_node,at_node,at_node,
				path_do=at_script_do,
				path_undo=at_script_dont,
				session,output,input
			)  		
		}
		########################################################################
        # make function reiterate ##############################################
        if(do_flow==(length(logfile$summary[,1])+1)){
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
		if(do_flow==(length(logfile$summary[,1])+2)){
			output$summa_html<<-renderText(enviMass:::summary_html(logfile$summary));
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character")
			output$measurements<<-DT::renderDataTable( # in case a node deletes "measurements" and the table output is still waiting
				measurements[,c("ID","Name","Type","Mode","Place","Date","Time","include","profiled","tag1","tag2","tag3","date_end","time_end","ID_2")]
			); 	
		}
        if(do_flow<(length(logfile$summary[,1])+3)){
			invalidateLater(1, session=NULL)
			cat("Calculating...");
			return("Calculating...")
		}else{		
			output$summa_html<<-renderText(enviMass:::summary_html(logfile$summary));
			isolate(init$b<-(init$b+1))
			if(any(ls()=="logfile")){stop("\n illegal logfile detected #2 in server_calculation.r!")}
			cat("Calculations completed, with a \n")
			time_diff<-(Sys.time()-time_start)
			print(time_diff)
			cat("\n")
			return("Calculations completed \n")
		}
        ########################################################################
      }else{
        cat("Project inconsistent\n");
		shinyjs:::info(say);
		return("Project inconsistent\n");
      }
    }
})
output$had_calculated<-renderText(paste(maincalc()))  
