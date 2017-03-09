# write variables to logfile ###################################################
# observe what has changed and reset Tasks_to_do accordingly ###################

observe({
    input$savepar;
    input$saveflow;
    if(	
		(exists("logfile")) & ((isolate(input$savepar)) || (isolate(input$saveflow)))
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
		# checking parameter validity ##########################################
		all_ok<-TRUE
		for(i in 1:length(logfile$parameters)){
			for_param<-names(logfile$parameters)[i]
			if(any(names(input)==for_param)){
				found_it<-try(eval(parse(text=paste("new_param<-","as.character(isolate(input$",for_param,"))",sep=""))))
				if(class(found_it)=="try-error"){	# present in UI?
					cat("\n Error for parameter input detected: ");cat(for_param);					
					shinyjs:::info(paste("Error on parameter input detected for ",for_param,". Comma instead of dot-seperated? Or vice versa?",sep=""))
					all_ok<-FALSE
				}else{
					found_it<-eval(parse(text=paste("new_param<-","as.character(isolate(input$",for_param,"))",sep="")))
					if(is.na(found_it) || length(found_it)==0){
						cat("\n Invalid parameter input detected: ");cat(for_param);					
						shinyjs:::info(paste("Invalid parameter input detected for ",for_param,". Comma instead of dot-seperated? Or vice versa?",sep=""))
						all_ok<-FALSE
					}
				}
			}		
		}
		######################################################################## 		
		
		if(all_ok){ # otherwise report problem to try() in server_calculations.r
			######################################################################## 						
			# parameter settings ###################################################
			cat("\n Checking parameter changes ... ")
			found<-c()
			# check standard parameters defined in function newproject()
			for(i in 1:length(logfile$parameters)){	
				if(names(logfile$parameters)[i]==""){next} # for any empty entries
				if(names(logfile$parameters)[i]=="external"){next} # for external parameters		
				old_param<-logfile$parameters[i]		
				for_param<-names(logfile$parameters)[i]
				eval(parse(text=paste("new_param<-","as.character(isolate(input$",for_param,"))",sep="")))
				if( length(new_param)>0 ){ # in shiny input list?
					if(!enviMass:::comp_list(old_param,new_param,as_pairs=FALSE)){
							# report changed parameter
							cat(paste("\n","Changed parameter ",for_param," from ",sep=""));
							cat("\n");print(old_param);cat("\n");cat(" to ");cat("\n");print(new_param);		
							# find affected files
							found<-c(found,
								paste("logfile$parameters$",names(logfile$parameters)[i],sep="")
							)
							logfile$parameters[i]<<-new_param
					}else{
						cat(paste(".",sep=""))
					}		
				}else{
					cat(paste("\n",for_param," not found in input list. ",sep=""))
				}
			}
			# check external parameters defined in workflow_parameters.r
			external_old<-logfile$parameters$external
			source(file="workflow_parameters.r")
			if(length(logfile$parameters$external)){
				for(i in 1:length(logfile$parameters$external)){	
					if(any(names(external_old)==names(logfile$parameters$external)[i])){
						this<-which(names(external_old)==names(logfile$parameters$external)[i])
						old_param<-external_old[[this]]
						new_param<-logfile$parameters$external[[i]]
						if(!enviMass:::comp_list(old_param,new_param,as_pairs=FALSE)){
							cat("*")
							found<-c(found,paste("logfile$parameters$external$",names(logfile$parameters$external)[i],sep=""))
						}
					}else{ # seems to be a new parameter
						found<-c(found,paste("logfile$parameters$external$",names(logfile$parameters$external)[i],sep=""))				
					}
				}
			}
			# any parameter changed ?
			if(any(found=="logfile$parameters$progressBar")){
				found<-found[found!="logfile$parameters$progressBar"]
			}
			if(length(found)){
				affected<-enviMass:::workflow_where(found)
				if(length(affected)>0){
					cat("\nAdapt settings affecting nodes: ")
					print(affected);cat("\n")
					for(i in 1:length(affected)){
						enviMass:::workflow_set(
							down=affected[i],
							check_node=TRUE, 	# only applied if affected step switched on - only in 
							single_file=TRUE,
							single_node=TRUE
						)
					}
				}
			}
			# update !single_file entry changes in measurements
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			output$measurements<<-DT::renderDataTable(
				measurements[,c("ID","Name","Type","Mode","Place","Date","Time","include","profiled","tag1","tag2","tag3","date_end","time_end","ID_2")]
			); 
			cat("Done checking of parameter changes.\n")
			######################################################################## 		

			########################################################################     
			# PW path ##############################################################
			logfile[[4]]<<-as.character(isolate(input$PWpath));
			########################################################################     

			######################################################################## 
			# screening - adducts ##################################################		
			at3<-logfile$adducts_pos
			logfile$adducts_pos<<-as.character(isolate(input$adducts_pos))
			at4<-logfile$adducts_pos
			at5<-logfile$adducts_neg
			logfile$adducts_neg<<-as.character(isolate(input$adducts_neg))
			at6<-logfile$adducts_neg
			if( any(is.na(match(at3,at4))) || any(is.na(match(at4,at3))) || any(is.na(match(at5,at6))) || any(is.na(match(at6,at5))) ){ 
				enviMass:::workflow_set(
					down="pattern",
					check_node=TRUE,
					single_file=TRUE,
					single_node=TRUE
				)
			}
			if(do_debug){cat("\n at_3")}
			######################################################################## 

			######################################################################## 
			# nontarget grouping - adducts #########################################		
			at3<-logfile$adducts_pos_group
			logfile$adducts_pos_group<<-as.character(isolate(input$adducts_pos_group))
			at4<-logfile$adducts_pos_group
			at5<-logfile$adducts_neg_group
			logfile$adducts_neg_group<<-as.character(isolate(input$adducts_neg_group))
			at6<-logfile$adducts_neg_group
			if( any(is.na(match(at3,at4))) || any(is.na(match(at4,at3))) || any(is.na(match(at5,at6))) || any(is.na(match(at6,at5))) ){ 
				enviMass:::workflow_set(
					down="adducts",
					check_node=TRUE,
					single_file=TRUE,
					single_node=TRUE
				)
			}
			if(do_debug){cat("\n at_4")}
			########################################################################
			
			######################################################################## 		
			# subtraction files ####################################################
			at1<-logfile$Positive_subtraction_files
			logfile$Positive_subtraction_files<<-c(isolate(input$files_pos_select_subtract),"FALSE")
			at2<-logfile$Positive_subtraction_files
			if(any(is.na(at2))){stop("\nThere was an issue reading out the new settings - maybe comma / dot separation was not fullfilled?")}		
			if(!enviMass:::comp_list(at1,at2,as_pairs=FALSE)){ # both steps take partly the same parameters! 
				enviMass:::workflow_set(
					down="blind",
					check_node=TRUE,
					single_file=TRUE,
					single_node=TRUE
				)
			}	
			at1<-logfile$Negative_subtraction_files
			logfile$Negative_subtraction_files<<-c(isolate(input$files_neg_select_subtract),"FALSE")
			at2<-logfile$Negative_subtraction_files
			if(any(is.na(at2))){stop("\nThere was an issue reading out the new settings - maybe comma / dot separation was not fullfilled?")}		
			if(!enviMass:::comp_list(at1,at2,as_pairs=FALSE)){ # both steps take partly the same parameters! 
				enviMass:::workflow_set(
					down="blind",
					check_node=TRUE,
					single_file=TRUE,
					single_node=TRUE
				)
			}		
			if(do_debug){cat("\n at_8b")}
			########################################################################  		  
			
			########################################################################
			# workflow settings ####################################################
			cat("\n Checking workflow changes ... ")
			found<-c()
			for(i in 1:length(logfile$workflow)){	
				if(names(logfile$workflow)[i]==""){next} # for any empty entries
				old_node<-logfile$workflow[i]
				for_node<-names(logfile$workflow)[i]
				eval(parse(text=paste("new_node<-","as.character(isolate(input$",for_node,"))",sep="")))
				if( length(new_node)>0 ){ # in shiny input list?
					if(!enviMass:::comp_list(old_node,new_node,as_pairs=TRUE)){
							# report changed parameter
							cat(paste("\n","Changed parameter ",for_node," from ",sep=""));
							cat("\n");print(old_node);cat("\n");cat(" to ");cat("\n");print(new_node);		
							# find affected files
							found<-c(found,for_node)
							logfile$workflow[i]<<-new_node
					}else{
						cat(paste(".",sep=""))
					}		
				}else{
					cat(paste("\n",for_node," not found in input list.\n",sep=""))
				}
			}
			# any parameter changed ?
			if(length(found)){
					cat("\nAdapt settings affecting nodes: ")
					print(found);cat("\n")
					for(i in 1:length(found)){
						enviMass:::workflow_set(
							down=found[i],
							check_node=FALSE,
							single_file=TRUE,
							single_node=FALSE
						) # do not change
					}
			}
			cat("Done checking of workflow changes.\n")
			########################################################################
		
			########################################################################
			save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
			cat("settings changed \n");
			output$dowhat<<-renderText("Project settings modified");
			########################################################################
			if(any(ls()=="logfile")){stop("\n illegal logfile detected #2 in server_variable_out.r!")}
			if(do_debug){cat("\n at_27")}
		}
    }

})
################################################################################
