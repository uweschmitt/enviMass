##############################################################################
# update compounds, measurements, etc ########################################
##############################################################################
  
# ADD IS #####################################################################
observe({ # update selectable adducts by ionization mode 
	input$ISadd_charge
	if(isolate(input$ISadd_charge)=="positive"){
		updateSelectInput(session, "ISadd_add", "Main adduct (+):", choices = c("FALSE",as.character(adducts[adducts[,6]=="positive",1])), selected="FALSE")
	}
	if(isolate(input$ISadd_charge)=="negative"){
		updateSelectInput(session, "ISadd_add", "Main adduct (-):", choices = c("FALSE",as.character(adducts[adducts[,6]=="negative",1])), selected="FALSE")
	}  
})
observe({
    input$AddIS
    if(input$AddIS){
		IS1<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character")
		IS2<-rep("FALSE",length=16)
		IS2[1]<-as.character(isolate(input$ISadd_ID))
		IS2[2]<-as.character(isolate(input$ISadd_name))
		IS2[3]<-as.character(isolate(input$ISadd_formula))
		IS2[4]<-as.character(isolate(input$ISadd_RT))
		if(isolate(input$ISadd_RTtol_use)){
			IS2[5]<-as.character(isolate(input$ISadd_RTtol))
		}
		IS2[6]<-as.character(isolate(input$ISadd_add))
		IS2[7]<-as.character(isolate(input$ISadd_charge))
		IS2[8]<-as.character(isolate(input$ISadd_use_recal))
		IS2[9]<-as.character(isolate(input$ISadd_use_screen))
		IS2[10]<-as.character(isolate(input$ISadd_rest_adduct))
		IS2[11]<-as.character(isolate(input$ISadd_remark))
		IS2[12]<-as.character(isolate(input$ISadd_tag1))
		IS2[13]<-as.character(isolate(input$ISadd_tag2))
		IS2[14]<-as.character(isolate(input$ISadd_tag3))
		if(isolate(input$ISadd_date)){
			IS2[15]<-as.character(isolate(input$ISadd_date_range[1]))
			IS2[16]<-as.character(isolate(input$ISadd_date_range[2]))
		}
		IS<-rbind(IS2,IS1);
		write.table(IS,file=file.path(logfile[[1]],"dataframes","IS.txt"),row.names=FALSE,sep="\t",quote=FALSE)
		rm(IS,IS1,IS2);
		#############################################################################
		# adjust task/workflow settings #############################################
		if( 
			(logfile$workflow[2]=="yes") || 
			(logfile$workflow[5]=="yes") ||
			(logfile$workflow[11]=="TRUE") ||
			(logfile$workflow[13]=="TRUE") ||
			(logfile$workflow[15]=="yes")
		){	# must rerun: 
			if(logfile$workflow[2]=="yes"){# must rerun: recal
				logfile$Tasks_to_redo[3]<-TRUE;
				logfile$Tasks_to_redo[3]<<-TRUE;				
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				measurements[,12]<-FALSE;
				write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);				
			}
			if(logfile$workflow[3]=="yes"){# must rerun: align
				logfile$Tasks_to_redo[5]<-TRUE;
				logfile$Tasks_to_redo[5]<<-TRUE;				
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				measurements[,13]<-FALSE;
				write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
			}
			# must not rerun: norm
			# must rerun: is_pattern / target_pattern
			logfile$Tasks_to_redo[8]<-TRUE;
			logfile$Tasks_to_redo[8]<<-TRUE;			
			if(logfile$workflow[9]=="yes"){# must rerun: profiled
				logfile$Tasks_to_redo[6]<-TRUE;
				logfile$Tasks_to_redo[6]<<-TRUE;				
			}
			if(logfile$workflow[10]=="yes"){# must rerun: trendblind
				logfile$Tasks_to_redo[7]<-TRUE;
				logfile$Tasks_to_redo[7]<<-TRUE;				
			}
			if(logfile$workflow[11]=="TRUE"){# must rerun: screen_IS_sam  
				logfile$Tasks_to_redo[10]<-TRUE;
				logfile$Tasks_to_redo[10]<<-TRUE;				
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				measurements[,17]<-FALSE;
				write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);												
			}
			if(logfile$workflow[12]=="TRUE"){# must rerun: screen_target_sam
				logfile$Tasks_to_redo[9]<-TRUE;
				logfile$Tasks_to_redo[9]<<-TRUE;				
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				measurements[,18]<-FALSE;
				write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);												
			}		
			if(logfile$workflow[13]=="TRUE"){# must rerun: screen_IS_comp
				logfile$Tasks_to_redo[12]<-TRUE;
				logfile$Tasks_to_redo[12]<<-TRUE;				
			}			
			if(logfile$workflow[14]=="TRUE"){# must rerun: screen_target_comp
				logfile$Tasks_to_redo[11]<-TRUE;
				logfile$Tasks_to_redo[11]<<-TRUE;				
			}
			if(logfile$workflow[15]=="yes"){# must rerun: profnorm - IS-based normalization
				logfile$Tasks_to_redo[15]<-TRUE;
				logfile$Tasks_to_redo[15]<<-TRUE;				
			}			
			if(logfile$workflow[16]=="yes"){# must rerun: homol
				logfile$Tasks_to_redo[14]<-TRUE;
				logfile$Tasks_to_redo[14]<<-TRUE;				
			}			
			if(logfile$workflow[17]=="yes"){# must rerun: massdef
				logfile$Tasks_to_redo[16]<-TRUE;
				logfile$Tasks_to_redo[16]<<-TRUE;				
			}			
		}	  
		#############################################################################			
		output$IS<<-DT::renderDataTable(read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character"));
		logfile[[2]][3:7]<-rep(TRUE,length(3:7));
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));      
		output$dowhat<-renderText("Added IS compound");
    }
})

# ADD IS LIST ################################################################
observe({
 	input$ISlist_path
	if(  (length(isolate(input$ISlist_path))) ){
		if( file.exists(as.character(isolate(input$ISlist_path[[4]]))) ){ 
			IS<-read.table(file=as.character(isolate(input$ISlist_path[[4]])),header=TRUE,sep="\t",colClasses = "character")
			write.table(IS,file=file.path(logfile[[1]],"dataframes","IS.txt"),row.names=FALSE,sep="\t",quote=FALSE)
			rm(IS)
			#############################################################################
			# adjust task/workflow settings #############################################
			if( 
				(logfile$workflow[2]=="yes") || 
				(logfile$workflow[5]=="yes") ||
				(logfile$workflow[11]=="TRUE") ||
				(logfile$workflow[13]=="TRUE") ||
				(logfile$workflow[15]=="yes")
			){	# must rerun: 
				if(logfile$workflow[2]=="yes"){# must rerun: recal
					logfile$Tasks_to_redo[3]<-TRUE;
					logfile$Tasks_to_redo[3]<<-TRUE;				
					measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
					measurements[,12]<-FALSE;
					write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);				
				}
				if(logfile$workflow[3]=="yes"){# must rerun: align
					logfile$Tasks_to_redo[5]<-TRUE;
					logfile$Tasks_to_redo[5]<<-TRUE;				
					measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
					measurements[,13]<-FALSE;
					write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
				}
				# must not rerun: norm
				# must rerun: is_pattern / target_pattern
				logfile$Tasks_to_redo[8]<-TRUE;
				logfile$Tasks_to_redo[8]<<-TRUE;			
				if(logfile$workflow[9]=="yes"){# must rerun: profiled
					logfile$Tasks_to_redo[6]<-TRUE;
					logfile$Tasks_to_redo[6]<<-TRUE;				
				}
				if(logfile$workflow[10]=="yes"){# must rerun: trendblind
					logfile$Tasks_to_redo[7]<-TRUE;
					logfile$Tasks_to_redo[7]<<-TRUE;				
				}
				if(logfile$workflow[11]=="TRUE"){# must rerun: screen_IS_sam  
					logfile$Tasks_to_redo[10]<-TRUE;
					logfile$Tasks_to_redo[10]<<-TRUE;				
					measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
					measurements[,17]<-FALSE;
					write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);												
				}
				if(logfile$workflow[12]=="TRUE"){# must rerun: screen_target_sam
					logfile$Tasks_to_redo[9]<-TRUE;
					logfile$Tasks_to_redo[9]<<-TRUE;				
					measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
					measurements[,18]<-FALSE;
					write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);												
				}		
				if(logfile$workflow[13]=="TRUE"){# must rerun: screen_IS_comp
					logfile$Tasks_to_redo[12]<-TRUE;
					logfile$Tasks_to_redo[12]<<-TRUE;				
				}			
				if(logfile$workflow[14]=="TRUE"){# must rerun: screen_target_comp
					logfile$Tasks_to_redo[11]<-TRUE;
					logfile$Tasks_to_redo[11]<<-TRUE;				
				}
				if(logfile$workflow[15]=="yes"){# must rerun: profnorm - IS-based normalization
					logfile$Tasks_to_redo[15]<-TRUE;
					logfile$Tasks_to_redo[15]<<-TRUE;				
				}			
				if(logfile$workflow[16]=="yes"){# must rerun: homol
					logfile$Tasks_to_redo[14]<-TRUE;
					logfile$Tasks_to_redo[14]<<-TRUE;				
				}			
				if(logfile$workflow[17]=="yes"){# must rerun: massdef
					logfile$Tasks_to_redo[16]<-TRUE;
					logfile$Tasks_to_redo[16]<<-TRUE;				
				}			
			}
			####################################################################
			save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));   
 			output$IS<<-DT::renderDataTable(read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character"));
			output$dowhat<-renderText("Added IS list");
		}
	}
})  
##############################################################################	
  
# DELETE IS ##################################################################
observe({
    input$DeleteIS
    if(input$DeleteIS){
		IS<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character")
		IS<-IS[IS$ID!=as.character(isolate(input$ISdelete_ID)),]
		write.table(IS,file=file.path(logfile[[1]],"dataframes","IS.txt"),row.names=FALSE,sep="\t",quote=FALSE)
		#############################################################################
		# adjust task/workflow settings #############################################
		if( 
				(logfile$workflow[2]=="yes") || 
				(logfile$workflow[5]=="yes") ||
				(logfile$workflow[11]=="TRUE") ||
				(logfile$workflow[13]=="TRUE") ||
				(logfile$workflow[15]=="yes")
		){	# must rerun: 
				if(logfile$workflow[2]=="yes"){# must rerun: recal
					logfile$Tasks_to_redo[3]<-TRUE;
					logfile$Tasks_to_redo[3]<<-TRUE;				
					measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
					measurements[,12]<-FALSE;
					write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);				
				}
				if(logfile$workflow[3]=="yes"){# must rerun: align
					logfile$Tasks_to_redo[5]<-TRUE;
					logfile$Tasks_to_redo[5]<<-TRUE;				
					measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
					measurements[,13]<-FALSE;
					write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
				}
				# must not rerun: norm
				# must rerun: is_pattern / target_pattern
				logfile$Tasks_to_redo[8]<-TRUE;
				logfile$Tasks_to_redo[8]<<-TRUE;			
				if(logfile$workflow[9]=="yes"){# must rerun: profiled
					logfile$Tasks_to_redo[6]<-TRUE;
					logfile$Tasks_to_redo[6]<<-TRUE;				
				}
				if(logfile$workflow[10]=="yes"){# must rerun: trendblind
					logfile$Tasks_to_redo[7]<-TRUE;
					logfile$Tasks_to_redo[7]<<-TRUE;				
				}
				if(logfile$workflow[11]=="TRUE"){# must rerun: screen_IS_sam  
					logfile$Tasks_to_redo[10]<-TRUE;
					logfile$Tasks_to_redo[10]<<-TRUE;				
					measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
					measurements[,17]<-FALSE;
					write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);												
				}
				if(logfile$workflow[12]=="TRUE"){# must rerun: screen_target_sam
					logfile$Tasks_to_redo[9]<-TRUE;
					logfile$Tasks_to_redo[9]<<-TRUE;				
					measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
					measurements[,18]<-FALSE;
					write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);												
				}		
				if(logfile$workflow[13]=="TRUE"){# must rerun: screen_IS_comp
					logfile$Tasks_to_redo[12]<-TRUE;
					logfile$Tasks_to_redo[12]<<-TRUE;				
				}			
				if(logfile$workflow[14]=="TRUE"){# must rerun: screen_target_comp
					logfile$Tasks_to_redo[11]<-TRUE;
					logfile$Tasks_to_redo[11]<<-TRUE;				
				}
				if(logfile$workflow[15]=="yes"){# must rerun: profnorm - IS-based normalization
					logfile$Tasks_to_redo[15]<-TRUE;
					logfile$Tasks_to_redo[15]<<-TRUE;				
				}			
				if(logfile$workflow[16]=="yes"){# must rerun: homol
					logfile$Tasks_to_redo[14]<-TRUE;
					logfile$Tasks_to_redo[14]<<-TRUE;				
				}			
				if(logfile$workflow[17]=="yes"){# must rerun: massdef
					logfile$Tasks_to_redo[16]<-TRUE;
					logfile$Tasks_to_redo[16]<<-TRUE;				
				}			
		}
		####################################################################
		output$IS<<-DT::renderDataTable(read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character"));
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));      
		output$dowhat<-renderText("Deleted compound");
    }
})
  
# ADD TARGET #################################################################
observe({ # update selectable adducts by ionization mode 
	input$targetsadd_charge
	if(isolate(input$targetsadd_charge)=="positive"){
		updateSelectInput(session, "targetsadd_add", "Main adduct (+):", choices = c("FALSE",as.character(adducts[adducts[,6]=="positive",1])), selected="FALSE")
	}
	if(isolate(input$targetsadd_charge)=="negative"){
		updateSelectInput(session, "targetsadd_add", "Main adduct (-):", choices = c("FALSE",as.character(adducts[adducts[,6]=="negative",1])), selected="FALSE")
	}
})
observe({
    input$Addtargets
    if(input$Addtargets){
		targets1<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character")
		targets2<-rep("FALSE",length=19)
		targets2[1]<-as.character(isolate(input$targetsadd_ID))
		targets2[2]<-as.character(isolate(input$targetsadd_name))
		targets2[3]<-as.character(isolate(input$targetsadd_formula))
		targets2[4]<-as.character(isolate(input$targetsadd_RT))
		if(isolate(input$targetsadd_RTtol_use)){
			targets2[5]<-as.character(isolate(input$targetsadd_RTtol))
		}
		#  targets2[6]<-  .... some IS ID
		targets2[7]<-as.character(isolate(input$targetsadd_add))
		targets2[8]<-as.character(isolate(input$targetsadd_charge))
		targets2[9]<-as.character(isolate(input$targetsadd_use_recal))
		targets2[10]<-as.character(isolate(input$targetsadd_use_screen))
		targets2[11]<-as.character(isolate(input$targetsadd_rest_adduct))
		targets2[12]<-as.character(isolate(input$targetsadd_remark))
		targets2[13]<-as.character(isolate(input$targetsadd_tag1))
		targets2[14]<-as.character(isolate(input$targetsadd_tag2))
		targets2[15]<-as.character(isolate(input$targetsadd_tag3))
		if(isolate(input$targetsadd_date)){
			targets2[16]<-as.character(isolate(input$targetsadd_date_range[1]))
			targets2[17]<-as.character(isolate(input$targetsadd_date_range[2]))
		}
		#  targets2[18]<-  .... intercept 
		#  targets2[19]<-  .... slope	
		targets<-rbind(targets2,targets1);
		write.table(targets,file=file.path(logfile[[1]],"dataframes","targets.txt"),row.names=FALSE,sep="\t",quote=FALSE)      
		rm(targets,targets1,targets2);
		#############################################################################
		# adjust task/workflow settings #############################################
		if( 
				(logfile$workflow[2]=="yes") || 
				(logfile$workflow[6]=="yes") ||
				(logfile$workflow[12]=="TRUE") ||
				(logfile$workflow[14]=="TRUE") 
		){	# must rerun: 
				if(logfile$workflow[2]=="yes"){# must rerun: recal
					logfile$Tasks_to_redo[3]<-TRUE;
					logfile$Tasks_to_redo[3]<<-TRUE;				
					measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
					measurements[,12]<-FALSE;
					write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);				
				}
				if(logfile$workflow[3]=="yes"){# must rerun: align
					logfile$Tasks_to_redo[5]<-TRUE;
					logfile$Tasks_to_redo[5]<<-TRUE;				
					measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
					measurements[,13]<-FALSE;
					write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
				}
				# must not rerun: norm
				# must rerun: is_pattern / target_pattern
				logfile$Tasks_to_redo[8]<-TRUE;
				logfile$Tasks_to_redo[8]<<-TRUE;			
				if(logfile$workflow[9]=="yes"){# must rerun: profiled
					logfile$Tasks_to_redo[6]<-TRUE;
					logfile$Tasks_to_redo[6]<<-TRUE;				
				}
				if(logfile$workflow[10]=="yes"){# must rerun: trendblind
					logfile$Tasks_to_redo[7]<-TRUE;
					logfile$Tasks_to_redo[7]<<-TRUE;				
				}
				if(logfile$workflow[11]=="TRUE"){# must rerun: screen_IS_sam  
					logfile$Tasks_to_redo[10]<-TRUE;
					logfile$Tasks_to_redo[10]<<-TRUE;				
					measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
					measurements[,17]<-FALSE;
					write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);												
				}
				if(logfile$workflow[12]=="TRUE"){# must rerun: screen_target_sam
					logfile$Tasks_to_redo[9]<-TRUE;
					logfile$Tasks_to_redo[9]<<-TRUE;				
					measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
					measurements[,18]<-FALSE;
					write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);												
				}		
				if(logfile$workflow[13]=="TRUE"){# must rerun: screen_IS_comp
					logfile$Tasks_to_redo[12]<-TRUE;
					logfile$Tasks_to_redo[12]<<-TRUE;				
				}			
				if(logfile$workflow[14]=="TRUE"){# must rerun: screen_target_comp
					logfile$Tasks_to_redo[11]<-TRUE;
					logfile$Tasks_to_redo[11]<<-TRUE;				
				}
				if(logfile$workflow[15]=="yes"){# must rerun: profnorm - IS-based normalization
					logfile$Tasks_to_redo[15]<-TRUE;
					logfile$Tasks_to_redo[15]<<-TRUE;				
				}			
				if(logfile$workflow[16]=="yes"){# must rerun: homol
					logfile$Tasks_to_redo[14]<-TRUE;
					logfile$Tasks_to_redo[14]<<-TRUE;				
				}			
				if(logfile$workflow[17]=="yes"){# must rerun: massdef
					logfile$Tasks_to_redo[16]<-TRUE;
					logfile$Tasks_to_redo[16]<<-TRUE;				
				}			
		}
		#############################################################################			
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));      
		output$targets<<-DT::renderDataTable(read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character"));      
		output$dowhat<-renderText("Added target compound");
    }
})

# ADD TARGET LIST ############################################################
observe({
 	input$targetlist_path
	if(  (length(isolate(input$targetlist_path))) ){
		if( file.exists(as.character(isolate(input$targetlist_path[[4]]))) ){ 
			targets<-read.table(file=as.character(isolate(input$targetlist_path[[4]])),header=TRUE,sep="\t",colClasses = "character")
			write.table(targets,file=file.path(logfile[[1]],"dataframes","targets.txt"),row.names=FALSE,sep="\t",quote=FALSE)
			rm(targets)
			#############################################################################
			# adjust task/workflow settings #############################################
			if( 
				(logfile$workflow[2]=="yes") || 
				(logfile$workflow[6]=="yes") ||
				(logfile$workflow[12]=="TRUE") ||
				(logfile$workflow[14]=="TRUE") 
			){	# must rerun: 
				if(logfile$workflow[2]=="yes"){# must rerun: recal
					logfile$Tasks_to_redo[3]<-TRUE;
					logfile$Tasks_to_redo[3]<<-TRUE;				
					measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
					measurements[,12]<-FALSE;
					write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);				
				}
				if(logfile$workflow[3]=="yes"){# must rerun: align
					logfile$Tasks_to_redo[5]<-TRUE;
					logfile$Tasks_to_redo[5]<<-TRUE;				
					measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
					measurements[,13]<-FALSE;
					write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
				}
				# must not rerun: norm
				# must rerun: is_pattern / target_pattern
				logfile$Tasks_to_redo[8]<-TRUE;
				logfile$Tasks_to_redo[8]<<-TRUE;			
				if(logfile$workflow[9]=="yes"){# must rerun: profiled
					logfile$Tasks_to_redo[6]<-TRUE;
					logfile$Tasks_to_redo[6]<<-TRUE;				
				}
				if(logfile$workflow[10]=="yes"){# must rerun: trendblind
					logfile$Tasks_to_redo[7]<-TRUE;
					logfile$Tasks_to_redo[7]<<-TRUE;				
				}
				if(logfile$workflow[11]=="TRUE"){# must rerun: screen_IS_sam  
					logfile$Tasks_to_redo[10]<-TRUE;
					logfile$Tasks_to_redo[10]<<-TRUE;				
					measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
					measurements[,17]<-FALSE;
					write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);												
				}
				if(logfile$workflow[12]=="TRUE"){# must rerun: screen_target_sam
					logfile$Tasks_to_redo[9]<-TRUE;
					logfile$Tasks_to_redo[9]<<-TRUE;				
					measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
					measurements[,18]<-FALSE;
					write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);												
				}		
				if(logfile$workflow[13]=="TRUE"){# must rerun: screen_IS_comp
					logfile$Tasks_to_redo[12]<-TRUE;
					logfile$Tasks_to_redo[12]<<-TRUE;				
				}			
				if(logfile$workflow[14]=="TRUE"){# must rerun: screen_target_comp
					logfile$Tasks_to_redo[11]<-TRUE;
					logfile$Tasks_to_redo[11]<<-TRUE;				
				}
				if(logfile$workflow[15]=="yes"){# must rerun: profnorm - IS-based normalization
					logfile$Tasks_to_redo[15]<-TRUE;
					logfile$Tasks_to_redo[15]<<-TRUE;				
				}			
				if(logfile$workflow[16]=="yes"){# must rerun: homol
					logfile$Tasks_to_redo[14]<-TRUE;
					logfile$Tasks_to_redo[14]<<-TRUE;				
				}			
				if(logfile$workflow[17]=="yes"){# must rerun: massdef
					logfile$Tasks_to_redo[16]<-TRUE;
					logfile$Tasks_to_redo[16]<<-TRUE;				
				}			
			}
			#############################################################################			
			save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));   
 			output$targets<<-DT::renderDataTable(read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character"));
			output$dowhat<-renderText("Added targets list");
		}
	}
})  
##############################################################################	
 
# DELETE TARGET ##############################################################
observe({
    input$Deletetargets
    if(input$Deletetargets){
		targets<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character")
		targets<-targets[targets$ID!=as.character(isolate(input$targetsdelete_ID)),]
		write.table(targets,file=file.path(logfile[[1]],"dataframes","targets.txt"),row.names=FALSE,sep="\t",quote=FALSE)      	  
		#############################################################################
		# adjust task/workflow settings #############################################
		if( 
				(logfile$workflow[2]=="yes") || 
				(logfile$workflow[6]=="yes") ||
				(logfile$workflow[12]=="TRUE") ||
				(logfile$workflow[14]=="TRUE") 
		){	# must rerun: 
				if(logfile$workflow[2]=="yes"){# must rerun: recal
					logfile$Tasks_to_redo[3]<-TRUE;
					logfile$Tasks_to_redo[3]<<-TRUE;				
					measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
					measurements[,12]<-FALSE;
					write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);				
				}
				if(logfile$workflow[3]=="yes"){# must rerun: align
					logfile$Tasks_to_redo[5]<-TRUE;
					logfile$Tasks_to_redo[5]<<-TRUE;				
					measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
					measurements[,13]<-FALSE;
					write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
				}
				# must not rerun: norm
				# must rerun: is_pattern / target_pattern
				logfile$Tasks_to_redo[8]<-TRUE;
				logfile$Tasks_to_redo[8]<<-TRUE;			
				if(logfile$workflow[9]=="yes"){# must rerun: profiled
					logfile$Tasks_to_redo[6]<-TRUE;
					logfile$Tasks_to_redo[6]<<-TRUE;				
				}
				if(logfile$workflow[10]=="yes"){# must rerun: trendblind
					logfile$Tasks_to_redo[7]<-TRUE;
					logfile$Tasks_to_redo[7]<<-TRUE;				
				}
				if(logfile$workflow[11]=="TRUE"){# must rerun: screen_IS_sam  
					logfile$Tasks_to_redo[10]<-TRUE;
					logfile$Tasks_to_redo[10]<<-TRUE;				
					measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
					measurements[,17]<-FALSE;
					write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);												
				}
				if(logfile$workflow[12]=="TRUE"){# must rerun: screen_target_sam
					logfile$Tasks_to_redo[9]<-TRUE;
					logfile$Tasks_to_redo[9]<<-TRUE;				
					measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
					measurements[,18]<-FALSE;
					write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);												
				}		
				if(logfile$workflow[13]=="TRUE"){# must rerun: screen_IS_comp
					logfile$Tasks_to_redo[12]<-TRUE;
					logfile$Tasks_to_redo[12]<<-TRUE;				
				}			
				if(logfile$workflow[14]=="TRUE"){# must rerun: screen_target_comp
					logfile$Tasks_to_redo[11]<-TRUE;
					logfile$Tasks_to_redo[11]<<-TRUE;				
				}
				if(logfile$workflow[15]=="yes"){# must rerun: profnorm - IS-based normalization
					logfile$Tasks_to_redo[15]<-TRUE;
					logfile$Tasks_to_redo[15]<<-TRUE;				
				}			
				if(logfile$workflow[16]=="yes"){# must rerun: homol
					logfile$Tasks_to_redo[14]<-TRUE;
					logfile$Tasks_to_redo[14]<<-TRUE;				
				}			
				if(logfile$workflow[17]=="yes"){# must rerun: massdef
					logfile$Tasks_to_redo[16]<-TRUE;
					logfile$Tasks_to_redo[16]<<-TRUE;				
				}			
		}
		#############################################################################			
		output$targets<<-DT::renderDataTable(read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character"));      	  
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));      
		output$dowhat<-renderText("Deleted compound");
    }
})
##############################################################################
  
# ADD MEASUREMENT ############################################################
addmeasu<-reactive({
	input$Measadd_path
	if(  (length(isolate(input$Measadd_path))) ){
		if( file.exists(as.character(isolate(input$Measadd_path[[4]]))) ){
			if(
			  substr(as.character(isolate(input$Measadd_path[[1]])),nchar(as.character(isolate(input$Measadd_path[[1]])))-3,nchar(as.character(isolate(input$Measadd_path[[1]]))))==".RAW" ||
			  substr(as.character(isolate(input$Measadd_path[[1]])),nchar(as.character(isolate(input$Measadd_path[[1]])))-3,nchar(as.character(isolate(input$Measadd_path[[1]]))))==".raw" ||
			  substr(as.character(isolate(input$Measadd_path[[1]])),nchar(as.character(isolate(input$Measadd_path[[1]])))-3,nchar(as.character(isolate(input$Measadd_path[[1]]))))==".Raw"
			){
				if( file.exists(file.path(logfile$PW)) ){
					measurements1<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
					nameit<-names(measurements1);
					measurements1<-measurements1[measurements1[,1]!="-",]		
					if(isolate(input$Measadd_ID_autom)=="yes"){
						newID<-getID(as.numeric(measurements1[,1]))
					}else{
						newID<-as.character(isolate(input$Measadd_ID))
					}
					file.copy(
					  from=isolate(input$Measadd_path[[4]]),
					  to=file.path(logfile[[1]],"files",paste(as.character(newID),".raw",sep="")),
					  overwrite=TRUE);
					PWfile(
					  infile=file.path(logfile[[1]],"files",paste(as.character(newID),".raw",sep="")),
					  file.path(logfile[[1]],"files"),
					  as.character(isolate(input$PWpath)),
					  notintern=FALSE,
					  use_format="mzXML");     				  
					file.remove(file.path(logfile[[1]],"files",paste(as.character(newID),".raw",sep="")))
					file.remove(isolate(input$Measadd_path[[4]]));
					if(  file.exists(file.path(logfile[[1]],"files",paste(newID,".mzXML",sep=""))) || file.exists(file.path(logfile[[1]],"files",paste(newID,".mzXML",sep="")))  ){
						measurements2<-c(
							as.character(newID),
							as.character(isolate(input$Measadd_name)),
							as.character(isolate(input$Measadd_type)),
							as.character(isolate(input$Measadd_mode)),
							as.character(isolate(input$Measadd_place)),
							as.character(isolate(input$Measadd_date)),
							as.character(isolate(input$Measadd_time)),
							as.character(isolate(input$Measadd_incl)),
							"TRUE","FALSE","FALSE","FALSE","FALSE","FALSE","FALSE",
							"FALSE","FALSE","FALSE","FALSE","FALSE","FALSE"
						);		  
						measurements3<-rbind(measurements2,measurements1);
						names(measurements3)<-nameit;
						measurements3<-enviMass:::convDate(measurements3);
						write.csv(measurements3,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
						rm(measurements1,measurements2,measurements3);
						output$measurements<<-DT::renderDataTable(read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character")); 
						#############################################################################
						# adjust task/workflow settings #############################################
						doit<-as.character(isolate(input$Measadd_incl))
						doit<<-as.character(isolate(input$Measadd_incl))
						if(doit=="TRUE"){
							if(logfile$workflow[1]=="yes"){# must rerun: qc
								logfile$Tasks_to_redo[2]<-TRUE;
								logfile$Tasks_to_redo[2]<<-TRUE;				
							}
							if(logfile$workflow[2]=="yes"){# must rerun: recal
								logfile$Tasks_to_redo[3]<-TRUE;
								logfile$Tasks_to_redo[3]<<-TRUE;				
							}
							if(logfile$workflow[3]=="yes"){# must rerun: align
								logfile$Tasks_to_redo[5]<-TRUE;
								logfile$Tasks_to_redo[5]<<-TRUE;				
								#measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
								#measurements[,13]<-FALSE;
								#write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
							}
							if(logfile$workflow[4]=="yes"){# must rerun: norm
								logfile$Tasks_to_redo[4]<-TRUE;
								logfile$Tasks_to_redo[4]<<-TRUE;				
							}
							# skip is_pattern
							# skip target_pattern
							if(logfile$workflow[7]=="TRUE"){# must rerun: component isotopologues
								logfile$Tasks_to_redo[13]<-TRUE;
								logfile$Tasks_to_redo[13]<<-TRUE;		
							}
							if(logfile$workflow[8]=="TRUE"){# must rerun: component adducts
								logfile$Tasks_to_redo[13]<-TRUE;
								logfile$Tasks_to_redo[13]<<-TRUE;	
							}
							if(logfile$workflow[9]=="yes"){# must rerun: profiled
								logfile$Tasks_to_redo[6]<-TRUE;
								logfile$Tasks_to_redo[6]<<-TRUE;				
							}
							if(logfile$workflow[10]=="yes"){# must rerun: trendblind
								logfile$Tasks_to_redo[7]<-TRUE;
								logfile$Tasks_to_redo[7]<<-TRUE;				
							}
							if(logfile$workflow[11]=="TRUE"){# must rerun: screen_IS_sam  
								logfile$Tasks_to_redo[10]<-TRUE;
								logfile$Tasks_to_redo[10]<<-TRUE;				
							}
							if(logfile$workflow[12]=="TRUE"){# must rerun: screen_target_sam
								logfile$Tasks_to_redo[9]<-TRUE;
								logfile$Tasks_to_redo[9]<<-TRUE;				
							}		
							if(logfile$workflow[13]=="TRUE"){# must rerun: screen_IS_comp
								logfile$Tasks_to_redo[12]<-TRUE;
								logfile$Tasks_to_redo[12]<<-TRUE;				
							}			
							if(logfile$workflow[14]=="TRUE"){# must rerun: screen_target_comp
								logfile$Tasks_to_redo[11]<-TRUE;
								logfile$Tasks_to_redo[11]<<-TRUE;				
							}
							if(logfile$workflow[15]=="yes"){# must rerun: profnorm - IS-based normalization
								logfile$Tasks_to_redo[15]<-TRUE;
								logfile$Tasks_to_redo[15]<<-TRUE;				
							}			
							if(logfile$workflow[16]=="yes"){# must rerun: homol
								logfile$Tasks_to_redo[14]<-TRUE;
								logfile$Tasks_to_redo[14]<<-TRUE;				
							}			
							if(logfile$workflow[17]=="yes"){# must rerun: massdef
								logfile$Tasks_to_redo[16]<-TRUE;
								logfile$Tasks_to_redo[16]<<-TRUE;				
							}			
						}
						#############################################################################			
						save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));      
						output$dowhat<-renderText("Measurement added");
						cat("Measurement added\n")
						return("Measurement added\n")
					}else{
						output$dowhat<-renderText("File corrupted? - upload failed!");
						cat("File corrupted? - upload failed!")
						return("File corrupted? - upload failed!")
					}
				}else{
					output$dowhat<-renderText("Path to PW MSConvert invalid");
					cat("Path to PW MSConvert invalid")
					return("Path to PW MSConvert invalid")
				}
			}else{ #ok
				measurements1<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				nameit<-names(measurements1);
				measurements1<-measurements1[measurements1[,1]!="-",]		
				if(isolate(input$Measadd_ID_autom)=="yes"){
					newID<-getID(as.numeric(measurements1[,1]))
				}else{
					newID<-as.character(isolate(input$Measadd_ID))
				}
				file.copy(
					from=isolate(input$Measadd_path[[4]]),
					to=file.path(logfile[[1]],"files",paste(as.character(newID),".mzXML",sep="")),
					overwrite=TRUE);
				file.remove(isolate(input$Measadd_path[[4]]));
				if( (file.exists(file.path(logfile[[1]],"files",paste(newID,".mzXML",sep="")))) || (file.exists(file.path(logfile[[1]],"files",paste(newID,".mzXML",sep=""))))){
					measurements2<-c(
						as.character(newID),
						as.character(isolate(input$Measadd_name)),
						as.character(isolate(input$Measadd_type)),
						as.character(isolate(input$Measadd_mode)),
						as.character(isolate(input$Measadd_place)),
						as.character(isolate(input$Measadd_date)),
						as.character(isolate(input$Measadd_time)),
						as.character(isolate(input$Measadd_incl)),
							"TRUE","FALSE","FALSE","FALSE","FALSE","FALSE","FALSE",
							"FALSE","FALSE","FALSE","FALSE","FALSE","FALSE"
					)
					measurements3<-rbind(measurements2,measurements1);
					names(measurements3)<-nameit;
					measurements3<-enviMass:::convDate(measurements3);
					write.csv(measurements3,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
					rm(measurements1,measurements2,measurements3);
					#############################################################################
					# adjust task/workflow settings #############################################
					doit<-as.character(isolate(input$Measadd_incl))
					doit<<-as.character(isolate(input$Measadd_incl))
					if(doit=="TRUE"){
							if(logfile$workflow[1]=="yes"){# must rerun: qc
								logfile$Tasks_to_redo[2]<-TRUE;
								logfile$Tasks_to_redo[2]<<-TRUE;				
							}
							if(logfile$workflow[2]=="yes"){# must rerun: recal
								logfile$Tasks_to_redo[3]<-TRUE;
								logfile$Tasks_to_redo[3]<<-TRUE;				
							}
							if(logfile$workflow[3]=="yes"){# must rerun: align
								logfile$Tasks_to_redo[5]<-TRUE;
								logfile$Tasks_to_redo[5]<<-TRUE;				
								#measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
								#measurements[,13]<-FALSE;
								#write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
							}
							if(logfile$workflow[4]=="yes"){# must rerun: norm
								logfile$Tasks_to_redo[4]<-TRUE;
								logfile$Tasks_to_redo[4]<<-TRUE;				
							}
							# skip is_pattern
							# skip target_pattern
							if(logfile$workflow[7]=="TRUE"){# must rerun: component isotopologues
								logfile$Tasks_to_redo[13]<-TRUE;
								logfile$Tasks_to_redo[13]<<-TRUE;		
							}
							if(logfile$workflow[8]=="TRUE"){# must rerun: component adducts
								logfile$Tasks_to_redo[13]<-TRUE;
								logfile$Tasks_to_redo[13]<<-TRUE;	
							}
							if(logfile$workflow[9]=="yes"){# must rerun: profiled
								logfile$Tasks_to_redo[6]<-TRUE;
								logfile$Tasks_to_redo[6]<<-TRUE;				
							}
							if(logfile$workflow[10]=="yes"){# must rerun: trendblind
								logfile$Tasks_to_redo[7]<-TRUE;
								logfile$Tasks_to_redo[7]<<-TRUE;				
							}
							if(logfile$workflow[11]=="TRUE"){# must rerun: screen_IS_sam  
								logfile$Tasks_to_redo[10]<-TRUE;
								logfile$Tasks_to_redo[10]<<-TRUE;				
							}
							if(logfile$workflow[12]=="TRUE"){# must rerun: screen_target_sam
								logfile$Tasks_to_redo[9]<-TRUE;
								logfile$Tasks_to_redo[9]<<-TRUE;				
							}		
							if(logfile$workflow[13]=="TRUE"){# must rerun: screen_IS_comp
								logfile$Tasks_to_redo[12]<-TRUE;
								logfile$Tasks_to_redo[12]<<-TRUE;				
							}			
							if(logfile$workflow[14]=="TRUE"){# must rerun: screen_target_comp
								logfile$Tasks_to_redo[11]<-TRUE;
								logfile$Tasks_to_redo[11]<<-TRUE;				
							}
							if(logfile$workflow[15]=="yes"){# must rerun: profnorm - IS-based normalization
								logfile$Tasks_to_redo[15]<-TRUE;
								logfile$Tasks_to_redo[15]<<-TRUE;				
							}			
							if(logfile$workflow[16]=="yes"){# must rerun: homol
								logfile$Tasks_to_redo[14]<-TRUE;
								logfile$Tasks_to_redo[14]<<-TRUE;				
							}			
							if(logfile$workflow[17]=="yes"){# must rerun: massdef
								logfile$Tasks_to_redo[16]<-TRUE;
								logfile$Tasks_to_redo[16]<<-TRUE;				
							}			
					}
					#############################################################################			
					output$measurements<<-DT::renderDataTable(read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character"));
					save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));      
					output$dowhat<-renderText("Measurement added");
					cat("Measurement added\n")
					return("Measurement added\n")
				}else{
					output$dowhat<-renderText("File corrupted? - upload failed!");
					cat("File corrupted? - upload failed!")
					return("File corrupted? - upload failed!")					
				}
			}
		}else{
			output$dowhat<-renderText("File must be reloaded");
			return("File must be reloaded")
		}
    } #ok
}) #ok
output$had_meas_added<-renderText(paste(addmeasu()))  
##############################################################################
  
# DELETE MEASUREMENT #########################################################
observe({
    input$Measdel
    if(input$Measdel){
      measurements1<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
      if(any(measurements1[,1]==as.character(isolate(input$Measdel_ID)))){
		# anything left? 
		if(length(measurements1[measurements1[,1]!=as.character(isolate(input$Measdel_ID)),1])==0){
			measurements1<-data.frame(c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("FALSE"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"));
			names(measurements1)<-c("ID","Name","Type","Mode","Place","Date","Time","include","copied?","picked?",
			"checked?","recal?","align?","norm?","feat?","comp?","IS_screen?","tar_screen?","tag1","tag2","tag3")
			adjustit<-"FALSE"
        }else{
		    measurements1<-measurements1[measurements1[,1]!=as.character(isolate(input$Measdel_ID)),]
			if(any(as.character(measurements1[,8]))=="TRUE"){
				adjustit<-"TRUE"
			}else{
				adjustit<-"FALSE"
			}
		}
		write.csv(measurements1,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
        file.remove(file.path(logfile[[1]],"files",paste(as.character(isolate(input$Measdel_ID)),".mzXML",sep="")))
		#############################################################################
		# adjust task/workflow settings #############################################
		if(adjustit=="TRUE"){
							if(logfile$workflow[1]=="yes"){# must rerun: qc
								logfile$Tasks_to_redo[2]<-TRUE;
								logfile$Tasks_to_redo[2]<<-TRUE;				
							}
							if(logfile$workflow[2]=="yes"){# must rerun: recal
								logfile$Tasks_to_redo[3]<-TRUE;
								logfile$Tasks_to_redo[3]<<-TRUE;				
							}
							if(logfile$workflow[3]=="yes"){# must rerun: align
								logfile$Tasks_to_redo[5]<-TRUE;
								logfile$Tasks_to_redo[5]<<-TRUE;				
								#measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
								#measurements[,13]<-FALSE;
								#write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
							}
							if(logfile$workflow[4]=="yes"){# must rerun: norm
								logfile$Tasks_to_redo[4]<-TRUE;
								logfile$Tasks_to_redo[4]<<-TRUE;				
							}
							# skip is_pattern
							# skip target_pattern
							if(logfile$workflow[7]=="TRUE"){# must rerun: component isotopologues
								logfile$Tasks_to_redo[13]<-TRUE;
								logfile$Tasks_to_redo[13]<<-TRUE;		
							}
							if(logfile$workflow[8]=="TRUE"){# must rerun: component adducts
								logfile$Tasks_to_redo[13]<-TRUE;
								logfile$Tasks_to_redo[13]<<-TRUE;	
							}
							if(logfile$workflow[9]=="yes"){# must rerun: profiled
								logfile$Tasks_to_redo[6]<-TRUE;
								logfile$Tasks_to_redo[6]<<-TRUE;				
							}
							if(logfile$workflow[10]=="yes"){# must rerun: trendblind
								logfile$Tasks_to_redo[7]<-TRUE;
								logfile$Tasks_to_redo[7]<<-TRUE;				
							}
							if(logfile$workflow[11]=="TRUE"){# must rerun: screen_IS_sam  
								logfile$Tasks_to_redo[10]<-TRUE;
								logfile$Tasks_to_redo[10]<<-TRUE;				
							}
							if(logfile$workflow[12]=="TRUE"){# must rerun: screen_target_sam
								logfile$Tasks_to_redo[9]<-TRUE;
								logfile$Tasks_to_redo[9]<<-TRUE;				
							}		
							if(logfile$workflow[13]=="TRUE"){# must rerun: screen_IS_comp
								logfile$Tasks_to_redo[12]<-TRUE;
								logfile$Tasks_to_redo[12]<<-TRUE;				
							}			
							if(logfile$workflow[14]=="TRUE"){# must rerun: screen_target_comp
								logfile$Tasks_to_redo[11]<-TRUE;
								logfile$Tasks_to_redo[11]<<-TRUE;				
							}
							if(logfile$workflow[15]=="yes"){# must rerun: profnorm - IS-based normalization
								logfile$Tasks_to_redo[15]<-TRUE;
								logfile$Tasks_to_redo[15]<<-TRUE;				
							}			
							if(logfile$workflow[16]=="yes"){# must rerun: homol
								logfile$Tasks_to_redo[14]<-TRUE;
								logfile$Tasks_to_redo[14]<<-TRUE;				
							}			
							if(logfile$workflow[17]=="yes"){# must rerun: massdef
								logfile$Tasks_to_redo[16]<-TRUE;
								logfile$Tasks_to_redo[16]<<-TRUE;				
							}			
		}	
        #############################################################################			
        output$measurements<<-DT::renderDataTable(read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character"));
        save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));      
        output$dowhat<-renderText("Measurement deleted");
      }else{
        output$dowhat<-renderText("Invalid ID");
      }
    }
})
##############################################################################

# IMPORT MEASUREMENTS ########################################################
observe({
    input$Import_project
    if(input$Import_project){
		cat("\n Importing project files ...")
        measurements_1<<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		file_in<<-as.character(isolate(input$import_pro_dir))
		measurements_2<<-read.csv(file=file.path(file_in,"dataframes","measurements"),colClasses = "character");
		if(any(measurements_2[,1]!="-")){
			for(i in 1:length(measurements_2[,1])){
				print(as.character(i))
				if(isolate(input$Merge_project)){ # avoid duplicates?
					if(
						any(
							( measurements_1[,3]==measurements_2[i,3] )&
							( measurements_1[,4]==measurements_2[i,4] )&				
							( measurements_1[,5]==measurements_2[i,5] )&
							( measurements_1[,6]==measurements_2[i,6] )&					
							( measurements_1[,7]==measurements_2[i,7] )
						)
					){
						print("skipped a duplicate")
						next;
					}
				}
				print(as.character(i))
				if(all(measurements_1[,1]!="-")){
					newID<-getID(as.numeric(measurements_1[,1]))
				}else{
					newID<-1			
				}
				if( # mzML to mzXML conversion required?
					file.exists(file.path(file_in,"files",paste(as.character(measurements_2[i,1]),".mzML",sep=""))) &
					!file.exists(file.path(file_in,"files",paste(as.character(measurements_2[i,1]),".mzXML",sep="")))
				){ 
					PWfile(
						file.path(file_in,"files",paste(as.character(measurements_2[i,1]),".mzML",sep="")),
						file.path(file_in,"files"),
						as.character(isolate(input$PWpath)),
						notintern=FALSE,
						use_format="mzXML");    
				}
				file.copy( # copy raw data 
					  from=file.path(file_in,"files",paste(as.character(measurements_2[i,1]),".mzXML",sep="")),
					  to=file.path(logfile[[1]],"files",paste(as.character(newID),".mzXML",sep="")),
					  overwrite=TRUE);
				file.copy(
					  from=file.path(file_in,"MSlist",as.character(measurements_2[i,1])),
					  to=file.path(logfile[[1]],"MSlist",as.character(newID)),
					  overwrite=TRUE);	  
				file.copy(
					  from=file.path(file_in,"peaklist",as.character(measurements_2[i,1])),
					  to=file.path(logfile[[1]],"peaklist",as.character(newID)),
					  overwrite=TRUE);
				measurements_1<-rbind(measurements_1,measurements_2[i,])	
				at<-length(measurements_1[,1])
				measurements_1[at,1]<-newID
				measurements_1<-measurements_1[measurements_1[,1]!="-",]
			}
			write.csv(measurements_1,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
			output$measurements<<-DT::renderDataTable(read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character")); 
			rm(measurements_1,measurements_2);
			cat(" done.")		
		}else{
			cat(" no files to import - project empty?.")
		}
	}
})
##############################################################################

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  