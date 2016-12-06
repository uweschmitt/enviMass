##############################################################################
# update compounds, measurements, etc ########################################
##############################################################################

##############################################################################   
# ADD IS #####################################################################
observe({ # update selectable adducts by ionization mode 
	input$ISadd_charge
	if(isolate(input$ISadd_charge)=="positive"){
		updateSelectInput(session, "ISadd_add", "Main adduct (+):", choices = c("FALSE",as.character(adducts[adducts[,6]=="positive",1])), selected="FALSE")
		updateSelectInput(session, "IS_quant_add", "Adduct for quantification:", choices = c("FALSE",as.character(adducts[adducts[,6]=="positive",1])), selected="FALSE")
	}
	if(isolate(input$ISadd_charge)=="negative"){
		updateSelectInput(session, "ISadd_add", "Main adduct (-):", choices = c("FALSE",as.character(adducts[adducts[,6]=="negative",1])), selected="FALSE")
		updateSelectInput(session, "IS_quant_add", "Adduct for calibration & quantification:", choices = c("FALSE",as.character(adducts[adducts[,6]=="negative",1])), selected="FALSE")
	}  
})
observe({
    input$AddIS
	if(any(objects(envir=as.environment(".GlobalEnv"))=="IS")){rm(IS,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="IS")){rm(IS)}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="IS_1")){rm(IS_1,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="IS_1")){rm(IS_1)}
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
		#if(isolate(input$ISadd_date)){
		#	IS2[15]<-as.character(isolate(input$ISadd_date_range[1]))
		#	IS2[16]<-as.character(isolate(input$ISadd_date_range[2]))
		#}
		IS2[15]<-"FALSE"
		IS2[16]<-"FALSE"		
		IS2[17]<-as.character(isolate(input$Lower_intensity_bound))
		IS2[18]<-as.character(isolate(input$Upper_intensity_bound))
		IS2[19]<-as.character(isolate(input$IS_quant_add))
		IS2[20]<-as.character(isolate(input$IS_quant_peak))		
		IS2[21]<-as.character(isolate(input$IS_quant_rule))		
		IS<-rbind(IS2,IS1);
		write.table(IS,file=file.path(logfile[[1]],"dataframes","IS.txt"),row.names=FALSE,sep="\t",quote=FALSE)
		rm(IS,IS1,IS2);
		#############################################################################
		# adjust task/workflow settings #############################################
		enviMass:::workflow_set(logfile,down="pattern")		  
		#############################################################################			
		output$IS<<-DT::renderDataTable(read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character"));
		logfile[[2]][3:7]<<-rep(TRUE,length(3:7));
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));      
		output$dowhat<-renderText("Added IS compound");
		if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_obs_Add.r!")}
		enviMass:::reset_selections(session)
    }
	if(any(objects(envir=as.environment(".GlobalEnv"))=="IS")){rm(IS,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="IS")){rm(IS)}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="IS_1")){rm(IS_1,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="IS_1")){rm(IS_1)}
})
############################################################################## 

############################################################################## 
# LOAD & MODIFY IS ###########################################################
observe({
    input$LoadIS
	if(any(objects(envir=as.environment(".GlobalEnv"))=="IS")){rm(IS,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="IS")){rm(IS)}
    if(input$LoadIS){
		IS<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character")
		this<-which(IS[,"ID"]==isolate(input$ISmodif_ID))
		if(length(this)>0){
			cat("\n loading IS specifications")
			updateTextInput(session,"ISadd_ID",value=IS[this,"ID"])
			updateTextInput(session,"ISadd_name",value=IS[this,"Name"])
			updateTextInput(session,"ISadd_formula",value=IS[this,"Formula"])
			updateSelectizeInput(session,"ISadd_charge", selected = IS[this,"ion_mode" ])					
			updateTextInput(session,"ISadd_RT", value=IS[this,"RT"]) 
			if(IS[this,"RT_tolerance"]!="FALSE"){
				updateCheckboxInput(session, "ISadd_RTtol_use", value=TRUE) 
				updateTextInput(session,"ISadd_RTtol", value = IS[this,"RT_tolerance"])
			}else{
				updateCheckboxInput(session, "ISadd_RTtol_use", value=FALSE)       
				updateTextInput(session,"ISadd_RTtol", value = "FALSE")				
			}
			updateSelectizeInput(session,"ISadd_add",selected = IS[this,"main_adduct"])
			updateCheckboxInput(session, "ISadd_rest_adduct", "Restrict screening to main adduct?",value = as.logical(IS[this,"restrict_adduct"]))
			updateCheckboxInput(session, "ISadd_use_recal", "To be used for m/z recalibration?", value = as.logical(IS[this,"use_for_recalibration"]))
			updateCheckboxInput(session, "ISadd_use_screen", "To be used for screening?", value = as.logical(IS[this,"use_for_screening"]))
			updateTextInput(session,"ISadd_remark", value = IS[this,"Remark"])
			updateTextInput(session,"ISadd_tag1", value = IS[this,"tag1"])
			updateTextInput(session,"ISadd_tag2", value = IS[this,"tag2"])
			updateTextInput(session,"ISadd_tag3", value = IS[this,"tag3"])
			updateSelectizeInput(session,"IS_quant_add", selected = IS[this,"Quant_adduct"])
			updateNumericInput(session,"IS_quant_peak",  value=IS[this,"Quant_peak"])
			updateTextInput(session,"Lower_intensity_bound", value = IS[this,"Lower_intensity_bound"])	
			updateTextInput(session,"Upper_intensity_bound", value = IS[this,"Upper_intensity_bound"])
			updateSelectizeInput(session,"IS_quant_rule", selected = IS[this,"Quant_rule"])							
		}else{
			showModal(modalDialog(
				title = "Loading compound specifications failed","No internal standard with this ID available.",
				easyClose = TRUE,footer = NULL
			))			
		}
		rm(IS)
	}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="IS")){rm(IS,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="IS")){rm(IS)}
})

observe({
    input$ModifIS
	if(any(objects(envir=as.environment(".GlobalEnv"))=="IS_mod")){rm(IS_mod,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="IS_mod")){rm(IS_mod)}	
	if(any(objects(envir=as.environment(".GlobalEnv"))=="targets")){rm(targets,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="targets")){rm(targets)}	
    if(input$ModifIS){
		cat("\n Saving IS modifications")	
		new_IS<-c(
			as.character(isolate(input$ISadd_ID)),
			as.character(isolate(input$ISadd_name)),
			as.character(isolate(input$ISadd_formula)),
			as.character(isolate(input$ISadd_RT)),
			as.character(isolate(input$ISadd_RTtol)),
			as.character(isolate(input$ISadd_add)),
			as.character(isolate(input$ISadd_charge)),
			as.character(isolate(input$ISadd_use_recal)),
			as.character(isolate(input$ISadd_use_screen)),
			as.character(isolate(input$ISadd_rest_adduct)),
			as.character(isolate(input$ISadd_remark)),
			as.character(isolate(input$ISadd_tag1)),
			as.character(isolate(input$ISadd_tag2)),
			as.character(isolate(input$ISadd_tag3)),
			"FALSE",
			"FALSE",
			as.character(isolate(input$Lower_intensity_bound)),
			as.character(isolate(input$Upper_intensity_bound)),
			as.character(isolate(input$IS_quant_add)),
			as.character(isolate(input$IS_quant_peak)),
			as.character(isolate(input$IS_quant_rule))			
		)
		IS_mod<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character")
		this<-which(IS_mod[,"ID"]==isolate(input$ISadd_ID))		
		if(length(this)>0){
			IS_mod[this,]<-new_IS
		}else{ # make new entry
			IS_mod<-rbind(IS_mod,new_IS)
		}
		# check before saving ...
		targets<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character")
		say<-enviMass:::check_compounds(
			intstand_check=IS_mod,
			targets,
			isotopes,
			adducts,
			logfile
		)
		if(say=="Project consistent"){
			write.table(IS_mod,file=file.path(logfile[[1]],"dataframes","IS.txt"),row.names=FALSE,sep="\t",quote=FALSE)
			enviMass:::workflow_set(logfile,down="pattern")
			output$IS<-DT::renderDataTable(read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character"));
			output$dowhat<-renderText("Modified an IS entry");
		}else{
			showModal(modalDialog(
				title = "Modification could not be saved.",
				paste("Reason: ",say,sep=""),
				easyClose = TRUE,footer = NULL
			 ))		
		}
		######################################################################
		enviMass:::reset_selections(session)
	}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="IS_mod")){rm(IS_mod,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="IS_mod")){rm(IS_mod)}	
	if(any(objects(envir=as.environment(".GlobalEnv"))=="targets")){rm(targets,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="targets")){rm(targets)}
})
############################################################################## 

############################################################################## 
# ADD IS LIST ################################################################
observe({
 	input$ISlist_path
	if(any(objects(envir=as.environment(".GlobalEnv"))=="IS_in")){rm(IS_in,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="IS_in")){rm(IS_in)}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="targets")){rm(targets,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="targets")){rm(targets)}	
	if(  (length(isolate(input$ISlist_path))) ){
		if( file.exists(as.character(isolate(input$ISlist_path[[4]]))) ){ 
			IS_in<-try({
				read.table(file=as.character(isolate(input$ISlist_path[[4]])),header=TRUE,sep="\t",colClasses = "character")
			})	
			if(class(IS_in)=="try-error"){
				showModal(modalDialog(
						title = "Selected file could not be loaded","Make sure to provide a tab-seperated .txt file for a compound list, and retry.",
						easyClose = TRUE,footer = NULL
					  ))				
			}else{		
				targets<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character")
				say<-enviMass:::check_compounds(
					intstand_check=IS_in,
					targets,
					isotopes,
					adducts,
					logfile
				)
				rm(targets)
				if(say=="Project consistent"){
					if(isolate(input$ISlist_save_copy)){ # make copy of old IS table
						at_time<-Sys.time()
						at_time<-gsub(" ","_", at_time)
						at_time<-gsub("-","", at_time)
						at_time<-gsub(":","_", at_time)
						IS<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character")
						write.table(IS,file=file.path(logfile[[1]],"dataframes",paste("IS_",at_time,".txt",sep="")),row.names=FALSE,sep="\t",quote=FALSE)
						rm(IS)
					}
					write.table(IS_in,file=file.path(logfile[[1]],"dataframes","IS.txt"),row.names=FALSE,sep="\t",quote=FALSE)
				}else{
					showModal(modalDialog(
							title = "Compound list not consistent",
							paste("Although the selected file could be loaded, the following issue arose: ",say,". Please revise.",sep=""),
							easyClose = TRUE,footer = NULL
						  ))					
				}
				rm(IS_in)
				#############################################################################
				# adjust task/workflow settings #############################################
				enviMass:::workflow_set(logfile,down="pattern")	
				####################################################################
				save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));   
				output$IS<-DT::renderDataTable(read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character"));
				output$dowhat<-renderText("Imported IS list");
			}
			if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_obs_Add.r!")}
			enviMass:::reset_selections(session)
		}
	}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="IS_in")){rm(IS_in,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="IS_in")){rm(IS_in)}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="targets")){rm(targets,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="targets")){rm(targets)}	
})  
##############################################################################	

############################################################################## 
# SAVE IS LIST ###############################################################
observe({
	input$download_IS
	if(any(objects(envir=as.environment(".GlobalEnv"))=="IS")){rm(IS,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="IS")){rm(IS)}
	if(length(isolate(input$download_IS))>0){
		IS<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character")
		out_list<-isolate(input$download_IS)
		write.table(IS,file=as.character(parseSavePath(getVolumes()(),out_list)[,3]),row.names=FALSE,sep="\t",quote=FALSE)
	}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="IS")){rm(IS,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="IS")){rm(IS)}
}) 
##############################################################################

##############################################################################   
# DELETE IS ##################################################################
observe({
    input$DeleteIS
	if(any(objects(envir=as.environment(".GlobalEnv"))=="IS")){rm(IS,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="IS")){rm(IS)}
    if(input$DeleteIS){
		IS<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character")
		IS<-IS[IS$ID!=as.character(isolate(input$ISdelete_ID)),]
		write.table(IS,file=file.path(logfile[[1]],"dataframes","IS.txt"),row.names=FALSE,sep="\t",quote=FALSE)
		#############################################################################
		# adjust task/workflow settings #############################################
		enviMass:::workflow_set(logfile,down="pattern")	
		####################################################################
		output$IS<-DT::renderDataTable(read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character"));
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));      
		output$dowhat<-renderText("Deleted compound");
		if(any(ls()=="logfile")){stop("illegal logfile detected #1 in server_obs_Add.r!")}
		enviMass:::reset_selections(session)
    }
	if(any(objects(envir=as.environment(".GlobalEnv"))=="IS")){rm(IS,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="IS")){rm(IS)}
})
############################################################################## 
 
 
 
##############################################################################  
# ADD TARGET #################################################################
observe({ # update selectable adducts by ionization mode 
	input$targetsadd_charge
	if(isolate(input$targetsadd_charge)=="positive"){
		updateSelectInput(session, "targetsadd_add", "Main adduct (+):", choices = c("FALSE",as.character(adducts[adducts[,6]=="positive",1])), selected="FALSE")
		updateSelectInput(session, "target_quant_add", "Adduct used for calibration & quantification:", choices = c("FALSE",as.character(adducts[adducts[,6]=="positive",1])), selected="FALSE")
	}
	if(isolate(input$targetsadd_charge)=="negative"){
		updateSelectInput(session, "targetsadd_add", "Main adduct (-):", choices = c("FALSE",as.character(adducts[adducts[,6]=="negative",1])), selected="FALSE")
		updateSelectInput(session, "target_quant_add", "Adduct used for calibration & quantification:", choices = c("FALSE",as.character(adducts[adducts[,6]=="negative",1])), selected="FALSE")
	}
})
observe({
    input$Addtargets
	if(any(objects(envir=as.environment(".GlobalEnv"))=="targets1")){rm(targets1,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="targets1")){rm(targets1)}	
	if(any(objects(envir=as.environment(".GlobalEnv"))=="targets")){rm(targets,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="targets")){rm(targets)}
    if(isolate(input$Addtargets)){
		if(verbose){cat("\n in Add")}
		targets1<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character")
		targets2<-rep("FALSE",length=19)
		targets2[1]<-as.character(isolate(input$targetsadd_ID))
		targets2[2]<-as.character(isolate(input$targetsadd_name))
		targets2[3]<-as.character(isolate(input$targetsadd_formula))
		targets2[4]<-as.character(isolate(input$targetsadd_RT))
		if(isolate(input$targetsadd_RTtol_use)){
			targets2[5]<-as.character(isolate(input$targetsadd_RTtol))
		}
		targets2[6]<-as.character(isolate(input$target_quant_peak))  
		targets2[7]<-as.character(isolate(input$targetsadd_add))
		targets2[8]<-as.character(isolate(input$targetsadd_charge))
		targets2[9]<-as.character(isolate(input$targetsadd_use_recal))
		targets2[10]<-as.character(isolate(input$targetsadd_use_screen))
		targets2[11]<-as.character(isolate(input$targetsadd_rest_adduct))
		targets2[12]<-as.character(isolate(input$targetsadd_remark))
		targets2[13]<-as.character(isolate(input$targetsadd_tag1))
		targets2[14]<-as.character(isolate(input$targetsadd_tag2))
		targets2[15]<-as.character(isolate(input$targetsadd_tag3))
		#if(isolate(input$targetsadd_date)){
		#	targets2[16]<-as.character(isolate(input$targetsadd_date_range[1]))
		#	targets2[17]<-as.character(isolate(input$targetsadd_date_range[2]))
		#}
		targets2[16]<-"FALSE"
		targets2[17]<-"FALSE"
		targets2[18]<-as.character(isolate(input$warn_1)) 
		targets2[19]<-as.character(isolate(input$warn_2))	
		targets2[20]<-as.character(isolate(input$target_quant_add)) 
		targets2[21]<-as.character(isolate(input$target_quant_peak))	
		targets2[22]<-as.character(isolate(input$target_quant_rule))				
		targets<-rbind(targets2,targets1);
		write.table(targets,file=file.path(logfile[[1]],"dataframes","targets.txt"),row.names=FALSE,sep="\t",quote=FALSE)      
		rm(targets,targets1,targets2);
		#############################################################################
		# adjust task/workflow settings #############################################
		enviMass:::workflow_set(logfile,down="pattern")	
		#############################################################################			
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));      
		output$targets<-DT::renderDataTable(read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character"));      
		output$dowhat<-renderText("Added target compound");
		if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_obs_Add.r!")}
		enviMass:::reset_selections(session)
    }
	if(any(objects(envir=as.environment(".GlobalEnv"))=="targets1")){rm(targets1,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="targets1")){rm(targets1)}	
	if(any(objects(envir=as.environment(".GlobalEnv"))=="targets")){rm(targets,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="targets")){rm(targets)}
})
############################################################################## 

############################################################################## 
# LOAD & MODIFY TARGETS ######################################################
observe({
    input$Loadtarget
	if(any(objects(envir=as.environment(".GlobalEnv"))=="targets")){rm(targets,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="targets")){rm(targets)}
    if(isolate(input$Loadtarget)){
		if(verbose){cat("\n in Load")}
		targets<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character")
		this<-which(targets[,"ID"]==isolate(input$targetmodif_ID))
		if(length(this)>0){
			updateTextInput(session,"targetsadd_ID",value=targets[this,"ID"])
			updateTextInput(session,"targetsadd_name",value=targets[this,"Name"])
			updateTextInput(session,"targetsadd_formula",value=targets[this,"Formula"])
			updateSelectizeInput(session,"targetsadd_charge", selected = targets[this,"ion_mode" ])					
			updateTextInput(session,"targetsadd_RT", value=targets[this,"RT"])
			if(targets[this,"RT_tolerance"]!="FALSE"){
				updateCheckboxInput(session, "targetsadd_RTtol_use", value=TRUE) 
				updateTextInput(session,"targetsadd_RTtol", value = targets[this,"RT_tolerance"])
			}else{
				updateCheckboxInput(session, "targetsadd_RTtol_use", value=FALSE)       
				updateTextInput(session,"targetsadd_RTtol", value = "FALSE")				
			}
			updateSelectizeInput(session,"targetsadd_add",selected = targets[this,"main_adduct"])
			updateCheckboxInput(session, "targetsadd_rest_adduct", "Restrict screening to main adduct?",value = as.logical(targets[this,"restrict_adduct"]))
			updateCheckboxInput(session, "targetsadd_use_recal", "To be used for m/z recalibration?", value = as.logical(targets[this,"use_for_recalibration"]))
			updateCheckboxInput(session, "targetsadd_use_screen", "To be used for screening?", value = as.logical(targets[this,"use_for_screening"]))
			updateTextInput(session,"targetsadd_remark", value = targets[this,"Remark"])
			updateTextInput(session,"targetsadd_tag1", value = targets[this,"tag1"])
			updateTextInput(session,"targetsadd_tag2", value = targets[this,"tag2"])
			updateTextInput(session,"targetsadd_tag3", value = targets[this,"tag3"])
			updateSelectizeInput(session,"target_quant_add", selected = targets[this,"Quant_adduct"])
			updateNumericInput(session,"target_quant_peak",  value=targets[this,"Quant_peak"])
			updateTextInput(session,"target_quant_ISID",  value=targets[this,"ID_internal_standard"])		
			updateTextInput(session,"warn_1", value = targets[this,"warn_1"])	
			updateTextInput(session,"warn_2", value = targets[this,"warn_2"])
			updateSelectizeInput(session,"target_quant_rule", selected = targets[this,"Quant_rule"])							
		}else{
			showModal(modalDialog(
				title = "Loading compound specifications failed","No target with this ID available.",
				easyClose = TRUE,footer = NULL
			))			
		}
		rm(targets)
	}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="targets")){rm(targets,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="targets")){rm(targets)}
})

observe({
    input$Modiftarget
	if(any(objects(envir=as.environment(".GlobalEnv"))=="targets")){rm(targets,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="targets")){rm(targets)}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="target_mod")){rm(target_mod,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="target_mod")){rm(target_mod)}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="new_target")){rm(new_target,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="new_target")){rm(new_target)}	
	if(any(objects(envir=as.environment(".GlobalEnv"))=="IS")){rm(IS,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="IS")){rm(IS)}	
    if(isolate(input$Modiftarget)){
		if(verbose){cat("\n in Modif")}	
		new_target<-c(
			as.character(isolate(input$targetsadd_ID)),
			as.character(isolate(input$targetsadd_name)),
			as.character(isolate(input$targetsadd_formula)),
			as.character(isolate(input$targetsadd_RT)),
			as.character(isolate(input$targetsadd_RTtol)),
			as.character(isolate(input$target_quant_ISID)),
			as.character(isolate(input$targetsadd_add)),
			as.character(isolate(input$targetsadd_charge)),
			as.character(isolate(input$targetsadd_use_recal)),
			as.character(isolate(input$targetsadd_use_screen)),
			as.character(isolate(input$targetsadd_rest_adduct)),
			as.character(isolate(input$targetsadd_remark)),
			as.character(isolate(input$targetsadd_tag1)),
			as.character(isolate(input$targetsadd_tag2)),
			as.character(isolate(input$targetsadd_tag3)),
			"FALSE",
			"FALSE",
			as.character(isolate(input$warn_1)),
			as.character(isolate(input$warn_2)),
			as.character(isolate(input$target_quant_add)),
			as.character(isolate(input$target_quant_peak)),
			as.character(isolate(input$target_quant_rule))			
		)
		target_mod<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character")
		this<-which(target_mod[,"ID"]==isolate(input$targetsadd_ID))		
		if(length(this)>0){
			target_mod[this,]<-new_target
		}else{ # make new entry
			target_mod<-rbind(target_mod,new_target)
		}
		# check before saving ...
		IS<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character")
		say<-enviMass:::check_compounds(
			intstand_check=IS,
			targets=target_mod,
			isotopes,
			adducts,
			logfile
		)
		if(say=="Project consistent"){
			write.table(target_mod,file=file.path(logfile[[1]],"dataframes","targets.txt"),row.names=FALSE,sep="\t",quote=FALSE)
			enviMass:::workflow_set(logfile,down="pattern")
			output$targets<-DT::renderDataTable(read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character"));
			output$dowhat<-renderText("Modified a target entry");
		}else{
			showModal(modalDialog(
				title = "Modification could not be saved.",
				paste("Reason: ",say,sep=""),
				easyClose = TRUE,footer = NULL
			 ))		
		}
		######################################################################
		enviMass:::reset_selections(session)
	}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="targets")){rm(targets,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="targets")){rm(targets)}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="target_mod")){rm(target_mod,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="target_mod")){rm(target_mod)}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="IS")){rm(IS,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="IS")){rm(IS)}	
})
############################################################################## 

############################################################################## 
# ADD TARGET LIST ############################################################
observe({
 	input$targetlist_path
	if(any(objects(envir=as.environment(".GlobalEnv"))=="targets")){rm(targets,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="targets")){rm(targets)}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="IS_current")){rm(IS_current,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="IS_current")){rm(IS_current)}	
	if(any(objects(envir=as.environment(".GlobalEnv"))=="target_in")){rm(target_in,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="target_in")){rm(target_in)}	
	if(  (length(isolate(input$targetlist_path))) ){
		if( file.exists(as.character(isolate(input$targetlist_path[[4]]))) ){ 
			target_in<-try({
				read.table(file=as.character(isolate(input$targetlist_path[[4]])),header=TRUE,sep="\t",colClasses = "character")
			})	
			if(class(target_in)=="try-error"){
				showModal(modalDialog(
						title = "Selected file could not be loaded","Make sure to provide a tab-seperated .txt file for a target compound list, and retry.",
						easyClose = TRUE,footer = NULL
					  ))				
			}else{		
				IS_current<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character")
				say<-enviMass:::check_compounds(
					intstand_check=IS_current,
					targets=target_in,
					isotopes,
					adducts,
					logfile
				)
				rm(IS_current)
				if(say=="Project consistent"){
					if(isolate(input$targetlist_save_copy)){ # make copy of old IS table
						at_time<-Sys.time()
						at_time<-gsub(" ","_", at_time)
						at_time<-gsub("-","", at_time)
						at_time<-gsub(":","_", at_time)
						targets<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character")
						write.table(targets,file=file.path(logfile[[1]],"dataframes",paste("targets_",at_time,".txt",sep="")),row.names=FALSE,sep="\t",quote=FALSE)
						rm(targets)
					}
					write.table(target_in,file=file.path(logfile[[1]],"dataframes","targets.txt"),row.names=FALSE,sep="\t",quote=FALSE)
				}else{
					showModal(modalDialog(
							title = "Compound list not consistent",
							paste("Although the selected file could be loaded, the following issue arose: ",say,". Please revise.",sep=""),
							easyClose = TRUE,footer = NULL
						  ))					
				}
				#############################################################################
				# adjust task/workflow settings #############################################
				enviMass:::workflow_set(logfile,down="pattern")	
				####################################################################
				save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));   
				output$targets<-DT::renderDataTable(read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character"));
				output$dowhat<-renderText("Imported target compound list");
			}
			if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_obs_Add.r!")}
		}
		enviMass:::reset_selections(session)
	}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="targets")){rm(targets,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="targets")){rm(targets)}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="IS_current")){rm(IS_current,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="IS_current")){rm(IS_current)}	
	if(any(objects(envir=as.environment(".GlobalEnv"))=="target_in")){rm(target_in,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="target_in")){rm(target_in)}	
})  
##############################################################################

############################################################################## 
# SAVE TARGET LIST ###########################################################
observe({
	input$download_target
	if(any(objects(envir=as.environment(".GlobalEnv"))=="targets")){rm(targets,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="targets")){rm(targets)}
	if(length(isolate(input$download_target))>0){
		targets<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character")
		out_list<-isolate(input$download_target)
		write.table(targets,file=as.character(parseSavePath(getVolumes()(),out_list)[,3]),row.names=FALSE,sep="\t",quote=FALSE)
		rm(targets)
	}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="targets")){rm(targets,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="targets")){rm(targets)}
}) 
##############################################################################

############################################################################## 
# DELETE TARGET LIST #########################################################
observe({
    input$Deletetargets
	if(any(objects(envir=as.environment(".GlobalEnv"))=="targets")){rm(targets,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="targets")){rm(targets)}
    if(input$Deletetargets){
		targets<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character")
		targets<-targets[targets$ID!=as.character(isolate(input$targetsdelete_ID)),]
		write.table(targets,file=file.path(logfile[[1]],"dataframes","targets.txt"),row.names=FALSE,sep="\t",quote=FALSE)      	  
		#############################################################################
		# adjust task/workflow settings #############################################
		enviMass:::workflow_set(logfile,down="pattern")	
		#############################################################################			
		output$targets<<-DT::renderDataTable(read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character"));      	  
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));      
		output$dowhat<-renderText("Deleted compound");
		if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_obs_Add.r!")}
		enviMass:::reset_selections(session)
    }
	if(any(objects(envir=as.environment(".GlobalEnv"))=="targets")){rm(targets,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="targets")){rm(targets)}
})
##############################################################################
 
 
 
 
 
 
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
					measurements1<-measurements1[measurements1[,"ID"]!="-",]		
					#if(isolate(input$Measadd_ID_autom)=="yes"){
						newID<-getID(as.numeric(measurements1[,"ID"]))
					#}else{
						#newID<-as.character(isolate(input$Measadd_ID))
					#}
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
						if(isolate(input$Measadd_type)=="calibration"){ # default for calibration files: no profiling!
							use_profiling<-"FALSE"
							start_date<-as.character(isolate(input$Measadd_cal_date1))
							start_time<-as.character(isolate(input$Measadd_cal_time1))
							tag1<-as.character(isolate(input$Measadd_tag1))
							tag2<-as.character(isolate(input$Measadd_tag2))
							tag3<-as.character(isolate(input$Measadd_tag1))
						}
						if(isolate(input$Measadd_type)=="sample" || isolate(input$Measadd_type)=="blank"){
							use_profiling<-as.character(isolate(input$Measadd_profiled))
							start_date<-as.character(isolate(input$Measadd_date))
							start_time<-as.character(isolate(input$Measadd_time))		
							tag1<-as.character(isolate(input$Measadd_tag1))						
							tag2<-"FALSE"
							tag3<-as.character(isolate(input$Measadd_tag1))
						}
						if(isolate(input$Measadd_type)=="spiked"){
							use_profiling<-"FALSE"
							start_date<-as.character(isolate(input$Measadd_recov_date)) # anything
							start_time<-as.character(isolate(input$Measadd_recov_time)) # anything	
							tag1<-FALSE
							tag2<-as.character(isolate(input$Measadd_spiked_tag2))
							tag3<-FALSE
						}
						measurements2<-c(
							as.character(newID),
							as.character(isolate(input$Measadd_name)),
							as.character(isolate(input$Measadd_type)),
							as.character(isolate(input$Measadd_mode)),
							as.character(isolate(input$Measadd_place)),
							start_date,start_time,
							as.character(isolate(input$Measadd_incl)),
							"TRUE","FALSE","FALSE","FALSE","FALSE","FALSE",
							use_profiling,
							"FALSE","FALSE","FALSE",
							tag1,tag2,tag3,
							as.character(isolate(input$Measadd_cal_date2)),							
							as.character(isolate(input$Measadd_cal_time2)),
							"FALSE","FALSE","FALSE","FALSE","FALSE",
							as.character(isolate(input$Measadd_ID2))							
						)
						measurements3<-rbind(measurements2,measurements1,stringsAsFactors=FALSE);
						names(measurements3)<-nameit;
						measurements3[,"Date"]<-enviMass:::convDate(measurements3[,"Date"]);
						measurements3[,"date_end"]<-enviMass:::convDate(measurements3[,"date_end"]);
						write.csv(measurements3,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
						rm(measurements1,measurements2,measurements3);
						measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character")
						output$measurements<<-DT::renderDataTable(
							measurements[,c("ID","Name","Type","Mode","Place","Date","Time","include","profiled","tag1","tag2","tag3","date_end","time_end","ID_2")]
						); 
						#############################################################################
						# adjust task/workflow settings #############################################
						doit<-as.character(isolate(input$Measadd_incl))
						doit<<-as.character(isolate(input$Measadd_incl))
						if(doit=="TRUE"){
							if( isolate(input$Measadd_type)!="calibration" ){ # exclude calibration
								enviMass:::workflow_set(logfile,down="peakpicking",single_file=TRUE,except="calibration")	
							}
							if( isolate(input$Measadd_type)=="calibration" ){ # still, do everything
								enviMass:::workflow_set(logfile,down="peakpicking",single_file=TRUE)	
								enviMass:::workflow_set(logfile,down="calibration",single_file=TRUE)									
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
				measurements1<-measurements1[measurements1[,"ID"]!="-",]		
				#if(isolate(input$Measadd_ID_autom)=="yes"){
					newID<-getID(as.numeric(measurements1[,"ID"]))
				#}else{
				#	newID<-as.character(isolate(input$Measadd_ID))
				#}
				file.copy(
					from=isolate(input$Measadd_path[[4]]),
					to=file.path(logfile[[1]],"files",paste(as.character(newID),".mzXML",sep="")),
					overwrite=TRUE);
				file.remove(isolate(input$Measadd_path[[4]]));
				if( (file.exists(file.path(logfile[[1]],"files",paste(newID,".mzXML",sep="")))) || (file.exists(file.path(logfile[[1]],"files",paste(newID,".mzXML",sep=""))))){
					if(isolate(input$Measadd_type)=="calibration"){ # default for calibration files: no profiling!
						use_profiling<-"FALSE"
						start_date<-as.character(isolate(input$Measadd_cal_date1))
						start_time<-as.character(isolate(input$Measadd_cal_time1))
						tag1<-as.character(isolate(input$Measadd_tag1))
						tag2<-as.character(isolate(input$Measadd_tag2))
						tag3<-as.character(isolate(input$Measadd_tag1))
					}
					if(isolate(input$Measadd_type)=="sample" || isolate(input$Measadd_type)=="blank"){
						use_profiling<-as.character(isolate(input$Measadd_profiled))
						start_date<-as.character(isolate(input$Measadd_date))
						start_time<-as.character(isolate(input$Measadd_time))
						tag1<-as.character(isolate(input$Measadd_tag1))						
						tag2<-"FALSE"
						tag3<-as.character(isolate(input$Measadd_tag1))
					}
					if(isolate(input$Measadd_type)=="spiked"){
						use_profiling<-"FALSE"
						start_date<-as.character(isolate(input$Measadd_recov_date)) # anything
						start_time<-as.character(isolate(input$Measadd_recov_time)) # anything	
						tag1<-FALSE
						tag2<-as.character(isolate(input$Measadd_spiked_tag2))
						tag3<-FALSE
					}
					measurements2<-c(
						as.character(newID),
						as.character(isolate(input$Measadd_name)),
						as.character(isolate(input$Measadd_type)),
						as.character(isolate(input$Measadd_mode)),
						as.character(isolate(input$Measadd_place)),
						start_date,start_time,
						as.character(isolate(input$Measadd_incl)),
						"TRUE","FALSE","FALSE","FALSE","FALSE","FALSE",
						use_profiling,
						"FALSE","FALSE","FALSE",
						tag1,tag2,tag3,
						as.character(isolate(input$Measadd_cal_date2)),							
						as.character(isolate(input$Measadd_cal_time2)),
						"FALSE","FALSE","FALSE","FALSE","FALSE",
						as.character(isolate(input$Measadd_ID2))
					)
					measurements3<-rbind(measurements2,measurements1,stringsAsFactors=FALSE);
					names(measurements3)<-nameit;
					measurements3[,"Date"]<<-enviMass:::convDate(measurements3[,"Date"]);
					measurements3[,"date_end"]<-enviMass:::convDate(measurements3[,"date_end"]);
					write.csv(measurements3,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
					rm(measurements1,measurements2,measurements3);
					#############################################################################
					# adjust task/workflow settings #############################################
					doit<-as.character(isolate(input$Measadd_incl))
					doit<<-as.character(isolate(input$Measadd_incl))
					if(doit=="TRUE"){
						if( isolate(input$Measadd_type)!="calibration" ){ # exclude calibration
							enviMass:::workflow_set(logfile,down="peakpicking",single_file=TRUE,except="calibration")	
						}
						if( isolate(input$Measadd_type)=="calibration" ){ # still, do everything
							enviMass:::workflow_set(logfile,down="peakpicking",single_file=TRUE)	
							enviMass:::workflow_set(logfile,down="calibration",single_file=TRUE)
						}						
					}
					#############################################################################			
					measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character")
					output$measurements<<-DT::renderDataTable(
						measurements[,c("ID","Name","Type","Mode","Place","Date","Time","include","profiled","tag1","tag2","tag3","date_end","time_end","ID_2")]
					); 
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
			#########################################################################			
			# subtraction files, positive: ##########################################
			measurements3<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			if(any( (measurements3[,"ID"]!="-") & (measurements3[,"Mode"]=="positive") & (measurements3[,"Type"]!="sample"))){
				IDs_pos<-measurements3[
					(measurements3[,"Mode"]=="positive") & (measurements3[,"Type"]!="sample")
				,"ID"]
				names_pos<-measurements3[
					(measurements3[,"Mode"]=="positive") & (measurements3[,"Type"]!="sample")
				,"Name"]
				IDs_pos<-paste(IDs_pos,names_pos,sep=" - ")
				if(any(logfile[[13]]!="FALSE")){
					select_pos<-logfile[[13]]
					select_pos<-select_pos[select_pos!="FALSE"]
					# include changes from file additions / removals
					select_pos<-select_pos[!is.na(match(select_pos,IDs_pos))]
					logfile[[13]]<<-c(select_pos,"FALSE")
				}else{
					select_pos<-NULL
				}
				updateCheckboxGroupInput(session,inputId="files_pos_select_subtract", label="", choices=IDs_pos, selected = select_pos)
			}
			# subtraction files, negative: ##########################################
			if(any( (measurements3[,"ID"]!="-") & (measurements3[,"Mode"]=="negative") & (measurements3[,"Type"]!="sample"))){
				IDs_neg<-measurements3[
					(measurements3[,"Mode"]=="negative") & (measurements3[,"Type"]!="sample")
				,"ID"]
				names_neg<-measurements3[
					(measurements3[,"Mode"]=="negative") & (measurements3[,"Type"]!="sample")
				,"Name"]
				IDs_neg<-paste(IDs_neg,names_pos,sep=" - ")
				if(any(logfile[[14]]!="FALSE")){
					select_neg<-logfile[[14]]
					select_neg<-select_neg[select_neg!="FALSE"]
					select_neg<-select_neg[!is.na(match(select_neg,IDs_neg))]
					logfile[[14]]<<-c(select_neg,"FALSE")
				}else{
					select_neg<-NULL
				}
				updateCheckboxGroupInput(session,inputId="files_neg_select_subtract", label="", choices=IDs_neg, selected = select_neg)
			}
			rm(measurements3)
			#########################################################################
			save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
			#########################################################################
		}else{
			output$dowhat<-renderText("File must be reloaded");
			return("File must be reloaded")
		}
    } #ok
	if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_obs_Add.r!")}
}) #ok
output$had_meas_added<-renderText(paste(addmeasu()))  
##############################################################################
  
############################################################################## 
# DELETE MEASUREMENT #########################################################
observe({
    input$Measdel
    if(input$Measdel){
      measurements1<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
      if(any(measurements1[,"ID"]==as.character(isolate(input$Measdel_ID)))){
		# anything left? 
		if(length(measurements1[measurements1[,"ID"]!=as.character(isolate(input$Measdel_ID)),"ID"])==0){
			measurements1<-data.frame(c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),
				c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),
				c("-"),c("-"),c("-"),c("-"),c("-"),c("-"));
			names(measurements1)<-c("ID","Name","Type","Mode","Place","Date","Time","include","copied","peakpicking",
			  "checked","recal","align","norm","profiled","LOD","IS_screen","tar_screen","tag1","tag2","tag3","date_end","time_end",
			  "isotopologues","adducts","homologues","EIC_correlation","blind","ID_2")
			adjustit<-"FALSE"
        }else{
			delete_type<-measurements1[measurements1[,"ID"]==as.character(isolate(input$Measdel_ID)),"Type"]
		    measurements1<-measurements1[measurements1[,"ID"]!=as.character(isolate(input$Measdel_ID)),]
			if(any(as.character(measurements1[,"include"])=="TRUE")){
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
			if(delete_type!="calibration"){ # exclude calibration
				enviMass:::workflow_set(logfile,down="peakpicking",single_file=TRUE,except="calibration")	
			}
			if(delete_type=="calibration"){ # still, do everything
				enviMass:::workflow_set(logfile,down="peakpicking",single_file=TRUE)
				enviMass:::workflow_set(logfile,down="calibration",single_file=TRUE)
			}
		}	
		#########################################################################			
		# subtraction files, positive: ##########################################
		measurements3<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		if(any( (measurements3[,"ID"]!="-") & (measurements3[,"Mode"]=="positive") & (measurements3[,"Type"]!="sample"))){
			IDs_pos<-measurements3[
				(measurements3[,"Mode"]=="positive") & (measurements3[,"Type"]!="sample")
			,"ID"]
			names_pos<-measurements3[
				(measurements3[,"Mode"]=="positive") & (measurements3[,"Type"]!="sample")
			,"Name"]
			IDs_pos<-paste(IDs_pos,names_pos,sep=" - ")
			if(any(logfile[[13]]!="FALSE")){
				select_pos<-logfile[[13]]
				select_pos<-select_pos[select_pos!="FALSE"]
				# include changes from file additions / removals
				select_pos<-select_pos[!is.na(match(select_pos,IDs_pos))]
				logfile[[13]]<<-c(select_pos,"FALSE")
			}else{
				select_pos<-NULL
			}
			updateCheckboxGroupInput(session,inputId="files_pos_select_subtract", label="", choices=IDs_pos, selected = select_pos)
		}
		# subtraction files, negative: ##########################################
		if(any( (measurements3[,"ID"]!="-") & (measurements3[,"Mode"]=="negative") & (measurements3[,"Type"]!="sample"))){
			IDs_neg<-measurements3[
				(measurements3[,"Mode"]=="negative") & (measurements3[,"Type"]!="sample")
			,"ID"]
			names_neg<-measurements3[
				(measurements3[,"Mode"]=="negative") & (measurements3[,"Type"]!="sample")
			,"Name"]
			IDs_neg<-paste(IDs_neg,names_pos,sep=" - ")
			if(any(logfile[[14]]!="FALSE")){
				select_neg<-logfile[[14]]
				select_neg<-select_neg[select_neg!="FALSE"]
				select_neg<-select_neg[!is.na(match(select_neg,IDs_neg))]
				logfile[[14]]<<-c(select_neg,"FALSE")
			}else{
				select_neg<-NULL
			}
			updateCheckboxGroupInput(session,inputId="files_neg_select_subtract", label="", choices=IDs_neg, selected = select_neg)
		}
		rm(measurements3)
		#############################################################################
		# delete results ############################################################
		if(file.exists(file.path(logfile$project_folder,"peaklist",as.character(isolate(input$Measdel_ID))))){
			file.remove(file.path(logfile[[1]],"files",paste(as.character(isolate(input$Measdel_ID)),".mzXML",sep="")))
		}
		if( file.exists( file.path(logfile[[1]],"results","LOD",paste("plot_LOD_",isolate(input$Measdel_ID),".png",sep="") ) ) ){
			file.remove(file.path(logfile[[1]],"results","LOD",paste("plot_LOD_",isolate(input$Measdel_ID),".png",sep="") ))
		}	
		if( file.exists( file.path(logfile[[1]],"pics",paste("peakhist_",isolate(input$Measdel_ID),sep="")) ) ){
			file.remove( file.path(logfile[[1]],"pics",paste("peakhist_",isolate(input$Measdel_ID),sep=""))  )
		}		
		if( file.exists( file.path(logfile[[1]],"pics",paste("recal_",isolate(input$Measdel_ID),sep="")) ) ){
			file.remove( file.path(logfile[[1]],"pics",paste("recal_",isolate(input$Measdel_ID),sep="")) )
		}		
		if( file.exists(file.path(logfile[[1]],"pics",paste("peakmzRT_",isolate(input$Measdel_ID),sep="")) ) ){
			file.remove(file.path(logfile[[1]],"pics",paste("peakmzRT_",isolate(input$Measdel_ID),sep="")) )
		}			
		if( file.exists(file.path(logfile[[1]],"results","componentization","adducts",paste(isolate(input$Measdel_ID),sep="")) ) ){
			file.remove(file.path(logfile[[1]],"results","componentization","adducts",paste(isolate(input$Measdel_ID),sep="")) )
		}			
		if( file.exists(file.path(logfile[[1]],"results","componentization","isotopologues",paste(isolate(input$Measdel_ID),sep="")) ) ){
			file.remove(file.path(logfile[[1]],"results","componentization","isotopologues",paste(isolate(input$Measdel_ID),sep="")) )
		}			
		if( file.exists(file.path(logfile[[1]],"results","componentization","EIC_corr",paste(isolate(input$Measdel_ID),sep="")) ) ){
			file.remove(file.path(logfile[[1]],"results","componentization","EIC_corr",paste(isolate(input$Measdel_ID),sep="")) )
		}			
		#############################################################################
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
        #############################################################################			
		measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character")
		output$measurements<<-DT::renderDataTable(
			measurements[,c("ID","Name","Type","Mode","Place","Date","Time","include","profiled","tag1","tag2","tag3","date_end","time_end","ID_2")]
		); 
        save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));      
        output$dowhat<-renderText("Measurement deleted");
      }else{
        output$dowhat<-renderText("Invalid ID");
      }
    }
	if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_obs_Add.r!")}
})
##############################################################################

############################################################################## 
# IMPORT MEASUREMENTS ########################################################
impproj<-reactive({
    input$Import_project
    if(input$Import_project){
		cat("\n Importing project files ...")
        measurements_1<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		file_in<<-as.character(isolate(input$import_pro_dir))
		measurements_2<-read.csv(file=file.path(file_in,"dataframes","measurements"),colClasses = "character");
		if(any(measurements_2[,1]!="-")){
			for(i in 1:length(measurements_2[,"ID"])){
				print(as.character(i))
				if(isolate(input$Merge_project)){ # avoid duplicates?
					if(
						any(
							( measurements_1[,"Type"]==measurements_2[i,"Type"] )&
							( measurements_1[,"Mode"]==measurements_2[i,"Mode"] )&				
							( measurements_1[,"Place"]==measurements_2[i,"Place"] )&
							( measurements_1[,"Date"]==measurements_2[i,"Date"] )&					
							( measurements_1[,"Time"]==measurements_2[i,"Time"] )
						)
					){
						print("skipped a duplicate")
						next;
					}
				}
				print(as.character(i))
				if(all(measurements_1[,"ID"]!="-")){
					newID<-getID(as.numeric(measurements_1[,"ID"]))
				}else{
					newID<-1			
				}
				if( # mzML to mzXML conversion required?
					file.exists(file.path(file_in,"files",paste(as.character(measurements_2[i,"ID"]),".mzML",sep=""))) &
					!file.exists(file.path(file_in,"files",paste(as.character(measurements_2[i,"ID"]),".mzXML",sep="")))
				){ 
					PWfile(
						file.path(file_in,"files",paste(as.character(measurements_2[i,"ID"]),".mzML",sep="")),
						file.path(file_in,"files"),
						as.character(isolate(input$PWpath)),
						notintern=FALSE,
						use_format="mzXML");    
				}
				file.copy( # copy raw data 
					  from=file.path(file_in,"files",paste(as.character(measurements_2[i,"ID"]),".mzXML",sep="")),
					  to=file.path(logfile[[1]],"files",paste(as.character(newID),".mzXML",sep="")),
					  overwrite=TRUE);
				file.copy(
					  from=file.path(file_in,"MSlist",as.character(measurements_2[i,"ID"])),
					  to=file.path(logfile[[1]],"MSlist",as.character(newID)),
					  overwrite=TRUE);	  
				file.copy(
					  from=file.path(file_in,"peaklist",as.character(measurements_2[i,"ID"])),
					  to=file.path(logfile[[1]],"peaklist",as.character(newID)),
					  overwrite=TRUE);
				measurements_1<-rbind(measurements_1,measurements_2[i,],stringsAsFactors=FALSE)	
				at<-length(measurements_1[,"ID"])
				measurements_1[at,"ID"]<-newID
				measurements_1<-measurements_1[measurements_1[,"ID"]!="-",]
			}
			write.csv(measurements_1,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character")
			output$measurements<<-DT::renderDataTable(
				measurements[,c("ID","Name","Type","Mode","Place","Date","Time","include","profiled","tag1","tag2","tag3","date_end","time_end","ID_2")]
			); 
			rm(measurements_1,measurements_2);
			enviMass:::workflow_set(logfile,down="peakpicking",single_file=TRUE)			
			#########################################################################			
			# subtraction files, positive: ##########################################
			measurements3<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			if(any( (measurements3[,"ID"]!="-") & (measurements3[,"Mode"]=="positive") & (measurements3[,"Type"]!="sample"))){
				IDs_pos<-measurements3[
					(measurements3[,"Mode"]=="positive") & (measurements3[,"Type"]!="sample")
				,"ID"]
				names_pos<-measurements3[
					(measurements3[,"Mode"]=="positive") & (measurements3[,"Type"]!="sample")
				,"Name"]
				IDs_pos<-paste(IDs_pos,names_pos,sep=" - ")
				if(any(logfile[[13]]!="FALSE")){
					select_pos<-logfile[[13]]
					select_pos<-select_pos[select_pos!="FALSE"]
					# include changes from file additions / removals
					select_pos<-select_pos[!is.na(match(select_pos,IDs_pos))]
					logfile[[13]]<<-c(select_pos,"FALSE")
				}else{
					select_pos<-NULL
				}
				updateCheckboxGroupInput(session,inputId="files_pos_select_subtract", label="", choices=IDs_pos, selected = select_pos)
			}
			# subtraction files, negative: ##########################################
			if(any( (measurements3[,"ID"]!="-") & (measurements3[,"Mode"]=="negative") & (measurements3[,"Type"]!="sample"))){
				IDs_neg<-measurements3[
					(measurements3[,"Mode"]=="negative") & (measurements3[,"Type"]!="sample")
				,"ID"]
				names_neg<-measurements3[
					(measurements3[,"Mode"]=="negative") & (measurements3[,"Type"]!="sample")
				,"Name"]
				IDs_neg<-paste(IDs_neg,names_pos,sep=" - ")
				if(any(logfile[[14]]!="FALSE")){
					select_neg<-logfile[[14]]
					select_neg<-select_neg[select_neg!="FALSE"]
					select_neg<-select_neg[!is.na(match(select_neg,IDs_neg))]
					logfile[[14]]<<-c(select_neg,"FALSE")
				}else{
					select_neg<-NULL
				}
				updateCheckboxGroupInput(session,inputId="files_neg_select_subtract", label="", choices=IDs_neg, selected = select_neg)
			}
			rm(measurements3)
			#########################################################################
			logfile$summary[1,2]<<-"TRUE"
			output$summa_html<<-renderText(enviMass:::summary_html(logfile$summary));
			save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
			#########################################################################			
			output$dowhat<-renderText("Files imported.");
			cat(" done.")		
			enviMass:::reset_selections(session)
		}else{
			cat(" no files to import - project empty?.")
			output$dowhat<-renderText("Failed import: no files.");
		}
	}
	if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_obs_Add.r!")}
})
output$had_import_project<-renderText(paste(impproj()))  
##############################################################################

############################################################################## 
# MODIFY SINGLE MEASUREMENT FILES ############################################  
# LOAD
observe({
	input$Modif_load
	if(isolate(input$Modif_load)){
		measurements3<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		atID<-as.character(isolate(input$Modif_ID))
		if(any(measurements3[,"ID"]==atID)){
			updateTextInput(session, "Modif_name",value = as.character(measurements3[measurements3[,"ID"]==atID,"Name"]))
			updateSelectInput(session,"Modif_type","Type:", choices = c("sample", "blank", "calibration", "spiked"), selected = as.character(measurements3[measurements3[,"ID"]==atID,"Type"]))	
			updateSelectInput(session, "Modif_mode", selected = as.character(measurements3[measurements3[,"ID"]==atID,"Mode"]))
			updateTextInput(session, "Modif_place",value = as.character(measurements3[measurements3[,"ID"]==atID,"Place"]))
			updateDateInput(session, "Modif_date", value = as.character(measurements3[measurements3[,"ID"]==atID,"Date"]))
			updateTextInput(session, "Modif_time",value = as.character(measurements3[measurements3[,"ID"]==atID,"Time"]))
			updateTextInput(session, "Modif_tag1",value = as.character(measurements3[measurements3[,"ID"]==atID,"tag1"]))
			updateTextInput(session, "Modif_tag2",value = as.character(measurements3[measurements3[,"ID"]==atID,"tag2"]))
			updateTextInput(session, "Modif_spiked_tag2",value = as.character(measurements3[measurements3[,"ID"]==atID,"tag2"]))			
			updateTextInput(session, "Modif_tag3",value = as.character(measurements3[measurements3[,"ID"]==atID,"tag3"]))
			updateSelectInput(session, "Modif_include", selected = as.character(measurements3[measurements3[,"ID"]==atID,"include"]))
			updateSelectInput(session, "Modif_profiled", selected = as.character(measurements3[measurements3[,"ID"]==atID,"profiled"]))
			updateDateInput(session, "Modif_cal_date1", value = as.character(measurements3[measurements3[,"ID"]==atID,"Date"]))
			updateTextInput(session, "Modif_cal_time1",value = as.character(measurements3[measurements3[,"ID"]==atID,"Time"]))
			updateDateInput(session, "Modif_cal_date2", value = as.character(measurements3[measurements3[,"ID"]==atID,"date_end"]))
			updateTextInput(session, "Modif_cal_time2",value = as.character(measurements3[measurements3[,"ID"]==atID,"time_end"]))
			updateDateInput(session, "Modif_recov_date", value = as.character(measurements3[measurements3[,"ID"]==atID,"Date"]))
			updateTextInput(session, "Modif_recov_time",value = as.character(measurements3[measurements3[,"ID"]==atID,"Time"]))				
			updateTextInput(session, "Modif_ID2",value = as.character(measurements3[measurements3[,"ID"]==atID,"ID_2"]))			
			output$dowhat<-renderText("Specifications loaded into mask.");
			cat("\n specifications loaded into mask")
			rm(measurements3)
		}else{
			updateTextInput(session,"Modif_name",value = "INVALID ID")		
			updateTextInput(session,"Modif_place",value = "INVALID ID")
			updateTextInput(session,"Modif_tag3",value = "INVALID ID")			
		}
	}
})  
# EXPORT
observe({
	input$Modif_export
	if(isolate(input$Modif_export)){
		measurements3<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		atID<-as.character(isolate(input$Modif_ID))
		if(any(measurements3[,1]==atID)){
			if(isolate(input$Modif_type)=="calibration"){ # default for calibration files: no profiling!
				use_profiling<-"FALSE"
				start_date<-as.character(isolate(input$Modif_cal_date1))
				start_time<-as.character(isolate(input$Modif_cal_time1))
				tag2<-as.character(isolate(input$Modif_tag2))
			}
			if(isolate(input$Modif_type)=="sample" || isolate(input$Modif_type)=="blank"){
				use_profiling<-as.character(isolate(input$Modif_profiled))
				start_date<-as.character(isolate(input$Modif_date))
				start_time<-as.character(isolate(input$Modif_time))	
				tag2<-as.character(isolate(input$Modif_tag2))				
			}
			if(isolate(input$Modif_type)=="spiked"){ 
				use_profiling<-"FALSE"
				start_date<-as.character(isolate(input$Modif_recov_date))
				start_time<-as.character(isolate(input$Modif_recov_time))
				tag2<-as.character(isolate(input$Modif_spiked_tag2))
			}			
			measurements3[measurements3[,"ID"]==atID,"Name"]<-as.character(isolate(input$Modif_name))
			measurements3[measurements3[,"ID"]==atID,"Type"]<-as.character(isolate(input$Modif_type))
			measurements3[measurements3[,"ID"]==atID,"Mode"]<-as.character(isolate(input$Modif_mode))
			measurements3[measurements3[,"ID"]==atID,"Place"]<-as.character(isolate(input$Modif_place))
			measurements3[measurements3[,"ID"]==atID,"Date"]<-start_date
			measurements3[measurements3[,"ID"]==atID,"Date"]<-enviMass:::convDate(measurements3[measurements3[,"ID"]==atID,"Date"]);
			measurements3[measurements3[,"ID"]==atID,"Time"]<-start_time	
			measurements3[measurements3[,"ID"]==atID,"tag1"]<-as.character(isolate(input$Modif_tag1))
			measurements3[measurements3[,"ID"]==atID,"tag2"]<-tag2
			measurements3[measurements3[,"ID"]==atID,"tag3"]<-as.character(isolate(input$Modif_tag3))	
			measurements3[measurements3[,"ID"]==atID,"include"]<-as.character(isolate(input$Modif_include))				
			measurements3[measurements3[,"ID"]==atID,"profiled"]<-use_profiling	
			measurements3[measurements3[,"ID"]==atID,"date_end"]<-as.character(isolate(input$Modif_cal_date2))
			measurements3[measurements3[,"ID"]==atID,"date_end"]<-enviMass:::convDate(measurements3[measurements3[,"ID"]==atID,"date_end"]);
			measurements3[measurements3[,"ID"]==atID,"time_end"]<-as.character(isolate(input$Modif_cal_time2))
			measurements3[measurements3[,"ID"]==atID,"ID_2"]<-as.character(isolate(input$Modif_ID2))
			write.csv(measurements3,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
			output$dowhat<-renderText("Specifications saved to file table.");
			cat("\n specifications exported from mask to file table")
			rm(measurements3)
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character")
			output$measurements<<-DT::renderDataTable(
				measurements[,c("ID","Name","Type","Mode","Place","Date","Time","include","profiled","tag1","tag2","tag3","date_end","time_end","ID_2")]
			);
			######################################################################			
			# subtraction files, positive: #######################################
			measurements3<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			if(any( (measurements3[,"ID"]!="-") & (measurements3[,"Mode"]=="positive") & (measurements3[,"Type"]!="sample"))){
				IDs_pos<-measurements3[
					(measurements3[,"Mode"]=="positive") & (measurements3[,"Type"]!="sample")
				,"ID"]
				names_pos<-measurements3[
					(measurements3[,"Mode"]=="positive") & (measurements3[,"Type"]!="sample")
				,"Name"]
				IDs_pos<-paste(IDs_pos,names_pos,sep=" - ")
				if(any(logfile[[13]]!="FALSE")){
					select_pos<-logfile[[13]]
					select_pos<-select_pos[select_pos!="FALSE"]
					# include changes from file additions / removals
					select_pos<-select_pos[!is.na(match(select_pos,IDs_pos))]
					logfile[[13]]<<-c(select_pos,"FALSE")
				}else{
					select_pos<-NULL
				}
				updateCheckboxGroupInput(session,inputId="files_pos_select_subtract", label="", choices=IDs_pos, selected = select_pos)
			}
			# subtraction files, negative: #######################################
			if(any( (measurements3[,"ID"]!="-") & (measurements3[,"Mode"]=="negative") & (measurements3[,"Type"]!="sample"))){
				IDs_neg<-measurements3[
					(measurements3[,"Mode"]=="negative") & (measurements3[,"Type"]!="sample")
				,"ID"]
				names_neg<-measurements3[
					(measurements3[,"Mode"]=="negative") & (measurements3[,"Type"]!="sample")
				,"Name"]
				IDs_neg<-paste(IDs_neg,names_neg,sep=" - ")
				if(any(logfile[[14]]!="FALSE")){
					select_neg<-logfile[[14]]
					select_neg<-select_neg[select_neg!="FALSE"]
					select_neg<-select_neg[!is.na(match(select_neg,IDs_neg))]
					logfile[[14]]<<-c(select_neg,"FALSE")
				}else{
					select_neg<-NULL
				}
				updateCheckboxGroupInput(session,inputId="files_neg_select_subtract", label="", choices=IDs_neg, selected = select_neg)
			}
			rm(measurements3)
			######################################################################
			# Adjust workflow ####################################################	
			if( isolate(input$Modif_type)!="calibration" ){ # exclude calibration
				enviMass:::workflow_set(logfile,down="peakpicking",single_file=TRUE,except="calibration")	
			}
			if( isolate(input$Modif_type)=="calibration" ){ # still, do everything
				enviMass:::workflow_set(logfile,down="peakpicking",single_file=TRUE)	
				enviMass:::workflow_set(logfile,down="calibration",single_file=TRUE)
			}							
			######################################################################
			output$summa_html<<-renderText(enviMass:::summary_html(logfile$summary));
			enviMass:::reset_selections(session)
			save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
			######################################################################
		}
	}
})  
##############################################################################  

############################################################################## 
# MODIFY CALIBRATION GROUP ###################################################
# LOAD
observe({
	input$Load_cal
	if(isolate(input$Load_cal)){
		measurements3<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		measurements3<-measurements3[measurements3[,"Mode"]==isolate(input$Modif_cal_mode),,drop=FALSE]
		#measurements3<-measurements3[measurements3[,4]=="positive",,drop=FALSE]		
		measurements3<-measurements3[measurements3[,"Type"]=="calibration",,drop=FALSE]
		measurements3<-measurements3[measurements3$tag2==isolate(input$Modif_cal_group),,drop=FALSE]
		if(length(measurements3[,"ID"])>0){
			updateDateInput(session, "Modif_calgroup_date1", value = as.character(measurements3[1,"Date"]))
			updateTextInput(session, "Modif_calgroup_time1",value = as.character(measurements3[1,"Time"]))
			updateDateInput(session, "Modif_calgroup_date2", value = as.character(measurements3[1,"date_end"]))
			updateTextInput(session, "Modif_calgroup_time2",value = as.character(measurements3[1,"time_end"]))			
			cat("\nCalibration group loaded")
			concen<-paste(measurements3$tag1[order(as.numeric(measurements3$tag1),decreasing=FALSE)],sep="",collapse=", ")
			text_out<-paste("Group with ",length(measurements3[,"ID"])," files selected, containing target concentrations of ",concen,".",sep="")
			if(length(unique(measurements3[,"Date"]))!=1){
				text_out<-paste(text_out,"WARNING: group files have different start dates. Should be corrected!", sep = '<br/>')
			}
			if(length(unique(measurements3[,"Time"]))!=1){
				text_out<-paste(text_out,"WARNING: group files have different start times. Should be corrected!", sep = '<br/>')			
			}
			if(length(unique(measurements3[,"date_end"]))!=1){
				text_out<-paste(text_out,"WARNING: group files have different end dates. Should be corrected!", sep = '<br/>')			
			}
			if(length(unique(measurements3[,"time_end"]))!=1){
				text_out<-paste(text_out,"WARNING: group files have different end times. Should be corrected!", sep = '<br/>')	
			}			
			output$Modif_cal_text_load<-renderText({text_out})
		}else{
			cat("\nInvalid calibration group")		
			output$Modif_cal_text_load<-renderText({"Invalid calibration group selected. Nothing to load."})
		}
	}
})
# MODIFY
observe({
	input$Change_cal
	if(isolate(input$Change_cal)){
		measurements3<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		for_those<-which( 
			(measurements3$tag2==isolate(input$Modif_cal_group)) &
			(measurements3$Mode==isolate(input$Modif_cal_mode))
		)
		cat(for_those)
		if(length(for_those)>0){
			measurements3[for_those,"Date"]<-as.character(isolate(input$Modif_calgroup_date1))
			measurements3[for_those,"Date"]<-enviMass:::convDate(measurements3[for_those,"Date"]);
			measurements3[for_those,"Time"]<-as.character(isolate(input$Modif_calgroup_time1))	
			measurements3[for_those,"date_end"]<-as.character(isolate(input$Modif_calgroup_date2))
			measurements3[for_those,"date_end"]<-enviMass:::convDate(measurements3[for_those,"date_end"]);
			measurements3[for_those,"time_end"]<-as.character(isolate(input$Modif_calgroup_time2))		
			any_include<-any(measurements3[for_those,"include"]=="TRUE")
			write.csv(measurements3,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
			rm(measurements3)
			any_calibrated<-FALSE
			if(isolate(input$Modif_cal_mode)=="positive"){ # check occurence in positive calibration models
				if(file.exists(file.path(logfile[[1]],"quantification",paste("cal_models_pos_",isolate(input$Modif_cal_group),sep="")))){
					any_calibrated<-TRUE;
				}
			}
			if(isolate(input$Modif_cal_mode)=="negative"){ # check occurence in negative calibration models
				if(file.exists(file.path(logfile[[1]],"quantification",paste("cal_models_neg_",isolate(input$Modif_cal_group),sep="")))){
					any_calibrated<-TRUE;
				}
			}			
			if( any_include & any_calibrated ){ # included & calibration models exist? Changed time period only affects quantification, calibration models remain the same
				enviMass:::workflow_set(down="quantification",check_node=TRUE,single_file=FALSE)	
				enviMass:::workflow_set(logfile,down="calibration",single_file=TRUE)	
				enviMass:::reset_selections(session)				
			}	
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character")
			output$measurements<<-DT::renderDataTable(
				measurements[,c("ID","Name","Type","Mode","Place","Date","Time","include","profiled","tag1","tag2","tag3","date_end","time_end","ID_2")]
			); 
			output$Modif_cal_text_load<-renderText({"Modified specifications saved."})
			cat("\n Changed calibration group specifications.")
		}else{
			output$Modif_cal_text_load<-renderText({"Invalid group to modify. Done nothing."})
			cat("\n Invalid group to modify. Done nothing.")		
		}
	}
})	
# COPY
observe({
	input$Copy_cal
	if(isolate(input$Copy_cal)){
		measurements3<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		# new group name valid?
		itsok<-TRUE
		cat("\nCheck modified group name:")
		if(any(measurements3[measurements3[,"Mode"]==isolate(input$Modif_cal_mode),,drop=FALSE]$tag2==isolate(input$Copy_cal_group))){
			itsok<-FALSE
			cat(" invalid.")
			output$Modif_cal_text_load<-renderText({"New specification cannot be saved. The chosen group name already exists!"})
		}else{
			if(grepl("_",isolate(input$Copy_cal_group))){
				output$Modif_cal_text_load<-renderText({"New specification cannot be saved; no underscores permitted."})
				itsok<-FALSE
				cat(" invalid.")			
			}else{
				cat(" ok.")
			}
		}
		for_those<-which( 
			(measurements3$tag2==isolate(input$Modif_cal_group)) &
			(measurements3$Mode==isolate(input$Modif_cal_mode))
		)
		if( itsok & (length(for_those)>0) ){ # = a valid new group name		
			measurements4<-measurements3[for_those,]
			measurements4[,"Date"]<-as.character(isolate(input$Modif_calgroup_date1))
			measurements4[,"Date"]<-enviMass:::convDate(measurements4[,"Date"]);
			measurements4[,"Time"]<-as.character(isolate(input$Modif_calgroup_time1))	
			measurements4[,"date_end"]<-as.character(isolate(input$Modif_calgroup_date2))
			measurements4[,"date_end"]<-enviMass:::convDate(measurements4[,"date_end"]);
			measurements4[,"time_end"]<-as.character(isolate(input$Modif_calgroup_time2))	
			measurements4[,"LOD"]<-"FALSE"	# redo LOD!
			measurements4$tag2<-rep(isolate(input$Copy_cal_group),length(measurements4$tag2))		
			for(i in 1:length(measurements4[,"ID"])){
				oldID<-measurements4[i,"ID"]
				# get new IDs!	
				newID<-enviMass:::getID(as.numeric(c(measurements3[,"ID"],measurements4[,"ID"]))) # here, measurements4 still partly contain duplicated, old IDs
				newID<-as.character(newID)
				measurements4[i,"ID"]<-newID
				# copy mzML-files with new IDs. must exist!
				if(file.exists(file.path(logfile$project_folder,"files",paste(oldID,".mzXML",sep="")))){
					file.copy(
						from=file.path(logfile$project_folder,"files",paste(oldID,".mzXML",sep="")),
						to=file.path(logfile$project_folder,"files",paste(newID,".mzXML",sep="")),	
						overwrite=TRUE
					)
					cat("\n   copied .mzXML file.")
				}else{
					stop("Copying calibration files.mzXML: missing file. Aborted. DEBUG YOUR PROJECT!")
				}
				# copy MSlists with new IDs, if existing!				
				if(file.exists(file.path(logfile$project_folder,"peaklist",oldID))){
					file.copy(
						from=file.path(logfile$project_folder,"peaklist",oldID),
						to=file.path(logfile$project_folder,"peaklist",newID),
						overwrite=TRUE
					)
					cat("\n   copied calibration peaklist.")
				}
				# copy peaklists with new IDs, if existing!
				if(file.exists(file.path(logfile$project_folder,"MSlist",oldID))){
					file.copy(
						from=file.path(logfile$project_folder,"MSlist",oldID),
						to=file.path(logfile$project_folder,"MSlist",newID),
						overwrite=TRUE
					)
					cat("\n   copied calibration MSlist.")
				}
				if(file.exists(file.path(logfile[[1]],"pics",paste("recal_",oldID,sep="")))){
					file.copy(				
						from=file.path(logfile[[1]],"pics",paste("recal_",oldID,sep="")),
						to=file.path(logfile[[1]],"pics",paste("recal_",newID,sep="")),
						overwrite=TRUE
					)
				}
				if(file.exists(file.path(logfile[[1]],"pics",paste("peakhist_",oldID,sep="")))){
					file.copy(				
						from=file.path(logfile[[1]],"pics",paste("peakhist_",oldID,sep="")),
						to=file.path(logfile[[1]],"pics",paste("peakhist_",newID,sep="")),
						overwrite=TRUE
					)
				}
				if(file.exists(file.path(logfile[[1]],"pics",paste("peakmzRT_",oldID,sep="")))){
					file.copy(				
						from=file.path(logfile[[1]],"pics",paste("peakmzRT_",oldID,sep="")),
						to=file.path(logfile[[1]],"pics",paste("peakmzRT_",newID,sep="")),
						overwrite=TRUE
					)
				}
				if(file.exists(file.path(logfile[[1]],"results","LOD",paste("plot_LOD_",oldID,".png",sep="")))){
					file.copy(				
						from=file.path(logfile[[1]],"results","LOD",paste("plot_LOD_",oldID,".png",sep="")),
						to=file.path(logfile[[1]],"results","LOD",paste("plot_LOD_",newID,".png",sep="")),
						overwrite=TRUE
					)
				}				
			}			
			measurements3<-rbind(measurements3,measurements4)
			write.csv(measurements3,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
			rm(measurements3,measurements4)
			# copy existing calibration models
			if(isolate(input$Modif_cal_mode)=="positive"){ # positive
				if(file.exists(file.path(logfile[[1]],"quantification",paste("cal_models_pos_",isolate(input$Modif_cal_group),sep="")))){
					source(file=file.path(logfile[[1]],"quantification",paste("cal_models_pos_",isolate(input$Modif_cal_group),sep="")),local=as.environment(".GlobalEnv"));					
					names(cal_models_pos)<<-isolate(input$Copy_cal_group) # rename!
					dump("cal_models_pos",file=file.path(logfile[[1]],"quantification",paste("cal_models_pos_",isolate(input$Copy_cal_group),sep="")),envir=as.environment(".GlobalEnv"));					
					rm(cal_models_pos,envir=as.environment(".GlobalEnv"))
				}
			}
			if(isolate(input$Modif_cal_mode)=="negative"){ # negative
				if(file.exists(file.path(logfile[[1]],"quantification",paste("cal_models_neg_",isolate(input$Modif_cal_group),sep="")))){		
					source(file=file.path(logfile[[1]],"quantification",paste("cal_models_neg_",isolate(input$Modif_cal_group),sep="")),local=as.environment(".GlobalEnv"));
					names(cal_models_neg)<<-isolate(input$Copy_cal_group) # rename!
					dump("cal_models_neg",file=file.path(logfile[[1]],"quantification",paste("cal_models_neg_",isolate(input$Copy_cal_group),sep="")),envir=as.environment(".GlobalEnv"));					
					rm(cal_models_neg,envir=as.environment(".GlobalEnv"))
				}
			}
			enviMass:::workflow_set(down="LOD",check_node=TRUE,single_file=TRUE)
			enviMass:::workflow_set(logfile,down="calibration",single_file=TRUE)
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character")
			output$measurements<<-DT::renderDataTable(
				measurements[,c("ID","Name","Type","Mode","Place","Date","Time","include","profiled","tag1","tag2","tag3","date_end","time_end","ID_2")]
			);
			cat("Calibration file set copied.")
			output$Modif_cal_text_load<-renderText({"Calibration file set copied."})
			enviMass:::reset_selections(session)
		}else{
			cat("Calibration file set not copied.")
			output$Modif_cal_text_load<-renderText({"Calibration file set not copied, such a group already exists and cannot be overwritten."})		
		}
	}
})	
# DELETE
observe({
	input$yes_delete_cal
	if(isolate(input$yes_delete_cal)){
		toggleModal(session,"Del_cal_confirm", toggle = "close")
		measurements3<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		for_those<-which( 
			(measurements3$tag2==isolate(input$Modif_cal_group)) &
			(measurements3$Mode==isolate(input$Modif_cal_mode))
		)
		if(length(for_those)>0){
			cat(for_those)
			any_include<-any(measurements3[for_those,"include"]=="TRUE")
			rem_IDs<-measurements3[for_those,"ID"]
			measurements3<-measurements3[-for_those,]
			write.csv(measurements3,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
			rm(measurements3)
			for(i in 1:length(rem_IDs)){ # remove MSlists & peaklist
				if(file.exists(file.path(logfile$project_folder,"files",paste(rem_IDs[i],".mzXML",sep="")))){
					file.remove(file.path(logfile$project_folder,"files",paste(rem_IDs[i],".mzXML",sep="")))
				}else{
					cat("\nMissing file.mzXML while deleting calibration files. Debug this project?")
				}
				if(file.exists(file.path(logfile$project_folder,"peaklist",rem_IDs[i]))){
					file.remove(file.path(logfile$project_folder,"peaklist",rem_IDs[i]))
				}
				if(file.exists(file.path(logfile$project_folder,"MSlist",rem_IDs[i]))){
					file.remove(file.path(logfile$project_folder,"MSlist",rem_IDs[i]))
				}	
				if(file.exists(file.path(logfile[[1]],"pics",paste("recal_",rem_IDs[i],sep="")))){
					file.remove(file.path(logfile[[1]],"pics",paste("recal_",rem_IDs[i],sep="")))
				}
				if(file.exists(file.path(logfile[[1]],"pics",paste("peakhist_",rem_IDs[i],sep="")))){
					file.remove(file.path(logfile[[1]],"pics",paste("peakhist_",rem_IDs[i],sep="")))
				}
				if(file.exists(file.path(logfile[[1]],"pics",paste("peakmzRT_",rem_IDs[i],sep="")))){
					file.remove(file.path(logfile[[1]],"pics",paste("peakmzRT_",rem_IDs[i],sep="")))
				}
				if(file.exists(file.path(logfile[[1]],"results","LOD",paste("plot_LOD_",rem_IDs[i],".png",sep="")))){
					file.remove(file.path(logfile[[1]],"results","LOD",paste("plot_LOD_",rem_IDs[i],".png",sep="")))
				}				
				if(file.exists(file.path(logfile[[1]],"results","componentization","adducts",rem_IDs[i]))){
					file.remove(file.path(logfile[[1]],"results","componentization","adducts",rem_IDs[i]))
				}				
				if(file.exists(file.path(logfile[[1]],"results","componentization","EIC_corr",rem_IDs[i]))){
					file.remove(file.path(logfile[[1]],"results","componentization","EIC_corr",rem_IDs[i]))
				}
				if(file.exists(file.path(logfile[[1]],"results","componentization","isotopologues",rem_IDs[i]))){
					file.remove(file.path(logfile[[1]],"results","componentization","isotopologues",rem_IDs[i]))
				}			
				cat("\n file removed")
			}
			if(isolate(input$Modif_cal_mode)=="positive"){ # positive
				if(file.exists(file.path(logfile[[1]],"quantification",paste("cal_models_pos_",isolate(input$Modif_cal_group),sep="")))){
					file.remove(file.path(logfile[[1]],"quantification",paste("cal_models_pos_",isolate(input$Modif_cal_group),sep="")))					
				}
			}
			if(isolate(input$Modif_cal_mode)=="negative"){ # negative
				if(file.exists(file.path(logfile[[1]],"quantification",paste("cal_models_neg_",isolate(input$Modif_cal_group),sep="")))){
					file.remove(file.path(logfile[[1]],"quantification",paste("cal_models_neg_",isolate(input$Modif_cal_group),sep="")))					
				}
			}
			if( any_include ){ # included & calibration models exist?
				enviMass:::workflow_set(down="calibration",check_node=TRUE,single_file=FALSE)	# is this optional? after all, the sets are just removed ...
				enviMass:::workflow_set(down="quantification",check_node=TRUE,single_file=FALSE)	
				updateSelectInput(session,"Ion_mode_Cal",selected = "none") # stops, in combination with Tasks_to_redo, invalid selections in the calibration tab!
			}	
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character")
			output$measurements<<-DT::renderDataTable(
				measurements[,c("ID","Name","Type","Mode","Place","Date","Time","include","profiled","tag1","tag2","tag3","date_end","time_end","ID_2")]
			); 
			output$Modif_cal_text_load<-renderText({"Calibration group deleted."})
			enviMass:::reset_selections(session)
			cat("\n Calibration group deleted.")
		}else{
			output$Modif_cal_text_load<-renderText({"Invalid group to delete. Done nothing."})
			cat("\n Invalid group to delete. Done nothing.")		
		}
	}
})	
##############################################################################  

############################################################################## 
# BATCH IMPORT FOLDER ########################################################  
impfolder<-reactive({
	input$Import_file_folder
	if(isolate(input$Import_file_folder)){
		file_in<<-as.character(isolate(input$import_file_folder))
		getfiles<<-list.files(path=file_in)
		if(length(getfiles)>0){
			cat("\nStarting upload ...");
			many<-0;					
			for(i in 1:length(getfiles)){
				filepath<-file.path(file_in,getfiles[i])
				file_ending<-enviMass:::filetype(getfiles[i],check=TRUE)
				if(
					file.exists(filepath) & file_ending # in case of modifications meanwhile
				){
					cat(paste("\n   processing file # ",i,sep=""));
					file_ending<-enviMass:::filetype(getfiles[i],check=FALSE)
					measurements1<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
					nameit<-names(measurements1);
					measurements1<-measurements1[measurements1[,1]!="-",,drop=FALSE]		
					if(!isolate(input$Import_file_folder_overwrite)){ # skip if file of same name exists in project?
						if(any(measurements1[,"Name"]==getfiles[i])){
							cat(" - skipped.");
							output$dowhat<-renderText("File import - file skipped.");
							next;
						}
					}
					# define minimum available date
					if(length(measurements1[,"ID"])==0){
						at_date<<-as.character(isolate(input$Measadd_date))
					}else{
						all_dates<-measurements1[,"Date"]
						at_date<<-enviMass:::minDate(all_dates,get_min=FALSE)
						at_date<<-enviMass:::incrDate(Date=at_date,increment=1);#cat(paste("\n",at_date))	
					}
					newID<-as.character(getID(as.numeric(measurements1[,"ID"])))
					if(file_ending==".mzXML"){			
						file.copy(
							from=filepath,
							to=file.path(logfile[[1]],"files",paste(newID,".mzXML",sep="")),
							overwrite=TRUE)	
						if( file.exists(file.path(logfile[[1]],"files",paste(newID,".mzXML",sep=""))) ){ # check: copy completed?			
							mz1<-readMzXmlData:::readMzXmlFile(
								mzXmlFile=file.path(logfile[[1]],"files",paste(newID,".mzXML",sep="")),
								removeMetaData = FALSE,verbose = FALSE)
							ioniz<-mz1[[1]]$metaData$polarity											
							with_mode<-"unknown"
							if(ioniz=="+"){
								with_mode<-"positive"
							}
							if(ioniz=="-"){
								with_mode<-"negative"
							}							
							measurements2<-c(
								newID,
								getfiles[i], # name
								"sample",
								with_mode,
								"FALSE", # place
								at_date, # incremented in order of upload
								as.character("12:00:00"),
								"TRUE", # to be included?
								"TRUE","FALSE","FALSE","FALSE","FALSE","FALSE",
								"TRUE", # to be profiled?
								"FALSE","FALSE","FALSE",
								"FALSE","FALSE","FALSE",
								at_date,
								as.character("12:00:00"),
								"FALSE","FALSE","FALSE","FALSE","FALSE","FALSE"								
							)								
							measurements3<-rbind(measurements2,measurements1,stringsAsFactors=FALSE);
							names(measurements3)<-nameit;
							measurements3[,"Date"]<-enviMass:::convDate(measurements3[,"Date"]);
							measurements3[,"date_end"]<-enviMass:::convDate(measurements3[,"date_end"]);
							write.csv(measurements3,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
							rm(measurements1,measurements2,measurements3);
							#############################################################################			
							measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character")
							output$measurements<<-DT::renderDataTable(
								measurements[,c("ID","Name","Type","Mode","Place","Date","Time","include","profiled","tag1","tag2","tag3","date_end","time_end","ID_2")]
							);
							save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));      
							output$dowhat<-renderText("Files copied");
							cat(" - file copied")
							#############################################################################
							many<-(many+1)
						}else{
							output$dowhat<-renderText("File corrupted? - copying from folder failed!");
							cat("File corrupted? - copying from folder failed!")
							return("File corrupted? - copying from folder failed!")					
						}
					}else{ # ... then its a .raw file		
						if( file.exists(file.path(logfile$PW)) ){
							file.copy(
								from=filepath,
								to=file.path(logfile[[1]],"files",paste(newID,".raw",sep="")),
								overwrite=TRUE)			
							PWfile(
								infile=file.path(logfile[[1]],"files",paste(newID,".raw",sep="")),
								file.path(logfile[[1]],"files"),
								as.character(isolate(input$PWpath)),
								notintern=FALSE,
								use_format="mzXML");     				  
							file.remove(file.path(logfile[[1]],"files",paste(newID,".raw",sep="")))
							if( file.exists(file.path(logfile[[1]],"files",paste(newID,".mzXML",sep=""))) ){ # copy completed and conversion ok?			
								mz1<-readMzXmlData:::readMzXmlFile(
									mzXmlFile=file.path(logfile[[1]],"files",paste(newID,".mzXML",sep="")),
									removeMetaData = FALSE,verbose = FALSE)
								ioniz<-mz1[[1]]$metaData$polarity											
								with_mode<-"unknown"
								if(ioniz=="+"){
									with_mode<-"positive"
								}
								if(ioniz=="-"){
									with_mode<-"negative"
								}							
								measurements2<-c(
									newID,
									getfiles[i], # name
									"sample",
									with_mode,
									"FALSE", # place
									at_date, # incremented in order of upload
									as.character("12:00:00"),
									"TRUE", # to be included?
									"TRUE","FALSE","FALSE","FALSE","FALSE","FALSE",
									"TRUE", # to be profiled?
									"FALSE","FALSE","FALSE",
									"FALSE","FALSE","FALSE",
									at_date,
									as.character("12:00:00"),
									"FALSE","FALSE","FALSE","FALSE","FALSE","FALSE"								
								)	
								measurements3<-rbind(measurements2,measurements1,stringsAsFactors=FALSE);
								names(measurements3)<-nameit;
								measurements3[,"Date"]<-enviMass:::convDate(measurements3[,"Date"]);
								measurements3[,"date_end"]<-enviMass:::convDate(measurements3[,"date_end"]);
								write.csv(measurements3,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
								rm(measurements1,measurements2,measurements3);
								#############################################################################			
								measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character")
								output$measurements<<-DT::renderDataTable(
									measurements[,c("ID","Name","Type","Mode","Place","Date","Time","include","profiled","tag1","tag2","tag3","date_end","time_end","ID_2")]
								);
								save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));      
								output$dowhat<-renderText("Files copied");
								cat(" - file copied")
								#############################################################################
								many<-(many+1)
							}else{
								output$dowhat<-renderText("File corrupted? - copying from folder failed!");
								cat("File corrupted? - copying from folder failed!")
								return("File corrupted? - copying from folder failed!")					
							}
						}else{
							output$dowhat<-renderText(".RAW-file: path to PW MSConvert.exe invalid. Please correct in tabs Settings/General!");
							cat(".RAW-file: path to PW MSConvert.exe invalid. Please correct in tabs Settings/General!")
							return(".RAW-file: path to PW MSConvert.exe invalid. Please correct in tabs Settings/General!")
						}	
					}
				}
			}
			if(many>0){
				enviMass:::workflow_set(logfile,down="peakpicking",single_file=TRUE) 
				logfile$summary[1,2]<<-"TRUE"
				output$summa_html<<-renderText(enviMass:::summary_html(logfile$summary));
				save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
				output$dowhat<-renderText(paste(many,"files imported"))
				cat(paste("\n",many,"files imported"))
				enviMass:::reset_selections(session)
				return(paste(many,"files imported"))
			}else{
				cat("\nNo files imported")
				return(paste("No files imported"))			
			}
			############################ BAUSTELLE	
			# Default: files all imported as "sample" - no update of blind subtraction files thus implemented
			############################ BAUSTELLE	
		}else{
			cat("\nNothing to import: no files or wrong folder.")
			return("Nothing to import: no files or wrong folder.")
		}
	}
	if(any(ls()=="logfile")){stop("\n illegal logfile detected #2 in server_obs_Add.r!")}
}) #ok
output$had_import_folder<-renderText(paste(impfolder()))
##############################################################################  
 
############################################################################## 
# IMPORT PARAMETERS ##########################################################  
observe({
    input$Import_project_para
    if(input$Import_project_para){
		cat("\n Importing project parameters ...")
		logfile_here<<-logfile;
		file_in<-as.character(isolate(input$import_pro_dir_paras))
		load(file.path(file_in,"logfile.emp"),envir=as.environment(".GlobalEnv"))
		logfile_other<<-logfile;
		rm(logfile,envir=as.environment(".GlobalEnv"))
		logfile<<-logfile_here
		if(logfile$version==logfile_other$version){
			logfile[[4]]<<-logfile_other[[4]]
			logfile$parameters<<-logfile_other$parameters
			logfile$adducts_pos<<-logfile_other$adducts_pos		
			logfile$adducts_neg<<-logfile_other$adducts_neg		 
			logfile$isotopes<<-logfile_other$isotopes		
			rm(logfile_other,logfile_here,envir=as.environment(".GlobalEnv"))
			save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp")); 
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			source("server_variables_in.R", local=TRUE)
			output$dowhat<<-renderText("Parameters imported.");
			enviMass:::reset_selections(session)
			cat(" done. \n")
		}else{
			output$dowhat<<-renderText("Parameters failed: incompatible enviMass versions.");
			cat(" failed. \n")		
		}
	}
	if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_obs_Add.r!")}
}) 
############################################################################## 
 
##############################################################################
# PLOT FILE OVERVIEW #########################################################
ranges_overview <- reactiveValues(x = NULL, y = NULL)

observe({
	input$Measadd_path
	input$Measdel
	input$Import_project
	input$Modif_export
	input$Change_cal
	input$Copy_cal
	input$yes_delete_cal
	output$file_overview <- renderPlot({
		enviMass:::plot_measurements(logfile,ranges_overview)
	})
})

observeEvent(input$file_overview_dblclick, {
    brush <- input$file_overview_brush
    if (!is.null(brush)) {
		ranges_overview$x <- c(brush$xmin, brush$xmax)
		ranges_overview$y <- NULL#c(brush$ymin, brush$ymax)
    } else {
		ranges_overview$x <- NULL
		ranges_overview$y <- NULL
    }
})

observeEvent(input$file_overview_brush, {
		brush <- input$file_overview_brush
		if (!is.null(brush)) {
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");	
			dated<-measurements[,"Date"]
			timed<-measurements[,"Time"]
			datetime<-c()
			for(i in 1:length(timed)){
				datetime<-c(datetime,paste(dated[i],timed[i],"CET",sep=" "))
			}
			atPOSIX<-as.POSIXct(datetime);	
			#############################################
			# positive, samples			
			these<-measurements[
				(measurements[,"Mode"]=="positive") &
				(measurements[,"Type"]=="sample") &
				(atPOSIX>=as.POSIXct(brush$xmin, origin = "1970-01-01")) &
				(atPOSIX<=as.POSIXct(brush$xmax, origin = "1970-01-01"))
			,"ID"]
			output$info_files_pos_samp <- renderText({
				if(length(these)>0){paste("<font color=\"black\"> Sample IDs: ",paste(these,collapse=", ")," </font>",sep="")}else{paste("<font color=\"black\"> No sample files selected </font>",sep="")}
			})
			# positive, blind	
			these2<-measurements[
				(measurements[,"Mode"]=="positive") &
				(measurements[,"Type"]=="blank") &
				(atPOSIX>=as.POSIXct(brush$xmin, origin = "1970-01-01")) &
				(atPOSIX<=as.POSIXct(brush$xmax, origin = "1970-01-01"))
			,"ID"]
			output$info_files_pos_blind<- renderText({
				if(length(these2)>0){paste("<font color=\"green\"> Blanks/blind IDs: ",paste(these2,collapse=", ")," </font>",sep="")}else{paste("<font color=\"green\"> No blind files selected </font>",sep="")}
			})
			# positive, calibration	
			if(any((measurements[,"Type"]=="calibration") & (measurements[,"Mode"]=="positive"))){
				these3<-which((measurements[,"Mode"]=="positive") &(measurements[,"Type"]=="calibration") )
				date_start<-measurements[these3,"Date"]
				date_end<-measurements[these3,"date_end"]
				time_start<-measurements[these3,"Time"]
				time_end<-measurements[these3,"time_end"]
				datetime_start<-c()
				datetime_end<-c()
				for(i in 1:length(time_start)){
					datetime_start<-c(datetime_start,paste(date_start[i],time_start[i],"CET",sep=" "))
					datetime_end<-c(datetime_end,paste(date_end[i],time_end[i],"CET",sep=" "))
				}
				atPOSIX_start<-as.POSIXct(datetime_start);	
				atPOSIX_end<-as.POSIXct(datetime_end);	
				these3<-these3[
					(atPOSIX_start<=as.POSIXct(brush$xmax, origin = "1970-01-01")) &
					(atPOSIX_end>=as.POSIXct(brush$xmin, origin = "1970-01-01"))
				]
				these4<-unique(measurements[these3,"tag2"])
				these3<-measurements[these3,"ID"]
			}else{
				these3<-c()
				these4<-c()
			}
			output$info_files_pos_cal<- renderText({
				if(length(these3)>0){
					paste("<font color=\"red\"> Calibration file IDs: ",paste(these3,collapse=", ")," </font>",sep="")
				}else{
					paste("<font color=\"red\"> No calibration files selected </font>",sep="")
				}
			})
			output$info_files_pos_calgroup<- renderText({
				if(length(these4)>0){
					paste("<font color=\"red\"> Calibration file groups: ",paste(these4,collapse=", ")," </font>",sep="")
				}else{
					paste("<font color=\"red\"> No calibration file groups selected </font>",sep="")
				}
			})		
			# positive, spiked	
			these5<-measurements[
				(measurements[,"Mode"]=="positive") &
				(measurements[,"Type"]=="spiked") &
				(atPOSIX>=as.POSIXct(brush$xmin, origin = "1970-01-01")) &
				(atPOSIX<=as.POSIXct(brush$xmax, origin = "1970-01-01"))
			,"ID"]
			output$info_files_pos_spiked<- renderText({
				if(length(these5)>0){paste("<font color=\"blue\"> Spiked file IDs: ",paste(these5,collapse=", ")," </font>",sep="")}else{paste("<font color=\"blue\"> No spiked files selected </font>",sep="")}
			})			
			#############################################
			# negative, samples			
			these6<-measurements[
				(measurements[,"Mode"]=="negative") &
				(measurements[,"Type"]=="sample") &
				(atPOSIX>=as.POSIXct(brush$xmin, origin = "1970-01-01")) &
				(atPOSIX<=as.POSIXct(brush$xmax, origin = "1970-01-01"))
			,"ID"]
			output$info_files_neg_samp <- renderText({
				if(length(these6)>0){paste("<font color=\"black\"> Sample IDs: ",paste(these6,collapse=", ")," </font>",sep="")}else{paste("<font color=\"black\"> No sample files selected </font>",sep="")}
			})
			# negative, blind	
			these7<-measurements[
				(measurements[,"Mode"]=="negative") &
				(measurements[,"Type"]=="blank") &
				(atPOSIX>=as.POSIXct(brush$xmin, origin = "1970-01-01")) &
				(atPOSIX<=as.POSIXct(brush$xmax, origin = "1970-01-01"))
			,"ID"]
			output$info_files_neg_blind<- renderText({
				if(length(these7)>0){paste("<font color=\"green\"> Blanks/blind IDs: ",paste(these7,collapse=", ")," </font>",sep="")}else{paste("<font color=\"green\"> No blind files selected </font>",sep="")}
			})
			# negative, calibration	
			if(any((measurements[,"Type"]=="calibration") & (measurements[,"Mode"]=="negative"))){
				these8<-which((measurements[,"Mode"]=="negative") &(measurements[,"Type"]=="calibration") )
				date_start<-measurements[these8,"Date"]
				date_end<-measurements[these8,"date_end"]
				time_start<-measurements[these8,"Time"]
				time_end<-measurements[these8,"time_end"]
				datetime_start<-c()
				datetime_end<-c()
				for(i in 1:length(time_start)){
					datetime_start<-c(datetime_start,paste(date_start[i],time_start[i],"CET",sep=" "))
					datetime_end<-c(datetime_end,paste(date_end[i],time_end[i],"CET",sep=" "))
				}
				atPOSIX_start<-as.POSIXct(datetime_start);	
				atPOSIX_end<-as.POSIXct(datetime_end);	
				these8<-these8[
					(atPOSIX_start<=as.POSIXct(brush$xmax, origin = "1970-01-01")) &
					(atPOSIX_end>=as.POSIXct(brush$xmin, origin = "1970-01-01"))
				]
				these9<-unique(measurements[these8,"tag2"])
				these8<-measurements[these8,"ID"]
			}else{
				these8<-c()
				these9<-c()
			}
			output$info_files_neg_cal<- renderText({
				if(length(these8)>0){
					paste("<font color=\"red\"> Calibration file IDs: ",paste(these8,collapse=", ")," </font>",sep="")
				}else{
					paste("<font color=\"red\"> No calibration files selected </font>",sep="")
				}
			})
			output$info_files_neg_calgroup<- renderText({
				if(length(these9)>0){
					paste("<font color=\"red\"> Calibration file groups: ",paste(these9,collapse=", ")," </font>",sep="")
				}else{
					paste("<font color=\"red\"> No calibration file groups selected </font>",sep="")
				}
			})		
			# negative, spiked	
			these10<-measurements[
				(measurements[,"Mode"]=="negative") &
				(measurements[,"Type"]=="spiked") &
				(atPOSIX>=as.POSIXct(brush$xmin, origin = "1970-01-01")) &
				(atPOSIX<=as.POSIXct(brush$xmax, origin = "1970-01-01"))
			,"ID"]
			output$info_files_neg_spiked<- renderText({
				if(length(these10)>0){paste("<font color=\"blue\"> Spiked file IDs: ",paste(these10,collapse=", ")," </font>",sep="")}else{paste("<font color=\"blue\"> No spiked files selected </font>",sep="")}
			})			
			#############################################		
			rm(measurements)
		}	
	})
	
##############################################################################
  
  
  
  
  
  
