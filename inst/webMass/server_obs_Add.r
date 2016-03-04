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
		IS2[17]<-as.character(isolate(input$Lower_intensity_bound))
		IS2[18]<-as.character(isolate(input$Upper_intensity_bound))
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
			enviMass:::workflow_set(logfile,down="pattern")	
			####################################################################
			save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));   
 			output$IS<<-DT::renderDataTable(read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character"));
			output$dowhat<-renderText("Added IS list");
			if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_obs_Add.r!")}
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
		enviMass:::workflow_set(logfile,down="pattern")	
		####################################################################
		output$IS<<-DT::renderDataTable(read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character"));
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));      
		output$dowhat<-renderText("Deleted compound");
		if(any(ls()=="logfile")){stop("illegal logfile detected #1 in server_obs_Add.r!")}
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
		targets2[18]<-as.character(isolate(input$warn_1)) 
		targets2[19]<-as.character(isolate(input$warn_2))	
		targets<-rbind(targets2,targets1);
		write.table(targets,file=file.path(logfile[[1]],"dataframes","targets.txt"),row.names=FALSE,sep="\t",quote=FALSE)      
		rm(targets,targets1,targets2);
		#############################################################################
		# adjust task/workflow settings #############################################
		enviMass:::workflow_set(logfile,down="pattern")	
		#############################################################################			
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));      
		output$targets<<-DT::renderDataTable(read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character"));      
		output$dowhat<-renderText("Added target compound");
		if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_obs_Add.r!")}
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
			workflow_set(logfile,down="pattern")	
			#############################################################################			
			save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));   
 			output$targets<<-DT::renderDataTable(read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character"));
			output$dowhat<-renderText("Added targets list");
			if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_obs_Add.r!")}
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
		enviMass:::workflow_set(logfile,down="pattern")	
		#############################################################################			
		output$targets<<-DT::renderDataTable(read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character"));      	  
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));      
		output$dowhat<-renderText("Deleted compound");
		if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_obs_Add.r!")}
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
							"TRUE","FALSE","FALSE","FALSE","FALSE","FALSE",
							as.character(isolate(input$Measadd_profiled)),
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
							enviMass:::workflow_set(logfile,down="peakpicking",single_file=TRUE)	
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
							"TRUE","FALSE","FALSE","FALSE","FALSE","FALSE",
							as.character(isolate(input$Measadd_profiled)),
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
						enviMass:::workflow_set(logfile,down="peakpicking",single_file=TRUE)	
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
			#########################################################################			
			# subtraction files, positive: ##########################################
			measurements3<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			if(any( (measurements3[,1]!="-") & (measurements3[,4]=="positive") & (measurements3[,3]!="sample"))){
				IDs_pos<-measurements3[
					(measurements3[,4]=="positive") & (measurements3[,3]!="sample")
				,1]
				names_pos<-measurements3[
					(measurements3[,4]=="positive") & (measurements3[,3]!="sample")
				,2]
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
			if(any( (measurements3[,1]!="-") & (measurements3[,4]=="negative") & (measurements3[,3]!="sample"))){
				IDs_neg<-measurements3[
					(measurements3[,4]=="negative") & (measurements3[,3]!="sample")
				,1]
				names_neg<-measurements3[
					(measurements3[,4]=="negative") & (measurements3[,3]!="sample")
				,2]
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
			"checked?","recal?","align?","norm?","profiled?","comp?","IS_screen?","tar_screen?","tag1","tag2","tag3")
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
			enviMass:::workflow_set(logfile,down="peakpicking",single_file=TRUE)	
		}	
		#########################################################################			
		# subtraction files, positive: ##########################################
		measurements3<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		if(any( (measurements3[,1]!="-") & (measurements3[,4]=="positive") & (measurements3[,3]!="sample"))){
			IDs_pos<-measurements3[
				(measurements3[,4]=="positive") & (measurements3[,3]!="sample")
			,1]
			names_pos<-measurements3[
				(measurements3[,4]=="positive") & (measurements3[,3]!="sample")
			,2]
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
		if(any( (measurements3[,1]!="-") & (measurements3[,4]=="negative") & (measurements3[,3]!="sample"))){
			IDs_neg<-measurements3[
				(measurements3[,4]=="negative") & (measurements3[,3]!="sample")
			,1]
			names_neg<-measurements3[
				(measurements3[,4]=="negative") & (measurements3[,3]!="sample")
			,2]
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
		#############################################################################
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
        #############################################################################			
        output$measurements<<-DT::renderDataTable(read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character"));
        save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));      
        output$dowhat<-renderText("Measurement deleted");
      }else{
        output$dowhat<-renderText("Invalid ID");
      }
    }
	if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_obs_Add.r!")}
})
##############################################################################

# IMPORT MEASUREMENTS ########################################################
observe({
    input$Import_project
    if(input$Import_project){
		cat("\n Importing project files ...")
        measurements_1<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		file_in<<-as.character(isolate(input$import_pro_dir))
		measurements_2<-read.csv(file=file.path(file_in,"dataframes","measurements"),colClasses = "character");
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
			enviMass:::workflow_set(logfile,down="peakpicking",single_file=TRUE) 
			#########################################################################			
			# subtraction files, positive: ##########################################
			measurements3<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			if(any( (measurements3[,1]!="-") & (measurements3[,4]=="positive") & (measurements3[,3]!="sample"))){
				IDs_pos<-measurements3[
					(measurements3[,4]=="positive") & (measurements3[,3]!="sample")
				,1]
				names_pos<-measurements3[
					(measurements3[,4]=="positive") & (measurements3[,3]!="sample")
				,2]
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
			if(any( (measurements3[,1]!="-") & (measurements3[,4]=="negative") & (measurements3[,3]!="sample"))){
				IDs_neg<-measurements3[
					(measurements3[,4]=="negative") & (measurements3[,3]!="sample")
				,1]
				names_neg<-measurements3[
					(measurements3[,4]=="negative") & (measurements3[,3]!="sample")
				,2]
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
			output$summa_html<<-renderText(enviMass:::summary_html(logfile$summary));
			save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
			#########################################################################			
			output$dowhat<-renderText("Files imported.");
			cat(" done.")		
		}else{
			cat(" no files to import - project empty?.")
			output$dowhat<-renderText("Failed import: no files.");
		}
	}
	if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_obs_Add.r!")}
})
##############################################################################

# MODIFY MEASUREMENTS ########################################################  
# LOAD
observe({
	input$Modif_load
	if(isolate(input$Modif_load)){
		measurements3<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		atID<-as.character(isolate(input$Modif_ID))
		if(any(measurements3[,1]==atID)){
			updateTextInput(session, "Modif_name",value = as.character(measurements3[measurements3[,1]==atID,2]))
			updateSelectInput(session,"Modif_type","Type:", choices = c("sample", "blank", "doted", "other"), selected = as.character(measurements3[measurements3[,1]==atID,3]))
			updateSelectInput(session, "Modif_incl", selected = as.character(measurements3[measurements3[,1]==atID,8]))		
			updateSelectInput(session, "Modif_mode", selected = as.character(measurements3[measurements3[,1]==atID,4]))
			updateTextInput(session, "Modif_place",value = as.character(measurements3[measurements3[,1]==atID,5]))
			updateDateInput(session, "Modif_date", value = as.character(measurements3[measurements3[,1]==atID,6]))
			updateTextInput(session, "Modif_time",value = as.character(measurements3[measurements3[,1]==atID,7]))		
			updateTextInput(session, "Modif_tag3",value = as.character(measurements3[measurements3[,1]==atID,21]))
			updateSelectInput(session, "Modif_profiled", selected = as.character(measurements3[measurements3[,1]==atID,15]))	
			output$dowhat<-renderText("Specifications loaded into mask.");
			cat("\n specifications loaded into mask")
			rm(measurements3)
		}else{
			updateTextInput(session, "Modif_name",value = "INVALID ID")		
			updateTextInput(session, "Modif_place",value = "INVALID ID")
			updateTextInput(session, "Modif_tag3",value = "INVALID ID")			
			
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
			measurements3[measurements3[,1]==atID,2]<-as.character(isolate(input$Modif_name))
			measurements3[measurements3[,1]==atID,3]<-as.character(isolate(input$Modif_type))
			measurements3[measurements3[,1]==atID,8]<-as.character(isolate(input$Modif_incl))
			measurements3[measurements3[,1]==atID,4]<-as.character(isolate(input$Modif_mode))
			measurements3[measurements3[,1]==atID,5]<-as.character(isolate(input$Modif_place))
			measurements3[measurements3[,1]==atID,6]<-as.character(isolate(input$Modif_date))
			measurements3[measurements3[,1]==atID,]<-enviMass:::convDate(measurements3[measurements3[,1]==atID,]);
			measurements3[measurements3[,1]==atID,7]<-as.character(isolate(input$Modif_time))		
			measurements3[measurements3[,1]==atID,21]<-as.character(isolate(input$Modif_tag3))	
			measurements3[measurements3[,1]==atID,15]<-as.character(isolate(input$Modif_profiled))	
			enviMass:::workflow_set(down="peakpicking",check_node=TRUE,single_file=TRUE)
			write.csv(measurements3,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
			output$dowhat<-renderText("Specifications saved to file table.");
			cat("\n specifications exported from mask to file table")
			rm(measurements3)
			output$measurements<<-DT::renderDataTable(read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character")); 
			######################################################################			
			# subtraction files, positive: #######################################
			measurements3<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			if(any( (measurements3[,1]!="-") & (measurements3[,4]=="positive") & (measurements3[,3]!="sample"))){
				IDs_pos<-measurements3[
					(measurements3[,4]=="positive") & (measurements3[,3]!="sample")
				,1]
				names_pos<-measurements3[
					(measurements3[,4]=="positive") & (measurements3[,3]!="sample")
				,2]
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
			if(any( (measurements3[,1]!="-") & (measurements3[,4]=="negative") & (measurements3[,3]!="sample"))){
				IDs_neg<-measurements3[
					(measurements3[,4]=="negative") & (measurements3[,3]!="sample")
				,1]
				names_neg<-measurements3[
					(measurements3[,4]=="negative") & (measurements3[,3]!="sample")
				,2]
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
			######################################################################
			output$summa_html<<-renderText(enviMass:::summary_html(logfile$summary));
			save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
			######################################################################
		}
	}
})  
##############################################################################  
  
############################################################################## 
# Import parameters ##########################################################  
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
		logfile[[4]]<<-logfile_other[[4]]
		logfile[[5]]<<-logfile_other[[5]]
		logfile[[7]]<<-logfile_other[[7]]		
		logfile[[8]]<<-logfile_other[[8]]		 
 		logfile[[9]]<<-logfile_other[[9]]		
		rm(logfile_other,logfile_here,envir=as.environment(".GlobalEnv"))
 		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp")); 
		measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		source("server_variables_in.R", local=TRUE)
		output$dowhat<<-renderText("Parameters imported.");
		cat(" done. \n")
	}
	if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_obs_Add.r!")}
}) 
############################################################################## 
 
  
  
  
  
  
  
  