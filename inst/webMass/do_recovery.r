
	#######################################################################################################################
	# delete old results ##################################################################################################
	those<-list(0)
	those[[1]]<-file.path(logfile[[1]],"quantification","target_recov_table_pos")			
	those[[2]]<-file.path(logfile[[1]],"quantification","target_recov_table_neg")	
	for(n in 1:length(those)){
		if(file.exists(those[[n]])){
			file.remove(those[[n]])
		}
	}	
	rm(those)
	measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");		
	######################################################################################################################
	
	######################################################################################################################	
	# POSITIVE ###########################################################################################################
	if(
		any(measurements[,"Mode"]=="positive" & measurements[,"Type"]=="spiked" & measurements[,"include"]=="TRUE") &
		file.exists(file.path(logfile[[1]],"quantification","target_quant_table_pos")) 
	){
		
		load(file.path(logfile[[1]],"quantification","target_quant_table_pos"))
		those_files<-measurements[(measurements[,"Mode"]=="positive" & measurements[,"Type"]=="spiked" & measurements[,"include"]=="TRUE"),,drop=FALSE]
		atdate<-those_files[,6]
		atdate<-as.Date(atdate);
		attime<-those_files[,7]
		attime<-as.difftime(attime);
		ord<-order(as.numeric(atdate),as.numeric(attime),as.numeric(those_files[,1]),decreasing=TRUE);
		those_files<-those_files[ord,]	
		if(logfile$parameters$recov_files_included!="FALSE"){
			if(logfile$parameters$recov_files_included<length(those_files[,1])){
				those_files<-those_files[1:logfile$parameters$recov_files_included,]		
			}
		}
		those_targets<-target_quant_table_pos[6:length(target_quant_table_pos[,1]),1:2]
		target_recov_table_pos<-matrix(nrow=(length(those_targets[,1])+4),ncol=(length(those_files[,1])+2),"")
		colnames(target_recov_table_pos)<-c("Target ID","Target name",those_files[,"ID"])
		rownames(target_recov_table_pos)<-c("Name","Type","Date","Time",those_targets[,1])
		target_recov_table_pos[1,]<-c("","",as.character(those_files[,"Name"]))
		target_recov_table_pos[2,]<-c("","",as.character(those_files[,"Type"]))
		target_recov_table_pos[3,]<-c("","",as.character(those_files[,"Date"]))
		target_recov_table_pos[4,]<-c("","",as.character(those_files[,"Time"]))		
		target_recov_table_pos[,1]<-c("","","","",those_targets[,1])
		target_recov_table_pos[,2]<-c("","","","",those_targets[,2])		
		##################################################################################################################
		for(i in 1:length(those_files[,"ID"])){
				from_ID<-those_files[i,"ID"]
				to_ID<-those_files[i,"tag2"]
				if(!any(measurements[measurements[,"Mode"]=="positive","ID"]==to_ID)){ # this should not happen anyway - included in check_project
					cat("\n WARNING: Missing relation for spiked file detected! Please revise");
					next;
				}
				if(!any(colnames(target_quant_table_pos)==from_ID)){
					next;
				}
				for(j in 5:length(target_recov_table_pos[,1])){
					target_ID<-target_recov_table_pos[j,1]
					from_quant<-target_quant_table_pos[
						target_quant_table_pos[,1]==target_ID,
						colnames(target_quant_table_pos)==from_ID
					]		
					if(grepl("!",from_quant)){next}
					to_quant<-target_quant_table_pos[
						target_quant_table_pos[,1]==target_ID,
						colnames(target_quant_table_pos)==to_ID
					]
					if(grepl("!",to_quant)){next}
					from_quant<-as.numeric(strsplit(from_quant,",")[[1]])
					to_quant<-as.numeric(strsplit(to_quant,",")[[1]])				
					recov<-c()
					for(n in 1:length(from_quant)){
						for(m in 1:length(to_quant)){
							recov<-c(recov,
								from_quant[n]-to_quant[m]
							)
						}
					}
					recov<-recov[recov>=0] # cannot be negatively concentrated!
					if(length(recov)==0){next}
					recov<-paste(as.character(recov),collapse=",")
					target_recov_table_pos[
						target_recov_table_pos[,1]==target_ID,
						colnames(target_recov_table_pos)==from_ID
					]<-recov				
				}			
		}
		##################################################################################################################
		save(target_recov_table_pos,file=file.path(logfile[[1]],"quantification","target_recov_table_pos"))
		rm(target_quant_table_pos,target_recov_table_pos)
	}
	######################################################################################################################

	######################################################################################################################	
	# NEGATIVE ###########################################################################################################
	if(
		any(measurements[,"Mode"]=="negative" & measurements[,"Type"]=="spiked" & measurements[,"include"]=="TRUE") &
		file.exists(file.path(logfile[[1]],"quantification","target_quant_table_neg")) 
	){
		
		load(file.path(logfile[[1]],"quantification","target_quant_table_neg"))
		those_files<-measurements[(measurements[,"Mode"]=="negative" & measurements[,"Type"]=="spiked" & measurements[,"include"]=="TRUE"),,drop=FALSE]
		atdate<-those_files[,6]
		atdate<-as.Date(atdate);
		attime<-those_files[,7]
		attime<-as.difftime(attime);
		ord<-order(as.numeric(atdate),as.numeric(attime),as.numeric(those_files[,1]),decreasing=TRUE);
		those_files<-those_files[ord,]	
		if(logfile$parameters$recov_files_included!="FALSE"){
			if(logfile$parameters$recov_files_included<length(those_files[,1])){
				those_files<-those_files[1:logfile$parameters$recov_files_included,]		
			}
		}
		those_targets<-target_quant_table_neg[6:length(target_quant_table_neg[,1]),1:2]
		target_recov_table_neg<-matrix(nrow=(length(those_targets[,1])+4),ncol=(length(those_files[,1])+2),"")
		colnames(target_recov_table_neg)<-c("Target ID","Target name",those_files[,"ID"])
		rownames(target_recov_table_neg)<-c("Name","Type","Date","Time",those_targets[,1])
		target_recov_table_neg[1,]<-c("","",as.character(those_files[,"Name"]))
		target_recov_table_neg[2,]<-c("","",as.character(those_files[,"Type"]))
		target_recov_table_neg[3,]<-c("","",as.character(those_files[,"Date"]))
		target_recov_table_neg[4,]<-c("","",as.character(those_files[,"Time"]))		
		target_recov_table_neg[,1]<-c("","","","",those_targets[,1])
		target_recov_table_neg[,2]<-c("","","","",those_targets[,2])		
		##################################################################################################################
		for(i in 1:length(those_files[,"ID"])){
				from_ID<-those_files[i,"ID"]
				to_ID<-those_files[i,"tag2"]
				if(!any(measurements[measurements[,"Mode"]=="negative","ID"]==to_ID)){ # this should not happen anyway - included in check_project
					cat("\n WARNING: Missing relation for spiked file detected! Please revise");
					next;
				}
				if(!any(colnames(target_quant_table_neg)==from_ID)){
					next;
				}
				for(j in 5:length(target_recov_table_neg[,1])){
					target_ID<-target_recov_table_neg[j,1]
					from_quant<-target_quant_table_neg[
						target_quant_table_neg[,1]==target_ID,
						colnames(target_quant_table_neg)==from_ID
					]		
					if(grepl("!",from_quant)){next}
					to_quant<-target_quant_table_neg[
						target_quant_table_neg[,1]==target_ID,
						colnames(target_quant_table_neg)==to_ID
					]
					if(grepl("!",to_quant)){next}
					from_quant<-as.numeric(strsplit(from_quant,",")[[1]])
					to_quant<-as.numeric(strsplit(to_quant,",")[[1]])				
					recov<-c()
					for(n in 1:length(from_quant)){
						for(m in 1:length(to_quant)){
							recov<-c(recov,
								from_quant[n]-to_quant[m]
							)
						}
					}
					recov<-recov[recov>=0] # cannot be negatively concentrated!
					if(length(recov)==0){next}
					recov<-paste(as.character(recov),collapse=",")
					target_recov_table_neg[
						target_recov_table_neg[,1]==target_ID,
						colnames(target_recov_table_neg)==from_ID
					]<-recov				
				}				
		}
		##################################################################################################################
		save(target_recov_table_neg,file=file.path(logfile[[1]],"quantification","target_recov_table_neg"))
		rm(target_quant_table_neg,target_recov_table_neg)
	}
	######################################################################################################################
	
	######################################################################################################################
	rm(measurements)
	
	
	