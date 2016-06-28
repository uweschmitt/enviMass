#' @title Check enviMass project
#'
#' @description \code{} runs several cconsistency checks on a enviMass project
#'
#' @param logfile enviMass project logfile
#' @param isotopes Isotope list
#' @param adducts Adducts list from package enviPat
#' @param skipcheck Logical. Should project check be skipped?
#'
#' @return Character string with either specific error message or a "Project consistent" message.
#' 
#' @details enviMass workflow function; run before further calculations are started in the workflow.
#' 

checkproject<-function(isotopes,adducts,skipcheck=FALSE,...){
  say<-"Project consistent"
  if(skipcheck){
	return(say);
  }
  if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in check_project.r!")}
  ###############################################################################
  # wrong upstream "must" executions? ###########################################
  must<-logfile[[12]]
  for(i in 1:length(must[1,])){
	for(j in 1:length(must[,i])){	
		if(must[j,i]==-1){
			if(logfile$workflow[names(logfile$workflow)==colnames(must)[i]]=="yes"){
				if(logfile$workflow[names(logfile$workflow)==rownames(must)[j]]=="yes"){
					say<-paste("wokflow step",names(logfile$workflow)[i],"excludes",names(logfile$workflow)[j],"- adapt workflow settings!")
				}
			}
		}		
	}
  }
  ##############################################################################
  # directories available? ##################################################### 
  if(!file.exists(file.path(logfile[[1]],"files"))){say<-"files directory missing!"}
  if(!file.exists(file.path(logfile[[1]],"MSlist"))){say<-"MSlist directory missing!"}  
  if(!file.exists(file.path(logfile[[1]],"features"))){say<-"features directory missing!"}
  if(!file.exists(file.path(logfile[[1]],"results"))){say<-"results directory missing!"}
	if(!file.exists(file.path(logfile[[1]],"results","screening"))){say<-"results/screening directory missing!"} 
	if(!file.exists(file.path(logfile[[1]],"quantification"))){say<-"results/quantification directory missing!"} 
	if(!file.exists(file.path(logfile[[1]],"results","LOD"))){say<-"results/LOD directory missing!"} 
	if(!file.exists(file.path(logfile[[1]],"results","recalibration"))){say<-"results/recalibration directory missing!"} 	
  if(!file.exists(file.path(logfile[[1]],"dataframes"))){say<-"dataframes directory missing!"}
  if(!file.exists(file.path(logfile[[1]],"pics"))){say<-"pics directory missing!"}
  if(!file.exists(file.path(logfile[[1]],"exports"))){say<-"exports directory missing!"}  
  ##############################################################################
  # compounds available & ok? ##################################################
	intstand_check<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character",blank.lines.skip=TRUE);
	if(length(intstand_check$Formula[intstand_check$Formula!="-"])==0){
		cat("No internal standards available \n") 
	}else{
		if(
			!all(names(intstand_check)==c("ID","Name","Formula","RT","RT_tolerance","main_adduct","ion_mode",
			"use_for_recalibration","use_for_screening", "restrict_adduct","Remark","tag1","tag2","tag3",
			"from","to","Lower_intensity_bound","Upper_intensity_bound","Quant_adduct","Quant_peak"))
		){
			say<-paste("Incorrect or missing columns in internal standard compound table - compare to tables in a new project for correct entries.")
		}
		if(any(grepl(",",intstand_check[,names(intstand_check)=="Lower_intensity_bound"],fixed=TRUE))){
			say<-"Correct IS table: use .-separator, no commas for bounds."
		}
		if(any(grepl(",",intstand_check[,names(intstand_check)=="Upper_intensity_bound"],fixed=TRUE))){
			say<-"Correct IS table: use .-separator, no commas for bounds."
		}
		if(any(duplicated(intstand_check[,1]))){
			those<-intstand_check[duplicated(intstand_check[,1]),1]
			those<-as.character(those)
			those<-paste(those,collapse=",")
			say<-paste("Duplicated IS IDs found:",those,". Please correct!")
		}  
		if(TRUE){ # ID check per positive OR negative - not run
			if(any(intstand_check[,names(intstand_check)=="ion_mode"]=="positive")){
				if(any(duplicated(intstand_check[intstand_check[,names(intstand_check)=="ion_mode"]=="positive",1]))){
					those<-intstand_check[duplicated(intstand_check[intstand_check[,names(intstand_check)=="ion_mode"]=="positive",1]),1]
					those<-as.character(those)
					those<-paste(those,collapse=",")
					say<-paste("Duplicated IS IDs found (positive ion_mode):",those,". Please correct!")
				}  
			}
			if(any(intstand_check[,names(intstand_check)=="ion_mode"]=="negative")){
				if(any(duplicated(intstand_check[intstand_check[,names(intstand_check)=="ion_mode"]=="negative",1]))){
					those<-intstand_check[duplicated(intstand_check[intstand_check[,names(intstand_check)=="ion_mode"]=="negative",1]),1]
					those<-as.character(those)
					those<-paste(those,collapse=",")
					say<-paste("Duplicated IS IDs found (negative ion_mode):",those,". Please correct!")
				}  
			}
		}
		checked<-enviPat::check_chemform(isotopes, intstand_check[,3])
		if(any(checked[,1])){
			print(intstand_check[intstand_check[,1]==TRUE,1])
			say<-"Invalid internal standard formula detected ... abort. Check R console for concerned target IDs"
		}else{ # save corrected formulas!
			intstand_check[,3]<-checked[,2]
			write.table(intstand_check,file=file.path(logfile[[1]],"dataframes","IS.txt"),row.names=FALSE,sep="\t",quote=FALSE)
		}
		if(any(!charmatch(intstand_check[,7],c("positive","negative"),nomatch=FALSE))){say<-"Invalid IS column 7 entry (not positive/negative)"}
		if(any(!charmatch(intstand_check[intstand_check[,7]=="positive",6],c(as.character(adducts[adducts[,6]=="positive",1]),"FALSE"),nomatch=FALSE))){
			wrongadduct<-intstand_check[intstand_check[,7]=="positive",6][!charmatch(intstand_check[intstand_check[,7]=="positive",6],c(as.character(adducts[adducts[,6]=="positive",1]),"FALSE"),nomatch=FALSE)]
			wrongadduct<-wrongadduct[1]
			say<-paste("IS column 6: wrong adduct",wrongadduct)
		}
		if(any(!charmatch(intstand_check[intstand_check[,7]=="negative",6],c(as.character(adducts[adducts[,6]=="negative",1]),"FALSE"),nomatch=FALSE))){
			wrongadduct<-intstand_check[intstand_check[,7]=="negative",6][!charmatch(intstand_check[intstand_check[,7]=="negative",6],c(as.character(adducts[adducts[,6]=="negative",1]),"FALSE"),nomatch=FALSE)]
			wrongadduct<-wrongadduct[1]
			say<-paste("IS column 6: wrong adduct",wrongadduct)
		}
		if(any(is.na(as.numeric(intstand_check[,4])))){say<-"IS column 4 not numeric"}
		if(any(is.na(as.numeric(intstand_check[intstand_check[,5]!="FALSE",5])))){say<-"IS compound table column 5 not numeric or an empty row is present!"}
		if(any(!charmatch(intstand_check[,8],c("TRUE","FALSE"),nomatch=FALSE))){say<-"IS compound table column 8 not logical or an empty row is present!"}
		if(any(!charmatch(intstand_check[,9],c("TRUE","FALSE"),nomatch=FALSE))){say<-"IS compound table column 9 not logical or an empty row is present!"}
		if(any(!(charmatch(intstand_check[,10],c("TRUE","FALSE"),nomatch=FALSE))) ){say<-"IS compound table column 10 not logical or an empty row is present!"}
		rm(checked,intstand_check)
	}
	targets_check<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character",blank.lines.skip=TRUE);
	if(length(targets_check$Formula[targets_check$Formula!="-"])==0){
		  cat("No internal targets available \n")    
	}else{
		if(
			!all(names(targets_check)==c("ID","Name","Formula","RT","RT_tolerance","ID_internal_standard","main_adduct","ion_mode",
			"use_for_recalibration","use_for_screening","restrict_adduct","Remark","tag1","tag2","tag3","from","to","warn_1","warn_2","Quant_adduct","Quant_peak"))
		){
			say<-paste("Incorrect or missing columns in target compound table - compare to tables from a new project for correct entries.")
		}
		if(any(grepl(",",targets_check[,names(targets_check)=="warn_1"],fixed=TRUE))){
			say<-"Correct targets table: use .-separator, no commas for warn_1."
		}
		if(any(grepl(",",targets_check[,names(targets_check)=="warn_2"],fixed=TRUE))){
			say<-"Correct targets table: use .-separator, no commas for warn_2."
		}		
		if(any(duplicated(targets_check[,1]))){
			those<-targets_check[duplicated(targets_check[,1]),1]
			those<-as.character(those)
			those<-paste(those,collapse=",")
			say<-paste("Duplicated target IDs found:",those,". Please correct!")		
		}  
		checked<-enviPat::check_chemform(isotopes, targets_check[,3])
		if(any(checked[,1])){
			print(targets_check[checked[,1]==TRUE,1])
			say<-"Invalid target formula detected ... abort. Check R console for concerned target IDs"
		}else{
			targets_check[,3]<-checked[,2]
			write.table(targets_check,file=file.path(logfile[[1]],"dataframes","targets.txt"),row.names=FALSE,sep="\t",quote=FALSE)        
		}  
		if(any(!charmatch(targets_check[,8],c("positive","negative"),nomatch=FALSE))){say<-"Invalid target column 8 entry (not positive/negative)"}
		if(any(!charmatch(targets_check[targets_check[,8]=="positive",7],c(as.character(adducts[adducts[,6]=="positive",1]),"FALSE"),nomatch=FALSE))){
			wrongadduct<-targets_check[targets_check[,8]=="positive",7][!charmatch(targets_check[targets_check[,8]=="positive",7],c(as.character(adducts[adducts[,6]=="positive",1]),"FALSE"),nomatch=FALSE)]
			wrongadduct<-wrongadduct[1]
			say<-paste("target column 7: wrong adduct",wrongadduct)
		}
		if(any(!charmatch(targets_check[targets_check[,8]=="negative",7],c(as.character(adducts[adducts[,6]=="negative",1]),"FALSE"),nomatch=FALSE))){
			wrongadduct<-targets_check[targets_check[,8]=="negative",7][!charmatch(targets_check[targets_check[,8]=="negative",7],c(as.character(adducts[adducts[,6]=="negative",1]),"FALSE"),nomatch=FALSE)]		
			wrongadduct<-wrongadduct[1]
			say<-paste("target column 7: wrong adduct",wrongadduct)
		}
		if(any(is.na(as.numeric(targets_check[,4])))){say<-"targets column 4 not numeric"}
		if(any(is.na(as.numeric(targets_check[targets_check[,5]!="FALSE",5])))){say<-"target compound table column 5 not numeric or an empty row is present!"}
		if(any(!charmatch(targets_check[,9],c("TRUE","FALSE"),nomatch=FALSE))){say<-"target compound table column 9 not logical or an empty row is present!"}
		if(any(!charmatch(targets_check[,10],c("TRUE","FALSE"),nomatch=FALSE))){say<-"target compound table column 10 not logical or an empty row is present!"}
		if(any(!(charmatch(targets_check[,11],c("TRUE","FALSE"),nomatch=FALSE))) ){say<-"target compound table column 11 not logical or an empty row is present!"}	
		rm(checked,targets_check)  
	}   
	##############################################################################
	# check compounds for calibration & quantification ###########################
	if(
		(logfile$workflow[names(logfile$workflow)=="calibration"])=="yes" ||
		(logfile$workflow[names(logfile$workflow)=="quantification"])=="yes"	
	){
		targets_check<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character",blank.lines.skip=TRUE);
		intstand_check<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character",blank.lines.skip=TRUE);
		# check if all relations / adducts are correct
		for(i in 1:length(targets_check[,1])){
			if(targets_check[i,6]!="FALSE"){
				ID_target_missing<-targets_check[i,1];
				ID_IS_missing<-targets_check[i,6];
				target_quan_adduct<-targets_check[i,names(targets_check)=="Quant_adduct"]				
				target_avail_adduct<-targets_check[i,names(targets_check)=="main_adduct"]
				if(targets_check[i,names(targets_check)=="restrict_adduct"]=="FALSE"){
					target_avail_adduct<-c(target_avail_adduct,logfile[[7]],logfile[[8]])
				}
				if( (!any(target_avail_adduct==target_quan_adduct)) | (target_quan_adduct=="FALSE") ){ # does the target quantification adduct exist?
					say<-paste(
					"Adduct used in quantification for target compound with ID ",
					ID_target_missing," not valid. Either it is set to FALSE and shouldn`t or the required adduct is not included in the adduct selection of this target! Please adjust.",sep="")
				}
				found_standard<-TRUE
				if( !any(intstand_check[ intstand_check[,7]==targets_check[i,8],1]==targets_check[i,6]) ){ # does the internal standard exist?		
					say<-paste(
					"Quantification/calibration problem: internal standard with ID ",
					ID_IS_missing," for target compound with ID ",ID_target_missing,
					" not found at this ionization mode. Check compound target table, column ID_internal_standard.",
					sep="")
					found_standard<-FALSE
				}
				if(found_standard){ # does the standard quantification adduct exist?
					j<-(intstand_check[,1]==ID_IS_missing)
					IS_quan_adduct<-intstand_check[j,names(intstand_check)=="Quant_adduct"]				
					IS_avail_adduct<-intstand_check[j,names(intstand_check)=="main_adduct"]	
					if(intstand_check[j,names(intstand_check)=="restrict_adduct"]=="FALSE"){
						IS_avail_adduct<-c(IS_avail_adduct,logfile[[7]],logfile[[8]])
					}
					if( (!any(IS_avail_adduct==IS_quan_adduct)) | (IS_quan_adduct=="FALSE") ){ # does the target quantification adduct exist?
						say<-paste(
						"Adduct used in quantification for internal standard with ID ",
						ID_IS_missing," not valid. Either it is set to FALSE and shouldn`t or the required adduct is not included in the adduct selection of this internal standard! Please adjust.",sep="")
					}
				}
			}
		}	
	}
  # enough compounds for recalibration available? ##############################
  if(logfile[[6]][2]=="yes"){
    if(logfile[[5]][30]=="Internal standards"){
 	  IS<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
	  if(length(IS[IS[,8]=="TRUE",1])<10){
        say<-"Not enough internal standards available for mass recalibration ... abort."    
      }
    }
    if(logfile[[5]][30]=="Target compounds"){
	  targets<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");
	  if(length(targets[targets[,9]=="TRUE",1])<10){
        say<-"Not enough target compounds available for mass recalibration ... abort."    
      }
    }
    if(logfile[[5]][30]=="both"){
  	  IS<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");	  
	  a<-length(IS[IS[,8]=="TRUE",1])
 	  targets<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");  
	  b<-length(targets[targets[,9]=="TRUE",1])
      if((a<10)||(b<10)){
        say<-"Not enough target compounds + internal standards available for mass recalibration ... abort."    
      }
    }
  }
  ##############################################################################
  # parameters ok? #############################################################
  # (1) on trend time lags #####################################################
  if(logfile$workflow[names(logfile$workflow)=="trendblind"]=="yes"){
	  lags<-as.numeric(strsplit(as.character(logfile[[5]][34]),",")[[1]])
	  if(any(is.na(lags))){say<-"Invalid trend lags - have you used comma separated numerics?"}
	  measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	  if(!any(measurements[,8]=="TRUE")){
		say<-"No file included into workflow; nothing to be calculated."
		return(say)
	  }
	  measurements<-measurements[measurements[,8]=="TRUE",]
	  dated<-measurements[,6]
	  timed<-measurements[,7]
	  datetime<-c()
	  if (length(timed) > 0) {
		for(i in 1:length(timed)){
			datetime<-c(datetime,paste(dated[i],timed[i],"CET",sep=" "))
		}
		atPOSIX<-as.POSIXct(datetime);
		atPOSIX<-as.numeric(atPOSIX)
		if(min(lags)>(((max(atPOSIX)-min(atPOSIX))/(24*60*60))+1)){say<-"Trend lags longer than time span of the measurements ... abort"}
	  }
	  rm(lags); ##############################################################################
	  if(
		(logfile$workflow[names(logfile$workflow)=="profiling"]=="yes") &
		( !any(measurements[,names(measurements)=="profiled"]=="TRUE") )
	  ){
			say<-"Workflow option profiling enabled, but no file included for profiling."
	  }
  }
  # data sets ok? ##############################################################
  filed<-list.files(file.path(logfile[[1]],"files"))
  if(!length(filed)){say<-"No measurements available!"}
  measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character")
  if(!all(!duplicated(measurements[,1]))){say<-"Duplicated measurements IDs. Revise!"}
  #measurements<-measurements[!duplicated(measurements[,1]),]
  if(any(is.na(as.numeric(measurements[,1])))){say<-"Non-numeric measurements IDs. Revise!"}  
  measurements_ID<-measurements[,1]
  measurements_ID<-paste(measurements_ID,".mzXML",sep="")
  if(any(match(measurements_ID,filed,nomatch=0)==0)){say<-paste("Missing mzXML file - file corrupted? Compare project file folder for consistency!",sep="")}
  ##############################################################################
  # progress bar? ##############################################################
  if(interactive() && !.Platform$OS.type == "windows" && .Platform$GUI == "Rgui" && logfile[[5]][21]=="TRUE"){
	say<-"Disable the progress bar in the Settings General Tab; works only under Windows OS"
  }
  ##############################################################################
  # blind peak subtraction enabled, but no blind files selected? ###############
  if(
	(logfile$workflow[names(logfile$workflow)=="blinds"]=="yes") &
	(logfile$parameters[[85]]==FALSE) &
	(logfile$parameters[[86]]==FALSE) &
	(logfile$parameters[[87]]==FALSE) &
	(logfile$parameters[[88]]==FALSE) 
  ){
	say<-"Blind detection enabled but blind settings disabled? Please adjust."
  }
  ##############################################################################
  if(any(ls()=="logfile")){stop("\n illegal logfile detected #2 in check_project.r!")}
  return(say);
  ##############################################################################  
}
