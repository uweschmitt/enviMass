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

checkproject<-function(isotopes,adducts,skipcheck=FALSE){
  say<-"Project consistent"
  if(skipcheck){
	return(say);
  }

  if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in check_project.r!")}
  ##############################################################################
  # directories available? ##################################################### 
  if(!file.exists(file.path(logfile[[1]],"files"))){say<-"files directory missing!"}
  if(!file.exists(file.path(logfile[[1]],"MSlist"))){say<-"MSlist directory missing!"}  
  if(!file.exists(file.path(logfile[[1]],"features"))){say<-"features directory missing!"}
  if(!file.exists(file.path(logfile[[1]],"results"))){say<-"results directory missing!"}
	if(!file.exists(file.path(logfile[[1]],"results","screening"))){say<-"results/screening directory missing!"} 
	if(!file.exists(file.path(logfile[[1]],"results","quantification"))){say<-"results/quantification directory missing!"} 
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
		if(any(duplicated(intstand_check[,1]))){say<-"Duplicated IS IDs found ... abort."}  
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
		if(any(is.na(as.numeric(intstand_check[intstand_check[,5]!="FALSE",5])))){say<-"IS compound table column 5 not numeric"}
		if(any(!charmatch(intstand_check[,8],c("TRUE","FALSE"),nomatch=FALSE))){say<-"IS compound table column 8 not logical"}
		if(any(!charmatch(intstand_check[,9],c("TRUE","FALSE"),nomatch=FALSE))){say<-"IS compound table column 9 not logical"}
		if(any(!(charmatch(intstand_check[,10],c("TRUE","FALSE"),nomatch=FALSE))) ){say<-"IS compound table column 10 not logical"}
		rm(checked,intstand_check)
	}
	targets_check<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character",blank.lines.skip=TRUE);
	if(length(targets_check$Formula[targets_check$Formula!="-"])==0){
		  cat("No internal targets available \n")    
	}else{
		if(any(duplicated(targets_check[,1]))){say<-"Duplicated target IDs found ... abort."}  
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
		if(any(is.na(as.numeric(targets_check[,4])))){say<-"IS column 4 not numeric"}
		if(any(is.na(as.numeric(targets_check[targets_check[,5]!="FALSE",5])))){say<-"target compound table column 5 not numeric"}
		if(any(!charmatch(targets_check[,9],c("TRUE","FALSE"),nomatch=FALSE))){say<-"target compound table column 9 not logical"}
		if(any(!charmatch(targets_check[,10],c("TRUE","FALSE"),nomatch=FALSE))){say<-"target compound table column 10 not logical"}
		if(any(!(charmatch(targets_check[,11],c("TRUE","FALSE"),nomatch=FALSE))) ){say<-"target compound table column 11 not logical"}	
		rm(checked,targets_check)  
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
  if(logfile[[2]][[7]]){
	  lags<-as.numeric(strsplit(as.character(logfile[[5]][34]),",")[[1]])
	  if(any(is.na(lags))){say<-"Invalid trend lags - have you used comma separated numerics?"}
	  measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
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
  }
  # data sets ok? ##############################################################
  filed<-list.files(file.path(logfile[[1]],"files"))
  if(!length(filed)){say<-"No measurements available!"}
  measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
  if(!all(!duplicated(measurements[,1]))){say<-"Duplicated measurements IDs. Revise!"}
  measurements<-measurements[,1]
  measurements<-paste(measurements,".mzXML",sep="")
  if(any(match(measurements,filed,nomatch=0)==0)){say<-paste("Missing mzXML file - file corrupted? Compare project mzML folder for consistency!",sep="")}
  ##############################################################################
  # progress bar? ##############################################################
  if(interactive() && !.Platform$OS.type == "windows" && .Platform$GUI == "Rgui" && logfile[[5]][21]=="TRUE"){
	say<-"Disable the progress bar; works only under Windows OS"
  }
  ##############################################################################
  if(any(ls()=="logfile")){stop("\n illegal logfile detected #2 in check_project.r!")}
  return(say);
  
}
