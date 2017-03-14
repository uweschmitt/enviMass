#' @title Check enviMass project
#'
#' @description \code{} runs several cconsistency checks on a enviMass project
#'
#' @param intstand_check 
#' @param
#' @param
#' @param
#' @param
#' @param
#' @param
#'
#' @return Character string with either specific error message or a "Project consistent" message.
#' 
#' @details enviMass workflow function; run before further calculations are started in the workflow.
#' 

check_compounds<-function(intstand_check,targets_check,isotopes,adducts,logfile,write_tables=FALSE,check_relat=TRUE){

	say<-"Project consistent"
	##############################################################################
	# compounds available & ok? ##################################################
	if(length(intstand_check$Formula[intstand_check$Formula!="-"])==0){
		cat("No internal standards available \n") 
	}else{
		if(
			!all(names(intstand_check)==c("ID","Name","Formula","RT","RT_tolerance","main_adduct","ion_mode",
			"use_for_recalibration","use_for_screening", "restrict_adduct","Remark","tag1","tag2","tag3",
			"from","to","Lower_intensity_bound","Upper_intensity_bound","Quant_adduct","Quant_peak","Quant_rule"))
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
		# check columns ###########################################################
		# ID with underscore? #####################################################
		if(any(grepl("_",intstand_check[,"ID"]))){
			say<-"Please revise: do not use underscores in internal standard IDs."
		}
		# main adduct correct? ####################################################
		if(any(!charmatch(intstand_check[intstand_check[,"ion_mode"]=="positive","main_adduct"],c(as.character(adducts[adducts[,6]=="positive",1]),"FALSE"),nomatch=FALSE))){
			wrongadduct<-intstand_check[intstand_check[,"ion_mode"]=="positive","main_adduct"][!charmatch(intstand_check[intstand_check[,"ion_mode"]=="positive","main_adduct"],c(as.character(adducts[adducts[,6]=="positive",1]),"FALSE"),nomatch=FALSE)]
			wrongadduct<-wrongadduct[1]
			say<-paste("IS main_adduct (column 6, positive ion_mode): wrong main_adduct:",wrongadduct,". Pleas revise!")
		}
		if(any(!charmatch(intstand_check[intstand_check[,"ion_mode"]=="negative","main_adduct"],c(as.character(adducts[adducts[,6]=="negative",1]),"FALSE"),nomatch=FALSE))){
			wrongadduct<-intstand_check[intstand_check[,"ion_mode"]=="negative","main_adduct"][!charmatch(intstand_check[intstand_check[,"ion_mode"]=="negative","main_adduct"],c(as.character(adducts[adducts[,6]=="negative",1]),"FALSE"),nomatch=FALSE)]
			wrongadduct<-wrongadduct[1]
			say<-paste("IS main_adduct (column 6, negative ion_mode): wrong main_adduct:",wrongadduct,". Pleas revise!")
		}
		# ionization mode correct? ################################################
		if(any(!charmatch(intstand_check[,"ion_mode"],c("positive","negative"),nomatch=FALSE))){
			say<-"Invalid IS ion_mode (column 7) entry (must be positive or negative)"
		}
		# RT numeric? #############################################################
		got<-which(is.na(as.numeric(intstand_check[,"RT"])))
		if(length(got)>0){
			say<-paste("Invalid RT (column 4) for internal standards with ID ",intstand_check[got,"ID"],". Please revise!",sep="")
		}
		# RT tolerance numeric? ###################################################
		got<-which(is.na(as.numeric(intstand_check[intstand_check[,"RT_tolerance"]!="FALSE","RT_tolerance"])))
		if(length(got)>0){
			say<-paste("Invalid RT_tolerance (column 5) for internal standards with ID ",intstand_check[got,"ID"],". must be either FALSE or numeric. Please revise!",sep="")
		}
		# use_for_recalibration correct? ########################################## 
		got<-which((intstand_check[,"use_for_recalibration"]!="FALSE")&(intstand_check[,"use_for_recalibration"]!="TRUE"))
		if(length(got)>0){
			say<-paste("Invalid use_for_recalibration for internal standards with ID ",intstand_check[got,"ID"],". must be either FALSE or TRUE. Please revise!",sep="")
		}		
		# use_for_screening correct? ################################names(##############
		got<-which((intstand_check[,"use_for_screening"]!="FALSE")&(intstand_check[,"use_for_screening"]!="TRUE"))
		if(length(got)>0){
			say<-paste("Invalid use_for_screening for internal standards with ID ",intstand_check[got,"ID"],". must be either FALSE or TRUE. Please revise!",sep="")
		}		
		# restrict_adduct correct? ################################################
		got<-which((intstand_check[,"restrict_adduct"]!="FALSE")&(intstand_check[,"restrict_adduct"]!="TRUE"))
		if(length(got)>0){
			say<-paste("Invalid restrict_adduct for internal standards with ID ",intstand_check[got,"ID"],". must be either FALSE or TRUE. Please revise!",sep="")
		}		
		# check from 
		# check to
		# check Lower_intensity_bound ##############################################
		got<-which(is.na(as.numeric(intstand_check[,"Lower_intensity_bound"])))
		if(length(got)>0){
			say<-paste("Invalid Lower_intensity_bound for internal standards with ID ",intstand_check[got,"ID"],". Must be numeric (i.e., 0 if not specified). Please revise!",sep="")
		}
		# check Upper_intensity_bound ##############################################
		got<-which(is.na(as.numeric(intstand_check[,"Upper_intensity_bound"])))
		if(length(got)>0){
			say<-paste("Invalid Upper_intensity_bound for internal standards with ID ",intstand_check[got,"ID"],". Must be numeric or Inf if not specified. Please revise!",sep="")
		}
		# quantization adduct correct? #############################################
		if(any(!charmatch(intstand_check[intstand_check[,"ion_mode"]=="positive","Quant_adduct"],c(as.character(adducts[adducts[,6]=="positive",1]),"FALSE"),nomatch=FALSE))){
			wrongadduct<-intstand_check[intstand_check[,"ion_mode"]=="positive","Quant_adduct"][!charmatch(intstand_check[intstand_check[,"ion_mode"]=="positive","Quant_adduct"],c(as.character(adducts[adducts[,6]=="positive",1]),"FALSE"),nomatch=FALSE)]
			wrongadduct<-wrongadduct[1]
			say<-paste("IS main_adduct (column 19, positive ion_mode): wrong Quant_adduct:",wrongadduct,". Pleas revise!")
		}
		if(any(!charmatch(intstand_check[intstand_check[,"ion_mode"]=="negative","Quant_adduct"],c(as.character(adducts[adducts[,6]=="negative",1]),"FALSE"),nomatch=FALSE))){
			wrongadduct<-intstand_check[intstand_check[,"ion_mode"]=="negative","Quant_adduct"][!charmatch(intstand_check[intstand_check[,"ion_mode"]=="negative","Quant_adduct"],c(as.character(adducts[adducts[,6]=="negative",1]),"FALSE"),nomatch=FALSE)]
			wrongadduct<-wrongadduct[1]
			say<-paste("IS main_adduct (column 19, negative ion_mode): wrong Quant_adduct:",wrongadduct,". Pleas revise!")
		}
		# check Quant_peak ########################################################
		got<-which(is.na(as.numeric(intstand_check[intstand_check[,"Quant_peak"]!="FALSE","Quant_peak"])))
		if(length(got)>0){
			say<-paste("Invalid Quant_peak for internal standards with ID ",intstand_check[got,"ID"],". must be either FALSE or an integer. Please revise!",sep="")
		}		
		# check Quant_rule ########################################################
		got<-which(any(!charmatch(intstand_check[,"Quant_rule"],c("most intense peak","closest RT","closest m/z","FALSE"))))
		if(length(got)>0){		
			say<-paste("Invalid Quant_rule for internal standards with ID ",intstand_check[got,"ID"],". must be either most 'intense peak', 'closest RT' or 'closest m/z'. Please revise!",sep="")
		}
		###########################################################################
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
			say<-paste("Invalid molecular formula detected for internal standard with ID:",intstand_check[checked[,1]==TRUE,"ID"])
		}else{ # save corrected formulas!
			if(write_tables){	
				intstand_check[,3]<-checked[,2]
				write.table(intstand_check,file=file.path(logfile[[1]],"dataframes","IS.txt"),row.names=FALSE,sep="\t",quote=FALSE)
			}
		}
		rm(checked)
	}
	if(length(targets_check$Formula[targets_check$Formula!="-"])==0){
		  cat("No targets available \n")    
	}else{
		if(
			!all(names(targets_check)==c("ID","Name","Formula","RT","RT_tolerance","ID_internal_standard","main_adduct","ion_mode",
			"use_for_recalibration","use_for_screening","restrict_adduct","Remark","tag1","tag2","tag3","from","to","warn_1","warn_2",
			"Quant_adduct","Quant_peak","Quant_rule"))
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
		# check columns ###########################################################
		# ID with underscore? #####################################################
		# check columns ###########################################################
		# ID with underscore? #####################################################
		if(any(grepl("_",targets_check[,"ID"]))){
			say<-"Please revise: do not use underscores for target IDs."
		}
		# main adduct correct? ####################################################
		if(any(!charmatch(targets_check[targets_check[,"ion_mode"]=="positive","main_adduct"],c(as.character(adducts[adducts[,6]=="positive",1]),"FALSE"),nomatch=FALSE))){
			wrongadduct<-targets_check[targets_check[,"ion_mode"]=="positive","main_adduct"][!charmatch(targets_check[targets_check[,"ion_mode"]=="positive","main_adduct"],c(as.character(adducts[adducts[,6]=="positive",1]),"FALSE"),nomatch=FALSE)]
			wrongadduct<-wrongadduct[1]
			say<-paste("target main_adduct (column 7, positive ion_mode): wrong main_adduct:",wrongadduct,". Pleas revise!")
		}
		if(any(!charmatch(targets_check[targets_check[,"ion_mode"]=="negative","main_adduct"],c(as.character(adducts[adducts[,6]=="negative",1]),"FALSE"),nomatch=FALSE))){
			wrongadduct<-targets_check[targets_check[,"ion_mode"]=="negative","main_adduct"][!charmatch(targets_check[targets_check[,"ion_mode"]=="negative","main_adduct"],c(as.character(adducts[adducts[,6]=="negative",1]),"FALSE"),nomatch=FALSE)]
			wrongadduct<-wrongadduct[1]
			say<-paste("target main_adduct (column 7, negative ion_mode): wrong main_adduct:",wrongadduct,". Pleas revise!")
		}
		# ionization mode correct? ################################################
		if(any(!charmatch(targets_check[,"ion_mode"],c("positive","negative"),nomatch=FALSE))){
			say<-"Invalid target ion_mode (column 8) entry (must be positive or negative)"
		}
		# RT numeric? #############################################################
		got<-which(is.na(as.numeric(targets_check[,"RT"])))
		if(length(got)>0){
			say<-paste("Invalid RT (column 4) for target with ID ",targets_check[got,"ID"],". Please revise!",sep="")
		}
		# RT tolerance numeric? ###################################################
		got<-which(is.na(as.numeric(targets_check[targets_check[,"RT_tolerance"]!="FALSE","RT_tolerance"])))
		if(length(got)>0){
			say<-paste("Invalid RT_tolerance (column 5) for target with ID ",targets_check[got,"ID"],". must be either FALSE or numeric. Please revise!",sep="")
		}
		# use_for_recalibration correct? ########################################## 
		got<-which((targets_check[,"use_for_recalibration"]!="FALSE")&(targets_check[,"use_for_recalibration"]!="TRUE"))
		if(length(got)>0){
			say<-paste("Invalid use_for_recalibration for target with ID ",targets_check[got,"ID"],". must be either FALSE or TRUE. Please revise!",sep="")
		}		
		# use_for_screening correct? ################################names(##############
		got<-which((targets_check[,"use_for_screening"]!="FALSE")&(targets_check[,"use_for_screening"]!="TRUE"))
		if(length(got)>0){
			say<-paste("Invalid use_for_screening for target with ID ",targets_check[got,"ID"],". must be either FALSE or TRUE. Please revise!",sep="")
		}		
		# restrict_adduct correct? ################################################
		got<-which((targets_check[,"restrict_adduct"]!="FALSE")&(targets_check[,"restrict_adduct"]!="TRUE"))
		if(length(got)>0){
			say<-paste("Invalid restrict_adduct for target with ID ",targets_check[got,"ID"],". must be either FALSE or TRUE. Please revise!",sep="")
		}		
		# check from 
		# check to
		# quantization adduct correct? #############################################
		if(any(!charmatch(targets_check[targets_check[,"ion_mode"]=="positive","Quant_adduct"],c(as.character(adducts[adducts[,6]=="positive",1]),"FALSE"),nomatch=FALSE))){
			wrongadduct<-targets_check[targets_check[,"ion_mode"]=="positive","Quant_adduct"][!charmatch(targets_check[targets_check[,"ion_mode"]=="positive","Quant_adduct"],c(as.character(adducts[adducts[,6]=="positive",1]),"FALSE"),nomatch=FALSE)]
			wrongadduct<-wrongadduct[1]
			say<-paste("target main_adduct (column 20, positive ion_mode): wrong Quant_adduct:",wrongadduct,". Pleas revise!")
		}
		if(any(!charmatch(targets_check[targets_check[,"ion_mode"]=="negative","Quant_adduct"],c(as.character(adducts[adducts[,6]=="negative",1]),"FALSE"),nomatch=FALSE))){
			wrongadduct<-targets_check[targets_check[,"ion_mode"]=="negative","Quant_adduct"][!charmatch(targets_check[targets_check[,"ion_mode"]=="negative","Quant_adduct"],c(as.character(adducts[adducts[,6]=="negative",1]),"FALSE"),nomatch=FALSE)]
			wrongadduct<-wrongadduct[1]
			say<-paste("target main_adduct (column 20, negative ion_mode): wrong Quant_adduct:",wrongadduct,". Pleas revise!")
		}
		# check Quant_peak ########################################################
		got<-which(is.na(as.numeric(targets_check[targets_check[,"Quant_peak"]!="FALSE","Quant_peak"])))
		if(length(got)>0){
			say<-paste("Invalid Quant_peak for internal standards with ID ",targets_check[got,"ID"],". must be either FALSE or an integer. Please revise!",sep="")
		}		
		# check Quant_rule ########################################################
		got<-which(any(!charmatch(targets_check[,"Quant_rule"],c("most intense peak","closest RT","closest m/z","FALSE"))))
		if(length(got)>0){		
			say<-paste("Invalid Quant_rule for internal standards with ID ",targets_check[got,"ID"],". must be either most 'intense peak', 'closest RT' or 'closest m/z'. Please revise!",sep="")
		}
		# check concentration warn levels #########################################
		if(logfile$workflow[names(logfile$workflow)=="quantification"]=="yes"){
			this<-which(is.na(as.numeric(targets_check$warn_1[targets_check$warn_1!="FALSE"])))
			if(length(this)>0){
				say<-paste(
					"Invalid entries for concentration warn level #1 found in target table rows #",paste(as.character(this),collapse=", "),
					". Entries must be FALSE or numbers, please revise.",sep=""
				)	
			}
			this<-which(is.na(as.numeric(targets_check$warn_2[targets_check$warn_2!="FALSE"])))
			if(length(this)>0){
				say<-paste(
					"Invalid entries for concentration warn level #2 found in target table rows #",paste(as.character(this),collapse=", "),
					". Entries must be FALSE or numbers, please revise.",sep=""
				)	
			}
		} 
		###########################################################################
		# formulas make sense? ####################################################
		checked<-enviPat::check_chemform(isotopes, targets_check[,3])
		if(any(checked[,1])){
			say<-paste("Invalid molecular formula detected for target with ID:",targets_check[checked[,1]==TRUE,"ID"])
		}else{
			if(write_tables){
				targets_check[,3]<-checked[,2]
				write.table(targets_check,file=file.path(logfile[[1]],"dataframes","targets.txt"),row.names=FALSE,sep="\t",quote=FALSE)       
			}			
		}
		rm(checked) 		
	}   
	##############################################################################
	# check compounds and files for calibration & quantification #################
	if(
		(
			(logfile$workflow[names(logfile$workflow)=="calibration"])=="yes" ||
			(logfile$workflow[names(logfile$workflow)=="quantification"])=="yes"
		)&(check_relat)
	){
		# check if all relations / adducts are correct
		for(i in 1:length(targets_check[,"ID"])){
			if(targets_check[i,"ID_internal_standard"]!="FALSE"){ # ISTD entry triggers check - does it also trigger the quantification?
				ID_target_missing<-targets_check[i,"ID"];
				ID_IS_missing<-targets_check[i,"ID_internal_standard"];
				target_quan_adduct<-targets_check[i,names(targets_check)=="Quant_adduct"]	
				target_avail_adduct<-targets_check[i,names(targets_check)=="main_adduct"]
				if(targets_check[i,names(targets_check)=="restrict_adduct"]=="FALSE"){
					if(targets_check[i,"ion_mode"]=="positive"){
						target_avail_adduct<-c(target_avail_adduct,logfile[["adducts_pos"]])
					}else{
						target_avail_adduct<-c(target_avail_adduct,logfile[["adducts_neg"]])
					}
				}
				# does the target quantification adduct exist?
				if( (!any(target_avail_adduct==target_quan_adduct)) | (target_quan_adduct=="FALSE") ){ 
					say<-paste(
					"Adduct used in quantification for target compound with ID ",
					ID_target_missing," not valid. Either it is set to FALSE and shouldn`t or the required adduct is not included in the adduct selection of this target! Please adjust.",sep="")
				}
				# does the internal standard exist?
				found_standard<-TRUE
				if( !any(intstand_check[ intstand_check[,"ion_mode"]==targets_check[i,"ion_mode"],1]==targets_check[i,"ID_internal_standard"]) ){ 		
					say<-paste(
					"Quantification/calibration problem: internal standard with ID ",
					ID_IS_missing," for target compound with ID ",ID_target_missing,
					" not found at this ionization mode. Check compound target table, column ID_internal_standard.",
					sep="")
					found_standard<-FALSE
				}
				# does the standard quantification adduct exist?
				if(found_standard){ 
					j<-(intstand_check[,1]==ID_IS_missing)
					IS_quan_adduct<-intstand_check[j,names(intstand_check)=="Quant_adduct"]				
					IS_avail_adduct<-intstand_check[j,names(intstand_check)=="main_adduct"]	
					if(intstand_check[j,names(intstand_check)=="restrict_adduct"]=="FALSE"){
						IS_avail_adduct<-c(IS_avail_adduct,logfile[[7]],logfile[[8]])
					}
					# does the target quantification adduct exist?
					if( (!any(IS_avail_adduct==IS_quan_adduct)) | (IS_quan_adduct=="FALSE") ){ 
						say<-paste(
						"Adduct used in quantification for internal standard with ID ",
						ID_IS_missing," not valid. Either it is set to FALSE and shouldn`t or the required adduct is not included in the adduct selection of this internal standard! Please adjust.",sep="")
					}
				}
			}
		}	
	}
	##############################################################################
	return(say);
	##############################################################################  
	
}
