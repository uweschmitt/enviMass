#' @title Start new enviMass workflow project
#'
#' @description \code{newproject} initializes an enviMass project
#'
#' @param pro_name Character string of project name.
#' @param pro_dir Character string. Path to a valid folder that will contains the project folder.
#' @param IS Dataframe containing internal standard (IS) compounds
#' @param targets Dataframe containing target and suspect compounds
#'
#' @details enviMass workflow function. Creates a project folder, various subfolders, a logfile, dummy outputs, etc
#' 

newproject<-function(pro_name,pro_dir,IS,targets){

  ##############################################################################
  # checks #####################################################################
  if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in newproject.r!")}
  if(grepl("\\",pro_dir,fixed=TRUE)){ # project directory established? 
	pro_dir<-gsub("\\",.Platform$file.sep,pro_dir,fixed=TRUE)
  }
  ##############################################################################
  # create project #############################################################
  if(!file.exists(file.path(pro_dir,pro_name))){ # should have been created in check_path
	dir.create(file.path(pro_dir,pro_name),recursive=TRUE)              # project folder
  }
  dir.create(file.path(pro_dir,pro_name,"files"),recursive=TRUE)      # subfolder
  dir.create(file.path(pro_dir,pro_name,"MSlist"),recursive=TRUE)     # subfolder
  dir.create(file.path(pro_dir,pro_name,"peaklist"),recursive=TRUE)   # subfolder  
  dir.create(file.path(pro_dir,pro_name,"features"),recursive=TRUE)   # subfolder
  dir.create(file.path(pro_dir,pro_name,"results"),recursive=TRUE)    # subfolder
	dir.create(file.path(pro_dir,pro_name,"results","screening"),recursive=TRUE)    	# subfolder  
	dir.create(file.path(pro_dir,pro_name,"results","LOD"),recursive=TRUE)    	# subfolder  
	dir.create(file.path(pro_dir,pro_name,"results","recalibration"),recursive=TRUE)   # subfolder 
	dir.create(file.path(pro_dir,pro_name,"results","componentization"),recursive=TRUE)   # subfolder 
	dir.create(file.path(pro_dir,pro_name,"results","componentization","adducts"),recursive=TRUE)   # subfolder 
	dir.create(file.path(pro_dir,pro_name,"results","componentization","isotopologues"),recursive=TRUE)   # subfolder 	
	dir.create(file.path(pro_dir,pro_name,"results","componentization","EIC_corr"),recursive=TRUE)   # subfolder 
	dir.create(file.path(pro_dir,pro_name,"results","componentization","homologues"),recursive=TRUE)   # subfolder 	
  dir.create(file.path(pro_dir,pro_name,"dataframes"),recursive=TRUE) # subfolder
  dir.create(file.path(pro_dir,pro_name,"pics"),recursive=TRUE)       # subfolder
  dir.create(file.path(pro_dir,pro_name,"exports"),recursive=TRUE)    # subfolder
  dir.create(file.path(pro_dir,pro_name,"quantification"),recursive=TRUE)    # subfolder  
  ##############################################################################
  # write compound tables ######################################################  
  write.table(IS,file=file.path(pro_dir,pro_name,"dataframes","IS.txt"),row.names=FALSE,sep="\t",quote=FALSE)      	  	
  write.table(targets,file=file.path(pro_dir,pro_name,"dataframes","targets.txt"),row.names=FALSE,sep="\t",quote=FALSE)      	  
  # write measurement table #################################################### 
  measurements<-data.frame(c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),
    c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),
	c("-"),c("-"),c("-"),c("-"),c("-"),c("-"));
  names(measurements)<-c("ID","Name","Type","Mode","Place","Date","Time","include","copied","peakpicking",
  "checked","recal","align","norm","profiled","LOD","IS_screen","tar_screen","tag1","tag2","tag3","date_end","time_end",
  "isotopologues","adducts","homologues","EIC_correlation","blind","ID_2")
  write.csv(measurements,file=file.path(pro_dir,pro_name,"dataframes","measurements"),row.names=FALSE)
  ##############################################################################
  # create & save a logfile ####################################################
	workflow_depend<-read.table(
		file="workflow_depend"		
	)
	workflow_depend<-as.matrix(workflow_depend)
	workflow_must<-read.table(
		file="workflow_must"			
	)
	workflow_must<-as.matrix(workflow_must)
	############################################################################
	# some checks on matrices ##################################################
	if(length(colnames(workflow_depend))!=length(rownames(workflow_depend))){
		stop("\n Invalid matrix workflow_depend: matrix not quadratic; please revise!")
	}
	if(length(colnames(workflow_must))!=length(rownames(workflow_must))){
		stop("\n Invalid matrix workflow_must: matrix not quadratic; please revise!")
	}	
	if(any(is.na(match(colnames(workflow_depend),rownames(workflow_depend))))){
		stop("\n Invalid matrix workflow_depend: unequal row vs. column names!")	
	}
	if(any(is.na(match(colnames(workflow_must),rownames(workflow_must))))){
		stop("\n Invalid matrix workflow_must: unequal row vs. column names!")	
	}	
	if(any(duplicated(colnames(workflow_must)))){
		stop("\n Invalid matrix workflow_must: duplicated node names!")		
	}	
	if(any(duplicated(colnames(workflow_depend)))){
		stop("\n Invalid matrix workflow_depend: duplicated node names!")		
	}
	if(any(is.na(match(colnames(workflow_depend),colnames(workflow_must))))){
		stop("\n Different nodes in workflow_depend vs. workflow_must: revise!")	
	}	
	if(any(is.na(match(colnames(workflow_must),colnames(workflow_depend))))){
		stop("\n Different nodes in workflow_depend vs. workflow_must: revise!")	
	}	
	############################################################################
	# some checks on script availability #######################################
	files<-list.files()
	for(i in 1:length(colnames(workflow_depend))){
		if(!any(files==paste("do_",colnames(workflow_depend)[i],".r",sep=""))){
			stop(paste("Missing do_ scripts.r for node ",colnames(workflow_depend)[i],sep=""))
		}
		if(!any(files==paste("dont_",colnames(workflow_depend)[i],".r",sep=""))){
			stop(paste("Missing dont_ scripts.r for node ",colnames(workflow_depend)[i],sep=""))
		}
	}
	############################################################################
	logfile<-list(0);
    # folder name ##############################################################
    logfile[[1]]<-file.path(pro_dir,pro_name);
    names(logfile)[1]<-c("project_folder")
    # nodes - what MUST be (re)done? ###########################################
	logfile[[2]]<-rep(FALSE,length(colnames(workflow_must)));
	names(logfile[[2]])<-colnames(workflow_must)
	names(logfile)[2]<-c("Tasks_to_redo");    
	# summary project status ###################################################
    tasks<-names(logfile[[2]]) # based on above Tasks_to_redo
    doneit<-rep(FALSE,length(tasks))
    summar<-data.frame(tasks,doneit,stringsAsFactors = FALSE)
    names(summar)<-c("Tasks","Done?")
    logfile[[3]]<-summar
    names(logfile)[3]<-c("summary")
    # ProteoWizard MSConvert path ##############################################
    logfile[[4]]<-"C:/Program Files/ProteoWizard/ProteoWizard 3.0.5140/msconvert.exe"
    names(logfile)[4]<-c("PW MSconvert path")
    # Parameters settings ######################################################
		# order of entries not a sequence! #####################################
		logfile[[5]]<-list(0)
		names(logfile)[5]<-c("parameters")
		# peak picking ###########################################################
		logfile$parameters$peak_MSlevel<-"1"; 	
		logfile$parameters$peak_drtgap<-"300"; 
		logfile$parameters$peak_dmzdens<-"3.5"; 
		logfile$parameters$peak_minpeak<-"40"; 		
		logfile$parameters$peak_drtsmall2<-"20"; 		
		logfile$parameters$peak_drtfill<-"10"; 		
		logfile$parameters$peak_drtdens2<-"120";
		logfile$parameters$peak_minint_log10<-"4"; 	 
		logfile$parameters$peak_SN<-"5"; 	
		logfile$parameters$peak_SB<-"2"; 	
		logfile$parameters$peak_recurs<-"3"; 
		logfile$parameters$peak_ended<-"1"; 	
		logfile$parameters$peak_weight<-"1"; 	
		logfile$parameters$peak_maxint_log10<-"6.5"; 	   
		logfile$parameters$peak_perc_cut<-"0"; 	
		logfile$parameters$peak_which_intensity<-"maximum"
		# show progbar? ########################################################
		logfile$parameters$progressBar<-"FALSE";	 
		# isotope patterns #####################################################
		logfile$parameters$resolution<-"Elite_R240000@400";
		# recalibration ########################################################
		logfile$parameters$recal_use<-"Internal standards"; 	
		logfile$parameters$recal_dmz<-"3"; 					
		logfile$parameters$recal_ppm<-"TRUE"; 				
		logfile$parameters$recal_drt<-"30"; 					
		logfile$parameters$recal_maxdmz<-"30";						
		# replicate intersection ################################################
		logfile$parameters$replicate_dmz<-"3";						
		logfile$parameters$replicate_ppm<-"TRUE";						
		logfile$parameters$replicate_recalib<-"FALSE";					
		logfile$parameters$replicate_delRT<-"30";					
		logfile$parameters$replicate_IS_dInt<-"5";						
		# trend detection ######################################################
		logfile$parameters$notrend<-"TRUE";		
		logfile$parameters$trend_lags<-"4,7,14"; 	
		logfile$parameters$trend_threshold<-"3";			
		# blind subtraction ####################################################	
		logfile$parameters$trend_blind<-"yes";				
		logfile$parameters$blind_threshold<-"100";			
		logfile$parameters$blind_dmz<-"3";			
		logfile$parameters$blind_ppm<-"TRUE";			
		logfile$parameters$blind_drt<-"30";				
		logfile$parameters$subtract_pos_bydate<-"FALSE";		
		logfile$parameters$subtract_pos_byfile<-"FALSE";		
		logfile$parameters$subtract_neg_bydate<-"FALSE";		
		logfile$parameters$subtract_neg_byfile<-"FALSE";		
		logfile$parameters$blind_omit<-"no";			
		# profiling ############################################################
		logfile$parameters$prof_maxfiles<-"100";		
		logfile$parameters$upto_file<-"FALSE";		
		logfile$parameters$prof_dmz<-"3";		
		logfile$parameters$prof_ppm<-"TRUE";		
		logfile$parameters$prof_drt<-"60";			
		logfile$parameters$prof_select<-"FALSE";		
		logfile$parameters$replicates_prof<-"yes";		
		# IS screening #########################################################
		logfile$parameters$IS_drt1<-"30"; 			# RT tolerance of peaks in sample relative to their expected RT [s]
		logfile$parameters$IS_drt2<-"10"; 			# RT tolerance of peaks within an isotope pattern [s]
		logfile$parameters$IS_dmz<-"3";				# m/z tolerance ...
		logfile$parameters$IS_ppm<-"TRUE";			# ... given in pppm?
		logfile$parameters$IS_inttol<-"30";			# Intensity tolerance %
		logfile$parameters$IS_intcut<-"50000";		# Lower intensity threhold
		logfile$parameters$IS_w1<-"0.8";    		# Matching score
		logfile$parameters$screen_IS_cutit<-"FALSE";    	# Cut off match combiantions below matching score?	
		logfile$parameters$screen_IS_maxonly<-"FALSE";    	# Screen only most intense isotopologue peak?		
		# target screening #####################################################
		logfile$parameters$tar_drt1<-"30"; 		# RT tolerance of peaks in sample relative to their expected RT [s]
		logfile$parameters$tar_drt2<-"10"; 		# RT tolerance of peaks within an isotope pattern [s]
		logfile$parameters$tar_dmz<-"3";		# m/z tolerance ...
		logfile$parameters$tar_ppm<-"TRUE";		# ... given in pppm?
		logfile$parameters$tar_inttol<-"30";	# Intensity tolerance %
		logfile$parameters$tar_intcut<-"50000";	# Lower intensity threhold
		logfile$parameters$tar_w1<-"0.8";    	# Matching score	
		logfile$parameters$screen_target_cutit<-"FALSE";    	# Cut off match combiantions below matching score?		
		logfile$parameters$screen_target_maxonly<-"FALSE";    	# Screen only most intense isotopologue peak?			
		# IS-based normalization ###############################################
		logfile$parameters$ISnorm_percfiles<-"90";		# Minimum percentage of files covered by each IS profile %
		logfile$parameters$ISnorm_numbIS<-"15";			# Minimum number of IS profiles
		logfile$parameters$ISnorm_medblank<-"FALSE";	# Show median deviation of blank/blind profiles?
		logfile$parameters$ISnorm_usesubblank<-"TRUE";	# Use subsampling
		logfile$parameters$ISnorm_numblank<-"100";		# Number of blank/blind profiles in subsample
		logfile$parameters$ISnorm_medsam<-"FALSE";		# Show median deviation of sample (i.e., non-blank) profiles?
		logfile$parameters$ISnorm_usesubsam<-"TRUE";	# Use subsampling
		logfile$parameters$ISnorm_numsam<-"100";		# Number of sample profiles in subsample
		logfile$parameters$ISnorm_score<-"0.8";			# Screening threshold
		# subtraction ##########################################################
		logfile$parameters$subtr_IS<-"yes"; 		
		logfile$parameters$subtr_target<-"yes"; 		
		logfile$parameters$subtr_blind<-"yes"; 			
		logfile$parameters$subtr_spiked<-"yes"; 		
		# add custom parameters ################################################
		source(file="workflow_parameters.r",local=TRUE)
	if(any(duplicated(names(logfile$parameters)))){stop("Duplicated parameter names found - revise!")}	
	# Workflow settings ########################################################
    logfile$workflow<-0    # based on above Tasks_to_redo
    names(logfile)[6]<-c("workflow")
	for(i in 1:length(names(logfile[[2]]))){
		# use simple initial workflow settings
		if(any(names(logfile[[2]])[i]==c("peakpicking","LOD","profiling","IS_screen","target_screen"))){
			logfile$workflow[i]<-"yes"; 
		}else{
			logfile$workflow[i]<-"no"; 		
		}
		names(logfile$workflow)[i]<-names(logfile[[2]])[i]
	}
	################################################################################################
	# define matrix of downstream workflow dependencies (==1) and ##################################
	# recalculations of previous steps if their results are overwritten (==2), e.g. IS_subtr or ####
	# target_subtr or target screening result tables/lists for quantification ######################
	# requires only a definition of direct ones - indirect ones will be in workflow_set.r ##########
	# dependencies must simply go after their parent node ########################################## 
	# order here actually irrelevant, because calculation order set with workflow_schedule #########
	logfile[[11]]<-workflow_depend
	names(logfile)[11]<-"workflow_depend"	
	################################################################################################
	# define upstream workflow "musts", i.e., upstream nodes on which`s execution a node ###########
	# depends. 0 = not dependent. 1 = dependent. -1 = MUST NOT be executed (not yet further implemented)  	
	logfile[[12]]<-workflow_must
	names(logfile)[12]<-"workflow_must"	
	################################################################################################	
	# reorder summary into workflow ################################################################
	schedule<-enviMass:::workflow_schedule(logfile$workflow_depend,logfile$workflow_must)
	if(!is.data.frame(schedule)){stop(schedule)}
	set_order<-match(schedule[,1],logfile$summary[,1])
	logfile$summary<-logfile$summary[set_order,]	
	################################################################################################
    # positive adducts #########################################################
    logfile[[7]]<-0   
    names(logfile)[7]<-c("adducts_pos")
    logfile[[7]]<-"M+H";
    # negative adducts #########################################################
    logfile[[8]]<-0   
    names(logfile)[8]<-c("adducts_neg")
    logfile[[8]]<-"M-H";
    # isotopes #################################################################      
    logfile[[9]]<-"";
    names(logfile)[9]<-c("isotopes")
	# enviMass version number ##################################################
    logfile[[10]]<-3.117
    names(logfile)[10]<-c("version")   
	# subtraction files ########################################################
	logfile[[13]]<-"FALSE"
	names(logfile)[13]<-"Positive_subtraction_files"
	logfile[[14]]<-"FALSE"
	names(logfile)[14]<-"Negative_subtraction_files"
	# calibration model ########################################################
	cal_models_pos<-list()
	save(cal_models_pos,file=file.path(logfile$project_folder,"quantification","cal_models_pos"));	
	cal_models_neg<-list()
	save(cal_models_neg,file=file.path(logfile$project_folder,"quantification","cal_models_neg"));	
    # measurement data.frame ###################################################
	save(logfile,file=file.path(pro_dir,pro_name,"logfile.emp"));  
	rm(logfile)
	##############################################################################
	return(file.path(pro_dir,pro_name,"logfile.emp"));
  
}
