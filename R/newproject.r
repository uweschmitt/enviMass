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
  if(grepl("\\",pro_dir,fixed=TRUE)){
	pro_dir<-gsub("\\",.Platform$file.sep,pro_dir,fixed=TRUE)
  }
  if(file.exists(file.path(pro_dir,pro_name))){
	return(FALSE);
  }
  if(!file.create(file.path(pro_dir,pro_name),showWarnings = FALSE)){
    return(FALSE); 
  }else{
	file.remove(file.path(pro_dir,pro_name))
  }
  ##############################################################################
  # create project #############################################################
  dir.create(file.path(pro_dir,pro_name),recursive=TRUE)              # project folder
  dir.create(file.path(pro_dir,pro_name,"files"),recursive=TRUE)      # subfolder
  dir.create(file.path(pro_dir,pro_name,"MSlist"),recursive=TRUE)     # subfolder
  dir.create(file.path(pro_dir,pro_name,"peaklist"),recursive=TRUE)   # subfolder  
  dir.create(file.path(pro_dir,pro_name,"features"),recursive=TRUE)   # subfolder
  dir.create(file.path(pro_dir,pro_name,"results"),recursive=TRUE)    # subfolder
  dir.create(file.path(pro_dir,pro_name,"dataframes"),recursive=TRUE) # subfolder
  dir.create(file.path(pro_dir,pro_name,"pics"),recursive=TRUE)       # subfolder
  dir.create(file.path(pro_dir,pro_name,"exports"),recursive=TRUE)    # subfolder
  ##############################################################################
  # write compound tables ######################################################  
  write.table(IS,file=file.path(pro_dir,pro_name,"dataframes","IS.txt"),row.names=FALSE,sep="\t",quote=FALSE)      	  	
  write.table(targets,file=file.path(pro_dir,pro_name,"dataframes","targets.txt"),row.names=FALSE,sep="\t",quote=FALSE)      	  
  # write measurement table #################################################### 
  measurements<-data.frame(c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),
    c("-"),c("FALSE"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"));
  names(measurements)<-c("ID","Name","Type","Mode","Place","Date","Time","include","copied?","picked?",
  "checked?","recal?","align?","norm?","feat?","comp?","IS_screen?","tar_screen?","tag1","tag2","tag3")
  write.csv(measurements,file=file.path(pro_dir,pro_name,"dataframes","measurements"),row.names=FALSE)
  ##############################################################################
  # create & save a logfile ####################################################
  logfile<-list(0);
    # folder name ##############################################################
    logfile[[1]]<-file.path(pro_dir,pro_name);
    names(logfile)[1]<-c("project_folder")
    # what MUST be done? ####################################################### 
	logfile[[2]]<-rep(FALSE,16);
	names(logfile[[2]])<-c(
		"peakpick","QC","recal","normalize","allign","profiling","trendblind","pattern",
		"screen_tar_file","screen_IS_file","screen_tar_comp","screen_IS_comp",
		"comp","homol","norm_prof","mass_defect"
	)	
    names(logfile)[2]<-c("Tasks_to_redo"); 
    # summary project status ###################################################
    tasks<-c(
		"Data available?",
		"Peak pick?",
		"QC?",
		"Isotope pattern?",
		"m/z recal.?",
		"Alligned?",
		"Intensity norm.?",
		"Profiled?",
		"IS norm.?",
		"Trend+Blind?",
		"",#"IS file-screen?",
		"",#"Target file-screen?",
		"",#"Compon?",
		"",#"IS comp-screen?",
		"",#"Target comp-screen?",
		"",#"Homologues?",
		""#"Mass defect?"
	 )
    doneit<-c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE)
    summar<-data.frame(tasks,doneit)
    names(summar)<-c("Tasks","Done?")
    logfile[[3]]<-summar
    names(logfile)[3]<-c("summary")
    # ProteoWizard MSConvert path ##############################################
    logfile[[4]]<-"C:\\Program Files\\ProteoWizard\\ProteoWizard 3.0.5140\\msconvert"
    names(logfile)[4]<-c("PW MSconvert path")
    # Parameters settings ######################################################
		# order of entries not a sequence! #####################################
		logfile[[5]]<-list(0)
		names(logfile)[5]<-c("parameters")
		# peak picking ###########################################################
		logfile$parameters[[1]]<-"1"; 	names(logfile$parameters)[1]<-"peak_MSlevel"     
		logfile$parameters[[2]]<-"300"; names(logfile$parameters)[2]<-"peak_drtgap"     
		logfile$parameters[[3]]<-"3.5"; names(logfile$parameters)[3]<-"peak_dmzdens"   	
		logfile$parameters[[4]]<-"4"; 	names(logfile$parameters)[4]<-"peak_minpeak"   	
		logfile$parameters$parameters<-"20"; 	names(logfile$parameters)[5]<-"peak_drtsmall2"   	
		logfile$parameters[[6]]<-"10"; 	names(logfile$parameters)[6]<-"peak_drtfill"  	
		logfile$parameters[[7]]<-"120"; names(logfile$parameters)[7]<-"peak_drtdens2"   
		logfile$parameters[[8]]<-"4"; 	names(logfile$parameters)[8]<-"peak_minint_log10"   
		logfile$parameters[[9]]<-"5"; 	names(logfile$parameters)[9]<-"peak_SN"  
		logfile$parameters[[10]]<-"2"; 	names(logfile$parameters)[10]<-"peak_SB"   
		logfile$parameters[[11]]<-"3"; 	names(logfile$parameters)[11]<-"peak_recurs"  
		logfile$parameters[[12]]<-"1"; 	names(logfile$parameters)[12]<-"peak_ended"  
		logfile$parameters[[13]]<-"1"; 	names(logfile$parameters)[13]<-"peak_weight"
		logfile$parameters[[14]]<-"6.5"; 	names(logfile$parameters)[14]<-"peak_maxint_log10"    
		# show progbar? ########################################################
		logfile$parameters[[21]]<-"FALSE";	names(logfile$parameters)[21]<-"progressBar"    
		# isotope patterns #####################################################
		logfile$parameters[[22]]<-"Elite/R240000@400";names(logfile$parameters)[22]<-"resolution"
		# recalibration ########################################################
		logfile$parameters[[30]]<-"Internal standards"; 	names(logfile$parameters)[30]<-"recal_use" 
		logfile$parameters[[31]]<-"3"; 						names(logfile$parameters)[31]<-"recal_dmz" 
		logfile$parameters[[32]]<-"TRUE"; 					names(logfile$parameters)[32]<-"recal_ppm" 
		logfile$parameters[[33]]<-"30"; 					names(logfile$parameters)[33]<-"recal_drt" 
		# trend detection ######################################################
		logfile$parameters[[29]]<-"TRUE";		names(logfile$parameters)[29]<-"notrend"
		logfile$parameters[[34]]<-"4,7,14"; 	names(logfile$parameters)[34]<-"trend_lags" 
		logfile$parameters[[35]]<-"3";			names(logfile$parameters)[35]<-"trend_threshold"
		# blind subtraction ####################################################
		logfile$parameters[[36]]<-"yes";		names(logfile$parameters)[36]<-"blind_do"	
		logfile$parameters[[37]]<-"100";		names(logfile$parameters)[37]<-"blind_threshold"	
		# profiling ############################################################
		logfile$parameters[[38]]<-"100";		names(logfile$parameters)[38]<-"prof_maxfiles"	
		logfile$parameters[[80]]<-"FALSE";		names(logfile$parameters)[80]<-"upto_file"
		logfile$parameters[[39]]<-"3";			names(logfile$parameters)[39]<-"prof_dmz"
		logfile$parameters[[40]]<-"TRUE";		names(logfile$parameters)[40]<-"prof_ppm"
		logfile$parameters[[41]]<-"60";			names(logfile$parameters)[41]<-"prof_drt"
		# IS screening #########################################################
		logfile$parameters[[42]]<-"30"; 		names(logfile$parameters)[42]<-"IS_drt1"	# RT tolerance of peaks in sample relative to their expected RT [s]
		logfile$parameters[[43]]<-"10"; 		names(logfile$parameters)[43]<-"IS_drt2"	# RT tolerance of peaks within an isotope pattern [s]
		logfile$parameters[[44]]<-"200"; 		names(logfile$parameters)[44]<-"IS_drt3"# RT tolerance of peaks in blank/blind relative to their expected RT [s]
		logfile$parameters[[45]]<-"3";			names(logfile$parameters)[45]<-"IS_dmz"# m/z tolerance ...
		logfile$parameters[[46]]<-"TRUE";		names(logfile$parameters)[46]<-"IS_ppm"# ... given in pppm?
		logfile$parameters[[47]]<-"30";			names(logfile$parameters)[47]<-"IS_inttol"# Intensity tolerance %
		logfile$parameters[[48]]<-"5E4";		names(logfile$parameters)[48]<-"IS_intcut"	# Lower intensity threhold
		logfile$parameters[[49]]<-"0.8";    	names(logfile$parameters)[49]<-"IS_w1" # Score weight for mass matching
		logfile$parameters[[50]]<-"0.2";		names(logfile$parameters)[50]<-"IS_w2"	# Score weight for relative intensity matching
		logfile$parameters[[51]]<-"0.0";		names(logfile$parameters)[51]<-"IS_w3"	# Score weight for occurrence in blank/blind		
		# target screening #####################################################
		logfile$parameters[[55]]<-"30"; 		names(logfile$parameters)[55]<-"tar_drt1"	# RT tolerance of peaks in sample relative to their expected RT [s]
		logfile$parameters[[56]]<-"10"; 		names(logfile$parameters)[56]<-"tar_drt2"	# RT tolerance of peaks within an isotope pattern [s]
		logfile$parameters[[57]]<-"200"; 		names(logfile$parameters)[57]<-"tar_drt3"# RT tolerance of peaks in blank/blind relative to their expected RT [s]
		logfile$parameters[[58]]<-"3";			names(logfile$parameters)[58]<-"tar_dmz"# m/z tolerance ...
		logfile$parameters[[59]]<-"TRUE";		names(logfile$parameters)[59]<-"tar_ppm"# ... given in pppm?
		logfile$parameters[[60]]<-"30";			names(logfile$parameters)[60]<-"tar_inttol"# Intensity tolerance %
		logfile$parameters[[61]]<-"5E4";		names(logfile$parameters)[61]<-"tar_intcut"	# Lower intensity threhold
		logfile$parameters[[62]]<-"0.8";    	names(logfile$parameters)[62]<-"tar_w1" # Score weight for mass matching
		logfile$parameters[[63]]<-"0.2";		names(logfile$parameters)[63]<-"tar_w2"	# Score weight for relative intensity matching
		logfile$parameters[[64]]<-"0.0";		names(logfile$parameters)[64]<-"tar_w3"	# Score weight for occurrence in blank/blind		
		# IS-based normalization ###############################################
		logfile$parameters[[70]]<-"90";		names(logfile$parameters)[70]<-"ISnorm_percfiles"# Minimum percentage of files covered by each IS profile %
		logfile$parameters[[71]]<-"15";		names(logfile$parameters)[71]<-"ISnorm_numbIS"# Minimum number of IS profiles
		logfile$parameters[[72]]<-"FALSE";	names(logfile$parameters)[72]<-"ISnorm_medblank"# Show median deviation of blank/blind profiles?
		logfile$parameters[[73]]<-"TRUE";	names(logfile$parameters)[73]<-"ISnorm_usesubblank"# Use subsampling
		logfile$parameters[[74]]<-"100";	names(logfile$parameters)[74]<-"ISnorm_numblank"	# Number of blank/blind profiles in subsample
		logfile$parameters[[75]]<-"FALSE";	names(logfile$parameters)[75]<-"ISnorm_medsam"# Show median deviation of sample (i.e., non-blank) profiles?
		logfile$parameters[[76]]<-"TRUE";	names(logfile$parameters)[76]<-"ISnorm_usesubsam"# Use subsampling
		logfile$parameters[[77]]<-"100";	names(logfile$parameters)[77]<-"ISnorm_numsam"	# Number of sample profiles in subsample
		logfile$parameters[[78]]<-".8";		names(logfile$parameters)[78]<-"ISnorm_score"# Screening threshold
		# Componentization #####################################################
		
 
		# Homologue series detection ###########################################

		  
    # Workflow settings ########################################################
    logfile$workflow<-0   
    names(logfile)[6]<-c("workflow")
		logfile$workflow[1]<-"yes"; 	names(logfile$workflow)[1]<-"qc" 
		logfile$workflow[2]<-"yes"; 	names(logfile$workflow)[2]<-"recal" 
		logfile$workflow[3]<-"yes"; 	names(logfile$workflow)[3]<-"align" 
		logfile$workflow[4]<-"yes"; 	names(logfile$workflow)[4]<-"norm" 
		logfile$workflow[5]<-"yes"; 	names(logfile$workflow)[5]<-"is_pattern" 
		logfile$workflow[6]<-"yes"; 	names(logfile$workflow)[6]<-"target_pattern" 
		logfile$workflow[7]<-"TRUE"; 	names(logfile$workflow)[7]<-"Comp_isotop" 
		logfile$workflow[8]<-"TRUE"; 	names(logfile$workflow)[8]<-"Comp_add" 
		logfile$workflow[9]<-"yes"; 	names(logfile$workflow)[9]<-"profiled" 
		logfile$workflow[10]<-"yes"; 	names(logfile$workflow)[10]<-"trenddetect"     
		logfile$workflow[11]<-"TRUE"; 	names(logfile$workflow)[11]<-"screen_IS_sam" 
		logfile$workflow[12]<-"TRUE"; 	names(logfile$workflow)[12]<-"screen_target_sam" 
		logfile$workflow[13]<-"TRUE"; 	names(logfile$workflow)[13]<-"screen_IS_comp" 
		logfile$workflow[14]<-"TRUE"; 	names(logfile$workflow)[14]<-"screen_target_comp" 		
		logfile$workflow[15]<-"yes"; 	names(logfile$workflow)[15]<-"profnorm" 
		logfile$workflow[16]<-"yes"; 	names(logfile$workflow)[16]<-"homol" 
		logfile$workflow[17]<-"yes"; 	names(logfile$workflow)[17]<-"massdef" 	  
		
		
		
    # positive adducts #########################################################
    logfile[[7]]<-0   
    names(logfile)[7]<-c("adducts_pos")
    logfile[[7]]<-"M+H";

    # negative adducts #########################################################
    logfile[[8]]<-0   
    names(logfile)[8]<-c("adducts_neg")
    logfile[[8]]<-"M-H";
    
    # isotopes #################################################################      
    logfile[[9]]<-0   
    names(logfile)[9]<-c("isotopes")
    logfile[[9]]<-"";
	# enviMass version number ##################################################
    logfile[[10]]<-2.2   
    names(logfile)[9]<-c("version")   
    # measurement data.frame ###################################################
	save(logfile,file=file.path(pro_dir,pro_name,"logfile.emp"));  
  ##############################################################################
  return(file.path(pro_dir,pro_name,"logfile.emp"));
}