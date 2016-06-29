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
	dir.create(file.path(pro_dir,pro_name,"results","screening"),recursive=TRUE)    	# subfolder  
	dir.create(file.path(pro_dir,pro_name,"results","LOD"),recursive=TRUE)    	# subfolder  
	dir.create(file.path(pro_dir,pro_name,"results","recalibration"),recursive=TRUE)   # subfolder 
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
    c("-"),c("FALSE"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"));
  names(measurements)<-c("ID","Name","Type","Mode","Place","Date","Time","include","copied","picked",
  "checked","recal","align","norm","profiled","LOD","IS_screen","tar_screen","tag1","tag2","tag3","date_end","time_end")
  write.csv(measurements,file=file.path(pro_dir,pro_name,"dataframes","measurements"),row.names=FALSE)
  ##############################################################################
  # create & save a logfile ####################################################
  logfile<-list(0);
    # folder name ##############################################################
    logfile[[1]]<-file.path(pro_dir,pro_name);
    names(logfile)[1]<-c("project_folder")
    # what MUST be done? ####################################################### 
	logfile[[2]]<-rep(FALSE,22);
	names(logfile[[2]])<-c(
		"peakpicking","qc","recal","norm","align","profiling","trendblind","pattern",
		"replicates","IS_screen","target_screen","LOD","quantification","blinds","IS_normaliz","IS_subtr","target_subtr",
		"blind_subtr","calibration","isotopologues","adducts","homologues"
	)	
    names(logfile)[2]<-c("Tasks_to_redo"); 
    # summary project status ###################################################
    tasks<-c(
		"Data available?",
		"peakpicking",
		"qc",
		"pattern",
		"recal",
		"align",
		"norm",
		"profiling",
		"IS_normaliz",
		"trendblind",
		"replicates",
		"IS_screen",
		"target_screen",
		"LOD",
		"quantification",
		"blinds",
		"IS_subtr",
		"target_subtr",
		"blind_subtr",
		"calibration",
		"isotopologues",
		"adducts",
		"homologues"
	 )
    doneit<-rep(FALSE,length(tasks))
    summar<-data.frame(tasks,doneit)
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
		logfile$parameters[[79]]<-"30";						names(logfile$parameters)[79]<-"recal_maxdmz"
		# replicate intersection ################################################
		logfile$parameters[[15]]<-"3";						names(logfile$parameters)[15]<-"replicate_dmz"
		logfile$parameters[[16]]<-"TRUE";					names(logfile$parameters)[16]<-"replicate_ppm"		
		logfile$parameters[[17]]<-"FALSE";					names(logfile$parameters)[17]<-"replicate_recalib"		
		logfile$parameters[[18]]<-"30";						names(logfile$parameters)[18]<-"replicate_delRT"	
		logfile$parameters[[19]]<-"5";						names(logfile$parameters)[19]<-"replicate_IS_dInt"
		# trend detection ######################################################
		logfile$parameters[[29]]<-"TRUE";		names(logfile$parameters)[29]<-"notrend"
		logfile$parameters[[34]]<-"4,7,14"; 	names(logfile$parameters)[34]<-"trend_lags" 
		logfile$parameters[[35]]<-"3";			names(logfile$parameters)[35]<-"trend_threshold"
		# blind subtraction ####################################################	
		logfile$parameters[[36]]<-"yes";		names(logfile$parameters)[36]<-"trend_blind"			
		logfile$parameters[[37]]<-"100";		names(logfile$parameters)[37]<-"blind_threshold"	
		logfile$parameters[[82]]<-"3";			names(logfile$parameters)[82]<-"blind_dmz"
		logfile$parameters[[83]]<-"TRUE";		names(logfile$parameters)[83]<-"blind_ppm"		
		logfile$parameters[[84]]<-"30";			names(logfile$parameters)[84]<-"blind_drt"		
		logfile$parameters[[85]]<-"FALSE";		names(logfile$parameters)[85]<-"subtract_pos_bydate"
		logfile$parameters[[86]]<-"FALSE";		names(logfile$parameters)[86]<-"subtract_pos_byfile"
		logfile$parameters[[87]]<-"FALSE";		names(logfile$parameters)[87]<-"subtract_neg_bydate"
		logfile$parameters[[88]]<-"FALSE";		names(logfile$parameters)[88]<-"subtract_neg_byfile"	
		logfile$parameters[[89]]<-"no";			names(logfile$parameters)[89]<-"blind_omit"	
		# profiling ############################################################
		logfile$parameters[[38]]<-"100";		names(logfile$parameters)[38]<-"prof_maxfiles"	
		logfile$parameters[[80]]<-"FALSE";		names(logfile$parameters)[80]<-"upto_file"
		logfile$parameters[[39]]<-"3";			names(logfile$parameters)[39]<-"prof_dmz"
		logfile$parameters[[40]]<-"TRUE";		names(logfile$parameters)[40]<-"prof_ppm"
		logfile$parameters[[41]]<-"60";			names(logfile$parameters)[41]<-"prof_drt"
		logfile$parameters[[90]]<-"FALSE";		names(logfile$parameters)[90]<-"prof_select"
		logfile$parameters[[91]]<-"yes";		names(logfile$parameters)[91]<-"replicates_prof"
		# IS screening #########################################################
		logfile$parameters[[42]]<-"30"; 		names(logfile$parameters)[42]<-"IS_drt1"	# RT tolerance of peaks in sample relative to their expected RT [s]
		logfile$parameters[[43]]<-"10"; 		names(logfile$parameters)[43]<-"IS_drt2"	# RT tolerance of peaks within an isotope pattern [s]
		logfile$parameters[[45]]<-"3";			names(logfile$parameters)[45]<-"IS_dmz"# m/z tolerance ...
		logfile$parameters[[46]]<-"TRUE";		names(logfile$parameters)[46]<-"IS_ppm"# ... given in pppm?
		logfile$parameters[[47]]<-"30";			names(logfile$parameters)[47]<-"IS_inttol"# Intensity tolerance %
		logfile$parameters[[48]]<-"5E4";		names(logfile$parameters)[48]<-"IS_intcut"	# Lower intensity threhold
		logfile$parameters[[49]]<-"0.8";    	names(logfile$parameters)[49]<-"IS_w1" # Matching score
		logfile$parameters[[50]]<-"FALSE";    	names(logfile$parameters)[50]<-"screen_IS_cutit" # Cut off match combiantions below matching score?	
		logfile$parameters[[51]]<-"FALSE";    	names(logfile$parameters)[51]<-"screen_IS_maxonly" # Screen only most intense isotopologue peak?		
		# target screening #####################################################
		logfile$parameters[[55]]<-"30"; 		names(logfile$parameters)[55]<-"tar_drt1"	# RT tolerance of peaks in sample relative to their expected RT [s]
		logfile$parameters[[56]]<-"10"; 		names(logfile$parameters)[56]<-"tar_drt2"	# RT tolerance of peaks within an isotope pattern [s]
		logfile$parameters[[58]]<-"3";			names(logfile$parameters)[58]<-"tar_dmz"# m/z tolerance ...
		logfile$parameters[[59]]<-"TRUE";		names(logfile$parameters)[59]<-"tar_ppm"# ... given in pppm?
		logfile$parameters[[60]]<-"30";			names(logfile$parameters)[60]<-"tar_inttol"# Intensity tolerance %
		logfile$parameters[[61]]<-"5E4";		names(logfile$parameters)[61]<-"tar_intcut"	# Lower intensity threhold
		logfile$parameters[[62]]<-"0.8";    	names(logfile$parameters)[62]<-"tar_w1" # Matching score	
		logfile$parameters[[63]]<-"FALSE";    	names(logfile$parameters)[63]<-"screen_target_cutit" # Cut off match combiantions below matching score?		
		logfile$parameters[[64]]<-"FALSE";    	names(logfile$parameters)[64]<-"screen_target_maxonly" # Screen only most intense isotopologue peak?			
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
	logfile$workflow[5]<-"yes"; 	names(logfile$workflow)[5]<-"pattern" 
	logfile$workflow[6]<-"yes"; 	names(logfile$workflow)[6]<-"LOD" 
	logfile$workflow[7]<-"yes"; 	names(logfile$workflow)[7]<-"peakpicking" 
	logfile$workflow[8]<-"yes"; 	names(logfile$workflow)[8]<-"quantification" 
	logfile$workflow[9]<-"yes"; 	names(logfile$workflow)[9]<-"profiling" 
	logfile$workflow[10]<-"yes"; 	names(logfile$workflow)[10]<-"trendblind"     
	logfile$workflow[11]<-"yes"; 	names(logfile$workflow)[11]<-"IS_screen" 
	logfile$workflow[12]<-"yes"; 	names(logfile$workflow)[12]<-"target_screen" 
	logfile$workflow[13]<-"yes"; 	names(logfile$workflow)[13]<-"replicates" 
	logfile$workflow[14]<-"yes"; 	names(logfile$workflow)[14]<-"blinds" 		
	logfile$workflow[15]<-"yes"; 	names(logfile$workflow)[15]<-"IS_normaliz" 
	logfile$workflow[16]<-"yes"; 	names(logfile$workflow)[16]<-"IS_subtr" 
	logfile$workflow[17]<-"yes"; 	names(logfile$workflow)[17]<-"target_subtr" 
	logfile$workflow[18]<-"yes"; 	names(logfile$workflow)[18]<-"blind_subtr" 
	logfile$workflow[19]<-"yes"; 	names(logfile$workflow)[19]<-"calibration" 	
	logfile$workflow[20]<-"yes"; 	names(logfile$workflow)[20]<-"isotopologues"	
	logfile$workflow[21]<-"yes"; 	names(logfile$workflow)[21]<-"adducts"	
	logfile$workflow[22]<-"yes"; 	names(logfile$workflow)[22]<-"homologues"	
	################################################################################################
	# define matrix of downstream workflow dependencies and ########################################
	# recalculations of previous steps if their results are overwritten, e.g. IS_subtr or ##########
	# target_subtr #################################################################################
	# requires only a definition of direct ones - inderect ones will be in workflow_set.r ##########
	# below specified in a row-wise fashion (but stored columnwise): ###############################
	# define workflow order of logfile$Tasks_to_redo by server.calculation.r #######################
	# dependencies must simply go after their parent node ########################################## 
	# order here actually irrelevant, because calculation order set in server_calculation  #########	
	work_names<-names(logfile$Tasks_to_redo)[1:22]
	depend<-matrix(ncol=length(work_names),nrow=length(work_names),0)
	colnames(depend)<-work_names
	rownames(depend)<-work_names					# peakpicking	qc	recal	norm	align	profiling	trendblind	pattern		replicates	IS_screen	target_screeen	LOD		calibration	quantification	blinds IS_normaliz	IS_subtr	target_subtr	blind_subtr	isotopologues	adducts homologues
	depend[,colnames(depend)=="peakpicking"]<-		c(0,			1,	1,		1,		1,		1,			1,			0,			1,			1,			1,				1,		1,			1,				1,		1,			1,			1,				1,			0,				0,		0)
	depend[,colnames(depend)=="qc"]<-				c(0,			0,	1,		1,		1,		1,			1,			0,			1,			1,			1,				1,		1,			1,				1,		1,			1,			1,				1,			0,				0,		0)
	depend[,colnames(depend)=="pattern"]<-			c(0,			0,	1,		0,		0,		0,			0,			0,			0,			1,			1,				0,		1,			1,				0,		1,			1,			1,				0,			0,				0,		0)
	depend[,colnames(depend)=="recal"]<-			c(0,			0,	0,		0,		0,		1,			1,			0,			1,			1,			1,				0,		1,			1,				0,		1,			1,			1,				1,			0,				0,		0)
	depend[,colnames(depend)=="align"]<-			c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		1,			1,				0,		0,			0,			0,				1,			0,				0,		0)
	depend[,colnames(depend)=="norm"]<-				c(0,			0,	0,		0,		0,		1,			1,			0,			0,			1,			1,				1,		1,			1,				1,		1,			0,			0,				1,			0,				0,		0)
	depend[,colnames(depend)=="blinds"]<-			c(0,			0,	0,		0,		0,		1,			1,			0,			0,			1,			1,				0,		1,			1,				0,		1,			1,			1,				1,			0,				0,		0)
	depend[,colnames(depend)=="replicates"]<-		c(0,			0,	0,		0,		0,		1,			1,			0,			0,			1,			1,				1,		1,			1,				0,		1,			1,			1,				1,			0,				0,		0)
	depend[,colnames(depend)=="profiling"]<-		c(0,			0,	0,		0,		0,		0,			1,			0,			0,			1,			1,				0,		0,			1,				0,		1,			1,			1,				1,			0,				0,		0)
	depend[,colnames(depend)=="IS_screen"]<-		c(0,			0,	0,		0,		0,		1,			0,			0,			0,			0,			0,				0,		0,			1,				0,		1,			1,			0,				0,			0,				0,		0)
	depend[,colnames(depend)=="target_screen"]<-	c(0,			0,	0,		0,		0,		1,			0,			0,			0,			0,			0,				0,		0,			1,				0,		1,			0,			1,				0,			0,				0,		0)
	depend[,colnames(depend)=="IS_normaliz"]<-		c(0,			0,	0,		0,		0,		0,			1,			0,			0,			0,			0,				0,		0,			0,				0,		1,			0,			0,				0,			0,				0,		0)
	depend[,colnames(depend)=="trendblind"]<-		c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	depend[,colnames(depend)=="LOD"]<-				c(0,			0,	0,		0,		0,		0,			0,			0,			0,			1,			1,				0,		1,			0,				0,		1,			1,			1,				0,			0,				0,		0)
	depend[,colnames(depend)=="calibration"]<-		c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	depend[,colnames(depend)=="quantification"]<-	c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	depend[,colnames(depend)=="IS_subtr"]<-			c(0,			0,	0,		0,		0,		1,			1,			0,			0,			0,			0,				0,		0,			0,				0,		0,			0,			0,				1,			0,				0,		0)
	depend[,colnames(depend)=="target_subtr"]<-		c(0,			0,	0,		0,		0,		1,			1,			0,			0,			0,			0,				0,		0,			0,				0,		0,			0,			0,				1,			0,				0,		0)
	depend[,colnames(depend)=="blind_subtr"]<-		c(0,			0,	0,		0,		0,		1,			1,			0,			0,			0,			0,				0,		0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	depend[,colnames(depend)=="isotopologues"]<-	c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	depend[,colnames(depend)=="adducts"]<-			c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	depend[,colnames(depend)=="homologues"]<-		c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	logfile[[11]]<-depend
	names(logfile)[11]<-"workflow_depend"
	################################################################################################
	# define upstream workflow "musts", i.e., upstream nodes on which`s execution a node ###########
	# depends. 0 = not dependent. 1 = dependent. -1 = MUST NOT be executed ######################### 
	must<-matrix(ncol=length(work_names),nrow=length(work_names),0)
	colnames(must)<-work_names
	rownames(must)<-work_names					# peakpicking	qc	recal	norm	align	profiling	trendblind	pattern		replicates	IS_screen	target_screeen	LOD		calibration		quantification	blinds IS_normaliz	IS_subtr	target_subtr	blind_subtr	isotopologues	adducts homologues
	must[,colnames(must)=="peakpicking"]<-		c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,				0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="qc"]<-				c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,				0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="pattern"]<-			c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,				0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="recal"]<-			c(1,			0,	0,		0,		0,		0,			0,			1,			0,			0,			0,				0,		0,				0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="align"]<-			c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,				0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="norm"]<-				c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,				0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="blinds"]<-			c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,				0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="replicates"]<-		c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,				0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="profiling"]<-		c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,				0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="IS_screen"]<-		c(1,			0,	0,		0,		0,		0,			0,			1,			0,			0,			0,				0,		0,				0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="target_screen"]<-	c(1,			0,	0,		0,		0,		0,			0,			1,			0,			0,			0,				0,		0,				0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="IS_normaliz"]<-		c(1,			0,	0,		0,		0,		0,			0,			1,			0,			0,			0,				0,		0,				0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="trendblind"]<-		c(1,			0,	0,		0,		0,		1,			0,			0,			0,			0,			0,				0,		0,				0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="LOD"]<-				c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,				0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="calibration"]<-		c(1,			0,	0,		0,		0,		0,			0,			1,			0,			0,			0,				0,		0,				0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="quantification"]<-	c(1,			0,	0,		0,		0,		0,			0,			1,			0,			1,			1,				0,		0,				0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="IS_subtr"]<-			c(1,			0,	0,		0,		0,		1,			0,			1,			0,			1,			1,				0,		0,				0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="target_subtr"]<-		c(1,			0,	0,		0,		0,		1,			0,			1,			0,			1,			1,				0,		0,				0,				0,		0,			0,			0,				0,			0,				0,		0)	
	must[,colnames(must)=="blind_subtr"]<-		c(1,			0,	0,		0,		0,		1,			0,			0,			0,			0,			0,				0,		0,				0,				1,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="isotopologues"]<-	c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,				0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="adducts"]<-			c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,				0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="homologues"]<-		c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,				0,				0,		0,			0,			0,				0,			0,				0,		0)
	logfile[[12]]<-must
	names(logfile)[12]<-"workflow_must"	
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
    logfile[[9]]<-0   
    names(logfile)[9]<-c("isotopes")
    logfile[[9]]<-"";
	# enviMass version number ##################################################
    logfile[[10]]<-3.101 
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