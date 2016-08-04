# include updates - if older projects are reopened
if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_updates.r!")}
#stop("\n\nMaintenance work; enviMass will be back in a couple of hours! Please update again later.")
#if(logfile[[10]]<3.100){
#if(logfile[[10]]<3.102){
if(TRUE){
	
	cat("\n Updating to version 3.100 ...")
	################################################################################################
	# create missing folder
	if(!file.exists(file.path(logfile$project_folder,"results","screening"))
	){
		dir.create(file.path(logfile$project_folder,"results","screening"),recursive=TRUE)    	# subfolder  
	}
	if(!file.exists(file.path(logfile$project_folder,"quantification"))
	){	
		dir.create(file.path(logfile$project_folder,"quantification"),recursive=TRUE)   # subfolder 
	}
	if(!file.exists(file.path(logfile$project_folder,"results","LOD"))
	){
		dir.create(file.path(logfile$project_folder,"results","LOD"),recursive=TRUE)    	# subfolder  
	}
	if(!file.exists(file.path(logfile$project_folder,"results","recalibration"))
	){
		dir.create(file.path(logfile$project_folder,"results","recalibration"),recursive=TRUE)    	# subfolder  
	}

	# another column in peaklists for the replicates!
	IDs<-list.files(file.path(logfile[[1]],"peaklist"))
	if(length(IDs)>0){
		for(i in 1:length(IDs)){
			load(file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])),envir=as.environment(".GlobalEnv"),verbose=FALSE);
			if(any(colnames(peaklist)=="keep")){break} # ok, has been done before
			keep<-rep(1,length(peaklist[,1])) # 1 == TRUE
			peaklist<-cbind(peaklist,keep)
			colnames(peaklist)[15]<-"keep_2";
			save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])))
			rm(peaklist)
		}
	}
	# extend logfile$summary
	first_ext<-FALSE;
	if(!any(logfile$summary[,1]=="replicates")){
		logfile$summary[,1]<<-as.character(logfile$summary[,1])
		logfile$summary[,2]<<-as.character(logfile$summary[,2])
		logfile$summary[11,1]<<-"replicates"
		logfile$summary[11,2]<<-"FALSE"	
		first_ext<-TRUE;
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
	}
	if(!any(logfile$summary[,1]=="IS_screen")){
		logfile$summary[,1]<<-as.character(logfile$summary[,1])
		logfile$summary[,2]<<-as.character(logfile$summary[,2])
		logfile$summary[12,1]<<-"IS_screen"
		logfile$summary[12,2]<<-"FALSE"	
		first_ext<-TRUE;
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
	}
	if(!any(logfile$summary[,1]=="target_screen")){
		logfile$summary[,1]<<-as.character(logfile$summary[,1])
		logfile$summary[,2]<<-as.character(logfile$summary[,2])
		logfile$summary[13,1]<<-"target_screen"
		logfile$summary[13,2]<<-"FALSE"	
		first_ext<-TRUE;
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
	}
	if(!any(logfile$summary[,1]=="LOD")){
		logfile$summary[,1]<<-as.character(logfile$summary[,1])
		logfile$summary[,2]<<-as.character(logfile$summary[,2])
		logfile$summary[14,1]<<-"LOD"
		logfile$summary[14,2]<<-"FALSE"	
		first_ext<-TRUE;
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
	}
	if(!any(logfile$summary[,1]=="quantification")){
		logfile$summary[,1]<<-as.character(logfile$summary[,1])
		logfile$summary[,2]<<-as.character(logfile$summary[,2])
		logfile$summary[15,1]<<-"quantification"
		logfile$summary[15,2]<<-"FALSE"	
		first_ext<-TRUE;
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
	}
	if(!any(logfile$summary[,1]=="blinds")){
		logfile$summary[,1]<<-as.character(logfile$summary[,1])
		logfile$summary[,2]<<-as.character(logfile$summary[,2])
		logfile$summary[16,1]<<-"blinds"
		logfile$summary[16,2]<<-"FALSE"	
		logfile$parameters[[82]]<<-"3";			names(logfile$parameters)[82]<<-"blind_dmz"
		logfile$parameters[[83]]<<-"TRUE";		names(logfile$parameters)[83]<<-"blind_ppm"		
		logfile$parameters[[84]]<<-"30";			names(logfile$parameters)[84]<<-"blind_drt"	
		first_ext<-TRUE;
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
	}
	if(!any(logfile$summary[,1]=="IS_subtr")){	
		logfile$summary[17,1]<<-"IS_subtr"
		logfile$summary[17,2]<<-"FALSE"	
		first_ext<-TRUE;
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))		
	}
	if(!any(logfile$summary[,1]=="target_subtr")){	
		logfile$summary<<-rbind(logfile$summary,c("",""))
		logfile$summary[18,1]<<-"target_subtr"
		logfile$summary[18,2]<<-"FALSE"	
		first_ext<-TRUE;
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
	}
	if(!any(logfile$summary[,1]=="blind_subtr")){	
		logfile$summary<<-rbind(logfile$summary,c("",""))
		logfile$summary[19,1]<<-"blind_subtr"
		logfile$summary[19,2]<<-"FALSE"	
		first_ext<-TRUE;
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
	}	
	if(first_ext){
		logfile$summary[,1]<<-c(
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
			"blind_subtr"
		)
	}
	
	# insert missing parameters
	if(!any(names(logfile$parameters)=="replicate_dmz")){
		logfile$parameters[[15]]<<-"3";names(logfile$parameters)[15]<<-"replicate_dmz"
		logfile$parameters[[16]]<<-"TRUE";names(logfile$parameters)[16]<<-"replicate_ppm"		
		logfile$parameters[[17]]<<-"FALSE";names(logfile$parameters)[17]<<-"replicate_recalib"		
		logfile$parameters[[18]]<<-"30";names(logfile$parameters)[18]<<-"replicate_delRT"	
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
	}
	if(!any(names(logfile$parameters)=="screen_IS_cutit")){
		logfile$parameters[[50]]<<-"FALSE";    	names(logfile$parameters)[50]<<-"screen_IS_cutit" # Cut off match combiantions below matching score?		
		logfile$parameters[[63]]<<-"FALSE";    	names(logfile$parameters)[63]<<-"screen_target_cutit" # Cut off match combiantions below matching score?			
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
	}
	if(!any(names(logfile$parameters)=="peak_perc_cut")){	
		logfile$parameters[[92]]<<-"0"; 	names(logfile$parameters)[92]<<-"peak_perc_cut"  
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
	}
	# logfile$Tasks_to_redo ##################################################################
	if(length(logfile[[2]])<23){
		logfile[[2]]<<-rep("TRUE",23)
	}
	names(logfile[[2]])<<-c(
		"peakpicking","qc","recal","norm","align","profiling","trendblind","pattern",
		"replicates","IS_screen","target_screen","LOD","calibration","recovery","quantification","blinds","IS_normaliz","IS_subtr","target_subtr",
		"blind_subtr","isotopologues","adducts","homologues"
	)	
	names(logfile)[2]<<-c("Tasks_to_redo"); 
	save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
	load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
	# logfile$workflow ######################################################################
	if(!any(names(logfile$workflow)=="IS_screen")){
		logfile$workflow[11]<<-"yes"; 	names(logfile$workflow)[11]<<-"IS_screen" 
	}
	if(!any(names(logfile$workflow)=="target_screen")){
		logfile$workflow[12]<<-"yes"; 	names(logfile$workflow)[12]<<-"target_screen" 
	}
	if(!any(names(logfile$workflow)=="replicates")){
		logfile$workflow[13]<<-"yes"; 	names(logfile$workflow)[13]<<-"replicates" 
	}
	if(!any(names(logfile$workflow)=="LOD")){
		logfile$workflow[6]<<-"yes"; 	names(logfile$workflow)[6]<<-"LOD" 
	}
	if(!any(names(logfile$workflow)=="quantification")){
		logfile$workflow[8]<<-"yes"; 	names(logfile$workflow)[8]<<-"quantification" 
	}
	if(!any(names(logfile$workflow)=="blinds")){
		logfile$workflow[14]<<-"yes"; 	names(logfile$workflow)[14]<<-"blinds" 
	}
	if(any(names(logfile$workflow)=="profnorm")){
		logfile$workflow[15]<<-"yes"; 	names(logfile$workflow)[15]<<-"IS_normaliz" 
	}	
	if(any(names(logfile$workflow)=="profiled")){
		logfile$workflow[15]<<-"yes"; 	names(logfile$workflow)[9]<<-"profiling" 
	}	
	logfile$workflow[5]<<-"yes"; 	names(logfile$workflow)[5]<<-"pattern" 
	logfile$workflow[7]<<-"yes"; 	names(logfile$workflow)[7]<<-"peakpicking" 		
									names(logfile$workflow)[10]<<-"trendblind" 		
	
	if(!any(names(logfile$workflow)=="IS_subtr")){
		logfile$workflow[16]<<-"yes"; 	names(logfile$workflow)[16]<<-"IS_subtr" 
	}
	if(!any(names(logfile$workflow)=="target_subtr")){
		logfile$workflow[17]<<-"yes"; 	names(logfile$workflow)[17]<<-"target_subtr" 
	}	
	if(!any(names(logfile$workflow)=="blind_subtr")){
		logfile$workflow[18]<<-"yes"; 	names(logfile$workflow)[18]<<-"blind_subtr" 
	}	
	# enforce a pattern recalculation - zero abundance debug
	if(TRUE){
		redo_pattern<-FALSE
		if( file.exists(file.path(logfile[[1]],"results","pattern_pos_IS")) ){	
			load(file=file.path(logfile[[1]],"results","pattern_pos_IS"),envir=as.environment(".GlobalEnv"));
			for(i in 1:length(pattern_pos_IS)){
				if( any(pattern_pos_IS[[i]][,2]==0) ){
					redo_pattern<-TRUE
				}
			}
			rm(pattern_pos_IS,envir=as.environment(".GlobalEnv"))
		}
		if( file.exists(file.path(logfile[[1]],"results","pattern_neg_IS")) ){	
			load(file=file.path(logfile[[1]],"results","pattern_neg_IS"),envir=as.environment(".GlobalEnv"));
			for(i in 1:length(pattern_neg_IS)){
				if( any(pattern_neg_IS[[i]][,2]==0) ){
					redo_pattern<-TRUE
				}
			}
			rm(pattern_neg_IS,envir=as.environment(".GlobalEnv"))
		}
		if( file.exists(file.path(logfile[[1]],"results","pattern_pos_target")) ){	
			load(file=file.path(logfile[[1]],"results","pattern_pos_target"),envir=as.environment(".GlobalEnv"));
			for(i in 1:length(pattern_pos_target)){
				if( any(pattern_pos_target[[i]][,2]==0) ){
					redo_pattern<-TRUE
				}
			}
			rm(pattern_pos_target,envir=as.environment(".GlobalEnv"))
		}
		if( file.exists(file.path(logfile[[1]],"results","pattern_neg_target")) ){	
			load(file=file.path(logfile[[1]],"results","pattern_neg_target"),envir=as.environment(".GlobalEnv"));
			for(i in 1:length(pattern_neg_target)){
				if( any(pattern_neg_target[[i]][,2]==0) ){
					redo_pattern<-TRUE
				}
			}
			rm(pattern_neg_target,envir=as.environment(".GlobalEnv"))
		}
		if(redo_pattern){
			logfile$workflow[names(logfile$workflow)=="pattern"]<<-"yes"; 	
			logfile$summary[logfile$summary[,1]=="pattern",2]<<-"TRUE";
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="pattern"]<<-TRUE
		}
	}
	################################################################################################
	# define matrix of downstream workflow dependencies ############################################
	# (or with recalculations of previous steps if previous results need to be again written to, ###
	# e.g. IS_subtr or target_subtr) ###############################################################
	# requires only a definition of direct ones - indirect ones will be in workflow_set.r ##########
	# below specified in a row-wise fashion (but stored and retrieved columnwise): #################
	# define workflow order of logfile$Tasks_to_redo by server.calculation.r #######################
	# dependencies must simply go after their parent node ########################################## 
	# order here actually irrelevant, because calculation order set in server_calculation  #########	
	work_names<-names(logfile$Tasks_to_redo)[1:23]
	depend<-matrix(ncol=length(work_names),nrow=length(work_names),0)
	colnames(depend)<-work_names
	rownames(depend)<-work_names					# peakpicking	qc	recal	norm	align	profiling	trendblind	pattern		replicates	IS_screen	target_screeen	LOD		calibration	recovery	quantification	blinds IS_normaliz	IS_subtr	target_subtr	blind_subtr	isotopologues	adducts homologues
	depend[,colnames(depend)=="peakpicking"]<-		c(0,			1,	1,		1,		1,		1,			1,			0,			1,			1,			1,				1,		1,			0,			1,				1,		1,			1,			1,				1,			0,				0,		0)
	depend[,colnames(depend)=="qc"]<-				c(0,			0,	1,		1,		1,		1,			1,			0,			1,			1,			1,				1,		1,			0,			1,				1,		1,			1,			1,				1,			0,				0,		0)
	depend[,colnames(depend)=="pattern"]<-			c(0,			0,	1,		0,		0,		0,			0,			0,			0,			1,			1,				0,		1,			0,			1,				0,		1,			1,			1,				0,			0,				0,		0)
	depend[,colnames(depend)=="recal"]<-			c(0,			0,	0,		0,		0,		1,			1,			0,			1,			1,			1,				0,		1,			0,			1,				0,		1,			1,			1,				1,			0,				0,		0)
	depend[,colnames(depend)=="align"]<-			c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		1,			0,			1,				0,		0,			0,			0,				1,			0,				0,		0)
	depend[,colnames(depend)=="norm"]<-				c(0,			0,	0,		0,		0,		1,			1,			0,			0,			1,			1,				1,		1,			0,			1,				1,		1,			0,			0,				1,			0,				0,		0)
	depend[,colnames(depend)=="blinds"]<-			c(0,			0,	0,		0,		0,		1,			1,			0,			0,			1,			1,				0,		1,			0,			1,				0,		1,			1,			1,				1,			0,				0,		0)
	depend[,colnames(depend)=="replicates"]<-		c(0,			0,	0,		0,		0,		1,			1,			0,			0,			1,			1,				1,		1,			0,			1,				0,		1,			1,			1,				1,			0,				0,		0)
	depend[,colnames(depend)=="profiling"]<-		c(0,			0,	0,		0,		0,		0,			1,			0,			0,			1,			1,				0,		0,			0,			1,				0,		1,			1,			1,				1,			0,				0,		0)
	depend[,colnames(depend)=="IS_screen"]<-		c(0,			0,	0,		0,		0,		1,			0,			0,			0,			0,			0,				0,		0,			0,			1,				0,		1,			1,			0,				0,			0,				0,		0)
	depend[,colnames(depend)=="target_screen"]<-	c(0,			0,	0,		0,		0,		1,			0,			0,			0,			0,			0,				0,		0,			0,			1,				0,		1,			0,			1,				0,			0,				0,		0)
	depend[,colnames(depend)=="IS_normaliz"]<-		c(0,			0,	0,		0,		0,		0,			1,			0,			0,			0,			0,				0,		0,			0,			0,				0,		1,			0,			0,				0,			0,				0,		0)
	depend[,colnames(depend)=="trendblind"]<-		c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	depend[,colnames(depend)=="LOD"]<-				c(0,			0,	0,		0,		0,		0,			0,			0,			0,			1,			1,				0,		1,			0,			0,				0,		1,			1,			1,				0,			0,				0,		0)
	depend[,colnames(depend)=="calibration"]<-		c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	depend[,colnames(depend)=="recovery"]<-			c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	depend[,colnames(depend)=="quantification"]<-	c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			1,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	depend[,colnames(depend)=="IS_subtr"]<-			c(0,			0,	0,		0,		0,		1,			1,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				1,			0,				0,		0)
	depend[,colnames(depend)=="target_subtr"]<-		c(0,			0,	0,		0,		0,		1,			1,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				1,			0,				0,		0)
	depend[,colnames(depend)=="blind_subtr"]<-		c(0,			0,	0,		0,		0,		1,			1,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	depend[,colnames(depend)=="isotopologues"]<-	c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	depend[,colnames(depend)=="adducts"]<-			c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	depend[,colnames(depend)=="homologues"]<-		c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	logfile[[11]]<<-depend
	names(logfile)[11]<<-"workflow_depend"
	################################################################################################
	# define upstream workflow "musts", i.e., upstream nodes on which`s execution a node ###########
	# depends. 0 = not dependent. 1 = dependent. -1 = MUST NOT be executed ######################### 
	must<-matrix(ncol=length(work_names),nrow=length(work_names),0)
	colnames(must)<-work_names
	rownames(must)<-work_names					# peakpicking	qc	recal	norm	align	profiling	trendblind	pattern		replicates	IS_screen	target_screeen	LOD		calibration	recovery	quantification	blinds IS_normaliz	IS_subtr	target_subtr	blind_subtr	isotopologues	adducts homologues
	must[,colnames(must)=="peakpicking"]<-		c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="qc"]<-				c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="pattern"]<-			c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="recal"]<-			c(1,			0,	0,		0,		0,		0,			0,			1,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="align"]<-			c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="norm"]<-				c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="blinds"]<-			c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="replicates"]<-		c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="profiling"]<-		c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="IS_screen"]<-		c(1,			0,	0,		0,		0,		0,			0,			1,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="target_screen"]<-	c(1,			0,	0,		0,		0,		0,			0,			1,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="IS_normaliz"]<-		c(1,			0,	0,		0,		0,		0,			0,			1,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="trendblind"]<-		c(1,			0,	0,		0,		0,		1,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="LOD"]<-				c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="calibration"]<-		c(1,			0,	0,		0,		0,		0,			0,			1,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="recovery"]<-			c(1,			0,	0,		0,		0,		0,			0,			1,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="quantification"]<-	c(1,			0,	0,		0,		0,		0,			0,			1,			0,			1,			1,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="IS_subtr"]<-			c(1,			0,	0,		0,		0,		1,			0,			1,			0,			1,			1,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="target_subtr"]<-		c(1,			0,	0,		0,		0,		1,			0,			1,			0,			1,			1,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)	
	must[,colnames(must)=="blind_subtr"]<-		c(1,			0,	0,		0,		0,		1,			0,			0,			0,			0,			0,				0,		0,			0,			0,				1,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="isotopologues"]<-	c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="adducts"]<-			c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	must[,colnames(must)=="homologues"]<-		c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
	logfile[[12]]<<-must
	names(logfile)[12]<<-"workflow_must"	
	################################################################################################		
	logfile[[10]]<<-3.100
	names(logfile)[10]<<-"version"
	################################################################################################
	save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
	load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv")) 

}

if(logfile[[10]]<3.101){	

	cat("\n Updating to version 3.101 ...")
	################################################################################################
	# create missing folder
	if(!file.exists(file.path(logfile$project_folder,"quantification"))){
		dir.create(file.path(logfile$project_folder,"quantification"),recursive=TRUE) # subfolder
	}
 	################################################################################################
	# updating columns in IS compound table ########################################################
	intstand<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
	if(any(names(intstand)=="Lower intensity bound")){# remove misnomer
		intstand<-intstand[,!(names(intstand)=="Lower intensity bound")]
		intstand<-intstand[,!(names(intstand)=="Lower intensity bound")]
	}
	if(any(names(intstand)=="Lower.intensity.bound")){# remove misnomer
		intstand<-intstand[,!(names(intstand)=="Lower.intensity.bound")]
		intstand<-intstand[,!(names(intstand)=="Upper.intensity.bound")]
	}
	if(!any(names(intstand)=="Lower_intensity_bound")){
		names_1<-names(intstand)
		intstand<-cbind(intstand,rep(0,length(intstand[,1])),rep(Inf,length(intstand[,1])))
		names(intstand)<-c(names_1,"Lower_intensity_bound","Upper_intensity_bound")
	}	
	if(any(names(intstand)=="Use_for_screening")){
		names(intstand)[names(intstand)=="Use_for_screening"]<-"use_for_screening"
	}	
	if(any(names(intstand)=="Use_for_recalibration")){
		names(intstand)[names(intstand)=="Use_for_recalibration"]<-"use_for_recalibration"
	}		
	write.table(intstand,file=file.path(logfile[[1]],"dataframes","IS.txt"),row.names=FALSE,sep="\t",quote=FALSE)
	rm(intstand)
	################################################################################################
	# insert missing parameters ####################################################################
	if(!any(names(logfile$parameters)=="subtract_pos_bydate")){
		logfile$parameters[[85]]<<-"FALSE";		names(logfile$parameters)[85]<<-"subtract_pos_bydate"
		logfile$parameters[[86]]<<-"FALSE";		names(logfile$parameters)[86]<<-"subtract_pos_byfile"
		logfile$parameters[[87]]<<-"FALSE";		names(logfile$parameters)[87]<<-"subtract_neg_bydate"
		logfile$parameters[[88]]<<-"FALSE";		names(logfile$parameters)[88]<<-"subtract_neg_byfile"	
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
	}
	if(!any(names(logfile$parameters)=="blind_omit")){
		logfile$parameters[[89]]<<-"no"; 		names(logfile$parameters)[89]<<-"blind_omit"	
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
	}		
	if(!any(names(logfile$parameters)=="prof_select")){
		logfile$parameters[[90]]<<-"FALSE"; 	names(logfile$parameters)[90]<<-"prof_select"
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
	}	
	if(!any(names(logfile$parameters)=="trend_blind")){	
		logfile$parameters[[36]]<<-"yes"; 		names(logfile$parameters)[36]<<-"trend_blind"		
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))	
	}
	if(!any(names(logfile$parameters)=="replicate_IS_dInt")){	
		logfile$parameters[[19]]<<-"5";		names(logfile$parameters)[19]<<-"replicate_IS_dInt"	
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))	
	}				
	if(!any(names(logfile$parameters)=="replicates_prof")){	
		logfile$parameters[[91]]<<-"yes";		names(logfile$parameters)[91]<<-"replicates_prof"
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))	
	}
	if(!any(names(logfile$parameters)=="screen_IS_maxonly")){	
		logfile$parameters[[64]]<<-"FALSE";    	names(logfile$parameters)[64]<<-"screen_target_maxonly" # Screen only most intense isotopologue peak?	
		logfile$parameters[[51]]<<-"FALSE";    	names(logfile$parameters)[51]<<-"screen_IS_maxonly" # Screen only most intense isotopologue peak?		
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))	
	}
	################################################################################################
	# updating columns in targets compound table ###################################################
	targets<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");
	if(!any(names(targets)=="warn_1")){
		names_1<-names(targets)
		targets<-cbind(targets,rep("FALSE",length(targets[,1])),rep("FALSE",length(targets[,1])))
		names(targets)<-c(names_1,"warn_1","warn_2")
	}
	if(any(names(targets)=="Use_for_screening")){
		names(targets)[names(targets)=="Use_for_screening"]<-"use_for_screening"
	}	
	if(any(names(targets)=="Use_for_recalibration")){
		names(targets)[names(targets)=="Use_for_recalibration"]<-"use_for_recalibration"
	}		
	if(any(names(targets)=="intensity_warn_1")){	# remove misnomer
		targets<-targets[,!(names(targets)=="intensity_warn_1")]
		targets<-targets[,!(names(targets)=="intensity_warn_2")]		
	}
	if(any(names(targets)=="RT.tolerance")){	# remove typo
		names(targets)[names(targets)=="RT.tolerance"]<-"RT_tolerance"
	}	
	if(any(names(targets)=="intercept")){	# remove obsolete column
		targets<-targets[,names(targets)!="intercept"]
	}		
	if(any(names(targets)=="slope")){	# remove obsolete column
		targets<-targets[,names(targets)!="slope"]
	}		
	write.table(targets,file=file.path(logfile[[1]],"dataframes","targets.txt"),row.names=FALSE,sep="\t",quote=FALSE)
	rm(targets)
	################################################################################################
	# another column in peaklists for the blind subtraction ########################################
	IDs<-list.files(file.path(logfile[[1]],"peaklist"))
	if(length(IDs)>0){
		for(i in 1:length(IDs)){
			load(file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])),envir=as.environment(".GlobalEnv"),verbose=FALSE);
			if(any(colnames(peaklist)=="keep_2")){break} # ok, has been done before
			keep_2<-rep(1,length(peaklist[,1])) # 1 == TRUE
			peaklist<-cbind(peaklist,keep_2)
			colnames(peaklist)[16]<-"keep_2";
			save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])))
			rm(peaklist)
		}
	}
	################################################################################################	
	# store subtraction files ######################################################################
	if(!any(names(logfile)=="Positive_subtraction_files")){
		logfile[[13]]<<-"FALSE"
		names(logfile)[13]<<-"Positive_subtraction_files"
	}
	if(!any(names(logfile)=="Negative_subtraction_files")){
		logfile[[14]]<<-"FALSE"
		names(logfile)[14]<<-"Negative_subtraction_files"
	}
	################################################################################################	
	logfile[[10]]<<-3.101
	names(logfile)[10]<<-"version"
	################################################################################################		
	save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
	load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv")) 
}

if(logfile[[10]]<3.102){

	cat("\n Updating to version 3.102 ...")
	################################################################################################
	# updating columns in IS compound table ########################################################
	intstand<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
	if(!any(names(intstand)=="Quant_adduct")){
		names_1<-names(intstand)
		quant_add<-rep("",length(intstand[,1]))
		quant_add[intstand[,names(intstand)=="ion_mode"]=="positive"]<-"M+H"
		quant_add[intstand[,names(intstand)=="ion_mode"]=="negative"]<-"M-H"		
		intstand<-cbind(intstand,quant_add,rep("1",length(intstand[,1])))
		names(intstand)<-c(names_1,"Quant_adduct","Quant_peak")
	}		
	write.table(intstand,file=file.path(logfile[[1]],"dataframes","IS.txt"),row.names=FALSE,sep="\t",quote=FALSE)
	rm(intstand)
	################################################################################################
	# updating columns in target compound table ####################################################
	targets<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");
	if(!any(names(targets)=="Quant_adduct")){
		names_1<-names(targets)
		quant_add<-rep("",length(targets[,1]))
		quant_add[targets[,names(targets)=="ion_mode"]=="positive"]<-"M+H"
		quant_add[targets[,names(targets)=="ion_mode"]=="negative"]<-"M-H"		
		targets<-cbind(targets,quant_add,rep("1",length(targets[,1])))
		names(targets)<-c(names_1,"Quant_adduct","Quant_peak")
	}		
	write.table(targets,file=file.path(logfile[[1]],"dataframes","targets.txt"),row.names=FALSE,sep="\t",quote=FALSE)
	rm(targets)
	################################################################################################
	if(!any(names(logfile$parameters)=="recal_maxdmz")){	
		logfile$parameters[[79]]<<-"30";		names(logfile$parameters)[79]<<-"recal_maxdmz"	
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))	
	}				
	################################################################################################	
	# extend logfile$summary
	if(!any(logfile$summary[,1]=="calibration")){
		logfile$summary[,1]<<-as.character(logfile$summary[,1])
		logfile$summary[,2]<<-as.character(logfile$summary[,2])
		logfile$summary[20,1]<<-"calibration"
		logfile$summary[20,2]<<-"FALSE"	
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
	}
	if(!any(logfile$summary[,1]=="isotopologues")){
		logfile$summary[,1]<<-as.character(logfile$summary[,1])
		logfile$summary[,2]<<-as.character(logfile$summary[,2])
		logfile$summary[21,1]<<-"isotopologues"
		logfile$summary[21,2]<<-"FALSE"	
		logfile$summary[22,1]<<-"adducts"
		logfile$summary[22,2]<<-"FALSE"			
		logfile$summary[23,1]<<-"homologues"
		logfile$summary[23,2]<<-"FALSE"			
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
	}	
	if(!any(logfile$summary[,1]=="recovery")){	
		logfile$summary[24,1]<<-"recovery"
		logfile$summary[24,2]<<-"FALSE"	
		first_ext<-TRUE;
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))		
	}
	################################################################################################
	# extend logfile$workflow
	if(!any(names(logfile$workflow)=="calibration")){	
		logfile$workflow[19]<<-"yes"; 	names(logfile$workflow)[19]<<-"calibration" 		
	}
	if(!any(names(logfile$workflow)=="isotopologues")){		
		logfile$workflow[20]<<-"yes"; 	names(logfile$workflow)[20]<<-"isotopologues"	
	}	
	if(!any(names(logfile$workflow)=="adducts")){		
		logfile$workflow[21]<<-"yes"; 	names(logfile$workflow)[21]<<-"adducts"	
	}	
	if(!any(names(logfile$workflow)=="homologues")){
		logfile$workflow[22]<<-"yes"; 	names(logfile$workflow)[22]<<-"homologues"
	}
	if(!any(names(logfile$workflow)=="recovery")){
		logfile$workflow[23]<<-"yes"; 	names(logfile$workflow)[23]<<-"recovery" 
	}		
	################################################################################################	
	# modify measurements table ####################################################################
	measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	if( !any(names(measurements)=="date_end") ){
		measurements<-cbind(
			measurements,
			rep("2019-06-08",length(measurements[,1])),
			rep("12:00:00",length(measurements[,1]))
		)
	}
	if( any(names(measurements)=="feat.") ){
		names(measurements)[names(measurements)=="feat."]<-"profiled";
		measurements[,names(measurements)=="profiled"]<-"TRUE";
	}
	names(measurements)<-c("ID","Name","Type","Mode","Place","Date","Time","include","copied","picked",
	"checked","recal","align","norm","profiled","LOD","IS_screen","tar_screen","tag1","tag2","tag3","date_end","time_end")
	write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
	rm(measurements)
	################################################################################################	
	logfile[[10]]<<-3.102
	names(logfile)[10]<<-"version"
	################################################################################################		
	save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
	load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv")) 

}

if(any(ls()=="logfile")){stop("\n illegal logfile detected #2 in server_updates.r!")}




















