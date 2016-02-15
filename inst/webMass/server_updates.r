# include updates - if older projects are reopened
if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_updates.r!")}




if(logfile[[10]]<3.100){	
	
	# create missing folder
	if(!file.exists(file.path(logfile$project_folder,"results","screening"))
	){
		dir.create(file.path(logfile$project_folder,"results","screening"),recursive=TRUE)    	# subfolder  
	}
	if(!file.exists(file.path(logfile$project_folder,"results","quantification"))
	){	
		dir.create(file.path(logfile$project_folder,"results","quantification"),recursive=TRUE)   # subfolder 
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
	for(i in 1:length(IDs)){
		load(file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])),envir=as.environment(".GlobalEnv"),verbose=FALSE);
		if(any(colnames(peaklist)=="keep")){break} # ok, has been done before
		keep<-rep(1,length(peaklist[,1])) # 1 == TRUE
		peaklist<-cbind(peaklist,keep)
		save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])))
		rm(peaklist)
	}

	# extend logfile$summary
	if(!any(logfile$summary[,1]=="Replicate filter")){
		logfile$summary[,1]<<-as.character(logfile$summary[,1])
		logfile$summary[,2]<<-as.character(logfile$summary[,2])
		logfile$summary[11,1]<<-"Replicate filter"
		logfile$summary[11,2]<<-"FALSE"	
		logfile$summary
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
	}

	if(!any(logfile$summary[,1]=="IS screening")){
		logfile$summary[,1]<<-as.character(logfile$summary[,1])
		logfile$summary[,2]<<-as.character(logfile$summary[,2])
		logfile$summary[12,1]<<-"IS screening"
		logfile$summary[12,2]<<-"FALSE"	
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
	}

	if(!any(logfile$summary[,1]=="Target screening")){
		logfile$summary[,1]<<-as.character(logfile$summary[,1])
		logfile$summary[,2]<<-as.character(logfile$summary[,2])
		logfile$summary[13,1]<<-"Target screening"
		logfile$summary[13,2]<<-"FALSE"	
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
	}

	if(!any(logfile$summary[,1]=="LOD")){
		logfile$summary[,1]<<-as.character(logfile$summary[,1])
		logfile$summary[,2]<<-as.character(logfile$summary[,2])
		logfile$summary[14,1]<<-"LOD"
		logfile$summary[14,2]<<-"FALSE"	
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
	}

	if(!any(logfile$summary[,1]=="quantification")){
		logfile$summary[,1]<<-as.character(logfile$summary[,1])
		logfile$summary[,2]<<-as.character(logfile$summary[,2])
		logfile$summary[15,1]<<-"quantification"
		logfile$summary[15,2]<<-"FALSE"	
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
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
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


	# logfile$Tasks_to_redo

	names(logfile[[2]])<<-c(
		"peakpick","QC","recal","normalize","allign","profiling","trendblind","pattern",
		"replicates","IS_screen","target_screen","LOD","quantification","blinds","norm_prof","-"
	)	
	 names(logfile)[2]<<-c("Tasks_to_redo"); 

	# logfile$workflow
	if(!any(names(logfile$workflow)=="screen_IS")){
		logfile$workflow[11]<<-"yes"; 	names(logfile$workflow)[11]<<-"screen_IS" 
	}
	if(!any(names(logfile$workflow)=="screen_target")){
		logfile$workflow[12]<<-"yes"; 	names(logfile$workflow)[12]<<-"screen_target" 
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

	logfile$workflow[5]<<-"yes"; 	names(logfile$workflow)[5]<<-"pattern" 
	logfile$workflow[7]<<-"yes"; 	names(logfile$workflow)[7]<<-"peakpicking" 		
	logfile$workflow[16]<<-"yes"; 	names(logfile$workflow)[16]<<-"-" 
	logfile$workflow[17]<<-"yes"; 	names(logfile$workflow)[17]<<-"-" 


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
			logfile$summary[logfile$summary[,1]=="Isotope pattern?",2]<<-"TRUE";
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="pattern"]<<-TRUE
		}
	}

	################################################################################################	
	logfile[[10]]<<-3.100
	names(logfile)[10]<<-"version"
	################################################################################################
	# define matrix of downstream workflow dependencies ############################################
	# requires only a definition of direct ones - inderect ones will be in workflow_set.r ##########
	# below specified in a row-wise fashion (but stored columnwise): ###############################
	# define workflow order of logfile$Tasks_to_redo by server.calculation.r #######################
	# dependencies must simply go after their parent node ########################################## 
	# order here actually irrelevant, because calculation order set in server_calculation  #########	
	work_order<-c(1,2,8,3,5,4,14,9,6,12,10,11,15,13,7)
	work_names<-(names(logfile$Tasks_to_redo)[work_order])
	depend<-matrix(ncol=length(work_names),nrow=length(work_names),0)
	colnames(depend)<-work_names
	rownames(depend)<-work_names					# peakpick	QC	pattern	recal	allign	normalize	replicates	profiling	IS_screen	target_screen	norm_prof	trendblind	LOD		quantification	blinds
	depend[,colnames(depend)=="peakpick"]<-			c(0,		1,	0,		1,		1,		1,			1,			1,			1,			1,				1,			1,			1,		1,				1)
	depend[,colnames(depend)=="QC"]<-				c(0,		0,	0,		1,		1,		0,			1,			1,			1,			1,				1,			1,			1,		1,				1)
	depend[,colnames(depend)=="pattern"]<-			c(0,		0,	0,		1,		1,		0,			0,			0,			1,			1,				1,			0,			0,		1,				0)
	depend[,colnames(depend)=="recal"]<-			c(0,		0,	0,		0,		1,		0,			1,			1,			1,			1,				1,			0,			0,		1,				1)
	depend[,colnames(depend)=="allign"]<-			c(0,		0,	0,		0,		0,		0,			1,			1,			1,			1,				1,			0,			0,		0,				1)
	depend[,colnames(depend)=="normalize"]<-		c(0,		0,	0,		0,		0,		0,			0,			1,			1,			1,				0,			1,			1,		1,				1)
	depend[,colnames(depend)=="blinds"]<-			c(0,		0,	0,		0,		0,		0,			1,			1,			1,			1,				1,			1,			1,		1,				1)
	depend[,colnames(depend)=="replicates"]<-		c(0,		0,	0,		0,		0,		0,			0,			1,			1,			1,				1,			1,			1,		1,				0)
	depend[,colnames(depend)=="profiling"]<-		c(0,		0,	0,		0,		0,		0,			0,			0,			1,			1,				1,			1,			0,		1,				0)
	depend[,colnames(depend)=="IS_screen"]<-		c(0,		0,	0,		0,		0,		0,			0,			0,			0,			1,				1,			0,			0,		1,				0)
	depend[,colnames(depend)=="target_screen"]<-	c(0,		0,	0,		0,		0,		0,			0,			0,			0,			0,				0,			0,			0,		1,				0)
	depend[,colnames(depend)=="norm_prof"]<-		c(0,		0,	0,		0,		0,		0,			0,			0,			0,			0,				0,			1,			0,		0,				0)
	depend[,colnames(depend)=="trendblind"]<-		c(0,		0,	0,		0,		0,		0,			0,			0,			0,			0,				0,			0,			0,		1,				0)
	depend[,colnames(depend)=="LOD"]<-				c(0,		0,	0,		0,		0,		0,			0,			0,			1,			1,				1,			1,			0,		1,				0)
	depend[,colnames(depend)=="quantification"]<-	c(0,		0,	0,		0,		0,		0,			0,			0,			0,			0,				0,			1,			0,		0,				0)	
	logfile[[11]]<<-depend
	names(logfile)[11]<<-"workflow_depend"
	################################################################################################
	# define upstream workflow "musts", i.e., upstream nodes on which`s execution a node ###########
	# depends. 0 = not dependent. 1 = dependent. -1 = MUST NOT be executed ######################### 
	must<-matrix(ncol=length(work_names),nrow=length(work_names),0)
	colnames(must)<-work_names
	rownames(must)<-work_names					# peakpick	QC	pattern	recal	allign	normalize	replicates	profiling	IS_screen	target_screen	norm_prof	trendblind	LOD		quantification	blinds
	must[,colnames(must)=="peakpick"]<-			c(0,		0,	0,		0,		0,		0,			0,			0,			0,			0,				0,			0,			0,		0,				0)
	must[,colnames(must)=="QC"]<-				c(1,		0,	0,		0,		0,		0,			0,			0,			0,			0,				0,			0,			0,		0,				0)
	must[,colnames(must)=="pattern"]<-			c(0,		0,	0,		0,		0,		0,			0,			0,			0,			0,				0,			0,			0,		0,				0)
	must[,colnames(must)=="recal"]<-			c(1,		0,	1,		0,		0,		0,			0,			0,			0,			0,				0,			0,			0,		0,				0)
	must[,colnames(must)=="allign"]<-			c(1,		0,	0,		0,		0,		0,			0,			0,			0,			0,				0,			0,			0,		0,				0)
	must[,colnames(must)=="normalize"]<-		c(1,		0,	0,		0,		0,		0,			0,			0,			0,			0,				0,			0,			0,		0,				0)
	must[,colnames(must)=="blinds"]<-			c(1,		0,	0,		0,		0,		0,			0,			0,			0,			0,				0,			0,			0,		0,				0)
	must[,colnames(must)=="replicates"]<-		c(1,		0,	0,		0,		0,		0,			0,			0,			0,			0,				0,			0,			0,		0,				0)
	must[,colnames(must)=="profiling"]<-		c(1,		0,	0,		0,		0,		0,			0,			0,			0,			0,				0,			0,			0,		0,				0)
	must[,colnames(must)=="IS_screen"]<-		c(1,		0,	1,		0,		0,		0,			0,			0,			0,			0,				0,			0,			0,		0,				0)
	must[,colnames(must)=="target_screen"]<-	c(1,		0,	0,		0,		0,		0,			0,			0,			0,			0,				0,			0,			0,		0,				0)
	must[,colnames(must)=="norm_prof"]<-		c(1,		0,	1,		0,		0,		0,			0,			1,			1,			0,				0,			0,			0,		0,				0)
	must[,colnames(must)=="trendblind"]<-		c(1,		0,	0,		0,		0,		0,			0,			1,			0,			0,				0,			0,			0,		0,				0)
	must[,colnames(must)=="LOD"]<-				c(1,		0,	0,		0,		0,		0,			0,			0,			0,			0,				0,			0,			0,		0,				0)
	must[,colnames(must)=="quantification"]<-	c(1,		0,	1,		0,		0,		0,			0,			0,			1,			1,				0,			0,			0,		0,				0)
	logfile[[12]]<<-must
	names(logfile)[12]<<-"workflow_must"	
	################################################################################################	
	save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
	load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv")) 

}






if(any(ls()=="logfile")){stop("\n illegal logfile detected #2 in server_updates.r!")}
