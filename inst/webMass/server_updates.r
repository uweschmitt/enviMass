# include updates - if older projects are reopened
if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_updates.r!")}




#if(logfile[[10]]<3.100){	
if(TRUE){
	
	cat("\n Updating to version 3.1 ...")
	################################################################################################
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
	if(!any(logfile$summary[,1]=="replicates")){
		logfile$summary[,1]<<-as.character(logfile$summary[,1])
		logfile$summary[,2]<<-as.character(logfile$summary[,2])
		logfile$summary[11,1]<<-"replicates"
		logfile$summary[11,2]<<-"FALSE"	
		logfile$summary
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
	}

	if(!any(logfile$summary[,1]=="IS_screen")){
		logfile$summary[,1]<<-as.character(logfile$summary[,1])
		logfile$summary[,2]<<-as.character(logfile$summary[,2])
		logfile$summary[12,1]<<-"IS_screen"
		logfile$summary[12,2]<<-"FALSE"	
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
	}

	if(!any(logfile$summary[,1]=="target_screen")){
		logfile$summary[,1]<<-as.character(logfile$summary[,1])
		logfile$summary[,2]<<-as.character(logfile$summary[,2])
		logfile$summary[13,1]<<-"target_screen"
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


	# logfile$Tasks_to_redo ##################################################################
	names(logfile[[2]])<<-c(
		"peakpicking","qc","recal","norm","align","profiling","trendblind","pattern",
		"replicates","IS_screen","target_screen","LOD","quantification","blinds","IS_normaliz","-"
	)	
	 names(logfile)[2]<<-c("Tasks_to_redo"); 

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
	logfile$workflow[16]<<-"yes"; 	names(logfile$workflow)[16]<<-"-" 
	logfile$workflow[17]<<-"yes"; 	names(logfile$workflow)[17]<<-"-" 

	# logfile$summary #################################################################
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
		""
	 )

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
	logfile[[10]]<<-3.100
	names(logfile)[10]<<-"version"
	################################################################################################
	# define matrix of downstream workflow dependencies ############################################
	# requires only a definition of direct ones - indirect ones will be in workflow_set.r ##########
	# below specified in a row-wise fashion (but stored and retrieved columnwise): #################
	# define workflow order of logfile$Tasks_to_redo by server.calculation.r #######################
	# dependencies must simply go after their parent node ########################################## 
	# order here actually irrelevant, because calculation order set in server_calculation  #########	
	work_names<-names(logfile$Tasks_to_redo)[1:15]
	depend<-matrix(ncol=length(work_names),nrow=length(work_names),0)
	colnames(depend)<-work_names
	rownames(depend)<-work_names					# peakpicking	qc	recal	norm	align	profiling	trendblind	pattern		replicates	IS_screen	target_screeen	LOD		quantification	blinds IS_normaliz
	depend[,colnames(depend)=="peakpicking"]<-		c(0,			1,	1,		1,		1,		1,			1,			1,			1,			1,			1,				1,		1,				1,		1)
	depend[,colnames(depend)=="qc"]<-				c(0,			0,	1,		1,		1,		1,			1,			0,			1,			1,			1,				1,		1,				1,		1)
	depend[,colnames(depend)=="pattern"]<-			c(0,			0,	1,		0,		0,		0,			0,			0,			0,			1,			1,				0,		1,				0,		1)
	depend[,colnames(depend)=="recal"]<-			c(0,			0,	0,		0,		0,		1,			1,			0,			1,			1,			1,				0,		1,				0,		1)
	depend[,colnames(depend)=="align"]<-			c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,				0,		0)
	depend[,colnames(depend)=="norm"]<-				c(0,			0,	0,		0,		0,		1,			1,			0,			0,			1,			1,				1,		1,				1,		1)
	depend[,colnames(depend)=="blinds"]<-			c(0,			0,	0,		0,		0,		1,			1,			0,			0,			1,			1,				0,		1,				0,		1)
	depend[,colnames(depend)=="replicates"]<-		c(0,			0,	0,		0,		0,		1,			1,			0,			0,			1,			1,				1,		1,				0,		1)
	depend[,colnames(depend)=="profiling"]<-		c(0,			0,	0,		0,		0,		0,			1,			0,			0,			1,			1,				0,		1,				0,		1)
	depend[,colnames(depend)=="IS_screen"]<-		c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		1,				0,		1)
	depend[,colnames(depend)=="target_screen"]<-	c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		1,				0,		1)
	depend[,colnames(depend)=="IS_normaliz"]<-		c(0,			0,	0,		0,		0,		0,			1,			0,			0,			0,			0,				0,		0,				0,		1)
	depend[,colnames(depend)=="trendblind"]<-		c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,				0,		0)
	depend[,colnames(depend)=="LOD"]<-				c(0,			0,	0,		0,		0,		0,			0,			0,			0,			1,			1,				0,		1,				0,		1)
	depend[,colnames(depend)=="quantification"]<-	c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,				0,		0)	
	logfile[[11]]<<-depend
	names(logfile)[11]<<-"workflow_depend"
	################################################################################################
	# define upstream workflow "musts", i.e., upstream nodes on which`s execution a node ###########
	# depends. 0 = not dependent. 1 = dependent. -1 = MUST NOT be executed ######################### 
	must<-matrix(ncol=length(work_names),nrow=length(work_names),0)
	colnames(must)<-work_names
	rownames(must)<-work_names					# peakpicking	qc	recal	norm	align	profiling	trendblind	pattern		replicates	IS_screen	target_screeen	LOD		quantification	blinds IS_normaliz
	must[,colnames(must)=="peakpicking"]<-		c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,				0,		0)
	must[,colnames(must)=="qc"]<-				c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,				0,		0)
	must[,colnames(must)=="pattern"]<-			c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,				0,		0)
	must[,colnames(must)=="recal"]<-			c(1,			0,	0,		0,		0,		0,			0,			1,			0,			0,			0,				0,		0,				0,		0)
	must[,colnames(must)=="align"]<-			c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,				0,		0)
	must[,colnames(must)=="norm"]<-				c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,				0,		0)
	must[,colnames(must)=="blinds"]<-			c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,				0,		0)
	must[,colnames(must)=="replicates"]<-		c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,				0,		0)
	must[,colnames(must)=="profiling"]<-		c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,				0,		0)
	must[,colnames(must)=="IS_screen"]<-		c(1,			0,	0,		0,		0,		0,			0,			1,			0,			0,			0,				0,		0,				0,		0)
	must[,colnames(must)=="target_screen"]<-	c(1,			0,	0,		0,		0,		0,			0,			1,			0,			0,			0,				0,		0,				0,		0)
	must[,colnames(must)=="IS_normaliz"]<-		c(1,			0,	0,		0,		0,		0,			0,			1,			0,			0,			0,				0,		0,				0,		0)
	must[,colnames(must)=="trendblind"]<-		c(1,			0,	0,		0,		0,		1,			0,			0,			0,			0,			0,				0,		0,				0,		0)
	must[,colnames(must)=="LOD"]<-				c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,				0,		0)
	must[,colnames(must)=="quantification"]<-	c(1,			0,	0,		0,		0,		0,			0,			1,			0,			1,			1,				0,		0,				0,		0)
	logfile[[12]]<<-must
	names(logfile)[12]<<-"workflow_must"	
	################################################################################################	
	save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
	load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv")) 

}

if(logfile[[10]]<3.101){	

	cat("\n Updating to version 3.101 ...")
	################################################################################################
	# create missing folder
	if(!file.exists(file.path(logfile$project_folder,"quantification"))
	){
		dir.create(file.path(logfile$project_folder,"quantification"),recursive=TRUE)    	# subfolder  
	}
	################################################################################################	
	# updating columns in IS compound table ########################################################
 	intstand<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
	if(!any(names(intstand)=="Lower intensity bound")){
		names_1<-names(intstand)
		intstand<-cbind(intstand,rep(-Inf,length(intstand[,1])),rep(Inf,length(intstand[,1])))
		names(intstand)<-c(names_1,"Lower intensity bound","Upper intensity bound")
		write.table(intstand,file=file.path(logfile[[1]],"dataframes","IS.txt"),row.names=FALSE,sep="\t",quote=FALSE) 
	}
	rm(intstand)
	################################################################################################	
	# updating columns in targets compound table ################################################### 	
 	targets<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");	
	if(!any(names(targets)=="intensity_warn_1")){	
		names_1<-names(targets)
		targets<-cbind(targets,rep("FALSE",length(targets[,1])),rep("FALSE",length(targets[,1])))		
		names(targets)<-c(names_1,"intensity_warn_1","intensity_warn_2")
		write.table(targets,file=file.path(logfile[[1]],"dataframes","targets.txt"),row.names=FALSE,sep="\t",quote=FALSE)	
	}
	rm(targets)
	################################################################################################	
	logfile[[10]]<<-3.101
	names(logfile)[10]<<-"version"
	save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
	load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv")) 
}

if(any(ls()=="logfile")){stop("\n illegal logfile detected #2 in server_updates.r!")}
