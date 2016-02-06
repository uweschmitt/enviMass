# include updates - if older projects are reopened
if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_updates.r!")}

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
	"replicates","IS_screen","target_screen","LOD","quantification","-","norm_prof","-"
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
logfile$workflow[5]<<-"yes"; 	names(logfile$workflow)[5]<<-"pattern" 
logfile$workflow[7]<<-"yes"; 	names(logfile$workflow)[7]<<-"peakpicking" 
logfile$workflow[14]<<-"yes"; 	names(logfile$workflow)[14]<<-"-" 		
logfile$workflow[16]<<-"yes"; 	names(logfile$workflow)[16]<<-"-" 
logfile$workflow[17]<<-"yes"; 	names(logfile$workflow)[17]<<-"-" 

save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv")) 

if(any(ls()=="logfile")){stop("\n illegal logfile detected #2 in server_updates.r!")}



