# include updates - if older projects are reopened

if(!file.exists(file.path(logfile$project_folder,"results","screening"))
){
	dir.create(file.path(logfile$project_folder,"results","screening"),recursive=TRUE)    	# subfolder  
}

if(!file.exists(file.path(logfile$project_folder,"results","quantification"))
){	
	dir.create(file.path(logfile$project_folder,"results","quantification"),recursive=TRUE)   # subfolder 
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
	logfile$summary[,1]<-as.character(logfile$summary[,1])
	logfile$summary[,2]<-as.character(logfile$summary[,2])
	logfile$summary[11,1]<-"Replicate filter"
	logfile$summary[11,2]<-"FALSE"	
	logfile$summary
	save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
	load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}

if(!any(logfile$summary[,1]=="IS screening")){
	logfile$summary[,1]<-as.character(logfile$summary[,1])
	logfile$summary[,2]<-as.character(logfile$summary[,2])
	logfile$summary[12,1]<-"IS screening"
	logfile$summary[12,2]<-"FALSE"	
	save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
	load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}

if(!any(logfile$summary[,1]=="Target screening")){
	logfile$summary[,1]<-as.character(logfile$summary[,1])
	logfile$summary[,2]<-as.character(logfile$summary[,2])
	logfile$summary[13,1]<-"Target screening"
	logfile$summary[13,2]<-"FALSE"	
	save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
	load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}













