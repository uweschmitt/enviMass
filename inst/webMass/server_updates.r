# include updates - if older projects are reopened

if(!file.exists(file.path(logfile$project_folder,"results","screening"))
){
	dir.create(file.path(logfile$project_folder,"results","screening"),recursive=TRUE)    	# subfolder  
}

if(!file.exists(file.path(logfile$project_folder,"results","quantification"))
){	
	dir.create(file.path(logfile$project_folder,"results","quantification"),recursive=TRUE)   # subfolder 
}





