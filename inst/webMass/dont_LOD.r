####################################################################
# delete old LOD gams & results ####################################
those<-list.files(file.path(logfile$project_folder,"results","LOD"))
if(length(those)>0){
	for(i in 1:length(those)){
		file.remove(file.path(logfile$project_folder,"results","LOD",those[i]))
	}
}

measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
measurements[,16]<-FALSE;
write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);




