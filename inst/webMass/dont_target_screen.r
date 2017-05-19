########################################################################################################	
# clean ! ##############################################################################################
those<-list.files(file.path(logfile$project_folder,"results","screening"))
if(length(those)>0){
	for(i in 1:length(those)){
		if(grepl("target",those[i])){ # distinguish from IS
			file.remove(file.path(logfile$project_folder,"results","screening",those[i]))
		}
	}
}
########################################################################################################	
	







