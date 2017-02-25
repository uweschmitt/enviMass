
	#######################################################################################################################
	# delete old results ##################################################################################################
	those<-list(0)
	those[[1]]<-file.path(logfile[[1]],"quantification","target_recov_table_pos")			
	those[[2]]<-file.path(logfile[[1]],"quantification","target_recov_table_neg")	
	for(n in 1:length(those)){
		if(file.exists(those[[n]])){
			file.remove(those[[n]])
		}
	}	
	rm(those)
	measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");		
	######################################################################################################################
	
