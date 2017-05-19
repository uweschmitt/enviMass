
	########################################################################################
	# remove SOME old results ##############################################################
	those<-list(0)
	those[[1]]<-file.path(logfile[[1]],"quantification","target_quant_table_pos")		
	those[[2]]<-file.path(logfile[[1]],"quantification","target_quant_table_pos_warn")				
	those[[3]]<-file.path(logfile[[1]],"quantification","target_quant_table_neg")			
	those[[4]]<-file.path(logfile[[1]],"quantification","target_quant_table_neg_warn")	
	for(n in 1:length(those)){
		if(file.exists(those[[n]])){
			file.remove(those[[n]])
		}
	}	
	rm(those)
	########################################################################################
	





