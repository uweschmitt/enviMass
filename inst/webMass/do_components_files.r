###############################################################################################
# Run a file-wise componentization ############################################################
###############################################################################################

	###########################################################################################
	do_isot<-(logfile$workflow[names(logfile$workflow)=="isotopologues"]=="yes")
	do_addu<-(logfile$workflow[names(logfile$workflow)=="adducts"]=="yes")
	do_homol<-(logfile$workflow[names(logfile$workflow)=="homologues"]=="yes")
	if( do_isot | do_addu ){ # homol alone not sufficient to run nontarget::combine
		measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		for(b in 1:length(measurements[,"ID"])){
			if( 
				(measurements[b,"include"]=="TRUE") & 			# included?
				(measurements[b,"components_files"]=="FALSE")  	# not yet done
			){ 
			
				##########################################################################
				# exclude files that do not end up in profiles ###########################
				if( (mute(logfile$parameters$prof_select=="TRUE")) & (measurements$profiled[b]=="FALSE") ){
					cat("\n Skip file - not included in profile building.");next;
				}
				##########################################################################
				cat(paste("\n Doing ",as.character(b)," of ",as.character(length(measurements[,"ID"]))," files: ",sep=""))	
				for_file<-measurements[b,"ID"]
				##########################################################################
				# get isotopologue grouping results ######################################
				if(
					do_isot & file.exists(file.path(logfile[[1]],"results","componentization","isotopologues",paste("full",for_file,sep="_")))
				){
					load(file.path(logfile[[1]],"results","componentization","isotopologues",paste("full",for_file,sep="_")))
					cat("load isot. pattern - ")
				}else{
					pattern<-FALSE
				}
				##########################################################################
				# get adduct grouping results ############################################
				if(
					do_addu & file.exists(file.path(logfile[[1]],"results","componentization","adducts",paste("full",for_file,sep="_")))
				){
					load(file.path(logfile[[1]],"results","componentization","adducts",paste("full",for_file,sep="_")))
					cat("load adduct groups - ")
				}else{
					adduct<-FALSE
				}
				##########################################################################
				# get adduct grouping results ############################################
				if(
					do_homol & file.exists(file.path(logfile[[1]],"results","componentization","homologues",paste("full",for_file,sep="_")))
				){
					load(file.path(logfile[[1]],"results","componentization","homologues",paste("full",for_file,sep="_")))
					cat("load homologues - ")
				}else{
					homol<-FALSE
				}
				##########################################################################
				# build components #######################################################
				cat("combine: ")
				component<-nontarget::combine(
					pattern, 
					adduct, 
					homol, 
					rules = c(FALSE, FALSE, FALSE), 
					dont = FALSE)	
				save(component,file=(file.path(logfile[[1]],"results","componentization","components",paste(for_file))))
				##########################################################################	
				rm(pattern,adduct,homol,component)
				measurements[b,"components_files"]<-"TRUE"
				write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
				##########################################################################	
				#cat("done.")
			}
		}
		rm(measurements)
	}else{
		# clean component folder ###############################################################
		
# BAUSTELLE
		
	}
	############################################################################################



	

