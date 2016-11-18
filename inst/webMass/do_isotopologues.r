######################################################################################################################
# Run isotopologue grouping, filewise ################################################################################
###################################################################################################################### 	

    measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	cat("Isotopologue grouping: ")
	load(file.path(logfile[[1]],"results","componentization","isotopologues","quantiz"))
	#save(quantiz,file=file.path(logfile[[1]],"results","componentization","isotopologues","quantiz"))
	if(quantiz[[1]]!=logfile$parameters$resolution){
		cat("\n WARNING: seems the quantized data for isotopologue grouping does NOT MATCH your selected resolution! Please resolve this issue.")
	}
	for(b in 1:length(measurements[,"ID"])){
		if( 
			(measurements[b,names(measurements)=="include"]=="TRUE") & 			# included?
			(measurements[b,names(measurements)=="isotopologues"]=="FALSE")  	# not yet done
		){ 
		
			if(
				(measurements[b,names(measurements)=="profiled"]=="FALSE") & 	
				(logfile$parameters$prof_select=="TRUE")
			){	next }	
			##########################################################################
			# exclude files that do not end up in profiles ###########################
			if( (mute(logfile$parameters$prof_select=="TRUE")) & (measurements$profiled[b]=="FALSE") ){
				cat("\n Skip file - not included in profile building.");next;
			}
			##########################################################################
			cat(paste("\n Doing ",as.character(b)," of ",as.character(length(measurements[,"ID"]))," files: ",sep=""))		
			##########################################################################
			# LOAD FILES & REMOVE OLD RESULTS ########################################
			cat("loading - ")
			for_file<-measurements[b,"ID"]
			if( file.exists(file.path(logfile[[1]],"results","componentization","isotopologues",for_file) ) ){
				file.remove(file.path(logfile[[1]],"results","componentization","isotopologues",for_file) )
			}			
			load(file=file.path(logfile[[1]],"peaklist",as.character(for_file))); # Peaklist  
			peaklist<-peaklist[order(peaklist[,10],decreasing=FALSE),] # match with IDs - for saving pattern; IDs are retrieved for pairs seperately
			if((logfile$workflow[names(logfile$workflow)=="EIC_correlation"]=="yes") & FALSE){ # load EIC correlation results - removed, also in EIC -> isot. depends matrix!
				if(file.exists(file.path(logfile[[1]],"results","componentization","EIC_corr",for_file))){
					load(file.path(logfile[[1]],"results","componentization","EIC_corr",for_file))
					exclude<-EIC_pairs[
						EIC_pairs[,4]<logfile$parameters$external$EICor_mincor # ==2000 indicates not enough overlap, despite enough EIC signals present
					,1:2,drop=FALSE]
					rm(EIC_pairs)
					if(length(exclude[,1])==0){
						exclude<-FALSE
					}else{
						cat(" with exclusion -")
					}
				}else{
					exclude<-FALSE
				}
			}else{
				exclude<-FALSE
			}
			##########################################################################	
			cat("grouping - ")
			peaklist2<-as.data.frame(peaklist[,c(12,13,14)])				
			#pattern<-enviMass:::pattern.search3(
			relat<-try(
				enviMass:::pattern.search3(
					peaklist2,
					quantiz,
					mztol=logfile$parameters$external$isotop_mztol,
					ppm=logfile$parameters$external$isotop_ppm,
					inttol=logfile$parameters$external$isotop_inttol,
					rttol=logfile$parameters$external$isotop_rttol,
					use_isotopes=FALSE,
					use_charges=logfile$parameters$external$isotop_use_charges,
					use_marker=TRUE,
					quick=TRUE,
					isotopes,
					exclude
				)
			)
			if(class(relat)=="try-error"){
				cat("\n Isotopologue detection failed - adpat parameters?");
				next;
			}				
			if(length(relat[,1])==0){
				cat("\n No adduct relations detected");
				next;
			}
			Isot_pairs<-relat
			those<-(Isot_pairs[,1]>Isot_pairs[,2])
			if(any(those)){
				Isot_pairs[those,]<-Isot_pairs[those,c(2,1)]
			}
			Isot_pairs<-Isot_pairs[order(Isot_pairs[,1],Isot_pairs[,2],decreasing=FALSE),]
			# insert peak IDs - peaklist may have been reordered before! ##############
			Isot_pairs[,1]<-peaklist[Isot_pairs[,1],"peak_ID"]
			Isot_pairs[,2]<-peaklist[Isot_pairs[,2],"peak_ID"]				
			save(Isot_pairs,file=(file.path(logfile[[1]],"results","componentization","isotopologues",paste(for_file,sep="_"))))
			##########################################################################	
			measurements[b,names(measurements)=="isotopologues"]<-"TRUE"
			write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
			rm(peaklist,peaklist2)			
			cat("done.")
			##########################################################################		
		}else{
			cat("\n Isotopologues detected before.")
		}
	}
	rm(quantiz,measurements)

