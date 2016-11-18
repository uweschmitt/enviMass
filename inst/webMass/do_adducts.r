# Run adduct grouping, filewise

    measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	cat("Adduct grouping: ")
	for(b in 1:length(measurements[,"ID"])){
		if( 
			(measurements[b,names(measurements)=="include"]=="TRUE") & 			# included?
			(measurements[b,names(measurements)=="adducts"]=="FALSE")  	# not yet done
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
			if( file.exists(file.path(logfile[[1]],"results","componentization","adducts",for_file) ) ){
				file.remove(file.path(logfile[[1]],"results","componentization","adducts",for_file) )
			}			
			# Peaklist
			load(file=file.path(logfile[[1]],"peaklist",as.character(for_file)));   
			peaklist<-peaklist[order(peaklist[,10],decreasing=FALSE),] # match with IDs - for saving pattern; IDs are retrieved for pairs seperately
			# EIC pairs
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
						cat(" with exclusion - ")
					}
				}else{
					exclude<-FALSE
				}
			}else{
				exclude<-FALSE
			}
			##########################################################################	
			cat("grouping - ")			
			if(measurements[b,names(measurements)=="Mode"]=="positive"){
				with_adducts<-logfile$parameters$external$adducts_pos
				with_mode<-"positive"
			}else{
				with_adducts<-logfile$parameters$external$adducts_neg
				with_mode<-"negative"
			}
			if(length(with_adducts)<2){
				cat("\n Not enough adducts for this ionization mode specified - skipped ...")
				next;
			}			
			peaklist2<-as.data.frame(peaklist[,c(12,13,14)])	
			relat<-try(
				enviMass:::adduct.search2( # dont name it "adducts" -> conflict
					peaklist2, 
					adducts, 
					rttol = logfile$parameters$external$adducts_rttol, 
					mztol = logfile$parameters$external$adducts_mztol,
					ppm = logfile$parameters$external$adducts_ppm, 
					use_adducts = with_adducts, 
					ion_mode = with_mode,
					exclude
				)
			)
			if(class(relat)=="try-error"){
				cat("\n Adduct detection failed - adpat parameters?");
				next;
			}				
			if(length(relat[,1])==0){
				cat("\n No adduct relations detected");
				next;
			}
			Adduct_pairs<-relat[,c(1,2)]
			those<-(Adduct_pairs[,1]>Adduct_pairs[,2])
			if(any(those)){
				Adduct_pairs[those,]<-Adduct_pairs[those,c(2,1)]
			}
			Adduct_pairs<-Adduct_pairs[order(Adduct_pairs[,1],Adduct_pairs[,2],decreasing=FALSE),]
			Adduct_pairs[,1]<-peaklist[Adduct_pairs[,1],"peak_ID"]	# although sorted above - just to be save
			Adduct_pairs[,2]<-peaklist[Adduct_pairs[,2],"peak_ID"]	# although sorted above - just to be save				
			save(Adduct_pairs,file=(file.path(logfile[[1]],"results","componentization","adducts",paste(for_file,sep=""))))
			rm(peaklist,peaklist2,those,Adduct_pairs,with_adducts,with_mode,for_file)
			##########################################################################	
			measurements[b,names(measurements)=="adducts"]<-"TRUE"
			write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
			cat("done.")
			##########################################################################		
		}else{
			cat("\n Adducts extracted before.")
		}
	}
	rm(measurements)
	