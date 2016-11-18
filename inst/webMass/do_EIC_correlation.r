# Run EIC correlation, filewise

    measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	EICor_plot_it<-FALSE # plot intermediate EIC pairs
	do_cor<-TRUE		# omit correlation
	cat("EIC correlations: ")
	for(b in 1:length(measurements[,"ID"])){
		if( 
			(measurements[b,names(measurements)=="include"]=="TRUE") & 			# included?
			(measurements[b,names(measurements)=="EIC_correlation"]=="FALSE")  	# not yet done
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
			for_file<-measurements[b,1]
			if( file.exists(file.path(logfile[[1]],"results","componentization","EIC_corr",for_file) ) ){
				file.remove(file.path(logfile[[1]],"results","componentization","EIC_corr",for_file) )
			}			
			# Peaklist
			load(file=file.path(logfile[[1]],"peaklist",as.character(for_file)));   
			# MSlist
			load(file=file.path(logfile[[1]],"MSlist",as.character(for_file)));   
			##########################################################################	
			# EXTRACT CANDIDATE PAIRS ################################################
			cat("pairing - ")
			ord<-order(peaklist[,14],decreasing=FALSE)
			peaklist<-peaklist[ord,,drop=FALSE]
			#use_peak<-MSlist[[7]][peaklist[,10],3]>=1#logfile$parameters$external$EICor_minpeaks # too few data points - included above; but interferes with do_componentization.r strategy!
			at<-1
			paired<-matrix(ncol=4,nrow=1E7,0)
			len<-dim(paired)[1]
			for(i in 1:(length(peaklist[,1])-1)){
				#if(!use_peak[i]){next} # enough data points?
				for(j in (i+1):length(peaklist[,1])){
					#if(!use_peak[j]){next} # enough data points?
					RTdif<-(peaklist[j,14]-peaklist[i,14])
					if(RTdif<=logfile$parameters$external$EICor_delRT){
						paired[at,1]<-i
						paired[at,2]<-j
						paired[at,3]<-RTdif		
						at<-(at+1)
						if(at>len){
							cat(".")
							paired2<-matrix(ncol=4,nrow=1E7,0)
							paired<-rbind(paired,paired2)
							len<-dim(paired)[1]
						}
					}else{
						break;
					}
				}
			}
			paired<-paired[1:(at-1),,drop=FALSE]
			if(length(paired[,1])==0){
				cat("nothing found; aborted.");next;
			}
			###########################################################################
			# EIC CORRELATION #########################################################
			cat("correlating - ")
			skipped<-0
			for(i in 1:length(paired[,1])){			
				if(i>1){ # no reload
					if(paired[i,1]!=paired[i-1,1]){
						PEAK_ID_1<-peaklist[paired[i,1],"peak_ID"]
						start_1<-MSlist[[7]][PEAK_ID_1,1]
						end_1<-MSlist[[7]][PEAK_ID_1,2]
						del_1<-MSlist[[7]][PEAK_ID_1,3]
					}	
				}else{
					PEAK_ID_1<-peaklist[paired[i,1],10]
					start_1<-MSlist[[7]][PEAK_ID_1,1]
					end_1<-MSlist[[7]][PEAK_ID_1,2]
					del_1<-MSlist[[7]][PEAK_ID_1,3]
				}
				if(del_1<logfile$parameters$external$EICor_minpeaks){ # not enough EIC signals, to be removed later
					paired[i,4]<-1000
					next;
				}
				PEAK_ID_2<-peaklist[paired[i,2],10]
				start_2<-MSlist[[7]][PEAK_ID_2,1]
				end_2<-MSlist[[7]][PEAK_ID_2,2]
				del_2<-MSlist[[7]][PEAK_ID_2,3]
				if(del_2<logfile$parameters$external$EICor_minpeaks){ # not enough EIC signals, to be removed later
					paired[i,4]<-1000
					next;
				}				
				####################################################################
				# extract data, along shorter data set & correlate #################
				if(del_1<=del_2){
					those<-match(MSlist[[4]][[2]][start_1:end_1,3],MSlist[[4]][[2]][start_2:end_2,3])
					if(sum(!is.na(those))<logfile$parameters$external$EICor_minpeaks){ # not enough EIC signal OVERLAP
						#paired[i,4]<-(-2000)
						skipped<-(skipped+1);
						next;
					}
					if(do_cor){		
						paired[i,4]<-cor(
							MSlist[[4]][[2]][start_1:end_1,2][!is.na(those)],
							MSlist[[4]][[2]][start_2:end_2,2][those[!is.na(those)]],
							method="spearman"
						)			
					}
				}else{
					those<-match(MSlist[[4]][[2]][start_2:end_2,3],MSlist[[4]][[2]][start_1:end_1,3])
					if(sum(!is.na(those))<logfile$parameters$external$EICor_minpeaks){ # not enough EIC signal OVERLAP
						#paired[i,4]<-(-2000)						
						skipped<-(skipped+1);
						next;
					}
					if(do_cor){			
						paired[i,4]<-cor(
							MSlist[[4]][[2]][start_2:end_2,2][!is.na(those)],
							MSlist[[4]][[2]][start_1:end_1,2][those[!is.na(those)]],
							method="spearman"
						)
					}
				}
				####################################################################
				# plotit ! #########################################################
				if(EICor_plot_it){
					if(paired[i,4]>=EICor_mincor){
						split.screen(c(2,1))
						xlim<-c(
							min(c(MSlist[[4]][[2]][start_1:end_1,3],MSlist[[4]][[2]][start_2:end_2,3])),
							max(c(MSlist[[4]][[2]][start_1:end_1,3],MSlist[[4]][[2]][start_2:end_2,3]))
						)
						par(mar=c(3.5,3.5,.1,.1))
						screen(1)
							plot(
								MSlist[[4]][[2]][start_1:end_1,3],
								MSlist[[4]][[2]][start_1:end_1,2],
								type="h",xlim=xlim,ylab="Intensity",xlab=""
							)
						screen(2)
						par(mar=c(3.5,3.5,.1,.1))
							plot(
								MSlist[[4]][[2]][start_2:end_2,3],
								MSlist[[4]][[2]][start_2:end_2,2],
								type="h",xlim=xlim,ylab="Intensity",xlab="RT"
							)
						close.screen(n, all.screens = TRUE)
						Sys.sleep(.3)
					}
				}
				#####################################################################
			}
			##########################################################################	
			# FILTER & SAVE RESULTS ##################################################
			cat("filtering: ")
			EIC_pairs<-paired#[paired[,4]>=logfile$parameters$external$EICor_mincor,,drop=FALSE]
			EIC_pairs[,1]<-peaklist[EIC_pairs[,1],10] # insert peak ID! peaklist has been resorted here!
			EIC_pairs[,2]<-peaklist[EIC_pairs[,2],10] # insert peak ID! peaklist has been resorted here!			
			#EIC_pairs<-EIC_pairs[EIC_pairs[,4]!=1000,,drop=FALSE] # no! dont drop them; artefacts in co_componentization otherwise!
			if(length(EIC_pairs[,1])==0){
				cat("nothing found; aborted.");next;
			}
			# re-arrange: values increasing per row and towards bottoms over columns 
			those<-(EIC_pairs[,1]>EIC_pairs[,2])
			if(any(those)){
				EIC_pairs[those,]<-EIC_pairs[those,c(2,1,3,4),drop=FALSE]
			}
			EIC_pairs<-EIC_pairs[order(EIC_pairs[,1],EIC_pairs[,2],decreasing=FALSE),,drop=FALSE]
			save(EIC_pairs,file=file.path(logfile[[1]],"results","componentization","EIC_corr",for_file))
			cat(as.character(round((sum(EIC_pairs[,4]>=logfile$parameters$external$EICor_mincor)/length(paired[,1])*100),digits=2)));
			cat("% = ");
			cat(as.character(length(EIC_pairs[,1])));
			cat(" pairs corrrelated - ");
			rm(paired,EIC_pairs,peaklist,MSlist);
			##########################################################################	
			measurements[b,names(measurements)=="EIC_correlation"]<-"TRUE"
			write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
			cat("done.")
			##########################################################################		
		}else{
			cat("\n EIC correlation done before.")
		}
	}
	rm(measurements)

