if(
	(logfile$workflow[15]=="yes" && logfile$summary[9,2]=="FALSE") || (logfile$Tasks_to_redo[15]=="TRUE") 
){

    measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	if( (any(measurements[,4]=="positive") & file.exists(file.path(as.character(logfile[[1]]),"results","pattern_pos_IS"))) ){
		if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="peaklist")){rm(peaklist)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_pos")){rm(profileList_pos,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="profileList_pos")){rm(profileList_pos)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="pattern_pos_IS")){rm(pattern_pos_IS,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="pattern_pos_IS")){rm(pattern_pos_IS)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="patternRT_pos_IS")){rm(patternRT_pos_IS,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="patternRT_pos_IS")){rm(patternRT_pos_IS)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="patternDelRT_pos_IS")){rm(patternDelRT_pos_IS,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="patternDelRT_pos_IS")){rm(patternDelRT_pos_IS)}
		load(file=file.path(as.character(logfile[[1]]),"results","profileList_pos"),envir=as.environment(".GlobalEnv"));	
		load(file=file.path(logfile[[1]],"results","pattern_pos_IS"),envir=as.environment(".GlobalEnv"));
		load(file=file.path(logfile[[1]],"results","patternRT_pos_IS"),envir=as.environment(".GlobalEnv"));
		load(file=file.path(logfile[[1]],"results","patternDelRT_pos_IS"),envir=as.environment(".GlobalEnv"));
		peaks<-profileList_pos[[7]];
		peaklist<-peaks[,c(14,16,15)];
		treshscore<-as.numeric(logfile$parameters[78])
		we1=(1-as.numeric(logfile$parameters[49]))
		we2=(1-we1)
		cat("- screening:");
		screened<-screening(	peaklist, 
								blanklist=FALSE, 
								pattern_pos_IS, 
								RT = patternRT_pos_IS,
								dmz = as.numeric(logfile$parameters[45]), 
								ppm = logfile$parameters[46], 
								dRT = patternDelRT_pos_IS, 
								dRTwithin = as.numeric(logfile$parameters[43]), 
								dRTblank = FALSE, 
								dInt = as.numeric(logfile$parameters[47]), 
								Intcut = as.numeric(logfile$parameters[48]), 
								w1=we1, 
								w2=we2, 	
								w3=0	
		)
		# set matrix to sort & store data from a profile ###################################
		atPOSIX<-profileList_pos[[3]];
		sampletype<-profileList_pos[[9]];
		sampleID<-profileList_pos[[4]];
		atdate<-c();
		attime<-c();
		for(i in 1:length(atPOSIX)){
				  atdate<-c(atdate, strsplit(atPOSIX[i]," ")[[1]][1]);
				  attime<-c(attime, strsplit(atPOSIX[i]," ")[[1]][2]);
		}
		attime<-as.difftime(attime);
		atdate<-as.Date(atdate);
		ord<-order(as.numeric(atdate),as.numeric(attime),sampleID);
		atPOSIXsort<-atPOSIX[ord];
		atdate<-atdate[ord];
		attime<-attime[ord];
		sampleID<-sampleID[ord];
		sampletype<-sampletype[ord];
		# to retrieve data
		timeset<-matrix(nrow=length(atPOSIX),ncol=5,0);
		for(i in 1:length(sampleID)){
			  if(sampletype[i]=="sample"){
				timeset[i,2]<-as.numeric(sampleID[i]);
			  }
			  if(sampletype[i]=="blank"){
				timeset[i,3]<-as.numeric(sampleID[i]);
			  }
		}
		# screen IS intensity profiles #####################################################
		cat("- on IS profiles");
		min_count<-floor(length(profileList_pos[[4]])*as.numeric(logfile$parameters[[70]])/100);
		lis_delint_IS<-list(0)
		lis_median_IS<-list(0)
		stillin<-rep(TRUE,length(peaks[,1]))
		for(p in 1:length(timeset[,1])){
			lis_delint_IS[[p]]<-numeric(0)
			lis_median_IS[[p]]<-numeric(0)
		}
		for(i in 1:length(screened)){
			if(length(screened[[i]])!=1){
				for(j in 1:length(screened[[i]][,1])){
					if(  !(as.character(screened[[i]][j,9]))=="NaN"  ){
						if( as.numeric(as.character(screened[[i]][j,9]))>=treshscore ){	
							hits<-strsplit(as.character(screened[[i]][j,2]),"/")[[1]];
							hits<-hits[hits!="-"];
							hits<-as.numeric(hits);
							for(b in 1:length(hits)){
								if( peaks[hits[b],3]>min_count ){
									profID<-as.numeric(peaks[hits[b],4])
									timeset[,4:5]<-0;
									timeset[,c(4,5)] <-.Call("fill_timeset",
														as.numeric(timeset),
														as.numeric(profileList_pos[[2]][(profileList_pos[[7]][profileList_pos[[7]][,4]==profID,1]:profileList_pos[[7]][profileList_pos[[7]][,4]==profID,2]),6]), # sampleIDs
														as.numeric(profileList_pos[[2]][(profileList_pos[[7]][profileList_pos[[7]][,4]==profID,1]:profileList_pos[[7]][profileList_pos[[7]][,4]==profID,2]),2]), # intensities
														as.integer(length(timeset[,1])),
														PACKAGE="enviMass"
									)
									timeset[timeset[,4]!=0,4]<-log10(timeset[timeset[,4]!=0,4])
									timeset[timeset[,5]!=0,5]<-log10(timeset[timeset[,5]!=0,5])
									median_int<-c(timeset[timeset[,4]!=0,4],timeset[timeset[,5]!=0,5])
									median_int<-median(median_int);
									for(k in 1:length(timeset[,1])){
										if(timeset[k,5]!=0){
											lis_delint_IS[[k]]<-c(lis_delint_IS[[k]],(timeset[k,5]-median_int))
											lis_median_IS[[k]]<-c(lis_median_IS[[k]],median_int)									
											next;
										}
										if(timeset[k,4]!=0){
											lis_delint_IS[[k]]<-c(lis_delint_IS[[k]],(timeset[k,4]-median_int))
											lis_median_IS[[k]]<-c(lis_median_IS[[k]],median_int)	
										}
									}
									stillin[hits[b]]<-FALSE									
								}	
							}
						}
					}	
				}
			}
		}
		# screen other profiles ###############################################################
		if( (logfile$parameters[72]=="TRUE" || logfile$parameters[75]=="TRUE") ){
			cat("- on blank / non-blank profiles");
			peaks<-peaks[stillin,];
			peaks<-peaks[sample(1:length(peaks[,1]),length(peaks[,1]), replace = FALSE),];
			if(logfile$parameters[72]=="TRUE"){ # use blank
				count_b<-0
				if( logfile$parameters[73]=="TRUE" ){ # use subsampling
					max_count_b<-as.numeric(logfile$parameters[74])
				}else{
					max_count_b<-length(peaks[,1])
				}
			}else{
				count_b<-0
				max_count_b<-0
			}
			if(logfile$parameters[75]=="TRUE"){ # use nonblank
				count_nb<-0
				if( logfile$parameters[76]=="TRUE" ){ # use subsampling
					max_count_nb<-as.numeric(logfile$parameters[77])
				}else{
					max_count_nb<-length(peaks[,1])
				}
			}else{
				count_nb<-0
				max_count_nb<-0
			}
			lis_delint_nb<-list(0)
			lis_median_nb<-list(0)
			lis_delint_b<-list(0)
			lis_median_b<-list(0)			
			stillin<-rep(TRUE,length(peaks[,1]))
			for(p in 1:length(timeset[,1])){
				lis_delint_nb[[p]]<-numeric(0)
				lis_median_nb[[p]]<-numeric(0)
				lis_delint_b[[p]]<-numeric(0)
				lis_median_b[[p]]<-numeric(0)
			}
			for(i in 1:length(peaks[,1])){		
				profID<-as.numeric(peaks[i,4])
				timeset[,4:5]<-0;
				timeset[,c(4,5)] <-.Call("fill_timeset",
					as.numeric(timeset),
					as.numeric(profileList_pos[[2]][(profileList_pos[[7]][profileList_pos[[7]][,4]==profID,1]:profileList_pos[[7]][profileList_pos[[7]][,4]==profID,2]),6]), # sampleIDs
					as.numeric(profileList_pos[[2]][(profileList_pos[[7]][profileList_pos[[7]][,4]==profID,1]:profileList_pos[[7]][profileList_pos[[7]][,4]==profID,2]),2]), # intensities
					as.integer(length(timeset[,1])),
					PACKAGE="enviMass"
				)
				if(	 # for blind #############################################3
					(logfile$parameters[72]=="TRUE") &&
					(any(timeset[,5]>0)) &&
					(count_b<max_count_b)
				){ 
					timeset[timeset[,4]!=0,4]<-log10(timeset[timeset[,4]!=0,4])
					timeset[timeset[,5]!=0,5]<-log10(timeset[timeset[,5]!=0,5])
					median_int<-c(timeset[timeset[,4]!=0,4],timeset[timeset[,5]!=0,5])
					median_int<-median(median_int);
					for(k in 1:length(timeset[,1])){
						if(timeset[k,5]!=0){
							lis_delint_b[[k]]<-c(lis_delint_b[[k]],(timeset[k,5]-median_int))
							lis_median_b[[k]]<-c(lis_median_b[[k]],median_int)									
							next;
						}
						if(timeset[k,4]!=0){
							lis_delint_b[[k]]<-c(lis_delint_b[[k]],(timeset[k,4]-median_int))
							lis_median_b[[k]]<-c(lis_median_b[[k]],median_int)	
						}
					}
					count_b<-(count_b+1);
					next;
				}
				if(	 # for non-blind ##########################################
					(logfile$parameters[75]=="TRUE") &&
					(any(timeset[,4]>0)) &&
					(count_nb<max_count_nb)
				){ 
					timeset[timeset[,4]!=0,4]<-log10(timeset[timeset[,4]!=0,4])
					median_int<-median(timeset[timeset[,4]!=0,4]);
					for(k in 1:length(timeset[,1])){
						if(timeset[k,4]!=0){
							lis_delint_nb[[k]]<-c(lis_delint_nb[[k]],(timeset[k,4]-median_int))
							lis_median_nb[[k]]<-c(lis_median_nb[[k]],median_int)									
						}
					}
					count_nb<-(count_nb+1);
				}
				if( # enough subsampling?
					(count_b>=max_count_b) &&
					(count_nb>=max_count_nb) 
				){
					break
				}
			}
		}
		# correct intensities - replace in profileList, recalculate mean_int ###############
		cat("- correcting intensities");
		corfac<-c()
		for(k in 1:length(lis_delint_IS)){
			if( length(lis_delint_IS[[k]])>as.numeric(logfile$parameters[71]) ){
				corfac<-c( corfac,10^median(lis_delint_IS[[k]]) )
			}else{			
				corfac<-c(corfac,1)
			}
		}
		corr_intens <- .Call(	"correct_intens",
								as.numeric(corfac),	  # correction factor
								as.integer(sampleID),       
								as.numeric(profileList_pos[[2]][,2]), # intensities
								as.integer(profileList_pos[[2]][,6]),  
								PACKAGE="enviMass"
							)
		profileList_pos[[2]][,2]<-corr_intens
		for(k in 1:length(profileList_pos[[7]][,8])){
			profileList_pos[[7]][k,16]<-mean(profileList_pos[[2]][(profileList_pos[[7]][k,1]:profileList_pos[[7]][k,2]),2])
		}
		save(profileList_pos,file=file.path(as.character(logfile[[1]]),"results","profileList_pos"));	
		# derive plots #####################################################################
		path=file.path(logfile[[1]],"pics","profnorm_pos")
		png(filename = path, bg = "white", width = 2000, height = 480)    	
			ylimit_del=c(0,0)
			count_IS=c(10000,0)
			count_b=c(10000,0)
			count_nb=c(10000,0)
			for( k in 1:length(lis_delint_IS) ){
				if(length(lis_delint_IS[[k]])>0){ # on IS
					if( min(lis_delint_IS[[k]])<ylimit_del[1] ){
						ylimit_del[1]<-min(lis_delint_IS[[k]])
					}
					if( max(lis_delint_IS[[k]])>ylimit_del[2] ){
						ylimit_del[2]<-max(lis_delint_IS[[k]])
					}
				}
				if(length(lis_delint_IS[[k]])<count_IS[1]){
					count_IS[1]<-length(lis_delint_IS[[k]])
				}
				if(length(lis_delint_IS[[k]])>count_IS[2]){
					count_IS[2]<-length(lis_delint_IS[[k]])
				}
				if( logfile$parameters[72]=="TRUE" ){ # on blank
					if(length(lis_delint_b[[k]])>0){ 
						if( median(lis_delint_b[[k]])<ylimit_del[1] ){
							ylimit_del[1]<-median(lis_delint_b[[k]])
						}
						if( median(lis_delint_b[[k]])>ylimit_del[2] ){
							ylimit_del[2]<-median(lis_delint_b[[k]])
						}
					}
					if(length(lis_delint_b[[k]])<count_b[1]){
						count_b[1]<-length(lis_delint_b[[k]])
					}
					if(length(lis_delint_b[[k]])>count_b[2]){
						count_b[2]<-length(lis_delint_b[[k]])
					}
				}
				if( logfile$parameters[75]=="TRUE" ){ # on non-blank
					if(length(lis_delint_nb[[k]])>0){ 
						if( median(lis_delint_nb[[k]])<ylimit_del[1] ){
							ylimit_del[1]<-median(lis_delint_nb[[k]])
						}
						if( median(lis_delint_nb[[k]])>ylimit_del[2] ){
							ylimit_del[2]<-median(lis_delint_nb[[k]])
						}
					}
					if(length(lis_delint_nb[[k]])<count_nb[1]){
						count_nb[1]<-length(lis_delint_nb[[k]])
					}
					if(length(lis_delint_nb[[k]])>count_nb[2]){
						count_nb[2]<-length(lis_delint_nb[[k]])
					}
				}					
			}				
			par(mar=c(5,4.5,1,10))
			plot.new()
			plot.window(xlim=c(-1,length(timeset[,1])+1),ylim=c(ylimit_del[1]-0.3,ylimit_del[2]))
			box();
			axis(side=1,at=seq(1,length(lis_delint_IS),1),labels=sampleID,las=2,cex.axis=1)
			axis(side=2,cex.axis=1.3);
			title(xlab="Temporal sequence of file IDs",ylab="Deviation from median log10 intensity",cex.lab=1.5)
			abline(h=0,col="red")	
			for(k in 1:length(lis_delint_IS)){
				if(timeset[k,3]!=0){
					abline(v=k,col="orange")
				}
				points( rep(k,length(lis_delint_IS[[k]])),lis_delint_IS[[k]],pch=19,cex=0.7,col="lightgrey" )
			}
			if( (logfile$parameters[72]=="TRUE") ){	
				for(k in 1:length(lis_delint_b)){
					if(length(lis_delint_b[[k]])>0){
						points(k,median(lis_delint_b[[k]]),pch=21,cex=1.7,bg="blue")
					}	
				}
			}
			if( (logfile$parameters[75]=="TRUE") ){	
				for(k in 1:length(lis_delint_nb)){
					if(length(lis_delint_nb[[k]])>0){ 
						points(k,median(lis_delint_nb[[k]]),pch=21,cex=1.7,bg="green3")
					}
				}
			}
			for(k in 1:length(lis_delint_IS)){
				points( k, median(lis_delint_IS[[k]]) , pch=21,cex=1.7, bg="red")
			}			
			plot.window(xlim=c(0,10),ylim=c(0,10))
			legend(-0.2,3,
				pch=c(21,19,21,21,19),
				pt.cex=1.7,cex=1.3,bty="n",
				legend=c(	"IS median deviation",
							"IS single profile deviation",
							"Blank median deviation",
							"Non-blank median deviation",
							"Blank file"),
				pt.bg=c("red","lightgrey","blue","green3","white"),
				col=c("black","lightgrey","black","black","white")
			)
			lines(x=c(-0.2,-0.1),y=c(.45,.45),col="orange")
		dev.off() 				
		exprprofnorm_pos<-list(src=file.path(logfile[[1]],"pics","profnorm_pos"))
		if(isolate(input$Ion_mode)=="positive"){
			output$profnorm<-renderImage(exprprofnorm_pos, deleteFile = FALSE)
		}
		# on peak counts 
		path=file.path(logfile[[1]],"pics","profcount_pos")
		png(filename = path, bg = "white", width = 2000, height = 480)    	
			par(mar=c(5,4.5,1,10))
			plot.new()
			plot.window(xlim=c(-1,length(timeset[,1])+1),ylim=c(count_IS[1]-1,count_IS[2]+1))	
			axis(side=1,at=seq(1,length(lis_delint_IS),1),labels=sampleID,las=2,cex.axis=1)
			axis(2,col="blue",col.ticks="red",col.axis="red",cex.axis=1.3);
			box()
			title(xlab="Temporal sequence of file IDs",ylab="",cex.lab=1.5)
			for(k in 1:length(lis_delint_IS)){
				if(timeset[k,3]!=0){
					abline(v=k,col="orange")
				}
			}
			countit<-c()	
			for(k in 1:length(lis_delint_IS)){
				countit<-c(countit,length(lis_delint_IS[[k]]))	
			}
			lines(countit,col="red",lwd=2)
			mtext("Number of IS peaks", side = 2, line =2.5, col="red",cex=1.5)			
			if( logfile$parameters[72]=="TRUE" ){	
				plot.window( xlim=c(-1,length(timeset[,1])+1),ylim=c(count_b[1]-1,count_b[2]+1) )	
				countit<-c()	
				for( k in 1:length(lis_delint_IS) ){
					countit<-c(countit,length(lis_delint_b[[k]]))	
				}
				lines(countit,col="blue",lwd=2)
				axis(4,col="blue",col.ticks="blue",col.axis="blue",cex.axis=1.3)
				mtext("Number of blank peaks", side = 4, line =2.3, col="blue", cex=1.5)
			}
			if( logfile$parameters[75]=="TRUE" ){	
				plot.window( xlim=c(-1,length(timeset[,1])+1),ylim=c(count_nb[1]-1,count_nb[2]+1) )	
				countit<-c()	
				for( k in 1:length(lis_delint_IS) ){
					countit<-c(countit,length(lis_delint_nb[[k]]))	
				}
				lines(countit,col="green3",lwd=2)
				axis(4,col="green3",col.ticks="green3",col.axis="green3",line=4.5,cex.axis=1.3)
				mtext("Number of non-blank peaks", side = 4, line =6.8,col="green3", cex=1.5)	
			}
		dev.off() 
		exprprofcount_pos<-list(src=file.path(logfile[[1]],"pics","profcount_pos"))
		if(isolate(input$Ion_mode)=="positive"){
			output$profcount<-renderImage(exprprofcount_pos, deleteFile = FALSE)
		}
		####################################################################################
	}
	########################################################################################
	if( any(measurements[,4]=="negative") & file.exists(file.path(as.character(logfile[[1]]),"results","pattern_neg_IS")) ){
		if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="peaklist")){rm(peaklist)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_neg")){rm(profileList_neg,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="profileList_neg")){rm(profileList_neg)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="pattern_neg_IS")){rm(pattern_neg_IS,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="pattern_neg_IS")){rm(pattern_neg_IS)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="patternRT_neg_IS")){rm(patternRT_neg_IS,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="patternRT_neg_IS")){rm(patternRT_neg_IS)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="patternDelRT_neg_IS")){rm(patternDelRT_neg_IS,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="patternDelRT_neg_IS")){rm(patternDelRT_neg_IS)}
		load(file=file.path(as.character(logfile[[1]]),"results","profileList_neg"),envir=as.environment(".GlobalEnv"));	
		load(file=file.path(logfile[[1]],"results","pattern_neg_IS"),envir=as.environment(".GlobalEnv"));
		load(file=file.path(logfile[[1]],"results","patternRT_neg_IS"),envir=as.environment(".GlobalEnv"));
		load(file=file.path(logfile[[1]],"results","patternDelRT_neg_IS"),envir=as.environment(".GlobalEnv"));
		peaks<-profileList_neg[[7]];
		peaklist<-peaks[,c(14,16,15)];
		treshscore<-as.numeric(logfile$parameters[78])
		we1=(1-as.numeric(logfile$parameters[49]))
		we2=(1-we1)
		cat("- screening:");
		screened<-screening(	peaklist, 
								blanklist=FALSE, 
								pattern_neg_IS, 
								RT = patternRT_neg_IS,
								dmz = as.numeric(logfile$parameters[45]), 
								ppm = logfile$parameters[46], 
								dRT = patternDelRT_neg_IS, 
								dRTwithin = as.numeric(logfile$parameters[43]), 
								dRTblank = FALSE, 
								dInt = as.numeric(logfile$parameters[47]), 
								Intcut = as.numeric(logfile$parameters[48]), 
								w1=we1, 
								w2=we2, 	
								w3=0	
		)
		# set matrix to sort & store data from a profile ###################################
		cat("- on IS profiles")
		atPOSIX<-profileList_neg[[3]];
		sampletype<-profileList_neg[[9]];
		sampleID<-profileList_neg[[4]];
		atdate<-c();
		attime<-c();
		for(i in 1:length(atPOSIX)){
				  atdate<-c(atdate, strsplit(atPOSIX[i]," ")[[1]][1]);
				  attime<-c(attime, strsplit(atPOSIX[i]," ")[[1]][2]);
		}
		attime<-as.difftime(attime);
		atdate<-as.Date(atdate);
		ord<-order(as.numeric(atdate),as.numeric(attime),sampleID);
		atPOSIXsort<-atPOSIX[ord];
		atdate<-atdate[ord];
		attime<-attime[ord];
		sampleID<-sampleID[ord];
		sampletype<-sampletype[ord];
		# to retrieve data
		timeset<-matrix(nrow=length(atPOSIX),ncol=5,0);
		for(i in 1:length(sampleID)){
			  if(sampletype[i]=="sample"){
				timeset[i,2]<-as.numeric(sampleID[i]);
			  }
			  if(sampletype[i]=="blank"){
				timeset[i,3]<-as.numeric(sampleID[i]);
			  }
		}
		# screen IS intensity profiles #####################################################
		min_count<-floor(length(profileList_neg[[4]])*as.numeric(logfile$parameters[[70]])/100);
		lis_delint_IS<-list(0)
		lis_median_IS<-list(0)
		stillin<-rep(TRUE,length(peaks[,1]))
		for(p in 1:length(timeset[,1])){
			lis_delint_IS[[p]]<-numeric(0)
			lis_median_IS[[p]]<-numeric(0)
		}
		for(i in 1:length(screened)){
			if(length(screened[[i]])!=1){
				for(j in 1:length(screened[[i]][,1])){
					if(  !(as.character(screened[[i]][j,9]))=="NaN"  ){
						if( as.numeric(as.character(screened[[i]][j,9]))>=treshscore ){
							hits<-strsplit(as.character(screened[[i]][j,2]),"/")[[1]]
							hits<-hits[hits!="-"];
							hits<-as.numeric(hits);
							for(b in 1:length(hits)){
								if( peaks[hits[b],3]>min_count ){
									profID<-as.numeric(peaks[hits[b],4])
									timeset[,4:5]<-0;
									timeset[,c(4,5)] <-.Call("fill_timeset",
														as.numeric(timeset),
														as.numeric(profileList_neg[[2]][(profileList_neg[[7]][profileList_neg[[7]][,4]==profID,1]:profileList_neg[[7]][profileList_neg[[7]][,4]==profID,2]),6]), # sampleIDs
														as.numeric(profileList_neg[[2]][(profileList_neg[[7]][profileList_neg[[7]][,4]==profID,1]:profileList_neg[[7]][profileList_neg[[7]][,4]==profID,2]),2]), # intensities
														as.integer(length(timeset[,1])),
														PACKAGE="enviMass"
									)
									timeset[timeset[,4]!=0,4]<-log10(timeset[timeset[,4]!=0,4])
									timeset[timeset[,5]!=0,5]<-log10(timeset[timeset[,5]!=0,5])
									median_int<-c(timeset[timeset[,4]!=0,4],timeset[timeset[,5]!=0,5])
									median_int<-median(median_int);
									for(k in 1:length(timeset[,1])){
										if(timeset[k,5]!=0){
											lis_delint_IS[[k]]<-c(lis_delint_IS[[k]],(timeset[k,5]-median_int))
											lis_median_IS[[k]]<-c(lis_median_IS[[k]],median_int)									
											next;
										}
										if(timeset[k,4]!=0){
											lis_delint_IS[[k]]<-c(lis_delint_IS[[k]],(timeset[k,4]-median_int))
											lis_median_IS[[k]]<-c(lis_median_IS[[k]],median_int)	
										}
									}
									stillin[hits[b]]<-FALSE									
								}	
							}
						}
					}	
				}
			}
		}
		# screen other profiles ###############################################################
		if( logfile$parameters[72]=="TRUE" || logfile$parameters[75]=="TRUE" ){
			cat("- on blank / non-blank profiles");
			peaks<-peaks[stillin,];
			peaks<-peaks[sample(1:length(peaks[,1]),length(peaks[,1]), replace = FALSE),];
			if(logfile$parameters[72]=="TRUE"){ # use blank
				count_b<-0
				if( logfile$parameters[73]=="TRUE" ){ # use subsampling
					max_count_b<-as.numeric(logfile$parameters[74])
				}else{
					max_count_b<-length(peaks[,1])
				}
			}else{
				count_b<-0
				max_count_b<-0
			}
			if(logfile$parameters[75]=="TRUE"){ # use nonblank
				count_nb<-0
				if( logfile$parameters[76]=="TRUE" ){ # use subsampling
					max_count_nb<-as.numeric(logfile$parameters[77])
				}else{
					max_count_nb<-length(peaks[,1])
				}
			}else{
				count_nb<-0
				max_count_nb<-0
			}
			lis_delint_nb<-list(0)
			lis_median_nb<-list(0)
			lis_delint_b<-list(0)
			lis_median_b<-list(0)			
			stillin<-rep(TRUE,length(peaks[,1]))
			for(p in 1:length(timeset[,1])){
				lis_delint_nb[[p]]<-numeric(0)
				lis_median_nb[[p]]<-numeric(0)
				lis_delint_b[[p]]<-numeric(0)
				lis_median_b[[p]]<-numeric(0)
			}
			for(i in 1:length(peaks[,1])){		
				profID<-as.numeric(peaks[i,4])
				timeset[,4:5]<-0;
				timeset[,c(4,5)] <-.Call("fill_timeset",
					as.numeric(timeset),
					as.numeric(profileList_neg[[2]][(profileList_neg[[7]][profileList_neg[[7]][,4]==profID,1]:profileList_neg[[7]][profileList_neg[[7]][,4]==profID,2]),6]), # sampleIDs
					as.numeric(profileList_neg[[2]][(profileList_neg[[7]][profileList_neg[[7]][,4]==profID,1]:profileList_neg[[7]][profileList_neg[[7]][,4]==profID,2]),2]), # intensities
					as.integer(length(timeset[,1])),
					PACKAGE="enviMass"
				)
				if(	 # for blind #############################################3
					(logfile$parameters[72]=="TRUE") &&
					(any(timeset[,5]>0)) &&
					(count_b<max_count_b)
				){ 
					timeset[timeset[,4]!=0,4]<-log10(timeset[timeset[,4]!=0,4])
					timeset[timeset[,5]!=0,5]<-log10(timeset[timeset[,5]!=0,5])
					median_int<-c(timeset[timeset[,4]!=0,4],timeset[timeset[,5]!=0,5])
					median_int<-median(median_int);
					for(k in 1:length(timeset[,1])){
						if(timeset[k,5]!=0){
							lis_delint_b[[k]]<-c(lis_delint_b[[k]],(timeset[k,5]-median_int))
							lis_median_b[[k]]<-c(lis_median_b[[k]],median_int)									
							next;
						}
						if(timeset[k,4]!=0){
							lis_delint_b[[k]]<-c(lis_delint_b[[k]],(timeset[k,4]-median_int))
							lis_median_b[[k]]<-c(lis_median_b[[k]],median_int)	
						}
					}
					count_b<-(count_b+1);
					next;
				}
				if(	 # for non-blind ##########################################
					(logfile$parameters[75]=="TRUE") &&
					(any(timeset[,4]>0)) &&
					(count_nb<max_count_nb)
				){ 
					timeset[timeset[,4]!=0,4]<-log10(timeset[timeset[,4]!=0,4])
					median_int<-median(timeset[timeset[,4]!=0,4]);
					for(k in 1:length(timeset[,1])){
						if(timeset[k,4]!=0){
							lis_delint_nb[[k]]<-c(lis_delint_nb[[k]],(timeset[k,4]-median_int))
							lis_median_nb[[k]]<-c(lis_median_nb[[k]],median_int)									
						}
					}
					count_nb<-(count_nb+1);
				}
				if( # enough subsampling?
					(count_b>=max_count_b) &&
					(count_nb>=max_count_nb) 
				){
					break
				}
			}
		}
		# correct intensities - replace in profileList, recalculate mean_int ###############
		cat("- correcting intensities");
		corfac<-c()
		for(k in 1:length(lis_delint_IS)){
			if( length(lis_delint_IS[[k]])>as.numeric(logfile$parameters[71]) ){
				corfac<-c( corfac,10^median(lis_delint_IS[[k]]) )
			}else{			
				corfac<-c(corfac,1)
			}
		}
		corr_intens <- .Call(	"correct_intens",
								as.numeric(corfac),	  # correction factor
								as.integer(sampleID),       
								as.numeric(profileList_neg[[2]][,2]), # intensities
								as.integer(profileList_neg[[2]][,6]),  
								PACKAGE="enviMass"
							)
		profileList_neg[[2]][,2]<-corr_intens
		for(k in 1:length(profileList_neg[[7]][,8])){
			profileList_neg[[7]][k,16]<-mean(profileList_neg[[2]][(profileList_neg[[7]][k,1]:profileList_neg[[7]][k,2]),2])
		}
		save(profileList_neg,file=file.path(as.character(logfile[[1]]),"results","profileList_neg"));	
		# derive plots #####################################################################
		path=file.path(logfile[[1]],"pics","profnorm_neg")
		png(filename = path, bg = "white", width = 2000, height = 480)    	
			ylimit_del=c(0,0)
			count_IS=c(10000,0)
			count_b=c(10000,0)
			count_nb=c(10000,0)
			for( k in 1:length(lis_delint_IS) ){
				if(length(lis_delint_IS[[k]])>0){ # on IS
					if( min(lis_delint_IS[[k]])<ylimit_del[1] ){
						ylimit_del[1]<-min(lis_delint_IS[[k]])
					}
					if( max(lis_delint_IS[[k]])>ylimit_del[2] ){
						ylimit_del[2]<-max(lis_delint_IS[[k]])
					}
				}
				if(length(lis_delint_IS[[k]])<count_IS[1]){
					count_IS[1]<-length(lis_delint_IS[[k]])
				}
				if(length(lis_delint_IS[[k]])>count_IS[2]){
					count_IS[2]<-length(lis_delint_IS[[k]])
				}
				
				if( (logfile$parameters[72]=="TRUE") ){ # on blank
					if(length(lis_delint_b[[k]])>0){ 
						if( median(lis_delint_b[[k]])<ylimit_del[1] ){
							ylimit_del[1]<-median(lis_delint_b[[k]])
						}
						if( median(lis_delint_b[[k]])>ylimit_del[2] ){
							ylimit_del[2]<-median(lis_delint_b[[k]])
						}
					}
					if(length(lis_delint_b[[k]])<count_b[1]){
						count_b[1]<-length(lis_delint_b[[k]])
					}
					if(length(lis_delint_b[[k]])>count_b[2]){
						count_b[2]<-length(lis_delint_b[[k]])
					}
				}
				
				if( (logfile$parameters[75]=="TRUE") ){ # on non-blank
					if(length(lis_delint_nb[[k]])>0){ 
						if( median(lis_delint_nb[[k]])<ylimit_del[1] ){
							ylimit_del[1]<-median(lis_delint_nb[[k]])
						}
						if( median(lis_delint_nb[[k]])>ylimit_del[2] ){
							ylimit_del[2]<-median(lis_delint_nb[[k]])
						}
					}
					if(length(lis_delint_nb[[k]])<count_nb[1]){
						count_nb[1]<-length(lis_delint_nb[[k]])
					}
					if(length(lis_delint_nb[[k]])>count_nb[2]){
						count_nb[2]<-length(lis_delint_nb[[k]])
					}
				}					
			}				
			par(mar=c(5,4.5,1,10))
			plot.new()
			plot.window(xlim=c(-1,length(timeset[,1])+1),ylim=c(ylimit_del[1]-0.3,ylimit_del[2]))
			box();
			axis(side=1,at=seq(1,length(lis_delint_IS),1),labels=sampleID,las=2,cex.axis=1)
			axis(side=2,cex.axis=1.3);
			title(xlab="Temporal sequence of file IDs",ylab="Deviation from median log10 intensity",cex.lab=1.5)
			abline(h=0,col="red")	
			for(k in 1:length(lis_delint_IS)){
				if(timeset[k,3]!=0){
					abline(v=k,col="orange")
				}
				points( rep(k,length(lis_delint_IS[[k]])),lis_delint_IS[[k]],pch=19,cex=0.7,col="lightgrey" )
			}
			if( (logfile$parameters[72]=="TRUE") ){	
				for(k in 1:length(lis_delint_b)){
					if(length(lis_delint_b[[k]])>0){
						points(k,median(lis_delint_b[[k]]),pch=21,cex=1.7,bg="blue")
					}	
				}
			}
			if( (logfile$parameters[75]=="TRUE") ){	
				for(k in 1:length(lis_delint_nb)){
					if(length(lis_delint_nb[[k]])>0){ 
						points(k,median(lis_delint_nb[[k]]),pch=21,cex=1.7,bg="green3")
					}
				}
			}
			for(k in 1:length(lis_delint_IS)){
				points( k, median(lis_delint_IS[[k]]) , pch=21,cex=1.7, bg="red")
			}			
			plot.window(xlim=c(0,10),ylim=c(0,10))
			legend(-0.2,3,
				pch=c(21,19,21,21,19),
				pt.cex=1.7,cex=1.3,bty="n",
				legend=c(	"IS median deviation",
							"IS single profile deviation",
							"Blank median deviation",
							"Non-blank median deviation",
							"Blank file"),
				pt.bg=c("red","lightgrey","blue","green3","white"),
				col=c("black","lightgrey","black","black","white")
			)
			lines(x=c(-0.2,-0.1),y=c(.45,.45),col="orange")
		dev.off() 				
		exprprofnorm_neg<-list(src=file.path(logfile[[1]],"pics","profnorm_neg"))
		if(isolate(input$Ion_mode)=="negative"){
			output$profnorm<-renderImage(exprprofnorm_neg, deleteFile = FALSE)
		}
		# on peak counts 
		path=file.path(logfile[[1]],"pics","profcount_neg")
		png(filename = path, bg = "white", width = 2000, height = 480)    	
			par(mar=c(5,4.5,1,10))
			plot.new()
			plot.window(xlim=c(-1,length(timeset[,1])+1),ylim=c(count_IS[1]-1,count_IS[2]+1))	
			axis(side=1,at=seq(1,length(lis_delint_IS),1),labels=sampleID,las=2,cex.axis=1)
			axis(2,col="blue",col.ticks="red",col.axis="red",cex.axis=1.3);
			box()
			title(xlab="Temporal sequence of file IDs",ylab="",cex.lab=1.5)
			for(k in 1:length(lis_delint_IS)){
				if(timeset[k,3]!=0){
					abline(v=k,col="orange")
				}
			}
			countit<-c()	
			for(k in 1:length(lis_delint_IS)){
				countit<-c(countit,length(lis_delint_IS[[k]]))	
			}
			lines(countit,col="red",lwd=2)
			mtext("Number of IS peaks", side = 2, line =2.5, col="red",cex=1.5)			
			if( logfile$parameters[72]=="TRUE" ){	
				plot.window( xlim=c(-1,length(timeset[,1])+1),ylim=c(count_b[1]-1,count_b[2]+1) )	
				countit<-c()	
				for(k in 1:length(lis_delint_IS)){
					countit<-c(countit,length(lis_delint_b[[k]]))	
				}
				lines(countit,col="blue",lwd=2)
				axis(4,col="blue",col.ticks="blue",col.axis="blue",cex.axis=1.3)
				mtext("Number of blank peaks", side = 4, line =2.3, col="blue", cex=1.5)
			}
			if( logfile$parameters[75]=="TRUE" ){	
				plot.window( xlim=c(-1,length(timeset[,1])+1),ylim=c(count_nb[1]-1,count_nb[2]+1) )	
				countit<-c()	
				for(k in 1:length(lis_delint_IS)){
					countit<-c(countit,length(lis_delint_nb[[k]]))	
				}
				lines(countit,col="green3",lwd=2)
				axis(4,col="green3",col.ticks="green3",col.axis="green3",line=4.5,cex.axis=1.3)
				mtext("Number of non-blank peaks", side = 4, line =6.8,col="green3", cex=1.5)	
			}
		dev.off() 
		exprprofcount_neg<-list(src=file.path(logfile[[1]],"pics","profcount_neg"))
		if(isolate(input$Ion_mode)=="negative"){
			output$profcount<-renderImage(exprprofcount_neg, deleteFile = FALSE)
		}
		####################################################################################

	}
	########################################################################################
	logfile$summary[9,2]<<-"TRUE";
    logfile$summary[9,2]<-"TRUE";
	logfile$Tasks_to_redo[15]<-"FALSE";
	logfile$Tasks_to_redo[15]<<-"FALSE";
    save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
	summa[9,2]<<-"done"
    summa[9,2]<-"done"
	output$summar<<-renderTable(summa)
    cat("IS-based profile normalization done \n");
    output$dowhat<<-renderText("IS-normalization done ... wait")

}else{

	if(logfile$workflow[15]=="no"){
		logfile$summary[9,2]<<-"FALSE";
		logfile$summary[9,2]<-"FALSE";
		path=file.path(logfile[[1]],"pics","profnorm_pos")
			png(filename = path, bg = "white")
			plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
			dev.off()
		    expr30p<-list(src=path)
			output$profnorm_pos<-renderImage(expr30p, deleteFile = FALSE)		
		path=file.path(logfile[[1]],"pics","profcount_pos")
			png(filename = path, bg = "white")
			plot.new();plot.window(xlim=c(1,1),ylim=c(1,1));#box();text(1,1,label="not available",cex=1.5,col="darkred")
			dev.off()
		    expr31p<-list(src=path)
			output$profcount_pos<-renderImage(expr31p, deleteFile = FALSE)
		path=file.path(logfile[[1]],"pics","profnorm_neg")
			png(filename = path, bg = "white")
			plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
			dev.off()
		    expr30n<-list(src=path)
			output$profnorm_neg<-renderImage(expr30n, deleteFile = FALSE)		
		path=file.path(logfile[[1]],"pics","profcount_neg")
			png(filename = path, bg = "white")
			plot.new();plot.window(xlim=c(1,1),ylim=c(1,1));#box();text(1,1,label="not available",cex=1.5,col="darkred")
			dev.off()
		    expr31n<-list(src=path)
			output$profcount_neg<-renderImage(expr31n, deleteFile = FALSE)		
	}
	logfile$Tasks_to_redo[15]<-"FALSE";
	logfile$Tasks_to_redo[15]<<-"FALSE";
    save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
    summa[9,2]<-"skipped";
    summa[9,2]<<-"skipped";
	output$summar<<-renderTable(summa)
    cat("IS-based profile normalization skipped \n");
    output$dowhat<<-renderText("IS-normalization skipped ... wait")

}