
    ############################################################################
    # retrieve monoisotopic masses for IS ######################################
    ############################################################################
    measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
    measurements_incl<-measurements[measurements[,"include"]=="TRUE",,drop=FALSE]
	
    leng<-length(measurements_incl[,1])
	mz_pos<-c();
	RT_pos<-c();
	mz_neg<-c();
	RT_neg<-c();
    if(any(measurements_incl[,"Mode"]=="positive")){ # positive ##########################
      if(logfile$parameters$recal_use=="Internal standards"){
        if(file.exists(file.path(logfile[[1]],"results","intmass_pos_IS"))){
			if(any(objects(envir=as.environment(".GlobalEnv"))=="intmass_pos_IS")){rm(intmass_pos_IS,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="intmass_pos_IS")){rm(intmass_pos_IS)}
			load(file=file.path(logfile[[1]],"results","intmass_pos_IS"),envir=as.environment(".GlobalEnv"))                    
			mz_pos<-c(mz_pos,intmass_pos_IS[,1]);
			RT_pos<-c(RT_pos,intmass_pos_IS[,2]);
		}else{cat("\n IS recalibration masses not found ... recalibration skipped!")}
	  }
      if(logfile$parameters$recal_use=="Target compounds"){
        if(file.exists(file.path(logfile[[1]],"results","intmass_pos_target"))){	  
			if(any(objects(envir=as.environment(".GlobalEnv"))=="intmass_pos_target")){rm(intmass_pos_target,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="intmass_pos_target")){rm(intmass_pos_target)}
			load(file=file.path(logfile[[1]],"results","intmass_pos_target"),envir=as.environment(".GlobalEnv"))                    
			mz_pos<-c(mz_pos,intmass_pos_target[,1]);
			RT_pos<-c(RT_pos,intmass_pos_target[,2]);
		}else{cat("\n Target recalibration masses not found ... recalibration skipped!")}      
	  }
      if(logfile$parameters$recal_use=="both"){
        if(file.exists(file.path(logfile[[1]],"results","intmass_pos_IS"))){	  
			if(any(objects(envir=as.environment(".GlobalEnv"))=="intmass_pos_IS")){rm(intmass_pos_IS,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="intmass_pos_IS")){rm(intmass_pos_IS)}
			load(file=file.path(logfile[[1]],"results","intmass_pos_IS"),envir=as.environment(".GlobalEnv"))                    
			mz_pos<-c(mz_pos,intmass_pos_IS[,1]);
			RT_pos<-c(RT_pos,intmass_pos_IS[,2]);
		}else{cat("\n IS recalibration masses not found ... recalibration skipped?")}		
        if(file.exists(file.path(logfile[[1]],"results","intmass_pos_target"))){	
			if(any(objects(envir=as.environment(".GlobalEnv"))=="intmass_pos_target")){rm(intmass_pos_target,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="intmass_pos_target")){rm(intmass_pos_target)}		
			load(file=file.path(logfile[[1]],"results","intmass_pos_target"),envir=as.environment(".GlobalEnv"))                    
			mz_pos<-c(mz_pos,intmass_pos_target[,1]);
			RT_pos<-c(RT_pos,intmass_pos_target[,2]);
		}else{cat("\n Target recalibration masses not found ... recalibration skipped?")}			
	  }
	  mz_pos<-c(as.numeric(as.character(mz_pos)));
      RT_pos<-c(as.numeric(as.character(RT_pos)));
    }
    if(any(measurements_incl[,"Mode"]=="negative")){ # negative ##########################
      if(logfile$parameters$recal_use=="Internal standards"){
        if(file.exists(file.path(logfile[[1]],"results","intmass_neg_IS"))){	
			if(any(objects(envir=as.environment(".GlobalEnv"))=="intmass_neg_IS")){rm(intmass_neg_IS,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="intmass_neg_IS")){rm(intmass_neg_IS)}				
			load(file=file.path(logfile[[1]],"results","intmass_neg_IS"),envir=as.environment(".GlobalEnv"))                    
			mz_neg<-c(mz_neg,intmass_neg_IS[,1]);
			RT_neg<-c(RT_neg,intmass_neg_IS[,2]);
		}else{cat("\n IS recalibration masses not found ... recalibration skipped!")}		
      }
      if(logfile$parameters$recal_use=="Target compounds"){
        if(file.exists(file.path(logfile[[1]],"results","intmass_neg_target"))){		  
			if(any(objects(envir=as.environment(".GlobalEnv"))=="intmass_neg_target")){rm(intmass_neg_target,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="intmass_neg_target")){rm(intmass_neg_target)}				
			load(file=file.path(logfile[[1]],"results","intmass_neg_target"),envir=as.environment(".GlobalEnv"))                    
			mz_neg<-c(mz_neg,intmass_neg_target[,1]);
			RT_neg<-c(RT_neg,intmass_neg_target[,2]);
		}else{cat("\n target recalibration masses not found ... recalibration skipped!")}   		
      }
      if(logfile$parameters$recal_use=="both"){
        if(file.exists(file.path(logfile[[1]],"results","intmass_neg_IS"))){		  
			if(any(objects(envir=as.environment(".GlobalEnv"))=="intmass_neg_IS")){rm(intmass_neg_IS,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="intmass_neg_IS")){rm(intmass_neg_IS)}				
			load(file=file.path(logfile[[1]],"results","intmass_neg_IS"),envir=as.environment(".GlobalEnv"))                    
			mz_neg<-c(mz_neg,intmass_neg_IS[,1]);
			RT_neg<-c(RT_neg,intmass_neg_IS[,2]);
		}else{cat("\n IS recalibration masses not found ... recalibration skipped?")}	
        if(file.exists(file.path(logfile[[1]],"results","intmass_neg_target"))){		
			if(any(objects(envir=as.environment(".GlobalEnv"))=="intmass_neg_target")){rm(intmass_neg_target,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="intmass_neg_target")){rm(intmass_neg_target)}				
			load(file=file.path(logfile[[1]],"results","intmass_neg_target"),envir=as.environment(".GlobalEnv"))                    
			mz_neg<-c(mz_neg,intmass_neg_target[,1]);
			RT_neg<-c(RT_neg,intmass_neg_target[,2]);
		}else{cat("\n target recalibration masses not found ... recalibration skipped?")}  		
      }
	  mz_neg<-c(as.numeric(as.character(mz_neg)));
      RT_neg<-c(as.numeric(as.character(RT_neg)));
    }
    for( i in 1:length(measurements_incl[,"ID"]) ){  
		if( (measurements_incl[i,"include"]=="TRUE")&(measurements_incl[i,"recal"]=="FALSE")  ){	  
		  if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
		  if(any(objects()=="peaklist")){rm(peaklist)}
		  load(file=file.path(logfile[[1]],"peaklist","",as.character(measurements_incl[i,"ID"])),envir=as.environment(".GlobalEnv"));   
		  if( (measurements_incl[i,"Mode"]=="positive") & (measurements_incl[i,"include"]=="TRUE") ){  
			if( length(mz_pos)>0 ){
				  peak_recal<-recalib(
						peaklist=peaklist[,c("m/z","max_int","RT")],
						mz=mz_pos,
						tolmz=as.numeric(logfile$parameters$recal_dmz),
						ppm=as.character(logfile$parameters$recal_ppm),
						ret=RT_pos,
						tolret=as.numeric(logfile$parameters$recal_drt),
						what="mass",
						one=TRUE,
						knot=5,
						plotit=TRUE,
						path_1=file.path(logfile[[1]],"pics",paste("recal_",as.character(measurements_incl[i,"ID"]),sep="")),
						path_2=file.path(logfile[[1]],"results","recalibration",paste("recal_gam_",as.character(measurements_incl[i,"ID"]),sep="")),
						plot_ppm=c(2,5,10),
						max_recal=as.numeric(logfile$parameters$recal_maxdmz)
					)
					if(length(peak_recal)>1){
					  peaklist[,c(12,13,14)]<-peak_recal
					  save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(measurements_incl[i,"ID"])));   
					  measurements_incl[i,"recal"]<-"FALSE";
					  cat(paste("Mass recalibration ",i," of ",leng," done.\n",sep=""))            
					}else{
					  cat(paste(peak_recal," \n",sep=""))
					  cat(paste("Mass recalibration ",i," of ",leng," - infeasible.\n",sep=""))        
					}
					measurements[
						measurements[,"ID"]==measurements_incl[i,"ID"]
					,"recal"]<-TRUE;
			}else{		
				png(filename = file.path(logfile[[1]],"pics",paste("recal_",as.character(measurements_incl[i,"ID"]),sep="")), bg = "white", width = 1100, height= 300)
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,"No recalibration \n compounds listed\n(positive ionization).",cex=1)
				dev.off();
				measurements[
					measurements[,"ID"]==measurements_incl[i,1]
				,"recal"]<-FALSE;
			}
		  }
		  if( (measurements_incl[i,"Mode"]=="negative")  & (measurements_incl[i,"include"]=="TRUE") ){
			if(length(mz_neg)>0){
				peak_recal<-recalib(
					peaklist=peaklist[,c(1,4,5)],
					mz=mz_neg,
					tolmz=as.numeric(logfile$parameters$recal_dmz),
					ppm=as.character(logfile$parameters$recal_ppm),
					ret=RT_neg,
					tolret=as.numeric(logfile$parameters$recal_drt),
					what="mass",
					one=TRUE,
					knot=5,
					plotit=TRUE,
					path_1=file.path(logfile[[1]],"pics",paste("recal_",as.character(measurements_incl[i,"ID"]),sep="")),
					path_2=file.path(logfile[[1]],"results","recalibration",paste("recal_gam_",as.character(measurements_incl[i,"ID"]),sep="")),
					plot_ppm=c(2,5,10),
					max_recal=as.numeric(logfile$parameters$recal_maxdmz)					  
				)
				if(length(peak_recal)>1){
				  peaklist[,c(12,13,14)]<-peak_recal
				  save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(measurements_incl[i,"ID"])));   
				  measurements_incl[i,"recal"]<-"FALSE";
				  cat(paste("Mass recalibration ",i," of ",leng," done.\n",sep=""))
				}else{
				  cat(paste(peak_recal," \n",sep=""))
				  cat(paste("Mass recalibration ",i," of ",leng," - infeasible.\n",sep=""))
				}
				measurements[
						measurements[,"ID"]==measurements_incl[i,"ID"]
				,"recal"]<-TRUE;
			}else{		
				png(filename = file.path(logfile[[1]],"pics",paste("recal_",as.character(measurements_incl[i,"ID"]),sep="")), bg = "white", width = 1100, height= 300)
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,"No recalibration \n compounds listed\n(negative ionization).",cex=1)
				dev.off();
				measurements[
					measurements[,"ID"]==measurements_incl[i,"ID"]
				,"recal"]<-FALSE;
			}	
		  }
		}  
    }
	write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);

 

