
    ############################################################################
    # retrieve monoisotopic masses for IS ######################################
    ############################################################################
    measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
    measurements<-measurements[measurements[,8]=="TRUE",]
    leng<-length(measurements[,1])
	mz_pos<-c();
	RT_pos<-c();
	mz_neg<-c();
	RT_neg<-c();
    if(any(measurements[,4]=="positive")){ # positive ##########################
      if(logfile$parameters[30]=="Internal standards"){
        if(file.exists(file.path(logfile[[1]],"results","intmass_pos_IS"))){
			if(any(objects(envir=as.environment(".GlobalEnv"))=="intmass_pos_IS")){rm(intmass_pos_IS,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="intmass_pos_IS")){rm(intmass_pos_IS)}
			load(file=file.path(logfile[[1]],"results","intmass_pos_IS"),envir=as.environment(".GlobalEnv"))                    
			mz_pos<-c(mz_pos,intmass_pos_IS[,1]);
			RT_pos<-c(RT_pos,intmass_pos_IS[,2]);
		}else{cat("\n IS recalibration masses not found ... recalibration skipped!")}
	  }
      if(logfile$parameters[30]=="Target compounds"){
        if(file.exists(file.path(logfile[[1]],"results","intmass_pos_target"))){	  
			if(any(objects(envir=as.environment(".GlobalEnv"))=="intmass_pos_target")){rm(intmass_pos_target,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="intmass_pos_target")){rm(intmass_pos_target)}
			load(file=file.path(logfile[[1]],"results","intmass_pos_target"),envir=as.environment(".GlobalEnv"))                    
			mz_pos<-c(mz_pos,intmass_pos_target[,1]);
			RT_pos<-c(RT_pos,intmass_pos_target[,2]);
		}else{cat("\n target recalibration masses not found ... recalibration skipped!")}      
	  }
      if(logfile$parameters[30]=="both"){
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
		}else{cat("\n IS target masses not found ... recalibration skipped?")}			
	  }
	  mz_pos<-c(as.numeric(as.character(mz_pos)));
      RT_pos<-c(as.numeric(as.character(RT_pos)));
    }
    if(any(measurements[,4]=="negative")){ # negative ##########################
      if(logfile$parameters[30]=="Internal standards"){
        if(file.exists(file.path(logfile[[1]],"results","intmass_neg_IS"))){	
			if(any(objects(envir=as.environment(".GlobalEnv"))=="intmass_neg_IS")){rm(intmass_neg_IS,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="intmass_neg_IS")){rm(intmass_neg_IS)}				
			load(file=file.path(logfile[[1]],"results","intmass_neg_IS"),envir=as.environment(".GlobalEnv"))                    
			mz_neg<-c(mz_neg,intmass_neg_IS[,1]);
			RT_neg<-c(RT_neg,intmass_neg_IS[,2]);
		}else{cat("\n IS recalibration masses not found ... recalibration skipped!")}		
      }
      if(logfile$parameters[30]=="Target compounds"){
        if(file.exists(file.path(logfile[[1]],"results","intmass_neg_target"))){		  
			if(any(objects(envir=as.environment(".GlobalEnv"))=="intmass_neg_target")){rm(intmass_neg_target,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="intmass_neg_target")){rm(intmass_neg_target)}				
			load(file=file.path(logfile[[1]],"results","intmass_neg_target"),envir=as.environment(".GlobalEnv"))                    
			mz_neg<-c(mz_neg,intmass_neg_target[,1]);
			RT_neg<-c(RT_neg,intmass_neg_target[,2]);
		}else{cat("\n target recalibration masses not found ... recalibration skipped!")}   		
      }
      if(logfile$parameters[30]=="both"){
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
    for( i in 1:length(measurements[,1]) ){  
	if( (measurements[i,8]=="TRUE")&(measurements[i,12]=="FALSE")  ){	  
	  if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
	  if(any(objects()=="peaklist")){rm(peaklist)}
      load(file=file.path(logfile[[1]],"peaklist","",as.character(measurements[i,1])),envir=as.environment(".GlobalEnv"));   
      if( (measurements[i,4]=="positive") & (measurements[i,8]=="TRUE") ){  
		if( length(mz_pos)>0 ){
			  peak_recal<-recalib(
					  peaklist=peaklist[,c(1,4,5)],
					  mz=mz_pos,
					  tolmz=as.numeric(logfile$parameters[31]),
					  ppm=as.character(logfile$parameters[32]),
					  ret=RT_pos,
					  tolret=as.numeric(logfile$parameters[33]),
					  what="mass",
					  one=TRUE,
					  knot=5,
					  plotit=TRUE,
					  path_1=file.path(logfile[[1]],"pics",paste("recal_",as.character(measurements[i,1]),sep="")),
					  path_2=file.path(logfile[[1]],"results","recalibration",paste("recal_gam_",as.character(measurements[i,1]),sep=""))
				)
				if(length(peak_recal)>1){
				  peaklist[,c(12,13,14)]<-peak_recal
				  save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(measurements[i,1])));   
				  measurements[i,12]<-"FALSE";
				  cat(paste("Mass recalibration ",i," of ",leng," done.\n",sep=""))            
				}else{
				  cat(paste(peak_recal," \n",sep=""))
				  cat(paste("Mass recalibration ",i," of ",leng," - infeasible.\n",sep=""))        
				}
				measurements[i,12]<-TRUE;
		}else{		
			png(filename = file.path(logfile[[1]],"pics",paste("recal_",as.character(measurements[i,1]),sep="")), bg = "white", width = 1100, height= 300)
			plot.new()
			plot.window(xlim=c(0,1),ylim=c(0,1))
			text(0.5,0.5,"No recalibration \n compounds listed\n(positive ionization).",cex=1)
			dev.off();
			measurements[i,12]<-FALSE;
		}
	  }
      if( (measurements[i,4]=="negative")  & (measurements[i,8]=="TRUE") ){
		if(length(mz_neg)>0){
			peak_recal<-recalib(
				  peaklist=peaklist[,c(1,4,5)],
				  mz=mz_neg,
				  tolmz=as.numeric(logfile$parameters[31]),
				  ppm=as.character(logfile$parameters[32]),
				  ret=RT_neg,
				  tolret=as.numeric(logfile$parameters[33]),
				  what="mass",
				  one=TRUE,
				  knot=5,
				  plotit=TRUE,
				  path_1=file.path(logfile[[1]],"pics",paste("recal_",as.character(measurements[i,1]),sep="")),
				  path_2=file.path(logfile[[1]],"results","recalibration",paste("recal_gam_",as.character(measurements[i,1]),sep=""))
			)
			if(length(peak_recal)>1){
			  peaklist[,c(12,13,14)]<-peak_recal
			  save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(measurements[i,1])));   
			  measurements[i,12]<-"FALSE";
			  cat(paste("Mass recalibration ",i," of ",leng," done.\n",sep=""))
			}else{
			  cat(paste(peak_recal," \n",sep=""))
			  cat(paste("Mass recalibration ",i," of ",leng," - infeasible.\n",sep=""))
			}
			measurements[i,12]<-TRUE;
		}else{		
			png(filename = file.path(logfile[[1]],"pics",paste("recal_",as.character(measurements[i,1]),sep="")), bg = "white", width = 1100, height= 300)
			plot.new()
			plot.window(xlim=c(0,1),ylim=c(0,1))
			text(0.5,0.5,"No recalibration \n compounds listed\n(negative ionization).",cex=1)
			dev.off();
			measurements[i,12]<-FALSE;
		}	
	  }
	}  
    }
	write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
    measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	output$measurements<-DT::renderDataTable(read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character")); 

 

