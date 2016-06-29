# Peak picking ##################################################################

	output$dowhat<-renderText("Peak picking ... please wait");
	if(any(search()=="package:nlme")){detach(package:nlme,force=TRUE);addit<-TRUE}else{addit<-FALSE}
    measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
    leng<-length(measurements[,1]);         
    for(i in 1:leng){ 
        # (measurement included & not yet picked) OR (peakpick forced) 
            if( (measurements[i,8]=="TRUE")&&(measurements[i,10]=="FALSE") ){
				cat(paste("\n    Peak picking sample ",as.character(i)," of ",as.character(leng),": "));                        
				MSlist<-enviPick::readMSdata(
					file.path(logfile[[1]],"files",paste(as.character(measurements[i,1]),".mzXML",sep="")),
					MSlevel=logfile$parameters[1],  # MSlevel
					progbar=logfile$parameters[21], # progbar
					minRT=FALSE,maxRT=FALSE,minmz=FALSE,maxmz=FALSE,
					ion_mode=measurements[i,4]
				);
				cat(" data read -");  				
				if(logfile$parameters[92]>"0"){
					len1<-length(MSlist[[4]][[2]][,2])
					MSlist[[4]][[2]]<-MSlist[[4]][[2]][
						MSlist[[4]][[2]][,2]>=quantile(MSlist[[4]][[2]][,2],(as.numeric(logfile$parameters[92])/100))
					,]
					len2<-length(MSlist[[4]][[2]][,2])
					cat(paste(" ",as.character(len1-len2),"of",as.character(len1),"data points discarded -"))
				}
				MSlist<-enviPick::mzagglom(
					MSlist,
					((as.numeric(logfile$parameters[[3]])*2)+1),
					ppm=TRUE,  
					as.numeric(logfile$parameters[[2]]),
					as.numeric(logfile$parameters[[4]]),
					10^as.numeric(logfile$parameters[[14]]),
					progbar=as.logical(logfile$parameters[21])
				);
				cat(" partitioned -");
				MSlist<-enviPick::mzclust(      
					MSlist,
					as.numeric(logfile$parameters[[3]]),
					ppm=TRUE,
					60,
					as.numeric(logfile$parameters[[4]]),
					10^as.numeric(logfile$parameters[[14]]),
					progbar=as.logical(logfile$parameters[21]),
					merged=TRUE,from=FALSE,to=FALSE
				);
				cat(" clustered -");
				MSlist<-enviPick::mzpick(       
					MSlist,
					as.numeric(logfile$parameters[[4]]),  # minpeak
					as.numeric(logfile$parameters[[5]]),  # drtsmall2
					as.numeric(logfile$parameters[[6]]),  # drtfill                
					as.numeric(logfile$parameters[[7]]),  # drtdens2
					as.numeric(logfile$parameters[[11]]), # recurs
					as.numeric(logfile$parameters[[13]]), # weight
					as.numeric(logfile$parameters[[10]]), # SB
					as.numeric(logfile$parameters[[9]]),  # SN               
					10^as.numeric(logfile$parameters[[8]]),  # minint
					10^as.numeric(logfile$parameters[[14]]), # maxint
					as.numeric(logfile$parameters[[12]]), # ended
					progbar=logfile$parameters[21],
					from=FALSE,to=FALSE
				);
				#MSlist[[8]]<-cbind(MSlist[[8]],MSlist[[8]][,c(1,4,5)]);
				save(MSlist,file=file.path(logfile[[1]],"MSlist",as.character(measurements[i,1])));   
				peaklist<-MSlist[[8]];
				peaklist<-cbind(peaklist,peaklist[,c(1,4,5)])
				colnames(peaklist)[12]<-"m/z_corr";
				colnames(peaklist)[13]<-"sum_int_corr";
				colnames(peaklist)[14]<-"RT_corr";      
				keep<-rep(1,length(peaklist[,1])) # 1 == TRUE
				peaklist<-cbind(peaklist,keep,keep) # replicate & blind indicators
				colnames(peaklist)[15]<-"keep";
				colnames(peaklist)[16]<-"keep_2";
				save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(measurements[i,1])));   
				cat(" plot -");  
				path=file.path(logfile[[1]],"pics",paste("peakhist_",as.character(measurements[i,1]),sep=""))
				png(filename = path, bg = "white")    
				hist(log10(MSlist[[4]][[2]][,2]),breaks=200,xlab="log10(Intensity)",main="All measurements (black) vs. measurements in peaks (red)")
				hist(log10(MSlist[[4]][[2]][MSlist[[4]][[2]][,7]!=0,2]),breaks=200,col="red",add=TRUE)        
				dev.off() 
				path=file.path(logfile[[1]],"pics",paste("peakmzRT_",as.character(measurements[i,1]),sep=""))
				png(filename = path, bg = "white")    
				plot(MSlist[[8]][,1],MSlist[[8]][,5],xlab="m/z",ylab="RT",pch=19,cex=0.3,main="Picked peaks")
				dev.off()
				if(any(objects(envir=as.environment(".GlobalEnv"))=="MSlist")){rm(MSlist,envir=as.environment(".GlobalEnv"))}
				if(any(objects()=="MSlist")){rm(MSlist)}
				measurements[i,10]<-TRUE;
				output$measurements<-DT::renderDataTable(read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character")); 
				write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
				cat(" picked -"); 
				formulator<-cbind(
					peaklist[,c(1,3)],rep(0,length(peaklist[,c(1)])),rep(0,length(peaklist[,c(1)])),
					peaklist[,c(5)]/60,rep(0,length(peaklist[,c(1)])),peaklist[,c(6)]/60,rep(0,length(peaklist[,c(1)])),
					peaklist[,c(7)]/60,rep(0,length(peaklist[,c(1)])))
				colnames(formulator)<-c(
					"Centroid m/z","Intensity","Peak S/N","Scan Number","RT (min.)","Start Scan Number",
					"Start RT (min.)","End Scan Number","End RT (min.)","Mass Chromatogram S/N")
				write.table(formulator,file=file.path(logfile[[1]],"exports",paste(as.character(measurements[i,2]),".txt",sep="")),row.names=FALSE,sep="\t",quote=FALSE);				
				if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
				if(any(objects()=="peaklist")){rm(peaklist)}
				cat(" exported. ");
            }
    }
    measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	output$measurements<-DT::renderDataTable(read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character")); 
	if(addit){library(nlme)}
    cat("Peak picking completed \n"); 	  
    updateSelectInput(session, "sel_meas_ID", label = "Select file by ID:", choices =  c("none",as.character(measurements[,1])), selected = "none")


