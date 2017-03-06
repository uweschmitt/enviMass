# Peak picking ##################################################################

	output$dowhat<-renderText("Peak picking ... please wait");
	if(any(search()=="package:nlme")){detach(package:nlme,force=TRUE);addit<-TRUE}else{addit<-FALSE}
    measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
    leng<-length(measurements[,"ID"]);         
    for(i in 1:leng){ 
        # (measurement included & not yet picked) OR (peakpick forced) 
            if( (measurements[i,"include"]=="TRUE")&&(measurements[i,"peakpicking"]=="FALSE") ){
				cat(paste("\n    Peak picking sample ",as.character(i)," of ",as.character(leng),": "));    
				if(logfile$parameters$cut_RT=="TRUE"){
					use_minRT<-as.numeric(logfile$parameters$cut_RT_min)
					use_maxRT<-as.numeric(logfile$parameters$cut_RT_max)
					cat("(filter RT range)")
				}else{
					use_minRT<-FALSE
					use_maxRT<-FALSE				
				}		
				if(logfile$parameters$cut_mass=="TRUE"){
					use_minmass<-as.numeric(logfile$parameters$cut_mass_min)
					use_maxmass<-as.numeric(logfile$parameters$cut_mass_max)
					cat("(filter mass range)")
				}else{
					use_minmass<-FALSE
					use_maxmass<-FALSE				
				}				
				MSlist<-enviPick::readMSdata(
					file.path(logfile[[1]],"files",paste(as.character(measurements[i,1]),".mzXML",sep="")),
					MSlevel=logfile$parameters$peak_MSlevel,  # MSlevel
					progbar=logfile$parameters$progressBar, # progbar
					minRT=use_minRT,
					maxRT=use_maxRT,
					minmz=use_minmass,
					maxmz=use_maxmass,
					ion_mode=measurements[i,"Mode"]
				);
				if(any(MSlist[[4]][[2]][,2]==0)){
					cat("\n do_peakpicking: zero intensities found - resolve issue before proceding.")
					MSlist[[4]][[2]]<-MSlist[[4]][[2]][MSlist[[4]][[2]][,2]!=0,,drop=FALSE]
				}				
				cat(" data read -");  				
				if(as.numeric(logfile$parameters$peak_perc_cut)>0){
					len1<-length(MSlist[[4]][[2]][,2])
					MSlist[[4]][[2]]<-MSlist[[4]][[2]][
						MSlist[[4]][[2]][,2]>=quantile(MSlist[[4]][[2]][,2],(as.numeric(logfile$parameters$peak_perc_cut)/100))
					,]
					len2<-length(MSlist[[4]][[2]][,2])
					cat(paste(" ",as.character(len1-len2),"of",as.character(len1),"data points discarded -"))
				}
				MSlist<-enviPick::mzagglom(
					MSlist,
					((as.numeric(logfile$parameters$peak_dmzdens)*2)+1),
					ppm=TRUE,  
					as.numeric(logfile$parameters$peak_drtgap),
					as.numeric(logfile$parameters$peak_minpeak),
					10^as.numeric(logfile$parameters$peak_maxint_log10),
					progbar=as.logical(logfile$parameters$progressBar)
				);
				cat(" partitioned -");
				MSlist<-enviPick::mzclust(      
					MSlist,
					as.numeric(logfile$parameters$peak_dmzdens),
					ppm=TRUE,
					60,
					as.numeric(logfile$parameters$peak_minpeak),
					10^as.numeric(logfile$parameters$peak_maxint_log10),
					progbar=as.logical(logfile$parameters$progressBar),
					merged=TRUE,from=FALSE,to=FALSE
				);
				cat(" clustered -");
				MSlist<-enviPick::mzpick(       
					MSlist,
					as.numeric(logfile$parameters$peak_minpeak),  # minpeak
					as.numeric(logfile$parameters$peak_drtsmall2),  # drtsmall2
					as.numeric(logfile$parameters$peak_drtfill),  # drtfill                
					as.numeric(logfile$parameters$peak_drtdens2),  # drtdens2
					as.numeric(logfile$parameters$peak_recurs), # recurs
					as.numeric(logfile$parameters$peak_weight), # weight
					as.numeric(logfile$parameters$peak_SB), # SB
					as.numeric(logfile$parameters$peak_SN),  # SN               
					10^as.numeric(logfile$parameters$peak_minint_log10),  # minint
					10^as.numeric(logfile$parameters$peak_maxint_log10), # maxint
					as.numeric(logfile$parameters$peak_ended), # ended
					progbar=logfile$parameters$progressBar,
					from=FALSE,to=FALSE
				);
				if(any(MSlist[[8]][,3]==0)){stop("\n do_peakpicking: zero intensities found - resolve issue before proceding.")}
				#MSlist[[8]]<-cbind(MSlist[[8]],MSlist[[8]][,c(1,4,5)]);
				save(MSlist,file=file.path(logfile[[1]],"MSlist",as.character(measurements[i,"ID"])));   
				peaklist<-MSlist[[8]];
				if(length(peaklist)==0){
					stop("No peaks picked - wrong parameters (e.g., intensity thresholds too high)?")
				}
				peaklist<-cbind(peaklist,peaklist[,1],rep(0,length(peaklist[,4])),peaklist[,5])
				colnames(peaklist)[12]<-"m/z_corr";
				colnames(peaklist)[13]<-"int_corr";
				colnames(peaklist)[14]<-"RT_corr";      
				keep<-rep(1,length(peaklist[,1])) # 1 == TRUE
				peaklist<-cbind(peaklist,keep,keep) # replicate & blind indicators
				colnames(peaklist)[15]<-"keep";
				colnames(peaklist)[16]<-"keep_2";
				save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(measurements[i,"ID"])));   
				cat(" plotted -");  
				path=file.path(logfile[[1]],"pics",paste("peakhist_",as.character(measurements[i,"ID"]),sep=""))
				png(filename = path, bg = "white")    
				a<-hist(log10(MSlist[[4]][[2]][,2]),breaks=200,plot=FALSE)
				hist(log10(MSlist[[4]][[2]][,2]),breaks=a$breaks,
					xlim=c(min(a$breaks),max(a$breaks)),ylim=c(0,max(a$counts)),
					xlab="log10(Intensity)",main="All data points (white histrogram), those in peaks (red) and their count ratio (blue points)",cex.main=.8)
				b<-hist(log10(MSlist[[4]][[2]][MSlist[[4]][[2]][,7]!=0,2]),breaks=a$breaks,col="red",add=TRUE)        
				atfrac<-(b$counts/a$counts)
				plot.window(xlim=c(min(a$breaks),max(a$breaks)),ylim=c(0,1))
				points(a$mids[!is.na(atfrac)],atfrac[!is.na(atfrac)],col="blue",cex=.5,pch=19)
				axis(4)
				dev.off() 				
				path=file.path(logfile[[1]],"pics",paste("peakmzRT_",as.character(measurements[i,1]),sep=""))
				png(filename = path, bg = "white")    
				plot(MSlist[[8]][,1],MSlist[[8]][,5],xlab="m/z",ylab="RT",pch=19,cex=0.3,main="Picked peaks")
				dev.off()
				if(any(objects(envir=as.environment(".GlobalEnv"))=="MSlist")){rm(MSlist,envir=as.environment(".GlobalEnv"))}
				if(any(objects()=="MSlist")){rm(MSlist)}
				measurements[i,"peakpicking"]<-TRUE;
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
	if(addit){library(nlme)}
    cat("Peak picking completed \n"); 	  
    

