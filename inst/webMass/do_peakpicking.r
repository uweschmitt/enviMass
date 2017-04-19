# Peak picking ##################################################################

	output$dowhat<-renderText("Peak picking ... please wait");
	if(any(search()=="package:nlme")){detach(package:nlme,force=TRUE);addit<-TRUE}else{addit<-FALSE}
    measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
    leng<-dim(measurements)[1];         
    for(i in 1:leng){ 
        # (measurement included & not yet picked) OR (peakpick forced) 
            if( (measurements[i,"include"]=="TRUE") & (measurements[i,"peakpicking"]=="FALSE") ){

				##################################################################
				cat(paste("\n    Peak picking sample ",as.character(i)," of ",as.character(leng),": "));    
				if(logfile$parameters$cut_RT=="TRUE"){
					use_minRT<-(as.numeric(logfile$parameters$cut_RT_min)*60)
					use_maxRT<-(as.numeric(logfile$parameters$cut_RT_max)*60)
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
				##################################################################
				if(logfile$parameters$is_example=="FALSE"){
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
					cat(" file read -"); 
					##############################################################
					if(logfile$parameters$peak_estimate=="TRUE"){
						use_peak_perc_cut<-0
						estim_values<-try({enviMass:::dens_filter(MSlist,plotit=FALSE,n=2000,m=5)},silent=TRUE)
						if(class(estim_values)!="try-error"){
							use_peak_dmzdens<-estim_values[[1]]
							use_peak_minint_log10<-estim_values[[2]]
							if(as.numeric(logfile$parameters$peak_maxint_log10)<use_peak_minint_log10){
								use_peak_maxint_log10<-log10(max(MSlist[["Scans"]][[2]][,"intensity"])+1)
							}else{
								use_peak_maxint_log10<-as.numeric(logfile$parameters$peak_maxint_log10)
							} 
							len1<-dim(MSlist[["Scans"]][[2]])[1]
							MSlist[["Scans"]][[2]]<-MSlist[["Scans"]][[2]][
								MSlist[["Scans"]][[2]][,"intensity"]>=(10^use_peak_minint_log10)
							,]
							len2<-dim(MSlist[["Scans"]][[2]])[1]
							cat(paste(" ",as.character(len1-len2),"of",as.character(len1),"data points discarded by absolute threshold -"))
							cat(" data filtered -"); 
						}else{
							use_peak_dmzdens<-as.numeric(logfile$parameters$peak_dmzdens)
							use_peak_minint_log10<-as.numeric(logfile$parameters$peak_minint_log10)
							use_peak_maxint_log10<-as.numeric(logfile$parameters$peak_maxint_log10)					
							use_peak_perc_cut<-as.numeric(logfile$parameters$peak_perc_cut)
							cat(" no data filtering possible -"); 
						}	
					}else{
						use_peak_dmzdens<-as.numeric(logfile$parameters$peak_dmzdens)
						use_peak_minint_log10<-as.numeric(logfile$parameters$peak_minint_log10)
						use_peak_maxint_log10<-as.numeric(logfile$parameters$peak_maxint_log10)					
						use_peak_perc_cut<-as.numeric(logfile$parameters$peak_perc_cut)
						if(use_peak_perc_cut>0){ # not to be used with filtering estimates - that uses absolute threshold intensities
							len1<-dim(MSlist[["Scans"]][[2]])[1]
							MSlist[["Scans"]][[2]]<-MSlist[["Scans"]][[2]][
								MSlist[["Scans"]][[2]][,"intensity"]>=quantile(MSlist[["Scans"]][[2]][,"intensity"],(use_peak_perc_cut/100))
							,]
							len2<-dim(MSlist[["Scans"]][[2]])[1]
							cat(paste(" ",as.character(len1-len2),"of",as.character(len1),"data points discarded by fraction -"))
						}
					}	
					##############################################################
				}else{ # no mzXML files for example projects -> use MSlist
					load(file=file.path(logfile[[1]],"MSlist",as.character(measurements[i,"ID"]))); 
					MSlist[["Scans"]][[2]]<-MSlist[["Scans"]][[2]][ # re-set
						order(MSlist[["Scans"]][[2]][,"RT"],decreasing=FALSE)
					,]
					MSlist[["Scans"]][[2]][,"partID"]<-0
					MSlist[["Scans"]][[2]][,"clustID"]<-0					
					MSlist[["Scans"]][[2]][,"peakID"]<-0
					use_peak_dmzdens<-as.numeric(logfile$parameters$peak_dmzdens)
					use_peak_minint_log10<-0
					use_peak_maxint_log10<-as.numeric(logfile$parameters$peak_maxint_log10)					
					use_peak_perc_cut<-0
					cat("example MSlist loaded -"); 
				}
				##################################################################
				if(any(MSlist[["Scans"]][[2]][,"intensity"]==0)){
					cat("\n Note in do_peakpicking: zero intensities found and discarded.")
					MSlist[["Scans"]][[2]]<-MSlist[["Scans"]][[2]][MSlist[["Scans"]][[2]][,"intensity"]!=0,,drop=FALSE]
				}	
				MSlist<-enviPick::mzagglom(
					MSlist,
					((use_peak_dmzdens*2)+1),
					ppm=TRUE,  
					as.numeric(logfile$parameters$peak_drtgap),
					as.numeric(logfile$parameters$peak_minpeak),
					10^use_peak_maxint_log10,
					progbar=as.logical(logfile$parameters$progressBar)
				);
				cat(" partitioned -");
				##################################################################
				MSlist<-enviPick::mzclust(      
					MSlist,
					use_peak_dmzdens,
					ppm=TRUE,
					60,
					as.numeric(logfile$parameters$peak_minpeak),
					10^use_peak_maxint_log10,
					progbar=as.logical(logfile$parameters$progressBar),
					merged=TRUE,from=FALSE,to=FALSE
				);
				cat(" clustered -");
				##################################################################
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
					10^use_peak_minint_log10,  # minint
					10^use_peak_maxint_log10, # maxint
					as.numeric(logfile$parameters$peak_ended), # ended
					progbar=logfile$parameters$progressBar,
					from=FALSE,to=FALSE
				);
				if(any(MSlist[["Peaklist"]][,3]==0)){stop("\n do_peakpicking: zero intensities found - resolve issue before proceding.")}
				save(MSlist,file=file.path(logfile[[1]],"MSlist",as.character(measurements[i,"ID"])));   
				peaklist<-MSlist[["Peaklist"]];
				if(length(peaklist)==0){
					stop("No peaks picked - wrong parameters (e.g., intensity thresholds too high)?")
				}
				##################################################################
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
				a<-hist(log10(MSlist[["Scans"]][[2]][,2]),breaks=200,plot=FALSE)
				hist(log10(MSlist[["Scans"]][[2]][,2]),breaks=a$breaks,
					xlim=c(min(a$breaks),max(a$breaks)),ylim=c(0,max(a$counts)),
					xlab="log10(Intensity)",main="All data points (white histrogram), those in peaks (red) and their count ratio (blue points)",cex.main=.8)
				b<-hist(log10(MSlist[["Scans"]][[2]][MSlist[["Scans"]][[2]][,7]!=0,2]),breaks=a$breaks,col="red",add=TRUE)        
				atfrac<-(b$counts/a$counts)
				plot.window(xlim=c(min(a$breaks),max(a$breaks)),ylim=c(0,1))
				points(a$mids[!is.na(atfrac)],atfrac[!is.na(atfrac)],col="blue",cex=.5,pch=19)
				axis(4)
				dev.off() 				
				path=file.path(logfile[[1]],"pics",paste("peakmzRT_",as.character(measurements[i,1]),sep=""))
				png(filename = path, bg = "white")    
				plot(MSlist[["Peaklist"]][,1],MSlist[["Peaklist"]][,5],xlab="m/z",ylab="RT",pch=19,cex=0.3,main="Picked peaks")
				dev.off()
				if(any(objects(envir=as.environment(".GlobalEnv"))=="MSlist")){rm(MSlist,envir=as.environment(".GlobalEnv"))}
				if(any(objects()=="MSlist")){rm(MSlist)}
				measurements[i,"peakpicking"]<-TRUE;
				write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
				cat(" picked."); 
				##################################################################

            }
    }
	if(addit){library(nlme)}
    cat("Peak picking completed \n"); 	  
    

