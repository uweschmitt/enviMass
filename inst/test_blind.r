########################################################################
# On blind/blank subtraction ###########################################
cat("\n Checking blind subtraction ...")
sam<-1000
if(
	(logfile$workflow[names(logfile$workflow)=="blind"])=="yes"
){
	cat(" included ...")
	####################################################################
	ppm<-logfile$parameters$blind_ppm
	dmz<-as.numeric(logfile$parameters$blind_dmz)
	dRT<-as.numeric(logfile$parameters$blind_drt)
	int_ratio<-as.numeric(logfile$parameters$blind_threshold)
	####################################################################
	measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	measurements<-measurements[measurements[,"include"]=="TRUE",]
	
	if(
		(any(measurements$Type=="sample")) & (any(measurements$Type=="blank"))
	){
		cat(" files available ...")
		atdate<-(measurements[,"Date"])
		attime<-(measurements[,"Time"])
		attime2<-as.difftime(attime);
		atdate<-as.Date(atdate);
		sampleID<-(measurements[,"ID"])
		ord<-order(as.numeric(atdate),as.numeric(attime2),sampleID);	
		measurements<-measurements[ord,]
		measurements_sam<-measurements[measurements$Type=="sample",]
		IDsam<-measurements_sam$ID
		measurements_blank<-measurements[measurements$Type=="blank",]
		IDblank<-measurements_blank$ID		
		for(i in 1:length(IDsam)){
			cat("\n - ")
			############################################################
			load(file=file.path(logfile[[1]],"peaklist",as.character(IDsam[i])),envir=as.environment(".GlobalEnv"),verbose=FALSE);
			peaks_sam<-peaklist;rm(peaklist)
			############################################################
			# check peaks found in blinds = must be present in #########
			# ANY blind measurement	####################################
			if(any(peaks_sam[,colnames(peaks_sam)=="keep_2"]==0)){
				those<-which(peaks_sam[,colnames(peaks_sam)=="keep_2"]==0)
				if(length(those)>sam){
					those<-sample(those,sam,replace=FALSE)
				}
				found<-rep(FALSE,length(those))
				for(k in 1:length(IDblank)){
					if(measurements_sam[measurements_sam$ID==IDsam[i],"Mode"]!=measurements_blank[measurements_blank$ID==IDblank[k],"Mode"]){next}
					cat("*")
					load(file=file.path(logfile[[1]],"peaklist",as.character(IDblank[k])),envir=as.environment(".GlobalEnv"),verbose=FALSE);
					peaks_blank<-peaklist;rm(peaklist)
					for(j in 1:length(those)){
						MASS<-peaks_sam[those[j],colnames(peaks_sam)=="m/z_corr"]
						RT<-peaks_sam[those[j],colnames(peaks_sam)=="RT_corr"]
						INT<-peaks_sam[those[j],colnames(peaks_sam)=="sum_int_corr"]
						if(ppm){
							mz_low<-(MASS-2*(dmz*MASS/1E6))
							mz_up<-(MASS+2*(dmz*MASS/1E6))
						}else{
							mz_low<-(MASS-2*(dmz/1000))
							mz_up<-(MASS+2*(dmz/1000))							
						}
						RT_low<-(RT-dRT)
						RT_up<-(RT+dRT)
						int_low<-(INT/int_ratio)
						got<-length(peaks_blank[
							(peaks_blank[,colnames(peaks_blank)=="m/z_corr"]>=mz_low) &
							(peaks_blank[,colnames(peaks_blank)=="m/z_corr"]<=mz_up) &
							(peaks_blank[,colnames(peaks_blank)=="RT_corr"]>=RT_low) &
							(peaks_blank[,colnames(peaks_blank)=="RT_corr"]<=RT_up) &							
							(peaks_blank[,colnames(peaks_blank)=="sum_int_corr"]>=int_low) 						
						,1])
						if(got>0){found[j]<-TRUE}
					}
				}
				if(any(!found)){cat("\nWARNING: wrongly non-subtracted peak(s) found!")}
			}
			############################################################
			# check peaks not found in blinds ##########################
			if(any(peaks_sam[,colnames(peaks_sam)=="keep_2"]==1)){
				those<-which(peaks_sam[,colnames(peaks_sam)=="keep_2"]==1)
				if(length(those)>sam){
					those<-sample(those,sam,replace=FALSE)
				}
				found<-rep(FALSE,length(those))
				for(k in 1:length(IDblank)){	
					if(measurements_sam[measurements_sam$ID==IDsam[i],"Mode"]!=measurements_blank[measurements_blank$ID==IDblank[k],"Mode"]){next} # same mode?
					skip<-TRUE	
					selec_pos<-logfile$Positive_subtraction_files
					selec_pos<-selec_pos[selec_pos!="FALSE"]
					subID<-strsplit(selec_pos," - ")
					subID<-lapply(subID, `[[`, 1)
					subID<-unlist(subID)
					if( # selective
						(measurements_sam[measurements_sam$ID==IDsam[i],"Mode"]=="positive")&
						(logfile$parameters$subtract_pos_byfile=="TRUE") &
						(any(subID==IDblank[k]))
					){
						skip<-FALSE
					}
					if( # by date/time order
						((measurements_sam[measurements_sam$ID==IDsam[i],"Mode"]=="positive")&(logfile$parameters$subtract_pos_bydate=="TRUE")) ||
						((measurements_sam[measurements_sam$ID==IDsam[i],"Mode"]=="negative")&(logfile$parameters$subtract_neg_bydate=="TRUE")) 
					){					
						at<-which(measurements[,"ID"]==IDsam[i])
						if(at>1){
							for(m in (at-1):1){ # backward
								if((measurements[m,"Type"]=="blank") & (measurements[m,"Mode"]==measurements_sam[measurements_sam$ID==IDsam[i],"Mode"])){
									break;
								}
							}
							if(measurements[m,"ID"]==IDblank[k]){
								skip<-FALSE
							}
						}						
					}
					if(skip){next}		
					cat(".")
					load(file=file.path(logfile[[1]],"peaklist",as.character(IDblank[k])),envir=as.environment(".GlobalEnv"),verbose=FALSE);
					peaks_blank<-peaklist;rm(peaklist)
					for(j in 1:length(those)){
						MASS<-peaks_sam[those[j],colnames(peaks_sam)=="m/z_corr"]
						RT<-peaks_sam[those[j],colnames(peaks_sam)=="RT_corr"]
						INT<-peaks_sam[those[j],colnames(peaks_sam)=="sum_int_corr"]
						if(ppm){
							mz_low<-(MASS-2*(dmz*MASS/1E6))
							mz_up<-(MASS+2*(dmz*MASS/1E6))
						}else{
							mz_low<-(MASS-2*(dmz/1000))
							mz_up<-(MASS+2*(dmz/1000))							
						}
						RT_low<-(RT-dRT)
						RT_up<-(RT+dRT)
						int_low<-(INT/int_ratio)
						got<-length(peaks_blank[
							(peaks_blank[,colnames(peaks_blank)=="m/z_corr"]>=mz_low) &
							(peaks_blank[,colnames(peaks_blank)=="m/z_corr"]<=mz_up) &
							(peaks_blank[,colnames(peaks_blank)=="RT_corr"]>=RT_low) &
							(peaks_blank[,colnames(peaks_blank)=="RT_corr"]<=RT_up) &							
							(peaks_blank[,colnames(peaks_blank)=="sum_int_corr"]>=int_low) 						
						,1])
						if(got>0){found[j]<-TRUE}
					}
				}				
				if(any(found)){cat("\nWARNING: wrongly subtracted peak(s) found!")}				
			}
			############################################################	
			rm(peaks_sam)
			############################################################
		}
	}
	rm(measurements)
	####################################################################
}
cat("\n Checking blind subtraction ... done.")
########################################################################


