measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
measurements<-measurements[measurements[,8]=="TRUE",]
IDs<-list.files(file.path(logfile[[1]],"peaklist"))
filetype<-(measurements[,3])
ionmode<-(measurements[,4])
atdate<-(measurements[,6])
attime<-(measurements[,7])
attime2<-as.difftime(attime);
atdate<-as.Date(atdate);
sampleID<-(measurements[,1])
ord<-order(as.numeric(atdate),as.numeric(attime2),sampleID);
ppm<-logfile$parameters[[83]]
dmz<-as.numeric(logfile$parameters[[82]])
dRT<-as.numeric(logfile$parameters[[84]])
int_ratio<-as.numeric(logfile$parameters[[37]])

if(FALSE){ # debug parameters
	ppm<-TRUE
	dmz<-3
	dRT<-30
	int_ratio<-10
}

# clean old entries #####################################################################################################
if(length(IDs)>0){
for(i in 1:length(IDs)){
	if(any(measurements[,1]==IDs[i])){
		if(filetype[measurements[,1]==IDs[i]]=="sample"){
			load(file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])),envir=as.environment(".GlobalEnv"),verbose=FALSE);
			keep_2<-rep(1,length(peaklist[,1])) # 1 == TRUE
			peaklist[,colnames(peaklist)=="keep_2"]<-keep_2
			save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])))
			rm(peaklist)
		}
	}else{
		cat("\n Orphaned peaklist detected - from an older workflow run?")
	}
}
}

#######################################################################################################
# run last blank by date & time subtraction ###########################################################
if((logfile$parameters[[85]]=="TRUE") || (logfile$parameters[[87]]=="TRUE")){
	blank_ID_last<-"FALSE"
	for(i in 2:length(ord)){ # can skip first file
		if((logfile$parameters[[85]]=="FALSE") & (ionmode[ord[i]]=="positive")){next}
		if((logfile$parameters[[87]]=="FALSE") & (ionmode[ord[i]]=="negative")){next}	
		if(filetype[ord[i]]=="sample"){
			sam_ID<-sampleID[ord[i]]
			found_blank<-FALSE
			for(j in (i-1):1){ # backward
				if((filetype[ord[j]]=="blank") & (ionmode[ord[i]]==ionmode[ord[j]])){
					blank_ID<-sampleID[ord[j]]
					found_blank<-TRUE
					break;
				}
			}
			if(!found_blank){
				next;
			}
			if(blank_ID!=blank_ID_last){ # only load, not reload
				load(file=file.path(logfile[[1]],"peaklist",as.character(blank_ID)),verbose=FALSE);
				peaks_blank<-peaklist[,c(12,13,14)]
				rm(peaklist)
			}
			load(file=file.path(logfile[[1]],"peaklist",as.character(sam_ID)),verbose=FALSE);
			peaks_sample<-peaklist[,c(12,13,14)]		
			getit <- search_peak( 
				peaklist=peaks_blank, 
				mz=peaks_sample[,1], 
				dmz=dmz, # precheck for profiles
				ppm=ppm, 
				RT=peaks_sample[,3], 
				dRT=dRT,
				onlymax=TRUE,
				int_ratio=int_ratio,
				int=peaks_sample[,2],
				get_matches=FALSE
			)	
			peaklist[getit=="TRUE",colnames(peaklist)=="keep_2"]<-0
			save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(sam_ID)))
			blank_ID_last<-blank_ID
			cat(paste("\n",
				round(sum(peaklist[,colnames(peaklist)=="keep_2"]==0)/length(peaklist[,1])*100,digits=1),
				" % of ",
				length(peaklist[,1]),
				" peaks blind filtered (files ",
				sam_ID," vs. ",blank_ID,", ",ionmode[ord[i]],", by date & time)."
			,sep=""))
			rm(peaklist);
		}
	}
}
#######################################################################################################

#######################################################################################################
# run the selective subtraction #######################################################################
# positive ############################################################################################
if( (logfile$parameters[[86]]=="TRUE") & any(logfile$Positive_subtraction_files!="FALSE") ){
	selec_pos<-logfile$Positive_subtraction_files
	selec_pos<-selec_pos[selec_pos!="FALSE"]
	for(i in 1:length(IDs)){
		if(any(measurements[,1]==IDs[i])){
			if( filetype[measurements[,1]==IDs[i]]=="sample" &  ionmode[measurements[,1]==IDs[i]]=="positive" ){
				load(file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])),verbose=FALSE);
				sam_peaklist<-peaklist;rm(peaklist);
				for(j in 1:length(selec_pos)){
					subID<-strsplit(selec_pos[j]," - ")[[1]][1]
					load(file=file.path(logfile[[1]],"peaklist",as.character(subID)),verbose=FALSE)
					peaks_blank<-peaklist[,c(12,13,14)];rm(peaklist);
					getit <- search_peak( 
						peaklist=peaks_blank, 
						mz=sam_peaklist[,12], 
						dmz=dmz, # precheck for profiles
						ppm=ppm, 
						RT=sam_peaklist[,14], 
						dRT=dRT,
						onlymax=TRUE,
						int_ratio=int_ratio,
						int=sam_peaklist[,13],
						get_matches=FALSE
					)	
					sam_peaklist[getit=="TRUE",colnames(sam_peaklist)=="keep_2"]<-0
					rm(peaks_blank)
				}
				peaklist<-sam_peaklist
				save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])))
				cat(paste("\n",
					round(sum(peaklist[,colnames(peaklist)=="keep_2"]==0)/length(peaklist[,1])*100,digits=1),
					" % of ",
					length(peaklist[,1]),
					" peaks blind filtered (selective, file ",
					as.character(IDs[i]),"). "
				,sep=""))
				rm(peaklist,sam_peaklist);
			}
		}else{
			cat("\n Orphaned peaklist detected - from an older workflow run?")
		}
	}
}
#######################################################################################################

#######################################################################################################
# run the selective subtraction #######################################################################
# negative ############################################################################################
if( (logfile$parameters[[88]]=="TRUE") & any(logfile$Negative_subtraction_files!="FALSE") ){
	selec_neg<-logfile$Negative_subtraction_files
	selec_neg<-selec_neg[selec_neg!="FALSE"]
	for(i in 1:length(IDs)){
		if(any(measurements[,1]==IDs[i])){	
			if(filetype[measurements[,1]==IDs[i]]=="sample" &  ionmode[measurements[,1]==IDs[i]]=="negative"){
				load(file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])),verbose=FALSE);
				sam_peaklist<-peaklist;rm(peaklist);
				for(j in 1:length(selec_neg)){
					subID<-strsplit(selec_neg[j]," - ")[[1]][1]
					load(file=file.path(logfile[[1]],"peaklist",as.character(subID)),verbose=FALSE)
					peaks_blank<-peaklist[,c(12,13,14)];rm(peaklist);
					getit <- search_peak( 
						peaklist=peaks_blank, 
						mz=sam_peaklist[,12], 
						dmz=dmz, # precheck for profiles
						ppm=ppm, 
						RT=sam_peaklist[,14], 
						dRT=dRT,
						onlymax=TRUE,
						int_ratio=int_ratio,
						int=sam_peaklist[,13],
						get_matches=FALSE
					)	
					sam_peaklist[getit=="TRUE",colnames(sam_peaklist)=="keep_2"]<-0
					rm(peaks_blank)
				}
				peaklist<-sam_peaklist
				save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])))
				cat(paste("\n",
					round(sum(peaklist[,colnames(peaklist)=="keep_2"]==0)/length(peaklist[,1])*100,digits=1),
					" % of ",
					length(peaklist[,1]),
					" peaks blind filtered (selective, file ",
					as.character(IDs[i]),"). "
				,sep=""))
				rm(peaklist,sam_peaklist);
			}
		}else{
			cat("\n Orphaned peaklist detected - from an older workflow run?")
		}		
	}
}
#######################################################################################################





