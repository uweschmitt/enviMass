
measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
filetype<-(measurements[,3])
ionmode<-(measurements[,4])
atdate<-(measurements[,6])
attime<-(measurements[,7])
attime2<-as.difftime(attime);
atdate<-as.Date(atdate);
sampleID<-(measurements[,1])
ord<-order(as.numeric(atdate),as.numeric(attime2),sampleID);
ppm<-logfile$parameters[[83]]
dmz<-logfile$parameters[[82]]
dRT<-logfile$parameters[[84]]
int_ratio<-logfile$parameters[[37]]

if(FALSE){ # debug parameters
	ppm<-TRUE
	dmz<-3
	dRT<-30
	int_ratio<-10
}

blank_ID_last<-"FALSE"
for(i in 2:length(ord)){ # can skip first file
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
			load(file=file.path(logfile[[1]],"peaklist",as.character(blank_ID)),envir=as.environment(".GlobalEnv"),verbose=FALSE);
			peaks_blank<-peaklist[,c(12,13,14)]
			rm(peaklist)
		}
		load(file=file.path(logfile[[1]],"peaklist",as.character(sam_ID)),envir=as.environment(".GlobalEnv"),verbose=FALSE);
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
		#save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(sam_ID)))
		blank_ID_last<-blank_ID
		cat(paste("\n",
			round(sum(peaklist[,colnames(peaklist)=="keep_2"]==0)/length(peaklist[,1])*100,digits=1),
			" % of ",
			length(peaklist[,1]),
			" peaks blind filtered (files ",
			sam_ID," vs. ",blank_ID,", ",ionmode[ord[i]],")."
		,sep=""))
		rm(peaklist);
	}

}











