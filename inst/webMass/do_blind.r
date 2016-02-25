measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
filetype<-(measurements[,3])
atdate<-(measurements[,6])
attime<-(measurements[,7])
attime2<-as.difftime(attime);
atdate<-as.Date(atdate);
sampleID<-(measurements[,1])
ord<-order(as.numeric(atdate),as.numeric(attime2),sampleID);

for(i in 2:length(ord)){ # can skip first file
	if(filetype[ord[i]]=="sample"){
		sam_ID<-sampleID[ord[i]]
		found_blank<-FALSE
		for(j in (i-1):1){ # backward
			if(filetype[ord[j]]=="blank"){
				blank_ID<-filetype[ord[j]]
				found_blank<-TRUE
			}
		}
		if(!found_blank){
			next;
		}

		

	}
}


