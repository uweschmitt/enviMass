use_int<-"max_int"
#use_int<-"sum_int"

if(
	( logfile$workflow[4]=="yes" && logfile$summary[7,2]=="FALSE" ) || ( logfile$Tasks_to_redo[4]=="TRUE" ) 
){

	measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	measurements<-measurements[measurements[,8]=="TRUE",]
	# positive #################################################################
	if(any(measurements[,4]=="positive")){
		measurements_pos<-measurements[measurements[,4]=="positive",]
		leng<-length(measurements_pos[,8])
		meanint<-c();
		maxint<-c();
		minint<-c();
		for(i in 1:leng){
			if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="peaklist")){rm(peaklist)}
			load(file=file.path(logfile[[1]],"peaklist",as.character(measurements_pos[i,1])),envir=as.environment(".GlobalEnv"));
			if(use_int=="max_int"){
				meanint<-c(meanint,median(peaklist[,3]));
				maxint<-c(maxint,max(peaklist[,3]));
				minint<-c(minint,min(peaklist[,3]));
			}else{
				meanint<-c(meanint,median(peaklist[,4]));
				maxint<-c(maxint,max(peaklist[,4]));
				minint<-c(minint,min(peaklist[,4]));
			}
			rm(peaklist,envir=as.environment(".GlobalEnv"))
		}
		atmean<-median(meanint)
		png(filename = file.path(logfile[[1]],"pics","int_distr_pos"), bg = "white",width = 680, height = 480)
		plot.new();
		plot.window(xlim=c(0,(leng+1)),ylim=c(log10(min(minint)),log10(max(maxint))))
		title(xlab="Measurement sequence (w/o outliers)",ylab="log10 Intensity",main="sample:green / blank:blue / other:grey")
		box();axis(1);axis(2);
		for(i in 1:leng){
			if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="peaklist")){rm(peaklist)}
			load(file=file.path(logfile[[1]],"peaklist",as.character(measurements_pos[i,1])),envir=as.environment(".GlobalEnv"));
			doneit<-FALSE
			if(measurements_pos[i,3]=="sample"){
				abline(h=log10(atmean),col="red",lty=2)
				if(use_int=="max_int"){
					boxplot(log10(peaklist[,3]),add=TRUE,at=i,cex=0.5,col="green")
				}else{
					boxplot(log10(peaklist[,4]),add=TRUE,at=i,cex=0.5,col="green")				
				}
				doneit<-TRUE
			}
			if(measurements_pos[i,3]=="blank"){
				abline(h=log10(atmean),col="red",lty=2)
				if(use_int=="max_int"){
					boxplot(log10(peaklist[,3]),add=TRUE,at=i,cex=0.5,col="blue")
				}else{
					boxplot(log10(peaklist[,4]),add=TRUE,at=i,cex=0.5,col="blue")				
				}
				doneit<-TRUE
			}
			if(doneit==FALSE){
				abline(h=log10(atmean),col="red",lty=2)
				if(use_int=="max_int"){
					boxplot(log10(peaklist[,3]),add=TRUE,at=i,cex=0.5,col="lightgrey")
				}else{
					boxplot(log10(peaklist[,4]),add=TRUE,at=i,cex=0.5,col="lightgrey")				
				}				
				doneit<-TRUE
			}
			if(use_int=="max_int"){
				peaklist[,13]<-c(peaklist[,3]/median(peaklist[,3])*atmean)
			}else{
				peaklist[,13]<-c(peaklist[,4]/median(peaklist[,4])*atmean)			
			}
			save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(measurements_pos[i,1])))
			rm(peaklist,envir=as.environment(".GlobalEnv"))
		}
		dev.off()
		expr1p<-list(src=file.path(logfile[[1]],"pics","int_distr_pos"))
		output$pic_int_distr_pos<-renderImage(expr1p, deleteFile = FALSE)
	}	
	# negative #################################################################	
	if(any(measurements[,4]=="negative")){
		measurements_neg<-measurements[measurements[,4]=="negative",]
		leng<-length(measurements_neg[,8])
		meanint<-c();
		maxint<-c();
		minint<-c();
		for(i in 1:leng){
			if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="peaklist")){rm(peaklist)}
			load(file=file.path(logfile[[1]],"peaklist",as.character(measurements_neg[i,1])),envir=as.environment(".GlobalEnv"));
			if(use_int=="max_int"){
				meanint<-c(meanint,median(peaklist[,3]));
				maxint<-c(maxint,max(peaklist[,3]));
				minint<-c(minint,min(peaklist[,3]));
			}else{	
				meanint<-c(meanint,median(peaklist[,4]));
				maxint<-c(maxint,max(peaklist[,4]));
				minint<-c(minint,min(peaklist[,4]));			
			}	
			rm(peaklist,envir=as.environment(".GlobalEnv"))
		}
		atmean<-median(meanint)
		png(filename = file.path(logfile[[1]],"pics","int_distr_neg"), bg = "white",width = 680, height = 480)
		plot.new();
		plot.window(xlim=c(0,(leng+1)),ylim=c(log10(min(minint)),log10(max(maxint))))
		title(xlab="Measurement sequence (w/o outliers)",ylab="log10 Intensity",main="sample:green / blank:blue / other:grey")
		box();axis(1);axis(2);
		for(i in 1:leng){
			if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="peaklist")){rm(peaklist)}
			load(file=file.path(logfile[[1]],"peaklist",as.character(measurements_neg[i,1])),envir=as.environment(".GlobalEnv"));
			doneit<-FALSE
			if(measurements_neg[i,3]=="sample"){
				abline(h=log10(atmean),col="red",lty=2)
				if(use_int=="max_int"){
					boxplot(log10(peaklist[,3]),add=TRUE,at=i,cex=0.5,col="green")
				}else{
					boxplot(log10(peaklist[,4]),add=TRUE,at=i,cex=0.5,col="green")
				}
				doneit<-TRUE
			}
			if(measurements_neg[i,3]=="blank"){
				abline(h=log10(atmean),col="red",lty=2)
				if(use_int=="max_int"){
					boxplot(log10(peaklist[,3]),add=TRUE,at=i,cex=0.5,col="blue")
				}else{
					boxplot(log10(peaklist[,4]),add=TRUE,at=i,cex=0.5,col="blue")
				}
				doneit<-TRUE
			}
			if(doneit==FALSE){
				abline(h=log10(atmean),col="red",lty=2)
				if(use_int=="max_int"){
					boxplot(log10(peaklist[,3]),add=TRUE,at=i,cex=0.5,col="lightgrey")
				}else{
					boxplot(log10(peaklist[,4]),add=TRUE,at=i,cex=0.5,col="lightgrey")
				}
				doneit<-TRUE
			}
			if(use_int=="max_int"){
				peaklist[,13]<-c(peaklist[,3]/median(peaklist[,3])*atmean)			
			}else{
				peaklist[,13]<-c(peaklist[,4]/median(peaklist[,4])*atmean)
			}
			save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(measurements_neg[i,1])))
			if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="peaklist")){rm(peaklist)}
		}
		dev.off()
		expr1n<-list(src=file.path(logfile[[1]],"pics","int_distr_neg"))
		output$pic_int_distr_neg<-renderImage(expr1n, deleteFile = FALSE)
	}	
	############################################################################
	logfile$summary[7,2]<<-"TRUE";
    logfile$summary[7,2]<-"TRUE";
	logfile$Tasks_to_redo[4]<-"FALSE";
	logfile$Tasks_to_redo[4]<<-"FALSE";
    save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
	summa[7,2]<-"done"
	summa[7,2]<<-"done"
	output$summar<<-renderTable(summa)
    cat("Normalization done \n");
    output$dowhat<<-renderText("Normalization done ... wait")

}else{
		
		# re-insert non-normalized intensities #######################################3
        if( logfile$workflow[4]=="no" && logfile$summary[7,2]=="TRUE" ){
			cat("Resetting intensities \n");
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			measurements<-measurements[measurements[,8]=="TRUE",]
			leng<-length(measurements[,8])
			for(i in 1:leng){
				if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
				if(any(objects()=="peaklist")){rm(peaklist)}
				load(file=file.path(logfile[[1]],"peaklist",as.character(measurements[i,1])),envir=as.environment(".GlobalEnv"));
				if(use_int=="max_int"){
					peaklist[,13]<-peaklist[,3];
				}else{
					peaklist[,13]<-peaklist[,4];				
				}
				save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(measurements[i,1])))
			}
		}
		if(logfile$workflow[4]=="no"){
			logfile$summary[7,2]<<-"FALSE";
			logfile$summary[7,2]<-"FALSE";
			path=file.path(logfile[[1]],"pics","int_distr_pos")
				png(filename = path, bg = "white")
				plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
				dev.off()
				expr3p<-list(src=path)
				output$pic_int_distr_pos<-renderImage(expr3p, deleteFile = FALSE)
				path=file.path(logfile[[1]],"pics","int_distr_neg")
			png(filename = path, bg = "white")
				plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
				dev.off()
				expr3n<-list(src=path)
				output$pic_int_distr_neg<-renderImage(expr3n, deleteFile = FALSE)
		}
		logfile$Tasks_to_redo[4]<-"FALSE";
	    logfile$Tasks_to_redo[4]<<-"FALSE";
        save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
        summa[7,2]<-"skipped";
        summa[7,2]<<-"skipped";
		output$summar<<-renderTable(summa)
        cat("Normalization skipped \n");
        output$dowhat<<-renderText("Normalization skipped ... wait")

}