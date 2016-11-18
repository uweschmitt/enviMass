
use_int<-logfile$parameters$peak_which_intensity

			cat("Resetting intensities \n");
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			measurements<-measurements[measurements[,"include"]=="TRUE",,drop=FALSE]
			leng<-length(measurements[,"include"])
			for(i in 1:leng){
				if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
				if(any(objects()=="peaklist")){rm(peaklist)}
				if(!file.exists(file.path(logfile[[1]],"peaklist",as.character(measurements[i,"ID"])))){
					next;
				}
				load(file=file.path(logfile[[1]],"peaklist",as.character(measurements[i,"ID"])),envir=as.environment(".GlobalEnv"));
				if(use_int=="maximum"){
					peaklist[,13]<-peaklist[,3];
				}else{
					peaklist[,13]<-peaklist[,4];				
				}
				save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(measurements[i,"ID"])))
			}
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
			rm(measurements)
