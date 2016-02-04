
path=file.path(logfile[[1]],"pics","plotQCa_pos")
				png(filename = path, bg = "white")
				plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
				dev.off()
				expr1p<-list(src=path)
				output$plotQCa_pos<-renderImage(expr1p, deleteFile = FALSE)	
path=file.path(logfile[[1]],"pics","plotQCb_pos")
				png(filename = path, bg = "white")
				plot.new();plot.window(xlim=c(1,1),ylim=c(1,1));#box();text(1,1,label="not available",cex=1.5,col="darkred")
				dev.off()
				expr2p<-list(src=path)
				output$plotQCb_pos<-renderImage(expr2p, deleteFile = FALSE)
path=file.path(logfile[[1]],"pics","plotQCa_neg")
				png(filename = path, bg = "white")
				plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
				dev.off()
				expr1n<-list(src=path)
				output$plotQCa_neg<-renderImage(expr1n, deleteFile = FALSE)	
path=file.path(logfile[[1]],"pics","plotQCb_neg")
				png(filename = path, bg = "white")
				plot.new();plot.window(xlim=c(1,1),ylim=c(1,1));#box();text(1,1,label="not available",cex=1.5,col="darkred")
				dev.off()
				expr2n<-list(src=path)
				output$plotQCb_neg<-renderImage(expr2n, deleteFile = FALSE)
logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="QC"]<-"FALSE";
logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="QC"]<<-"FALSE";
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp")); 
measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
measurements[,11]<-"skip";
write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);          
measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");






