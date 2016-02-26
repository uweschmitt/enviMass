
recal_gams<-list.files(file.path(logfile[[1]],"results","recalibration"))
if(length(recal_gams)>0){
	for(i in 1:length(recal_gams)){
		file.remove(file.path(logfile[[1]],"results","recalibration",recal_gams[i])) 
	}
}

path=file.path(logfile[[1]],"pics","recal_none")
png(filename = path, bg = "white")
plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
dev.off()
exprrec<-list(src=path)
output$recal_pic<-renderImage(exprrec, deleteFile = FALSE);		
output$peakhist_pic<-renderImage(exprrec, deleteFile = FALSE);
output$peakmzRT_pic<-renderImage(exprrec, deleteFile = FALSE);		

IDs<-list.files(file.path(logfile[[1]],"peaklist"))
if(length(IDs)>0){
for(i in 1:length(IDs)){
	load(file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])),verbose=FALSE);
	peaklist[,12]<-peaklist[,1]		
	save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])))
	rm(peaklist)
}
}

