##############################################################################
# Export profile list selection ##############################################
##############################################################################
exported_1<-reactive({
    input$expo_profiles
	if(isolate(init$a)=="TRUE"){
		if( 
			isolate(input$expo_profiles) &
			any(objects(envir=as.environment(".GlobalEnv"))=="profileList") & 
			any(objects(envir=as.environment(".GlobalEnv"))=="profpeaks") 	& 
			any(objects(envir=as.environment(".GlobalEnv"))=="profpeaks2")	
		){
			export<-enviMass:::export_profile_list(
				profileList,	
				profpeaks2,
				progbar=FALSE
			)
			getnames<-c("Date_time","sample_ID","blank_ID")
			leng<-((length(export[1,])-3)/2)
			for(i in 1:leng){
				getnames<-c(getnames,as.character(i),"-")
			}
			pathed<-file.path(logfile[[1]],"exports","profile_list.txt")
			write.table(export,file=pathed,col.names = getnames,row.names=FALSE)
			return("Latest export finished")	
		}else{
			return("Latest export failed")
		}
	}
})	
output$expo1<-renderText(paste(exported_1()))
##############################################################################

##############################################################################
# Export peaklist selection ##################################################
##############################################################################
exported_2<-reactive({
    input$expo_peaklist
	peaklistID<-isolate(input$sel_meas)
	if(isolate(init$a)=="TRUE"){
		if(
			(peaklistID!=0)&(file.exists(file.path(logfile[[1]],"peaklist",as.character(peaklistID))))
		){
			load(file=file.path(logfile[[1]],"peaklist",as.character(peaklistID)),envir=as.environment(".GlobalEnv"),verbose=FALSE);
			peaklist<-peaklist[
				peaklist[,colnames(peaklist)=="keep"]==1
			,,drop=FALSE]
			peaklist<-peaklist[
				peaklist[,colnames(peaklist)=="keep_2"]==1
			,,drop=FALSE]		
			peaklist<-peaklist[,c(12,13,14),drop=FALSE]
			write.csv(
				peaklist,
				file=file.path(logfile[[1]],"exports",paste("peaklist",as.character(peaklistID),".csv",sep="")),
				row.names=FALSE
			);
			rm(peaklist)	
			return("Latest export finished")	
		}else{
			return("No export conducted")
		}
	}
})	
output$expo2<-renderText(paste(exported_2()))
##############################################################################




























