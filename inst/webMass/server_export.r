##############################################################################
# Export peaklist selection ##################################################
##############################################################################
exported_1<-reactive({
    input$expo_profiles
    if( 
		(isolate(init$a)=="TRUE") & 
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
})	
output$expo1<-renderText(paste(exported_1()))
##############################################################################
