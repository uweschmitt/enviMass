if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_obs_screening.r!")}

###########################################################################################################
# SPECIFY IONIZATION MODE #################################################################################
observe({ 
	input$Ion_mode_Cal 
	init$b
	if(isolate(init$a)=="TRUE"){
		if(isolate(input$Ion_mode_Cal)=="positive"){
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			measurements<-measurements[measurements[,3]=="calibration",,drop=FALSE]
			measurements<-measurements[measurements[,4]=="positive",,drop=FALSE]
			if(length(measurements[,1])>0){
				those<-unique(measurements$tag2)
				if(all(those!="FALSE")){
					those<-c("none",those)
					updateSelectInput(session,"Cal_file_set","Specify calibration file group",choices = those, selected = those[1])
				}else{ # trigger warning
					cat("all calibration groups must have a tag2 other than FALSE!")
				}
			}
		}
		if(isolate(input$Ion_mode_Cal)=="negative"){
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			measurements<-measurements[measurements[,3]=="calibration",,drop=FALSE]
			measurements<-measurements[measurements[,4]=="negative",,drop=FALSE]
			if(length(measurements[,1])>0){
				those<-unique(measurements$tag2)
				if(all(those!="FALSE")){
					those<-c("none",those)
					updateSelectInput(session,"Cal_file_set","Specify calibration file group",choices = those, selected = those[1])
				}else{ # trigger warning
					cat("all calibration groups must have a tag2 other than FALSE!")
				}
			}
		}
	}
})
###########################################################################################################

###########################################################################################################
# SPECIFY CALIBRATION GROUP ###############################################################################
observe({ 
	input$Cal_file_set
	init$b
	if(isolate(init$a)=="TRUE"){
		if(isolate(input$Cal_file_set)!="none"){
		


		}
	}
})
###########################################################################################################

  
 if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_obs_screening.r!")}