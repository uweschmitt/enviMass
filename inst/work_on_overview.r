




ui <- fluidPage(
		helpText("The below plot indicates available files as dots at their respective date and time, listed over the different file categories 
		and seperately for each of the two ion modes."),
		helpText("Select a time period via mouse brush to list the file IDs for it. Double-click into the area to zoom in; double-click again to zoom out."),
		plotOutput("file_overview", 
			dblclick = "file_overview_dblclick",
			brush = brushOpts(
				id = "file_overview_brush",
				resetOnNew = TRUE,
				direction = "x",
				fill="red"
			),
			height = "600px"
		),
		HTML('<hr noshade="noshade" />'),
		HTML('<br><p><font size="5"> + </font><font size="3"> Selected file IDs, positive ionization:</p></font>'),
		htmlOutput("info_files_pos_samp"),htmlOutput("info_files_pos_blind"),htmlOutput("info_files_pos_cal"),htmlOutput("info_files_pos_calgroup"),htmlOutput("info_files_pos_spiked"),
		HTML('<hr noshade="noshade" />'),		
		HTML('<br><p><font size="5"> - </font><font size="3"> Selected file IDs, negative ionization:</p></font>'),
		htmlOutput("info_files_neg_samp"),htmlOutput("info_files_neg_blind"),htmlOutput("info_files_neg_cal"),htmlOutput("info_files_neg_calgroup"),htmlOutput("info_files_neg_spiked")
)

server <- function(input, output) {

	ranges_overview <- reactiveValues(x = NULL, y = NULL)
	output$file_overview <- renderPlot({
		enviMass:::plot_measurements(logfile,ranges_overview)
	})
	####################################################
	observeEvent(input$file_overview_dblclick, {
		brush <- input$file_overview_brush
		if (!is.null(brush)) {
			ranges_overview$x <- c(brush$xmin, brush$xmax)
			ranges_overview$y <- NULL#c(brush$ymin, brush$ymax)
		} else {
			ranges_overview$x <- NULL
			ranges_overview$y <- NULL
		}
	})
	####################################################
	observeEvent(input$file_overview_brush, {
		brush <- input$file_overview_brush
		if (!is.null(brush)) {
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");	
			dated<-measurements[,"Date"]
			timed<-measurements[,"Time"]
			datetime<-c()
			for(i in 1:length(timed)){
				datetime<-c(datetime,paste(dated[i],timed[i],"CET",sep=" "))
			}
			atPOSIX<-as.POSIXct(datetime);	
			#############################################
			# positive, samples			
			these<-measurements[
				(measurements[,"Mode"]=="positive") &
				(measurements[,"Type"]=="sample") &
				(atPOSIX>=as.POSIXct(brush$xmin, origin = "1970-01-01")) &
				(atPOSIX<=as.POSIXct(brush$xmax, origin = "1970-01-01"))
			,"ID"]
			output$info_files_pos_samp <- renderText({
				if(length(these)>0){paste("<font color=\"black\"> Sample IDs: ",paste(these,collapse=", ")," </font>",sep="")}else{paste("<font color=\"black\"> No sample files selected </font>",sep="")}
			})
			# positive, blind	
			these2<-measurements[
				(measurements[,"Mode"]=="positive") &
				(measurements[,"Type"]=="blank") &
				(atPOSIX>=as.POSIXct(brush$xmin, origin = "1970-01-01")) &
				(atPOSIX<=as.POSIXct(brush$xmax, origin = "1970-01-01"))
			,"ID"]
			output$info_files_pos_blind<- renderText({
				if(length(these2)>0){paste("<font color=\"green\"> Blanks/blind IDs: ",paste(these2,collapse=", ")," </font>",sep="")}else{paste("<font color=\"green\"> No blind files selected </font>",sep="")}
			})
			# positive, calibration	
			if(any((measurements[,"Type"]=="calibration") & (measurements[,"Mode"]=="positive"))){
				these3<-which((measurements[,"Mode"]=="positive") &(measurements[,"Type"]=="calibration") )
				date_start<-measurements[these3,"Date"]
				date_end<-measurements[these3,"date_end"]
				time_start<-measurements[these3,"Time"]
				time_end<-measurements[these3,"time_end"]
				datetime_start<-c()
				datetime_end<-c()
				for(i in 1:length(time_start)){
					datetime_start<-c(datetime_start,paste(date_start[i],time_start[i],"CET",sep=" "))
					datetime_end<-c(datetime_end,paste(date_end[i],time_end[i],"CET",sep=" "))
				}
				atPOSIX_start<-as.POSIXct(datetime_start);	
				atPOSIX_end<-as.POSIXct(datetime_end);	
				these3<-these3[
					(atPOSIX_start<=as.POSIXct(brush$xmax, origin = "1970-01-01")) &
					(atPOSIX_end>=as.POSIXct(brush$xmin, origin = "1970-01-01"))
				]
				these4<-unique(measurements[these3,"tag2"])
				these3<-measurements[these3,"ID"]
			}else{
				these3<-c()
				these4<-c()
			}
			output$info_files_pos_cal<- renderText({
				if(length(these3)>0){
					paste("<font color=\"red\"> Calibration file IDs: ",paste(these3,collapse=", ")," </font>",sep="")
				}else{
					paste("<font color=\"red\"> No calibration files selected </font>",sep="")
				}
			})
			output$info_files_pos_calgroup<- renderText({
				if(length(these4)>0){
					paste("<font color=\"red\"> Calibration file groups: ",paste(these4,collapse=", ")," </font>",sep="")
				}else{
					paste("<font color=\"red\"> No calibration file groups selected </font>",sep="")
				}
			})		
			# positive, spiked	
			these5<-measurements[
				(measurements[,"Mode"]=="positive") &
				(measurements[,"Type"]=="spiked") &
				(atPOSIX>=as.POSIXct(brush$xmin, origin = "1970-01-01")) &
				(atPOSIX<=as.POSIXct(brush$xmax, origin = "1970-01-01"))
			,"ID"]
			output$info_files_pos_spiked<- renderText({
				if(length(these5)>0){paste("<font color=\"blue\"> Spiked file IDs: ",paste(these5,collapse=", ")," </font>",sep="")}else{paste("<font color=\"blue\"> No spiked files selected </font>",sep="")}
			})			
			#############################################
			# negative, samples			
			these6<-measurements[
				(measurements[,"Mode"]=="negative") &
				(measurements[,"Type"]=="sample") &
				(atPOSIX>=as.POSIXct(brush$xmin, origin = "1970-01-01")) &
				(atPOSIX<=as.POSIXct(brush$xmax, origin = "1970-01-01"))
			,"ID"]
			output$info_files_neg_samp <- renderText({
				if(length(these6)>0){paste("<font color=\"black\"> Sample IDs: ",paste(these6,collapse=", ")," </font>",sep="")}else{paste("<font color=\"black\"> No sample files selected </font>",sep="")}
			})
			# negative, blind	
			these7<-measurements[
				(measurements[,"Mode"]=="negative") &
				(measurements[,"Type"]=="blank") &
				(atPOSIX>=as.POSIXct(brush$xmin, origin = "1970-01-01")) &
				(atPOSIX<=as.POSIXct(brush$xmax, origin = "1970-01-01"))
			,"ID"]
			output$info_files_neg_blind<- renderText({
				if(length(these7)>0){paste("<font color=\"green\"> Blanks/blind IDs: ",paste(these7,collapse=", ")," </font>",sep="")}else{paste("<font color=\"green\"> No blind files selected </font>",sep="")}
			})
			# negative, calibration	
			if(any((measurements[,"Type"]=="calibration") & (measurements[,"Mode"]=="negative"))){
				these8<-which((measurements[,"Mode"]=="negative") &(measurements[,"Type"]=="calibration") )
				date_start<-measurements[these8,"Date"]
				date_end<-measurements[these8,"date_end"]
				time_start<-measurements[these8,"Time"]
				time_end<-measurements[these8,"time_end"]
				datetime_start<-c()
				datetime_end<-c()
				for(i in 1:length(time_start)){
					datetime_start<-c(datetime_start,paste(date_start[i],time_start[i],"CET",sep=" "))
					datetime_end<-c(datetime_end,paste(date_end[i],time_end[i],"CET",sep=" "))
				}
				atPOSIX_start<-as.POSIXct(datetime_start);	
				atPOSIX_end<-as.POSIXct(datetime_end);	
				these8<-these8[
					(atPOSIX_start<=as.POSIXct(brush$xmax, origin = "1970-01-01")) &
					(atPOSIX_end>=as.POSIXct(brush$xmin, origin = "1970-01-01"))
				]
				these9<-unique(measurements[these8,"tag2"])
				these8<-measurements[these8,"ID"]
			}else{
				these8<-c()
				these9<-c()
			}
			output$info_files_neg_cal<- renderText({
				if(length(these8)>0){
					paste("<font color=\"red\"> Calibration file IDs: ",paste(these8,collapse=", ")," </font>",sep="")
				}else{
					paste("<font color=\"red\"> No calibration files selected </font>",sep="")
				}
			})
			output$info_files_neg_calgroup<- renderText({
				if(length(these9)>0){
					paste("<font color=\"red\"> Calibration file groups: ",paste(these9,collapse=", ")," </font>",sep="")
				}else{
					paste("<font color=\"red\"> No calibration file groups selected </font>",sep="")
				}
			})		
			# negative, spiked	
			these10<-measurements[
				(measurements[,"Mode"]=="negative") &
				(measurements[,"Type"]=="spiked") &
				(atPOSIX>=as.POSIXct(brush$xmin, origin = "1970-01-01")) &
				(atPOSIX<=as.POSIXct(brush$xmax, origin = "1970-01-01"))
			,"ID"]
			output$info_files_neg_spiked<- renderText({
				if(length(these10)>0){paste("<font color=\"blue\"> Spiked file IDs: ",paste(these10,collapse=", ")," </font>",sep="")}else{paste("<font color=\"blue\"> No spiked files selected </font>",sep="")}
			})			
			#############################################			

			rm(measurements)
		}	
	})
	####################################################
	
}


shinyApp(ui, server)




