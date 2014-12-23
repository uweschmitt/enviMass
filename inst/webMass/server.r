options(shiny.maxRequestSize=2000*1024^2)

shinyServer(function(input, output, session){
################################################################################
################################################################################

  ##############################################################################
  # load data ##################################################################  
  if(!any(objects(envir=as.environment(".GlobalEnv"))=="isotopes")){data(isotopes,package="enviPat",envir=as.environment(".GlobalEnv"))}
  if(!any(objects(envir=as.environment(".GlobalEnv"))=="adducts")){data(adducts,package="enviPat",envir=as.environment(".GlobalEnv"))}
  if(!any(objects(envir=as.environment(".GlobalEnv"))=="resolution_list")){data(resolution_list,package="enviPat",envir=as.environment(".GlobalEnv"))}
  ##############################################################################
  # define variables, inputs, outputs - if not in server.startup.R #############
  tried<-try(getVolumes()(),silent=FALSE)
  if(!inherits(tried,"try-error")){
	shinyFileChoose(input,"pro_dir3", session=session, roots=getVolumes(), filetypes=c("emp") )
  }else{
	createAlert(session,inputId = "alert_4", alertId="a4", title = NULL, message="logfile select disabled, used folder path input",type = "alarm",append=FALSE,block=TRUE,dismiss=TRUE)
  }
  output$textit<-renderText("Waiting...")
  output$dowhat<-renderText("Open")
  output$isotable<-renderTable(isotopes)
  updateCheckboxGroupInput(session, "adducts_pos", "Positive ions:", choices =  as.character(adducts[adducts[,6]=="positive",1]),selected=as.character(adducts[adducts[,6]=="positive",1][1]))
  updateCheckboxGroupInput(session, "adducts_neg", "Negative ions:", choices =  as.character(adducts[adducts[,6]=="negative",1]),selected=as.character(adducts[adducts[,6]=="negative",1][1]))               
  updateSelectInput(session, "resolution", "Instrument resolution:", choices =  names(resolution_list), selected= (names(resolution_list)[1]))                   
  init<-reactiveValues() 	# reactive value to indicate ...
  init$a<-"FALSE"  			# ... if/when a project is opened
  init$b<-1					# ... increments each time a calculation is run or ion mode selected; trigger observers to refresh results
  logobusy<-list(src="circ.gif");
  output$logobusy1<-renderImage(logobusy, deleteFile = FALSE);
  output$logobusy2<-renderImage(logobusy, deleteFile = FALSE);
  output$logobusy3<-renderImage(logobusy, deleteFile = FALSE);  
  peakrot<-list(src="peakrot.gif");
  output$peakrot<-renderImage(peakrot, deleteFile = FALSE); 
  ##############################################################################
  # start projects or load them - get their individual logfile #################
  # also load parameter settings ###############################################
  source("server_startup.R", local=TRUE)
  ##############################################################################  
  # add compounds, measurements, etc ###########################################
  source("server_obs_Add.R", local=TRUE)
  ##############################################################################  
  # observe request for results ################################################
  source("server_obs_res_meas.R", local=TRUE)
  ##############################################################################  
  # observe parameter changes ##################################################  
  source("server_variables_out.R", local=TRUE)
  ##############################################################################  
  # run calculations ###########################################################
  source("server_calculation.R", local=TRUE)
  observe({ # Restart enviMass
    input$Restart
    if(input$Restart){
      output$textit<<-renderText("Waiting...")
	  source("server_cleaner.R", local=TRUE);		
      cat("Restart\n")
    }
  })
  observe({ # Exit enviMass
    input$Exit
    if(input$Exit){
		source("server_cleaner.R", local=TRUE);		
		stopApp(returnValue="Quit enviMass browser session")
    }
  })
  ##############################################################################

################################################################################
################################################################################
})
                                                      