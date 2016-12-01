options(shiny.maxRequestSize=2000*1024^2)
verbose<-TRUE

shinyServer(function(input, output, session){
################################################################################
################################################################################

  cat("\n I run in ");print(environment())
  in_envir<<-environmentName(environment())
  if(any(ls()=="logfile")){stop("\n illegal logfile detected in server.r #1")}  
  ##############################################################################
  # load data ##################################################################  
  if(!any(objects(envir=as.environment(".GlobalEnv"))=="isotopes")){data(isotopes,package="enviPat",envir=as.environment(".GlobalEnv"))}
  if(!any(objects(envir=as.environment(".GlobalEnv"))=="adducts")){data(adducts,package="enviPat",envir=as.environment(".GlobalEnv"))}
  if(!any(objects(envir=as.environment(".GlobalEnv"))=="resolution_list")){data(resolution_list,package="enviPat",envir=as.environment(".GlobalEnv"))}
  if(any(names(resolution_list)=="Elite/R240000@400")){
	shinyjs:::info("library enviPat is not up to date - or you have loaded an old workspace containing old enviPat resolution data lists. Update enviPat and clean your workspace before continuing with enviMass!");
  }	
  ##############################################################################
  # define variables, inputs, outputs - if not in server.startup.R #############
  tried<-try(getVolumes()(),silent=FALSE)
  if(!inherits(tried,"try-error")){
	shinyFileChoose(input, "pro_dir3", session=session, roots=getVolumes(), filetypes=c("emp"), updateFreq = 30000)
	shinyFileSave(input, "download_IS", roots=getVolumes(), updateFreq=30000)
	shinyFileSave(input, "download_target", roots=getVolumes(), updateFreq=30000)	
  }else{
	createAlert(session,anchorId = "alert_4", alertId="a4", title = NULL, content="logfile select disabled, used folder path input",style = "alarm",append=FALSE,dismiss=TRUE)
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
  source("server_startup.r", local=TRUE) 
  ##############################################################################  
  # add compounds, measurements, etc ###########################################
  source("server_obs_Add.r", local=TRUE)
  ##############################################################################  
  # observe request for results ################################################
  source("server_obs_res_meas.r", local=TRUE)
  ##############################################################################  
  # observe parameter changes ##################################################  
  source("server_variables_out.r", local=TRUE)
  ##############################################################################  
  # observe export requests ####################################################  
  source("server_export.r", local=TRUE)
  ############################################################################## 
  # output screening results ###################################################
  source("server_obs_screening.r", local=TRUE)  
  ############################################################################## 
  # observe calibration sets ###################################################
  source("server_obs_calibration.r", local=TRUE)  
  # output network js ##########################################################
  source("server_force.r", local=TRUE)    
  ##############################################################################  
  # run calculations ###########################################################
  source("server_calculation.r", local=TRUE)
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
                                                      