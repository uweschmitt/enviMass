sidebarPanel(
  tags$head(
        tags$style(type="text/css", "select { width: 200px; }"),
        tags$style(type="text/css", "textarea { max-width: 185px; }"),
        tags$style(type="text/css", ".jslider { max-width: 200px; }"),
        tags$style(type='text/css', ".well { max-width: 230px; }"),
        tags$style(type='text/css', ".span4 { max-width: 350px; }")
      ),
    verbatimTextOutput("textit"),
    # start panel ##############################################################
    conditionalPanel(
		condition = "output.textit == 'Waiting...'",
		bsAlert("failed_new"),bsAlert("failed_open"),
		# (1) to start a new project #############################################
		tags$h4("Start new project"),
		textInput("pro_name", "New project name:", value = "new_project_name"),
		textInput("pro_dir", "Project directory:", value = "D:\\Users\\"),	  
		bsPopover("pro_dir", 
			title = "Path to folder that will contain the new project.",
			content = "The new project folder will be generated automatically, using the project name.", 
			placement = "right", trigger = "hover"),
		bsActionButton("newit","Start",style="success"),
		# (2) to open an existing project ########################################
		tags$h4("Open existing project"),
		helpText("Type path to project folder ..."),
		textInput("pro_dir2", "", value = "C:\\...\\old_project_name"),
		bsPopover("pro_dir2", 
			title = "Insert full path, including the project folder, but excluding the logfile.emp.",
			content = "Using your OS explorer, you may navigate into your project folder and copy/paste the full path.", 
			placement = "right", trigger = "hover"),
		helpText("... or browse for its logfile.emp:"),
		bsAlert("alert_4"),
		shinyFilesButton(id="pro_dir3", label="logfile select", 
			title="Please select a logfile.emp in the project folder", FALSE),
		bsPopover("pro_dir3", 
			title = "Select logfile.emp",
			content = "... located inside the project folder of the desired project.", 
			placement = "right", trigger = "hover"),	
		helpText(""),				  
		bsActionButton("openit","Open",style="success"),
		textOutput("had_opened")
    ),
    # action panel #############################################################
    conditionalPanel(
      condition = "output.textit != 'Waiting...'",
	  bsAlert("reset"),
      helpText("Current state:"),
      verbatimTextOutput("dowhat"),
      helpText("Finished tasks:"),
      tableOutput("summar"),
	  HTML('<hr noshade="noshade" />'),
      bsActionButton("Calc","Calculate",style="danger"),
	  bsPopover("Calc", 
		title = "Start new project (re)calculation.",
		content = "The current settings for parameters and workflow steps will be used. Calculation results will be displayed in the results tabset.", 
		placement = "right", trigger = "hover"),
	  textOutput("had_calculated"),
	  HTML('<hr noshade="noshade" />'),
      actionButton("Restart","Back"),
	  bsPopover("Restart", 
		title = "Return to start page ...",content = "... to start a new project or open an existing one.", 
		placement = "right", trigger = "hover"),
      actionButton("Exit","Exit")
    ),
	HTML('<font color="white">')
)
