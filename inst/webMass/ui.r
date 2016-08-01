widget_style <-
  "display: inline-block;
  vertical-align: text-top;
  padding: 7px;
  border: solid;
  border-width: 2px;
  border-radius: 4px;
  border-color: #ACA;
  background-color:#cccccc;"

widget_style2 <-
  "display: inline-block;
  vertical-align: text-top;
  padding: 7px;
  border: solid;
  border-width: 2px;
  border-radius: 4px;
  border-color: #ACA;
  background-color:#FFFFFF;"

widget_style3 <-
  "display: inline-block;
  vertical-align: text-top;
  padding: 7px;
  border: solid;
  border-width: 2px;
  border-radius: 4px;
  border-color: #FFFFFF;
  background-color:#FFFFFF;"

widget_style4 <-
  "display: inline-block;
  vertical-align: text-top;
  align: center;
  padding: 7px;
  border: solid;
  border-width: 2px;
  border-radius: 4px;
  border-color: #FFFFFF;
  background-color:#FFFFFF;"

widget_style5 <-
  "display: inline-block;
  vertical-align: text-top;
  padding: 7px;
  border: solid;
  border-width: 3px;
  border-radius: 4px;
  border-color: darkgrey;
  background-color: #cccccc;"

widget_style6 <-
  "max-width: 600px;"  
 
shinyUI(
	fluidPage(
	################################################################################
	################################################################################
	
		##############################################################################
		conditionalPanel( 
			condition = "output.textit == 'Waiting...'", 
			titlePanel(
				HTML('
					<p style="background-color: darkgreen">
					<font color="#FFFFFF" size="6">
					&nbsp enviMass v3.1
					</font><br/></p>'),
					windowTitle="enviMass v3.1"
				)
		),
		sidebarLayout(
			##########################################################################
			#sidebarPanel( # included in sourced script - check why!
			source("ui_sidebar.R", local=TRUE)$value,
			#)  
			##########################################################################
			mainPanel(
				tags$head(	
					tags$style(HTML("
						li a{color: black; background-color: darkgrey}; 
						.tabs-above > .nav > li[class=active] > a {
							background-color: #870000;
							color: #FFFFFF;};	
					"))
				),
				useShinyjs(),  # Set up shinyjs
				HTML('</font>'),
				source("ui_busy.R", local=TRUE)$value,  
				HTML('</font>'),
				source("ui_mainPanel_startup.R", local=TRUE)$value,
				HTML('</font>'),	
				source("ui_mainPanel.R", local=TRUE)$value 
			, style = "float:left")
			##########################################################################  
		) 
	################################################################################
	################################################################################
	)
)





