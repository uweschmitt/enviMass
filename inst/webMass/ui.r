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
  
shinyUI(
pageWithSidebar(
################################################################################
################################################################################

  ##############################################################################
  headerPanel(
	HTML('
	<p style="background-color:grey">
	<font color="#FFFFFF" size="5">
	&nbsp enviMass v2.0
	</font><br/></p>'),
	windowTitle="enviMass"
  ),
  ##############################################################################
  source("ui_sidebar.R", local=TRUE),   
  ##############################################################################
  mainPanel(
    HTML('</font>'),
	source("ui_busy.R", local=TRUE),  
    HTML('</font>'),
	source("ui_mainPanel_startup.R", local=TRUE),
	HTML('</font>'),	
	source("ui_mainPanel.R", local=TRUE) 
  )
  ##############################################################################
   
################################################################################
################################################################################
)
)





