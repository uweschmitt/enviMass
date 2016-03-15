
conditionalPanel( 
		condition = "output.textit == 'Waiting...'",
		div(style = widget_style4,	  
			tags$h4("Welcome to enviMass version 3.1"),
			tags$h5("The new workflow for processing high-resolution LC-MS data and temporal trend detection"),
			imageOutput("peakrot", width="10%",height="10%")
		),
        HTML(
			'<p>
			<br/>
			<br/>
			Workflow features:
			<li> Data conversion & centroidization (utilizes ProteoWizard <a href="http://proteowizard.sourceforge.net/" target="_blank">msconvert</a>)
			<li> Data import (.mzXML)
			<li> EIC extraction			
			<li> Peak picking
			<li> Mass recalibration
			<!--	<li> RT Alignment 		--> 
			<li> Quality control
			<li> Profile time series extraction or processing of replicate files 
			<!-- <li> Replicate sample merging --> 
			<li> Intensity normalization: global & profile-based
			<li> Trend detection
			<li> Blind/blank subtraction	
			<li> Profile mass estimation	
			<!--<li> Nontarget peak grouping to components 	--> 
			<!--<li> Isotope pattern calculation --> 
			<!--<li> Target, suspect and internal standard screening: sample-wise & component-wise --> 
			<!--<li> Homologue series detection --> 
			<!--<li> Mass defect analysis --> 
			<li> Handles different ionization modes in one project
			<br/>
			<br/>
			To start a new enviMass project, choose a project name, give a valid project directory and press Start. 
			Existing projects cannot be overwritten.
			<br/>
			To open an existing enviMass project, insert (copypaste) the path to your enviMass project folder (contains
			a logfile.emp) and press Open.
			</p> 
			'
		),
		div(style = widget_style2,	  
			HTML('<p><font color="green">
			Interested in isotope pattern calculation, profile convolution and centroidization? Then check our interactive
			<a href="http://www.envipat.eawag.ch/" target="_blank">enviPat</a> online tool. 
			</font>
			</p>')
		),
		div(style = widget_style2,	  
			HTML('<p><font color="green">
			Have a look at your raw data or run peak picking independently of enviMass? Then check the new 
			<a href="http://cran.r-project.org/web/packages/enviPick/index.html" target="_blank">enviPick R package</a>! 
			</font>
			</p>')
		),
	    HTML('<font color="white">')
)

