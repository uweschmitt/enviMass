
conditionalPanel( 
		condition = "output.textit == 'Waiting...'",
		div(style = widget_style4,	  
			tags$h4("Welcome to enviMass version 3.128"),
			tags$h5("The new workflow for processing high-resolution LC-MS data and temporal trend detection"),
			imageOutput("peakrot", width="10%",height="10%")
		),
        HTML(
			'<p>
			<br/>
			<br/>
			Workflow features:
			<li> Data import (.mzXML, centroided data)
			<li> Direct data conversion & centroidization for Thermo .raw files (utilizes ProteoWizard <a href="http://proteowizard.sourceforge.net/" target="_blank">msconvert</a>)
			<li> EIC extraction			
			<li> Peak picking
			<li> Data visualization
			<li> Mass recalibration
			<!--	<li> RT Alignment 		--> 
			<li> Quality control
			<li> Profile time series extraction or processing of replicate files 
			<li> Replicate sample intersection
			<li> Intensity normalization: median- & ISTD-based
			<li> Trend detection
			<li> Blind/blank subtraction	
			<li> Profile mass estimation	
			<li> Nontarget peak grouping to components 
			<li> Isotope pattern calculation  
			<li> Target, suspect and internal standard screening
			<li> Homologue series detection
			<li> Calibration, quantification and recovery steps
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
			Interested in isotope pattern calculation, profile convolution and centroidization? Then check the interactive
			<a href="http://www.envipat.eawag.ch/" target="_blank">enviPat online tool</a>. 
			</font>
			</p>')
		),
		div(style = widget_style2,	  
			HTML('<p><font color="green">
			Run a homologue series detection separately? Then check the new 
			<a href="http://www.envihomolog.eawag.ch/" target="_blank">online tool</a>! 
			</font>
			</p>')
		),
	    HTML('<font color="white">')
)

