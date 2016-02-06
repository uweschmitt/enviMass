	  
    conditionalPanel( 
		condition = "output.textit != 'Waiting...'",  
		tags$head(
			tags$style(type = "text/css", "li a{color: black; background-color: darkgrey}") 
		),		
		tags$style(HTML("
			.tabs-above > .nav > li[class=active] > a {
			background-color: #870000;
			color: #FFFFFF;}")
		),
		tags$h5(""),
		bsAlert("alert_3"),
		tabsetPanel(
		#navbarPage("",
        ########################################################################
        # MEASUREMENTS #########################################################
        ########################################################################
		tabPanel("Files",
   		  div(style = widget_style,
            tags$h4("Add LC-HRMS files"),
            helpText("To add a new file, set the below specifications accordingly and press Add"),
            textInput("Measadd_ID", "Numeric ID:", value = "123"),
			radioButtons("Measadd_ID_autom", "...or assign ID automatically?", c("yes"="yes","no"="no")),
            textInput("Measadd_name", "Name:", value = "Sample 1"),
            selectInput("Measadd_type", "Type:", choices = c("sample", "blank", "doted", "other")),
            selectInput("Measadd_incl", "Include?", choices = c("TRUE","FALSE")),
            selectInput("Measadd_mode", "Choose ionization mode:", choices = c("positive", "negative")),
            textInput("Measadd_place", "Place:", value = "Rhine"),
			dateInput("Measadd_date", "Date", value = NULL, min = NULL,max = NULL, format = "yyyy-mm-dd", startview = "month",weekstart = 0, language = "en"),
            textInput("Measadd_time", "Time:(HH:MM:SS)", value = "12:00:00"),
            fileInput("Measadd_path", "Select centroided .mzXML file:", multiple = FALSE, accept = c(".mzXML",".raw")),
			bsPopover("Measadd_path", 
				title = "WARNING",
				content = "Files must be centroided. Check Package enviPick whether peaks can be picked properly from your files. Tested with Orbitrap files only.
				If reload fails, press cancel in the file upload window first", 
				placement = "right", trigger = "hover"),
			textOutput("had_meas_added")
		  ),
          div(style = widget_style,
            tags$h4("Delete file"),
            textInput("Measdel_ID", "ID:", value = "123"),
            actionButton("Measdel","Remove")
          ),
		  div(style = widget_style,
            tags$h4("Import project (files only)"),
			textInput("import_pro_dir", "", value = "C:\\...\\old_project_name"),
			bsPopover("import_pro_dir", 
				title = "Insert full path, including the project folder, but excluding the logfile.emp.",
				content = "Using your OS explorer, you may navigate into your project folder and copy/paste the full path.", 
				placement = "right", trigger = "hover"),
			checkboxInput("Merge_project", "Omit duplicates?", FALSE),
			bsPopover("Merge_project", 
				title = "File duplicate handling",
				content = "A file with the same type, time, date, ionization and place as one which already exists will not be imported.", 
				placement = "right", trigger = "hover"),
            actionButton("Import_project","Import")
          ),
          helpText(""),
          DT::dataTableOutput("measurements")
        ),
        ########################################################################
        # Compounds ############################################################
        ########################################################################
        tabPanel("Compounds",
		 tabsetPanel(
            tabPanel("Internal standards",
                  div(style = widget_style,
                    tags$h4("Add internal standard compound"),
                    helpText("To add a new internal standard, fill out the below form and press Add"),
                    textInput("ISadd_ID", "Unique ID:", value = "123_XYZ"),          
                    textInput("ISadd_name", "Name:", value = "CompoundX"),
                    textInput("ISadd_formula", "Formula:", value = "C10H12O10"),          
                    textInput("ISadd_RT", "Retention time (RT) [min]:", value = "7.52"),    
					checkboxInput("ISadd_RTtol_use", "Use specific +/- RT tolerance [min]:", FALSE),                    
                    textInput("ISadd_RTtol","", value = "2"),          
                    selectInput(inputId="ISadd_charge", label="Ionization mode", choices=c("positive","negative"), selected = "positive", multiple = FALSE),
					popify(
						selectInput(inputId="ISadd_add", label="Main adduct:", choices= "FALSE", selected = "FALSE", multiple = FALSE),
							title = "Adduct definition for a compound",
							content = 	"A compound-specific adduct can be defined here; non-specific adducts are defined in the Settings tab. Unless the below restriction is enabled, the compound-specific adduct is used alongside the non-specific ones. Seperate compound entries have to be made when including more than one compound-specific adduct.", 
							placement = "right", trigger = "hover"),		
					checkboxInput("ISadd_rest_adduct", "Restrict screening to main adduct?", FALSE),       
					checkboxInput("ISadd_use_recal", "To be used for m/z recalibration?", TRUE),
					checkboxInput("ISadd_use_screen", "To be used for screening?", TRUE),
					textInput("ISadd_remark", "remark", value = "none"),
                    textInput("ISadd_tag1", "tag1", value = "none"),
                    textInput("ISadd_tag2", "tag2", value = "none"),
                    textInput("ISadd_tag3", "tag3", value = "none"),			
					checkboxInput("ISadd_date", "Restrict to temporal range?", FALSE),	
					dateRangeInput("ISadd_date_range", label="", start = NULL, end = NULL,
						min = NULL, max = NULL, format = "yyyy-mm-dd",
						startview = "month", weekstart = 0, language = "en",
						separator = " to "),
					actionButton("AddIS","Add")
                  ),
                  div(style = widget_style,
                    tags$h4("Remove internal standard compound"),
                    helpText("To delete a compound from the list, type in its ID and press Delete"),
                    textInput("ISdelete_ID", "ID for deletion:", value = "123_XYZ"),          
                    actionButton("DeleteIS","Delete")
                  ),
				  div(style = widget_style2,
                    tags$h4("Import compound list"),
					fileInput("ISlist_path", "Select IS.txt file, e.g. from dataframes folder of another project", multiple = FALSE, accept = c(".txt"))
                  ),
                helpText(""),
                DT::dataTableOutput("IS")
            ),
            tabPanel("Targets", 
                  div(style = widget_style,
                    tags$h4("Add target or suspect compound"),
                    helpText("To add a new target, fill out the below form and press Add."), #For suspects, define RT tolerances <br> spannning the full elution time range and set its RT to the middle time point."),
                    textInput("targetsadd_ID", "Unique ID:", value = "123_XYZ"),          
                    textInput("targetsadd_name", "Name:", value = "CompoundY"),
                    textInput("targetsadd_formula", "Formula:", value = "C10H12O10"),          
                    textInput("targetsadd_RT", "Retention time (RT) [min]:", value = "7.52"),    
					checkboxInput("targetsadd_RTtol_use", "Use specific +/- RT tolerance [min]:", FALSE),                    
                    textInput("targetsadd_RTtol","", value = "2"),          
                    selectInput(inputId="targetsadd_charge", label="Ionization mode", choices=c("positive","negative"), selected = "positive", multiple = FALSE),
					popify(
						selectInput(inputId="targetsadd_add", label="Main adduct:", choices= "FALSE", selected = "FALSE", multiple = FALSE),
							title = "Adduct definition for a compound",
							content = 	"A compound-specific adduct can be defined here; non-specific adducts are defined in the Settings tab. Unless the below restriction is enabled, the compound-specific adduct is used alongside the non-specific ones. Seperate compound entries have to be made when including more than one compound-specific adduct.", 
							placement = "right", trigger = "hover"),		
					checkboxInput("targetsadd_rest_adduct", "Restrict screening to main adduct?", FALSE),                    
					checkboxInput("targetsadd_use_recal", "To be used for m/z recalibration?", TRUE),
					checkboxInput("targetsadd_use_screen", "To be used for screening?", TRUE),
					textInput("targetsadd_remark", "remark", value = "none"),
                    textInput("targetsadd_tag1", "tag1", value = "none"),
                    textInput("targetsadd_tag2", "tag2", value = "none"),
                    textInput("targetsadd_tag3", "tag3", value = "none"),			
					checkboxInput("targetsadd_date", "Restrict to temporal range?", FALSE),	
					dateRangeInput("targetsadd_date_range", label="", start = NULL, end = NULL,
						min = NULL, max = NULL, format = "yyyy-mm-dd",
						startview = "month", weekstart = 0, language = "en",
						separator = " to "),
                    actionButton("Addtargets","Add")
                  ),
                  div(style = widget_style,
                    tags$h4("Remove target compound"),
                    helpText("To delete a compound from the list, type in its ID and press Delete"),
                    textInput("targetsdelete_ID", "ID for deletion:", value = "123_XYZ"),          
                    actionButton("Deletetargets","Delete")
                  ),
                  div(style = widget_style2,
                    tags$h4("Import compound list"),
					fileInput("targetlist_path", "Select targets.txt file, e.g. from dataframes folder of another project", multiple = FALSE, accept = c(".txt"))
                  ),
                helpText(""),
                DT::dataTableOutput("targets")
            )
          )  
        ),
        ########################################################################
        # WORKFLOW OPTIONS #####################################################
        ########################################################################
        tabPanel("Workflow options",
			tags$h5("Apply settings to project?"), 
			bsButton("saveflow","Apply",style="warning"),
			bsAlert("alert_1"),
				# block 1 ######################################################
				HTML('<hr noshade="noshade" />'),
				HTML('<p style="background-color:darkgrey"; align="center"> <font color="#FFFFFF"> File upload </font></p> '),
				HTML('<p style="background-color:darkgrey"; align="center"> <font color="#FFFFFF"> Peak picking </font></p> '),
				HTML('<p style="background-color:darkgrey"; align="center"> <font color="#FFFFFF"> Quality control filter </font></p> '),
					radioButtons("qc", "Include?", c("yes"="yes","no"="no")),          
				HTML('<hr noshade="noshade" />'),
				HTML('<h1 align="center"> &#x21e9; </h1> '),
				# block 2 ######################################################
				HTML('<p style="background-color:darkblue"; align="center"> <font color="#FFFFFF"> Mass recalibration </font></p> '),
					radioButtons("recal", "Include?", c("yes"="yes","no"="no")),                
				HTML('<p style="background-color:darkblue"; align="center"> <font color="#FFFFFF"> Median intensity normalization </font></p> '),
					radioButtons("intnorm", "Include?", c("yes"="yes","no"="no")),
				#HTML('<p style="background-color:darkgreen"; align="center"> <font color="#FFFFFF"> RT alignment </font></p> '),
				#radioButtons("RTalign", "Include?", c("yes"="yes","no"="no")),     
				HTML('<p style="background-color:darkblue"; align="center"> <font color="#FFFFFF"> Replicate filter </font></p> '),				
					radioButtons("replicates", "Include? ", c("yes"="yes","no"="no")),
				HTML('<hr noshade="noshade" />'),
				HTML('<h1 align="center"> &#x21e9; </h1> '),                     
				# block 3 ######################################################
				HTML('<p style="background-color:darkgreen"; align="center"> <font color="#FFFFFF"> Profile extraction </font></p>'),
					radioButtons("profiled", "Include? ", c("yes"="yes","no"="no")),
				HTML('<p style="background-color:darkgreen"; align="center"> <font color="#FFFFFF"> LOD interpolation </font></p>'),
					radioButtons("LOD_interpol", "Include? ", c("yes"="yes","no"="no")),					
				HTML('<p style="background-color:darkgreen"; align="center"> <font color="#FFFFFF"> Compound screening </font></p> '),
					checkboxInput("screen_IS", "Screen internal standards?", TRUE),
					checkboxInput("screen_target", "Screen targets/suspects?", TRUE),	
				HTML('<hr noshade="noshade" />'),
				HTML('<h1 align="center"> &#x21e9; </h1> '),  					
				# block 4 ######################################################
				HTML('<p style="background-color:darkred"; align="center"> <font color="#FFFFFF"> Normalization using IS-profiles </font></p> '),
					radioButtons("profnorm", "Include? ", c("yes"="yes","no"="no")),
				HTML('<p style="background-color:darkred"; align="center"> <font color="#FFFFFF"> Quantification </font></p> '),
					radioButtons("quantif", "Include? ", c("yes"="yes","no"="no")),
				HTML('<p style="background-color:darkred"; align="center"> <font color="#FFFFFF"> Trend detection & blind subtraction </font></p> '),
					radioButtons("trenddetect", "Include? ", c("yes"="yes","no"="no")),
				HTML('<hr noshade="noshade" />')                  
				#HTML('<p style="background-color:darkred"; align="center"> <font color="#FFFFFF"> Peak grouping (componentization) </font></p>'),
				#	checkboxInput("Comp_isotop", "Group isotopologue peaks?", TRUE),
				#	checkboxInput("Comp_add", "Group adduct peaks?", TRUE),
				#HTML('<p style="background-color:darkred"; align="center"> <font color="#FFFFFF"> Component-wise screening </font></p> '),
				#	checkboxInput("screen_IS_comp", "Screen internal standards?", TRUE),
				#	checkboxInput("screen_target_comp", "Screen internal targets/sustpects?", TRUE),
				#HTML('<hr noshade="noshade" />'),
				#HTML('<h1 align="center"> &#x21e9; </h1> '),                     
				# block 5 ######################################################
				#HTML('<p style="background-color:orange"; align="center"> <font color="#FFFFFF"> Homologue series detection </font></p> '),
					#radioButtons("homol", "Include? ", c("yes"="yes","no"="no")),	
				#HTML('<p style="background-color:orange"; align="center"> <font color="#FFFFFF"> Mass defect calculation </font></p> '),
					#radioButtons("massdef", "Include? ", c("yes"="yes","no"="no")),	
				#HTML('<hr noshade="noshade" />')
				################################################################
        ),
        ########################################################################
        # PARAMETER SETTINGS ###################################################
        ########################################################################
        tabPanel("Settings",     
		  tags$h5("Apply settings to project?"), 
		  bsButton("savepar","Apply",style="warning"),
		  bsAlert("alert_2"),
          HTML('<hr noshade="noshade" />'),
          tabsetPanel(
            # PEAK PICKING #####################################################
            tabPanel("Peak picking",
				div(style = widget_style2,
					tags$h5("Peak definition"), 
					numericInput("peak_minpeak", "Minimum number of measurements per peak ...", 4),
					sliderInput("peak_drtsmall2", "... within a given RT window [s]", min = 1, max = 200, value = 20, step= 0.1),
					sliderInput("peak_drtfill", "Maximum RT gap length to be interpolated [s]", min = 0, max = 60, value = 10, step= 0.1),
					sliderInput("peak_drtdens2", "Peak definition - Maximum RT length of a single peak", min = 10, max = 1500, value = 120, step= 0.1),
					sliderInput("peak_minint", "Minimum log10(intensity) threshold", min = 0, max = 10, value = 4, step= .1),
					numericInput("peak_SN", "Minimum Signal/Noise", 5),
					numericInput("peak_SB", "Minimum Signal/Base", 2),
					numericInput("peak_recurs", "Maximum possible number of peaks within a single EIC", 3)
				),
				div(style = widget_style2,
					tags$h5("EIC partitioning & clustering"),
					sliderInput("peak_drtgap", "Maximum retention time gap in an EIC", min = 20, max = 1500, value = 300, step= 1),
					sliderInput("peak_dmzdens", "Maximum m/z deviation of a measurement from its EIC mean [ppm]", min = 1, max = 20, value = 3.5, step= 0.1)       				
				),
				div(style = widget_style,
					tags$h5("Advanced"),
					numericInput("peak_ended", "How often can a peak detection fail to end the recursion? - peak picking", 1),
					numericInput("peak_weight", "Weight for assigning measurements to a peak - peak picking", 1),
					sliderInput("peak_maxint", "Upper log10(intensity) safety threshold", min = 0, max = 15, value = 6.7, step= .1)
				),
              tags$h4(""),
              tags$h4("")              
            ),
            # ADDUCTS ##########################################################
            tabPanel("Adducts",
                  helpText("Select adducts for the calculation of centroid masses of IS & targets isotope patterns.\n Used for recalibration and screening."),
                  div(style = widget_style3,checkboxGroupInput("adducts_pos", "Positive ions:", "none")),
                  div(style = widget_style3,checkboxGroupInput("adducts_neg", "Negative ions:", "none"))
            ),
			# RESOLUTION #######################################################
            tabPanel("Resolution",
                  div(style = widget_style3,selectInput("resolution", "Instrument resolution:", "none",width='600px')),
				  imageOutput("plot_resolution", height="auto")	  
            ),	
            # RECALIBRATION ####################################################
            tabPanel("Recalibration",
              tags$h5("Mass recalibration"),
              div(
                style = widget_style3,
                selectInput("recal_what", "Reference compounds:", c("Internal standards","Target compounds","both"),"Internal standards",multiple=FALSE),                
                numericInput("recal_dmz", "m/z tolerance ...", 3),                
                selectInput("recal_ppm", "... given in:", choices = c("ppm"="TRUE","absolute"="FALSE"), "TRUE"),
                numericInput("recal_drt", "RT tolerance [s]", 30)   
              )
            ),
			# REPLICATES #######################################################
            tabPanel("Replicates",
				tags$h5("Replicate files are defined (i.e., grouped) by the tag3 entry (not FALSE) in the measurements table"),
				numericInput("replicate_dmz", "m/z tolerance ...", 3),                
				selectInput("replicate_ppm", "... given in:", choices = c("ppm"="TRUE","absolute"="FALSE"), "TRUE"),	
				#selectInput("replicate_recalib", "... and corrected by recalibration results (if available)", choices = c("TRUE"="TRUE","FALSE"="FALSE"), "FALSE"),	
				numericInput("replicate_delRT", "RT tolerance of a compound peaks across replicate samples [s]", 30)
			),	
            # ALLIGNMENT #######################################################
            #tabPanel("Allignment",
            #  tags$h5("RT allignment")
            #),
            # SCREENING ########################################################
            tabPanel("Screening",
              tabsetPanel(
                tabPanel("IS",
					div(style = widget_style2,
						tags$h5("Retention time"),
						numericInput("screen_IS_delRT", "RT tolerance of peaks in sample relative to their expected RT [s]", 30),
						numericInput("screen_IS_dRTwithin", "RT tolerance of peaks within an isotope pattern [s]", 50),
						numericInput("screen_IS_dRTblank", "RT tolerance of peaks in blank/blind relative to their expected RT [s]", 200)
					),
				  	div(style = widget_style2,
						tags$h5("m/z"),
						numericInput("screen_IS_dmz", "m/z tolerance ...", 3),                
						selectInput("screen_IS_ppm", "... given in:", choices = c("ppm"="TRUE","absolute"="FALSE"), "TRUE")
 					),
					div(style = widget_style2,
						tags$h5("Intensity"),
						sliderInput("screen_IS_dInt", "Intensity tolerance %", min = 0, max = 100, value = 30, step= .2),
						numericInput("screen_IS_Intcut", "Lower intensity threhold", 5E4)                
					),
					div(style = widget_style2,
						tags$h5("Scoring"),
						numericInput("screen_IS_w1", "Score weight for mass matching", 0.8),                
						numericInput("screen_IS_w2", "Score weight for relative intensity matching", 0.2),                
						numericInput("screen_IS_w3", "Score weight for occurrence in blank/blind", 0.0)                
					)
                ),
                tabPanel("Targets & Suspects",
					div(style = widget_style2,
						tags$h5("Retention time"),
						numericInput("screen_target_delRT", "RT tolerance of peaks in sample relative to their expected RT [s]", 30),
						numericInput("screen_target_dRTwithin", "RT tolerance of peaks within an isotope pattern [s]", 50),
						numericInput("screen_target_dRTblank", "RT tolerance of peaks in blank/blind relative to their expected RT [s]", 200)
					),
				  	div(style = widget_style2,
						tags$h5("m/z"),
						numericInput("screen_target_dmz", "m/z tolerance ...", 3),                
						selectInput("screen_target_ppm", "... given in:", choices = c("ppm"="TRUE","absolute"="FALSE"), "TRUE")
 					),
					div(style = widget_style2,
						tags$h5("Intensity"),
						sliderInput("screen_target_dInt", "Intensity tolerance %", min = 0, max = 100, value = 30, step= .2),
						numericInput("screen_target_Intcut", "Lower intensity threhold", 5E4)                
					),
					div(style = widget_style,
						tags$h5("Scoring"),
						numericInput("screen_target_w1", "Score weight for mass matching", 0.8),                
						numericInput("screen_target_w2", "Score weight for relative intensity matching", 0.2),                
						numericInput("screen_target_w3", "Score weight for occurrence in blank/blind", 0.0)                
					)
                )
              )
            ),
            # INTENSITY NORMALIZATION ##########################################
            tabPanel("Normalization",
				tags$h5("Intensity normalization based on IS-profiles"),
				sliderInput("profnorm_cover_files", "Minimum of files covered by each IS profile (%)", min = 0, max = 100, value = 90, step= 1),
				sliderInput("profnorm_threshold", "Screening threshold", min = 0, max = 1, value = 0.8, step= .01),
				numericInput("profnorm_cover_isccount", "Minimum number of IS profiles", 15),
				HTML('<hr noshade="noshade" />'),
				checkboxInput("profnorm_use_blank", "Show median deviation of blank/blind profiles?", FALSE),
				checkboxInput("profnorm_use_blank_sample", "Use subsampling", FALSE),
				numericInput("profnorm_use_blank_samplecount", "Number of blank/blind profiles in subsample", 100),
				HTML('<hr noshade="noshade" />'),
				checkboxInput("profnorm_use_nonblank", "Show median deviation of sample (i.e., non-blank) profiles?", FALSE),	
				checkboxInput("profnorm_use_nonblank_sample", "Use subsampling", FALSE),
				numericInput("profnorm_use_nonblank_samplecount", "Number of sample profiles in subsample", 100)
			),
            # PROFILING ########################################################
            tabPanel("Profiling",
				tags$h5("Profile extraction"),
				sliderInput("prof_sets", "Maximum number of newest samples to be processed", min = 50, max = 3000, value = 100, step= 1),
				numericInput("prof_dmz", "Peak deviation within profiles: m/z tolerance ...", 3),                
                selectInput("prof_ppm", "... given in:", choices = c("ppm"="TRUE","absolute"="FALSE"), "TRUE"),
                numericInput("prof_drt", "Peak deviation within profiles: RT tolerance [s]", 60)            ),
            # TREND ############################################################
            tabPanel("Trends",
				tags$h5("Trend detection:"),
				textInput("trend_lags", "Time lags of trends [days], comma-separated:", value = "4,7,14"),
				numericInput("trend_thres", "Trend vs. mean+variance intensity threshold:", 3),
				selectInput("notrend", "Do not show global trend - instead, report it as maximum intensity above blind", choices = c("TRUE"="TRUE","FALSE"="FALSE"), "FALSE")
            ),
            # BLIND #############################################################
            tabPanel("Blind",
				tags$h5("Blind subtraction:"),
				selectInput("blind_do", "Run a blind subtraction...", choices = c("yes"="yes","no"="no"), "yes"),
				numericInput("blind_fold", "...if intensity ratio sample/blind <", 100)           
            ),			
            # GENERAL SETTINGS #################################################
            tabPanel("General",
				div(style = widget_style3,
					textInput("PWpath", "Path to Proteowizard MSConvert (use / and include .exe)", value = "C:/Program Files/ProteoWizard/ProteoWizard 3.0.5140/msconvert.exe")
				),
				div(style = widget_style,
					tags$h5("Debug tools"),
					selectInput("progressbar", "Show progress bars (Windows OS only)", choices = c("TRUE","FALSE"), selected="FALSE"),
					selectInput("do_project_check", "Skip the project check before calculations?", choices = c("TRUE","FALSE"), selected="FALSE"),					
					textInput("upto_file", "Up to file with ID:", value = "FALSE"),
					tags$h6("Reset project without peak picking:"),
					bsButton("reset_1","Reset",style="danger"),
					tags$h6("Reset project with peak picking:"),
					bsButton("reset_2","Reset",style="danger")
				)
			)			
          )
        ),
        ########################################################################
        # RESULTS ##############################################################
        ########################################################################
        tabPanel("Results", 	
			tabsetPanel( 
				tabPanel("Profiles",
					div(style = widget_style5,
						textOutput("had_ion"),	
						selectInput("Ion_mode", label=NULL, c("positive","negative"), selected = ("positive"), multiple = FALSE)
					),
					HTML('<hr noshade="noshade" />'),  
					#navbarPage("", 
					tabsetPanel( 
						tabPanel("Summary",										
								tags$h5("Filter profile list:"),
								div(style = widget_style3,numericInput("filterProf_minmass", "Minimum m/z:", 0)),
								div(style = widget_style3,numericInput("filterProf_maxmass", "Maximum m/z:", 3000)),
								div(style = widget_style3,numericInput("filterProf_minrt", "Minimum RT [s]:", 0)),
								div(style = widget_style3,numericInput("filterProf_maxrt", "Maximum RT [s]:", 100000)),			
								div(style = widget_style3,radioButtons("filterProf_meanblind", "Use mean above blind?", c("no"="no","yes"="yes"))),
								bsPopover("filterProf_meanblind", 
									title = "Replicates, not time series ...",
									content = "Get profiles with mean sample intensity x times above mean blank intensities (set Sort profile list by: maximum or mean intensity; x = to be set in blind settings panel). Useful if your files are replicates and not a time sequences.", 
									placement = "top", trigger = "hover"),
								div(style = widget_style3,radioButtons("filterProf_notblind", "Not in blind?", c("no"="no","yes"="yes"))),
								div(style = widget_style3,selectInput("filterProf_sort", "Sort profile list by:", 
									choices = c("ID","mean m/z","mean RT","maximum intensity","mean intensity","global trend intensity","current trend intensity"), selected="current trend intensity")),
								div(style = widget_style3,numericInput("filterProf_count", "Restrict list size:", 500)),
								HTML('<hr noshade="noshade" />'),
								bsCollapse(multiple = FALSE, open = "col2", id = "collapse2",
									bsCollapsePanel("Profile statistics", 
										div(style = widget_style3,tags$h6("Number of peaks:"), textOutput("atprof1")),
										div(style = widget_style3,tags$h6("Number of profiles:"), textOutput("atprof2")),
										div(style = widget_style3,tags$h6("...containing blind peaks:"), textOutput("atprof3")),
										div(style = widget_style3,tags$h6("...with a past trend:"), textOutput("atprof4")),
										div(style = widget_style3,tags$h6("...with a current trend:"), textOutput("atprof5")),
										value="test4"),							
									bsCollapsePanel("Intensity histogram", 
										imageOutput("profilehisto", height="auto"),
										value="test5"),
									bsCollapsePanel("Profile list", 		
										tableOutput("allproftable"),
										value="test5")
								),
								div(style = widget_style3,
									bsButton("expo_profiles","Export filtered profile list",style="info"),
									textOutput("expo1"),
									bsPopover("expo_profiles", 
										title = "Export above filtered profiles",
										content = "Time-sorted peak intensities of profiles and their mean mass & RT are exported as profiles.txt to the export folder of this project. WARNING: restrict list size to avoid lengthy exports!", 
										placement = "right", trigger = "hover"))				
										
						),
						tabPanel("Newest trends",
								tags$h5("Comparison of current vs. global trends by profile ID"),
								imageOutput("boxprofile", height="auto"),
								HTML(
									'<p>
									</br></br>
									The above boxplot (grey) shows the intensity distributions of all trends of concern, listing the IDs, mean masses (m/z) and mean retention time (RT)
									of the profiles with the most intense trends on the right.
									Colored points are used to elucidate the current trend intensities from the latest input file. Note that a file is only included if 
									surviving the quality check (QC). The red dots signify profiles with intensities in the outlier range of the global (past and latest) trend intensities; 
									green dots symbolize those below. 
									</br></br>
									For more information on the individual profiles, navigate to the Single Profile tab and use the profile ID to extract further profile information.
									</p>'
								)
						),
						tabPanel("Single Profile",
								tags$h5("Extraction of individual profiles"),					
								HTML('<p> Enter the ID of a profile to extract relevant information. Profile IDs are listed both in the Summary tab and the Newest trends tab. 
								Alternatively, sort and filter the profile list in the Summary tab and choose an entry number to show a listed profile. </p>'),							
								div(style = widget_style3,numericInput("profID", "profile ID:", 0)),
								div(style = widget_style3,numericInput("profentry", "Entry # in (filtered, sorted) profile list:", 0)),
								div(style = widget_style3,radioButtons("prof_log", "Logarithmic intensity?", c("no"="no","yes"="yes"))),
								div(style = widget_style3,textOutput("peak_number")),
								imageOutput("timeprofile", height="auto"),
								bsCollapse(multiple = FALSE, open = "col1", id = "collapse1",
									bsCollapsePanel("Profile EICs & Peak viewer", 
										div(style = widget_style3,numericInput("profpeakID", "Peak entry #:", min=0, 0)),
										bsPopover("profpeakID", title = "View extracted chromatograms & peaks of the selected profile (sample files only).",
											content = "Select peaks in the order listed in the profile peak table, i.e. from latest to oldest file.", placement = "right", trigger = "hover"),
										div(style = widget_style3,textOutput("prof_peak_text")),
										HTML('<hr noshade="noshade" />'),
										imageOutput("profile_position", height="auto"),
										imageOutput("profile_EIC", height="auto"),
										value="test1"),
									bsCollapsePanel("Profile mass estimation",
										div(style = widget_style3,
											bsButton("dens_mass","Get mass estimates",style="success"),
											numericInput("boot_size", "Size of bootstrap sample:",min=10, 200),
											radioButtons("use_weight", "Weight by intensity?", c("no"="no","yes"="yes"))
										),
										tags$h5("m/z estimate (dark blue):"),
										textOutput("prof_mass"),
										imageOutput("massdens", height="auto"),
										imageOutput("massint", height="auto"),
										value="test3"),
									bsCollapsePanel("Profile peak table", 
										DT::dataTableOutput("oneproftable"),
										value="test2")
								)
						),
						tabPanel("Normalization",            
								imageOutput("profnorm", height="auto"),
								imageOutput("profcount", height="auto")
						)#,
						#id="navbar_prof",inverse=FALSE,collapsible=TRUE,fluid=TRUE
					)
				),
				tabPanel("Quality control",
					tabsetPanel(
						tabPanel("Positive ionization ",
							tags$h5("Quantile distribution of peak intensities:"),           
							imageOutput("plotQCa_pos", height="auto"),
							tags$h5("Outliers:"),
							imageOutput("plotQCb_pos", height="auto"),
							tags$h5("Intensity distribution:"),                    
							imageOutput("pic_int_distr_pos", width = "100%", height = "250px")
						),
						tabPanel("Negative ionization ",
							tags$h5("Quantile distribution of peak intensities:"),           
							imageOutput("plotQCa_neg", height="auto"),
							tags$h5("Outliers:"),
							imageOutput("plotQCb_neg", height="auto"),
							tags$h5("Intensity distribution:"),                    
							imageOutput("pic_int_distr_neg", width = "100%", height = "250px")
						)
					)	
				),
                tabPanel("EIC & Peaks",
					div(style = widget_style3,selectInput("sel_meas_ID", "Select file ID:", choices = c("none"), "none")),
					div(style = widget_style3,numericInput("sel_peak_ID", "Select peak ID:", 0)),
					imageOutput("EIC1", height="auto"),
					imageOutput("EIC2", height="auto"),
					imageOutput("EIC3", height="auto")
                ),
                tabPanel("Processing",            
					selectInput("sel_meas", "Select file ID:", choices = c("none"), "none"),
					imageOutput("recal_pic", height="auto"),
					imageOutput("peakhist_pic", height="auto"),
					imageOutput("peakmzRT_pic", height="auto")
                ),
				
				
				tabPanel("Compound screening",
					tabsetPanel(
						tabPanel("Positive ionization",
							selectInput(inputId="Pos_compound_select",label="",choices=c("Target compounds","Internal standards"), 
								selected = "Target compounds", multiple = FALSE),
			
				
							fluidRow(
								column(9, DT::dataTableOutput('Table_IS_screening_pos')),
								column(3, verbatimTextOutput('Table_IS_screening_pos_row'))
							)
						
						
						),
						tabPanel("Negative ionization",
							selectInput(inputId="Neg_compound_select",label="",choices=c("Target compounds","Internal standards"), selected = "Target compounds", multiple = FALSE)
	
						
						)
					)	
				)

				
            )
        ),
        ########################################################################
        # HELP #################################################################
        ########################################################################
        tabPanel("Manual",
			tags$iframe(style="height:800px; width:110%;", src="manual.pdf")
		),
        ########################################################################
        # ABOUT ################################################################
        ########################################################################
        tabPanel("About",
		   tags$h5("Citing enviMass"),
		   HTML('<p>
			Loos, M., Ruff, M., Singer, H., 2013. enviMass v3.2 - Software workflow for the monitoring of temporal micropollutant dynamics using LC-HRMS data
			</p> '),
		   tags$h5("Contact, maintainer:"),
		   HTML('<p> Martin Loos, Martin.Loos@eawag.ch </p> '),
		   tags$h5("License:"),
		   HTML('<p>GPL-2 </p> ')
        )
        ########################################################################
      ),
	  HTML('<font color="white">') # camouflage of TRUE from sourcing
    )

