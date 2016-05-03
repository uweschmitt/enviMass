	  
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
			HTML('<hr noshade="noshade" />'),
				bsCollapse(multiple = FALSE, open = "files_open", id = "files",
					bsCollapsePanel("Add LC-HRMS file", 		
						helpText("To add a new file, set the below specifications accordingly and select it."),
						HTML('<hr noshade="noshade" />'),												
						fluidRow(
							column(width = 5, numericInput("Measadd_ID", "Numeric ID:", value = "123")),
							column(width = 5, radioButtons("Measadd_ID_autom", "...or assign ID automatically (highly recommended)?", c("yes"="yes","no"="no")))							
						),
						HTML('<hr noshade="noshade" />'),
						fluidRow(
							column(width = 5, textInput("Measadd_name", "Name:", value = "Sample 1")),
							column(width = 5, selectInput("Measadd_type", "Type:", choices = c("sample", "blank", "doted", "other"))),
							column(width = 5, selectInput("Measadd_incl", "Include?", choices = c("TRUE","FALSE"))),
							column(width = 5, selectInput("Measadd_mode", "Choose ionization mode:", choices = c("positive", "negative")))	
							
						),
						HTML('<hr noshade="noshade" />'),
						fluidRow(
							column(width = 5,textInput("Measadd_place", "Place:", value = "Rhine")),		
							column(width = 5,dateInput("Measadd_date", "Date", value = NULL, min = NULL,max = NULL, format = "yyyy-mm-dd", startview = "month",weekstart = 0, language = "en")),	
							column(width = 5,textInput("Measadd_time", "Time:(HH:MM:SS)", value = "12:00:00")),							
							column(width = 5,textInput("Measadd_tag3", "Replicate group", value = "FALSE"))	
						),
						HTML('<hr noshade="noshade" />'),
						fluidRow(
							column(width = 5, selectInput("Measadd_profiled", "Use for profiling?", choices = c("TRUE","FALSE"), selected = TRUE))
						),						
						HTML('<hr noshade="noshade" />'),
						div(style = widget_style,
							fileInput("Measadd_path", "Select centroided .mzXML file:", multiple = FALSE, accept = c(".mzXML",".raw")),
							bsPopover("Measadd_path", 
								title = "WARNING",
								content = "Files must be centroided. Check Package enviPick whether peaks can be picked properly from your files. Tested with Orbitrap files only.
								If reload fails, press cancel in the file upload window first", 
								placement = "right", trigger = "hover"),
							textOutput("had_meas_added")		
						)
					),
					bsCollapsePanel("Delete LC-HRMS file", 		
						tags$h5("Delete file by its unique ID from the below file table"),
						textInput("Measdel_ID", "ID:", value = "123"),
						actionButton("Measdel","Remove")		
					),				
					bsCollapsePanel("Modify a file specification", 
						fluidRow(
							column(width = 4, helpText("Load settings of a file into below mask by its ID, modify and then export the new settings into the main table. Modifications make a full recalculation default.")),
							column(width = 3, textInput("Modif_ID", "ID:", value = "123"), actionButton("Modif_load","Load"))
						),
						HTML('<hr noshade="noshade" />'),
						fluidRow(
							column(width = 5, textInput("Modif_name", "Name:", value = "Sample 1")),
							column(width = 5, selectInput("Modif_type", "Type:", choices = c("sample", "blank", "doted", "other"))),
							column(width = 5, selectInput("Modif_mode", "Choose ionization mode:", choices = c("positive", "negative")))	
							
						),
						HTML('<hr noshade="noshade" />'),
						fluidRow(
							column(width = 5,textInput("Modif_place", "Place:", value = "Rhine")),		
							column(width = 5,dateInput("Modif_date", "Date", value = NULL, min = NULL,max = NULL, format = "yyyy-mm-dd", startview = "month",weekstart = 0, language = "en")),	
							column(width = 5,textInput("Modif_time", "Time:(HH:MM:SS)", value = "12:00:00")),							
							column(width = 5,textInput("Modif_tag3", "Replicate group", value = "FALSE"))	
						),
						HTML('<hr noshade="noshade" />'),
						fluidRow(
							column(width = 5, selectInput("Modif_include","Include in workflow?",choices = c("TRUE","FALSE"),selected="TRUE")),
							column(width = 5, selectInput("Modif_profiled","Use for profiling?",choices = c("TRUE","FALSE"),selected="TRUE"))
						),	
						HTML('<hr noshade="noshade" />'),
						actionButton("Modif_export","Save")
					),
					bsCollapsePanel("Batch upload from folder", 	
						tags$h5("Read in batches of files (.mzXML) from a folder; file specifications will be guessed and can later be modified above. 
						File dates are set in the sequence of files provided in the folder."),	
						textInput("import_file_folder", "Folder path:", value = "C:\\...\\folder"),
						bsPopover("import_file_folder", 
							title = "Insert full path, including the project folder, but excluding the logfile.emp.",
							content = "Using your OS explorer, you may navigate into your project folder and copy/paste the full path.", 
							placement = "right", trigger = "hover"),
						checkboxInput("Import_file_folder_overwrite", "Import files with same name?", FALSE),
						bsPopover("Import_file_folder_overwrite", 
							title = "File import to existing files",
							content = "If a file in the folder with the same filename as one already existing in the project is found, should it be imported (new ID assigned)?", 
							placement = "right", trigger = "hover"),
						actionButton("Import_file_folder","Import"),
						HTML('<hr noshade="noshade" />'),
						textOutput("had_import_folder")						
					),
					bsCollapsePanel("Import files from another project", 		
						tags$h5("Select project folder to import files from:"),
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
						actionButton("Import_project","Import"),		
						HTML('<hr noshade="noshade" />'),
						textOutput("had_import_project")	
					)					
				),	
			HTML('<hr noshade="noshade" />'),
			DT::dataTableOutput("measurements")	
	    ),
        ########################################################################
        # Compounds ############################################################
        ########################################################################
        tabPanel("Compounds",
		 tabsetPanel(
            tabPanel("Internal standards",
				HTML('<hr noshade="noshade" />'),
				bsCollapse(multiple = FALSE, open = "IS_comp_open", id = "IS_comp",
					bsCollapsePanel("Add internal standard compound", 
						fluidRow(
							column(width = 5, helpText("To add a new internal standard, fill out the below form and press Add") ),
							column(width = 2, offset = 0.3, actionButton("AddIS","Add"))
						),
						HTML('<hr noshade="noshade" />'),
						div(style = widget_style,
							textInput("ISadd_ID", "Unique ID:", value = "123_XYZ"),          
							textInput("ISadd_name", "Name:", value = "CompoundX"),
							textInput("ISadd_formula", "Formula:", value = "C10H12O10")
						),
						div(style = widget_style,
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
							checkboxInput("ISadd_use_screen", "To be used for screening?", TRUE)
						),
						div(style = widget_style,
							textInput("ISadd_remark", "remark", value = "none"),
							textInput("ISadd_tag1", "tag1", value = "none"),
							textInput("ISadd_tag2", "tag2", value = "none"),
							textInput("ISadd_tag3", "tag3", value = "none")
						),
						div(style = widget_style,							
							checkboxInput("ISadd_date", "Restrict to temporal range?", FALSE),	
							dateRangeInput("ISadd_date_range", label="", start = NULL, end = NULL,
								min = NULL, max = NULL, format = "yyyy-mm-dd",
								startview = "month", weekstart = 0, language = "en",
								separator = " to "),
							textInput("Lower_intensity_bound", "Screening: lower intensity bound (.-separated)", value = "0"),						
							textInput("Upper_intensity_bound", "Screening: upper intensity bound (.-separated)", value = "Inf")
						)
					),
					bsCollapsePanel("Remove internal standard compound", 					
						helpText("To delete a compound from the list, type in its ID and press Delete"),
						textInput("ISdelete_ID", "ID for deletion:", value = "123_XYZ"),          
						actionButton("DeleteIS","Delete")					
					),
					bsCollapsePanel("Import compound list", 					
						fileInput("ISlist_path", "Select IS.txt file, e.g. from dataframes folder of another project", multiple = FALSE, accept = c(".txt"))
					),
					bsCollapsePanel("Modify in external editor", 					
						HTML('
							<p><font>
								The below compound table can be assembled and modified in external text editors or Excel and then imported via the above import step.
								To do so, use file IS.txt from the dataframe folder of a new enviMass project as a template. 
								Mind the character encoding (e.g., ANSI) of .txt files when modifying.
							</font></p>
							<p style="background-color:darkred"; align="center"> <font color="#FFFFFF"> 
								Any such modifications must strictly adhere to the following rules to avoid frustration:
							</font></p> 
							<ol>
								<li>Headers: strictly use file IS.txt generated in a new enviMass project (dataframes folder) as a template. No empty spaces in header names permitted (e.g. ion_mode NOT ion mode).</li>
								<li>Input format: text file (.txt), tab delimited.</li>
								<li>Compound names: no special signs permitted. Use big and small letters, numbers, underscores, hyphen, brackets and empty spaces - and absolutely nothing else.</li>
								<li>No empty columns. If you are not sure what to fill in, use what is given  in the template IS.txt file.</li>
								<li>Absolutely NO duplicated IDs.</li>
								<li>Do not delete columns; their number, order and content types are all fixed.</li>
								<li>No uncompleted entries per compounds.</li>
								<li>Numeric entries with decimal points: dot-separated.</li>	
							</ol>
							<p style="background-color:darkgreen"; align="center"> <font color="#FFFFFF"> 
								Help for some column contents:
							</font></p> 
							<ol>
								<li>ID: numbers and characters permitted; no empty spaces or special characters. Unique, absolutely NO duplicates permitted.</li>
								<li>RT_tolerance: FALSE or a compound-specific retention time tolerance. Overwrites the one set as standard value in tab Settings/Screening/IS.</li>
								<li>main_adduct: FALSE or name of a special adduct to be used for this compound entry. Valid adduct names can be found in tab Settings/Adduct.</li>								
								<li>restrict_adduct: TRUE or FALSE. Only use the main_adduct (if specified) for this compound and ignore the ones specified in tab Settings/Adduct?</li>
							</ol>
						')		
					)
				),
                HTML('<hr noshade="noshade" />'),
                DT::dataTableOutput("IS")
            ),
            tabPanel("Targets", 
				HTML('<hr noshade="noshade" />'),
				bsCollapse(multiple = FALSE, open = "target_comp_open", id = "target_comp",
					bsCollapsePanel("Add target compound", 
						fluidRow(
							column(width = 5, helpText("To add a new target, fill out the below form and press Add.") ),
							column(width = 2, offset = 0.3,  actionButton("Addtargets","Add"))
						),			
						div(style = widget_style,
							textInput("targetsadd_ID", "Unique ID:", value = "123_XYZ"),          
							textInput("targetsadd_name", "Name:", value = "CompoundY"),
							textInput("targetsadd_formula", "Formula:", value = "C10H12O10")          
						),
						div(style = widget_style,						
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
							checkboxInput("targetsadd_use_screen", "To be used for screening?", TRUE)
						),
						div(style = widget_style,							
							textInput("targetsadd_remark", "remark", value = "none"),
							textInput("targetsadd_tag1", "tag1", value = "none"),
							textInput("targetsadd_tag2", "tag2", value = "none"),
							textInput("targetsadd_tag3", "tag3", value = "none")
						),
						div(style = widget_style,
							checkboxInput("targetsadd_date", "Restrict to temporal range?", FALSE),	
							dateRangeInput("targetsadd_date_range", label="", start = NULL, end = NULL,
								min = NULL, max = NULL, format = "yyyy-mm-dd",
								startview = "month", weekstart = 0, language = "en",
								separator = " to "),
							textInput("warn_1", "First concentration warn level (.-separated)", value = "FALSE"),							
							textInput("warn_2", "Second concentration warn level (.-separated)", value = "FALSE")						
						)
					),
					bsCollapsePanel("Remove target compound", 
						helpText("To delete a compound from the list, type in its ID and press Delete"),
						textInput("targetsdelete_ID", "ID for deletion:", value = "123_XYZ"),          
						actionButton("Deletetargets","Delete")					
					),
					bsCollapsePanel("Import compound list", 					
						fileInput("targetlist_path", "Select targets.txt file, e.g. from dataframes folder of another project", multiple = FALSE, accept = c(".txt"))
					),
					bsCollapsePanel("Modify in external editor", 					
						HTML('
							<p><font>
								The below compound table can be assembled and modified in external text editors or Excel and then imported via the above import step.
								To do so, use file targets.txt from the dataframe folder of a new enviMass project as a template.
								Mind the character encoding (e.g., ANSI) of .txt files when modifying.
							</font></p>
							<p style="background-color:darkred"; align="center"> <font color="#FFFFFF"> 
								Any such modifications must strictly adhere to the following rules to avoid frustration:
							</font></p> 
							<ol>
								<li>Headers: strictly use file targets.txt generated in a new enviMass project (dataframes folder) as a template. No empty spaces in header names permitted (e.g. ion_mode NOT ion mode).</li>
								<li>Input format: text file (.txt), tab delimited.</li>
								<li>Compound names: no special signs permitted. Use big and small letters, numbers, underscores, hyphen, brackets and empty spaces - and absolutely nothing else.</li>
								<li>No empty columns. If you are not sure what to fill in, use what is given  in the template targets.txt file.</li>
								<li>Absolutely NO duplicated IDs.</li>
								<li>Do not delete columns; their number, order and content types are all fixed.</li>
								<li>No uncompleted entries per compounds.</li>
								<li>Numeric entries with decimal points: dot-separated.</li>								
							</ol>
							<p style="background-color:darkgreen"; align="center"> <font color="#FFFFFF"> 
								Help for some column contents:
							</font></p> 
							<ol>
								<li>ID: numbers and characters permitted; no empty spaces or special characters. Unique, absolutely NO duplicates permitted.</li>
								<li>RT_tolerance: FALSE or a compound-specific retention time tolerance. Overwrites the one set as standard value in tab Settings/Screening/targets.</li>
								<li>ID_internal_standard: unique ID of a internal standard compound to be used for quantification. Set to FALSE otherwise.</li>
								<li>main_adduct: FALSE or name of a special adduct to be used for this compound entry. Valid adduct names can be found in tab Settings/Adduct.</li>								
								<li>restrict_adduct: TRUE or FALSE. Only use the main_adduct (if specified) for this compound and ignore the ones specified in tab Settings/Adduct?</li>
							</ol>
						')		
					)
				),
				HTML('<hr noshade="noshade" />'),  
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
			HTML('<hr noshade="noshade" />'),
			tags$p(align="justify","Below is an illustrative network graph of dependencies between steps in the workflow. Enabled steps are shown in 
				dark blue; disabled ones in light blue; interactive and zoomable. When modifying workflow steps, parameters or inputs, enviMass 
				dynamically adjusts and minimizes all required recalculations via their relevant dependencies. These recalculations will ultimately be 
				enforced when pressing the left-sided Calculate button."),
				networkD3:::forceNetworkOutput("force_workflow", width = 1500, height = 400),
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
					fluidRow(
						column(width = 2, radioButtons("recal", "Include?", c("yes"="yes","no"="no")) ),
						column(width = 10, offset = 0.3,
							tags$p(align="justify","Theoretical masses of internal standard compounds are used to correct systematic offsets in measured masses.")
						)
					),
				HTML('<p style="background-color:darkblue"; align="center"> <font color="#FFFFFF"> Median intensity normalization </font></p> '),
					radioButtons("intnorm", "Include?", c("yes"="yes","no"="no")),
				#HTML('<p style="background-color:darkgreen"; align="center"> <font color="#FFFFFF"> RT alignment </font></p> '),
				#radioButtons("RTalign", "Include?", c("yes"="yes","no"="no")),  
				HTML('<p style="background-color:darkblue"; align="center"> <font color="#FFFFFF"> Blank / blind peak detection </font></p> '),				
					fluidRow(
						column(width = 2, radioButtons("blind_filter", "Detect?", c("yes"="yes","no"="no")) ),
						column(width = 2, radioButtons("blind_omit", "Remove?", c("yes"="yes","no"="no")) ),
						column(width = 8, offset = 0.3,
							tags$p(align="justify","Detects sample peaks which are also present in blind/blank files. Check Settings Blind Tab for file selection. 
							Optionally, affected peaks can early be removed before downstream processing. Recommended for e.g. effect on compound screening - otherwise
							later removal (cp. red steps) recommended.")
						)
					),				
				HTML('<p style="background-color:darkblue"; align="center"> <font color="#FFFFFF"> Replicate filter </font></p> '),				
					fluidRow(
						column(width = 2, radioButtons("replicates", "Include?", c("yes"="yes","no"="no")) ),
						column(width = 2, radioButtons("replicates_prof", "Use in profiling?", c("yes"="yes","no"="no")) ),
						column(width = 8, offset = 0.3,
							tags$p(align="justify","Filters out picked peaks which are not ubiquitously present in a set of measurements, within tolerances specified in the Settings/Replicates tab.
							Typically, such a set would be composed of replicate measurements. 
							A set can be defined by a joint string other than FALSE (e.g. replicates_A) in the tag_3 column of the file table, as assembled in the Files tab.
							Optionally, profiles can be extracted in the replicates first, and later merged over all files.")
						)
					),
				HTML('<hr noshade="noshade" />'),
				HTML('<h1 align="center"> &#x21e9; </h1> '),                     
				# block 3 ######################################################
				HTML('<p style="background-color:darkgreen"; align="center"> <font color="#FFFFFF"> Profile extraction </font></p>'),
					fluidRow(
						column(width = 2, radioButtons("profiled", "Include? ", c("yes"="yes","no"="no")) ),
						column(width = 10, offset = 0.3,
							tags$p(align="justify","An intensity descent assorts peaks into profiles. Uses a fixed retention time tolerance window and
							an adaptive mass tolerance window. Only peaks of files marked for profiling are used (column profiled in the files table / File Tab).")
						)
					),					
				HTML('<p style="background-color:darkgreen"; align="center"> <font color="#FFFFFF"> LOD interpolation </font></p>'),
					fluidRow(
						column(width = 2, radioButtons("LOD_interpol", "Include? ", c("yes"="yes","no"="no")) ),
						column(width = 10, offset = 0.3,
							tags$p(align="justify","For each measurement, estimates a RT-dependent intensity threshold below which peaks are not expected to get picked. 
							Can be used for the below compound screening.")
						)
					),
				HTML('<p style="background-color:darkgreen"; align="center"> <font color="#FFFFFF"> Compound screening </font></p> '),
					fluidRow(
						column(width = 3, 
							radioButtons("screen_IS", "Screen internal standards?", c("yes"="yes","no"="no"))
						),
						column(width = 3, 
							radioButtons("screen_target", "Screen targets/suspects?", c("yes"="yes","no"="no"))
						),						
						column(width = 6, offset = 0.3,
							tags$p(align="justify","Uses the LOD thresholds estimated in the above step. If the LOD interpolation is not run, a fixed intensity threshold as specified in the
							Settings/Screening tabs is used. All peaks of files marked for inclusion and profiling are used (column profiled in the files table / File Tab) - unless peaks are
							removed explicitly in the early blind detection (see blue workflow steps).")
						)
					),	
				HTML('<p style="background-color:darkgreen"; align="center"> <font color="#FFFFFF"> Quantification </font></p> '),
					radioButtons("quantif", "Include? ", c("yes"="yes","no"="no")),					
				HTML('<hr noshade="noshade" />'),
				HTML('<h1 align="center"> &#x21e9; </h1> '),  					
				# block 4 ######################################################
				HTML('<p style="background-color:darkred"; align="center"> <font color="#FFFFFF"> Normalization using IS-profiles </font></p> '),
					fluidRow(
						column(width = 2, radioButtons("profnorm", "Include? ", c("yes"="yes","no"="no")) ),
						column(width = 10, offset = 0.3,
							tags$p(align="justify","Relies on the above internal standard screening and profile extraction (see green steps). Intensities of picked peaks in each measurements 
							are normalized by the median deviation over all internal standards from their individual median intensity taken over their individual profiles. Thus,
							this approach differs from the less reliable Median intensity normalization above, which can be skipped in this case (see blue steps).")
						)
					),					
				HTML('<p style="background-color:darkred"; align="center"> <font color="#FFFFFF"> Compound subtraction </font></p> '),
					fluidRow(
						column(width = 3, 
							radioButtons("subtr_IS", "Subtract internal standards?", c("yes"="yes","no"="no"))
						),
						column(width = 3, 
							radioButtons("subtr_target", "Subtract targets/suspects?", c("yes"="yes","no"="no"))
						),						
						column(width = 6, offset = 0.3,
							tags$p(align="justify","Run a profile recalculation omitting compound peaks belonging to matches >= the cutoff score defined in the Settings/Screening Tab.
							The original compound screening results remain visible, but affected peaks will not be part of the final profiles.")
						)
					),	
				HTML('<p style="background-color:darkred"; align="center"> <font color="#FFFFFF"> Blind peak subtraction </font></p> '),
					fluidRow(
						column(width = 2, 
							radioButtons("subtr_blind", "Subtract?", c("yes"="yes","no"="no"))
						),	
						column(width = 10, offset = 0.3,
							tags$p(align="justify","Subtract peaks from profiles which have also been detected in blind/blank samples but not removed yet (cp. blind detection of blue steps)?")
						)
					),	
				HTML('<p style="background-color:darkred"; align="center"> <font color="#FFFFFF"> Trend detection </font></p> '),
					fluidRow(
						column(width = 2, radioButtons("trenddetect", "Include? ", c("yes"="yes","no"="no")) ),
						column(width = 10, offset = 0.3,
							tags$p(align="justify","Contains a separate blind detection step. Herein, intensities of blind/blank peaks are interpolated over the time series.
							This interpolation and subtraction is only applicable if the separate blind filter step is disabled (see above blue steps).")
						)
					),							
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
				numericInput("replicate_dmz", "+/- m/z tolerance ...", 3),                
				selectInput("replicate_ppm", "... given in:", choices = c("ppm"="TRUE","absolute"="FALSE"), "TRUE"),	
				#selectInput("replicate_recalib", "... and corrected by recalibration results (if available)", choices = c("TRUE"="TRUE","FALSE"="FALSE"), "FALSE"),	
				numericInput("replicate_delRT", "RT tolerance window of peaks caused by the same analyte across replicate samples [s]", 30),
				sliderInput("replicate_IS_dInt", "Intensity tolerance %", min = 0, max = 100, value = 30, step= .2)
			),	
            # ALLIGNMENT #######################################################
            #tabPanel("Alignment",
            #  tags$h5("RT alignment")
            #),
            # SCREENING ########################################################
            tabPanel("Screening",
              tabsetPanel(
                tabPanel("IS",
					div(style = widget_style2,
						tags$h5("Retention time"),
						numericInput("screen_IS_delRT", "RT tolerance of peaks in sample relative to their expected RT [s]", 30),
						numericInput("screen_IS_dRTwithin", "RT tolerance of peaks within an isotope pattern [s]", 50)
					),
				  	div(style = widget_style2,
						tags$h5("Mass"),
						numericInput("screen_IS_dmz", "m/z tolerance ...", 3),                
						selectInput("screen_IS_ppm", "... given in:", choices = c("ppm"="TRUE","absolute"="FALSE"), "TRUE")
 					),
					div(style = widget_style2,
						tags$h5("Intensity"),
						sliderInput("screen_IS_dInt", "Intensity tolerance %", min = 0, max = 100, value = 30, step= .2),
						numericInput("screen_IS_Intcut", "Lower intensity threshold (if LOD interpolation disabled)", 5E4)                
					),
					div(style = widget_style2,
						tags$h5("Scoring"),
						numericInput("screen_IS_w1", "Cutoff score [0,1]", 0.8)              
					)
                ),
                tabPanel("Targets & Suspects",
					div(style = widget_style2,
						tags$h5("Retention time"),
						numericInput("screen_target_delRT", "RT tolerance of peaks in sample relative to their expected RT [s]", 30),
						numericInput("screen_target_dRTwithin", "RT tolerance of peaks within an isotope pattern [s]", 50)
					),
				  	div(style = widget_style2,
						tags$h5("Mass"),
						numericInput("screen_target_dmz", "m/z tolerance ...", 3),                
						selectInput("screen_target_ppm", "... given in:", choices = c("ppm"="TRUE","absolute"="FALSE"), "TRUE")
 					),
					div(style = widget_style2,
						tags$h5("Intensity"),
						sliderInput("screen_target_dInt", "Intensity tolerance %", min = 0, max = 100, value = 30, step= .2),
						numericInput("screen_target_Intcut", "Lower intensity threhold", 5E4)                
					),
					div(style = widget_style2,
						tags$h5("Scoring"),
						numericInput("screen_target_w1", "Cutoff score [0,1]", 0.8)           
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
                numericInput("prof_drt", "Peak deviation within profiles: RT tolerance [s]", 60),     
				HTML('<hr noshade="noshade" />'),				
				selectInput("prof_select", "Omit files with table entry profiled=FALSE from profiling?", choices = c("TRUE"="TRUE","FALSE"="FALSE"), selected="FALSE")			
			),
            # TREND ############################################################
            tabPanel("Trends",
				tags$h5("Trend detection:"),
				textInput("trend_lags", "Time lags of trends [days], comma-separated:", value = "4,7,14"),
				numericInput("trend_thres", "Trend vs. mean+variance intensity threshold:", 3),
				HTML('<hr noshade="noshade" />'),
				radioButtons("trend_blind", "Additional blind interpolation and subtraction per profile?", c("yes"="yes","no"="no")),
				HTML('<hr noshade="noshade" />'),
				selectInput("notrend", "Do not show global trend - instead, report it as maximum intensity above blind", choices = c("TRUE"="TRUE","FALSE"="FALSE"), "FALSE")
            ),
            # BLIND #############################################################
            tabPanel("Blind",
				tags$h5("Blind subtraction:"),
				numericInput("blind_fold", "Intensity threshold ratio sample/blind <", 100),
				numericInput("blind_dmz", "Mass uncertainty (+/-) ...", 3), 
                selectInput("blind_ppm", "... given in:", choices = c("ppm"="TRUE","absolute"="FALSE"), "TRUE"),				
                numericInput("blind_drt", "RT tolerance  [s]", 60),       
				HTML('<hr noshade="noshade" />'),
				tags$h5("Positive ionization mode"),
				checkboxInput("subtract_pos_bydate", "Subtract with the next blank/blind file preceding a sample by its date & time?", FALSE),				
				checkboxInput("subtract_pos_byfile", "Additional non-sample files to subtract each sample file with (i.e. not preceding by date only), choose file ID:", FALSE),				
				checkboxGroupInput("files_pos_select_subtract", label="", choices=c("FALSE"), selected = NULL),
				HTML('<hr noshade="noshade" />'),
				tags$h5("Negative ionization mode"),
				checkboxInput("subtract_neg_bydate", "Subtract with the next blank/blind file preceding a sample by its date & time?", FALSE),				
				checkboxInput("subtract_neg_byfile", "Additional non-sample files to subtract each sample file with (i.e. not preceding by date only), choose file ID:", FALSE),							
				checkboxGroupInput("files_neg_select_subtract", label="", choices=c("FALSE"), selected = NULL)
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
			),
			# IMPORT PARAMETERS FROM ANOTHER PROJECT ###########################
            tabPanel("Import",			
				tags$h5("Import all parameter settings from another project, excluding blind file selection."),
						textInput("import_pro_dir_paras", "", value = "C:\\...\\other_project_name"),
						bsPopover("import_pro_dir_paras", 
							title = "Insert full path, including the project folder, but excluding the logfile.emp.",
							content = "Using your OS explorer, you may navigate into your project folder and copy/paste the full path.", 
							placement = "right", trigger = "hover"),
						actionButton("Import_project_para","Import")		
			)
          )
        ),
        ########################################################################
        # RESULTS ##############################################################
        ########################################################################
        tabPanel("Results", 	
			tabsetPanel( 
				######################################################################################################################
				tabPanel("Profiles",
					div(style = widget_style5,
						textOutput("had_ion"),	
						selectInput("Ion_mode", label=NULL, c("positive","negative"), selected = ("positive"), multiple = FALSE)
					),
					HTML('<hr noshade="noshade" />'),  
					#navbarPage("", 
					tabsetPanel( 

### Baustelle						
						tabPanel("Summary",										
								tags$h5("Filter profile list:"),
								div(style = widget_style3,numericInput("filterProf_minmass", "Minimum m/z:", 0)),
								div(style = widget_style3,numericInput("filterProf_maxmass", "Maximum m/z:", 3000)),
								div(style = widget_style3,numericInput("filterProf_minrt", "Minimum RT [s]:", 0)),
								div(style = widget_style3,numericInput("filterProf_maxrt", "Maximum RT [s]:", 100000)),			
								div(style = widget_style3,radioButtons("filterProf_meanblind", "Use mean above blind?", c("no"="no","yes"="yes"))),
								bsPopover("filterProf_meanblind", 
									title = "Replicates, not time series ...",
									content = "Get profiles with mean sample intensity x times above mean blank intensities (set Sort profile list by: maximum or mean intensity; x = to be set in blind settings panel). Useful if all your sample files are replicates and not a time sequences.", 
									placement = "top", trigger = "hover"),
								div(style = widget_style3,radioButtons("filterProf_notblind", "Not in blind?", c("no"="no","yes"="yes"))),
								div(style = widget_style3,selectInput("filterProf_sort", "Sort profile list by:", 
									choices = c("ID","mean m/z","mean RT","maximum intensity","mean intensity","global trend intensity","current trend intensity"), selected="current trend intensity")),
								div(style = widget_style3,numericInput("filterProf_count", "Restrict list size:", 500)),
								conditionalPanel( # IS filter				
										condition = "input.screen_IS == 'yes'",
										tags$h5("IS compounds filter:"),										
										HTML('<hr noshade="noshade" />')
								),
								conditionalPanel( # target filter				
										condition = "input.screen_target == 'yes'",	
										tags$h5("Target compounds filter:"),	
										HTML('<hr noshade="noshade" />')
								),															
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
### Baustelle	
						),
						tabPanel("Latest trends",
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
				######################################################################################################################
				tabPanel("Quality control",
					tabsetPanel(
						tabPanel("Positive ionization ",
							tags$h5("Quantile distribution of peak intensities:"),           
							imageOutput("plotQCa_pos", height="auto"),
							tags$h5("Outliers:"),
							imageOutput("plotQCb_pos", height="auto"),
							tags$h5("Intensity distribution for median intensity normalization:"),                    
							imageOutput("pic_int_distr_pos", width = "100%", height = "250px")
						),
						tabPanel("Negative ionization ",
							tags$h5("Quantile distribution of peak intensities:"),           
							imageOutput("plotQCa_neg", height="auto"),
							tags$h5("Outliers:"),
							imageOutput("plotQCb_neg", height="auto"),
							tags$h5("Intensity distribution for median intensity normalization:"),                    
							imageOutput("pic_int_distr_neg", width = "100%", height = "250px")
						)
					)	
				),
				######################################################################################################################
                tabPanel("EIC & Peaks",
					div(style = widget_style3,selectInput("sel_meas_ID", "Select file ID:", choices = c("none"), "none")),
					div(style = widget_style3,numericInput("sel_peak_ID", "Select peak ID:", 0)),
					imageOutput("EIC1", height="auto"),
					imageOutput("EIC2", height="auto"),
					imageOutput("EIC3", height="auto")
                ),
                tabPanel("Processing",            
					selectInput("sel_meas", "Select file ID:", choices = c("none"), "none"),
					HTML('<hr noshade="noshade" />'),
					fluidRow(										
						column(4,tags$h5("File name:"),textOutput('file_proc_name') ),
						column(4,tags$h5("File type: "),textOutput('file_proc_type')),										
						column(4,tags$h5("Ionization mode: "),textOutput('file_proc_mode'))										
					),									
					HTML('<hr noshade="noshade" />'),
					fluidRow(										
						column(4,tags$h5("Number of peaks: "),textOutput('file_peak_number') ),
						column(4,tags$h5("% of peaks affected by blind filter: "),textOutput('file_blind_rem')),										
						column(4,tags$h5("% of peaks removed by replicate filter: "),textOutput('file_repl_rem'))										
					),					
					HTML('<hr noshade="noshade" />'),
					imageOutput("recal_pic", height="auto"),
					HTML('<hr noshade="noshade" />'),
					imageOutput("peakhist_pic", height="auto"),
					HTML('<hr noshade="noshade" />'),
					imageOutput("peakmzRT_pic", height="auto"),
					HTML('<hr noshade="noshade" />'),
					imageOutput("LOD_pic", height="auto"),
					HTML('<hr noshade="noshade" />')
				),
				######################################################################################################################
				tabPanel("Compound screening",
					tabsetPanel(
						tabPanel("Positive ionization",
							fluidRow(
								column(4, selectInput(inputId="Pos_compound_select",label="",choices=c("Target compounds","Internal standards"), 
									selected = "Target compounds", multiple = FALSE)),
								column(4, selectInput(inputId="screen_pos_summarize", label="", choices = c("Show all adducts"="yes","Collapse adducts"="no"), "yes"))
							),
							bsCollapse(multiple = FALSE, open = NULL, id = "collapse_screen_pos_one",
								bsCollapsePanel(title="Pattern match for selected compound", #style="info",
									textOutput('screening_details_comp_pos'),
									plotOutput("plot_pattern_pos")
								),
								bsCollapsePanel(title="Characteristics for selected compound",
									textOutput('screening_details_comp_pos2'),
									HTML('<hr noshade="noshade" />'),
									fluidRow(										
										column(4, selectInput(inputId="selec_pos_x",label="x axis",
											choices=c("m/z","RT","Intensity","Date&time","Type","Place","Conz."),selected = "m/z", multiple = FALSE)),
										column(4, selectInput(inputId="selec_pos_y",label="y axis",
											choices=c("m/z","RT","Intensity","Date&time","Type","Place","Conz."),selected = "RT", multiple = FALSE)),									
										column(4, radioButtons("selec_pos_log_rat", "Log intensity?", c("yes"="yes","no"="no"),inline=TRUE))
									),
									HTML('<hr noshade="noshade" />'),						
									plotOutput("plot_selec_dist_pos"),
									conditionalPanel(				
										condition = "input.Pos_compound_select == 'Internal standards'",					
										HTML('<hr noshade="noshade" />'),
										fluidRow(										
											column(4,tags$p(align="justify","Adopt a new log intensity range for the future screening of this selected compound-adduct?")),
											column(3,numericInput("screen_int_pos_low", "Lower bound", 0,step=0.1)),
											column(3,numericInput("screen_int_pos_up", "Upper bound", 10,step=0.1)),
											column(2,bsButton("save_int_pos"," Adopt",style="warning",icon=icon("bookmark")))
										)
									)
								),
								bsCollapsePanel(title="Screening table for selected compound", 
									textOutput('screening_details_comp_pos3'),
									HTML('<hr noshade="noshade" />'),
									DT::dataTableOutput('Table_screening_selected_pos')
								)
							),
							HTML('<hr noshade="noshade" />'),
							tags$p(align="justify","The below sample and blank matches give the number of files with matches above the cutoff score, 
							with multiple matches per file above this cutoff merged."),
							DT::dataTableOutput('Table_screening_pos'),
							HTML('<hr noshade="noshade" />'),
							bsCollapse(multiple = FALSE, open = NULL, id = "collapse_screen_pos_all",
								bsCollapsePanel(title="Summary plots",
									fluidRow(
										column(width = 4, offset = 0.6,
												tags$p(align="justify","Characteristics of all signal peaks for screened compounds from the above table.")
										),
										column(4, selectInput(inputId="Summ_pos_x",label="x axis",choices=c("m/z","RT","log Intensity","m/z deviation [ppm]","RT deviation","Time sequence"),selected = "m/z", multiple = FALSE)),
										column(4, selectInput(inputId="Summ_pos_y",label="y axis",choices=c("m/z","RT","log Intensity","m/z deviation [ppm]","RT deviation","Time sequence"),selected = "RT", multiple = FALSE))									
									),
									plotOutput("plot_pattern_distrib_pos"),
									HTML('<hr noshade="noshade" />'),
									fluidRow(										
										column(4,tags$p(align="justify","Ratios of sample vs. blank intensities for screened compounds from the above table.")),
										column(4,radioButtons("screen_pos_log_rat", "Log scale?", c("yes"="yes","no"="no"),inline=TRUE))
									),
									plotOutput("plot_aboveBlank_pos",height = 250)
								)
							)	
						),
						tabPanel("Negative ionization",
							fluidRow(
								column(4, selectInput(inputId="Neg_compound_select",label="",choices=c("Target compounds","Internal standards"), 
									selected = "Target compounds", multiple = FALSE)),
								column(4, selectInput(inputId="screen_neg_summarize", label="", choices = c("Show all adducts"="yes","Collapse adducts"="no"), "yes"))
							),
							bsCollapse(multiple = FALSE, open = NULL, id = "collapse_screen_pos_one",
								bsCollapsePanel(title="Pattern match for selected compound", #style="info",
									textOutput('screening_details_comp_neg'),
									plotOutput("plot_pattern_neg")
								),
								bsCollapsePanel(title="Characteristics for selected compound",
									textOutput('screening_details_comp_neg2'),
									HTML('<hr noshade="noshade" />'),
									fluidRow(										
										column(4, selectInput(inputId="selec_neg_x",label="x axis",
											choices=c("m/z","RT","Intensity","Date&time","Type","Place","Conz."),selected = "m/z", multiple = FALSE)),
										column(4, selectInput(inputId="selec_neg_y",label="y axis",
											choices=c("m/z","RT","Intensity","Date&time","Type","Place","Conz."),selected = "RT", multiple = FALSE)),									
										column(4, radioButtons("selec_neg_log_rat", "Log intensity?", c("yes"="yes","no"="no"),inline=TRUE))
									),
									HTML('<hr noshade="noshade" />'),						
									plotOutput("plot_selec_dist_neg"),
									conditionalPanel(				
										condition = "input.Neg_compound_select == 'Internal standards'",					
										HTML('<hr noshade="noshade" />'),
										fluidRow(										
											column(4,tags$p(align="justify","Adopt a new log intensity range for the future screening of this selected compound-adduct?")),
											column(3,numericInput("screen_int_neg_low", "Lower bound", 0,step=0.1)),
											column(3,numericInput("screen_int_neg_up", "Upper bound", 10,step=0.1)),
											column(2,bsButton("save_int_neg"," Adopt",style="warning",icon=icon("bookmark")))
										)
									)
								),
								bsCollapsePanel(title="Screening table for selected compound", 
									textOutput('screening_details_comp_neg3'),
									HTML('<hr noshade="noshade" />'),
									DT::dataTableOutput('Table_screening_selected_neg')
								)
							),
							HTML('<hr noshade="noshade" />'),
							tags$p(align="justify","The below sample and blank matches give the number of files with matches above the cutoff score, 
							with multiple matches per file above this cutoff merged."),
							DT::dataTableOutput('Table_screening_neg'),
							HTML('<hr noshade="noshade" />'),
							bsCollapse(multiple = FALSE, open = NULL, id = "collapse_screen_neg_all",
								bsCollapsePanel(title="Summary plots",
									fluidRow(
										column(width = 4, offset = 0.6,
												tags$p(align="justify","Characteristics of all signal peaks for screened compounds from the above table.")
										),
										column(4, selectInput(inputId="Summ_neg_x",label="x axis",choices=c("m/z","RT","log Intensity","m/z deviation [ppm]","RT deviation","Time sequence"),selected = "m/z", multiple = FALSE)),
										column(4, selectInput(inputId="Summ_neg_y",label="y axis",choices=c("m/z","RT","log Intensity","m/z deviation [ppm]","RT deviation","Time sequence"),selected = "RT", multiple = FALSE))									
									),
									plotOutput("plot_pattern_distrib_neg"),
									HTML('<hr noshade="noshade" />'),
									fluidRow(										
										column(4,tags$p(align="justify","Ratios of sample vs. blank intensities for screened compounds from the above table.")),
										column(4,radioButtons("screen_neg_log_rat", "Log scale?", c("yes"="yes","no"="no"),inline=TRUE))
									),
									plotOutput("plot_aboveBlank_neg",height = 250)
								)
							)	
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
			Loos, M., Ruff, M., Singer, H., 2013. enviMass v3.1 - Software workflow for the monitoring of temporal micropollutant dynamics using LC-HRMS data
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

