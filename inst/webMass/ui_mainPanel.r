	  
    conditionalPanel( 
		condition = "output.textit != 'Waiting...'",  	
		tags$h5(""),
		bsAlert("alert_3"),
		tabsetPanel(
		#navbarPage("",
        ########################################################################
        # MEASUREMENTS #########################################################
        ########################################################################
		tabPanel("Files",
			HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/inputs/files.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; Check help for details</a></p>'),	
				bsCollapse(multiple = FALSE, open = "files_open", id = "files",
					# ADD FILE #################################################
					bsCollapsePanel("Add LC-MS file", 		
						helpText("To add a new file.mzXML, set its specifications below and upload it."),
						HTML('<hr noshade="noshade" />'),
						fluidRow(
							column(width = 5, textInput("Measadd_name", "Name:", value = "File XY")),
							column(width = 5, selectInput("Measadd_type", "Type:", choices = c("sample", "blank", "calibration", "spiked"))),
							column(width = 5, selectInput("Measadd_incl", "Include?", choices = c("TRUE","FALSE"))),
							column(width = 5, selectInput("Measadd_mode", "Choose ionization mode:", choices = c("positive", "negative")))	
						),
						conditionalPanel(
							condition = "input.Measadd_type != 'calibration' & input.Measadd_type != 'spiked'",						
							HTML('<hr noshade="noshade" />'),
							fluidRow(
								column(width = 5,textInput("Measadd_place", "Place:", value = "Rhine")),		
								column(width = 5,dateInput("Measadd_date", "Date:", value = NULL, min = NULL,max = NULL, format = "yyyy-mm-dd", startview = "month",weekstart = 0, language = "en")),	
								column(width = 5,textInput("Measadd_time", "Time (HH:MM:SS):", value = "12:00:00")),							
								column(width = 5,
									conditionalPanel(
										condition = "input.Measadd_type == 'samples'",								
											textInput("Measadd_tag3", "Replicate group (tag3):", value = "FALSE")
									)
								),
								column(width = 5,textInput("Measadd_ID2", "Custom ID:", value = "FALSE")),
								column(width = 5, selectInput("Measadd_profiled", "Use for profiling (if Settings/Profiling/Omit adjusted to do so)?", choices = c("TRUE","FALSE"), selected = "TRUE"))
							)			
						),	
						conditionalPanel(
							condition = "input.Measadd_type == 'calibration'", 
							HTML('<hr noshade="noshade" />'),
							fluidRow(
								column(width = 5,textInput("Measadd_tag1", "Concentration (no units; tag1)", value = "FALSE")),
								column(width = 5,textInput("Measadd_tag2", "Name of calibration file set (required; tag2)", value = "Group A")),
								column(width = 5,dateInput("Measadd_cal_date1", "Date start", value = NULL, min = NULL,max = NULL, format = "yyyy-mm-dd", startview = "month",weekstart = 0, language = "en")),
								column(width = 5,textInput("Measadd_cal_time1", "Time start (HH:MM:SS)", value = "12:00:00")),
								column(width = 5,dateInput("Measadd_cal_date2", "Date end", value="2018-01-01", min = NULL,max = NULL, format = "yyyy-mm-dd", startview = "month",weekstart = 0, language = "en")),
								column(width = 5,textInput("Measadd_cal_time2", "Time end (HH:MM:SS)", value = "12:00:00"))							
							)
						),
						conditionalPanel(
							condition = "input.Measadd_type == 'spiked'", 
							HTML('<hr noshade="noshade" />'),
							fluidRow(
								column(width = 5,textInput("Measadd_spiked_tag2", "ID of reference file (tag2)", value = "FALSE")),
								column(width = 5,dateInput("Measadd_recov_date", "Date start", value = NULL, min = NULL,max = NULL, format = "yyyy-mm-dd", startview = "month",weekstart = 0, language = "en")),
								column(width = 5,textInput("Measadd_recov_time", "Time start (HH:MM:SS)", value = "12:00:00"))
							)
						),						
						HTML('<hr noshade="noshade" />'),
						div(style = widget_style,
							fileInput("Measadd_path", "Upload centroided .mzXML file:", multiple = FALSE, accept = c(".mzXML",".raw")),
							bsPopover("Measadd_path", 
								title = "WARNING",
								content = "Files must be centroided. Check Package enviPick whether peaks can be picked properly from your files. Tested with Orbitrap files only.
								If reload fails, press cancel in the file upload window first", 
								placement = "right", trigger = "hover"),
							textOutput("had_meas_added")		
						)
					),
					# DELETE FILE ##################################################
					bsCollapsePanel("Delete LC-MS file", 		
						tags$h5("Delete file by its unique ID from the below file table"),
						textInput("Measdel_ID", "ID:", value = "123"),
						bsButton("Measdel","Remove",style="primary")		
					),				
					# MODIFY FILE ##################################################
					bsCollapsePanel("Modify specifications for a single file", 
						fluidRow(
							column(width = 4, helpText("Load settings of a file into below mask by its ID, modify and then save the new settings into the main table. Modifications make a full recalculation default.")),
							column(width = 3, textInput("Modif_ID", "ID:", value = "123"), bsButton("Modif_load","Load",style="primary"))
						),
						HTML('<hr noshade="noshade" />'),
						fluidRow(
							column(width = 5, textInput("Modif_name", "Name:", value = "Sample 1")),
							column(width = 5, selectInput("Modif_type", "Type:", choices = c("sample", "blank", "calibration", "spiked"))),
							column(width = 5, selectInput("Modif_include","Include in workflow?",choices = c("TRUE","FALSE"),selected="TRUE")),
							column(width = 5, selectInput("Modif_mode", "Choose ionization mode:", choices = c("positive", "negative")))	
						),
						conditionalPanel(
							condition = "input.Modif_type == 'sample' | input.Modif_type == 'blank'", 						
							HTML('<hr noshade="noshade" />'),
							fluidRow(
								column(width = 5,textInput("Modif_place", "Place:", value = "Rhine")),										
								column(width = 5,dateInput("Modif_date", "Date", value = NULL, min = NULL,max = NULL, format = "yyyy-mm-dd", startview = "month",weekstart = 0, language = "en")),	
								column(width = 5,textInput("Modif_time", "Time:(HH:MM:SS)", value = "12:00:00")),
								column(width = 5,textInput("Modif_tag3", "Replicate group (tag3)", value = "FALSE")),
								column(width = 5,textInput("Modif_ID2", "Custom ID", value = "FALSE")),
								column(width = 5, selectInput("Modif_profiled","Use for profiling (if Settings/Profiling/Omit adjusted to do so)?",choices = c("TRUE","FALSE"),selected="TRUE"))								
							)
						),
						conditionalPanel(
							condition = "input.Modif_type == 'calibration'", 
							HTML('<hr noshade="noshade" />'),
							fluidRow(
								column(width = 5,textInput("Modif_tag1", "Concentration (no units; tag1)", value = "0")),
								column(width = 5,textInput("Modif_tag2", "Calibration group (tag2)", value = "FALSE")),
								column(width = 5,dateInput("Modif_cal_date1", "Date start", value = NULL, min = NULL,max = NULL, format = "yyyy-mm-dd", startview = "month",weekstart = 0, language = "en")),
								column(width = 5,textInput("Modif_cal_time1", "Time start (HH:MM:SS)", value = "12:00:00")),
								column(width = 5,dateInput("Modif_cal_date2", "Date end", value = NULL, min = NULL,max = NULL, format = "yyyy-mm-dd", startview = "month",weekstart = 0, language = "en")),
								column(width = 5,textInput("Modif_cal_time2", "Time end (HH:MM:SS)", value = "12:00:00"))							
							)
						),
						conditionalPanel(
							condition = "input.Modif_type == 'spiked'", 
							HTML('<hr noshade="noshade" />'),
							fluidRow(
								column(width = 5,textInput("Modif_spiked_tag2", "ID of file to subtract from (tag2)", value = "FALSE")),
								column(width = 5,dateInput("Modif_recov_date", "Date", value = NULL, min = NULL,max = NULL, format = "yyyy-mm-dd", startview = "month",weekstart = 0, language = "en")),	
								column(width = 5,textInput("Modif_recov_time", "Time:(HH:MM:SS)", value = "12:00:00"))
							)
						),						
						HTML('<hr noshade="noshade" />'),
						bsButton("Modif_export","Save",style="primary")
					),
					# MODIFY CALIBRATION GROUP ######################################
					bsCollapsePanel("Modify, copy or delete a calibration group", 
						fluidRow(
							column(width = 4, helpText("Load a calibration group to be modified, deleted or copied below. To do so, select the ionization mode, insert the calibration group name (tag2 in the file table) and press Load.")),
							column(width = 3, selectInput("Modif_cal_mode", "Choose ionization mode:", choices = c("positive", "negative"))),
							column(width = 3, textInput("Modif_cal_group", "Calibration group (tag2):", value = "existinggroup"), bsButton("Load_cal","Load",style="primary"))
						),
						HTML('<hr noshade="noshade" />'),
						htmlOutput('Modif_cal_text_load'),
						HTML('<hr noshade="noshade" />'),
						helpText("Modify the specifications for all files of the above loaded calibration group and press Save to make the changes permanent"),
						fluidRow(
							column(width = 5,dateInput("Modif_calgroup_date1", "Date start", value = NULL, min = NULL,max = NULL, format = "yyyy-mm-dd", startview = "month",weekstart = 0, language = "en")),
							column(width = 5,textInput("Modif_calgroup_time1", "Time start (HH:MM:SS)", value = "12:00:00")),
							column(width = 5,dateInput("Modif_calgroup_date2", "Date end", value = NULL, min = NULL,max = NULL, format = "yyyy-mm-dd", startview = "month",weekstart = 0, language = "en")),
							column(width = 5,textInput("Modif_calgroup_time2", "Time end (HH:MM:SS)", value = "12:00:00"))
						),
						bsButton("Change_cal","Save",style="primary"),					
						HTML('<hr noshade="noshade" />'),
						helpText("Copy the calibration group to create a new calibration group with a different name. This will directly incorporate the above (changed) specifications for the new group."),
						textInput("Copy_cal_group", "New calibration group name (no underscores permitted!):", value = "newgroup"),
						bsButton("Copy_cal","Copy",style="primary"),	
						HTML('<hr noshade="noshade" />'),
						helpText("Delete all the files and any associated contents of the loaded calibration group."),
						bsButton("Del_cal","Delete",style="primary"),
						bsModal("Del_cal_confirm", "Sure about deleting the selected calibration file set?", "Del_cal", size = "small",
							bsButton("yes_delete_cal","Yes",style="warning")
						)
					),		
					# BATCH UPLOAD ##################################################
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
						bsButton("Import_file_folder","Import",style="primary"),
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
						bsButton("Import_project","Import",style="primary"),		
						HTML('<hr noshade="noshade" />'),
						textOutput("had_import_project")	
					),
					bsCollapsePanel("File overview", 
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
						HTML('<font size="5"> + </font><font size="3"> Last selected file IDs, positive ionization:</font>'),
						htmlOutput("info_files_pos_samp"),htmlOutput("info_files_pos_blind"),htmlOutput("info_files_pos_cal"),htmlOutput("info_files_pos_calgroup"),htmlOutput("info_files_pos_spiked"),
						HTML('<hr noshade="noshade" />'),		
						HTML('<font size="5"> - </font><font size="3"> Last selected file IDs, negative ionization:</font>'),
						htmlOutput("info_files_neg_samp"),htmlOutput("info_files_neg_blind"),htmlOutput("info_files_neg_cal"),htmlOutput("info_files_neg_calgroup"),htmlOutput("info_files_neg_spiked")
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
				HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/inputs/IS.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; Check help for details</a></p>'),	
				bsCollapse(multiple = FALSE, open = "IS_comp_open", id = "IS_comp",
					bsCollapsePanel("Add or modify an internal standard compound", 
						fluidRow(
							column(width = 5, helpText("To add a new internal standard, fill out the below form and press Add") ),
							column(width = 2, offset = 0.3, bsButton("AddIS","Add",style="primary"))
						),
						HTML('<hr noshade="noshade" />'),
						fluidRow(
							column(width = 5, helpText("To modify a compound insert its ID into the field to the right and press Load,
								which updates the below form. Modifications in the form can then be saved with the Modify button.
								If the compound ID in the form is changed, a new entry to the compound table with this ID will be made.
								Beware: if another compound with this new ID exists, it will be replaced in the table.") ),
							column(width = 4, 
								textInput("ISmodif_ID", "Compound ID:", value = "123_XYZ"),
								bsButton("LoadIS","Load",style="primary")),
							column(width = 3, offset = 0.3, bsButton("ModifIS","Modify",style="primary"))
						),
						HTML('<hr style="border: .8px solid " />'),
						fluidRow(
							column(width = 2, textInput("ISadd_ID", "Unique ID:", value = "123_XYZ")),  
							column(width = 8, helpText("Must be unique, numbers and characters permitted; no empty spaces or special characters or underscores permitted.") )							
						),
						HTML('<hr style="border: .8px solid " />'),
						fluidRow(         
							column(width = 3, textInput("ISadd_name", "Name:", value = "CompoundX")),
							column(width = 3, selectInput("ISadd_charge", label="Ionization mode", choices=c("positive","negative"), selected = "positive", multiple = FALSE))									
						),
						HTML('<hr style="border: .8px solid " />'),
						fluidRow(
							column(width = 3, textInput("ISadd_formula", "Formula:", value = "C10H12O10")),
							column(width = 8, helpText("Must consist of an upper case letter, possibly followed by lower case letters; to refer to individual isotopes (e.g., from isotope 
									labelling of a molecule, e.g., N5 vs. [15]N2N3), square brackets may precede the capital letter. Any other symbols which may be part of a chemical formula (e.g., charges (+), 
									dashes, asterisks, ...) are not permitted.") )										
						),
						HTML('<hr noshade="noshade" />'),
						fluidRow(
							column(width = 3, textInput("ISadd_RT", "Retention time (RT) [min]:", value = "7.52")),    
							column(width = 3, checkboxInput("ISadd_RTtol_use", "Use specific +/- RT tolerance [min]:", FALSE)),                    
							column(width = 3, textInput("ISadd_RTtol","", value = "2"))
						),
						HTML('<hr noshade="noshade" />'),
						fluidRow(						
							column(width = 4, 
								selectInput("ISadd_add", label="Main adduct:", choices= "FALSE", selected = "FALSE", multiple = FALSE),
								checkboxInput("ISadd_rest_adduct", "Restrict screening to main adduct?", FALSE)),
							column(width = 8, helpText("A compound-specific adduct can be defined here; general adducts to be considered for all compounds can be defined in the Settings/Adducts tab. 	
														Unless the below restriction, the compound-specific adduct is used alongside the general ones. 
														Seperate compound entries have to be made when including more than one compound-specific adduct."))							
						),
						HTML('<hr noshade="noshade" />'),
						fluidRow(	
							column(width = 4, checkboxInput("ISadd_use_recal", "To be used for m/z recalibration?", TRUE)),
							column(width = 4, checkboxInput("ISadd_use_screen", "To be used for screening?", TRUE))
						),
						HTML('<hr noshade="noshade" />'),
						fluidRow(
							column(width = 4, textInput("ISadd_remark", "remark", value = "none")),
							column(width = 4, textInput("ISadd_tag1", "tag1", value = "none")),
							column(width = 4, textInput("ISadd_tag2", "tag2", value = "none")),
							column(width = 4, textInput("ISadd_tag3", "tag3", value = "none"))
						),
						HTML('<hr noshade="noshade" />'),
						#div(style = widget_style,							
						#	checkboxInput("ISadd_date", "Restrict to temporal range?", FALSE),	
						#	dateRangeInput("ISadd_date_range", label="", start = NULL, end = NULL,
						#		min = NULL, max = NULL, format = "yyyy-mm-dd",
						#		startview = "month", weekstart = 0, language = "en",
						#		separator = " to ")
						#),
						fluidRow(
							column(width = 12, tags$h4("Quantification settings")),
							column(width = 4,
								helpText("Used adduct:"),
								selectInput("IS_quant_add", label=NULL, choices= "FALSE", selected = "FALSE", multiple = FALSE)
							),		
							column(width = 5,
								helpText("Used centroid peak (sorted by mass, 1=monoisotopic):"),
								numericInput("IS_quant_peak", label=NULL, 1)
							),
							column(width = 8,
								fluidRow(
									column(width = 5,
										HTML('<span class="help-block">Lower log<sub>10</sub> intensity bound of centroid peak used for quantification (.-separated):</span>'),
										textInput("Lower_intensity_bound", label=NULL, value = "0")),						
									column(width = 5,
										HTML('<span class="help-block">Upper log<sub>10</sub> intensity bound of centroid peak used for quantification (.-separated):</span>'),
										textInput("Upper_intensity_bound", label=NULL, value = "Inf"))
								)
							),
							column(width = 5,
								helpText("To rank several centroid peaks to quantify with, use:"),
								selectInput("IS_quant_rule", label=NULL, choices=c("most intense peak","closest RT","closest m/z"), selected = "most intense peak",width='400px'))							
						)
					),
					bsCollapsePanel("Remove internal standard compound", 					
						helpText("To delete a compound from the list, type in its ID and press Delete"),
						textInput("ISdelete_ID", "ID for deletion:", value = "123_XYZ"),          
						bsButton("DeleteIS","Delete",style="primary")					
					),
					bsCollapsePanel("Import / export internal standard compound list", 		
						helpText("Import IS compund list.txt file, e.g. from the dataframes folder of another project (where it can be found as IS.txt):"),					
						fileInput("ISlist_path", NULL, multiple = FALSE, accept = c(".txt")),
						checkboxInput("ISlist_save_copy", "Save a copy of the current IS compound table below if it is to be replaced by the import?", TRUE),
						bsPopover("ISlist_save_copy", 
							title = "Safety copy",
							content = "The copy can be found in the dataframes folder of your project, named as IS_date_time. You may reload this backup later to undo a compound import.", 
							placement = "right", trigger = "hover"),
						HTML('<hr noshade="noshade" />'),
						helpText("Export and save the below IS table as .txt file. The latter can again be reloaded after modifications, using the above import."),
						shinySaveButton(id="download_IS", label="Save", title="Save below IS table", filetype=list(text='txt'), buttonType = "default", class = NULL)	
					),
					bsCollapsePanel("Modify in external editor", 					
						HTML('
							<p><font>
								The below compound table can be exported as .txt file and modified in external text editors, OpenOffice Calc or Excel and then again imported. 
								Export and import MUST be done with the above Import / export functions to check for correctness of the tables. The .txt files are tab-separated when exported and must be so for import.
								Mind the character encoding (e.g., ANSI) of .txt files when modifying and ensure your imported compound table does not contain empty rows, especially at the table end.
							</font></p>
							<p style="background-color:darkred"; align="center"> <font color="#FFFFFF"> 
								Any external modifications must strictly adhere to the following rules:
							</font></p> 
							<ol>
								<li>Input format: text file (.txt), tab delimited.</li>
								<li>Column headers: DO NOT MODIFY.</li>
								<li>No empty columns. If you are not sure what to fill in, simply stay with what is set in the exported file.</li>
								<li>Absolutely NO duplicated IDs or IDs containing underscores.</li>
								<li>Do not delete columns; their number, order and content types are all fixed for good reasons.</li>
								<li>No uncompleted entries per compounds,i.e., row.</li>
								<li>Numeric entries with decimal points: dot-separated.</li>	
							</ol>
							<p style="background-color:darkgreen"; align="center"> <font color="#FFFFFF"> 
								Help for column contents of the internal standard table:
							</font></p> 
							<ol>
								<li><b>ID</b>: must be unique, numbers and characters permitted; no empty spaces or special characters. Absolutely NO duplicate IDs and NO IDs with underscores permitted.</li>
								<li><b>Name</b>: (= compound name) no special signs permitted. Use big and small letters, numbers, underscores, hyphen, brackets and empty spaces - and absolutely nothing else.</li>
								<li><b>Formula</b>: molecular formula of compounds: must consist of an upper case letter, possibly followed by lower case letters; to refer to individual isotopes (e.g., from isotope 
									labelling of a molecule, e.g., N5 vs. [15]N2N3), square brackets may precede the capital letter. Any other symbols which may be part of a chemical formula (e.g., charges (+), 
									dashes, asterisks, ...) are not permitted. </li>
								<li><b>RT</b>: retention time, IN MINUTES. Must be specified. The value will be converted to seconds during workflow usage.</li>
								<li><b>RT_tolerance</b>: FALSE or a compound-specific retention time tolerance given IN MINUTES. Overwrites the one set as standard value in tab Settings/Screening/IS. 
									The RT_tolerance will be automatically converted to seconds during workflow usage.</li>
								<li><b>main_adduct</b>: FALSE or name of a special adduct to be used for this compound entry. Valid adduct names can be found in tab Settings/Adduct.</li>								
								<li><b>ion_mod</b>e: positive or negative. If a compound is to be screened in both modes, two entries (table rows) are required.</li>
								<li><b>use_for_recalibration</b>: TRUE or FALSE.</li>
								<li><b>use_for_screening</b>: TRUE or FALSE.</li>
								<li><b>restrict_adduct</b>: TRUE or FALSE. 
									TRUE: only use the <b>main_adduct</b> (must then be specified) for this compound and ignore the ones specified in tab Settings/Adduct? 
									FALSE: the adducts specified in tab Settings/Adduct and the one set in column <b>main_adduct</b> (if specified) are all considered for this compound.</li>
								<li><b>Remark</b>: Character string for your remark on this compound (no tabs.</li>
								<li><b>tag1</b>: Character string for specification of the compound (no tabs), e.g. pharmaceutical.</li>	
								<li><b>tag2</b>: Character string for further specifications (no tabs).</li>
								<li><b>tag3</b>: Character string for further specifications (no tabs).</li>
								<li>from: ignore.</li>
								<li>to: ignore.</li>
								<li><b>Lower_intensity_bound</b>: peaks with log<sub>10</sub> intensities below this bound are not used for quantification. Set to 0 to omit; can be interactively set for each IS compound after screening, in the Results tab.</li>	
								<li><b>Upper_intensity_bound</b>: peaks with log<sub>10</sub> intensities above this bound are not used for quantification. Set to Inf to omit; can be interactively set for each IS compound after screening, in the Results tab.</li>
								<li><b>Quant_adduct</b>: ESI adduct to be used for quantification purposes. Valid adduct names can be found among those selected in tab Settings/Adduct or must be in agreement with <b>main_adduct</b>. </li>
								<li><b>Quant_peak</b>: Integer number refering to the centroid peak (ordered by increasing mass; 1 = monoisotopic) of the adduct to be used for quantification purposes. </li>
								<li><b>Quant_rule</b>: Any ONE of exactly these choices "most intense peak","closest RT" or "closest m/z". Used for quantification if several peaks are available. 
									The one quantified concentration after applying the rule will be listed first in the quantification tables for targets linked to this internal standard.</li>
							</ol>
						')		
					)
				),
                HTML('<hr noshade="noshade" />'),
                DT::dataTableOutput("IS")
            ),
            tabPanel("Targets", 
				HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/inputs/targets.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; Check help for details</a></p>'),	
				bsCollapse(multiple = FALSE, open = "target_comp_open", id = "target_comp",
					bsCollapsePanel("Add or modify a target compound", 
						fluidRow(
							column(width = 5, helpText("To add a new target, fill out the below form and press Add.") ),
							column(width = 2, offset = 0.3, bsButton("Addtargets","Add",style="primary"))
						),
						HTML('<hr noshade="noshade" />'),
						fluidRow(
							column(width = 5, helpText("To modify a compound insert its ID into the field to the right and press Load,
								which updates the below form. Modifications in the form can then be saved with the Modify button.
								If the compound ID in the form is changed, a new entry to the compound table with this ID will be made.
								Beware: if another compound with this new ID already exists, it will be replaced in the table.") ),
							column(width = 4, 
								textInput("targetmodif_ID", "Compound ID:", value = "123_XYZ"),
								bsButton("Loadtarget","Load",style="primary")),
							column(width = 3, offset = 0.3, bsButton("Modiftarget","Modify",style="primary"))
						),
						HTML('<hr style="border: .8px solid " />'),
						fluidRow(
							column(width = 2, textInput("targetsadd_ID", "Unique ID:", value = "123_XYZ")),          
							column(width = 8, helpText("Must be unique, numbers and characters permitted; no empty spaces or special characters or underscores permitted."))	
						),
						HTML('<hr style="border: .8px solid " />'),
						fluidRow(         
							column(width = 3, textInput("targetsadd_name", "Name:", value = "CompoundX")),
							column(width = 3, selectInput("targetsadd_charge", label="Ionization mode", choices=c("positive","negative"), selected = "positive", multiple = FALSE))									
						),						
						HTML('<hr style="border: .8px solid " />'),
						fluidRow(
							column(width = 3, textInput("targetsadd_formula", "Formula:", value = "C10H12O10")),
							column(width = 8, helpText("Must consist of an upper case letter, possibly followed by lower case letters; to refer to individual isotopes (e.g., from isotope 
									labelling of a molecule, e.g., N5 vs. [15]N2N3), square brackets may precede the capital letter. Any other symbols which may be part of a chemical formula (e.g., charges (+), 
									dashes, asterisks, ...) are not permitted."))		
						),						
						HTML('<hr noshade="noshade" />'),
						fluidRow(
							column(width = 3, textInput("targetsadd_RT", "Retention time (RT) [min]:", value = "7.52")),    
							column(width = 3, checkboxInput("targetsadd_RTtol_use", "Use specific +/- RT tolerance [min]:", FALSE)),                    
							column(width = 3, textInput("targetsadd_RTtol","", value = "2"))
						),
						HTML('<hr noshade="noshade" />'),
						fluidRow(						
							column(width = 4, 
								selectInput("targetsadd_add", label="Main adduct:", choices= "FALSE", selected = "FALSE", multiple = FALSE),
								checkboxInput("targetsadd_rest_adduct", "Restrict screening to main adduct?", FALSE)),
							column(width = 8, helpText("A compound-specific adduct can be defined here; general adducts to be considered for all compounds can be defined in the Settings/Adducts tab. 	
														Unless the below restriction, the compound-specific adduct is used alongside the general ones. 
														Seperate compound entries have to be made when including more than one compound-specific adduct."))							
						),
						HTML('<hr noshade="noshade" />'),
						fluidRow(	
							column(width = 4, checkboxInput("targetsadd_use_recal", "To be used for m/z recalibration?", TRUE)),
							column(width = 4, checkboxInput("targetsadd_use_screen", "To be used for screening?", TRUE))
						),
						HTML('<hr noshade="noshade" />'),
						fluidRow(
							column(width = 4, textInput("targetsadd_remark", "remark", value = "none")),
							column(width = 4, textInput("targetsadd_tag1", "tag1", value = "none")),
							column(width = 4, textInput("targetsadd_tag2", "tag2", value = "none")),
							column(width = 4, textInput("targetsadd_tag3", "tag3", value = "none"))
						),
						HTML('<hr noshade="noshade" />'),
						#div(style = widget_style,							
						#	checkboxInput("targetsadd_date", "Restrict to temporal range?", FALSE),	
						#	dateRangeInput("targetsadd_date_range", label="", start = NULL, end = NULL,
						#		min = NULL, max = NULL, format = "yyyy-mm-dd",
						#		startview = "month", weekstart = 0, language = "en",
						#		separator = " to ")
						#),
						fluidRow(
							column(width = 12, tags$h4("Quantification settings")),
							column(width = 4,
								helpText("Used adduct:"),
								selectInput("target_quant_add", label=NULL, choices= "FALSE", selected = "FALSE", multiple = FALSE)
							),		
							column(width = 5,
								helpText("Used centroid peak (sorted by mass, 1=monoisotopic):"),
								numericInput("target_quant_peak", label=NULL, 1)
							),
							column(width = 5,
								helpText("ID of internal standard used for calibration & quantification for this target:"),
								textInput("target_quant_ISID", label=NULL, value = "FALSE")
							),							
							column(width = 8,
								fluidRow(
									column(width = 5,
										helpText("First concentration warn level (.-separated):"),
										textInput("warn_1", label=NULL, value = "FALSE")),						
									column(width = 5,
										helpText("Second concentration warn level (.-separated):"),
										textInput("warn_2", label=NULL, value = "FALSE"))
								)
							),
							column(width = 5,
								helpText("To rank several centroid peaks to quantify with, use:"),
								selectInput("target_quant_rule", label=NULL, choices=c("most intense peak","closest RT","closest m/z"), selected = "most intense peak",width='400px'))							
						)
					),
					bsCollapsePanel("Remove target compound", 
						helpText("To delete a compound from the list, type in its ID and press Delete"),
						textInput("targetsdelete_ID", "ID for deletion:", value = "123_XYZ"),          
						bsButton("Deletetargets","Delete",style="primary")					
					),
					bsCollapsePanel("Import / export target compound list", 		
						helpText("Import target compound list.txt file, e.g. from the dataframes folder of another project (where it can be found as targets.txt):"),					
						fileInput("targetlist_path", NULL, multiple = FALSE, accept = c(".txt")),
						checkboxInput("targetlist_save_copy", "Save a copy of the current target compound table below if it is to be replaced by the import?", TRUE),
						bsPopover("targetlist_save_copy", 
							title = "Safety copy",
							content = "The copy can be found in the dataframes folder of your project, named as targets_date_time. You may reload this backup later to undo a compound import.", 
							placement = "right", trigger = "hover"),
						HTML('<hr noshade="noshade" />'),
						helpText("Export and save the below target compound table as .txt file. The latter can again be reloaded after modifications, using the above import."),
						shinySaveButton(id="download_target", label="Save", title="Save below target compound table", filetype=list(text='txt'), buttonType = "default", class = NULL)	
					),
					bsCollapsePanel("Modify in external editor", 					
						HTML('
							<p><font>
								The below compound table can be exported as .txt file and modified in external text editors, OpenOffice Calc or Excel and then again imported. 
								Export and import MUST be done with the above Import / export functions to check for correctness of the tables. The .txt files are tab-separated when exported and must be so for import.
								Mind the character encoding (e.g., ANSI) of .txt files when modifying and ensure your imported compound table does not contain empty rows, especially at the table end.
							</font></p>
							<p style="background-color:darkred"; align="center"> <font color="#FFFFFF"> 
								Any external modifications must strictly adhere to the following rules:
							</font></p> 
							<ol>
								<li>Input format: text file (.txt), tab delimited.</li>
								<li>Column headers: DO NOT MODIFY.</li>
								<li>No empty columns. If you are not sure what to fill in, simply stay with what is set in the exported file.</li>
								<li>Absolutely NO duplicated IDs or IDs containing underscores.</li>
								<li>Do not delete columns; their number, order and content types are all fixed for good reasons.</li>
								<li>No uncompleted entries per compounds,i.e., row.</li>
								<li>Numeric entries with decimal points: dot-separated.</li>
								<li>If you want to make links to the internal standard table via column <b>ID_internal_standard</b>, the internal standard
								must exists in the according table. Thus, better start setting up the internal standard table first.</li>								
							</ol>
							<p style="background-color:darkgreen"; align="center"> <font color="#FFFFFF"> 
								Help for target table column contents:
							</font></p> 
							<ol>
								<li><b>ID</b>: must be unique, numbers and characters permitted; no empty spaces or special characters. Absolutely NO duplicate IDs and NO IDs with underscores permitted.</li>
								<li><b>Name</b>: (= compound name) no special signs permitted. Use big and small letters, numbers, underscores, hyphen, brackets and empty spaces - and absolutely nothing else.</li>
								<li><b>Formula</b>: molecular formula of compounds: must consist of an upper case letter, possibly followed by lower case letters; to refer to individual isotopes (e.g., from isotope 
									labelling of a molecule, e.g., N5 vs. [15]N2N3), square brackets may precede the capital letter. Any other symbols which may be part of a chemical formula (e.g., charges (+), 
									dashes, asterisks, ...) are not permitted. </li>
								<li><b>RT</b>: retention time, IN MINUTES. Must be specified. The value will be converted to seconds during workflow usage.</li>
								<li><b>ID_internal_standard</b>: FALSE or a valid ID of an internal standards to be found in the CURRENT internal standard table. Only targets linked with this ID will be considered for quantification.</li>
								<li><b>RT_tolerance</b>: FALSE or a compound-specific retention time tolerance given IN MINUTES. Overwrites the one set as standard value in tab Settings/Screening/IS. 
									The RT_tolerance will be automatically converted to seconds during workflow usage.</li>
								<li><b>main_adduct</b>: FALSE or name of a special adduct to be used for this compound entry. Valid adduct names can be found in tab Settings/Adduct.</li>								
								<li><b>ion_mod</b>e: positive or negative. If a compound is to be screened in both modes, two entries (table rows) are required.</li>
								<li><b>use_for_recalibration</b>: TRUE or FALSE.</li>
								<li><b>use_for_screening</b>: TRUE or FALSE.</li>
								<li><b>restrict_adduct</b>: TRUE or FALSE. 
									TRUE: only use the <b>main_adduct</b> (must then be specified) for this compound and ignore the ones specified in tab Settings/Adduct? 
									FALSE: the adducts specified in tab Settings/Adduct and the one set in column <b>main_adduct</b> (if specified) are all considered for this compound.</li>
								<li><b>Remark</b>: Character string for your remark on this compound (no tabs.</li>
								<li><b>tag1</b>: Character string for specification of the compound (no tabs), e.g. pharmaceutical.</li>	
								<li><b>tag2</b>: Character string for further specifications (no tabs).</li>
								<li><b>tag3</b>: Character string for further specifications (no tabs).</li>
								<li>from: ignore.</li>
								<li>to: ignore.</li>
								<li><b>warn_1</b>: First concentration warn level, without units. Quantified targets above this level are marked orange in the quantification table. </li>
								<li><b>warn_2</b>: Second concentration warn level, without units. Quantified targets above this level are marked red in the quantification table. </li>
								<li><b>Quant_adduct</b>: ESI adduct to be used for quantification purposes. Valid adduct names can be found among those selected in tab Settings/Adduct or must be in agreement with <b>main_adduct</b>. </li>
								<li><b>Quant_peak</b>: Integer number refering to the centroid peak (ordered by increasing mass; 1 = monoisotopic) of the adduct to be used for quantification purposes. </li>
								<li><b>Quant_rule</b>: Any ONE of exactly these choices "most intense peak","closest RT" or "closest m/z". Used for quantification if several peaks are available. 
									The one quantified concentration after applying the rule will be listed first in the quantification tables.</li>
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
			tags$h5("Apply workflow settings to project?"), 
			bsButton("saveflow","Apply",style="warning"),
			bsAlert("alert_1"),
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
					radioButtons("norm", "Include?", c("yes"="yes","no"="no")),
				#HTML('<p style="background-color:darkgreen"; align="center"> <font color="#FFFFFF"> RT alignment </font></p> '),
				#radioButtons("RTalign", "Include?", c("yes"="yes","no"="no")),  
				HTML('<p style="background-color:darkblue"; align="center"> <font color="#FFFFFF"> Blank / blind peak detection </font></p> '),				
					fluidRow(
						column(width = 2, radioButtons("blind", "Detect?", c("yes"="yes","no"="no")) ),
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
						column(width = 2, radioButtons("profiling", "Include? ", c("yes"="yes","no"="no")) ),
						column(width = 10, offset = 0.3,
							tags$p(align="justify","An intensity descent assorts peaks into profiles. Uses a fixed retention time tolerance window and
							an adaptive mass tolerance window. Only peaks of files marked for profiling are used (column profiled in the files table / File Tab).")
						)
					),					
				HTML('<p style="background-color:darkgreen"; align="center"> <font color="#FFFFFF"> LOD interpolation </font></p>'),
					fluidRow(
						column(width = 2, radioButtons("LOD", "Include? ", c("yes"="yes","no"="no")) ),
						column(width = 10, offset = 0.3,
							tags$p(align="justify","For each measurement, estimates a RT-dependent intensity threshold below which peaks are not expected to get picked. 
							Can be used for the below compound screening.")
						)
					),
				HTML('<p style="background-color:darkgreen"; align="center"> <font color="#FFFFFF"> Compound screening </font></p> '),
					fluidRow(
						column(width = 3, 
							radioButtons("IS_screen", "Screen internal standards?", c("yes"="yes","no"="no"))
						),
						column(width = 3, 
							radioButtons("target_screen", "Screen target compounds?", c("yes"="yes","no"="no"))
						),						
						column(width = 6, offset = 0.3,
							tags$p(align="justify","Uses the LOD thresholds estimated in the above step. If the LOD interpolation is not run, a fixed intensity threshold as specified in the
							Settings/Screening tabs is used. All peaks of files marked for inclusion and profiling are used (column profiled in the files table / File Tab) - unless peaks are
							removed explicitly in the early blind detection (see blue workflow steps).")
						)
					),	
				HTML('<hr noshade="noshade" />'),
				HTML('<h1 align="center"> &#x21e9; </h1> '),  					
				# block 4 ######################################################					
				HTML('<p style="background-color:orange"; align="center"> <font color="#FFFFFF"> Calibration </font></p> '),
					fluidRow(
						column(width = 2, radioButtons("calibration", "Include? ", c("yes"="yes","no"="no")) ),
						column(width = 10, offset = 0.3,
							tags$p(align="justify","This step screens for calibration sets of target and internal standard compound peaks, using
							the provided calibration files. The sets can be used in the Calibration tab to establish specific calibration models (curves)
							for quantification. If selected, the extraction of these calibration peaks will be affected by the above mass recalibration, 
							replicate intersection, blind subtraction and LOD interpolation steps. Once you have established the desired calibration models,
							you can deselect this step, until new calibration files have been added."),
							HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/quantification.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="right">&#8594; More info.</a></p>')	
						)
					),
				HTML('<p style="background-color:orange"; align="center"> <font color="#FFFFFF"> Quantification </font></p> '),
					fluidRow(
						column(width = 2, radioButtons("quantification", "Include? ", c("yes"="yes","no"="no")) ),
						column(width = 10, offset = 0.3,
							tags$p(align="justify","Based on the calibration models from the above step, an estimation of target compound concentrations
							from their intensity ratios to their individual internal standard compounds is derived."),
							HTML('<a href="http://www.looscomputing.ch/eng/enviMass/topics/quantification.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="right">&#8594; More info.</a>')	
						)
					),
				HTML('<p style="background-color:orange"; align="center"> <font color="#FFFFFF"> Recovery </font></p> '),
					fluidRow(
						column(width = 2, radioButtons("recovery", "Include? ", c("yes"="yes","no"="no")) ),
						column(width = 10, offset = 0.3,
							tags$p(align="justify","Following the above calibration and quantification steps, a concentration recovery of spiked target compounds is calculated.
							Requires upload of spiked files."),
							HTML('<a href="http://www.looscomputing.ch/eng/enviMass/topics/quantification.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="right">&#8594; More info.</a>')	
						)
					),
				HTML('<hr noshade="noshade" />'),
				HTML('<h1 align="center"> &#x21e9; </h1> '),  		
				# block X ######################################################					
				HTML('<p style="background-color:black"; align="center"> <font color="#FFFFFF"> EIC correlation </font></p> '),
					fluidRow(
						#column(width = 2, radioButtons("EIC_correlation", "Include?", c("yes"="yes","no"="no"))),
						column(width = 10, offset = 0.3,tags$p(align="justify","Under construction"))
					),	
				HTML('<p style="background-color:black"; align="center"> <font color="#FFFFFF"> Isotopologue grouping </font></p> '),
					fluidRow(
						#column(width = 2, radioButtons("isotopologues", "Include?", c("yes"="yes","no"="no"))),
						column(width = 10, offset = 0.3,tags$p(align="justify","Under construction"))
					),	
				HTML('<p style="background-color:black"; align="center"> <font color="#FFFFFF"> Adduct grouping </font></p> '),
					fluidRow(
						#column(width = 2, radioButtons("adducts", "Include?", c("yes"="yes","no"="no"))),
						column(width = 10, offset = 0.3,tags$p(align="justify","Under construction"))
					),	
				HTML('<p style="background-color:black"; align="center"> <font color="#FFFFFF"> Homologue series detection </font></p> '),
					fluidRow(
						#column(width = 2, radioButtons("homologues", "Include?", c("yes"="yes","no"="no"))),
						column(width = 10, offset = 0.3,tags$p(align="justify","Under construction"))
					),	
				HTML('<hr noshade="noshade" />'),
				HTML('<h1 align="center"> &#x21e9; </h1> '),  				
				# block 5 ######################################################
				HTML('<p style="background-color:darkred"; align="center"> <font color="#FFFFFF"> Intensity normalization using IS-profiles </font></p> '),
					fluidRow(
						column(width = 2, radioButtons("IS_normaliz", "Include? ", c("yes"="yes","no"="no")) ),
						column(width = 10, offset = 0.3,
							tags$p(align="justify","Relies on the above profile extraction and internal standard screening (green steps). Intensities of picked peaks in each measurement 
							are normalized by the median deviation all internal standards have in the measurement from their individual median profile intensity.
							This approach can replace the less reliable Median intensity normalization above (blue steps). Internal standards must have been spiked at constant concentrations.")
						)
					),					
				HTML('<p style="background-color:darkred"; align="center"> <font color="#FFFFFF"> Profile filtering </font></p> '),
					fluidRow(
						column(width = 2, 
							radioButtons("subtr", "Include?", c("yes"="yes","no"="no"),inline=TRUE)
						),
						column(width = 10, offset = 0.3,
							tags$p(align="justify","Run a profile recalculation omitting...")						
						)
					),	
					HTML('<hr noshade="noshade" />'), 
					fluidRow(
						column(width = 6, offset = 0.3,
							tags$p(align="justify","... peaks belonging to matches >= the cutoff score defined in the Settings/Screening Tab.
							The original compound screening results remain visible, but affected peaks will not be part of the final profiles.")
						),
						column(width = 3, 
							radioButtons("subtr_IS", "Subtract internal standard peaks?", c("yes"="yes","no"="no"),inline=TRUE)
						),
						column(width = 3, 
							radioButtons("subtr_target", "Subtract target compound peaks?", c("yes"="yes","no"="no"),inline=TRUE)
						)
					),	
					HTML('<hr noshade="noshade" />') ,
					fluidRow(
						column(width = 10, offset = 0.3,
							tags$p(align="justify","... peaks which have also been detected in blind/blank samples but not removed yet (cp. blind detection of blue steps)?")
						),
						column(width = 2, 
							radioButtons("subtr_blind", "Subtract?", c("yes"="yes","no"="no"),inline=TRUE)
						)
					),	
					HTML('<hr noshade="noshade" />') ,
					fluidRow(
						column(width = 10, offset = 0.3,
							tags$p(align="justify","... peaks from spiked files?")
						),
						column(width = 2, 
							radioButtons("subtr_spiked", "Subtract?", c("yes"="yes","no"="no"),inline=TRUE)
						)
					),	
				HTML('<p style="background-color:darkred"; align="center"> <font color="#FFFFFF"> Trend detection </font></p> '),
					fluidRow(
						column(width = 2, radioButtons("trendblind", "Include? ", c("yes"="yes","no"="no")) ),
						column(width = 10, offset = 0.3,
							tags$p(align="justify","Depending on settings, this can contain a separate blind detection step. Herein, intensities of blind/blank peaks are interpolated over the time series.
							This interpolation and subtraction is only applicable if the separate blind filter step is disabled (see above blue steps and the preceding red step).")
						)
					),							            
				#HTML('<p style="background-color:darkred"; align="center"> <font color="#FFFFFF"> Peak grouping (componentization) </font></p>'),
				#	checkboxInput("Comp_isotop", "Group isotopologue peaks?", TRUE),
				HTML('<p style="background-color:darkred"; align="center"> <font color="#FFFFFF"> Componentization </font></p> '),
					fluidRow(
						#column(width = 2, radioButtons("components", "Include?", c("yes"="yes","no"="no"))),
						column(width = 10, offset = 0.3,tags$p(align="justify","Under construction"))
					),
				HTML('<hr noshade="noshade" />') 
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
		  tags$h5("Apply paramter settings to project?"), 
		  bsButton("savepar","Apply",style="warning"),
		  bsAlert("alert_2"),
		  HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/inputs/parameters.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; Check help on parameters</a></p>'),	
          tabsetPanel(
            # PEAK PICKING #####################################################
            tabPanel("Peak picking",
				div(style = widget_style2,
					tags$h5("Peak definition"), 
					numericInput("peak_minpeak", "Minimum number of measurements per peak ...", 4),
					sliderInput("peak_drtsmall2", "... within a given RT window [s]", min = 1, max = 200, value = 20, step= 0.1),
					sliderInput("peak_drtfill", "Maximum RT gap length to be interpolated [s]", min = 0, max = 60, value = 10, step= 0.1),
					sliderInput("peak_drtdens2", "Peak definition - Maximum RT length of a single peak", min = 10, max = 1500, value = 120, step= 0.1),
					sliderInput("peak_minint_log10", "Minimum log10(intensity) threshold", min = 0, max = 10, value = 4, step= .1),
					numericInput("peak_SN", "Minimum Signal/Noise", 5),
					numericInput("peak_SB", "Minimum Signal/Base", 2),
					numericInput("peak_recurs", "Maximum possible number of peaks within a single EIC", 3)
				),
				div(style = widget_style2,
					tags$h5("EIC partitioning & clustering"),
					sliderInput("peak_drtgap", "Maximum retention time gap in an EIC", min = 20, max = 1500, value = 300, step= 1),
					sliderInput("peak_dmzdens", "Maximum m/z deviation of a measurement from its EIC mean [ppm]", min = 1, max = 100, value = 3.5, step= 0.1)       				
				),
				div(style = widget_style,
					tags$h5("Advanced"),
					numericInput("peak_ended", "How often can a peak detection fail to end the recursion? - peak picking", 1),
					numericInput("peak_weight", "Weight for assigning measurements to a peak - peak picking", 1),
					numericInput("peak_maxint_log10", "Upper log10(intensity) safety threshold", 6.7),					
					#sliderInput("peak_maxint_log10", "Upper log10(intensity) safety threshold", min = 0, max = 15, value = 6.7, step= .1),
					sliderInput("peak_perc_cut", "Percentage of low-intense data points to discard", min = 0, max = 100, value = 0, step= .1)
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
                  div(style = widget_style3,selectInput("resolution", label="Instrument resolution:", choices="none", selected="none", width='600px')),
				  imageOutput("plot_resolution", height="auto")	  
            ),	
            # RECALIBRATION ####################################################
            tabPanel("Recalibration",
              tags$h5("Mass recalibration"),
              div(
                style = widget_style3,
                selectInput("recal_use", "Reference compounds:", c("Internal standards","Target compounds","both"),"Internal standards",multiple=FALSE),                
                numericInput("recal_dmz", "m/z tolerance ...", 3),            
                numericInput("recal_maxdmz", "Maximum allowable m/z correction ...", 30),  				
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
				numericInput("replicate_IS_dInt", "Intensity tolerance X (log 10 scale, 10^X):", 9)
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
						numericInput("IS_drt1", "RT tolerance of peaks relative to their expected RT [s]", 30),
						numericInput("IS_drt2", "RT tolerance of peaks within an isotope pattern [s]", 50)
					),
				  	div(style = widget_style2,
						tags$h5("Mass"),
						numericInput("IS_dmz", "m/z tolerance ...", 3),                
						selectInput("IS_ppm", "... given in:", choices = c("ppm"="TRUE","absolute"="FALSE"), "TRUE")
 					),
					div(style = widget_style2,
						tags$h5("Intensity"),
						sliderInput("IS_inttol", "Intensity tolerance %", min = 0, max = 100, value = 30, step= .2),
						numericInput("IS_intcut", "Lower intensity threshold (if LOD interpolation disabled)", 5E4)                
					),
					div(style = widget_style2,
						tags$h5("Scoring"),
						numericInput("IS_w1", "Cutoff score [0,1]", 0.8),         
						HTML('<h1 align="center"> &#x21f3; </h1> '),						
						selectInput("screen_IS_cutit", "Exclude matches below cutoff score?", choices = c("TRUE"="TRUE","FALSE"="FALSE"), "FALSE"),
						HTML('<h1 align="center"> &#x21f3; </h1> '),
						selectInput("screen_IS_maxonly", "Screen only most intense isotopologue peak?", choices = c("TRUE"="TRUE","FALSE"="FALSE"), "FALSE")
					)
                ),
                tabPanel("Targets & Suspects",
					div(style = widget_style2,
						tags$h5("Retention time"),
						numericInput("tar_drt1", "RT tolerance of peaks in sample relative to their expected RT [s]", 30),
						numericInput("tar_drt2", "RT tolerance of peaks within an isotope pattern [s]", 50)
					),
				  	div(style = widget_style2,
						tags$h5("Mass"),
						numericInput("tar_dmz", "m/z tolerance ...", 3),                
						selectInput("tar_ppm", "... given in:", choices = c("ppm"="TRUE","absolute"="FALSE"), "TRUE")
 					),
					div(style = widget_style2,
						tags$h5("Intensity"),
						sliderInput("tar_inttol", "Intensity tolerance %", min = 0, max = 100, value = 30, step= .2),
						numericInput("tar_intcut", "Lower intensity threshold", 5E4)                
					),
					div(style = widget_style2,
						tags$h5("Scoring"),
						numericInput("tar_w1", "Cutoff score [0,1]", 0.8),
						HTML('<h1 align="center"> &#x21f3; </h1> '),
						selectInput("screen_target_cutit", "Exclude matches below cutoff score?", choices = c("TRUE"="TRUE","FALSE"="FALSE"), "FALSE"),
						HTML('<h1 align="center"> &#x21f3; </h1> '),
						selectInput("screen_target_maxonly", "Screen only most intense isotopologue peak?", choices = c("TRUE"="TRUE","FALSE"="FALSE"), "FALSE")					
					)
                )
              )
            ),
			# QUANTIFICATION & RECOVERY ########################################
			tabPanel("Quantification",
				HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/quantification.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; Check help for details & parameter descriptions.</a></p>'),	
				numericInput("quant_files_included", "Number of latest file to include in the quantification:", 30),
				numericInput("recov_files_included", "Number of latest spiked files to include for recovery estimation:", 30)			
			),
            # INTENSITY NORMALIZATION ##########################################
            tabPanel("Normalization",
				tags$h5("Intensity normalization based on IS-profiles"),
				sliderInput("ISnorm_percfiles", "Minimum of files covered by each IS profile (%)", min = 0, max = 100, value = 90, step= 1),
				sliderInput("ISnorm_score", "Screening threshold", min = 0, max = 1, value = 0.8, step= .01),
				numericInput("ISnorm_numbIS", "Minimum number of IS profiles", 15),
				HTML('<hr noshade="noshade" />'),
				checkboxInput("ISnorm_medblank", "Show median deviation of blank/blind profiles?", FALSE),
				checkboxInput("ISnorm_usesubblank", "Use subsampling", FALSE),
				numericInput("ISnorm_numblank", "Number of blank/blind profiles in subsample", 100),
				HTML('<hr noshade="noshade" />'),
				checkboxInput("ISnorm_medsam", "Show median deviation of sample (i.e., non-blank) profiles?", FALSE),	
				checkboxInput("ISnorm_usesubsam", "Use subsampling", FALSE),
				numericInput("ISnorm_numsam", "Number of sample profiles in subsample", 100)
			),
            # PROFILING ########################################################
            tabPanel("Profiling",
				tags$h5("Profile extraction"),
				sliderInput("prof_maxfiles", "Maximum number of newest files to be processed (by date/time) ", min = 50, max = 3000, value = 100, step= 1),
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
				numericInput("trend_threshold", "Trend vs. mean+variance intensity threshold:", 3),
				HTML('<hr noshade="noshade" />'),
				radioButtons("trend_blind", "Additional blind interpolation and subtraction per profile?", c("yes"="yes","no"="no")),
				HTML('<hr noshade="noshade" />'),
				selectInput("notrend", "Do not show global trend - instead, report it as maximum intensity above blind", choices = c("TRUE"="TRUE","FALSE"="FALSE"), "FALSE")
            ),
            # BLIND #############################################################
            tabPanel("Blind",
				tags$h5("Blind subtraction:"),
				numericInput("blind_threshold", "Intensity threshold ratio sample/blind <", 100),
				numericInput("blind_dmz", "Mass uncertainty (+/-) ...", 3), 
                selectInput("blind_ppm", "... given in:", choices = c("ppm"="TRUE","absolute"="FALSE"), "TRUE"),				
                numericInput("blind_drt", "RT tolerance [s]:", 60),       
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
					selectInput("progressBar", "Show progress bars (Windows OS only)", choices = c("TRUE","FALSE"), selected="FALSE"),
					selectInput("do_project_check", "Skip the project check before calculations?", choices = c("TRUE","FALSE"), selected="FALSE"),					
					selectInput("ignore_large_files", "Ignore .mzXML and MSlist files during check?", choices = c("TRUE","FALSE"), selected="FALSE"),					
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
        # CALIBRATION ##########################################################
        ########################################################################
        tabPanel("Calibration", 	
			HTML('<a href="http://www.looscomputing.ch/eng/enviMass/topics/quantification.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; More info.</a>'),
			tabsetPanel( 
				tabPanel("Create",
					helpText("Select the ionization mode to load the available calibration file sets into the below selection."),
					selectInput("Ion_mode_Cal", label="Ionization mode", c("none","positive","negative"), selected = ("none"), multiple = FALSE),	
					HTML('<hr noshade="noshade" />'),
					HTML('<h1 align="center"> &#x21e9; </h1> '),
					conditionalPanel(
						condition = "input.Ion_mode_Cal != 'none'", 	
						selectInput("Cal_file_set",label="Specify calibration file set",choices=c("none"),selected = "none", multiple = FALSE),
						conditionalPanel(
							condition = "input.Cal_file_set != 'none' & input.Ion_mode_Cal != 'none'", 						
							bsButton("Cal_file_set_delete", label="Delete all models for this calibration file set?", icon = icon("fire-extinguisher"))
						),
						HTML('<hr noshade="noshade" />'),
						HTML('<h1 align="center"> &#x21e9; </h1> '),
						conditionalPanel(
							condition = "input.Cal_file_set != 'none' & input.Ion_mode_Cal != 'none'", 
							helpText("Select Target compounds which are linked to an Internal standard (column 6 of the target compound table) to (resume) work on their individual calibration models below. 
							The screened compounds can also be viewed in the Results/Compound screening tab, choosing the calibration files. Target - Internal standard links can be modified either directly 
							in the target compound table (above tab Compound -> Targets) or with the below orange button Use for quantification"),
							fluidRow(
								column(3, selectInput("Cal_target_name",label="Target name",choices=c("none"),selected = "none", multiple = FALSE)),							
								column(3, selectInput("Cal_target_ID",label="Target ID",choices=c("none"),selected = "none", multiple = FALSE))
							),
							fluidRow(
								column(3, selectInput("Cal_IS_name",label="Internal standard name",choices=c("none"),selected = "none", multiple = FALSE)),
								column(3, selectInput("Cal_IS_ID",label="Internal standard ID",choices=c("none"),selected = "none", multiple = FALSE))
							),
							textOutput('number_missing_models'),
							HTML('<hr noshade="noshade" />'),			
							helpText("Proceed through entries in the target compound list:"),
							fluidRow(
								column(1, bsButton("Cal_first", label="", icon = icon("step-backward"))),
									bsTooltip("Cal_first", title="Jump to first Target", placement = "top", trigger = "hover", options = NULL),	
								column(1, bsButton("Cal_previous", label="", icon = icon("arrow-left"))),
									bsTooltip("Cal_previous", title="Jump to preceeding (or last, if beginning of table is reached) Target", placement = "top", trigger = "hover", options = NULL),
								column(1, bsButton("Cal_next", label="", icon = icon("arrow-right"))),		
									bsTooltip("Cal_next", title="Jump to next (or first, if end of table is reached) Target", placement = "top", trigger = "hover", options = NULL),	
								column(1, bsButton("Cal_next_missing", label="", icon = icon("arrow-right"), style="info")),		
									bsTooltip("Cal_next_missing", title="Jump to next Target with missing calibration model", placement = "top", trigger = "hover", options = NULL),	
								column(1, bsButton("Cal_missing", label="", icon = icon("map-pin"), style="info")),	
									bsTooltip("Cal_missing", title="Jump to first Target with missing calibration model", placement = "top", trigger = "hover", options = NULL),
								column(1, bsButton("Cal_last", label="", icon = icon("step-forward"))),
									bsTooltip("Cal_last", title="Jump to last Target", placement = "top", trigger = "hover", options = NULL)							
							),							
							HTML('<hr noshade="noshade" />'),
							HTML('<h1 align="center"> &#x21e9; </h1> '),
							bsButton("save_Cal","Save/replace model",style="success"),
							bsPopover("save_Cal", 
								title = "Save the selected calibration model for later modification or to quantify with? Overwrites any existing model.",
								content = "The model is only used for quantification if the Internal standard is linked to the selected Target compound in column 6 of the target compound table. You may modify this default link with the orange button on the right.", 
								placement = "bottom", trigger = "hover"),							
							bsButton("use_Cal","Use for quantification",style="warning"),
							bsPopover("use_Cal", 
								title = "Make entry to target compound list",
								content = "Save the selected Internal standard as the one to quantify with (= entry into column 6 of the target compound table)?", 
								placement = "bottom", trigger = "hover"),							
							bsButton("remove_Cal","Remove model",style="danger"),
							bsPopover("remove_Cal", 
								title = "Delete any existing model for the compound selection.",
								content = "Target compounds without a calibration model for their linked Internal standard  (= entry into column 6 of the target compound table) will not be quantified.", 
								placement = "bottom", trigger = "hover"),							
							bsButton("reload_Cal","Reload data",style="info"),	
							bsPopover("reload_Cal", 
								title = "Refresh below plot and table for the selected compound pair of Target and Internal standard ...",
								content = "... e.g., after you have changed the intensity bounds of the Internal standard in the tab <em> Results -> Screening& -> Ionization -> Internal standards (show all adducts) -> Characteristics for selected compound. </em>", 
								placement = "bottom", trigger = "hover"),							
							conditionalPanel(
								condition = "input.Cal_file_set != 'none' & input.Ion_mode_Cal != 'none' & input.Cal_target_ID != 'none' & input.Cal_target_name != 'none' & input.Cal_IS_ID != 'none' & input.Cal_IS_name != 'none'  ", 
								HTML('<hr noshade="noshade" />'),
								div(style = widget_style6,
									plotOutput("cal_plot", 
										dblclick = "cal_plot_dblclick",
										brush = brushOpts(
										  id = "cal_plot_brush",
										  resetOnNew = TRUE
										)
									)
								),
								HTML('<hr noshade="noshade" />'),
								fluidRow(
									column(4,selectInput("cal_model", "Select calibration model (red in above plot)", choices = c("linear","quadratic"), "linear")),
									column(4,checkboxInput("cal_model_0intercept", "Force 0-intercept?",  width = NULL)),
									column(4,checkboxInput("cal_model_weight", "Weight by inverse target intensity?",  width = NULL))
								),
								textOutput('cal_model_summary'),
								HTML('<hr noshade="noshade" />'),
								fluidRow(
									column(4,
										checkboxInput("cal_model_bound_low", "Set a lower bound for the intensity ratio?",  width = NULL),
										numericInput("cal_model_bound_low_value", "Lower bound log10:", 0)
									),
									column(4,
										checkboxInput("cal_model_bound_up", "Set an upper bound for the intensity ratio?",  width = NULL),
										numericInput("cal_model_bound_up_value", "Upper bound log10:", 20)
									)									
								),						
								HTML('<hr noshade="noshade" />'),
								helpText("Click into the below table rows to select and deselect data points for the above calibration model (red line):"),
								fluidRow(
									column(12,DT::dataTableOutput('cal_table'))
								),
								HTML('<hr noshade="noshade" />'),		
								numericInput("use_precision", "Numeric match precision", 0.01),
								bsPopover("use_precision", 
									title = "Precision with which intensity ratios at a given concentration from a calibration model are matched to those in the above table.",
									content = "Only applied if a calibration model exists. Rows in the table for which no matches exist are automatically deselected.", 
									placement = "right", trigger = "hover")
							)
						)
					)
				)#,	
				#tabPanel("Import",	
				#	HTML('<hr noshade="noshade" />'),
				#	helpText("Import calibration files from another project (WARNING - replaces calibration files in current project):"),
				#	textInput("import_pro_dir_Cal", "", value = "C:\\...\\other_project_name"),
				#	bsPopover("import_pro_dir_Cal", 
				#		title = "Insert full path, including the project folder, but excluding the logfile.emp.",
				#		content = "Using your OS explorer, you may navigate into your project folder and copy/paste the full path.", 
				#		placement = "right", trigger = "hover"),
				#	bsButton("Import_project_Cal","Import",style="info")	
				#)
			)
        ),
        ########################################################################
        # RESULTS ##############################################################
        ########################################################################
        tabPanel("Results", 	
			tabsetPanel( 
				######################################################################################################################
                tabPanel("Data viewer",
					div(style = widget_style3,numericInput("sel_meas_ID", "Type in file ID:", 0)),
					HTML('<hr noshade="noshade" />'),
					conditionalPanel(			
						condition = "input.sel_meas_ID != 0",					
						bsCollapse(multiple = FALSE, open = NULL, id = "collapse_screen_pos_one",
							bsCollapsePanel(title="All picked peaks & raw data", 
								navbarPage("Settings:",
									tabPanel("Data",
										fluidRow(
											column(3,checkboxInput("peaks_mz_RT_use_peaks", "Plot peaks?", TRUE)),
											column(4,checkboxInput("peaks_mz_RT_use_raw", "Show raw data (as density grid >1E5 data points)?", FALSE)),
											column(3,checkboxInput("peaks_mz_RT_use_IDs", "Display peak IDs?", FALSE))								
										)
									),
									tabPanel("View",
										fluidRow(
											column(5,checkboxInput('showPanel1', 'Show marginal intensity distributions (for <1E5 raw data points)', FALSE)),
											column(5,checkboxInput('showPanel2', 'Show interactive 3D plot (for <1E5 raw data points)', FALSE))
										)
									),
									tabPanel("Search & Filter",
										fluidRow(
											column(4,checkboxInput("peaks_mz_RT_use_window", "Add a blue search window centered at coordinates:", FALSE)),									
											column(2,	
												textInput("peaks_mz_RT_use_window_mass", "m/z",  "216.101"),
												helpText( a("Calculate a mass?", href="http://www.envipat.eawag.ch/index.php",target="_blank"))
											),
											column(2,textInput("peaks_mz_RT_use_window_RT", "RT", "500")),
											column(3,textInput("peaks_mz_RT_use_window_RT_tol", "RT tolerance", "60"))
										),
										fluidRow(
											column(4,checkboxInput("peaks_mz_RT_use_bar", "Add a blue ppm bar to search window?", FALSE)),									
											column(3,textInput("peaks_mz_RT_use_bar_value", "ppm", "10"))
										),
										HTML('<hr noshade="noshade" />'),
										fluidRow(
											column(6,sliderInput("plot_filter_intensity", "Filter intensity range (log10):", min=0, max=8, value=c(0,8), width = '100%', dragRange=TRUE))
										),
										HTML('<hr noshade="noshade" />'),
										helpText("If included in the workflow:"),							
										fluidRow(
											column(4,checkboxInput("plot_filter_blind", "Include blind subtraction?", FALSE)),									
											column(4,checkboxInput("plot_filter_replicates", "Include replicate filter?", FALSE))
										),
										HTML('<hr noshade="noshade" />')
									),
									tabPanel("Hide",HTML('<hr noshade="noshade" />')),
									inverse=FALSE
								),			
								#HTML('<hr noshade="noshade" />'),	
								#fluidRow(
								#	column(2,bsButton("peaks_mz_RT_zoom_out","Zoom out",style="info"))			
								
								#),									
								plotOutput("plot_peaks_mz_RT", 
									width = "100%", height = "650px",
									dblclick = "plot_peaks_mz_RT_dblclick",
									#click = "plot_peaks_mz_RT_click",
									hover = "plot_peaks_mz_RT_hover",
									brush = brushOpts(
										id = "plot_peaks_mz_RT_brush",
										resetOnNew = TRUE
									)
								),
								HTML('<p><font>
									<span style="color:black"><b>&bull;</b></span> Picked peaks <span style="color:gray"><b>&bull;</b></span> Raw data centroids <span style="color:red"><b>&bull;</b></span> Peak centroids 							
								</p></font>'),
								conditionalPanel(condition = 'input.showPanel1',
									plotOutput("plot_peaks_mz_int", width = "100%", height = "320px"),
									plotOutput("plot_peaks_RT_int", width = "100%", height = "320px")
								),
								conditionalPanel(condition = 'input.showPanel2',
									plotly:::plotlyOutput("plot_peaks_3D", width = "100%", height = "800px")
								)
							),
							bsCollapsePanel(title="Individual EICs and peaks", 
								div(style = widget_style3,numericInput("sel_peak_ID", "Specify peak ID:", 0)),
								HTML('<hr noshade="noshade" />'),								
								plotOutput("EIC1", 
									dblclick = "plot_EIC1_dblclick",
									brush = brushOpts(
										id = "plot_EIC1_brush",
										resetOnNew = TRUE
										)
								),
								plotOutput("EIC2", 
									dblclick = "plot_EIC2_dblclick",
									brush = brushOpts(
										id = "plot_EIC2_brush",
										resetOnNew = TRUE
										)
								),
								plotOutput("EIC3", 
									dblclick = "plot_EIC3_dblclick",
									brush = brushOpts(
										id = "plot_EIC3_brush",
										resetOnNew = TRUE
										)
								)						
							)
						)
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
				tabPanel("Processing",            
					numericInput("sel_meas", "Type in file ID:", 0),
					conditionalPanel(			
						condition = "output.dowhat != 'Invalid ID chosen to view processing results.'",
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
						imageOutput("LOD_pic", height="auto"),
						HTML('<hr noshade="noshade" />'),
						div(style = widget_style3,
							bsButton("expo_peaklist","Export peaklist in .csv format",style="info"),
							textOutput("expo2"),
							bsPopover("expo_peaklist", 
								title = "Export peaklist of above selected file",
								content = "Export as peaklist.csv to the export folder of this project, with three columns of mass, intensity and RT. Peaklists are affected by blind and replicate filters contained in the workflow.", 
								placement = "right", trigger = "hover")),
						HTML('<hr noshade="noshade" />')		
					)
				),
				######################################################################################################################
				tabPanel("Screening & quantification",
					tabsetPanel(
						tabPanel("Positive ionization",
							fluidRow(
								column(3, selectInput("Pos_compound_select",label="",choices=c("Choose","Target compounds","Internal standards","Quantification","Recovery","File-wise counts"), 
									selected = "Choose", multiple = FALSE)),
								conditionalPanel(			
									condition = "input.Pos_compound_select == 'Internal standards' || input.Pos_compound_select == 'Target compounds'",										
										column(3, selectInput("screen_pos_summarize", label="", choices = c("Show all adducts"="yes","Collapse adducts"="no"), "yes")),
										column(4, selectInput("Pos_type_select",label="",choices=c("Sample/blind files","Calibration files"), 
											selected = "Non-calibration files", multiple = FALSE))			
								)							
							),
							conditionalPanel(			
								condition = "input.Pos_compound_select == 'Internal standards' || input.Pos_compound_select == 'Target compounds'",	
									conditionalPanel(			
									condition = "input.screen_pos_summarize == 'yes'",	
									bsCollapse(multiple = FALSE, open = NULL, id = "collapse_screen_pos_one",
										bsCollapsePanel(title="Pattern match for selected compound", #style="info",
											textOutput('screening_details_comp_pos'),
											plotOutput("plot_pattern_pos", 
												dblclick = "plot_pattern_pos_dblclick",
												brush = brushOpts(
													id = "plot_pattern_pos_brush",
													resetOnNew = TRUE
												)
											)
										),
										bsCollapsePanel(title="Characteristics for selected compound",
											textOutput('screening_details_comp_pos2'),
											HTML('<hr noshade="noshade" />'),
											fluidRow(										
												column(4, selectInput("selec_pos_x",label="x axis",
													choices=c("m/z","RT","Intensity","Date&time","Type","Place","Conz."),selected = "m/z", multiple = FALSE)),
												column(4, selectInput("selec_pos_y",label="y axis",
													choices=c("m/z","RT","Intensity","Date&time","Type","Place","Conz."),selected = "Intensity", multiple = FALSE)),									
												column(4, radioButtons("selec_pos_log_rat", "Log intensity?", c("yes"="yes","no"="no"),inline=TRUE))
											),
											HTML('<hr noshade="noshade" />'),						
											plotOutput("plot_selec_dist_pos"),
											HTML('<hr noshade="noshade" />'),
											textOutput('info_IS_bounds_pos'),
											conditionalPanel(				
												condition = "(input.Pos_compound_select == 'Internal standards') & (output.info_IS_bounds_pos != 'Compound/adduct not used for quantification')",					
													fluidRow(										
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
											column(4, selectInput("Summ_pos_x",label="x axis",choices=c("m/z","Measured RT","log Intensity","m/z deviation [ppm]","RT deviation within","Time sequence","Expected RT"),selected = "m/z", multiple = FALSE)),
											column(4, selectInput("Summ_pos_y",label="y axis",choices=c("m/z","Measured RT","log Intensity","m/z deviation [ppm]","RT deviation within","Time sequence","Expected RT"),selected = "RT", multiple = FALSE))									
										),
										plotOutput("plot_pattern_distrib_pos"),
										HTML('<hr noshade="noshade" />'),
										fluidRow(										
											column(4,textOutput('count_aboveBlank_pos')),
											column(4,radioButtons("screen_pos_log_rat", "Log scale?", c("yes"="yes","no"="no"),inline=TRUE))
										),
										plotOutput("plot_aboveBlank_pos",height = 250)
									)
								)	
							),		
							conditionalPanel(			
								condition = "input.Pos_compound_select == 'Quantification'",							
								HTML('<hr noshade="noshade" />'),
								DT::dataTableOutput('target_quant_table_pos')
							),		
							conditionalPanel(			
								condition = "input.Pos_compound_select == 'Recovery'",							
								HTML('<hr noshade="noshade" />'),
								DT::dataTableOutput('target_recov_table_pos')
							),												
							conditionalPanel(			
								condition = "input.Pos_compound_select == 'File-wise counts'",	
								tags$p(align="justify","The below table lists the number of compounds which have been positively screened above the cutoff score per file. Matches for different adducts of the same compound are counted separately."),
								HTML('<hr noshade="noshade" />'),
								DT::dataTableOutput('count_file_compound_pos'),
								HTML('<hr noshade="noshade" />')
							)						
						),
						tabPanel("Negative ionization",
							fluidRow(
								column(3, selectInput("Neg_compound_select",label="",choices=c("Choose","Target compounds","Internal standards","Quantification","Recovery","File-wise counts"), 
									selected = "Choose", multiple = FALSE)),
								conditionalPanel(			
									condition = "input.Neg_compound_select == 'Internal standards' || input.Neg_compound_select == 'Target compounds'",										
										column(3, selectInput("screen_neg_summarize", label="", choices = c("Show all adducts"="yes","Collapse adducts"="no"), "yes")),
										column(4, selectInput("Neg_type_select",label="",choices=c("Sample/blind files","Calibration files"), 
											selected = "Non-calibration files", multiple = FALSE))			
								)							
							),
							conditionalPanel(			
								condition = "input.Neg_compound_select == 'Internal standards' || input.Neg_compound_select == 'Target compounds'",	
								conditionalPanel(			
								condition = "input.screen_neg_summarize == 'yes'",	
									bsCollapse(multiple = FALSE, open = NULL, id = "collapse_screen_neg_one",
										bsCollapsePanel(title="Pattern match for selected compound", #style="info",
											textOutput('screening_details_comp_neg'),
											plotOutput("plot_pattern_neg", 
												dblclick = "plot_pattern_neg_dblclick",
												brush = brushOpts(
													id = "plot_pattern_neg_brush",
													resetOnNew = TRUE
												)
											)
										),
										bsCollapsePanel(title="Characteristics for selected compound",
											textOutput('screening_details_comp_neg2'),
											HTML('<hr noshade="noshade" />'),
											fluidRow(										
												column(4, selectInput("selec_neg_x",label="x axis",
													choices=c("m/z","RT","Intensity","Date&time","Type","Place","Conz."),selected = "m/z", multiple = FALSE)),
												column(4, selectInput("selec_neg_y",label="y axis",
													choices=c("m/z","RT","Intensity","Date&time","Type","Place","Conz."),selected = "Intensity", multiple = FALSE)),									
												column(4, radioButtons("selec_neg_log_rat", "Log intensity?", c("yes"="yes","no"="no"),inline=TRUE))
											),
											HTML('<hr noshade="noshade" />'),						
											plotOutput("plot_selec_dist_neg"),
											HTML('<hr noshade="noshade" />'),
											textOutput('info_IS_bounds_neg'),
											conditionalPanel(				
												condition = "(input.Neg_compound_select == 'Internal standards') & (output.info_IS_bounds_neg != 'Compound/adduct not used for quantification')",					
													fluidRow(										
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
											column(4, selectInput("Summ_neg_x",label="x axis",choices=c("m/z","Measured RT","log Intensity","m/z deviation [ppm]","RT deviation within","Time sequence","Expected RT"),selected = "m/z", multiple = FALSE)),
											column(4, selectInput("Summ_neg_y",label="y axis",choices=c("m/z","Measured RT","log Intensity","m/z deviation [ppm]","RT deviation within","Time sequence","Expected RT"),selected = "RT", multiple = FALSE))									
										),
										plotOutput("plot_pattern_distrib_neg"),
										HTML('<hr noshade="noshade" />'),
										fluidRow(										
											column(4,textOutput('count_aboveBlank_neg')),
											column(4,radioButtons("screen_neg_log_rat", "Log scale?", c("yes"="yes","no"="no"),inline=TRUE))
										),
										plotOutput("plot_aboveBlank_neg",height = 250)
									)
								)	
							),		
							conditionalPanel(			
								condition = "input.Neg_compound_select == 'Quantification'",							
								HTML('<hr noshade="noshade" />'),
								DT::dataTableOutput('target_quant_table_neg')
							),								
							conditionalPanel(			
								condition = "input.Neg_compound_select == 'Recovery'",							
								HTML('<hr noshade="noshade" />'),
								DT::dataTableOutput('target_recov_table_neg')
							),	
							conditionalPanel(			
								condition = "input.Neg_compound_select == 'File-wise counts'",	
								tags$p(align="justify","The below table lists the number of compounds which have been positively screened above the cutoff score per file. Matches for different adducts of the same compound are counted separately."),
								HTML('<hr noshade="noshade" />'),
								DT::dataTableOutput('count_file_compound_neg'),
								HTML('<hr noshade="noshade" />')
							)						
						)

					)	
				),
				######################################################################################################################
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
									content = "Get profiles with mean sample intensity x times above mean blank intensities (set Sort profile list by: maximum or mean intensity; x = to be set in blind settings panel). Useful if all your sample files are replicates and not a time sequences.", 
									placement = "top", trigger = "hover"),
								div(style = widget_style3,radioButtons("filterProf_notblind", "Not in blind?", c("no"="no","yes"="yes"))),
								div(style = widget_style3,selectInput("filterProf_sort", "Sort profile list by:", 
									choices = c("ID","mean m/z","mean RT","maximum intensity","mean intensity","global trend intensity","current trend intensity","total peak number"), selected="current trend intensity")),
								div(style = widget_style3,numericInput("filterProf_count", "Restrict list size:", 500)),
								#conditionalPanel( # IS filter				
								#		condition = "input.screen_IS == 'yes'",
								#		tags$h5("IS compounds filter:"),										
								#		HTML('<hr noshade="noshade" />')
								#),
								#conditionalPanel( # target filter				
								#		condition = "input.screen_target == 'yes'",	
								#		tags$h5("Target compounds filter:"),	
								#		HTML('<hr noshade="noshade" />')
								#),															
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
								plotOutput("timeprofile", 
									dblclick = "timeprofile_dblclick",
									brush = brushOpts(
									  id = "timeprofile_brush",
									  resetOnNew = TRUE
									)
								),
								bsCollapse(multiple = FALSE, open = "col1", id = "collapse1",
									bsCollapsePanel("Profile EICs & Peak viewer", 
										div(style = widget_style3,numericInput("profpeakID", "Peak entry #:", min=0, 0)),
										bsPopover("profpeakID", title = "View extracted chromatograms & peaks of the selected profile (sample files only).",
											content = "Select peaks in the order listed in the profile peak table, i.e. from latest to oldest file.", placement = "right", trigger = "hover"),
										div(style = widget_style3,textOutput("prof_peak_text")),
										HTML('<hr noshade="noshade" />'),				
										plotOutput("profile_position"),										
										plotOutput("profile_EIC"),											
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
				)		
            )				
        ),
        ########################################################################
        # HELP & ABOUT #########################################################
        ########################################################################
        tabPanel("Help & About",
			tags$h4("Help"),
			helpText( a("For further help, instructions, topics, examples and requests please visit the enviMass website", href="http://www.envimass.ch/",target="_blank")),
			HTML('<hr noshade="noshade" />'),							
			tags$h4("Citing enviMass"),
			HTML('<p> Martin Loos, 2016, enviMass version 3.1. Zenodo. http://doi.org/10.5281/zenodo.48164 </p> '),
			HTML('<hr noshade="noshade" />'),	
			tags$h4("Contact, author, maintainer:"),
			helpText( a("Martin Loos, mloos@looscomputing.ch", href="http://looscomputing.ch/eng/contact.htm",target="_blank") ),
			tags$h4("Contributors:"),
			helpText("Steffen Ruppe, Matthias Ruff, Jan Mazacek, Heinz Singer"),
			HTML('<hr noshade="noshade" />'),	
			tags$h4("License enviMass version 3.117 :"),
			helpText( a("creative commons Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0)", href="https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode", target="_blank") ),
			helpText("enviMass is distributed in the hope that it will be useful and was tested in different settings, but WITHOUT ANY WARRANTY (see license). 
				TERMS: You must give appropriate credit, indicate changes, not redistribute modifications and not use material for commercial purposes. ")	
		)
        ########################################################################
      ),
	  HTML('<font color="white">') # hide TRUE from sourcing
    )
