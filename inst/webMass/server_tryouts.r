
if(FALSE){			
if(file.exists(file.path(logfile$project_folder,"results","screening","IS_screening_pos"))){			
	load(file=file.path(as.character(logfile$project_folder),"results","screening","IS_screening_pos"),envir=as.environment(".GlobalEnv"), verbose=TRUE);

	addResourcePath("project22", file.path(as.character(logfile$project_folder)))
	
	dir.create(file.path(path.package("shiny", quiet = FALSE),"www","screening_pics"), showWarnings = FALSE, recursive = TRUE, mode = "0777")
	copy_files<-list.files(file.path(as.character(logfile$project_folder),"results","screening"))
	if(length(copy_files)>0){	for(i in 1:length(copy_files)){
		file.copy(	
				from=file.path(as.character(logfile$project_folder),"results","screening",copy_files[i]),
				to=file.path(path.package("shiny", quiet = FALSE),"www","screening_pics",copy_files[i]),				
				overwrite = TRUE, 
				recursive = FALSE,
				copy.mode = TRUE
		)
	}}
	

	#output$Table_IS_screening_pos <- DT::renderDataTable({
	#	DT::datatable(IS_screening_pos, escape = FALSE,selection = 'single') # HERE
	#})
	output$Table_IS_screening_pos <- DT::renderDataTable({
			DT::datatable(IS_screening_pos, escape = FALSE,selection = 'single') %>% 
				formatStyle(
					'Score',
					background = styleColorBar(c(0,1), 'lightgreen'),
					#backgroundSize = '1 0',
					#backgroundRepeat = 'no-repeat',
					backgroundPosition = 'left'
		)
		
		},
		server = TRUE
	)
	
	used_rows<-c() # ISSUE with selection! s requires a deselection before a row can be loaded anew
	output$Table_IS_screening_pos_row = renderPrint({
		s = input$Table_IS_screening_pos_rows_selected
		#s<-s[length(s)]
		if (length(s)) {
			cat('These rows were selected:\n\n')
			#cat(IS_screening_pos[s,3], sep = ', ')
			cat(s, sep = ', ')
		}
	})
		
}			
}			
			
			
			
			
			
			
			
			
			
			
			
			
			
unload(pkg = "D:\\Users\\uchemadmin\\Desktop\\MS\\R_packages\\enviMass\\enviMass\\enviMass")
clean_dll(pkg = "D:\\Users\\uchemadmin\\Desktop\\MS\\R_packages\\enviMass\\enviMass\\enviMass")
compile_dll(pkg = "D:\\Users\\uchemadmin\\Desktop\\MS\\R_packages\\enviMass\\enviMass\\enviMass", quiet = FALSE)
load_all(pkg = "D:\\Users\\uchemadmin\\Desktop\\MS\\R_packages\\enviMass\\enviMass\\enviMass")
#set_path(c(old_path,get_path()))
