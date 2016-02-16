if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_obs_screening.r!")}

##############################################################################
# update screening results ###################################################
##############################################################################
# Dont understand the environments in which results_screen_IS_pos ends up! ###
# POSITIVE IONIZATION ########################################################
observe({ 
	input$Pos_compound_select  
	input$screen_pos_summarize
	init$b
	if(isolate(init$a)=="TRUE"){
		if( isolate(input$Pos_compound_select=="Target compounds") ){
			cat("\n Looking at positive targets")
			if( file.exists(file=file.path(logfile$project_folder,"results","screening","results_screen_target_pos")) ){
				load(file=file.path(logfile$project_folder,"results","screening","results_screen_target_pos"))
				if( isolate(input$screen_pos_summarize=="yes") ){
					results_screen_target_pos<<-results_screen_target_pos[[1]]
				}else{
					results_screen_target_pos<<-results_screen_target_pos[[2]]
				}
				output$Table_screening_pos <- DT::renderDataTable({
					DT::datatable(results_screen_target_pos, escape = FALSE,selection = 'single',rownames=FALSE) %>% 
						formatStyle(
							'Max. sample score',
							background = styleColorBar(c(0,1), 'lightgreen'),
							backgroundPosition = 'left'
						)
				},server = TRUE)
				rm(results_screen_target_pos)	
			}else{	
				output$Table_screening_pos <- DT::renderDataTable({
					DT::datatable(as.data.frame(cbind("")),selection = 'single',rownames=FALSE,colnames="No target screening results available")
				},server = TRUE)		
			}
		}
		if( isolate(input$Pos_compound_select=="Internal standards") ){
			cat("\n Looking at positive standards")
			if( file.exists(file=file.path(logfile$project_folder,"results","screening","results_screen_IS_pos")) ){
				load(file=file.path(logfile$project_folder,"results","screening","results_screen_IS_pos"))
				if( isolate(input$screen_pos_summarize=="yes") ){
					results_screen_IS_pos<<-results_screen_IS_pos[[1]]
				}else{
					results_screen_IS_pos<<-results_screen_IS_pos[[2]]
				}
				output$Table_screening_pos <- DT::renderDataTable({
					DT::datatable(results_screen_IS_pos, escape = FALSE,selection = 'single',rownames=FALSE) %>% 
						formatStyle(
							'Max. sample score',
							background = styleColorBar(c(0,1), 'lightgreen'),
							backgroundPosition = 'left'
						)
				},server = TRUE)
				rm(results_screen_IS_pos)	
			}else{	
				output$Table_screening_pos <- DT::renderDataTable({
					DT::datatable(as.data.frame(cbind("")),selection = 'single',rownames=FALSE,colnames="No internal standard screening results available")
				},server = TRUE)		
			}
		}
	}
})  
# NEGATIVE IONIZATION ########################################################
observe({ 
	input$Neg_compound_select  
	input$screen_neg_summarize
	init$b
	if(isolate(init$a)=="TRUE"){
		if( isolate(input$Neg_compound_select=="Target compounds") ){
			cat("\n Looking at negative targets")
			if( file.exists(file=file.path(logfile$project_folder,"results","screening","results_screen_target_neg")) ){
				load(file=file.path(logfile$project_folder,"results","screening","results_screen_target_neg"))
				if( isolate(input$screen_neg_summarize=="yes") ){
					results_screen_target_neg<<-results_screen_target_neg[[1]]
				}else{
					results_screen_target_neg<<-results_screen_target_neg[[2]]
				}
				output$Table_screening_neg <- DT::renderDataTable({
					DT::datatable(results_screen_target_neg, escape = FALSE,selection = 'single',rownames=FALSE) %>% 
						formatStyle(
							'Max. sample score',
							background = styleColorBar(c(0,1), 'lightgreen'),
							backgroundPosition = 'left'
						)
				},server = TRUE)
				rm(results_screen_target_neg)	
			}else{	
				output$Table_screening_neg <- DT::renderDataTable({
					DT::datatable(as.data.frame(cbind("")),selection = 'single',rownames=FALSE,colnames="No target screening results available")
				},server = TRUE)		
			}
		}
		if( isolate(input$Neg_compound_select=="Internal standards") ){
			cat("\n Looking at negative standards")
			if( file.exists(file=file.path(logfile$project_folder,"results","screening","results_screen_IS_neg")) ){
				load(file=file.path(logfile$project_folder,"results","screening","results_screen_IS_neg"))
				if( isolate(input$screen_neg_summarize=="yes") ){
					results_screen_IS_neg<<-results_screen_IS_neg[[1]]
				}else{
					results_screen_IS_neg<<-results_screen_IS_neg[[2]]
				}
				output$Table_screening_neg <- DT::renderDataTable({
					DT::datatable(results_screen_IS_neg, escape = FALSE,selection = 'single',rownames=FALSE) %>% 
						formatStyle(
							'Max. sample score',
							background = styleColorBar(c(0,1), 'lightgreen'),
							backgroundPosition = 'left'
						)
				},server = TRUE)
				rm(results_screen_IS_neg)	
			}else{	
				output$Table_screening_neg <- DT::renderDataTable({
					DT::datatable(as.data.frame(cbind("")),selection = 'single',rownames=FALSE,colnames="No internal standard screening results available")
				},server = TRUE)		
			}
		}
	}
})  
############################################################################## 

##############################################################################
# Plot screening results #####################################################
##############################################################################
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
############################################################################## 




















if(any(ls()=="logfile")){stop("\n illegal logfile detected #2 in server_obs_screening.r!")}