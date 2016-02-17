if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_obs_screening.r!")}

##############################################################################
# update screening results ###################################################
##############################################################################
# Dont understand the environments:: in which table_screening ends up! #######
# POSITIVE IONIZATION ########################################################
observe({ 
	input$Pos_compound_select  
	input$screen_pos_summarize
	init$b
	if(isolate(init$a)=="TRUE"){
		found_table<-FALSE
		if( isolate(input$Pos_compound_select=="Target compounds") ){
			cat("\n Looking at positive targets_selec")
			if( file.exists(file=file.path(logfile$project_folder,"results","screening","results_screen_target_pos")) ){
				load(file=file.path(logfile$project_folder,"results","screening","results_screen_target_pos"))
				if( isolate(input$screen_pos_summarize=="yes") ){
					results_screen_pos<-results_screen_target_pos[[1]]
				}else{
					results_screen_pos<-results_screen_target_pos[[2]]
				}
				table_screening<-DT::datatable(results_screen_pos, escape = FALSE,selection = 'single') %>% 
					formatStyle('Max. sample score',background = styleColorBar(c(0,1), 'lightgreen'),backgroundPosition = 'right')
				output$Table_screening_pos <- DT::renderDataTable({table_screening},server = TRUE)
				rm(results_screen_target_pos)
				found_table<-TRUE
			}else{	
				output$Table_screening_pos <- DT::renderDataTable({
					DT::datatable(as.data.frame(cbind("")),selection = 'single',rownames=FALSE,colnames="No target screening results available")
				},server = TRUE)		
			}
		}
		if( isolate(input$Pos_compound_select=="Internal standards") ){
			cat("\n Looking at positive standards_selec")
			if( file.exists(file=file.path(logfile$project_folder,"results","screening","results_screen_IS_pos")) ){
				load(file=file.path(logfile$project_folder,"results","screening","results_screen_IS_pos"))
				if( isolate(input$screen_pos_summarize=="yes") ){
					results_screen_pos<-results_screen_IS_pos[[1]]
				}else{
					results_screen_pos<-results_screen_IS_pos[[2]]
				}
				table_screening<-DT::datatable(results_screen_pos, escape = FALSE,selection = 'single') %>% 
					formatStyle('Max. sample score',background = styleColorBar(c(0,1), 'lightgreen'),backgroundPosition = 'right')
				output$Table_screening_pos <- DT::renderDataTable({table_screening},server = TRUE)
				rm(results_screen_IS_pos)	
				found_table<-TRUE
			}else{	
				output$Table_screening_pos <- DT::renderDataTable({
					DT::datatable(as.data.frame(cbind("")),selection = 'single',rownames=FALSE,colnames="No internal standard screening results available")
				},server = TRUE)		
			}
		}
		if(found_table){
		
			output$Table_screening_pos_row <- renderPrint({
				
				s<-input$Table_screening_pos_rows_selected
				
				#s<-s[is.na(match(s,s_pos))]
				s_pos<<-c(s)
				
				if (length(s)) {
					cat('These rows were selected:\n\n')
					#cat(results_screen_pos[s,5], sep = ', ')
					cat(s, sep = ', ')
				}
			
			
			})
		}
		
		
	}
})  
# NEGATIVE IONIZATION ########################################################
observe({ 
	input$Neg_compound_select  
	input$screen_neg_summarize
	init$b
	if(isolate(init$a)=="TRUE"){
		if( isolate(input$Neg_compound_select=="Target compounds_selec") ){
			cat("\n Looking at negative targets")
			if( file.exists(file=file.path(logfile$project_folder,"results","screening","results_screen_target_neg")) ){
				load(file=file.path(logfile$project_folder,"results","screening","results_screen_target_neg"))
				if( isolate(input$screen_neg_summarize=="yes") ){
					results_screen<-results_screen_target_neg[[1]]
				}else{
					results_screen<-results_screen_target_neg[[2]]
				}
				table_screening<-DT::datatable(results_screen, escape = FALSE,selection = 'single',rownames=FALSE) %>% 
					formatStyle('Max. sample score',background = styleColorBar(c(0,1), 'lightgreen'),backgroundPosition = 'right')
				output$Table_screening_neg <- DT::renderDataTable({table_screening},server = TRUE)
				rm(results_screen_target_neg)	
			}else{	
				output$Table_screening_neg <- DT::renderDataTable({
					DT::datatable(as.data.frame(cbind("")),selection = 'single',rownames=FALSE,colnames="No target screening results available")
				},server = TRUE)		
			}
		}
		if( isolate(input$Neg_compound_select=="Internal standards_selec") ){
			cat("\n Looking at negative standards")
			if( file.exists(file=file.path(logfile$project_folder,"results","screening","results_screen_IS_neg")) ){
				load(file=file.path(logfile$project_folder,"results","screening","results_screen_IS_neg"))
				if( isolate(input$screen_neg_summarize=="yes") ){
					results_screen<-results_screen_IS_neg[[1]]
				}else{
					results_screen<-results_screen_IS_neg[[2]]
				}
				table_screening<-DT::datatable(results_screen, escape = FALSE,selection = 'single',rownames=FALSE) %>% 
					formatStyle('Max. sample score',background = styleColorBar(c(0,1), 'lightgreen'),backgroundPosition = 'right')
				output$Table_screening_neg <- DT::renderDataTable({table_screening},server = TRUE)
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


if(any(ls()=="logfile")){stop("\n illegal logfile detected #2 in server_obs_screening.r!")}