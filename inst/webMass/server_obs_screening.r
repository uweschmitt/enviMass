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
				load(file=file.path(as.character(logfile[[1]]),"results","profileList_pos"),envir=as.environment(".GlobalEnv"));	
				load(file=file.path(logfile[[1]],"results","pattern_pos_target"),envir=as.environment(".GlobalEnv"));
				pattern<-pattern_pos_target;rm(pattern_pos_target,envir=as.environment(".GlobalEnv"));
				load(file=file.path(logfile[[1]],"results","patternRT_pos_target"),envir=as.environment(".GlobalEnv"));
				pattern_RT<-patternRT_pos_target;rm(patternRT_pos_target,envir=as.environment(".GlobalEnv"));
				load(file=file.path(logfile[[1]],"results","patternDelRT_pos_target"),envir=as.environment(".GlobalEnv"));
				pattern_delRT<-patternDelRT_pos_target;rm(patternDelRT_pos_target,envir=as.environment(".GlobalEnv"));
				load(file=file.path(logfile$project_folder,"results","screening","res_target_pos_screen"))
				res_pos_screen<-res_target_pos_screen;rm(res_target_pos_screen);
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				compound_table<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");				
			}else{	
				output$Table_screening_pos <- DT::renderDataTable({
					DT::datatable(as.data.frame(cbind("")),selection = 'single',rownames=FALSE,colnames="No target screening results available")
				},server = TRUE)	
				found_table<-FALSE				
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
				load(file=file.path(as.character(logfile[[1]]),"results","profileList_pos"),envir=as.environment(".GlobalEnv"));	
				load(file=file.path(logfile[[1]],"results","pattern_pos_IS"),envir=as.environment(".GlobalEnv"));
				pattern<-pattern_pos_IS;rm(pattern_pos_IS,envir=as.environment(".GlobalEnv"));
				load(file=file.path(logfile[[1]],"results","patternRT_pos_IS"),envir=as.environment(".GlobalEnv"));
				pattern_RT<-patternRT_pos_IS;rm(patternRT_pos_IS,envir=as.environment(".GlobalEnv"));
				load(file=file.path(logfile[[1]],"results","patternDelRT_pos_IS"),envir=as.environment(".GlobalEnv"));
				pattern_delRT<-patternDelRT_pos_IS;rm(patternDelRT_pos_IS,envir=as.environment(".GlobalEnv"));
				load(file=file.path(logfile$project_folder,"results","screening","res_IS_pos_screen"))
				res_pos_screen<-res_IS_pos_screen;rm(res_IS_pos_screen);
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				compound_table<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
			}else{	
				output$Table_screening_pos <- DT::renderDataTable({
					DT::datatable(as.data.frame(cbind("")),selection = 'single',rownames=FALSE,colnames="No internal standard screening results available")
				},server = TRUE)	
				found_table<-FALSE				
			}
		}
		if(found_table){
			# name selected compound & adduct
			named_compound <- renderText({
				s<-input$Table_screening_pos_row_last_clicked
				if (length(s)) {	paste(results_screen_pos[s,2],results_screen_pos[s,3],
				compound_table[compound_table[,1]==results_screen_pos[s,1],3]
				,sep=", ")	}
			})
			output$screening_details_comp_pos <-named_compound
			output$screening_details_comp_pos2<-named_compound
			output$screening_details_comp_pos3<-named_compound
			# plot pattern & matches
			output$plot_pattern <- renderPlot({
				s<-input$Table_screening_pos_row_last_clicked
				if (length(s) & isolate(input$screen_pos_summarize=="yes")) {	
					use_comp<-(	grepl(paste(results_screen_pos[s,1],"_",sep=""),names(pattern),fixed=TRUE)&
								grepl(paste("_",results_screen_pos[s,3],"_",sep=""),names(pattern),fixed=TRUE))
					pattern_sel<-pattern[use_comp][[1]]
					res_pos_screen_sel<-res_pos_screen[use_comp][[1]]
					int_tol<-as.numeric(logfile$parameters$IS_inttol)
					ylim_up<-(int_tol+100+(int_tol*0.1))
					plot(pattern_sel[,1],pattern_sel[,2],type="h",lwd=3,col="red",xlab="m/z",ylab="Rescaled intensity",ylim=c(0,ylim_up))
					if(length(res_pos_screen_sel)>0){
						for(i in 1:length(res_pos_screen_sel)){
							if(length(res_pos_screen_sel[[i]])>0){
								for(j in 1:length(res_pos_screen_sel[[i]])){
									found_matches<-res_pos_screen_sel[[i]][[j]]$Peaks
									x<-profileList_pos[[2]][found_matches[,2],1]
									y<-(profileList_pos[[2]][found_matches[,2],2]*res_pos_screen_sel[[i]][[j]][[6]])
									y<-y[order(x)]
									x<-x[order(x)]
									lines(x=x,y=y,col="grey")										
									points(x,y,col="darkgreen",cex=2)
								}
							}
						}
					}
					points(pattern_sel[,1],pattern_sel[,2],type="h",lwd=3,col="red",xlab="m/z",ylab="Rescaled intensity",ylim=c(0,110))
					plot.window(xlim=c(0,1),ylim=c(0,1))
					legend(0.8,1,legend=c("Theoretical pattern","Matches","Co-occurrences"),fill=c("red","darkgreen","grey"),border=c("red","darkgreen","grey"))
				}
			})
			# make Table over samples
			output$Table_screening_selected_pos<-DT::renderDataTable({
				s<-input$Table_screening_pos_row_last_clicked
				if (length(s) & isolate(input$screen_pos_summarize=="yes")) {
					use_comp<-(	grepl(paste(results_screen_pos[s,1],"_",sep=""),names(pattern),fixed=TRUE)&
								grepl(paste("_",results_screen_pos[s,3],"_",sep=""),names(pattern),fixed=TRUE) )
					res_pos_screen_sel<-res_pos_screen[use_comp][[1]]
					which_where<-c();which_peaks<-c();sample_type<-c();score_1<-c();score_2<-c();delppm<-c();delRT<-c();inte<-c();
					with_peaks<-c();with_file<-c();with_s<-c();
					
					if(length(res_pos_screen_sel)>0){
						for(i in 1:length(res_pos_screen_sel)){
							if(length(res_pos_screen_sel[[i]])>0){
								for(j in 1:length(res_pos_screen_sel[[i]])){
									which_where<-c(which_where,measurements[i,1]);
									sample_type<-c(sample_type,measurements[i,3]);
									which_peaks<-c(which_peaks,paste(res_pos_screen_sel[[i]][[j]]$Peaks[,1],collapse=", "))
									score_1<-c(score_1,round(res_pos_screen_sel[[i]][[j]]$score_1,digits=2));
									score_2<-c(score_2,round(res_pos_screen_sel[[i]][[j]]$score_2,digits=2));
									delppm<-c(delppm,paste(as.character(round(res_pos_screen_sel[[i]][[j]][[4]],digits=2)),collapse=", "));
									#delRT<-c(delRT,paste(as.character(round(res_pos_screen_sel[[i]][[j]][[4]],digits=2)),collapse=", "));
									found_matches<-res_pos_screen_sel[[i]][[j]]$Peaks
									with_peaks<-c(with_peaks,paste(as.character(found_matches[,2]),collapse=", "));
									delRT<-c(delRT,
										paste(as.character(round(profileList_pos[[2]][found_matches[,2],3],digits=2)),collapse=", ")
									)
									inte<-c(inte,
										paste(as.character(round(log10(profileList_pos[[2]][found_matches[,2],2]),digits=2)),collapse=", ")									
									)
									with_file<-c(with_file,i)
									with_s<-c(with_s,s)
								}
							}
						}
					}
					DT::datatable(
						as.data.frame(cbind(which_where,sample_type,which_peaks,score_1,score_2,delppm,delRT,inte,
						with_peaks,with_file,with_s),row.names = NULL,stringsAsFactors=FALSE),
						rownames = FALSE, colnames=c("File ID","File type","Pattern matches","Score > LOD","Score < LOD",
							"m/z deviation (ppm)","RT","log Intensity","Peak IDs","m","i")
					)
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