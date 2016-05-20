observe({
    input$newit
    if(isolate(input$newit)){
		if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_startup.r!")}
		# clean workspace from any previous workflow results #####################
		source("server_cleaner.R", local=TRUE);		  
		# create alert ###########################################################
		createAlert(session, anchorId="alert_3", alertId="a3", 
				title = "Getting started?", 
				content="<p> Using the below panels, (1) add <span style='color: black'> Files</span>, (2) optionally add <span style='color: black'> Compounds</span>, (3) select the <span style='color: black'> 
					Workflow options </span> and (4) adjust the parameter <span style='color: black'> Settings </span> (make sure your data are suitable for enviPick peak picking, best by trying the seperate 
					<a href='http://cran.r-project.org/web/packages/enviPick/index.html' target='_blank'> package</a> beforehand). Then, press <span style='color: red'> Calculate </span> (a sidebar message will 
					tell you if your project is consistent) and wait for the <span style='color: black'> Results </span> to appear (depending on the number of files, this may take a while; especially if peak 
					picking has not been run yet).
					<p>You can exit your project and reopen it later to add new files and compounds or to change settings. enviMass will then adjust recalculations to avoid unnecessary workflow 
					steps (e.g. to not redo peak picking for all files if only a few new have been added).</p>
					",
				style="info",append=FALSE,dismiss=TRUE)
		##########################################################################
		if(!exists("IS",envir=as.environment(".GlobalEnv"))){data(IS,package="enviMass")}
		if(!exists("targets",envir=as.environment(".GlobalEnv"))){data(targets,package="enviMass")}
		logfile_path<-enviMass:::newproject(isolate(input$pro_name),isolate(input$pro_dir),IS,targets);
		if(logfile_path!="FALSE"){
			output$textit<-renderText(as.character(logfile_path));
			load(logfile_path,envir=as.environment(".GlobalEnv"));
			output$summa_html<-renderText(enviMass:::summary_html(logfile$summary));
			output$dowhat<-renderText("Started new project");
			output$IS<-DT:::renderDataTable(read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character"));
			output$targets<-DT:::renderDataTable(read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character"));      
			measurements<-read.csv(file=file.path(logfile$project_folder,"dataframes","measurements"),colClasses = "character")
			output$measurements<-DT:::renderDataTable(read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character"))      		
			# SET DUMMY RESULTS ####################################################
			# (1) Peak picking #####################################################
			path=file.path(logfile$project_folder,"pics","EIC1");
				png(filename = path, bg = "white", width = 1100, height= 300)
				plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=2)
				dev.off();
				expr_peak<-list(src=file.path(logfile[[1]],"pics","EIC1"));
				output$EIC1<-renderImage(expr_peak, deleteFile = FALSE);
				output$EIC2<-renderImage(expr_peak, deleteFile = FALSE);
			# (1) QC ###############################################################
			path=file.path(logfile$project_folder,"pics","plotQCa_pos")
				png(filename = path, bg = "white")
				plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
				dev.off()
				expr1p<-list(src=path)
				output$plotQCa_pos<-renderImage(expr1p, deleteFile = FALSE)	
			path=file.path(logfile$project_folder,"pics","plotQCb_pos")
				png(filename = path, bg = "white")
				plot.new();plot.window(xlim=c(1,1),ylim=c(1,1));#box();text(1,1,label="not available",cex=1.5,col="darkred")
				dev.off()
				expr2p<-list(src=path)
				output$plotQCb_pos<-renderImage(expr2p, deleteFile = FALSE)
			path=file.path(logfile$project_folder,"pics","plotQCa_neg")
				png(filename = path, bg = "white")
				plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
				dev.off()
				expr1n<-list(src=path)
				output$plotQCa_neg<-renderImage(expr1n, deleteFile = FALSE)	
			path=file.path(logfile$project_folder,"pics","plotQCb_neg")
				png(filename = path, bg = "white")
				plot.new();plot.window(xlim=c(1,1),ylim=c(1,1));#box();text(1,1,label="not available",cex=1.5,col="darkred")
				dev.off()
				expr2n<-list(src=path)
				output$plotQCb_neg<-renderImage(expr2n, deleteFile = FALSE)
			# (3) Normalization ####################################################
			path=file.path(logfile$project_folder,"pics","int_distr_pos")
				png(filename = path, bg = "white")
				plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
				dev.off()
				expr3p<-list(src=path)
				output$pic_int_distr_pos<-renderImage(expr3p, deleteFile = FALSE)
			path=file.path(logfile$project_folder,"pics","int_distr_neg")
				png(filename = path, bg = "white")
				plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
				dev.off()
				expr3n<-list(src=path)
				output$pic_int_distr_neg<-renderImage(expr3n, deleteFile = FALSE)
			# (4) Recalibration ####################################################
			path=file.path(logfile[[1]],"pics","recal_none")
				png(filename = path, bg = "white")
				plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
				dev.off()
				exprrec<-list(src=path)
				output$recal_pic<-renderImage(exprrec, deleteFile = FALSE);		
				output$peakhist_pic<-renderImage(exprrec, deleteFile = FALSE);
				output$peakmzRT_pic<-renderImage(exprrec, deleteFile = FALSE);	
			# (Y) IS normalization #################################################
			path=file.path(logfile$project_folder,"pics","profnorm_pos")
				png(filename = path, bg = "white")
				plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
				dev.off()
				expr30p<-list(src=path)
				output$profnorm_pos<-renderImage(expr30p, deleteFile = FALSE)		
			path=file.path(logfile$project_folder,"pics","profcount_pos")
				png(filename = path, bg = "white")
				plot.new();plot.window(xlim=c(1,1),ylim=c(1,1));#box();text(1,1,label="not available",cex=1.5,col="darkred")
				dev.off()
				expr31p<-list(src=path)
				output$profcount_pos<-renderImage(expr31p, deleteFile = FALSE)
			path=file.path(logfile$project_folder,"pics","profnorm_neg")
				png(filename = path, bg = "white")
				plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
				dev.off()
				expr30n<-list(src=path)
				output$profnorm_neg<-renderImage(expr30n, deleteFile = FALSE)		
			path=file.path(logfile$project_folder,"pics","profcount_neg")
				png(filename = path, bg = "white")
				plot.new();plot.window(xlim=c(1,1),ylim=c(1,1));#box();text(1,1,label="not available",cex=1.5,col="darkred")
				dev.off()
				expr31n<-list(src=path)
				output$profcount_neg<-renderImage(expr31n, deleteFile = FALSE)		
			 

			# (X) Profiling, trends, blind ########################################
			path=file.path(logfile$project_folder,"pics","boxprofile_pos")
				png(filename = path, bg = "white")
				plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
				dev.off()
				expr4p<-list(src=path)
				output$boxprofile<-renderImage(expr4p, deleteFile = FALSE)
			path=file.path(logfile$project_folder,"pics","boxprofile_neg")
				png(filename = path, bg = "white")
				plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
				dev.off()
				expr4n<-list(src=path)
				output$boxprofile<-renderImage(expr4n, deleteFile = FALSE)

			# Parse selectable isotopes ###########################################
			elements<-unique(as.character(isotopes[1:295,1]))
			elements<-elements[order(elements)]
			isotopos<-c()
			for(k in 1:length(elements)){
				if(length(isotopes[isotopes[,1]==elements[k],1])>1){
					elem_isos<-isotopes[isotopes[,1]==elements[k],]
					elem_isos<-as.character(elem_isos[-1,2])
					isotopos<-c(isotopos,elem_isos)
				}
			}
			#updateCheckboxGroupInput(session, "isos", "Select relevant isotopes:", choices = as.character(isotopos),selected=c("13C","34S","81Br","37Cl"))               
			########################################################################
			if(any(ls()=="logfile")){stop("\n illegal logfile detected #1b in server_startup.r!")}
			source("server_variables_in.R", local=TRUE)
			########################################################################
			isolate(init$a<-"TRUE")
			cat("Started a new project\n");
		}else{
			createAlert(session, anchorId="failed_new", alertId = "failed_new_id", title = "Invalid project path", 
				content = "Project already exists or the specified path is invalid.", 
				style = "danger", dismiss = TRUE, append = FALSE)
			cat("Invalid - project already exists or path invalid \n")
		}
    }
	if(any(ls()=="logfile")){stop("\n illegal logfile detected #2 in server_startup.r!")}
})
  
#  observe({ # didnt work with the busy message, remember?
maincalc2<-reactive({
    input$openit
    if(isolate(input$openit)){
		#closeAlert(session, alertId="a3")
		if(any(ls()=="logfile")){stop("\n illegal logfile detected #3 in server_startup.r!")}
		# clean workspace from any previous workflow results #####################
		source("server_cleaner.R", local=TRUE);	  
		##########################################################################
		cat("observed openit")
		file_in<-NA
		try({
			file_in<-as.character(parseFilePaths(getVolumes()(),isolate(input$pro_dir3))[1,4])
		})
		if(is.na(file_in)){ # take string input, format properly
			#cat("\n is NA! \n");cat(file_in);
			file_in<-as.character(isolate(input$pro_dir2))
			if(grepl("\\",file_in,fixed=TRUE)){
				file_in<-gsub("\\",.Platform$file.sep,file_in,fixed=TRUE)
			}
		}else{ # take file selection input - should be properly formatted already
			#cat("\n is not NA! \n");cat(file_in);
			file_in<-strsplit(file_in,.Platform$file.sep)[[1]]			
			file_in<-file_in[-length(file_in)]
			file_in<-paste0(file_in,collapse=.Platform$file.sep)
		}
		
		
		if(file.exists(file.path(file_in,"logfile.emp"))){
			load(file.path(file_in,"logfile.emp"),envir=as.environment(".GlobalEnv"))
			#######################################################################			
			logfile$project_folder<<-as.character(file_in);
			# include version updates #############################################
			source("server_updates.R", local=TRUE);	 
			save(logfile,file=file.path(file_in,"logfile.emp"));
			output$textit<-renderText(logfile$project_folder);
			output$summa_html<-renderText(enviMass:::summary_html(logfile$summary));
			output$dowhat<-renderText("Opened existing project");
			output$IS<-DT::renderDataTable(read.table(file=file.path(logfile$project_folder,"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character"));
			output$targets<-DT::renderDataTable(read.table(file=file.path(logfile$project_folder,"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character"));              
			measurements<-read.csv(file=file.path(logfile$project_folder,"dataframes","measurements"),colClasses = "character")
			output$measurements<-DT::renderDataTable(measurements)        
			# RETRIEVE RESULTS #####################################################
			# (1) Peak picking & preprocessing #####################################
			path=file.path(logfile$project_folder,"pics","EIC1");
				png(filename = path, bg = "white", width = 1100, height= 300)
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,"nothing selected \n or not available",cex=2)
				dev.off();
				expr_peak<-list(src=file.path(logfile$project_folder,"pics","EIC1",sep=""));
				output$EIC1<-renderImage(expr_peak, deleteFile = FALSE);
				output$EIC2<-renderImage(expr_peak, deleteFile = FALSE);
			# (2) QC ###############################################################		
			if(file.exists(file.path(logfile$project_folder,"pics","plotQCa_pos"))){
			  expr1p<-list(src=file.path(logfile$project_folder,"pics","plotQCa_pos"))
			  output$plotQCa_pos<-renderImage(expr1p, deleteFile = FALSE)
			}
			if(file.exists(file.path(logfile$project_folder,"pics","plotQCb_pos"))){
			  expr2p<-list(src=file.path(logfile$project_folder,"pics","plotQCb_pos"))
			  output$plotQCb_pos<-renderImage(expr2p, deleteFile = FALSE)
			}
		   if(file.exists(file.path(logfile$project_folder,"pics","plotQCa_neg"))){
			  expr1n<-list(src=file.path(logfile$project_folder,"pics","plotQCa_neg"))
			  output$plotQCa_neg<-renderImage(expr1n, deleteFile = FALSE)
			}
			if(file.exists(file.path(logfile$project_folder,"pics","plotQCb_neg"))){
			  expr2n<-list(src=file.path(logfile$project_folder,"pics","plotQCb_neg"))
			  output$plotQCb_neg<-renderImage(expr2n, deleteFile = FALSE)
			}
			# (3) Normalization ####################################################
			if(file.exists(file.path(logfile$project_folder,"pics","int_distr_pos"))){
			  expr3p<-list(src=file.path(logfile$project_folder,"pics","int_distr_pos"))
			  output$pic_int_distr_pos<-renderImage(expr3p, deleteFile = FALSE)
			}
			if(file.exists(file.path(logfile$project_folder,"pics","int_distr_neg"))){
			  expr3n<-list(src=file.path(logfile$project_folder,"pics","int_distr_neg"))
			  output$pic_int_distr_neg<-renderImage(expr3n, deleteFile = FALSE)
			}
			# Recalibration, sample peaks ########################################## 
			path=file.path(logfile$project_folder,"pics","recal_none")
				png(filename = path, bg = "white")
				plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
				dev.off()
				exprrec<-list(src=path)
				output$recal_pic<-renderImage(exprrec, deleteFile = FALSE);		
				output$peakhist_pic<-renderImage(exprrec, deleteFile = FALSE);
				output$peakmzRT_pic<-renderImage(exprrec, deleteFile = FALSE);
			
			
			# (4) Available measurements ###########################################
			# if(logfile$workflow[2]=="yes"){	
				updateSelectInput(session, "sel_meas", label = "Select file by ID:", choices =  c("none",as.character(measurements[,1])), selected = "none")
			# }
			updateSelectInput(session, "sel_meas_ID", label = "Select file by ID:", choices =  c("none",as.character(measurements[,1])), selected = "none")
			# (5) RT Alignment #####################################################
						
			# (6) IS-Normalization #################################################
			if(file.exists(file.path(logfile$project_folder,"pics","profnorm_pos"))){
				if(isolate(input$Ion_mode)=="positive"){
					exprprofnorm_pos<-list(src=file.path(logfile$project_folder,"pics","profnorm_pos"))
					output$profnorm<-renderImage(exprprofnorm_pos, deleteFile = FALSE)
				}
			}
			if(file.exists(file.path(logfile$project_folder,"pics","profcount_pos"))){
				if(isolate(input$Ion_mode)=="positive"){
					exprprofcount_pos<-list(src=file.path(logfile$project_folder,"pics","profcount_pos"))
					output$profcount<-renderImage(exprprofcount_pos, deleteFile = FALSE)
				}
			}
			if(file.exists(file.path(logfile$project_folder,"pics","profnorm_neg"))){
				if(isolate(input$Ion_mode)=="negative"){
					exprprofnorm_neg<-list(src=file.path(logfile$project_folder,"pics","profnorm_neg"))
					output$profnorm<-renderImage(exprprofnorm_neg, deleteFile = FALSE)
				}
			}
			if(file.exists(file.path(logfile$project_folder,"pics","profcount_neg"))){
				if(isolate(input$Ion_mode)=="negative"){
					exprprofcount_neg<-list(src=file.path(logfile$project_folder,"pics","profcount_neg"))
					output$profcount<-renderImage(exprprofcount_neg, deleteFile = FALSE)
				}
			}			
			# (X) Profiling, trends, blind #########################################		
			if(file.exists(file.path(logfile$project_folder,"pics","boxprofile_pos"))){
				if(isolate(input$Ion_mode)=="positive"){
					expr4p<-list(src=file.path(logfile$project_folder,"pics","boxprofile_pos"))
					output$boxprofile<-renderImage(expr4p, deleteFile = FALSE)
				}
			}
			if(file.exists(file.path(logfile$project_folder,"pics","boxprofile_neg"))){
				if(isolate(input$Ion_mode)=="negative"){
					expr4n<-list(src=file.path(logfile$project_folder,"pics","boxprofile_neg"))
					output$boxprofile<-renderImage(expr4n, deleteFile = FALSE)
				}
			}
			if(file.exists(file.path(logfile$project_folder,"results","profileList_pos"))){
				if(isolate(input$Ion_mode)=="positive"){
					load(file=file.path(as.character(logfile$project_folder),"results","profileList_pos"),envir=as.environment(".GlobalEnv"), verbose=TRUE);
					assign("profileList",profileList_pos,envir=as.environment(".GlobalEnv"));
				}
			}	
			if(file.exists(file.path(logfile$project_folder,"results","profileList_neg"))){
				if(isolate(input$Ion_mode)=="negative"){
					load(file=file.path(as.character(logfile$project_folder),"results","profileList_neg"),envir=as.environment(".GlobalEnv"), verbose=TRUE);
					assign("profileList",profileList_neg,envir=as.environment(".GlobalEnv"));
				}
			}	
			if(file.exists(file.path(logfile$project_folder,"results","profpeaks_pos"))){
				if(isolate(input$Ion_mode)=="positive"){
					load(file=file.path(as.character(logfile$project_folder),"results","profpeaks_pos"),envir=as.environment(".GlobalEnv"), verbose=TRUE);
					assign("profpeaks",profpeaks_pos,envir=as.environment(".GlobalEnv"));
				}
			}
			if(file.exists(file.path(logfile$project_folder,"results","profpeaks_neg"))){
				if(isolate(input$Ion_mode)=="negative"){
					load(file=file.path(as.character(logfile$project_folder),"results","profpeaks_neg"),envir=as.environment(".GlobalEnv"), verbose=TRUE);
					assign("profpeaks",profpeaks_neg,envir=as.environment(".GlobalEnv"));
				}
			}
			

			if(file.exists(file.path(logfile$project_folder,"pics","profilehisto.png"))){ 
				expr6<-list(src=file.path(logfile$project_folder,"pics","profilehisto.png"))
				output$profilehisto<-renderImage(expr6, deleteFile = FALSE)
			}		
			# Parse selectable isotopes ###########################################
			elements<-unique(as.character(isotopes[1:295,1]))
			elements<-elements[order(elements)]
			isotopos<-c()
			for(k in 1:length(elements)){
				if(length(isotopes[isotopes[,1]==elements[k],1])>1){
					elem_isos<-isotopes[isotopes[,1]==elements[k],]
					elem_isos<-as.character(elem_isos[-1,2])
					isotopos<-c(isotopos,elem_isos)
				}
			}
			#updateCheckboxGroupInput(session, "isos", "Select relevant isotopes:", choices = as.character(isotopos),selected=c("13C","34S","81Br","37Cl"))               		
			########################################################################  
			# screening results - show only target results by default ##############
			if( isolate(input$Pos_compound_select=="Target compounds") ){
				cat("Looking at positive targets_startup \n")
				if( file.exists(file=file.path(logfile$project_folder,"results","screening","results_screen_target_pos")) ){
					load(file=file.path(logfile$project_folder,"results","screening","results_screen_target_pos"))
					screen_dev_pos<-results_screen_target_pos[[3]]
					rat_sam_blank_pos<-results_screen_target_pos[[1]][,10,drop=FALSE]
					if( isolate(input$screen_pos_summarize=="yes") ){
						results_screen<-results_screen_target_pos[[1]]
					}else{
						results_screen<-results_screen_target_pos[[2]]
					}
					output$Table_screening_pos <- DT::renderDataTable({
						DT::datatable(results_screen, escape = FALSE,selection = 'single',rownames=FALSE) %>% 
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
			if( isolate(input$Pos_compound_select=="Target compounds") ){
				cat("Looking at negative targets_startup \n")
				if( file.exists(file=file.path(logfile$project_folder,"results","screening","results_screen_target_neg")) ){
					load(file=file.path(logfile$project_folder,"results","screening","results_screen_target_neg"))
					screen_dev_neg<-results_screen_target_neg[[3]]
					rat_sam_blank_neg<-results_screen_target_neg[[1]][,10,drop=FALSE]
					if( isolate(input$screen_neg_summarize=="yes") ){
						results_screen<-results_screen_target_neg[[1]]
					}else{
						results_screen<-results_screen_target_neg[[2]]
					}
					output$Table_screening_neg <- DT::renderDataTable({
						DT::datatable(results_screen, escape = FALSE,selection = 'single',rownames=FALSE) %>% 
							formatStyle(
								'Max. sample score',
								background = styleColorBar(c(0,1), 'lightgreen'),
								backgroundPosition = 'left'
							)					
					},server = TRUE)
					rm(results_screen_target_neg,table_screening)	
				}else{	
					output$Table_screening_neg <- DT::renderDataTable({
						DT::datatable(as.data.frame(cbind("")),selection = 'single',rownames=FALSE,colnames="No target screening results available")
					},server = TRUE)		
				}
			}
			###########################################################################	
			cat(objects())
			if(any(ls()=="logfile")){stop("\n illegal logfile detected #3b in server_startup.r!")}
			source("server_variables_in.R", local=TRUE)
			if(isolate(init$a=="FALSE")){
				isolate(init$a<-"TRUE")
			}else{
				isolate(init$b<-(init$b+1));cat(" - ")
			}
			cat("\nProject opened\n")
			return("Project available\n")
		}else{
			createAlert(session, anchorId="failed_open", alertId = "failed_open_id", title = "Invalid project path", 
				content = "Project with specified path does not exist!", 
				style = "danger", dismiss = TRUE, append = FALSE)
			cat("Invalid - project already exists or path invalid \n")
			cat("Invalid project!\n")
			return("Project invalid\n")
		}
    }
	if(any(ls()=="logfile")){stop("\n illegal logfile detected #3 in server_startup.r!")}
  })
output$had_opened<-renderText(paste(maincalc2())) 




