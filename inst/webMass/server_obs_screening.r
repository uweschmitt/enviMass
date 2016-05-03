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
			if( 
				(file.exists(file=file.path(logfile$project_folder,"results","screening","results_screen_target_pos"))) #&
				#(logfile$workflow[names(logfile$workflow)=="target_subtr"]!="yes")
			){			
				load(file=file.path(logfile$project_folder,"results","screening","results_screen_target_pos"))
				screen_dev_pos<-results_screen_target_pos[[3]] # contains sample vs. blank intensity ratios
				screen_dev_pos[,2]<-log10(screen_dev_pos[,2])
				rat_sam_blank_pos<-results_screen_target_pos[[1]][,10,drop=FALSE]
				if( isolate(input$screen_pos_summarize=="yes") ){
					results_screen_pos<-results_screen_target_pos[[1]]
				}else{
					results_screen_pos<-results_screen_target_pos[[2]]
				}
				table_screening_pos<-DT::datatable(results_screen_pos, escape = FALSE,selection = 'single') %>% 
					formatStyle('Max. sample score',background = styleColorBar(c(0,1), 'lightgreen'),backgroundPosition = 'right')
				output$Table_screening_pos <- DT::renderDataTable({table_screening_pos},server = TRUE)
				rm(results_screen_target_pos)
				found_table<-TRUE
				load(file=file.path(as.character(logfile[[1]]),"results","profileList_pos_copy"),envir=as.environment(".GlobalEnv"));	
				load(file=file.path(logfile[[1]],"results","pattern_pos_target"),envir=as.environment(".GlobalEnv"));
				pattern_pos<-pattern_pos_target;rm(pattern_pos_target,envir=as.environment(".GlobalEnv"));
				patt_pos_ID<-rep("",length(pattern_pos))
				patt_pos_add<-rep("",length(pattern_pos))
				for(i in 1:length(pattern_pos)){
					patt_pos_ID[i]<-strsplit(names(pattern_pos[i]),"_",fixed=TRUE)[[1]][1]
					patt_pos_add[i]<-strsplit(names(pattern_pos[i]),"_",fixed=TRUE)[[1]][2]
				}
				load(file=file.path(logfile[[1]],"results","patternRT_pos_target"),envir=as.environment(".GlobalEnv"));
				pattern_RT_pos<-patternRT_pos_target;rm(patternRT_pos_target,envir=as.environment(".GlobalEnv"));
				load(file=file.path(logfile[[1]],"results","patternDelRT_pos_target"),envir=as.environment(".GlobalEnv"));
				pattern_delRT_pos<-patternDelRT_pos_target;rm(patternDelRT_pos_target,envir=as.environment(".GlobalEnv"));
				load(file=file.path(logfile$project_folder,"results","screening","res_target_pos_screen"))
				res_pos_screen<-res_target_pos_screen;rm(res_target_pos_screen);
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				compound_table<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");	
				# For IS without Conz.:
				updateSelectInput(session,inputId="selec_pos_x",label="x axis",choices=c("m/z","RT","Intensity","Date&time","Type","Place","Conz."),selected = "m/z")
				updateSelectInput(session,inputId="selec_pos_y",label="y axis",choices=c("m/z","RT","Intensity","Date&time","Type","Place","Conz."),selected = "RT")
			}else{	
				output$Table_screening_pos <- DT::renderDataTable({
					DT::datatable(as.data.frame(cbind("")),selection = 'single',rownames=FALSE,colnames="No target screening results available")
				},server = TRUE)	
				found_table<-FALSE				
			}
		}
		if( isolate(input$Pos_compound_select=="Internal standards") ){
			cat("\n Looking at positive standards_selec")
			if( 
				(file.exists(file=file.path(logfile$project_folder,"results","screening","results_screen_IS_pos"))) #&
				#(logfile$workflow[names(logfile$workflow)=="IS_subtr"]!="yes")
			){			
				load(file=file.path(logfile$project_folder,"results","screening","results_screen_IS_pos"))
				screen_dev_pos<-results_screen_IS_pos[[3]] # contains sample vs. blank intensity ratios		
				screen_dev_pos[,2]<-log10(screen_dev_pos[,2])
				rat_sam_blank_pos<-results_screen_IS_pos[[1]][,10,drop=FALSE]				
				if( isolate(input$screen_pos_summarize=="yes") ){
					results_screen_pos<-results_screen_IS_pos[[1]]
				}else{
					results_screen_pos<-results_screen_IS_pos[[2]]
				}
				table_screening_pos<-DT::datatable(results_screen_pos, escape = FALSE,selection = 'single') %>% 
					formatStyle('Max. sample score',background = styleColorBar(c(0,1), 'lightgreen'),backgroundPosition = 'right')
				output$Table_screening_pos <- DT::renderDataTable({table_screening_pos},server = TRUE)
				rm(results_screen_IS_pos)	
				found_table<-TRUE
				load(file=file.path(as.character(logfile[[1]]),"results","profileList_pos_copy"),envir=as.environment(".GlobalEnv"));	
				load(file=file.path(logfile[[1]],"results","pattern_pos_IS"),envir=as.environment(".GlobalEnv"));
				pattern_pos<-pattern_pos_IS;rm(pattern_pos_IS,envir=as.environment(".GlobalEnv"));
				patt_pos_ID<-rep("",length(pattern_pos))
				patt_pos_add<-rep("",length(pattern_pos))
				for(i in 1:length(pattern_pos)){
					patt_pos_ID[i]<-strsplit(names(pattern_pos[i]),"_",fixed=TRUE)[[1]][1]
					patt_pos_add[i]<-strsplit(names(pattern_pos[i]),"_",fixed=TRUE)[[1]][2]
				}
				load(file=file.path(logfile[[1]],"results","patternRT_pos_IS"),envir=as.environment(".GlobalEnv"));
				pattern_RT_pos<-patternRT_pos_IS;rm(patternRT_pos_IS,envir=as.environment(".GlobalEnv"));
				load(file=file.path(logfile[[1]],"results","patternDelRT_pos_IS"),envir=as.environment(".GlobalEnv"));
				pattern_delRT_pos<-patternDelRT_pos_IS;rm(patternDelRT_pos_IS,envir=as.environment(".GlobalEnv"));
				load(file=file.path(logfile$project_folder,"results","screening","res_IS_pos_screen"))
				res_pos_screen<-res_IS_pos_screen;rm(res_IS_pos_screen);
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				compound_table<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
				# For IS without Conz.:
				updateSelectInput(session,inputId="selec_pos_x",label="x axis",choices=c("m/z","RT","Intensity","Date&time","Type","Place"),selected = "m/z")
				updateSelectInput(session,inputId="selec_pos_y",label="y axis",choices=c("m/z","RT","Intensity","Date&time","Type","Place"),selected = "RT")	
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
			output$plot_pattern_pos <- renderPlot({
				s<-input$Table_screening_pos_row_last_clicked
				if (length(s) & isolate(input$screen_pos_summarize=="yes")) {	
					use_comp<-(
						(patt_pos_ID==as.character(results_screen_pos[s,1])) & (patt_pos_add==as.character(results_screen_pos[s,3]))
					)
					if(sum(use_comp)>1){stop("Report this issue for a debug on server_obs_screening.r #1")}															
					pattern_sel<-pattern_pos[use_comp][[1]]
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
				}else{
					plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(.5,.5,labels="No compound selected or adducts collapsed",col="red",cex=1.6)				
				}
			})
			# make Table over samples
			output$Table_screening_selected_pos<-DT::renderDataTable({
				s<-input$Table_screening_pos_row_last_clicked
				if (length(s) & isolate(input$screen_pos_summarize=="yes")) {
					use_comp<-(
						(patt_pos_ID==as.character(results_screen_pos[s,1])) & (patt_pos_add==as.character(results_screen_pos[s,3]))
					)
					if(sum(use_comp)>1){stop("Report this issue for a debug on server_obs_screening.r #1")}	
					res_pos_screen_sel<-res_pos_screen[use_comp][[1]]
					which_where<-c();which_peaks<-c();sample_type<-c();score_1<-c();score_2<-c();delppm<-c();delRT<-c();inte<-c();
					with_peaks<-c();#with_file<-c();with_s<-c();
					IDs<-as.numeric(measurements[,1]) 
					if(length(res_pos_screen_sel)>0){
						for(i in 1:length(res_pos_screen_sel)){
							if(length(res_pos_screen_sel[[i]])>0){
								for(j in 1:length(res_pos_screen_sel[[i]])){
									which_where<-c(which_where,measurements[IDs==i,1]);
									sample_type<-c(sample_type,measurements[IDs==i,3]);
									which_peaks<-c(which_peaks,paste(res_pos_screen_sel[[i]][[j]]$Peaks[,1],collapse=", "))
									score_1<-c(score_1,round(res_pos_screen_sel[[i]][[j]]$score_1,digits=2));
									score_2<-c(score_2,round(res_pos_screen_sel[[i]][[j]]$score_2,digits=2));
									delppm<-c(delppm,paste(as.character(round(res_pos_screen_sel[[i]][[j]][[4]],digits=2)),collapse=", "));
									#delRT<-c(delRT,paste(as.character(round(res_pos_screen_sel[[i]][[j]][[4]],digits=2)),collapse=", "));
									with_peaks<-c(with_peaks,paste(as.character(
										profileList_pos[[2]][ # insert peak IDs - original $Peaks refer to an entry only!
											res_pos_screen_sel[[i]][[j]]$Peaks[,2],4
										]
									),collapse=", "));
									delRT<-c(delRT,
										paste(as.character(round(res_pos_screen_sel[[i]][[j]]$RT,digits=2)),collapse=", ")
									)
									inte<-c(inte,
										paste(as.character(round(log10(res_pos_screen_sel[[i]][[j]]$Intensity),digits=2)),collapse=", ")									
									)
									#with_file<-c(with_file,i)
									#with_s<-c(with_s,s)
								}
							}
						}
					}
					DT::datatable(
						as.data.frame(cbind(which_where,sample_type,which_peaks,score_1,score_2,delppm,delRT,inte,with_peaks),
						#with_file,with_s),
						row.names = NULL,stringsAsFactors=FALSE),
						rownames = FALSE, colnames=c("File ID","File type","Pattern matches","Score > LOD","Score < LOD",
							"m/z deviation (ppm)","RT","log Intensity","Peak IDs")#,"m","i")
					)
				}else{
					DT::datatable(as.data.frame(cbind("")),selection = 'single',rownames=FALSE,colnames="No compound selected or adducts collapsed")
				}
			},server = TRUE)
			# initialize intensity range for selected internal standard
			observe({ 
				s<-input$Table_screening_pos_row_last_clicked
				if (length(s) & isolate(input$screen_pos_summarize=="yes") & isolate(input$Pos_compound_select=="Internal standards")) {
					updateNumericInput(session, "screen_int_pos_low", "Lower bound", value = 5,step=0.1)   
					updateNumericInput(session, "screen_int_pos_up", "Upper bound", value = 6,step=0.1) 
				}	
			})
			# export intensity range for selected internal standard
			observe({
				input$save_int_pos
				s<-isolate(input$Table_screening_pos_row_last_clicked)
				if (isolate(input$save_int_pos) & length(s) & isolate(input$screen_pos_summarize=="yes") & isolate(input$Pos_compound_select=="Internal standards")) {
				
				
				}
			})
			output$plot_selec_dist_pos <- renderPlot({
				s<-input$Table_screening_pos_row_last_clicked
				input$selec_pos_log_rat
				input$selec_pos_x
				input$selec_pos_y
				if (length(s) & isolate(input$screen_pos_summarize=="yes")) {			
					use_comp<-(
						(patt_pos_ID==as.character(results_screen_pos[s,1])) & (patt_pos_add==as.character(results_screen_pos[s,3]))
					)
					if(sum(use_comp)>1){stop("Report this issue for a debug on server_obs_screening.r #1")}	
					pattern_sel<-pattern_pos[use_comp][[1]]
					res_pos_screen_sel<-res_pos_screen[use_comp][[1]]
					cut_score<-as.numeric(logfile$parameters$IS_w1)	
					IDs<-as.numeric(measurements[,1])
					# extract relevant data for the compound - adduct ##########################
					if(length(res_pos_screen_sel)>0){
						mass<-c();inte<-c();RT<-c();cutit<-c();atdate<-c();attime<-c();placed<-c();typed<-c();
						for(i in 1:length(res_pos_screen_sel)){
							if(length(res_pos_screen_sel[[i]])>0){
								for(j in 1:length(res_pos_screen_sel[[i]])){
									mass<-c(mass,res_pos_screen_sel[[i]][[j]][[7]])
									inte<-c(inte,res_pos_screen_sel[[i]][[j]][[8]])
									RT<-c(RT,res_pos_screen_sel[[i]][[j]][[9]])
									local_score<-0
									if(!is.na(res_pos_screen_sel[[i]][[j]]$score_1)){
										local_score<-(res_pos_screen_sel[[i]][[j]]$score_1)
									}
									if( (local_score>=1) || (is.na(res_pos_screen_sel[[i]][[j]]$score_1)) ){
										if(!is.na(res_pos_screen_sel[[i]][[j]]$score_2)){
											local_score<-(local_score+res_pos_screen_sel[[i]][[j]]$score_2)
										}
									}									
									if(local_score>=cut_score){
										cutit<-c(cutit,1);
									}else{
										cutit<-c(cutit,0);
									}
									lengi<-length(res_pos_screen_sel[[i]][[j]][[7]])
									placed<-c(placed,rep(measurements[IDs==i,5],lengi))
									typed<-c(typed,rep(measurements[IDs==i,3],lengi))
									atdate<-c(atdate,rep(measurements[IDs==i,6],lengi))
									attime<-c(attime,rep(measurements[IDs==i,7],lengi))								
								}
							}
						}
						timed<-as.POSIXct(paste(atdate,attime,"CET",sep=" "))
						timed2<-pretty(timed)
						timelimit<-c(min(timed),max(timed))		
						placed<-as.factor(placed)
						typed<-as.factor(typed)	
						int_pos_sel_low<-(10^input$screen_int_pos_low)
						int_pos_sel_up<-(10^input$screen_int_pos_up)
						# plot ######################################################################
						if(isolate(input$selec_pos_log_rat=="yes")){
							inte<-log10(inte)
							if(int_pos_sel_low!=FALSE){int_pos_sel_low<-log10(int_pos_sel_low)}
							if(int_pos_sel_up!=FALSE){int_pos_sel_up<-log10(int_pos_sel_up)}						
						}
						if(isolate(input$selec_pos_x=="m/z")){sel_pos_x<-mass;xlab_sel<-"m/z"}
						if(isolate(input$selec_pos_x=="RT")){sel_pos_x<-RT;xlab_sel<-"RT"}	
						if(isolate(input$selec_pos_x=="Intensity")){sel_pos_x<-inte;
							if(isolate(input$selec_pos_log_rat=="yes")){xlab_sel<-"log Intensity"}else{xlab_sel<-"Intensity"}
						}
						if(isolate(input$selec_pos_x=="Date&time")){sel_pos_x<-timed;xlab_sel<-"Time"}
						if(isolate(input$selec_pos_x=="Type")){sel_pos_x<-typed;xlab_sel<-"Type"}
						if(isolate(input$selec_pos_x=="Place")){sel_pos_x<-typed;xlab_sel<-"Place"}
						
						if(isolate(input$selec_pos_y=="m/z")){sel_pos_y<-mass;ylab_sel<-"m/z"}
						if(isolate(input$selec_pos_y=="RT")){sel_pos_y<-RT;ylab_sel<-"RT"}	
						if(isolate(input$selec_pos_y=="Intensity")){sel_pos_y<-inte;
							if(isolate(input$selec_pos_log_rat=="yes")){ylab_sel<-"log Intensity"}else{ylab_sel<-"Intensity"}
						}
						if(isolate(input$selec_pos_y=="Date&time")){sel_pos_y<-timed;ylab_sel<-"Time"}
						if(isolate(input$selec_pos_y=="Type")){sel_pos_y<-typed;ylab_sel<-"Type"}
						if(isolate(input$selec_pos_y=="Place")){sel_pos_y<-typed;ylab_sel<-"Place"}												
						par(mar=c(4,4,.8,.5))
						plot(sel_pos_x,sel_pos_y,pch=19,cex=.7,
							xlab=xlab_sel,ylab=ylab_sel,col="darkgrey"
						)
						if( 
							(isolate(input$selec_pos_x=="Intensity")||isolate(input$selec_pos_y=="Intensity")) & 
							int_pos_sel_low!=FALSE &
							int_pos_sel_up!=FALSE &
							isolate(input$selec_pos_x!="Type") &
							isolate(input$selec_pos_x!="Place") &
							isolate(input$Pos_compound_select=="Internal standards")
						){
							if(isolate(input$selec_pos_x=="Intensity")){
								rec_x1<-int_pos_sel_low
								rec_x2<-int_pos_sel_up						
							}else{
								if(input$selec_pos_x=="Date&time"){
									rec_x1<-min(timed)
									rec_x2<-max(timed)
								}else{
									rec_x1<-0
									rec_x2<-2*max(sel_pos_x)
								}
							}
							if(isolate(input$selec_pos_y=="Intensity")){
								rec_y1<-int_pos_sel_low
								rec_y2<-int_pos_sel_up						
							}else{
								if(input$selec_pos_y=="Date&time"){
									rec_y1<-min(timed)
									rec_y2<-max(timed)
								}else{
									rec_y1<-0
									rec_y2<-2*max(sel_pos_x)
								}
							}
							rect(rec_x1,rec_y1,rec_x2,rec_y2,border = NA,col="orange3")					
						}
						if(isolate(input$selec_pos_x=="m/z")){abline(v=pattern_sel[,1],col="red")}
						if(isolate(input$selec_pos_y=="m/z")){abline(h=pattern_sel[,1],col="red")}
						points(sel_pos_x[cutit==0],sel_pos_y[cutit==0],pch=19,cex=.7,col="darkgrey");
						points(sel_pos_x[cutit==1],sel_pos_y[cutit==1],pch=19,cex=.7,col="black");					
						box();
					}else{
						plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(.5,.5,labels="No screening matches for this compound",col="red",cex=1.6)
					}
				}else{
					plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(.5,.5,labels="No compound selected or adducts collapsed",col="red",cex=1.6)
				}
			},width = "auto", height = "auto")
			output$plot_pattern_distrib_pos <- renderPlot({
				if(length(screen_dev_pos)>0){
					use_x<-input$Summ_pos_x
					use_y<-input$Summ_pos_y
					if(use_x!=use_y){
						par(mar=c(4,4,.5,.5))
						plot(
							screen_dev_pos[,colnames(screen_dev_pos)==use_x],
							screen_dev_pos[,colnames(screen_dev_pos)==use_y],
							pch=19,cex=.3,xlab=use_x,ylab=use_y,col="lightgrey"
						)
						points(
							screen_dev_pos[screen_dev_pos[,6]==1,colnames(screen_dev_pos)==use_x],
							screen_dev_pos[screen_dev_pos[,6]==1,colnames(screen_dev_pos)==use_y],
							pch=19,cex=.4,xlab=use_x,ylab=use_y,col="black"
						)					
						if(use_x=="m/z deviation [ppm]" | use_x=="RT deviation"){abline(v=0,col="red")}
						if(use_y=="m/z deviation [ppm]" | use_y=="RT deviation"){abline(h=0,col="red")}
						plot.window(xlim=c(0,1),ylim=c(0,1))
						legend(0.9,1,title="Cutoff score",legend=c("below","above"),fill=c("lightgrey","black"),border=c("lightgrey","black"))
					}else{
						par(mar=c(4,4,.5,.5))
						plot.new()
						plot.window(xlim=c(min(screen_dev_pos[,colnames(screen_dev_pos)==use_x]),max(screen_dev_pos[,colnames(screen_dev_pos)==use_x])),ylim=c(0,3.3))
						boxplot(screen_dev_pos[screen_dev_pos[,6]==1,colnames(screen_dev_pos)==use_x],
						horizontal=TRUE,xlab=use_x,width=1.3,at=2,add=TRUE)				
						boxplot(screen_dev_pos[screen_dev_pos[,6]==0,colnames(screen_dev_pos)==use_x],
						horizontal=TRUE,xlab=use_x,width=1.3,at=1,add=TRUE,col="grey")			
						plot.window(xlim=c(0,1),ylim=c(0,1))
						legend(0.8,1,title="Cutoff score",legend=c("above","below"),fill=c("white","lightgrey","black"),border=c("black","lightgrey"))
					}
				}
			})			
			output$plot_aboveBlank_pos <- renderPlot({
				if(any(rat_sam_blank_pos>0)){
					par(mar=c(4,4,.5,.5))
					if(input$screen_pos_log_rat=="yes"){
						boxplot(log10(rat_sam_blank_pos[rat_sam_blank_pos>0]),
						horizontal=TRUE,xlab="log10 intensity ratio",width=1.3)
					}else{
						boxplot((rat_sam_blank_pos[rat_sam_blank_pos>0]),
						horizontal=TRUE,xlab="Intensity ratio",width=1.3)
					}
				}else{
					plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(.5,.5,labels="No ratios available",col="red",cex=1.6)
				}
			},width = "auto", height = 250)
			
		} # if(found_table)
	} # if init$a
})  
##############################################################################

##############################################################################
# NEGATIVE IONIZATION ########################################################
observe({ 
	input$Neg_compound_select  
	input$screen_neg_summarize
	init$b
	if(isolate(init$a)=="TRUE"){
		found_table<-FALSE
		if( isolate(input$Neg_compound_select=="Target compounds") ){
			cat("\n Looking at negative targets_selec")
			if( 				
				(file.exists(file=file.path(logfile$project_folder,"results","screening","results_screen_target_neg"))) #&
				#(logfile$workflow[names(logfile$workflow)=="target_subtr"]!="yes")
			){
				load(file=file.path(logfile$project_folder,"results","screening","results_screen_target_neg"))
				screen_dev_neg<-results_screen_target_neg[[3]] # contains sample vs. blank intensity ratios
				screen_dev_neg[,2]<-log10(screen_dev_neg[,2])
				rat_sam_blank_neg<-results_screen_target_neg[[1]][,10,drop=FALSE]
				if( isolate(input$screen_neg_summarize=="yes") ){
					results_screen_neg<-results_screen_target_neg[[1]]
				}else{
					results_screen_neg<-results_screen_target_neg[[2]]
				}
				table_screening_neg<-DT::datatable(results_screen_neg, escape = FALSE,selection = 'single') %>% 
					formatStyle('Max. sample score',background = styleColorBar(c(0,1), 'lightgreen'),backgroundPosition = 'right')
				output$Table_screening_neg <- DT::renderDataTable({table_screening_neg},server = TRUE)
				rm(results_screen_target_neg)
				found_table<-TRUE
				load(file=file.path(as.character(logfile[[1]]),"results","profileList_neg_copy"),envir=as.environment(".GlobalEnv"));	
				load(file=file.path(logfile[[1]],"results","pattern_neg_target"),envir=as.environment(".GlobalEnv"));
				pattern_neg<-pattern_neg_target;rm(pattern_neg_target,envir=as.environment(".GlobalEnv"));
				patt_neg_ID<-rep("",length(pattern_neg))
				patt_neg_add<-rep("",length(pattern_neg))
				for(i in 1:length(pattern_neg)){
					patt_neg_ID[i]<-strsplit(names(pattern_neg[i]),"_",fixed=TRUE)[[1]][1]
					patt_neg_add[i]<-strsplit(names(pattern_neg[i]),"_",fixed=TRUE)[[1]][2]
				}
				load(file=file.path(logfile[[1]],"results","patternRT_neg_target"),envir=as.environment(".GlobalEnv"));
				pattern_RT_neg<-patternRT_neg_target;rm(patternRT_neg_target,envir=as.environment(".GlobalEnv"));
				load(file=file.path(logfile[[1]],"results","patternDelRT_neg_target"),envir=as.environment(".GlobalEnv"));
				pattern_delRT_neg<-patternDelRT_neg_target;rm(patternDelRT_neg_target,envir=as.environment(".GlobalEnv"));
				load(file=file.path(logfile$project_folder,"results","screening","res_target_neg_screen"))
				res_neg_screen<-res_target_neg_screen;rm(res_target_neg_screen);
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				compound_table<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");	
				# For IS without Conz.:
				updateSelectInput(session,inputId="selec_neg_x",label="x axis",choices=c("m/z","RT","Intensity","Date&time","Type","Place","Conz."),selected = "m/z")
				updateSelectInput(session,inputId="selec_neg_y",label="y axis",choices=c("m/z","RT","Intensity","Date&time","Type","Place","Conz."),selected = "RT")
			}else{	
				output$Table_screening_neg <- DT::renderDataTable({
					DT::datatable(as.data.frame(cbind("")),selection = 'single',rownames=FALSE,colnames="No target screening results available")
				},server = TRUE)	
				found_table<-FALSE				
			}
		}
		if( isolate(input$Neg_compound_select=="Internal standards") ){
			cat("\n Looking at negative standards_selec")
			if(
				(file.exists(file=file.path(logfile$project_folder,"results","screening","results_screen_IS_neg"))) #&
				#(logfile$workflow[names(logfile$workflow)=="IS_subtr"]!="yes")			
			){
				load(file=file.path(logfile$project_folder,"results","screening","results_screen_IS_neg"))
				screen_dev_neg<-results_screen_IS_neg[[3]] # contains sample vs. blank intensity ratios		
				screen_dev_neg[,2]<-log10(screen_dev_neg[,2])
				rat_sam_blank_neg<-results_screen_IS_neg[[1]][,10,drop=FALSE]				
				if( isolate(input$screen_neg_summarize=="yes") ){
					results_screen_neg<-results_screen_IS_neg[[1]]
				}else{
					results_screen_neg<-results_screen_IS_neg[[2]]
				}
				table_screening_neg<-DT::datatable(results_screen_neg, escape = FALSE,selection = 'single') %>% 
					formatStyle('Max. sample score',background = styleColorBar(c(0,1), 'lightgreen'),backgroundPosition = 'right')
				output$Table_screening_neg <- DT::renderDataTable({table_screening_neg},server = TRUE)
				rm(results_screen_IS_neg)	
				found_table<-TRUE
				load(file=file.path(as.character(logfile[[1]]),"results","profileList_neg_copy"),envir=as.environment(".GlobalEnv"));	
				load(file=file.path(logfile[[1]],"results","pattern_neg_IS"),envir=as.environment(".GlobalEnv"));
				pattern_neg<-pattern_neg_IS;rm(pattern_neg_IS,envir=as.environment(".GlobalEnv"));
				patt_neg_ID<-rep("",length(pattern_neg))
				patt_neg_add<-rep("",length(pattern_neg))
				for(i in 1:length(pattern_neg)){
					patt_neg_ID[i]<-strsplit(names(pattern_neg[i]),"_",fixed=TRUE)[[1]][1]
					patt_neg_add[i]<-strsplit(names(pattern_neg[i]),"_",fixed=TRUE)[[1]][2]
				}
				load(file=file.path(logfile[[1]],"results","patternRT_neg_IS"),envir=as.environment(".GlobalEnv"));
				pattern_RT_neg<-patternRT_neg_IS;rm(patternRT_neg_IS,envir=as.environment(".GlobalEnv"));
				load(file=file.path(logfile[[1]],"results","patternDelRT_neg_IS"),envir=as.environment(".GlobalEnv"));
				pattern_delRT_neg<-patternDelRT_neg_IS;rm(patternDelRT_neg_IS,envir=as.environment(".GlobalEnv"));
				load(file=file.path(logfile$project_folder,"results","screening","res_IS_neg_screen"))
				res_neg_screen<-res_IS_neg_screen;rm(res_IS_neg_screen);
				measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				compound_table<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
				# For IS without Conz.:
				updateSelectInput(session,inputId="selec_neg_x",label="x axis",choices=c("m/z","RT","Intensity","Date&time","Type","Place"),selected = "m/z")
				updateSelectInput(session,inputId="selec_neg_y",label="y axis",choices=c("m/z","RT","Intensity","Date&time","Type","Place"),selected = "RT")	
			}else{	
				output$Table_screening_neg <- DT::renderDataTable({
					DT::datatable(as.data.frame(cbind("")),selection = 'single',rownames=FALSE,colnames="No internal standard screening results available")
				},server = TRUE)	
				found_table<-FALSE				
			}
		}
		if(found_table){
			# name selected compound & adduct
			named_compound <- renderText({
				s<-input$Table_screening_neg_row_last_clicked
				if (length(s)) {	paste(results_screen_neg[s,2],results_screen_neg[s,3],
				compound_table[compound_table[,1]==results_screen_neg[s,1],3]
				,sep=", ")	}
			})
			output$screening_details_comp_neg <-named_compound
			output$screening_details_comp_neg2<-named_compound
			output$screening_details_comp_neg3<-named_compound
			# plot pattern & matches
			output$plot_pattern_neg <- renderPlot({
				s<-input$Table_screening_neg_row_last_clicked
				if (length(s) & isolate(input$screen_neg_summarize=="yes")) {	
					use_comp<-(
						(patt_neg_ID==as.character(results_screen_neg[s,1])) & (patt_neg_add==as.character(results_screen_neg[s,3]))
					)
					if(sum(use_comp)>1){stop("Report this issue for a debug on server_obs_screening.r #1")}	
					pattern_sel<-pattern_neg[use_comp][[1]]
					res_neg_screen_sel<-res_neg_screen[use_comp][[1]]
					int_tol<-as.numeric(logfile$parameters$IS_inttol)
					ylim_up<-(int_tol+100+(int_tol*0.1))
					plot(pattern_sel[,1],pattern_sel[,2],type="h",lwd=3,col="red",xlab="m/z",ylab="Rescaled intensity",ylim=c(0,ylim_up))
					if(length(res_neg_screen_sel)>0){
						for(i in 1:length(res_neg_screen_sel)){
							if(length(res_neg_screen_sel[[i]])>0){
								for(j in 1:length(res_neg_screen_sel[[i]])){
									found_matches<-res_neg_screen_sel[[i]][[j]]$Peaks
									x<-profileList_neg[[2]][found_matches[,2],1]
									y<-(profileList_neg[[2]][found_matches[,2],2]*res_neg_screen_sel[[i]][[j]][[6]])
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
				}else{
					plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(.5,.5,labels="No compound selected or adducts collapsed",col="red",cex=1.6)				
				}
			})
			# make Table over samples
			output$Table_screening_selected_neg<-DT::renderDataTable({
				s<-input$Table_screening_neg_row_last_clicked
				if(length(s) & isolate(input$screen_neg_summarize=="yes")){
					use_comp<-(
						(patt_neg_ID==as.character(results_screen_neg[s,1])) & (patt_neg_add==as.character(results_screen_neg[s,3]))
					)
					if(sum(use_comp)>1){stop("Report this issue for a debug on server_obs_screening.r #1")}							
					res_neg_screen_sel<-res_neg_screen[use_comp][[1]]
					which_where<-c();which_peaks<-c();sample_type<-c();score_1<-c();score_2<-c();delppm<-c();delRT<-c();inte<-c();
					with_peaks<-c();#with_file<-c();with_s<-c();
					IDs<-as.numeric(measurements[,1]) 
					if(length(res_neg_screen_sel)>0){
						for(i in 1:length(res_neg_screen_sel)){
							if(length(res_neg_screen_sel[[i]])>0){
								for(j in 1:length(res_neg_screen_sel[[i]])){
									which_where<-c(which_where,measurements[IDs==i,1]);
									sample_type<-c(sample_type,measurements[IDs==i,3]);
									which_peaks<-c(which_peaks,paste(res_neg_screen_sel[[i]][[j]]$Peaks[,1],collapse=", "))
									score_1<-c(score_1,round(res_neg_screen_sel[[i]][[j]]$score_1,digits=2));
									score_2<-c(score_2,round(res_neg_screen_sel[[i]][[j]]$score_2,digits=2));
									delppm<-c(delppm,paste(as.character(round(res_neg_screen_sel[[i]][[j]][[4]],digits=2)),collapse=", "));
									#delRT<-c(delRT,paste(as.character(round(res_neg_screen_sel[[i]][[j]][[4]],digits=2)),collapse=", "));
									with_peaks<-c(with_peaks,paste(as.character(
										profileList_neg[[2]][ # insert peak IDs - original $Peaks refer to an entry only!
											res_neg_screen_sel[[i]][[j]]$Peaks[,2],4
										]
									),collapse=", "));
									delRT<-c(delRT,
										paste(as.character(round(res_neg_screen_sel[[i]][[j]]$RT,digits=2)),collapse=", ")
									)
									inte<-c(inte,
										paste(as.character(round(log10(res_neg_screen_sel[[i]][[j]]$Intensity),digits=2)),collapse=", ")									
									)
									#with_file<-c(with_file,i)
									#with_s<-c(with_s,s)
								}
							}
						}
					}
					DT::datatable(
						as.data.frame(cbind(which_where,sample_type,which_peaks,score_1,score_2,delppm,delRT,inte,with_peaks),
						#with_file,with_s),
						row.names = NULL,stringsAsFactors=FALSE),
						rownames = FALSE, colnames=c("File ID","File type","Pattern matches","Score > LOD","Score < LOD",
							"m/z deviation (ppm)","RT","log Intensity","Peak IDs")#,"m","i")
					)
				}else{
					DT::datatable(as.data.frame(cbind("")),selection = 'single',rownames=FALSE,colnames="No compound selected or adducts collapsed")
				}
			},server = TRUE)

			# initialize intensity range for selected internal standard
			observe({ 
				s<-input$Table_screening_neg_row_last_clicked
				if (length(s) & isolate(input$screen_neg_summarize=="yes") & isolate(input$Neg_compound_select=="Internal standards")) {
					updateNumericInput(session, "screen_int_neg_low", "Lower bound", value = 5,step=0.1)   
					updateNumericInput(session, "screen_int_neg_up", "Upper bound", value = 6,step=0.1) 
# BAUSTELLE					
				}	
			})
			# export intensity range for selected internal standard
			observe({
				input$save_int_neg
				s<-isolate(input$Table_screening_neg_row_last_clicked)
				if (isolate(input$save_int_neg) & length(s) & isolate(input$screen_neg_summarize=="yes") & isolate(input$Neg_compound_select=="Internal standards")) {
				
				
				}
			})
			output$plot_selec_dist_neg <- renderPlot({
				s<-input$Table_screening_neg_row_last_clicked
				input$selec_neg_log_rat
				input$selec_neg_x
				input$selec_neg_y
				if (length(s) & isolate(input$screen_neg_summarize=="yes")) {			
					use_comp<-(
						(patt_neg_ID==as.character(results_screen_neg[s,1])) & (patt_neg_add==as.character(results_screen_neg[s,3]))
					)
					if(sum(use_comp)>1){stop("Report this issue for a debug on server_obs_screening.r #1")}	
					pattern_sel<-pattern_neg[use_comp][[1]]
					res_neg_screen_sel<-res_neg_screen[use_comp][[1]]
					cut_score<-as.numeric(logfile$parameters$IS_w1)	
					IDs<-as.numeric(measurements[,1])
					# extract relevant data for the compound - adduct ##########################
					if(length(res_neg_screen_sel)>0){
						mass<-c();inte<-c();RT<-c();cutit<-c();atdate<-c();attime<-c();placed<-c();typed<-c();
						for(i in 1:length(res_neg_screen_sel)){
							if(length(res_neg_screen_sel[[i]])>0){
								for(j in 1:length(res_neg_screen_sel[[i]])){
									mass<-c(mass,res_neg_screen_sel[[i]][[j]][[7]])
									inte<-c(inte,res_neg_screen_sel[[i]][[j]][[8]])
									RT<-c(RT,res_neg_screen_sel[[i]][[j]][[9]])
									local_score<-0
									if(!is.na(res_neg_screen_sel[[i]][[j]]$score_1)){
										local_score<-(res_neg_screen_sel[[i]][[j]]$score_1)
									}
									if( (local_score>=1) || (is.na(res_neg_screen_sel[[i]][[j]]$score_1)) ){
										if(!is.na(res_neg_screen_sel[[i]][[j]]$score_2)){
											local_score<-(local_score+res_neg_screen_sel[[i]][[j]]$score_2)
										}
									}									
									if(local_score>=cut_score){
										cutit<-c(cutit,1);
									}else{
										cutit<-c(cutit,0);
									}
									lengi<-length(res_neg_screen_sel[[i]][[j]][[7]])
									placed<-c(placed,rep(measurements[IDs==i,5],lengi))
									typed<-c(typed,rep(measurements[IDs==i,3],lengi))
									atdate<-c(atdate,rep(measurements[IDs==i,6],lengi))
									attime<-c(attime,rep(measurements[IDs==i,7],lengi))								
								}
							}
						}
						timed<-as.POSIXct(paste(atdate,attime,"CET",sep=" "))
						timed2<-pretty(timed)
						timelimit<-c(min(timed),max(timed))		
						placed<-as.factor(placed)
						typed<-as.factor(typed)	
						int_neg_sel_low<-(10^input$screen_int_neg_low)
						int_neg_sel_up<-(10^input$screen_int_neg_up)
						# plot ######################################################################
						if(isolate(input$selec_neg_log_rat=="yes")){
							inte<-log10(inte)
							if(int_neg_sel_low!=FALSE){int_neg_sel_low<-log10(int_neg_sel_low)}
							if(int_neg_sel_up!=FALSE){int_neg_sel_up<-log10(int_neg_sel_up)}						
						}
						if(isolate(input$selec_neg_x=="m/z")){sel_neg_x<-mass;xlab_sel<-"m/z"}
						if(isolate(input$selec_neg_x=="RT")){sel_neg_x<-RT;xlab_sel<-"RT"}	
						if(isolate(input$selec_neg_x=="Intensity")){sel_neg_x<-inte;
							if(isolate(input$selec_neg_log_rat=="yes")){xlab_sel<-"log Intensity"}else{xlab_sel<-"Intensity"}
						}
						if(isolate(input$selec_neg_x=="Date&time")){sel_neg_x<-timed;xlab_sel<-"Time"}
						if(isolate(input$selec_neg_x=="Type")){sel_neg_x<-typed;xlab_sel<-"Type"}
						if(isolate(input$selec_neg_x=="Place")){sel_neg_x<-typed;xlab_sel<-"Place"}
						
						if(isolate(input$selec_neg_y=="m/z")){sel_neg_y<-mass;ylab_sel<-"m/z"}
						if(isolate(input$selec_neg_y=="RT")){sel_neg_y<-RT;ylab_sel<-"RT"}	
						if(isolate(input$selec_neg_y=="Intensity")){sel_neg_y<-inte;
							if(isolate(input$selec_neg_log_rat=="yes")){ylab_sel<-"log Intensity"}else{ylab_sel<-"Intensity"}
						}
						if(isolate(input$selec_neg_y=="Date&time")){sel_neg_y<-timed;ylab_sel<-"Time"}
						if(isolate(input$selec_neg_y=="Type")){sel_neg_y<-typed;ylab_sel<-"Type"}
						if(isolate(input$selec_neg_y=="Place")){sel_neg_y<-typed;ylab_sel<-"Place"}												
						par(mar=c(4,4,.8,.5))
						plot(sel_neg_x,sel_neg_y,pch=19,cex=.7,
							xlab=xlab_sel,ylab=ylab_sel,col="darkgrey"
						)
						if( 
							(isolate(input$selec_neg_x=="Intensity")||isolate(input$selec_neg_y=="Intensity")) & 
							int_neg_sel_low!=FALSE &
							int_neg_sel_up!=FALSE &
							isolate(input$selec_neg_x!="Type") &
							isolate(input$selec_neg_x!="Place") &
							isolate(input$Neg_compound_select=="Internal standards")
						){
							if(isolate(input$selec_neg_x=="Intensity")){
								rec_x1<-int_neg_sel_low
								rec_x2<-int_neg_sel_up						
							}else{
								if(input$selec_neg_x=="Date&time"){
									rec_x1<-min(timed)
									rec_x2<-max(timed)
								}else{
									rec_x1<-0
									rec_x2<-2*max(sel_neg_x)
								}
							}
							if(isolate(input$selec_neg_y=="Intensity")){
								rec_y1<-int_neg_sel_low
								rec_y2<-int_neg_sel_up						
							}else{
								if(input$selec_neg_y=="Date&time"){
									rec_y1<-min(timed)
									rec_y2<-max(timed)
								}else{
									rec_y1<-0
									rec_y2<-2*max(sel_neg_x)
								}
							}
							rect(rec_x1,rec_y1,rec_x2,rec_y2,border = NA,col="orange3")					
						}
						if(isolate(input$selec_neg_x=="m/z")){abline(v=pattern_sel[,1],col="red")}
						if(isolate(input$selec_neg_y=="m/z")){abline(h=pattern_sel[,1],col="red")}
						points(sel_neg_x[cutit==0],sel_neg_y[cutit==0],pch=19,cex=.7,col="darkgrey");
						points(sel_neg_x[cutit==1],sel_neg_y[cutit==1],pch=19,cex=.7,col="black");					
						box();
					}else{
						plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(.5,.5,labels="No screening matches for this compound",col="red",cex=1.6)
					}
				}else{
					plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(.5,.5,labels="No compound selected or adducts collapsed",col="red",cex=1.6)
				}
			},width = "auto", height = "auto")

			output$plot_pattern_distrib_neg <- renderPlot({
				if(length(screen_dev_neg)>0){
					use_x<-input$Summ_neg_x
					use_y<-input$Summ_neg_y
					if(use_x!=use_y){
						par(mar=c(4,4,.5,.5))
						plot(
							screen_dev_neg[,colnames(screen_dev_neg)==use_x],
							screen_dev_neg[,colnames(screen_dev_neg)==use_y],
							pch=19,cex=.3,xlab=use_x,ylab=use_y,col="lightgrey"
						)
						points(
							screen_dev_neg[screen_dev_neg[,6]==1,colnames(screen_dev_neg)==use_x],
							screen_dev_neg[screen_dev_neg[,6]==1,colnames(screen_dev_neg)==use_y],
							pch=19,cex=.4,xlab=use_x,ylab=use_y,col="black"
						)					
						if(use_x=="m/z deviation [ppm]" | use_x=="RT deviation"){abline(v=0,col="red")}
						if(use_y=="m/z deviation [ppm]" | use_y=="RT deviation"){abline(h=0,col="red")}
						plot.window(xlim=c(0,1),ylim=c(0,1))
						legend(0.9,1,title="Cutoff score",legend=c("below","above"),fill=c("lightgrey","black"),border=c("lightgrey","black"))
					}else{
						par(mar=c(4,4,.5,.5))
						plot.new()
						plot.window(xlim=c(min(screen_dev_neg[,colnames(screen_dev_neg)==use_x]),max(screen_dev_neg[,colnames(screen_dev_neg)==use_x])),ylim=c(0,3.3))
						boxplot(screen_dev_neg[screen_dev_neg[,6]==1,colnames(screen_dev_neg)==use_x],
						horizontal=TRUE,xlab=use_x,width=1.3,at=2,add=TRUE)				
						boxplot(screen_dev_neg[screen_dev_neg[,6]==0,colnames(screen_dev_neg)==use_x],
						horizontal=TRUE,xlab=use_x,width=1.3,at=1,add=TRUE,col="grey")			
						plot.window(xlim=c(0,1),ylim=c(0,1))
						legend(0.8,1,title="Cutoff score",legend=c("above","below"),fill=c("white","lightgrey","black"),border=c("black","lightgrey"))
					}
				}
			})			
			output$plot_aboveBlank_neg <- renderPlot({
				if(any(rat_sam_blank_neg>0)){
					par(mar=c(4,4,.5,.5))
					if(input$screen_neg_log_rat=="yes"){
						boxplot(log10(rat_sam_blank_neg[rat_sam_blank_neg>0]),
						horizontal=TRUE,xlab="log10 intensity ratio",width=1.3)
					}else{
						boxplot((rat_sam_blank_neg[rat_sam_blank_neg>0]),
						horizontal=TRUE,xlab="Intensity ratio",width=1.3)
					}
				}else{
					plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(.5,.5,labels="No ratios available",col="red",cex=1.6)
				}
			},width = "auto", height = 250)
			
		} # if(found_table)
	} # if init$a
})  
##############################################################################

if(any(ls()=="logfile")){stop("\n illegal logfile detected #2 in server_obs_screening.r!")}




