##############################################################################
# observe components outputs #################################################
##############################################################################
if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_obs_screening.r!")}
ee	<-	reactiveValues() # reactive value ...
ee$entry<-0
verbose<-TRUE

observe({ # - A
	input$sel_meas_comp 
	if(verbose){cat("\n in Comp_A")}
	if(isolate(init$a)=="TRUE"){
		do_isot<-(logfile$workflow[names(logfile$workflow)=="isotopologues"]=="yes")
		do_addu<-(logfile$workflow[names(logfile$workflow)=="adducts"]=="yes")
		do_homol<-(logfile$workflow[names(logfile$workflow)=="homologues"]=="yes")	
		do_EIC<-(logfile$workflow[names(logfile$workflow)=="EIC_correlation"]=="yes")
		if(
			file.exists(file.path(logfile[[1]],"results","componentization","components",isolate(input$sel_meas_comp))) &			
			isolate(input$sel_meas_comp)!="" & # ... finds emtpy folder otherwise
			( do_isot | do_addu | do_homol )
		){
			######################################################################		
			if(do_isot){ # update selection entries for the atom bound estimations
				updateSelectInput(session,inputId="atom_bound_addpeaks",
					choices=c("(a) peaks in same isotopologue group","(b) all peaks with similar RT"),
					selected="(a) peaks in same isotopologue group"
				)			
			}else{
				updateSelectInput(session,inputId="atom_bound_addpeaks",
					choices=c("(b) all peaks with similar RT"),
					selected="(b) all peaks with similar RT"
				)					
			}
			######################################################################
			# load componentization results ######################################
			load(file.path(logfile[[1]],"results","componentization","components",isolate(input$sel_meas_comp)),envir=as.environment(".GlobalEnv"))
			if(verbose){cat("\n in Comp_A_1")}
			cat("\n Loaded file")
			# output components summary table ####################################
			if(!is.null(dim(component[["pattern peak list"]]))){
				num_peaks_all<-dim(component[["pattern peak list"]])[1]
			}else{
				num_peaks_all<-dim(component[["adduct peak list"]])[1]
			}
			num_comp<-dim(component[[1]])[1]
			reduc<-round((num_peaks_all/num_comp),digits=2)
			num_isot_peaks<-rep(0,num_comp)
			num_adduc_peaks<-rep(0,num_comp)
			for(i in 1:num_comp){
				if(component[[1]][i,3]!="-"){
					num_isot_peaks[i]<-length(strsplit(component[[1]][i,3],",")[[1]])
				}
				if(component[[1]][i,5]!="-"){
					num_adduc_peaks[i]<-length(strsplit(component[[1]][i,5],",")[[1]])
				}
			}
			min2_size_comp<-round((sum((num_isot_peaks+num_adduc_peaks)>1)/num_comp),digits=2)
			median_size_comp<-round(mean(num_isot_peaks+num_adduc_peaks),digits=2)
			max_size_comp<-max(num_isot_peaks+num_adduc_peaks)
			output$num_peaks_all<-renderText(paste("Number of peaks: ",as.character(num_peaks_all),sep=""))
			output$num_comp<-renderText(paste("Number of components: ",as.character(num_comp),sep=""))
			output$reduc<-renderText(paste("Reduction factor: ",as.character(reduc),sep=""))
			output$min2_size_comp<-renderText(paste("Fraction of components with min. 2 peaks: ",as.character(min2_size_comp),sep=""))
			output$median_size_comp<-renderText(paste("Mean number of peaks per component: ",as.character(median_size_comp),sep=""))
			output$max_size_comp<-renderText(paste("Max number of peaks in a component: ",as.character(max_size_comp),sep=""))
			# output component table #############################################
			comp_table<-component[[1]][,c(1,15,13,14,16,3,5,6,7,11,12),drop=FALSE]
			comp_table[,2]<-round(comp_table[,2],digits=1)
			comp_table[,4]<-round(comp_table[,4],digits=5)			
			output$comp_table <- DT::renderDataTable(
				DT::datatable(
					comp_table,
					colnames=c(
						"Component ID",
						"Max. peak intens.","Max. peak ID","Max. peak m/z","Max. peak RT",
						"ID(s) isot. peaks","ID(s) adduct peaks","ID(s) homol. series","ID(s) interfering peaks",
						"Isot. peaks adducts","Adduct peak adducts"
					),
					rownames=FALSE,
					extensions = c('Buttons'),
					options = list(
						lengthMenu = c(100,200,400),
						ordering=F,
						dom = 'Bfrtip',
						buttons = c('excel')#buttons = c('excel', 'pdf', 'print', 'csv'),
					)
				),
				server = FALSE
			)
			# output text summary ################################################
			if((length(component[["pattern peak list"]])>1) & (do_isot)){found_isos<-TRUE}else{found_isos<-FALSE}
			if((length(component[["adduct peak list"]])>1) & (do_addu)){found_addu<-TRUE}else{found_addu<-FALSE}
			if((length(component[["homologue list"]])>1) & (do_homol)){found_homol<-TRUE}else{found_homol<-FALSE}
			get_comp_state<-"Available components for this file:"
			if(found_isos & found_addu & found_homol){
				get_comp_state<-paste(get_comp_state," isotopologues, adducts and homologues",sep="")
			}
			if(found_isos & found_addu & !found_homol){
				get_comp_state<-paste(get_comp_state," isotopologues and adducts",sep="")
			}
			if(found_isos & !found_addu & found_homol){
				get_comp_state<-paste(get_comp_state," isotopologues and homologues",sep="")
			}
			if(!found_isos & found_addu & found_homol){
				get_comp_state<-paste(get_comp_state," adducts and homologues",sep="")
			}
			if(found_isos & !found_addu & !found_homol){
				get_comp_state<-paste(get_comp_state," only isotopologues",sep="")
			}
			if(!found_isos & found_addu & !found_homol){
				get_comp_state<-paste(get_comp_state," only adducts",sep="")
			}
			if(!found_isos & !found_addu & found_homol){
				get_comp_state<-paste(get_comp_state," only homologues homologues",sep="")
			}
			output$sel_meas_comp_state<-renderText(get_comp_state)
		}else{
			updateSelectInput(session,inputId="atom_bound_addpeaks",
				choices=c("(b) all peaks with similar RT"),
				selected="(b) all peaks with similar RT"
			)								
			if(!file.exists(file.path(logfile[[1]],"peaklist",isolate(input$sel_meas_comp)))){
				output$sel_meas_comp_state<-renderText("Invalid file ID")			
			}else{
				output$sel_meas_comp_state<-renderText("No componentization results for this file available")
			}
			if(verbose){cat("\n in Comp_A_2")}
		}
		# load homologue results #############################################
		if(
			file.exists(file.path(logfile[[1]],"results","componentization","homologues",paste("full",isolate(input$sel_meas_comp),sep="_"))) &
			do_homol 
		){
			if(verbose){cat("\n in Comp_A_3")}
			load(file.path(logfile[[1]],"results","componentization","homologues",paste("full",isolate(input$sel_meas_comp),sep="_")),envir=as.environment(".GlobalEnv"))			
			# output homol. series table #####################################
			output$homol_table <- DT::renderDataTable(
				datatable(
					cbind(homol[[3]][,c(1,2,3)],round(homol[[3]][,4],digits=1)),
					colnames=c("Series ID","Peak IDs","m/z difference","RT difference [s]"),
					rownames=FALSE
				)
			)
			# output homol. series plot ######################################
			output$homol_plot <- renderPlot({					
				nontarget:::plothomol(homol,
					xlim = FALSE, ylim = FALSE
				)
			},res=100)	
		}
		######################################################################
		#if(
		#	file.exists(file.path(logfile[[1]],"results","componentization","EIC_corr",as.character(isolate(input$sel_meas_comp)))) &
		#	do_EIC 
		#){
		#if(FALSE){
		#	if(verbose){cat("\n in Comp_A_4")}
			#load(file.path(logfile[[1]],"results","componentization","EIC_corr",as.character(isolate(input$sel_meas_comp)))) 
			#load(file=file.path(logfile[[1]],"peaklist",as.character(isolate(input$sel_meas_comp))));   
			#EIC_pairs			
		#}
		#}
	}
	##########################################################################
	
})	

observe({ # - B
	input$sel_meas_comp_peak 
	if(verbose){cat("\n in Comp_B")}
	if(!is.na(as.numeric(isolate(input$sel_meas_comp_peak)))){
	if(isolate(init$a)=="TRUE" & as.numeric(isolate(input$sel_meas_comp_peak))>0){
		if(
			file.exists(file.path(logfile[[1]],"results","componentization","components",isolate(input$sel_meas_comp))) &
			(isolate(input$sel_meas_comp_peak)>0)
		){
			# search in isotop. peaks
			that<-which(!is.na(unlist(lapply(strsplit(component[[1]][,3],","), match, x=as.numeric(isolate(input$sel_meas_comp_peak))))))
			# search in adduct peaks
			if(length(that)==0){
				that<-which(!is.na(unlist(lapply(strsplit(component[[1]][,5],","), match, x=as.numeric(isolate(input$sel_meas_comp_peak))))))
			}
			# search in interfering peaks
			if(length(that)==0){
				that<-which(!is.na(unlist(lapply(strsplit(component[[1]][,7],","), match, x=as.numeric(isolate(input$sel_meas_comp_peak))))))
			}
			if(length(that)==1){
				if(verbose){cat("\n in Comp_B_1")}
				updateNumericInput(session,"sel_meas_comp_comp",value=that)
			}else{
				if(verbose){cat("\n in Comp_B_2")}
				cat("\n Invalid peak selected!")
				updateNumericInput(session,"sel_meas_comp_comp",value=0)
			}
		}
	}
	}
})	

observe({ # - C
	input$sel_meas_comp_comp 
	if(verbose){cat("\n in Comp_C")}
	if(!is.na(as.numeric(isolate(input$sel_meas_comp_comp)))){
	if(isolate(init$a)=="TRUE" & as.numeric(isolate(input$sel_meas_comp_comp))){
		if(
			file.exists(file.path(logfile[[1]],"results","componentization","components",isolate(input$sel_meas_comp))) &
			(isolate(input$sel_meas_comp_comp)>0)
		){
			if(
				any(component[[1]][,1]==as.numeric(isolate(input$sel_meas_comp_comp)))
			){
				if(verbose){cat("\n in Comp_C_1")}
				ee$entry<-as.numeric(isolate(input$sel_meas_comp_comp))	
			}else{
				if(verbose){cat("\n in Comp_C_2")}
				cat("\n Invalid component selected!")
				updateNumericInput(session,"sel_meas_comp_comp",value=0)
				updateNumericInput(session,"sel_meas_comp_peak",value=0)
			}			
		}
	}
	}
})	

observe({ # - D: generate outputs
	ee$entry 
	if(isolate(ee$entry)>0){
		if(verbose){cat("\n in Comp_D_1")}
		got_comp<-enviMass:::plotcomp_parts(component, compoID=as.numeric(isolate(ee$entry)), what="check")
		if(got_comp){
			output$found_compo<-renderText("")
			# output spectrum
			output$comp_plot_spec <- renderPlot({	
				enviMass:::plotcomp_parts(component, compoID=as.numeric(isolate(ee$entry)), what="spec")
			},res=110)				
			# output circular plot
			output$comp_plot_circ <- renderPlot({	
				enviMass:::plotcomp_parts(component, compoID=as.numeric(isolate(ee$entry)), what="circ")
			},res=110)				
			# output tables 
			comp_table<<-enviMass:::plotcomp_parts(component, compoID=as.numeric(isolate(ee$entry)), what="table")
			inser<-rep("",length(comp_table$relations[,2]))
			do_these<-which(grepl("<->",comp_table$relations[,2]))
			if(length(do_these)>0){
				inser[do_these]<-paste("Same isotopologues of different adducts:",comp_table$relations[do_these,2],sep=" ")
			}
			do_these<-which(grepl(",z=",comp_table$relations[,2]))
			if(length(do_these)>0){					
				inser[do_these]<-"Different isotopologues of the same adduct"
			}	
			output$comp_table_a <- DT::renderDataTable(
					datatable(
					cbind(
						unlist(lapply(strsplit(comp_table$relations[,1],"-"), `[[`,1)),
						unlist(lapply(strsplit(comp_table$relations[,1],"-"), `[[`,2)),
						inser,
						round(comp_table$relations[,3],digits=3)					
					),
					colnames=c("First peaks","Second peak","Links","Intensity ratio")
				)
			)
			output$comp_table_b <- DT::renderDataTable(
				datatable(
					cbind(
						round(comp_table$a[,1],digits=1),
						round(comp_table$a[,2],digits=5),
						format(comp_table$a[,3],scientific=TRUE,digits=2),
						round(comp_table$a[,4],digits=2)						
					),
					colnames=c("Peak ID","m/z","Intensity","RT [s]")
				)
			)
			if(comp_table[[4]]=="Not part of a homologue series"){
				output$comp_table_c  <- DT::renderDataTable(
					datatable(
						data.frame("Peaks in selected component are not part of any homologue series"),
						colnames="No results"
					)
				)
			}else{
				these<-unique(strsplit(comp_table[[4]],"/")[[1]][-1])
				output$comp_table_c  <- DT::renderDataTable(
					datatable(
						cbind(
							homol[[3]][these,1],
							round(homol[[3]][these,3],digits=4),
							round(homol[[3]][these,4],digits=1)
						),
						colnames=c("Series ID","m/z difference","RT difference [s]")
					)
				)
			}
			output$comp_table_d <- DT::renderDataTable(
				datatable(
					cbind(
						round(comp_table[[2]][,1],digits=1),
						round(comp_table[[2]][,2],digits=5),
						format(comp_table[[2]][,3],scientific=TRUE,digits=2),
						round(comp_table[[2]][,4],digits=2)						
					),							
					colnames=c("Peak ID","m/z","Intensity","RT [s]")	
				)
			)
		}else{
			output$found_compo<-renderText("The selected component contains only one peak.")
			output$comp_plot_spec <- renderPlot({	
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(.5,.5,labels="The selected component contains only one peak.")
			},res=110)				
			# output circular plot
			output$comp_plot_circ <- renderPlot({	
				plot.new()
			},res=110)				
			# output tables 
			output$comp_table_a <- DT::renderDataTable(
				datatable(data.frame("No results"),colnames="No results")
			)
			output$comp_table_b <- DT::renderDataTable(
				datatable(data.frame("No results"),colnames="No results")
			)
			output$comp_table_c  <- DT::renderDataTable(
				datatable(data.frame("No results"),colnames="No results")
			)
			output$comp_table_d <- DT::renderDataTable(
				datatable(data.frame("No results"),colnames="No results")
			)		
		}	
		######################################################################
		
	}
})

##############################################################################
  
  
##############################################################################
# Calculate atom bounds ######################################################
# Update shift bounds from element election
output$atom_bounds_that <- renderUI({
	if(length(input$atom_bounds_this)){
		cat("\n Am observing")
		those_elements <- input$atom_bounds_this
		lapply(those_elements, function(i){
			numericInput(inputId = paste0("ppm_", i),label = paste0("Maximum shifts for ", i),value = 30, width='200px')
		})
	}else{
		cat("\n Am observing nothing")		
	}
})

observe({ # - F: generate outputs
	input$atom_bound_peak 
	input$sel_meas_comp
	input$atom_bound_addpeaks
	input$atom_bounds_calculate
	if(isolate(init$a)=="TRUE"){
		
		if(file.exists(file.path(logfile[[1]],"peaklist",isolate(input$sel_meas_comp)))){
			found_peaklist<-TRUE	
		}else{
			found_peaklist<-FALSE
		}
		if(!is.na(isolate(input$atom_bound_peak))>0){
			if(	
				(isolate(input$atom_bound_peak)>0) & 
				(isolate(input$sel_meas_comp)>0) &
				found_peaklist
			){
				if(verbose){cat("\n in Atoms_1")}
				######################################################################
				# get additional peaks ###############################################
				if(isolate(input$atom_bound_addpeaks)=="(b) all peaks with similar RT"){ # load peaklist 
					if(verbose){cat("\n in Atoms_2")}		
					load(file.path(logfile[[1]],"peaklist",isolate(input$sel_meas_comp)),envir=as.environment(".GlobalEnv"))
					at_peak<<-which(peaklist[,"peak_ID"]==isolate(input$atom_bound_peak))
					if(length(at_peak)>0){
						atom_peaks<<-peaklist[
							(peaklist[,"m/z_corr"]>peaklist[at_peak,"m/z_corr"]) &
							((abs(peaklist[,"RT_corr"]-peaklist[at_peak,"RT_corr"]))<=logfile$parameters$isotop_rttol) &
							(abs(peaklist[,"m/z_corr"]-peaklist[at_peak,"m/z_corr"])<=10)
						,,drop=FALSE]
						atom_peaks<<-rbind(
							peaklist[at_peak,],
							atom_peaks)
					}	
				}else{ # get peaks from components - only selectable if selectInput adapted accordingly
					if(verbose){cat("\n in Atoms_3")}	
					at_peak<<-which(!is.na(unlist(lapply(strsplit(component[[1]][,3],","), match, x=as.numeric(isolate(input$atom_bound_peak))))))
					# search in adduct peaks
					if(length(at_peak)==0){
						at_peak<<-which(!is.na(unlist(lapply(strsplit(component[[1]][,5],","), match, x=as.numeric(isolate(input$atom_bound_peak))))))
					}
					# search in interfering peaks
					if(length(at_peak)==0){
						at_peak<<-which(!is.na(unlist(lapply(strsplit(component[[1]][,7],","), match, x=as.numeric(isolate(input$atom_bound_peak))))))
					}				
					if(length(at_peak)>0){	# at_peak = at which component		
						get_peaks<<-c()
						if(component[[1]][at_peak,3]!="-"){
							get_peaks<<-c(get_peaks,strsplit(component[[1]][at_peak,3],",")[[1]])
						}
						if(component[[1]][at_peak,5]!="-"){
							get_peaks<<-c(get_peaks,strsplit(component[[1]][at_peak,5],",")[[1]])
						}			
						if(component[[1]][at_peak,7]!="-"){
							get_peaks<<-c(get_peaks,strsplit(component[[1]][at_peak,7],",")[[1]])
						}			
						get_peaks<<-as.numeric(get_peaks)
						peaklist<<-component[[2]][get_peaks,1:4,drop=FALSE]
						at_peak<<-which(peaklist[,"peak ID"]==isolate(input$atom_bound_peak))
						atom_peaks<<-peaklist[
							(peaklist[,"m/z_corr"]>peaklist[at_peak,"m/z_corr"]) &
							((abs(peaklist[,"RT_corr"]-peaklist[at_peak,"RT_corr"]))<=logfile$parameters$isotop_rttol) &
							(abs(peaklist[,"m/z_corr"]-peaklist[at_peak,"m/z_corr"])<=10)
						,1:3,drop=FALSE]
						atom_peaks<<-rbind(
							peaklist[at_peak,1:3],
							atom_peaks)
					}else{
						cat("\n Invalid peak selected - not found among components.")
					}
				}
				######################################################################
				# get LOD ############################################################
				if(length(at_peak)>0){
					if(
						(logfile$workflow[names(logfile$workflow)=="LOD"]=="yes") &
						(file.exists(file=file.path(logfile$project_folder,"results","LOD","LOD_splined")))
					){
						if(verbose){cat("\n in Atoms_4")}
						load(file=file.path(logfile$project_folder,"results","LOD","LOD_splined"));
						with_model<-which(names(LOD_splined)==paste("LOD_",as.character(isolate(input$sel_meas_comp)),sep=""))			
						#with_model<-which(names(LOD_splined)==paste("LOD_","1163",sep=""))			
						if(length(with_model)>0){						
							use_LOD<<-10^(predict(LOD_splined[[with_model]],atom_peaks[1,"RT_corr"])$y[[1]])	
						}else{
							cat("\n Shouldn`t there be a LOD spline (case A)? Could not find it!")
							use_LOD<<-logfile$parameters$tar_intcut			
						}					
					}else{
						if(verbose){cat("\n in Atoms_5")}
						cat("\n Shouldn`t there be a LOD spline (case B)? Could not find it!")
						use_LOD<<-logfile$parameters$tar_intcut
					}
					##################################################################
					# plot peaks #####################################################
					output$atom_bound_plot_peak <- renderPlot({	
						plot(
							atom_peaks[,"m/z_corr"],atom_peaks[,"int_corr"],
							xlab="m/z",ylab="Intensity",type="h",col="darkgreen",lwd=2,
							ylim=c(0,max(atom_peaks[,"int_corr"])))
						points(
							atom_peaks[1,"m/z_corr"],atom_peaks[1,"int_corr"],
							type="h",col="red",lwd=2)
						abline(h=use_LOD,col="gray")
					},res=110)				
					##################################################################		
					# run estimate function ##########################################
					elements<<-isolate(input$atom_bounds_this)
					charges<<-c(1,2,3,4)
					if(length(elements)>0){
						dmz<<-c()	
						lapply(elements,function(i){
							x<-paste0("ppm_",i)
							dmz<<-c(dmz,
									input[[x]]
							)
						})
						atom_counts<<-try({					
							enviMass:::atoms(
								masses=atom_peaks[,"m/z_corr"],
								intensities=atom_peaks[,"int_corr"],
								elements,
								dmz,
								ppm=TRUE,
								charges,
								isotopes,
								int_cut=use_LOD,
								inttol=0.2,
								use_C=as.logical(isolate(input$atom_bound_wcarbon)),
								must_peak=FALSE
							)
						})
					}		
					if(class(atom_counts)!="try-error"){
						if(verbose){cat("\n in Atoms_6")}
						output$atom_count_table <- DT::renderDataTable(			
							datatable(as.data.frame(cbind(charges,atom_counts),row.names = NULL),
								colnames=c("Charge z",elements),rownames=FALSE)
						)
					}else{
						output$atom_count_table <- DT::renderDataTable(
							datatable(data.frame("No results- sth went wrong - debug?"),colnames="No results")
						)
						cat("\n Atom count estimation failed - debug!")
					}
					##################################################################			
				}else{
					output$atom_bound_plot_peak <- renderPlot({	
						plot.new()
						plot.window(xlim=c(0,1),ylim=c(0,1))
						text(0.5,0.5,labels="Invalid peak ID")
					},res=110)			
				}	
				######################################################################		
			}else{
				output$atom_bound_plot_peak <- renderPlot({	
					plot.new()
					plot.window(xlim=c(0,1),ylim=c(0,1))
					text(0.5,0.5,labels="No peaks available")
				},res=110)			
			}
		}
	}
})

##############################################################################
  
  
  
  
  
  
  
  
  
  
  
  
