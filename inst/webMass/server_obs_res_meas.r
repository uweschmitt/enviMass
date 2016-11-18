##############################################################################
# update results for individual measurements #################################
##############################################################################
observe({
    input$sel_meas
	if(!is.na(isolate(input$sel_meas))){
    if(isolate(input$sel_meas)!=0){
		##########################################################################	
		measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		if(any(measurements[,"ID"]==as.character(isolate(input$sel_meas)))){
			output$file_proc_name<-renderText(measurements[measurements[,"ID"]==as.character(isolate(input$sel_meas)),"Name"])
			output$file_proc_type<-renderText(measurements[measurements[,"ID"]==as.character(isolate(input$sel_meas)),"Type"])
			output$file_proc_mode<-renderText(measurements[measurements[,"ID"]==as.character(isolate(input$sel_meas)),"Mode"])		
			# peaklist info ##########################################################
			if(file.exists(file.path(logfile$project_folder,"peaklist",as.character(isolate(input$sel_meas))))){
				load(file=file.path(logfile$project_folder,"peaklist",as.character(isolate(input$sel_meas))),envir=as.environment(".GlobalEnv"),verbose=FALSE);			
				output$file_peak_number<-renderText(as.character(length(peaklist[,1])));	
				blind_rem<-round(
					(sum(peaklist[,colnames(peaklist)=="keep_2"]==0))/length(peaklist[,1])*100
				,digits=3)
				output$file_blind_rem<-renderText(as.character(blind_rem));
				repl_rem<-round(
					(sum(peaklist[,colnames(peaklist)=="keep"]==0))/length(peaklist[,1])*100
				,digits=3)
				output$file_repl_rem<-renderText(as.character(repl_rem));
				#rm(peaklist,envir=as.environment(".GlobalEnv")) # wtf?
			}else{
				cat("\n no peaklist for processing view found")
			}
			##########################################################################		
			pics<-list.files(file.path(logfile[[1]],"pics"))
			# recalibration ##########################################################
			if(
				any(pics==paste("recal_",as.character(isolate(input$sel_meas)),sep="")) & (logfile$workflow[names(logfile$workflow)=="recal"]=="yes")
			){
				expr1<-list(src=file.path(logfile[[1]],"pics",paste("recal_",as.character(isolate(input$sel_meas)),sep="")))
				output$recal_pic<-renderImage(expr1, deleteFile = FALSE)
				cat("\n Found recal_pic")
			}else{
				cat("\n Not found recal_pic")			
			}
			##########################################################################
			# LOD  ###################################################################
			if( file.exists( file.path(logfile[[1]],"results","LOD",paste("plot_LOD_",as.character(isolate(input$sel_meas)),".png",sep="") ) ) ){
				expr_LOD<-list( src=file.path(logfile[[1]],"results","LOD",paste("plot_LOD_",as.character(isolate(input$sel_meas)),".png",sep="")) )
				output$LOD_pic<-renderImage(expr_LOD, deleteFile = FALSE)	
			}else{
				cat("\n LOD pic file not found")
			}
			##########################################################################
			output$dowhat<-renderText("Processing per file viewed.");	
		}else{
			output$dowhat<-renderText("Invalid ID chosen to view processing results.");		
		}
    }
	}
})
##############################################################################

##############################################################################
# retrieve peak information ##################################################
##############################################################################
ranges_peaks_mz_RT <- reactiveValues(x = NULL, y = NULL)
refresh_plot<-reactiveValues()
refresh_plot$a<-1
observe({
	refresh_plot$a
    input$sel_meas_ID
	input$peaks_mz_RT_use_peaks
	input$peaks_mz_RT_use_raw
	input$peaks_mz_RT_use_IDs
	input$peaks_mz_RT_use_window
	input$peaks_mz_RT_use_window_mass
	input$peaks_mz_RT_use_window_RT
	input$peaks_mz_RT_use_bar	
	input$peaks_mz_RT_use_bar_value	
	input$peaks_mz_RT_use_window_RT_tol
	input$plot_filter_intensity
	input$plot_filter_blind
	input$plot_filter_replicates
	if(isolate(init$a)=="TRUE"){
		if(!is.na(isolate(input$sel_meas_ID))){ 
			if(!any(objects(envir=as.environment(".GlobalEnv"))=="atit")){ # atit -> dont load the same MSlist twice = too slow
				assign("atit",0,envir=as.environment(".GlobalEnv"))
			}	
			if(	
				file.exists(file.path(logfile[[1]],"MSlist",as.character(isolate(input$sel_meas_ID)))) &
				file.exists(file.path(logfile[[1]],"peaklist",as.character(isolate(input$sel_meas_ID)))) 			
			){
				if(verbose){cat("\n Plotting")}
				if((isolate(input$sel_meas_ID)!=atit)){
					if(any(objects(envir=as.environment(".GlobalEnv"))=="MSlist")){rm(MSlist,envir=as.environment(".GlobalEnv"))}
					if(any(objects()=="MSlist")){rm(MSlist)}				
					load(file.path(logfile[[1]],"MSlist",as.character(isolate(input$sel_meas_ID))), envir=as.environment(".GlobalEnv"))
					load(file.path(logfile[[1]],"peaklist",as.character(isolate(input$sel_meas_ID))), envir=as.environment(".GlobalEnv"))
					if(verbose){cat("\nafter:\n "); print(gc());cat("\n MSlist file loaded")}
					assign("atit",isolate(input$sel_meas_ID),envir=as.environment(".GlobalEnv"))
					max_int<-max(max(MSlist[["Peaklist"]][,"max_int"]),max(MSlist[["Scans"]][[2]][,"intensity"]))
					min_int<-min(min(MSlist[["Peaklist"]][,"max_int"]),min(MSlist[["Scans"]][[2]][,"intensity"]))
					updateSliderInput(session,"plot_filter_intensity", min=round(log10(min_int),digits=1), max=round(log10(max_int),digits=1), value=c(log10(min_int),log10(max_int)))
				}
				# prepare plotting information #############################################################
				if(!is.null(isolate(ranges_peaks_mz_RT$x))){
					x_lim<-isolate(ranges_peaks_mz_RT$x)
				}else{
					x_lim<-c(min(MSlist[["Scans"]][[2]][,"m/z"]),max(MSlist[["Scans"]][[2]][,"m/z"]))
				}
				if(!is.null(isolate(ranges_peaks_mz_RT$y))){
					y_lim<-isolate(ranges_peaks_mz_RT$y)
				}else{
					y_lim<-c(min(MSlist[["Scans"]][[2]][,"RT"]),max(MSlist[["Scans"]][[2]][,"RT"]))
				}		
				use_these<-which(
					(peaklist[,"m/z"]>=x_lim[1]) &
					(peaklist[,"m/z"]<=x_lim[2]) &
					(peaklist[,"RT"]>=y_lim[1]) &
					(peaklist[,"RT"]<=y_lim[2]) &
					(peaklist[,"max_int"]>=isolate(10^input$plot_filter_intensity[1])) &
					(peaklist[,"max_int"]<=isolate(10^input$plot_filter_intensity[2])) 
				)
				if((length(use_these)>0) & isolate(input$plot_filter_blind)){
					use_these<-use_these[peaklist[use_these,"keep_2"]==1]
				}
				if((length(use_these)>0) & isolate(input$plot_filter_replicates)){
					use_these<-use_these[peaklist[use_these,"keep"]==1]
				}				
				if(length(use_these)>0){
					z_lim<-c(min(peaklist[use_these,"max_int"]),max(peaklist[use_these,"max_int"]))
				}else{
					z_lim<-c(min(peaklist[,"max_int"]),max(peaklist[,"max_int"]))
				}
				if(isolate(input$peaks_mz_RT_use_raw)){ # must be prerequisite to evaluate those, below
					if(is.null(isolate(ranges_peaks_mz_RT$x))){
						x_min<-min(MSlist[["Scans"]][[2]][,"m/z"])
						x_max<-max(MSlist[["Scans"]][[2]][,"m/z"])							
					}else{
						x_min<-min(MSlist[["Scans"]][[2]][
							MSlist[["Scans"]][[2]][,"m/z"]>=isolate(ranges_peaks_mz_RT$x[1])
						,"m/z"])
						x_max<-max(MSlist[["Scans"]][[2]][
							MSlist[["Scans"]][[2]][,"m/z"]<=isolate(ranges_peaks_mz_RT$x[2])
						,"m/z"])	
					}
					if(is.null(isolate(ranges_peaks_mz_RT$y))){
						y_min<-min(MSlist[["Scans"]][[2]][,"RT"])
						y_max<-max(MSlist[["Scans"]][[2]][,"RT"])							
					}else{
						y_min<-min(MSlist[["Scans"]][[2]][
							MSlist[["Scans"]][[2]][,"RT"]>=isolate(ranges_peaks_mz_RT$y[1])
						,"RT"])
						y_max<-max(MSlist[["Scans"]][[2]][
							MSlist[["Scans"]][[2]][,"RT"]<=isolate(ranges_peaks_mz_RT$y[2])
						,"RT"])	
					}					
					those<-which(
						(MSlist[["Scans"]][[2]][,"m/z"]>=x_min) &
						(MSlist[["Scans"]][[2]][,"m/z"]<=x_max) &
						(MSlist[["Scans"]][[2]][,"RT"]>=y_min) &
						(MSlist[["Scans"]][[2]][,"RT"]<=y_max) &
						(MSlist[["Scans"]][[2]][,"intensity"]>=isolate(10^input$plot_filter_intensity[1])) &
						(MSlist[["Scans"]][[2]][,"intensity"]<=isolate(10^input$plot_filter_intensity[2]))						
					)			
					if(length(those)<=1E6){
						colorit<-rep("gray",length(those))
						colorit[MSlist[["Scans"]][[2]][those,"peakID"]!=0]<-"red"
					}
				}else{
					those<-c()
				}
				####################################################################################			
				output$plot_peaks_mz_RT <- renderPlot({
					par(mar=c(4, 4, 3, .1))
					plot.new()
					plot.window(xlim=x_lim,ylim=y_lim)
					title(
						xlab="m/z [Th]", ylab="RT",
						main="Draw rectangles and double-click into them to zoom in, double-click again to zoom fully out. Bottom plots adapt accordingly.",cex.main=.75
					)
					box();axis(1);axis(2);
					# add raw data ? #########################################################
					if(isolate(input$peaks_mz_RT_use_raw)){
						if(length(those)<=1E5 & length(those)>0){
							points(
								MSlist[["Scans"]][[2]][those,"m/z"],
								MSlist[["Scans"]][[2]][those,"RT"],
								pch=19,cex=.5,col=colorit
							)
						}else{
							smoothScatter(
								x=MSlist[["Scans"]][[2]][those,"m/z"], 
								y=MSlist[["Scans"]][[2]][those,"RT"],
								 colramp = colorRampPalette(c("white", "red")),
								nbin = 200, add = TRUE
							)
						}
					}
					# add picked peaks ? ######################################################
					if(isolate(input$peaks_mz_RT_use_peaks)){		
						if(length(use_these)>0){
							if(!isolate(input$peaks_mz_RT_use_IDs)){
								colorit<-"black"
							}else{
								colorit<-"darkgrey"
							}
							points(
								peaklist[use_these,"m/z"],
								peaklist[use_these,"RT"],						
								pch=21,cex=.8,col=colorit
							)
						}
					}
					# add peak IDs? ###########################################################
					if(isolate(input$peaks_mz_RT_use_IDs)){
						if(length(use_these)>0){
							text(						
								peaklist[use_these,"m/z"],
								peaklist[use_these,"RT"], 
								labels = as.character(peaklist[use_these,"peak_ID"]),
								pos = NULL, col	= "darkred", cex=.6
							)
						}					
					}					
					# add search window ? ######################################################
					if(isolate(input$peaks_mz_RT_use_window)){	
						del_mz<-((x_lim[2]-x_lim[1])/15)
						at_mz<-isolate(input$peaks_mz_RT_use_window_mass)
						at_RT<-isolate(input$peaks_mz_RT_use_window_RT)
						rect(
							xleft=(at_mz-del_mz), 
							ybottom=(at_RT-isolate(input$peaks_mz_RT_use_window_RT_tol)), 
							xright=(at_mz+del_mz), 
							ytop=(at_RT+isolate(input$peaks_mz_RT_use_window_RT_tol)),
							col=NULL, border="blue",lwd=2)
						if(isolate(input$peaks_mz_RT_use_bar)){		
							del_ppm<-(at_mz*isolate(input$peaks_mz_RT_use_bar_value)/1E6)
							lines(
								x=c((at_mz-del_ppm),(at_mz+del_ppm)),
								y=c(at_RT,at_RT),
								col="blue",lwd=2)
						}	
					}
					############################################################################	
				}, res = 100, execOnResize=TRUE)
				####################################################################################
				output$plot_peaks_mz_int <- renderPlot({		
						par(mar=c(4, 4, .8, .2))
						plot.new()
						plot.window(xlim=x_lim,ylim=c(0,z_lim[2]))
						title(xlab="m/z [Th]", ylab="Intensity")
						box();axis(1);axis(2);
						# add raw data ? #########################################################			
						if(isolate(input$peaks_mz_RT_use_raw)){			
							if(length(those)<=1E5 & length(those)>0){
								points(
									MSlist[["Scans"]][[2]][those,"m/z"],
									MSlist[["Scans"]][[2]][those,"intensity"],
									type="h",pch=19,cex=.5,col=colorit
								)
							}						
						}
						# add picked peaks ? #####################################################
						if(isolate(input$peaks_mz_RT_use_peaks)){						
							if(length(use_these)>0){
								if(!isolate(input$peaks_mz_RT_use_IDs)){
									colorit<-"black"
								}else{
									colorit<-"darkgrey"
								}
								points(
									peaklist[use_these,"m/z"],
									peaklist[use_these,"max_int"],						
									pch=21,cex=.8,col=colorit,type="h"
								)					
							}
						}
						# add ppm bar? #############################################################
						if(isolate(input$peaks_mz_RT_use_window)){	
							if(isolate(input$peaks_mz_RT_use_bar)){	
								at_mz<-isolate(input$peaks_mz_RT_use_window_mass)							
								del_ppm<-(at_mz*isolate(input$peaks_mz_RT_use_bar_value)/1E6)
								lines(
									x=c((at_mz-del_ppm),(at_mz+del_ppm)),
									y=c(0.5*z_lim[2],.5*z_lim[2]),
									col="blue",lwd=2)
							}	
						}						
				}, res = 100, execOnResize=TRUE)	
				####################################################################################
				output$plot_peaks_RT_int <- renderPlot({		
						par(mar=c(4, 4, .8, .2))
						plot.new()
						plot.window(xlim=y_lim,ylim=c(0,z_lim[2]))
						title(xlab="RT", ylab="Intensity")
						box();axis(1);axis(2);
						# add raw data ? #########################################################			
						if(isolate(input$peaks_mz_RT_use_raw)){			
							if(length(those)<=1E5 & length(those)>0){
								points(
									MSlist[["Scans"]][[2]][those,"RT"],
									MSlist[["Scans"]][[2]][those,"intensity"],
									type="h",pch=19,cex=.5,col=colorit
								)
							}						
						}
						# add picked peaks ? #####################################################
						if(isolate(input$peaks_mz_RT_use_peaks)){						
							if(length(use_these)>0){
								if(!isolate(input$peaks_mz_RT_use_IDs)){
									colorit<-"black"
								}else{
									colorit<-"darkgrey"
								}
								points(
									peaklist[use_these,"RT"],
									peaklist[use_these,"max_int"],						
									pch=21,cex=.8,col=colorit,type="h"
								)					
							}
						}
						# add RT window? ##########################################################
						if(isolate(input$peaks_mz_RT_use_window)){	
								at_RT<-isolate(input$peaks_mz_RT_use_window_RT)
								lines(
									x=c(
										(at_RT-isolate(input$peaks_mz_RT_use_window_RT_tol)),
										(at_RT+isolate(input$peaks_mz_RT_use_window_RT_tol))
									),
									y=c(.5*z_lim[2],.5*z_lim[2]),
									col="blue",lwd=2)
							}	
				}, res = 100, execOnResize=TRUE)	
				####################################################################################		
				# plotly output ####################################################################
				# only peaks? ######################################################################
				if(isolate(input$peaks_mz_RT_use_peaks) & !isolate(input$peaks_mz_RT_use_raw)){
					if(length(use_these)<=1E5 & length(use_these)>0){
						sub_peaks<-as.data.frame(peaklist[use_these,c("m/z","RT","max_int")])
						names(sub_peaks)<-c("m_z","RT","Intensity")				
						sub_peaks[,"Intensity"]<-(sub_peaks[,"Intensity"]/2)
						output$plot_peaks_3D <- renderPlotly({
							p <- plotly:::plot_ly(type="scatter3d",mode="markers",showlegend=TRUE)%>%
							plotly:::add_trace(p,
								x = ~m_z, y = ~RT, z = ~Intensity,
								data = sub_peaks,
								color=I("black"),
								size = I(1),
								name = "",
								error_z=list(
									color="black",
									thickness=0,
									symmetric = TRUE, 
									type = "data" ,
									array = sub_peaks$Intensity
								)
							)	
						})					
					}else{
						output$plot_peaks_3D <- renderPlotly({plotly:::plot_ly(type="scatter3d",mode="markers",showlegend=TRUE)})
					}
				}
				# only raw data? ###########################################################						
				if(!isolate(input$peaks_mz_RT_use_peaks) & isolate(input$peaks_mz_RT_use_raw)){
					if(length(those)<=1E5 & length(those)>0){  # implies number of peaks is lower, too
						sub_MSlist<-as.data.frame(MSlist[["Scans"]][[2]][those,c("m/z","RT","intensity","peakID")])
						names(sub_MSlist)<-c("m_z","RT","Intensity","peakID")		
						sub_MSlist[,"Intensity"]<-(sub_MSlist[,"Intensity"]/2)
						output$plot_peaks_3D <- renderPlotly({
							p <- plotly:::plot_ly(type="scatter3d",mode="markers",showlegend=TRUE)%>%
							plotly:::add_trace(p,
								x = ~m_z, y = ~RT, z = ~Intensity,
								data = sub_MSlist[sub_MSlist[,"peakID"]==0,],
								color=I("gray"),
								size = I(1),
								name = "",
								error_z=list(
									color="gray",
									thickness=0,
									symmetric = TRUE, 
									type = "data" ,
									array = sub_MSlist[sub_MSlist[,"peakID"]==0,]$Intensity
								)
							)%>%	
							plotly:::add_trace(p,
								x = ~m_z, y = ~RT, z = ~Intensity,
								data = sub_MSlist[sub_MSlist[,"peakID"]!=0,],
								color=I("red"),
								size = I(1),
								name = "",
								error_z=list(
									color="red",
									thickness=0,
									symmetric = TRUE, 
									type = "data" ,
									array = sub_MSlist[sub_MSlist[,"peakID"]!=0,]$Intensity
								)
							)
						})					
					}else{
						output$plot_peaks_3D <- renderPlotly({plotly:::plot_ly(type="scatter3d",mode="markers",showlegend=TRUE)})					
					}
				}					
				# peaks & raw data? ########################################################
				if( isolate(input$peaks_mz_RT_use_peaks) & isolate(input$peaks_mz_RT_use_raw)){
					if(length(those)<=1E5 & length(those)>0){ # implies number of peaks is lower, too
						sub_MSlist<-as.data.frame(MSlist[["Scans"]][[2]][those,c("m/z","RT","intensity","peakID")])
						names(sub_MSlist)<-c("m_z","RT","Intensity","peakID")		
						sub_MSlist[,"Intensity"]<-(sub_MSlist[,"Intensity"]/2)
						sub_peaks<-as.data.frame(peaklist[use_these,c("m/z","RT","max_int")])
						names(sub_peaks)<-c("m_z","RT","Intensity")				
						sub_peaks[,"Intensity"]<-(sub_peaks[,"Intensity"]/2)
						output$plot_peaks_3D <- renderPlotly({
							p <- plotly:::plot_ly(type="scatter3d",mode="markers",showlegend=TRUE)%>%
							plotly:::add_trace(p,
								x = ~m_z, y = ~RT, z = ~Intensity,
								data = sub_MSlist[sub_MSlist[,"peakID"]==0,],
								color=I("gray"),
								size = I(1),
								name = "",
								error_z=list(
									color="gray",
									thickness=0,
									symmetric = TRUE, 
									type = "data" ,
									array = sub_MSlist[sub_MSlist[,"peakID"]==0,]$Intensity
								)
							)%>%	
							plotly:::add_trace(p,
								x = ~m_z, y = ~RT, z = ~Intensity,
								data = sub_MSlist[sub_MSlist[,"peakID"]!=0,],
								color=I("red"),
								size = I(1),
								name = "",
								error_z=list(
									color="red",
									thickness=0,
									symmetric = TRUE, 
									type = "data" ,
									array = sub_MSlist[sub_MSlist[,"peakID"]!=0,]$Intensity
								)
							)%>%
							plotly:::add_trace(p,
								x = ~m_z, y = ~RT, z = ~Intensity,
								data = sub_peaks,
								color=I("black"),
								size = I(1),
								name = "",
								error_z=list(
									color="black",
									thickness=0,
									symmetric = TRUE, 
									type = "data" ,
									array = sub_peaks$Intensity
								)
							)	
						})							
					}else{
						output$plot_peaks_3D <- renderPlotly({plotly:::plot_ly(type="scatter3d",mode="markers",showlegend=TRUE)})					
					}
				}						
				# nothing? ##################################################################
				if( !isolate(input$peaks_mz_RT_use_peaks) & !isolate(input$peaks_mz_RT_use_raw)){
					output$plot_peaks_3D <- renderPlotly({plotly:::plot_ly(type="scatter3d",mode="markers",showlegend=TRUE)})	
				}
				#############################################################################
			}else{
				output$plot_peaks_mz_RT <- renderPlot({})
				output$plot_peaks_mz_int <- renderPlot({})
				output$plot_peaks_RT_int <- renderPlot({})
				output$plot_peaks_3D <- renderPlotly({plotly:::plot_ly(type="scatter3d",mode="markers",showlegend=TRUE)})	
			}
		}
	}
})

# When a double-click happens, check if there's a brush on the plot.
# If so, zoom to the brush bounds; if not, reset the zoom.
observeEvent(input$plot_peaks_mz_RT_dblclick, { 
    brush <- input$plot_peaks_mz_RT_brush
    if (!is.null(brush)) {
		isolate(ranges_peaks_mz_RT$x <- c(brush$xmin, brush$xmax))
		isolate(ranges_peaks_mz_RT$y <- c(brush$ymin, brush$ymax))
    } else {
		isolate(ranges_peaks_mz_RT$x <- NULL)
		isolate(ranges_peaks_mz_RT$y <- NULL)
    }
	refresh_plot$a<-(refresh_plot$a+1)
	cat("\n Zooming with brush")
})

observe({
	input$peaks_mz_RT_zoom_out
	if(isolate(init$a)=="TRUE"){
		if(!is.na(isolate(input$sel_meas_ID))){ 
			if(!is.null(isolate(ranges_peaks_mz_RT$x))){
				cat("\n Zooming out on X")
				old_range<-abs(isolate(ranges_peaks_mz_RT$x[2]-ranges_peaks_mz_RT$x[1]))
				isolate(ranges_peaks_mz_RT$x[1]<-ranges_peaks_mz_RT$x[1]-.5*old_range)
				isolate(ranges_peaks_mz_RT$x[2]<-ranges_peaks_mz_RT$x[2]+.5*old_range)
			}
			
			
		refresh_plot$a<-(refresh_plot$a+1)		
		}
	}
})

observe({
    input$sel_meas_ID
	input$sel_peak_ID
	if(!is.na(isolate(input$sel_meas_ID))){ # if user deletes entry!
		if(isolate(input$sel_meas_ID)!=0){
			if(atit==isolate(input$sel_meas_ID)){ # above MSlist upload worked?
				if( !is.na(isolate(input$sel_peak_ID)) & 
					(isolate(input$sel_peak_ID)!=0) & 
					any(MSlist[[8]][,10]==isolate(input$sel_peak_ID)) &
					any(objects(envir=as.environment(".GlobalEnv"))=="MSlist")
				){
					EIC_ID<<-unique(MSlist[[8]][MSlist[[8]][,10]==isolate(input$sel_peak_ID),9]);
					peakit<<-MSlist[[4]][[2]][c(MSlist[[7]][isolate(input$sel_peak_ID),1]:MSlist[[7]][isolate(input$sel_peak_ID),2]),]			
					if(length(peakit)>7){
						EICit<<-MSlist[[4]][[2]][c(MSlist[[6]][EIC_ID,1]:MSlist[[6]][EIC_ID,2]),]
						output$EIC1 <- renderPlot({
							if(length(EICit)>7){
								plot(EICit[,3],EICit[,2],type="h",col="darkgrey",xlab="RT",ylab="Intensity",main="EIC (grey) & Peak (red)",xlim=c(min(MSlist[[4]][[1]]),max(MSlist[[4]][[1]])))
							}else{
								plot(EICit[3],EICit[2],type="h",col="darkgrey",xlab="RT",ylab="Intensity",main="EIC (grey) & Peak (red)")
							}
							if(length(peakit)>7){	
								points(peakit[,3],peakit[,2],type="h",col="red",lwd=2)
							}else{
								points(peakit[3],peakit[2],type="h",col="red",lwd=2)				
							}
						})	
						output$EIC2 <- renderPlot({
							if(length(EICit)>7){
								plot(EICit[,3],EICit[,2],type="h",col="darkgrey",xlab="RT",ylab="Intensity")
							}else{
								plot(EICit[3],EICit[2],type="h",col="darkgrey",xlab="RT",ylab="Intensity")
							}
							if(length(peakit)>7){	
								points(peakit[,3],peakit[,2],type="h",col="red",lwd=2)
							}else{
								points(peakit[3],peakit[2],type="h",col="red",lwd=2)				
							}
						})	
						output$EIC3 <- renderPlot({
							if(length(EICit)>7){
								plot(EICit[,3],EICit[,1],pch=19,col="darkgrey",xlab="RT",ylab="m/z")			
							}else{
								plot(EICit[3],EICit[1],pch=19,col="darkgrey",xlab="RT",ylab="m/z")
							}
							if(length(peakit)>7){	
								points(peakit[,3],peakit[,1],pch=19,col="red",cex=1.5)
							}else{
								points(peakit[,3],peakit[,1],pch=19,col="red",cex=1.5)				
							}
						})	
						cat("\n EIC & peak extracted")
					}else{
						cat("\n Peak based on single measurement - plotting skipped.")
					}
				}	
			}else{
				output$EIC1 <- renderPlot({plot.new()})
				output$EIC2 <- renderPlot({plot.new()})
				output$EIC3 <- renderPlot({plot.new()})
			}
		}
	}
})
##############################################################################

##############################################################################
# update results for changes in ion mode selection ###########################
##############################################################################
maincalc3<-reactive({
	input$Ion_mode
	if( (isolate(init$a)=="TRUE") & (isolate(input$Ion_mode)=="positive") ){
		exprprofnorm_pos<-list(src=file.path(logfile[[1]],"pics","profnorm_pos"))
		output$profnorm<-renderImage(exprprofnorm_pos, deleteFile = FALSE)
		exprprofcount_pos<-list(src=file.path(logfile[[1]],"pics","profcount_pos"))
		output$profcount<-renderImage(exprprofcount_pos, deleteFile = FALSE)
		if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList")){ rm(profileList, inherits = TRUE) }
		if(file.exists(file.path(logfile[[1]],"results","profileList_pos"))){
			if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_pos")){rm(profileList_pos,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="profileList_pos")){rm(profileList_pos)}				
			load(file=file.path(as.character(logfile[[1]]),"results","profileList_pos"),envir=as.environment(".GlobalEnv"), verbose=TRUE);
			assign("profileList",profileList_pos,envir=as.environment(".GlobalEnv"));
		}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="profpeaks")){ rm(profpeaks, inherits = TRUE) }
		if(file.exists(file.path(logfile[[1]],"results","profpeaks_pos"))){
			if(any(objects(envir=as.environment(".GlobalEnv"))=="profpeaks_pos")){rm(profpeaks_pos,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="profpeaks_pos")){rm(profpeaks_pos)}				
			load(file=file.path(as.character(logfile[[1]]),"results","profpeaks_pos"),envir=as.environment(".GlobalEnv"), verbose=TRUE);
			assign("profpeaks",profpeaks_pos,envir=as.environment(".GlobalEnv"));
		}
		expr4p<-list(src=file.path(logfile[[1]],"pics","boxprofile_pos"))
		output$boxprofile<-renderImage(expr4p, deleteFile = FALSE)		
		isolate(init$b<<-(init$b+1))
		if(any(objects()=="profileList")){stop("illegal profpeaks2 found, #1");}
		if(any(objects()=="profpeaks")){stop("illegal profpeaks found, #1");}
		return("Select ionization (switch to negative):\n")
	}
	if( (isolate(init$a)=="TRUE") &  (isolate(input$Ion_mode)=="negative") ){
		exprprofnorm_neg<-list(src=file.path(logfile[[1]],"pics","profnorm_neg"))
		output$profnorm<-renderImage(exprprofnorm_neg, deleteFile = FALSE)
		exprprofcount_neg<-list(src=file.path(logfile[[1]],"pics","profcount_neg"))
		output$profcount<-renderImage(exprprofcount_neg, deleteFile = FALSE)
		if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList")){ rm(profileList, inherits = TRUE) }
		if(file.exists(file.path(logfile[[1]],"results","profileList_neg"))){
			if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_neg")){rm(profileList_neg,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="profileList_neg")){rm(profileList_neg)}				
			load(file=file.path(as.character(logfile[[1]]),"results","profileList_neg"),envir=as.environment(".GlobalEnv"), verbose=TRUE);
			assign("profileList",profileList_neg,envir=as.environment(".GlobalEnv"));
		}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="profpeaks")){ rm(profpeaks, inherits = TRUE) }
		if(file.exists(file.path(logfile[[1]],"results","profpeaks_neg"))){
			if(any(objects(envir=as.environment(".GlobalEnv"))=="profpeaks_neg")){rm(profpeaks_neg,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="profpeaks_neg")){rm(profpeaks_neg)}				
			load(file=file.path(as.character(logfile[[1]]),"results","profpeaks_neg"),envir=as.environment(".GlobalEnv"), verbose=TRUE);
			assign("profpeaks",profpeaks_neg,envir=as.environment(".GlobalEnv"));
		}
		expr4n<-list(src=file.path(logfile[[1]],"pics","boxprofile_neg"))
		output$boxprofile<-renderImage(expr4n, deleteFile = FALSE)	
		isolate(init$b<<-(init$b+1))
		if(any(objects()=="profileList")){stop("illegal profpeaks2 found, #2");}
		if(any(objects()=="profpeaks")){stop("illegal profpeaks found, #2");}
		return("Select ionization (switch to positive):\n")	
	}
})
output$had_ion<-renderText(paste(maincalc3())) 
##############################################################################

##############################################################################
# Sort and filter the profile list ###########################################
##############################################################################
maincalc6<-reactive({
	init$a # in project?
	init$b # number of calculations - update profpeaks2 after each such
	#cat(" \n ... ");print(isolate(init$b));cat(" ... ")
	input$filterProf_maxmass
	input$filterProf_minmass
	input$filterProf_minrt
	input$filterProf_maxrt
	input$filterProf_meanblind
	input$filterProf_notblind
	input$filterProf_sort
	input$filterProf_count
    if( 
		(isolate(init$a)=="TRUE") & 
		(any(objects(envir=as.environment(".GlobalEnv"))=="profpeaks")) & 
		!is.na(isolate(input$filterProf_minmass)) & 
		!is.na(isolate(input$filterProf_maxmass)) & 
		!is.na(isolate(input$filterProf_minrt)) & 
		!is.na(isolate(input$filterProf_maxrt)) 
	){
		cat("\n profilepeaks filtered and sorted")		
		if(any(objects(envir=as.environment(".GlobalEnv"))=="profpeaks2")){rm(profpeaks2,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="profileList")){stop("illegal profpeaks2 found, #3");}
		if(any(objects()=="profpeaks")){stop("illegal profpeaks found, #3");}
		if(any(objects()=="profpeaks2")){stop("illegal profpeaks2 found, #3");}
		assign("profpeaks2",profpeaks,envir=as.environment(".GlobalEnv"));
		if( length(profpeaks2)>13 ){profpeaks2<<-profpeaks2[profpeaks2[,1]>=isolate(input$filterProf_minmass),,drop = FALSE]}else{ if( length(profpeaks2)==13 ){ profpeaks2<<-profpeaks2[profpeaks2[1]>=isolate(input$filterProf_minmass),drop = FALSE] }}
		if( length(profpeaks2)>13 ){profpeaks2<<-profpeaks2[profpeaks2[,1]<=isolate(input$filterProf_maxmass),,drop = FALSE]}else{ if( length(profpeaks2)==13 ){  profpeaks2<<-profpeaks2[profpeaks2[1]<=isolate(input$filterProf_maxmass),drop = FALSE] }}
		if( length(profpeaks2)>13 ){profpeaks2<<-profpeaks2[profpeaks2[,3]>=isolate(input$filterProf_minrt),,drop = FALSE]}else{ if( length(profpeaks2)==13 ){  profpeaks2<<-profpeaks2[profpeaks2[3]>=isolate(input$filterProf_minrt),drop = FALSE] }}
		if( length(profpeaks2)>13 ){profpeaks2<<-profpeaks2[profpeaks2[,3]<=isolate(input$filterProf_maxrt),,drop = FALSE]}else{ if( length(profpeaks2)==13 ){  profpeaks2<<-profpeaks2[profpeaks2[3]<=isolate(input$filterProf_maxrt),drop = FALSE] }}
		if( length(profpeaks2)>13 ){
			if(isolate(input$filterProf_meanblind)=="yes"){
				profpeaks2<<-profpeaks2[(profpeaks2[,6]==1),,drop = FALSE] # above blind OR not in blind, = profileList[[7]][k,9]
			}
		}else{
			if( length(profpeaks2)==13 ){
				profpeaks2<<-profpeaks2[(profpeaks2[,6]==1)]
			}
		}
		if( length(profpeaks2)>13 ){
			if(isolate(input$filterProf_notblind)=="yes"){
				profpeaks2<<-profpeaks2[profpeaks2[,5]==0,,drop = FALSE] # not in blind, = profileList[[7]][k,8]
			}
		}else{
			if( length(profpeaks2)==13 ){
				profpeaks2<<-profpeaks2[(profpeaks2[,5]==0)]
			}		
		}
		if( length(profpeaks2)>13 ){
			if(isolate(input$filterProf_sort)=="ID"){
				profpeaks2<<-profpeaks2[order(profpeaks2[,10],decreasing=FALSE),]
			}
			if(isolate(input$filterProf_sort)=="mean m/z"){
				profpeaks2<<-profpeaks2[order(profpeaks2[,1],decreasing=FALSE),]
			}
			if(isolate(input$filterProf_sort)=="mean RT"){
				profpeaks2<<-profpeaks2[order(profpeaks2[,3],decreasing=FALSE),]
			}
			if(isolate(input$filterProf_sort)=="maximum intensity"){
				profpeaks2<<-profpeaks2[order(profpeaks2[,4],decreasing=TRUE),]
			}
			if(isolate(input$filterProf_sort)=="total peak number"){
				profpeaks2<<-profpeaks2[order(profpeaks2[,11],profpeaks2[,4],decreasing=TRUE),]
			}			
			if(isolate(input$filterProf_sort)=="mean intensity"){
				profpeaks2<<-profpeaks2[order(profpeaks2[,2],decreasing=TRUE),]
			}
			if(isolate(input$filterProf_sort)=="global trend intensity"){
				profpeaks2<<-profpeaks2[order(profpeaks2[,12],decreasing=TRUE),]
				profpeaks2<<-profpeaks2[profpeaks2[,12]!=0,]
			}
			if(isolate(input$filterProf_sort)=="current trend intensity"){
				profpeaks2<<-profpeaks2[order(profpeaks2[,13],decreasing=TRUE),]
				profpeaks2<<-profpeaks2[profpeaks2[,13]!=0,]
			}
		}
		if(length(profpeaks2)!=0){
			if(length(profpeaks2)==13){ # contains just one row?
				if(
					profpeaks2[1]>=isolate(input$filterProf_minmass) &
					profpeaks2[1]<=isolate(input$filterProf_maxmass) &
					profpeaks2[3]>=isolate(input$filterProf_minrt) &
					profpeaks2[3]<=isolate(input$filterProf_maxrt) 
				){ 
					if( ((isolate(input$filterProf_meanblind)=="yes") & (profpeaks2[6]==1)) || (isolate(input$filterProf_meanblind)=="no") ){
						if( ((isolate(input$filterProf_notblind)=="yes") & (profpeaks2[5]==1)) || (isolate(input$filterProf_notblind)=="no") ){
							output$allproftable<-renderTable(profpeaks2)
							#output$atprof1<-renderText({ "1" })
							output$atprof2<-renderText({ "1" })
							output$atprof3<-renderText({ "1" })
							output$atprof4<-renderText({ "1" })
							output$atprof5<-renderText({ "1" })			
							path=file.path(logfile[[1]],"pics","profilehisto.png");
							png(filename = path, bg = "white", width = 1100);
								plot.new()
								plot.window(xlim=c(0,1),ylim=c(0,1))
								text(0.5,0.5,labels="0 profiles left for these filter settings",cex=1.8,col="red")
							dev.off();
							expr6<-list(src=file.path(logfile[[1]],"pics","profilehisto.png"));
							output$profilehisto<-renderImage(expr6, deleteFile = FALSE);
							updateNumericInput(session,"profID",value = 0);
							updateNumericInput(session,"profentry",value = 0);	
							return("1")
						}else{
							#output$atprof1<-renderText({ "0" })
							output$atprof2<-renderText({ "0" })
							output$atprof3<-renderText({ "0" })
							output$atprof4<-renderText({ "0" })
							output$atprof5<-renderText({ "0" })
							path=file.path(logfile[[1]],"pics","profilehisto.png");
							png(filename = path, bg = "white", width = 1100);
								plot.new()
								plot.window(xlim=c(0,1),ylim=c(0,1))
								text(0.5,0.5,labels="0 profiles left for these filter settings",cex=1.8,col="red")
							dev.off();
							expr6<-list(src=file.path(logfile[[1]],"pics","profilehisto.png"));
							output$profilehisto<-renderImage(expr6, deleteFile = FALSE);
							output$allproftable<-renderText("No profiles left")
							updateNumericInput(session,"profID",value = 0);
							updateNumericInput(session,"profentry",value = 0);
							return("0")
						}
					}else{
						#output$atprof1<-renderText({ "0" })
						output$atprof2<-renderText({ "0" })
						output$atprof3<-renderText({ "0" })
						output$atprof4<-renderText({ "0" })
						output$atprof5<-renderText({ "0" })
						path=file.path(logfile[[1]],"pics","profilehisto.png");
						png(filename = path, bg = "white", width = 1100);
							plot.new()
							plot.window(xlim=c(0,1),ylim=c(0,1))
							text(0.5,0.5,labels="0 profiles for these filter settings",cex=1.8,col="red")
						dev.off();
						expr6<-list(src=file.path(logfile[[1]],"pics","profilehisto.png"));
						output$profilehisto<-renderImage(expr6, deleteFile = FALSE);
						output$allproftable<-renderText("No profiles left")
						updateNumericInput(session,"profID",value = 0);
						updateNumericInput(session,"profentry",value = 0);
						return("0")
					}
				}else{
					#output$atprof1<-renderText({ "0" })
					output$atprof2<-renderText({ "0" })
					output$atprof3<-renderText({ "0" })
					output$atprof4<-renderText({ "0" })
					output$atprof5<-renderText({ "0" })
					path=file.path(logfile[[1]],"pics","profilehisto.png");
					png(filename = path, bg = "white", width = 1100);
						plot.new()
						plot.window(xlim=c(0,1),ylim=c(0,1))
						text(0.5,0.5,labels="0 profiles for these filter settings",cex=1.8,col="red")
					dev.off();
					expr6<-list(src=file.path(logfile[[1]],"pics","profilehisto.png"));
					output$profilehisto<-renderImage(expr6, deleteFile = FALSE);
					output$allproftable<-renderText("No profiles left")
					updateNumericInput(session,"profID",value = 0);
					updateNumericInput(session,"profentry",value = 0);
					return("0")
				}
			}else{
				if(length(profpeaks2)>13){
					# summary
					atit1<-sum(profpeaks2[,11]) 
					#output$atprof1<-renderText({ atit1 })
					atit2<-length(profpeaks2[,11])
					output$atprof2<-renderText({ atit2 })
					atit3<-length(profpeaks2[profpeaks2[,5]==1,11])
					output$atprof3<-renderText({ atit3 })
					atit4<-length(profpeaks2[profpeaks2[,6]==1,11])
					output$atprof4<-renderText({ atit4 })
					atit5<-length(profpeaks2[profpeaks2[,13]!=0,11])
					output$atprof5<-renderText({ atit5 })
					# intensity histogram
					path=file.path(logfile[[1]],"pics","profilehisto.png");
                    png(filename = path, bg = "white", width = 600);
                    plot_profiles_intensity_histograms(mean_intensities=profpeaks2[,2],
                                                       max_intensities=profpeaks2[,4],
                                                       past_incidents=profpeaks2[,12],
                                                       current_incidents=profpeaks2[,13]);
                    dev.off();
					expr6<-list(src=file.path(logfile[[1]],"pics","profilehisto.png"));
					output$profilehisto<-renderImage(expr6, deleteFile = FALSE);
					# table
					if( (length(profpeaks2[,1])>isolate(input$filterProf_count))  &  !is.na(isolate(input$filterProf_count)) ){
						profpeaks2<<-profpeaks2[1:isolate(input$filterProf_count),]
					}
					profpeaks2<<-as.data.frame(profpeaks2)
					profpeaks2[,1]<<-format(profpeaks2[,1],digits=8)
					profpeaks2[,2]<<-format(profpeaks2[,2],scientific=TRUE,digits=2)
					profpeaks2[,4]<<-format(profpeaks2[,4],scientific=TRUE,digits=2)
					profpeaks2[,5]<<-as.integer(profpeaks2[,5])
					profpeaks2[,6]<<-as.integer(profpeaks2[,6])
					profpeaks2[,7]<<-format(profpeaks2[,7],scientific=TRUE,digits=2)
					profpeaks2[,10]<<-as.integer(profpeaks2[,10])
					profpeaks2[,11]<<-as.integer(profpeaks2[,11])
					profpeaks2[,12]<<-format(profpeaks2[,12],scientific=TRUE,digits=2)
					profpeaks2[,13]<<-format(profpeaks2[,13],scientific=TRUE,digits=2)
					profpeaks3<<-profpeaks2[,c(10:13,1:9)]
					names(profpeaks3)<<-c("profile ID","number of peaks","global trend intensity","current trend intensity","mean m/z", "mean intensity", "mean RT", "maximum Intensity", "in blind?", "above blind?", "m/z variance", "minimum RT", "maximum RT")
					output$allproftable<-renderTable(profpeaks3)
					updateNumericInput(session,"profID",value = 0);
					updateNumericInput(session,"profentry",value = 0);
					return(as.character(atit1));
				}else{
					#output$atprof1<-renderText({ "0" })
					output$atprof2<-renderText({ "0" })
					output$atprof3<-renderText({ "0" })
					output$atprof4<-renderText({ "0" })
					output$atprof5<-renderText({ "0" })
					path=file.path(logfile[[1]],"pics","profilehisto.png");
					png(filename = path, bg = "white", width = 1100);
						plot.new()
						plot.window(xlim=c(0,1),ylim=c(0,1))
						text(0.5,0.5,labels="0 profiles for these filter settings",cex=1.8,col="red")
					dev.off();
					expr6<-list(src=file.path(logfile[[1]],"pics","profilehisto.png"));
					output$profilehisto<-renderImage(expr6, deleteFile = FALSE);
					output$allproftable<-renderText("No profiles left")
					updateNumericInput(session,"profID",value = 0);
					updateNumericInput(session,"profentry",value = 0);
					return("0")
				}
			}
		}else{
			#output$atprof1<-renderText({ "0" })
			output$atprof2<-renderText({ "0" })
			output$atprof3<-renderText({ "0" })
			output$atprof4<-renderText({ "0" })
			output$atprof5<-renderText({ "0" })
			path=file.path(logfile[[1]],"pics","profilehisto.png");
			png(filename = path, bg = "white", width = 1100);
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,labels="0 profiles for these filter settings",cex=1.8,col="red")
			dev.off()
			expr6<-list(src=file.path(logfile[[1]],"pics","profilehisto.png"))
			output$profilehisto<-renderImage(expr6, deleteFile = FALSE);
			output$allproftable<-renderText("No profiles left")
			updateNumericInput(session,"profID",value = 0)
			updateNumericInput(session,"profentry",value = 0)
			return("0")
		}
	}else{
		if( isolate(init$a)=="TRUE" ){
			cat("\n No profiles available\n")
			#output$atprof1<-renderText({"0"}) # now used as reactive output
			output$atprof2<-renderText({"0"})
			output$atprof3<-renderText({"0"})
			output$atprof4<-renderText({"0"})
			output$atprof5<-renderText({"0"})	
			output$allproftable<-renderText("No profiles available")
			path=file.path(logfile[[1]],"pics","profilehisto.png")
			png(filename = path, bg = "white", width = 1100)
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,labels="0 profiles for these filter settings",cex=1.8,col="red")
			dev.off()
			expr6<-list(src=file.path(logfile[[1]],"pics","profilehisto.png"))
			output$profilehisto<-renderImage(expr6, deleteFile = FALSE)
			updateNumericInput(session,"profID",value = 0)
			updateNumericInput(session,"profentry",value = 0)
			return("0")
		}
	}
	if(any(objects()=="profileList")){stop("illegal profpeaks2 found, #4");}
	if(any(objects()=="profpeaks")){stop("illegal profpeaks found, #4");}
	if(any(objects()=="profpeaks2")){stop("illegal profpeaks2 found, #4");}
	if(any(objects()=="profpeaks3")){stop("illegal profpeaks3 found, #4");}
})	
output$atprof1<-renderText(paste(maincalc6()))
output$peak_number<-renderText(paste(maincalc6())) 
##############################################################################

##############################################################################
# update results for individual profileIDs ###################################
##############################################################################
ranges_timeprofile <- reactiveValues(x = NULL, y = NULL)

observe({
    input$profID
	init$b
    if( (isolate(init$a)=="TRUE") &  
		(!is.na(isolate(input$profID))) & 
		(isolate(input$profID)!=0) & 
		any(objects(envir=as.environment(".GlobalEnv"))=="profileList") & 
		any(objects(envir=as.environment(".GlobalEnv"))=="profpeaks") 	& 
		any(objects(envir=as.environment(".GlobalEnv"))=="profpeaks2")
	){
		if(any(objects()=="profileList")){stop("illegal profileList found, #5");}
		if(any(objects()=="profpeaks")){stop("illegal profpeaks found, #5");}
		if(any(objects()=="profpeaks2")){stop("illegal profpeaks2 found, #5");}
		if(any(objects()=="profpeaks3")){stop("illegal profpeaks3 found, #5");}
		if(any(profileList[[7]][,4]==as.numeric(isolate(input$profID)))){
			cat("\n plotting profile with ID ");cat(as.numeric(isolate(input$profID)));
			if(logfile$parameters$trend_blind=="yes"){
				blindsubtract<-TRUE
			}else{
				blindsubtract<-FALSE
			}
			lagit<-as.numeric(strsplit(logfile$parameters$trend_lags,",")[[1]])		
			if(isolate(input$prof_log)=="yes"){
				logscaled<-TRUE
			}else{
				logscaled<-FALSE
			}
			output$timeprofile <- renderPlot({			
				assign("peakTable",plotaprofile(
					profileList,
					profileID=as.numeric(isolate(input$profID)),
					logint=logscaled,
					blindsub=blindsubtract,
					blindfold=as.numeric(logfile$parameters$blind_threshold),
					lags=lagit,
					threshold=as.numeric(logfile$parameters$trend_threshold),
					ranges_x=ranges_timeprofile$x,
					ranges_y=ranges_timeprofile$y,
					),envir=as.environment(".GlobalEnv")
				)	
				if((!is.null(ranges_timeprofile$x))||(!is.null(ranges_timeprofile$y))){
					mtext("Now zoomed in",side=3,col="gray")
				}
			})			
			output$oneproftable<-DT::renderDataTable(peakTable);
			updateNumericInput(session,"profpeakID",value = 0);
			path=file.path(logfile[[1]],"pics","massdens.png");
			png(filename = path, bg = "white", width = 550,height=200);			
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,labels="Waiting...",cex=1.8,col="red")
			dev.off();
			expr_massdens<-list(src=file.path(logfile[[1]],"pics","massdens.png"));
			output$massdens<-renderImage(expr_massdens, deleteFile = FALSE);				
			path=file.path(logfile[[1]],"pics","mass_int.png");
			png(filename = path, bg = "white", width = 550,height=300);	
				plot.new()
			dev.off();
			expr_mass_int<-list(src=file.path(logfile[[1]],"pics","mass_int.png"));
			output$massint<-renderImage(expr_mass_int, deleteFile = FALSE);	
		}else{
			output$timeprofile <- renderPlot({	
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,labels="Nothing to plot - invalid ID",cex=1.8,col="red")
			})
			output$oneproftable<-renderText("")
			updateNumericInput(session,"profpeakID",value = 0);		
			path=file.path(logfile[[1]],"pics","massdens.png");
			png(filename = path, bg = "white", width = 550,height=200);			
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,labels="Not available",cex=1.8,col="red")
			dev.off();
			expr_massdens<-list(src=file.path(logfile[[1]],"pics","massdens.png"));
			output$massdens<-renderImage(expr_massdens, deleteFile = FALSE);				
			path=file.path(logfile[[1]],"pics","mass_int.png");
			png(filename = path, bg = "white", width = 550,height=300);	
				plot.new()
			dev.off();
			expr_mass_int<-list(src=file.path(logfile[[1]],"pics","mass_int.png"));
			output$massint<-renderImage(expr_mass_int, deleteFile = FALSE);	
		}
	}else{
		if(isolate(init$a)=="TRUE"){
			output$timeprofile <- renderPlot({	
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,labels="Nothing to plot - invalid ID",cex=1.8,col="red")
			})			
			output$oneproftable<-renderText("")
			updateNumericInput(session,"profpeakID",value = 0);
			path=file.path(logfile[[1]],"pics","massdens.png");
			png(filename = path, bg = "white", width = 550,height=200);			
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,labels="Not available",cex=1.8,col="red")
			dev.off();
			expr_massdens<-list(src=file.path(logfile[[1]],"pics","massdens.png"));
			output$massdens<-renderImage(expr_massdens, deleteFile = FALSE);			
			path=file.path(logfile[[1]],"pics","mass_int.png");
			png(filename = path, bg = "white", width = 550,height=300);	
				plot.new()
			dev.off();
			expr_mass_int<-list(src=file.path(logfile[[1]],"pics","mass_int.png"));
			output$massint<-renderImage(expr_mass_int, deleteFile = FALSE);	
		}
	}
})	
##############################################################################

##############################################################################
# update results per profilepeak list entry index ############################
##############################################################################
observe({
    input$profentry
	init$b
    if(	(isolate(init$a)=="TRUE") &  
		!is.na(isolate(input$profentry)) & 
		(isolate(input$profentry)!=0) & 
		(any(objects(envir=as.environment(".GlobalEnv"))=="profileList")) & 
		any(objects(envir=as.environment(".GlobalEnv"))=="profpeaks") & 
		any(objects(envir=as.environment(".GlobalEnv"))=="profpeaks2")
	){
		if( (isolate(input$profentry)<=length(profpeaks2[,1])) & 
			(isolate(input$profentry)>0) & 
			(isolate(input$profentry)<=length(profpeaks2[,1])) 
		){ 
				if( any( profileList[[7]][,4]==as.numeric(profpeaks2[isolate(input$profentry),10]) ) ){
					updateNumericInput(session,"profID",value = as.numeric(as.character(profpeaks2[isolate(input$profentry),10])))				
				}
		}else{
			output$timeprofile <- renderPlot({	
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,labels="Invalid list entry",cex=1.8,col="red")
			})			
			output$oneproftable<-renderText("")
		}
	}else{
		if(isolate(init$a)=="TRUE"){
			output$timeprofile <- renderPlot({	
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,labels="Nothing to plot - invalid ID",cex=1.8,col="red")
			})			
			output$oneproftable<-renderText("")
		}
	}
})	
##############################################################################

##############################################################################
# PLOT ZOOM & TABLE SELECTION ################################################
# When a double-click happens, check if there's a brush on the plot.
# If so, zoom to the brush bounds; if not, reset the zoom.
observeEvent(input$timeprofile_dblclick, { # - N
	if(verbose){cat("\n in N")}
    brush <- input$timeprofile_brush
    if (!is.null(brush)) {
		ranges_timeprofile$x <- c(brush$xmin, brush$xmax)
		ranges_timeprofile$y <- c(brush$ymin, brush$ymax)
    } else {
		ranges_timeprofile$x <- NULL
		ranges_timeprofile$y <- NULL
    }
})
##############################################################################

##############################################################################
# get EICs for individual profiles ###########################################
##############################################################################
maincalc4<-reactive({
	input$profpeakID
	if( 	
		(isolate(init$a)=="TRUE") & 
		(!is.na(isolate(input$profID))) & 
		(isolate(input$profpeakID)>0) & 
		(any(objects(envir=as.environment(".GlobalEnv"))=="peakTable")) &
		(any(objects(envir=as.environment(".GlobalEnv"))=="profileList")) 
	){
		if( (isolate(input$profpeakID)<=length(peakTable[,1])) & (any(profileList[[7]][,4]==as.numeric(isolate(input$profID)))) ){
			# positioning plot ###############################################
			path=file.path(logfile[[1]],"pics","profile_position.png");
			png(filename = path, bg = "white", width = 1100,height=150);
				timed<-as.POSIXct(paste(peakTable[,1],peakTable[,2],sep=" "))
				par_old<-par(mar=c(1,1,1,1))				
				plot.new()
				plot.window(xlim=c(min(timed),max(timed)),ylim=c(0,max(max(as.numeric(as.character(peakTable[,4]))),max(as.numeric(as.character(peakTable[,6]))))))
				abline(v=timed[as.numeric(isolate(input$profpeakID))],col="darkgrey",lwd=5)
				box();
				points(timed[peakTable[,3]!=0],peakTable[peakTable[,3]!=0,4],type="l",col="darkgreen");
				points(timed[peakTable[,5]!=0],peakTable[peakTable[,5]!=0,6],type="l",col="red");
				par(par_old);
			dev.off();
			expr_profile_position<-list(src=file.path(logfile[[1]],"pics","profile_position.png"));
			output$profile_position<-renderImage(expr_profile_position, deleteFile = FALSE);		
			# EIC plot ########################################################
			if(any(objects(envir=as.environment(".GlobalEnv"))=="MSlist")){rm(MSlist,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="MSlist")){rm(MSlist)}				
			if(any(objects(envir=as.environment(".GlobalEnv"))=="EIC_ID")){rm(EIC_ID,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="EIC_ID")){rm(EIC_ID)}				
			fileID<-peakTable[,3]
			fileID[fileID=="0"]<-peakTable[peakTable[,3]=="0",5]
			if(any(peakTable[as.numeric(isolate(input$profpeakID)),7]!=0)){
				load(file.path(logfile[[1]],"MSlist",fileID[as.numeric(isolate(input$profpeakID))]), envir=as.environment(".GlobalEnv"))
				cat("\n MSlist loaded");		
				EIC_ID<<-unique(MSlist[[8]][MSlist[[8]][,10]==as.numeric(peakTable[as.numeric(isolate(input$profpeakID)),7]),9]);
				peakit<<-MSlist[[4]][[2]][c(MSlist[[7]][as.numeric(peakTable[as.numeric(isolate(input$profpeakID)),7]),1]:MSlist[[7]][as.numeric(peakTable[as.numeric(isolate(input$profpeakID)),7]),2]),]		
				#if(length(peakit)>7){
				EICit<<-MSlist[[4]][[2]][c(MSlist[[6]][EIC_ID,1]:MSlist[[6]][EIC_ID,2]),]
				path=file.path(logfile[[1]],"pics","profile_EIC.png");
				png(filename = path, bg = "white", width = 1100,height=300);
						par_old<-par(mar=c(2,2,1,1))							
						if(length(EICit)>7){
							plot(EICit[,3],EICit[,2],type="h",col="darkgrey",xlab="RT",ylab="Intensity",xlim=c(min(MSlist[[4]][[1]]),max(MSlist[[4]][[1]])))
						}else{
							plot(EICit[3],EICit[2],type="h",col="darkgrey",xlab="RT",ylab="Intensity")
						}
						if(length(peakit)>7){	
							points(peakit[,3],peakit[,2],type="h",col="red",lwd=2)
						}else{
							points(peakit[3],peakit[2],type="h",col="red",lwd=2)				
						}
						par(par_old);	
				dev.off();
				expr_profile_EIC<-list(src=file.path(logfile[[1]],"pics","profile_EIC.png"));
				output$profile_EIC<-renderImage(expr_profile_EIC, deleteFile = FALSE);	
				rm(EIC_ID,peakit,envir=as.environment(".GlobalEnv"))
				return(
					paste("= sample file ID: ",as.character(fileID[as.numeric(isolate(input$profpeakID))])," (",as.character(timed[as.numeric(isolate(input$profpeakID))]),")" )
				);		
				#}					
			}else{
				path=file.path(logfile[[1]],"pics","profile_EIC.png");
				png(filename = path, bg = "white", width = 1100,height=300);
					plot.new()
					plot.window(xlim=c(0,1),ylim=c(0,1))
					text(0.5,0.5,labels="No sample peak available",cex=1.8,col="red")
				dev.off();
				expr_profile_EIC<-list(src=file.path(logfile[[1]],"pics","profile_EIC.png"));
				output$profile_EIC<-renderImage(expr_profile_EIC, deleteFile = FALSE);		
				return("No peak or EIC available");
			}
		}else{
			path=file.path(logfile[[1]],"pics","profile_position.png");
			png(filename = path, bg = "white", width = 900,height=150);
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,labels="Selection out of range",cex=1.8,col="red")
			dev.off();
			expr_profile_position<-list(src=file.path(logfile[[1]],"pics","profile_position.png"));
			output$profile_position<-renderImage(expr_profile_position, deleteFile = FALSE);		
			path=file.path(logfile[[1]],"pics","profile_EIC.png");
			png(filename = path, bg = "white", width = 1100,height=300);
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,labels="No sample peak available",cex=1.8,col="red")
			dev.off();
			expr_profile_EIC<-list(src=file.path(logfile[[1]],"pics","profile_EIC.png"));
			output$profile_EIC<-renderImage(expr_profile_EIC, deleteFile = FALSE);		
			return("No peak or EIC available");
		}
	}else{
		if( (isolate(init$a)=="TRUE") ){
			path=file.path(logfile[[1]],"pics","profile_position.png");
			png(filename = path, bg = "white", width = 900,height=150);
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,labels="Nothing to plot",cex=1.8,col="red")
			dev.off();
			expr_profile_position<-list(src=file.path(logfile[[1]],"pics","profile_position.png"));
			output$profile_position<-renderImage(expr_profile_position, deleteFile = FALSE);
			path=file.path(logfile[[1]],"pics","profile_EIC.png");
			png(filename = path, bg = "white", width = 1100,height=300);
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,labels="No sample peak available",cex=1.8,col="red")
			dev.off();
			expr_profile_EIC<-list(src=file.path(logfile[[1]],"pics","profile_EIC.png"));
			output$profile_EIC<-renderImage(expr_profile_EIC, deleteFile = FALSE);		
			return("No peak or EIC available");
		}
	}
})
output$prof_peak_text<-renderText(paste(maincalc4())) 
##############################################################################
	
##############################################################################
# get mass estimates for individual profiles #################################
##############################################################################
maincalc5<-reactive({
	input$dens_mass
	if(	(isolate(input$dens_mass)) &
		(isolate(init$a)=="TRUE") & 
		(isolate(input$profID)!=0) & 
		(any(objects(envir=as.environment(".GlobalEnv"))=="profileList")) 
	){
		if(any(profileList[[7]][,4]==as.numeric(isolate(input$profID)))){
			######################################################################
			cat("\n kernel density ...")
			if(isolate(input$use_weight)=="yes"){use_weights<-TRUE}else{use_weights<-FALSE}
			path=file.path(logfile[[1]],"pics","massdens.png");
			png(filename = path, bg = "white", width = 550,height=400);	
				getmass<-mass_dens(
						profileList,
						profileID=as.numeric(isolate(input$profID)),
						bootstrap=TRUE,
						boot_size=as.numeric(isolate(input$boot_size)),
						use_weights
				)
			dev.off();
			expr_massdens<-list(src=file.path(logfile[[1]],"pics","massdens.png"));
			output$massdens<-renderImage(expr_massdens, deleteFile = FALSE);	
			path=file.path(logfile[[1]],"pics","mass_int.png");
			png(filename = path, bg = "white", width = 550,height=300);	
					mass_int(
						profileList,
						profileID=as.numeric(isolate(input$profID))
					)
			dev.off();
			expr_mass_int<-list(src=file.path(logfile[[1]],"pics","mass_int.png"));
			output$massint<-renderImage(expr_mass_int, deleteFile = FALSE);	
			return(getmass);
			######################################################################
		}else{
			return("Not available");
		}
	}else{
		######################################################################
		cat("\n kernel density ...")
		path=file.path(logfile[[1]],"pics","massdens.png");
		png(filename = path, bg = "white", width = 550,height=200);			
			plot.new()
			plot.window(xlim=c(0,1),ylim=c(0,1))
			text(0.5,0.5,labels="Not available",cex=1.8,col="red")
		dev.off();
		expr_massdens<-list(src=file.path(logfile[[1]],"pics","massdens.png"));
		output$massdens<-renderImage(expr_massdens, deleteFile = FALSE);		
		return("...");
		######################################################################
	}
})
output$prof_mass<-renderText(paste(maincalc5())) 
##############################################################################

##############################################################################
# Observe project reset buttons ##############################################
##############################################################################
observe({
    input$reset_1
    if( (isolate(init$a)=="TRUE") & isolate(input$reset_1) ){

		if(any(ls()=="logfile")){stop(paste("\n illegal logfile detected in server_obs_res_mean.r #1"))}
		logfile$Tasks_to_redo<<-replace(logfile$Tasks_to_redo,-1,TRUE)
		logfile$Tasks_to_redo<<-replace(logfile$Tasks_to_redo,1,FALSE)
		measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		measurements[,c("checked","recal","align","norm", "LOD","IS_screen","tar_screen","isotopologues","adducts","homologues","EIC_correlation","blind")]<-"FALSE"
		write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
		createAlert(session,anchorId = "reset", alertId="reset1", title = NULL, content="Project reset w/o peak picking",style = "warning",append=FALSE,dismiss=TRUE)
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		cat("\nReset without peak picking \n")
	}
})
observe({
    input$reset_2
    if( (isolate(init$a)=="TRUE") & isolate(input$reset_2) ){
		if(any(ls()=="logfile")){stop(paste("\n illegal logfile detected in server_obs_res_mean.r #1"))}		
		logfile$Tasks_to_redo<<-replace(logfile$Tasks_to_redo,,TRUE)
		measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		if(TRUE){
			measurements[,c("peakpicking","checked","recal","align","norm", "LOD","IS_screen","tar_screen","isotopologues","adducts","homologues","EIC_correlation","blind")]<-"FALSE"
		}else{ # in-script switch: only redo peak-picking for replicate files
			measurements[measurements[,names(measurements)=="tag3"]!=FALSE,c("peakpicking","checked","recal","align","norm", "LOD","IS_screen","tar_screen","isotopologues","adducts","homologues","EIC_correlation","blind")]<-FALSE
		}
		write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
		rm(measurements)
		# delete all peaklists
		those<-list.files(file.path(logfile$project_folder,"peaklist"))
		if(length(those)>0){
			for(i in 1:length(those)){
				file.remove(file.path(logfile$project_folder,"peaklist",those[i]))
			}
		}
		createAlert(session,anchorId = "reset", alertId="reset2", title = NULL, content="Project reset",style = "warning",append=FALSE,dismiss=TRUE)
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		cat("\nTotal reset \n")
	}
})
##############################################################################

##############################################################################
# Observe resolution data sets ###############################################
##############################################################################
observe({
    input$resolution
	init$a
    if( (isolate(init$a)=="TRUE") & (isolate(input$resolution)!="none") ){
		#cat(input$resolution);cat("\n")
		path=file.path(logfile[[1]],"pics","resolution")
		png(filename = path, bg = "white")
			that<-resolution_list[names(resolution_list) == as.character(input$resolution)][[1]]
			plot(that[,1],that[,2],pch=19,cex=0.5,xlab="m/z",ylab="Resolution")		
		dev.off()
		exprres<-list(src=file.path(logfile[[1]],"pics","resolution"))
		output$plot_resolution<-renderImage(exprres, deleteFile = FALSE)	
	}
})
##############################################################################










