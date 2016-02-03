##############################################################################
# update results for individual measurements #################################
##############################################################################
observe({
    input$sel_meas
    if(isolate(input$sel_meas)!="none"){
    pics<-list.files(file.path(logfile[[1]],"pics"))
	# recalibration ##########################################################
    if(
        any(pics==paste("recal_",isolate(input$sel_meas),sep="")) & (logfile$workflow[2]=="yes")
    ){
        expr1<-list(src=file.path(logfile[[1]],"pics",paste("recal_",isolate(input$sel_meas),sep="")))
        output$recal_pic<-renderImage(expr1, deleteFile = FALSE)
	}
    ##########################################################################
    # peak picking ###########################################################
    if(
        any(pics==paste("peakhist_",isolate(input$sel_meas),sep=""))
    ){
        expr2<-list(src=file.path(logfile[[1]],"pics",paste("peakhist_",isolate(input$sel_meas),sep="")))
        output$peakhist_pic<-renderImage(expr2, deleteFile = FALSE)
    }
    if(
        any(pics==paste("peakmzRT_",isolate(input$sel_meas),sep=""))
    ){
        expr3<-list(src=file.path(logfile[[1]],"pics",paste("peakmzRT_",isolate(input$sel_meas),sep="")))
        output$peakmzRT_pic<-renderImage(expr3, deleteFile = FALSE)
    }
    ##########################################################################
    }
})
##############################################################################

##############################################################################
# retrieve peak information ##################################################
##############################################################################
observe({
    input$sel_meas_ID
	input$sel_peak_ID
	if(!any(objects(envir=as.environment(".GlobalEnv"))=="atit")){
		#atit<<-(0)
		assign("atit",0,envir=as.environment(".GlobalEnv"))
	}
    if(isolate(input$sel_meas_ID)!="none"){
		if( (file.exists(file.path(logfile[[1]],"MSlist",isolate(input$sel_meas_ID)))) & 
			(isolate(input$sel_meas_ID)!=atit) 
		){
				if(any(objects(envir=as.environment(".GlobalEnv"))=="MSlist")){rm(MSlist,envir=as.environment(".GlobalEnv"))}
				if(any(objects()=="MSlist")){rm(MSlist)}				
				load(file.path(logfile[[1]],"MSlist",isolate(input$sel_meas_ID)), envir=as.environment(".GlobalEnv"))
				cat("\nafter:\n "); print(gc());	
				cat("\n MSlist file loaded")
				#atit<<-isolate(input$sel_meas_ID)
				assign("atit",isolate(input$sel_meas_ID),envir=as.environment(".GlobalEnv"))
		}
		if( !is.na(isolate(input$sel_peak_ID)) & 
			(isolate(input$sel_peak_ID)!=0) & 
			any(MSlist[[8]][,10]==isolate(input$sel_peak_ID)) &
			any(objects(envir=as.environment(".GlobalEnv"))=="MSlist")
		){
			EIC_ID<<-unique(MSlist[[8]][MSlist[[8]][,10]==isolate(input$sel_peak_ID),9]);
			peakit<<-MSlist[[4]][[2]][c(MSlist[[7]][isolate(input$sel_peak_ID),1]:MSlist[[7]][isolate(input$sel_peak_ID),2]),]			
			if(length(peakit)>7){
				EICit<<-MSlist[[4]][[2]][c(MSlist[[6]][EIC_ID,1]:MSlist[[6]][EIC_ID,2]),]
				path=file.path(logfile[[1]],"pics","EIC1");
				png(filename = path, bg = "white", width = 1100, height= 300)
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
				dev.off();
				expr6<-list(src=file.path(logfile[[1]],"pics","EIC1"));
				output$EIC1<-renderImage(expr6, deleteFile = FALSE);
				path=file.path(logfile[[1]],"pics","EIC2");
				png(filename = path, bg = "white", width = 1100, height= 300);
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
				dev.off();
				expr7<-list(src=file.path(logfile[[1]],"pics","EIC2"));
				output$EIC2<-renderImage(expr7, deleteFile = FALSE);
				path=file.path(logfile[[1]],"pics","EIC3");
				png(filename = path, bg = "white", width = 1100, height= 300);
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
				dev.off();
				expr8<-list(src=file.path(logfile[[1]],"pics","EIC3"));
				output$EIC3<-renderImage(expr8, deleteFile = FALSE);
				cat("\n EIC & peak extracted")
			}else{
				cat("\n Peak based on single measurement - plotting skipped.")
			}
		}	
    }
})
##############################################################################

##############################################################################
# update results for changes ion mode selection ##############################
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
    if( (isolate(init$a)=="TRUE") & 
		(any(objects(envir=as.environment(".GlobalEnv"))=="profpeaks")) & 
		!is.na(isolate(input$filterProf_minmass)) & 
		!is.na(isolate(input$filterProf_maxmass)) & 
		!is.na(isolate(input$filterProf_minrt)) & 
		!is.na(isolate(input$filterProf_maxrt)) 
	){
		cat("\n profilepeaks filtered and sorted")		
		if(any(objects(envir=as.environment(".GlobalEnv"))=="profpeaks2")){rm(profpeaks2,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="profpeaks2")){rm(profpeaks2)}
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
		if(length(profpeaks2)>13){
			if(isolate(input$filterProf_notblind)=="yes"){
				profpeaks2<<-profpeaks2[profpeaks2[,5]==0,,drop = FALSE] # not in blind, = profileList[[7]][k,8]
			}
		}else{
			if( length(profpeaks2)==13 ){
				profpeaks2<<-profpeaks2[(profpeaks2[,5]==0)]
			}		
		}
		if(length(profpeaks2)>13){
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
					profpeaks2<<-profpeaks2[,c(10:13,1:9)]
					names(profpeaks2)<<-c("profile ID","number of peaks","global trend intensity","current trend intensity","mean m/z", "mean intensity", "mean RT", "maximum Intensity", "in blind?", "above blind?", "m/z variance", "minimum RT", "maximum RT")
					output$allproftable<-renderTable(profpeaks2)
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
			dev.off();
			expr6<-list(src=file.path(logfile[[1]],"pics","profilehisto.png"));
			output$profilehisto<-renderImage(expr6, deleteFile = FALSE);
			output$allproftable<-renderText("No profiles left")
			updateNumericInput(session,"profID",value = 0);
			updateNumericInput(session,"profentry",value = 0);
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
			path=file.path(logfile[[1]],"pics","profilehisto.png");
			png(filename = path, bg = "white", width = 1100);
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,labels="0 profiles for these filter settings",cex=1.8,col="red")
			dev.off();
			expr6<-list(src=file.path(logfile[[1]],"pics","profilehisto.png"));
			output$profilehisto<-renderImage(expr6, deleteFile = FALSE);
			updateNumericInput(session,"profID",value = 0);
			updateNumericInput(session,"profentry",value = 0);
			return("0")
		}
	}
})	
output$atprof1<-renderText(paste(maincalc6()))
output$peak_number<-renderText(paste(maincalc6())) 
##############################################################################

##############################################################################
# update results for individual profileIDs ###################################
##############################################################################
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
		if(any(profileList[[7]][,4]==as.numeric(isolate(input$profID)))){
			if(logfile$parameters[[36]]=="yes"){
				blindsubtract<-TRUE
			}else{
				blindsubtract<-FALSE
			}
			lagit<-as.numeric(strsplit(logfile$parameters[[34]],",")[[1]])		
			if(isolate(input$prof_log)=="yes"){
				logscaled<-TRUE
			}else{
				logscaled<-FALSE
			}
			path=file.path(logfile[[1]],"pics","timeprofile");
			png(filename = path, bg = "white", width = 1100);
			assign("peakTable",plotaprofile(
				profileList,
				profileID=as.numeric(isolate(input$profID)),
				logint=logscaled,
				blindsub=blindsubtract,
				blindfold=as.numeric(logfile$parameters[[37]]),
				lags=lagit,
				threshold=as.numeric(logfile$parameters[[35]])
				),envir=as.environment(".GlobalEnv")
			)
			dev.off();
			expr5<-list(src=file.path(logfile[[1]],"pics","timeprofile"));
			output$timeprofile<-renderImage(expr5, deleteFile = FALSE);
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
			path=file.path(logfile[[1]],"pics","timeprofile");
			png(filename = path, bg = "white", width = 1100);
			plot.new()
			plot.window(xlim=c(0,1),ylim=c(0,1))
			text(0.5,0.5,labels="Nothing to plot - invalid ID",cex=1.8,col="red")
			dev.off();
			expr5<-list(src=file.path(logfile[[1]],"pics","timeprofile"));
			output$timeprofile<-renderImage(expr5, deleteFile = FALSE);
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
			path=file.path(logfile[[1]],"pics","timeprofile");
			png(filename = path, bg = "white", width = 1100);
			plot.new()
			plot.window(xlim=c(0,1),ylim=c(0,1))
			text(0.5,0.5,labels="Nothing to plot - invalid ID",cex=1.8,col="red")
			dev.off();
			expr5<-list(src=file.path(logfile[[1]],"pics","timeprofile"));
			output$timeprofile<-renderImage(expr5, deleteFile = FALSE);
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
				if(any(profileList[[7]][,4]==as.numeric(profpeaks2[isolate(input$profentry),1]))){
					updateNumericInput(session,"profID",value = as.numeric(as.character(profpeaks2[isolate(input$profentry),1])))				
				}
		}else{
			path=file.path(logfile[[1]],"pics","timeprofile");
			png(filename = path, bg = "white", width = 1100);
			plot.new()
			plot.window(xlim=c(0,1),ylim=c(0,1))
			text(0.5,0.5,labels="Invalid list entry",cex=1.8,col="red")
			dev.off();
			expr5<-list(src=file.path(logfile[[1]],"pics","timeprofile"));
			output$timeprofile<-renderImage(expr5, deleteFile = FALSE);
			output$oneproftable<-renderText("")
		}
	}else{
		if(isolate(init$a)=="TRUE"){
			path=file.path(logfile[[1]],"pics","timeprofile");
			png(filename = path, bg = "white", width = 1100);
			plot.new()
			plot.window(xlim=c(0,1),ylim=c(0,1))
			text(0.5,0.5,labels="Nothing to plot - invalid ID",cex=1.8,col="red")
			dev.off();
			expr5<-list(src=file.path(logfile[[1]],"pics","timeprofile"));
			output$timeprofile<-renderImage(expr5, deleteFile = FALSE);
			output$oneproftable<-renderText("")
		}
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
		logfile$Tasks_to_redo<-replace(logfile$Tasks_to_redo,-1,TRUE)
		logfile$Tasks_to_redo<<-replace(logfile$Tasks_to_redo,-1,TRUE)
		logfile$Tasks_to_redo<-replace(logfile$Tasks_to_redo,1,FALSE)
		logfile$Tasks_to_redo<<-replace(logfile$Tasks_to_redo,1,FALSE)
		measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		measurements[,c(11,12,13,14)]<-"FALSE"
		write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
		createAlert(session,anchorId = "reset", alertId="reset1", title = NULL, content="Project reset w/o peak picking",style = "warning",append=FALSE,dismiss=TRUE)
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		cat("\nReset without peak picking \n")
	}
})
observe({
    input$reset_2
    if( (isolate(init$a)=="TRUE") & isolate(input$reset_2) ){
		logfile$Tasks_to_redo<-replace(logfile$Tasks_to_redo,,TRUE)
		logfile$Tasks_to_redo<<-replace(logfile$Tasks_to_redo,,TRUE)
		measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		measurements[,c(11,12,13,14)]<-"FALSE"
		write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
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

