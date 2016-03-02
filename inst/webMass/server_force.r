output$force_workflow <- networkD3:::renderForceNetwork({

	input$saveflow;
	if(isolate(init$a)=="TRUE"){
	
		depends<-logfile$workflow_depend
		src<-c()
		target<-c()
		val<-c()
		for(i in 1:length(depends[,1])){
			if( logfile$workflow[names(logfile$workflow)==rownames(depends)[i]]=="yes" ){	
				for(j in 1:length(depends[1,])){
					if(depends[i,j]==1 & i!=j & 
					logfile$workflow[names(logfile$workflow)==colnames(depends)[j]]=="yes"
					){
						src<-c(src,(i-1))
						target<-c(target,(j-1))
						val<-c(val,length(depends[,1])-sum(depends[j,]))
					}
				}
			}
		}
		MisLinks<-data.frame(src,target,val)
		names(MisLinks)<-c("source","target","value")
		MisLinks<-MisLinks[order(MisLinks[,1],MisLinks[,2]),]
		MisNodes<-data.frame(rownames(depends),rep(0,length(rownames(depends))),rep(2,length(rownames(depends))),stringsAsFactors =FALSE)
		names(MisNodes)<-c("name","group","size")
		for(i in 1:length(logfile$workflow)){
			if(logfile$workflow[i]=="yes"){
				MisNodes[MisNodes[,1]==names(logfile$workflow)[i],2]<-3
			}else{
				MisNodes[MisNodes[,1]==names(logfile$workflow)[i],2]<-4
			}
		}
		
		MisNodes[MisNodes[,1]=="peakpicking",1]<-"Peak picking"
		MisNodes[MisNodes[,1]=="qc",1]<-"Quality control"
		MisNodes[MisNodes[,1]=="recal",1]<-"Mass recalibration"
		MisNodes[MisNodes[,1]=="profiling",1]<-"Profile extraction"		
		MisNodes[MisNodes[,1]=="trendblind",1]<-"Trend detection"
		MisNodes[MisNodes[,1]=="pattern",1]<-"Isotope pattern calculations"
		MisNodes[MisNodes[,1]=="replicates",1]<-"Replicate intersection"
		MisNodes[MisNodes[,1]=="IS_screen",1]<-"Internal standard screening"
		MisNodes[MisNodes[,1]=="target_screen",1]<-"Target screening"
		MisNodes[MisNodes[,1]=="IS_normaliz",1]<-"Internal standard intensity normalization"
		MisNodes[MisNodes[,1]=="norm",1]<-"Median intensity normalization"
		MisNodes[MisNodes[,1]=="blinds",1]<-"Blind subtraction"		
		MisNodes[MisNodes[,1]=="quantification",1]<-"Quantification"				
		MisNodes[MisNodes[,1]=="IS_subtr",1]<-"IS subtraction"	
		MisNodes[MisNodes[,1]=="target_subtr",1]<-"target subtraction"	
		
		networkD3:::forceNetwork(Links = MisLinks, Nodes = MisNodes,
					Source = "source", Target = "target",
					NodeID = "name", Group = "group", 
					zoom=TRUE, opacity = 0.8)
					
					
	}
	
})