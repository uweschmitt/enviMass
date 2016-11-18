
path<-"D:/Users/uchemadmin/Desktop/MS/R_packages/enviMass_devel/R"
files<-list.files(path)

path<-"D:/Users/uchemadmin/Desktop/MS/R_packages/enviMass_devel/inst/webMass"
files<-list.files(path)



for(i in 1:length(files)){
	if(	TRUE 
		#|| (grepl("do_",files[i])) ||
		#(grepl("dont_",files[i])) ||
		#(grepl("server_",files[i]))		
	){
		textit<-readLines(
			file.path(path,files[i])
		)
		##################################################################################
		found<-1
		while(found>0){
			found<-0
			for(j in 1:length(textit)){
				if(grepl("logfile$parameters[",textit[j], fixed = TRUE)){
					found<-(found+1)
					at<-gregexpr(pattern="logfile$parameters[",	textit[j], fixed = TRUE) # get starting position
					from<-at[[1]][1]
					if(substr(textit[j],from+19,from+19)=="["){
						brack<-2
					}else{
						brack<-1
					}
					to<-(from+19-1)
					do_brack<-brack
					while(do_brack>0){ # find end of brackets - end position
						to<-(to+1)
						if(substr(textit[j],to,to)=="]"){
							do_brack<-(do_brack-1)
						}
					}
					to_replace_old<-substr(textit[j],from,to)
					if(brack==1){
						number<-substr(textit[j],(from+19),(to-1))
					}else{
						number<-substr(textit[j],(from+20),(to-2))				
					}
					if(grepl(",",number, fixed = TRUE)){
						found<-(found-1)
						next;
					}
					if(grepl(":",number, fixed = TRUE)){
						found<-(found-1)
						next;
					}
					param_name<-names(logfile$parameters)[as.numeric(number)]	
					to_replace_new<-paste("logfile$parameters$",param_name,sep="")
					textit[j]<-gsub(to_replace_old,to_replace_new,textit[j],fixed = TRUE)
					cat(paste(to_replace_old,"with",to_replace_new));cat("\n")
				}
			}
		}
		##################################################################################
		writeLines(textit,
			file.path(path,files[i])
		)
	}
	#stop()
}






