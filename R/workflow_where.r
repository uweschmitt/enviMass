#' @title Find parameters in workflow scripts
#'
#' @export
#'
#' @description Schedules workflow nodes
#'
#' @param strings Vector of character strings with parameter names
#' @param check_mute Check if parameters are muted 
#' 
#' @details enviMass workflow function, detects parameter usage in workflow scripts; returns names affected nodes
#' 

workflow_where<-function(strings,check_mute=TRUE){

	##########################################################################################
	affected<-c()
	files<-list.files()
	for(i in 1:length(files)){
		if(	
			(grepl("do_",files[i])) || (grepl("dont_",files[i])) # workflow scripts		
		){
			textit<-readLines(
				file.path(files[i])
			)
			##################################################################################
			gotit<-FALSE
			for(j in 1:length(textit)){
				for(k in 1:length(strings)){
					if(grepl(strings[k],textit[j], fixed = TRUE)){
						if(check_mute){
							all_muted<-TRUE
							in_here<-gregexpr(strings[k],textit[j])[[1]]
							for(n in 1:length(in_here)){
								if(in_here[n]<6){
									all_muted<-FALSE;
									break;
								}
								if(substr(textit[j],(in_here[n]-5),(in_here[n]-2))!="mute"){
									all_muted<-FALSE
									break;
								}		
							}
							if(all_muted){next} # skip because parametes are muted
						}
						if(grepl("do_",files[i])){
							this<-strsplit(files[i],"do_")[[1]][2]
						}else{
							this<-strsplit(files[i],"dont_")[[1]][2]		
						}
						this<-strsplit(this,".",fixed=TRUE)[[1]][1]
						affected<-c(affected,this)
						gotit<-TRUE;
						break;
					}
				}
				if(gotit){
					break # continue with next file
				}
			}
			##################################################################################	
		}
	}
	##########################################################################################
	affected<-unique(affected)
	return(affected)	
	##########################################################################################
	
}



