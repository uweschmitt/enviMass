#' @title Pattern decomposition function 
#'
#' @description Pattern decomposition function 
#'
#' @param cent_peak_mat
#' @param pattern
#' @param profileList
#' @param LOD
#' @param RT_tol_inside
#' @param int_tol
#' @param score_cut 
#' @param plotit
#' 
#' @details enviMass workflow function
#' 

recomb_score<-function(cent_peak_mat,pattern_compound,profileList,LOD,RT_tol_inside,int_tol,score_cut,plotit=FALSE){
			
	#######################################################################
	if(!is.matrix(cent_peak_mat)){stop("cent_peak_mat must be a matrix")}
	if(!is.matrix(pattern_compound)){stop("pattern_compound must be a matrix")}
	if((!is.numeric(LOD))||(length(LOD)>1)){stop("LOD must be numeric")}
	if((!is.numeric(RT_tol_inside))||(length(RT_tol_inside)>1)){stop("LOD must be numeric")}
	if((!is.numeric(int_tol))||(length(int_tol)>1)){stop("LOD must be numeric")}	
	#######################################################################		
	results<-list(0)
	at_results<-1
	checked<-TRUE
	check_nodes<-list()
	check_nodes[[1]]<-cent_peak_mat # initialize with full set
	check_nodes_index<-list()
	check_nodes_index[[1]]<-length(cent_peak_mat[,1])
	while(checked){
		new_nodes<-list()
		new_nodes_index<-list()
		at_new_nodes<-1
		checked<-FALSE
		for(k in 1:length(check_nodes)){
			if(plotit){	
				rescale<-weighted.mean(
					x=(pattern_compound[check_nodes[[k]][,1],2]/profileList[[2]][check_nodes[[k]][,2],2]),
				    w=( profileList[[2]][check_nodes[[k]][,2],2] / (int_tol/100*profileList[[2]][check_nodes[[k]][,2],2]) )
				)
				plot(
					log10(pattern_compound[,2]/rescale),
					log10(pattern_compound[,2]/rescale),
					pch=19,col="lightgray",cex=1.5,
					xlab="Theoretical intensity",
					ylab="Measured intensity"
				)
				abline(h=log10(LOD),col="red")
				abline(v=log10(LOD),col="red")
				abline(0,1,col="lightgrey",lwd=2)
				points(
					log10(pattern_compound[check_nodes[[k]][,1],2]/rescale),
					log10(profileList[[2]][check_nodes[[k]][,2],2]),
					pch=19,col="darkgreen",cex=1
				)
			}			
			if( 
				check_plaus(
					cent_peak_combi=check_nodes[[k]],
					pattern_centro=pattern_compound,
					profileList,
					RT_tol_inside,
					int_tol
				) 
			){
				results[[at_results]]<-list()
				results[[at_results]][[1]]<-check_nodes[[k]]
				# score intensities ...
				rescale<-weighted.mean(
					x=(pattern_compound[check_nodes[[k]][,1],2]/profileList[[2]][check_nodes[[k]][,2],2]),
				    w=( profileList[[2]][check_nodes[[k]][,2],2] / (int_tol/100*profileList[[2]][check_nodes[[k]][,2],2]) )
				)
				above_LOD<-((pattern_compound[,2]/rescale)>LOD)
				# ... measured "above LOD threshold:"
				if(any(above_LOD)){
					score1<-(
						sum(pattern_compound[check_nodes[[k]][,1],2][above_LOD[check_nodes[[k]][,1]]])/sum(pattern_compound[above_LOD,2])
					)
					score1<-round(score1,digits=4)
				}else{
					score1<-NA
				}
				results[[at_results]][[2]]<-score1
				# ... measured "below LOD threshold:"
				if(any(!above_LOD)){
					score2<-(
						sum(pattern_compound[check_nodes[[k]][,1],2][!above_LOD[check_nodes[[k]][,1]]]) / sum(pattern_compound[!above_LOD,2])
					)
				}else{
					score2<-NA
				}
				results[[at_results]][[3]]<-score2
				results[[at_results]][[4]]<-((pattern_compound[check_nodes[[k]][,1],1]-profileList[[2]][check_nodes[[k]][,2],1])/mean(pattern_compound[check_nodes[[k]][,1],1])*1E6)				
				results[[at_results]][[5]]<-(mean(profileList[[2]][check_nodes[[k]][,2],3])-profileList[[2]][check_nodes[[k]][,2],3])
				results[[at_results]][[6]]<-rescale
				names(results[[at_results]])<-c("Peaks","score_1","score_2","ppm deviation","RT deviation from mean","rescale factor")
				at_results<-(at_results+1)
				if(plotit){box(col="green",lwd=5);title(main=paste(score1,score2,sep=" - "));Sys.sleep(3);}
				}else{
				# maker smaller combinations by omission of one (centroid,peak)
				if(check_nodes_index[[k]]>0){ 
					# nothing to inherit - 
					# - this combination is either part of another larger one or
					# - has been build from low-combining that one
					len<-(length(check_nodes[[k]][,1]):1)
					len<-(len[1:check_nodes_index[[k]]])
					len<-rev(len)
					for(z in 1:check_nodes_index[[k]]){
						new_nodes[[at_new_nodes]]<-check_nodes[[k]][-(len[z]),,drop=FALSE]
						new_nodes_index[[at_new_nodes]]<-(check_nodes_index[[k]]-z)
						at_new_nodes<-(at_new_nodes+1)
					}
					checked<-TRUE
				}
				if(plotit){box(col="red",lwd=5);Sys.sleep(2);}
			}
		}
		check_nodes<-new_nodes
		check_nodes_index<-new_nodes_index
	}
	return(results)
	#######################################################################
	
}
