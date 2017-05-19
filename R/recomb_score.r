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
#' @param use_score_cut Logical.
#' @param score_cut Numeric.
#' @param plot_it Logical.
#' @param verbose Logical.
#' @param RT_seperate Logical.
#' 
#' @details enviMass workflow function
#' 

recomb_score<-function(
		cent_peak_mat,
		pattern_compound,
		profileList,
		LOD,
		RT_tol_inside,
		int_tol,
		use_score_cut=FALSE,
		score_cut=0,		
		plot_it=FALSE,
		verbose=FALSE,
		RT_seperate=FALSE
	){
			
	#######################################################################
	if(!is.matrix(cent_peak_mat)){stop("cent_peak_mat must be a matrix")}
	if(!is.matrix(pattern_compound)){stop("pattern_compound must be a matrix")}
	if((!is.numeric(LOD))||(length(LOD)>1)){stop("LOD must be numeric")}
	if((!is.numeric(RT_tol_inside))||(length(RT_tol_inside)>1)){stop("LOD must be numeric")}
	if((!is.numeric(int_tol))||(length(int_tol)>1)){stop("LOD must be numeric")}	
	if(any(profileList[[2]][,2]==0)){stop("Data with zereo intensity found - debug")}
	#######################################################################		
	results<-list()
	at_results<-1
	checked<-TRUE
	check_nodes<-list()
	check_nodes_index<-list()
	# Pre-decompose by RT gaps, especially for suspect screening 
	if((length(cent_peak_mat[,2])>1) & RT_seperate){
		cent_peak_mat<-cent_peak_mat[
			order(profileList[[2]][cent_peak_mat[,2],3],decreasing=FALSE)
		,,drop=FALSE]	
		at_RT<-profileList[[2]][cent_peak_mat[,2],3]
		in_node<-rep(1,length(cent_peak_mat[,1]))
		init_node<-1
		for(i in 2:length(cent_peak_mat[,2])){
			if((at_RT[i]-at_RT[i-1])>RT_tol_inside){
				init_node<-(init_node+1)
				in_node[i]<-init_node
			}
		}
		if(init_node>1){
			for(i in 1:max(init_node)){
				check_nodes[[i]]<-cent_peak_mat[in_node==i,,drop=FALSE]
				check_nodes[[i]]<-(check_nodes[[i]][order(check_nodes[[i]][,1],decreasing=FALSE),,drop=FALSE])
				check_nodes_index[[i]]<-sum(in_node==i)
			}
		}else{
			check_nodes[[1]]<-cent_peak_mat # initialize with full set
			check_nodes_index[[1]]<-length(cent_peak_mat[,1])
		}
	}else{
		check_nodes[[1]]<-cent_peak_mat # initialize with full set
		check_nodes_index[[1]]<-length(cent_peak_mat[,1])
	}
	while(checked){
		if(verbose){cat("\n new round:")}
		new_nodes<-list()
		new_nodes_index<-list()
		at_new_nodes<-1
		checked<-FALSE
		for(k in 1:length(check_nodes)){
			if(verbose){cat("\n  :: "); cat(paste(check_nodes[[k]][,2],collapse=","))}
			if(plot_it){	
				rescale<-weighted.mean(
					x=(pattern_compound[check_nodes[[k]][,1],2]/profileList[[2]][check_nodes[[k]][,2],2]),
				    w=(profileList[[2]][check_nodes[[k]][,2],2] / (int_tol/100*profileList[[2]][check_nodes[[k]][,2],2]) )
				)
				plot(
					log10(pattern_compound[,2]/rescale),
					log10(pattern_compound[,2]/rescale),
					pch=19,col="lightgray",cex=1.5,
					xlab="Theoretical intensity",ylab="Measured intensity",
					xlim=c(min(log10(pattern_compound[,2]/rescale)),max(log10(pattern_compound[,2]/rescale))),
					ylim=c(min(log10(profileList[[2]][check_nodes[[k]][,2],2])),max(log10(profileList[[2]][check_nodes[[k]][,2],2])))
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
			# score intensities ...
			rescale<-weighted.mean(
				x=(pattern_compound[check_nodes[[k]][,1],2]/profileList[[2]][check_nodes[[k]][,2],2]),
				w=(profileList[[2]][check_nodes[[k]][,2],2]/(int_tol/100*profileList[[2]][check_nodes[[k]][,2],2]) )
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
			if(use_score_cut=="TRUE"){
				if(is.na(score1)){ 		# ... discarded if below LOD
					if(verbose){cat("- all below LOD.")}
					next;
				}
				if(score1<=score_cut){ 	# 
					if(verbose){cat("- discarded by score.")}				
					next;
				}				
			}
			contained<-FALSE
			if(length(results)>0){
				for(n in 1:length(results)){
					if( !any(is.na(match(check_nodes[[k]][,2],results[[n]][[1]][,2]))) ){
						contained<-TRUE;
						break; # skip - already contained in some larger combination
					}
				}
			}
			if(contained){
				if(verbose){cat("-redundant.")}
				next
			}			
			get_plaus<-check_plaus( # returns TRUE or FALSE
				cent_peak_combi=check_nodes[[k]],
				pattern_centro=pattern_compound,
				profileList,
				RT_tol_inside,
				int_tol
			) 
			if(get_plaus){
				if(verbose){cat("-TRUE")}
				if(length(results)>0){
					for(n in 1:length(results)){
						if( !any(is.na(match(check_nodes[[k]][,2],results[[n]][[1]][,2]))) ){
							contained<-TRUE;
							break; # skip - already contained in some larger combination
						}
					}
				}
				if(verbose){cat("-ok")}
				results[[at_results]]<-list()
				check_nodes[[k]]<-check_nodes[[k]][order(check_nodes[[k]][,1],decreasing=FALSE),,drop=FALSE]
				results[[at_results]][[1]]<-check_nodes[[k]]
				results[[at_results]][[2]]<-score1
				# ... measured "below LOD threshold:"
				if(any(!above_LOD)){
					score2<-(
						#sum(pattern_compound[check_nodes[[k]][,1],2][!above_LOD[check_nodes[[k]][,1]]]) / sum(pattern_compound[!above_LOD,2])
						sum(pattern_compound[check_nodes[[k]][,1],2][!above_LOD[check_nodes[[k]][,1]]]) / sum(pattern_compound[,2])
					)
				}else{
					score2<-NA
				}
				results[[at_results]][[3]]<-score2
				results[[at_results]][[4]]<-((pattern_compound[check_nodes[[k]][,1],1]-profileList[[2]][check_nodes[[k]][,2],1])/mean(pattern_compound[check_nodes[[k]][,1],1])*1E6)				
				results[[at_results]][[5]]<-(mean(profileList[[2]][check_nodes[[k]][,2],3])-profileList[[2]][check_nodes[[k]][,2],3])
				results[[at_results]][[6]]<-rescale
				results[[at_results]][[7]]<-profileList[[2]][check_nodes[[k]][,2],1]
				results[[at_results]][[8]]<-profileList[[2]][check_nodes[[k]][,2],2]
				results[[at_results]][[9]]<-profileList[[2]][check_nodes[[k]][,2],3]
				names(results[[at_results]])<-c("Peaks","score_1","score_2","ppm deviation","RT deviation from mean","rescale factor","m/z","Intensity","RT")
				at_results<-(at_results+1)
				if(plot_it){box(col="green",lwd=5);title(main=paste(score1,score2,k,sep=" - "));Sys.sleep(3);}
			}else{
				if(verbose){cat("-FALSE")}
				# maker smaller combinations by omission of one (centroid,peak)
				if(check_nodes_index[[k]]>0){ 
					# nothing to inherit - 
					# - this combination is either part of another larger one or
					# - has been build from low-combining that one
					len<-(length(check_nodes[[k]][,1]):1)
					len<-(len[1:check_nodes_index[[k]]])
					len<-rev(len)
					for(z in 1:check_nodes_index[[k]]){
						new_nodes[[at_new_nodes]]<-(check_nodes[[k]][-(len[z]),,drop=FALSE])
						#cat(len[z]);cat(" - ")
						new_nodes_index[[at_new_nodes]]<-(check_nodes_index[[k]]-z)
						#cat(check_nodes_index[[k]]-z);cat(" - ")
						at_new_nodes<-(at_new_nodes+1)
					}
					checked<-TRUE
				}
				if(plot_it){box(col="red",lwd=5);title(main=k);Sys.sleep(2);}
			}
		}
		check_nodes<-new_nodes
		check_nodes_index<-new_nodes_index
	
	} # while
	#######################################################################	
	# resort list - largest combinations go first #########################
	results2<-list()
	len<-rep(0,length(results))
	for(k in 1:length(results)){
		len[k]<-length(results[[k]][[1]][,1])
	}
	ord<-order(len)
	for(k in 1:length(results)){	
		results2[[k]]<-results[[ord[k]]]
	}
	results2<-results
	#######################################################################
	return(results2)
	#######################################################################
	
}
