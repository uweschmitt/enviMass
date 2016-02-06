#' @title Pattern decomposition function 
#'
#' @description Pattern decomposition function 
#'
#' @param cent_peak_mat
#' @param pattern
#' @param profileList
#' @param LOD
#' @param RT_tol_inside
#' 
#' @details enviMass workflow function
#' 

recomb_score<-function(cent_peak_mat,pattern,profileList,LOD,RT_tol_inside){
			
	#######################################################################
	if(!is.matrix(cent_peak_mat)){stop("cent_peak_mat must be a matrix")}
	if(!is.matrix(pattern)){stop("pattern must be a matrix")}
	if((!is.numeric(LOD))||(length(LOD)>1)){stop("LOD must be numeric")}
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
			if( check_plaus(check_nodes[[k]],profileList,RT_tol_inside) ){
				results[[at_results]]<-check_nodes[[k]]
				at_results<-(at_results+1)
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
			}
		}
		check_nodes<-new_nodes
		check_nodes_index<-new_nodes_index
	}
	return(results)
	#######################################################################
	
}
