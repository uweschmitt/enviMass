#' @title Execute workflow node
#'
#' @export
#'
#' @description Schedules workflow nodes
#'
#' @param depend logfile matrix downstream dependencies (1) & recalculations (2)
#' @param must logfile matrix upstream musts (1) & unimplemented donts (-1)
#' 
#' @details enviMass workflow function, returns ordered string with workflow steps
#' 

workflow_schedule<-function(depend,must){

	######################################################################################
	######################################################################################
	# insert must into depends - to be save ##############################################
	for(i in 1:length(must[1,])){
		if(sum(must[,i])>0){
			those<-rownames(must)[must[,i]!=0]
			for_that<-colnames(must)[i]
			for(j in 1:length(those)){ # mark dependency in must-matrix
				depend[rownames(depend)==for_that,colnames(depend)==those[j]]<-1
			}
		}
	}
	######################################################################################
	# insert upstream ==2 dependencies as downstream ==1 dependencies ####################
	# just to be save ####################################################################
	# ==2 only relevant for workflow_set, not for the scheduling #########################
	for(i in 1:length(depend[1,])){
		if(any(depend[,i]==2)){
			depend[i,
				depend[,i]==2
			]<-1
		}	
	}
	depend[depend==2]<-0
	######################################################################################
	# derive workflow order from depend ##################################################	
	do_for<-length(depend[1,])
	say<-"ok"
	node<-c()
	node_level<-c()
	for(i in 1:do_for){ # depend be able to derive order in that maximum number of steps ...
		at<-apply(depend,1,sum)
		do_these<-which(at==0)
		if(length(node)>0){
			do_these<-do_these[
				is.na(match(
					colnames(depend)[do_these],node
				))
			]
		}
		if(length(do_these)==0 & !all(at==0)){
			say<-"workflow_scheduler issue#1"
			break;
		}
		node<-c(node,colnames(depend)[do_these])
		node_level<-c(node_level,rep(i,length(do_these)))		
		# clean downstream dependencies
		for(j in 1:length(do_these)){
			depend[do_these[j],]<-0 # downstream dependencies
			depend[,do_these[j]]<-0 # upstream dependencies
		}
	}
	if(any(is.na(match(colnames(depend),node)))){say<-"workflow_scheduler issue#2"} # node missing
	if(any(depend>0) & say!="ok"){say<-"workflow_scheduler issue#3"} 				# node missing
	######################################################################################
	if(say=="ok"){
		return(
			data.frame(node,node_level,stringsAsFactors=TRUE)
		)
	}else{
		return(say)	
	}
}



