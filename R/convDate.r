#' @title Convert dates
#'
#' @description \code{convDate} streamlines date inputs for interlan workflow usage
#'
#' @param measurements enviMass file list
#' 
#' @details enviMass workflow function
#' 

convDate<-function(measurements){
	
	leng<-length(measurements[,1])
	getit<-c()
	for(i in 1:leng){
		dated<-as.character(measurements[i,6])
		done<-FALSE
		if((grepl("[.]",dated))){
			dated<-strsplit(dated,"[.]")[[1]];
			done<-TRUE;
		}
		if(!done){if((grepl("-",dated))){
			dated<-strsplit(dated,"-")[[1]];
			done<-TRUE;
		}}		
		ldated<-nchar(dated)
		if(ldated[1]==max(ldated)){
			Y<-dated[1]
			M<-dated[2]
			D<-dated[3]
		}else{
			Y<-dated[3]
			M<-dated[2]
			D<-dated[1]
		}
		that<-paste(Y,M,D,sep="-");
		getit<-c(getit,that);
	}
	measurements[,6]<-getit;
	return(measurements);
	
}