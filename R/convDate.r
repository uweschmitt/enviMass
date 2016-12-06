#' @title Convert dates
#'
#' @description \code{convDate} streamlines date inputs for workflow usage
#'
#' @param for_dates Vector with date strings
#' 
#' @details enviMass workflow function
#' 

convDate<-function(for_dates){
	
	#if(any(!is.character(for_dates))){stop("\n debug for_dates in convDate: non-character input found!")}
	leng<-length(for_dates)
	getit<-c()
	for(i in 1:leng){
		if(for_dates[i]!="FALSE"){
			dated<-for_dates[i]
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
		}else{
			getit<-c(getit,"FALSE");
		}
	}
	return(getit);
	
}

