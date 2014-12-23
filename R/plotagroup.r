#' @title Plot time series of component profiles
#'
#' @export
#'
#' @description \code{plotagroup} extracts all profiles part of a component with a certain ID and plots them.
#'
#' @param components List of components
#' @param ID Component ID
#' @param profileList A profile list.
#' @param logint Logical. Plot parameter.
#' @param blindsub Logical. Plot parameter.
#' @param blindfold Numeric. Plot parameter.
#' @param lags Numeric vector. Plot parameter.
#' @param threshold Numeric. Plot parameter.
#' @param add Logical. Plot parameter.
#' @param textit Logical. Plot parameter.
#' 
#' @details enviMass workflow plot function
#' 


plotagroup<-function(
	components,
	ID,
	profileList,
	logint=FALSE,
	blindsub=TRUE,
	blindfold=100,
	lags=c(5,14),
	threshold=3,
	add=FALSE,
	textit=FALSE	
){

	############################################################################################
	leng<-length(components[[1]][[ID]][,1])
	getit<-unique(as.numeric(c(components[[1]][[ID]][,c(1)],components[[1]][[ID]][,c(2)])))
	maxint<-c()
	for(i in 1:length(getit)){
		int<-as.numeric(profileList[[2]][(profileList[[7]][profileList[[7]][,4]==getit[i],1]:profileList[[7]][profileList[[7]][,4]==getit[i],2]),2])
		maxint<-c(maxint,max(int))
	}
	getit<-getit[order(maxint,decreasing=TRUE)]
	plotaprofile(
		profileList,
		profileID=getit[1],
		logint,
		blindsub,
		blindfold,
		lags,
		threshold,
		add=FALSE,
		textit=FALSE
	)
	if(length(getit)>2){
		for(i in 2:length(getit)){
			plotaprofile(
				profileList,
				profileID=getit[i],
				logint=FALSE,
				blindsub=TRUE,
				blindfold,
				lags,
				threshold,
				add=TRUE,
				textit=FALSE,
				colorit="blue"
			)
		}
	}
	############################################################################################
	
}