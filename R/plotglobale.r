#' @title Plot all time series from a list of profiles.
#'
#' @export
#'
#' @description \code{plotglobal} plots all or a subset of most intense time series contained in \code{plofileList}
#'
#' @param profileList A profile list.
#' @param n_global Numeric. How many global most intense profiles to include?
#' @param n_present Numeric. How many present most intense profiles to include?
#' @param use_present Logical. 
#' @param highlight Logical.
#' 
#' @details enviMass workflow function
#' 


plotglobal<-function(
	profileList,
	n_global=10,
	n_present=10,
	use_present=FALSE,
	highlight=FALSE
){

	##################################################################################
	getit<-seq(1,length(profileList[[7]][,1]),1)
	getit_global<-getit[profileList[[7]][,6]!=0]
	getit_global<-getit_global[order(profileList[[7]][getit_global,6],decreasing=TRUE)]	
	getit_global<-getit_global[1:n_global]
	if(use_present=="FALSE"){
		getit_present<-getit[profileList[[7]][,5]!=0]
		getit_present<-getit_present[order(profileList[[7]][getit_present,5],decreasing=TRUE)]
		getit_present<-getit_present[1:n_present]
	}else{
		getit_present<-use_present
	}
	##################################################################################
	# global #########################################################################
	plotaprofile(
		profileList,
		getit_global[1],
		logint=FALSE,
		blindsub=FALSE,
		blindfold=100,
		lags=c(5,14),
		threshold=3,
		add=FALSE,
		textit=FALSE,
		simple=TRUE,
		colorit="darkgrey",
		use_lwd=1
	)
	if(n_global>1){
	for(i in 2:n_global){
		plotaprofile(
			profileList,
			getit_global[i],
			logint=FALSE,
			blindsub=FALSE,
			blindfold=100,
			lags=c(5,14),
			threshold=3,
			add=TRUE,
			textit=FALSE,
			simple=TRUE,
			colorit="darkgrey",
			use_lwd=1
		)
	}
	}
	##################################################################################
	# highlight ######################################################################
	if(highlight[1]!=FALSE){
	for(i in 1:length(highlight)){
		plotaprofile(
			profileList,
			highlight[i],
			logint=FALSE,
			blindsub=FALSE,
			blindfold=100,
			lags=c(5,14),
			threshold=3,
			add=TRUE,
			textit=FALSE,
			simple=TRUE,
			colorit="red",
			use_lwd=2
		)
	}
	}
	##################################################################################
	# present ########################################################################
	if(highlight[1]==FALSE){
	for(i in 1:n_present){
		plotaprofile(
			profileList,
			getit_present[i],
			logint=FALSE,
			blindsub=FALSE,
			blindfold=100,
			lags=c(5,14),
			threshold=3,
			add=TRUE,
			textit=FALSE,
			simple=TRUE,
			colorit="red",
			use_lwd=1
		)
	}
	}
	##################################################################################

}