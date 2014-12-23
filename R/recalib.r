#' @title Function for m/z and RT recalibration.
#'
#' @export
#'
#' @description Function for m/z and RT recalibration.
#'
#' @param peaklist Matrix or data.frame of sample peaks. 3 columns: m/z, intensity, RT.
#' @param mz Numeric vector. Masses to be matched within large tolerances of \code{tolmz}.
#' @param tolmz Numeric. +/- m/z tolerance (precision)
#' @param ppm Logical. \code{tolmz} given in ppm?
#' @param ret Numeric vector. Retention times to be matched; units equal to those of the input files
#' @param tolret Numeric. RT tolerance for matches with \code{peaklist}; units equal to those of the input files
#' @param what Character strings "mass" or "ret".
#' @param one Logical. Only use recalibration \code{mz, ret} that can be matched unambiguously.
#' @param knot Integer. Number of spline knots.
#' @param plotit Logical. Produce recalibration plot?
#' @param path_1 Logical \code{FALSE} or character string. If not \code{FALSE}, filepath to output the plot
#' @param path_2 Logical \code{FALSE} or character string. Filepath for saving GAM model.
#' @param stopit Logical. Triggers a full R error (\code{TRUE}) or returns an error message string (\code{FALSE}) if failing. 
#'
#' @return Recalibrated \code{peaklist}.
#' 
#' @details enviMass workflow function. \code{mz} and \code{ret} must be of same length and specify the points to recalibrate with.
#' A minimum of 15 matches must be found to conduct the recalibration.
#'


recalib<-function(
  peaklist,
  mz,
  tolmz=10,
  ppm=TRUE,
  ret,
  tolret,
  what="mass",
  one=TRUE,
  knot=5,
  plotit=FALSE,
  path_1=FALSE,
  path_2=FALSE,  
  stopit=FALSE
  ){


  ##############################################################################
  # search for concerned peaks #################################################
  if(what!="mass" && what!="ret"){stop("what what?")};
  peaks<-search_peak(peaklist,mz,tolmz,ppm,ret,tolret);
  # collect concurring peaks ###################################################
  getit1<-c();  # expected mz
  getit2<-c();  # observed peaklist
  for(i in 1:length(peaks)){
    if(peaks[i]!="FALSE"){
      put<-as.numeric(strsplit(peaks[i],"/")[[1]]);
      if(what=="mass"){
        if(one==TRUE & length(put)==1){
          getit1<-c(getit1,rep(mz[i],length(put)));  # expected
          getit2<-c(getit2,peaklist[put,1]);         # observed
        }else{
          getit1<-c(getit1,rep(mz[i],length(put)));  # expected
          getit2<-c(getit2,peaklist[put,1]);         # observed
        }
      }else{ # what=="ret"
        if(one==TRUE & length(put)==1){
          getit1<-c(getit1,rep(ret[i],length(put))); # expected
          getit2<-c(getit2,peaklist[put,3]);         # observed
        }else{
          getit1<-c(getit1,rep(ret[i],length(put))); # expected
          getit2<-c(getit2,peaklist[put,3]);         # observed
        }
      }
    }
  }
  getit3<-c(getit1-getit2); # observed-expected
  ##############################################################################
  # stop if too few data found! ################################################
  if(length(getit3)<15){
		if(path_1!="FALSE"){png(filename = path_1, bg = "white")}    
		plot.new();
		plot.window(xlim=c(1,1),ylim=c(1,1));
		box();
		text(1,1,label="not available",cex=1.5,col="darkred") 
		if(path_1!="FALSE"){dev.off()}   
		if(!stopit){
			return("Too few data points for fit!\n");
		}else{
			stop("Too few data points for fit!\n")
		}
  }
  ##############################################################################
  # train gam model ############################################################
  that<-data.frame(getit2,getit3)
  names(that)<-c("obs","delta")
  #attach(that)
  model<-mgcv::gam(delta~s(obs,bs="ts",k=knot),data=that);
  if(plotit==TRUE){
    if(what=="mass"){    
      if(path_1!="FALSE"){png(filename = path_1, bg = "white")}
      plot(getit2,getit3,pch=19,cex=0.5,xlab="m/z",ylab="Expected m/z - observed m/z",main="Recalibration results");
      abline(h=0,col="red");
      points(getit2[order(getit2)],predict(model)[order(getit2)],col="red",type="l",lwd=2);
      #points(getit2,getit3-predict(model),col="green",pch=19,cex=0.7);      
      if(path_1!="FALSE"){dev.off()}   
    }else{
      if(path_1!="FALSE"){png(filename = path_1, bg = "white")} 
      plot(getit2,getit3,pch=19,cex=0.5,xlab="Retention time",ylab="Expected RT - observed RT",main="Recalibration results");
      abline(h=0,col="red");
      points(getit2[order(getit2)],predict(model)[order(getit2)],col="red",type="l",lwd=2);
      #points(getit2,getit3-predict(model),col="green",pch=19,cex=0.7);
      if(path_1!="FALSE"){dev.off()}      
    }
  }
  ##############################################################################
  # predict -> recalibrate peaklist ############################################
  if(what=="mass"){
      that<-data.frame("obs"=peaklist[,1],"delta"=peaklist[,1]);
      pred2<-mgcv::predict.gam(model,newdata=that);
      newpeaks<-peaklist;
      newpeaks[,1]<-c(peaklist[,1]+pred2);
  }else{
      that<-data.frame("obs"=peaklist[,3],"delta"=peaklist[,3]);
      pred2<-mgcv::predict.gam(model,newdata=that);
      newpeaks<-peaklist;
      newpeaks[,3]<-c(peaklist[,3]+pred2);
  }
  if(path_2!="FALSE"){
	save(model,file=path_2)
  }
  return(newpeaks)
  ##############################################################################
}
