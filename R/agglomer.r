#' @title Partition the peak lists pooled by \code{\link{startprofiles}}. 
#'
#' @export
#'
#' @description \code{agglomer} partitions all peaks from various files into unassociated regions.
#'
#' @param profileList A profile list.
#' @param dmass Numeric. m/z gap size
#' @param ppm Logical. \code{dmass} given in ppm?
#' @param dret Numeric. RT gap size; units equal to those of the input files
#'
#' @return Updated profile list
#' 
#' @details enviMass workflow function
#' 
#' @seealso  \code{startprofiles}, \code{partcluster} 

agglomer<-function(
	profileList,
	dmass=2,
	ppm=TRUE,
	dret=60
){

  ##############################################################################
  if(!is.loaded("agglom")){stop(".Call to agglom failed; aborted.")}
  if(!is.loaded("indexed")){stop(".Call to indexed failed; aborted.")}
  if(!is.numeric(dmass)){stop("dmass must be numeric; aborted.")}
  if(!is.numeric(dret)){stop("dret must be numeric; aborted.")}
  if(!is.logical(ppm)){stop("ppm must be logical; aborted.")}
  ##############################################################################
  if(ppm){ppm2<-1}else{ppm2<-0}
  part <- .Call("agglom",
    as.numeric(profileList[[2]][,1]),
    as.numeric(profileList[[2]][,3]),
    as.integer(profileList[[2]][,6]),
    as.integer(ppm2),
    as.numeric((dmass*2)+1),
    as.numeric(dret),
    PACKAGE="enviMass" 
  )
  profileList[[2]][,7]<-part;
  profileList[[2]]<-profileList[[2]][order(profileList[[2]][,7],decreasing=FALSE),]
  maxit<-max(profileList[[2]][,7])
  index<- .Call("indexed",
    as.integer(profileList[[2]][,7]),
    as.integer(maxit),
    as.integer(3),
    PACKAGE="enviMass" 
  )
  index<-index[index[,1]!=0,]
  colnames(index)<-c("start_ID","end_ID","number_peaks")
  profileList[[6]]<-index
  profileList[[1]][[2]]<-TRUE;
  ##############################################################################
  return(profileList)
  
}


