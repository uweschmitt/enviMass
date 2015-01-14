#' @title Plot intensities histogram of a profile
#'
#' @export
#'
#' @description \code{plot_profile_intensities_histogram} plots the intensities histogram of given profiles

#'
#' @param mean_intensities vector of mean intensities in profiles
#' @param max_intensities vector of max intensities in profiles
#' @param past_incidentss vector of global trands in profiles
#' @param current_incidentss vector of current trends in profiles
#' 
#' @details enviMass workflow plot function
#' 



plot_profiles_intensity_histograms <- function(mean_intensities, max_intensities, past_incidents,
                                               current_incidents) {

    .gtz <- function(what)
    {
        return(what[what > 0]);
    }

    mean_intensities <- .gtz(mean_intensities);
    max_intensities <- .gtz(max_intensities);
    past_incidents <- .gtz(past_incidents);
    current_incidents <- .gtz(current_incidents);

    if(length(mean_intensities)>0){    # MEAN INTENSIties
        aa<-c(hist(log10(mean_intensities), breaks=100, plot=FALSE), 0)
        aa<-max(aa$counts)
    }else{
        aa<-(0)
    }
    if(length(max_intensities)>0){    # MAX INTENSIties
        bb<-c(hist(log10(max_intensities), breaks=100, plot=FALSE), 0)
        bb<-max(bb$counts)
    }else{
        bb<-(0)
    }
    if(length(past_incidents)>0){  # PAST incidents
        cc<-c(hist(log10(past_incidents), breaks=100, plot=FALSE), 0)
        cc<-max(cc$counts)
    }else{
        cc<-(0)
    }
    if(length(current_incidents)>0){   # CURRENT incidents
        dd<-c(hist(log10(current_incidents), breaks=100, plot=FALSE), 0)
        dd<-max(dd$counts)
    }else{
        dd<-(0)
    }
    aaa<-max(aa, bb, cc, dd)
    if(aaa>0){
        if (aa != 0) {
            hist(log10(mean_intensities), breaks=100, xlab="log10 Intensities", border="darkgrey",
                 col=NULL, ylim=c(0, aaa), main="Intensities distributions of profiles")
        }
        if (bb != 0) {
            hist(log10(max_intensities), breaks=100, add=TRUE, border="darkgreen", col=NULL)
        }
        if(cc != 0) {
            hist(log10(past_incidents), breaks=100, add=TRUE, border="red", col=NULL)
        }
        if(dd != 0) {
            hist(log10(current_incidents), breaks=100, add=TRUE, border="blue", col=NULL)
            rug(log10(current_incidents), col="blue", quiet=TRUE)
        }
        plot.window(xlim=c(0, 10), ylim=c(0, 10))
        text(7, 9, labels="Mean intensities", col="darkgrey", pos=4)
        text(7, 8.5, labels="Maximum intensities", col="darkgreen", pos=4)
        text(7, 8, labels="Global trend intensities", col="red", pos=4)
        text(7, 7.5, labels="Current trend intensities", col="blue", pos=4)
    }else{
        plot.new()
        plot.window(xlim=c(0, 1), ylim=c(0, 1))
        text(0.5, 0.5, labels="histogram infeasible", cex=1.8, col="red")
    }
}

