#' @title enviMass Workflow GUI
#'
#' @export
#'
#' @description \code{webMass()} starts the browser user interface.
#' 
#' @details Requires a browser other than Internet Explorer. Has only been tested on Firefox and Google Chrome; 
#' make sure to set the relevant browser as your default browser.
#' 



webMass<-function(){
  shiny::runApp(system.file('webMass', package='enviMass'))
}

