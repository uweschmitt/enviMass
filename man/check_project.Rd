% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/check_project.r
\name{check_project}
\alias{check_project}
\title{Check enviMass project}
\usage{
check_project(logfile, isotopes, adducts)
}
\arguments{
\item{logfile}{enviMass project logfile}

\item{isotopes}{Isotope list}

\item{adducts}{Adducts list from package enviPat}
}
\value{
Character string with either specific error message or a "Project consistent" message.
}
\description{
\code{} runs several cconsistency checks on a enviMass project
}
\details{
enviMass workflow function; run before further calculations are started in the workflow.
}

