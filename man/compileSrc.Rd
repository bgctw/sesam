\name{compileSrc}
\alias{compileSrc}
\title{compileSrc}
\description{Compile all c files in src directory to dll}
\usage{compileSrc(options = "", unlink = TRUE, dynFilenameLocal = file.path("src", 
    paste(dllBaseName, .Platform$dynlib.ext, sep = "")), dllBaseName = getPackageId(), 
    snowfall = FALSE)}
\arguments{
  \item{options}{further options to R CMD SHLIB}
  \item{unlink}{unlink dll before}
  \item{dynFilenameLocal}{the name of the dll}
  \item{dllBaseName}{the base name of the dll, may be provided for performance reasons.}
  \item{snowfall}{if TRUE, then dll is unlinked and relinked in all snowfall subprocesses}
}



\author{Thomas Wutzler}





