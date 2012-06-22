\name{put2mdiCodeDir}
\alias{put2mdiCodeDir}
\title{put2mdiCodeDir}
\description{Copy the source-Code, package zip to mdi-code directory.}
\usage{put2mdiCodeDir(mdiBaseDir = "/Net/Groups/C-Side/MDI", mdiCodeDir = file.path(mdiBaseDir, 
    "_code/R"), doRebuild = TRUE, pInfo = getPackageDesc())}
\arguments{
  \item{mdiBaseDir}{The path to the mdi-root directory}
  \item{mdiCodeDir}{the relative path to the code directory}
  \item{doRebuild}{if TRUE then \code{\link{buildPkg}} is invoked}
  \item{pInfo}{overwrite when requiring different arguments to \code{\link{getPackageDesc}}}
}
\details{When done with local development, archive the code on the MDI code directory
Take care - the existing version in the MDI code directory will be overwritten.}


\author{Thomas Wutzler}





