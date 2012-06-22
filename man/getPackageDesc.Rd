\name{getPackageDesc}
\alias{getPackageDesc}
\title{getPackageDesc}
\description{Parse the DESCRIPTION file}
\usage{getPackageDesc(userInstallDir = "c:/tmp")}
\arguments{
  \item{userInstallDir}{The directory of local user libraries (if the R-library directory is not writeable)}
}
\details{Assumes that current working directory is the top level of the package,
and that the DESCRIPTION file resides in this current directory}
\value{list with entries
\item{pkg}{the package id}
\item{pkgVersion}{the package version}
\item{installDir}{the directory, where the package will be installed to}
\item{htmlRoot}{the directory where the generated html documentation pages reside.}}

\author{Thomas Wutzler}



\seealso{\link{twDev}}
\examples{
if( FALSE ){
	getPackageDesc()
}
}
