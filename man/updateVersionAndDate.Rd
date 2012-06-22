\name{updateVersionAndDate}
\alias{updateVersionAndDate}
\title{updateVersionAndDate}
\description{update the version and the date in the DESCRIPTION file }
\usage{updateVersionAndDate(versionUpdate = c(svnHead = "svnHead", svnCur = "svnCur", 
    inc = "inc"), packageDir = ".")}
\arguments{
  \item{versionUpdate}{method of updating the version\describe{
\item{svnHead}{replace number after last dot by the subversion head revision +1}
\item{svnCur}{replace number after last dot by highest subversion revision number of the current workspace}
\item{inc}{increment the number after the last dot}
}}
  \item{packageDir}{the root directory of the package (default current dir)}
}
\details{Replaces the part after last dot with}
\value{character vector: the version and the revisions }

\author{Thomas Wutzler}



\seealso{\link{twDev}}
\examples{
#mtrace(updateVersionAndDate)
if( FALSE ){
	# take care. This will change the DESCRIPTION file, may revert afterwards
	updateVersionAndDate(versionUpdate="svnCur")
	{tmp <- scan('DESCRIPTION', what = character(), sep = '\n'); tmp[ grep("Version:",tmp)]}
	updateVersionAndDate(versionUpdate="inc")
	{tmp <- scan('DESCRIPTION', what = character(), sep = '\n'); tmp[ grep("Version:",tmp)]}
	updateVersionAndDate()
	{tmp <- scan('DESCRIPTION', what = character(), sep = '\n'); tmp[ grep("Version:",tmp)]}
}
}
