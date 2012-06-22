\name{genRd}
\alias{genRd}
\title{genRd}
\description{Generate the Rd files from inlinedoc comments}
\usage{genRd(genHtml = TRUE, showInBrowser = genHtml)}
\arguments{
  \item{genHtml}{if TRUE, the package is compiled and html-Files are generated}
  \item{showInBrowser}{if TRUE then package is installed, and generated documentation is shown in a brower window}
}
\details{The man directory is cleared first. Then inlinedocs is executed. 
Third all Rd-files from inst/genData are copied to the man dir, overwriting the generated files.
The third step is useful in providing a more detailed package.Rd and description of datasets.}


\author{Thomas Wutzler}



\seealso{\link{twDev}}
\examples{
if( FALSE ){	# do not execute on installation
	genInlinedocsRd()
}
}
