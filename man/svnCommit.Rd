\name{svnCommit}
\alias{svnCommit}
\title{svnCommit}
\description{running svn commit -m msg}
\usage{svnCommit(msg, updateDesc = TRUE)}
\arguments{
  \item{msg}{commit message}
  \item{updateDesc}{if TRUE, update date and version of the DESCRIPTION file (see \code{\link{updateVersionAndDate}})}
}



\author{Thomas Wutzler}




\examples{
if( FALSE ){
	svnCommit("Testing svnCommit")
}
}
