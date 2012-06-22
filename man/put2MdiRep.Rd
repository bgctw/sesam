\name{put2MdiRep}
\alias{put2MdiRep}
\title{put2MdiRep}
\description{publish the generated html package description on MDI internal pages.}
\usage{put2MdiRep(mdiBaseDir = "/Net/Groups/C-Side/MDI", mdiCodeHtmlDir = file.path(mdiBaseDir, 
    "public_html/code_doc"), pInfo = getPackageDesc(), doGenRd = TRUE, 
    doPutCode = TRUE, ...)}
\arguments{
  \item{mdiBaseDir}{The path to the mdi-root directory}
  \item{mdiCodeHtmlDir}{the relative path to the code directory}
  \item{pInfo}{overwrite when requiring different arguments to \code{\link{getPackageDesc}}}
  \item{doGenRd}{if TRUE then genRd() is invoked, which generates the html pages.}
  \item{doPutCode}{if TRUE then \code{\link{put2mdiCodeDir}} is called}
  \item{\dots}{further arguments to \code{\link{put2mdiCodeDir}}, such as userInstallDir}
}
\details{Access as http://www.bgc-jena.mpg.de/bgc-mdi/code_doc/<pkg>/html/00Index.html

Assumes that the package has been installed. This is done if \code{genRd=TRUE}.

Best run from R inside a linux shell:
When run from windows, must still set the permissions in a unix shell: chmod 775 <dir>.
Sometimes Windows also has problems in accessing the html-dir.}


\author{Thomas Wutzler}





