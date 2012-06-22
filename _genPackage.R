#------------ generating package tw.DEMC
# inlcuding generation of Rd files from inline-docs
# install.packages("inlinedocs")
pkg<-"twDev"

.tmp.loadPackages <- function(){
	# loading libraries, sourcing the code and loading data
	# usually done on startup library(MyPackage)
	
	# this code uses several packages that need to be installed once
	# install.packages(c("RUnit","inlinedocs","R.methodsS3", "abind","twMisc"), repos = c("http://R-Forge.R-project.org","http://cran.rakanu.com/"), dep = TRUE)
	# in case that you use not the current R-version, you need to download the sources tarball from https://r-forge.r-project.org/R/?group_id=887
	#   upack the sources, and issue from a shell from the folder above extracted folder twMisc "R CMD INSTALL twMisc"
	
	library(RUnit)
	library(abind)
	library(R.methodsS3)
	library(debug)
	
	tmp <- sapply(Sys.glob(file.path("R","*.R")), source)
	#data( list=twStripFileExt(basename(Sys.glob(file.path("data","*.RData")))))
}

.tmp.installInlineDocs <- function(){
	prevWd <- setwd("..")
	system(	paste("R CMD INSTALL --html inlinedocs/pkg/inlinedocs", sep="") )
	setwd(prevWd)
}
.tmp.inlinedocs <- function(){
	# generate documentation
	
	# generate RD Files
	pkg<-"twDev"
	library(inlinedocs)
	unlink( file.path("man","*.Rd") )
	#prevWd <- setwd("..")
	#tryCatch( package.skeleton.dx(pkg), finally=setwd(prevWd))
	package.skeleton.dx(".")
	try(file.copy( Sys.glob(file.path("inst","genData","*.Rd")), "man", overwrite=TRUE ), silent=TRUE)	# copy descriptions of data
	
	# generate the HTML  files
	prevWd <- setwd("..")
	system(	paste("R CMD INSTALL --html ",pkg, sep="") )
	setwd(prevWd)
	
	
	# show in Browser
	htmlRoot <- file.path( system.file(package = pkg), "html" )
	html_viewer <- function(path) {
		browser <- getOption("browser")
		if(is.null(browser) && .Platform$OS.type == "windows")
			shell.exec(chartr("/", "\\", path))
		else browseURL(paste("file://", URLencode(path), sep=""))
	}
	html_viewer(file.path(htmlRoot,paste0(pkg,"-package.html")) )
	
	# copy to the generated html into working directory
	#file.copy( htmlRoot, ".", recursive=TRUE)
	
	updateVersionAndDate()
	
}

.tmp.UnitTests <- function(){
	# executing Unit Tests before actual installation
	test_dir("inst/tests")
	
	# let R check package consistency
	copy2clip( cmd <- paste("R CMD check --no-manual ",pkg, sep="") )
	prevWd <- setwd("..");	tryCatch( system(cmd), finally=setwd(prevWd))
}





