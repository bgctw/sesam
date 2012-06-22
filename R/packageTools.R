getPackageId <- function(
	### Read the package id from the DESCRIPTION file in the current working directory.
	...	##<< further arguments to \code{\link{getPackageDesc}}
){
	##seealso<< \link{twDev}
	getPackageDesc(...)$pkg
	### Scalar String: the package identifier
}

getPackageDesc <- function(
	### Parse the DESCRIPTION file
	userInstallDir = "c:/tmp"	##<< The directory of local user libraries (if the R-library directory is not writeable)
){
	##seealso<< \link{twDev}
	#
	##details<<
	## Assumes that current working directory is the top level of the package,
	## and that the DESCRIPTION file resides in this current directory
	installDir <- if( file.access(.libPaths()[1],2)==0 ) .libPaths()[1] else  
			userInstallDir          # if standard libDir is not writeable, this gives the directory where the package will be installed
	desc <- readLines("DESCRIPTION")
	pkg <- str_trim(gsub("(^Package:\\s+)(.+)$", "\\2", desc[grep("^Package:",desc)], perl=TRUE))
	##value<< list with entries
	list(	
		pkg = pkg	##<< the package id
		,pkgVersion = gsub("(^Version:\\s+)(.+)$", "\\2", str_trim(desc[grep("^Version:",desc)]), perl=TRUE)	##<< the package version
		,installDir = installDir	##<< the directory, where the package will be installed to
		,htmlRoot = file.path( installDir, pkg, "html" )  ##<< the directory where the generated html documentation pages reside.
	)
	##end<<
}
attr(getPackageDesc,"ex") <- function(){
	if( FALSE ){
		getPackageDesc()
	}
}

compileSrc <- function(
	### Compile all c files in src directory to dll
	options = ""		##<< further options to R CMD SHLIB
	,unlink=TRUE		##<< unlink dll before
	,dynFilenameLocal = file.path("src",paste(dllBaseName, .Platform$dynlib.ext, sep = "")) ##<< the name of the dll
	,dllBaseName = getPackageId()	##<< the base name of the dll, may be provided for performance reasons.
	,snowfall=FALSE		##<< if TRUE, then dll is unlinked and relinked in all snowfall subprocesses 
){
	if( unlink ){
		try(dyn.unload(dynFilenameLocal), silent=TRUE)
		if( snowfall ){
			sfExport("dynFilenameLocal")
			tmp <- try( sfClusterEval( dyn.unload(dynFilenameLocal)),  silent=TRUE)
		}
	}
	cFileNames <- Sys.glob(file.path("src","*.c"))
	#unlink("src/*.o")
	system(tmp <- paste("R CMD SHLIB -o ",dynFilenameLocal,paste(cFileNames, collapse=" ")))
	if( unlink ){
		dyn.load(dynFilenameLocal)
		if( snowfall ){
			sfClusterEval( dyn.load(dynFilenameLocal) )
		}
	}
}

cleanSrc <- function(
	### Delete all object files and dll in src directory 
	unlink=TRUE		##<< unlink dll before
	,dynFilenameLocal = file.path("src",paste(dllBaseName, .Platform$dynlib.ext, sep = "")) ##<< the name of the dll
	,dllBaseName = getPackageId()	##<< the base name of the dll, may be provided for performance reasons.
	,snowfall=FALSE		##<< if TRUE, then dll is unlinked and relinked in all snowfall subprocesses 
){
	if( unlink ){
		tmp <- try( dyn.unload(dynFilenameLocal), silent=TRUE )
		if( snowfall ){
			sfExport("dynFilenameLocal")
			tmp <- try( sfClusterEval( dyn.unload(dynFilenameLocal)),  silent=TRUE)
		}
	}
	unlink("src/*.o")
	unlink(dynFilenameLocal)
}

genRd <- function(
	### Generate the Rd files from inlinedoc comments
	genHtml=TRUE			##<< if TRUE, the package is compiled and html-Files are generated			
	,showInBrowser = genHtml	##<< if TRUE then package is installed, and generated documentation is shown in a brower window
){
	##seealso<< \link{twDev}
	#
	##details<<  
	## The man directory is cleared first. Then inlinedocs is executed. 
	## Third all Rd-files from inst/genData are copied to the man dir, overwriting the generated files.
	## The third step is useful in providing a more detailed package.Rd and description of datasets.
	pkg <- getPackageId()
	# generate RD Files
	if( require(inlinedocs) ){
		unlink( file.path("man","*.Rd") )
		#prevWd <- setwd("..")
		#tryCatch( package.skeleton.dx(pkg), finally=setwd(prevWd))
		package.skeleton.dx(".")
		try(file.copy( Sys.glob(file.path("inst","genData","*.Rd")), "man", overwrite=TRUE ), silent=TRUE)	# copy descriptions of data
		if( genHtml ){
			# generate the HTML  files
			prevWd <- setwd("..")
			system(	paste("R CMD INSTALL --html ",pkg, sep="") )
			setwd(prevWd)
			# show in Browser
		}
		if( showInBrowser ){
			htmlRoot <- file.path( system.file(package = pkg), "html" )
			html_viewer <- function(path) {
				browser <- getOption("browser")
				if(is.null(browser) && .Platform$OS.type == "windows")
					shell.exec(chartr("/", "\\", path))
				else browseURL(paste("file://", URLencode(path), sep=""))
			}
			html_viewer(file.path(htmlRoot,paste(pkg,"-package.html",sep="")) )
		}
	}else{ warning("could not load package inlinedocs")}
}
attr(genRd,"ex") <- function(){
	if( FALSE ){	# do not execute on installation
		genInlinedocsRd()
	}
}

runTests <- function(
	### running unit tests in default directory.
){
	test_dir("inst/tests")
}

runRCheck <- function(
	### running R CMD check
	options="--no-manual"	##<< further options
	,pkg = getPackageId()	##<< package name, default parse DESCRIPTION
){
	cmd <- paste("R CMD check",options,pkg)
	prevWd <- setwd("..");	tryCatch( system(cmd), finally=setwd(prevWd))
}
attr(runRCheck,"ex") <- function(){
	if( FALSE ){
		runRCheck()
	}
}

installPkg <- function(
	### running R CMD check
	options="--html"		##<< further options
	,pkg = getPackageId()	##<< package name, default parse DESCRIPTION
){
	cmd <- paste("R CMD INSTALL",options,pkg)
	prevWd <- setwd("..");	tryCatch( system(cmd), finally=setwd(prevWd))
}
attr(installPkg,"ex") <- function(){
	if( FALSE ){
		installPkg()
	}
}

svnUpdate <- function(
	### running svn update
){
	cmd <- paste("svn update")
	system(cmd)
}
attr(svnUpdate,"ex") <- function(){
	if( FALSE ){
		svnUpdate()
	}
}

svnCommit <- function(
	### running svn commit -m msg
	msg					##<< commit message
	,updateDesc=TRUE	##<< if TRUE, update date and version of the DESCRIPTION file (see \code{\link{updateVersionAndDate}})
){
	if( updateDesc )
		updateVersionAndDate()
	cmd <- paste("svn commit -m '",msg,"'")
	system(cmd)
}
attr(svnCommit,"ex") <- function(){
	if( FALSE ){
		svnCommit("Testing svnCommit")
	}
}

loadLibs <- function(
	### loading all the libraries indicated in the depends tag of the DESCRIPTION file
){
	description.content <- scan('DESCRIPTION', what = character(), sep = '\n')
	lineDep <- description.content[grep('Depends:', description.content)]
	deps0 <- unlist(lapply( str_split( sub('[Dd]epends:(.*)$',"\\1" ,   lineDep), "," ), str_trim ))
	deps <- {tmp <- grep("R ",deps0); if( length(tmp) ) deps0[-tmp] else deps0 }	# remove the dependence on R
	tmp <- sapply( deps, function(dep){ 
			#print(dep)
			require(dep, character.only = TRUE) 
		})
	### boolean vector indicatin if package was loaded successfully
}

loadDll <- function(
	### Load DLL (Shared object) from src directory
	dynFilenameLocal = file.path("src",paste(pkg, .Platform$dynlib.ext, sep = ""))
	,pkg = getPackageId()	##<< the name of the package, may be provided for performance reasons.
){
	if( file.exists(dynFilenameLocal) ){
		dyn.load(dynFilenameLocal)
	}else{
		warning(paste("dynLib",dynFilenameLocal,"not found."))
	}
}

sourceR <- function(
	### source all the files in the R subdirectory.
){
	tmp <- sapply(Sys.glob(file.path("R","*.R")), source)		# load all the R code in file in dir R
}
attr(sourceR,"ex") <- function(){
	if( FALSE ){
		svnCommit("Testing svnCommit")
	}
}

loadPkg <- function(
	### \code{\link{loadLibs}}, \code{\link{loadDll}}, and \code{\link{sourceR}}
){
	loadLibs()
	loadDll()
	sourceR()
}

reloadPkg <- function(
	### Reloading a package
	pkg		##<< name of the package (scalar character string)
	,force=TRUE	##<< logical: should a package be detached even though other attached packages depend on it?
	,...	##<< other arguments to \code{\link{detach}}
){
	detach(paste0("package:",pkg), unload="TRUE", character.only = TRUE, force=force,...)
	library(pkg, character.only = TRUE)
}

buildPkg <- function(
	### create the zip file of the package, that can be installed from R-GUI.
	pkg = getPackageId()
){
	##details<< 
	## Run from windows it will be a zip file 
	## while run from a nix system it will be a tgz file relating to the architecture.
	prevWd <- setwd("..")
	tryCatch(
		system(	paste("R CMD INSTALL --build ",pkg, sep="") )
		,finally=setwd(prevWd)
	)
}

put2mdiCodeDir <- function(
	### Copy the source-Code, package zip to mdi-code directory.
	mdiBaseDir="/Net/Groups/C-Side/MDI"		##<< The path to the mdi-root directory
	,mdiCodeDir = file.path(mdiBaseDir,"_code/R") ##<< the relative path to the code directory
	,doRebuild=TRUE		##<< if TRUE then \code{\link{buildPkg}} is invoked
	,pInfo=getPackageDesc()	##<< overwrite when requiring different arguments to \code{\link{getPackageDesc}}
){
	##details<< 
	## When done with local development, archive the code on the MDI code directory
	## Take care - the existing version in the MDI code directory will be overwritten.
	pkg <- pInfo$pkg
	pkgVersion <- pInfo$pkgVersion
	#
	destDir <- file.path(mdiCodeDir,pkg)  	# might need to change mounted path e.g. m:
	unlink(destDir, recursive=TRUE) 
	try(dir.create(destDir, mode = "0755"))		# mode ignored on windows, need to set by hand, best on unix eg. dialog chmod 775 <dir>
	file.copy( file.path("..", pkg), mdiCodeDir, overwrite=TRUE, recursive=TRUE)
	#
	# rebuild the zip file 
	if( doRebuild ) buildPkg(pkg=pkg)
	# 	
	#  copy also the buildt zip and/or tgz files generated above
	file.copy( Sys.glob( file.path("..",paste(pkg,'_',pkgVersion,'*',sep=""))), mdiCodeDir )
}
#trace(put2mdiCodeDir,recover)	#untrace(put2mdiCodeDir)

put2MdiRep <- function(
	### Publish the generated html package description on MDI internal pages.
	mdiBaseDir="/Net/Groups/C-Side/MDI"		##<< The path to the mdi-root directory
	,mdiCodeHtmlDir = file.path(mdiBaseDir,"public_html/code_doc") ##<< the relative path to the code directory
	,pInfo=getPackageDesc()	##<< overwrite when requiring different arguments to \code{\link{getPackageDesc}}
	,doGenRd=TRUE			##<< if TRUE then genRd() is invoked, which generates the html pages.
	,doPutCode=TRUE		##<< if TRUE then \code{\link{put2mdiCodeDir}} is called 
	,...				##<< further arguments to \code{\link{put2mdiCodeDir}}, such as userInstallDir 
){
	##details<<
	## Access as http://www.bgc-jena.mpg.de/bgc-mdi/code_doc/<pkg>/html/00Index.html
	## 
	## Assumes that the package has been installed. This is done if \code{genRd=TRUE}.
	##
	## Best run from R inside a linux shell:
	## When run from windows, must still set the permissions in a unix shell: chmod 775 <dir>.
	## Sometimes Windows also has problems in accessing the html-dir.
	#
	pkg <- pInfo$pkg
	if( doPutCode ) put2mdiCodeDir(mdiBaseDir=mdiBaseDir, pInfo=pInfo, ...)
	if( doGenRd ) genRd(showInBrowser=FALSE)
	#
	#Sys.glob(file.path(htmlRoot,"*"))
	destDir <- file.path(mdiCodeHtmlDir,pkg)		# might need to change mounted path e.g. m:
	#unlink(destDir, recursive=TRUE) 
	unlink(Sys.glob(file.path(destDir,"*")), recursive=TRUE)
	try(dir.create(destDir, mode = "0775"))		# mode ignored on windows, need to set by hand, best on unix eg. dialog chmod 775 <dir>
	#
	# copy the hmtl description generated above
	pkgRoot <- system.file(package = pkg)
	htmlRoot <- file.path( pkgRoot, "html" )
	file.copy( htmlRoot, destDir, recursive=TRUE)
	#Sys.glob( file.path(htmlRoot,"*") )
	file.copy( file.path(pkgRoot,"DESCRIPTION"), destDir)
	try( system( paste("chmod 755 ",destDir,"/DESCRIPTION",sep="")), silent=TRUE) 
	# copy the pdf description generated above
	#file.copy( file.path("..",paste(pkg,".pdf",sep="")), destDir )
	# copy the buildt zip and/or tgz files generated above
	#file.copy( Sys.glob( file.path("..",paste(pkg,'_',pInfo$pkgVersion,'*',sep=""))), destDir )
}
#trace(put2MdiRep,recover) 	#untrace(put2MdiRep)






