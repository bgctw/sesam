.getTopSvnRevision <- function(
	### use svnversion to get the lasted used revision number
){
	revisions  <- system('svnversion -n .', intern =TRUE)
	if( 0 != length(grep(":", revisions)) ){
		warning(".getTopSvnRevision: mixed revisions, commit first")
	}
	topRevision <- sub("^[^:]*:","",revisions) # only keep last of mixed revisions
	if( 0 != length(grep("M", topRevision)) ){
		warning(".getTopSvnRevision: modified revision, maybe only version number changed, maybe need to commit first")
	}
	topRevision
}
attr(.getTopSvnRevision,"ex") <- function(){
	if( FALSE ){
		.getTopSvnRevision()	# svn need to be installed 			
	}
}

.getHeadSvnRevision <- function(
### use svnversion to get the lasted used revision number
){
	revisions  <- system('svn info --revision HEAD', intern =TRUE)
	sub("Revision: ?(.*)$","\\1",revisions[grep("Revision:",revisions)])
}
attr(.getHeadSvnRevision,"ex") <- function(){
	if( FALSE ){
		.getHeadSvnRevision()	# svn need to be installed 			
	}
}

updateVersionAndDate <- function(
	### update the version and the date in the DESCRIPTION file 
	versionUpdate=c(		##<< method of updating the version
		##describe<<
		svnHead="svnHead"	##<< replace number after last dot by the subversion head revision +1 
		,svnCur="svnCur"	##<< replace number after last dot by highest subversion revision number of the current workspace
		,inc="inc"			##<< increment the number after the last dot
		) 
		##end<<
	,packageDir="."			##<< the root directory of the package (default current dir)
){
	##seealso<< \link{twDev}
	
	##details<<
	## Replaces the part after last dot with
	prevWd <- setwd(packageDir)
	tryCatch({
		description.content <- scan('DESCRIPTION', what = character(), sep = '\n')
		version.line <- grep('Version:', description.content)
		#
		if( version.line ){
			versionUpdate <- match.arg(versionUpdate) 
			version.number <- switch( versionUpdate,
				"svnCur" = as.integer(sub("[[:alpha:]]","",.getTopSvnRevision() )) # remove M and S
				,"svnHead" = 1+as.integer(.getHeadSvnRevision()) 
				,"inc" = 1+as.integer(sub('.*\\.([^.]*$)',"\\1" ,   description.content[version.line]))
				, stop(paste("unknown versionUpdate",versionUpdate))
			)
			if( is.finite(version.number) ){
				newLine <- sub('\\.[^.]*$', paste('.',version.number,sep=""),   description.content[version.line])
				description.content[version.line] <- newLine
			}
		} else{
			warning("no Version: line found in DESCRIPTION file.")
		}
		#
		date.line <- grep('Date:', description.content)
		if( date.line ){
			newLineDate <- sub("^Date:.*$", paste("Date: ",format(Sys.time(),"%Y-%m-%d"),sep=""),   paste(description.content[date.line]))
			description.content[date.line] <- newLineDate
		}
		write(description.content, file = 'DESCRIPTION', sep = '\n')
	}, finally=setwd(prevWd) )
	version.number
	### character vector: the version and the revisions 
}
attr(updateVersionAndDate,"ex") <- function(){
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

