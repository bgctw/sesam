.tmp.f <- function(
### dummy, holds code that is not to be executed
){
	install.packages("R.utils") #GString
	install.packages("R.oo")
	install.packages("debug")
	install.packages("deSolve")
	install.packages("nlme")
	install.packages("snowfall") #parallel
	install.packages("sensitivity")
	install.packages("mnormt")
	install.packages("mvtnorm")
	install.packages("abind")
	install.packages("coda")
	install.packages("abind")
	install.packages("inlinedocs",repos="http://r-forge.r-project.org")
	install.packages("RUnit")
	install.packages("ggplot2")
	library(inlinedocs)			#package.skeleton.dx
	library(debug)
}

twWin <- function(
	### Creates a new window with default size and par.
	width=4.6,height=3.2,pointsize=10, record=TRUE
	,xpos=40,ypos=-80
	,...
){
	##seealso<< \link{twDev}, \link{dev.set}
	cdev = dev.cur()
	windows(width=width,height=height,pointsize=pointsize, record=record, xpos=xpos, ypos=ypos, ...)	# nice to get 4 plot in powerpoint screen
	par( las=1 )					#also y axis labels horizontal
	par(mar=c(2.0,3.3,0,0)+0.3 )  #margins
	par(tck=0.02 )				#axe-tick length inside plots             
	par(mgp=c(1.1,0.2,0) )  #positioning of axis title, axis labels, axis
	#par( mfrow=c(1,2) )	#tow plots in window
	cdev
	### the device that was active before opening the new window
}

str3 <- function(
	### Compactly Display the Structure with 3 nesting levels and without attributes
	object	##<< the object to display
	,...	##<< further arguments to str		
	##seealso<< \code{\link{copy2clip}}, \link{twDev}
){
	str( object, max.level=3,  give.attr=FALSE, ...)
	### result of \code{\link{str}( object, max.level=3,  give.attr=FALSE, ...)}
}

copy2clip <- function(
	### copies argument to the clipboard
	x, col.names=NA, row.names=FALSE, quote=FALSE, sep="\t", ... 
){	
	##seealso<< \link{twDev}
	
	
	##<<seealso \code{\link{write.table}}
	write.table(x,"clipboard",sep=sep, row.names=row.names, col.names=row.names, quote=quote, ...)
	x
	### x again
}  



twDf2wikiTable <- function(
	### Mediawiki-Code for table of given data.frame. 
	ds		##<< data.frame
	,tableProps='style="float: right;"'	##<< additional text inlcuded in the header
	,copy2clip=TRUE		##<< if TRUE then attempts to copy to the (Windows)clipboard
){
	##seealso<< \code{\link{copy2clip}}, \link{twDev}
	
	rows <- paste(by(ds, 1:nrow(ds), function(row){ paste("|-\n|",paste(row[1,],collapse=" || ")) }),collapse="\n")
	heading <-  paste("!",paste(colnames(ds),collapse=" !! "))
	tableProps <- 'style="float: right;"'
	#library(R.utils)
	#tmp <- as.character(GString('{| class="wikitable" frame="hsides" ${tableProps} \n${heading}\n${rows}\n|}'))
	tmp <- paste('{| class="wikitable" frame="hsides" ',tableProps,' \n',heading,'\n',rows,'\n|}',sep="")
	if( copy2clip) copy2clip(tmp)
	tmp
	### String, Side-effect: copied to clipboard
}





