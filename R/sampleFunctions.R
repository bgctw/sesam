sampleFunc <- function(
	### Adding two vectors
	a.n.v		##<< numeric vector: first term
	,b.n.v		##<< numeric vector: second term
){
	##seealso<< \code{\link{myPackageId}}
	#
	##details<< 
	## If the two vectors are of different length, the shorter one
	## is recycled.
	#
	a.n.v + b.n.v 
	### The componentwise sum of the two arguments
}
attr(sampleFunc,"ex") <- function(){
	# This example appears in the generated help.
	sampleFunc(1:10, 2)
}
