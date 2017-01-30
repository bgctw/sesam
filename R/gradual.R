getGradualWeight <- function(
		### get a gradual weigth with near 0 for x < bound, near 1 for x > bound, logistic in between
		x
		,bound=0.1
){
	plogis(x, scale=bound/5.3)	#5.3 the 0.5percent quantile
}
attr(getGradualWeight,"ex") <- function(){
	x <- seq(-0.11,0.11,length.out=61)
	w <- getGradualWeight(x)
	plot( w ~ x)
}
