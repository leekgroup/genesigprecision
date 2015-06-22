# This is a convenience function to take a set of simulation results from
# a BatchJobs run and output them as a table with bias, variance, and percent gain

make_table <- function(y){
	outmat <- do.call(rbind, y)
	means <- matrix(colMeans(outmat), 4, 3, byrow=T)
	vars <- apply(outmat, 2, var)
	col <- rot <- vector("numeric", 4)
	for(i in 1:4){
		idx <- i + 2*(i-1)
		rot[i] <- (vars[idx] - vars[idx + 1])/vars[idx]
		col[i] <- (vars[idx] - vars[idx + 2])/vars[idx]
	}

	vars <- matrix(vars, 4, 3, byrow=T)

	final <- cbind(means[,1], vars[,1], means[,2], vars[,2], means[,3], vars[,3], rot, col)

	xtable(final, digits=5)
}
