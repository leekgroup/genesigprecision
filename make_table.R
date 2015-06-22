
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
	final[,"rot"] <- final[,"rot"] * 100
	final[,"col"] <- final[,"col"] * 100
	rownames(final) <- c("$W_{-ER}$", "$W_C$", "$W_G$", "$W_{CG}$")
	xt <- xtable(final, digits=c(0,4,5,4,5,4,5,2,2))
	print(xt, sanitize.rownames.function = identity)
}
