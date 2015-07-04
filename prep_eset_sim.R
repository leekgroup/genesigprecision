# This function prepares ExpressionSet data for simulation run
# It relabels the relevant data columns and makes MammaPrint predictions using the genefu package

prep_eset_sim <- function(dat){
	library(genefu)
	# Need annotation to make MammaPrint predictions via genefu package
	colnames(fData(dat)) <- c("EntrezGene.ID", "Symbol")
	pd <- pData(dat)
	pd$mammaprint <- gene70(t(exprs(dat)), fData(dat), do.mapping=T)$risk
	# Indicators for grade
	pd$g2ind <- ifelse(pd$grade == 2, 1, 0)
	pd$g3ind <- ifelse(pd$grade == 3, 1, 0)
	pd$y <- pd$surv5
	pd
}
