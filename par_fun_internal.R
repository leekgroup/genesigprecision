#'
#' Parallelized Adjusted Estimator Simulation Function
#'
#' This function constructs a simulated dataset from the data saves in genesigprecision_data.Rda
#' and computes the adjusted and unadjusted treatment effect estimates. This is done by splitting the
#' original dataset into training and testing sets, building prediction models on the training data,
#' and applying them to resampled sets of the test data (nrun resamples).
#'
#' @param seed An integer random seed
#' @param rand (logical) T indicates that you would like to permute all labels and remove corellation in the dataset. F indicates that you would like to permute records only and retain relationships between covariates
#' @param nrun Number of test set itreations run
#'
#' @export
#'
#' @return A 100 x 15 matrix with results for W_{-age}, W_C, W_G, W_C;W_G, W_CG (described in paper) over 100 simulations

par_fun_internal <- function(input, rand=F){

	#set.seed(seed)

	# If rand==F, we are doing a normal resampling. If rand==T, we are perumuting labels	
	if(!rand){
		sampfun <- function(pd){ pd[sample(nrow(pd), replace=TRUE),]}
	} else {
		sampfun <- function(pd){ytmp <- sample(pd$y, replace=TRUE); pd <- pd[sample(nrow(pd), replace=TRUE),]; pd$y <- ytmp; pd}
        }

	options(warn=-1)
	library(rpart)
	library(BatchJobs)
	source("functions.R")
	source("internal_looper_nomod.R")
	load("genesigprecision_data.Rda")

	# We break the data up into 25% "training", 75% testing, and only use the testing data
	pdidx <- sample(1:nrow(pd), nrow(pd)/4, replace=F)
        pdtrain <- pd[pdidx,]
        pdtest <- pd[-pdidx,]

	reg <- makeRegistry(id=input$id, seed=input$seed)	
	ids <- batchMap(reg, internal_looper_nomod, (input$seed+1):(input$seed+100), more.args=list(pdtest=pdtest, rand=rand))

	done <- submitJobs(reg, wait=function(retries) 100, max.retries=10)

	# To not step on the registry's toes while it's initializing
	Sys.sleep(5)

	jobs_done <- FALSE
	while(!jobs_done){
		jobs_done <- waitForJobs(reg,progressbar=FALSE)
	}

	do.call(rbind, loadResults(reg))
}
