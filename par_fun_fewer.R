#'
#' Parallelized Adjusted Estimator Simulation Function
#'
#' This function constructs a simulated dataset from the data saves in genesigprecision_data.Rda
#' and computes the adjusted and unadjusted treatment effect estimates. This process is repeated
#' 100 times. This function is called by BatchJobs 100 times to complete and aggregate results
#' of 10,000 simulations of the described form. We split data into training and test sets and
#' resample from the test data. We get approximations of precision gain from the resampled test data
#' using covariates that are predictions from the models built on the training data.
#'
#' @param seed An integer random seed
#' @param rand (logical) T indicates that you would like to permute all labels and remove corellation in the dataset. F indicates that you would like to permute records only and retain relationships between covariates
#'
#' @export
#'
#' @return A 100 x 15 matrix with results for W_{-age}, W_C, W_G, W_C;W_G, W_CG (described in paper) over 100 simulations

par_fun_fewer <- function(seed,rand=F){

	# If rand==F, we are doing a normal resampling. If rand==T, we are perumuting labels	
	if(!rand){
		sampfun <- function(pd){ pd[sample(nrow(pd), replace=TRUE),]}
	} else {
		sampfun <- function(pd){ytmp <- sample(pd$y, replace=TRUE); pd <- pd[sample(nrow(pd), replace=TRUE),]; pd$y <- ytmp; pd}
	}

	options(warn=-1)
	source("functions.R")
	load("genesigprecision_data.Rda")

	result <- matrix(NA,100,12)
        i <- 1

        while(i <= 100){

                pdtmp <- sampfun(pd) # We sample the data as per the scheme selected

                pdtmp$trt <- rbinom(nrow(pdtmp), 1, 0.5) # Exogenous treatment assignment

		# Here, we only adjust for tumor size and ER status as clinincal covariates
		out_clin_no_er <- tryCatch(run_analysis(pdtmp, c("Characteristics.TumorSize")), error=function(e) e)
                out_clin_er <- tryCatch(run_analysis(pdtmp, c("Characteristics.TumorSize", "Factor.Value.ER.status")), error=function(e) e)
                out_gen <- tryCatch(run_analysis(pdtmp, c("Factor.Value.MammaPrint.prediction")), error=function(e) e)
                out_cg <- tryCatch(run_analysis(pdtmp, c("Characteristics.TumorSize", "Factor.Value.ER.status","Factor.Value.MammaPrint.prediction")), error=function(e) e)

                if(any(unlist(lapply(list(out_clin_no_er, out_clin_er, out_gen, out_cg), inherits, "error")))){

                        print("I made an error...check why.")
                        next
                }

                result[i,] <- c(unlist(out_clin_no_er), unlist(out_clin_er), unlist(out_gen), unlist(out_cg))
                i <- i+1
        }

        result
}

