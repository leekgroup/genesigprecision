#'
#' Parallelized Adjusted Estimator Simulation Function
#'
#' This function constructs a simulated dataset from the data saves in genesigprecision_data.Rda
#' and computes the adjusted and unadjusted treatment effect estimates. This process is repeated
#' 100 times. This function is called by BatchJobs 100 times to complete and aggregate results
#' of 10,000 simulations of the described form.
#'
#' @param seed An integer random seed
#' @param rand (logical) T indicates that you would like to permute all labels and remove corellation in the dataset. F indicates that you would like to permute records only and retain relationships between covariates
#'
#' @export
#'
#' @return A 100 x 15 matrix with results for W_{-age}, W_C, W_G, W_C;W_G, W_CG (described in paper) over 100 simulations

par_fun <- function(seed, rand=F){

	set.seed(seed)

	# If rand==F, we are doing a normal resampling. If rand==T, we are perumuting labels	
	if(!rand){
		sampfun <- function(pd){ pd[sample(nrow(pd), replace=TRUE),]}
	} else {
		sampfun <- function(pd){as.data.frame(lapply(pd, sample, replace=F))}
	}

	options(warn=-1)
	library(rpart)
	source("functions.R")
	load("genesigprecision_data.Rda")

	result <- matrix(NA,100,15)
	i <- 1

	while(i <= 100){

		pdtmp <- sampfun(pd) # We sample the data as per the scheme selected

		pdtmp$trt <- rbinom(nrow(pd), 1, 0.5) # Random assignment of treatment

		# Clinical model excluding age
		mod_clin_no_age <- rpart(y~Characteristics.TumorSize + Factor.Value.ER.status + g2ind + g3ind, data=pdtmp)
            mod_clin_no_age <- prune(mod_clin_no_age, cp=mod_clin_no_age$cptable[which.min(mod_clin_no_age$cptable[,"xerror"]),"CP"])
            pred_clin_no_age <- predict(mod_clin_no_age)

		# Clinical model including ER status
		mod_clin_age <- rpart(y~Characteristics.Age + Characteristics.TumorSize + Factor.Value.ER.status + g2ind + g3ind, data=pdtmp)
		mod_clin_age <- prune(mod_clin_age, cp=mod_clin_age$cptable[which.min(mod_clin_age$cptable[,"xerror"]),"CP"])
		pred_clin_age <- predict(mod_clin_age)

		# MammaPrint model
		mod_gen <- rpart(y~Factor.Value.MammaPrint.prediction, data=pdtmp)
		mod_gen <- prune(mod_gen, cp=mod_gen$cptable[which.min(mod_gen$cptable[,"xerror"]),"CP"])
		pred_gen <- predict(mod_gen)

		
		# Model including all predictors (all clinical + MammaPrint)
		mod_cpg <- rpart(y~Characteristics.Age + Characteristics.TumorSize + Factor.Value.ER.status + g2ind + g3ind + Factor.Value.MammaPrint.prediction, data=pdtmp)
		mod_cpg <- prune(mod_cpg, cp=mod_cpg$cptable[which.min(mod_cpg$cptable[,"xerror"]),"CP"])
		pred_cpg <- predict(mod_cpg)

		pdtmp$pred_clin_no_age <- pred_clin_no_age
		pdtmp$pred_clin_age <- pred_clin_age
		pdtmp$pred_gen <- pred_gen
		pdtmp$pred_cpg <- pred_cpg

		out_clin_no_age <- tryCatch(run_analysis(pdtmp, c("pred_clin_no_age")), error=function(e) e)
		out_clin_age <- tryCatch(run_analysis(pdtmp, c("pred_clin_age")), error=function(e) e)
		out_gen <- tryCatch(run_analysis(pdtmp, c("pred_gen")), error=function(e) e)
		# c_g is if we include seperate clinical and genetic predictors; cg is if we have all covarites together in one predictor
		out_c_g <- tryCatch(run_analysis(pdtmp, c("pred_clin_age", "pred_gen")), error=function(e) e)
		out_cg <-  tryCatch(run_analysis(pdtmp, c("pred_cpg")), error=function(e) e)
		
		if(any(unlist(lapply(list(out_clin_no_age, out_clin_age, out_gen, out_c_g, out_cg), inherits, "error")))){

			print("I made an error...check why.")
			next
		}
		
		result[i,] <- c(unlist(out_clin_no_age), unlist(out_clin_age), unlist(out_gen), unlist(out_c_g), unlist(out_cg))
		i <- i+1
	}
	
	result
}
