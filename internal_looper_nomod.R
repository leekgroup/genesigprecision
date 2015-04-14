internal_looper_nomod <- function(seed, pdtest, rand){
	#set.seed(seed)

	# If rand==F, we are doing a normal resampling. If rand==T, we are perumuting labels	
	if(!rand){
		sampfun <- function(pd){ pd[sample(1:nrow(pd), size=296, replace=TRUE),]}
	} else {
		sampfun <- function(pd){ytmp <- sample(pd$y, replace=TRUE); pd <- pd[sample(nrow(pd), replace=TRUE),]; pd$y <- ytmp; pd}
        }

	library(rpart)
	source("functions.R")
	result <- matrix(NA,100,12)
	i <- 1

	while(i <= 100){

		pdtmp <- sampfun(pdtest) # We sample the data as per the scheme selected

		pdtmp$trt <- rbinom(nrow(pdtmp), 1, 0.5) # Random assignment of treatment

		out_clin_no_er <- tryCatch(run_analysis(pdtmp, c("Characteristics.Age", "Characteristics.TumorSize", "g2ind", "g3ind")), error=function(e) e)
                out_clin_er <- tryCatch(run_analysis(pdtmp, c("Characteristics.Age", "Characteristics.TumorSize", "Factor.Value.ER.status", "g2ind", "g3ind")), error=function(e) e)
                out_gen <- tryCatch(run_analysis(pdtmp, c("Factor.Value.MammaPrint.prediction")), error=function(e) e)
                out_cg <- tryCatch(run_analysis(pdtmp, c("Characteristics.Age", "Characteristics.TumorSize", "Factor.Value.ER.status", "g2ind", "g3ind", "Factor.Value.MammaPrint.prediction")), error=function(e) e)

                if(any(unlist(lapply(list(out_clin_no_er, out_clin_er, out_gen, out_cg), inherits, "error")))){

                        print("I made an error...check why.")
                        next
                }

                result[i,] <- c(unlist(out_clin_no_er), unlist(out_clin_er), unlist(out_gen), unlist(out_cg))
                i <- i+1
	}

	result
}

