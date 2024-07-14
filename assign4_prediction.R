### ONLY past your predict.y function here
# predict.y function
predict.y <- function(x) {
  load("fit_params.Rdata")
  y_pred <- x %*% beta.final
  return(y_pred)
}


