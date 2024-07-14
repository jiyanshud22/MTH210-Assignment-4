bridgeReg <- function(y, X, alpha, lambda, maxiter, tol) {
  dist <- tol + 1
  iter <- 0
  ncols <- ncol(X)
  curr <- matrix(1, ncol = 1, nrow = ncols)
  
  while (dist > tol) {
    iter <- iter + 1
    if (iter > maxiter) {
      stop("Maximum iterations reached")
    }
    prev <- curr
    mjs <- alpha * abs(curr)^(alpha - 2)
    # Add regularization term to ensure invertibility
    curr <- qr.solve(t(X) %*% X + lambda / alpha * diag(ncols) + 1e-6 * diag(ncols), t(X) %*% y)
    dist <- norm(prev - curr, "2")
  }
  return(curr)
}

data <- read.csv("assign4_train.csv")
y <- data[, 1]
X <- data[, -1]
n <- nrow(X)
p <- ncol(X)

X <- cbind(1, X)
X <- as.matrix(X)

lambda_values <- 10^(seq(2.6, 3, by = 0.01))
CV_err <- numeric(length = length(lambda_values))

XX <- t(X) %*% X

for (i in 1:length(lambda_values)) {
  lambda <- lambda_values[i]
  track_cv <- 0
  for (j in 1:n) {
    X_train <- X[-j, , drop = FALSE]
    y_train <- y[-j]
    
    beta_train <- bridgeReg(y_train, X_train, 2, lambda, 1e5, .01)
    
    y_pred <- X[j, ] %*% beta_train
    
    track_cv <- track_cv + (y[j] - y_pred)^2
  }
  CV_err[i] <- track_cv / n
}

chosen_lambda <- lambda_values[which.min(CV_err)]

beta_final <- solve(XX + chosen_lambda * diag(dim(XX)[1])) %*% t(X) %*% y

save(beta_final, file = "fit_params.Rdata")
