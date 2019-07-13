get_EF_lambda <- function(number_points, lambda_value, num_assets, sigma0, returns0) {
  # A function to get the EF based upon the following input.
  # An expansion of the simple example in
  # http://horicky.blogspot.co.uk/2013/01/optimization-in-r.html
  
  # Input: number_points (number of points we want in the EF)
  #        lambda_value (size of expected returns required)
  #        num_assets (number of assets)
  #        sigma (covariance matrix for the num_assets subset)
  #        returns0 (returns for the num_assets subset)
  # Ouput: solution - weights
  #        min cost / risk vector (the same here, as we are solving the 2-obj problem rather than the \lambda-formulation)
  
  options(digits=16) # Use more digits for calculation precision

  D.Matrix <- 2*sigma0*lambda_value
  d.Vector <- returns0*(1-lambda_value)
  
  sumX.Equality <- matrix(rep(1,num_assets), ncol=1) 
  # Constraints matrix A:
  A.Matrix <- cbind(sumX.Equality, diag(num_assets))
  risk.Vector <- rep(0,number_points)
  count <- 1 # Count from 1
  solutions <- matrix(nrow = number_points, ncol = num_assets)

  # This is the RHS of Ax >= b:
  b.Vector <- c(1, rep(0,num_assets))

  out <- solve.QP(Dmat=D.Matrix, dvec=d.Vector, Amat=A.Matrix, bvec=b.Vector, meq=1)
  soln = out$solution
  risk.Vector[count] <- t(soln) %*% sigma0 %*% soln
  solutions[count,] <- soln # record the asset allocation
  count <- count + 1
  returns.Number <- soln %*% returns0
  
  return(list("solutions" = solutions, "risks" = risk.Vector, "return" = returns.Number))
}
