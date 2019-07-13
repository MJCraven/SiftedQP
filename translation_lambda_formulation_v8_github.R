rm(list=ls()) # Clear workspace
library(openxlsx)
library(quadprog)
library(reshape2)
library(ggplot2)
library(raster)
library(plot.matrix)
library(RColorBrewer)
library(gplots)
setwd("C:/Users/mcraven1/Desktop/") # Working directory - Beasley datasets and Cesarone solutions should be in here.
source("get_EF_lambda.R")
source("sift_v3.R")
options(digits=20)
number_runs <- 1
time_comp <- matrix(nrow=number_runs, ncol=1)
time_sift <- matrix(nrow=number_runs, ncol=1)

# Parameters
sample_EFs <- 0 # = 1 if print samples of EFs, 0 otherwise
massets <- 3 # k-value
nlambda <- 200 # if we change this then we often obtain different numbers of cases
standard_datasets <- 1 # selects one of Beasley's (D.) datasets if == 1; else some other dataset
port_number <- 3 # Can change to 1, 2, 3, 4 or 5
vector_product_threshold <- 10^-18
lambda_0_vector <- c(0.342, 0.149, 0.426, 0.107, 0.115)

# Array pre-allocation limit - this should only really go up to max(icount_useful), but in practise we do not know what value that takes.
arr_size = 1000
allweights_size <- 100
do_pdf <- 1 # = 0 if no pdf, 1 if pdf

se0 <- as.numeric(Sys.time())
rn0 <- 1e8*(se0-floor(se0))

if ( standard_datasets == 1 ) {
  data_string <- paste0("Dataset (D", port_number, "), k = ", massets, ", l = ", nlambda, sep="")
  print(paste0("Using ", data_string, sep=""))
  table1 <- read.table(paste0("port", port_number, "-returns.txt", sep=""))
  mu_return_vector <- table1$V1
  sds <- table1$V2
  port_correlation <- read.table(paste0("port", port_number, "-corr.txt", sep=""))
  lambda_0 <- lambda_0_vector[port_number]
  
  lengthmu <- length(mu_return_vector)
  # Read in Cesarone dataset
  # From this data: # RetRisk??_SolVal.txt (column 1: return; columns 2-10: risk values when K=2,3,...,10; column 11: Unconstrained Markowitz risk) 
  Cesarone <- read.delim(paste0("RetRisk", lengthmu, "_SolVal.txt"), header=FALSE)
  
} else {
	# Use some other dataset - don't need sds as it is only used to compute the covariance matrix
	# So just need mu vector [mu_return_vector] and covariance matrix [sigma]
  alt_dataset <- 493 # 389, 493
  
  data_string <- paste0("Dataset (S&P-", alt_dataset, "), k = ", massets, ", l = ", nlambda, sep="")
  
  # Convert format
  mu_return_vector <- read.csv(paste0("returns-", alt_dataset, ".csv", sep=""), sep=",", header = FALSE)
	mu_return_vector <- as.numeric(as.vector(mu_return_vector))
	sigma1 <- read.delim(paste0("covariance matrix-", alt_dataset, ".txt", sep=""), header= FALSE)

	if ( alt_dataset == 389 ) {
	  lambda_0 <- 0.0000001
	  lambda_1 <- 1
	  load("389-unc_EF.RData")
	  uef_rets <- df1$Soln_y
	  uef_risks <- df1$Soln_x
	} else {
	  lambda_0 <- 0.489
	  lambda_1 <- 1
	  load("493-unc_EF.RData")
	  uef_rets <- df1$Soln_y
	  uef_risks <- df1$Soln_x
	}
}

number_assets <- length(mu_return_vector)

rdmtitle <- paste0(rn0, "-n=",  number_assets, "-k=", massets, "-nlambda=", nlambda)
if ( do_pdf == 1 ) {
  pdf(paste0(rdmtitle, ".pdf",sep=""), onefile=T, paper="A4r")
}

for ( run0 in 1:number_runs ) {
  print(paste0("Run ", run0, sep=""))
  
  Time1 <- proc.time()[3]
  
if ( standard_datasets == 1 ) {
  lambda_1 <- 1
} else {
  #
}
dlambda <- (lambda_1 - lambda_0)/(nlambda-1)

run_all_lambda <- 1

if ( standard_datasets == 1 ) {
	# Get correlations matrix
	corr_matrix <- matrix(nrow=number_assets, ncol=number_assets)
	for ( i in 1: length(port_correlation$V1) ) {
	  row1 <- port_correlation$V1[i]
	  col1 <- port_correlation$V2[i]
	  corr_matrix[row1,col1] <- port_correlation$V3[i]
	  corr_matrix[col1,row1] <- port_correlation$V3[i] # And the other half, so corr_matrix is symmetric
	}

	# Construct covariance matrix
	sigma <- matrix(nrow=number_assets, ncol=number_assets)
	for ( i in 1:number_assets ) {
	  for ( j in 1:number_assets ) {
	    sigma[i,j] <- corr_matrix[i,j]*sds[i]*sds[j]
	  }
	}
} else {
  # Put covariance matrix in proper format
    sigma <- matrix(nrow=number_assets, ncol=number_assets)
    for ( i in 1:number_assets ) {
      command1 <- paste0("sigma[,i] <- sigma1$V", i)
      eval(parse(text=command1))
      }
}

list_combination <- function(n, k, indx) {
  thislot_start <- 0
  thislot_end <- 0
  num_combs <- choose(n,k)
  j <- 0
  caselist <- rep(0,k) # k zeroes
  
  count <-0
  lastlot <- 0
  for ( kindex in 1:k ) {
    found <- 0
    k1 <- k - kindex
    while ( ( found == 0 ) & ( thislot_end <= num_combs ) ) {
      j <- j + 1
      count <- count + 1
      thislot_start <- thislot_end
      nextlot <- choose(n-j, k1)
      thislot_end <- thislot_start + nextlot
      if ( ( indx <= thislot_end ) & ( indx > thislot_start ) ) {
        thisdigit <- count
        found<- 1
        caselist[kindex] <- thisdigit
      }
    }
    thislot_end <- thislot_start
  }
  return(caselist)
}

Q2 <- sigma
returns2 <- mu_return_vector
nsamp <- length(returns2)
Q1 <- Q2

# Re-order in order of return [max = 1]:
nassets <- nsamp
bigret <- matrix(nrow=nassets, ncol=1, -10)
bigretindex <- matrix(nrow=nassets, ncol=1, 0)
iselect <- matrix(nrow=nassets, ncol=1, seq(1,nassets))
returns3 <- returns2
returns1 <- returns2

# Re-order returns vector
for ( i in 1:nassets ) {
  bigret[i] <- -10
  for ( j in i:nassets ){
    
    if ( bigret[i] < returns3[j] ) {
      bigret[i] <- returns3[j]
      bigretindex[i] <- j
    }
  }
  tempret <- returns3[bigretindex[i]]
  returns3[bigretindex[i]] <- returns3[i]
  returns3[i] <- tempret
  tempindex <- iselect[i]
  iselect[i] <- iselect[bigretindex[i]]
  iselect[bigretindex[i]] <- tempindex
}

# Re-order risk matrix
Q1<- matrix(nrow=nassets, ncol=nassets, 0)
for ( jcol in 1:nassets ) {
  for ( irow in 1:nassets ) {
    Q1[irow,jcol] <- Q2[iselect[irow], iselect[jcol]]
  }
}

n_possibles<- choose(nassets, massets) #number of possible combinations

# Initialise matrices - in same cases, lists may be used.
vol2 <- matrix(nrow=nlambda, ncol=arr_size, 0)
retn2 <- matrix(nrow=nlambda, ncol=arr_size,0)
cost <- matrix(nrow=nlambda, ncol=arr_size, 0)
mincost <- matrix(nrow=nlambda, ncol=1, 1000)  #
rundata <- matrix(nrow=nlambda, ncol=1, 0)     #
minrisk <- matrix(nrow=n_possibles, ncol=1, 0) #
maxrisk <- matrix(nrow=n_possibles, ncol=1, 0) #
maxret <- matrix(nrow=n_possibles, ncol=1, 0)  #
minret <- matrix(nrow=n_possibles, ncol=1, 0)  #
caselist_useful <- matrix(nrow=massets, ncol=arr_size, 0)

allweights <- array(dim = c(nlambda, massets,allweights_size),0) # array containing icount_useful weights matrices [each size nlambda*massets]  

massets1 <- massets
iselect <- seq(1,nassets)
icount_useful <- 0

for ( irun in 1:n_possibles ) { # Loop through all the possibilities
  if ( irun %% 10000 == 0 ){
    print(c("Case = ", irun, "/", n_possibles))
  }
  r<-matrix(nrow=massets, ncol=1, 0)
  
  # Compute unconstrained EF if needed
  if (irun > n_possibles) {
    massets <- nassets
    Q <- Q1
    returns <- returns1
  } else {
    # Identify assets that form sub-portfolio
    caselist <- list_combination(nassets, massets, irun)
    iselect <- caselist
    returns <- matrix(nrow=massets, ncol=1, 0)
    Q <- matrix(nrow=massets, ncol=massets, 0)
    
    # Write data into returns vector and Q matrix
    for( jcol in 1:massets ){
      returns[jcol] <- returns3[iselect[jcol]]
      for( irow in 1:massets ){
        Q[irow,jcol] <- Q1[iselect[irow], iselect[jcol]]
      }
    }
    
    # Check if `current` sub-EF would be dominated by previously-computed ones (analytical for k=2)
    if ( massets == 0 ) {
      x0<- -(Q[1,2]-Q[2,2])/(Q[1,1]+Q[2,2]-(2*Q[1,2]))
      if ((x0 > 0) & (x0 < 1)){
        minrisk[irun] <- ((Q[1,1]*Q[2,2])-(Q[1,2]*Q[1,2])) / ( Q[1,1]+Q[2,2]-2*Q[1,2] )
        minret[irun] <- (1/(Q[1,1]+Q[2,2]-(2*Q[1,2])))*(Q[1,1]*returns[2]-Q[1,2]*(returns[1]+returns[2])+Q[2,2]*returns[1])
        maxret[irun] <- max(returns)
      } else {
        minrisk[irun] <- min(Q[1,1],Q[2,2])
        if ( Q[1,1] > Q[2,2] ) {
          minret[irun] <- returns[2]
        } else {
          minret[irun] <- returns[1] 
        }
      }
      
      if ( returns[1] > returns[2]) {
        maxrisk[irun] <- Q[1,1]
        maxret[irun] <- returns[1]
      } else {
        maxrisk[irun] <- Q[2,2]
        maxret[irun] <- returns[2]
      } 
    } else { # if massets != 0 (was != 2)
      
      # Find minimum risk and the corresponding return:
      # We have the returns vector in descending numerical order
      # Here we are solving the Markowitz problem but only for one point on the EF: where $\lambda = 0$.
      # i.e., max return.
      
      #bound <- max(returns)
      l0 <- get_EF_lambda(number_points=1, lambda_value=lambda_1, num_assets=massets, sigma0=Q, returns0=returns )
      minrisk[irun] <- l0$risks
      soln <- l0$solutions
      minret[irun] <- l0$return
      thismaxret <- -1000
      for ( iasset in 1:massets ) {
        thisretn <- returns[iasset]
        if ( thisretn > thismaxret ) {
          thismaxret <- thisretn
          thismaxrisk <- Q[iasset,iasset]
          thisasset <- iasset
        }
      }
      
      # Find max return and corresponding risk
      # assumes max return is first element
        maxrisk[irun] <- thismaxrisk
        maxret[irun] <- thismaxret
    }
  } # end else irun > n_possibles
  
  r <- returns
  
  if ( irun > 1 ) {
    run_all_lambda <- 0
    
    # Checking for domination - general case
    ilambda <- 0
    while ( run_all_lambda != 1 && ilambda < nlambda-1 ) {
      ilambda <- ilambda + 1
      # Check min risk point: vol2 and retn2 contain `current` overall EF
      vec1x <- vol2[ilambda, rundata[ilambda]]-minrisk[irun]
      vec1y <- retn2[ilambda, rundata[ilambda]]-minret[irun]
      vec2x <- vol2[ilambda+1, rundata[ilambda+1]]-minrisk[irun]
      vec2y <- retn2[ilambda+1, rundata[ilambda+1]]-minret[irun]
      crossprod1 <- (vec1x*vec2y)-(vec2x*vec1y)

      # Check max return point
      vec1x <- vol2[ilambda, rundata[ilambda]]-maxrisk[irun]
      vec1y <- retn2[ilambda, rundata[ilambda]]-maxret[irun]
      vec2x <- vol2[ilambda+1, rundata[ilambda+1]]-maxrisk[irun]
      vec2y <- retn2[ilambda+1, rundata[ilambda+1]]-maxret[irun]
      crossprod2 <- (vec1x*vec2y)-(vec2x*vec1y)

      # Check max some middle point (this is an approximation)
      # The numbers below may change based upon the dataset
      vec1x <- vol2[ilambda, rundata[ilambda]]-(.75*minrisk[irun]+.25*maxrisk[irun])
      vec1y <- retn2[ilambda, rundata[ilambda]]-(0.5*maxret[irun]+0.5*minret[irun])
      vec2x <- vol2[ilambda+1, rundata[ilambda+1]]-(.75*minrisk[irun]+.25*maxrisk[irun])
      vec2y <- retn2[ilambda+1, rundata[ilambda+1]]-(0.5*maxret[irun]+0.5*minret[irun])
      crossprod3 <- (vec1x*vec2y)-(vec2x*vec1y)
      
      volstor1 <- vol2[ilambda, rundata[ilambda]]
      volstor2 <- vol2[ilambda+1, rundata[ilambda+1]]
      retstor1 <- retn2[ilambda, rundata[ilambda]]
      retstor2 <- retn2[ilambda+1, rundata[ilambda+1]]

      if ( ( crossprod1 < -vector_product_threshold ) | ( crossprod2 < -vector_product_threshold ) | ( crossprod3 < -vector_product_threshold ) ) {
        run_all_lambda <- 1
      } else {
        #print("Dominated!")
      }
    } # end of while
  }
  
  if ( sample_EFs == 1 ) {
    # Identify case and plot EF to give a WW-representation
    if ( ( caselist[1] %% 10 == 1 ) && ( caselist[2] %% 5 == 1 ) ) {
      run_all_lambda <- -1 # 1
    }
  }
  
  if ( irun == 1 | abs( run_all_lambda ) == 1 ) { 
    # We have found a non-dominated sub-EF - get full sub-EF for all lambda_i by QP
    icount_useful <- icount_useful+1

    if ( icount_useful > arr_size ) {
      print("Adding another column.")
      vol2 <- cbind(vol2, rep(0, nlambda))
      retn2 <- cbind(retn2, rep(0, nlambda))
      cost <- cbind(cost, rep(0, nlambda))
      empty <- array(0, dim=c(nlambda,massets))
      allweights <- abind(allweights, empty, along=3) # add extra 'columns' [are we able to do this?]
      caselist_useful <-cbind(caselist_useful, rep(0, massets))
    }
    
    lambda <- lambda_0
    weights <- matrix(nrow=nlambda, ncol=massets, 0) # We can use these if we wish
    
    lowest_return = minret[irun]

    for ( i in 1:nlambda ) {
      # Do QP procedure for each point on the full sub-EF
      l0 <- get_EF_lambda(number_points=1, lambda_value=lambda_0+dlambda*(i-1), num_assets=massets, sigma0=Q, returns0=returns )
      soln <- l0$solutions
      vol2[i, icount_useful] <- l0$risks
      retn2[i, icount_useful] <- l0$return #soln
      cost[i, icount_useful] <- (lambda*vol2[i, icount_useful]) - ((1-lambda)*retn2[i, icount_useful])

      # In this case the same. Usually we would be minimising $\lambda r - (1-\lambda)\overline{R}$, but in this case we are minimising risk r.
      weights[i, ] <- l0$solutions # Row i of weights 'vector' - do not currently do anything with these.
      
      allweights[i, ,icount_useful] <- l0$solutions
      
      # Mincost[i] is the `current' min cost for lambda_i
      lambda <- lambda + dlambda # Go to the next lambda value - rewritten in terms of returns below.
      if ( cost[i, icount_useful] < mincost[i] ) {
        mincost[i] <- cost[i, icount_useful]
        rundata[i] <- icount_useful
      }
    } # end for for i
  } # end of if irun
  
  # Plots
  
  if ( irun <= n_possibles ) {
    # Generates plots of sub-EF
    if ( abs(run_all_lambda) == 1 | irun == 1 ) {
      print( c( "irun = ", irun, "case = ", caselist, "count = ", icount_useful ) )
      caselist_useful[ ,icount_useful] <-caselist 
      tempc = caselist
     } 
  } else {
    # Sequence of plots
    print( c("irun = ", irun, "unconstrained") )
     }
} # end of for irun

massets <- massets1

# Basic statistics
avcost <- mean(mincost)
retvec <- matrix(nrow=icount_useful*nlambda, ncol=1, 0)
riskvec <- matrix(nrow=icount_useful*nlambda, ncol=1, 0)
iveccount<- 0

# Collate data from each sub-EF to export to sieving code
for ( j in 1:icount_useful ) {
  for ( i in 1:nlambda ) {
    iveccount <- iveccount+1
    retvec[iveccount] <- retn2[i,j]
    riskvec[iveccount] <- vol2[i,j]
  }
}

Time2 <- proc.time()[3]
print(c("Computation time: ", Time2-Time1))

vols0 <- vol2[1:nlambda,1:icount_useful]
rets0 <- retn2[1:lambda,1:icount_useful]
plot(diag(Q1),returns3[1:nassets],col="black",type="p",xlim = c(min(vols0),max(max(vols0),max(diag(Q1)))), ylim = c(min(returns3[1:nassets]),max(rets0)),xlab="Volatility",ylab="Return",main="")
legend_text <- 1:icount_useful
for ( i in 1:icount_useful ) {
  lines(vol2[ ,i], retn2[ ,i], col="red", lwd=2)
}
green_CCEF=matrix(nrow = nlambda,ncol = 2)
for ( i in 1:nlambda ) {
  green_CCEF[i,1] <- vol2[i,rundata[i]]  # x co-ords
  green_CCEF[i,2] <- retn2[i,rundata[i]] # y co-ords
}
points(green_CCEF, pch = 19, col="green", bg="green")

# Set ggplot2 title format to centred
theme_update(plot.title = element_text(hjust = 0.5))
cf <- c(retn2[1:nlambda,1:icount_useful])
df <- data.frame( cf, c(vol2[1:nlambda,1:icount_useful]))
df2 <- melt(df, id.vars=1)
g2 <- ggplot(data=df2, aes(y=cf, x=value, color="red"))
g2 <- g2 + geom_point()
g2 <- g2 + labs(title=paste0("Non-Dom'd Sub-EFs: Pf Ret vs Vol [", icount_useful, " EFs]", sep=""), x="Volatility", y="Return")
green_CCEF <- data.frame(green_CCEF)
green_CCEF <- melt(green_CCEF, id.vars=1)
g2 <- g2 + geom_point(data=green_CCEF, aes(x=X1, y=value, color="green"))
g2 <- g2 + theme(legend.position="none")
print(g2)

# Now pass all the points from all sub-EFs to the sifting algorithm
print(c("Sifting ", icount_useful, " EFs."))
non_dom_EF <- sift(x_array=riskvec, y_array=retvec)

nnondom <- length(non_dom_EF$Risks)  # last entry spurious??
nondom_weights <- matrix(nrow = nnondom ,ncol=nassets)

#print(c(non_dom_EF$Indexes,non_dom_EF$Returns))
for ( k in 1:nnondom ) {
  global_index <- non_dom_EF$Indexes[k] 
  row_index <- (global_index-1)%/%nlambda+1 # which sub-EF
  col_index <- (global_index-1)%%nlambda+1  # which lambda?
  for (asset in 1:massets){
  nondom_weights[k,caselist_useful[asset,row_index]] <- allweights[col_index,asset,row_index]
  }
}

                          
if ( run0 == 1 ) { ###
  write.xlsx(non_dom_EF, file=paste(rn0, "n", nassets, "_k", massets, "_lambda", nlambda, "_CCEF.xlsx", sep=""))
  saveRDS(non_dom_EF, file=paste(rn0, "n", nassets, "_k", massets, "_lambda", nlambda, "_CCEF.rds", sep=""))
} ###

print(c("Ratio: ", icount_useful, "/", irun, "=", icount_useful/irun))
Time3 <- round(proc.time()[3] - Time1,digits=2)
print(c("Total time taken: ", Time3))

time_sift[run0] <- Time3-Time2+Time1
time_comp[run0] <- Time2-Time1

print(c("Sifting time: ", time_sift[run0]))
print(c("Computation time: ", time_comp[run0]))

# Plot with Cesarone CCEF if we have a standard dataset
# Close up
if ( standard_datasets == 1 ) {
  plot(non_dom_EF$Risks, non_dom_EF$Returns, type="n", xlab="Risk", ylab="Return", main="")
  points(non_dom_EF$Risks, non_dom_EF$Returns, col="red", lwd=1)
  lines(Cesarone[, massets], Cesarone[,1], col="blue", lwd=3)
  legend("bottomright", legend=c(paste0("NDCCEF (D", port_number, ")", sep=""), paste0("Cesarone CCEF (k=", massets, ")", sep="")), col=c("red", "blue"), lty=c(1,1,1), cex=0.8)
}

# Proper scale #######
if ( standard_datasets == 1 ) {
  points(non_dom_EF$Risks, non_dom_EF$Returns, col="red", lwd=1)
  lines(Cesarone[, massets], Cesarone[,1], col="blue", lwd=3)
  legend("bottomright", legend=c(paste0("NDCCEF (D", port_number, ")", sep=""), paste0("Cesarone CCEF (k=", massets, ")", sep="")), col=c("red", "blue"), lty=c(1,1,1), cex=0.8)
}

if ( standard_datasets == 0 ) {
  plot(non_dom_EF$Risks, non_dom_EF$Returns, type="n", xlim=c(min(uef_risks),max(uef_risks)), ylim=c(min(uef_rets),max(uef_rets)),       xlab="Risk", ylab="Return", main="")
  points(non_dom_EF$Risks, non_dom_EF$Returns, col="red", lwd=1)
  lines(uef_risks, uef_rets, col="blue", lwd=3)
  legend("bottomright", legend=c(paste0("NDCCEF (", alt_dataset, ")", sep=""), paste0("UEF", sep="")), col=c("red", "blue"), lty=c(1,1,1), cex=0.8)
}

par(mar=c(5.1, 4.1, 4.1, 4.1))
plot(nondom_weights,col=rainbow(5, start=0, end=1), key=list(cex.axis=1))

# Plot restricted heatmap
# nnondom = number of non-dominated sub-EFs
# Compute maximal index of non-zero weight

max_index <- 0
for ( i in 1:nnondom ) {
   max_index_comp <- max(which(!is.na(nondom_weights[i,]))) # max index of non-NA weight on row i
   if ( max_index_comp > max_index ) {
     max_index <- max_index_comp
   }
}

lastbit <- nnondom-20+1

# Another try
plot(rbind(nondom_weights[1:20,1:max_index],nondom_weights[lastbit:nnondom,1:max_index]),breaks=NULL, main="Weights of top and bottom 20 points on CEF by return", xlab = "Asset number (in order of return)", col=rainbow(6, start=0, end=1), key=list(cex.axis=1))

# Samples every 5
heatmap.2( nondom_weights, Rowv=FALSE, Colv=FALSE, dendrogram='none', density.info='none', trace='none', key=TRUE, symkey=FALSE,lhei = c(2,6), col=brewer.pal(11,"Spectral"), margins = c(5,10 ))
heatmap.2( rbind(nondom_weights[1:20,1:max_index],nondom_weights[lastbit:nnondom,1:max_index]), Rowv=FALSE, Colv=FALSE, dendrogram='none', density.info='none', trace='none', key=TRUE, symkey=FALSE,lhei = c(2,6), col=brewer.pal(11,"Spectral"), margins = c(5,10 ))
nondom_weights2 <- nondom_weights
nondom_weights2[nondom_weights2 < 10^-10] <- NA

# Samples every 10 points
if ( max_index > 10 ) {
labCol1 <- c(rep(NA,max_index))
labRow1 <- c(rep(NA,nnondom))
labCol1[seq(1,max_index,floor(max_index/10))] <- seq(1,max_index,floor(max_index/10))
labCol1[max_index] <- max_index
labRow1[seq(1,nnondom,floor(nnondom/100))] <- c(1,10*seq(1+floor(nnondom/100),nnondom,floor(nnondom/100)))
labRow1[nnondom] <- 10*nnondom
heatmap.2( rbind(nondom_weights2[seq(1,nnondom,10),1:max_index],nondom_weights2[nnondom,1:max_index]), Rowv=FALSE, Colv=FALSE, dendrogram='none', density.info='none', trace='none', key=TRUE,
  key.xtickfun=function() {
  cex <- par("cex")*par("cex.axis")
  side <- 1
  line <- 0
  col <- par("col.axis")
  mtext(0, side=side, at=0, adj=0, line=line, cex=cex, col=col)
  mtext(1, side=side, at=1, adj=1, line=line, cex=cex, col=col)
  return(list(labels=FALSE, tick=FALSE))
}, xlab="Rank of asset by return", ylab="CCEF point number",
labRow = labRow1,
labCol = labCol1,
col=brewer.pal(11,"Spectral"), #margins = c(5,10 )
)
}
} # end of run

print(c("Mean/SD sifting time: ", mean(time_sift), sd(time_sift)))
print(c("Mean/SD computation time: ", mean(time_comp), sd(time_comp)))

if ( run0 == 1 ) { ###
  write.xlsx(nondom_weights, file=paste(rn0, "-n", nassets, "_k", massets, "_lambda", nlambda, "_nd-weights.xlsx", sep=""))
} ###
dev.off()
dev.off()
