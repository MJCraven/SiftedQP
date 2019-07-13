sift <- function(x_array, y_array) {
  # Original code by D. Graham
  # Identifies non-dominated data 
  # input data in xarray(1:iend) and yarray(1:iend)
  # output data in xenv(1:istart), yenv(1:istart)
  # iend = istart when program finishes
  istart <- 0
  iend <- length(x_array) # was 500
  pass <- 0
  maxpass <- 1 # should just need single pass
  errcnt <- 1
  pass <- 0
  options(digits=16) # Use more digits for calculation precision - DIG originally used 8
  while ( ( errcnt  > 0 ) & ( pass < maxpass ) ) {
    pass <- pass + 1
    #print(c(' Pass number ', pass ))
    ienv <- matrix(nrow=(iend+1), ncol=1, 0)
    xenv <- matrix(nrow=(iend+1), ncol=1, 0)
    yenv <- matrix(nrow=(iend+1), ncol=1, 0)
    index_array <- c(1:(iend)) # initialise index array
    
    istart <- 0
    ymax <- -1
    while ( istart < iend ) {
      
      dom <- matrix(nrow=iend+1, ncol=1, 0)
      istart <- istart+1
      
      # find max y
      ymax = -1e20
      for ( iy in istart:iend ) {
        yi <- y_array[iy]
        if ( yi >= ymax ) {
          ymax <- yi
          xmax <- x_array[iy]
          imax <- iy
        }
      }
      yenv[istart] <- ymax
      xenv[istart] <- xmax
      ienv[istart] <- imax
      if(dom[imax] == 2)
        print (c(' Error imax ', imax, ' dom ', dom[imax]))
      end  
      
      # swap to get max at istart entry
      temp <- xmax
      x_array[imax] <- x_array[istart]
      x_array[istart] <- temp
      temp <- index_array[istart]
      index_array[istart] <- index_array[imax]
      index_array[imax] <- temp
      temp <- ymax
      y_array[imax] <- y_array[istart]
      y_array[istart] <- temp

      dom[imax] <- dom[istart]
      dom[istart] <- 0

      # find next non-dominated
      if ( iend > 1 ) { # trivial case to avoid the case where we have just one point (nothing to compare it with!)
        for ( i in (istart+1):iend ) {
          dy1 <- (y_array[i]) - (y_array[istart])
          dx1 <- (x_array[i]) - (x_array[istart])
          dy <- dy1
          dx <- dx1
          if( ( dy <= 0 ) & ( dx >= 0 ) ) {
            dom[i] <- 2
          }
        } # end of for 
      } # end of if
      ifinish <- iend

      # remove dominated points
      for ( i in (istart+1):(ifinish-1) ) {
        while ( ( dom[i] == 2 ) & ( iend > istart ) & ( i < iend ) ) {
          for ( j in i:(ifinish - 1) ) {
            y_array[j] <- y_array[j+1]
            x_array[j] <- x_array[j+1]
            index_array[j] <- index_array[j+1]
            dom[j] <- dom[j+1]
          }
          iend <- iend - 1
        }
      } # end of for
      
      if((dom[ifinish]==2) & ( iend > istart )){
        iend <- iend-1
      }
      icount <- 0
      
      x_plot <- matrix(nrow=iend, ncol=1, 0)
      y_plot <- matrix(nrow=iend, ncol=1, 0)
      for ( i in 1:iend ) {
        if ( dom[i] == 0 ) {
          icount <- icount+1
          x_plot[icount] <- x_array[i]
          y_plot[icount] <- y_array[i]   
        }
      }
    } # end of while istart < iend
    
    print(icount)
    errcnt <- 1
    for ( i in 1:( istart-1 ) ) {
      dx <- xenv[i+1]-xenv[i]
      dy <- yenv[i+1]-yenv[i]
      if ( ( dx > 0 ) | ( dy > 0 ) ) {
        print(c('error in ordering, row ', ( i+1 ), ' dx, dy ', dx, ' ', dy, ' final istart =', istart,' iend =', iend))
        errcnt <- errcnt + 1
      }
    }
    if ( errcnt == 1 ) {
      errcnt <- 0
      print('NO ERRORS - CONVERGED ')
    }
    
    cf <- c(yenv[1:istart])
    df <- data.frame( cf, c(xenv[1:istart]))
    df2 <- melt(df, id.vars=1)
    
    # so rets come first followed by vols
    istart0 <- istart
    g2 <- ggplot(data=df2, aes(y=cf, x=value, color="red"))
    g2 <- g2 + geom_point()
        g2 <- g2 + labs(title=paste0("Non-Dominated Sub-EFs: Portfolio Ret vs Vol\n(", toString(istart0), " Points, ", data_string, ")"), x="Volatility", y="Return")
        g2 <- g2 + theme(legend.position="none")
    print(g2)
  } # end of while (errcnt > 0 && pass < maxpass)
  return(list("Returns" = yenv[1:istart], "Risks" = xenv[1:istart], "Indexes" = index_array[1:istart]))
}
