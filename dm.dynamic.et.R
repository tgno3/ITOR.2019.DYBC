dm.dynamic.et <- function(xdata, ydata, zdata, finalz, rts = "crs", orientation = "i"){
  # Load library
  #library(lpSolveAPI)
  
  # Initial checks
  if(dim(xdata)[1] != dim(ydata)[1] | dim(xdata)[length(dim(xdata))] != dim(ydata)[length(dim(ydata))]) stop('Data must be balanced.')
  if(!(rts %in% c("crs", "vrs", "irs", "drs")))                                                         stop('rts must be "crs", "vrs", "irs", or "drs".')
  if(!(orientation %in% c("i", "o")))                                                                   stop('orientation must be "i" or "o".')
  
  # Parameters
  xdata    <- if(length(dim(xdata)) != 3)  array(xdata, c(dim(xdata)[1], 1, dim(xdata)[2])) else as.array(xdata)
  ydata    <- if(length(dim(ydata)) != 3)  array(ydata, c(dim(ydata)[1], 1, dim(ydata)[2])) else as.array(ydata)
  zdata    <- if(length(dim(zdata)) != 3)  array(zdata, c(dim(zdata)[1], 1, dim(zdata)[2])) else as.array(zdata)
  finalz   <- if(length(dim(finalz)) != 2) array(finalz, c(length(finalz), 1))              else as.array(finalz)
  initialz <- finalz + apply(zdata, 1, sum)
  n        <- dim(xdata)[1]
  m        <- dim(xdata)[2]
  s        <- dim(ydata)[2]
  b        <- dim(zdata)[2]
  t        <- dim(xdata)[3]
  
  # Data frames
  results.efficiency.s <- array(NA, dim = c(n, 1))
  results.lambda       <- array(NA, dim = c(n, n))
  results.efficiency.t <- array(NA, dim = c(n, t))
  results.xslack       <- array(NA, dim = c(n, m, t))
  results.yslack       <- array(NA, dim = c(n, s, t))
  results.zslack       <- array(NA, dim = c(n, b, t))
  results.islack       <- array(NA, dim = c(n, b, t))
  results.fslack       <- array(NA, dim = c(n, b, t))
  
  
  # Pointers
  p.eff <- n + 1
  p.xsl <- n + t + 1
  p.zsl <- n + t + m*t + 1
  p.ysl <- n + t + m*t + b*t + 1
  p.isl <- n + t + m*t + b*t + s*t + 1
  p.fsl <- n + t + m*t + b*t + s*t + b + 1
  p.end <- n + t + m*t + b*t + s*t + b + b + 1
  
  # LP
  for(j in 1:n){
    # Declare LP
    lp.it <- make.lp(0, n + t + m*t + b*t + s*t + b + b) # lambda + eff + xslack + zslack + yslack + initialZslack + finalZslack
    
    # Set objective
    if(orientation == "i") set.objfn(lp.it, c(rep( 1/t, t)), indices = c(p.eff:(p.xsl - 1)))
    if(orientation == "o") set.objfn(lp.it, c(rep(-1/t, t)), indices = c(p.eff:(p.xsl - 1)))
    
    # RTS
    if(rts == "vrs") add.constraint(lp.it, c(rep(1, n)), indices = c(1:n), "=", 1)
    if(rts == "crs") set.constr.type(lp.it, 0, 1)
    if(rts == "irs") add.constraint(lp.it, c(rep(1, n)), indices = c(1:n), ">=", 1)
    if(rts == "drs") add.constraint(lp.it, c(rep(1, n)), indices = c(1:n), "<=", 1)
    
    # LP
    for(k in 1:t){
      
      # Input constraint
      for(i in 1:m){
        if(orientation == "i"){
          add.constraint(lp.it, c(xdata[, i, k], -xdata[j, i, k], 1), indices = c(1:n, n + k, p.xsl - 1 + m*(k - 1) + i), "=", 0)
        }else{
          add.constraint(lp.it, c(xdata[, i, k], 1), indices = c(1:n, p.xsl - 1 + m*(k - 1) + i), "=", xdata[j, i, k])
        }
      }
      
      # Stock constraint
      for(i in 1:b){
        if(orientation == "i"){
          add.constraint(lp.it, c(zdata[, i, k], -zdata[j, i, k], 1), indices = c(1:n, n + k, p.zsl - 1 + b*(k - 1) + i), "=", 0)
        }else{
          add.constraint(lp.it, c(zdata[, i, k], 1), indices = c(1:n, p.zsl - 1 + b*(k - 1) + i), "=", zdata[j, i, k])
        }
      }
      
      # Output constraint
      for(r in 1:s){
        if(orientation == "i"){
          add.constraint(lp.it, c(ydata[, r, k], -1), indices = c(1:n, p.ysl - 1 + s*(k - 1) + r), "=", ydata[j, r, k])
        }else{
          add.constraint(lp.it, c(ydata[, r, k], -ydata[j, r, k], -1), indices = c(1:n, n + k, p.ysl - 1 + s*(k - 1) + r), "=", 0)
        }
      }
    }
    
    # Initial stock constraints
    for(i in 1:b){
      add.constraint(lp.it, c(initialz[, i], 1), indices = c(1:n, p.isl - 1 + i), "=", initialz[j, i])
    }
    
    # Final stock constraints
    for(i in 1:b){
      add.constraint(lp.it, c(finalz[, i], -1), indices = c(1:n, p.fsl - 1 + i), "=", finalz[j, i])
    }
    
    # Bounds
    if(orientation == "i"){
      set.bounds(lp.it, upper = c(rep(Inf, n), rep(   1, t), rep(Inf, p.end - p.xsl)))  
      set.bounds(lp.it, lower = c(rep(0,   n), rep(-Inf, t), rep(  0, p.end - p.xsl)))
    }else{
      set.bounds(lp.it, lower = c(rep(0,   n), rep(   1, t), rep(  0, p.end - p.xsl)))
    }
    
    # Solve
    solve.lpExtPtr(lp.it)
    
    # Get results
    results.efficiency.s[j]   <- abs(get.objective(lp.it))
    temp.p                    <- get.variables(lp.it)
    results.efficiency.t[j, ] <- temp.p[p.eff:(p.xsl - 1)]
    
    # Stage II
    # Link previous solutions
    for(k in 1:t){add.constraint(lp.it, c(1), indices = c(p.eff + k - 1), "=", results.efficiency.t[j, k])}
    
    # Slack sum max
    set.objfn(lp.it, c(rep(-1, (p.end - p.xsl))), indices = c(p.xsl:(p.end - 1)))
    
    # Solve
    solve.lpExtPtr(lp.it)
    
    # Get results
    temp.s              <- get.variables(lp.it)
    results.lambda[j, ] <- temp.s[1:n]
    results.xslack[j,,] <- array(temp.s[p.xsl:(p.zsl - 1)], c(m, t))
    results.zslack[j,,] <- array(temp.s[p.zsl:(p.ysl - 1)], c(b, t))
    results.yslack[j,,] <- array(temp.s[p.ysl:(p.isl - 1)], c(s, t))
    results.islack[j,,] <- array(temp.s[p.isl:(p.fsl - 1)], c(b, t))
    results.fslack[j,,] <- array(temp.s[p.fsl:(p.end - 1)], c(b, t))
  }
  
  # Store results
  results <- list(eff.s  = results.efficiency.s, 
                  eff.t  = results.efficiency.t, 
                  lambda = results.lambda, 
                  xslack = results.xslack, 
                  yslack = results.yslack, 
                  zslack = results.zslack, 
                  islack = results.islack, 
                  fslack = results.fslack)
  return(results)
}
