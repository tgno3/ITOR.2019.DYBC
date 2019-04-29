dm.dynamic.ll <- function(xdata, ydata, zdata, budget, rts = "crs", orientation = "i"){
  # Load library
  #library(lpSolveAPI)
  
  # Initial checks
  #if(!(3 %in% c(length(dim(xdata)), length(dim(ydata)))))              stop('Data must be 3-dimensional.')
  #if(dim(xdata)[length(dim(xdata))] != dim(ydata)[length(dim(ydata))]) stop('Data must be balanced.')
  if(!(rts %in% c("crs", "vrs", "irs", "drs")))                        stop('rts must be "crs", "vrs", "irs", or "drs".')
  if(!(orientation %in% c("i", "o", "n")))                             stop('orientation must be "i", "o", or "n".')
  
  # Parameters
  xdata  <- if(length(dim(xdata))  != 3) array( xdata,  c(dim(xdata)[1], 1, dim(xdata)[2])) else as.array(xdata)
  ydata  <- if(length(dim(ydata))  != 3) array( ydata,  c(dim(ydata)[1], 1, dim(ydata)[2])) else as.array(ydata)
  zdata  <- if(length(dim(zdata))  != 3) array( zdata,  c(dim(zdata)[1], 1, dim(zdata)[2])) else as.array(zdata)
  budget <- if(length(dim(budget)) != 2) array(budget, c(length(budget), 1))                else as.array(budget)
  n      <- dim(xdata)[1]
  m      <- dim(xdata)[2]
  s      <- dim(ydata)[2]
  b      <- dim(zdata)[2]
  t      <- dim(xdata)[3]
  
  # Budget available at each T
  zt     <- array(budget, c(n, t))
  for(i in 2:t){zt[, i] <- zt[, i - 1] - zdata[,, (i - 1), drop = F]}
  
  # Data frames
  results.lambda       <- array(NA, dim = c(n, n, t))
  results.efficiency.s <- array(NA, dim = c(n, 1))
  results.efficiency.t <- array(NA, dim = c(n, t))
  results.xslack       <- array(NA, dim = c(n, m, t))
  results.zslack       <- array(NA, dim = c(n, b, t))
  results.yslack       <- array(NA, dim = c(n, s, t))
  results.aslack       <- array(NA, dim = c(n, 1, t))
  
  # Pointers
  p.eff <- n*t + 1
  p.xsl <- n*t + t + 1
  p.zsl <- n*t + t + m*t + 1
  p.ysl <- n*t + t + m*t + b*t + 1
  p.asl <- n*t + t + m*t + b*t + s*t + 1
  p.end <- n*t + t + m*t + b*t + s*t + b*t + 1
  
  # LP
  for(j in 1:n){
    # Declare LP
    lp.ba <- make.lp(0, n*t + t + m*t + b*t + s*t + b*t) # lambda + eff + xslack + zslack + yslack + aslack
    
    # Set objective
    if(orientation == "i") set.objfn(lp.ba, c(rep( 1/t, t)), indices = c(p.eff:(p.xsl - 1)))
    if(orientation == "o") set.objfn(lp.ba, c(rep(-1/t, t)), indices = c(p.eff:(p.xsl - 1)))
    
    # RTS
    if(rts == "vrs") for(k in 1:t){add.constraint (lp.ba, c(rep(1, n)), indices = c(((k-1)*n + 1):((k-1)*n + n)), "=", 1)}
    if(rts == "crs") for(k in 1:t){set.constr.type(lp.ba, 0, k)}
    if(rts == "irs") for(k in 1:t){add.constraint (lp.ba, c(rep(1, n)), indices = c(((k-1)*n + 1):((k-1)*n + n)), ">=", 1)}
    if(rts == "drs") for(k in 1:t){add.constraint (lp.ba, c(rep(1, n)), indices = c(((k-1)*n + 1):((k-1)*n + n)), "<=", 1)}
    
    # LP
    for(k in 1:t){
      
      # Input constraint
      for(i in 1:m){
        if(orientation == "i"){
          add.constraint(lp.ba, c(xdata[, i, k], -xdata[j, i, k], 1), indices = c(((k-1)*n + 1):((k-1)*n + n), n*t + k, p.xsl - 1 + m*(k - 1) + i), "=", 0)
        }else{
          add.constraint(lp.ba, c(xdata[, i, k], 1), indices = c(((k-1)*n + 1):((k-1)*n + n), p.xsl - 1 + m*(k - 1) + i), "=", xdata[j, i, k])
        }
      }
      
      # Budget-spent constraint
      for(i in 1:b){
        if(orientation == "i"){
          add.constraint(lp.ba, c(zdata[, i, k], -zdata[j, i, k], 1), indices = c(((k-1)*n + 1):((k-1)*n + n), n*t + k, p.zsl - 1 + b*(k - 1) + i), "=", 0)
        }else{
          add.constraint(lp.ba, c(zdata[, i, k], 1), indices = c(((k-1)*n + 1):((k-1)*n + n), p.zsl - 1 + b*(k - 1) + i), "=", zdata[j, i, k])
        }
      }
      
      # Output constraint
      for(r in 1:s){
        if(orientation == "i"){
          add.constraint(lp.ba, c(ydata[, r, k], -1), indices = c(((k-1)*n + 1):((k-1)*n + n), p.ysl - 1 + s*(k - 1) + r), "=", ydata[j, r, k])
        }else{
          add.constraint(lp.ba, c(ydata[, r, k], -ydata[j, r, k], -1), indices = c(((k-1)*n + 1):((k-1)*n + n), n*t + k, p.ysl - 1 + s*(k - 1) + r), "=", 0)
        }
      }
      
      # Budget-available constraints
      add.constraint(lp.ba, c(zt[, k], 1), indices = c(((k-1)*n + 1):((k-1)*n + n), p.asl - 1 + k), "=", zt[j, k])
    }
    
    # Bounds
    set.bounds(lp.ba, lower = c(rep(0, n*t), rep(-Inf, t), rep(0, p.end - p.xsl)))
    
    # Solve
    solve.lpExtPtr(lp.ba)
    
    # Get results
    results.efficiency.s[j]  <- abs(get.objective(lp.ba))
    temp.p                   <- get.variables(lp.ba)
    results.lambda[j,,]      <- array(temp.p[1:(n*t)], c(n, t))
    results.efficiency.t[j,] <- temp.p[p.eff:(p.xsl - 1)]
    results.xslack[j,,]      <- array(temp.p[p.xsl:(p.zsl - 1)], c(m, t))
    results.zslack[j,,]      <- array(temp.p[p.zsl:(p.ysl - 1)], c(b, t))
    results.yslack[j,,]      <- array(temp.p[p.ysl:(p.asl - 1)], c(s, t))
    results.aslack[j,,]      <- array(temp.p[p.asl:(p.end - 1)], c(1, t))
    
    # Stage II
    # Link previous solutions
    add.constraint(lp.ba, rep(1, t), indices = c(p.eff:(p.xsl - 1)), "=", results.efficiency.t[j])
    
    # slack sum max
    set.objfn(lp.ba, c(rep(-1, (p.end - p.xsl))), indices = c(p.xsl:(p.end - 1)))
    
    # solve
    solve.lpExtPtr(lp.ba)
    
    # Get results
    temp.s              <- get.variables(lp.ba)
    results.lambda[j,,] <- array(temp.s[1:(n*t)],           c(n, t))
    results.xslack[j,,] <- array(temp.s[p.xsl:(p.zsl - 1)], c(m, t))
    results.zslack[j,,] <- array(temp.s[p.zsl:(p.ysl - 1)], c(b, t))
    results.yslack[j,,] <- array(temp.s[p.ysl:(p.asl - 1)], c(s, t))
    results.aslack[j,,] <- array(temp.s[p.asl:(p.end - 1)], c(b, t))
  }
  
  # Store results
  results <- list(eff.s  = results.efficiency.s, 
                  eff.t  = results.efficiency.t, 
                  lambda = results.lambda, 
                  xslack = results.xslack, 
                  yslack = results.yslack, 
                  zslack = results.zslack,
                  aslack = results.aslack)
  return(results)
}
