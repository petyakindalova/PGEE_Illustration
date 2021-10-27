# Updated functions from package 'binarySimCLF_1.0.tar.gz'
# code was originally written for logit link - we use log link
# 
gen.dataPP <- function(beta, nc, cl.size=2:6, p, rho){
  cl.size <- round(runif(nc, cl.size[1], cl.size[length(cl.size)]))
  N <- sum(cl.size)
  X1i <- rep(rbinom(nc, 1, p), times=cl.size)
  j <- cl.size
  if(length(unique(cl.size)) == 1) tij <- c(sapply(j, function(x) 0.2*(1:x)))
  if(length(unique(cl.size)) > 1) tij <- unlist(sapply(j, function(x) 0.2*(1:x)))
  intercept <- rep(1, N)
  dat <- cbind(intercept, X1i, obstime=tij)
  a <- exp(apply(dat, 1 , function(r) r%*% beta))
  #pi <- a/(1+a)
  pi <- a #changed by PK
  id <- rep(1:nc, times=j)
  dd <- data.frame(id, pi)
  d.pi <- split(dd, id)
  R <- lapply(cl.size, function(x) xch(x, rho))
  V <- list()
  for(i in 1:nc){
    #print(d.pi[[i]]$pi)
    #print(R[[i]])
    V[[i]] <- cor2var(R[[i]], d.pi[[i]]$pi)
    #print(V[[i]])
  }
  B <- lapply(V, function(Vi) allReg(Vi));
  # Checks CLF compatibility.
  clf.compat <- NULL
  for(i in 1:nc){
    clf.compat[i] <- blrchk(d.pi[[i]][,2], V[[i]])
  }
  y <- list()
  for(i in 1:nc){
    if(clf.compat[i]){
      y[[i]] = mbsclf(1, d.pi[[i]][,2], B[[i]])$y
    }
    if(clf.compat[i]==FALSE){
      y[[i]] = rep(NA, cl.size[i])
    }
  }
  #print(y)
  yij <- unlist(y)
  data <- data.frame(id=id, yij= yij, intercept= intercept,
                     X1i=X1i, obstime=tij)
  return(data)
}

cor2var <- function(r, mu){
  p <- dim(r)[1];
  #print(r)
  d <- sqrt(diag(mu));
  #print(d)
  V <- d %*% r %*% d;
  #print(V)
  rownames(V) <- NULL;
  return(V);
}

