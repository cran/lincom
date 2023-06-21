wmse <- function(mk, n1, r=1, w=2, contract=0.8, lbmdis=TRUE)
{if(!lbmdis) mk <- -mk
 if(r <= 0) stop("Relative cost of FP to FN needs to be positive")

 a <- dim(mk)
 n0 <- a[1]-n1
 p <- a[2]
 M1 <- cbind(rep(-1,n1),mk[1:n1,])
 M0 <- cbind(rep(-1,n0),mk[-(1:n1),])

 init <- unname(glm(c(rep(1,n1),rep(0,n0)) ~ mk,family=binomial())$coef)
 init <- init/sum(abs(init[-1]))
 init[1] <- -init[1]
 btk <- init
 opt_coef <- btk
 opt_obj <- wmse_objcon(M1,M0,btk,r,w)$obj
 init_obj <- opt_obj

 lc <- mk %*% btk[-1]
 a <- sort(lc[1:n1])
 epi <- max(a[-1]-a[-n1])
 b <- sort(lc[-(1:n1)])
 epi <- max(epi,b[-1]-b[-n0])

 btk[1] <- btk[1]-epi
 fstepi <- TRUE
 diff <- 1
 df <- 1

 while(fstepi || df > 1e-3 || diff > 1e-3)
 {if(! fstepi) {
    btk <- result$coef
    oldoutobj <- outobj
    epi <- epi*contract}
    
  first <- TRUE
  reduct <- 1
  while(first || reduct > 1e-3)
  {if(! first) {
     btk <- result$coef
     oldobj <- obj}
   dok <- r*colSums(M0[M0 %*% btk - epi > 0,,drop=FALSE])/(epi*n0)+
          colSums(M1[M1 %*% btk > 0,,drop=FALSE])/(epi*n1)+
	  c(0,w*(1+r)*sign(btk[-1]))

   result <- wmse_core(M1,M0,n1,n0,p,r,epi,btk,dok)
   obj <- wmse_objcon(M1,M0,result$coef,r,w,epi)$obj
   if(! first) reduct <- oldobj-obj
   first <- FALSE
  }
  outobj <- wmse_objcon(M1,M0,result$coef,r,w)$obj
  df <- abs(outobj-obj)
  if(outobj < opt_obj) {
    opt_obj <- outobj
    opt_coef <- result$coef}
  if(! fstepi) diff <- abs(outobj-oldoutobj)
  fstepi <- FALSE
 }
 opt_coef <- opt_coef/sum(abs(opt_coef[-1]))
 opt_obj <- wmse_objcon(M1,M0,opt_coef,r,w)$obj
 list(coef=as.vector(opt_coef)[-1], obj=opt_obj, threshold=opt_coef[1],
      init_coef=init[-1], init_obj=init_obj+1, init_threshold=init[1])
}

wmse_core <- function(M1,M0,n1,n0,p,r,epi,btk,dok)
{ prob <- list()

  prob$sense <- "min"

  prob$c <- c(-dok,dok[-1],rep(1/(epi*n1),n1),rep(r/(epi*n0),n0),
              rep(0,n1+n0))

  M <- rbind(M1,M0)[,-1]
  A <- as.matrix.csr(M)
  A <- cbind(as.matrix.csr(rep(-1,n1+n0)),A,-A,-as(n1+n0,"matrix.diag.csr"),
             as(n1+n0,"matrix.diag.csr"))
  A <- rbind(A,
   	 as.matrix.csr(c(0,rep(1,2*p),rep(0,2*(n1+n0))),1,1+2*(p+n1+n0)))
  prob$A <- as(A,"CsparseMatrix")

  prob$bc <- rbind(c(rep(epi,n1),rep(0,n0),-Inf),
                   c(rep(epi,n1),rep(0,n0),1))

  prob$bx <- rbind(c(-Inf,rep(0,2*(p+n1+n0))),
                   rep(Inf,1+2*(p+n1+n0)))

  rslt <- Rmosek::mosek(prob,opts = list(verbose = 0))
  if (rslt$response$code != 0) stop(sprintf("Mosek optimization unsuccessful: %s (%d),", rslt$response$msg, rslt$response$code))
  f <- rslt$sol$itr$xx
  coef <- c(f[1],f[2:(p+1)]-f[(p+2):(2*p+1)])

  list(coef=coef)
}

wmse_objcon <- function(M1,M0,beta,r,w,epi=0)
{n1 <- dim(M1)[1]
 n0 <- dim(M0)[1]
 obj <- ifelse(epi == 0,
               r*mean(as.numeric(M0 %*% beta >= 0))+
               mean(as.numeric(M1 %*% beta < 0)),
               r*mean(pmax(M0 %*% beta,0)-pmax(M0 %*% beta - epi,0))/epi
	       +1-mean(pmax(M1 %*% beta,0)-pmax(M1 %*% beta - epi,0))/epi+
               -w*(1+r)*sum(abs(beta[-1])))
 list(obj=obj,con=sum(abs(beta[-1]))-1)
}
