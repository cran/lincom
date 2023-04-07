eum <- function(mk, n1, s0, w=2, grdpt=10, contract=0.8, fixsens=TRUE, lbmdis=TRUE)
{if(!lbmdis) mk <- -mk
  if(!fixsens) {
    mk <- -rbind(mk[-(1:n1),],mk[1:n1,])
    n1 <- dim(mk)[1]-n1
  }
 if(s0 <= 0 || s0 >= 1) stop("Control level out of permissible bound (0,1)")

 a <- dim(mk)
 n0 <- a[1]-n1
 if(min(n1,n0) < 2) stop("# cases or controls needs to be >=2")
 p <- a[2]
 M1 <- cbind(rep(-1,n1),mk[1:n1,])
 M0 <- cbind(rep(-1,n0),mk[-(1:n1),])

 if(grdpt == 0) {
   init <- unname(glm(c(rep(1,n1),rep(0,n0)) ~ mk,family=binomial())$coef)
   init[1] <- sort(mk[1:n1,] %*% init[-1])[floor(n1*(1-s0))+1]
   init <- init/sum(abs(init[-1]))}
 else init <- grid(mk,n1,s0,grdpt)$coef
 btk <- init
 opt_coef <- btk
 opt_obj <- objcon(M1,M0,btk,s0,w)$obj
 init_sp <- 1-opt_obj

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
   dok <- colSums(M0[M0 %*% btk - epi > 0,,drop=FALSE])/(epi*n0)+c(0,w*sign(btk[-1]))
   ck <-  sum(pmax(M1 %*% btk,0))/(epi*n1)
   dck <- colSums(M1[M1 %*% btk > 0,,drop=FALSE])/(epi*n1)

   result <- ccp_core(M1,M0,n1,n0,p,s0,epi,btk,dok,ck,dck)
   obj <- objcon(M1,M0,result$coef,s0,w,epi)$obj
   if(! first) reduct <- oldobj-obj
   first <- FALSE
  }
  outobj <- objcon(M1,M0,result$coef,s0,w)$obj
  df <- abs(outobj-obj)
  if(outobj < opt_obj) {
    opt_obj <- outobj
    opt_coef <- result$coef}
  if(! fstepi) diff <- abs(outobj-oldoutobj)
  fstepi <- FALSE
 }
 opt_coef[-1] <- opt_coef[-1]/sum(abs(opt_coef[-1]))
 opt_coef[1] <- sort(mk[1:n1,] %*% opt_coef[-1])[floor(n1*(1-s0))+1]
 opt_obj <- mean(as.integer(mk[-(1:n1),] %*% opt_coef[-1] >= opt_coef[1]))
 list(coef=as.vector(opt_coef)[-1], hs=1-opt_obj, threshold=opt_coef[1],
      init_coef=init[-1], init_hs=init_sp, init_threshold=init[1])
}

grid <- function(mark,n1,sn,grdpt)
{a <- dim(mark)
 n <- a[1]
 p <- a[2]
 est <- .Fortran("gridsch",as.double(mark),as.integer(n),as.integer(p),
                 as.integer(n1),as.double(sn),as.integer(grdpt),
		 coef=double(p+1),sp=double(1),integer(p),integer(n1),
		 double(n))
 list(sp=est$sp, coef=est$coef)
}

ccp_core <- function(M1,M0,n1,n0,p,s0,epi,btk,dok,ck,dck)
{ prob <- list()

  prob$sense <- "min"

  prob$c <- c(-dok,dok[-1],rep(0,n1),rep(1/(epi*n0),n0),
              rep(0,n1+n0))

  M <- rbind(M1,M0)[,-1]
  A <- as.matrix.csr(M)
  A <- cbind(as.matrix.csr(rep(-1,n1+n0)),A,-A,-as(n1+n0,"matrix.diag.csr"),
             as(n1+n0,"matrix.diag.csr"))
  A <- rbind(A,
         as.matrix.csr(c(-dck,dck[-1],rep(1/(epi*n1),n1),
	               rep(0,n1+2*n0)),1,1+2*(p+n1+n0)),
   	 as.matrix.csr(c(0,rep(1,2*p),rep(0,2*(n1+n0))),1,1+2*(p+n1+n0)))
  prob$A <- as(A,"CsparseMatrix")

  prob$bc <- rbind(c(rep(epi,n1),rep(0,n0),-Inf,-Inf),
                   c(rep(epi,n1),rep(0,n0),-s0+ck-sum(dck*btk),1))

  prob$bx <- rbind(c(-Inf,rep(0,2*(p+n1+n0))),
                   rep(Inf,1+2*(p+n1+n0)))

  r <- Rmosek::mosek(prob,opts = list(verbose = 0))
  if (r$response$code != 0) stop(sprintf("Mosek optimization unsuccessful: %s (%d),", r$response$msg, r$response$code))
  f <- r$sol$itr$xx
  coef <- c(f[1],f[2:(p+1)]-f[(p+2):(2*p+1)])

  list(coef=coef)
}

objcon <- function(M1,M0,beta,s0,w,epi=0)
{n1 <- dim(M1)[1]
 n0 <- dim(M0)[1]
 if(epi == 0)
   {obj <- mean(as.numeric(M0 %*% beta >= 0))
    con <- s0-mean(as.numeric(M1 %*% beta >= 0))}
 else {
    obj <- mean(pmax(M0 %*% beta,0)-pmax(M0 %*% beta - epi,0))/epi+
           w*(1-sum(abs(beta[-1])))
    con <- s0-mean(pmax(M1 %*% beta,0)-pmax(M1 %*% beta - epi,0))/epi}
 list(obj=obj,con=c(con,sum(abs(beta[-1]))-1))
}
