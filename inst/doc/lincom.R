## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----install, eval=FALSE, message=FALSE, warning=FALSE------------------------
#  install.packages("lincom")

## ----simulation, eval=TRUE, message=FALSE, warning=FALSE----------------------
library(lincom)
## simulate 3 biomarkers for 100 cases and 100 controls
n1 <- 100
n0 <- 100
mk <- rbind(matrix(rnorm(3*n1),ncol=3),matrix(rnorm(3*n0),ncol=3))
mk[1:n1,1] <- mk[1:n1,1]/sqrt(2)+1
mk[1:n1,2] <- mk[1:n1,2]*sqrt(2)+1

## ----eum, message=FALSE, warning=FALSE----------------------------------------
## The following two lines are commented out - require installation of 'MOSEK' to run
#lcom1 <- eum(mk, n1=n1, s0=0.95, grdpt=0)
#lcom1

## ----wmse, message=FALSE, warning=FALSE---------------------------------------
## default relative weight r=1.
## Require installation of 'MOSEK' to run
## The following two lines are commented out - require installation of 'MOSEK' to run
#lcom2 <- wmse(mk, n1=n1)
#lcom2

