\name{eum}
\alias{eum}
\title{Empirical performance optimization for specificity (or sensitivity) at
  a controlled sensitivity (or specificity) level}
\description{Linear combination of multiple biomarkers}
\usage{eum(mk, n1, s0, w=2, grdpt=10, contract=0.8, fixsens=TRUE, lbmdis=TRUE)}
\arguments{
  \item{mk}{biomarker values of cases followed by controls, with each row
  containing multiple markers from an individual.}
  \item{n1}{size of cases.}
  \item{s0}{controlled level of sensitivity or specificity.}
  \item{w}{weight for l1 norm of combination coefficient in the objective
  function (w>1 guarantees sound asymptotic properties).}
  \item{grdpt}{number of grid points in coarse grid search for initial value;
  if grdpt=0, use logistic regression instead.}
  \item{contract}{reduction factor in the sequence of approximation parameters
  for indicator function.}
  \item{fixsens}{fixing sensitivity if True, and specificity otherwise.}
  \item{lbmdis}{larger biomarker value is more associated with cases if True,
  and controls otherwise.}
}
\value{
  \item{coef}{estimated combination coefficient, with unity l1 norm.}
  \item{hs}{empirical estimate of specificity at controlled
  sensitivity, or vice versa.}
  \item{threshold}{estimated threshold.}
  \item{init_coef}{initial combination coefficient, with unity l1 norm.}
  \item{init_hs}{initial specificity at controlled sensitivity, or vice versa.}
  \item{init_threshold}{estimated threshold for the initial combination
  coefficient.}
}
\examples{
## simulate 3 biomarkers for 100 cases and 100 controls
mk <- rbind(matrix(rnorm(300),ncol=3),matrix(rnorm(300),ncol=3))
mk[1:100,1] <- mk[1:100,1]/sqrt(2)+1
mk[1:100,2] <- mk[1:100,2]*sqrt(2)+1

## linear combination to empirically maximize specificity at controlled 0.95
## sensitivity
## Require installation of 'MOSEK' to run
\dontrun{
lcom <- eum(mk, 100, 0.95, grdpt=0)
}
}
\references{Huang and Sanda (2022). Linear biomarker combination for
constrained classification. \emph{The Annals of Statistics} 50, 2793--2815}
\author{Yijian Huang}
