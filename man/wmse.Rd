\name{wmse}
\alias{wmse}
\title{Weighted Manski's maximum score estimator}
\description{empirical minimization of averaged false positive rate and false negative rate}
\usage{wmse(mk, n1, r=1, w=2, contract=0.8, lbmdis=TRUE)}
\arguments{
  \item{mk}{biomarker values of cases followed by controls, with each row
  containing multiple markers from an individual.}
  \item{n1}{size of cases.}
  \item{r}{weight of false positive rate relative to false negative rate.}
  \item{w}{weight for l1 norm of combination coefficient in the objective
  function (w>1 guarantees sound asymptotic properties).}
  \item{contract}{reduction factor in the sequence of approximation parameters
  for indicator function.}
  \item{lbmdis}{larger biomarker value is more associated with cases if True,
  and controls otherwise.}
}
\value{
  \item{coef}{estimated combination coefficient, with unity l1 norm.}
  \item{obj}{empirical objective function: r * false positive rate + false negative rate.}
  \item{threshold}{estimated threshold.}
  \item{init_coef}{initial combination coefficient from logistic regression,
  with unity l1 norm.}
  \item{init_obj}{empirical objective function for the initial combination
  coefficient.}
  \item{init_threshold}{estimated threshold for the initial combination
  coefficient.}
}
\examples{
## simulate 3 biomarkers for 100 cases and 100 controls
mk <- rbind(matrix(rnorm(300),ncol=3),matrix(rnorm(300),ncol=3))
mk[1:100,1] <- mk[1:100,1]/sqrt(2)+1
mk[1:100,2] <- mk[1:100,2]*sqrt(2)+1

## linear combination to empirically minimize averaged false positive rate and
## false negative rate
## Require installation of 'MOSEK' to run
\dontrun{
lcom <- wmse(mk, 100)
}
}
\references{Huang and Sanda (2022). Linear biomarker combination for
constrained classification. \emph{The Annals of Statistics} 50, 2793--2815}
\author{Yijian Huang}
