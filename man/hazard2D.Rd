\name{hazard2D}
\alias{hazard2D}
\title{
Local Linear Estimator of the Two-Dimensional Marker-Dependent Hazard. 
}
\description{
Local linear estimator of the marker-dependent hazard introduced by Nielsen (1998). The implementation considers two dimensions only: a one-dimensional marker and time (duration). It assumes aggregated data in the form of occurrences and exposure. The bandwidth can be provided or estimated using  cross-validation (see details below).
}
\usage{
hazard2D(t.grid, z.grid, o.zt, e.zt, bs.grid, cv=FALSE)
}
\arguments{
  \item{t.grid}{
a vector of \code{M} grid points for the time dimension.
}
 \item{z.grid}{
a vector of \code{M} grid points for the marker dimension.
}
 \item{o.zt}{
a matrix of occurrences (\code{M} times \code{M}).
}
 \item{e.zt}{
a matrix of exposures (\code{M} times \code{M}).
}
 \item{bs.grid}{
a matrix with a grid of 2-dimensional bandwidths (by rows).
}
 \item{cv}{
logical, if \code{cv=TRUE} (default) bandwidth is estimated
by cross-validation.
}
}
\details{
The marker-dependent hazard local linear estimator was introduced by Nielsen (1998) and bandwidth selection provided by Gámiz et al. (2013). It assumes the Aalen multiplicative intensity model. The function implements such estimator in the case of a one-dimensional marker only (first dimension). Data are assumed to be aggregated in the form of occurrences and exposure (see Gámiz et al. 2013) in a two-dimensional grid of values (z,t).

The estimator involves on a two-dimensional kernel function and a two-dimensional bandwidth. The implemented kernel is multiplicative with \code{K(u)=(3003/2048)*(1-(u)^2)^6)*(abs(u)<1)}. The bandwidth can be provided in the argument \code{bs.grid}  in the form of a matrix (\code{2} times \code{1}). Data-driven badwidth selection  is also supported. If \code{cv=TRUE} (default) then the bandwith is estimated using cross-validation from a grid of \code{nb} two-dimensional bandwidths provided in \code{bs.grid} (\code{nb} times \code{2}).

This marker-dependent hazard local linear estimator assumes that full information is available in the form of occurrences and exposures. In a emergent or developing pandemic, this estimator will most likely be infeasible from typically available data. Thus, the estimator will  be computed after the necessary information is estimated through the function \code{hazard2Dmiss()}.  See Gámiz et al. (2022,2024a,b) for more details.
}
\value{
\item{hi.zt}{A (\code{M} times \code{M}) matrix with the estimated hazard evaluated at the grid points.}
\item{bcv}{A two-dimensional vector with the bandwidth used to compute the estimator (estimated by cross-validation if \code{cv=TRUE}).}
}
\references{
Gámiz, M.L., Janys, L., Martínez-Miranda,  M.D. and Nielsen, J.P.  (2013). Bandwidth selection in marker
dependent kernel hazard estimation. Computational Statistics & Data Analysis,  68, 155–169.

Gámiz, M.L., Mammen, E., Martínez-Miranda, M.D. and Nielsen, J.P.  (2022). Missing link survival analysis with applications to available pandemic data. Computational Statistics & Data Analysis, 169, 107405.

Gámiz, M.L., Mammen, E., Martínez-Miranda, M.D. and Nielsen, J.P. (2024a). Low quality exposure and point processes with a view to the first phase of a pandemic. arXiv:2308.09918.

Gámiz, M.L., Mammen, E., Martínez-Miranda, M.D. and Nielsen, J.P. (2024b). Monitoring a developing pandemic with available data. arXiv:2308.09919.

Nielsen, J.P. (1998). Marker dependent kernel estimation from local linear estimation.  Scandinavian Actuarial Journal, 2, 113-124.
}
\author{
M.L. Gámiz, E. Mammen, M.D. Martínez-Miranda and J.P. Nielsen.
}
\note{
Infeasible estimator to be evaluated through the function \code{hazard2Dmiss()}.
}

\seealso{
\code{\link{hazard2Dmiss}}, \code{\link{rate2Dmiss}}
}

\examples{
## 1. Define a true 2D-hazard evaluated at M*M grid points
M<-100
alpha<-dbeta((1:M)/(M+1),2,2)/(M+1)
alpha.zt<-matrix(NA,M,M)
for(z in 1:M) for (t in 1:(M-z+1)) alpha.zt[z,t]<-alpha[z]*alpha[t]

## 2. Simulate data from the true hazard (Aalen multiplicative model)
N<-10000 # sample size
set.seed(1)
##  simulate new arrivals
Ei.new<-hist(sort(runif(n=N,max=M)),breaks=0:M,plot=FALSE)$counts
##  simulate matrices of exposure (e.zt) and occurrences (o.zt)
o.zt<-e.zt<-matrix(0,M,M)
e.zt[,1]<-as.integer(Ei.new)
for(z in 1:M)  for (t in 1:(M-z+1)){
 if(e.zt[z,t]>0) o.zt[z,t]<-rbinom(1,as.integer(e.zt[z,t]),alpha.zt[z,t]) 
    else o.zt[z,t]<-0
 if(t<(M-z+1)) e.zt[z,(t+1)]<-e.zt[z,t]-o.zt[z,t]
 }

## 3. Estimate the 2d-hazard with fixed bandwidth
bs.grid<-t(c(M/2,M/2))
alpha.estim<-hazard2D(1:M,1:M,o.zt,e.zt,bs.grid,cv=FALSE)

## 4. Compare true and estimated hazards
persp(1:M,1:M,alpha.zt,main='True hazard',theta=60,xlab='z',ylab='t',zlab='')
persp(1:M,1:M,alpha.estim$hi.zt,main='Estimated hazard',theta=60,
      xlab='z',ylab='t',zlab='')

}

\keyword{ survival }
\keyword{ nonparametric }
