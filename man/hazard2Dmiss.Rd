\name{hazard2Dmiss}
\alias{hazard2Dmiss}
\title{
Local Linear Estimator of the Two-Dimensional Marker-Dependent Hazard from Missing-Survival-Link Data. 
}
\description{
Local linear estimator of the marker-dependent hazard in the case of missing-survival-link data. Hazard is assumed to have two dimensions: a one-dimensional marker (typically the notification date) and time (duration). It is assumed the situation of observing aggregated data in the form of occurrences and exposures. The missing-survival link problem means that the duration is not directly observed. The estimator follows from an iterative algorithm where, at each step, full information including duration is estimated first and then, the local linear estimator is computed evaluating the function \code{hazard2D}.

}
\usage{
hazard2Dmiss(t.grid, z.grid, Oi1.z, Oi2.z, Ei.z, bs.grid,
    cv=TRUE, epsilon=1e-4, max.ite=50)
}
\arguments{
  \item{t.grid}{
a vector of \code{M} grid points for the time dimension.
}
 \item{z.grid}{
a vector of \code{M} grid points for the marker dimension.
}
 \item{Oi1.z}{
a vector of length \code{M} with the number of deaths notified each day in \code{z.grid}.}
 \item{Oi2.z}{
a vector of length \code{M} with the number of recoveries notified each day in \code{z.grid}.}
 \item{Ei.z}{
a vector of length \code{M} with the total number of people in the hospital each day in \code{z.grid}.. %, arrivals at any day from 1:z).
}
 \item{bs.grid}{
a matrix with a grid of 2-dimensional bandwidths (by rows).
}
 \item{cv}{
logical, if \code{cv=TRUE} (default) bandwidth is estimated by cross-validation.
}
 \item{epsilon}{ 
a numeric value with the tolerance in the iterative algorithm. Default value is \code{epsilon=1e-4}.}
 \item{max.ite}{ 
a integer value with the maximum number of iterations in the iterative algorithm. Default value is \code{max.ite=50}.
}
}

\details{
Hazard is assumed having two dimensions: a one-dimensional marker  and time. It is assumed the situation of observing aggregated data in the form of occurrences and exposure, as  the function \code{hazard2D} does. The difference is that the time  dimension (duration) is not observed. The estimator follows from an iterative algorithm where, at each step, full information, estimating the duration, is constructed first and then the local linear estimator is computed by evaluating the function \code{hazard2D}. 

Missing-link-survival data means that duration is not observed directly. Each day (z), we get information on number of people remaining in hospital (exposure, Ei.z), the number of deaths (Oi1.z) and the number of recoveries (Oi2.z) on that day.
}

\value{
\item{hi.zt}{a (\code{M} times \code{M}) matrix with the estimated hazard (death+recovery) evaluated at the grid points.}
\item{hi1.zt}{a (\code{M} times \code{M}) matrix with the estimated hazard of deaths evaluated at the grid points.}
\item{hi2.zt}{a (\code{M} times \code{M}) matrix with the estimated hazard of recoveries evaluated at the grid points.}
\item{bcv}{a two-dimensional vector with the bandwidth used to compute the estimator (estimated by cross-validation if \code{cv=TRUE}).}
 \item{tol}{a numeric value with the achieved tolerance value in the algorithm.}
 \item{it}{an integer with the number of iterations performed in the algorithm.}
}
\references{
Gámiz, M.L., Mammen, E., Martínez-Miranda, M.D. and Nielsen, J.P.  (2022). Missing link survival analysis with applications to available pandemic data. Computational Statistics & Data Analysis, 169, 107405.

Gámiz, M.L., Mammen, E., Martínez-Miranda, M.D. and Nielsen, J.P. (2024a). Low quality exposure and point processes with a view to the first phase of a pandemic. arXiv:2308.09918.

Gámiz, M.L., Mammen, E., Martínez-Miranda, M.D. and Nielsen, J.P. (2024b). Monitoring a developing pandemic with available data. arXiv:2308.09919.
}

\author{
M.L. Gámiz, E. Mammen, M.D. Martínez-Miranda and J.P. Nielsen.
}

\seealso{
\code{\link{hazard2D}}, \code{\link{medtime}}, \code{\link{poutcome}}
}

\examples{
\donttest{
## Survival analysis of duration in covid-19 data
data('covid')
Ei.z<-covid$Hospi   # exposure for survival analysis
Oi1.z<-covid$Death  # deaths
Oi2.z<-covid$Recov  # recoveries
## compute incremental values
Oi1.z<-diff(Oi1.z)
Oi2.z<-diff(Oi2.z)
Ei.z<-Ei.z[-1]     # exposure is cumulative
M<-length(Ei.z)
t.grid<-z.grid<-1:M
## notification date (marker)
ddates<-covid$Date

## Hazard estimate with a fixed bandwidth
bs<-t(c(150,150))
res.h<-hazard2Dmiss(t.grid,z.grid,Oi1.z,Oi2.z,Ei.z,bs.grid=bs,cv=FALSE)
hi1.zt<-res.h$hi1.zt ## 2D-hazard of recoveries
hi2.zt<-res.h$hi2.zt ## 2D-hazard of deaths

## Plot of hazard of deaths on several dates
z1<-c(13,44,105,197)
zdates<-ddates[z1] ; nz<-length(z1)
t.min<-35  # consider a maximum duration of 35 days for the plots
ti<-1:t.min ; n0<-length(ti)

yy<-range(hi1.zt[z1,1:n0],na.rm=TRUE)
yy[1]<-yy[1]-.0003;yy[2]<-yy[2]+.03
plot(ti,hi1.zt[nz,1:n0], main='Deaths',type='l',
  xlab='Time (days) from admission',ylab='Hazard',ylim=yy)
for(i in 2:nz) lines(ti,hi1.zt[z1[i-1],1:n0],lwd=3,col=i,lty=i)
legend('topright',c('Date of admision', as.character(zdates)),
  lty=c(NA,2:nz,1),lwd=c(NA,rep(3,nz-1),2),col=c(NA,2:nz,1),bty='n')

## Same for recoveries 
yy<-range(hi2.zt[z1,1:n0],na.rm=TRUE)
yy[2]<-yy[2]+.05
plot(ti,hi2.zt[nz,1:n0], main='Recoveries',type='l',
    xlab='Time (days) from admission',ylab='Hazard',ylim=yy)
for(i in 2:nz) lines(ti,hi2.zt[z1[i-1],1:n0],lwd=3,col=i,lty=i)
    legend('topright',c('Date of admision',as.character(zdates)),
    lty=c(NA,2:nz,1),lwd=c(NA,rep(3,nz-1),2),col=c(NA,2:nz,1),bty='n')
}
}


\keyword{ survival }
\keyword{ nonparametric }
