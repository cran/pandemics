\name{rate2Dmiss}
\alias{rate2Dmiss}
\title{
Local Linear Estimator of the Two-Dimensional Infection (or Hospitalization) Rate from Missing-Survival-Link Data.
}
\description{
Local linear estimator of the infection or the hospitalization rate in the case of missing-survival-link data (Gámiz et al. 2024a). The rate is assumed to have two dimensions: a one-dimensional marker (typically the notification date) and time (duration). It is assumed the situation of observing aggregated data in the form of occurrences and exposures. The missing-survival link problem means that the duration is not directly observed. The estimator follows from an iterative algorithm where, at each step, full information including duration is estimated first and then, the local linear estimator is computed evaluating the function \code{hazard2D}.

}
\usage{
rate2Dmiss(t.grid, z.grid, Oi.z, Ei.z1, bs.grid,
    cv=TRUE, epsilon=1e-4, max.ite=50)
}
\arguments{
  \item{t.grid}{
a vector of \code{M} grid points for the time dimension.
}
 \item{z.grid}{
a vector of \code{M} grid points for the marker dimension.
}
 \item{Oi.z}{
a vector of length \code{M} with the number of people tested positive (infection rate) or with the number of hospitalizations (hospitalization rate) notified  ech day in \code{z.grid}.
}
 \item{Ei.z1}{
a vector of length \code{M} with the number of people tested positive   notified  each day in \code{z.grid}.
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
an integer value with the maximum number of iterations in the iterative algorithm. Default value is \code{max.ite=50}.
}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
Hazard is assumed having two dimensions: a one-dimensional marker (typically the notification date) and time (duration). It is assumed the situation of observing aggregated data in the form of occurrences and exposure, as  the function \code{hazard2D} does. The difference is that the time dimension (duration) is not directly observed. The estimator follows from an iterative algorithm where, at each step, full information, estimating the duration, is constructed first and then the local linear estimator is computed evaluating the function \code{hazard2D}. See more details on the local linear hazard estimator in the documentation of \code{\link{hazard2D}}.

To estimate the infection rate we assume that each day (z), we have information on the number of people tested positive. The vector of occurrences (\code{Oi.z})  is the observed number of people tested positive each day, removing the first \code{d} days. Here \code{d} is typically 1 (an infected person might infect one day after being tested positive), see the examples below. The vector of exposure (\code{Ei.z}) is the observed number of people tested positive each day, removing the last \code{d} days.

To estimate the hospitalization rate we assume that each day (z), we have information on  the number of new hospitalizations (\code{Oi.z}), and the number of people tested positive (\code{Ei.z})  each day.
}

\value{
\item{hi.zt}{a (\code{M} times \code{M}) matrix with the estimated rate evaluated at the grid points.}

\item{bcv}{a two-dimensional vector with the bandwidth used to compute the estimator (estimated by cross-validation if \code{cv=TRUE}).}

 \item{tol}{a numeric value with the achieved tolerance value in the algorithm.}

\item{it}{an integer with the number of iterations performed in the algorithm.}
}

\references{
Gámiz, M.L., Mammen, E., Martínez-Miranda, M.D. and Nielsen, J.P. (2024a). Low quality exposure and point processes with a view to the first phase of a pandemic. arXiv:2308.09918.

Gámiz, M.L., Mammen, E., Martínez-Miranda, M.D. and Nielsen, J.P. (2024b). Monitoring a developing pandemic with available data. arXiv:2308.09919.
}

\author{
M.L. Gámiz, E. Mammen, M.D. Martínez-Miranda and J.P. Nielsen.
}

\seealso{
\code{\link{hazard2D}}, \code{\link{forecasting}}
}

\examples{
## Analysis of the infection and hospitalization processes
## data are (cumulative) number of hospitalizations and positive tested

data('covid')
## We remove the first 56 rows (no data on testing until 13th May)
covid2<-covid[-c(1:56),]
M2<-nrow(covid2)

## 1. Rate of infection
Ei.new<-covid2$Posit
delay<-1;M2<-M2-delay
Oi.z<-Ei.new[-(1:delay)]; Ei.z1<-Ei.new[1:M2]
t.grid<-z.grid<-1:M2
bs<-t(c(5,10))
RInf<-rate2Dmiss(t.grid,z.grid,Oi.z,Ei.z1,bs.grid=bs,cv=FALSE,
       epsilon=1e-4,max.ite=50)
hi.zt<-RInf$hi.zt # the estimated infection rate

## Plot the estimated infection rate at few notification days
ddates<-covid2$Date
z1<-c(19,49,80,141)
nz<-length(z1)
zdates<-ddates[z1]
t.min<-min(M2-z1+1)
ti<-1:t.min;n0<-length(ti)

## displaced curves
alphas.I<-matrix(NA,M2,M2) # upper-triangular matrix
for(j in 1:M2) alphas.I[j,j:M2]<-hi.zt[j,1:(M2-j+1)]
alphas<-alphas.I[z1,-(1:18)]
M2a<-ncol(alphas)
yy<-c(0,max(alphas,na.rm=TRUE)+0.1)

plot(1:M2a,alphas[nz,],type='l',ylab='',xlab='Date of notification',
      main= 'Dynamic rate of infection',ylim=yy,lwd=2,lty=1,xaxt='n')
axis(1,at=z1-18,labels=c(zdates))
for (i in 2:nz) lines(1:M2a,alphas[i-1,],lwd=3,col=i,lty=i)
legend('topleft',c('Starting on date:',as.character(zdates)),
     lty=c(NA,2:nz,1),lwd=c(NA,rep(3,nz-1),2),col=c(NA,2:nz,1),bty='n')

## 2. Rate of hospitalization
Hi<-covid2$Hospi ; Hi<-Hi[-1]
Ri<-covid2$Recov ; Ri<-diff(Ri)
Di<-covid2$Death ; Di<-diff(Di)
M2<-length(Di)
## New hospitalizations are Hi-Ri-Di
newHi<-Hi[-1]-(Hi[-M2]-Ri[-M2]-Di[-M2])
newHi<-c(Hi[1],newHi)
newHi[newHi<0]<-0; # possible inconsistency in the data
Oi.z<-as.integer(newHi)
Ei.z1<-covid2$Posit
t.grid<-z.grid<-1:M2
bs<-t(c(20,10))
RHosp<-rate2Dmiss(t.grid,z.grid,Oi.z,Ei.z1,bs.grid=bs,cv=FALSE,
     epsilon=1e-4,max.ite=50)
hi.zt<-RHosp$hi.zt # the estimated rate

## Plot the estimated rate at few notification days
z1<-c(19,49,80,141)
nz<-length(z1)
zdates<-ddates[z1]
t.min<-min(M2-z1+1)
ti<-1:t.min;n0<-length(ti)
## displaced curves
alphas.I<-matrix(NA,M2,M2) # upper-triangular matrix
for(j in 1:M2) alphas.I[j,j:M2]<-hi.zt[j,1:(M2-j+1)]
alphas<-alphas.I[z1,-(1:18)]
M2a<-ncol(alphas)
yy<-c(0,max(alphas,na.rm=TRUE)+0.01)
plot(1:M2a,alphas[nz,],type='l',ylab='',xlab='Date of notification',
     main= 'Dynamic rate of hospitalization',ylim=yy,lwd=2,xaxt='n')
axis(1,at=z1-18,labels=c(zdates))
for (i in 2:nz) lines(1:M2a,alphas[i-1,],lwd=3,col=i,lty=i)
legend('topright',c('Starting on date:',as.character(zdates)),
    lty=c(NA,2:nz,1),lwd=c(NA,rep(3,nz-1),2),col=c(NA,2:nz,1),bty='n')
}


\keyword{ survival }
\keyword{ nonparametric }
