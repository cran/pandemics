\name{poutcome}
\alias{poutcome}
\title{
Probability of Outcome by Cause Specific. 
}
\description{
From the two-dimensional estimated hazards of deaths and recoveries, the probability that a person, who has been in hospital for a number of days, leaves the hospital  alive or death, depending on the date of admission.

}
\usage{
poutcome(hi1.zt,hi2.zt,z1)
}
\arguments{
  \item{hi1.zt}{
a matrix with the estimated hazard of deaths  (\code{M} times \code{M}). 
}
  \item{hi2.zt}{
a matrix with the estimated hazard of recoveries  (\code{M} times \code{M}). 
}
 \item{z1}{
(optional) a vector of indexes between 1 and \code{M} indicating the admission days to evaluate the probabilities. If missing then \code{z1<-c(seq(1,M-1,by=2),M-1)}.
}
}

\value{
\item{alive.zt}{ a matrix (\code{M} times \code{M}) with the computed probabilities of leaving the hospital alive (each column corresponds to a day in \code{z1}).}
\item{death.zt}{ a matrix (\code{M} times \code{M}) with the computed probabilities of dying in hospital (each column corresponds to a day in \code{z1}).}
}

\references{
Gámiz, M.L., Mammen, E., Martínez-Miranda, M.D. and Nielsen, J.P. (2024b). Monitoring a developing pandemic with available data. arXiv:2308.09919.
}
\author{
M.L. Gámiz, E. Mammen, M.D. Martínez-Miranda and J.P. Nielsen.
}

\note{
Evaluate the function \code{hazard2Dmiss} to create the estimated hazards for the arguments \code{hi1.zt} and \code{hi2.zt}. 
}

\seealso{
\code{\link{hazard2Dmiss}}
}

\examples{
\donttest{
data('covid')
Ei.z<-covid$Hospi   # exposure for survival analysis
Oi1.z<-covid$Death  # deaths
Oi2.z<-covid$Recov  # recoveries
# compute incremental values
Oi1.z<-diff(Oi1.z)
Oi2.z<-diff(Oi2.z)
Ei.z<-Ei.z[-1]     # exposure is left as cumulative
M<-length(Ei.z)
t.grid<-z.grid<-1:M
# notification date (marker)
ddates<-covid$Date

## First compute the estimated hazard
bs<-t(c(150,150))
res.h<-hazard2Dmiss(t.grid,z.grid,Oi1.z,Oi2.z,Ei.z,bs.grid=bs,cv=FALSE)
hi1.zt<-res.h$hi1.zt # 2D-hazard of deaths
hi2.zt<-res.h$hi2.zt # 2D-hazard of recoveries

## Now the probabilities at few values of the marker (admission dates)
z1<-c(13,44,105,197)
zdates<-ddates[z1]
nz<-length(z1)
t.min<-35  # maximum duration is 35
ti<-1:t.min;n0<-length(ti)
res<-poutcome(hi1.zt,hi2.zt,z1)
alive.zt<-res$alive.zt
death.zt<-res$death.zt

# Cause= recovery
plot(ti,alive.zt[1:n0,1],ylim=c(0,1),lwd=2,type='l',
     main='Probability to get out alive',
     ylab='',xlab='Time from admission (days)')
for(i in 2:nz) lines(ti,alive.zt[1:n0,i-1],lwd=3,col=i,lty=i)
legend('bottom',legend=zdates,lty=c(2:nz,1),
  lwd=c(rep(3,nz-1),2),col=c(2:nz,1),bty='n')

# Cause= death
plot(ti,death.zt[1:n0,1],ylim=c(0,1),lwd=2,type='l',
     main='Probability of dying in the hospital',
     ylab='',xlab='Time from admission (days)')
for(i in 2:nz) lines(ti,death.zt[1:n0,i-1],lwd=3,col=i,lty=i)
legend('top',legend=zdates,lty=c(2:nz,1),lwd=c(rep(3,nz-1),2),
       col=c(2:nz,1),bty='n')
}
}

\keyword{ survival }
\keyword{ nonparametric }
