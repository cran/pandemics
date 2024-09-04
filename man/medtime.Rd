\name{medtime}
\alias{medtime}
\title{
Median of the Time Spent in Hospital by Date of Admission. 
}
\description{
From the two-dimensional estimated hazard of death/recovery, the median of the time spent in hospital is computed  depending on the date of admission.

}
\usage{
medtime(hi.zt,z1)
}
\arguments{
  \item{hi.zt}{
a matrix with the estimated hazard of death+recovery  (\code{M} times \code{M}), obtained  evaluating the function  \code{hazard2Dmiss}.
}
 \item{z1}{
(optional) a vector of indexes between 1 and \code{M} indicating the admission days to evaluate the median. If missing then \code{z1<-c(seq(1,M-1,by=2),M-1)}.
}
}
\value{
A vector with the computed median times for each day in \code{z1}.}

\references{
Gámiz, M.L., Mammen, E., Martínez-Miranda, M.D. and Nielsen, J.P. (2024b). Monitoring a developing pandemic with available data. arXiv:2308.09919.
}
\author{
M.L. Gámiz, E. Mammen, M.D. Martínez-Miranda and J.P. Nielsen.
}

\note{
Evaluate the function \code{hazard2Dmiss} to create the estimated hazard for the argument \code{hi.zt}. 
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
hi.zt<-res.h$hi.zt # =hi2.zt+hi1.zt (two possible outcomes)

## Now the median time at few values of the marker (admission dates)
z1<-c(seq(1,M-1,by=30),M-1)
nz<-length(z1)
res<-medtime(hi.zt,z1)

plot(z1,res,ylab='days',xaxt = "n",type='p',pch=16,
  xlim=range(z1), xlab='Date of admission',
  main='Median time from admission to exit (death+recovery)')
axis(1,at=z1,labels=ddates[z1],cex=1.2)
}
}

\keyword{ survival }
\keyword{ nonparametric }
