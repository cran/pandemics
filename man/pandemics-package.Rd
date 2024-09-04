\name{pandemics-package}
\alias{pandemics-package}
\alias{pandemics}
\docType{package}
\title{
Monitoring a Developing Pandemic with Available Data.
}
\description{
Full dynamic system to describe and forecast the spread
and the severity of a developing pandemic based on typically available data. The system involves three different transitions: infection-infection, infection-hospitalization and hospitalization-death/recovery. The intensities of these transitions are dynamic and estimated using nonparametric local linear estimators.

Estimation is performed from aggregated data consisting of number of infections (positive tests), hospitalizations, deaths and recoveries notified each day.

The package can be used to provide forecasts of new infections, hospitalizations and deaths, as well as typical survival indicators such as the median time from admission in hospital to recovery or death, or the probability that a patient who has been in hospital for a number of days can leave it alive.

}
\details{
\tabular{ll}{
Package: \tab pandemics\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2024-30-08\cr
License: \tab GPL-2\cr
}
}
\author{
M.L. Gamiz, E. Mammen, M.D. Martinez-Miranda and J.P. Nielsen

Maintainer: Maria Dolores Martinez-Miranda <mmiranda@ugr.es>
}
\references{
Gámiz, M.L., Mammen, E., Martínez-Miranda, M.D. and Nielsen, J.P.  (2022). Missing link survival analysis with applications to available pandemic data. Computational Statistics & Data Analysis, 169, 107405.

Gámiz, M.L., Mammen, E., Martínez-Miranda, M.D. and Nielsen, J.P. (2024a). Low quality exposure and point processes with a view to the first phase of a pandemic. arXiv:2308.09918.

Gámiz, M.L., Mammen, E., Martínez-Miranda, M.D. and Nielsen, J.P. (2024b). Monitoring a developing pandemic with available data. arXiv:2308.09919.
}
\keyword{ package }
\keyword{ nonparametric }
\keyword{ survival }

\examples{
## Forecasting new infections in October 2020 from data up to 30-Sep
data('covid')
## We remove the first 56 rows  (no data on testing until 13th May)
covid2<-covid[-c(1:56),]
M2<-nrow(covid2)
Hi<-covid2$Hospi
Hi<-Hi[-1]
Ri<-covid2$Recov
Ri<-diff(Ri)
Di<-covid2$Death
Di<-diff(Di)
Pi<-covid2$Posit
M2<-length(Di)
# New hospitalizations are Hi-Ri-Di
newHi<-Hi[-1]-(Hi[-M2]-Ri[-M2]-Di[-M2])
newHi<-as.integer(c(Hi[1],newHi))
newHi[newHi<0]<-0
ddates<-covid2$Date

## 1. First estimate the infection rate
Ms<-141   # up to 30-Sep
Ei.new<-Pi[1:Ms]
delay<-1;Msd<-Ms-delay
Oi.z<-Ei.new[-(1:delay)]
Ei.z1<-Ei.new[1:Msd];
t.grid<-z.grid<-1:Msd
bs<-t(c(5,10))
RInf<-rate2Dmiss(t.grid,z.grid,Oi.z,Ei.z1,bs.grid=bs,cv=FALSE)
RoInf<-RInf$hi.zt

## 2. Forecasting now with a given infection indicator
Cval<-1.5
period<-32  # forecasts up to 1st November, 32 days
fore<-forecasting(Cval,period,RoInf,Pz=Pi[1:Ms],
    newHz=newHi[1:Ms],Hz=Hi[1:Ms],Dz=Di[1:Ms],Rz=Ri[1:Ms])
Pz.pred<-fore$Pz.pred

## 3. Plot forecasts and compare with observed values
## (future values are shown for predictions validation)
Pz.obs<-Pi
yy<-range(Pz.pred,Pz.obs,na.rm=TRUE)
plot(1:(Ms-1),Pz.obs[3:(Ms+1)],ylab='',xlab='Date of notification',
     main='Forecasts of new positives in October 2020',pch=19,
     ylim=yy,xaxt='n',xlim=c(1,Ms+1+period))
oat<-c(1,17,32,47,62,78,93,109,124,139,154,170)
olab<-ddates[oat]
axis(1,at=oat,labels=olab)
points(Ms:(Ms+period-2),Pz.obs[(Ms+2):(Ms+period)],col=1,pch=1)
lines(Ms:(Ms+period-2),Pz.pred[-1],col=2,lty=2,lwd=3)
legend('topleft',c('Data: daily number of new positives until 30-Sep',
    paste('Forecasts with Cval=',round(Cval,2)),
    'True numbers of new positives in October'),
    pch=c(19,NA,1),col=c(1,2,1),lty=c(NA,2,NA),lwd=c(NA,3,NA),bty='n')

}
