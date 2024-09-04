\name{forecasting}
\alias{forecasting}
\title{
Forecasting infections, hospitalizations and deaths in hospital.}
\description{
Forecasting infections, hospitalizations and deaths in hospital. From estimated two-dimensional rates of infections and hospilizations and estimated hazards of deaths and recoveries, new infections, hospitalizations and deaths are predicted in a forecasting period.
}
\usage{
forecasting(Cval=1, period, RoInf, RoHosp, RoDeath, RoRec,
            Pz, newHz, Hz, Dz, Rz)

}
\arguments{
  \item{Cval}{
a numeric value with the infection indicator. Default value is \code{Cval=1}.
}
 \item{period}{
an integer value with the number of days to forecast.
}
 \item{RoInf}{
a matrix with the estimated rate of infection (\code{M-1} times \code{M-1}) estimated from \code{M} days, evaluating the function  \code{rate2Dmiss}.
}
 \item{RoHosp}{ (optional)
a matrix with the estimated rate of hospitalization (\code{M} times \code{M}), obtained evaluating the function  \code{rate2Dmiss}.
}
 \item{RoDeath}{ (optional)
a matrix with the estimated hazard of deaths  (\code{M} times \code{M}), obtained  evaluating the function  \code{hazard2Dmiss}.
}
 \item{RoRec}{ (optional)
a matrix with the estimated hazard of recoveries  (\code{M} times \code{M}), obtained  evaluating the function  \code{hazard2Dmiss}.
}
 \item{Pz}{
a vector with the observed (past) number of people tested positive each day (length \code{M}).  }
\item{newHz}{
a vector with the observed number of new hospitalizations each day (length \code{M}).  }
\item{Hz}{
a vector with the observed total number people in   hospitals each day (length \code{M}).  }
\item{Dz}{
a vector with the observed number of deaths each day (length \code{M}).  }
\item{Rz}{
a vector with the observed number of recoveries each day (length \code{M}).  }
}

\details{
To create estimated rates and hazards for the arguments evaluate the functions \code{hazard2Dmiss} and \code{rate2Dmiss}. If the argument \code{RoHosp} is missing then only infections are predicted. If any of the argument \code{RoDeath} or \code{RoRec} is missing  then deaths (recoveries) are not predicted.

The infection indicator \code{Cval} is typically provided by experts. At the more recent observation (t), it indicates whether the future (t+h, h > 0) will be different or equal to the immediate past. If \code{Cval=1} then we can forecast the immediate future based on the immediate past. There might be periods where little is happening (and the indicator might have a tendency
to increase slowly) and there might be few but very important change-points where measures (e.g. lock down) are introduced to minimize future infections (and the indicator might drop dramatically in a matter of days). See Gámiz et al. (2024a) for more details and the relationship to the well-known reproduction number.
}

\value{
\item{Pz.fitted}{a vector with the fitted values for the number of infections in the past. }
\item{Pz.pred}{a vector with the predicted values for the number of infections in the forecasting period. }
\item{newHz.fitted}{a vector with the fitted values for the number of hospitalizations in the past. }
\item{newHz.pred}{a vector with the predicted values for the number of hospitalizations in the forecasting period. }
\item{Dz.fitted}{a vector with the fitted values for the number of deaths in the past. }
\item{Dz.pred}{a vector with the predicted values for the number of deaths in the forecasting period. }
\item{Rz.fitted}{a vector with the fitted values for the number of recoveries in the past. }
\item{Rz.pred}{a vector with the predicted values for the number of recoveries in the forecasting period. }
}


\references{
Gámiz, M.L., Mammen, E., Martínez-Miranda, M.D. and Nielsen, J.P. (2024a). Low quality exposure and point processes with a view to the first phase of a pandemic. arXiv:2308.09918.

Gámiz, M.L., Mammen, E., Martínez-Miranda, M.D. and Nielsen, J.P. (2024b). Monitoring a developing pandemic with available data. arXiv:2308.09919.
}

\author{
M.L. Gámiz, E. Mammen, M.D. Martínez-Miranda and J.P. Nielsen.
}

\seealso{
\code{\link{hazard2Dmiss}}, \code{\link{rate2Dmiss}}
}

\examples{
## Forecasting October 2020 from past data
data('covid')
## We remove the first 56 rows(no data on testing until 13th May)
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

## 1. First we estimate the death and recovery hazards,
## and the infection and hospitalization rates.
## Estimation using data up to 30-Sep-2020
Ms<-141   # up to 30-Sep

## 1.1. Infection rate
Ei.new<-Pi[1:Ms]
delay<-1;Msd<-Ms-delay
Oi.z<-Ei.new[-(1:delay)]
Ei.z1<-Ei.new[1:Msd];
t.grid<-z.grid<-1:Msd
bs<-t(c(5,10))
RInf<-rate2Dmiss(t.grid,z.grid,Oi.z,Ei.z1,
     bs.grid=bs,cv=FALSE)
RoInf<-RInf$hi.zt

## 1.2. Hospitalization rate
Oi.z<-newHi[1:Ms]
# The exposure are the positives data
Ei.z1<-Pi[1:Ms]
t.grid<-z.grid<-1:Ms
bs<-t(c(150,150))
RHosp<-rate2Dmiss(t.grid,z.grid,Oi.z,Ei.z1,
     bs.grid=bs,cv=FALSE)
RoHosp<-matrix(as.numeric(RHosp$hi.zt),Ms,Ms)

## 1.3. Hazards of deaths and recoveries
Oi1.z<-Di[1:Ms]  # deaths
Oi2.z<-Ri[1:Ms]   # recoveries
Ei.z<-Hi[1:Ms]    # exposure is left as cumulative
t.grid<-z.grid<-1:Ms
bs<-t(c(150,40))
res.h<-hazard2Dmiss(t.grid,z.grid,Oi1.z,Oi2.z,Ei.z,
     bs.grid=bs,cv=FALSE)
RoDeath<-as.matrix(res.h$hi1.zt,Ms,Ms) ## 2D-hazard of deaths
RoRec<-as.matrix(res.h$hi2.zt,Ms,Ms) ## 2D-hazard of recoveries

## 2. Forecasting with a given infection indicator
Cval<-1.5
period<-32  # forecasts up to 1st November, 32 days
fore<-forecasting(Cval,period,RoInf,RoHosp,RoDeath,RoRec,
    Pi[1:Ms],newHi[1:Ms],Hi[1:Ms],Di[1:Ms],Ri[1:Ms])
Hz.pred<-fore$newHz.pred
Pz.pred<-fore$Pz.pred
Dz.pred<-fore$Dz.pred

## 3. Plot forecasts and compare with observed values
## (future values are shown for predictions validation)
Pz.obs<-Pi
Hz.obs<-newHi
Dz.obs<-Di

## Graph with the new infections up to 30-Sep and forecasts
yy<-range(Pz.pred,Pz.obs,na.rm=TRUE)
plot(1:(Ms-1),Pz.obs[3:(Ms+1)],ylab='',xlab='Date of notification',
     main='Forecasts of new positives in October 2020',pch=19,
     ylim=yy,xaxt='n',xlim=c(1,Ms+1+period))
oat<-c(1,17,32,47,62,78,93,109,124,139,154,170)
olab<-ddates[oat]
axis(1,at=oat,labels=olab)
## forecasts start on 30-sep but we only plot from 1-October
points(Ms:(Ms+period-2),Pz.obs[(Ms+2):(Ms+period)],col=1,pch=1)
lines(Ms:(Ms+period-2),Pz.pred[-1],col=2,lty=2,lwd=3)
legend('topleft',c('Data: daily number of new positives until 30-Sep',
    paste('Forecasts with Cval=',round(Cval,2)),
    'True numbers of new positives in October'),
    pch=c(19,NA,1),col=c(1,2,1),lty=c(NA,2,NA),lwd=c(NA,3,NA),bty='n')

## Graph with the new hospitalizations up to 30-Sep and forecasts
ylim1<-range(Hz.pred,Hz.obs[-(1:2)],na.rm=TRUE)
plot(1:(Ms-1),Hz.obs[3:(Ms+1)],ylab='',xlab='Date of admission',
     main='Forecasts of new hospitalizations in October 2020',
     pch=19,
     ylim=ylim1,xaxt='n',xlim=c(1,Ms+1+period))
oat<-c(1,17,32,47,62,78,93,109,124,139,154,170)
olab<-ddates[oat]
axis(1,at=oat,labels=olab)
points(Ms:(Ms+period-2),Hz.obs[(Ms+2):(Ms+period)])
lines(Ms:(Ms+period-2),Hz.pred[-1],col=2,lty=2,lwd=2)
legend('topleft',c('Data: daily number of new hospitalizations until 30-Sep',
    paste('Forecasts with Cval=',round(Cval,2)),
    'True numbers of new hospitalizations in October'),
    pch=c(19,NA,1),col=c(1,2,1),lty=c(NA,2,NA),lwd=c(NA,3,NA),bty='n')


## Graph with deaths up to 30-Sep and forecasts
ylim1<-range(Dz.pred,Dz.obs[-(1:2)],na.rm=TRUE)
plot(1:(Ms-1),Dz.obs[3:(Ms+1)],ylab='',xlab='Date of admission',
     main='Forecasts of deaths in October 2020',
     pch=19,ylim=ylim1,xaxt='n',xlim=c(1,Ms+1+period))
oat<-c(1,17,32,47,62,78,93,109,124,139,154,170)
olab<-ddates[oat]
axis(1,at=oat,labels=olab)
points(Ms:(Ms+period-2),Dz.obs[(Ms+2):(Ms+period)])
lines(Ms:(Ms+period-2),Dz.pred[-1],col=2,lty=2,lwd=2)
legend('topleft',c('Data: daily number of deaths until 30-Sep',
    paste('Forecasts with Cval=',round(Cval,2)),
    'True numbers of deaths in October'),
    pch=c(19,NA,1),col=c(1,2,1),lty=c(NA,2,NA),
    lwd=c(NA,3,NA),bty='n')

}

\keyword{ survival }
\keyword{ nonparametric }
