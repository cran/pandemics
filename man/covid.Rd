\name{covid}
\alias{covid}
\docType{data}
\title{
COVID-19 Outbreak Data in France
}
\description{
Daily data of COVID-19 cases (positive-tested),  hospitalizations, deaths and recovered published in the open platform for French public data during the period 2020-03-18 to 2020-11-01.
}
\usage{data("covid")}
\format{
  Data frame with 229 observations and 5 variables:

[,1] \code{Date} Notification date

[,2] \code{Hospi} Daily total number of persons in hospital

[,3] \code{Death} Daily total number of deaths (in hospital)

[,4] \code{Recov} Daily total number of persons discharged from hospital

[,5] \code{Posit} Daily total number of persons tested positive
}

\source{
Data were downloaded from \url{https://www.data.gouv.fr/es/datasets/} on 5th of January 2022.
}
\references{
Gámiz, M.L., Mammen, E., Martínez-Miranda, M.D. and Nielsen, J.P. (2024a). Low quality exposure and point processes with a view to the first phase of a pandemic. arXiv:2308.09918.

Gámiz, M.L., Mammen, E., Martínez-Miranda, M.D. and Nielsen, J.P. (2024b). Monitoring a developing pandemic with available data. arXiv:2308.09919.
}

\examples{
data("covid")
Hi<-covid$Hospi
Hi<-Hi[-1]
Ri<-covid$Recov
Ri<-diff(Ri)
Di<-covid$Death
Di<-diff(Di)
M<-length(Di)
## New hospitalizations are Hi-Ri-Di
newHi<-Hi[-1]-(Hi[-M]-Ri[-M]-Di[-M])
newHi<-c(Hi[1],newHi)
newHi[newHi<0]<-0 # (inconsistencies in the data)
plot(covid$Date[-1],newHi,xlab='Notification date',ylab='Number of persons',
main='New hospitalizations',type='l')
}
\keyword{datasets}
