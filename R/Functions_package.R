################################################################################
##### Hazard/rates estimation
################################################################################

####
## 2-dimensional local linear hazard with full information
####
hazard2D<-function(t.grid,z.grid,o.zt,e.zt,bs.grid,cv=FALSE)
{
  O.tz<-t(o.zt)
  E.tz<-t(e.zt)

  # Univariate kernel: sextic kernel
  K<-function(u) { return(((3003/2048)*(1-(u)^2)^6)*(abs(u)<1))}

  ## 1. Compute kernel evaluations and moments of the kernel
  ## 1.1. Kernel evaluations at dimension t
  M<-length(t.grid)
  Tt<-matrix(rep(t.grid, times=M),nrow = M, ncol = M,byrow=FALSE)
  # Tt is an MxM matrix with the values of the grid
  Tx<-Tt-t(Tt)
  nb<-nrow(bs.grid)
  K.bt<-array(0,dim=c(M,M,nb))
  # each element of K.bt is the Kernel evaluated at ut below
  ut<-array(0,dim=c(M,M,nb))
  for(b in 1:nb)
  {
    ut[,,b]=Tx/bs.grid[b,2]
    K.bt[,,b]<-apply(ut[,,b],1:2,K)
    #this works and returns the value of the Kernel function at each individual value of the array ut
    K.bt[,,b]<-K.bt[,,b]/(bs.grid[b,2])
  }

  ## 1.2. Dimension z (evaluation points are the same but bandwidth not necessarily)
  if (any(t.grid!=z.grid)){
    Zt<-matrix(rep(z.grid, times=M),nrow = M, ncol = M,byrow=FALSE)
    Zx<-Zt-t(Zt)
  } else Zx<-Tx
  K.bz<-array(0,dim=c(M,M,nb))
  ut<-array(0,dim=c(M,M,nb))
  for(b in 1:nb)
  {
    ut[,,b]=Zx/bs.grid[b,1]
    K.bz[,,b]<-apply(ut[,,b],1:2,K)
    #this works and returns the value of the Kernel function at each individual value of the array ut
    K.bz[,,b]<-K.bz[,,b]/(bs.grid[b,1])
  }
  rm(ut,Tt)

  ## 1.3. Moment calculations involved in the estimator
  c0<-c10<-c11<-d00<-d01<-d11<-det.D<-den<-array(0,dim=c(M,M,nb))
  for (b in 1:nb)
  {
    c0[,,b]<-(K.bt[,,b]) %*% (E.tz) %*%t( K.bz[,,b])
    c10[,,b]<-(K.bt[,,b] * Tx) %*% (E.tz) %*%t(K.bz[,,b])
    c11[,,b]<-(K.bt[,,b]) %*% (E.tz) %*% t(K.bz[,,b]*Zx)
    d00[,,b]<-(K.bt[,,b] * (Tx^2)) %*% (E.tz) %*% t(K.bz[,,b])
    d01[,,b]<-(K.bt[,,b] * Tx) %*% (E.tz) %*% t(K.bz[,,b]*Zx)
    d11[,,b]<-(K.bt[,,b]) %*% (E.tz) %*% t(K.bz[,,b]*(Zx^2))
    det.D[,,b]<- d00[,,b]*d11[,,b]-d01[,,b]^2
    den[,,b]<-det.D[,,b]*c0[,,b] -
      ( c10[,,b]*(d11[,,b]*c10[,,b]-d01[,,b]*c11[,,b])
        + c11[,,b]*(d00[,,b]*c11[,,b] - d01[,,b]*c10[,,b]) )
  }

  ## 2. Create a function to compute the LL hazard from results above
  ##    and for one single bandwidth
  hazard.b<-function(K.bt,K.bz,Tx,Zx,b,O.tz,c10,c11,d00,d01,d11,det.D,den)
  {
    num1<-((K.bt[,,b]) %*% (O.tz)  %*%  t(K.bz[,,b])) *(det.D[,,b])
    num2<-( (K.bt[,,b] * Tx) %*% (O.tz)  %*%  t(K.bz[,,b]))*(d11[,,b]*c10[,,b] - d01[,,b]*c11[,,b])
    num3<-( (K.bt[,,b]) %*% (O.tz)  %*%  t(K.bz[,,b]*Zx))*(d00[,,b]*c11[,,b] - d01[,,b]*c10[,,b])
    num<-num1-num2-num3
    estim<-num/den[,,b] ### alpha(t,z)
    estim[den[,,b]==0]<-NA
    estim<-t(estim);    ## alpha(z,t)
    # First dimension is z= start and second is t = duration
    ## Moreover 1 <= z <= M; and 1 <= t <= M-z+1=z2 (end)
    ## so, we cannot estimate beyond z2, that is,
    ##  estim.zt[z,t]=0 for t in [M-z+1, M]
    M<-nrow(estim)
    for (z in 2:M){for (t in (M-z+2):M){estim[z,t]<-NA}}
    ## return a vector with dimension Mt*Mz=M^2
    estim<-as.vector(estim)
    estim[estim<0]<-NA
    return(estim)
  }


  ## 3. Compute now the local linear estimator with possible CV-bandwidth
  if(cv==TRUE & nb>1)
  {
    ## 3.1. The leave-one-out (loo) LL-hazard estimator
    # Create a function similar to hazard.b but leave-one-out
    hazard.loo.b<-function(K.bt,K.bz,Tx,Zx,b,O.tz,c10,c11,d00,d01,d11,
                           det.D,den,M)
    {
      estim.ij<-double(M*M)
      ij<-l<-0
      for (i in 1:M)
      {
        for (j in 1:M)
        {
          l<-l+1
          if(j>M-i+1){
            estim.ij[l]<-NA
          } else {
            if (O.tz[i,j]>0) {
              Oij.prev<-O.tz[i,j];O.tz[i,j]<- O.tz[i,j]-1;ij=1
            }
            num1<-( (K.bt[i,,b]) %*% (O.tz)  %*%  (K.bz[j,,b])) *(det.D[i,j,b])
            num2<-( (K.bt[i,,b] * Tx[i,]) %*% (O.tz)  %*%  K.bz[j,,b])*(d11[i,j,b]*c10[i,j,b] - d01[i,j,b]*c11[i,j,b])
            num3<-( (K.bt[i,,b]) %*% (O.tz)  %*%  (K.bz[j,,b]*Zx[j,]) )*(d00[i,j,b]*c11[i,j,b] - d01[i,j,b]*c10[i,j,b])
            num<-num1-num2-num3
            estim<-num/den[i,j,b]; estim[den[i,j,b]==0]<-NA
            estim.ij[l]<-estim
            if (ij>0) {O.tz[i,j]<-Oij.prev; ij=0}
          }
        }
      }
      return(estim.ij)
    }

    # Finally evaluate the above function above for each bandwidth in the grid
    estim.loo.zt.bs<-sapply(1:nb, function(b){
      return(hazard.loo.b(K.bt,K.bz,Tx,Zx,b,O.tz,c10,
                          c11,d00,d01,d11,det.D,den,M))})
    # estim.loo.zt.bs is a matrix with Mt*Mz rows and nb columns

    ## 3.2. The usual estimator (not loo)
    estim.zt.bs<-sapply(1:nb, function(b){
      return(hazard.b(K.bt,K.bz,Tx,Zx,b,O.tz,c10,
                      c11,d00,d01,d11,det.D,den))})
    ## estim.zt.bs is a matrix with Mt*Mz rows and nb columns

    ## 3.3. Compute the CV-banwidth choice with the above matrices
    vec.E.zt<-as.double(e.zt)
    vec.O.zt<-as.double(o.zt)
    #  The cross-validation function
    CV.b<-function(b)
    {
      dif1.b<-sum((estim.zt.bs[,b])^2 ,na.rm=TRUE)
      dif2.b<-sum(estim.loo.zt.bs[,b]* (vec.O.zt/vec.E.zt),na.rm=TRUE)
      cv.b<-dif1.b-2*dif2.b
      if (cv.b==Inf | cv.b==0) cv.b<-NA
      return(cv.b)
    }
    cv.values <- sapply(1:nb, CV.b)
    bcv<-which.min(cv.values)
    hi.zt<-estim.zt.bs[,bcv]
    hi.zt<-matrix(hi.zt,M,M,byrow=FALSE)
    hi.zt[which(hi.zt<0)]<-NA
  } else {
    estim.zt.b<-hazard.b(K.bt,K.bz,Tx,Zx,b=1,O.tz,c10,
                         c11,d00,d01,d11,det.D,den)
    hi.zt<-matrix(as.double(estim.zt.b),M,M,byrow=FALSE)
    hi.zt[which(hi.zt<0)]<-NA
    bcv<-1
  }

  return(list(hi.zt=hi.zt,bcv=t(bs.grid[bcv,])))
}

####
## Our algorithm to estimate the survival hazards from truncated data
####
hazard2Dmiss<-function(t.grid,z.grid,Oi1.z,Oi2.z,Ei.z,bs.grid,cv=TRUE,
                         epsilon=1e-4,max.ite=50)
{
  ## sample information is as follows:
  # Each day z, we get information on number of people remaining in hospital: Ei.z
  # as well as deaths (Oi1.z) and recoveries (Oi2.z) on that day (z)
  # Oi1.z= number of deaths notified on the day z;
  # Oi2.z=number of recoveries notified on the day z
  # Ei.z= number of people in the hospital on the day z, they have arrived at any day from 1:z;

  Oi.z<-Oi1.z+Oi2.z  # number of recoveries+deaths
  M<-length(Oi.z)    # total number of days considered in the sample



  #################################################################################
  # Construction of ocurrences and exposures
  #################################################################################

  ## Step 1. Initialization: construct them from an initial guess (exponential)

  Ei.zt<-Oi.zt<-matrix(0,M,M)  # NOT OBSERVED! only the sums by rows are!!
  #Oi.zt[z,t]<- number of subjects that leave the hospital (die or recover) the day z and have duration in hospital equal t
  ## stay in hospital from the day z-t+1 to the day z;
  ## columns are duration  (t) and  rows are the reporting day (z)
  #Ei.zt[z,t]<- number of subjects staying in the hospital on the day z and have duration in hospital equal t
  # they have arrived at any day in the interval (z-t+1,z) and still are in hospital on the day z
  # columns denote duration (t) and  rows are the reporting day (z)
  Ei.new<-Ei.z[-1]-(Ei.z[-M]-Oi.z[-M])
  Ei.new<-c(Ei.z[1],Ei.new);
  Ei.new[Ei.new<0]<-0   #to overcome certain data incongruences
  Ei.zt[,1]<-as.integer(Ei.new) # new arrivals each day:


  # dimension z=notification day; dimension t=duration
  for(z in 1:M)
  {
    for(t in 1:(M-z+1)) ## we fill in the matrices following the diagonal-track
    {
      if( (Ei.zt[(z+t-1),t]>0) & (Oi.z[z+t-1]>0) )
      {
        Oi.zt[(z+t-1),t]<-Ei.zt[(z+t-1),t]/Oi.z[(z+t-1)]
      } else {
        Oi.zt[(z+t-1),t]<-0
      }
      if(t<(M-z+1)){
        Ei.zt[(z+t),(t+1)]<-Ei.zt[(z+t-1),t]-Oi.zt[(z+t-1),t]
      }
    }
  }


  ## To estimate Oi.zt, define:
  #  q(z,t) = O(z,t)/(sum_d' O(z,t'))
  #  q.zt is the density of occurrences by duration t
  #  conditioned to notification day z
  #  From q.zt we create:
  #  Oi.zt<- q(z,t)*Oi.z[z]

  total.Oz<-rowSums(Oi.zt)
  q.zt<-Oi.zt/total.Oz
  q.zt[which(is.na(q.zt))]<-0
  Oi.zt<-q.zt*Oi.z;

  ## To estimate Ei.zt, define:
  #  h(z,t) = E(z,t')/(sum_d E(z,t'))
  #  h(z,t) is the density of exposure by duration t conditioned
  #  to notification day z
  #  From h(z,t) we construct:
  #  Ei.zt<- sum_z h(z,t)*Ei.z[z]
  total.Ez<-rowSums(Ei.zt) ## = Ei.z #the exposure by date (z)
  h.zt<-Ei.zt/total.Ez
  h.zt[is.na(h.zt)]<-0
  Ei.zt<-h.zt*Ei.z;

  Oi.zt[Ei.zt==0]<-0

  ## Rearrange the matrices to compute local linear estimators
  ## each row is marked by the date the subjects enter the hospital
  ## the matrices are triangular with no element below the secondary diagonal
  oi.zt<-matrix(0,M,M)
  ei.zt<-matrix(0,M,M)
  for(z in 1:M)
  {
    for(t in 1:(M-z+1))
    {
      ei.zt[z,t]<-Ei.zt[(z+t-1),t]
      oi.zt[z,t]<-Oi.zt[(z+t-1),t]
    }
  }

  est0<-hazard2D(t.grid,z.grid,o.zt=oi.zt,e.zt=ei.zt,bs.grid,cv=TRUE)
  hi.zt.0<-est0$hi.zt
  bcv0<-est0$bcv

  ## Step 2. Iterations until convergence or stopping criteria

  tol<-1 #initial tolerance: max(abs(alphai.zt-hi.zt)/alphai.zt,na.rm=T)
  it<-0
  while((tol>epsilon) & (it<max.ite))
  {
    it<-it+1
    hi.zt<-hi.zt.0
    bcv<-bcv0
    # Repeat step 1 but with estimated hazard for duration
    Oi.zt<-Ei.zt<-Si.zt<-pi.zt<-S0i.zt<-matrix(0,M,M)
    for(z in 1:M){
      Si.zt[z,]<-exp(-cumsum(hi.zt[z,]))
      S0i.zt[z,]<-c(1,Si.zt[z,1:(M-1)]);
      pi.zt[z,]<- 1-Si.zt[z,]/S0i.zt[z,];kk<-which(S0i.zt[z,]==0)
      pi.zt[z,kk]<-1;pi.zt[z,is.na(pi.zt[z,])==T]<-1
    }
    Ei.zt[,1]<-as.integer(Ei.new)

    for(z in 1:M)
    {
      for(t in 1:(M-z+1))
      {
        if((Ei.zt[(z+t-1),t]>0)&(is.na(hi.zt[z,t])==FALSE)){
          Oi.zt[(z+t-1),t]<-Ei.zt[(z+t-1),t]*pi.zt[z,t]
        } else Oi.zt[(z+t-1),t]<-0
        if(t<(M-z+1)){Ei.zt[(z+t),(t+1)]<-Ei.zt[(z+t-1),t]-Oi.zt[(z+t-1),t]}
      }
    }

    total.Oz<-rowSums(Oi.zt) # = Oi
    total.Ez<-rowSums(Ei.zt) # = Ei
    ####

    q.zt<-Oi.zt/total.Oz
    h.zt<-Ei.zt/total.Ez
    q.zt[which(is.na(q.zt))]<-0;h.zt[which(is.na(h.zt))]<-0

    Oi.zt<-q.zt*Oi.z
    Ei.zt<-h.zt*Ei.z
    Oi.zt[Ei.zt==0]<-0

    oi.zt<-matrix(0,M,M)
    ei.zt<-matrix(0,M,M)

    for(z in 1:M)
    {
      for(t in 1:(M-z+1))
      {
        ei.zt[z,t]<-Ei.zt[(z+t-1),t]
        oi.zt[z,t]<-Oi.zt[(z+t-1),t]
      }
    }

    if (it<5){
      est<-hazard2D(t.grid,z.grid,o.zt=oi.zt,e.zt=ei.zt,bs.grid,cv=TRUE)
      hi.zt<-est$hi.zt
      bcv<-est$bcv
    } else {
      est<-hazard2D(t.grid,z.grid,o.zt=oi.zt,e.zt=ei.zt,
                       bs.grid=bcv,cv=FALSE) ## corrected 6th Feb 2024
      hi.zt<-est$hi.zt
    }

    tol<-(sum( (hi.zt.0-hi.zt)^2,na.rm=T))/(sum(hi.zt.0^2,na.rm=T)+1e-6)
    hi.zt.0<-hi.zt
    message('Iteration ',it, '. Tolerance=',tol)

  }

  ## Step 3 (final): estimate hazard for deaths and recoveries separately

  Oi1.zt<-q.zt*Oi1.z
  Oi2.zt<-q.zt*Oi2.z
  oi1.zt<-oi2.zt<-ei.zt<-matrix(0,M,M)
  for(z in 1:M)
  {
    for(t in 1:(M-z+1))
    {
      ei.zt[z,t]<-Ei.zt[(z+t-1),t]
      oi1.zt[z,t]<-Oi1.zt[(z+t-1),t]
      oi2.zt[z,t]<-Oi2.zt[(z+t-1),t]
    }
  }
  ## hazard estimate for deaths
  if (sum(Oi1.z)==0){ hi1.zt<-NA
  } else {
    est1<-hazard2D(t.grid,z.grid,o.zt=oi1.zt,e.zt=ei.zt,bs.grid,cv=TRUE)
    hi1.zt<-est1$hi.zt
    bcv1<-est1$bcv
  }

  ## estimate of recovery hazard:
  if (sum(Oi2.z)==0){hi2.zt<-NA
  } else {
    est2<-hazard2D(t.grid,z.grid,o.zt=oi2.zt,e.zt=ei.zt,bs.grid,cv=TRUE)
    hi2.zt<-est2$hi.zt
    bcv2<-est2$bcv
  }





  result<-list(hi.zt=hi.zt,hi1.zt=hi1.zt,hi2.zt=hi2.zt,
               bcv=c(bcv1,bcv2),tol=tol,it=it)
          ## it may return also the last generated occurrences and exposure
               # , o.zt=oi.zt,o1.zt=oi1.zt,o2.zt=oi2.zt,e.zt=ei.zt)

  return(result)
}


#   ## From the hazard we can compute other measures of interest :
#   ## Medians of the total time until outcome depending on time (z1)

# hi<-hi2.zt+hi1.zt # 2 possible outcomes
medtime<-function(hi.zt,z1)
{
  hi.zt[is.na(hi.zt)]<-0
  M<-nrow(hi.zt)
  if (missing(z1)) z1<-c(seq(1,M-1,by=2),M-1)
  S.list<-c()
  n<-length(z1)
  for(i in 1:n)
  {
    zi<-z1[i]
    hi<-hi.zt[zi,]
    Si<-list(cumprod(1-hi))
    S.list<-c(S.list,Si)
  }
  v.med<-double(n)
  times<-1:M
  for(i in 1:n)
  {
    ii<-which(S.list[[i]]<0.5)
    if (length(ii)==0) v.med[i]<-NA else v.med[i]<-times[min(ii)]
  }
  return(v.med)
}


#   ## 2. Probability of outcome (alive or death) depending on time (z1)
 
poutcome<-function(hi1.zt,hi2.zt,z1)
{
  M<-nrow(hi1.zt)
  if (missing(z1)) z1<-c(seq(1,M-1,by=2),M-1)
  alive<-function(td,hid,hir)
  {
    Md<-length(td)
    hid[is.na(hid)]<-0
    hir[is.na(hir)]<-0
    Si.all<-cumprod(1-(hid+hir))
    n<-length(Si.all)
    Si.0<-c(1,Si.all[-n])
    prob.r<-hir*Si.all
    p.alive<-cumsum(prob.r[n:1])
    p.alive<-p.alive[n:1]/Si.0
    p.death<-1-p.alive
    res<-list(p.alive=p.alive[1:Md],p.death=p.death[1:Md])
    return(res)
  }

  nz<-length(z1)
  alive.zt<-death.zt<-matrix(NA,M,nz)
  for(j in 1:nz)
  {
    ti<-1:(M-z1[j]+1)
    n<-length(ti)
    probs.j<-alive(td=ti,hid=hi1.zt[z1[j],1:n],hir=hi2.zt[z1[j],1:n])
    alive.zt[1:(M-z1[j]+1),j]<-probs.j$p.alive
    death.zt[1:(M-z1[j]+1),j]<-probs.j$p.death
  }
  result<-list(alive.zt=alive.zt, death.zt=death.zt)
  return(result)
}

## Our algorithm to estimate the rate of infection from truncated data
## The function below is a version of hazard2D_trunc() used to compute
## the estimation of the hospitalization rate and infection rate
## supplied arguments (Oi.z, Ei.z1) changing depending on which of the two we want
rate2Dmiss<-function(t.grid,z.grid,Oi.z,Ei.z1,bs.grid,cv=TRUE,
                       epsilon=1e-4,max.ite=50)
{
  # Oi.z= number of hospitalized notified on the day z;
  # Ei.z1= number new positive tested on the day z;
  # This is the constant exposure for being hospitalized
  ## The counting process we are interested is:
  ## N_z(t) = number of hospitalized among people that were positive
  ##on the day z-t+1
  ## This process has intensity lambda_z(t)= alpha_z(t)*Ei.zt[z,1]
  ## In our previous works: lambda_z(t)=alpha_z(t)*Ei.zt[z+t-1,t],
  ## with Ei.zt[z+t-1,t]=Ei.zt[z+t-2,t-1]-Oi.zt[z+t-2,t-1],
  ## in other words, the exposure is updated by removing the occurrences each day.
  ## Now the exposure is constant and equal to the number of new positive on the day z, z=1,2,..M
  ## We need modify this step of the old algorithm: Ei.zt[z+t-1,t]<-Ei.zt[z,1], for all t=1,2,...,M-z+1
  ## Finally, remember we are doing time-dependent hazards, being the marker variable the notification date= z


  M<-length(Oi.z)   # total number of days considered in the sample
  #################################################################################
  # Construction of ocurrences and exposures
  #################################################################################

  Oi.zt<-Ei.zt<-matrix(0,M,M) # NOT OBSERVED! only the sums by rows are!!

  ## Step 1. Initialization: construct them from an initial guess (exponential)

  for(z in 1:M) for(t in 1:(M-z+1)) Ei.zt[(z+t-1),t]<-Ei.z1[z]

  for(z in 1:M)
  {
    for(t in 1:(M-z+1))
    {
      if((Ei.zt[(z+t-1),t]>0)&(Oi.z[z+t-1]>0)) {                                   ### deterministic version!!
        Oi.zt[(z+t-1),t]<-Ei.zt[(z+t-1),t]/Oi.z[(z+t-1)]
      } else Oi.zt[(z+t-1),t]<-0
    }
  }

  total.Oz<-rowSums(Oi.zt)
  q.zt<-Oi.zt/total.Oz
  q.zt[which(is.na(q.zt))]<-0
  Oi.zt<-q.zt*Oi.z
  Oi.zt[Ei.zt==0]<-0

  oi.zt<-ei.zt<-matrix(0,M,M)

  for(z in 1:M)
  {
    for(t in 1:(M-z+1))
    {
      ei.zt[z,t]<-Ei.zt[(z+t-1),t]
      oi.zt[z,t]<-Oi.zt[(z+t-1),t]
    }
  }

  estim<-hazard2D(t.grid,z.grid,o.zt=oi.zt,e.zt=ei.zt,
                     bs.grid,cv=cv)
  hi.zt.0<-estim$hi.zt
  bcv.0<-estim$bcv

  ## Step 2. Iterations until convergence or stopping criteria
  tol<-1 #initial tolerance
  it<-0
  message('Running the algorithm, please be patient.')
  while((tol>epsilon) & (it<max.ite))
  {
    it<-it+1
    hi.zt<-hi.zt.0
    bcv<-bcv.0
    pi.zt<-Oi.zt<-matrix(0,M,M)

    for(z in 1:M)
    {
      pi.zt[z,]<-hi.zt[z,]*exp(-hi.zt[z,])
      for(t in 1:(M-z+1))
      {
        if((Ei.zt[(z+t-1),t]>0)&(is.na(pi.zt[z,t])==FALSE))
        {
          Oi.zt[(z+t-1),t]<-Ei.zt[(z+t-1),t]*pi.zt[z,t]
        }else Oi.zt[(z+t-1),t]<-0
      }
    }
    ### totals by row (should be Oi).
    total.Oz<-rowSums(Oi.zt) # = Oi
    q.zt<-Oi.zt/total.Oz
    q.zt[which(is.na(q.zt))]<-0
    #### estimated 2d-dimensional occurrences
    Oi.zt<-q.zt*Oi.z
    Oi.zt[Ei.zt==0]<-0

    ### Obtain the information of occurrences and exposure in terms of marker (z1) before passing on to csda13-code for hazard estimation:
    oi.zt<-matrix(0,M,M)

    for(z in 1:M) for(t in 1:(M-z+1)) oi.zt[z,t]<-Oi.zt[(z+t-1),t]

    if (it<5){
      estim<-hazard2D(t.grid,z.grid,o.zt=oi.zt,e.zt=ei.zt,
                       bs.grid,cv=cv)
    } else {
      estim<-hazard2D(t.grid,z.grid,o.zt=oi.zt,e.zt=ei.zt,
                         bs.grid=bcv,cv=FALSE)
    }
    hi.zt<-estim$hi.zt
    bcv<-estim$bcv

    tol<-(sum( (hi.zt.0-hi.zt)^2,na.rm=T))/(sum(hi.zt.0^2,na.rm=T)+1e-6)
    hi.zt.0<-hi.zt
    message('Iteration ',it, '. Tolerance=',tol)
  }

  result<-list(hi.zt=hi.zt,bcv=bcv,tol=tol,it=it,
  # return also the last generated occurrences and exposure
  o.zt=oi.zt,e.zt=ei.zt)

  return(result)
}

################################################################################
##### FORECASTING
################################################################################


## Function to computed forecasts of the number of new infected
## and new hospitalizations
forecasting<-function(Cval=1,period,RoInf,RoHosp,RoDeath,RoRec,
                      Pz,newHz,Hz,Dz,Rz)
{
  M<-nrow(RoInf)
  M1<-M+1
  Pz.obs<-Pz[1:M1]
  Hz.obs<-Hz[1:M1]
  newHz.obs<-newHz[1:M1]
  Dz.obs<-Dz[1:M1]
  Rz.obs<-Rz[1:M1]
  
  ############################
  #### Step 1: Forecasting for number of positives every day in the interval [(M+1)+1,(M+1)+period],
  #### the data for RoInf: (P_1,P_2,.....,P_{M},H_{M+1})
  alphas.I<-matrix(0,M,M) ### the maximum time-length (duration=t) we estimate the rate of infection is M
  for(j in 1:M){alphas.I[j,j:M]<-RoInf[j,1:(M-j+1)]} ## this is an upper-triangular matrix
  Pz.fitted<-colSums(Pz.obs[1:M]*alphas.I,na.rm=T)
  
  late.estim.I<-alphas.I[,M] ##
  # rate increases linearly up to Cval times
  v.C<-1+ ((Cval-1)/period) * (1:period)
  Aux2<-Aux1<-matrix(0,period+M,period)
  for(j in 1:period){Aux1[(j+1):(j+M),j]<-late.estim.I}
  for(k in 1:period){Aux2[(k+1):(k+M),k]<-v.C[k]}
  fc.alphas.I<-Aux1*Aux2 ### alphas.I
  
  Pz.pred<-Pz
  
  z<-1
  while(z<=period)
  {
    new.Pz.pred<-sum(Pz.pred*fc.alphas.I[1:length(Pz.pred),z],na.rm=T)
    Pz.pred<-c(Pz.pred,new.Pz.pred)
    z<-z+1
  }
  ############################
  
  if (missing(RoHosp)) {
    newHz.fitted<-rep(NA,M1);newHz.pred<-rep(NA,period)
    Dz.fitted<-rep(NA,M1);Dz.pred<-rep(NA,period)
    Rz.fitted<-rep(NA,M1);Rz.pred<-rep(NA,period)
    } else {   
  
      ### Step 2: Forecasting for number of NEW hospitalizations every day in the interval [(M+1)+1,(M+1)+period]
    ### Hospitalization rate is estimated for a maximum time-length of M1=M+1
    alphas.H<-matrix(0,M1,M1) ### the maximum time-length (duration=t) we estimate the rate of infection is M
    for(j in 1:M1){alphas.H[j,j:M1]<-RoHosp[j,1:(M1-j+1)]} ## this is an upper-triangular matrix
    newHz.fitted<-colSums(Pz.obs*alphas.H,na.rm=T)
    
    late.estim.H<-alphas.H[,M1] 
    
    fc.alphas.H<-matrix(0,M1+period,period) 
    
    for(j in 1:period){fc.alphas.H[(j+1):(j+M1),j]<-late.estim.H} ### we extrapolate the latest estimation of RoH with C=1
    
    ### To predict hospitalizations in the period [(M+1)+1,(M+1)+period] we need the predictions of infected
    ### just obtained above: Pz.pred which is a vector of dimension (M+1)+period
    newHz.pred<-colSums(Pz.pred[1:(M1+period)]*fc.alphas.H,na.rm=T)
    
    ############################
    if (missing(RoDeath)==FALSE & missing(RoRec)==FALSE) {
      ### Step 3: Forecasting for number of deaths and recoveries every day in the interval [(M+1)+1,(M+1)+period]
    ### Death-rate and Rec-rate are estimated for a maximum time-length of M1=M+1
    
    # Build matrices of alphas.D and alphas.R in the complete time period: (1,M+1+period)
    # the observed period first:
    alphas.D<-alphas.R<-matrix(0,M1,M1) ### the maximum time-length (duration=t) we estimate the rate of infection is M
    for(z in 1:M1){
      alphas.D[z,z:M1]<-RoDeath[z,1:(M1-z+1)]
      alphas.R[z,z:M1]<-RoRec[z,1:(M1-z+1)]
      } ## this is an upper-triangular matrix
    
    
    # now the extrapolation to (M+1, M+1+period):  
    late.estim.D<-alphas.D[,M1] 
    late.estim.R<-alphas.R[,M1]
    fc.alphas.D<-matrix(0,M1+period,period) 
    fc.alphas.R<-matrix(0,M1+period,period)
    for(z in 1:period){
      fc.alphas.D[(z+1):(z+M1),z]<-late.estim.D
      fc.alphas.R[(z+1):(z+M1),z]<-late.estim.R} ### we extrapolate the latest estimation of RoDeath and RoRec
    
    # all together:                                         
    zeros<-matrix(0,period,M1)
    alphas.pred.D<-rbind(alphas.D,zeros)
    alphas.pred.D<-cbind(alphas.pred.D,fc.alphas.D)
    alphas.pred.R<-rbind(alphas.R,zeros)
    alphas.pred.R<-cbind(alphas.pred.R,fc.alphas.R)
    
    D.zt<-R.zt<-H.zt<-matrix(0,M1+period,M1+period)
    diag(H.zt)<-c(newHz.fitted,newHz.pred)#[1:(M1+period)];
    diag(D.zt)<-diag(H.zt)*diag(alphas.pred.D)
    diag(R.zt)<-diag(H.zt)*diag(alphas.pred.R)
    
    for(z in 1:(M1+period-1)){
      for(t in (z+1):(M1+period)){
        H.zt[z,t]<-H.zt[z,t-1]-D.zt[z,t-1]-R.zt[z,t-1]
        D.zt[z,t]<-H.zt[z,t]*alphas.pred.D[z,t]
        R.zt[z,t]<-H.zt[z,t]*alphas.pred.R[z,t]
      }
    }
    Dz.fitted<-colSums(D.zt,na.rm=T)##
    Rz.fitted<-colSums(R.zt,na.rm=T)
    ###%%%
    
    ### Step 4: Forecasting for number of people hospitalized every day 
    Hz.fitted<-colSums(H.zt,na.rm=T)
    }
    }  
  
  return(list(Pz.fitted=Pz.fitted,
              Pz.pred=Pz.pred[((M1)+1):(M+1+period)],
              newHz.fitted=newHz.fitted,
              newHz.pred=newHz.pred,
              Dz.fitted=Dz.fitted[1:M1],
              Dz.pred=Dz.fitted[((M1)+1):(M+1+period)],
              Rz.fitted=Rz.fitted[1:M1],
              Rz.pred=Rz.fitted[((M1)+1):(M1+period)]))
}

