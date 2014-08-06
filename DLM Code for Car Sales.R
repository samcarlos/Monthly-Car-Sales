data=read.csv("/merged.data.csv")

data=data[,-which(colSums(is.na(data))>0)]
data=data[,-which(colSums(data==0)>0)]
data=data[,-c(1,31,32)]

library(dlm)
library(FKF)

#function that takes a vector and block repeats it n times
repBdiag=function(x,n){
  Bdiag=x
  for(i in 2:n){
    Bdiag=bdiag(Bdiag, x)
  }
  return(Bdiag)
}


y <- ts(data, frequency = 1)
k=dim(y)[2]
#should be set much higher, but due to slow speed only used 3000 for prototype
MC <- 1000
TT <- nrow(y)
nAhead=12

#setting up arrays to keep output
gibbsTheta <- array(0, dim=c(TT+1,2*k+k*11, MC-1)) 
gibbsV <- array(0, dim=c(k,k, MC))
gibbsWmu <- array(0, dim=c(k,k, MC)) 
gibbsWbeta <- array(0, dim=c(k,k, MC)) 
gibbsWseas <- array(0, dim=c(11*k,11*k, MC)) 
gibbsThetaPredicted <- array(0, dim=c(nAhead,k, MC-1)) 


seasonal=c(1,rep(0,10))
## prior hyperparameters
delta0 <- 100
delta1 <- 100;
delta2 <- 50
delta3=100

V0 <- (delta0-2) *diag(c(rep(10^2, k)))
Wmu0 <- (delta1-2) * diag(rep(100,k))
Wbeta0 <- (delta2 -2) * diag(rep(2,k))
Wseas0=(delta3-2) * diag(rep(seasonal,k)*10)


mod <- dlm(FF = cbind(matrix(c(1,0),nrow=1) %x% diag(k), t(repBdiag(seasonal,k))),
           V = diag(k),
           GG = bdiag(matrix(c(1,0,1,1),2,2) %x% diag(k),repBdiag(dlmModSeas(12)$GG,k)),
           W = bdiag(diag(k), diag(k), repBdiag(diag(seasonal),k)),
           m0 = c(y[1,],rep(0,k), rep(0,11*k)),
           C0 = diag(x = 1e7, nrow = (2*k+k*11)))
# starting values 

mod$V <- gibbsV[,,1] <- V0/(delta0-2)
gibbsWmu[,,1] <- Wmu0/(delta1-2)
gibbsWbeta[,,1] <- Wbeta0/(delta2-2)
gibbsWseas[,,1] <- Wseas0/(delta3-2)

mod$W <- bdiag(gibbsWmu[,,1], gibbsWbeta[,,1],gibbsWseas[,,1])

ct=rep(0,k)
dt=rep(0, (2*k+k*11))
P0=as.matrix(diag(mod$C0))

y1=t(matrix(as.double(y),dim(y)[1], dim(y)[2]))

seasonal.loc=which(diag(Wseas0)>0)+2*k


library(FKF)

ct1=as.matrix(rep(0,k))
FF=mod$FF
GG=mod$GG

library(MASS)
sh=2+TT/2
lp=matrix(0,108,29)
for(it in 1: (MC-1)){
  # kalman filter stage
  modelFKF=fkf(a0=mod$m0,P0=mod$C0, dt=rep(0, dim(mod$W)[1]), ct=ct1, Tt=mod$GG, Zt=mod$FF, HHt=mod$W,GGt=mod$V, yt=(y1))
  #FFBS algorithm
  theta=matrix(0,109,dim(mod$GG))
  theta[109,]=modelFKF$att[,108]+mvrnorm(n=1, rep(0,377),modelFKF$Ptt[,,108])
  for(i in 108:1){
    h=modelFKF$att[,i]+modelFKF$Ptt[,,i]%*%t(mod$GG)%*%solve(modelFKF$Pt[,,i+1])%*% (theta[i+1,]-modelFKF$at[,i+1])
    H=modelFKF$Ptt[,,i]-modelFKF$Ptt[,,i]%*%t(mod$GG)%*%solve(modelFKF$Pt[,,i+1])%*%mod$GG%*%modelFKF$Ptt[,,i]
    
    H[which(abs(H)<1e-3)]=0
    theta[i,]=mvrnorm(n=1, h,H)
  }
  #predict future states
  ytsp <- tsp(mod$m0)
  p <- length(mod$m0)
  m <- nrow(mod$FF)
  a <- rbind(modelFKF$att[,108], matrix(0, nAhead, p))
  R <- vector("list", nAhead + 1)
  R[[1]] <- modelFKF$Ptt[,,108]
  f <- matrix(0, nAhead, m)
  Q <- vector("list", nAhead)
  for (t in 1:nAhead) {
    a[t + 1, ] <- mod$GG %*% a[t, ]
    R[[t + 1]] <- mod$GG %*% R[[t]] %*% t(mod$GG) + mod$W
    f[t, ] <- mod$FF %*% a[t + 1, ]
    Q[[t]] <- mod$FF %*% R[[t + 1]] %*% t(mod$FF) + mod$V
  }
  
  #update parameters and save outpute
  gibbsThetaPredicted[,,it]=f 
  gibbsTheta[,,it]=theta
  # update V
  S <- crossprod(y-theta[-1,] %*% t(mod$FF)) + V0
  gibbsV[,,it+1] <- solve(rwishart(df=delta0+1+TT,p=k,Sigma=solve(S)))
  mod$V <- gibbsV[,,it+1]
  # update Wmu and Wbeta
  theta.center <- theta[-1,]-(theta[-(TT),]  %*% t(mod$GG))
  SS1 <- crossprod(theta.center)[1:k,1:k] + Wmu0
  SS2 <- crossprod(theta.center)[(k+1):(2*k),(k+1):(2*k)] + Wbeta0
  SS3=.001+diag(crossprod(theta.center))[seasonal.loc]/2
  psi <- rgamma(k, shape = sh, rate = SS3)
  gibbsWmu[,,it+1] <- solve(rwishart(df=delta1+1+TT, Sigma=solve(SS1)))
  gibbsWbeta[,,it+1] <- solve(rwishart(df=delta2+1+TT, Sigma=solve(SS2)))
  temp.mat=matrix(0,k*11,k*11)
  diag(temp.mat)[seasonal.loc-2*k]=1/psi
  gibbsWseas[,,it+1]=temp.mat
  
  mod$W <- bdiag(gibbsWmu[,,it+1], gibbsWbeta[,,it+1],gibbsWseas[,,it+1])     
}


