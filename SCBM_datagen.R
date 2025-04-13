################################################################################
# Module: Datagen
# Description: This script is part of the SCBM framework and implements
# core functionality related to datagen.
################################################################################


# ====================
# Load Required Libraries
# ====================
library(grf)

# --------------------------------------------------
# Function: transOUT
# Purpose : [Add description of what the function does]
# Arguments: [Add argument descriptions if needed]
# --------------------------------------------------
transOUT<- function(Y,Z,X){
# --- Construct a Data Frame ---
  datq=data.frame(Z,X)
  mdl=regression_forest(X=X,Y=Z)
# --- Construct a Data Frame ---
  XX=data.frame(X)
  Zpre=predict(mdl,newdata=XX)
  Zpre[which(Zpre>0.90)]=0.90
  Zpre[which(Zpre<0.10)]=0.10
  Zpe=Z*Y/Zpre+(1-Z)*(-Y)/(1-Zpre)
# --- Return Final Output ---
  return(Zpe)
}

# --------------------------------------------------
# Function: transOUT.rct
# Purpose : [Add description of what the function does]
# Arguments: [Add argument descriptions if needed]
# --------------------------------------------------
transOUT.rct<- function(Y,Z,X){
# --- Construct a Data Frame ---
  dat=data.frame(Z,X)
# --- Construct a Data Frame ---
  XX=data.frame(X)
  Zpre=0.5
  Zpe=Z*Y/Zpre+(1-Z)*(-Y)/(1-Zpre)
# --- Return Final Output ---
  return(Zpe)
}
pmpv=function(X){
  for (i in 1:length(X)) {
    print(paste0("[[",i,"]]"))
    print(mean(X[[i]]))
    print(var(X[[i]]))
  }
}

# --------------------------------------------------
# Function: fnum_tau
# Purpose : [Add description of what the function does]
# Arguments: [Add argument descriptions if needed]
# --------------------------------------------------
fnum_tau<- function(X){
  # mx=c(3,2,3,2,3,2,3,2,3,2,3,2)
  # tx=c(5,6,7,8,9,10,5,6,7,8,9,10)
  
  f1 <- rep(0,nrow(X))
  
  f2 <- -3*X[,1]-3*X[,3]
  
  f3 <- f2
  
  f4 <- 0
  
  f5 <- 0
  
  f6 <- X[,1]+2*X[,3]+X[,5]+2*X[,7]+X[,9]#new
  
  f7 <- 5*as.numeric(X[,1]>1)*as.numeric(X[,7] > 0) + 5*as.numeric(X[,5] > 1)-5#new
  
  f8 <- X[,1]*X[,5] + 2*X[,3]*X[,7] + 2*X[,2]*X[,4] -2 #new
  
  f9 <- 2*sin(pi*X[,1]*X[,7]) + 2*(X[,3]+0.5)^2 + X[,9] -2 #new
  
  f10 <- 1/(1+exp(-X[,1]))+2*X[,3]^2+2*X[,9]-2#new
  
  f11 <- 0
  
  return (list(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11))
}



fnum_mu<- fnum_tau
 

# --------------------------------------------------
# Function: Xgen
# Purpose : [Add description of what the function does]
# Arguments: [Add argument descriptions if needed]
# --------------------------------------------------
Xgen<- function(n,p){
  di=1:p
  n=n
  ji <- seq(from=1,to=length(di),by=2)
  ou <- setdiff(di,ji)
  X=matrix(numeric(n*length(di)),nrow = n)
  for (i in 1:length(di)) {
    if(i%in%ji){
      X[,i]=rnorm(n)
    }
    else{
      X[,i]=rbinom(n,1,0.5)
    }
  }
# --- Return Final Output ---
return(X)
}

# --------------------------------------------------
# Function: datagen
# Purpose : [Add description of what the function does]
# Arguments: [Add argument descriptions if needed]
# --------------------------------------------------
datagen<- function(mu,tau,type,X,noise=1){
  n=nrow(X)
  if (type==1) {
    t=rbinom(n,1,0.5)
  }else{
    t=rbinom(n, 1, 0.3 + 0.2 * (X[, 1] > 1)+ 0.2 * (X[, 3] > 1))

  }
  
  y=mu+(t-0.5)*tau+rnorm(n,0,noise)
  

  if (type==1) {#response
# --- Construct a Data Frame ---
    dat=data.frame(transOUT.rct(y,t,X),t,tau,X,y)
  }else{
# --- Construct a Data Frame ---
    dat=data.frame(transOUT(y,t,X),t,tau,X,y)
  }

# --- Return Final Output ---
return(dat)
}