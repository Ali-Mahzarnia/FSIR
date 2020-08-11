#3rd example of betas that are not at origin
#EXAMPLE WITH DUIFFERENT INTIAL BETAS
library(fda)
library(matlib)
# Generate N=500 browinain motion, with p=100 observation
# inastead of time in (0,1) we take time between (0,100)
N=500
p=1000
X=matrix(0, ncol=N, nrow=p)
for(n in 1:N){
  X[, n]=cumsum(rnorm(p))/sqrt(p)
}
# AS WELL as generating the error term for 500 scalar response Y where error~N(0,0.1)
eps=matrix(0, ncol=N, nrow=1)
for(n in 1:N){
  eps[, n]=rnorm(1,0,0.1)
}
#true beta 1 and beta 2
beta1 = function(t){
  return((2*t-1)^3)}
beta2 = function(t){
  return(cos(pi*(2*t-1)))}

#evaulate true beta 1 and beta 2 on (0,1) , basically on each of p=1000 observation point
b1=matrix(0, ncol=1, nrow=p)
b2=matrix(0, ncol=1, nrow=p)
for(i in 0:p){
  j=i
  b1[i]=beta1(j/p)
  b2[i]=beta2(j/p)
}


#inner product of X and beta 1 and 2 and constructing Y
# Y=e(beta1.X)+e(|beta2*x|)+epsilon
Y=matrix(0, ncol=N, nrow=1)
#these two are inner products of x and beta1 and beta2
Xb1=matrix(0, ncol=N, nrow=1)
Xb2=matrix(0, ncol=N, nrow=1)
for(n in 1:N){
  Xb1[n]=(X[,n]%*%b1)/p
  Xb2[n]=(X[,n]%*%b2)/p
}
Y=exp(Xb1)+exp(abs(Xb2))+eps
# make a matrix
X.obs = X[(1:100)*p/100,]
tt=(1:100)/100
nbasisno=30

FSIRresult=FSIR(Y,X.obs,5,basisno= nbasisno,tt,3)

#eigenvalues 
#FSIRresult$eigvalue
#proportion of variability 
FSIRresult$eigvalue/sum(FSIRresult$eigvalue)
#recreate betas
B.basis=create.bspline.basis(rangeval=c(0,1), nbasis=nbasisno)
betahat=FSIRresult$betahat
par(mfrow=c(1,2))
f=fd(betahat[,1], B.basis)
plot(eval.fd(tt, f),ylim=c(-2,2))
lines(y=b1,x=(1:p)*100/p)
lines(y=-b1,x=(1:p)*100/p, col='red')

f=fd(betahat[,2], B.basis)
plot(eval.fd(tt, f),ylim=c(-2,2))
lines(y=b2,x=(1:p)*100/p)
lines(y=-b2,x=(1:p)*100/p, col='red')

# projected results :
par(mfrow=c(1,2))
Xresult=FSIRresult$projected
plot(Xb1,Xresult[1,])
plot(Xb2,Xresult[2,])

#install.packages("scatterplot3d") # Install
library("scatterplot3d") # load
par(mfrow=c(1,2))
scatterplot3d(Xb1,Xb2,Y, highlight.3d=TRUE, angle = 200)
scatterplot3d(Xresult[1,],Xresult[2,],Y, highlight.3d=TRUE, angle = 200)

par(mfrow=c(2,2))
plot(Xresult[1,],Y)
plot(Xb1,Y)

plot(Xresult[2,],Y)
plot(Xb2,Y)
