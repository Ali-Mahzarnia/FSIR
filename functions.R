
# to test through after generating data
# $s=5; basisno=25 ; X=X.obs ; kn=2
 

FSIR=function(Y,X,s,basisno,tt,kn){
  p=dim(X)[1]
  N=dim(X)[2]
  B.basis=create.bspline.basis(rangeval=c(0,1), nbasis=basisno)
  X.fd=smooth.basis(tt,y=X, fdParobj=B.basis)
  Xcoef=X.fd$fd$coefs
  Xcoefc=Xcoef-rowMeans(Xcoef)
  Gram=bsplinepen(B.basis, 0)
  L= Xcoefc%*%t(Xcoefc)%*%Gram/N
  #FPCA on L 
  #first psudue sqrt of gram :
  Geig=eig(Gram); Geigv=Geig$values ; Geigvec=Geig$vectors
  ind= (abs(Geigv)>10^(-15))
  #gram^0.5
  Ghalf= Geigvec[,ind] %*% diag(Geigv[ind]^(0.5)) %*% t(Geigvec[,ind])
  #Gram ^-0.5
  Ghalfneg= Geigvec[,ind] %*% diag(Geigv[ind]^(-0.5)) %*% t(Geigvec[,ind])
  #now use for pca
  Lpceig=eig(Ghalf%*%Xcoefc%*%t(Xcoefc)%*%Ghalf/N)
  Leigf=Ghalfneg%*%(Lpceig$vectors);  Leigv=Lpceig$values;
  #inverse of operator L
  #ind= (abs(Leigv)>10^(-15))
  Linv= Leigf[,1:kn] %*% diag(Leigv[1:kn]^(-1)) %*% t(Leigf[,1:kn])
  
  
  
#now bind to order based on Y and find E(X|Y)
  YX=rbind(Y,X)
  YX=YX[,order(YX[1,])]
  XY.fd=smooth.basis(tt,y=YX[-1,], fdParobj=B.basis)
  XYcoef=XY.fd$fd$coefs
  EXYcoef=apply(array(XYcoef, c(basisno, N/s, s)), c(1, 3), mean)
  EXYcoefc=EXYcoef-rowMeans(EXYcoef)
  Le=EXYcoefc%*%t(EXYcoefc)%*%Gram/s

  F=Linv%*%Le
  #FPCA on F
  Feig=eig(sqrt(Gram)%*%F%*%t(F)%*%sqrt(Gram))
  beta=Ghalfneg%*%Feig$vectors;  Feigv=Feig$values
  # z is 25*500 each row is 500 data projected to new dimension from biggest to smallest eigenvecs
  Z=t(beta)%*%Gram%*%Xcoef
  return(list(betahat=beta, projected=Z, eigvalue=Feigv))
}

