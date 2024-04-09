###########################
### Example 1 tensor mean regression
###########################
library(rTensor)
library(Matrix)
library(glmnet)
library(quantreg)

n1=n2=n3=150
rank=c(2,2,2)
dim=c(5,5,5)
n=n1+n2+n3
K=3


########### generate X y B in tensor mean regression #################
set.seed(2022)
e = rnorm(n,mean = 0, sd = 1)
X<-as.tensor(array(rnorm(prod(c(dim[1], dim[2], dim[3], n)),0,1), dim = c(dim[1], dim[2], dim[3], n)))
X1=X[,,,1:n1];X2=X[,,,(n1+1):(n1+n2)];X3=X[,,,(n1+n2+1):(n1+n2+n3)]


################### generate parameter B 123 #########################
G1=as.tensor(array(rnorm(prod(rank),0,1),dim = rank))
G2=as.tensor(array(rnorm(prod(rank),0,1),dim = rank))
G3=as.tensor(array(rnorm(prod(rank),0,1),dim = rank))

U1 <- qr.Q(qr(matrix(rnorm(dim[1]*rank[1],0,1),dim[1],rank[1])))
U2 <- qr.Q(qr(matrix(rnorm(dim[2]*rank[2],0,1),dim[2],rank[2])))
U3 <- qr.Q(qr(matrix(rnorm(dim[3]*rank[3],0,1),dim[3],rank[3])))

ls <- list(  'mat1' = U1,'mat2' = U2,'mat3' = U3  )
B1 <- ttl(G1,ls,ms=c(1,2,3))
B2 <- ttl(G2,ls,ms=c(1,2,3))
B3 <- ttl(G3,ls,ms=c(1,2,3))


y1=vector()#first group y_1
for (i in 1:n1){y1[i]=innerProd(X[,,,i],B1)+e[i]}
y2=vector()#second group y_2
for (i in (n1+1):(n1+n2)){y2[i-n1]=innerProd(X[,,,i],B2)+e[i]}
y3=vector()#third group y_3
for (i in (n1+n2+1):(n1+n2+n3)){y3[i-(n1+n2)]=innerProd(X[,,,i],B3)+e[i]}

vec_X1=matrix(0,n1,dim[1]*dim[2]*dim[3]) 
for (i in 1:n1){vec_X1[i,]=t(rTensor::vec(X[,,,i]))}
vec_X2=matrix(0,n2,dim[1]*dim[2]*dim[3])
for (i in 1:n2){vec_X2[i,]=t(rTensor::vec(X[,,,i+n1]))}
vec_X3=matrix(0,n3,dim[1]*dim[2]*dim[3])
for (i in 1:n3){vec_X3[i,]=t(rTensor::vec(X[,,,i+n1+n2]))}



######### select rank in tensor mean regression ###################
rankgrid = as.matrix( expand.grid(c(2,3,4),c(2,3,4),c(2,3,4)) )
BIC = vector()
for (i in 1:nrow(rankgrid)) {
  ranktemp = rankgrid[i,]
  rank_BIC = try(main_m(RANK=ranktemp,X,y1,y2,y3,vec_X1,vec_X2,vec_X3))
  if("try-error" %in%class(rank_BIC)){cat("error when rank =", ranktemp);next}
  BIC[i] = rank_BIC$BIC
  print(i)
}
rankhat = rankgrid[which.min(BIC),]



########  estimate in tensor mean regression  ################
fit_ours = main_m(RANK=rankhat,X,y1,y2,y3,vec_X1,vec_X2,vec_X3)


#estimation error of tensor mean regression
c(fnorm(fit_ours$b1-B1)/fnorm(B1),fnorm(fit_ours$b2-B2)/fnorm(B2),fnorm(fit_ours$b3-B3)/fnorm(B3))








###########################
### Example 2 tensor quantile regression
###########################
library(rTensor)
library(Matrix)
library(glmnet)
library(quantreg)


n1=n2=n3=150
rank=c(2,2,2)
tau=0.5
dim=c(5,5,5)
n=n1+n2+n3
K=3



########### generate X y B in tensor quantile regression #################
set.seed(2022)

q_tau = qnorm(tau,0,1)
e = rnorm(n,mean = -q_tau, sd = 1) 

X<-as.tensor(array(rnorm(prod(c(dim[1], dim[2], dim[3], n)),0,1), dim = c(dim[1], dim[2], dim[3], n)))
X1=X[,,,1:n1];X2=X[,,,(n1+1):(n1+n2)];X3=X[,,,(n1+n2+1):(n1+n2+n3)]


#################### parameter B 123 #######################
G1=as.tensor(array(rnorm(prod(rank),0,1),dim = rank))
G2=as.tensor(array(rnorm(prod(rank),0,1),dim = rank))
G3=as.tensor(array(rnorm(prod(rank),0,1),dim = rank))

U1 <- qr.Q(qr(matrix(rnorm(dim[1]*rank[1],0,1),dim[1],rank[1])))
U2 <- qr.Q(qr(matrix(rnorm(dim[2]*rank[2],0,1),dim[2],rank[2])))
U3 <- qr.Q(qr(matrix(rnorm(dim[3]*rank[3],0,1),dim[3],rank[3])))

ls <- list(  'mat1' = U1,'mat2' = U2,'mat3' = U3  )
B1 <- ttl(G1,ls,ms=c(1,2,3))
B2 <- ttl(G2,ls,ms=c(1,2,3))
B3 <- ttl(G3,ls,ms=c(1,2,3))


y1=vector()#first group y_1
for (i in 1:n1){y1[i]=innerProd(X[,,,i],B1)+e[i]}
y2=vector()#second group y_2
for (i in (n1+1):(n1+n2)){y2[i-n1]=innerProd(X[,,,i],B2)+e[i]}
y3=vector()#third group y_3
for (i in (n1+n2+1):(n1+n2+n3)){y3[i-(n1+n2)]=innerProd(X[,,,i],B3)+e[i]}

vec_X1=matrix(0,n1,dim[1]*dim[2]*dim[3])
for (i in 1:n1){vec_X1[i,]=t(rTensor::vec(X[,,,i]))}
vec_X2=matrix(0,n2,dim[1]*dim[2]*dim[3])
for (i in 1:n2){vec_X2[i,]=t(rTensor::vec(X[,,,i+n1]))}
vec_X3=matrix(0,n3,dim[1]*dim[2]*dim[3])
for (i in 1:n3){vec_X3[i,]=t(rTensor::vec(X[,,,i+n1+n2]))}


######### select rank in tensor quantile regression  ###################
rankgrid = as.matrix( expand.grid(c(2,3,4),c(2,3,4),c(2,3,4)) )
BIC = vector()
for (i in 1:nrow(rankgrid)) {
  ranktemp = rankgrid[i,]
  rank_BIC = try(main_q(RANK=ranktemp,X,y1,y2,y3,vec_X1,vec_X2,vec_X3))
  if("try-error" %in%class(rank_BIC)){cat("error when rank =", ranktemp);next}
  BIC[i] = rank_BIC$BIC
  print(i)
}
rankhat = rankgrid[which.min(BIC),]



###################  estimate in tensor quantile regression ###################
fit_ours = main_q(RANK=rankhat,X,y1,y2,y3,vec_X1,vec_X2,vec_X3)

#estimation error of tensor quantile regression
c(fnorm(fit_ours$b1-B1)/fnorm(B1),fnorm(fit_ours$b2-B2)/fnorm(B2),fnorm(fit_ours$b3-B3)/fnorm(B3))



