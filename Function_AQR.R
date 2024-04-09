#################################################
###  AQR(alternating quantile regression), estimate B1 B2 B3
#################################################
library(rTensor)
library(Matrix)
library(glmnet)
library(quantreg)


Int_q <- function(RANK,y1,y2,y3,vec_X1,vec_X2,vec_X3){
  ### RANK:predetermined Tucker rank
  ### y1 y2 y3: three groups' response vector
  ### vec_X1 vec_X2 vec_X3: each row is vectorized X
  
  
  B1_int <- as.tensor(array( rnorm(prod(dim),0,1) ,dim = dim))
  B2_int <- as.tensor(array( rnorm(prod(dim),0,1) ,dim = dim))
  B3_int <- as.tensor(array( rnorm(prod(dim),0,1) ,dim = dim))
  ###### HOSVD to get U1 U2 U3 G1 G2 G3 #######
  hosvd_B1 = hosvd(B1_int,ranks=RANK)
  hosvd_B2 = hosvd(B2_int,ranks=RANK)
  hosvd_B3 = hosvd(B3_int,ranks=RANK)
  
  G1_int=hosvd_B1$Z
  G2_int=hosvd_B2$Z
  G3_int=hosvd_B3$Z
  
  U1_1_int=hosvd_B1$U[[1]]
  U2_1_int=hosvd_B1$U[[2]]
  U3_1_int=hosvd_B1$U[[3]]
  
  U1_2_int=hosvd_B2$U[[1]]
  U2_2_int=hosvd_B2$U[[2]]
  U3_2_int=hosvd_B2$U[[3]]
  
  U1_3_int=hosvd_B3$U[[1]]
  U2_3_int=hosvd_B3$U[[2]]
  U3_3_int=hosvd_B3$U[[3]]
  
  U1_int=(U1_1_int+ U1_2_int +U1_3_int)/3
  U2_int=(U2_1_int+ U2_2_int +U2_3_int)/3
  U3_int=(U3_1_int+ U3_2_int +U3_3_int)/3
  
  rtn = list(B1_int=B1_int,B2_int=B2_int,B3_int=B3_int,G1_int=G1_int,G2_int=G2_int,G3_int=G3_int,
             U1_int=U1_int,U2_int=U2_int,U3_int=U3_int)
  return(rtn)
}



main_q <- function(RANK, X,y1,y2,y3,vec_X1,vec_X2,vec_X3,it_max=20,ep=0.01){
  ### RANK:predetermined Tucker rank
  ### X: 4-order tensor X[,,,i] is i-th sample
  ### y1 y2 y3: three groups' response vector
  ### vec_X1 vec_X2 vec_X3: each row is vectorized X
  ### it_max: maximum number of iterations
  ### ep: convergence threshold
  
  #########initial estimator######
  int = Int_q(RANK,y1,y2,y3,vec_X1,vec_X2,vec_X3)
  g1_old = int$G1_int
  g2_old = int$G2_int
  g3_old = int$G3_int
  u1_old = int$U1_int
  u2_old = int$U2_int
  u3_old = int$U3_int
  
  Loss=c()
  k=1
  while (k<=it_max) {
    u_list_123_old <- list('mat1'=u1_old,'mat2'=u2_old,'mat3'=u3_old)
    b1_old <- ttl(g1_old,u_list_123_old,ms=c(1,2,3))
    b2_old <- ttl(g2_old,u_list_123_old,ms=c(1,2,3))
    b3_old <- ttl(g3_old,u_list_123_old,ms=c(1,2,3))
    
    ######## update U1 U2 U3  create new X ########
    #create new x_u1 to update U1
    row1_u1 <- matrix(0,n1,prod(dim))
    for (i in 1:n1){row1_u1[i,]=t(vec(k_unfold(X[,,,i],1)))}
    row1_x_u1 <- row1_u1%*%kronecker(kronecker(u3_old,u2_old)%*%
                                       base::t(k_unfold(g1_old,1)@data),diag(dim[1]))
    
    row2_u1 <- matrix(0,n2,prod(dim))
    for (i in 1:n2){row2_u1[i,]=t(vec(k_unfold(X[,,,i+n1],1)))}
    row2_x_u1 <- row2_u1%*%kronecker(kronecker(u3_old,u2_old)%*%
                                       base::t(k_unfold(g2_old,1)@data),diag(dim[1]))
    
    row3_u1 <- matrix(0,n3,prod(dim))
    for (i in 1:n3){row3_u1[i,]=t(vec(k_unfold(X[,,,i+n1+n2],1)))}
    row3_x_u1 <- row3_u1%*%kronecker(kronecker(u3_old,u2_old)%*%
                                       base::t(k_unfold(g3_old,1)@data),diag(dim[1]))
    
    x_u1 <- rbind(row1_x_u1 , row2_x_u1 , row3_x_u1)
    beta_u1 <- rq(c(y1,y2,y3)~x_u1, tau = tau)$coefficients[-1]
    u1_new <-  matrix(beta_u1,dim[1],RANK[1])
    u1_new <- apply(u1_new,2,FUN = function(x) x/sqrt(sum(x^2)))
    
    #create new x_u2 to update U2
    row1_u2 <- matrix(0,n1,prod(dim))
    for (i in 1:n1){row1_u2[i,]=t(vec(k_unfold(X[,,,i],2)))}
    row1_x_u2 <- row1_u2%*%kronecker(kronecker(u3_old,u1_new)%*%
                                       base::t(k_unfold(g1_old,2)@data),diag(dim[2]))
    
    row2_u2 <- matrix(0,n2,prod(dim))
    for (i in 1:n2){row2_u2[i,]=t(vec(k_unfold(X[,,,i+n1],2)))}
    row2_x_u2 <- row2_u2%*%kronecker(kronecker(u3_old,u1_new)%*%
                                       base::t(k_unfold(g2_old,2)@data),diag(dim[2]))
    
    row3_u2 <- matrix(0,n3,prod(dim))
    for (i in 1:n3){row3_u2[i,]=t(vec(k_unfold(X[,,,i+n1+n2],2)))}
    row3_x_u2 <- row3_u2%*%kronecker(kronecker(u3_old,u1_new)%*%
                                       base::t(k_unfold(g3_old,2)@data),diag(dim[2]))
    
    x_u2 <- rbind(row1_x_u2 , row2_x_u2 , row3_x_u2)
    beta_u2 <- rq(c(y1,y2,y3)~x_u2, tau = tau)$coefficients[-1]
    u2_new <-  matrix(beta_u2,dim[2],RANK[2])
    u2_new <- apply(u2_new,2,FUN = function(x) x/sqrt(sum(x^2)))
    
    #create new x_u3 to update U3
    row1_u3 <- matrix(0,n1,prod(dim))
    for (i in 1:n1){row1_u3[i,]=t(vec(k_unfold(X[,,,i],3)))}
    row1_x_u3 <- row1_u3%*%kronecker(kronecker(u2_new,u1_new)%*%
                                       base::t(k_unfold(g1_old,3)@data),diag(dim[3]))
    
    row2_u3 <- matrix(0,n2,prod(dim))
    for (i in 1:n2){row2_u3[i,]=t(vec(k_unfold(X[,,,i+n1],3)))}
    row2_x_u3 <- row2_u3%*%kronecker(kronecker(u2_new,u1_new)%*%
                                       base::t(k_unfold(g2_old,3)@data),diag(dim[3]))
    
    row3_u3 <- matrix(0,n3,prod(dim))
    for (i in 1:n3){row3_u3[i,]=t(vec(k_unfold(X[,,,i+n1+n2],3)))}
    row3_x_u3 <- row3_u3%*%kronecker(kronecker(u2_new,u1_new)%*%
                                       base::t(k_unfold(g3_old,3)@data),diag(dim[3]))
    
    x_u3 <- rbind(row1_x_u3 , row2_x_u3 , row3_x_u3)
    beta_u3 <- rq(c(y1,y2,y3)~x_u3, tau = tau)$coefficients[-1]
    u3_new <-  matrix(beta_u3,dim[3],RANK[3])
    # u3_new <- qr.Q(qr(matrix(beta_u3,dim[3],R3)))
    u3_new <- apply(u3_new,2,FUN = function(x) x/sqrt(sum(x^2)))
    
    
    ############ create y_mid_1,2,3 and  x_mid_1,2,3,then update G1 G2 G3 ######
    u_list_new <- list('mat1'=u3_new,'mat2'=u2_new,'mat3'=u1_new)
    x_mid_1 <- vec_X1%*%kronecker_list(u_list_new)
    x_mid_2 <- vec_X2%*%kronecker_list(u_list_new)
    x_mid_3 <- vec_X3%*%kronecker_list(u_list_new)
    g1_new<-as.tensor(array(rq(y1~x_mid_1, tau = tau)$coefficients[-1],dim=RANK))
    g2_new<-as.tensor(array(rq(y2~x_mid_2, tau = tau)$coefficients[-1],dim=RANK))
    g3_new<-as.tensor(array(rq(y3~x_mid_3, tau = tau)$coefficients[-1],dim=RANK))
    
    
    #############  update B1 B2 B3 ######
    u_list_123 <- list('mat1'=u1_new,'mat2'=u2_new,'mat3'=u3_new)
    b1_new <- ttl(g1_new,u_list_123,ms=c(1,2,3))
    b2_new <- ttl(g2_new,u_list_123,ms=c(1,2,3))
    b3_new <- ttl(g3_new,u_list_123,ms=c(1,2,3))
    
    
    ###### compute loss #############
    res = c( (y1-vec_X1%*%vec(b1_new)), (y2-vec_X2%*%vec(b2_new)), (y3-vec_X3%*%vec(b3_new)) )
    Loss_new = mean(res*(tau - (res<0)))
    
    ########### compute BIC ########
    df = K*prod(RANK)+sum(RANK*(dim-RANK))
    BIC = n*log(Loss_new)+(df/2)*log(n)
    
    
    
    if ((fnorm(b1_new-b1_old)+
         fnorm(b2_new-b2_old)+
         fnorm(b3_new-b3_old))< ep) {
      result <- list(k=k,u1=u1_new,u2=u2_new,u3=u3_new,
                     g1=g1_new,g2=g2_new,g3=g3_new,
                     b1=b1_new,b2=b2_new,b3=b3_new,BIC = BIC)
      return(result)
    }
    u1_old=u1_new;u2_old=u2_new;u3_old=u3_new
    g1_old=g1_new;g2_old=g2_new;g3_old=g3_new
    
    Loss[k]=Loss_new
    k=k+1  
  }
  result <- list(k=k,u1=u1_new,u2=u2_new,u3=u3_new,
                 g1=g1_new,g2=g2_new,g3=g3_new,
                 b1=b1_new,b2=b2_new,b3=b3_new,BIC = BIC)
  return(result)
  
}









