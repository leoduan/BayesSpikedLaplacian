require("rstiefel")

rcbmf <-
  function(A,B,C,X, n_tries=1E4)
  {
    #simulate from the matrix bmf distribution as described in Hoff(2009) 
    #this is one Gibbs step, and must be used iteratively
    
    #warning - do not start X at the eigenvectors of A:
    #The relative weights on A then can become infinity. 
    #Instead, you can start close to A, e.g. 
    # X<-rmf.matrix(UA[,1:R]*m)
    
    m<-dim(X)[1] ;  R<-dim(X)[2]
    x_bak = X[,1]
    
    
    if(m>R)
    {
      
      for(r in sample( c(2:R)))
      {
        N<-NullC(X[,-r])
        An<-B[r,r]*t(N)%*%(A)%*%N ; cn<- t(N)%*%C[,r]
        X[,r]<-N%*%rbmf.vector.gibbs(An,cn,t(N)%*%X[,r])
      }
  
    }
    
    if(m==R)
    {
      for(s in seq(1,R,length=R))
      {
        r<-sort(sample(c(2:n),2))
        N<-NullC( X[,-r]  )
        An<- t(N)%*%A%*%N
        Cn<- t(N)%*%C[,r]
        X[,r]<-N%*%rbmf.O2(An,B[r,r],Cn)
      }
    }
    
    #update the first column so that it stays positive
    {
      r = 1
      N<-NullC(X[,-r])
      An<-B[r,r]*t(N)%*%(A)%*%N ; cn<- t(N)%*%C[,r]
      X[,r]<-N%*%rbmf.vector.gibbs(An,cn,t(N)%*%X[,r])
      
      while( any(is.na(X[,r]))){
        X[,r]<-N%*%rbmf.vector.gibbs(An,cn,t(N)%*%X[,r])
      }
      counter =0
      while(min(X[,1])<0 & counter< n_tries)
      {
        r= 1
        X[,r]<-N%*%rbmf.vector.gibbs(An,cn,t(N)%*%X[,r])
        counter = counter+1
      }
      
      if(min(X[,1])<0){
        X[,1]<- x_bak
      }
    }
    
    return (X)
  }


# B<- matrix(rnorm(100),10)
# B<- (B+t(B))/2 + diag(1,10)*5
# eigen(B)$values
# 
# 
# 
# A<- matrix(rnorm(400),20)
# A<- (A+t(A))/2+ diag(1,20)*6
# eigen(A)$values
# 
# C<- matrix(rnorm(200),20)
# 
# C[,1]<- 1
# X<- qr.Q(qr(C))
# X[,1]<- abs(X[,1])
# 
# X<- r.constrained.bmf.matrix.gibbs(A,B,C,X,1E4)
# X[,1]
