require("rstiefel")

rcbing <-
  function(A,B,X, n_tries=1E4)
  {
    #simulate from the matrix bmf distribution as described in Hoff(2009) 
    #this is one Gibbs step, and must be used iteratively
    
    ### assumes B is a diagonal matrix with *decreasing* entries 
    x_bak = X[,1]
    
    
    m<-dim(X)[1] ;  R<-dim(X)[2]
    if(m>R)
    {
      for(r in sample( c(2:R)))
      {
        N<-NullC(X[,-r])
        An<-B[r,r]*t(N)%*%(A)%*%N 
        X[,r]<-N%*%rbing.vector.gibbs(An,t(N)%*%X[,r])
      }
    }
    
    #If m=R then the fc of one vector given all the others is 
    #just +-1 times the vector in the null space. In this case, 
    #the matrix needs to be updated at least two columns at a 
    #time. 
    if(m==R)
    {
      for(s in c(2:R))
      {
        r<-sort( sample(c(2:R),2) )
        N<-NullC( X[,-r]  )
        An<- t(N)%*%A%*%N
        #X[,r]<-N%*%rbing.O2(An,B[r,r]) 
        X[,r]<-N%*%rbing.Op(An,B[r,r]) 
      }
    }
    
    {
      r = 1
      N<-NullC(X[,-r])
      An<-B[r,r]*t(N)%*%(A)%*%N 
      X[,r]<-N%*%rbing.vector.gibbs(An,t(N)%*%X[,r])
      
      while( any(is.na(X[,r]))){
        X[,r]<-N%*%rbing.vector.gibbs(An,t(N)%*%X[,r])
      }
      
      counter =0
      while(min(X[,1])<0 & counter< n_tries)
      {
        r= 1
        X[,r]<-N%*%rbing.vector.gibbs(An,t(N)%*%X[,r])
        counter = counter+1
      }
      
      if(min(X[,1])<0){
        X[,1]<- x_bak
      }
    }
    
    
    X
  }
  
  
