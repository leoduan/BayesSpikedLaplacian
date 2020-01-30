setwd("~/git/BayesLaplacian/")
source("sim1.R")

K = 10
theta = 1
Lam = abs(rnorm(K,0,0.01))
Lam[0] = 0
I_n = diag(1,n)
I_K = diag(1,K)


Omega = diag(1,K)
Omega[1,1]=0

eta = rep(1,K)

sigma_eps = 1E-3
sigma_lambda1 = 0.1**2
sigma_lambda2 = 0.1**2

sigma_theta = 0.1
mu_theta=1

a_sigma_eps = 2
b_sigma_eps = 1E-8

sigma_mu = 1

svdL= svd(L)
uL = svdL$u
dL = svdL$d

idx = rev(c((n-K+1):n))
Q<- uL[,idx]

Q<- Q+rnorm(length(Q))*0.001
Q<- qr.Q(qr(Q))
if(min(Q[,1])<0){
  Q[,1]= -Q[,1]
}

Lam = dL[idx]
theta = mean(dL[-idx])
E = L - Q%*%(Lam-I_K*theta)%*%t(Q)
diag(E)<- 0

sigma_eps = sum(E**2)/2/(n*(n-1)/2)
sigma_eps

require("truncnorm")
source("rcbing.R")

for (i in (1:100)){
  # update Q
  par1 = (I_n * 2.0001 -L)/sigma_eps/2
  par2 = diag(theta- Lam)

  
  for(k in 1:5){
    Q<- tryCatch({
      rcbing(par1, par2 ,Q, n_tries = 1E3)
      # rbing.matrix.gibbs(par1,par2,Q)
    },error = function(e) {
      print("err")
      Q
    },
    finally = {Q}
    )
  }
  
  
  {
    # update Lambda
    Var1 = 1/(1/(2* sigma_eps) + 1/sigma_lambda1)
    M1 = Var1 * (   diag(t(Q)%*% (L)%*%Q)/
                      (2* sigma_eps) )
    
    Var2 = 1/(1/(2* sigma_eps) + 1/sigma_lambda2)
    M2 = Var2 * ( diag(t(Q)%*% (L)%*%Q)/
                    (2* sigma_eps) +  theta/sigma_lambda2)
    
    
    Lam1 = rtruncnorm(K,mean=M1, sd=sqrt(Var1),a = 0, b = 2)
    Lam2 = rtruncnorm(K,mean=M2, sd=sqrt(Var2),a = 0, b = 2)
    # 
    # Lam[1] = 0
    Lam<- eta*Lam1 + (1-eta)*Lam2
    
    Lam[1]<- 1E-5
    prob_eta <- cbind( dtruncnorm(Lam, mean=0, sd=sqrt(sigma_lambda1),a = 0, b = 2),
                       dtruncnorm(Lam, mean=mu_theta, sd=sqrt(sigma_lambda2),a = 0, b = 2))
    eta<- (runif(K)< ((prob_eta/rowSums(prob_eta))[,1]))*1
    # 
    Lam[1]<- 0
    eta[1]<-2
    
    # Update Sigma_lambda1
    sigma_lambda1<- 1/rgamma(1, sum(eta==1)/2+ 2, rate = sum( (Lam[eta==1])**2)/2+ 1E-1)
    sigma_lambda2<- 1/rgamma(1, sum(eta==0)/2+ 2, rate = sum( (Lam[eta==0]-mu_theta)**2)/2+ 1E-1)
    
    
    # Update Theta
    
    Var = 1/ (  (n-K)/2/sigma_eps + 1/sigma_theta)
    M10 = (sum(diag(L)) - sum(diag(t(Q)%*%L%*%Q)))/
      2/sigma_eps + mu_theta/sigma_theta
    M1 = Var * M10
    
    theta = rtruncnorm(1, mean=M1, sd= sqrt(Var), a= 0,b= 2)
  }
  
  # Update L
  M1 = diag(Q%*%diag(Lam-theta)%*%t(Q)) + theta
  Var = 2* sigma_eps
  diag(L) = rnorm(n, M1,sqrt(Var))
  
  # Update sigma_eps
  E = L - Q%*%(diag(Lam)-I_K*theta)%*%t(Q)
  diag(E)<- 0
  sigma_eps = 1/rgamma(1, shape =  n*(n-1)/2/2 ,
                       rate = sum(E**2)/2/2 )
  
  print(sigma_eps)
}

Lam
theta

plot(Q[,2],Q[,3])


