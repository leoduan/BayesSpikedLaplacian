rm(list=ls())
ns = c(10, 20, 30)
ns_cumsum = cumsum(ns)
n = sum(ns)

A = matrix(0,n,n)

for (i in 1:length(ns_cumsum)){
  if(i ==1){
    idx1 = 1
  }else{
    idx1 = ns_cumsum[i-1]+1
  }
  idx2 = ns_cumsum[i]
  
  n0 = idx2 -idx1+1
  
  mat = matrix(runif(n0**2, 0,1),n0)
  mat = (t(mat)+mat)/2
  #     mat = (mat<0.5)
  A[idx1:idx2,idx1:idx2]= mat

}

A0 = A
# add noise

mat = matrix(runif(n**2, 0,0.01),n,n)
mat = (t(mat)+mat)/2

A = A + mat
A = abs(A)

diag(A)<- 0
# Add articulation points

A[9,10] = A[10,9] =1
A[29,30] = A[30,29] =1

# Add noisy articulation points

A[9,12] = A[12,9] =1
A[29,35] = A[35,29] =1

# Plot the adjacency

heatmap(A)


# plotmat(A)

# Degree

D = diag(rowSums(A))
Dinv = solve(D)
DinvHalf = Dinv**(0.5)


# Normalized Laplacian

L = DinvHalf%*% (D-A) %*% DinvHalf

