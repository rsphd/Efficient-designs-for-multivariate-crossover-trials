library(MASS)
library(fastmatrix) #For Kronecker Product

p=3
t=3
n=18

rho_11=seq(0.01,0.99,by=0.01)

J=matrix(rep(1,p^2),nrow=p,ncol=p)
one_vec = rep(1,n^2)
J_n = matrix(one_vec,nrow=n,ncol=n)
identity_n = diag(n)
H_n = identity_n - J_n/n
one_vec_t = rep(1,t^2)
J_t = matrix(one_vec_t,nrow=t,ncol=t)
identity_t = diag(t)
H_t = identity_t - J_t/t
psi=matrix(c(0,0,0,1,0,0,0,1,0),nrow=p,ncol=p,byrow = T)

T_d1 = matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3,byrow=TRUE)
T_d2 = matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3,byrow=TRUE)
T_d3 = matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3,byrow=TRUE)
T_d4 = matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3,byrow=TRUE)
T_d5 = matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3,byrow=TRUE)
T_d6 = matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3,byrow=TRUE)
T_d7 = matrix(c(0,0,1,1,0,0,0,1,0),nrow=3,ncol=3,byrow=TRUE)
T_d8 = matrix(c(0,0,1,1,0,0,0,1,0),nrow=3,ncol=3,byrow=TRUE)
T_d9 = matrix(c(0,0,1,1,0,0,0,1,0),nrow=3,ncol=3,byrow=TRUE)
T_d10 = matrix(c(0,0,1,1,0,0,0,1,0),nrow=3,ncol=3,byrow=TRUE)
T_d11 = matrix(c(0,0,1,1,0,0,0,1,0),nrow=3,ncol=3,byrow=TRUE)
T_d12 = matrix(c(0,0,1,1,0,0,0,1,0),nrow=3,ncol=3,byrow=TRUE)
T_d13 = matrix(c(0,1,0,0,0,1,1,0,0),nrow=3,ncol=3,byrow=TRUE)
T_d14 = matrix(c(0,1,0,0,0,1,1,0,0),nrow=3,ncol=3,byrow=TRUE)
T_d15 = matrix(c(0,1,0,0,0,1,1,0,0),nrow=3,ncol=3,byrow=TRUE)
T_d16 = matrix(c(0,1,0,0,0,1,1,0,0),nrow=3,ncol=3,byrow=TRUE)
T_d17 = matrix(c(0,1,0,0,0,1,1,0,0),nrow=3,ncol=3,byrow=TRUE)
T_d18 = matrix(c(0,1,0,0,0,1,1,0,0),nrow=3,ncol=3,byrow=TRUE)

T_d1_orthogonal=matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3,byrow=TRUE)

#For Mat(0.5)

efficiency=c()

for(j in 1:length(rho_11))
{
  V=matrix(c(1,rho_11[j],rho_11[j]^2,rho_11[j],1,rho_11[j],rho_11[j]^2,rho_11[j],1),nrow=3,ncol=3)
  
  V_star=solve(V) - solve(V)%*%J%*%solve(V)/sum(rowSums(solve(V)))
  
  C_11 = t(T_d1)%*%V_star%*%T_d1 + t(T_d2)%*%V_star%*%T_d2 + t(T_d3)%*%V_star%*%T_d3 + t(T_d4)%*%V_star%*%T_d4 + t(T_d5)%*%V_star%*%T_d5 + t(T_d6)%*%V_star%*%T_d6 + t(T_d7)%*%V_star%*%T_d7 + t(T_d8)%*%V_star%*%T_d8 + t(T_d9)%*%V_star%*%T_d9 + t(T_d10)%*%V_star%*%T_d10 + t(T_d11)%*%V_star%*%T_d11 + t(T_d12)%*%V_star%*%T_d12 + t(T_d13)%*%V_star%*%T_d13 + t(T_d14)%*%V_star%*%T_d14 + t(T_d15)%*%V_star%*%T_d15 + t(T_d16)%*%V_star%*%T_d16 + t(T_d17)%*%V_star%*%T_d17 + t(T_d18)%*%V_star%*%T_d18
  C_12 = t(T_d1)%*%V_star%*%psi%*%T_d1 + t(T_d2)%*%V_star%*%psi%*%T_d2 + t(T_d3)%*%V_star%*%psi%*%T_d3 + t(T_d4)%*%V_star%*%psi%*%T_d4 + t(T_d5)%*%V_star%*%psi%*%T_d5 + t(T_d6)%*%V_star%*%psi%*%T_d6 + t(T_d7)%*%V_star%*%psi%*%T_d7 + t(T_d8)%*%V_star%*%psi%*%T_d8 + t(T_d9)%*%V_star%*%psi%*%T_d9 + t(T_d10)%*%V_star%*%psi%*%T_d10 + t(T_d11)%*%V_star%*%psi%*%T_d11 + t(T_d12)%*%V_star%*%psi%*%T_d12 + t(T_d13)%*%V_star%*%psi%*%T_d13 + t(T_d14)%*%V_star%*%psi%*%T_d14 + t(T_d15)%*%V_star%*%psi%*%T_d15 + t(T_d16)%*%V_star%*%psi%*%T_d16 + t(T_d17)%*%V_star%*%psi%*%T_d17 + t(T_d18)%*%V_star%*%psi%*%T_d18
  C_21 = t(C_12)
  C_22 = t(T_d1)%*%t(psi)%*%V_star%*%psi%*%T_d1 + t(T_d2)%*%t(psi)%*%V_star%*%psi%*%T_d2 + t(T_d3)%*%t(psi)%*%V_star%*%psi%*%T_d3 + t(T_d4)%*%t(psi)%*%V_star%*%psi%*%T_d4 + t(T_d5)%*%t(psi)%*%V_star%*%psi%*%T_d5 + t(T_d6)%*%t(psi)%*%V_star%*%psi%*%T_d6 + t(T_d7)%*%t(psi)%*%V_star%*%psi%*%T_d7 + t(T_d8)%*%t(psi)%*%V_star%*%psi%*%T_d8 + t(T_d9)%*%t(psi)%*%V_star%*%psi%*%T_d9 + t(T_d10)%*%t(psi)%*%V_star%*%psi%*%T_d10 + t(T_d11)%*%t(psi)%*%V_star%*%psi%*%T_d11 + t(T_d12)%*%t(psi)%*%V_star%*%psi%*%T_d12 + t(T_d13)%*%t(psi)%*%V_star%*%psi%*%T_d13 + t(T_d14)%*%t(psi)%*%V_star%*%psi%*%T_d14 + t(T_d15)%*%t(psi)%*%V_star%*%psi%*%T_d15 + t(T_d16)%*%t(psi)%*%V_star%*%psi%*%T_d16 + t(T_d17)%*%t(psi)%*%V_star%*%psi%*%T_d17 + t(T_d18)%*%t(psi)%*%V_star%*%psi%*%T_d18
  
  C=C_11 - C_12%*%ginv(C_22)%*%C_21
  
  e_11=sum(diag(t(T_d1_orthogonal)%*%V_star%*%T_d1_orthogonal))
  e_12=sum(diag(t(T_d1_orthogonal)%*%V_star%*%psi%*%T_d1_orthogonal))
  e_22=sum(diag(t(T_d1_orthogonal)%*%t(psi)%*%V_star%*%psi%*%T_d1_orthogonal))-V_star[1,1]/t
  
  E=n/(t-1)*matrix(c(e_11,e_12,e_12,e_22), nrow=2, ncol=2, byrow=TRUE)
  
  C_orthogonal=(det(E)/e_22)*H_t
  
  efficiency[j]=sum(diag(C))/sum(diag(C_orthogonal))
  
  print(j)
}

data=data.frame(efficiency)
names(data)=c("Efficiency")

#Calling library writexl
library(writexl)
#Storing dataframe in given path as excel file
write_xlsx(data,"D:/Code for Manuscript 2 New 2/Mat_half_prop_example_p3t3.xlsx")

#For Mat(1.5)

efficiency=c()

for(j in 1:length(rho_11))
{
  V=matrix(c(1,(1-log(rho_11[j]))*rho_11[j],(1-2*log(rho_11[j]))*rho_11[j]^2,(1-log(rho_11[j]))*rho_11[j],1,(1-log(rho_11[j]))*rho_11[j],(1-2*log(rho_11[j]))*rho_11[j]^2,(1-log(rho_11[j]))*rho_11[j],1),nrow=p,ncol=p)
  
  V_star=solve(V) - solve(V)%*%J%*%solve(V)/sum(rowSums(solve(V)))
  
  C_11 = t(T_d1)%*%V_star%*%T_d1 + t(T_d2)%*%V_star%*%T_d2 + t(T_d3)%*%V_star%*%T_d3 + t(T_d4)%*%V_star%*%T_d4 + t(T_d5)%*%V_star%*%T_d5 + t(T_d6)%*%V_star%*%T_d6 + t(T_d7)%*%V_star%*%T_d7 + t(T_d8)%*%V_star%*%T_d8 + t(T_d9)%*%V_star%*%T_d9 + t(T_d10)%*%V_star%*%T_d10 + t(T_d11)%*%V_star%*%T_d11 + t(T_d12)%*%V_star%*%T_d12 + t(T_d13)%*%V_star%*%T_d13 + t(T_d14)%*%V_star%*%T_d14 + t(T_d15)%*%V_star%*%T_d15 + t(T_d16)%*%V_star%*%T_d16 + t(T_d17)%*%V_star%*%T_d17 + t(T_d18)%*%V_star%*%T_d18
  C_12 = t(T_d1)%*%V_star%*%psi%*%T_d1 + t(T_d2)%*%V_star%*%psi%*%T_d2 + t(T_d3)%*%V_star%*%psi%*%T_d3 + t(T_d4)%*%V_star%*%psi%*%T_d4 + t(T_d5)%*%V_star%*%psi%*%T_d5 + t(T_d6)%*%V_star%*%psi%*%T_d6 + t(T_d7)%*%V_star%*%psi%*%T_d7 + t(T_d8)%*%V_star%*%psi%*%T_d8 + t(T_d9)%*%V_star%*%psi%*%T_d9 + t(T_d10)%*%V_star%*%psi%*%T_d10 + t(T_d11)%*%V_star%*%psi%*%T_d11 + t(T_d12)%*%V_star%*%psi%*%T_d12 + t(T_d13)%*%V_star%*%psi%*%T_d13 + t(T_d14)%*%V_star%*%psi%*%T_d14 + t(T_d15)%*%V_star%*%psi%*%T_d15 + t(T_d16)%*%V_star%*%psi%*%T_d16 + t(T_d17)%*%V_star%*%psi%*%T_d17 + t(T_d18)%*%V_star%*%psi%*%T_d18
  C_21 = t(C_12)
  C_22 = t(T_d1)%*%t(psi)%*%V_star%*%psi%*%T_d1 + t(T_d2)%*%t(psi)%*%V_star%*%psi%*%T_d2 + t(T_d3)%*%t(psi)%*%V_star%*%psi%*%T_d3 + t(T_d4)%*%t(psi)%*%V_star%*%psi%*%T_d4 + t(T_d5)%*%t(psi)%*%V_star%*%psi%*%T_d5 + t(T_d6)%*%t(psi)%*%V_star%*%psi%*%T_d6 + t(T_d7)%*%t(psi)%*%V_star%*%psi%*%T_d7 + t(T_d8)%*%t(psi)%*%V_star%*%psi%*%T_d8 + t(T_d9)%*%t(psi)%*%V_star%*%psi%*%T_d9 + t(T_d10)%*%t(psi)%*%V_star%*%psi%*%T_d10 + t(T_d11)%*%t(psi)%*%V_star%*%psi%*%T_d11 + t(T_d12)%*%t(psi)%*%V_star%*%psi%*%T_d12 + t(T_d13)%*%t(psi)%*%V_star%*%psi%*%T_d13 + t(T_d14)%*%t(psi)%*%V_star%*%psi%*%T_d14 + t(T_d15)%*%t(psi)%*%V_star%*%psi%*%T_d15 + t(T_d16)%*%t(psi)%*%V_star%*%psi%*%T_d16 + t(T_d17)%*%t(psi)%*%V_star%*%psi%*%T_d17 + t(T_d18)%*%t(psi)%*%V_star%*%psi%*%T_d18
  
  C=C_11 - C_12%*%ginv(C_22)%*%C_21
  
  e_11=sum(diag(t(T_d1_orthogonal)%*%V_star%*%T_d1_orthogonal))
  e_12=sum(diag(t(T_d1_orthogonal)%*%V_star%*%psi%*%T_d1_orthogonal))
  e_22=sum(diag(t(T_d1_orthogonal)%*%t(psi)%*%V_star%*%psi%*%T_d1_orthogonal))-V_star[1,1]/t
  
  E=n/(t-1)*matrix(c(e_11,e_12,e_12,e_22), nrow=2, ncol=2, byrow=TRUE)
  
  C_orthogonal=(det(E)/e_22)*H_t
  
  efficiency[j]=sum(diag(C))/sum(diag(C_orthogonal))
  
  print(j)
}

data=data.frame(efficiency)
names(data)=c("Efficiency")

#Calling library writexl
library(writexl)
#Storing dataframe in given path as excel file
write_xlsx(data,"D:/Code for Manuscript 2 New 2/Mat_3half_prop_example_p3t3.xlsx")

#For Mat(infinity)

efficiency=c()

for(j in 1:length(rho_11))
{
  V=matrix(c(1,rho_11[j],rho_11[j]^4,rho_11[j],1,rho_11[j],rho_11[j]^4,rho_11[j],1),nrow=3,ncol=3)
  
  V_star=solve(V) - solve(V)%*%J%*%solve(V)/sum(rowSums(solve(V)))
  
  C_11 = t(T_d1)%*%V_star%*%T_d1 + t(T_d2)%*%V_star%*%T_d2 + t(T_d3)%*%V_star%*%T_d3 + t(T_d4)%*%V_star%*%T_d4 + t(T_d5)%*%V_star%*%T_d5 + t(T_d6)%*%V_star%*%T_d6 + t(T_d7)%*%V_star%*%T_d7 + t(T_d8)%*%V_star%*%T_d8 + t(T_d9)%*%V_star%*%T_d9 + t(T_d10)%*%V_star%*%T_d10 + t(T_d11)%*%V_star%*%T_d11 + t(T_d12)%*%V_star%*%T_d12 + t(T_d13)%*%V_star%*%T_d13 + t(T_d14)%*%V_star%*%T_d14 + t(T_d15)%*%V_star%*%T_d15 + t(T_d16)%*%V_star%*%T_d16 + t(T_d17)%*%V_star%*%T_d17 + t(T_d18)%*%V_star%*%T_d18
  C_12 = t(T_d1)%*%V_star%*%psi%*%T_d1 + t(T_d2)%*%V_star%*%psi%*%T_d2 + t(T_d3)%*%V_star%*%psi%*%T_d3 + t(T_d4)%*%V_star%*%psi%*%T_d4 + t(T_d5)%*%V_star%*%psi%*%T_d5 + t(T_d6)%*%V_star%*%psi%*%T_d6 + t(T_d7)%*%V_star%*%psi%*%T_d7 + t(T_d8)%*%V_star%*%psi%*%T_d8 + t(T_d9)%*%V_star%*%psi%*%T_d9 + t(T_d10)%*%V_star%*%psi%*%T_d10 + t(T_d11)%*%V_star%*%psi%*%T_d11 + t(T_d12)%*%V_star%*%psi%*%T_d12 + t(T_d13)%*%V_star%*%psi%*%T_d13 + t(T_d14)%*%V_star%*%psi%*%T_d14 + t(T_d15)%*%V_star%*%psi%*%T_d15 + t(T_d16)%*%V_star%*%psi%*%T_d16 + t(T_d17)%*%V_star%*%psi%*%T_d17 + t(T_d18)%*%V_star%*%psi%*%T_d18
  C_21 = t(C_12)
  C_22 = t(T_d1)%*%t(psi)%*%V_star%*%psi%*%T_d1 + t(T_d2)%*%t(psi)%*%V_star%*%psi%*%T_d2 + t(T_d3)%*%t(psi)%*%V_star%*%psi%*%T_d3 + t(T_d4)%*%t(psi)%*%V_star%*%psi%*%T_d4 + t(T_d5)%*%t(psi)%*%V_star%*%psi%*%T_d5 + t(T_d6)%*%t(psi)%*%V_star%*%psi%*%T_d6 + t(T_d7)%*%t(psi)%*%V_star%*%psi%*%T_d7 + t(T_d8)%*%t(psi)%*%V_star%*%psi%*%T_d8 + t(T_d9)%*%t(psi)%*%V_star%*%psi%*%T_d9 + t(T_d10)%*%t(psi)%*%V_star%*%psi%*%T_d10 + t(T_d11)%*%t(psi)%*%V_star%*%psi%*%T_d11 + t(T_d12)%*%t(psi)%*%V_star%*%psi%*%T_d12 + t(T_d13)%*%t(psi)%*%V_star%*%psi%*%T_d13 + t(T_d14)%*%t(psi)%*%V_star%*%psi%*%T_d14 + t(T_d15)%*%t(psi)%*%V_star%*%psi%*%T_d15 + t(T_d16)%*%t(psi)%*%V_star%*%psi%*%T_d16 + t(T_d17)%*%t(psi)%*%V_star%*%psi%*%T_d17 + t(T_d18)%*%t(psi)%*%V_star%*%psi%*%T_d18
  
  C=C_11 - C_12%*%ginv(C_22)%*%C_21
  
  e_11=sum(diag(t(T_d1_orthogonal)%*%V_star%*%T_d1_orthogonal))
  e_12=sum(diag(t(T_d1_orthogonal)%*%V_star%*%psi%*%T_d1_orthogonal))
  e_22=sum(diag(t(T_d1_orthogonal)%*%t(psi)%*%V_star%*%psi%*%T_d1_orthogonal))-V_star[1,1]/t
  
  E=n/(t-1)*matrix(c(e_11,e_12,e_12,e_22), nrow=2, ncol=2, byrow=TRUE)
  
  C_orthogonal=(det(E)/e_22)*H_t
  
  efficiency[j]=sum(diag(C))/sum(diag(C_orthogonal))
  
  print(j)
}

data=data.frame(efficiency)
names(data)=c("Efficiency")

#Calling library writexl
library(writexl)
#Storing dataframe in given path as excel file
write_xlsx(data,"D:/Code for Manuscript 2 New 2/Mat_inf_prop_example_p3t3.xlsx")

library(readxl)
Mat_half_prop_ex_p3t3=read_excel("D:/Code for Manuscript 2 New 2/Mat_half_prop_example_p3t3.xlsx")
Mat_3half_prop_ex_p3t3=read_excel("D:/Code for Manuscript 2 New 2/Mat_3half_prop_example_p3t3.xlsx")
Mat_inf_prop_ex_p3t3=read_excel("D:/Code for Manuscript 2 New 2/Mat_inf_prop_example_p3t3.xlsx")

max(Mat_half_prop_ex_p3t3$Efficiency)
max(Mat_3half_prop_ex_p3t3$Efficiency)
max(Mat_inf_prop_ex_p3t3$Efficiency)


