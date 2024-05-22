#V1:Mat(infinity) VR: Mat(1.5)
library(MASS)
library(raster)

p=3
t=3
n=6
identity_p=diag(p)
J=matrix(rep(1,p^2),nrow=p,ncol=p)
V=matrix(c(0,0,0,1,0,0,0,1,0),nrow=p,ncol=p,byrow = T)
H=identity_p-J/t
one_vec = rep(1,n^2)
J_n = matrix(one_vec,nrow=n,ncol=n)
identity_n = diag(n)
H_n = identity_n - J_n/n
one_vec_t = rep(1,t^2)
J_t = matrix(one_vec_t,nrow=t,ncol=t)
identity_t = diag(t)
H_t = identity_t - J_t/t

rho=c(seq(-0.999,-0.001,by=0.001),seq(0.001, 0.999, by = 0.001))
rho_star=c()
sigma_12=c()

tr_c=c()
tr_c_orthogonal=c()
tr_c_balanceduniform=c()
tr_c_uniform=c()

rho_11 = c(seq(0.01, 0.99,by = 0.01))

rd1=c()
rd2=c()

min_relative_diff1=c()
min_relative_diff2=c()

max_relative_diff1=c()
max_relative_diff2=c()

for(j in 1:length(rho_11))
{
  for (i in 1:length(rho))
  {
    sigma_11 = 5 #Any positive value is allowed
    sigma_22= sigma_11
    
    V1=sigma_11*matrix(c(1,rho_11[j],rho_11[j]^4,rho_11[j],1,rho_11[j],rho_11[j]^4,rho_11[j],1),nrow=p,ncol=p)
    #V1 is V_1
    
    V2=matrix(c(1,(1-log(rho_11[j]))*rho_11[j],(1-2*log(rho_11[j]))*rho_11[j]^2,(1-log(rho_11[j]))*rho_11[j],1,(1-log(rho_11[j]))*rho_11[j],(1-2*log(rho_11[j]))*rho_11[j]^2,(1-log(rho_11[j]))*rho_11[j],1),nrow=p,ncol=p)
    #V2 is V_R
    
    rho_star[i] = rho[i]*sqrt(sigma_22/sigma_11)
    
    V1_star=solve(V1) - solve(V1)%*%J%*%solve(V1)/sum(rowSums(solve(V1)))
    
    V2_star=solve(V2) - solve(V2)%*%J%*%solve(V2)/sum(rowSums(solve(V2)))
    
    sigma_12[i] = sigma_11*(1-rho[i]^2)
    a1 = sum(diag(V1_star)) + ((1+rho_star[i]^2)/sigma_12[i])*sum(diag(V2_star))
    b1 = sum(diag(V1_star%*%V)) + ((1+rho_star[i]^2)/sigma_12[i])*sum(diag(V2_star%*%V))
    c1 = sum(diag(H%*%t(V)%*%V1_star%*%V)) + ((1+rho_star[i]^2)/sigma_12[i])*sum(diag(H%*%t(V)%*%V2_star%*%V))
    tr_c[i] = n*(a1-((b1^2)/c1)) #Upperbound of trace
    
    S1=V1_star + rho_star[i]^2/sigma_12[i]*V2_star
    S2=rho_star[i]/sigma_12[i]*V2_star
    S3=t(S2)
    S4=1/sigma_12[i]*V2_star
    
    #For relative difference for orthogonal array
    rd1[i] = n*((1+rho_star[i]^2) *(sum(diag(V1_star%*%V))*sum(diag(H%*%t(V)%*%V2_star%*%V))-sum(diag(V2_star%*%V))*sum(diag(H%*%t(V)%*%V1_star%*%V)))^2/((sigma_12[i]*sum(diag(H%*%t(V)%*%V1_star%*%V)) + (1+rho_star[i]^2)*sum(diag(H%*%t(V)%*%V2_star%*%V)))*sum(diag(H%*%t(V)%*%V1_star%*%V))*sum(diag(H%*%t(V)%*%V2_star%*%V))))/tr_c[i] *100
    
    T_d1 = matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3,byrow=T)
    T_d2 = matrix(c(0,1,0,0,0,1,1,0,0),nrow=3,ncol=3,byrow=T)
    T_d3 = matrix(c(0,0,1,1,0,0,0,1,0),nrow=3,ncol=3,byrow=T)
    T_d4 = matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3,byrow=T)
    T_d5 = matrix(c(0,1,0,0,0,1,1,0,0),nrow=3,ncol=3,byrow=T)
    T_d6 = matrix(c(0,0,1,1,0,0,0,1,0),nrow=3,ncol=3,byrow=T)
    
    C_11_11 = t(T_d1)%*%S1%*%T_d1 + t(T_d2)%*%S1%*%T_d2 + t(T_d3)%*%S1%*%T_d3 + t(T_d4)%*%S1%*%T_d4 + t(T_d5)%*%S1%*%T_d5 + t(T_d6)%*%S1%*%T_d6 
    C_11_12 = -(t(T_d1)%*%S2%*%T_d1 + t(T_d2)%*%S2%*%T_d2 + t(T_d3)%*%S2%*%T_d3 + t(T_d4)%*%S2%*%T_d4 + t(T_d5)%*%S2%*%T_d5 + t(T_d6)%*%S2%*%T_d6 )
    C_11_21 = t(C_11_12)
    C_11_22 = t(T_d1)%*%S4%*%T_d1 + t(T_d2)%*%S4%*%T_d2 + t(T_d3)%*%S4%*%T_d3 + t(T_d4)%*%S4%*%T_d4 + t(T_d5)%*%S4%*%T_d5 + t(T_d6)%*%S4%*%T_d6
    
    C_12_11 = (t(T_d1)%*%S1%*%V%*%T_d1 + t(T_d2)%*%S1%*%V%*%T_d2 + t(T_d3)%*%S1%*%V%*%T_d3 + t(T_d4)%*%S1%*%V%*%T_d4 + t(T_d5)%*%S1%*%V%*%T_d5 + t(T_d6)%*%S1%*%V%*%T_d6 )%*%H
    C_12_12 = -(t(T_d1)%*%S2%*%V%*%T_d1 + t(T_d2)%*%S2%*%V%*%T_d2 + t(T_d3)%*%S2%*%V%*%T_d3 + t(T_d4)%*%S2%*%V%*%T_d4 + t(T_d5)%*%S2%*%V%*%T_d5 + t(T_d6)%*%S2%*%V%*%T_d6 )%*%H
    C_12_21 = -(t(T_d1)%*%S3%*%V%*%T_d1 + t(T_d2)%*%S3%*%V%*%T_d2 + t(T_d3)%*%S3%*%V%*%T_d3 + t(T_d4)%*%S3%*%V%*%T_d4 + t(T_d5)%*%S3%*%V%*%T_d5 + t(T_d6)%*%S3%*%V%*%T_d6 )%*%H
    C_12_22 = (t(T_d1)%*%S4%*%V%*%T_d1 + t(T_d2)%*%S4%*%V%*%T_d2 + t(T_d3)%*%S4%*%V%*%T_d3 + t(T_d4)%*%S4%*%V%*%T_d4 + t(T_d5)%*%S4%*%V%*%T_d5 + t(T_d6)%*%S4%*%V%*%T_d6 )%*%H
    
    C_22_11 = H%*%(t(T_d1)%*%t(V)%*%S1%*%V%*%T_d1 + t(T_d2)%*%t(V)%*%S1%*%V%*%T_d2 + t(T_d3)%*%t(V)%*%S1%*%V%*%T_d3 + t(T_d4)%*%t(V)%*%S1%*%V%*%T_d4 + t(T_d5)%*%t(V)%*%S1%*%V%*%T_d5 + t(T_d6)%*%t(V)%*%S1%*%V%*%T_d6 )%*%H
    C_22_12 = -(H%*%(t(T_d1)%*%t(V)%*%S2%*%V%*%T_d1 + t(T_d2)%*%t(V)%*%S2%*%V%*%T_d2 + t(T_d3)%*%t(V)%*%S2%*%V%*%T_d3 + t(T_d4)%*%t(V)%*%S2%*%V%*%T_d4 + t(T_d5)%*%t(V)%*%S2%*%V%*%T_d5 + t(T_d6)%*%t(V)%*%S2%*%V%*%T_d6 )%*%H)
    C_22_21 = t(C_22_12)
    C_22_22 = H%*%(t(T_d1)%*%t(V)%*%S4%*%V%*%T_d1 + t(T_d2)%*%t(V)%*%S4%*%V%*%T_d2 + t(T_d3)%*%t(V)%*%S4%*%V%*%T_d3 + t(T_d4)%*%t(V)%*%S4%*%V%*%T_d4 + t(T_d5)%*%t(V)%*%S4%*%V%*%T_d5 + t(T_d6)%*%t(V)%*%S4%*%V%*%T_d6 )%*%H
    
    C_11 = rbind(cbind(C_11_11,C_11_12), cbind(C_11_21,C_11_22))
    C_12 = rbind(cbind(C_12_11,C_12_12), cbind(C_12_21,C_12_22))
    C_21 = t(C_12)
    C_22 = rbind(cbind(C_22_11,C_22_12), cbind(C_22_21,C_22_22))
    
    C = C_11 - C_12 %*% ginv(C_22) %*% C_21
    tr_c_uniform[i] = sum(diag(C))
    
    #For relative difference for uniform design
    rd2[i] = (tr_c[i] - tr_c_uniform[i])/tr_c[i]*100
  }
  
  min_relative_diff1[j]=min(rd1)
  min_relative_diff2[j]=min(rd2)
  
  max_relative_diff1[j]=max(rd1)
  max_relative_diff2[j]=max(rd2)
  
  print(j)  
}


data=data.frame(min_relative_diff1,max_relative_diff1,min_relative_diff2,max_relative_diff2)
names(data)=c("Minimum_Orthogonal","Maximum_Orthogonal","Minimum_Uniform","Maximum_Uniform")

#Calling library writexl
library(writexl)
#Storing dataframe in given path as excel file
write_xlsx(data,"D:/Code for Manuscript 2 New 2/Mat_inf_Mat_3half_p3t3.xlsx")

#Calling library readxl
library(readxl)

Mat_inf_Mat_3half_p3t3=read_excel("D:/Code for Manuscript 2 New 2/Mat_inf_Mat_3half_p3t3.xlsx")

#Plot for Minimum

plot(rho_11,Mat_inf_Mat_3half_p3t3$Minimum_Orthogonal,ylim = range(c(Mat_inf_Mat_3half_p3t3$Minimum_Orthogonal,100)),xlab=expression('r'[(1)]), ylab="Minimum relative difference (in %)",type="l",lty="solid",lwd=3)
par(new=TRUE)
lines(rho_11,Mat_inf_Mat_3half_p3t3$Minimum_Uniform,lty="dotted",lwd=3)

legend(x = "right", lty = c("solid","dotted"), lwd=c(3,3), text.font = 4, col= c("black"),text.col = "black", legend=c(expression('d'^"*"),expression('d'[1])))

#Plot for Maximum

plot(rho_11,Mat_inf_Mat_3half_p3t3$Maximum_Orthogonal,ylim = range(c(Mat_inf_Mat_3half_p3t3$Maximum_Orthogonal,100)),xlab=expression('r'[(1)]), ylab="Maximum relative difference (in %)",type="l",lty="solid",lwd=3)
par(new=TRUE)
lines(rho_11,Mat_inf_Mat_3half_p3t3$Maximum_Uniform,lty="dotted",lwd=3)

legend(x = "right", lty = c("solid","dotted"), lwd=c(3,3), text.font = 4, col= c("black"),text.col = "black", legend=c(expression('d'^"*"),expression('d'[1])))



