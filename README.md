# Fitting-Gamma-Distribution-to-Auto-Claims-Insurance-Policies-
To estimate the MOM /MLE estimators  of the Gamma Distribution fit, determing the best in fitting the auto-claim amount using their MSE's and Bias 

set.seed(222)

sple_claim_amt<-sample(x=AutoClaims$PAID,size =200,replace =T) # sample of 200 claim amt

head(sple_claim_amt)

# Distribution of the Auto-claim
hist(sple_claim_amt,freq = F,main = NULL,xlab = "sample claim Amount ")   #histogram with density plot

# Image 1
<img width="409" height="300" alt="image" src="https://github.com/user-attachments/assets/1b122e7e-d49d-4c60-9ea5-dcfc9fc6cec3" />


## Method of moment Gamma(alpha,Beta)
mom.gamma<-function(x)
{
  mu<-mean(x)
  varx<-var(x)
  
  #solve for shape and rate 
  Alpha<-mu^2/varx
  Beta<-mu/varx
  
  return(c(alpha=Alpha,beta=Beta)) 
}

Mom.est<-mom.gamma(sple_claim_amt)

Mom.est[[1]]  #  moment shape;
Mom.est[[2]]  #  moment rate 

## Maximum Likelihood Gamma(alpha,Beta)
loglik<-function(x,par)                              
{
  Alpha<-par[1];Beta<-par[2]
  llik<-sum(log(((Beta^Alpha)/(gamma(Alpha)))*
                  x^(Alpha-1)*exp(-Beta*x)))
  return(-llik)
}  #log likelihood function

MLE.est<-optim(par=c(Mom.est[[1]],Mom.est[[2]]), fn=loglik,
               method ="Nelder-Mead",x=sple_claim_amt)       #optimization 

MLE.est$par[1]  #mle shape; MLE.est$par[2]                   #mle shape




# Fiting the mom  and mle gamma curve on the sample 

curve(dgamma(x,shape=Mom.est[[1]],rate=Mom.est[[2]], log = FALSE),
      add = TRUE, col = "red",lty=2, lwd = 2)                             # plot the curve of the mom

curve(dgamma(x,shape=MLE.est$par[1],rate=MLE.est$par[2], log = FALSE),
      add = TRUE, col = "darkblue",lty=1, lwd = 2)                        #plot the curve of the mle 

legend("topright", legend=c("mom_gamma","mle_gamma"), lty=c(2,1), 
       col = c("red","darkblue"),lwd=2)                                    # legend for the plot 
# Image 2
  <img width="445" height="356" alt="image" src="https://github.com/user-attachments/assets/f2e5762f-da8b-44b2-9ab1-aa65807bfec1" />     


# Bootstrap
R<-1000   #number of iteration 
n=200     #sample size

bootstrapdata_1<-sapply(1:R,                                                    #bootstrap sampling 
                 function(i)sample(x=sple_claim_amt,size=200,replace =T)) 
head(bootstrapdata_1)
Estimators_1<-function(Data,par)   #function for estimates 
{
  Mom.est<-mom.gamma(Data)
  MLE.est<-optim(par=c(Mom.est[[1]],Mom.est[[2]]), fn=loglik,
                 method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
                            "Brent"),x=Data)
  Tab<-c(MOM=Mom.est,MLE=MLE.est$par)
  names(Tab)<-c("MOM.alpha","MOM.beta",
                "MLE.alpha","MLE.beta")
  return(Tab)
}


XTY<-data.frame(t(sapply(1:R,function(i)Estimators_1(bootstrapdata_1[,i],
                 c(Mom.est[[1]],Mom.est[[2]])))))

head(XTY)

# Funtion of Bias 
bias.fun_1<-function(theta.hat_1,theta_1)
{
  mean(theta.hat_1)-theta_1
}

Alpha.hat_1<-XTY[,c(1,3)]

Beta.hat_1<-XTY[,-c(1,3)]
##################
# Computing the Bias 
Bias.alpha_1<-apply(Alpha.hat_1,2,
                  function(y)bias.fun_1(theta.hat_1=y,theta_1=Mom.est[[1]]))


Bias.beta_1<-apply(Beta.hat_1,2,
                 function(y)bias.fun_1(theta.hat_1=y,theta_1=Mom.est[[2]]))                    
# Funtion of MSE
MSE.fun<-function(theta.hat_1,theta_1)
{
  mean((theta.hat_1-theta_1)^2)
}


################## 
# Computing the MSE  
MSE.alpha<-apply(Alpha.hat_1,2,
                 function(y)MSE.fun(theta.hat_1=y,theta_1=Mom.est[[1]]))


MSE.beta<-apply(Beta.hat_1,2,
                function(y)MSE.fun(theta.hat_1=y,theta_1=Mom.est[[2]])) 

# Results of the MSE and Bias

result<-cbind(Bias.alpha_1,Bias.beta_1,MSE.alpha,MSE.beta)
print(result)
