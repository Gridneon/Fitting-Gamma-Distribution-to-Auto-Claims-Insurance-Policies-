## Take the autoclaim data in package insurance 
## sample 200 observation
## fit the gamma distribution to your sample using the method of moment and mle
## Superimpose the density functions on histogram of your sample
## Comment on your results


## Resample observations of the size n=200
# from your sample( so called Bootstap)
## Find the MOM and MLE
## Repeat this 1000 times
## Find the MSE and Bias for the MOM and MLE estimators







set.seed(222)
sple_claim_amt<-sample(x=AutoClaims$PAID,size =200,replace =T) # sample of 200 claim amt
head(sple_claim_amt)
hist(sple_claim_amt,freq = F,main = NULL,xlab = "sample claim Amount ")   #histogram with density plot 


## using method of moment Gamma(alpha,Beta)
# alpha= mean(x)/variance(x)     # beta= mean^2(x)/variance(x)
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

Mom.est[[1]]  #  moment shape 
Mom.est[[2]]  #  moment rate 


#log likelihood function
loglik<-function(x,par)
{
  Alpha<-par[1];Beta<-par[2]
  llik<-sum(log(((Beta^Alpha)/(gamma(Alpha)))*
                  x^(Alpha-1)*exp(-Beta*x)))
  return(-llik)
} 

MLE.est<-optim(par=c(Mom.est[[1]],Mom.est[[2]]), fn=loglik,
               method ="Nelder-Mead",x=sple_claim_amt)

MLE.est$par[1]  #mle shape 
MLE.est$par[2]  #mle shape


#fiting the moment gamma curve on the sample 
curve(dgamma(x,shape=Mom.est[[1]],rate=Mom.est[[2]], log = FALSE),
      add = TRUE, col = "red",lty=2, lwd = 2)

#fiting the mle gamma curve on the sample 
curve(dgamma(x,shape=MLE.est$par[1],rate=MLE.est$par[2], log = FALSE),
      add = TRUE, col = "darkblue",lty=1, lwd = 2)



legend("topright", legend=c("mom_gamma","mle_gamma"), lty=c(2,1), 
       col = c("red","darkblue"),lwd=2)



R<-1000
n=200


bootstrapdata_1<-sapply(1:R,
                        function(i)sample(x=sple_claim_amt,size=200,replace =T))
head(bootstrapdata_1)
Estimators_1<-function(Data,par)
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

# funtion of bias 
bias.fun_1<-function(theta.hat_1,theta_1)
{
  mean(theta.hat_1)-theta_1
}

################

Alpha.hat_1<-XTY[,c(1,3)]

Beta.hat_1<-XTY[,-c(1,3)]
##################
#computing the bias of the 
Bias.alpha_1<-apply(Alpha.hat_1,2,
                    function(y)bias.fun_1(theta.hat_1=y,theta_1=Mom.est[[1]]))


Bias.beta_1<-apply(Beta.hat_1,2,
                   function(y)bias.fun_1(theta.hat_1=y,theta_1=Mom.est[[2]]))                    
# funtion of bias 
MSE.fun<-function(theta.hat_1,theta_1)
{
  mean((theta.hat_1-theta_1)^2)
}


################## 
#computing the MSE  
MSE.alpha<-apply(Alpha.hat_1,2,
                 function(y)MSE.fun(theta.hat_1=y,theta_1=Mom.est[[1]]))


MSE.beta<-apply(Beta.hat_1,2,
                function(y)MSE.fun(theta.hat_1=y,theta_1=Mom.est[[2]])) 

####results of the MSE and Bias

result<-cbind(Bias.alpha_1,Bias.beta_1,MSE.alpha,MSE.beta)
print(result)