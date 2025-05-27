# log-linear zero-inflated generalized Poisson (ZIGP) INGARCHX* Model with delay parameters 
# rainfall and temperature
set.seed(199933)
library(mvtnorm)
library(HMMpa)
library(readxl)

nosim  = 1
result = NULL
S1	 = matrix(0,nosim,5)
S2	 = matrix(0,nosim,5)
S3	 = matrix(0,nosim,5)
S4	 = matrix(0,nosim,5)
S5	 = matrix(0,nosim,5)
S6   = matrix(0,nosim,5)
S7   = matrix(0,nosim,5)
S8   = matrix(0,nosim,5)
S9   = matrix(0,nosim,5)

nob=780
Mean_resid=matrix(0,1,nob)
Mean_resid_square=matrix(0,1,nob)
Mean_yt_predict=matrix(0,1,nob)
Median_yt_predict=matrix(0,1,nob)
U_yt_predict=matrix(0,1,nob)
L_yt_predict=matrix(0,1,nob)
DIC=matrix(0,nosim,1)

r1	 = NULL
r2	 = NULL
r3     	= NULL
r4     	= NULL
r5     	= NULL
r6     	= NULL
r7     	= NULL
r8     	= NULL
r9     	= NULL
####Generate ZIGP Distribution Function####

#Generate random variable x from ZIGP distribution
cdf0<-function(x,lam,rho)
{
  if(x==0){
    return(rho+(1-rho)*exp(-lam))
  }
}

cdf1<-function(x,lam,psi,rho){
  v<-seq(1,x,by=1)
  return(rho+(1-rho)*exp(-lam)+sum((1-rho)*lam*((lam+psi*v)^(v-1))*exp(-lam-psi*v)/factorial(v)))
}

pois.cdf<-function(x,lam,rho,psi){
  if(x==0){
    sum=cdf0(0,lam,rho)
  }else{
    sum=cdf1(x,lam,psi,rho)
  }
  return(sum)
}

# Creating a function for ZIGP model since not available sa R

r.zigp <- function(n, lam,rho,psi)
{
  U <- runif(n)
  X <- rep(0,n)
  # loop through each uniform
  for(i in 1:n)
  {
    # first check if you are in the first interval
    if(U[i] < pois.cdf(0,lam,rho,psi))
    {
      X[i] <- 0
    } else
    {
      # while loop to determine which subinterval,I, you are in
      # terminated when B = TRUE
      B = FALSE
      I = 0
      while(B == FALSE)
      {
        # the interval to check
        int <- c( pois.cdf(I, lam,rho,psi), pois.cdf(I+1,lam,rho,psi))
        # see if the uniform is in that interval
        if( (U[i] > int[1]) & (U[i] < int[2]) )
        {
          # if so, quit the while loop and store the value
          X[i] <- I+1
          B = TRUE
        } else
        {
          # If not, continue the while loop and increase I by 1
          I=I+1
        }
      }
    }
  }
  return(X)
}


#for the delay parameter (since no available function sa mode in R)
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}



# Likelihood sa ZIGP model
liko<-function(nob,yt,xt1,xt2,alpha0,alpha1,beta,rho,psi,om1,om2,d1,d2,d0){
  liko=0
  
  for (t in (d0+1):nob)
  { 
    
    
    nu[t]=alpha0+ alpha1*log(yt[t-1]+1) + beta*nu[t-1]
    nu_s[t] = ((1-psi)/(1-rho))*nu[t]+om1*xt1[t-d1]+om2*xt2[t-d2]
    lam_s[t]=exp(nu_s[t])
    
    
    if(yt[t]==0){
      liko=liko+log(rho+(1-rho)*exp(-(lam_s[t])))
    }
    else{
      liko= liko+log(1-rho)+log(lam_s[t])+(yt[t]-1)*log(lam_s[t]+psi*yt[t])-(lam_s[t])-psi*yt[t]
    }
    
  }
  return(liko)
}


library(readxl)
iligan_data <- readxl::read_xlsx(path="C:\\Users\\ACER\\Downloads\\ThesisIT!\\IC\\Iligan_Dengue_Cases.xlsx", sheet=1) 



for (isi in 1:nosim)
{
  #yt<- iligan_data$ILIGAN[1:780]  #Data1
  
  #************  (for TS plot labeling)
  # Define the time variable
  #time <- seq(2010, 2024, by = 2)  # Adjust for skipping by 2 years
  
  # Plot the time series with the time labels
  #plot.ts(Yt, xaxt = "n", ylab = "Yt", xlab = "Time")  # Add axis labels
  #axis(1, at = seq(1, nob, length.out = length(time)), labels = time)  # Adjust for skipping by 2 years
  
  # Add title
  #title(main = "Data 1: Dengue Cases in Iligan City", line = 2, col.main = "black", cex.main = 0.8)
  #************  
  
  plot.ts(yt)
  
  d0=3
  
  xt1 <- iligan_data$rainfall
  xt2 <- iligan_data$temperature
  
  
  
  
  
  # plot(ts(xt))
  #****************************************************************New Part
  # *******************Random walk Metropolis
  # Setting up starting values
  
  alpha0	=	0.01
  alpha1	=	0.01
  beta	=	0.01
  om1=0.01
  om2=0.01
  rho=0.01
  psi=0.01
  d1=1
  d2=1
  lam_s=NULL
  nu=NULL
  nu[1:(d0)]	= yt[1]
  nu_o=NULL
  nu_o[1:(d0)]	= yt[1]
  nu_s	=	NULL
  nu_star=NULL
  star=NULL
  star1=NULL
  yt_pred=NULL
  m=NULL
  v=NULL
  res=NULL
  
  count	=0
  count1=0
  count2=0
  count3=0
  count4=0
  
  
  ### set stepsize 
  step_al	   = 	0.1
  step_c     =  0.05
  step_psi   =  0.05
  step_rho   =  0.05
  
  #Hyperparameters
  a_1=3
  a_2=2
  c_1=2
  c_2=1
  
  
  M 	= 	20000
  ind	=   	8000
  draws = matrix(0,M,9)
  d=matrix(0,M,1) #DIC
  resid=matrix(0,M,nob)
  resid_square=matrix(0,M,nob)
  yt_predict=matrix(0,M,nob)
  
  for (i in 1:M){
    
    lik9=0
    lik10=0
    lognor1=0
    lognor2=0
    lognor3=0
    lognor4=0
    lognor5=0
    lognor6=0
    lognor7=0
    lognor8=0
    
    old=c(alpha0,alpha1,beta)
    old1=psi
    old2=c(om1,om2)
    
    
    
    if(i  <= ind)
    {
      repeat 
      {
        star[1]	=  alpha0 +  rnorm(1,0,1)*step_al
        star[2]	=  alpha1 +  rnorm(1,0,1)*step_al
        star[3]	=  beta  +   rnorm(1,0,1)*step_al
        if ( (star[1]>0)&(star[2] >0 )&( star[3] >0) & (star[2]+ star[3]<1) ) {break}
      }
    }else
      repeat
      {
        star=rmvnorm(1,IK_mean,IK_cov)
        if ( (star[1]>0)&(star[2] >0 )&( star[3] >0) & (star[2]+ star[3]<1) ) {break}
        
        lognor2=dmvnorm(star,IK_mean,IK_cov,log=T)
        lognor1=dmvnorm(old,IK_mean,IK_cov,log=T)
      }
    
    lik1<-liko(nob,yt,xt1,xt2,alpha0,alpha1,beta,rho,psi,om1,om2,d1,d2,d0)
    lik2<-liko(nob,yt,xt1,xt2,star[1],star[2],star[3],rho,psi,om1,om2,d1,d2,d0)
    
    lik = lognor1-lognor2+lik2- lik1
    u 	=	 runif(1)
    if (log(u)<lik) 
    {
      alpha0  =  star[1]
      alpha1 	=  star[2]
      beta    =  star[3]
      count	  =  count+1
    }
    
    
    if(i  <= ind)
    { 
      repeat
      {
        psi_star = psi +rnorm(1,0,1)*step_psi
        if((psi_star>=0& psi_star<1)){break} 
      }
    }else
      repeat
      {
        psi_star = rnorm(1,0,1)*step_psi
        if((psi_star>=0& psi_star<1)){break} 
        
        lognor4=dnorm(psi_star,IK_mean1,IK_cov1,log=T)
        lognor3=dnorm(old1,IK_mean1,IK_cov1,log=T)
      }
    
    
    
    
    lik3<-liko(nob,yt,xt1,xt2,alpha0,alpha1,beta,rho,psi,om1,om2,d1,d2,d0)
    lik4<-liko(nob,yt,xt1,xt2,alpha0,alpha1,beta,rho,psi_star,om1,om2,d1,d2,d0)
    
    lik = lognor3-lognor4+lik4-lik3
    u 	=	 runif(1)
    if (log(u)<lik) 
    {
      
      psi 	=  psi_star
      count1	= count1+1
    }
    
    repeat
    {
      rho_star = rho +rnorm(1,0,1)*step_rho
      if((rho_star>0& rho_star<1)){break} 
    }
    
    lik5<-liko(nob,yt,xt1,xt2,alpha0,alpha1,beta,rho,psi,om1,om2,d1,d2,d0)
    lik6<-liko(nob,yt,xt1,xt2,alpha0,alpha1,beta,rho_star,psi,om1,om2,d1,d2,d0)
    
    lik = lik6-lik5
    u 	=	 runif(1)
    if (log(u)<lik) 
    {
      
      rho 	=  rho_star
      count2	= count2+1
    }
    
    if(i  <= ind)
    {
      repeat 
      {
        star1[1]	= om1 + rnorm(1,0,1)*step_c
        star1[2]	= om2 + rnorm(1,0,1)*step_c
        if ( (star1[1]>0)&(star1[2]>0)) {break}
      }
    }else
      repeat
      {
        star1=rmvnorm(1,IK_mean3,IK_cov3)
        if ( (star1[1]>0)&(star1[2]>0)) {break}
        
        lognor6=dmvnorm(star1,IK_mean3,IK_cov3,log=T)
        lognor5=dmvnorm(old2,IK_mean3,IK_cov3,log=T)
      }
    
    lik7<-liko(nob,yt,xt1,xt2,alpha0,alpha1,beta,rho,psi,om1,om2,d1,d2,d0)+(d1-1)*log(om1)-d2*om1+(d1-1)*log(om2)-d2*om2
    lik8<-liko(nob,yt,xt1,xt2,alpha0,alpha1,beta,rho,psi,star1[1],star1[2],d1,d2,d0)+(d1-1)*log(star1[1])-d2*star1[1]+(d1-1)*log(star1[2])-d2*star1[2]
    
    
    lik = lognor5-lognor6+lik8-lik7
    u 	=	 runif(1)
    if (log(u)<lik) 
    {
      om1     =  star1[1]
      om2     =  star1[2]
      count3	= count3+1
    }
    
    
    
    
    lik_b=NULL
    lik_bm=NULL
    sum_lik=0
    prob=NULL
    
    
    for(j in 1:d0){
      
      lik_b[j]=liko(nob,yt,xt1,xt2,alpha0,alpha1,beta,rho,psi,om1,om2,j,d2,d0)
      
    }
    
    max_lik=max(lik_b[1],lik_b[2],lik_b[3])
    
    for(j in 1:d0){
      
      lik_bm[j]=exp(lik_b[j]-max_lik)
      
      sum_lik=sum_lik+lik_bm[j]
      
    }
    
    for(j in 1:d0){
      
      prob[j]=lik_bm[j]/sum_lik
      
    }
    
    cum_prob=cumsum(prob)
    
    U <- runif(1)
    d1 <- rep(0,1)
    
    if(U < cum_prob[1])
    {
      d1 <- 1
    } else
    {
      # while loop to determine which subinterval,I, you are in
      # terminated when B = TRUE
      B = FALSE
      I = 1
      while(B == FALSE)
      {
        # the interval to check
        int <- c(cum_prob[I], cum_prob[I+1])
        # see if the uniform is in that interval
        if( (U > int[1]) & (U < int[2]) )
        {
          # if so, quit the while loop and store the value
          d1 <- I+1
          B = TRUE
        } else
        {
          # If not, continue the while loop and increase I by 1
          I=I+1
        }
      }
    }
    
    lik_b2=NULL
    lik_bm2=NULL
    sum_lik2=0
    prob2=NULL
    
    
    for(j in 1:d0){
      
      lik_b2[j]=liko(nob,yt,xt1,xt2,alpha0,alpha1,beta,rho,psi,om1,om2,d1,j,d0)
      
    }
    
    max_lik=max(lik_b2[1],lik_b2[2],lik_b2[3])
    
    for(j in 1:d0){
      
      lik_bm2[j]=exp(lik_b2[j]-max_lik)
      
      sum_lik2=sum_lik2+lik_bm2[j]
      
    }
    
    for(j in 1:d0){
      
      prob2[j]=lik_bm2[j]/sum_lik2
      
    }
    
    cum_prob=cumsum(prob2)
    
    U <- runif(1)
    d2 <- rep(0,1)
    
    if(U < cum_prob[1])
    {
      d2 <- 1
    } else
    {
      # while loop to determine which subinterval,I, you are in
      # terminated when B = TRUE
      B = FALSE
      I = 1
      while(B == FALSE)
      {
        # the interval to check
        int <- c(cum_prob[I], cum_prob[I+1])
        # see if the uniform is in that interval
        if( (U > int[1]) & (U < int[2]) )
        {
          # if so, quit the while loop and store the value
          d2 <- I+1
          B = TRUE
        } else
        {
          # If not, continue the while loop and increase I by 1
          I=I+1
        }
      } 
    }
    
    draws[i,] = c(alpha0,alpha1,beta,psi,rho,om1,om2,d1,d2)
    
    
    m= NULL
    v=NULL
    res=NULL
    yt_pred= NULL
    
    for (t in (d0+1):nob)
    {
      nu_o[t]=alpha0+ alpha1*log(yt[t-1]+1) + beta*nu_o[t-1]
      nu_s[t] = ((1-psi)/(1-rho))*nu_o[t]+om1*xt1[t-d1]+om2*xt2[t-d2]
      lam_s[t]=exp(nu_s[t])
      
      m[t]= (1 - rho)/((1 - psi) * lam_s[t]) #ZIGP mean formula
      v[t]=(1-rho)*((rho*lam_s[t]^2)/(1-psi)^2+(lam_s[t]/(1-psi)^3))
      
      res[t] = (yt[t]-m[t])/sqrt(v[t])
    }
    
    resid[i,]=res
    resid_square[i,]=res^2
    
    for (t in (d0+1):nob)
    {
      
      nu_o[t]=alpha0+ alpha1*log(yt[t-1]+1) + beta*nu_o[t-1]
      nu_s[t] = ((1-psi)/(1-rho))*nu_o[t]+om1*xt1[t-d1]+om2*xt2[t-d2]
      lam_s[t]=exp(nu_s[t])
      
      yt_pred[t] =r.zigp(1, lam_s[t],rho,psi)
      
    }
    yt_predict[i,]=yt_pred
    
    
    for (t in (d0+1):nob)
    { 
      nu_o[t]=alpha0+ alpha1*log(yt[t-1]+1) + beta*nu_o[t-1]
      nu_s[t] = ((1-psi)/(1-rho))*nu_o[t]+om1*xt1[t-d1]+om2*xt2[t-d2]
      lam_s[t]=exp(nu_s[t])
      
      if(yt[t]==0){
        lik9=lik9-2*log(rho+(1-rho)*exp(-(lam_s[t])))
      }
      else{
        lik9= lik9-2*log(1-rho)-2*log(lam_s[t])-2*(yt[t]-1)*log(lam_s[t]+psi*yt[t])+2*(lam_s[t])+2*psi*yt[t]+2*lgamma(yt[t]+1)
      }
      
    }
    #d=matrix(0,M,length(lik9))
    d[i,]=c(lik9)
    
    
    
    if(i%%2000 == 0)
    {cat("************  SIMU AT: ",isi,"\n")
      cat(i,"\n")
      cat("alpha0, alpha1,beta",alpha0, alpha1,beta,"\n")         
      cat("accept. rate",100*count /i,"\n")
      cat("psi",psi, "\n")         
      cat("accept. rate",100*count1 /i,"\n")
      cat("rho",rho, "\n")         
      cat("accept. rate",100*count2 /i,"\n")
      cat("c1,c2",om1,om2, "\n")         
      cat("accept. rate",100*count3 /i,"\n")
      cat("d1","d2",d1,d2, "\n")         
      # cat("accept. rate",100*count4 /i,"\n")
    }
    if(i==ind)
    {
      #IK_mean=c(mean(draws[1001:ind,1]),mean(draws[1001:ind,2]),mean(draws[1001:ind,3]))
      #IK_cov=cov(draws[1001:ind,1:3])
      #IK_mean1=c(mean(draws[1001:ind,6]),mean(draws[1001:ind,7]))
      #IK_cov1=cov(draws[1001:ind,6:7])
      
      IK_mean=c(mean(draws[1001:ind,1]),mean(draws[1001:ind,2]),mean(draws[1001:ind,3]))
      IK_cov=cov(draws[1001:ind,1:3])
      IK_mean1=mean(draws[1001:ind,4])
      IK_cov1=var(draws[1001:ind,4])
      #IK_mean2=mean(draws[1001:ind,5])
      #IK_cov2=var(draws[1001:ind,5])
      IK_mean3=c(mean(draws[1001:ind,6]),mean(draws[1001:ind,7]))
      IK_cov3=cov(draws[1001:ind,6:7])
      
      
    }
  }
  
  MCMC=(ind+1):M
  alpha0_mean=mean(draws[MCMC,1])
  alpha1_mean=mean(draws[MCMC,2])
  beta_mean=mean(draws[MCMC,3])
  psi_mean=mean(draws[MCMC,4])
  rho_mean=mean(draws[MCMC,5])
  om1_mean=mean(draws[MCMC,6])
  om2_mean=mean(draws[MCMC,7])
  d1_mode=getmode(draws[MCMC,8])
  d2_mode=getmode(draws[MCMC,9])
  
  
  for (t in (d0+1):nob)
  { 
    nu_o[t]=alpha0_mean+ alpha1_mean*log(yt[t-1]+1) + beta_mean*nu_o[t-1]
    nu_s[t] = ((1-psi_mean)/(1-rho_mean))*nu_o[t]+om1_mean*xt1[t-d1_mode]+om2_mean*xt2[t-d2_mode]
    lam_s[t]=exp(nu_s[t])
    if(yt[t]==0){
      lik10=lik10-2*log(rho_mean+(1-rho_mean)*exp(-(lam_s[t])))
    }
    else{
      lik10= lik10-2*log(1-rho_mean)-2*log(lam_s[t])-2*(yt[t]-1)*log(lam_s[t]+psi_mean*yt[t])+2*(lam_s[t])+2*psi_mean*yt[t]+2*lgamma(yt[t]+1)
    } 
  }
  
  DIC[isi, 1] = 2*(mean(d[MCMC,1]))-lik10
  
  
  ############################# PLOT
  names		 = c(list(expression(alpha0)),list(expression(alpha1)),list(expression(beta)), list(expression(psi)),list(expression(rho)),list(bquote(omega[0])),list(bquote(omega[1])),"d1","d2")
  if (isi==1) 
  {
    par(mai=c(0.6,0.4,0.5,0.4),mfrow=c(4,4))
    for (i in 1:7){
      ts.plot(draws[8001:M,i],xlab="iterations",ylab="",main=names[i])
      acf(draws[8001:M,i],main=names[i])
    }
  }
  ######################################################################
  #############################
  S1[isi,1]=mean(draws[MCMC,1])
  S1[isi,2]=median(draws[MCMC,1])
  S1[isi,3]=sd(draws[MCMC,1])
  S1[isi,4]=quantile(draws[MCMC,1],0.025)
  S1[isi,5]=quantile(draws[MCMC,1],0.975)
  
  S2[isi,1]=mean(draws[MCMC,2])
  S2[isi,2]=median(draws[MCMC,2])
  S2[isi,3]=sd(draws[MCMC,2])
  S2[isi,4]=quantile(draws[MCMC,2],0.025)
  S2[isi,5]=quantile(draws[MCMC,2],0.975)
  
  S3[isi,1]=mean(draws[MCMC,3])
  S3[isi,2]=median(draws[MCMC,3])
  S3[isi,3]=sd(draws[MCMC,3])
  S3[isi,4]=quantile(draws[MCMC,3],0.025)
  S3[isi,5]=quantile(draws[MCMC,3],0.975)
  
  S4[isi,1]=mean(draws[MCMC,4])
  S4[isi,2]=median(draws[MCMC,4])
  S4[isi,3]=sd(draws[MCMC,4])
  S4[isi,4]=quantile(draws[MCMC,4],0.025)
  S4[isi,5]=quantile(draws[MCMC,4],0.975)
  
  S5[isi,1]=mean(draws[MCMC,5])
  S5[isi,2]=median(draws[MCMC,5])
  S5[isi,3]=sd(draws[MCMC,5])
  S5[isi,4]=quantile(draws[MCMC,5],0.025)
  S5[isi,5]=quantile(draws[MCMC,5],0.975)
  
  S6[isi,1]=mean(draws[MCMC,6])
  S6[isi,2]=median(draws[MCMC,6])
  S6[isi,3]=sd(draws[MCMC,6])
  S6[isi,4]=quantile(draws[MCMC,6],0.025)
  S6[isi,5]=quantile(draws[MCMC,6],0.975)
  
  S7[isi,1]=mean(draws[MCMC,7])
  S7[isi,2]=median(draws[MCMC,7])
  S7[isi,3]=sd(draws[MCMC,7])
  S7[isi,4]=quantile(draws[MCMC,7],0.025)
  S7[isi,5]=quantile(draws[MCMC,7],0.975)
  
  S8[isi,1]=getmode(draws[MCMC,8])
  S9[isi,1]=getmode(draws[MCMC,9])
  
  Mean_resid=apply(resid[MCMC,-1:-3],2,mean)
  Mean_resid_square=apply(resid_square[MCMC,-1:-3],2,mean)
  Mean_yt_predict=apply(yt_predict[MCMC,-1:-3],2,mean)
  
  
  #_______________________________________________________________
#To save in excel
  #  MCMCfile<-cbind(S1,S2,S3,S4,S5,S6,S7,S8,S9)
  
 # colnames(MCMCfile)<-c("mean_alpha0","median_alpha0","std_alpha0","P025_alpha0","P975_alpha0",
 #                       "mean_alpha1","median_alpha1","std_alpha1","P025_alpha1","P975_alpha1",
  #                      "mean_beta","median_beta","std_beta","P025_beta","P975_beta",
   #                     "mean_psi","median_psi","std_psi","P025_psi","P975_psi",
    #                    "mean_rho","median_rho","std_rho","P025_rho","P975_rho",
     #                   "mean_om1","median_om1","std_om1","P025_om1","P975_om1",
      #                  "mean_om2","median_om2","std_om2","P025_om2","P975_om2",
       #                 "mode_d1","NAN","NAN","NAN","NAN",
        #                "mode_d2","NAN","NAN","NAN","NAN")
  ### save si2.out
  #write.csv(MCMCfile,"C:/Users/ACER/Downloads/ZI GROUP CODES/2cov_old_real.csv")
#______________________________________________________________________________________

  }

r1=apply(S1,2,mean)
r2=apply(S2,2,mean)
r3=apply(S3,2,mean)
r4=apply(S4,2,mean)
r5=apply(S5,2,mean)
r6=apply(S6,2,mean)
r7=apply(S7,2,mean)
r8=apply(S8,2,mean)
r9=apply(S9,2,mean)
result=round(rbind(r1,r2,r3,r4,r5,r6,r7,r8,r9),4)
colnames(result) <- c("mean","median","std","P025","P975")
rownames(result) <- names
result
acceptrate=count/M
acceptrate1=count1/M
acceptrate2=count2/M
acceptrate3=count3/M

acceptrate
acceptrate1
acceptrate2
acceptrate3


D=apply(DIC,2,mean)
#D


#---------------For Posterior Odds Ratio use

# Estimate log marginal likelihood using harmonic mean in log-space
#harmonic_mean_log <- function(logL_vector) {
#  max_logL <- max(logL_vector)
#  return(max_logL - log(mean(exp(logL_vector - max_logL))))
#}

#logL_raintemp <- -0.5 * d[(ind+1):M, 1]  # Convert -2*logL to logL
#mlog_raintemp <- harmonic_mean_log(logL_rain&temp)  # Save log marginal likelihood

# Show the result
#cat("Log Marginal Likelihood (Rainfall & Temp model):", mlog_raintemp, "\n")
#--------------------------------------


#par(mai=c(0.8,0.4,0.6,0.4),mfrow=c(1,3))
#plot(ts(Mean_resid),main="Standardized Residual Plot")
#abline(h=0,col="red")
#acf(Mean_resid,main="ACF of Standardized Residual")
#acf(Mean_resid_square,main="ACF of Standardized Squared Residual")



library(tidyr)
library(dplyr)
library(ggplot2)

newdata <- data.frame(yt[-1:-3], Mean_yt_predict)
newdata$year = seq(from = 1, to = nrow(newdata))
names(newdata) <- c("Observed", "Predicted", "year")

plot = newdata %>% 
  pivot_longer(cols = 'Observed':'Predicted', names_to = "var", values_to = "val") %>%
  ggplot(aes(x = year, y = val, color = var, group = var, linetype = var)) +
  geom_line(size = 0.75) +
  scale_color_manual(name = "", values = c("blue", "red")) +
  scale_linetype_manual(name = "", values = c(2, 1)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab("") +
  ylab("") +
  ggtitle("Data 1: Iligan - Rainfall & Temperature") +
  theme(plot.title = element_text(hjust = 0.5,size = 20, face="bold",margin=margin(0,0,30,0)))+
  theme(axis.text.x  = element_text(size=19,colour = "black"))+
  theme(axis.text.y  = element_text(size=19 ,colour="black"))+
  theme(axis.line = element_line(size=0.5, colour = "black"))+
  theme(legend.title =  element_blank())+
  theme(legend.text = element_text(size=16))+
  theme(legend.position = c(0.1,0.9))+
  scale_x_continuous(breaks = c(1,53,105,157,209,261,313,365,417,469,521,573,625,677,729), #PLUS 53
                     labels = c("2010", "2011", "2012", "2013", "2014", "2015","2016","2017","2018","2019","2020","2021", "2022","2023","2024"))


print(plot)