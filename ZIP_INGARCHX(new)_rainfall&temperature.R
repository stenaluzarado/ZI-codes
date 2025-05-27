# zero inflated Poisson (ZIP) INGARCHX* Model with delay parameters - rainfall & temperature
set.seed(199933)
library(mvtnorm)
library(HMMpa)
nosim  = 1 #no. of simulation
result = NULL #nagcreate og object nga ddto ibutang ang result
S1	 = matrix(0,nosim,5) #create matrix para sa tanan results
S2	 = matrix(0,nosim,5) #S1,S2,... kay kung pila kabook parameters
S3	 = matrix(0,nosim,5) 
S4	 = matrix(0,nosim,5)
S5	 = matrix(0,nosim,5)
S6	 = matrix(0,nosim,5)
S7	 = matrix(0,nosim,5)
S8	 = matrix(0,nosim,5)

nob=780
Mean_resid= matrix(0,1,nob)
Mean_resid_square= matrix(0,1,nob)
Mean_yt_predict= matrix(0,1,nob)
U_yt_predict=matrix(0,1,nob) #upper 95% cf
L_yt_predict=matrix(0,1,nob) #lower 95% cf
DIC=matrix(0,nosim,1)

r1	 = NULL
r2	 = NULL
r3   = NULL
r4   = NULL
r5   = NULL
r6   = NULL
r7   = NULL
r8   = NULL


# getting mode (para sa delay parameter)
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# likelihood zero-inflated Poisson
liko <- function(nob, yt, xt1,xt2, alpha0, alpha1, beta, rho, omega1, omega2, d1,d2, d0) {
  liko = 0
  lam <- numeric(nob)
  lam_s <- numeric(nob)
  
  for (t in (d0 + 1):nob) {
    lam[t] = alpha0 + alpha1 * yt[t - 1] + beta * lam[t - 1]
    lam_s[t] = lam[t] + omega1 * xt1[t - d1] + omega2 * xt2[t - d2]
    liko = liko + (yt[t]) * log(lam_s[t]) - (lam_s[t])
    
    if(yt[t]==0){
      liko=liko+log(rho+(1-rho)*exp(-(lam_s[t])))
    }
    else{
      liko= liko+log(1-rho)+log(lam_s[t])+(yt[t])
    }
    
  }
  
  return(liko)
}


library(readxl)
iligan_data <- readxl::read_xlsx(path="C:\\Users\\ACER\\Downloads\\dummy.xlsx") 



for (isi in 1:nosim)
{
  #yt<- iligan_data$ILIGAN[1:780]    #Data1
  
  plot.ts(yt)
  
  d0=3
  
  xt1 <- iligan_data$rainfall
  xt2 <- iligan_data$temperature
  
  # *********Random Walk Metropolis
  # Setting up starting values
  
  alpha0 = 0.01
  alpha1 = 0.01
  beta = 0.01
  rho = 0.01
  omega1 = 0.01
  omega2 = 0.01
  d1 = 1
  d2=1
  lam = NULL
  lam[1:(d0+1)] = yt[1]
  lam_s = NULL
  lam_s[1:(d0+1)] = yt[1]
  star = NULL
  star1 =NULL
  yt_pred=NULL
  m=NULL # mean
  v=NULL # variance
  res=NULL # residual
  
  
  count = 0 #ang mag count kung nadawat ba o dli sa acceptance rate
  count1 = 0
  count2 = 0
  count3 = 0
  count4 = 0
  
  
  ### step size
  step_a0=0.08
  step_al = 0.1
  step_c = 0.05
  step_rho=0.05
  
  # hyperparameters
  c_1 = 2
  c_2 =1
  
  M 	= 	20000
  ind	=   	8000
  draws = matrix(0,M,8)
  d=matrix(0,M,1)
  resid=matrix(0,M,nob)
  resid_square=matrix(0,M,nob)
  yt_predict=matrix(0,M,nob)
  
  
  ###########
  
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
    old1=c(omega1,omega2)
    
    
    if(i  <= ind)
    {
      repeat 
      {
        star[1]	= alpha0 +  rnorm(1,0,1)*step_a0
        star[2]	=  alpha1 +  rnorm(1,0,1)*step_al
        star[3]		=  beta  +   rnorm(1,0,1)*step_al
        if(star[1]>0& star[2]>0& star[3]>0& (star[2]+star[3]) <1) {break}
      }
    }else
      repeat
      {
        star=rmvnorm(1,IK_mean,IK_cov)
        if(star[1]>0& star[2]>0& star[3]>0& (star[2]+star[3]) <1) {break}
        
        lognor2=dmvnorm(star,IK_mean,IK_cov,log=T)
        lognor1=dmvnorm(old,IK_mean,IK_cov,log=T)
      }
    
    lik1<-liko(nob, yt, xt1,xt2, alpha0, alpha1, beta, rho, omega1, omega2, d1, d2,d0)
    lik2<-liko(nob,yt,xt1,xt2,star[1],star[2],star[3], rho, omega1, omega2, d1, d2, d0)
    
    lik = lognor1-lognor2+lik2- lik1
    u 	=	 runif(1)
    if (log(u)<lik) 
    {
      alpha0  =  star[1]
      alpha1 	=  star[2]
      beta    =  star[3]
      count	= count+1
    }
    
    repeat
    {
      rho_star = rho +rnorm(1,0,1)*step_rho
      if((rho_star>0& rho_star<1)){break} 
    }
    
    lik5<-liko(nob,yt,xt1,xt2,alpha0,alpha1,beta,rho,omega1,omega2, d1, d2,d0)
    lik6<-liko(nob,yt,xt1,xt2,alpha0,alpha1,beta,rho_star,omega1,omega2, d1, d2,d0)
    
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
        star1[1]	= omega1 + rnorm(1,0,1)*step_c
        star1[2]	= omega2 + rnorm(1,0,1)*step_c
        if ( (star1[1]>0)&(star1[2]>0)) {break}
      }
    }else
      repeat
      {
        star1=rmvnorm(1,IK_mean1,IK_cov1)
        if ( (star1[1]>0)&(star1[2]>0)) {break}
        
        lognor4=dmvnorm(star1,IK_mean1,IK_cov1,log=T)
        lognor3=dmvnorm(old1,IK_mean1,IK_cov1,log=T)
      }
    
    lik7<-liko(nob,yt,xt1,xt2,alpha0,alpha1,beta,rho,omega1,omega2,d1,d2,d0)+(c_1-1)*log(omega1)- c_2*omega1 +(c_1-1)*log(omega2)- c_2*omega2
    lik8<-liko(nob,yt,xt1,xt2,alpha0,alpha1,beta,rho,star1[1],star1[2],d1,d2,d0)+(c_1-1)*log(star1[1]) - c_2*star1[1] +(c_1-1)*log(star1[2]) - c_2*star1[2]
    
    
    lik = lognor3-lognor4+lik8-lik7
    u 	=	 runif(1)
    if (log(u)<lik) 
    {
      omega1     =  star1[1]
      omega2     =  star1[2]
      count3	= count3+1
    }
    
    
    lik_b=NULL
    lik_bm=NULL
    sum_lik=0
    prob=NULL
    
    
    for(j in 1:d0){
      
      lik_b[j]=liko(nob,yt,xt1,xt2,alpha0,alpha1,beta,rho,omega1,omega2,j,d2,d0)
      
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
      
      lik_b2[j]=liko(nob,yt,xt1,xt2,alpha0,alpha1,beta,rho,omega1,omega2,d1,j,d0)
      
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
    
    
    draws[i,] = c(alpha0, alpha1, beta, rho, omega1, omega2, d1, d2)
    
    
    m= NULL
    v=NULL
    res=NULL
    yt_pred= NULL
    
    ################################################################################    
    
    
    for (t in (d0+1):nob)
    {
      lam[t] = alpha0 + alpha1 * yt[t - 1] + beta * lam[t - 1]
      lam_s[t] = lam[t] + omega1 * xt1[t - d1] + omega2 * xt2[t - d2]
      
      m[t]= (1-rho)*lam_s[t]  #ZIP mean formula
      v[t]= (1 - rho) * lam_s[t] * (1 + rho * lam_s[t]) #ZIP variance formula
      res[t]=(yt[t]-m[t])/sqrt(v[t])
      
    }
    
    resid[i,]=res
    resid_square[i,]=res^2
    
    
    for (t in (d0+1):nob)
    {
      lam[t] = alpha0 + alpha1 * yt[t - 1] + beta * lam[t - 1]
      lam_s[t] = lam[t] + omega1 * xt1[t - d1] + omega2 * xt2[t - d2]
      
      yt_pred[t] = rpois(1, lambda = lam_s[t])
      
    }
    yt_predict[i,]=yt_pred
    
    
    
    for(t in (d0+1):nob)
    {
      lam[t] = alpha0 + alpha1 * yt[t - 1] + beta * lam[t - 1]
      lam_s[t] = lam[t] + omega1 * xt1[t - d1] + omega2 * xt2[t - d2]
      
      
      if(yt[t]==0){
        lik9=lik9-2*log(rho+(1-rho)*exp(-(lam_s[t])))
      }
      else{
        lik9= lik9-2*log(1-rho)- 2*(yt[t]) * log(lam_s[t]) + 2*(lam_s[t])+ 2*lgamma(yt[t]+1)
      }
      
    }
    
    d[i,] = c(lik9)
    
    if(i%%2000 == 0)
    {cat("************  SIMU AT: ",isi,"\n")
      cat(i,"\n")
      cat("alpha0, alpha1,beta",alpha0, alpha1,beta,"\n")         
      cat("accept. rate",100*count /i,"\n")
      cat("rho",rho, "\n")         
      cat("accept. rate",100*count2 /i,"\n")
      cat("omega1, omega2", omega1,omega2, "\n")         
      cat("accept. rate",100*count3 /i,"\n")
      cat("d1","d2",d1,d2, "\n")         
      # cat("accept. rate",100*count4 /i,"\n")
      
    }
    if(i==ind)
    {
      IK_mean=c(mean(draws[1001:ind,1]),mean(draws[1001:ind,2]),mean(draws[1001:ind,3])) 
      IK_cov=cov(draws[1001:ind,1:3])       #alpha_0, alpha_1, beta1
      IK_mean1=c(mean(draws[1001:ind,5]),mean(draws[1001:ind,6]))
      IK_cov1=cov(draws[1001:ind,5:6])      #omega1, omega2
      
      
    }
  }
  
  # MCMC 
  MCMC= (ind+1):M
  alpha0_mean =mean(draws[MCMC,1])
  alpha1_mean =mean(draws[MCMC,2])
  beta_mean =mean(draws[MCMC,3])
  rho_mean =mean(draws[MCMC,4])
  omega1_mean=mean(draws[MCMC,5])
  omega2_mean=mean(draws[MCMC,6])
  d1_mode =getmode(draws[MCMC,7])
  d2_mode =getmode(draws[MCMC,8])
  
  
  for(t in (d0+1):nob)
  {
    lam[t] = alpha0_mean + alpha1_mean * yt[t - 1] + beta_mean * lam[t - 1]
    lam_s[t] = lam[t] + omega1_mean * xt1[t - d1_mode] + omega2_mean * xt2[t - d2_mode]
    
    
    if(yt[t]==0){
      lik10=lik10-2*log(rho_mean+(1-rho_mean)*exp(-(lam_s[t])))
    }
    else{
      lik10= lik10-2*log(1-rho)- 2*(yt[t]) * log(lam_s[t]) + 2*(lam_s[t])+ 2*lgamma(yt[t]+1)
    } 
    
  }
  
  DIC[isi,1]=2*(mean(d[MCMC,1]))-lik10
  
  
  
  ############################# PLOT
  names		 = c(list(bquote(alpha[0])),list(bquote(alpha[1])),list(bquote(beta[1])), list(bquote(rho)),list(bquote(omega[1])),list(bquote(omega[2])), "d1", "d2")
  if (isi==1) 
  {
    par(mai=c(0.6,0.4,0.5,0.4),mfrow=c(4,4))
    for (i in 1:6){
      ts.plot(draws[8001:M,i],xlab="iterations",ylab="",main=names[i])
      acf(draws[8001:M,i],main=names[i])
      
    }
  }
  ######################################################################
  #############################
  MCMC=(ind+1):M
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
  
  S7[isi,1]=getmode(draws[MCMC,7])
  S8[isi,1]=getmode(draws[MCMC,8])
  
  Mean_resid = apply(resid[MCMC,-1:-3],2,mean)
  Mean_resid_square = apply(resid_square[MCMC,-1:-3],2,mean)
  Mean_yt_predict = apply(yt_predict[MCMC,-1:-3],2,mean)
  
}

# To save in excel
#_________________________________________________________
#MCMCfile<-cbind(S1,S2,S3,S4,S5)

#colnames(MCMCfile)<-c("mean_alpha0","median_alpha0","std_alpha0","P025_alpha0","P975_alpha0",
#                     "mean_alpha1","median_alpha1","std_alpha1","P025_alpha1","P975_alpha1",
#                    "mean_beta","median_beta","std_beta","P025_beta","P975_beta", 
#                   "mean_omega1","median_omega1","std_omega1","P025_omega1","P975_omega1",
#                  "mode_d1","NAN","NAN","NAN","NAN")
### save si2.out
#write.csv(MCMCfile,"C:\\Users\\63960\\Desktop\\POISSON\\p- 1 covariate\\cdo-rainfall.csv")
# write.csv(MCMCfile,"C:/Users/jolia/Downloads/ZIGP_model1_00_model2_rainfall_temp.r")
#_____________________________________________________________________


r1=apply(S1,2,mean)
r2=apply(S2,2,mean)
r3=apply(S3,2,mean)
r4=apply(S4,2,mean)
r5=apply(S5,2,mean)
r6=apply(S6,2,mean)
r7=apply(S7,2,mean)
r8=apply(S8,2,mean)

result=round(rbind(r1,r2,r3,r4,r5,r6,r7,r8),4)
colnames(result) <- c("mean","median","std","P025","P975")
rownames(result) <- names
result
acceptrate=count/M
acceptrate2=count2/M
acceptrate3=count3/M

acceptrate
acceptrate2
acceptrate3


#############################################################################

D=apply(DIC,2,mean)
D


#---------------For Posterior Odds Ratio use

# Estimate log marginal likelihood using harmonic mean in log-space
#harmonic_mean_log <- function(logL_vector) {
#  max_logL <- max(logL_vector)
#  return(max_logL - log(mean(exp(logL_vector - max_logL))))
#}

#logL_raintemp <- -0.5 * d[(ind+1):M, 1]  # Convert -2*logL to logL
#mlog_raintemp <- harmonic_mean_log(logL_raintemp)  # Save log marginal likelihood

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
  scale_x_continuous(breaks = c(1,53,105,157,209,261,313,365,417,469,521,573,625,677,729), #PLUS 52
                     labels = c("2010", "2011", "2012", "2013", "2014", "2015","2016","2017","2018","2019","2020","2021", "2022","2023","2024"))


print(plot)