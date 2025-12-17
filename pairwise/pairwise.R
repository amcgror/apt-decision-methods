library('boot')
library('dplyr')
library(doParallel)
library(foreach)

tsup <-c(6,7,8)
T1 <-c(rep(0.0497,3))
thresh1 <- c(rep(0.969609375,3))
thresh2 <- c(rep(0.0695703125,3))
pow <- c(rep(0, 3))
ess <- c(rep(0,3))

#CHANGE
name <- ".csv"

for(m in 1:3){
  
  n_treats_start <- 3 #number of treatments
  min_n <- 50 #min number of patients per arm needed before interim analysis
  max_n <- 200 #max number of patients per arm needed for trial conclusion 
  thresh_sup <-thresh1[m]
  thresh_fut <-thresh2[m]
  mu0 <- 5
  mu1 <- tsup[m]
  n0 <- 50
  alpha0<-n0/2
  beta0 <- (9*n0)/2
  treats <- c("T1", "T2", "T3")
  recruit<- 50 #recruit rate
  nsamps<-10000
  sup_treat <- "T3"
  n0a<-5
  alpha_a <-n0a/2
  beta_a <- (9*n0a)/2
  delta <- 0
  
simulation <- function(n_treats,
                       min_n,
                       max_n,
                       thresh_sup,
                       thresh_fut,
                       mu0,
                       alpha0,
                       beta0,
                       n0, 
                       treats,
                       recruit,
                       delta){
  n_treats <- n_treats_start
  recruit_rate <- matrix(recruit,1, n_treats); colnames(recruit_rate) <- treats
  Result = array(dim=c(2,n_treats), 
                 dimnames = list(c("Sample.size", "Result"),treats))
  ind <- matrix(0,1, n_treats); colnames(ind) <- treats
  sup_ind<-0
  
  
  ##simulate prior data
  t<-matrix(NA, nrow=min_n, ncol=n_treats-1)
  #n_samples <- matrix(0,1,n_treats);colnames(n_samples) <- treats
  
  
  tau = rgamma(1, alpha0, beta0)
  mu = rnorm(1, mu0, sqrt(1/(tau*n0)))
  for(i in 1:n_treats-1){
    for(j in 1:min_n){
      t[j,i] = rnorm(1, mu, 1/sqrt(tau))
      #n_samples[i] <- length(t[,i])
    }}
  mu_1<- rnorm(1, mu1, sqrt(1/(tau*n0)))
  t<-cbind(t, rnorm(min_n, mu_1, 1/sqrt(tau)))
  ntot <- min_n
  
  
  
  
  while(TRUE){
    if(ntot >= min_n){
      ntot<-nrow(t)
      post_mu <- matrix(NA, nrow=nsamps, ncol=n_treats)
      n_new <- t(cbind(rep(0, n_treats))); colnames(n_new) <- treats 
      ##get posterior data
      mu_pos <- c(rep(0, n_treats))
      beta_pos <- c(rep(0, n_treats))
      alpha_pos <-c(rep(0, n_treats))
      df <- c(rep(0, n_treats))
      scale <- c(rep(0, n_treats))
      location<-c(rep(0, n_treats))
      #data comes from priors then samples 10000 from the posterior distribution
      #https://en.wikipedia.org/wiki/Conjugate_prior
      for(i in 1:n_treats){
        mu_pos[i]=((n0a*mu0)+(ntot*mean(t[,i])))/(n0a+ntot)
        beta_pos[i]=beta_a+.5*(sum((t[,i]-mean(t[,i]))^2))+
          (ntot*n0a)/(ntot+n0a)*.5*((mean(t[,i])-mu0)^2)
        alpha_pos[i] = alpha_a+ntot/2
        n_pos= n0a+ntot 
        df[i]=2*alpha_pos[i]
        location[i]=mu_pos[i]
        scale[i]=sqrt(beta_pos[i]/(n_pos*alpha_pos[i]))
        post_mu[,i] = rt(nsamps, df[i])*scale[i]+mu_pos[i]           
      }
  
      
    
      ranking<-c(1:n_treats)
      
      #PAIR-WISE COMPARISON
      
      #Get pairwise probability of  treatmenti > treatmentj
      A <- matrix(0, ncol=n_treats, nrow=n_treats)
      
      for(i in 1:n_treats){
        for(j in 1:n_treats){
          A[i,j] = sum(post_mu[,i] > post_mu[,j]+delta)/nsamps
        }
      }
      
      #compare it to matrix of 0.95 
      B <-matrix(data = c(rep(thresh_sup, n_treats^2)), ncol=n_treats, nrow=n_treats)
      C <- matrix(data = c(rep(thresh_fut, n_treats^2)), ncol=n_treats, nrow=n_treats)
      
      #check for superiority
      D <- A >= B
      
      sup <- matrix(apply(D, 1, sum),1, n_treats); colnames(sup) <- treats
      
      
      #check for inferiority
      E <- A < C
      
      infe <-matrix(apply(E, 1, sum),1, n_treats); colnames(infe) <- treats
      
      ##################################
      for(i in 1:n_treats){
        # n.samples <- as.vector(n_samples[,treats[i]])
        ind.state <- as.vector(ind[,treats[i]])
        
        if(ntot >= min_n){
          if(ind.state == 0){
            sup.trt <- as.vector(sup[,treats[i]])
            inf.trt <-as.vector(infe[,treats[i]]) 
            
            if(sup.trt == n_treats-1){
              Result["Sample.size", treats[i]] <- ntot
              Result["Result", treats[i]] <- "Superiority"
              ind[,treats[i]] = 1
            }
            
            if(inf.trt == n_treats){
              Result["Sample.size", treats[i]] <- ntot
              Result["Result", treats[i]] <- "Inferiority"
              ind[,treats[i]] = 2
            }
            
            if(sup.trt < n_treats-1 && inf.trt < n_treats){
              Result["Sample.size", treats[i]] <- ntot
              Result["Result", treats[i]] <- "Nothing"
              ind[,treats[i]] = 0
            }
          }
        }
        if(ntot >= max_n){
          Result["Sample.size", treats[i]] <- ntot
          # Result["Result", treats[i]] <- "Inconclusive"
          ind[,treats[i]] = 1
        }
        
        
        
        
      }
      
      if(sup_treat %in% ind && ind[,sup_treat]==2){
        sup_ind<-1
      }
      #if(sum(ind==1) >0){break}
      
      
      
      if(sum(ind==2) > 0){
        index <- which(ind==2)
        n_treats <- n_treats - length(index)
        treats <- treats[-index]
        recruit_rate <- recruit_rate[-index]
        recruit_rate <-matrix(recruit_rate,1,n_treats);colnames(recruit_rate) <- treats
        n_new <- n_new[-index]
        n_new <-matrix(n_new,1,n_treats);colnames(n_new) <- treats
        ind <- ind[,-index]
        ind <- matrix(ind,1,n_treats);colnames(ind) <- treats
        t <- t[,-index]
        
        # post.mu <- post.mu[,-index]
        
        
        for(j in 1:n_treats){
          n_new[,treats[j]] <- recruit_rate[,treats[j]]
        }
        n <- n_new + ntot
        
      }
      if(sum(ind==2) == 0){
        
        n_new <- t(cbind(rep(0, n_treats))); colnames(n_new) <- treats
        for(j in 1:n_treats){
          n_new[,treats[j]] <- recruit_rate[,treats[j]]}
        n <- n_new + ntot
        
      }
    }
    
    
    n_new <- n-min_n
    
    if(length(n) == 1){
      Result["Sample.size", treats] <- ntot
      Result["Result", treats] <- "Superiority"
      
      {break}
    }
    if(sum(ind==1) >0){break}
    
    if(sup_ind==0){
      tnew<-matrix(NA, nrow=recruit_rate, ncol=n_treats-1)
      
      
      for(i in 1:n_treats-1){
        for(j in 1:recruit){
          
          tnew[j,i] =  rnorm(1, mu, 1/sqrt(tau))
          
        }}
      tnew<-cbind(tnew, rnorm(recruit_rate[1], mu_1, 1/sqrt(tau)))
      t<-rbind(t, tnew)
    }
    
    if(sup_ind==1){
      tnew<-matrix(NA, nrow=recruit_rate, ncol=n_treats)
      
      
      for(i in 1:n_treats){
        for(j in 1:recruit){
          
          tnew[j,i] =  rnorm(1, mu, 1/sqrt(tau))
          
        }}
      t<-rbind(t, tnew)
    }
  }
  return(Result)
}




#parallel computing code
n_sim <- 10000
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
registerDoParallel(cl)
Sim_Res <- foreach(k = 1:n_sim, .combine = rbind,
                   .packages = c("boot", "dplyr")
) %dopar% { 
  set.seed(k + 123)
  simulation(n_treats=n_treats,
             min_n=min_n,
             max_n=max_n,
             thresh_sup=thresh_sup,
             thresh_fut=thresh_fut,
             mu0=mu0,
             alpha0=alpha0,
             beta0=beta0,
             n0=n0,
             treats=treats,
             recruit=recruit,
             delta=delta) 
}
stopCluster(cl)


results <- as.data.frame(Sim_Res)

#CHANGE 
pow[m] <- sum(results[sup_treat]=="Superiority")/n_sim

nums <- seq(from=1, to = 20000, by=2)
ess[m]<-mean(as.numeric(unlist(results[nums,])))

}

df <-data.frame("tsup"=tsup, "T1"=T1, "Thresh_sup"=thresh1, "Thresh_fut"=thresh2,
                "Power"=pow, "Ess"=ess)

write.csv(df, name)
