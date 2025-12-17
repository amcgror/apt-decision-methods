library(Rcpp)
library(nloptr)
library('boot')
library('dplyr')
library(doParallel)
library(foreach)
library(lme4)

sim_fun <- function(x){
  
  n_treats_start <- 8 #number of treatments
  min_n <- 50 #min number of patients per arm needed before interim analysis
  max_n <- 200 #max number of patients per arm needed for trial conclusion 
  thresh_sup <- x[1]
  thresh_fut <-x[2]
  mu0 <- 5
  n0 <- 50
  alpha0<-n0/2
  beta0 <- (9*n0)/2
  treats <- c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8")
  recruit<- 50 #recruit rate
  nsamps<-10000
  n0a<-5
  alpha_a <-n0a/2
  beta_a <- (9*n0a)/2
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
                         recruit){
    n_treats <- n_treats_start
    recruit_rate <- matrix(recruit,1, n_treats); colnames(recruit_rate) <- treats
    Result = array(dim=c(2,n_treats), 
                   dimnames = list(c("Sample.size", "Result"),treats))
    ind <- matrix(0,1, n_treats); colnames(ind) <- treats
    
    
    ##simulate prior data
    t<-matrix(NA, nrow=min_n, ncol=n_treats)
    #n_samples <- matrix(0,1,n_treats);colnames(n_samples) <- treats
    
    tau = rgamma(1, alpha0, beta0)
    mu = rnorm(1, mu0, sqrt(1/(tau*n0)))
    for(i in 1:n_treats){
      for(j in 1:min_n){
        t[j,i] = rnorm(1, mu, 1/sqrt(tau))
        #n_samples[i] <- length(t[,i])
      }}
    
    ntot <- min_n
    
    
    
    
    while(TRUE){
      if(ntot >= min_n){
        ntot <- nrow(t)
        post_mu <- matrix(NA, nrow=nsamps, ncol=n_treats)
        n_new <- t(cbind(rep(0, n_treats))); colnames(n_new) <- treats 
        ##get posterior data
        mu_pos <- c(rep(0, n_treats))
        beta_pos <- c(rep(0, n_treats))
        alpha_pos <-c(rep(0, n_treats))
        scale <- c(rep(0, n_treats))
        df <- c(rep(0, n_treats))
        location<-c(rep(0, n_treats))
        #data comes from priors then samples 10000 from the posterior distribution
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
        # ntot<-nrow(t)
        
        #mean_rank
        ranking<-c(1:n_treats)
        
        cppFunction('NumericMatrix rankRows(NumericMatrix post_mu, NumericVector ranking) {
  int nsamps = post_mu.nrow();
  int n_treats = post_mu.ncol();
  NumericMatrix ranks(nsamps, n_treats);
  
  for(int k = 0; k < nsamps; k++) {
    NumericVector x = -post_mu.row(k);
    NumericVector s = clone(x);
    std::sort(s.begin(), s.end());
    IntegerVector index(n_treats);
    NumericVector ranked(n_treats);
    
    for(int i = 0; i < n_treats; i++) {
      
      for(int j = 0; j < n_treats; j++) {
        if(x[j] == s[i]) {
          index[j] = i+1;
          
        }
      }
    }
    
  

    
    ranks.row(k) = index;
  }
  
  return ranks;
}')
        
        ranks <- rankRows(post_mu=post_mu, ranking=ranking)
        mean_rank <- matrix(apply(t(ranks), 1, mean),1, n_treats); colnames(mean_rank) <- treats
        
        for(i in 1:n_treats){
          # n.samples <- as.vector(n_samples[,treats[i]])
          ind.state <- as.vector(ind[,treats[i]])
          
          if(ntot >= min_n){
            if(ind.state == 0){
              mean_rank_ <- as.numeric(mean_rank[,treats[i]])
              
              if(mean_rank_ < thresh_sup){
                Result["Sample.size", treats[i]] <- ntot
                Result["Result", treats[i]] <- "Superiority"
                ind[,treats[i]] = 1
              }
              
              if(mean_rank_ > thresh_fut){
                Result["Sample.size", treats[i]] <- ntot
                Result["Result", treats[i]] <- "Inferiority"
                ind[,treats[i]] = 2
              }
              
              if(mean_rank_ >= thresh_sup && mean_rank_ <= thresh_fut){
                Result["Sample.size", treats[i]] <- ntot
                Result["Result", treats[i]] <- "Nothing"
                ind[,treats[i]] = 0
              }
              
              
            } ###add scenario where if its the last treatment its superior
          }
          if(ntot >= max_n){
            Result["Sample.size", treats[i]] <- ntot
            # Result["Result", treats[i]] <- "Inconclusive"
            ind[,treats[i]] = 1
          }
          
          
          
          
        }
        
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
            n_new[,treats[j]] <- recruit_rate[,treats[j]]}
          n <- n_new + ntot
          
        }else{
          
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
      
      tnew<-matrix(NA, nrow=recruit_rate, ncol=n_treats)
      
      
      for(i in 1:n_treats){
        for(j in 1:recruit){
          
          tnew[j,i] =  rnorm(1, mu, 1/sqrt(tau))
          
        }}
      t<-rbind(t, tnew)
    }
    return(Result)
  }
  
  
  
  
  
  
  #parallel computing code
  n_sim <- 10000
  ncores <- detectCores() - 1
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  Sim_Res <- foreach(k = 1:n_sim, .combine = rbind,
                     .packages = c("boot", "dplyr", "Rcpp")
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
               recruit=recruit) 
  }
  stopCluster(cl)
  
  
  results <- as.data.frame(Sim_Res)
  
  t1<- sum(results=="Superiority")/n_sim
  
  n_treats_start <- 8 #number of treatments
  min_n <- 50 #min number of patients per arm needed before interim analysis
  max_n <- 200 #max number of patients per arm needed for trial conclusion 
  #thresh_sup <-thresh1[m]
  #thresh_fut <-thresh2[m]
  thresh_sup <- x[1]
  thresh_fut <-x[2]
  mu0 <- 5
  #mu1 <- tsup[m]
  mu1<-6
  n0 <- 50
  alpha0<-n0/2
  beta0 <- (9*n0)/2
  treats <- c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8")
  recruit<- 50 #recruit rate
  nsamps<-10000
  sup_treat <- "T8"
  n0a<-5
  alpha_a <-n0a/2
  beta_a <- (9*n0a)/2
  
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
                         recruit){
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
        ntot <- nrow(t)
        post_mu <- matrix(NA, nrow=nsamps, ncol=n_treats)
        n_new <- t(cbind(rep(0, n_treats))); colnames(n_new) <- treats 
        ##get posterior data
        mu_pos <- c(rep(0, n_treats))
        beta_pos <- c(rep(0, n_treats))
        alpha_pos <-c(rep(0, n_treats))
        scale <- c(rep(0, n_treats))
        df <- c(rep(0, n_treats))
        location<-c(rep(0, n_treats))
        #data comes from priors then samples 10000 from the posterior distribution
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
        # ntot<-nrow(t)
        
        
        #mean_rank
        ranking<-c(1:n_treats)
        
        cppFunction('NumericMatrix rankRows(NumericMatrix post_mu, NumericVector ranking) {
  int nsamps = post_mu.nrow();
  int n_treats = post_mu.ncol();
  NumericMatrix ranks(nsamps, n_treats);
  
  for(int k = 0; k < nsamps; k++) {
    NumericVector x = -post_mu.row(k);
    NumericVector s = clone(x);
    std::sort(s.begin(), s.end());
    IntegerVector index(n_treats);
    NumericVector ranked(n_treats);
    
    for(int i = 0; i < n_treats; i++) {
      
      for(int j = 0; j < n_treats; j++) {
        if(x[j] == s[i]) {
          index[j] = i+1;
          
        }
      }
    }
    
  

    
    ranks.row(k) = index;
  }
  
  return ranks;
}')
        
        ranks <- rankRows(post_mu=post_mu, ranking=ranking)
        mean_rank <- matrix(apply(t(ranks), 1, mean),1, n_treats); colnames(mean_rank) <- treats
        
        for(i in 1:n_treats){
          # n.samples <- as.vector(n_samples[,treats[i]])
          ind.state <- as.vector(ind[,treats[i]])
          
          if(ntot >= min_n){
            if(ind.state == 0){
              mean_rank_ <- as.numeric(mean_rank[,treats[i]])
              
              if(mean_rank_ < thresh_sup){
                Result["Sample.size", treats[i]] <- ntot
                Result["Result", treats[i]] <- "Superiority"
                ind[,treats[i]] = 1
              }
              
              if(mean_rank_ > thresh_fut){
                Result["Sample.size", treats[i]] <- ntot
                Result["Result", treats[i]] <- "Inferiority"
                ind[,treats[i]] = 2
              }
              
              if(mean_rank_ >= thresh_sup && mean_rank_ <= thresh_fut){
                Result["Sample.size", treats[i]] <- ntot
                Result["Result", treats[i]] <- "Nothing"
                ind[,treats[i]] = 0
              }
              
              
            } ###add scenario where if its the last treatment its superior
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
                     .packages = c("boot", "dplyr", "Rcpp")
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
               recruit=recruit) 
  }
  stopCluster(cl)
  
  results <- as.data.frame(Sim_Res)
  9
  pow <- sum(results[sup_treat]=="Superiority")/n_sim
  
  if(t1 > 0.05){
    power = 0
  }
  
  if(t1 <= 0.05){
    power=-pow
  }
  
  return(power)
}


p<-c(1.1, 7.9)
out <- Nelder_Mead(sim_fun, p, control=list(maxfun=110))

write.csv(out, file=".csv")

