#### Figure1 (moran'I weight plot) ####
library(dplyr)
pm2.5_dat <- read.csv(file = "C:\\Users\\USER\\Desktop\\논문 코드 공유\\pm2.5_dat2.csv" ,header=T,fileEncoding = 'euc-kr')
name <- pm2.5_dat[,4]
s <- pm2.5_dat[,1:2]
s1 <- cbind(s[,1]*88.8,s[,2]*111)
pm2.5_dat <- pm2.5_dat[,5:ncol(pm2.5_dat)]
weights_mat <- as.matrix(dist(s1))
weights <- (weights_mat < 50 & weights_mat > 0)

library(igraph)
w <- graph_from_adjacency_matrix(
  weights
)

set.seed(69) # 69
plot(w, vertex.size=0.01,edge.arrow.size=0.01,
     vertex.label=name,vertex.label.cex=0.8)

#### moran'I p-value ####
library(ape)
moran <- c()
for(i in 1:ncol(pm2.5_dat)){
  moran[i] <- Moran.I(pm2.5_dat[,i], weights, scaled = FALSE, na.rm = FALSE,
                      alternative = "greater")$p.value
}
moran


#### Figure2 (time plot) ####
pm2.5_dat %>% as.matrix %>% apply(.,2,mean) %>% 
  plot(main='mean',xlab='Time(month)',ylab='value',type='l')
abline(v=c(9,21,33),col='red',lty=2)
legend('topright',legend='1월',lty=2,col='red')
pm2.5_dat %>% as.matrix %>% apply(.,2,sd) %>% 
  plot(main='sd',xlab='Time(month)',ylab='value',type='l')
abline(v=c(9,21,33),col='red',lty=2)
legend('topright',legend='1월',lty=2,col='red')
acf(pm2.5_dat %>% as.matrix %>% apply(.,2,mean),main='acf of mean')
acf(pm2.5_dat %>% as.matrix %>% apply(.,2,sd),main='acf of sd')

#### Figure3 (simulation) ####
library(nimble)
library(dplyr)
library(truncnorm)
# GP kernel
expcov <- nimbleFunction(
  run = function(dists_abs = double(2), dists_squ = double(2),
                 c0 = double(0), c1 = double(0),
                 c2 = double(0), c3 = double(0)) {
    returnType(double(2))
    n1 <- dim(dists_abs)[1]
    n2 <- dim(dists_abs)[2]
    result <- matrix(nrow = n1, ncol = n2, init = FALSE)
    for(i in 1:n1)
      for(j in 1:n2)
        result[i, j] <- c0*c1^dists_squ[i,j]+
      c2*c3^(sin(pi*dists_abs[i,j]/12))^2
    for(i in 1:n1)
      result[i, i] <- result[i, i] + 1e-6
    return(result)
  })
cExpcov <- compileNimble(expcov)

X=1
nmiss=0
loglik <- matrix(0,X,9)
correct_rate <- matrix(0,X,9)
for(x in 1:X){
  N=50
  Time = 9
  K=10
  J=5
  s=matrix(0,N,Time)
  for(t in 1:Time){
    s[,t] <- seq(0.4,20,length=N)
  }
  y=matrix(0,N,Time)
  for(t in 1:Time){
    y[,t][s[,t] > 5+5/8*(t-1) & s[,t] < 10+5/8*(t-1)] <- rnorm(sum(s[,t] > 5+5/8*(t-1) & s[,t] < 10+5/8*(t-1)),mean=20,sd=1)
    y[,t][s[,t] < 5+5/8*(t-1) | s[,t] > 10+5/8*(t-1)] <- rnorm(sum(s[,t] < 5+5/8*(t-1) | s[,t] > 10+5/8*(t-1)),mean=1,sd=1)
  }
  
  Time <- ncol(y)
  
  # model define
  dists_abs <- as.matrix(dist(c(1:(Time+nmiss))))
  dists_squ <- as.matrix(dist(c(1:(Time+nmiss))))^2
  s_tilde <- seq(4,20,length=J)
  code <- nimbleCode({
    for(k in 1:K){
      # mu ~ GP
      mu_0[k] ~ dnorm(mean=0, var=1000)
      c0_mu[k] ~ dinvgamma(shape=1e-3,scale=1e-3)
      c1_mu[k] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
      c2_mu[k] ~ dinvgamma(shape=1e-3,scale=1e-3)
      c3_mu[k] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
      cov_GP_mu[k,1:Time, 1:Time] <- expcov(dists_abs=dists_abs[1:Time,1:Time],
                                            dists_squ=dists_squ[1:Time,1:Time],
                                            c0=c0_mu[k],c1=c1_mu[k],
                                            c2=c2_mu[k],c3=c3_mu[k])
      mu_mu[k,1:Time] <- mu_0[k]*ones[1:Time]
      mu[k,1:Time] ~ dmnorm(mu_mu[k,1:Time],cov=cov_GP_mu[k,1:Time,1:Time])
      
      
      # beta0 ~ GP
      beta0_0[k] ~ dnorm(mean=0, var=1000)
      c0_beta0[k] ~ dinvgamma(shape=1e-3,scale=1e-3)
      c1_beta0[k] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
      c2_beta0[k] ~ dinvgamma(shape=1e-3,scale=1e-3)
      c3_beta0[k] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
      cov_GP_beta0[k,1:Time, 1:Time] <- expcov(dists_abs=dists_abs[1:Time,1:Time],
                                               dists_squ=dists_squ[1:Time,1:Time],
                                               c0=c0_beta0[k],c1=c1_beta0[k],
                                               c2=c2_beta0[k],c3=c3_beta0[k])
      mu_beta0[k,1:Time] <- beta0_0[k]*ones[1:Time]
      beta0[k,1:Time] ~ dmnorm(mu_beta0[k,1:Time],
                               cov=cov_GP_beta0[k,1:Time,1:Time])
      
      # beta1 ~ GP
      for(j in 1:J){
        beta1_0[k,j] ~ dnorm(mean=0, var=1000)
        c0_beta1[k,j] ~ dinvgamma(shape=1e-3,scale=1e-3)
        c1_beta1[k,j] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
        c2_beta1[k,j] ~ dinvgamma(shape=1e-3,scale=1e-3)
        c3_beta1[k,j] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
        cov_GP_beta1[k,j,1:Time, 1:Time] <- expcov(dists_abs=dists_abs[1:Time,1:Time],
                                                   dists_squ=dists_squ[1:Time,1:Time],
                                                   c0=c0_beta1[k,j],c1=c1_beta1[k,j],
                                                   c2=c2_beta1[k,j],c3=c3_beta1[k,j])
        mu_beta1[k,j,1:Time] <- beta1_0[k,j]*ones[1:Time]
        beta1[k,j,1:Time] ~ dmnorm(mu_beta1[k,j,1:Time],
                                   cov=cov_GP_beta1[k,j,1:Time,1:Time])
      }
      psi[k] ~ dunif(min=0.05,max=5)
    }
    
    # sigma2 ~ GP
    sigma2_0 ~ dnorm(mean=0, var=1000)
    c0_sigma2 ~ dinvgamma(shape=1e-3,scale=1e-3)
    c1_sigma2 ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
    c2_sigma2 ~ dinvgamma(shape=1e-3,scale=1e-3)
    c3_sigma2 ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
    cov_GP_sigma2[1:Time, 1:Time] <- expcov(dists_abs=dists_abs[1:Time,1:Time],
                                            dists_squ=dists_squ[1:Time,1:Time],
                                            c0=c0_sigma2,c1=c1_sigma2,
                                            c2=c2_sigma2,c3=c3_sigma2)
    mu_sigma2[1:Time] <- sigma2_0*ones[1:Time]
    log(sigma2[1:Time]) ~ dmnorm(mu_sigma2[1:Time],
                                 cov=cov_GP_sigma2[1:Time,1:Time])
    
    for(i in 1:N){           # LSBP
      for(k in 1:K){
        for(j in 1:J){
          kernel[i,j,k] <- exp(-sqrt(sum((s[i,1]-s_tilde[j])^2))/psi[k])
        }
      }
    }
    for(i in 1:N){
      for(k in 1:K){
        for(t in 1:Time){
          for(j in 1:J){
            g[i,j,k,t] <- beta1[k,j,t]*kernel[i,j,k]
          }
          sum_g[i,k,t] <- sum(g[i,1:J,k,t]) + beta0[k,t]
          p[i,k,t] <- 1/(1+exp(-sum_g[i,k,t]))
        }
      }
    }
    for(i in 1:N){
      for(t in 1:Time){
        w[i,1:K,t] <- stick_breaking(p[i,1:(K-1),t])
        z[i,t] ~ dcat(w[i,1:K,t])
        y[i,t] ~ dnorm(mean=mu[z[i,t],t],
                       sd=sqrt(sigma2[t])) ## likelihood
      }
    }
  })
  
  data=list(y=y)
  constants=list(s=s,s_tilde=s_tilde,dists_abs=dists_abs,dists_squ=dists_squ,
                 K=K,N=N,Time=Time,J=J,ones=rep(1,Time))
  
  # runMCMC
  
  inits=list(mu=matrix(rnorm(K*Time,mean=0,sd=1),K,Time),
             log_sigma2=rnorm(Time,mean=0,sd=1),
             beta0=matrix(rnorm(K*Time,mean=0,sd=1),K,Time),
             beta1=array(rnorm(K*J*Time,mean=0,sd=1),c(K,J,Time)),
             
             mu_0=rnorm(K,mean=0,sd=1),
             sigma2_0=rnorm(1,mean=0,sd=1),
             beta0_0=rnorm(K,mean=0,sd=1),
             beta1_0=matrix(rnorm(K*J,mean=0,sd=1),K,J),
             
             c0_mu=rinvgamma(K,shape=10,scale=10),
             c1_mu=rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             c2_mu=rinvgamma(K,shape=10,scale=10),
             c3_mu=rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             
             c0_sigma2=rinvgamma(1,shape=10,scale=10),
             c1_sigma2=rtruncnorm(1,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             c2_sigma2=rinvgamma(1,shape=10,scale=10),
             c3_sigma2=rtruncnorm(1,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             
             c0_beta0=rinvgamma(K,shape=10,scale=10),
             c1_beta0=rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             c2_beta0=rinvgamma(K,shape=10,scale=10),
             c3_beta0=rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             
             c0_beta1=matrix(rinvgamma(K*J,shape=10,scale=10),K,J),
             c1_beta1=matrix(rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),K,J),
             c2_beta1=matrix(rinvgamma(K*J,shape=10,scale=10),K,J),
             c3_beta1=matrix(rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),K,J),
             
             z=matrix(sample(1:5,N*Time,replace=T),N,Time),
             psi=runif(K,min=0.05,max=5))
  
  monitors=c('mu','sigma2','beta0','beta1',
             'mu_0','sigma2_0','beta0_0','beta1_0',
             'c0_mu','c1_mu','c2_mu','c3_mu',
             'c0_sigma2','c1_sigma2','c2_sigma2','c3_sigma2',
             'c0_beta0','c1_beta0','c2_beta0','c3_beta0',
             'c0_beta1','c1_beta1','c2_beta1','c3_beta1',
             'z','psi','p','w')
  niter=10000
  nburnin=5000
  model <- nimbleModel(code = code, data = data, constants = constants, inits = inits)
  ## Ensure we have the nodes needed to simulate new datasets
  dataNodes <- model$getNodeNames(dataOnly = TRUE)
  parentNodes <- model$getParents(dataNodes, stochOnly = TRUE)  # `getParents` is new in nimble 0.11.0
  ## Ensure we have both data nodes and deterministic intermediates (e.g., lifted nodes)
  simNodes <- model$getDependencies(parentNodes, self = FALSE)
  cmodel  <- compileNimble(model)
  mcmc    <- buildMCMC(model, monitors = parentNodes)
  cmcmc   <- compileNimble(mcmc, project = model)
  samples <- runMCMC(cmcmc, niter = niter, nburnin = nburnin, thin = 1)
  ppSamplerNF <- nimbleFunction(
    setup = function(model, mcmc) {
      dataNodes <- model$getNodeNames(dataOnly = TRUE)
      parentNodes <- model$getParents(dataNodes, stochOnly = TRUE)
      cat("Stochastic parents of data are:", paste(parentNodes, collapse = ','), ".\n")
      simNodes <- model$getDependencies(parentNodes, self = FALSE)
      vars <- mcmc$mvSamples$getVarNames()  # need ordering of variables in mvSamples / samples matrix
      cat("Using posterior samples of:", paste(vars, collapse = ','), ".\n")
      n <- length(model$expandNodeNames(dataNodes, returnScalarComponents = TRUE))
    },
    run = function(samples = double(2)) {
      nSamp <- dim(samples)[1]
      ppSamples <- matrix(nrow = nSamp, ncol = n)
      for(i in 1:nSamp) {
        values(model, vars) <<- samples[i, ]
        model$simulate(simNodes, includeData = TRUE)
        ppSamples[i, ] <- values(model, dataNodes)
      }
      returnType(double(2))
      return(ppSamples)
    })
  ## Create the sampler for this model and this MCMC.
  ppSampler <- ppSamplerNF(model, mcmc)
  cppSampler <- compileNimble(ppSampler, project = model)
  system.time(ppSamples_via_nf <- cppSampler$run(samples))
  y_hat <- ppSamples_via_nf %>% apply(.,2,median) %>% matrix(.,nrow=50)
  sigma_hat <- ppSamples_via_nf %>% apply(.,2,sd) %>% matrix(.,nrow=50)
  dim(ppSamples_via_nf)
  
  sigma2_samples <- exp(samples[,grep('sigma2',colnames(samples))])
  sigma2 <- apply(sigma2_samples,2,median)
  
  par(mfrow=c(3,3))
  for(t in 1:3){
    plot(s[,t],y[,t],ylim=c(-2,30),xlab='',ylab='')
    par(new=T)
    plot(s[,1],y_hat[,t],ylim=c(-2,30),col='red',xlab='',ylab='')
    par(new=F)
    loglik[x,t] <- mean(dnorm(x = y[,t],mean = y_hat[,t],sd=sigma2[t],log=T))
  }
  legend("topright",legend=c("data","estimated"),col=c("black","red"),
         border="white",cex=0.6,pch=c(1,3))
  for(t in 4:Time){
    plot(s[,t],y[,t],ylim=c(-2,30),xlab='',ylab='')
    par(new=T)
    plot(s[,1],y_hat[,t],ylim=c(-2,30),col='red',xlab='',ylab='')
    par(new=F)
    loglik[x,t] <- mean(dnorm(x = y[,t],mean = y_hat[,t],sd=sigma2[t],log=T))
  }
  par(mfrow=c(1,1))
  
  true_group=matrix(0,N,Time)
  for(t in 1:Time){
    true_group[,t] <-
      (s[,t] > 5+5/8*(t-1) & s[,t] < 10+5/8*(t-1))
  }
  for(t in 1:Time){
    correct_rate[x,t] <-
      sum(true_group[,t] ==
            (y_hat[,t] > 10))/N
  }
}
correct_rate %>% mean
loglik %>% mean



#### simulation nmiss = 1 #### 
library(nimble)
library(dplyr)
library(truncnorm)
# GP kernel
expcov <- nimbleFunction(
  run = function(dists_abs = double(2), dists_squ = double(2),
                 c0 = double(0), c1 = double(0),
                 c2 = double(0), c3 = double(0)) {
    returnType(double(2))
    n1 <- dim(dists_abs)[1]
    n2 <- dim(dists_abs)[2]
    result <- matrix(nrow = n1, ncol = n2, init = FALSE)
    for(i in 1:n1)
      for(j in 1:n2)
        result[i, j] <- c0*c1^dists_squ[i,j]+
      c2*c3^(sin(pi*dists_abs[i,j]/12))^2
    for(i in 1:n1)
      result[i, i] <- result[i, i] + 1e-6
    return(result)
  })
cExpcov <- compileNimble(expcov)
c=4  # number of kernel params

X=100
nmiss=1
loglik <- matrix(0,X,nmiss)
correct_rate <- matrix(0,X,nmiss)
for(x in 1:X){
  N=50
  Time = 9
  K=10
  J=5
  s=matrix(0,N,Time)
  for(t in 1:Time){
    s[,t] <- seq(0.4,20,length=N)
  }
  y=matrix(0,N,Time)
  for(t in 1:Time){
    y[,t][s[,t] > 5+5/8*(t-1) & s[,t] < 10+5/8*(t-1)] <- rnorm(sum(s[,t] > 5+5/8*(t-1) & s[,t] < 10+5/8*(t-1)),mean=20,sd=1)
    y[,t][s[,t] < 5+5/8*(t-1) | s[,t] > 10+5/8*(t-1)] <- rnorm(sum(s[,t] < 5+5/8*(t-1) | s[,t] > 10+5/8*(t-1)),mean=1,sd=1)
  }
  
  new_time = sample(2:8,size=nmiss,replace=FALSE)
  new_time = sort(new_time)
  train_time <- (1:Time)[-new_time]
  s_true <- s[,new_time]
  s_true <- s_true %>% as.matrix
  y_true <- y[,new_time]
  y_true <- y_true %>% as.matrix
  s <- s[,-new_time]
  y <- y[,-new_time]
  Time <- ncol(y)
  nmiss2=Time+nmiss
  
  # model define
  dists_abs <- as.matrix(dist(c(1:(Time+nmiss))[-new_time]))
  dists_squ <- as.matrix(dist(c(1:(Time+nmiss))[-new_time]))^2
  s_tilde <- seq(4,20,length=J)
  code <- nimbleCode({
    for(k in 1:K){
      # mu ~ GP
      mu_0[k] ~ dnorm(mean=0, var=1000)
      c0_mu[k] ~ dinvgamma(shape=1e-3,scale=1e-3)
      c1_mu[k] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
      c2_mu[k] ~ dinvgamma(shape=1e-3,scale=1e-3)
      c3_mu[k] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
      cov_GP_mu[k,1:Time, 1:Time] <- expcov(dists_abs=dists_abs[1:Time,1:Time],
                                            dists_squ=dists_squ[1:Time,1:Time],
                                            c0=c0_mu[k],c1=c1_mu[k],
                                            c2=c2_mu[k],c3=c3_mu[k])
      mu_mu[k,1:Time] <- mu_0[k]*ones[1:Time]
      mu[k,1:Time] ~ dmnorm(mu_mu[k,1:Time],cov=cov_GP_mu[k,1:Time,1:Time])
      
      
      # beta0 ~ GP
      beta0_0[k] ~ dnorm(mean=0, var=1000)
      c0_beta0[k] ~ dinvgamma(shape=1e-3,scale=1e-3)
      c1_beta0[k] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
      c2_beta0[k] ~ dinvgamma(shape=1e-3,scale=1e-3)
      c3_beta0[k] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
      cov_GP_beta0[k,1:Time, 1:Time] <- expcov(dists_abs=dists_abs[1:Time,1:Time],
                                               dists_squ=dists_squ[1:Time,1:Time],
                                               c0=c0_beta0[k],c1=c1_beta0[k],
                                               c2=c2_beta0[k],c3=c3_beta0[k])
      mu_beta0[k,1:Time] <- beta0_0[k]*ones[1:Time]
      beta0[k,1:Time] ~ dmnorm(mu_beta0[k,1:Time],
                               cov=cov_GP_beta0[k,1:Time,1:Time])
      
      # beta1 ~ GP
      for(j in 1:J){
        beta1_0[k,j] ~ dnorm(mean=0, var=1000)
        c0_beta1[k,j] ~ dinvgamma(shape=1e-3,scale=1e-3)
        c1_beta1[k,j] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
        c2_beta1[k,j] ~ dinvgamma(shape=1e-3,scale=1e-3)
        c3_beta1[k,j] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
        cov_GP_beta1[k,j,1:Time, 1:Time] <- expcov(dists_abs=dists_abs[1:Time,1:Time],
                                                   dists_squ=dists_squ[1:Time,1:Time],
                                                   c0=c0_beta1[k,j],c1=c1_beta1[k,j],
                                                   c2=c2_beta1[k,j],c3=c3_beta1[k,j])
        mu_beta1[k,j,1:Time] <- beta1_0[k,j]*ones[1:Time]
        beta1[k,j,1:Time] ~ dmnorm(mu_beta1[k,j,1:Time],
                                   cov=cov_GP_beta1[k,j,1:Time,1:Time])
      }
      psi[k] ~ dunif(min=0.05,max=5)
    }
    
    # sigma2 ~ GP
    sigma2_0 ~ dnorm(mean=0, var=1000)
    c0_sigma2 ~ dinvgamma(shape=1e-3,scale=1e-3)
    c1_sigma2 ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
    c2_sigma2 ~ dinvgamma(shape=1e-3,scale=1e-3)
    c3_sigma2 ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
    cov_GP_sigma2[1:Time, 1:Time] <- expcov(dists_abs=dists_abs[1:Time,1:Time],
                                            dists_squ=dists_squ[1:Time,1:Time],
                                            c0=c0_sigma2,c1=c1_sigma2,
                                            c2=c2_sigma2,c3=c3_sigma2)
    mu_sigma2[1:Time] <- sigma2_0*ones[1:Time]
    log(sigma2[1:Time]) ~ dmnorm(mu_sigma2[1:Time],
                                 cov=cov_GP_sigma2[1:Time,1:Time])
    
    for(i in 1:N){           # LSBP
      for(k in 1:K){
        for(j in 1:J){
          kernel[i,j,k] <- exp(-sqrt(sum((s[i,1]-s_tilde[j])^2))/psi[k])
        }
      }
    }
    for(i in 1:N){
      for(k in 1:K){
        for(t in 1:Time){
          for(j in 1:J){
            g[i,j,k,t] <- beta1[k,j,t]*kernel[i,j,k]
          }
          sum_g[i,k,t] <- sum(g[i,1:J,k,t]) + beta0[k,t]
          p[i,k,t] <- 1/(1+exp(-sum_g[i,k,t]))
        }
      }
    }
    for(i in 1:N){
      for(t in 1:Time){
        w[i,1:K,t] <- stick_breaking(p[i,1:(K-1),t])
        z[i,t] ~ dcat(w[i,1:K,t])
        y[i,t] ~ dnorm(mean=mu[z[i,t],t],
                       sd=sqrt(sigma2[t])) ## likelihood
      }
    }
  })
  
  data=list(y=y)
  constants=list(s=s,s_tilde=s_tilde,dists_abs=dists_abs,dists_squ=dists_squ,
                 K=K,N=N,Time=Time,J=J,ones=rep(1,Time))
  
  # runMCMC
  
  inits=list(mu=matrix(rnorm(K*Time,mean=0,sd=1),K,Time),
             log_sigma2=rnorm(Time,mean=0,sd=1),
             beta0=matrix(rnorm(K*Time,mean=0,sd=1),K,Time),
             beta1=array(rnorm(K*J*Time,mean=0,sd=1),c(K,J,Time)),
             
             mu_0=rnorm(K,mean=0,sd=1),
             sigma2_0=rnorm(1,mean=0,sd=1),
             beta0_0=rnorm(K,mean=0,sd=1),
             beta1_0=matrix(rnorm(K*J,mean=0,sd=1),K,J),
             
             c0_mu=rinvgamma(K,shape=10,scale=10),
             c1_mu=rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             c2_mu=rinvgamma(K,shape=10,scale=10),
             c3_mu=rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             
             c0_sigma2=rinvgamma(1,shape=10,scale=10),
             c1_sigma2=rtruncnorm(1,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             c2_sigma2=rinvgamma(1,shape=10,scale=10),
             c3_sigma2=rtruncnorm(1,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             
             c0_beta0=rinvgamma(K,shape=10,scale=10),
             c1_beta0=rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             c2_beta0=rinvgamma(K,shape=10,scale=10),
             c3_beta0=rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             
             c0_beta1=matrix(rinvgamma(K*J,shape=10,scale=10),K,J),
             c1_beta1=matrix(rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),K,J),
             c2_beta1=matrix(rinvgamma(K*J,shape=10,scale=10),K,J),
             c3_beta1=matrix(rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),K,J),
             
             z=matrix(sample(1:5,N*Time,replace=T),N,Time),
             psi=runif(K,min=0.05,max=5))
  
  monitors=c('mu','sigma2','beta0','beta1',
             'mu_0','sigma2_0','beta0_0','beta1_0',
             'c0_mu','c1_mu','c2_mu','c3_mu',
             'c0_sigma2','c1_sigma2','c2_sigma2','c3_sigma2',
             'c0_beta0','c1_beta0','c2_beta0','c3_beta0',
             'c0_beta1','c1_beta1','c2_beta1','c3_beta1',
             'z','psi','p','w')
  niter=10000
  nburnin=5000
  nmcmc <- niter-nburnin
  simulMCMC <- nimbleMCMC(code=code,constants = constants,data = data,
                          inits=inits,
                          monitors=monitors,niter=niter,nburnin=nburnin,thin=1)
  
  # seperate samples
  mu_samples <- simulMCMC[,grep('mu',colnames(simulMCMC))[(c*K+1):(c*K+K*Time)]]
  mu <- apply(mu_samples,2,median) %>% matrix(.,K,Time)
  sigma2_samples <- simulMCMC[,grep('sigma2',colnames(simulMCMC))[(c+1):(c+Time)]]
  sigma2 <- apply(sigma2_samples,2,median)
  beta0samples <- simulMCMC[,grep('beta0',colnames(simulMCMC))[1:(K*Time)]]
  beta0 <- apply(beta0samples,2,median) %>% matrix(.,K,Time)
  beta1samples <- simulMCMC[,grep('beta1',colnames(simulMCMC))[1:(K*J*Time)]]
  beta1 <- apply(beta1samples,2,median) %>% array(.,c(K,J,Time))
  
  mu_0_samples <- simulMCMC[,grep('mu_0',colnames(simulMCMC))]
  mu_0 <- apply(mu_0_samples,2,median)
  sigma2_0_samples <- simulMCMC[,grep('sigma2_0',colnames(simulMCMC))]
  sigma2_0 <- median(sigma2_0_samples)
  beta0_0_samples <- simulMCMC[,grep('beta0_0',colnames(simulMCMC))]
  beta0_0 <- apply(beta0_0_samples,2,median)
  beta1_0_samples <- simulMCMC[,grep('beta1_0',colnames(simulMCMC))]
  beta1_0 <- apply(beta1_0_samples,2,median) %>% matrix(.,K,J)
  
  c0_mu_samples <- simulMCMC[,grep('c0_mu',colnames(simulMCMC))]
  c0_mu <- apply(c0_mu_samples,2,median)
  c1_mu_samples <- simulMCMC[,grep('c1_mu',colnames(simulMCMC))]
  c1_mu <- apply(c1_mu_samples,2,median)
  c2_mu_samples <- simulMCMC[,grep('c2_mu',colnames(simulMCMC))]
  c2_mu <- apply(c2_mu_samples,2,median)
  c3_mu_samples <- simulMCMC[,grep('c3_mu',colnames(simulMCMC))]
  c3_mu <- apply(c3_mu_samples,2,median)
  
  c0_sigma2_samples <- simulMCMC[,grep('c0_sigma2',colnames(simulMCMC))]
  c0_sigma2 <- median(c0_sigma2_samples)
  c1_sigma2_samples <- simulMCMC[,grep('c1_sigma2',colnames(simulMCMC))]
  c1_sigma2 <- median(c1_sigma2_samples)
  c2_sigma2_samples <- simulMCMC[,grep('c2_sigma2',colnames(simulMCMC))]
  c2_sigma2 <- median(c2_sigma2_samples)
  c3_sigma2_samples <- simulMCMC[,grep('c3_sigma2',colnames(simulMCMC))]
  c3_sigma2 <- median(c3_sigma2_samples)
  
  c0_beta0_samples <- simulMCMC[,grep('c0_beta0',colnames(simulMCMC))]
  c0_beta0 <- apply(c0_beta0_samples,2,median)
  c1_beta0_samples <- simulMCMC[,grep('c1_beta0',colnames(simulMCMC))]
  c1_beta0 <- apply(c1_beta0_samples,2,median)
  c2_beta0_samples <- simulMCMC[,grep('c2_beta0',colnames(simulMCMC))]
  c2_beta0 <- apply(c2_beta0_samples,2,median)
  c3_beta0_samples <- simulMCMC[,grep('c3_beta0',colnames(simulMCMC))]
  c3_beta0 <- apply(c3_beta0_samples,2,median)
  
  c0_beta1_samples <- simulMCMC[,grep('c0_beta1',colnames(simulMCMC))]
  c0_beta1 <- apply(c0_beta1_samples,2,median) %>% matrix(.,K,J)
  c1_beta1_samples <- simulMCMC[,grep('c1_beta1',colnames(simulMCMC))]
  c1_beta1 <- apply(c1_beta1_samples,2,median) %>% matrix(.,K,J)
  c2_beta1_samples <- simulMCMC[,grep('c2_beta1',colnames(simulMCMC))]
  c2_beta1 <- apply(c2_beta1_samples,2,median) %>% matrix(.,K,J)
  c3_beta1_samples <- simulMCMC[,grep('c3_beta1',colnames(simulMCMC))]
  c3_beta1 <- apply(c3_beta1_samples,2,median) %>% matrix(.,K,J)
  
  z_samples <- simulMCMC[,grep('z',colnames(simulMCMC))]
  z <- apply(z_samples,2,function(x) table(x) %>% which.max %>% names) %>%
    as.integer %>% matrix(.,N,Time)
  psi_samples <- simulMCMC[,grep('psi',colnames(simulMCMC))]
  psi <- apply(psi_samples,2,median)
  
  y_pred1 <- matrix(0,N,nmiss)
  sigma2_pred1 <- c()
  for(new in 1:nmiss){
    times9_1_ind <- new_time[new]-1*(new-1)
    times1_1_ind <- new_time[new]-1*new
    times9_1 <- train_time[times9_1_ind]
    times1_1 <- train_time[times1_1_ind]
    beta1_pred_1 <- (beta1[,,times1_1_ind]*(times9_1-new_time[new]) +
                       beta1[,,times9_1_ind]*(new_time[new]-times1_1))/(times9_1-times1_1)
    
    beta0_pred_1 <- (beta0[,times1_1_ind]*(times9_1-new_time[new]) +
                       beta0[,times9_1_ind]*(new_time[new]-times1_1))/(times9_1-times1_1)
    
    mu_pred1 <- (mu[,times1_1_ind]*(times9_1-new_time[new]) +
                   mu[,times9_1_ind]*(new_time[new]-times1_1))/(times9_1-times1_1)
    
    kernel <- g <- array(0,c(N,K,J))
    sum_g <- p <- w <- matrix(0,N,K)
    for(i in 1:N){
      for(k in 1:K){
        for(j in 1:J){
          kernel[i,k,j] <- exp(-sqrt(sum((s[i,1]-s_tilde[j])^2))/psi[k])
        }
      }
    }
    for(i in 1:N){
      for(k in 1:K){
        for(j in 1:J){
          g[i,k,j] <- beta1_pred_1[k,j] * kernel[i,k,j]
        }
      }
    }
    for(i in 1:N){
      for(k in 1:K){
        sum_g[i,k] <- sum(g[i,k,1:J]) + beta0_pred_1[k]
        p[i,k] <- 1/(1+exp(-sum_g[i,k]))
      }
      w[i,] <- stick_breaking(p[i,][1:(K-1)])
    }
    # y_pred1 <- c()
    for(i in 1:N){
      y_pred1[i,new] <- sum(w[i,] * mu_pred1)
    }
    sigma2_pred1[new] <- (sigma2[times1_1_ind]*(times9_1-new_time[new]) +
                            sigma2[times9_1_ind]*(new_time[new]-times1_1))/(times9_1-times1_1)
    
  }
  par(mfrow=c(1,nmiss))
  
  for(t in 1:nmiss){
    plot(s[,1],y_true[,t],ylim=c(-1,30),xlab='',ylab='',main=new_time[t])
    par(new=T)
    plot(s[,1],y_pred1[,t],ylim=c(-1,30),xlab='',ylab='',col='red',pch=3)
    par(new=F)
    loglik[x,t] <- mean(dnorm(x = y_true[,t],mean = y_pred1[,t],sd=sqrt(sigma2_pred1[t]),log=T))
  }
  par(mfrow=c(1,1))
  true_group=matrix(0,N,nmiss)
  for(t in 1:nmiss){
    true_group[,t] <-
      (s_true[,t] > 5+5/8*(new_time[t]-1) & s_true[,t] < 10+5/8*(new_time[t]-1))
    
  }
  for(t in 1:nmiss){
    correct_rate[x,t] <-
      sum(true_group[,t] ==
            (y_pred1[,t] > 10))/N
    
  }
}
loglik %>% mean
correct_rate %>% mean

#### simulation nmiss = 2 #### 
library(nimble)
library(dplyr)
library(truncnorm)
# GP kernel
expcov <- nimbleFunction(
  run = function(dists_abs = double(2), dists_squ = double(2),
                 c0 = double(0), c1 = double(0),
                 c2 = double(0), c3 = double(0)) {
    returnType(double(2))
    n1 <- dim(dists_abs)[1]
    n2 <- dim(dists_abs)[2]
    result <- matrix(nrow = n1, ncol = n2, init = FALSE)
    for(i in 1:n1)
      for(j in 1:n2)
        result[i, j] <- c0*c1^dists_squ[i,j]+
      c2*c3^(sin(pi*dists_abs[i,j]/12))^2
    for(i in 1:n1)
      result[i, i] <- result[i, i] + 1e-6
    return(result)
  })
cExpcov <- compileNimble(expcov)
c=4  # number of kernel params

X=100
nmiss=2
loglik <- matrix(0,X,nmiss)
correct_rate <- matrix(0,X,nmiss)
for(x in 1:X){
  N=50
  Time = 9
  K=10
  J=5
  s=matrix(0,N,Time)
  for(t in 1:Time){
    s[,t] <- seq(0.4,20,length=N)
  }
  y=matrix(0,N,Time)
  for(t in 1:Time){
    y[,t][s[,t] > 5+5/8*(t-1) & s[,t] < 10+5/8*(t-1)] <- rnorm(sum(s[,t] > 5+5/8*(t-1) & s[,t] < 10+5/8*(t-1)),mean=20,sd=1)
    y[,t][s[,t] < 5+5/8*(t-1) | s[,t] > 10+5/8*(t-1)] <- rnorm(sum(s[,t] < 5+5/8*(t-1) | s[,t] > 10+5/8*(t-1)),mean=1,sd=1)
  }
  
  new_time = sample(2:8,size=nmiss,replace=FALSE)
  new_time = sort(new_time)
  train_time <- (1:Time)[-new_time]
  s_true <- s[,new_time]
  s_true <- s_true %>% as.matrix
  y_true <- y[,new_time]
  y_true <- y_true %>% as.matrix
  s <- s[,-new_time]
  y <- y[,-new_time]
  Time <- ncol(y)
  nmiss2=Time+nmiss
  
  # model define
  dists_abs <- as.matrix(dist(c(1:(Time+nmiss))[-new_time]))
  dists_squ <- as.matrix(dist(c(1:(Time+nmiss))[-new_time]))^2
  s_tilde <- seq(4,20,length=J)
  code <- nimbleCode({
    for(k in 1:K){
      # mu ~ GP
      mu_0[k] ~ dnorm(mean=0, var=1000)
      c0_mu[k] ~ dinvgamma(shape=1e-3,scale=1e-3)
      c1_mu[k] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
      c2_mu[k] ~ dinvgamma(shape=1e-3,scale=1e-3)
      c3_mu[k] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
      cov_GP_mu[k,1:Time, 1:Time] <- expcov(dists_abs=dists_abs[1:Time,1:Time],
                                            dists_squ=dists_squ[1:Time,1:Time],
                                            c0=c0_mu[k],c1=c1_mu[k],
                                            c2=c2_mu[k],c3=c3_mu[k])
      mu_mu[k,1:Time] <- mu_0[k]*ones[1:Time]
      mu[k,1:Time] ~ dmnorm(mu_mu[k,1:Time],cov=cov_GP_mu[k,1:Time,1:Time])
      
      
      # beta0 ~ GP
      beta0_0[k] ~ dnorm(mean=0, var=1000)
      c0_beta0[k] ~ dinvgamma(shape=1e-3,scale=1e-3)
      c1_beta0[k] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
      c2_beta0[k] ~ dinvgamma(shape=1e-3,scale=1e-3)
      c3_beta0[k] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
      cov_GP_beta0[k,1:Time, 1:Time] <- expcov(dists_abs=dists_abs[1:Time,1:Time],
                                               dists_squ=dists_squ[1:Time,1:Time],
                                               c0=c0_beta0[k],c1=c1_beta0[k],
                                               c2=c2_beta0[k],c3=c3_beta0[k])
      mu_beta0[k,1:Time] <- beta0_0[k]*ones[1:Time]
      beta0[k,1:Time] ~ dmnorm(mu_beta0[k,1:Time],
                               cov=cov_GP_beta0[k,1:Time,1:Time])
      
      # beta1 ~ GP
      for(j in 1:J){
        beta1_0[k,j] ~ dnorm(mean=0, var=1000)
        c0_beta1[k,j] ~ dinvgamma(shape=1e-3,scale=1e-3)
        c1_beta1[k,j] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
        c2_beta1[k,j] ~ dinvgamma(shape=1e-3,scale=1e-3)
        c3_beta1[k,j] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
        cov_GP_beta1[k,j,1:Time, 1:Time] <- expcov(dists_abs=dists_abs[1:Time,1:Time],
                                                   dists_squ=dists_squ[1:Time,1:Time],
                                                   c0=c0_beta1[k,j],c1=c1_beta1[k,j],
                                                   c2=c2_beta1[k,j],c3=c3_beta1[k,j])
        mu_beta1[k,j,1:Time] <- beta1_0[k,j]*ones[1:Time]
        beta1[k,j,1:Time] ~ dmnorm(mu_beta1[k,j,1:Time],
                                   cov=cov_GP_beta1[k,j,1:Time,1:Time])
      }
      psi[k] ~ dunif(min=0.05,max=5)
    }
    
    # sigma2 ~ GP
    sigma2_0 ~ dnorm(mean=0, var=1000)
    c0_sigma2 ~ dinvgamma(shape=1e-3,scale=1e-3)
    c1_sigma2 ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
    c2_sigma2 ~ dinvgamma(shape=1e-3,scale=1e-3)
    c3_sigma2 ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
    cov_GP_sigma2[1:Time, 1:Time] <- expcov(dists_abs=dists_abs[1:Time,1:Time],
                                            dists_squ=dists_squ[1:Time,1:Time],
                                            c0=c0_sigma2,c1=c1_sigma2,
                                            c2=c2_sigma2,c3=c3_sigma2)
    mu_sigma2[1:Time] <- sigma2_0*ones[1:Time]
    log(sigma2[1:Time]) ~ dmnorm(mu_sigma2[1:Time],
                                 cov=cov_GP_sigma2[1:Time,1:Time])
    
    for(i in 1:N){           # LSBP
      for(k in 1:K){
        for(j in 1:J){
          kernel[i,j,k] <- exp(-sqrt(sum((s[i,1]-s_tilde[j])^2))/psi[k])
        }
      }
    }
    for(i in 1:N){
      for(k in 1:K){
        for(t in 1:Time){
          for(j in 1:J){
            g[i,j,k,t] <- beta1[k,j,t]*kernel[i,j,k]
          }
          sum_g[i,k,t] <- sum(g[i,1:J,k,t]) + beta0[k,t]
          p[i,k,t] <- 1/(1+exp(-sum_g[i,k,t]))
        }
      }
    }
    for(i in 1:N){
      for(t in 1:Time){
        w[i,1:K,t] <- stick_breaking(p[i,1:(K-1),t])
        z[i,t] ~ dcat(w[i,1:K,t])
        y[i,t] ~ dnorm(mean=mu[z[i,t],t],
                       sd=sqrt(sigma2[t])) ## likelihood
      }
    }
  })
  
  data=list(y=y)
  constants=list(s=s,s_tilde=s_tilde,dists_abs=dists_abs,dists_squ=dists_squ,
                 K=K,N=N,Time=Time,J=J,ones=rep(1,Time))
  
  # runMCMC
  
  inits=list(mu=matrix(rnorm(K*Time,mean=0,sd=1),K,Time),
             log_sigma2=rnorm(Time,mean=0,sd=1),
             beta0=matrix(rnorm(K*Time,mean=0,sd=1),K,Time),
             beta1=array(rnorm(K*J*Time,mean=0,sd=1),c(K,J,Time)),
             
             mu_0=rnorm(K,mean=0,sd=1),
             sigma2_0=rnorm(1,mean=0,sd=1),
             beta0_0=rnorm(K,mean=0,sd=1),
             beta1_0=matrix(rnorm(K*J,mean=0,sd=1),K,J),
             
             c0_mu=rinvgamma(K,shape=10,scale=10),
             c1_mu=rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             c2_mu=rinvgamma(K,shape=10,scale=10),
             c3_mu=rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             
             c0_sigma2=rinvgamma(1,shape=10,scale=10),
             c1_sigma2=rtruncnorm(1,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             c2_sigma2=rinvgamma(1,shape=10,scale=10),
             c3_sigma2=rtruncnorm(1,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             
             c0_beta0=rinvgamma(K,shape=10,scale=10),
             c1_beta0=rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             c2_beta0=rinvgamma(K,shape=10,scale=10),
             c3_beta0=rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             
             c0_beta1=matrix(rinvgamma(K*J,shape=10,scale=10),K,J),
             c1_beta1=matrix(rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),K,J),
             c2_beta1=matrix(rinvgamma(K*J,shape=10,scale=10),K,J),
             c3_beta1=matrix(rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),K,J),
             
             z=matrix(sample(1:5,N*Time,replace=T),N,Time),
             psi=runif(K,min=0.05,max=5))
  
  monitors=c('mu','sigma2','beta0','beta1',
             'mu_0','sigma2_0','beta0_0','beta1_0',
             'c0_mu','c1_mu','c2_mu','c3_mu',
             'c0_sigma2','c1_sigma2','c2_sigma2','c3_sigma2',
             'c0_beta0','c1_beta0','c2_beta0','c3_beta0',
             'c0_beta1','c1_beta1','c2_beta1','c3_beta1',
             'z','psi','p','w')
  niter=10000
  nburnin=5000
  nmcmc <- niter-nburnin
  simulMCMC <- nimbleMCMC(code=code,constants = constants,data = data,
                          inits=inits,
                          monitors=monitors,niter=niter,nburnin=nburnin,thin=1)
  
  # seperate samples
  mu_samples <- simulMCMC[,grep('mu',colnames(simulMCMC))[(c*K+1):(c*K+K*Time)]]
  mu <- apply(mu_samples,2,median) %>% matrix(.,K,Time)
  sigma2_samples <- simulMCMC[,grep('sigma2',colnames(simulMCMC))[(c+1):(c+Time)]]
  sigma2 <- apply(sigma2_samples,2,median)
  beta0samples <- simulMCMC[,grep('beta0',colnames(simulMCMC))[1:(K*Time)]]
  beta0 <- apply(beta0samples,2,median) %>% matrix(.,K,Time)
  beta1samples <- simulMCMC[,grep('beta1',colnames(simulMCMC))[1:(K*J*Time)]]
  beta1 <- apply(beta1samples,2,median) %>% array(.,c(K,J,Time))
  
  mu_0_samples <- simulMCMC[,grep('mu_0',colnames(simulMCMC))]
  mu_0 <- apply(mu_0_samples,2,median)
  sigma2_0_samples <- simulMCMC[,grep('sigma2_0',colnames(simulMCMC))]
  sigma2_0 <- median(sigma2_0_samples)
  beta0_0_samples <- simulMCMC[,grep('beta0_0',colnames(simulMCMC))]
  beta0_0 <- apply(beta0_0_samples,2,median)
  beta1_0_samples <- simulMCMC[,grep('beta1_0',colnames(simulMCMC))]
  beta1_0 <- apply(beta1_0_samples,2,median) %>% matrix(.,K,J)
  
  c0_mu_samples <- simulMCMC[,grep('c0_mu',colnames(simulMCMC))]
  c0_mu <- apply(c0_mu_samples,2,median)
  c1_mu_samples <- simulMCMC[,grep('c1_mu',colnames(simulMCMC))]
  c1_mu <- apply(c1_mu_samples,2,median)
  c2_mu_samples <- simulMCMC[,grep('c2_mu',colnames(simulMCMC))]
  c2_mu <- apply(c2_mu_samples,2,median)
  c3_mu_samples <- simulMCMC[,grep('c3_mu',colnames(simulMCMC))]
  c3_mu <- apply(c3_mu_samples,2,median)
  
  c0_sigma2_samples <- simulMCMC[,grep('c0_sigma2',colnames(simulMCMC))]
  c0_sigma2 <- median(c0_sigma2_samples)
  c1_sigma2_samples <- simulMCMC[,grep('c1_sigma2',colnames(simulMCMC))]
  c1_sigma2 <- median(c1_sigma2_samples)
  c2_sigma2_samples <- simulMCMC[,grep('c2_sigma2',colnames(simulMCMC))]
  c2_sigma2 <- median(c2_sigma2_samples)
  c3_sigma2_samples <- simulMCMC[,grep('c3_sigma2',colnames(simulMCMC))]
  c3_sigma2 <- median(c3_sigma2_samples)
  
  c0_beta0_samples <- simulMCMC[,grep('c0_beta0',colnames(simulMCMC))]
  c0_beta0 <- apply(c0_beta0_samples,2,median)
  c1_beta0_samples <- simulMCMC[,grep('c1_beta0',colnames(simulMCMC))]
  c1_beta0 <- apply(c1_beta0_samples,2,median)
  c2_beta0_samples <- simulMCMC[,grep('c2_beta0',colnames(simulMCMC))]
  c2_beta0 <- apply(c2_beta0_samples,2,median)
  c3_beta0_samples <- simulMCMC[,grep('c3_beta0',colnames(simulMCMC))]
  c3_beta0 <- apply(c3_beta0_samples,2,median)
  
  c0_beta1_samples <- simulMCMC[,grep('c0_beta1',colnames(simulMCMC))]
  c0_beta1 <- apply(c0_beta1_samples,2,median) %>% matrix(.,K,J)
  c1_beta1_samples <- simulMCMC[,grep('c1_beta1',colnames(simulMCMC))]
  c1_beta1 <- apply(c1_beta1_samples,2,median) %>% matrix(.,K,J)
  c2_beta1_samples <- simulMCMC[,grep('c2_beta1',colnames(simulMCMC))]
  c2_beta1 <- apply(c2_beta1_samples,2,median) %>% matrix(.,K,J)
  c3_beta1_samples <- simulMCMC[,grep('c3_beta1',colnames(simulMCMC))]
  c3_beta1 <- apply(c3_beta1_samples,2,median) %>% matrix(.,K,J)
  
  z_samples <- simulMCMC[,grep('z',colnames(simulMCMC))]
  z <- apply(z_samples,2,function(x) table(x) %>% which.max %>% names) %>%
    as.integer %>% matrix(.,N,Time)
  psi_samples <- simulMCMC[,grep('psi',colnames(simulMCMC))]
  psi <- apply(psi_samples,2,median)
  
  y_pred1 <- matrix(0,N,nmiss)
  sigma2_pred1 <- c()
  for(new in 1:nmiss){
    times9_1_ind <- new_time[new]-1*(new-1)
    times1_1_ind <- new_time[new]-1*new
    times9_1 <- train_time[times9_1_ind]
    times1_1 <- train_time[times1_1_ind]
    beta1_pred_1 <- (beta1[,,times1_1_ind]*(times9_1-new_time[new]) +
                       beta1[,,times9_1_ind]*(new_time[new]-times1_1))/(times9_1-times1_1)
    
    beta0_pred_1 <- (beta0[,times1_1_ind]*(times9_1-new_time[new]) +
                       beta0[,times9_1_ind]*(new_time[new]-times1_1))/(times9_1-times1_1)
    
    mu_pred1 <- (mu[,times1_1_ind]*(times9_1-new_time[new]) +
                   mu[,times9_1_ind]*(new_time[new]-times1_1))/(times9_1-times1_1)
    
    kernel <- g <- array(0,c(N,K,J))
    sum_g <- p <- w <- matrix(0,N,K)
    for(i in 1:N){
      for(k in 1:K){
        for(j in 1:J){
          kernel[i,k,j] <- exp(-sqrt(sum((s[i,1]-s_tilde[j])^2))/psi[k])
        }
      }
    }
    for(i in 1:N){
      for(k in 1:K){
        for(j in 1:J){
          g[i,k,j] <- beta1_pred_1[k,j] * kernel[i,k,j]
        }
      }
    }
    for(i in 1:N){
      for(k in 1:K){
        sum_g[i,k] <- sum(g[i,k,1:J]) + beta0_pred_1[k]
        p[i,k] <- 1/(1+exp(-sum_g[i,k]))
      }
      w[i,] <- stick_breaking(p[i,][1:(K-1)])
    }
    # y_pred1 <- c()
    for(i in 1:N){
      y_pred1[i,new] <- sum(w[i,] * mu_pred1)
    }
    sigma2_pred1[new] <- (sigma2[times1_1_ind]*(times9_1-new_time[new]) +
                            sigma2[times9_1_ind]*(new_time[new]-times1_1))/(times9_1-times1_1)
    
  }
  par(mfrow=c(1,nmiss))
  
  for(t in 1:nmiss){
    plot(s[,1],y_true[,t],ylim=c(-1,30),xlab='',ylab='',main=new_time[t])
    par(new=T)
    plot(s[,1],y_pred1[,t],ylim=c(-1,30),xlab='',ylab='',col='red',pch=3)
    par(new=F)
    loglik[x,t] <- mean(dnorm(x = y_true[,t],mean = y_pred1[,t],sd=sqrt(sigma2_pred1[t]),log=T))
  }
  par(mfrow=c(1,1))
  true_group=matrix(0,N,nmiss)
  for(t in 1:nmiss){
    true_group[,t] <-
      (s_true[,t] > 5+5/8*(new_time[t]-1) & s_true[,t] < 10+5/8*(new_time[t]-1))
    
  }
  for(t in 1:nmiss){
    correct_rate[x,t] <-
      sum(true_group[,t] ==
            (y_pred1[,t] > 10))/N
    
  }
}
loglik %>% mean
correct_rate %>% mean
#### simulation nmiss = 3 #### 
library(nimble)
library(dplyr)
library(truncnorm)
# GP kernel
expcov <- nimbleFunction(
  run = function(dists_abs = double(2), dists_squ = double(2),
                 c0 = double(0), c1 = double(0),
                 c2 = double(0), c3 = double(0)) {
    returnType(double(2))
    n1 <- dim(dists_abs)[1]
    n2 <- dim(dists_abs)[2]
    result <- matrix(nrow = n1, ncol = n2, init = FALSE)
    for(i in 1:n1)
      for(j in 1:n2)
        result[i, j] <- c0*c1^dists_squ[i,j]+
      c2*c3^(sin(pi*dists_abs[i,j]/12))^2
    for(i in 1:n1)
      result[i, i] <- result[i, i] + 1e-6
    return(result)
  })
cExpcov <- compileNimble(expcov)
c=4  # number of kernel params

X=100
nmiss=3
loglik <- matrix(0,X,nmiss)
correct_rate <- matrix(0,X,nmiss)
for(x in 1:X){
  N=50
  Time = 9
  K=10
  J=5
  s=matrix(0,N,Time)
  for(t in 1:Time){
    s[,t] <- seq(0.4,20,length=N)
  }
  y=matrix(0,N,Time)
  for(t in 1:Time){
    y[,t][s[,t] > 5+5/8*(t-1) & s[,t] < 10+5/8*(t-1)] <- rnorm(sum(s[,t] > 5+5/8*(t-1) & s[,t] < 10+5/8*(t-1)),mean=20,sd=1)
    y[,t][s[,t] < 5+5/8*(t-1) | s[,t] > 10+5/8*(t-1)] <- rnorm(sum(s[,t] < 5+5/8*(t-1) | s[,t] > 10+5/8*(t-1)),mean=1,sd=1)
  }
  
  new_time = sample(2:8,size=nmiss,replace=FALSE)
  new_time = sort(new_time)
  train_time <- (1:Time)[-new_time]
  s_true <- s[,new_time]
  s_true <- s_true %>% as.matrix
  y_true <- y[,new_time]
  y_true <- y_true %>% as.matrix
  s <- s[,-new_time]
  y <- y[,-new_time]
  Time <- ncol(y)
  nmiss2=Time+nmiss
  
  # model define
  dists_abs <- as.matrix(dist(c(1:(Time+nmiss))[-new_time]))
  dists_squ <- as.matrix(dist(c(1:(Time+nmiss))[-new_time]))^2
  s_tilde <- seq(4,20,length=J)
  code <- nimbleCode({
    for(k in 1:K){
      # mu ~ GP
      mu_0[k] ~ dnorm(mean=0, var=1000)
      c0_mu[k] ~ dinvgamma(shape=1e-3,scale=1e-3)
      c1_mu[k] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
      c2_mu[k] ~ dinvgamma(shape=1e-3,scale=1e-3)
      c3_mu[k] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
      cov_GP_mu[k,1:Time, 1:Time] <- expcov(dists_abs=dists_abs[1:Time,1:Time],
                                            dists_squ=dists_squ[1:Time,1:Time],
                                            c0=c0_mu[k],c1=c1_mu[k],
                                            c2=c2_mu[k],c3=c3_mu[k])
      mu_mu[k,1:Time] <- mu_0[k]*ones[1:Time]
      mu[k,1:Time] ~ dmnorm(mu_mu[k,1:Time],cov=cov_GP_mu[k,1:Time,1:Time])
      
      
      # beta0 ~ GP
      beta0_0[k] ~ dnorm(mean=0, var=1000)
      c0_beta0[k] ~ dinvgamma(shape=1e-3,scale=1e-3)
      c1_beta0[k] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
      c2_beta0[k] ~ dinvgamma(shape=1e-3,scale=1e-3)
      c3_beta0[k] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
      cov_GP_beta0[k,1:Time, 1:Time] <- expcov(dists_abs=dists_abs[1:Time,1:Time],
                                               dists_squ=dists_squ[1:Time,1:Time],
                                               c0=c0_beta0[k],c1=c1_beta0[k],
                                               c2=c2_beta0[k],c3=c3_beta0[k])
      mu_beta0[k,1:Time] <- beta0_0[k]*ones[1:Time]
      beta0[k,1:Time] ~ dmnorm(mu_beta0[k,1:Time],
                               cov=cov_GP_beta0[k,1:Time,1:Time])
      
      # beta1 ~ GP
      for(j in 1:J){
        beta1_0[k,j] ~ dnorm(mean=0, var=1000)
        c0_beta1[k,j] ~ dinvgamma(shape=1e-3,scale=1e-3)
        c1_beta1[k,j] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
        c2_beta1[k,j] ~ dinvgamma(shape=1e-3,scale=1e-3)
        c3_beta1[k,j] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
        cov_GP_beta1[k,j,1:Time, 1:Time] <- expcov(dists_abs=dists_abs[1:Time,1:Time],
                                                   dists_squ=dists_squ[1:Time,1:Time],
                                                   c0=c0_beta1[k,j],c1=c1_beta1[k,j],
                                                   c2=c2_beta1[k,j],c3=c3_beta1[k,j])
        mu_beta1[k,j,1:Time] <- beta1_0[k,j]*ones[1:Time]
        beta1[k,j,1:Time] ~ dmnorm(mu_beta1[k,j,1:Time],
                                   cov=cov_GP_beta1[k,j,1:Time,1:Time])
      }
      psi[k] ~ dunif(min=0.05,max=5)
    }
    
    # sigma2 ~ GP
    sigma2_0 ~ dnorm(mean=0, var=1000)
    c0_sigma2 ~ dinvgamma(shape=1e-3,scale=1e-3)
    c1_sigma2 ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
    c2_sigma2 ~ dinvgamma(shape=1e-3,scale=1e-3)
    c3_sigma2 ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
    cov_GP_sigma2[1:Time, 1:Time] <- expcov(dists_abs=dists_abs[1:Time,1:Time],
                                            dists_squ=dists_squ[1:Time,1:Time],
                                            c0=c0_sigma2,c1=c1_sigma2,
                                            c2=c2_sigma2,c3=c3_sigma2)
    mu_sigma2[1:Time] <- sigma2_0*ones[1:Time]
    log(sigma2[1:Time]) ~ dmnorm(mu_sigma2[1:Time],
                                 cov=cov_GP_sigma2[1:Time,1:Time])
    
    for(i in 1:N){           # LSBP
      for(k in 1:K){
        for(j in 1:J){
          kernel[i,j,k] <- exp(-sqrt(sum((s[i,1]-s_tilde[j])^2))/psi[k])
        }
      }
    }
    for(i in 1:N){
      for(k in 1:K){
        for(t in 1:Time){
          for(j in 1:J){
            g[i,j,k,t] <- beta1[k,j,t]*kernel[i,j,k]
          }
          sum_g[i,k,t] <- sum(g[i,1:J,k,t]) + beta0[k,t]
          p[i,k,t] <- 1/(1+exp(-sum_g[i,k,t]))
        }
      }
    }
    for(i in 1:N){
      for(t in 1:Time){
        w[i,1:K,t] <- stick_breaking(p[i,1:(K-1),t])
        z[i,t] ~ dcat(w[i,1:K,t])
        y[i,t] ~ dnorm(mean=mu[z[i,t],t],
                       sd=sqrt(sigma2[t])) ## likelihood
      }
    }
  })
  
  data=list(y=y)
  constants=list(s=s,s_tilde=s_tilde,dists_abs=dists_abs,dists_squ=dists_squ,
                 K=K,N=N,Time=Time,J=J,ones=rep(1,Time))
  
  # runMCMC
  
  inits=list(mu=matrix(rnorm(K*Time,mean=0,sd=1),K,Time),
             log_sigma2=rnorm(Time,mean=0,sd=1),
             beta0=matrix(rnorm(K*Time,mean=0,sd=1),K,Time),
             beta1=array(rnorm(K*J*Time,mean=0,sd=1),c(K,J,Time)),
             
             mu_0=rnorm(K,mean=0,sd=1),
             sigma2_0=rnorm(1,mean=0,sd=1),
             beta0_0=rnorm(K,mean=0,sd=1),
             beta1_0=matrix(rnorm(K*J,mean=0,sd=1),K,J),
             
             c0_mu=rinvgamma(K,shape=10,scale=10),
             c1_mu=rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             c2_mu=rinvgamma(K,shape=10,scale=10),
             c3_mu=rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             
             c0_sigma2=rinvgamma(1,shape=10,scale=10),
             c1_sigma2=rtruncnorm(1,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             c2_sigma2=rinvgamma(1,shape=10,scale=10),
             c3_sigma2=rtruncnorm(1,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             
             c0_beta0=rinvgamma(K,shape=10,scale=10),
             c1_beta0=rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             c2_beta0=rinvgamma(K,shape=10,scale=10),
             c3_beta0=rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             
             c0_beta1=matrix(rinvgamma(K*J,shape=10,scale=10),K,J),
             c1_beta1=matrix(rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),K,J),
             c2_beta1=matrix(rinvgamma(K*J,shape=10,scale=10),K,J),
             c3_beta1=matrix(rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),K,J),
             
             z=matrix(sample(1:5,N*Time,replace=T),N,Time),
             psi=runif(K,min=0.05,max=5))
  
  monitors=c('mu','sigma2','beta0','beta1',
             'mu_0','sigma2_0','beta0_0','beta1_0',
             'c0_mu','c1_mu','c2_mu','c3_mu',
             'c0_sigma2','c1_sigma2','c2_sigma2','c3_sigma2',
             'c0_beta0','c1_beta0','c2_beta0','c3_beta0',
             'c0_beta1','c1_beta1','c2_beta1','c3_beta1',
             'z','psi','p','w')
  niter=10000
  nburnin=5000
  nmcmc <- niter-nburnin
  simulMCMC <- nimbleMCMC(code=code,constants = constants,data = data,
                          inits=inits,
                          monitors=monitors,niter=niter,nburnin=nburnin,thin=1)
  
  # seperate samples
  mu_samples <- simulMCMC[,grep('mu',colnames(simulMCMC))[(c*K+1):(c*K+K*Time)]]
  mu <- apply(mu_samples,2,median) %>% matrix(.,K,Time)
  sigma2_samples <- simulMCMC[,grep('sigma2',colnames(simulMCMC))[(c+1):(c+Time)]]
  sigma2 <- apply(sigma2_samples,2,median)
  beta0samples <- simulMCMC[,grep('beta0',colnames(simulMCMC))[1:(K*Time)]]
  beta0 <- apply(beta0samples,2,median) %>% matrix(.,K,Time)
  beta1samples <- simulMCMC[,grep('beta1',colnames(simulMCMC))[1:(K*J*Time)]]
  beta1 <- apply(beta1samples,2,median) %>% array(.,c(K,J,Time))
  
  mu_0_samples <- simulMCMC[,grep('mu_0',colnames(simulMCMC))]
  mu_0 <- apply(mu_0_samples,2,median)
  sigma2_0_samples <- simulMCMC[,grep('sigma2_0',colnames(simulMCMC))]
  sigma2_0 <- median(sigma2_0_samples)
  beta0_0_samples <- simulMCMC[,grep('beta0_0',colnames(simulMCMC))]
  beta0_0 <- apply(beta0_0_samples,2,median)
  beta1_0_samples <- simulMCMC[,grep('beta1_0',colnames(simulMCMC))]
  beta1_0 <- apply(beta1_0_samples,2,median) %>% matrix(.,K,J)
  
  c0_mu_samples <- simulMCMC[,grep('c0_mu',colnames(simulMCMC))]
  c0_mu <- apply(c0_mu_samples,2,median)
  c1_mu_samples <- simulMCMC[,grep('c1_mu',colnames(simulMCMC))]
  c1_mu <- apply(c1_mu_samples,2,median)
  c2_mu_samples <- simulMCMC[,grep('c2_mu',colnames(simulMCMC))]
  c2_mu <- apply(c2_mu_samples,2,median)
  c3_mu_samples <- simulMCMC[,grep('c3_mu',colnames(simulMCMC))]
  c3_mu <- apply(c3_mu_samples,2,median)
  
  c0_sigma2_samples <- simulMCMC[,grep('c0_sigma2',colnames(simulMCMC))]
  c0_sigma2 <- median(c0_sigma2_samples)
  c1_sigma2_samples <- simulMCMC[,grep('c1_sigma2',colnames(simulMCMC))]
  c1_sigma2 <- median(c1_sigma2_samples)
  c2_sigma2_samples <- simulMCMC[,grep('c2_sigma2',colnames(simulMCMC))]
  c2_sigma2 <- median(c2_sigma2_samples)
  c3_sigma2_samples <- simulMCMC[,grep('c3_sigma2',colnames(simulMCMC))]
  c3_sigma2 <- median(c3_sigma2_samples)
  
  c0_beta0_samples <- simulMCMC[,grep('c0_beta0',colnames(simulMCMC))]
  c0_beta0 <- apply(c0_beta0_samples,2,median)
  c1_beta0_samples <- simulMCMC[,grep('c1_beta0',colnames(simulMCMC))]
  c1_beta0 <- apply(c1_beta0_samples,2,median)
  c2_beta0_samples <- simulMCMC[,grep('c2_beta0',colnames(simulMCMC))]
  c2_beta0 <- apply(c2_beta0_samples,2,median)
  c3_beta0_samples <- simulMCMC[,grep('c3_beta0',colnames(simulMCMC))]
  c3_beta0 <- apply(c3_beta0_samples,2,median)
  
  c0_beta1_samples <- simulMCMC[,grep('c0_beta1',colnames(simulMCMC))]
  c0_beta1 <- apply(c0_beta1_samples,2,median) %>% matrix(.,K,J)
  c1_beta1_samples <- simulMCMC[,grep('c1_beta1',colnames(simulMCMC))]
  c1_beta1 <- apply(c1_beta1_samples,2,median) %>% matrix(.,K,J)
  c2_beta1_samples <- simulMCMC[,grep('c2_beta1',colnames(simulMCMC))]
  c2_beta1 <- apply(c2_beta1_samples,2,median) %>% matrix(.,K,J)
  c3_beta1_samples <- simulMCMC[,grep('c3_beta1',colnames(simulMCMC))]
  c3_beta1 <- apply(c3_beta1_samples,2,median) %>% matrix(.,K,J)
  
  z_samples <- simulMCMC[,grep('z',colnames(simulMCMC))]
  z <- apply(z_samples,2,function(x) table(x) %>% which.max %>% names) %>%
    as.integer %>% matrix(.,N,Time)
  psi_samples <- simulMCMC[,grep('psi',colnames(simulMCMC))]
  psi <- apply(psi_samples,2,median)
  
  y_pred1 <- matrix(0,N,nmiss)
  sigma2_pred1 <- c()
  for(new in 1:nmiss){
    times9_1_ind <- new_time[new]-1*(new-1)
    times1_1_ind <- new_time[new]-1*new
    times9_1 <- train_time[times9_1_ind]
    times1_1 <- train_time[times1_1_ind]
    beta1_pred_1 <- (beta1[,,times1_1_ind]*(times9_1-new_time[new]) +
                       beta1[,,times9_1_ind]*(new_time[new]-times1_1))/(times9_1-times1_1)
    
    beta0_pred_1 <- (beta0[,times1_1_ind]*(times9_1-new_time[new]) +
                       beta0[,times9_1_ind]*(new_time[new]-times1_1))/(times9_1-times1_1)
    
    mu_pred1 <- (mu[,times1_1_ind]*(times9_1-new_time[new]) +
                   mu[,times9_1_ind]*(new_time[new]-times1_1))/(times9_1-times1_1)
    
    kernel <- g <- array(0,c(N,K,J))
    sum_g <- p <- w <- matrix(0,N,K)
    for(i in 1:N){
      for(k in 1:K){
        for(j in 1:J){
          kernel[i,k,j] <- exp(-sqrt(sum((s[i,1]-s_tilde[j])^2))/psi[k])
        }
      }
    }
    for(i in 1:N){
      for(k in 1:K){
        for(j in 1:J){
          g[i,k,j] <- beta1_pred_1[k,j] * kernel[i,k,j]
        }
      }
    }
    for(i in 1:N){
      for(k in 1:K){
        sum_g[i,k] <- sum(g[i,k,1:J]) + beta0_pred_1[k]
        p[i,k] <- 1/(1+exp(-sum_g[i,k]))
      }
      w[i,] <- stick_breaking(p[i,][1:(K-1)])
    }
    # y_pred1 <- c()
    for(i in 1:N){
      y_pred1[i,new] <- sum(w[i,] * mu_pred1)
    }
    sigma2_pred1[new] <- (sigma2[times1_1_ind]*(times9_1-new_time[new]) +
                            sigma2[times9_1_ind]*(new_time[new]-times1_1))/(times9_1-times1_1)
    
  }
  par(mfrow=c(1,nmiss))
  
  for(t in 1:nmiss){
    plot(s[,1],y_true[,t],ylim=c(-1,30),xlab='',ylab='',main=new_time[t])
    par(new=T)
    plot(s[,1],y_pred1[,t],ylim=c(-1,30),xlab='',ylab='',col='red',pch=3)
    par(new=F)
    loglik[x,t] <- mean(dnorm(x = y_true[,t],mean = y_pred1[,t],sd=sqrt(sigma2_pred1[t]),log=T))
  }
  par(mfrow=c(1,1))
  true_group=matrix(0,N,nmiss)
  for(t in 1:nmiss){
    true_group[,t] <-
      (s_true[,t] > 5+5/8*(new_time[t]-1) & s_true[,t] < 10+5/8*(new_time[t]-1))
    
  }
  for(t in 1:nmiss){
    correct_rate[x,t] <-
      sum(true_group[,t] ==
            (y_pred1[,t] > 10))/N
    
  }
}
loglik %>% mean
correct_rate %>% mean
#### simulation nmiss = 4 #### 
library(nimble)
library(dplyr)
library(truncnorm)
# GP kernel
expcov <- nimbleFunction(
  run = function(dists_abs = double(2), dists_squ = double(2),
                 c0 = double(0), c1 = double(0),
                 c2 = double(0), c3 = double(0)) {
    returnType(double(2))
    n1 <- dim(dists_abs)[1]
    n2 <- dim(dists_abs)[2]
    result <- matrix(nrow = n1, ncol = n2, init = FALSE)
    for(i in 1:n1)
      for(j in 1:n2)
        result[i, j] <- c0*c1^dists_squ[i,j]+
      c2*c3^(sin(pi*dists_abs[i,j]/12))^2
    for(i in 1:n1)
      result[i, i] <- result[i, i] + 1e-6
    return(result)
  })
cExpcov <- compileNimble(expcov)
c=4  # number of kernel params

X=100
nmiss=4
loglik <- matrix(0,X,nmiss)
correct_rate <- matrix(0,X,nmiss)
for(x in 1:X){
  N=50
  Time = 9
  K=10
  J=5
  s=matrix(0,N,Time)
  for(t in 1:Time){
    s[,t] <- seq(0.4,20,length=N)
  }
  y=matrix(0,N,Time)
  for(t in 1:Time){
    y[,t][s[,t] > 5+5/8*(t-1) & s[,t] < 10+5/8*(t-1)] <- rnorm(sum(s[,t] > 5+5/8*(t-1) & s[,t] < 10+5/8*(t-1)),mean=20,sd=1)
    y[,t][s[,t] < 5+5/8*(t-1) | s[,t] > 10+5/8*(t-1)] <- rnorm(sum(s[,t] < 5+5/8*(t-1) | s[,t] > 10+5/8*(t-1)),mean=1,sd=1)
  }
  
  new_time = sample(2:8,size=nmiss,replace=FALSE)
  new_time = sort(new_time)
  train_time <- (1:Time)[-new_time]
  s_true <- s[,new_time]
  s_true <- s_true %>% as.matrix
  y_true <- y[,new_time]
  y_true <- y_true %>% as.matrix
  s <- s[,-new_time]
  y <- y[,-new_time]
  Time <- ncol(y)
  nmiss2=Time+nmiss
  
  # model define
  dists_abs <- as.matrix(dist(c(1:(Time+nmiss))[-new_time]))
  dists_squ <- as.matrix(dist(c(1:(Time+nmiss))[-new_time]))^2
  s_tilde <- seq(4,20,length=J)
  code <- nimbleCode({
    for(k in 1:K){
      # mu ~ GP
      mu_0[k] ~ dnorm(mean=0, var=1000)
      c0_mu[k] ~ dinvgamma(shape=1e-3,scale=1e-3)
      c1_mu[k] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
      c2_mu[k] ~ dinvgamma(shape=1e-3,scale=1e-3)
      c3_mu[k] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
      cov_GP_mu[k,1:Time, 1:Time] <- expcov(dists_abs=dists_abs[1:Time,1:Time],
                                            dists_squ=dists_squ[1:Time,1:Time],
                                            c0=c0_mu[k],c1=c1_mu[k],
                                            c2=c2_mu[k],c3=c3_mu[k])
      mu_mu[k,1:Time] <- mu_0[k]*ones[1:Time]
      mu[k,1:Time] ~ dmnorm(mu_mu[k,1:Time],cov=cov_GP_mu[k,1:Time,1:Time])
      
      
      # beta0 ~ GP
      beta0_0[k] ~ dnorm(mean=0, var=1000)
      c0_beta0[k] ~ dinvgamma(shape=1e-3,scale=1e-3)
      c1_beta0[k] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
      c2_beta0[k] ~ dinvgamma(shape=1e-3,scale=1e-3)
      c3_beta0[k] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
      cov_GP_beta0[k,1:Time, 1:Time] <- expcov(dists_abs=dists_abs[1:Time,1:Time],
                                               dists_squ=dists_squ[1:Time,1:Time],
                                               c0=c0_beta0[k],c1=c1_beta0[k],
                                               c2=c2_beta0[k],c3=c3_beta0[k])
      mu_beta0[k,1:Time] <- beta0_0[k]*ones[1:Time]
      beta0[k,1:Time] ~ dmnorm(mu_beta0[k,1:Time],
                               cov=cov_GP_beta0[k,1:Time,1:Time])
      
      # beta1 ~ GP
      for(j in 1:J){
        beta1_0[k,j] ~ dnorm(mean=0, var=1000)
        c0_beta1[k,j] ~ dinvgamma(shape=1e-3,scale=1e-3)
        c1_beta1[k,j] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
        c2_beta1[k,j] ~ dinvgamma(shape=1e-3,scale=1e-3)
        c3_beta1[k,j] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
        cov_GP_beta1[k,j,1:Time, 1:Time] <- expcov(dists_abs=dists_abs[1:Time,1:Time],
                                                   dists_squ=dists_squ[1:Time,1:Time],
                                                   c0=c0_beta1[k,j],c1=c1_beta1[k,j],
                                                   c2=c2_beta1[k,j],c3=c3_beta1[k,j])
        mu_beta1[k,j,1:Time] <- beta1_0[k,j]*ones[1:Time]
        beta1[k,j,1:Time] ~ dmnorm(mu_beta1[k,j,1:Time],
                                   cov=cov_GP_beta1[k,j,1:Time,1:Time])
      }
      psi[k] ~ dunif(min=0.05,max=5)
    }
    
    # sigma2 ~ GP
    sigma2_0 ~ dnorm(mean=0, var=1000)
    c0_sigma2 ~ dinvgamma(shape=1e-3,scale=1e-3)
    c1_sigma2 ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
    c2_sigma2 ~ dinvgamma(shape=1e-3,scale=1e-3)
    c3_sigma2 ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
    cov_GP_sigma2[1:Time, 1:Time] <- expcov(dists_abs=dists_abs[1:Time,1:Time],
                                            dists_squ=dists_squ[1:Time,1:Time],
                                            c0=c0_sigma2,c1=c1_sigma2,
                                            c2=c2_sigma2,c3=c3_sigma2)
    mu_sigma2[1:Time] <- sigma2_0*ones[1:Time]
    log(sigma2[1:Time]) ~ dmnorm(mu_sigma2[1:Time],
                                 cov=cov_GP_sigma2[1:Time,1:Time])
    
    for(i in 1:N){           # LSBP
      for(k in 1:K){
        for(j in 1:J){
          kernel[i,j,k] <- exp(-sqrt(sum((s[i,1]-s_tilde[j])^2))/psi[k])
        }
      }
    }
    for(i in 1:N){
      for(k in 1:K){
        for(t in 1:Time){
          for(j in 1:J){
            g[i,j,k,t] <- beta1[k,j,t]*kernel[i,j,k]
          }
          sum_g[i,k,t] <- sum(g[i,1:J,k,t]) + beta0[k,t]
          p[i,k,t] <- 1/(1+exp(-sum_g[i,k,t]))
        }
      }
    }
    for(i in 1:N){
      for(t in 1:Time){
        w[i,1:K,t] <- stick_breaking(p[i,1:(K-1),t])
        z[i,t] ~ dcat(w[i,1:K,t])
        y[i,t] ~ dnorm(mean=mu[z[i,t],t],
                       sd=sqrt(sigma2[t])) ## likelihood
      }
    }
  })
  
  data=list(y=y)
  constants=list(s=s,s_tilde=s_tilde,dists_abs=dists_abs,dists_squ=dists_squ,
                 K=K,N=N,Time=Time,J=J,ones=rep(1,Time))
  
  # runMCMC
  
  inits=list(mu=matrix(rnorm(K*Time,mean=0,sd=1),K,Time),
             log_sigma2=rnorm(Time,mean=0,sd=1),
             beta0=matrix(rnorm(K*Time,mean=0,sd=1),K,Time),
             beta1=array(rnorm(K*J*Time,mean=0,sd=1),c(K,J,Time)),
             
             mu_0=rnorm(K,mean=0,sd=1),
             sigma2_0=rnorm(1,mean=0,sd=1),
             beta0_0=rnorm(K,mean=0,sd=1),
             beta1_0=matrix(rnorm(K*J,mean=0,sd=1),K,J),
             
             c0_mu=rinvgamma(K,shape=10,scale=10),
             c1_mu=rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             c2_mu=rinvgamma(K,shape=10,scale=10),
             c3_mu=rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             
             c0_sigma2=rinvgamma(1,shape=10,scale=10),
             c1_sigma2=rtruncnorm(1,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             c2_sigma2=rinvgamma(1,shape=10,scale=10),
             c3_sigma2=rtruncnorm(1,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             
             c0_beta0=rinvgamma(K,shape=10,scale=10),
             c1_beta0=rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             c2_beta0=rinvgamma(K,shape=10,scale=10),
             c3_beta0=rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             
             c0_beta1=matrix(rinvgamma(K*J,shape=10,scale=10),K,J),
             c1_beta1=matrix(rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),K,J),
             c2_beta1=matrix(rinvgamma(K*J,shape=10,scale=10),K,J),
             c3_beta1=matrix(rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),K,J),
             
             z=matrix(sample(1:5,N*Time,replace=T),N,Time),
             psi=runif(K,min=0.05,max=5))
  
  monitors=c('mu','sigma2','beta0','beta1',
             'mu_0','sigma2_0','beta0_0','beta1_0',
             'c0_mu','c1_mu','c2_mu','c3_mu',
             'c0_sigma2','c1_sigma2','c2_sigma2','c3_sigma2',
             'c0_beta0','c1_beta0','c2_beta0','c3_beta0',
             'c0_beta1','c1_beta1','c2_beta1','c3_beta1',
             'z','psi','p','w')
  niter=10000
  nburnin=5000
  nmcmc <- niter-nburnin
  simulMCMC <- nimbleMCMC(code=code,constants = constants,data = data,
                          inits=inits,
                          monitors=monitors,niter=niter,nburnin=nburnin,thin=1)
  
  # seperate samples
  mu_samples <- simulMCMC[,grep('mu',colnames(simulMCMC))[(c*K+1):(c*K+K*Time)]]
  mu <- apply(mu_samples,2,median) %>% matrix(.,K,Time)
  sigma2_samples <- simulMCMC[,grep('sigma2',colnames(simulMCMC))[(c+1):(c+Time)]]
  sigma2 <- apply(sigma2_samples,2,median)
  beta0samples <- simulMCMC[,grep('beta0',colnames(simulMCMC))[1:(K*Time)]]
  beta0 <- apply(beta0samples,2,median) %>% matrix(.,K,Time)
  beta1samples <- simulMCMC[,grep('beta1',colnames(simulMCMC))[1:(K*J*Time)]]
  beta1 <- apply(beta1samples,2,median) %>% array(.,c(K,J,Time))
  
  mu_0_samples <- simulMCMC[,grep('mu_0',colnames(simulMCMC))]
  mu_0 <- apply(mu_0_samples,2,median)
  sigma2_0_samples <- simulMCMC[,grep('sigma2_0',colnames(simulMCMC))]
  sigma2_0 <- median(sigma2_0_samples)
  beta0_0_samples <- simulMCMC[,grep('beta0_0',colnames(simulMCMC))]
  beta0_0 <- apply(beta0_0_samples,2,median)
  beta1_0_samples <- simulMCMC[,grep('beta1_0',colnames(simulMCMC))]
  beta1_0 <- apply(beta1_0_samples,2,median) %>% matrix(.,K,J)
  
  c0_mu_samples <- simulMCMC[,grep('c0_mu',colnames(simulMCMC))]
  c0_mu <- apply(c0_mu_samples,2,median)
  c1_mu_samples <- simulMCMC[,grep('c1_mu',colnames(simulMCMC))]
  c1_mu <- apply(c1_mu_samples,2,median)
  c2_mu_samples <- simulMCMC[,grep('c2_mu',colnames(simulMCMC))]
  c2_mu <- apply(c2_mu_samples,2,median)
  c3_mu_samples <- simulMCMC[,grep('c3_mu',colnames(simulMCMC))]
  c3_mu <- apply(c3_mu_samples,2,median)
  
  c0_sigma2_samples <- simulMCMC[,grep('c0_sigma2',colnames(simulMCMC))]
  c0_sigma2 <- median(c0_sigma2_samples)
  c1_sigma2_samples <- simulMCMC[,grep('c1_sigma2',colnames(simulMCMC))]
  c1_sigma2 <- median(c1_sigma2_samples)
  c2_sigma2_samples <- simulMCMC[,grep('c2_sigma2',colnames(simulMCMC))]
  c2_sigma2 <- median(c2_sigma2_samples)
  c3_sigma2_samples <- simulMCMC[,grep('c3_sigma2',colnames(simulMCMC))]
  c3_sigma2 <- median(c3_sigma2_samples)
  
  c0_beta0_samples <- simulMCMC[,grep('c0_beta0',colnames(simulMCMC))]
  c0_beta0 <- apply(c0_beta0_samples,2,median)
  c1_beta0_samples <- simulMCMC[,grep('c1_beta0',colnames(simulMCMC))]
  c1_beta0 <- apply(c1_beta0_samples,2,median)
  c2_beta0_samples <- simulMCMC[,grep('c2_beta0',colnames(simulMCMC))]
  c2_beta0 <- apply(c2_beta0_samples,2,median)
  c3_beta0_samples <- simulMCMC[,grep('c3_beta0',colnames(simulMCMC))]
  c3_beta0 <- apply(c3_beta0_samples,2,median)
  
  c0_beta1_samples <- simulMCMC[,grep('c0_beta1',colnames(simulMCMC))]
  c0_beta1 <- apply(c0_beta1_samples,2,median) %>% matrix(.,K,J)
  c1_beta1_samples <- simulMCMC[,grep('c1_beta1',colnames(simulMCMC))]
  c1_beta1 <- apply(c1_beta1_samples,2,median) %>% matrix(.,K,J)
  c2_beta1_samples <- simulMCMC[,grep('c2_beta1',colnames(simulMCMC))]
  c2_beta1 <- apply(c2_beta1_samples,2,median) %>% matrix(.,K,J)
  c3_beta1_samples <- simulMCMC[,grep('c3_beta1',colnames(simulMCMC))]
  c3_beta1 <- apply(c3_beta1_samples,2,median) %>% matrix(.,K,J)
  
  z_samples <- simulMCMC[,grep('z',colnames(simulMCMC))]
  z <- apply(z_samples,2,function(x) table(x) %>% which.max %>% names) %>%
    as.integer %>% matrix(.,N,Time)
  psi_samples <- simulMCMC[,grep('psi',colnames(simulMCMC))]
  psi <- apply(psi_samples,2,median)
  
  y_pred1 <- matrix(0,N,nmiss)
  sigma2_pred1 <- c()
  for(new in 1:nmiss){
    times9_1_ind <- new_time[new]-1*(new-1)
    times1_1_ind <- new_time[new]-1*new
    times9_1 <- train_time[times9_1_ind]
    times1_1 <- train_time[times1_1_ind]
    beta1_pred_1 <- (beta1[,,times1_1_ind]*(times9_1-new_time[new]) +
                       beta1[,,times9_1_ind]*(new_time[new]-times1_1))/(times9_1-times1_1)
    
    beta0_pred_1 <- (beta0[,times1_1_ind]*(times9_1-new_time[new]) +
                       beta0[,times9_1_ind]*(new_time[new]-times1_1))/(times9_1-times1_1)
    
    mu_pred1 <- (mu[,times1_1_ind]*(times9_1-new_time[new]) +
                   mu[,times9_1_ind]*(new_time[new]-times1_1))/(times9_1-times1_1)
    
    kernel <- g <- array(0,c(N,K,J))
    sum_g <- p <- w <- matrix(0,N,K)
    for(i in 1:N){
      for(k in 1:K){
        for(j in 1:J){
          kernel[i,k,j] <- exp(-sqrt(sum((s[i,1]-s_tilde[j])^2))/psi[k])
        }
      }
    }
    for(i in 1:N){
      for(k in 1:K){
        for(j in 1:J){
          g[i,k,j] <- beta1_pred_1[k,j] * kernel[i,k,j]
        }
      }
    }
    for(i in 1:N){
      for(k in 1:K){
        sum_g[i,k] <- sum(g[i,k,1:J]) + beta0_pred_1[k]
        p[i,k] <- 1/(1+exp(-sum_g[i,k]))
      }
      w[i,] <- stick_breaking(p[i,][1:(K-1)])
    }
    # y_pred1 <- c()
    for(i in 1:N){
      y_pred1[i,new] <- sum(w[i,] * mu_pred1)
    }
    sigma2_pred1[new] <- (sigma2[times1_1_ind]*(times9_1-new_time[new]) +
                            sigma2[times9_1_ind]*(new_time[new]-times1_1))/(times9_1-times1_1)
    
  }
  par(mfrow=c(1,nmiss))
  
  for(t in 1:nmiss){
    plot(s[,1],y_true[,t],ylim=c(-1,30),xlab='',ylab='',main=new_time[t])
    par(new=T)
    plot(s[,1],y_pred1[,t],ylim=c(-1,30),xlab='',ylab='',col='red',pch=3)
    par(new=F)
    loglik[x,t] <- mean(dnorm(x = y_true[,t],mean = y_pred1[,t],sd=sqrt(sigma2_pred1[t]),log=T))
  }
  par(mfrow=c(1,1))
  true_group=matrix(0,N,nmiss)
  for(t in 1:nmiss){
    true_group[,t] <-
      (s_true[,t] > 5+5/8*(new_time[t]-1) & s_true[,t] < 10+5/8*(new_time[t]-1))
    
  }
  for(t in 1:nmiss){
    correct_rate[x,t] <-
      sum(true_group[,t] ==
            (y_pred1[,t] > 10))/N
    
  }
}
loglik %>% mean
correct_rate %>% mean
#### simulation nmiss = 5 #### 
library(nimble)
library(dplyr)
library(truncnorm)
# GP kernel
expcov <- nimbleFunction(
  run = function(dists_abs = double(2), dists_squ = double(2),
                 c0 = double(0), c1 = double(0),
                 c2 = double(0), c3 = double(0)) {
    returnType(double(2))
    n1 <- dim(dists_abs)[1]
    n2 <- dim(dists_abs)[2]
    result <- matrix(nrow = n1, ncol = n2, init = FALSE)
    for(i in 1:n1)
      for(j in 1:n2)
        result[i, j] <- c0*c1^dists_squ[i,j]+
      c2*c3^(sin(pi*dists_abs[i,j]/12))^2
    for(i in 1:n1)
      result[i, i] <- result[i, i] + 1e-6
    return(result)
  })
cExpcov <- compileNimble(expcov)
c=4  # number of kernel params

X=100
nmiss=5
loglik <- matrix(0,X,nmiss)
correct_rate <- matrix(0,X,nmiss)
for(x in 1:X){
  N=50
  Time = 9
  K=10
  J=5
  s=matrix(0,N,Time)
  for(t in 1:Time){
    s[,t] <- seq(0.4,20,length=N)
  }
  y=matrix(0,N,Time)
  for(t in 1:Time){
    y[,t][s[,t] > 5+5/8*(t-1) & s[,t] < 10+5/8*(t-1)] <- rnorm(sum(s[,t] > 5+5/8*(t-1) & s[,t] < 10+5/8*(t-1)),mean=20,sd=1)
    y[,t][s[,t] < 5+5/8*(t-1) | s[,t] > 10+5/8*(t-1)] <- rnorm(sum(s[,t] < 5+5/8*(t-1) | s[,t] > 10+5/8*(t-1)),mean=1,sd=1)
  }
  
  new_time = sample(2:8,size=nmiss,replace=FALSE)
  new_time = sort(new_time)
  train_time <- (1:Time)[-new_time]
  s_true <- s[,new_time]
  s_true <- s_true %>% as.matrix
  y_true <- y[,new_time]
  y_true <- y_true %>% as.matrix
  s <- s[,-new_time]
  y <- y[,-new_time]
  Time <- ncol(y)
  nmiss2=Time+nmiss
  
  # model define
  dists_abs <- as.matrix(dist(c(1:(Time+nmiss))[-new_time]))
  dists_squ <- as.matrix(dist(c(1:(Time+nmiss))[-new_time]))^2
  s_tilde <- seq(4,20,length=J)
  code <- nimbleCode({
    for(k in 1:K){
      # mu ~ GP
      mu_0[k] ~ dnorm(mean=0, var=1000)
      c0_mu[k] ~ dinvgamma(shape=1e-3,scale=1e-3)
      c1_mu[k] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
      c2_mu[k] ~ dinvgamma(shape=1e-3,scale=1e-3)
      c3_mu[k] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
      cov_GP_mu[k,1:Time, 1:Time] <- expcov(dists_abs=dists_abs[1:Time,1:Time],
                                            dists_squ=dists_squ[1:Time,1:Time],
                                            c0=c0_mu[k],c1=c1_mu[k],
                                            c2=c2_mu[k],c3=c3_mu[k])
      mu_mu[k,1:Time] <- mu_0[k]*ones[1:Time]
      mu[k,1:Time] ~ dmnorm(mu_mu[k,1:Time],cov=cov_GP_mu[k,1:Time,1:Time])
      
      
      # beta0 ~ GP
      beta0_0[k] ~ dnorm(mean=0, var=1000)
      c0_beta0[k] ~ dinvgamma(shape=1e-3,scale=1e-3)
      c1_beta0[k] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
      c2_beta0[k] ~ dinvgamma(shape=1e-3,scale=1e-3)
      c3_beta0[k] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
      cov_GP_beta0[k,1:Time, 1:Time] <- expcov(dists_abs=dists_abs[1:Time,1:Time],
                                               dists_squ=dists_squ[1:Time,1:Time],
                                               c0=c0_beta0[k],c1=c1_beta0[k],
                                               c2=c2_beta0[k],c3=c3_beta0[k])
      mu_beta0[k,1:Time] <- beta0_0[k]*ones[1:Time]
      beta0[k,1:Time] ~ dmnorm(mu_beta0[k,1:Time],
                               cov=cov_GP_beta0[k,1:Time,1:Time])
      
      # beta1 ~ GP
      for(j in 1:J){
        beta1_0[k,j] ~ dnorm(mean=0, var=1000)
        c0_beta1[k,j] ~ dinvgamma(shape=1e-3,scale=1e-3)
        c1_beta1[k,j] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
        c2_beta1[k,j] ~ dinvgamma(shape=1e-3,scale=1e-3)
        c3_beta1[k,j] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
        cov_GP_beta1[k,j,1:Time, 1:Time] <- expcov(dists_abs=dists_abs[1:Time,1:Time],
                                                   dists_squ=dists_squ[1:Time,1:Time],
                                                   c0=c0_beta1[k,j],c1=c1_beta1[k,j],
                                                   c2=c2_beta1[k,j],c3=c3_beta1[k,j])
        mu_beta1[k,j,1:Time] <- beta1_0[k,j]*ones[1:Time]
        beta1[k,j,1:Time] ~ dmnorm(mu_beta1[k,j,1:Time],
                                   cov=cov_GP_beta1[k,j,1:Time,1:Time])
      }
      psi[k] ~ dunif(min=0.05,max=5)
    }
    
    # sigma2 ~ GP
    sigma2_0 ~ dnorm(mean=0, var=1000)
    c0_sigma2 ~ dinvgamma(shape=1e-3,scale=1e-3)
    c1_sigma2 ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
    c2_sigma2 ~ dinvgamma(shape=1e-3,scale=1e-3)
    c3_sigma2 ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
    cov_GP_sigma2[1:Time, 1:Time] <- expcov(dists_abs=dists_abs[1:Time,1:Time],
                                            dists_squ=dists_squ[1:Time,1:Time],
                                            c0=c0_sigma2,c1=c1_sigma2,
                                            c2=c2_sigma2,c3=c3_sigma2)
    mu_sigma2[1:Time] <- sigma2_0*ones[1:Time]
    log(sigma2[1:Time]) ~ dmnorm(mu_sigma2[1:Time],
                                 cov=cov_GP_sigma2[1:Time,1:Time])
    
    for(i in 1:N){           # LSBP
      for(k in 1:K){
        for(j in 1:J){
          kernel[i,j,k] <- exp(-sqrt(sum((s[i,1]-s_tilde[j])^2))/psi[k])
        }
      }
    }
    for(i in 1:N){
      for(k in 1:K){
        for(t in 1:Time){
          for(j in 1:J){
            g[i,j,k,t] <- beta1[k,j,t]*kernel[i,j,k]
          }
          sum_g[i,k,t] <- sum(g[i,1:J,k,t]) + beta0[k,t]
          p[i,k,t] <- 1/(1+exp(-sum_g[i,k,t]))
        }
      }
    }
    for(i in 1:N){
      for(t in 1:Time){
        w[i,1:K,t] <- stick_breaking(p[i,1:(K-1),t])
        z[i,t] ~ dcat(w[i,1:K,t])
        y[i,t] ~ dnorm(mean=mu[z[i,t],t],
                       sd=sqrt(sigma2[t])) ## likelihood
      }
    }
  })
  
  data=list(y=y)
  constants=list(s=s,s_tilde=s_tilde,dists_abs=dists_abs,dists_squ=dists_squ,
                 K=K,N=N,Time=Time,J=J,ones=rep(1,Time))
  
  # runMCMC
  
  inits=list(mu=matrix(rnorm(K*Time,mean=0,sd=1),K,Time),
             log_sigma2=rnorm(Time,mean=0,sd=1),
             beta0=matrix(rnorm(K*Time,mean=0,sd=1),K,Time),
             beta1=array(rnorm(K*J*Time,mean=0,sd=1),c(K,J,Time)),
             
             mu_0=rnorm(K,mean=0,sd=1),
             sigma2_0=rnorm(1,mean=0,sd=1),
             beta0_0=rnorm(K,mean=0,sd=1),
             beta1_0=matrix(rnorm(K*J,mean=0,sd=1),K,J),
             
             c0_mu=rinvgamma(K,shape=10,scale=10),
             c1_mu=rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             c2_mu=rinvgamma(K,shape=10,scale=10),
             c3_mu=rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             
             c0_sigma2=rinvgamma(1,shape=10,scale=10),
             c1_sigma2=rtruncnorm(1,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             c2_sigma2=rinvgamma(1,shape=10,scale=10),
             c3_sigma2=rtruncnorm(1,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             
             c0_beta0=rinvgamma(K,shape=10,scale=10),
             c1_beta0=rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             c2_beta0=rinvgamma(K,shape=10,scale=10),
             c3_beta0=rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
             
             c0_beta1=matrix(rinvgamma(K*J,shape=10,scale=10),K,J),
             c1_beta1=matrix(rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),K,J),
             c2_beta1=matrix(rinvgamma(K*J,shape=10,scale=10),K,J),
             c3_beta1=matrix(rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),K,J),
             
             z=matrix(sample(1:5,N*Time,replace=T),N,Time),
             psi=runif(K,min=0.05,max=5))
  
  monitors=c('mu','sigma2','beta0','beta1',
             'mu_0','sigma2_0','beta0_0','beta1_0',
             'c0_mu','c1_mu','c2_mu','c3_mu',
             'c0_sigma2','c1_sigma2','c2_sigma2','c3_sigma2',
             'c0_beta0','c1_beta0','c2_beta0','c3_beta0',
             'c0_beta1','c1_beta1','c2_beta1','c3_beta1',
             'z','psi','p','w')
  niter=10000
  nburnin=5000
  nmcmc <- niter-nburnin
  simulMCMC <- nimbleMCMC(code=code,constants = constants,data = data,
                          inits=inits,
                          monitors=monitors,niter=niter,nburnin=nburnin,thin=1)
  
  # seperate samples
  mu_samples <- simulMCMC[,grep('mu',colnames(simulMCMC))[(c*K+1):(c*K+K*Time)]]
  mu <- apply(mu_samples,2,median) %>% matrix(.,K,Time)
  sigma2_samples <- simulMCMC[,grep('sigma2',colnames(simulMCMC))[(c+1):(c+Time)]]
  sigma2 <- apply(sigma2_samples,2,median)
  beta0samples <- simulMCMC[,grep('beta0',colnames(simulMCMC))[1:(K*Time)]]
  beta0 <- apply(beta0samples,2,median) %>% matrix(.,K,Time)
  beta1samples <- simulMCMC[,grep('beta1',colnames(simulMCMC))[1:(K*J*Time)]]
  beta1 <- apply(beta1samples,2,median) %>% array(.,c(K,J,Time))
  
  mu_0_samples <- simulMCMC[,grep('mu_0',colnames(simulMCMC))]
  mu_0 <- apply(mu_0_samples,2,median)
  sigma2_0_samples <- simulMCMC[,grep('sigma2_0',colnames(simulMCMC))]
  sigma2_0 <- median(sigma2_0_samples)
  beta0_0_samples <- simulMCMC[,grep('beta0_0',colnames(simulMCMC))]
  beta0_0 <- apply(beta0_0_samples,2,median)
  beta1_0_samples <- simulMCMC[,grep('beta1_0',colnames(simulMCMC))]
  beta1_0 <- apply(beta1_0_samples,2,median) %>% matrix(.,K,J)
  
  c0_mu_samples <- simulMCMC[,grep('c0_mu',colnames(simulMCMC))]
  c0_mu <- apply(c0_mu_samples,2,median)
  c1_mu_samples <- simulMCMC[,grep('c1_mu',colnames(simulMCMC))]
  c1_mu <- apply(c1_mu_samples,2,median)
  c2_mu_samples <- simulMCMC[,grep('c2_mu',colnames(simulMCMC))]
  c2_mu <- apply(c2_mu_samples,2,median)
  c3_mu_samples <- simulMCMC[,grep('c3_mu',colnames(simulMCMC))]
  c3_mu <- apply(c3_mu_samples,2,median)
  
  c0_sigma2_samples <- simulMCMC[,grep('c0_sigma2',colnames(simulMCMC))]
  c0_sigma2 <- median(c0_sigma2_samples)
  c1_sigma2_samples <- simulMCMC[,grep('c1_sigma2',colnames(simulMCMC))]
  c1_sigma2 <- median(c1_sigma2_samples)
  c2_sigma2_samples <- simulMCMC[,grep('c2_sigma2',colnames(simulMCMC))]
  c2_sigma2 <- median(c2_sigma2_samples)
  c3_sigma2_samples <- simulMCMC[,grep('c3_sigma2',colnames(simulMCMC))]
  c3_sigma2 <- median(c3_sigma2_samples)
  
  c0_beta0_samples <- simulMCMC[,grep('c0_beta0',colnames(simulMCMC))]
  c0_beta0 <- apply(c0_beta0_samples,2,median)
  c1_beta0_samples <- simulMCMC[,grep('c1_beta0',colnames(simulMCMC))]
  c1_beta0 <- apply(c1_beta0_samples,2,median)
  c2_beta0_samples <- simulMCMC[,grep('c2_beta0',colnames(simulMCMC))]
  c2_beta0 <- apply(c2_beta0_samples,2,median)
  c3_beta0_samples <- simulMCMC[,grep('c3_beta0',colnames(simulMCMC))]
  c3_beta0 <- apply(c3_beta0_samples,2,median)
  
  c0_beta1_samples <- simulMCMC[,grep('c0_beta1',colnames(simulMCMC))]
  c0_beta1 <- apply(c0_beta1_samples,2,median) %>% matrix(.,K,J)
  c1_beta1_samples <- simulMCMC[,grep('c1_beta1',colnames(simulMCMC))]
  c1_beta1 <- apply(c1_beta1_samples,2,median) %>% matrix(.,K,J)
  c2_beta1_samples <- simulMCMC[,grep('c2_beta1',colnames(simulMCMC))]
  c2_beta1 <- apply(c2_beta1_samples,2,median) %>% matrix(.,K,J)
  c3_beta1_samples <- simulMCMC[,grep('c3_beta1',colnames(simulMCMC))]
  c3_beta1 <- apply(c3_beta1_samples,2,median) %>% matrix(.,K,J)
  
  z_samples <- simulMCMC[,grep('z',colnames(simulMCMC))]
  z <- apply(z_samples,2,function(x) table(x) %>% which.max %>% names) %>%
    as.integer %>% matrix(.,N,Time)
  psi_samples <- simulMCMC[,grep('psi',colnames(simulMCMC))]
  psi <- apply(psi_samples,2,median)
  
  y_pred1 <- matrix(0,N,nmiss)
  sigma2_pred1 <- c()
  for(new in 1:nmiss){
    times9_1_ind <- new_time[new]-1*(new-1)
    times1_1_ind <- new_time[new]-1*new
    times9_1 <- train_time[times9_1_ind]
    times1_1 <- train_time[times1_1_ind]
    beta1_pred_1 <- (beta1[,,times1_1_ind]*(times9_1-new_time[new]) +
                       beta1[,,times9_1_ind]*(new_time[new]-times1_1))/(times9_1-times1_1)
    
    beta0_pred_1 <- (beta0[,times1_1_ind]*(times9_1-new_time[new]) +
                       beta0[,times9_1_ind]*(new_time[new]-times1_1))/(times9_1-times1_1)
    
    mu_pred1 <- (mu[,times1_1_ind]*(times9_1-new_time[new]) +
                   mu[,times9_1_ind]*(new_time[new]-times1_1))/(times9_1-times1_1)
    
    kernel <- g <- array(0,c(N,K,J))
    sum_g <- p <- w <- matrix(0,N,K)
    for(i in 1:N){
      for(k in 1:K){
        for(j in 1:J){
          kernel[i,k,j] <- exp(-sqrt(sum((s[i,1]-s_tilde[j])^2))/psi[k])
        }
      }
    }
    for(i in 1:N){
      for(k in 1:K){
        for(j in 1:J){
          g[i,k,j] <- beta1_pred_1[k,j] * kernel[i,k,j]
        }
      }
    }
    for(i in 1:N){
      for(k in 1:K){
        sum_g[i,k] <- sum(g[i,k,1:J]) + beta0_pred_1[k]
        p[i,k] <- 1/(1+exp(-sum_g[i,k]))
      }
      w[i,] <- stick_breaking(p[i,][1:(K-1)])
    }
    # y_pred1 <- c()
    for(i in 1:N){
      y_pred1[i,new] <- sum(w[i,] * mu_pred1)
    }
    sigma2_pred1[new] <- (sigma2[times1_1_ind]*(times9_1-new_time[new]) +
                            sigma2[times9_1_ind]*(new_time[new]-times1_1))/(times9_1-times1_1)
    
  }
  par(mfrow=c(1,nmiss))
  
  for(t in 1:nmiss){
    plot(s[,1],y_true[,t],ylim=c(-1,30),xlab='',ylab='',main=new_time[t])
    par(new=T)
    plot(s[,1],y_pred1[,t],ylim=c(-1,30),xlab='',ylab='',col='red',pch=3)
    par(new=F)
    loglik[x,t] <- mean(dnorm(x = y_true[,t],mean = y_pred1[,t],sd=sqrt(sigma2_pred1[t]),log=T))
  }
  par(mfrow=c(1,1))
  true_group=matrix(0,N,nmiss)
  for(t in 1:nmiss){
    true_group[,t] <-
      (s_true[,t] > 5+5/8*(new_time[t]-1) & s_true[,t] < 10+5/8*(new_time[t]-1))
    
  }
  for(t in 1:nmiss){
    correct_rate[x,t] <-
      sum(true_group[,t] ==
            (y_pred1[,t] > 10))/N
    
  }
}
loglik %>% mean
correct_rate %>% mean


#### Figure4 (kernel center map) ####
library(stringr)
library(progress)
library(dplyr)
library(ggplot2)
library(rgdal)      ## retired package...

map <- readOGR('SIG_20221119/sig.shp',encoding = 'euc-kr')
new_map <- fortify(map,region = 'SIG_CD')
kernel_center <- read.csv('misemise/kernel_center2.csv',header=T,fileEncoding = 'euc-kr')
colnames(kernel_center) <- c('X','long','lat')
coord <- data.frame(utmk.long=kernel_center$long, utmk.lat=kernel_center$lat)
convertCoordSystem <- function(long, lat, from.crs, to.crs){
  xy <- data.frame(long=long, lat=lat)
  coordinates(xy) <- ~long+lat
  
  from.crs <- CRS(from.crs)
  from.coordinates <- SpatialPoints(xy, proj4string=from.crs)
  
  to.crs <- CRS(to.crs)
  changed <- as.data.frame(SpatialPoints(spTransform(from.coordinates, to.crs)))
  names(changed) <- c("long2", "lat2")
  
  return(changed)
}
to.crs = "+proj=tmerc +lat_0=38 +lon_0=127.5 +k=0.9996 +x_0=1000000 +y_0=2000000 +ellps=GRS80 +units=m +no_defs"
from.crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
coord <- cbind(coord, convertCoordSystem(coord$utmk.long, coord$utmk.lat, from.crs, to.crs))
kernel_center <- cbind(kernel_center,coord[,3:4])
map_korea_ggplot <- new_map %>% 
  ggplot(aes(x=long, y=lat, group = group))
map_korea_ggplot+
  geom_polygon(fill='white', color='black')+
  coord_quickmap()+
  geom_point(data=kernel_center,aes(x=long2,y=lat2),size=3,col='red',group=11110.1)


#### compare other bayesian method ####
pm2.5_dat <- read.csv(file = 'misemise/pm2.5_dat_result16.csv',header=T)
library(ggplot2)
library(dplyr)
library(spTimer)
colnames(pm2.5_dat)
Time=43
pm2.5_dat[,4+Time+2]
library(verification)
num <- 1
crps_my <-crps_bms<- c()

for(num in 1:(Time-2)){
  crps_my[num] <- verification::crps(pm2.5_dat[,(Time+5-num)],
                                     pm2.5_dat[,c((2*Time+5-num),(3*Time+5-num))])$CRPS
}
for(num in 1:Time){
  crps_bms[num] <- verification::crps(pm2.5_dat2[,(Time+5-num)],
                                      pm2.5_dat2[,c((2*Time+5-num),(3*Time+5-num))])$CRPS
}
mean(crps_my)
mean(crps_bms)

spT.validation(unlist(c(pm2.5_dat[,c(5:(5+Time-1))])),
               unlist(c(pm2.5_dat[,Time+c(5:(5+Time-1))])))  ### fitted my
spT.validation(unlist(c(pm2.5_dat2[,c(5:(5+Time-1))])),
               unlist(c(pm2.5_dat2[,Time+c(5:(5+Time-1))])))




#### PM2.5 data analysis ####
## data generate ##
library(nimble)
library(dplyr)
library(truncnorm)
## GP kernel ##
expcov <- nimbleFunction(
  run = function(dists_abs = double(2), dists_squ = double(2),
                 c0 = double(0), c1 = double(0),
                 c2 = double(0), c3 = double(0)) {
    returnType(double(2))
    n1 <- dim(dists_abs)[1]
    n2 <- dim(dists_abs)[2]
    result <- matrix(nrow = n1, ncol = n2, init = FALSE)
    for(i in 1:n1)
      for(j in 1:n2)
        result[i, j] <- c0*c1^dists_squ[i,j]+
      c2*c3^(sin(pi*dists_abs[i,j]/12))^2
    for(i in 1:n1)
      result[i, i] <- result[i, i] + 1e-6
    return(result)
  })
cExpcov <- compileNimble(expcov)
c=4  # number of kernel params
nmiss=0
N=162

pm2.5_dat <- read.csv('misemise/pm2.5_dat2.csv',heade=T,fileEncoding = 'euc-kr')

colnames(pm2.5_dat)
pm2.5_dat <- arrange(pm2.5_dat,CTPRV_NM)
pm2.5_dat$CTPRV_NM
# View(pm2.5_dat)
y_full <- pm2.5_dat[1:N,5:ncol(pm2.5_dat)]
s <- pm2.5_dat[1:N,1:2]
y <- pm2.5_dat[1:N,5:(ncol(pm2.5_dat)-nmiss)]
y <- as.matrix(y)
mean_y <- mean(y)
sd_y <- sd(y)
y <- (y-mean(y))/sd(y)

Time = ncol(y)
K=20
s_tilde <- read.csv('misemise/kernel_center2.csv',header=T,
                    fileEncoding = 'euc-kr',row.names=1)
J <- nrow(s_tilde)

## model define ##
dists_abs <- as.matrix(dist(1:Time))
dists_squ <- as.matrix(dist(1:Time))^2
code <- nimbleCode({
  for(k in 1:K){
    # mu ~ GP
    mu_0[k] ~ dnorm(mean=0, var=1000)
    c0_mu[k] ~ dinvgamma(shape=1e-3,scale=1e-3)
    c1_mu[k] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
    c2_mu[k] ~ dinvgamma(shape=1e-3,scale=1e-3)
    c3_mu[k] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
    cov_GP_mu[k,1:Time, 1:Time] <- expcov(dists_abs=dists_abs[1:Time,1:Time],
                                          dists_squ=dists_squ[1:Time,1:Time],
                                          c0=c0_mu[k],c1=c1_mu[k],
                                          c2=c2_mu[k],c3=c3_mu[k])
    mu_mu[k,1:Time] <- mu_0[k]*ones[1:Time]
    mu[k,1:Time] ~ dmnorm(mu_mu[k,1:Time],cov=cov_GP_mu[k,1:Time,1:Time])


    # beta0 ~ GP
    beta0_0[k] ~ dnorm(mean=0, var=1000)
    c0_beta0[k] ~ dinvgamma(shape=1e-3,scale=1e-3)
    c1_beta0[k] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
    c2_beta0[k] ~ dinvgamma(shape=1e-3,scale=1e-3)
    c3_beta0[k] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
    cov_GP_beta0[k,1:Time, 1:Time] <- expcov(dists_abs=dists_abs[1:Time,1:Time],
                                             dists_squ=dists_squ[1:Time,1:Time],
                                             c0=c0_beta0[k],c1=c1_beta0[k],
                                             c2=c2_beta0[k],c3=c3_beta0[k])
    mu_beta0[k,1:Time] <- beta0_0[k]*ones[1:Time]
    beta0[k,1:Time] ~ dmnorm(mu_beta0[k,1:Time],
                             cov=cov_GP_beta0[k,1:Time,1:Time])

    # beta1 ~ GP
    for(j in 1:J){
      beta1_0[k,j] ~ dnorm(mean=0, var=1000)
      c0_beta1[k,j] ~ dinvgamma(shape=1e-3,scale=1e-3)
      c1_beta1[k,j] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
      c2_beta1[k,j] ~ dinvgamma(shape=1e-3,scale=1e-3)
      c3_beta1[k,j] ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
      cov_GP_beta1[k,j,1:Time, 1:Time] <- expcov(dists_abs=dists_abs[1:Time,1:Time],
                                                 dists_squ=dists_squ[1:Time,1:Time],
                                                 c0=c0_beta1[k,j],c1=c1_beta1[k,j],
                                                 c2=c2_beta1[k,j],c3=c3_beta1[k,j])
      mu_beta1[k,j,1:Time] <- beta1_0[k,j]*ones[1:Time]
      beta1[k,j,1:Time] ~ dmnorm(mu_beta1[k,j,1:Time],
                                 cov=cov_GP_beta1[k,j,1:Time,1:Time])
    }
    psi[k] ~ dunif(min=0.05,max=5)
  }

  # sigma2 ~ GP
  sigma2_0 ~ dnorm(mean=0, var=1000)
  c0_sigma2 ~ dinvgamma(shape=1e-3,scale=1e-3)
  c1_sigma2 ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
  c2_sigma2 ~ dinvgamma(shape=1e-3,scale=1e-3)
  c3_sigma2 ~ T(dnorm(mean=1/2,sd=1),1e-6,1-1e-6)
  cov_GP_sigma2[1:Time, 1:Time] <- expcov(dists_abs=dists_abs[1:Time,1:Time],
                                          dists_squ=dists_squ[1:Time,1:Time],
                                          c0=c0_sigma2,c1=c1_sigma2,
                                          c2=c2_sigma2,c3=c3_sigma2)
  mu_sigma2[1:Time] <- sigma2_0*ones[1:Time]
  log(sigma2[1:Time]) ~ dmnorm(mu_sigma2[1:Time],
                               cov=cov_GP_sigma2[1:Time,1:Time])

  for(i in 1:N){           # LSBP
    for(k in 1:K){
      for(j in 1:J){
        kernel[i,j,k] <- exp(-sqrt(sum((s[i,]-s_tilde[j,])^2))/psi[k])
      }
    }
  }
  for(i in 1:N){
    for(k in 1:K){
      for(t in 1:Time){
        for(j in 1:J){
          g[i,j,k,t] <- beta1[k,j,t]*kernel[i,j,k]
        }
        sum_g[i,k,t] <- sum(g[i,1:J,k,t]) + beta0[k,t]
        p[i,k,t] <- 1/(1+exp(-sum_g[i,k,t]))
      }
    }
  }
  for(i in 1:N){
    for(t in 1:Time){
      w[i,1:K,t] <- stick_breaking(p[i,1:(K-1),t])
      z[i,t] ~ dcat(w[i,1:K,t])
      y[i,t] ~ dnorm(mean=mu[z[i,t],t],
                     sd=sqrt(sigma2[t])) ## likelihood
    }
  }
})

data=list(y=y)
constants=list(s=s,s_tilde=s_tilde,dists_abs=dists_abs,dists_squ=dists_squ,
               K=K,N=N,Time=Time,J=J,ones=rep(1,Time))

## runMCMC ##

inits=list(mu=matrix(rnorm(K*Time,mean=0,sd=1),K,Time),
           log_sigma2=rnorm(Time,mean=0,sd=1),
           beta0=matrix(rnorm(K*Time,mean=0,sd=1),K,Time),
           beta1=array(rnorm(K*J*Time,mean=0,sd=1),c(K,J,Time)),

           mu_0=rnorm(K,mean=0,sd=1),
           sigma2_0=rnorm(1,mean=0,sd=1),
           beta0_0=rnorm(K,mean=0,sd=1),
           beta1_0=matrix(rnorm(K*J,mean=0,sd=1),K,J),

           c0_mu=rinvgamma(K,shape=10,scale=10),
           c1_mu=rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
           c2_mu=rinvgamma(K,shape=10,scale=10),
           c3_mu=rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),

           c0_sigma2=rinvgamma(1,shape=10,scale=10),
           c1_sigma2=rtruncnorm(1,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
           c2_sigma2=rinvgamma(1,shape=10,scale=10),
           c3_sigma2=rtruncnorm(1,a=1e-6,b=1-1e-6,mean=1/2,sd=1),

           c0_beta0=rinvgamma(K,shape=10,scale=10),
           c1_beta0=rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),
           c2_beta0=rinvgamma(K,shape=10,scale=10),
           c3_beta0=rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),

           c0_beta1=matrix(rinvgamma(K*J,shape=10,scale=10),K,J),
           c1_beta1=matrix(rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),K,J),
           c2_beta1=matrix(rinvgamma(K*J,shape=10,scale=10),K,J),
           c3_beta1=matrix(rtruncnorm(K,a=1e-6,b=1-1e-6,mean=1/2,sd=1),K,J),

           z=matrix(sample(1:5,N*Time,replace=T),N,Time),
           psi=runif(K,min=0.05,max=5))

monitors=c('mu','sigma2','beta0','beta1',
           'mu_0','sigma2_0','beta0_0','beta1_0',
           'c0_mu','c1_mu','c2_mu','c3_mu',
           'c0_sigma2','c1_sigma2','c2_sigma2','c3_sigma2',
           'c0_beta0','c1_beta0','c2_beta0','c3_beta0',
           'c0_beta1','c1_beta1','c2_beta1','c3_beta1',
           'z','psi','p','w')
niter=30000
nburnin=20000
model <- nimbleModel(code = code, data = data,
                     constants = constants, inits = inits)
## Ensure we have the nodes needed to simulate new datasets
dataNodes <- model$getNodeNames(dataOnly = TRUE)
parentNodes <- model$getParents(dataNodes, stochOnly = TRUE)  # `getParents` is new in nimble 0.11.0
## Ensure we have both data nodes and deterministic intermediates (e.g., lifted nodes)
simNodes <- model$getDependencies(parentNodes, self = FALSE)
cmodel  <- compileNimble(model)
mcmc    <- buildMCMC(model, monitors = parentNodes)
cmcmc   <- compileNimble(mcmc, project = model)
samples <- runMCMC(cmcmc, niter = niter, nburnin = nburnin, thin=1)
ppSamplerNF <- nimbleFunction(
  setup = function(model, mcmc) {
    dataNodes <- model$getNodeNames(dataOnly = TRUE)
    parentNodes <- model$getParents(dataNodes, stochOnly = TRUE)
    cat("Stochastic parents of data are:", paste(parentNodes, collapse = ','), ".\n")
    simNodes <- model$getDependencies(parentNodes, self = FALSE)
    vars <- mcmc$mvSamples$getVarNames()  # need ordering of variables in mvSamples / samples matrix
    cat("Using posterior samples of:", paste(vars, collapse = ','), ".\n")
    n <- length(model$expandNodeNames(dataNodes, returnScalarComponents = TRUE))
  },
  run = function(samples = double(2)) {
    nSamp <- dim(samples)[1]
    ppSamples <- matrix(nrow = nSamp, ncol = n)
    for(i in 1:nSamp) {
      values(model, vars) <<- samples[i, ]
      model$simulate(simNodes, includeData = TRUE)
      ppSamples[i, ] <- values(model, dataNodes)
    }
    returnType(double(2))
    return(ppSamples)
  })
## Create the sampler for this model and this MCMC.
ppSampler <- ppSamplerNF(model, mcmc)
cppSampler <- compileNimble(ppSampler, project = model)
## Check ordering of variables is same in 'vars' and in 'samples'.
colnames(samples)
# identical(colnames(samples), model$expandNodeNames(mcmc$mvSamples$getVarNames()))
# set.seed(1)
system.time(ppSamples_via_nf <- cppSampler$run(samples))
y_hat <- ppSamples_via_nf %>% apply(.,2,median) %>% matrix(.,nrow=162)
y_hat <- y_hat*sd_y+mean_y
sigma2_samples <- exp(samples[,grep('sigma2',colnames(samples))])
sigma2 <- apply(sigma2_samples,2,median)*sd_y^2
sigma2 <- matrix(sigma2,nrow=162,ncol=Time,byrow=T)
colnames(sigma2) <- paste0('sigma2_',1:43)
par(mfrow=c(2,4))
loglik=c()
for(n in 1:Time){
  plot(y_full[,n],ylim=c(-1,50),xlab='',ylab='',main=n)
  par(new=T)
  plot(y_hat[,n],ylim=c(-1,50),xlab='',ylab='',col='red',pch=3)
  par(new=F)
  loglik[n] <- mean(dnorm(x = y_full[,n],mean = y_hat[,n],
                          sd=sqrt(sigma2[,n]),log=T))
}
loglik %>% mean
par(mfrow=c(1,1))
sqrt(sigma2[1,]) %>% plot(main='sd',xlab='Time',ylab='value',type='l')
abline(v=c(9,21,33,45,57,69),col='red')
colnames(y_hat) <- paste0('mu',1:Time)

# saving file
pm2.5_dat_result <- cbind(pm2.5_dat,y_hat,sigma2)
# write.csv(pm2.5_dat_result,file='pm2.5_dat_result16.csv',row.names = F)
z_samples <- samples[,grep('z',colnames(samples))]
z <- apply(z_samples,2,function(x) table(x) %>% which.max %>% names) %>% as.integer %>% matrix(.,N,Time)
# write.csv(z,file='z.csv',row.names = F)

#### Figure5 (PM2.5 plot) ####
pm2.5_dat <- read.csv(file = 'misemise/pm2.5_dat_result16.csv',header=T)
z <- read.csv('misemise/z.csv',header=T)
colnames(z) <- paste0('z',1:43)
pm2.5_dat3 <- cbind(pm2.5_dat[,1:(43+4)],z)

library(ggplot2)
library(dplyr)
library(spTimer)
library(dplyr)
center <- read.csv('center.csv',header=T)[,-1]
center_arr <- arrange(center,CTPRV_NM,SGNG_NM)
pm2.5_dat_arr <- arrange(pm2.5_dat,CTPRV_NM,SGNG_NM)

pm2.5_dat_arr <- pm2.5_dat_arr[c(1:19,rep(20,3),21:nrow(pm2.5_dat_arr)),] ## 162 -> 250 인접지역으로 채우기
rownames(pm2.5_dat_arr) <- NULL
pm2.5_dat_arr <- pm2.5_dat_arr[c(1:31,rep(32,3),rep(33,4),34:nrow(pm2.5_dat_arr)),]
rownames(pm2.5_dat_arr) <- NULL
pm2.5_dat_arr <- pm2.5_dat_arr[c(1:39,rep(40,2),41,rep(42,2),43:nrow(pm2.5_dat_arr)),]
rownames(pm2.5_dat_arr) <- NULL
pm2.5_dat_arr <- pm2.5_dat_arr[c(1:49,rep(50,3),51:nrow(pm2.5_dat_arr)),]
rownames(pm2.5_dat_arr) <- NULL
pm2.5_dat_arr <- pm2.5_dat_arr[c(1:72,rep(73,5),74:nrow(pm2.5_dat_arr)),]
rownames(pm2.5_dat_arr) <- NULL
pm2.5_dat_arr <- pm2.5_dat_arr[c(1:104,rep(105,2),106:nrow(pm2.5_dat_arr)),]
rownames(pm2.5_dat_arr) <- NULL
pm2.5_dat_arr <- pm2.5_dat_arr[c(1:106,rep(107,5),rep(108,8),109:nrow(pm2.5_dat_arr)),]
rownames(pm2.5_dat_arr) <- NULL
pm2.5_dat_arr <- pm2.5_dat_arr[c(1:119,rep(120,5),rep(121,16),122:nrow(pm2.5_dat_arr)),]
rownames(pm2.5_dat_arr) <- NULL
pm2.5_dat_arr <- pm2.5_dat_arr[c(1:140,rep(141,25),142,rep(143,5),144:nrow(pm2.5_dat_arr)),]
rownames(pm2.5_dat_arr) <- NULL
pm2.5_dat_arr <- pm2.5_dat_arr[c(1:171,rep(172,10),173:nrow(pm2.5_dat_arr)),]
rownames(pm2.5_dat_arr) <- NULL
pm2.5_dat_arr <- pm2.5_dat_arr[c(1:214,rep(215,2),216:nrow(pm2.5_dat_arr)),]
rownames(pm2.5_dat_arr) <- NULL
pm2.5_dat_arr <- pm2.5_dat_arr[c(1:231,rep(232,2),233:nrow(pm2.5_dat_arr)),]
rownames(pm2.5_dat_arr) <- NULL
pm2.5_dat_arr <- pm2.5_dat_arr[c(1:245,rep(246,4),247),]
rownames(pm2.5_dat_arr) <- NULL

center_dat <- cbind(center_arr,pm2.5_dat_arr[,c(1:2,5:(2*Time+4))])
colnames(center_dat)

colnames(center_dat)[2:3] <- c('X','Y')
center_dat <- center_dat %>% arrange(.,id)

library(stringr)
library(progress)
library(dplyr)
library(ggplot2)
library(rgdal)
map <- readOGR('SIG_20221119/sig.shp',encoding = 'euc-kr')
new_map <- fortify(map,region = 'SIG_CD')
center_dat_list <- list()
center_dat %>% dim
center_dat %>% colnames
center_dat_mean <- cbind(center_dat[,1:7],rowMeans(center_dat[,8:(Time+7)]))
colnames(center_dat_mean)[8] <- '초미세먼지농도'
center_dat %>% dim
for(i in 1:(2*Time)){
  center_dat_list[[i]] <- center_dat[,c(1:7,7+i)]
  colnames(center_dat_list[[i]])[8] <- '초미세먼지농도'
}
P_merge_list=list()
pb <- progress_bar$new(total=(Time))
for(i in 33:35){
  P_merge_list[[i]] <- merge(new_map,center_dat_list[[i]],by='id')
  colnames(P_merge_list[[i]])[2:3] <- c('X좌표','Y좌표')
}
pb <- progress_bar$new(total=Time)
for(i in (Time+33):(Time+35)){
  P_merge_list[[i]] <- merge(new_map,center_dat_list[[i]],by='id')
  colnames(P_merge_list[[i]])[2:3] <- c('X좌표','Y좌표')
}
i=33; i=34; i=35
ggplot() +
  geom_polygon(data = P_merge_list[[i]],
               aes(x=X좌표, y=Y좌표, group=group, fill = 초미세먼지농도), color='blue')+
  scale_fill_gradient(low = "white", high = "black", space = "Lab",
                      guide = "colourbar",limits =c(5,40))+
  theme_bw() + labs(title = str_c("전국 초미세먼지 지도 (20",(i+3)%/%12+22-Time%/%12,"년 "
                                  ,ceiling((i+4)%%12.01),'월)'))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5))
ggplot() +
  geom_polygon(data = P_merge_list[[i+Time]],
               aes(x=X좌표, y=Y좌표, group=group, fill = 초미세먼지농도), color='blue')+
  scale_fill_gradient(low = "white", high = "black", space = "Lab",
                      guide = "colourbar",limits =c(5,40))+
  theme_bw() + labs(title = str_c("추정된 전국 초미세먼지 지도 (20",(i+3)%/%12+22-Time%/%12,"년 "
                                  ,ceiling((i+4)%%12.01),'월)'))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5))
ggplot() +
  geom_polygon(data = P_merge_list[[i+Time]],
               aes(x=X좌표, y=Y좌표, group=group, fill = 초미세먼지농도), color='blue')+
  scale_fill_gradient(low = "white", high = "black", space = "Lab",
                      guide = "colourbar",limits =c(5,40))+
  theme_bw() + labs(title = str_c("추정된 전국 초미세먼지 지도 (20",(i+3)%/%12+22-Time%/%12,"년 "
                                  ,ceiling((i+4)%%12.01),'월)'))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5))

aaa <- c('A','B','C','D','E','F','G','H','I','J','K','L')
P_merge_list2[[33]]$그룹번호 <- aaa[P_merge_list2[[33]]$그룹번호]
P_merge_list2[[34]]$그룹번호 <- aaa[P_merge_list2[[34]]$그룹번호]
P_merge_list2[[35]]$그룹번호 <- aaa[P_merge_list2[[35]]$그룹번호]
color <- c('red','green','yellow','cyan2',
           'orange','pink'); i =33
color <- c('green','yellow','cyan2',
           'orange','pink'); i=34
color <- c('hotpink2','green','yellow','cyan2',
           'green4','orange'); i=35
ggplot() +
  geom_polygon(data = P_merge_list2[[i]],
               aes(x=X좌표, y=Y좌표, group=group, fill = 그룹번호), color='blue')+
  scale_fill_manual(values=color)+
  # scale_fill_gradient(low = "white", high = "black", space = "Lab",
  #                     guide = "colourbar",limits =c(1,11))+
  theme_bw() + labs(title = str_c("군집화된 전국 초미세먼지 지도 (20",(i+3)%/%12+22-Time%/%12,"년 "
                                  ,ceiling((i+4)%%12.01),'월)'))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5))






#### fitting other models ####
library(dplyr)
library(spTimer)
library(bmstdr)

f2 <- y ~ t
tn <- Time
sn <- 162
valids <- c(1, 5, 10)
validt <- sort(sample(1:Time,20))
vrows <- getvalidrows(sn=sn, tn=tn, valids=valids, validt=validt)
M3 <- Bsptime(package="spTimer",formula=f2, data=pm2.5, 
              coordtype="utm", coords=4:5,
              validrows=vrows, model="GP", report=5,
              mchoice=F, scale.transform = "NONE")

M5 <- Bsptime(package="inla", model="AR", formula=f2, data=pm2.5, 
              coordtype="utm", coords=4:5, scale.transform = "NONE", 
              mchoice=T, validrows=vrows)


M6 <- Bsptime(package="spTDyn", model="GP", formula=f2, data=pm2.5, 
              coordtype="utm", coords=4:5, scale.transform = "NONE", n.report=2)

M7 <- Bsptime(package="spBayes",  formula=f2, data=pm2.5, 
              prior.sigma2=c(2, 25),
              prior.tau2 =c(2, 25),
              prior.sigma.eta =c(2, 0.001),
              coordtype="utm", 
              coords=4:5, scale.transform = "SQRT", 
              mchoice=T,  N=5000,  n.report=200)

modfit3 <- M3$fit
modfit5 <- M5$fit
modfit6 <- M6$fit
modfit7 <- M7$fit

fitall3 <- data.frame(modfit3$fitted)
fitall5 <- data.frame(modfit5$fitted)
fitall6 <- data.frame(modfit6$fitted)
fitall7 <- data.frame(modfit7$fitted)

fitted3_mu <- fitall3[,1] %>% matrix(.,nrow=162,byrow=T)
fitted3_mu <- fitted3_mu*sd_y+mean_y
fitted3_sd <- fitall3[,2] %>% matrix(.,nrow=162,byrow=T)
fitted3_sd <- fitted3_sd*sd_y

fitted5_mu <- modfit5$summary.fitted.values[1:nrow(pm2.5),1] %>% matrix(.,nrow=162,byrow=T)
fitted5_mu <- fitted5_mu*sd_y+mean_y
fitted5_sd <- modfit5$summary.fitted.values[1:nrow(pm2.5),2] %>% matrix(.,nrow=162,byrow=T)
fitted5_sd <- fitted5_sd*sd_y

fitted6_mu <- fitall6[,1] %>% matrix(.,nrow=162,byrow=T)
fitted6_mu <- fitted6_mu*sd_y+mean_y
fitted6_sd <- fitall6[,2] %>% matrix(.,nrow=162,byrow=T)
fitted6_sd <- fitted6_sd*sd_y

fitted7_mu <- fitall7[,1] %>% matrix(.,nrow=162,byrow=T)
fitted7_mu <- fitted7_mu*sd_y+mean_y
fitted7_sd <- fitall7[,2] %>% matrix(.,nrow=162,byrow=T)
fitted7_sd <- fitted7_sd*sd_y

colnames(fitted3_mu) <- paste0('mu',1:43)
colnames(fitted3_sd) <- paste0('sd',1:43)

colnames(fitted5_mu) <- paste0('mu',1:43)
colnames(fitted5_sd) <- paste0('sd',1:43)

colnames(fitted6_mu) <- paste0('mu',1:43)
colnames(fitted6_sd) <- paste0('sd',1:43)

colnames(fitted7_mu) <- paste0('mu',1:43)
colnames(fitted7_sd) <- paste0('sd',1:43)

pm2.5_dat <- read.csv('pm2.5_dat2.csv',header=T,fileEncoding = 'euc-kr')
pm2.5_dat1 <- cbind(pm2.5_dat,fitted3_mu,fitted3_sd)
pm2.5_dat2 <- cbind(pm2.5_dat,fitted5_mu,fitted5_sd)
pm2.5_dat3 <- cbind(pm2.5_dat,fitted6_mu,fitted6_sd)
pm2.5_dat4 <- cbind(pm2.5_dat,fitted7_mu,fitted7_sd)

#### compare other model ####
colnames(pm2.5_dat)
Time=43
pm2.5_dat[,4+Time+2]
library(verification)
num <- 1
crps_my <-crps_bms1<-crps_bms2<-crps_bms3<-crps_bms4<- c()

for(num in 1:(Time-2)){
  crps_my[num] <- verification::crps(pm2.5_dat[,(Time+5-num)],
                                     pm2.5_dat[,c((2*Time+5-num),(3*Time+5-num))])$CRPS
}

for(num in 1:Time){
  crps_bms1[num] <- verification::crps(pm2.5_dat1[,(Time+5-num)],
                                      pm2.5_dat1[,c((2*Time+5-num),(3*Time+5-num))])$CRPS
}
for(num in 1:Time){
  crps_bms2[num] <- verification::crps(pm2.5_dat2[,(Time+5-num)],
                                      pm2.5_dat2[,c((2*Time+5-num),(3*Time+5-num))])$CRPS
}
for(num in 1:Time){
  crps_bms3[num] <- verification::crps(pm2.5_dat3[,(Time+5-num)],
                                       pm2.5_dat3[,c((2*Time+5-num),(3*Time+5-num))])$CRPS
}
for(num in 1:Time){
  crps_bms4[num] <- verification::crps(pm2.5_dat4[,(Time+5-num)],
                                       pm2.5_dat4[,c((2*Time+5-num),(3*Time+5-num))])$CRPS
}
mean(crps_my)
mean(crps_bms1)
mean(crps_bms2)
mean(crps_bms3)
mean(crps_bms4)

### fitted my
spT.validation(unlist(c(pm2.5_dat[,c(5:(5+Time-1))])),
               unlist(c(pm2.5_dat[,Time+c(5:(5+Time-1))])))

### fitted other packages
spT.validation(unlist(c(pm2.5_dat1[,c(5:(5+Time-1))])),
               unlist(c(pm2.5_dat1[,Time+c(5:(5+Time-1))]))) spT.validation(unlist(c(pm2.5_dat2[,c(5:(5+Time-1))])),
               unlist(c(pm2.5_dat2[,Time+c(5:(5+Time-1))])))
spT.validation(unlist(c(pm2.5_dat3[,c(5:(5+Time-1))])),
               unlist(c(pm2.5_dat3[,Time+c(5:(5+Time-1))])))
spT.validation(unlist(c(pm2.5_dat4[,c(5:(5+Time-1))])),
               unlist(c(pm2.5_dat4[,Time+c(5:(5+Time-1))])))