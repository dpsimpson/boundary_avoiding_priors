setwd("~/Documents/boundary_avoiding_priors/")
library(rstan)
library(xtable)
options(mc.cores=parallel::detectCores())

y <- c(28,  8, -3,  7, -1,  1, 18, 12)
sigma <- c(15, 10, 16, 11,  9, 11, 10, 18)
standata <- list(y =y, sigma = sigma,J=length(y))
mod <- stan_model("schools_ln.stan")

means <- c(0,0.5,1,1.5,2,2.5,3,4,5,6,7,8,9,10)
standata$tau_sigma <- 1.0
curve(dlnorm(x,10,1),from = 0, to=5,log="y",ylim = c(1e-36,1))
for (mean in means){
  curve(dlnorm(x,mean, 1), add=TRUE,log="y")
}

output= list()
for (mean in means){
  standata$tau_mu <- mean
  

  fit <- sampling(mod, data = standata, control= list(adapt_delta=0.8),refresh=0)
  divs = sum(rstan:::sampler_param_vector(fit,"divergent__"))
  energies_by_chain <- rstan:::sampler_param_matrix(fit, "energy__")
  EBFMIs <- apply(energies_by_chain, 2, function(x) {
    numer <- sum(diff(x)^2)/length(x)
    denom <- var(x)
    numer/denom
  })
  low_BFMI <- any(EBFMIs < 0.2)
  
  
  tmp = list(mu = mean, prior_median= exp(mean), one_pc = qlnorm(p=0.01,meanlog = mean, sdlog=1),
            divergences = divs, low_BFMI=low_BFMI, tau_low =summary(fit)$summary["tau",4], 
            tau_med =summary(fit)$summary["tau",6],
            tau_up = summary(fit)$summary["tau",8]
 )
  output = rbind(output,tmp)
  

}
rownames(output) <- c()

xx= xtable(output,display = c("d","f","f","f","d","s","f","f","f"))
xx
print(xx,type = "html",include.rownames = FALSE)

  mod2 = stan_model("eight_schools_ncp.stan")
  y <- c(28,  8, -3,  7, -1,  1, 18, 12)
  sigma <- c(15, 10, 16, 11,  9, 11, 10, 18)
  standata2 <- list(y =y, sigma = sigma,J=length(y))
  
  fit2=sampling(mod2, data = standata2, control= list(adapt_delta=0.9999999),refresh=0)
  
  ## Inverse gamma
  
  y <- c(28,  8, -3,  7, -1,  1, 18, 12)
  sigma <- c(15, 10, 16, 11,  9, 11, 10, 18)
  standata <- list(y =y, sigma = sigma,J=length(y))
  mod <- stan_model("schools_ig.stan")
  
  standata$tau_shape <- 2
  rates <- rev(c(100,50,40,30,20,10,5,2,1,0.5))
  output= list()
  for (rate in rates){
    standata$tau_rate <- rate
    
    fit <- sampling(mod, data = standata, control= list(adapt_delta=0.8),refresh=0)
    divs = sum(rstan:::sampler_param_vector(fit,"divergent__"))
    energies_by_chain <- rstan:::sampler_param_matrix(fit, "energy__")
    EBFMIs <- apply(energies_by_chain, 2, function(x) {
      numer <- sum(diff(x)^2)/length(x)
      denom <- var(x)
      numer/denom
    })
    low_BFMI <- any(EBFMIs < 0.2)
    
    
    tmp = list(beta = rate, prior_mean = rate, one_pc = 1/qgamma(p=0.991,shape=2,rate=rate),
               divergences = divs, low_BFMI=low_BFMI, tau_low =summary(fit)$summary["tau",4], 
               tau_med =summary(fit)$summary["tau",6],
               tau_up = summary(fit)$summary["tau",8]
    )
    output = rbind(output,tmp)
    
    
  }
  rownames(output) <- c()
  
  xx= xtable(output,display = c("d","f","f","f","d","s","f","f","f"))
  xx
  print(xx,type = "html",include.rownames = FALSE)
