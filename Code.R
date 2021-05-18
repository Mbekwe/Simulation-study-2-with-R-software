####################     Simulation to evaluate the dependence on prevalence of the selected measures     #################### 

####################     Simulation of correlated binary outcome using the method of Lunn and Davies:     #################### 
####################                   A note on generating correlated binary variables                   #################### 

library(lme4)
library(ICC)
library(foreach)
library(doParallel)
library(openxlsx)
library(ggplot2)
library(cowplot)
library(minerva)


##########     est_anova: a function to estimate the intraclass correlation coefficient (ICC)     ##########
##########                        using the analysis of variance estimator                        ##########

###   Arguments   ###
# g: a group variable
# out: the outcome
# data: a dataframe containing g and out

###   Value   ###
# Rho: the estimated intraclass correlation coefficient

est_anova <- function(g, out, data)
{
  # Sample size
  N <- nrow(data)
  
  l <- list(g = substitute(g), out = substitute(out))
  d <- data.frame(g = eval(expr = l$g, envir = data), out=eval(expr = l$out, envir = data))
  
  # Number of clusters
  k <- length(levels(d$g))
  
  # Clusters size
  size <- rowSums(table(d$g,d$out))
  
  # Number of successes in clusters
  zi <- table(d$g,d$out)[, 2]
  
  nA <- (1/(k-1))*(N - sum(size^2)/N)
  msb <- (1/(k-1))*((sum((zi^2)/size)) - ((sum(zi))^2)/N)
  msw <- (1/(N-k))*(sum(zi) - (sum((zi^2)/size)))
  res <- (msb - msw)/(msb + (nA - 1)*msw)
  
  return(list(Rho = res))
}


##########     make_pairs: a function to constitute pairs     ##########
###   Arguments   ##
# dat : a dataframe of 2 columns

###   Value   ###
# interm_res: Pairs constituted

make_pair <- function(x)
{
  n <- length(x)
  
  res <- matrix(nrow=n*(n-1), ncol=2)
  res[,1] <- rep(x, each = n-1)
  
  interm_res <- c()
  for (i in 1:n) interm_res <- c(interm_res,x[-i])
  res[,2] <- interm_res

  res
}

make_pairs <- function(dat)
{
  clust_siz <- table(dat[,1])
  
  interm_res <- c()
  for (i in levels(dat[,1]))
  {
    interm_res <- rbind(interm_res , make_pair(x = dat[dat[,1]==i, 2]))
  }
  interm_res <- data.frame(x1=interm_res[,1] , x2 =interm_res[,2])
  interm_res
}


##########    approx_TCC.1: a function to approximate the tetrachoric correlation coefficient (TCC)     ##########

###   Arguments   ###
# mat: a 2*2 table
# delta: the stop's criterion

###   Value   ###
# rtet: the TCC

approx_TCC.1 <- function(mat, delta)
{
  n  <- sum(mat)
  a <- mat[1,1]; b <- mat[1,2]; c <- mat[2,1]; d <- mat[2,2]
  
  if (b==0 & c==0) r_final <- 1
  
  else if (a==0 & d==0) r_final <- -1
  
  else {
    
    if (qnorm((((a+c)-(b+d))/n + 1 )/2)==0 & qnorm((((a+b)-(c+d))/n + 1 )/2)==0) r_final <- cos(pi*b/(a+b))
    
    else
    {
      if(a==0) {a <- 0.5}
      else if(b==0) {b <- 0.5}
      else if(c==0) {c <- 0.5}
      else if(d==0) {d <- 0.5}
      
      n  <- a + b + c + d
      h <- qnorm((((a+c)-(b+d))/n + 1 )/2)
      k <- qnorm((((a+b)-(c+d))/n + 1 )/2) 
      H <- dnorm(h)
      K <- dnorm(k)
      eps <- (a*d - b*c)/(n^2*H*K)
      
      # First terms
      u <- c(); u[1] <- 1; u[2] <- h*k; u[3] <- (h^2-1)*(k^2-1); u[4] <- h*k*(h^2-3)*(k^2-3)
      
      for (i in 5:150) 
      {
        n <- i - 2
        u[i] <- n*(2*n - 1 - h^2 - k^2)*u[i-2] - n*(n - 1)*(n - 2)^2*u[i-4] + h*k*(u[i-1] + n*(n - 1)*u[i-3])
      } 
      
      coefs <- c(); r <- c(); s <- c()
      
      coefs[1] <- u[1]
      roots <- polyroot(c(-eps,coefs))
      
      r_int  <- Re(roots[which(Mod(roots)==abs(Re(roots)))])
      
      #dr: relative difference
      diff <- 1
      j <- 0
      
      while (diff > delta & j < 150)
      {
        r_final <- r_int
        
        j <- j+1
        
        coefs[j+1] <- u[j+1]/factorial(j+1)
        
        roots <- polyroot(c(-eps,coefs))
        
        real_roots <- Re(roots[which(Mod(roots)==abs(Re(roots)))])
        
        if (length(real_roots)!=0) {r_int <- sign(real_roots)[(which.min(abs(real_roots)))]*min(abs(real_roots))} 
        
        else {r_int <- r_final}
        
        diff <- abs(r_int - r_final)
      }
    }
  }
  return(list(rtet = r_final))
}


##########     sum_f_int: a function to calculate the binary ICC knowing the continuous ICC     ##########

f_int <- function(x, h)
{
  (1/sqrt(1-x^2))*exp(-h^2/(1+x))
}

###   Arguments   ###
# r_tet: the continuous intraclass correlation coefficient
# n: the number of intervals
# h: the threshold of dichotomization
# p: the prevalence

###   Value   ###
# res: the estimated binary ICC

sum_f_int <- function(r_tet, n, h, p)
{ 
  x0 <- 0; s <- 0
  delta <- (r_tet-0)/n
  
  for (i in 1:(n-1))
  {
    x0 <- x0 + delta
    s <- s + f_int(x=x0, h)
  }
  
  s <- delta*((f_int(x=0, h) + f_int(x=r_tet, h))/2 + s)
  
  res <- (1/(2*pi*p*(1-p)))*s
  
  return(res)
}


##########     sum_f_int: a function to calculate the continuous ICC knowing the binary ICC     ##########

###   Arguments   ###
# rhobin: the binary intraclass correlation coefficient
# n: the number of intervals
# h: the threshold of dichotomization
# p: the prevalence

###   Value   ###
# res: the estimated continuous ICC

rec_sum_f_int <- function(rhobin, n, h, p)
{
  f <- function(x, n, h, p) sum_f_int(x, n, h, p)  - rhobin
  res <- uniroot(f = f, interval = c(0,1), n = n, h = h, p = p)$root
  res
}


##########     Estimation of the selected measures     ##########

###   Arguments   ###
# out: the outcome
# g: the group variable
# dat: data
# nsim: the number of simulations

###  Values   ###
# vpc1, vpc2, vpc3, vpc4: the variance partition coefficients
# mor: the median odds ratio

f.int <- function(out, g, dat, nsim)
{
  l <- list(g = substitute(g), out = substitute(out))
  d <- data.frame(g = eval(expr = l$g, envir = dat), out=eval(expr = l$out, envir = dat))
  
  # Fit models
  # m1 <- glmer(out ~ 1 + 1|g, data = d, family = binomial(link = "logit"), nAGQ = 20)
  
  m1 <- glmer(out ~ 1 + 1|g, data = d, family = binomial(link = "logit"), nAGQ = 20,
              control = glmerControl(check.conv.grad = .makeCC("stop", tol = 2e-3, relTol = NULL)))
  
  ## VPC: method A
  # estimate of intercept
  beta0 <- as.numeric(fixef(m1))
  
  # estimate of the between group variance
  sigma_carre_u0 <- as.numeric(VarCorr(m1)[[1]])
  
  # estimate of the prevalence
  p1i_ij <- prop.table(table(d$out))[2]
  
  # estimate of the level one variance
  v11 <- p1i_ij*(1 - p1i_ij)
  
  # estimate of the level two variance
  v12 <- (sigma_carre_u0 *p1i_ij^2)/(1 + exp(beta0))^2
  
  # estimate of the variance partition coefficient
  vpc1 <- v12/(v12 + v11)
  
  
  ## VPC: method B
  # simulate 
  u0j <- rnorm(n = nsim, mean = 0, sd = sqrt(as.numeric(VarCorr(m1)[[1]])))
  
  # estimate of intercept
  beta0 <- as.numeric(fixef(m1))
  
  # estimate of the prevalence
  p2i_ij <- exp(beta0 + u0j)/(1 + exp(beta0 + u0j))
  
  # estimate of the level one variance
  v21 <- mean(p2i_ij*(1-p2i_ij))
  
  # estimate of the level two variance
  v22 <- var(p2i_ij)
  
  # estimate of the variance partition coefficient
  vpc2 <- v22/(v22 + v21)
  
  
  ## VPC: method C
  
  # estimate of the variance partition coefficient
  vpc3 <- ICCest(x = g, y = as.numeric(out), data = d)$ICC
  
  if(vpc3 < 0) {vpc3 <- 0}
  
  ## VPC: method D
  # the level one variance is 1
  
  # estimate of the level two variance
  v42 <- as.numeric(VarCorr(m1)[[1]])
  
  # estimate of the variance partition coefficient
  vpc4 <- v42/(v42 + pi^2/3)
  
  
  ## MOR
  # the level two variance is sigma_carre_u0
  
  # estimate of the variance partition coefficient
  mor <- exp(qnorm(0.75)*sqrt(2*sigma_carre_u0))
  
  return(list(vpc1 = vpc1, vpc2 = vpc2, vpc3 = vpc3, vpc4 = vpc4, mor = mor))
}

###   Arguments   ###
# k: the number of clusters
# m: the average size of clusters
# v: the variance of cluster sizes
# rho: the intraclass correlation coefficient
# p: the prevalence

###   Values   ###
# p_est
# vpc1_est, vpc2_est, vpc3_est, vpc4_est: the estimated variance partition coefficients
# mor_est: the estimated median odds ratio
# tcc_est: the estimated TCC using original formula
# neg_tcc_est: an output egal to 1 if the estimated TCC is less than 0 and egal to 0 otherwise
# rho_est: the estimated ICC
# neg_rho_est: an output egal to 1 if the estimated ICC is less than 0 and egal to 0 otherwise
# tcc_kirk_est: the estimated TCC using Kirk's formula

sim.1 <- function(k, m, v, rho, p)
{
  # variable cluster sizes
  clust_siz <- rnbinom(n = k, size = m^2/(v-m), mu = m)
  
  # delete empty clusters
  i <- 0
  while (sum(is.element(el = c(0,1), set = clust_siz)) >= 1 & i < 1000)
  {
    clust_siz <- rnbinom(n = k, size = m^2/(v-m), mu = m)
    i <- i + 1
  }
  
  # sample size
  n <- sum(clust_siz)
  
  
  # Simulation of correlated binary outcome using the method of Lunn and Davies
  # model: Xij = (1-Uij)Yij + UijZi
  # Zi ~ Binom(1,pi); Yij ~ Binom(1,pi); Uij ~ Binom(1,sqrt(rho))
  
  Z <- rep(rbinom(k, 1, p), clust_siz)
  Y <- rbinom(n, 1, p)
  U <- rbinom(n, 1, sqrt(rho))
  
  X <- (1-U)*Y + U*Z
  
  while(length(unique(X)) == 1)
  {
    Z <- rep(rbinom(k, 1, p), clust_siz)
    Y <- rbinom(n, 1, p)
    U <- rbinom(n, 1, sqrt(rho))
    
    X <- (1-U)*Y + U*Z
  }
  
  # estimated prevalence
  p_est <- as.numeric(prop.table(table(X))[2])
  
  # variable "cluster"
  clust <- factor(rep(1:k, clust_siz))
  
  data_sim <- data.frame(Outcome = X, Cluster = clust)
  
  res <- f.int(out = Outcome, g = Cluster, dat = data_sim, nsim = 5000)
  
  vpc1_est <- res$vpc1
  vpc2_est <- res$vpc2
  vpc3_est <- res$vpc3
  vpc4_est <- res$vpc4
  
  mor_est <- res$mor
  
  
  # estimation of the tetrachoric correlation coefficient using binary data
  tcc_est <- approx_TCC.1(mat = table(make_pairs(data_sim[ , c("Cluster","Outcome")])), delta = 0.00001)$rtet
  
  neg_tcc_est <- ifelse(tcc_est < 0, 1, 0)
  if(tcc_est < 0) {tcc_est <- 0}
  
  # estimation of rho using the Anova estimator
  rho_est <- est_anova(g = Cluster, out = Outcome, data = data_sim)$Rho
  
  neg_rho_est <- ifelse(rho_est < 0, 1, 0)
  if(rho_est < 0) {rho_est <- 0}
  
  # estimation of the tetrachoric correlation coefficient using ICC anova estimation
  tcc_kirk_est <- rec_sum_f_int(rhobin = rho_est, n = 10, h = qnorm(p), p = p)
  
  return(list(p_est = p_est, vpc1_est = vpc1_est, vpc2_est = vpc2_est, vpc3_est = vpc3_est, vpc4_est = vpc4_est, mor_est = mor_est,
              tcc_est = tcc_est, neg_tcc_est = neg_tcc_est, rho_est = rho_est, neg_rho_est = neg_rho_est, tcc_kirk_est = tcc_kirk_est))
  
}

##########                                            Data generation                                            ##########

nsimul <- 10000
ncores <- 8

##########################################     With rho = 0.01     ##########################################


prevalence <- seq(from = 0.01, to = 0.99, length.out = 20)
s <- length(prevalence)

Bin_ICC_001 <- c()
for (i in 1:s) Bin_ICC_001[i] <- sum_f_int(r_tet = 0.01, n = 100, h = qnorm(1-prevalence[i]), p = prevalence[i])

res <- c(); RES <- array(data = NA, dim = c(nsimul, 20, 11), dimnames = list(NULL, NULL, c("Mean_p", "Mean_vpc1", "Mean_vpc2", "Mean_vpc3", 
                                                                                           "Mean_vpc4", "Mean_mor", "Mean_tcc", "mean_neg_tcc_est", "Mean_rho_est", "mean_neg_rho_est",
                                                                                           "Mean_tcc_kirk_est")))

for (j in 1:s)
{
  registerDoParallel(cores = ncores)
  res <- foreach(1:nsimul, .packages = c("lme4", "ICC")) %dopar% sim.1(k = 10, m = 25, v = 225, rho = Bin_ICC_001[j], p = prevalence[j])
  
  RES[, j , 1] <- sapply(res, function(a){a$p_est})
  RES[, j , 2] <- sapply(res, function(a){a$vpc1_est})
  RES[, j , 3] <- sapply(res, function(a){a$vpc2_est})
  RES[, j , 4] <- sapply(res, function(a){a$vpc3_est})
  RES[, j , 5] <- sapply(res, function(a){a$vpc4_est})
  RES[, j , 6] <- sapply(res, function(a){a$mor_est})
  RES[, j , 7] <- sapply(res, function(a){a$rho_est})
  RES[, j , 8] <- sapply(res, function(a){a$tcc_est})
  RES[, j , 9] <- sapply(res, function(a){a$tcc_kirk_est})
  RES[, j , 10] <- sapply(res, function(a){a$neg_tcc_est})
  RES[, j , 11] <- sapply(res, function(a){a$neg_rho_est})
  
  stopImplicitCluster()
}

wb <- createWorkbook()
addWorksheet(wb, "Mean_p")
addWorksheet(wb, "Mean_vpc1")
addWorksheet(wb, "Mean_vpc2")
addWorksheet(wb, "Mean_vpc3")
addWorksheet(wb, "Mean_vpc4")
addWorksheet(wb, "Mean_mor")
addWorksheet(wb, "Mean_tcc")
addWorksheet(wb, "mean_neg_tcc_est")
addWorksheet(wb, "Mean_rho_est")
addWorksheet(wb, "mean_neg_rho_est")
addWorksheet(wb, "Mean_tcc_kirk_est")

writeDataTable(wb, "Mean_p", x = as.data.frame(RES[, , 1]))
writeDataTable(wb, "Mean_vpc1", x = as.data.frame(RES[, , 2]))
writeDataTable(wb, "Mean_vpc2", x = as.data.frame(RES[, , 3]))
writeDataTable(wb, "Mean_vpc3", x = as.data.frame(RES[, , 4]))
writeDataTable(wb, "Mean_vpc4", x = as.data.frame(RES[, , 5]))
writeDataTable(wb, "Mean_mor", x = as.data.frame(RES[, , 6]))
writeDataTable(wb, "Mean_tcc", x = as.data.frame(RES[, , 7]))
writeDataTable(wb, "mean_neg_tcc_est", x = as.data.frame(RES[, , 8]))
writeDataTable(wb, "Mean_rho_est", x = as.data.frame(RES[, , 9]))
writeDataTable(wb, "mean_neg_rho_est", x = as.data.frame(RES[, , 10]))
writeDataTable(wb, "Mean_tcc_kirk_est", x = as.data.frame(RES[, , 11]))

saveWorkbook(wb, file = "10000and10and25and225and001.xlsx")

##########################################

res <- c(); RES <- array(data = NA, dim = c(nsimul, 20, 11), dimnames = list(NULL, NULL, c("Mean_p", "Mean_vpc1", "Mean_vpc2", "Mean_vpc3", 
                                                                                           "Mean_vpc4", "Mean_mor", "Mean_tcc", "mean_neg_tcc_est", "Mean_rho_est", "mean_neg_rho_est",
                                                                                           "Mean_tcc_kirk_est")))

for (j in 1:s)
{
  registerDoParallel(cores = ncores)
  res <- foreach(1:nsimul, .packages = c("lme4", "ICC")) %dopar% sim.1(k = 20, m = 25, v = 225, rho = Bin_ICC_001[j], p = prevalence[j])
  
  RES[, j , 1] <- sapply(res, function(a){a$p_est})
  RES[, j , 2] <- sapply(res, function(a){a$vpc1_est})
  RES[, j , 3] <- sapply(res, function(a){a$vpc2_est})
  RES[, j , 4] <- sapply(res, function(a){a$vpc3_est})
  RES[, j , 5] <- sapply(res, function(a){a$vpc4_est})
  RES[, j , 6] <- sapply(res, function(a){a$mor_est})
  RES[, j , 7] <- sapply(res, function(a){a$rho_est})
  RES[, j , 8] <- sapply(res, function(a){a$tcc_est})
  RES[, j , 9] <- sapply(res, function(a){a$tcc_kirk_est})
  RES[, j , 10] <- sapply(res, function(a){a$neg_tcc_est})
  RES[, j , 11] <- sapply(res, function(a){a$neg_rho_est})
  
  stopImplicitCluster()
}

wb <- createWorkbook()
addWorksheet(wb, "Mean_p")
addWorksheet(wb, "Mean_vpc1")
addWorksheet(wb, "Mean_vpc2")
addWorksheet(wb, "Mean_vpc3")
addWorksheet(wb, "Mean_vpc4")
addWorksheet(wb, "Mean_mor")
addWorksheet(wb, "Mean_tcc")
addWorksheet(wb, "mean_neg_tcc_est")
addWorksheet(wb, "Mean_rho_est")
addWorksheet(wb, "mean_neg_rho_est")
addWorksheet(wb, "Mean_tcc_kirk_est")

writeDataTable(wb, "Mean_p", x = as.data.frame(RES[, , 1]))
writeDataTable(wb, "Mean_vpc1", x = as.data.frame(RES[, , 2]))
writeDataTable(wb, "Mean_vpc2", x = as.data.frame(RES[, , 3]))
writeDataTable(wb, "Mean_vpc3", x = as.data.frame(RES[, , 4]))
writeDataTable(wb, "Mean_vpc4", x = as.data.frame(RES[, , 5]))
writeDataTable(wb, "Mean_mor", x = as.data.frame(RES[, , 6]))
writeDataTable(wb, "Mean_tcc", x = as.data.frame(RES[, , 7]))
writeDataTable(wb, "mean_neg_tcc_est", x = as.data.frame(RES[, , 8]))
writeDataTable(wb, "Mean_rho_est", x = as.data.frame(RES[, , 9]))
writeDataTable(wb, "mean_neg_rho_est", x = as.data.frame(RES[, , 10]))
writeDataTable(wb, "Mean_tcc_kirk_est", x = as.data.frame(RES[, , 11]))

saveWorkbook(wb, file = "10000and20and25and225and001.xlsx")

##########################################

res <- c(); RES <- array(data = NA, dim = c(nsimul, 20, 11), dimnames = list(NULL, NULL, c("Mean_p", "Mean_vpc1", "Mean_vpc2", "Mean_vpc3", 
                                                                                           "Mean_vpc4", "Mean_mor", "Mean_tcc", "mean_neg_tcc_est", "Mean_rho_est", "mean_neg_rho_est",
                                                                                           "Mean_tcc_kirk_est")))

for (j in 1:s)
{
  registerDoParallel(cores = ncores)
  res <- foreach(1:nsimul, .packages = c("lme4", "ICC")) %dopar% sim.1(k = 50, m = 25, v = 225, rho = Bin_ICC_001[j], p = prevalence[j])
  
  RES[, j , 1] <- sapply(res, function(a){a$p_est})
  RES[, j , 2] <- sapply(res, function(a){a$vpc1_est})
  RES[, j , 3] <- sapply(res, function(a){a$vpc2_est})
  RES[, j , 4] <- sapply(res, function(a){a$vpc3_est})
  RES[, j , 5] <- sapply(res, function(a){a$vpc4_est})
  RES[, j , 6] <- sapply(res, function(a){a$mor_est})
  RES[, j , 7] <- sapply(res, function(a){a$rho_est})
  RES[, j , 8] <- sapply(res, function(a){a$tcc_est})
  RES[, j , 9] <- sapply(res, function(a){a$tcc_kirk_est})
  RES[, j , 10] <- sapply(res, function(a){a$neg_tcc_est})
  RES[, j , 11] <- sapply(res, function(a){a$neg_rho_est})
  
  stopImplicitCluster()
}

wb <- createWorkbook()
addWorksheet(wb, "Mean_p")
addWorksheet(wb, "Mean_vpc1")
addWorksheet(wb, "Mean_vpc2")
addWorksheet(wb, "Mean_vpc3")
addWorksheet(wb, "Mean_vpc4")
addWorksheet(wb, "Mean_mor")
addWorksheet(wb, "Mean_tcc")
addWorksheet(wb, "mean_neg_tcc_est")
addWorksheet(wb, "Mean_rho_est")
addWorksheet(wb, "mean_neg_rho_est")
addWorksheet(wb, "Mean_tcc_kirk_est")

writeDataTable(wb, "Mean_p", x = as.data.frame(RES[, , 1]))
writeDataTable(wb, "Mean_vpc1", x = as.data.frame(RES[, , 2]))
writeDataTable(wb, "Mean_vpc2", x = as.data.frame(RES[, , 3]))
writeDataTable(wb, "Mean_vpc3", x = as.data.frame(RES[, , 4]))
writeDataTable(wb, "Mean_vpc4", x = as.data.frame(RES[, , 5]))
writeDataTable(wb, "Mean_mor", x = as.data.frame(RES[, , 6]))
writeDataTable(wb, "Mean_tcc", x = as.data.frame(RES[, , 7]))
writeDataTable(wb, "mean_neg_tcc_est", x = as.data.frame(RES[, , 8]))
writeDataTable(wb, "Mean_rho_est", x = as.data.frame(RES[, , 9]))
writeDataTable(wb, "mean_neg_rho_est", x = as.data.frame(RES[, , 10]))
writeDataTable(wb, "Mean_tcc_kirk_est", x = as.data.frame(RES[, , 11]))

saveWorkbook(wb, file = "10000and50and25and225and001.xlsx")


##########################################     With rho = 0.05     ##########################################

Bin_ICC_005 <- c()
for (i in 1:s) Bin_ICC_005[i] <- sum_f_int(r_tet = 0.05, n = 100, h = qnorm(1-prevalence[i]), p = prevalence[i])

res <- c(); RES <- array(data = NA, dim = c(nsimul, 20, 11), dimnames = list(NULL, NULL, c("Mean_p", "Mean_vpc1", "Mean_vpc2", "Mean_vpc3", 
                                                                                           "Mean_vpc4", "Mean_mor", "Mean_tcc", "mean_neg_tcc_est", "Mean_rho_est", "mean_neg_rho_est",
                                                                                           "Mean_tcc_kirk_est")))

for (j in 1:s)
{
  registerDoParallel(cores = ncores)
  res <- foreach(1:nsimul, .packages = c("lme4", "ICC")) %dopar% sim.1(k = 10, m = 25, v = 225, rho = Bin_ICC_005[j], p = prevalence[j])
  
  RES[, j , 1] <- sapply(res, function(a){a$p_est})
  RES[, j , 2] <- sapply(res, function(a){a$vpc1_est})
  RES[, j , 3] <- sapply(res, function(a){a$vpc2_est})
  RES[, j , 4] <- sapply(res, function(a){a$vpc3_est})
  RES[, j , 5] <- sapply(res, function(a){a$vpc4_est})
  RES[, j , 6] <- sapply(res, function(a){a$mor_est})
  RES[, j , 7] <- sapply(res, function(a){a$rho_est})
  RES[, j , 8] <- sapply(res, function(a){a$tcc_est})
  RES[, j , 9] <- sapply(res, function(a){a$tcc_kirk_est})
  RES[, j , 10] <- sapply(res, function(a){a$neg_tcc_est})
  RES[, j , 11] <- sapply(res, function(a){a$neg_rho_est})
  
  stopImplicitCluster()
}

wb <- createWorkbook()
addWorksheet(wb, "Mean_p")
addWorksheet(wb, "Mean_vpc1")
addWorksheet(wb, "Mean_vpc2")
addWorksheet(wb, "Mean_vpc3")
addWorksheet(wb, "Mean_vpc4")
addWorksheet(wb, "Mean_mor")
addWorksheet(wb, "Mean_tcc")
addWorksheet(wb, "mean_neg_tcc_est")
addWorksheet(wb, "Mean_rho_est")
addWorksheet(wb, "mean_neg_rho_est")
addWorksheet(wb, "Mean_tcc_kirk_est")

writeDataTable(wb, "Mean_p", x = as.data.frame(RES[, , 1]))
writeDataTable(wb, "Mean_vpc1", x = as.data.frame(RES[, , 2]))
writeDataTable(wb, "Mean_vpc2", x = as.data.frame(RES[, , 3]))
writeDataTable(wb, "Mean_vpc3", x = as.data.frame(RES[, , 4]))
writeDataTable(wb, "Mean_vpc4", x = as.data.frame(RES[, , 5]))
writeDataTable(wb, "Mean_mor", x = as.data.frame(RES[, , 6]))
writeDataTable(wb, "Mean_tcc", x = as.data.frame(RES[, , 7]))
writeDataTable(wb, "mean_neg_tcc_est", x = as.data.frame(RES[, , 8]))
writeDataTable(wb, "Mean_rho_est", x = as.data.frame(RES[, , 9]))
writeDataTable(wb, "mean_neg_rho_est", x = as.data.frame(RES[, , 10]))
writeDataTable(wb, "Mean_tcc_kirk_est", x = as.data.frame(RES[, , 11]))

saveWorkbook(wb, file = "10000and10and25and225and005.xlsx")

##########################################

res <- c(); RES <- array(data = NA, dim = c(nsimul, 20, 11), dimnames = list(NULL, NULL, c("Mean_p", "Mean_vpc1", "Mean_vpc2", "Mean_vpc3", 
                                                                                           "Mean_vpc4", "Mean_mor", "Mean_tcc", "mean_neg_tcc_est", "Mean_rho_est", "mean_neg_rho_est",
                                                                                           "Mean_tcc_kirk_est")))

for (j in 1:s)
{
  registerDoParallel(cores = ncores)
  res <- foreach(1:nsimul, .packages = c("lme4", "ICC")) %dopar% sim.1(k = 20, m = 25, v = 225, rho = Bin_ICC_005[j], p = prevalence[j])
  
  RES[, j , 1] <- sapply(res, function(a){a$p_est})
  RES[, j , 2] <- sapply(res, function(a){a$vpc1_est})
  RES[, j , 3] <- sapply(res, function(a){a$vpc2_est})
  RES[, j , 4] <- sapply(res, function(a){a$vpc3_est})
  RES[, j , 5] <- sapply(res, function(a){a$vpc4_est})
  RES[, j , 6] <- sapply(res, function(a){a$mor_est})
  RES[, j , 7] <- sapply(res, function(a){a$rho_est})
  RES[, j , 8] <- sapply(res, function(a){a$tcc_est})
  RES[, j , 9] <- sapply(res, function(a){a$tcc_kirk_est})
  RES[, j , 10] <- sapply(res, function(a){a$neg_tcc_est})
  RES[, j , 11] <- sapply(res, function(a){a$neg_rho_est})
  
  stopImplicitCluster()
}

wb <- createWorkbook()
addWorksheet(wb, "Mean_p")
addWorksheet(wb, "Mean_vpc1")
addWorksheet(wb, "Mean_vpc2")
addWorksheet(wb, "Mean_vpc3")
addWorksheet(wb, "Mean_vpc4")
addWorksheet(wb, "Mean_mor")
addWorksheet(wb, "Mean_tcc")
addWorksheet(wb, "mean_neg_tcc_est")
addWorksheet(wb, "Mean_rho_est")
addWorksheet(wb, "mean_neg_rho_est")
addWorksheet(wb, "Mean_tcc_kirk_est")

writeDataTable(wb, "Mean_p", x = as.data.frame(RES[, , 1]))
writeDataTable(wb, "Mean_vpc1", x = as.data.frame(RES[, , 2]))
writeDataTable(wb, "Mean_vpc2", x = as.data.frame(RES[, , 3]))
writeDataTable(wb, "Mean_vpc3", x = as.data.frame(RES[, , 4]))
writeDataTable(wb, "Mean_vpc4", x = as.data.frame(RES[, , 5]))
writeDataTable(wb, "Mean_mor", x = as.data.frame(RES[, , 6]))
writeDataTable(wb, "Mean_tcc", x = as.data.frame(RES[, , 7]))
writeDataTable(wb, "mean_neg_tcc_est", x = as.data.frame(RES[, , 8]))
writeDataTable(wb, "Mean_rho_est", x = as.data.frame(RES[, , 9]))
writeDataTable(wb, "mean_neg_rho_est", x = as.data.frame(RES[, , 10]))
writeDataTable(wb, "Mean_tcc_kirk_est", x = as.data.frame(RES[, , 11]))

saveWorkbook(wb, file = "10000and20and25and225and005.xlsx")

##########################################

res <- c(); RES <- array(data = NA, dim = c(nsimul, 20, 11), dimnames = list(NULL, NULL, c("Mean_p", "Mean_vpc1", "Mean_vpc2", "Mean_vpc3", 
                                                                                           "Mean_vpc4", "Mean_mor", "Mean_tcc", "mean_neg_tcc_est", "Mean_rho_est", "mean_neg_rho_est",
                                                                                           "Mean_tcc_kirk_est")))

for (j in 1:s)
{
  registerDoParallel(cores = ncores)
  res <- foreach(1:nsimul, .packages = c("lme4", "ICC")) %dopar% sim.1(k = 50, m = 25, v = 225, rho = Bin_ICC_005[j], p = prevalence[j])
  
  RES[, j , 1] <- sapply(res, function(a){a$p_est})
  RES[, j , 2] <- sapply(res, function(a){a$vpc1_est})
  RES[, j , 3] <- sapply(res, function(a){a$vpc2_est})
  RES[, j , 4] <- sapply(res, function(a){a$vpc3_est})
  RES[, j , 5] <- sapply(res, function(a){a$vpc4_est})
  RES[, j , 6] <- sapply(res, function(a){a$mor_est})
  RES[, j , 7] <- sapply(res, function(a){a$rho_est})
  RES[, j , 8] <- sapply(res, function(a){a$tcc_est})
  RES[, j , 9] <- sapply(res, function(a){a$tcc_kirk_est})
  RES[, j , 10] <- sapply(res, function(a){a$neg_tcc_est})
  RES[, j , 11] <- sapply(res, function(a){a$neg_rho_est})
  
  stopImplicitCluster()
}

wb <- createWorkbook()
addWorksheet(wb, "Mean_p")
addWorksheet(wb, "Mean_vpc1")
addWorksheet(wb, "Mean_vpc2")
addWorksheet(wb, "Mean_vpc3")
addWorksheet(wb, "Mean_vpc4")
addWorksheet(wb, "Mean_mor")
addWorksheet(wb, "Mean_tcc")
addWorksheet(wb, "mean_neg_tcc_est")
addWorksheet(wb, "Mean_rho_est")
addWorksheet(wb, "mean_neg_rho_est")
addWorksheet(wb, "Mean_tcc_kirk_est")

writeDataTable(wb, "Mean_p", x = as.data.frame(RES[, , 1]))
writeDataTable(wb, "Mean_vpc1", x = as.data.frame(RES[, , 2]))
writeDataTable(wb, "Mean_vpc2", x = as.data.frame(RES[, , 3]))
writeDataTable(wb, "Mean_vpc3", x = as.data.frame(RES[, , 4]))
writeDataTable(wb, "Mean_vpc4", x = as.data.frame(RES[, , 5]))
writeDataTable(wb, "Mean_mor", x = as.data.frame(RES[, , 6]))
writeDataTable(wb, "Mean_tcc", x = as.data.frame(RES[, , 7]))
writeDataTable(wb, "mean_neg_tcc_est", x = as.data.frame(RES[, , 8]))
writeDataTable(wb, "Mean_rho_est", x = as.data.frame(RES[, , 9]))
writeDataTable(wb, "mean_neg_rho_est", x = as.data.frame(RES[, , 10]))
writeDataTable(wb, "Mean_tcc_kirk_est", x = as.data.frame(RES[, , 11]))

saveWorkbook(wb, file = "10000and50and25and225and005.xlsx")


##########################################     With rho = 0.3     ##########################################

Bin_ICC_03 <- c()
for (i in 1:s) Bin_ICC_03[i] <- sum_f_int(r_tet = 0.3, n = 100, h = qnorm(1-prevalence[i]), p = prevalence[i])

res <- c(); RES <- array(data = NA, dim = c(nsimul, 20, 11), dimnames = list(NULL, NULL, c("Mean_p", "Mean_vpc1", "Mean_vpc2", "Mean_vpc3", 
                                                                                           "Mean_vpc4", "Mean_mor", "Mean_tcc", "mean_neg_tcc_est", "Mean_rho_est", "mean_neg_rho_est",
                                                                                           "Mean_tcc_kirk_est")))

for (j in 1:s)
{
  registerDoParallel(cores = ncores)
  res <- foreach(1:nsimul, .packages = c("lme4", "ICC")) %dopar% sim.1(k = 10, m = 25, v = 225, rho = Bin_ICC_03[j], p = prevalence[j])
  
  RES[, j , 1] <- sapply(res, function(a){a$p_est})
  RES[, j , 2] <- sapply(res, function(a){a$vpc1_est})
  RES[, j , 3] <- sapply(res, function(a){a$vpc2_est})
  RES[, j , 4] <- sapply(res, function(a){a$vpc3_est})
  RES[, j , 5] <- sapply(res, function(a){a$vpc4_est})
  RES[, j , 6] <- sapply(res, function(a){a$mor_est})
  RES[, j , 7] <- sapply(res, function(a){a$rho_est})
  RES[, j , 8] <- sapply(res, function(a){a$tcc_est})
  RES[, j , 9] <- sapply(res, function(a){a$tcc_kirk_est})
  RES[, j , 10] <- sapply(res, function(a){a$neg_tcc_est})
  RES[, j , 11] <- sapply(res, function(a){a$neg_rho_est})
  
  stopImplicitCluster()
}

wb <- createWorkbook()
addWorksheet(wb, "Mean_p")
addWorksheet(wb, "Mean_vpc1")
addWorksheet(wb, "Mean_vpc2")
addWorksheet(wb, "Mean_vpc3")
addWorksheet(wb, "Mean_vpc4")
addWorksheet(wb, "Mean_mor")
addWorksheet(wb, "Mean_tcc")
addWorksheet(wb, "mean_neg_tcc_est")
addWorksheet(wb, "Mean_rho_est")
addWorksheet(wb, "mean_neg_rho_est")
addWorksheet(wb, "Mean_tcc_kirk_est")

writeDataTable(wb, "Mean_p", x = as.data.frame(RES[, , 1]))
writeDataTable(wb, "Mean_vpc1", x = as.data.frame(RES[, , 2]))
writeDataTable(wb, "Mean_vpc2", x = as.data.frame(RES[, , 3]))
writeDataTable(wb, "Mean_vpc3", x = as.data.frame(RES[, , 4]))
writeDataTable(wb, "Mean_vpc4", x = as.data.frame(RES[, , 5]))
writeDataTable(wb, "Mean_mor", x = as.data.frame(RES[, , 6]))
writeDataTable(wb, "Mean_tcc", x = as.data.frame(RES[, , 7]))
writeDataTable(wb, "mean_neg_tcc_est", x = as.data.frame(RES[, , 8]))
writeDataTable(wb, "Mean_rho_est", x = as.data.frame(RES[, , 9]))
writeDataTable(wb, "mean_neg_rho_est", x = as.data.frame(RES[, , 10]))
writeDataTable(wb, "Mean_tcc_kirk_est", x = as.data.frame(RES[, , 11]))

saveWorkbook(wb, file = "10000and10and25and225and03.xlsx")

##########################################

res <- c(); RES <- array(data = NA, dim = c(nsimul, 20, 11), dimnames = list(NULL, NULL, c("Mean_p", "Mean_vpc1", "Mean_vpc2", "Mean_vpc3", 
                                                                                           "Mean_vpc4", "Mean_mor", "Mean_tcc", "mean_neg_tcc_est", "Mean_rho_est", "mean_neg_rho_est",
                                                                                           "Mean_tcc_kirk_est")))

for (j in 1:s)
{
  registerDoParallel(cores = ncores)
  res <- foreach(1:nsimul, .packages = c("lme4", "ICC")) %dopar% sim.1(k = 20, m = 25, v = 225, rho = Bin_ICC_03[j], p = prevalence[j])
  
  RES[, j , 1] <- sapply(res, function(a){a$p_est})
  RES[, j , 2] <- sapply(res, function(a){a$vpc1_est})
  RES[, j , 3] <- sapply(res, function(a){a$vpc2_est})
  RES[, j , 4] <- sapply(res, function(a){a$vpc3_est})
  RES[, j , 5] <- sapply(res, function(a){a$vpc4_est})
  RES[, j , 6] <- sapply(res, function(a){a$mor_est})
  RES[, j , 7] <- sapply(res, function(a){a$rho_est})
  RES[, j , 8] <- sapply(res, function(a){a$tcc_est})
  RES[, j , 9] <- sapply(res, function(a){a$tcc_kirk_est})
  RES[, j , 10] <- sapply(res, function(a){a$neg_tcc_est})
  RES[, j , 11] <- sapply(res, function(a){a$neg_rho_est})
  
  stopImplicitCluster()
}

wb <- createWorkbook()
addWorksheet(wb, "Mean_p")
addWorksheet(wb, "Mean_vpc1")
addWorksheet(wb, "Mean_vpc2")
addWorksheet(wb, "Mean_vpc3")
addWorksheet(wb, "Mean_vpc4")
addWorksheet(wb, "Mean_mor")
addWorksheet(wb, "Mean_tcc")
addWorksheet(wb, "mean_neg_tcc_est")
addWorksheet(wb, "Mean_rho_est")
addWorksheet(wb, "mean_neg_rho_est")
addWorksheet(wb, "Mean_tcc_kirk_est")

writeDataTable(wb, "Mean_p", x = as.data.frame(RES[, , 1]))
writeDataTable(wb, "Mean_vpc1", x = as.data.frame(RES[, , 2]))
writeDataTable(wb, "Mean_vpc2", x = as.data.frame(RES[, , 3]))
writeDataTable(wb, "Mean_vpc3", x = as.data.frame(RES[, , 4]))
writeDataTable(wb, "Mean_vpc4", x = as.data.frame(RES[, , 5]))
writeDataTable(wb, "Mean_mor", x = as.data.frame(RES[, , 6]))
writeDataTable(wb, "Mean_tcc", x = as.data.frame(RES[, , 7]))
writeDataTable(wb, "mean_neg_tcc_est", x = as.data.frame(RES[, , 8]))
writeDataTable(wb, "Mean_rho_est", x = as.data.frame(RES[, , 9]))
writeDataTable(wb, "mean_neg_rho_est", x = as.data.frame(RES[, , 10]))
writeDataTable(wb, "Mean_tcc_kirk_est", x = as.data.frame(RES[, , 11]))

saveWorkbook(wb, file = "10000and20and25and225and03.xlsx")

##########################################

res <- c(); RES <- array(data = NA, dim = c(nsimul, 20, 11), dimnames = list(NULL, NULL, c("Mean_p", "Mean_vpc1", "Mean_vpc2", "Mean_vpc3", 
                                                                                           "Mean_vpc4", "Mean_mor", "Mean_tcc", "mean_neg_tcc_est", "Mean_rho_est", "mean_neg_rho_est",
                                                                                           "Mean_tcc_kirk_est")))

for (j in 1:s)
{
  registerDoParallel(cores = ncores)
  res <- foreach(1:nsimul, .packages = c("lme4", "ICC")) %dopar% sim.1(k = 50, m = 25, v = 225, rho = Bin_ICC_03[j], p = prevalence[j])
  
  RES[, j , 1] <- sapply(res, function(a){a$p_est})
  RES[, j , 2] <- sapply(res, function(a){a$vpc1_est})
  RES[, j , 3] <- sapply(res, function(a){a$vpc2_est})
  RES[, j , 4] <- sapply(res, function(a){a$vpc3_est})
  RES[, j , 5] <- sapply(res, function(a){a$vpc4_est})
  RES[, j , 6] <- sapply(res, function(a){a$mor_est})
  RES[, j , 7] <- sapply(res, function(a){a$rho_est})
  RES[, j , 8] <- sapply(res, function(a){a$tcc_est})
  RES[, j , 9] <- sapply(res, function(a){a$tcc_kirk_est})
  RES[, j , 10] <- sapply(res, function(a){a$neg_tcc_est})
  RES[, j , 11] <- sapply(res, function(a){a$neg_rho_est})
  
  stopImplicitCluster()
}

wb <- createWorkbook()
addWorksheet(wb, "Mean_p")
addWorksheet(wb, "Mean_vpc1")
addWorksheet(wb, "Mean_vpc2")
addWorksheet(wb, "Mean_vpc3")
addWorksheet(wb, "Mean_vpc4")
addWorksheet(wb, "Mean_mor")
addWorksheet(wb, "Mean_tcc")
addWorksheet(wb, "mean_neg_tcc_est")
addWorksheet(wb, "Mean_rho_est")
addWorksheet(wb, "mean_neg_rho_est")
addWorksheet(wb, "Mean_tcc_kirk_est")

writeDataTable(wb, "Mean_p", x = as.data.frame(RES[, , 1]))
writeDataTable(wb, "Mean_vpc1", x = as.data.frame(RES[, , 2]))
writeDataTable(wb, "Mean_vpc2", x = as.data.frame(RES[, , 3]))
writeDataTable(wb, "Mean_vpc3", x = as.data.frame(RES[, , 4]))
writeDataTable(wb, "Mean_vpc4", x = as.data.frame(RES[, , 5]))
writeDataTable(wb, "Mean_mor", x = as.data.frame(RES[, , 6]))
writeDataTable(wb, "Mean_tcc", x = as.data.frame(RES[, , 7]))
writeDataTable(wb, "mean_neg_tcc_est", x = as.data.frame(RES[, , 8]))
writeDataTable(wb, "Mean_rho_est", x = as.data.frame(RES[, , 9]))
writeDataTable(wb, "mean_neg_rho_est", x = as.data.frame(RES[, , 10]))
writeDataTable(wb, "Mean_tcc_kirk_est", x = as.data.frame(RES[, , 11]))

saveWorkbook(wb, file = "10000and50and25and225and03.xlsx")


##########                                            Graphic representation                                            ##########

alpha <- c(seq(from = 1, to = 2, length.out = 300),  seq(from = 2.001, to = 5, length.out = 300), 
           seq(from = 5.001, to = 10, length.out = 200), seq(from = 5.001, to = 30, length.out = 300),
           seq(from = 30.001, to = 100, length.out = 600), seq(from = 100.001, to = 1000, length.out = 300), rep(x = 1, times=2000))

beta <- c(rep(x = 1, times = 2000), seq(from = 1, to = 2, length.out = 300),  seq(from = 2.001, to = 5, length.out = 300), 
          seq(from = 5.001, to = 10, length.out = 200), seq(from = 5.001, to = 30, length.out = 300),
          seq(from = 30.001, to = 100, length.out = 600), seq(from = 100.001, to = 1000, length.out = 300)) 

x <- alpha/(alpha + beta)
y <- alpha*beta/((alpha + beta)^2*(alpha + beta + 1))
z <- y/(x*(1-x))

Max_ICC <- data.frame(Prevalence = x, Rho = z)


Res <- matrix(NA, 20, 9)
for (i in 1:9)
{ 
  d <- read.xlsx(xlsxFile = "10000and10and25and225and001.xlsx",
                 sheet = i, startRow = 2, colNames = FALSE)
  Res[ ,i] <- colMeans(d)
}

P_est <- read.xlsx(xlsxFile = "10000and10and25and225and001.xlsx",
                   sheet = 1, startRow = 2, colNames = FALSE)
Rho_est <- read.xlsx(xlsxFile = "10000and10and25and225and001.xlsx",
                     sheet = 7, startRow = 2, colNames = FALSE)

Rd_est <- c(); Inf_rho_max  <- c(); Prev_est <- c()
for (j in 1:s)
{  
  rho_max <- Max_ICC[which.min(abs(Max_ICC$Prevalence - prevalence[j])), ]$Rho
  ind_inf <- which(Rho_est[ ,j] <= rho_max)
  prev_est <- mean(P_est[ ,j][ind_inf])
  rho_inf_est <- mean(Rho_est[ ,j][ind_inf])
  rd_est <- 100*(rho_max - rho_inf_est)/rho_max
  inf_rho_max <- 100*mean(Rho_est[ ,j] <= rho_max)
  
  Rd_est <- c(Rd_est, rd_est)
  Inf_rho_max <- c(Inf_rho_max, inf_rho_max)
  Prev_est <- c(Prev_est, prev_est)
}


d10000and10and25and225and001 <- data.frame(Mean_p = Res[ ,1], Mean_vpc1 = Res[,2], Mean_vpc2 = Res[,3], Mean_vpc3 = Res[,4], 
                                           Mean_vpc4 =  Res[,5], Mean_mor = Res[,6], Mean_rho_est = Res[,7], Mean_tcc = Res[,8],
                                           Mean_tcc_kirk_est = Res[,9], Mean_prev_est = Prev_est, Mean_rd_est = Rd_est,
                                           Mean_inf_rho_max = Inf_rho_max)


Res <- matrix(NA, 20, 9)
for (i in 1:9)
{ 
  d <- read.xlsx(xlsxFile = "10000and20and25and225and001.xlsx",
                 sheet = i, startRow = 2, colNames = FALSE)
  Res[ ,i] <- colMeans(d)
}

P_est <- read.xlsx(xlsxFile = "10000and20and25and225and001.xlsx",
                   sheet = 1, startRow = 2, colNames = FALSE)
Rho_est <- read.xlsx(xlsxFile = "10000and20and25and225and001.xlsx",
                     sheet = 7, startRow = 2, colNames = FALSE)

Rd_est <- c(); Inf_rho_max  <- c(); Prev_est <- c()
for (j in 1:s)
{  
  rho_max <- Max_ICC[which.min(abs(Max_ICC$Prevalence - prevalence[j])), ]$Rho
  ind_inf <- which(Rho_est[ ,j] <= rho_max)
  prev_est <- mean(P_est[ ,j][ind_inf])
  rho_inf_est <- mean(Rho_est[ ,j][ind_inf])
  rd_est <- 100*(rho_max - rho_inf_est)/rho_max
  inf_rho_max <- 100*mean(Rho_est[ ,j] <= rho_max)
  
  Rd_est <- c(Rd_est, rd_est)
  Inf_rho_max <- c(Inf_rho_max, inf_rho_max)
  Prev_est <- c(Prev_est, prev_est)
}


d10000and20and25and225and001 <- data.frame(Mean_p = Res[ ,1], Mean_vpc1 = Res[,2], Mean_vpc2 = Res[,3], Mean_vpc3 = Res[,4], 
                                           Mean_vpc4 =  Res[,5], Mean_mor = Res[,6], Mean_rho_est = Res[,7], Mean_tcc = Res[,8],
                                           Mean_tcc_kirk_est = Res[,9], Mean_prev_est = Prev_est, Mean_rd_est = Rd_est,
                                           Mean_inf_rho_max = Inf_rho_max)


Res <- matrix(NA, 20, 9)
for (i in 1:9)
{ 
  d <- read.xlsx(xlsxFile = "10000and50and25and225and001.xlsx",
                 sheet = i, startRow = 2, colNames = FALSE)
  Res[ ,i] <- colMeans(d)
}

P_est <- read.xlsx(xlsxFile = "10000and50and25and225and001.xlsx",
                   sheet = 1, startRow = 2, colNames = FALSE)
Rho_est <- read.xlsx(xlsxFile = "10000and50and25and225and001.xlsx",
                     sheet = 7, startRow = 2, colNames = FALSE)

Rd_est <- c(); Inf_rho_max  <- c(); Prev_est <- c()
for (j in 1:s)
{  
  rho_max <- Max_ICC[which.min(abs(Max_ICC$Prevalence - prevalence[j])), ]$Rho
  ind_inf <- which(Rho_est[ ,j] <= rho_max)
  prev_est <- mean(P_est[ ,j][ind_inf])
  rho_inf_est <- mean(Rho_est[ ,j][ind_inf])
  rd_est <- 100*(rho_max - rho_inf_est)/rho_max
  inf_rho_max <- 100*mean(Rho_est[ ,j] <= rho_max)
  
  Rd_est <- c(Rd_est, rd_est)
  Inf_rho_max <- c(Inf_rho_max, inf_rho_max)
  Prev_est <- c(Prev_est, prev_est)
}


d10000and50and25and225and001 <- data.frame(Mean_p = Res[ ,1], Mean_vpc1 = Res[,2], Mean_vpc2 = Res[,3], Mean_vpc3 = Res[,4], 
                                           Mean_vpc4 =  Res[,5], Mean_mor = Res[,6], Mean_rho_est = Res[,7], Mean_tcc = Res[,8],
                                           Mean_tcc_kirk_est = Res[,9], Mean_prev_est = Prev_est, Mean_rd_est = Rd_est,
                                           Mean_inf_rho_max = Inf_rho_max)

Res <- matrix(NA, 20, 9)
for (i in 1:9)
{ 
  d <- read.xlsx(xlsxFile = "10000and10and25and225and005.xlsx",
                 sheet = i, startRow = 2, colNames = FALSE)
  Res[ ,i] <- colMeans(d)
}

P_est <- read.xlsx(xlsxFile = "10000and10and25and225and005.xlsx",
                   sheet = 1, startRow = 2, colNames = FALSE)
Rho_est <- read.xlsx(xlsxFile = "10000and10and25and225and005.xlsx",
                     sheet = 7, startRow = 2, colNames = FALSE)

Rd_est <- c(); Inf_rho_max  <- c(); Prev_est <- c()
for (j in 1:s)
{  
  rho_max <- Max_ICC[which.min(abs(Max_ICC$Prevalence - prevalence[j])), ]$Rho
  ind_inf <- which(Rho_est[ ,j] <= rho_max)
  prev_est <- mean(P_est[ ,j][ind_inf])
  rho_inf_est <- mean(Rho_est[ ,j][ind_inf])
  rd_est <- 100*(rho_max - rho_inf_est)/rho_max
  inf_rho_max <- 100*mean(Rho_est[ ,j] <= rho_max)
  
  Rd_est <- c(Rd_est, rd_est)
  Inf_rho_max <- c(Inf_rho_max, inf_rho_max)
  Prev_est <- c(Prev_est, prev_est)
}


d10000and10and25and225and005 <- data.frame(Mean_p = Res[ ,1], Mean_vpc1 = Res[,2], Mean_vpc2 = Res[,3], Mean_vpc3 = Res[,4], 
                                           Mean_vpc4 =  Res[,5], Mean_mor = Res[,6], Mean_rho_est = Res[,7], Mean_tcc = Res[,8],
                                           Mean_tcc_kirk_est = Res[,9], Mean_prev_est = Prev_est, Mean_rd_est = Rd_est,
                                           Mean_inf_rho_max = Inf_rho_max)


Res <- matrix(NA, 20, 9)
for (i in 1:9)
{ 
  d <- read.xlsx(xlsxFile = "10000and20and25and225and005.xlsx",
                 sheet = i, startRow = 2, colNames = FALSE)
  Res[ ,i] <- colMeans(d)
}

P_est <- read.xlsx(xlsxFile = "10000and20and25and225and005.xlsx",
                   sheet = 1, startRow = 2, colNames = FALSE)
Rho_est <- read.xlsx(xlsxFile = "10000and20and25and225and005.xlsx",
                     sheet = 7, startRow = 2, colNames = FALSE)

Rd_est <- c(); Inf_rho_max  <- c(); Prev_est <- c()
for (j in 1:s)
{  
  rho_max <- Max_ICC[which.min(abs(Max_ICC$Prevalence - prevalence[j])), ]$Rho
  ind_inf <- which(Rho_est[ ,j] <= rho_max)
  prev_est <- mean(P_est[ ,j][ind_inf])
  rho_inf_est <- mean(Rho_est[ ,j][ind_inf])
  rd_est <- 100*(rho_max - rho_inf_est)/rho_max
  inf_rho_max <- 100*mean(Rho_est[ ,j] <= rho_max)
  
  Rd_est <- c(Rd_est, rd_est)
  Inf_rho_max <- c(Inf_rho_max, inf_rho_max)
  Prev_est <- c(Prev_est, prev_est)
}


d10000and20and25and225and005 <- data.frame(Mean_p = Res[ ,1], Mean_vpc1 = Res[,2], Mean_vpc2 = Res[,3], Mean_vpc3 = Res[,4], 
                                           Mean_vpc4 =  Res[,5], Mean_mor = Res[,6], Mean_rho_est = Res[,7], Mean_tcc = Res[,8],
                                           Mean_tcc_kirk_est = Res[,9], Mean_prev_est = Prev_est, Mean_rd_est = Rd_est,
                                           Mean_inf_rho_max = Inf_rho_max)


Res <- matrix(NA, 20, 9)
for (i in 1:9)
{ 
  d <- read.xlsx(xlsxFile = "10000and50and25and225and005.xlsx",
                 sheet = i, startRow = 2, colNames = FALSE)
  Res[ ,i] <- colMeans(d)
}

P_est <- read.xlsx(xlsxFile = "10000and50and25and225and005.xlsx",
                   sheet = 1, startRow = 2, colNames = FALSE)
Rho_est <- read.xlsx(xlsxFile = "10000and50and25and225and005.xlsx",
                     sheet = 7, startRow = 2, colNames = FALSE)

Rd_est <- c(); Inf_rho_max  <- c(); Prev_est <- c()
for (j in 1:s)
{  
  rho_max <- Max_ICC[which.min(abs(Max_ICC$Prevalence - prevalence[j])), ]$Rho
  ind_inf <- which(Rho_est[ ,j] <= rho_max)
  prev_est <- mean(P_est[ ,j][ind_inf])
  rho_inf_est <- mean(Rho_est[ ,j][ind_inf])
  rd_est <- 100*(rho_max - rho_inf_est)/rho_max
  inf_rho_max <- 100*mean(Rho_est[ ,j] <= rho_max)
  
  Rd_est <- c(Rd_est, rd_est)
  Inf_rho_max <- c(Inf_rho_max, inf_rho_max)
  Prev_est <- c(Prev_est, prev_est)
}


d10000and50and25and225and005 <- data.frame(Mean_p = Res[ ,1], Mean_vpc1 = Res[,2], Mean_vpc2 = Res[,3], Mean_vpc3 = Res[,4], 
                                           Mean_vpc4 =  Res[,5], Mean_mor = Res[,6], Mean_rho_est = Res[,7], Mean_tcc = Res[,8],
                                           Mean_tcc_kirk_est = Res[,9], Mean_prev_est = Prev_est, Mean_rd_est = Rd_est,
                                           Mean_inf_rho_max = Inf_rho_max)


Res <- matrix(NA, 20, 9)
for (i in 1:9)
{ 
  d <- read.xlsx(xlsxFile = "10000and10and25and225and03.xlsx", sheet = i, startRow = 2, colNames = FALSE)
  Res[ ,i] <- colMeans(d)
}

P_est <- read.xlsx(xlsxFile = "10000and10and25and225and03.xlsx", sheet = 1, startRow = 2, colNames = FALSE)

Rho_est <- read.xlsx(xlsxFile = "10000and10and25and225and03.xlsx", sheet = 7, startRow = 2, colNames = FALSE)

Rd_est <- c(); Inf_rho_max  <- c(); Prev_est <- c()
for (j in 1:s)
{  
  rho_max <- Max_ICC[which.min(abs(Max_ICC$Prevalence - prevalence[j])), ]$Rho
  ind_inf <- which(Rho_est[ ,j] <= rho_max)
  prev_est <- mean(P_est[ ,j][ind_inf])
  rho_inf_est <- mean(Rho_est[ ,j][ind_inf])
  rd_est <- 100*(rho_max - rho_inf_est)/rho_max
  inf_rho_max <- 100*mean(Rho_est[ ,j] <= rho_max)
  
  Rd_est <- c(Rd_est, rd_est)
  Inf_rho_max <- c(Inf_rho_max, inf_rho_max)
  Prev_est <- c(Prev_est, prev_est)
}


d10000and10and25and225and03 <- data.frame(Mean_p = Res[ ,1], Mean_vpc1 = Res[,2], Mean_vpc2 = Res[,3], Mean_vpc3 = Res[,4], 
                                          Mean_vpc4 =  Res[,5], Mean_mor = Res[,6], Mean_rho_est = Res[,7], Mean_tcc = Res[,8],
                                          Mean_tcc_kirk_est = Res[,9], Mean_prev_est = Prev_est, Mean_rd_est = Rd_est,
                                          Mean_inf_rho_max = Inf_rho_max)


Res <- matrix(NA, 20, 9)
for (i in 1:9)
{ 
  d <- read.xlsx(xlsxFile = "10000and20and25and225and03.xlsx",
                 sheet = i, startRow = 2, colNames = FALSE)
  Res[ ,i] <- colMeans(d)
}

P_est <- read.xlsx(xlsxFile = "10000and20and25and225and03.xlsx",
                   sheet = 1, startRow = 2, colNames = FALSE)
Rho_est <- read.xlsx(xlsxFile = "10000and20and25and225and03.xlsx",
                     sheet = 7, startRow = 2, colNames = FALSE)

Rd_est <- c(); Inf_rho_max  <- c(); Prev_est <- c()
for (j in 1:s)
{  
  rho_max <- Max_ICC[which.min(abs(Max_ICC$Prevalence - prevalence[j])), ]$Rho
  ind_inf <- which(Rho_est[ ,j] <= rho_max)
  prev_est <- mean(P_est[ ,j][ind_inf])
  rho_inf_est <- mean(Rho_est[ ,j][ind_inf])
  rd_est <- 100*(rho_max - rho_inf_est)/rho_max
  inf_rho_max <- 100*mean(Rho_est[ ,j] <= rho_max)
  
  Rd_est <- c(Rd_est, rd_est)
  Inf_rho_max <- c(Inf_rho_max, inf_rho_max)
  Prev_est <- c(Prev_est, prev_est)
}


d10000and20and25and225and03 <- data.frame(Mean_p = Res[ ,1], Mean_vpc1 = Res[,2], Mean_vpc2 = Res[,3], Mean_vpc3 = Res[,4], 
                                          Mean_vpc4 =  Res[,5], Mean_mor = Res[,6], Mean_rho_est = Res[,7], Mean_tcc = Res[,8],
                                          Mean_tcc_kirk_est = Res[,9], Mean_prev_est = Prev_est, Mean_rd_est = Rd_est,
                                          Mean_inf_rho_max = Inf_rho_max)


Res <- matrix(NA, 20, 9)
for (i in 1:9)
{ 
  d <- read.xlsx(xlsxFile = "10000and50and25and225and03.xlsx",
                 sheet = i, startRow = 2, colNames = FALSE)
  Res[ ,i] <- colMeans(d)
}

P_est <- read.xlsx(xlsxFile = "10000and50and25and225and03.xlsx",
                   sheet = 1, startRow = 2, colNames = FALSE)
Rho_est <- read.xlsx(xlsxFile = "10000and50and25and225and03.xlsx",
                     sheet = 7, startRow = 2, colNames = FALSE)

Rd_est <- c(); Inf_rho_max  <- c(); Prev_est <- c()
for (j in 1:s)
{  
  rho_max <- Max_ICC[which.min(abs(Max_ICC$Prevalence - prevalence[j])), ]$Rho
  ind_inf <- which(Rho_est[ ,j] <= rho_max)
  prev_est <- mean(P_est[ ,j][ind_inf])
  rho_inf_est <- mean(Rho_est[ ,j][ind_inf])
  rd_est <- 100*(rho_max - rho_inf_est)/rho_max
  inf_rho_max <- 100*mean(Rho_est[ ,j] <= rho_max)
  
  Rd_est <- c(Rd_est, rd_est)
  Inf_rho_max <- c(Inf_rho_max, inf_rho_max)
  Prev_est <- c(Prev_est, prev_est)
}


d10000and50and25and225and03 <- data.frame(Mean_p = Res[ ,1], Mean_vpc1 = Res[,2], Mean_vpc2 = Res[,3], Mean_vpc3 = Res[,4], 
                                          Mean_vpc4 =  Res[,5], Mean_mor = Res[,6], Mean_rho_est = Res[,7], Mean_tcc = Res[,8],
                                          Mean_tcc_kirk_est = Res[,9], Mean_prev_est = Prev_est, Mean_rd_est = Rd_est,
                                          Mean_inf_rho_max = Inf_rho_max)


d10000and10and25and225and001$Mean_p <- seq(from = 0.01, to = 0.99, length.out = 20)
d10000and20and25and225and001$Mean_p <- seq(from = 0.01, to = 0.99, length.out = 20)
d10000and50and25and225and001$Mean_p <- seq(from = 0.01, to = 0.99, length.out = 20)
d10000and10and25and225and005$Mean_p <- seq(from = 0.01, to = 0.99, length.out = 20)
d10000and20and25and225and005$Mean_p <- seq(from = 0.01, to = 0.99, length.out = 20)
d10000and50and25and225and005$Mean_p <- seq(from = 0.01, to = 0.99, length.out = 20)
d10000and10and25and225and03$Mean_p <- seq(from = 0.01, to = 0.99, length.out = 20)
d10000and20and25and225and03$Mean_p <- seq(from = 0.01, to = 0.99, length.out = 20)
d10000and50and25and225and03$Mean_p <- seq(from = 0.01, to = 0.99, length.out = 20)

##########     Graphic representation for VPCs     ##########

p1 <- ggplot(data = d10000and10and25and225and001) +
  geom_line(aes(x = Mean_p, y = Mean_vpc1, col = "Estimated VPC1"), size = 1.5) +
  geom_line(aes(x = Mean_p, y = Mean_vpc2, col = "Estimated VPC2"), size = 1.5) +
  geom_line(aes(x = Mean_p, y = Mean_vpc3, col = "Estimated VPC3"), size = 1.5) + 
  geom_line(aes(x = Mean_p, y = Mean_vpc4, col = "Estimated VPC4"), size = 1.5) +
  xlab("Prevalence") + ylab("VPC") + scale_colour_manual(name = "", values = c("#00BA38", "#9683EC", "#619CFF", "#F8766D"),
                                                         breaks = c("Estimated VPC1","Estimated VPC2","Estimated VPC3", "Estimated VPC4"), 
                                                         labels = c(expression(VPC[1]), expression(VPC[2]), expression(VPC[3]), expression(VPC[4]))) 

p2 <- ggplot(data = d10000and20and25and225and001) +
  geom_line(aes(x = Mean_p, y = Mean_vpc1), size = 1.5, col = "#00BA38") +
  geom_line(aes(x = Mean_p, y = Mean_vpc2), size = 1.5, col = "#9683EC") +
  geom_line(aes(x = Mean_p, y = Mean_vpc3), size = 1.5, col = "#619CFF") + 
  geom_line(aes(x = Mean_p, y = Mean_vpc4), size = 1.5, col = "#F8766D") +
  xlab("Prevalence") + ylab("VPC") 

p3 <- ggplot(data = d10000and50and25and225and001) +
  geom_line(aes(x = Mean_p, y = Mean_vpc1), size = 1.5, col = "#00BA38") +
  geom_line(aes(x = Mean_p, y = Mean_vpc2), size = 1.5, col = "#9683EC") +
  geom_line(aes(x = Mean_p, y = Mean_vpc3), size = 1.5, col = "#619CFF") + 
  geom_line(aes(x = Mean_p, y = Mean_vpc4), size = 1.5, col = "#F8766D")  +
  xlab("Prevalence") + ylab("VPC") 

p4 <- ggplot(data = d10000and10and25and225and005) +
  geom_line(aes(x = Mean_p, y = Mean_vpc1), size = 1.5, col = "#00BA38") +
  geom_line(aes(x = Mean_p, y = Mean_vpc2), size = 1.5, col = "#9683EC") +
  geom_line(aes(x = Mean_p, y = Mean_vpc3), size = 1.5, col = "#619CFF") + 
  geom_line(aes(x = Mean_p, y = Mean_vpc4), size = 1.5, col = "#F8766D") +
  xlab("Prevalence") + ylab("VPC") 

p5 <- ggplot(data = d10000and20and25and225and005) +
  geom_line(aes(x = Mean_p, y = Mean_vpc1), size = 1.5, col = "#00BA38") +
  geom_line(aes(x = Mean_p, y = Mean_vpc2), size = 1.5, col = "#9683EC") +
  geom_line(aes(x = Mean_p, y = Mean_vpc3), size = 1.5, col = "#619CFF") + 
  geom_line(aes(x = Mean_p, y = Mean_vpc4), size = 1.5, col = "#F8766D") +
  xlab("Prevalence") + ylab("VPC") 

p6 <- ggplot(data = d10000and50and25and225and005) +
  geom_line(aes(x = Mean_p, y = Mean_vpc1), size = 1.5, col = "#00BA38") +
  geom_line(aes(x = Mean_p, y = Mean_vpc2), size = 1.5, col = "#9683EC") +
  geom_line(aes(x = Mean_p, y = Mean_vpc3), size = 1.5, col = "#619CFF") + 
  geom_line(aes(x = Mean_p, y = Mean_vpc4), size = 1.5, col = "#F8766D") +
  xlab("Prevalence") + ylab("VPC") 

p7 <- ggplot(data = d10000and10and25and225and03) +
  geom_line(aes(x = Mean_p, y = Mean_vpc1), size = 1.5, col = "#00BA38") +
  geom_line(aes(x = Mean_p, y = Mean_vpc2), size = 1.5, col = "#9683EC") +
  geom_line(aes(x = Mean_p, y = Mean_vpc3), size = 1.5, col = "#619CFF") + 
  geom_line(aes(x = Mean_p, y = Mean_vpc4), size = 1.5, col = "#F8766D") +
  xlab("Prevalence") + ylab("VPC") 

p8 <- ggplot(data = d10000and20and25and225and03) +
  geom_line(aes(x = Mean_p, y = Mean_vpc1), size = 1.5, col = "#00BA38") +
  geom_line(aes(x = Mean_p, y = Mean_vpc2), size = 1.5, col = "#9683EC") +
  geom_line(aes(x = Mean_p, y = Mean_vpc3), size = 1.5, col = "#619CFF") + 
  geom_line(aes(x = Mean_p, y = Mean_vpc4), size = 1.5, col = "#F8766D") +
  xlab("Prevalence") + ylab("VPC") 


p9 <- ggplot(data = d10000and50and25and225and03) +
  geom_line(aes(x = Mean_p, y = Mean_vpc1), size = 1.5, col = "#00BA38") +
  geom_line(aes(x = Mean_p, y = Mean_vpc2), size = 1.5, col = "#9683EC") +
  geom_line(aes(x = Mean_p, y = Mean_vpc3), size = 1.5, col = "#619CFF") + 
  geom_line(aes(x = Mean_p, y = Mean_vpc4), size = 1.5, col = "#F8766D") +
  xlab("Prevalence") + ylab("VPC") 


p10 <- p1 + theme(legend.position = "none")

p11 <- get_legend(p1 + theme(legend.position=c(0.60,0.6), legend.justification = "center"))

plot_grid(p11, NULL, NULL, NULL, NULL, p10, p2, p3, NULL, p4, p5, p6, NULL, p7, p8, p9, ncol = 4, nrow = 4, 
          rel_widths = c(1, 4, 4, 4), rel_heights = c(1, 2, 2, 2), 
          labels = c("", "k = 10", "k = 20", "k = 50", "(A)", "", "", "", "(B)", "", "", "", "(C)", "", "",""), 
          label_size = 20, label_x = 0.55,  label_y = 0.65, hjust = 0.5)


##########     Graphic representation for MOR     ##########

p1 <- ggplot(data = d10000and10and25and225and001) + geom_line(aes(x = Mean_p, y = Mean_mor), size = 1.5) +
  xlab("Prevalence") + ylab("Median odds ratio") + scale_x_continuous(breaks=c(0.1, 0.3, 0.5, 0.7, 0.9)) 

p2 <- ggplot(data = d10000and20and25and225and001) + geom_line(aes(x = Mean_p, y = Mean_mor), size = 1.5) +
  xlab("Prevalence") + ylab("Median odds ratio") + scale_x_continuous(breaks=c(0.1, 0.3, 0.5, 0.7, 0.9)) 

p3 <- ggplot(data = d10000and50and25and225and001) + geom_line(aes(x = Mean_p, y = Mean_mor), size = 1.5) +
  xlab("Prevalence") + ylab("Median odds ratio") 

p4 <- ggplot(data = d10000and10and25and225and005) + geom_line(aes(x = Mean_p, y = Mean_mor), size = 1.5) +
  xlab("Prevalence") + ylab("Median odds ratio") + scale_x_continuous(breaks=c(0.1, 0.3, 0.5, 0.7, 0.9))


p5 <- ggplot(data = d10000and20and25and225and005) + geom_line(aes(x = Mean_p, y = Mean_mor), size = 1.5) +
  xlab("Prevalence") + ylab("Median odds ratio") + scale_x_continuous(breaks=c(0.1, 0.3, 0.5, 0.7, 0.9)) 

p6 <- ggplot(data = d10000and50and25and225and005) + geom_line(aes(x = Mean_p, y = Mean_mor), size = 1.5) +
  xlab("Prevalence") + ylab("Median odds ratio") + scale_x_continuous(breaks=c(0.1, 0.3, 0.5, 0.7, 0.9)) 

p7 <- ggplot(data = d10000and10and25and225and03) + geom_line(aes(x = Mean_p, y = Mean_mor), size = 1.5) +
  xlab("Prevalence") + ylab("Median odds ratio") + scale_x_continuous(breaks=c(0.1, 0.3, 0.5, 0.7, 0.9)) 

p8 <- ggplot(data = d10000and20and25and225and03) + geom_line(aes(x = Mean_p, y = Mean_mor), size = 1.5) +
  xlab("Prevalence") + ylab("Median odds ratio") + scale_x_continuous(breaks=c(0.1, 0.3, 0.5, 0.7, 0.9)) 

p9 <- ggplot(data = d10000and50and25and225and03) + geom_line(aes(x = Mean_p, y = Mean_mor), size = 1.5) +
  xlab("Prevalence") + ylab("Median odds ratio") + scale_x_continuous(breaks=c(0.1, 0.3, 0.5, 0.7, 0.9))    

plot_grid(NULL, NULL, NULL, NULL, NULL, p1, p2, p3, NULL, p4, p5, p6, NULL, p7, p8, p9, ncol = 4, nrow = 4, 
          rel_widths = c(0.4, 2, 2, 2), rel_heights = c(0.4, 2, 2, 2), labels = c("", "k = 10", "k = 20", "k = 50", 
                                                                                  "(A)", "", "", "", "(B)", "", "", "", "(C)", "", "",""), label_size = 20, label_x = 0.55,  label_y = 0.65, hjust=0.5)



##########     Graphic representation for TCCs     ##########

p1 <- ggplot(data = d10000and10and25and225and001) +
  geom_line(aes(x = Mean_p, y = Mean_rho_est, col = "Estimated binary ICC"), size = 1.5) +
  geom_line(aes(x = Mean_p, y = Mean_tcc_kirk_est, col = "Estimated TCC (Kirk formula)"), size = 1.5) +
  geom_line(aes(x = Mean_p, y = Mean_tcc, col = "Estimated TCC (Original formula)"), size = 1.5) +
  geom_hline(aes(yintercept = 0.01, col = "Theoretical continuous ICC"), size = 1.5) +
  xlab("Prevalence") + ylab("") + scale_colour_manual(name="", values=c("#00BA38", "#9683EC", "#619CFF", "#F8766D"))

p2 <- ggplot(data = d10000and20and25and225and001) + 
  geom_line(aes(x = Mean_p, y = Mean_rho_est), size = 1.5, col = "#00BA38") +
  geom_line(aes(x = Mean_p, y = Mean_tcc_kirk_est), size = 1.5, col = "#9683EC") +
  geom_line(aes(x = Mean_p, y = Mean_tcc), size = 1.5, col = "#619CFF") +
  geom_hline(aes(yintercept = 0.01), size = 1.5, col = "#F8766D") + xlab("Prevalence") + ylab("") 

p3 <- ggplot(data = d10000and50and25and225and001) + 
  geom_line(aes(x = Mean_p, y = Mean_rho_est), size = 1.5, col = "#00BA38") +
  geom_line(aes(x = Mean_p, y = Mean_tcc_kirk_est), size = 1.5, col = "#9683EC") +
  geom_line(aes(x = Mean_p, y = Mean_tcc), size = 1.5, col = "#619CFF") +
  geom_hline(aes(yintercept = 0.01), size = 1.5, col = "#F8766D") + xlab("Prevalence") + ylab("") 

p4 <- ggplot(data = d10000and10and25and225and005) + 
  geom_line(aes(x = Mean_p, y = Mean_rho_est), size = 1.5, col = "#00BA38") +
  geom_line(aes(x = Mean_p, y = Mean_tcc_kirk_est), size = 1.5, col = "#9683EC") +
  geom_line(aes(x = Mean_p, y = Mean_tcc), size = 1.5, col = "#619CFF") +
  geom_hline(aes(yintercept = 0.05), size = 1.5, col = "#F8766D") + xlab("Prevalence") + ylab("")  

p5 <- ggplot(data = d10000and20and25and225and005) + 
  geom_line(aes(x = Mean_p, y = Mean_rho_est), size = 1.5, col = "#00BA38") +
  geom_line(aes(x = Mean_p, y = Mean_tcc_kirk_est), size = 1.5, col = "#9683EC") +
  geom_line(aes(x = Mean_p, y = Mean_tcc), size = 1.5, col = "#619CFF") +
  geom_hline(aes(yintercept = 0.05), size = 1.5, col = "#F8766D") + xlab("Prevalence") + ylab("")  

p6 <- ggplot(data = d10000and50and25and225and005) + 
  geom_line(aes(x = Mean_p, y = Mean_rho_est), size = 1.5, col = "#00BA38") +
  geom_line(aes(x = Mean_p, y = Mean_tcc_kirk_est), size = 1.5, col = "#9683EC") +
  geom_line(aes(x = Mean_p, y = Mean_tcc), size = 1.5, col = "#619CFF") +
  geom_hline(aes(yintercept = 0.05), size = 1.5, col = "#F8766D") + xlab("Prevalence") + ylab("") 
p7 <- ggplot(data = d10000and10and25and225and03) + 
  geom_line(aes(x = Mean_p, y = Mean_rho_est), size = 1.5, col = "#00BA38") +
  geom_line(aes(x = Mean_p, y = Mean_tcc_kirk_est), size = 1.5, col = "#9683EC") +
  geom_line(aes(x = Mean_p, y = Mean_tcc), size = 1.5, col = "#619CFF") +
  geom_hline(aes(yintercept = 0.3), size = 1.5, col = "#F8766D") + xlab("Prevalence") + ylab("") 

p8 <- ggplot(data = d10000and20and25and225and03) + 
  geom_line(aes(x = Mean_p, y = Mean_rho_est), size = 1.5, col = "#00BA38") +
  geom_line(aes(x = Mean_p, y = Mean_tcc_kirk_est), size = 1.5, col = "#9683EC") +
  geom_line(aes(x = Mean_p, y = Mean_tcc), size = 1.5, col = "#619CFF") +
  geom_hline(aes(yintercept = 0.3), size = 1.5, col = "#F8766D") + xlab("Prevalence") + ylab("") 

p9 <- ggplot(data = d10000and50and25and225and03) + 
  geom_line(aes(x = Mean_p, y = Mean_rho_est), size = 1.5, col = "#00BA38") +
  geom_line(aes(x = Mean_p, y = Mean_tcc_kirk_est), size = 1.5, col = "#9683EC") +
  geom_line(aes(x = Mean_p, y = Mean_tcc), size = 1.5, col = "#619CFF") +
  geom_hline(aes(yintercept = 0.3), size = 1.5, col = "#F8766D") + xlab("Prevalence") + ylab("")  

p10 <- p1 + theme(legend.position = "none")

p11 <- get_legend(p1 + theme(legend.position=c(0.98,0.6), legend.justification = "center"))

plot_grid(p11, NULL, NULL, NULL, NULL, p10, p2, p3, NULL, p4, p5, p6, NULL, p7, p8, p9, ncol = 4, nrow = 4, 
          rel_widths = c(1, 2, 2, 2), rel_heights = c(1, 2, 2, 2), labels=c("", "k = 10", "k = 20", "k = 50", "(A)", 
                                                                            "", "", "", "(B)", "", "", "", "(C)", "", "",""), label_size = 20, label_x = 0.7,  label_y = 0.65, hjust=0.5)



##########     Graphic representation for the relative deviation to the theoretical maximum ICC value     ##########

p1 <- ggplot(data = d10000and10and25and225and001) + geom_line(aes(x = Mean_prev_est, y = Mean_rd_est), size = 1.5) + 
  geom_point(aes(x = Mean_prev_est, y = Mean_inf_rho_max), col = "red", pch = 15, size = 1.5) + 
  xlab("Prevalence") + ylab("Distance to the \n max ICC value (%)") 

p2 <- ggplot(data = d10000and20and25and225and001)  + geom_line(aes(x = Mean_prev_est, y = Mean_rd_est), size = 1.5) + 
  geom_point(aes(x = Mean_prev_est, y = Mean_inf_rho_max), col = "red", pch = 15, size = 1.5) + 
  xlab("Prevalence") + ylab("Distance to the \n max ICC value (%)") 

p3 <- ggplot(data = d10000and50and25and225and001)  + geom_line(aes(x = Mean_prev_est, y = Mean_rd_est), size = 1.5) + 
  geom_point(aes(x = Mean_prev_est, y = Mean_inf_rho_max), col = "red", pch = 15, size = 1.5) + 
  xlab("Prevalence") + ylab("Distance to the \n max ICC value (%)") 

p4 <- ggplot(data = d10000and10and25and225and005)  + geom_line(aes(x = Mean_prev_est, y = Mean_rd_est), size = 1.5) + 
  geom_point(aes(x = Mean_prev_est, y = Mean_inf_rho_max), col = "red", pch = 15, size = 1.5) + 
  xlab("Prevalence") + ylab("Distance to the \n max ICC value (%)") 

p5 <- ggplot(data = d10000and20and25and225and005)  + geom_line(aes(x = Mean_prev_est, y = Mean_rd_est), size = 1.5) + 
  geom_point(aes(x = Mean_prev_est, y = Mean_inf_rho_max), col = "red", pch = 15, size = 1.5) + 
  xlab("Prevalence") + ylab("Distance to the \n max ICC value (%)") 

p6 <- ggplot(data = d10000and50and25and225and005)  + geom_line(aes(x = Mean_prev_est, y = Mean_rd_est), size = 1.5) + 
  geom_point(aes(x = Mean_prev_est, y = Mean_inf_rho_max), col = "red", pch = 15, size = 1.5) + 
  xlab("Prevalence") + ylab("Distance to the \n max ICC value (%)") 

p7 <- ggplot(data = d10000and10and25and225and03)  + geom_line(aes(x = Mean_prev_est, y = Mean_rd_est), size = 1.5) + 
  geom_point(aes(x = Mean_prev_est, y = Mean_inf_rho_max), col = "red", pch = 15, size = 1.5) + 
  xlab("Prevalence") + ylab("Distance to the \n max ICC value (%)") 

p8 <- ggplot(data = d10000and20and25and225and03)  + geom_line(aes(x = Mean_prev_est, y = Mean_rd_est), size = 1.5) + 
  geom_point(aes(x = Mean_prev_est, y = Mean_inf_rho_max), col = "red", pch = 15, size = 1.5) + 
  xlab("Prevalence") + ylab("Distance to the \n max ICC value (%)") 

p9 <- ggplot(data = d10000and50and25and225and03)  + geom_line(aes(x = Mean_prev_est, y = Mean_rd_est), size = 1.5) + 
  geom_point(aes(x = Mean_prev_est, y = Mean_inf_rho_max), col = "red", pch = 15, size = 1.5) + 
  xlab("Prevalence") + ylab("Distance to the \n max ICC value (%)") 

plot_grid(NULL, NULL, NULL, NULL, NULL, p1, p2, p3, NULL, p4, p5, p6, NULL, p7, p8, p9, ncol = 4, nrow = 4, 
          rel_widths = c(0.4, 2, 2, 2), rel_heights = c(0.4, 2, 2, 2), labels = c("", "k = 10", "k = 20", "k = 50", 
                                                                                  "(A)", "", "", "", "(B)", "", "", "", "(C)", "", "",""), label_size = 20, label_x = 0.55,  label_y = 0.65, hjust=0.5)

##########     Proportion of null between-cluster variance     ##########

f0 <- function(x) length(which(x == 0))*100/length(x)

RESULTS <- data.frame(Theo_P = round(seq(from = 0.01, to = 0.99, length.out = 20), 2),
                      S1 = apply(X = read.xlsx(xlsxFile = "10000and10and25and225and001.xlsx",
                                               sheet = 2, startRow = 2, colNames = FALSE), MARGIN = 2, FUN = f0),
                      S2 = apply(X = read.xlsx(xlsxFile = "10000and20and25and225and001.xlsx",
                                               sheet = 2, startRow = 2, colNames = FALSE), MARGIN = 2, FUN = f0),
                      S3 = apply(X = read.xlsx(xlsxFile = "10000and50and25and225and001.xlsx",
                                               sheet = 2, startRow = 2, colNames = FALSE), MARGIN = 2, FUN = f0),
                      S4 = apply(X = read.xlsx(xlsxFile = "10000and10and25and225and005.xlsx",
                                               sheet = 2, startRow = 2, colNames = FALSE), MARGIN = 2, FUN = f0),
                      S5 = apply(X = read.xlsx(xlsxFile = "10000and20and25and225and005.xlsx",
                                               sheet = 2, startRow = 2, colNames = FALSE), MARGIN = 2, FUN = f0),
                      S6 = apply(X = read.xlsx(xlsxFile = "10000and50and25and225and005.xlsx",
                                               sheet = 2, startRow = 2, colNames = FALSE), MARGIN = 2, FUN = f0),
                      S7 = apply(X = read.xlsx(xlsxFile = "10000and10and25and225and03.xlsx",
                                                 sheet = 2, startRow = 2, colNames = FALSE), MARGIN = 2, FUN = f0),
                      S8 = apply(X = read.xlsx(xlsxFile = "10000and20and25and225and03.xlsx",
                                               sheet = 2, startRow = 2, colNames = FALSE), MARGIN = 2, FUN = f0),
                      S9 = apply(X = read.xlsx(xlsxFile = "10000and50and25and225and03.xlsx",
                                               sheet = 2, startRow = 2, colNames = FALSE), MARGIN = 2, FUN = f0))
RESULTS



##########                               Maximal Information Coefficient calculation                               ##########

f_mic <- function(dat)
{
  cor_rho_est <- round(mine(dat$Mean_p, dat$Mean_rho_est)$MIC, digits = 2)
  cor_vpc1 <- round(mine(dat$Mean_p, dat$Mean_vpc1)$MIC, digits = 2)
  cor_vpc2 <- round(mine(dat$Mean_p, dat$Mean_vpc2)$MIC, digits = 2)
  cor_vpc4 <- round(mine(dat$Mean_p, dat$Mean_vpc4)$MIC, digits = 2)
  cor_mor <- round(mine(dat$Mean_p, dat$Mean_mor)$MIC, digits = 2)
  cor_tcc_o <- round(mine(dat$Mean_p, dat$Mean_tcc)$MIC, digits = 2)
  cor_tcc_k <- round(mine(dat$Mean_p, dat$Mean_tcc_kirk_est)$MIC, digits = 2)
  cor_rd <- round(mine(dat$Mean_p, dat$Mean_rd_est)$MIC, digits = 2)
  
  return(list(cor_rho_est = cor_rho_est, cor_vpc1 = cor_vpc1, cor_vpc2 = cor_vpc2, cor_vpc4 = cor_vpc4, cor_mor = cor_mor,
              cor_tcc_o = cor_tcc_o, cor_tcc_k = cor_tcc_k, cor_rd = cor_rd))
}

RES <- matrix(NA, nrow = 9, ncol = 8)
RES[1, ] <- unlist(f_mic(dat = d10000and10and25and225and001))
RES[2, ] <- unlist(f_mic(dat = d10000and20and25and225and001))
RES[3, ] <- unlist(f_mic(dat = d10000and50and25and225and001))
RES[4, ] <- unlist(f_mic(dat = d10000and10and25and225and005))
RES[5, ] <- unlist(f_mic(dat = d10000and20and25and225and005))
RES[6, ] <- unlist(f_mic(dat = d10000and50and25and225and005))
RES[7, ] <- unlist(f_mic(dat = d10000and10and25and225and03))
RES[8, ] <- unlist(f_mic(dat = d10000and20and25and225and03))
RES[9, ] <- unlist(f_mic(dat = d10000and50and25and225and03))
RES


#########              Caluculation of intracluster correlation measures on Pythagore data              ##########

f.meas <- function(out, g, dat, est_p, nsim)
{
  l <- list(g = substitute(g), out = substitute(out))
  d <- data.frame(g = eval(expr = l$g, envir = dat), out=eval(expr = l$out, envir = dat))
  
  # Fit models
  # m1 <- glmer(out ~ 1 + 1|g, data = d, family = binomial(link = "logit"), nAGQ = 20)
  
  m1 <- glmer(out ~ 1 + 1|g, data = d, family = binomial(link = "logit"), nAGQ = 20,
              control = glmerControl(check.conv.grad = .makeCC("stop", tol = 2e-3, relTol = NULL)))
  
  ## VPC: method A
  # estimate of intercept
  beta0 <- as.numeric(fixef(m1))
  
  # estimate of the between group variance
  sigma_carre_u0 <- as.numeric(VarCorr(m1)[[1]])
  
  # estimate of the prevalence
  p1i_ij <- prop.table(table(d$out))[2]
  
  # estimate of the level one variance
  v11 <- p1i_ij*(1 - p1i_ij)
  
  # estimate of the level two variance
  v12 <- (sigma_carre_u0 *p1i_ij^2)/(1 + exp(beta0))^2
  
  # estimate of the variance partition coefficient
  vpc1 <- v12/(v12 + v11)
  
  
  ## VPC: method B
  # simulate 
  u0j <- rnorm(n = nsim, mean = 0, sd = sqrt(as.numeric(VarCorr(m1)[[1]])))
  
  # estimate of intercept
  beta0 <- as.numeric(fixef(m1))
  
  # estimate of the prevalence
  p2i_ij <- exp(beta0 + u0j)/(1 + exp(beta0 + u0j))
  
  # estimate of the level one variance
  v21 <- mean(p2i_ij*(1-p2i_ij))
  
  # estimate of the level two variance
  v22 <- var(p2i_ij)
  
  # estimate of the variance partition coefficient
  vpc2 <- v22/(v22 + v21)
  
  
  ## VPC: method C
  
  # estimate of the variance partition coefficient
  vpc3 <- ICCest(x = g, y = as.numeric(out), data = d)$ICC
  
  if(vpc3 < 0) {vpc3 <- 0}
  
  ## VPC: method D
  # the level one variance is 1
  
  # estimate of the level two variance
  v42 <- as.numeric(VarCorr(m1)[[1]])
  
  # estimate of the variance partition coefficient
  vpc4 <- v42/(v42 + pi^2/3)
  
  
  ## MOR
  # the level two variance is sigma_carre_u0
  
  # estimate of the median odds ratio
  mor <- exp(qnorm(0.75)*sqrt(2*sigma_carre_u0))
  
  tcc_original <- approx_TCC.1(mat = table(make_pairs(d[ , c("g","out")])), delta = 0.00001)$rtet
  
  rho_est <- est_anova(g = g, out = out, data = d)$Rho
  
  tcc_kirk <- rec_sum_f_int(rhobin = est_anova(g = g, out = out, data = d)$Rho, n = 10, h = qnorm(est_p), p = est_p)
  
  rho_max <- Max_ICC[which.min(abs(Max_ICC$Prevalence - est_p)), ]$Rho
  rd_est <- 100*(rho_max - rho_est)/rho_max
  
  return(list(vpc1 = vpc1, vpc2 = vpc2, vpc3 = vpc3, vpc4 = vpc4, mor = mor,
              tcc_original = tcc_original, rho_est = rho_est, tcc_kirk = tcc_kirk, rd_est = rd_est))
}


set.seed(1234)
f.meas(out = femme_nal, g = num_centre, dat = data_control, est_p = 0.459, nsim = 5000)
f.meas(out = femme_nal, g = num_centre, dat = data_intervention, est_p = 0.540, nsim = 5000)


#########  Confidence interval

B <- 1000
RES <- matrix(NA, nrow = 1000, ncol = 9)
for (b in 1:B)
{
  d_boot <- c()
  clust_sp <- unique(sample(unique(data_control$num_centre), replace = T))
  
  for (i in 1:length(clust_sp))
  {
    res <- d_control[data_control$num_centre == clust_sp[i], ]
    d_boot <- rbind(d_boot, res)
  }
  RES[b, ] <- as.numeric(f.meas(out = femme_nal, g = factor(num_centre), est_p = as.numeric(prop.table(table(d_boot$femme_nal))[2]),
                                dat = d_boot, nsim = 5000))
}

quantile(RES[, 1], probs = 0.025); quantile(RES[, 1], probs = 0.975)
quantile(RES[, 2], probs = 0.025); quantile(RES[, 2], probs = 0.975)
quantile(RES[, 3], probs = 0.025); quantile(RES[, 3], probs = 0.975)
quantile(RES[, 4], probs = 0.025); quantile(RES[, 4], probs = 0.975)
quantile(RES[, 5], probs = 0.025); quantile(RES[, 5], probs = 0.975)
quantile(RES[, 6], probs = 0.025); quantile(RES[, 6], probs = 0.975)
quantile(RES[, 7], probs = 0.025); quantile(RES[, 7], probs = 0.975)
quantile(RES[, 8], probs = 0.025); quantile(RES[, 8], probs = 0.975)
quantile(RES[, 9], probs = 0.025); quantile(RES[, 9], probs = 0.975)


RES <- matrix(NA, nrow = 1000, ncol = 9)
for (b in 1:B)
{
  d_boot <- c()
  clust_sp <- unique(sample(unique(data_intervention$num_centre), replace = T))
  
  for (i in 1:length(clust_sp))
  {
    res <- d_intervention[data_intervention$num_centre == clust_sp[i], ]
    d_boot <- rbind(d_boot, res)
  }
  RES[b, ] <- as.numeric(f.meas(out = femme_nal, g = factor(num_centre), est_p = as.numeric(prop.table(table(d_boot$femme_nal))[2]),
                                dat = d_boot, nsim = 5000))
}

quantile(RES[, 1], probs = 0.025); quantile(RES[, 1], probs = 0.975)
quantile(RES[, 2], probs = 0.025); quantile(RES[, 2], probs = 0.975)
quantile(RES[, 3], probs = 0.025); quantile(RES[, 3], probs = 0.975)
quantile(RES[, 4], probs = 0.025); quantile(RES[, 4], probs = 0.975)
quantile(RES[, 5], probs = 0.025); quantile(RES[, 5], probs = 0.975)
quantile(RES[, 6], probs = 0.025); quantile(RES[, 6], probs = 0.975)
quantile(RES[, 7], probs = 0.025); quantile(RES[, 7], probs = 0.975)
quantile(RES[, 8], probs = 0.025); quantile(RES[, 8], probs = 0.975)
quantile(RES[, 9], probs = 0.025); quantile(RES[, 9], probs = 0.975)


####################            Simulation of correlated continuous outcome and dichotmisation            #################### 


###   Arguments   ###
# k: the number of clusters
# m: the average size of clusters
# v: the variance of cluster sizes
# rho: the intraclass correlation coefficient
# threshold: the threshold for dichotomisation

###   Values   ###
# p_est
# vpc1_est, vpc2_est, vpc3_est, vpc4_est: the estimated variance partition coefficients
# mor_est: the estimated median odds ratio
# tcc_est: the estimated TCC using original formula
# neg_tcc_est: an output egal to 1 if the estimated TCC is less than 0 and egal to 0 otherwise
# rho_est: the estimated ICC
# neg_rho_est: an output egal to 1 if the estimated ICC is less than 0 and egal to 0 otherwise
# tcc_kirk_est: the estimated TCC using Kirk's formula

sim.2 <- function(k, m, v, rho, threshold)
{
  # variable cluster sizes
  clust_siz <- rnbinom(n = k, size = m^2/(v-m), mu = m)
  
  # delete empty clusters
  i <- 0
  while (sum(is.element(el = c(0,1), set = clust_siz)) >= 1 & i < 1000)
  {
    clust_siz <- rnbinom(n = k, size = m^2/(v-m), mu = m)
    i <- i + 1
  }
  
  # sample size
  n <- sum(clust_siz)
  
  ## model: yij = mu + alphai + epsij
  ## mu = 0; alphai ~ rnorm(0,ss), epsij ~ logis(0, 1-ss)
  ## rho = ss/(ss+1-ss) = ss
  
  ss <- rho
  
  # generate alphai, i = 1,...,k
  alpha <- rep(rnorm(n = k, mean = 0, sd = sqrt(ss)), clust_siz)
  
  # generate epsilonij
  epsilon <- rlogis(n = n, location = 0, scale = sqrt(1-ss))
  
  y <- alpha + epsilon
  
  while (min(y) >= min(threshold) | max(y) <= max(threshold))
  {
    # generate alphai, i = 1,...,k
    alpha <- rep(rnorm(n = k, mean = 0, sd = sqrt(ss)), clust_siz)
    
    # generate epsilonij
    epsilon <- rlogis(n = n, location = 0, scale = sqrt(1-ss))
    
    y <- alpha + epsilon
  }
  
  # variable "cluster"
  clust <- factor(rep(1:k, clust_siz))
  
  th <- length(threshold)
  p_est <- c(); vpc1_est <- c(); vpc2_est <- c(); vpc3_est <- c(); vpc4_est <- c(); mor_est <- c()
  rho_est <- c(); tcc_est <- c(); tcc_kirk_est <- c()
  neg_tcc_est <- c(); neg_rho_est <- c()
  
  for (i in 1:th)
  {
    ybin <- factor(ifelse(y > threshold[i], 1, 0))
    
    data_sim <- data.frame(Outcome = ybin, Cluster = clust)
    
    p_est[i]  <- as.numeric(prop.table(table(data_sim$Outcome))[2]) 
    
    res <- f.int(out = Outcome, g = Cluster, dat = data_sim, nsim = 5000)
    
    vpc1_est[i] <- res$vpc1
    vpc2_est[i] <- res$vpc2
    vpc3_est[i] <- res$vpc3
    vpc4_est[i] <- res$vpc4
    mor_est[i]  <- res$mor
    
    # estimation of the TCC using original formula
    tcc_est[i] <- approx_TCC.1(mat = table(make_pairs(data_sim[ , c("Cluster","Outcome")])), delta = 0.00001)$rtet
    
    neg_tcc_est[i] <- ifelse(tcc_est[i] < 0, 1, 0)
    if(tcc_est[i] < 0) {tcc_est[i] <- 0}
    
    # estimation of rho using the Analysis of variance estimator
    rho_est[i] <- est_anova(g = Cluster, out = Outcome, data = data_sim)$Rho
    
    neg_rho_est[i] <- ifelse(rho_est[i]  < 0, 1, 0)
    if(rho_est[i]  < 0) {rho_est[i]  <- 0}
    
    # estimation of the TCC using Kirk's formula
    tcc_kirk_est[i] <- rec_sum_f_int(rhobin = rho_est[i], n = 10, h = threshold[i], p = pnorm(threshold[i]))
    
  }
  
  return(list(p_est = p_est, vpc1_est = vpc1_est, vpc2_est = vpc2_est, vpc3_est = vpc3_est, vpc4_est = vpc4_est, mor_est = mor_est,
              tcc_est = tcc_est, neg_tcc_est = neg_tcc_est, rho_est = rho_est, neg_rho_est = neg_rho_est, tcc_kirk_est = tcc_kirk_est))
  
}


##########                                            Data generation                                            ##########

##########################################     With rho = 0.01     ##########################################

ncores <- 8
nfor <- 10
nsimul <- 1000

Res <- list()
for (j in 1:nfor)
{ 
  registerDoParallel(cores = ncores)
  res <- foreach(i = 1:nsimul, .packages = "lme4", .errorhandling = "remove") %dopar% sim.2(k = 10, m = 25, v = 225, rho = 0.01, 
                                                                                            threshold = seq(from = 5, to = -5, length.out = 20))
  Res  <- c(Res, res)
  stopImplicitCluster()
}

RES <- array(data = NA, dim = c(nsimul*nfor, 20, 11), dimnames = list(NULL, NULL, c("Mean_p", "Mean_vpc1", "Mean_vpc2", "Mean_vpc3", 
                                                                                    "Mean_vpc4", "Mean_mor", "Mean_tcc", "mean_neg_tcc_est", "Mean_rho_est", "mean_neg_rho_est",
                                                                                    "Mean_tcc_kirk_est")))

RES[, , 1] <- t(sapply(Res, function(a){a$p_est}))
RES[, , 2] <- t(sapply(Res, function(a){a$vpc1_est}))
RES[, , 3] <- t(sapply(Res, function(a){a$vpc2_est}))
RES[, , 4] <- t(sapply(Res, function(a){a$vpc3_est}))
RES[, , 5] <- t(sapply(Res, function(a){a$vpc4_est}))
RES[, , 6] <- t(sapply(Res, function(a){a$mor_est}))
RES[, , 7] <- t(sapply(Res, function(a){a$rho_est}))
RES[, , 8] <- t(sapply(Res, function(a){a$tcc_est}))
RES[, , 9] <- t(sapply(Res, function(a){a$tcc_kirk_est}))
RES[, , 10] <- t(sapply(Res, function(a){a$neg_tcc_est}))
RES[, , 11] <- t(sapply(Res, function(a){a$neg_rho_est}))


wb <- createWorkbook()
addWorksheet(wb, "Mean_p")
addWorksheet(wb, "Mean_vpc1")
addWorksheet(wb, "Mean_vpc2")
addWorksheet(wb, "Mean_vpc3")
addWorksheet(wb, "Mean_vpc4")
addWorksheet(wb, "Mean_mor")
addWorksheet(wb, "Mean_tcc")
addWorksheet(wb, "mean_neg_tcc_est")
addWorksheet(wb, "Mean_rho_est")
addWorksheet(wb, "mean_neg_rho_est")
addWorksheet(wb, "Mean_tcc_kirk_est")

writeDataTable(wb, "Mean_p", x = as.data.frame(RES[, , 1]))
writeDataTable(wb, "Mean_vpc1", x = as.data.frame(RES[, , 2]))
writeDataTable(wb, "Mean_vpc2", x = as.data.frame(RES[, , 3]))
writeDataTable(wb, "Mean_vpc3", x = as.data.frame(RES[, , 4]))
writeDataTable(wb, "Mean_vpc4", x = as.data.frame(RES[, , 5]))
writeDataTable(wb, "Mean_mor", x = as.data.frame(RES[, , 6]))
writeDataTable(wb, "Mean_tcc", x = as.data.frame(RES[, , 7]))
writeDataTable(wb, "mean_neg_tcc_est", x = as.data.frame(RES[, , 8]))
writeDataTable(wb, "Mean_rho_est", x = as.data.frame(RES[, , 9]))
writeDataTable(wb, "mean_neg_rho_est", x = as.data.frame(RES[, , 10]))
writeDataTable(wb, "Mean_tcc_kirk_est", x = as.data.frame(RES[, , 11]))

saveWorkbook(wb, file = "L10000and10and25and225and001.xlsx")



##########################################     With rho = 0.05     ##########################################

Res <- list()
for (j in 1:nfor)
{ 
  registerDoParallel(cores = ncores)
  res <- foreach(i = 1:nsimul, .packages = "lme4", .errorhandling = "remove") %dopar% sim.2(k = 10, m = 25, v = 225, rho = 0.05, 
                                                                                            threshold = seq(from = 5, to = -5, length.out = 20))
  Res  <- c(Res, res)
  stopImplicitCluster()
}

RES <- array(data = NA, dim = c(nsimul*nfor, 20, 11), dimnames = list(NULL, NULL, c("Mean_p", "Mean_vpc1", "Mean_vpc2", "Mean_vpc3", 
                                                                                    "Mean_vpc4", "Mean_mor", "Mean_tcc", "mean_neg_tcc_est", "Mean_rho_est", "mean_neg_rho_est",
                                                                                    "Mean_tcc_kirk_est")))

RES[, , 1] <- t(sapply(Res, function(a){a$p_est}))
RES[, , 2] <- t(sapply(Res, function(a){a$vpc1_est}))
RES[, , 3] <- t(sapply(Res, function(a){a$vpc2_est}))
RES[, , 4] <- t(sapply(Res, function(a){a$vpc3_est}))
RES[, , 5] <- t(sapply(Res, function(a){a$vpc4_est}))
RES[, , 6] <- t(sapply(Res, function(a){a$mor_est}))
RES[, , 7] <- t(sapply(Res, function(a){a$rho_est}))
RES[, , 8] <- t(sapply(Res, function(a){a$tcc_est}))
RES[, , 9] <- t(sapply(Res, function(a){a$tcc_kirk_est}))
RES[, , 10] <- t(sapply(Res, function(a){a$neg_tcc_est}))
RES[, , 11] <- t(sapply(Res, function(a){a$neg_rho_est}))


wb <- createWorkbook()
addWorksheet(wb, "Mean_p")
addWorksheet(wb, "Mean_vpc1")
addWorksheet(wb, "Mean_vpc2")
addWorksheet(wb, "Mean_vpc3")
addWorksheet(wb, "Mean_vpc4")
addWorksheet(wb, "Mean_mor")
addWorksheet(wb, "Mean_tcc")
addWorksheet(wb, "mean_neg_tcc_est")
addWorksheet(wb, "Mean_rho_est")
addWorksheet(wb, "mean_neg_rho_est")
addWorksheet(wb, "Mean_tcc_kirk_est")

writeDataTable(wb, "Mean_p", x = as.data.frame(RES[, , 1]))
writeDataTable(wb, "Mean_vpc1", x = as.data.frame(RES[, , 2]))
writeDataTable(wb, "Mean_vpc2", x = as.data.frame(RES[, , 3]))
writeDataTable(wb, "Mean_vpc3", x = as.data.frame(RES[, , 4]))
writeDataTable(wb, "Mean_vpc4", x = as.data.frame(RES[, , 5]))
writeDataTable(wb, "Mean_mor", x = as.data.frame(RES[, , 6]))
writeDataTable(wb, "Mean_tcc", x = as.data.frame(RES[, , 7]))
writeDataTable(wb, "mean_neg_tcc_est", x = as.data.frame(RES[, , 8]))
writeDataTable(wb, "Mean_rho_est", x = as.data.frame(RES[, , 9]))
writeDataTable(wb, "mean_neg_rho_est", x = as.data.frame(RES[, , 10]))
writeDataTable(wb, "Mean_tcc_kirk_est", x = as.data.frame(RES[, , 11]))

saveWorkbook(wb, file = "L10000and10and25and225and005.xlsx")


##########################################     With rho = 0.3     ##########################################

Res <- list()
for (j in 1:nfor)
{ 
  registerDoParallel(cores = ncores)
  res <- foreach(i = 1:nsimul, .packages = "lme4", .errorhandling = "remove") %dopar% sim.2(k = 10, m = 25, v = 225, rho = 0.3, 
                                                                                            threshold = seq(from = 5, to = -5, length.out = 20))
  Res  <- c(Res, res)
  stopImplicitCluster()
}

length(Res)

RES <- array(data = NA, dim = c(nsimul*nfor, 20, 11), dimnames = list(NULL, NULL, c("Mean_p", "Mean_vpc1", "Mean_vpc2", "Mean_vpc3", 
                                                                                    "Mean_vpc4", "Mean_mor", "Mean_tcc", "mean_neg_tcc_est", "Mean_rho_est", "mean_neg_rho_est",
                                                                                    "Mean_tcc_kirk_est")))

RES[, , 1] <- t(sapply(Res, function(a){a$p_est}))
RES[, , 2] <- t(sapply(Res, function(a){a$vpc1_est}))
RES[, , 3] <- t(sapply(Res, function(a){a$vpc2_est}))
RES[, , 4] <- t(sapply(Res, function(a){a$vpc3_est}))
RES[, , 5] <- t(sapply(Res, function(a){a$vpc4_est}))
RES[, , 6] <- t(sapply(Res, function(a){a$mor_est}))
RES[, , 7] <- t(sapply(Res, function(a){a$rho_est}))
RES[, , 8] <- t(sapply(Res, function(a){a$tcc_est}))
RES[, , 9] <- t(sapply(Res, function(a){a$tcc_kirk_est}))
RES[, , 10] <- t(sapply(Res, function(a){a$neg_tcc_est}))
RES[, , 11] <- t(sapply(Res, function(a){a$neg_rho_est}))


wb <- createWorkbook()
addWorksheet(wb, "Mean_p")
addWorksheet(wb, "Mean_vpc1")
addWorksheet(wb, "Mean_vpc2")
addWorksheet(wb, "Mean_vpc3")
addWorksheet(wb, "Mean_vpc4")
addWorksheet(wb, "Mean_mor")
addWorksheet(wb, "Mean_tcc")
addWorksheet(wb, "mean_neg_tcc_est")
addWorksheet(wb, "Mean_rho_est")
addWorksheet(wb, "mean_neg_rho_est")
addWorksheet(wb, "Mean_tcc_kirk_est")

writeDataTable(wb, "Mean_p", x = as.data.frame(RES[, , 1]))
writeDataTable(wb, "Mean_vpc1", x = as.data.frame(RES[, , 2]))
writeDataTable(wb, "Mean_vpc2", x = as.data.frame(RES[, , 3]))
writeDataTable(wb, "Mean_vpc3", x = as.data.frame(RES[, , 4]))
writeDataTable(wb, "Mean_vpc4", x = as.data.frame(RES[, , 5]))
writeDataTable(wb, "Mean_mor", x = as.data.frame(RES[, , 6]))
writeDataTable(wb, "Mean_tcc", x = as.data.frame(RES[, , 7]))
writeDataTable(wb, "mean_neg_tcc_est", x = as.data.frame(RES[, , 8]))
writeDataTable(wb, "Mean_rho_est", x = as.data.frame(RES[, , 9]))
writeDataTable(wb, "mean_neg_rho_est", x = as.data.frame(RES[, , 10]))
writeDataTable(wb, "Mean_tcc_kirk_est", x = as.data.frame(RES[, , 11]))

saveWorkbook(wb, file = "L10000and10and25and225and03.xlsx")


##########                                            Graphic representation                                            ##########

Res <- matrix(NA, 20, 9)
for (i in 1:9)
{ 
  d <- read.xlsx(xlsxFile = "L10000and10and25and225and001.xlsx",
                 sheet = i, startRow = 2, colNames = FALSE)
  Res[ ,i] <- colMeans(d)
}

P_est <- read.xlsx(xlsxFile = "L10000and10and25and225and001.xlsx",
                   sheet = 1, startRow = 2, colNames = FALSE)
Rho_est <- read.xlsx(xlsxFile = "L10000and10and25and225and001.xlsx",
                     sheet = 7, startRow = 2, colNames = FALSE)

Prev <- Res[ , 1]
Rd_est <- c(); Inf_rho_max  <- c(); Prev_est <- c()
for (j in 1:20)
{  
  rho_max <- Max_ICC[which.min(abs(Max_ICC$Prevalence - Prev[j])), ]$Rho
  ind_inf <- which(Rho_est[ ,j] <= rho_max)
  prev_est <- mean(P_est[ ,j][ind_inf])
  rho_inf_est <- mean(Rho_est[ ,j][ind_inf])
  rd_est <- 100*(rho_max - rho_inf_est)/rho_max
  inf_rho_max <- 100*mean(Rho_est[ ,j] <= rho_max)
  
  Rd_est <- c(Rd_est, rd_est)
  Inf_rho_max <- c(Inf_rho_max, inf_rho_max)
  Prev_est <- c(Prev_est, prev_est)
}


d10000and10and25and225and001 <- data.frame(Mean_p = Res[,1], Mean_vpc1 = Res[,2], Mean_vpc2 = Res[,3], Mean_vpc3 = Res[,4], 
                                           Mean_vpc4 =  Res[,5], Mean_mor = Res[,6], Mean_rho_est = Res[,7], Mean_tcc = Res[,8],
                                           Mean_tcc_kirk_est = Res[,9], Mean2_p = Prev_est, Mean_rd_est = Rd_est, Mean_inf_rho_max = Inf_rho_max)


Res <- matrix(NA, 20, 9)
for (i in 1:9)
{ 
  d <- read.xlsx(xlsxFile = "L10000and10and25and225and005.xlsx",
                 sheet = i, startRow = 2, colNames = FALSE)
  Res[ ,i] <- colMeans(d)
}


P_est <- read.xlsx(xlsxFile = "L10000and10and25and225and005.xlsx",
                   sheet = 1, startRow = 2, colNames = FALSE)
Rho_est <- read.xlsx(xlsxFile = "L10000and10and25and225and005.xlsx",
                     sheet = 7, startRow = 2, colNames = FALSE)

Prev <- Res[ , 1]
Rd_est <- c(); Inf_rho_max  <- c(); Prev_est <- c()
for (j in 1:20)
{  
  rho_max <- Max_ICC[which.min(abs(Max_ICC$Prevalence - Prev[j])), ]$Rho
  ind_inf <- which(Rho_est[ ,j] <= rho_max)
  prev_est <- mean(P_est[ ,j][ind_inf])
  rho_inf_est <- mean(Rho_est[ ,j][ind_inf])
  rd_est <- 100*(rho_max - rho_inf_est)/rho_max
  inf_rho_max <- 100*mean(Rho_est[ ,j] <= rho_max)
  
  Rd_est <- c(Rd_est, rd_est)
  Inf_rho_max <- c(Inf_rho_max, inf_rho_max)
  Prev_est <- c(Prev_est, prev_est)
}


d10000and10and25and225and005 <- data.frame(Mean_p = Res[,1], Mean_vpc1 = Res[,2], Mean_vpc2 = Res[,3], Mean_vpc3 = Res[,4], 
                                           Mean_vpc4 =  Res[,5], Mean_mor = Res[,6], Mean_rho_est = Res[,7], Mean_tcc = Res[,8],
                                           Mean_tcc_kirk_est = Res[,9], Mean2_p = Prev_est, Mean_rd_est = Rd_est, Mean_inf_rho_max = Inf_rho_max)


Res <- matrix(NA, 20, 9)
for (i in 1:9)
{ 
  d <- read.xlsx(xlsxFile = "L10000and10and25and225and03.xlsx",
                 sheet = i, startRow = 2, colNames = FALSE)
  Res[ ,i] <- colMeans(d)
}
P_est <- read.xlsx(xlsxFile = "L0000and10and25and225and03.xlsx",
                   sheet = 1, startRow = 2, colNames = FALSE)
Rho_est <- read.xlsx(xlsxFile = "L10000and10and25and225and03.xlsx",
                     sheet = 7, startRow = 2, colNames = FALSE)

Prev <- Res[ , 1]
Rd_est <- c(); Inf_rho_max  <- c(); Prev_est <- c()
for (j in 1:20)
{  
  rho_max <- Max_ICC[which.min(abs(Max_ICC$Prevalence - Prev[j])), ]$Rho
  ind_inf <- which(Rho_est[ ,j] <= rho_max)
  prev_est <- mean(P_est[ ,j][ind_inf])
  rho_inf_est <- mean(Rho_est[ ,j][ind_inf])
  rd_est <- 100*(rho_max - rho_inf_est)/rho_max
  inf_rho_max <- 100*mean(Rho_est[ ,j] <= rho_max)
  
  Rd_est <- c(Rd_est, rd_est)
  Inf_rho_max <- c(Inf_rho_max, inf_rho_max)
  Prev_est <- c(Prev_est, prev_est)
}


d10000and10and25and225and03 <- data.frame(Mean_p = Res[,1], Mean_vpc1 = Res[,2], Mean_vpc2 = Res[,3], Mean_vpc3 = Res[,4], 
                                          Mean_vpc4 =  Res[,5], Mean_mor = Res[,6], Mean_rho_est = Res[,7], Mean_tcc = Res[,8],
                                          Mean_tcc_kirk_est = Res[,9], Mean2_p = Prev_est, Mean_rd_est = Rd_est, Mean_inf_rho_max = Inf_rho_max)


##########     Graphic representation for VPCs     ##########

p1 <- ggplot(data = d10000and10and25and225and001) +
  geom_line(aes(x = Mean_p, y = Mean_vpc1, col = "Estimated VPC1"), size = 1) +
  geom_line(aes(x = Mean_p, y = Mean_vpc2, col = "Estimated VPC2"), size = 1) +
  geom_line(aes(x = Mean_p, y = Mean_vpc3, col = "Estimated VPC3"), size = 1) + 
  geom_line(aes(x = Mean_p, y = Mean_vpc4, col = "Estimated VPC4"), size = 1) +
  xlab("Prevalence") + ylab("VPC") + scale_colour_manual(name = "", values = c("#00BA38", "#9683EC", "#619CFF", "#F8766D"),
                                                         breaks = c("Estimated VPC1","Estimated VPC2","Estimated VPC3", "Estimated VPC4"), 
                                                         labels = c(expression(VPC[1]), expression(VPC[2]), expression(VPC[3]), expression(VPC[4]))) 

p2 <- ggplot(data = d10000and10and25and225and005) +
  geom_line(aes(x = Mean_p, y = Mean_vpc1), size = 1, col = "#00BA38") +
  geom_line(aes(x = Mean_p, y = Mean_vpc2), size = 1, col = "#9683EC") +
  geom_line(aes(x = Mean_p, y = Mean_vpc3), size = 1, col = "#619CFF") + 
  geom_line(aes(x = Mean_p, y = Mean_vpc4), size = 1, col = "#F8766D") 

p3 <- ggplot(data = d10000and10and25and225and03) +
  geom_line(aes(x = Mean_p, y = Mean_vpc1), size = 1, col = "#00BA38") +
  geom_line(aes(x = Mean_p, y = Mean_vpc2), size = 1, col = "#9683EC") +
  geom_line(aes(x = Mean_p, y = Mean_vpc3), size = 1, col = "#619CFF") + 
  geom_line(aes(x = Mean_p, y = Mean_vpc4), size = 1, col = "#F8766D") 

p10 <- p1 + theme(legend.position = "none")

p11 <- get_legend(p1 + theme(legend.position=c(0.60,0.6), legend.justification = "center"))

plot_grid(p11, NULL, NULL, p10, NULL, p2, NULL, p3, ncol = 2, nrow = 4, 
          rel_widths = c(1, 4), rel_heights = c(2.1, 4, 4, 4), 
          labels = c("",  "k = 10", "(A)", "", "(B)", "", "(C)", ""), 
          label_size = 10, label_x = 0.55,  label_y = 0.65, hjust = 0.5)


##########     Graphic representation for MOR     ##########

p1 <- ggplot(data = d10000and10and25and225and001) + geom_line(aes(x = Mean_p, y = Mean_mor), size = 1) +
  xlab("Prevalence") + ylab("Median odds ratio") + scale_x_continuous(breaks=c(0.1, 0.3, 0.5, 0.7, 0.9)) 

p2 <- ggplot(data = d10000and10and25and225and005) + geom_line(aes(x = Mean_p, y = Mean_mor), size = 1) +
  xlab("Prevalence") + ylab("Median odds ratio") + scale_x_continuous(breaks=c(0.1, 0.3, 0.5, 0.7, 0.9)) 

p3 <- ggplot(data = d10000and10and25and225and03) + geom_line(aes(x = Mean_p, y = Mean_mor), size = 1) +
  xlab("Prevalence") + ylab("Median odds ratio") + scale_x_continuous(breaks=c(0.1, 0.3, 0.5, 0.7, 0.9)) 

plot_grid(NULL, NULL, NULL, p1, NULL, p2, NULL, p3, ncol = 2, nrow = 4, 
          rel_widths = c(1, 4), rel_heights = c(1, 4, 4, 4), 
          labels = c("",  "k = 10", "(A)", "", "(B)", "", "(C)", ""), 
          label_size = 10, label_x = 0.55,  label_y = 0.65, hjust = 0.5)

##########     Graphic representation for TCCs     ##########

p1 <- ggplot(data = d10000and10and25and225and001) +
  geom_line(aes(x = Mean_p, y = Mean_rho_est, col = "Estimated binary ICC"), size = 1) +
  geom_line(aes(x = Mean_p, y = Mean_tcc_kirk_est, col = "Estimated TCC (Kirk formula)"), size = 1) +
  geom_line(aes(x = Mean_p, y = Mean_tcc, col = "Estimated TCC (Original formula)"), size = 1) +
  geom_hline(aes(yintercept = 0.01, col = "Theoretical continuous ICC"), size = 1) +
  xlab("Prevalence") + ylab("") + scale_colour_manual(name="", values=c("#00BA38", "#9683EC", "#619CFF", "#F8766D")) 

p2 <- ggplot(data = d10000and10and25and225and005) + 
  geom_line(aes(x = Mean_p, y = Mean_rho_est), size = 1, col = "#00BA38") +
  geom_line(aes(x = Mean_p, y = Mean_tcc_kirk_est), size = 1, col = "#9683EC") +
  geom_line(aes(x = Mean_p, y = Mean_tcc), size = 1, col = "#619CFF") +
  geom_hline(aes(yintercept = 0.05), size = 1, col = "#F8766D") + xlab("Prevalence") + ylab("") 

p3 <- ggplot(data = d10000and10and25and225and03) + 
  geom_line(aes(x = Mean_p, y = Mean_rho_est), size = 1, col = "#00BA38") +
  geom_line(aes(x = Mean_p, y = Mean_tcc_kirk_est), size = 1, col = "#9683EC") +
  geom_line(aes(x = Mean_p, y = Mean_tcc), size = 1, col = "#619CFF") +
  geom_hline(aes(yintercept = 0.3), size = 1, col = "#F8766D") + xlab("Prevalence") + ylab("") 
  
p10 <- p1 + theme(legend.position = "none")

p11 <- get_legend(p1 + theme(legend.position=c(0.98, 0.6), legend.justification = "center"))

plot_grid(p11, NULL, NULL, p10, NULL, p2, NULL, p3, ncol = 2, nrow = 4, 
          rel_widths = c(2, 4), rel_heights = c(2, 4, 4, 4), 
          labels = c("",  "k = 10", "(A)", "", "(B)", "", "(C)", ""), 
          label_size = 10, label_x = 0.55,  label_y = 0.65, hjust = 0.5)


##########     Graphic representation for the relative deviation to the theoretical maximum ICC value     ##########

p1 <- ggplot(data = d10000and10and25and225and001) + geom_line(aes(x = Mean2_p, y = Mean_rd_est), size = 1) + 
  geom_point(aes(x = Mean2_p, y = Mean_inf_rho_max), col = "red", pch = 15, size = 1.5)  +
  xlab("Prevalence") + ylab("Distance to the \n max ICC value (%)") 

p2 <- ggplot(data = d10000and10and25and225and005)  + geom_line(aes(x = Mean2_p, y = Mean_rd_est), size = 1) + 
  geom_point(aes(x = Mean2_p, y = Mean_inf_rho_max), col = "red", pch = 15, size = 1.5) + 
  xlab("Prevalence") + ylab("Distance to the \n max ICC value (%)") 

p3 <- ggplot(data = d10000and10and25and225and03)  + geom_line(aes(x = Mean2_p, y = Mean_rd_est), size = 1) + 
  geom_point(aes(x = Mean2_p, y = Mean_inf_rho_max), col = "red", pch = 15, size = 1.5) +
  xlab("Prevalence") + ylab("Distance to the \n max ICC value (%)")  

plot_grid(NULL, NULL, NULL, p1, NULL, p2, NULL, p3, ncol = 2, nrow = 4, 
          rel_widths = c(0.4, 2), rel_heights = c(0.5, 2, 2, 2), 
          labels = c("",  "k = 10", "(A)", "", "(B)", "", "(C)", ""), 
          label_size = 10, label_x = 0.55,  label_y = 0.65, hjust = 0.5)


##########                               Maximal Information Coefficient calculation                               ##########

RES <- matrix(NA, nrow = 3, ncol = 8)
RES[1, ] <- unlist(f_mic(dat = d10000and10and25and225and001))
RES[2, ] <- unlist(f_mic(dat = d10000and10and25and225and005))
RES[3, ] <- unlist(f_mic(dat = d10000and10and25and225and03))
RES
