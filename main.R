setwd("") # \set working directory 
source("functions.R")

#Load needed libraries
library(latex2exp)
library(tidyverse)
library(readxl)
library(dplyr)
library(timeDate)
library(ggplot2)



############################### Load the dataset ############################### 
set_opt <- read_excel("Datasets.xlsx", 
                     col_types= c("text", "date", "date", "numeric", "text", 
                                 "numeric", "numeric", "numeric", 
                                 "numeric", "numeric"))

# Characteristics of the option set
names(set_opt) # Functions to get the column names 
str(set_opt) # Compactly display the internal structure 
class(set_opt) # Class 
dim(set_opt) # Retrieve the dimension 
glimpse(set_opt) # Get a glimpse 



############################### Clean the dataset ##############################
# Aplly the ds_cleaningMC function to clean the dataset, by removing the 
# options that do not satisfy the Merton's contraints 
clean1 <- ds_cleaningMC(dataset = set_opt)  


############################### Volatility smile ###############################
# Extract the merge_set, which is composed of the call_set and the put_set 
# (recall we only have calls in this dataset)
clean_set <- clean1$merge_set 

# Values extracted from clean_set
S0  <- clean_set$underlying_price # Underlying price 
ttm <- clean_set$ttm/252 # Time to maturity
K   <- clean_set$strike # Strike price
r   <- clean_set$interest_rate # Interest rate
opt_price  <- clean_set$option_price # Option price
TypeOption <- clean_set$option_type # Option type

ImpliedVol <- numeric(length = length(opt_price)) 
set.seed(123)
# Use of a for-cycle to minimize the implied volatility for each option with 
# optim, given initial guess of 0.2 and method L-BFGS-B.
for(i in c(1:length(opt_price))){
  dummy <- optim(par = 0.2, fn = impliedVola_optimFUN, lower = 0.001,
                 upper = 1, method ="L-BFGS-B", S0 = S0[i], K = K[i],
                 ttm = ttm[i], r = r[i], TypeOption = TypeOption[i],
                 mkt_price = opt_price[i]) 
  # the if else statement tells if there is a successful completion 
  # (convergence = 0). 
  if(dummy$convergence==0){ cat("\n", "okay")
  }else{
    cat("\n objective function ", c(i, dummy$convergence))
  }
  ImpliedVol[i] <- dummy$par
} 

# Add to the clean set the column with implied volatilities
set_wIV <- clean_set %>% mutate(imp_vol = ImpliedVol) 

# As all the calls have same expiration date, plot the volatility smile, 
# graph of the implied volatilities and strike prices 
plot(set_wIV$strike, set_wIV$imp_vol, ylab = "IV", xlab = "K") 



################################### Strategy ###################################
S0 <- set_wIV$underlying_price[1] # Extract the initial price
K1 <- 1.05*S0 # Define K1
K2 <- 0.95*S0 # Define K2

# Interpolation to get implied volatility for the strike wanted
interpolation <- approx(x = set_wIV$strike, y = set_wIV$imp_vol, xout = K1) 

sigma <- interpolation$y[1]
sigma 

ttm=21/252 # Same maturity as options in dataset
r <- set_wIV$interest_rate[1] # Extract the interest rate from the dataset

price_call1 <- call_priceBS(S0=S0, K=K1, r=r, ttm=ttm, sigma=sigma) 
price_call1 # First call price
price_call2 <- call_priceBS(S0=S0, K=K2, r=r, ttm=ttm, sigma=sigma) 
price_call2 # Second call price

price_put <- price_call2 - S0 + K2*exp(-r*ttm) # Apply the put-call parity
# to the call with strike price K2 in order to find the price of the put
# option with strike price K2

price_strategy <- price_call1 + price_put 
price_strategy # Strategy price = sum of the prices of the 2 options

premium1<- price_call1 # Premium paid for the call option at time t0
premium2 <- price_put # Premium paid for the put option at time t0

S_T <- seq(45, 80, 0.01)

# Apply the function to find the strategy profit and see its plot 
strategy_profit <- strategy(S_T = S_T, K1=K1, K2=K2, premium1=premium1,
                            premium2=premium2, plt = TRUE)

# Negative profit when K2 - premium2 - premium1 < S_T < K1 + premium1 + premium2 
# ==> when 58.1 < S_T < 66.66
b1 <- (log((K1+premium1+premium2)/S0) - (r-0.5*sigma^2)*(ttm))/(sigma*sqrt(ttm))
b2 <- (log((K2-premium2-premium1)/S0) - (r-0.5*sigma^2)*(ttm))/(sigma*sqrt(ttm))
prob_loss <- pnorm(b1)-pnorm(b2) 
prob_loss



########################### Calibration of the model ###########################
# Calibration of the volatility using as distance the relative mean squared 
# error (RMSE) 
res_RMSEcal <- optim(par = 0.19, fn = calibrationRMSE_optimFUN, lower = 0.01,
                     method ="L-BFGS-B",
                     S0 = S0, K = K,
                     ttm = ttm,
                     r = r, mkt_price = opt_price, TypeOption = 
                       unique(TypeOption))

res_RMSEcal$par



############################ Monte Carlo simulation ############################
S0 <- set_wIV$underlying_price[1] # Extract the underlying price
sigma <- res_RMSEcal$par #  Extract the calibrated volatility 
r <- set_wIV$interest_rate[1] # Extract the interest rate
ttm <- 30/252 # Time to maturity
N<-6 # Delta = ttm/N ==> N = ttm/delta = (30/252)/(5/252) = 6

# Black & Scholes call option price using Monte Carlo simulation for different 
# number of simulations
set.seed(123)
BS_call_MC10 <- BS_call_MC(S0 = S0, sigma=sigma, ttm=ttm, r=r, nsim=10)

set.seed(123)
BS_call_MC50 <- BS_call_MC(S0 = S0, sigma=sigma, ttm=ttm, r=r, nsim=50)

set.seed(123)
BS_call_MC100 <- BS_call_MC(S0 = S0, sigma=sigma, ttm=ttm, r=r, nsim=100)

set.seed(123)
BS_call_MC150 <- BS_call_MC(S0 = S0, sigma=sigma, ttm=ttm, r=r, nsim=150)

set.seed(123)
BS_call_MC200 <- BS_call_MC(S0 = S0, sigma=sigma, ttm=ttm, r=r, nsim=200)

set.seed(123)
BS_call_MC250 <- BS_call_MC(S0 = S0, sigma=sigma, ttm=ttm, r=r, nsim=250)

set.seed(123)
BS_call_MC500 <- BS_call_MC(S0 = S0, sigma=sigma, ttm=ttm, r=r, nsim=500)

set.seed(123)
BS_call_MC1000 <- BS_call_MC(S0 = S0, sigma=sigma, ttm=ttm, r=r, nsim=1000)

set.seed(123)
BS_call_MC100000 <- BS_call_MC(S0 = S0, sigma=sigma, ttm=ttm, r=r, nsim=100000)

set.seed(999)
m <- c(10, 50, 100, 150, 200, 250, 500, 1000) #Vector of numbers of simulations 
p1 <- NULL 
err <- NULL
nM <- length(m)
repl <- 100
mat <- matrix(NA, repl, nM)

# Construction of the box-plot graph associated to the number of simulations 
# in the vector m
for (k in 1:nM) {
  tmp <- numeric(repl)
  for (i in 1:repl) tmp[i] <- BS_call_MC(S0 = S0,sigma = sigma, r = r, ttm=ttm, 
                                         nsim = m[k])$MCprice
  mat[, k] <- tmp
  p1 <- c(p1, mean(tmp))
  err <- c(err, sd(tmp))
}
colnames(mat) <- m

minP <- min(p1 - err)
maxP <- max(p1 + err)
plot(m, p1, type = "n", ylim = c(minP, maxP), axes = F, 
     ylab = "MC price", 
     xlab = "MC simulations")
lines(m, p1 + err, col = "blue")
lines(m, p1 - err, col = "blue")
axis(2)
axis(1, m)
boxplot(mat, add = TRUE, at = m, boxwex = 15, col = "orange",axes = F)
points(m, p1, col = "blue", lwd = 3, lty = 3)


