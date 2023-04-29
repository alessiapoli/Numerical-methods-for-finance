
########################## Cleaning the dataset ################################
ds_cleaningMC <- function(dataset) {
  
  # Definition of option to be considered either call or put. 
  cond_type_call <- (dataset$option_type == "C" | dataset$option_type == "c" |
                       dataset$option_type == "Call" | 
                       dataset$option_type == "call")
  
  cond_type_put  <- (dataset$option_type == "P" | dataset$option_type == "p" |
                       dataset$option_type == "Put" | 
                       dataset$option_type == "put")
  
  cond_type <- cond_type_call | cond_type_put
  stopifnot("The type of option is wrong/missing" = cond_type == TRUE)
 
  splitting_call <- cond_type_call
  ds_callsplit <- dataset[splitting_call, ]
  
  # Definition of Merton's contraints for a call option
  cond_call <- ds_callsplit$underlying_price >= ds_callsplit$option_price   &
    ds_callsplit$option_price >= pmax(ds_callsplit$underlying_price   -
                                        ds_callsplit$strike *
                                        exp(-ds_callsplit$interest_rate * 
                                              (ds_callsplit$ttm/252)), 0)
  
  ds_call <- ds_callsplit[cond_call,] # Call options that satisfy the 
                                      # Merton's contraints
  
  nec <- abs(dim(ds_call)[1] - dim(ds_callsplit)[1])
  cat("Number of eliminated calls:", nec, fill = TRUE)
  
  percc <- (nrow(set_opt) - nec) / nrow(set_opt) * 100
  cat("Percentage of the call options in the cleaned dataset:", 
      percc, fill = TRUE, "%")
  
  
  splitting_put <- cond_type_put
  ds_putsplit <- dataset[splitting_put, ]
  
  # Definition of Merton's contraints for a put option
  cond_put <- ds_putsplit$strike *
    exp(-ds_putsplit$interest_rate *
          (ds_putsplit$ttm/252)) >= ds_putsplit$option_price &
    ds_putsplit$option_price >= pmax(ds_putsplit$strike *
                                       exp(-ds_putsplit$interest_rate * 
                                             (ds_putsplit$ttm/252)) -
                                       ds_putsplit$underlying_price , 0)
  
  ds_put <- ds_putsplit[cond_put,] # Put options that satisfy the 
                                   # Merton's contraints
  
  nep <- abs(dim(ds_put)[1] - dim(ds_putsplit)[1]) 
  cat("Number of eliminated puts:", nep, fill = TRUE)
  
  percp <- (nrow(set_opt) - nep) / nrow(set_opt) * 100
  cat("Percentage of the put options in the cleaned dataset", percp, 
      fill = TRUE, "%")
  
  return(list(call_set = ds_call, put_set = ds_put, 
              merge_set = rbind(ds_call, ds_put)))
}



############################ Volatility smile ##################################
option_priceBS <- function(S0, K, r, ttm, sigma, TypeOption){
  
  # Arguments of the Gaussian CDF in the Black & Scholes call option price
  d1 <- ( log(S0/K) + (r + 0.5 * sigma^2) * ttm)/(sigma * sqrt(ttm))
  d2 <- d1 - sigma * sqrt(ttm)
  
  # Black & Scholes call option prices
  if(TypeOption == "C" | TypeOption == "Call" | TypeOption == "call"){
    opt_price <- S0 * pnorm(d1) - K * exp(-r * ttm) * pnorm(d2)
  
    # Black & Scholes put option prices  
  }else if (TypeOption == "P" | TypeOption == "Put" | TypeOption == "put"){
    opt_price <- K * exp(-r * ttm) * pnorm(-d2) - S0 * pnorm(-d1)
    
  }else{
    stop("Type of the option is not defined")
    
  }
  return(opt_price)
}

# Formula to find the implied volatility using Black & Scholes prices
impliedVola_optimFUN <- function(S0, K, r, ttm, par, TypeOption, mkt_price){

    res <- (option_priceBS(S0 = S0, K = K, ttm = ttm, r = r, sigma = par,
                         TypeOption = TypeOption) - mkt_price)^2

    return(res)
}



################################### Strategy ###################################
# Function to calculate the call options Black & Scholes prices
call_priceBS <- function(S0, K, r, ttm, sigma){
  
  d1 <- ( log(S0/K) + (r + 0.5 * sigma^2) * ttm)/(sigma * sqrt(ttm))
  d2 <- d1 - sigma * sqrt(ttm)
  
  call_price <- S0 * pnorm(d1) - K * exp(-r * ttm) * pnorm(d2)
  
  return(call_price)
}
  
# Function to find the strategy profit, which also returns the plot
strategy <- function(S_T, K1, K2, premium1, premium2, plt = TRUE) {
  
  res_call <- S_T - K1 - premium1 
  res_put <- K2 - S_T - premium2
  
  # Profit of a long position in a call with strike price K1 
  profit_long_call <- pmax(-premium1, res_call) 
  # Profit of a long position in a put with strike price K2
  profit_long_put <- pmax(-premium2, res_put)
  # Profit of the strangle strategy 
  strategy_profit <- rowSums(cbind(profit_long_call, profit_long_put))
  
  dt <- data.frame(cbind(S_T, profit_long_call, profit_long_put, 
                         strategy_profit))
  
  # Construction of the plot
  if(plt){
    p <- ggplot(dt, aes(x = S_T)) + theme_light() +
      geom_line(aes(y = profit_long_call, color = "Long Call",
                    linetype = "Long Call"), lwd= 0.8) +
      geom_line(aes(y = profit_long_put, color="Long Put",
                    linetype = "Long Put"), lwd= 0.8)+
      geom_line(aes(y= strategy_profit, color = "Strategy",
                    linetype = "Strategy"), lwd= 0.8) +
      xlab(TeX("$S_T$")) + ylab("Profit") +geom_hline(yintercept= 0, 
                                  color = "black", lwd= 0.8 ) +
      scale_linetype_manual(name="Legend",
                            breaks = c("Long Call", "Long Put", "Strategy"),
                            values= c("dashed", "dashed", "solid")) +
      scale_color_manual(name="Legend",
                         breaks = c("Long Call", "Long Put", "Strategy"),
                         values = c("red", "blue", "green")) +
      theme(legend.position="top")
    print(p)
  }
  
  return(strategy_profit)
}



########################## Calibration function ################################
# Function to calibrate the volatility using as distance the relative mean 
# squared error (RMSE) 
calibrationRMSE_optimFUN <- function(S0, K, r, ttm, par, TypeOption, mkt_price){
  
  error <- ((option_priceBS(S0 = S0, K = K, ttm = ttm, r = r, sigma = par,
                            TypeOption = TypeOption) - mkt_price)/mkt_price)^2
  
  return(mean(error))
}



######################### Monte Carlo simulation ###############################
# GBM simulation of the underlying prices at different times using the for-loop
simGBM_forloop <- function(S0, ttm, sigma, r, nsim, N=N){ 
  
  S0_initial <- rep(S0, nsim)
  
  deltat <- (ttm)/N
  
  path_gbm <- matrix(NA, nsim, N+1) 
  
  path_gbm[, 1] <- S0_initial
  
  for(t in c(2:(N + 1))){
    # path-gbm -> underlying prices
    path_gbm[, t] <- path_gbm[, t-1] * exp((r - 0.5*sigma^2)*(deltat) +
                                             sigma*sqrt(deltat)*rnorm(nsim))
  }
  
  time_grid <- seq(0 , ttm, by = deltat)
  
  return(list(grid = time_grid, paths = path_gbm))
}



# Function to define the Black & Scholes call option prices using 
# Monte-Carlo simulation
BS_call_MC <- function(S0, sigma, ttm, r, nsim){
  
  sim_us <- simGBM_forloop(S0 = S0, ttm = ttm,sigma = sigma, r = r, 
                           nsim = nsim, N = N)$paths
  
  # Definition of the final payoff
  ST <- sim_us[, ncol(sim_us)] 
  S_prod <- apply(sim_us, 1, function(x) exp(sum(log(x)))) #or simply prod(x)
  sim_final_payoff <- pmax(ST-(S_prod)^(1/(N+1)), 0) 
  
  # Construction of the confidence interval given alpha = 95%
  discounted_final_payoff <- exp(-r*(ttm))*sim_final_payoff # call prices
  variance_ofMC_foraCall <- var(discounted_final_payoff)
  term1 <- 1.96*sqrt(variance_ofMC_foraCall)/sqrt(nsim) 
  
  MC_price <- (1/nsim)*sum(discounted_final_payoff) #mean of all the call prices
  
  return(list(MCprice = MC_price, LB= MC_price-term1, UB = MC_price+term1))
}

