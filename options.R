###################################################################################################
# Black-Scholes-Metron formula for price of a European option.
###################################################################################################
# Inputs :-
# S is price of the underlying stock
# K is strike price of the option
# tau is time remaining to maturity of the option, as a year fraction
# vol is volatility, Black-Scholes-Metron assumes flat vol across strikes
# r is the risk-free rate, defaults to 4.6% (long-term average Fed rate)
# d is the annual dividend rate, defaults to 0
# is_call indicates whether the option is a call or a put (use is_call = FALSE for put option)
###################################################################################################
# Outputs :-
# price is the Black-Scholes price of the option 
###################################################################################################

opt_bs_price <- function(S, K, tau, vol, r=0.046, d=0, is_call=TRUE) {
  if (tau == 0) {
    # Case when option is at maturity
    if (is_call) {
      return(pmax(S - K, 0))
    } else {
      return(pmax(K - S, 0))
    }
  }
  d1 <- (log(S/K) + (r-d + 0.5*vol^2)*tau) / (vol*sqrt(tau))
  d2 <-  d1 - vol*sqrt(tau)
  if(is_call) {
    price <- S*exp(-d*tau)*pnorm(d1) - K*exp(-r*tau)*pnorm(d2)
    }
  else {
    price <- -S*exp(-d*tau)*pnorm(-d1) + K*exp(-r*tau)*pnorm(-d2)
    }
  return(price)
  }


###################################################################################################
# Delta of an option using finite forward difference on Black-Scholes formula
###################################################################################################
# Additional Inputs :-
# dS is the size of forward step, defaults to 1 cent (or equivalent of 1/100 units of currency)
###################################################################################################
# Outputs :-
# delta is the option price's sensitivity to S
###################################################################################################

opt_get_delta <- function(S, K, tau, vol, r=0.046, d=0, is_call=TRUE, dS=0.01) {
  S_plus <- S + dS
  price <- opt_bs_price(S, K, tau, vol, r, d, is_call)
  price_plus <- opt_bs_price(S_plus, K, tau, vol, r, d, is_call)
  delta <- (price_plus - price) / dS
  return(delta)
}


###################################################################################################
# Gamma of an option using finite forward difference on Black-Scholes formula
###################################################################################################
# Additional Inputs :-
# dS is the size of forward step, defaults to 0.5 cent (or equivalent of 1/200 units of currency)
###################################################################################################
# Outputs :-
# gamma is the option delta's sensitivity to S
###################################################################################################

opt_get_gamma <- function(S, K, tau, vol, r=0.046, d=0, is_call=TRUE, dS=0.005) {
  S_plus <- S + dS
  S_plusplus <- S + dS + dS
  price <- opt_bs_price(S, K, tau, vol, r, d, is_call)
  price_plus <- opt_bs_price(S_plus, K, tau, vol, r, d, is_call)
  price_plusplus <- opt_bs_price(S_plusplus, K, tau, vol, r, d, is_call)
  gamma <- (price_plusplus - 2*price_plus + price) / dS^2
  return(gamma)
}


###################################################################################################
# Vega of an option using finite forward difference on Black-Scholes formula
###################################################################################################
# Additional Inputs :-
# dvol is the size of forward step, defaults to 0.001 (0.1% volatility)
###################################################################################################
# Outputs :-
# vega is the option price's sensitivity to vol
###################################################################################################

opt_get_vega <- function(S, K, tau, vol, r=0.046, d=0, is_call=TRUE, dvol=0.001) {
  vol_plus <- vol + dvol
  price <- opt_bs_price(S, K, tau, vol, r, d, is_call)
  price_plus <- opt_bs_price(S, K, tau, vol_plus, r, d, is_call)
  vega <- (price_plus - price) / dvol
  return(vega)
}


###################################################################################################
# Theta of an option using finite forward difference on Black-Scholes formula
###################################################################################################
# Additional Inputs :-
# dtau is the size of forward step, defaults to 1/365 (1 day)
###################################################################################################
# Outputs :-
# theta is the option price's sensitivity to the passage of time
###################################################################################################

opt_get_theta <- function(S, K, tau, vol, r=0.046, d=0, is_call=TRUE, dtau=1/365) {
  tau_plus <- tau + dtau
  price <- opt_bs_price(S, K, tau, vol, r, d, is_call)
  price_yesterday <- opt_bs_price(S, K, tau_plus, vol, r, d, is_call)
  theta <- (price - price_yesterday) / dtau
  return(theta)
}

