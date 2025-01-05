# Option Greeks at different levels of Moneyness and Maturity

How do option values and greeks evolve as the contract approaches maturity?
The Black-Scholes option pricing framework defines the relationship between implied volatility of the underlying, the time remaining to expiry, and the value of the option.

Using an option pricing formula and numerical finite differencing to calculate partial derivatives, it is possible to compute the greeks and study how they evolve as time progresses.

In this project, I visualize payoffs and greeks for a single European option, using the ggplot2 package.