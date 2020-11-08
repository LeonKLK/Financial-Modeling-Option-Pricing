#Introduction

This is a final year projet in financial modeling under the Physics Department in The Chinese University of Hong Kong(CUHK).

This projet is about option pricing using numerical method(explicit, implicit, Crank–Nicolson, Monte Carlo), for vanilla option and spread option.

The output contains the option price, Delta and Gamma for different trading strategies. Notice that this code isn't practical in real market since we have different assumptions, for example different parameters are fixed and we asuumped the market is efficient(stock price follows a complete random walk). 

#Remarks

For the explicit, implicit and Crank Nicolson finite difference method, notice that for the interested range, the error could be huge at the begining and end, where normally won't be reached in the market(5 s.d.).

For Monte Carlo method, this method is the most accurate but the most time consuming, thus one can always choose the desire range, for his best use of time. Notice that for the random number from a Guassian distribution, we use the Box-Muller method.
