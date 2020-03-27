# portfolio_sharpe_ratio

*Info:
This Matlab package optimizes portfolio weights based on Sharpe ratio,
average total return, and volatility (standard deviation of returns). The
basic idea is to provide a directory as input 'csv_dir'. This folder
should contain .csv files of historical data for each ticker desired to 
comprise part of a portfolio. The function will then return an optimized 
weighting scheme based on high Sharpe ratio, high total return, and low 
volatility (you can provide a custom weighting for each criterion). It 
will also output historical performance data for the portfolio, an 
alternative equally-weighted portfolio, and each individual ticker. 
Furthermore, you can choose to plot a matrix showing correlations among 
the individual assets, as well as a 3D scatter plot showing where the 
optimized portfolio falls among other randomly generated portfolios on 
the dimensions of Sharpe ratio, mean annual return, and volatility (SD of
returns). The function was designed using data from Yahoo Finance
(https://finance.yahoo.com/) but should work with other data sources
provided the formatting is similar. Try optimizing the sample portfolio 
included (with data from 1998-2018) to gain a better understanding of use
cases.

*Cautionary note: your results may be skewed if your data does not go back
sufficiently far and includes only one portion of a market cycle (e.g.,
all bull market, no recessions). To gain insight into performance over
the whole market cycle, try to include data going back to 2008 or 2000, 
if not longer. When this is not possible, note that assets with high
volatility and high annual returns during a bull market will often be the
same assets that sustain the largest losses during economic downturns.
In the case of major recessions these losses can sometimes exceed -40% 
in a calendar year.

*Disclaimer: This open-source research tool is not intended to provide 
  investment advice. It is intended only for informational purposes, and 
  the user is not recommended to use the tool to make actual investment 
  decisions. Seek a duly licensed professional for investment advice.
  
##
If you find portfolio_sharpe_ratio useful and would like to support its continued development, feel free to send a cup of coffee! :) <br><br>
[![paypal](https://www.paypalobjects.com/en_US/i/btn/btn_donateCC_LG.gif)](https://paypal.me/ElliotLayden?locale.x=en_US)
