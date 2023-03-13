\name{RolWinWavCor-package}
\alias{RolWinWavCor-package}
\alias{RolWinWavCor}
\docType{package}
\title{Estimate rolling window wavelet correlation between two time series
}
\description{
'RolWinWavCor' estimates and plots as a heat map rolling window wavelet correlation (RWWC) coefficients statistically significant (within the 95\% CI) between two regular (evenly spaced) time series. 'RolWinWavCor' is based on the method proposed by Polanco-Martínez et al. (2018), Physica A 490, 1211-1227 (<doi:10.1016/j.physa.2017.08.065>). The RWWC was designed for financial time series analysis, but other types of data (e.g., climatic, ecological, geological, etc.) can be analysed. 'RolWinWavCor' contains two functions: (1) \code{\link{estim_RWWC}} that is used to estimate the rolling window wavelet correlation coefficients and their statistical significance, and (2) \code{\link{plot_estim_RWWC}} that plots the time series under analysis and a heat map of the rolling window correlation coefficients that are statistically significant. The functions included in 'RolWinWavCor' are highly flexible since these contains some parameters to personalize the time series under study and the heat maps of the rolling window wavelet correlations. Moreover, we have also included a data set named \code{EU_stock_markets} that contain nine EU stock market indices (closing prices) to exemplify the use of the functions contained in 'RolWinWavCor'. 
}
\details{
\tabular{ll}{
Package: \tab RolWinWavCor \cr
Type: \tab Package\cr
Version: \tab 0.4\cr
Date: \tab 2023-03-09\cr
License: \tab GPL (>= 2) \cr
LazyLoad: \tab yes\cr
}
RolWinWavCor package contains two functions: (1) \code{\link{estim_RWWC}} that estimates the rolling window wavelet correlation (WLMC) coefficients and their 95\% confidence intervals (95\% CI); (2) \code{\link{plot_estim_RWWC}} that plots the time series under analysis and the rolling window correlation coefficients that are statistical significant (within the 95\% CI) as a heat map, which are obtained via the function \code{\link{estim_RWWC}}. 
}
\note{
 Dependencies: \emph{waveslim}
}
\author{
Josué M. Polanco-Martínez (a.k.a. jomopo).\cr
Excellence Unit GECOS, IME, Universidad de Salamanca, Salamanca, SPAIN.\cr
BC3 - Basque Centre for Climate Change, Leioa, SPAIN.\cr 
Web1: \url{https://scholar.google.es/citations?user=8djLIhcAAAAJ&hl=en}.\cr
Web2: \url{https://www.researchgate.net/profile/Josue-Polanco-Martinez}.\cr
Email: \email{josue.m.polanco@gmail.com} \cr
Acknowledgement:\cr 
We acknowledge to the Excellence Unit GECOS (grant reference number CLU-2019-03), Universidad de Salamanca for its funding support. 
}  
\references{
Polanco-Martínez, J. M., Fernández-Macho, J., Neumann, M. B., & Faria, S. H. (2018). A pre-crisis vs. crisis analysis of peripheral EU stock markets by means of wavelet transform and a nonlinear causality test. Physica A: Statistical Mechanics and its Applications, 490, 1211-1227. <URL: \doi{10.1016/j.physa.2017.08.065}>. 
}
