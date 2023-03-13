\name{plot_estim_RWWC}
\alias{plot_estim_RWWC}
\alias{plt_estimation_RWWC}
\title{Plot the rolling window wavelet correlation coefficients 
}
\description{The \code{\link{plot_estim_RWWC}} function plots the rolling window wavelet correlation (RWWC) coefficients that are statistically significant between two regular time series as a heat map and also plots the time series under study. The function is based on the work of Polanco-Martínez et al. (2018). \code{\link{plot_estim_RWWC}} is fed by the function \code{\link{estim_RWWC}} that is contained in our \R package 'RolWinWavCor'. 
}
\usage{
plot_estim_RWWC(inputdata, DATES="null", wavcorinput, Wname, J, W, 
  Align="center", vartsX="X", vartsY="Y", coltsX="black", 
  coltsY="blue", CEXAXIS=1, CEXLAB=1)
}
\arguments{
\item{inputdata}{ 
A matrix of three columns: the first one is the time (regular or evenly spaced) and the other two columns are the variables under study. This data set is the same used in the function \code{\link{estim_RWWC}}. 
}
\item{DATES}{
This optional parameter contains the times of the time series under study. If this parameter is not provided it is computed using the \code{inputdata}.
}
\item{wavcorinput}{
This parameter contains the output of the function \code{\link{estim_RWWC}} and consists of a multidimensional matrix with three columns: the first one contains the correlation coefficients and the second and third columns are the left and right bound of the 95\% confidence interval. 
}
\item{Wname}{
Name of the wavelet filter used in the wavelet transform (MODWT) decomposition and must be the same that was used with the function \code{\link{estim_RWWC}}. 
}
\item{J}{
The maximum level of the MODWT decomposition and must be the same used with the function \code{\link{estim_RWWC}}. 
}
\item{W}{
The window-length or size of the window used when the rolling window wavelet correlation coefficients are estimated and this must have the same value that was used in \code{\link{estim_RWWC}}.
}
\item{Align}{
This is used to align the rolling object and must be the same as used in the function \code{\link{estim_RWWC}}.
}
\item{vartsX, vartsY}{
Names of the first (e.g. ``X'') and the second (e.g. ``Y'') variable under study.
}
\item{coltsX, coltsY}{ 
The colors used to plot the first and second variable. By default the colors are black and blue for the first and second variable, respectively. 
}
\item{CEXAXIS}{ 
This parameter is used to plot the size of the X and Y axes. Its default value is 1.
}
\item{CEXLAB}{
This parameter is used to plot the size of the X-axis and Y-axis labels. Its default value is 1.
}
}
\details{
The \code{\link{plot_estim_RWWC}} function plots the time series under analysis and the rolling window wavelet correlation coefficients that are statistically significant (within the 95\% CI) as a heat map. This function is also based on the work of Polanco-Martínez et al. (2018).
}
\value{
Output: a multi-plot displayed via screen containing the time series under scrutiny and a heat map of the rolling window wavelet correlation coefficients that are statistically significant. 
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
\examples{
# We reproduce Figure 2 presented in Polanco-Martínez et al. (2018). 
datPIGS     <- EU_stock_markets 
sindatePIGS <- datPIGS[-1]    
sindatePIGS <- sindatePIGS[c(1:5, 8)]
lrdatPIGS   <- apply(log(sindatePIGS), 2, diff)
lrDATES     <- as.Date(datPIGS[,1][-1])
tsdatPIGS   <- ts(lrdatPIGS, start=1, freq=1)
Nnam        <- dim(tsdatPIGS)[2]
lrdatPIGS   <- lrdatPIGS[,1:Nnam]
inputdata   <- tsdatPIGS[,c(2,5)]

Wname <- "la8"  
J     <- 4    
W     <- 241 
Align <- "center"

rwwc <- estim_RWWC(inputdata, Wname, J, W, Align=Align)
wavcor.output <- rwwc 
DATES         <- lrDATES
plot_estim_RWWC(inputdata, DATES=DATES, wavcor.output, Wname, J, W, 
                  Align=Align, CEXAXIS=1.2)
}
\keyword{plot_estimate_the_RWWC}
\keyword{plot_estimate_the_Rolling_Window_Wavelet_Correlation} 
