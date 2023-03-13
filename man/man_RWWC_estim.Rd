\name{estim_RWWC}
\alias{estim_RWWC}
\alias{estimation_RWWC}
\title{Estimate rolling window wavelet correlation coefficients  
}
\description{The \code{\link{estim_RWWC}} function estimates the rolling window wavelet correlation (RWWC) coefficients between two regular or evenly spaced time series. This function is based on the work of Polanco-Martínez et al. (2018). However, we use the functions \code{modwt} and \code{wave.correlation} from the \R package 'waveslim' (Whitcher et al. 2000, Whitcher 2022) to decompose via MODWT the two time series under study and to estimate the wavelet correlation between these two time series. 
}
\usage{
estim_RWWC(inputdata, Wname="la8", J, W, Align="center", Scale=TRUE)
}
\arguments{
\item{inputdata}{ 
A matrix of three columns: the first one is the time (regular or evenly spaced) and the other two columns are the variables under study. 
}
\item{Wname}{
Name of the wavelet filter used in the wavelet transform (MODWT) decomposition. There are several wavelet filters to be used, but we use by default the Daubechies orthonormal compactly supported wavelet of length L=8, that is, ``\code{la8}''. A relatively long wavelet filter (e.g. LA(4) or LA(8)) is adequate to analyse non-stationary time series but also correlation structures that are not stationary, as was pointed out by Fernández-Macho (2018) and Polanco-Martínez et al. (2020). 
}
\item{J}{
The maximum level of the MODWT decomposition. It is recommended to use \eqn{log_{2}(N) - 3}, where \eqn{N} is the number of rows or elements of \code{inputdata} (Fernández-Macho 2018, Polanco-Martínez et al. 2020). 
}
\item{W}{
The window-length or size of the window used when the rolling window wavelet correlation coefficients are estimated. The optimal window-length is a key question when the RWWC is estimated and the answer is closely related with the time-scales of the phenomenon under study. It is advisable to use one window-length around the target time-scale, not too short and not too long, and always test several window-lengths to analyse its sensibility (Polanco-Martínez et al. (2018)). 
}
\item{Align}{
This is used to align the rolling object. There are three options: \code{left}, \code{center} and \code{right}, but the \code{center} option is used by default to ensure that variations in the correlations are aligned with the variations in the relationships of the variables under study, rather than being shifted to the left or right (Polanco-Martínez 2019, Polanco-Martínez 2020). However, this implies that the window-length must be an odd number.
}
\item{Scale}{
Scale (by default is ``TRUE''; ``FALSE''  otherwise) is used to ``normalize'' or ``standardize'' the variables under analysis. It is advisable to ''normalize/standardize'' the time series under study to have them in the same scales. 
}
}
\details{
The \code{\link{estim_RWWC}} function estimates the rolling window wavelet correlation coefficients and their statistical significance between two time series sampled on identical time points for a given window-length. The \code{\link{estim_RWWC}} function is based on the work of Polanco-Martínez et al. (2018), but this function uses the functions \code{modwt} and \code{wave.correlation} from the \R package 'waveslim' package (Whitcher et al. 2000, Whitcher et al. 2022) to decompose via MODWT the two time series under analysis and to compute the wavelet correlation coefficients between these two series. 
}
\value{
Outputs:
A multidimensional matrix containing three columns: the rolling window wavelet correlation coefficients and their respective 95\% confidence interval in the second (lower bound) and third (upper bound) columns. 
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
Polanco-Martínez, J. M., Fernández-Macho, J., Neumann, M. B., & Faria, S. H. (2018). A pre-crisis vs. crisis analysis of peripheral EU stock markets by means of wavelet transform and a nonlinear causality test. Physica A: Statistical Mechanics and its Applications, 490, 1211-1227. <URL: \doi{10.1016/j.physa.2017.08.065}>. \cr 

Fernández-Macho, J. (2018). Time-localized wavelet multiple regression and correlation. Physica A: Statistical Mechanics and its Applications, 492, 1226-1238. <URL: \doi{10.1016/j.physa.2017.11.050}>. \cr 

Whitcher, B., Guttorp, P., & Percival, D. B. (2000). Wavelet analysis of covariance with application to atmospheric time series. Journal of Geophysical Research: Atmospheres, 105(D11), 14941-14962. <URL: \doi{10.1029/2000JD900110}>. \cr

Gencay, R., F. Selcuk and B. Whitcher (2001) An Introduction to Wavelets and Other Filtering Methods in Finance and Economics, Academic Press. \cr 

Whitcher, B. (2022). waveslim: Basic Wavelet Routines for One-, Two-, and Three-Dimensional Signal Processing.  <URL: https://cran.r-project.org/package=waveslim>, R package version 1.8.4.
}
\examples{
# We reproduce Figure 2 presented in Polanco-Martínez et al. (2018). 
datPIGS     <- EU_stock_markets 
sindatePIGS <- datPIGS[-1] # To remove dates 
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
}
\keyword{estimate_the_RWWC}
\keyword{estimate_the_Rolling_Window_Wavelet_Correlation} 
