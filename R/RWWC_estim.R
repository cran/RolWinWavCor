########################################################################
# "estim_RWWC" function from the R package RolWinWavCor  	       #
########################################################################
# The function "estim_RWWC" estimates the rolling window wavelet       #
# correlation coefficients between two regular time series. The method #
# is based on Polanco-Martínez et al (2018), Physica A 490, 211-1227,  #
# https://doi.org/10.1016/j.physa.2017.08.065 			       #
########################################################################

########################################################################
#   Copyright (C) 2023 by Josué M. Polanco-Martínez                    #
#   This file/code is part of the R RolWinWavCor package               #
########################################################################
#								     
#   RolWinWavCor is free software:  
#   you can redistribute it and/or modify it under the terms of the GNU 
#   General Public License as published by the Free Software 
#   Foundation, either version 3 of the License, or (at your option) 
#   any later version.
#
#   RolWinWavCor is distributed in the hope that it will be 
#   useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
#   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with RolWinWavCor If not, see <http://www.gnu.org/licenses/>.
#
########################################################################

#----------------------------------------------------------------------#
#                Description of arguments (INPUTS):                    #
#----------------------------------------------------------------------#
# 1. inputdata: It is a matrix of 3 columns (time and the two variables#
# 		under study).					       # 
# 2. Wname:     This is the name of the wavelet filter used in the     # 
#               MODWT decomposition of the time series under study.    #
# 3. J: 	It is the maximum level of the MODWT decomposition and #
#  		this must be an integer number.			       #
# 4. W:		Window-length or number of elements of the window where# 
# 		the rolling wavelet correlations are estimated.	       #
# 5. Align:     This is used to align the rolling object. There are    # 
#  		three options: “left”, “center” and “right,” but the   #
#	        “center” option is used by default to ensure that      # 
#	        variations in the correlations are aligned with the    #
# 	   	variations in the relationships of the variables under #
#               study, but the window-lengths (W) must be odd.	       #
# 6. Scale: 	This is used to normalise or standardise the variables #
# 		under analysis (the default option is “TRUE”; otherwise#
#		it is “FALSE”). 				       #
#----------------------------------------------------------------------#

 
# Function: start
 estim_RWWC <- function(inputdata, Wname="la8", J, W, Align="center", Scale=TRUE) { 

 #------------------------------------------------------------------------#
 # Check 1: inputdata MUST contain three columns: 
 #          time, variable X and variable Y 
 if (dim(inputdata)[2] != 2 | is.ts(inputdata) != TRUE) {
  stop("\n W A R N I N G: The input data MUST be a ts object with 
   2 columns (first and second variable under analysis). Thank you for 
  using RolWinWavCor package. \n ")
 }

 # Check 2: the times MUST be regular/evenly spaced - no gaps! 
 Deltat <- diff(time(inputdata))  # Deltat is the temporal resolution! 
 if (length(unique(Deltat)) != 1) {
  stop("\n W A R N I N G: The input data must be regular (evenly spaced 
   time). Otherwise, please, consider to address this drawback before 
   using RolWinWavCor. Thank you so much for using RolWinWavCor package. \n")
 }

 # Check 3: W MUST be odd if Align="center" 
 if (W %% 2 == 0 & Align == "center" | W < 3) { 
  stop(paste("\n W A R N I N G: W is EVEN and Align has been defined 
   as `center' or W is < 3. Thank you for using RolWinWavCor package. \n"))
  }

 # Check 4: scaling data: [X_i - mean(X_i)]/sd(X-i), mean=0 & sd=1
 if(isTRUE(Scale)) { 
  inputdata <- apply(inputdata, 2, scale)
 }
 # ----------------------------------------------------------------------- 


 NL   <- dim(inputdata)[1]
 NlW  <- NL - W + 1 
 w    <- W
 wavcor.inputdata <- array(NA, c(J, 3, NL)) 

 if (Align == "left")   { id <- 1 }
 if (Align == "center") { id <- 1; NlW <- NlW; w <- W} 
 if (Align == "right")  { id <- 1; NlW <- NlW; w <- W - 1} 


  for (i in id:NlW) { 
    modwt.tsdatX   <- modwt(inputdata[i:w,1], Wname, n.levels=J)
    modwt.tsdatY   <- modwt(inputdata[i:w,2], Wname, n.levels=J)
    bw.modwtsdatX  <- brick.wall(modwt.tsdatX, Wname)
    bw.modwtsdatY  <- brick.wall(modwt.tsdatY, Wname)
    wavcormodwtsXY <- wave.correlation(bw.modwtsdatX, bw.modwtsdatY, N=NL)
    if (Align == "right") {
     wavcor.inputdata[,,i+W-1] <- as.matrix(wavcormodwtsXY[-(J+1),])
    }
    if (Align == "center") {
     wavcor.inputdata[,,i+ceiling(W/2)-1] <- as.matrix(wavcormodwtsXY[-(J+1),])
    }
    if (Align == "left") {
     wavcor.inputdata[,,i] <- as.matrix(wavcormodwtsXY[-(J+1),])
    }
    w <- w + 1 
  } 
 return(wavcor.inputdata)

 } 
# Function: end  

