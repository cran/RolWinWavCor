########################################################################
# "plot_estim_RWWC" function from the R package RolWinWavCor  	       #
########################################################################
# The function "plot_estim_RWWC" estimates the rolling window wavelet  #
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
# 		under study) that is used in the function "estim_RWWC".#
# 2. DATES:     This optional parameter contains the times of the time #
#               series under study. If this par. is not provided, this #
#               is computed using inputdata. 		               #
# 3. wavcorinput: This parameter contain the output of the function    #
#                 "estim_RWWC" and it is a list of 3 columns:	       # 
#                 correlation coefficients and the left and right CI.  #
# 4. Wname:     This is the name of the wavelet filter used in the     # 
#               MODWT decomposition of the time series under study.    #
# 5. J: 	It is the maximum level of the MODWT decomposition and #
#  		this must be an integer number.			       #
# 6. W:		Window-length or number of elements of the window where# 
# 		the rolling wavelet correlations are estimated.	       #
# 7. Align:     this is used to align the rolling object. There are    # 
#  		three options: “left”, “center” and “right,” but the   #
#	        “center” option is used by default to ensure that      # 
#	        variations in the correlations are aligned with the    #
# 	   	variations in the relationships of the variables under #
#               study, but the window-lengths (W) must be odd.	       #
# 8. vartsX:    Name of the first variable, e.g. "X". 		       # 
# 9. vartsY:    Name of the second variable, e.g. "Y". 		       # 
# 10. coltsX    Color used to plot the first variable. 		       #
# 11. coltsY:   Color used to plot the second variable.		       #
# 12. CEXAXIS:  This parameter is used to plot the size of the X and   #
# 		and Y axes. Its default value is 1. 		       #
# 13. CEXLAB:   This parameter is used to plot the size of the X-axis  #
# 		and Y-axis labels. Its default value is 1.	       #
#----------------------------------------------------------------------#
 

# Function: start
 plot_estim_RWWC <- function(inputdata, DATES="null", wavcorinput,   
                        Wname, J, W, Align="center", vartsX="X", 
                        vartsY="Y", coltsX="black", coltsY="blue",  
                        CEXAXIS=1, CEXLAB=1) { 

 wavcorinput <- wavcorinput

 NL   <- dim(inputdata)[1]
 Lup  <- W/2 + 1 
 Lin  <- NL - W/2 
 RealDates <- DATES[Lup:Lin]

 mgpv <- c(3,1,0)
 LWav <- dim(wavcorinput)[4]
 xx   <- RealDates; yy <- 1:J; 

 pal     <- colorRampPalette(c("blue", "cyan", "green", "white", 
    	    "yellow", "orange", "red"))
 ncolors <- 200
 breaks  <- seq(-1,1,,ncolors+1) 
 colors  <- pal(ncolors)
 # vertical colorbar 
 levs     <- breaks[-1] - diff(breaks)/2
 Lrangcol <- levs 
 
 #:: To make the barcolor
 rangecol <- (0:(ncolors-1)) / (ncolors-1)
 byseq    <- floor(ncolors/10)
 atlab    <- seq(2, length(rangecol), by=byseq)
 Lrangcol <- rangecol[atlab]
 #####################################
 to3Dp    <- as.matrix(wavcorinput[1:J,1,])
 #####################################
 # Testing if correlation values are within the CI 
 CIi      <- as.matrix(wavcorinput[1:J,2,]) 
 CIs      <- as.matrix(wavcorinput[1:J,3,]) 
 NW <- dim(wavcorinput)[3] 
 for (j in 1:J) { 
  id <- which(CIi[j,] >= to3Dp[j,] & to3Dp[j,] >= CIs[j,]) 
  to3Dp[j,id] <- NA 
 } 
 #####################################
 mat.zero <- to3Dp
 #####################################
 rangev   <- seq(min(to3Dp, na.rm=TRUE), max(to3Dp, na.rm=TRUE), 
              length.out=ncolors)
 rangebar <- matrix(rangev,nrow=1,ncol=ncolors,byrow=TRUE)
 Lrangbar <- pal(ncolors)
 Lrangbar <- rangebar[atlab]

 ################################# 
 oldpar <- par(no.readonly = TRUE)
 on.exit(par(oldpar))
 par(mar=c(2.35, 4.95, 2.2, 3.5) + 0.1) 
 layout(matrix(c(1,2,3), 3, 1, byrow=FALSE), heights=c(2, 3.5, 0.65))
 ################################# 
 
 # To generate D1, D2,..., DJ
 Ds <- paste("D", 1:J, sep="")

 # Lab names in case names are not defined! 
 labynames <- colnames(inputdata)
 vartsX    <- labynames[1]
 vartsY    <- labynames[2]

 # Plot time series 
 plot(DATES, inputdata[,1], t="l", xaxs="i", col=coltsX, las=1, 
  xlab="Time", ylab="", yaxt="n", cex.axis=CEXAXIS, 
  main=paste(vartsX, " vs. ", vartsY, sep="")) 
 points(DATES, inputdata[,2], t="l", xaxs="i", col=coltsY, yaxt="n") 
 axis(2, at=pretty(inputdata[,1]), labels=pretty(inputdata[,1]), 
  col.lab=coltsX, las=1, cex.axis=CEXAXIS)
 mtext(2, text=vartsX, col=coltsX, line=2.5, cex=CEXLAB, las=3)
 axis(4, at=pretty(inputdata[,2]), labels=pretty(inputdata[,2]), 
  col.lab=coltsY, col.axis=coltsY, las=1, cex.axis=CEXAXIS)
 mtext(4, text=vartsY, col=coltsY, line=2.5, cex=CEXLAB, las=3)
 mtext(1, text="Time", line=2.75, cex=CEXLAB)
 # Plot heat maps of wavelet rolling window correlations 
 image(DATES, yy, z=t(mat.zero), col=colors, xlab=" ", cex.main=0.5, 
  breaks=breaks, ylab="Wavelet scale", yaxt="n", 
  cex.lab=1.6*CEXLAB, cex.axis=CEXAXIS, cex.main=1.2, main="") 
 axis(2, at=1:J, labels=Ds[1:J], las=2, cex.axis=CEXAXIS, cex.lab=CEXLAB)
 # Plot color bar 
 image(z=t(rangebar), axes=FALSE, col=colors, frame.plot=TRUE, 
  breaks=breaks, yaxt="n", xaxt="n")
 mtext(1, text="Wavelet correlation coefficients", line=-2.5, 
  cex=0.75*CEXLAB)
 axis(1, at=Lrangcol, labels=round(Lrangbar, digits=2), las=1,  
  cex.axis=0.85*CEXAXIS, gap.axis=0.35)

 } 
# Function: end  

 
 
