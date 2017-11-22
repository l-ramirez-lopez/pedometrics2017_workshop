#------------------------------- Info ------------------------------------------
# Description: Splice correction for systematic shits observed in the data
#              at specific wavelengths/wavenumbers
# 
# Inputs:      sensorRepetitions.txt
#
# Authors:     Leo Ramirez-Lopez & Alexandre Wadoux
#              ramirez-lopez.l@buchi.com; alexandre.wadoux@wur.nl 
#
# Date:        Jun 2017
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set the language of R to English
Sys.setenv(language = "EN")

# Call the required packages 
# prospectr
require(prospectr)


# USER: specifiy working directy
wd <- "C:/Users/raml/Documents/pedometrics2017"

# R: Set the working directory
setwd(wd)

# USER: specifiy the input files (including the subdirectory 
# that is not specified in the working directy)
inputfile1 <- "data/sensorRepetitions.txt"

# R: read the data
sdata <- read.table(inputfile1,
                    header = TRUE, 
                    check.names = FALSE, 
                    sep ="\t")

# USER: indicate here the name of the first wavelength/wavenumber 
# of the spectra  as it appears in the data table
firstw <- 350


# R: extract in one object only the spectra from the "data" table...
spc <- sdata[,which(colnames(sdata) == firstw):ncol(sdata)]

# R: remove from "data" spectra...
sdata <- sdata[,-c(which(colnames(sdata) == firstw):ncol(sdata)), drop = FALSE]

# R: put back the spectra in the "data" object but as a 
# sub-element of "data"
sdata$spc <-spc

# R: remove the spc object since it is already a sub-element of "data"
rm(spc)

# R: extract from the column names of the spectra sub-element 
# the vector of wavelengths/wavenumbers
wavs <- colnames(sdata$spc)

# R: Since the "wavs" vector is a character string vector
# we will need to transform it to a numeric vector
# NOTE that the names of the columns of the spectra 
# must be writen only with numbers (otherwise they will not be
# correctly converted from characters to numbers)
wavs <- as.numeric(wavs)


# R: plot the spectra...
# USER: before plotting provide names of the axes 
# of the plot
xax <- "Wavelength, nm"
yax <- "Reflectance"


matplot(x = wavs, y = t(sdata$spc),
        xlab = xax,
        ylab = yax,
        ylim = c(0, 1),
        type = "l",
        lty = 1,
        col = rgb(red = 1, green = 0, blue = 0, alpha = 0.3),
        main = "Spectra with shifts")


# Note the "shifts" in the plot of the spectra at 1000 and 1830 nm. 
# This is due to artefacts of the spectrometer. 
# This can be corrected with the spliceCorrection function 
# of the prospectr package

# USER: Indicate the exact wavelengths at which the shits are
# located
sshifts <- c(1000, 1830)

# R: Correct the spectral "shifts" using the spliceCorrection function 
sdata$spc_corrected <- spliceCorrection(X = sdata$spc, 
                                       wav = wavs, 
                                       splice = sshifts)


# R: plot the corrected spectra...                          
matplot(x = wavs, y = t(sdata$spc_corrected),
        xlab = xax,
        ylab = yax,
        ylim = c(0, 1),
        type = "l",
        lty = 1,
        col = rgb(red = 0, green = 0, blue = 1, alpha = 0.3),
        main = "Spectra corrected")
