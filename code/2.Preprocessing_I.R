#------------------------------- Info ------------------------------------------
# Description: Spectral transformation from reflectance to absorbance 
#              Spectral de-noising
#              Spectral resampling
#              Spectral differentiation
# 
# Inputs:      noisy_spectra.txt
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
inputfile1 <- "data/noisy_spectra.txt"

# read the data
sdata <- read.table(inputfile1, 
                    header = TRUE, 
                    check.names = FALSE, 
                    sep ="\t")

# USER: indicate here the name of the first wavelength/wavenumber 
# of the spectra  as it appears in the data table
firstw <- 350


# R: extract in one object only the spectra from the "data" table...
spc <- sdata[,which(colnames(sdata) == firstw):ncol(sdata)]


# --- 1. Transformation from reflectance to absorbance ----
# R: Transform the reflectance spectra into absorbance spectra...
# log(1/R) or -log(R)
spc <- -log(spc)

# R: remove from "data" spectra...
sdata <- sdata[,-c(which(colnames(sdata) == firstw):ncol(sdata)), drop = FALSE]

# R: put back the spectra in the "data" object but as a 
# sub-element of "data"
sdata$spc <-spc

# R: remove the spc object since it is already a sub-element of "data"
rm(spc)


# --- 2. Noise removal ----
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
yax <- "Absorbance"


matplot(x = wavs, y = t(sdata$spc),
        xlab = xax,
        ylab = yax,
        type = "l",
        lty = 1,
        col = rgb(red = 1, green = 0, blue = 0, alpha = 0.3),
        main = "Noisy spectra")
grid()

# Note the random noise along all the spectra


# We present two options

# First option: Reduce the noise by applying a moving average window 
#   using the movav function 
#   USER: First define the window size (e.g. 11 bands)
swindow_ma <- 11

#   R: Apply the movav function to the spectra
#   Note that the 5 first and last wavelengths are lost in the process
#   (i.e. the vector of wavelengths is shorter now because of the window size
#   of 11 bands)
# The number of wavelengths lost at the begining and 
# at the the of the spectra is (window size  - 1)/2
sdata$spc_ma <- movav(X = sdata$spc, w = swindow_ma)  

#   R: extract the vector of remaining wavelengths
wavs_ma <- colnames(sdata$spc_ma)
wavs_ma <- as.numeric(wavs_ma)



# Second option: Reduce the noise by applying the Savitzky-Golay filter 
#   (savitzkyGolay function) 
#   USER: First define the window size (e.g. 11 bands)...
swindow_sg <- 11
#   USER: then, define the order of the polynomial 
#   to fit the  points within the window (e.g. 2 )
poly_sg <- 2


#   R: Apply the savitzkyGolay function to the spectra
#   Note that the 5 first and last wavelengths are lost in the process
#   (i.e. the vector of wavelengths is shorter now because of the window size
#   of 11 bands)
sdata$spc_sg <- savitzkyGolay(X = sdata$spc, m = 0, p = poly_sg, w = swindow_sg) 

#   R: extract the vector of remaining wavelengths
wavs_sg <- colnames(sdata$spc_sg)
wavs_sg <- as.numeric(wavs_sg)


# R: plot the de-noised spectra resulting from applying the moving average...                          
matplot(x = wavs_ma, y = t(sdata$spc_ma),
        xlab = xax,
        ylab = yax,
        type = "l",
        lty = 1,
        col = rgb(red = 0.5, green = 0, blue = 0.5, alpha = 0.3),
        main = "De-noised spectra (moving average)")
grid()

# R: plot the de-noised spectra resulting from applying the Savitzky-Golay filter...                          
matplot(x = wavs_sg, y = t(sdata$spc_sg),
        xlab = xax,
        ylab = yax,
        type = "l",
        lty = 1,
        col = rgb(red = 0, green = 0.5, blue = 0.5, alpha = 0.3),
        main = "De-noised spectra (Savitzky-Golay)")
grid()

