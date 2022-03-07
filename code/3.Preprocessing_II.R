#------------------------------- Info ------------------------------------------
# Description: Spectral differentiation
# 
# Inputs:      spectra_soil_profile_Colombia.txt
#
# Authors:     Leo Ramirez-Lopez & Alexandre Wadoux
#              ramirez-lopez.l@buchi.com; alexandre.wadoux@wur.nl 
#
# Date:        Jun 2017
#
# Actualization: Melissa Lis-Gutierrez, Tatiana Moreno & Leo Ramirez-Lopez
#               mlisg@unal.edu.co; tmorenom@unal.edu.co; 
#               ramirez-lopez.l@buchi.com
#
# Date:        Feb 2022
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set the language of R to English
Sys.setenv(language = "EN")

# Call the required packages 
# prospectr
require(prospectr)

# USER: specify working directory
wd <- "C:/Users/raml/Documents/pedometrics2017"

# R: Set the working directory
setwd(wd)

# USER: specify the input files (including the subdirectory 
# that is not specified in the working directory)
inputfile_1 <- "data/spectra_soil_profile_Colombia.txt"

# read the data
s_data <- read.table(inputfile_1, 
                    header = TRUE, 
                    check.names = FALSE, 
                    sep ="\t")

# USER: indicate here the name of the first wavelength/wavenumber 
# of the spectra  as it appears in the data table
firstw <- 350


# R: extract in one object only the spectra from the "data" table...
spc <- s_data[,which(colnames(s_data) == firstw):ncol(s_data)]


# --- 1. Transformation from reflectance to absorbance ----
# R: Transform the reflectance spectra into absorbance spectra...
# log(1/R) or -log(R)
spc <- -log(spc)

# R: remove from "data" spectra...
s_data <- s_data[,-c(which(colnames(s_data) == firstw):ncol(s_data)), drop = FALSE]

# R: put back the spectra in the "data" object but as a 
# sub-element of "data"
s_data$spc <-spc

# R: remove the spc object since it is already a sub-element of "data"
rm(spc)

# R: extract from the column names of the spectra sub-element 
# the vector of wavelengths/wavenumbers
wavs <- colnames(s_data$spc)


# R: Since the "wavs" vector is a character string vector
# we will need to transform it to a numeric vector
# NOTE that the names of the columns of the spectra 
# must be written only with numbers (otherwise they will not be
# correctly converted from characters to numbers)
wavs <- as.numeric(wavs)


# R: plot the spectra...
# USER: before plotting provide names of the axes 
# of the plot
xax <- "Wavelength, nm"
yax <- "Absorbance"


matplot(x = wavs, y = t(s_data$spc),
        xlab = xax,
        ylab = yax,
        type = "l",
        lty = 1,
        col = rgb(red = 1, green = 0, blue = 0, alpha = 0.3),
        main = "Spectra of a soil profile")
grid()

# --- 2. Spectra differentiation ----
# Since differentiation amplifies spectral noise we can combine 
# smoothing and differentiation.
# The savitzkyGolay function can be used to apply differentiation
# It smooths out the coefficients of the polynomial in the moving window
# so that the noise is not dramatically amplified

# --- 2.1 First derivative ----
# USER: First define the differentiation order (e.g. first derivative)
difforder_a <- 1

# USER: Define the window size (e.g. 5 bands)...
swindow_sg <- 5

# USER: then, define the order of the polynomial 
#   to fit the  points within the window (e.g. 2 )
poly_sg <- 2

# R: Apply the savitzkyGolay function to the spectra
# Note that the 2 first and last wavelengths are lost in the process
# (i.e. the vector of wavelengths is shorter now because of the window size
# of 11 bands. 
# The number of wavelengths lost at the beginning and 
# at the the of the spectra is (window size  - 1)/2
s_data$spc_sgdiff_1st <- savitzkyGolay(X = s_data$spc, m = difforder_a, p = poly_sg, w = swindow_sg) 

#   R: extract the vector of remaining wavelengths
wavs_sgdiff <- colnames(s_data$spc_sgdiff_1st)
wavs_sgdiff <- as.numeric(wavs_sgdiff)

# R: plot the spectra...
# USER: before plotting provide name of the y axis 
# of the plot (since it is not absorbance anymore)
yaxder_a <- "1st Derivarive"

# R: plot the derivatives of the spectra...                          
matplot(x = wavs_sgdiff, y = t(s_data$spc_sgdiff_1st),
        xlab = xax,
        ylab = yaxder_a,
        type = "l",
        lty = 1,
        col = rgb(red = 0, green = 0.5, blue = 0.5, alpha = 0.3),
        main = "First Savitzky-Golay derivative")
grid()


# --- 2.2 Second derivative ----
# USER: First define the differentiation order (e.g. first derivative)
difforder_b <- 2

# R: Apply the savitzkyGolay function to the spectra
# Note that the 2 first and last wavelengths are lost in the process
# (i.e. the vector of wavelengths is shorter now because of the window size
# of 11 bands. 
# The number of wavelengths lost at the beginning and 
# at the the of the spectra is (window size  - 1)/2
s_data$spc_sgdiff_2nd <- savitzkyGolay(X = s_data$spc, m = difforder_b, p = poly_sg, w = swindow_sg) 

# R: plot the spectra...
# USER: before plotting provide name of the y axis 
# of the plot (since it is not absorbance anymore)
yaxder_b <- "2nd Derivarive"

# R: plot the derivatives of the spectra...                          
matplot(x = wavs_sgdiff, y = t(s_data$spc_sgdiff_2nd),
        xlab = xax,
        ylab = yaxder_b,
        type = "l",
        lty = 1,
        col = rgb(red = 0, green = 0.5, blue = 0.5, alpha = 0.3),
        main = "Second Savitzky-Golay derivative")
grid()
