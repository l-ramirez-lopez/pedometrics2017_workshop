#------------------------------- Info ------------------------------------------
# Description: Modelling spectra
#
# Inputs:      br_spectra_3.txt
#
# Authors:     Leo Ramirez-Lopez, Melissa Lis-Gutierrez & Tatiana Moreno
#              ramirez-lopez.l@buchi.com;
#
# Date:        Dic 2021
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
# resemble
require(resemble)
# caret
require(caret)

# USER: specify the input files (including the subdirectory
# that is not specified in the working directory)
inputfile_1 <- "data/br_spectra_3.txt"

# R: read the data
s_data <- read.table(inputfile_1,
                     header = TRUE,
                     check.names = FALSE,
                     sep = "\t"
)

# USER: indicate here the name of the first wavelength/wavenumber
# of the spectra  as it appears in the data table
first_w <- 401

# R: extract in one object only the spectra from the "data" table...
spc <- as.matrix(s_data[, which(colnames(s_data) == first_w):ncol(s_data)])

# R: remove from "data" spectra...
s_data <- s_data[, -c(which(colnames(s_data) == first_w):ncol(s_data)), drop = FALSE]

# R: put back the spectra in the "data" object but as a
# sub-element of "data"
s_data$spc <- spc

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

# --- 1. Preprocessing ----
# --- 1.1. Splice correction ----
# Note the "shifts" in the plot of the spectra at 1000 and 1830 nm. 
# This is due to artefacts of the spectrometer. 
# This can be corrected with the spliceCorrection function 
# of the prospectr package

# USER: Indicate the exact wavelengths at which the shits are
# located
s_shifts <- c(997, 1833, 1837)

# R: Correct the spectral "shifts" using the spliceCorrection function 
s_data$spc_corrected <- spliceCorrection(X = s_data$spc, 
                                         wav = wavs, 
                                         splice = s_shifts
)


# R: plot the corrected spectra...                          
matplot(x = wavs, y = t(s_data$spc_corrected),
        xlab = xax,
        ylab = yax,
        ylim = c(0, 1),
        type = "l",
        lty = 1,
        col = rgb(red = 0, green = 0, blue = 1, alpha = 0.3),
        main = "Spectra corrected")



# --- 1.2. SNV ----
# R: apply the Standard Normal Variate on the absorbance spectra
# in order to account for scattering effects
s_data$spc_snv <- standardNormalVariate(X = -log(s_data$spc_corrected))

# --- 1.3. Spectra differentiation ----
# Since differentiation amplifies spectral noise we can combine
# smoothing and differentiation.
# The savitzkyGolay function can be used to apply differentiation
# It smooths out the coefficients of the polynomial in the moving window
# so that the noise is not dramatically amplified

# --- 1.3.1 First derivative ----
# USER: First define the differentiation order (e.g. first derivative)
difforder_a <- 1

# USER: Define the window size (e.g. 7 bands)...
swindow_sg <- 7

# USER: then, define the order of the polynomial
#   to fit the  points within the window (e.g. 1 )
poly_sg <- 1

# R: Apply the savitzkyGolay function to the spectra
# Note that the 2 first and last wavelengths are lost in the process
# (i.e. the vector of wavelengths is shorter now because of the window size
# of 11 bands.
# The number of wavelengths lost at the begining and
# at the the of the spectra is (window size  - 1)/2 en este caso la funci?n de orden 0
s_data$spc_snv_d1 <- savitzkyGolay(
  X = s_data$spc_snv, m = difforder_a, p = poly_sg, w = swindow_sg
)

#   R: extract the vector of remaining wavelengths
wavs_sg_diff <- colnames(s_data$spc_snv_d1)
wavs_sg_diff <- as.numeric(wavs_sg_diff)

# R: plot the spectra...
# USER: before plotting provide name of the y axis
# of the plot (since it is not absorbance anymore)

xlab <- "Wavelength, nm"
ylab <- "1st Derivarive(SNV(Absorbance))"

# R: plot the derivatives of the spectra...

matplot(x = wavs_sg_diff, y = t(s_data$spc_snv_d1),
        xlab = xlab,
        ylab = ylab,
        type = "l",
        lty = 1,
        col = rgb(red = 0, green = 0.5, blue = 0.5, alpha = 0.3),
        main = "First Savitzky-Golay derivative"
)
grid()

# --- 2. Modelling ----
control <- trainControl(method = "LGOCV",
                        number = 50,
                        p = 0.75
)

grid_pls <- data.frame(ncomp = 1:12)

model_pls <- train(x = s_data$spc_snv_d1,
                   y = s_data$Clay,
                   method = "pls",
                   trControl = control,
                   tuneGrid = grid_pls
)

model_pls$finalModel$coefficients[ , 1, 6]
plot(model_pls$finalModel$coefficients[ , 1, 6], type = "l")

#is used to input new samples and validate with the predictive model.
predict(object = model_pls, newdata = s_data$spc_snv_d1)
