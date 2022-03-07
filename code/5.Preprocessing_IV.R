#------------------------------- Info ------------------------------------------
# Description: Scatter corrections using:
#              - Multiplicative Scatter Correction (msc)
#              - Standard Normal Variate (snv)
# 
# Inputs:      br_spectra_1.txt
#              br_spectra_2.txt
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

# pls
require(pls)

# USER: specify working directory
wd <- "C:/Users/raml/Documents/pedometrics2017"

# R: Set the working directory
setwd(wd)

# USER: specify the input files (including the subdirectory 
# that is not specified in the working directory)
inputfile_1 <- "data/br_spectra_1.txt"

# R: read the data
s_data <- read.table(inputfile_1, 
                    header = TRUE, 
                    check.names = FALSE, 
                    sep ="\t")

# USER: indicate here the name of the first wavelength/wavenumber 
# of the spectra  as it appears in the data table
firstw <- 352

# R: extract in one object only the spectra from the "data" table...
spc <- as.matrix(s_data[,which(colnames(s_data) == firstw):ncol(s_data)])

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
yax <- "Reflectance"

matplot(x = wavs, y = t(s_data$spc),
        xlab = xax,
        ylab = yax,
        type = "l",
        lty = 1,
        col = rgb(red = 1, green = 0, blue = 0, alpha = 0.3),
        main = "Soil spectra")
grid()

#---- 1. Continuum Removal ----
# First of all note that the technique of continuum removal
# can only be applied to either reflectance or absorbance spectra

# The continuumRemoval function of the prospectr package for applying 
# this method
# Here we remove the continuum of a set of reflectance spectra

# USER: Indicate wheter you are working with reflectance 
# or absorbance spectra. Use "R" to indicate reflectance  
# and "A" to indicate absorbance
tp <- "R"

# R: Remove the continuum from the spectra
s_data$spc_cr <- continuumRemoval(X = s_data$spc, wav = wavs, type = tp)


# R: plot the msc spectra of the old and new data
# USER: before plotting provide the new names of the y axis 
# of the plot
yax_cr <- "Reflectance (continuum removed)"


# USER: define the color of the spectra
colo <- rgb(red = 0, green = 0, blue = 0, alpha = 0.3)

# R: create the plots
matplot(x = wavs, y = t(s_data$spc_cr),
        xlab = xax,
        ylab = yax_cr,
        type = "l",
        lty = 1,
        col = colo,
        main = "Continuum Removal")
grid()



#---- 1.2 Standard Normal Variate ----
# The standardNormalVariate function from the prospectr package can 
# be used to apply this type of correction

# R: apply the Standard Normal Variate on the spectra
s_data$spc_snv <- standardNormalVariate(X = s_data$spc)
 
# R: apply the Standard Normal Variate on the new spectra 
# in this case we do not need any information form the old data
# since the transformation for each spectrum is independent from
# any other data (it operates row-wise)
new_data$spc_snv <- standardNormalVariate(X = new_data$spc)
  
# R: plot the snv spectra of the old and new data
# USER: before plotting provide the new names of the y axis 
# of the plot
yax_snv <- "snv(Reflectance)"


# R: create the plots
matplot(x = wavs, y = t(s_data$spc_snv),
        xlab = xax,
        ylab = yax_snv,
        type = "l",
        lty = 1,
        col = colo,
        main = "snv(spectra)")
grid()
matlines(x = wavs, y = t(new_data$spc_snv),
         lty = 1,
         col = colnw)
legend("bottomright", 
       legend = c("Old", "New"), 
       col = c(colo, colnw), 
       lty = 1, cex = 1, box.lty = 3, ncol = 2,
       box.col = rgb(1,1,1,0), bg = rgb(1,1,1,0))

#---- 1. Centering and scaling the spectra ----
# The scale function can be used for applying either 
# centering and scaling or both. When both are applied, 
# centering takes place before scaling

# Example 1
# USER: Define what you need:
# e.g. Only center the snv corrected spectral data
centering <- TRUE
scaling <- FALSE

# R: center the data
s_data$spc_snv_cnt <- scale(s_data$spc_snv, center = centering, scale = scaling)

# R: plot the data

# R: plot the snv spectra of the old and new data
# USER: before plotting provide the new names of the y axis 
# of the plot
yax_snv <- "Centred snv(Reflectance)"

# R: create the plots
matplot(x = wavs, y = t(s_data$spc_snv_cnt),
        xlab = xax,
        ylab = yax_snv,
        type = "l",
        lty = 1,
        col = colo,
        main = "Centred snv(spectra)")
grid()

# Example 2
# USER: Define what you need:
# e.g. Perform both centering and scaling the snv corrected spectral data
centering <- TRUE
scaling <- TRUE

# R: center the data
s_data$spc_snv_cs <- scale(s_data$spc_snv, center = centering, scale = scaling)


# R: plot the data

# R: plot the snv spectra of the old and new data
# USER: before plotting provide the new names of the y axis 
# of the plot
yax_snv <- "Centred and scaled snv(Reflectance)"

# R: create the plots
matplot(x = wavs, y = t(s_data$spc_snv_cs),
        xlab = xax,
        ylab = yax_snv,
        type = "l",
        lty = 1,
        col = colo,
        main = "Centred and scaled snv(spectra)")
grid()

