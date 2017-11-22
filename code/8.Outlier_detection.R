#------------------------------- Info ------------------------------------------
# Description: Outlier detection
# 
# Inputs:      world_data_3644_samples.txt
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
# resemble
require(resemble)


# USER: specifiy working directy
wd <- "C:/Users/raml/Documents/pedometrics2017"

# R: Set the working directory
setwd(wd)

# USER: specifiy the input files (including the subdirectory 
# that is not specified in the working directy)
inputfile1 <- "data/world_data_3644_samples.txt"

# R: read the data
sdata <- read.table(inputfile1, 
                    header = TRUE, 
                    check.names = FALSE, 
                    sep ="\t")

# USER: indicate here the name of the first wavelength/wavenumber 
# of the spectra  as it appears in the data table
firstw <- 350

# R: extract in one object only the spectra from the "data" table...
spc <- as.matrix(sdata[,which(colnames(sdata) == firstw):ncol(sdata)])

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
        type = "l",
        lty = 1,
        col = rgb(red = 0, green = 0.4, blue = 0.8, alpha = 0.05),
        main = "Soil spectra from the ICRAF world spectral library")
grid()


# R: apply the Standard Normal Variate on the spectra
# in order to account for scattering effects
sdata$spc_snv <- standardNormalVariate(X = sdata$spc)


# The aim of outlier detection is to identify atypical spectra (samples) 
# that are different from the average of the population of samples. Therefore,
# the distance of each sample in the dataset to the average spectrum 
# is used to assess whether the sample is an outlier or not.

# In NIR spectroscopy the outlier detection is usually carried out in the 
# principal compüonent space of the data. Therefore, it is nessary to apply this
# technique to the spectral data before the distances are computed


# R: compute the PCs of the spectra
# USER: indicate the maximum amount of cummalative variance explained
# that needs to be retained in the PCs
maxexplvar <- 0.99

pcspectra <- pcProjection(Xr = sdata$spc_snv, 
                          pcSelection = list("cumvar", maxexplvar), 
                          center = TRUE, scaled = FALSE)

# R: Calculate the average of the PC socres
worldspccenter <- colMeans(pcspectra$scores)

# R: Since the result of the "colMeans" function is a vector
# it is necssary to reformat it to a matrix of 1 row
worldspccenter <- t(as.matrix(worldspccenter))


# R: Calculate the dissimilarity between each sample and the 
# average spectrum (in the PC space by using the Mahalanobis distance)
wmahald <- fDiss(Xr =  pcspectra$scores, 
                 X2 =  worldspccenter, 
                 method = "mahalanobis", 
                 center = FALSE, scaled = FALSE)


# How to decide whether a sample is an outlier or not...
# Most methods use "arbitraty" dissimilarity  limits
# over which a sample is considered as an outlier

# For example, samples with a Mahalanobis dissimilarity larger 
# than 1 can be considered as outliers

# R: plot the dissimilarity scores vs the index of the sample 
plot(wmahald,
     pch = 16,
     col = rgb(red = 0, green = 0.4, blue = 0.8, alpha = 0.5),
     ylab = "Dissimilarity score")

# R: add a horizontal line to better vizualize the samples 
# with Mahalanobis dissimilarity scores larger than 1
# (1 is the standard arbitrary threshold)
abline(h = 1, col = "red")

# R: identify the samples with Mahalanobis dissimilarity 
# scores larger than 1
wmahald > 1

# R: Obtain the indices of the outliers
indx_out_m <- which(wmahald > 1)

# R: whoh many outliers?
length(indx_out_m)

# R: plot the first two PCs along the the identified outliers
plot(x = pcspectra$scores[,1],
     y = pcspectra$scores[,2],
     xlab = "PC 1",
     ylab = "PC 2",
     pch = 16,
     col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5),
     main = "Mahalanobis outliers")
grid()
points(x = pcspectra$scores[indx_out_m,1],
       y = pcspectra$scores[indx_out_m,2],
       pch = 16,
       col = "red")


# Another method is the H distance which is 
# equivalent to the Mahalanobis distance in the PC space, 
# i.e H can be computed from it
# The H distance is computed as follows: 
# The the root squared of the product between the squared Mahalanobis distance 
# and the number of PCs
hs <- (wmahald^2 * pcspectra$n.components)^0.5

# R: plot the h scores vs the index of the sample 
plot(hs,
     pch = 16,
     col = rgb(red = 0, green = 0.8, blue = 0.3, alpha = 0.5),
     ylab = "H")

# R: add a horizontal line to better vizualize the samples 
# with H scores larger than 3
# (3 is the standard arbitrary threshold for H used in the literature)
abline(h = 3, col = "red")

# R: identify the samples with H scores larger than 3
hs > 3

# R: Obtain the indices of the outliers
indx_out_h <- which(hs > 3)

# R: whoh many outliers?
length(indx_out_h)

# R: plot the first two PCs along the the identified outliers
plot(x = pcspectra$scores[,1],
     y = pcspectra$scores[,2],
     xlab = "PC 1",
     ylab = "PC 2",
     pch = 16,
     col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5),
     main = "H outliers")
grid()
points(x = pcspectra$scores[indx_out_h,1],
       y = pcspectra$scores[indx_out_h,2],
       pch = 16,
       col = "red")


