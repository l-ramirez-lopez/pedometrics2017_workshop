#------------------------------- Info ------------------------------------------
# Description: Measuring the similarity/dissimilarity between spectra 
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


# R: apply the Standard Normal Variate on the absorbance spectra
# in order to account for scattering effects
sdata$spc_snv <- standardNormalVariate(X = -log(sdata$spc))


# R: plot the spectra...
# USER: before plotting provide names of the axes 
# of the plot
xax <- "Wavelength, nm"
yax <- "snv(Absrobance)"

matplot(x = wavs, y = t(sdata$spc_snv),
        xlab = xax,
        ylab = yax,
        type = "l",
        lty = 1,
        col = rgb(red = 0, green = 0.4, blue = 0.8, alpha = 0.5),
        main = "Soil spectra")
grid()

#---- 1. Mahalanobis distance computed on the principal components (PC) space ----
## First of all we need to compute the PCs of the spectra
# for doing so we can use the pcProjection function of the resemble package

# R: compute the PCs of the spectra
# USER: indicate the maximum amount of cummalative variance explained
# that needs to be retained in the PCs
maxexplvar <- 0.99

pcspectra <- pcProjection(Xr = sdata$spc_snv, 
                          pcSelection = list("cumvar", maxexplvar), 
                          center = TRUE, scaled = FALSE)

# R: get a summary of the PC object created (Compactly Display the Structure of the object)
str(pcspectra)

# R: get just the names of the sub-objects in the PC object created 
names(pcspectra)

# R: plot the first two scores of the PCs 
plot(x = pcspectra$scores[,1],
     y = pcspectra$scores[,2],
     xlab = "PC 1",
     ylab = "PC 2",
     type = "p",
     pch = 16,
     col = rgb(red = 0, green = 0.4, blue = 0.8, alpha = 0.2),
     main = "Score plot")
grid()

# To obtain the pairwise Mahalanobis distances we can use the 
# fDiss of the resemble package

# R: compute the pairwise Mahalanobis distances
md1 <- fDiss(Xr = pcspectra$scores, 
             X2 = pcspectra$scores, 
             method = "mahalanobis", 
             center = FALSE, scaled = FALSE)

# Since the Mahalanobis distance is equivalent to the Euclidean distance 
# after the score variables are centered and standardized/scaled to unit variance,
# We can also apply the Euclidean distance algorithm to the centered and scaled PC scores

# R: compute the pairwise Euclidean distances on the entered and scaled PC scores
md2 <- fDiss(Xr = pcspectra$scores, 
             X2 = pcspectra$scores, 
             method = "euclid", 
             center = TRUE, scaled = TRUE)

# Compare the first three pairwise distances for md1 and md2
md1[1:3,1:3]
md2[1:3,1:3]


# Let's say you just need to compare the dissimilarity matrix between 
# all the samples in your spectral data set and the first three
# spectra in the same dataset

#R: compute the dissimilarity matrix
md_first_three <- fDiss(Xr = pcspectra$scores, 
                        X2 = pcspectra$scores[1:3,], 
                        method = "euclid", 
                        center = TRUE, scaled = TRUE)

#---- 2. Correlation dissimilarity ----
# for computing the pairwise correlation dissimilarity between spectra 
# you can use the correlation function of R

# R: compute the correlation dissimilarity 
scorrelation <- cor(t(pcspectra$scores))
sc1 <- (1 - scorrelation)/2

# R: alternatively you can also use the corDiss from the 
# resemble package
sc2 <- corDiss(Xr = pcspectra$scores, 
               X2 = pcspectra$scores, 
               center = FALSE, scaled = FALSE)

# Compare the first three pairwise distances for sc1 and sc2
sc1[1:3,1:3]
sc2[1:3,1:3]

#---- 3. Correlation dissimilarity ----
# for computing the pairwise correlation dissimilarity between spectra 
# you can use the correlation function of R

# R: compute the correlation dissimilarity 
scorrelation <- cor(t(sdata$spc_snv))
cd1 <- (1 - scorrelation)/2

# R: alternatively you can also use the corDiss in the 
# resemble package
cd2 <- corDiss(Xr = sdata$spc_snv, 
               X2 = sdata$spc_snv, 
               center = FALSE, scaled = FALSE)

# Compare the first three pairwise distances for sc1 and sc2
cd1[1:3,1:3]
cd2[1:3,1:3]
round(x = cd2[1:3,1:3], digits = 14)

# The correlation dissimilarity???? can also be applied by using a moving window approach

# USER: set a window size (the data points that will cover the window)
w4cd <- 51 

# R: compute the moving window correlation dissimilatity 
# (Note: it can be computationally demanding)
mwcd <- corDiss(Xr = sdata$spc_snv, 
                X2 = sdata$spc_snv, 
                ws = w4cd,
                center = FALSE, scaled = FALSE)

# Check the first three dissimilatiy scores
mwcd[1:3,1:3]
round(x = mwcd[1:3,1:3], digits = 14)

#---- 3. Cosine dissimilarity (a.k.a. Spectral Angle Mapper)----
# for computing the pairwise correlation dissimilarity between spectra 
# you can use the fDiss function in the resemble package

# R: compute the pairwise Spectral Angle Mapper dissimilarity 
samd <- fDiss(Xr = sdata$spc_snv, 
              X2 = sdata$spc_snv, 
              method = "cosine", 
              center = FALSE, scaled = FALSE)

# Check the first three dissimilatiy scores
samd[1:3,1:3]



