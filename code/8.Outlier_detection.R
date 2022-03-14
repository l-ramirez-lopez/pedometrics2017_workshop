#------------------------------- Info ------------------------------------------
# Description: Outlier detection
# 
# Inputs:      world_data_3644_samples.txt
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
# resemble
require(resemble)


# USER: specify working directory
wd <- "C:/Users/raml/Documents/pedometrics2017"

# R: Set the working directory
setwd(wd)

# USER: specify the input files (including the subdirectory 
# that is not specified in the working directory)
inputfile_1 <- "data/world_data_3644_samples.txt"

# R: read the data
s_data <- read.table(inputfile_1, 
                    header = TRUE, 
                    check.names = FALSE, 
                    sep ="\t")

# USER: indicate here the name of the first wavelength/wavenumber 
# of the spectra  as it appears in the data table
first_w <- 350

# R: extract in one object only the spectra from the "data" table...
spc <- as.matrix(s_data[,which(colnames(s_data) == first_w):ncol(s_data)])

# R: remove from "data" spectra...
s_data <- s_data[,-c(which(colnames(s_data) == first_w):ncol(s_data)), drop = FALSE]

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
        col = rgb(red = 0, green = 0.4, blue = 0.8, alpha = 0.05),
        main = "Soil spectra from the ICRAF world spectral library")
grid()


# R: apply the Standard Normal Variate on the spectra
# in order to account for scattering effects
s_data$spc_snv <- standardNormalVariate(X = s_data$spc)


# The aim of outlier detection is to identify atypical spectra (samples) 
# that are different from the average of the population of samples. Therefore,
# the distance of each sample in the dataset to the average spectrum 
# is used to assess whether the sample is an outlier or not.

# In NIR spectroscopy the outlier detection is usually carried out in the 
# principal component space of the data. Therefore, it is necessary to apply this
# technique to the spectral data before the distances are computed


# R: compute the PCs of the spectra
# USER: indicate the maximum amount of cumulative variance explained
# that needs to be retained in the PCs
maxexplvar <- 0.99

pc_spectra <- pc_projection(Xr = s_data$spc_snv, 
                            pc_selection = list("cumvar", maxexplvar), 
                            center = TRUE, scaled = FALSE)

# R: Calculate the average of the PC scores
world_spc_center <- colMeans(pc_spectra$scores)

# R: Since the result of the "colMeans" function is a vector
# it is necessary to reformat it to a matrix of 1 row
world_spc_center <- t(as.matrix(world_spc_center))


# R: Calculate the dissimilarity between each sample and the 
# average spectrum (in the PC space by using the Mahalanobis distance)
w_mahald <- f_diss(Xr =  pc_spectra$scores, 
                   Xu =  world_spc_center, 
                   diss_method = "mahalanobis", 
                   center = TRUE, scale = FALSE)


# How to decide whether a sample is an outlier or not...
# Most methods use "arbitraty" dissimilarity  limits
# over which a sample is considered as an outlier

# For example, samples with a Mahalanobis dissimilarity larger 
# than 1 can be considered as outliers

# R: plot the dissimilarity scores vs the index of the sample 
plot(w_mahald,
     pch = 16,
     col = rgb(red = 0, green = 0.4, blue = 0.8, alpha = 0.5),
     ylab = "Dissimilarity score")

# R: add a horizontal line to better visualize the samples 
# with Mahalanobis dissimilarity scores larger than 1
# (1 is the standard arbitrary threshold)
abline(h = 1, col = "red")

# R: identify the samples with Mahalanobis dissimilarity 
# scores larger than 1
w_mahald > 1

# R: Obtain the indices of the outliers
indx_out_m <- which(w_mahald > 1)

# R: who many outliers?
length(indx_out_m)

# R: plot the first two PCs along the the identified outliers
plot(x = pc_spectra$scores[,1],
     y = pc_spectra$scores[,2],
     xlab = "PC 1",
     ylab = "PC 2",
     pch = 16,
     col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5),
     main = "Mahalanobis outliers")
grid()
points(x = pc_spectra$scores[indx_out_m,1],
       y = pc_spectra$scores[indx_out_m,2],
       pch = 16,
       col = "red")


# Another method is the H distance which is 
# equivalent to the Mahalanobis distance in the PC space, 
# i.e H can be computed from it
# The H distance is computed as follows: 
# The the root squared of the product between the squared Mahalanobis distance 
# and the number of PCs
hs <- (w_mahald^2 * pc_spectra$n_components)^0.5

# R: plot the h scores vs the index of the sample 
plot(
  hs,
  pch = 16,
  col = rgb(red = 0, green = 0.8, blue = 0.3, alpha = 0.5),
  ylab = "H"
)

# R: add a horizontal line to better visualize the samples 
# with H scores larger than 3
# (3 is the standard arbitrary threshold for H used in the literature)
abline(h = 3, col = "red")

# R: identify the samples with H scores larger than 3
hs > 3

# R: Obtain the indices of the outliers
indx_out_h <- which(hs > 3)

# R: who many outliers?
length(indx_out_h)

# R: plot the first two PCs along the the identified outliers
plot(x = pc_spectra$scores[,1],
     y = pc_spectra$scores[,2],
     xlab = "PC 1",
     ylab = "PC 2",
     pch = 16,
     col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5),
     main = "H outliers")
grid()
points(x = pc_spectra$scores[indx_out_h, 1],
       y = pc_spectra$scores[indx_out_h, 2],
       pch = 16,
       col = "red")

