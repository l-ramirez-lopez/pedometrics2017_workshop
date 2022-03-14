#------------------------------- Info ------------------------------------------
# Description: Measuring the similarity/dissimilarity between spectra
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


# USER: specify  working directory
wd <- "C:/Users/raml/Documents/pedometrics2017"

# R: Set the working directory
setwd(wd)

# USER: specify  the input files (including the subdirectory
# that is not specified in the working directory)
inputfile_1 <- "data/world_data_3644_samples.txt"

# R: read the data
s_data <- read.table(inputfile_1,
                    header = TRUE,
                    check.names = FALSE,
                    sep = "\t"
)

# USER: indicate here the name of the first wavelength/wavenumber
# of the spectra  as it appears in the data table
first_w <- 350

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


# R: apply the Standard Normal Variate on the absorbance spectra
# in order to account for scattering effects
s_data$spc_snv <- standardNormalVariate(X = -log(s_data$spc))


# R: plot the spectra...
# USER: before plotting provide names of the axes
# of the plot
xax <- "Wavelength, nm"
yax <- "snv(Absrobance)"

matplot(x = wavs, y = t(s_data$spc_snv),
        xlab = xax,
        ylab = yax,
        type = "l",
        lty = 1,
        col = rgb(red = 0, green = 0.4, blue = 0.8, alpha = 0.5),
        main = "Soil spectra"
)

grid()

#---- 1. Mahalanobis distance computed on the principal components (PC) space ----
## First of all we need to compute the PCs of the spectra
# for doing so we can use the pc_projection function of the resemble package

# R: compute the PCs of the spectra
# USER: indicate the maximum amount of cumulative variance explained
# that needs to be retained in the PCs
max_explained_var <- 0.99

pc_spectra <- pc_projection(Xr = s_data$spc_snv,
                            pc_selection = list("cumvar", max_explained_var),
                            center = TRUE, scale = FALSE
)

# R: get a summary of the PC object created (Compactly Display the Structure of the object)
str(pc_spectra)

# R: get just the names of the sub-objects in the PC object created
names(pc_spectra)

# R: plot the first two scores of the PCs
plot(x = pc_spectra$scores[, 1],
     y = pc_spectra$scores[, 2],
     xlab = "PC 1",
     ylab = "PC 2",
     type = "p",
     pch = 16,
     col = rgb(red = 0, green = 0.4, blue = 0.8, alpha = 0.2),
     main = "Score plot"
)
grid()

# To obtain the pairwise Mahalanobis distances we can use the
# f_diss function of the resemble package

# R: compute the pairwise Mahalanobis distances
md_1 <- f_diss(
  Xr = pc_spectra$scores,
  Xu = NULL,
  diss_method = "mahalanobis",
  center = TRUE
)


# Since the Mahalanobis distance is equivalent to the Euclidean distance
# after the score variables are centered and standardized/scaled to unit variance,
# We can also apply the Euclidean distance algorithm to the centered and scale PC scores

# R: compute the pairwise Euclidean distances on the entered and scale PC scores
md_2 <- dissimilarity(
  Xr = pc_spectra$scores,
  Xu = NULL,
  diss_method = "euclid",
  center = TRUE, scale = TRUE
)

# Compare the first three pairwise distances for md_1 and md_2
md_1$dissimilarity[1:3, 1:3]
md_2$dissimilarity[1:3, 1:3]


# Let's say you just need to compare the dissimilarity matrix between
# all the samples in your spectral data set and the first three
# spectra in the same dataset

# R: compute the dissimilarity matrix
md_first_three <- dissimilarity(Xr = pc_spectra$scores,
                                Xu = pc_spectra$scores[1:3, ],
                                diss_method = "euclid",
                                center = TRUE, scale = TRUE
)

#---- 2. Correlation dissimilarity ----
# for computing the pairwise correlation dissimilarity between spectra
# you can use the correlation function of R

# R: compute the correlation dissimilarity
s_correlation <- cor(t(pc_spectra$scores))
sc_1 <- (1 - s_correlation) / 2

# R: alternatively you can also use the dissimilarity function from the
# resemble package
sc_2 <- dissimilarity(Xr = pc_spectra$scores,
                      Xu = pc_spectra$scores,
                      diss_method = "cor", 
                      center = FALSE, scale = FALSE
)

# Compare the first three pairwise distances for sc_1 and sc_2
sc_1[1:3, 1:3]
sc_2$dissimilarity[1:3, 1:3]

#---- 3. Correlation dissimilarity ----
# for computing the pairwise correlation dissimilarity between spectra
# you can use the correlation function of R

# R: compute the correlation dissimilarity
s_correlation <- cor(t(s_data$spc_snv))
cd_1 <- (1 - s_correlation) / 2

# R: alternatively you can also use the dissimilarity function in the
# resemble package
cd_2 <- dissimilarity(Xr = s_data$spc_snv,
                      Xu = s_data$spc_snv,
                      diss_method = "cor",
                      center = FALSE, scale = FALSE
)

# Compare the first three pairwise distances for sc_1 and sc_2
cd_1[1:3, 1:3]
cd_2$dissimilarity[1:3, 1:3]
round(x = cd_2$dissimilarity[1:3, 1:3], digits = 14)

# The correlation dissimilarity???? can also be applied by using a moving window approach

# USER: set a window size (the data points that will cover the window)
w_4 <- 51

# R: compute the moving window correlation dissimilarity
# (Note: it can be computationally demanding)
mw_cd <- dissimilarity(Xr = s_data$spc_snv,
                       Xu = s_data$spc_snv,
                       diss_method = "cor",
                       ws = w_4,
                       center = FALSE, scale = FALSE
)

# Check the first three dissimilarity scores
mw_cd$dissimilarity[1:3, 1:3]
round(x = mw_cd$dissimilarity[1:3, 1:3], digits = 14)

#---- 3. Cosine dissimilarity (a.k.a. Spectral Angle Mapper)----
# for computing the pairwise correlation dissimilarity between spectra
# you can use the dissimilarity function in the resemble package

# R: compute the pairwise Spectral Angle Mapper dissimilarity
sa_md <- dissimilarity(Xr = s_data$spc_snv,
                       Xu = s_data$spc_snv,
                       diss_method = "cosine",
                       center = FALSE, scale = FALSE
)

# Check the first three dissimilarity scores
sa_md$dissimilarity[1:3, 1:3]
