#------------------------------- Info ------------------------------------------
# Description: Assesing the similarity/dissimilarity between spectra 
#              based on the the Mahalanobis distance computed on the 
#              principal component space
# 
# Inputs:      br_spectra_3.txt
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
inputfile1 <- "data/br_spectra_3.txt"

# R: read the data
sdata <- read.table(inputfile1, 
                    header = TRUE, 
                    check.names = FALSE, 
                    sep ="\t")

# USER: indicate here the name of the first wavelength/wavenumber 
# of the spectra  as it appears in the data table
firstw <- 401

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


# R: get a summary of the sdata object created (Compactly Display the Structure of the object)
str(sdata)

# This is a dataset that includes soil data for particle size distribution (sand, silt and clay),
# and vis-NIR spectra
# In this example we will use the particle size distribution data as side information to
# help us to asses the dissimilatity measures that we get from the spectra
# In this approach the side information means information about one variable which is available 
# for a group of samples. It is assumed that there is a correlation (or at least an indirect 
# or secondary correlation) between this side information and soil spectra. In other words, 
# this approach is based on the assumption that the similarity measures between the spectra 
# of a given group of soil samples should be able to reflect their similarity also in terms 
# of the side information (e.g. compositional similarity).



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


#---- 1. Random selection of samples ----

# For this example randomly split the data into two substes
# (select 25% of the data)
# R: Obtain the number of samples in the data
ns <- nrow(sdata)

# USER: indicate the percentaje that you want to select
pd <- 0.25

# R: compute the integer number corresponding to the 25% of the samples
nsamples2select <- round(x = ns * pd, digits = 0)

# R: compute the vector of sample indexes
orig_indexes <- 1:ns
 
# R: Randomly select the samples
# before the selection use a random initialization number (seed). 
# This allows you to reproduce your random selection. For this,
# the set.seed function of R can be used. You just need to provide any 
# number you want. In this example we use Leo's birthday (in case you 
# want to send gifts to Switzerland on this date!) ;)
# Note that if you want to reproduce the random slection on the same vector,
# you need to use the same random initialization number right after the selection
# in conjunction with the set.seed function
set.seed(241180) 
sel_sample_indx <- sample(x = orig_indexes, size = nsamples2select)

# R: now split the data
# first the indexes randomly found
xu_data <- sdata[sel_sample_indx,]
# Then the rest (i.e. the ones that were not randomly found)
xr_data <- sdata[-sel_sample_indx,]

# Let's assume that the side information is available only for the 
# xr_data (along with the spectral data) and that the  xr_data only
# contains the spectral information (i.e. the side information is not available 
# for this set)


#---- 2. Mahalanobis distance computed on the principal components (PC) space ----
## First of all we need to compute the PCs of the spectra
# for doing so we can use the pcProjection function of the resemble package

# Since principal component anylsisis for the spectra does not require
# any side information for its execution, then we can apply the principal componets
# to the combined spectra (xr_data  and xu_data)

# R: combine the spectral data of xr_data  and xu_data (by rows) into one single dataset
# using the rbind function of R
comb_x <- rbind(xr_data$spc_snv, xu_data$spc_snv)

# R: compute the PCs of the combined spectra
# USER: indicate the maximum amount of cummalative variance explained
# that needs to be retained in the PCs
maxexplvar <- 0.99

pcspectra <- pcProjection(Xr = comb_x, 
                          pcSelection = list("cumvar", maxexplvar), 
                          center = TRUE, scaled = FALSE)

# R: get a summary of the PC object created (Compactly Display the Structure of the object)
str(pcspectra)

# R: get just the names of the sub-objects in the PC object created 
names(pcspectra)

# The first 237 rows of the score matrix correspond to the socres of the
# spectral data of xr_data and the remaining ones correspond to the socres of the 
# xu_data



# R: plot the first two scores of the PCs of the xr_data
plot(x = pcspectra$scores[1:nrow(xr_data), 1],
     y = pcspectra$scores[1:nrow(xr_data), 2],
     xlab = "PC 1",
     ylab = "PC 2",
     type = "p",
     pch = 16,
     col = rgb(red = 0, green = 0.4, blue = 0.8, alpha = 0.5),
     main = "Score plot")
grid()

# R: Add to the above plot the first two scores of the PCs of the xu_data
points(x = pcspectra$scores[-c(1:nrow(xr_data)), 1],
       y = pcspectra$scores[-c(1:nrow(xr_data)),2],
       xlab = "PC 1",
       ylab = "PC 2",
       pch = 16,
       col = rgb(red = 0.8, green = 0.4, blue = 0, alpha = 0.5))

# Since the xr_data is the dataset that contains the side information
# we execute the side information analysis on the dissimilarity matrix 
# of the socres of xr_data (i.e. first 237 rows of  pcspectra$scores)

# Step 1: A distance matrix is derived from the spectral matrix X. 
# This spectral matrix has a side information Y.

# R: First we compute the pairwise Mahalanobis distances
# To obtain the pairwise Mahalanobis distances we can use the 
# fDiss of the resemble package
md_xr <- fDiss(Xr =  pcspectra$scores[1:nrow(xr_data), ], 
               X2 =  pcspectra$scores[1:nrow(xr_data), ], 
               method = "mahalanobis", 
               center = FALSE, scaled = FALSE)

# md_xr is the dissimilarity matrix of the spectra in xr_data
# computed in its PC space

# Step 2. By using the distance matrix, for each sample in X select its closet 
# (most spectrally similar) sample.

# USER: First create an empty object to store the results of the indices of 
# nearest neighbors 
nearest_n <- NULL
# R: use a for loop to iterate over each of the columns in the 
# md_xr matrix
for(i in 1:nrow(md_xr)){
  # get the order of the indices of the 
  # values from the smallest to the largest
  # and select only the second one (the nearest neighbour)
  # Since the dissimilarity of a sample to itself it is zero
  # then the first one in the ordered list does not indicate the 
  # nearest neighbor but the sample itself
  
  # R: Temporarily store the results of the ith iteration in an object (nn_i)
  nn_i <- order(x = md_xr[,i])[2]
  
  # R: concatenate the results (store them in the nearest_n object)
  nearest_n <- c(nearest_n, nn_i)
}

# The resulting object is now a vector containing the 
# indices of the nearest neighbors of each of the 247 samples 
# in xr_data (e.g. the nearest neighbor of the first sample in xr_data
# is the sample number 181)
nearest_n

# Alternatively you can use an internal function of the resemble package 
# to identify the nearest neighbors from a squared symetric matrix of 
# dissimilarities. In this case a matrix with one column is returned 
nearest_n2 <- resemble:::which_min(X = md_xr, cores = 1)

nearest_n2

# Step 3. Compare the side information of each sample in X to the side 
# information of its corresponding closest sample. Evaluate the statistics 
# of these comparisons for assessing the reliability of the distance metric
# algorithm

# Let's say that the side information is the clay content ("Clay" column in xr_data)
# These are the clay content values of each sample in xr_data
xr_data[,"Clay"]

# R: Get the clay content values of the nearest neighbors found in xr_data for 
# each sample in xr_data
# This is done by using the nearest_n vector which contains the indices 
# of the nearest neighbors. The xr_data can be reordered according to these indices
xr_data[nearest_n,"Clay"]

# You can compare the samples and its nearest neighbors 
# in terms of their clay content values. For doing so, you can use the root mean
# square of differences (RMSD)

# R: first compute the vector of squared differences...
sqdiff <- (xr_data[,"Clay"] - xr_data[nearest_n,"Clay"])^2

# R: ... then compute the average...
mean_sqdiff <- mean(sqdiff)

# R: ... and finally compute the RMSD
rmsd <- mean_sqdiff^0.5

# R: You can also compute the correlation coefficient between the
# the clay content values of the samples and the caly contetnt 
# values of their nearest neighbors 
r <- cor(xr_data[,"Clay"], xr_data[nearest_n,"Clay"])

# You can also use a graphical comparison of the side 
# information (clay content values) of the samples and 
# their nearest neighbors 
plot(x = xr_data[,"Clay"], 
     y = xr_data[nearest_n,"Clay"], 
     xlab = "Clay content of the samples, %",
     ylab = "Clay content of the nearest neighbors, %",
     pch = 16,
     col = rgb(red = 1, green = 0.2, blue = 0.2, alpha = 0.5))
grid()

# A GOOD SPECTRAL DISSIMILARITY MEASURE SHOULD ALSO 
# BE ABLE TO REPRESENT THE COMPOSIIONAL DISSIMILARITY!
# To select a given spectral dissimilarity matrix (e.g. computed by different methods)
# you can choose the one that maximizes the compositional similarity between samples
# (e.g. that minimizes the RMSD)


# The simEval fucntion from the resemble package can be used to asses the 
# reliability of the dissimilarity matrices

# R: Use the simEval function to get the RMSD and the r
d_eval <- simEval(d = md_xr, sideInf = xr_data$Clay)

# R: get a summary of the d_eval object created with the simEval fucntion 
# (Compactly Display the Structure of the object)
str(d_eval)

# R: get just the names of the sub-objects in the simEval object created 
names(d_eval)

# R: get the evaluation results
d_eval$eval

# R: get the results of the nearest neighbor search 
d_eval$firstNN

#---- 3. Other examples with the Mahalanobis distance computed on the PCs ----

# This approach can be used to select the optimal number of PCs for dissimilarity computations
# You can evaluate for example a series of dissimilarity matrices computed from
# a sequence of different number of PCs retained

# R: define the sequence of the total number of PCs to be retained
# (e.g. from 2 to 20)
max_pcs <- 2:20

# USER: First create an empty object to store the results of the evaluation 
# results for each dissimilarity matrix computed on the diffrent number of PCs
eval_results <- NULL

# R: use a for loop to iterate over each of the columns in the 
# md_xr matrix
for(i in max_pcs){
  # R: Temporarily store the pcProjection results of the ith iteration 
  # in an object (i_pcspectra)
  i_pcspectra <- pcProjection(Xr = comb_x, 
                              pcSelection = list("manual", i), 
                              center = TRUE, scaled = FALSE)
  
  # R: Temporarily store the fDiss results of the ith iteration 
  # in an object (i_md_xr)
  i_md_xr <- fDiss(Xr =  i_pcspectra$scores[1:nrow(xr_data), ], 
                   X2 =  i_pcspectra$scores[1:nrow(xr_data), ], 
                   method = "mahalanobis", 
                   center = FALSE, scaled = FALSE)
  
  # R: Temporarily store the results of the simEval of the ith iteration 
  # in an object (nn_i)
  i_eval <- simEval(d = i_md_xr, sideInf = xr_data$Clay)

  # R: concatenate the results by rows (store them in the eval_results object)  
  eval_results <- rbind(eval_results, i_eval$eval)
}
  
# Check the results (in this case the optimal number of PCs is 7)
eval_results

# R: plot the RMSD results...
# USER: before plotting provide names of the axes 
# of the plot
plot(x = max_pcs, 
     y = eval_results$rmsd,
     xlab = "Principal component",
     ylab = "RMSD",
     type = "p",
     pch = 16,
     col = rgb(red = 0, green = 0.4, blue = 0.8, alpha = 0.5))
grid()



# The above code/procedures to select the optimal number of PCs 
# is summarized in one single function in the resemble package
# The orthoDiss function offers the possibility to compute the 
# Mahalanobis distance in the PC space by retaining the optimal 
# number of PCs that minimize the RMSD

# Note that the side information is provided in the Yr argument of this 
# function (which is the side information corresponding to the spectra in Xr)

# R: compute the Mahalanobis distance in the PC space suing the 
# side information approach ("opc" method which is specified in the pcSelection 
# argument). In this case we evaluate a sequence of dissimilarity matrices
# from 2 to 20
opcd <- orthoDiss(Xr = xr_data$spc_snv, 
                  X2 = xu_data$spc_snv, 
                  Yr = xr_data$Clay, 
                  pcSelection = list("opc", 20), 
                  method = "pca", 
                  center = TRUE, scaled = FALSE)

# R: get a summary of the opcd object created with the orthoDiss fucntion 
# (Compactly Display the Structure of the object)
str(opcd)

# R: get just the names of the sub-objects in the opcd object created 
names(opcd)

# R: get the optimal numer of PCs retained
opcd$n.components

# R: get the dissimilarity matrix
opcd$dissimilarity
