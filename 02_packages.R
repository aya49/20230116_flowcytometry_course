# install generic packages from CRAN using the install.packages() function.
install.packages("devtools") # helps install packages
install.packages("stringr") # nice string manipulation functions
install.packages("umap") # UMAP dimensionality reduction function
install.packages("plyr") # data handling functions
install.packages("dplyr") # more data handling functions

# if the BiocManager package has not been installed, install it.
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# install flow cytometry (and other biostatistic) packages from Bioconductor 
# using the "install" function from the BiocManager package.
BiocManager::install("flowCore")
BiocManager::install("flowWorkspace")
BiocManager::install("PeacoQC")
BiocManager::install("flowDensity")
BiocManager::install("FlowSOM")