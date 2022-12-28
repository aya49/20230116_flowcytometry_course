source("utils.R")

# packages
install_pkgs(c(
  'devtools', # helps install packages
  'pracma', # practical numerical math functions
  'stringr', # nice string manipulation functions
  'umap', # UMAP dimensionality reduction
  'plyr', # data frame handling functions
  'dplyr'
), repo="cran")

install_pkgs(c(
  # flow cytometry data type handling packages
  'flowCore',
  'flowWorkspace', # devtools::install_github("RGLab/flowWorkspace", ref="trunk")
  'flowDensity', # Bioc; automated gating and plotting
  'FlowSOM',
  'PeacoQC' # Bioc; QA
), repo="bioc")
