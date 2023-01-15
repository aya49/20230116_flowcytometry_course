## preprocessing fcs file
## author: alice yue
## input: raw fcs file
## output: preprocessed fcs file

# load packages
library("flowCore")
library("PeacoQC")

# path to fcs file 
fcs_path <- system.file("extdata", "111.fcs", package="PeacoQC")
# fcs_path <- "/home/user/folder/file.fcs"

# load fcs file
f <- flowCore::read.FCS(fcs_path)

# explore fcs file
# flowCore::exprs(f) # raw cell x marker matrix
head(flowCore::exprs(f))
dim(flowCore::exprs(f))
flowCore::markernames(f) # marker names (excluding morphology columns)
f@parameters@data


## 1.1 compensate ####
spillover <- flowCore::keyword(f)$SPILL
f <- flowCore::compensate(f, spillover=spillover)

## 1.2 logicle transform ####
transformList <- flowCore::estimateLogicle(f, channels=colnames(spillover))
f <- flowWorkspace::transform(f, transformList)

## 1.3 cleaning; see res_dir for plot ####
## parameter: Mean Absolute Deviation (MAD) distance (decrease = less strict)
## parameter: IsolationTree (IT) (decrease = more strict)
channels <- c(1, 3, 5:14, 18, 21)
res_dir <- "/home/maruko/projects/temp" # where to save PeacoQC plot
fmr <- PeacoQC::RemoveMargins(f, channels=channels, output="full")
pQC <- PeacoQC::PeacoQC(fmr[["flowframe"]], channels=channels,
                        plot=TRUE, save_fcs=FALSE, report=FALSE,
                        output_directory=res_dir)

# final preprocessed fcs!
fc <- pQC[["FinalFF"]] 

# clean memory by removing variables
rm("fmr") # remove variable
rm("pQC")
gc() # remove removed variables from memory