## preprocessing fcs file
## author: alice yue
## input: raw fcs file (multiple)
## output: preprocessed fcs file

# specify path to folder containing all .fcs files you want to analyze
folder_path <- "path/to/folder"

# list .fcs files in folder
fcs_paths <- list.files(folder_path, full.names=TRUE)
fcs_paths <- fcs_paths[grepl("[.]fcs$", fcs_paths, ignore.case=TRUE)]

# old and new folder name
folder_split <- strsplit(folder_path, "/")
folder_old <- folder_split[[1]][length(folder_split[[1]])]
folder_clean <- "clean"

# loop through each .fcs file path
for (fcs_path in fcs_paths) {
    
    # load fcs file
    f <- flowCore::read.FCS(fcs_path)
    
    ## 1.1 compensate ####
    spillover <- flowCore::keyword(f)$SPILL
    # spillover <- read.csv("path", header=TRUE, row.names=TRUE)
    # flowCore::keyword(f)$SPILL <- spillover
    f <- flowCore::compensate(f, spillover=spillover)
    
    ## 1.2 logicle transform ####
    transformList <- flowCore::estimateLogicle(f, channels=colnames(spillover))
    f <- flowWorkspace::transform(f, transformList)
    
    # let's look at the Time vs FSC-A plot to see the flow of cells
    flowDensity::plotDens(f, channels=c("Time", "FSC-A"))
    
    ## 1.3 cleaning; see res_dir for plot ####
    ## parameter: Mean Absolute Deviation (MAD) distance (decrease = less strict)
    ## parameter: IsolationTree (IT) (decrease = more strict)
    fmr <- PeacoQC::RemoveMargins(f, channels=grep("[FS]SC", colnames(flowCore::exprs(f))), output="full")
    pQC <- PeacoQC::PeacoQC(fmr[["flowframe"]], channels=colnames(spillover),
                            plot=FALSE, save_fcs=FALSE, report=FALSE)
    
    fcs_path_new <- gsub(folder_old, folder_clean, fcs_path)
    flowCore::write.FCS(f, fcs_path_new)
}