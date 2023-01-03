## gate fcs file
## author: alice yue
## input: raw fcs file
## output: gated flowworkspace

# path to raw fcs file
fcs_path <- "/home/maruko/projects/sangerP2.fcs"

## load fcs file
f <- flowCore::read.FCS(fcs_path)

## 1.1 compensate
spillover <- flowCore::keyword(f)$SPILL
f <- flowCore::compensate(f, spillover=spillover)

## 1.2 logicle transform
transformList <- flowCore::estimateLogicle(f, channels=colnames(spillover))
f <- flowWorkspace::transform(f, transformList)

## 1.3 cleaning; see res_path for plot
fmr <- PeacoQC::RemoveMargins(f, channels=c("FSC-A", "SSC-A"), output="full")
pQC <- PeacoQC::PeacoQC(fmr[["flowframe"]], channels=colnames(spillover),
                        plot=FALSE, save_fcs=FALSE, report=FALSE)
fc <- pQC[["FinalFF"]]

# clean memory by removing variables
rm("fmr")
rm("pQC")
gc()


## 2 GATING: cell population identification based on a gating strategy ####

# initialize gating set (containing 1 file)
fs <- flowCore::flowSet(list(f))
# fs <- flowCore::flowSet(list(f1, f2, f3)) # can add more files if needed
gs <- flowWorkspace::GatingSet(fs)

## 2.1 gating singlets ####

# get threshold gates
gate_singlets <- flowDensity::deGate(f, channel=c("SSC-A"), tinypeak.removal=0.95, upper=TRUE, use.upper=TRUE, )
gate_singlets_low <- flowDensity::deGate(f, channel=c("SSC-A"), use.percentile=TRUE, percentile=0.0001)

# gate
singlets_high <- flowDensity::flowDensity(
    f, channels=c("FSC-A", "SSC-A"), 
    position=c(NA, F), gates=c(NA, gate_singlets))
singlets <- flowDensity::flowDensity(
    singlets_high, channels=c("FSC-A", "SSC-A"), 
    position=c(NA, T), gates=c(NA, gate_singlets_low))

# plot
flowDensity::plotDens(f, channels=c("FSC-A","SSC-A"), main="All", xlim=c(0, 275000), ylim=c(0, 275000))
lines(singlets@filter)
abline(h=c(gate_singlets, gate_singlets_low), lty="dashed", col="red")
