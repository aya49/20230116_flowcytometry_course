## gate fcs file with clustering
## author: alice yue
## input: raw fcs file
## output: gated flowworkspace

# load packages
library("flowCore")
library("flowWorkspace") # also imports ggcyto and ggplot2 packages
library("Rphenograph") # also imports igraph package
library("FlowSOM")

# directory to save results in
res_dir <- "/home/maruko/projects/gating"

# set seed for randomness
set.seed(4)

# load gating set
gs <- flowWorkspace::load_gs(paste0(res_dir, "/gs"))
cpops <- flowWorkspace::gs_get_leaf_nodes(gs)

# prepare colours for plots!
cpopcol <- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", 
             "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")


## 2 GATING: cell population identification based on a gating strategy ####

## 2.15 gating b2bcells > preB/MZB/folB ####
# CD21: high > MZ > fol > pre > low
ff <- flowWorkspace::cytoframe_to_flowFrame(
    as(flowWorkspace::gs_pop_get_data(gs, "b2bcells"), Class="list")[[1]] )

# get flowframe > 2D frame
ff2D <- flowCore::exprs(ff)[,c("BV421-A","PE-A")]

# specify number of clusters
n_clust <- 3

# try phenograph
pc_ <- Rphenograph::Rphenograph(ff2D, 50)
pc <- as.numeric(igraph::membership(pc_[[2]]))

# consolidate phenograph clusters with k-medoid
clustu <- unique(pc)
clustmeds <- matrix(0, ncol=2, nrow=length(clustu))
colnames(clustmeds) <- c("BV421-A", "PE-A")
for (i in clustu) {
    clusti <- pc==i
    clustmeds[i,] <- c(median(ff2D[clusti,1]), median(ff2D[clusti,2]))
}
pkc_ <- cluster::pam(clustmeds, k=n_clust)
pkc <- pkc_$clustering[pc]

# try flowsom
fc_ <- FlowSOM::FlowSOM(ff2D, nClus=n_clust)
fc <- as.numeric(FlowSOM::GetMetaclusters(fc_))


## plot! ####
ff2Dd <- data.frame(
    x=ff2D[,1], y=ff2D[,2], 
    phenograph=pc, phenograph_kmedoids=pkc, flowsom=fc)

gp2p <- ggplot2::ggplot(ff2Dd, ggplot2::aes(x=x, y=y, colour=factor(phenograph))) + ggplot2::geom_point(size=0.1)
gp2k <- ggplot2::ggplot(ff2Dd, ggplot2::aes(x=x, y=y, colour=factor(phenograph_kmedoids))) + ggplot2::geom_point(size=0.1)
gp2f <- ggplot2::ggplot(ff2Dd, ggplot2::aes(x=x, y=y, colour=factor(flowsom))) + ggplot2::geom_point(size=0.1)

gp2p
gp2k
gp2f

# compare with original gating
ggcyto::autoplot(gs[[1]], "folbcells", bins=100)


# # nonconvex hull function outdated: devtools::install_github("andrewzm/INLA")
# # register gates into gs
# cluster <- 1 # choose a cluster
# nchull <- INLA::inla.nonconvex.hull(ff2D[clusts==cluster,, drop=FALSE])
# gate <- flowCore::polygonGate(.gate=nchull)
# node <- flowWorkspace::gs_pop_add(gs, gate, name="folbcells", parent="b2bcells")
# flowWorkspace::recompute(gs)


## save gating set! ####
# flowWorkspace::save_gs(gs, path=paste0(res_dir, "/gs2"))
# gs <- flowWorkspace::load_gs(paste0(res_dir, "/gs"))
