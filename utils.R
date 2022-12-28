## input: a list of packages
## output: installs given packages
install_pkgs_bioc <- function(pkgs2) {
  # install bioconductor packages
  pkgs3 <- setdiff(pkgs2, installed.packages())
  if (length(pkgs3) > 0) {
    devtools::install_bioc(pkgs3, verbose=F)
  }
}
install_pkgs_cran <- function(pkgs1) {
  # install cran packages
  pkgs2 <- setdiff(pkgs1, installed.packages())
  if (length(pkgs2) > 0) {
    devtools::install_cran(pkgs2, verbose=F)
  }
}
install_pkgs <- function(pkgs, repo=NULL) {
  # install devtools to install packages from github
  if (!"devtools"%in%installed.packages()) {
    install.packages("devtools", verbose=F)
  }
  
  # install github packages
  pkgsgh <- grepl("[/]", pkgs)
  if (any(pkgsgh)) {
    pkgsghp <- pkgs[pkgsgh]
    pkgsghn <- unlist( lapply(strsplit(pkgsghp, "/"), function(x) x[length(x)]) )
    pkgs1 <- pkgsghp[!pkgsghn %in% installed.packages()]
    if (length(pkgs1) > 0) {
      for (pkgs1_ in pkgs1) {
        devtools::install_github(pkgs1, verbose=F)
      }
    }
    suppressWarnings(sapply(pkgsghn, require, character.only=TRUE))
  }
  
  if (is.null(repo)) {
    install_pkgs_cran(pkgs[!pkgsgh])
    install_pkgs_bioc(pkgs[!pkgsgh])
  } else if (repo=="bioc") {
    install_pkgs_bioc(pkgs[!pkgsgh])
  } else if (repo=="cran") {
    install_pkgs_cran(pkgs[!pkgsgh])
  }
  
  suppressWarnings(sapply(pkgs[!pkgsgh], require, character.only=TRUE))
}

## input: start time and message
## output: duration since start time
time_output <- function(start, msg="") {
  start <- as.POSIXct(start)
  end <- Sys.time()
  time_elapsed <- difftime(end, start, units="secs")
  print(paste0(msg, ifelse(msg == "", "", ": "),
               tstr(start), "-", tstr(end), " > ", tstr(time_elapsed)))
}
tstr <- function(time) format(.POSIXct(time), "%H:%M:%S")

## input: exprs matrix
## output: density plot
denscolour <- function(exprs,coli,colp) {
  if (missing(coli)) {
    dat <- exprs[,1:2]
  } else {
    dat <- exprs[,coli]
  }
  if (missing(colp))
    colp <- colorRampPalette(c("blue", "turquoise", 
                               "green", "yellow", "orange", "red"))
  require(grDevices)
  col <- densCols(dat, colramp=colp) #grDevices
  return(col)
}

## from the brinkmanlab
## input: 2D matrix
## output: 2D density + colours for each point
densCols <- function(x, y=NULL, nbin=c(750,750), bandwidth=NULL, colramp=colorRampPalette(c("blue", "turquoise", "green", "yellow", "orange", "red"))) {
  xy <- xy.coords(x, y, setLab=FALSE)
  select <- is.finite(xy$x) & is.finite(xy$y)
  X <- cbind(xy$x, xy$y)[select, ]
  
  if (is.null(bandwidth)) {
    bandwidth <- diff(apply(X, 2, stats::quantile, probs=c(0.05, 0.95), na.rm=TRUE, names=FALSE))/25
    bandwidth[bandwidth == 0] <- 1
  }
  if (!is.numeric(bandwidth)) stop("'bandwidth' must be numeric")
  if (any(bandwidth <= 0)) stop("'bandwidth' must be positive")
  
  map <- KernSmooth::bkde2D(X, bandwidth=bandwidth, gridsize=nbin)
  mkBreaks <- function(u) u - diff(range(u))/(length(u) - 1)/2
  xbin <- cut(X[, 1], mkBreaks(map$x1), labels=FALSE)
  ybin <- cut(X[, 2], mkBreaks(map$x2), labels=FALSE)
  dens <- map$fhat[cbind(xbin, ybin)]
  dens[is.na(dens)] <- 0
  colpal <- cut(dens, length(dens), labels=FALSE)
  cols <- rep(NA_character_, length(select))
  cols[select] <- colramp(length(dens))[colpal]
  
  return(list(cols=cols, map=map))
}
