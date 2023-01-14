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

## rotate 2D frame
## input: 2D matrix and angle
## output: rotated 2D matrix
rotate_data <- function(data, theta=pi/2 - atan(lm(data.new[,1] ~ data.new[,2])$coefficients[2])) {
  data %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2,byrow=T)
}