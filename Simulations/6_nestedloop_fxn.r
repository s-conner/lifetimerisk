nestedloop <- function(x,
                varnames, sign=rep(1, length(varnames)),
                varlabels=NULL){
  ##
  if (!inherits(x, "data.frame"))
    stop("Argument 'x' must be a data.frame.")
  ##
  mo <- matrix(sign,
               nrow=dim(x)[[1]], ncol=length(varnames),
               byrow=TRUE)
  xo <- x[,varnames]
  ##
  ## Re-ordering:
  res <- x[do.call(order, mo*xo),]
  ##
  attr(res, "varnames")  <- varnames
  attr(res, "varlabels") <- varlabels
  ##
  class(res) <- c("nestedloop", class(res))
  ##
  res
}



lines.nestedloop <- function(x,
                             varnames=attributes(x)$varnames,
                             varlabels=attributes(x)$varlabels,
                             which="v",
                             col=if (which=="r") "#999999" else "black",
                             ymin.refline, ymax.refline,
                             cex.ref=0.9,
                             log=TRUE,
                             ...){
  ##
  nvar <- length(varnames)
  ##
  if (length(col)==1)
    col <- rep(col, nvar)
  ##
  if (which=="v"){
    ##
    ## Vertical lines
    ##
    nlen <- rep(NA, nvar)
    ##
    for (i in 1:nvar)
      nlen[i] <- length(unique(x[,varnames[i]]))
    ##
    cnlen <- cumprod(nlen)
    ##
    for (i in (nvar-1):1)
      abline(v=cnlen[nvar]*(0:cnlen[i])/cnlen[i]+1,
             col=col[i])
  }
  else if (which=="r"){
    ##
    ## Reference lines
    ##
    if (is.null(varlabels))
      varlabels <- varnames
    ##
    labels.varnames <- rep("", nvar)
    ##
    for (i in 1:length(varnames)){
      if (is.factor(x[,varnames[i]]))
        varvals <- unique(x[,varnames[i]])
      else{
        varvals <- format(unique(x[,varnames[i]]))
        varvals <- sub("^[[:space:]]*(.*?)[[:space:]]*$",
                       "\\1",
                       varvals,
                       perl=TRUE)
      }
      ##
      labels.varnames[i] <- paste(varlabels[i],
                                  " (",
                                  paste(varvals,
                                        collapse=", "),
                                  ")", sep="")
    }
    ##
    if (log){
      ymax <- log(ymax.refline)
      ymin <- log(ymin.refline)
    }
    else{
      ymax <- ymax.refline
      ymin <- ymin.refline
    }
    ##
    distance <- (ymax-ymin)/nvar
    ##
    ypos <- ymax-0.2*distance-(1/nvar)*(0:(nvar-1))*(ymax-ymin)
    ypos.ref.max <- ypos-0.20*distance
    ypos.ref.min <- ypos-0.75*distance
    ##
    if (log){
      ypos <- exp(ypos)
      ypos.ref.max <- exp(ypos.ref.max)
      ypos.ref.min <- exp(ypos.ref.min)
    }
    ##
    for (i in 1:nvar){
      ##print(c(ypos.ref.max[i], ypos.ref.min[i]))
      ##
      ##abline(h=ypos[i], lwd=1, col="red")
      ##abline(h=ypos.ref.max[i], lwd=1, col="blue")
      ##abline(h=ypos.ref.min[i], lwd=1, col="green")
      text(1, ypos[i], labels.varnames[i], adj=0, cex=cex.ref)
      ##
      xvar <- x[,varnames[i]]
      if (is.factor(xvar))
          xvar <- as.numeric(xvar)
      xvar <- ypos.ref.min[i] +
        (xvar-min(xvar))/(max(xvar)-min(xvar))*(ypos.ref.max[i]-ypos.ref.min[i])
      
      xvar <- c(xvar, tail(xvar,1))
      lines(xvar, col=col[i], type="s", lwd=1)
    }
  }
  ##
  invisible(NULL)
}


panel.nestedloop <- function(x,
                             varnames=attributes(x)$varnames,
                             varlabels=attributes(x)$varlabels,
                             which="v",
                             col=if (which=="r") "#999999" else "black",
                             ymin.refline, ymax.refline,
                             cex.ref=0.9,
                             log=TRUE,
                             ...){
  ##
  nvar <- length(varnames)
  ##
  if (length(col)==1)
    col <- rep(col, nvar)
  ##
  if (which=="v"){
    ##
    ## Vertical lines
    ##
    nlen <- rep(NA, nvar)
    ##
    for (i in 1:nvar)
      nlen[i] <- length(unique(x[,varnames[i]]))
    ##
    cnlen <- cumprod(nlen)
    ##
    for (i in (nvar-1):1)
      abline(v=cnlen[nvar]*(0:cnlen[i])/cnlen[i]+1,
             col=col[i])
  }
  else if (which=="r"){
    ##
    ## Reference lines
    ##
    if (is.null(varlabels))
      varlabels <- varnames
    ##
    labels.varnames <- rep("", nvar)
    ##
    for (i in 1:length(varnames)){
      if (is.factor(x[,varnames[i]]))
        varvals <- unique(x[,varnames[i]])
      else{
        varvals <- format(unique(x[,varnames[i]]))
        varvals <- sub("^[[:space:]]*(.*?)[[:space:]]*$",
                       "\\1",
                       varvals,
                       perl=TRUE)
      }
      ##
      labels.varnames[i] <- paste(varlabels[i],
                                  " (",
                                  paste(varvals,
                                        collapse=", "),
                                  ")", sep="")
    }
    ##
    if (log){
      ymax <- log(ymax.refline)
      ymin <- log(ymin.refline)
    }
    else{
      ymax <- ymax.refline
      ymin <- ymin.refline
    }
    ##
    distance <- (ymax-ymin)/nvar
    ##
    ypos <- ymax-0.2*distance-(1/nvar)*(0:(nvar-1))*(ymax-ymin)
    ypos.ref.max <- ypos-0.20*distance
    ypos.ref.min <- ypos-0.75*distance
    ##
    if (log){
      ypos <- exp(ypos)
      ypos.ref.max <- exp(ypos.ref.max)
      ypos.ref.min <- exp(ypos.ref.min)
    }
    ##
    for (i in 1:nvar){
      ##print(c(ypos.ref.max[i], ypos.ref.min[i]))
      ##
      ##abline(h=ypos[i], lwd=1, col="red")
      ##abline(h=ypos.ref.max[i], lwd=1, col="blue")
      ##abline(h=ypos.ref.min[i], lwd=1, col="green")
      ltext(1, ypos[i], labels.varnames[i], adj=c(0, 0.5), cex=cex.ref)
      ##
      xvar <- x[,varnames[i]]
      if (is.factor(xvar))
          xvar <- as.numeric(xvar)
      xvar <- ypos.ref.min[i] +
        (xvar-min(xvar))/(max(xvar)-min(xvar))*
          (ypos.ref.max[i]-ypos.ref.min[i])
      
      xvar <- c(xvar, tail(xvar,1))
      llines(xvar, col=col[i], type="s", lwd=1)
    }
  }
  ##
  invisible(NULL)
}
