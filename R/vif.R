# Author: Babak Naimi, naimi.b@gmail.com
# Date :  Oct. 2019
# Last update: September 2023
# Version 1.7
# Licence GPL v3

.vif <- function(.dd) {
  z<-rep(NA,ncol(.dd))
  names(z) <- colnames(.dd)
  for (i in 1:ncol(.dd)) {
    z[i] <-  1 / (1 - summary(lm(.dd[,i]~.,data=.dd[-i]))$r.squared)
  }
  return(z)
}


.vif2 <- function(y,w) {
  z<-rep(NA,length(w))
  names(z) <- colnames(y)[w]
  for (i in 1:length(w)) {
    z[i] <-  1/(1-summary(lm(as.formula(paste(colnames(y)[w[i]],"~.",sep='')),data=y))$r.squared)
  }
  return(z)
}

.maxCor <- function(k){
  k <- abs(k)
  n <- nrow(k)
  for (i in 1:n) k[i:n,i] <- NA
  w <- which.max(k)
  c(rownames(k)[((w%/%nrow(k))+1)],colnames(k)[w%%nrow(k)])
}
#---
# it checks the correlation and if both of the variables in the pair with max-Cor are in keep, it 
# goes to next order of maximum correlation until a pair with at least one of them not in keep is selected
.maxCor2 <- function(k,keep){
  k <- abs(k)
  n <- nrow(k)
  for (i in 1:n) k[i:n,i] <- NA
  #w <- which.max(k)
  o <- order(k,decreasing = TRUE)
  w <- o[1]
  vn <- c(rownames(k)[((w%/%nrow(k))+1)],colnames(k)[w%%nrow(k)])
  if (all(vn %in% keep)) {
    LOOP <-TRUE
    j <- 2
    while(LOOP) {
      w <- o[j]
      vn <- c(rownames(k)[((w%/%nrow(k))+1)],colnames(k)[w%%nrow(k)])
      if (all(vn %in% keep) && j < length(o)) j <- j+1
      else LOOP <- FALSE
    }
    return (list(variables=vn,warning_level=j))
  } else {
    return(vn)
  }
  
}


.minCor <- function(k){
  k <- abs(k)
  rr<-c();cc<-c();co<-c()
  for (c in 1:(ncol(k)-1)) {
    for (r in (c+1):nrow(k)){
      rr<-c(rr,rownames(k)[r]);cc<-c(cc,colnames(k)[c])
      co <- c(co,k[r,c])
    }
  }
  w <- which.min(co)
  c(rr[w],cc[w])
}
#-----------
# Builds the transitive exclusion chains from an exclusionLog data.frame.
# Returns a named list: each name is a final retained variable; the value is a
# character vector of excluded predictors attributed to it, ordered from the
# most recently excluded (direct pairing partner) to the most distally excluded.
# Example: list(Bio7 = c("Bio10", "Bio5")) means Bio7 beat Bio10 which had
# previously beaten Bio5, so the chain reads Bio7 <- Bio10 <- Bio5.
.buildChains <- function(log) {
  if (is.null(log) || nrow(log) == 0L) return(list())
  accumulated <- list()
  for (i in seq_len(nrow(log))) {
    ex <- log$excluded[i]
    kp <- log$correlated_with[i]
    # Collect the chain of ex: ex itself plus anything ex had previously won
    ex_chain <- if (!is.null(accumulated[[ex]])) c(ex, accumulated[[ex]]) else ex
    # ex is now truly removed; clear its accumulated record
    accumulated[[ex]] <- NULL
    # Transfer ex_chain to the winner (kp)
    if (is.null(accumulated[[kp]])) {
      accumulated[[kp]] <- ex_chain
    } else {
      accumulated[[kp]] <- c(accumulated[[kp]], ex_chain)
    }
  }
  accumulated
}
#-----------------
if (!isGeneric("vif")) {
  setGeneric("vif", function(x,size, ...)
    standardGeneric("vif"))
}  
setMethod('vif', signature(x='RasterStackBrick'),
          function(x, size) {
            if (nlayers(x) == 1) stop("The Raster object should have at least two layers")
            if (missing(size)) {
              if (ncell(x) < 7000) size <- NULL
              else size <- 5000
            }
            
            if (is.null(size)) x <- as.data.frame(x,na.rm=TRUE)
            else x <- sampleRandom(x,size,na.rm=TRUE)
            
            x <- na.omit(x)
            v <- .vif(x)
            
            data.frame(Variables=names(v),VIF=as.vector(v))
          }
)
#-----
setMethod('vif', signature(x='SpatRaster'),
          function(x, size) {
            if (nlyr(x) == 1) stop("The Raster object should have at least two layers")
            if (missing(size)) {
              if (ncell(x) < 7000) size <- NULL
              else size <- 5000
            }
            
            if (is.null(size)) x <- as.data.frame(x,na.rm=TRUE)
            else x <- spatSample(x,size,na.rm=TRUE)
            
            x <- na.omit(x)
            v <- .vif(x)
            
            data.frame(Variables=names(v),VIF=as.vector(v))
          }
)
setMethod('vif', signature(x='data.frame'),
          function(x, size) {
            if (ncol(x) == 1) stop("At least two variables are needed to quantify vif")
            x <- na.omit(x)
            
            if (missing(size)) {
              if (nrow(x) < 7000) size <- NULL
              else size <- 5000
            }
            if (!is.null(size)) {
              if(nrow(x) > size) x <- x[sample(1:nrow(x),size),]
            }
            v <- .vif(x)
            data.frame(Variables=names(v),VIF=as.vector(v))
          }
)

setMethod('vif', signature(x='matrix'),
          function(x, size) {
            if (ncol(x) == 1) stop("At least two variables are needed to quantify vif")
            
            x <- na.omit(x)
            
            if (missing(size)) {
              if (nrow(x) < 7000) size <- NULL
              else size <- 5000
            }
            if (!is.null(size)) {
              if(nrow(x) > size) x <- x[sample(1:nrow(x),size),]
            }
            
            v <- .vif(x)
            data.frame(Variables=names(v),VIF=as.vector(v))
          }
)

if (!isGeneric("vifcor")) {
  setGeneric("vifcor", function(x, th= 0.9,keep=NULL, size, method='pearson',...)
    standardGeneric("vifcor"))
}  
setMethod('vifcor', signature(x='RasterStackBrick'),
          function(x, th=0.9, keep=NULL, size, method='pearson') {
            if (nlayers(x) == 1) stop("The Raster object should have at least two layers")
            if (missing(method) || !method %in% c('pearson','kendall','spearman')) method <- 'pearson'
            
            if (missing(size)) {
              if (ncell(x) < 7000) size <- NULL
              else size <- 5000
            }
            
            if (missing(keep)) keep <- NULL
            #-----------
            
            if (is.null(size)) x <- as.data.frame(x,na.rm=TRUE)
            else x <- sampleRandom(x,size,na.rm=TRUE)
            
            x <- na.omit(x)
            
            vn <- colnames(x)
            
            if (!is.null(keep)) {
              if (is.numeric(keep)) {
                .w <- keep %in% c(1:length(vn))
                if (all(!.w)) stop(paste0('the values in keep are out of range; should be between 1 and ',ncol(x),' for your dataset...!'))
                if (any(!.w)) {
                  keep <- keep[.w]
                  warning('some of layer numbers in keep are out of range (should be between 1 and nlyr(x)), so they are ignored!')
                }
                keep <- vn[keep]
              } else if (is.character(keep)) {
                .w <- keep %in% colnames(x)
                if (all(!.w)) stop('None of the variable names in keep are available in x dataset!')
                if (any(!.w)) {
                  keep <- keep[.w]
                  warning('some of the variable names in keep are not available in x (check if it is typos), so they are ignored!')
                }
              }
            }
            #----------------
            .warn <- FALSE
            
            LOOP <- TRUE
            n <- new("VIF")
            n@variables <- colnames(x)
            exc <- c()
            excLog <- data.frame(step=integer(), excluded=character(), correlated_with=character(),
                                 correlation=numeric(), vif_excluded=numeric(), vif_retained=numeric(),
                                 stringsAsFactors=FALSE)
            step_i <- 0L
            if (is.null(keep)) {
              while (LOOP) {
                xcor <- abs(cor(x, method=method))
                mx <- .maxCor(xcor)
                if (xcor[mx[1],mx[2]] >= th) {
                  step_i <- step_i + 1L
                  w1 <- which(colnames(xcor) == mx[1])
                  w2 <- which(rownames(xcor) == mx[2])
                  v <- .vif2(x,c(w1,w2))
                  ex <- mx[which.max(v[mx])]
                  kp <- mx[mx != ex]
                  excLog <- rbind(excLog, data.frame(step=step_i, excluded=ex, correlated_with=kp,
                                   correlation=xcor[mx[1],mx[2]], vif_excluded=v[ex], vif_retained=v[kp],
                                   stringsAsFactors=FALSE))
                  exc <- c(exc,ex)
                  x <- x[,-which(colnames(x) == ex)]
                } else LOOP <- FALSE
              }
            } else {
              
              while (LOOP) {
                xcor <- abs(cor(x, method=method))
                mx <- .maxCor2(xcor,keep=keep)
                if (is.list(mx)) {
                  .warn <- TRUE
                  mx <- mx[[1]]
                }
                if (xcor[mx[1],mx[2]] >= th) {
                  step_i <- step_i + 1L
                  w1 <- which(colnames(xcor) == mx[1])
                  w2 <- which(rownames(xcor) == mx[2])
                  if (any(mx %in% keep)) {
                    ex <- mx[!mx %in% keep]
                    if (length(ex) == 1L) {
                      kp <- mx[mx != ex]
                      v_pair <- .vif2(x, c(w1, w2))
                      excLog <- rbind(excLog, data.frame(step=step_i, excluded=ex, correlated_with=kp,
                                       correlation=xcor[mx[1],mx[2]], vif_excluded=v_pair[ex], vif_retained=v_pair[kp],
                                       stringsAsFactors=FALSE))
                    }
                  } else {
                    v <- .vif2(x,c(w1,w2))
                    ex <- mx[which.max(v[mx])]
                    kp <- mx[mx != ex]
                    excLog <- rbind(excLog, data.frame(step=step_i, excluded=ex, correlated_with=kp,
                                     correlation=xcor[mx[1],mx[2]], vif_excluded=v[ex], vif_retained=v[kp],
                                     stringsAsFactors=FALSE))
                  }
                  
                  exc <- c(exc,ex)
                  x <- x[,-which(colnames(x) == ex)]
                } else LOOP <- FALSE
              }
            }
            
            if (.warn) {
              if (length(keep) > 2) warning('At least two variables specified in "keep" are strongly correlated (i.e., are subjected to the collinearity issue)!')
              else warning('The two variables specified in "keep" are strongly correlated (i.e., are subjected to the collinearity issue)!')
            }
            #---
            if (length(exc) > 0) n@excluded <- exc
            n@exclusionLog <- excLog
            n@chains <- .buildChains(excLog)
            v <- .vif(x)
            n@corMatrix <- cor(x, method=method)
            n@results <- data.frame(Variables=names(v),VIF=as.vector(v))
            n
          }
)
#-------
setMethod('vifcor', signature(x='SpatRaster'),
          function(x, th=0.9, keep=NULL, size, method='pearson') {
            if (nlyr(x) == 1) stop("The Raster object should have at least two layers")
            if (missing(method) || !method %in% c('pearson','kendall','spearman')) method <- 'pearson'
            
            if (missing(size)) {
              if (ncell(x) < 7000) size <- NULL
              else size <- 5000
            }
            
            if (missing(keep)) keep <- NULL
            #-----------
            
            if (is.null(size)) x <- as.data.frame(x,na.rm=TRUE)
            else x <- spatSample(x,size,na.rm=TRUE)
            
            x <- na.omit(x)
            
            vn <- colnames(x)
            
            if (!is.null(keep)) {
              if (is.numeric(keep)) {
                .w <- keep %in% c(1:length(vn))
                if (all(!.w)) stop(paste0('the values in keep are out of range; should be between 1 and ',ncol(x),' for your dataset...!'))
                if (any(!.w)) {
                  keep <- keep[.w]
                  warning('some of layer numbers in keep are out of range (should be between 1 and nlyr(x)), so they are ignored!')
                }
                keep <- vn[keep]
              } else if (is.character(keep)) {
                .w <- keep %in% colnames(x)
                if (all(!.w)) stop('None of the variable names in keep are available in x dataset!')
                if (any(!.w)) {
                  keep <- keep[.w]
                  warning('some of the variable names in keep are not available in x (check if it is typos), so they are ignored!')
                }
              }
            }
            #----------------
            
            .warn <- FALSE
            
            LOOP <- TRUE
            n <- new("VIF")
            n@variables <- colnames(x)
            exc <- c()
            excLog <- data.frame(step=integer(), excluded=character(), correlated_with=character(),
                                 correlation=numeric(), vif_excluded=numeric(), vif_retained=numeric(),
                                 stringsAsFactors=FALSE)
            step_i <- 0L
            if (is.null(keep)) {
              while (LOOP) {
                xcor <- abs(cor(x, method=method))
                mx <- .maxCor(xcor)
                if (xcor[mx[1],mx[2]] >= th) {
                  step_i <- step_i + 1L
                  w1 <- which(colnames(xcor) == mx[1])
                  w2 <- which(rownames(xcor) == mx[2])
                  v <- .vif2(x,c(w1,w2))
                  ex <- mx[which.max(v[mx])]
                  kp <- mx[mx != ex]
                  excLog <- rbind(excLog, data.frame(step=step_i, excluded=ex, correlated_with=kp,
                                   correlation=xcor[mx[1],mx[2]], vif_excluded=v[ex], vif_retained=v[kp],
                                   stringsAsFactors=FALSE))
                  exc <- c(exc,ex)
                  x <- x[,-which(colnames(x) == ex)]
                } else LOOP <- FALSE
              }
            } else {
              
              while (LOOP) {
                xcor <- abs(cor(x, method=method))
                mx <- .maxCor2(xcor,keep=keep)
                if (is.list(mx)) {
                  .warn <- TRUE
                  mx <- mx[[1]]
                }
                if (xcor[mx[1],mx[2]] >= th) {
                  step_i <- step_i + 1L
                  w1 <- which(colnames(xcor) == mx[1])
                  w2 <- which(rownames(xcor) == mx[2])
                  if (any(mx %in% keep)) {
                    ex <- mx[!mx %in% keep]
                    if (length(ex) == 1L) {
                      kp <- mx[mx != ex]
                      v_pair <- .vif2(x, c(w1, w2))
                      excLog <- rbind(excLog, data.frame(step=step_i, excluded=ex, correlated_with=kp,
                                       correlation=xcor[mx[1],mx[2]], vif_excluded=v_pair[ex], vif_retained=v_pair[kp],
                                       stringsAsFactors=FALSE))
                    }
                  } else {
                    v <- .vif2(x,c(w1,w2))
                    ex <- mx[which.max(v[mx])]
                    kp <- mx[mx != ex]
                    excLog <- rbind(excLog, data.frame(step=step_i, excluded=ex, correlated_with=kp,
                                     correlation=xcor[mx[1],mx[2]], vif_excluded=v[ex], vif_retained=v[kp],
                                     stringsAsFactors=FALSE))
                  }
                  
                  exc <- c(exc,ex)
                  x <- x[,-which(colnames(x) == ex)]
                } else LOOP <- FALSE
              }
            }
            
            if (.warn) {
              if (length(keep) > 2) warning('At least two variables specified in "keep" are strongly correlated (i.e., are subjected to the collinearity issue)!')
              else warning('The two variables specified in "keep" are strongly correlated (i.e., are subjected to the collinearity issue)!')
            }
            #---
            if (length(exc) > 0) n@excluded <- exc
            n@exclusionLog <- excLog
            n@chains <- .buildChains(excLog)
            v <- .vif(x)
            n@corMatrix <- cor(x, method=method)
            n@results <- data.frame(Variables=names(v),VIF=as.vector(v))
            n
          }
)
#----------
setMethod('vifcor', signature(x='data.frame'),
          function(x, th=0.9, keep=NULL, size, method='pearson') {
            if (ncol(x) == 1) stop("The data.frame should have at least two columns")
            if (missing(method) || !method %in% c('pearson','kendall','spearman')) method <- 'pearson'
            
            x <- na.omit(x)
            
            if (missing(size)) {
              if (nrow(x) < 7000) size <- NULL
              else size <- 5000
            }
            
            if (missing(keep)) keep <- NULL
            #-----------
            
            if (!is.null(size)) {
              if(nrow(x) > size) x <- x[sample(1:nrow(x),size),]
            }
            #----
            vn <- colnames(x)
            
            if (!is.null(keep)) {
              if (is.numeric(keep)) {
                .w <- keep %in% c(1:length(vn))
                if (all(!.w)) stop(paste0('the values in keep are out of range; should be between 1 and ',ncol(x),' for your dataset...!'))
                if (any(!.w)) {
                  keep <- keep[.w]
                  warning('some of layer numbers in keep are out of range (should be between 1 and ncol(x)), so they are ignored!')
                }
                keep <- vn[keep]
              } else if (is.character(keep)) {
                .w <- keep %in% colnames(x)
                if (all(!.w)) stop('None of the variable names in keep are available in x!')
                if (any(!.w)) {
                  keep <- keep[.w]
                  warning('some of the variable names in keep are not available in x (check if it is typos), so they are ignored!')
                }
              }
            }
            #----------------
            .warn <- FALSE
            
            LOOP <- TRUE
            n <- new("VIF")
            n@variables <- colnames(x)
            exc <- c()
            excLog <- data.frame(step=integer(), excluded=character(), correlated_with=character(),
                                 correlation=numeric(), vif_excluded=numeric(), vif_retained=numeric(),
                                 stringsAsFactors=FALSE)
            step_i <- 0L
            if (is.null(keep)) {
              while (LOOP) {
                xcor <- abs(cor(x, method=method))
                mx <- .maxCor(xcor)
                if (xcor[mx[1],mx[2]] >= th) {
                  step_i <- step_i + 1L
                  w1 <- which(colnames(xcor) == mx[1])
                  w2 <- which(rownames(xcor) == mx[2])
                  v <- .vif2(x,c(w1,w2))
                  ex <- mx[which.max(v[mx])]
                  kp <- mx[mx != ex]
                  excLog <- rbind(excLog, data.frame(step=step_i, excluded=ex, correlated_with=kp,
                                   correlation=xcor[mx[1],mx[2]], vif_excluded=v[ex], vif_retained=v[kp],
                                   stringsAsFactors=FALSE))
                  exc <- c(exc,ex)
                  x <- x[,-which(colnames(x) == ex)]
                } else LOOP <- FALSE
              }
            } else {
              
              while (LOOP) {
                xcor <- abs(cor(x, method=method))
                mx <- .maxCor2(xcor,keep=keep)
                if (is.list(mx)) {
                  .warn <- TRUE
                  mx <- mx[[1]]
                }
                if (xcor[mx[1],mx[2]] >= th) {
                  step_i <- step_i + 1L
                  w1 <- which(colnames(xcor) == mx[1])
                  w2 <- which(rownames(xcor) == mx[2])
                  if (any(mx %in% keep)) {
                    ex <- mx[!mx %in% keep]
                    if (length(ex) == 1L) {
                      kp <- mx[mx != ex]
                      v_pair <- .vif2(x, c(w1, w2))
                      excLog <- rbind(excLog, data.frame(step=step_i, excluded=ex, correlated_with=kp,
                                       correlation=xcor[mx[1],mx[2]], vif_excluded=v_pair[ex], vif_retained=v_pair[kp],
                                       stringsAsFactors=FALSE))
                    }
                  } else {
                    v <- .vif2(x,c(w1,w2))
                    ex <- mx[which.max(v[mx])]
                    kp <- mx[mx != ex]
                    excLog <- rbind(excLog, data.frame(step=step_i, excluded=ex, correlated_with=kp,
                                     correlation=xcor[mx[1],mx[2]], vif_excluded=v[ex], vif_retained=v[kp],
                                     stringsAsFactors=FALSE))
                  }
                  
                  exc <- c(exc,ex)
                  x <- x[,-which(colnames(x) == ex)]
                } else LOOP <- FALSE
              }
            }
            
            if (.warn) {
              if (length(keep) > 2) warning('At least two variables specified in "keep" are strongly correlated (i.e., are subjected to the collinearity issue)!')
              else warning('The two variables specified in "keep" are strongly correlated (i.e., are subjected to the collinearity issue)!')
            }
            #---
            if (length(exc) > 0) n@excluded <- exc
            n@exclusionLog <- excLog
            n@chains <- .buildChains(excLog)
            v <- .vif(x)
            n@corMatrix <- cor(x, method=method)
            n@results <- data.frame(Variables=names(v),VIF=as.vector(v))
            n
          }
)
#--------

setMethod('vifcor', signature(x='matrix'),
          function(x, th=0.9, keep=NULL, size, method='pearson') {
            if (ncol(x) == 1) stop("The matrix should have at least two columns")
            if (missing(method) || !method %in% c('pearson','kendall','spearman')) method <- 'pearson'
            
            x <- as.data.frame(x)
            x <- na.omit(x)
            
            if (missing(size)) {
              if (nrow(x) < 7000) size <- NULL
              else size <- 5000
            }
            
            if (missing(keep)) keep <- NULL
            #-----------
            
            if (!is.null(size)) {
              if(nrow(x) > size) x <- x[sample(1:nrow(x),size),]
            }
            #----
            vn <- colnames(x)
            
            if (!is.null(keep)) {
              if (is.numeric(keep)) {
                .w <- keep %in% c(1:length(vn))
                if (all(!.w)) stop(paste0('the values in keep are out of range; should be between 1 and ',ncol(x),' for your dataset...!'))
                if (any(!.w)) {
                  keep <- keep[.w]
                  warning('some of layer numbers in keep are out of range (should be between 1 and ncol(x)), so they are ignored!')
                }
                keep <- vn[keep]
              } else if (is.character(keep)) {
                .w <- keep %in% colnames(x)
                if (all(!.w)) stop('None of the variable names in keep are available in x!')
                if (any(!.w)) {
                  keep <- keep[.w]
                  warning('some of the variable names in keep are not available in x (check if it is typos), so they are ignored!')
                }
              }
            }
            #----------------
            .warn <- FALSE
            
            LOOP <- TRUE
            n <- new("VIF")
            n@variables <- colnames(x)
            exc <- c()
            excLog <- data.frame(step=integer(), excluded=character(), correlated_with=character(),
                                 correlation=numeric(), vif_excluded=numeric(), vif_retained=numeric(),
                                 stringsAsFactors=FALSE)
            step_i <- 0L
            if (is.null(keep)) {
              while (LOOP) {
                xcor <- abs(cor(x, method=method))
                mx <- .maxCor(xcor)
                if (xcor[mx[1],mx[2]] >= th) {
                  step_i <- step_i + 1L
                  w1 <- which(colnames(xcor) == mx[1])
                  w2 <- which(rownames(xcor) == mx[2])
                  v <- .vif2(x,c(w1,w2))
                  ex <- mx[which.max(v[mx])]
                  kp <- mx[mx != ex]
                  excLog <- rbind(excLog, data.frame(step=step_i, excluded=ex, correlated_with=kp,
                                   correlation=xcor[mx[1],mx[2]], vif_excluded=v[ex], vif_retained=v[kp],
                                   stringsAsFactors=FALSE))
                  exc <- c(exc,ex)
                  x <- x[,-which(colnames(x) == ex)]
                } else LOOP <- FALSE
              }
            } else {
              
              while (LOOP) {
                xcor <- abs(cor(x, method=method))
                mx <- .maxCor2(xcor,keep=keep)
                if (is.list(mx)) {
                  .warn <- TRUE
                  mx <- mx[[1]]
                }
                if (xcor[mx[1],mx[2]] >= th) {
                  step_i <- step_i + 1L
                  w1 <- which(colnames(xcor) == mx[1])
                  w2 <- which(rownames(xcor) == mx[2])
                  if (any(mx %in% keep)) {
                    ex <- mx[!mx %in% keep]
                    if (length(ex) == 1L) {
                      kp <- mx[mx != ex]
                      v_pair <- .vif2(x, c(w1, w2))
                      excLog <- rbind(excLog, data.frame(step=step_i, excluded=ex, correlated_with=kp,
                                       correlation=xcor[mx[1],mx[2]], vif_excluded=v_pair[ex], vif_retained=v_pair[kp],
                                       stringsAsFactors=FALSE))
                    }
                  } else {
                    v <- .vif2(x,c(w1,w2))
                    ex <- mx[which.max(v[mx])]
                    kp <- mx[mx != ex]
                    excLog <- rbind(excLog, data.frame(step=step_i, excluded=ex, correlated_with=kp,
                                     correlation=xcor[mx[1],mx[2]], vif_excluded=v[ex], vif_retained=v[kp],
                                     stringsAsFactors=FALSE))
                  }
                  
                  exc <- c(exc,ex)
                  x <- x[,-which(colnames(x) == ex)]
                } else LOOP <- FALSE
              }
            }
            
            if (.warn) {
              if (length(keep) > 2) warning('At least two variables specified in "keep" are strongly correlated (i.e., are subjected to the collinearity issue)!')
              else warning('The two variables specified in "keep" are strongly correlated (i.e., are subjected to the collinearity issue)!')
            }
            #---
            if (length(exc) > 0) n@excluded <- exc
            n@exclusionLog <- excLog
            n@chains <- .buildChains(excLog)
            v <- .vif(x)
            n@corMatrix <- cor(x, method=method)
            n@results <- data.frame(Variables=names(v),VIF=as.vector(v))
            n
          }
)
########################################

if (!isGeneric("vifstep")) {
  setGeneric("vifstep", function(x, th= 10,keep=NULL, size, method='pearson',...)
    standardGeneric("vifstep"))
}

setMethod('vifstep', signature(x='RasterStackBrick'),
          function(x, th=10, keep=NULL, size, method='pearson') {
            if (nlayers(x) == 1) stop("The Raster object should have at least two layers!")
            if (missing(method) || !method %in% c('pearson','kendall','spearman')) method <- 'pearson'
            
            if (missing(size)) {
              if (ncell(x) < 7000) size <- NULL
              else size <- 5000
            }
            
            if (missing(keep)) keep <- NULL
            #-----------
            
            if (is.null(size)) x <- as.data.frame(x,na.rm=TRUE)
            else x <- sampleRandom(x,size,na.rm=TRUE)
            
            x <- na.omit(x)
            
            vn <- colnames(x)
            
            if (!is.null(keep)) {
              if (is.numeric(keep)) {
                .w <- keep %in% c(1:length(vn))
                if (all(!.w)) stop(paste0('the values in keep are out of range; should be between 1 and ',ncol(x),' for your dataset...!'))
                if (any(!.w)) {
                  keep <- keep[.w]
                  warning('some of layer numbers in keep are out of range (should be between 1 and nlyr(x)), so they are ignored!')
                }
                keep <- vn[keep]
              } else if (is.character(keep)) {
                .w <- keep %in% colnames(x)
                if (all(!.w)) stop('None of the variable names in keep are available in x dataset!')
                if (any(!.w)) {
                  keep <- keep[.w]
                  warning('some of the variable names in keep are not available in x (check if it is typos), so they are ignored!')
                }
              }
            }
            #----------------
            LOOP <- TRUE
            n <- new("VIF")
            n@variables <- colnames(x)
            exc <- c()
            excLog <- data.frame(step=integer(), excluded=character(), correlated_with=character(),
                                 correlation=numeric(), vif_excluded=numeric(), vif_retained=numeric(),
                                 stringsAsFactors=FALSE)
            step_i <- 0L
            if (is.null(keep)) {
              while (LOOP) {
                v <- .vif(x)
                if (v[which.max(v)] >= th) {
                  step_i <- step_i + 1L
                  ex <- names(v[which.max(v)])
                  xcor_step <- abs(cor(x, method=method))
                  diag(xcor_step) <- NA
                  kp <- names(which.max(xcor_step[ex, ]))
                  excLog <- rbind(excLog, data.frame(step=step_i, excluded=ex, correlated_with=kp,
                                   correlation=xcor_step[ex, kp], vif_excluded=v[ex], vif_retained=NA_real_,
                                   stringsAsFactors=FALSE))
                  exc <- c(exc,ex)
                  x <- x[,-which(colnames(x) == ex)]
                } else LOOP=FALSE
              }
            } else {
              while (LOOP) {
                v <- .vif(x)
                v <- v[!names(v) %in% keep]
                if (v[which.max(v)] >= th) {
                  step_i <- step_i + 1L
                  ex <- names(v[which.max(v)])
                  xcor_step <- abs(cor(x, method=method))
                  diag(xcor_step) <- NA
                  kp <- names(which.max(xcor_step[ex, ]))
                  excLog <- rbind(excLog, data.frame(step=step_i, excluded=ex, correlated_with=kp,
                                   correlation=xcor_step[ex, kp], vif_excluded=v[ex], vif_retained=NA_real_,
                                   stringsAsFactors=FALSE))
                  exc <- c(exc,ex)
                  x <- x[,-which(colnames(x) == ex)]
                } else LOOP=FALSE
              }
            }
            #---
            if (length(exc) > 0) n@excluded <- exc
            n@exclusionLog <- excLog
            n@chains <- .buildChains(excLog)
            v <- .vif(x)
            n@corMatrix <- cor(x, method=method)
            n@results <- data.frame(Variables=names(v),VIF=as.vector(v))
            n
          }
)

setMethod('vifstep', signature(x='data.frame'),
          function(x, th=10, keep=NULL, size, method='pearson') {
            if (ncol(x) == 1) stop("The data.frame should have at least two variables!")
            if (missing(method) || !method %in% c('pearson','kendall','spearman')) method <- 'pearson'
            
            x <- na.omit(x)
            
            if (missing(size)) {
              if (nrow(x) < 6000) size <- NULL
              else size <- 5000
            }
            
            if (missing(keep)) keep <- NULL
            #-----------
            if (!is.null(size)) {
              if(nrow(x) > size) x <- x[sample(1:nrow(x),size),]
            }
            #----
            vn <- colnames(x)
            
            if (!is.null(keep)) {
              if (is.numeric(keep)) {
                .w <- keep %in% c(1:length(vn))
                if (all(!.w)) stop(paste0('the values in keep are out of range; should be between 1 and ',ncol(x),' for your dataset...!'))
                if (any(!.w)) {
                  keep <- keep[.w]
                  warning('some of layer numbers in keep are out of range (should be between 1 and ncol(x)), so they are ignored!')
                }
                keep <- vn[keep]
              } else if (is.character(keep)) {
                .w <- keep %in% colnames(x)
                if (all(!.w)) stop('None of the variable names in keep are available in x!')
                if (any(!.w)) {
                  keep <- keep[.w]
                  warning('some of the variable names in keep are not available in x (check if it is typos), so they are ignored!')
                }
              }
            }
            #----------------
            
            LOOP <- TRUE
            n <- new("VIF")
            n@variables <- colnames(x)
            exc <- c()
            excLog <- data.frame(step=integer(), excluded=character(), correlated_with=character(),
                                 correlation=numeric(), vif_excluded=numeric(), vif_retained=numeric(),
                                 stringsAsFactors=FALSE)
            step_i <- 0L
            if (is.null(keep)) {
              while (LOOP) {
                v <- .vif(x)
                if (v[which.max(v)] >= th) {
                  step_i <- step_i + 1L
                  ex <- names(v[which.max(v)])
                  xcor_step <- abs(cor(x, method=method))
                  diag(xcor_step) <- NA
                  kp <- names(which.max(xcor_step[ex, ]))
                  excLog <- rbind(excLog, data.frame(step=step_i, excluded=ex, correlated_with=kp,
                                   correlation=xcor_step[ex, kp], vif_excluded=v[ex], vif_retained=NA_real_,
                                   stringsAsFactors=FALSE))
                  exc <- c(exc,ex)
                  x <- x[,-which(colnames(x) == ex)]
                } else LOOP=FALSE
              }
            } else {
              while (LOOP) {
                v <- .vif(x)
                v <- v[!names(v) %in% keep]
                if (v[which.max(v)] >= th) {
                  step_i <- step_i + 1L
                  ex <- names(v[which.max(v)])
                  xcor_step <- abs(cor(x, method=method))
                  diag(xcor_step) <- NA
                  kp <- names(which.max(xcor_step[ex, ]))
                  excLog <- rbind(excLog, data.frame(step=step_i, excluded=ex, correlated_with=kp,
                                   correlation=xcor_step[ex, kp], vif_excluded=v[ex], vif_retained=NA_real_,
                                   stringsAsFactors=FALSE))
                  exc <- c(exc,ex)
                  x <- x[,-which(colnames(x) == ex)]
                } else LOOP=FALSE
              }
            }
            if (length(exc) > 0) n@excluded <- exc
            n@exclusionLog <- excLog
            n@chains <- .buildChains(excLog)
            v <- .vif(x)
            n@corMatrix <- cor(x, method=method)
            n@results <- data.frame(Variables=names(v),VIF=as.vector(v))
            n
          }
)

setMethod('vifstep', signature(x='matrix'),
          function(x, th=10, keep=NULL, size, method='pearson') {
            if (ncol(x) == 1) stop("The matrix should have at least two columns (variables)!")
            if (missing(method) || !method %in% c('pearson','kendall','spearman')) method <- 'pearson'
            
            x <- as.data.frame(x)
            x <- na.omit(x)
            
            if (missing(size)) {
              if (nrow(x) < 6000) size <- NULL
              else size <- 5000
            }
            
            if (missing(keep)) keep <- NULL
            #-----------
            
            if (!is.null(size)) {
              if(nrow(x) > size) x <- x[sample(1:nrow(x),size),]
            }
            #----
            vn <- colnames(x)
            
            if (!is.null(keep)) {
              if (is.numeric(keep)) {
                .w <- keep %in% c(1:length(vn))
                if (all(!.w)) stop(paste0('the values in keep are out of range; should be between 1 and ',ncol(x),' for your dataset...!'))
                if (any(!.w)) {
                  keep <- keep[.w]
                  warning('some of layer numbers in keep are out of range (should be between 1 and ncol(x)), so they are ignored!')
                }
                keep <- vn[keep]
              } else if (is.character(keep)) {
                .w <- keep %in% colnames(x)
                if (all(!.w)) stop('None of the variable names in keep are available in x!')
                if (any(!.w)) {
                  keep <- keep[.w]
                  warning('some of the variable names in keep are not available in x (check if it is typos), so they are ignored!')
                }
              }
            }
            #----------------
            LOOP <- TRUE
            
            n <- new("VIF")
            n@variables <- colnames(x)
            exc <- c()
            excLog <- data.frame(step=integer(), excluded=character(), correlated_with=character(),
                                 correlation=numeric(), vif_excluded=numeric(), vif_retained=numeric(),
                                 stringsAsFactors=FALSE)
            step_i <- 0L
            if (is.null(keep)) {
              while (LOOP) {
                v <- .vif(x)
                if (v[which.max(v)] >= th) {
                  step_i <- step_i + 1L
                  ex <- names(v[which.max(v)])
                  xcor_step <- abs(cor(x, method=method))
                  diag(xcor_step) <- NA
                  kp <- names(which.max(xcor_step[ex, ]))
                  excLog <- rbind(excLog, data.frame(step=step_i, excluded=ex, correlated_with=kp,
                                   correlation=xcor_step[ex, kp], vif_excluded=v[ex], vif_retained=NA_real_,
                                   stringsAsFactors=FALSE))
                  exc <- c(exc,ex)
                  x <- x[,-which(colnames(x) == ex)]
                } else LOOP=FALSE
              }
            } else {
              while (LOOP) {
                v <- .vif(x)
                v <- v[!names(v) %in% keep]
                if (v[which.max(v)] >= th) {
                  step_i <- step_i + 1L
                  ex <- names(v[which.max(v)])
                  xcor_step <- abs(cor(x, method=method))
                  diag(xcor_step) <- NA
                  kp <- names(which.max(xcor_step[ex, ]))
                  excLog <- rbind(excLog, data.frame(step=step_i, excluded=ex, correlated_with=kp,
                                   correlation=xcor_step[ex, kp], vif_excluded=v[ex], vif_retained=NA_real_,
                                   stringsAsFactors=FALSE))
                  exc <- c(exc,ex)
                  x <- x[,-which(colnames(x) == ex)]
                } else LOOP=FALSE
              }
            }
            if (length(exc) > 0) n@excluded <- exc
            n@exclusionLog <- excLog
            n@chains <- .buildChains(excLog)
            v <- .vif(x)
            n@corMatrix <- cor(x, method=method)
            n@results <- data.frame(Variables=names(v),VIF=as.vector(v))
            n
          }
)
#------------


setMethod('vifstep', signature(x='SpatRaster'),
          function(x, th=10, keep=NULL, size, method='pearson') {
            if (nlyr(x) == 1) stop("The Raster object should have at least two layers!")
            if (missing(method) || !method %in% c('pearson','kendall','spearman')) method <- 'pearson'
            if (missing(size)) {
              if (ncell(x) < 7000) size <- NULL
              else size <- 5000
            }
            
            if (missing(keep)) keep <- NULL
            #-----------
            
            if (is.null(size)) x <- as.data.frame(x,na.rm=TRUE)
            else x <- spatSample(x,size,na.rm=TRUE)
            
            x <- na.omit(x)
            
            vn <- colnames(x)
            
            if (!is.null(keep)) {
              if (is.numeric(keep)) {
                .w <- keep %in% c(1:length(vn))
                if (all(!.w)) stop(paste0('the values in keep are out of range; should be between 1 and ',ncol(x),' for your dataset...!'))
                if (any(!.w)) {
                  keep <- keep[.w]
                  warning('some of layer numbers in keep are out of range (should be between 1 and nlyr(x)), so they are ignored!')
                }
              keep <- vn[keep]
              } else if (is.character(keep)) {
                .w <- keep %in% colnames(x)
                if (all(!.w)) stop('None of the variable names in keep are available in x dataset!')
                if (any(!.w)) {
                  keep <- keep[.w]
                  warning('some of the variable names in keep are not available in x (check if it is typos), so they are ignored!')
                }
              }
            }
            #----------------
            LOOP <- TRUE
            n <- new("VIF")
            n@variables <- colnames(x)
            exc <- c()
            excLog <- data.frame(step=integer(), excluded=character(), correlated_with=character(),
                                 correlation=numeric(), vif_excluded=numeric(), vif_retained=numeric(),
                                 stringsAsFactors=FALSE)
            step_i <- 0L
            if (is.null(keep)) {
              while (LOOP) {
                v <- .vif(x)
                if (v[which.max(v)] >= th) {
                  step_i <- step_i + 1L
                  ex <- names(v[which.max(v)])
                  xcor_step <- abs(cor(x, method=method))
                  diag(xcor_step) <- NA
                  kp <- names(which.max(xcor_step[ex, ]))
                  excLog <- rbind(excLog, data.frame(step=step_i, excluded=ex, correlated_with=kp,
                                   correlation=xcor_step[ex, kp], vif_excluded=v[ex], vif_retained=NA_real_,
                                   stringsAsFactors=FALSE))
                  exc <- c(exc,ex)
                  x <- x[,-which(colnames(x) == ex)]
                } else LOOP=FALSE
              }
            } else {
              while (LOOP) {
                v <- .vif(x)
                v <- v[!names(v) %in% keep]
                if (v[which.max(v)] >= th) {
                  step_i <- step_i + 1L
                  ex <- names(v[which.max(v)])
                  xcor_step <- abs(cor(x, method=method))
                  diag(xcor_step) <- NA
                  kp <- names(which.max(xcor_step[ex, ]))
                  excLog <- rbind(excLog, data.frame(step=step_i, excluded=ex, correlated_with=kp,
                                   correlation=xcor_step[ex, kp], vif_excluded=v[ex], vif_retained=NA_real_,
                                   stringsAsFactors=FALSE))
                  exc <- c(exc,ex)
                  x <- x[,-which(colnames(x) == ex)]
                } else LOOP=FALSE
              }
            }
            
            if (length(exc) > 0) n@excluded <- exc
            n@exclusionLog <- excLog
            n@chains <- .buildChains(excLog)
            v <- .vif(x)
            n@corMatrix <- cor(x, method=method)
            n@results <- data.frame(Variables=names(v),VIF=as.vector(v))
            n
          }
)



