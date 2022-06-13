##R ROUTINES FOR THE R PACKAGE DLNM (C) ANTONIO GASPARRINI 2017##
#
cyclic <-
  function(x, df=NULL, knots=NULL, degree=3, intercept=FALSE, Boundary.knots=NULL) {
    #
    nx <- names(x)
    x <- as.vector(x)
    range <- range(x,na.rm=TRUE)
    bound <- Boundary.knots
    nax <- is.na(x)
    if(nas <- any(nax)) x <- x[!nax]
    if((degree <- as.integer(degree)) < 1) stop("'degree' must be integer >= 1")
    #
    # DEFINE BOUNDARY KNOTS
    if(is.null(bound)) {
      bound <- range
    } else {
      if(length(bound)!=2||diff(bound)<=0)
        stop(("boundary knots must be two increasing values"))
      if(min(x)<bound[1]||max(x)>bound[2])
        stop("values must be within boundary knots")
    }
    #
    # DEFINE INTERNAL KNOTS AND DF
    if(is.null(knots)) {
      if(is.null(df)) df <- max(1,degree-1+intercept)
      nik <- df-intercept
      if(nik<degree-1) stop("basis dimension too small for b-spline degree")
      dx <- diff(range)/(nik+1)
      knots <- seq(range[1]+dx, range[2]-dx, length=nik)
    } else {
      if(min(knots)<=bound[1]||max(knots)>=bound[2])
        stop("internal knots must be within boundary knots")
      df <- length(knots)+2+intercept
    }
    #
    #  TRANSFORMATION AND INTERCEPT
    basis <- cSplineDes(x,c(bound[1],knots,bound[2]),degree+1)
    if(!intercept) basis <- basis[,-1L,drop=FALSE]
    #
    ################################################################################
    #
    # RE-INSERT MISSING
    if(nas) {
      nmat <- matrix(NA,length(nax),ncol(basis))
      nmat[!nax,] <- basis
      basis <- nmat
    }
    #
    # NAMES AND ATTRIBUTES
    dimnames(basis) <- list(nx,seq(ncol(basis)))
    attributes(basis) <- c(attributes(basis),list(df=df,knots=knots,degree=degree,
                                                  intercept=intercept,Boundary.knots=bound))
    #
    class(basis) <- c("cyclic","matrix")
    #
    return(basis)
  }
