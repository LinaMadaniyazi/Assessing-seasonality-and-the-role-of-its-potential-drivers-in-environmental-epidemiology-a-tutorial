##################################################################################################################
# FUNCTION TO ESTIMATE PEAK OF SEASONALITY
#ADAPTED FROM FUNCTION TO ESTIMATE MINIMUM OF AN EXPOSURE-RESPONSE FUNCTION FROM A FITTED MODEL
#Tobias, Armstrong, Gasparrini. Epidemiol, 2017.
##################################################################################################################

##################################################################################################################
#
findmax <- function(basis,model=NULL,coef=NULL,vcov=NULL,at=NULL,from=NULL,to=NULL,by=NULL,sim=FALSE,nsim=5000) {

  ################################################################################################################
  # ARGUMENTS:
  # - basis: A SPLINE OR OTHER BASIS FOR AN EXPOSURE x CREATED BY DLNM FUNCTION CROSSBASIS OR ONEBASIS
  # - model: THE modelTED MODEL
  # - coef AND vcov: COEF AND VCOV FOR basis IF model IS NOT PROVIDED
  # - at: A NUMERIC VECTOR OF x VALUES OVER WHICH THE MAXIMUM IS SOUGHT
  # OR
  # - from, to: RANGE OF x VALUES OVER WHICH THE MAXIMUM IS SOUGHT.
  # - by: INCREMENT OF THE SEQUENCES x VALUES OVER WHICH THE MAXIMUM IS SOUGHT
  # - sim: IF BOOTSTRAP SIMULATION SAMPLES SHOULD BE RETURNED
  # - nsim: NUMBER OF SIMULATION SAMPLES
  ################################################################################################################

  ################################################################################################################
  # CHECK AND DEFINE BASIS
  if(!any(class(basis)%in%c("crossbasis","cyclic")))
    stop("the first argument must be an object of class 'crossbasis' or 'cyclic'")

  # SET COEF, VCOV CLASS AND LINK
  if(!is.null(model)) {
    model.class <- class(model)
    coef <- dlnm:::getcoef(model,model.class)
    ind <- grep("spline.season",names(coef))
    coef <- coef[ind]
    vcov <- dlnm:::getvcov(model,model.class)[ind,ind,drop=FALSE]
    model.link <- dlnm:::getlink(model,model.class)
  } else model.class <- NA
  #
  # CHECK
  if(length(coef)!=ncol(basis) || length(coef)!=dim(vcov)[1] || any(is.na(coef)) || any(is.na(vcov)))
    stop("model or coef/vcov not consistent with basis")

  ##calculate logrr, se##
  logrr<- matrix(NA, 1, 366 )
  logrr_se<- matrix(NA, 1, 366 )

  Trange <- 1:366
  bvar<- basis

  l <- length(Trange)
  for(j in 1:l){ logrr[1,j] <-  t(bvar[j,]) %*% as.numeric( coef)  }
  # standard error of logrr
  for(j in 1:l){ logrr_se[1,j] <- sqrt( as.numeric( t(bvar[j,]) %*% vcov %*% (bvar[j,]) ) ) }

  ################################################################################################################
  # FIND THE MAXIMUM
  max <-  apply(logrr, 1, which.max )


  ################################################################################################################
  # SIMULATIONS
  #
  if(sim) {
    # SIMULATE COEFFICIENTS
    k <- length(coef)
    eigen <- eigen(vcov)
    set.seed(123)
    X <- matrix(rnorm(length(coef)*nsim),nsim)
    coefsim <- coef + eigen$vectors %*% diag(sqrt(eigen$values),k) %*% t(X)
    # COMPUTE MAXIMUM
    maxsim <- apply(coefsim,2,function(coefi) {
      for(j in 1:l){ logrr[1,j] <-  t(bvar[j,]) %*% as.numeric( coefi)  }
      return(apply(logrr, 1, which.max ))
    })


    recenter_doy<-function(doy,oldcenter,newcenter=182,daysinyear=366) {
      return( 1+(doy-oldcenter+newcenter-1)%%daysinyear)
    }

    # RECENTER PEAK
    recenter_doy<-function(doy,oldcenter,newcenter=182,daysinyear=366) {
      return( 1+(doy-oldcenter+newcenter-1)%%daysinyear)
    }

    peak<-recenter_doy(maxsim,oldcenter = max,newcenter = 182)
    peak.low<- recenter_doy(quantile(peak,c(2.5)/100),oldcenter=182,newcenter=max)
    peak.up<- recenter_doy(quantile(peak,c(97.5)/100),oldcenter=182,newcenter=max)

    if(sim)
      (warning("WARNING: This function assumes that the seasonality curve has a single peak with (for the CI) a symmetric sampling distribution!"))
    else (NA)

    est_ci<-round(c(peak=max,peak=peak.low,peak=peak.up),digits = 0)
    names(est_ci)[2:3]<-c("peak.low","peak.up")
  }
  #
  ################################################################################################################
  #

  #
#  recenter_sim<-recenter_doy(sim,oldcenter=max,newcenter=182)
  res <- if(sim)  est_ci else max
#  peak.low<- recenter_doy(quantile(recenter_sim,c(2.5)/100),oldcenter=182,newcenter=max)
 # peak.up<- recenter_doy(quantile(recenter_sim,c(97.5)/100),oldcenter=182,newcenter=max)
  #
#
  return(res)

}

# END OF FUNCTION ################################################################################################



