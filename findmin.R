##################################################################################################################
# FUNCTION TO ESTIMATE TROUGH OF SEASONALITY
#ADAPTED FROM FUNCTION TO ESTIMATE MINIMUM OF AN EXPOSURE-RESPONSE FUNCTION FROM A FITTED MODEL
#Tobias, Armstrong, Gasparrini. Epidemiol, 2017.
##################################################################################################################

##################################################################################################################
#
findmin <- function(basis,model=NULL,coef=NULL,vcov=NULL,at=NULL,from=NULL,to=NULL,by=NULL,sim=FALSE,nsim=5000) {

  ################################################################################################################
  # ARGUMENTS:
  # - basis: A SPLINE OR OTHER BASIS FOR AN EXPOSURE x CREATED BY DLNM FUNCTION CROSSBASIS OR ONEBASIS
  # - model: THE modelTED MODEL
  # - coef AND vcov: COEF AND VCOV FOR basis IF model IS NOT PROVIDED
  # - at: A NUMERIC VECTOR OF x VALUES OVER WHICH THE MINIMUM IS SOUGHT
  # OR
  # - from, to: RANGE OF x VALUES OVER WHICH THE MINIMUM IS SOUGHT.
  # - by: INCREMENT OF THE SEQUENCES x VALUES OVER WHICH THE MINIMUM IS SOUGHT
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
  # FIND THE MINIMUM
  min <-  apply(logrr, 1, which.min )

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
    # COMPUTE MINIMUM
    minsim <- apply(coefsim,2,function(coefi) {
      for(j in 1:l){ logrr[1,j] <-  t(bvar[j,]) %*% as.numeric( coefi)  }
      return(apply(logrr, 1, which.min ))
    })

    # RECENTER TROUGH
    recenter_doy<-function(doy,oldcenter,newcenter=182,daysinyear=366) {
      return( 1+(doy-oldcenter+newcenter-1)%%daysinyear)
    }

    trough<-recenter_doy(minsim,oldcenter = min,newcenter = 182)
    trough.low<- recenter_doy(quantile(trough,c(2.5)/100),oldcenter=182,newcenter=min)
    trough.up<- recenter_doy(quantile(trough,c(97.5)/100),oldcenter=182,newcenter=min)

    if(sim)
      (warning("WARNING: This function assumes that the seasonality curve has a single trough with (for the CI) a symmetric sampling distribution!"))
    else (NA)

    est_ci<-round(c(trough=min,trough=trough.low,trough=trough.up),digits = 0)
    names(est_ci)[2:3]<-c("trough.low","trough.up")
  }
  #
  ################################################################################################################
  #
  res <- if(sim) est_ci else min
  #
  return(res)
}

# END OF FUNCTION ################################################################################################



