###
### ADAPTED FROM FUNCTION DEVELOPED BY Antonio Gasparrini 2014
# https://github.com/gasparrini/2014_gasparrini_BMCmrm_Rcodedata/blob/master/attrdl.R
################################################################################
# FUNCTION FOR COMPUTING ATTRIBUTABLE MEASURES FROM DLNM
#   REQUIRES dlnm VERSION 2.0.9 AND ON
################################################################################
#
# DISCLAIMER:
#   THE CODE COMPOSING THIS FUNCTION HAS NOT BEEN SYSTEMATICALLY TESTED. THE
#   PRESENCE OF BUGS CANNOT BE RULED OUT. ALSO, ALTHOUGH WRITTEN GENERICALLY
#   FOR WORKING IN DIFFERENT SCENARIOS AND DATA, THE FUNCTION HAS NOT BEEN
#   TESTED IN CONTEXTS DIFFERENT THAN THE EXAMPLE INCLUDED IN THE PAPER.
#   IT IS RESPONSIBILITY OF THE USER TO CHECK THE RELIABILITY OF THE RESULTS IN
#   DIFFERENT APPLICATIONS.

################################################################################
# SEE THE PDF WITH A DETAILED DOCUMENTATION AT www.ag-myresearch.com
#
#   - x: AN EXPOSURE VECTOR OR (ONLY FOR dir="back") A MATRIX OF LAGGED EXPOSURES
#   - basis: THE CROSS-BASIS COMPUTED FROM x
#   - cases: THE CASES VECTOR OR (ONLY FOR dir="forw") THE MATRIX OF FUTURE CASES
#   - model: THE FITTED MODEL
#   - coef AND vcov: COEF AND VCOV FOR basis IF model IS NOT PROVIDED
#   - type: EITHER "an" OR "af" FOR ATTRIBUTABLE NUMBER OR FRACTION
#   - tot: IF TRUE, THE TOTAL ATTRIBUTABLE RISK IS COMPUTED
#   - sim: IF SIMULATION SAMPLES SHOULD BE RETURNED. ONLY FOR tot=TRUE
#   - nsim: NUMBER OF SIMULATION SAMPLES
################################################################################
attrs <- function(x,basis,data,model=NULL,coef=NULL,vcov=NULL,type,
 tot=TRUE,sim=FALSE,nsim=5000) {
################################################################################
  type <- match.arg(type,c("an","af"))

  # cases:
  cases<-matrix(by(data$death,data$seasonal,FUN = mean))
################################################################################
#
  # EXTRACT COEF AND VCOV IF MODEL IS PROVIDED

################################################################################
#
  # CHECK AND DEFINE BASIS
  if(!any(class(basis)%in%c("cyclic")))
    stop("the first argument must be an object of class 'cyclic'")

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
#
  ## for centering##
  maxcen<- apply(logrr, 1, which.min )
  bvar_cen<- bvar[maxcen,]

  # centered logrr
  l <- length(Trange)
  for(j in 1:l){ logrr[1,j] <- as.numeric( t(bvar[j,]-bvar_cen) %*% as.numeric( coef) ) }
  # standard error of logrr
  for(j in 1:l){ logrr_se[1,j] <- sqrt( as.numeric( t(bvar[j,]-bvar_cen) %*% vcov%*% (bvar[j,]-bvar_cen) ) ) }

################################################################################
#
  # COMPUTE AF AND AN
  af <- 1-exp(-logrr)
  an <- af[1,]*cases[,1]
#
  # TOTAL
  if(tot) {
    isna <- is.na(an)
    an <- sum(an,na.rm=T)
    af <- an/sum(cases[!isna])
  }
#
################################################################################
#
  # EMPIRICAL CONFIDENCE INTERVALS
  if(!tot && sim) {
    sim <- FALSE
    warning("simulation samples only returned for tot=T")
  }
  if(sim) {
    # SAMPLE COEF
    k <- length(coef)
    eigen <- eigen(vcov)
    X <- matrix(rnorm(length(coef)*nsim),nsim)
    coefsim <- coef + eigen$vectors %*% diag(sqrt(eigen$values),k) %*% t(X)
    # RUN THE LOOP
    # pre_afsim <- (1 - exp(- Xpredall %*% coefsim)) * cases # a matrix
    # afsim <- colSums(pre_afsim,na.rm=TRUE) / sum(cases[!isna],na.rm=TRUE)
    logrri<- matrix(NA, 1, 366 )
    afsim <- apply(coefsim,2, function(coefi) {
      for(j in 1:l){ logrri[1,j] <-  t(bvar[j,]) %*% as.numeric( coefi)  }
      ## for centering##
      maxceni<- apply(logrri, 1, which.min )
      bvar_ceni<- bvar[maxceni,]
      # centered logrr
      l <- length(Trange)
      for(j in 1:l){ logrri[1,j] <- as.numeric( t(bvar[j,]-bvar_ceni) %*% as.numeric( coefi) ) }
      # COMPUTE AF AND AN
      afi <- 1-exp(-logrri)
      ani <- afi[1,]*cases[,1]
      #
      # TOTAL
      if(tot) {
        isna <- is.na(ani)
        ani <- sum(ani,na.rm=T)
        afi <- ani/sum(cases[!isna])
      }
      sum(ani[!is.na(ani)])/sum(cases[!is.na(ani)])
    })
    #ansim <- afsim*den
    ansim <- afsim
  }
#
################################################################################
#    minsim <- apply(coefsim,2,function(coefi) {

  res <- if(sim) {
    if(type=="an") ansim else afsim
  } else {
    if(type=="an") an else af
  }
#
  return(res)
}

#
