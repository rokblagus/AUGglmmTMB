


#' Matrix power operator
#'
#' Computes the matrix raised to a given power using the eigen decomposition.
#' This is an internal helper function and is not exported.
#'
#' @param x A square numeric matrix.
#' @param n Numeric exponent.
#'
#' @return A numeric matrix representing `x^n`.
#' @keywords internal
#' @examples
#' \dontrun{
#' A <- matrix(c(2,1,1,2),2,2)
#' A %^% 2
#' }

"%^%" <- function(x, n)   with(eigen(x), vectors %*% (values^n * t(vectors)))


#' Compute fixed-effects penalty strength
#'
#' Computes the penalty strength for the fixed effects as suggested by Košuta et al.
#'
#' @param data A list containing the model data. Must include the element
#' \code{X} for the fixed effects design matrix. Optionally it also includes \code{M} for the number of trials and  \code{W} for the weights; if missing they are set to 1 internally for all observations.
#'
#' @return A numeric value equal to \eqn{\sqrt{p/n}}, where \eqn{p} is the number
#' of fixed-effect parameters and \eqn{n} is the effective sample size.
#'
#' @details
#' The penalty strength is defined as \eqn{\sqrt{p/n}}, following the recommendation
#' of Košuta et al., where \eqn{p} is the number of columns of \code{X} and \eqn{n}
#' is the effective sample size.
#'
#' @examples
#' data(birds)
#' Y<-birds$parasites
#' X <- model.matrix(parasites~migration+food, birds)
#' Z1<-model.matrix(parasites~migration, birds)
#' Z2<-model.matrix(parasites~1, birds)
#' grouping1<-as.numeric(as.factor(birds$Phylogenetic.Tomi))
#' grouping2<-as.numeric(as.factor(birds$species))
#' xdf<-list(Y=Y,X=X,Z=list(Z1,Z2),
#'   grouping=list(grouping1,grouping2),
#'   M=rep(1,nrow(X)),W=rep(1,nrow(X)))
#'
#' mv_multiplier(xdf)
#'
#' @export


mv_multiplier <- function(data) {
  if (is.null(data$X)) stop("The suplied data are not in correct format, X is missing.")
 # if (is.null(data$M)|is.null(data$W)) stop("The suplied data are not in correct format.")
  if (is.null(data$M)) data$M<-rep(1,nrow(data$X))
  if (is.null(data$W)) data$W<-rep(1,nrow(data$X))

  n <- sum(data$M*data$W)
  p <- ncol(data$X)
  return( sqrt(p / n))
}

#' Construct model formula from data list
#'
#' Builds a \code{formula} object consistent with the assumed structure of the
#' input \code{data}. The function supports both single and multiple random-effects
#' terms.
#'
#' @param data A list containing model components, including \code{Y}, \code{M},
#' \code{X}, \code{Z}, and \code{grouping}. If multiple random effects are present,
#' \code{Z} and \code{grouping} must be lists.
#'
#' @return A \code{formula} object suitable for fitting a GLMM.
#'
#' @details
#' The response is constructed as \code{cbind(Y, M - Y)} and the fixed-effects
#' component is specified as \code{-1 + X}. Random-effects terms are added as
#' \code{(-1 + Z_i | grouping_i)} for each random-effects component.
#'
#' @keywords internal

make_formula<-function(data){
  t1<-"cbind(Y,M-Y)~-1+X"
  if (!inherits(data$Z,"list")) t2<-"(-1+Z1|grouping1)" else {
    t2<-rep(NA,length(data$Z))
    for (i in 1:length(data$Z)){
      t2[i]<-paste0("(-1+Z",i,"|grouping",i,")")
    }
    t2<-paste(t2,collapse="+")

  }
  formula(paste(c(t1,t2),collapse = "+"))

}

#' Restructure input data for model fitting
#'
#' Converts the input \code{data} list into a standardized format required for
#' model fitting. In particular, random-effects components are renamed and
#' expanded to match the expected naming convention.
#'
#' @param data A list containing model components, including \code{Y}, \code{M},
#' \code{W}, \code{X}, \code{Z}, and \code{grouping}. The elements \code{Z} and
#' \code{grouping} can be either single objects or lists for multiple random effects.
#'
#' @return A list with elements \code{Y}, \code{M}, \code{W}, \code{X}, and
#' appropriately named random-effects components \code{Z1, Z2, ...} and
#' \code{grouping1, grouping2, ...}.
#'
#' @details
#' If a single random-effects term is provided, it is stored as \code{Z1} and
#' \code{grouping1}. If multiple random effects are present, each component is
#' assigned a unique name following the pattern \code{Z_i} and \code{grouping_i}.
#'
#' @keywords internal


restructure_data<-function(data){
  data_rest<-list(Y=data$Y,M=data$M,W=data$W,X=data$X)
  if (!inherits(data$Z,"list")) {
    data_rest$Z1<-data$Z
    data_rest$grouping1<-data$grouping
  } else {
    for (i in 1:length(data$Z)){
      nms_Z<-paste0("Z",i)
      nms_gr<-paste0("grouping",i)
      data_rest[[nms_Z]]<-data$Z[[i]]
      data_rest[[nms_gr]]<-data$grouping[[i]]
    }
  }
  data_rest
}


#' Construct augmented data for fixed-effects penalization
#'
#' Creates an augmented dataset by adding pseudo-observations corresponding to
#' the fixed-effects penalty.
#'
#' @param data A list containing the model data with elements \code{Y}, \code{M},
#' \code{X}, \code{W}, and random-effects components \code{Z1, Z2, ...}
#' with corresponding grouping variables \code{grouping1, grouping2, ...}.
#'
#' @param cfe Numeric. Strength of the fixed-effects penalty.
#'
#' @param fix_beta Numeric vector of fixed-effects coefficients used to compute
#' the linear predictor.
#'
#' @param link_fun Character. Link function. Supported options are
#' \code{"logit"}, \code{"probit"}, \code{"cloglog"}, \code{"cauchit"},
#' and \code{"loglog"}.
#'
#' @return A list with the same structure as \code{data}, containing the augmented
#' response \code{Y}, total counts \code{M}, design matrix \code{X}, weights
#' \code{W}, and corresponding random-effects components. The augmented dataset
#' includes two additional pseudo-observations per original observation.
#'
#' @details
#' The function constructs pseudo-data corresponding to the fixed-effects penalty
#' using a data-augmentation approach.
#'
#' @examples
#' data(birds)
#' Y<-birds$parasites
#' X <- model.matrix(parasites~migration+food, birds)
#' Z1<-model.matrix(parasites~migration, birds)
#' Z2<-model.matrix(parasites~1, birds)
#' grouping1<-as.numeric(as.factor(birds$Phylogenetic.Tomi))
#' grouping2<-as.numeric(as.factor(birds$species))
#' xdf<-list(Y=Y,X=X,Z=list(Z1,Z2),
#'   grouping=list(grouping1,grouping2),
#'   M=rep(1,nrow(X)),W=rep(1,nrow(X)))
#'
#' xdft<-AUGglmmTMB:::restructure_data(xdf)
#' pd_fix<-make_pd_fix(xdft,mv_multiplier(xdft),rep(0,ncol(xdft$X)),"logit")
#'
#' lapply(xdft,function(x) if(is.null(dim(x))) length(x) else dim(x) )
#' lapply(pd_fix,function(x) if(is.null(dim(x))) length(x) else dim(x) )
#'
#' @export


make_pd_fix<-function(data,cfe,fix_beta,link_fun){

  get_mu <- switch(link_fun,
                   logit = function(x) {
                     pn<-plogis(x)
                     pn[pn<1e-12]<-1e-12
                     pn[pn>(1-1e-12)]<-1-1e-12
                     pn
                     },
                   probit = function(x) {
                     pn<-pnorm(x)
                     pn[pn<1e-12]<-1e-12
                     pn[pn>(1-1e-12)]<-1-1e-12
                     pn
                     },
                   cloglog = function(x) {
                     pn<-1 - exp(-exp(x))
                     pn[pn<1e-12]<-1e-12
                     pn[pn>(1-1e-12)]<-1-1e-12
                     pn
                     },
                   cauchit = function(x) {
                     pn<-atan(x)/pi+0.5
                     pn[pn<1e-12]<-1e-12
                     pn[pn>(1-1e-12)]<-1-1e-12
                     pn
                     },
                   loglog = function(x) {
                     pn<-exp(-exp(-x))
                     pn[pn<1e-12]<-1e-12
                     pn[pn>(1-1e-12)]<-1-1e-12
                     pn
                     },
                   stop("Unsupported link function")
  )

  get_omega<-switch(link_fun,
                    logit = function(x) {
                      mu <- plogis(x)
                    mu*(1-mu)},
                    probit = function(x) {
                      pn<-pnorm(x)
                      pn[pn<1e-12]<-1e-12
                      pn[pn>(1-1e-12)]<-1-1e-12

                      dnorm(x)**2/(pn*(1-pn))
                    },
                    cloglog = function(x) {
                      exp(2*x)/(exp(exp(x))-1)
                    },
                    cauchit = function(x) {
                     1/( (1+x**2)**2*(pi**2/4-atan(x)**2))
                    },
                    loglog = function(x) {
                      exp(-2*x)/(exp(exp(-x))-1)
                    },
                        stop("Unsupported link function")
  )

  get_d_primeomega <- switch(link_fun,
                             logit = function(x) {
                               mu <- plogis(x)
                               1 - 2 * mu
                             },

                             probit = function(x) {
                               phi <- dnorm(x)
                               phi[phi<1e-12]<-1e-12
                               Phi <- pnorm(x)
                               -(x * Phi * (1 - Phi)) / phi
                             },

                             cloglog = function(x) {
                              # ((1 - exp(x)) * exp(x - exp(x)) * (exp(exp(x)) - 1)) / exp(2 * x)
                             (1-exp(-exp(x)))*(1-exp(x))/exp(x)
                               },

                             loglog = function(x) {
                               (1-exp(-exp(-x)))*(exp(-x)-1)/exp(-x)
                             },

                             cauchit = function(x) {
                               -(2 * x / pi) * (pi^2 / 4 - atan(x)^2)
                             },

                             stop("Unsupported link function")
  )


  X<-rbind(data$X,data$X,data$X)








  lp<-data$X%*%matrix(fix_beta,ncol=1)


  omega<-c(get_omega(lp))

  W<-diag(data$W*data$M*omega)
  Xw <- (W %^% 0.5) %*% data$X
  qrXw <- qr(Xw)
  H <- tcrossprod(qr.Q(qrXw))
  hi<-diag(H)

  W<-c(data$W,hi*cfe/data$M,hi*cfe/data$M)
mud<-get_mu(lp)
  #dp<-data$M*get_d_primeomega(lp)+2*get_mu(lp)-1
  dp<-data$M*get_d_primeomega(lp)+2*mud-1
  #y1<-data$Y+data$M*(get_d_primeomega(lp)+2*get_mu(lp)-1)
  #y2<-data$M-data$Y+data$M*(get_d_primeomega(lp)+2*get_mu(lp)-1)
  y1<-data$Y+dp
  y2<-data$M-data$Y+dp

  #ksi<-function(x){ #x is lp!
  #  p1<--data$M*(get_d_primeomega(x)+2*get_mu(x)-1)/get_mu(x)
  #  p2<-data$M*(get_d_primeomega(x)+2*get_mu(x)-1)/(1-get_mu(x))
  #  pp<-cbind(p1,p2)

   # apply(pp,1,max)

  #}
  #ksi<-function(x){ #x is lp!
   # p1<--dp/get_mu(x)
  #  p2<-dp/(1-get_mu(x))
    p1<--dp/mud
    p2<-dp/(1-mud)
    pp<-cbind(p1,p2)

   ksid<- apply(pp,1,max)

  #}

  #y1<-y1+ksi(lp)*get_mu(lp)
  y1<-y1+ksid*mud
  y1[y1 < 0] <- 0
  #y2<-y2+ksi(lp)*get_mu(lp)
  y2<-y2+ksid*mud
  y2[y2 < 0] <- 0
  #m<-data$M+ksi(lp)
  m<-data$M+ksid

  Y<-c(data$Y,y1,y2)
  M<-c(data$M,m,m)


  pd1<-list(Y=Y,M=M,X=X,W=W)



  for (i in 1:sum(grepl("Z",names(data)))){
    nmZ<-paste0("Z",i)
    nmgr<-paste0("grouping",i)
    Z0<-matrix(0,nrow=nrow(data[[nmZ]]),ncol=ncol(data[[nmZ]]))
    Z<-rbind(data[[nmZ]],Z0,Z0)

    grouping<-c(data[[nmgr]],data[[nmgr]],data[[nmgr]])
    pd1[[nmZ]]<-Z
    pd1[[nmgr]]<-grouping
  }

  pd1


}




#' Create pseudo-data for random-effects penalty
#'
#' Constructs pseudo-observations corresponding to the random-effects penalty
#' using an eigen-decomposition approach.
#'
#' @param psi Numeric matrix. The random-effects covariance (or precision) matrix
#'   to be penalized.
#' @param nu Numeric. Degrees-of-freedom parameter controlling the strength of the
#'   penalty.
#' @param const Numeric. Large constant used as the total count for pseudo-observations
#'   (defaults to 1e8).
#' @param param Character. Determines whether \code{psi} is interpreted as a
#'   \code{"precision"} or \code{"variance"} matrix. Must be one of \code{"precision"} or
#'   \code{"variance"}.
#' @param inv_link_fun Function. Inverse link function applied to the pseudo-data. Defaults
#'   to the logit inverse: \code{function(x) 1/(1+exp(-x))}.
#'
#' @return A list containing a single element \code{data}, which itself is a list with
#'   the following components:
#'   \describe{
#'     \item{\code{Y}}{Vector of pseudo-responses for the random-effects penalty.}
#'     \item{\code{grouping}}{Integer vector indicating the grouping of each pseudo-observation.}
#'     \item{\code{nn}}{Vector of total counts for the pseudo-observations.}
#'     \item{\code{Z}}{Random-effects design matrix corresponding to the pseudo-data.}
#'   }
#'
#' @details
#' The function performs an eigen-decomposition of \code{psi} (scaled according
#' to \code{nu} and \code{param}) to generate pseudo-data that encode the
#' random-effects penalty. The pseudo-data
#' can then be combined with the original data for penalized GLMM fitting.
#'
#' @examples
#' \dontrun{
#' psi <- matrix(c(1,0.5,0.5,1), ncol=2)
#' pseudo <- make_pseudo_data_rand(psi, nu=5)
#' str(pseudo)
#' }
#'
#' @export


make_pseudo_data_rand<-function(psi,nu,const=1e8,param="precision",inv_link_fun=function(x) 1/(1+exp(-x))){
  if (is.null(match.arg(param,c("precision","variance")))) stop("param needs to be one of: precision,variance")

  q<-ncol(psi)
  if (param=="precision") cc<-(nu-q-1)/q
  if (param=="variance") cc<-(nu+q+1)/q

  cc<-max(c(floor(cc),1))
  true<-psi/cc
  ee<-eigen(true,TRUE)
  ui<-list()
  for (j in 1:q){
    ui[[j]]<-sqrt(ee$values[j])*ee$vectors[,j]
  }

  pi<-list()

  for (j in 1:length(ui)){
    I<-diag(rep(1,length(ui[[j]])))
    pi[[j]]<-inv_link_fun(I%*%matrix(ui[[j]],ncol=1))
  }
  Y<-unlist(pi)

  id<-rep(1:q,each=q)
  n<-rep(const,length(id))

  Zi<-matrix(0,ncol=q,nrow=q)

  for (j in 1:q){
    Zi[j,j]<-1
  }
  for (j in 1:q){
    if (j==1) Z=Zi else Z<-rbind(Z,Zi)
  }

  fact<-cc
  if (fact>1){
    Y<-rep(Y,fact)
    n<-rep(n,fact)
    id<-rep(1:(q*fact),each=q)
    for (j in 1:(q*fact)){
      if (j==1) Z=Zi else Z<-rbind(Z,Zi)
    }
  }

  data0<-list(Y=Y,grouping=id,nn=n,Z=Z)


  list(data=data0)
}


#' Create augmented dataset for random-effects penalty
#'
#' Constructs an augmented dataset by adding pseudo-data corresponding to a
#' random-effects penalty for a selected random effect. This allows penalized
#' estimation of random effects.
#'
#' @param data List. Original dataset with elements \code{Y}, \code{M}, \code{W},
#'   \code{X}, and random-effects components \code{Z} and \code{grouping}.
#'   \code{Z} and \code{grouping} may be lists if multiple random effects are present.
#' @param nu Numeric vector or list. Degrees-of-freedom parameters for the random-effects
#'   penalties, one per random effect.
#' @param psi List of numeric matrices. Random-effects covariance (or variance) matrices,
#'   one per random effect.
#' @param link_fun Character. The inverse link function to use. Supported values:
#'   \code{"logit"}, \code{"probit"}, \code{"cloglog"}, \code{"cauchit"}, \code{"loglog"}.
#' @param which_re Integer. Index of the random effect for which the pseudo-data are generated
#'   (defaults to 1).
#'
#' @return A list representing the augmented dataset with the following elements:
#'   \describe{
#'     \item{\code{Y}}{Response vector including original and pseudo-observations.}
#'     \item{\code{M}}{Total counts vector including original and pseudo-observations.}
#'     \item{\code{X}}{Design matrix for fixed effects (original and zero-padded rows).}
#'     \item{\code{W}}{Observation weights including pseudo-data weights.}
#'     \item{\code{Z1, Z2, ...}}{Random-effects design matrices including pseudo-data.}
#'     \item{\code{grouping1, grouping2, ...}}{Grouping vectors including pseudo-data.}
#'   }
#'
#' @details
#' The function generates pseudo-observations for a selected random effect using
#' an eigen-decomposition of the corresponding \code{psi} matrix. The pseudo-data
#' are then combined with the original dataset, zero-padding other random-effects
#' design matrices to maintain the correct dimensions. This augmented dataset can
#' be directly used for penalized GLMM fitting.
#'
#' @examples
#' data(birds)
#' Y<-birds$parasites
#' X <- model.matrix(parasites~migration+food, birds)
#' Z1<-model.matrix(parasites~migration, birds)
#' Z2<-model.matrix(parasites~1, birds)
#' grouping1<-as.numeric(as.factor(birds$Phylogenetic.Tomi))
#' grouping2<-as.numeric(as.factor(birds$species))
#' xdf<-list(Y=Y,X=X,Z=list(Z1,Z2),
#'   grouping=list(grouping1,grouping2),
#'   M=rep(1,nrow(X)),W=rep(1,nrow(X)))
#'
#' xdft<-AUGglmmTMB:::restructure_data(xdf)
#' pd_re_1<-make_pd_re(xdft, nu=list(3), psi=list(diag(1,2,2)),
#'   link_fun = "logit", which_re = 1)
#'
#' lapply(pd_re_1,function(x) if(is.null(dim(x))) length(x) else dim(x) )
#'
#' pd_re_2<-make_pd_re(xdft, nu=list(NULL,1), psi=list(NULL,matrix(1,1,1)),
#'   link_fun = "logit", which_re = 2)
#'
#' lapply(pd_re_2,function(x) if(is.null(dim(x))) length(x) else dim(x) )
#'
#' @export


make_pd_re<-function(data,nu,psi,link_fun,which_re=1){

    link_f <- switch(link_fun,
                       logit = function(x) plogis(x),
                       probit = function(x) pnorm(x),
                       cloglog = function(x) 1 - exp(-exp(x)),
                       cauchit = function(x) atan(x)/pi+0.5,
                       loglog = function(x) exp(-exp(-x)),
                       stop("Unsupported link function")
  )
  pd1<-make_pseudo_data_rand(psi[[which_re]],nu[[which_re]],const=1e8,param="variance",inv_link_fun=link_f)



  Y<-c(data$Y,pd1$data$Y)
  X0<-matrix(0,nrow=length(pd1$data$Y),ncol=ncol(data$X))
  X<-rbind(data$X,X0)
  M<-c(data$M,rep(1,length(pd1$data$Y)))


  W<-c(data$W,pd1$data$nn)

  pd2<-list(Y=Y,M=M,X=X,W=W)


  for (i in 1:sum(grepl("Z",names(data)))){
    nmZ<-paste0("Z",i)
    nmgr<-paste0("grouping",i)
    if (i==which_re){
      Z<-rbind(data[[nmZ]],pd1$data$Z)

      grouping<-c(data[[nmgr]],max(data[[nmgr]])+pd1$data$grouping)

    } else {
    Z0<-matrix(0,nrow=nrow(pd1$data$Z),ncol=ncol(data[[nmZ]]))
    Z<-rbind(data[[nmZ]],Z0)

    grouping<-c(data[[nmgr]],max(data[[nmgr]])+pd1$data$grouping)
    }

    pd2[[nmZ]]<-Z
    pd2[[nmgr]]<-grouping
  }


  pd2
}


#' @title Deprecated: Create pseudo-data for GLM (logit only)
#' @description Internal function for creating pseudo-data for a GLM using
#'   the fixed-effects penalty. Only works for the logit link. This function
#'   is deprecated and should not be used in new code.
#'
#' @param data List. Original dataset with elements \code{Y}, \code{M}, \code{X},
#'   \code{W}, and \code{o}.
#' @param cfe Numeric. Fixed-effects penalty strength.
#' @param fix_beta Numeric vector. Current fixed-effects coefficient estimates.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{\code{Y}}{Weighted response vector including pseudo-data.}
#'     \item{\code{M}}{Weighted total counts vector including pseudo-data.}
#'     \item{\code{X}}{Design matrix for fixed effects including pseudo-data.}
#'     \item{\code{O}}{Offset vector including pseudo-data.}
#'   }
#' @keywords internal


make_pd_fix_glm_logit<-function(data,cfe,fix_beta){

  Y<-c(data$Y,data$Y,data$M-data$Y)
  X<-rbind(data$X,data$X,data$X)
  M<-c(data$M,data$M,data$M)

  pi0<-c(plogis(data$X%*%matrix(fix_beta,ncol=1)))
  W<-diag(data$W*data$M*pi0*(1-pi0))
  Xw <- (W %^% 0.5) %*% data$X
  qrXw <- qr(Xw)
  H <- tcrossprod(qr.Q(qrXw))
  hi<-diag(H)

  W<-c(data$W,hi*cfe/data$M,hi*cfe/data$M)

  #O<-c(data$o,rep(0,2*nrow(data$X)))
  O<-c(data$o,data$o,data$o)

  list(Y=Y*W,M=M*W,X=X,O=O)

}

#' @title Deprecated: Penalized GLM via Data Augmentation (logit only)
#' @description Internal function implementing a penalized GLM via data
#'   augmentation. Only works for the logit link. This approach is no longer
#'   used in the main workflow.
#'
#' @param data List. Original dataset with elements \code{Y}, \code{M}, \code{X},
#'   \code{W}, and \code{o}.
#' @param cfe Numeric. Fixed-effects penalty strength (as in Kosuta et al.).
#' @param maxIter Integer. Maximum number of iterations for the iterative GLM fitting. Default is 15.
#' @param tol Numeric. Convergence tolerance for coefficient updates. Default is 1e-9.
#' @param link_fun Character. Link function to use. Only "logit" is supported.
#' @param save_coef Logical. If TRUE, returns coefficient paths across iterations. Default is FALSE.
#'
#' @return If \code{save_coef = FALSE}, returns a \code{glm} object. If
#'   \code{save_coef = TRUE}, returns a list with elements:
#'   \describe{
#'     \item{\code{fit}}{The fitted \code{glm} object.}
#'     \item{\code{coefs}}{Matrix of coefficient estimates across iterations.}
#'   }
#'
#' @keywords internal

my_pen_glm<-function(data,cfe,maxIter=15,tol=1e-9,link_fun="logit",save_coef=FALSE){
  if (link_fun!="logit") stop("Only logit link is supported")

  beta_fix<-rep(0,ncol(data$X))
  flag=TRUE
  ii=0
  while(flag==TRUE&ii<maxIter){
  ii=ii+1
  xdf<-make_pd_fix_glm_logit(data,cfe,beta_fix)
  fit<-glm(cbind(Y,M-Y)~-1+X+offset(O),data=xdf,family=binomial(link=link_fun))
  beta_fix<-fit$coefficients
  if (ii>1){
    est<-fit$coefficients

    if (max(abs(est-est0))<tol) flag=FALSE

  }
  est0<-fit$coefficients
  if (save_coef){
    if (ii==1) coefs<-est0 else coefs<-rbind(coefs,est0)
  }
  }
  if (save_coef) return(list(fit,coefs)) else return(fit)

}





#' Fit a penalized binomial mixed model via iterative algorithms
#'
#' Fits a penalized binomial mixed model (GLMM) using Algorithm 1 or Algorithm 2 presented in Košuta et al.
#'
#' @param data A list containing the model data with elements:
#' \describe{
#'   \item{Y}{Response vector.}
#'   \item{M}{Number of trials. Can be missing in which case it is internally set to 1 for all observations.}
#'   \item{W}{Likelihood, frequency or precision weights. Can be missing in which case it is internally set to 1 for all observations.}
#'   \item{X}{Design matrix for fixed effects.}
#'   \item{Z}{Random-effects design matrix or a list of matrices for multiple
#'   random-effects terms.}
#'   \item{grouping}{Grouping variable or a list of grouping variables corresponding
#'   to each random-effects term.}
#' }
#' If multiple random effects are present, \code{Z} and \code{grouping} must be
#' lists of equal length.
#'
#' @param cfe Numeric. Strength of the penalty applied to the fixed effects
#' (as in Košuta et al.).
#'
#' @param nu,psi Penalty parameters for the random effects. If \code{NULL},
#' no random-effects penalty is applied and only the fixed-effects penalty is used.
#' If provided, these have to be lists, that can be of the same length as \code{Z}, in which case
#' the penalty is applied to all random effects. If shorter, the penalty is applied
#' only to the first \code{length(nu)} random-effects terms.
#'
#' @param fit_pGLM Logical. If \code{TRUE}, uses Algorithm 2; if \code{FALSE}, uses Algorithm 1. Default is \code{FALSE}.
#'
#' @param maxiter Integer. Maximum number of iterations. Default is \code{50}
#'
#' @param tol Numeric. Convergence tolerance based on changes in the parameter
#' vector (fixed effects and variance components). Default is \code{1e-6}.
#'
#' @param link_fun Character. Link function to use. Supported options include
#' \code{"logit"}, \code{"probit"}, \code{"cloglog"}, \code{"cauchit"}, and
#' \code{"loglog"}. Note that \code{"cauchit"} and \code{"loglog"} are not
#' available in \pkg{glmmTMB}, although they are supported internally. Default is \code{"logit"}.
#'
#' @param save_coef Logical. If \code{TRUE}, stores the parameter estimates
#' (fixed effects and variance components) at each iteration. Default is \code{FALSE}.
#'
#' @param inter_iter Integer. Number of internal iterations used by
#' \pkg{glmmTMB}. Default is \code{1e3}.
#'
#' @param use_previous Logical. If \code{TRUE}, uses parameter estimates from the
#' previous iteration as starting values. Default is \code{TRUE}.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{fit}{The fitted penalized GLMM object.}
#'   \item{loglik}{The unpenalized marginal log-likelihood.}
#'   \item{coefs}{Matrix of GLMM parameter estimates by iteration.}
#'   \item{coefs_glm}{Matrix of GLM parameter estimates by iteration (only relevant
#'   when \code{fit_pGLM = TRUE}).}
#' }
#'
#' @details
#' This function is used internally by \code{\link{AUGglmmTMB}} and is not meant to be used directly (it can be, if used correctly, however (much) faster).
#'
#' The function implements penalized likelihood estimation for binomial mixed model with optional
#' penalties on both fixed and random effects. Two algorithms are available:
#' Algorithm 1 and Algorithm 2. The single-iteration estimator from Košuta et el. is obtained by specifying \code{fit_pGLM = TRUE,maxiter=1}.
#'
#' @seealso \code{\link{AUGglmmTMB}}, \code{\link{get_psi}}, \code{\link{get_data_plot_cloglik}}
#'
#' @examples
#' \dontrun{
#'
#'
#' data(birds)
#' Y<-birds$parasites
#' X <- model.matrix(parasites~migration+food, birds)
#' Z1<-model.matrix(parasites~migration, birds)
#' Z2<-model.matrix(parasites~1, birds)
#' grouping1<-as.numeric(as.factor(birds$Phylogenetic.Tomi))
#' grouping2<-as.numeric(as.factor(birds$species))
#' xdf<-list(Y=Y,X=X,Z=list(Z1,Z2),
#'   grouping=list(grouping1,grouping2),
#'   M=rep(1,nrow(X)),W=rep(1,nrow(X)))
#'
#' #MPL-A1
#'
#' fit_mpl1<-mpl_fitter(xdf,mv_multiplier(xdf),nu=list(3),psi=list(diag(1,2,2)),
#'   fit_pGLM=FALSE,maxiter=50,tol=1e-6,link_fun="logit",
#'     save_coef=TRUE,inter_iter=1e3,use_previous=TRUE)
#' summary(fit_mpl1$fit)
#' fit_mpl1$coefs
#'
#' #MPL-A2
#'
#' fit_mpl2<-mpl_fitter(xdf,mv_multiplier(xdf),nu=list(3),psi=list(diag(1,2,2)),
#'   fit_pGLM=TRUE,maxiter=50,tol=1e-6,link_fun="logit",
#'     save_coef=TRUE,inter_iter=1e3,use_previous=TRUE)
#' summary(fit_mpl2$fit)
#' fit_mpl2$coefs
#'
#'#' #MPL-A2(1)
#'
#' fit_mpl21<-mpl_fitter(xdf,mv_multiplier(xdf),nu=list(3),psi=list(diag(1,2,2)),
#'   fit_pGLM=TRUE,maxiter=1,tol=1e-6,link_fun="logit",
#'     save_coef=TRUE,inter_iter=1e3,use_previous=TRUE)
#' summary(fit_mpl21$fit)
#' fit_mpl21$coefs
#'
#' #MPL-A1: penalty on bith REs
#'
#' fit_mpl1_2repen<-mpl_fitter(xdf,mv_multiplier(xdf),nu=list(3,1),psi=list(diag(1,2,2),diag(1,1,1)),
#'   fit_pGLM=FALSE,maxiter=50,tol=1e-6,link_fun="logit",
#'     save_coef=TRUE,inter_iter=1e3,use_previous=TRUE)
#' summary(fit_mpl1_2repen$fit)
#' fit_mpl1_2repen$coefs
#' }
#' @export


mpl_fitter<-function(data,cfe,nu,psi,fit_pGLM=FALSE,maxiter=50,tol=1e-6,link_fun="logit",save_coef=FALSE,inter_iter=1e3,use_previous=TRUE){
  if (link_fun=="cauchit") stop("Cauchit link is currently not supported by glmmTMB")
  if (link_fun=="loglog") stop("Loglog link is currently not supported by glmmTMB")
  if (is.null(data$Y)) stop("data must be a list, element Y is missing")
  if (is.null(data$M)) data$M<-rep(1,length(data$Y)) #stop("data must be a list, element M is missing")
  if (is.null(data$X)) stop("data must be a list, element X is missing")
  if (is.null(data$Z)) stop("data must be a list, element Z is missing")
  if (is.null(data$W)) data$W<-rep(1,length(data$Y)) #stop("data must be a list, element W is missing")
  if (is.null(data$grouping)) stop("data must be a list, element grouping is missing")
if (!is.null(nu)) if (!is.list(nu)) stop("Arguments nu and psi need to be lists.")
  formula_fit<-make_formula(data)
  data<-restructure_data(data)
  if (cfe==0) maxiter<-1
  flag=TRUE
  ii=0
  while (flag==TRUE&ii<(maxiter)){
  ii=ii+1
  if (fit_pGLM) {
    if (ii==1) data$o<-rep(0,length(data$Y))
    data2<-data
    data2$Y<-data2$Y*data2$W
    data2$M<-data2$M*data2$W
    #firth<-glm(cbind(Y,M-Y)~-1+X+offset(o),data=data2,family=binomial(link=link_fun),method="brglmFit",type="MPL_Jeffreys",a=cfe)
    firth<-glm(cbind(Y,M-Y)~-1+X+offset(o),data=data2,family=binomial(link=link_fun),method=brglm2::brglmFit,type="MPL_Jeffreys",a=cfe)
    #firth<-my_pen_glm(data2,cfe,maxIter=15,tol=1e-9,link_fun="logit",save_coef=FALSE)

    beta_fix<-firth$coefficients
   # beta_fix[beta_fix>1e5]<-1e5
  #  beta_fix[beta_fix<(-1e5)]<--1e5
    if (save_coef) {if (ii==1) coefs_glm<-beta_fix else coefs_glm<-rbind(coefs_glm,beta_fix)}
  } else {

 if (ii==1) beta_fix<-rep(0,ncol(data$X))
  }

  xdf<-make_pd_fix(data,cfe,beta_fix,link_fun)

  if (!is.null(nu)){
    for (ll in 1:length(nu)){
    xdf<-make_pd_re(xdf,nu,psi,link_fun,ll)
    }
  }

  xdf$Y<-xdf$Y*xdf$W
  xdf$M<-xdf$M*xdf$W


 #if (use_previous){
#  if (ii==1){
 #   gmm.fit<-glmmTMB(formula_fit, data = xdf, family = binomial(link=link_fun),control=glmmTMBControl(optCtrl=list(iter.max=inter_iter,eval.max=inter_iter),profile = FALSE,collect = FALSE), se = FALSE, verbose = FALSE, doFit = TRUE, REML = FALSE)
#  } else {
 #   gmm.fit<-glmmTMB(formula_fit, data = xdf,
  #                   start=list(beta=gmm.fit$fit$par[names(gmm.fit$fit$par)%in%"beta"],theta=gmm.fit$fit$par[names(gmm.fit$fit$par)%in%"theta"]),
   #                  family = binomial(link=link_fun),control=glmmTMBControl(optCtrl=list(iter.max=inter_iter,eval.max=inter_iter),profile = FALSE,collect = FALSE), se = FALSE, verbose = FALSE, doFit = TRUE, REML = FALSE)

  #}
 #} else {
#   gmm.fit<-glmmTMB(formula_fit, data = xdf, family = binomial(link=link_fun),control=glmmTMBControl(optCtrl=list(iter.max=inter_iter,eval.max=inter_iter),profile = FALSE,collect = FALSE), se = FALSE, verbose = FALSE, doFit = TRUE, REML = FALSE)

 #}

  start_list <- NULL
  if (use_previous && ii > 1) {
    start_list <- list(
      beta  = gmm.fit$fit$par[names(gmm.fit$fit$par) %in% "beta"],
      theta = gmm.fit$fit$par[names(gmm.fit$fit$par) %in% "theta"]
    )
  }

  # Fit the model
  gmm.fit <- glmmTMB(
    formula_fit,
    data = xdf,
    family = binomial(link = link_fun),
    start = start_list,
    control = glmmTMBControl(
      optCtrl = list(iter.max = inter_iter, eval.max = inter_iter),
      profile = FALSE,
      collect = FALSE
    ),
    se = FALSE,
    verbose = FALSE,
    doFit = TRUE,
    REML = FALSE
  )



  if (fit_pGLM){
    pi1<-predict(gmm.fit,newdata=data)
    pi0<-predict(gmm.fit,newdata=data,re.form = NA)
    data$o<-pi1-pi0
  } else {
    beta_fix<-fixef(gmm.fit)$cond
  }


  if (ii>1){
  est<-gmm.fit$fit$par

    if (max(abs(est-est0))<tol) flag=FALSE

  }
  est0<-gmm.fit$fit$par
  if (save_coef){
    if (ii==1) coefs<-est0 else coefs<-rbind(coefs,est0)
  }

  }#end while



    if (save_coef&!fit_pGLM) coefs_glm<-NA
    if (!save_coef) {coefs_glm<-NA;coefs<-NA}





  #gmm.fit
    #refit to ogrig + pd Z

    if (!is.null(nu)){


      xdfa<-make_pd_re(data,nu,psi,link_fun)

      gmm.fit2<-glmmTMB(gmm.fit$call$formula,data=xdfa,family=binomial(link=link_fun),
                        start=list(beta=gmm.fit$fit$par[names(gmm.fit$fit$par)%in%"beta"],theta=gmm.fit$fit$par[names(gmm.fit$fit$par)%in%"theta"]),
                        control = glmmTMBControl(optCtrl = list(iter.max=0, eval.max=0),rank_check ="skip",conv_check="skip")
      )





    }


    #refit to orig data only
    gmm.fit3<-glmmTMB(gmm.fit$call$formula,data=data,family=binomial(link=link_fun),
                      start=list(beta=gmm.fit$fit$par[names(gmm.fit$fit$par)%in%"beta"],theta=gmm.fit$fit$par[names(gmm.fit$fit$par)%in%"theta"]),
                      control = glmmTMBControl(optCtrl = list(iter.max=0, eval.max=0),rank_check ="skip",conv_check="skip")
    )



    if (!is.null(nu))  list(fit=gmm.fit2,loglik=-gmm.fit3$fit$objective,coefs=coefs,coefs_glm=coefs_glm) else list(fit=gmm.fit3,loglik=-gmm.fit3$fit$objective,coefs=coefs,coefs_glm=coefs_glm)




}









#' Estimate penalty parameters for random effects
#'
#' Computes the penalty parameters \eqn{\tau} and \eqn{\Psi} for the random-effects
#' penalty, following the approach of Košuta et al. The penalty is assumed to be
#' applied only to the first random-effects term.
#'
#' @param data A list containing the model data. See \code{mpl_fitter} for details.
#' @param cfe Numeric. Strength of the penalty applied to the fixed effects.
#' @param nu Penalty parameter for the random effects. See \code{mpl_fitter}
#' for details.
#' @param fit_pGLM Logical. Indicates whether Algorithm 2 (\code{TRUE}) or Algorithm 1 (\code{FALSE}) is used. Default is \code{FALSE}.
#' @param maxiter Integer. Maximum number of iterations. Default is \code{50}.
#' @param tol Numeric. Convergence tolerance. Default is \code{1e-6}.
#' @param link_fun Character. Link function. Default is \code{"logit"}.
#' @param save_coef Logical. Whether to store coefficient paths. Default is \code{FALSE}.
#' @param inter_iter Integer. Number of internal iterations used by \pkg{glmmTMB}. Default is \code{1e3}.
#' @param use_previous Logical. Whether to use previous estimates as starting values. Default is \code{TRUE}.
#'
#' @return A list with elements:
#' \describe{
#'   \item{tau}{Estimated value of \eqn{\tau}.}
#'   \item{psi}{Estimated value of \eqn{\Psi}.}
#' }
#'
#' @details
#' This function is used internally by \code{\link{AUGglmmTMB}} and is not meant to be used directly (it can be, if used correctly, however (much) faster).
#'
#' The procedure assumes that the random-effects penalty is imposed only on the
#' first random-effects component. The estimates are obtained using an iterative
#' procedure consistent with the penalized likelihood framework described in
#' Košuta et al.
#'
#' @seealso \code{\link{AUGglmmTMB}}, \code{\link{get_data_plot_cloglik}}, \code{\link{mpl_fitter}}
#'
#' @examples
#' \dontrun{
#'
#' data(birds)
#' Y<-birds$parasites
#' X <- model.matrix(parasites~migration+food, birds)
#' Z1<-model.matrix(parasites~migration, birds)
#' Z2<-model.matrix(parasites~1, birds)
#' grouping1<-as.numeric(as.factor(birds$Phylogenetic.Tomi))
#' grouping2<-as.numeric(as.factor(birds$species))
#' xdf<-list(Y=Y,X=X,Z=list(Z1,Z2),
#'   grouping=list(grouping1,grouping2),
#'   M=rep(1,nrow(X)),W=rep(1,nrow(X)))
#'
#' psi_opt<-get_psi(xdf,mv_multiplier(xdf),nu=list(3),
#'   fit_pGLM=TRUE,maxiter=1,tol=1e-6,link_fun="logit",
#'     save_coef=FALSE,inter_iter=1e3,use_previous=FALSE)
#' psi_opt
#'
#' }
#' @export


get_psi<-function(data,cfe,nu,fit_pGLM=FALSE,maxiter=50,tol=1e-6,link_fun="logit",save_coef=FALSE,inter_iter=1e3,use_previous=FALSE){

  fit0<-mpl_fitter(data,cfe=cfe,nu=NULL,psi=NULL,fit_pGLM=fit_pGLM,maxiter=maxiter,tol=tol,link_fun=link_fun,save_coef = FALSE,use_previous=use_previous)

  D_est<-VarCorr(fit0$fit)$cond[[1]]


  get_cond_lik<-function(model,Y){

    p_hat<-predict(model,type="response")[1:length(Y)]

    p1<-log(p_hat)
    p2<-log(1-p_hat)

    p1[p1<(-1e5)]<-(-1e5)
    p2[p2<(-1e5)]<-(-1e5)

    sum(  Y*p1+(1-Y)*p2  )

  }





  get_lim<-function(model,Y,X){


    p_hat<-predict(model,type="response")[1:length(Y)]

    dl_deta<-Y-p_hat


    grad_beta <- t(X) %*% dl_deta


    vcov_beta <- vcov(model)$cond
    if (any(is.na(vcov_beta))) vcov_beta<-4*diag(1,ncol(X),ncol(X))
    if (any(diag(vcov_beta<0))) vcov_beta<-4*diag(1,ncol(X),ncol(X)) #res from Kos
    sqrt(t(grad_beta)%*%vcov_beta%*%grad_beta)


  }



  q<-nrow(D_est)
  tau=0
  ee<-eigen(D_est)

  lm<-1 #identity as the shrink target
  li<-ee$values+tau*(lm-ee$values)
  if (q>1) psi0<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q else psi0<-ee$vectors%*%diag(li,1,1)%*%t(ee$vectors)*3*q




  fit_rok1<-mpl_fitter(data,cfe=cfe,nu=nu,psi=list(psi0),fit_pGLM=fit_pGLM,maxiter=maxiter,tol=tol,link_fun=link_fun,save_coef = FALSE,use_previous=use_previous)


  lik0<-get_cond_lik(fit_rok1$fit,data$Y)





  lim<-get_lim(fit_rok1$fit,data$Y,data$X) #1SE rule




  tau_finder_maxi<-function(tau){
    li<-ee$values+tau*(lm-ee$values)
    if (q>1) psi0<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q else psi0<-ee$vectors%*%diag(li,1,1)%*%t(ee$vectors)*3*q


    fit_rok<-mpl_fitter(data,cfe=cfe,nu=nu,psi=list(psi0),fit_pGLM=fit_pGLM,maxiter=maxiter,tol=tol,link_fun=link_fun,save_coef = FALSE,use_previous=use_previous)
    -get_cond_lik(fit_rok$fit,data$Y)

  }

  tau_finder<-function(tau){
    li<-ee$values+tau*(lm-ee$values)
    if (q>1) psi0<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q else psi0<-ee$vectors%*%diag(li,1,1)%*%t(ee$vectors)*3*q


    fit_rok<-mpl_fitter(data,cfe=cfe,nu=nu,psi=list(psi0),fit_pGLM=fit_pGLM,maxiter=maxiter,tol=tol,link_fun=link_fun,save_coef = FALSE,use_previous=use_previous)
    abs(get_cond_lik(fit_rok$fit,data$Y)-lik0)-lim

  }

  tol_min<-.Machine$double.eps^0.2
  tol_max<-.Machine$double.eps^0.2

  opti<-optimise(tau_finder_maxi,c(0,1))

  if (opti$minimum<tol_min){
    uni<-  try(uniroot(tau_finder,lower=0,upper = 1),silent=TRUE)
    if (inherits(uni, "try-error")) {
      opt_tau <- 1
    } else {
      opt_tau <- uni$root
    }
  } else {
    if ((-opti$objective-lik0)<lim) opt_tau<-opti$minimum else {
      uni<-  uniroot(tau_finder,lower=tol_min,upper = opti$minimum)
      opt_tau<-uni$root
    }
  }
  if (opt_tau>(1-tol_max)) opt_tau<-1-tol_max
  if (opt_tau<tol_min) opt_tau<-tol_min


  tau1<-opt_tau

  li<-ee$values+opt_tau*(lm-ee$values)
  if (q>1) psi0<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q else psi0<-ee$vectors%*%diag(li,1,1)%*%t(ee$vectors)*3*q

  list(tau=tau1,psi=psi0)

}



#' Conditional log-likelihood profile for random-effects shrinkage
#'
#' Computes the conditional log-likelihood over a grid
#' of shrinkage parameters (\eqn{\tau}) applied to the random-effects covariance
#' matrix. The function constructs a sequence of covariance matrices using
#' eigenvalue shrinkage toward the identity matrix, refits the model for each
#' value of \eqn{\tau}, and evaluates the corresponding conditional log-likelihood.
#' A standard-error-based threshold is also returned for model selection purposes
#' (e.g., the 1-SE rule).
#'
#' @param data A list containing the data components as required by
#'   \code{mpl_fitter}, including \code{Y}, \code{M}, \code{W}, \code{X},
#'   and random-effects structures (\code{Z}, \code{grouping}).
#' @param cfe Numeric. Strength of the fixed-effects penalty.
#' @param nu Numeric. Parameter controlling the random-effects penalty.
#' @param n_taus Integer. Number of grid points for the shrinkage parameter
#'   \eqn{\tau} in the interval [0, 1].
#' @param fit_pGLM Logical. If \code{TRUE}, uses Algorithm 2; otherwise uses Algorithm 1.
#' @param maxiter Integer. Maximum number of iterations for the fitting procedure.
#' @param tol Numeric. Convergence tolerance for parameter estimates.
#' @param link_fun Character. Link function to use. Options include
#'   \code{"logit"}, \code{"probit"}, \code{"cloglog"}, \code{"cauchit"},
#'   and \code{"loglog"}.
#' @param save_coef Logical. Whether to store coefficient paths across iterations.
#' @param inter_iter Integer. Number of internal iterations used by
#'   \code{glmmTMB}.
#' @param use_previous Logical. Whether to use previous estimates as starting
#'   values in iterative updates.
#'
#' @return A list with components:
#' \describe{
#'   \item{taus}{Numeric vector of shrinkage parameter values.}
#'   \item{cloglik}{Numeric vector of corresponding conditional log-likelihood values.}
#'   \item{se}{Numeric value representing a standard-error-based threshold.}
#' }
#'
#' @details
#' This function is used internally by \code{\link{AUGglmmTMB}} and is not meant to be used directly (it can be, if used correctly, however (much) faster).
#'
#' The procedure starts by fitting a model without a random-effects penalty to obtain
#' an initial estimate of the covariance matrix. A sequence of covariance matrices is
#' then generated by shrinking the eigenvalues toward the identity matrix. For each
#' value of \eqn{\tau}, the model is refitted and the conditional log-likelihood is
#' evaluated.
#'
#' The reported standard error is based on a quadratic approximation using the
#' gradient of the conditional log-likelihood and the estimated covariance matrix
#' of the fixed effects, and can be used in conjunction with the 1-SE rule for
#' selecting \eqn{\tau}.
#'
#' @seealso \code{\link{AUGglmmTMB}}, \code{\link{get_psi}}, \code{\link{mpl_fitter}}
#'
#' @examples
#' \dontrun{
#'
#' data(birds)
#' Y<-birds$parasites
#' X <- model.matrix(parasites~migration+food, birds)
#' Z1<-model.matrix(parasites~migration, birds)
#' Z2<-model.matrix(parasites~1, birds)
#' grouping1<-as.numeric(as.factor(birds$Phylogenetic.Tomi))
#' grouping2<-as.numeric(as.factor(birds$species))
#' xdf<-list(Y=Y,X=X,Z=list(Z1,Z2),
#'   grouping=list(grouping1,grouping2),
#'   M=rep(1,nrow(X)),W=rep(1,nrow(X)))
#'
#' cloglik<-get_data_plot_cloglik(xdf,mv_multiplier(xdf),nu=list(3),n_taus=20,
#'   fit_pGLM=TRUE,maxiter=1,tol=1e-6,link_fun="logit",
#'     save_coef=FALSE,inter_iter=1e3,use_previous=FALSE)
#'
#' psi_opt<-get_psi(xdf,mv_multiplier(xdf),nu=list(3),
#'   fit_pGLM=TRUE,maxiter=1,tol=1e-6,link_fun="logit",
#'     save_coef=FALSE,inter_iter=1e3,use_previous=FALSE)
#'
#' plot(cloglik$taus,cloglik$cloglik,type="l",xlab="tau",ylab="conditional log-likelihood")
#' abline(h=cloglik$cloglik[1]+cloglik$se,lty=3,col="red")
#' abline(v=psi_opt$tau,lty=3,col="blue")
#'
#' }
#' @export

get_data_plot_cloglik<-function(data,cfe,nu,n_taus,fit_pGLM=FALSE,maxiter=50,tol=1e-6,link_fun="logit",save_coef=FALSE,inter_iter=1e3,use_previous=FALSE){

  fit0<-mpl_fitter(data,cfe=cfe,nu=NULL,psi=NULL,fit_pGLM=fit_pGLM,maxiter=maxiter,tol=tol,link_fun=link_fun,save_coef = FALSE,use_previous=use_previous)

  D_est<-VarCorr(fit0$fit)$cond[[1]]


  get_cond_lik<-function(model,Y){

    p_hat<-predict(model,type="response")[1:length(Y)]

    p1<-log(p_hat)
    p2<-log(1-p_hat)

    p1[p1<(-1e5)]<-(-1e5)
    p2[p2<(-1e5)]<-(-1e5)

    sum(  Y*p1+(1-Y)*p2  )

  }




  get_lim<-function(model,Y,X){


    p_hat<-predict(model,type="response")[1:length(Y)]

    dl_deta<-Y-p_hat


    grad_beta <- t(X) %*% dl_deta


    vcov_beta <- vcov(model)$cond
    if (any(is.na(vcov_beta))) vcov_beta<-4*diag(1,ncol(X),ncol(X))
    if (any(diag(vcov_beta<0))) vcov_beta<-4*diag(1,ncol(X),ncol(X)) #res from Kos
    sqrt(t(grad_beta)%*%vcov_beta%*%grad_beta)


  }



  q<-nrow(D_est)
  tau=0
  ee<-eigen(D_est)

  lm<-1 #identity as the shrink target
  li<-ee$values+tau*(lm-ee$values)
  if (q>1) psi0<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q else psi0<-ee$vectors%*%diag(li,1,1)%*%t(ee$vectors)*3*q




  fit_rok1<-mpl_fitter(data,cfe=cfe,nu=nu,psi=list(psi0),fit_pGLM=fit_pGLM,maxiter=maxiter,tol=tol,link_fun=link_fun,save_coef = FALSE,use_previous=use_previous)


  lik0<-get_cond_lik(fit_rok1$fit,data$Y)





  lim<-get_lim(fit_rok1$fit,data$Y,data$X) #1SE rule




  tau_finder_maxi<-function(tau){
    li<-ee$values+tau*(lm-ee$values)
    if (q>1) psi0<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q else psi0<-ee$vectors%*%diag(li,1,1)%*%t(ee$vectors)*3*q


    fit_rok<-mpl_fitter(data,cfe=cfe,nu=nu,psi=list(psi0),fit_pGLM=fit_pGLM,maxiter=maxiter,tol=tol,link_fun=link_fun,save_coef = FALSE,use_previous=use_previous)
    -get_cond_lik(fit_rok$fit,data$Y)

  }

seq_tau<-seq(from=0,to=1,length.out=n_taus)
seq_tau<-seq_tau[-1]
c_lik<-rep(NA,length(seq_tau))
zz=0
for (ii in seq_tau) {
  zz=zz+1
  c_lik[zz]<-tau_finder_maxi(ii)
}

seq_tau<-c(0,seq_tau)
c_lik<-c(lik0,-c_lik)

list(taus=seq_tau,cloglik=c_lik,se=lim)

}



#' Control options for AUGglmmTMB fitting
#'
#' Create a list of control parameters to customize the fitting behavior of
#' \code{\link{AUGglmmTMB}}.
#'
#' @param fit_pGLM Logical. If \code{TRUE}, uses the penalized GLM step (Algorithm 2
#'   in Košuta et al.), otherwise it uses Algorithm 1. Default is \code{FALSE}.
#' @param maxiter Integer. Maximum number of iterations for the penalized likelihood
#'   fitting procedure. Default is 50.
#' @param tol Numeric. Convergence tolerance based on changes in parameter estimates.
#'   Default is 1e-6.
#' @param save_coef Logical. If \code{TRUE}, stores coefficient estimates at each
#'   iteration. Default is \code{FALSE}.
#' @param inter_iter Integer. Number of internal iterations used by
#'   \pkg{glmmTMB}. Default is 1000.
#' @param use_previous Logical. If \code{TRUE}, uses estimates from the previous
#'   iteration as starting values. Default is \code{FALSE}.
#'
#' @return A named list containing all specified control parameters.
#'
#' @examples
#'
#' #Algorithm 1 from Košuta et al.
#' ctrl <- AUGglmmTMBControl(fit_pGLM = FALSE)
#'
#' #Algorithm 2 from Košuta et al.
#' ctrl2 <- AUGglmmTMBControl(fit_pGLM = TRUE)
#'
#' #single-iteration approximation to MPL from Košuta et al.
#' ctrl3 <- AUGglmmTMBControl(fit_pGLM = TRUE,maxiter=1)
#'
#' @export

AUGglmmTMBControl <- function(fit_pGLM   = FALSE,
                              maxiter    = 50,
                              tol        = 1e-6,
                              save_coef  = FALSE,
                              inter_iter = 1e3,
                              use_previous = FALSE) {

  namedList(fit_pGLM,maxiter,tol,save_coef,inter_iter,use_previous)
}

#' Create penalty options for AUGglmmTMB
#'
#' Constructs a list of penalty parameters to be used with \code{\link{AUGglmmTMB}}.
#' This includes penalties for fixed effects and optional penalties for random effects.
#'
#' @param cfe Numeric. Strength of penalty on fixed effects. Default is \code{NULL},
#'   in which case \code{AUGglmmTMB} computes it internally via \code{\link{mv_multiplier}} as suggested by Košuta et al.
#'   Setting \code{cfe = 0} disables the fixed-effects penalty.
#' @param autrepen Logical. If \code{TRUE}, automatically estimates the random-effects
#'       penalty parameters \eqn{\tau} and \eqn{\Psi} using the procedure described in Košuta et al (see also \code{\link{get_psi}}). Default is \code{FALSE}.
#'       When \code{TRUE}, the penalty is only applied to the 1st random effect.
#' @param nu Optional numeric vector specifying random-effects penalty parameters. Ignored when \code{autrepen=TRUE}. Default is \code{NULL} in which case the random effects are not penalized.
#' @param psi Optional list of random-effects penalty matrices. Ignored when \code{autrepen=TRUE}. Ignored when \code{nu=NULL}; when \code{nu} is not \code{NULL}, \code{psi} has to be the list of the same length as \code{nu} containing matrices of appropriate dimensions. Default is \code{NULL}.
#' @param plot Logical flag for plotting the conditional likelihood. Ignored when \code{autrepen=FALSE}. When \code{TRUE}, the computation time can be substantially increased. Default is \code{FALSE}.
#' @param ntaus Numeric, number of different values of \eqn{\tau} to be used when plotting the conditional likelihood. Ignored when \code{plot=FALSE}. Using a large value can substantially increase the computation time. Default is \code{50}.
#'
#' @return A named list containing all specified penalty parameters.
#'
#' @details
#' This function does not perform any fitting itself. It simply returns a structured list
#' that specifies how \code{\link{AUGglmmTMB}} should apply penalties to fixed and random effects.
#' Use this function to easily create the \code{penOpt} argument for \code{\link{AUGglmmTMB}}.
#' Setting \code{autrepen=FALSE} and \code{nu=NULL} turns off the penalty on the random effects covariance matrices.
#'
#' @examples
#' # Default penalty on the fixed effects, no penalty on the random effects
#' pen <- AUGglmmTMBPenalty()
#'
#' # Fixed-effects penalty strength set to 0.5, data-driven random-effects penalty parameters
#' pen2 <- AUGglmmTMBPenalty(cfe = 0.5, autrepen = TRUE)
#'
#' # No fixed-effects penalty, data-driven random-effects penalty parameters
#' pen3 <- AUGglmmTMBPenalty(cfe = 0, autrepen = TRUE)
#'
#' # Turn off both penalties (the same as glmmTMB)
#' pen4 <- AUGglmmTMBPenalty(cfe = 0)
#'
#' @export


AUGglmmTMBPenalty<-function(cfe=NULL,
                            autrepen=FALSE,
                            nu=NULL,
                            psi=NULL,
                            plot=FALSE,ntaus=50){

  namedList(cfe,autrepen,nu,psi,plot,ntaus)
}



#' Fit a penalized binomial mixed model
#'
#' Fits a penalized binomial generalized linear mixed model with optional
#' penalties on both fixed and random effects using the approach of Košuta et al.
#'
#' @param formula A \code{\link{formula}} object specifying the fixed and random
#'   effects in standard R formula syntax. Random effects should be of the form
#'   \code{(expr | factor)}.
#' @param data A \code{data.frame} containing all variables used in \code{formula}. When missing values are present, the observations with any missing value are removed as in \code{\link{na.omit}}.
#' @param weights Optional. Column name, as in \code{\link{glmmTMB}} in \code{data} specifying
#'   observation weights. Default is \code{NULL}, in which case all weights are 1.
#' @param link Character. Link function to use for the binomial response. Default is
#'   \code{"logit"}.
#' @param penOpt List of penalty options, created by \code{\link{AUGglmmTMBPenalty}}:
#'   \describe{
#'     \item{\code{cfe}}{Numeric. Strength of penalty on fixed effects. If \code{NULL} (the default),
#'       internally computed via \code{\link{mv_multiplier}} as suggested by Košuta et al. Setting \code{cfe = 0} disables fixed-effect penalty.}
#'     \item{\code{autrepen}}{Logical. If \code{TRUE}, automatically estimates the random-effects
#'       penalty parameters \eqn{\tau} and \eqn{\Psi} using the procedure described in Košuta et al (see also \code{\link{get_psi}}). When \code{TRUE}, the penalty is only applied to the 1st random effect. Default is \code{FALSE}.}
#'     \item{\code{nu}}{Optional list specifying random-effects penalty parameters. Ignored when \code{autrepen=TRUE}. Can be of the same length as the number of specified random effects in which case the penalty is applied to all random effects; when shorter, the penalty is applied only to the first \code{length(nu)} random effects. Default is \code{NULL} in which case the random effects are not penalized.}
#'     \item{\code{psi}}{Optional list of random-effects penalty matrices. Ignored when \code{autrepen=TRUE}. Ignored when \code{nu=NULL}; when \code{nu} is not \code{NULL}, \code{psi} has to be the list of the same length as \code{nu} containing matrices of appropriate dimensions. Default is \code{NULL}.}
#'   \item{\code{plot}}{Logical flag for plotting the conditional likelihood. Ignored when \code{autrepen=FALSE}. When \code{TRUE}, the computation time can be substantially increased. Default is \code{FALSE}.}
#'   \item{\code{ntaus}}{Numeric, number of different values of \eqn{\tau} to be used when plotting the conditional likelihood. Ignored when \code{plot=FALSE}. Using a large value can substantially increase the computation time. Default is \code{50}.}
#'   }
#' @param control List of control parameters created via \code{\link{AUGglmmTMBControl}}:
#' \describe{
#' \item{\code{fit_pGLM}}{Logical. If \code{TRUE}, uses the penalized GLM step (Algorithm 2
#'   in Košuta et al.), otherwise it uses Algorithm 1. Default is \code{FALSE}.}
#' \item{\code{maxiter}}{Integer. Maximum number of iterations for the penalized likelihood
#'   fitting procedure. Default is 50.}
#' \item{\code{tol}}{Numeric. Convergence tolerance based on changes in parameter estimates.
#'   Default is 1e-6.}
#' \item{\code{save_coef}}{Logical. If \code{TRUE}, stores coefficient estimates at each
#'   iteration. Default is \code{FALSE}.}
#' \item{\code{inter_iter}}{Integer. Number of internal iterations used by
#'   \pkg{glmmTMB}. Default is 1000.}
#' \item{\code{use_previous}}{Logical. If \code{TRUE}, uses estimates from the previous
#'   iteration as starting values. Default is \code{FALSE}.}
#' }
#'
#' @return A list with elements:
#' \describe{
#'   \item{fit}{A list with the elements:
#'   \describe{
#'   \item{fit}{The fitted penalized GLMM object of class \code{\link{glmmTMB}}.}
#'   \item{loglik}{The unpenalized marginal log-likelihood.}
#'   \item{coefs}{Matrix of GLMM parameter estimates by iteration.}
#'   \item{coefs_glm}{Matrix of GLM parameter estimates by iteration (only relevant
#'   when \code{fit_pGLM = TRUE}).}
#' }}
#'   \item{optre}{List with elements \code{opt_tau} and \code{opt_psi} corresponding
#'     to estimated random-effects penalty parameters if \code{autrepen = TRUE}, otherwise \code{NULL}.}
#' }
#'
#' @details
#' Setting \code{autrepen=TRUE} uses the data-driven procedure proposed by Košuta et al. to determine the penalty parameters; the parameter \code{nu} is set to \eqn{2q-1} internally, any other value supplied in \code{nu} is ignored.
#' If \code{autrepen=FALSE} and \code{nu=NULL}, no random-effects penalty is applied and only the fixed-effects penalty is used.
#' When \code{autrepen=FALSE}, if provided, \code{nu} and \code{psi}, can be of the same length as the number of specified random effects, in which case
#' the penalty is applied to all random effects. If shorter, the penalty is applied
#' only to the first \code{length(nu)} random-effects terms. When specified, \code{nu} and \code{psi} need to be of the same length: each element of \code{nu} and \code{psi} are the penalty parameters \eqn{\tau} and \eqn{\Psi}, respectively, for the corresponding random effect. When specified, the elements of \code{psi} need to be matrices of appropriate dimensions.
#'
#'
#' @seealso \code{\link{AUGglmmTMBPenalty}},\code{\link{AUGglmmTMBControl}}, \code{\link{get_psi}}, \code{\link{mpl_fitter}}
#'
#' @examples
#' \dontrun{
#' data(birds)
#'
#' #default fixed-effects penalty, no penalty on the random effect;
#' #single-iteration approximation to MPL
#'
#' fit <- AUGglmmTMB(
#'   cbind(parasites, 1 - parasites) ~ migration + food + (1 | species),
#'   data = birds,
#'   weights = NULL,
#'   link = "logit",
#'   penOpt = AUGglmmTMBPenalty(),
#'   control = AUGglmmTMBControl(fit_pGLM = TRUE, maxiter = 1)
#' )
#' summary(fit$fit$fit)
#'
#' #the same
#'
#' fit2 <- AUGglmmTMB(
#'   parasites ~ migration + food + (1 | species),
#'   data = birds,
#'   weights = NULL,
#'   link = "logit",
#'   penOpt = AUGglmmTMBPenalty(),
#'   control = AUGglmmTMBControl(fit_pGLM = TRUE, maxiter = 1)
#' )
#' summary(fit2$fit$fit)
#'
#' #default penalty on the fixed effects, penalty on both random effects;
#' #single-iteration approximation to MPL
#'
#' fit3 <- AUGglmmTMB(
#'   parasites ~ migration + food + (migration|phylogenetic) + (1 | species),
#'   data = birds,
#'   weights = NULL,
#'   link = "logit",
#'   penOpt = AUGglmmTMBPenalty(nu=list(3,1),psi=list(diag(1,2,2),matrix(1,1,1))),
#'   control = AUGglmmTMBControl(fit_pGLM = TRUE, maxiter = 1)
#' )
#' summary(fit3$fit$fit)
#'
#' #default penalty on the fixed effects, data-driven penalty on phylogenetic, no penalty on species;
#' #single-iteration approximation to MPL
#'
#' fit4 <- AUGglmmTMB(
#'   parasites ~ migration + food + (migration|phylogenetic) + (1 | species),
#'   data = birds,
#'   weights = NULL,
#'   link = "logit",
#'   penOpt = AUGglmmTMBPenalty(autrepen=TRUE,plot=TRUE,ntaus =20),
#'   control = AUGglmmTMBControl(fit_pGLM = TRUE, maxiter = 1)
#' )
#'
#' fit4$optre
#' summary(fit4$fit$fit)
#' }
#'
#' @export

AUGglmmTMB<-function(formula,data,weights=NULL,link="logit",
                     penOpt=AUGglmmTMBPenalty(),
                     control=AUGglmmTMBControl()){


  data<-na.omit(data)
  formula<-as.formula(formula)



 # if (is.null(link)) link_fun<-"logit" else link_fun<-link
  link_fun<-link

  fix_formula <- nobars(formula)




  random_form <- findbars(formula)


  re_list <- NULL
  j=1
  for(i in random_form){
    random_parts <- i
    bar <- random_parts[[1]]
    lhs <- random_parts[[2]]
    rhs <- random_parts[[3]]

    if(length(as.character(rhs))>1) stop("Only random effects of type (REexpr1 | factor1) + (REexpr2 | factor2) + ... are supported. Expressions of type (REexpr | factor1:factor2) or (REexpr | factor1/factor2) are not supported.")

    res <- list(expr=as.formula(paste0("~", deparse(lhs))),
                gr=as.character(rhs))

    re_list[[j]] <- res
    j=j+1

  }

  rand_formula <- re_list



  X <- model.matrix(fix_formula, data)

  Z<-list()
  grouping<-list()
  for(i in 1:length(rand_formula)){
    Z[[i]]<- model.matrix(rand_formula[[i]][[1]], data)
    grouping[[i]]<- as.numeric(factor(data[, rand_formula[[i]][[2]]]))
  }





  # build model frame with weights evaluated in data




  mf <- stats::model.frame(
    fix_formula,
    data = data)#,weights = weights)



  # extract response as evaluated by model.frame
  Yraw <- stats::model.response(mf)

  # Case 1: matrix response (cbind)
  if (is.matrix(Yraw)) {
    if (ncol(Yraw) != 2) {
      stop("cbind response must have exactly two columns")
    }

    Y <- Yraw[, 1]
    M <- rowSums(Yraw)


  } else {
    Y <- Yraw
    M <- rep(1, length(Y))


  }

  #if (is.null(weights)){
  #  W <- rep(1, length(Y))
  #} else {

  #if (is.numeric(weights)){
  #  if (length(weights)!=nrow(data)) stop("The specified weights are not correct length.") else W<-weights
  #} else {
  #  if (is.character(weights)) W<-data[,weights] else {
  #    w_expr <- substitute(weights)
  #    W <- eval(w_expr, envir = data, enclos = parent.frame())
  #    }
  #}
  #}
  #W<- model.weights(mf)
  #if (is.null(W)) W <- rep(1, length(Y))
  mc <- match.call()

  # keep only the relevant arguments for model.frame
  mf <- mc[c(1L, match(c("data", "weights"), names(mc), 0L))]
  #mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")  # replace function name with model.frame

  # evaluate model.frame in the formula environment
  #fr <- eval(mf, envir = environment(formula))
  fr <- eval(mf, envir = parent.frame())
  # number of observations
  nobs <- nrow(fr)

  # extract weights safely
  W <- as.vector(model.weights(fr))
  if (is.null(W)) W <- rep(1, nobs)

  data_mpl<-namedList(Y,M,W,X,Z,grouping)

  if (is.null(penOpt$cfe)) penOpt$cfe<-mv_multiplier(data_mpl)

  if (penOpt$autrepen){
    q<-ncol(data_mpl$Z[[1]])
    psi_opt<-suppressMessages(suppressWarnings(get_psi(data=data_mpl,cfe=penOpt$cfe,nu=list(2*q-1),fit_pGLM=control$fit_pGLM,maxiter=control$maxiter,
                                                       tol=control$tol,link_fun=link_fun,
                     save_coef=FALSE,inter_iter=control$inter_iter,use_previous=control$use_previous)))

    opt_tau<-psi_opt$tau
    opt_psi<-psi_opt$psi

    fit<-suppressMessages(suppressWarnings(mpl_fitter(data=data_mpl,
                    cfe=penOpt$cfe,nu=list(2*q-1),psi=list(opt_psi),
                    fit_pGLM=control$fit_pGLM,maxiter=control$maxiter,tol=control$tol,link_fun=link_fun,
                    save_coef=control$save_coef,inter_iter=control$inter_iter,use_previous=control$use_previous)))
if (penOpt$plot){
  psi_taus<-suppressMessages(suppressWarnings(get_data_plot_cloglik(data=data_mpl,cfe=penOpt$cfe,nu=list(2*q-1),n_taus=penOpt$ntaus,fit_pGLM=control$fit_pGLM,maxiter=control$maxiter,
                                  tol=control$tol,link_fun=link_fun,save_coef=FALSE,inter_iter=control$inter_iter,use_previous=control$use_previous)))

  plot(psi_taus$taus,psi_taus$cloglik,type="l",xlab="tau",ylab="conditional log-likelihood")
  if (psi_taus$cloglik[1]<psi_taus$cloglik[2])  lm<-psi_taus$cloglik[1]+psi_taus$se else psi_taus$cloglik[1]-psi_taus$se
  abline(h=lm,lty=3,col="red")
   abline(v=psi_opt$tau,lty=3,col="blue")

}
  } else {
    opt_tau<-NULL
    opt_psi<-NULL
    fit<-suppressMessages(suppressWarnings(mpl_fitter(data=data_mpl,
                    cfe=penOpt$cfe,nu=penOpt$nu,psi=penOpt$psi,
                    fit_pGLM=control$fit_pGLM,maxiter=control$maxiter,tol=control$tol,link_fun=link_fun,
                    save_coef=control$save_coef,inter_iter=control$inter_iter,use_previous=control$use_previous)))
  }

  list(fit=fit,optre=list(tau=opt_tau,psi=opt_psi))

}
