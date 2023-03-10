#' @title Peuso Maximum Likelihood (PML) Estimator of SAR model with Partial Network Data
#' @description \code{qmleSAR} implements a maximum likelihood-based  estimator of the linear-in-mean SAR model when only the linking probabilities are available or can be estimated.
#' @param formula object of class \link[stats]{formula}: a symbolic description of the model. The `formula` should be as for example \code{y ~ x1 + x2 | gy | gx1 + gx2}
#' where `y` is the endogenous vector, the listed variables before the pipe, `x1`, `x2` are the individual exogenous variables, `gy` is the average of `y` among friends, and
#' `gx1`, `gx2` are the contextual observed variables. If `gy` is observed and `gx1`, `gx2` are not, the formula should be
#' \code{y ~ x1 + x2 | gy}. If `gy` is not observed and `gx1`, `gx2` are, the formula should be \code{y ~ x1 + x2 || gx1 + gx2}. If `gy`, `gx1`, and `gx2` are not observed, the 
#' the formula should simply be \code{y ~ x1 + x2}.
#' @param  contextual logical; if true, this means that all individual variables will be set as contextual variables. In contrast \code{\link{mcmcSAR}},
#' `formula` as `y ~ x1 + x2` and `contextual` as `TRUE` is not equivalent to set formula as `y ~ x1 + x2 || gx1 + gx2`. `formula = y ~ x1 + x2` means that `gy`, `gx1`, and `gx2` 
#' are not observed and `contextual = TRUE` means that the estimated model includes contextual effects.
#' @param pos.gx vector of integer indicating the positions of the `X = [x1, x2, ...]` columns corresponding to observed `GX`. If `formula` is given as `y ~ x1 + x2 || gx2`,
#' then the corresponding position is 2 because `gx2` corresponds to the second variable in `X`.  If  `formula` is given as `y ~ x1 + x2 || gx2 + gx1`, then `pos.gx` should be `c(2, 1)`,
#' because the first observed `gx` corresponds to the second variable in `X` and the second observed `gx` correspond to the first variable in `X`. If `pos.gx` is missing, then 
#' the program assumes that the observed contextual variables are in the same order than that of `X`.
#' @param fixed.effects logical; if true, group heterogeneity is included as fixed effects.
#' @param dnetwork a list, where the m-th elements is the matrix of link probability in the m-th sub-network. 
#' @param R the number of simulation `R`.
#' @param opt.ctr a list a controls to paste into \link[optimize]{stat}. The peer effect is the solution of a nonlinear equation having zero on the right-hand side. 
#' In practice, the program minimizes the square of the left-hand side using \link[optimize]{stat} and verifies that the minimum is zero.#' In practice, the program minimizes the square of the left-hand-side using \link[optimize]{stat} and verifies that the minimum is zero.
#' @param cond.var logical; if true the estimator variance conditional on `dnetwork` will be computed.
#' @param `print.proc` (logical) indicates if the optimization process should be printed step by step.
#' @param data optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If missing, the variables are taken from \code{environment(formula)}, typically the environment from which `qmleSAR` is called.
#' @return A list consisting of:
#'     \item{n.group}{number of groups.}
#'     \item{N}{vector of each group size.}
#'     \item{time}{elapsed time to run the qml estimator.}
#'     \item{estimates}{vector of estimated parameters.}
#'     \item{formula}{input value of `formula`.}
#'     \item{contextual}{input value of `contextual`.}
#'     \item{fixed.effects}{input value of `fixed.effects`.}
#'     \item{details}{other details of the model.}
#' @importFrom stats optimize
#' @examples 
#' \donttest{
#' # Number of groups
#' M        <- 100
#' # size of each group
#' N        <- rep(30,M)
#' # covariates
#' X        <- cbind(rnorm(sum(N),0,5),rpois(sum(N),7))
#' # network formation model parameter
#' rho      <- c(-0.8, 0.2, -0.1)
#' # individual effects
#' beta     <- c(2, 1, 1.5, 5, -3)
#' # endogenous effects
#' alpha    <- 0.4
#' # std-dev errors
#' se       <- 1
#' # network
#' tmp      <- c(0, cumsum(N))
#' X1l      <- lapply(1:M, function(x) X[c(tmp[x] + 1):tmp[x+1],1])
#' X2l      <- lapply(1:M, function(x) X[c(tmp[x] + 1):tmp[x+1],2])
#' dist.net <- function(x, y) abs(x - y)
#' X1.mat   <- lapply(1:M, function(m) {
#'   matrix(kronecker(X1l[[m]], X1l[[m]], FUN = dist.net), N[m])})
#' X2.mat   <- lapply(1:M, function(m) {
#'   matrix(kronecker(X2l[[m]], X2l[[m]], FUN = dist.net), N[m])})
#' Xnet     <- as.matrix(cbind("Const" = 1,
#'                             "dX1"   = mat.to.vec(X1.mat),
#'                             "dX2"   = mat.to.vec(X2.mat)))
#' ynet     <- Xnet %*% rho
#' ynet     <- c(1*((ynet + rlogis(length(ynet))) > 0))
#' G0       <- vec.to.mat(ynet, N, normalise = FALSE)
#' # normalise
#' G0norm   <- norm.network(G0)
#' # Matrix GX
#' GX       <- peer.avg(G0norm, X)
#' # simulate dependent variable use an external package
#' y        <- CDatanet::simsar(~ X, contextual = TRUE, Glist = G0norm,
#'                              theta = c(alpha, beta, se))
#' Gy       <- y$Gy
#' y        <- y$y
#' # build dataset
#' dataset           <- as.data.frame(cbind(y, X, Gy, GX))
#' colnames(dataset) <- c("y","X1","X2", "Gy", "GX1", "GX2")
#' nNet      <- nrow(Xnet) # network formation model sample size
#' Aobs      <- sample(1:nNet, round(0.3*nNet)) # We observed 30%
#' # We can estimate rho using the gml function from the stats package
#' logestim  <- glm(ynet[Aobs] ~ -1 + Xnet[Aobs,], family = binomial(link = "logit"))
#' slogestim <- summary(logestim)
#' rho.est   <- logestim$coefficients
#' rho.var   <- slogestim$cov.unscaled # we also need the covariance of the estimator
#' 
#' d.logit     <- lapply(1:M, function(x) {
#'   out       <- 1/(1 + exp(-rho.est[1] - rho.est[2]*X1.mat[[x]] -
#'                             rho.est[3]*X2.mat[[x]]))
#'   diag(out) <- 0
#'   out})
#' qmle.logit   <- qmleSAR(y ~ X1 + X2, dnetwork = d.logit, contextual = TRUE,
#'                       qmle.ctr  = list(R = 100L, print = TRUE), data = dataset)
#' summary(qmle.logit, dnetwork = d.logit, data = dataset)
#' }
#' @export
qmleSAR <- function(formula,
                   contextual    = FALSE,
                   pos.gx        = NULL,
                   fixed.effects = FALSE,
                   dnetwork,
                   R             = 5L,
                   opt.ctr       = list(interval = c(-1, 1), tol = .Machine$double.eps^0.25),
                   print.proc    = TRUE,
                   data){
  
  t1           <- Sys.time()
  if(!exists(".Random.seed")) set.seed(0)
  seed         <- .Random.seed 
  S            <- 1L
  
  # data
  env.f        <- environment(formula)
  f.t.data     <- formula.to.data.smm(formula = formula, data = data, fixed.effects = fixed.effects) 
  formula      <- f.t.data$formula; environment(formula) <- env.f
  X            <- f.t.data$X
  y            <- f.t.data$y
  GX1          <- f.t.data$GX
  Gy           <- f.t.data$Gy
  Gyobs        <- !is.null(Gy)
  GX1obs       <- !is.null(GX1)
  col.x        <- colnames(X)
  col.gx1      <- colnames(GX1)
  intercept    <- ("(Intercept)" %in% col.x)
  Kx           <- length(col.x) 
  Kx1          <- length(col.gx1) 
  Kx2          <- Kx - intercept - Kx1
  if(Kx1 > 0 & is.null(pos.gx)){
    pos.gx     <- 1:Kx1
  }
  if(Kx1 != length(pos.gx)) stop("length(pos.gx) not equal to the number of observed gx")
  if(Kx1 > 0){
    if((max(pos.gx) + intercept) > Kx) stop("max(pos.gx) greater than the number of columns in X")
  }
  if(GX1obs & !contextual){
    warning("contextual is set fixed but contextual variables are observed")
  }
  X1           <- X[,pos.gx + intercept, drop = FALSE]
  X2           <- NULL
  if(is.null(pos.gx)){
    X2         <- X[,((1 + intercept):Kx), drop = FALSE]
  } else {
    X2         <- X[,((1 + intercept):Kx)[-pos.gx], drop = FALSE]
  }
  if(Kx1 == 0){
    GX1        <- X1
  }
  
  #sizes
  N            <- sapply(dnetwork, nrow)
  M            <- length(N)
  Nsum         <- sum(N)
  Ncum         <- c(0, cumsum(N))
  Ilist        <- lapply(N, diag)
  
  # type
  part         <- NULL
  if((Gyobs & !contextual) | (Gyobs & contextual & Kx2 == 0)){
    part       <- 0
  }
  if(Gyobs & contextual & Kx2 != 0){
    part       <- 2
  }
  if(!Gyobs){
    part       <- 1
  }
  model        <- paste0(part, "_", as.numeric(contextual))

  # Others
  ctr.optim    <- list(R = R, distr = dnetwork, Ilist = Ilist, y = y, Gy = Gy, X = X, X1 = X1, GX1 = GX1, 
                     X2 = X2, Kx = Kx, Kx1 = Kx1, Kx2 = Kx2, M = M, Ncum = Ncum, FE = fixed.effects, 
                     seed = seed, print.proc = print.proc, S = S)
  
  # optimization
  arg.optim    <- c(list(f = get(paste0("qmleopt", model))), opt.ctr, ctr.optim)
  out.optim    <- do.call(what = optimize, args = arg.optim) 
  alpha        <- out.optim$minimum 
  foc          <- out.optim$objective
  
  # compute beta and sigma
  bse          <- do.call(what = get(paste0("qmlebeta", model)), args = c(list(alpha = alpha), ctr.optim))
  beta         <- c(bse$beta)
  se2          <- bse$se2
  
  col.gx       <- NULL
  if(contextual) {
    if(Kx2 == 0){
      col.gx   <- col.gx1
    } else {
      col.gx   <- c(col.gx1, paste0("G: ", colnames(X2)))
    }
  }
  
  theta        <- c(beta, alpha, se2)
  names(theta) <- c(col.x, col.gx, "peer effects", "se2")
  
  # time
  t2           <- Sys.time()
  timer        <- as.numeric(difftime(t2, t1, units = "secs")) 
  nhours       <- floor(timer/3600)
  nminutes     <- floor((timer-3600*nhours)/60)%%60
  nseconds     <- timer-3600*nhours-60*nminutes
  
  # details
  infos        <- list(optimize     = out.optim,
                       Gyobs        = Gyobs, 
                       Kx           = Kx,
                       Kx1          = Kx1,
                       Kx2          = Kx2,
                       time         = c("Elapsed times (seconds)" = timer),
                       seed         = seed)
  
  
  if (print.proc > 0) {
    cat("Elapsed time: ", nhours, " HH ", nminutes, " mm ", round(nseconds), " ss \n \n")
  }
  
  out          <- list(n.group       = M,
                       N             = N,
                       estimate      = theta,
                       objective     = foc,
                       formula       = formula,
                       contextual    = contextual,
                       fixed.effects = fixed.effects,
                       details       = infos)
  class(out)   <-"qmleSAR"
  out
}

#' @title Summarizing the PML Estimator of SAR model with Partial Network Data
#' @description Summary and print methods for the class `qmlwSAR`.
#' @param object an object of class "qmleSAR", output of the function \code{\link{qmleSAR}}.
#' @param x an object of class "summary.smmSAR" or "smmSAR", output of the functions \code{\link{summary.qmleSAR}} or
#' \code{\link{qmleSAR}}.
#' @param .fun,.args are used to simulate from the distribution of `dnetwork`. `.fun` is the simulator function
#' where `.args` is a list of its arguments. Typically `do.call(.fun, .args)` is supposed to simulate one `dnetwork` from
#' the distribution.
#' @param sim the number of simulations of `dnetwork`.
#' @param ncores the number of cores to be used for the simulation. Use a lot of cores for fast simulations.
#' @param dnetwork a list, where the m-th elements is the matrix of link probability in the m-th sub-network. 
#' @param data optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If missing, the variables are taken from \code{environment(formula)}, typically the environment from which `smmSAR` is called.
#' @param ... further arguments passed to or from other methods.
#' @return A list consisting of:
#'     \item{n.group}{number of groups.}
#'     \item{N}{vector of each group size.}
#'     \item{estimates}{vector of estimated parameters.}
#'     \item{formula}{input value of `formula`.}
#'     \item{contextual}{input value of `contextual`.}
#'     \item{fixed.effects}{input value of `fixed.effects`.}
#'     \item{details}{other details of the model.}
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach "%dopar%"
#' @importFrom doRNG "%dorng%"
#' @importFrom stats var
#' @export
"summary.qmleSAR" <- function(object, 
                             .fun, 
                             .args, 
                             sim    = 30,
                             ncores = 1,
                             dnetwork, 
                             data,
                             ...){
  # stopifnot(inherits(object, "smmSAR"))
  # details       <- object$details
  # derM          <- details$av.grad.m
  # aveMM         <- details$`av.m%*%t(m)`
  # aveM          <- details$av.m
  # 
  # 
  # out       <- list(n.group       = object$n.group,
  #                   N             = object$N,
  #                   estimates     = object$estimates,
  #                   cov           = varcov,
  #                   formula       = object$formula,
  #                   contextual    = object$contextual,
  #                   fixed.effects = object$fixed.effects,
  #                   smm.ctr       = object$smm.ctr,
  #                   details       = details,
  #                   ...           = ...)
  # class(out) <- "summary.qmleSAR"
  # out
}


fSIGMA <- function(.fun, .args, fmvzeta, Afmvzeta, M) {
  Afmvzeta$distr  <- do.call(.fun, .args)
  tmp   <- do.call(fmvzeta, Afmvzeta)
  sumMM <- tmp$sumMM
  sumM  <- tmp$sumM
  list(VZ = sumMM - sumM %*% t(sumM)/M, EZ = sumM)
}

#' @rdname summary.qmleSAR
#' @export
"print.summary.qmleSAR"  <- function(x, ...) {
  stopifnot(inherits(x, "summary.qmleSAR"))
  
  M          <- x$n.group
  N          <- x$N
  estimates  <- x$estimates
  sumN       <- sum(N)
  tmp        <- fcoefficients(x$estimates, sqrt(diag(x$cov)))
  out_print  <- tmp$out_print
  out        <- tmp$out
  out_print  <- c(list(out_print), x[-(1:9)], list(...))
  
  nR         <- x$smm.ctr$R
  nS         <- 1#x$smm.ctr$S
  nT         <- 1#x$smm.ctr$T
  
  cat("Simulated Method of Moments estimation of SAR model", "\n\n")
  cat("Formula = ", Reduce(paste, deparse(x$formula)), "\n\n", sep = "")
  cat("Contextual effects: ", ifelse(x$contextual, "Yes", "No"), "\n", sep = "")
  cat("Fixed effects: ", ifelse(x$fixed.effects, "Yes", "No"), "\n\n", sep = "")
  cat("Network details\n")
  cat("GX ", ifelse(x$details$GXobs, "Observed", "Not Observed"), "\n", sep = "")
  cat("Gy ", ifelse(x$details$Gyobs, "Observed", "Not Observed"), "\n", sep = "")
  cat("Number of groups: ", M, "\n", sep = "")
  cat("Sample size     : ", sum(N), "\n\n", sep = "")

  
  #cat("Simulation tuning parameters:\n")
  #cat("R = ", nR, ifelse(!is.null(nS), paste0(", S = ", nS), ""), ifelse(!is.null(nT), paste0(", T = ", nT), "")  , "\n\n", sep = "")
  cat("Simulation settings\n")
  cat("R = ", nR, "\n", sep = "")
  cat("Smoother : ", x$details$smoother, sep = "")
  if(x$details$smoother)cat(" (h = ", x$details$h, ")", sep = "")
  cat("\n\n")
  
  cat("Coefficients:\n")
  do.call("print", out_print)
  cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  
  out        <- list(n.group       = x$n.group,
                     N             = x$N,
                     estimates     = x$estimates,
                     cov           = x$cov,
                     coefficients  = out,
                     formula       = x$formula,
                     contextual    = x$contextual,
                     fixed.effects = x$fixed.effects,
                     smm.ctr       = x$smm.ctr)
  class(out) <- "print.summary.qmleSAR"
  invisible(out)
}

#' @rdname summary.qmleSAR
#' @export
"print.qmleSAR"  <- function(x,
                            dnetwork, 
                            .fun, 
                            .args, 
                            sim    = NULL,
                            ncores = 1,
                            data, 
                            ...) {
  stopifnot(inherits(x, "qmleSAR"))
  print(summary(object   = x,
                dnetwork = dnetwork, 
                .fun     = .fun, 
                .args    = .args, 
                sim      = sim,
                ncores   = ncores,
                data     = data, ...))
}