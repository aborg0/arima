package com.mind_era.arima

/*
from: https://svn.r-project.org/R/trunk/src/library/stats/R/arima.R
arima <- function(x, order = c(0L, 0L, 0L),
                  seasonal = list(order = c(0L, 0L, 0L), period = NA),
                  xreg = NULL, include.mean = TRUE,
                  transform.pars = TRUE, fixed = NULL, init = NULL,
                  method = c("CSS-ML", "ML", "CSS"), n.cond,
                  SSinit = c("Gardner1980", "Rossignol2011"),
                  optim.method = "BFGS",
                  optim.control = list(), kappa = 1e6)
{
    "%+%" <- function(a, b) .Call(C_TSconv, a, b)

    SSinit <- match.arg(SSinit)
    SS.G <- SSinit == "Gardner1980"
    ## helper of armafn(), called by optim()
    upARIMA <- function(mod, phi, theta)
    {
        p <- length(phi); q <- length(theta)
        mod$phi <- phi; mod$theta <- theta
        r <- max(p, q + 1L)
        if(p > 0) mod$T[1L:p, 1L] <- phi
	if(r > 1L)
	    mod$Pn[1L:r, 1L:r] <-
		if(SS.G) .Call(C_getQ0, phi, theta)
		else .Call(C_getQ0bis, phi, theta, tol = 0)# tol=0: less checking
	else
	    mod$Pn[1L, 1L] <- if (p > 0) 1/(1 - phi^2) else 1
        mod$a[] <- 0
        mod
    }

    arimaSS <- function(y, mod)
    {
        ## next call changes mod components a, P, Pn so beware!
        .Call(C_ARIMA_Like, y, mod, 0L, TRUE)
    }

    ## the objective function called by optim()
    armafn <- function(p, trans)
    {
        par <- coef
        par[mask] <- p
        trarma <- .Call(C_ARIMA_transPars, par, arma, trans)
	if(is.null(Z <- tryCatch(upARIMA(mod, trarma[[1L]], trarma[[2L]]),
				 error = function(e) NULL)))
	    return(.Machine$double.xmax)# bad parameters giving error, e.g. in solve(.)
        if(ncxreg > 0) x <- x - xreg %*% par[narma + (1L:ncxreg)]
        ## next call changes Z components a, P, Pn so beware!
        res <- .Call(C_ARIMA_Like, x, Z, 0L, FALSE)
        s2 <- res[1L]/res[3L]
        0.5*(log(s2) + res[2L]/res[3L])
    }

    armaCSS <- function(p)
    {
        par <- as.double(fixed)
        par[mask] <- p
        trarma <- .Call(C_ARIMA_transPars, par, arma, FALSE)
        if(ncxreg > 0) x <- x - xreg %*% par[narma + (1L:ncxreg)]
        res <- .Call(C_ARIMA_CSS, x, arma, trarma[[1L]], trarma[[2L]],
                     as.integer(ncond), FALSE)
        0.5 * log(res)
    }

    arCheck <- function(ar)
    {
        p <- max(which(c(1, -ar) != 0)) - 1
        if(!p) return(TRUE)
        all(Mod(polyroot(c(1, -ar[1L:p]))) > 1)
    }

    maInvert <- function(ma)
    {
        ## polyroot can't cope with leading zero.
        q <- length(ma)
        q0 <- max(which(c(1,ma) != 0)) - 1L
        if(!q0) return(ma)
        roots <- polyroot(c(1, ma[1L:q0]))
        ind <- Mod(roots) < 1
        if(all(!ind)) return(ma)
        if(q0 == 1) return(c(1/ma[1L], rep.int(0, q - q0)))
        roots[ind] <- 1/roots[ind]
        x <- 1
        for (r in roots) x <- c(x, 0) - c(0, x)/r
        c(Re(x[-1L]), rep.int(0, q - q0))
    }

    series <- deparse(substitute(x))
    if(NCOL(x) > 1L)
        stop("only implemented for univariate time series")
    method <- match.arg(method)

    x <- as.ts(x)
    if(!is.numeric(x))
        stop("'x' must be numeric")
    storage.mode(x) <- "double"  # a precaution
    dim(x) <- NULL
    n <- length(x)

    if(!missing(order))
        if(!is.numeric(order) || length(order) != 3L || any(order < 0))
            stop("'order' must be a non-negative numeric vector of length 3")
    if(!missing(seasonal))
        if(is.list(seasonal)) {
            if(is.null(seasonal$order))
                stop("'seasonal' must be a list with component 'order'")
            if(!is.numeric(seasonal$order) || length(seasonal$order) != 3L
               || any(seasonal$order < 0L))
                stop("'seasonal$order' must be a non-negative numeric vector of length 3")
        } else if(is.numeric(order)) {
            if(length(order) == 3L) seasonal <- list(order=seasonal)
            else ("'seasonal' is of the wrong length")
        } else stop("'seasonal' must be a list with component 'order'")

    if (is.null(seasonal$period) || is.na(seasonal$period)
        ||seasonal$period == 0) seasonal$period <- frequency(x)
    arma <- as.integer(c(order[-2L], seasonal$order[-2L], seasonal$period,
                         order[2L], seasonal$order[2L]))
    narma <- sum(arma[1L:4L])

    xtsp <- tsp(x)
    tsp(x) <- NULL
    Delta <- 1.
    for(i in seq_len(order[2L])) Delta <- Delta %+% c(1., -1.)
    for(i in seq_len(seasonal$order[2L]))
        Delta <- Delta %+% c(1, rep.int(0, seasonal$period-1), -1)
    Delta <- - Delta[-1L]
    nd <- order[2L] + seasonal$order[2L]
    n.used <- sum(!is.na(x)) - length(Delta)
    if (is.null(xreg)) {
        ncxreg <- 0L
    } else {
        nmxreg <- deparse(substitute(xreg))
        if (NROW(xreg) != n) stop("lengths of 'x' and 'xreg' do not match")
        ncxreg <- NCOL(xreg)
        xreg <- as.matrix(xreg)
        storage.mode(xreg) <- "double"
    }
    class(xreg) <- NULL
    if (ncxreg > 0L && is.null(colnames(xreg)))
        colnames(xreg) <-
            if(ncxreg == 1L) nmxreg else paste0(nmxreg, 1L:ncxreg)
    if (include.mean && (nd == 0L)) {
        xreg <- cbind(intercept = rep(1, n), xreg = xreg)
        ncxreg <- ncxreg + 1L
    }
    if(method == "CSS-ML") {
        anyna <- anyNA(x)
        if(ncxreg) anyna <- anyna || anyNA(xreg)
        if(anyna) method <- "ML"
    }

    if (method == "CSS" || method == "CSS-ML") {
        ncond <- order[2L] + seasonal$order[2L] * seasonal$period
        ncond1 <- order[1L] + seasonal$period * seasonal$order[1L]
	ncond <- ncond + if(!missing(n.cond)) max(n.cond, ncond1) else ncond1
    } else ncond <- 0

    if (is.null(fixed)) fixed <- rep(NA_real_, narma + ncxreg)
    else if(length(fixed) != narma + ncxreg) stop("wrong length for 'fixed'")
    mask <- is.na(fixed)
##    if(!any(mask)) stop("all parameters were fixed")
    no.optim <- !any(mask)
    if(no.optim) transform.pars <- FALSE
    if(transform.pars) {
        ind <- arma[1L] + arma[2L] + seq_len(arma[3L])
        if (any(!mask[seq_len(arma[1L])]) || any(!mask[ind])) {
            warning("some AR parameters were fixed: setting transform.pars = FALSE")
            transform.pars <- FALSE
        }
    }
    init0 <- rep.int(0, narma)
    parscale <- rep(1, narma)
    if (ncxreg) {
        cn <- colnames(xreg)
        orig.xreg <- (ncxreg == 1L) || any(!mask[narma + 1L:ncxreg])
        if (!orig.xreg) {
            S <- svd(na.omit(xreg))
            xreg <- xreg %*% S$v
        }
        dx <- x
        dxreg <- xreg
        if(order[2L] > 0L) {
            dx <- diff(dx, 1L, order[2L])
            dxreg <- diff(dxreg, 1L, order[2L])
        }
        if(seasonal$period > 1L & seasonal$order[2L] > 0) {
            dx <- diff(dx, seasonal$period, seasonal$order[2L])
            dxreg <- diff(dxreg, seasonal$period, seasonal$order[2L])
        }
        fit <- if(length(dx) > ncol(dxreg))
            lm(dx ~ dxreg - 1, na.action = na.omit)
        else list(rank = 0L)
        if(fit$rank == 0L) {
            ## Degenerate model. Proceed anyway so as not to break old code
            fit <- lm(x ~ xreg - 1, na.action = na.omit)
        }
        isna <- is.na(x) | apply(xreg, 1L, anyNA)
        n.used <- sum(!isna) - length(Delta)
        init0 <- c(init0, coef(fit))
        ses <- summary(fit)$coefficients[, 2L]
        parscale <- c(parscale, 10 * ses)
    }
    if (n.used <= 0) stop("too few non-missing observations")

    if(!is.null(init)) {
        if(length(init) != length(init0))
            stop("'init' is of the wrong length")
        if(any(ind <- is.na(init))) init[ind] <- init0[ind]
        if(method == "ML") {
            ## check stationarity
            if(arma[1L] > 0)
                if(!arCheck(init[1L:arma[1L]]))
                    stop("non-stationary AR part")
            if(arma[3L] > 0)
                if(!arCheck(init[sum(arma[1L:2L]) + 1L:arma[3L]]))
                    stop("non-stationary seasonal AR part")
            if(transform.pars)
                init <- .Call(C_ARIMA_Invtrans, as.double(init), arma)
        }
    } else init <- init0

    coef <- as.double(fixed)
    if(!("parscale" %in% names(optim.control)))
       optim.control$parscale <- parscale[mask]

    if(method == "CSS") {
        res <- if(no.optim)
            list(convergence=0L, par=numeric(), value=armaCSS(numeric()))
        else
            optim(init[mask], armaCSS,  method = optim.method, hessian = TRUE,
                  control = optim.control)
        if(res$convergence > 0)
            warning(gettextf("possible convergence problem: optim gave code = %d",
                             res$convergence), domain = NA)
        coef[mask] <- res$par
        ## set model for predictions
        trarma <- .Call(C_ARIMA_transPars, coef, arma, FALSE)
	mod <- makeARIMA(trarma[[1L]], trarma[[2L]], Delta, kappa, SSinit)
        if(ncxreg > 0) x <- x - xreg %*% coef[narma + (1L:ncxreg)]
        arimaSS(x, mod)
        val <- .Call(C_ARIMA_CSS, x, arma, trarma[[1L]], trarma[[2L]],
                     as.integer(ncond), TRUE)
        sigma2 <- val[[1L]]
        var <- if(no.optim) numeric() else solve(res$hessian * n.used)
    } else {
        if(method == "CSS-ML") {
            res <- if(no.optim)
                list(convergence=0L, par=numeric(), value=armaCSS(numeric()))
            else
                optim(init[mask], armaCSS,  method = optim.method,
                      hessian = FALSE, control = optim.control)
            if(res$convergence == 0) init[mask] <- res$par
            ## check stationarity
            if(arma[1L] > 0)
                if(!arCheck(init[1L:arma[1L]]))
                    stop("non-stationary AR part from CSS")
            if(arma[3L] > 0)
                if(!arCheck(init[sum(arma[1L:2L]) + 1L:arma[3L]]))
                    stop("non-stationary seasonal AR part from CSS")
            ncond <- 0L
        }
        if(transform.pars) {
            init <- .Call(C_ARIMA_Invtrans, init, arma)
            ## enforce invertibility
            if(arma[2L] > 0) {
                ind <- arma[1L] + 1L:arma[2L]
                init[ind] <- maInvert(init[ind])
            }
            if(arma[4L] > 0) {
                ind <- sum(arma[1L:3L]) + 1L:arma[4L]
                init[ind] <- maInvert(init[ind])
            }
        }
        trarma <- .Call(C_ARIMA_transPars, init, arma, transform.pars)
	mod <- makeARIMA(trarma[[1L]], trarma[[2L]], Delta, kappa, SSinit)
        res <- if(no.optim)
            list(convergence = 0, par = numeric(),
                 value = armafn(numeric(), as.logical(transform.pars)))
        else
            optim(init[mask], armafn, method = optim.method,
                  hessian = TRUE, control = optim.control,
                  trans = as.logical(transform.pars))
        if(res$convergence > 0)
            warning(gettextf("possible convergence problem: optim gave code = %d",
                             res$convergence), domain = NA)
        coef[mask] <- res$par
        if(transform.pars) {
            ## enforce invertibility
            if(arma[2L] > 0L) {
                ind <- arma[1L] + 1L:arma[2L]
                if(all(mask[ind]))
                    coef[ind] <- maInvert(coef[ind])
            }
            if(arma[4L] > 0L) {
                ind <- sum(arma[1L:3L]) + 1L:arma[4L]
                if(all(mask[ind]))
                    coef[ind] <- maInvert(coef[ind])
            }
            if(any(coef[mask] != res$par))  {  # need to re-fit
                oldcode <- res$convergence
                res <- optim(coef[mask], armafn, method = optim.method,
                             hessian = TRUE,
                             control = list(maxit = 0L,
                             parscale = optim.control$parscale),
                             trans = TRUE)
                res$convergence <- oldcode
                coef[mask] <- res$par
            }
            ## do it this way to ensure hessian was computed inside
            ## stationarity region
            A <- .Call(C_ARIMA_Gradtrans, as.double(coef), arma)
            A <- A[mask, mask]
	    var <- crossprod(A, solve(res$hessian * n.used, A))
            coef <- .Call(C_ARIMA_undoPars, coef, arma)
        } else var <- if(no.optim) numeric() else solve(res$hessian * n.used)
        trarma <- .Call(C_ARIMA_transPars, coef, arma, FALSE)
	mod <- makeARIMA(trarma[[1L]], trarma[[2L]], Delta, kappa, SSinit)
        val <- if(ncxreg > 0L)
            arimaSS(x - xreg %*% coef[narma + (1L:ncxreg)], mod)
        else arimaSS(x, mod)
        sigma2 <- val[[1L]][1L]/n.used
    }
    value <- 2 * n.used * res$value + n.used + n.used * log(2 * pi)
    aic <- if(method != "CSS") value + 2*sum(mask) + 2 else NA
    nm <- NULL
    if (arma[1L] > 0L) nm <- c(nm, paste0("ar", 1L:arma[1L]))
    if (arma[2L] > 0L) nm <- c(nm, paste0("ma", 1L:arma[2L]))
    if (arma[3L] > 0L) nm <- c(nm, paste0("sar", 1L:arma[3L]))
    if (arma[4L] > 0L) nm <- c(nm, paste0("sma", 1L:arma[4L]))
    if (ncxreg > 0L) {
        nm <- c(nm, cn)
        if(!orig.xreg) {
            ind <- narma + 1L:ncxreg
            coef[ind] <- S$v %*% coef[ind]
            A <- diag(narma + ncxreg)
            A[ind, ind] <- S$v
            A <- A[mask, mask]
            var <- A %*% var %*% t(A)
        }
    }
    names(coef) <- nm
    if(!no.optim) dimnames(var) <- list(nm[mask], nm[mask])
    resid <- val[[2L]]
    tsp(resid) <- xtsp
    class(resid) <- "ts"
    structure(list(coef = coef, sigma2 = sigma2, var.coef = var, mask = mask,
		   loglik = -0.5 * value, aic = aic, arma = arma,
		   residuals = resid, call = match.call(), series = series,
		   code = res$convergence, n.cond = ncond, nobs = n.used,
		   model = mod),
	      class = "Arima")
}


print.Arima <-
    function (x, digits = max(3L, getOption("digits") - 3L), se = TRUE, ...)
{
    cat("\nCall:", deparse(x$call, width.cutoff = 75L), "", sep = "\n")
    if (length(x$coef)) {
        cat("Coefficients:\n")
        coef <- round(x$coef, digits = digits)
        ## use NROW as if all coefs are fixed there are no var.coef's
        if (se && NROW(x$var.coef)) {
            ses <- rep.int(0, length(coef))
            ses[x$mask] <- round(sqrt(diag(x$var.coef)), digits = digits)
            coef <- matrix(coef, 1L, dimnames = list(NULL, names(coef)))
            coef <- rbind(coef, s.e. = ses)
        }
        print.default(coef, print.gap = 2)
    }
    cm <- x$call$method
    if(is.null(cm) || cm != "CSS")
        cat("\nsigma^2 estimated as ", format(x$sigma2, digits = digits),
            ":  log likelihood = ", format(round(x$loglik, 2L)),
            ",  aic = ", format(round(x$aic, 2L)), "\n", sep = "")
    else
        cat("\nsigma^2 estimated as ",
            format(x$sigma2, digits = digits),
            ":  part log likelihood = ", format(round(x$loglik,2)),
            "\n", sep = "")
    invisible(x)
}


predict.Arima <-
    function (object, n.ahead = 1L, newxreg = NULL, se.fit = TRUE, ...)
{
    myNCOL <- function(x) if (is.null(x)) 0 else NCOL(x)
    rsd <- object$residuals
    xr <- object$call$xreg
    xreg <- if (!is.null(xr)) eval.parent(xr) else NULL
    ncxreg <- myNCOL(xreg)
    if (myNCOL(newxreg) != ncxreg)
        stop("'xreg' and 'newxreg' have different numbers of columns")
    class(xreg) <- NULL
    xtsp <- tsp(rsd)
    n <- length(rsd)
    arma <- object$arma
    coefs <- object$coef
    narma <- sum(arma[1L:4L])
    if (length(coefs) > narma) {
        if (names(coefs)[narma + 1L] == "intercept") {
            xreg <- cbind(intercept = rep(1, n), xreg)
            newxreg <- cbind(intercept = rep(1, n.ahead), newxreg)
            ncxreg <- ncxreg + 1L
        }
        xm <- if(narma == 0) drop(as.matrix(newxreg) %*% coefs)
        else drop(as.matrix(newxreg) %*% coefs[-(1L:narma)])
    }
    else xm <- 0
    if (arma[2L] > 0L) {
        ma <- coefs[arma[1L] + 1L:arma[2L]]
        if (any(Mod(polyroot(c(1, ma))) < 1))
            warning("MA part of model is not invertible")
    }
    if (arma[4L] > 0L) {
        ma <- coefs[sum(arma[1L:3L]) + 1L:arma[4L]]
        if (any(Mod(polyroot(c(1, ma))) < 1))
            warning("seasonal MA part of model is not invertible")
    }
    z <- KalmanForecast(n.ahead, object$model)
    pred <- ts(z[[1L]] + xm, start = xtsp[2L] + deltat(rsd),
               frequency = xtsp[3L])
    if (se.fit) {
        se <- ts(sqrt(z[[2L]] * object$sigma2),
                 start = xtsp[2L] + deltat(rsd),
                 frequency = xtsp[3L])
        list(pred=pred, se=se)
    }
    else pred
}


makeARIMA <- function(phi, theta, Delta, kappa = 1e6,
		      SSinit = c("Gardner1980", "Rossignol2011"),
		      tol = .Machine$double.eps)
{
    if(anyNA(phi))   warning(gettextf("NAs in '%s'", "phi"), domain=NA)
    if(anyNA(theta)) warning(gettextf("NAs in '%s'", "theta"), domain=NA)
    p <- length(phi); q <- length(theta)
    r <- max(p, q + 1L); d <- length(Delta)
    rd <- r + d
    Z <- c(1., rep.int(0, r-1L), Delta)
    T <- matrix(0., rd, rd)
    if(p > 0) T[1L:p, 1L] <- phi
    if(r > 1L) {
        ind <- 2:r
        T[cbind(ind-1L, ind)] <- 1
    }
    if(d > 0L) {
        T[r+1L, ] <- Z
        if(d > 1L) {
            ind <- r + 2:d
            T[cbind(ind, ind-1)] <- 1
        }
    }
    if(q < r - 1L) theta <- c(theta, rep.int(0, r-1L-q))
    R <- c(1, theta, rep.int(0, d))
    V <- R %o% R
    h <- 0.
    a <- rep(0., rd)
    Pn <- P <- matrix(0., rd, rd)
    if(r > 1L)
        Pn[1L:r, 1L:r] <- switch(match.arg(SSinit),
                                 "Gardner1980" = .Call(C_getQ0, phi, theta),
                                 "Rossignol2011" = .Call(C_getQ0bis, phi, theta, tol),
                                 stop("invalid 'SSinit'"))
    else Pn[1L, 1L] <- if(p > 0) 1/(1 - phi^2) else 1
    if(d > 0L) Pn[cbind(r+1L:d, r+1L:d)] <- kappa
    list(phi=phi, theta=theta, Delta=Delta, Z=Z, a=a, P=P, T=T, V=V,
         h=h, Pn=Pn)
}

coef.Arima <- function (object, ...) object$coef

vcov.Arima <- function (object, ...) object$var.coef

logLik.Arima <- function (object, ...) {
    res <- if(is.na(object$aic)) NA
    else structure(object$loglik, df = sum(object$mask) + 1, nobs = object$nobs)
    class(res) <- "logLik"
    res
}
https://svn.r-project.org/R/trunk/src/library/stats/src/arima.c

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdlib.h> // for abs
#include <string.h>

#include <R.h>
#include "ts.h"
#include "statsR.h" // for getListElement

#ifndef max
#define max(a,b) ((a < b)?(b):(a))
#endif
#ifndef min
#define min(a,b) ((a < b)?(a):(b))
#endif


/*
  KalmanLike, internal to StructTS:
  .Call(C_KalmanLike, y, mod$Z, mod$a, mod$P, mod$T, mod$V, mod$h, mod$Pn,
        nit, FALSE, update)
  KalmanRun:
  .Call(C_KalmanLike, y, mod$Z, mod$a, mod$P, mod$T, mod$V, mod$h, mod$Pn,
        nit, TRUE, update)
*/

/* y vector length n of observations
   Z vector length p for observation equation y_t = Za_t +  eps_t
   a vector length p of initial state
   P p x p matrix for initial state uncertainty (contemparaneous)
   T  p x p transition matrix
   V  p x p = RQR'
   h = var(eps_t)
   anew used for a[t|t-1]
   Pnew used for P[t|t -1]
   M used for M = P[t|t -1]Z

   op is FALSE for KalmanLike, TRUE for KalmanRun.
   The latter computes residuals and states and has
   a more elaborate return value.

   Almost no checking here!
 */

SEXP
KalmanLike(SEXP sy, SEXP mod, SEXP sUP, SEXP op, SEXP update)
{
    int lop = asLogical(op);
    mod = PROTECT(duplicate(mod));

    SEXP sZ = getListElement(mod, "Z"), sa = getListElement(mod, "a"),
	sP = getListElement(mod, "P"), sT = getListElement(mod, "T"),
	sV = getListElement(mod, "V"), sh = getListElement(mod, "h"),
	sPn = getListElement(mod, "Pn");

    if (TYPEOF(sy) != REALSXP || TYPEOF(sZ) != REALSXP ||
	TYPEOF(sa) != REALSXP || TYPEOF(sP) != REALSXP ||
	TYPEOF(sPn) != REALSXP ||
	TYPEOF(sT) != REALSXP || TYPEOF(sV) != REALSXP)
	error(_("invalid argument type"));

    int n = LENGTH(sy), p = LENGTH(sa);
    double *y = REAL(sy), *Z = REAL(sZ), *T = REAL(sT), *V = REAL(sV),
	*P = REAL(sP), *a = REAL(sa), *Pnew = REAL(sPn), h = asReal(sh);

    double *anew = (double *) R_alloc(p, sizeof(double));
    double *M = (double *) R_alloc(p, sizeof(double));
    double *mm = (double *) R_alloc(p * p, sizeof(double));
    // These are only used if(lop), but avoid -Wall trouble
    SEXP ans = R_NilValue, resid = R_NilValue, states = R_NilValue;
    if(lop) {
	PROTECT(ans = allocVector(VECSXP, 3));
	SET_VECTOR_ELT(ans, 1, resid = allocVector(REALSXP, n));
	SET_VECTOR_ELT(ans, 2, states = allocMatrix(REALSXP, n, p));
	SEXP nm = PROTECT(allocVector(STRSXP, 3));
	SET_STRING_ELT(nm, 0, mkChar("values"));
	SET_STRING_ELT(nm, 1, mkChar("resid"));
	SET_STRING_ELT(nm, 2, mkChar("states"));
	setAttrib(ans, R_NamesSymbol, nm);
	UNPROTECT(1);
    }

    double sumlog = 0.0, ssq = 0.0;
    int nu = 0;
    for (int l = 0; l < n; l++) {
	for (int i = 0; i < p; i++) {
	    double tmp = 0.0;
	    for (int k = 0; k < p; k++)
		tmp += T[i + p * k] * a[k];
	    anew[i] = tmp;
	}
	if (l > asInteger(sUP)) {
	    for (int i = 0; i < p; i++)
		for (int j = 0; j < p; j++) {
		    double tmp = 0.0;
		    for (int k = 0; k < p; k++)
			tmp += T[i + p * k] * P[k + p * j];
		    mm[i + p * j] = tmp;
		}
	    for (int i = 0; i < p; i++)
		for (int j = 0; j < p; j++) {
		    double tmp = V[i + p * j];
		    for (int k = 0; k < p; k++)
			tmp += mm[i + p * k] * T[j + p * k];
		    Pnew[i + p * j] = tmp;
		}
	}
	if (!ISNAN(y[l])) {
	    nu++;
	    double *rr = NULL /* -Wall */;
	    if(lop) rr = REAL(resid);
	    double resid0 = y[l];
	    for (int i = 0; i < p; i++)
		resid0 -= Z[i] * anew[i];
	    double gain = h;
	    for (int i = 0; i < p; i++) {
		double tmp = 0.0;
		for (int j = 0; j < p; j++)
		    tmp += Pnew[i + j * p] * Z[j];
		M[i] = tmp;
		gain += Z[i] * M[i];
	    }
	    ssq += resid0 * resid0 / gain;
	    if(lop) rr[l] = resid0 / sqrt(gain);
	    sumlog += log(gain);
	    for (int i = 0; i < p; i++)
		a[i] = anew[i] + M[i] * resid0 / gain;
	    for (int i = 0; i < p; i++)
		for (int j = 0; j < p; j++)
		    P[i + j * p] = Pnew[i + j * p] - M[i] * M[j] / gain;
	} else {
	    double *rr = NULL /* -Wall */;
	    if(lop) rr = REAL(resid);
	    for (int i = 0; i < p; i++)
		a[i] = anew[i];
	    for (int i = 0; i < p * p; i++)
		P[i] = Pnew[i];
	    if(lop) rr[l] = NA_REAL;
	}
	if(lop) {
	    double *rs = REAL(states);
	    for (int j = 0; j < p; j++) rs[l + n*j] = a[j];
	}
    }

    SEXP res = PROTECT(allocVector(REALSXP, 2));
    REAL(res)[0] = ssq/nu; REAL(res)[1] = sumlog/nu;
    if(lop) {
	SET_VECTOR_ELT(ans, 0, res);
	if(asLogical(update)) setAttrib(ans, install("mod"), mod);
	UNPROTECT(3);
	return ans;
    } else {
	if(asLogical(update)) setAttrib(res, install("mod"), mod);
	UNPROTECT(2);
	return res;
    }
}

SEXP
KalmanSmooth(SEXP sy, SEXP mod, SEXP sUP)
{
    SEXP sZ = getListElement(mod, "Z"), sa = getListElement(mod, "a"),
	sP = getListElement(mod, "P"), sT = getListElement(mod, "T"),
	sV = getListElement(mod, "V"), sh = getListElement(mod, "h"),
	sPn = getListElement(mod, "Pn");

    if (TYPEOF(sy) != REALSXP || TYPEOF(sZ) != REALSXP ||
	TYPEOF(sa) != REALSXP || TYPEOF(sP) != REALSXP ||
	TYPEOF(sT) != REALSXP || TYPEOF(sV) != REALSXP)
	error(_("invalid argument type"));

    SEXP ssa, ssP, ssPn, res, states = R_NilValue, sN;
    int n = LENGTH(sy), p = LENGTH(sa);
    double *y = REAL(sy), *Z = REAL(sZ), *a, *P,
	*T = REAL(sT), *V = REAL(sV), h = asReal(sh), *Pnew;
    double *at, *rt, *Pt, *gains, *resids, *Mt, *L, gn, *Nt;
    Rboolean var = TRUE;

    PROTECT(ssa = duplicate(sa)); a = REAL(ssa);
    PROTECT(ssP = duplicate(sP)); P = REAL(ssP);
    PROTECT(ssPn = duplicate(sPn)); Pnew = REAL(ssPn);

    PROTECT(res = allocVector(VECSXP, 2));
    SEXP nm = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(nm, 0, mkChar("smooth"));
    SET_STRING_ELT(nm, 1, mkChar("var"));
    setAttrib(res, R_NamesSymbol, nm);
    UNPROTECT(1);
    SET_VECTOR_ELT(res, 0, states = allocMatrix(REALSXP, n, p));
    at = REAL(states);
    SET_VECTOR_ELT(res, 1, sN = allocVector(REALSXP, n*p*p));
    Nt = REAL(sN);

    double *anew, *mm, *M;
    anew = (double *) R_alloc(p, sizeof(double));
    M = (double *) R_alloc(p, sizeof(double));
    mm = (double *) R_alloc(p * p, sizeof(double));

    Pt = (double *) R_alloc(n * p * p, sizeof(double));
    gains = (double *) R_alloc(n, sizeof(double));
    resids = (double *) R_alloc(n, sizeof(double));
    Mt = (double *) R_alloc(n * p, sizeof(double));
    L = (double *) R_alloc(p * p, sizeof(double));

    for (int l = 0; l < n; l++) {
	for (int i = 0; i < p; i++) {
	    double tmp = 0.0;
	    for (int k = 0; k < p; k++)
		tmp += T[i + p * k] * a[k];
	    anew[i] = tmp;
	}
	if (l > asInteger(sUP)) {
	    for (int i = 0; i < p; i++)
		for (int j = 0; j < p; j++) {
		    double tmp = 0.0;
		    for (int k = 0; k < p; k++)
			tmp += T[i + p * k] * P[k + p * j];
		    mm[i + p * j] = tmp;
		}
	    for (int i = 0; i < p; i++)
		for (int j = 0; j < p; j++) {
		    double tmp = V[i + p * j];
		    for (int k = 0; k < p; k++)
			tmp += mm[i + p * k] * T[j + p * k];
		    Pnew[i + p * j] = tmp;
		}
	}
	for (int i = 0; i < p; i++) at[l + n*i] = anew[i];
	for (int i = 0; i < p*p; i++) Pt[l + n*i] = Pnew[i];
	if (!ISNAN(y[l])) {
	    double resid0 = y[l];
	    for (int i = 0; i < p; i++)
		resid0 -= Z[i] * anew[i];
	    double gain = h;
	    for (int i = 0; i < p; i++) {
		double tmp = 0.0;
		for (int j = 0; j < p; j++)
		    tmp += Pnew[i + j * p] * Z[j];
		Mt[l + n*i] = M[i] = tmp;
		gain += Z[i] * M[i];
	    }
	    gains[l] = gain;
	    resids[l] = resid0;
	    for (int i = 0; i < p; i++)
		a[i] = anew[i] + M[i] * resid0 / gain;
	    for (int i = 0; i < p; i++)
		for (int j = 0; j < p; j++)
		    P[i + j * p] = Pnew[i + j * p] - M[i] * M[j] / gain;
	} else {
	    for (int i = 0; i < p; i++) {
		a[i] = anew[i];
		Mt[l + n * i] = 0.0;
	    }
	    for (int i = 0; i < p * p; i++)
		P[i] = Pnew[i];
	    gains[l] = NA_REAL;
	    resids[l] = NA_REAL;
	}
    }

    /* rt stores r_{t-1} */
    rt = (double *) R_alloc(n * p, sizeof(double));
    for (int l = n - 1; l >= 0; l--) {
	if (!ISNAN(gains[l])) {
	    gn = 1/gains[l];
	    for (int i = 0; i < p; i++)
		rt[l + n * i] = Z[i] * resids[l] * gn;
	} else {
	    for (int i = 0; i < p; i++) rt[l + n * i] = 0.0;
	    gn = 0.0;
	}

	if (var) {
	    for (int i = 0; i < p; i++)
		for (int j = 0; j < p; j++)
		    Nt[l + n*i + n*p*j] = Z[i] * Z[j] * gn;
	}

	if (l < n - 1) {
	    /* compute r_{t-1} */
	    for (int i = 0; i < p; i++)
		for (int j = 0; j < p; j++)
		    mm[i + p * j] = ((i==j) ? 1:0) - Mt[l + n * i] * Z[j] * gn;
	    for (int i = 0; i < p; i++)
		for (int j = 0; j < p; j++) {
		    double tmp = 0.0;
		    for (int k = 0; k < p; k++)
			tmp += T[i + p * k] * mm[k + p * j];
		    L[i + p * j] = tmp;
		}
	    for (int i = 0; i < p; i++) {
		double tmp = 0.0;
		for (int j = 0; j < p; j++)
		    tmp += L[j + p * i] * rt[l + 1 + n * j];
		rt[l + n * i] += tmp;
	    }
	    if(var) { /* compute N_{t-1} */
		for (int i = 0; i < p; i++)
		    for (int j = 0; j < p; j++) {
			double tmp = 0.0;
			for (int k = 0; k < p; k++)
			    tmp += L[k + p * i] * Nt[l + 1 + n*k + n*p*j];
			mm[i + p * j] = tmp;
		    }
		for (int i = 0; i < p; i++)
		    for (int j = 0; j < p; j++) {
			double tmp = 0.0;
			for (int k = 0; k < p; k++)
			    tmp += mm[i + p * k] * L[k + p * j];
			Nt[l + n*i + n*p*j] += tmp;
		    }
	    }
	}

	for (int i = 0; i < p; i++) {
	    double tmp = 0.0;
	    for (int j = 0; j < p; j++)
		tmp += Pt[l + n*i + n*p*j] * rt[l + n * j];
	    at[l + n*i] += tmp;
	}
    }
    if (var)
	for (int l = 0; l < n; l++) {
	    for (int i = 0; i < p; i++)
		for (int j = 0; j < p; j++) {
		    double tmp = 0.0;
		    for (int k = 0; k < p; k++)
			tmp += Pt[l + n*i + n*p*k] * Nt[l + n*k + n*p*j];
		    mm[i + p * j] = tmp;
		}
	    for (int i = 0; i < p; i++)
		for (int j = 0; j < p; j++) {
		    double tmp = Pt[l + n*i + n*p*j];
		    for (int k = 0; k < p; k++)
			tmp -= mm[i + p * k] * Pt[l + n*k + n*p*j];
		    Nt[l + n*i + n*p*j] = tmp;
		}
	}
    UNPROTECT(4);
    return res;
}


SEXP
KalmanFore(SEXP nahead, SEXP mod, SEXP update)
{
    mod = PROTECT(duplicate(mod));
    SEXP sZ = getListElement(mod, "Z"), sa = getListElement(mod, "a"),
	sP = getListElement(mod, "P"), sT = getListElement(mod, "T"),
	sV = getListElement(mod, "V"), sh = getListElement(mod, "h");

    if (TYPEOF(sZ) != REALSXP ||
	TYPEOF(sa) != REALSXP || TYPEOF(sP) != REALSXP ||
	TYPEOF(sT) != REALSXP || TYPEOF(sV) != REALSXP)
	error(_("invalid argument type"));

    int  n = asInteger(nahead), p = LENGTH(sa);
    double *Z = REAL(sZ), *a = REAL(sa), *P = REAL(sP), *T = REAL(sT),
	*V = REAL(sV), h = asReal(sh);
    double *mm, *anew, *Pnew;

    anew = (double *) R_alloc(p, sizeof(double));
    Pnew = (double *) R_alloc(p * p, sizeof(double));
    mm = (double *) R_alloc(p * p, sizeof(double));
    SEXP res, forecasts, se;
    PROTECT(res = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(res, 0, forecasts = allocVector(REALSXP, n));
    SET_VECTOR_ELT(res, 1, se = allocVector(REALSXP, n));
    {
	SEXP nm = PROTECT(allocVector(STRSXP, 2));
	SET_STRING_ELT(nm, 0, mkChar("pred"));
	SET_STRING_ELT(nm, 1, mkChar("var"));
	setAttrib(res, R_NamesSymbol, nm);
	UNPROTECT(1);
    }
    for (int l = 0; l < n; l++) {
	double fc = 0.0;
	for (int i = 0; i < p; i++) {
	    double tmp = 0.0;
	    for (int k = 0; k < p; k++)
		tmp += T[i + p * k] * a[k];
	    anew[i] = tmp;
	    fc += tmp * Z[i];
	}
	for (int i = 0; i < p; i++)
	    a[i] = anew[i];
	REAL(forecasts)[l] = fc;

	for (int i = 0; i < p; i++)
	    for (int j = 0; j < p; j++) {
		double tmp = 0.0;
		for (int k = 0; k < p; k++)
		    tmp += T[i + p * k] * P[k + p * j];
		mm[i + p * j] = tmp;
	    }
	for (int i = 0; i < p; i++)
	    for (int j = 0; j < p; j++) {
		double tmp = V[i + p * j];
		for (int k = 0; k < p; k++)
		    tmp += mm[i + p * k] * T[j + p * k];
		Pnew[i + p * j] = tmp;
	    }
	double tmp = h;
	for (int i = 0; i < p; i++)
	    for (int j = 0; j < p; j++) {
		P[i + j * p] = Pnew[i + j * p];
		tmp += Z[i] * Z[j] * P[i + j * p];
	    }
	REAL(se)[l] = tmp;
    }
    if(asLogical(update)) setAttrib(res, install("mod"), mod);
    UNPROTECT(2);
    return res;
}


static void partrans(int p, double *raw, double *new)
{
    int j, k;
    double a, work[100];

    if(p > 100) error(_("can only transform 100 pars in arima0"));

    /* Step one: map (-Inf, Inf) to (-1, 1) via tanh
       The parameters are now the pacf phi_{kk} */
    for(j = 0; j < p; j++) work[j] = new[j] = tanh(raw[j]);
    /* Step two: run the Durbin-Levinson recursions to find phi_{j.},
       j = 2, ..., p and phi_{p.} are the autoregression coefficients */
    for(j = 1; j < p; j++) {
	a = new[j];
	for(k = 0; k < j; k++)
	    work[k] -= a * new[j - k - 1];
	for(k = 0; k < j; k++) new[k] = work[k];
    }
}

SEXP ARIMA_undoPars(SEXP sin, SEXP sarma)
{
    int *arma = INTEGER(sarma), mp = arma[0], mq = arma[1], msp = arma[2],
	v, n = LENGTH(sin);
    double *params, *in = REAL(sin);
    SEXP res = allocVector(REALSXP, n);

    params = REAL(res);
    for (int i = 0; i < n; i++) params[i] = in[i];
    if (mp > 0) partrans(mp, in, params);
    v = mp + mq;
    if (msp > 0) partrans(msp, in + v, params + v);
    return res;
}


SEXP ARIMA_transPars(SEXP sin, SEXP sarma, SEXP strans)
{
    int *arma = INTEGER(sarma), trans = asLogical(strans);
    int mp = arma[0], mq = arma[1], msp = arma[2], msq = arma[3],
	ns = arma[4], i, j, p = mp + ns * msp, q = mq + ns * msq, v;
    double *in = REAL(sin), *params = REAL(sin), *phi, *theta;
    SEXP res, sPhi, sTheta;

    PROTECT(res = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(res, 0, sPhi = allocVector(REALSXP, p));
    SET_VECTOR_ELT(res, 1, sTheta = allocVector(REALSXP, q));
    phi = REAL(sPhi);
    theta = REAL(sTheta);
    if (trans) {
	int n = mp + mq + msp + msq;

	params = (double *) R_alloc(n, sizeof(double));
	for (i = 0; i < n; i++) params[i] = in[i];
	if (mp > 0) partrans(mp, in, params);
	v = mp + mq;
	if (msp > 0) partrans(msp, in + v, params + v);
    }
    if (ns > 0) {
	/* expand out seasonal ARMA models */
	for (i = 0; i < mp; i++) phi[i] = params[i];
	for (i = 0; i < mq; i++) theta[i] = params[i + mp];
	for (i = mp; i < p; i++) phi[i] = 0.0;
	for (i = mq; i < q; i++) theta[i] = 0.0;
	for (j = 0; j < msp; j++) {
	    phi[(j + 1) * ns - 1] += params[j + mp + mq];
	    for (i = 0; i < mp; i++)
		phi[(j + 1) * ns + i] -= params[i] * params[j + mp + mq];
	}
	for (j = 0; j < msq; j++) {
	    theta[(j + 1) * ns - 1] += params[j + mp + mq + msp];
	    for (i = 0; i < mq; i++)
		theta[(j + 1) * ns + i] += params[i + mp] *
		    params[j + mp + mq + msp];
	}
    } else {
	for (i = 0; i < mp; i++) phi[i] = params[i];
	for (i = 0; i < mq; i++) theta[i] = params[i + mp];
    }
    UNPROTECT(1);
    return res;
}

#if !defined(atanh) && defined(HAVE_DECL_ATANH) && !HAVE_DECL_ATANH
extern double atanh(double x);
#endif
static void invpartrans(int p, double *phi, double *new)
{
    int j, k;
    double a, work[100];

    if(p > 100) error(_("can only transform 100 pars in arima0"));

    for(j = 0; j < p; j++) work[j] = new[j] = phi[j];
    /* Run the Durbin-Levinson recursions backwards
       to find the PACF phi_{j.} from the autoregression coefficients */
    for(j = p - 1; j > 0; j--) {
	a = new[j];
	for(k = 0; k < j; k++)
	    work[k]  = (new[k] + a * new[j - k - 1]) / (1 - a * a);
	for(k = 0; k < j; k++) new[k] = work[k];
    }
    for(j = 0; j < p; j++) new[j] = atanh(new[j]);
}

SEXP ARIMA_Invtrans(SEXP in, SEXP sarma)
{
    int *arma = INTEGER(sarma), mp = arma[0], mq = arma[1], msp = arma[2],
	i, v, n = LENGTH(in);
    SEXP y = allocVector(REALSXP, n);
    double *raw = REAL(in), *new = REAL(y);

    for(i = 0; i < n; i++) new[i] = raw[i];
    if (mp > 0) invpartrans(mp, raw, new);
    v = mp + mq;
    if (msp > 0) invpartrans(msp, raw + v, new + v);
    return y;
}

#define eps 1e-3
SEXP ARIMA_Gradtrans(SEXP in, SEXP sarma)
{
    int *arma = INTEGER(sarma), mp = arma[0], mq = arma[1], msp = arma[2],
	n = LENGTH(in);
    SEXP y = allocMatrix(REALSXP, n, n);
    double *raw = REAL(in), *A = REAL(y), w1[100], w2[100], w3[100];

    for (int i = 0; i < n; i++)
	for (int j = 0; j < n; j++)
	    A[i + j*n] = (i == j);
    if(mp > 0) {
	for (int i = 0; i < mp; i++) w1[i] = raw[i];
	partrans(mp, w1, w2);
	for (int i = 0; i < mp; i++) {
	    w1[i] += eps;
	    partrans(mp, w1, w3);
	    for (int j = 0; j < mp; j++) A[i + j*n] = (w3[j] - w2[j])/eps;
	    w1[i] -= eps;
	}
    }
    if(msp > 0) {
	int v = mp + mq;
	for (int i = 0; i < msp; i++) w1[i] = raw[i + v];
	partrans(msp, w1, w2);
	for(int i = 0; i < msp; i++) {
	    w1[i] += eps;
	    partrans(msp, w1, w3);
	    for(int j = 0; j < msp; j++)
		A[i + v + (j+v)*n] = (w3[j] - w2[j])/eps;
	    w1[i] -= eps;
	}
    }
    return y;
}


SEXP
ARIMA_Like(SEXP sy, SEXP mod, SEXP sUP, SEXP giveResid)
{
    SEXP sPhi = getListElement(mod, "phi"),
	sTheta = getListElement(mod, "theta"),
	sDelta = getListElement(mod, "Delta"),
	sa = getListElement(mod, "a"),
	sP = getListElement(mod, "P"),
	sPn = getListElement(mod, "Pn");

    if (TYPEOF(sPhi) != REALSXP || TYPEOF(sTheta) != REALSXP ||
	TYPEOF(sDelta) != REALSXP || TYPEOF(sa) != REALSXP ||
	TYPEOF(sP) != REALSXP || TYPEOF(sPn) != REALSXP)
	error(_("invalid argument type"));

    SEXP res, nres, sResid = R_NilValue;
    int n = LENGTH(sy), rd = LENGTH(sa), p = LENGTH(sPhi),
	q = LENGTH(sTheta), d = LENGTH(sDelta), r = rd - d;
    double *y = REAL(sy), *a = REAL(sa), *P = REAL(sP), *Pnew = REAL(sPn);
    double *phi = REAL(sPhi), *theta = REAL(sTheta), *delta = REAL(sDelta);
    double sumlog = 0.0, ssq = 0, *anew, *mm = NULL, *M;
    int nu = 0;
    Rboolean useResid = asLogical(giveResid);
    double *rsResid = NULL /* -Wall */;

    anew = (double *) R_alloc(rd, sizeof(double));
    M = (double *) R_alloc(rd, sizeof(double));
    if (d > 0) mm = (double *) R_alloc(rd * rd, sizeof(double));

    if (useResid) {
	PROTECT(sResid = allocVector(REALSXP, n));
	rsResid = REAL(sResid);
    }

    for (int l = 0; l < n; l++) {
	for (int i = 0; i < r; i++) {
	    double tmp = (i < r - 1) ? a[i + 1] : 0.0;
	    if (i < p) tmp += phi[i] * a[0];
	    anew[i] = tmp;
	}
	if (d > 0) {
	    for (int i = r + 1; i < rd; i++) anew[i] = a[i - 1];
	    double tmp = a[0];
	    for (int i = 0; i < d; i++) tmp += delta[i] * a[r + i];
	    anew[r] = tmp;
	}
	if (l > asInteger(sUP)) {
	    if (d == 0) {
		for (int i = 0; i < r; i++) {
		    double vi = 0.0;
		    if (i == 0) vi = 1.0; else if (i - 1 < q) vi = theta[i - 1];
		    for (int j = 0; j < r; j++) {
			double tmp = 0.0;
			if (j == 0) tmp = vi; else if (j - 1 < q) tmp = vi * theta[j - 1];
			if (i < p && j < p) tmp += phi[i] * phi[j] * P[0];
			if (i < r - 1 && j < r - 1) tmp += P[i + 1 + r * (j + 1)];
			if (i < p && j < r - 1) tmp += phi[i] * P[j + 1];
			if (j < p && i < r - 1) tmp += phi[j] * P[i + 1];
			Pnew[i + r * j] = tmp;
		    }
		}
	    } else {
		/* mm = TP */
		for (int i = 0; i < r; i++)
		    for (int j = 0; j < rd; j++) {
			double tmp = 0.0;
			if (i < p) tmp += phi[i] * P[rd * j];
			if (i < r - 1) tmp += P[i + 1 + rd * j];
			mm[i + rd * j] = tmp;
		    }
		for (int j = 0; j < rd; j++) {
		    double tmp = P[rd * j];
		    for (int k = 0; k < d; k++)
			tmp += delta[k] * P[r + k + rd * j];
		    mm[r + rd * j] = tmp;
		}
		for (int i = 1; i < d; i++)
		    for (int j = 0; j < rd; j++)
			mm[r + i + rd * j] = P[r + i - 1 + rd * j];

		/* Pnew = mmT' */
		for (int i = 0; i < r; i++)
		    for (int j = 0; j < rd; j++) {
			double tmp = 0.0;
			if (i < p) tmp += phi[i] * mm[j];
			if (i < r - 1) tmp += mm[rd * (i + 1) + j];
			Pnew[j + rd * i] = tmp;
		    }
		for (int j = 0; j < rd; j++) {
		    double tmp = mm[j];
		    for (int k = 0; k < d; k++)
			tmp += delta[k] * mm[rd * (r + k) + j];
		    Pnew[rd * r + j] = tmp;
		}
		for (int i = 1; i < d; i++)
		    for (int j = 0; j < rd; j++)
			Pnew[rd * (r + i) + j] = mm[rd * (r + i - 1) + j];
		/* Pnew <- Pnew + (1 theta) %o% (1 theta) */
		for (int i = 0; i <= q; i++) {
		    double vi = (i == 0) ? 1. : theta[i - 1];
		    for (int j = 0; j <= q; j++)
			Pnew[i + rd * j] += vi * ((j == 0) ? 1. : theta[j - 1]);
		}
	    }
	}
	if (!ISNAN(y[l])) {
	    double resid = y[l] - anew[0];
	    for (int i = 0; i < d; i++)
		resid -= delta[i] * anew[r + i];

	    for (int i = 0; i < rd; i++) {
		double tmp = Pnew[i];
		for (int j = 0; j < d; j++)
		    tmp += Pnew[i + (r + j) * rd] * delta[j];
		M[i] = tmp;
	    }

	    double gain = M[0];
	    for (int j = 0; j < d; j++) gain += delta[j] * M[r + j];
	    if(gain < 1e4) {
		nu++;
		ssq += resid * resid / gain;
		sumlog += log(gain);
	    }
	    if (useResid) rsResid[l] = resid / sqrt(gain);
	    for (int i = 0; i < rd; i++)
		a[i] = anew[i] + M[i] * resid / gain;
	    for (int i = 0; i < rd; i++)
		for (int j = 0; j < rd; j++)
		    P[i + j * rd] = Pnew[i + j * rd] - M[i] * M[j] / gain;
	} else {
	    for (int i = 0; i < rd; i++) a[i] = anew[i];
	    for (int i = 0; i < rd * rd; i++) P[i] = Pnew[i];
	    if (useResid) rsResid[l] = NA_REAL;
	}
    }

    if (useResid) {
	PROTECT(res = allocVector(VECSXP, 3));
	SET_VECTOR_ELT(res, 0, nres = allocVector(REALSXP, 3));
	REAL(nres)[0] = ssq;
	REAL(nres)[1] = sumlog;
	REAL(nres)[2] = (double) nu;
	SET_VECTOR_ELT(res, 1, sResid);
	UNPROTECT(2);
	return res;
    } else {
	nres = allocVector(REALSXP, 3);
	REAL(nres)[0] = ssq;
	REAL(nres)[1] = sumlog;
	REAL(nres)[2] = (double) nu;
	return nres;
    }
}

/* do differencing here */
/* arma is p, q, sp, sq, ns, d, sd */
SEXP
ARIMA_CSS(SEXP sy, SEXP sarma, SEXP sPhi, SEXP sTheta,
	  SEXP sncond, SEXP giveResid)
{
    SEXP res, sResid = R_NilValue;
    double ssq = 0.0, *y = REAL(sy), tmp;
    double *phi = REAL(sPhi), *theta = REAL(sTheta), *w, *resid;
    int n = LENGTH(sy), *arma = INTEGER(sarma), p = LENGTH(sPhi),
	q = LENGTH(sTheta), ncond = asInteger(sncond);
    int ns, nu = 0;
    Rboolean useResid = asLogical(giveResid);

    w = (double *) R_alloc(n, sizeof(double));
    for (int l = 0; l < n; l++) w[l] = y[l];
    for (int i = 0; i < arma[5]; i++)
	for (int l = n - 1; l > 0; l--) w[l] -= w[l - 1];
    ns = arma[4];
    for (int i = 0; i < arma[6]; i++)
	for (int l = n - 1; l >= ns; l--) w[l] -= w[l - ns];

    PROTECT(sResid = allocVector(REALSXP, n));
    resid = REAL(sResid);
    if (useResid) for (int l = 0; l < ncond; l++) resid[l] = 0;

    for (int l = ncond; l < n; l++) {
	tmp = w[l];
	for (int j = 0; j < p; j++) tmp -= phi[j] * w[l - j - 1];
	for (int j = 0; j < min(l - ncond, q); j++)
	    tmp -= theta[j] * resid[l - j - 1];
	resid[l] = tmp;
	if (!ISNAN(tmp)) {
	    nu++;
	    ssq += tmp * tmp;
	}
    }
    if (useResid) {
	PROTECT(res = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(res, 0, ScalarReal(ssq / (double) (nu)));
	SET_VECTOR_ELT(res, 1, sResid);
	UNPROTECT(2);
	return res;
    } else {
	UNPROTECT(1);
	return ScalarReal(ssq / (double) (nu));
    }
}

SEXP TSconv(SEXP a, SEXP b)
{
    int na, nb, nab;
    SEXP ab;
    double *ra, *rb, *rab;

    PROTECT(a = coerceVector(a, REALSXP));
    PROTECT(b = coerceVector(b, REALSXP));
    na = LENGTH(a);
    nb = LENGTH(b);
    nab = na + nb - 1;
    PROTECT(ab = allocVector(REALSXP, nab));
    ra = REAL(a); rb = REAL(b); rab = REAL(ab);
    for (int i = 0; i < nab; i++) rab[i] = 0.0;
    for (int i = 0; i < na; i++)
	for (int j = 0; j < nb; j++)
	    rab[i + j] += ra[i] * rb[j];
    UNPROTECT(3);
    return (ab);
}

/* based on code from AS154 */

static void
inclu2(size_t np, double *xnext, double *xrow, double ynext,
       double *d, double *rbar, double *thetab)
{
    double cbar, sbar, di, xi, xk, rbthis, dpi;
    size_t i, k, ithisr;

/*   This subroutine updates d, rbar, thetab by the inclusion
     of xnext and ynext. */

    for (i = 0; i < np; i++) xrow[i] = xnext[i];

    for (ithisr = 0, i = 0; i < np; i++) {
	if (xrow[i] != 0.0) {
	    xi = xrow[i];
	    di = d[i];
	    dpi = di + xi * xi;
	    d[i] = dpi;
	    cbar = di / dpi;
	    sbar = xi / dpi;
	    for (k = i + 1; k < np; k++) {
		xk = xrow[k];
		rbthis = rbar[ithisr];
		xrow[k] = xk - xi * rbthis;
		rbar[ithisr++] = cbar * rbthis + sbar * xk;
	    }
	    xk = ynext;
	    ynext = xk - xi * thetab[i];
	    thetab[i] = cbar * thetab[i] + sbar * xk;
	    if (di == 0.0) return;
	} else
	    ithisr = ithisr + np - i - 1;
    }
}

#ifdef DEBUG_Q0bis
# include <R_ext/Print.h>
  double chk_V(double v[], char* nm, int jj, int len) {
    // len = length(<vector>)  <==> index must be in  {0, len-1}
    if(jj < 0 || jj >= len)
	REprintf(" %s[%2d]\n", nm, jj);
    return(v[jj]);
  }
#endif

/*
  Matwey V. Kornilov's implementation of algorithm by
  Dr. Raphael Rossignol
  See https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=14682 for details.
*/
SEXP getQ0bis(SEXP sPhi, SEXP sTheta, SEXP sTol)
{
    SEXP res;
    int p = LENGTH(sPhi), q = LENGTH(sTheta);
    double *phi = REAL(sPhi), *theta = REAL(sTheta); // tol = REAL(sTol)[0];

    int i,j, r = max(p, q + 1);

    /* Final result is block product
     *   Q0 = A1 SX A1^T + A1 SXZ A2^T + (A1 SXZ A2^T)^T + A2 A2^T ,
     * where A1 [i,j] = phi[i+j],
     *       A2 [i,j] = ttheta[i+j],  and SX, SXZ are defined below */
    PROTECT(res = allocMatrix(REALSXP, r, r));
    double *P = REAL(res);

    /* Clean P */
    Memzero(P, r*r);

#ifdef DEBUG_Q0bis
#define _ttheta(j) chk_V(ttheta, "ttheta", j, q+1)// was  r
#define _tphi(j)   chk_V(tphi,   "tphi",   j, p+1)
#define _rrz(j)    chk_V(rrz,    "rrz",    j, q)
#else
#define _ttheta(j) ttheta[j]
#define _tphi(j) tphi[j]
#define _rrz(j)  rrz [j]
#endif

    double *ttheta = (double *) R_alloc(q + 1, sizeof(double));
    /* Init ttheta = c(1, theta) */
    ttheta[0] = 1.;
    for (i = 1; i < q + 1; ++i) ttheta[i] = theta[i - 1];

    if( p > 0 ) {
	int r2 = max(p + q, p + 1);
	SEXP sgam = PROTECT(allocMatrix(REALSXP, r2, r2)),
	    sg = PROTECT(allocVector(REALSXP, r2));
	double *gam = REAL(sgam);
	double *g = REAL(sg);
	double *tphi = (double *) R_alloc(p + 1, sizeof(double));
	/* Init tphi = c(1, -phi) */
	tphi[0] = 1.;
	for (i = 1; i < p + 1; ++i) tphi[i] = -phi[i - 1];

    /* Compute the autocovariance function of U, the AR part of X */

    /* Gam := C1 + C2 ; initialize */
	Memzero(gam, r2*r2);

    /* C1[E] */
	for (j = 0; j < r2; ++j)
	    for (i = j; i < r2 && i - j < p + 1; ++i)
		gam[j*r2 + i] += _tphi(i-j);

    /* C2[E] */
	for (i = 0; i < r2; ++i)
	    for (j = 1; j < r2 && i + j < p + 1; ++j)
		gam[j*r2 + i] += _tphi(i+j);

    /* Initialize g = (1 0 0 .... 0) */
	g[0] = 1.;
	for (i = 1; i < r2; ++i)
	    g[i] = 0.;

    /* rU = solve(Gam, g)  -> solve.default() -> .Internal(La_solve, .,)
     * --> fiddling with R-objects -> C and then F77_CALL(.) of dgesv, dlange, dgecon
     * FIXME: call these directly here, possibly even use 'info' instead of error(.)
     * e.g., in case of exact singularity.
     */
	SEXP callS = PROTECT(lang4(install("solve.default"), sgam, sg, sTol)),
	    su = PROTECT(eval(callS, R_BaseEnv));
	double *u = REAL(su);
    /* SX = A SU A^T */
    /* A[i,j]  = ttheta[j-i] */
    /* SU[i,j] = u[abs(i-j)] */
    /* Q0 += ( A1 SX A1^T == A1 A SU A^T A1^T) */
	// (relying on good compiler optimization here:)
	for (i = 0; i < r; ++i)
	    for (j = i; j < r; ++j)
		for (int k = 0; i + k < p; ++k)
		    for (int L = k; L - k < q + 1; ++L)
			for (int m = 0; j + m < p; ++m)
			    for (int n = m; n - m < q + 1; ++n)
				P[r*i + j] += phi[i + k] * phi[j + m] *
				    _ttheta(L - k) * _ttheta(n - m) * u[abs(L - n)];
	UNPROTECT(4);
    /* Compute correlation matrix between X and Z */
    /* forwardsolve(C1, g) */
    /* C[i,j] = tphi[i-j] */
    /* g[i] = _ttheta(i) */
	double *rrz = (double *) R_alloc(q, sizeof(double));
	if(q > 0) {
	    for (i = 0; i < q; ++i) {
		rrz[i] = _ttheta(i);
		for (j = max(0, i - p); j < i; ++j)
		    rrz[i] -= _rrz(j) * _tphi(i-j);
	    }
	}

    /* Q0 += A1 SXZ A2^T + (A1 SXZ A2^T)^T */
    /* SXZ[i,j] = rrz[j-i-1], j > 0 */
	for (i = 0; i < r; ++i)
	    for (j = i; j < r; ++j) {
		int k, L;
		for (k = 0; i + k < p; ++k)
		    for (L = k+1; j + L < q + 1; ++L)
			P[r*i + j] += phi[i + k] * _ttheta(j + L) * _rrz(L - k - 1);
		for (k = 0; j + k < p; ++k)
		    for (L = k+1; i + L < q + 1; ++L)
			P[r*i + j] += phi[j + k] * _ttheta(i + L) * _rrz(L - k - 1);
	    }
    } // end if(p > 0)

    /* Q0 += A2 A2^T */
    for (i = 0; i < r; ++i)
	for (j = i; j < r; ++j)
	    for (int k = 0; j + k < q + 1; ++k)
		 P[r*i + j] += _ttheta(i + k) * _ttheta(j + k);

    /* Symmetrize result */
    for (i = 0; i < r; ++i)
	for (j = i+1; j < r; ++j)
	    P[r*j + i] = P[r*i + j];

    UNPROTECT(1);
    return res;
}

SEXP getQ0(SEXP sPhi, SEXP sTheta)
{
    SEXP res;
    int  p = LENGTH(sPhi), q = LENGTH(sTheta);
    double *phi = REAL(sPhi), *theta = REAL(sTheta);

    /* thetab[np], xnext[np], xrow[np].  rbar[rbar] */
    /* NB: nrbar could overflow */
    int r = max(p, q + 1);
    size_t np = r * (r + 1) / 2, nrbar = np * (np - 1) / 2, npr, npr1;
    size_t indi, indj, indn, i, j, ithisr, ind, ind1, ind2, im, jm;


    /* This is the limit using an int index.  We could use
       size_t and get more on a 64-bit system,
       but there seems no practical need. */
    if(r > 350) error(_("maximum supported lag is 350"));
    double *xnext, *xrow, *rbar, *thetab, *V;
    xnext = (double *) R_alloc(np, sizeof(double));
    xrow = (double *) R_alloc(np, sizeof(double));
    rbar = (double *) R_alloc(nrbar, sizeof(double));
    thetab = (double *) R_alloc(np, sizeof(double));
    V = (double *) R_alloc(np, sizeof(double));
    for (ind = 0, j = 0; j < r; j++) {
	double vj = 0.0;
	if (j == 0) vj = 1.0; else if (j - 1 < q) vj = theta[j - 1];
	for (i = j; i < r; i++) {
	    double vi = 0.0;
	    if (i == 0) vi = 1.0; else if (i - 1 < q) vi = theta[i - 1];
	    V[ind++] = vi * vj;
	}
    }

    PROTECT(res = allocMatrix(REALSXP, r, r));
    double *P = REAL(res);

    if (r == 1) {
	if (p == 0) P[0] = 1.0; // PR#16419
	else P[0] = 1.0 / (1.0 - phi[0] * phi[0]);
	UNPROTECT(1);
	return res;
    }
    if (p > 0) {
/*      The set of equations s * vec(P0) = vec(v) is solved for
	vec(P0).  s is generated row by row in the array xnext.  The
	order of elements in P is changed, so as to bring more leading
	zeros into the rows of s. */

	for (i = 0; i < nrbar; i++) rbar[i] = 0.0;
	for (i = 0; i < np; i++) {
	    P[i] = 0.0;
	    thetab[i] = 0.0;
	    xnext[i] = 0.0;
	}
	ind = 0;
	ind1 = -1;
	npr = np - r;
	npr1 = npr + 1;
	indj = npr;
	ind2 = npr - 1;
	for (j = 0; j < r; j++) {
	    double phij = (j < p) ? phi[j] : 0.0;
	    xnext[indj++] = 0.0;
	    indi = npr1 + j;
	    for (i = j; i < r; i++) {
		double ynext = V[ind++];
		double phii = (i < p) ? phi[i] : 0.0;
		if (j != r - 1) {
		    xnext[indj] = -phii;
		    if (i != r - 1) {
			xnext[indi] -= phij;
			xnext[++ind1] = -1.0;
		    }
		}
		xnext[npr] = -phii * phij;
		if (++ind2 >= np) ind2 = 0;
		xnext[ind2] += 1.0;
		inclu2(np, xnext, xrow, ynext, P, rbar, thetab);
		xnext[ind2] = 0.0;
		if (i != r - 1) {
		    xnext[indi++] = 0.0;
		    xnext[ind1] = 0.0;
		}
	    }
	}

	ithisr = nrbar - 1;
	im = np - 1;
	for (i = 0; i < np; i++) {
	    double bi = thetab[im];
	    for (jm = np - 1, j = 0; j < i; j++)
		bi -= rbar[ithisr--] * P[jm--];
	    P[im--] = bi;
	}

/*        now re-order p. */

	ind = npr;
	for (i = 0; i < r; i++) xnext[i] = P[ind++];
	ind = np - 1;
	ind1 = npr - 1;
	for (i = 0; i < npr; i++) P[ind--] = P[ind1--];
	for (i = 0; i < r; i++) P[i] = xnext[i];
    } else {

/* P0 is obtained by backsubstitution for a moving average process. */

	indn = np;
	ind = np;
	for (i = 0; i < r; i++)
	    for (j = 0; j <= i; j++) {
		--ind;
		P[ind] = V[ind];
		if (j != 0) P[ind] += P[--indn];
	    }
    }
    /* now unpack to a full matrix */
    for (i = r - 1, ind = np; i > 0; i--)
	for (j = r - 1; j >= i; j--)
	    P[r * i + j] = P[--ind];
    for (i = 0; i < r - 1; i++)
	for (j = i + 1; j < r; j++)
	    P[i + r * j] = P[j + r * i];
    UNPROTECT(1);
    return res;
}
*/

//import algebra.ring.{Ring}
import cats.data._
import com.mind_era.arima.OrdinaryLeastSquares.XY
import scalin.algos.RankFactorization
import spire.implicits._
import spire.math._
import scalin.{Mat, Pivot, Vec}
import scalin.syntax.assign._
import scalin.mutable._
import scalin.mutable.dense._
import scribe.Logging
import spire.algebra._
import spire.math.poly.RootFinder

import scala.collection.{breakOut, mutable}
import scala.collection.immutable.BitSet
import scala.math.ScalaNumber
import scala.reflect.ClassTag

case class Order(p: Natural, d: Natural, q: Natural) {
  def _1: Natural = p
  def _2: Natural = d
  def _3: Natural = q
}

object Order {
  def apply(): Order = Order(Natural.zero, Natural.zero, Natural.zero)
}

case class Seasonal(order: Order, period: Option[Natural])

sealed trait ArimaLearningMethod

object ArimaLearningMethod {

  case object CssMl extends ArimaLearningMethod

  case object Css extends ArimaLearningMethod

  case object Ml extends ArimaLearningMethod

}

sealed trait SsInit

object SsInit {

  case object Gardner1980 extends SsInit

  case object Rossignol2011 extends SsInit

}

sealed trait OptimizationMethod

object OptimizationMethod {

  case object BFGS extends OptimizationMethod

}

import com.mind_era.arima.ArimaLearningMethod._
import com.mind_era.arima.SsInit._
import com.mind_era.arima.OptimizationMethod._

case class ArimaInternal[V: Eq : Field](phi: IndexedSeq[V], theta: IndexedSeq[V], t: Mat[V], pn: Mat[V],
                                        a: IndexedSeq[V])

case class ArimaResult[V: Eq : Field](coefficients: IndexedSeq[V], sigma2: V, varCoef: V, mask: BitSet,
                                      logLikelihood: V, aic: V, arma: IndexedSeq[Natural], residuals: IndexedSeq[V],
                                      series: IndexedSeq[V], convergenceCode: Int, numConditionals: Natural,
                                      numObservations: Natural, model: ArimaInternal[V])

/**
  * Created by aborg on 20/05/2017.
  */
case class Arima[@specialized(Double) V: Eq : Field : Pivot : ClassTag : RootFinder: spire.algebra.Order: Trig,
ErrV <: ScalaNumber : Field : NRoot](
                                      x: Vec[V], order: Order = Order(),
                                      seasonalParam: Option[Seasonal] = Some(Seasonal(Order(), None)),
                                      var xReg: Option[Mat[V]] = None, includeMean: Boolean = true,
                                      var transformPars: Boolean = true, var fixed: Option[IndexedSeq[Option[V]]] = None,
                                      init: Option[IndexedSeq[V]] = None,
                                      method: ArimaLearningMethod = CssMl, nCondOpt: Option[Natural],
                                      ssInit: SsInit = Rossignol2011,
                                      optimizationMethod: OptimizationMethod = BFGS,
                                      var optimizationControl: Map[String, Any] = Map.empty, kappa: Double = 1e6)(
  implicit val conv: Convert[V, ErrV]) extends Logging {

  import Arima._

  val series: DenseVec[V] = x.toVec
  val n: Int = series.length
  val seasonal: Seasonal = seasonalValue(seasonalParam)
  private val seasonalPeriod = seasonal.period.get
  val arma: IndexedSeq[Natural] = Vector(order.p, order.q, seasonal.order.p, seasonal.order.q, seasonalPeriod,
    order.d, seasonal.order.d)
  val Seq(armaR1/*p*/, armaR2/*q*/, armaR3/*sp*/, armaR4/*sq*/, armaR5/*sper*/, armaR6/*d*/, armaR7/*sd*/) = arma
  val nArma: Natural = arma.take(4).reduce(_ + _)
  val `1, -1`: Vec[V] = oneMinusOne
  val `1, seasonalPeriod-1 0s, -1`: DenseVec[V] = Vec.ones[V](1) + Vec.zeros[V](seasonalPeriod.toInt - 1) + Vec[V](
    Ring.negate(Ring.one[V]))
  val delta: Vec[V] = -Vec.fromSeq((Natural.one to seasonal.order.d).foldLeft(
    (Natural.one to order.d).foldLeft(Vec[V](Ring.one[V]))((acc, _) => convolutionVec(acc, `1, -1`)))(
    (acc, _) => convolutionVec(acc, `1, seasonalPeriod-1 0s, -1`)
  ).toIndexedSeq.tail)
  val nd: Natural = order.d + seasonal.order.d
  // Assuming all data is present
  var nUsed: Int = x.length - delta.length
  xReg.foreach(m => require(m.nRows == x.length, s"Wrong length for xReg: ${m.nRows} for x: ${x.length}"))
  val ncxreg: Int = xReg.map(mat => mat.nCols).getOrElse(0)
  if (includeMean && nd == 0) {
    xReg = Some(if (xReg.isEmpty)
      Mat.tabulate(n, 1)((_, _) => Ring.one[V])
    else
      Mat.tabulate(n, 1)((_, _) => Ring.one[V]).horzcat(xReg.get))
  }
  // No missing values in x, so no need to change CssMl to Ml in that case

  val nCond: Natural = method match {
    case Css | CssMl =>
      val nCond1 = order.p + seasonalPeriod * seasonal.order.p
      order.d + seasonal.order.d * seasonalPeriod + (nCond1 max nCondOpt.getOrElse(Natural.zero))
    case Ml => Natural.zero
  }
  fixed.foreach(seq => require(seq.length == nArma + ncxreg, s"Wrong length for fix: ${seq.length}, should be ${nArma + ncxreg}"))
  fixed = fixed.orElse(Some(Vector.fill((nArma + ncxreg).toInt)(None)))
  val mask: IndexedSeq[Boolean] = fixed.get.map(_.nonEmpty)
  val noOptim: Boolean = mask.exists(!_)
  if (noOptim) transformPars = false
  var ind: IndexedSeq[Natural] = _
  if (transformPars) {
    ind = ((arma(0) + arma(1)).toInt until arma.take(3).reduce(_ + _).toInt).map(n => Natural(n))
    if (mask.take(arma(0).toInt).exists(v => v) || ind.view.map(i =>mask(i.toInt)).exists(v => v)) {
      logger.warn("Some AR parameters are fixed, turning off transformPars")
      transformPars = false
    }
  }
  // Using after initialization init0 as init!
  var init0: Vec[V] = Vec.tabulate(nArma.toInt)(_ => Ring.zero[V])
  var parScale: Vec[ErrV] = Vec.tabulate(nArma.toInt)(_ => Ring.one[ErrV])


  if (ncxreg > 0) {
    val origXreg: Boolean = ncxreg == 1 || (nArma.toInt until (nArma.toInt + ncxreg)).exists(i => mask(i))
    if (!origXreg) {
      val SVD(_, _, v) = svd(xReg.get)
      xReg = xReg.map(mat => mat * v)
    }
    var dx: Vec[V] = x
    var dxReg: Mat[V] = xReg.get
    if (order.d > Natural.zero) {
      dx = diff(dx, Natural.one, order.d)
      dxReg = diff(dxReg, Natural.one, order.d)
    }
    if (seasonalPeriod > Natural.one && seasonal.order.d > Natural.zero) {
      dx = diff(dx, seasonalPeriod, seasonal.order.d)
      dxReg = diff(dxReg, seasonalPeriod, seasonal.order.d)
    }
    val fit: OrdinaryLeastSquares.Result[V, OrdinaryLeastSquares.CoefficientErrors[V, ErrV]] = (if (dx.length > dxReg.nCols) {
      val ols = OrdinaryLeastSquaresWithoutIntercept.olsWithErrors[V, ErrV](dx.toIndexedSeq.zip(dxReg.rowSeq).map{case (y, xs) => XY(xs.toIndexedSeq, y)})
      if (ols.beta.exists(c => c.value != Ring.zero[V])) Some(ols) else None
    } else {
      xReg.map(mat => OrdinaryLeastSquaresWithoutIntercept.olsWithErrors[V, ErrV](x.toIndexedSeq.zip(mat.rowSeq).map{case (y, xs) => XY(xs.toIndexedSeq, y)}))
    }).getOrElse(OrdinaryLeastSquaresWithoutIntercept.olsWithErrors[V, ErrV](x.toIndexedSeq.zip(xReg.get.rowSeq).map{case (y, xs) => XY(xs.toIndexedSeq, y)}))
    val sumNAs = 0
    nUsed = x.length - sumNAs - delta.length
    init0 = init0 ++ Vec(fit.beta.map(_.value): _*)
    // single value
    val ses10: Vec[ErrV] = Vec(fit.beta.map(coeff => coeff.stdError * 10d): _*)
    parScale = parScale ++ ses10
    //dxReg.inverse
  }
  if (nUsed <= 0) throw new IllegalArgumentException("Too few observations")
  init.foreach(initV => {
    if (initV.length != init0.length)
      throw new IllegalStateException(s"`init` has wrong length: ${initV.length} vs ${init0.length}")
    // Missing values -not supported- in init should come from init0

    method match {
      case Ml =>
        if (arma(0) > Natural.zero && !arCheck(initV.take(arma(0).toInt)))
          throw new IllegalArgumentException("Non-stationary AR part")
        if (armaR3 > Natural.zero && !arCheck(initV.slice((armaR1 + armaR2).toInt, (armaR1 + armaR2 + armaR3).toInt)))
          throw new IllegalArgumentException("Non-stationary seasional AR part")
        if (transformPars) {
          init0 = arimaInvTrans(initV, arma)
        }
      case _ => // Nothing to do
    }
  })
  def armaCss(p: IndexedSeq[V]): (ErrV, IndexedSeq[V]) = {
    val par = fixed.getOrElse(IndexedSeq.empty).toBuffer
    var idx = -1
    par.indices.foreach(i => if (mask(i)) {
      par.update(i, Some(p({idx += 1; idx})))
    })
    val trArma = transPars(par.map(_.get).toIndexedSeq, arma, trans = false)
    if (ncxreg > 0) {
      // TODO convolution for Mat?
      x -= convolutionVec(xReg.get.rowSeq.head, Vec(par.slice(nArma.toInt, nArma.toInt + ncxreg).map(_.get): _*))
    }
    ???
  }
  // init0 is to be used from now on
  private[this] var coef: IndexedSeq[Option[V]] = fixed.get
  if (optimizationControl.contains("parscale")) {
    optimizationControl = optimizationControl.updated("parscale", mask)
  }
  if (method == ArimaLearningMethod.Css) {

  } else { // CSS ML or ML
    if (method == ArimaLearningMethod.CssMl) {

    }

  }

  def arCheck(ar: IndexedSeq[V]): Boolean = {
    val p = ar.lastIndexWhere(_ != Ring.zero[V])
    p == 0 || {
      val coeffs: Map[Int, V] = (ar.take(p).zipWithIndex.map{case (v, i) => (i + 1) -> -v} :+ (0 -> Ring.one[V]))(
        breakOut)
      // TODO find complex roots and check their length
      Polynomial(coeffs).roots.forall(_ > Ring.one[V])
    }
  }
}

object Arima {
  import spire.syntax._
  def atanh[V: Trig: Field](x: V): V = log((Ring.one[V] + x) / (Ring.one[V] - x)) / (Ring.one[V] + Ring.one[V])

  private def oneMinusOne[V: Ring]: Vec[V] = Vec(Ring.one, Ring.negate(Rig.one))
  private def seasonalValue(seasonalParam: Option[Seasonal]): Seasonal = {
    seasonalParam.map(s => s.copy(period = Some(s.period.map(p =>
      if (p == Natural.zero)
        Natural.one
      else
        p).getOrElse(Natural.one)))).getOrElse(Seasonal(Order(), Some(Natural.one)))
  }

  private def convolution[V:Rig](as: IndexedSeq[V], bs: IndexedSeq[V]): IndexedSeq[V] = {
    val res = mutable.Buffer.fill[V](as.length + bs.length - 1)(Rig.zero[V])
    for (i <- as.indices) {
      for (j <- bs.indices) {
        res(i + j) += as(i) * bs(j)
      }
    }
    res.toIndexedSeq
  }

  private def convolutionVec[V: Rig](as: Vec[V], bs: Vec[V]): Vec[V] = {
    val res: DenseVec[V] = scalin.mutable.DenseVec.fillConstant(as.length + bs.length - 1)(Rig.zero)
    for (i <- 0 until as.length)
      for (j <- 0 until bs.length)
        res.set(i + j,  as(i) * bs(j))
    res
  }

  case class SVD[V: Eq: Field](u: Mat[V], d: Vec[V]/*Diag[V] would be better*/, v: Mat[V])

  // TODO implement
  def svd[V: Eq: Field](matrix: Mat[V]): SVD[V] = ???
  def diff[V: Eq: Field](dx: Vec[V], lag: Natural, d: Natural): Vec[V] = {
    require(lag > Natural.zero)
    require(d > Natural.zero)
    d match {
      case Natural.one => Vec.fromSeq((0 until dx.length - lag.toInt).map(i => dx(i + lag.toInt) - dx(i)))
      case _ => diff(diff(dx, lag, d - Natural.one), lag, Natural.one)
    }
  }
  def diff[V: Eq: Field](dx: Mat[V], lag: Natural, d: Natural): Mat[V] = {
    dx.colSeq.map(diff(_, lag, d).toColMat).reduce((col1, col2) => col1.horzcat(col2))
  }

  def invParTrans[V: Field: Trig](p: Natural, phi: IndexedSeq[V]): scala.IndexedSeq[V] = {
//    import spire.syntax.ring._
    require(p > Natural.zero)
    val res = phi.toBuffer
    val work = phi.toBuffer
    var a: V = Ring.zero[V]
    for(j <- Range.inclusive(p.toInt - 1, 1, -1)) {
      a = res(j)
      for(k <- 0 until j) {
        work(k) = (res(k) + a * res(j - k - 1)) / (1 - a * a)
      }
      for(k <- 0 until j) res(k) = work(k)
    }
    for (j <- 0 until p.toInt) {
      res.update(j, atanh(res(j)))
    }
    res.toIndexedSeq
  }

  case class PhiTheta[V](phi: IndexedSeq[V], theta: IndexedSeq[V])

  def transPars[V: Field: Trig](in: IndexedSeq[V], arma: IndexedSeq[Natural], trans: Boolean): PhiTheta[V] = {
    import _root_.scala.collection.mutable.{IndexedSeq => MutIndexedSeq}
    val IndexedSeq(mp, mq, msp, msq, ns, _*) = arma
    val p = mp + ns * msp
    val q = mq + ns * msq
    var params: IndexedSeq[V] = in.map(v => v)
    if (trans) {
      val n = mp + mq + msp + msq
      if (mp > Natural.zero) {
        params = parTrans(mp, in) ++ params.drop(mp.toInt)
      }
      val v = (mp + mq).toInt
      if (msp > Natural.zero) {
        params = params.take(v) ++ parTrans(msp, in.drop(v)) ++ params.drop(v + msp.toInt)
      }
    }
    if (ns > Natural.zero) {
      val phi = MutIndexedSeq(params.take(mp.toInt) ++ Seq.fill((p - mp).toInt)(Field.zero[V]): _*)
      val theta = MutIndexedSeq(params.slice(mp.toInt, (mp + mq).toInt) ++ Seq.fill((q- mq).toInt)(Field.zero[V]): _*)
      for (j <- 0 until msp.toInt) {
        phi.updated((j + 1) * ns.toInt - 1, phi((j + 1) * ns.toInt - 1) + params(j +(mp + mq).toInt))
        for (i <- 0 until mp.toInt) {
          phi.updated((j + 1) * ns.toInt + i, phi((j + 1) * ns.toInt + i) - params(i) * params(j + (mp + mq).toInt))
        }
      }
      for (j <- 0 until msq.toInt) {
        theta.updated((j + 1) * ns.toInt - 1, theta((j + 1) * ns.toInt - 1)+ params(j + (mp + mq + msp).toInt))
        for (i <- 0 until mq.toInt)
          theta.updated((j + 1) * ns.toInt +i, theta((j + 1) * ns.toInt +i) + params(i + mp.toInt) * params(j + (mp + mq + msp).toInt))

      }
      PhiTheta(phi = IndexedSeq(phi: _*), theta = IndexedSeq(theta: _*))
    } else {
      PhiTheta(phi = params.take(mp.toInt), theta = params.slice(mp.toInt, (mp + mq).toInt))
    }
  }

  def parTrans[V: Field: Trig](p: Natural, raw: IndexedSeq[V]): IndexedSeq[V] = {
    val res = raw.toBuffer
    for (j <- 0 until p.toInt) res(j) = tanh(res(j))
    val work = mutable.Buffer[V]()
    res.copyToBuffer(work)
    for (j <- 1 until p.toInt) {
      val a = res(j)
      for (k <- 0 until j) {
        work(k) -= a * res(j - k - 1)
      }
      for (k <- 0 until j) {
        res(k) = work(k)
      }
    }
    res.toIndexedSeq
  }

  def undoPars[V: Field: Trig](inp: IndexedSeq[V], arma: IndexedSeq[Natural]): Vec[V] = {
    val IndexedSeq(p, q, sp, _*) = arma
    var res = inp
    if (p > Natural.zero) {
      res = parTrans(p, inp)
    }
    val v = (p + q).toInt
    if (sp > Natural.zero) {
      res = res.take(v) ++ parTrans(sp, res.drop(v))
    }
    Vec(res: _*)
  }

  def arimaInvTrans[V: Field: Trig](initV: IndexedSeq[V], arma: IndexedSeq[Natural]): Vec[V] = {
    val IndexedSeq(p, q, sp, _*) = arma
    var res = initV
    if (p > 0) {
      res = invParTrans(p, initV)
    }
    val v = (p + q).toInt
    if (sp > 0) {
      res = res.take(v) ++ invParTrans(sp, initV.drop(v))
    }
    Vec(res: _*)
  }
}
