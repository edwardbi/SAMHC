#' Cut a sequence of numbers into bins with equal numbers in each bin
#'
#' @param x vector of values based on which cut points will be determined
#' @param y vector of values to be cut, default to be the same as \code{x}
#' @param breaks number of cut points
#' @param keep.inx indices of y that will be categorized as 1 or the largest bin
#'     even if their values are out of range of x
#'
#' @export
#'
rweCut <- function(x, y = x, breaks = 5, keep.inx = NULL) {
  cuts    <- quantile(x, seq(0, 1,length = breaks+1));
  cuts[1] <- cuts[1] - 0.001;
  rst     <- rep(NA, length(y));
  for (i in 2:length(cuts)) {
    inx      <- which(y > cuts[i-1] & y <= cuts[i]);
    rst[inx] <- i-1;
  }
  
  if (!is.null(keep.inx)) {
    inx <- which(y[keep.inx] <= cuts[1]);
    if (0 < length(inx)) {
      rst[keep.inx[inx]] <- 1;
    }
    
    inx <- which(y[keep.inx] > cuts[length(cuts)]);
    if (0 < length(inx)) {
      rst[keep.inx[inx]] <- length(cuts) - 1;
    }
  }
  
  rst
}


#' Summarize PS scores
#'
#' Get number of subjects and the distances of PS distributions for each PS
#' strata
#'
#' @inheritParams rwe_ps
#'
#' @param data_withps A class \code{RWE_DWITHPS} list. See \code{\link{rwe_ps}}.
#' @param min_n0 threshold for number of external subjects, below which the
#'     external data in the current stratum will be ignored by setting the PS
#'     distance to 0. Default value 10.
#' @param ... Parameters for \code{get_distance}, e.g., \code{metric} with
#'     options such as overlapping area (\code{ovl}).
#'
#' @return A class \code{RWE_PSDIST} dataframe with columns
#'   \itemize{
#' \item{Strata}{Index of stratum. 0 represents the overall information}
#' \item{N0,N1}{Number of subjects in group 0 and 1}
#' \item{N00, N10}{Number of arm 0 subjects in group 0 and 1, when arm exists}
#' \item{N01, N11}{Number of arm 1 subjects in group 0 and 1, when arm exists}
#' \item{Dist} Distance}
#'
#'
#' @examples
#'
#' dta_ps <- rwe_ps(ex_dta,
#'                  v_covs = paste("V", 1:7, sep = ""),
#'                  v_grp = "Group",
#'                  cur_grp_level = "current")
#'
#'
#' rwe_ps_dist(dta_ps, metric = "ovl")
#'
#'
#' @export
#'
rwe_ps_dist <- function(data_withps, min_n0 = 10,
                        trt_arm_level = NULL,
                        ...) {
  
  f_narm <- function(inx, dataps) {
    if (is.null(dataps[["_arm_"]]))
      return(c(length(inx), 0, 0))
    n0 <- length(which(0 == dataps[inx, "_arm_"]))
    n1 <- length(which(1 == dataps[inx, "_arm_"]))
    c(length(inx), n0, n1)
  }
  
  stopifnot(inherits(data_withps,
                     what = get_rwe_class("DWITHPS")))
  
  dataps   <- data_withps$data;
  nstrata  <- data_withps$nstrata;
  rst      <- NULL;
  for (i in seq_len(nstrata)) {
    inx_ps0 <- i == dataps[["_strata_"]] & 0 == dataps[["_grp_"]]
    inx_ps1 <- i == dataps[["_strata_"]] & 1 == dataps[["_grp_"]];
    n0_01   <- f_narm(which(inx_ps0), dataps)
    n1_01   <- f_narm(which(inx_ps1), dataps)
    
    if (!is.null(trt_arm_level) &
        !is.null(dataps[["_arm_"]])) {
      inx_ps0 <- inx_ps0 & trt_arm_level == dataps[["_arm_"]]
      inx_ps1 <- inx_ps1 & trt_arm_level == dataps[["_arm_"]]
    }
    
    ps0 <- dataps[which(inx_ps0), "_ps_"];
    ps1 <- dataps[which(inx_ps1), "_ps_"];
    
    if (0 == length(ps0) | 0 == length(ps1))
      warning("No samples in strata");
    
    if (any(is.na(c(ps0, ps1))))
      warning("NA found in propensity scores in a strata");
    
    if (length(ps0) < min_n0) {
      warning("Not enough data in the external data
                     in the current stratum.
                     External data ignored.");
      cur_dist <- 0;
    } else {
      cur_dist <- get_distance(ps0, ps1, ...)
    }
    
    rst <- rbind(rst, c(i, n0_01, n1_01, cur_dist))
  }
  
  ## overall
  inx_tot_ps0 <- which(0 == dataps[["_grp_"]])
  inx_tot_ps1 <- which(1 == dataps[["_grp_"]])
  n0_tot_01   <- f_narm(inx_tot_ps0, dataps)
  n1_tot_01   <- f_narm(inx_tot_ps1, dataps)
  
  ps0         <- dataps[inx_tot_ps0, "_ps_"];
  ps1         <- dataps[inx_tot_ps1, "_ps_"];
  all_dist    <- get_distance(ps0, ps1, ...)
  rst         <- rbind(rst, c(0, n0_tot_01, n1_tot_01, all_dist))
  
  ##return
  colnames(rst) <- c("Strata",
                     "N0", "N00", "N01",
                     "N1", "N10", "N11",
                     "Dist")
  rst           <- data.frame(rst)
  class(rst)    <- append(get_rwe_class("PSDIST"),
                          class(rst))
  
  rst
}

#' Compute distance from F0 to F1
#'
#' @param type type of distances. ovl: overlapping coefficient, kl:
#'     1/(1+Kullback-Leibler divergence)
#' @param n.bins number of bins for KL computation
#' @param epsilon small integer for Dirichlet smoothing
#'
#' @return a vector with the number of samples in group 0, the number of samples
#'     in group 1, and 1/(1+KL divergence) from group 0 to group 1 when type is
#'     kl, or the overlapping coefficient when type is ovl
#' @export
#'
rweDist <- function(sample.F0, sample.F1, n.bins = 10, type = c("ovl", "kl"), epsilon = 10^-6) {
  
  type     <- match.arg(type);
  
  smps     <- c(sample.F0, sample.F1);
  n0       <- length(sample.F0);
  n1       <- length(sample.F1);
  
  if (0 == n0 | 0 == n1)
    return(c(n0, n1, NA));
  
  if (1 == length(unique(smps))) {
    cut.smps <- rep(1, n0+n1)
    n.bins   <- 1;
    warning("Distributions for computing distances are degenerate.",
            call. = FALSE);
  } else {
    cut.smps <- rweCut(smps, breaks = n.bins);
  }
  
  rst <- 0;
  for (j in 1:n.bins) {
    n0.j <- length(which(j == cut.smps[1:n0]));
    n1.j <- length(which(j == cut.smps[(n0+1):(n0+n1)]));
    
    rst  <- rst + switch(type,
                         kl = {ep0  <- (n0.j+epsilon)/(n0 + epsilon * n.bins);
                         ep1  <- (n1.j+epsilon)/(n1 + epsilon * n.bins);
                         ep1 * log(ep1/ep0)},
                         ovl = min(n0.j/n0, n1.j/n1));
  }
  
  if ("kl" == type)
    rst <- 1/(1+rst);
  
  rst;
}



#' Distance between two distributions
#'
#' Get balance by different metrics
#'
#' @param cov0 Vector (or matrix for \code{mhb}) of samples from the first
#'     distribution
#' @param cov1 Vector (or matrix for \code{mhb}) of samples from the second
#'     distribution
#' @param metric Metrics of distances with options
#' \describe{ \item{ovl}{Overlapping area}
#'     \item{ksd}{Kullback-Leibler distance} \item{std}{Standardized difference
#'     in mean} \item{abb}{Absolute difference in mean} \item{ley}{Levy
#'     distance} \item{mhb}{Mahalanobis distance}
#' }
#'
#' @return A real value of the distance
#'
#' @examples
#'
#' x <- rnorm(100,  mean = 0, sd = 1)
#' y <- rnorm(1000, mean = 1, sd = 2)
#' get_distance(x, y, "ovl")
#' get_distance(x, y, "abd")
#'
#' @export
#'
get_distance <- function(cov0, cov1,
                         metric = c("ovl", "ksd", "std", "abd",
                                    "ley", "mhb")) {
  metric <- match.arg(metric)
  switch(metric,
         std = {
           s <- sqrt((var(cov1) + var(cov0)) / 2)
           abs(mean(cov1) - mean(cov0)) / s
         },
         abd = abs(mean(cov0) - mean(cov1)),
         ovl = metric_ovl(cov0, cov1),
         ksd = metric_ksd(cov0, cov1),
         ley = metric_ley(cov0, cov1),
         mhb = metric_mhb(cov0, cov1)
  )
}


## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
##             PRIVATE FUNCTIONS
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------

## overlapping coefficient
metric_ovl <- function(cov0, cov1) {
  cov <- c(cov0, cov1)
  if (length(unique(cov)) <= 10) {
    all_x <- c(rep(0, length(cov0)),
               rep(1, length(cov1)))
    pt    <- apply(prop.table(table(cov, all_x), 2),
                   1, min)
    
    return(sum(pt))
  }
  
  mn <- min(cov) * 1.25;
  mx <- max(cov) * 1.25;
  
  f1 <- approxfun(density(cov1, from = mn, to = mx,
                          bw = "nrd"));
  f0 <- approxfun(density(cov0, from = mn, to = mx,
                          bw = "nrd"));
  
  fn <- function(x)
    pmin(f1(x), f0(x))
  
  s <- try(integrate(fn, lower = mn, upper = mx,
                     subdivisions = 500)$value)
  
  ifelse(inherits(s, "try-error"),
         NA,
         s)
}

## K-S distance
metric_ksd <- function(cov0, cov1) {
  cov    <- c(cov0, cov1);
  cdf_1  <- ecdf(cov1);
  cdf_0  <- ecdf(cov0);
  1/max(abs(cdf_1(cov) - cdf_0(cov)))
}

## Levy distance
metric_ley <- function(cov0, cov1) {
  cov   <- c(cov0, cov1);
  cdf_1 <- ecdf(cov1);
  cdf_0 <- ecdf(cov0);
  e     <- max(abs(cdf_1(cov) - cdf_0(cov)))
  
  if (length(unique(cov)) <= 10)
    return(e)
  
  x     <- seq(min(cov), max(cov), length.out = 1000);
  check <- all(cdf_0(x - e) - e <= cdf_1(x) &
                 cdf_1(x) <= cdf_0(x + e) + e)
  
  while (check) {
    e <- e - .01
    check <- all(cdf_0(x - e) - e <= cdf_1(x) &
                   cdf_1(x) <= cdf_0(x + e) + e)
  }
  
  1/e
}

## mahalanobis balance
##
## covs should be a reduced datset that contains only those covariates
## that will be used for calculating Mahalanobis balance, for example,
## covs = dat[,1:6]
## trt should be the exposure variable,
## for example, trt=dat$X
##
metric_mhb <- function(cov0, cov1) {
  cov01 <- rbind(cov0, cov1)
  sinv  <- solve(cov(cov01))
  x0    <- colMeans(cov0)
  x1    <- colMeans(cov1)
  
  rst <- sum((t(x1 - x0) %*% sinv) * (x1 - x0))
  1/rst
}


#' Generate frequency table for factor columns
#'
#' @return a vector with the number of samples in group 0, the number of samples
#'     in group 1, and the KL divergence from group 0 to group 1
#' @export
#'
rweFreqTbl <- function(data, var.groupby, vars = NULL) {
  
  if (is.null(vars))
    vars <- colnames(data);
  
  rst <- NULL;
  for (v in vars) {
    if (!is.factor(data[[v]]))
      next;
    
    cur.freq <- data %>% count_(c(var.groupby, v)) %>%
      group_by_(.dots = var.groupby) %>%
      mutate(Sum = sum(n), Freq = n/sum(n)) %>%
      mutate_if(is.factor, as.character) %>%
      mutate(Cov = v) %>%
      rename_(Value = v);
    
    rst <- rbind(rst, data.frame(cur.freq));
  }
  
  rst
}

#' Split A into Bins
#'
#' Split A into bins to minimize the number difference between two arms
#'
#' @param ns1.trt numbers of treatment arm patients in each stratum
#' @param ns1.ctl numbers of control arm patients in each stratum
#' @param ns0     numbers of external patients in each stratum
#' @param A       number of target patients to borrow
#'
#' @export
#'
rweEvenLmbdS <- function(ns1.trt, ns1.ctl, ns0, A, init.lmbds = NULL) {
  
  f.target <- function(lmbds) {
    cl  <- c(lmbds, A - sum(lmbds));
    rst <- sum((ns1.trt - ns1.ctl - cl)^2);
    rst
  }
  
  f.gradient <- function(lmbds) {
    g <- -2*(ns1.trt[-NS] - ns1.ctl[-NS] - lmbds);
  }
  
  NS <- length(ns1.trt);
  
  stopifnot(NS == length(ns1.ctl) &
              NS == length(ns0));
  
  
  ## there is only one stratum
  if (1 == NS) {
    return(A);
  }
  
  ## multiple strata
  
  ## restrictions lmb > 0; sum(lmb) < A; lmb < ns0;
  ui <- rbind(diag(NS-1),
              -1 * diag(NS-1),
              rep(-1, NS-1),
              rep(1, NS-1));
  
  ci <- c(rep(0, NS-1),
          -ns0[-NS],
          -A,
          A - ns0[NS]);
  
  if (is.null(init.lmbds))
    init.lmbds <- rep(A/NS, NS-1);
  
  rst <- constrOptim(theta = init.lmbds,
                     f     = f.target,
                     grad  = f.gradient,
                     ui, ci, mu = 1e-04, control = list(),
                     outer.iterations = 100, outer.eps = 1e-05,
                     hessian = FALSE);
  
  c(rst$par, A - sum(rst$par))
}

#' Get weights
#'
#' @param A target number of subjects to be borrowed
#' @param m.lambda method to split A. rs: by overlapping coefficient; even: by
#'     minimizing trt and control imbalance in numbers
#'
#' @return power parameter before standardization
#'
#' @export
#'
rweGetLambda <- function(A, rs = NULL, ns1.trt = NULL, ns1.ctl = NULL, ns0,
                         m.lambda = c("rs", "even", "inverse"), ...) {
  m.lambda <- match.arg(m.lambda);
  
  rst <- switch(m.lambda,
                rs      = apply(cbind(ns0, A * rs/sum(rs)), 1, min),
                even    = rweEvenLmbdS(ns1.trt, ns1.ctl, ns0, A, ...),
                inverse = {mrs <- 1/(1-rs);
                apply(cbind(ns0, A * mrs/sum(mrs)), 1, min)})
  rst
}


#' Get weights
#'
#' @param A target number of subjects to be borrowed
#' @param m.lambda method to split A. rs: by overlapping coefficient; even: by
#'     minimizing trt and control imbalance in numbers
#'
#' @return power parameter before standardization
#'
#' @export
#'
rweGpsLambda <- function(A, ps_dist) {
  ps_dist <- ps_dist[which(ps_dist$Strata > 0), ]
  inx     <- seq(4, ncol(ps_dist), by = 2)
  
  ## average distance
  avg   <- apply(ps_dist[, inx], 1, mean)
  avg_a <- A * avg / sum(avg)
  
  ## lambdas
  lambdas <- apply(cbind(avg_a, ps_dist[, inx]),
                   1,
                   function(x) {
                     x[-1] / sum(x[-1]) * x[1]
                   })
  
  ## minimum
  rst <- apply(cbind(as.vector(t(lambdas)),
                     as.vector(data.matrix(ps_dist[, inx-1]))
  ),
  1,
  min)
  
  matrix(rst, nrow = nrow(ps_dist))
}


#' Summary statistics
#' l
#'
#'
#'
#'
#' @export
#'
rweSummary <- function(cur.m, cur.var, true.theta,  cur.ci = NULL) {
  rst      <- c(thetahat = cur.m,
                thetavar = cur.var,
                bias     = cur.m - true.theta,
                mse      = (cur.m - true.theta)^2)
  
  if (!is.null(cur.ci)) {
    range.ci <- range(cur.ci)
    rst <- c(rst,
             width    = range.ci[2] - range.ci[1],
             cover    = true.theta >= range.ci[1] & true.theta <= range.ci[2],
             lb       = as.numeric(range.ci[1]),
             ub       = as.numeric(range.ci[2]))
  }
  rst
}


#' Combining simulation results
#'
#' @param lst.rst List of simulation results. Each element represents a
#'     replication
#'
#' @export
#'
rweSimuCombine <- function(lst.rst, fun = mean, ignore.error = TRUE, ...) {
  if (ignore.error) {
    err.inx <- NULL;
    for (i in 1:length(lst.rst)) {
      if ("try-error" == class(lst.rst[[i]]))
        err.inx <- c(err.inx, i);
    }
    
    if (!is.null(err.inx))
      lst.rst <- lst.rst[-err.inx];
  }
  
  nreps       <- length(lst.rst);
  rep1        <- lst.rst[[1]];
  lst.combine <- rep(list(NULL), length(rep1));
  for (i in 1:nreps) {
    for (j in 1:length(rep1)) {
      cur.value        <- lst.rst[[i]][[j]];
      lst.combine[[j]] <- rbind(lst.combine[[j]],
                                as.vector(as.matrix(cur.value)));
    }
  }
  
  for (j in 1:length(lst.combine)) {
    cur.rst           <- apply(lst.combine[[j]], 2, fun, ...);
    dim(cur.rst)      <- dim(as.matrix(rep1[[j]]));
    dimnames(cur.rst) <- dimnames(as.matrix(rep1[[j]]));
    lst.combine[[j]]  <- cur.rst;
  }
  
  names(lst.combine) <- names(rep1);
  lst.combine
}



## --------------------------------------------------------------------------------
## --------------------------------------------------------------------------------
##             PRIVATE FUNCTIONS
## --------------------------------------------------------------------------------
## --------------------------------------------------------------------------------

get.rwe.class <- function(c.str = c("DWITHPS", "PSDIST", "D_GPS", "GPSDIST")) {
  c.str <- match.arg(c.str);
  switch(c.str,
         DWITHPS  = "RWE_DWITHPS",
         PSDIST   = "RWE_PSDIST",
         D_GPS    = "RWE_D_GPS",
         GPSDIST  = "RWE_GPSDIST")
}


make.global <- function(alist, dest.env='.GlobalEnv') {
  for (i in 1:length(alist)) {
    assign(names(alist[i]), alist[[i]], dest.env );
  }
}

expit <- function(x) {
  ex <- exp(x);
  ex/(1+ex);
}

get.xbeta <- function(covX, regCoeff) {
  if (length(regCoeff) > 0 &
      length(regCoeff) != ncol(covX)) {
    stop("Number of coefficients does not match with the design matrix.")
  }
  
  apply(covX, 1, function(x) {
    sum(x * regCoeff)}
  )
}

get.covmat <- function(StDevCovar, corrCovar) {
  n.x      <- length(StDevCovar);
  Vars     <- StDevCovar*StDevCovar;
  CovarMat <- matrix(NA, n.x, n.x);
  for (i in 1:n.x) {
    CovarMat[i,i] <- Vars[i];
    for (j in i:n.x) {
      if (j == i) {
        CovarMat[i,i] <- Vars[i];
        next;
      }
      CovarMat[i, j] <- corrCovar*StDevCovar[i]*StDevCovar[j];
      CovarMat[j, i] <- CovarMat[i, j];
    }
  }
  
  CovarMat
}

## cut covariates into categories
get_cov_cat <- function(covX, breaks = NULL) {
  f.cut <- function(x, bs) {
    if (is.null(bs))
      return(x);
    
    ibs <- sort(unique(c(-Inf, bs, Inf)));
    rst <- as.numeric(cut(x, breaks = ibs)) - 1;
    factor(rst, levels = 0:length(bs))
  }
  
  if (is.null(breaks))
    return(covX);
  
  if (is.numeric(breaks)) {
    rst <- apply(covX, 1, f.cut, breaks);
    rst <- t(rst);
  } else if (is.list(breaks)) {
    rst <- covX;
    for (i in 1:min(ncol(covX), length(breaks))) {
      rst[,i] <- f.cut(covX[,i],
                       breaks[[i]]);
    }
  }
  
  rst
}



#' Call STAN models
#'
#'
#' @param chains STAN parameter. Number of Markov chainsm
#' @param iter STAN parameter. Number of iterations
#' @param warmup STAN parameter. Number of burnin.
#' @param control STAN parameter. See \code{rstan::stan} for details.
#' @param ... other options to call STAN sampling such as \code{thin},
#'     \code{algorithm}. See \code{rstan::sampling} for details.#'
#'
#'
#' @export
#'
rweSTAN <- function(lst.data, stan.mdl = "powerp",
                    chains = 1, iter = 120000, warmup = 20000,
                    control = list(adapt_delta=0.95), ...) {

  
  # if(stan.mdl=="powerps"){
  #   #model=stan.models[[1]]
  #   model=stan_model(model_code=powerps,model_name="powerps");
  # }else if(stan.mdl=="powerp"){
  #   #model=stan.models[[2]]
  #   model<- stan_model(model_code=powerp,model_name="powerp");
  # }else if(stan.mdl=="powerpsbinary"){
  #   #model=stan.models[[3]]
  #   model<- stan_model(model_code=powerpsbinary,model_name="powerpsbinary");
  # }
  stan.rst <- sampling(model,
                              data    = lst.data,
                              chains  = chains,
                              iter    = iter,
                              warmup  = warmup,
                              control = control,refresh=0,
                              ...);

  
  stan.rst;
}



#' Get Posterior for all stratum
#'
#' @param data class DWITHPS data frame
#' @param type type of outcomes
#' @param A    target number of subjects to be borrowed
#' @param RS   parameters for dirichelet prior
#' @param Fix.RS whether treat RS as fixed or the prior mean of vs
#' @param ...  extra parameters for calling function \code{\link{rweSTAN}}
#'
#' @export
#'
rwePsPowDrawPost <- function(data, type = c("continuous", "binary"),
                             A = 0, RS = NULL, Fix.RS = FALSE,
                             v.outcome = "Y",  ...) {
  
  stopifnot(v.outcome %in% colnames(data));
  type <- match.arg(type);
  
  ## prepare stan data
  data   <- data[!is.na(data[["_strata_"]]),];
  S      <- max(data[["_strata_"]]);
  stan.d <- NULL;
  
  Y1     <- NULL;
  INX1   <- NULL;
  for (i in 1:S) {
    cur.d1 <- data[data[["_strata_"]] == i & data[["_grp_"]] == 1, v.outcome];
    cur.d0 <- data[data[["_strata_"]] == i & data[["_grp_"]] == 0, v.outcome];
    
    if (0 == length(cur.d1) | 0 == length(cur.d0)) {
      stop(paste("Stratum ", i, " contains no subjects from group 1 or 0", sep = ""));
    }
    
    cur.n1 <- length(cur.d1);
    cur.d  <- c(N0 = length(cur.d0), YBAR0 = mean(cur.d0), SD0   = sd(cur.d0),
                N1 = cur.n1,         YBAR1 = mean(cur.d1), YSUM1 = sum(cur.d1));
    
    stan.d <- rbind(stan.d, cur.d);
    Y1     <- c(Y1, cur.d1);
    INX1   <- c(INX1, rep(i, length = cur.n1));
  }
  
  if (is.null(RS))
    RS <- rep(1/S, S);
  
  lst.data  <- list(S     = S,
                    A     = A,
                    RS    = RS,
                    FIXVS = as.numeric(Fix.RS),
                    N0    = stan.d[,"N0"],
                    N1    = stan.d[,"N1"],
                    YBAR0 = stan.d[,"YBAR0"],
                    SD0   = stan.d[,"SD0"]);
  
  ## sampling
  if ("continuous" == type) {
   # stan.mdl  <-  "powerps";
    stan.mdl  <- ifelse(1 == S,  "powerp", "powerps");
    lst.data <- c(lst.data,
                  list(TN1  = length(Y1),
                       Y1   = Y1,
                       INX1 = INX1));
    rst.post  <- rweSTAN(lst.data = lst.data, stan.mdl = stan.mdl, ...);
    rst.theta <- rstan::extract(rst.post, pars = "theta")$theta;
  } else {
    lst.data <- c(lst.data,
                  list(YBAR1 = as.numeric(stan.d[,"YBAR1"]),
                       YSUM1 = as.numeric(stan.d[,"YSUM1"])));
    
    if (1 < S) {
      rst.post  <- rweSTAN(lst.data = lst.data, stan.mdl = "powerpsbinary", ...);
      rst.theta <- rstan::extract(rst.post, pars = "theta")$theta;
    } else {
      rst.post  <- NULL;
      rst.theta <- with(lst.data, rbeta(2000, YSUM1+A*YBAR0+1, N1-YSUM1+A*(1-YBAR0)+1));
    }
  }
  
  ## return
  list(post.theta = rst.theta,
       stan.rst   = rst.post);
}

#' Summary Posterior theta
#'
#'
#' @param post.theta posterior samples from STAN
#' @param true.theta true value of theta
#' @param quants     quantiles
#' @param weights    number of subjects in each strata as the weight
#'
#' @export
#'
rweSummaryPost <- function(post.theta, true.theta, quants = c(0.025, 0.975), weights=1) {
  
  ##cur.post      <- apply(post.theta, 2, function(x) {mean(weights*x)});
  cur.post        <- post.theta;
  cur.m           <- mean(cur.post);
  cur.var         <- var(cur.post);
  cur.ci          <- quantile(cur.post, quants);
  ##post.var        <- apply(post.theta, 1, var);
  ##names(post.var) <- paste("var", 1:length(post.var), sep = "");
  
  rst <- rweSummary(cur.m, cur.var, cur.ci, true.theta);
  rst
  
}


#' Get Prior Only for all stratum
#'
#' @param data class DWITHPS data frame
#' @param type type of outcomes
#' @param A    target number of subjects to be borrowed
#' @param RS   parameters for dirichelet prior
#' @param ...  extra parameters for calling function \code{\link{rweSTAN}}
#'
#' @export
#'
rwePsPowDrawPrior <- function(data, type = c("continuous", "binary"), A = 0, RS = NULL,
                              v.outcome = "Y",  iter = 2000, ...) {
  
  stopifnot(v.outcome %in% colnames(data));
  type <- match.arg(type);
  
  ## prepare stan data
  data   <- data[!is.na(data[["_strata_"]]),];
  S      <- max(data[["_strata_"]]);
  
  if (is.null(RS))
    RS <- rep(1/S, S);
  
  rst.theta <- NULL;
  for (i in 1:S) {
    cur.d0 <- data[data[["_strata_"]] == i & data[["_grp_"]] == 0, v.outcome];
    
    if (0 == length(cur.d0)) {
      warning(paste("Stratum ", i, " contains no subjects from group 0", sep = ""));
      next;
    }
    
    N0    <- length(cur.d0);
    YBAR0 <- mean(cur.d0);
    SD0   <- sd(cur.d0);
    ALPHA <- min(1, A*RS[i]/sum(RS)/N0);
    
    cat(A, ALPHA, N0, RS[i], "\n");
    
    ## sampling
    if ("continuous" == type) {
      rst.post  <- rweSTAN(lst.data = list(N0    = N0,
                                           YBAR0 = YBAR0,
                                           SD0   = SD0,
                                           ALPHA = ALPHA),
                           stan.mdl = "prior",
                           iter = iter,
                           ...);
      cur.theta <- rstan::extract(rst.post, pars = "theta")$theta;
    } else {
      cur.theta <- rbeta(iter, ALPHA*YBAR0*N0+1, ALPHA*(1-YBAR0)*N0+1);
    }
    
    rst.theta <- rbind(rst.theta, cbind(i, cur.theta));
  }
  
  ## return
  colnames(rst.theta) <- c("Stratum", "Samples");
  rst.theta           <- data.frame(rst.theta);
  rst.theta
}




#' Get propensity scores
#'
#' @param ... parameters to get propensity scores
#' @param d1.arm define strata based on a specific arm. Ignored if NULL.
#'
#' @export
#'
rwePS <- function(data, ps.fml = NULL,
                  v.grp   = "group",
                  v.arm   = "arm",
                  v.covs  = "V1",
                  d1.arm  = NULL,
                  d1.grp  = 1,
                  nstrata = 5, ...) {
  
  dnames <- colnames(data);
  stopifnot(v.grp %in% dnames);
  
  ## generate formula
  if (is.null(ps.fml))
    ps.fml <- as.formula(paste(v.grp, "~",
                               paste(v.covs, collapse = "+"),
                               sep = ""))
  
  ## d1 index will be kept in the results
  d1.inx   <- d1.grp == data[[v.grp]];
  keep.inx <- which(d1.inx);
  
  ## for 2-arm studies only
  if (!is.null(d1.arm))
    d1.inx <- d1.inx & d1.arm == data[[v.arm]];
  
  ## get ps
  all.ps  <- get_ps(data, ps.fml = ps.fml, ...);
  D1.ps   <- all.ps[which(d1.inx)];
  
  
  ## add columns to data
  grp     <- rep(1, nrow(data));
  grp[which(data[[v.grp]] != d1.grp)] <- 0;
  
  data[["_ps_"]]     <- all.ps;
  data[["_grp_"]]    <- grp;
  data[["_arm_"]]    <- data[[v.arm]];
  
  
  ## stratification
  if (nstrata > 0) {
    strata  <- rweCut(D1.ps, all.ps, breaks = nstrata, keep.inx = keep.inx);
    data[["_strata_"]] <- strata;
  }
  
  ## return
  rst <- list(data    = data,
              ps.fml  = ps.fml,
              nstrata = nstrata);
  class(rst) <- get.rwe.class("DWITHPS");
  
  rst
}

#' Get generalized propensity scores
#'
#' @param ... parameters to get propensity scores
#'
#' @export
#'
rweGPS <- function(data, ps.fml = NULL,
                   v.grp   = "group",
                   v.covs  = "V1",
                   nstrata = 5, ...) {
  
  ## likelihood
  f_ll <- function(beta) {
    xbeta <- apply(d_mat, 1,
                   function(x) sum(x * beta[-(1:nd1)]))
    
    exb <- sapply(xbeta,
                  function(x) {
                    tx <- beta[1:nd1] + x
                    tx <- 1 + sum(exp(tx))
                    log(tx)
                  })
    
    rst <- sum(beta[1:nd1] * n_j)
    rst <- rst + sum(xbeta[inx_j])
    rst <- rst - sum(exb)
    rst
  }
  
  f_gradient <- function(beta) {
    xbeta <- apply(d_mat, 1,
                   function(x) sum(x * beta[-(1:nd1)]))
    
    exb <- sapply(xbeta,
                  function(x) {
                    tx <- beta[1:nd1] + x
                    tx <- exp(tx)
                  })
    s_exb  <- apply(exb, 2, sum)
    s_1exb <- sapply(s_exb, function(x) {1 / (1+x)})
    
    g <- numeric(length(beta))
    for (i in 1:nd1) {
      g[i] <- n_j[i] - sum(s_1exb * exb[i,])
    }
    
    for (i in 1:nx) {
      g[i + nd1] <- sum(d_mat[inx_j, i]) - sum(s_1exb * s_exb * d_mat[, i])
    }
    
    g
  }
  
  dnames <- colnames(data);
  stopifnot(v.grp %in% dnames);
  
  nd  <- max(data[[v.grp]])
  nd1 <- nd - 1
  n_j <- NULL
  for (i in 2:nd) {
    n_j <- c(n_j, sum(data[[v.grp]] == i))
  }
  inx_j <- which(data[[v.grp]] > 1)
  
  ## design matrix
  if (is.null(ps.fml))
    ps.fml <- as.formula(paste(v.grp, "~", paste(v.covs, collapse = "+"), sep = ""))
  d_mat <- model.matrix(ps.fml, data)[, -1]
  
  ## mle
  ## mle0 <- optim(rep(0, nd1 + NC), f_ll,
  ##               method = "Nelder-Mead",
  ##               control = list(fnscale = -1, ...))
  mle <- optim(rep(0, nd1 + ncol(d_mat)),
               f_ll, gr = f_gradient,
               method = "BFGS", control = list(fnscale = -1, ...))
  
  ps_par <- mle$par
  
  ## gps and strata
  gps    <- apply(d_mat, 1,
                  function(x) sum(x * mle$par[-(1:nd1)]))
  
  strata <- rweCut(gps[which(1 == data[[v.grp]])], gps, breaks = nstrata)
  
  ## return
  data[["_gps_"]]    <- gps
  data[["_strata_"]] <- strata
  data[["_grp_"]]    <- data[[v.grp]]
  
  rst <- list(data    = data,
              ps.fml  = ps.fml,
              nd      = nd,
              nstrata = nstrata,
              mle     = mle)
  
  class(rst) <- get.rwe.class("D_GPS");
  rst
}


#' Get number of subjects and the distances of PS distributions for each PS
#' strata
#'
#' @param min.n0 threshold for N0, below which the external data in the
#'     current stratum will be ignored by setting the PS distance to 0
#'
#' @param d1.arm calculate distance based on a specific arm. Ignored if NULL.
#'
#' @export
#'
rwePSDist <- function(data.withps, n.bins = 10, min.n0 = 10,
                      type = c("ovl", "kl"), d1.arm = NULL, ...) {
  
  f.narm <- function(inx, dataps) {
    
    if (is.null(dataps[["_arm_"]]))
      return(c(length(inx), 0,0));
    
    n0 <- length(which(0 == dataps[inx, "_arm_"]));
    n1 <- length(which(1 == dataps[inx, "_arm_"]));
    
    c(length(inx), n0 ,n1);
  }
  
  stopifnot(inherits(data.withps,
                     what = get.rwe.class("DWITHPS")));
  
  type <- match.arg(type);
  
  dataps   <- data.withps$data;
  nstrata  <- data.withps$nstrata;
  rst      <- NULL;
  for (i in 1:nstrata) {
    
    inx.ps0 <- i == dataps[["_strata_"]] & 0 == dataps[["_grp_"]];
    inx.ps1 <- i == dataps[["_strata_"]] & 1 == dataps[["_grp_"]];
    n0.01   <- f.narm(which(inx.ps0), dataps);
    n1.01   <- f.narm(which(inx.ps1), dataps);
    
    if (!is.null(d1.arm) & !is.null(dataps[["_arm_"]])) {
      inx.ps0 <- inx.ps0 & d1.arm == dataps[["_arm_"]];
      inx.ps1 <- inx.ps1 & d1.arm == dataps[["_arm_"]];
    }
    
    ps0 <- dataps[which(inx.ps0), "_ps_"];
    ps1 <- dataps[which(inx.ps1), "_ps_"];
    
    if (0 == length(ps0) | 0 == length(ps1))
      warning("No samples in strata");
    
    if (any(is.na(c(ps0, ps1))))
      warning("NA found in propensity scores in a strata");
    
    if (length(ps0) < min.n0) {
      warning("Not enough data in the external data in the current stratum.
                     External data ignored.");
      cur.dist <- 0;
    } else {
      cur.dist <- rweDist(ps0, ps1, n.bins = n.bins, type = type, ...);
    }
    
    rst <- rbind(rst, c(i, n0.01, n1.01, cur.dist));
  }
  
  ## overall
  inx.tot.ps0 <- which(0 == dataps[["_grp_"]]);
  inx.tot.ps1 <- which(1 == dataps[["_grp_"]]);
  n0.tot.01   <- f.narm(inx.tot.ps0, dataps);
  n1.tot.01   <- f.narm(inx.tot.ps1, dataps);
  
  ps0         <- dataps[inx.tot.ps0, "_ps_"];
  ps1         <- dataps[inx.tot.ps1, "_ps_"];
  all.dist    <- rweDist(ps0, ps1, n.bins = nstrata*n.bins, type = type, ...);
  rst         <- rbind(rst, c(0, n0.tot.01, n1.tot.01, all.dist));
  
  
  colnames(rst) <- c("Strata", "N0", "N00", "N01", "N1", "N10", "N11", "Dist");
  rst           <- data.frame(rst);
  class(rst)    <- append(get.rwe.class("PSDIST"), class(rst));
  
  rst
}

#' Get number of subjects and the distances of PS distributions for each PS
#' strata for multiple data sources
#'
#' @param min.n0 threshold for N0, below which the external data in the
#'     current stratum will be ignored by setting the PS distance to 0
#'
#' @param d1.arm calculate distance based on a specific arm. Ignored if NULL.
#'
#' @export
#'
rweGpsDist <- function(data.gps, n.bins = 10, min.n0 = 10, type = c("ovl", "kl"),  ...) {
  
  stopifnot(inherits(data.gps,
                     what = get.rwe.class("D_GPS")));
  
  type     <- match.arg(type)
  dataps   <- data.gps$data
  nstrata  <- data.gps$nstrata
  nd       <- data.gps$nd
  
  rst      <- NULL
  for (i in 1:nstrata) {
    inx.ps1 <- i == dataps[["_strata_"]] & 1 == dataps[["_grp_"]]
    ps1     <- dataps[which(inx.ps1), "_gps_"];
    if (0 == length(ps1))
      warning(paste("No samples in strata", i, "in current study"))
    
    dist <- length(ps1)
    for (j in 2:nd) {
      inx.psj <- i == dataps[["_strata_"]] & j == dataps[["_grp_"]]
      psj     <- dataps[which(inx.psj), "_gps_"]
      
      if (0 == length(psj))
        warning(paste("No samples in strata", i, "in Study", j))
      
      if (length(psj) < min.n0) {
        warning("Not enough data in the external data in
                         the current stratum.")
        cur_dist <- 0
      } else {
        cur_dist <- rweDist(psj, ps1, n.bins = n.bins, type = type, ...)
      }
      
      dist <- c(dist, length(psj), cur_dist)
    }
    rst <- rbind(rst, c(i, dist))
  }
  
  ##overall
  inx.tot.ps1 <- which(1 == dataps[["_grp_"]]);
  ps1         <- dataps[inx.tot.ps1, "_gps_"];
  
  dist <- length(ps1)
  for (j in 2:nd) {
    inx.tot.psj <- which(j == dataps[["_grp_"]])
    psj         <- dataps[inx.tot.psj, "_gps_"]
    cur_dist    <- rweDist(psj, ps1, n.bins = nstrata * n.bins,
                           type = type, ...)
    dist        <- c(dist, length(psj), cur_dist)
  }
  rst <- rbind(rst, c(0, dist));
  
  
  colnames(rst) <- c("Strata", "N1",
                     paste(rep(c("N", "Dist"), nd - 1),
                           rep(2:nd, each = 2), sep = ""))
  rst           <- data.frame(rst);
  class(rst)    <- append(get.rwe.class("GPSDIST"), class(rst));
  
  rst
}


#' Get the actual power term in the power prior
#'
#' @param psdist A RWE_PSDIST type object
#' @param a power term
#' @param overall.ess ratio of overall added number of patients to N1
#' @param adjust.size whether adjust for sizes in group 0 and 1 in the power
#'     term
#' @param adjust.dist whether adjust for distance in ps scores in the power term
#'
#' @export
#'
rweGetPowerA <- function(psdist, a = NULL, overall.ess = 0.3,
                         adjust.size = TRUE, adjust.dist = TRUE) {
  
  stopifnot(inherits(psdist, what = get.rwe.class("PSDIST")));
  
  ## compute a
  if (is.null(a)) {
    stopifnot(1 == adjust.size);
    stopifnot(overall.ess >= 0);
    
    if (1 == adjust.dist) {
      a <- overall.ess / mean(psdist$Dist);
    } else {
      a <- overall.ess;
    }
  }
  
  
  ## compute as, power term for each strata
  rst <- rep(a, nrow(psdist));
  if (1 == adjust.size) {
    rst <- rst * psdist$N1 / psdist$N0;
  }
  
  if (1 == adjust.dist) {
    rst <- rst * psdist$Dist;
  }
  
  list(a  = a,
       as = rst);
}

#' Match on PS
#'
#' Match patients in RWD with patients in current study based on PS
#'
#' @param dta_cur current study data
#' @param dta_ext external data source data
#'
#' @export
#'
rwe_match_ps <- function(dta_cur, dta_ext, n_match = 3, ps.fml = NULL,
                         v.covs  = "V1") {
  
  n_cur <- nrow(dta_cur)
  n_ext <- nrow(dta_ext)
  ratio <- min(n_match, floor(n_ext / n_cur))
  
  dta_cur$grp_tmp <- 1
  dta_ext$grp_tmp <- 0
  dta_ext         <- dta_ext[, colnames(dta_cur)]
  
  dta    <- rbind(dta_cur, dta_ext)
  dta_ps <- rwePS(data = dta, v.grp = "grp_tmp", v.covs = v.covs, nstrata = 0)
  
  ps        <- dta_ps$data[["_ps_"]]
  target    <- ps[1:n_cur]
  candidate <- ps[-(1:n_cur)]
  
  ## match
  rst <- cMatch(target, candidate, ratio = ratio)[]
  rst <- rst[1:(ratio * n_cur)] + 1
  
  cbind(pid      = rep(1:n_cur, each = ratio),
        match_id = rst)
}

## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
##             PRIVATE FUNCTIONS
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------

## compute propensity scores
get_ps <- function(dta, ps.fml, type = c("logistic", "randomforest"),
                   ntree = 5000,
                   ..., grp = NULL, ps.cov = NULL) {
  
  type <- match.arg(type);
  
  ## generate formula
  if (is.null(ps.fml))
    ps.fml <- as.formula(paste(grp, "~", paste(ps.cov, collapse="+"),
                               sep=""));
  
  ## identify grp if passed from formula
  grp <- all.vars(ps.fml)[1];
  
  ## fit model
  switch(type,
         logistic = {
           glm.fit <- glm(ps.fml, family=binomial, data=dta, ...);
           est.ps <- glm.fit$fitted;
         },
         randomforest = {
           dta[[grp]] <- as.factor(dta[[grp]]);
           rf.fit     <- randomForest(ps.fml, data = dta,
                                      ntree = ntree, ...);
           est.ps     <- predict(rf.fit, type = "prob")[,2];
         });
  est.ps
}

#' Get unbalance in baseline covariates
#'
#' @param diff If TRUE, get the difference in covariates between groups.
#'     Otherwise, get covariates for each group separately
#'
#' @param var.group Column name in pts that corresponds to treat group
#'
#' @inheritParams simupara
#'
#' @return If diff is TRUE, return a dataset with columns V and Diff. Otherwise,
#'     return a dataset with columns V, Z and Value. In both cases, column V
#'     represents the name of covariates.
#'
#' @export
#'
rweUnbalance <- function(nPat, ..., pts = NULL, covs = NULL, diff = TRUE,
                         var.group = "Z", cov.pattern = "^V[0-9]+$") {
  if (is.null(pts)) {
    pts <- rweSimuTwoArm(nPat, ...);
  }
  
  ##unbalance
  inx.0 <- which(0 == pts[, var.group]);
  if (is.null(covs)) {
    c.xy  <- colnames(pts);
    c.xy  <- c.xy[grep(cov.pattern, c.xy)];
  } else {
    c.xy <- covs;
  }
  
  unb   <- NULL;
  for (i in 1:length(c.xy)) {
    x0     <- sample(pts[inx.0,  c.xy[i]], size = nPat, replace = TRUE);
    x1     <- sample(pts[-inx.0, c.xy[i]], size = nPat, replace = TRUE);
    
    if (diff) {
      x.diff <- x1 - x0;
      unb    <- rbind(unb, data.frame(V=c.xy[i], Diff=x1 - x0));
    } else {
      unb    <- rbind(unb, data.frame(V=c.xy[i], Z=1, Value=x1));
      unb    <- rbind(unb, data.frame(V=c.xy[i], Z=0, Value=x0));
    }
  }
  
  ## make group column factor if it exists
  if (!diff) {
    unb$Z <- as.factor(unb$Z)
  }
  
  unb
}


#' Get balance by different metrics
#'
#' @param cov0 covariates from group 0
#' @param cov1 covariates from group 1
#'
#' @export
#'
rweBalMetric <- function(cov0, cov1, metric = c("std", "abd", "ovl", "ksd",
                                                "ley", "mhb")) {
  metric <- match.arg(metric);
  switch(metric,
         std = {
           s <- sqrt((var(cov1) + var(cov0)) / 2)
           abs(mean(cov1) - mean(cov0)) / s;
         },
         abd = abs(mean(cov0) - mean(cov1)),
         ovl = metric.ovl(cov0, cov1),
         ksd = metric.ksd(cov0, cov1),
         ley = metric.ley(cov0, cov1),
         mhb = metric.mhb(cov0, cov1)
  )
}

## --------------------------------------------------------------------------------
## --------------------------------------------------------------------------------
##             PRIVATE FUNCTIONS
## --------------------------------------------------------------------------------
## --------------------------------------------------------------------------------

## overlapping coefficient
metric.ovl <- function(cov0, cov1) {
  cov <- c(cov0, cov1);
  if (length(unique(cov)) <= 10) {
    all.x <- c(rep(0, length(cov0)), rep(1, length(cov1)));
    pt    <- apply(prop.table(table(cov, all.x), 2), 1, min);
    
    ## reversed to measure imbalance
    return(1-sum(pt))
  }
  
  mn <- min(cov) * 1.25;
  mx <- max(cov) * 1.25;
  f1 <- approxfun(density(cov1, from = mn, to = mx, bw = "nrd"));
  f0 <- approxfun(density(cov0, from = mn, to = mx, bw = "nrd"));
  
  fn <- function(x)
    pmin(f1(x), f0(x))
  
  s <- try(integrate(fn, lower = mn, upper = mx,
                     subdivisions = 500)$value)
  ## Reverse: measure imbalance
  ifelse(inherits(s, "try-error"), NA, 1-s)
}

## K-S distance
metric.ksd <- function(cov0, cov1) {
  cov <- c(cov0, cov1);
  F1  <- ecdf(cov1);
  F0  <- ecdf(cov0);
  max(abs(F1(cov) - F0(cov)));
}

## Levy distance
metric.ley <- function(cov0, cov1) {
  cov <- c(cov0, cov1);
  F1  <- ecdf(cov1);
  F0  <- ecdf(cov0);
  e   <- max(abs(F1(cov) - F0(cov)));
  
  if (length(unique(cov)) <= 10)
    return(e)
  
  x     <- seq(min(cov), max(cov), length.out=1000);
  check <- all(F0(x-e) - e <= F1(x) & F1(x) <= F0(x+e) + e)
  
  while (check) {
    e <- e-.01
    check <- all(F0(x-e) - e <= F1(x) & F1(x) <= F0(x+e) + e);
  }
  
  e
}

#' mahalanobis balance
#'
#' covs should be a reduced datset that contains only those covariates
#' that will be used for calculating Mahalanobis balance, for example,
#' covs = dat[,1:6]
#' trt should be the exposure variable, for example, trt=dat$X
metric.mhb <- function(cov0, cov1) {
  S    <- rbind(cov0, cov1)
  Sinv <- solve(cov(S))
  x0   <- colMeans(S[1:length(cov0), ])
  x1   <- colMeans(S[-(1:length(cov0)), ])
  
  sum((t(x1 - x0) %*% Sinv) * (x1 - x0))
}



##-----------------------------------------------------------------------------
##     FUNCTIONS RELATED TO USING MARGINAL STATISTICS TO
##     1) SAMPLE PATIENTS
##     2) CALCULATE LIKELIHOOD
##     3)
##-----------------------------------------------------------------------------

#' Sample patient based on summary statistics
#'
#' @param target_stats target summary statistics
#' @param dta_ext external data
#' @param method sampling methods
#' @param n     number of patients to be drawn in each try
#' @param weights weights based on likelihood
#'
#' @export
#'
rwe_margin_sample <- function(target_stats, dta_ext,
                              method = c("genetic", "sa",
                                         "naive", "weighted", "ps"),
                              n_min = 300, max_steps = 10000, epsilon = NULL,
                              seed = NULL, ...) {
  
  method <- match.arg(method)
  
  if (!is.null(seed))
    set.seed(seed)
  
  rst <- tkCallFun(c("private_margin_sample_", method),
                   target_stats, dta_ext, n_min, max_steps, epsilon,
                   ...)
  
  ## return
  best_selected <- rst$best_selected
  best_val      <- rst$best_val
  
  list(selected = dta_ext[1 == best_selected, ],
       ext_ind  = best_selected,
       distance = best_val$distance,
       diff     = best_val$diff,
       stats    = best_val$stats,
       target   = target_stats)
}


#' Simulate covariates based on summary statistics
#'
#'
#' @return data frame with n patients, each column represents a covariate and
#'     each row represents a patient
#'
#' @export
#'
rwe_margin_simu <- function(target_stats, n = 500) {
  rst <- lapply(target_stats, private_margin_simu)
  rst <- data.frame(rst)
}

#' Calculate log-likelihood based on summary statistics
#'
#'
#' @param y matrix of data, each column represents a covariate and each row
#'     represents a patient
#' @return data frame for all y
#'
#' @export
#'
rwe_margin_ll <- function(target_stats, y) {
  rst <- NULL
  for (i in seq_len(length(target_stats))) {
    cur_stats <- target_stats[[i]]
    cur_y     <- y[[names(target_stats)[i]]]
    cur_ll    <- private_margin_ll(cur_stats, cur_y)
    rst       <- cbind(rst, cur_ll)
  }
  
  ll      <- apply(rst, 1, sum)
  weights <- exp(ll)
  weights <- weights / sum(weights)
  
  cbind(ll = ll, weights = weights)
}


#' Extract summary statistics
#'
#' Extract summary statistics from a data frame based on an existing summary
#' statistics list
#'
#' @return list of summary statistics
#'
#' @export
#'
rwe_extract_stats <- function(target_stats, y) {
  rst <- list()
  for (i in seq_len(length(target_stats))) {
    cur_v   <- names(target_stats)[i]
    cur_y   <- y[[cur_v]]
    cur_s   <- target_stats[[i]]
    
    cur_rst <- switch(
      cur_s$type,
      discrete   = tkExtractStats(cur_y, xlev = cur_s$values),
      continuous = tkExtractStats(cur_y, type = "continuous"),
      quants     = tkExtractStats(cur_y, type = "quants",
                                  quants = cur_s$quants[, 1]))
    rst[[cur_v]] <- cur_rst
  }
  rst
}

#' Calculate differences in summary statistics
#'
#' @return vector of differences in summary statistics
#'
#' @export
#'
rwe_margin_stat_diff <- function(target_stats, y, type = c("max", "sum"),
                                 epsilon = NULL, ...) {
  
  type        <- match.arg(type)
  target_stats_y <- rwe_extract_stats(target_stats, y)
  
  ## all differences
  rst <- NULL
  for (i in seq_len(length(target_stats))) {
    cur_diff <- private_stat_diff(target_stats[[i]],
                                  target_stats_y[[i]])
    rst      <- c(rst, cur_diff)
  }
  
  ## distance summary
  dist <- switch(type,
                 max = max(abs(rst)),
                 sum = sum(abs(rst)),
                 9999)
  
  ## acceptable or not
  yn <- FALSE
  if (!is.null(epsilon))
    yn <- dist < epsilon
  
  list(stats      = target_stats_y,
       diff       = rst,
       distance   = dist,
       acceptable = yn)
}



## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
##             PRIVATE FUNCTIONS
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------

## Simulate a covariate based on summary statistics
##
##
## @return column vector with n patients
##
private_margin_simu <- function(target_stats, n = 500) {
  ## sample based on quantiles
  fsmp_quants <- function(n, range, quants) {
    qmat <- private_quantile_mat(range, quants)
    rstq <- sample(seq_len(nrow(qmat)),
                   size = n, replace = TRUE, prob = qmat[, 1])
    rst  <- sapply(rstq, function(x) {
      runif(1, min = qmat[x, 2], max = qmat[x, 3])
    }, simplify = TRUE)
    
    rst
  }
  
  rst <- switch(target_stats$type,
                discrete = {
                  smps <- sample(seq_len(length(target_stats$values)),
                                 size    = n,
                                 replace = TRUE,
                                 prob    = target_stats$probs)
                  factor(target_stats$values[smps])
                },
                continuous = rnorm(n,
                                   target_stats$mean,
                                   target_stats$sd),
                quants = fsmp_quants(n,
                                     target_stats$range,
                                     target_stats$quants),
                NULL)
  rst
}

## Calculate log likelihood
##
## Calculate likelihood of a covariate based on its summary statistics
##
## @param y vector of data
## @param target_stats summary statistics
##
## @return column vector with n patients
##
private_margin_ll <- function(target_stats, y) {
  ## ll based proportions
  f_dis <- function(y, values, probs) {
    lp <- log(probs)
    sapply(y, function(x) {
      lp[which(x == values)]
    }, simplify = TRUE)
  }
  
  rst <- switch(target_stats$type,
                discrete   = f_dis(y, target_stats$values, target_stats$probs),
                continuous = dnorm(y, target_stats$mean, target_stats$sd,
                                   log = TRUE),
                NULL)
  rst
}

## Quantile matrix
##
## Get matrix data of quantiles
##
## @return matrix with 3 columns. 1: probabilities, 2: lower bound,
## 3: upper bound
##
private_quantile_mat <- function(range, quants) {
  qmat  <- NULL
  lastq <- 0
  lastx <- range[1]
  for (i in seq_len(nrow(quants))) {
    curq <- quants[i, 1]
    curx <- quants[i, 2]
    
    cmat  <- c(curq - lastq, lastx, curx)
    qmat  <- rbind(qmat, cmat)
    
    lastq <- curq
    lastx <- curx
  }
  
  cmat <- c(1 - lastq, lastx, range[2])
  qmat <- rbind(qmat, cmat)
}

## Difference in summary statistics
private_stat_diff <- function(stat0, stat1) {
  
  stopifnot(stat0$type == stat1$type)
  
  rst <- switch(stat0$type,
                discrete = {
                  stopifnot(identical(stat0$values,
                                      stat1$values))
                  diff <- stat0$probs - stat1$probs
                },
                continuous = {
                  c(stat0$mean - stat1$mean,
                    stat0$sd   - stat1$sd)
                },
                quants = {
                  stopifnot(identical(stat0$quants[, 1],
                                      stat1$quants[, 1]))
                  stat0$quants[, 2] - stat1$quant[, 2]
                },
                ... = NULL)
}

## A tool function for getting distance
f_dist <- function(target_stats, dta_ext, ind_selected, ...) {
  cur_sel <- which(1 == ind_selected)
  
  if (0 == length(cur_sel))
    return(NULL)
  
  rst_dist <- rwe_margin_stat_diff(target_stats,
                                   y = dta_ext[cur_sel, ],
                                   ...)
  rst_dist
}

## Naive Sampling from RWD
private_margin_sample_naive <- function(target_stats, dta_ext, n_min = 300,
                                        max_steps = 10000, ...,
                                        weighted = FALSE, ps = FALSE) {
  
  f_select <- function() {
    cur_n   <- sample(n_min:n_max, 1)
    
    if (!ps) {
      cur_inx <- sample(1:n_max, cur_n, replace = FALSE,
                        prob = weights)
    } else {
      simu_target <- rwe_margin_simu(n = cur_n, target_stats)
      cur_inx     <- rwe_match_ps(simu_target,
                                  dta_ext,
                                  n_match = 1,
                                  v.covs  = names(target_stats))
    }
    
    cur_selected           <- rep(0, n_max)
    cur_selected[cur_inx]  <- 1
    cur_val                <- f_dist(target_stats, dta_ext, cur_selected,
                                     ...)
    
    list(selected = cur_selected,
         val      = cur_val)
  }
  
  ## total patients
  n_max <- nrow(dta_ext)
  
  ## set weights
  if (weighted) {
    weights <- rwe_margin_ll(target_stats, dta_ext)[, "weights"]
  } else {
    weights <- rep(1/n_max, n_max)
  }
  
  ## initiate
  cur_rst       <- f_select()
  best_selected <- cur_rst$selected
  best_val      <- cur_rst$val
  
  k <- 1
  while (k < max_steps & !best_val$acceptable) {
    cur_rst <- f_select()
    
    if (best_val$distance > cur_rst$val$distance) {
      best_selected <- cur_rst$selected
      best_val      <- cur_rst$val
      
      ## cat("-----", k, sum(best_selected), best_val$distance, "\n")
    }
    
    k <- k + 1
  }
  
  list(best_selected = best_selected,
       best_val      = best_val)
}

## Weighted Sampling from RWD
private_margin_sample_weighted <- function(...) {
  private_margin_sample_naive(..., weighted = TRUE, ps = FALSE)
}


## PS Sampling from RWD
private_margin_sample_ps <- function(...) {
  private_margin_sample_naive(..., weighted = FALSE, ps = TRUE)
}

## Simulated Annealing for Sampling from RWD
private_margin_sample_sa <- function(target_stats, dta_ext, n_min = 300,
                                     max_steps = 10000, ..., alpha = 0.1,
                                     weighted = TRUE, swap_size = 5,
                                     p_extend = 0.5) {
  
  ## select samples from available or selected pool
  f_swap <- function(ind_selected) {
    ## extend or not
    to_extend <- rbinom(1, 1, p_extend)
    
    if (0 == to_extend | all(1 == ind_selected)) {
      ## remove current selected
      inx <- which(1 == ind_selected)
      rst <- sample(inx, size = min(swap_size, length(inx) - n_min),
                    prob = 1 - weights[inx])
      ind_selected[rst] <- 0
    } else {
      ## add to selected patients
      inx <- which(0 == ind_selected)
      rst <- sample(inx, size = min(swap_size, length(inx)),
                    prob = weights[inx])
      ind_selected[rst] <- 1
    }
    ind_selected
    
  }
  
  ## total patients
  n_max <- nrow(dta_ext)
  
  ## set weights
  if (weighted) {
    weights <- rwe_margin_ll(target_stats, dta_ext)[, "weights"]
  } else {
    weights <- rep(1/nrow(dta_ext), nrow(dta_ext))
  }
  
  ## indicator of selected patients
  cur_n                 <- sample(n_min:n_max, 1)
  cur_inx               <- sample(1:n_max, cur_n, replace = FALSE,
                                  prob = weights)
  cur_selected          <- rep(0, n_max)
  cur_selected[cur_inx] <- 1
  cur_val               <- f_dist(target_stats, dta_ext, cur_selected, ...)
  
  ## initiate best values
  best_selected <- cur_selected
  best_val      <- cur_val
  
  k <- 1
  while (k < max_steps & !best_val$acceptable) {
    temp_T        <- alpha * k / max_steps
    next_selected <- f_swap(cur_selected)
    next_val      <- f_dist(target_stats, dta_ext, next_selected, ...)
    
    diff_val      <- next_val$distance - cur_val$distance
    p_update      <- exp( - diff_val / temp_T)
    
    if (p_update >= runif(1)) {
      cur_selected <- next_selected
      cur_val      <- next_val
    }
    
    if (best_val$distance > cur_val$distance) {
      best_selected <- cur_selected
      best_val      <- cur_val
      
      ## cat("-----", k, sum(best_selected), best_val$distance, "\n")
    }
    
    k <- k + 1
  }
  
  list(best_selected = best_selected,
       best_val      = best_val)
}

## Genetic algorithm for sampling
private_margin_sample_genetic <- function(target_stats, dta_ext, n_min = 300,
                                          max_steps = 10000, ..., monitor = FALSE) {
  
  fitness <- function(ind_selected, ...) {
    rst <- f_dist(target_stats, dta_ext, ind_selected, ...)
    return(-rst$distance)
  }
  
  rst <- ga(type = "binary",
            fitness = fitness, nBits = nrow(dta_ext),
            keepBest = TRUE, maxiter = max_steps,
            monitor = monitor, ...)
  
  best_selected <- as.numeric(rst@solution)
  best_val      <- f_dist(target_stats, dta_ext, best_selected, ...)
  
  list(best_selected = best_selected,
       best_val      = best_val)
}

## Random forest for Sampling from RWD
private_margin_sample_rf <- function(target_stats, dta_ext, n_min = 300,
                                     max_steps = 10000, ...) {
  
  ## current distance
  f_dis <- function(ind_selected) {
    cur_sel  <- which(1 == ind_selected)
    rst_dist <- rwe_margin_stat_diff(target_stats,
                                     y = dta_ext[cur_sel, ],
                                     type = type)
    rst_dist
  }
  
  ## initial random seed
  if (!is.null(seed))
    set.seed(seed)
  
  ## set weights
  if (weighted) {
    weights <- rwe_margin_ll(target_stats, dta_ext)[, "weights"]
  } else {
    weights <- rep(1/nrow(dta_ext), nrow(dta_ext))
  }
  
  ## indicator of selected patients
  cur_selected <- rep(1, nrow(dta_ext))
  cur_val      <- f_dis(cur_selected)
  best_selected <- cur_selected
  best_val      <- cur_val
  
  ind_tree <- rep(1, nrow(dta_ext))
  k        <- 1
  while (sum(cur_selected) > n_min
         & best_val$distance > epsilon
         & k < max_steps) {
    cur_size <- sample(1:swap_size, 1)
    inx_cand <- which(1 == cur_selected)
    cur_inx  <- sample(inx_cand,
                       size = cur_size,
                       prob = 1 - weights[inx_cand])
    next_selected          <- cur_selected
    next_selected[cur_inx] <- 0
    next_val               <- f_dis(next_selected)
    diff_val               <- next_val$distance - cur_val$distance
    
    if (diff_val < 0) {
      cur_selected <- next_selected
      cur_val      <- next_val
      
      if (best_val$distance > cur_val$distance) {
        best_selected <- cur_selected
        best_val      <- cur_val
        
        cat("iter", k, "distance", best_val$distance, "\n")
      }
    }
    
    ## keep node
    ind_tree[cur_inx] <- 0
    ## update steps
    k                 <- k + 1
  }
}



#' Plot unbalance in covariates
#'
#'
#' @export
#'
rwePlotUnbalance <- function(data.unb,
                             var.x     = "Diff",
                             var.group = NULL,
                             xlim      = NULL,
                             ylim      = NULL,
                             title     = "",
                             f.grid    = formula("V~Study"),
                             adjust    = 1) {
  
  if (is.null(var.group)) {
    rst <- ggplot(data.unb, aes_string(x=var.x)) +
      geom_density(alpha = 0.4, fill = "gray", na.rm = TRUE, adjust = adjust) +
      geom_vline(xintercept = 0, linetype="dashed", col = "red");
  } else {
    rst <- ggplot(data.unb,
                  aes_string(x     = var.x,
                             group = var.group,
                             color = var.group,
                             linetype = var.group)) +
      geom_density(alpha = 0.2, fill = "gray", na.rm = TRUE, adjust = adjust, );
  }
  
  rst <- rst + labs(x = "", y="", title=title) +
    theme_bw() +
    theme(strip.background = element_blank(),
          panel.grid       = element_blank(),
          panel.border     = element_blank(),
          axis.text.y      = element_blank(),
          axis.text.x      = element_text(size=9),
          axis.ticks.y     = element_blank(),
          panel.spacing    = unit(0, "lines"))+
    facet_grid(f.grid);
  
  if (!is.null(xlim))
    rst <- rst + scale_x_continuous(limits = xlim, breaks = 0);
  if (!is.null(ylim))
    rst <- rst + scale_y_continuous(limits = ylim);
  
  rst
}


#' Plot PS distributions
#'
#'
#' @method plot RWE_DWITHPS
#'
#' @export
#'
plot.RWE_DWITHPS <- function(x, type = c("ps", "balance"), ...) {
  type <- match.arg(type);
  
  switch(type,
         ps = plotRwePs(x, ...),
         balance = plotRweBalance(x, ...));
}


#' add a method summary2
#'
#'
#' @export
#'
summary2 <- function (x, ...) {
  UseMethod("summary2", x)
}

#' Summary unbalance metrics
#'
#'
#' @export
#'
#' @method summary2 RWE_DWITHPS
#'
summary2.RWE_DWITHPS <- function(x, v.cov = NULL, label.cov = v.cov, ...) {
  if (is.null(v.cov))
    v.cov <- all.vars(x$ps.fml)[-1];
  
  if (is.null(label.cov))
    label.cov <- v.cov;
  
  nstrata      <- x$nstrata;
  dtaps        <- x$data;
  dtaps$Strata <- dtaps[["_strata_"]];
  dtaps$Group  <- dtaps[["_grp_"]];
  
  rst <- NULL;
  for (iv in 1:length(v.cov)) {
    v        <- v.cov[iv];
    cur.cov  <- dtaps[[v]];
    if (is.factor(cur.cov)) {
      cur.cov <- as.numeric(cur.cov);
    }
    
    for (stra in 1:max(dtaps$Strata, na.rm=TRUE)) {
      cur.cov0 <- cur.cov[which(0 == dtaps$Group & stra == dtaps$Strata)];
      cur.cov1 <- cur.cov[which(1 == dtaps$Group & stra == dtaps$Strata)];
      
      cur.met  <- rweBalMetric(cur.cov0, cur.cov1, ...);
      rst      <- rbind(rst,
                        data.frame(Covariate = label.cov[iv],
                                   Strata    = stra,
                                   Metric    = cur.met)
      );
    }
    
    ##overall
    cur.cov0 <- cur.cov[which(0 == dtaps$Group & !is.na(dtaps$Strata))];
    cur.cov1 <- cur.cov[which(1 == dtaps$Group & !is.na(dtaps$Strata))];
    cur.met  <- rweBalMetric(cur.cov0, cur.cov1, ...);
    rst      <- rbind(rst,
                      data.frame(Covariate = label.cov[iv],
                                 Strata    = 0,
                                 Metric    = cur.met)
    );
    
  }
  
  rst
}

## --------------------------------------------------------------------------------
## --------------------------------------------------------------------------------
##             PRIVATE FUNCTIONS
## --------------------------------------------------------------------------------
## --------------------------------------------------------------------------------

## plot density of propensity score
plotRwePs <- function(data.withps, overall.inc = TRUE, add.text = TRUE,
                      facet.scales = "free_y", ...) {
  stopifnot(inherits(data.withps,
                     what = get.rwe.class("DWITHPS")));
  
  pskl        <- rwePSDist(data.withps, ...);
  nstrata     <- data.withps$nstrata;
  dtaps       <- data.withps$data;
  
  pskl$Strata <- as.factor(c(paste("Stratum ", 1:nstrata, sep = ""), "Overall"));
  xlim        <- range(dtaps[which(!is.na(dtaps[["_strata_"]])),"_ps_"], na.rm = TRUE);
  
  all.data <- NULL;
  for (i in 1:nstrata) {
    cur.sub  <- dtaps[which(i == dtaps[["_strata_"]]),];
    cur.data <- data.frame(Strata = paste("Stratum ", i, sep = ""),
                           Ps     = cur.sub[["_ps_"]],
                           Group  = cur.sub[["_grp_"]]);
    all.data <- rbind(all.data, cur.data);
  }
  
  if (overall.inc) {
    cur.data <-  data.frame(Strata = "Overall",
                            Ps     = dtaps[["_ps_"]],
                            Group  = dtaps[["_grp_"]]);
    
    all.data <- rbind(all.data, cur.data);
  } else {
    pskl <- pskl %>% filter(Strata != "Overall");
  }
  
  all.data$Group <- as.factor(all.data$Group);
  rst <- ggplot(data = all.data, aes(x = Ps)) +
    geom_density(alpha = 0.2,
                 aes(group = Group,
                     fill  = Group,
                     linetype = Group),
                 trim  = TRUE,
                 na.rm = TRUE) +
    labs(x = "Propensity Score", y = "Density") +
    scale_y_continuous(breaks = NULL) +
    scale_x_continuous(limits = xlim) +
    scale_fill_manual(values=c("gray20", "gray80")) +
    theme_bw() +
    theme(strip.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(colour = "black"),
          panel.spacing = unit(0, "lines")) +
    facet_grid(Strata ~ ., scales = facet.scales);
  
  if (add.text) {
    rst <- rst +
      geom_text(x = Inf, y = Inf, hjust = 1, vjust = 1,
                aes(label = paste('n0=', N0,
                                  ", n1=", N1,
                                  ", OVL=", format(Dist, digits = 3),
                                  sep = "")),
                data = pskl, size = 4);
  }
  rst
}

plot.balance.fac <- function(dtaps, v, overall.inc = TRUE) {
  cur.d <- rweFreqTbl(dtaps, var.groupby = c("Strata", "Group"), vars = v);
  cur.d <- cur.d %>% dplyr::filter(!is.na(Strata));
  cur.d$Strata <- paste("Stratum ", cur.d$Strata, sep = "")
  if (overall.inc) {
    cur.overall <- rweFreqTbl(dtaps, var.groupby = "Group", vars = v);
    cur.overall$Strata <- "Overall";
    cur.d <- rbind(cur.d, cur.overall);
  }
  cur.d$Group <- as.factor(cur.d$Group);
  cur.d$Value <- as.factor(cur.d$Value);
  
  rst <- ggplot(data = cur.d, aes(x = Value, y = Freq)) +
    geom_bar(alpha = 0.4,
             stat = "identity",
             position = "dodge",
             color = "black",
             aes(group = Group,
                 fill  = Group)) +
    scale_fill_manual(values=c("gray20", "gray80")) +
    scale_y_continuous(breaks = NULL, limits = c(0,1)) +
    labs(x = "", y = "") +
    facet_grid(Strata ~ .);
  rst
}

plot.balance.cont <- function(dtaps, v, nstrata, overall.inc = TRUE, facet.scales = "free_y") {
  cur.d <- NULL;
  for (i in 1:nstrata) {
    cur.sub      <- dtaps[which(i == dtaps[["_strata_"]]),];
    cur.v        <- data.frame(Cov    = v,
                               Value  = cur.sub[[v]],
                               Group  = cur.sub[["_grp_"]]);
    cur.v$Strata <- paste("Stratum ", i, sep = "");
    cur.d        <- rbind(cur.d, cur.v);
  }
  
  if (overall.inc) {
    cur.sub      <- dtaps;
    cur.v        <- data.frame(Cov    = v,
                               Value  = cur.sub[[v]],
                               Group  = cur.sub[["_grp_"]]);
    cur.v$Strata <- paste("Overall");
    cur.d        <- rbind(cur.d, cur.v);
  }
  cur.d$Group <- as.factor(cur.d$Group);
  
  rst <- ggplot(data = cur.d, aes(x = Value)) +
    geom_density(alpha = 0.2,
                 aes(group = Group,
                     fill  = Group,
                     linetype = Group),
                 na.rm = TRUE) +
    scale_y_continuous(breaks = NULL) +
    scale_fill_manual(values=c("gray20", "white")) +
    labs(x = "", y = "") +
    facet_grid(Strata ~ ., scales = facet.scales);
  rst
}


plotRweBalance <- function(data.withps, overall.inc = TRUE, v.cov = NULL,
                           facet.scales = "free_y", label.cov = v.cov, legend.width = 0.08,
                           ...) {
  
  if (is.null(v.cov))
    v.cov <- all.vars(data.withps$ps.fml)[-1];
  
  if (is.null(label.cov))
    label.cov <- v.cov;
  
  nstrata      <- data.withps$nstrata;
  dtaps        <- data.withps$data;
  dtaps$Strata <- dtaps[["_strata_"]];
  dtaps$Group  <- dtaps[["_grp_"]];
  
  rst <- list();
  for (v in v.cov) {
    if (is.factor(dtaps[[v]])) {
      cur.p <- plot.balance.fac(dtaps, v, overall.inc = overall.inc);
    } else {
      cur.p <- plot.balance.cont(dtaps, v, nstrata = nstrata,
                                 overall.inc = overall.inc, facet.scales = facet.scales);
    }
    cur.p <- cur.p +
      labs(title = label.cov[v == v.cov]) +
      theme_bw() +
      theme(strip.background = element_blank(),
            strip.placement  = "right",
            strip.text       = element_blank(),
            panel.grid       = element_blank(),
            panel.border     = element_blank(),
            panel.spacing    = unit(0, "lines"),
            plot.title = element_text(hjust = 0.5),
            legend.position  = "none",
            plot.margin      = unit(c(1,0,1,-0.5), "lines"));
    
    rst[[v]] <- cur.p;
  }
  
  rst[[length(rst)]] <- rst[[length(rst)]] +
    theme(strip.text = element_text(size=8),
          legend.position = "right")
  
  rst$nrow        <- 1;
  rst$rel_widths <- c(rep(1, length(v.cov)-1),
                      1+legend.width*length(v.cov));
  do.call(plot_grid, rst);
}




#' Simulate covariates following a mixture of multivariate normal distribution
#'
#' @inheritParams simupara
#'
#' @export
#'
rweSimuCov <- function(nPat, muCov, sdCov, corCov, mix.phi = 1,
                       seed = NULL,
                       cov.breaks = NULL) {
  
  f.cur <- function(x, i) {
    if (is.array(x)) {
      rst <- x[min(i, nrow(x)),];
    } else {
      rst <- x;
    }
    
    rst
  }
  
  stopifnot(is.numeric(mix.phi) | any(mix.phi < 0));
  stopifnot(nPat > 0);
  
  if (!is.null(seed))
    set.seed(seed);
  
  n.pts <- rmultinom(1, nPat, mix.phi);
  cov.x <- NULL;
  for (i in 1:length(mix.phi)) {
    
    if (0 == n.pts[i])
      next;
    
    cur.mu  <- f.cur(muCov, i);
    cur.sd  <- f.cur(sdCov, i);
    cur.cor <- corCov[min(i, length(corCov))];
    cur.x   <- rmvnorm(n.pts[i],
                       mean  = cur.mu,
                       sigma = get.covmat(cur.sd, cur.cor));
    cov.x   <- rbind(cov.x, cur.x);
  }
  
  colnames(cov.x) <- paste("V", 1:ncol(cov.x), sep = "");
  cov.x           <- data.frame(cov.x);
  cov.x           <- get_cov_cat(cov.x, cov.breaks);
}

#' Simulate X multiplied by Beta
#'
#' @inheritParams simupara
#' @param ... Parameters for simulating covariates by function
#'     \code{\link{rweSimuCov}}
#'
#'
#' @export
#'
rweXBeta <- function(..., regCoeff, cov.x = NULL, fmla = NULL) {
  
  stopifnot(inherits(fmla, "formula") | is.null(fmla));
  
  if (is.null(cov.x))
    cov.x <- rweSimuCov(...);
  
  if (is.null(fmla)) {
    fmla <- formula(paste("~",
                          paste(colnames(cov.x), collapse = "+")));
  }
  
  d.matrix <- model.matrix(fmla, cov.x);
  xbeta    <- get.xbeta(d.matrix, regCoeff);
  xbeta
}


#' Compute standard error of the random error
#'
#' @inheritParams simupara
#' @inheritParams rweGetYSig
#'
#' @return mean of xbeta and standard error or the random error term
#'
#' @export
#'
rweGetYSig <- function(..., nPat=500000, xbeta = NULL, sig2Ratio = 1) {
  if (is.null(xbeta))
    xbeta   <- rweXBeta(nPat=nPat, ...);
  
  v.xbeta <- var(xbeta);
  ysig    <- sqrt(v.xbeta * sig2Ratio);
  
  c(mean(xbeta), ysig);
}

#' Get intercept for a binary outcome.
#'
#' The binary outcome may be an outcome or a treatment assignment.
#'
#' @inheritParams simupara
#'
#' @param ... Parameters for simulating covariates by function
#'     \code{\link{rweXBeta}} and \code{\link{rweSimuCov}}
#'
#' @return standard error or the random error term
#'
#' @export
#'
rweGetBinInt <- function(..., regCoeff, nPat=500000, xbeta = NULL, bin.mu = 0.5) {
  ## fill in 0 for intercept temporarily
  if (is.null(xbeta))
    ey <- rweXBeta(nPat, regCoeff = c(0, regCoeff), ...);
  
  fx <- function(b0, bmu) {
    logp <- (b0 + ey) - log(1 + exp(b0+ey));
    m    <- mean(exp(logp));
    rst  <- abs(m - bmu);
  }
  
  mey <- max(abs(ey));
  rst <- NULL;
  
  for (i in 1:length(bin.mu)) {
    cur.rst <- optimize(fx, c(-100 - mey, 100 + mey),
                        bmu = bin.mu[i])$minimum
    rst     <- c(rst, cur.rst);
  }
  
  rst
}

#' Simulate random errors
#'
#' @inheritParams simupara
#'
#' @export
#'
rweSimuError <- function(nPat,
                         error.type = c("normal", "skewed"),
                         ysig   = 1,
                         skew.n = NULL,skew.p = NULL, skew.noise = 0.0001,
                         ...) {
  
  type <- match.arg(error.type);
  rst <- switch(type,
                normal = {rnorm(nPat, 0, ysig)},
                skewed = {
                  mu        <- skew.n * (1-skew.p) / skew.p;
                  va        <- skew.n * (1-skew.p) / skew.p^2;
                  noise.sig <- skew.noise;
                  rst       <- rnbinom(nPat, skew.n, skew.p);
                  rst       <- rst - mu + rnorm(nPat, 0, noise.sig);
                  rst       <- rst/sqrt(va + noise.sig^2)*ysig;
                });
  
  rst
}

#' Simulate outcomes for a single arm study
#'
#' @inheritParams simupara
#'
#' @export
#'
rweSimuSingleArm <- function(nPat, muCov, sdCov, corCov, regCoeff, mix.phi = 1,
                             cov.breaks = NULL,
                             fmla = NULL,
                             type = c("continuous", "binary"),
                             ysig = NULL, sig2Ratio=1, b0 = NULL, bin.mu = 0.5,
                             ...) {
  
  type  <- match.arg(type);
  COV.X <- rweSimuCov(nPat       = nPat,
                      muCov      = muCov,
                      sdCov      = sdCov,
                      corCov     = corCov,
                      mix.phi    = mix.phi,
                      cov.breaks = cov.breaks);
  ##simulate Y
  if ("continuous" == type) {
    ## epsilon
    if (is.null(ysig)) {
      ysig <- rweGetYSig(muCov      = muCov,
                         sdCov      = sdCov,
                         corCov     = corCov,
                         mix.phi    = mix.phi,
                         cov.breaks = cov.breaks,
                         regCoeff   = regCoeff,
                         sig2Ratio  = sig2Ratio,
                         fmla       = fmla)[2];
    }
    
    XBETA   <- rweXBeta(cov.x = COV.X, regCoeff = regCoeff, fmla = fmla);
    EPSILON <- rweSimuError(nPat, ysig = ysig,...);
    Y       <- XBETA + EPSILON;
  } else if ("binary" == type) {
    if (is.null(b0)) {
      stopifnot(!is.null(bin.mu));
      b0  <- rweGetBinInt(bin.mu,
                          muCov      = muCov,
                          sdCov      = sdCov,
                          corCov     = corCov,
                          regCoeff   = regCoeff,
                          mix.phi    = mix.phi,
                          cov.breaks = cov.breaks,
                          fmla       = fmla);
    }
    regCoeff <- c(b0, regCoeff);
    XBETA    <- rweXBeta(cov.x = COV.X, regCoeff = regCoeff, fmla = fmla);
    Y        <- rbinom(nPat, 1, expit(XBETA));
  }
  
  ##return
  Data           <- cbind(1:nPat, Y, COV.X);
  colnames(Data) <- c("pid", "Y", paste("V", 1:ncol(COV.X), sep=""));
  data.frame(Data);
}

#' Simulate continuous outcomes for a two arm study
#'
#' @inheritParams simupara
#'
#' @export
#'
rweSimuTwoArm <- function(nPat, muCov, sdCov, corCov, trt.effect = 0,
                          regCoeff.y, regCoeff.z=0,
                          mix.phi = 1, cov.breaks = NULL,
                          fmla.y = NULL, fmla.z = NULL, ysig = NULL,
                          b0 = NULL, z1.p = 0.5,
                          sig2Ratio = 2,  ..., do.simu=TRUE) {
  
  ##treatment assignment
  if (is.null(b0)) {
    b0  <- rweGetBinInt(bin.mu   = z1.p,
                        muCov    = muCov,
                        sdCov    = sdCov,
                        corCov   = corCov,
                        regCoeff = regCoeff.z,
                        mix.phi  = mix.phi,
                        fmla     = fmla.z);
  }
  
  if (is.null(ysig)) {
    ysig <- rweGetYSig(muCov      = muCov,
                       sdCov      = sdCov,
                       corCov     = corCov,
                       mix.phi    = mix.phi,
                       cov.breaks = cov.breaks,
                       regCoeff   = regCoeff.y,
                       sig2Ratio  = sig2Ratio,
                       fmla       = fmla.y)[2];
  }
  
  simu.data <- NULL;
  if (do.simu) {
    ##covariates
    COV.X   <- rweSimuCov(nPat       = nPat,
                          muCov      = muCov,
                          sdCov      = sdCov,
                          corCov     = corCov,
                          mix.phi    = mix.phi,
                          cov.breaks = cov.breaks);
    
    if (identical(0, regCoeff.z)) {
      Z  <- rbinom(nPat, 1, z1.p);
    } else {
      xbeta.z <- rweXBeta(cov.x    = COV.X,
                          regCoeff = regCoeff.z,
                          fmla     = fmla.z);
      Z       <- rbinom(nPat, 1, expit(b0 + xbeta.z));
    }
    
    xbeta.y <- rweXBeta(cov.x    = COV.X,
                        regCoeff = regCoeff.y,
                        fmla     = fmla.y);
    
    epsilon.y <- rweSimuError(nPat, ysig = ysig, ...);
    Y         <- Z * trt.effect + xbeta.y + epsilon.y;
    simu.data <- cbind(pid=1:nPat, Y=Y, Z=Z, COV.X);
  }
  
  list(true.effect = trt.effect,
       simu.data   = simu.data,
       b0ysig      = c(b0 = b0, ysig = ysig));
}

#' Simulate data from an existing dataset
#'
#'
#' @inheritParams simupara
#'
#' @export
#'
rweSimuFromTrial <- function(nPat, trial.data, group = "A", outcome = "Y",
                             with.replacement = TRUE, seed = NULL,
                             permute = TRUE, f.subset = NULL,
                             permute.trteffect = 0, permute.interaction = 0,
                             simu.group = group, simu.outcome = outcome) {
  if (!is.null(seed))
    set.seed(seed);
  
  if (1 == length(nPat)) {
    ## set the same sample size for the two groups
    nPat <- rep(nPat,2);
  }
  
  if (!permute) {
    arms   <- unique(trial.data[[group]]);
    rst    <- NULL;
    mean.y <- NULL;
    for (i in 1:length(arms)) {
      cur.d   <- trial.data[arms[i] == trial.data[[group]],];
      cur.n   <- nPat[min(i, length(nPat))];
      
      stopifnot(with.replacement | nrow(cur.d) > cur.n);
      
      mean.y   <- c(mean.y, mean(cur.d[[outcome]]));
      cur.smp  <- cur.d[sample(1:nrow(cur.d), cur.n, replace=with.replacement), ];
      
      cur.smp[[simu.group]]   <- i - 1;
      cur.smp[[simu.outcome]] <- cur.d[[outcome]];
      
      rst <- rbind(rst, cur.smp);
    }
    
    trt.effect <- mean.y;
    simu.data  <- rst;
  } else {
    cur.d <- trial.data;
    cur.n <- sum(nPat);
    
    stopifnot(with.replacement | nrow(cur.d) > cur.n);
    
    smp.inx <- sample(1:nrow(cur.d), cur.n, replace=with.replacement);
    cur.smp <- cur.d[smp.inx, ];
    grps    <- NULL;
    for (i in 1:length(nPat)) {
      grps <- c(grps, rep(i-1, nPat[i]));
    }
    cur.smp[[simu.group]]   <- grps;
    cur.smp[[simu.outcome]] <- cur.smp[[outcome]];
    
    ##introduce main effect to the last group
    inx.last <- which(max(grps) == grps);
    cur.smp[inx.last, simu.outcome] <- permute.trteffect + cur.smp[inx.last, simu.outcome];
    
    ##introduce interaction effect
    if (is.function(f.subset) & permute.interaction != 0) {
      all.subgrp <- f.subset(cur.d);
      stopifnot(all(all.subgrp %in% c(0,1)));
      
      egbi   <- mean(all.subgrp);
      subgrp <- all.subgrp[smp.inx];
      cur.smp[inx.last, simu.outcome] <- cur.smp[inx.last, simu.outcome] +
        permute.interaction * (subgrp[inx.last] - egbi);
    }
    
    trt.effect <- permute.trteffect;
    simu.data  <- cur.smp;
  }
  
  ## randomize the order of simulated patients
  list(true.effect = trt.effect,
       simu.data   = simu.data[sample(1:nrow(simu.data)),]);
}


