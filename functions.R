###################################################################################################
######################### Script with additional functions ########################################
###################################################################################################

#' Gaussian mixture modelling for the verification of reference intervals
#'
#' @param x Numeric vector of positive laboratory values without missing values.
#'   The length of \code{x} should be at least 200 observations.
#' @param lognormal Logical; if \code{TRUE}, a lognormal transformation is applied
#'   before clustering. Default is \code{FALSE}.
#' @param targets Optional numeric vector specifying target values. Default is \code{NULL}.
#' @param plot.it Logical; if \code{TRUE}, diagnostic plots are generated.
#'   Default is \code{TRUE}.
#' @param add.boxplot Logical; if \code{TRUE}, a boxplot is added to the histogram.
#'   Default is \code{TRUE}.
#' @param plot.legend Logical; if \code{TRUE}, a legend is added to the plot.
#'   Default is \code{TRUE}.
#' @param pos.legend Character string specifying the legend position, passed to
#'   \code{legend()}. Default is \code{"topright"}.
#' @param plot.bic Logical; if \code{TRUE}, the Bayesian Information Criterion (BIC)
#'   is plotted for model comparison. Default is \code{FALSE}.
#' @param main Character string specifying the main title of the plot.
#' @param xlab Character string specifying the x-axis label.
#' @param hist.bins Integer; number of bins used for the histogram.
#'   Default is 50.
#' @param model Character string specifying the variance model used for clustering.
#'   If \code{NULL}, the model is selected automatically. Possible values include
#'   \code{"E"} (equal variance) and \code{"V"} (variable variance).
#' @param n.cluster Integer or integer vector specifying the number of clusters.
#'   If \code{NULL}, the optimal number of clusters is selected automatically.
#' @param apply.rounding Logical; if \code{TRUE}, results are rounded.
#'   Default is \code{TRUE}.
#' @param digits Integer specifying the number of decimal places used for rounding.
lab_mclust <- function(x, lognormal = FALSE, targets = NULL,
                       plot.it = TRUE, add.boxplot = TRUE, 
                       plot.legend = TRUE, pos.legend = "topright",
                       plot.bic = FALSE,
                       main = "", xlab = "",
                       hist.bins = 50, model = NULL, n.cluster = NULL,
                       apply.rounding = TRUE, digits = NULL){

  if(lognormal){
    xx <- log(x)
  }else{
    xx <- x
  }
  if(is.null(model)){
    mc <- Mclust(xx, G = n.cluster)
  }else{
    mc <- Mclust(xx, modelNames = model, G = n.cluster)
  }
  for (i in 1 : mc$G) {
    if(is.na(mc$parameters$variance$sigmasq[i])){
      mc$parameters$variance$sigmasq[i] <- mc$parameters$variance$sigmasq[1]
    }
  }
  res.tab <- data.frame(matrix(NA, mc$G, 5))
  for(i in 1 : mc$G){
    res.tab[i, 3] <- round(mc$parameters$pro[i] * 100, 1)
    res.tab[i, 4] <- round(mc$parameters$mean[i], 3)
    res.tab[i, 5] <- round(sqrt(mc$parameters$variance$sigmasq[i]), 3)
  }
  if(lognormal){
    colnames(res.tab) <- c("ll", "ul", "percent", "meanlog", "sdlog")
    for (i in 1 : mc$G) {
      res.tab[i, 1] <- qlnorm(0.025, mc$parameters$mean[i], sqrt(mc$parameters$variance$sigmasq[i]))
      res.tab[i, 2] <- qlnorm(0.975, mc$parameters$mean[i], sqrt(mc$parameters$variance$sigmasq[i]))
    }
  }else{
    colnames(res.tab) <- c("ll", "ul", "percent", "mean", "sd")
    for (i in 1 : mc$G) {
      res.tab[i, 1] <- qnorm(0.025, mc$parameters$mean[i], sqrt(mc$parameters$variance$sigmasq[i]))
      res.tab[i, 2] <- qnorm(0.975, mc$parameters$mean[i], sqrt(mc$parameters$variance$sigmasq[i]))
    }
  }
  if(is.null(digits)){
    digits <- 2 - floor(log10(median(xx)))
    if(digits < 0){digits <- 0}
  }
  if(apply.rounding){
    res.tab[, 1] <- round(res.tab[, 1], digits)
    res.tab[, 2] <- round(res.tab[, 2], digits)
  }
  
  if(plot.it){
    col = rainbow(9)
    d <- density(x)
    y.max <- max(d$y) * 1.1
    if(add.boxplot){y.max <- y.max  * 1.4}
    if(length(x) > 200){
      breaks <- seq(0.9 * min(x), 1.1 * max(x), length.out = hist.bins)
    }else{
      breaks <- "Sturges"
    }
    hist(x, freq = FALSE, 
         breaks = breaks,
         col = "white", border = "grey", 
         ylim = c(0, y.max),yaxt = "n",
         main = main, xlab = xlab, ylab = "")
    box()
    lines(d, lty = 2)
    if(add.boxplot){
      boxplot(x, horizontal = TRUE, at = y.max * 0.9, boxwex = y.max * 0.1, add = TRUE)
    }
    if(lognormal){
      for (i in 1 : nrow(res.tab)) {
        curve(dlnorm(x, meanlog = mc$parameters$mean[i], sdlog = sqrt(mc$parameters$variance$sigmasq[i])) * mc$parameters$pro[i],
              from = min(x), to = max(x), lwd = 2, col = col[i], add = TRUE)
      }
    }else{
      for (i in 1 : nrow(res.tab)) {
        curve(dnorm(x, mean = mc$parameters$mean[i], sd = sqrt(mc$parameters$variance$sigmasq[i])) * mc$parameters$pro[i],
              from = min(x), to = max(x), lwd = 2, col = col[i], add = TRUE)
      }
    }
    if(!is.null(targets)){
      lines(rep(targets[1], 2), c(0, y.max * 0.8), lty = 2)
      lines(rep(targets[2], 2), c(0, y.max * 0.8), lty = 2)
      text(targets, rep(y.max * 0.85, 2), targets)
    }
    if(plot.legend){
      legend(pos.legend, 
             paste0(round(res.tab$ll, digits), "-", round(res.tab$ul, digits),
                   " (", res.tab$percent, "%)"),
             lwd = 2, col = col[1 : nrow(res.tab)], cex = 0.8, bty = "n")
    }
    if(plot.bic){
      plot(mc$BIC)
    }
  }
  if(is.null(n.cluster)){n.c <- mc$G} else {
    n.c <- paste(mc$G, "from", deparse(n.cluster)) }
  return(list(n.cluster = noquote(n.c), 
              stats = res.tab, BIC = mc$BIC))
}

#' Computes the standard Box-Cox transformation.
#'
#' @param x Data to be transformed.
#' @param lambda The parameter of the Box-Cox transformation.
#'
#' @return The Box-Cox transformed data.
box.cox.trans <- function(x,lambda=1){
  if (lambda==0){
    return(log(x))
  }else{
    return((x^lambda - 1)/lambda)
  }
}

#' Computes the standard inverse Box-Cox transformation.
#'
#' @param x Data to be transformed.
#' @param lambda The parameter of the (inverse) Box-Cox transformation.
#'
#' @return The inverse Box-Cox transformed data.
box.cox.inv.trans <- function(x,lambda=1){
  if (lambda==0){
    return(exp(x))
  }else{
    return((x*lambda + 1)^(1/lambda))
  }
}

#' Computes a normal-approximation confidence interval for a quantile.
#'
#' @param p The target quantile (probability), e.g. 0.975.
#' @param alpha Significance level for the confidence interval.
#' @param n Sample size.
#'
#' @return A numeric vector of length 2 containing the lower and upper
#'   bounds of the confidence interval for the quantile.
i.norm <- function(p=0.975,alpha=0.1,n=120){
  half.width <- abs(qnorm(alpha/2)*sqrt(p*(1-p))/(dnorm(qnorm(p))*sqrt(n)))
  centre.point <- qnorm(p)
  return(c(centre.point-half.width,centre.point+half.width))
}

#' Transforms lower and upper limits using the Box-Cox transformation.
#'
#' @param ll Lower limit.
#' @param ul Upper limit.
#' @param lambda Parameter of the Box-Cox transformation.
#'
#' @return A numeric vector containing the transformed lower and upper
#'   limits.
compute.mu.sigma <- function(ll,ul,lambda=0){
  ll.t <- box.cox.trans(ll,lambda=lambda)
  ul.t <- box.cox.trans(ul,lambda=lambda)
  
  return(list(lower = ll.t, upper = ul.t))
}

#' Computes VeRUS-type confidence limits based on a Box-Cox transformed scale.
#'
#' The function transforms the provided lower and upper limits to a Box-Cox
#' scale, estimates the corresponding normal mean and standard deviation,
#' constructs normal-approximation confidence intervals for selected
#' quantiles, and transforms the results back to the original scale.
#'
#' @param ll Lower limit.
#' @param ul Upper limit.
#' @param lambda Parameter of the Box-Cox transformation.
#' @param delta Optional shift added after back-transformation.
#' @param p Numeric vector of probabilities for the lower and upper quantiles.
#' @param alpha Significance level for the confidence intervals.
#' @param n Sample size.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{lower.lim.low}{Lower bound of the CI for the lower quantile.}
#'   \item{lower.lim.upp}{Upper bound of the CI for the lower quantile.}
#'   \item{upper.lim.low}{Lower bound of the CI for the upper quantile.}
#'   \item{upper.lim.upp}{Upper bound of the CI for the upper quantile.}
#' }
verus.limits <- function(ll, ul, lambda, delta=0, p=c(0.025,0.975), alpha=0.1, n=120){
  ll.t <- box.cox.trans(ll,lambda=lambda)
  ul.t <- box.cox.trans(ul,lambda=lambda)
  
  mu <- (ll.t + ul.t)/2
  sigma <- (ul.t - ll.t)/(qnorm(0.975)-qnorm(0.025)) #factor 3.92
  
  ll <- box.cox.inv.trans(i.norm(p=p[1],alpha=alpha,n=n)*sigma + mu,lambda=lambda) + delta 
  ul <- box.cox.inv.trans(i.norm(p=p[2],alpha=alpha,n=n)*sigma + mu,lambda=lambda) + delta 
  
  return(list(lower.lim.low=ll[1],lower.lim.upp=ll[2],upper.lim.low=ul[1],upper.lim.upp=ul[2]))
}

#' Computes the mid-range coefficient of variation.
#'
#' @param ll Lower limit or lower value.
#' @param ul Upper limit or upper value.
#'
#' @return A numeric value representing the relative spread between the upper and lower values.
mocov <- function(ll, ul){
  return((ul - ll) / (ul + ll))
}

### Modified for VeRUS from the reflimR package ###
interpretation_VeRUS <- function(limits, targets, lambda){
  if(length(limits) != 2){stop("Limits must be a vector with length 2.")}
  if(limits[1] >= limits[2]){stop("The upper limit must be greater than the lower limit.")}
  if(limits[1] <= 0 | limits[2] <= 0){stop("Only positive limit values allowed.")}
  if(length(targets) != 2){stop("targets must be a vector with length 2.")}
  if(targets[1] >= targets[2]){stop("The upper target limit must be greater than the lower target limit.")}
  if(targets[1] <= 0 | targets[2] <= 0){stop("Only positive target values allowed.")}
  
  tol.lim <- unlist(verus.limits(limits[1], limits[2], lambda = lambda))
  col.lim <- c(lower.limit = rgb(0.7, 0.7, 0.7, 0.5), upper.limit = rgb(0.7, 0.7, 0.7, 0.5))
  dev.lim <- c(lower.limit = NA, upper.limit = NA)
  
  tol.tar <- unlist(verus.limits(targets[1], targets[2], lambda = lambda))
  col.tar <- col.lim
  if(limits[1] >= tol.tar[1] & limits[1] <= tol.tar[2]){
    col.lim[1] <-  rgb(0, 1, 0, 0.5)
  } else {
    if((tol.lim[1] >= tol.tar[1] & tol.lim[1] <= tol.tar[2]) |
       (tol.lim[2] >= tol.tar[1] & tol.lim[2] <= tol.tar[2])){
      col.lim[1] <- rgb(1, 1, 0, 0.5)
    } else {col.lim[1] <- rgb(1, 0, 0, 0.5)}
  }
  if(limits[2] >= tol.tar[3] & limits[2] <= tol.tar[4]){
    col.lim[2] <-  rgb(0, 1, 0, 0.5)
  } else {
    if((tol.lim[3] >= tol.tar[3] & tol.lim[3] <= tol.tar[4]) |
       (tol.lim[4] >= tol.tar[3] &  tol.lim[4] <= tol.tar[4])){
      col.lim[2] <- rgb(1, 1, 0, 0.5)
    } else {col.lim[2] <- rgb(1, 0, 0, 0.5)}
  }
  dev.lim[1 : 2] <- c("within tolerance", "within tolerance")
  for(i in 1 : 2){
    if(col.lim[i] == rgb(1, 1, 0, 0.5)){
      if(limits[i] < targets[i]){
        dev.lim[i] <- "slightly decreased"
      } else {
        dev.lim[i] <- "slightly increased"
      }
    }
    if(col.lim[i] == rgb(1, 0, 0, 0.5)){
      if(limits[i] < targets[i]){
        dev.lim[i] <- "markedly decreased"
      } else {
        dev.lim[i] <- "markedly increased"
      }
    }
  }
  return(list(tol.lim = tol.lim, tol.tar = tol.tar,
              col.lim = col.lim, col.tar = col.tar,
              dev.lim = dev.lim))
}

reflim_VeRUS <- function(x, lognormal = NULL, targets = NULL,
                   perc.trunc = 2.5, n.min = 200, apply.rounding = TRUE,
                   plot.it = TRUE, plot.all = FALSE, print.n = TRUE,
                   main = "reference limits", xlab = "x", lambda){
  
  xx <- na.omit(x)
  
  result <- list(stats = c(mean = NA, sd = NA, n.total = NA, n.trunc = NA),
                 lognormal = lognormal,
                 limits = c(lower.lim = NA, upper.lim = NA,
                            lower.lim.low = NA, lower.lim.upp = NA,
                            upper.lim.low = NA, upper.lim.upp = NA),
                 targets = c(lower.lim = NA, upper.lim = NA,
                             lower.lim.low = NA, lower.lim.up = NA,
                             upper.lim.low = NA, upper.lim.upp = NA),
                 perc.norm = NA,
                 confidence.int = c(lower.lim.low = NA, lower.lim.upp = NA,
                                    upper.lim.low = NA, upper.lim.upp = NA,
                                    n = NA),
                 interpretation = c(lower.limit = NA, upper.limit = NA),
                 remarks = NA)
  
  if(!is.numeric(xx)){
    warning("x must be numeric. Non-numeric values removed.")
    
    xx <- as.numeric(xx)
    xx <- na.omit(xx)
    
    result$remarks <- "Non-numeric values removed"
  }
  if(min(xx) <= 0){
    warning("Only positive values allowed. values <= 0 removed.")
    xx <- xx[xx > 0]
    result$remarks <- "Values <= 0 removed"
  }
  if(!is.null(targets)){
    
    targets <- na.omit(as.numeric(targets))
    
    if(length(targets) != 2){
      warning("targets must be a vector with length 2. NA not allowed. Targets removed.")
      targets = NULL
      result$remarks <- "Unsuitable target values removed"
    }
  }
  if(!is.null(targets)){
    if(is.na(targets[1]) | is.na(targets[2])){
      warning("Targets must be numeric. NA not allowed. Targets removed.")
      targets = NULL
      result$remarks <- "Unsuitable target values removed"
    }
  }
  if(!is.null(targets)){
    if(targets[1] >= targets[2]){
      warning("The upper target limit must be greater than the lower target limit. Targets removed. ")
      targets = NULL
      result$remarks <- "Unsuitable target values removed"
    }
  }
  if(!is.null(targets)){
    if(targets[1] <= 0 | targets[2] <= 0){
      warning("Only positive target values allowed. Targets removed.")
      targets = NULL
      result$remarks <- "Unsuitable target values removed"
    }
  }
  n <- length(xx)
  if(n < 40){
    warning(paste("n = ", n, ". The absolute minimum for reference limit estimation is 40. NAs returned."))
    result$stats[3] <- n
    result$remarks <- "Total n < 40"
    return(result)
  }
  if(n < n.min){
    warning(paste("n =", n, "where a minimum of", n.min, "is required. n.min has been set to 40 at a potential loss of accuracy."))
    result$stats[3] <- n
    result$remarks <- "Attention: low n."
    n.min <- 40
  }
  
  digits <- adjust_digits(median(xx))$digits
  if(is.null(lognormal)){
    plot.logtype <- TRUE
    lognormal <- lognorm(xx, plot.it = FALSE)$lognormal
  } else {
    plot.logtype <- FALSE
  }
  
  res.lognorm <- lognorm(xx, plot.it = FALSE)
  res.trunc <- iboxplot(xx, lognormal = lognormal,
                        perc.trunc = perc.trunc,
                        apply.rounding = apply.rounding,
                        plot.it = FALSE)
  n.trunc <- length(res.trunc$trunc)
  if(n.trunc < 40){
    warning(paste("n = ", n.trunc, "after truncation. The absolute minimum for reference limit estimation is 40. NAs returned."))
    result$stats[3] <- n
    result$stats[4] <- n.trunc
    result$remarks <- "n < 40 after truncation."
    return(result)
  }
  if(n.trunc < n.min){
    warning(paste("n.trunc =", n.trunc, "where a minimum of", n.min, "is required. n.min has been set to 40 at a potential loss of accuracy."))
    result$stats[3] <- n
    result$stats[4] <- n.trunc
    result$remarks <- "Low n after truncation."
    n.min <- 40
  }
  res.qq <- truncated_qqplot(res.trunc$trunc, lognormal = lognormal,
                             perc.trunc = perc.trunc, n.min = n.min,
                             apply.rounding = apply.rounding, plot.it = FALSE)$result
  res.ci <- conf_int95(n = n,
                       lower.limit = as.numeric(res.qq[3]),
                       upper.limit = as.numeric(res.qq[4]),
                       lognormal = lognormal, apply.rounding = apply.rounding)
  if(res.qq[3] > 0){
    res.pu <- unlist(verus.limits(as.numeric(res.qq[3]), as.numeric(res.qq[4]), lambda = lambda))
  } else {
    warning("Estimated lower limit <= 0. No tolerance limits calculated. No graphics produced.")
    res.pu <- rep(NA, 4)
    targets = NULL
    result$remarks <- "Lower limit <= 0"
  }
  
  res.lim <- c(as.numeric(res.qq[3 : 4]), as.numeric(res.pu))
  names(res.lim) <- c("lower.lim", "upper.lim", "lower.lim.low", "lower.lim.upp", "upper.lim.low", "upper.lim.upp")
  if(apply.rounding){res.lim <- round(res.lim, digits)}
  res.tar <- c(lower.lim = NA, upper.lim = NA,
               lower.lim.low = NA, lower.lim.up = NA,
               upper.lim.low = NA, upper.lim.upp = NA)
  dev.lim <- c(lower.limit = NA, upper.limit = NA)
  if(!is.null(targets)){
    ip <- interpretation_VeRUS(res.lim[1 : 2], targets, lambda)
    res.tar[1 : 2] <- targets
    res.tar[3 : 6] <- ip$tol.tar
    if(apply.rounding){res.tar <- round(res.tar, digits)}
    dev.lim <- ip$dev.lim
  }
  if(res.qq[3] > 0){
    if(plot.all){plot.it <- TRUE}
    if(plot.all){
      oldpar <- par(mfrow = c(2, 2))
      on.exit(par(oldpar))
    }
    if(plot.it){
      rh <- ri_hist_VeRUS(xx, lognormal = lognormal, stats = res.qq[1 : 2],
                    limits = res.qq[3 : 4], targets = targets,
                    perc.norm = res.trunc$perc.norm,
                    main = main, xlab = xlab, lambda = lambda)
      if(print.n){
        legend("topright", legend = paste("n =", n.trunc, "after truncation"),
               bty = "n", cex = 0.75)
      }
    }
    if(plot.all){
      lognorm(xx, main = "Step 1: Bowley skewness", xlab = "", plot.logtype = plot.logtype)
      iboxplot(xx, lognormal = lognormal, perc.trunc = perc.trunc,
               apply.rounding = apply.rounding,
               main = "Step 2: iBoxplot", xlab = "")
      truncated_qqplot(res.trunc$trunc, lognormal = lognormal,
                       perc.trunc = perc.trunc, n.min = n.min,
                       apply.rounding = apply.rounding,
                       main = "Step 3: Q-Q plot", xlab = "", ylab = "")
    }
  }
  
  result$stats = c(res.qq[1 : 2], n.total = n, n.trunc = n.trunc)
  result$lognormal = lognormal
  result$limits = res.lim
  result$targets = res.tar
  result$perc.norm = res.trunc$perc.norm
  result$confidence.int = res.ci[1 : 4]
  result$interpretation = dev.lim
  
  return(result)
}

ri_hist_VeRUS <- function(x, lognormal, stats, limits, perc.norm,
                    targets = NULL, remove.extremes = TRUE,
                    main = "reflim", xlab = "x", lambda){
  
  xx <- na.omit(x)
  if(!is.numeric(xx)){stop("x must be numeric.")}
  if(min(xx) <= 0){stop("x must be a vector of positive numbers.")}
  if(length(stats) != 2){stop("stats must be a vector with length 2 containing mean (or meanlog) and sd (or sdlog).")}
  if(length(limits) != 2){stop("limits must be a vector with length 2.")}
  if(is.numeric(limits)){
    if(limits[1] >= limits[2]){stop("The upper limit must be greater than the lower limit.")}
    if(limits[1] <= 0 | limits[2] <= 0){stop("Only positive limit values allowed.")}
  }
  if(!is.null(targets)){
    if(length(targets) != 2){stop("targets must be a vector with length 2.")}
    if(targets[1] >= targets[2]){stop("The upper target limit must be greater than the lower target limit.")}
    if(targets[1] <= 0 | targets[2] <= 0){stop("Only positive limit values allowed.")}
  }
  
  digits <- adjust_digits(median(xx))$digits
  n <- length(xx)
  if(n < 40){{stop(paste0("n = ", n, ". The absolute minimum for reference limit estimation is 40."))}}
  
  if(remove.extremes){xx <- xx[xx <= median(xx) + 8 * IQR(xx)]}
  if(n < 200){breaks <- "Sturges"} else {
    difference <- max(xx) - min(xx)
    from = min(xx) - 0.1 * difference
    if(from < 0) {from <- 0}
    to = max(xx) + 0.1 * difference
    if(to < limits[2]) { to <- limits[2] + 0.1 * difference}
    breaks <- seq(from = from, to = to, by = (limits[2] - limits[1]) / 10)
  }
  d <- density(xx)
  d.max <- max(d$y)
  
  if(lognormal){
    d1 <- dlnorm(d$x, stats[1], stats[2])
  } else {
    d1 <- dnorm(d$x, stats[1], stats[2])
  }
  
  tol.lim <- unlist(verus.limits(limits[1], limits[2], lambda = lambda))
  col.lim <- c(rgb(0.7, 0.7, 0.7, 0.5), rgb(0.7, 0.7, 0.7, 0.5))
  dev.lim <- c(lower.limit = NA, upper.limit = NA)
  
  if(!is.null(targets)){
    ip <- interpretation_VeRUS(limits, targets, lambda)
    tol.lim <- ip$tol.lim
    tol.tar <- ip$tol.tar
    col.lim <- ip$col.lim
    col.tar <- ip$col.tar
    dev.lim <- ip$dev.lim
  }
  
  hist(xx, freq = FALSE, breaks = breaks, yaxt = "n",
       ylim = c(0, max(d.max, max(d1 * perc.norm / 100))),
       col = "white", border = "grey",
       main = main, xlab = xlab, ylab = "")
  
  if(!is.null(targets)){
    rect(tol.tar[1], 0, tol.tar[2], d.max * 0.8, col = col.tar[1], border = NA)
    rect(tol.tar[3], 0, tol.tar[4], d.max * 0.8, col = col.tar[2], border = NA)
  }
  rect(tol.lim[1], 0, tol.lim[2], d.max * 0.8, col = col.lim[1], border = NA)
  rect(tol.lim[3], 0, tol.lim[4], d.max * 0.8, col = col.lim[2], border = NA)
  
  lines(d, lty = 3)
  lines(d$x, d1 * perc.norm / 100, lwd = 2, col = "blue")
  d2 <- d$y - d1
  d2[d2 < 0] <- 0
  lines(d$x, d2, col = "red", lty = 2)
  
  lines(rep(limits[1], 2), c(0, d.max * 0.8), lty = 2)
  lines(rep(limits[2], 2), c(0, d.max * 0.8), lty = 2)
  text(limits, rep(d.max * 0.85, 2), round(limits, digits))
  if(!is.null(targets)){
    lines(rep(targets[1], 2), c(0, d.max * 0.8))
    lines(rep(targets[2], 2), c(0, d.max * 0.8))
    text(targets[1 : 2], rep(d.max * 0.92, 2),
         round(targets[1 : 2], digits), col = "grey")
  }
  return(list(lognormal = lognormal, percent_normal = perc.norm, interpretation = dev.lim))
}