source("utils.R")
source("basis.R")
source("exampleRun.R")
source("derivatives.R")

library(MASS)
library(microbenchmark)
library(zeallot)

smacofSetup <- function(delta, w, xold) {
  w <- w / sum(w)
  v <- mmatrix(w)
  dold <- dist(xold)
  delta <- delta / sqrt(sum(w * (delta ^ 2)))
  lbd <- sum(w * delta * dold) / sum(w * dold ^ 2)
  xold <- lbd * xold
  dold <- lbd * dold
  sold <- sum(w * (delta - dold) ^ 2) / 2
  return(list(delta, w, v, xold, dold, sold))
}

sSmacofConfs <- function(delta, w, xold, sold, epsconf, memory, verbose) {
  guttlist <- list(xold)
  repeat {
    xnew <- guttman(delta, xold, w = w)
    chng <- rms(xold - xnew)
    dnew <- dist(xnew)
    snew <- sum(w * (delta - dnew) ^ 2) / 2
    guttlist <- c(guttlist, list(xnew))
    if (verbose > 1L) {
      cat(" confs ", formatC(c(sold, snew, chng), digits = 15, format = "f"), "\n")
    }
    if ((chng < epsconf) || (length(guttlist) == memory)) {
      break
    }
    xold <- xnew
    sold <- snew
  }
  return(list(xnew, guttlist))
}

sSmacofSteps <-
  function(delta, w, v, xold, guttlist, otmax, epsstep, verbose) {
    m <- length(guttlist)
    told <- ei(m, m)
    dold <- dist(xold)
    sold <- sum(w * (delta - dold) ^ 2) / 2
    bold <- mmatrix(w * (delta / dold))
    vcon <- listCrossprod(guttlist, v)
    print(eigen(vcon)$values)
    otel <- 1
    repeat {
      bcon <- listCrossprod(guttlist, bold)
      print(eigen(bcon)$values)
      tnew <- drop(ginv(vcon) %*% bcon %*% told)
      chng <- rms(told - tnew)
      xnew <- listSum(guttlist, tnew)
      dnew <- dist(xnew)
      snew <- sum(w * (delta - dnew) ^ 2) / 2
      if (verbose > 1L) {
        cat(" steps ", formatC(c(sold, snew, chng), digits = 15, format = "f"), "\n")
      }
      bnew <- mmatrix(w * delta / dnew)
      if ((otel == otmax) || (chng < epsstep)) {
        break
      }
      otel <- otel + 1
      bold <- bnew
      told <- tnew
      sold <- snew
    }
    return(list(tnew, xnew, dnew, snew, bcon, vcon))
  }

sSmacofStaps <- function() {
  
}

sSmacof <-
  function(delta,
           p = 2,
           xold = torgerson(delta, p),
           w = oneDist(attr(delta, "Size")),
           itmax = 10000L,
           otmax = 1L,
           memory = 2L,
           epsover = 1e-6,
           epsconf = 1e-6,
           epsstep = 1e-6,
           verbose = 0L) {
    c(delta, w, v, xold, dold, sold) %<-%
      smacofSetup(delta, w, xold)
    itel <- 1
    repeat {
      c(xnew, guttlist) %<-%
        sSmacofConfs(delta, w, xold, sold, epsconf, memory, verbose)
      c(tnew, xnew, dnew, snew, bcon, vcon) %<-%
        sSmacofSteps(delta, w, xold = xnew, v, guttlist, otmax, epsstep, verbose)
      chng <- rms(xold - xnew)
      if (verbose > 0L) {
        cat(
          "itel ",
          formatC(itel, format = "d"),
          "snew ",
          formatC(snew, digits = 10, format = "f"),
          "diff ",
          formatC(sold - snew, digits = 10, format = "f"),
          "chng ",
          formatC(chng, digits = 10, format = "f"),
          "\n"
        )
      }
      if ((chng < epsover) || (itel == itmax)) {
        break
      }
      itel <- itel + 1
      xold <- xnew
      sold <- snew
    }
    return(list(
      xmat = xnew,
      loss = snew,
      thet = tnew,
      itel = itel,
      bcon = bcon,
      vcon = vcon
      ))
  }
