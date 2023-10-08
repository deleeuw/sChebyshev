
vSmacof <-
  function(delta,
           mbas,
           w = oneDist(attr(delta, "Size")),
           told = rep(1, length(mbas)),
           itmax = 1000,
           eps1 = 1e-15,
           eps2 = 1e-15,
           verbose = TRUE) {
    w <- w / sum(w)
    delta <- delta / sqrt(sum(w * (delta ^ 2)))
    m <- length(mbas)
    v <- mmatrix(w)
    vv <- matrix(0, m, m)
    for (i in 1:m) {
      for (j in 1:m) {
        vv[i, j] <- sum(mbas[[i]] * (v %*% mbas[[j]]))
      }
    }
    vi <- ginv(vv)
    zold <- listSum(mbas, told)
    dold <- dist(zold)
    sold <- sum(w * (delta - dold) ^ 2) / 2.0
    itel <- 1
    repeat {
      bold <- mmatrix(w * (dlmat / dold))
      bb <- matrix(0, m, m)
      for (i in 1:m) {
        for (j in 1:m) {
          bb[i, j] <- sum(mbas[[i]] * bold %*% mbas[[j]])
        }
      }
      tnew <- vi %*% bb %*% told
      znew <- listSum(mbas, tnew)
      dnew <- dist(znew)
      snew <- sum(w * (delta - dnew) ^ 2) / 2.0
      chng <- sqrt(sum((told - tnew) ^ 2))
      if (verbose) {
        cat(
          "itel ",
          formatC(itel, digits = 4, format = "d"),
          "sold ",
          formatC(sold, digits = 15, format = "f"),
          "snew ",
          formatC(snew, digits = 15, format = "f"),
          "chng ",
          formatC(chng, digits = 15, format = "f"),
          "\n"
        )
      }
      if ((itel == itmax) || (((sold - snew) < eps1) && (chng < eps2))) {
        break
      }
      told <- tnew 
      zold <- znew
      sold <- snew
      dold <- dnew
      itel <- itel + 1
    }
    return(list(
      c = tnew,
      z = znew,
      d = dnew,
      s = snew,
      itel = itel
    ))
  }
