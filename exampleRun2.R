exampleRun2 <- function(delta) {
  args <- list(
    list(delta, otmax = 1L, memory = 1L),
    list(delta, otmax = 5L, memory = 1L),
    list(delta, otmax = 10L, memory = 1L),
    list(delta, otmax = 100L, memory = 1L),
    list(delta, otmax = 1L, memory = 2L),
    list(delta, otmax = 5L, memory = 2L),
    list(delta, otmax = 10L, memory = 2L),
    list(delta, otmax = 100L, memory = 2L),
    list(delta, otmax = 1L, memory = 3L),
    list(delta, otmax = 5L, memory = 3L),
    list(delta, otmax = 10L, memory = 3L),
    list(delta, otmax = 100L, memory = 3L),
    list(delta, otmax = 1L, memory = 4L),
    list(delta, otmax = 5L, memory = 4L),
    list(delta, otmax = 10L, memory = 4L),
    list(delta, otmax = 100L, memory = 4L)
  )
  n <- length(args)
  run <- as.list(1:n)
  trun <- as.list(1:n)
  for (i in 1:n) {
    run[[i]] <- do.call("sSmacof", args[[i]])
    trun[[i]] <-
      suppressWarnings(fivenum(microbenchmark(do.call(
        "sSmacof", args[[i]]
      ), times = 11L)$time / (10 ^ 6)))
  }
  results <- matrix(0, n, 3)
  row.names(results) <-
    c(
      "otmax = 1,    memory = 1",
      "otmax = 5,    memory = 1",
      "otmax = 10,   memory = 1",
      "otmax = 100,  memory = 1",
      "otmax = 1,    memory = 2",
      "otmax = 5,    memory = 2",
      "otmax = 10,   memory = 2",
      "otmax = 100,  memory = 2",
      "otmax = 1,    memory = 3",
      "otmax = 5,    memory = 3",
      "otmax = 10,   memory = 3",
      "otmax = 100,  memory = 3",
      "otmax = 1,    memory = 4",
      "otmax = 5,    memory = 4",
      "otmax = 10,   memory = 4",
      "otmax = 100,  memory = 4"
    )
  colnames(results) <-
    c("itel", "stress", "time")
  for (i in 1:n) {
    results[i,] <-
      c(run[[i]]$itel,
        run[[i]]$loss,
        trun[[i]][[3]])
  }
  return(list(results = results, run = run))
}
