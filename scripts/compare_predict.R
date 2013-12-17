#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 2) {
    print("Usage: [predicted_wins] [actual_wins]")
    q()
}

explode <- function(v, h) {
  result <- rep(0, sum(h))
  index <- 1
  for(i in 1:length(v)) {
    if(h[i] > 0) {
      for(j in 1:h[i]) {
        result[index] <- v[i]
        index <- index + 1
      }
    }
  }
  return(result)
}


predicted <- read.table(argv[1], header=T)
actual <- read.table(argv[2], header=F)


m <- ncol(actual)

pindex <- 2 #skip first row
error <- rep(0,ncol(predicted) - 1)
games <- rep(0,ncol(predicted) - 1)
for(i in 1:m) {
  if(i < m) {
    for(j in (i+1):m) {
      games[pindex - 1] <- actual[i,j] + actual[j,i]
      p <- median(explode(predicted[,1] - 0.5, predicted[,pindex]))
      error[pindex - 1] <- abs(actual[i,j] - p)
      paste(actual[i,j],  p, "\n")
      pindex <- pindex + 1
    }
  }
}

cat(paste("games:", sum(games), "\n"))
cat(paste("correct:", sum(games) - sum(error), "\n"))
cat(paste("incorrect:", sum(error), "\n"))
