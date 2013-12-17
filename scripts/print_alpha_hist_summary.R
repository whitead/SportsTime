#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 3) {
    print("Usage: [alpha_hist] [team_labels] [square_win_matrix]")
    q()
}

hist <- read.table(argv[1], header=T)
labels <- read.table(argv[2])

m <- ncol(hist) - 1

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

team.strengths <- rep(0, nrow(labels))
names(team.strengths) <- labels[,1]

team.data <- matrix(0, nrow=sum(hist[,2]), ncol=nrow(labels))
colnames(team.data) <- labels[,1]

#explode the histograms
for(i in 2:ncol(hist)) {
  v <- explode(hist[,1], hist[,i])
  team.data[,i - 1] <- v
  team.strengths[i - 1] <- median(v)
}

#make win loss records
wins.matrix <- read.table(argv[3])
wins <- rep(0, m)
losses <- rep(0, m)
for(i in 1:m) {
  wins[i] <- sum(wins.matrix[i,])
  losses[i] <- sum(wins.matrix[,i])
}


#sort the teams and make labels
team.strengths.labels <- rep("", nrow(labels))
for(i in 1:m) {
  index <- order(team.strengths, decreasing=TRUE)[i]
  cat(paste(i, labels[index,1], team.strengths[index], "\n"))
  team.strengths.labels[i] <- paste(labels[index,1],
                                    paste("(",
                                          wins[index], " - ",
                                          losses[index],")",
                                          sep=""),
                                    format(i, width=2))
}

cairo_pdf("ordered.pdf", width=7, height=7)
par(family="LMSans10", fg="dark gray", mar=c(4, 15, 3, 2))
palette <- colorRampPalette(c("blue", "red"))
boxplot(team.data[,order(team.strengths)], horizontal=T, las=1, xlab=expression(alpha), col=palette(m), border="black", outline=F,
        names=rev(team.strengths.labels))
graphics.off()
