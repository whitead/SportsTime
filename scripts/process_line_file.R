#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 2) {
    print("Usage: [input csv] [beta-binomial input file] [start game] [end game]")
    q()
}

d <- read.csv(argv[1])

if(length(argv) == 2) {
    start <- 1
    end <- nrow(d)
} else if(length(argv) == 3) {
    start <- 1
    end <- as.integer(argv[3])
} else {
    start <- as.integer(argv[3])
    end <- as.integer(argv[4])
}

games <- matrix(ncol=32, nrow=32, 0)
row.names(games) <- levels(d$Home.Team)
colnames(games) <- levels(d$Home.Team)
for(i in start:end) {

        
    if(d$Home.Score[i] > d$Visitor.Score[i])
        games[d$Home.Team[i], d$Visitor[i]] <- games[d$Home.Team[i], d$Visitor[i]] + 1
    else
        games[d$Visitor[i], d$Home.Team[i]] <- games[d$Visitor[i], d$Home.Team[i]] + 1
}
write.table(games, file=paste(argv[2], "_wins.txt", sep=""), row.names=FALSE, col.names=FALSE)

cat("Wins:")
print(apply(games, 1, sum))

games <- matrix(ncol=32, nrow=32, 0)
for(i in end:nrow(d)) {
  games[d$Home.Team[i], d$Visitor[i]] <- games[d$Home.Team[i], d$Visitor[i]] + 1
  games[d$Visitor[i], d$Home.Team[i]] <- games[d$Visitor[i], d$Home.Team[i]] + 1 
}
if(nrow(d) - end > 1)
  write.table(games, file=paste(argv[2], "_remaining.txt", sep=""), row.names=FALSE, col.names=FALSE)
