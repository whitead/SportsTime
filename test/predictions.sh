#! /bin/bash

N=162

for i in `seq 100 5 $(($N - 1))`; do
#    echo "Training with $i"
    ../scripts/process_line_file.R ~/Dropbox/SportsData/nfllines/nfl2013lines.csv train $i
#    echo "Testing on $i through $N"
    ../scripts/process_line_file.R ~/Dropbox/SportsData/nfllines/nfl2013lines.csv test $((i+1)) $((i+14))
    ../src/beta-binomial -w train_wins.txt -m 32 -n 2500 -p test_games.txt 
    ../scripts/compare_predict.R predicted_wins.txt test_wins.txt
done;
