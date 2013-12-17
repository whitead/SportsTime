#! /bin/bash

N=25

for i in `seq 10 $(($N - 1))`; do
    ../scripts/process_line_file.R ~/Dropbox/SportsData/nfllines/nfl2013lines.csv train $i
 ../scripts/process_line_file.R ~/Dropbox/SportsData/nfllines/nfl2013lines.csv test $i $(($N - $i))
 ../src/beta-binomial -w train_wins.txt -m 32 -n 200 -p train_remaining.txt 
 ../scripts/compare_predict.R predicted_wins.txt test_wins.txt
done;