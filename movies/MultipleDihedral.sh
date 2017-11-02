#!/bin/bash

touch output.csv # file with the end result

files=$(find . -name '*.xyz' -exec echo {} \;)
for file in $files; do # replace with your "for every file loop"
  touch tempfile1
  touch tempfile2
  g++ -std=c++14 DihedralTool.cpp && ./a.out $file $1 $2 $3 $4 | tr ' ' '\n' > tempfile1 # program.sh is your software, we pipe it into tr to split it into a column
  paste -d ',' tempfile1 output.csv > tempfile2 # paste the column with a delimeter to a tempfile
  mv tempfile2 output.csv # clean up the tempfiles on every iteration
  rm tempfile1
  echo $file
done