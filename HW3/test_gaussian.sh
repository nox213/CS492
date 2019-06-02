#!/bin/sh

echo 'n = 1000'
./gaussian.out 1000 > t1
echo 'n = 3000'
./gaussian.out 3000 > t2
echo 'n = 5000'
./gaussian.out 5000 > t3

rm -f result_gaussian.txt
touch result_gaussian.txt
cat t1 t2 t3 > result_gaussian.txt
rm t1
rm t2
rm t3
