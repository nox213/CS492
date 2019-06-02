#!/bin/sh

echo 'n = 1000'
./dense.out 1000> t1
echo 'n = 3000'
./dense.out 3000 > t2
echo 'n = 5000'
./dense.out 5000 > t3

rm -f result.txt
touch result.txt
cat t1 t2 t3 > result.txt
rm t1
rm t2
rm t3
