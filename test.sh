#!/bin/sh

echo 'n = 1000 p = 1'
./dense.out 1000 1 > t1_1
echo 'n = 1000 p = 3'
./dense.out 1000 3 > t1_3
echo 'n = 1000 p = 6'
./dense.out 1000 6 > t1_6
echo 'n = 3000 p = 1'
./dense.out 3000 1 > t2_1
echo 'n = 3000 p = 3'
./dense.out 3000 3 > t2_3
echo 'n = 3000 p = 6'
./dense.out 3000 6 > t2_6
echo 'n = 5000 p = 1'
./dense.out 5000 1 > t3_1
echo 'n = 5000 p = 3'
./dense.out 5000 3 > t3_3
echo 'n = 5000 p = 6'
./dense.out 5000 6 > t3_6

rm -f result.txt
touch result.txt
cat t1_1 t1_3 t1_6 t2_1 t2_3 t2_6 t3_1 t3_3 t3_6 > result.txt
rm t1_*
rm t2_*
rm t3_*
