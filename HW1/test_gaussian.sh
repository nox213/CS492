#!/bin/sh

echo 'n = 1000 p = 2' > t1_1
./gaussian.out 1000 2 > t1_1
echo '----------------------'
echo 'n = 1000 p = 4'
./gaussian.out 1000 4 > t1_3
echo '----------------------'
echo 'n = 1000 p = 6'
./gaussian.out 1000 6 > t1_6
echo '----------------------' 
echo 'n = 3000 p = 2'
./gaussian.out 3000 2 > t2_1
echo '----------------------'
echo 'n = 3000 p = 4'
./gaussian.out 3000 4 > t2_3
echo '----------------------'
echo 'n = 3000 p = 6'
./gaussian.out 3000 6 > t2_6
echo '----------------------'
echo 'n = 5000 p = 2'
./gaussian.out 5000 2 > t3_1
echo '----------------------'
echo 'n = 5000 p = 4'
./gaussian.out 5000 4 > t3_3
echo '----------------------'
echo 'n = 5000 p = 6'
./gaussian.out 5000 6 > t3_6

rm -f result_gaussian.txt
touch result_gaussian.txt
cat t1_1 t1_3 t1_6 t2_1 t2_3 t2_6 t3_1 t3_3 t3_6 > result_gaussian.txt
rm t1_*
rm t2_*
rm t3_*
