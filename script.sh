#/bin/bash

./sparse.out ./matrix/2cubes_sphere.mtx 2048 $1 > t1
echo 'done'
./sparse.out ./matrix/cage12.mtx 1024 $1 > t2
echo 'done'
./sparse.out ./matrix/consph.mtx 2048 $1 > t3
echo 'done'
./sparse.out ./matrix/cop20k_A.mtx 2048 $1 > t4
echo 'done'
./sparse.out ./matrix/filter3D.mtx 2048 $1 > t5
echo 'done'
./sparse.out ./matrix/hood.mtx 1024 $1 > t6
echo 'done'
./sparse.out ./matrix/m133-b3.mtx 1024 $1 > t7
echo 'done'
./sparse.out ./matrix/mac_econ_fwd500.mtx 1024 $1 > t8
echo 'done'
./sparse.out ./matrix/majorbasis.mtx 1024 $1 > t9
echo 'done'
./sparse.out ./matrix/mario002.mtx 512 $1 > t10
echo 'done'
./sparse.out ./matrix/mc2depi.mtx 512 $1 > t11
echo 'done'
./sparse.out ./matrix/offshore.mtx 1024 $1 > t12
echo 'done'
./sparse.out ./matrix/patents_main.mtx 1024 $1 > t13
echo 'done'
./sparse.out ./matrix/pdb1HYS.mtx 4096 $1 > t14
echo 'done'
./sparse.out ./matrix/poisson3Da.mtx 16384 $1 > t15
echo 'done'
./sparse.out ./matrix/pwtk.mtx 1024 $1 > t16
echo 'done'
./sparse.out ./matrix/rma10.mtx 4096 $1 > t17
echo 'done'
./sparse.out ./matrix/scircuit.mtx 1024 $1 > t18
echo 'done'
./sparse.out ./matrix/shipsec1.mtx 1024 $1 > t19
echo 'done'
./sparse.out ./matrix/webbase-1M.mtx 256 $1 > t20
echo 'done'

cat t1 t2 t3 t4 t5 t6 t7 t8 t9 t10 t11 t12 t13 t14 t15 t16 t17 t18 t19 t20 > sparse_$1.txt
rm t1 t2 t3 t4 t5 t6 t7 t8 t9 t10 t11 t12 t13 t14 t15 t16 t17 t18 t19 t20

