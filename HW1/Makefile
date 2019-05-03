all : dense.out gaussian.out sparse.out

dense.out : dense.o
	gcc -O3 -o dense.out -pthread -lm $<

dense.o : dense.c
	gcc -c -O3 -lm $<

gaussian.out : gaussian.o
	gcc -O3 -o gaussian.out -pthread $< -lm

gaussian.o : gaussian.c
	gcc -c -O3 $< -lm

sparse.out : sparse.o mmreader.o mmreader.hpp
	g++ -O3 -o sparse.out -pthread -lm -std=c++11 $< mmreader.o

sparse.o : sparse.cpp mmreader.hpp
	g++ -c -O3 $< -std=c++11
mmreader.o : mmreader.cpp mmreader.hpp
	g++ -c $< -std=c++11

clean :
	rm -f dense.o dense.out
	rm -f gaussian.o gaussian.out
	rm -f sparse.o sparse.out
	rm -f mmreader.o mmreader.out
	

