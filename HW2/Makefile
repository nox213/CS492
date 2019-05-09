OPTION = -O2 -fopenmp

all : dense.out gaussian.out sparse.out

dense.out : dense.o
	gcc $(OPTION) -o dense.out -lm  -g $<

dense.o : dense.c
	gcc $(OPTION) -c -lm  -g $<

gaussian.out : gaussian.o
	gcc $(OPTION) -o gaussian.out $< -lm

gaussian.o : gaussian.c
	gcc $(OPTION) -c $< -lm

sparse.out : sparse.o mmreader.o mmreader.hpp
	g++ $(OPTION) -o sparse.out -std=c++11 $< mmreader.o

sparse.o : sparse.cpp mmreader.hpp
	g++ $(OPTION) -c $< -std=c++11
mmreader.o : mmreader.cpp mmreader.hpp
	g++ $(OPTION) -c $< -std=c++11

clean :
	rm -f dense.o dense.out
	rm -f gaussian.o gaussian.out
	rm -f sparse.o sparse.out
	rm -f mmreader.o mmreader.out
	

