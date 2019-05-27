OPTION =
CC = nvcc

all : dense.out gaussian.out sparse.out

dense.out : dense.o
	$(CC) $(OPTION) -o dense.out -lm $<

dense.o : dense.cu
	$(CC) $(OPTION) -c -lm  $<

gaussian.out : gaussian.o
	$(CC) $(OPTION) -o gaussian.out $< -lm

gaussian.o : gaussian.cu
	$(CC) $(OPTION) -c $< -lm

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
	

