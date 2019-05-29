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
	$(CC) $(OPTION) -o sparse.out $< mmreader.o

sparse.o : sparse.cu mmreader.hpp
	$(CC) $(OPTION) --std=c++11 -c $< 
mmreader.o : mmreader.cpp mmreader.hpp
	g++ $(OPTION) --std=c++11 -c $< 

clean :
	rm -f dense.o dense.out
	rm -f gaussian.o gaussian.out
	rm -f sparse.o sparse.out
	rm -f mmreader.o mmreader.out
	

