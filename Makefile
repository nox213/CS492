all : dense.out gaussian.out

dense.out : dense.o
	gcc -o dense.out -pthread -lm $<

dense.o : dense.c
	gcc -c $<

gaussian.out : gaussian.o
	gcc -o gaussian.out -pthread -lm $<

gaussian.o : gaussian.c
	gcc -c $<

clean :
	rm -f dense.o dense.out
	rm -f gaussian.o gaussian.out
	

