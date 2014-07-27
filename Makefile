all: sphericalinterpolate

CFLAGS=-g
CC=g++ -g

stripack.o: stripack.cc
	$(CC) -c stripack.cc

ssrfpack.o: ssrfpack.cc
	$(CC) -c ssrfpack.cc

SphericalSplineInterpolator.o: SphericalSplineInterpolator.cc SphericalSplineInterpolator.hh
	$(CC) -c SphericalSplineInterpolator.cc

main.o: main.cc
	$(CC) -c main.cc

sphericalinterpolate: stripack.o ssrfpack.o main.o SphericalSplineInterpolator.o 
	$(CC) -o sphericalinterpolate stripack.o ssrfpack.o SphericalSplineInterpolator.o  main.o -lf2c -lm -lstdc++ -g

clean:
	rm -f *.o sphericalinterpolate
