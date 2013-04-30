# Makefile: CSE391 Final Project
# Geometric Multigrid

CC:=gcc

multigrid: jacobi.o main.o assemble.o utils.o IO.o
	$(CC) jacobi.o main.o assemble.o utils.o IO.o -fopenmp -std=c++11 -o multigrid

jacobi.o: jacobi.C jacobi.h
	$(CC) jacobi.C -c -fopenmp -std=c++11

assemble.o: assemble.C assemble.h utils.C utils.h
	$(CC) assemble.C -c -fopenmp -std=c++11

utils.o: utils.C utils.h
	$(CC) utils.C -c -fopenmp -std=c++11

IO.o: IO.C IO.h
	$(CC) IO.C -c -fopenmp -std=c++11

main.o: main.C
	$(CC) main.C -c -fopenmp -std=c++11

run:
	@ export OMP_NUM_THREADS=1
	@./multigrid

clean:
	@ rm *.o multigrid

run2:
	@ export OMP_NUM_THREADS=2
	@./mis_shared

