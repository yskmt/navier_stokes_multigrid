# Makefile: CSE391 Final Project
# Geometric Multigrid

multigrid: jacobi.o main.o
	g++ jacobi.o main.o -fopenmp -o multigrid

jacobi.o: jacobi.C jacobi.h
	g++ jacobi.C -c -fopenmp

main.o: main.C
	g++ main.C -c -fopenmp

run:
	@ export OMP_NUM_THREADS=1
	@./multigrid

run2:
	@ export OMP_NUM_THREADS=2
	@./mis_shared

