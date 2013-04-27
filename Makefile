# Makefile: CSE391 Final Project
# Geometric Multigrid

multigrid: jacobi.o main.o assemble.o utils.o
	g++ jacobi.o main.o assemble.o utils.o -fopenmp -o multigrid

jacobi.o: jacobi.C jacobi.h
	g++ jacobi.C -c -fopenmp

assemble.o: assemble.C assemble.h utils.C utils.h
	g++ assemble.C -c -fopenmp

utils.o: utils.C utils.h
	g++ utils.C -c -fopenmp

main.o: main.C
	g++ main.C -c -fopenmp

run:
	@ export OMP_NUM_THREADS=1
	@./multigrid

clean:
	@ rm *.o multigrid

run2:
	@ export OMP_NUM_THREADS=2
	@./mis_shared

