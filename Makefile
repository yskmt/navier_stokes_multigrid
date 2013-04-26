# Makefile: CSE391 Final Project
# Geometric Multigrid

multigrid: main.C
	g++ main.C -fopenmp -o multigrid

run:
	@ export OMP_NUM_THREADS=1
	@./multigrid

run2:
	@ export OMP_NUM_THREADS=2
	@./mis_shared

