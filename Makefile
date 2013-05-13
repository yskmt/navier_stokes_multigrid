# Makefile: CSE391 Final Project
# Geometric Multigrid

CC=g++
CFLAGS=-fopenmp  -std=c++11
LDFLAGS=-fopenmp  -std=c++11
SOURCES=main.C assemble.C  IO.C  jacobi.C utils.C  v_cycle.C advection.C viscosity.C msort.C pressure.C
OBJ=$(SOURCES:.C=.o)
EXE=multigrid

all: $(SOURCES) $(EXE)

$(EXE): $(OBJ) 
	$(CC) $(LDFLAGS) $(OBJ) -o $@

.C.o:
	$(CC) $(CFLAGS) -c $< -o $@

jacobi.o:jacobi.h
main.o:jacobi.h assemble.h utils.h IO.h v_cycle.h advection.h viscosity.h
viscosity.o: viscosity.h
pressure.o: pressure.h
IO.o: IO.h advection.h
v_cycle.o: v_cycle.h pressure.h

run:
	@ export OMP_NUM_THREADS=1
	@./multigrid 1 1 12 12 12

clean:
	@ rm *.o multigrid

