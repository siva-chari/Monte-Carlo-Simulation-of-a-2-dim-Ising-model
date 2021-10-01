SRC = include_mods.f90 include_defs.f90 main.f90
#SRC = $(wildcard *.f90)
OBJ = $(SRC:.f90=.o)
EXE = run_ising_2d.exe

CC =       gfortran
CCFLAGS =  -O3 -std=f2008
LDFLAGS =       #-lblas -latlas -llapack
RM =	   rm -rf

$(EXE): $(OBJ)
	$(CC) -o $(EXE) $(OBJ) $(CCFLAGS) $(LDFLAGS)

%.o: %.f90
	$(CC) -c $(CCFLAGS) $*.f90 $(CCFLAGS) $(LDFLAGS)

.PHONY: clean

clean:
	$(RM) $(EXE) $(OBJ) *.mod

