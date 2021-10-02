# 2d-Ising-model-Fortran-90
A Fortran 90 code for the 2-dim Ising model

Explanation:

main.f90 : --> is the main program, all the necessary tasks to be performed, are written as subroutines in the file named "include_defs.f90"
and, "include_mods.f90" contains the modules defined that are necessary for the data moving across the routines/ or defining a kind of global data, that every subroutine can see. 

To run/ execute the code:

1. keep all the code files that are 3 fortran files and also the Makefile, in the same directory.
2. To compile & create the executable, open a terminal(hope you're using Linux). and go to that path of the folder where you kept all the above files. 
3. type "make" and hit enter.
4. It will create an executable for you (by default, it will create "run_ising_2d.exe"), which can run the program.
5. now, just type "./run_ising_2d.exe" in the terminal, it will start running your code.
6. Congratulations, your code is running now, it will also create a few output data files, from this simulation, which you may need to use for the purpose of plotting and any further analysis.
