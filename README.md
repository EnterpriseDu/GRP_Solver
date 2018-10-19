***A Brief Instruction***
  1. Enter the directory 'Riemann_solver/'.
  2. Run 'make'.
  3. Include the header 'Riemann_solver/inc/Riemann_solver.h' in the file call the solvers.
  4. Link to the library 'Riemann_solver/Riemann_solver.a' when compiling your programm.


In this package we mainly implement the GRP solver for the two-dimensional Euler equations of the polytropic gases. The codes are in the file './src/system_GRP_solver.c'. The definition of the API can be also found in the file 'Riemann_solver/inc/Riemann_solver.h'.


They rely on the Riemann solver. Here we implement an exact Riemann solver in './src/Riemann_solver_exact.c'.

The classical HLL solver is also implemented in './src/HLL.c'.

GRP sovlers for some particular scalar equations are implemented in './src/scalar_GRP_solver.c'.

