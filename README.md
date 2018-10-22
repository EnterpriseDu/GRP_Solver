A Brief Instruction
===

How to use
---

  1. Enter the directory 'GRP_solver/'.
  2. Run 'make'.
  3. Include the header 'GRP_solver/inc/Riemann_solver.h' in the file call the solvers.
  4. Link to the library 'GRP_solver/GRP_Solver.a' when compiling your programm.

introduction
---

In this package we mainly implement the GRP solver for the two-dimensional Euler equations of the polytropic gases. The codes are in the file './src/system_GRP_solver.c'. The definition of the API can be also found in the file 'GRP_solver/inc/Riemann_solver.h'.


They rely on the Riemann solver. Here we implement an exact Riemann solver in './src/Riemann_solver_exact.c'.

The classical HLL solver is also implemented in './src/HLL.c'.

GRP sovlers for some particular scalar equations are implemented in './src/scalar_GRP_solver.c'.


Numerical Experiments
---

The non-linear GRP solver (by disabling the acoustic case) is tested by the transported of the sine wave and the isentropic flow.
The numerical errors and orders of the latter problem are

|  m  | L_1 error | L_1 order | L_\infty error | L_\infty order |
|:---:|:---------:|:---------:|:--------------:|:--------------:|
|  20 |  9.85e-5  |           |     3.86e-4    |                |
|  40 |  2.86e-5  |    1.79   |     1.30e-4    |       1.58     |
|  80 |  7.52e-2  |    1.93   |     3.46e-5    |       1.90     |
| 160 |  1.93e-5  |    1.96   |     8.92e-6    |       1.96     |
| 320 |  4.86e-6  |    1.99   |     2.23e-6    |       2.00     |
| 640 |  1.23e-6  |    1.98   |     6.03e-7    |       1.89     |
|1280 |  3.14e-7  |    1.97   |     1.87e-7    |       1.69     |

The acoustic case codes are also tested by the same experiments and the results are almost the same which is reasonable since the flow is smooth.