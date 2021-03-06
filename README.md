A Brief Instruction to *GRP_Solver* -- an implementation of the GRP solver
===

A quick start
---

  1. Cast to the directory 'GRP_solver/'.
  2. Run 'make'.
  3. Include the header 'GRP_solver/inc/Riemann_solver.h' in the source file where the solvers are called.
  4. Link to the archive 'GRP_solver/GRP_Solver.a' when compiling your programm.



Introduction
---

In this package we mainly implement the GRP solver for the two-dimensional Euler equations of the polytropic gases. The codes are in the file './src/system_GRP_solver.c'. The definition of the API can be also found in the file 'GRP_solver/inc/Riemann_solver.h'.


They rely on the Riemann solver. Here we implement an exact Riemann solver in './src/Riemann_solver_exact.c'.

The classical HLL solver is also implemented in './src/HLL.c'.

GRP sovlers for some particular scalar equations are implemented in './src/scalar_GRP_solver.c'.



Interface of the GRP solver 'Euler_GRP_solver'
---

Since the GRP solver is the key part of the present program package, we briefly introduce the interface of 'linear_GRP_solver_Edir_Q1D' in 'src/system_Euler_solver.c'.

Firstly, please refer to the mathematical instruction provided in the package and references therein for basic information of generalized Riemann problem of hyperbolic conservation laws.

The output
-

The first two parameters of 'linear_GRP_solver_Edir_Q1D' are the output of the GRP solver.

1. 'wave_speed' contains two double variable which are the speeds of the left- and right-most waves, respectively. For example, if the left wave is a rarefaction wave, wave_speed[0] is $u_L-c_L$. If the left wave is a shock, wave_speed[0] is the shock speed.
2. 'out' contains the Riemann solution and the instantaneous directional derivative of the solution at the initial discontinuity. The Riemann solution and the instantaneous derivative are stored in 'out.VAR' and 'out.DER', respectively. The GRP solver can calculate the instantaneous directional derivative along an arbitrary direction. This parameter is a 'EulerPack' object. 'EulerPack' is a structure defined in 'inc/Riemann_solver.h' containing the fluid state. Details of 'EulerPack' will be given below.

The input
-

The following parameters are the input of the GRP solver.
1. The third and fourth parameters 'lambda_x' and 'lambda_y' defines the direction of the directional derivative. For a fixed cell interface, they are 0.
2. 'n_trans' will be introduced after the definition of 'EulerPack' is given in the next subsection.
3. The sixth and seventh parameter 'wL' and 'wR' are the initial fluid states on the left and right sides of the initial discontinuity, respectively. As 'out', they are also of the type 'EulerPack', the definition of which will be given below.
4. The last parameter 'para' contains some configuration parameters:
  **a**. 'eps' is epsilon;
  **b**. 'tol' is the tolerance for the exact Rieman solver;
  **c**. 'N' is the maximam number of iterations in the Riemann solver;
  **d**. 'geo_factor' is the geometric factor of the quasi 1-D flow, whose value should be 0 for planar flow;
  **e**. 'radius' and 'nDim' are currenly abandoned.

The structure 'EulerPack'
-

The definition of 'EulerPack' is given in 'inc/Riemann_solver.h'. It is designed to contain the fluid state.
1. The structure member 'gamma' is the heat ratio of the fluid.
2. Structure members 'VAR', 'DER' and 'TGT' contains the reconstructed limiting values of the fluid state and its derivatives along the normal and tangential directions of the initial discontinuity, respectively. They are of the type 'EulerVar', the definition of which will be given later.

The structure 'EulerVar' contains four *double* variables and an array of *double*.

The four *double* variables 'rho', 'u', 'v' and 'p' represent the density, normal and tangential components of the velocity and the pressure, respectively. So, 'wL.VAR.rho' is the limiting value of the density on the left side of the initial discontinuity and 'wR.TGT.p' is the limiting value of the tangential derivative of the pressure on the right side of the initial discontinuity.

The array 'trans' contains the additional variables following transport equations in addition to the original Euler equations. They can either follow

$\phi_t+(u,v)\cdot\nabla\phi=0$, or $(\rho\phi)_t+\nabla\cdot(\rho(u,v)\phi)=0$.

The length of the array 'trans' in 'EulerVAR' is given by the fifth parameter 'n_trans'.

**Remark:** The heat ratios in 'wL' and 'wR' are not necessarily the same. This allows the current version of the GRP solver to deal with multi-component fluid with varying heat ratio.

**Remark:** The Riemann solver 'Euler_Godunov_solver' has the same interface as that of 'linear_GRP_solver_Edir_Q1D'.

Numerical experiments
---

The non-linear GRP solver (by disabling the acoustic case) is tested by the transported of the sine wave.

|  m  | L_1 error | L_1 order | L_\infty error | L_\infty order |
|:---:|:---------:|:---------:|:--------------:|:--------------:|
|  20 |  9.85e-6  |           |     3.86e-5    |                |
|  40 |  2.86e-6  |    1.79   |     1.30e-5    |       1.58     |
|  80 |  7.52e-7  |    1.93   |     3.46e-6    |       1.90     |
| 160 |  1.92e-7  |    1.97   |     8.87e-7    |       1.97     |
| 320 |  4.86e-8  |    1.98   |     2.24e-7    |       1.99     |
| 640 |  1.22e-8  |    1.99   |     5.60e-8    |       2.00     |
|1280 |  3.06e-9  |    2.00   |     1.39e-8    |       2.00     |
|2560 |  7.64e-10 |    2.00   |     3.14e-9    |       2.03     |

The acoustic case codes are also tested by the same experiments and the results are almost the same which is reasonable since the flow is smooth.

Applied in the two-stage fourth-order GRP scheme, the numerical results for the same problme are listed in the folowing table.

|  m  | L_1 error | L_1 order | L_\infty error | L_\infty order |
|:---:|:---------:|:---------:|:--------------:|:--------------:|
|  20 |  4.35e-5  |           |     7.01e-5    |                |
|  40 |  1.37e-6  |   5.00    |     2.15e-6    |        5.03    |
|  80 |  4.28e-8  |   5.00    |     6.72e-8    |        5.00    |
| 160 |  1.34e-9  |   5.00    |     2.11e-9    |        5.00    |
| 320 |  4.26e-11 |   4.98    |     6.69e-10   |        4.98    |
| 640 |  1.40e-12 |   4.93    |     2.22e-11   |        4.91    |
|1280 |  5.06e-13 |   4.79    |     1.13e-13   |        4.29    |

