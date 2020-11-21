#define N_EULER 5

typedef struct{
  double gamma, rho, u, v, p;
  double d_gamma, d_rho, d_u, d_v, d_p;
  double t_gamma, t_rho, t_u, t_v, t_p;
} RSboundary;

typedef struct{
  double rho, u, v, p; //primitive variables
  double * trans;//transported variables
} EulerVar;
typedef struct{
  double gamma, d_gamma;
  EulerVar VAR, DER, TGT;
} EulerPack;

typedef struct{
  double eps, tol;
  int N;
  double radius;
  int nDim;
  double geo_factor;
} RSparameters;

void linear_GRP_solver_Edir_Q1D
(double *wave_speed, EulerPack *out, double lambda_x, double lambda_y, int n_trans, EulerPack const *wL, EulerPack const *wR, RSparameters *para);
void Euler_GRP_solver
(double wave_speed[2], EulerPack *out, double lambda_x, double lambda_y, int n_trans, EulerPack const *wL, EulerPack const *wR, RSparameters *para);

/*
 * GRP solvers for the burgers equation and the advection equation
 */
void burgers_GRP_solver
(double D[], double U[], double lambda, double rho_L, double rho_R, double d_rho_L, double d_rho_R, double eps);
void advection_GRP_solver
(double D[], double U[], double a, double lambda, double rho_L, double rho_R, double d_rho_L, double d_rho_R, double eps);


/*
 * the exact Riemann solver for the Euler equaitons of polytropic gases
 */
double Riemann_solver_exact(double * U_star, double * P_star, double gammaL, double gammaR, double U_l, double U_r, double P_l, double P_r, double c_l, double c_r, double dist, int * CRW, double eps, double tol, int N);

/*
 * some linearized Riemann solvers
 */
void HLL(double wave_speed[2], double F[5], double U[5], double lambda, RSboundary *wL, RSboundary *wR, RSparameters *para);

void HLL_consv(double wave_speed[2], double F[5], double U[5], double lambda, RSboundary *wL, RSboundary *wR, RSparameters *para);

/*
 *GRP and ADER solvers
 *They rely on the Riemann solver declared above.
 */
void Euler_GRP_solver(double wave_speed[2], EulerPack *out, double lambda_x, double lambda_y, int n_trans, EulerPack const *wL, EulerPack const *wR, RSparameters *para);

void vacuum(double wave_speed[2], EulerPack *out, double lambda_x, double lambda_y, int n_trans, EulerPack const *wL, EulerPack const *wR, RSparameters *para);

void vacuum_Godunov(double wave_speed[2], EulerPack *out, double lambda_x, double lambda_y, int n_trans, EulerPack const *wL, EulerPack const *wR, RSparameters *para);


void Euler_GRP_left(double wave_speed[2], EulerPack *res, double lambda_x, double lambda_y, double u_star, double acc, double dudy, EulerPack *wL, RSparameters *para);

void Euler_Godunov_left(double wave_speed[2], EulerPack *res, double u_star, EulerPack *wL, RSparameters *para);

void vacuum_left(double wave_speed[2], EulerPack *res, double u_star, double acc, EulerPack *wL, RSparameters *para);




void copy_EulerVar(EulerVar *dest, EulerVar *source, int n);
void copy_EulerPack(EulerPack *dest, EulerPack *source, int n_trans);
