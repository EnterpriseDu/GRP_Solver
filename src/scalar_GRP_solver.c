

void burgers_GRP_solver
(double D[], double U[], double lambda, double rho_L, double rho_R, double d_rho_L, double d_rho_R, double eps)
{
  double sigma = 0.5*(rho_L + rho_R);
  double v_L = rho_L - lambda, v_R = rho_R - lambda;


  if(v_R > -eps)
    if(v_L < eps)  //transonic rarefaction
    {
      U[0] = lambda;
      D[0] = 0.0;
    }
    else  //wave propagates to the right
    {
      U[0] = rho_L;
      D[0] = (lambda - U[0]) * d_rho_L;
    }
  else
    if(v_L < eps)  //wave propagates to the left
    {
      U[0] = rho_R;
      D[0] = (lambda - U[0]) * d_rho_R;
    }
    else  //shock
      if(sigma > 0.0)
      {
	U[0] = rho_L;
	D[0] = (lambda - U[0]) * d_rho_L;
      }
      else
      {
	U[0] = rho_R;
	D[0] = (lambda - U[0]) * d_rho_R;
      }
}

void advection_GRP_solver
(double D[], double U[], double a, double lambda, double rho_L, double rho_R, double d_rho_L, double d_rho_R, double eps)
{
  double sigma = 0.5*(rho_L + rho_R);
  double u = a-lambda;


  if(u > 0.0)
  {
    U[0] = rho_L;
    D[0] = -a * d_rho_L;
  }
  else
  {
    U[0] = rho_R;
    D[0] = -a * d_rho_R;
  }
}
