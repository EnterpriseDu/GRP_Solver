#include <math.h>
#include <stdio.h>

#include "Riemann_solver.h"



void ADER2_solver_P2P
(double wave_speed[2], double D[4], double U[4], double lambda, double gamma, double eps, double tol,
 double rho_L, double u_L, double v_L, double p_L,
 double rho_R, double u_R, double v_R, double p_R,
 double d_rho_L, double d_u_L, double d_v_L, double d_p_L,
 double d_rho_R, double d_u_R, double d_v_R, double d_p_R,
 double t_rho_L, double t_u_L, double t_v_L, double t_p_L,
 double t_rho_R, double t_u_R, double t_v_R, double t_p_R)
{
  double c_L, c_R, cricit;
  int CRW[2];
  double dist, PHIxL, PHIxR, PHIyL, PHIyR;
  double u_star, p_star, rho_star_L, rho_star_R, c_star_L, c_star_R, c_star, c_square;

  double speed_L, speed_R, zeta = (gamma-1.0)/(gamma+1.0);


  c_L = sqrt(gamma * p_L / rho_L);
  c_R = sqrt(gamma * p_R / rho_R);

  cricit= u_R - u_L - 2.0*(c_R+c_L)/(gamma-1.0);
  if((rho_L < eps) || (p_L < eps) || (rho_R < eps) || (p_R < eps) || (cricit > -eps))
  {
    vacuum_P2P(wave_speed, D, U, cricit, c_L, c_R, lambda, gamma, eps,
	       rho_L, u_L, v_L, p_L, rho_R, u_R, v_R, p_R,
	       d_rho_L, d_u_L, d_v_L, d_p_L, d_rho_R, d_u_R, d_v_R, d_p_R,
	       t_rho_L, t_u_L, t_v_L, t_p_L, t_rho_R, t_u_R, t_v_R, t_p_R);
    return;
  }

  double g0, g1, g2, g3, g4;


  dist = sqrt((rho_L-rho_R)*(rho_L-rho_R) + (u_L-u_R)*(u_L-u_R) + (v_L-v_R)*(v_L-v_R) + (p_L-p_R)*(p_L-p_R));
//=========Solve the Riemann problem==========
  if(dist < eps)
  {
    wave_speed[0] = u_L-c_L;
    wave_speed[1] = u_R+c_R;
    CRW[0] = 0; CRW[1] = 0;
    u_star = 0.5*(u_R+u_L);
    if(u_star > 0.0)
      p_star = p_L;
    else
      p_star = p_R;
    c_star_L = c_L;
    c_star_R = c_R;
    rho_star_L = rho_L;
    rho_star_R = rho_R;
  }
  else
  {
    Riemann_solver_exact(&u_star, &p_star, gamma, gamma, u_L, u_R, p_L, p_R, c_L, c_R, 1.0, CRW, eps, tol, 500);

    if(CRW[0]){
      rho_star_L = rho_L*pow(p_star/p_L, 1.0/gamma);
        c_star_L =   c_L*pow(p_star/p_L, 0.5*(gamma-1.0)/gamma);
      wave_speed[0] =   u_L - c_L;
    }  else{
      rho_star_L = rho_L*(p_star+zeta*p_L)/(p_L+zeta*p_star);
        c_star_L = sqrt(gamma * p_star / rho_star_L);
      wave_speed[0] = u_L - c_L*sqrt(0.5*((gamma+1.0)*(p_star/p_L) + (gamma-1.0))/gamma);
    }
    if(CRW[1]){
      rho_star_R = rho_R*pow(p_star/p_R,1.0/gamma);
        c_star_R =   c_R*pow(p_star/p_R, 0.5*(gamma-1.0)/gamma);
      wave_speed[1] =   u_R + c_R;
    }  else{
      rho_star_R = rho_R*(p_star+zeta*p_R)/(p_R+zeta*p_star);
        c_star_R = sqrt(gamma * p_star / rho_star_R);
      wave_speed[1] = u_R + c_R*sqrt(0.5*((gamma+1.0)*(p_star/p_R) + (gamma-1.0))/gamma);
    }
  }

//=========Compute the state in the star region==========
  if(u_star > lambda)
  {
    U[2] = v_L;
    if(wave_speed[0] > lambda)
    {
      U[0] = rho_L;
      U[1] = u_L;
      U[3] = p_L;
      c_star = c_L;
      c_square = c_star*c_star;
    }
    else if(CRW[0] && (u_star-c_star_L>lambda))
    {
      U[1] = zeta*(u_L+2.0*(c_L+lambda)/(gamma-1.0));
      c_star = U[1] - lambda;
      c_square = c_star*c_star;
      U[3] = pow(c_star/c_L, 2.0*gamma/(gamma-1.0)) * p_L;
      U[0] = gamma*U[3]/c_square;
    }
    else
    {
      U[0] = rho_star_L;
      U[1] = u_star;
      U[3] = p_star;
      c_star = c_star_L;
      c_square = c_star*c_star;
    }
  }
  else
  {
    U[2] = v_R;
    if(wave_speed[1] < lambda)
    {
      U[0] = rho_R;
      U[1] = u_R;
      U[3] = p_R;
      c_star = c_R;
      c_square = c_star*c_star;
    }
    else if(CRW[0] && (u_star+c_star_L)<lambda)
    {
      U[1] = zeta*(u_R-2.0*(c_R-lambda)/(gamma-1.0));
      c_star = lambda-U[1];
      c_square = c_star*c_star;
      U[3] = pow(c_star/c_R, 2.0*gamma/(gamma-1.0)) * p_R;
      U[0] = gamma*U[3]/c_square;
    }
    else
    {
      U[0] = rho_star_R;
      U[1] = u_star;
      U[3] = p_star;
      c_star = c_star_R;
      c_square = c_star*c_star;
    }
  }

//=========Compute the direvatives in the star region==========
  if(U[1] > lambda)
  {
    if(U[1]-c_star > lambda)
    {
      D[0] = -U[1]*d_rho_L - U[0]*d_u_L;
      D[1] = -U[1]*  d_u_L - d_p_L/U[0];
      D[2] = -U[1]*  d_v_L;
      D[3] = -U[1]*  d_p_L - U[0]*d_u_L*c_square;

      g0 = U[2]*t_rho_L + U[0]*t_v_L;
      g1 = U[2]*  t_u_L;
      g2 = U[2]*  t_v_L + t_p_L/U[0];
      g3 = U[2]*  t_p_L + U[0]*t_v_L*c_square;
    }
    else
    {
      PHIxL = d_u_L + d_p_L/(U[0]*c_star);
      PHIxR = d_u_R - d_p_R/(U[0]*c_star);
      PHIyL = t_u_L + t_p_L/(U[0]*c_star);
      PHIyR = t_u_R - t_p_R/(U[0]*c_star);

      D[1] = -0.5*((U[1]+c_star)*PHIxL + (U[1]-c_star)*PHIxR);
      D[3] = -0.5*((U[1]+c_star)*PHIxL - (U[1]-c_star)*PHIxR)*U[0]*c_star;
      D[2] = -U[1]*d_v_L;
      D[0] = (D[3] + U[1]*(d_p_L - c_square*d_rho_L))/c_square;

      g1 = 0.5*U[2]*(PHIyL + PHIyR);
      g3 = 0.5*(U[2]*(PHIyL - PHIyR) + 2.0*c_star*t_v_L)*U[0]*c_star;
      g2 = U[2]*t_v_L + 0.5*c_star*(PHIyL - PHIyR);
      g0 = (g3 - U[2]*(t_p_L - c_square*t_rho_L))/c_square;
      /*
      D[1] = -0.5*(((U[1]+c_star)*PHIxL + (U[1]-c_star)*PHIxR) + U[2]*(PHIyL + PHIyR));
      D[3] = -0.5*(((U[1]+c_star)*PHIxL - (U[1]-c_star)*PHIxR) + U[2]*(PHIyL - PHIyR) + 2.0*c_star*t_v_L)*U[0]*c_star;
      D[2] = -U[1]*d_v_L - U[2]*t_v_L - 0.5*c_star*(PHIyL - PHIyR);
      D[0] = (D[3] + U[1]*(d_p_L - c_square*d_rho_L) + U[2]*(t_p_L - c_square*t_rho_L))/c_square;
      //*/
    }
  }
  else
  {
    if(U[1]+c_star < lambda)
    {
      D[0] = -U[1]*d_rho_R - U[0]*d_u_R;
      D[1] = -U[1]*  d_u_R - d_p_R/U[0];
      D[2] = -U[1]*  d_v_R;
      D[3] = -U[1]*  d_p_R - U[0]*d_u_R*c_square;

      g0 = U[2]*t_rho_R + U[0]*t_v_R;
      g1 = U[2]*  t_u_R;
      g2 = U[2]*  t_v_R + t_p_R/U[0];
      g3 = U[2]*  t_p_R + U[0]*t_v_R*c_square;
    }
    else
    {
      PHIxL = d_u_L + d_p_L/(U[0]*c_star);
      PHIxR = d_u_R - d_p_R/(U[0]*c_star);
      PHIyL = t_u_L + t_p_L/(U[0]*c_star);
      PHIyR = t_u_R - t_p_R/(U[0]*c_star);

      D[1] = -0.5*((U[1]+c_star)*PHIxL + (U[1]-c_star)*PHIxR);
      D[3] = -0.5*((U[1]+c_star)*PHIxL - (U[1]-c_star)*PHIxR)*U[0]*c_star;
      D[2] = -U[1]*d_v_R;
      D[0] = (D[3] + U[1]*(d_p_R - c_square*d_rho_R))/c_square;

      g1 = 0.5*U[2]*(PHIyL + PHIyR);
      g3 = 0.5*(U[2]*(PHIyL - PHIyR) + 2.0*c_star*t_v_R)*U[0]*c_star;
      g2 = U[2]*t_v_R + 0.5*c_star*(PHIyL - PHIyR);
      g0 = (g3 - U[2]*(t_p_R - c_square*t_rho_R))/c_square;
      /*
      D[1] = -0.5*(((U[1]+c_star)*PHIxL + (U[1]-c_star)*PHIxR) + U[2]*(PHIyL + PHIyR));
      D[3] = -0.5*(((U[1]+c_star)*PHIxL - (U[1]-c_star)*PHIxR) + U[2]*(PHIyL - PHIyR) + 2.0*c_star*t_v_R)*U[0]*c_star;
      D[2] = -U[1]*d_v_R - U[2]*t_v_R - 0.5*c_star*(PHIyL - PHIyR);
      D[0] = (D[3] + U[1]*(d_p_R - c_square*d_rho_R) + U[2]*(t_p_R - c_square*t_rho_R))/c_square;
      //*/
    }
  }

  D[0] = D[0] - g0;
  D[1] = D[1] - g1;
  D[2] = D[2] - g2;
  D[3] = D[3] - g3;

  /*
  if(U[1] > lambda)
  {
    if(U[1]-c_star > lambda)
    {
      D[0] = -U[1]*d_rho_L - U[2]*t_rho_L - U[0]*(d_u_L+t_v_L);
      D[1] = -U[1]*  d_u_L - U[2]*  t_u_L - d_p_L/U[0];
      D[2] = -U[1]*  d_v_L - U[2]*  t_v_L - t_p_L/U[0];
      D[3] = -U[1]*  d_p_L - U[2]*  t_p_L - U[0]*(d_u_L+t_v_L)*c_square;
    }
    else
    {
      PHIxL = d_u_L + d_p_L/(U[0]*c_star);
      PHIxR = d_u_R - d_p_R/(U[0]*c_star);
      PHIyL = t_u_L + t_p_L/(U[0]*c_star);
      PHIyR = t_u_R - t_p_R/(U[0]*c_star);
      D[1] = -0.5*(((U[1]+c_star)*PHIxL + (U[1]-c_star)*PHIxR) + U[2]*(PHIyL + PHIyR));
      D[3] = -0.5*(((U[1]+c_star)*PHIxL - (U[1]-c_star)*PHIxR) + U[2]*(PHIyL - PHIyR) + 2.0*c_star*t_v_L)*U[0]*c_star;
      D[2] = -U[1]*d_v_L - U[2]*t_v_L - 0.5*c_star*(PHIyL - PHIyR);
      D[0] = (D[3] + U[1]*(d_p_L - c_square*d_rho_L) + U[2]*(t_p_L - c_square*t_rho_L))/c_square;
    }
  }
  else
  {
    if(U[1]+c_star < lambda)
    {
      D[0] = -U[1]*d_rho_R - U[2]*t_rho_R - U[0]*(d_u_R+t_v_R);
      D[1] = -U[1]*  d_u_R - U[2]*  t_u_R - d_p_R/U[0];
      D[2] = -U[1]*  d_v_R - U[2]*  t_v_R - t_p_R/U[0];
      D[3] = -U[1]*  d_p_R - U[2]*  t_p_R - U[0]*(d_u_R+t_v_R)*c_square;
    }
    else
    {
      PHIxL = d_u_L + d_p_L/(U[0]*c_star);
      PHIxR = d_u_R - d_p_R/(U[0]*c_star);
      PHIyL = t_u_L + t_p_L/(U[0]*c_star);
      PHIyR = t_u_R - t_p_R/(U[0]*c_star);
      D[1] = -0.5*(((U[1]+c_star)*PHIxL + (U[1]-c_star)*PHIxR) + U[2]*(PHIyL + PHIyR));
      D[3] = -0.5*(((U[1]+c_star)*PHIxL - (U[1]-c_star)*PHIxR) + U[2]*(PHIyL - PHIyR) + 2.0*c_star*t_v_R)*U[0]*c_star;
      D[2] = -U[1]*d_v_R - U[2]*t_v_R - 0.5*c_star*(PHIyL - PHIyR);
      D[0] = (D[3] + U[1]*(d_p_R - c_square*d_rho_R) + U[2]*(t_p_R - c_square*t_rho_R))/c_square;
    }
  }
  */
}



void ADER2_solver_C2P
(double wave_speed[2], double D[4], double U[4], double lambda, double gamma, double eps, double tol,
 double rho_L, double u_L, double v_L, double p_L,
 double rho_R, double u_R, double v_R, double p_R,
 double d_rho_L, double d_u_L, double d_v_L, double d_p_L,
 double d_rho_R, double d_u_R, double d_v_R, double d_p_R,
 double sL1, double sL2, double sL3, double sL4,
 double sR1, double sR2, double sR3, double sR4)
{
  double c_L, c_R, cricit;
  int CRW[2];
  double dist, PHIxL, PHIxR, PHIyL, PHIyR;
  double u_star, p_star, rho_star_L, rho_star_R, c_star_L, c_star_R, c_star, c_square;

  double speed_L, speed_R, zeta = (gamma-1.0)/(gamma+1.0);



  c_L = sqrt(gamma * p_L / rho_L);
  c_R = sqrt(gamma * p_R / rho_R);
  cricit= u_R - u_L - 2.0*(c_R+c_L)/(gamma-1.0);
  if((rho_L < eps) || (p_L < eps) || (rho_R < eps) || (p_R < eps) || (cricit > -eps))
  {
    vacuum_C2P(wave_speed, D, U, cricit, c_L, c_R, lambda, gamma, eps,
	       rho_L, u_L, v_L, p_L, rho_R, u_R, v_R, p_R,
	       d_rho_L, d_u_L, d_v_L, d_p_L, d_rho_R, d_u_R, d_v_R, d_p_R,
	       sL1, sL2, sL3, sL4, sR1, sR2, sR3, sR4);
    return;
  }



  double g00, g0, g1, g2, g3, g4, sD1, sD2, sD3, sD4, PROD, V, H_star;


  dist = sqrt((rho_L-rho_R)*(rho_L-rho_R) + (u_L-u_R)*(u_L-u_R) + (v_L-v_R)*(v_L-v_R) + (p_L-p_R)*(p_L-p_R));
//=========Solve the Riemann problem==========
  if(dist < eps)
  {
    wave_speed[0] = u_L-c_L;
    wave_speed[1] = u_R+c_R;
    CRW[0] = 0; CRW[1] = 0;
    u_star = 0.5*(u_R+u_L);
    if(u_star > 0.0)
      p_star = p_L;
    else
      p_star = p_R;
    c_star_L = c_L;
    c_star_R = c_R;
    rho_star_L = rho_L;
    rho_star_R = rho_R;
  }
  else
  {
    Riemann_solver_exact(&u_star, &p_star, gamma, gamma, u_L, u_R, p_L, p_R, c_L, c_R, 1.0, CRW, eps, tol, 500);

    if(CRW[0]){
      rho_star_L = rho_L*pow(p_star/p_L, 1.0/gamma);
        c_star_L =   c_L*pow(p_star/p_L, 0.5*(gamma-1.0)/gamma);
      wave_speed[0] =   u_L - c_L;
    }  else{
      rho_star_L = rho_L*(p_star+zeta*p_L)/(p_L+zeta*p_star);
        c_star_L = sqrt(gamma * p_star / rho_star_L);
      wave_speed[0] = u_L - c_L*sqrt(0.5*((gamma+1.0)*(p_star/p_L) + (gamma-1.0))/gamma);
    }
    if(CRW[1]){
      rho_star_R = rho_R*pow(p_star/p_R,1.0/gamma);
        c_star_R =   c_R*pow(p_star/p_R, 0.5*(gamma-1.0)/gamma);
      wave_speed[1] =   u_R + c_R;
    }  else{
      rho_star_R = rho_R*(p_star+zeta*p_R)/(p_R+zeta*p_star);
        c_star_R = sqrt(gamma * p_star / rho_star_R);
      wave_speed[1] = u_R + c_R*sqrt(0.5*((gamma+1.0)*(p_star/p_R) + (gamma-1.0))/gamma);
    }
  }

//=========Compute the state in the star region==========
  if(u_star > lambda)
  {
    U[2] = v_L;
    if(wave_speed[0] > lambda)
    {
      U[0] = rho_L;
      U[1] = u_L;
      U[3] = p_L;
      c_star = c_L;
      c_square = c_star*c_star;
    }
    else if(CRW[0] && (u_star-c_star_L>lambda))
    {
      U[1] = zeta*(u_L+2.0*(c_L+lambda)/(gamma-1.0));
      c_star = U[1] - lambda;
      c_square = c_star*c_star;
      U[3] = pow(c_star/c_L, 2.0*gamma/(gamma-1.0)) * p_L;
      U[0] = gamma*U[3]/c_square;
    }
    else
    {
      U[0] = rho_star_L;
      U[1] = u_star;
      U[3] = p_star;
      c_star = c_star_L;
      c_square = c_star*c_star;
    }
  }
  else
  {
    U[2] = v_R;
    if(wave_speed[1] < lambda)
    {
      U[0] = rho_R;
      U[1] = u_R;
      U[3] = p_R;
      c_star = c_R;
      c_square = c_star*c_star;
    }
    else if(CRW[0] && (u_star+c_star_L)<lambda)
    {
      U[1] = zeta*(u_R-2.0*(c_R-lambda)/(gamma-1.0));
      c_star = lambda-U[1];
      c_square = c_star*c_star;
      U[3] = pow(c_star/c_R, 2.0*gamma/(gamma-1.0)) * p_R;
      U[0] = gamma*U[3]/c_square;
    }
    else
    {
      U[0] = rho_star_R;
      U[1] = u_star;
      U[3] = p_star;
      c_star = c_star_R;
      c_square = c_star*c_star;
    }
  }

//=========Compute the direvatives in the star region==========
  V = U[1]*U[1] + U[2]*U[2];
  if(U[1] > lambda)
  {
    if(U[1]-c_star > lambda)
    {
      D[0] = -U[1]*d_rho_L - U[0]*d_u_L;
      D[1] = -U[1]*  d_u_L - d_p_L/U[0];
      D[2] = -U[1]*  d_v_L;
      D[3] = -U[1]*  d_p_L - U[0]*d_u_L*c_square;
    }
    else
    {
      PHIxL = d_u_L + d_p_L/(U[0]*c_star);
      PHIxR = d_u_R - d_p_R/(U[0]*c_star);
      D[1] = -0.5*((U[1]+c_star)*PHIxL + (U[1]-c_star)*PHIxR);
      D[3] = -0.5*((U[1]+c_star)*PHIxL - (U[1]-c_star)*PHIxR)*U[0]*c_star;
      D[2] = -U[1]*d_v_L;
      D[0] = (D[3] + U[1]*(d_p_L - c_square*d_rho_L))/c_square;


      sD1 = sR1 - sL1;
      sD2 = sR2 - sL2;
      sD3 = sR3 - sL3;
      sD4 = sR4 - sL4;
      H_star = 0.5*V + c_square/(gamma-1.0);
      PROD = sD1*0.5*V - sD2*U[1] - sD3*U[2] + sD4;
      PROD = PROD*(gamma-1.0) + c_star*(U[1]*sD1 - sD2);
      PROD = 0.5*PROD/c_square;
      sL1 = sL1 + PROD;
      sL2 = sL2 + PROD*(U[1]-c_star);
      sL3 = sL3 + PROD*U[2];
      sL4 = sL4 + PROD*(H_star-U[1]*c_star);
    }
    g1 = sL3;
    g2 = U[2]*(sL2-U[1]*sL1)/U[0];
    g0 = 0.5*V*sL1 - U[1]*sL2 - U[2]*sL3 + sL4;
    g00 = U[2]*sL1-sL3;
    g3 = ((gamma-1.0)*g0 - U[2]*g00) / U[0];
    g4 = (gamma-1.0)*U[2]*g0 - c_square*g00;
    D[0] = D[0] - g1;
    D[1] = D[1] - g2;
    D[2] = D[2] - g3;
    D[3] = D[3] - g4;
  }
  else
  {
    if(U[1]+c_star < lambda)
    {
      D[0] = -U[1]*d_rho_R - U[0]*d_u_R;
      D[1] = -U[1]*  d_u_R - d_p_R/U[0];
      D[2] = -U[1]*  d_v_R;
      D[3] = -U[1]*  d_p_R - U[0]*d_u_R*c_square;
    }
    else
    {
      PHIxL = d_u_L + d_p_L/(U[0]*c_star);
      PHIxR = d_u_R - d_p_R/(U[0]*c_star);
      D[1] = -0.5*((U[1]+c_star)*PHIxL + (U[1]-c_star)*PHIxR);
      D[3] = -0.5*((U[1]+c_star)*PHIxL - (U[1]-c_star)*PHIxR)*U[0]*c_star;
      D[2] = -U[1]*d_v_R;
      D[0] = (D[3] + U[1]*(d_p_R - c_square*d_rho_R))/c_square;


      sD1 = sL1 - sR1;
      sD2 = sL2 - sR2;
      sD3 = sL3 - sR3;
      sD4 = sL4 - sR4;
      H_star = 0.5*V + c_square/(gamma-1.0);
      PROD = sD1*0.5*V - sD2*U[1] - sD3*U[2] + sD4;
      PROD = PROD*(gamma-1.0) - c_star*(U[1]*sD1 - sD2);
      PROD = 0.5*PROD/c_square;
      sR1 = sR1 + PROD;
      sR2 = sR2 + PROD*(U[1]+c_star);
      sR3 = sR3 + PROD*U[2];
      sR4 = sR4 + PROD*(H_star+U[1]*c_star);
    }
    g1 = sR3;
    g2 = U[2]*(sR2-U[1]*sR1)/U[0];
    g0 = 0.5*V*sR1 - U[1]*sR2 - U[2]*sR3 + sR4;
    g00 = U[2]*sR1-sR3;
    g3 = ((gamma-1.0)*g0 - U[2]*g00) / U[0];
    g4 = (gamma-1.0)*U[2]*g0 - c_square*g00;
    D[0] = D[0] - g1;
    D[1] = D[1] - g2;
    D[2] = D[2] - g3;
    D[3] = D[3] - g4;
  }
}




void ADER2_solver_Bpre
(double wave_speed[2], double D[4], double U[4], double lambda, double gamma, double eps, double tol,
 double rho_L, double u_L, double v_L, double p_L,
 double rho_R, double u_R, double v_R, double p_R,
 double d_rho_L, double d_u_L, double d_v_L, double d_p_L,
 double d_rho_R, double d_u_R, double d_v_R, double d_p_R,
 double t_rho_L, double t_u_L, double t_v_L, double t_p_L,
 double t_rho_R, double t_u_R, double t_v_R, double t_p_R)
{
  double c_L, c_R, cricit;
  int CRW[2];
  double dist, PHIxL, PHIxR;
  double u_star, p_star, rho_star_L, rho_star_R, c_star_L, c_star_R, c_star, c_square;
  double rho_star, RIEL, RIER, v_star;
  
  double speed_L, speed_R, zeta = (gamma-1.0)/(gamma+1.0);


  double hL1, hL2, hL3, hL4, hR1, hR2, hR3, hR4, hD1, hD2, hD3, hD4, PROD, V, H_star;

  hL1 = t_rho_L*v_L + rho_L*t_v_L;
  hL2 = hL1*u_L + rho_L*v_L*t_u_L;
  hL3 = hL1*v_L + rho_L*v_L*t_v_L + t_p_L;
  hL4 = (0.5*t_rho_L*(u_L*u_L+v_L*v_L) + rho_L*(u_L*t_u_L+v_L*t_v_L) + t_p_L) * v_L;
  hL4 = hL4 + (0.5*rho_L*(u_L*u_L+v_L*v_L)+p_L) * t_v_L;

  hR1 = t_rho_R*v_R + rho_R*t_v_R;
  hR2 = hR1*u_R + rho_R*v_R*t_u_R;
  hR3 = hR1*v_R + rho_R*v_R*t_v_R + t_p_R;
  hR4 = (0.5*t_rho_R*(u_R*u_R+v_R*v_R) + rho_R*(u_R*t_u_R+v_R*t_v_R) + t_p_R) * v_R;
  hR4 = hR4 + (0.5*rho_R*(u_R*u_R+v_R*v_R)+p_R) * t_v_R;

  c_L = sqrt(gamma * p_L / rho_L);
  c_R = sqrt(gamma * p_R / rho_R);
  cricit= u_R - u_L - 2.0*(c_R+c_L)/(gamma-1.0);
  if((rho_L < eps) || (p_L < eps) || (rho_R < eps) || (p_R < eps) || (cricit > -eps))
  {
    vacuum_Bpre(wave_speed, D, U, cricit, c_L, c_R, lambda, gamma, eps,
		rho_L, u_L, v_L, p_L, rho_R, u_R, v_R, p_R,
		d_rho_L, d_u_L, d_v_L, d_p_L, d_rho_R, d_u_R, d_v_R, d_p_R,
		hL1, hL2, hL3, hL4, hR1, hR2, hR3, hR4);
    return;
  }


  
  dist = sqrt((rho_L-rho_R)*(rho_L-rho_R) + (u_L-u_R)*(u_L-u_R) + (v_L-v_R)*(v_L-v_R) + (p_L-p_R)*(p_L-p_R));

  
//=========Solving the Riemann problem==========
  if(dist < eps)
  {
    wave_speed[0] = u_L-c_L;
    wave_speed[1] = u_R+c_R;
    CRW[0] = 0; CRW[1] = 0;
    u_star = 0.5*(u_R+u_L);
    if(u_star > 0.0)
      p_star = p_L;
    else
      p_star = p_R;
    c_star_L = c_L;
    c_star_R = c_R;
    rho_star_L = rho_L;
    rho_star_R = rho_R;
  }
  else
  {
    Riemann_solver_exact(&u_star, &p_star, gamma, gamma, u_L, u_R, p_L, p_R, c_L, c_R, 1.0, CRW, eps, tol, 500);

    if(CRW[0]){
      rho_star_L = rho_L*pow(p_star/p_L, 1.0/gamma);
        c_star_L =   c_L*pow(p_star/p_L, 0.5*(gamma-1.0)/gamma);
      wave_speed[0] =   u_L - c_L;
    }  else{
      rho_star_L = rho_L*(p_star+zeta*p_L)/(p_L+zeta*p_star);
        c_star_L = sqrt(gamma * p_star / rho_star_L);
      wave_speed[0] = u_L - c_L*sqrt(0.5*((gamma+1.0)*(p_star/p_L) + (gamma-1.0))/gamma);
    }
    if(CRW[1]){
      rho_star_R = rho_R*pow(p_star/p_R,1.0/gamma);
        c_star_R =   c_R*pow(p_star/p_R, 0.5*(gamma-1.0)/gamma);
      wave_speed[1] =   u_R + c_R;
    }  else{
      rho_star_R = rho_R*(p_star+zeta*p_R)/(p_R+zeta*p_star);
        c_star_R = sqrt(gamma * p_star / rho_star_R);
      wave_speed[1] = u_R + c_R*sqrt(0.5*((gamma+1.0)*(p_star/p_R) + (gamma-1.0))/gamma);
    }
  }


//=========Computing the state in the star region==========
  if(u_star > lambda)
  {
    U[2] = v_L;
    if(wave_speed[0] > lambda)
    {
      U[0] = rho_L;
      U[1] = u_L;
      U[3] = p_L;
      c_star = c_L;
      c_square = c_star*c_star;
    }
    else if(CRW[0] && (u_star-c_star_L>lambda))
    {
      U[1] = zeta*(u_L+2.0*(c_L+lambda)/(gamma-1.0));
      c_star = U[1] - lambda;
      c_square = c_star*c_star;
      U[3] = pow(c_star/c_L, 2.0*gamma/(gamma-1.0)) * p_L;
      U[0] = gamma*U[3]/c_square;
    }
    else
    {
      U[0] = rho_star_L;
      U[1] = u_star;
      U[3] = p_star;
      c_star = c_star_L;
      c_square = c_star*c_star;
    }
  }
  else
  {
    U[2] = v_R;
    if(wave_speed[1] < lambda)
    {
      U[0] = rho_R;
      U[1] = u_R;
      U[3] = p_R;
      c_star = c_R;
      c_square = c_star*c_star;
    }
    else if(CRW[0] && (u_star+c_star_L)<lambda)
    {
      U[1] = zeta*(u_R-2.0*(c_R-lambda)/(gamma-1.0));
      c_star = lambda-U[1];
      c_square = c_star*c_star;
      U[3] = pow(c_star/c_R, 2.0*gamma/(gamma-1.0)) * p_R;
      U[0] = gamma*U[3]/c_square;
    }
    else
    {
      U[0] = rho_star_R;
      U[1] = u_star;
      U[3] = p_star;
      c_star = c_star_R;
      c_square = c_star*c_star;
    }
  }



//=========Computing direvatives in the star region==========
  V = U[1]*U[1] + U[2]*U[2];
  if(U[1] > lambda)
  {
    if(U[1]-c_star > lambda)
    {
      D[0] = -U[1]*d_rho_L - U[0]*d_u_L;
      D[1] = -U[1]*  d_u_L - d_p_L/U[0];
      D[2] = -U[1]*  d_v_L;
      D[3] = -U[1]*  d_p_L - U[0]*d_u_L*c_square;
    }
    else
    {
      PHIxL = d_u_L + d_p_L/(U[0]*c_star);
      PHIxR = d_u_R - d_p_R/(U[0]*c_star);
      D[1] = -0.5*((U[1]+c_star)*PHIxL + (U[1]-c_star)*PHIxR);
      D[3] = -0.5*((U[1]+c_star)*PHIxL - (U[1]-c_star)*PHIxR)*U[0]*c_star;
      D[2] = -U[1]*d_v_L;
      D[0] = (D[3] + U[1]*(d_p_L - c_square*d_rho_L))/c_square;

      
      hD1 = hR1 - hL1;
      hD2 = hR2 - hL2;
      hD3 = hR3 - hL3;
      hD4 = hR4 - hL4;
      H_star = 0.5*V + c_square/(gamma-1.0);
      PROD = hD1*0.5*V - hD2*U[1] - hD3*U[2] + hD4;
      PROD = PROD*(gamma-1.0) + c_star*(U[1]*hD1 - hD2);
      PROD = 0.5*PROD/c_square;
      hL1 = hL1 + PROD;
      hL2 = hL2 + PROD*(U[1]-c_star);
      hL3 = hL3 + PROD*U[2];
      hL4 = hL4 + PROD*(H_star-U[1]*c_star);
    }
    
    hL2 = (hL2 - U[1]*hL1)/U[0];
    hL3 = (hL3 - U[2]*hL1)/U[0];
    hL4 = (gamma-1.0)*(hL4 - (U[1]*hL2+U[2]*hL3) + 0.5*hL1*V);
    D[0] = D[0] - hL1;
    D[1] = D[1] - hL2;
    D[2] = D[2] - hL3;
    D[3] = D[3] - hL4;
  }
  else
  {
    if(U[1]+c_star < lambda)
    {
      D[0] = -U[1]*d_rho_R - U[0]*d_u_R;
      D[1] = -U[1]*  d_u_R - d_p_R/U[0];
      D[2] = -U[1]*  d_v_R;
      D[3] = -U[1]*  d_p_R - U[0]*d_u_R*c_square;
    }
    else
    {
      PHIxL = d_u_L + d_p_L/(U[0]*c_star);
      PHIxR = d_u_R - d_p_R/(U[0]*c_star);
      D[1] = -0.5*((U[1]+c_star)*PHIxL + (U[1]-c_star)*PHIxR);
      D[3] = -0.5*((U[1]+c_star)*PHIxL - (U[1]-c_star)*PHIxR)*U[0]*c_star;
      D[2] = -U[1]*d_v_R;
      D[0] = (D[3] + U[1]*(d_p_R - c_square*d_rho_R))/c_square;
      
      hD1 = hL1 - hR1;
      hD2 = hL2 - hR2;
      hD3 = hL3 - hR3;
      hD4 = hL4 - hR4;
      H_star = 0.5*V + c_square/(gamma-1.0);
      PROD = hD1*0.5*V - hD2*U[1] - hD3*U[2] + hD4;
      PROD = PROD*(gamma-1.0) - c_star*(U[1]*hD1 - hD2);
      PROD = 0.5*PROD/c_square;
      hR1 = hR1 + PROD;
      hR2 = hR2 + PROD*(U[1]+c_star);
      hR3 = hR3 + PROD*U[1];
      hR4 = hR4 + PROD*(H_star+U[1]*c_star);
    }
  
    hR2 = (hR2 - U[1]*hR1)/U[0];
    hR3 = (hR3 - U[2]*hR1)/U[0];
    hR4 = (gamma-1.0)*(hR4 - (U[1]*hR2+U[2]*hR3) + 0.5*hR1*V);
    D[0] = D[0] - hR1;
    D[1] = D[1] - hR2;
    D[2] = D[2] - hR3;
    D[3] = D[3] - hR4;
  }


}
