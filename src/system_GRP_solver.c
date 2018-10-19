/*
 * This is the GRP solver for the two-dimensional Euler equations of
 * polytropic gases:
 *                    w_t + f(w)_x + g(w)_y = 0
 * where
 *     / rho   \      / rho*u     \      / rho*v     \
 * w = | rho*u |  f = | rho*u^2+p |  g = | rho*u*v   |                  (1)
 *     | rho*v |      | rho*u*v   |      | rho*v^2+p |
 *     \ E     / ,    \ u*(E+p)   / ,    \ v*(E+p)   / ,
 * E=0.5*rho*(u^2+v^2) + p/(gamma-1).
 *
 * The outputs are:
 * wave_speed : wave_speed[0] is the speed of the left wave
 *              wave_speed[1] is the speed of the right wave
 *              they are used in some circumstances
 *
 * The inputs are:
 * /        \
 * |lambda_x| : the direction along which we solve the GRP
 * |lambda_y|   in the case of fixed meshes, lambda is 0
 * \        /
 *   gamma    : is the heat ratio for the polytropic gases
 *   eps      : is the constant epsilon
 *  n_trans   : number of transported variables subject to
 *              additional transport equations
 *   (wL,wR)  : Initial states
 *    para    : parameters
 * +--------------------------------------------------------------+
 * | The codes for calculations of derivatives of the transported |
 * | variables are to be added.                                   |
 * +--------------------------------------------------------------+
 *
 * /   rho_L \
 * |     u_L | : is the left state
 * |     v_L |
 * \     p_L /
 * / d_rho_L \
 * |   d_u_L | : is the normal derivatives of the left state
 * |   d_v_L |
 * \   d_p_L /
 * / t_rho_L \
 * |   t_u_L | : is the tangential derivatives of the left state
 * |   t_v_L |
 * \   t_p_L /
 * /   rho_R \
 * |     u_R | : is the right state
 * |     v_R |
 * \     p_R /
 * / d_rho_R \
 * |   d_u_R | : is the normal derivatives of the right state
 * |   d_v_R |
 * \   d_p_R /
 * / t_rho_R \
 * |   t_u_R | : is the tangential derivatives of the right state
 * |   t_v_R |
 * \   t_p_R /
 *
 * In the one-dimensional computations, just set tangential derivatives
 * to be zeros.
 *
 * References:
 * [1] M. Ben-Artzi and J. Falcovitz, Generalized Riemann Problems in
 *     Computational Fluid Dynamics, Cambridge University Press, Cam-
 *     bridge, 2003.
 * [2] M. Ben-Artzi, J. Li and G. Warnecke, A direct GRP scheme for
 *     compressible fluid flows, J. Comput. Phys., 218 (2006), 19-43.
 */



#include <math.h>
#include <stdio.h>

#include "Riemann_solver.h"


void linear_GRP_solver_Edir_Q1D
(double *wave_speed, EulerPack *out, double lambda_x, double lambda_y, int n_trans, EulerPack const *wL, EulerPack const *wR, RSparameters *para)
{
  double const eps = para->eps;
  double const tol = para->tol;
  int const N = para->N;
  double const geo_factor = para->geo_factor;

  double const gammaL  = wL->gamma;
  double const rho_L   = wL->VAR.rho;
  double const u_L     = wL->VAR.u;
  double const v_L     = wL->VAR.v;
  double const p_L     = wL->VAR.p;
  double const d_rho_L = wL->DER.rho;
  double const d_u_L   = wL->DER.u;
  double const d_v_L   = wL->DER.v;
  double const d_p_L   = wL->DER.p;
  double const t_rho_L = wL->TGT.rho;
  double const t_u_L   = wL->TGT.u;
  double const t_v_L   = wL->TGT.v;
  double const t_p_L   = wL->TGT.p;
  double const gammaR  = wR->gamma;
  double const rho_R   = wR->VAR.rho;
  double const u_R     = wR->VAR.u;
  double const v_R     = wR->VAR.v;
  double const p_R     = wR->VAR.p;
  double const d_rho_R = wR->DER.rho;
  double const d_u_R   = wR->DER.u;
  double const d_v_R   = wR->DER.v;
  double const d_p_R   = wR->DER.p;
  double const t_rho_R = wR->TGT.rho;
  double const t_u_R   = wR->TGT.u;
  double const t_v_R   = wR->TGT.v;
  double const t_p_R   = wR->TGT.p;


  int CRW[2]={0,0};
  double dist;
  double c_L, c_R, C, c_frac;

  double d_Phi, d_Psi, TdS, VAR;
  double D_rho, D_u, D_v, D_p, D_z, D_phi, T_rho, T_u, T_v, T_p, T_z, T_phi; 
  double u_star, p_star, rho_star_L, rho_star_R, c_star_L, c_star_R;

  double H1, H2, H3;
  double a_L, b_L, d_L, a_R, b_R, d_R, detA;
  double L_u, L_p, L_rho;

  double u_t_mat, p_t_mat;
  double shk_u_s, shk_u_L, shk_u_R;
  
  const double zetaL = (gammaL-1.0)/(gammaL+1.0);
  const double zetaR = (gammaR-1.0)/(gammaR+1.0);
 
  double rho_x;
  double g_rho, g_u, g_p, f;

  double speed_L, speed_R;

  c_L = sqrt(gammaL * p_L / rho_L);
  c_R = sqrt(gammaR * p_R / rho_R);

  dist = (u_L-u_R)*(u_L-u_R) + (p_L-p_R)*(p_L-p_R);
  //=========acoustic case==========
  Riemann_solver_exact(&u_star, &p_star, gammaL, gammaR, u_L, u_R, p_L, p_R, c_L, c_R, 1.0, CRW, eps, tol, 500);
  if(CRW[0])
    {
      rho_star_L = rho_L*pow(p_star/p_L, 1.0/gammaL);
      c_star_L = c_L*pow(p_star/p_L, 0.5*(gammaL-1.0)/gammaL);
      speed_L = u_L - c_L;
    }
  else
    {
      rho_star_L = rho_L*(p_star+zetaL*p_L)/(p_L+zetaL*p_star);
      c_star_L = sqrt(gammaL * p_star / rho_star_L);
      speed_L = u_L - c_L*sqrt(0.5*((gammaL+1.0)*(p_star/p_L) + (gammaL-1.0))/gammaL);
    }
  if(CRW[1])
    {
      rho_star_R = rho_R*pow(p_star/p_R,1.0/gammaR);
      c_star_R = c_R*pow(p_star/p_R, 0.5*(gammaR-1.0)/gammaR);
      speed_R = u_R + c_R;
    }
  else
    {
      rho_star_R = rho_R*(p_star+zetaR*p_R)/(p_R+zetaR*p_star);
      c_star_R = sqrt(gammaR * p_star / rho_star_R);
      speed_R = u_R + c_R*sqrt(0.5*((gammaR+1.0)*(p_star/p_R) + (gammaR-1.0))/gammaR);
    }
  wave_speed[0] = speed_L;
  wave_speed[1] = speed_R;


  if(dist < eps)
    {
      wave_speed[0] = speed_L;
      wave_speed[1] = speed_R;

      if(speed_L > lambda_x) //the direction is on the left side of all the three waves
	{
	  out->VAR.rho = rho_L;
	  out->VAR.u =   u_L;
	  out->VAR.v =   v_L;
	  out->VAR.p =   p_L;
	  out->DER.rho = -(u_L-lambda_x)*d_rho_L - (v_L-lambda_y)*t_rho_L - rho_L*(d_u_L+t_v_L);
	  out->DER.u = -(u_L-lambda_x)*d_u_L   - (v_L-lambda_y)*t_u_L   - d_p_L/rho_L;
	  out->DER.v = -(u_L-lambda_x)*d_v_L   - (v_L-lambda_y)*t_v_L   - t_p_L/rho_L;
	  out->DER.p = -(u_L-lambda_x)*d_p_L   - (v_L-lambda_y)*t_p_L   - rho_L*c_L*c_L*(d_u_L+t_v_L) ;
	}
      else if(speed_R < lambda_x) //the direction is on the right side of all the three waves
	{
	  out->VAR.rho = rho_R;
	  out->VAR.u =   u_R;
	  out->VAR.v =   v_R;
	  out->VAR.p =   p_R;
	  out->DER.rho = -(u_R-lambda_x)*d_rho_R - (v_R-lambda_y)*t_rho_R - rho_R*(d_u_R+t_v_R);
	  out->DER.u = -(u_R-lambda_x)*d_u_R   - (v_R-lambda_y)*t_u_R   - d_p_R/rho_R;
	  out->DER.v = -(u_R-lambda_x)*d_v_R   - (v_R-lambda_y)*t_v_R   - t_p_R/rho_R;
	  out->DER.p = -(u_R-lambda_x)*d_p_R   - (v_R-lambda_y)*t_p_R   - rho_R*c_R*c_R*(d_u_R+t_v_R);
	}
      else
	{
	  if(CRW[0] && ((u_star-c_star_L) > lambda_x)) // the direction is in a 1-CRW
	    {
	      out->VAR.u = zetaL*(u_L+2.0*(c_L+lambda_x)/(gammaL-1.0));
	      C = out->VAR.u - lambda_x;
	      out->VAR.p = pow(C/c_L, 2.0*gammaL/(gammaL-1.0)) * p_L;
	      out->VAR.rho = gammaL*out->VAR.p/C/C;
	      out->VAR.v = v_L;
	    }
	  else if(CRW[1] && ((u_star+c_star_R) < lambda_x)) // the direction is in a 3-CRW
	    {
	      out->VAR.u = zetaR*(u_R-2.0*(c_R-lambda_x)/(gammaR-1.0));
	      C = lambda_x-out->VAR.u;
	      out->VAR.p = pow(C/c_R, 2.0*gammaR/(gammaR-1.0)) * p_R;
	      out->VAR.rho = gammaR*out->VAR.p/C/C;
	      out->VAR.v = v_R;
	    }	
	  else if(u_star > lambda_x) //the direction is between the 1-wave and the contact discontinuety
	    {
	      out->VAR.rho = rho_star_L;
	      out->VAR.u =   u_star;
	      out->VAR.v =        v_L;
	      out->VAR.p =   p_star;
	      C    =   c_star_L;
	    }
	  else //the direction is between the contact discontinuety and the 3-wave
	    {
	      out->VAR.rho = rho_star_R;
	      out->VAR.u =   u_star;
	      out->VAR.v =        v_R;
	      out->VAR.p =   p_star;
	      C    =   c_star_R;
	    }			

	  D_p = 0.5*((d_u_L*(out->VAR.rho*C) + d_p_L) - (d_u_R*(out->VAR.rho*C) - d_p_R));			
	  T_p = 0.5*((t_u_L*(out->VAR.rho*C) + t_p_L) - (t_u_R*(out->VAR.rho*C) - t_p_R));
	  D_u = 0.5*(d_u_L + d_p_L/(out->VAR.rho*C) + d_u_R - d_p_R/(out->VAR.rho*C));			
	  T_u = 0.5*(t_u_L + t_p_L/(out->VAR.rho*C) + t_u_R - t_p_R/(out->VAR.rho*C));			
	  if(u_star > lambda_x)
	    {
	      D_v = d_v_L;
	      T_v = t_v_L;
	      D_rho = d_rho_L - d_p_L/(C*C) + D_p/(C*C);
	      T_rho = t_rho_L - t_p_L/(C*C) + T_p/(C*C);				
	    }
	  else
	    {
	      D_v = d_v_R;
	      T_v = t_v_R;
	      D_rho = d_rho_R - d_p_R/(C*C) + D_p/(C*C);
	      T_rho = t_rho_R - t_p_R/(C*C) + T_p/(C*C);
	    }
	  out->DER.rho = -(out->VAR.u-lambda_x)*D_rho - (out->VAR.v-lambda_y)*T_rho - out->VAR.rho*(D_u+T_v);
	  out->DER.u = -(out->VAR.u-lambda_x)*D_u   - (out->VAR.v-lambda_y)*T_u   - D_p/out->VAR.rho;
	  out->DER.v = -(out->VAR.u-lambda_x)*D_v   - (out->VAR.v-lambda_y)*T_v   - T_p/out->VAR.rho;
	  out->DER.p = -(out->VAR.u-lambda_x)*D_p   - (out->VAR.v-lambda_y)*T_p   - out->VAR.rho*C*C*(D_u+T_v);	
	}
      return;
    }

  //=========non-acoustic case==========
  //------trivial case------
  if(speed_L > lambda_x) //the direction is on the left side of all the three waves
    {
      out->VAR.rho = rho_L;
      out->VAR.u =   u_L;
      out->VAR.v =   v_L;
      out->VAR.p =   p_L;
      out->DER.rho = -(u_L-lambda_x)*d_rho_L - (v_L-lambda_y)*t_rho_L - rho_L*(d_u_L+t_v_L);
      out->DER.u = -(u_L-lambda_x)*d_u_L   - (v_L-lambda_y)*t_u_L   - d_p_L/rho_L;
      out->DER.v = -(u_L-lambda_x)*d_v_L   - (v_L-lambda_y)*t_v_L   - t_p_L/rho_L;
      out->DER.p = -(u_L-lambda_x)*d_p_L   - (v_L-lambda_y)*t_p_L   - rho_L*c_L*c_L*(d_u_L+t_v_L);
    }
  else if(speed_R < lambda_x) //the direction is on the right side of all the three waves
    {
      out->VAR.rho = rho_R;
      out->VAR.u =   u_R;
      out->VAR.v =   v_R;
      out->VAR.p =   p_R;
      out->DER.rho = -(u_R-lambda_x)*d_rho_R - (v_R-lambda_y)*t_rho_R - rho_R*(d_u_R+t_v_R);
      out->DER.u = -(u_R-lambda_x)*d_u_R   - (v_R-lambda_y)*t_u_R   - d_p_R/rho_R;
      out->DER.v = -(u_R-lambda_x)*d_v_R   - (v_R-lambda_y)*t_v_R   - t_p_R/rho_R;
      out->DER.p = -(u_R-lambda_x)*d_p_R   - (v_R-lambda_y)*t_p_R   - rho_R*c_R*c_R*(d_u_R+t_v_R);
    }
  else//----non-trivial case----
    {
      if(CRW[0] && ((u_star-c_star_L) > lambda_x)) // the direction is in a 1-CRW
	{
	  out->VAR.u = zetaL*(u_L+2.0*(c_L+lambda_x)/(gammaL-1.0));
	  C = out->VAR.u - lambda_x;
	  out->VAR.p = pow(C/c_L, 2.0*gammaL/(gammaL-1.0)) * p_L;
	  out->VAR.rho = gammaL*out->VAR.p/C/C;
	  out->VAR.v = v_L;

	  c_frac = C/c_L;
	  TdS = (d_p_L - d_rho_L*c_L*c_L)/(gammaL-1.0)/rho_L;
	  d_Psi = d_u_L + (gammaL*d_p_L/c_L - c_L*d_rho_L)/(gammaL-1.0)/rho_L;

	  out->DER.u = ((1.0+zetaL)*pow(c_frac, 0.5/zetaL) + zetaL*pow(c_frac, (1.0+zetaL)/zetaL));
	  out->DER.u = out->DER.u/(1.0+2.0*zetaL) * TdS;
	  out->DER.u = out->DER.u - c_L*pow(c_frac, 0.5/zetaL)*d_Psi;
	  out->DER.p = out->VAR.rho*(out->VAR.u - lambda_x)*out->DER.u;

	  out->DER.rho = out->VAR.rho*(out->VAR.u - lambda_x)*pow(c_frac, (1.0+zetaL)/zetaL)*TdS*(gammaL-1.0);
	  out->DER.rho = (out->DER.rho + out->DER.p) / C/C;

	  out->DER.v = -(out->VAR.u - lambda_x)*d_v_L*out->VAR.rho/rho_L;
	}
      else if(CRW[1] && ((u_star+c_star_R) < lambda_x)) // the direction is in a 3-CRW
	{
	  out->VAR.u = zetaR*(u_R-2.0*(c_R-lambda_x)/(gammaR-1.0));
	  C = lambda_x-out->VAR.u;
	  out->VAR.p = pow(C/c_R, 2.0*gammaR/(gammaR-1.0)) * p_R;
	  out->VAR.rho = gammaR*out->VAR.p/C/C;
	  out->VAR.v = v_R;
					
	  c_frac = C/c_R;
	  TdS = (d_p_R - d_rho_R*c_R*c_R)/(gammaR-1.0)/rho_R;
	  d_Phi = d_u_R - (gammaR*d_p_R/c_R - c_R*d_rho_R)/(gammaR-1.0)/rho_R;

	  out->DER.u = ((1.0+zetaR)*pow(c_frac, 0.5/zetaR) + zetaR*pow(c_frac, (1.0+zetaR)/zetaR));
	  out->DER.u = out->DER.u/(1.0+2.0*zetaR) * TdS;
	  out->DER.u = out->DER.u + c_R*pow(c_frac, 0.5/zetaR)*d_Phi;
	  out->DER.p = out->VAR.rho*(out->VAR.u-lambda_x)*out->DER.u;

	  out->DER.rho = out->VAR.rho*(out->VAR.u-lambda_x)*pow(c_frac, (1.0+zetaR)/zetaR)*TdS*(gammaR-1.0);
	  out->DER.rho = (out->DER.rho + out->DER.p) / C/C;

	  out->DER.v = -(out->VAR.u-lambda_x)*d_v_R*out->VAR.rho/rho_R;
	}
      else//--non-sonic case--
	{
	  if(u_star < lambda_x) //the direction is between the contact discontinuety and the 3-wave
	    {
	      out->VAR.rho = rho_star_R;
	      out->VAR.u =   u_star;
	      out->VAR.v =   v_R;
	      out->VAR.p =   p_star;
	      C = c_star_R;
	    }
	  else //the direction is between the 1-wave and the contact discontinuety
	    {
	      out->VAR.rho = rho_star_L;
	      out->VAR.u =   u_star;
	      out->VAR.v =   v_L;
	      out->VAR.p =   p_star;
	      C = c_star_L;
	    }

	  //determine a_L, b_L and d_L
	  if(CRW[0]) //the 1-wave is a CRW
	    {
	      a_L = 1.0;
	      b_L = 1.0 / rho_star_L / c_star_L;
	      c_frac = c_star_L/c_L;
	      TdS = (d_p_L - d_rho_L*c_L*c_L)/(gammaL-1.0)/rho_L;
	      d_Psi = d_u_L + (gammaL*d_p_L/c_L - c_L*d_rho_L)/(gammaL-1.0)/rho_L;
	      d_L = ((1.0+zetaL)*pow(c_frac, 0.5/zetaL) + zetaL*pow(c_frac, (1.0+zetaL)/zetaL));
	      d_L = d_L/(1.0+2.0*zetaL) * TdS;
	      d_L = d_L - c_L*pow(c_frac, 0.5/zetaL) * d_Psi;
	    }
	  else //the 1-wave is a shock
	    {
	      shk_u_s = -sqrt(0.5*((gammaL+1.0)*p_L   +(gammaL-1.0)*p_star)/rho_star_L);
	      shk_u_L = -sqrt(0.5*((gammaL+1.0)*p_star+(gammaL-1.0)*p_L   )/rho_L);

	      VAR = sqrt((1-zetaL)/(rho_L*(p_star+zetaL*p_L)));

	      H1 =  0.5*VAR * (p_star+(1.0+2.0*zetaL)*p_L)/(p_star+zetaL*p_L);
	      H2 = -0.5*VAR * ((2.0+zetaL)*p_star + zetaL*p_L)/(p_star+zetaL*p_L);
	      H3 = -0.5*VAR * (p_star-p_L) / rho_L;

	      L_p = -1.0/rho_L - shk_u_L*H2;
	      L_u = shk_u_L + rho_L*(c_L*c_L*H2 + H3);
	      L_rho = -shk_u_L * H3;

	      a_L = 1.0 - rho_star_L* shk_u_s * H1;
	      b_L = -shk_u_s/(rho_star_L*c_star_L*c_star_L)+ H1;
	      d_L = L_rho*d_rho_L + L_u*d_u_L + L_p*d_p_L;
	    }
	  //determine a_R, b_R and d_R
	  if(CRW[1]) //the 3-wave is a CRW
	    {
	      a_R = 1.0;
	      b_R = -1.0 / rho_star_R / c_star_R;
	      c_frac = c_star_R/c_R;
	      TdS = (d_p_R - d_rho_R*c_R*c_R)/(gammaR-1.0)/rho_R;
	      d_Phi = d_u_R - (gammaR*d_p_R/c_R - c_R*d_rho_R)/(gammaR-1.0)/rho_R;
	      d_R = ((1.0+zetaR)*pow(c_frac, 0.5/zetaR) + zetaR*pow(c_frac, (1.0+zetaR)/zetaR));
	      d_R = d_R/(1.0+2.0*zetaR) * TdS;
	      d_R = d_R + c_R*pow(c_frac, 0.5/zetaR)*d_Phi;
	    }
	  else //the 3-wave is a shock
	    {
	      shk_u_s = sqrt(0.5*((gammaR+1.0)*p_R   + (gammaR-1.0)*p_star)/rho_star_R);
	      shk_u_R = sqrt(0.5*((gammaR+1.0)*p_star+ (gammaR-1.0)*p_R   )/rho_R);

	      VAR  = sqrt((1.0-zetaR)/(rho_R*(p_star+zetaR*p_R)));

	      H1 = 0.5* VAR * (p_star+(1+2.0*zetaR)*p_R)/(p_star+zetaR*p_R);
	      H2 = -0.5*VAR * ((2.0+zetaR)*p_star+zetaR*p_R)/(p_star+zetaR*p_R);
	      H3 = -0.5*(p_star-p_R)* VAR /rho_R;

	      L_p = -1.0/rho_R + shk_u_R*H2;
	      L_u = shk_u_R - rho_R*(c_R*c_R*H2 + H3);
	      L_rho = shk_u_R * H3;

	      a_R = 1.0 +rho_star_R* shk_u_s * H1;
	      b_R = -(shk_u_s/(rho_star_R*c_star_R*c_star_R) + H1);
	      d_R = L_rho*d_rho_R + L_u*d_u_R + L_p*d_p_R;
	    }

	  detA = a_L*b_R - b_L*a_R;
	  u_t_mat = (b_R*d_L - b_L*d_R)/detA;
	  p_t_mat = (a_L*d_R - a_R*d_L)/detA;

	  //already total D!
	  out->DER.u = u_t_mat + (u_star-lambda_x)/out->VAR.rho/C/C * p_t_mat;
	  out->DER.p = p_t_mat + (u_star-lambda_x)*out->VAR.rho * u_t_mat;
	
	  if(u_star < lambda_x) //the direction is between the contact discontinuety and the 3-wave
	    {
	      if(CRW[1]) //the 3-wave is a CRW
		{
		  //already total D!
		  out->DER.rho = rho_star_R*(u_star-lambda_x)*pow(c_star_R/c_R, (1.0+zetaR)/zetaR)*(d_p_R - d_rho_R*c_R*c_R)/rho_R;
		  out->DER.rho = (out->DER.rho + out->DER.p) / c_star_R/c_star_R;

		  out->DER.v = -out->VAR.u*d_v_R*out->VAR.rho/rho_R;
		  out->DER.v = out->DER.v + lambda_x*d_v_R;
		}
	      else //the 3-wave is a shock
		{
		  shk_u_s = sqrt(0.5*((gammaR+1.0)*p_R   + (gammaR-1.0)*p_star)/rho_star_R);
		  shk_u_R = sqrt(0.5*((gammaR+1.0)*p_star+ (gammaR-1.0)*p_R   )/rho_R);

		  VAR = p_R + zetaR*p_star;
		  H1 = rho_R * p_R    * (1.0 - zetaR*zetaR) / VAR/VAR;
		  H2 = rho_R * p_star * (zetaR*zetaR - 1.0) / VAR/VAR;
		  H3 = (p_star + zetaR*p_R)/VAR;

		  L_rho = shk_u_R * H3 * d_rho_R;
		  L_u = -rho_R * (H2*c_R*c_R + H3) * d_u_R;
		  L_p = H2 * shk_u_R * d_p_R;

		  out->DER.rho = ((u_star+shk_u_s)/c_star_R/c_star_R - u_star*H1)*p_t_mat + rho_star_R*u_star*shk_u_s*H1*u_t_mat;
		  out->DER.rho = (out->DER.rho - u_star*(L_p+L_rho+L_u)) / shk_u_s;

		  f = shk_u_R*(H2*d_p_R + H3*d_rho_R) - rho_R*(H2*c_R*c_R+H3)*d_u_R;
		  rho_x = (f + H1*(p_t_mat - rho_star_R*shk_u_s*u_t_mat) - out->DER.rho) / (shk_u_R+u_R);//shk_spd;
		  //out->DER.rho = out->DER.rho + lambda_x*rho_x; (Dr. X.Lei's comments)

		  out->DER.v = -out->VAR.u * shk_u_R * d_v_R / shk_u_s;
		  out->DER.v = out->DER.v + lambda_x*d_v_R;
		}
	    }
	  else //the direction is between the 1-wave and the contact discontinuety
	    {
	      if(CRW[0]) //the 1-wave is a CRW
		{
		  //already total D!
		  out->DER.rho = rho_star_L*(u_star-lambda_x)*pow(c_star_L/c_L, (1.0+zetaL)/zetaL)*(d_p_L - d_rho_L*c_L*c_L)/rho_L;
		  out->DER.rho = (out->DER.rho + out->DER.p) / c_star_L/c_star_L;

		  out->DER.v = -out->VAR.u*d_v_L*out->VAR.rho/rho_L;
		  out->DER.v = out->DER.v + lambda_x*d_v_L;
		}
	      else //the 1-wave is a shock
		{
		  shk_u_s = -sqrt(0.5*((gammaL+1.0)*p_L   +(gammaL-1.0)*p_star)/rho_star_L);
		  shk_u_L = -sqrt(0.5*((gammaL+1.0)*p_star+(gammaL-1.0)*p_L   )/rho_L);

		  VAR = p_L + zetaL*p_star;

		  H1 = rho_L * p_L    * (1.0 - zetaL*zetaL) / VAR/VAR;
		  H2 = rho_L * p_star * (zetaL*zetaL - 1.0) / VAR/VAR;
		  H3 = (p_star + zetaL*p_L)/VAR;

		  L_rho = shk_u_L * H3 * d_rho_L;
		  L_u = -rho_L*(H2*c_L*c_L + H3) * d_u_L;
		  L_p = H2 * shk_u_L * d_p_L;

		  out->DER.rho = ((u_star+shk_u_s)/c_star_L/c_star_L - H1*u_star)*p_t_mat + rho_star_L*u_star*shk_u_s*H1*u_t_mat;
		  out->DER.rho = (out->DER.rho - u_star*(L_p+L_rho+L_u))/ shk_u_s;

		  f = shk_u_L*(H2*d_p_L + H3*d_rho_L) - rho_L*(H2*c_L*c_L+H3)*d_u_L;
		  rho_x = (f + H1*(p_t_mat - rho_star_L*shk_u_s*u_t_mat) - out->DER.rho) / (shk_u_L+u_L);
		  //out->DER.rho = out->DER.rho + lambda_x*rho_x; (Dr. X.Lei's comments)

		  out->DER.v = -out->VAR.u * shk_u_L * d_v_L / shk_u_s;
		  out->DER.v = out->DER.v + lambda_x*d_v_L;
		}
	    }
	  //--end of non-sonic case--
	}
      T_p = 0.5*((t_u_L*(out->VAR.rho*C) + t_p_L) - (t_u_R*(out->VAR.rho*C) - t_p_R));
      T_u = 0.5*(t_u_L + t_p_L/(out->VAR.rho*C) + t_u_R - t_p_R/(out->VAR.rho*C));
      if (u_star > lambda_x)
	{
	  T_rho = t_rho_L - t_p_L/(C*C) + T_p/(C*C);
	  out->DER.rho = out->DER.rho - (out->VAR.v-lambda_y)*T_rho - out->VAR.rho*t_v_L;
	  out->DER.u = out->DER.u - (out->VAR.v-lambda_y)*T_u;
	  out->DER.v = out->DER.v - (out->VAR.v-lambda_y)*t_v_L - T_p/out->VAR.rho;
	  out->DER.p = out->DER.p - (out->VAR.v-lambda_y)*T_p   - out->VAR.rho*C*C*t_v_L;							
	}
      else
	{
	  T_rho = t_rho_R - t_p_R/(C*C) + T_p/(C*C);
	  out->DER.rho = out->DER.rho - (out->VAR.v-lambda_y)*T_rho - out->VAR.rho*t_v_R;
	  out->DER.u = out->DER.u - (out->VAR.v-lambda_y)*T_u;
	  out->DER.v = out->DER.v - (out->VAR.v-lambda_y)*t_v_R - T_p/out->VAR.rho;
	  out->DER.p = out->DER.p - (out->VAR.v-lambda_y)*T_p   - out->VAR.rho*C*C*t_v_R;
	}
      //----end of non-trivial case----
    }
}




void Euler_Godunov_solver(double wave_speed[2], EulerPack *out, double lambda_x, double lambda_y, int n_trans, EulerPack const *wL, EulerPack const *wR, RSparameters *para)
{
  double const eps = para->eps;
  double const tol = para->tol;
  int const N = para->N;

  double const gammaL  = wL->gamma;
  double const rho_L   = wL->VAR.rho;
  double const u_L     = wL->VAR.u;
  double const v_L     = wL->VAR.v;
  double const p_L     = wL->VAR.p;
  double const gammaR  = wR->gamma;
  double const rho_R   = wR->VAR.rho;
  double const u_R     = wR->VAR.u;
  double const v_R     = wR->VAR.v;
  double const p_R     = wR->VAR.p;

  double const c_L = sqrt(gammaL * p_L / rho_L);
  double const c_R = sqrt(gammaR * p_R / rho_R);

  double const gammaL_1 = gammaL-1.0;
  double const gammaR_1 = gammaR-1.0;
  double const _over_gammaL_1 = 1.0/gammaL_1;
  double const _over_gammaR_1 = 1.0/gammaR_1;
  //double const nuL = 1.0/(gammaL-1.0);
  //double const nuR = 1.0/(gammaR-1.0);
  double const muL = (gammaL-1.0)/(gammaL+1.0);
  double const muR = (gammaR-1.0)/(gammaR+1.0);
  double const muL_2 = muL*muL;
  double const muR_2 = muR*muR;



  double dist;
  double c_frac, cricit, c_star, c_square;
  int CRW[2];
  double speed_L, speed_R;
  double d_Phi, d_Psi, TdS, VAR, RIEL, RIER, PHIxL, PHIxR, PHIyL, PHIyR;
  double u_star, p_star, rho_star_L, rho_star_R, c_star_L, c_star_R;

  double PI, H1, H2, H3;
  double a_L, b_L, d_L, a_R, b_R, d_R, detA;
  double L_u, L_p, L_rho;

  double u_t_mat, p_t_mat, D0_u_tau, D0_p_tau;
  double T_rho, T_u, T_v, T_p;

  double shk_spd, shk_u_s, shk_u_L, shk_u_R;
  //double C, rho_x;
  double g_rho, g_u, g_p, f;

  int it_trans;




  cricit= (u_R - 2.0*_over_gammaR_1*c_R) - (u_L + 2.0*_over_gammaL_1*c_L);
  if((rho_L < eps) || (p_L < eps) || (rho_R < eps) || (p_R < eps) || (cricit > -eps))
  {
    vacuum(wave_speed, out, lambda_x, lambda_y, n_trans, wL, wR, para);
    return;
  }



  //dist = (rho_L-rho_R)*(rho_L-rho_R) + (u_L-u_R)*(u_L-u_R) + (p_L-p_R)*(p_L-p_R);
  dist = (u_L-u_R)*(u_L-u_R) + (p_L-p_R)*(p_L-p_R);
//=========acoustic case==========
  if(dist < eps)
  {
    wave_speed[0] = u_L-c_L;
    wave_speed[1] = u_R+c_R;

    rho_star_L = rho_L;
    rho_star_R = rho_R;
      c_star_L =   c_L;
      c_star_R =   c_R;
      u_star = 0.5*(u_R+u_L);
      p_star = 0.5*(p_R+p_L);

    if(lambda_x < u_star)
    {
      c_star = c_L;
      c_square = c_star*c_star;
      out->gamma   = gammaL;
      out->VAR.rho = rho_star_L;
      out->VAR.u   =   u_star;
      out->VAR.v   =   v_L;
      out->VAR.p   =   p_star;
      for(it_trans = 0; it_trans < n_trans; ++it_trans)
	out->VAR.trans[it_trans] = wL->VAR.trans[it_trans];
    }
    else
    {
      c_star = c_R;
      c_square = c_star*c_star;
      out->gamma   = gammaR;
      out->VAR.rho = rho_star_R;
      out->VAR.u   =   u_star;
      out->VAR.v   =   v_R;
      out->VAR.p   =   p_star;
      for(it_trans = 0; it_trans < n_trans; ++it_trans)
	out->VAR.trans[it_trans] = wR->VAR.trans[it_trans];
    }
    return;
  }


  Riemann_solver_exact(&u_star, &p_star, gammaL, gammaR, u_L, u_R, p_L, p_R, c_L, c_R, 1.0, CRW, eps, tol, N);
  if(fabs(u_star) < eps)
    u_star = 0.0;
  if(CRW[0]){
    rho_star_L = rho_L*pow(p_star/p_L, 1.0/gammaL);
      c_star_L =   c_L*pow(p_star/p_L, 0.5*(gammaL-1.0)/gammaL);
       speed_L =   u_L - c_L;
  }  else{
    rho_star_L = rho_L*(p_star+muL*p_L)/(p_L+muL*p_star);
      c_star_L = sqrt(gammaL * p_star / rho_star_L);
       speed_L = u_L - c_L*sqrt(0.5*((gammaL+1.0)*(p_star/p_L) + (gammaL-1.0))/gammaL);
  }
  if(CRW[1]){
    rho_star_R = rho_R*pow(p_star/p_R,1.0/gammaR);
      c_star_R =   c_R*pow(p_star/p_R, 0.5*(gammaR-1.0)/gammaR);
       speed_R =   u_R + c_R;
  }  else{
    rho_star_R = rho_R*(p_star+muR*p_R)/(p_R+muR*p_star);
      c_star_R = sqrt(gammaR * p_star / rho_star_R);
       speed_R = u_R + c_R*sqrt(0.5*((gammaR+1.0)*(p_star/p_R) + (gammaR-1.0))/gammaR);
  }
  wave_speed[0] = speed_L;
  wave_speed[1] = speed_R;

  //------trivial case------
  if(speed_L > lambda_x) //the direction is on the left side of all the three waves
  {
    c_star = c_L;
    c_square = c_star*c_star;

    out->gamma   = gammaL;
    out->VAR.rho = rho_L;
    out->VAR.u   =   u_L;
    out->VAR.v   =   v_L;
    out->VAR.p   =   p_L;

    for(it_trans = 0; it_trans < n_trans; ++it_trans)
      out->VAR.trans[it_trans] = wL->VAR.trans[it_trans];
  }
  else if(speed_R < lambda_x) //the direction is on the right side of all the three waves
  {
    c_star = c_R;
    c_square = c_star*c_star;

    out->gamma   = gammaR;
    out->VAR.rho = rho_R;
    out->VAR.u   = u_R;
    out->VAR.v   = v_R;
    out->VAR.p   = p_R;

    for(it_trans = 0; it_trans < n_trans; ++it_trans)
      out->VAR.trans[it_trans] = wR->VAR.trans[it_trans];
  }
  else//----non-trivial case----
  {
    if((CRW[0]) && ((u_star-c_star_L) > lambda_x)) // the direction is in a 1-CRW
    {
      out->gamma   = gammaL;
      //out->VAR.u   = muL*(u_L+2.0*(c_L+lambda_x)/(gammaL-1.0));
      out->VAR.u   = muL*(u_L+2.0*_over_gammaL_1*(c_L+lambda_x));
      c_star = out->VAR.u - lambda_x;
      c_square = c_star*c_star;
      out->VAR.p   = pow(c_star/c_L, 2.0*gammaL*_over_gammaL_1) * p_L;
      out->VAR.rho = gammaL*out->VAR.p / c_square;
      out->VAR.v   = v_L;

      for(it_trans = 0; it_trans < n_trans; ++it_trans)
	out->VAR.trans[it_trans] = wL->VAR.trans[it_trans];
    }
    else if((CRW[1]) && ((u_star+c_star_R) < lambda_x)) // the direction is in a 3-CRW
    {
      out->gamma   = gammaR;
      //out->VAR.u   = muR*(u_R-2.0*(c_R-lambda)/(gammaR-1.0));
      out->VAR.u   = muR*(u_R-2.0*_over_gammaR_1*(c_R-lambda_x));
      c_star = lambda_x - out->VAR.u;
      c_square = c_star*c_star;
      out->VAR.p   = pow(c_star/c_R, 2.0*gammaR*_over_gammaR_1) * p_R;
      out->VAR.rho = gammaR*out->VAR.p  /c_square;
      out->VAR.v   = v_R;

      for(it_trans = 0; it_trans < n_trans; ++it_trans)
	out->VAR.trans[it_trans] = wR->VAR.trans[it_trans];
    }
    else//--non-sonic case--
    {
      if(u_star < lambda_x) //the direction is between the contact discontinuety and the 3-wave
      {
	out->gamma   = gammaR;
	out->VAR.rho = rho_star_R;
	out->VAR.u   =   u_star;
	out->VAR.v   =   v_R;
	out->VAR.p   =   p_star;
	for(it_trans = 0; it_trans < n_trans; ++it_trans)
	  out->VAR.trans[it_trans] = wR->VAR.trans[it_trans];
      }
      else //the direction is between the 1-wave and the contact discontinuety
      {
	out->gamma   = gammaL;
	out->VAR.rho = rho_star_L;
	out->VAR.u   =   u_star;
	out->VAR.v   =   v_L;
	out->VAR.p   =   p_star;
	for(it_trans = 0; it_trans < n_trans; ++it_trans)
	  out->VAR.trans[it_trans] = wL->VAR.trans[it_trans];
      }
    }
  }
}
