#include <math.h>
#include <stdio.h>

#include "Riemann_solver.h"


void Euler_GRP_solver(double wave_speed[2], EulerPack *out, double lambda_x, double lambda_y, int n_trans, EulerPack const *wL, EulerPack const *wR, RSparameters *para)
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
  double d_Phi, d_Psi, TdS, VAR, RIEL, RIER, PHIxR, PHIyR, PHIxL, PHIyL, ETAxL, ETAxR, ETAyL, ETAyR;
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


  //double t_U_s1, t_U_s2, t_U_s3, t_U_s4, t_gm_star;
  //double t_U_D1, t_U_D2, t_U_D3, t_U_D4;
  //double g1, g2, g3, g4, g5, temp1, temp2, V, H_star, PROD;
  /* double const t_U_L1 = t_rho_L; */
  /* double const t_U_L2 = t_rho_L*u_L + rho_L*t_u_L; */
  /* double const t_U_L3 = t_rho_L*v_L + rho_L*t_v_L; */
  /* double const t_U_L4 = 0.5*(t_U_L2*u_L + rho_L*u_L*t_u_L + t_U_L3*v_L + rho_L*v_L*t_v_L) + t_p_L*_over_gammaR_1; */
  /* double const t_U_R1 = t_rho_R; */
  /* double const t_U_R2 = t_rho_R*u_R + rho_R*t_u_R; */
  /* double const t_U_R3 = t_rho_R*v_R + rho_R*t_v_R; */
  /* double const t_U_R4 = 0.5*(t_U_R2*u_R + rho_R*u_R*t_u_R + t_U_R3*v_R + rho_R*v_R*t_v_R) + t_p_R*_over_gammaR_1; */
  /* double t_U_s1, t_U_s2, t_U_s3, t_U_s4, t_gm_star; */
  /* double t_U_D1, t_U_D2, t_U_D3, t_U_D4; */





  cricit= (u_R - 2.0*_over_gammaR_1*c_R) - (u_L + 2.0*_over_gammaL_1*c_L);
  if((rho_L < eps) || (p_L < eps) || (rho_R < eps) || (p_R < eps) || (cricit > -eps))
  {
    vacuum(wave_speed, out, lambda_x, lambda_y, n_trans, wL, wR, para);
    return;
  }



  Riemann_solver_exact(&u_star, &p_star, gammaL, gammaR, u_L, u_R, p_L, p_R, c_L, c_R, 1.0, CRW, eps, tol, 500);
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



  //dist = (rho_L-rho_R)*(rho_L-rho_R) + (u_L-u_R)*(u_L-u_R) + (p_L-p_R)*(p_L-p_R);
  dist = (u_L-u_R)*(u_L-u_R) + (p_L-p_R)*(p_L-p_R);
//=========acoustic case==========
  if(0)//(dist < eps)
  {
    wave_speed[0] = u_L-c_L;
    wave_speed[1] = u_R+c_R;

    PHIxR = d_u_R - d_p_R/(rho_R*c_R);
    PHIxL = d_u_L + d_p_L/(rho_L*c_L);
    PHIyR = t_u_R - t_p_R/(rho_R*c_R);
    PHIyL = t_u_L + t_p_L/(rho_L*c_L);

    ETAxL = d_rho_L - d_p_L/(c_L*c_L);
    ETAxR = d_rho_R - d_p_R/(c_R*c_R);
    ETAyL = t_rho_L - t_p_L/(c_L*c_L);
    ETAyR = t_rho_R - t_p_R/(c_R*c_R);

    if(lambda_x < u_star)
    {
      if(u_star-c_L < lambda_x)
      {
	c_star = c_star_L;
	c_square = c_star*c_star;
	out->gamma   = gammaL;
	out->VAR.rho = rho_star_L;
	out->VAR.u   =   u_star;
	out->VAR.v   =   v_L;
	out->VAR.p   =   p_star;

	out->DER.u = -0.5*( (out->VAR.u-c_star-lambda_x)*PHIxR + (out->VAR.u+c_star-lambda_x)*PHIxL + (out->VAR.v-lambda_y)*( PHIyR+PHIyL));
	out->DER.p = -0.5*(-(out->VAR.u-c_star-lambda_x)*PHIxR + (out->VAR.u+c_star-lambda_x)*PHIxL + (out->VAR.v-lambda_y)*(-PHIyR+PHIyL) + 2.0*c_star*t_v_L) * out->VAR.rho*c_star;
	out->DER.p   = out->DER.p - geo_factor*out->VAR.rho*out->VAR.u*c_square;/* geometric effect */
	out->DER.rho = -(out->VAR.u-lambda_x)*ETAxL - (out->VAR.v-lambda_y)*ETAyL + out->DER.p/c_square;
	out->DER.v = -(out->VAR.u-lambda_x)*d_v_L - (out->VAR.v-lambda_y)*t_v_L - t_p_L/out->VAR.rho;


	/* out->DER.u   = -0.5 * (((u_L+c_L)*PHIxL + (u_R-c_R)*PHIxR) + v_L*(PHIyL + PHIyR)); */
	/* out->DER.p   = -0.5 * (((u_L+c_L)*PHIxL - (u_R-c_R)*PHIxR) + v_L*(PHIyL - PHIyR) + 2.0*c_L*t_v_L) * rho_L*c_L; */
	/* out->DER.rho = (out->DER.p + out->DER.u*(d_p_L - c_square*d_rho_L) + out->DER.v*(t_p_L - c_square*t_rho_L)) / c_square; */
	/* out->DER.v   = -u_star*d_v_L - v_L*t_v_L - 0.5*c_L*(PHIyL - PHIyR); */
	/* out->DER.rho += lambda_x * d_rho_L + lambda_y * t_rho_L; */
	/* out->DER.u   += lambda_x * d_u_L   + lambda_y * t_u_L; */
	/* out->DER.v   += lambda_x * d_v_L   + lambda_y * t_v_L; */
	/* out->DER.p   += lambda_x * d_p_L   + lambda_y * t_p_L; */
      }
      else
      {
	c_star = c_L;
	c_square = c_star*c_star;
	out->gamma   = gammaL;
	out->VAR.rho = rho_L;
	out->VAR.u   =   u_L;
	out->VAR.v   =   v_L;
	out->VAR.p   =   p_L;

	out->DER.rho = (lambda_x-out->VAR.u) *d_rho_L + (lambda_y-out->VAR.v) *t_rho_L - out->VAR.rho*(d_u_L+t_v_L);
	out->DER.u   = (lambda_x-out->VAR.u) *  d_u_L + (lambda_y-out->VAR.v) *  t_u_L - d_p_L/out->VAR.rho;
	out->DER.v   = (lambda_x-out->VAR.u) *  d_v_L + (lambda_y-out->VAR.v) *  t_v_L - t_p_L/out->VAR.rho;
	out->DER.p   = (lambda_x-out->VAR.u) *  d_p_L + (lambda_y-out->VAR.v) *  t_p_L - out->VAR.rho*(d_u_L+t_v_L)*c_square;
	/* geometric effect */
	out->DER.p   = out->DER.p   - geo_factor*out->VAR.rho*out->VAR.u*c_square;
	out->DER.rho = out->DER.rho - geo_factor*out->VAR.rho*out->VAR.u;
      }
      for(it_trans = 0; it_trans < n_trans; ++it_trans)
      {
	out->VAR.trans[it_trans] = wL->VAR.trans[it_trans];
	out->DER.trans[it_trans] = (lambda_x-wL->VAR.u)*wL->DER.trans[it_trans]*(out->VAR.rho/rho_L) + (lambda_y-wL->VAR.v)*wL->TGT.trans[it_trans];
      }
    }
    else
    {
      if(lambda_x < u_star+c_R)
      {
	c_star = c_star_R;
	c_square = c_star*c_star;
	out->gamma   = gammaR;
	out->VAR.rho = rho_star_R;
	out->VAR.u   =   u_star;
	out->VAR.v   =   v_R;
	out->VAR.p   =   p_star;

	out->DER.u = -0.5*( (out->VAR.u-c_star-lambda_x)*PHIxR + (out->VAR.u+c_star-lambda_x)*PHIxL + (out->VAR.v-lambda_y)*( PHIyR+PHIyL));
	out->DER.p = -0.5*(-(out->VAR.u-c_star-lambda_x)*PHIxR + (out->VAR.u+c_star-lambda_x)*PHIxL + (out->VAR.v-lambda_y)*(-PHIyR+PHIyL) + 2.0*c_star*t_v_R) * out->VAR.rho*c_star;
	out->DER.p   = out->DER.p - geo_factor*out->VAR.rho*out->VAR.u*c_square;/* geometric effect */
	out->DER.rho = -(out->VAR.u-lambda_x)*ETAxR - (out->VAR.v-lambda_y)*ETAyR + out->DER.p/c_square;
	out->DER.v = -(out->VAR.u-lambda_x)*d_v_R - (out->VAR.v-lambda_y)*t_v_R - t_p_R/out->VAR.rho;

	/* out->DER.u   = -0.5 * (((u_L+c_L)*PHIxL + (u_R-c_R)*PHIxR) + v_R*(PHIyL + PHIyR)); */
	/* out->DER.p   = -0.5 * (((u_L+c_L)*PHIxL - (u_R-c_R)*PHIxR) + v_R*(PHIyL - PHIyR) + 2.0*c_R*t_v_R) * rho_R*c_R; */
	/* out->DER.rho = (out->DER.p   + out->DER.u*(d_p_R - c_square*d_rho_R) + out->DER.v*(t_p_R - c_square*t_rho_R)) / c_square; */
	/* out->DER.v   = -u_star*d_v_R - v_R*t_v_R - 0.5*c_R*(PHIyL - PHIyR); */
	/* //D[4] = -u_R*d_gm_R - v_R*t_gm_R; */
	/* out->DER.rho += lambda_x * d_rho_R + lambda_y * t_rho_R; */
	/* out->DER.u   += lambda_x * d_u_R   + lambda_y * t_u_R; */
	/* out->DER.v   += lambda_x * d_v_R   + lambda_y * t_v_R; */
	/* out->DER.p   += lambda_x * d_p_R   + lambda_y * t_p_R; */
      }
      else
      {
	c_star = c_R;
	c_square = c_star*c_star;
	out->gamma   = gammaR;
	out->VAR.rho = rho_R;
	out->VAR.u   =   u_R;
	out->VAR.v   =   v_R;
	out->VAR.p   =   p_R;

	out->DER.rho = (lambda_x-out->VAR.u) *d_rho_R + (lambda_y-out->VAR.v) *t_rho_R - out->VAR.rho*(d_u_R+t_v_R);
	out->DER.u   = (lambda_x-out->VAR.u) *  d_u_R + (lambda_y-out->VAR.v) *  t_u_R - d_p_R/out->VAR.rho;
	out->DER.v   = (lambda_x-out->VAR.u) *  d_v_R + (lambda_y-out->VAR.v) *  t_v_R - t_p_R/out->VAR.rho;
	out->DER.p   = (lambda_x-out->VAR.u) *  d_p_R + (lambda_y-out->VAR.v) *  t_p_R - out->VAR.rho*(d_u_R+t_v_R)*c_square;
	/* geometric effect */
	out->DER.p   = out->DER.p   - geo_factor*out->VAR.rho*out->VAR.u*c_square;
	out->DER.rho = out->DER.rho - geo_factor*out->VAR.rho*out->VAR.u;
      }
      for(it_trans = 0; it_trans < n_trans; ++it_trans)
      {
	out->VAR.trans[it_trans] = wR->VAR.trans[it_trans];
	out->DER.trans[it_trans] = (lambda_x-wR->VAR.u)*wR->DER.trans[it_trans]*(out->VAR.rho/rho_R) + (lambda_y-wR->VAR.v)*wR->TGT.trans[it_trans];
      }
    }
    return;
  }



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

    out->DER.rho = (lambda_x-out->VAR.u) *d_rho_L + (lambda_y-out->VAR.v) *t_rho_L - out->VAR.rho*(d_u_L+t_v_L);
    out->DER.u   = (lambda_x-out->VAR.u) *  d_u_L + (lambda_y-out->VAR.v) *  t_u_L - d_p_L/out->VAR.rho;
    out->DER.v   = (lambda_x-out->VAR.u) *  d_v_L + (lambda_y-out->VAR.v) *  t_v_L - t_p_L/out->VAR.rho;
    out->DER.p   = (lambda_x-out->VAR.u) *  d_p_L + (lambda_y-out->VAR.v) *  t_p_L - out->VAR.rho*(d_u_L+t_v_L)*c_square;
    /* geometric effect */
    out->DER.p   = out->DER.p   - geo_factor*out->VAR.rho*out->VAR.u*c_square;
    out->DER.rho = out->DER.rho - geo_factor*out->VAR.rho*out->VAR.u;

    for(it_trans = 0; it_trans < n_trans; ++it_trans)
    {
      out->VAR.trans[it_trans] = wL->VAR.trans[it_trans];
      out->DER.trans[it_trans] = (lambda_x-wL->VAR.u) * wL->DER.trans[it_trans] + (lambda_y-wL->VAR.v) * wL->TGT.trans[it_trans];
    }
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

    out->DER.rho = (lambda_x-out->VAR.u) *d_rho_R + (lambda_y-out->VAR.v) *t_rho_R - out->VAR.rho*(d_u_R+t_v_R);
    out->DER.u   = (lambda_x-out->VAR.u) *  d_u_R + (lambda_y-out->VAR.v) *  t_u_R - d_p_R/out->VAR.rho;
    out->DER.v   = (lambda_x-out->VAR.u) *  d_v_R + (lambda_y-out->VAR.v) *  t_v_R - t_p_R/out->VAR.rho;
    out->DER.p   = (lambda_x-out->VAR.u) *  d_p_R + (lambda_y-out->VAR.v) *  t_p_R - out->VAR.rho*(d_u_R+t_v_R)*c_square;
    /* geometric effect */
    out->DER.p   = out->DER.p   - geo_factor*out->VAR.rho*out->VAR.u*c_square;
    out->DER.rho = out->DER.rho - geo_factor*out->VAR.rho*out->VAR.u;

    for(it_trans = 0; it_trans < n_trans; ++it_trans)
    {
      out->VAR.trans[it_trans] = wR->VAR.trans[it_trans];
      out->DER.trans[it_trans] = (lambda_x-wR->VAR.u) * wR->DER.trans[it_trans] + (lambda_y-wR->VAR.v) * wR->TGT.trans[it_trans];
    }
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

      c_frac = c_star/c_L;
      //TdS = (d_p_L - d_rho_L*c_L*c_L)/(gammaL-1.0)/rho_L;
      //d_Psi = d_u_L + (gammaL*d_p_L/c_L - c_L*d_rho_L)/(gammaL-1.0)/rho_L;
      TdS = _over_gammaL_1*(d_p_L - d_rho_L*c_L*c_L)/rho_L;
      d_Psi = d_u_L + _over_gammaL_1*(gammaL*d_p_L/c_L - c_L*d_rho_L)/rho_L;

      out->DER.u = ( (1.0+muL)*pow(c_frac, 0.5/muL) + muL*pow(c_frac, (1.0+muL)/muL) )/(1.0+2.0*muL);
      out->DER.u*= TdS;
      out->DER.u-= c_L*pow(c_frac, 0.5/muL)*d_Psi;
      out->DER.u-= c_L*t_v_L/(2.0*muL-1.0)*(c_frac*muL+pow(c_frac, 0.5/muL)*(muL-1.0));//tangential effect
      out->DER.u+= (lambda_y-out->VAR.v)*t_u_L;//tangential effect
      out->DER.u+= out->VAR.u*t_v_L;//tangential effect

      out->DER.p = out->VAR.rho*(out->VAR.u-lambda_x)*out->DER.u;
      out->DER.p+= (lambda_y-out->VAR.v)*t_p_L;//tangential effect

      out->DER.rho = out->VAR.rho*(out->VAR.u-lambda_x)*pow(c_frac, (1.0+muL)/muL)*(d_p_L - d_rho_L*c_L*c_L)/rho_L;
      out->DER.rho = (out->DER.rho + out->DER.p) /c_square;
      out->DER.rho+= (lambda_y-out->VAR.v)*t_rho_L;//tangential effect

      /* geometric effect */
      out->DER.p   = out->DER.p   - geo_factor*out->VAR.rho*out->VAR.u*c_square;
      out->DER.rho = out->DER.rho - geo_factor*out->VAR.rho*out->VAR.u;

      //out->DER.v   = -out->VAR.u*d_v_L*out->VAR.rho/rho_L;
      out->DER.v = (lambda_x-out->VAR.u)*d_v_L*(out->VAR.rho/rho_L);
      out->DER.v+= (lambda_y-out->VAR.v)*t_v_L;
      out->DER.v-= ( ((muL-1.0)*pow(c_frac, 2.0/muL-1.0)-1.0)/(muL-2.0)/out->VAR.rho ) * t_p_L;//tangential effect

      for(it_trans = 0; it_trans < n_trans; ++it_trans)
      {
	out->VAR.trans[it_trans] = wL->VAR.trans[it_trans];
	out->DER.trans[it_trans] = (lambda_x-out->VAR.u)*wL->DER.trans[it_trans]*(out->VAR.rho/rho_L) + (lambda_y-out->VAR.v)*wL->TGT.trans[it_trans];                 
      }                                                                             // ^
                                                                                    // | what's the meaning
                                                                                    // | of this
      //linear calculation:
      //PHIyL = t_u_L + t_p_L/(U[0]out->VAR.rho*c_star);
      //PHIyR = t_u_R - t_p_R/(U[0]out->VAR.rho*c_star);
      //g2 = 0.5* out->VAR.v  *(PHIyL + PHIyR);
      //g4 = 0.5*(out->VAR.v  *(PHIyL - PHIyR) + 2.0*c_star*t_v_L)*U[0]out->VAR.rho*c_star;
      //g3 = out->VAR.v  *t_v_L + 0.5*c_star*(PHIyL - PHIyR);  //error found by LEI corrected
      //g1 = (g4 - out->VAR.v  *(t_p_L - c_square*t_rho_L))/c_square;
      //g5 = out->VAR.v  *t_gm_L;
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

      c_frac = c_star/c_R;
      //TdS = (d_p_R - d_rho_R*c_R*c_R)/(gammaR-1.0)/rho_R;
      //d_Phi = d_u_R - (gammaR*d_p_R/c_R - c_R*d_rho_R)/(gammaR-1.0)/rho_R;
      TdS = (d_p_R - d_rho_R*c_R*c_R)/(gammaR-1.0)/rho_R;
      d_Phi = d_u_R - (gammaR*d_p_R/c_R - c_R*d_rho_R)/(gammaR-1.0)/rho_R;

      out->DER.u = ( (1.0+muR)*pow(c_frac, 0.5/muR) + muR*pow(c_frac, (1.0+muR)/muR) )/(1.0+2.0*muR);
      out->DER.u*= TdS;
      out->DER.u+= c_R*pow(c_frac, 0.5/muR)*d_Phi;
      out->DER.u+= c_R*t_v_R/(2.0*muR-1.0)*(c_frac*(muR-1.0)+pow(c_frac, 0.5/muR)*muR);//tangential effect
      out->DER.u+= (lambda_y-out->VAR.v)*t_u_R;//tangential effect
      out->DER.u+= out->VAR.u*t_v_R;//tangential effect

      out->DER.p = out->VAR.rho*(out->VAR.u-lambda_x)*out->DER.u;
      out->DER.p+= (lambda_y-out->VAR.v)*t_p_R;//tangential effect

      out->DER.rho = out->VAR.rho*(out->VAR.u-lambda_x)*pow(c_frac, (1.0+muR)/muR)*(d_p_R - d_rho_R*c_R*c_R)/rho_R;
      out->DER.rho = (out->DER.rho + out->DER.p) /c_square;
      out->DER.rho += (lambda_y-out->VAR.v)*t_rho_R;//tangential effect

      /* geometric effect */
      out->DER.p   = out->DER.p   - geo_factor*out->VAR.rho*out->VAR.u*c_square;
      out->DER.rho = out->DER.rho - geo_factor*out->VAR.rho*out->VAR.u;

      //out->DER.v   = -out->VAR.u  *d_v_R*out->VAR.rho/rho_R;
      out->DER.v = (lambda_x-out->VAR.u)*d_v_R*out->VAR.rho/rho_R;
      out->DER.v-= (lambda_y-out->VAR.v)*t_v_R;//tangential effect
      out->DER.v-= ( ((muR-1.0)*pow(c_frac, 2.0/muR-1.0)-1.0)/(muR-2.0)/out->VAR.rho ) * t_p_R;//tangential effect

      //D[4] = -(out->VAR.u   - lambda)*d_gm_R;
      for(it_trans = 0; it_trans < n_trans; ++it_trans)
      {
	out->VAR.trans[it_trans] = wR->VAR.trans[it_trans];
	out->DER.trans[it_trans] = (lambda_x-out->VAR.u)*wR->DER.trans[it_trans]*(out->VAR.rho/rho_R) + (lambda_y-out->VAR.v)*wR->TGT.trans[it_trans];                                              // ^
      }                                                                                    // | what's the meaning
                                                                                           // | of this

      //PHIyL = t_u_L + t_p_L/(U[0]out->VAR.rho*c_star);
      //PHIyR = t_u_R - t_p_R/(U[0]out->VAR.rho*c_star);
      //g2 = 0.5*out->VAR.v  *(PHIyL + PHIyR);
      //g4 = 0.5*(out->VAR.v  *(PHIyL - PHIyR) + 2.0*c_star*t_v_R)*U[0]out->VAR.rho*c_star;
      //g3 = out->VAR.v  *t_v_R + 0.5*c_star*(PHIyL - PHIyR);  //error found by LEI corrected
      //g1 = (g4 - out->VAR.v  *(t_p_R - c_square*t_rho_R))/c_square;
      //g5 = out->VAR.v  *t_gm_R;
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

	c_star = c_star_R;
	c_square = c_star*c_star;
	//V = out->VAR.u*out->VAR.u + out->VAR.v*out->VAR.v;
	T_v = t_v_R;
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

	c_star = c_star_L;
	c_square = c_star*c_star;
	//V = out->VAR.u*out->VAR.u + out->VAR.v*out->VAR.v;
	T_v = t_v_L;
      }
      //determine a_L, b_L and d_L
      if(CRW[0]) //the 1-wave is a CRW
      {
	a_L = 1.0;
        b_L = 1.0 / rho_star_L / c_star_L;
	c_frac = c_star_L/c_L;
	TdS = (d_p_L - d_rho_L*c_L*c_L)/(gammaL-1.0)/rho_L;
	d_Psi = d_u_L + (gammaL*d_p_L/c_L - c_L*d_rho_L)/(gammaL-1.0)/rho_L;
	d_L = ( (1.0+muL)*pow(c_frac, 0.5/muL) + muL*pow(c_frac, (1.0+muL)/muL) )/(1.0+2.0*muL);
	d_L = d_L * TdS;
	d_L = d_L - c_L*pow(c_frac, 0.5/muL) * d_Psi;
      }
      else //the 1-wave is a shock
      {
	shk_u_s = -sqrt(0.5*((gammaL+1.0)*p_L   +(gammaL-1.0)*p_star)/rho_star_L);//1-shock speed w.r.t. u_star
	shk_u_L = -sqrt(0.5*((gammaL+1.0)*p_star+(gammaL-1.0)*p_L   )/rho_L);     //1-shock speed w.r.t. u_L

	VAR = sqrt((1-muL)/(rho_L*(p_star+muL*p_L)));
	H1 =  0.5* VAR * (p_star+(1.0+2.0*muL)*p_L)/(p_star+muL*p_L);
	H2 = -0.5*VAR * ((2.0+muL)*p_star + muL*p_L)/(p_star+muL*p_L);
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
	d_R = ( (1.0+muR)*pow(c_frac, 0.5/muR) + muR*pow(c_frac, (1.0+muR)/muR) )/(1.0+2.0*muR);
	d_R = d_R * TdS;
	d_R = d_R + c_R*pow(c_frac, 0.5/muR)*d_Phi;
      }
      else //the 3-wave is a shock
      {
	shk_u_s = sqrt(0.5*((gammaR+1.0)*p_R   + (gammaR-1.0)*p_star)/rho_star_R);//1-shock speed w.r.t. u_star
	shk_u_R = sqrt(0.5*((gammaR+1.0)*p_star+ (gammaR-1.0)*p_R   )/rho_R);     //1-shock speed w.r.t. u_R

	VAR  = sqrt((1.0-muR)/(rho_R*(p_star+muR*p_R)));

	H1 = 0.5* VAR * (p_star+(1+2.0*muR)*p_R)/(p_star+muR*p_R);
	H2 = -0.5*VAR * ((2.0+muR)*p_star+muR*p_R)/(p_star+muR*p_R);
	H3 = -0.5*(p_star-p_R)* VAR /rho_R;

	L_p = -1.0/rho_R + shk_u_R*H2;
	L_u = shk_u_R - rho_R*(c_R*c_R*H2 + H3);
	L_rho = shk_u_R * H3;

	a_R = 1.0 +rho_star_R* shk_u_s * H1;
	b_R = -(shk_u_s/(rho_star_R*c_star_R*c_star_R) + H1);
	d_R = L_rho*d_rho_R + L_u*d_u_R + L_p*d_p_R;
      }


      T_p = 0.5*(t_p_L+t_p_R + (t_u_L-t_u_R)*(out->VAR.rho*c_star));//tangential effect
      T_u = 0.5*(t_u_L+t_u_R + (t_p_L-t_p_R)/(out->VAR.rho*c_star));//tangential effect
      d_L = d_L - a_L*v_L*T_u - b_L*v_L*T_p;
      d_R = d_R - a_R*v_R*T_u - b_R*v_R*T_p;

      detA = a_L*b_R - b_L*a_R;
      u_t_mat = (b_R*d_L - b_L*d_R)/detA;
      p_t_mat = (a_L*d_R - a_R*d_L)/detA;
      D0_p_tau = p_t_mat + out->VAR.v*T_p;
      D0_u_tau = u_t_mat + out->VAR.v*T_u;

      //already total D!
      out->DER.u = u_t_mat + ((out->VAR.u-lambda_x)/out->VAR.rho/c_square) * D0_p_tau;
      out->DER.p = p_t_mat + ((out->VAR.u-lambda_x)*out->VAR.rho)          * D0_u_tau;
      //out->DER.u   = u_t_mat + (u_star-lambda)*p_t_mat/rho_star_R/c_square;
      //out->DER.p   = p_t_mat + rho_star_R*(u_star-lambda) * u_t_mat;
      //out->DER.u   = u_t_mat + (u_star-lambda)*p_t_mat/rho_star_L/c_square;
      //out->DER.p   = p_t_mat + rho_star_L*(u_star-lambda) * u_t_mat;

      if(u_star < lambda_x) //the direction is between the contact discontinuety and the 3-wave
      {
	//out->VAR.rho = rho_star_R;
	//out->VAR.u   =   u_star;
	//out->VAR.v   =   v_R;
	//out->VAR.p   =   p_star;
	//out->gamma   = gammaR;
	//c_star = c_star_R;
	//c_square = c_star*c_star;
	//V = out->VAR.u*out->VAR.u + out->VAR.v*out->VAR.v;
	//D[4] = -(u_star-lambda)*d_gm_R;

	T_rho = t_rho_R + (-t_p_R+T_p)/c_square;//tangential effect
	if(CRW[1]) //the 3-wave is a CRW
	{
	  //already total D!
	  out->DER.rho = out->VAR.rho*(out->VAR.u-lambda_x)*pow(c_star_R/c_R, (1.0+muR)/muR)*(d_p_R - d_rho_R*c_R*c_R)/rho_R;
	  out->DER.rho = (out->DER.rho + out->DER.p + out->VAR.v*T_p) / c_square - (out->VAR.v-lambda_y)*T_rho;

	  out->DER.v = -out->VAR.u*d_v_R*(out->VAR.rho/rho_R);
	  out->DER.v+= lambda_x*d_v_R;
	  out->DER.v+= (out->VAR.u/c_star*((muR-1.0)*pow(c_frac, 2.0/muR-1.0)-1.0)/(muR-2.0)/out->VAR.rho) * t_p_R;//nonlinear tangential effect?
	  out->DER.v-= ((out->VAR.u/c_star+1.0)/out->VAR.rho) * T_p;//nonlinear tangential effect?

	  for(it_trans = 0; it_trans < n_trans; ++it_trans)
	    out->DER.trans[it_trans] = (lambda_x-out->VAR.u*(out->VAR.rho/rho_R))*wR->DER.trans[it_trans] + (lambda_y-out->VAR.v)*wR->TGT.trans[it_trans];
	}                                                       // ^ why moves in
	else //the 3-wave is a shock
	{
	  shk_u_s = sqrt(0.5*((gammaR+1.0)*p_R   + (gammaR-1.0)*p_star)/rho_star_R);
	  shk_u_R = sqrt(0.5*((gammaR+1.0)*p_star+ (gammaR-1.0)*p_R   )/rho_R);

	  VAR = p_R + muR*p_star;
	  H1 = rho_R * p_R    * (1.0 - muR_2) / (VAR*VAR);
	  H2 = rho_R * p_star * (muR_2 - 1.0) / (VAR*VAR);
	  H3 = (p_star + muR*p_R)/VAR;

	  //old code:
	  //out->DER.rho = ((u_star+shk_u_s)/c_star_R/c_star_R - u_star*H1)*p_t_mat + rho_star_R*u_star*shk_u_s*H1*u_t_mat;
	  //out->DER.rho = (out->DER.rho - u_star*(L_p+L_rho+L_u)) / shk_u_s;

	  //L_rho = shk_u_R * H3 * d_rho_R;
	  //L_u = -rho_R * (H2*c_R*c_R + H3) * d_u_R;
	  //L_p = H2 * shk_u_R * d_p_R;
	  //L_v = -t_v_R*(rho_R*c_R*c_R*H2+rho_R*H3);
	  f = shk_u_R*(H2*d_p_R + H3*d_rho_R) - rho_R*(H2*c_R*c_R+H3)*(d_u_R+t_v_R);/*t_v_R is the tangential effect
										      f = L_rho+L_p+L_u+L_v;*/
	  g_p = ( ((shk_u_R+u_R)-lambda_x)/c_square - (out->VAR.u-lambda_x)*H1 ) * D0_p_tau;//replace p_t_mat by D0_p_tau to incorporate the tangential effect
	  g_u = ( (out->VAR.u-lambda_x)*out->VAR.rho*shk_u_s*H1 ) * D0_u_tau;//replace u_t_mat by D0_u_tau to incorporate the tangential effect
	  out->DER.rho = ((out->VAR.u-lambda_x)*f - g_p - g_u) / (-shk_u_s);
	  out->DER.rho+= (lambda_y-out->VAR.v)*T_rho;//tangential effect

	  //old code:
	  //D[2]out->DER.v   = out->VAR.u   * shk_u_R * d_v_R / shk_u_s;
	  //D[2]out->DER.v   = -D[2]out->DER.v   + lambda*d_v_R;
	  out->DER.v = -(out->VAR.u*(shk_u_R*d_v_R-t_p_R/rho_R) + (shk_u_s+out->VAR.u)*T_p/out->VAR.rho) / shk_u_s;
	  out->DER.v+= lambda_x*d_v_R;
	  out->DER.v+= (lambda_y-out->VAR.v)*t_v_R;

	  for(it_trans = 0; it_trans < n_trans; ++it_trans)
	    out->DER.trans[it_trans] = (lambda_x-out->VAR.u*(shk_u_R/shk_u_s))*wR->DER.trans[it_trans] + (lambda_y-out->VAR.v)*wR->TGT.trans[it_trans];
	}

	/*
	//(Uy)* = (Uy)_R - ([(Uy)_R - (Uy)_L].L4)R4
	t_U_D1 = t_U_L1 - t_U_R1;
	t_U_D2 = t_U_L2 - t_U_R2;
	t_U_D3 = t_U_L3 - t_U_R3;
	t_U_D4 = t_U_L4 - t_U_R4;
	H_star = 0.5*V + c_square/(out->gamma  -1.0);
	PROD = t_U_D1*0.5*V - t_U_D2*out->VAR.u   - t_U_D3*out->VAR.v   + t_U_D4;
	PROD = PROD*(out->gamma  -1.0) - c_star*(out->VAR.u  *t_U_D1 - t_U_D2);
	PROD = 0.5*PROD/c_square;
	t_U_s1 = t_U_R1 + PROD;
	t_U_s2 = t_U_R2 + PROD*(out->VAR.u  +c_star);
	t_U_s3 = t_U_R3 + PROD*out->VAR.v  ;
	t_U_s4 = t_U_R4 + PROD*(H_star+out->VAR.u  *c_star);
	t_gm_star = t_gm_R;
	//(dU/dV).B.(Uy)*
	g1 = t_U_s3;
	g2 = out->VAR.v  *(-out->VAR.u  *t_U_s1+t_U_s2)/U[0]out->VAR.rho;
	temp1 = 0.5*V*t_U_s1 - out->VAR.u  *t_U_s2 - out->VAR.v  *t_U_s3 + t_U_s4;
	temp2 = -out->VAR.v  *t_U_s1+t_U_s3;
	g3 = ((out->gamma  -1.0)*temp1 + out->VAR.v  *temp2) / U[0]out->VAR.rho;
	g4 = (out->gamma  -1.0)*out->VAR.v  *temp1 + c_square*temp2;
	g5 = out->VAR.v   * t_gm_star;
	*/
      }
      else //the direction is between the 1-wave and the contact discontinuety
      {
	//out->VAR.rho = rho_star_L;
	//out->VAR.u   =   u_star;
	//out->VAR.v   =   v_L;
	//out->VAR.p   =   p_star;
	//out->gamma   = gammaL;
	//c_star = c_star_L;
	//c_square = c_star*c_star;
	//V = out->VAR.u*out->VAR.u + out->VAR.v*out->VAR.v;
	//D[4] = -(u_star-lambda)*d_gm_L;

	T_rho = t_rho_L + (-t_p_L+T_p)/c_square;//tangential effect

	if(CRW[0]) //the 1-wave is a CRW
	{
	  //already total D!
	  out->DER.rho = out->VAR.rho*(out->VAR.u-lambda_x)*pow(c_star_L/c_L, (1.0+muL)/muL)*(d_p_L - d_rho_L*c_L*c_L)/rho_L;
	  out->DER.rho = (out->DER.rho + out->DER.p + out->VAR.v*T_p) / c_square - (out->VAR.v-lambda_y)*T_rho;

	  out->DER.v = -out->VAR.u*d_v_L*(out->VAR.rho/rho_L);
	  out->DER.v+= lambda_x*d_v_L;
	  out->DER.v-= (out->VAR.u/c_star*((muL-1.0)*pow(c_frac, 2.0/muL-1.0)-1.0)/(muL-2.0)/out->VAR.rho) * t_p_L;//nonlinear tangential effect?
	  out->DER.v+= ((out->VAR.u/c_star-1.0)/out->VAR.rho) * T_p;//nonlinear tangential effect?

	  for(it_trans = 0; it_trans < n_trans; ++it_trans)
	    out->DER.trans[it_trans] = (lambda_x-out->VAR.u*(out->VAR.rho/rho_L))*wL->DER.trans[it_trans] + (lambda_y-out->VAR.v)*wL->TGT.trans[it_trans];
	}                                                       // ^ why moves in
	else //the 1-wave is a shock
	{
	  shk_u_s = -sqrt(0.5*((gammaL+1.0)*p_L   +(gammaL-1.0)*p_star)/rho_star_L);
	  shk_u_L = -sqrt(0.5*((gammaL+1.0)*p_star+(gammaL-1.0)*p_L   )/rho_L);

	  VAR = p_L + muL*p_star;
	  H1 = rho_L * p_L    * (1.0 - muL_2) / (VAR*VAR);
	  H2 = rho_L * p_star * (muL_2 - 1.0) / (VAR*VAR);
	  H3 = (p_star + muL*p_L)/VAR;

	  //old code:
	  //out->DER.rho = ((u_star+shk_u_s)/c_star_L/c_star_L - H1*u_star)*p_t_mat + rho_star_L*u_star*shk_u_s*H1*u_t_mat;
	  //out->DER.rho = (out->DER.rho - u_star*(L_p+L_rho+L_u))/ shk_u_s;

	  //L_rho = shk_u_L * H3 * d_rho_L;
	  //L_u = -rho_L*(H2*c_L*c_L + H3) * d_u_L;
	  //L_p = H2 * shk_u_L * d_p_L;
	  //L_v = -t_v_L*(rho_L*c_L*c_L*H2+rho_L*H3);
	  f = shk_u_L*(H2*d_p_L + H3*d_rho_L) - rho_L*(H2*c_L*c_L+H3)*(d_u_L+t_v_L);/*t_v_R is the tangential effect
										      f = L_rho+L_p+L_u+L_v;*/
	  g_p = ( ((shk_u_L+u_L)-lambda_x)/c_square - (out->VAR.u-lambda_x)*H1 ) * D0_p_tau;//replace p_t_mat by D0_p_tau to incorporate the tangential effect
	  g_u = ( (out->VAR.u-lambda_x)*out->VAR.rho*shk_u_s*H1 ) * D0_u_tau;//replace u_t_mat by D0_u_tau to incorporate the tangential effect
	  out->DER.rho = ((out->VAR.u-lambda_x)*f - g_p - g_u) / (-shk_u_s);
	  out->DER.rho+= (lambda_y-out->VAR.v)*T_rho;//tangential effect

	  //old code:
	  //D[2]out->DER.v   = out->VAR.u   * shk_u_L * d_v_L / shk_u_s;
	  //D[2]out->DER.v   = -D[2]out->DER.v   + lambda*d_v_L;
	  out->DER.v = -(out->VAR.u*(shk_u_L*d_v_L-t_p_L/rho_L) + (shk_u_s+out->VAR.u)*T_p/out->VAR.rho) / shk_u_s;
	  out->DER.v+= lambda_x*d_v_L;
	  out->DER.v+= (lambda_y-out->VAR.v)*t_v_L;

	  for(it_trans = 0; it_trans < n_trans; ++it_trans)
	    out->DER.trans[it_trans] = (lambda_x-out->VAR.u*(shk_u_L/shk_u_s))*wL->DER.trans[it_trans] + (lambda_y-out->VAR.v)*wL->TGT.trans[it_trans];
	}

	/*
	//(Uy)* = (Uy)_L + ([(Uy)_R - (Uy)_L].L1)R1
	t_U_D1 = t_U_R1 - t_U_L1;
	t_U_D2 = t_U_R2 - t_U_L2;
	t_U_D3 = t_U_R3 - t_U_L3;
	t_U_D4 = t_U_R4 - t_U_L4;
	H_star = 0.5*V + c_square/(out->gamma  -1.0);
	PROD = t_U_D1*0.5*V - t_U_D2*out->VAR.u   - t_U_D3*out->VAR.v   + t_U_D4;
	PROD = PROD*(out->gamma  -1.0) + c_star*(out->VAR.u  *t_U_D1 - t_U_D2);
	PROD = 0.5*PROD/c_square;
	t_U_s1 = t_U_L1 + PROD;
	t_U_s2 = t_U_L2 + PROD*(out->VAR.u  -c_star);
	t_U_s3 = t_U_L3 + PROD*out->VAR.v  ;
	t_U_s4 = t_U_L4 + PROD*(H_star-out->VAR.u  *c_star);
	t_gm_star = t_gm_L;
	//(dU/dV).B.(Uy)*
	g1 = t_U_s3;
	g2 = out->VAR.v  *(-out->VAR.u  *t_U_s1+t_U_s2)/U[0]out->VAR.rho;
	temp1 = 0.5*V*t_U_s1 - out->VAR.u  *t_U_s2 - out->VAR.v  *t_U_s3 + t_U_s4;
	temp2 = -out->VAR.v  *t_U_s1+t_U_s3;
	g3 = ((out->gamma  -1.0)*temp1 + out->VAR.v  *temp2) / U[0]out->VAR.rho;
	g4 = (out->gamma  -1.0)*out->VAR.v  *temp1 + c_square*temp2;
	g5 = out->VAR.v   * t_gm_star;
	*/
      }
    //--end of non-sonic case--
    }
    out->DER.u += lambda_y*T_u;
    out->DER.p += lambda_y*T_p;

    /* geometric effect */
    out->DER.p   = out->DER.p   - geo_factor*out->VAR.rho*out->VAR.u*c_square;
    out->DER.rho = out->DER.rho - geo_factor*out->VAR.rho*out->VAR.u;

    //out->DER.rho = out->DER.rho - g1;
    //D[1]out->DER.u   = D[1]out->DER.u   - g2;
    //D[2]out->DER.v   = D[2]out->DER.v   - g3;
    //out->DER.p   = out->DER.p   - g4;
    //D[4] = D[4] - g5;
  //----end of non-trivial case----
  }
}
