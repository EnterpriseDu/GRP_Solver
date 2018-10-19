#include <math.h>

#include "Riemann_solver.h"


void vacuum_Godunov(double wave_speed[2], EulerPack *out, double lambda_x, double lambda_y, int n_trans, EulerPack const *wL, EulerPack const *wR, RSparameters *para)
{
  double const eps = para->eps;

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


  double TdS, d_Psi, d_Phi, RIE_SL, RIE_RR, PHIyL, PHIyR;
  double c_star_L, c_star_R, C, c_frac, c_star, c_square;
  double g1, g2, g3, g4, g5;
  int it_trans;


  if(((rho_R < eps) || (p_R < eps)) && ((rho_L < eps) || (p_L < eps)))
  {
    wave_speed[0] = 0.0;
    wave_speed[1] = 0.0;
    out->gamma = 0.0;
    out->VAR.rho = 0.0;
    out->VAR.u   = 0.0;
    out->VAR.v   = 0.0;
    out->VAR.p   = 0.0;
    for(it_trans = 0; it_trans < n_trans; ++it_trans)
      out->VAR.trans[it_trans] = 0.0;
  }
  else if((rho_R < eps) || (p_R < eps))//right dry,the left wave is a Rarefaction
  {
    wave_speed[0] = u_L - c_L;
    wave_speed[1] = 0.0;
    RIE_SL = u_L + 2.0*_over_gammaL_1*c_L;//Riemann Invariant

    if(lambda_x < wave_speed[0])//x lie to the left of fan-region of
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
    else if(lambda_x < RIE_SL)//x lie in the  fan-region: SONIC CASE
    {
      out->gamma   = gammaL;
      out->VAR.u   = muL*(u_L+2.0*_over_gammaL_1*(c_L+lambda_x));
      c_star = out->VAR.u - lambda_x;
      c_square = c_star*c_star;
      out->VAR.p   = pow(c_star/c_L, 2.0*gammaL*_over_gammaL_1) * p_L;
      out->VAR.rho = gammaL*out->VAR.p / c_square;
      out->VAR.v   = v_L;
      for(it_trans = 0; it_trans < n_trans; ++it_trans)
	out->VAR.trans[it_trans] = wL->VAR.trans[it_trans];
    }
    else//x lie at the right of fan-region 
    {
      out->gamma = 0.0;
      out->VAR.rho = 0.0;
      out->VAR.u   = 0.0;
      out->VAR.v   = 0.0;
      out->VAR.p   = 0.0;
      for(it_trans = 0; it_trans < n_trans; ++it_trans)
	out->VAR.trans[it_trans] = 0.0;
    }
  }
  else if((rho_L < eps) || (p_L < eps))//left dry,the right wave is a Rarefaction
  {
    wave_speed[0] = 0.0;
    wave_speed[1] = u_R + c_R;
    RIE_RR = u_R - 2.0*_over_gammaR_1*c_R;//Riemann Invariant

    if(wave_speed[1] < lambda_x)
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
    else if(RIE_RR < lambda_x)//sonic case
    {
      out->gamma   = gammaR;
      out->VAR.u   = muR*(u_R-2.0*_over_gammaR_1*(c_R-lambda_x));
      c_star = lambda_x - out->VAR.u;
      c_square = c_star*c_star;
      out->VAR.p   = pow(c_star/c_R, 2.0*gammaR*_over_gammaR_1) * p_R;
      out->VAR.rho = gammaR*out->VAR.p  /c_square;
      out->VAR.v   = v_R;
      for(it_trans = 0; it_trans < n_trans; ++it_trans)
	out->VAR.trans[it_trans] = wR->VAR.trans[it_trans];
    }
    else
    {
      out->gamma = 0.0;
      out->VAR.rho = 0.0;
      out->VAR.u   = 0.0;
      out->VAR.v   = 0.0;
      out->VAR.p   = 0.0;
      for(it_trans = 0; it_trans < n_trans; ++it_trans)
	out->VAR.trans[it_trans] = 0.0;
    }
  }
  else//middle vacuum
  {
    wave_speed[0] = u_L - c_L;
    wave_speed[1] = u_R + c_R;
    RIE_SL = u_L + 2.0*_over_gammaL_1*c_L;
    RIE_RR = u_R - 2.0*_over_gammaR_1*c_R;

    if(lambda_x < wave_speed[0])
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
    else if(lambda_x < RIE_SL)
    {
      out->gamma   = gammaL;
      out->VAR.u   = muL*(u_L+2.0*_over_gammaL_1*(c_L+lambda_x));
      c_star = out->VAR.u - lambda_x;
      c_square = c_star*c_star;
      out->VAR.p   = pow(c_star/c_L, 2.0*gammaL*_over_gammaL_1) * p_L;
      out->VAR.rho = gammaL*out->VAR.p / c_square;
      out->VAR.v   = v_L;
      for(it_trans = 0; it_trans < n_trans; ++it_trans)
	out->VAR.trans[it_trans] = wL->VAR.trans[it_trans];
    }
    else if(wave_speed[1] < lambda_x)
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
    else if(RIE_RR < lambda_x)
    {
      out->gamma   = gammaR;
      out->VAR.u   = muR*(u_R-2.0*_over_gammaR_1*(c_R-lambda_x));
      c_star = lambda_x - out->VAR.u;
      c_square = c_star*c_star;
      out->VAR.p   = pow(c_star/c_R, 2.0*gammaR*_over_gammaR_1) * p_R;
      out->VAR.rho = gammaR*out->VAR.p  /c_square;
      out->VAR.v   = v_R;
      for(it_trans = 0; it_trans < n_trans; ++it_trans)
	out->VAR.trans[it_trans] = wR->VAR.trans[it_trans];
    }
    else
    {
      out->gamma = 0.0;
      out->VAR.rho = 0.0;
      out->VAR.u   = 0.0;
      out->VAR.v   = 0.0;
      out->VAR.p   = 0.0;
      for(it_trans = 0; it_trans < n_trans; ++it_trans)
	out->VAR.trans[it_trans] = 0.0;
    }
  }
}


void vacuum(double wave_speed[2], EulerPack *out, double lambda_x, double lambda_y, int n_trans, EulerPack const *wL, EulerPack const *wR, RSparameters *para)
{
  double const eps = para->eps;

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


  double TdS, d_Psi, d_Phi, RIE_SL, RIE_RR, PHIyL, PHIyR;
  double c_star_L, c_star_R, C, c_frac, c_star, c_square;
  double g1, g2, g3, g4, g5;
  int it_trans;


  if(((rho_R < eps) || (p_R < eps)) && ((rho_L < eps) || (p_L < eps)))
  {
    wave_speed[0] = 0.0;
    wave_speed[1] = 0.0;
    out->gamma = 0.0;
    out->VAR.rho = 0.0;
    out->VAR.u   = 0.0;
    out->VAR.v   = 0.0;
    out->VAR.p   = 0.0;
    for(it_trans = 0; it_trans < n_trans; ++it_trans)
      out->VAR.trans[it_trans] = 0.0;
    out->DER.rho = 0.0;
    out->DER.u   = 0.0;
    out->DER.v   = 0.0;
    out->DER.p   = 0.0;
    for(it_trans = 0; it_trans < n_trans; ++it_trans)
      out->DER.trans[it_trans] = 0.0;
  }
  else if((rho_R < eps) || (p_R < eps))//right dry,the left wave is a Rarefaction
  {
    wave_speed[0] = u_L - c_L;
    wave_speed[1] = 0.0;
    RIE_SL = u_L + 2.0*_over_gammaL_1*c_L;//Riemann Invariant

    if(lambda_x < wave_speed[0])//x lie to the left of fan-region of
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

      for(it_trans = 0; it_trans < n_trans; ++it_trans)
      {
	out->VAR.trans[it_trans] = wL->VAR.trans[it_trans];
	out->DER.trans[it_trans] = (lambda_x-wL->VAR.u) * wL->DER.trans[it_trans] + (lambda_y-wL->VAR.v) * wL->TGT.trans[it_trans];
      }
    }
    else if(lambda_x < RIE_SL)//x lie in the  fan-region: SONIC CASE
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

      //out->DER.v   = -out->VAR.u*d_v_L*out->VAR.rho/rho_L;
      out->DER.v = (lambda_x-out->VAR.u)*d_v_L*(out->VAR.rho/rho_L);
      out->DER.v+= (lambda_y-out->VAR.v)*t_v_L;
      out->DER.v-= ( ((muL-1.0)*pow(c_frac, 2.0/muL-1.0)-1.0)/(muL-2.0)/out->VAR.rho ) * t_p_L;//tangential effect

      for(it_trans = 0; it_trans < n_trans; ++it_trans)
      {
	out->VAR.trans[it_trans] = wL->VAR.trans[it_trans];
	out->DER.trans[it_trans] = (lambda_x-out->VAR.u)*wL->DER.trans[it_trans]*(out->VAR.rho/rho_L) + (lambda_y-out->VAR.v)*wL->TGT.trans[it_trans];                                              // ^
      }
      /*
      out->VAR.u   = muL*(u_L+2.0*(c_L+lambda)*_over_gammaL_1);
      c_star       = out->VAR.u - lambda;
      c_square     = c_star*c_star;
      out->VAR.p   = pow(c_star/c_L, 2.0*gammaL*_over_gammaL_1) * p_L;
      out->VAR.rho = gammaL * out->VAR.p / c_square;
      out->VAR.v   = v_L;
      out->gamma   = gammaL;

      c_frac = c_star/c_L;
      TdS = (_over_gammaL_1*(d_p_L - d_rho_L*c_L*c_L))/rho_L;
      d_Psi = d_u_L + (_over_gammaL_1*(gammaL*d_p_L/c_L - c_L*d_rho_L))/rho_L;

      D[1] = ( (1.0+muL)*pow(c_frac, 0.5/muL) + muL*pow(c_frac, (1.0+muL)/muL) )/(1.0+2.0*muL);
      D[1] = D[1] * TdS;
      D[1] = D[1] - c_L*pow(c_frac, 0.5/muL)*d_Psi;
      D[3] = out->VAR.rho*(out->VAR.u-lambda)*D[1];

      D[0] = out->VAR.rho*(out->VAR.u-lambda)*pow(c_frac, (1.0+muL)/muL)*(d_p_L - d_rho_L*c_L*c_L)/rho_L;
      D[0] = (D[0] + D[3]) /c_square;
      D[2] = -out->VAR.u*d_v_L*out->VAR.rho/rho_L;

      PHIyL = t_u_L + t_p_L/(out->VAR.rho*c_star);
      PHIyR = t_u_R - t_p_R/(out->VAR.rho*c_star);
      g2 = 0.5*out->VAR.v*(PHIyL + PHIyR);
      g4 = 0.5*(out->VAR.v*(PHIyL - PHIyR) + 2.0*c_star*t_v_L)*out->VAR.rho*c_star;
      g3 = out->VAR.v*t_v_L - 0.5*c_star*(PHIyL - PHIyR);
      g1 = (g4 - out->VAR.v*(t_p_L - c_square*t_rho_L))/c_square;
      g5 = out->VAR.v*t_gm_L;
      D[0] = D[0] - g1;
      D[1] = D[1] - g2;
      D[2] = D[2] - g3;
      D[3] = D[3] - g4;
      D[4] = D[4] - g5;
      */
    }
    else//x lie at the right of fan-region 
    {
      out->gamma = 0.0;
      out->VAR.rho = 0.0;
      out->VAR.u   = 0.0;
      out->VAR.v   = 0.0;
      out->VAR.p   = 0.0;
      for(it_trans = 0; it_trans < n_trans; ++it_trans)
	out->VAR.trans[it_trans] = 0.0;
      out->DER.rho = 0.0;
      out->DER.u   = 0.0;
      out->DER.v   = 0.0;
      out->DER.p   = 0.0;
      for(it_trans = 0; it_trans < n_trans; ++it_trans)
	out->DER.trans[it_trans] = 0.0;
    }
  }
  else if((rho_L < eps) || (p_L < eps))//left dry,the right wave is a Rarefaction
  {
    wave_speed[0] = 0.0;
    wave_speed[1] = u_R + c_R;
    RIE_RR = u_R - 2.0*_over_gammaR_1*c_R;//Riemann Invariant

    if(wave_speed[1] < lambda_x)
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

      for(it_trans = 0; it_trans < n_trans; ++it_trans)
      {
	out->VAR.trans[it_trans] = wR->VAR.trans[it_trans];
	out->DER.trans[it_trans] = (lambda_x-wR->VAR.u) * wR->DER.trans[it_trans] + (lambda_y-wR->VAR.v) * wR->TGT.trans[it_trans];
      }
      /*
      c_star = c_R;
      c_square = c_star*c_star;
      out->VAR.rho = rho_R;
      out->VAR.u =   u_R;
      out->VAR.v =   v_R;
      out->VAR.p =   p_R;
      out->gamma = gammaR;

      D[0] = -out->VAR.u*d_rho_R - out->VAR.v*t_rho_R - out->VAR.rho*(d_u_R+t_v_R);
      D[1] = -out->VAR.u*  d_u_R - out->VAR.v*  t_u_R - d_p_R/out->VAR.rho;
      D[2] = -out->VAR.u*  d_v_R - out->VAR.v*  t_v_R - t_p_R/out->VAR.rho;
      D[3] = -out->VAR.u*  d_p_R - out->VAR.v*  t_p_R - out->VAR.rho*(d_u_R+t_v_R)*c_square;
      D[4] = -out->VAR.u* d_gm_R - out->VAR.v*  t_gm_R;

      D[0] = D[0] + lambda * d_rho_R;
      D[1] = D[1] + lambda * d_u_R;
      D[2] = D[2] + lambda * d_v_R;
      D[3] = D[3] + lambda * d_p_R;
      D[4] = D[4] + lambda * d_gm_R;
      */
    }
    else if(RIE_RR < lambda_x)//sonic case
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

      //out->DER.v   = -out->VAR.u  *d_v_R*out->VAR.rho/rho_R;
      out->DER.v = (lambda_x-out->VAR.u)*d_v_R*out->VAR.rho/rho_R;
      out->DER.v-= (lambda_y-out->VAR.v)*t_v_R;//tangential effect
      out->DER.v-= ( ((muR-1.0)*pow(c_frac, 2.0/muR-1.0)-1.0)/(muR-2.0)/out->VAR.rho ) * t_p_R;//tangential effect

      //D[4] = -(out->VAR.u   - lambda)*d_gm_R;
      for(it_trans = 0; it_trans < n_trans; ++it_trans)
      {
	out->VAR.trans[it_trans] = wR->VAR.trans[it_trans];
	out->DER.trans[it_trans] = (lambda_x-out->VAR.u)*wR->DER.trans[it_trans]*(out->VAR.rho/rho_R) + (lambda_y-out->VAR.v)*wR->TGT.trans[it_trans];                                              // ^
      }
      /*
      out->VAR.u = muR*(u_R-2.0*(c_R-lambda)/(gammaR-1.0));
      c_star = lambda-out->VAR.u;
      c_square = c_star*c_star;
      out->VAR.p = pow(c_star/c_R, 2.0*gammaR/(gammaR-1.0)) * p_R;
      out->VAR.rho = gammaR*out->VAR.p/c_square;
      out->VAR.v = v_R;
      out->gamma = gammaR;

      c_frac = c_star/c_R;
      TdS = (d_p_R - d_rho_R*c_R*c_R)/(gammaR-1.0)/rho_R;
      d_Phi = d_u_R - (gammaR*d_p_R/c_R - c_R*d_rho_R)/(gammaR-1.0)/rho_R;

      D[1] = ( (1.0+muR)*pow(c_frac, 0.5/muR) + muR*pow(c_frac, (1.0+muR)/muR))/(1.0+2.0*muR);
      D[1] = D[1] * TdS;
      D[1] = D[1] + c_R*pow(c_frac, 0.5/muR)*d_Phi;
      D[3] = out->VAR.rho*(out->VAR.u-lambda)*D[1];

      D[0] = out->VAR.rho*(out->VAR.u-lambda)*pow(c_frac, (1.0+muR)/muR)*(d_p_R - d_rho_R*c_R*c_R)/rho_R;
      D[0] = (D[0] + D[3]) /c_square;
      D[2] = -out->VAR.u*d_v_L*out->VAR.rho/rho_R;

      PHIyL = t_u_L + t_p_L/(out->VAR.rho*c_star);
      PHIyR = t_u_R - t_p_R/(out->VAR.rho*c_star);
      g2 = 0.5*out->VAR.v*(PHIyL + PHIyR);
      g4 = 0.5*(out->VAR.v*(PHIyL - PHIyR) + 2.0*c_star*t_v_R)*out->VAR.rho*c_star;
      g3 = out->VAR.v*t_v_R - 0.5*c_star*(PHIyL - PHIyR);
      g1 = (g4 - out->VAR.v*(t_p_R - c_square*t_rho_R))/c_square;
      g5 = out->VAR.v*t_gm_R;
      D[0] = D[0] - g1;
      D[1] = D[1] - g2;
      D[2] = D[2] - g3;
      D[3] = D[3] - g4;
      D[4] = D[4] - g5;
      */
    }
    else
    {
      out->gamma = 0.0;
      out->VAR.rho = 0.0;
      out->VAR.u   = 0.0;
      out->VAR.v   = 0.0;
      out->VAR.p   = 0.0;
      for(it_trans = 0; it_trans < n_trans; ++it_trans)
	out->VAR.trans[it_trans] = 0.0;
      out->DER.rho = 0.0;
      out->DER.u   = 0.0;
      out->DER.v   = 0.0;
      out->DER.p   = 0.0;
      for(it_trans = 0; it_trans < n_trans; ++it_trans)
	out->DER.trans[it_trans] = 0.0;
    }
  }
  else//middle vacuum
  {
    wave_speed[0] = u_L - c_L;
    wave_speed[1] = u_R + c_R;
    RIE_SL = u_L + 2.0*_over_gammaL_1*c_L;
    RIE_RR = u_R - 2.0*_over_gammaR_1*c_R;

    if(lambda_x < wave_speed[0])
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

      for(it_trans = 0; it_trans < n_trans; ++it_trans)
      {
	out->VAR.trans[it_trans] = wL->VAR.trans[it_trans];
	out->DER.trans[it_trans] = (lambda_x-wL->VAR.u) * wL->DER.trans[it_trans] + (lambda_y-wL->VAR.v) * wL->TGT.trans[it_trans];
      }
    }
    else if(lambda_x < RIE_SL)
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

      //out->DER.v   = -out->VAR.u*d_v_L*out->VAR.rho/rho_L;
      out->DER.v = (lambda_x-out->VAR.u)*d_v_L*(out->VAR.rho/rho_L);
      out->DER.v+= (lambda_y-out->VAR.v)*t_v_L;
      out->DER.v-= ( ((muL-1.0)*pow(c_frac, 2.0/muL-1.0)-1.0)/(muL-2.0)/out->VAR.rho ) * t_p_L;//tangential effect

      for(it_trans = 0; it_trans < n_trans; ++it_trans)
      {
	out->VAR.trans[it_trans] = wL->VAR.trans[it_trans];
	out->DER.trans[it_trans] = (lambda_x-out->VAR.u)*wL->DER.trans[it_trans]*(out->VAR.rho/rho_L) + (lambda_y-out->VAR.v)*wL->TGT.trans[it_trans];                                              // ^
      }
    }
    else if(wave_speed[1] < lambda_x)
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

      for(it_trans = 0; it_trans < n_trans; ++it_trans)
      {
	out->VAR.trans[it_trans] = wR->VAR.trans[it_trans];
	out->DER.trans[it_trans] = (lambda_x-wR->VAR.u) * wR->DER.trans[it_trans] + (lambda_y-wR->VAR.v) * wR->TGT.trans[it_trans];
      }
    }
    else if(RIE_RR < lambda_x)
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

      //out->DER.v   = -out->VAR.u  *d_v_R*out->VAR.rho/rho_R;
      out->DER.v = (lambda_x-out->VAR.u)*d_v_R*out->VAR.rho/rho_R;
      out->DER.v-= (lambda_y-out->VAR.v)*t_v_R;//tangential effect
      out->DER.v-= ( ((muR-1.0)*pow(c_frac, 2.0/muR-1.0)-1.0)/(muR-2.0)/out->VAR.rho ) * t_p_R;//tangential effect

      //D[4] = -(out->VAR.u   - lambda)*d_gm_R;
      for(it_trans = 0; it_trans < n_trans; ++it_trans)
      {
	out->VAR.trans[it_trans] = wR->VAR.trans[it_trans];
	out->DER.trans[it_trans] = (lambda_x-out->VAR.u)*wR->DER.trans[it_trans]*(out->VAR.rho/rho_R) + (lambda_y-out->VAR.v)*wR->TGT.trans[it_trans];                                              // ^
      }
    }
    else
    {
      out->gamma = 0.0;
      out->VAR.rho = 0.0;
      out->VAR.u   = 0.0;
      out->VAR.v   = 0.0;
      out->VAR.p   = 0.0;
      for(it_trans = 0; it_trans < n_trans; ++it_trans)
	out->VAR.trans[it_trans] = 0.0;
      out->DER.rho = 0.0;
      out->DER.u   = 0.0;
      out->DER.v   = 0.0;
      out->DER.p   = 0.0;
      for(it_trans = 0; it_trans < n_trans; ++it_trans)
	out->DER.trans[it_trans] = 0.0;
    }
  }
}
