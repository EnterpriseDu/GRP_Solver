#include <math.h>

#include "Riemann_solver.h"


/*
 * In this solver we take in primitive variables
 */
void HLL(double wave_speed[2], double F[5], double U[5], double lambda, RSboundary *wL, RSboundary *wR, RSparameters *para)
{
  double SL, SR, SD, UL[5], UR[5], FL[5], FR[5];
  double const gammaL = wL->gamma;
  double const rho_L  = wL->rho;
  double const u_L    = wL->u;
  double const v_L    = wL->v;
  double const p_L    = wL->p;
  double const gammaR = wR->gamma;
  double const rho_R  = wR->rho;
  double const u_R    = wR->u;
  double const v_R    = wR->v;
  double const p_R    = wR->p;
  double const gammaL_1 = gammaL-1.0;
  double const gammaR_1 = gammaR-1.0;
  double const _over_gammaL_1 = 1.0/gammaL_1;
  double const _over_gammaR_1 = 1.0/gammaR_1;

  double const c_L = sqrt(gammaL * p_L / rho_L);
  double const c_R = sqrt(gammaR * p_R / rho_R);

  SL = u_L - c_L;
  if(u_R-c_R < SL)
    SL = u_R-c_R;
  SR = u_R + c_R;
  if(u_L+c_L > SR)
    SR = u_L+c_L;
  SD = 1.0/(SR-SL);
  wave_speed[0] = SL;
  wave_speed[1] = SR;

  UL[0] = rho_L;
  UL[1] = rho_L*u_L;
  UL[2] = rho_L*v_L;
  UL[3] = 0.5*rho_L*(u_L*u_L+v_L*v_L) + p_L*_over_gammaL_1;
  UL[4] = rho_L/(gammaL-1.0);
  UR[0] = rho_R;
  UR[1] = rho_R*u_R;
  UR[2] = rho_R*v_R;
  UR[3] = 0.5*rho_R*(u_R*u_R+v_R*v_R) + p_R*_over_gammaR_1;
  UR[4] = rho_R/(gammaR-1.0);
  FL[0] = UL[1];
  FL[1] = UL[1]*u_L + p_L;
  FL[2] = UL[2]*u_L;
  FL[3] = u_L*(UL[3]+p_L);
  FL[4] = UL[4]*u_L;
  FR[0] = UR[1];
  FR[1] = UR[1]*u_R + p_R;
  FR[2] = UR[2]*u_R;
  FR[3] = u_R*(UR[3]+p_R);
  FR[4] = UR[4]*u_R;

  if(SL > 0.0)
  {
    U[0] = UL[0];
    U[1] = UL[1];
    U[2] = UL[2];
    U[3] = UL[3];
    U[4] = UL[4];
    F[0] = FL[0];
    F[1] = FL[1];
    F[2] = FL[2];
    F[3] = FL[3];
    F[4] = FL[4];
  }
  else if(SR < 0.0)
  {
    U[0] = UR[0];
    U[1] = UR[1];
    U[2] = UR[2];
    U[3] = UR[3];
    U[4] = UR[4];
    F[0] = FR[0];
    F[1] = FR[1];
    F[2] = FR[2];
    F[3] = FR[3];
    F[4] = FR[4];
  }
  else
  {
    U[0] = (SR*UR[0]-SL*UL[0]+FL[0]-FR[0])*SD;
    U[1] = (SR*UR[1]-SL*UL[1]+FL[1]-FR[1])*SD;
    U[2] = (SR*UR[2]-SL*UL[2]+FL[2]-FR[2])*SD;
    U[3] = (SR*UR[3]-SL*UL[3]+FL[3]-FR[3])*SD;
    U[4] = (SR*UR[4]-SL*UL[4]+FL[4]-FR[4])*SD;
    F[0] = (SR*FL[0]-SL*FR[0]+SL*SR*(UR[0]-UL[0]))*SD;
    F[1] = (SR*FL[1]-SL*FR[1]+SL*SR*(UR[1]-UL[1]))*SD;
    F[2] = (SR*FL[2]-SL*FR[2]+SL*SR*(UR[2]-UL[2]))*SD;
    F[3] = (SR*FL[3]-SL*FR[3]+SL*SR*(UR[3]-UL[3]))*SD;
    F[4] = (SR*FL[4]-SL*FR[4]+SL*SR*(UR[4]-UL[4]))*SD;
  }

  /* if(U[1] > lambda) */
  /*   U[4] = gammaL; */
  /* else */
  /*   U[4] = gammaR; */
}


/*
 * In this solver we take in conservative variables
 */
void HLL_consv(double wave_speed[2], double F[5], double U[5], double lambda, RSboundary *wL, RSboundary *wR, RSparameters *para)
{
  double SL, SR, Sstar, SD;
  double aL, aR, aD;
  double UL[5], UR[5], FL[5], FR[5];

  double const gammaL = wL->gamma;
  double const rho_L  = wL->rho;
  double const mom_L  = wL->u;
  double const tan_L  = wL->v;
  double const ene_L  = wL->p;
  double const gammaR = wR->gamma;
  double const rho_R  = wR->rho;
  double const mom_R  = wR->u;
  double const tan_R  = wR->v;
  double const ene_R  = wR->p;
  double const gammaL_1 = gammaL-1.0;
  double const gammaR_1 = gammaR-1.0;
  double const _over_gammaL_1 = 1.0/gammaL_1;
  double const _over_gammaR_1 = 1.0/gammaR_1;

  double const u_L = mom_L/rho_L;
  double const v_L = tan_L/rho_L;
  double const p_L = gammaL_1*(ene_L-0.5*mom_L*u_L);
  double const u_R = mom_R/rho_R;
  double const v_R = tan_R/rho_R;
  double const p_R = gammaR_1*(ene_R-0.5*mom_R*u_R);
  double const c_L = sqrt(gammaL * p_L / rho_L);
  double const c_R = sqrt(gammaR * p_R / rho_R);

  SL = u_L - c_L;
  if(u_R-c_R < SL)
    SL = u_R-c_R;
  SR = u_R + c_R;
  if(u_L+c_L > SR)
    SR = u_L+c_L;
  SD = 1.0/(SR-SL);
  wave_speed[0] = SL;
  wave_speed[1] = SR;

  UL[0] = rho_L;
  UL[1] = mom_L;
  UL[2] = tan_L;
  UL[3] = ene_L;
  UL[4] = rho_L/(gammaL-1.0);
  UR[0] = rho_R;
  UR[1] = mom_R;
  UR[2] = tan_R;
  UR[3] = ene_R;
  UR[4] = rho_R/(gammaR-1.0);
  FL[0] = UL[1];
  FL[1] = UL[1]*u_L + p_L;
  FL[2] = UL[2]*u_L;
  FL[3] = (UL[3]+p_L)*u_L;
  FL[4] = UL[4]*u_L;
  FR[0] = UR[1];
  FR[1] = UR[1]*u_R + p_R;
  FR[2] = UR[2]*u_R;
  FR[3] = (UR[3]+p_R)*u_R;
  FR[4] = UR[4]*u_R;

  if(SL > 0.0)
  {
    U[0] = UL[0];
    U[1] = UL[1];
    U[2] = UL[2];
    U[3] = UL[3];
    U[4] = UL[4];
    F[0] = FL[0];
    F[1] = FL[1];
    F[2] = FL[2];
    F[3] = FL[3];
    F[4] = FL[4];
  }
  else if(SR < 0.0)
  {
    U[0] = UR[0];
    U[1] = UR[1];
    U[2] = UR[2];
    U[3] = UR[3];
    U[4] = UR[4];
    F[0] = FR[0];
    F[1] = FR[1];
    F[2] = FR[2];
    F[3] = FR[3];
    F[4] = FR[4];
  }
  else
  {
    U[0] = (SR*UR[0]-SL*UL[0]+FL[0]-FR[0])*SD;
    U[1] = (SR*UR[1]-SL*UL[1]+FL[1]-FR[1])*SD;
    U[2] = (SR*UR[2]-SL*UL[2]+FL[2]-FR[2])*SD;
    U[3] = (SR*UR[3]-SL*UL[3]+FL[3]-FR[3])*SD;
    U[4] = (SR*UR[4]-SL*UL[4]+FL[4]-FR[4])*SD;
    F[0] = (SR*FL[0]-SL*FR[0]+SL*SR*(UR[0]-UL[0]))*SD;
    F[1] = (SR*FL[1]-SL*FR[1]+SL*SR*(UR[1]-UL[1]))*SD;
    F[2] = (SR*FL[2]-SL*FR[2]+SL*SR*(UR[2]-UL[2]))*SD;
    F[3] = (SR*FL[3]-SL*FR[3]+SL*SR*(UR[3]-UL[3]))*SD;
    F[4] = (SR*FL[4]-SL*FR[4]+SL*SR*(UR[4]-UL[4]))*SD;
  }

  /* if(U[1] > lambda) */
  /*   U[4] = gammaL; */
  /* else */
  /*   U[4] = gammaR; */
}
