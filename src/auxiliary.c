
#include "Riemann_solver.h"


void copy_EulerVar(EulerVar *dest, EulerVar *source, int n)
{
  dest->rho = source->rho;
  dest->u   = source->u;
  dest->v   = source->v;
  dest->p   = source->p;

  int k;
  for(k = 0; k < n; ++k)
    dest->trans[k] = source->trans[k];
}
void copy_EulerPack(EulerPack *dest, EulerPack *source, int n_trans)
{
  dest->gamma   = source->gamma;
  dest->d_gamma = source->d_gamma;

  copy_EulerVar(&(dest->VAR), &(source->VAR), n_trans);
  copy_EulerVar(&(dest->DER), &(source->DER), n_trans);
  copy_EulerVar(&(dest->TGT), &(source->TGT), n_trans);
}
