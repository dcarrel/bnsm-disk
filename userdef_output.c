#include "pluto.h"
extern double UNIT_G;
/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{
  int i, j, k;
  double ***bern;
  double ***prs, ***rho, ***vx1, ***vx2, ***vx3, *x1;
  double gpot, vsq, rhoe;
    
  prs = d->Vc[PRS];
  rho = d->Vc[RHO];
  vx1 = d->Vc[VX1];
  vx2 = d->Vc[VX2];
  vx3 = d->Vc[VX3];

  x1  = grid->xgc[IDIR]; 
    
  bern = GetUserVar("bern");
  
  DOM_LOOP(k,j,i){
    vsq  = vx1[k][j][i]*vx1[k][j][i] + vx2[k][j][i]*vx2[k][j][i] + vx3[k][j][i]*vx3[k][j][i];
    gpot = -UNIT_G*g_inputParam[MBH]/x1[i];
    rhoe = prs[k][j][i]/rho[k][j][i]/(g_gamma-1);
    
    bern[k][j][i] = 0.5*vsq + gpot + rhoe + prs[k][j][i]/rho[k][j][i]; 
  }
}
/* ************************************************************* */
void ChangeOutputVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
  Image *image;

#if PARTICLES
  //SetOutputVar ("energy",PARTICLES_FLT_OUTPUT, NO);
//  SetOutputVar ("x1",    PARTICLES_FLT_OUTPUT, NO);
  //SetOutputVar ("vx1",   PARTICLES_FLT_OUTPUT, NO);
#endif

}





