/* /////////////////////////////////////////////////////////////////// */
/*! \file  
 *  \brief Specification of explicit first and second viscosity coefficients*/
/* /////////////////////////////////////////////////////////////////// */
#include "pluto.h"
/* ************************************************************************** */
void Visc_nu(double *v, double x1, double x2, double x3,
                        double *nu1, double *nu2)
/*! 
 *
 *  \param [in]      v  pointer to data array containing cell-centered quantities
 *  \param [in]      x1 real, coordinate value 
 *  \param [in]      x2 real, coordinate value 
 *  \param [in]      x3 real, coordinate value 
 *  \param [in, out] nu1  pointer to first viscous coefficient
 *  \param [in, out] nu2  pointer to second viscous coefficient
 *
 *  \return This function has no return value.
 * ************************************************************************** */
{
//  printf("Visocisty called\n");
  if(g_time > g_inputParam[TSTART]){//v[RHO] > 1.e-4){
    *nu1 = g_inputParam[ALPHA] * v[PRS] * pow(x1, 1.5)/v[RHO]/sqrt(g_inputParam[MBH]*CONST_G);
  }
  else{
    *nu1 = 0.0;
  }
  *nu2 = 0.0;
}
