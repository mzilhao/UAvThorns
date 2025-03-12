
#include <stdio.h>
#include <math.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"


void KerrQI(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT i,j,k;

  const CCTK_REAL delta2 = mass*mass - spin*spin;
  const CCTK_REAL delta  = sqrt(delta2);

  CCTK_REAL x1, y1, z1;
  CCTK_REAL RR, RR2, rBL, costh, costh2, rho2, rho, psi4, psi2, psi1, sigma, hh, alpha, alpha2;
  CCTK_REAL RRrBL, sinth2, sinth, aux, pert, R0pert2;
  CCTK_REAL Athph, ARph, HF, HE, Axx, Axy, Axz, Ayy, Ayz;

  CCTK_REAL R_x, R_y, R_z;
  CCTK_REAL sinth2ph_x, sinth2ph_y, sinthth_x, sinthth_y, sinthth_z, sinthx_th, sinthy_th, sinthz_th;


  for (k = 0; k < cctk_lsh[2]; ++k) {
    for (j = 0; j < cctk_lsh[1]; ++j) {
      for (i = 0; i < cctk_lsh[0]; ++i) {

        const CCTK_INT ind = CCTK_GFINDEX3D (cctkGH, i, j, k);

        x1     = x[ind] - x0;
        y1     = y[ind] - y0;
        z1     = z[ind] - z0;

        RR2    = x1*x1 + y1*y1 + z1*z1;
        RR     = sqrt(RR2);

        R0pert2 = (RR - R0pert)*(RR - R0pert);

        costh  = z1/RR;
        costh2 = costh*costh;
        sinth2 = 1. - costh2;
        sinth  = sqrt(sinth2);

        R_x    = x1/RR;
        R_y    = y1/RR;
        R_z    = z1/RR;

        sinth2ph_x = -y1/RR2;
        sinth2ph_y =  x1/RR2;

        sinthth_x  = z1*x1/(RR*RR2); 
        sinthth_y  = z1*y1/(RR*RR2); 
        sinthth_z  = -sinth2/RR; 

        sinthx_th  = x1 * costh;
        sinthy_th  = y1 * costh;
        sinthz_th  = -RR * sinth2;


        rBL    = RR + mass + 0.25*delta2 / RR;   // Boyer-Lindquist coordinate r

        RRrBL  = RR2 + RR*mass + 0.25*delta2;

        rho2   = rBL*rBL + spin*spin * costh2;
        rho    = sqrt(rho2);

        // sigma = (2.*mass*rBL) / rho2;
        sigma  = (2.*mass*RRrBL) * RR / (RRrBL*RRrBL + RR2*spin*spin * costh2);

        hh     = (1 + sigma) / (RRrBL*RRrBL + RR2*spin*spin * costh2);

        psi4   = rho2 / RR2;
        psi2   = sqrt(psi4);
        psi1   = sqrt(psi2);


        // non-axisymmetric perturbation
        pert = 1. + AA * (x1*x1 - y1*y1)/(mass*mass) * exp( -2.*R0pert2/delta2 );

        // 3-metric
        gxx[ind] = psi4/pert * ( 1. + spin*spin * hh * y1*y1 );
        gxy[ind] = - psi4/pert * spin*spin * hh * x1*y1;
        gxz[ind] = 0;
        gyy[ind] = psi4/pert * ( 1. + spin*spin * hh * x1*x1 );
        gyz[ind] = 0;
        gzz[ind] = psi4/pert;


        alpha  = (RR + 0.5*delta)*(RR - 0.5*delta) / RR *
                 1. / sqrt(rBL*rBL + spin*spin * ( 1. + sigma*sinth2));
        alpha2 = alpha*alpha;

        HF     = - spin*spin*spin * alpha * sigma/rho * costh;  // we are dividing by sinth2
        Athph  = HF / RR;                                       // we are dividing by sinth

        aux    =  rho2 * (rBL*rBL - spin*spin) + 2.*rBL*rBL * (rBL*rBL + spin*spin);

        HE     = spin*mass * aux / (rho*rho*rho) * 
                 1. / sqrt(rBL*rBL + spin*spin * ( 1. + sigma*sinth2));

        ARph   = HE / RR2;                                       // we are dividing by sinth2


        Axx = 2.*ARph *  R_x * sinth2ph_x                     +  2.*Athph *  sinthth_x * sinth2ph_x;
        Axy =    ARph * (R_x * sinth2ph_y + R_y * sinth2ph_x) +     Athph * (sinthth_x * sinth2ph_y + sinthth_y * sinth2ph_x);
        Axz =    ARph *                     R_z * sinth2ph_x  +     Athph *                           sinthth_z * sinth2ph_x;
        Ayy = 2.*ARph *  R_y * sinth2ph_y                     +  2.*Athph *  sinthth_y * sinth2ph_y;
        Ayz =    ARph *                     R_z * sinth2ph_y  +     Athph *                           sinthth_z * sinth2ph_y;

        kxx[ind] = Axx / psi2;
        kxy[ind] = Axy / psi2;
        kxz[ind] = Axz / psi2;
        kyy[ind] = Ayy / psi2;
        kyz[ind] = Ayz / psi2;
        kzz[ind] = 0.;

        // lapse
        if ( CCTK_EQUALS(initial_lapse, "psi^n") ) {
          alp[ind] = pow(psi1, initial_lapse_psi_exponent);
        }

      } /* for i */
    }   /* for j */
  }     /* for k */

}
