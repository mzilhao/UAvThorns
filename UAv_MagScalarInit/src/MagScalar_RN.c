#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"


void MagScalar_RN(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int imin[3], imax[3];

  for (int d = 0; d < 3; ++ d) {
    imin[d] = 0;
    imax[d] = cctk_lsh[d];
  }


  for (int i = imin[0]; i < imax[0]; ++i) {
    for (int j = imin[1]; j < imax[1]; ++j) {
      for (int k = imin[2]; k < imax[2]; ++k) {

        const int ind = CCTK_GFINDEX3D (cctkGH, i, j, k);

        const CCTK_REAL x1 = x[ind] - x0;
        const CCTK_REAL y1 = y[ind] - y0;
        const CCTK_REAL z1 = z[ind] - z0;

        const CCTK_REAL R2 = x1*x1 + y1*y1 + z1*z1;
        const CCTK_REAL R  = sqrt(R2);
        const CCTK_REAL R3 = R2 * R;

        const CCTK_REAL psi2 = pow(1 + 0.5 * par_m / R, 2) - 0.25 * pow(par_q / R, 2);
        const CCTK_REAL psi1 = sqrt(psi2);

        gxx[ind] = pow(psi2, 2);
        gxy[ind] = 0;
        gxz[ind] = 0;
        gyy[ind] = pow(psi2, 2);
        gyz[ind] = 0;
        gzz[ind] = pow(psi2, 2);

        kxx[ind] = 0;
        kxy[ind] = 0;
        kxz[ind] = 0;
        kyy[ind] = 0;
        kyz[ind] = 0;
        kzz[ind] = 0;


        /* EMG terms */
        Ex[ind]  = par_q * x1/R3 / pow(psi2, 3);
        Ey[ind]  = par_q * y1/R3 / pow(psi2, 3);
        Ez[ind]  = par_q * z1/R3 / pow(psi2, 3);

        /* account here for the different normalization of the MagScalarEvolve
           thorn when compared with ProcaEvolve */
        Ex[ind] *= 1/sqrt(4*M_PI);
        Ey[ind] *= 1/sqrt(4*M_PI);
        Ez[ind] *= 1/sqrt(4*M_PI);

        Ax[ind] = 0;
        Ay[ind] = 0;
        Az[ind] = 0;

        Aphi[ind]  = 0;

        Zeta[ind]  = 0;


        /* scalar terms */
        phi1[ind]  = 0;
        phi2[ind]  = 0;
        Kphi1[ind] = 0;
        Kphi2[ind] = 0;

        /* lapse */
        if ( CCTK_EQUALS(initial_lapse, "psi^n") ) {
          alp[ind] = pow(psi1, initial_lapse_psi_exponent);
        }

      }
    }
  }

}
