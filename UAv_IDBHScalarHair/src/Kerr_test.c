
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "util_Table.h"

#define SMALL (1.e-9)


void UAv_Kerr_test(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dxsq = CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(0);
  const CCTK_REAL dysq = CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(1);

  /* let's create arrays for the F1, F2, F0, phi0, W functions */
  const CCTK_INT N_points = cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; // total points

  CCTK_REAL *F1, *F2, *F0, *phi0, *W;
  CCTK_REAL *dW_dR, *dW_dth, *d2W_dth2;

  F1   = (CCTK_REAL *) malloc(N_points * sizeof(CCTK_REAL));
  F2   = (CCTK_REAL *) malloc(N_points * sizeof(CCTK_REAL));
  F0   = (CCTK_REAL *) malloc(N_points * sizeof(CCTK_REAL));
  phi0 = (CCTK_REAL *) malloc(N_points * sizeof(CCTK_REAL));
  W    = (CCTK_REAL *) malloc(N_points * sizeof(CCTK_REAL));

  dW_dR    = (CCTK_REAL *) malloc(N_points * sizeof(CCTK_REAL));
  dW_dth   = (CCTK_REAL *) malloc(N_points * sizeof(CCTK_REAL));
  d2W_dth2 = (CCTK_REAL *) malloc(N_points * sizeof(CCTK_REAL));

  for (int k = 0; k < cctk_lsh[2]; ++k) {
    for (int j = 0; j < cctk_lsh[1]; ++j) {
      for (int i = 0; i < cctk_lsh[0]; ++i) {

        const CCTK_INT ind  = CCTK_GFINDEX3D (cctkGH, i, j, k);

        const CCTK_REAL x1  = x[ind] - x0;
        const CCTK_REAL y1  = y[ind] - y0;
        const CCTK_REAL z1  = z[ind] - z0;

        const CCTK_REAL rho2 = x1*x1 + y1*y1;
        const CCTK_REAL rho  = sqrt(rho2);

        const CCTK_REAL RR2 = x1*x1 + y1*y1 + z1*z1;
        const CCTK_REAL RR  = sqrt(RR2);

        const CCTK_REAL costh  = z1/RR;

        /* note that there are divisions by RR in the following expressions.
           divisions by zero should be avoided by choosing a non-zero value for
           z0 (for instance) */

        F1[ind] = (1.0/2.0)*log(pow(-1 + 16*ct/(RR*pow(4 + rH/RR, 2)), 2) + 256*ct*(ct - rH)*pow(z1, 2)/(pow(RR, 4)*pow(4 + rH/RR, 4)));

        F2[ind] = -F1[ind] + (1.0/2.0)*log(pow(pow(-1 + 16*ct/(RR*pow(4 + rH/RR, 2)), 2) + 256*ct*(ct - rH)/(pow(RR, 2)*pow(4 + rH/RR, 4)), 2) - 256*ct*pow(4*RR - rH, 2)*(ct - rH)*rho2/(pow(RR, 4)*pow(4 + rH/RR, 4)*pow(4*RR + rH, 2)));

        F0[ind] = -F2[ind];

        W[ind] = 4096*sqrt(ct*(ct - rH))*(-1 + 16*ct/(RR*pow(4 + rH/RR, 2)))*(2*ct - rH)*exp(-2*F1[ind] - 2*F2[ind])/(pow(RR, 3)*pow(4 + rH/RR, 6));

        dW_dR[ind]    =  4096*pow(RR, 2)*sqrt(ct*(ct - rH))*(4*RR - rH)*(4*RR + rH)*(2*ct - rH)*(196608*pow(RR, 8) - 786432*pow(RR, 7)*ct + 393216*pow(RR, 7)*rH + 1310720*pow(RR, 6)*pow(ct, 2) - 1310720*pow(RR, 6)*ct*rH + 344064*pow(RR, 6)*pow(rH, 2) - 1048576*pow(RR, 5)*pow(ct, 3) + 1572864*pow(RR, 5)*pow(ct, 2)*rH - 868352*pow(RR, 5)*ct*pow(rH, 2) + 172032*pow(RR, 5)*pow(rH, 3) + 262144*pow(RR, 4)*pow(ct, 4) - 524288*pow(RR, 4)*pow(ct, 3)*rH + 557056*pow(RR, 4)*pow(ct, 2)*pow(rH, 2) - 65536*pow(RR, 4)*pow(ct, 2)*pow(rho, 2) - 294912*pow(RR, 4)*ct*pow(rH, 3) + 65536*pow(RR, 4)*ct*rH*pow(rho, 2) + 53760*pow(RR, 4)*pow(rH, 4) - 65536*pow(RR, 3)*pow(ct, 3)*pow(rH, 2) + 131072*pow(RR, 3)*pow(ct, 3)*pow(rho, 2) + 98304*pow(RR, 3)*pow(ct, 2)*pow(rH, 3) - 196608*pow(RR, 3)*pow(ct, 2)*rH*pow(rho, 2) - 54272*pow(RR, 3)*ct*pow(rH, 4) + 65536*pow(RR, 3)*ct*pow(rH, 2)*pow(rho, 2) + 10752*pow(RR, 3)*pow(rH, 5) + 5120*pow(RR, 2)*pow(ct, 2)*pow(rH, 4) - 24576*pow(RR, 2)*pow(ct, 2)*pow(rH, 2)*pow(rho, 2) - 5120*pow(RR, 2)*ct*pow(rH, 5) + 24576*pow(RR, 2)*ct*pow(rH, 3)*pow(rho, 2) + 1344*pow(RR, 2)*pow(rH, 6) + 8192*RR*pow(ct, 3)*pow(rH, 2)*pow(rho, 2) - 12288*RR*pow(ct, 2)*pow(rH, 3)*pow(rho, 2) - 192*RR*ct*pow(rH, 6) + 4096*RR*ct*pow(rH, 4)*pow(rho, 2) + 96*RR*pow(rH, 7) - 256*pow(ct, 2)*pow(rH, 4)*pow(rho, 2) + 256*ct*pow(rH, 5)*pow(rho, 2) + 3*pow(rH, 8))/pow(65536*pow(RR, 8) - 262144*pow(RR, 7)*ct + 131072*pow(RR, 7)*rH + 524288*pow(RR, 6)*pow(ct, 2) - 524288*pow(RR, 6)*ct*rH + 114688*pow(RR, 6)*pow(rH, 2) - 524288*pow(RR, 5)*pow(ct, 3) + 786432*pow(RR, 5)*pow(ct, 2)*rH - 376832*pow(RR, 5)*ct*pow(rH, 2) + 57344*pow(RR, 5)*pow(rH, 3) + 262144*pow(RR, 4)*pow(ct, 4) - 524288*pow(RR, 4)*pow(ct, 3)*rH + 393216*pow(RR, 4)*pow(ct, 2)*pow(rH, 2) - 65536*pow(RR, 4)*pow(ct, 2)*pow(rho, 2) - 131072*pow(RR, 4)*ct*pow(rH, 3) + 65536*pow(RR, 4)*ct*rH*pow(rho, 2) + 17920*pow(RR, 4)*pow(rH, 4) - 32768*pow(RR, 3)*pow(ct, 3)*pow(rH, 2) + 49152*pow(RR, 3)*pow(ct, 2)*pow(rH, 3) - 23552*pow(RR, 3)*ct*pow(rH, 4) + 3584*pow(RR, 3)*pow(rH, 5) + 2048*pow(RR, 2)*pow(ct, 2)*pow(rH, 4) + 8192*pow(RR, 2)*pow(ct, 2)*pow(rH, 2)*pow(rho, 2) - 2048*pow(RR, 2)*ct*pow(rH, 5) - 8192*pow(RR, 2)*ct*pow(rH, 3)*pow(rho, 2) + 448*pow(RR, 2)*pow(rH, 6) - 64*RR*ct*pow(rH, 6) + 32*RR*pow(rH, 7) - 256*pow(ct, 2)*pow(rH, 4)*pow(rho, 2) + 256*ct*pow(rH, 5)*pow(rho, 2) + pow(rH, 8), 2);

        dW_dth[ind]   = -2097152*pow(RR, 4)*ct*rho*sqrt(ct*(ct - rH))*pow(4*RR - rH, 2)*pow(4*RR + rH, 2)*(ct - rH)*(2*ct - rH)*(16*pow(RR, 2) - 16*RR*ct + 8*RR*rH + pow(rH, 2))*costh/pow(65536*pow(RR, 8) - 262144*pow(RR, 7)*ct + 131072*pow(RR, 7)*rH + 524288*pow(RR, 6)*pow(ct, 2) - 524288*pow(RR, 6)*ct*rH + 114688*pow(RR, 6)*pow(rH, 2) - 524288*pow(RR, 5)*pow(ct, 3) + 786432*pow(RR, 5)*pow(ct, 2)*rH - 376832*pow(RR, 5)*ct*pow(rH, 2) + 57344*pow(RR, 5)*pow(rH, 3) + 262144*pow(RR, 4)*pow(ct, 4) - 524288*pow(RR, 4)*pow(ct, 3)*rH + 393216*pow(RR, 4)*pow(ct, 2)*pow(rH, 2) - 65536*pow(RR, 4)*pow(ct, 2)*pow(rho, 2) - 131072*pow(RR, 4)*ct*pow(rH, 3) + 65536*pow(RR, 4)*ct*rH*pow(rho, 2) + 17920*pow(RR, 4)*pow(rH, 4) - 32768*pow(RR, 3)*pow(ct, 3)*pow(rH, 2) + 49152*pow(RR, 3)*pow(ct, 2)*pow(rH, 3) - 23552*pow(RR, 3)*ct*pow(rH, 4) + 3584*pow(RR, 3)*pow(rH, 5) + 2048*pow(RR, 2)*pow(ct, 2)*pow(rH, 4) + 8192*pow(RR, 2)*pow(ct, 2)*pow(rH, 2)*pow(rho, 2) - 2048*pow(RR, 2)*ct*pow(rH, 5) - 8192*pow(RR, 2)*ct*pow(rH, 3)*pow(rho, 2) + 448*pow(RR, 2)*pow(rH, 6) - 64*RR*ct*pow(rH, 6) + 32*RR*pow(rH, 7) - 256*pow(ct, 2)*pow(rH, 4)*pow(rho, 2) + 256*ct*pow(rH, 5)*pow(rho, 2) + pow(rH, 8), 2);

        d2W_dth2[ind] =  -2097152*pow(RR, 5)*ct*sqrt(ct*(ct - rH))*pow(4*RR - rH, 2)*pow(4*RR + rH, 2)*(ct - rH)*(2*ct - rH)*(16*pow(RR, 2) - 16*RR*ct + 8*RR*rH + pow(rH, 2))*(65536*pow(RR, 8)*pow(costh, 2) - 262144*pow(RR, 7)*ct*pow(costh, 2) + 131072*pow(RR, 7)*rH*pow(costh, 2) + 524288*pow(RR, 6)*pow(ct, 2)*pow(costh, 2) - 524288*pow(RR, 6)*ct*rH*pow(costh, 2) + 114688*pow(RR, 6)*pow(rH, 2)*pow(costh, 2) - 65536*pow(RR, 6)*pow(rho, 2) - 524288*pow(RR, 5)*pow(ct, 3)*pow(costh, 2) + 786432*pow(RR, 5)*pow(ct, 2)*rH*pow(costh, 2) - 376832*pow(RR, 5)*ct*pow(rH, 2)*pow(costh, 2) + 262144*pow(RR, 5)*ct*pow(rho, 2) + 57344*pow(RR, 5)*pow(rH, 3)*pow(costh, 2) - 131072*pow(RR, 5)*rH*pow(rho, 2) + 262144*pow(RR, 4)*pow(ct, 4)*pow(costh, 2) - 524288*pow(RR, 4)*pow(ct, 3)*rH*pow(costh, 2) + 393216*pow(RR, 4)*pow(ct, 2)*pow(rH, 2)*pow(costh, 2) + 196608*pow(RR, 4)*pow(ct, 2)*pow(rho, 2)*pow(costh, 2) - 524288*pow(RR, 4)*pow(ct, 2)*pow(rho, 2) - 131072*pow(RR, 4)*ct*pow(rH, 3)*pow(costh, 2) - 196608*pow(RR, 4)*ct*rH*pow(rho, 2)*pow(costh, 2) + 524288*pow(RR, 4)*ct*rH*pow(rho, 2) + 17920*pow(RR, 4)*pow(rH, 4)*pow(costh, 2) - 114688*pow(RR, 4)*pow(rH, 2)*pow(rho, 2) - 32768*pow(RR, 3)*pow(ct, 3)*pow(rH, 2)*pow(costh, 2) + 524288*pow(RR, 3)*pow(ct, 3)*pow(rho, 2) + 49152*pow(RR, 3)*pow(ct, 2)*pow(rH, 3)*pow(costh, 2) - 786432*pow(RR, 3)*pow(ct, 2)*rH*pow(rho, 2) - 23552*pow(RR, 3)*ct*pow(rH, 4)*pow(costh, 2) + 376832*pow(RR, 3)*ct*pow(rH, 2)*pow(rho, 2) + 3584*pow(RR, 3)*pow(rH, 5)*pow(costh, 2) - 57344*pow(RR, 3)*pow(rH, 3)*pow(rho, 2) - 262144*pow(RR, 2)*pow(ct, 4)*pow(rho, 2) + 524288*pow(RR, 2)*pow(ct, 3)*rH*pow(rho, 2) + 2048*pow(RR, 2)*pow(ct, 2)*pow(rH, 4)*pow(costh, 2) - 24576*pow(RR, 2)*pow(ct, 2)*pow(rH, 2)*pow(rho, 2)*pow(costh, 2) - 393216*pow(RR, 2)*pow(ct, 2)*pow(rH, 2)*pow(rho, 2) + 65536*pow(RR, 2)*pow(ct, 2)*pow(rho, 4) - 2048*pow(RR, 2)*ct*pow(rH, 5)*pow(costh, 2) + 24576*pow(RR, 2)*ct*pow(rH, 3)*pow(rho, 2)*pow(costh, 2) + 131072*pow(RR, 2)*ct*pow(rH, 3)*pow(rho, 2) - 65536*pow(RR, 2)*ct*rH*pow(rho, 4) + 448*pow(RR, 2)*pow(rH, 6)*pow(costh, 2) - 17920*pow(RR, 2)*pow(rH, 4)*pow(rho, 2) + 32768*RR*pow(ct, 3)*pow(rH, 2)*pow(rho, 2) - 49152*RR*pow(ct, 2)*pow(rH, 3)*pow(rho, 2) - 64*RR*ct*pow(rH, 6)*pow(costh, 2) + 23552*RR*ct*pow(rH, 4)*pow(rho, 2) + 32*RR*pow(rH, 7)*pow(costh, 2) - 3584*RR*pow(rH, 5)*pow(rho, 2) + 768*pow(ct, 2)*pow(rH, 4)*pow(rho, 2)*pow(costh, 2) - 2048*pow(ct, 2)*pow(rH, 4)*pow(rho, 2) - 8192*pow(ct, 2)*pow(rH, 2)*pow(rho, 4) - 768*ct*pow(rH, 5)*pow(rho, 2)*pow(costh, 2) + 2048*ct*pow(rH, 5)*pow(rho, 2) + 8192*ct*pow(rH, 3)*pow(rho, 4) + pow(rH, 8)*pow(costh, 2) - 448*pow(rH, 6)*pow(rho, 2) + 64*ct*pow(rH, 6)*pow(rho, 2)/RR - 32*pow(rH, 7)*pow(rho, 2)/RR + 256*pow(ct, 2)*pow(rH, 4)*pow(rho, 4)/pow(RR, 2) - 256*ct*pow(rH, 5)*pow(rho, 4)/pow(RR, 2) - pow(rH, 8)*pow(rho, 2)/pow(RR, 2))/pow(65536*pow(RR, 8) - 262144*pow(RR, 7)*ct + 131072*pow(RR, 7)*rH + 524288*pow(RR, 6)*pow(ct, 2) - 524288*pow(RR, 6)*ct*rH + 114688*pow(RR, 6)*pow(rH, 2) - 524288*pow(RR, 5)*pow(ct, 3) + 786432*pow(RR, 5)*pow(ct, 2)*rH - 376832*pow(RR, 5)*ct*pow(rH, 2) + 57344*pow(RR, 5)*pow(rH, 3) + 262144*pow(RR, 4)*pow(ct, 4) - 524288*pow(RR, 4)*pow(ct, 3)*rH + 393216*pow(RR, 4)*pow(ct, 2)*pow(rH, 2) - 65536*pow(RR, 4)*pow(ct, 2)*pow(rho, 2) - 131072*pow(RR, 4)*ct*pow(rH, 3) + 65536*pow(RR, 4)*ct*rH*pow(rho, 2) + 17920*pow(RR, 4)*pow(rH, 4) - 32768*pow(RR, 3)*pow(ct, 3)*pow(rH, 2) + 49152*pow(RR, 3)*pow(ct, 2)*pow(rH, 3) - 23552*pow(RR, 3)*ct*pow(rH, 4) + 3584*pow(RR, 3)*pow(rH, 5) + 2048*pow(RR, 2)*pow(ct, 2)*pow(rH, 4) + 8192*pow(RR, 2)*pow(ct, 2)*pow(rH, 2)*pow(rho, 2) - 2048*pow(RR, 2)*ct*pow(rH, 5) - 8192*pow(RR, 2)*ct*pow(rH, 3)*pow(rho, 2) + 448*pow(RR, 2)*pow(rH, 6) - 64*RR*ct*pow(rH, 6) + 32*RR*pow(rH, 7) - 256*pow(ct, 2)*pow(rH, 4)*pow(rho, 2) + 256*ct*pow(rH, 5)*pow(rho, 2) + pow(rH, 8), 3);

        phi0[ind] = 0.;

      }
    }
  }

  /* printf("F1 = %g\n", F1[0]); */
  /* printf("F2 = %g\n", F2[0]); */
  /* printf("F0 = %g\n", F0[0]); */
  /* printf("phi0 = %g\n", phi0[0]); */
  /* printf("W = %g\n", W[0]); */

  /* now we finally write the metric and all 3+1 quantities. first we write the
     3-metric, lapse and scalar fields */

  const CCTK_REAL tt = cctk_time;

  for (int k = 0; k < cctk_lsh[2]; ++k) {
    for (int j = 0; j < cctk_lsh[1]; ++j) {
      for (int i = 0; i < cctk_lsh[0]; ++i) {

        const CCTK_INT ind  = CCTK_GFINDEX3D (cctkGH, i, j, k);

        const CCTK_REAL x1  = x[ind] - x0;
        const CCTK_REAL y1  = y[ind] - y0;
        const CCTK_REAL z1  = z[ind] - z0;

        const CCTK_REAL RR2 = x1*x1 + y1*y1 + z1*z1;
        const CCTK_REAL RR  = sqrt(RR2);

        const CCTK_REAL ph = atan2(y1, x1);

        const CCTK_REAL cosph  = cos(ph);
        const CCTK_REAL sinph  = sin(ph);

        const CCTK_REAL cosmph = cos(mm*ph);
        const CCTK_REAL sinmph = sin(mm*ph);

        /* note the division by RR in the following. divisions by zero should be
           avoided by choosing a non-zero value for z0 (for instance) */
        const CCTK_REAL aux  = 1. + 0.25 * rH/RR;
        const CCTK_REAL aux4 = aux * aux * aux * aux;
        const CCTK_REAL psi4 = exp(2. * F1[ind]) * aux4;
        const CCTK_REAL psi2 = sqrt(psi4);
        const CCTK_REAL psi1 = sqrt(psi2);

        const CCTK_REAL h_rho2 = exp(2. * (F2[ind] - F1[ind])) - 1.;

        // 3-metric
        gxx[ind] = psi4 * (1. + h_rho2 * sinph * sinph);
        gxy[ind] = -psi4 * h_rho2 * sinph * cosph;
        gxz[ind] = 0;
        gyy[ind] = psi4 * (1. + h_rho2 * cosph * cosph);
        gyz[ind] = 0;
        gzz[ind] = psi4;

        CCTK_REAL alph = exp(F0[ind]) * (RR - 0.25*rH) / (RR + 0.25*rH);
        // FIXME
        if (alph < SMALL)
          alph = SMALL;

        // lapse
        if (CCTK_EQUALS(initial_lapse, "psi^n"))
          alp[ind] = pow(psi1, initial_lapse_psi_exponent);
        else if (CCTK_EQUALS(initial_lapse, "Kerr_test"))
          alp[ind] = alph;

        // shift
        if (CCTK_EQUALS(initial_shift, "Kerr_test")) {
          betax[ind] =  W[ind] * y1;
          betay[ind] = -W[ind] * x1;
          betaz[ind] =  0.;
        }


        // add perturbation
        CCTK_REAL phi0_l = phi0[ind];
        phi0_l *= 1. + pert_A * exp( -0.5*RR2/(pert_Rmax*pert_Rmax) )
                              * sin(2.*M_PI * RR / pert_lambda);

        const CCTK_REAL omega = mm * OmegaH;

        // scalar fields
        phi1[ind]  = phi0_l * (cos(omega * tt) * cosmph + sin(omega * tt) * sinmph);
        phi2[ind]  = phi0_l * (cos(omega * tt) * sinmph - sin(omega * tt) * cosmph);
        Kphi1[ind] = 0.5 * mm * (W[ind] - OmegaH) / alph * phi2[ind];
        Kphi2[ind] = 0.5 * mm * (OmegaH - W[ind]) / alph * phi1[ind];

      } /* for i */
    }   /* for j */
  }     /* for k */


  /* now we write the extrinsic curvature. since there are derivatives, we can
     only loop through the interior points. */

  for (int k = cctk_nghostzones[2]; k < cctk_lsh[2] - cctk_nghostzones[2]; k++)
    for (int j = cctk_nghostzones[1]; j < cctk_lsh[1] - cctk_nghostzones[1]; j++)
      for (int i = cctk_nghostzones[0]; i < cctk_lsh[0] - cctk_nghostzones[0]; i++) {

        const CCTK_INT ind  = CCTK_GFINDEX3D(cctkGH, i, j, k);

        const CCTK_REAL x1  = x[ind] - x0;
        const CCTK_REAL y1  = y[ind] - y0;
        const CCTK_REAL z1  = z[ind] - z0;

        const CCTK_REAL RR2 = x1*x1 + y1*y1 + z1*z1;
        const CCTK_REAL RR  = sqrt(RR2);

        const CCTK_REAL rho2 = x1*x1 + y1*y1;
        const CCTK_REAL rho  = sqrt(rho2);

        const CCTK_REAL costh  = z1/RR;
        const CCTK_REAL costh2 = costh*costh;
        const CCTK_REAL sinth2 = 1. - costh2;
        const CCTK_REAL sinth  = sqrt(sinth2);

        // R dW/dR
        CCTK_REAL RdWdR = RR * dW_dR[ind];

        // R sin(th) dW/dth
        CCTK_REAL RsinthdWdth = RR * sinth * dW_dth[ind];

        /*
        const CCTK_REAL R_x = x1/RR;   // dR/dx
        const CCTK_REAL R_y = y1/RR;   // dR/dy
        const CCTK_REAL R_z = z1/RR;   // dR/dz
        */

        const CCTK_REAL sinth2ph_x = -y1/RR2; // sin(th)^2 dphi/dx
        const CCTK_REAL sinth2ph_y =  x1/RR2; // sin(th)^2 dphi/dy

        const CCTK_REAL Rsinthth_x  = z1*x1/RR2; // R sin(th) dth/dx
        const CCTK_REAL Rsinthth_y  = z1*y1/RR2; // R sin(th) dth/dy
        const CCTK_REAL Rsinthth_z  = -sinth2;   // R sin(th) dth/dz

        CCTK_REAL alph = exp(F0[ind]) * (RR - 0.25*rH) / (RR + 0.25*rH);
        // FIXME
        if (alph < SMALL)
          alph = SMALL;

        const CCTK_REAL aux  = 1. + 0.25 * rH/RR;
        const CCTK_REAL aux4 = aux * aux * aux * aux;

        /*
          KRph/(R sin(th)^2) = - 1/2 exp(2F2) (1 + rH/(4R))^4 / alpha R dW/dR
          Kthph/(R sin(th))  = - 1/2 exp(2F2) (1 + rH/(4R))^4 / alpha R sin(th) dW/dth

          note that in the expressions below we absorb one R factor in the terms
          RdWdR and RsinthdWdth
        */
        const CCTK_REAL KRph_o_Rsinth2 = -0.5 * exp(2. * F2[ind]) * aux4 / alph * RdWdR;

        CCTK_REAL RsinthdWdth_o_sinth2;
        // if at the rho = 0 axis we need to regularize the division by sin(th)^2
        if (rho < sqrt(dxsq + dysq) * 0.25) {
          // R d^2W/dth^2
          CCTK_REAL Rd2th_W = RR * d2W_dth2[ind];

          RsinthdWdth_o_sinth2 = Rd2th_W;
        }
        else {
          RsinthdWdth_o_sinth2 = RsinthdWdth / sinth2;
        }

        const CCTK_REAL Kthph_o_Rsinth3 = -0.5 * exp(2. * F2[ind]) * aux4 / alph * RsinthdWdth_o_sinth2;


        kxx[ind] = 2.*KRph_o_Rsinth2 *  x1 * sinth2ph_x                     +  2.*Kthph_o_Rsinth3 *  Rsinthth_x * sinth2ph_x;
        kxy[ind] =    KRph_o_Rsinth2 * (x1 * sinth2ph_y + y1 * sinth2ph_x)  +     Kthph_o_Rsinth3 * (Rsinthth_x * sinth2ph_y + Rsinthth_y * sinth2ph_x);
        kxz[ind] =    KRph_o_Rsinth2 *                    z1 * sinth2ph_x   +     Kthph_o_Rsinth3 *                            Rsinthth_z * sinth2ph_x;
        kyy[ind] = 2.*KRph_o_Rsinth2 *  y1 * sinth2ph_y                     +  2.*Kthph_o_Rsinth3 *  Rsinthth_y * sinth2ph_y;
        kyz[ind] =    KRph_o_Rsinth2 *                    z1 * sinth2ph_y   +     Kthph_o_Rsinth3 *                            Rsinthth_z * sinth2ph_y;
        kzz[ind] = 0.;

      }

  /* now we use the function ExtrapolateGammas, from NewRad thorn, to
     extrapolate the extrinsic curvature to the outer boundary points. inner
     boundaries will be filled after synchronisation. */

  ExtrapolateGammas(cctkGH, kxx);
  ExtrapolateGammas(cctkGH, kxy);
  ExtrapolateGammas(cctkGH, kxz);
  ExtrapolateGammas(cctkGH, kyy);
  ExtrapolateGammas(cctkGH, kyz);
  ExtrapolateGammas(cctkGH, kzz);

  free(F1); free(F2); free(F0); free(phi0); free(W);
  free(dW_dR); free(dW_dth); free(d2W_dth2);

  return;
}
