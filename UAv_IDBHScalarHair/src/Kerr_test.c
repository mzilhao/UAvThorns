
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

static void apply_jacobian(CCTK_REAL dvar[3], CCTK_REAL jac[3][3])
{
  CCTK_REAL xdvar[3];

  for (int a = 0; a < 3; a++)
    xdvar[a] = 0.0;

  for (int a = 0; a < 3; a++)
    for (int b = 0; b < 3; b++) {
      xdvar[a] += dvar[b] * jac[b][a];
    }

  for (int a = 0; a < 3; a++)
    dvar[a] = xdvar[a];

  return;
}

void UAv_Kerr_test(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL dx12 = 12 * CCTK_DELTA_SPACE(0);
  CCTK_REAL dy12 = 12 * CCTK_DELTA_SPACE(1);
  CCTK_REAL dz12 = 12 * CCTK_DELTA_SPACE(2);

  bool use_jacobian = false;
  if (CCTK_IsFunctionAliased("MultiPatch_GetDomainSpecification"))
    use_jacobian = true;

  const CCTK_REAL *lJ11 =
    use_jacobian ? CCTK_VarDataPtr(cctkGH, 0, "Coordinates::J11") : NULL;
  const CCTK_REAL *lJ12 =
    use_jacobian ? CCTK_VarDataPtr(cctkGH, 0, "Coordinates::J12") : NULL;
  const CCTK_REAL *lJ13 =
    use_jacobian ? CCTK_VarDataPtr(cctkGH, 0, "Coordinates::J13") : NULL;
  const CCTK_REAL *lJ21 =
    use_jacobian ? CCTK_VarDataPtr(cctkGH, 0, "Coordinates::J21") : NULL;
  const CCTK_REAL *lJ22 =
    use_jacobian ? CCTK_VarDataPtr(cctkGH, 0, "Coordinates::J22") : NULL;
  const CCTK_REAL *lJ23 =
    use_jacobian ? CCTK_VarDataPtr(cctkGH, 0, "Coordinates::J23") : NULL;
  const CCTK_REAL *lJ31 =
    use_jacobian ? CCTK_VarDataPtr(cctkGH, 0, "Coordinates::J31") : NULL;
  const CCTK_REAL *lJ32 =
    use_jacobian ? CCTK_VarDataPtr(cctkGH, 0, "Coordinates::J32") : NULL;
  const CCTK_REAL *lJ33 =
    use_jacobian ? CCTK_VarDataPtr(cctkGH, 0, "Coordinates::J33") : NULL;


  /* let's create arrays for the F1, F2, F0, phi0, W functions */
  const CCTK_INT N_points = cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; // total points

  CCTK_REAL *F1, *F2, *F0, *phi0, *W;

  F1   = (CCTK_REAL *) malloc(N_points * sizeof(CCTK_REAL));
  F2   = (CCTK_REAL *) malloc(N_points * sizeof(CCTK_REAL));
  F0   = (CCTK_REAL *) malloc(N_points * sizeof(CCTK_REAL));
  phi0 = (CCTK_REAL *) malloc(N_points * sizeof(CCTK_REAL));
  W    = (CCTK_REAL *) malloc(N_points * sizeof(CCTK_REAL));

  for (int k = 0; k < cctk_lsh[2]; ++k) {
    for (int j = 0; j < cctk_lsh[1]; ++j) {
      for (int i = 0; i < cctk_lsh[0]; ++i) {

        const CCTK_INT ind  = CCTK_GFINDEX3D (cctkGH, i, j, k);

        const CCTK_REAL x1  = x[ind] - x0;
        const CCTK_REAL y1  = y[ind] - y0;
        const CCTK_REAL z1  = z[ind] - z0;

        const CCTK_REAL RR2 = x1*x1 + y1*y1 + z1*z1;

        CCTK_REAL RR  = sqrt(RR2);
        // prevent divisions from 0
        if (RR < SMALL)
          RR = SMALL;

        // from (quasi-)isotropic coordinate R to the metric coordinate r
        /* const CCTK_REAL rr = RR * (1. + 0.25 * rH / RR) * (1. + 0.25 * rH / RR); */

        const CCTK_REAL th = acos( z1/RR );

        F1[ind] = (1.0/2.0)*log(pow(-1 + 16*ct/(RR*pow(4 + rH/RR, 2)), 2) + 256*ct*(ct - rH)*pow(cos(th), 2)/(pow(RR, 2)*pow(4 + rH/RR, 4)));

        F2[ind] = -F1[ind] + (1.0/2.0)*log(pow(pow(-1 + 16*ct/(RR*pow(4 + rH/RR, 2)), 2) + 256*ct*(ct - rH)/(pow(RR, 2)*pow(4 + rH/RR, 4)), 2) - 256*ct*pow(4*RR - rH, 2)*(ct - rH)*pow(sin(th), 2)/(pow(RR, 2)*pow(4 + rH/RR, 4)*pow(4*RR + rH, 2)));

        F0[ind] = -F2[ind];

        W[ind] = 4096*sqrt(ct*(ct - rH))*(-1 + 16*ct/(RR*pow(4 + rH/RR, 2)))*(2*ct - rH)*exp(-2*F1[ind] - 2*F2[ind])/(pow(RR, 3)*pow(4 + rH/RR, 6));

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

        CCTK_REAL RR2 = x1*x1 + y1*y1 + z1*z1;
        if (RR2 < pow(SMALL, 2))
          RR2 = pow(SMALL, 2);
        const CCTK_REAL RR  = sqrt(RR2);

        CCTK_REAL rho2 = x1*x1 + y1*y1;
        if (rho2 < pow(SMALL, 2))
          rho2 = pow(SMALL, 2);
        const CCTK_REAL rho  = sqrt(rho2);

        const CCTK_REAL cosph  = x1/rho;
        const CCTK_REAL sinph  = y1/rho;

        CCTK_REAL cosmph = 0.;
        CCTK_REAL sinmph = 0.;
        if (mm == 0) {
          cosmph = 1.;
          sinmph = 0.;
        } else if (mm == 1) {
          cosmph = cosph;
          sinmph = sinph;
        } else if (mm == 2) {
          cosmph = cosph * cosph - sinph * sinph;
          sinmph = 2. * cosph * sinph;
        } else if (mm == 3) {
          cosmph = cosph * cosph * cosph - 3. * cosph * sinph * sinph;
          sinmph = 3. * cosph * cosph * sinph - sinph * sinph * sinph;
        } else {
          CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "mm = %d not implemented yet! Aborting.", mm);
        }

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

        const CCTK_INT ind     = CCTK_GFINDEX3D(cctkGH, i, j, k);

        const CCTK_INT indim1  = CCTK_GFINDEX3D(cctkGH, i-1, j, k);
        const CCTK_INT indip1  = CCTK_GFINDEX3D(cctkGH, i+1, j, k);
        const CCTK_INT indim2  = CCTK_GFINDEX3D(cctkGH, i-2, j, k);
        const CCTK_INT indip2  = CCTK_GFINDEX3D(cctkGH, i+2, j, k);

        const CCTK_INT indjm1  = CCTK_GFINDEX3D(cctkGH, i, j-1, k);
        const CCTK_INT indjp1  = CCTK_GFINDEX3D(cctkGH, i, j+1, k);
        const CCTK_INT indjm2  = CCTK_GFINDEX3D(cctkGH, i, j-2, k);
        const CCTK_INT indjp2  = CCTK_GFINDEX3D(cctkGH, i, j+2, k);

        const CCTK_INT indkm1  = CCTK_GFINDEX3D(cctkGH, i, j, k-1);
        const CCTK_INT indkp1  = CCTK_GFINDEX3D(cctkGH, i, j, k+1);
        const CCTK_INT indkm2  = CCTK_GFINDEX3D(cctkGH, i, j, k-2);
        const CCTK_INT indkp2  = CCTK_GFINDEX3D(cctkGH, i, j, k+2);

        const CCTK_REAL x1  = x[ind] - x0;
        const CCTK_REAL y1  = y[ind] - y0;
        const CCTK_REAL z1  = z[ind] - z0;

        CCTK_REAL RR2 = x1*x1 + y1*y1 + z1*z1;
        if (RR2 < pow(SMALL, 2))
          RR2 = pow(SMALL, 2);
        const CCTK_REAL RR  = sqrt(RR2);

        CCTK_REAL rho2 = x1*x1 + y1*y1;
        if (rho2 < pow(SMALL, 2))
          rho2 = pow(SMALL, 2);
        /* const CCTK_REAL rho  = sqrt(rho2); */

        /* get the jacobian to take derivatives. 0 is x, 1 is y and 2 is z */
        CCTK_REAL jac[3][3];
        if (use_jacobian) {
          jac[0][0] = lJ11[ind];
          jac[0][1] = lJ12[ind];
          jac[0][2] = lJ13[ind];
          jac[1][0] = lJ21[ind];
          jac[1][1] = lJ22[ind];
          jac[1][2] = lJ23[ind];
          jac[2][0] = lJ31[ind];
          jac[2][1] = lJ32[ind];
          jac[2][2] = lJ33[ind];
        } else {
          jac[0][0] = 1.0;
          jac[1][1] = 1.0;
          jac[2][2] = 1.0;
          jac[0][1] = 0.0;
          jac[0][2] = 0.0;
          jac[1][0] = 0.0;
          jac[1][2] = 0.0;
          jac[2][0] = 0.0;
          jac[2][1] = 0.0;
        }

        // 4th-order accurate first derivatives of the W function
        CCTK_REAL d1_W[3];

        // dW/dx
        d1_W[0] = (   -W[indip2] + 8*W[indip1]
                   - 8*W[indim1] +   W[indim2] ) / dx12;

        // dW/dy
        d1_W[1] = (   -W[indjp2] + 8*W[indjp1]
                   - 8*W[indjm1] +   W[indjm2] ) / dy12;

        // dW/dz
        d1_W[2] = (   -W[indkp2] + 8*W[indkp1]
                   - 8*W[indkm1] +   W[indkm2] ) / dz12;

        if (use_jacobian) {
          apply_jacobian(d1_W, jac);
        }

        // R dW/dR
        CCTK_REAL RdWdR = x1 * d1_W[0] + y1 * d1_W[1] + z1 * d1_W[2];

        // R sin(th) dW/dth
        CCTK_REAL RsinthdWdth = x1 * z1 * d1_W[0] + y1 * z1 * d1_W[1]
                              - (x1 * x1 + y1 * y1) * d1_W[2];

        const CCTK_REAL R_x = x1/RR;   // dR/dx
        const CCTK_REAL R_y = y1/RR;   // dR/dy
        const CCTK_REAL R_z = z1/RR;   // dR/dz

        const CCTK_REAL costh  = z1/RR;
        const CCTK_REAL costh2 = costh*costh;
        const CCTK_REAL sinth2 = 1. - costh2;

        const CCTK_REAL ph_x = -y1/rho2; // dphi/dx
        const CCTK_REAL ph_y =  x1/rho2; // dphi/dy

        const CCTK_REAL sinth2ph_x = -y1/RR2; // sin(th)^2 dphi/dx
        const CCTK_REAL sinth2ph_y =  x1/RR2; // sin(th)^2 dphi/dy

        const CCTK_REAL sinthth_x  = z1*x1/(RR*RR2); // sin(th) dth/dx
        const CCTK_REAL sinthth_y  = z1*y1/(RR*RR2); // sin(th) dth/dy
        const CCTK_REAL sinthth_z  = -sinth2/RR;     // sin(th) dth/dz

        CCTK_REAL alph = exp(F0[ind]) * (RR - 0.25*rH) / (RR + 0.25*rH);
        if (alph < SMALL)
          alph = SMALL;

        const CCTK_REAL aux  = 1. + 0.25 * rH/RR;
        const CCTK_REAL aux4 = aux * aux * aux * aux;

        // KRph/sin(th)^2 = - 1/2 R^2 exp(2F2) (1 + rH/(4R))^4 / alpha  dW/dR
        // Kthph/sin(th)  = - 1/2 R^2 exp(2F2) (1 + rH/(4R))^4 / alpha sin(th) dW/dth
        CCTK_REAL KRph_o_sinth2 = -0.5 * RR * exp(2. * F2[ind]) * aux4 / alph * RdWdR;
        CCTK_REAL Kthph_o_sinth = -0.5 * RR * exp(2. * F2[ind]) * aux4 / alph * RsinthdWdth;

        kxx[ind] = 2.*KRph_o_sinth2 *  R_x * sinth2ph_x                     +  2.*Kthph_o_sinth *  sinthth_x * ph_x;
        kxy[ind] =    KRph_o_sinth2 * (R_x * sinth2ph_y + R_y * sinth2ph_x) +     Kthph_o_sinth * (sinthth_x * ph_y + sinthth_y * ph_x);
        kxz[ind] =    KRph_o_sinth2 *                     R_z * sinth2ph_x  +     Kthph_o_sinth *                     sinthth_z * ph_x;
        kyy[ind] = 2.*KRph_o_sinth2 *  R_y * sinth2ph_y                     +  2.*Kthph_o_sinth *  sinthth_y * ph_y;
        kyz[ind] =    KRph_o_sinth2 *                     R_z * sinth2ph_y  +     Kthph_o_sinth *                     sinthth_z * ph_y;
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
}
