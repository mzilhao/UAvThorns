
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

void UAv_ID_read_data(CCTK_INT *, CCTK_INT *, CCTK_REAL [], CCTK_REAL [],
                   CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], CCTK_REAL []);


void UAv_IDBHScalarHair(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dxsq = CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(0);
  const CCTK_REAL dysq = CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(1);
  const CCTK_REAL dzsq = CCTK_DELTA_SPACE(2)*CCTK_DELTA_SPACE(2);

  CCTK_INT NF;      // NF will be the actual size of the arrays
  CCTK_INT NX;      // NX will be the number of X points
  CCTK_INT Ntheta;  // Ntheta will be the number of theta points

  CCTK_REAL *Xtmp, *thtmp, *F1_in, *F2_in, *F0_in, *phi0_in, *W_in;
  Xtmp     = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));
  thtmp    = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));
  F1_in    = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));
  F2_in    = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));
  F0_in    = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));
  phi0_in  = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));
  W_in     = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));

  // we get the data from the input file
  UAv_ID_read_data(&NF, &NX, Xtmp, thtmp, F1_in, F2_in, F0_in, phi0_in, W_in);

  Ntheta = NF/NX;

  CCTK_VInfo(CCTK_THORNSTRING, "NX     = %d", NX);
  CCTK_VInfo(CCTK_THORNSTRING, "Ntheta = %d", Ntheta);
  CCTK_VInfo(CCTK_THORNSTRING, "NF     = %d", NF);

  // now we create arrays with the X and theta coordinates
  CCTK_REAL X[NX], theta[Ntheta];
  for (int i = 0; i < NX; i++) {
    X[i]     = Xtmp[i];
    /* printf("X[%3d] = %lf\n", i, X[i]); */
  }
  for (int i = 0; i < Ntheta; i++) {
    theta[i] = thtmp[i*NX];
    /* printf("theta[%3d] = %lf\n", i, theta[i]); */
  }

  // the spacing in each coordinate is
  const CCTK_REAL dX     = (X[NX-1] - X[0])/(NX-1);
  const CCTK_REAL dtheta = (theta[Ntheta-1] - theta[0])/(Ntheta-1);

  /* printf("dX     = %e\n", dX); */
  /* printf("dtheta = %e\n", dtheta); */

  // make sure spacing is uniform in the provided grid
  for (int i = 1; i < NX; i++) {
    if (fabs(X[i] - X[i - 1] - dX) > SMALL) {
      printf("i = %d\n", i);
      CCTK_WARN(0, "X grid is not uniformly spaced. Aborting.");
    }
  }
  for (int j = 1; j < Ntheta; j++) {
    if (fabs(theta[j] - theta[j - 1] - dtheta) > SMALL)
      CCTK_WARN(0, "theta grid is not uniformly spaced. Aborting.");
  }


  CCTK_REAL *dW_dr_in, *dW_dth_in, *d2W_dth2_in, *d2W_drth_in;
  dW_dr_in    = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  dW_dth_in   = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  d2W_dth2_in = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));
  d2W_drth_in = (CCTK_REAL *) malloc(NF * sizeof(CCTK_REAL));

  // take the derivatives of the input data

  const CCTK_REAL oodX       = 1. / dX;
  const CCTK_REAL oodXsq     = oodX * oodX;
  const CCTK_REAL oodX12     = 1. / (12. * dX);
  const CCTK_REAL oodth12    = 1. / (12. * dtheta);
  const CCTK_REAL oodthsq12  = 1. / (12. * dtheta * dtheta);
  const CCTK_REAL oodXdth4   = 1. / (4.  * dX * dtheta);
  const CCTK_REAL oodXdth144 = 1. / (144. * dX * dtheta);

  for (int j = 0; j < Ntheta; j++) {
    for (int i = 0; i < NX; i++) {

      const CCTK_INT ind    = i + j*NX;

      const CCTK_INT indim1 = i-1 + j*NX;
      const CCTK_INT indip1 = i+1 + j*NX;
      const CCTK_INT indim2 = i-2 + j*NX;
      const CCTK_INT indip2 = i+2 + j*NX;
      const CCTK_INT indip3 = i+3 + j*NX;

      CCTK_INT indim2jm1 = i-2 + (j-1)*NX;
      CCTK_INT indim2jm2 = i-2 + (j-2)*NX;
      CCTK_INT indim2jp1 = i-2 + (j+1)*NX;
      CCTK_INT indim2jp2 = i-2 + (j+2)*NX;

      CCTK_INT indim1jm1 = i-1 + (j-1)*NX;
      CCTK_INT indim1jm2 = i-1 + (j-2)*NX;
      CCTK_INT indim1jp1 = i-1 + (j+1)*NX;
      CCTK_INT indim1jp2 = i-1 + (j+2)*NX;

      CCTK_INT indjm1 = i + (j-1)*NX;
      CCTK_INT indjm2 = i + (j-2)*NX;
      CCTK_INT indjp1 = i + (j+1)*NX;
      CCTK_INT indjp2 = i + (j+2)*NX;

      CCTK_INT indip1jm1 = i+1 + (j-1)*NX;
      CCTK_INT indip1jm2 = i+1 + (j-2)*NX;
      CCTK_INT indip1jp1 = i+1 + (j+1)*NX;
      CCTK_INT indip1jp2 = i+1 + (j+2)*NX;

      CCTK_INT indip2jm1 = i+2 + (j-1)*NX;
      CCTK_INT indip2jm2 = i+2 + (j-2)*NX;
      CCTK_INT indip2jp1 = i+2 + (j+1)*NX;
      CCTK_INT indip2jp2 = i+2 + (j+2)*NX;

      /* let's use the fact that the solution is axi-symmetric (and that
         theta[0] = 0) for the boundary points in j */
      if (j == 0) {
        indjm1    = i   + (j+1)*NX;
        indjm2    = i   + (j+2)*NX;
        indim1jm1 = i-1 + (j+1)*NX;
        indim2jm2 = i-2 + (j+2)*NX;
        indip1jm1 = i+1 + (j+1)*NX;
        indip2jm2 = i+2 + (j+2)*NX;
      } else if (j == 1) {
        indjm1    = i   + (j-1)*NX;
        indjm2    = i   + (j+1)*NX;
        indim1jm1 = i-1 + (j-1)*NX;
        indim2jm2 = i-2 + (j+1)*NX;
        indip1jm1 = i+1 + (j-1)*NX;
        indip2jm2 = i+2 + (j+1)*NX;
      } else if (j == Ntheta - 2) {
        indjp1    = i   + (j+1)*NX;
        indjp2    = i   + (j-1)*NX;
        indim1jp1 = i-1 + (j+1)*NX;
        indim2jp2 = i-2 + (j-1)*NX;
        indip1jp1 = i+1 + (j+1)*NX;
        indip2jp2 = i+2 + (j-1)*NX;
      } else if (j == Ntheta - 1) {
        indjp1    = i   + (j-1)*NX;
        indjp2    = i   + (j-2)*NX;
        indim1jp1 = i-1 + (j-1)*NX;
        indim2jp2 = i-2 + (j-2)*NX;
        indip1jp1 = i+1 + (j-1)*NX;
        indip2jp2 = i+2 + (j-2)*NX;
      }

      const CCTK_REAL lX = X[i];
      /* const CCTK_REAL lth = theta[j]; */

      /* Xtmp[idx]; */
      /* thtmp[idx]; */
      /* printf("X[%3d] = %lf\n", i, X[i]); */


      // 4th order accurate stencils
      const CCTK_REAL W_th = (-W_in[indjp2] + 8 * W_in[indjp1] - 8 * W_in[indjm1] + W_in[indjm2]) *
        oodth12;
      const CCTK_REAL W_thth =  (  -W_in[indjp2] + 16 * W_in[indjp1] - 30 * W_in[ind]
                             + 16 * W_in[indjm1] -      W_in[indjm2] ) * oodthsq12;

      CCTK_REAL W_X, W_Xth;

      if (i == 1 || i == NX - 2) {
        // 1st derivative with 2nd order accuracy (central stencils)
        W_X = (-W_in[indim1] + W_in[indip1]) * 0.5 * oodX;

        W_Xth  = ( W_in[indip1jp1] - W_in[indip1jm1] - W_in[indim1jp1] + W_in[indim1jm1] ) * oodXdth4;

      } else if (i == 0) {
        /* this point is X == 0, r == rH, R == rH/4. dW_dX goes to zero here. but
           since we're interested in dW_dr, and since drxdr diverges (here), we
           will use L'Hopital's rule. for that, we will write instead the 2nd
           derivative */

        // 2nd derivative with 2nd order accuracy (forward stencils)
        W_X = (2*W_in[ind] - 5*W_in[indip1] + 4*W_in[indip2] - W_in[indip3]) * oodXsq;

        // mixed (1st) derivatives with 2nd order accuracy (central stencils in j and forward in i)
        W_Xth = ( - 3*W_in[indjp1] + 3*W_in[indjm1] + 4*W_in[indip1jp1] - 4*W_in[indip1jm1]
                  - W_in[indip2jp1] + W_in[indip2jm1] ) * oodXdth4;

      } else if (i == NX - 1) {
        /* last radial point */

        // 1st derivative with 2nd order accuracy (backward stencils)
        W_X = (W_in[indim2] - 4*W_in[indim1] + 3*W_in[ind]) * 0.5 * oodX;
        W_Xth = 0.; // we don't actually use this variable at large r, so just
                    // set it to zero

      } else {
        // 4th order accurate stencils
        W_X    = (-W_in[indip2] + 8 * W_in[indip1] - 8 * W_in[indim1] + W_in[indim2]) * oodX12;
        W_Xth  = (
            -W_in[indim2jp2] +  8*W_in[indim1jp2] -  8*W_in[indip1jp2] +   W_in[indip2jp2]
         + 8*W_in[indim2jp1] - 64*W_in[indim1jp1] + 64*W_in[indip1jp1] - 8*W_in[indip2jp1]
         - 8*W_in[indim2jm1] + 64*W_in[indim1jm1] - 64*W_in[indip1jm1] + 8*W_in[indip2jm1]
         +   W_in[indim2jm2] -  8*W_in[indim1jm2] +  8*W_in[indip1jm2] -   W_in[indip2jm2] ) * oodXdth144;
      }

      // from the X coordinate used in the input files to the x coordinate
      const CCTK_REAL rx = C0*lX/(1. - lX);
      // from the x coordinate to the metric coordinate r
      // const CCTK_REAL rr = sqrt(rH*rH + rx*rx);

      // corresponding derivatives
      const CCTK_REAL dXdrx = 1./(C0 + rx) - rx/((C0 + rx)*(C0 + rx));

      CCTK_REAL drxdr;
      if (i == 0) { // rx == 0 (X == 0)
        drxdr = sqrt(rH*rH + rx*rx);
      } else {
        drxdr = sqrt(rH*rH + rx*rx)/rx;
      }
      const CCTK_REAL dXdr = dXdrx * drxdr;

      dW_dr_in[ind]    = dXdr * W_X;
      d2W_drth_in[ind] = dXdr * W_Xth;

      dW_dth_in[ind]   = W_th;
      d2W_dth2_in[ind] = W_thth;
    }
  }


  /* now we need to interpolate onto the actual grid points. first let's store
     the grid points themselves in the coordinates (X, theta). */
  const CCTK_INT N_interp_points = cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; // total points

  CCTK_REAL *X_g, *theta_g;
  X_g     = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  theta_g = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));

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
        // TODO: manter?
        if (RR < SMALL)
          RR = SMALL;

        // from (quasi-)isotropic coordinate R to the metric coordinate r
        const CCTK_REAL rr = RR * (1. + 0.25 * rH / RR) * (1. + 0.25 * rH / RR);

        // from the metric coordinate r to the x coordinate
        const CCTK_REAL rx = sqrt(rr*rr - rH*rH);

        // and finally to the X radial coordinate (used in input files)
        const CCTK_REAL lX = rx / (C0 + rx);

        CCTK_REAL ltheta = acos( z1/RR );
        if (ltheta > 0.5*M_PI)    // symmetry along the equatorial plane
          ltheta = M_PI - ltheta;

        X_g[ind]     = lX;
        theta_g[ind] = ltheta;
      }
    }
  }

  /* now for the interpolation */

  // TODO: add these
  CCTK_REAL *dW_dr, *dW_dth, *d2W_dth2, *d2W_drth;
  dW_dr    = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  dW_dth   = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  d2W_dth2 = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  d2W_drth = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));


  const CCTK_INT N_dims  = 2;   // 2-D interpolation

  const CCTK_INT N_input_arrays  = 5;
  const CCTK_INT N_output_arrays = 5;

  /* origin and stride of the input coordinates. with this Cactus reconstructs
     the whole X and theta array. */
  CCTK_REAL origin[N_dims];
  CCTK_REAL delta [N_dims];
  origin[0] = X[0];  origin[1] = theta[0];
  delta[0]  = dX;    delta[1]  = dtheta;

  /* points onto which we want to interpolate, ie, the grid points themselves in
     (X, theta) coordinates (computed above) */
  const void *interp_coords[N_dims];
  interp_coords[0] = (const void *) X_g;
  interp_coords[1] = (const void *) theta_g;


  /* input arrays */
  const void *input_arrays[N_input_arrays];
  CCTK_INT input_array_type_codes[N_input_arrays];
  CCTK_INT input_array_dims[N_dims];
  input_array_dims[0] = NX;
  input_array_dims[1] = Ntheta;

  input_array_type_codes[0] = CCTK_VARIABLE_REAL;
  input_array_type_codes[1] = CCTK_VARIABLE_REAL;
  input_array_type_codes[2] = CCTK_VARIABLE_REAL;
  input_array_type_codes[3] = CCTK_VARIABLE_REAL;
  input_array_type_codes[4] = CCTK_VARIABLE_REAL;

  /* Cactus stores and expects arrays in Fortran order, that is, faster in the
     first index. this is compatible with our input file, where the X coordinate
     is faster. */
  input_arrays[0] = (const void *) F1_in;
  input_arrays[1] = (const void *) F2_in;
  input_arrays[2] = (const void *) F0_in;
  input_arrays[3] = (const void *) phi0_in;
  input_arrays[4] = (const void *) W_in;


  /* output arrays */
  void *output_arrays[N_output_arrays];
  CCTK_INT output_array_type_codes[N_output_arrays];
  CCTK_REAL *F1, *F2, *F0, *phi0, *W;

  F1   = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  F2   = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  F0   = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  phi0 = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
  W    = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));

  output_array_type_codes[0] = CCTK_VARIABLE_REAL;
  output_array_type_codes[1] = CCTK_VARIABLE_REAL;
  output_array_type_codes[2] = CCTK_VARIABLE_REAL;
  output_array_type_codes[3] = CCTK_VARIABLE_REAL;
  output_array_type_codes[4] = CCTK_VARIABLE_REAL;

  output_arrays[0] = (void *) F1;
  output_arrays[1] = (void *) F2;
  output_arrays[2] = (void *) F0;
  output_arrays[3] = (void *) phi0;
  output_arrays[4] = (void *) W;


  /* handle and settings for the interpolation routine */
  int operator_handle, param_table_handle;
  operator_handle    = CCTK_InterpHandle("Lagrange polynomial interpolation");
  param_table_handle = Util_TableCreateFromString("order=4 boundary_extrapolation_tolerance={0.1 1.0 0.05 0.05}");

  CCTK_INFO("Interpolating result...");

  /* do the actual interpolation, and check for error returns */
  int status = CCTK_InterpLocalUniform(N_dims, operator_handle,
                                       param_table_handle,
                                       origin, delta,
                                       N_interp_points,
                                       CCTK_VARIABLE_REAL,
                                       interp_coords,
                                       N_input_arrays, input_array_dims,
                                       input_array_type_codes,
                                       input_arrays,
                                       N_output_arrays, output_array_type_codes,
                                       output_arrays);
  if (status < 0) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
    "interpolation screwed up!");
  }

  free(X_g); free(theta_g);
  free(Xtmp); free(thtmp);
  free(F1_in); free(F2_in); free(F0_in); free(phi0_in); free(W_in);

  free(dW_dr_in); free(dW_dth_in); free(d2W_dth2_in); free(d2W_drth_in);


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
        // TODO: manter?
        if (RR2 < pow(SMALL, 2))
          RR2 = pow(SMALL, 2);
        const CCTK_REAL RR  = sqrt(RR2);

        const CCTK_REAL rho2 = x1*x1 + y1*y1;
        const CCTK_REAL rho  = sqrt(rho2);

        const CCTK_REAL costh  = z1/RR;
        const CCTK_REAL costh2 = costh*costh;
        const CCTK_REAL sinth2 = 1. - costh2;
        const CCTK_REAL sinth  = sqrt(sinth2);

        /*
        const CCTK_REAL R_x = x1/RR;   // dR/dx
        const CCTK_REAL R_y = y1/RR;   // dR/dy
        const CCTK_REAL R_z = z1/RR;   // dR/dz
        */

        const CCTK_REAL sinth2ph_x = -y1/RR2; // sin(th)^2 dphi/dx
        const CCTK_REAL sinth2ph_y =  x1/RR2; // sin(th)^2 dphi/dy

        const CCTK_REAL R2sinth2ph_x = -y1;  // R^2 sin(th)^2 dphi/dx
        const CCTK_REAL R2sinth2ph_y =  x1;  // R^2 sin(th)^2 dphi/dy

        const CCTK_REAL Rsinthth_x  = z1*x1/RR2; // R sin(th) dth/dx
        const CCTK_REAL Rsinthth_y  = z1*y1/RR2; // R sin(th) dth/dy
        const CCTK_REAL Rsinthth_z  = -sinth2;   // R sin(th) dth/dz

        const CCTK_REAL ph = atan2(y1, x1);

        const CCTK_REAL cosph  = cos(ph);
        const CCTK_REAL sinph  = sin(ph);

        const CCTK_REAL cosmph = cos(mm*ph);
        const CCTK_REAL sinmph = sin(mm*ph);

        /* note the division by RR in the following. divisions by zero should be
           avoided by choosing a non-zero value for z0 (for instance) */
        const CCTK_REAL aux  = 1. + 0.25 * rH/RR;
        const CCTK_REAL aux2 = aux  * aux;
        const CCTK_REAL aux4 = aux2 * aux2;
        const CCTK_REAL aux5 = aux4 * aux;
        const CCTK_REAL aux6 = aux4 * aux2;
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

        /*
          KRph/(R sin(th)^2)  = - 1/2 exp(2F2-F0) (1 + rH/(4R))^6 R dW/dr
          Kthph/(R sin(th))^3 = - 1/2 exp(2F2-F0) (1 + rH/(4R))^5 dW/dth / sin(th) 1/(R - rH/4)
        */

        // KRph/(R sin(th)^2)
        const CCTK_REAL KRph_o_Rsinth2 = -0.5 * exp(2. * F2[ind] - F0[ind]) * aux6 * RR * dW_dr[ind];

        const CCTK_REAL den = RR - 0.25 * rH;

        // dW/dth / sin(th) 1/(R - rH/4)
        CCTK_REAL dWdth_o_sinth_den;

        // if at the rho = 0 axis we need to regularize the division by sin(th)
        if (rho < sqrt(dxsq + dysq) * 0.25)
          dWdth_o_sinth_den = d2W_dth2[ind] / den;
        // if at R ~ rH/4 we need to regularize the division by R - rH/4
        else if ( fabs(den) < sqrt(dxsq + dysq + dzsq) * 0.125 )
          dWdth_o_sinth_den = (1. - 0.25*0.25 * rH*rH / RR2) * d2W_drth[ind] / sinth;
        else
          dWdth_o_sinth_den = dW_dth[ind] / (den * sinth);

        // Kthph/(R sin(th))^3
        const CCTK_REAL Kthph_o_R3sinth3 = -0.5 * exp(2. * F2[ind] - F0[ind]) * aux5 * dWdth_o_sinth_den;


        // extrinsic curvature
        kxx[ind] = 2.*KRph_o_Rsinth2 *  x1 * sinth2ph_x                     +  2.*Kthph_o_R3sinth3 *  Rsinthth_x * R2sinth2ph_x;
        kxy[ind] =    KRph_o_Rsinth2 * (x1 * sinth2ph_y + y1 * sinth2ph_x)  +     Kthph_o_R3sinth3 * (Rsinthth_x * R2sinth2ph_y + Rsinthth_y * R2sinth2ph_x);
        kxz[ind] =    KRph_o_Rsinth2 *                    z1 * sinth2ph_x   +     Kthph_o_R3sinth3 *                              Rsinthth_z * R2sinth2ph_x;
        kyy[ind] = 2.*KRph_o_Rsinth2 *  y1 * sinth2ph_y                     +  2.*Kthph_o_R3sinth3 *  Rsinthth_y * R2sinth2ph_y;
        kyz[ind] =    KRph_o_Rsinth2 *                    z1 * sinth2ph_y   +     Kthph_o_R3sinth3 *                              Rsinthth_z * R2sinth2ph_y;
        kzz[ind] = 0.;


        CCTK_REAL alph = exp(F0[ind]) * (RR - 0.25*rH) / (RR + 0.25*rH);
        if (alph < SMALL)
          alph = SMALL;

        // lapse
        if (CCTK_EQUALS(initial_lapse, "psi^n"))
          alp[ind] = pow(psi1, initial_lapse_psi_exponent);
        else if (CCTK_EQUALS(initial_lapse, "HairyBH"))
          alp[ind] = alph;

        // shift
        if (CCTK_EQUALS(initial_shift, "HairyBH")) {
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

  free(F1); free(F2); free(F0); free(phi0); free(W);

  free(dW_dr); free(dW_dth); free(d2W_dth2); free(d2W_drth);

  return;
}
