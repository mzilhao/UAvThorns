#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

void MagScalarBase_Zero(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int const np = cctk_ash[0] * cctk_ash[1] * cctk_ash[2];

#pragma omp parallel for
  for (int i = 0; i < np; ++i) {
    Ex  [i] = 0.0;
    Ey  [i] = 0.0;
    Ez  [i] = 0.0;
    Ax  [i] = 0.0;
    Ay  [i] = 0.0;
    Az  [i] = 0.0;
    Aphi[i] = 0.0;
    Zeta[i] = 0.0;
    phi1[i] = 0.0;
    phi2[i] = 0.0;
    Kphi1[i]= 0.0;
    Kphi2[i]= 0.0;
  }

  if (CCTK_EQUALS(initial_data_setup_method, "init_some_levels") ||
      CCTK_EQUALS(initial_data_setup_method, "init_single_level")) {
    /* do nothing */
  } else if (CCTK_EQUALS(initial_data_setup_method, "init_all_levels")) {
      if (CCTK_ActiveTimeLevels(cctkGH, "MagScalarBase::Zeta") >= 2) {
#pragma omp parallel for
          for (int i = 0; i < np; ++i) {
              Ex_p  [i] = 0.0;
              Ey_p  [i] = 0.0;
              Ez_p  [i] = 0.0;
              Ax_p  [i] = 0.0;
              Ay_p  [i] = 0.0;
              Az_p  [i] = 0.0;
              Aphi_p[i] = 0.0;
              Zeta_p[i] = 0.0;
              phi1_p[i] = 0.0;
              phi2_p[i] = 0.0;
              Kphi1_p[i]= 0.0;
              Kphi2_p[i]= 0.0;
          }
      }

      if (CCTK_ActiveTimeLevels(cctkGH, "MagScalarBase::Zeta") >= 3) {
#pragma omp parallel for
          for (int i = 0; i < np; ++i) {
              Ex_p_p  [i] = 0.0;
              Ey_p_p  [i] = 0.0;
              Ez_p_p  [i] = 0.0;
              Ax_p_p  [i] = 0.0;
              Ay_p_p  [i] = 0.0;
              Az_p_p  [i] = 0.0;
              Aphi_p_p[i] = 0.0;
              Zeta_p_p[i] = 0.0;
              phi1_p_p[i] = 0.0;
              phi2_p_p[i] = 0.0;
              Kphi1_p_p[i]= 0.0;
              Kphi2_p_p[i]= 0.0;
          }
      }
  } else {
    CCTK_WARN(
        CCTK_WARN_ABORT,
        "Unsupported parameter value for InitBase::initial_data_setup_method");
  }
}
