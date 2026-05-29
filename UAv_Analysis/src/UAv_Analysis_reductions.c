#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void UAv_Analysis_IntegrateVol(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (do_analysis_every <= 0) {
    return;
  }

  if (cctk_iteration % do_analysis_every != 0) {
    return;
  }

  enum { num_fields = 17, num_out_vals = 1 };

  // Helper to have a homogeneous interface for reduction output to cactus variables
  CCTK_REAL out_vals[num_fields * num_out_vals];

  CCTK_INT reduction_handle = CCTK_ReductionHandle("sum");
  if (reduction_handle < 0) {
    CCTK_WARN(0, "Could not obtain a handle for sum reduction");
    return;
  }

  // All variables
  const char *const varnames[num_fields] = {
      // 0
      "UAv_Analysis::dE_gf_volume",
      // 1
      "UAv_Analysis::dJx_gf_volume",
      "UAv_Analysis::dJy_gf_volume",
      "UAv_Analysis::dJz_gf_volume",
      // 4
      "UAv_Analysis::dIxx_gf_volume",
      "UAv_Analysis::dIxy_gf_volume",
      "UAv_Analysis::dIxz_gf_volume",
      "UAv_Analysis::dIyy_gf_volume",
      "UAv_Analysis::dIyz_gf_volume",
      "UAv_Analysis::dIzz_gf_volume",
      // 10
      "UAv_Analysis::drho_gf_volume",
      "UAv_Analysis::dCoMx_gf_volume",
      "UAv_Analysis::dCoMy_gf_volume",
      "UAv_Analysis::dCoMz_gf_volume",
      // 14
      "UAv_Analysis::dpx_gf_volume",
      "UAv_Analysis::dpy_gf_volume",
      "UAv_Analysis::dpz_gf_volume"};


  // Get IDs of all variables
  CCTK_INT varid[num_fields];
  for (int i = 0; i < num_fields; ++i) {
    varid[i] = CCTK_VarIndex(varnames[i]);
    if (varid[i] < 0) {
      CCTK_VWARN(0, "Could not get index to grid array %s", varnames[i]);
    }
  }

  // REDUCTION
  if (*early_CoM == 0) {
    // If no early CoM, we reduce all GFs
    CCTK_INT ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, num_out_vals,
        CCTK_VARIABLE_REAL, out_vals, num_fields,
        varid[0],  // E
        varid[1],  varid[2],  varid[3],  // J_i
        varid[4],  varid[5],  varid[6],  varid[7],  varid[8],  varid[9],  // I_ij
        varid[10], varid[11], varid[12], varid[13], // rho and CoM^i
        varid[14], varid[15], varid[16]);// p_i
  
    if (ierr < 0) {
      CCTK_WARN(0, "Error while reducing the auxiliary XX_gf_volume grid functions.");
    }
  } else {
    // If early CoM, we don't recompute it
    const CCTK_INT num_reduced_fields = num_fields - 4; // we reduce 4 fewer variables (rho and CoM^i)
    
    // We use the reduced_vals array with the correct size, then affect output_vals properly
    CCTK_REAL reduced_vals[num_reduced_fields * num_out_vals];

    CCTK_INT ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, num_out_vals,
        CCTK_VARIABLE_REAL, reduced_vals, num_reduced_fields,
        varid[0],  // E
        varid[1],  varid[2],  varid[3],  // J_i
        varid[4],  varid[5],  varid[6],  varid[7],  varid[8],  varid[9],  // I_ij
// SKIP varid[10], varid[11], varid[12], varid[13], // rho and CoM^i
        varid[14], varid[15], varid[16]);// p_i

    if (ierr < 0) {
      CCTK_WARN(0, "Error while reducing the auxiliary XX_gf_volume grid functions.");
    }

    // E
    out_vals[0]  = reduced_vals[0];
    // J_i
    out_vals[1]  = reduced_vals[1];
    out_vals[2]  = reduced_vals[2];
    out_vals[3]  = reduced_vals[3];
    // I_ij
    out_vals[4]  = reduced_vals[4];
    out_vals[5]  = reduced_vals[5];
    out_vals[6]  = reduced_vals[6];
    out_vals[7]  = reduced_vals[7];
    out_vals[8]  = reduced_vals[8];
    out_vals[9]  = reduced_vals[9];
    // p_i
    out_vals[14] = reduced_vals[10];
    out_vals[15] = reduced_vals[11];
    out_vals[16] = reduced_vals[12];
  } // end if early CoM for reduction

  
  // Affect output values to Cactus variables  
  *total_energy = out_vals[0];
  
  *total_angular_momentum_x = out_vals[1];
  *total_angular_momentum_y = out_vals[2];
  *total_angular_momentum_z = out_vals[3];
  
  *Ixx = out_vals[4];
  *Ixy = out_vals[5];
  *Ixz = out_vals[6];
  *Iyy = out_vals[7];
  *Iyz = out_vals[8];
  *Izz = out_vals[9];

  if (*early_CoM == 0) {
    *center_of_mass_x = out_vals[11] / out_vals[10];
    *center_of_mass_y = out_vals[12] / out_vals[10];
    *center_of_mass_z = out_vals[13] / out_vals[10];
  }

  *linear_momentum_x = out_vals[14];
  *linear_momentum_y = out_vals[15];
  *linear_momentum_z = out_vals[16];

  /*
    Note that the integrated values just obtained do *not*, by default,
    take into account any grid symmetries that may be present.
    Use parameter `symmetry` for this. The desired symmetry factors may
    need to be implemented in UAv_Analysis_InitSymmetryFactors. In that
    case, they will be applied to the integrands directly.
    Having the results without symmetry factors applied can be useful,
    but the user should be aware of this choice.
  */
}

///////////////////////////////////////////////////////
// Special function when we need to compute the center of mass earlier than tracking.
// We only compute the GFs needed for the center of mass.
///////////////////////////////////////////////////////


void UAv_Analysis_early_CoM_reduce(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (*early_CoM == 0) {
    return;
  }

  if (do_analysis_every <= 0) {
    return;
  }

  if (cctk_iteration % do_analysis_every != 0) {
    return;
  }

  enum { num_reduced_fields = 4, num_out_vals = 1 };

  // Output array for reduction results
  CCTK_REAL out_vals[num_reduced_fields * num_out_vals];

  CCTK_INT reduction_handle = CCTK_ReductionHandle("sum");
  if (reduction_handle < 0) {
    CCTK_WARN(0, "Could not obtain a handle for sum reduction");
    return;
  }

  const char *const varnames[num_reduced_fields] = {
      "UAv_Analysis::drho_gf_volume",
      "UAv_Analysis::dCoMx_gf_volume",
      "UAv_Analysis::dCoMy_gf_volume",
      "UAv_Analysis::dCoMz_gf_volume"};

  
  // Get IDs of all variables
  CCTK_INT varid[num_reduced_fields];
  for (int i = 0; i < num_reduced_fields; ++i) {
    varid[i] = CCTK_VarIndex(varnames[i]);
    if (varid[i] < 0) {
      CCTK_VWARN(0, "Could not get index to grid array %s", varnames[i]);
    }
  }

  CCTK_INT ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, num_out_vals,
      CCTK_VARIABLE_REAL, out_vals, num_reduced_fields,
      varid[0], varid[1], varid[2], varid[3]);

  if (ierr < 0) {
    CCTK_WARN(0, "Error while reducing the auxiliary XX_gf_volume grid functions.");
    return;
  }

  *center_of_mass_x = out_vals[1] / out_vals[0];
  *center_of_mass_y = out_vals[2] / out_vals[0];
  *center_of_mass_z = out_vals[3] / out_vals[0];
}