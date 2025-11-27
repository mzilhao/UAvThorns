#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

// Initialize auxiliary members
void UAv_Initialization (CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;


  // -------------------------------
  // MULTIPATCH
  // -------------------------------

  // Trick to check multipatch usage without having to inherit the thorn
  *is_multipatch = CCTK_IsFunctionAliased("MultiPatch_GetDomainSpecification");
  
  // We need to check the volume form in the Initialization function and not in ParamCheck,
  // because the value of the parameter Coordinates::store_volume_form alone is not enough
  *use_volume_form = 0; // default: don't use volume form


  if (*is_multipatch > 0) {
    // If the simulation uses multipatch, make sure that the volume form variable is stored
    // WARNING: Not all patch systems use the volume form... 
    //          We issue the warning, but don't abort, to keep the abstraction layer and avoid splitting cases.
    CCTK_INT* volume_form_state_ptr = CCTK_VarDataPtr(cctkGH, 0, "Coordinates::volume_form_state");
    //CCTK_VINFO("volume_form_state = %d", *volume_form_state_ptr);
    
    if (volume_form_state_ptr != NULL) {
      // -1 is the "error" return value of ParameterGetType() (see repos/flesh/src/main/Parameters.c and repos/flesh/src/include/cctk_Parameter.h)
      CCTK_INT type = -1; 

      const CCTK_INT* store_volume_form_ptr = CCTK_ParameterGet("store_volume_form", "Coordinates", &type);
      if (store_volume_form_ptr == NULL || type != PARAMETER_BOOLEAN) {
        CCTK_ERROR("Problem acquiring pointer to Coordinates::store_volume_form parameter.");
      }
      // Coordinates::store_volume_form = yes
      else if (*store_volume_form_ptr) {
        // No volume form actually stored
        // TODO/WARNING: the variable might be uninitialized in that case? 
        //               It indeed happens that it is polluted and has 1 in memory, which causes the wrong behavior>>>
        if (*volume_form_state_ptr != 1) {
          CCTK_WARN(1, "The patch system that you are using does not store the volume form, although you set Coordinates::store_volume_form = yes. "     
                       "UAv_Analysis thorn will not be able to use the volume form for integrations, so results may be incorrect.");
        }
        // else, we're good to go
        else {
          CCTK_INFO("Using volume form from multipatch system for integrations in UAv_Analysis.");
          *use_volume_form = 1;
        }
      }
      // Coordinates::store_volume_form = no
      else {
        CCTK_WARN(1, "You are using a multipatch system, but you set Coordinates::store_volume_form = no. "   
                     "UAv_Analysis thorn will not be able to use the volume form for integrations, so results may be incorrect." 
                     "If your patch system allows volume form computation/storage, please set Coordinates::store_volume_form = yes.");
      }
    }
    else {
      CCTK_ERROR("Problem acquiring pointer to Coordinates::volume_form_state variable.");
    }
  }


  // -------------------------------
  // ORIGIN TRACKING
  // -------------------------------

  // Initialize origin tracking if needed
  if (track_origin_from_grid_scalar) {
    // Get the index of variables. It's not supposed to change during the simulation (I think).
    // Validity of parameters should have been checked in ParamCheck.
    // Seems a bit redundant to do this affectation here and not in ParamCheck, but more in the logic.
    
    // x source
    *origin_from_grid_scalar_index_x = CCTK_VarIndex (track_origin_source_x);
    // y source
    *origin_from_grid_scalar_index_y = CCTK_VarIndex (track_origin_source_y);
    // z source
    *origin_from_grid_scalar_index_z = CCTK_VarIndex (track_origin_source_z);

    CCTK_VINFO("Tracking origin used in the analysis with grid scalars.");
    CCTK_VINFO("x0 = %s", track_origin_source_x);
    CCTK_VINFO("y0 = %s", track_origin_source_y);
    CCTK_VINFO("z0 = %s", track_origin_source_z);
  }
  else { // no tracking from grid scalar
    // origin_from_grid_scalar_index not allocated in schedule.ccl in that case

    // We can already initialize the coordinates
    *x0 = origin_x;
    *y0 = origin_y;
    *z0 = origin_z;

    CCTK_VINFO("Using fixed origin in the analysis.");
    CCTK_VINFO("x0 = %g", *x0);
    CCTK_VINFO("y0 = %g", *y0);
    CCTK_VINFO("z0 = %g", *z0);
  }
}