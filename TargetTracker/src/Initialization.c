#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include "TargetTracker.h"

void TargetTracker_Initialization(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Initialize targets
  for (CCTK_INT itarget = 0; itarget < nmax_targets; itarget++) {
    // Get the index of variables. It's not supposed to change during the simulation (I think).
    // Validity of parameters should have been checked in ParamCheck.
    // Since validity has been checked, we can use the negative return value for "fixed" targets.
    // WARNING: Since the parameters are steerable, make sure the rest is robust
    
    
    // x source
    target_id_x[itarget] = CCTK_VarIndex (target_x[itarget]);
    // y source
    target_id_y[itarget] = CCTK_VarIndex (target_y[itarget]);
    // z source
    target_id_z[itarget] = CCTK_VarIndex (target_z[itarget]);

    // initial position
    target_loc_x[itarget] = initial_x[itarget];
    target_loc_y[itarget] = initial_y[itarget];
    target_loc_z[itarget] = initial_z[itarget];

    // Set initial active status of the target based on the current parameters.
    TargetActivationCondition(CCTK_PASS_CTOC, itarget);
    
    if (track[itarget]) {
      
      CCTK_VINFO("Initialized %s target %d with sources:", is_active[itarget] ? "active" : "inactive", itarget);
      
      // Check if fixed target in each dimension

      if (target_id_x[itarget] >= 0) {
        CCTK_VINFO("x: %s", target_x[itarget]);
      } else {
        CCTK_VINFO("x = %g (fixed)", initial_x[itarget]);
      }

      if (target_id_y[itarget] >= 0) {
        CCTK_VINFO("y: %s", target_y[itarget]);
      } else {
        CCTK_VINFO("y = %g (fixed)", initial_y[itarget]);
      }

      if (target_id_z[itarget] >= 0) {
        CCTK_VINFO("z: %s", target_z[itarget]);
      } else {
        CCTK_VINFO("z = %g (fixed)", initial_z[itarget]);
      }

    } // end if track

  } // for itarget
}
