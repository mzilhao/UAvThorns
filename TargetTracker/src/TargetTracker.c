#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cctk_Functions.h>

///////////////////////////////////////////////////////////////////////////////
// Declare functions to avoid warnings
////////////////////////////////////////////////////////////////////////////////

// Main function of the TargetTracker thorn
void TargetTracker_SetSurfaces(CCTK_ARGUMENTS);

// Function to update the is_active status of the target based on the current parameters.
CCTK_INT UpdateTargetStatus(CCTK_ARGUMENTS, CCTK_INT itarget);

// Wrapper function to change target for one dimension (x, y, z) if it changed
void TargetChangeOneDim (CCTK_ARGUMENTS, CCTK_INT itarget, 
                            CCTK_INT *const ptr_current_id, const char* tar_name, 
                            const char* dim_name);

                            // Helper function to print error message and terminate. Deactivates the concerned target.
void TerminateHelper(CCTK_ARGUMENTS, const char* message, CCTK_INT itarget);


///////////////////////////////////////////////////////////////////////////////
// This function updates the is_active status of the target based on the current parameters.
// It performs checks to see if the target has changed, and if the new one is valid (since it's always steerable).
///////////////////////////////////////////////////////////////////////////////
CCTK_INT UpdateTargetStatus(CCTK_ARGUMENTS, CCTK_INT itarget) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    // Check if the target is tracked
    is_active[itarget] = ( 
        track[itarget]
        // &&  track_every[itarget] > 0     // Should be satisfied by construction
        &&  (cctk_time >= start_tracking_after_time[itarget])
        &&  (cctk_time <= stop_tracking_after_time[itarget])
    );
    
    if (is_active[itarget]) {
        // Check if target has changed and is valid
        TargetChangeOneDim(CCTK_PASS_CTOC, itarget, &target_id_x[itarget], target_x[itarget], "x");
        TargetChangeOneDim(CCTK_PASS_CTOC, itarget, &target_id_y[itarget], target_y[itarget], "y");
        TargetChangeOneDim(CCTK_PASS_CTOC, itarget, &target_id_z[itarget], target_z[itarget], "z");
    }
    
    return is_active[itarget];
}


///////////////////////////////////////////////////////////////////////////////
// Wrapper function to change target for one dimension (x, y, z) if it changed
///////////////////////////////////////////////////////////////////////////////
void TargetChangeOneDim (CCTK_ARGUMENTS, CCTK_INT itarget, 
                            CCTK_INT *const ptr_current_id, const char* tar_name, 
                            const char* dim_name) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // "New" target
    CCTK_INT tar_id = CCTK_VarIndex(tar_name);
    if (tar_id != *ptr_current_id) { // Target has changed
    
        // Check if new target exists
        if (tar_id < 0) { // No such variable
            char error_message [1000]; 
            sprintf(error_message, "Error while getting ID for target %s variable %s", dim_name, tar_name);
            TerminateHelper(CCTK_PASS_CTOC, error_message, itarget);
            return;
        }
        
        // Check if new target is a grid scalar
        if (CCTK_GroupTypeFromVarI(tar_id) != CCTK_SCALAR) {
            char error_message [1000];
            sprintf(error_message, "Target %s variable %s is not a grid scalar.", dim_name, tar_name);
            TerminateHelper(CCTK_PASS_CTOC, error_message, itarget);
            return;
        }
        
        // Target is valid: update target_id
        if (verbose) {
            CCTK_VINFO("At iteration %d (simulation time %g), target %d %s component was successfully changed from '%s' to '%s'", 
                cctk_iteration, cctk_time,
                itarget, dim_name,
                CCTK_VarName(*ptr_current_id),
                tar_name
            );
        }
        *ptr_current_id = tar_id;
    }
}


///////////////////////////////////////////////////////////////////////////////
// This is the main function of the TargetTracker thorn. 
// It loops over all targets and performs tracking for those that are active during the current iteration.
///////////////////////////////////////////////////////////////////////////////
void TargetTracker_SetSurfaces(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    for (CCTK_INT itarget=0 ; itarget < nmax_targets ; itarget++) {

        // Perform tracking for this target
        // WARNING: Bypass update if we're not on a tracking iteration
        // track_every should be > 0 by construction
        if (cctk_iteration % track_every[itarget] == 0 && UpdateTargetStatus(CCTK_PASS_CTOC, itarget)) { // process iteration

            // Acquire pointer to target variable
            CCTK_REAL *target_loc_x_ptr = (CCTK_REAL *) CCTK_VarDataPtrI(cctkGH, 0, target_id_x[itarget]);
            CCTK_REAL *target_loc_y_ptr = (CCTK_REAL *) CCTK_VarDataPtrI(cctkGH, 0, target_id_y[itarget]);
            CCTK_REAL *target_loc_z_ptr = (CCTK_REAL *) CCTK_VarDataPtrI(cctkGH, 0, target_id_z[itarget]);

            CCTK_INT target_err = 0;
            if (target_loc_x_ptr == NULL) {
                char error_message [1000];
                sprintf(error_message, "Error while acquiring pointer to target %d x variable %s", itarget, target_x[itarget]);
                TerminateHelper(CCTK_PASS_CTOC, error_message, itarget);
                ++target_err;
            }
            else {
                *target_loc_x = *target_loc_x_ptr;
            }

            if (target_loc_y_ptr == NULL) {
                char error_message [1000];
                sprintf(error_message, "Error while acquiring pointer to target %d y variable %s", itarget, target_y[itarget]);
                TerminateHelper(CCTK_PASS_CTOC, error_message, itarget);
                ++target_err;
            }
            else {
                *target_loc_y = *target_loc_y_ptr;
            }

            if (target_loc_z_ptr == NULL) {
                char error_message [1000];
                sprintf(error_message, "Error while acquiring pointer to target %d z variable %s", itarget, target_z[itarget]);
                TerminateHelper(CCTK_PASS_CTOC, error_message, itarget);
                ++target_err;
            }
            else {
                *target_loc_z = *target_loc_z_ptr;
            }

            if (target_err > 0) {
                continue;
            }

            // Update target position
            if (which_surface_to_store_info[itarget] != -1) {
                int sn = which_surface_to_store_info[itarget];

                sf_centroid_x[sn] = *target_loc_x;
                sf_centroid_y[sn] = *target_loc_y;
                sf_centroid_z[sn] = *target_loc_z;

                sf_active[sn] = 1;
                sf_valid[sn] = 1;

                if (verbose) {
                    CCTK_VINFO("Setting spherical surface %d centroid from target #%d to (%g,%g,%g)",
                                sn, itarget, 
                                *target_loc_x, *target_loc_y, *target_loc_z);
                }
            }
        } //end if process iteration
    } // end for loop over targets
}


///////////////////////////////////////////////////////////////////////////////
// Helper function to print error message and terminate. Deactivates the concerned target.
///////////////////////////////////////////////////////////////////////////////
inline void TerminateHelper(CCTK_ARGUMENTS, const char* message, CCTK_INT itarget) {
    DECLARE_CCTK_ARGUMENTS;

    CCTK_VWARN(1, "%s", message);
    CCTK_VWARN(1, "Deactivating target %d and triggering termination at iteration %d (simulation time %g).", itarget, cctk_iteration, cctk_time);
    CCTK_TerminateNext(cctkGH);

    is_active[itarget] = 0;
    return;
}