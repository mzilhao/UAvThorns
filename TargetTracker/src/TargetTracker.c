#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cctk_Functions.h>
#include "TargetTracker.h"

///////////////////////////////////////////////////////////////////////////////
// Declare functions to avoid warnings
///////////////////////////////////////////////////////////////////////////////

// Main function of the TargetTracker thorn
void TargetTracker_SetSurfaces(CCTK_ARGUMENTS);

// Function to update the is_active status of the target based on the current parameters.
CCTK_INT UpdateTargetStatus(CCTK_ARGUMENTS, CCTK_INT itarget);

// Helper function to check a change of target for one dimension (x, y, z)
void TargetChangeOneDim (CCTK_ARGUMENTS, struct TargetInfoBundleOneDim bundle);

// Helper function to get data for pointer in one dimension (x, y, z)
CCTK_INT TargetGetDataOneDim (CCTK_ARGUMENTS, const struct TargetInfoBundleOneDim bundle, CCTK_REAL *const ptr_target_loc);


///////////////////////////////////////////////////////////////////////////////
// This function updates the is_active status of the target based on the current parameters.
// It calls the function that checks parameters and updates the various internal flags.
// It performs checks to see if the target has changed, and if the new one is valid (since it's always steerable).
///////////////////////////////////////////////////////////////////////////////
CCTK_INT UpdateTargetStatus(CCTK_ARGUMENTS, CCTK_INT itarget) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    // ----- Register old flag values to check for changes -----

    // Check if the target was tracked previously
    CCTK_INT was_active = is_active[itarget];

    // Check if the tracker was passively following a surface
    CCTK_INT was_loc_from_surface = is_loc_from_surface[itarget];

    // Check if the tracker was tracking the opposite of the target previously
    CCTK_INT was_opposite = is_opposite[itarget];

    // Check if the tracker position was adjusted previously
    CCTK_INT was_adjusted = is_adjusted[itarget];

    // Dynamic activation of adjustment will deactivate tracking the opposite of a target.
    // The interface quantities (including flags) will be updated below in TargetActivationCondition.
    if (was_opposite && !was_adjusted && tracker_adjust[itarget]) {
        // Info
        if (verbose) {
            CCTK_VINFO("At iteration %d (simulation time %g), tracker %d adjustment was dynamically activated, "
                "which deactivates tracking the opposite of the target.", 
                cctk_iteration, cctk_time, itarget);
        }
        // Set track_opposite parameter, will be effective at the next DECLARE_CCTK_PARAMETERS (in TargetActivationCondition)
        char param_name[100];
        sprintf(param_name, "track_opposite[%d]", itarget);
        CCTK_ParameterSet(param_name, "TargetTracker", "no");
    }


    // ----- Check if the target is tracked now -----
    TargetActivationCondition(CCTK_PASS_CTOC, itarget);

    // Info
    if (verbose) {
        if (is_loc_from_surface[itarget] != was_loc_from_surface) {
            CCTK_VINFO("At iteration %d (simulation time %g), tracker %d changed to %s mode.", 
                cctk_iteration, cctk_time, itarget,
                is_loc_from_surface[itarget] ? "'surface to tracker'" : "'tracker to surface'");
        }

        if (is_active[itarget] != was_active) {
            CCTK_VINFO("At iteration %d (simulation time %g), tracker %d was %sactivated.", 
                cctk_iteration, cctk_time, itarget, is_active[itarget] ? "" : "de");
        }

        if (is_adjusted[itarget] != was_adjusted) {
            CCTK_VINFO("At iteration %d (simulation time %g), tracker %d adjustment was %sactivated.", 
                cctk_iteration, cctk_time, itarget, is_adjusted[itarget] ? "" : "de");

            if (is_adjusted[itarget]) {
                CCTK_VINFO("Adjustment parameters for tracker %d: fac_x = %g, ori_x = %g; fac_y = %g, ori_y = %g; fac_z = %g, ori_z = %g.", 
                    itarget, adj_fac_x[itarget], adj_ori_x[itarget], adj_fac_y[itarget], adj_ori_y[itarget], adj_fac_z[itarget], adj_ori_z[itarget]);
            }
        }

        // Put this message after the adjustment message, since dynamically turning on adjustment
        // shuts down tracking the opposite.
        if (is_opposite[itarget] != was_opposite) {
            CCTK_VINFO("At iteration %d (simulation time %g), tracker %d is %s tracking the opposite of the target.", 
                cctk_iteration, cctk_time, itarget, is_opposite[itarget] ? "now" : "no longer");
        }
    } // end if verbose
    

    // ----- Check if target has changed and is valid -----
    if (is_active[itarget]) {
        struct TargetInfoBundleOneDim bundle_x = {itarget, &target_id_x[itarget], target_x[itarget], "x"};
        TargetChangeOneDim(CCTK_PASS_CTOC, bundle_x);
        struct TargetInfoBundleOneDim bundle_y = {itarget, &target_id_y[itarget], target_y[itarget], "y"};
        TargetChangeOneDim(CCTK_PASS_CTOC, bundle_y);
        struct TargetInfoBundleOneDim bundle_z = {itarget, &target_id_z[itarget], target_z[itarget], "z"};
        TargetChangeOneDim(CCTK_PASS_CTOC, bundle_z);
    }
    
    return is_active[itarget];
}


///////////////////////////////////////////////////////////////////////////////
// Helper function to check a change of target for one dimension (x, y, z)
///////////////////////////////////////////////////////////////////////////////
void TargetChangeOneDim (CCTK_ARGUMENTS, struct TargetInfoBundleOneDim bundle) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // Unpack bundle
    const CCTK_INT itarget    = bundle.itarget;
    CCTK_INT* ptr_current_id  = bundle.ptr_current_id;
    const char* tgt_name      = bundle.tgt_name;
    const char* dim_name      = bundle.dim_name;

    // "New" target
    CCTK_INT tgt_id = CCTK_VarIndex(tgt_name);
    if (tgt_id != *ptr_current_id) { // Target has changed
    
        // Check if target is fixed
        if (!CCTK_Equals(tgt_name, "")) {
            
            // Check if new target exists
            if (tgt_id < 0) { // No such variable
                char error_message [1000]; 
                sprintf(error_message, "Could not get index of target %d %s variable %s", itarget, dim_name, tgt_name);
                TerminateHelper(CCTK_PASS_CTOC, error_message, itarget);
                return;
            }
            
            // Check if new target is a grid scalar
            if (CCTK_GroupTypeFromVarI(tgt_id) != CCTK_SCALAR) {
                char error_message [1000];
                sprintf(error_message, "Target %d %s variable %s is not a grid scalar.", itarget, dim_name, tgt_name);
                TerminateHelper(CCTK_PASS_CTOC, error_message, itarget);
                return;
            }
        }
        
        // Target is valid: update target_id
        if (verbose) {
            CCTK_VINFO("At iteration %d (simulation time %g), target %d %s component was successfully changed from '%s' to '%s'", 
                cctk_iteration, cctk_time,
                itarget, dim_name,
                CCTK_VarName(*ptr_current_id),
                tgt_name
            );
        }
        *ptr_current_id = tgt_id;
    }
}


////////////////////////////////////////////////////////////////////////////////
// Helper function to get data for pointer in one dimension (x, y, z)
// Returns 0 if no error, 1 if error (and deactivates target and triggers termination).
////////////////////////////////////////////////////////////////////////////////
CCTK_INT TargetGetDataOneDim (CCTK_ARGUMENTS, const struct TargetInfoBundleOneDim bundle, CCTK_REAL *const ptr_target_loc) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // Unpack bundle
    const CCTK_INT itarget    = bundle.itarget;
    const CCTK_INT tgt_id     =*bundle.ptr_current_id;
    const char* tgt_name      = bundle.tgt_name;
    const char* dim_name      = bundle.dim_name;

    // A negative index here should mean fixed tracker (and not an error):
    // keep current position in that case
    if (tgt_id >= 0) {
        // Acquire pointer to target variable
        CCTK_REAL *ptr_new_loc = (CCTK_REAL *) CCTK_VarDataPtrI(cctkGH, 0, tgt_id);
        
        if (ptr_new_loc == NULL) {
            char error_message [1000];
            sprintf(error_message, "Error while acquiring pointer to target %d %s variable %s", itarget, dim_name, tgt_name);
            TerminateHelper(CCTK_PASS_CTOC, error_message, itarget);
            return 1;
        }
        else {
            *ptr_target_loc = *ptr_new_loc;
        }
    }
    return 0;
}


///////////////////////////////////////////////////////////////////////////////
// This is the main function of the TargetTracker thorn. 
// It loops over all trackers and performs tracking for those that are active during the current iteration.
///////////////////////////////////////////////////////////////////////////////
void TargetTracker_SetSurfaces(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    for (CCTK_INT itarget=0 ; itarget < nmax_targets ; itarget++) {

        // Perform tracking for this target
        // WARNING: Bypass update if we're not on a tracking iteration
        // track_every should be > 0 by construction
        if (cctk_iteration % track_every[itarget] == 0) { // process iteration
            // WARNING: UpdateTargetStatus needs to be called first to update the flags!
            if (UpdateTargetStatus(CCTK_PASS_CTOC, itarget)) { // active tracker in tracker to surface mode
            
                // Normal case: get target location from variables
                // In surface to tracker mode, the tracker is considered inactive

                CCTK_INT target_err = 0;
                // Get target position for each dimension and check for errors
                struct TargetInfoBundleOneDim bundle_x = {itarget, &target_id_x[itarget], target_x[itarget], "x"};
                target_err += TargetGetDataOneDim(CCTK_PASS_CTOC, bundle_x, &target_loc_x[itarget]);
                struct TargetInfoBundleOneDim bundle_y = {itarget, &target_id_y[itarget], target_y[itarget], "y"};
                target_err += TargetGetDataOneDim(CCTK_PASS_CTOC, bundle_y, &target_loc_y[itarget]);
                struct TargetInfoBundleOneDim bundle_z = {itarget, &target_id_z[itarget], target_z[itarget], "z"};
                target_err += TargetGetDataOneDim(CCTK_PASS_CTOC, bundle_z, &target_loc_z[itarget]);

                // An error will trigger termination at the end of the time step.
                if (target_err > 0) {
                    continue;
                }

                // Adjust location if needed (not applied to fixed tracker)
                // (We need the if statement because of dynamic steerability of the flags)
                if (is_opposite[itarget] || is_adjusted[itarget]) {
                    if (*bundle_x.ptr_current_id >= 0) {
                        target_loc_x[itarget] = adj_fac_x[itarget] * target_loc_x[itarget] + adj_ori_x[itarget];
                    }
                    if (*bundle_y.ptr_current_id >= 0) {
                        target_loc_y[itarget] = adj_fac_y[itarget] * target_loc_y[itarget] + adj_ori_y[itarget];
                    }
                    if (*bundle_z.ptr_current_id >= 0) {
                        target_loc_z[itarget] = adj_fac_z[itarget] * target_loc_z[itarget] + adj_ori_z[itarget];
                    }
                } // end if adjustment
                
                // Update spherical surface with tracker position
                if (which_surface_to_store_info[itarget] != -1) {
                    int sn = which_surface_to_store_info[itarget];

                    sf_centroid_x[sn] = target_loc_x[itarget];
                    sf_centroid_y[sn] = target_loc_y[itarget];
                    sf_centroid_z[sn] = target_loc_z[itarget];

                    sf_active[sn] = 1;
                    sf_valid[sn]  = 1;

                    if (verbose) {
                        CCTK_VINFO("Setting spherical surface %d centroid from tracker #%d to (%g,%g,%g)",
                                    sn, itarget, 
                                    target_loc_x[itarget], target_loc_y[itarget], target_loc_z[itarget]);
                    }
                }
            
            } else if (is_tracked[itarget] && is_loc_from_surface[itarget]) { // tracked target in surface to tracker mode
                // If loc_from_surface_mode is on, we don't try to get the target location from the variables,
                // but directly from the surface. This can be useful for continuity when a horizon is formed for instance.
                // The case where which_surface_to_store_info is invalid should be caught before.
                
                // TODO: Check if surface is active / is valid ?
                target_loc_x[itarget] = sf_centroid_x[which_surface_to_store_info[itarget]];
                target_loc_y[itarget] = sf_centroid_y[which_surface_to_store_info[itarget]];
                target_loc_z[itarget] = sf_centroid_z[which_surface_to_store_info[itarget]];

                if (verbose) {
                    CCTK_VINFO("Setting tracker %d location from surface %d to (%g,%g,%g)",
                                itarget, which_surface_to_store_info[itarget], 
                                target_loc_x[itarget], target_loc_y[itarget], target_loc_z[itarget]);
                }
            } // else we don't do anything
        } //end if process iteration
    } // end for loop over targets
}