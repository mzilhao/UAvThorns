#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cctk_Functions.h>

/*
 * Contains helper functions for the TargetTracker thorn.
 */

////////////////////////////////////////////////////////////////////////////////
// Helper struc to pass arguments to functions
////////////////////////////////////////////////////////////////////////////////
struct TargetInfoBundleOneDim {
    CCTK_INT itarget;
    CCTK_INT* ptr_current_id;
    const char* tgt_name;
    const char* dim_name;
};

/////////////////////////////////////////////////////////////////////////////////////////
// Helper function to print an error message and trigger termination after the time step.
// Deactivates the concerned target.
/////////////////////////////////////////////////////////////////////////////////////////
static inline void TerminateHelper(CCTK_ARGUMENTS, const char* message, CCTK_INT itarget) {
    DECLARE_CCTK_ARGUMENTS;

    CCTK_VWARN(1, "%s", message);
    CCTK_VWARN(1, "Deactivating target %d and triggering termination at iteration %d (simulation time %g).", itarget, cctk_iteration, cctk_time);
    CCTK_TerminateNext(cctkGH);

    is_active[itarget] = 0;
    return;
}

////////////////////////////////////////////////////////////////////////////////
// Helper function to update is_loc_from_surface flag.
// We also check validity of the surface index.
////////////////////////////////////////////////////////////////////////////////

static inline void UpdateIsLocFromSurface(CCTK_ARGUMENTS, CCTK_INT itarget) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // Update value of internal is_loc_from_surface flag based on loc_from_surface_mode parameter (which is STEERABLE=ALWAYS)
    is_loc_from_surface[itarget] = loc_from_surface_mode[itarget];

    if (is_loc_from_surface[itarget]) {
        // Not sure the second can happen in here because the relevant parameters
        // are not always steerable, so they pass through ParamCheck.
        // The first one can happen when steering loc_to_surface_mode from no to yes.
        // NOTE: This will be issued and terminate even if the target is not tracked nor active
        if (which_surface_to_store_info[itarget] == -1 || which_surface_to_store_info[itarget] >= nsurfaces) {
            char error_message [1000];
            sprintf(error_message,  "For target %d, 'loc_from_surface_mode' is yes "
                                    "but 'which_surface_to_store_info = %d' is invalid.",
                                    itarget, which_surface_to_store_info[itarget]);
        
            TerminateHelper(CCTK_PASS_CTOC, error_message, itarget);

        }
    }
    return;
}

////////////////////////////////////////////////////////////////////////////////
// Helper function to update is_adjusted flag.
////////////////////////////////////////////////////////////////////////////////

static inline void UpdateAdjustment(CCTK_ARGUMENTS, CCTK_INT itarget) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    is_adjusted[itarget] = tracker_adjust[itarget];
    if (is_adjusted[itarget]) {
        adj_fac_x[itarget] = tracker_x_adjust_factor[itarget];
        adj_fac_y[itarget] = tracker_y_adjust_factor[itarget];
        adj_fac_z[itarget] = tracker_z_adjust_factor[itarget];
        adj_ori_x[itarget] = tracker_x_adjust_origin[itarget];
        adj_ori_y[itarget] = tracker_y_adjust_origin[itarget];
        adj_ori_z[itarget] = tracker_z_adjust_origin[itarget];
    }
    return;
}


////////////////////////////////////////////////////////////////////////////////
// Helper function for factorization and consistency.
// Gives the condition for a target to be active based on the current parameters
// and updates the is_active status of the target accordingly.
////////////////////////////////////////////////////////////////////////////////

static inline void TargetActivationCondition(CCTK_ARGUMENTS, CCTK_INT itarget) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // Update value of internal is_tracked flag based on track parameter (which is STEERABLE=ALWAYS)
    is_tracked[itarget] = track[itarget];

    // Update value of is_loc_from_surface flag based on loc_from_surface_mode parameter (which is STEERABLE=ALWAYS)
    // WARNING: We consider that the target is inactive in that mode, but the tracker still has to follow the surface!
    UpdateIsLocFromSurface(CCTK_PASS_CTOC, itarget);

    // Update value of is_adjusted flag based on tracker_adjust parameter (which is STEERABLE=ALWAYS)
    UpdateAdjustment(CCTK_PASS_CTOC, itarget);

    is_active[itarget] = ( 
        is_tracked[itarget]
        && !loc_from_surface_mode[itarget]
        // &&  track_every[itarget] > 0     // Should be satisfied by construction
        &&  (cctk_time >= start_tracking_after_time[itarget])
        &&  (cctk_time <= stop_tracking_after_time[itarget])
    );
    return;
}