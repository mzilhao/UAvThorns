#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include "TargetTracker.h"

///////////////////////////////////////////////////////////////////////////////
// First initialization of one target
///////////////////////////////////////////////////////////////////////////////
void InitializeOneTarget (CCTK_ARGUMENTS, CCTK_INT itarget) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // Get the index of variables. It's not supposed to change during the simulation (I think).
    // Validity of parameters should have been checked in ParamCheck.
    // Since validity has been checked, a negative value here should mean fixed tracker (empty string "^$").
    // WARNING: Since the parameters are steerable, make sure the rest is robust
    
    
    // x source
    target_id_x[itarget] = CCTK_VarIndex (target_x[itarget]);
    // y source
    target_id_y[itarget] = CCTK_VarIndex (target_y[itarget]);
    // z source
    target_id_z[itarget] = CCTK_VarIndex (target_z[itarget]);

    // Initialize all adjustment parameters
    is_adjusted[itarget] = tracker_adjust[itarget];
    adj_fac_x[itarget]   = tracker_x_adjust_factor[itarget];
    adj_fac_y[itarget]   = tracker_y_adjust_factor[itarget];
    adj_fac_z[itarget]   = tracker_z_adjust_factor[itarget];
    adj_ori_x[itarget]   = tracker_x_adjust_origin[itarget];
    adj_ori_y[itarget]   = tracker_y_adjust_origin[itarget];
    adj_ori_z[itarget]   = tracker_z_adjust_origin[itarget];

    // Initial position
    // Will be overriden with the correct values at ANALYSIS if needed 
    // WARNING: if fixed target and adjusted target (not so sensible), this will be reapeated at t=0
    target_loc_x[itarget] = adj_fac_x[itarget] * initial_x[itarget] + adj_ori_x[itarget];
    target_loc_y[itarget] = adj_fac_y[itarget] * initial_y[itarget] + adj_ori_y[itarget];
    target_loc_z[itarget] = adj_fac_z[itarget] * initial_z[itarget] + adj_ori_z[itarget];
    
    // Set initial value of is_loc_from_surface flag.
    // ParamCheck prevents a wrong surface index here.
    // TargetActivationCondition will repeat it, 
    // but doing it once here allows to avoid triggering the change message.
    is_loc_from_surface[itarget] = loc_from_surface_mode[itarget];
    
    // Set initial active status of the target based on the current parameters.
    TargetActivationCondition(CCTK_PASS_CTOC, itarget);
    
    if (is_tracked[itarget]) {
        
        if (which_surface_to_store_info[itarget] != -1) {
            char* from_surf_str = loc_from_surface_mode[itarget] ? 
                                "'surface to tracker'" : 
                                "'tracker to surface'" ;
            CCTK_VINFO("Tracker %d associated with surface %d is in %s mode.", 
                        itarget, which_surface_to_store_info[itarget], from_surf_str);
        }
        
        if (is_adjusted[itarget]) {
            CCTK_VINFO("Adjustment is activated for tracker %d with parameters: fac_x = %g, ori_x = %g; fac_y = %g, ori_y = %g; fac_z = %g, ori_z = %g.", 
                itarget, adj_fac_x[itarget], adj_ori_x[itarget], adj_fac_y[itarget], adj_ori_y[itarget], adj_fac_z[itarget], adj_ori_z[itarget]);
        }
        
        CCTK_VINFO("Initialized %s tracker %d with sources:", is_active[itarget] ? "active" : "inactive", itarget);

        // Check if fixed target in each dimension

        if (target_id_x[itarget] >= 0) {
            CCTK_VINFO("x: %s", target_x[itarget]);
        } else {
            CCTK_VINFO("x = %g (fixed)", target_loc_x[itarget]);
        }

        if (target_id_y[itarget] >= 0) {
            CCTK_VINFO("y: %s", target_y[itarget]);
        } else {
            CCTK_VINFO("y = %g (fixed)", target_loc_y[itarget]);
        }

        if (target_id_z[itarget] >= 0) {
            CCTK_VINFO("z: %s", target_z[itarget]);
        } else {
            CCTK_VINFO("z = %g (fixed)", target_loc_z[itarget]);
        }

    } // end if track
}

///////////////////////////////////////////////////////////////////////////////
// Helper function to recover the name of the target variable for one dimension (x, y, z) from its checkpointed ID.
///////////////////////////////////////////////////////////////////////////////
void RecoverOneTargetNameOneDim (CCTK_ARGUMENTS, const struct TargetInfoBundleOneDim bundle) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // Unpack bundle
    const CCTK_INT itarget          = bundle.itarget;
    const CCTK_INT* ptr_current_id  = bundle.ptr_current_id;
    const char* tgt_name            = bundle.tgt_name;
    const char* dim_name            = bundle.dim_name;

    const CCTK_INT recovery_ID      = CCTK_VarIndex (tgt_name);
    char checkpointed_Name[150];
    if (*ptr_current_id >= 0) {
        sprintf(checkpointed_Name, "%s::%s", CCTK_ImpFromVarI(*ptr_current_id), CCTK_VarName(*ptr_current_id));
    }
    else {
        sprintf(checkpointed_Name, "%s", "");
    }

    // Output info
    if (recovery_ID != *ptr_current_id || !CCTK_Equals(tgt_name, checkpointed_Name)) {
        char message[1000];
        sprintf(message, "At recovery of target %d, parameter 'target_%s' is '%s' (ID: %d),\n"
                         "    but the variable with checkpointed index 'target_id_%s[%d]' = %d is '%s'.\n"
                         "    Setting target_%s[%d] = %s.",
                         itarget, dim_name, tgt_name, recovery_ID,
                         dim_name, itarget, *ptr_current_id, checkpointed_Name,
                         dim_name, itarget, checkpointed_Name);
        if (verbose) {
            CCTK_VINFO("%s", message);
        }
        else {
            CCTK_VWARN(CCTK_WARN_COMPLAIN, "%s", message);
        }
    } //end if mismatch

    // Set parameter to checkpointed value
    // A recovered negative index should mean a fixed target
    // The ParameterSet will take effect at the next DECLARE_CCTK_PARAMETERS
    char param_name[100];
    sprintf(param_name, "target_%s[%d]", dim_name, itarget);
    CCTK_ParameterSet(param_name, "TargetTracker", (*ptr_current_id >= 0) ? checkpointed_Name : "");
}

///////////////////////////////////////////////////////////////////////////////
// (Partial) recovery of parameters for one target (after restart from checkpoint)
// This is not completely robust for all parameters and all scenarios.
///////////////////////////////////////////////////////////////////////////////
void RecoverOneTarget (CCTK_ARGUMENTS, CCTK_INT itarget) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    // The ParameterSet will take effect at the next DECLARE_CCTK_PARAMETERS

    char param_name[100];

    /* Parameter: track 
     * Most important to recover since it controls the tracking and is always steerable.
     */
    // Info
    if (track[itarget] != is_tracked[itarget]) {
        char message[1000];
        sprintf(message, "At recovery of target %d, parameter 'track' is '%s', but the internal flag 'is_tracked' is '%s'.\n"
                         "    Setting track[%d] = %s.",
                         itarget, track[itarget] ? "yes" : "no", is_tracked[itarget] ? "yes" : "no",
                         itarget, is_tracked[itarget] ? "yes" : "no");
        if (verbose) {
            CCTK_VINFO("%s", message);
        }
        else {
            CCTK_VWARN(CCTK_WARN_COMPLAIN, "%s", message);
        }
    }
    // Set
    sprintf(param_name, "track[%d]", itarget);
    CCTK_ParameterSet(param_name, "TargetTracker", is_tracked[itarget] ? "yes" : "no");

    //////////////////////////////////////////////
    
    /* Parameter: loc_from_surface_mode
     * Controls the tracking mode and is always steerable.
     */
    // Info
    if (loc_from_surface_mode[itarget] != is_loc_from_surface[itarget]) {
        char message[1000];
        sprintf(message, "At recovery of target %d, parameter 'loc_from_surface_mode' is '%s', but the internal flag 'is_loc_from_surface' is '%s'.\n"
                         "    Setting loc_from_surface_mode[%d] = %s.",
                         itarget, loc_from_surface_mode[itarget] ? "yes" : "no", is_loc_from_surface[itarget] ? "yes" : "no",
                         itarget, is_loc_from_surface[itarget] ? "yes" : "no");
        if (verbose) {
            CCTK_VINFO("%s", message);
        }
        else {
            CCTK_VWARN(CCTK_WARN_COMPLAIN, "%s", message);
        }
    }
    // Set
    sprintf(param_name, "loc_from_surface_mode[%d]", itarget);
    CCTK_ParameterSet(param_name, "TargetTracker", is_loc_from_surface[itarget] ? "yes" : "no");
    
    //////////////////////////////////////////////
    
    /* Parameter: tracker_adjust
     * Controls the adjustment mode and is always steerable.
     */
    // Info
    if (tracker_adjust[itarget] != is_adjusted[itarget]) {
        char message[1000];
        sprintf(message, "At recovery of tracker %d, parameter 'tracker_adjust' is '%s', but the internal flag 'is_adjusted' is '%s'.\n"
                         "    Setting tracker_adjust[%d] = %s.",
                         itarget, tracker_adjust[itarget] ? "yes" : "no", is_adjusted[itarget] ? "yes" : "no",
                         itarget, is_adjusted[itarget] ? "yes" : "no");
        if (verbose) {
            CCTK_VINFO("%s", message);
        }
        else {
            CCTK_VWARN(CCTK_WARN_COMPLAIN, "%s", message);
        }
    }
    // Set
    sprintf(param_name, "tracker_adjust[%d]", itarget);
    CCTK_ParameterSet(param_name, "TargetTracker", is_adjusted[itarget] ? "yes" : "no");

    //////////////////////////////////////////////

    /* Parameters: target variables
     * We recover the names from their checkpointed IDs.
     * NOTE: If these parameters become not always steerable, this step is probably unnecessary.
     */
    const struct TargetInfoBundleOneDim bundle_x = {itarget, &target_id_x[itarget], target_x[itarget], "x"};
    RecoverOneTargetNameOneDim(CCTK_PASS_CTOC, bundle_x);
    const struct TargetInfoBundleOneDim bundle_y = {itarget, &target_id_y[itarget], target_y[itarget], "y"};
    RecoverOneTargetNameOneDim(CCTK_PASS_CTOC, bundle_y);
    const struct TargetInfoBundleOneDim bundle_z = {itarget, &target_id_z[itarget], target_z[itarget], "z"};
    RecoverOneTargetNameOneDim(CCTK_PASS_CTOC, bundle_z);

    //////////////////////////////////////////////

    /* Parameters that we don't try to recover:
     * - force_params_at_recovery: STEERABLE=RECOVER
     * - adjustments parameters (fac and ori): STEERABLE=RECOVER for simplicity for now, taken from the param file as such
     * - track_every: STEERABLE=RECOVER
     * - start/stop_tracking_after_time: STEERABLE=RECOVER
     * - which_surface_to_store_info: STEERABLE=RECOVER
     * - initial_x/y/z: STEERABLE=RECOVER, and they are ignored anyway if force_params_at_recovery = no
     * - verbose: STEERABLE=ALWAYS, but it's not critical (and the user should know if they try to steer it mid-run)
     */



}

///////////////////////////////////////////////////////////////////////////////
// Routine called at the first initialization of the simulation (i.e. not recovery)
///////////////////////////////////////////////////////////////////////////////
void TargetTracker_Initialization(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    for (CCTK_INT itarget = 0; itarget < nmax_targets; itarget++) {
        InitializeOneTarget(CCTK_PASS_CTOC, itarget);
    }
}

///////////////////////////////////////////////////////////////////////////////
// Routine called at recovery
///////////////////////////////////////////////////////////////////////////////
void TargetTracker_Recovery(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    /* Fundamental ambiguity: have parameters changed during the run and the parameter is reused identically,
     * or has the user purposefully changed the parameter file?
     * Having a fully consistent and robust architecture seems involved, so we have to make assumptions and concessions.
     *
     * By default, we want to keep the recovered parameters, but without having to fully duplicate parameters into the interface.
     * (I don't fully understand if we can get the checkpointed values of (steerable) parameters and store them somehow.)
     * The user can force a re-read of parameters at recovery (for those that we try to infer).
     */
    
    for (CCTK_INT itarget = 0; itarget < nmax_targets; itarget++) {
        if (force_params_at_recovery[itarget]) {
            if (verbose) {
                CCTK_VINFO("Forcing re-initialization of target %d from parameter file at recovery.", itarget);
            }
            InitializeOneTarget(CCTK_PASS_CTOC, itarget);
        }
        else {
            RecoverOneTarget(CCTK_PASS_CTOC, itarget);
        }
    }
}
