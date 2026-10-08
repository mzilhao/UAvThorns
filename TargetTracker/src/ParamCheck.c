#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cctk_Functions.h>

void TargetTracker_ParamCheck(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // -------------------------------
    // CENTER OF MASS CALCULATION
    // -------------------------------

    // See if the calculation could be requested at all
    CCTK_INT want_com = 0;
    for (CCTK_INT itarget = 0; itarget < nmax_targets; itarget++) {
        if (compute_com[itarget]) {
            want_com = 1;
        }
    }

    if (want_com && com_every <= 0) {
        CCTK_VPARAMWARN("'com_every' is <= 0, but some 'compute_com' parameter is set to yes: the center of mass will not be computed.");
    }

    // For multipatch / Llama grids, we need the volume form computed in Coordinates.
    // This is activated by the parameter Coordinates::store_volume_form = yes.
    // (adapted from UAv_Analysis)
    if (want_com && com_every > 0) {
        // Trick to check multipatch usage without having to inherit the thorn
        if( CCTK_IsFunctionAliased("MultiPatch_GetDomainSpecification") ) {

            // -1 is the "error" return value of ParameterGetType() (see the flesh)
            CCTK_INT type = -1;
            const CCTK_INT* store_volume_form_ptr = CCTK_ParameterGet("store_volume_form", "Coordinates", &type);

            if (store_volume_form_ptr == NULL || type != PARAMETER_BOOLEAN) {
                CCTK_ERROR("Problem acquiring pointer to Coordinates::store_volume_form parameter.");
            }
            else if (!(*store_volume_form_ptr)) {
                CCTK_VPARAMWARN("You are using a multipatch system, and the center of mass is computed, but you set Coordinates::store_volume_form = no. "
                                "TargetTracker will not be able to use the volume form for the integration. "
                                "Please set Coordinates::store_volume_form = yes.");
            }
            // else, we're good to go
        }
    }

    // Check the excision parameter (runtime check in the computation:
    // it is fatal there if the mask type is not registered by e.g. AHFinderDirect.)
    // WARNING: CCTK_IsFunctionAliased is not usable here (SpaceMask functions are
    //          plain C symbols, not schedule aliases)

    // -------------------------------
    // TRACKING
    // -------------------------------

    for (CCTK_INT itarget = 0; itarget < nmax_targets; itarget++) {
        
        if (loc_from_surface_mode[itarget] && which_surface_to_store_info[itarget] == -1) {
            CCTK_VPARAMWARN("For target %d, 'loc_from_surface_mode' is yes but no surface is specified in 'which_surface_to_store_info'.", itarget);
        }

        if (which_surface_to_store_info[itarget] != -1 && which_surface_to_store_info[itarget] >= nsurfaces) {
            CCTK_VPARAMWARN("For target %d, surface index is greater than the number of spherical surfaces.", itarget);
        }
        
        // Soft warning for the user if they set both tracker_adjust and track_opposite to yes.
        if (tracker_adjust[itarget] && track_opposite[itarget]) {
            CCTK_VWARN(CCTK_WARN_COMPLAIN, "For target %d, both 'tracker_adjust' and 'track_opposite' are set to yes."
                            "'track_opposite' takes precedence and 'tracker_adjust' will be ignored.", itarget);
        }
        
        // --------------------------------
        // CHECK TARGET VARIABLES
        // --------------------------------
        CCTK_INT index;
        // x
        index = CCTK_VarIndex (target_x[itarget]);
        if (!CCTK_Equals(target_x[itarget], "")) {
            if (index < 0) {
                CCTK_VPARAMWARN("Could not get index of target %d x variable: %s.", itarget, target_x[itarget]);
            }
            else if (CCTK_GroupTypeFromVarI(index) != CCTK_SCALAR) {
                CCTK_VPARAMWARN("Target %d x variable: %s is not a grid scalar.", itarget, target_x[itarget]);
            }
        }
        // y
        index = CCTK_VarIndex (target_y[itarget]);
        if (!CCTK_Equals(target_y[itarget], "")) {
            if (index < 0) {
                CCTK_VPARAMWARN("Could not get index of target %d y variable: %s.", itarget, target_y[itarget]);
            }
            else if (CCTK_GroupTypeFromVarI(index) != CCTK_SCALAR) {
                CCTK_VPARAMWARN("Target %d y variable: %s is not a grid scalar.", itarget, target_y[itarget]);
            }
        }
        // z
        index = CCTK_VarIndex (target_z[itarget]);
        if (!CCTK_Equals(target_z[itarget], "")) {
            if (index < 0) {
                CCTK_VPARAMWARN("Could not get index of target %d z variable: %s.", itarget, target_z[itarget]);
            }
            else if (CCTK_GroupTypeFromVarI(index) != CCTK_SCALAR) {
                CCTK_VPARAMWARN("Target %d z variable: %s is not a grid scalar.", itarget, target_z[itarget]);
            }
        }
    }
}
