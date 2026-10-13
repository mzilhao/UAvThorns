#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cctk_Functions.h>

#include <math.h>

#include "SpaceMask.h"

/*
 * Self-contained calculation of the center of mass of the matter distribution.
 *
 * The calculation is done for one target at a time, with a single set of
 * integrand grid functions which are recomputed for each target (a scheme
 * similar to the one of VolumeIntegrals_GRMHD).
 * The integrands are reduced with a sum reduction, and the center of mass is
 * obtained by dividing the reduced x^i-weighted integrands by the reduced
 * density integrand.
 */


///////////////////////////////////////////////////////////////////////////////
// Initialization of flags and quantities for the center of mass calculation
///////////////////////////////////////////////////////////////////////////////

void TargetTracker_ComFlags(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // -------------------------------
    // MULTIPATCH
    // -------------------------------

    // Setting the flags has to be done here, cannot be done in PARAMCHECK

    // Trick to check multipatch usage without having to inherit the thorn
    *is_multipatch = CCTK_IsFunctionAliased("MultiPatch_GetDomainSpecification");

    // default (non Llama grid): don't use volume form
    *use_volume_form = 0;

    // TODO: If there's ever a flag in Coordinates that indicates that the volume form is trivial
    //       (e.g. is_cartesian), we could check it and set *use_volume_form = 0 in that case.

    if (*is_multipatch > 0) {
        // If the simulation uses multipatch, make sure that the volume form variable is stored.
        // Supposedly, this should be handled correctly in Coordinates,
        // so this should be redundant with the ParamCheck on store_volume_form.
        CCTK_INT* volume_form_state_ptr = CCTK_VarDataPtr(cctkGH, 0, "Coordinates::volume_form_state");

        if (volume_form_state_ptr != NULL) {
            // No volume form actually stored (should not happen if ParamCheck passed)
            if (*volume_form_state_ptr != 1) {
                CCTK_ERROR("The patch system that you are using does not store the volume form, although you set Coordinates::store_volume_form = yes. "
                           "TargetTracker thorn will not be able to use the volume form for the center of mass integration. "
                           "There is likely a problem in thorn Coordinates.");
            }
            // else, we're good to go
            else {
                CCTK_INFO("Using volume form from multipatch system for the center of mass integration in TargetTracker.");
                *use_volume_form = 1;
            }
        }
        else {
            CCTK_ERROR("Problem acquiring pointer to Coordinates::volume_form_state variable.");
        }
    }


    // -------------------------------
    // SYMMETRY FACTORS
    // -------------------------------

    // Factors to apply to the integrands, depending on the symmetry of the problem.
    // Since the integrand grid functions are private to this thorn, there is no
    // need to declare parities on them through `SetCartSymVN` or similar functions.
    // This is handled internally by the thorn.

    if (CCTK_Equals(com_symmetry, "none")) {
        *sym_factor_drho  = 1;
        *sym_factor_dCoMx = 1;
        *sym_factor_dCoMy = 1;
        *sym_factor_dCoMz = 1;
    }
    else if (CCTK_Equals(com_symmetry, "bitant_x")) {
        *sym_factor_drho  = 2;
        *sym_factor_dCoMx = 0;
        *sym_factor_dCoMy = 2;
        *sym_factor_dCoMz = 2;
    }
    else if (CCTK_Equals(com_symmetry, "bitant_y")) {
        *sym_factor_drho  = 2;
        *sym_factor_dCoMx = 2;
        *sym_factor_dCoMy = 0;
        *sym_factor_dCoMz = 2;
    }
    else if (CCTK_Equals(com_symmetry, "bitant_z")) {
        *sym_factor_drho  = 2;
        *sym_factor_dCoMx = 2;
        *sym_factor_dCoMy = 2;
        *sym_factor_dCoMz = 0;
    }
    else {
        // Should not happen if this function is properly kept consistent with the parameters.
        CCTK_VWARN(CCTK_WARN_ALERT, "Unknown symmetry for the center of mass integrands: '%s'. Defaulting to 'com_symmetry=none'.", com_symmetry);
        *sym_factor_drho  = 1;
        *sym_factor_dCoMx = 1;
        *sym_factor_dCoMy = 1;
        *sym_factor_dCoMz = 1;
    }
}


///////////////////////////////////////////////////////////////////////////////
// Counter management
///////////////////////////////////////////////////////////////////////////////

void TargetTracker_ComCounterZero(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    *com_counter = 0;
}

void TargetTracker_ComCounterZeroRec(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    *com_counter = 0;
}

void TargetTracker_ComCounter(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // Do the calculation only if asked for and on the right iteration
    if (com_every <= 0 || cctk_iteration % com_every != 0) {
        // Make sure the counter is zero: it gates the WHILE loop of the group
        *com_counter = 0;
        return;
    }

    // Register the index of the last (highest-index) target with the calculation
    // activated, so that the WHILE group processes the computed targets from the
    // last one to the first one (iterating on '*com_counter' as 1 + target index).
    *com_counter = 0;
    for (CCTK_INT itarget = nmax_targets - 1 ; itarget >= 0 ; itarget--) {
        if (is_com_computed[itarget]) {
            *com_counter = itarget + 1;
            break;
        }
    }
}

void TargetTracker_ComDecrementCounter(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    (*com_counter)--;
}


///////////////////////////////////////////////////////////////////////////////
// Computation of the integrands for one target
///////////////////////////////////////////////////////////////////////////////

void TargetTracker_ComIntegrand(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // Decode the current target from the counter: the next target to process is
    // given by (*com_counter - 1), which is the highest-index target with the
    // calculation still to be done (see TargetTracker_ComCounter)
    const CCTK_INT itarget = *com_counter - 1;

    // Inactive target: nothing to do.
    // NOTE: The counter iterates over all target indices from the highest computed
    //       one down to 0, so targets in the gaps (not computed) are expected here.
    if (itarget < 0 || !is_com_computed[itarget]) {
        return;
    }


    // Get the mask type and state bits if we need to excise the horizon

    CCTK_INT type_bits     = -1;
    CCTK_INT state_outside = -1;

    if (excise_horizon) {
        type_bits     = SpaceMask_GetTypeBits("mask");
        state_outside = SpaceMask_GetStateBits("mask", "outside");

        if (type_bits < 0) {
            CCTK_ERROR("Thorn AHFinderDirect (or a similar thorn setting the mask) not activated, but excise_horizon requires it.");
        }

        if (state_outside < 0) {
            CCTK_ERROR("Error in obtaining the 'outside' state bits of the mask from SpaceMask.");
        }
    }


    // Volume element: with multipatch we need to multiply by the volume form,
    // else just use the standard Cartesian one

    CCTK_REAL dV_cart = cctk_delta_space[0] * cctk_delta_space[1] * cctk_delta_space[2];

    CCTK_REAL* volume_form = NULL;
    if (*use_volume_form > 0) {
        volume_form = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, "Coordinates::volume_form");
        if (volume_form == NULL) {
            CCTK_ERROR("Problem acquiring pointer to Coordinates::volume_form variable.");
        }
    }


    // Zero out the integrands before filling them
    // WARNING: This zeroes out everything, including the ghost zones,
    //          which are not included in the reductions anyway

    for (int k = 0 ; k < cctk_lsh[2] ; k++) {
    for (int j = 0 ; j < cctk_lsh[1] ; j++) {
    for (int i = 0 ; i < cctk_lsh[0] ; i++) {
        const CCTK_INT ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
        com_integrand_rho[ind] = 0;
        com_integrand_x[ind]   = 0;
        com_integrand_y[ind]   = 0;
        com_integrand_z[ind]   = 0;
    }
    }
    }

    /*
    * Compute the integrands
    *
    * These loops cover the interior of the local box, i.e. they are the C-translation
    * (0-based, half-open) of the Fortran loop
    *     do k = 1+cctk_nghostzones(3), cctk_lsh(3)-cctk_nghostzones(3)
    * (1-based, inclusive, covering C indices [nghostzones, lsh-nghostzones-1]).
    * This includes the grid layer exactly on a symmetry plane (e.g. z=0 for a bitant box),
    * which the reduction weights with 1/2 — necessary for the correct symmetry-normalized
    * integrals.
    *
    * Note that these loops will also exclude at least one layer of points on the physical boundary
    * (see https://lists.einsteintoolkit.org/pipermail/users/2024-September/009465.html),
    * but as of now we choose to live with that.
    *
    *
    * The integrands are zeroed outside/inside spheres (i.e. inside a shell if both are set) as
    *      com_compute_outside_radius < r <= com_compute_inside_radius 
    * with the (<=, <) conventions of VolumeIntegrals_GRMHD.
    * The spheres are centered at the (previous) location of the tracker, the point being
    * to be able to feed a new center of mass to the tracker.
    */

    const CCTK_REAL x_c   = target_loc_x[itarget];
    const CCTK_REAL y_c   = target_loc_y[itarget];
    const CCTK_REAL z_c   = target_loc_z[itarget];
    const CCTK_REAL r_max = com_compute_inside_radius[itarget];
    const CCTK_REAL r_min = com_compute_outside_radius[itarget];

    if (verbose) {
        CCTK_VINFO("CoM pass for target %d: center = (%g,%g,%g), inside_radius = %g, outside_radius = %g",
                   itarget, (double)x_c, (double)y_c, (double)z_c, (double)r_max, (double)r_min);
    }

    for (int k = cctk_nghostzones[2] ; k < cctk_lsh[2] - cctk_nghostzones[2] ; k++) {
    for (int j = cctk_nghostzones[1] ; j < cctk_lsh[1] - cctk_nghostzones[1] ; j++) {
    for (int i = cctk_nghostzones[0] ; i < cctk_lsh[0] - cctk_nghostzones[0] ; i++) {

        const CCTK_INT ind = CCTK_GFINDEX3D(cctkGH, i, j, k);


        // ----- Excision -----
        // Checking if inside the horizon, if asking for it to be excised
        // No need to compute the rest if inside the horizon

        int docalc = 1;
        if (excise_horizon) {
            docalc = SpaceMask_CheckState(space_mask, ind, "mask", "outside");
        }

        // ----- Volume restriction -----
        const CCTK_REAL dx_c = x[ind] - x_c;
        const CCTK_REAL dy_c = y[ind] - y_c;
        const CCTK_REAL dz_c = z[ind] - z_c;
        const CCTK_REAL r2_c = dx_c*dx_c + dy_c*dy_c + dz_c*dz_c;
        docalc = docalc && (r_max  <= 0 || r2_c <= r_max*r_max )
                         && (r_min <= 0 || r2_c >  r_min*r_min);


        // ----- Metric determinant -----

        const CCTK_REAL gxx_l = gxx[ind];
        const CCTK_REAL gxy_l = gxy[ind];
        const CCTK_REAL gxz_l = gxz[ind];
        const CCTK_REAL gyy_l = gyy[ind];
        const CCTK_REAL gyz_l = gyz[ind];
        const CCTK_REAL gzz_l = gzz[ind];

        const CCTK_REAL detgd =       gxx_l * gyy_l * gzz_l
                                + 2 * gxy_l * gxz_l * gyz_l
                                -     gxx_l * gyz_l * gyz_l
                                -     gyy_l * gxz_l * gxz_l
                                -     gzz_l * gxy_l * gxy_l;


        // ----- Volume element -----

        const CCTK_REAL dV = (volume_form != NULL) ? volume_form[ind] : dV_cart;

        const CCTK_REAL vol3 = sqrt(detgd) * dV;


        // ----- Density -----
        // Eulerian energy density rho = n^mu n^nu T_mu_nu

        const CCTK_REAL alph    = alp[ind];
        const CCTK_REAL betax_l = betax[ind];
        const CCTK_REAL betay_l = betay[ind];
        const CCTK_REAL betaz_l = betaz[ind];

        CCTK_REAL rho = eTtt[ind];
        rho = rho - 2 * betax_l * eTtx[ind];
        rho = rho - 2 * betay_l * eTty[ind];
        rho = rho - 2 * betaz_l * eTtz[ind];
        rho = rho + betax_l * betax_l * eTxx[ind];
        rho = rho + betay_l * betay_l * eTyy[ind];
        rho = rho + betaz_l * betaz_l * eTzz[ind];
        rho = rho + 2 * betax_l * betay_l * eTxy[ind];
        rho = rho + 2 * betax_l * betaz_l * eTxz[ind];
        rho = rho + 2 * betay_l * betaz_l * eTyz[ind];
        rho = rho / (alph * alph);

        const CCTK_REAL rho_vol4 = rho * alph * vol3;


        // ----- Integrands -----
        // dCoM^i = rho * x^i * alpha * sqrt(detgd) dV
        // The division by the integral of the density is done after the reduction.
        // We don't use x1 = x - x_track (and add the tracking location back after
        // the reduction), for simplicity.

        if (docalc) {
            com_integrand_rho[ind] = rho_vol4 * (*sym_factor_drho);
            com_integrand_x[ind]   = x[ind] * rho_vol4 * (*sym_factor_dCoMx);
            com_integrand_y[ind]   = y[ind] * rho_vol4 * (*sym_factor_dCoMy);
            com_integrand_z[ind]   = z[ind] * rho_vol4 * (*sym_factor_dCoMz);
        }

    }
    }
    }
}


///////////////////////////////////////////////////////////////////////////////
// Reduction of the integrands, and setting of the center of mass scalars
///////////////////////////////////////////////////////////////////////////////

void TargetTracker_ComDoSum(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // Decode the current target from the counter (see TargetTracker_ComIntegrand)
    const CCTK_INT itarget = *com_counter - 1;

    // Inactive target: nothing to do (should not happen in normal circumstances)
    if (!is_com_computed[itarget]) {
        return;
    }


    // ----- Reduction -----

    enum { num_fields = 4, num_out_vals = 1 };

    // Helper to have a homogeneous interface for reduction output to cactus variables
    CCTK_REAL out_vals[num_fields * num_out_vals];

    CCTK_INT reduction_handle = CCTK_ReductionHandle("sum");
    if (reduction_handle < 0) {
        CCTK_WARN(0, "Could not obtain a handle for sum reduction");
        return;
    }

    const char *const varnames[num_fields] = {
        "TargetTracker::com_integrand_rho",
        "TargetTracker::com_integrand_x",
        "TargetTracker::com_integrand_y",
        "TargetTracker::com_integrand_z"};

    // Get IDs of all variables
    CCTK_INT varid[num_fields];
    for (int i = 0 ; i < num_fields ; i++) {
        varid[i] = CCTK_VarIndex(varnames[i]);
        if (varid[i] < 0) {
            CCTK_VWARN(CCTK_WARN_ABORT, "Could not get index to grid function %s", varnames[i]);
        }
    }

    CCTK_INT ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, num_out_vals,
        CCTK_VARIABLE_REAL, out_vals, num_fields,
        varid[0], varid[1], varid[2], varid[3]);

    if (ierr < 0) {
        CCTK_WARN(0, "Error while reducing the com_integrand grid functions.");
        return;
    }


    // ----- Output -----
    // Center of mass: division by the mass integral (with sanity checks)
    // Note that with symmetry factors, the numerator integrands for some components
    // are deliberately zeroed out, so the division by the (nonzero) mass integral
    // is still well-defined.

    const CCTK_REAL mass = out_vals[0];

    // Sanity checks: the mass integral should be positive (and finite) to make sense
    // (typical failure case: no matter in the integration volume, e.g. a vacuum
    // simulation or the matter having left the sphere)
    const int mass_ok = (mass > 0) && isfinite(mass);

    int com_ok = mass_ok;
    for (int i = 0 ; i < 3 ; i++) {
        com_ok = com_ok && isfinite(out_vals[1 + i]);
    }

    if (com_ok) {
        target_com_x[itarget]    = out_vals[1] / mass;
        target_com_y[itarget]    = out_vals[2] / mass;
        target_com_z[itarget]    = out_vals[3] / mass;
        target_com_mass[itarget] = mass;
    }
    // Failure: keep the previous values
    else {
        CCTK_VWARN(CCTK_WARN_COMPLAIN,
                   "Invalid center of mass integral for target %d (mass = %g). Keeping the previous values (x,y,z) = (%g,%g,%g).",
                   itarget, (double)mass,
                   (double)target_com_x[itarget], (double)target_com_y[itarget], (double)target_com_z[itarget]);
    }


    if (verbose) {
        CCTK_VINFO("Center of mass of target %d: (x,y,z) = (%g,%g,%g), mass = %g",
                   itarget,
                   (double)target_com_x[itarget], (double)target_com_y[itarget], (double)target_com_z[itarget],
                   (double)target_com_mass[itarget]);
    }
}