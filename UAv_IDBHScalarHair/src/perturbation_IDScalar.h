// TODO: Not quite sure which includes are superfluous, and if we really the CCTK parts, but play safe
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

/*
Helpers to compute perturbations to be added to quantities (scalar field, metric, ...).
The idea is to factorize the code a bit, to make it easier to change and complement,
without spoiling legibility of the main code.

Due to the nature of the loops, we act here on the individual values at a point, not on grid functions.
*/


////////////////////////////////////////////////////////////////////////////////
// Helper struc to pass the coordinates
////////////////////////////////////////////////////////////////////////////////
struct PertCoords {
    // const CCTK_REAL x;
    // const CCTK_REAL y;
    // const CCTK_REAL z;

    const CCTK_REAL rho;
    const CCTK_REAL R;

    // const CCTK_REAL th;
    const CCTK_REAL ph;
};

////////////////////////////////////////////////////////////////////////////////
// Middle/low level function that switches to the precise perturbation function,
// computed as a the product of a radial seed, and an angular part.
// Note that the angular part includes multiplicative radial factors for 
// regularity.
// The radial coordinate might be R or rho, depending on the cases.
////////////////////////////////////////////////////////////////////////////////
// TODO: As a possible future extension, have different topologies of 
//       perturbations, e.g. cylindrical, toroidal, spherical, ...
// TODO: If we add many different options, consider factoring out the actual
//       computations to dedicated functions.
//       (radial vs. angular; specific seeds; ...)
////////////////////////////////////////////////////////////////////////////////

static inline CCTK_REAL UAv_IDScalar__ComputePerturbation (CCTK_ARGUMENTS, const struct PertCoords pertCoords) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    CCTK_REAL radial, angular;

    // ----- Compute the radial part of the perturbation -----
    if      (CCTK_EQUALS(perturbation_radial_seed, "constant")) {
        
        radial = perturbation_amplitude;
    
    }
    else if (CCTK_EQUALS(perturbation_radial_seed, "gaussian")) {
        
        const CCTK_REAL A = perturbation_amplitude;
        const CCTK_REAL R0 = perturbation_R0;
        const CCTK_REAL sigma = perturbation_width;

        const CCTK_REAL arg = (pertCoords.R - R0) / sigma;
        radial = A * exp(- 0.5 * arg * arg);
    
    }
    else {
    
        CCTK_VERROR("Unknown PertRadialSeed = %s in UAv_IDScalar__ComputePerturbation.", perturbation_radial_seed);
    
    }

    // ----- Compute the angular part of the perturbation -----
    
    // For now, simple "cylindrical" azimuthal l=m mode, i.e. rho^m * cos(m phi).
    // This may not satisfy all symmetries and regularity conditions,
    // but the perturbation is supposed to be small, and it's a perturbation.
    // (NOTE: for m==0, this is actually spherical.)
    {

        const CCTK_INT m = perturbation_azimuthal_number;
        // NOTE: We add factors of mu so that the angular part is dimensionless.
        angular = m==0 ? 1 : pow(pertCoords.rho*mu, m) * cos(m * pertCoords.ph);
    
    }

    return radial * angular;
}

////////////////////////////////////////////////////////////////////////////////
// High level function to return the perturbed value at a given point,
// given the base value and the coordinates.
// It returns the value of the field after perturbation.
//
// It basically just switches on the perturbation type and redirects the call.
////////////////////////////////////////////////////////////////////////////////

static inline CCTK_REAL UAv_IDScalar_Perturb (CCTK_ARGUMENTS, const CCTK_REAL base_value, const struct PertCoords pertCoords) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // Redundancy in parameters, to force the user to be explicit and not use defaults.
    // Also, we probably shouldn't hit that here
    if (CCTK_EQUALS(perturbation_field, "none")) {
        return base_value;
    }
    
    CCTK_REAL result = base_value;

    // ----- Switch on the perturbation type -----
    if      (CCTK_EQUALS(perturbation_type, "none")) {
        // do nothing, return the base value
    }
    else if (CCTK_EQUALS(perturbation_type, "relative")) {
        result *= (1. + UAv_IDScalar__ComputePerturbation(CCTK_PASS_CTOC, pertCoords));
    }
    else if (CCTK_EQUALS(perturbation_type, "multiplicative")) {
        result *= UAv_IDScalar__ComputePerturbation(CCTK_PASS_CTOC, pertCoords);
    }
    else if (CCTK_EQUALS(perturbation_type, "additive")) {
        result += UAv_IDScalar__ComputePerturbation(CCTK_PASS_CTOC, pertCoords);
    }
    else {
        CCTK_VERROR("Unknown perturbation type = %s in UAv_IDScalar__Perturb.", perturbation_type);
    }

    return result;
}