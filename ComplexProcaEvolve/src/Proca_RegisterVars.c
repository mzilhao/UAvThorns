
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void Proca_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr = 0, group, rhs, var;

  // register evolution and rhs gridfunction groups with MoL

  /* metric and extrinsic curvature */
  group = CCTK_GroupIndex("ADMBase::lapse");
  ierr += MoLRegisterSaveAndRestoreGroup(group);
  group = CCTK_GroupIndex("ADMBase::shift");
  ierr += MoLRegisterSaveAndRestoreGroup(group);
  group = CCTK_GroupIndex("ADMBase::metric");
  ierr += MoLRegisterSaveAndRestoreGroup(group);
  group = CCTK_GroupIndex("ADMBase::curv");
  ierr += MoLRegisterSaveAndRestoreGroup(group);

  /* Ei and rhs_Ei */
  group = CCTK_GroupIndex("ProcaBase::E1i");
  rhs   = CCTK_GroupIndex("ProcaEvolve::rhs_E1i");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  /* Ai and rhs_Ai */
  group = CCTK_GroupIndex("ProcaBase::A1i");
  rhs   = CCTK_GroupIndex("ProcaEvolve::rhs_A1i");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  /* Aphi and rhs_Aphi */
  var   = CCTK_VarIndex("ProcaBase::Aphi1");
  rhs   = CCTK_VarIndex("ProcaEvolve::rhs_Aphi1");
  ierr += MoLRegisterEvolved(var, rhs);

  /* Zeta and rhs_Zeta */
  var   = CCTK_VarIndex("ProcaBase::Zeta1");
  rhs   = CCTK_VarIndex("ProcaEvolve::rhs_Zeta1");
  ierr += MoLRegisterEvolved(var, rhs);
  
  /* Ei and rhs_Ei */
  group = CCTK_GroupIndex("ProcaBase::E2i");
  rhs   = CCTK_GroupIndex("ProcaEvolve::rhs_E2i");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  /* Ai and rhs_Ai */
  group = CCTK_GroupIndex("ProcaBase::A2i");
  rhs   = CCTK_GroupIndex("ProcaEvolve::rhs_A2i");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  /* Aphi and rhs_Aphi */
  var   = CCTK_VarIndex("ProcaBase::Aphi2");
  rhs   = CCTK_VarIndex("ProcaEvolve::rhs_Aphi2");
  ierr += MoLRegisterEvolved(var, rhs);

  /* Zeta and rhs_Zeta */
  var   = CCTK_VarIndex("ProcaBase::Zeta2");
  rhs   = CCTK_VarIndex("ProcaEvolve::rhs_Zeta2");
  ierr += MoLRegisterEvolved(var, rhs);

  if (ierr) CCTK_ERROR("Problems registering with MoL");

}
