
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"


subroutine MagScalar_calc_rhs( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  ! Fundamental variables
  CCTK_REAL                alph, beta(3)
  CCTK_REAL                hh(3,3), hu(3,3), trk, dethh, ch
  CCTK_REAL                lE(3), lA(3), lAphi, lZeta
  CCTK_REAL                lphi1, lphi2, lKphi1, lKphi2

  ! First derivatives
  CCTK_REAL                d1_alph(3), d1_beta(3,3)
  CCTK_REAL                d1_hh(3,3,3), d1_ch(3)
  CCTK_REAL                d1_lE(3,3), d1_lA(3,3), d1_lZeta(3), d1_lAphi(3)
  CCTK_REAL                d1_lphi1(3), d1_lphi2(3), d1_lKphi1(3), d1_lKphi2(3)

  ! Second derivatives
  CCTK_REAL                d2_lA(3,3,3)
  CCTK_REAL                d2_lphi1(3,3), d2_lphi2(3,3)

  ! Advection derivatives
  CCTK_REAL                ad1_lE(3), ad1_lA(3), ad1_lZeta, ad1_lAphi
  CCTK_REAL                ad1_lphi1, ad1_lphi2, ad1_lKphi1, ad1_lKphi2
  CCTK_REAL                d1_f(3)   ! Place holder for the advection derivs

  ! Auxiliary variables
  CCTK_REAL                cf1(3,3,3), cf2(3,3,3)
  CCTK_REAL                tr_dalp_dphi1, tr_cd2_phi1, tr_dch_dphi1
  CCTK_REAL                tr_dalp_dphi2, tr_cd2_phi2, tr_dch_dphi2
  CCTK_REAL                tr_A_A, tr_A_dphi1, tr_A_dphi2
  CCTK_REAL                rho_e, J_i(3)

  ! Covaraint derivatives
  CCTK_REAL                cd_lA(3,3), cd_dA(3,3,3)
  CCTK_REAL                cd2_lphi1(3,3), cd2_lphi2(3,3)

  ! Right hand sides
  CCTK_REAL                rhs_lE(3), rhs_lA(3), rhs_lZeta, rhs_lAphi
  CCTK_REAL                rhs_lphi1, rhs_lphi2 
  CCTK_REAL                rhs_lKphi1, rhs_lKphi2

  ! Misc variables
  CCTK_REAL                dx12, dy12, dz12, dxsq12, dysq12, dzsq12,         &
                           dxdy144, dxdz144, dydz144
  CCTK_REAL                odx60, ody60, odz60, odxsq180, odysq180, odzsq180,&
                           odxdy3600, odxdz3600, odydz3600
  CCTK_INT                 i, j, k
  CCTK_INT                 di, dj, dk
  CCTK_REAL, parameter ::  one = 1
  CCTK_INT                 a, b, c, m, n


  dx12 = 12*CCTK_DELTA_SPACE(1)
  dy12 = 12*CCTK_DELTA_SPACE(2)
  dz12 = 12*CCTK_DELTA_SPACE(3)

  dxsq12 = 12*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(1)
  dysq12 = 12*CCTK_DELTA_SPACE(2)*CCTK_DELTA_SPACE(2)
  dzsq12 = 12*CCTK_DELTA_SPACE(3)*CCTK_DELTA_SPACE(3)

  dxdy144 = 144*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2)
  dxdz144 = 144*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(3)
  dydz144 = 144*CCTK_DELTA_SPACE(2)*CCTK_DELTA_SPACE(3)


  odx60 = 1 / (60 * CCTK_DELTA_SPACE(1))
  ody60 = 1 / (60 * CCTK_DELTA_SPACE(2))
  odz60 = 1 / (60 * CCTK_DELTA_SPACE(3))

  odxsq180 = 1 / (180*CCTK_DELTA_SPACE(1)**2)
  odysq180 = 1 / (180*CCTK_DELTA_SPACE(2)**2)
  odzsq180 = 1 / (180*CCTK_DELTA_SPACE(3)**2)

  odxdy3600 = 1 / (3600*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2))
  odxdz3600 = 1 / (3600*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(3))
  odydz3600 = 1 / (3600*CCTK_DELTA_SPACE(2)*CCTK_DELTA_SPACE(3))


  ! convert ADM variables to BSSN-like ones
  call MagScalar_adm2bssn(CCTK_PASS_FTOF)


  !$OMP PARALLEL DO COLLAPSE(3) &
  !$OMP PRIVATE(alph, beta, hh, hu, trk, dethh, ch,&
  !$OMP lE, lA, lAphi, lZeta,&
  !$OMP lphi1, lphi2, lKphi1, lKphi2,&
  !$OMP d1_alph, d1_beta, d1_hh, d1_ch,&
  !$OMP d1_lE, d1_lA, d1_lZeta, d1_lAphi,&
  !$OMP d2_lA, ad1_lE, ad1_lA, ad1_lZeta, ad1_lAphi,&
  !$OMP d1_lphi1, d1_lphi2, d1_lKphi1, d1_lKphi2,&
  !$OMP d2_lphi1, d2_lphi2, ad1_lphi1, ad1_lphi2, ad1_lKphi1, ad1_lKphi2,&
  !$OMP d1_f, cf1, cf2, cd_lA, cd_dA, cd2_lphi1, cd2_lphi2,&
  !$OMP tr_dalp_dphi1, tr_cd2_phi1, tr_dch_dphi1,&
  !$OMP tr_dalp_dphi2, tr_cd2_phi2, tr_dch_dphi2,&
  !$OMP rho_e, J_i, tr_A_A, tr_A_dphi1, tr_A_dphi2,&
  !$OMP rhs_lE, rhs_lA, rhs_lZeta, rhs_lAphi,&
  !$OMP rhs_lphi1, rhs_lphi2, rhs_lKphi1, rhs_lKphi2,&
  !$OMP i, j, k,&
  !$OMP di, dj, dk,&
  !$OMP a, b, c, m, n)
  do k = 1+cctk_nghostzones(3), cctk_lsh(3)-cctk_nghostzones(3)
  do j = 1+cctk_nghostzones(2), cctk_lsh(2)-cctk_nghostzones(2)
  do i = 1+cctk_nghostzones(1), cctk_lsh(1)-cctk_nghostzones(1)

    !------------ Get local variables ----------
    ch        = chi(i,j,k)

    hh(1,1)   = hxx(i,j,k)
    hh(1,2)   = hxy(i,j,k)
    hh(1,3)   = hxz(i,j,k)
    hh(2,2)   = hyy(i,j,k)
    hh(2,3)   = hyz(i,j,k)
    hh(3,3)   = hzz(i,j,k)
    hh(2,1)   = hh(1,2)
    hh(3,1)   = hh(1,3)
    hh(3,2)   = hh(2,3)

    trk       = tracek(i,j,k)

    alph      = alp(i,j,k)

    beta(1)   = betax(i,j,k)
    beta(2)   = betay(i,j,k)
    beta(3)   = betaz(i,j,k)

    lE(1)     = Ex(i,j,k)
    lE(2)     = Ey(i,j,k)
    lE(3)     = Ez(i,j,k)

    lA(1)     = Ax(i,j,k)
    lA(2)     = Ay(i,j,k)
    lA(3)     = Az(i,j,k)

    lZeta     = Zeta(i,j,k)
    lAphi     = Aphi(i,j,k)

    lphi1   = phi1(i,j,k)
    lphi2   = phi2(i,j,k)

    lKphi1  = Kphi1(i,j,k)
    lKphi2  = Kphi2(i,j,k)
    !-------------------------------------------


    !------------ Invert 3-metric ----------------
    ! NOTE: deth = 1 by construction, but that is not satisfied numerically
    dethh =       hh(1,1) * hh(2,2) * hh(3,3)                              &
            + 2 * hh(1,2) * hh(1,3) * hh(2,3)                              &
            -     hh(1,1) * hh(2,3) ** 2                                   &
            -     hh(2,2) * hh(1,3) ** 2                                   &
            -     hh(3,3) * hh(1,2) ** 2
    hu(1,1) = (hh(2,2) * hh(3,3) - hh(2,3) ** 2     ) / dethh
    hu(2,2) = (hh(1,1) * hh(3,3) - hh(1,3) ** 2     ) / dethh
    hu(3,3) = (hh(1,1) * hh(2,2) - hh(1,2) ** 2     ) / dethh
    hu(1,2) = (hh(1,3) * hh(2,3) - hh(1,2) * hh(3,3)) / dethh
    hu(1,3) = (hh(1,2) * hh(2,3) - hh(1,3) * hh(2,2)) / dethh
    hu(2,3) = (hh(1,3) * hh(1,2) - hh(2,3) * hh(1,1)) / dethh
    hu(2,1) = hu(1,2)
    hu(3,1) = hu(1,3)
    hu(3,2) = hu(2,3)
    !-------------------------------------------


    if (derivs_order == 4) then

      !------------- Centered 1st derivatives -----------

      ! d1_ch(3)
      d1_ch(1) = (   -chi(i+2,j,k) + 8*chi(i+1,j,k)                          &
                  - 8*chi(i-1,j,k) +   chi(i-2,j,k) ) / dx12

      d1_ch(2) = (   -chi(i,j+2,k) + 8*chi(i,j+1,k)                          &
                  - 8*chi(i,j-1,k) +   chi(i,j-2,k) ) / dy12

      d1_ch(3) = (   -chi(i,j,k+2) + 8*chi(i,j,k+1)                          &
                  - 8*chi(i,j,k-1) +   chi(i,j,k-2) ) / dz12


      ! d1_hh(3,3,3)
      d1_hh(1,1,1) = (   -hxx(i+2,j,k) + 8*hxx(i+1,j,k)                      &
                      - 8*hxx(i-1,j,k) +   hxx(i-2,j,k) ) / dx12
      d1_hh(1,2,1) = (   -hxy(i+2,j,k) + 8*hxy(i+1,j,k)                      &
                      - 8*hxy(i-1,j,k) +   hxy(i-2,j,k) ) / dx12
      d1_hh(1,3,1) = (   -hxz(i+2,j,k) + 8*hxz(i+1,j,k)                      &
                      - 8*hxz(i-1,j,k) +   hxz(i-2,j,k) ) / dx12
      d1_hh(2,2,1) = (   -hyy(i+2,j,k) + 8*hyy(i+1,j,k)                      &
                      - 8*hyy(i-1,j,k) +   hyy(i-2,j,k) ) / dx12
      d1_hh(2,3,1) = (   -hyz(i+2,j,k) + 8*hyz(i+1,j,k)                      &
                      - 8*hyz(i-1,j,k) +   hyz(i-2,j,k) ) / dx12
      d1_hh(3,3,1) = (   -hzz(i+2,j,k) + 8*hzz(i+1,j,k)                      &
                      - 8*hzz(i-1,j,k) +   hzz(i-2,j,k) ) / dx12

      d1_hh(1,1,2) = (   -hxx(i,j+2,k) + 8*hxx(i,j+1,k)                      &
                      - 8*hxx(i,j-1,k) +   hxx(i,j-2,k) ) / dy12
      d1_hh(1,2,2) = (   -hxy(i,j+2,k) + 8*hxy(i,j+1,k)                      &
                      - 8*hxy(i,j-1,k) +   hxy(i,j-2,k) ) / dy12
      d1_hh(1,3,2) = (   -hxz(i,j+2,k) + 8*hxz(i,j+1,k)                      &
                      - 8*hxz(i,j-1,k) +   hxz(i,j-2,k) ) / dy12
      d1_hh(2,2,2) = (   -hyy(i,j+2,k) + 8*hyy(i,j+1,k)                      &
                      - 8*hyy(i,j-1,k) +   hyy(i,j-2,k) ) / dy12
      d1_hh(2,3,2) = (   -hyz(i,j+2,k) + 8*hyz(i,j+1,k)                      &
                      - 8*hyz(i,j-1,k) +   hyz(i,j-2,k) ) / dy12
      d1_hh(3,3,2) = (   -hzz(i,j+2,k) + 8*hzz(i,j+1,k)                      &
                      - 8*hzz(i,j-1,k) +   hzz(i,j-2,k) ) / dy12

      d1_hh(1,1,3) = (   -hxx(i,j,k+2) + 8*hxx(i,j,k+1)                      &
                      - 8*hxx(i,j,k-1) +   hxx(i,j,k-2) ) / dz12
      d1_hh(1,2,3) = (   -hxy(i,j,k+2) + 8*hxy(i,j,k+1)                      &
                      - 8*hxy(i,j,k-1) +   hxy(i,j,k-2) ) / dz12
      d1_hh(1,3,3) = (   -hxz(i,j,k+2) + 8*hxz(i,j,k+1)                      &
                      - 8*hxz(i,j,k-1) +   hxz(i,j,k-2) ) / dz12
      d1_hh(2,2,3) = (   -hyy(i,j,k+2) + 8*hyy(i,j,k+1)                      &
                      - 8*hyy(i,j,k-1) +   hyy(i,j,k-2) ) / dz12
      d1_hh(2,3,3) = (   -hyz(i,j,k+2) + 8*hyz(i,j,k+1)                      &
                      - 8*hyz(i,j,k-1) +   hyz(i,j,k-2) ) / dz12
      d1_hh(3,3,3) = (   -hzz(i,j,k+2) + 8*hzz(i,j,k+1)                      &
                      - 8*hzz(i,j,k-1) +   hzz(i,j,k-2) ) / dz12

      d1_hh(2,1,:) = d1_hh(1,2,:)
      d1_hh(3,1,:) = d1_hh(1,3,:)
      d1_hh(3,2,:) = d1_hh(2,3,:)


      ! d1_alph(3)
      d1_alph(1) = (   -alp(i+2,j,k) + 8*alp(i+1,j,k)                        &
                    - 8*alp(i-1,j,k) +   alp(i-2,j,k) ) / dx12

      d1_alph(2) = (   -alp(i,j+2,k) + 8*alp(i,j+1,k)                        &
                    - 8*alp(i,j-1,k) +   alp(i,j-2,k) ) / dy12

      d1_alph(3) = (   -alp(i,j,k+2) + 8*alp(i,j,k+1)                        &
                    - 8*alp(i,j,k-1) +   alp(i,j,k-2) ) / dz12

      ! d1_beta (3,3)
      d1_beta(1,1)  = (   -betax(i+2,j,k) + 8*betax(i+1,j,k)                 &
                       - 8*betax(i-1,j,k) +   betax(i-2,j,k) ) / dx12
      d1_beta(2,1)  = (   -betay(i+2,j,k) + 8*betay(i+1,j,k)                 &
                       - 8*betay(i-1,j,k) +   betay(i-2,j,k) ) / dx12
      d1_beta(3,1)  = (   -betaz(i+2,j,k) + 8*betaz(i+1,j,k)                 &
                       - 8*betaz(i-1,j,k) +   betaz(i-2,j,k) ) / dx12

      d1_beta(1,2)  = (   -betax(i,j+2,k) + 8*betax(i,j+1,k)                 &
                       - 8*betax(i,j-1,k) +   betax(i,j-2,k) ) / dy12
      d1_beta(2,2)  = (   -betay(i,j+2,k) + 8*betay(i,j+1,k)                 &
                       - 8*betay(i,j-1,k) +   betay(i,j-2,k) ) / dy12
      d1_beta(3,2)  = (   -betaz(i,j+2,k) + 8*betaz(i,j+1,k)                 &
                       - 8*betaz(i,j-1,k) +   betaz(i,j-2,k) ) / dy12

      d1_beta(1,3)  = (   -betax(i,j,k+2) + 8*betax(i,j,k+1)                 &
                       - 8*betax(i,j,k-1) +   betax(i,j,k-2) ) / dz12
      d1_beta(2,3)  = (   -betay(i,j,k+2) + 8*betay(i,j,k+1)                 &
                       - 8*betay(i,j,k-1) +   betay(i,j,k-2) ) / dz12
      d1_beta(3,3)  = (   -betaz(i,j,k+2) + 8*betaz(i,j,k+1)                 &
                       - 8*betaz(i,j,k-1) +   betaz(i,j,k-2) ) / dz12

      ! d1_lE(3,3)
      d1_lE(1,1) = (   -Ex(i+2,j,k) + 8*Ex(i+1,j,k)               &
                    - 8*Ex(i-1,j,k) +   Ex(i-2,j,k) ) / dx12
      d1_lE(2,1) = (   -Ey(i+2,j,k) + 8*Ey(i+1,j,k)               &
                    - 8*Ey(i-1,j,k) +   Ey(i-2,j,k) ) / dx12
      d1_lE(3,1) = (   -Ez(i+2,j,k) + 8*Ez(i+1,j,k)               &
                    - 8*Ez(i-1,j,k) +   Ez(i-2,j,k) ) / dx12

      d1_lE(1,2) = (   -Ex(i,j+2,k) + 8*Ex(i,j+1,k)               &
                    - 8*Ex(i,j-1,k) +   Ex(i,j-2,k) ) / dy12
      d1_lE(2,2) = (   -Ey(i,j+2,k) + 8*Ey(i,j+1,k)               &
                    - 8*Ey(i,j-1,k) +   Ey(i,j-2,k) ) / dy12
      d1_lE(3,2) = (   -Ez(i,j+2,k) + 8*Ez(i,j+1,k)               &
                    - 8*Ez(i,j-1,k) +   Ez(i,j-2,k) ) / dy12

      d1_lE(1,3) = (   -Ex(i,j,k+2) + 8*Ex(i,j,k+1)               &
                    - 8*Ex(i,j,k-1) +   Ex(i,j,k-2) ) / dz12
      d1_lE(2,3) = (   -Ey(i,j,k+2) + 8*Ey(i,j,k+1)               &
                    - 8*Ey(i,j,k-1) +   Ey(i,j,k-2) ) / dz12
      d1_lE(3,3) = (   -Ez(i,j,k+2) + 8*Ez(i,j,k+1)               &
                    - 8*Ez(i,j,k-1) +   Ez(i,j,k-2) ) / dz12

      ! d1_lA(3,3)
      d1_lA(1,1) = (   -Ax(i+2,j,k) + 8*Ax(i+1,j,k)               &
                    - 8*Ax(i-1,j,k) +   Ax(i-2,j,k) ) / dx12
      d1_lA(2,1) = (   -Ay(i+2,j,k) + 8*Ay(i+1,j,k)               &
                    - 8*Ay(i-1,j,k) +   Ay(i-2,j,k) ) / dx12
      d1_lA(3,1) = (   -Az(i+2,j,k) + 8*Az(i+1,j,k)               &
                    - 8*Az(i-1,j,k) +   Az(i-2,j,k) ) / dx12

      d1_lA(1,2) = (   -Ax(i,j+2,k) + 8*Ax(i,j+1,k)               &
                    - 8*Ax(i,j-1,k) +   Ax(i,j-2,k) ) / dy12
      d1_lA(2,2) = (   -Ay(i,j+2,k) + 8*Ay(i,j+1,k)               &
                    - 8*Ay(i,j-1,k) +   Ay(i,j-2,k) ) / dy12
      d1_lA(3,2) = (   -Az(i,j+2,k) + 8*Az(i,j+1,k)               &
                    - 8*Az(i,j-1,k) +   Az(i,j-2,k) ) / dy12

      d1_lA(1,3) = (   -Ax(i,j,k+2) + 8*Ax(i,j,k+1)               &
                    - 8*Ax(i,j,k-1) +   Ax(i,j,k-2) ) / dz12
      d1_lA(2,3) = (   -Ay(i,j,k+2) + 8*Ay(i,j,k+1)               &
                    - 8*Ay(i,j,k-1) +   Ay(i,j,k-2) ) / dz12
      d1_lA(3,3) = (   -Az(i,j,k+2) + 8*Az(i,j,k+1)               &
                    - 8*Az(i,j,k-1) +   Az(i,j,k-2) ) / dz12

      ! d1_lZeta(3)
      d1_lZeta(1) = (   -Zeta(i+2,j,k) + 8*Zeta(i+1,j,k)                        &
                     - 8*Zeta(i-1,j,k) +   Zeta(i-2,j,k) ) / dx12

      d1_lZeta(2) = (   -Zeta(i,j+2,k) + 8*Zeta(i,j+1,k)                        &
                     - 8*Zeta(i,j-1,k) +   Zeta(i,j-2,k) ) / dy12

      d1_lZeta(3) = (   -Zeta(i,j,k+2) + 8*Zeta(i,j,k+1)                        &
                     - 8*Zeta(i,j,k-1) +   Zeta(i,j,k-2) ) / dz12

      ! d1_lAphi(3)
      d1_lAphi(1) = (   -Aphi(i+2,j,k) + 8*Aphi(i+1,j,k)                        &
                     - 8*Aphi(i-1,j,k) +   Aphi(i-2,j,k) ) / dx12

      d1_lAphi(2) = (   -Aphi(i,j+2,k) + 8*Aphi(i,j+1,k)                        &
                     - 8*Aphi(i,j-1,k) +   Aphi(i,j-2,k) ) / dy12

      d1_lAphi(3) = (   -Aphi(i,j,k+2) + 8*Aphi(i,j,k+1)                        &
                     - 8*Aphi(i,j,k-1) +   Aphi(i,j,k-2) ) / dz12

      ! d1_lphi1(3)
      d1_lphi1(1) = (   -phi1(i+2,j,k) + 8*phi1(i+1,j,k)                          &
                    - 8*phi1(i-1,j,k) +   phi1(i-2,j,k) ) / dx12

      d1_lphi1(2) = (   -phi1(i,j+2,k) + 8*phi1(i,j+1,k)                          &
                    - 8*phi1(i,j-1,k) +   phi1(i,j-2,k) ) / dy12

      d1_lphi1(3) = (   -phi1(i,j,k+2) + 8*phi1(i,j,k+1)                          &
                    - 8*phi1(i,j,k-1) +   phi1(i,j,k-2) ) / dz12

      ! d1_lphi2(3)
      d1_lphi2(1) = (   -phi2(i+2,j,k) + 8*phi2(i+1,j,k)                          &
                    - 8*phi2(i-1,j,k) +   phi2(i-2,j,k) ) / dx12

      d1_lphi2(2) = (   -phi2(i,j+2,k) + 8*phi2(i,j+1,k)                          &
                    - 8*phi2(i,j-1,k) +   phi2(i,j-2,k) ) / dy12

      d1_lphi2(3) = (   -phi2(i,j,k+2) + 8*phi2(i,j,k+1)                          &
                    - 8*phi2(i,j,k-1) +   phi2(i,j,k-2) ) / dz12

      ! d1_lKphi1(3)
      d1_lKphi1(1) = (   -Kphi1(i+2,j,k) + 8*Kphi1(i+1,j,k)                          &
                      - 8*Kphi1(i-1,j,k) +   Kphi1(i-2,j,k) ) / dx12

      d1_lKphi1(2) = (   -Kphi1(i,j+2,k) + 8*Kphi1(i,j+1,k)                          &
                      - 8*Kphi1(i,j-1,k) +   Kphi1(i,j-2,k) ) / dy12

      d1_lKphi1(3) = (   -Kphi1(i,j,k+2) + 8*Kphi1(i,j,k+1)                          &
                      - 8*Kphi1(i,j,k-1) +   Kphi1(i,j,k-2) ) / dz12

      ! d1_lKphi2(3)
      d1_lKphi2(1) = (   -Kphi2(i+2,j,k) + 8*Kphi2(i+1,j,k)                          &
                      - 8*Kphi2(i-1,j,k) +   Kphi2(i-2,j,k) ) / dx12

      d1_lKphi2(2) = (   -Kphi2(i,j+2,k) + 8*Kphi2(i,j+1,k)                          &
                      - 8*Kphi2(i,j-1,k) +   Kphi2(i,j-2,k) ) / dy12

      d1_lKphi2(3) = (   -Kphi2(i,j,k+2) + 8*Kphi2(i,j,k+1)                          &
                      - 8*Kphi2(i,j,k-1) +   Kphi2(i,j,k-2) ) / dz12

      !--------------------------------------------------


      !------------- Centered 2nd derivatives -----------

      ! d2_lA(3,3,3)
      d2_lA(1,1,1) = (  -Ax(i+2,j,k) + 16*Ax(i+1,j,k) - 30*Ax(i,j,k) &
                    + 16*Ax(i-1,j,k) -    Ax(i-2,j,k) ) / dxsq12
      d2_lA(2,1,1) = (  -Ay(i+2,j,k) + 16*Ay(i+1,j,k) - 30*Ay(i,j,k) &
                    + 16*Ay(i-1,j,k) -    Ay(i-2,j,k) ) / dxsq12
      d2_lA(3,1,1) = (  -Az(i+2,j,k) + 16*Az(i+1,j,k) - 30*Az(i,j,k) &
                    + 16*Az(i-1,j,k) -    Az(i-2,j,k) ) / dxsq12

      d2_lA(1,2,2) = (  -Ax(i,j+2,k) + 16*Ax(i,j+1,k) - 30*Ax(i,j,k) &
                    + 16*Ax(i,j-1,k) -    Ax(i,j-2,k) ) / dysq12
      d2_lA(2,2,2) = (  -Ay(i,j+2,k) + 16*Ay(i,j+1,k) - 30*Ay(i,j,k) &
                    + 16*Ay(i,j-1,k) -    Ay(i,j-2,k) ) / dysq12
      d2_lA(3,2,2) = (  -Az(i,j+2,k) + 16*Az(i,j+1,k) - 30*Az(i,j,k) &
                    + 16*Az(i,j-1,k) -    Az(i,j-2,k) ) / dysq12

      d2_lA(1,3,3) = (  -Ax(i,j,k+2) + 16*Ax(i,j,k+1) - 30*Ax(i,j,k) &
                    + 16*Ax(i,j,k-1) -    Ax(i,j,k-2) ) / dzsq12
      d2_lA(2,3,3) = (  -Ay(i,j,k+2) + 16*Ay(i,j,k+1) - 30*Ay(i,j,k) &
                    + 16*Ay(i,j,k-1) -    Ay(i,j,k-2) ) / dzsq12
      d2_lA(3,3,3) = (  -Az(i,j,k+2) + 16*Az(i,j,k+1) - 30*Az(i,j,k) &
                    + 16*Az(i,j,k-1) -    Az(i,j,k-2) ) / dzsq12

      d2_lA(1,1,2) = (  -Ax(i-2,j+2,k) +  8*Ax(i-1,j+2,k) -  8*Ax(i+1,j+2,k) +   Ax(i+2,j+2,k) &
                     + 8*Ax(i-2,j+1,k) - 64*Ax(i-1,j+1,k) + 64*Ax(i+1,j+1,k) - 8*Ax(i+2,j+1,k) &
                     - 8*Ax(i-2,j-1,k) + 64*Ax(i-1,j-1,k) - 64*Ax(i+1,j-1,k) + 8*Ax(i+2,j-1,k) &
                     +   Ax(i-2,j-2,k) -  8*Ax(i-1,j-2,k) +  8*Ax(i+1,j-2,k) -   Ax(i+2,j-2,k) ) / dxdy144
      d2_lA(2,1,2) = (  -Ay(i-2,j+2,k) +  8*Ay(i-1,j+2,k) -  8*Ay(i+1,j+2,k) +   Ay(i+2,j+2,k) &
                     + 8*Ay(i-2,j+1,k) - 64*Ay(i-1,j+1,k) + 64*Ay(i+1,j+1,k) - 8*Ay(i+2,j+1,k) &
                     - 8*Ay(i-2,j-1,k) + 64*Ay(i-1,j-1,k) - 64*Ay(i+1,j-1,k) + 8*Ay(i+2,j-1,k) &
                     +   Ay(i-2,j-2,k) -  8*Ay(i-1,j-2,k) +  8*Ay(i+1,j-2,k) -   Ay(i+2,j-2,k) ) / dxdy144
      d2_lA(3,1,2) = (  -Az(i-2,j+2,k) +  8*Az(i-1,j+2,k) -  8*Az(i+1,j+2,k) +   Az(i+2,j+2,k) &
                     + 8*Az(i-2,j+1,k) - 64*Az(i-1,j+1,k) + 64*Az(i+1,j+1,k) - 8*Az(i+2,j+1,k) &
                     - 8*Az(i-2,j-1,k) + 64*Az(i-1,j-1,k) - 64*Az(i+1,j-1,k) + 8*Az(i+2,j-1,k) &
                     +   Az(i-2,j-2,k) -  8*Az(i-1,j-2,k) +  8*Az(i+1,j-2,k) -   Az(i+2,j-2,k) ) / dxdy144

      d2_lA(1,1,3) = (  -Ax(i-2,j,k+2) +  8*Ax(i-1,j,k+2) -  8*Ax(i+1,j,k+2) +   Ax(i+2,j,k+2) &
                     + 8*Ax(i-2,j,k+1) - 64*Ax(i-1,j,k+1) + 64*Ax(i+1,j,k+1) - 8*Ax(i+2,j,k+1) &
                     - 8*Ax(i-2,j,k-1) + 64*Ax(i-1,j,k-1) - 64*Ax(i+1,j,k-1) + 8*Ax(i+2,j,k-1) &
                     +   Ax(i-2,j,k-2) -  8*Ax(i-1,j,k-2) +  8*Ax(i+1,j,k-2) -   Ax(i+2,j,k-2) ) / dxdz144
      d2_lA(2,1,3) = (  -Ay(i-2,j,k+2) +  8*Ay(i-1,j,k+2) -  8*Ay(i+1,j,k+2) +   Ay(i+2,j,k+2) &
                     + 8*Ay(i-2,j,k+1) - 64*Ay(i-1,j,k+1) + 64*Ay(i+1,j,k+1) - 8*Ay(i+2,j,k+1) &
                     - 8*Ay(i-2,j,k-1) + 64*Ay(i-1,j,k-1) - 64*Ay(i+1,j,k-1) + 8*Ay(i+2,j,k-1) &
                     +   Ay(i-2,j,k-2) -  8*Ay(i-1,j,k-2) +  8*Ay(i+1,j,k-2) -   Ay(i+2,j,k-2) ) / dxdz144
      d2_lA(3,1,3) = (  -Az(i-2,j,k+2) +  8*Az(i-1,j,k+2) -  8*Az(i+1,j,k+2) +   Az(i+2,j,k+2) &
                     + 8*Az(i-2,j,k+1) - 64*Az(i-1,j,k+1) + 64*Az(i+1,j,k+1) - 8*Az(i+2,j,k+1) &
                     - 8*Az(i-2,j,k-1) + 64*Az(i-1,j,k-1) - 64*Az(i+1,j,k-1) + 8*Az(i+2,j,k-1) &
                     +   Az(i-2,j,k-2) -  8*Az(i-1,j,k-2) +  8*Az(i+1,j,k-2) -   Az(i+2,j,k-2) ) / dxdz144

      d2_lA(1,2,3) = (  -Ax(i,j-2,k+2) +  8*Ax(i,j-1,k+2) -  8*Ax(i,j+1,k+2) +   Ax(i,j+2,k+2) &
                     + 8*Ax(i,j-2,k+1) - 64*Ax(i,j-1,k+1) + 64*Ax(i,j+1,k+1) - 8*Ax(i,j+2,k+1) &
                     - 8*Ax(i,j-2,k-1) + 64*Ax(i,j-1,k-1) - 64*Ax(i,j+1,k-1) + 8*Ax(i,j+2,k-1) &
                     +   Ax(i,j-2,k-2) -  8*Ax(i,j-1,k-2) +  8*Ax(i,j+1,k-2) -   Ax(i,j+2,k-2) ) / dydz144
      d2_lA(2,2,3) = (  -Ay(i,j-2,k+2) +  8*Ay(i,j-1,k+2) -  8*Ay(i,j+1,k+2) +   Ay(i,j+2,k+2) &
                     + 8*Ay(i,j-2,k+1) - 64*Ay(i,j-1,k+1) + 64*Ay(i,j+1,k+1) - 8*Ay(i,j+2,k+1) &
                     - 8*Ay(i,j-2,k-1) + 64*Ay(i,j-1,k-1) - 64*Ay(i,j+1,k-1) + 8*Ay(i,j+2,k-1) &
                     +   Ay(i,j-2,k-2) -  8*Ay(i,j-1,k-2) +  8*Ay(i,j+1,k-2) -   Ay(i,j+2,k-2) ) / dydz144
      d2_lA(3,2,3) = (  -Az(i,j-2,k+2) +  8*Az(i,j-1,k+2) -  8*Az(i,j+1,k+2) +   Az(i,j+2,k+2) &
                     + 8*Az(i,j-2,k+1) - 64*Az(i,j-1,k+1) + 64*Az(i,j+1,k+1) - 8*Az(i,j+2,k+1) &
                     - 8*Az(i,j-2,k-1) + 64*Az(i,j-1,k-1) - 64*Az(i,j+1,k-1) + 8*Az(i,j+2,k-1) &
                     +   Az(i,j-2,k-2) -  8*Az(i,j-1,k-2) +  8*Az(i,j+1,k-2) -   Az(i,j+2,k-2) ) / dydz144

      d2_lA(:,2,1) = d2_lA(:,1,2)
      d2_lA(:,3,1) = d2_lA(:,1,3)
      d2_lA(:,3,2) = d2_lA(:,2,3)

      ! d2_lphi1(3,3)
      d2_lphi1(1,1) = (   -phi1(i+2,j,k) + 16*phi1(i+1,j,k) - 30*phi1(i,j,k)     &
                      + 16*phi1(i-1,j,k) -    phi1(i-2,j,k) ) / dxsq12

      d2_lphi1(2,2) = (   -phi1(i,j+2,k) + 16*phi1(i,j+1,k) - 30*phi1(i,j,k)     &
                      + 16*phi1(i,j-1,k) -    phi1(i,j-2,k) ) / dysq12

      d2_lphi1(3,3) = (   -phi1(i,j,k+2) + 16*phi1(i,j,k+1) - 30*phi1(i,j,k)     &
                      + 16*phi1(i,j,k-1) -    phi1(i,j,k-2) ) / dzsq12

      d2_lphi1(1,2) = (   -phi1(i-2,j+2,k) +  8*phi1(i-1,j+2,k) -  8*phi1(i+1,j+2,k) +   phi1(i+2,j+2,k)   &
                      + 8*phi1(i-2,j+1,k) - 64*phi1(i-1,j+1,k) + 64*phi1(i+1,j+1,k) - 8*phi1(i+2,j+1,k)   &
                      - 8*phi1(i-2,j-1,k) + 64*phi1(i-1,j-1,k) - 64*phi1(i+1,j-1,k) + 8*phi1(i+2,j-1,k)   &
                      +   phi1(i-2,j-2,k) -  8*phi1(i-1,j-2,k) +  8*phi1(i+1,j-2,k) -   phi1(i+2,j-2,k) ) / dxdy144

      d2_lphi1(1,3) = (   -phi1(i-2,j,k+2) +  8*phi1(i-1,j,k+2) -  8*phi1(i+1,j,k+2) +   phi1(i+2,j,k+2)   &
                      + 8*phi1(i-2,j,k+1) - 64*phi1(i-1,j,k+1) + 64*phi1(i+1,j,k+1) - 8*phi1(i+2,j,k+1)   &
                      - 8*phi1(i-2,j,k-1) + 64*phi1(i-1,j,k-1) - 64*phi1(i+1,j,k-1) + 8*phi1(i+2,j,k-1)   &
                      +   phi1(i-2,j,k-2) -  8*phi1(i-1,j,k-2) +  8*phi1(i+1,j,k-2) -   phi1(i+2,j,k-2) ) / dxdz144

      d2_lphi1(2,3) = (   -phi1(i,j-2,k+2) +  8*phi1(i,j-1,k+2) -  8*phi1(i,j+1,k+2) +   phi1(i,j+2,k+2)   &
                      + 8*phi1(i,j-2,k+1) - 64*phi1(i,j-1,k+1) + 64*phi1(i,j+1,k+1) - 8*phi1(i,j+2,k+1)   &
                      - 8*phi1(i,j-2,k-1) + 64*phi1(i,j-1,k-1) - 64*phi1(i,j+1,k-1) + 8*phi1(i,j+2,k-1)   &
                      +   phi1(i,j-2,k-2) -  8*phi1(i,j-1,k-2) +  8*phi1(i,j+1,k-2) -   phi1(i,j+2,k-2) ) / dydz144

      d2_lphi1(2,1) = d2_lphi1(1,2)
      d2_lphi1(3,1) = d2_lphi1(1,3)
      d2_lphi1(3,2) = d2_lphi1(2,3)

      ! d2_lphi2(3,3)
      d2_lphi2(1,1) = (   -phi2(i+2,j,k) + 16*phi2(i+1,j,k) - 30*phi2(i,j,k)     &
                      + 16*phi2(i-1,j,k) -    phi2(i-2,j,k) ) / dxsq12

      d2_lphi2(2,2) = (   -phi2(i,j+2,k) + 16*phi2(i,j+1,k) - 30*phi2(i,j,k)     &
                      + 16*phi2(i,j-1,k) -    phi2(i,j-2,k) ) / dysq12

      d2_lphi2(3,3) = (   -phi2(i,j,k+2) + 16*phi2(i,j,k+1) - 30*phi2(i,j,k)     &
                      + 16*phi2(i,j,k-1) -    phi2(i,j,k-2) ) / dzsq12

      d2_lphi2(1,2) = (   -phi2(i-2,j+2,k) +  8*phi2(i-1,j+2,k) -  8*phi2(i+1,j+2,k) +   phi2(i+2,j+2,k)   &
                      + 8*phi2(i-2,j+1,k) - 64*phi2(i-1,j+1,k) + 64*phi2(i+1,j+1,k) - 8*phi2(i+2,j+1,k)   &
                      - 8*phi2(i-2,j-1,k) + 64*phi2(i-1,j-1,k) - 64*phi2(i+1,j-1,k) + 8*phi2(i+2,j-1,k)   &
                      +   phi2(i-2,j-2,k) -  8*phi2(i-1,j-2,k) +  8*phi2(i+1,j-2,k) -   phi2(i+2,j-2,k) ) / dxdy144

      d2_lphi2(1,3) = (   -phi2(i-2,j,k+2) +  8*phi2(i-1,j,k+2) -  8*phi2(i+1,j,k+2) +   phi2(i+2,j,k+2)   &
                      + 8*phi2(i-2,j,k+1) - 64*phi2(i-1,j,k+1) + 64*phi2(i+1,j,k+1) - 8*phi2(i+2,j,k+1)   &
                      - 8*phi2(i-2,j,k-1) + 64*phi2(i-1,j,k-1) - 64*phi2(i+1,j,k-1) + 8*phi2(i+2,j,k-1)   &
                      +   phi2(i-2,j,k-2) -  8*phi2(i-1,j,k-2) +  8*phi2(i+1,j,k-2) -   phi2(i+2,j,k-2) ) / dxdz144

      d2_lphi2(2,3) = (   -phi2(i,j-2,k+2) +  8*phi2(i,j-1,k+2) -  8*phi2(i,j+1,k+2) +   phi2(i,j+2,k+2)   &
                      + 8*phi2(i,j-2,k+1) - 64*phi2(i,j-1,k+1) + 64*phi2(i,j+1,k+1) - 8*phi2(i,j+2,k+1)   &
                      - 8*phi2(i,j-2,k-1) + 64*phi2(i,j-1,k-1) - 64*phi2(i,j+1,k-1) + 8*phi2(i,j+2,k-1)   &
                      +   phi2(i,j-2,k-2) -  8*phi2(i,j-1,k-2) +  8*phi2(i,j+1,k-2) -   phi2(i,j+2,k-2) ) / dydz144

      d2_lphi2(2,1) = d2_lphi2(1,2)
      d2_lphi2(3,1) = d2_lphi2(1,3)
      d2_lphi2(3,2) = d2_lphi2(2,3)

      !------------ Advection derivatives --------
      if( use_advection_stencils /= 0 ) then

        di = int( sign( one, beta(1) ) )
        dj = int( sign( one, beta(2) ) )
        dk = int( sign( one, beta(3) ) )

        ! ad1_lE(3)
        d1_f(1) = di * ( -3*Ex(i-di,j,k) - 10*Ex(i,j,k) + 18*Ex(i+di,j,k) &
                        - 6*Ex(i+2*di,j,k) + Ex(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3*Ex(i,j-dj,k) - 10*Ex(i,j,k) + 18*Ex(i,j+dj,k) &
                        - 6*Ex(i,j+2*dj,k) + Ex(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3*Ex(i,j,k-dk) - 10*Ex(i,j,k) + 18*Ex(i,j,k+dk) &
                        - 6*Ex(i,j,k+2*dk) + Ex(i,j,k+3*dk)) / dz12
        ad1_lE(1) = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        d1_f(1) = di * ( -3*Ey(i-di,j,k) - 10*Ey(i,j,k) + 18*Ey(i+di,j,k) &
                        - 6*Ey(i+2*di,j,k) + Ey(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3*Ey(i,j-dj,k) - 10*Ey(i,j,k) + 18*Ey(i,j+dj,k) &
                        - 6*Ey(i,j+2*dj,k) + Ey(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3*Ey(i,j,k-dk) - 10*Ey(i,j,k) + 18*Ey(i,j,k+dk) &
                        - 6*Ey(i,j,k+2*dk) + Ey(i,j,k+3*dk)) / dz12
        ad1_lE(2) = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        d1_f(1) = di * ( -3*Ez(i-di,j,k) - 10*Ez(i,j,k) + 18*Ez(i+di,j,k) &
                        - 6*Ez(i+2*di,j,k) + Ez(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3*Ez(i,j-dj,k) - 10*Ez(i,j,k) + 18*Ez(i,j+dj,k) &
                        - 6*Ez(i,j+2*dj,k) + Ez(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3*Ez(i,j,k-dk) - 10*Ez(i,j,k) + 18*Ez(i,j,k+dk) &
                        - 6*Ez(i,j,k+2*dk) + Ez(i,j,k+3*dk)) / dz12
        ad1_lE(3) = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        ! ad1_lA(3)
        d1_f(1) = di * ( -3*Ax(i-di,j,k) - 10*Ax(i,j,k) + 18*Ax(i+di,j,k) &
                        - 6*Ax(i+2*di,j,k) + Ax(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3*Ax(i,j-dj,k) - 10*Ax(i,j,k) + 18*Ax(i,j+dj,k) &
                        - 6*Ax(i,j+2*dj,k) + Ax(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3*Ax(i,j,k-dk) - 10*Ax(i,j,k) + 18*Ax(i,j,k+dk) &
                        - 6*Ax(i,j,k+2*dk) + Ax(i,j,k+3*dk)) / dz12
        ad1_lA(1) = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        d1_f(1) = di * ( -3*Ay(i-di,j,k) - 10*Ay(i,j,k) + 18*Ay(i+di,j,k) &
                        - 6*Ay(i+2*di,j,k) + Ay(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3*Ay(i,j-dj,k) - 10*Ay(i,j,k) + 18*Ay(i,j+dj,k) &
                        - 6*Ay(i,j+2*dj,k) + Ay(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3*Ay(i,j,k-dk) - 10*Ay(i,j,k) + 18*Ay(i,j,k+dk) &
                        - 6*Ay(i,j,k+2*dk) + Ay(i,j,k+3*dk)) / dz12
        ad1_lA(2) = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        d1_f(1) = di * ( -3*Az(i-di,j,k) - 10*Az(i,j,k) + 18*Az(i+di,j,k) &
                        - 6*Az(i+2*di,j,k) + Az(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3*Az(i,j-dj,k) - 10*Az(i,j,k) + 18*Az(i,j+dj,k) &
                        - 6*Az(i,j+2*dj,k) + Az(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3*Az(i,j,k-dk) - 10*Az(i,j,k) + 18*Az(i,j,k+dk) &
                        - 6*Az(i,j,k+2*dk) + Az(i,j,k+3*dk)) / dz12
        ad1_lA(3) = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        ! ad1_lZeta
        d1_f(1) = di * ( -3*Zeta(i-di,j,k) - 10*Zeta(i,j,k) + 18*Zeta(i+di,j,k)  &
                        - 6*Zeta(i+2*di,j,k) + Zeta(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3*Zeta(i,j-dj,k) - 10*Zeta(i,j,k) + 18*Zeta(i,j+dj,k)  &
                        - 6*Zeta(i,j+2*dj,k) + Zeta(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3*Zeta(i,j,k-dk) - 10*Zeta(i,j,k) + 18*Zeta(i,j,k+dk)  &
                        - 6*Zeta(i,j,k+2*dk) + Zeta(i,j,k+3*dk)) / dz12
        ad1_lZeta = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        ! ad1_lAphi
        d1_f(1) = di * ( -3*Aphi(i-di,j,k) - 10*Aphi(i,j,k) + 18*Aphi(i+di,j,k)  &
                        - 6*Aphi(i+2*di,j,k) + Aphi(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3*Aphi(i,j-dj,k) - 10*Aphi(i,j,k) + 18*Aphi(i,j+dj,k)  &
                        - 6*Aphi(i,j+2*dj,k) + Aphi(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3*Aphi(i,j,k-dk) - 10*Aphi(i,j,k) + 18*Aphi(i,j,k+dk)  &
                        - 6*Aphi(i,j,k+2*dk) + Aphi(i,j,k+3*dk)) / dz12
        ad1_lAphi = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        ! ad1_lphi1
        d1_f(1) = di * ( -3*phi1(i-di,j,k) - 10*phi1(i,j,k) + 18*phi1(i+di,j,k)   &
                        - 6*phi1(i+2*di,j,k)  + phi1(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3*phi1(i,j-dj,k) - 10*phi1(i,j,k) + 18*phi1(i,j+dj,k)   &
                        - 6*phi1(i,j+2*dj,k)  + phi1(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3*phi1(i,j,k-dk) - 10*phi1(i,j,k) + 18*phi1(i,j,k+dk)   &
                        - 6*phi1(i,j,k+2*dk)  + phi1(i,j,k+3*dk)) / dz12
        ad1_lphi1 = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        ! ad1_lKphi1
        d1_f(1) = di * ( -3*Kphi1(i-di,j,k) - 10*Kphi1(i,j,k) + 18*Kphi1(i+di,j,k)   &
                        - 6*Kphi1(i+2*di,j,k)  + Kphi1(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3*Kphi1(i,j-dj,k) - 10*Kphi1(i,j,k) + 18*Kphi1(i,j+dj,k)   &
                        - 6*Kphi1(i,j+2*dj,k)  + Kphi1(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3*Kphi1(i,j,k-dk) - 10*Kphi1(i,j,k) + 18*Kphi1(i,j,k+dk)   &
                        - 6*Kphi1(i,j,k+2*dk)  + Kphi1(i,j,k+3*dk)) / dz12
        ad1_lKphi1 = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        ! ad1_lphi2
        d1_f(1) = di * ( -3*phi2(i-di,j,k) - 10*phi2(i,j,k) + 18*phi2(i+di,j,k)   &
                        - 6*phi2(i+2*di,j,k)  + phi2(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3*phi2(i,j-dj,k) - 10*phi2(i,j,k) + 18*phi2(i,j+dj,k)   &
                        - 6*phi2(i,j+2*dj,k)  + phi2(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3*phi2(i,j,k-dk) - 10*phi2(i,j,k) + 18*phi2(i,j,k+dk)   &
                        - 6*phi2(i,j,k+2*dk)  + phi2(i,j,k+3*dk)) / dz12
        ad1_lphi2 = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        ! ad1_lKphi2
        d1_f(1) = di * ( -3*Kphi2(i-di,j,k) - 10*Kphi2(i,j,k) + 18*Kphi2(i+di,j,k)   &
                        - 6*Kphi2(i+2*di,j,k)  + Kphi2(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3*Kphi2(i,j-dj,k) - 10*Kphi2(i,j,k) + 18*Kphi2(i,j+dj,k)   &
                        - 6*Kphi2(i,j+2*dj,k)  + Kphi2(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3*Kphi2(i,j,k-dk) - 10*Kphi2(i,j,k) + 18*Kphi2(i,j,k+dk)   &
                        - 6*Kphi2(i,j,k+2*dk)  + Kphi2(i,j,k+3*dk)) / dz12
        ad1_lKphi2 = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

      else

        ! ad1_lE(3)
        ad1_lE(1) = beta(1)*d1_lE(1,1) + beta(2)*d1_lE(1,2) + beta(3)*d1_lE(1,3)
        ad1_lE(2) = beta(1)*d1_lE(2,1) + beta(2)*d1_lE(2,2) + beta(3)*d1_lE(2,3)
        ad1_lE(3) = beta(1)*d1_lE(3,1) + beta(2)*d1_lE(3,2) + beta(3)*d1_lE(3,3)

        ! ad1_lA(3)
        ad1_lA(1) = beta(1)*d1_lA(1,1) + beta(2)*d1_lA(1,2) + beta(3)*d1_lA(1,3)
        ad1_lA(2) = beta(1)*d1_lA(2,1) + beta(2)*d1_lA(2,2) + beta(3)*d1_lA(2,3)
        ad1_lA(3) = beta(1)*d1_lA(3,1) + beta(2)*d1_lA(3,2) + beta(3)*d1_lA(3,3)

        ! ad1_lZeta
        ad1_lZeta = beta(1)*d1_lZeta(1) + beta(2)*d1_lZeta(2) + beta(3)*d1_lZeta(3)

        ! ad1_lAphi
        ad1_lAphi = beta(1)*d1_lAphi(1) + beta(2)*d1_lAphi(2) + beta(3)*d1_lAphi(3)

        ! ad1_lphi1
        ad1_lphi1 = beta(1)*d1_lphi1(1) + beta(2)*d1_lphi1(2) + beta(3)*d1_lphi1(3)

        ! ad1_lKphi1
        ad1_lKphi1 = beta(1)*d1_lKphi1(1) + beta(2)*d1_lKphi1(2) + beta(3)*d1_lKphi1(3)

        ! ad1_lphi2
        ad1_lphi2 = beta(1)*d1_lphi2(1) + beta(2)*d1_lphi2(2) + beta(3)*d1_lphi2(3)

        ! ad1_lKphi2
        ad1_lKphi2 = beta(1)*d1_lKphi2(1) + beta(2)*d1_lKphi2(2) + beta(3)*d1_lKphi2(3)

      end if
      !-------------------------------------------

    else if (derivs_order == 6) then

      !------------ Centered 1st derivatives -----
      ! d1_ch(3)
      d1_ch(1) = (  chi(i+3,j,k) - 9*chi(i+2,j,k) + 45*chi(i+1,j,k)          &
                  - chi(i-3,j,k) + 9*chi(i-2,j,k) - 45*chi(i-1,j,k) ) * odx60

      d1_ch(2) = (  chi(i,j+3,k) - 9*chi(i,j+2,k) + 45*chi(i,j+1,k)          &
                  - chi(i,j-3,k) + 9*chi(i,j-2,k) - 45*chi(i,j-1,k) ) * ody60

      d1_ch(3) = (  chi(i,j,k+3) - 9*chi(i,j,k+2) + 45*chi(i,j,k+1)          &
                  - chi(i,j,k-3) + 9*chi(i,j,k-2) - 45*chi(i,j,k-1) ) * odz60

      ! d1_hh(3,3,3)
      d1_hh(1,1,1) = (  hxx(i+3,j,k) - 9*hxx(i+2,j,k) + 45*hxx(i+1,j,k)      &
                      - hxx(i-3,j,k) + 9*hxx(i-2,j,k) - 45*hxx(i-1,j,k) ) * odx60
      d1_hh(1,2,1) = (  hxy(i+3,j,k) - 9*hxy(i+2,j,k) + 45*hxy(i+1,j,k)      &
                      - hxy(i-3,j,k) + 9*hxy(i-2,j,k) - 45*hxy(i-1,j,k) ) * odx60
      d1_hh(1,3,1) = (  hxz(i+3,j,k) - 9*hxz(i+2,j,k) + 45*hxz(i+1,j,k)      &
                      - hxz(i-3,j,k) + 9*hxz(i-2,j,k) - 45*hxz(i-1,j,k) ) * odx60
      d1_hh(2,2,1) = (  hyy(i+3,j,k) - 9*hyy(i+2,j,k) + 45*hyy(i+1,j,k)      &
                      - hyy(i-3,j,k) + 9*hyy(i-2,j,k) - 45*hyy(i-1,j,k) ) * odx60
      d1_hh(2,3,1) = (  hyz(i+3,j,k) - 9*hyz(i+2,j,k) + 45*hyz(i+1,j,k)      &
                      - hyz(i-3,j,k) + 9*hyz(i-2,j,k) - 45*hyz(i-1,j,k) ) * odx60
      d1_hh(3,3,1) = (  hzz(i+3,j,k) - 9*hzz(i+2,j,k) + 45*hzz(i+1,j,k)      &
                      - hzz(i-3,j,k) + 9*hzz(i-2,j,k) - 45*hzz(i-1,j,k) ) * odx60

      d1_hh(1,1,2) = (  hxx(i,j+3,k) - 9*hxx(i,j+2,k) + 45*hxx(i,j+1,k)      &
                      - hxx(i,j-3,k) + 9*hxx(i,j-2,k) - 45*hxx(i,j-1,k) ) * ody60
      d1_hh(1,2,2) = (  hxy(i,j+3,k) - 9*hxy(i,j+2,k) + 45*hxy(i,j+1,k)      &
                      - hxy(i,j-3,k) + 9*hxy(i,j-2,k) - 45*hxy(i,j-1,k) ) * ody60
      d1_hh(1,3,2) = (  hxz(i,j+3,k) - 9*hxz(i,j+2,k) + 45*hxz(i,j+1,k)      &
                      - hxz(i,j-3,k) + 9*hxz(i,j-2,k) - 45*hxz(i,j-1,k) ) * ody60
      d1_hh(2,2,2) = (  hyy(i,j+3,k) - 9*hyy(i,j+2,k) + 45*hyy(i,j+1,k)      &
                      - hyy(i,j-3,k) + 9*hyy(i,j-2,k) - 45*hyy(i,j-1,k) ) * ody60
      d1_hh(2,3,2) = (  hyz(i,j+3,k) - 9*hyz(i,j+2,k) + 45*hyz(i,j+1,k)      &
                      - hyz(i,j-3,k) + 9*hyz(i,j-2,k) - 45*hyz(i,j-1,k) ) * ody60
      d1_hh(3,3,2) = (  hzz(i,j+3,k) - 9*hzz(i,j+2,k) + 45*hzz(i,j+1,k)      &
                      - hzz(i,j-3,k) + 9*hzz(i,j-2,k) - 45*hzz(i,j-1,k) ) * ody60

      d1_hh(1,1,3) = (  hxx(i,j,k+3) - 9*hxx(i,j,k+2) + 45*hxx(i,j,k+1)      &
                      - hxx(i,j,k-3) + 9*hxx(i,j,k-2) - 45*hxx(i,j,k-1) ) * odz60
      d1_hh(1,2,3) = (  hxy(i,j,k+3) - 9*hxy(i,j,k+2) + 45*hxy(i,j,k+1)      &
                      - hxy(i,j,k-3) + 9*hxy(i,j,k-2) - 45*hxy(i,j,k-1) ) * odz60
      d1_hh(1,3,3) = (  hxz(i,j,k+3) - 9*hxz(i,j,k+2) + 45*hxz(i,j,k+1)      &
                      - hxz(i,j,k-3) + 9*hxz(i,j,k-2) - 45*hxz(i,j,k-1) ) * odz60
      d1_hh(2,2,3) = (  hyy(i,j,k+3) - 9*hyy(i,j,k+2) + 45*hyy(i,j,k+1)      &
                      - hyy(i,j,k-3) + 9*hyy(i,j,k-2) - 45*hyy(i,j,k-1) ) * odz60
      d1_hh(2,3,3) = (  hyz(i,j,k+3) - 9*hyz(i,j,k+2) + 45*hyz(i,j,k+1)      &
                      - hyz(i,j,k-3) + 9*hyz(i,j,k-2) - 45*hyz(i,j,k-1) ) * odz60
      d1_hh(3,3,3) = (  hzz(i,j,k+3) - 9*hzz(i,j,k+2) + 45*hzz(i,j,k+1)      &
                      - hzz(i,j,k-3) + 9*hzz(i,j,k-2) - 45*hzz(i,j,k-1) ) * odz60

      d1_hh(2,1,:) = d1_hh(1,2,:)
      d1_hh(3,1,:) = d1_hh(1,3,:)
      d1_hh(3,2,:) = d1_hh(2,3,:)

      ! d1_alph(3)
      d1_alph(1) = (  alp(i+3,j,k) - 9*alp(i+2,j,k) + 45*alp(i+1,j,k) &
                    - alp(i-3,j,k) + 9*alp(i-2,j,k) - 45*alp(i-1,j,k) ) * odx60

      d1_alph(2) = (  alp(i,j+3,k) - 9*alp(i,j+2,k) + 45*alp(i,j+1,k) &
                    - alp(i,j-3,k) + 9*alp(i,j-2,k) - 45*alp(i,j-1,k) ) * ody60

      d1_alph(3) = (  alp(i,j,k+3) - 9*alp(i,j,k+2) + 45*alp(i,j,k+1) &
                    - alp(i,j,k-3) + 9*alp(i,j,k-2) - 45*alp(i,j,k-1) ) * odz60

      ! d1_beta(3,3)
      d1_beta(1,1) = (  betax(i+3,j,k) - 9*betax(i+2,j,k) + 45*betax(i+1,j,k) &
                      - betax(i-3,j,k) + 9*betax(i-2,j,k) - 45*betax(i-1,j,k) ) * odx60
      d1_beta(2,1) = (  betay(i+3,j,k) - 9*betay(i+2,j,k) + 45*betay(i+1,j,k) &
                      - betay(i-3,j,k) + 9*betay(i-2,j,k) - 45*betay(i-1,j,k) ) * odx60
      d1_beta(3,1) = (  betaz(i+3,j,k) - 9*betaz(i+2,j,k) + 45*betaz(i+1,j,k) &
                      - betaz(i-3,j,k) + 9*betaz(i-2,j,k) - 45*betaz(i-1,j,k) ) * odx60

      d1_beta(1,2) = (  betax(i,j+3,k) - 9*betax(i,j+2,k) + 45*betax(i,j+1,k) &
                      - betax(i,j-3,k) + 9*betax(i,j-2,k) - 45*betax(i,j-1,k) ) * ody60
      d1_beta(2,2) = (  betay(i,j+3,k) - 9*betay(i,j+2,k) + 45*betay(i,j+1,k) &
                      - betay(i,j-3,k) + 9*betay(i,j-2,k) - 45*betay(i,j-1,k) ) * ody60
      d1_beta(3,2) = (  betaz(i,j+3,k) - 9*betaz(i,j+2,k) + 45*betaz(i,j+1,k) &
                      - betaz(i,j-3,k) + 9*betaz(i,j-2,k) - 45*betaz(i,j-1,k) ) * ody60

      d1_beta(1,3) = (  betax(i,j,k+3) - 9*betax(i,j,k+2) + 45*betax(i,j,k+1) &
                      - betax(i,j,k-3) + 9*betax(i,j,k-2) - 45*betax(i,j,k-1) ) * odz60
      d1_beta(2,3) = (  betay(i,j,k+3) - 9*betay(i,j,k+2) + 45*betay(i,j,k+1) &
                      - betay(i,j,k-3) + 9*betay(i,j,k-2) - 45*betay(i,j,k-1) ) * odz60
      d1_beta(3,3) = (  betaz(i,j,k+3) - 9*betaz(i,j,k+2) + 45*betaz(i,j,k+1) &
                      - betaz(i,j,k-3) + 9*betaz(i,j,k-2) - 45*betaz(i,j,k-1) ) * odz60

      ! d1_lE(3,3)
      d1_lE(1,1) = (  Ex(i+3,j,k) - 9*Ex(i+2,j,k) + 45*Ex(i+1,j,k) &
                    - Ex(i-3,j,k) + 9*Ex(i-2,j,k) - 45*Ex(i-1,j,k) ) * odx60
      d1_lE(2,1) = (  Ey(i+3,j,k) - 9*Ey(i+2,j,k) + 45*Ey(i+1,j,k) &
                    - Ey(i-3,j,k) + 9*Ey(i-2,j,k) - 45*Ey(i-1,j,k) ) * odx60
      d1_lE(3,1) = (  Ez(i+3,j,k) - 9*Ez(i+2,j,k) + 45*Ez(i+1,j,k) &
                    - Ez(i-3,j,k) + 9*Ez(i-2,j,k) - 45*Ez(i-1,j,k) ) * odx60

      d1_lE(1,2) = (  Ex(i,j+3,k) - 9*Ex(i,j+2,k) + 45*Ex(i,j+1,k) &
                    - Ex(i,j-3,k) + 9*Ex(i,j-2,k) - 45*Ex(i,j-1,k) ) * ody60
      d1_lE(2,2) = (  Ey(i,j+3,k) - 9*Ey(i,j+2,k) + 45*Ey(i,j+1,k) &
                    - Ey(i,j-3,k) + 9*Ey(i,j-2,k) - 45*Ey(i,j-1,k) ) * ody60
      d1_lE(3,2) = (  Ez(i,j+3,k) - 9*Ez(i,j+2,k) + 45*Ez(i,j+1,k) &
                    - Ez(i,j-3,k) + 9*Ez(i,j-2,k) - 45*Ez(i,j-1,k) ) * ody60

      d1_lE(1,3) = (  Ex(i,j,k+3) - 9*Ex(i,j,k+2) + 45*Ex(i,j,k+1) &
                    - Ex(i,j,k-3) + 9*Ex(i,j,k-2) - 45*Ex(i,j,k-1) ) * odz60
      d1_lE(2,3) = (  Ey(i,j,k+3) - 9*Ey(i,j,k+2) + 45*Ey(i,j,k+1) &
                    - Ey(i,j,k-3) + 9*Ey(i,j,k-2) - 45*Ey(i,j,k-1) ) * odz60
      d1_lE(3,3) = (  Ez(i,j,k+3) - 9*Ez(i,j,k+2) + 45*Ez(i,j,k+1) &
                    - Ez(i,j,k-3) + 9*Ez(i,j,k-2) - 45*Ez(i,j,k-1) ) * odz60

      ! d1_lA(3,3)
      d1_lA(1,1) = (  Ax(i+3,j,k) - 9*Ax(i+2,j,k) + 45*Ax(i+1,j,k) &
                    - Ax(i-3,j,k) + 9*Ax(i-2,j,k) - 45*Ax(i-1,j,k) ) * odx60
      d1_lA(2,1) = (  Ay(i+3,j,k) - 9*Ay(i+2,j,k) + 45*Ay(i+1,j,k) &
                    - Ay(i-3,j,k) + 9*Ay(i-2,j,k) - 45*Ay(i-1,j,k) ) * odx60
      d1_lA(3,1) = (  Az(i+3,j,k) - 9*Az(i+2,j,k) + 45*Az(i+1,j,k) &
                    - Az(i-3,j,k) + 9*Az(i-2,j,k) - 45*Az(i-1,j,k) ) * odx60

      d1_lA(1,2) = (  Ax(i,j+3,k) - 9*Ax(i,j+2,k) + 45*Ax(i,j+1,k) &
                    - Ax(i,j-3,k) + 9*Ax(i,j-2,k) - 45*Ax(i,j-1,k) ) * ody60
      d1_lA(2,2) = (  Ay(i,j+3,k) - 9*Ay(i,j+2,k) + 45*Ay(i,j+1,k) &
                    - Ay(i,j-3,k) + 9*Ay(i,j-2,k) - 45*Ay(i,j-1,k) ) * ody60
      d1_lA(3,2) = (  Az(i,j+3,k) - 9*Az(i,j+2,k) + 45*Az(i,j+1,k) &
                    - Az(i,j-3,k) + 9*Az(i,j-2,k) - 45*Az(i,j-1,k) ) * ody60

      d1_lA(1,3) = (  Ax(i,j,k+3) - 9*Ax(i,j,k+2) + 45*Ax(i,j,k+1) &
                    - Ax(i,j,k-3) + 9*Ax(i,j,k-2) - 45*Ax(i,j,k-1) ) * odz60
      d1_lA(2,3) = (  Ay(i,j,k+3) - 9*Ay(i,j,k+2) + 45*Ay(i,j,k+1) &
                    - Ay(i,j,k-3) + 9*Ay(i,j,k-2) - 45*Ay(i,j,k-1) ) * odz60
      d1_lA(3,3) = (  Az(i,j,k+3) - 9*Az(i,j,k+2) + 45*Az(i,j,k+1) &
                    - Az(i,j,k-3) + 9*Az(i,j,k-2) - 45*Az(i,j,k-1) ) * odz60

      ! d1_lZeta(3)
      d1_lZeta(1) = (  Zeta(i+3,j,k) - 9*Zeta(i+2,j,k) + 45*Zeta(i+1,j,k) &
                     - Zeta(i-3,j,k) + 9*Zeta(i-2,j,k) - 45*Zeta(i-1,j,k) ) * odx60

      d1_lZeta(2) = (  Zeta(i,j+3,k) - 9*Zeta(i,j+2,k) + 45*Zeta(i,j+1,k) &
                     - Zeta(i,j-3,k) + 9*Zeta(i,j-2,k) - 45*Zeta(i,j-1,k) ) * ody60

      d1_lZeta(3) = (  Zeta(i,j,k+3) - 9*Zeta(i,j,k+2) + 45*Zeta(i,j,k+1) &
                     - Zeta(i,j,k-3) + 9*Zeta(i,j,k-2) - 45*Zeta(i,j,k-1) ) * odz60

      ! d1_lAphi(3)
      d1_lAphi(1) = (  Aphi(i+3,j,k) - 9*Aphi(i+2,j,k) + 45*Aphi(i+1,j,k) &
                     - Aphi(i-3,j,k) + 9*Aphi(i-2,j,k) - 45*Aphi(i-1,j,k) ) * odx60

      d1_lAphi(2) = (  Aphi(i,j+3,k) - 9*Aphi(i,j+2,k) + 45*Aphi(i,j+1,k) &
                     - Aphi(i,j-3,k) + 9*Aphi(i,j-2,k) - 45*Aphi(i,j-1,k) ) * ody60

      d1_lAphi(3) = (  Aphi(i,j,k+3) - 9*Aphi(i,j,k+2) + 45*Aphi(i,j,k+1) &
                     - Aphi(i,j,k-3) + 9*Aphi(i,j,k-2) - 45*Aphi(i,j,k-1) ) * odz60

      !--------------------------------------------------


      !------------- Centered 2nd derivatives -----------

      ! d2_A(3,3,3)
      d2_lA(1,1,1) = (  2*Ax(i+3,j,k) - 27*Ax(i+2,j,k) + 270*Ax(i+1,j,k) - 490*Ax(i,j,k)&
                      + 2*Ax(i-3,j,k) - 27*Ax(i-2,j,k) + 270*Ax(i-1,j,k) ) * odxsq180
      d2_lA(2,1,1) = (  2*Ay(i+3,j,k) - 27*Ay(i+2,j,k) + 270*Ay(i+1,j,k) - 490*Ay(i,j,k)&
                      + 2*Ay(i-3,j,k) - 27*Ay(i-2,j,k) + 270*Ay(i-1,j,k) ) * odxsq180
      d2_lA(3,1,1) = (  2*Az(i+3,j,k) - 27*Az(i+2,j,k) + 270*Az(i+1,j,k) - 490*Az(i,j,k)&
                      + 2*Az(i-3,j,k) - 27*Az(i-2,j,k) + 270*Az(i-1,j,k) ) * odxsq180

      d2_lA(1,2,2) = (  2*Ax(i,j+3,k) - 27*Ax(i,j+2,k) + 270*Ax(i,j+1,k) - 490*Ax(i,j,k)&
                      + 2*Ax(i,j-3,k) - 27*Ax(i,j-2,k) + 270*Ax(i,j-1,k) ) * odysq180
      d2_lA(2,2,2) = (  2*Ay(i,j+3,k) - 27*Ay(i,j+2,k) + 270*Ay(i,j+1,k) - 490*Ay(i,j,k)&
                      + 2*Ay(i,j-3,k) - 27*Ay(i,j-2,k) + 270*Ay(i,j-1,k) ) * odysq180
      d2_lA(3,2,2) = (  2*Az(i,j+3,k) - 27*Az(i,j+2,k) + 270*Az(i,j+1,k) - 490*Az(i,j,k)&
                      + 2*Az(i,j-3,k) - 27*Az(i,j-2,k) + 270*Az(i,j-1,k) ) * odysq180

      d2_lA(1,3,3) = (  2*Ax(i,j,k+3) - 27*Ax(i,j,k+2) + 270*Ax(i,j,k+1) - 490*Ax(i,j,k)&
                      + 2*Ax(i,j,k-3) - 27*Ax(i,j,k-2) + 270*Ax(i,j,k-1) ) * odzsq180
      d2_lA(2,3,3) = (  2*Ay(i,j,k+3) - 27*Ay(i,j,k+2) + 270*Ay(i,j,k+1) - 490*Ay(i,j,k)&
                      + 2*Ay(i,j,k-3) - 27*Ay(i,j,k-2) + 270*Ay(i,j,k-1) ) * odzsq180
      d2_lA(3,3,3) = (  2*Az(i,j,k+3) - 27*Az(i,j,k+2) + 270*Az(i,j,k+1) - 490*Az(i,j,k)&
                      + 2*Az(i,j,k-3) - 27*Az(i,j,k-2) + 270*Az(i,j,k-1) ) * odzsq180

      d2_lA(1,1,2) = (    -Ax(i-3,j+3,k) +   9*Ax(i-2,j+3,k) -   45*Ax(i-1,j+3,k) +   45*Ax(i+1,j+3,k) -   9*Ax(i+2,j+3,k) +    Ax(i+3,j+3,k) &
                      +  9*Ax(i-3,j+2,k) -  81*Ax(i-2,j+2,k) +  405*Ax(i-1,j+2,k) -  405*Ax(i+1,j+2,k) +  81*Ax(i+2,j+2,k) -  9*Ax(i+3,j+2,k) &
                      - 45*Ax(i-3,j+1,k) + 405*Ax(i-2,j+1,k) - 2025*Ax(i-1,j+1,k) + 2025*Ax(i+1,j+1,k) - 405*Ax(i+2,j+1,k) + 45*Ax(i+3,j+1,k) &
                      + 45*Ax(i-3,j-1,k) - 405*Ax(i-2,j-1,k) + 2025*Ax(i-1,j-1,k) - 2025*Ax(i+1,j-1,k) + 405*Ax(i+2,j-1,k) - 45*Ax(i+3,j-1,k) &
                      -  9*Ax(i-3,j-2,k) +  81*Ax(i-2,j-2,k) -  405*Ax(i-1,j-2,k) +  405*Ax(i+1,j-2,k) -  81*Ax(i+2,j-2,k) +  9*Ax(i+3,j-2,k) &
                      +    Ax(i-3,j-3,k) -   9*Ax(i-2,j-3,k) +   45*Ax(i-1,j-3,k) -   45*Ax(i+1,j-3,k) +   9*Ax(i+2,j-3,k) -    Ax(i+3,j-3,k) ) * odxdy3600
      d2_lA(2,1,2) = (    -Ay(i-3,j+3,k) +   9*Ay(i-2,j+3,k) -   45*Ay(i-1,j+3,k) +   45*Ay(i+1,j+3,k) -   9*Ay(i+2,j+3,k) +    Ay(i+3,j+3,k) &
                      +  9*Ay(i-3,j+2,k) -  81*Ay(i-2,j+2,k) +  405*Ay(i-1,j+2,k) -  405*Ay(i+1,j+2,k) +  81*Ay(i+2,j+2,k) -  9*Ay(i+3,j+2,k) &
                      - 45*Ay(i-3,j+1,k) + 405*Ay(i-2,j+1,k) - 2025*Ay(i-1,j+1,k) + 2025*Ay(i+1,j+1,k) - 405*Ay(i+2,j+1,k) + 45*Ay(i+3,j+1,k) &
                      + 45*Ay(i-3,j-1,k) - 405*Ay(i-2,j-1,k) + 2025*Ay(i-1,j-1,k) - 2025*Ay(i+1,j-1,k) + 405*Ay(i+2,j-1,k) - 45*Ay(i+3,j-1,k) &
                      -  9*Ay(i-3,j-2,k) +  81*Ay(i-2,j-2,k) -  405*Ay(i-1,j-2,k) +  405*Ay(i+1,j-2,k) -  81*Ay(i+2,j-2,k) +  9*Ay(i+3,j-2,k) &
                      +    Ay(i-3,j-3,k) -   9*Ay(i-2,j-3,k) +   45*Ay(i-1,j-3,k) -   45*Ay(i+1,j-3,k) +   9*Ay(i+2,j-3,k) -    Ay(i+3,j-3,k) ) * odxdy3600
      d2_lA(3,1,2) = (    -Az(i-3,j+3,k) +   9*Az(i-2,j+3,k) -   45*Az(i-1,j+3,k) +   45*Az(i+1,j+3,k) -   9*Az(i+2,j+3,k) +    Az(i+3,j+3,k) &
                      +  9*Az(i-3,j+2,k) -  81*Az(i-2,j+2,k) +  405*Az(i-1,j+2,k) -  405*Az(i+1,j+2,k) +  81*Az(i+2,j+2,k) -  9*Az(i+3,j+2,k) &
                      - 45*Az(i-3,j+1,k) + 405*Az(i-2,j+1,k) - 2025*Az(i-1,j+1,k) + 2025*Az(i+1,j+1,k) - 405*Az(i+2,j+1,k) + 45*Az(i+3,j+1,k) &
                      + 45*Az(i-3,j-1,k) - 405*Az(i-2,j-1,k) + 2025*Az(i-1,j-1,k) - 2025*Az(i+1,j-1,k) + 405*Az(i+2,j-1,k) - 45*Az(i+3,j-1,k) &
                      -  9*Az(i-3,j-2,k) +  81*Az(i-2,j-2,k) -  405*Az(i-1,j-2,k) +  405*Az(i+1,j-2,k) -  81*Az(i+2,j-2,k) +  9*Az(i+3,j-2,k) &
                      +    Az(i-3,j-3,k) -   9*Az(i-2,j-3,k) +   45*Az(i-1,j-3,k) -   45*Az(i+1,j-3,k) +   9*Az(i+2,j-3,k) -    Az(i+3,j-3,k) ) * odxdy3600


      d2_lA(1,1,3) = (    -Ax(i-3,j,k+3) +   9*Ax(i-2,j,k+3) -   45*Ax(i-1,j,k+3) +   45*Ax(i+1,j,k+3) -   9*Ax(i+2,j,k+3) +    Ax(i+3,j,k+3) &
                      +  9*Ax(i-3,j,k+2) -  81*Ax(i-2,j,k+2) +  405*Ax(i-1,j,k+2) -  405*Ax(i+1,j,k+2) +  81*Ax(i+2,j,k+2) -  9*Ax(i+3,j,k+2) &
                      - 45*Ax(i-3,j,k+1) + 405*Ax(i-2,j,k+1) - 2025*Ax(i-1,j,k+1) + 2025*Ax(i+1,j,k+1) - 405*Ax(i+2,j,k+1) + 45*Ax(i+3,j,k+1) &
                      + 45*Ax(i-3,j,k-1) - 405*Ax(i-2,j,k-1) + 2025*Ax(i-1,j,k-1) - 2025*Ax(i+1,j,k-1) + 405*Ax(i+2,j,k-1) - 45*Ax(i+3,j,k-1) &
                      -  9*Ax(i-3,j,k-2) +  81*Ax(i-2,j,k-2) -  405*Ax(i-1,j,k-2) +  405*Ax(i+1,j,k-2) -  81*Ax(i+2,j,k-2) +  9*Ax(i+3,j,k-2) &
                      +    Ax(i-3,j,k-3) -   9*Ax(i-2,j,k-3) +   45*Ax(i-1,j,k-3) -   45*Ax(i+1,j,k-3) +   9*Ax(i+2,j,k-3) -    Ax(i+3,j,k-3) ) * odxdz3600
      d2_lA(2,1,3) = (    -Ay(i-3,j,k+3) +   9*Ay(i-2,j,k+3) -   45*Ay(i-1,j,k+3) +   45*Ay(i+1,j,k+3) -   9*Ay(i+2,j,k+3) +    Ay(i+3,j,k+3) &
                      +  9*Ay(i-3,j,k+2) -  81*Ay(i-2,j,k+2) +  405*Ay(i-1,j,k+2) -  405*Ay(i+1,j,k+2) +  81*Ay(i+2,j,k+2) -  9*Ay(i+3,j,k+2) &
                      - 45*Ay(i-3,j,k+1) + 405*Ay(i-2,j,k+1) - 2025*Ay(i-1,j,k+1) + 2025*Ay(i+1,j,k+1) - 405*Ay(i+2,j,k+1) + 45*Ay(i+3,j,k+1) &
                      + 45*Ay(i-3,j,k-1) - 405*Ay(i-2,j,k-1) + 2025*Ay(i-1,j,k-1) - 2025*Ay(i+1,j,k-1) + 405*Ay(i+2,j,k-1) - 45*Ay(i+3,j,k-1) &
                      -  9*Ay(i-3,j,k-2) +  81*Ay(i-2,j,k-2) -  405*Ay(i-1,j,k-2) +  405*Ay(i+1,j,k-2) -  81*Ay(i+2,j,k-2) +  9*Ay(i+3,j,k-2) &
                      +    Ay(i-3,j,k-3) -   9*Ay(i-2,j,k-3) +   45*Ay(i-1,j,k-3) -   45*Ay(i+1,j,k-3) +   9*Ay(i+2,j,k-3) -    Ay(i+3,j,k-3) ) * odxdz3600
      d2_lA(3,1,3) = (    -Az(i-3,j,k+3) +   9*Az(i-2,j,k+3) -   45*Az(i-1,j,k+3) +   45*Az(i+1,j,k+3) -   9*Az(i+2,j,k+3) +    Az(i+3,j,k+3) &
                      +  9*Az(i-3,j,k+2) -  81*Az(i-2,j,k+2) +  405*Az(i-1,j,k+2) -  405*Az(i+1,j,k+2) +  81*Az(i+2,j,k+2) -  9*Az(i+3,j,k+2) &
                      - 45*Az(i-3,j,k+1) + 405*Az(i-2,j,k+1) - 2025*Az(i-1,j,k+1) + 2025*Az(i+1,j,k+1) - 405*Az(i+2,j,k+1) + 45*Az(i+3,j,k+1) &
                      + 45*Az(i-3,j,k-1) - 405*Az(i-2,j,k-1) + 2025*Az(i-1,j,k-1) - 2025*Az(i+1,j,k-1) + 405*Az(i+2,j,k-1) - 45*Az(i+3,j,k-1) &
                      -  9*Az(i-3,j,k-2) +  81*Az(i-2,j,k-2) -  405*Az(i-1,j,k-2) +  405*Az(i+1,j,k-2) -  81*Az(i+2,j,k-2) +  9*Az(i+3,j,k-2) &
                      +    Az(i-3,j,k-3) -   9*Az(i-2,j,k-3) +   45*Az(i-1,j,k-3) -   45*Az(i+1,j,k-3) +   9*Az(i+2,j,k-3) -    Az(i+3,j,k-3) ) * odxdz3600

      d2_lA(1,2,3) = (    -Ax(i,j-3,k+3) +   9*Ax(i,j-2,k+3) -   45*Ax(i,j-1,k+3) +   45*Ax(i,j+1,k+3) -   9*Ax(i,j+2,k+3) +    Ax(i,j+3,k+3) &
                      +  9*Ax(i,j-3,k+2) -  81*Ax(i,j-2,k+2) +  405*Ax(i,j-1,k+2) -  405*Ax(i,j+1,k+2) +  81*Ax(i,j+2,k+2) -  9*Ax(i,j+3,k+2) &
                      - 45*Ax(i,j-3,k+1) + 405*Ax(i,j-2,k+1) - 2025*Ax(i,j-1,k+1) + 2025*Ax(i,j+1,k+1) - 405*Ax(i,j+2,k+1) + 45*Ax(i,j+3,k+1) &
                      + 45*Ax(i,j-3,k-1) - 405*Ax(i,j-2,k-1) + 2025*Ax(i,j-1,k-1) - 2025*Ax(i,j+1,k-1) + 405*Ax(i,j+2,k-1) - 45*Ax(i,j+3,k-1) &
                      -  9*Ax(i,j-3,k-2) +  81*Ax(i,j-2,k-2) -  405*Ax(i,j-1,k-2) +  405*Ax(i,j+1,k-2) -  81*Ax(i,j+2,k-2) +  9*Ax(i,j+3,k-2) &
                      +    Ax(i,j-3,k-3) -   9*Ax(i,j-2,k-3) +   45*Ax(i,j-1,k-3) -   45*Ax(i,j+1,k-3) +   9*Ax(i,j+2,k-3) -    Ax(i,j+3,k-3) ) * odydz3600
      d2_lA(2,2,3) = (    -Ay(i,j-3,k+3) +   9*Ay(i,j-2,k+3) -   45*Ay(i,j-1,k+3) +   45*Ay(i,j+1,k+3) -   9*Ay(i,j+2,k+3) +    Ay(i,j+3,k+3) &
                      +  9*Ay(i,j-3,k+2) -  81*Ay(i,j-2,k+2) +  405*Ay(i,j-1,k+2) -  405*Ay(i,j+1,k+2) +  81*Ay(i,j+2,k+2) -  9*Ay(i,j+3,k+2) &
                      - 45*Ay(i,j-3,k+1) + 405*Ay(i,j-2,k+1) - 2025*Ay(i,j-1,k+1) + 2025*Ay(i,j+1,k+1) - 405*Ay(i,j+2,k+1) + 45*Ay(i,j+3,k+1) &
                      + 45*Ay(i,j-3,k-1) - 405*Ay(i,j-2,k-1) + 2025*Ay(i,j-1,k-1) - 2025*Ay(i,j+1,k-1) + 405*Ay(i,j+2,k-1) - 45*Ay(i,j+3,k-1) &
                      -  9*Ay(i,j-3,k-2) +  81*Ay(i,j-2,k-2) -  405*Ay(i,j-1,k-2) +  405*Ay(i,j+1,k-2) -  81*Ay(i,j+2,k-2) +  9*Ay(i,j+3,k-2) &
                     +    Ay(i,j-3,k-3) -   9*Ay(i,j-2,k-3) +   45*Ay(i,j-1,k-3) -   45*Ay(i,j+1,k-3) +   9*Ay(i,j+2,k-3) -    Ay(i,j+3,k-3) ) * odydz3600
      d2_lA(3,2,3) = (    -Az(i,j-3,k+3) +   9*Az(i,j-2,k+3) -   45*Az(i,j-1,k+3) +   45*Az(i,j+1,k+3) -   9*Az(i,j+2,k+3) +    Az(i,j+3,k+3) &
                      +  9*Az(i,j-3,k+2) -  81*Az(i,j-2,k+2) +  405*Az(i,j-1,k+2) -  405*Az(i,j+1,k+2) +  81*Az(i,j+2,k+2) -  9*Az(i,j+3,k+2) &
                      - 45*Az(i,j-3,k+1) + 405*Az(i,j-2,k+1) - 2025*Az(i,j-1,k+1) + 2025*Az(i,j+1,k+1) - 405*Az(i,j+2,k+1) + 45*Az(i,j+3,k+1) &
                      + 45*Az(i,j-3,k-1) - 405*Az(i,j-2,k-1) + 2025*Az(i,j-1,k-1) - 2025*Az(i,j+1,k-1) + 405*Az(i,j+2,k-1) - 45*Az(i,j+3,k-1) &
                      -  9*Az(i,j-3,k-2) +  81*Az(i,j-2,k-2) -  405*Az(i,j-1,k-2) +  405*Az(i,j+1,k-2) -  81*Az(i,j+2,k-2) +  9*Az(i,j+3,k-2) &
                      +    Az(i,j-3,k-3) -   9*Az(i,j-2,k-3) +   45*Az(i,j-1,k-3) -   45*Az(i,j+1,k-3) +   9*Az(i,j+2,k-3) -    Az(i,j+3,k-3) ) * odydz3600

      d2_lA(:,2,1) = d2_lA(:,1,2)
      d2_lA(:,3,1) = d2_lA(:,1,3)
      d2_lA(:,3,2) = d2_lA(:,2,3)


      !------------ Advection derivatives --------
      if( use_advection_stencils /= 0 ) then

        di = int( sign( one, beta(1) ) )
        dj = int( sign( one, beta(2) ) )
        dk = int( sign( one, beta(3) ) )

        ! ad1_lE(3)
        d1_f(1) = di * (   2*Ex(i-2*di,j,k) - 24*Ex(i-di,j,k) - 35*Ex(i,j,k) + 80*Ex(i+di,j,k) &
                        - 30*Ex(i+2*di,j,k) + 8*Ex(i+3*di,j,k) - Ex(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * (   2*Ex(i,j-2*dj,k) - 24*Ex(i,j-dj,k) - 35*Ex(i,j,k) + 80*Ex(i,j+dj,k) &
                        - 30*Ex(i,j+2*dj,k) + 8*Ex(i,j+3*dj,k) - Ex(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * (   2*Ex(i,j,k-2*dk) - 24*Ex(i,j,k-dk) - 35*Ex(i,j,k) + 80*Ex(i,j,k+dk) &
                        - 30*Ex(i,j,k+2*dk) + 8*Ex(i,j,k+3*dk) - Ex(i,j,k+4*dk) ) * odz60
        ad1_lE(1) = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        d1_f(1) = di * (   2*Ey(i-2*di,j,k) - 24*Ey(i-di,j,k) - 35*Ey(i,j,k) + 80*Ey(i+di,j,k) &
                        - 30*Ey(i+2*di,j,k) +  8*Ey(i+3*di,j,k) -  Ey(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * (   2*Ey(i,j-2*dj,k) - 24*Ey(i,j-dj,k) - 35*Ey(i,j,k) + 80*Ey(i,j+dj,k) &
                        - 30*Ey(i,j+2*dj,k) +  8*Ey(i,j+3*dj,k) -  Ey(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * (   2*Ey(i,j,k-2*dk) - 24*Ey(i,j,k-dk) - 35*Ey(i,j,k) + 80*Ey(i,j,k+dk) &
                        - 30*Ey(i,j,k+2*dk) +  8*Ey(i,j,k+3*dk) -  Ey(i,j,k+4*dk) ) * odz60
        ad1_lE(2) = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        d1_f(1) = di * (   2*Ez(i-2*di,j,k) - 24*Ez(i-di,j,k) - 35*Ez(i,j,k) + 80*Ez(i+di,j,k) &
                        - 30*Ez(i+2*di,j,k) +  8*Ez(i+3*di,j,k) -  Ez(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * (   2*Ez(i,j-2*dj,k) - 24*Ez(i,j-dj,k) - 35*Ez(i,j,k) + 80*Ez(i,j+dj,k) &
                        - 30*Ez(i,j+2*dj,k) +  8*Ez(i,j+3*dj,k) -  Ez(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * (   2*Ez(i,j,k-2*dk) - 24*Ez(i,j,k-dk) - 35*Ez(i,j,k) + 80*Ez(i,j,k+dk) &
                        - 30*Ez(i,j,k+2*dk) +  8*Ez(i,j,k+3*dk) -  Ez(i,j,k+4*dk) ) * odz60
        ad1_lE(3) = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        ! ad1_lA(3)
        d1_f(1) = di * (   2*Ax(i-2*di,j,k) - 24*Ax(i-di,j,k) - 35*Ax(i,j,k) + 80*Ax(i+di,j,k) &
                        - 30*Ax(i+2*di,j,k) + 8*Ax(i+3*di,j,k) - Ax(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * (   2*Ax(i,j-2*dj,k) - 24*Ax(i,j-dj,k) - 35*Ax(i,j,k) + 80*Ax(i,j+dj,k) &
                        - 30*Ax(i,j+2*dj,k) + 8*Ax(i,j+3*dj,k) - Ax(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * (   2*Ax(i,j,k-2*dk) - 24*Ax(i,j,k-dk) - 35*Ax(i,j,k) + 80*Ax(i,j,k+dk) &
                        - 30*Ax(i,j,k+2*dk) + 8*Ax(i,j,k+3*dk) - Ax(i,j,k+4*dk) ) * odz60
        ad1_lA(1) = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        d1_f(1) = di * (   2*Ay(i-2*di,j,k) - 24*Ay(i-di,j,k) - 35*Ay(i,j,k) + 80*Ay(i+di,j,k) &
                        - 30*Ay(i+2*di,j,k) +  8*Ay(i+3*di,j,k) -  Ay(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * (   2*Ay(i,j-2*dj,k) - 24*Ay(i,j-dj,k) - 35*Ay(i,j,k) + 80*Ay(i,j+dj,k) &
                        - 30*Ay(i,j+2*dj,k) +  8*Ay(i,j+3*dj,k) -  Ay(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * (   2*Ay(i,j,k-2*dk) - 24*Ay(i,j,k-dk) - 35*Ay(i,j,k) + 80*Ay(i,j,k+dk) &
                        - 30*Ay(i,j,k+2*dk) +  8*Ay(i,j,k+3*dk) -  Ay(i,j,k+4*dk) ) * odz60
        ad1_lA(2) = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        d1_f(1) = di * (   2*Az(i-2*di,j,k) - 24*Az(i-di,j,k) - 35*Az(i,j,k) + 80*Az(i+di,j,k) &
                        - 30*Az(i+2*di,j,k) +  8*Az(i+3*di,j,k) -  Az(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * (   2*Az(i,j-2*dj,k) - 24*Az(i,j-dj,k) - 35*Az(i,j,k) + 80*Az(i,j+dj,k) &
                        - 30*Az(i,j+2*dj,k) +  8*Az(i,j+3*dj,k) -  Az(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * (   2*Az(i,j,k-2*dk) - 24*Az(i,j,k-dk) - 35*Az(i,j,k) + 80*Az(i,j,k+dk) &
                        - 30*Az(i,j,k+2*dk) +  8*Az(i,j,k+3*dk) -  Az(i,j,k+4*dk) ) * odz60
        ad1_lA(3) = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        ! ad1_lZeta
        d1_f(1) = di * (   2*Zeta(i-2*di,j,k) - 24*Zeta(i-di,j,k) - 35*Zeta(i,j,k) + 80*Zeta(i+di,j,k) &
                        - 30*Zeta(i+2*di,j,k) +  8*Zeta(i+3*di,j,k) -  Zeta(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * (   2*Zeta(i,j-2*dj,k) - 24*Zeta(i,j-dj,k) - 35*Zeta(i,j,k) + 80*Zeta(i,j+dj,k) &
                        - 30*Zeta(i,j+2*dj,k) +  8*Zeta(i,j+3*dj,k) -  Zeta(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * (   2*Zeta(i,j,k-2*dk) - 24*Zeta(i,j,k-dk) - 35*Zeta(i,j,k) + 80*Zeta(i,j,k+dk) &
                        - 30*Zeta(i,j,k+2*dk) +  8*Zeta(i,j,k+3*dk) -  Zeta(i,j,k+4*dk) ) * odz60
        ad1_lZeta = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        ! ad1_lAphi
        d1_f(1) = di * (   2*Aphi(i-2*di,j,k) - 24*Aphi(i-di,j,k) - 35*Aphi(i,j,k) + 80*Aphi(i+di,j,k) &
                        - 30*Aphi(i+2*di,j,k) +  8*Aphi(i+3*di,j,k) -  Aphi(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * (   2*Aphi(i,j-2*dj,k) - 24*Aphi(i,j-dj,k) - 35*Aphi(i,j,k) + 80*Aphi(i,j+dj,k) &
                        - 30*Aphi(i,j+2*dj,k) +  8*Aphi(i,j+3*dj,k) -  Aphi(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * (   2*Aphi(i,j,k-2*dk) - 24*Aphi(i,j,k-dk) - 35*Aphi(i,j,k) + 80*Aphi(i,j,k+dk) &
                        - 30*Aphi(i,j,k+2*dk) +  8*Aphi(i,j,k+3*dk) -  Aphi(i,j,k+4*dk) ) * odz60
        ad1_lAphi = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

      else

        ! ad1_lE(3)
        ad1_lE(1) = beta(1)*d1_lE(1,1) + beta(2)*d1_lE(1,2) + beta(3)*d1_lE(1,3)
        ad1_lE(2) = beta(1)*d1_lE(2,1) + beta(2)*d1_lE(2,2) + beta(3)*d1_lE(2,3)
        ad1_lE(3) = beta(1)*d1_lE(3,1) + beta(2)*d1_lE(3,2) + beta(3)*d1_lE(3,3)

        ! ad1_lA(3)
        ad1_lA(1) = beta(1)*d1_lA(1,1) + beta(2)*d1_lA(1,2) + beta(3)*d1_lA(1,3)
        ad1_lA(2) = beta(1)*d1_lA(2,1) + beta(2)*d1_lA(2,2) + beta(3)*d1_lA(2,3)
        ad1_lA(3) = beta(1)*d1_lA(3,1) + beta(2)*d1_lA(3,2) + beta(3)*d1_lA(3,3)

        ! ad1_lZeta
        ad1_lZeta = beta(1)*d1_lZeta(1) + beta(2)*d1_lZeta(2) + beta(3)*d1_lZeta(3)

        ! ad1_lAphi
        ad1_lAphi = beta(1)*d1_lAphi(1) + beta(2)*d1_lAphi(2) + beta(3)*d1_lAphi(3)

      end if
      !-------------------------------------------

    else
      call CCTK_WARN(0, "derivs_order not yet implemented.")
    end if

    !------------ Christoffel symbols ----------
    cf1 = 0
    do a = 1, 3
      do b = 1, 3
        do c = b, 3
          cf1(a,b,c) = 0.5d0 * (d1_hh(a,b,c) + d1_hh(a,c,b) - d1_hh(b,c,a))
        end do
      end do
    end do
    cf1(:,2,1) = cf1(:,1,2)
    cf1(:,3,1) = cf1(:,1,3)
    cf1(:,3,2) = cf1(:,2,3)

    cf2 = 0
    do a = 1, 3
      do b = 1, 3
        do c = b, 3
          do m = 1, 3
            cf2(a,b,c) = cf2(a,b,c) + hu(a,m) * cf1(m,b,c)
          end do
        end do
      end do
    end do
    cf2(:,2,1) = cf2(:,1,2)
    cf2(:,3,1) = cf2(:,1,3)
    cf2(:,3,2) = cf2(:,2,3)
    !-------------------------------------------


    !------------ Covariant derivatives --------
    cd_lA = d1_lA
    cd2_lphi1  = d2_lphi1
    cd2_lphi2  = d2_lphi2
    do a = 1, 3
      do b = 1, 3
        do m = 1, 3
          cd_lA(a,b) = cd_lA(a,b) - cf2(m,a,b) * lA(m)
          cd2_lphi1(a,b)  = cd2_lphi1(a,b) - cf2(m,a,b) * d1_lphi1(m)
          cd2_lphi2(a,b)  = cd2_lphi2(a,b) - cf2(m,a,b) * d1_lphi2(m)
       end do
      end do
   end do

   cd_dA = d2_lA
   do a = 1, 3
      do b = 1, 3
        do c = 1, 3
          do m = 1, 3
            cd_dA(a,b,c) = cd_dA(a,b,c) - cf2(m,a,c) * d1_lA(m,b)         &
                                        - cf2(m,b,c) * d1_lA(a,m)
           end do
         end do
       end do
    end do
    !-------------------------------------------


    !------------ Definition of Maxwell souces  --------
    rho_e = -2 * q * (2 * lphi1 * lKphi2 - 2 * lphi2 * lKphi1 + q * lAphi * (lphi1*lphi1 + lphi2*lphi2))
    do a = 1, 3
       J_i(a)=-2*q*(lphi1 * d1_lphi2(a) - lphi2 * d1_lphi1(a) + q * lA(a) * (lphi1*lphi1 + lphi2*lphi2))
    end do
    !-------------------------------------------


    !--------- Evolution of E, A, Aphi         ----------

    ! rhs_lE
    rhs_lE  = ad1_lE

    do a = 1, 3
      do m = 1, 3
        rhs_lE(a) = rhs_lE(a) - lE(m) * d1_beta(a,m)
      end do
    end do

    rhs_lE  = rhs_lE + alph * trk * lE

    do a = 1, 3
       do b = 1, 3
          rhs_lE(a) = rhs_lE(a) + alph * ch**conf_fac_exponent * hu(a,b) * d1_lZeta(b)      &
                                + alph * mu_V*mu_V * ch**conf_fac_exponent * hu(a,b) * lA(b)&
                                - alph * ch**conf_fac_exponent * hu(a,b) * J_i(b)  
          do c = 1, 3
             do m = 1, 3
                rhs_lE(a) = rhs_lE(a)                                                       &
                     + ch**(2*conf_fac_exponent) * hu(a,b) * hu(c,m) * d1_alph(m)           &
                                                 * ( d1_lA(c,b) - d1_lA(b,c) )              &
                     + alph * ch**(2*conf_fac_exponent) * hu(a,b) * hu(c,m)                 &
                                                 * ( cd_dA(m,b,c) - cd_dA(b,m,c) )          &
                     + 0.50d0*conf_fac_exponent * alph * ch**(2*conf_fac_exponent - 1)      &
                        * hu(a,b) * hu(c,m) * ( d1_lA(m,b)*d1_ch(c) - d1_lA(b,c)*d1_ch(m) )
             end do
          end do
       end do
    end do


    ! rhs_lA
    rhs_lA  = ad1_lA

    do a = 1, 3
      do m = 1, 3
        rhs_lA(a) = rhs_lA(a) + lA(m) * d1_beta(m,a)
      end do
    end do

    rhs_lA  = rhs_lA - alph * d1_lAphi - lAphi * d1_alph

    do a = 1,3
       do b = 1,3
          rhs_lA(a) = rhs_lA(a) - alph * hh(a,b) * lE(b) * ch**(-conf_fac_exponent)
       end do
    end do


    ! rhs_lAphi
    rhs_lAphi  = ad1_lAphi  + alph * trk * lAphi - alph * lZeta

    do a = 1, 3
        do b = 1, 3
           rhs_lAphi = rhs_lAphi                                                     &
                     + 0.5d0 * alph * conf_fac_exponent * ch**(conf_fac_exponent-1)  &
                             * hu(a,b) * lA(a) * d1_ch(b)                            &
                     - alph * ch**conf_fac_exponent * hu(a,b) * cd_lA(a,b)           &
                     - ch**conf_fac_exponent * hu(a,b) * lA(a) * d1_alph(b)
        end do
    end do

    
    !--------- Evolution of phi, Kphi ----------

    ! rhs_lphi1, rhs_lphi2
    rhs_lphi1  = ad1_lphi1
    rhs_lphi2  = ad1_lphi2

    rhs_lphi1 = rhs_lphi1 - 2.0d0 * alph * lKphi1

    rhs_lphi2 = rhs_lphi2 - 2.0d0 * alph * lKphi2


    ! rhs_lKphi1, rhs_lKphi2
    rhs_lKphi1 = ad1_lKphi1
    rhs_lKphi2 = ad1_lKphi2

    tr_dalp_dphi1 = 0
    tr_dalp_dphi2 = 0
    tr_cd2_phi1   = 0
    tr_cd2_phi2   = 0
    tr_dch_dphi1  = 0
    tr_dch_dphi2  = 0
    tr_A_A        = 0
    tr_A_dphi1    = 0
    tr_A_dphi2    = 0
    do a = 1, 3
      do b = 1, 3
        tr_dalp_dphi1 = tr_dalp_dphi1 + hu(a,b) * d1_alph(a) * d1_lphi1(b)
        tr_dalp_dphi2 = tr_dalp_dphi2 + hu(a,b) * d1_alph(a) * d1_lphi2(b)
        tr_cd2_phi1   = tr_cd2_phi1   + hu(a,b) * cd2_lphi1(a,b)
        tr_cd2_phi2   = tr_cd2_phi2   + hu(a,b) * cd2_lphi2(a,b)
        tr_dch_dphi1  = tr_dch_dphi1  + hu(a,b) * d1_ch(a) * d1_lphi1(b)
        tr_dch_dphi2  = tr_dch_dphi2  + hu(a,b) * d1_ch(a) * d1_lphi2(b)
        tr_A_A        = tr_A_A        + hu(a,b) * lA(a) * lA(b)
        tr_A_dphi1    = tr_A_dphi1    + hu(a,b) * lA(a) * d1_lphi1(b)
        tr_A_dphi2    = tr_A_dphi2    + hu(a,b) * lA(a) * d1_lphi2(b)
      end do
    end do

    rhs_lKphi1 = rhs_lKphi1 - 0.5d0 * ch * tr_dalp_dphi1                                                      &
                 + 0.5d0 * alph * ( - ch * tr_cd2_phi1 + 0.5d0 * tr_dch_dphi1                                 &
                                    + mu*mu * lphi1 * ( 1 - 8*V_lambda*( lphi1*lphi1 + lphi2*lphi2 )          &
                                    + 12*V_lambda*V_lambda                                                    &
                                    * ( lphi1*lphi1 + lphi2*lphi2 ) * ( lphi1*lphi1 + lphi2*lphi2 ) )         &
                                    + 2 * trk * lKphi1 )                                                      &
                 + 0.5d0 * q*q * alph * ( ch * tr_A_A - lAphi*lAphi) * lphi1                                  &
                 - q * alph * ( - ch * tr_A_dphi2 + 2 * lAphi * lKphi2 )                   

    rhs_lKphi2 = rhs_lKphi2 - 0.5d0 * ch * tr_dalp_dphi2                                                      &
                 + 0.5d0 * alph * ( - ch * tr_cd2_phi2 + 0.5d0 * tr_dch_dphi2                                 &
                                    + mu*mu * lphi2 * ( 1 - 8*V_lambda*( lphi1*lphi1 + lphi2*lphi2 )          &
                                    + 12*V_lambda*V_lambda                                                    &
                                    * ( lphi1*lphi1 + lphi2*lphi2 ) * ( lphi1*lphi1 + lphi2*lphi2 ) )         &
                                    + 2 * trk * lKphi2 )                                                      &
                 + 0.5d0 * q*q * alph * ( ch * tr_A_A - lAphi*lAphi) * lphi2                                  &
                 - q * alph * ( ch * tr_A_dphi1 - 2 * lAphi * lKphi1 )                  


    !--------- Evolution of Zeta ----------
    
    ! rhs_lZeta
    rhs_lZeta = ad1_lZeta - alph * kappa * lZeta + alph * mu_V*mu_V * lAphi - alph * rho_e

    do a = 1, 3
       rhs_lZeta = rhs_lZeta + alph * d1_lE(a,a)                         &
                 - 1.5d0 * conf_fac_exponent * alph * lE(a) * d1_ch(a) / ch
    end do

    rhs_lZeta = Zeta_Omega_fac * rhs_lZeta

    !-------------------------------------------

    ! if( abs(y(i,j,k)) < 1.0d-05 .and. abs(z(i,j,k)) < 1.0d-05 ) then
    !    write(*,*) 'i, j, k    = ', i, j, k
    !    write(*,*) 'x          = ', x(i,j,k)
    !    write(*,*) 'Ex         = ', lE(1)
    !    write(*,*) 'rhs_lE     = ', rhs_lE
    !    write(*,*) 'rhs_lA     = ', rhs_lA
    !    write(*,*) 'rhs_lAphi  = ', rhs_lAphi
    !    write(*,*) 'rhs_lZeta  = ', rhs_lZeta
    !    call flush(6)
    ! end if


    !-------- Write to the gridfunctions -------

    rhs_Ex(i,j,k) = rhs_lE(1)
    rhs_Ey(i,j,k) = rhs_lE(2)
    rhs_Ez(i,j,k) = rhs_lE(3)

    rhs_Ax(i,j,k) = rhs_lA(1)
    rhs_Ay(i,j,k) = rhs_lA(2)
    rhs_Az(i,j,k) = rhs_lA(3)

    rhs_Zeta(i,j,k) = rhs_lZeta

    rhs_Aphi(i,j,k) = rhs_lAphi

    rhs_phi1(i,j,k)  = rhs_lphi1
    rhs_phi2(i,j,k)  = rhs_lphi2
    rhs_Kphi1(i,j,k) = rhs_lKphi1
    rhs_Kphi2(i,j,k) = rhs_lKphi2

  end do
  end do
  end do
  !$OMP END PARALLEL DO

end subroutine MagScalar_calc_rhs
!
!=============================================================================
!

subroutine MagScalar_calc_rhs_bdry( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL, parameter :: one  = 1.0
  CCTK_REAL, parameter :: zero = 0.0
  CCTK_INT ierr

  ierr = NewRad_Apply(cctkGH, Ex, rhs_Ex, E0(1), one, n_E(1))
  ierr = NewRad_Apply(cctkGH, Ey, rhs_Ey, E0(2), one, n_E(2))
  ierr = NewRad_Apply(cctkGH, Ez, rhs_Ez, E0(3), one, n_E(3))

  ierr = NewRad_Apply(cctkGH, Ax, rhs_Ax, A0(1), one, n_A(1))
  ierr = NewRad_Apply(cctkGH, Ay, rhs_Ay, A0(2), one, n_A(2))
  ierr = NewRad_Apply(cctkGH, Az, rhs_Az, A0(3), one, n_A(3))

  ierr = NewRad_Apply(cctkGH, Aphi, rhs_Aphi, Aphi0, one, n_Aphi)
  ierr = NewRad_Apply(cctkGH, Zeta, rhs_Zeta, zero, one, n_Zeta)

  ierr = NewRad_Apply(cctkGH, phi1, rhs_phi1, phi1_0, one, n_phi1)
  ierr = NewRad_Apply(cctkGH, phi2, rhs_phi2, phi2_0, one, n_phi2)

  ierr = NewRad_Apply(cctkGH, Kphi1, rhs_Kphi1, Kphi1_0, one, n_Kphi1)
  ierr = NewRad_Apply(cctkGH, Kphi2, rhs_Kphi2, Kphi2_0, one, n_Kphi2)

end subroutine MagScalar_calc_rhs_bdry
