#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

subroutine MagScalar_calc_Tmunu( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  ! Fundamental variables
  CCTK_REAL                alph, beta(3)
  CCTK_REAL                gg(3,3), gu(3,3), detgg
  CCTK_REAL                lE(3), lB(3), lA(3), lAphi
  CCTK_REAL                Ed(3), Bd(3)
  CCTK_REAL                Tab(4,4)

  ! First derivatives
  CCTK_REAL                d1_lA(3,3)

  ! Auxiliary variables
  CCTK_REAL                eps_lc_d(3,3,3), eps_lc_u(3,3,3)

  ! Matter variables
  CCTK_REAL                srcE, srcjdi(3), srcSij(3,3)

  ! Misc variables
  CCTK_REAL                dx12, dy12, dz12
  CCTK_REAL                odx60, ody60, odz60
  CCTK_REAL                aux

  CCTK_REAL, parameter ::  one  = 1
  CCTK_REAL, parameter ::  pi   = acos(-one)
  CCTK_REAL, parameter ::  pi4  = 4*pi
  CCTK_REAL, parameter ::  pi8  = 8*pi
  CCTK_INT                 i, j, k
  CCTK_INT                 a, b, c, m

  dx12 = 12*CCTK_DELTA_SPACE(1)
  dy12 = 12*CCTK_DELTA_SPACE(2)
  dz12 = 12*CCTK_DELTA_SPACE(3)

  odx60 = 1 / (60 * CCTK_DELTA_SPACE(1))
  ody60 = 1 / (60 * CCTK_DELTA_SPACE(2))
  odz60 = 1 / (60 * CCTK_DELTA_SPACE(3))

  !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(k,j,i,a,b,c,m,aux,&
  !$OMP                                 alph,beta,&
  !$OMP                                 gg,gu,detgg,&
  !$OMP                                 lE,lB,lA,lAphi,&
  !$OMP                                 Ed,Bd,Tab,&
  !$OMP                                 eps_lc_d,eps_lc_u,&
  !$OMP                                 d1_lA, &
  !$OMP                                 srcE, srcjdi, srcSij)
  do k = 1+cctk_nghostzones(3), cctk_lsh(3)-cctk_nghostzones(3)
     do j = 1+cctk_nghostzones(2), cctk_lsh(2)-cctk_nghostzones(2)
        do i = 1+cctk_nghostzones(1), cctk_lsh(1)-cctk_nghostzones(1)

           !------------ Get local variables ----------

           alph      = alp(i,j,k)

           beta(1)   = betax(i,j,k)
           beta(2)   = betay(i,j,k)
           beta(3)   = betaz(i,j,k)

           gg(1,1)   = gxx(i,j,k)
           gg(1,2)   = gxy(i,j,k)
           gg(1,3)   = gxz(i,j,k)
           gg(2,2)   = gyy(i,j,k)
           gg(2,3)   = gyz(i,j,k)
           gg(3,3)   = gzz(i,j,k)
           gg(2,1)   = gg(1,2)
           gg(3,1)   = gg(1,3)
           gg(3,2)   = gg(2,3)

           lE(1)     = Ex(i,j,k)
           lE(2)     = Ey(i,j,k)
           lE(3)     = Ez(i,j,k)

           lA(1)     = Ax(i,j,k)
           lA(2)     = Ay(i,j,k)
           lA(3)     = Az(i,j,k)

           lAphi     = Aphi(i,j,k)

           Ed(1)     = gg(1,1) * lE(1) + gg(1,2) * lE(2) + gg(1,3) * lE(3)
           Ed(2)     = gg(2,1) * lE(1) + gg(2,2) * lE(2) + gg(2,3) * lE(3)
           Ed(3)     = gg(3,1) * lE(1) + gg(3,2) * lE(2) + gg(3,3) * lE(3)


           !------------ Invert 3-metric ----------------
           detgg   =     gg(1,1) * gg(2,2) * gg(3,3)                              &
                   + 2 * gg(1,2) * gg(1,3) * gg(2,3)                              &
                   -     gg(1,1) * gg(2,3) ** 2                                   &
                   -     gg(2,2) * gg(1,3) ** 2                                   &
                   -     gg(3,3) * gg(1,2) ** 2

           gu(1,1) = (gg(2,2) * gg(3,3) - gg(2,3) ** 2     ) / detgg
           gu(2,2) = (gg(1,1) * gg(3,3) - gg(1,3) ** 2     ) / detgg
           gu(3,3) = (gg(1,1) * gg(2,2) - gg(1,2) ** 2     ) / detgg
           gu(1,2) = (gg(1,3) * gg(2,3) - gg(1,2) * gg(3,3)) / detgg
           gu(1,3) = (gg(1,2) * gg(2,3) - gg(1,3) * gg(2,2)) / detgg
           gu(2,3) = (gg(1,3) * gg(1,2) - gg(2,3) * gg(1,1)) / detgg
           gu(2,1) = gu(1,2)
           gu(3,1) = gu(1,3)
           gu(3,2) = gu(2,3)
           !-------------------------------------------------


           !------------- Centered 1st derivatives ----------

           if (derivs_order == 4) then

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

           else if (derivs_order == 6) then

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

           else
             call CCTK_WARN(0, "derivs_order not yet implemented.")
           end if


           !------------ Levi-Civita tensor ----------
           eps_lc_u        =  0
           eps_lc_u(1,2,3) =  1
           eps_lc_u(2,3,1) =  1
           eps_lc_u(3,1,2) =  1
           eps_lc_u(3,2,1) = -1
           eps_lc_u(2,1,3) = -1
           eps_lc_u(1,3,2) = -1
           eps_lc_u = eps_lc_u / sqrt(detgg)
           eps_lc_d = eps_lc_u * detgg
           !------------------------------------------


           ! magnetic field B (here used as an auxiliary variable)
           lB = 0
           do a = 1, 3
             do b = 1, 3
               do m = 1, 3
                 lB(a) = lB(a) + eps_lc_u(a,b,m) * d1_lA(m,b)
               end do
             end do
           end do

           Bd(1) = gg(1,1) * lB(1) + gg(1,2) * lB(2) + gg(1,3) * lB(3)
           Bd(2) = gg(2,1) * lB(1) + gg(2,2) * lB(2) + gg(2,3) * lB(3)
           Bd(3) = gg(3,1) * lB(1) + gg(3,2) * lB(2) + gg(3,3) * lB(3)
           !-------------------------------------------

           !------------ Matter terms -----------------
           !
           ! mu = 0, 1, 2, 3; i,a = 1,2,3
           !
           ! n_mu = (-alph, 0, 0, 0)
           ! n^mu = (1, -betax, -betay, -betaz)/alph
           !
           ! rho = n^mu n^nu T_{mu nu}
           !     = (T_{00} - 2 beta^i T_{i0} + beta^i beta^j T_{ij})/alph^2
           !
           ! j_a = -h_a^mu n^nu T_{mu nu}
           !     = -(T_{a 0} - beta^j T_{a j})/alph
           !
           ! S_{a b} = h_{a mu} h_{b nu} T^{mu nu} = T_{a b}


           ! srcE = rho = n^mu n^nu T_{mu nu}
           srcE = mu*mu * lAphi*lAphi
           do a = 1, 3
              do b = 1, 3
                 srcE = srcE + ( lE(a) * lE(b) + lB(a) * lB(b) ) * gg(a,b)    &
                             + mu*mu * lA(a) * lA(b) * gu(a,b)
              end do
           end do
           srcE = srcE / pi8

           !srcjdi = j_a
           srcjdi = mu*mu * lAphi * lA
           do a = 1, 3
              do b = 1, 3
                 do m = 1, 3
                    srcjdi(a) = srcjdi(a) + eps_lc_d(a,b,m) * lE(b) * lB(m)
                 end do
              end do
           end do
           srcjdi = srcjdi / pi4


           ! srcSij = S_{a b}

           aux = 0
           do a = 1, 3
              do b = 1, 3
                 aux = aux + ( lE(a) * lE(b) + lB(a) * lB(b) ) * gg(a,b)    &
                           - mu*mu * lA(a) * lA(b) * gu(a,b)
              end do
           end do

           srcSij = 0.5 * (aux + mu*mu * lAphi*lAphi) * gg
           do a = 1, 3
              do b = 1, 3
                 srcSij(a,b) = srcSij(a,b) - Ed(a) * Ed(b) - Bd(a) * Bd(b)  &
                             + mu*mu * lA(a) * lA(b)
              end do
           end do
           srcSij = srcSij / pi4

           !------------------------------------------


           ! now to fill in the stress-energy tensor. note that we use Tab(4,4)
           ! for T_{0 0}
           !
           ! T_{a b} = S_{a b}
           !
           ! T_{0 0} = alph^2 rho - 2 alph beta^a j_a + beta^a beta^b S_{a b}
           !
           ! T_{0 a} = -alph j_a + beta^b S_{a b}

           Tab(1:3,1:3) = srcSij(1:3,1:3)

           Tab(1:3,4) = -alph * srcjdi(1:3)
           do b = 1, 3
              Tab(1:3,4) = Tab(1:3,4) + beta(b) * srcSij(1:3,b)
           end do
           Tab(4,1:3) = Tab(1:3,4)

           Tab(4,4) = alph**2 * srcE
           do a = 1, 3
              Tab(4,4) = Tab(4,4) - 2 * alph * beta(a) * srcjdi(a)
              do b = 1, 3
                 Tab(4,4) = Tab(4,4) + beta(a) * beta(b) * srcSij(a,b)
              end do
           end do

           ! and finally store it in the Tmunu variables
           eTtt(i,j,k) = eTtt(i,j,k) + Tab(4,4)
           eTtx(i,j,k) = eTtx(i,j,k) + Tab(4,1)
           eTty(i,j,k) = eTty(i,j,k) + Tab(4,2)
           eTtz(i,j,k) = eTtz(i,j,k) + Tab(4,3)
           eTxx(i,j,k) = eTxx(i,j,k) + Tab(1,1)
           eTxy(i,j,k) = eTxy(i,j,k) + Tab(1,2)
           eTxz(i,j,k) = eTxz(i,j,k) + Tab(1,3)
           eTyy(i,j,k) = eTyy(i,j,k) + Tab(2,2)
           eTyz(i,j,k) = eTyz(i,j,k) + Tab(2,3)
           eTzz(i,j,k) = eTzz(i,j,k) + Tab(3,3)

        end do
     end do
  end do
  !$OMP END PARALLEL DO

end subroutine MagScalar_calc_Tmunu
