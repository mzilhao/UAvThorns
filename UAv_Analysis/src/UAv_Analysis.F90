
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#include "SpaceMask.h"

subroutine UAv_Analysis_gfs( CCTK_ARGUMENTS )
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL alph, beta(3), Tab(4,4)
  CCTK_REAL gd(3,3), gu(3,3), detgd
  
  CCTK_REAL aux, S, rho
  CCTK_REAL mom(3)
  
  CCTK_REAL x0, y0, z0, x1, y1, z1, xtmp, ytmp, ztmp
  CCTK_INT gs_index_tmp_x, gs_index_tmp_y, gs_index_tmp_z
  
  CCTK_POINTER ori_ptr_x, ori_ptr_y, ori_ptr_z
  pointer (ori_ptr_x, xtmp)
  pointer (ori_ptr_y, ytmp)
  pointer (ori_ptr_z, ztmp)

  character*200 grid_scalar
  CCTK_INT grid_scalar_len
  
  CCTK_INT  i, j, k, m, n

  CCTK_INT type_bits, state_outside

  logical docalc

  type_bits     = -1
  state_outside = -1
  
  if (do_analysis_every .le. 0) then
     return
  end if
  
  if (MOD(cctk_iteration, do_analysis_every) .ne. 0 ) then
     return
  endif

  if (excise_horizon /= 0) then

     call SpaceMask_GetTypeBits(type_bits, "mask")
     call SpaceMask_GetStateBits(state_outside, "mask", "outside")

     if (type_bits < 0) then
        call CCTK_WARN(0, "Thorn AHFinderDirect not activated, but excise_horizon requires it.")
     end if

     if (state_outside < 0) then
        call CCTK_WARN(0, "Error in obtaining AHFinderDirect GetStateBits")
     end if

  end if

  ! Coordinates of the origin used in angular momentum, quadrupole
  if (track_origin_from_grid_scalar /= 0) then
      ! Get the index of variables. This is to check if they change.
      ! If they change, there's no point in using them, just use CCTK_VarDataPtr with the name directly
      ! If they don't change, we can store them at initialization and don't bother with CCTK_STRING in Fortran
      write(*,*) 'Getting variable indices in UAv_Analysis'
      ! x
      call CCTK_FortranString(grid_scalar_len, track_origin_source_x, grid_scalar)
      call CCTK_VarIndex(gs_index_tmp_x, grid_scalar(1:grid_scalar_len))
      write(*,*) 'Index of x source: ', gs_index_tmp_x
      ! y
      call CCTK_FortranString(grid_scalar_len, track_origin_source_y, grid_scalar)
      call CCTK_VarIndex(gs_index_tmp_y, grid_scalar(1:grid_scalar_len))
      write(*,*) 'Index of y source: ', gs_index_tmp_y
      ! z
      call CCTK_FortranString(grid_scalar_len, track_origin_source_z, grid_scalar)
      call CCTK_VarIndex(gs_index_tmp_z, grid_scalar(1:grid_scalar_len))
      write(*,*) 'Index of z source: ', gs_index_tmp_z

      call CCTK_VarDataPtrI(ori_ptr_x, cctkGH, 0, gs_index_tmp_x)
      call CCTK_VarDataPtrI(ori_ptr_y, cctkGH, 0, gs_index_tmp_y)
      call CCTK_VarDataPtrI(ori_ptr_z, cctkGH, 0, gs_index_tmp_z)
      
      x0 = xtmp
      y0 = ytmp
      z0 = ztmp
  else
      x0 = origin_x
      y0 = origin_y
      z0 = origin_z
  end if

!   write(*,*) 'Tracking origin in UAv_Analysis'
!   write(*,*) 'x0 = ', x0
!   write(*,*) 'y0 = ', y0
!   write(*,*) 'z0 = ', z0

  dE_gf_volume   = 0
  dJx_gf_volume  = 0
  dJy_gf_volume  = 0
  dJz_gf_volume  = 0
  dIxx_gf_volume = 0
  dIxy_gf_volume = 0
  dIxz_gf_volume = 0
  dIyy_gf_volume = 0
  dIyz_gf_volume = 0
  dIzz_gf_volume = 0
  if (compute_density_rho == 1) then
     density_rho    = 0
  end if
  if (compute_density_p == 1) then
     density_px     = 0
     density_py     = 0
     density_pz     = 0
  end if

  do k = 1+cctk_nghostzones(3), cctk_lsh(3)-cctk_nghostzones(3)
  do j = 1+cctk_nghostzones(2), cctk_lsh(2)-cctk_nghostzones(2)
  do i = 1+cctk_nghostzones(1), cctk_lsh(1)-cctk_nghostzones(1)

    ! checking if outside the horizon, if asking for it to be excised
    docalc = .true.
    if (excise_horizon /= 0) then
       if (.not. SpaceMask_CheckStateBitsF90(space_mask, i, j, k, type_bits, state_outside)) then
          docalc = .false.
       end if
    end if

    ! if inside the horizon, no need to compute the rest (continue with the next
    ! iteration of the do loop)
    if (.not. docalc) cycle

    !--------------Get local variables ----------
    gd(1,1) = gxx(i,j,k)
    gd(1,2) = gxy(i,j,k)
    gd(1,3) = gxz(i,j,k)
    gd(2,2) = gyy(i,j,k)
    gd(2,3) = gyz(i,j,k)
    gd(3,3) = gzz(i,j,k)
    gd(2,1) = gd(1,2)
    gd(3,1) = gd(1,3)
    gd(3,2) = gd(2,3)

    alph    = alp(i,j,k)

    beta(1) = betax(i,j,k)
    beta(2) = betay(i,j,k)
    beta(3) = betaz(i,j,k)

    x1      = x(i,j,k) - x0
    y1      = y(i,j,k) - y0
    z1      = z(i,j,k) - z0

    ! stress-energy tensor variables
    Tab = 0
    if (stress_energy_state /= 0) then
       Tab(4,4) = eTtt(i,j,k)
       Tab(4,1) = eTtx(i,j,k)
       Tab(4,2) = eTty(i,j,k)
       Tab(4,3) = eTtz(i,j,k)
       Tab(1,1) = eTxx(i,j,k)
       Tab(1,2) = eTxy(i,j,k)
       Tab(1,3) = eTxz(i,j,k)
       Tab(2,2) = eTyy(i,j,k)
       Tab(2,3) = eTyz(i,j,k)
       Tab(3,3) = eTzz(i,j,k)
       Tab(1,4) = Tab(4,1)
       Tab(2,4) = Tab(4,2)
       Tab(3,4) = Tab(4,3)
       Tab(2,1) = Tab(1,2)
       Tab(3,1) = Tab(1,3)
       Tab(3,2) = Tab(2,3)
    end if
    !--------------------------------------------


    !-------------- Invert metric ---------------
    detgd =       gd(1,1) * gd(2,2) * gd(3,3)                                &
            + 2 * gd(1,2) * gd(1,3) * gd(2,3)                                &
            -     gd(1,1) * gd(2,3) ** 2                                     &
            -     gd(2,2) * gd(1,3) ** 2                                     &
            -     gd(3,3) * gd(1,2) ** 2
    gu(1,1) = (gd(2,2) * gd(3,3) - gd(2,3) ** 2     ) / detgd
    gu(2,2) = (gd(1,1) * gd(3,3) - gd(1,3) ** 2     ) / detgd
    gu(3,3) = (gd(1,1) * gd(2,2) - gd(1,2) ** 2     ) / detgd
    gu(1,2) = (gd(1,3) * gd(2,3) - gd(1,2) * gd(3,3)) / detgd
    gu(1,3) = (gd(1,2) * gd(2,3) - gd(1,3) * gd(2,2)) / detgd
    gu(2,3) = (gd(1,3) * gd(1,2) - gd(2,3) * gd(1,1)) / detgd
    gu(2,1) = gu(1,2)
    gu(3,1) = gu(1,3)
    gu(3,2) = gu(2,3)
    !--------------------------------------------

    ! Eulerian energy density
    rho = Tab(4,4)
    do m = 1, 3
       rho = rho - 2 * beta(m) * Tab(m,4)
       do n = 1, 3
          rho = rho + beta(m) * beta(n) * Tab(m,n)
       end do
    end do
    rho = rho / ( alph * alph )

    if (compute_density_rho == 1) then
      density_rho(i,j,k) = rho
    end if
    
    ! momentum density
    do n = 1, 3
      mom(n) = Tab(4,n)
      do m = 1, 3
         mom(n) = mom(n) - beta(m) * Tab(m,n)
      end do
      mom(n) = - mom(n) / alph
    end do
   
    if (compute_density_p == 1) then
      density_px(i,j,k) = mom(1)
      density_py(i,j,k) = mom(2)
      density_pz(i,j,k) = mom(3)
    end if
   

    aux = 0
    do m = 1, 3
      do n = 1, 3
        aux = aux + beta(m) * beta(n) * Tab(m,n)
      end do
    end do
    aux = (Tab(4,4) - aux) / alph

    S = 0
    do m = 1, 3
       do n = 1, 3
          S = S + gu(m,n) * Tab(m,n)
       end do
    end do


    ! dE_gf_volume = (alpha h^ij T_ij + T_tt / alpha - beta^i beta^j T_ij / alpha) sqrt(detgd)

    dE_gf_volume(i,j,k)   = (alph * S + aux) * sqrt(detgd)

    ! dJz = (-y p_x + x p_y) sqrt(detgd)        + permutations
    dJz_gf_volume(i,j,k)  = (-y1 * mom(1) + x1 * mom(2)) * sqrt(detgd)
    dJx_gf_volume(i,j,k)  = (-z1 * mom(2) + y1 * mom(3)) * sqrt(detgd)
    dJy_gf_volume(i,j,k)  = (-x1 * mom(3) + z1 * mom(1)) * sqrt(detgd)
    
    ! dI_ij = rho * x^i x^j * alpha * sqrt(detgd)
    dIxx_gf_volume(i,j,k) = alph * rho * x1 * x1 * sqrt(detgd)
    dIxy_gf_volume(i,j,k) = alph * rho * x1 * y1 * sqrt(detgd)
    dIxz_gf_volume(i,j,k) = alph * rho * x1 * z1 * sqrt(detgd)
    dIyy_gf_volume(i,j,k) = alph * rho * y1 * y1 * sqrt(detgd)
    dIyz_gf_volume(i,j,k) = alph * rho * y1 * z1 * sqrt(detgd)
    dIzz_gf_volume(i,j,k) = alph * rho * z1 * z1 * sqrt(detgd)

  end do
  end do
  end do
end subroutine UAv_Analysis_gfs


subroutine UAv_Analysis_IntegrateVol( CCTK_ARGUMENTS )
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT ierr
  CCTK_INT reduction_handle, varid

  CCTK_REAL E_int
  CCTK_REAL Jx_int, Jy_int, Jz_int
  CCTK_REAL Ixx_int, Ixy_int, Ixz_int, Iyy_int, Iyz_int, Izz_int

  if (do_analysis_every .le. 0) then
     return
  end if

  if (MOD(cctk_iteration, do_analysis_every) .ne. 0 ) then
     return
  endif
  
  call CCTK_ReductionHandle(reduction_handle, 'sum')
  if (reduction_handle < 0) then
     call CCTK_WARN(0, 'Could not obtain a handle for sum reduction')
  end if

  ! energy

  ! get index to the integration array
  call CCTK_VarIndex(varid, 'UAv_Analysis::dE_gf_volume')
  if (varid < 0) then
     call CCTK_WARN(0, 'Could not get index to grid array dE_gf_volume')
  end if

  ! do a sum over all processors
  call CCTK_Reduce(ierr, cctkGH, -1, reduction_handle, 1, CCTK_VARIABLE_REAL, &
       E_int, 1, varid)
  if (ierr < 0) then
     call CCTK_WARN(0, 'Error while reducing dE_gf_volume')
  end if

  ! TODO: is there a way of doing all components at once?

  ! Jx
  call CCTK_VarIndex(varid, 'UAv_Analysis::dJx_gf_volume')
  if (varid < 0) then
     call CCTK_WARN(0, 'Could not get index to grid array dJx_gf_volume')
  end if
  call CCTK_Reduce(ierr, cctkGH, -1, reduction_handle, 1, CCTK_VARIABLE_REAL, &
       Jx_int, 1, varid)
  if (ierr < 0) then
     call CCTK_WARN(0, 'Error while reducing dJx_gf_volume')
  end if
  
  ! Jy
  call CCTK_VarIndex(varid, 'UAv_Analysis::dJy_gf_volume')
  if (varid < 0) then
     call CCTK_WARN(0, 'Could not get index to grid array dJy_gf_volume')
  end if
  call CCTK_Reduce(ierr, cctkGH, -1, reduction_handle, 1, CCTK_VARIABLE_REAL, &
       Jy_int, 1, varid)
  if (ierr < 0) then
     call CCTK_WARN(0, 'Error while reducing dJy_gf_volume')
  end if
  
  ! Jz
  call CCTK_VarIndex(varid, 'UAv_Analysis::dJz_gf_volume')
  if (varid < 0) then
     call CCTK_WARN(0, 'Could not get index to grid array dJz_gf_volume')
  end if
  call CCTK_Reduce(ierr, cctkGH, -1, reduction_handle, 1, CCTK_VARIABLE_REAL, &
       Jz_int, 1, varid)
  if (ierr < 0) then
     call CCTK_WARN(0, 'Error while reducing dJz_gf_volume')
  end if
  
  ! Ixx
  call CCTK_VarIndex(varid, 'UAv_Analysis::dIxx_gf_volume')
  if (varid < 0) then
     call CCTK_WARN(0, 'Could not get index to grid array dIxx_gf_volume')
  end if
  call CCTK_Reduce(ierr, cctkGH, -1, reduction_handle, 1, CCTK_VARIABLE_REAL, &
       Ixx_int, 1, varid)
  if (ierr < 0) then
     call CCTK_WARN(0, 'Error while reducing dIxx_gf_volume')
  end if

  ! Ixy
  call CCTK_VarIndex(varid, 'UAv_Analysis::dIxy_gf_volume')
  if (varid < 0) then
     call CCTK_WARN(0, 'Could not get index to grid array dIxy_gf_volume')
  end if
  call CCTK_Reduce(ierr, cctkGH, -1, reduction_handle, 1, CCTK_VARIABLE_REAL, &
       Ixy_int, 1, varid)
  if (ierr < 0) then
     call CCTK_WARN(0, 'Error while reducing dIxy_gf_volume')
  end if

  ! Ixz
  call CCTK_VarIndex(varid, 'UAv_Analysis::dIxz_gf_volume')
  if (varid < 0) then
     call CCTK_WARN(0, 'Could not get index to grid array dIxz_gf_volume')
  end if
  call CCTK_Reduce(ierr, cctkGH, -1, reduction_handle, 1, CCTK_VARIABLE_REAL, &
       Ixz_int, 1, varid)
  if (ierr < 0) then
     call CCTK_WARN(0, 'Error while reducing dIxz_gf_volume')
  end if

  ! Iyy
  call CCTK_VarIndex(varid, 'UAv_Analysis::dIyy_gf_volume')
  if (varid < 0) then
     call CCTK_WARN(0, 'Could not get index to grid array dIyy_gf_volume')
  end if
  call CCTK_Reduce(ierr, cctkGH, -1, reduction_handle, 1, CCTK_VARIABLE_REAL, &
       Iyy_int, 1, varid)
  if (ierr < 0) then
     call CCTK_WARN(0, 'Error while reducing dIyy_gf_volume')
  end if

  ! Iyz
  call CCTK_VarIndex(varid, 'UAv_Analysis::dIyz_gf_volume')
  if (varid < 0) then
     call CCTK_WARN(0, 'Could not get index to grid array dIyz_gf_volume')
  end if
  call CCTK_Reduce(ierr, cctkGH, -1, reduction_handle, 1, CCTK_VARIABLE_REAL, &
       Iyz_int, 1, varid)
  if (ierr < 0) then
     call CCTK_WARN(0, 'Error while reducing dIyz_gf_volume')
  end if

  ! Izz
  call CCTK_VarIndex(varid, 'UAv_Analysis::dIzz_gf_volume')
  if (varid < 0) then
     call CCTK_WARN(0, 'Could not get index to grid array dIzz_gf_volume')
  end if
  call CCTK_Reduce(ierr, cctkGH, -1, reduction_handle, 1, CCTK_VARIABLE_REAL, &
       Izz_int, 1, varid)
  if (ierr < 0) then
     call CCTK_WARN(0, 'Error while reducing dIzz_gf_volume')
  end if


  ! the multiplication with the volume element needs to be done here
  total_energy = E_int * cctk_delta_space(1) * cctk_delta_space(2) * cctk_delta_space(3)
  
  total_angular_momentum_x = Jx_int * cctk_delta_space(1) * cctk_delta_space(2) * cctk_delta_space(3)
  total_angular_momentum_y = Jy_int * cctk_delta_space(1) * cctk_delta_space(2) * cctk_delta_space(3)
  total_angular_momentum_z = Jz_int * cctk_delta_space(1) * cctk_delta_space(2) * cctk_delta_space(3)

  Ixx = Ixx_int * cctk_delta_space(1) * cctk_delta_space(2) * cctk_delta_space(3)
  Ixy = Ixy_int * cctk_delta_space(1) * cctk_delta_space(2) * cctk_delta_space(3)
  Ixz = Ixz_int * cctk_delta_space(1) * cctk_delta_space(2) * cctk_delta_space(3)
  Iyy = Iyy_int * cctk_delta_space(1) * cctk_delta_space(2) * cctk_delta_space(3)
  Iyz = Iyz_int * cctk_delta_space(1) * cctk_delta_space(2) * cctk_delta_space(3)
  Izz = Izz_int * cctk_delta_space(1) * cctk_delta_space(2) * cctk_delta_space(3)

  ! write(*,*) 'total_energy = ', total_energy

  ! note that the integrated values just obtained do *not* take into account
  ! any grid symmetries that may be present. so one needs to remember to
  ! multiply by the appropriate factor (ie, 8 in the case of using just an
  ! octant, etc) a posteriori.

end subroutine UAv_Analysis_IntegrateVol
