!> Accelerations due to the Coriolis force and momentum advection
module MOM_CoriolisAdv

! This file is part of MOM6. See LICENSE.md for the license.

!> \author Robert Hallberg, April 1994 - June 2002

use MOM_diag_mediator, only : post_data, query_averaging_enabled, diag_ctrl
use MOM_diag_mediator, only : post_product_u, post_product_sum_u
use MOM_diag_mediator, only : post_product_v, post_product_sum_v
use MOM_diag_mediator, only : register_diag_field, safe_alloc_ptr, time_type
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_open_boundary, only : ocean_OBC_type, OBC_DIRECTION_E, OBC_DIRECTION_W
use MOM_open_boundary, only : OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_string_functions, only : uppercase
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : accel_diag_ptrs, porous_barrier_type
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_wave_interface, only : wave_parameters_CS

implicit none ; private

public CorAdCalc, CoriolisAdv_init, CoriolisAdv_end

#include <MOM_memory.h>

!> Control structure for mom_coriolisadv
type, public :: CoriolisAdv_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  integer :: Coriolis_Scheme !< Selects the discretization for the Coriolis terms.
                             !! Valid values are:
                             !! - SADOURNY75_ENERGY - Sadourny, 1975
                             !! - ARAKAWA_HSU90     - Arakawa & Hsu, 1990, Energy & non-div. Enstrophy
                             !! - ROBUST_ENSTRO     - Pseudo-enstrophy scheme
                             !! - SADOURNY75_ENSTRO - Sadourny, JAS 1975, Enstrophy
                             !! - ARAKAWA_LAMB81    - Arakawa & Lamb, MWR 1981, Energy & Enstrophy
                             !! - ARAKAWA_LAMB_BLEND - A blend of Arakawa & Lamb with Arakawa & Hsu and Sadourny energy.
                             !! The default, SADOURNY75_ENERGY, is the safest choice then the
                             !! deformation radius is poorly resolved.
  integer :: KE_Scheme       !< KE_SCHEME selects the discretization for
                             !! the kinetic energy. Valid values are:
                             !!  KE_ARAKAWA, KE_SIMPLE_GUDONOV, KE_GUDONOV
  integer :: PV_Adv_Scheme   !< PV_ADV_SCHEME selects the discretization for PV advection
                             !! Valid values are:
                             !! - PV_ADV_CENTERED - centered (aka Sadourny, 75)
                             !! - PV_ADV_UPWIND1  - upwind, first order
  real    :: F_eff_max_blend !< The factor by which the maximum effective Coriolis
                             !! acceleration from any point can be increased when
                             !! blending different discretizations with the
                             !! ARAKAWA_LAMB_BLEND Coriolis scheme [nondim].
                             !! This must be greater than 2.0, and is 4.0 by default.
  real    :: wt_lin_blend    !< A weighting value beyond which the blending between
                             !! Sadourny and Arakawa & Hsu goes linearly to 0 [nondim].
                             !! This must be between 1 and 1e-15, often 1/8.
  logical :: no_slip         !< If true, no slip boundary conditions are used.
                             !! Otherwise free slip boundary conditions are assumed.
                             !! The implementation of the free slip boundary
                             !! conditions on a C-grid is much cleaner than the
                             !! no slip boundary conditions. The use of free slip
                             !! b.c.s is strongly encouraged. The no slip b.c.s
                             !! are not implemented with the biharmonic viscosity.
  logical :: bound_Coriolis  !< If true, the Coriolis terms at u points are
                             !! bounded by the four estimates of (f+rv)v from the
                             !! four neighboring v points, and similarly at v
                             !! points.  This option would have no effect on the
                             !! SADOURNY75_ENERGY scheme if it were possible to
                             !! use centered difference thickness fluxes.
  logical :: Coriolis_En_Dis !< If CORIOLIS_EN_DIS is defined, two estimates of
                             !! the thickness fluxes are used to estimate the
                             !! Coriolis term, and the one that dissipates energy
                             !! relative to the other one is used.  This is only
                             !! available at present if Coriolis scheme is
                             !! SADOURNY75_ENERGY.
  logical :: weno_velocity_smooth !! If true, use velocity to compute the weighting for WENO
  logical :: UP3_use_limiter !! If true, use flux limiter when UP3 scheme is called
  type(time_type), pointer :: Time !< A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the timing of diagnostic output.
  !>@{ Diagnostic IDs
  integer :: id_rv = -1, id_PV = -1, id_gKEu = -1, id_gKEv = -1
  integer :: id_rvxu = -1, id_rvxv = -1
  ! integer :: id_hf_gKEu    = -1, id_hf_gKEv    = -1
  integer :: id_hf_gKEu_2d = -1, id_hf_gKEv_2d = -1
  integer :: id_intz_gKEu_2d = -1, id_intz_gKEv_2d = -1
  ! integer :: id_hf_rvxu    = -1, id_hf_rvxv    = -1
  integer :: id_hf_rvxu_2d = -1, id_hf_rvxv_2d = -1
  integer :: id_h_gKEu = -1, id_h_gKEv = -1
  integer :: id_h_rvxu = -1, id_h_rvxv = -1
  integer :: id_intz_rvxu_2d = -1, id_intz_rvxv_2d = -1
  integer :: id_CAuS = -1, id_CAvS = -1
  !>@}
end type CoriolisAdv_CS

!>@{ Enumeration values for Coriolis_Scheme
integer, parameter :: SADOURNY75_ENERGY = 1
integer, parameter :: ARAKAWA_HSU90     = 2
integer, parameter :: ROBUST_ENSTRO     = 3
integer, parameter :: SADOURNY75_ENSTRO = 4
integer, parameter :: ARAKAWA_LAMB81    = 5
integer, parameter :: AL_BLEND          = 6
integer, parameter :: SECOND_ENSTRO = 7
integer, parameter :: UP3_ENSTRO = 8
integer, parameter :: wenovi_7th_ENSTRO = 9
integer, parameter :: wenovi_5th_ENSTRO = 15
integer, parameter :: wenovi_7th_split = 16
integer, parameter :: UP3_split = 17
integer, parameter :: UP3_PV_ENSTRO = 18
integer, parameter :: CEN4_ENSTRO = 19
character*(20), parameter :: SADOURNY75_ENERGY_STRING = "SADOURNY75_ENERGY"
character*(20), parameter :: ARAKAWA_HSU_STRING = "ARAKAWA_HSU90"
character*(20), parameter :: ROBUST_ENSTRO_STRING = "ROBUST_ENSTRO"
character*(20), parameter :: SADOURNY75_ENSTRO_STRING = "SADOURNY75_ENSTRO"
character*(20), parameter :: ARAKAWA_LAMB_STRING = "ARAKAWA_LAMB81"
character*(20), parameter :: AL_BLEND_STRING = "ARAKAWA_LAMB_BLEND"
character*(20), parameter :: SECOND_ENSTRO_STRING = "SECOND_ENSTRO"
character*(20), parameter :: UP3_ENSTRO_STRING = "UP3_ENSTRO"
character*(20), parameter :: UP3_SPLIT_STRING = "UP3_SPLIT"
character*(20), parameter :: UP3_PV_ENSTRO_STRING = "UP3_PV_ENSTRO"
character*(20), parameter :: WENOVI_5TH_ENSTRO_STRING = "WENOVI_5TH_ENSTRO"
character*(20), parameter :: WENOVI_7TH_ENSTRO_STRING = "WENOVI_7TH_ENSTRO"
character*(20), parameter :: WENOVI_7TH_SPLIT_STRING = "WENOVI_7TH_SPLIT"
character*(20), parameter :: CEN4_ENSTRO_STRING = "CEN4_ENSTRO"
!>@}
!>@{ Enumeration values for KE_Scheme
integer, parameter :: KE_ARAKAWA        = 10
integer, parameter :: KE_SIMPLE_GUDONOV = 11
integer, parameter :: KE_GUDONOV        = 12
integer, parameter :: KE_UP3            = 13
integer, parameter :: KE_wenovi_7th     = 14
character*(20), parameter :: KE_ARAKAWA_STRING = "KE_ARAKAWA"
character*(20), parameter :: KE_SIMPLE_GUDONOV_STRING = "KE_SIMPLE_GUDONOV"
character*(20), parameter :: KE_GUDONOV_STRING = "KE_GUDONOV"
character*(20), parameter :: KE_UP3_STRING = "KE_UP3"
character*(20), parameter :: KE_WENOVI_7TH_STRING = "KE_WENOVI_7TH"
!>@}
!>@{ Enumeration values for PV_Adv_Scheme
integer, parameter :: PV_ADV_CENTERED   = 21
integer, parameter :: PV_ADV_UPWIND1    = 22
character*(20), parameter :: PV_ADV_CENTERED_STRING = "PV_ADV_CENTERED"
character*(20), parameter :: PV_ADV_UPWIND1_STRING = "PV_ADV_UPWIND1"
!>@}

contains

!> Calculates the Coriolis and momentum advection contributions to the acceleration.
subroutine CorAdCalc(u, v, h, uh, vh, CAu, CAv, OBC, AD, G, GV, US, CS, pbv, Waves)
  type(ocean_grid_type),                      intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV !< Vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)    :: u  !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)    :: v  !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)    :: uh !< Zonal transport u*h*dy
                                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)    :: vh !< Meridional transport v*h*dx
                                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(out)   :: CAu !< Zonal acceleration due to Coriolis
                                                                  !! and momentum advection [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(out)   :: CAv !< Meridional acceleration due to Coriolis
                                                                  !! and momentum advection [L T-2 ~> m s-2].
  type(ocean_OBC_type),                       pointer       :: OBC !< Open boundary control structure
  type(accel_diag_ptrs),                      intent(inout) :: AD  !< Storage for acceleration diagnostics
  type(unit_scale_type),                      intent(in)    :: US  !< A dimensional unit scaling type
  type(CoriolisAdv_CS),                       intent(in)    :: CS  !< Control structure for MOM_CoriolisAdv
  type(porous_barrier_type),                  intent(in)    :: pbv !< porous barrier fractional cell metrics
  type(Wave_parameters_CS),         optional, pointer       :: Waves !< An optional pointer to Stokes drift CS

  ! Local variables
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    q, &        ! Layer potential vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1].
    qS, &       ! Layer Stokes vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1].
    Ih_q, &     ! The inverse of thickness interpolated to q points [H-1 ~> m-1 or m2 kg-1].
    Area_q      ! The sum of the ocean areas at the 4 adjacent thickness points [L2 ~> m2].

  real, dimension(SZIB_(G),SZJ_(G)) :: &
    a, b, c, d  ! a, b, c, & d are combinations of the potential vorticities
                ! surrounding an h grid point.  At small scales, a = q/4,
                ! b = q/4, etc.  All are in [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1],
                ! and use the indexing of the corresponding u point.

  real, dimension(SZI_(G),SZJ_(G)) :: &
    Area_h, &   ! The ocean area at h points [L2 ~> m2].  Area_h is used to find the
                ! average thickness in the denominator of q.  0 for land points.
    KE          ! Kinetic energy per unit mass [L2 T-2 ~> m2 s-2], KE = (u^2 + v^2)/2.
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    hArea_u, &  ! The cell area weighted thickness interpolated to u points
                ! times the effective areas [H L2 ~> m3 or kg].
    KEx, &      ! The zonal gradient of Kinetic energy per unit mass [L T-2 ~> m s-2],
                ! KEx = d/dx KE.
    uh_center   ! Transport based on arithmetic mean h at u-points [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJB_(G)) :: &
    hArea_v, &  ! The cell area weighted thickness interpolated to v points
                ! times the effective areas [H L2 ~> m3 or kg].
    KEy, &      ! The meridional gradient of Kinetic energy per unit mass [L T-2 ~> m s-2],
                ! KEy = d/dy KE.
    vh_center   ! Transport based on arithmetic mean h at v-points [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJ_(G)) :: &
    uh_min, uh_max, &   ! The smallest and largest estimates of the zonal volume fluxes through
                        ! the faces (i.e. u*h*dy) [H L2 T-1 ~> m3 s-1 or kg s-1]
    vh_min, vh_max, &   ! The smallest and largest estimates of the meridional volume fluxes through
                        ! the faces (i.e. v*h*dx) [H L2 T-1 ~> m3 s-1 or kg s-1]
    ep_u, ep_v  ! Additional pseudo-Coriolis terms in the Arakawa and Lamb
                ! discretization [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1].
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    dvdx, dudy, & ! Contributions to the circulation around q-points [L2 T-1 ~> m2 s-1]
    dvSdx, duSdy, & ! idem. for Stokes drift [L2 T-1 ~> m2 s-1]
    rel_vort, & ! Relative vorticity at q-points [T-1 ~> s-1].
    abs_vort, & ! Absolute vorticity at q-points [T-1 ~> s-1].
    stk_vort, & ! Stokes vorticity at q-points [T-1 ~> s-1].
    q2          ! Relative vorticity over thickness [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1].
  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)) :: &
    PV, &       ! A diagnostic array of the potential vorticities [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1].
    RV          ! A diagnostic array of the relative vorticities [T-1 ~> s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: CAuS ! Stokes contribution to CAu [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: CAvS ! Stokes contribution to CAv [L T-2 ~> m s-2]
  real :: fv1, fv2, fv3, fv4   ! (f+rv)*v at the 4 points surrounding a u points[L T-2 ~> m s-2]
  real :: fu1, fu2, fu3, fu4   ! -(f+rv)*u at the 4 points surrounding a v point [L T-2 ~> m s-2]
  real :: max_fv, max_fu       ! The maximum of the neighboring Coriolis accelerations [L T-2 ~> m s-2]
  real :: min_fv, min_fu       ! The minimum of the neighboring Coriolis accelerations [L T-2 ~> m s-2]

  real, parameter :: C1_12 = 1.0 / 12.0 ! C1_12 = 1/12 [nondim]
  real, parameter :: C1_24 = 1.0 / 24.0 ! C1_24 = 1/24 [nondim]
  real :: max_Ihq, min_Ihq       ! The maximum and minimum of the nearby Ihq [H-1 ~> m-1 or m2 kg-1].
  real :: hArea_q                ! The sum of area times thickness of the cells
                                 ! surrounding a q point [H L2 ~> m3 or kg].
  real :: vol_neglect            ! A volume so small that is expected to be
                                 ! lost in roundoff [H L2 ~> m3 or kg].
  real :: temp1, temp2           ! Temporary variables [L2 T-2 ~> m2 s-2].
  real :: eps_vel                ! A tiny, positive velocity [L T-1 ~> m s-1].

  real :: uhc, vhc               ! Centered estimates of uh and vh [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: uhm, vhm               ! The input estimates of uh and vh [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: c1, c2, c3, slope      ! Nondimensional parameters for the Coriolis limiter scheme [nondim]

  real :: Fe_m2         ! Temporary variable associated with the ARAKAWA_LAMB_BLEND scheme [nondim]
  real :: rat_lin       ! Temporary variable associated with the ARAKAWA_LAMB_BLEND scheme [nondim]
  real :: rat_m1        ! The ratio of the maximum neighboring inverse thickness
                        ! to the minimum inverse thickness minus 1 [nondim]. rat_m1 >= 0.
  real :: AL_wt         ! The relative weight of the Arakawa & Lamb scheme to the
                        ! Arakawa & Hsu scheme [nondim], between 0 and 1.
  real :: Sad_wt        ! The relative weight of the Sadourny energy scheme to
                        ! the other two with the ARAKAWA_LAMB_BLEND scheme [nondim],
                        ! between 0 and 1.

  real :: Heff1, Heff2  ! Temporary effective H at U or V points [H ~> m or kg m-2].
  real :: Heff3, Heff4  ! Temporary effective H at U or V points [H ~> m or kg m-2].
  real :: h_tiny        ! A very small thickness [H ~> m or kg m-2].
  real :: UHeff, VHeff  ! More temporary variables [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: QUHeff,QVHeff ! More temporary variables [H L2 T-2 ~> m3 s-2 or kg s-2].
  integer :: i, j, k, n, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: seventh_order, fifth_order, third_order, second_order ! Order of accuracy for the WENO calculations
  logical :: Stokes_VF
  real :: u_v, v_u, q_v, q_u ! u_v is the u velocity at v point, v_u is the v velocity at u point 
  real :: u_q1, u_q2, u_q3, u_q4, u_q5, u_q6, u_q7, u_q8 ! is the v velocity at q points
  real :: v_q1, v_q2, v_q3, v_q4, v_q5, v_q6, v_q7, v_q8 ! is the u velocity at q points
  real :: abs_vort_u, abs_vort_v  ! absolute vorticity at u and v points
  real :: rel_vort_u, rel_vort_v  ! relative vorticity at u and v points
  real :: f_u, f_v      ! Coriolis coefficient at u and v points
  real :: fv, fu        ! f*v at u point and -f*u at v point
  real :: fq1, fq2, fq3, fq4 !f/h at u point using Arakawa and Hsu reconstruction
!  real :: theta, psi    ! temperory variables for UP3 limiter

! To work, the following fields must be set outside of the usual
! is to ie range before this subroutine is called:
!   v(is-1:ie+2,js-1:je+1), u(is-1:ie+1,js-1:je+2), h(is-1:ie+2,js-1:je+2),
!   uh(is-1,ie,js:je+1) and vh(is:ie+1,js-1:je).

  if (.not.CS%initialized) call MOM_error(FATAL, &
         "MOM_CoriolisAdv: Module must be initialized before it is used.")

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB ; nz = GV%ke
  vol_neglect = GV%H_subroundoff * (1e-4 * US%m_to_L)**2
  eps_vel = 1.0e-10*US%m_s_to_L_T
  h_tiny = GV%Angstrom_H  ! Perhaps this should be set to h_neglect instead.

  !$OMP parallel do default(private) shared(Isq,Ieq,Jsq,Jeq,G,Area_h)
  do j=Jsq-1,Jeq+2 ; do I=Isq-1,Ieq+2
    Area_h(i,j) = G%mask2dT(i,j) * G%areaT(i,j)
  enddo ; enddo
  if (associated(OBC)) then ; do n=1,OBC%number_of_segments
    if (.not. OBC%segment(n)%on_pe) cycle
    I = OBC%segment(n)%HI%IsdB ; J = OBC%segment(n)%HI%JsdB
    if (OBC%segment(n)%is_N_or_S .and. (J >= Jsq-1) .and. (J <= Jeq+1)) then
      do i = max(Isq-1,OBC%segment(n)%HI%isd), min(Ieq+2,OBC%segment(n)%HI%ied)
        if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
          Area_h(i,j+1) = Area_h(i,j)
        else ! (OBC%segment(n)%direction == OBC_DIRECTION_S)
          Area_h(i,j) = Area_h(i,j+1)
        endif
      enddo
    elseif (OBC%segment(n)%is_E_or_W .and. (I >= Isq-1) .and. (I <= Ieq+1)) then
      do j = max(Jsq-1,OBC%segment(n)%HI%jsd), min(Jeq+2,OBC%segment(n)%HI%jed)
        if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
          Area_h(i+1,j) = Area_h(i,j)
        else ! (OBC%segment(n)%direction == OBC_DIRECTION_W)
          Area_h(i,j) = Area_h(i+1,j)
        endif
      enddo
    endif
  enddo ; endif
  !$OMP parallel do default(private) shared(Isq,Ieq,Jsq,Jeq,G,Area_h,Area_q)
  do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
    Area_q(i,j) = (Area_h(i,j) + Area_h(i+1,j+1)) + &
                  (Area_h(i+1,j) + Area_h(i,j+1))
  enddo ; enddo

  Stokes_VF = .false.
  if (present(Waves)) then ; if (associated(Waves)) then
    Stokes_VF = Waves%Stokes_VF
  endif ; endif

  !$OMP parallel do default(private) shared(u,v,h,uh,vh,CAu,CAv,G,GV,CS,AD,Area_h,Area_q,&
  !$OMP                        RV,PV,is,ie,js,je,Isq,Ieq,Jsq,Jeq,nz,vol_neglect,h_tiny,OBC,eps_vel, &
  !$OMP                        pbv, Stokes_VF)
  do k=1,nz

    ! Here the second order accurate layer potential vorticities, q,
    ! are calculated.  hq is  second order accurate in space.  Relative
    ! vorticity is second order accurate everywhere with free slip b.c.s,
    ! but only first order accurate at boundaries with no slip b.c.s.
    ! First calculate the contributions to the circulation around the q-point.
    if (Stokes_VF) then
      if (CS%id_CAuS>0 .or. CS%id_CAvS>0) then
        do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
          dvSdx(I,J) = (-Waves%us_y(i+1,J,k)*G%dyCv(i+1,J)) - &
                       (-Waves%us_y(i,J,k)*G%dyCv(i,J))
          duSdy(I,J) = (-Waves%us_x(I,j+1,k)*G%dxCu(I,j+1)) - &
                       (-Waves%us_x(I,j,k)*G%dxCu(I,j))
        enddo; enddo
      endif
      if (.not. Waves%Passive_Stokes_VF) then
        do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
          dvdx(I,J) = ((v(i+1,J,k)-Waves%us_y(i+1,J,k))*G%dyCv(i+1,J)) - &
                      ((v(i,J,k)-Waves%us_y(i,J,k))*G%dyCv(i,J))
          dudy(I,J) = ((u(I,j+1,k)-Waves%us_x(I,j+1,k))*G%dxCu(I,j+1)) - &
                      ((u(I,j,k)-Waves%us_x(I,j,k))*G%dxCu(I,j))
        enddo; enddo
      else
        do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
          dvdx(I,J) = (v(i+1,J,k)*G%dyCv(i+1,J)) - (v(i,J,k)*G%dyCv(i,J))
          dudy(I,J) = (u(I,j+1,k)*G%dxCu(I,j+1)) - (u(I,j,k)*G%dxCu(I,j))
        enddo; enddo
      endif
    else
      do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
        dvdx(I,J) = (v(i+1,J,k)*G%dyCv(i+1,J)) - (v(i,J,k)*G%dyCv(i,J))
        dudy(I,J) = (u(I,j+1,k)*G%dxCu(I,j+1)) - (u(I,j,k)*G%dxCu(I,j))
      enddo; enddo
    endif
    do J=Jsq-1,Jeq+1 ; do i=Isq-1,Ieq+2
      hArea_v(i,J) = 0.5*((Area_h(i,j) * h(i,j,k)) + (Area_h(i,j+1) * h(i,j+1,k)))
    enddo ; enddo
    do j=Jsq-1,Jeq+2 ; do I=Isq-1,Ieq+1
      hArea_u(I,j) = 0.5*((Area_h(i,j) * h(i,j,k)) + (Area_h(i+1,j) * h(i+1,j,k)))
    enddo ; enddo

    if (CS%Coriolis_En_Dis) then
      do j=Jsq,Jeq+1 ; do I=is-1,ie
        uh_center(I,j) = 0.5 * ((G%dy_Cu(I,j)*pbv%por_face_areaU(I,j,k)) * u(I,j,k)) * (h(i,j,k) + h(i+1,j,k))
      enddo ; enddo
      do J=js-1,je ; do i=Isq,Ieq+1
        vh_center(i,J) = 0.5 * ((G%dx_Cv(i,J)*pbv%por_face_areaV(i,J,k)) * v(i,J,k)) * (h(i,j,k) + h(i,j+1,k))
      enddo ; enddo
    endif

    ! Adjust circulation components to relative vorticity and thickness projected onto
    ! velocity points on open boundaries.
    if (associated(OBC)) then ; do n=1,OBC%number_of_segments
      if (.not. OBC%segment(n)%on_pe) cycle
      I = OBC%segment(n)%HI%IsdB ; J = OBC%segment(n)%HI%JsdB
      if (OBC%segment(n)%is_N_or_S .and. (J >= Jsq-1) .and. (J <= Jeq+1)) then
        if (OBC%zero_vorticity) then ; do I=OBC%segment(n)%HI%IsdB,OBC%segment(n)%HI%IedB
          dvdx(I,J) = 0. ; dudy(I,J) = 0.
        enddo ; endif
        if (OBC%freeslip_vorticity) then ; do I=OBC%segment(n)%HI%IsdB,OBC%segment(n)%HI%IedB
          dudy(I,J) = 0.
        enddo ; endif
        if (OBC%computed_vorticity) then ; do I=OBC%segment(n)%HI%IsdB,OBC%segment(n)%HI%IedB
          if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
            dudy(I,J) = 2.0*(OBC%segment(n)%tangential_vel(I,J,k) - u(I,j,k))*G%dxCu(I,j)
          else ! (OBC%segment(n)%direction == OBC_DIRECTION_S)
            dudy(I,J) = 2.0*(u(I,j+1,k) - OBC%segment(n)%tangential_vel(I,J,k))*G%dxCu(I,j+1)
          endif
        enddo ; endif
        if (OBC%specified_vorticity) then ; do I=OBC%segment(n)%HI%IsdB,OBC%segment(n)%HI%IedB
          if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
            dudy(I,J) = OBC%segment(n)%tangential_grad(I,J,k)*G%dxCu(I,j)*G%dyBu(I,J)
          else ! (OBC%segment(n)%direction == OBC_DIRECTION_S)
            dudy(I,J) = OBC%segment(n)%tangential_grad(I,J,k)*G%dxCu(I,j+1)*G%dyBu(I,J)
          endif
        enddo ; endif

        ! Project thicknesses across OBC points with a no-gradient condition.
        do i = max(Isq-1,OBC%segment(n)%HI%isd), min(Ieq+2,OBC%segment(n)%HI%ied)
          if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
            hArea_v(i,J) = 0.5 * (Area_h(i,j) + Area_h(i,j+1)) * h(i,j,k)
          else ! (OBC%segment(n)%direction == OBC_DIRECTION_S)
            hArea_v(i,J) = 0.5 * (Area_h(i,j) + Area_h(i,j+1)) * h(i,j+1,k)
          endif
        enddo

        if (CS%Coriolis_En_Dis) then
          do i = max(Isq-1,OBC%segment(n)%HI%isd), min(Ieq+2,OBC%segment(n)%HI%ied)
            if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
              vh_center(i,J) = (G%dx_Cv(i,J)*pbv%por_face_areaV(i,J,k)) * v(i,J,k) * h(i,j,k)
            else ! (OBC%segment(n)%direction == OBC_DIRECTION_S)
              vh_center(i,J) = (G%dx_Cv(i,J)*pbv%por_face_areaV(i,J,k)) * v(i,J,k) * h(i,j+1,k)
            endif
          enddo
        endif
      elseif (OBC%segment(n)%is_E_or_W .and. (I >= Isq-1) .and. (I <= Ieq+1)) then
        if (OBC%zero_vorticity) then ; do J=OBC%segment(n)%HI%JsdB,OBC%segment(n)%HI%JedB
          dvdx(I,J) = 0. ; dudy(I,J) = 0.
        enddo ; endif
        if (OBC%freeslip_vorticity) then ; do J=OBC%segment(n)%HI%JsdB,OBC%segment(n)%HI%JedB
          dvdx(I,J) = 0.
        enddo ; endif
        if (OBC%computed_vorticity) then ; do J=OBC%segment(n)%HI%JsdB,OBC%segment(n)%HI%JedB
          if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
            dvdx(I,J) = 2.0*(OBC%segment(n)%tangential_vel(I,J,k) - v(i,J,k))*G%dyCv(i,J)
          else ! (OBC%segment(n)%direction == OBC_DIRECTION_W)
            dvdx(I,J) = 2.0*(v(i+1,J,k) - OBC%segment(n)%tangential_vel(I,J,k))*G%dyCv(i+1,J)
          endif
        enddo ; endif
        if (OBC%specified_vorticity) then ; do J=OBC%segment(n)%HI%JsdB,OBC%segment(n)%HI%JedB
          if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
            dvdx(I,J) = OBC%segment(n)%tangential_grad(I,J,k)*G%dyCv(i,J)*G%dxBu(I,J)
          else ! (OBC%segment(n)%direction == OBC_DIRECTION_W)
            dvdx(I,J) = OBC%segment(n)%tangential_grad(I,J,k)*G%dyCv(i+1,J)*G%dxBu(I,J)
          endif
        enddo ; endif

        ! Project thicknesses across OBC points with a no-gradient condition.
        do j = max(Jsq-1,OBC%segment(n)%HI%jsd), min(Jeq+2,OBC%segment(n)%HI%jed)
          if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
            hArea_u(I,j) = 0.5*(Area_h(i,j) + Area_h(i+1,j)) * h(i,j,k)
          else ! (OBC%segment(n)%direction == OBC_DIRECTION_W)
            hArea_u(I,j) = 0.5*(Area_h(i,j) + Area_h(i+1,j)) * h(i+1,j,k)
          endif
        enddo
        if (CS%Coriolis_En_Dis) then
          do j = max(Jsq-1,OBC%segment(n)%HI%jsd), min(Jeq+2,OBC%segment(n)%HI%jed)
            if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
              uh_center(I,j) = (G%dy_Cu(I,j)*pbv%por_face_areaU(I,j,k)) * u(I,j,k) * h(i,j,k)
            else ! (OBC%segment(n)%direction == OBC_DIRECTION_W)
              uh_center(I,j) = (G%dy_Cu(I,j)*pbv%por_face_areaU(I,j,k)) * u(I,j,k) * h(i+1,j,k)
            endif
          enddo
        endif
      endif
    enddo ; endif

    if (associated(OBC)) then ; do n=1,OBC%number_of_segments
      if (.not. OBC%segment(n)%on_pe) cycle
      ! Now project thicknesses across cell-corner points in the OBCs.  The two
      ! projections have to occur in sequence and can not be combined easily.
      I = OBC%segment(n)%HI%IsdB ; J = OBC%segment(n)%HI%JsdB
      if (OBC%segment(n)%is_N_or_S .and. (J >= Jsq-1) .and. (J <= Jeq+1)) then
        do I = max(Isq-1,OBC%segment(n)%HI%IsdB), min(Ieq+1,OBC%segment(n)%HI%IedB)
          if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
            if (Area_h(i,j) + Area_h(i+1,j) > 0.0) then
              hArea_u(I,j+1) = hArea_u(I,j) * ((Area_h(i,j+1) + Area_h(i+1,j+1)) / &
                                               (Area_h(i,j) + Area_h(i+1,j)))
            else ; hArea_u(I,j+1) = 0.0 ; endif
          else ! (OBC%segment(n)%direction == OBC_DIRECTION_S)
            if (Area_h(i,j+1) + Area_h(i+1,j+1) > 0.0) then
              hArea_u(I,j) = hArea_u(I,j+1) * ((Area_h(i,j) + Area_h(i+1,j)) / &
                                               (Area_h(i,j+1) + Area_h(i+1,j+1)))
            else ; hArea_u(I,j) = 0.0 ; endif
          endif
        enddo
      elseif (OBC%segment(n)%is_E_or_W .and. (I >= Isq-1) .and. (I <= Ieq+1)) then
        do J = max(Jsq-1,OBC%segment(n)%HI%JsdB), min(Jeq+1,OBC%segment(n)%HI%JedB)
          if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
            if (Area_h(i,j) + Area_h(i,j+1) > 0.0) then
              hArea_v(i+1,J) = hArea_v(i,J) * ((Area_h(i+1,j) + Area_h(i+1,j+1)) / &
                                               (Area_h(i,j) + Area_h(i,j+1)))
            else ; hArea_v(i+1,J) = 0.0 ; endif
          else ! (OBC%segment(n)%direction == OBC_DIRECTION_W)
            hArea_v(i,J) = 0.5 * (Area_h(i,j) + Area_h(i,j+1)) * h(i,j+1,k)
            if (Area_h(i+1,j) + Area_h(i+1,j+1) > 0.0) then
              hArea_v(i,J) = hArea_v(i+1,J) * ((Area_h(i,j) + Area_h(i,j+1)) / &
                                               (Area_h(i+1,j) + Area_h(i+1,j+1)))
            else ; hArea_v(i,J) = 0.0 ; endif
          endif
        enddo
      endif
    enddo ; endif

    if (CS%no_slip) then
      do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
        rel_vort(I,J) = (2.0 - G%mask2dBu(I,J)) * (dvdx(I,J) - dudy(I,J)) * G%IareaBu(I,J)
      enddo; enddo
      if (Stokes_VF) then
        if (CS%id_CAuS>0 .or. CS%id_CAvS>0) then
          do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
            stk_vort(I,J) = (2.0 - G%mask2dBu(I,J)) * (dvSdx(I,J) - duSdy(I,J)) * G%IareaBu(I,J)
          enddo; enddo
        endif
      endif
    else
      do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
        rel_vort(I,J) = G%mask2dBu(I,J) * (dvdx(I,J) - dudy(I,J)) * G%IareaBu(I,J)
      enddo; enddo
      if (Stokes_VF) then
        if (CS%id_CAuS>0 .or. CS%id_CAvS>0) then
          do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
            stk_vort(I,J) = (2.0 - G%mask2dBu(I,J)) * (dvSdx(I,J) - duSdy(I,J)) * G%IareaBu(I,J)
          enddo; enddo
        endif
      endif
    endif

    do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
      abs_vort(I,J) = G%CoriolisBu(I,J) + rel_vort(I,J)
    enddo ; enddo

    do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
      hArea_q = (hArea_u(I,j) + hArea_u(I,j+1)) + (hArea_v(i,J) + hArea_v(i+1,J))
      Ih_q(I,J) = Area_q(I,J) / (hArea_q + vol_neglect)
      q(I,J) = abs_vort(I,J) * Ih_q(I,J)
    enddo; enddo

    if (Stokes_VF) then
      if (CS%id_CAuS>0 .or. CS%id_CAvS>0) then
        do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
          qS(I,J) = stk_vort(I,J) * Ih_q(I,J)
        enddo; enddo
      endif
    endif

    if (CS%id_rv > 0) then
      do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
        RV(I,J,k) = rel_vort(I,J)
      enddo ; enddo
    endif

    if (CS%id_PV > 0) then
      do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
        PV(I,J,k) = q(I,J)
      enddo ; enddo
    endif

    if (associated(AD%rv_x_v) .or. associated(AD%rv_x_u)) then
      do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
        q2(I,J) = rel_vort(I,J) * Ih_q(I,J)
      enddo ; enddo
    endif

    !   a, b, c, and d are combinations of neighboring potential
    ! vorticities which form the Arakawa and Hsu vorticity advection
    ! scheme.  All are defined at u grid points.

    if (CS%Coriolis_Scheme == ARAKAWA_HSU90) then
      do j=Jsq,Jeq+1
        do I=is-1,Ieq
          a(I,j) = (q(I,J) + (q(I+1,J) + q(I,J-1))) * C1_12
          d(I,j) = ((q(I,J) + q(I+1,J-1)) + q(I,J-1)) * C1_12
        enddo
        do I=Isq,Ieq
          b(I,j) = (q(I,J) + (q(I-1,J) + q(I,J-1))) * C1_12
          c(I,j) = ((q(I,J) + q(I-1,J-1)) + q(I,J-1)) * C1_12
        enddo
      enddo
    elseif (CS%Coriolis_Scheme == ARAKAWA_LAMB81) then
      do j=Jsq,Jeq+1 ; do I=Isq,Ieq+1
        a(I-1,j) = (2.0*(q(I,J) + q(I-1,J-1)) + (q(I-1,J) + q(I,J-1))) * C1_24
        d(I-1,j) = ((q(I,j) + q(I-1,J-1)) + 2.0*(q(I-1,J) + q(I,J-1))) * C1_24
        b(I,j) =   ((q(I,J) + q(I-1,J-1)) + 2.0*(q(I-1,J) + q(I,J-1))) * C1_24
        c(I,j) =   (2.0*(q(I,J) + q(I-1,J-1)) + (q(I-1,J) + q(I,J-1))) * C1_24
        ep_u(i,j) = ((q(I,J) - q(I-1,J-1)) + (q(I-1,J) - q(I,J-1))) * C1_24
        ep_v(i,j) = (-(q(I,J) - q(I-1,J-1)) + (q(I-1,J) - q(I,J-1))) * C1_24
      enddo ; enddo
    elseif (CS%Coriolis_Scheme == AL_BLEND) then
      Fe_m2 = CS%F_eff_max_blend - 2.0
      rat_lin = 1.5 * Fe_m2 / max(CS%wt_lin_blend, 1.0e-16)

      ! This allows the code to always give Sadourny Energy
      if (CS%F_eff_max_blend <= 2.0) then ; Fe_m2 = -1. ; rat_lin = -1.0 ; endif

      do j=Jsq,Jeq+1 ; do I=Isq,Ieq+1
        min_Ihq = MIN(Ih_q(I-1,J-1), Ih_q(I,J-1), Ih_q(I-1,J), Ih_q(I,J))
        max_Ihq = MAX(Ih_q(I-1,J-1), Ih_q(I,J-1), Ih_q(I-1,J), Ih_q(I,J))
        rat_m1 = 1.0e15
        if (max_Ihq < 1.0e15*min_Ihq) rat_m1 = max_Ihq / min_Ihq - 1.0
        ! The weights used here are designed to keep the effective Coriolis
        ! acceleration from any one point on its neighbors within a factor
        ! of F_eff_max.  The minimum permitted value is 2 (the factor for
        ! Sadourny's energy conserving scheme).

        ! Determine the relative weights of Arakawa & Lamb vs. Arakawa and Hsu.
        if (rat_m1 <= Fe_m2) then ; AL_wt = 1.0
        elseif (rat_m1 < 1.5*Fe_m2) then ; AL_wt = 3.0*Fe_m2 / rat_m1 - 2.0
        else ; AL_wt = 0.0 ; endif

        ! Determine the relative weights of Sadourny Energy vs. the other two.
        if (rat_m1 <= 1.5*Fe_m2) then ; Sad_wt = 0.0
        elseif (rat_m1 <= rat_lin) then
          Sad_wt = 1.0 - (1.5*Fe_m2) / rat_m1
        elseif (rat_m1 < 2.0*rat_lin) then
          Sad_wt = 1.0 - (CS%wt_lin_blend / rat_lin) * (rat_m1 - 2.0*rat_lin)
        else ; Sad_wt = 1.0 ; endif

        a(I-1,j) = Sad_wt * 0.25 * q(I-1,J) + (1.0 - Sad_wt) * &
                   ( ((2.0-AL_wt)* q(I-1,J) + AL_wt*q(I,J-1)) + &
                      2.0 * (q(I,J) + q(I-1,J-1)) ) * C1_24
        d(I-1,j) = Sad_wt * 0.25 * q(I-1,J-1) + (1.0 - Sad_wt) * &
                   ( ((2.0-AL_wt)* q(I-1,J-1) + AL_wt*q(I,J)) + &
                      2.0 * (q(I-1,J) + q(I,J-1)) ) * C1_24
        b(I,j) =   Sad_wt * 0.25 * q(I,J) + (1.0 - Sad_wt) * &
                   ( ((2.0-AL_wt)* q(I,J) + AL_wt*q(I-1,J-1)) + &
                      2.0 * (q(I-1,J) + q(I,J-1)) ) * C1_24
        c(I,j) =   Sad_wt * 0.25 * q(I,J-1) + (1.0 - Sad_wt) * &
                   ( ((2.0-AL_wt)* q(I,J-1) + AL_wt*q(I-1,J)) + &
                      2.0 * (q(I,J) + q(I-1,J-1)) ) * C1_24
        ep_u(i,j) = AL_wt  * ((q(I,J) - q(I-1,J-1)) + (q(I-1,J) - q(I,J-1))) * C1_24
        ep_v(i,j) = AL_wt * (-(q(I,J) - q(I-1,J-1)) + (q(I-1,J) - q(I,J-1))) * C1_24
      enddo ; enddo
    endif

    if (CS%Coriolis_En_Dis) then
    !  c1 = 1.0-1.5*RANGE ; c2 = 1.0-RANGE ; c3 = 2.0 ; slope = 0.5
      c1 = 1.0-1.5*0.5 ; c2 = 1.0-0.5 ; c3 = 2.0 ; slope = 0.5

      do j=Jsq,Jeq+1 ; do I=is-1,ie
        uhc = uh_center(I,j)
        uhm = uh(I,j,k)
        ! This sometimes matters with some types of open boundary conditions.
        if (G%dy_Cu(I,j) == 0.0) uhc = uhm

        if (abs(uhc) < 0.1*abs(uhm)) then
          uhm = 10.0*uhc
        elseif (abs(uhc) > c1*abs(uhm)) then
          if (abs(uhc) < c2*abs(uhm)) then ; uhc = (3.0*uhc+(1.0-c2*3.0)*uhm)
          elseif (abs(uhc) <= c3*abs(uhm)) then ; uhc = uhm
          else ; uhc = slope*uhc+(1.0-c3*slope)*uhm
          endif
        endif

        if (uhc > uhm) then
          uh_min(I,j) = uhm ; uh_max(I,j) = uhc
        else
          uh_max(I,j) = uhm ; uh_min(I,j) = uhc
        endif
      enddo ; enddo
      do J=js-1,je ; do i=Isq,Ieq+1
        vhc = vh_center(i,J)
        vhm = vh(i,J,k)
        ! This sometimes matters with some types of open boundary conditions.
        if (G%dx_Cv(i,J) == 0.0) vhc = vhm

        if (abs(vhc) < 0.1*abs(vhm)) then
          vhm = 10.0*vhc
        elseif (abs(vhc) > c1*abs(vhm)) then
          if (abs(vhc) < c2*abs(vhm)) then ; vhc = (3.0*vhc+(1.0-c2*3.0)*vhm)
          elseif (abs(vhc) <= c3*abs(vhm)) then ; vhc = vhm
          else ; vhc = slope*vhc+(1.0-c3*slope)*vhm
          endif
        endif

        if (vhc > vhm) then
          vh_min(i,J) = vhm ; vh_max(i,J) = vhc
        else
          vh_max(i,J) = vhm ; vh_min(i,J) = vhc
        endif
      enddo ; enddo
    endif

    ! Calculate KE and the gradient of KE
    call gradKE(u, v, h, KE, KEx, KEy, k, OBC, G, GV, US, CS)

    ! Calculate the tendencies of zonal velocity due to the Coriolis
    ! force and momentum advection.  On a Cartesian grid, this is
    !     CAu =  q * vh - d(KE)/dx.
    if (CS%Coriolis_Scheme == wenovi_7th_ENSTRO) then
      do j=js,je ; do I=Isq,Ieq    
        v_u = 0.25 * ((v(i+1,J,k) + v(i,J,k)) + (v(i,J-1,k) + v(i+1,J-1,k)))

        u_q1 = (uh(I,j-4,k) + uh(I,j-3,k)) * 0.5 
        u_q2 = (uh(I,j-3,k) + uh(I,j-2,k)) * 0.5
        u_q3 = (uh(I,j-2,k) + uh(I,j-1,k)) * 0.5
        u_q4 = (uh(I,j-1,k) + uh(I,j  ,k)) * 0.5
        u_q5 = (uh(I,j  ,k) + uh(I,j+1,k)) * 0.5
        u_q6 = (uh(I,j+1,k) + uh(I,j+2,k)) * 0.5
        u_q7 = (uh(I,j+2,k) + uh(I,j+3,k)) * 0.5
        u_q8 = (uh(I,j+3,k) + uh(I,j+4,k)) * 0.5

        v_q1 = (vh(i+1,J-4,k) + vh(i,J-4,k)) * 0.5
        v_q2 = (vh(i+1,J-3,k) + vh(i,J-3,k)) * 0.5
        v_q3 = (vh(i+1,J-2,k) + vh(i,J-2,k)) * 0.5
        v_q4 = (vh(i+1,J-1,k) + vh(i,J-1,k)) * 0.5
        v_q5 = (vh(i+1,J  ,k) + vh(i,J  ,k)) * 0.5
        v_q6 = (vh(i+1,J+1,k) + vh(i,J+1,k)) * 0.5
        v_q7 = (vh(i+1,J+2,k) + vh(i,J+2,k)) * 0.5
        v_q8 = (vh(i+1,J+3,k) + vh(i,J+3,k)) * 0.5
      
        third_order = (G%mask2dCu(I,j-2) * G%mask2dCu(I,j-1) * G%mask2dCu(I,j) * & 
                       G%mask2dCu(I,j+1) * G%mask2dCu(I,j+2))

        fifth_order   = third_order * G%mask2dCu(I,j-3) * G%mask2dCu(I,j+3) 
        seventh_order = fifth_order * G%mask2dCu(I,j-4) * G%mask2dCu(I,j-4)

        ! compute the masking to make sure that inland values are not used
        if (seventh_order == 1) then
            ! all values are valid, we use seventh order reconstruction
            call weno_seven_reconstruction(abs_vort(I,J-4),abs_vort(I,J-3),abs_vort(I,J-2),abs_vort(I,J-1), & 
                                           abs_vort(I,J)  ,abs_vort(I,J+1),abs_vort(I,J+2),abs_vort(I,J+3), & 
                                           u_q1, u_q2, u_q3, u_q4, u_q5, u_q6, u_q7, u_q8, v_u, q_u, CS%weno_velocity_smooth)

            ! call weno_seven_reconstruction(abs_vort(I,J-4),abs_vort(I,J-3),abs_vort(I,J-2),abs_vort(I,J-1), & 
            !                                abs_vort(I,J)  ,abs_vort(I,J+1),abs_vort(I,J+2),abs_vort(I,J+3), & 
            !                                v_q1, v_q2, v_q3, v_q4, v_q5, v_q6, v_q7, v_q8, v_u, q_uv, CS%weno_velocity_smooth)

            ! q_u = 0.5 * (q_uu + q_uv)
        elseif (fifth_order == 1) then
            ! all values are valid, we use fifth order reconstruction
            call weno_five_reconstruction(abs_vort(I,J-3),abs_vort(I,J-2),abs_vort(I,J-1), &
                                          abs_vort(I,J),  abs_vort(I,J+1),abs_vort(I,J+2), & 
                                          u_q2, u_q3, u_q4, u_q5, u_q6, u_q7, v_u, q_u, CS%weno_velocity_smooth)

        elseif (third_order == 1) then
            ! only the middle values are valid, we use third order reconstruction
            call weno_three_reconstruction(abs_vort(I,J-2),abs_vort(I,J-1),abs_vort(I,J),abs_vort(I,J+1), & 
                                           u_q3, u_q4, u_q5, u_q6, v_u, q_u, CS%weno_velocity_smooth)
        else ! Upwind first order
            if (v_u>0.) then
                q_u = abs_vort(I,J-1)
            else
                q_u = abs_vort(I,J)
            endif
        endif

        CAu(I,j,k) = (q_u * v_u) 
      enddo ; enddo            
    elseif (CS%Coriolis_Scheme == wenovi_7th_split) then
      do j=js,je ; do I=Isq,Ieq    
        v_u = 0.25 * ((v(i+1,J,k) + v(i,J,k)) + (v(i,J-1,k) + v(i+1,J-1,k)))

        u_q1 = (uh(I,j-4,k) + uh(I,j-3,k)) * 0.5 
        u_q2 = (uh(I,j-3,k) + uh(I,j-2,k)) * 0.5
        u_q3 = (uh(I,j-2,k) + uh(I,j-1,k)) * 0.5
        u_q4 = (uh(I,j-1,k) + uh(I,j  ,k)) * 0.5
        u_q5 = (uh(I,j  ,k) + uh(I,j+1,k)) * 0.5
        u_q6 = (uh(I,j+1,k) + uh(I,j+2,k)) * 0.5
        u_q7 = (uh(I,j+2,k) + uh(I,j+3,k)) * 0.5
        u_q8 = (uh(I,j+3,k) + uh(I,j+4,k)) * 0.5

      !  v_q1 = (vh(i+1,J-4,k) + vh(i,J-4,k)) * 0.5
      !  v_q2 = (vh(i+1,J-3,k) + vh(i,J-3,k)) * 0.5
      !  v_q3 = (vh(i+1,J-2,k) + vh(i,J-2,k)) * 0.5
      !  v_q4 = (vh(i+1,J-1,k) + vh(i,J-1,k)) * 0.5
      !  v_q5 = (vh(i+1,J  ,k) + vh(i,J  ,k)) * 0.5
      !  v_q6 = (vh(i+1,J+1,k) + vh(i,J+1,k)) * 0.5
      !  v_q7 = (vh(i+1,J+2,k) + vh(i,J+2,k)) * 0.5
      !  v_q8 = (vh(i+1,J+3,k) + vh(i,J+3,k)) * 0.5
      
        third_order = (G%mask2dCu(I,j-2) * G%mask2dCu(I,j-1) * G%mask2dCu(I,j) * & 
                       G%mask2dCu(I,j+1) * G%mask2dCu(I,j+2))

        fifth_order   = third_order * G%mask2dCu(I,j-3) * G%mask2dCu(I,j+3) 
        seventh_order = fifth_order * G%mask2dCu(I,j-4) * G%mask2dCu(I,j-4)

        ! compute the masking to make sure that inland values are not used
        if (seventh_order == 1) then
            ! all values are valid, we use seventh order reconstruction
            call weno_seven_reconstruction(rel_vort(I,J-4),rel_vort(I,J-3),rel_vort(I,J-2),rel_vort(I,J-1), & 
                                           rel_vort(I,J)  ,rel_vort(I,J+1),rel_vort(I,J+2),rel_vort(I,J+3), & 
                                           u_q1, u_q2, u_q3, u_q4, u_q5, u_q6, u_q7, u_q8, v_u, q_u, CS%weno_velocity_smooth)

            ! call weno_seven_reconstruction(abs_vort(I,J-4),abs_vort(I,J-3),abs_vort(I,J-2),abs_vort(I,J-1), & 
            !                                abs_vort(I,J)  ,abs_vort(I,J+1),abs_vort(I,J+2),abs_vort(I,J+3), & 
            !                                v_q1, v_q2, v_q3, v_q4, v_q5, v_q6, v_q7, v_q8, v_u, q_uv, CS%weno_velocity_smooth)

            ! q_u = 0.5 * (q_uu + q_uv)
        elseif (fifth_order == 1) then
            ! all values are valid, we use fifth order reconstruction
            call weno_five_reconstruction(rel_vort(I,J-3),rel_vort(I,J-2),rel_vort(I,J-1), &
                                          rel_vort(I,J),  rel_vort(I,J+1),rel_vort(I,J+2), & 
                                          u_q2, u_q3, u_q4, u_q5, u_q6, u_q7, v_u, q_u, CS%weno_velocity_smooth)

        elseif (third_order == 1) then
            ! only the middle values are valid, we use third order reconstruction
            call weno_three_reconstruction(rel_vort(I,J-2),rel_vort(I,J-1),rel_vort(I,J),rel_vort(I,J+1), & 
                                           u_q3, u_q4, u_q5, u_q6, v_u, q_u, CS%weno_velocity_smooth)
        else ! Upwind first order
            if (v_u>0.) then
                q_u = rel_vort(I,J-1)
            else
                q_u = rel_vort(I,J)
            endif
        endif

!        fv = 0.25 * &
!          ((G%CoriolisBu(I,J) * (v(i+1,J,k) + v(i,J,k))) + &
!           (G%CoriolisBu(I,J-1) * (v(i,J-1,k) + v(i+1,J-1,k))))  
!        fv = 0.25 * G%IdxCu(I,j) * &
!          ((G%CoriolisBu(I,J) * Ih_q(I,J) * (vh(i+1,J,k) + vh(i,J,k))) + &
!           (G%CoriolisBu(I,J-1) * Ih_q(I,J-1) * (vh(i,J-1,k) + vh(i+1,J-1,k))))  
        fq1 = (G%CoriolisBu(I,J)*Ih_q(I,J) + (G%CoriolisBu(I+1,J)*Ih_q(I+1,J) + & 
               G%CoriolisBu(I,J-1)*Ih_q(I,J-1))) * C1_12
        fq2 = (G%CoriolisBu(I,J)*Ih_q(I,J) + (G%CoriolisBu(I-1,J)*Ih_q(I-1,J) + &
               G%CoriolisBu(I,J-1)*Ih_q(I,J-1))) * C1_12
        fq3 = ((G%CoriolisBu(I,J)*Ih_q(I,J) + G%CoriolisBu(I-1,J-1)*Ih_q(I-1,J-1)) + &
                G%CoriolisBu(I,J-1)*Ih_q(I,J-1)) * C1_12
        fq4 = ((G%CoriolisBu(I,J)*Ih_q(I,J) + G%CoriolisBu(I+1,J-1)*Ih_q(I+1,J-1)) + &
                G%CoriolisBu(I,J-1)*Ih_q(I,J-1)) * C1_12
        fv = G%IdxCu(I,j) * &
          (fq1*vh(i+1,J,k) + fq2*vh(i,J,k) + fq3*vh(i,J-1,k) + fq4*vh(i+1,J-1,k))

        CAu(I,j,k) = (q_u * v_u) + fv
      enddo ; enddo            
    elseif (CS%Coriolis_Scheme == wenovi_5th_ENSTRO) then
      do j=js,je ; do I=Isq,Ieq
        v_u = 0.25 * ((v(i+1,J,k) + v(i,J,k)) + (v(i,J-1,k) + v(i+1,J-1,k)))

        v_q2 = (v(i+1,J-3,k) + v(i,J-3,k)) * 0.5
        v_q3 = (v(i+1,J-2,k) + v(i,J-2,k)) * 0.5
        v_q4 = (v(i+1,J-1,k) + v(i,J-1,k)) * 0.5
        v_q5 = (v(i+1,J  ,k) + v(i,J  ,k)) * 0.5
        v_q6 = (v(i+1,J+1,k) + v(i,J+1,k)) * 0.5
        v_q7 = (v(i+1,J+2,k) + v(i,J+2,k)) * 0.5
!
!        third_order = (G%mask2dCv(i+1,J-2) * G%mask2dCv(i+1,J-1) * G%mask2dCv(i+1,J) * G%mask2dCv(i+1,J+1) * &
!                       G%mask2dCv(i,J-2)   * G%mask2dCv(i,J-1)   * G%mask2dCv(i,J)   * G%mask2dCv(i,J+1))
!
!        fifth_order   = third_order * G%mask2dCv(i+1,J-3) * G%mask2dCv(i,J-3) * G%mask2dCv(i+1,J+2) * G%mask2dCv(i,J+2)
!
!        ! compute the masking to make sure that inland values are not used
!        if (fifth_order == 1) then
!            ! all values are valid, we use fifth order reconstruction
        call weno_five_reconstruction(rel_vort(I,J-3),rel_vort(I,J-2),rel_vort(I,J-1), &
                                      rel_vort(I,J),rel_vort(I,J+1),rel_vort(I,J+2),  &
                                      v_q2, v_q3, v_q4, v_q5, v_q6, v_q7, v_u, q_u, CS%weno_velocity_smooth)

!        elseif (third_order == 1) then
!            ! only the middle values are valid, we use third order reconstruction
!            call weno_three_reconstruction(abs_vort(I,J-2),abs_vort(I,J-1),abs_vort(I,J),abs_vort(I,J+1), &
!                                           v_q3, v_q4, v_q5, v_q6, v_u, q_u, CS%weno_velocity_smooth)
!        else ! Upwind first order
!            if (v_u>0.) then
!                q_u = abs_vort(I,J-1)
!            else
!                q_u = abs_vort(I,J)
!            endif
!        endif

        f_u = 0.5 * (G%CoriolisBu(I,J) + G%CoriolisBu(I,J-1))
        CAu(I,j,k) = ((q_u + f_u) * v_u)
      enddo ; enddo

    elseif (CS%Coriolis_Scheme == UP3_ENSTRO) then
      do j=js,je ; do I=Isq,Ieq
        v_u = 0.25 * ((v(i+1,J,k) + v(i,J,k)) + (v(i,J-1,k) + v(i+1,J-1,k)))
        call UP3_limiter_reconstruction(abs_vort(I,J-2), abs_vort(I,J-1),&
                abs_vort(I,J), abs_vort(I,J+1), v_u, abs_vort_u)
        CAu(I,j,k) = (abs_vort_u) * v_u 
      enddo ; enddo
    elseif (CS%Coriolis_Scheme == UP3_PV_ENSTRO) then
      do j=js,je ; do I=Isq,Ieq
        v_u = 0.25*G%IdxCu(I,j)*((vh(i+1,J,k) + vh(i,J,k)) + (vh(i,J-1,k) + vh(i+1,J-1,k)))
        call UP3_limiter_reconstruction(q(I,J-2), q(I,J-1),&
                q(I,J), q(I,J+1), v_u, q_u)
        CAu(I,j,k) = (q_u) * v_u 
      enddo ; enddo
    elseif (CS%Coriolis_Scheme == UP3_split) then
      do j=js,je ; do I=Isq,Ieq
        v_u = 0.25 * ((v(i+1,J,k) + v(i,J,k)) + (v(i,J-1,k) + v(i+1,J-1,k)))
        call UP3_limiter_reconstruction(rel_vort(I,J-2), rel_vort(I,J-1),&
                rel_vort(I,J), rel_vort(I,J+1), v_u, rel_vort_u)

!        fv = 0.25 * &
!          ((G%CoriolisBu(I,J) * (v(i+1,J,k) + v(i,J,k))) + &
!           (G%CoriolisBu(I,J-1) * (v(i,J-1,k) + v(i+1,J-1,k))))       
!        fv = 0.25 * G%IdxCu(I,j) * &
!          ((G%CoriolisBu(I,J) * Ih_q(I,J) * (vh(i+1,J,k) + vh(i,J,k))) + &
!           (G%CoriolisBu(I,J-1) * Ih_q(I,J-1) * (vh(i,J-1,k) + vh(i+1,J-1,k))))  
        fq1 = (G%CoriolisBu(I,J)*Ih_q(I,J) + (G%CoriolisBu(I+1,J)*Ih_q(I+1,J) + & 
               G%CoriolisBu(I,J-1)*Ih_q(I,J-1))) * C1_12
        fq2 = (G%CoriolisBu(I,J)*Ih_q(I,J) + (G%CoriolisBu(I-1,J)*Ih_q(I-1,J) + &
               G%CoriolisBu(I,J-1)*Ih_q(I,J-1))) * C1_12
        fq3 = ((G%CoriolisBu(I,J)*Ih_q(I,J) + G%CoriolisBu(I-1,J-1)*Ih_q(I-1,J-1)) + &
                G%CoriolisBu(I,J-1)*Ih_q(I,J-1)) * C1_12
        fq4 = ((G%CoriolisBu(I,J)*Ih_q(I,J) + G%CoriolisBu(I+1,J-1)*Ih_q(I+1,J-1)) + &
                G%CoriolisBu(I,J-1)*Ih_q(I,J-1)) * C1_12
        fv = G%IdxCu(I,j) * &
          (fq1*vh(i+1,J,k) + fq2*vh(i,J,k) + fq3*vh(i,J-1,k) + fq4*vh(i+1,J-1,k))
        CAu(I,j,k) = (rel_vort_u) * v_u + fv
      enddo ; enddo
    elseif (CS%Coriolis_Scheme == CEN4_ENSTRO) then
      do j=js,je ; do I=Isq,Ieq
        v_u = 0.25 * ((v(i+1,J,k) + v(i,J,k)) + (v(i,J-1,k) + v(i+1,J-1,k)))
        call CEN4_reconstruction(abs_vort(I,J-2), abs_vort(I,J-1),&
                abs_vort(I,J), abs_vort(I,J+1), abs_vort_u)
        CAu(I,j,k) = (abs_vort_u) * v_u
      enddo ; enddo
    endif

    if (CS%bound_Coriolis) then
      do j=js,je ; do I=Isq,Ieq
        fv1 = abs_vort(I,J) * v(i+1,J,k)
        fv2 = abs_vort(I,J) * v(i,J,k)
        fv3 = abs_vort(I,J-1) * v(i+1,J-1,k)
        fv4 = abs_vort(I,J-1) * v(i,J-1,k)

        max_fv = max(fv1, fv2, fv3, fv4)
        min_fv = min(fv1, fv2, fv3, fv4)

        CAu(I,j,k) = min(CAu(I,j,k), max_fv)
        CAu(I,j,k) = max(CAu(I,j,k), min_fv)
      enddo ; enddo
    endif

    ! Term - d(KE)/dx.
    do j=js,je ; do I=Isq,Ieq
      CAu(I,j,k) = CAu(I,j,k) - KEx(I,j)
    enddo ; enddo

    if (associated(AD%gradKEu)) then
      do j=js,je ; do I=Isq,Ieq
        AD%gradKEu(I,j,k) = -KEx(I,j)
      enddo ; enddo
    endif

    ! Calculate the tendencies of meridional velocity due to the Coriolis
    ! force and momentum advection.  On a Cartesian grid, this is
    !     CAv = - q * uh - d(KE)/dy.
    if (CS%Coriolis_Scheme == wenovi_7th_ENSTRO) then
      do J=Jsq,Jeq ; do i=is,ie
        u_v = 0.25* ((u(I-1,j,k) + u(I-1,j+1,k)) + (u(I,j,k) + u(I,j+1,k)))

        u_q1 = (uh(I-4,j+1,k) + uh(I-4,j,k)) * 0.5
        u_q2 = (uh(I-3,j+1,k) + uh(I-3,j,k)) * 0.5
        u_q3 = (uh(I-2,j+1,k) + uh(I-2,j,k)) * 0.5
        u_q4 = (uh(I-1,j+1,k) + uh(I-1,j,k)) * 0.5
        u_q5 = (uh(I  ,j+1,k) + uh(I  ,j,k)) * 0.5
        u_q6 = (uh(I+1,j+1,k) + uh(I+1,j,k)) * 0.5
        u_q7 = (uh(I+2,j+1,k) + uh(I+2,j,k)) * 0.5
        u_q8 = (uh(I+3,j+1,k) + uh(I+3,j,k)) * 0.5
        
        v_q1 = (vh(i-4,J,k) + vh(i-3,J,k)) * 0.5
        v_q2 = (vh(i-3,J,k) + vh(i-2,J,k)) * 0.5
        v_q3 = (vh(i-2,J,k) + vh(i-1,J,k)) * 0.5
        v_q4 = (vh(i-1,J,k) + vh(i  ,J,k)) * 0.5
        v_q5 = (vh(i  ,J,k) + vh(i+1,J,k)) * 0.5
        v_q6 = (vh(i+1,J,k) + vh(i+2,J,k)) * 0.5
        v_q7 = (vh(i+2,J,k) + vh(i+3,J,k)) * 0.5
        v_q8 = (vh(i+3,J,k) + vh(i+4,J,k)) * 0.5

        third_order = (G%mask2dCv(i-2,J) * G%mask2dCv(i-1,J) * G%mask2dCv(i,J) * G%mask2dCv(i+1,J) * & 
                       G%mask2dCv(i+2,J))

        fifth_order   = third_order * G%mask2dCv(i-3,J) * G%mask2dCv(i+3,J) 
        seventh_order = fifth_order * G%mask2dCv(i-4,J) * G%mask2dCv(i+4,J) 

        ! compute the masking to make sure that inland values are not used
        if (seventh_order == 1) then
            ! all values are valid, we use seventh order reconstruction
            call weno_seven_reconstruction(abs_vort(I-4,J),abs_vort(I-3,J),abs_vort(I-2,J),abs_vort(I-1,J), & 
                                           abs_vort(I,J)  ,abs_vort(I+1,J),abs_vort(I+2,J),abs_vort(I+3,J), & 
                                           v_q1, v_q2, v_q3, v_q4, v_q5, v_q6, v_q7, v_q8, u_v, q_v, CS%weno_velocity_smooth)

            ! all values are valid, we use seventh order reconstruction
            ! call weno_seven_reconstruction(abs_vort(I-4,J),abs_vort(I-3,J),abs_vort(I-2,J),abs_vort(I-1,J), & 
            !                                abs_vort(I,J)  ,abs_vort(I+1,J),abs_vort(I+2,J),abs_vort(I+3,J), & 
            !                                u_q1, u_q2, u_q3, u_q4, u_q5, u_q6, u_q7, u_q8, u_v, q_v, CS%weno_velocity_smooth)

        elseif (fifth_order == 1) then
            ! all values are valid, we use fifth order reconstruction
            call weno_five_reconstruction(abs_vort(I-3,J),abs_vort(I-2,J),abs_vort(I-1,J), &
                                          abs_vort(I,J),abs_vort(I+1,J),abs_vort(I+2,J), & 
                                          v_q2, v_q3, v_q4, v_q5, v_q6, v_q7, u_v, q_v, CS%weno_velocity_smooth)

        elseif (third_order == 1) then
            ! only the middle values are valid, we use third order reconstruction
                call weno_three_reconstruction(abs_vort(I-2,J),abs_vort(I-1,J),abs_vort(I,J),abs_vort(I+1,J), & 
                                               v_q3, v_q4, v_q5, v_q6, u_v, q_v, CS%weno_velocity_smooth)
        else ! Upwind first order!
            if (u_v>0.) then
                q_v = abs_vort(I-1,J)
            else
                q_v = abs_vort(I,J)
            endif
        endif

        CAv(i,J,k) = - (q_v * u_v)
      enddo ; enddo
    elseif (CS%Coriolis_Scheme == wenovi_7th_split) then
      do J=Jsq,Jeq ; do i=is,ie
        u_v = 0.25* ((u(I-1,j,k) + u(I-1,j+1,k)) + (u(I,j,k) + u(I,j+1,k)))

        u_q1 = (uh(I-4,j+1,k) + uh(I-4,j,k)) * 0.5
        u_q2 = (uh(I-3,j+1,k) + uh(I-3,j,k)) * 0.5
        u_q3 = (uh(I-2,j+1,k) + uh(I-2,j,k)) * 0.5
        u_q4 = (uh(I-1,j+1,k) + uh(I-1,j,k)) * 0.5
        u_q5 = (uh(I  ,j+1,k) + uh(I  ,j,k)) * 0.5
        u_q6 = (uh(I+1,j+1,k) + uh(I+1,j,k)) * 0.5
        u_q7 = (uh(I+2,j+1,k) + uh(I+2,j,k)) * 0.5
        u_q8 = (uh(I+3,j+1,k) + uh(I+3,j,k)) * 0.5
        
        v_q1 = (vh(i-4,J,k) + vh(i-3,J,k)) * 0.5
        v_q2 = (vh(i-3,J,k) + vh(i-2,J,k)) * 0.5
        v_q3 = (vh(i-2,J,k) + vh(i-1,J,k)) * 0.5
        v_q4 = (vh(i-1,J,k) + vh(i  ,J,k)) * 0.5
        v_q5 = (vh(i  ,J,k) + vh(i+1,J,k)) * 0.5
        v_q6 = (vh(i+1,J,k) + vh(i+2,J,k)) * 0.5
        v_q7 = (vh(i+2,J,k) + vh(i+3,J,k)) * 0.5
        v_q8 = (vh(i+3,J,k) + vh(i+4,J,k)) * 0.5

        third_order = (G%mask2dCv(i-2,J) * G%mask2dCv(i-1,J) * G%mask2dCv(i,J) * G%mask2dCv(i+1,J) * & 
                       G%mask2dCv(i+2,J))

        fifth_order   = third_order * G%mask2dCv(i-3,J) * G%mask2dCv(i+3,J) 
        seventh_order = fifth_order * G%mask2dCv(i-4,J) * G%mask2dCv(i+4,J) 

        ! compute the masking to make sure that inland values are not used
        if (seventh_order == 1) then
            ! all values are valid, we use seventh order reconstruction
            call weno_seven_reconstruction(rel_vort(I-4,J),rel_vort(I-3,J),rel_vort(I-2,J),rel_vort(I-1,J), & 
                                           rel_vort(I,J)  ,rel_vort(I+1,J),rel_vort(I+2,J),rel_vort(I+3,J), & 
                                           v_q1, v_q2, v_q3, v_q4, v_q5, v_q6, v_q7, v_q8, u_v, q_v, CS%weno_velocity_smooth)

            ! all values are valid, we use seventh order reconstruction
            ! call weno_seven_reconstruction(abs_vort(I-4,J),abs_vort(I-3,J),abs_vort(I-2,J),abs_vort(I-1,J), & 
            !                                abs_vort(I,J)  ,abs_vort(I+1,J),abs_vort(I+2,J),abs_vort(I+3,J), & 
            !                                u_q1, u_q2, u_q3, u_q4, u_q5, u_q6, u_q7, u_q8, u_v, q_v, CS%weno_velocity_smooth)

        elseif (fifth_order == 1) then
            ! all values are valid, we use fifth order reconstruction
            call weno_five_reconstruction(rel_vort(I-3,J),rel_vort(I-2,J),rel_vort(I-1,J), &
                                          rel_vort(I,J),rel_vort(I+1,J),rel_vort(I+2,J), & 
                                          v_q2, v_q3, v_q4, v_q5, v_q6, v_q7, u_v, q_v, CS%weno_velocity_smooth)

        elseif (third_order == 1) then
            ! only the middle values are valid, we use third order reconstruction
                call weno_three_reconstruction(rel_vort(I-2,J),rel_vort(I-1,J),rel_vort(I,J),rel_vort(I+1,J), & 
                                               v_q3, v_q4, v_q5, v_q6, u_v, q_v, CS%weno_velocity_smooth)
        else ! Upwind first order!
            if (u_v>0.) then
                q_v = rel_vort(I-1,J)
            else
                q_v = rel_vort(I,J)
            endif
        endif

!        fu = - 0.25* &
!            ((G%CoriolisBu(I-1,J)*(u(I-1,j,k) + u(I-1,j+1,k))) + &
!             (G%CoriolisBu(I,J)*(u(I,j,k) + u(I,j+1,k)))) 
!        fu = - 0.25 * G%IdyCv(i,J) * &
!            ((G%CoriolisBu(I-1,J)*Ih_q(I-1,J)*(uh(I-1,j,k) + uh(I-1,j+1,k))) + &
!             (G%CoriolisBu(I,J)*Ih_q(I,J)*(uh(I,j,k) + uh(I,j+1,k))))         
        fq1 = (G%CoriolisBu(I,J)*Ih_q(I,J) + (G%CoriolisBu(I-1,J)*Ih_q(I-1,J) + & 
               G%CoriolisBu(I-1,J-1)*Ih_q(I-1,J-1))) * C1_12
        fq2 = (G%CoriolisBu(I,J)*Ih_q(I,J) + (G%CoriolisBu(I-1,J)*Ih_q(I-1,J) + &
               G%CoriolisBu(I,J-1)*Ih_q(I,J-1))) * C1_12
        fq3 = ((G%CoriolisBu(I,J)*Ih_q(I,J) + G%CoriolisBu(I-1,J)*Ih_q(I-1,J)) + &
                G%CoriolisBu(I,J+1)*Ih_q(I,J+1)) * C1_12
        fq4 = ((G%CoriolisBu(I,J)*Ih_q(I,J) + G%CoriolisBu(I-1,J)*Ih_q(I-1,J)) + &
                G%CoriolisBu(I-1,J+1)*Ih_q(I-1,J+1)) * C1_12
        fu = - G%IdyCv(i,J) * &
          (fq1*uh(I-1,j,k) + fq2*uh(I,j,k) + fq3*uh(I,j+1,k) + fq4*uh(I-1,j+1,k))

        CAv(i,J,k) = - (q_v * u_v) + fu
      enddo ; enddo
    elseif (CS%Coriolis_Scheme == wenovi_5th_ENSTRO) then
      do J=Jsq,Jeq ; do i=is,ie
        u_v = 0.25* ((u(I-1,j,k) + u(I-1,j+1,k)) + (u(I,j,k) + u(I,j+1,k)))

        u_q2 = (u(I-3,j+1,k) + u(I-3,j,k)) * 0.5
        u_q3 = (u(I-2,j+1,k) + u(I-2,j,k)) * 0.5
        u_q4 = (u(I-1,j+1,k) + u(I-1,j,k)) * 0.5
        u_q5 = (u(I  ,j+1,k) + u(I  ,j,k)) * 0.5
        u_q6 = (u(I+1,j+1,k) + u(I+1,j,k)) * 0.5
        u_q7 = (u(I+2,j+1,k) + u(I+2,j,k)) * 0.5
!
!        third_order = (G%mask2dCu(I-2,j+1) * G%mask2dCu(I-1,j+1) * G%mask2dCu(I,j+1) * G%mask2dCu(I+1,j+1) * &
!                       G%mask2dCu(I-2,j)   * G%mask2dCu(I-1,j)   * G%mask2dCu(I,j)   * G%mask2dCu(I+1,j))
!
!        fifth_order   = third_order * G%mask2dCu(I-3,j+1) * G%mask2dCu(I-3,j) * G%mask2dCu(I+2,j+1) * G%mask2dCu(I+2,j)

!        ! compute the masking to make sure that inland values are not used
!
!        if (fifth_order == 1) then
!            ! all values are valid, we use fifth order reconstruction
        call weno_five_reconstruction(rel_vort(I-3,J),rel_vort(I-2,J),rel_vort(I-1,J), &
                                      rel_vort(I,J),rel_vort(I+1,J),rel_vort(I+2,J), &
                                      u_q2, u_q3, u_q4, u_q5, u_q6, u_q7, u_v, q_v, CS%weno_velocity_smooth)
!        elseif (third_order == 1) then
!            ! only the middle values are valid, we use third order reconstruction
!                call weno_three_reconstruction(abs_vort(I-2,J),abs_vort(I-1,J),abs_vort(I,J),abs_vort(I+1,J), &
!                                               u_q3, u_q4, u_q5, u_q6, u_v, q_v, CS%weno_velocity_smooth)
!        else ! Upwind first order!
!        if (u_v>0.) then
!          q_v = abs_vort(I-1,J)
!        else
!          q_v = abs_vort(I,J)
!        endif
!        endif

        f_v = 0.5 * (G%CoriolisBu(I-1,J) + G%CoriolisBu(I,J))
        CAv(i,J,k) = - ((q_v + f_v) * u_v)
      enddo ; enddo      
    elseif (CS%Coriolis_Scheme == UP3_ENSTRO) then
      do J=Jsq,Jeq ; do i=is,ie
        u_v = 0.25* ((u(I-1,j,k) + u(I-1,j+1,k)) + (u(I,j,k) + u(I,j+1,k)))
        call UP3_limiter_reconstruction(abs_vort(I-2,J), abs_vort(I-1,J),&
                abs_vort(I,J), abs_vort(I+1,J), u_v, abs_vort_v) 
        CAv(i,J,k) = - (abs_vort_v) * u_v
      enddo ; enddo
    elseif (CS%Coriolis_Scheme == UP3_PV_ENSTRO) then
      do J=Jsq,Jeq ; do i=is,ie
        u_v = 0.25*G%IdyCv(i,J)*((uh(I-1,j,k) + uh(I-1,j+1,k)) + (uh(I,j,k) + uh(I,j+1,k)))
        call UP3_limiter_reconstruction(q(I-2,J), q(I-1,J),&
                q(I,J), q(I+1,J), u_v, q_v) 
        CAv(i,J,k) = - (q_v) * u_v
      enddo ; enddo
    elseif (CS%Coriolis_Scheme == UP3_split) then
      do J=Jsq,Jeq ; do i=is,ie
        u_v = 0.25* ((u(I-1,j,k) + u(I-1,j+1,k)) + (u(I,j,k) + u(I,j+1,k)))
        call UP3_limiter_reconstruction(rel_vort(I-2,J), rel_vort(I-1,J),&
                rel_vort(I,J), rel_vort(I+1,J), u_v, rel_vort_v) 
!        fu = - 0.25* &
!            ((G%CoriolisBu(I-1,J)*(u(I-1,j,k) + u(I-1,j+1,k))) + &
!             (G%CoriolisBu(I,J)*(u(I,j,k) + u(I,j+1,k))))         
!        fu = - 0.25 * G%IdyCv(i,J) *&
!            ((G%CoriolisBu(I-1,J)*Ih_q(I-1,J)*(uh(I-1,j,k) + uh(I-1,j+1,k))) + &
!             (G%CoriolisBu(I,J)*Ih_q(I,J)*(uh(I,j,k) + uh(I,j+1,k))))         
        fq1 = (G%CoriolisBu(I,J)*Ih_q(I,J) + (G%CoriolisBu(I-1,J)*Ih_q(I-1,J) + & 
               G%CoriolisBu(I-1,J-1)*Ih_q(I-1,J-1))) * C1_12
        fq2 = (G%CoriolisBu(I,J)*Ih_q(I,J) + (G%CoriolisBu(I-1,J)*Ih_q(I-1,J) + &
               G%CoriolisBu(I,J-1)*Ih_q(I,J-1))) * C1_12
        fq3 = ((G%CoriolisBu(I,J)*Ih_q(I,J) + G%CoriolisBu(I-1,J)*Ih_q(I-1,J)) + &
                G%CoriolisBu(I,J+1)*Ih_q(I,J+1)) * C1_12
        fq4 = ((G%CoriolisBu(I,J)*Ih_q(I,J) + G%CoriolisBu(I-1,J)*Ih_q(I-1,J)) + &
                G%CoriolisBu(I-1,J+1)*Ih_q(I-1,J+1)) * C1_12
        fu = - G%IdyCv(i,J) * &
          (fq1*uh(I-1,j,k) + fq2*uh(I,j,k) + fq3*uh(I,j+1,k) + fq4*uh(I-1,j+1,k))
        CAv(i,J,k) = - (rel_vort_v) * u_v + fu
      enddo ; enddo
    elseif (CS%Coriolis_Scheme == CEN4_ENSTRO) then
      do J=Jsq,Jeq ; do i=is,ie
        u_v = 0.25* ((u(I-1,j,k) + u(I-1,j+1,k)) + (u(I,j,k) + u(I,j+1,k)))
        call CEN4_reconstruction(abs_vort(I-2,J), abs_vort(I-1,J),&
                abs_vort(I,J), abs_vort(I+1,J), abs_vort_v) 
        CAv(i,J,k) = - (abs_vort_v) * u_v
      enddo ; enddo
    endif
    ! Add in the additonal terms with Arakawa & Lamb.
    if ((CS%Coriolis_Scheme == ARAKAWA_LAMB81) .or. &
        (CS%Coriolis_Scheme == AL_BLEND)) then ; do J=Jsq,Jeq ; do i=is,ie
      CAv(i,J,k) = CAv(i,J,k) + &
            ((ep_v(i,j)*vh(i,J-1,k)) - (ep_v(i,j+1)*vh(i,J+1,k))) * G%IdyCv(i,J)
    enddo ; enddo ; endif

    if (Stokes_VF) then
      if (CS%id_CAuS>0 .or. CS%id_CAvS>0) then
        ! Computing the diagnostic Stokes contribution to CAv
        do J=Jsq,Jeq ; do i=is,ie
          CAvS(I,j,k) = 0.25 * &
                ((qS(I,J) * (uh(I,j+1,k) + uh(I,j,k))) + &
                 (qS(I,J-1) * (uh(I-1,j,k) + uh(I-1,j+1,k)))) * G%IdyCv(i,J)
        enddo; enddo
      endif
    endif

    if (CS%bound_Coriolis) then
      do J=Jsq,Jeq ; do i=is,ie
        fu1 = -abs_vort(I,J) * u(I,j+1,k)
        fu2 = -abs_vort(I,J) * u(I,j,k)
        fu3 = -abs_vort(I-1,J) * u(I-1,j+1,k)
        fu4 = -abs_vort(I-1,J) * u(I-1,j,k)

        max_fu = max(fu1, fu2, fu3, fu4)
        min_fu = min(fu1, fu2, fu3, fu4)

        CAv(I,j,k) = min(CAv(I,j,k), max_fu)
        CAv(I,j,k) = max(CAv(I,j,k), min_fu)
      enddo ; enddo
    endif    

    ! Term - d(KE)/dy.
    do J=Jsq,Jeq ; do i=is,ie
      CAv(i,J,k) = CAv(i,J,k) - KEy(i,J)
    enddo ; enddo
    if (associated(AD%gradKEv)) then
      do J=Jsq,Jeq ; do i=is,ie
        AD%gradKEv(i,J,k) = -KEy(i,J)
      enddo ; enddo
    endif

    if (associated(AD%rv_x_u) .or. associated(AD%rv_x_v)) then
      ! Calculate the Coriolis-like acceleration due to relative vorticity.
      if (CS%Coriolis_Scheme == SADOURNY75_ENERGY) then
        if (associated(AD%rv_x_u)) then
          do J=Jsq,Jeq ; do i=is,ie
            AD%rv_x_u(i,J,k) = - 0.25* &
              ((q2(I-1,j)*(uh(I-1,j,k) + uh(I-1,j+1,k))) + &
               (q2(I,j)*(uh(I,j,k) + uh(I,j+1,k)))) * G%IdyCv(i,J)
          enddo ; enddo
        endif

        if (associated(AD%rv_x_v)) then
          do j=js,je ; do I=Isq,Ieq
            AD%rv_x_v(I,j,k) = 0.25 * &
              ((q2(I,j) * (vh(i+1,J,k) + vh(i,J,k))) + &
               (q2(I,j-1) * (vh(i,J-1,k) + vh(i+1,J-1,k)))) * G%IdxCu(I,j)
          enddo ; enddo
        endif
      else
        if (associated(AD%rv_x_u)) then
          do J=Jsq,Jeq ; do i=is,ie
            AD%rv_x_u(i,J,k) = -G%IdyCv(i,J) * C1_12 * &
              (((((q2(I,J) + q2(I-1,J-1)) + q2(I-1,J)) * uh(I-1,j,k)) + &
                (((q2(I-1,J) + q2(I,J+1)) + q2(I,J)) * uh(I,j+1,k))) + &
               ((((q2(I-1,J) + q2(I,J-1)) + q2(I,J)) * uh(I,j,k))+ &
                (((q2(I,J) + q2(I-1,J+1)) + q2(I-1,J)) * uh(I-1,j+1,k))))
          enddo ; enddo
        endif

        if (associated(AD%rv_x_v)) then
          do j=js,je ; do I=Isq,Ieq
            AD%rv_x_v(I,j,k) = G%IdxCu(I,j) * C1_12 * &
              (((((q2(I+1,J) + q2(I,J-1)) + q2(I,J)) * vh(i+1,J,k)) + &
                (((q2(I-1,J-1) + q2(I,J)) + q2(I,J-1)) * vh(i,J-1,k))) + &
               ((((q2(I-1,J) + q2(I,J-1)) + q2(I,J)) * vh(i,J,k)) + &
                (((q2(I+1,J-1) + q2(I,J)) + q2(I,J-1)) * vh(i+1,J-1,k))))
          enddo ; enddo
        endif
      endif
    endif

  enddo ! k-loop.

  ! Here the various Coriolis-related derived quantities are offered for averaging.
  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_rv > 0) call post_data(CS%id_rv, RV, CS%diag)
    if (CS%id_PV > 0) call post_data(CS%id_PV, PV, CS%diag)
    if (CS%id_gKEu>0) call post_data(CS%id_gKEu, AD%gradKEu, CS%diag)
    if (CS%id_gKEv>0) call post_data(CS%id_gKEv, AD%gradKEv, CS%diag)
    if (CS%id_rvxu > 0) call post_data(CS%id_rvxu, AD%rv_x_u, CS%diag)
    if (CS%id_rvxv > 0) call post_data(CS%id_rvxv, AD%rv_x_v, CS%diag)
    if (Stokes_VF) then
      if (CS%id_CAuS > 0) call post_data(CS%id_CAuS, CAuS, CS%diag)
      if (CS%id_CAvS > 0) call post_data(CS%id_CAvS, CAvS, CS%diag)
    endif

    ! Diagnostics for terms multiplied by fractional thicknesses

    ! 3D diagnostics hf_gKEu etc. are commented because there is no clarity on proper remapping grid option.
    ! The code is retained for debugging purposes in the future.
    ! if (CS%id_hf_gKEu > 0) call post_product_u(CS%id_hf_gKEu, AD%gradKEu, AD%diag_hfrac_u, G, nz, CS%diag)
    ! if (CS%id_hf_gKEv > 0) call post_product_v(CS%id_hf_gKEv, AD%gradKEv, AD%diag_hfrac_v, G, nz, CS%diag)
    ! if (CS%id_hf_rvxv > 0) call post_product_u(CS%id_hf_rvxv, AD%rv_x_v, AD%diag_hfrac_u, G, nz, CS%diag)
    ! if (CS%id_hf_rvxu > 0) call post_product_v(CS%id_hf_rvxu, AD%rv_x_u, AD%diag_hfrac_v, G, nz, CS%diag)

    if (CS%id_hf_gKEu_2d > 0) call post_product_sum_u(CS%id_hf_gKEu_2d, AD%gradKEu, AD%diag_hfrac_u, G, nz, CS%diag)
    if (CS%id_hf_gKEv_2d > 0) call post_product_sum_v(CS%id_hf_gKEv_2d, AD%gradKEv, AD%diag_hfrac_v, G, nz, CS%diag)
    if (CS%id_intz_gKEu_2d > 0) call post_product_sum_u(CS%id_intz_gKEu_2d, AD%gradKEu, AD%diag_hu, G, nz, CS%diag)
    if (CS%id_intz_gKEv_2d > 0) call post_product_sum_v(CS%id_intz_gKEv_2d, AD%gradKEv, AD%diag_hv, G, nz, CS%diag)

    if (CS%id_hf_rvxv_2d > 0) call post_product_sum_u(CS%id_hf_rvxv_2d, AD%rv_x_v, AD%diag_hfrac_u, G, nz, CS%diag)
    if (CS%id_hf_rvxu_2d > 0) call post_product_sum_v(CS%id_hf_rvxu_2d, AD%rv_x_u, AD%diag_hfrac_v, G, nz, CS%diag)

    if (CS%id_h_gKEu > 0) call post_product_u(CS%id_h_gKEu, AD%gradKEu, AD%diag_hu, G, nz, CS%diag)
    if (CS%id_h_gKEv > 0) call post_product_v(CS%id_h_gKEv, AD%gradKEv, AD%diag_hv, G, nz, CS%diag)
    if (CS%id_h_rvxv > 0) call post_product_u(CS%id_h_rvxv, AD%rv_x_v, AD%diag_hu, G, nz, CS%diag)
    if (CS%id_h_rvxu > 0) call post_product_v(CS%id_h_rvxu, AD%rv_x_u, AD%diag_hv, G, nz, CS%diag)

    if (CS%id_intz_rvxv_2d > 0) call post_product_sum_u(CS%id_intz_rvxv_2d, AD%rv_x_v, AD%diag_hu, G, nz, CS%diag)
    if (CS%id_intz_rvxu_2d > 0) call post_product_sum_v(CS%id_intz_rvxu_2d, AD%rv_x_u, AD%diag_hv, G, nz, CS%diag)
  endif

end subroutine CorAdCalc

subroutine UP3_limiter_reconstruction(q1,q2,q3,q4,u,qr)
  real, intent(in)    :: q1, q2, q3, q4
  real, intent(in)    :: u
  real, intent(inout) :: qr
  real                :: theta, psi

  if (u>0.) then
    theta = (q2 - q1)/(q3 - q2 + 1e-20)
    psi = max(0., min(1., 1/3. + 1/6.*theta, theta)) ! limiter introduced by Koren (1993)
    qr = q2 + psi*(q3 - q2)
  else
    theta = (q4 - q3)/(q3 - q2 + 1e-20)
    psi = max(0., min(1., 1/3. + 1/6.*theta, theta))
    qr = q3 + psi*(q2 - q3)
  endif  

end subroutine UP3_limiter_reconstruction

subroutine CEN4_reconstruction(q1,q2,q3,q4,qr)
  real, intent(in)    :: q1, q2, q3, q4
  real, intent(inout) :: qr

  qr = (-q1 + 7*q2 + 7*q3 - q4)/12.

end subroutine CEN4_reconstruction

subroutine weno_three_reconstruction(q1, q2, q3, q4, u1, u2, u3, u4, u, qr, velocity_smoothing)
    real, intent(in)    :: q1, q2, q3, q4
    real, intent(in)    :: u1, u2, u3, u4
    real, intent(in)    :: u
    logical, intent(in) :: velocity_smoothing
    real, intent(inout) :: qr
    real :: c0, c1
    real :: b0, b1
    real :: tau, w0, w1
    real :: s
    
    if (u>0.) then
      call weno_three_reconstruction_0(q2, q3, c0)
      call weno_three_reconstruction_1(q1, q2, c1) 
      if (velocity_smoothing == .true.) then
        call weno_three_weight(u2, u3, b0)
        call weno_three_weight(u1, u2, b1)
      else
        call weno_three_weight(q2, q3, b0)
        call weno_three_weight(q1, q2, b1)
      endif
    else
      call weno_three_reconstruction_0(q3, q2, c0)
      call weno_three_reconstruction_1(q4, q3, c1)
      if (velocity_smoothing == .true.) then
        call weno_three_weight(u3, u2, b0)
        call weno_three_weight(u4, u3, b1)
      else
        call weno_three_weight(q3, q2, b0)
        call weno_three_weight(q4, q3, b1)
      endif
    endif
  
    tau = abs(b0-b1)
    w0  = 2./3. * (1 + (tau / (b0 + 1e-20))**2)
    w1  = 1./3. * (1 + (tau / (b1 + 1e-20))**2)
  
    s = 1. / (w0 + w1)
    w0 = w0 * s
    w1 = w1 * s
  
    qr = w0 * c0 + w1 * c1
  
end subroutine weno_three_reconstruction

subroutine weno_three_weight(q0, q1, w0)
    real, intent(in) :: q0,q1
    real, intent(inout) :: w0
  
    w0 = q0 * q0 - 2 * q0 * q1 + q1 * q1
    
end subroutine weno_three_weight

subroutine weno_three_reconstruction_0(q0, q1, w0)
    real, intent(in) :: q0,q1
    real, intent(inout) :: w0
  
    w0 = (q0 + q1) * 0.5

end subroutine weno_three_reconstruction_0

subroutine weno_three_reconstruction_1(q0, q1, w0)
    real, intent(in) :: q0,q1
    real, intent(inout) :: w0
  
    w0 = (- q0 + 3 * q1) * 0.5

end subroutine weno_three_reconstruction_1

subroutine weno_five_reconstruction(q1, q2, q3, q4, q5, q6, u1, u2, u3, u4, u5, u6, u, qr, velocity_smoothing)
    real, intent(in)    :: q1, q2, q3, q4, q5, q6
    real, intent(in)    :: u1, u2, u3, u4, u5, u6
    real, intent(in)    :: u
    logical, intent(in) :: velocity_smoothing
    real, intent(inout) :: qr
    real :: c0, c1, c2
    real :: b0, b1, b2
    real :: tau, w0, w1, w2
    real :: s
  
    if (u>0.) then
      call weno_five_reconstruction_0(q3, q4, q5, c0)
      call weno_five_reconstruction_1(q2, q3, q4, c1) 
      call weno_five_reconstruction_2(q1, q2, q3, c2)
      if (velocity_smoothing == .true.) then
        call weno_five_weight_0(u3, u4, u5, b0)
        call weno_five_weight_1(u2, u3, u4, b1)
        call weno_five_weight_2(u1, u2, u3, b2)
      else
        call weno_five_weight_0(q3, q4, q5, b0)
        call weno_five_weight_1(q2, q3, q4, b1)
        call weno_five_weight_2(q1, q2, q3, b2)
      endif
    else
      call weno_five_reconstruction_0(q4, q3, q2, c0)
      call weno_five_reconstruction_1(q5, q4, q3, c1)
      call weno_five_reconstruction_2(q6, q5, q4, c2)
      if (velocity_smoothing == .true.) then
        call weno_five_weight_0(u4, u3, u2, b0)
        call weno_five_weight_1(u5, u4, u3, b1)
        call weno_five_weight_2(u6, u5, u4, b2)
      else
        call weno_five_weight_0(q4, q3, q2, b0)
        call weno_five_weight_1(q5, q4, q3, b1)
        call weno_five_weight_2(q6, q5, q4, b2)
      endif
    endif
  
    tau = abs(b0 - b2)
    w0  = 3./10. * (1 + (tau / (b0 + 1e-20))**2)
    w1  = 3./5.  * (1 + (tau / (b1 + 1e-20))**2)
    w2  = 1./10. * (1 + (tau / (b2 + 1e-20))**2)
  
    s = 1. / (w0 + w1 + w2)
    w0 = w0 * s
    w1 = w1 * s
    w2 = w2 * s
  
    qr = w0 * c0 + w1 * c1 + w2 * c2
  
end subroutine weno_five_reconstruction

subroutine weno_five_weight_0(q0, q1, q2, w0)
  real, intent(in) :: q0,q1,q2
  real, intent(inout) :: w0

  w0 = q0 * (10 * q0 - 31 * q1 + 11 * q2) + q1 * (25 * q1 - 19 * q2) + 4 * q2 * q2
  
end subroutine weno_five_weight_0

subroutine weno_five_weight_1(q0, q1, q2, w1)
  real, intent(in) :: q0,q1,q2
  real, intent(inout) :: w1

  w1 = q0 * (4 * q0 - 13 * q1 + 5 * q2) + q1 * (13 * q1 - 13 * q2) + 4 * q2 * q2

end subroutine weno_five_weight_1

subroutine weno_five_weight_2(q0, q1, q2, w2)
  real, intent(in) :: q0, q1, q2
  real, intent(inout) :: w2

  w2 = q0 * (4 * q0 - 19 * q1 + 11 * q2) + q1 * (25 * q1 - 31 * q2) + 10 * q2 * q2

end subroutine weno_five_weight_2

subroutine weno_five_reconstruction_0(q0, q1, q2, p0)
  real, intent(in) :: q0,q1,q2
  real, intent(inout) :: p0

  p0 = (2*q0 + 5*q1 - q2) / 6.

end subroutine weno_five_reconstruction_0

subroutine weno_five_reconstruction_1(q0, q1, q2, p1)
  real, intent(in) :: q0,q1,q2
  real, intent(inout) :: p1

  p1 = (-q0 + 5*q1 + 2*q2) / 6.

end subroutine weno_five_reconstruction_1

subroutine weno_five_reconstruction_2(q0, q1, q2, p2)
  real, intent(in) :: q0,q1,q2
  real, intent(inout) :: p2

  p2 = (2*q0 - 7*q1 + 11*q2) / 6.

end subroutine weno_five_reconstruction_2

subroutine weno_seven_reconstruction(q1, q2, q3, q4, q5, q6, q7, q8, u1, u2, u3, u4, u5, u6, u7, u8, u, qr, velocity_smoothing)
  real, intent(in)    :: q1, q2, q3, q4, q5, q6, q7, q8
  real, intent(in)    :: u1, u2, u3, u4, u5, u6, u7, u8
  real, intent(in)    :: u
  logical, intent(in) :: velocity_smoothing
  real, intent(inout) :: qr
  real :: c0, c1, c2, c3
  real :: b0, b1, b2, b3
  real :: tau, w0, w1, w2, w3
  real :: s

  if (u>0.) then
    call weno_seven_reconstruction_0(q4, q5, q6, q7, c0)
    call weno_seven_reconstruction_1(q3, q4, q5, q6, c1) 
    call weno_seven_reconstruction_2(q2, q3, q4, q5, c2)
    call weno_seven_reconstruction_3(q1, q2, q3, q4, c3)
    if (velocity_smoothing == .true.) then
      call weno_seven_weight_0(u4, u5, u6, u7, b0)
      call weno_seven_weight_1(u3, u4, u5, u6, b1)
      call weno_seven_weight_2(u2, u3, u4, u5, b2)
      call weno_seven_weight_3(u1, u2, u3, u4, b3)
    else
      call weno_seven_weight_0(q4, q5, q6, q7, b0)
      call weno_seven_weight_1(q3, q4, q5, q6, b1)
      call weno_seven_weight_2(q2, q3, q4, q5, b2)
      call weno_seven_weight_3(q1, q2, q3, q4, b3)
    endif
  else
    call weno_seven_reconstruction_0(q5, q4, q3, q2, c0)
    call weno_seven_reconstruction_1(q6, q5, q4, q3, c1)
    call weno_seven_reconstruction_2(q7, q6, q5, q4, c2)
    call weno_seven_reconstruction_3(q8, q7, q6, q5, c3)
    if (velocity_smoothing == .true.) then
      call weno_seven_weight_0(u5, u4, u3, u2, b0)
      call weno_seven_weight_1(u6, u5, u4, u3, b1)
      call weno_seven_weight_2(u7, u6, u5, u4, b2)
      call weno_seven_weight_3(u8, u7, u6, u5, b3)
    else
      call weno_seven_weight_0(q5, q4, q3, q2, b0)
      call weno_seven_weight_1(q6, q5, q4, q3, b1)
      call weno_seven_weight_2(q7, q6, q5, q4, b2)
      call weno_seven_weight_3(q8, q7, q6, q5, b3)
    endif
  endif

  tau = abs(b0 + 3 * b1 - 3 * b2 - b3)
  w0  = 4./35.  * (1 + (tau / (b0 + 1e-20))**2)
  w1  = 18./35. * (1 + (tau / (b1 + 1e-20))**2)
  w2  = 12./35. * (1 + (tau / (b2 + 1e-20))**2)
  w3  = 1./35.  * (1 + (tau / (b3 + 1e-20))**2)

  s = 1. / (w0 + w1 + w2 + w3)
  w0 = w0 * s
  w1 = w1 * s
  w2 = w2 * s
  w3 = w3 * s

  qr = w0 * c0 + w1 * c1 + w2 * c2 + w3 * c3

end subroutine weno_seven_reconstruction

subroutine weno_seven_weight_0(q0, q1, q2, q3, w0)
  real, intent(in) :: q0,q1,q2,q3
  real, intent(inout) :: w0

  w0 = q0 * (2.107 * q0 - 9.402 * q1 + 7.042 * q2 - 1.854 * q3) + q1 * (11.003 * q1 - 17.246 * q2 + 4.642 * q3) + &
       q2 * (7.043 * q2 - 3.882 * q3) + 0.547 * q3 * q3
       
end subroutine weno_seven_weight_0

subroutine weno_seven_weight_1(q0, q1, q2, q3, w1)
  real, intent(in) :: q0,q1,q2,q3
  real, intent(inout) :: w1

  w1 = q0 * (0.547 * q0 - 2.522 * q1 + 1.922 * q2 - 0.494 * q3) + q1 * (3.443 * q1 - 5.966 * q2 + 1.602 * q3) + &
       q2 * (2.843 * q2 - 1.642 * q3) + 0.267 * q3 * q3
       
end subroutine weno_seven_weight_1
    
subroutine weno_seven_weight_2(q0, q1, q2, q3, w2)
  real, intent(in) :: q0,q1,q2,q3
  real, intent(inout) :: w2

  w2 = q0 * (0.267 * q0 - 1.642 * q1 + 1.602 * q2 - 0.494 * q3) + q1 * (2.843 * q1 - 5.966 * q2 + 1.922 * q3) + &
       q2 * (3.443 * q2 - 2.522 * q3) + 0.547 * q3 * q3
       
end subroutine weno_seven_weight_2

subroutine weno_seven_weight_3(q0, q1, q2, q3, w3)
  real, intent(in) :: q0,q1,q2,q3
  real, intent(inout) :: w3

  w3 = q0 * (0.547  * q0 - 3.882 * q1 + 4.642 * q2 - 1.854 * q3) + q1 * (7.043 * q1 - 17.246 * q2 + 7.042 * q3) + &
       q2 * (11.003 * q2 - 9.402 * q3) + 2.107 * q3 * q3
       
end subroutine weno_seven_weight_3

subroutine weno_seven_reconstruction_0(q0, q1, q2, q3, p0)
  real, intent(in) :: q0,q1,q2,q3
  real, intent(inout) :: p0

  p0 = (6*q0 + 26*q1 - 10*q2 + 2*q3) / 24.

end subroutine weno_seven_reconstruction_0

subroutine weno_seven_reconstruction_1(q0, q1, q2, q3, p1)
  real, intent(in) :: q0,q1,q2,q3
  real, intent(inout) :: p1

  p1 = (-2*q0 + 14*q1 + 14*q2 - 2*q3) / 24.

end subroutine weno_seven_reconstruction_1

subroutine weno_seven_reconstruction_2(q0, q1, q2, q3, p2)
  real, intent(in) :: q0,q1,q2,q3
  real, intent(inout) :: p2

  p2 = (2*q0 - 10*q1 + 26*q2 + 6*q3) / 24.

end subroutine weno_seven_reconstruction_2

subroutine weno_seven_reconstruction_3(q0, q1, q2, q3, p3)
  real, intent(in) :: q0,q1,q2,q3
  real, intent(inout) :: p3

  p3 = (-6*q0 + 26*q1 - 46*q2 + 50*q3) / 24.

end subroutine weno_seven_reconstruction_3

!> Calculates the acceleration due to the gradient of kinetic energy.
subroutine gradKE(u, v, h, KE, KEx, KEy, k, OBC, G, GV, US, CS)
  type(ocean_grid_type),                      intent(in)  :: G   !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)  :: GV  !< Vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)  :: u   !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)  :: v   !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h   !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G) ,SZJ_(G) ),         intent(out) :: KE  !< Kinetic energy per unit mass [L2 T-2 ~> m2 s-2]
  real, dimension(SZIB_(G),SZJ_(G) ),         intent(out) :: KEx !< Zonal acceleration due to kinetic
                                                                 !! energy gradient [L T-2 ~> m s-2]
  real, dimension(SZI_(G) ,SZJB_(G)),         intent(out) :: KEy !< Meridional acceleration due to kinetic
                                                                 !! energy gradient [L T-2 ~> m s-2]
  integer,                                    intent(in)  :: k   !< Layer number to calculate for
  type(ocean_OBC_type),                       pointer     :: OBC !< Open boundary control structure
  type(unit_scale_type),                      intent(in)  :: US  !< A dimensional unit scaling type
  type(CoriolisAdv_CS),                       intent(in)  :: CS  !< Control structure for MOM_CoriolisAdv
  ! Local variables
  real :: um, un, up, uq, vm, vn, vp, vq         ! Temporary variables [L T-1 ~> m s-1].
  real :: um2, up2, vm2, vp2     ! Temporary variables [L2 T-2 ~> m2 s-2].
  real :: um2a, up2a, vm2a, vp2a ! Temporary variables [L4 T-2 ~> m4 s-2].
  real :: third_order_u, third_order_v  ! Product of mask values to determine the boundary
  real :: fifth_order_u, fifth_order_v  ! Product of mask values to determine the boundary
  real :: seventh_order_u, seventh_order_v  ! Product of mask values to determine the boundary
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, n

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB


  ! Calculate KE (Kinetic energy for use in the -grad(KE) acceleration term).
  if (CS%KE_Scheme == KE_ARAKAWA) then
    ! The following calculation of Kinetic energy includes the metric terms
    ! identified in Arakawa & Lamb 1982 as important for KE conservation.  It
    ! also includes the possibility of partially-blocked tracer cell faces.
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      KE(i,j) = ( ( (G%areaCu( I ,j)*(u( I ,j,k)*u( I ,j,k))) + &
                    (G%areaCu(I-1,j)*(u(I-1,j,k)*u(I-1,j,k))) ) + &
                  ( (G%areaCv(i, J )*(v(i, J ,k)*v(i, J ,k))) + &
                    (G%areaCv(i,J-1)*(v(i,J-1,k)*v(i,J-1,k))) ) )*0.25*G%IareaT(i,j)
    enddo ; enddo
  elseif (CS%KE_Scheme == KE_SIMPLE_GUDONOV) then
    ! The following discretization of KE is based on the one-dimensional Gudonov
    ! scheme which does not take into account any geometric factors
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      up = 0.5*( u(I-1,j,k) + ABS( u(I-1,j,k) ) ) ; up2 = up*up
      um = 0.5*( u( I ,j,k) - ABS( u( I ,j,k) ) ) ; um2 = um*um
      vp = 0.5*( v(i,J-1,k) + ABS( v(i,J-1,k) ) ) ; vp2 = vp*vp
      vm = 0.5*( v(i, J ,k) - ABS( v(i, J ,k) ) ) ; vm2 = vm*vm
      KE(i,j) = ( max(up2,um2) + max(vp2,vm2) ) *0.5
    enddo ; enddo
  elseif (CS%KE_Scheme == KE_GUDONOV) then
    ! The following discretization of KE is based on the one-dimensional Gudonov
    ! scheme but has been adapted to take horizontal grid factors into account
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      up = 0.5*( u(I-1,j,k) + ABS( u(I-1,j,k) ) ) ; up2a = up*up*G%areaCu(I-1,j)
      um = 0.5*( u( I ,j,k) - ABS( u( I ,j,k) ) ) ; um2a = um*um*G%areaCu( I ,j)
      vp = 0.5*( v(i,J-1,k) + ABS( v(i,J-1,k) ) ) ; vp2a = vp*vp*G%areaCv(i,J-1)
      vm = 0.5*( v(i, J ,k) - ABS( v(i, J ,k) ) ) ; vm2a = vm*vm*G%areaCv(i, J )
      KE(i,j) = ( max(um2a,up2a) + max(vm2a,vp2a) )*0.5*G%IareaT(i,j)
    enddo ; enddo
  elseif (CS%KE_Scheme == KE_UP3) then
    ! The following discretization of KE is based on the one-dimensional third-order
    ! upwind scheme which does not take horizontal grid factors into account
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    ! compute the masking to make sure that inland values are not used
      third_order_u = (G%mask2dCu(I-2,j) * G%mask2dCu(I-1,j)* &
                     G%mask2dCu(I,j) * G%mask2dCu(I+1,j))

      if (third_order_u == 1) then
        up = (-u(I-2,j,k) + 7*u(I-1,j,k) + 7*u(I,j,k) - u(I+1,j,k))/12.
        call UP3_limiter_reconstruction(u(I-2,j,k), u(I-1,j,k),&
                u(I,j,k), u(I+1,j,k), up, um)
!        if (up>0.) then
!          um = (-u(I-2,j,k) + 5*u(I-1,j,k) + 2*u(I,j,k))/6.
!        elseif (up<0.) then
!          um = (2*u(I-1,j,k) + 5*u(I,j,k) - u(I+1,j,k))/6.
!        else
!          um = up
!        endif
      else
        up = (u(I-1,j,k) + u(I,j,k))*0.5
        if (up>0.) then
          um = u(I-1,j,k)
        elseif (up<0.) then
          um = u(I,j,k)
        else
          um = up
        endif
      endif

      third_order_v = (G%mask2dCv(i,J-2) * G%mask2dCv(i,J-1)* &
                     G%mask2dCv(i,J) * G%mask2dCv(i,J+1))  
      if (third_order_v ==1) then
        vp = (-v(i,J-2,k) + 7*v(i,J-1,k) + 7*v(i,J,k) - v(i,J+1,k))/12.
        call UP3_limiter_reconstruction(v(i,J-2,k), v(i,J-1,k),&
                v(i,J,k), v(i,J+1,k), vp, vm)
!        if (vp>0.) then
!          vm = (-v(i,J-2,k) + 5*v(i,J-1,k) + 2*v(i,J,k))/6.
!        elseif (vp<0.) then
!          vm = (2*v(i,J-1,k) + 5*v(i,J,k) - v(i,J+1,k))/6.
!        else
!          vm = vp
!        endif
      else
        vp = (v(i,J-1,k) + v(i,J,k))*0.5
        if (vp>0.) then
          vm = v(i,J-1,k)
        elseif (vp<0.) then
          vm = v(i,J,k)
        else
          vm = vp
        endif
      endif

      KE(i,j) = ( um*um + vm*vm )*0.5
    enddo ; enddo
  elseif (CS%KE_Scheme == KE_wenovi_7th) then
    ! The following discretization of KE is based on the one-dimensional third-order
    ! upwind scheme which does not take horizontal grid factors into account
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    ! compute the masking to make sure that inland values are not used

      third_order_u = (G%mask2dCu(I-2,j) * G%mask2dCu(I-1,j) * G%mask2dCu(I,j) * & 
                     G%mask2dCu(I+1,j)) 

      fifth_order_u   = third_order_u * G%mask2dCu(I-3,j) * G%mask2dCu(I+2,j) 
      seventh_order_u = fifth_order_u * G%mask2dCu(I-4,j) * G%mask2dCu(I+3,j)*0

      ! compute the masking to make sure that inland values are not used
      if (seventh_order_u == 1) then
        ! all values are valid, we use seventh order reconstruction
        call weno_seven_reconstruction(u(I-4,j,k),u(I-3,j,k),u(I-2,j,k),u(I-1,j,k), & 
                                       u(I,j,k),u(I+1,j,k),u(I+2,j,k),u(I+3,j,k), & 
                                       1., 1., 1., 1., 1., 1., 1., 1., 1., up, CS%weno_velocity_smooth)
        call weno_seven_reconstruction(u(I-4,j,k),u(I-3,j,k),u(I-2,j,k),u(I-1,j,k), & 
                                       u(I,j,k),u(I+1,j,k),u(I+2,j,k),u(I+3,j,k), & 
                                       1., 1., 1., 1., 1., 1., 1., 1., -1., uq, CS%weno_velocity_smooth)
      elseif (fifth_order_u == 1) then
        ! all values are valid, we use fifth order reconstruction
        call weno_five_reconstruction(u(I-3,j,k),u(I-2,j,k),u(I-1,j,k), &
                                      u(I,j,k),  u(I+1,j,k),u(I+2,j,k), & 
                                       1., 1., 1., 1., 1., 1., 1., up, CS%weno_velocity_smooth)
        call weno_five_reconstruction(u(I-3,j,k),u(I-2,j,k),u(I-1,j,k), &
                                      u(I,j,k),  u(I+1,j,k),u(I+2,j,k), & 
                                       1., 1., 1., 1., 1., 1., -1., uq, CS%weno_velocity_smooth)
      elseif (third_order_u == 1) then
        ! only the middle values are valid, we use third order reconstruction
        call weno_three_reconstruction(u(I-2,j,k),u(I-1,j,k),u(I,j,k),u(I+1,j,k), & 
                                       1., 1., 1., 1., 1., up, CS%weno_velocity_smooth)
        call weno_three_reconstruction(u(I-2,j,k),u(I-1,j,k),u(I,j,k),u(I+1,j,k), & 
                                       1., 1., 1., 1., -1., uq, CS%weno_velocity_smooth)
      else ! Upwind first order
        up = u(I-1,j,k)
        uq = u(I,j,k)
      endif
      un = (up + uq)*0.5
      !um = (SIGN(1., un)*0.5*up - SIGN(1., un)*0.5*uq) + un ! um becomes up if un is 0
      if (un>0.) then
        um = up
      elseif (un<0.) then
        um = uq
      else
        um = 0.
      endif


    ! compute the masking to make sure that inland values are not used

      third_order_v = (G%mask2dCv(i,J-2) * G%mask2dCv(i,J-1) * G%mask2dCv(i,J) * & 
                     G%mask2dCv(i,J+1)) 

      fifth_order_v   = third_order_v * G%mask2dCv(i,J-3) * G%mask2dCv(i,J+2) 
      seventh_order_v = fifth_order_v * G%mask2dCv(i,J-4) * G%mask2dCv(i,J+3)*0

      ! compute the masking to make sure that inland values are not used
      if (seventh_order_v == 1) then
        ! all values are valid, we use seventh order reconstruction
        call weno_seven_reconstruction(v(i,J-4,k),v(i,J-3,k),v(i,J-2,k),v(i,J-1,k), & 
                                       v(i,J,k),v(i,J+1,k),v(i,J+2,k),v(i,J+3,k), & 
                                       1., 1., 1., 1., 1., 1., 1., 1., 1., vp, CS%weno_velocity_smooth)
        call weno_seven_reconstruction(v(i,J-4,k),v(i,J-3,k),v(i,J-2,k),v(i,J-1,k), & 
                                       v(i,J,k),v(i,J+1,k),v(i,J+2,k),v(i,J+3,k), & 
                                       1., 1., 1., 1., 1., 1., 1., 1., -1., vq, CS%weno_velocity_smooth)
      elseif (fifth_order_v == 1) then
        ! all values are valid, we use fifth order reconstruction
        call weno_five_reconstruction(v(i,J-3,k),v(i,J-2,k),v(i,J-1,k), &
                                      v(i,J,k),  v(i,J+1,k),v(i,J+2,k), & 
                                       1., 1., 1., 1., 1., 1., 1., vp, CS%weno_velocity_smooth)
        call weno_five_reconstruction(v(i,J-3,k),v(i,J-2,k),v(i,J-1,k), &
                                      v(i,J,k),  v(i,J+1,k),v(i,J+2,k), & 
                                       1., 1., 1., 1., 1., 1., -1., vq, CS%weno_velocity_smooth)
      elseif (third_order_v == 1) then
        ! only the middle values are valid, we use third order reconstruction
        call weno_three_reconstruction(v(i,J-2,k),v(i,J-1,k),v(i,J,k),v(i,J+1,k), & 
                                       1., 1., 1., 1., 1., vp, CS%weno_velocity_smooth)
        call weno_three_reconstruction(v(i,J-2,k),v(i,J-1,k),v(i,J,k),v(i,J+1,k), & 
                                       1., 1., 1., 1., -1., vq, CS%weno_velocity_smooth)
      else ! Upwind first order
        vp = v(i,J-1,k)
        vq = v(i,J,k)
      endif
      vn = (vp + vq)*0.5
     ! vm = (SIGN(1., vn)*0.5*vp - SIGN(1., vn)*0.5*vq) + vn ! vm becomes vp if vn is 0
      if (vn>0.) then
        vm = vp
      elseif (vn<0.) then
        vm = vq
      else
        vm = 0.
      endif

      KE(i,j) = ( um*um + vm*vm )*0.5
    enddo ; enddo
    

  endif

  ! Term - d(KE)/dx.
  do j=js,je ; do I=Isq,Ieq
    KEx(I,j) = (KE(i+1,j) - KE(i,j)) * G%IdxCu(I,j)
  enddo ; enddo

  ! Term - d(KE)/dy.
  do J=Jsq,Jeq ; do i=is,ie
    KEy(i,J) = (KE(i,j+1) - KE(i,j)) * G%IdyCv(i,J)
  enddo ; enddo

  if (associated(OBC)) then
    do n=1,OBC%number_of_segments
      if (OBC%segment(n)%is_N_or_S) then
        do i=OBC%segment(n)%HI%isd,OBC%segment(n)%HI%ied
          KEy(i,OBC%segment(n)%HI%JsdB) = 0.
        enddo
      elseif (OBC%segment(n)%is_E_or_W) then
        do j=OBC%segment(n)%HI%jsd,OBC%segment(n)%HI%jed
          KEx(OBC%segment(n)%HI%IsdB,j) = 0.
        enddo
      endif
    enddo
  endif

end subroutine gradKE

!> Initializes the control structure for MOM_CoriolisAdv
subroutine CoriolisAdv_init(Time, G, GV, US, param_file, diag, AD, CS)
  type(time_type), target, intent(in)    :: Time !< Current model time
  type(ocean_grid_type),   intent(in)    :: G    !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< Vertical grid structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Runtime parameter handles
  type(diag_ctrl), target, intent(inout) :: diag !< Diagnostics control structure
  type(accel_diag_ptrs),   target, intent(inout) :: AD !< Storage for acceleration diagnostics
  type(CoriolisAdv_CS),    intent(inout) :: CS   !< Control structure for MOM_CoriolisAdv
  ! Local variables
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_CoriolisAdv" ! This module's name.
  character(len=20)  :: tmpstr
  character(len=400) :: mesg
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = GV%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  CS%initialized = .true.
  CS%diag => diag ; CS%Time => Time

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "NOSLIP", CS%no_slip, &
                 "If true, no slip boundary conditions are used; otherwise "//&
                 "free slip boundary conditions are assumed. The "//&
                 "implementation of the free slip BCs on a C-grid is much "//&
                 "cleaner than the no slip BCs. The use of free slip BCs "//&
                 "is strongly encouraged, and no slip BCs are not used with "//&
                 "the biharmonic viscosity.", default=.false.)

  call get_param(param_file, mdl, "CORIOLIS_EN_DIS", CS%Coriolis_En_Dis, &
                 "If true, two estimates of the thickness fluxes are used "//&
                 "to estimate the Coriolis term, and the one that "//&
                 "dissipates energy relative to the other one is used.", &
                 default=.false.)

  ! Set %Coriolis_Scheme
  ! (Select the baseline discretization for the Coriolis term)
  call get_param(param_file, mdl, "CORIOLIS_SCHEME", tmpstr, &
                 "CORIOLIS_SCHEME selects the discretization for the "//&
                 "Coriolis terms. Valid values are: \n"//&
                 "\t SADOURNY75_ENERGY - Sadourny, 1975; energy cons. \n"//&
                 "\t ARAKAWA_HSU90     - Arakawa & Hsu, 1990 \n"//&
                 "\t SADOURNY75_ENSTRO - Sadourny, 1975; enstrophy cons. \n"//&
                 "\t ARAKAWA_LAMB81    - Arakawa & Lamb, 1981; En. + Enst.\n"//&
                 "\t ARAKAWA_LAMB_BLEND - A blend of Arakawa & Lamb with \n"//&
                 "\t                      Arakawa & Hsu and Sadourny energy \n"//&
                 "\t SECOND_ENSTRO - 2nd-order enstrophy cons. \n"//&
                 "\t UP3_ENSTRO - 3rd-order enstrophy cons. \n"//&
                 "\t UP3_PV_ENSTRO - 3rd-order PV enstrophy cons. \n"//&
                 "\t UP3_SPLIT - 3rd-order enstrophy cons. \n"//&
                 "\t CEN4_ENSTRO - 4th-order enstrophy cons. \n"//&
                 "\t WENOVI_5TH_ENSTRO - 5th-order enstrophy cons. \n"//&
                 "\t WENOVI_7TH_SPLIT - 7th-order enstrophy cons. \n"//&
                 "\t WENOVI_7TH_ENSTRO - 7th-order enstrophy cons.",&

                 default=SADOURNY75_ENERGY_STRING)
  tmpstr = uppercase(tmpstr)
  select case (tmpstr)
    case (SADOURNY75_ENERGY_STRING)
      CS%Coriolis_Scheme = SADOURNY75_ENERGY
    case (ARAKAWA_HSU_STRING)
      CS%Coriolis_Scheme = ARAKAWA_HSU90
    case (SADOURNY75_ENSTRO_STRING)
      CS%Coriolis_Scheme = SADOURNY75_ENSTRO
    case (ARAKAWA_LAMB_STRING)
      CS%Coriolis_Scheme = ARAKAWA_LAMB81
    case (AL_BLEND_STRING)
      CS%Coriolis_Scheme = AL_BLEND
    case (ROBUST_ENSTRO_STRING)
      CS%Coriolis_Scheme = ROBUST_ENSTRO
      CS%Coriolis_En_Dis = .false.
    case (SECOND_ENSTRO_STRING)
      CS%Coriolis_Scheme = SECOND_ENSTRO
    case (UP3_ENSTRO_STRING)
      CS%Coriolis_Scheme = UP3_ENSTRO
    case (UP3_PV_ENSTRO_STRING)
      CS%Coriolis_Scheme = UP3_PV_ENSTRO
    case (UP3_SPLIT_STRING)
      CS%Coriolis_Scheme = UP3_split
    case (CEN4_ENSTRO_STRING)
      CS%Coriolis_Scheme = CEN4_ENSTRO
    case (WENOVI_7TH_ENSTRO_STRING)
      CS%Coriolis_Scheme = wenovi_7th_ENSTRO      
    case (WENOVI_7TH_SPLIT_STRING)
      CS%Coriolis_Scheme = wenovi_7th_split      
    case (WENOVI_5TH_ENSTRO_STRING)
      CS%Coriolis_Scheme = wenovi_5th_ENSTRO
    case default
      call MOM_mesg('CoriolisAdv_init: Coriolis_Scheme ="'//trim(tmpstr)//'"', 0)
      call MOM_error(FATAL, "CoriolisAdv_init: Unrecognized setting "// &
            "#define CORIOLIS_SCHEME "//trim(tmpstr)//" found in input file.")
  end select
  if (CS%Coriolis_Scheme == wenovi_7th_ENSTRO .or. &
      CS%Coriolis_Scheme == wenovi_7th_split .or. &
      CS%Coriolis_Scheme == wenovi_5th_ENSTRO) then
    call get_param(param_file, mdl, "WENO_VELOCITY_SMOOTH", CS%weno_velocity_smooth, &
            "If true, use velocity to compute weighting for WENO. ", &
                  default=.false.)
  endif
  if (CS%Coriolis_Scheme == AL_BLEND) then
    call get_param(param_file, mdl, "CORIOLIS_BLEND_WT_LIN", CS%wt_lin_blend, &
                 "A weighting value for the ratio of inverse thicknesses, "//&
                 "beyond which the blending between Sadourny Energy and "//&
                 "Arakawa & Hsu goes linearly to 0 when CORIOLIS_SCHEME "//&
                 "is ARAWAKA_LAMB_BLEND. This must be between 1 and 1e-16.", &
                 units="nondim", default=0.125)
    call get_param(param_file, mdl, "CORIOLIS_BLEND_F_EFF_MAX", CS%F_eff_max_blend, &
                 "The factor by which the maximum effective Coriolis "//&
                 "acceleration from any point can be increased when "//&
                 "blending different discretizations with the "//&
                 "ARAKAWA_LAMB_BLEND Coriolis scheme.  This must be "//&
                 "greater than 2.0 (the max value for Sadourny energy).", &
                 units="nondim", default=4.0)
    CS%wt_lin_blend = min(1.0, max(CS%wt_lin_blend,1e-16))
    if (CS%F_eff_max_blend < 2.0) call MOM_error(WARNING, "CoriolisAdv_init: "//&
           "CORIOLIS_BLEND_F_EFF_MAX should be at least 2.")
  endif

  mesg = "If true, the Coriolis terms at u-points are bounded by "//&
         "the four estimates of (f+rv)v from the four neighboring "//&
         "v-points, and similarly at v-points."
  if (CS%Coriolis_En_Dis .and. (CS%Coriolis_Scheme == SADOURNY75_ENERGY)) then
    mesg = trim(mesg)//"  This option is "//&
                 "always effectively false with CORIOLIS_EN_DIS defined and "//&
                 "CORIOLIS_SCHEME set to "//trim(SADOURNY75_ENERGY_STRING)//"."
  else
    mesg = trim(mesg)//"  This option would "//&
                 "have no effect on the SADOURNY Coriolis scheme if it "//&
                 "were possible to use centered difference thickness fluxes."
  endif
  call get_param(param_file, mdl, "BOUND_CORIOLIS", CS%bound_Coriolis, mesg, &
                 default=.false.)
  if ((CS%Coriolis_En_Dis .and. (CS%Coriolis_Scheme == SADOURNY75_ENERGY)) .or. &
      (CS%Coriolis_Scheme == ROBUST_ENSTRO)) CS%bound_Coriolis = .false.

  ! Set KE_Scheme (selects discretization of KE)
  call get_param(param_file, mdl, "KE_SCHEME", tmpstr, &
                 "KE_SCHEME selects the discretization for acceleration "//&
                 "due to the kinetic energy gradient. Valid values are: \n"//&
                 "\t KE_ARAKAWA, KE_SIMPLE_GUDONOV, KE_GUDONOV, KE_UP3, &
                 KE_WENOVI_7TH", &
                 default=KE_ARAKAWA_STRING)
  tmpstr = uppercase(tmpstr)
  select case (tmpstr)
    case (KE_ARAKAWA_STRING); CS%KE_Scheme = KE_ARAKAWA
    case (KE_SIMPLE_GUDONOV_STRING); CS%KE_Scheme = KE_SIMPLE_GUDONOV
    case (KE_GUDONOV_STRING); CS%KE_Scheme = KE_GUDONOV
    case (KE_UP3_STRING); CS%KE_Scheme = KE_UP3
    case (KE_WENOVI_7TH_STRING); CS%KE_Scheme = KE_wenovi_7th
    case default
      call MOM_mesg('CoriolisAdv_init: KE_Scheme ="'//trim(tmpstr)//'"', 0)
      call MOM_error(FATAL, "CoriolisAdv_init: "// &
               "#define KE_SCHEME "//trim(tmpstr)//" in input file is invalid.")
  end select

  if (CS%Coriolis_Scheme == UP3_ENSTRO .or. &
      CS%Coriolis_Scheme == UP3_PV_ENSTRO .or. CS%Coriolis_Scheme == UP3_split .or. &
      CS%KE_Scheme == KE_UP3) then
    call get_param(param_file, mdl, "UP3_USE_LIMITER", CS%UP3_use_limiter, &
            "If true, use flux limiter for UP3 scheme ", &
                  default=.false.)
  endif

  ! Set PV_Adv_Scheme (selects discretization of PV advection)
  call get_param(param_file, mdl, "PV_ADV_SCHEME", tmpstr, &
                 "PV_ADV_SCHEME selects the discretization for PV "//&
                 "advection. Valid values are: \n"//&
                 "\t PV_ADV_CENTERED - centered (aka Sadourny, 75) \n"//&
                 "\t PV_ADV_UPWIND1  - upwind, first order", &
                 default=PV_ADV_CENTERED_STRING)
  select case (uppercase(tmpstr))
    case (PV_ADV_CENTERED_STRING)
      CS%PV_Adv_Scheme = PV_ADV_CENTERED
    case (PV_ADV_UPWIND1_STRING)
      CS%PV_Adv_Scheme = PV_ADV_UPWIND1
    case default
      call MOM_mesg('CoriolisAdv_init: PV_Adv_Scheme ="'//trim(tmpstr)//'"', 0)
      call MOM_error(FATAL, "CoriolisAdv_init: "// &
                     "#DEFINE PV_ADV_SCHEME in input file is invalid.")
  end select

  CS%id_rv = register_diag_field('ocean_model', 'RV', diag%axesBL, Time, &
     'Relative Vorticity', 's-1', conversion=US%s_to_T)

  CS%id_PV = register_diag_field('ocean_model', 'PV', diag%axesBL, Time, &
     'Potential Vorticity', 'm-1 s-1', conversion=GV%m_to_H*US%s_to_T)

  CS%id_gKEu = register_diag_field('ocean_model', 'gKEu', diag%axesCuL, Time, &
     'Zonal Acceleration from Grad. Kinetic Energy', 'm s-2', conversion=US%L_T2_to_m_s2)

  CS%id_gKEv = register_diag_field('ocean_model', 'gKEv', diag%axesCvL, Time, &
     'Meridional Acceleration from Grad. Kinetic Energy', 'm s-2', conversion=US%L_T2_to_m_s2)

  CS%id_rvxu = register_diag_field('ocean_model', 'rvxu', diag%axesCvL, Time, &
     'Meridional Acceleration from Relative Vorticity', 'm s-2', conversion=US%L_T2_to_m_s2)

  CS%id_rvxv = register_diag_field('ocean_model', 'rvxv', diag%axesCuL, Time, &
     'Zonal Acceleration from Relative Vorticity', 'm s-2', conversion=US%L_T2_to_m_s2)

  CS%id_CAuS = register_diag_field('ocean_model', 'CAu_Stokes', diag%axesCuL, Time, &
     'Zonal Acceleration from Stokes Vorticity', 'm s-2', conversion=US%L_T2_to_m_s2)
  ! add to AD

  CS%id_CAvS = register_diag_field('ocean_model', 'CAv_Stokes', diag%axesCvL, Time, &
     'Meridional Acceleration from Stokes Vorticity', 'm s-2', conversion=US%L_T2_to_m_s2)
  ! add to AD

  !CS%id_hf_gKEu = register_diag_field('ocean_model', 'hf_gKEu', diag%axesCuL, Time, &
  !   'Fractional Thickness-weighted Zonal Acceleration from Grad. Kinetic Energy', &
  !   'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  CS%id_hf_gKEu_2d = register_diag_field('ocean_model', 'hf_gKEu_2d', diag%axesCu1, Time, &
     'Depth-sum Fractional Thickness-weighted Zonal Acceleration from Grad. Kinetic Energy', &
     'm s-2', conversion=US%L_T2_to_m_s2)

  !CS%id_hf_gKEv = register_diag_field('ocean_model', 'hf_gKEv', diag%axesCvL, Time, &
  !   'Fractional Thickness-weighted Meridional Acceleration from Grad. Kinetic Energy', &
  !   'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  CS%id_hf_gKEv_2d = register_diag_field('ocean_model', 'hf_gKEv_2d', diag%axesCv1, Time, &
     'Depth-sum Fractional Thickness-weighted Meridional Acceleration from Grad. Kinetic Energy', &
     'm s-2', conversion=US%L_T2_to_m_s2)

  CS%id_h_gKEu = register_diag_field('ocean_model', 'h_gKEu', diag%axesCuL, Time, &
     'Thickness Multiplied Zonal Acceleration from Grad. Kinetic Energy', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  CS%id_intz_gKEu_2d = register_diag_field('ocean_model', 'intz_gKEu_2d', diag%axesCu1, Time, &
     'Depth-integral of Zonal Acceleration from Grad. Kinetic Energy', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)

  CS%id_h_gKEv = register_diag_field('ocean_model', 'h_gKEv', diag%axesCvL, Time, &
     'Thickness Multiplied Meridional Acceleration from Grad. Kinetic Energy', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  CS%id_intz_gKEv_2d = register_diag_field('ocean_model', 'intz_gKEv_2d', diag%axesCv1, Time, &
     'Depth-integral of Meridional Acceleration from Grad. Kinetic Energy', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)

  !CS%id_hf_rvxu = register_diag_field('ocean_model', 'hf_rvxu', diag%axesCvL, Time, &
  !   'Fractional Thickness-weighted Meridional Acceleration from Relative Vorticity', &
  !   'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  CS%id_hf_rvxu_2d = register_diag_field('ocean_model', 'hf_rvxu_2d', diag%axesCv1, Time, &
     'Depth-sum Fractional Thickness-weighted Meridional Acceleration from Relative Vorticity', &
     'm s-2', conversion=US%L_T2_to_m_s2)

  !CS%id_hf_rvxv = register_diag_field('ocean_model', 'hf_rvxv', diag%axesCuL, Time, &
  !   'Fractional Thickness-weighted Zonal Acceleration from Relative Vorticity', &
  !   'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  CS%id_hf_rvxv_2d = register_diag_field('ocean_model', 'hf_rvxv_2d', diag%axesCu1, Time, &
     'Depth-sum Fractional Thickness-weighted Zonal Acceleration from Relative Vorticity', &
     'm s-2', conversion=US%L_T2_to_m_s2)

  CS%id_h_rvxu = register_diag_field('ocean_model', 'h_rvxu', diag%axesCvL, Time, &
     'Thickness Multiplied Meridional Acceleration from Relative Vorticity', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  CS%id_intz_rvxu_2d = register_diag_field('ocean_model', 'intz_rvxu_2d', diag%axesCv1, Time, &
     'Depth-integral of Meridional Acceleration from Relative Vorticity', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)

  CS%id_h_rvxv = register_diag_field('ocean_model', 'h_rvxv', diag%axesCuL, Time, &
     'Thickness Multiplied Zonal Acceleration from Relative Vorticity', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  CS%id_intz_rvxv_2d = register_diag_field('ocean_model', 'intz_rvxv_2d', diag%axesCu1, Time, &
     'Depth-integral of Fractional Thickness-weighted Zonal Acceleration from Relative Vorticity', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)

  ! Allocate memory needed for the diagnostics that have been enabled.
  if ((CS%id_gKEu > 0) .or. (CS%id_hf_gKEu_2d > 0) .or. &
    ! (CS%id_hf_gKEu > 0) .or. &
      (CS%id_h_gKEu > 0) .or. (CS%id_intz_gKEu_2d > 0)) then
    call safe_alloc_ptr(AD%gradKEu, IsdB, IedB, jsd, jed, nz)
  endif
  if ((CS%id_gKEv > 0) .or. (CS%id_hf_gKEv_2d > 0) .or. &
    ! (CS%id_hf_gKEv > 0) .or. &
      (CS%id_h_gKEv > 0) .or. (CS%id_intz_gKEv_2d > 0)) then
    call safe_alloc_ptr(AD%gradKEv, isd, ied, JsdB, JedB, nz)
  endif
  if ((CS%id_rvxu > 0) .or. (CS%id_hf_rvxu_2d > 0) .or. &
    ! (CS%id_hf_rvxu > 0) .or. &
      (CS%id_h_rvxu > 0) .or. (CS%id_intz_rvxu_2d > 0)) then
    call safe_alloc_ptr(AD%rv_x_u, isd, ied, JsdB, JedB, nz)
  endif
  if ((CS%id_rvxv > 0) .or. (CS%id_hf_rvxv_2d > 0) .or. &
    ! (CS%id_hf_rvxv > 0) .or. &
      (CS%id_h_rvxv > 0) .or. (CS%id_intz_rvxv_2d > 0)) then
    call safe_alloc_ptr(AD%rv_x_v, IsdB, IedB, jsd, jed, nz)
  endif

  if ((CS%id_hf_gKEv_2d > 0) .or. &
    ! (CS%id_hf_gKEv > 0) .or. (CS%id_hf_rvxu > 0) .or. &
      (CS%id_hf_rvxu_2d > 0)) then
    call safe_alloc_ptr(AD%diag_hfrac_v, isd, ied, JsdB, JedB, nz)
  endif
  if ((CS%id_hf_gKEu_2d > 0) .or. &
    ! (CS%id_hf_gKEu > 0) .or. (CS%id_hf_rvxv > 0) .or. &
      (CS%id_hf_rvxv_2d > 0)) then
    call safe_alloc_ptr(AD%diag_hfrac_u, IsdB, IedB, jsd, jed, nz)
  endif
  if ((CS%id_h_gKEu > 0) .or. (CS%id_intz_gKEu_2d > 0) .or. &
      (CS%id_h_rvxv > 0) .or. (CS%id_intz_rvxv_2d > 0)) then
    call safe_alloc_ptr(AD%diag_hu, IsdB, IedB, jsd, jed, nz)
  endif
  if ((CS%id_h_gKEv > 0) .or. (CS%id_intz_gKEv_2d > 0) .or. &
      (CS%id_h_rvxu > 0) .or. (CS%id_intz_rvxu_2d > 0)) then
    call safe_alloc_ptr(AD%diag_hv, isd, ied, JsdB, JedB, nz)
  endif

end subroutine CoriolisAdv_init

!> Destructor for coriolisadv_cs
subroutine CoriolisAdv_end(CS)
  type(CoriolisAdv_CS), intent(inout) :: CS !< Control structure for MOM_CoriolisAdv
end subroutine CoriolisAdv_end

!> \namespace mom_coriolisadv
!!
!! This file contains the subroutine that calculates the time
!! derivatives of the velocities due to Coriolis acceleration and
!! momentum advection.  This subroutine uses either a vorticity
!! advection scheme from Arakawa and Hsu, Mon. Wea. Rev. 1990, or
!! Sadourny's (JAS 1975) energy conserving scheme.  Both have been
!! modified to use general orthogonal coordinates as described in
!! Arakawa and Lamb, Mon. Wea. Rev. 1981.  Both schemes are second
!! order accurate, and allow for vanishingly small layer thicknesses.
!! The Arakawa and Hsu scheme globally conserves both total energy
!! and potential enstrophy in the limit of nondivergent flow.
!! Sadourny's energy conserving scheme conserves energy if the flow
!! is nondivergent or centered difference thickness fluxes are used.
!!
!! A small fragment of the grid is shown below:
!! \verbatim
!!
!!    j+1  x ^ x ^ x   At x:  q, CoriolisBu
!!    j+1  > o > o >   At ^:  v, CAv, vh
!!    j    x ^ x ^ x   At >:  u, CAu, uh, a, b, c, d
!!    j    > o > o >   At o:  h, KE
!!    j-1  x ^ x ^ x
!!        i-1  i  i+1  At x & ^:
!!           i  i+1    At > & o:
!! \endverbatim
!!
!! The boundaries always run through q grid points (x).

end module MOM_CoriolisAdv
