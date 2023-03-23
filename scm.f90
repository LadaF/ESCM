module kinds
  use iso_fortran_env, only: real64
  
  integer, parameter :: rp = real64 
end module

module grid_variables
  use kinds

  implicit none
  
  integer :: nz, nzw
  
  real(rp) :: dz
  
  type fields
    real(rp), dimension(:), allocatable :: u, v, w, theta
  
    real(rp), dimension(:), allocatable :: e, eps, theta_eps
  
    real(rp), dimension(:), allocatable :: up_wp, vp_wp, thetap_wp, wp_eps, Km, Kh
  contains
    procedure :: fields_sum, fields_multiply, fields_assign
    generic :: operator(+) => fields_sum
    generic :: operator(*) => fields_multiply
    generic :: assignment(=) => fields_assign
    procedure :: allocate => allocate_grid_variables
  end type
  
  real(rp), dimension(:), allocatable :: z, zw
  
 
contains

  subroutine init_grid
    integer :: k
    real(rp) :: z_top, z_surface, dz
    
    z_surface = 0
    z_top = 3000
    
    nz = 150
    nzw = nz-1
  
    dz = (z_top - z_surface) / nz
    
    allocate(z(1:nz))
    allocate(zw(0:nzw+1))
    
    do k = 0, nzw+1
      zw(k) = z_surface + dz * k
    end do
    
    do k = 1, nz
      z(k) = (zw(k) + zw(k-1)) / 2
    end do
    
  end subroutine


  subroutine allocate_grid_variables(u)
    class(fields) :: u
  
 
    allocate(u%u(1:nz))
    allocate(u%v(1:nz))
    allocate(u%w(1:nz))
    allocate(u%theta(1:nz))
    allocate(u%e(0:nz+1))
    allocate(u%eps(0:nz+1))
    allocate(u%up_wp(0:nz+1))
    allocate(u%vp_wp(0:nz+1))
    allocate(u%thetap_wp(0:nz+1))
    allocate(u%Km(0:nz+1))
    allocate(u%Kh(0:nz+1))
  
  end subroutine
  
  function fields_sum(a, b) result(c)
    type(fields) :: c
    class(fields), intent(in) :: a, b
    c%u = a%u + b%u
    c%v = a%v + b%v
    c%theta = a%theta + b%theta
    c%e = a%e + b%e
    c%eps = a%eps + b%eps
  end function
  
  function fields_multiply(a, b) result(c)
    type(fields) :: c
    class(fields), intent(in) :: a
    real(rp), intent(in) :: b
    c%u = a%u * b
    c%v = a%v * b
    c%theta = a%theta * b
    c%e = a%e * b
    c%eps = a%eps * b
  end function
  
  subroutine fields_assign(a, b)
    class(fields), intent(inout) :: a
    class(fields), intent(in) :: b
    a%u = b%u
    a%v = b%v
    a%theta = b%theta
    a%e = b%e
    a%eps = b%eps
    where(a%eps(1:nzw)<1e-6_rp) a%eps(1:nzw)=1e-6_rp
    where(a%e(1:nzw)<1e-4_rp) a%e(1:nzw)=1e-4_rp
  end subroutine
end module


module global_variables
  use kinds
  
  real(rp), parameter :: pi = acos(-1._rp)
  
  real(rp) :: f_Coriolis = 1.114e-4_rp
  
  real(rp) :: g = 9.81_rp
  
  real(rp), dimension(:), allocatable :: u_g, v_g
  
  real(rp) :: z0 = 0.1
  
contains

  subroutine init_geostrophic_windspeed
    use grid_variables
  
    allocate(u_g(1:nz), v_g(1:nz))
    u_g = 12
    v_g = 0
  end subroutine

end module 



module boundary_conditions
  use kinds
  use global_variables
  use grid_variables
  implicit none
  
  
  public set_boundary_conditions
  public init_boundary_conditions

  !else theta is fixed
  logical :: theta_flux_fixed = .true.

  
  real(rp) :: theta_surf = 300._rp
  real(rp) :: thetap_wp_surf = 0.05
  
  real(rp) :: u_star = 0
  real(rp) :: up_wp_surf, vp_wp_surf
    

contains

  subroutine init_boundary_conditions(theta_flux_fixed_arg, theta_s, thetap_wp_s)
    logical, intent(in) :: theta_flux_fixed_arg
    real(rp), intent(in) :: theta_s, thetap_wp_s
    theta_flux_fixed = theta_flux_fixed_arg
    theta_surf = theta_s
    thetap_wp_surf = thetap_wp_s
  end subroutine
  
  
  subroutine set_boundary_conditions(u, theta_s, thetap_wp_s)
    type(fields) :: u
    real(rp), optional :: theta_s, thetap_wp_s
    
    if ((theta_flux_fixed) .and. (present(thetap_wp_s))) then
      thetap_wp_surf = thetap_wp_s
    else if ((.not.theta_flux_fixed) .and. (present(theta_s))) then
      theta_surf = theta_s
    end if
    
    call compute_surface_fluxes(u)

    u%up_wp(0) = up_wp_surf
    u%vp_wp(0) = vp_wp_surf
    u%thetap_wp(0) = thetap_wp_surf

    u%up_wp(nz) = 0
    u%vp_wp(nz) = 0
    u%thetap_wp(nz) = 0

    u%e(0) = 0
    u%eps(0) = 0
    u%Km(0) = 0
    u%Kh(0) = 0
    
    u%e(nzw+1) = 0
    u%eps(nzw+1) = 0
    u%Km(nzw+1) = 0
    u%Kh(nzw+1) = 0
  end subroutine
  
  
  subroutine compute_surface_fluxes(u)
    type(fields) :: u
    real(rp) :: u_mag
    
    u_mag = hypot(u%u(1), u%v(1))
    
    call compute_u_star_theta_flux(u_mag, u%theta(1), u_star, thetap_wp_surf)

    up_wp_surf = - u%u(1)/u_mag * u_star**2
    vp_wp_surf = - u%v(1)/u_mag * u_star**2
  end subroutine
  

  subroutine compute_u_star_theta_flux(u_mag, theta, u_star, theta_flux)
    real(rp), intent(in) :: u_mag, theta
    real(rp), intent(inout) :: u_star, theta_flux
    if (theta_flux_fixed) then
    
      call fixed_flux_u_star(u_mag, theta, u_star, theta_flux)
      
    else
    
      call fixed_theta_u_star(u_mag, theta, u_star, theta_flux)
    
    end if
  end subroutine
  
  subroutine fixed_flux_u_star(u_mag, theta, u_star, theta_flux)
    real(rp), intent(in) :: u_mag, theta
    real(rp), intent(out) :: u_star, theta_flux
    call WM_MO_Flux_ustar(u_mag, z(1)-zw(0), u_star, z0, theta_flux, theta, g)
  end subroutine
  
  subroutine fixed_theta_u_star(u_mag, theta, u_star, theta_flux)
    real(rp), intent(in) :: u_mag, theta
    real(rp), intent(out) :: u_star, theta_flux
    call WM_MO_Dirichlet_ustar_tfl(u_star, theta_flux, z(1)-zw(0), &
                                   z0, z0, u_mag, theta, theta_surf)
  end subroutine
  
  
  !the following taken from ELMM
  
  
  subroutine WM_MO_Flux_ustar(vel,dist,ustar,z0,temperature_flux,temperature_ref,grav_acc)
    real(rp),intent(inout) :: ustar
    real(rp),parameter  :: eps = 1e-3
    real(rp),intent(in) :: vel,dist,z0,temperature_flux
    real(rp),intent(in) :: temperature_ref,grav_acc
    real(rp) :: ustar2,zL,zL2,Psi
    integer :: i

    if (dist<=z0) then

      stop "first model layer must be above z0"
      
    else

     i = 0

     Psi = 0
     
     do
      i = i+1
      ustar2 = ustar

      ustar = ustar + (vel*0.4_rp/(log(max((dist/z0)-Psi,1E-5_rp))) - ustar) / 2
      
      if (ustar<1E-4) then
       zL = -10000*sign(temperature_flux, 1.0_rp)
      else
       zL2 = Obukhov_zL(ustar,temperature_flux, temperature_ref, grav_acc, dist)
       zL = zL+(zL2-zL)/2
      end if

      Psi = PsiM_MO(zL)

      if  (abs(ustar-ustar2)/max(abs(ustar),1.e-3_rp)<eps) exit

      if (i>=50) then
                  ustar = 0
                  exit
      end if

     end do

    end if

  end subroutine WM_MO_Flux_ustar
  
  pure real(rp) function Obukhov_zL(ustar,temperature_flux,tempref,g,z)
    real(rp),intent(in):: ustar,temperature_flux,tempref,g,z

    Obukhov_zL = z*(0.4_rp*(g/tempref)*temperature_flux)/(-ustar**3)
  end function Obukhov_zL
  
  pure real(rp) function PsiM_MO(zeta)
    real(rp),intent(in):: zeta
    real(rp) :: x

    if (zeta<0) then
      x = (1-15._rp*zeta)**(1/4._rp)
      PsiM_MO = log(((1+x**2)/2._rp)*((1+x)/2._rp)**2)-2._rp*atan(x)+pi/2
    else
      PsiM_MO = - 4.8_rp * zeta
    end if
  end function PsiM_MO


  pure real(rp) function PsiH_MO(zeta)
    real(rp),intent(in):: zeta
    real(rp) :: x

    if (zeta<0) then
      x = (1-15._rp*zeta)**(1/4._rp)
      PsiH_MO = 2._rp*log((1+x**2)/2._rp)
    else
      PsiH_MO = - 7.8_rp * zeta
    end if
  end function PsiH_MO
  
  subroutine WM_MO_Dirichlet_ustar_tfl(ustar,temperature_flux,dist,z0,z0H,vel,temp,temp_surf)
    real(rp),intent(inout) :: ustar,temperature_flux
    real(rp),intent(in) :: vel,dist,z0,z0H,temp,temp_surf
    real(rp),parameter :: eps = 1e-3_rp
    real(rp),parameter :: k_U = 0.4_rp
    real(rp),parameter :: k_T = 0.4_rp
    real(rp) :: zL, zL0, psi_m, psi_h, tempdif
    integer :: i

    tempdif = temp - temp_surf
    
    ustar = vel*k_U/log(dist/z0)

   
    if (dist<=z0) then

      stop "first model layer must be above z0"

    else
       psi_m = 0
       psi_h = 0
       zL0 = -10000._rp
    
       i = 0
       do
         i = i+1
         
         !L does not contain the Karman constant!
         zL =  Obukhov_zL(ustar,temperature_flux,temp,g,dist)
         
         psi_m = PsiM_MO(zL)
         psi_h = PsiH_MO(zL)
         
         ustar = vel * k_U / (log(dist/z0) - psi_m)
         temperature_flux = - tempdif * ustar * k_T / (log(dist/z0H) - psi_h)
         
         if  (i>1 .and. abs(zL-zL0)/max(abs(zL),1.e-4_rp)<eps) exit
         if (i>=50.or.zL>10000) exit
                 
         zl0 = zL
       end do

       if (i>=50.or.zL>10000) then
         ustar = 0
         temperature_flux = 0
       end if
    end if
  end subroutine WM_MO_Dirichlet_ustar_tfl

end module



module initial_conditions
  use kinds
  use grid_variables
  use global_variables
  implicit none

  real(rp) :: theta_mixed = 290
  real(rp) :: gamma_free = 0.004_rp
  real(rp) :: mixed_top = 400
  real(rp) :: delta_theta = 0.01_rp
  
contains

  subroutine set_initial_conditions(u)
    type(fields), intent(inout) :: u
    integer :: k

    u%u = u_g
    u%v = v_g
    
    do k = 1, nz
      if (z(k) < mixed_top) then
        u%theta(k) = theta_mixed
      else
        u%theta(k) = theta_mixed + delta_theta + gamma_free * (zw(k) - mixed_top)
      end if
    end do
  
    do k = 1, nzw
      if (z(k) < mixed_top) then
        u%e(k) = 0.01_rp
        u%eps(k) = 0.0001_rp
      else
        u%e(k) = 0.01_rp
        u%eps(k) = 0.0001_rp
      end if
    end do
  end subroutine

end module




module turbulence_closure
  use kinds
  use grid_variables
  use global_variables
  use boundary_conditions
  
  implicit none
  
! Duynkerke and Driedonks (1987); Stubley and Rooney (1986)
   real,parameter :: c1 = 1.35_rp, c2 = 0.09_rp, c3 = 1.44_rp, c4 = 1.92_rp, c5 = 0.77_rp, c_counter = 5.0
contains

  subroutine compute_turbulence_K(in)
    type(fields), intent(inout) :: in     !only `in`, diagnostics
    integer :: k
    real(rp) :: dth_dz, du_dz, dv_dz, Ri, alpha, shear, w_star
    real(rp) :: pbl_height, gamma_counter
    integer :: pbl_k
    
    !compute the countergradient term
    if (thetap_wp_surf>0) then
      !find out PBL height to apply the countergradient term there and to compute w_star
      do k = 1, nzw
        shear = (du_dz**2 + dv_dz**2)
        if (shear < 1e-5_rp) then
          Ri = 1
        else
          Ri = g/((in%theta(k+1)+in%theta(k))/2) * dth_dz / (du_dz**2 + dv_dz**2)
        end if
        if (Ri>0.25_rp) then
          pbl_height = zw(k)
          pbl_k = k
        end if          
      end do
      
      w_star = (g/in%theta(1)*thetap_wp_surf*pbl_height)
      gamma_counter = c_counter * thetap_wp_surf / w_star / pbl_height
    else
      ! we may consider computing the SBL heigh too, for diagnostics
      gamma_counter = 0
      pbl_k = 0
      pbl_height = 0
    end if
    
    !the main turbulent fluxes
    do k = 1, nzw
      
      in%Km(k) = c2 * in%e(k)**2 / in%eps(k)
      
      
      !TODO: virtual temperature
      dth_dz = (in%theta(k+1) - in%theta(k)) / (z(k+1) - z(k))
      du_dz = (in%u(k+1) - in%u(k)) / (z(k+1) - z(k))
      dv_dz = (in%v(k+1) - in%v(k)) / (z(k+1) - z(k))

      shear = (du_dz**2 + dv_dz**2)
      if (shear < 1e-12_rp.and.(in%theta(k+1)>in%theta(k))) then
        Ri = 1
      else if (shear < 1e-12_rp.and.(in%theta(k+1)<=in%theta(k))) then
        Ri = -10000
      else
        Ri = g/((in%theta(k+1)+in%theta(k))/2) * dth_dz / (du_dz**2 + dv_dz**2)
      end if
       
      if (Ri<0.16_rp) then
        alpha = 1.318_rp * (0.2231_rp - Ri) / (0.2341_rp - Ri)
      else
        alpha = 1.318_rp * (0.2231_rp - 0.16_rp) / (0.2341_rp - 0.16_rp)
      end if
      
      in%Kh(k) = alpha * in%Km(k)
      
      in%up_wp(k) = - in%Km(k) * du_dz
      in%vp_wp(k) = - in%Km(k) * dv_dz
      in%up_wp(k) = - in%Km(k) * du_dz
      
      if (k<=pbl_k) then
        in%thetap_wp(k) = - in%Kh(k) * (dth_dz - gamma_counter)
      else
        in%thetap_wp(k) = - in%Kh(k) * (dth_dz)
      end if
    end do
    
  end subroutine
  

end module


module prognostic_equations
  use kinds
  use grid_variables
  use global_variables
  
  implicit none
  
contains

  subroutine compute_dt_turbulence(out, in)
    use turbulence_closure
    type(fields), intent(inout) :: out
    type(fields), intent(in) :: in
    
    real(rp) :: dth_dz, du_dz, dv_dz, Kh_Pbuoy, Km_shear
    integer :: k
    !values are increments dy_dt * delta_t
    
    associate(u=>in%u, v=>in%v, theta=>in%theta, e=>in%e, eps=>in%eps, Km=>in%Km, Kh=>in%Kh)

    do k = 1, nzw
      dth_dz = (theta(k+1) - theta(k)) / (z(k+1) - z(k))
      du_dz = (u(k+1) - u(k)) / (z(k+1) - z(k))
      dv_dz = (v(k+1) - v(k)) / (z(k+1) - z(k))
      
      Kh_Pbuoy = Kh(k) * g / ((theta(k+1)+theta(k))/2) *  dth_dz
      Km_shear = Km(k) * (du_dz**2 + dv_dz**2)
      
      out%eps(k) =   c3 * eps(k) / e(k) * max( Km_shear -  Kh_Pbuoy, Km_shear)  &
                    - c4 * eps(k)**2 / e(k) &
                    + c5 * (  (Km(k+1)+Km(k)) * (eps(k+1)-eps(k)) / 2 / (zw(k+1)-zw(k))     &
                            - (Km(k)+Km(k-1)) * (eps(k)-eps(k-1)) / 2 / (zw(k)-zw(k-1))   ) &
                        / (z(k+1) - z(k))           
      
      out%e(k) =   Km_shear - Kh_Pbuoy &
                 + c1 * (  (Km(k+1)+Km(k)) * (e(k+1)-e(k)) / 2 / (zw(k+1)-zw(k))     &
                         - (Km(k)+Km(k-1)) * (e(k)-e(k-1)) / 2 / (zw(k)-zw(k-1))   ) &
                       / (z(k+1) - z(k)) &
                 - eps(k)
    end do
    
    end associate
  
  end subroutine
  
  subroutine compute_dt_flow(out, in)
    type(fields), intent(inout) :: out
    type(fields), intent(in) :: in
    integer :: k
    
    do k = 1, nz
      out%u(k) =   f_Coriolis * (in%v(k) - v_g(k)) - (in%up_wp(k) - in%up_wp(k-1)) / (zw(k) - zw(k-1))
      out%v(k) = - f_Coriolis * (in%u(k) - u_g(k)) - (in%vp_wp(k) - in%vp_wp(k-1)) / (zw(k) - zw(k-1))
      out%theta(k) =                       - (in%thetap_wp(k) - in%thetap_wp(k-1)) / (zw(k) - zw(k-1))
    end do
  end subroutine
end module

module saving
  use kinds
  use grid_variables
  
  implicit none
contains

  subroutine save_profiles(u, time)
    type(fields) :: u
    real(rp) :: time
    integer :: iu, k
    
    open(newunit=iu, file="vars_p.txt")
    write(iu,'(*(g0,1x))') "# k", "z", "u", "v" ,"theta"
    do k = 1, nz
      write(iu,*) k, z(k), u%u(k), u%v(k), u%theta(k)
    end do
    close(iu)

    open(newunit=iu, file="vars_w.txt")
    write(iu,'(*(g0,1x))') "# k", "zw", "e", "eps", "u'w'", "theta'w'", "Km", "Kh"
    do k = 1, nzw
      write(iu,*) k, zw(k), u%e(k), u%eps(k), u%up_wp(k), u%thetap_wp(k), u%Km(k), u%Kh(k)
    end do
    close(iu)
  end subroutine

end module

include "save_nc.f90"

program scm
  use kinds
  use global_variables
  use boundary_conditions
  use initial_conditions
  use boundary_conditions
  use turbulence_closure
  use prognostic_equations
  use saving
  use save_nc
  
  implicit none
  
  type(fields) :: u, du, u2, du2
  
  real(rp) :: time, delta_t, start_time, end_time, save_period
  
  integer :: itime, iter, max_it = 1e5, save_period_time_step
  
  start_time = 0._rp
  end_time = 24*3600
  delta_t = 0.5_rp
  
  save_period = 60._rp
  save_period_time_step = nint(save_period / delta_t)
  
  time = start_time

  
  max_it = nint((end_time-start_time) / delta_t)
  
  call init_grid
  
  call u%allocate
  call u2%allocate
  call du%allocate
  call du2%allocate
  
  call init_geostrophic_windspeed
  
  call init_boundary_conditions(theta_flux_fixed_arg=.false., theta_s=280._rp, thetap_wp_s=0._rp)
  
  call set_initial_conditions(u)
  
  call save_nc_init
  
  call save_nc_profile(u, time)
  
  do itime = 1, max_it
    call set_boundary_conditions(u, theta_s = theta_func(time))
    
    call compute_turbulence_K(u)
    call compute_dt_turbulence(du, u)
    call compute_dt_flow(du, u)
    
    u2 = u + du * delta_t
        
    do iter = 1, 4
      call set_boundary_conditions(u2, theta_s = theta_func(time+delta_t))
      
      call compute_turbulence_K(u2)
      call compute_dt_turbulence(du2, u2)
      call compute_dt_flow(du2, u2)

      u2 = u + (du + du2) * (delta_t/2)
    end do
    
    time = time + delta_t
    write(*,'(2x,f0.2)') time
    
    u = u2

    if (mod(itime, save_period_time_step)==0) call save_nc_profile(u, time)
  end do
  
   
!   call save_profiles(u, time)
  call save_nc_finalize
  
contains

  function theta_flux_func(time) result(res)
    real(rp) :: res
    real(rp), intent(in) :: time
    real(rp), parameter :: period = 24 * 3600
    res = 0.1 * sin(2 * pi * time / period)
    
    if (res<0) res = res/40
  end function

  function theta_func(time) result(res)
    real(rp) :: res
    real(rp), intent(in) :: time
    real(rp), parameter :: period = 24 * 3600
    res = 290 + 5 * sin(2 * pi * time / period)
    
  end function

end program scm

