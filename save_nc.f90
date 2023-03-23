module save_nc
  use kinds
  use grid_variables
  use netcdf
  
  implicit none 
  
  private
  
  public :: save_nc_init, save_nc_finalize, save_nc_profile
  
  integer :: id, dims(3)
  integer :: dim_time, dim_z, dim_zw
  integer :: vars(10)
  integer :: n_profile = 0
  
  character, parameter :: nl = new_line(' ')

contains

  subroutine save_nc_init
    integer :: err, k, idx
    
    idx = 1
    
    err = nf90_create('scm_result.nc', NF90_CLOBBER, id)

    err = nf90_def_dim(id, 'time', NF90_UNLIMITED, dims(1))

    err = nf90_def_dim(id, 'z', nz, dims(2))

    err = nf90_def_dim(id, 'zw', nzw+1, dims(3))

    err = nf90_def_var(id, 'time', NF90_REAL, dims(1:1), dim_time)
    
    err = nf90_def_var(id, 'z', NF90_REAL, dims(2:2), dim_z)
    
    err = nf90_def_var(id, 'zw', NF90_REAL, dims(3:3), dim_zw)
    
    err = nf90_put_att(id, dim_time, 'units', 's')
    err = nf90_put_att(id, dim_time, 'long_name', 'time')

    err = nf90_put_att(id, dim_time, '_FillValue', NF90_FILL_REAL)
    err = nf90_put_att(id, dim_z, '_FillValue', NF90_FILL_REAL)
    err = nf90_put_att(id, dim_zw, '_FillValue', NF90_FILL_REAL)

    err = nf90_put_att(id, NF90_GLOBAL, '_FillValue', NF90_FILL_REAL)
    err = nf90_put_att(id, NF90_GLOBAL, 'Model', 'Educational Single Column Model'//nl// &
                                                'https://')

    err = nf90_put_att(id, NF90_GLOBAL, 'Turbulence scheme', &
                                        'e-epsilon'//nl// &
                                        'following Zhang, Wang, Xue (2020) https://10.1175/MWR-D-19-0084.1')
    call def_z('u', 'm/s')
    call def_z('v', 'm/s')
    call def_z('theta', 'K')
    
    call def_zw('e', 'm^2/s^2')
    call def_zw('eps', 'm^2/s^3')
    call def_zw("u'w'", 'm^2/s^2')
    call def_zw("v'w'", 'm^2/s^2')
    call def_zw("theta'w'", 'm^2/s^2')
    call def_zw("Km", 'm^2/s')
    call def_zw("Kh", 'm^2/s')
    
    err = nf90_enddef(id)
   
    err = nf90_put_var(id, dim_z, [(z(k), k=1, nz)] )
    err = nf90_put_var(id, dim_zw, [(zw(k), k=0, nzw)] )
    
  contains
    subroutine def_z(name, units, long_name)
      character(*) :: name, units
      character(*), optional :: long_name
      call def(name, units, [2,1], long_name)
    end subroutine

    subroutine def_zw(name, units, long_name)
      character(*) :: name, units
      character(*), optional :: long_name
      call def(name, units, [3,1], long_name)
    end subroutine
    
    subroutine def(name, units, ds, long_name)
      character(*) :: name, units
      integer :: ds(:)
      character(*), optional :: long_name
      
      
      err = nf90_def_var(id, name, NF90_REAL, dims(ds), vars(idx))
    
      err = nf90_put_att(id, vars(idx), 'units', units)
      err = nf90_put_att(id, vars(idx), '_FillValue', NF90_FILL_REAL)
      if (present(long_name)) then
        err = nf90_put_att(id, vars(idx), 'long_name', long_name)
      else
        err = nf90_put_att(id, vars(idx), 'long_name', name)
      end if
      idx = idx + 1
    end subroutine
    
  end subroutine
  
  subroutine save_nc_finalize
    integer :: err
  
    err = nf90_close(id)
  end subroutine
  
  
  subroutine save_nc_profile(u, time)
    type(fields), intent(in) :: u
    real(rp), intent(in) :: time
    integer :: err, idx
    n_profile = n_profile + 1

    err = nf90_put_var(id, dim_time, [time], start=[n_profile], count=[1])
    idx = 1    
    call put_z(u%u)
    call put_z(u%v)
    call put_z(u%theta)
    call put_zw(u%e)
    call put_zw(u%eps)
    call put_zw(u%up_wp)
    call put_zw(u%vp_wp)
    call put_zw(u%thetap_wp)
    call put_zw(u%Km)
    call put_zw(u%Kh)
  contains
    subroutine put_z(arr)
      real(rp) :: arr(:)
      err = nf90_put_var(id, vars(idx), arr, start=[1,n_profile], count=[nz,1])
      idx = idx + 1
    end subroutine
    subroutine put_zw(arr)
      real(rp) :: arr(:)
      err = nf90_put_var(id, vars(idx), arr, start=[1,n_profile], count=[nzw+1,1])
      idx = idx + 1
    end subroutine
  end subroutine
end module
