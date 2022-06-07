module ant_prj_mod
  implicit none
  real(8),allocatable :: lat(:),lon(:)

  contains

    subroutine ant_prj(nphi,dphi,beamwidth,lat_ct,lon_ct,satzen,satazi,satran)
      use axis
      use math
      use footprint
      integer,intent(in)::nphi
      real(8),intent(in)::dphi,beamwidth,lat_ct,lon_ct,satzen,satazi,satran
      real(8) :: los_unit_ecr(3),satpos_ecr(3),ortho_unit_ecr(3),phi
      real(8) :: fovVector(3),curFovVector(3),geolat,geolon        
      integer::iphi 

      if (allocated(lat)) deallocate(lat)
      if (allocated(lon)) deallocate(lon)
      allocate(lat(nphi),lon(nphi))


      call determine_los_satpos_ecr(lat_ct,lon_ct,satran,satazi,satzen, &
                                    beamwidth, &
                                    los_unit_ecr,satpos_ecr,ortho_unit_ecr)

      fovVector = rotate_vector(ortho_unit_ecr,beamwidth/2.0,los_unit_ecr)
      do iphi=1,nphi        
        phi=dble(iphi-1)*dphi
        curFovVector = rotate_vector(los_unit_ecr,phi,fovVector)
        call compute_geolocation(satpos_ecr, curFovVector,EQUAT_RAD,FLATFAC,&
                                 geolat,geolon)
        lat(iphi)=geolat !geodetic radians
        lon(iphi)=geolon !radians
      enddo


      return
    end subroutine ant_prj
end module ant_prj_mod
