!======================
!Created on 10/21/2020
!@author: Jun Zhou @UMD
!======================

module ant_prj_mod
  implicit none
  real(8),allocatable :: lat(:,:),lon(:,:)

  contains

    subroutine ant_prj(ntheta,nphi,dtheta,dphi,lat_ct,lon_ct,satzen,satazi,satran)
      use axis
      use math
      use footprint
      integer,intent(in)::nphi,ntheta
      real(8),intent(in)::dphi,dtheta,lat_ct,lon_ct,satzen,satazi,satran
      real(8) :: los_unit_ecr(3),satpos_ecr(3),ortho_unit_ecr(3)
      real(8) :: fovVector(3),curFovVector(3),geolat,geolon,theta,phi    
      integer::iphi,itheta


      if (allocated(lat)) deallocate(lat)
      if (allocated(lon)) deallocate(lon)
      allocate(lat(ntheta,nphi),lon(ntheta,nphi))


      call determine_los_satpos_ecr(lat_ct,lon_ct,satran,satazi,satzen, &
                                    los_unit_ecr,satpos_ecr,ortho_unit_ecr)

      do itheta=1,ntheta
        theta=dble(itheta-1)*dtheta
        fovVector = rotate_vector(ortho_unit_ecr,theta,los_unit_ecr)
        do iphi=1,nphi        
          phi=dble(iphi-1)*dphi
          curFovVector = rotate_vector(los_unit_ecr,phi,fovVector)
          call compute_geolocation(satpos_ecr, curFovVector,EQUAT_RAD,FLATFAC,&
                                   geolat,geolon)
          lat(itheta,iphi)=geolat 
          lon(itheta,iphi)=geolon
        enddo
      enddo


      return
    end subroutine ant_prj
end module ant_prj_mod
