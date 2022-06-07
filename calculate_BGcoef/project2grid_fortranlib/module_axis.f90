
module axis
   use constants_m
   implicit none
   public::determine_deltathetaphi_prime
contains

  subroutine determine_deltathetaphi_prime(beamwidth, &
                                   satzen_fovcenter,satalt_fovcenter,&
                                   ntheta_prime,nphi_prime,delta_theta_prime,delta_phi_prime)
    implicit none
    real(kind=8),intent(in)::beamwidth,satzen_fovcenter,satalt_fovcenter
    real(kind=8),intent(out)::delta_theta_prime,delta_phi_prime
    integer(4),intent(out)::ntheta_prime,nphi_prime
    real(kind=8)::satzen,satalt,radiom_half_beamwidth_rad,theta_prime_max,delta_theta_prime_max,&
                  delta_phi_prime_max,delta_theta_prime_equiv,delta_phi_prime_equiv


    radiom_half_beamwidth_rad=beamwidth/2d0
    theta_prime_max=2.5d0*radiom_half_beamwidth_rad
    delta_theta_prime_max=0.1*deg2rad !theta_prime_max/100d0
    delta_phi_prime_max=delta_theta_prime_max/sin(theta_prime_max)  


    satzen=satzen_fovcenter+theta_prime_max
    
    satalt=satalt_fovcenter*cos(satzen_fovcenter)/cos(satzen)

    delta_theta_prime_equiv=(delx)*1000d0*cos(satzen)/satalt

    delta_phi_prime_equiv=(delx)*1000d0*cos(satzen)/(satalt*sin(theta_prime_max))


    delta_theta_prime         = delta_theta_prime_equiv
    delta_phi_prime           = delta_phi_prime_equiv
    if(delta_theta_prime .gt. delta_theta_prime_max) &
       delta_theta_prime       = delta_theta_prime_max
    if(delta_phi_prime .gt. delta_phi_prime_max) &
       delta_phi_prime         = delta_phi_prime_max
    ntheta_prime              = floor(theta_prime_max/delta_theta_prime)+1
    nphi_prime                = floor(2.*pi/delta_phi_prime)+1
    delta_theta_prime         = theta_prime_max/float(ntheta_prime-1)
    delta_phi_prime           = 2.*pi/float(nphi_prime-1)

    return
  end subroutine determine_deltaThetaPhi_prime


   !calculate the distance between two points on the spherical surface
   !This function is from ATMS-VIIRS-COLLO-V1/Module-commfunc.f90
   function arc_length(lat1,lon1,lat2,lon2)
      implicit none
      real(8) ::  arc_length
      real(8),intent(in) :: lat1,lat2,lon1,lon2 !geocentric radian
      real(8) :: BB

      BB=SIN(lat1)*SIN(lat2)+COS(lat1)*COS(lat2)*COS(lon1-lon2)
      IF(BB.GT.1.0d0) BB=1.0d0
      IF(BB.LT.-1.0d0) BB=-1.0d0
      BB=ACOS(BB)
      arc_length=EQUAT_RAD*BB/1000.0d0  !(km)

    return

  end function

  subroutine calculate_grid_index(lat_gedrad,lon_rad,grid_origin_lat_cenrad,grid_origin_lon_rad,&
                                  numx,numy,nx_foot,ny_foot) 
    implicit none
    real(8),intent(in)::lat_gedrad,lon_rad
    real(kind=4),intent(in)::grid_origin_lat_cenrad,grid_origin_lon_rad
    integer(kind=2),intent(in)::numx,numy
    integer(kind=2),intent(out)::nx_foot,ny_foot
    real(8)::lat_cenrad,dx,dy

    lat_cenrad=latitude_transform(lat_gedrad,FLATFAC,1)
    dx=arc_length(lat_cenrad,dble(grid_origin_lon_rad),lat_cenrad,lon_rad)
    dy=arc_length(dble(grid_origin_lat_cenrad),lon_rad,lat_cenrad,lon_rad)
    nx_foot=nint(dx/delx)+1
    ny_foot=nint(dy/dely)+1

    if(nx_foot>numx .or. ny_foot>numy .or. nx_foot<=0 .or. ny_foot<=0)then
      print *,'calculate_grid_index: nx_foot,ny_foot are invalid!'
      print *,lon_rad,grid_origin_lon_rad,dx,delx
      print *,lat_cenrad,grid_origin_lat_cenrad,dy,dely
      print *,'numx=',numx
      print *,'numy=',numy
      print *,'nx_foot=',nx_foot
      print *,'ny_foot=',ny_foot
      stop
    endif

    return
  end subroutine calculate_grid_index


!! Description:
!!
!!*********   Given a Earth equatorial radius, the geocentric latitude, and
!!*********   the flattening factor, find the geocentric range.
!!
!! Usage:
!!
!!   [range] = geocentric_range( latitude, equatorialRadius, flatteningFactor )
!!
!! Input:
!!
!! latitude    : Geocentric latitude  of location 1, in degrees (geographical coordinates).
!!     equatorialRadius
!!    : Earth equaorial radius in meters ( WGS84 is 6378137.0 meters).
!!     flatteningFactor
!!     : Flattening factor of the Earth, WGS84 value is 1/298.257223560
!!
!! Output:
!!
!!     range       :  Geocetric distance at the given geocentric latitude.
!! 
!! Algorithm: 
!!     x**2/a**2+y**2/b**2=1  x=r*cos(lat) y=r*sin(lat)  =>
!!     r=ab/sqrt((a*sin(lat))**2+(b*cos(lat))**2)
!! 
!***********************function geocentric_range*********************************
 ! input is lat in degree 
 function geocentric_range(geocentricLatitude, equatorialRadius,flatteningFactor)

    real(8)::geocentricLatitude, equatorialRadius,flatteningFactor  
    real(8)::geocentric_range,polarRadius,lat
   
    polarRadius = (1.0e0 - flatteningFactor) * equatorialRadius
    lat = geocentricLatitude * deg2rad
    !geocentric_range = sqrt((equatorialRadius * cos(lat))**2 + (polarRadius * sin(lat))**2)
    geocentric_range = equatorialRadius*polarRadius/ &  !modify: correct the algorithm
                       sqrt((equatorialRadius*sin(lat))**2+(polarRadius*cos(lat))**2)

    return                      
  end function geocentric_range

!input is lat in radian:
  function geocentric_range_rad(geocentricLatitude, equatorialRadius,flatteningFactor)
    
    real(8)::geocentricLatitude, equatorialRadius,flatteningFactor  
    real(8)::polarRadius,lat,geocentric_range_rad
    
    polarRadius = (1.0e0 - flatteningFactor) * equatorialRadius
    lat = geocentricLatitude 
    !geocentric_range = sqrt((equatorialRadius * cos(lat))**2 + (polarRadius * sin(lat))**2)
    geocentric_range_rad = equatorialRadius*polarRadius/ &  !modify: correct the algorithm
                       sqrt((equatorialRadius*sin(lat))**2+(polarRadius*cos(lat))**2)

    return                      
   
  end function geocentric_range_rad
!*********************************************************************************************

!***********************function latitude_transform*********************************************
!input: latRef: geocentric lat or geodetic lat (rad)
!flatteningFactor: flattening (a-b)/a
!option: 1 geodetic => geocentric
!        2          <=
!latitude_transform: transformed result (rad)
  function latitude_transform(latRef, flatteningFactor, option)
    real(8)::latitude_transform,latRef,flatteningFactor
    integer::option
    real(8)latIn,ff,latout,latitudeOutput

    latIn = latRef
    ff = 1.0 - flatteningFactor
    if( option ==  1) then
      latOut = atan( ((ff)**2) * tan(latIn)) 
      latitudeOutput = latOut 
    else 
      latOut = atan((1.0/ (ff**2)) * tan(latIn)) 
      latitudeOutput = latOut 
    end if

    latitude_transform=latitudeOutput
    return
  end function latitude_transform
!****************************************************************************************

!**********************function compute_geolocation**************************************
  !  The basic equations are:
  !  P + lambda LOS = G  ( P = satellite position, LOS = line of sight, G =geolocation point).
  !  and lambda is the slant range.
  !  Gx^2 / a^2 + Gy^2 / a^2 + Gz^2 / c^2 = 1    Earth ellipsoid equation, a =equatorial radius.
  !  c is polar radius ( c = ( 1-f) * a) where f is flattening factor.
  
  subroutine compute_geolocation(satellitePosition,lineOfSight,earthRadius,flatFact,&
                                geoLat,geoLon)

    real(8),intent(in)::lineOfSight(3),satellitePosition(3),flatFact,earthRadius
    real(8),intent(out)::geoLat,geoLon
    real(8)::geolocationPoint(3),polarRadius,termA,termB,termC,radical
    real(8)::slantRange,slantRange1,slantRange2
    real(8)::magXY,geocentricLat
    
    polarRadius = earthRadius * ( 1.0 - flatFact)
    termA = ((lineOfSight(1) / earthRadius)**2) + &
           ((lineOfSight(2) / earthRadius)**2) + &
           ((lineOfSight(3) / polarRadius)**2)

    termB =  ( satellitePosition(1) * lineOfSight(1) / (earthRadius**2)) + &
            ( satellitePosition(2) * lineOfSight(2) / (earthRadius**2)) + &
            ( satellitePosition(3) * lineOfSight(3) / (polarRadius**2))  
    termB = termB * 2.0
 
    termC =  (satellitePosition(1)/earthRadius)**2 + &
             (satellitePosition(2)/earthRadius)**2 + &
             (satellitePosition(3)/polarRadius)**2 - 1.0

   radical = termB**2 - (4.0 * termA * termC)
   !print *,'radical=',radical
   if( radical .lt. 0.0) then 
      !  The line of sight does not intercept the Earth ellipsoid.
      geolat=rmissing !jun10/5
      geolon=rmissing !jun10/5
      return  !jun10/5
   endif 
   if( radical .eq. 0.0) then 
      !  The line of sight  intercepts the Earth ellipsoid tangentially.
      slantRange = -termB / (2.0 * termA)
      geolocationPoint = satellitePosition + slantRange * lineOfSight
   endif 
   if( radical .gt.  0.0) then 
      !  The line of sight intercepts the Earth ellipsoid at 2 point, the solution
      !  is the shorter slant range.
      slantRange1 = (-termB - sqrt(radical))  / (2.0 * termA)
      slantRange2 = (-termB + sqrt(radical))  / (2.0 * termA)
      slantRange = min(slantRange1, slantRange2)
      geolocationPoint = satellitePosition + slantRange * lineOfSight
   endif
   geoLon = atan2(geolocationPoint(2),geolocationPoint(1))
   magXY = sqrt( geolocationPoint(1)**2 +  geolocationPoint(2)**2) 
   geocentricLat = atan2( geolocationPoint(3), magXY)
   geoLat = latitude_transform( geocentricLat, flatFact, 2)
   return
  
   end subroutine compute_geolocation
!*****************************************************************************

!****************function rotate_vector***************************************
function rotate_vector(rotationAxis, angle, oldVector)
  use math
  real(8)::rotate_vector(3),rotationAxis(3),oldVector(3)
  real(8)::angle,cosAngle,sinAngle
  real(8)::firstTerm(3),secondTerm(3),thirdTerm(3)
  
  cosAngle = cos(angle)        
  sinAngle = sin(angle)        

  firstTerm = oldVector * cosAngle 
  secondTerm =  cross( oldVector, rotationAxis) * sinAngle 
  thirdTerm = (dot_product(oldVector, rotationAxis) * rotationAxis) * (1.0-cosAngle)
  rotate_vector = firstTerm - secondTerm + thirdTerm 
  return 
end function rotate_vector
end module
