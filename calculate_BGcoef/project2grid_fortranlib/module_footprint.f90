
module footprint
   use axis
   use math
   use constants_m
   implicit none 
   
   public :: determine_los_satpos_ecr

contains

  subroutine determine_los_satpos_ecr(geodeticLat, longitude,&
                                      range1, azimuth, zenith,&
                                      los_unit_ecr,satpos_ecr,ortho_unit_ecr)
    real(8)::geodeticLat, longitude,range1, azimuth, zenith,&
             los_unit_ecr(3),satpos_ecr(3),ortho_unit_ecr(3)
    real(8)::lonRad,azimuthRad,zenithRad,geodeticLatRad
    real(8)::geocentricLatRad,geocentricRange,geocentricTempX
    real(8)::geocentricPoint(3),losPointX(3),losPointX1(3),losPointX2(3),&
             losPointX3(3),losPointX4(3),rotationVector(3)
    real(8)::geocentricLat,useLonRad,orthoVectorLOS(3),&
             satellitePosition(3),rtemp(3)

    
    lonRad = longitude  !* deg2rad
    azimuthRad = azimuth !* deg2rad
    zenithRad = zenith   !* deg2rad
    geodeticLatRad = geodeticLat !* deg2rad

    !  Find the geocentric latitude.
    !  1 = geodetic to geocentric, 2 = geocentric to geodetic.
    geocentricLatRad = latitude_transform( geodeticLatRad, FLATFAC, 1)
    geocentricLat = geocentricLatRad * rad2deg

    geocentricRange  = geocentric_range(geocentricLat, EQUAT_RAD, FLATFAC)

    geocentricTempX = geocentricRange * cos(geocentricLatRad)
    geocentricPoint(3) = geocentricRange * sin(geocentricLatRad)
    geocentricPoint(1) = geocentricRange * cos(geocentricLatRad) * cos( lonRad)
    geocentricPoint(2) = geocentricRange * cos(geocentricLatRad) * sin( lonRad)


    !  Find the satellite line-of-sight ( look vector) assuming
    !  that the geolocation point is (1, 0, 0). The next step will be
    !  to perform 2 rotation to find the LOS in the geocentric frame
    !  or Earth Centered Reference (ECR).
    !  Here, the X axis is outward at longitude zero, Z is toward the North pole
    !  and Y is orthogonal to X and Z.
    losPointX = (/ -1.0e0, 0.0e0, 0.0e0/)

    rotationVector = (/ 0.0e0, -1.0e0, 0.0e0/)
    losPointX1 = rotate_vector( rotationVector, zenithRad, losPointX);
    rotationVector = (/-1.0e0, 0.0e0, 0.0e0/)
    losPointX2 = rotate_vector( rotationVector, azimuthRad, losPointX1);
    !  Now find the LOS in the geocentric frame. This is done by
    !  performing 2 rotations: 
    !   1) geodetic latitude around [0,-1,0],
    !   2) longitude around [ 0, 0, 1].

    rotationVector = (/ 0.0e0, -1.0e0, 0.0e0/)
    losPointX3 = rotate_vector( rotationVector, geodeticLatRad, losPointX2)
    rotationVector = (/ 0.0e0, 0.0e0, 1.0e0/)

    if (lonRad .LE.  0.0) then
       useLonRad = lonRad + 2.0*pi
    else
       useLonRad = lonRad
    end if
    losPointX4 = rotate_vector( rotationVector, useLonRad, losPointX3)


    !------------------
    satellitePosition = geocentricPoint - (losPointX4 * range1)
    rtemp=(/0.0e0, 0.0e0, 1.0e0/)
    orthoVectorLOS = cross(losPointX4,rtemp)
    orthoVectorLOS = orthoVectorLOS /(sqrt(dot_product(orthoVectorLOS,orthoVectorLOS)))

    !--- output ---- 
    los_unit_ecr=losPointX4
    satpos_ecr=satellitePosition
    ortho_unit_ecr=orthoVectorLOS

    return 
  end subroutine determine_los_satpos_ecr



end module
