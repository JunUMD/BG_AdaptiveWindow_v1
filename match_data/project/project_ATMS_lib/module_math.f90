module math
  use constants_m
  implicit none

contains
!********function findgen************************
function findgen(n)
  integer::i,n
  real(8)::findgen(n)

  do  i=1,n
    findgen(i)=float(i-1)
  end do
  return
end function findgen
!*************************************************
!**********function cross*************************
function cross(a, b)
  real(8) :: cross(3)
  real(8) :: a(3), b(3)

  cross(1) = a(2) * b(3) - a(3) * b(2)
  cross(2) = a(3) * b(1) - a(1) * b(3)
  cross(3) = a(1) * b(2) - a(2) * b(1)
  return
end function cross
!**********function cross2d*****************************
function cross2d(a, b)
  real(8) :: cross2d
  real(8) :: a(2), b(2)
  cross2d = a(1) * b(2) - a(2) * b(1)
  return
end function cross2d
!*********function norm2*******************************
function norm2(x)
  implicit none
  intrinsic :: dot_product,sqrt
  real(8),intent(in)::x(:)
  real(8) :: norm2
  norm2=sqrt(dot_product(x,x))
end function
!*******************************************************
function IfInCircle(rX,rY,rX_Circle,rY_Circle,n)
  real(8)::rX,rY,rDiff1(2),rDiff2(2),rTheta,rCosTheta,rSinTheta
  integer(4)::i,n
  real(8)::rX_Circle(n),rY_Circle(n),rNorm1,rNorm2
  logical:: IfInCircle

  rTheta=0.
  IfInCircle=.False.
  do i=1,n-1
    rDiff1(1)=rX_Circle(i)-rX
    rDiff1(2)=rY_Circle(i)-rY
    rDiff2(1)=rX_Circle(i+1)-rX
    rDiff2(2)=rY_Circle(i+1)-rY
    rNorm1=sqrt(rDiff1(1)**2+rDiff1(2)**2)
    rNorm2=sqrt(rDiff2(1)**2+rDiff2(2)**2)
    if ((rNorm1<zero) .or. (rNorm2<zero)) then
      rTheta=2*pi
      exit
     else
       rCosTheta=dot_product(rDiff1,rDiff2)/(rNorm1*rNorm2)
       rSinTheta=cross2d(rDiff1,rDiff2)/(rNorm1*rNorm2)
       if (rSinTheta<zero) then
         rTheta=rTheta-acos(rCosTheta)
       else
         rTheta=rTheta+acos(rCosTheta)
       end if
     end if
  end do
  if (abs(rTheta)<zero) then   !if rTheta=0, outside
    IfInCircle=.False. 
  else if (abs(rTheta)>6.0) then !if rTheta=2Pi,inside 
    IfInCircle=.True.
  end if
end function IfInCircle
!*************************************************

!!Fit polynomial:
!!y=a(1)+a(2)*(x-xmean)+a(3)*(x-xmean)^2+a(4)*(x-xmean)^3 ...
!!# of data: n
!!# of terms: m
!!data points: x(n),y(n)
!!fitting coef: a(m)
!!fitting error: dt1,dt2,dt3
!!This subroutine is from <<Fortran common use algorithm set>> 
!!version 2, Xu Shiliang, p341

SUBROUTINE HPIR1(X,Y,A,N,M,DT1,DT2,DT3)
implicit none
integer :: I,J,K,M,N
DIMENSION X(N),Y(N),A(M),S(20),T(20),B(20)
DOUBLE PRECISION X,Y,A,S,T,B,DT1,DT2,DT3,&
                 Z,D1,P,C,D2,G,Q,DT
DO 5 I=1,M
5 A(I)=0.0
IF (M.GT.N) M=N
IF (M.GT.20) M=20
Z=0.0
DO 10 I=1,N
10 Z=Z+X(I)/N
B(1)=1.0
D1=N
P=0.0
C=0.0
DO 20 I=1,N
  P=P+(X(I)-Z)
  C=C+Y(I)
20 CONTINUE
C=C/D1
P=P/D1
A(1)=C*B(1)
IF (M.GT.1) THEN
  T(2)=1.0
  T(1)=-P
  D2=0.0
  C=0.0
  G=0.0
  DO 30 I=1,N
    Q=X(I)-Z-P
    D2=D2+Q*Q
    C=Y(I)*Q+C
    G=(X(I)-Z)*Q*Q+G
30  CONTINUE

  C=C/D2
  P=G/D2
  Q=D2/D1
  D1=D2
  A(2)=C*T(2)
  A(1)=C*T(1)+A(1)
END IF
DO 100 J=3,M
  S(J)=T(J-1)
  S(J-1)=-P*T(J-1)+T(J-2)
  IF (J.GE.4) THEN
    DO 40 K=J-2,2,-1
40    S(K)=-P*T(K)+T(K-1)-Q*B(K)
  END IF
  S(1)=-P*T(1)-Q*B(1)
  D2=0.0
  C=0.0
  G=0.0
  DO 70 I=1,N
    Q=S(J)
    DO 60 K=J-1,1,-1
60    Q=Q*(X(I)-Z)+S(K)
    D2=D2+Q*Q
    C=Y(I)*Q+C
    G=(X(I)-Z)*Q*Q+G
70  CONTINUE
  C=C/D2
  P=G/D2
  Q=D2/D1
  D1=D2
  A(J)=C*S(J)
  T(J)=S(J)
  DO 80 K=J-1,1,-1
    A(K)=C*S(K)+A(K)
    B(K)=T(K)
    T(K)=S(K)
80  CONTINUE
100 CONTINUE
DT1=0.0
DT2=0.0
DT3=0.0
DO 120 I=1,N
  Q=A(M)
  DO 110 K=M-1,1,-1
110  Q=Q*(X(I)-Z)+A(K)
  DT=Q-Y(I)
  IF (ABS(DT).GT.DT3) DT3=ABS(DT)
  DT1=DT1+DT*DT
  DT2=DT2+ABS(DT)
120 CONTINUE
RETURN
END SUBROUTINE HPIR1


!===========================================================
!interpolate n variables of a point p between p1 
!and p2
!factor= distance from p  to p1
!       _________________________
!        distance from p2 to p1
!longitude should be first transformed from 0-+-pi to 0-2pi,
!then changed back after the interpolation
!author: Jun Zhou 5/19/2015
!============================================================
subroutine interpolate_1d(p1,p2,factor,p,longitude)

  real(8),intent(in) :: p1,p2,factor
  logical,intent(in),optional :: longitude
  real(8),intent(out) :: p
  
  real(8) :: lon(2)

  if(present(longitude))then
     lon(1)=p1
     lon(2)=p2
     if(lon(1)<zero) lon(1)=lon(1)+pi*2d0
     if(lon(2)<zero) lon(2)=lon(2)+pi*2d0
     p=lon(1)*(1d0-factor)+lon(2)*factor
     if(p>pi) p=p-pi*2d0
  else
     p=p1*(1d0-factor)+p2*factor
  endif

return
end subroutine interpolate_1d

subroutine interpolate_nd(n,p1,p2,factor,p)

  integer(4),intent(in) :: n
  real(8),intent(in) :: p1(n),p2(n),factor
  real(8),intent(out) :: p(n)
  
  integer(4) :: i

  do i=1,n
     p(i)=p1(i)*(1d0-factor)+p2(i)*factor
  enddo

return
end subroutine interpolate_nd

end module
