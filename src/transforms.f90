module transform
  implicit none

contains
  function Rinv(u,a)
    real,dimension(5,5) :: Rinv
    real :: u,a
    real :: gamma = 1.4

  ! Rinv(1,1) = ((gamma - 1)/4)*u*u/(a*a) + u/(2*a)
  ! Rinv(2,1) = 1 - ((gamma - 1)/2) * u*u/(a*a)
  ! Rinv(3,1) = ((gamma - 1)/4)*u*u/(a*a) - u/(2*a)
  !
  ! Rinv(1,2) = -1 * ( ((gamma - 1)/2)*u/(a*a) + 1/(2*a))
  ! Rinv(2,2) = (gamma - 1)*u/(a*a)
  ! Rinv(3,2) = -1 * ( ((gamma - 1)/2)*u/(a*a) - 1/(2*a))
  !
  ! Rinv(1,3) = ((gamma - 1)/(2*a*a))
  ! Rinv(2,3) = -1*((gamma - 1)/(a*a))
  ! Rinv(3,3) = ((gamma - 1)/(2*a*a))

    Rinv(1,:) = (/1,0,0,0,0/) !
    Rinv(2,:) = (/0,1,0,0,0/) !
    Rinv(3,:) = (/0,0,1,0,0/) !
    Rinv(4,:) = (/0,0,0,1,0/) !
    Rinv(5,:) = (/0,0,0,0,1/) !

  end function Rinv


  function Rcalc(u,a)
    real,dimension(5,5) :: Rcalc
    real :: gamma = 1.4
    real :: u,a

    ! Rcalc(1,1) = 1
    ! Rcalc(2,1) = u-a
    ! Rcalc(3,1) = a*a/(gamma - 1) + 0.5*u*u - u*a
    !
    ! Rcalc(1,2) = 1
    ! Rcalc(2,2) = u
    ! Rcalc(3,2) = 0.5*u*u
    !
    ! Rcalc(1,3) = 1
    ! Rcalc(2,3) = u+a
    ! Rcalc(3,3) = a*a/(gamma - 1) + 0.5*u*u + u*a

    Rcalc(1,:) = (/1,0,0,0,0/) !
    Rcalc(2,:) = (/0,1,0,0,0/) !
    Rcalc(3,:) = (/0,0,1,0,0/) !
    Rcalc(4,:) = (/0,0,0,1,0/) !
    Rcalc(5,:) = (/0,0,0,0,1/) !

  end function Rcalc
end module transform
