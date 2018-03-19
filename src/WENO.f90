subroutine WENO(lambda,F,q,R_i,Rinv_i,hp_i,hn_i)
  implicit none
  real,dimension(1,5) :: u,v,umm,um,up,upp,vmm,vm,vp,vpp,p0n,p1n,p2n,B0n,B1n
  real,dimension(1,5) :: B2n,alpha0n,alpha1n,alpha2n,alphasumn,w0n,w1n,w2n,hn
  real,dimension(1,5) :: p0p,p1p,p2p,B0p,B1p,B2p,alpha0p,alpha1p,alpha2p
  real,dimension(1,5) :: alphasump,w0p,w1p,w2p,hp,res,res1,res2,hp_i,hn_i
  real,dimension(6,5) :: F,q,u_i,v_i,g_i,w_i,temp
  real :: d0n,d0p,d1n,d1p,d2n,d2p,lambda,epsilon
  integer :: n_x,i
  real,dimension(5,5) :: R_i,Rinv_i

  ! do i=1,6
  !   w_i(i,1) = Rinv_i(1,1) * q(i,1) + Rinv_i(1,2) * q(i,2) + Rinv_i(1,3) * q(i,3)
  !   w_i(i,2) = Rinv_i(2,1) * q(i,1) + Rinv_i(2,2) * q(i,2) + Rinv_i(2,3) * q(i,3)
  !   w_i(i,3) = Rinv_i(3,1) * q(i,1) + Rinv_i(3,2) * q(i,2) + Rinv_i(3,3) * q(i,3)
  !
  !   g_i(i,1) = Rinv_i(1,1) * F(i,1) + Rinv_i(1,2) * F(i,2) + Rinv_i(1,3) * F(i,3)
  !   g_i(i,2) = Rinv_i(2,1) * F(i,1) + Rinv_i(2,2) * F(i,2) + Rinv_i(2,3) * F(i,3)
  !   g_i(i,3) = Rinv_i(3,1) * F(i,1) + Rinv_i(3,2) * F(i,2) + Rinv_i(3,3) * F(i,3)
  !
  ! end do
  w_i = q
  g_i = F


  v_i = 0.5 * (g_i + lambda*w_i)

  temp = 0.5 * (g_i - lambda*w_i)

  ! u_i = turnx(temp,6,-1)


  ! Right flux
  ! compute u_{i+1/2}^{+}
  vmm(1,:) = v_i(1,:)
  vm(1,:)  = v_i(2,:)
  v(1,:)   = v_i(3,:)
  vp(1,:)  = v_i(4,:)
  vpp(1,:) = v_i(5,:)

  ! Polynomials
  p0p = (2.0*vmm - 7.0*vm + 11.0*v)/6.0
  p1p = (-vm     + 5.0*v  + 2.0*vp)/6.0
  p2p = (2.0*v   + 5.0*vp - vpp   )/6.0

  ! Smoothness indicators
  B0p = 13.0/12.0*(vmm-2.0 * vm + v  )**2.0 + 1.0/4.0*(vmm - 4.0*vm + 3.0*v)**2.0
  B1p = 13.0/12.0*(vm -2.0 * v  + vp )**2.0 + 1.0/4.0*(vm  - vp        )**2.0
  B2p = 13.0/12.0*(v  -2.0 * vp + vpp)**2.0 + 1.0/4.0*(3.0*v - 4.0*vp + vpp)**2.0

  !constants
  d0p = 1.0/10
  d1p = 6.0/10
  d2p = 3.0/10
  epsilon = 1E-6

  !alpha weights
  alpha0p = d0p/(epsilon + B0p)**2.0
  alpha1p = d1p/(epsilon + B1p)**2.0
  alpha2p = d2p/(epsilon + B2p)**2.0
  alphasump = alpha0p + alpha1p + alpha2p

  ! ENO weights
  w0p = alpha0p/alphasump
  w1p = alpha1p/alphasump
  w2p = alpha2p/alphasump

  ! Flux at the boundary u_{i+1/2}^{+}
  hp = w0p*p0p + w1p*p1p + w2p*p2p


  ! Left Flux
  ! compute u_{i+1/2}^{-}
  umm(1,:) = temp(2,:)
  um(1,:)  = temp(3,:)
  u(1,:)   = temp(4,:)
  up(1,:)  = temp(5,:)
  upp(1,:) = temp(6,:)

  ! Polynomials
  p0n = ( -umm + 5.0*um + 2.0*u  )/6.0
  p1n = ( 2.0*um + 5.0*u  - up   )/6.0
  p2n = (11.0*u  - 7.0*up + 2.0*upp)/6.0

  ! Smoothness indicators
  B0n = 13.0/12.0*(umm-2.0*um+u  )**2.0 + 1.0/4.0*(umm - 4.0*um + 3.0*u)**2.0
  B1n = 13.0/12.0*(um -2.0*u +up )**2.0 + 1.0/4.0*(um  - up        )**2.0
  B2n = 13.0/12.0*(u  -2.0*up+upp)**2.0 + 1.0/4.0*(3.0*u - 4.0*up + upp)**2.0

  !constants
  d0n = 3.0/10
  d1n = 6.0/10
  d2n = 1.0/10

  ! Alpha Weights
  alpha0n = d0n/(epsilon + B0n)**2.0
  alpha1n = d1n/(epsilon + B1n)**2.0
  alpha2n = d2n/(epsilon + B2n)**2.0
  alphasumn = alpha0n + alpha1n + alpha2n

  ! ENO weights
  w0n = alpha0n/alphasumn
  w1n = alpha1n/alphasumn
  w2n = alpha2n/alphasumn

  ! Numerical Flux at boundary u_{i+1/2}^{-}
  hn = w0n*p0n + w1n*p1n + w2n*p2n


  ! hp_i(1,1) = R_i(1,1) * hp(1,1) + R_i(1,2) * hp(1,2) + R_i(1,3) * hp(1,3)
  ! hp_i(1,2) = R_i(2,1) * hp(1,1) + R_i(2,2) * hp(1,2) + R_i(2,3) * hp(1,3)
  ! hp_i(1,3) = R_i(3,1) * hp(1,1) + R_i(3,2) * hp(1,2) + R_i(3,3) * hp(1,3)
  !
  ! hn_i(1,1) = R_i(1,1) * hn(1,1) + R_i(1,2) * hn(1,2) + R_i(1,3) * hn(1,3)
  ! hn_i(1,2) = R_i(2,1) * hn(1,1) + R_i(2,2) * hn(1,2) + R_i(2,3) * hn(1,3)
  ! hn_i(1,3) = R_i(3,1) * hn(1,1) + R_i(3,2) * hn(1,2) + R_i(3,3) * hn(1,3)
  hp_i = hp
  hn_i = hn


  return
end subroutine WENO
