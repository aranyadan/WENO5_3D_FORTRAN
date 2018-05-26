subroutine build_flux(q,n_x,n_y,n_z,delx,dely,delz,Re,Pr,Suth,F,G,H,VISCOUS)
  implicit none
  integer :: n_x,n_y,n_z,delx,dely,delz,VISCOUS
  real :: gamma=1.4,Re,Suth,Pr
  real, dimension(n_x,n_y,n_z,5) :: q,F, G,H
  real, dimension(n_x,n_y,n_z) :: u,v,w,p,rho,E,a,tauxx,tauxy,tauxz,tauyy,tauyz,&
                                  tauzz,q_x,q_y,q_z,T
  real, dimension(n_x,n_y,n_z,4) :: vel

  call primitives(q,n_x,n_y,n_z,rho,u,v,w,E,p,a)
  T = p/rho

  if(VISCOUS==1) then
    vel(:,:,:,1) = u(:,:,:)
    vel(:,:,:,2) = v(:,:,:)
    vel(:,:,:,3) = w(:,:,:)
    vel(:,:,:,4) = T(:,:,:)
    call viscous_fluxes(vel,p,rho,n_x,n_y,n_z,delx,dely,delz,Suth,Re,Pr,tauxx,&
                        tauxy,tauxz,tauyy,tauyz,tauzz,q_x,q_y,q_z)
  else
    tauxx=0.0
    tauyy=0.0
    tauzz=0.0
    q_x=0.0
    q_y=0.0
    q_z=0.0
  endif

  F(:,:,:,1) = rho*u
  F(:,:,:,2) = rho*u*u + p + (- tauxx)*VISCOUS
  F(:,:,:,3) = rho*u*v + (- tauxy)*VISCOUS
  F(:,:,:,4) = rho*u*w + (- tauxz)*VISCOUS
  F(:,:,:,5) = rho*u*E + p*u + (- u*tauxx - v*tauxy - w*tauxz - q_x)*VISCOUS

  G(:,:,:,1) = rho*v
  G(:,:,:,2) = rho*v*u + (- tauxy)*VISCOUS
  G(:,:,:,3) = rho*v*v + p + (- tauyy)*VISCOUS
  G(:,:,:,4) = rho*v*w + (-tauyz)*VISCOUS
  G(:,:,:,5) = rho*v*E + p*v + (- u*tauxy - v*tauyy - w*tauyz - q_y)*VISCOUS

  H(:,:,:,1) = rho*w
  H(:,:,:,2) = rho*w*u + (- tauxz)*VISCOUS
  H(:,:,:,3) = rho*w*v + (- tauyz)*VISCOUS
  H(:,:,:,4) = rho*w*w + p + (- tauzz)*VISCOUS
  H(:,:,:,5) = rho*w*E + p*w + (- u*tauxz - v*tauyz - w*tauzz - q_z)*VISCOUS


  ! write(*,*)MAXVAL(MAXVAL(tauxx(:,:),1))
  ! write(*,*)MINVAL(MINVAL(rho(:,:),1))


end subroutine build_flux
