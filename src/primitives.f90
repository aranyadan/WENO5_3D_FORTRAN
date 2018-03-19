subroutine primitives(q,n_x,n_y,n_z,rho,u,v,w,E,p,a)
  implicit none
  integer :: n_x,n_y,n_z
  real,dimension(n_x,n_y,n_z,5) :: q
  real :: gamma = 1.4
  real,dimension(n_x,n_y,n_z) :: rho,u,v,w,E,p,a

  rho = q(:,:,:,1)
  u = q(:,:,:,2)/rho
  v = q(:,:,:,3)/rho
  w = q(:,:,:,4)/rho
  E = q(:,:,:,5)/rho
  p = (gamma-1.0)*rho*(E-0.5*(u*u+v*v+w*w))
  a = SQRT(gamma*p/rho)

end subroutine primitives
