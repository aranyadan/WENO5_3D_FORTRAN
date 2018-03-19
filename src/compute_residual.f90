subroutine compute_residual(q,qo,residual,n_x,n_y,n_z)
  implicit none
  real, dimension(n_x,n_y,n_z,5) :: q,qo
  real, dimension(n_x,n_y,n_z) :: u,v,E,p,a,u2,rho,w
  real, dimension(5) :: res
  real :: residual
  integer :: dim,n_x,n_y,n_z

  call primitives(q,n_x,n_y,n_z,rho,u,v,w,E,p,a)
  call primitives(qo,n_x,n_y,n_z,rho,u2,v,w,E,p,a)
  do dim=1,1
    res(dim) = SUM(ABS(u-u2))/(MAX(1,size(u)));
  end do
  residual = res(1);


end subroutine compute_residual
