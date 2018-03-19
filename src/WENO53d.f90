subroutine WENO53d(lambda,F,q,n_x,n_y,n_z,hp,hn,dir)
  use transform
  implicit none
  integer :: n_x,i,j,dir,k,n_y,n_z
  real,dimension(n_x,n_y,n_z,5) :: F,q,hp,hn,Ftemp,qtemp
  real,dimension(-1:n_x+3,-1:n_y+3,-1:n_z+3,5) :: Fnew, qnew
  real :: lambda
  real,dimension(1) :: rho,u,v,E,p,a,w
  real,dimension(6,5) :: F_i,q_i,g_i,v_i
  real,dimension(5,5) :: R_i,Rinv_i
  real,dimension(1,5) :: q_h,hpr,hnr

  call pad_fluxes(F,Fnew,n_x,n_y,n_z,5)
  call pad_fluxes(q,qnew,n_x,n_y,n_z,5)


  select case (dir)
    case(1)
      do k=1,n_z
        do j=1,n_y
          do i=1,n_x
            F_i(1:6,:) = Fnew(i-2:i+3,j,k,:)
            q_i(1:6,:) = qnew(i-2:i+3,j,k,:)

            q_h(1,:) = (qnew(i,j,k,:) + qnew(i+1,j,k,:))/2.0

            call primitives(q_h,1,1,1,rho,u,v,w,E,p,a)
            R_i = Rcalc(u(1),a(1))
            Rinv_i = Rinv(u(1),a(1))
            call WENO(lambda,F_i,q_i,R_i,Rinv_i,hpr,hnr)
            hp(i,j,k,:) = hpr(1,:)
            hn(i,j,k,:) = hnr(1,:)
          end do
        end do
      end do
    case(2)
      do k=1,n_z
        do j=1,n_y
          do i=1,n_x
            F_i(1:6,:) = Fnew(i,j-2:j+3,k,:)
            q_i(1:6,:) = qnew(i,j-2:j+3,k,:)

            q_h(1,:) = (qnew(i,j,k,:) + qnew(i,j+1,k,:))/2.0

            call primitives(q_h,1,1,1,rho,u,v,w,E,p,a)
            R_i = Rcalc(u(1),a(1))
            Rinv_i = Rinv(u(1),a(1))
            call WENO(lambda,F_i,q_i,R_i,Rinv_i,hpr,hnr)
            hp(i,j,k,:) = hpr(1,:)
            hn(i,j,k,:) = hnr(1,:)
          end do
        end do
      end do
    case(3)
      do k=1,n_z
        do j=1,n_y
          do i=1,n_x
            F_i(1:6,:) = Fnew(i,j,k-2:k+3,:)
            q_i(1:6,:) = qnew(i,j,k-2:k+3,:)

            q_h(1,:) = (qnew(i,j,k,:) + qnew(i,j,k+1,:))/2.0

            call primitives(q_h,1,1,1,rho,u,v,w,E,p,a)
            R_i = Rcalc(u(1),a(1))
            Rinv_i = Rinv(u(1),a(1))
            call WENO(lambda,F_i,q_i,R_i,Rinv_i,hpr,hnr)
            hp(i,j,k,:) = hpr(1,:)
            hn(i,j,k,:) = hnr(1,:)
          end do
        end do
      end do
  end select

  return
end subroutine WENO53d
