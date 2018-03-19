subroutine viscous_fluxes(vel,p,rho,n_x,n_y,n_z,delx,dely,delz,Suth,Re,Pr,tauxx,tauxy,tauxz,tauyy,tauyz,tauzz,q_x,q_y,q_z)
  implicit none
  real, dimension(n_x,n_y,n_z,4) :: vel     ! 1-xvel, 2-yvel, 3-zvel, 4-Temp
  real, dimension(n_x,n_y,n_z) :: p,rho,tauxx,tauxy,tauxz,tauyy,tauyz,tauzz,q_x,q_y,q_z
  real, dimension(-1:n_x+3,-1:n_y+3,-1:n_z+3,4) :: velnew
  real :: delx,dely,delz,Suth,R_gas_const = 287.0, T,meu,lambda,Re,Pr,Therm_const
  real :: u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,T_x,T_y,T_z
  real, dimension(1,4) :: temp_der,temp_der2
  real, dimension(6,4) :: temp,temp2
  real, parameter :: gamma=1.4
  integer :: n_x,n_y,n_z,i,j,k

  Therm_const = gamma / ((gamma-1.0)*Pr*Re)

  call pad_fluxes(vel,velnew,n_x,n_y,n_z,4)
  ! Calculate the derivatives in x and y directions
  do k=1,n_z
    do j=1,n_y
      do i =1,n_x

        temp(1:6,:) = velnew(i-2:i+3,j,k,:)
        call get_first_deriv(temp,delx,temp_der)
        u_x = temp_der(1,1); v_x = temp_der(1,2); w_x = temp_der(1,3); T_x = temp_der(1,4)

        temp(1:6,:) = velnew(i,j-2:j+3,k,:)
        call get_first_deriv(temp,dely,temp_der)
        u_y = temp_der(1,1); v_y = temp_der(1,2); w_y = temp_der(1,3); T_y = temp_der(1,4)

        temp(1:6,:) = velnew(i,j,k-2:k+3,:)
        call get_first_deriv(temp,delz,temp_der)
        u_z = temp_der(1,1); v_z = temp_der(1,2); w_z = temp_der(1,3); T_z = temp_der(1,4)

        T = vel(i,j,k,4)
        meu = (T**(3.0/2.0)) * ( (1+Suth)/(T+Suth) ) 
        lambda = (-2.0/3.0)*meu

        tauxx(i,j,k) = (lambda*(u_x+v_y+w_z) + 2*meu*u_x)/Re
        tauyy(i,j,k) = (lambda*(u_x+v_y+w_z) + 2*meu*v_y)/Re
        tauzz(i,j,k) = (lambda*(u_x+v_y+w_z) + 2*meu*w_z)/Re
        tauxy(i,j,k) = (meu*(u_y+v_x))/Re
        tauxz(i,j,k) = (meu*(u_z+w_x))/Re
        tauyz(i,j,k) = (meu*(w_y+v_z))/Re
        q_x(i,j,k) = T_x * Therm_const
        q_y(i,j,k) = T_y * Therm_const
        q_z(i,j,k) = T_z * Therm_const

      end do
    end do
end do

end subroutine viscous_fluxes


subroutine get_first_deriv(vel,delta,deriv)
  implicit none
  real, dimension(6,3) :: vel
  real :: delta
  real, dimension(1,3) :: deriv

  deriv(1,:) = ( (1.0/20.0)*vel(1,:) + (-0.5)*vel(2,:) + (-1.0/3.0)*vel(3,:) &
               +	vel(4,:) + (-0.25)*vel(5,:) + (1.0/30.0)*vel(6,:) )/delta

  ! deriv(1,:) = ( (1.0/12.0)*vel(1,:) -(2.0/3.0)*vel(2,:) +(2.0/3.0)*vel(4,:) -(1.0/12.0)*vel(5,:) )/(delta)


end subroutine get_first_deriv
