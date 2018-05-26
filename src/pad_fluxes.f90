subroutine pad_fluxes(F,Fnew,n_x,n_y,n_z,vars)
  implicit none
  real, dimension(n_x,n_y,n_z,vars) :: F
  real, dimension(-1:n_x+3,-1:n_y+3,-1:n_z+3,vars) :: Fnew
  integer :: n_x,n_y,n_z,vars

  Fnew(1:n_x,1:n_y,1:n_z,:) = F(1:n_x,1:n_y,1:n_z,:)



  ! !!! Zero order extrapolation
  ! Fnew(-1:0,:,:,:) = spread(Fnew(1,:,:,:),1,2)
  ! Fnew(:,-1:0,:,:) = spread(Fnew(:,1,:,:),2,2)
  Fnew(:,:,-1:0,:) = spread(Fnew(:,:,1,:),3,2)
  Fnew(n_x+1:n_x+3,:,:,:) = spread(Fnew(n_x,:,:,:),1,3)
  Fnew(:,n_y+1:n_y+3,:,:) = spread(Fnew(:,n_y,:,:),2,3)
  Fnew(:,:,n_z+1:n_z+3,:) = spread(Fnew(:,:,n_z,:),3,3)



  !! Periodic BC
  Fnew(-1:0,:,:,:) = Fnew(n_x-1:n_x,:,:,:)
  Fnew(:,-1:0,:,:) = Fnew(:,n_y-1:n_y,:,:)
  ! Fnew(:,:,-1:0,:) = Fnew(:,:,n_z-1:n_z,:)
  ! Fnew(n_x+1:n_x+3,:,:,:) = Fnew(1:3,:,:,:)
  ! Fnew(:,n_y+1:n_y+3,:,:) = Fnew(:,1:3,:,:)
  ! Fnew(:,:,n_z+1:n_z+3,:) = Fnew(:,:,1:3,:)



end subroutine pad_fluxes
