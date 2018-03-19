module saver
  implicit none
contains
  function save_data(q,x,y,z,n_x,n_y,n_z,t,id)
    integer :: n_x,n_y,n_z,id,i,j,k,save_data
    real :: t,gamma=1.4
    real,dimension(n_x,n_y,n_z,5) :: q
    real,dimension(n_x,n_y,n_z) :: u,p,rho,E,a,M,v,Temp,w
    real, dimension(n_x) :: x
    real, dimension(n_y) :: y
    real, dimension(n_z) :: z
    character(len=5) :: charI
    character(len=1024) :: fname

    rho = q(:,:,:,1)
    u = q(:,:,:,2)/rho
    v = q(:,:,:,3)/rho
    w = q(:,:,:,4)/rho
    E = q(:,:,:,5)/rho
    p = (gamma-1)*rho*(E-0.5*(u*u+v*v+w*w))
    a = SQRT(gamma*p/rho)
    ! M = u/a
    Temp = p/rho

    write(charI,'(I5)') id
    write(fname,'(a,i6.6,a,f8.4,a)') "./data/frame",id,"t_",t,".dat"


    open(unit = 100,file = fname)
    do k=1,n_z
      do j=1,n_y
        do i=1,n_x
          write(100,*)x(i),y(j),z(k),u(i,j,k),v(i,j,k),w(i,j,k),p(i,j,k),rho(i,j,k),Temp(i,j,k)
        end do
      end do
    end do
    close(100)
  save_data=1
  end function save_data

  function get_video(val)
    integer :: val,get_video
    character(len=1024) :: property,command
    select case (val)
      case(1)
        property='u_velocity'
      case(2)
        property= 'v_velocity'
      case(3)
        property= 'w_velocity'
      case(4)
        property='pressure'
      case(5)
        property='density'
      case(6)
        property='Temperature'
      case default
        property='y'
    end select
    write(command,'(a,a,a,a,a)') "avconv -i ""./plots/",trim(property),"%06d.png"" -r 30 ./plots/",trim(property),".mp4"
    call system(command)
    get_video=1
  end function get_video
end module saver
