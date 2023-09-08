module mpi_vars
  ! This module keeps track of the
  ! variables associated with MPI
  integer :: myid
  integer :: numprocs,ierr,tag
  integer :: pmax(3), m_3did(3), m_off(3), m_3dneigh(3,2)
  integer :: i_start(3),i_end(3)
end module mpi_vars

MODULE loc_arrays
  ! Local arrays
  ! Grids
  double precision,dimension(:),allocatable :: x,y,z
  ! Solution arrays
  double precision, dimension(:,:,:), allocatable :: u,up,uproj,force,uproj_old
  double precision, dimension(:,:,:), allocatable :: v,vp,vproj,vproj_old
  double precision, dimension(:,:,:,:), allocatable :: kv,ku
END MODULE loc_arrays

program wave3d
  use mpi_vars
  use mpi
  use loc_arrays
  implicit none
  double precision, parameter :: pi = acos(-1.d0)
  integer BC_array(3,2) ! For setting boundary conditions
  integer :: nx,ny,nz,nt
  DOUBLE PRECISION :: tend, omega,d_omega
  DOUBLE PRECISION :: x_min,x_max,y_min,y_max,z_min,z_max
  integer :: n_olp

  integer :: m_nx,m_ny,m_nz
  ! for splitting the number of gridpoints
  integer :: nn,m_n,reminder,dim
  ! Generic integers
  integer :: i,j,k,p1,p2,p3,ii1,ii2,ii3,dir
  integer :: it,it_helm
  !
  DOUBLE PRECISION ::hx,hy,hz,dt,t,hi2,alp_bc
  DOUBLE PRECISION ::res_loc,residual
  !
  integer :: status(MPI_STATUS_SIZE),n_dble
  CHARACTER(7) :: charit
  CHARACTER(7) :: charx
  CHARACTER(7) :: chary
  CHARACTER(7) :: charz
  integer :: vol2linidx
  !
  real :: start_time

  DOUBLE PRECISION :: tol
  logical :: print_x_slice
  integer :: i_x_slice

  vol2linidx(p1,p2,p3) = (p1-1 + pmax(1)*(p2-1)+pmax(2)*pmax(1)*(p3-1))


  ! Set up BC
  BC_array = 0
  BC_array(1,1) = 1 ! at x_min
  BC_array(1,2) = 1 ! at x_max
  BC_array(2,1) = 1 ! at y_min
  BC_array(2,2) = 0 ! at y_max
  BC_array(3,1) = 1 ! at z_min
  BC_array(3,2) = 1 ! at z_max

  x_min = -1.d0
  y_min = -1.d0
  z_min = -1.d0
  x_max =  1.d0
  y_max =  1.d0
  z_max =  1.d0

  call cpu_time(start_time)
  ! initialize the MPI world
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  ! if (myid .eq. 0) WRITE(*,*) "Starting Wave-3D with :", numprocs, " processors."
  ! if (myid .eq. 0) open(21,file='results.txt',status='unknown')

  n_olp = 1

  tol = 5.0d-5
  d_omega = 1.0d0
  omega = 20.0d0

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if(myid.eq.0) nx = max(ceiling(10.d0*omega),20)
  if (myid .eq. 0) WRITE(*,*) "Using nx = ", nx, " points."
  call MPI_Bcast(nx,1,MPI_INT,0,MPI_Comm_world,ierr)
  ny = nx
  nz = nx
  ! Distribute the processors in 3D
  call PROCDIST3D( nx, ny, nz, numprocs, pmax )
  ! my3did is in [1,pmax(1)] x [1,pmax(2)] x [1,pmax(3)]
  call MYMESH( myid, pmax, m_3did )
  !
  if ((nx .lt. pmax(1)) .or. (ny .lt. pmax(2)) .or. (nz .lt. pmax(3))) then
     if(myid == 0) write(*,*) 'Too few gridpoints!'
     call mpi_abort(MPI_COMM_WORLD,0,ierr)
  end if
  if (myid .eq. 0) WRITE(*,*) "The grid is distributed on :  ", pmax, " processors."
  ! neighbors for the MPI communication
  do dim = 1,3
     ii1 =  (2-dim)*(3-dim)/2
     ii2 = -(1-dim)*(3-dim)
     ii3 =  (1-dim)*(2-dim)/2
     do dir = -1,1,2
        if (((m_3did(dim).eq.1) .and. (dir .lt. 0)).or.&
             ((m_3did(dim).eq.pmax(dim)) .and. (dir .gt. 0))) then
           m_3dneigh(dim,1+(dir+1)/2) = vol2linidx(&
                m_3did(1)-dir*ii1*(pmax(1)-1),&
                m_3did(2)-dir*ii2*(pmax(2)-1),&
                m_3did(3)-dir*ii3*(pmax(3)-1))
        else
           m_3dneigh(dim,1+(dir+1)/2) = vol2linidx(&
                m_3did(1)+dir*ii1,&
                m_3did(2)+dir*ii2,&
                m_3did(3)+dir*ii3)
        end if
     end do
  end do
  ! Set BC for the communication.
  do dim = 1,3
     if (m_3did(dim).eq.1) m_3dneigh(dim,1) = MPI_PROC_NULL
     if (m_3did(dim).eq.pmax(dim)) m_3dneigh(dim,2) = MPI_PROC_NULL
  end do
  !
  ! Now we compute how many points per processor
  ! and where they start
  !
  hx = (x_max-x_min)/dble(nx)
  hy = (y_max-y_min)/dble(ny)
  hz = (z_max-z_min)/dble(nz)

  hi2 = 1.d0/hx**2
  !
  do dim = 1,3
     ! Number of gridpoints to be distributed.
     ! 0:nx
     if(dim.eq.1) nn = nx+1
     if(dim.eq.2) nn = ny+1
     if(dim.eq.3) nn = nz+1
     ! Eeach processor will have at least (floor(nn/pmax(1)) + 1) gridpoints
     ! The +1 is due to 0:m_nx numbering
     m_n = nn / pmax(dim)
     reminder = mod(nn,pmax(dim))
     if (m_3did(dim) .le. reminder ) then
        m_n = m_n+1
        m_off(dim) = (m_3did(dim)-1)*(m_n)
     else
        m_off(dim) = (reminder)*(m_n+1)+(m_3did(dim)-1-reminder)*(m_n)
     end if
     if(dim.eq.1) m_nx = m_n
     if(dim.eq.2) m_ny = m_n
     if(dim.eq.3) m_nz = m_n
  end do

  i_start = 1
  i_end(1) = m_nx
  i_end(2) = m_ny
  i_end(3) = m_nz
  ! Set BC
  do dim = 1,3
     if ((m_3did(dim).eq.1).and.(dim < 3)) i_start(dim) = 2
     if (m_3did(dim).eq.pmax(dim)) i_end(dim) = i_end(dim)-1
  end do

  ! Now allocate various arrays
  allocate(x(m_nx),y(m_ny),z(m_nz))
  allocate(u(-n_olp+1:m_nx+n_olp,-n_olp+1:m_ny+n_olp,-n_olp+1:m_nz+n_olp),&
       ku(m_nx,m_ny,m_nz,4),up(m_nx,m_ny,m_nz),force(m_nx,m_ny,m_nz))
  allocate(uproj(i_start(1):i_end(1),i_start(2):i_end(2),i_start(3):i_end(3)))
  allocate(uproj_old(i_start(1):i_end(1),i_start(2):i_end(2),i_start(3):i_end(3)))

  allocate(v(m_nx,m_ny,m_nz),kv(m_nx,m_ny,m_nz,4),vp(m_nx,m_ny,m_nz))
  allocate(vproj(i_start(1):i_end(1),i_start(2):i_end(2),i_start(3):i_end(3)))
  allocate(vproj_old(i_start(1):i_end(1),i_start(2):i_end(2),i_start(3):i_end(3)))

  ! Grids
  do i = 1,m_nx
     x(i) = x_min + hx*dble(m_off(1)+i-1)
  end do
  do i = 1,m_ny
     y(i) = y_min + hy*dble(m_off(2)+i-1)
  end do
  do i = 1,m_nz
     z(i) = z_min + hz*dble(m_off(3)+i-1)
  end do

  tend = 2.d0*pi/omega

  dt = 0.4d0*hx
  nt = ceiling(tend/dt)
  dt = tend/dble(nt)

  do k = i_start(3),i_end(3)
     do j = i_start(2),i_end(2)
        do i = i_start(1),i_end(1)
           force(i,j,k) = omega**3*exp(-36*omega**2*((x(i)-0.01d0)**2+(y(j)-0.012d0)**2+(z(k)-0.005d0)**2))
        end do
     end do
  end do

  uproj = 0.d0
  vproj = 0.d0
  it_helm = 0
  do while (it_helm < 10)
     !!
!!$      if(mod(it_helm,2).eq.0) then
!!$         sprod_loc = sum(r_t_0*v_0)
!!$         call MPI_Allreduce(sprod_loc,alpha_0,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!!$         alpha_0 = rho_0/alpha_0
!!$         alpha_1 = alpha_0
!!$         u_1 = u_0 -alpha_0*v_0
!!$      end if

     call evolve_and_project
     res_loc = sum((uproj-uproj_old)**2)
     call MPI_Allreduce(res_loc,residual,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     ! write(*,*) it_helm,sqrt(dble(it_helm+1))*abs(tau_0)
     if(myid.eq.0) write(*,*) it_helm,sqrt(residual)
     it_helm = it_helm + 1
  end do

  ! if(myid.eq.0) close(23)

  print_x_slice = .false.
  do i = 1,m_nx
     if ((x(i) .gt. 0.1d0) .and. (abs(x(i)-0.1d0) .lt. hx)) then
        print_x_slice = .true.
        i_x_slice = i
     end if
  end do
  if(print_x_slice) then
     WRITE(charit,"(I7.7)") 0
     WRITE(charx,"(I7.7)") m_3did(1)
     WRITE(chary,"(I7.7)") m_3did(2)
     WRITE(charz,"(I7.7)") m_3did(3)
     call printdble2d(uproj(i_x_slice,:,:),i_start(2),i_end(2),i_start(3),i_end(3),&
          "uxs"//charit//"_"//charx//"_"//chary//"_"//charz//".txt")
     call printdble2d(vproj(i_x_slice,:,:),i_start(2),i_end(2),i_start(3),i_end(3),&
          "vxs"//charit//"_"//charx//"_"//chary//"_"//charz//".txt")
  end if

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  deallocate(x,y,z)
  deallocate(u,ku,up,force)
  deallocate(v,kv,vp)
  deallocate(uproj,vproj)
  deallocate(uproj_old,vproj_old)

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call MPI_FINALIZE(ierr)

contains

  subroutine evolve_and_project

    uproj_old = uproj
    vproj_old = vproj
    u = 0.d0
    ! Set initial data
    do k = i_start(3),i_end(3)
       do j = i_start(2),i_end(2)
          do i = i_start(1),i_end(1)
             u(i,j,k) = uproj(i,j,k)
             v(i,j,k) = vproj(i,j,k)
          end do
       end do
    end do

    t = 0.d0
    ! Start integrating
    do k = i_start(3),i_end(3)
       do j = i_start(2),i_end(2)
          do i = i_start(1),i_end(1)
             uproj(i,j,k) = 0.5d0*dt*u(i,j,k)*(cos(omega*t)-0.25d0)
             vproj(i,j,k) = 0.5d0*dt*v(i,j,k)*(cos(omega*t)-0.25d0)
          end do
       end do
    end do

    alp_bc = -1.d0

    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    ! Now timestep with rk4
    do it = 1,nt
       t = dble(it-1)*dt

       do k = i_start(3),i_end(3)
          do j = i_start(2),i_end(2)
             do i = i_start(1),i_end(1)
                up(i,j,k) = u(i,j,k)
                vp(i,j,k) = v(i,j,k)
             end do
          end do
       end do

       call communicate_u
       ! Fill in BC
       call enforce_bc

       do k = i_start(3),i_end(3)
          do j = i_start(2),i_end(2)
             do i = i_start(1),i_end(1)
                kv(i,j,k,1) = (hi2*(-6.d0*u(i,j,k) + u(i+1,j,k) + u(i-1,j,k) &
                     + u(i,j+1,k) + u(i,j-1,k) + u(i,j,k+1) + u(i,j,k-1)) &
                     -cos(omega*t)*force(i,j,k))
                ku(i,j,k,1) = v(i,j,k)
             end do
          end do
       end do
       do k = i_start(3),i_end(3)
          do j = i_start(2),i_end(2)
             do i = i_start(1),i_end(1)
                v(i,j,k) = vp(i,j,k) + 0.5d0*dt*kv(i,j,k,1)
                u(i,j,k) = up(i,j,k) + 0.5d0*dt*ku(i,j,k,1)
             end do
          end do
       end do

       t = t + 0.5d0*dt

       call communicate_u
       ! Fill in BC
       call enforce_bc

       do k = i_start(3),i_end(3)
          do j = i_start(2),i_end(2)
             do i = i_start(1),i_end(1)
                kv(i,j,k,2) = (hi2*(-6.d0*u(i,j,k) + u(i+1,j,k) + u(i-1,j,k) &
                     + u(i,j+1,k) + u(i,j-1,k) + u(i,j,k+1) + u(i,j,k-1)) &
                     -cos(omega*t)*force(i,j,k))
                ku(i,j,k,2) = v(i,j,k)
             end do
          end do
       end do
       do k = i_start(3),i_end(3)
          do j = i_start(2),i_end(2)
             do i = i_start(1),i_end(1)
                v(i,j,k) = vp(i,j,k) + 0.5d0*dt*kv(i,j,k,2)
                u(i,j,k) = up(i,j,k) + 0.5d0*dt*ku(i,j,k,2)
             end do
          end do
       end do

       call communicate_u
       ! Fill in BC
       call enforce_bc

       do k = i_start(3),i_end(3)
          do j = i_start(2),i_end(2)
             do i = i_start(1),i_end(1)
                kv(i,j,k,3) = (hi2*(-6.d0*u(i,j,k) + u(i+1,j,k) + u(i-1,j,k) &
                     + u(i,j+1,k) + u(i,j-1,k) + u(i,j,k+1) + u(i,j,k-1)) &
                     -cos(omega*t)*force(i,j,k))
                ku(i,j,k,3) = v(i,j,k)
             end do
          end do
       end do
       do k = i_start(3),i_end(3)
          do j = i_start(2),i_end(2)
             do i = i_start(1),i_end(1)
                v(i,j,k) = vp(i,j,k) + dt*kv(i,j,k,3)
                u(i,j,k) = up(i,j,k) + dt*ku(i,j,k,3)
             end do
          end do
       end do

       t = t + 0.5d0*dt
       call communicate_u
       ! Fill in BC
       call enforce_bc

       do k = i_start(3),i_end(3)
          do j = i_start(2),i_end(2)
             do i = i_start(1),i_end(1)
                kv(i,j,k,4) = (hi2*(-6.d0*u(i,j,k) + u(i+1,j,k) + u(i-1,j,k) &
                     + u(i,j+1,k) + u(i,j-1,k) + u(i,j,k+1) + u(i,j,k-1)) &
                     -cos(omega*t)*force(i,j,k))
                ku(i,j,k,4) = v(i,j,k)
             end do
          end do
       end do

       do k = i_start(3),i_end(3)
          do j = i_start(2),i_end(2)
             do i = i_start(1),i_end(1)
                v(i,j,k) = vp(i,j,k) + (dt/6.d0)*(kv(i,j,k,1) + 2.d0*(kv(i,j,k,2) + kv(i,j,k,3)) + kv(i,j,k,4))
                u(i,j,k) = up(i,j,k) + (dt/6.d0)*(ku(i,j,k,1) + 2.d0*(ku(i,j,k,2) + ku(i,j,k,3)) + ku(i,j,k,4))
             end do
          end do
       end do
       ! Do some more projection
       do k = i_start(3),i_end(3)
          do j = i_start(2),i_end(2)
             do i = i_start(1),i_end(1)
                uproj(i,j,k) = uproj(i,j,k) + dt*u(i,j,k)*(cos(omega*t)-0.25d0)
                vproj(i,j,k) = vproj(i,j,k) + dt*v(i,j,k)*(cos(omega*t)-0.25d0)
             end do
          end do
       end do

       !     WRITE(charit,"(I7.7)") it
       !     call printdble2d(u(10,i_start(2):i_end(2),i_start(3):i_end(3)),i_start(2),i_end(2),i_start(3),i_end(3),"u"//charit//".txt")

    end do
    ! Remove the last half
    do k = i_start(3),i_end(3)
       do j = i_start(2),i_end(2)
          do i = i_start(1),i_end(1)
             uproj(i,j,k) = (2.d0/tend)*(uproj(i,j,k)-0.5d0*dt*u(i,j,k)*(cos(omega*tend)-0.25d0))
             vproj(i,j,k) = (2.d0/tend)*(vproj(i,j,k)-0.5d0*dt*v(i,j,k)*(cos(omega*tend)-0.25d0))
          end do
       end do
    end do
    call MPI_Barrier(MPI_COMM_WORLD,ierr)

  end subroutine evolve_and_project

  subroutine enforce_bc

    ! Fill in BC
    ! THE ZERO DIRICHLET ARE NEVER CHANGED...

    ! x = x_min and x_max
    if ((m_3did(1).eq.1).and.(bc_array(1,1).eq.1)) then
       i = i_start(1)
       do k = i_start(3),i_end(3)
          do j = i_start(2),i_end(2)
             u(i-1,j,k) = u(i+1,j,k) + alp_bc*2.d0*hx*v(i,j,k)
          end do
       end do
    end if

    if ((m_3did(3).eq.pmax(3)).and.(bc_array(1,2).eq.1)) then
       i = i_end(1)
       do k = i_start(3),i_end(3)
          do j = i_start(2),i_end(2)
             u(i+1,j,k) = u(i-1,j,k) + alp_bc*2.d0*hx*v(i,j,k)
          end do
       end do
    end if

    ! y = y_min and y_max
    if ((m_3did(2).eq.1).and.(bc_array(2,1).eq.1)) then
       j = i_start(2)
       do k = i_start(3),i_end(3)
          do i = i_start(1),i_end(1)
             u(i,j-1,k) = u(i,j+1,k) + alp_bc*2.d0*hx*v(i,j,k)
          end do
       end do
    end if

    if ((m_3did(3).eq.pmax(2)).and.(bc_array(2,2).eq.1)) then
       j = i_end(2)
       do k = i_start(3),i_end(3)
          do i = i_start(1),i_end(1)
             u(i,j+1,k) = u(i,j-1,k) + alp_bc*2.d0*hx*v(i,j,k)
          end do
       end do
    end if

    if ((m_3did(3).eq.1).and.(bc_array(3,1).eq.1)) then
       k = i_start(3)
       do j = i_start(2),i_end(2)
          do i = i_start(1),i_end(1)
             u(i,j,k-1) = u(i,j,k+1) + alp_bc*2.d0*hx*v(i,j,k)
          end do
       end do
    end if

    if ((m_3did(3).eq.pmax(3)).and.(bc_array(3,2).eq.1)) then
       k = i_end(3)
       do j = i_start(2),i_end(2)
          do i = i_start(1),i_end(1)
             u(i,j,k+1) = u(i,j,k-1) + alp_bc*2.d0*hx*v(i,j,k)
          end do
       end do
    end if



  end subroutine enforce_bc

  subroutine communicate_u
    ! Communicate the solution to friends and neighbours
    tag = 1
    n_dble = n_olp*(m_ny+2*n_olp)*(m_nz+2*n_olp)
    ! send to the neighbor to the right receive from left
    call MPI_Sendrecv(&
         u(m_nx-(n_olp-1):m_nx,-n_olp+1:m_ny+n_olp,-n_olp+1:m_nz+n_olp),n_dble,&
         MPI_DOUBLE_PRECISION,m_3dneigh(1,2),tag,&
         u(1-n_olp:0,-n_olp+1:m_ny+n_olp,-n_olp+1:m_nz+n_olp),n_dble,&
         MPI_DOUBLE_PRECISION,m_3dneigh(1,1),tag,MPI_COMM_WORLD,status,ierr)
    ! send to the neighbor to the left receive from right
    call MPI_Sendrecv(&
         u(1:n_olp,-n_olp+1:m_ny+n_olp,-n_olp+1:m_nz+n_olp),n_dble,&
         MPI_DOUBLE_PRECISION,m_3dneigh(1,1),tag,&
         u(m_nx+1:m_nx+n_olp,-n_olp+1:m_ny+n_olp,-n_olp+1:m_nz+n_olp),n_dble,&
         MPI_DOUBLE_PRECISION,m_3dneigh(1,2),tag,MPI_COMM_WORLD,status,ierr)

    n_dble = n_olp*(m_nx+2*n_olp)*(m_nz+2*n_olp)
    ! y-top -> y-bot
    call MPI_Sendrecv(&
         u(-n_olp+1:m_nx+n_olp,m_ny-(n_olp-1):m_ny,-n_olp+1:m_nz+n_olp),n_dble,&
         MPI_DOUBLE_PRECISION,m_3dneigh(2,2),tag,&
         u(-n_olp+1:m_nx+n_olp,1-n_olp:0,-n_olp+1:m_nz+n_olp),n_dble,&
         MPI_DOUBLE_PRECISION,m_3dneigh(2,1),tag,MPI_COMM_WORLD,status,ierr)
    ! y-bot -> y-top
    call MPI_Sendrecv(&
         u(-n_olp+1:m_nx+n_olp,1:n_olp,-n_olp+1:m_nz+n_olp),n_dble,&
         MPI_DOUBLE_PRECISION,m_3dneigh(2,1),tag,&
         u(-n_olp+1:m_nx+n_olp,m_ny+1:m_ny+n_olp,-n_olp+1:m_nz+n_olp),n_dble,&
         MPI_DOUBLE_PRECISION,m_3dneigh(2,2),tag,MPI_COMM_WORLD,status,ierr)

    n_dble = n_olp*(m_nx+2*n_olp)*(m_ny+2*n_olp)
    ! z-top -> z-bot
    call MPI_Sendrecv(&
         u(-n_olp+1:m_nx+n_olp,-n_olp+1:m_ny+n_olp,m_nz-(n_olp-1):m_nz),n_dble,&
         MPI_DOUBLE_PRECISION,m_3dneigh(3,2),tag,&
         u(-n_olp+1:m_nx+n_olp,-n_olp+1:m_ny+n_olp,1-n_olp:0),n_dble,&
         MPI_DOUBLE_PRECISION,m_3dneigh(3,1),tag,MPI_COMM_WORLD,status,ierr)
    ! z-bot -> z-top
    call MPI_Sendrecv(&
         u(-n_olp+1:m_nx+n_olp,-n_olp+1:m_ny+n_olp,1:n_olp),n_dble,&
         MPI_DOUBLE_PRECISION,m_3dneigh(3,1),tag,&
         u(-n_olp+1:m_nx+n_olp,-n_olp+1:m_ny+n_olp,m_nz+1:m_nz+n_olp),n_dble,&
         MPI_DOUBLE_PRECISION,m_3dneigh(3,2),tag,MPI_COMM_WORLD,status,ierr)
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
  end subroutine communicate_u

end program wave3d

subroutine printdble1d(u,nx1,nx2,str)
  implicit none
  integer, intent(in) :: nx1,nx2
  DOUBLE PRECISION, intent(in) :: u(nx1:nx2)
  character(len=*), intent(in) :: str
  integer :: i
  open(2,file=trim(str),status='unknown')
  do i=nx1,nx2,1
     write(2,fmt='(E24.16)',advance='no') u(i)
  end do
  write(2,'()')
  close(2)
end subroutine printdble1d

subroutine printdble2d(u,nx1,nx2,ny1,ny2,str)
  implicit none
  integer, intent(in) :: nx1,nx2,ny1,ny2
  DOUBLE PRECISION, intent(in) :: u(nx1:nx2,ny1:ny2)
  character(len=*), intent(in) :: str
  integer :: i,j
  open(2,file=trim(str),status='unknown')
  do j=ny1,ny2,1
     do i=nx1,nx2,1
        if(abs(u(i,j)) .lt. 1.0d-18 ) then
           write(2,fmt='(E24.16)',advance='no') 0.d0
        else
           write(2,fmt='(E24.16)',advance='no') u(i,j)
        end if
     end do
     write(2,'()')
  end do
  close(2)
end subroutine printdble2d
