module scalars_module
! HUMIDITY subroutines in place but not yet turned on !!
use types,only:rprec
use param 
use sim_param,only:u,v,w,theta,q,path,pcon,dudt,dvdt,dwdt
use bottombc !Includes patches subroutine
use sgsmodule,only:Nu_t,magS
implicit none
integer, parameter:: tag_counter = 200
logical, parameter:: SCALAR_DEBUG=.FALSE.
!!!!!!--------------------------------------------
! Part I of the scalar files - contains the basic subroutines
! Also look at scalars_module2.f90 for other subroutines !! 
! CONTAINS subroutines :
! theta_all_in_one - Performs the various steps for theta calculation
! humidity_all_in_one - Performs the various steps for humidity
! scalar_RHS_calc - computes the RHS side of the scalar evolution equation
! calcbeta - computes the buoyancy term for temperature
! step_scalar - time evolves the scalar
! obukhov - computes the obukhov similarity terms for use in scalars,wallstress and derivwall
! Authored by Vijayant Kumar
! Last modified - April 24, 2004
!!!!!!--------------------------------------------
!real(kind=rprec),dimension(ld,ny,nz):: scalar
!real(kind=rprec),dimension(ld,ny,nz):: theta,q ! theta and q specified in sim_param
 
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
real(kind=rprec),dimension(ld,ny,$lbz:nz):: beta_scal,Pr_
!real(kind=rprec),dimension(ld,ny,nz):: dsdx,dsdy,dsdz
real(kind=rprec),dimension(ld,ny,$lbz:nz):: dTdz,dqdz ! Only ones needed for output
! Might need to add x and y derivatives here in case they need to be outputted
! Right now they are in the "scalar"_all_in_one routines below !!
real(kind=rprec),dimension(ld,ny,$lbz:nz):: RHS_Tf,RHS_T,RHS_qf,RHS_q
real(kind=rprec), dimension(ld,ny,$lbz:nz):: sgs_t3,sgs_q3 !defines the surface sgs flux

!real(kind=rprec),dimension(:,:,:),allocatable:: scalar,q
!real(kind=rprec),parameter::g=9.81,inv_strength=0.008 !inversion strength (K/m)
real(kind=rprec),dimension(nx,ny)::L,wstar !defines obukhov length and convective vel scale, w_star
!real(kind=rprec),dimension(nx,ny)::z_os !T_s and q_s are defined in bottombc.f90
real(kind=rprec),dimension(ld,ny)::T_s_filtered !filtered T_s for calc of wT_s

!real(kind=rprec),dimension(nx,ny)::ustar_avg ! Defines the local u_star as calculated in obukhov

! Now define local u_star in bottombc.f90 and calculate in wallstress.f90 and use that value
! everywhere else
integer, parameter:: obukhov_output=0 !Controls whether obukhov variables are outputted by scalar_slice
integer, parameter:: wt_s_vector_dim1=no_days*86400/300+1
!integer, parameter:: wt_s_vector_dim1=floor(no_days*86400._rprec/300._rprec)+1

!real(kind=rprec),dimension(floor(no_days*86400._rprec/300._rprec)+1,1) :: wt_s_vector
real(kind=rprec),dimension(wt_s_vector_dim1,1) :: wt_s_vector
! Variables for heterogeneity analysis
! hetero_array_freqz = number of time steps equivalent to 20 seconds
!integer,parameter:: hetero_array_freqz=int(20/dt_dim),hetero_count_out=p_count
integer,parameter:: hetero_array_freqz=100,hetero_count_out=p_count
integer,save::time_ind


! Variables added for pollen
! Chamecki - 08/01/2006
REAL(kind=rprec),DIMENSION(ld,ny,$lbz:nz)  :: Kc_t	    ! 3D matrix of SGS diffusivity for pollen
REAL(kind=rprec),DIMENSION(ld,ny,$lbz:nz)  :: dPCondz	    ! Vertical derivatives (needed also for output)
REAL(kind=rprec),DIMENSION(ld,ny,$lbz:nz)  :: RHS_PConf,RHS_PCon ! RHS for PCon equation
REAL(kind=rprec),DIMENSION(ld,ny,$lbz:nz)  :: sgs_PCon3     ! Defines the sgs vertical flux
REAL(kind=rprec),DIMENSION(ld,ny,$lbz:nz)  :: res_PCon3     ! Defines the resolved vertical flux
REAL(kind=rprec),DIMENSION(ld,ny,$lbz:nz)  :: sgs_PCon1     ! Defines the sgs x-direction flux
REAL(kind=rprec),DIMENSION(ld,ny,$lbz:nz)  :: Cs2Sc         ! Dynamic coefficient for scalar equation

REAL(kind=rprec),DIMENSION(nx,ny)	   :: P_surf_flux   ! Surface pollen flux     
REAL(kind=rprec),DIMENSION(nx,ny)	   :: deposition    ! Surface pollen deposition (this is actually the net flux=deposition-source)
REAL(kind=rprec),DIMENSION(nx,ny)	   :: P_surf_flux_dep! Surface pollen flux for Cr=0 everywhere
REAL(kind=rprec),DIMENSION(nx,ny)	   :: Real_dep      ! This is the real deposition (using Cr=0 everywhere)
REAL(kind=rprec) 	  		   :: flux_out,flux_out_prev

! Variables for interpolation of velocity field
REAL(kind=rprec),DIMENSION(nx,nx)          :: matrix_x      ! Matrix for x direction
REAL(kind=rprec),DIMENSION(ny,ny)          :: matrix_y      ! Matrix for y direction
REAL(kind=rprec),DIMENSION(nx)             :: dvector_x     ! Derivative vector for x direction
REAL(kind=rprec),DIMENSION(ny)             :: dvector_y     ! Derivative vector for y direction

! Variables for pollen balance
REAL(kind=rprec):: released,airborne,deposited,gone

! Scaling factors for time-varying u* and Co
REAL(kind=rprec):: scale_us,scale_Co
REAL(kind=rprec),DIMENSION(3):: scale_Co3

!++++++++++++++++++++++++++++++++++
!DY Variables added by Di Yang
REAL(kind=rprec),dimension(ld,ny,$lbz:nz):: beta_pcon
!++++++++++++++++++++++++++++++++++




contains
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
subroutine theta_all_in_one
use topbc,only:sponge
use test_filtermodule
implicit none
real:: wt_s_current,dummy_t
integer :: jz,ios,counter
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
real(kind=rprec),dimension(ld,ny,$lbz:nz)::dTdx,dTdy
real(kind=rprec),dimension($lbz:nz)::sponge_theta


  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    wt_s_current=wt_s
  end if

  ! Right now set the Prandtl num matrix equal to a constant Prandtl
  Pr_=Pr 

  call filt_da(theta,dTdx,dTdy)

  if ((USE_MPI) .AND. (SCALAR_DEBUG)) then
    do jz=0,nz
      print *,'jz,theta,dTdx,dTdy,coord',jz,theta(4,4,jz),dTdx(4,4,jz),&
      dTdy(4,4,jz),coord
    end do
  elseif (SCALAR_DEBUG) then
    do jz=0,nz
      print *,'jz,theta,dTdx,dTdy',jz,theta(4,4,jz),dTdx(4,4,jz),dTdy(4,4,jz)
    end do
  end if

  call ddz_uv (dTdz,theta)

  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
    dTdz(:,:,Nz)=dTdz_top/T_scale*z_i ! Valid for temperature
  end if


$if ($MPI)
  !print *,'synchronizing in theta_all_in for coord = ',coord
  ! Need to synchronize w and dTdz across the processors for uvp node
  ! based computation in scalar_RHS_calc (calculation of RHS_m)
  call mpi_sendrecv (w(1, 1, 1), ld*ny, MPI_RPREC, down, tag_counter+1,  &
                     w(1, 1, nz), ld*ny, MPI_RPREC, up, tag_counter+1,   &
                     comm, status, ierr)   
  call mpi_sendrecv (dTdz(1, 1, 1), ld*ny, MPI_RPREC, down, tag_counter+3,  &
                     dTdz(1, 1, nz), ld*ny, MPI_RPREC, up, tag_counter+3,   &
                     comm, status, ierr)   

  ! Also need to synchronize Nu_t across the processors for 
  ! computation in scalar_RHS_calc (calculation of RHS_m)
  call mpi_sendrecv (Nu_t(1, 1, 1), ld*ny, MPI_RPREC, down, tag_counter+3,  &
                     Nu_t(1, 1, nz), ld*ny, MPI_RPREC, up, tag_counter+3,   &
                     comm, status, ierr)   
  call mpi_sendrecv (Nu_t(1, 1, nz-1), ld*ny, MPI_RPREC, up, tag_counter+4,  &
                     Nu_t(1, 1, 0), ld*ny, MPI_RPREC, down, tag_counter+4,   &
                     comm, status, ierr)
$endif


  if ((USE_MPI) .AND. (SCALAR_DEBUG)) then
    do jz=0,nz
      print *,'jz,dTdz,Nu_t,coord',jz,dTdz(4,4,jz),&
      Nu_t(4,4,jz),coord
    end do
  elseif (SCALAR_DEBUG) then
    do jz=0,nz
      print *,'jz,dTdz,Nu_t',jz,dTdz(4,4,jz),Nu_t(4,4,jz)
    end do
  end if

  if (S_FLAG) then
    RHS_Tf=RHS_T

    ! Perform test filtering of T_s for calculation of surf fluxes
    if ((jt_total .eq. SCAL_init) .and. (lbc .eq. 0) .and. (remote_homog_flag .eq. 0)) then
      print *,'T_s b4 filtering',sqrt(sum((T_s-sum(T_s)/float(nx*ny))**2))/float(nx*ny)
      T_s_filtered(1:nx,1:ny)=T_s
      call test_filter(T_s_filtered,G_test)
      T_s=T_s_filtered(1:nx,1:ny)
      print *,'T_s after filtering',sqrt(sum((T_s-sum(T_s)/float(nx*ny))**2))/float(nx*ny)
    end if
    if ((jt_total .eq. nsteps-1) .and. (lbc .eq. 0) .and. (remote_homog_flag .eq. 0)) then
      print *,'T_s after filtering',sqrt(sum((T_s-sum(T_s)/float(nx*ny))**2))/float(nx*ny)
    end if

    call scalar_RHS_calc(theta,dTdx,dTdy,dTdz,T_s,z_os,RHS_T,sgs_t3,wt_s_current)

    if (ubc==1 .and. damping_method==2) then !add the damping term to the scalar equation
      do jz=1,nz-1
        RHS_T(1:nx,1:ny,jz)=RHS_T(1:nx,1:ny,jz)-0.5_rprec*(sponge(jz)+sponge(jz+1))*&
                           (theta(1:nx,1:ny,jz)-sum(theta(1:nx,1:ny,jz))/(nx*ny))
      end do
    end if

    if ((USE_MPI) .AND. (SCALAR_DEBUG)) then
      do jz=0,nz
        print *,'jz,RHS_Tf@AFTER_RHS,RHS_T,coord',jz,RHS_Tf(4,4,jz),RHS_T(4,4,jz),coord
      end do
    elseif (SCALAR_DEBUG) then
      do jz=0,nz
        print *,'jz,RHS_Tf@AFTER_RHS,RHS_T',jz,RHS_Tf(4,4,jz),RHS_T(4,4,jz)
      end do
    end if

    ! Calculates the buoyancy term which gets added to the vertical momentum equation
    call calcbeta(theta) 

    if (SCALAR_DEBUG) then
      print *,'forcing beta_scal=0 to decouple'
      beta_scal=0._rprec
    end if

    if ((USE_MPI) .AND. (SCALAR_DEBUG)) then
      do jz=0,nz
        print *,'jz,sgs_t3,beta_scal,coord',jz,sgs_t3(4,4,jz),beta_scal(4,4,jz),coord
      end do
    elseif (SCALAR_DEBUG) then
      do jz=0,nz
        print *,'jz,sgs_t3,beta_scal',jz,sgs_t3(4,4,jz),beta_scal(4,4,jz)
      end do
    end if

    if (jt .eq. SCAL_INIT .and. (.not. initsc)) then
      RHS_Tf=RHS_T
    end if

    if ((USE_MPI) .AND. (SCALAR_DEBUG)) then
      do jz=0,nz
        print *,'jz,RHS_Tf,RHS_T,coord',jz,RHS_Tf(4,4,jz),RHS_T(4,4,jz),coord
      end do
    elseif (SCALAR_DEBUG) then
      do jz=0,nz
        print *,'jz,RHS_Tf,RHS_T',jz,RHS_Tf(4,4,jz),RHS_T(4,4,jz)
      end do
    end if

    call step_scalar(theta,RHS_T,RHS_Tf)
    
  end if



  if ((USE_MPI) .AND. (SCALAR_DEBUG)) then
    do jz=0,nz
      print *,'jz,theta@after step,coord',jz,theta(4,4,jz),coord
    end do
  elseif (SCALAR_DEBUG) then
    do jz=0,nz
      print *,'jz,theta@after step',jz,theta(4,4,jz)
    end do
  end if


$if ($MPI)
  call mpi_sendrecv (theta(1, 1, 1), ld*ny, MPI_RPREC, down, tag_counter+7,  &
                     theta(1, 1, nz), ld*ny, MPI_RPREC, up, tag_counter+7,   &
                     comm, status, ierr)   
  call mpi_sendrecv (theta(1, 1, nz-1), ld*ny, MPI_RPREC, up, tag_counter+8,  &
                     theta(1, 1, 0), ld*ny, MPI_RPREC, down, tag_counter+8,   &
                     comm, status, ierr)
$endif


  if ((USE_MPI) .AND. (SCALAR_DEBUG)) then
    do jz=0,nz
      print *,'jz,theta@after step&sync,coord',jz,theta(4,4,jz),coord
    end do
  elseif (SCALAR_DEBUG) then
    do jz=0,nz
      print *,'jz,theta@after step&sync',jz,theta(4,4,jz)
    end do
  end if

end subroutine theta_all_in_one


!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
!subroutine humidity_all_in_one
!! Not correct now.. DONT USE ....
!use sim_param
!!use scalars_module
!implicit none
!real(kind=rprec),dimension(ld,ny,nz)::dqdx,dqdy
!real:: wq_s_current
!!integer:: jt

!wq_s_current=wq_s

!Pr_=Pr !Right now set the Prandtl no matrix equal to a constant Prandtl
!! number as specified in param. could use the dynamic model ideal to compute Pr as well !!
!! The plan-averaged dynamic Prandtl number model is already coded. just need to put it in !!

!call filt_da(q,dqdx,dqdy)
!call ddz_uv (dqdz,q) ! Calculate vertical derivatives !! 
!dqdz(:,:,Nz)=0 ! Valid for humidity and other passive scalars

!if (S_FLAG) then
! RHS_q=RHS_q
!! q_s is the surface humidity - specified using patch_or_remote()
!call scalar_RHS_calc(q,dqdx,dqdy,dqdz,q_s,z_os,RHS_q,sgs_q3,wq_s_current)

!if (jt==1. .and. (.not. initu)) then
!RHS_qf=RHS_q
!end if

!call step_scalar(q,RHS_q,RHS_qf)
!end if

!end subroutine humidity_all_in_one



!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
subroutine calcbeta (scalar)
! This calculates the buoyancy term (beta_scal) to be added to the vertical
! momentum equation for temperature
! Authored by Vijayant Kumar
! Last updated April 14, 2004
implicit none
integer::i, j, k, jz_min
!real(kind=rprec),dimension(ld,ny,nz),intent(out)::beta_scal
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
real(kind=rprec),dimension(ld,ny,$lbz:nz),intent(in)::scalar
!DY Debugging1 by Di Yang
!DY real(kind=rprec),dimension(nz)::scalar_bar
real(kind=rprec),dimension($lbz:nz)::scalar_bar
!DY Debugging1 end here
real(kind=rprec)::g_hat,above, below

  !..Non-dimensionalize gravity
  g_hat=g*(z_i/(u_star**2))
  beta_scal=0._rprec

  !..Note Beta is stored on W nodes, but Theta is on UVP nodes
  !....We do not time-advance the ground nodes, so start at k=2
  ! VK: Inserted the averaging code inside this file itself
  ! rather than doing it in prof
  do k=$lbz,nz
    scalar_bar(k)=0.0    
    do j=1,ny
      do i=1,nx
        scalar_bar(k)=scalar_bar(k)+scalar(i,j,k)
      end do
   end do
   scalar_bar(k)=scalar_bar(k)/(nx*ny)
  end do

  !....We do not time-advance the ground nodes, so start at k=2
  !.. For the MPI case, this means that we start from jz=2 for
  !.. coord=0 and jz=1 otherwise... enable by an if statment

  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    jz_min = 2
  else
    jz_min = 1
  end if

  do k=jz_min,Nz-1
    do j=1,Ny
      do i=1,nx
        above=(scalar(i,j,k)-scalar_bar(k))/scalar_bar(k)
        below=(scalar(i,j,k-1)-scalar_bar(k-1))/scalar_bar(k-1)
        beta_scal(i,j,k)=g_hat*(above + below)/2._rprec
      end do
    end do     
  end do
  return
  
end subroutine calcbeta
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------



subroutine scalar_RHS_calc(scalar,dsdx,dsdy,dsdz,S_Surf,z_os,RHS,sgs_vert,surf_flux_current)

  use test_filtermodule

  implicit none

  $if ($MPI)
    $define $lbz 0
  $else
    $define $lbz 1
  $endif

  integer:: i, j, k, jz
  integer:: jz_min,ubc_jz
!  integer:: patch(nx,ny),patchnum(types)
  real:: surf_flux_current,crap2
!  real(kind=rprec),dimension(ld,ny,nz):: u,v,w - !No need as already invoked using sim_param
  real(kind=rprec),dimension(ld,ny,$lbz:nz):: dsdx,dsdy,dsdz
  real(kind=rprec),dimension(ld,ny,$lbz:nz):: RHS,temp
  real(kind=rprec),dimension(ld,ny,$lbz:nz):: scalar
  real(kind=rprec),dimension(ld,ny,$lbz:nz):: dtemp,sgs_vert
!  real,dimension(ld,ny,nz):: dtemp,s,txz,tyz,sgs_vert,Pr_,Nu_t
  real(kind=rprec),dimension(ld_big,ny2,$lbz:nz):: u_m,v_m,w_m,dsdx_m,dsdy_m,dsdz_m
  real(kind=rprec),dimension(ld_big,ny2,$lbz:nz):: RHS_m

  !cVK - changed the dimensions for RHS_m,u_m etc. to ld_big
  !cVK as otherwise it causes segmentation errors
  real(kind=rprec),dimension(nx,ny):: ustar_local,S_Surf,surf_flux,z_os
  real(kind=rprec),dimension(ld,ny):: scalar_node_1 ! test filtered and used for computing surface flux
!  real(kind=rprec),dimension(nx,ny):: psi_h,phi_h
!  real(kind=rprec),dimension(nx,ny):: psi_h,phi_h,ustar2
  real(kind=rprec),dimension (ptypes):: ustar,wt_s2
  character (64) :: fname_hetero

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for ocean LES with Stokes drift
real(kind=rprec),dimension($lbz:nz) :: ust
real (rprec) :: z
!DY End here
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for ocean LES with Stokes drift
!DY Start here
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if(OCEAN_FLAG .and. STOKES_FLAG) then
   do jz=$lbz,nz
      $if ($MPI)
      z = (coord*(nz-1) + jz - 0.5_rprec) * dz
      $else
      z=(real(jz)-0.5_rprec)*dz
      $endif
      ust(jz)=U_stokes*exp(-2._rprec*wavenm_w*z)
   enddo
endif
!+++++++++++++++
!DY End here
!+++++++++++++++


  if ((USE_MPI) .AND. (SCALAR_DEBUG)) then
    do jz=0,nz
      print *,'jz,RHS_T@entry in scal,RHS_m,coord',jz,RHS(4,4,jz),&
      RHS_m(4,4,jz),coord
    end do
  elseif (SCALAR_DEBUG) then
    do jz=0,nz
      print *,'jz,RHS_T@entry in scal,RHS_m',jz,RHS(4,4,jz),&
      RHS_m(4,4,jz)
    end do
  end if


  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    
    if (SCALAR_DEBUG) then
      do jz=1,nx
        print *,'jz,ustar_avg,coord',jz,ustar_avg(jz,4),coord
      end do
    elseif (SCALAR_DEBUG) then
      do jz=1,nx
        print *,'jz,ustar_avg',jz,ustar_avg(jz,4)
      end do
    end if

    ustar=0.

    if (patch_flag==1) then !if block 101

      if (spatial_flux_flag) then !SPATIAL_FLUX IF BLOCK! if block 102

        if (jt .eq. SCAL_init) then
          print *,'reading spatial flux data from file !!'
          open(unit=77,file=path//'spatial_flux.dat',status='unknown')
          do j=1,ny
            read(77,5168) (surf_flux(i,j),i=1,nx)
          end do
          surf_flux=surf_flux/u_star/T_scale !The readin surf_flux is dimensional - so non-dimensionalize
        else if(jt .GT. SCAL_init) then
          surf_flux=sgs_t3(1:nx,1:ny,1)
        end if !end for if loop for jt .eq. SCAL_init

        ustar_local=ustar_avg !set ustar as value computed in obukhov

        do j=1,ny
          do i=1,nx
            dsdz(i,j,1) =-phi_h(i,j)*surf_flux(i,j)/(ustar_local(i,j)*vonk*DZ/2._rprec) !set the gradient at the first point
          end do
        end do
!!VIJ!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! Need to work on this later to remove the above do loop
!          print *,'TEST1',dsdz(5:9,4:6,1)
!          dsdz(:,:,1) =-phi_h(:,:)*surf_flux(:,:)/(ustar_local(:,:)*vonk*DZ/2._rprec) !set the gradient at the first point
!          print *,'TEST2',dsdz(5:9,4:6,1)
!!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


5169    format(1400(E16.10))

      else    ! SPATIAL FLUX IF BLOCK CONTINUES !block 102 conts.

        do k=1,ptypes
          wt_s2(k)=(-1.)**(k+1)*surf_flux_current 
          !VK Creates patches of type ptypes with alternating signs of heat flux
        end do

!c This computes the average value of the scalar
!c S_surf refers to the value assigned to the surface according to
!c routine patches.f based on number of patches, and parameters in
!c dimen.h
! Also calculated is a local averaged value of u_star from the subgrid
! stress at the wall

        do j=1,ny
          do i=1,nx
            do k=1,ptypes
              if (patch(i,j)==k) then
!                ustar(patch(i,j))=ustar(patch(i,j))+(txz(i,j,1)**2+&
!                tyz(i,j,1)**2)**.25/patchnum(patch(i,j))
                ustar(patch(i,j))=ustar(patch(i,j))+ustar_avg(i,j)/patchnum(patch(i,j))
              end if
            end do
          end do
        end do


        ! distribute ustar value to patches
        do j=1,ny
          do i=1,nx
            do k=1,ptypes
              if (patch(i,j)==k) then
!                ustar2(i,j)=ustar(k)
                ustar_local(i,j)=ustar(k)
              end if
            end do
          end do
        end do

        ! Compute surface flux and dsdz at z=DZ/2
        do j=1,Ny
          do i=1,Nx

            ! lbc=1 is used for prescribing the surface flux while
            ! lbc=0 has been used to prescribe the temperature
            if (lbc==1.and.scalar(1,1,1)<2) then
              do k=1,ptypes
                if (patch(i,j)==k) then

                  ! The u_star is coming from dimen.h = Ug for coriolis and is not
                  ! the local u_star computed from stress at the surface.
                  surf_flux(i,j)=wt_s2(k)/T_scale/u_star

                end if
              end do
  
            else if (lbc==0.and.scalar(1,1,1)<2) then

              ustar_local=ustar_avg
              surf_flux(i,j)=(S_Surf(i,j)-scalar(i,j,1))*vonk*ustar_local(i,j)&
                             /(dlog(dz/(2._rprec*z_os(i,j)))-psi_h(i,j))

              !  equ. 4.28. in brutsaert
              ! Impose the gradient value at the first node.
              ! Note that ustar2 is a variable coming from sgs stresses at the wall,
              ! ustar2 has a single value for all the nodes in a single patch.
              ! phi_h is obtained from the routine obukhov.f, surf_flux is as computed above
              ! and vonk and dz are constants. Everything is in dimensionless form
            
	    end if

            ! Now we have the lowest dsdz on the UVP nodes all others on w nodes
            dsdz(i,j,1) =-phi_h(i,j)*surf_flux(i,j)/(ustar_local(i,j)*vonk*DZ/2._rprec)

          end do
        end do

      end if !end spatial heat flux flag !end if block 102
    
    
    elseif (remote_flag==1) then !if block 101 conts.
      
      if (lbc .eq. 0) then

        ustar_local=ustar_avg

        ! use variable scalar_node_1 to store the scalar field at node 1
        ! test filter this variable to bring the surface flux calculated later in a local
        ! formulation more in line with the average formulation as prescribed by
        ! Monin-Obukhov theory

        scalar_node_1=scalar(:,:,1)


        if (remote_flux_homog_flag .eq. 1) then
          if (mod(jt,100)==0) print *,'apply MO in avg sense: remote_flux homog_flag ON'
          crap2=sum(scalar_node_1(1:nx,1:ny))/real(nx*ny);
          scalar_node_1(1:nx,1:ny)=crap2
          ustar_local=sum(sqrt(u(1:nx,1:ny,1)**2+v(1:nx,1:ny,1)**2))/float(nx*ny)*vonK/&
                    (dlog(0.5_rprec*dz/exp(sum(dlog(zo(1:nx,1:ny)))/float(nx*ny)))-sum(psi_m(1:nx,1:ny))/float(nx*ny))
          if (mod(jt,1000)==0) then
            print *,'Testing remote flux homog calc:u*,theta_1,theta_sfc,psi_m,psi_h,phi_h,zo:(5:6,5)',&
              ustar_local(5:6,5),scalar_node_1(5:6,5),S_surf(5:6,5),psi_m(5:6,5),psi_h(5:6,5),phi_h(5:6,5),&
              zo(5:6,5)
          end if
        
	else

          ! Filtering only makes sense if the flux calculation is not based on mean values (remote_flux_homog_flga .ne. 1)
          call test_filter(scalar_node_1,G_test)

        end if


        do j=1,ny
          do i=1,nx
            surf_flux(i,j)=(S_Surf(i,j)-scalar_node_1(i,j))*vonk*ustar_local(i,j)&
        		   /(dlog(dz/(2._rprec*z_os(i,j)))-psi_h(i,j))

            dsdz(i,j,1) =-phi_h(i,j)*surf_flux(i,j)/(ustar_local(i,j)*vonk*DZ/2._rprec)
          end do
        end do

        if (jt==SCAL_init+1) then
          open(unit=56,file=path//'output/surf_flux_step1.txt',status="unknown",position="append")
          do j=1,ny
            write(56,5168) (surf_flux(i,j)*u_star*T_scale,i=1,nx)
          end do
          close(56)
        end if

 5168   format(1400(E14.5))
        
      else if ((lbc .eq. 1) .and. (spatial_flux_flag) ) then

        if (jt .eq. SCAL_init) then
          print *,'reading spatial flux data from file !!'
          open(unit=77,file='./spatial_flux.dat',status='unknown')
          do j=1,ny
            read(77,5168) (surf_flux(i,j),i=1,nx)
          end do
          if (remote_to_patch_flag) then
            print *,'mean surf flux BEFORE remote to patch: ',sum(surf_flux)/float(nx*ny) 
            print *,'Creating 2 patch spatial heat flux field from spatial flux data ...'
            call remote_to_patch(surf_flux,1) 
            print *,'mean surf flux AFTER remote to patch: ',sum(surf_flux)/float(nx*ny)
          else if (remote_homog_flag == 1) then
            print *,'Homogenizing the spatial flux field as remote_homog_flag == 1'
            surf_flux=sum(surf_flux)/float(nx*ny)
          end if
          print *,'surf flux read',surf_flux(2:4,2:4) 
          surf_flux=surf_flux/u_star/T_scale !The readin surf_flux is dimensional - so non-dimensionalize
        else if(jt .GT. SCAL_init) then
          surf_flux=sgs_t3(1:nx,1:ny,1)
        end if !end for if loop for jt .eq. SCAL_init

        ustar_local=ustar_avg !set ustar as value computed in obukhov

        do j=1,ny
          do i=1,nx
            dsdz(i,j,1) =-phi_h(i,j)*surf_flux(i,j)/(ustar_local(i,j)*vonk*DZ/2._rprec) !set the gradient at the first point
          end do
        end do

      end if

    end if !end if block 101


!!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Do the heterogeneity analysis
    if ((lbc .eq. 0) .AND. (jt_total .gt. SCAL_init)) then !ONLY FOR SURF BC = surface temperature !IF BLOCK 200

!      if (jt .eq. hetero_array_freqz) then !If Block 200.1 STARTS
      if (jt_total .eq. SCAL_init+1) then !If Block 200.1 STARTS
        print *,'Hetero vars: hetero_array_freqz,hetero_count_out = ',hetero_array_freqz,hetero_count_out
        if (mod(hetero_count_out,hetero_array_freqz) .NE. 0) then !BLOCK 200.1.1 STARTS
          print *,'Hetero count out not exactly divisible by hetero_array_freqz .. QUITTING' 
          print *,'Please change the values in scalars_module.f90'
          STOP  
          time_ind = 1
        end if  ! BLOCK 200.1.1 ENDS 
      end if ! BLOCK 200.1 ENDS
      
      if (mod(jt,hetero_array_freqz) .eq. 0) then !BLOCK 200.2 STARTS
        print *, 'Writing hetero fields out at jt = ',jt
        write (fname_hetero, '(a,i6.6,a)') path//'output/fields_3d/hetero_data',time_ind,'.bin'
        open(77,file=fname_hetero,form='unformatted',position='append')
        write(77) jt_total,real(scalar_node_1),real(ustar_local),&
                  real(psi_h),real(surf_flux),&
                  real(phi_h)
        close(77)
        if (mod(jt,hetero_count_out) .eq. 0) then !BLOCK 200.2.1 STARTS
          print *, 'Changing file name string counter at jt = ',jt
          time_ind = time_ind + 1
        end if !BLOCK 200.2.1 ENDS
      end if !BLOCK 200.2 ENDS
    end if !IF BLOCK 200 ENDS

!!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if ((USE_MPI) .AND. (SCALAR_DEBUG)) then
      do jz=0,nz
        print *,'jx,surf_flux,ustar_local,phi_h,dsdz@jz_1,coord',jz,surf_flux(jz,4),&
        ustar_local(jz,4),phi_h(jz,4),dsdz(jz,4,1),coord
      end do
    elseif (SCALAR_DEBUG) then
      do jz=0,nz
        print *,'jx,surf_flux,ustar_local,phi_h,dsdz@jz_1',jz,surf_flux(jz,4),&
        ustar_local(jz,4),phi_h(jz,4),dsdz(jz,4,1)
      end do
    end if

  end if !end for if (USE_MPI ...) block


  call dealias1(u,u_m)
  call dealias1(v,v_m)
  call dealias1(w,w_m)
  call dealias1(dsdx,dsdx_m)
  call dealias1(dsdy,dsdy_m)
  call dealias1(dsdz,dsdz_m)

  if ((USE_MPI) .AND. (SCALAR_DEBUG)) then
    do jz=0,nz
      print *,'jz,u_m,v_m,w_m,dsdx_m,dsdy_m,dsdz_m@b4 dealias@SCAL,coord',jz,u_m(4,4,jz),&
               v_m(4,4,jz),w_m(4,4,jz),dsdx_m(4,4,jz),dsdy_m(4,4,jz),dsdz_m(4,4,jz),coord
    end do
  elseif (SCALAR_DEBUG) then
    do jz=0,nz
      print *,'jz,RHS_T@b4 dealias@SCAL,RHS_m',jz,RHS(4,4,jz),&
              RHS_m(4,4,jz)
    end do
  end if

  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    jz_min = 2
  else
    jz_min = 1
  end if
 
!  if ((USE_MPI)) then

  ubc_jz = nz-1

!  else
!    ubc_jz=nz
!  end if


! Now compute the RHS term of the filtered scalar equation. 
! Note that this is the advection term with the scalar as 
! the diffusion term has been thrown away. This is done step 
! by step for each of the expanded arrays from dealias1 separately
! for the first node & last node AND the rest of the nodes.

! xxxxx ------Comments valid for MPI case only ---------XXXX
! For MPI case, all the nodes have fluid nodes (1:nz-1) except for
! near-the wall processes (2:nz-1 for w nodes and 1:nz-1 for uvp nodes)
! and the top nodes (1:nz)
! The following loop executes from 2:nz-1 for the near-wall process
! and 1:nz-1 for the other processes. Also, the subsequent 2 loops
! take care of the first node for the near-wall process (coord = 0)
! and the topmost node for the top process (coord = nproc-1).
! Note also that we need to make ghost data available for dTdz and
! w for the topmost node (jz=n) within each process and therefore,
! this synchronization (MPI_SENDRECV()) has been done in the subroutine
! theta_all_in_one ()
! xxxxx --------- MPI Comment block ends ------------------XXXX

!  do k=2,Nz-1
  do k=jz_min,nz-1
    do j=1,Ny2
      do i=1,Nx2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Modified by Di Yang for ocean LES with Stokes drift
!DY Start here
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY        RHS_m(i,j,k)=u_m(i,j,k)*dsdx_m(i,j,k)+v_m(i,j,k)*dsdy_m(i,j,k)&
!DY            +(w_m(i,j,k)*dsdz_m(i,j,k)+w_m(i,j,k+1)*dsdz_m(i,j,k+1))/2._rprec
         if(OCEAN_FLAG .and. STOKES_FLAG) then
            RHS_m(i,j,k)=(u_m(i,j,k)+ust(k))*dsdx_m(i,j,k)+v_m(i,j,k)*dsdy_m(i,j,k)&
                 +(w_m(i,j,k)*dsdz_m(i,j,k)+w_m(i,j,k+1)*dsdz_m(i,j,k+1))/2._rprec
         else
            RHS_m(i,j,k)=u_m(i,j,k)*dsdx_m(i,j,k)+v_m(i,j,k)*dsdy_m(i,j,k)&
                 +(w_m(i,j,k)*dsdz_m(i,j,k)+w_m(i,j,k+1)*dsdz_m(i,j,k+1))/2._rprec
         endif
!+++++++++++++++
!DY End here
!+++++++++++++++
      end do
    end do
  end do


  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    do j=1,Ny2
      do i=1,Nx2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Modified by Di Yang for ocean LES with Stokes drift
!DY Start here
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY        RHS_m(i,j,1)=u_m(i,j,1)*dsdx_m(i,j,1)+v_m(i,j,1)*dsdy_m(i,j,1)&
!DY                     +(0.5_rprec*w_m(i,j,2))*dsdz_m(i,j,2)
!DY!        RHS_m(i,j,Nz)=u_m(i,j,Nz)*dsdx_m(i,j,Nz)+v_m(i,j,Nz)*dsdy_m(i,j,Nz)
         if(OCEAN_FLAG .and. STOKES_FLAG) then
            RHS_m(i,j,1)=(u_m(i,j,1)+ust(1))*dsdx_m(i,j,1)+v_m(i,j,1)*dsdy_m(i,j,1)&
                 +(0.5_rprec*w_m(i,j,2))*dsdz_m(i,j,2)
         else
            RHS_m(i,j,1)=u_m(i,j,1)*dsdx_m(i,j,1)+v_m(i,j,1)*dsdy_m(i,j,1)&
                 +(0.5_rprec*w_m(i,j,2))*dsdz_m(i,j,2)
         endif
!+++++++++++++++
!DY End here
!+++++++++++++++
      end do
    end do
  end if
 
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
    do j=1,Ny2
      do i=1,Nx2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Modified by Di Yang for ocean LES with Stokes drift
!DY Start here
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY!        RHS_m(i,j,1)=u_m(i,j,1)*dsdx_m(i,j,1)+v_m(i,j,1)*dsdy_m(i,j,1)&
!DY!                    +(0.5_rprec*w_m(i,j,2))*dsdz_m(i,j,2)
!DY        RHS_m(i,j,Nz)=u_m(i,j,Nz)*dsdx_m(i,j,Nz)+v_m(i,j,Nz)*dsdy_m(i,j,Nz)
         if(OCEAN_FLAG .and. STOKES_FLAG) then
            RHS_m(i,j,Nz)=(u_m(i,j,Nz)+ust(Nz))*dsdx_m(i,j,Nz)+v_m(i,j,Nz)*dsdy_m(i,j,Nz)
         else
            RHS_m(i,j,Nz)=u_m(i,j,Nz)*dsdx_m(i,j,Nz)+v_m(i,j,Nz)*dsdy_m(i,j,Nz)
         endif
!+++++++++++++++
!DY End here
!+++++++++++++++
      end do
    end do
  end if

  if ((USE_MPI) .AND. (SCALAR_DEBUG)) then
    do jz=0,nz
      print *,'jz,RHS_T@b4 dealias@SCAL,RHS_m,coord',jz,RHS(4,4,jz),&
      RHS_m(4,4,jz),coord
    end do
  elseif (SCALAR_DEBUG) then
    do jz=0,nz
      print *,'jz,RHS_T@b4 dealias@SCAL,RHS_m',jz,RHS(4,4,jz),&
      RHS_m(4,4,jz)
    end do
  end if

  call dealias2(RHS,RHS_m)

  if ((USE_MPI) .AND. (SCALAR_DEBUG)) then
    do jz=0,nz
      print *,'jz,RHS_T@after dealias@SCAL,RHS_m,coord',jz,RHS(4,4,jz),&
      RHS_m(4,4,jz),coord
    end do
  elseif (SCALAR_DEBUG) then
    do jz=0,nz
      print *,'jz,RHS_T@after dealias@SCAL,RHS_m',jz,RHS(4,4,jz),&
      RHS_m(4,4,jz)
    end do
  end if

!c...Now building the SGS part of the RHS.
! Here the sgs_term for scalars is built up using Nu_t from sgs_stag_W.f
! and dividing it by the turbulent Prandtl # specified in dimen.h
!c....Note: Since we bring the Convective term to RHS its sign changes.
!c....Below "Temp" is used for SGS flux; its divergence is added to RHS
!VK.. Nu_t is on w nodes everywhere except at z=dz/2.. while
!VK dsdx is on uvp nodes.. so, interpolate Nu_t as we want temp to
!VK be on uvp nodes
! All this valid only till Pr_ is a constant..
! This will need a makeover once Pr_ becomes dynamic as well...

  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    ! at jz=1, both Nu_t and dsdx are on uvp nodes .. no need for interp
    do j=1,Ny
      do i=1,Nx
        temp(i,j,1)=(1./Pr_(i,j,1))*Nu_t(i,j,1)*dsdx(i,j,1)
      end do
    end do
  end if

!  do k=1,Nz
  do k=jz_min,ubc_jz
    do j=1,Ny
      do i=1,Nx
!        temp(i,j,k)=(1./Pr_(i,j,k))*Nu_t(i,j,k)*dsdx(i,j,k)
        temp(i,j,k)=(1./Pr_(i,j,k))*0.5_rprec*(Nu_t(i,j,k)+Nu_t(i,j,k+1))*dsdx(i,j,k)
      end do
    end do
  end do

  call DDX (dtemp, temp)  

  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    ! at jz=1, both Nu_t and dsdy are on uvp nodes .. no need for interp
    do j=1,Ny
      do i=1,Nx
        RHS(i,j,1) = (-1.*RHS(i,j,1) + dtemp(i,j,1))
        temp(i,j,1)=(1./Pr_(i,j,1))*Nu_t(i,j,1)*dsdy(i,j,1)
      end do
    end do
  end if

!  do k=1,ubc_jz
  do k=jz_min,ubc_jz
    ! Nu_t is on w nodes and dsdy is on uvp nodes for jz=2 to nz
    do j=1,Ny
      do i=1,Nx
        RHS(i,j,k) = (-1.*RHS(i,j,k) + dtemp(i,j,k))
        temp(i,j,k)=(1./Pr_(i,j,k))*0.5_rprec*(Nu_t(i,j,k)+Nu_t(i,j,k+1))*dsdy(i,j,k)
!        temp(i,j,k)=(1./Pr_(i,j,k))*Nu_t(i,j,k)*dsdy(i,j,k)
      end do
    end do
  end do

  if ((USE_MPI) .AND. (SCALAR_DEBUG)) then
    do jz=0,nz
      print *,'jz,RHS_T@after ddx@SCAL,coord',jz,RHS(4,4,jz),coord
    end do
  elseif (SCALAR_DEBUG) then
    do jz=0,nz
      print *,'jz,RHS_T@after ddx@SCAL',jz,RHS(4,4,jz)
    end do
  end if

  call DDY (dtemp, temp)   
  
  !...Use MO flux at wall for the scalar sgs term !
  ! Note that the total contribution to the scalar sgs term at
  ! the first node comes from the surface flux computed above from
  ! the specified heat flux, wt_s

  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    do j=1,Ny
      do i=1,Nx
        RHS(i,j,1) = RHS(i,j,1) + dtemp(i,j,1)
        temp(i,j,1) = -1.*surf_flux(i,j)
        sgs_vert(i,j,1) =surf_flux(i,j)
      end do
    end do
  end if


  ! Note sgs_vert is -1*temp because sgs_vert is modeled as -Nu_t*dsdz/Pr
  ! while temp is the same term but w/o the minus sign due to the additional
  ! minus outside the scalar fluctuation flux term in RHS
  ! need to run this loop nz due to the nature of the differenetiation in ddz_w

!  do k=2,Nz
!  do k=jz_min,ubc_jz
  do k=jz_min,nz
    do j=1,Ny
      do i=1,Nx
        RHS(i,j,k) = RHS(i,j,k) + dtemp(i,j,k)
!        temp(i,j,k)=(1./Pr_(i,j,k))*0.5*(Nu_t(i,j,k)+Nu_t(i,j,k-1))*dsdz(i,j,k)
        temp(i,j,k)=(1./Pr_(i,j,k))*Nu_t(i,j,k)*dsdz(i,j,k)
        sgs_vert(i,j,k)=-1.*temp(i,j,k)
      end do
    end do
  end do


  ! The SGS_z flux is on the W nodes, but DDZ_W will put it back on UVP nodes! 
  ! Also note that sgs_vert(i,j,k) influences the computations in 
  ! OBUKHOV.f and is not involved in any computations in this routine.
  ! sgs_t3(i,j,1) (<w'theta'> is used for computing wt at the surface in OBUKHOV)

  if ((USE_MPI) .AND. (SCALAR_DEBUG)) then
    do jz=0,nz
      print *,'jz,RHS_T@after ddy@SCAL,coord',jz,RHS(4,4,jz),coord
    end do
  elseif (SCALAR_DEBUG) then
    do jz=0,nz
      print *,'jz,RHS_T@after ddy@SCAL',jz,RHS(4,4,jz)
    end do
  end if

  call DDZ_w (dtemp, temp)

  do k=1,ubc_jz
    do j=1,Ny
      do i=1,Nx
        RHS(i,j,k) = RHS(i,j,k) + dtemp(i,j,k)
      end do
    end do
  end do

  if ((USE_MPI) .AND. (SCALAR_DEBUG)) then
    do jz=0,nz
      print *,'jz,RHS_T@after ddz@SCAL,coord',jz,RHS(4,4,jz),coord
    end do
  elseif (SCALAR_DEBUG) then
    do jz=0,nz
      print *,'jz,RHS_T@after ddz@SCAL',jz,RHS(4,4,jz)
    end do
  end if

!5167  format (110(1x,e14.5)) 

!  return 

end subroutine scalar_RHS_calc
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------




subroutine step_scalar(scalar,RHS_pre,RHS_post)
!subroutine step_scalar(scalar,RHS_T,RHS_Tf,wt_s_current)
implicit none
integer:: i,j,k
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
real(kind=rprec),dimension(ld,ny,$lbz:nz)::scalar, RHS_pre, RHS_post
!real(kind=rprec)::wt_s_current
!cVK - This routine moves the scalar field (scalar in this case)
!cVK - forward in time using the scalar from previous time step
!cVK - and the RHS terms from the previous two time steps 
!cVK - using second order Adams-Bashforth scheme

! Note that in the last staments within this file, we set the value
! of scalar at the topmost node based on prescribed bc (inv_strength)
! and so, we could save some computation by only performing
! the scalar computation till Nz-1 global node...

do k=1,nz-1
!do k=1,nz
      do j=1,ny
             do i=1,nx
                 scalar(i,j,k)= scalar(i,j,k)+dt*(1.5_rprec*RHS_pre(i,j,k)-0.5_rprec*RHS_post(i,j,k))
             end do
      end do
end do     

!VK Note that this subroutine was designed to be working with a set of scalars (incl.
!VK temperature and humidity and so, these boundary conditions as given below should
!VK be interpreted in the right context and not just for temperature
!VK For example, the first if block refers to a condition with humidity while the
!VK second and third statements are focussed to temperature

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
! if MPI - then clicks and is used for the process dealing wih the top nodes
! else in general is used for temp bc
! Note that L_z*nproc is the total vertical extent of the domain for the MPI and
! non-MPI cases ..(with MPI in place, we can not use L_z>z_i anymore)
     if ((L_z*nproc)>z_i) then ! for temperature and non-neutral case
!     if (scalar(2,2,2)<2._rprecand.L_z>z_i) then ! for temperature and non-neutral case
!         print *,'setting scalar value at topmost global node for coord = ',coord
!         scalar(:,:,Nz)=scalar(:,:,Nz-1)+inv_strength/T_scale*z_i*dz !ubc 
         scalar(:,:,Nz)=scalar(:,:,Nz-1)+dTdz_top/T_scale*z_i*dz !ubc 
! inv_strength - refers to the slope of the inversion ~ 0.003K/Km for temperature
     else ! for everything else - neutral and passive scalars (may be modified depending on need)
         scalar(:,:,Nz)=scalar(:,:,Nz-1)
     end if
end if

!return
end subroutine step_scalar



!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
subroutine obukhov  

use types,only:rprec
!use param
use sim_param,only:u,v,theta,path
use bottombc !Includes patches subroutine, phi_m,psi_m,phi_h,psi_h,ustar_avg
use test_filtermodule
!use scalars_module
implicit none
!integer,intent(in)::jt
integer:: jx,jy

real(kind=rprec), dimension(ld,ny):: wt_avg,wq_avg,theta_avg,u1,v1
real(kind=rprec), dimension(nx,ny):: x,zeta ! wstar, L already defined above
real(kind=rprec) g_,wt_,wq_,ustar_,theta_,L_,zo_,u_avg(nx,ny),fr,wstar_avg

real(kind=rprec),save:: obuk_L,obuk_ustar,obuk_phi_m,obuk_phi_h,obuk_psi_m   
real(kind=rprec),save:: obuk_wt_sfc,obuk_psi_h,obuk_zo,obuk_wstar   

  if (jt .LT. SCAL_init .OR. (.NOT. S_FLAG)) then  
    if (.not. initsc) psi_m=0._rprec !psi_m is being read in initial.f90
    
    phi_m=1._rprec  
    psi_h=0._rprec  
    phi_h=1._rprec 
    ustar_avg=sum(sqrt(u(1:nx,1:ny,1)**2+v(1:nx,1:ny,1)**2))/float(nx*ny)*vonK/&
              (dlog(0.5_rprec*dz/zo)-sum(psi_m(1:nx,1:ny))/float(nx*ny))
    L=0._rprec
    wstar=0._rprec
    return
  end if  

  ! Do the following only for jt .eq. SCAL_init
  ! This becomes curcial for the supercomputing runs as we need to break the
  ! simulation into smaller chunks and thereby need an accurate way to continue
  ! the simulation from the vel_sc.out file..
  ! therefore, added sgs_t3(:,:,1) i.e. the surface flux to the list of output variables
  ! in io.f90
  if (jt .EQ. SCAL_init) then
    obuk_L=0._rprec;obuk_ustar=0._rprec;obuk_wstar=0._rprec;obuk_phi_m=0._rprec
    obuk_phi_h=0._rprec;obuk_psi_m=0._rprec;obuk_psi_h=0._rprec;obuk_zo=0._rprec
    if (.not. initsc) then
      sgs_t3(:,:,1)=wt_s/u_star/T_scale
    end if
  end if

  !  nondimensionalize g
  g_=g/(u_star**2/z_i)

  theta_avg=theta(:,:,1) 
  wt_avg=sgs_t3(:,:,1) ! We need only the surface flux - defined by sgs
  zo_=exp(sum(dlog(zo(1:nx,1:ny)))/float(nx*ny))

  ! averages over x-y plane @ z = 1
  wt_=sum(sgs_t3(1:nx,1:ny,1))/float(nx*ny)
  ustar_=sum(sqrt(u(1:nx,1:ny,1)**2+v(1:nx,1:ny,1)**2))/float(nx*ny)*vonK/&
            (dlog(0.5_rprec*dz/zo_)-sum(psi_m(1:nx,1:ny))/float(nx*ny))
  theta_=sum(theta_avg(1:nx,1:ny))/float(nx*ny)


  if ((patch_flag==1 .and. num_patch==1) .OR. (OB_homog_flag)) then  
    do jx=1,nx
      do jy=1,ny
        wt_avg(jx,jy)=wt_
        ustar_avg(jx,jy)=ustar_
        theta_avg(jx,jy)=theta_
      end do
    end do
  else
    u1=u(:,:,1)
    v1=v(:,:,1)
    call test_filter(u1,G_test)
    call test_filter(v1,G_test)
    do jx=1,nx
      do jy=1,ny
        u_avg(jx,jy)=sqrt(u1(jx,jy)**2+v1(jx,jy)**2)
      end do
    end do
    ustar_avg(1:nx,:)=u_avg(:,:)*vonK/(dlog(0.5_rprec*dz/zo(:,:))-psi_m(:,:))
  end if


  ! Compute Obukhov Length  
  do jx=1,ny
    do jy=1,nx
  
      L(jx,jy)=-ustar_avg(jx,jy)**3/(vonk*g_/theta_avg(jx,jy)*wt_avg(jx,jy))

      ! w_star is defined as [(g/<T_0>)*<w'T'>*z_i]**(1/3) (refer to 
      ! Nieuwstadt et al., Turbulent Shear flows, 1991)
      ! Therefore, for our case, where we are computing the non-dimensional w_star,
      ! the formula transforms to [(g_nd/<T_0_nd>)*<w'T'>_nd]**(1/3) where the suffix
      ! _nd refers to being non-dimensionalized using Ug (coriolis velocity), T_scale (300K)
      ! and z_i
      wstar(jx,jy)=sign((g_/theta_avg(jx,jy)*abs(wt_avg(jx,jy)))**(1./3.),wt_avg(jx,jy))

      !    wstar(jx,jy)=sign((g_/theta_avg(jx,jy)*abs(wt_avg(jx,jy))*z_i)**(1./3.),wt_avg(jx,jy))
      ! The above is the earlier wrong formula where there has been this additional z_i which for
      ! the post-processing means division by 10 as the usual value of z_i=1000 which with cube root on
      ! w_star just becomes a multiplication by 10. So, in case you feel that the value of w_star is quite
      ! high and the data is dated before Dec. 7th, 2004, please make sure to divide by 10
      !c  for unstable conditions
      if ((L(jx,jy)<0._rprec) .and. (wt_avg(jx,jy) .ne. 0._rprec)) then
        x(jx,jy)=(1._rprec-16._rprec*dz/2._rprec/L(jx,jy))**.25_rprec
        psi_m(jx,jy)=2._rprec*dlog((1.+x(jx,jy))/2._rprec)+&
        dlog((1._rprec+x(jx,jy)**2)/2._rprec)-2._rprec*datan(x(jx,jy))+pi/2._rprec
        psi_h(jx,jy)=2._rprec*dlog((1._rprec+x(jx,jy)**2)/2._rprec)
        phi_m(jx,jy)=x(jx,jy)**(-1)
        phi_h(jx,jy)=x(jx,jy)**(-2)
      else if ((L(jx,jy)>0._rprec).and.(wt_avg(jx,jy).ne. 0._rprec)) then

        ! Implementing new formulations for phi and psi for stable case
        ! using Cheng & Brutsaert (2004): source - Brutsaert's book from
        ! Marc's Hydrology course
        ! the new relations are from the GABLS study
        zeta(jx,jy)=0.5_rprec*dz/L(jx,jy)

        ! %%%%%%%%%%%%%%%%%%%%%%%% Formulations used in GABLS (Good for z/L < 1)%%%%%%%%%%%%%%%%
        ! %%%%%%%%%%%%%%%%%%%%%%%% Also apeearing in Hogstorm, BLM, 1988)%%%%%%%%%%%%%%%%
        if ((jan_diurnal_run) .OR. (GABLS_diurnal_test)) then !USE MO functions for z/L < 1

          ! %%%%%%%%%%%%%%%%%%%%%%%% Formulations used in Cheng & Brutsaert (BLM,2005) %%%%%%%%%%%%
          ! X=z/L where z in present case is dz/2 and L is the Obukhov length
          ! The equations: phi_m(X) = 1+a*[(X+(X^b)*{(1+X^b)^(-1+1/b)})/(X+{(1+X^b)^(1/b)})];a=6.1;b=2.5
          ! The equations: phi_h(X) = 1+c*[(X+(X^d)*{(1+X^d)^(-1+1/d)})/(X+{(1+X^d)^(1/d)})];c=5.3;d=1.1
          ! The equations: psi_m(X) = -a*ln[X+{(1+X^b)^(1/b)}];a=6.1;b=2.5
          ! The equations: psi_h(X) = -c*ln[X+{(1+X^d)^(1/d)}];c=5.3;d=1.1
          ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          phi_m(jx,jy)=1._rprec+6.1_rprec*(zeta(jx,jy)+zeta(jx,jy)**2.5_rprec*&
                       ((1._rprec+zeta(jx,jy)**2.5_rprec)**(-1._rprec+1._rprec/2.5_rprec)))/(zeta(jx,jy)+&
                       ((1._rprec+zeta(jx,jy)**2.5_rprec)**(1._rprec/2.5_rprec)))

          psi_m(jx,jy)=-1._rprec*6.1_rprec*dlog(zeta(jx,jy)+((1._rprec+&
                        zeta(jx,jy)**2.5_rprec)**(1._rprec/2.5_rprec)))
            
          phi_h(jx,jy)=1._rprec+5.3_rprec*(zeta(jx,jy)+zeta(jx,jy)**1.1_rprec*&
                       ((1._rprec+zeta(jx,jy)**1.1_rprec)**(-1._rprec+1._rprec/1.1_rprec)))/(zeta(jx,jy)+&
                       ((1._rprec+zeta(jx,jy)**1.1_rprec)**(1._rprec/1.1_rprec)))

          psi_h(jx,jy)=-1._rprec*5.3_rprec*dlog(zeta(jx,jy)+((1._rprec+&
                         zeta(jx,jy)**1.1_rprec)**(1._rprec/1.1_rprec)))

          ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        else

          ! %%%%%%%%%%%%%%%%%%%%%%%% Formulations used in Brutsaert (1982) %%%%%%%%%%%%%%%%%%%%
          !  phi_m(jx,jy)=1._rprec+5.0_rprec*zeta(jx,jy)
          !  phi_h(jx,jy)=1._rprec+5.0_rprec*zeta(jx,jy)
          !  psi_m(jx,jy)=-1._rprec*5.0_rprec*zeta(jx,jy)
          !  psi_h(jx,jy)=-1._rprec*5.0_rprec*zeta(jx,jy)

          phi_m(jx,jy)=1._rprec+4.8_rprec*zeta(jx,jy)
          phi_h(jx,jy)=1._rprec+7.8_rprec*zeta(jx,jy)
          psi_m(jx,jy)=-1._rprec*4.8_rprec*zeta(jx,jy)
          psi_h(jx,jy)=-1._rprec*7.8_rprec*zeta(jx,jy)

        end if ! end for the jan_diurnal_run if block     
      
      else
         psi_m(jx,jy)=0._rprec
         psi_h(jx,jy)=0._rprec
         phi_m(jx,jy)=1._rprec
         phi_h(jx,jy)=1._rprec
      
      end if ! (Loop5 ends)

    end do
  end do

  L_=-(ustar_**3)/(vonk*(g_/theta_)*wt_)
  wstar_avg=sign((g_/theta_*abs(wt_))**(1./3.),wt_)

  write (6,7780) L_*z_i,ustar_*u_star,theta_*T_scale,(sum(T_s)/float(nx*ny))*T_scale,dz/2/L_,wt_*u_star*T_scale
  
7780 format ('L(m),ustar(m/s),theta_1(K),T_s(K),z/L,wt_s(Km/s):',(5(1x,E12.6),1x,E12.6))


  !C! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !-------------------- OUTPUT ------------------------------
  ! Output the heat flux time series to a file to be used later

  open (unit=47,file=path//'output/WT_sfc_tseries.out',status="unknown",position="append")
  write(47,5168) (jt_total+1)*dt,wt_
  close(47)
  !C! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if (((initsc) .AND. ((jt_total+1) .ge. SCAL_init)) .OR. &
     ((.not. initsc) .AND. ((jt_total+1) .ge. SCAL_init+1))) then
  
    fr=1._rprec/c_count
    obuk_L=obuk_L+fr*L_
    obuk_ustar=obuk_ustar+fr*ustar_
    obuk_wstar=obuk_wstar+fr*wstar_avg
    obuk_phi_m=obuk_phi_m+fr*sum(phi_m(:,:))/float(nx*ny)
    obuk_psi_m=obuk_psi_m+fr*sum(psi_m(:,:))/float(nx*ny)
    obuk_phi_h=obuk_phi_h+fr*sum(phi_h(:,:))/float(nx*ny)
    obuk_psi_h=obuk_psi_h+fr*sum(psi_h(:,:))/float(nx*ny)
    obuk_zo=obuk_zo+fr*zo_
    obuk_wt_sfc=obuk_wt_sfc+fr*wt_

    if (mod(jt,c_count)==0) then
      open (unit=47,file=path//'output/mo.out',status="unknown",position="append")
      write(47,5168) (jt_total+1)*dt,obuk_L,obuk_ustar,obuk_wstar,obuk_phi_m,obuk_psi_m,&
                     obuk_phi_h,obuk_psi_h,obuk_zo,obuk_wt_sfc
      close(47)

      obuk_L=0._rprec;obuk_ustar=0._rprec;obuk_wstar=0._rprec;obuk_phi_m=0._rprec
      obuk_phi_h=0._rprec;obuk_psi_m=0._rprec;obuk_psi_h=0._rprec;obuk_zo=0._rprec
      obuk_wt_sfc=0._rprec;
    end if
  end if

5168     format(E14.5,9(1x,E14.5))

  return

end subroutine obukhov 

!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------


















! Main routine for pollen time advance
! Chamecki - 08/01/2006

subroutine pollen_all_in_one

  use test_filtermodule

  implicit none

  integer:: i, j
  real(kind=rprec):: sfc_flux_current,dummy_t
  integer :: jz,ios,counter

$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

  real(kind=rprec),dimension(ld,ny,$lbz:nz)::dPCondx,dPCondy


  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then

    ! This smight be useful in the near future
    if (jan_diurnal_run) then !perform time variable heat flux
!      if (jt .eq. SCAL_init) then
!        counter=1
!        open(unit=56,file=path//'heatflux_diurnal_run.txt')
!        do
!          read(56,*,iostat=ios) wt_s_vector(counter,1)
!          if (ios/=0) exit
!          counter=counter+1
!        end do
!        print *,'length of wt_s_vector = ',size(wt_s_vector)
!      end if
!
!      ! note that jt_total is always jt-1 as it is updated only at the end of each step
!      ! and hence the formula below reflects that !!
!      if (GABLS_diurnal_test) then !read a continuous WT tseries for GABLS diurnal
!        wt_s_current=wt_s_vector(jt_total+1,1)
!      else ! for the actual jan diurnal run 
!        wt_s_current=wt_s_vector(floor(((jt_total)*dt*z_i/u_star)/300._rprec)+1,1)
!      end if
!
!      print '(A,F12.6,A,I3.3)','wt_s = ',wt_s_current&
!           ,' at time (hour) = ',floor(((jt_total+1)*dt*z_i/u_star)/3600)
    else    
      sfc_flux_current=sfc_flux
    end if
  end if


  ! Calculate horizontal derivatives
  call filt_da(PCon,dPCondx,dPCondy)
  
  
  ! Calculate vertical derivative
  call ddz_uv(dPCondz,PCon)  

  
$if ($MPI)
! Need to synchronize w and dPCondz across the processors for uvp node
! based computation in scalar_RHS_calc (calculation of RHS_m)
  call mpi_sendrecv (w(1, 1, 1), ld*ny, MPI_RPREC, down, tag_counter+1,  &
                     w(1, 1, nz), ld*ny, MPI_RPREC, up, tag_counter+1,   &
                     comm, status, ierr)   
  call mpi_sendrecv (w(1, 1, nz-1), ld*ny, MPI_RPREC, up, tag_counter+2,  &
                     w(1, 1, 0), ld*ny, MPI_RPREC, down, tag_counter+2,   &
                     comm, status, ierr)
  call mpi_sendrecv (dPCondz(1, 1, 1), ld*ny, MPI_RPREC, down, tag_counter+3,  &
                     dPCondz(1, 1, nz), ld*ny, MPI_RPREC, up, tag_counter+3,   &
                     comm, status, ierr)   
  call mpi_sendrecv (dPCondz(1, 1, nz-1), ld*ny, MPI_RPREC, up, tag_counter+4,  &
                     dPCondz(1, 1, 0), ld*ny, MPI_RPREC, down, tag_counter+4,   &
                     comm, status, ierr)

! Also need to synchronize Nu_t across the processors for 
! computation in scalar_RHS_calc (calculation of RHS_m)
  call mpi_sendrecv (Nu_t(1, 1, 1), ld*ny, MPI_RPREC, down, tag_counter+5,  &
                     Nu_t(1, 1, nz), ld*ny, MPI_RPREC, up, tag_counter+5,   &
                     comm, status, ierr)   
  call mpi_sendrecv (Nu_t(1, 1, nz-1), ld*ny, MPI_RPREC, up, tag_counter+6,  &
                     Nu_t(1, 1, 0), ld*ny, MPI_RPREC, down, tag_counter+6,   &
                     comm, status, ierr)
$endif



  ! Main loop to advance PCon in time
  IF (PCon_FLAG) THEN
  

    RHS_PConf=RHS_PCon


!    ! Perform test filtering of PCon_s for calculation of surf fluxes
!    if ((jt_total .eq. PCon_init) .and. (lbc .eq. 0) .and. (remote_homog_flag .eq. 0)) then
!      print *,'T_s b4 filtering',sqrt(sum((T_s-sum(T_s)/float(nx*ny))**2))/float(nx*ny)
!      T_s_filtered(1:nx,1:ny)=T_s
!      call test_filter(T_s_filtered,G_test)
!      T_s=T_s_filtered(1:nx,1:ny)
!      print *,'T_s after filtering',sqrt(sum((T_s-sum(T_s)/float(nx*ny))**2))/float(nx*ny)
!    end if
!    if ((jt_total .eq. nsteps-1) .and. (lbc .eq. 0) .and. (remote_homog_flag .eq. 0)) then
!      print *,'T_s after filtering',sqrt(sum((T_s-sum(T_s)/float(nx*ny))**2))/float(nx*ny)
!    end if


    ! Calculate RHS of equation
    call pollen_RHS_calc(PCon,dPCondx,dPCondy,dPCondz,PCon_s,zo_PCon,RHS_PCon,sgs_PCon3,sfc_flux_current)
    
    


    if ((jt .eq. PCon_init) .or. (pointsource .and. (jt .eq. ini_src))) then
      RHS_PConf=RHS_PCon
    end if

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for buoyancy term added to vertical momentum equation
!DY Start here
!!$    if(active_pcon) then
!!$       call calcbeta_pcon(PCon)
!!$    endif
!DY End here
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    call step_pollen(PCon,RHS_PCon,RHS_PConf)


!    ! Test to verify standard deviation of wind direction at the source point
!    OPEN(1,FILE=path//'output/source_wind.out',STATUS="unknown",POSITION="append")
!    WRITE(1,'(I8,3F12.6)')jt,u(xps,yps,zps),v(xps,yps,zps),w(xps,yps,zps)
!    CLOSE(1)


    ! Check mass balance
    
    dummy_t=(released-airborne-deposited-gone)

    ! Grains released (point source)
    if (pointsource .and. (jt_total >= ini_src)) then
      released=((jt_total-ini_src+1)*dt*Q_src*dx*dy*dz)
    end if
    if (pointsource .and. (jt_total >= end_src)) then
      released=((end_src-ini_src)*dt*Q_src*dx*dy*dz)
    end if
    
!    ! Grains airborne and deposited in the valid domain
!    airborne=(SUM(PCon(2:nx-x_relaxp-1,:,1:nz-1))*dx*dy*dz)
!    airborne=airborne+(SUM(PCon(1,:,1:nz-1))+SUM(PCon(nx-x_relaxp,:,1:nz-1)))*dx*dy*dz/2._rprec
!    deposited=(SUM(deposition(2:nx-x_relaxp-1,:))*dx*dy)
!    deposited=deposited+(SUM(deposition(1,:))+SUM(deposition(nx-x_relaxp,:)))*dx*dy/2._rprec
!
!    ! Integrated flux at x=(nx-x_relaxp)
!    gone=gone+(SUM(PCon(nx-x_relaxp,:,:)*u(nx-x_relaxp,:,:))*dt*dy*dz)

    
    ! Grains airborne and deposited in the valid domain
    airborne=(SUM(PCon(2:nx-x_relaxp-5,:,1:nz-1))*dx*dy*dz)
    airborne=airborne+(SUM(PCon(1,:,1:nz-1))+SUM(PCon(nx-x_relaxp-4,:,1:nz-1)))*dx*dy*dz/2._rprec
    deposited=(SUM(deposition(2:nx-x_relaxp-5,:))*dx*dy)
    deposited=deposited+(SUM(deposition(1,:))+SUM(deposition(nx-x_relaxp-4,:)))*dx*dy/2._rprec

    ! Integrated flux at x=(nx-x_relaxp)
    flux_out_prev=flux_out
    flux_out=(SUM(PCon(nx-x_relaxp-4,:,1:nz-1)*u(nx-x_relaxp-4,:,1:nz-1))*dt*dy*dz)

    flux_out=flux_out-(SUM((1._rprec/Sc)* &
                           0.5_rprec*(Nu_t(nx-x_relaxp-4,:,2:nz-1)+Nu_t(nx-x_relaxp-4,:,3:nz))* &
			   dPCondx(nx-x_relaxp-4,:,2:nz-1))*dt*dy*dz)

    flux_out=flux_out-(SUM((1._rprec/Sc)*Nu_t(nx-x_relaxp-4,:,1)*dPCondx(nx-x_relaxp-4,:,1))*dt*dy*dz)


    flux_out=flux_out-(SUM(PCon(1,:,1:nz-1)*u(1,:,1:nz-1))*dt*dy*dz)

    flux_out=flux_out+(SUM((1._rprec/Sc)* &
                           0.5_rprec*(Nu_t(1,:,2:nz-1)+Nu_t(1,:,3:nz))* &
			   dPCondx(1,:,2:nz-1))*dt*dy*dz)

    flux_out=flux_out+(SUM((1._rprec/Sc)*Nu_t(1,:,1)*dPCondx(1,:,1))*dt*dy*dz)

    gone=gone+1.5_rprec*flux_out-0.5_rprec*flux_out_prev



!    ! Grains airborne and deposited in the valid domain
!    airborne=(SUM(PCon(1:nx-x_relaxp-3,:,1:nz-1))*dx*dy*dz)
!    deposited=(SUM(deposition(1:nx-x_relaxp-3,:))*dx*dy)

!    ! Integrated flux at x=(nx-x_relaxp)
!    gone=gone+(SUM(PCon(nx-x_relaxp-2,:,:)*u(nx-x_relaxp-2,:,:))*dt*dy*dz)




! Use this if there is no inflow condition
    airborne=(SUM(PCon(1:nx,:,1:nz-1))*dx*dy*dz)
    deposited=(SUM(deposition(1:nx,:))*dx*dy)
    gone=0._rprec
        
    
    OPEN(1,FILE=path//'output/pollen_balance.out',STATUS="unknown",POSITION="append")
    WRITE(1,'(I8,6G18.8)')jt_total,(released*PCon_scale*z_i*z_i*z_i), &                                 ! Pollen grains released (#)
                             (airborne*PCon_scale*z_i*z_i*z_i), &                                 ! Pollen grains airborne (#)
			     (deposited*PCon_scale*z_i*z_i*z_i), &                                ! Pollen grains deposited (#)
			     (gone*PCon_scale*z_i*z_i*z_i), &                                     ! Pollen grains out of domain (#)
			     ((released-airborne-deposited-gone)*PCon_scale*z_i*z_i*z_i), &
                             ((released-airborne-deposited-gone)-dummy_t*PCon_scale*z_i*z_i*z_i)

!    WRITE(1,'(I8,26F18.8)')jt,released,airborne,deposited,gone,(released-airborne-deposited-gone), &
!                             SUM(PCon(1:nx,:,1)),SUM(PCon(1:nx,:,2)), &
!                             SUM(PCon(1:nx,:,3)),SUM(PCon(1:nx,:,4)), &
!                             SUM(PCon(1:nx,:,5)),SUM(PCon(1:nx,:,6)), &
!                             SUM(PCon(1:nx,:,28)),SUM(PCon(1:nx,:,29)), &
!                             SUM(PCon(1:nx,:,30)),SUM(PCon(1:nx,:,31)), &
!                             SUM(PCon(1:nx,:,32)),SUM(PCon(1:nx,:,33)), &
!			     SUM(dPCondz(1:nx,:,1)),SUM(dPCondz(1:nx,:,2)),&
!			     (SUM(RHS_PCon(1:nx,:,1))*dx*dy*dt), &
!			     (SUM(RHS_PCon(1:nx,:,2))*dx*dy*dt), &
!			     (SUM(RHS_PCon(1:nx,:,3))*dx*dy*dt),&
!			     (SUM(RHS_PCon(1:nx,:,nz))*dx*dy*dt), &
!			     (SUM(RHS_PCon(1:nx,:,nz-1))*dx*dy*dt), &
!			     (SUM(RHS_PCon(1:nx,:,nz-2))*dx*dy*dt), &
!			     (released-airborne-deposited-gone)-dummy_t
			     
!    PRINT*,'grains released ',jt,released
!    PRINT*,'airborne grains ',jt,airborne
!    PRINT*,'deposited grains',jt,deposited
!    PRINT*,'grains gone     ',jt,gone
!    PRINT*,'lost grains     ',jt,released-airborne-deposited-gone
    CLOSE(1)

    ! Output timeseries at fixed locations for statistics
    OPEN(1,FILE=path//'output/time_series.out',STATUS="unknown",POSITION="append")
    WRITE(1,'(I8,7G18.8)')jt_total,              &       ! time
                          PCon(xps,yps,zps),     &       ! source point
                          PCon(xps+8,yps,zps),   &       ! 
                          PCon(xps+8,yps-4,zps), &       ! 
                          PCon(xps+8,yps+4,zps), &       ! 
                          PCon(xps+8,yps,zps-4), &       ! 
                          PCon(xps+8,yps,zps+4), &       ! 
                          PCon(xps+16,yps,zps)           ! 


    CLOSE(1)

 
  end if
  

$if ($MPI)
  call mpi_sendrecv (PCon(1, 1, 1), ld*ny, MPI_RPREC, down, tag_counter+7,  &
                     PCon(1, 1, nz), ld*ny, MPI_RPREC, up, tag_counter+7,   &
                     comm, status, ierr)   
  call mpi_sendrecv (PCon(1, 1, nz-1), ld*ny, MPI_RPREC, up, tag_counter+8,  &
                     PCon(1, 1, 0), ld*ny, MPI_RPREC, down, tag_counter+8,   &
                     comm, status, ierr)
$endif


end subroutine pollen_all_in_one




!
!  Sorry, but it was too much work to use the scalar_RHS_calc...
!

subroutine pollen_RHS_calc(scalar,dsdx,dsdy,dsdz,S_Surf,z_os,RHS,sgs_vert,surf_flux_current)

  use test_filtermodule

  implicit none

  $if ($MPI)
    $define $lbz 0
  $else
    $define $lbz 1
  $endif

  integer:: i, j, k, jz
  integer:: jz_min,ubc_jz
  real(kind=rprec):: surf_flux_current,term1,dist_src
  real(kind=rprec),dimension(ld,ny,$lbz:nz):: dsdx,dsdy,dsdz
  real(kind=rprec),dimension(ld,ny,$lbz:nz):: RHS,temp
  real(kind=rprec),dimension(ld,ny,$lbz:nz):: scalar
  real(kind=rprec),dimension(ld,ny,$lbz:nz):: dtemp,sgs_vert
  real(kind=rprec),dimension(ld_big,ny2,$lbz:nz):: u_m,v_m,w_m,dsdx_m,dsdy_m,dsdz_m
  real(kind=rprec),dimension(ld_big,ny2,$lbz:nz):: RHS_m

  real(kind=rprec),dimension(nx,ny):: ustar_local,S_Surf,z_os
  real(kind=rprec),dimension(nx,ny):: surf_flux_m   ! previous surface flux for correct calculation of deposition
  real(kind=rprec),dimension(ld,ny):: scalar_node_1 ! test filtered and used for computing surface flux

  real(kind=rprec),dimension (ptypes):: ustar,wt_s2
  character (64) :: fname_hetero

  integer:: ncpu_source, ncpu_source_bl, ncpu_source_ab
  integer:: zps_local, zps_bl_local, zps_ab_local

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for ocean LES with Stokes drift
real(kind=rprec),dimension($lbz:nz) :: ust
real (rprec) :: z
!DY End here
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for ocean LES with Stokes drift
!DY Start here
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if(OCEAN_FLAG .and. STOKES_FLAG) then
   do jz=$lbz,nz
      $if ($MPI)
      z = (coord*(nz-1) + jz - 0.5_rprec) * dz
      $else
      z=(real(jz)-0.5_rprec)*dz
      $endif
      ust(jz)=U_stokes*exp(-2._rprec*wavenm_w*z)
   enddo
endif
!+++++++++++++++
!DY End here
!+++++++++++++++

  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    
    ustar=0.



!    do k=1,ptypes
!      wt_s2(k)=(-1.)**(k+1)*surf_flux_current 
!      !VK Creates patches of type ptypes with alternating signs of heat flux
!    end do

    do k=1,ptypes
      wt_s2(k)=surf_flux_current 
    end do

!c This computes the average value of the scalar
!c S_surf refers to the value assigned to the surface according to
!c routine patches.f based on number of patches, and parameters in
!c dimen.h
! Also calculated is a local averaged value of u_star from the subgrid
! stress at the wall

    do j=1,ny
      do i=1,nx
        do k=1,ptypes
          if (patch(i,j)==k) then
!            ustar(patch(i,j))=ustar(patch(i,j))+(txz(i,j,1)**2+&
!            tyz(i,j,1)**2)**.25/patchnum(patch(i,j))
            ustar(patch(i,j))=ustar(patch(i,j))+ustar_avg(i,j)/patchnum(patch(i,j))
          end if
        end do
      end do
    end do


    ! distribute ustar value to patches
    do j=1,ny
      do i=1,nx
        do k=1,ptypes
          if (patch(i,j)==k) then
!            ustar2(i,j)=ustar(k)
            ustar_local(i,j)=ustar(k)
          end if
        end do
      end do
    end do
    
    ! Stored previous surface flux to integrate deposition
    surf_flux_m=P_surf_flux

    ! Compute surface flux and dsdz at z=DZ/2
    do j=1,Ny
      do i=1,Nx

        ! lbc=1 is used for prescribing the surface flux while
        ! lbc=0 has been used to prescribe the temperature
        if (lbcp==1) then
          do k=1,ptypes
            if (patch(i,j)==k) then
!              P_surf_flux(i,j)=wt_s2(k)/PCon_scale/u_star
              P_surf_flux(i,j)=wt_s2(k)
              ! The u_star above is coming from dimen.h = Ug for coriolis and is not
              ! the local u_star computed from stress at the surface.
            end if
!            print *,i,j,wt_s2(k),P_surf_flux(i,j),T_scale,u_star
          end do
  
        else if (lbcp==0) then
!          ustar_local=sum(sqrt(u(1:nx,1:ny,1)**2+v(1:nx,1:ny,1)**2))/float(nx*ny)*vonK/&
!       	       (dlog(0.5_rprec*dz/zo)-sum(psi_m(1:nx,1:ny))/float(nx*ny))
          ustar_local=ustar_avg
!          print *,'using temp bc for patches !!'
!          P_surf_flux(i,j)=(S_Surf(i,j)-scalar(i,j,1))*vonk*ustar2(i,j)&
!       		  /(phi_h(i,j)*dlog(dz/(2._rprec*z_os(i,j))))

          
  	  ! For pollen
	  IF (settling) THEN
!	    term1=(dz/(2._rprec*z_os(i,j)))**(-settling_vel/(vonk*ustar_local(i,j)))
	    term1=((2._rprec*z_os(i,j))/dz)**(Sc_boun*settling_vel/(vonk*ustar_local(i,j)))
	    P_surf_flux(i,j)=((S_Surf(i,j)*term1-scalar(i,j,1))*settling_vel)/(1._rprec-term1)
  	  ELSE
	    P_surf_flux(i,j)=((S_Surf(i,j)-scalar(i,j,1))*vonk*ustar_local(i,j))/(Sc_boun*dlog(dz/(2._rprec*z_os(i,j))))
	  END IF
	  

!  equ. 4.28. in brutsaert
! Impose the gradient value at the first node.
! Note that ustar2 is a variable coming from sgs stresses at the wall,
! ustar2 has a single value for all the nodes in a single patch.
! phi_h is obtained from the routine obukhov.f, surf_flux is as computed above
! and vonk and dz are constants. Everything is in dimensionless form
            
    	end if

        ! Now we have the lowest dsdz on the UVP nodes all others on w nodes
!        dsdz(i,j,1) =-phi_h(i,j)*P_surf_flux(i,j)/(ustar2(i,j)*vonk*DZ/2._rprec)
!        dsdz(i,j,1) =-phi_h(i,j)*P_surf_flux(i,j)/(ustar_local(i,j)*vonk*DZ/2._rprec)

        



        ! Calculate vertical derivative at first vertical node
        ! Note: in the finite volume formualtion this is never used to update the concentration field
        IF (settling) THEN
          dsdz(i,j,1)=(-Sc_boun*P_surf_flux(i,j)-settling_vel*scalar(i,j,1))/(ustar_local(i,j)*vonk*DZ/2._rprec)
        ELSE
          dsdz(i,j,1)=-(Sc_boun*P_surf_flux(i,j))/(ustar_local(i,j)*vonk*dz/2._rprec)
        END IF
	

	! Do not allow positive fluxes for point source case
	! i.e. fluxes from the ground when PCon becomes negative due to Gibbs
	IF (pointsource .AND. P_surf_flux(i,j)>0) THEN
	  P_surf_flux(i,j)=0._rprec
	  dsdz(i,j,1)=0._rprec
	END IF
	
	
	! Calculate pollen deposition and source
        deposition(i,j)=deposition(i,j)+(-dt*(1.5_rprec*P_surf_flux(i,j)-0.5_rprec*surf_flux_m(i,j)))
	

      end do
    end do
    
    
!    PRINT*,MAXVAL(scalar),MINVAL(scalar)


  end if !end for if (USE_MPI ...) block
  
  
  
  
  call dealias1(u,u_m)
  call dealias1(v,v_m)
  call dealias1(w,w_m)
  call dealias1(dsdx,dsdx_m)
  call dealias1(dsdy,dsdy_m)
  call dealias1(dsdz,dsdz_m)



  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    jz_min = 2
  else
    jz_min = 1
  end if
 

  ubc_jz = nz-1


! Now compute the RHS term of the filtered scalar equation. 
! Note that this is the advection term with the scalar as 
! the diffusion term has been thrown away. This is done step 
! by step for each of the expanded arrays from dealias1 separately
! for the first node & last node AND the rest of the nodes.

! xxxxx ------Comments valid for MPI case only ---------XXXX
! For MPI case, all the nodes have fluid nodes (1:nz-1) except for
! near-the wall processes (2:nz-1 for w nodes and 1:nz-1 for uvp nodes)
! and the top nodes (1:nz)
! The following loop executes from 2:nz-1 for the near-wall process
! and 1:nz-1 for the other processes. Also, the subsequent 2 loops
! take care of the first node for the near-wall process (coord = 0)
! and the topmost node for the top process (coord = nproc-1).
! Note also that we need to make ghost data available for dTdz and
! w for the topmost node (jz=n) within each process and therefore,
! this synchronization (MPI_SENDRECV()) has been done in the subroutine
! theta_all_in_one ()
! xxxxx --------- MPI Comment block ends ------------------XXXX

  do k=jz_min,nz-1
    do j=1,Ny2
      do i=1,Nx2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Modified by Di Yang for ocean LES with Stokes drift
!DY Start here
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY  	RHS_m(i,j,k)=u_m(i,j,k)*dsdx_m(i,j,k)+v_m(i,j,k)*dsdy_m(i,j,k)&
!DY  	    +(w_m(i,j,k)*dsdz_m(i,j,k)+w_m(i,j,k+1)*dsdz_m(i,j,k+1))/2._rprec
         if(OCEAN_FLAG .and. STOKES_FLAG) then
            RHS_m(i,j,k)=(u_m(i,j,k)+ust(k))*dsdx_m(i,j,k)+v_m(i,j,k)*dsdy_m(i,j,k)&
                 +(w_m(i,j,k)*dsdz_m(i,j,k)+w_m(i,j,k+1)*dsdz_m(i,j,k+1))/2._rprec
         else
            RHS_m(i,j,k)=u_m(i,j,k)*dsdx_m(i,j,k)+v_m(i,j,k)*dsdy_m(i,j,k)&
                 +(w_m(i,j,k)*dsdz_m(i,j,k)+w_m(i,j,k+1)*dsdz_m(i,j,k+1))/2._rprec
         endif
!+++++++++++++++
!DY End here
!+++++++++++++++
      end do
    end do
  end do


  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    do j=1,Ny2
      do i=1,Nx2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Modified by Di Yang for ocean LES with Stokes drift
!DY Start here
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY  	RHS_m(i,j,1)=u_m(i,j,1)*dsdx_m(i,j,1)+v_m(i,j,1)*dsdy_m(i,j,1)&
!DY  		     +(0.5_rprec*w_m(i,j,2))*dsdz_m(i,j,2)
         if(OCEAN_FLAG .and. STOKES_FLAG) then
            RHS_m(i,j,1)=(u_m(i,j,1)+ust(1))*dsdx_m(i,j,1)+v_m(i,j,1)*dsdy_m(i,j,1)&
                 +(0.5_rprec*w_m(i,j,2))*dsdz_m(i,j,2)
         else
            RHS_m(i,j,1)=u_m(i,j,1)*dsdx_m(i,j,1)+v_m(i,j,1)*dsdy_m(i,j,1)&
                 +(0.5_rprec*w_m(i,j,2))*dsdz_m(i,j,2)
         endif
!+++++++++++++++
!DY End here
!+++++++++++++++
      end do
    end do
  end if
  
  
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
    do j=1,Ny2
      do i=1,Nx2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Modified by Di Yang for ocean LES with Stokes drift
!DY Start here
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY  	RHS_m(i,j,Nz)=u_m(i,j,Nz)*dsdx_m(i,j,Nz)+v_m(i,j,Nz)*dsdy_m(i,j,Nz)
         if(OCEAN_FLAG .and. STOKES_FLAG) then
            RHS_m(i,j,Nz)=(u_m(i,j,Nz)+ust(Nz))*dsdx_m(i,j,Nz)+v_m(i,j,Nz)*dsdy_m(i,j,Nz)
         else
            RHS_m(i,j,Nz)=u_m(i,j,Nz)*dsdx_m(i,j,Nz)+v_m(i,j,Nz)*dsdy_m(i,j,Nz)
         endif
!+++++++++++++++
!DY End here
!+++++++++++++++
      end do
    end do
  end if

!!$  do k=jz_min,nz-1
!!$    do j=1,Ny2
!!$      do i=1,Nx2
!!$  	RHS_m(i,j,k)=u_m(i,j,k)*dsdx_m(i,j,k)+v_m(i,j,k)*dsdy_m(i,j,k)&
!!$  	    +(w_m(i,j,k)*dsdz_m(i,j,k)+w_m(i,j,k+1)*dsdz_m(i,j,k+1))/2._rprec
!!$      end do
!!$    end do
!!$  end do
!!$
!!$
!!$  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
!!$    do j=1,Ny2
!!$      do i=1,Nx2
!!$  	RHS_m(i,j,1)=u_m(i,j,1)*dsdx_m(i,j,1)+v_m(i,j,1)*dsdy_m(i,j,1)&
!!$  		     +(0.5_rprec*w_m(i,j,2))*dsdz_m(i,j,2)
!!$      end do
!!$    end do
!!$  end if
!!$  
!!$  
!!$  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
!!$    do j=1,Ny2
!!$      do i=1,Nx2
!!$  	RHS_m(i,j,Nz)=u_m(i,j,Nz)*dsdx_m(i,j,Nz)+v_m(i,j,Nz)*dsdy_m(i,j,Nz)
!!$      end do
!!$    end do
!!$  end if
    



  call dealias2(RHS,RHS_m)


!c...Now building the SGS part of the RHS.
! Here the sgs_term for scalars is built up using Nu_t from sgs_stag_W.f
! and dividing it by the turbulent Prandtl # specified in dimen.h
!c....Note: Since we bring the Convective term to RHS its sign changes.
!c....Below "Temp" is used for SGS flux; its divergence is added to RHS
!VK.. Nu_t is on w nodes everywhere except at z=dz/2.. while
!VK dsdx is on uvp nodes.. so, interpolate Nu_t as we want temp to
!VK be on uvp nodes
! All this valid only till Pr_ is a constant..
! This will need a makeover once Pr_ becomes dynamic as well...

  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    ! at jz=1, both Nu_t and dsdx are on uvp nodes .. no need for interp
    do j=1,Ny
      do i=1,Nx
        temp(i,j,1)=(1._rprec/Sc)*Nu_t(i,j,1)*dsdx(i,j,1)
      end do
    end do
  end if

  do k=jz_min,ubc_jz
    do j=1,Ny
      do i=1,Nx
        temp(i,j,k)=(1._rprec/Sc)*0.5_rprec*(Nu_t(i,j,k)+Nu_t(i,j,k+1))*dsdx(i,j,k)
      end do
    end do
  end do

  call DDX (dtemp, temp)  

  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    ! at jz=1, both Nu_t and dsdy are on uvp nodes .. no need for interp
    do j=1,Ny
      do i=1,Nx
        RHS(i,j,1) = (-1._rprec*RHS(i,j,1) + dtemp(i,j,1))
        temp(i,j,1)=(1._rprec/Sc)*Nu_t(i,j,1)*dsdy(i,j,1)
      end do
    end do
  end if

  do k=jz_min,ubc_jz
    ! Nu_t is on w nodes and dsdy is on uvp nodes for jz=2 to nz
    do j=1,Ny
      do i=1,Nx
        RHS(i,j,k) = (-1._rprec*RHS(i,j,k) + dtemp(i,j,k))
        temp(i,j,k)=(1._rprec/Sc)*0.5_rprec*(Nu_t(i,j,k)+Nu_t(i,j,k+1))*dsdy(i,j,k)
      end do
    end do
  end do

  call DDY (dtemp, temp)   
  
  !...Use MO flux at wall for the scalar sgs term !
  ! Note that the total contribution to the scalar sgs term at
  ! the first node comes from the surface flux computed above from
  ! the specified heat flux, wt_s

  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    do j=1,Ny
      do i=1,Nx
        RHS(i,j,1)=RHS(i,j,1)+dtemp(i,j,1)
        temp(i,j,1)=-1._rprec*P_surf_flux(i,j)
        sgs_vert(i,j,1)=P_surf_flux(i,j)
      end do
    end do
  end if


  ! Note sgs_vert is -1*temp because sgs_vert is modeled as -Nu_t*dsdz/Pr
  ! while temp is the same term but w/o the minus sign due to the additional
  ! minus outside the scalar fluctuation flux term in RHS
  ! need to run this loop nz due to the nature of the differenetiation in ddz_w

  ! Modified to account for pollen settling velocity term
  ! Chamecki 10/12/2006
  IF (settling) THEN

    DO k=jz_min,nz
      DO j=1,Ny
        DO i=1,Nx
          
	  ! Add y-component of turbulent diffusion to RHS
	  RHS(i,j,k)=RHS(i,j,k)+dtemp(i,j,k)
	  
	  ! Calculate z-component of turbulent diffusion
          temp(i,j,k)=(1._rprec/Sc)*Nu_t(i,j,k)*dsdz(i,j,k)
	  
	  ! Add the gravitational settling part
          temp(i,j,k)=temp(i,j,k)+settling_vel*(scalar(i,j,k-1)+scalar(i,j,k))/2._rprec
	  
          ! Store vertical flux
	  sgs_vert(i,j,k)=-1._rprec*temp(i,j,k)

        END DO
      END DO
    END DO
    
  ELSE

    DO k=jz_min,nz
      DO j=1,Ny
        DO i=1,Nx
          RHS(i,j,k) = RHS(i,j,k) + dtemp(i,j,k)
          temp(i,j,k)=(1._rprec/Sc)*Nu_t(i,j,k)*dsdz(i,j,k)
          sgs_vert(i,j,k)=-1._rprec*temp(i,j,k)
        END DO
      END DO
    END DO
  
  END IF
  
  ! Impose zero vertical flux at top boundary
  temp(:,:,nz)=0._rprec
  sgs_vert(:,:,nz)=0._rprec


  ! The SGS_z flux is on the W nodes, but DDZ_W will put it back on UVP nodes! 
  ! Also note that sgs_vert(i,j,k) influences the computations in 
  ! OBUKHOV.f and is not involved in any computations in this routine.
  ! sgs_t3(i,j,1) (<w'theta'> is used for computing wt at the surface in OBUKHOV)

  call DDZ_w (dtemp, temp)

  do k=1,ubc_jz
    do j=1,Ny
      do i=1,Nx
        RHS(i,j,k)=RHS(i,j,k)+dtemp(i,j,k)
      end do
    end do
  end do
  

!  PRINT*,MAXVAL(RHS(1:nx,1:ny,1)),MAXLOC(RHS(1:nx,1:ny,1))
!  PRINT*,MAXVAL(scalar(1:nx,1:ny,1)),MAXLOC(scalar(1:nx,1:ny,1))
  
!  PRINT*,RHS(5,15,1),scalar(5,15,1)

!  PRINT*,scalar(4,15,1),scalar(5,15,1),scalar(6,15,1)
!  PRINT*,scalar(5,14,1),scalar(5,15,1),scalar(5,16,1)
!  PRINT*,scalar(5,15,1),scalar(5,15,2)
  
!  IF (jt_total==2) STOP

  
  ! Add source term if needed
  IF (pointsource .AND. (jt_total>=ini_src .AND. jt_total<end_src)) THEN
      
    if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then

    IF (single_point) THEN    
    ! Source in 1 point
!    RHS(xps,yps,zps)=RHS(xps,yps,zps)+Q_src



!    ! Divide source in 9 points (horizontal plane)
!    ! Center (=1/4)
!    RHS(xps,yps,zps)=RHS(xps,yps,zps)+Q_src/4
!    ! Sides (=1/8)
!    RHS(xps-1,yps,zps)=RHS(xps-1,yps,zps)+Q_src/8
!    RHS(xps+1,yps,zps)=RHS(xps+1,yps,zps)+Q_src/8
!    RHS(xps,yps-1,zps)=RHS(xps,yps-1,zps)+Q_src/8
!    RHS(xps,yps+1,zps)=RHS(xps,yps+1,zps)+Q_src/8
!    ! Corners (=1/16)
!    RHS(xps-1,yps-1,zps)=RHS(xps-1,yps-1,zps)+Q_src/16
!    RHS(xps-1,yps+1,zps)=RHS(xps-1,yps+1,zps)+Q_src/16
!    RHS(xps+1,yps-1,zps)=RHS(xps+1,yps-1,zps)+Q_src/16
!    RHS(xps+1,yps+1,zps)=RHS(xps+1,yps+1,zps)+Q_src/16




    ! Divide source in 9 points (horizontal plane)
    ! Center (=1/4)
    RHS(xps,yps,zps)=RHS(xps,yps,zps)+Q_src/4/2
    ! Sides (=1/8)
    RHS(xps-1,yps,zps)=RHS(xps-1,yps,zps)+Q_src/8/2
    RHS(xps+1,yps,zps)=RHS(xps+1,yps,zps)+Q_src/8/2
    RHS(xps,yps-1,zps)=RHS(xps,yps-1,zps)+Q_src/8/2
    RHS(xps,yps+1,zps)=RHS(xps,yps+1,zps)+Q_src/8/2
    ! Corners (=1/16)
    RHS(xps-1,yps-1,zps)=RHS(xps-1,yps-1,zps)+Q_src/16/2
    RHS(xps-1,yps+1,zps)=RHS(xps-1,yps+1,zps)+Q_src/16/2
    RHS(xps+1,yps-1,zps)=RHS(xps+1,yps-1,zps)+Q_src/16/2
    RHS(xps+1,yps+1,zps)=RHS(xps+1,yps+1,zps)+Q_src/16/2

    ! Divide source in 9 points (horizontal plane)
    ! Center (=1/4)
    RHS(xps,yps,zps-1)=RHS(xps,yps,zps-1)+Q_src/4/4
    ! Sides (=1/8)
    RHS(xps-1,yps,zps-1)=RHS(xps-1,yps,zps-1)+Q_src/8/4
    RHS(xps+1,yps,zps-1)=RHS(xps+1,yps,zps-1)+Q_src/8/4
    RHS(xps,yps-1,zps-1)=RHS(xps,yps-1,zps-1)+Q_src/8/4
    RHS(xps,yps+1,zps-1)=RHS(xps,yps+1,zps-1)+Q_src/8/4
    ! Corners (=1/16)
    RHS(xps-1,yps-1,zps-1)=RHS(xps-1,yps-1,zps-1)+Q_src/16/4
    RHS(xps-1,yps+1,zps-1)=RHS(xps-1,yps+1,zps-1)+Q_src/16/4
    RHS(xps+1,yps-1,zps-1)=RHS(xps+1,yps-1,zps-1)+Q_src/16/4
    RHS(xps+1,yps+1,zps-1)=RHS(xps+1,yps+1,zps-1)+Q_src/16/4

    ! Divide source in 9 points (horizontal plane)
    ! Center (=1/4)
    RHS(xps,yps,zps+1)=RHS(xps,yps,zps+1)+Q_src/4/4
    ! Sides (=1/8)
    RHS(xps-1,yps,zps+1)=RHS(xps-1,yps,zps+1)+Q_src/8/4
    RHS(xps+1,yps,zps+1)=RHS(xps+1,yps,zps+1)+Q_src/8/4
    RHS(xps,yps-1,zps+1)=RHS(xps,yps-1,zps+1)+Q_src/8/4
    RHS(xps,yps+1,zps+1)=RHS(xps,yps+1,zps+1)+Q_src/8/4 
    ! Corners (=1/16)
    RHS(xps-1,yps-1,zps+1)=RHS(xps-1,yps-1,zps+1)+Q_src/16/4
    RHS(xps-1,yps+1,zps+1)=RHS(xps-1,yps+1,zps+1)+Q_src/16/4
    RHS(xps+1,yps-1,zps+1)=RHS(xps+1,yps-1,zps+1)+Q_src/16/4
    RHS(xps+1,yps+1,zps+1)=RHS(xps+1,yps+1,zps+1)+Q_src/16/4

    ELSE

    DO j=1,ny
      DO i=1,nx
        RHS(i,j,zps)=RHS(i,j,zps)+Q_src
      END DO
    END DO 

    END IF

    end if

  END IF
  

end subroutine pollen_RHS_calc





! This is a simplified version of inflow_cond() in forcing.f90
! Chamecki 08/21/2006
subroutine inflow_pollen()
use types, only : rprec
use param, only : jx_p, ny, nz, pi, x_relaxp
use sim_param, only : PCon
implicit none

integer :: jx
integer :: i
integer :: jx_0

real (rprec) :: factor

!---------------------------------------------------------------------



  ! Initial point for damping
  jx_0 = jx_p - x_relaxp
  
  ! Changed from i=1 to i=0 to include last point in the domain
  do i = 0, x_relaxp-1
    jx = jx_p - i
    factor = 0.5_rprec * (1._rprec - cos (pi * real (i, rprec) / x_relaxp))
    PCon(jx, 1:ny, 1:nz) =  factor * PCon(jx_0, 1:ny, 1:nz)
  end do


end subroutine inflow_pollen




! This is a simplified version of step_scalar()
! Chamecki 10/04/2006
subroutine step_pollen(scalar,RHS_pre,RHS_post)
implicit none
integer:: i,j,k
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

  real(kind=rprec),dimension(ld,ny,$lbz:nz)::scalar, RHS_pre, RHS_post

  do k=1,nz-1
    do j=1,ny
      do i=1,nx
        scalar(i,j,k)= scalar(i,j,k)+dt*(1.5_rprec*RHS_pre(i,j,k)-0.5_rprec*RHS_post(i,j,k))
      end do
    end do
  end do     

  scalar(:,:,Nz)=scalar(:,:,Nz-1)

end subroutine step_pollen























! Here is the finite volumes discretization
! For now, all the routines are specific
! Afer everything is working, I should clean up all the pollen
! routines and put them together in a nice and organized final version.

! Main routine for pollen time advance
! Chamecki - 03/29/2007

SUBROUTINE pollen_all_in_one2

  use test_filtermodule

  IMPLICIT NONE

  INTEGER:: i, j
  INTEGER :: jz,ios,counter

  REAL(kind=rprec):: sim_time,data_time,data_time_prev
  REAL(kind=rprec):: data_ustar,data_ustar_prev
  REAL(kind=rprec),DIMENSION(3)     :: data_Co,data_Co_prev
  REAL(kind=rprec),DIMENSION(nx,ny) :: PCon_sfc_temp
  REAL(kind=rprec):: airborne_global,deposited_global

$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif



  ! Main loop to advance PCon in time
  IF (PCon_FLAG) THEN
  
    ! Store previous time step RHS for AB2 time integration
    RHS_PConf=RHS_PCon
    
    ! Calculate vertical derivative just as diagnostic variable
    ! Note: in the finite volume formualtion this is never used to update the concentration field
    CALL ddz_uv(dPCondz,PCon) 
    
$if ($MPI)
  !print *,'synchronizing in pollen _all_in_one2 for coord = ',coord  ! Need to synchronize w and dPCondz across the processors for uvp node
  ! based computation in scalar_RHS_calc (calculation of RHS_m)
  call mpi_sendrecv (w(1, 1, 1), ld*ny, MPI_RPREC, down, tag_counter+1,  &
                     w(1, 1, nz), ld*ny, MPI_RPREC, up, tag_counter+1,   &
                     comm, status, ierr)  
  call mpi_sendrecv (dudt(1, 1, 1), ld*ny, MPI_RPREC, down, tag_counter+2,  &
                     dudt(1, 1, nz), ld*ny, MPI_RPREC, up, tag_counter+2,   &
                     comm, status, ierr)
  call mpi_sendrecv (dvdt(1, 1, 1), ld*ny, MPI_RPREC, down, tag_counter+3,  &
                     dvdt(1, 1, nz), ld*ny, MPI_RPREC, up, tag_counter+3,   &
                     comm, status, ierr)
  call mpi_sendrecv (dwdt(1, 1, 1), ld*ny, MPI_RPREC, down, tag_counter+4,  &
                     dwdt(1, 1, nz), ld*ny, MPI_RPREC, up, tag_counter+4,   &
                     comm, status, ierr)
  ! Also need to synchronize Nu_t across the processors for   ! computation in scalar_RHS_calc (calculation of RHS_m)
  call mpi_sendrecv (Nu_t(1, 1, 1), ld*ny, MPI_RPREC, down, tag_counter+5,  &
                     Nu_t(1, 1, nz), ld*ny, MPI_RPREC, up, tag_counter+5,   &
                     comm, status, ierr)  
  call mpi_sendrecv (Nu_t(1, 1, nz-1), ld*ny, MPI_RPREC, up, tag_counter+6,  &
                     Nu_t(1, 1, 0), ld*ny, MPI_RPREC, down, tag_counter+6,   &
                     comm, status, ierr)
  call mpi_sendrecv (magS(1, 1, 1), ld*ny, MPI_RPREC, down, tag_counter+7,  &
                     magS(1, 1, nz), ld*ny, MPI_RPREC, up, tag_counter+7,   &
                     comm, status, ierr)
  call mpi_sendrecv (magS(1, 1, nz-1), ld*ny, MPI_RPREC, up, tag_counter+8,  &
                     magS(1, 1, 0), ld*ny, MPI_RPREC, down, tag_counter+8,   &
                     comm, status, ierr)
$endif
    
    ! This is to use time evolution of u* and Co for the simulation of the ragweed field experiment
    ! Note that scale_us is used in main.f90 to scale the pressure forcing
    IF (rag06) THEN
    
      ! Determine current time (in seconds)
      sim_time=jt_total*dt_dim
      
      ! Initialize time from file
      data_time=-1._rprec
         
      ! Open file and read appropriate values for ustar and Co    
      OPEN(1,FILE=path//'LES_bc.dat',STATUS='OLD',ACTION='READ')
      READ(1,*)
      READ(1,*)
      DO WHILE (data_time<=sim_time)
      
        ! Save previous values
        data_time_prev=data_time
        data_ustar_prev=data_ustar
        data_Co_prev(:)=data_Co(:)
        
        ! Read next value
        READ(1,*)data_time,data_ustar,data_Co(1),data_Co(2),data_Co(3)
      
      END DO
      
      ! Calculate sclaing factors using linear interpolation
      scale_us=data_ustar_prev+((sim_time-data_time_prev)/(data_time-data_time_prev))*(data_ustar-data_ustar_prev)
      scale_Co3(:)=data_Co_prev(:)+((sim_time-data_time_prev)/(data_time-data_time_prev))*(data_Co(:)-data_Co_prev(:))
     
      ! Get dimensionless values
      scale_us=scale_us/u_star
      scale_Co3(:)=scale_Co3(:)/PCon_scale
      
      ! Scale the surface concentration
      PCon_sfc_temp=0._rprec
      DO j=1,ny

        ! Linear interpolation in the y direction
        ! Very specific to the bc and grid for the ragfield simulation
        IF (j<=30) THEN
          scale_Co=scale_Co3(1)
        ELSE IF (j<=35) THEN
          scale_Co=(1._rprec-((j-30)*2-1)/10._rprec)*scale_Co3(1)+(((j-30)*2-1)/10._rprec)*scale_Co3(2)
        ELSE IF (j<=40) THEN
          scale_Co=(1._rprec-((j-35)*2-1)/10._rprec)*scale_Co3(2)+(((j-35)*2-1)/10._rprec)*scale_Co3(3)
        ELSE
          scale_Co=scale_Co3(3)
        END IF

        WHERE (PCon_s(:,j)/=0._rprec)
          PCon_sfc_temp(:,j)=scale_Co*PCon_s(:,j)
        END WHERE

!        DO i=1,nx
!          IF (PCon_s(i,j)/=0._rprec) PCon_sfc_temp(i,j)=scale_Co*PCon_s(i,j)
!        END DO

      END DO


!      OPEN(1,FILE=path//'output/check_Psrc.out',STATUS="unknown",POSITION="append")
!      DO i=1,nx
!        WRITE(1,'(51F18.8)')i*dx,(PCon_sfc_temp(i,j),j=1,ny)
!      END DO
!      CLOSE(1)


      OPEN(1,FILE=path//'output/scaling.out',STATUS="unknown",POSITION="append")
      WRITE(1,'(I8,2F16.6)')jt_total,(scale_us*u_star),(scale_Co*PCon_scale)
      CLOSE(1)

      
    ELSE
      
      IF (jt_total>=ini_src .AND. jt_total<end_src) THEN
        PCon_sfc_temp=PCon_s
      ELSE
        PCon_sfc_temp=0._rprec
      END IF
      
    END IF

      RHS_PConf=RHS_PCon

    ! Calculate RHS of equation
    ! The total flux out of the domain is also calculated inside the routine
!    CALL pollen_RHS_calc2(PCon,dPCondz,PCon_s,zo_PCon,RHS_PCon,sgs_PCon3,res_PCon3,sfc_flux)
    CALL pollen_RHS_calc2(PCon,dPCondz,PCon_sfc_temp,zo_PCon,RHS_PCon,sgs_PCon3,res_PCon3,sfc_flux,sgs_PCon1)
    
            
    ! For first time step, RHS(n-1)=RHS(n) for time integration
    IF (((jt .eq. PCon_init) .or. (pointsource .and. (jt .eq. ini_src))) &
       .and. (.not. initPCon)) THEN
      RHS_PConf=RHS_PCon
!++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for testing
!DY Start here
      print*, coord,"pointsource start"
!DY End here
!++++++++++++++++++++++++++++++++++++++
    END IF

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for buoyancy term added to vertical momentum equation
!DY Start here
!!$    if(active_pcon) then
!!$       call calcbeta_pcon(PCon)
!!$    endif
!DY End here
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! Time integration using AB2
    CALL step_pollen2(PCon,RHS_PCon,RHS_PConf)

$if ($MPI)
  call mpi_sendrecv (PCon(1, 1, 1), ld*ny, MPI_RPREC, down, tag_counter+7,  &
                     PCon(1, 1, nz), ld*ny, MPI_RPREC, up, tag_counter+7,   &
                     comm, status, ierr)
  call mpi_sendrecv (PCon(1, 1, nz-1), ld*ny, MPI_RPREC, up, tag_counter+8,  &
                     PCon(1, 1, 0), ld*ny, MPI_RPREC, down, tag_counter+8,   &
                     comm, status, ierr)
$endif

    ! Check mass balance

    ! Grains released (point source) - uses AB2
    IF (pointsource .and. (jt_total >= ini_src)) THEN
      released=1.5_rprec*((jt_total-ini_src+1)*dt*Q_src*dx*dy*dz) &
              -0.5_rprec*((jt_total-ini_src)*dt*Q_src*dx*dy*dz)
    END IF
    IF (pointsource .and. (jt_total >= end_src)) THEN
      released=((end_src-ini_src)*dt*Q_src*dx*dy*dz)
    END IF
    

    ! Grains airborne and deposited in the valid domain
    airborne=(SUM(PCon(1:nx,1:ny,1:nz-1))*dx*dy*dz)
    deposited=(SUM(deposition(1:nx,1:ny))*dx*dy)
    gone=gone+1.5_rprec*flux_out-0.5_rprec*flux_out_prev
    
    IF (periodicbcx .AND. periodicbcy) gone=0._rprec

$if ($MPI)

  call mpi_reduce (airborne,airborne_global,1,MPI_RPREC,MPI_SUM,0,comm,ierr)
  call mpi_reduce (deposited,deposited_global,1,MPI_RPREC,MPI_SUM,0,comm,ierr)
  if (rank == 0) then
    OPEN(1,FILE=path//'output/pollen_balance.out',STATUS="unknown",POSITION="append")
    WRITE(1,'(I8,6G18.8)')jt_total,(released), &                           ! Pollen grains released (#)
                             (airborne_global), &                          ! Pollen grains airborne (#)
                             (deposited_global), &                         ! Pollen grains deposited (#)
                             (gone), &                                     ! Pollen grains out of domain (#)
                             ((released-airborne-deposited-gone))


    CLOSE(1)
  endif

$else

! Use this if there is no inflow condition
!    airborne=(SUM(PCon(1:nx,:,1:nz-1))*dx*dy*dz)
!    deposited=(SUM(deposition(1:nx,:))*dx*dy)
!    gone=0._rprec
        
    
!    OPEN(1,FILE=path//'output/pollen_balance.out',STATUS="unknown",POSITION="append")
!    WRITE(1,'(I8,6G18.8)')jt_total,(released*PCon_scale*z_i*z_i*z_i), &                           ! Pollen grains released (#)
!                             (airborne*PCon_scale*z_i*z_i*z_i), &                                 ! Pollen grains airborne (#)
!			     (deposited*PCon_scale*z_i*z_i*z_i), &                                ! Pollen grains deposited (#)
!			     (gone*PCon_scale*z_i*z_i*z_i), &                                     ! Pollen grains out of domain (#)
!			     ((released-airborne-deposited-gone)*PCon_scale*z_i*z_i*z_i)
!
!
!    CLOSE(1)

    
    ! Dimensionless balance
    OPEN(1,FILE=path//'output/pollen_balance.out',STATUS="unknown",POSITION="append")
    WRITE(1,'(I8,6G18.8)')jt_total,(released), &                           ! Pollen grains released (#)
                             (airborne), &                                 ! Pollen grains airborne (#)
			     (deposited), &                                ! Pollen grains deposited (#)
			     (gone), &                                     ! Pollen grains out of domain (#)
			     ((released-airborne-deposited-gone))


    CLOSE(1)

$endif

  END IF
  

END SUBROUTINE pollen_all_in_one2




!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Add by Di Yang
!DY Calculate buoyancy term due to particle concentration
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine calcbeta_pcon (PCon)
! This calculates the buoyancy term due to oil concentration (beta_pcon) 
! to be added to the vertical momentum equation
! Authored by Di Yang

implicit none
integer::i, j, k, jz_min
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
real(kind=rprec),dimension(ld,ny,$lbz:nz),intent(in)::PCon
real(kind=rprec)::g_hat

  !..Non-dimensionalize gravity
  g_hat=g*(z_i/(u_star**2))
  beta_pcon=0._rprec

  do k=$lbz,Nz-1
    do j=1,Ny
      do i=1,nx
        beta_pcon(i,j,k)=g_hat*(densratio_pcon-1._rprec)*V_pcon*Pcon(i,j,k)
      end do
    end do     
  end do
  return
  
end subroutine calcbeta_pcon
!+++++++++++++++++++++++++++++++++++++++++
!DY Subroutine calcbeta_pcon end here
!+++++++++++++++++++++++++++++++++++++++++





!
!  This is the finite volumes version of RHS calculation
!  This routine include both QUICK and SMART schemes with periodic and inflow/outflow bc's.
!  The main reference for the discretization of the advection term is the paper by 
!    Waterson and Deconinck (1995) cited by Xie et al JoT (2004).
!

subroutine pollen_RHS_calc2(scalar,dsdz,S_Surf,z_os,RHS,sgs_vert,res_vert,surf_flux_current,sgs_x)


  use sim_param,only: dudx,dvdy,dwdz

  implicit none

  $if ($MPI)
    $define $lbz 0
  $else
    $define $lbz 1
  $endif

  integer:: i, j, k, jz
  integer:: jz_min,ubc_jz
  real(kind=rprec):: surf_flux_current,term1
  real(kind=rprec):: term_alpha,term_xi,term_tao,term_y,term_r,term_omega

  real(kind=rprec),dimension(ld,ny,$lbz:nz)  :: RHS
  real(kind=rprec),dimension(ld,ny,$lbz:nz)  :: scalar
  real(kind=rprec),dimension(ld,ny,$lbz:nz)  :: dsdz
  real(kind=rprec),dimension(ld,ny,$lbz:nz)  :: sgs_vert                     ! Store SGS vertical pollen flux
  real(kind=rprec),dimension(ld,ny,$lbz:nz)  :: res_vert                     ! Store resolved vertical pollen flux
  real(kind=rprec),dimension(ld,ny,$lbz:nz)  :: sgs_x                     ! Store SGS x-direction pollen flux
  
  real(kind=rprec)                           :: rx,ry,rz                     ! Gradient ratios in x, y and z
  real(kind=rprec),dimension(ld,ny,nz)       :: temp_int                     ! To help with interpolations
  REAL(kind=rprec),DIMENSION(nx+1,ny+1,nz)   :: scalar_x,scalar_y,scalar_z   ! Interpolated scalar in x, y and z directions
  REAL(kind=rprec),DIMENSION(ld,ny+1,nz)     :: u_int,v_int,w_int 	   ! Interpolated velocity field for convective term
  real(kind=rprec),dimension(ny,nz-1)        :: ghost_x0,ghost_xLx           ! Ghost nodes at x=-dx/2 and x=Lx+dx/2
  real(kind=rprec),dimension(nx,nz-1)        :: ghost_y0,ghost_yLy           ! Ghost nodes at y=-dy/2 and y=Ly+dy/2
  real(kind=rprec),dimension(ld,ny)          :: ghost_z0
  real(kind=rprec),dimension(ld,ny,$lbz:nz)       :: magS_int                    ! Interpolated |S|
  

  real(kind=rprec),dimension(nx,ny)          :: ustar_local,S_Surf,z_os

  real(kind=rprec),dimension(nx,ny)          :: surf_flux_m                  ! Previous surface flux for correct calculation of deposition
                                                                             ! (i.e. AB2 deposition integration)
  real(kind=rprec),dimension(nx,ny)          :: surf_flux_src_m              ! Previous surface flux for total flux

  real(kind=rprec),dimension(ld,ny):: scalar_node_1 ! test filtered and used for computing surface flux

  real(kind=rprec),dimension (ptypes)        :: ustar_patch                  ! u* averaged for each patch
  real(kind=rprec),dimension (ptypes)        :: sfc_flx_patch                ! Surface flux for each patch

  real(kind=rprec)                           :: delta                        ! Filter size

  integer:: ncpu_source
  real(kind=rprec):: zps_local 

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for ocean LES with Stokes drift
real(kind=rprec),dimension($lbz:nz) :: ust
real (rprec) :: z
!DY End here
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for ocean LES with Stokes drift
!DY Start here
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if(OCEAN_FLAG .and. STOKES_FLAG) then
   do jz=$lbz,nz
      $if ($MPI)
      z = (coord*(nz-1) + jz - 0.5_rprec) * dz
      $else
      z=(real(jz)-0.5_rprec)*dz
      $endif
      ust(jz)=U_stokes*exp(-2._rprec*wavenm_w*z)
   enddo
endif
!+++++++++++++++
!DY End here
!+++++++++++++++

  !
  ! Step 1 - Prepare bottom boundary condition
  !
  
  ! For now the imposed surface flux is the same for all patches
  DO k=1,ptypes
    sfc_flx_patch(k)=surf_flux_current 
  END DO

  ! Initialize patch-averaged u*
  ustar_patch=0._rprec

  ! Calculate patch-averaged u*
  DO j=1,ny
    DO i=1,nx
      DO k=1,ptypes
  	IF (patch(i,j)==k) THEN
  	  ustar_patch(patch(i,j))=ustar_patch(patch(i,j))+ustar_avg(i,j)/patchnum(patch(i,j))
  	END IF
      END DO
    END DO
  END DO


  ! Assign u* to each node based on patch type
  DO j=1,ny
    DO i=1,nx
      DO k=1,ptypes
  	IF (patch(i,j)==k) THEN
  	  ustar_local(i,j)=ustar_patch(k)
  	END IF
      END DO
    END DO
  END DO
  
  ! Store previous surface flux to integrate deposition
  surf_flux_m=P_surf_flux
  surf_flux_src_m=P_surf_flux_dep

  ! Compute surface flux and at z=dz/2 using similarity theory
  DO j=1,ny
    DO i=1,nx

      ! lbc=1 is used for prescribing the surface flux while
      IF (lbcp==1) THEN
  	
        ! Assign surface flux based on patch type
        DO k=1,ptypes
  	  IF (patch(i,j)==k) P_surf_flux(i,j)=sfc_flx_patch(k)
  	END DO
        
        ! If imposed flux, src and total are the same
        P_surf_flux_dep=P_surf_flux
  
      ! lbc=0 has been used to prescribe the temperature
      ELSE IF (lbcp==0) THEN

  	! Use average u* instead of local value based on patch type
        ustar_local=ustar_avg

        ! Calculate surface flux for neutral atmospheric stability
        IF (settling) THEN
!          term1=((2._rprec*z_os(i,j))/dz)**(Sc_boun*settling_vel/(vonk*ustar_local(i,j)))
          term1=((z_os(i,j)-d0(i,j))/(dz/2._rprec-d0(i,j)))**(Sc_boun*settling_vel/(vonk*ustar_local(i,j)))
          P_surf_flux(i,j)=((S_Surf(i,j)*term1-scalar(i,j,1))*settling_vel)/(1._rprec-term1)
          ! Calculate effects of atmospheric stability
          IF (wt_s .ne. 0._rprec) THEN
            term_alpha = Sc_boun*settling_vel/vonk/ustar_local(i,j)
            term_xi = (dz/2._rprec-d0(i,j))/L(i,j)
            term_omega = 1._rprec
            IF (L(i,j) .lt. 0.0_rprec) THEN ! For unstable atmosphere
              term_tao = 16._rprec*term_xi*(0.5_rprec-term_alpha)-(1._rprec+term_alpha)
              term_y = 2._rprec*term_alpha  &
                     / ((term_tao**2._rprec  &
                     - 4._rprec*term_alpha*16._rprec*term_xi  &
                     * (1._rprec+term_alpha-0.5_rprec))**0.5_rprec  &
                     - term_tao)
              term_r = term_y**2._rprec/term_alpha  &
                     + (1._rprec-term_y)**2._rprec  &
                     - 0.5_rprec*(16._rprec*term_xi)**2._rprec  &
                     / (1._rprec-16._rprec*term_xi*term_y)**2._rprec  &
                     * term_y**2._rprec / term_alpha  &
                     * (1._rprec-term_y)**2._rprec
              term_omega = (1._rprec+term_alpha)**(1._rprec+term_alpha-0.5_rprec)  &
                         * term_r**(-0.5_rprec)  &
                         * (term_y/term_alpha)**term_alpha  &
                         * (1._rprec-term_y)  &
                         * (1._rprec-16._rprec*term_xi*term_y)**(-0.5_rprec)
            END IF
            IF (L(i,j) .gt. 0._rprec) THEN ! For stable atmosphere
              term_omega = 1._rprec+5._rprec*(term_alpha/(term_alpha+1._rprec))*term_xi
            END IF
            P_surf_flux(i,j)=P_surf_flux(i,j)/term_omega
          END IF
        ELSE
          term1=dlog((dz/2._rprec-d0(i,j))/(z_os(i,j)-d0(i,j)))
          P_surf_flux(i,j)=((S_Surf(i,j)-scalar(i,j,1))*vonk*ustar_local(i,j))/(Sc_boun*term1)
        END IF

        ! Calculate real deposition surface flux for neutral atmospheric stability
        ! This is the surface flux if there is Cr=0 everywhere in the domain and is used to characterize the deposition
        IF (settling) THEN
!          term1=((2._rprec*z_os(i,j))/dz)**(Sc_boun*settling_vel/(vonk*ustar_local(i,j)))
          term1=((z_os(i,j)-d0(i,j))/(dz/2._rprec-d0(i,j)))**(Sc_boun*settling_vel/(vonk*ustar_local(i,j)))
          P_surf_flux_dep(i,j)=((0._rprec-scalar(i,j,1))*settling_vel)/(1._rprec-term1)
        ELSE
          term1=dlog((dz/2._rprec-d0(i,j))/(z_os(i,j)-d0(i,j)))
          P_surf_flux_dep(i,j)=((0._rprec-scalar(i,j,1))*vonk*ustar_local(i,j))/(Sc_boun*term1)
        END IF
        
      END IF


      ! Calculate vertical derivative at first vertical node
      ! Note: in the finite volume formualtion this is never used to update the concentration field
      IF ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
      IF (settling) THEN
        dsdz(i,j,1)=(-Sc_boun*P_surf_flux(i,j)-settling_vel*scalar(i,j,1))/(ustar_local(i,j)*vonk*(dz/2._rprec-d0(i,j)))
      ELSE
        dsdz(i,j,1)=-(Sc_boun*P_surf_flux(i,j))/(ustar_local(i,j)*vonk*(dz/2._rprec-d0(i,j)))
      END IF
      END IF


      ! Do not allow positive fluxes for point source case
      ! i.e. fluxes from the ground when PCon becomes negative due to Gibbs
      IF (pointsource .AND. P_surf_flux(i,j)>0) THEN
        P_surf_flux(i,j)=0._rprec
      END IF
      
      ! For first time step, RHS(n-1)=RHS(n) for time integration
      IF ((jt .eq. PCon_init) .or. (pointsource .and. (jt .eq. ini_src))) THEN
        surf_flux_m=P_surf_flux
        surf_flux_src_m=P_surf_flux_dep
      END IF      
      
      ! Calculate pollen deposition (time integration using AB2 scheme)
      deposition(i,j)=deposition(i,j)+(-dt*(1.5_rprec*P_surf_flux(i,j)-0.5_rprec*surf_flux_m(i,j)))

      ! Calculate pollen source
      Real_dep(i,j)=Real_dep(i,j)+(-dt*(1.5_rprec*P_surf_flux_dep(i,j)-0.5_rprec*surf_flux_src_m(i,j)))

    END DO
  END DO
  

  !
  ! Step 2 - Prepare boundary conditions (face values or ghost nodes)
  !


  IF (periodicbcx) THEN

    ! For periodic boundary condition in x and y
    ! The ghost nodes are just explicit enforcement of periodicity
    ghost_x0(1:ny,1:nz-1)=scalar(nx,1:ny,1:nz-1)
    ghost_xLx(1:ny,1:nz-1)=scalar(1,1:ny,1:nz-1)
    
    ! The face values have to be calculated in the interpolation scheme below
    
  ELSE

    ! For inflow condition/outflow
    ! The ghost nodes enforce zero derivatives
    ghost_xLx(1:ny,1:nz-1)=scalar(nx,1:ny,1:nz-1)
    
    ! The face value enforces zero concentration at the inlet
    scalar_x(1,1:ny,1:nz-1)=0._rprec

    ! The ghost node is also set to be zero
    ! The alternative approach would be to calculate its value based 
    ! on linear extrapolation for the advective term and 2nd order
    ! differences for the diffusion term
    ghost_x0(1:ny,1:nz-1)=0._rprec
    
    ! The other face values are set based on a "rough" approximation
    ! I have to rethink this approach
    scalar_x(nx+1,1:ny,1:nz-1)=scalar(nx,1:ny,1:nz-1)
  
  END IF


  IF (periodicbcy) THEN

    ! For periodic boundary condition in x and y
    ! The ghost nodes are just explicit enforcement of periodicity
    ghost_y0(1:nx,1:nz-1)=scalar(1:nx,ny,1:nz-1)
    ghost_yLy(1:nx,1:nz-1)=scalar(1:nx,1,1:nz-1)
    
    ! The face values have to be calculated in the interpolation scheme below
    
  ELSE

    ! For inflow condition/outflow
    ! The ghost nodes enforce zero derivatives
    ghost_y0(1:nx,1:nz-1)=scalar(1:nx,1,1:nz-1)
!    ghost_yLy(1:nx,1:nz-1)=scalar(1:nx,ny,1:nz-1)
    
    ! The other face values are set based on a "rough" approximation
    ! I have to rethink this approach
    scalar_y(1:nx,1,1:nz-1)=scalar(1:nx,1,1:nz-1)
!    scalar_y(1:nx,ny+1,1:nz-1)=scalar(1:nx,ny,1:nz-1)



    ! For ragwee field simulation modify b.c.
    scalar_y(1:nx,ny+1,1:nz-1)=0._rprec
    ghost_yLy(1:nx,1:nz-1)=0._rprec
  
  END IF

  ghost_z0(1:ld,1:ny)=0._rprec  

  !
  ! Step 3 - Interpolations
  !
  

  !
  ! 3a - u-component
  !  
  
  ! Assemble inverted matrix for x-direction
  matrix_x=0._rprec
  DO i=1,(nx-1)
    DO j=1,i
      matrix_x(i,j)=-(nx-i)*(1._rprec/nx)
    END DO
    DO j=(i+1),nx
      matrix_x(i,j)=i*(1._rprec/nx)
    END DO
  END DO
  matrix_x(nx,:)=1._rprec/nx
  
  
  ! Calculate u-component of velocity
  ! Loop over all yz-locations
  DO j=1,ny
    DO k=1,nz
  
      ! Assemble derivative vector
      dvector_x(1:nx-1)=dx*dudx(1:nx-1,j,k)
      dvector_x(nx)=SUM(u(1:nx,j,k))
  
      ! Calculate velocity field     
      DO i=1,nx
        u_int(i,j,k)=SUM(matrix_x(:,i)*dvector_x(:))
      END DO
      
    END DO
  END DO
  
  ! Use periodicity
  u_int(nx+1,1:ny,1:nz)=u_int(1,1:ny,1:nz)


  !
  ! 3b - v-component
  !  
  
  ! Assemble inverted matrix for y-direction
  matrix_y=0._rprec
  DO i=1,(ny-1)
    DO j=1,i
      matrix_y(i,j)=-(ny-i)*(1._rprec/ny)
    END DO
    DO j=(i+1),ny
      matrix_y(i,j)=i*(1._rprec/ny)
    END DO
  END DO
  matrix_y(ny,:)=1._rprec/ny
  
  
  ! Calculate v-component of velocity
  ! Loop over all xz-locations
  DO i=1,nx
    DO k=1,nz
  
      ! Assemble derivative vector
      dvector_y(1:ny-1)=dy*dvdy(i,1:ny-1,k)
      dvector_y(ny)=SUM(v(i,1:ny,k))
  
      ! Calculate velocity field     
      DO j=1,ny
        v_int(i,j,k)=SUM(matrix_y(:,j)*dvector_y(:))
      END DO
      
    END DO
  END DO
  
  ! Use periodicity
  v_int(1:nx,ny+1,1:nz)=v_int(1:nx,1,1:nz)


  !
  ! 3c - w-component
  !  
  
  ! w is already on the right position
  w_int(1:nx,1:ny,1:nz)=w(1:nx,1:ny,1:nz)
  
  
  
  ! This is the new interpolation for |S| - center of finite volumes
  ! (interpolation only in z)

  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
    ! To avoid problems in the upper boundary
    magS(1:nx,1:ny,nz)=magS(1:nx,1:ny,nz-1)
  end if
  
  magS_int(1:nx,1:ny,2:nz-1)=(magS(1:nx,1:ny,2:nz-1)+magS(1:nx,1:ny,3:nz))/2._rprec

  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    ! The first level does not need any interpolation
    magS_int(1:nx,1:ny,1)=magS(1:nx,1:ny,1)
  else
    magS_int(1:nx,1:ny,1)=(magS(1:nx,1:ny,1)+magS(1:nx,1:ny,$lbz))/2._rprec
  end if

$if($MPI)
  call mpi_sendrecv (magS_int(1, 1, 1), ld*ny, MPI_RPREC, down, tag_counter+5,  &
                     magS_int(1, 1, nz), ld*ny, MPI_RPREC, up, tag_counter+5,   &
                     comm, status, ierr)
  call mpi_sendrecv (magS_int(1, 1, nz-1), ld*ny, MPI_RPREC, up, tag_counter+6,  &
                     magS_int(1, 1, 0), ld*ny, MPI_RPREC, down, tag_counter+6,   &
                     comm, status, ierr)
$endif

    
  ! Add settling velocity to interpolated vertical velocity
  ! Note: w_int @ k=1 and k=nz is zero (imposed B.C.)
  !       there is no need to add the settling for these two locations since
  !       on the boundaries the settling is accounted for together with the 
  !       vertical turbulent diffusion (they add to zero at the top and to the
  !       theoretical equilibrium profile at the surface).
  IF (settling) THEN
    w_int(1:nx,1:ny,2:nz-1)=w_int(1:nx,1:ny,2:nz-1)-settling_vel
    if (USE_MPI .and. coord > 0) then
      w_int(1:nx,1:ny,1)=w_int(1:nx,1:ny,1)-settling_vel
    end if
    if (USE_MPI .and. coord < nproc-1) then
      w_int(1:nx,1:ny,nz)=w_int(1:nx,1:ny,nz)-settling_vel
    end if
    IF (PCon_acc) THEN
      ! Multiply by taup=ws/g (use dimensionless gravity)
      ! Add contribution to u_int (linear interpolation to faces)
      u_int(2:nx,1:ny,1:nz)=u_int(2:nx,1:ny,1:nz)-0.5D0*(dudt(1:nx-1,1:ny,1:nz)+dudt(2:nx,1:ny,1:nz))*settling_vel/(9.81D0/(u_star**2/z_i))
      u_int(1,1:ny,1:nz)=u_int(1,1:ny,1:nz)-0.5D0*(dudt(nx,1:ny,1:nz)+dudt(1,1:ny,1:nz))*settling_vel/(9.81D0/(u_star**2/z_i))
      u_int(nx+1,1:ny,1:nz)=u_int(nx+1,1:ny,1:nz)-0.5D0*(dudt(nx,1:ny,1:nz)+dudt(1,1:ny,1:nz))*settling_vel/(9.81D0/(u_star**2/z_i))

      ! Add contribution to v_int (linear interpolation to faces)
      v_int(1:nx,2:ny,1:nz)=v_int(1:nx,2:ny,1:nz)-0.5D0*(dvdt(1:nx,1:ny-1,1:nz)+dvdt(1:nx,2:ny,1:nz))*settling_vel/(9.81D0/(u_star**2/z_i))
      v_int(1:nx,1,1:nz)=v_int(1:nx,1,1:nz)-0.5D0*(dvdt(1:nx,ny,1:nz)+dvdt(1:nx,1,1:nz))*settling_vel/(9.81D0/(u_star**2/z_i))
      v_int(1:nx,ny+1,1:nz)=v_int(1:nx,ny+1,1:nz)-0.5D0*(dvdt(1:nx,ny,1:nz)+dvdt(1:nx,1,1:nz))*settling_vel/(9.81D0/(u_star**2/z_i))

      ! Add contribution to w_int (no interpolation needed - see note above about settling velocity)
      w_int(1:nx,1:ny,2:nz-1)=w_int(1:nx,1:ny,2:nz-1)-dwdt(1:nx,1:ny,2:nz-1)*settling_vel/(9.81D0/(u_star**2/z_i))
      if (USE_MPI .and. coord > 0) then
        w_int(1:nx,1:ny,1)=w_int(1:nx,1:ny,1)-dwdt(1:nx,1:ny,1)*settling_vel/(9.81D0/(u_star**2/z_i))
      end if
      if (USE_MPI .and. coord < nproc-1) then
        w_int(1:nx,1:ny,nz)=w_int(1:nx,1:ny,nz)-dwdt(1:nx,1:ny,nz)*settling_vel/(9.81D0/(u_star**2/z_i))
      end if
    END IF
  END IF
  

  ! Calculate the interpolated concentrations
  ! This is actually where QUICK, SMART, etc are different
  
  !
  ! QUICK scheme
  !
  
  IF (PCon_scheme==2) THEN
  
    ! Interpolated concentrations in x
    DO k=1,nz-1
      DO j=1,ny
              
        ! interior nodes
        DO i=3,nx-1
          ! Choose the appropriate upwind direction
          IF (u_int(i,j,k) >= 0._rprec) THEN  
            scalar_x(i,j,k)=(3._rprec*scalar(i,j,k)+6._rprec*scalar(i-1,j,k)-scalar(i-2,j,k))/8._rprec
          ELSE      
            scalar_x(i,j,k)=(3._rprec*scalar(i-1,j,k)+6._rprec*scalar(i,j,k)-scalar(i+1,j,k))/8._rprec
          END IF
        END DO
        
        ! nodes using BC
        
        ! Second node (i=2)
        ! Choose the appropriate upwind direction
        IF (u_int(2,j,k) >= 0._rprec) THEN  
          scalar_x(2,j,k)=(3._rprec*scalar(2,j,k)+6._rprec*scalar(1,j,k)-ghost_x0(j,k))/8._rprec
        ELSE    
          scalar_x(2,j,k)=(3._rprec*scalar(1,j,k)+6._rprec*scalar(2,j,k)-scalar(3,j,k))/8._rprec
        END IF

        ! Before last node (i=nx)
        ! Choose the appropriate upwind direction
        IF (u_int(nx,j,k) >= 0._rprec) THEN  
          scalar_x(nx,j,k)=(3._rprec*scalar(nx,j,k)+6._rprec*scalar(nx-1,j,k)-scalar(nx-2,j,k))/8._rprec
        ELSE    
          scalar_x(nx,j,k)=(3._rprec*scalar(nx-1,j,k)+6._rprec*scalar(nx,j,k)-ghost_xLx(j,k))/8._rprec
        END IF

      END DO        
    END DO
    
    
    
    ! Interpolated concentrations in y
    DO k=1,nz-1
      DO i=1,nx
              
        ! interior nodes
        DO j=3,ny-1       
          ! Choose the appropriate upwind direction
          IF (v_int(i,j,k) >= 0._rprec) THEN  
            scalar_y(i,j,k)=(3._rprec*scalar(i,j,k)+6._rprec*scalar(i,j-1,k)-scalar(i,j-2,k))/8._rprec
          ELSE      
            scalar_y(i,j,k)=(3._rprec*scalar(i,j-1,k)+6._rprec*scalar(i,j,k)-scalar(i,j+1,k))/8._rprec
          END IF
        END DO
        
        ! nodes using BC
        
        ! Second node (j=2)

        ! Choose the appropriate upwind direction
        IF (v_int(i,2,k) >= 0._rprec) THEN  
          scalar_y(i,2,k)=(3._rprec*scalar(i,2,k)+6._rprec*scalar(i,1,k)-ghost_y0(i,k))/8._rprec
        ELSE    
          scalar_y(i,2,k)=(3._rprec*scalar(i,1,k)+6._rprec*scalar(i,2,k)-scalar(i,3,k))/8._rprec
        END IF      

        ! Before last node (j=ny)
        ! Choose the appropriate upwind direction
        IF (v_int(i,ny,k) >= 0._rprec) THEN  
          scalar_y(i,ny,k)=(3._rprec*scalar(i,ny,k)+6._rprec*scalar(i,ny-1,k)-scalar(i,ny-2,k))/8._rprec
        ELSE    
          scalar_y(i,ny,k)=(3._rprec*scalar(i,ny-1,k)+6._rprec*scalar(i,ny,k)-ghost_yLy(i,k))/8._rprec
        END IF

      END DO        
    END DO


    ! Interpolated concentrations in z
    ! Note: interpolated values at k=1 and k=nz are not needed, since the fluxes are
    ! imposed through BC and w'=0 anyway
    DO k=3,nz-2
      DO j=1,ny
        DO i=1,nx        
        ! Choose the appropriate upwind direction
        IF (w_int(i,j,k) >= 0._rprec) THEN  
          scalar_z(i,j,k)=(3._rprec*scalar(i,j,k)+6._rprec*scalar(i,j,k-1)-scalar(i,j,k-2))/8._rprec
        ELSE      
          scalar_z(i,j,k)=(3._rprec*scalar(i,j,k-1)+6._rprec*scalar(i,j,k)-scalar(i,j,k+1))/8._rprec
        END IF
        END DO
      END DO        
    END DO
    
    ! For the nodes k=2 and k=nz-1
    DO j=1,ny
      DO i=1,nx
        
        ! k=2
        ! Choose the appropriate upwind direction
        IF (w_int(i,j,2) >= 0._rprec) THEN  
!          scalar_z(i,j,2)=(3._rprec*scalar(i,j,2)+6._rprec*scalar(i,j,1)-scalar(i,j,0))/8._rprec
        ! For now use a simple interpolation in this case
          if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
          scalar_z(i,j,2)=(scalar(i,j,2)+scalar(i,j,1))/2._rprec
          else
          scalar_z(i,j,2)=(3._rprec*scalar(i,j,2)+6._rprec*scalar(i,j,1)-scalar(i,j,$lbz))/8._rprec
          endif
        ELSE    
          scalar_z(i,j,2)=(3._rprec*scalar(i,j,1)+6._rprec*scalar(i,j,2)-scalar(i,j,3))/8._rprec
        END IF

        ! k=nz-1
        ! Choose the appropriate upwind direction
        IF (w_int(i,j,nz-1) >= 0._rprec) THEN  
          scalar_z(i,j,nz-1)=(3._rprec*scalar(i,j,nz-1)+6._rprec*scalar(i,j,nz-2)-scalar(i,j,nz-3))/8._rprec
        ELSE
          if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then    
!          scalar_z(i,j,nz-1)=(3._rprec*scalar(i,j,nz-2)+6._rprec*scalar(i,j,nz-1)-scalar(i,j,nz))/8._rprec
          scalar_z(i,j,nz-1)=(3._rprec*scalar(i,j,nz-2)+6._rprec*scalar(i,j,nz-1)-0._rprec)/8._rprec
          else
          scalar_z(i,j,nz-1)=(3._rprec*scalar(i,j,nz-2)+6._rprec*scalar(i,j,nz-1)-scalar(i,j,nz))/8._rprec
          end if
        END IF

      END DO
    END DO      



    ! Compute face nodes for periodic boundary conditions
    IF (periodicbcx) THEN
    
      DO k=1,nz-1
        DO j=1,ny
                        
          ! First node (i=1)
	  ! Choose the appropriate upwind direction
          IF (u_int(1,j,k) >= 0._rprec) THEN  
            scalar_x(1,j,k)=(3._rprec*scalar(1,j,k)+6._rprec*scalar(nx,j,k)-scalar(nx-1,j,k))/8._rprec
          ELSE      
            scalar_x(1,j,k)=(3._rprec*scalar(nx,j,k)+6._rprec*scalar(1,j,k)-scalar(2,j,k))/8._rprec
          END IF
          
          ! Last node (i=nx+1)
 	  ! Choose the appropriate upwind direction
          IF (u_int(nx+1,j,k) >= 0._rprec) THEN  
            scalar_x(nx+1,j,k)=(3._rprec*scalar(1,j,k)+6._rprec*scalar(nx,j,k)-scalar(nx-1,j,k))/8._rprec
          ELSE	  
            scalar_x(nx+1,j,k)=(3._rprec*scalar(nx,j,k)+6._rprec*scalar(1,j,k)-scalar(2,j,k))/8._rprec
          END IF

        END DO        
      END DO

    END IF


    IF (periodicbcy) THEN

      DO k=1,nz-1
        DO i=1,nx

          ! First node (j=1)
          ! Choose the appropriate upwind direction
          IF (v_int(i,1,k) >= 0._rprec) THEN  
            scalar_y(i,1,k)=(3._rprec*scalar(i,1,k)+6._rprec*scalar(i,ny,k)-scalar(i,ny-1,k))/8._rprec
          ELSE	
            scalar_y(i,1,k)=(3._rprec*scalar(i,ny,k)+6._rprec*scalar(i,1,k)-scalar(i,2,k))/8._rprec
          END IF      

          ! Last node (j=ny+1)
          ! Choose the appropriate upwind direction
          IF (v_int(i,ny+1,k) >= 0._rprec) THEN  
            scalar_y(i,ny+1,k)=(3._rprec*scalar(i,1,k)+6._rprec*scalar(i,ny,k)-scalar(i,ny-1,k))/8._rprec
          ELSE	
            scalar_y(i,ny+1,k)=(3._rprec*scalar(i,ny,k)+6._rprec*scalar(i,1,k)-scalar(i,2,k))/8._rprec
          END IF

        END DO        
      END DO

    END IF

        

  !
  ! SMART scheme
  !  
    
  ELSE IF (PCon_scheme==3) THEN
  
    ! Interpolated concentrations in x
    DO k=1,nz-1
      DO j=1,ny
              
        ! interior nodes
        DO i=3,nx-1
          ! Choose the appropriate upwind direction
          IF (u_int(i,j,k) >= 0._rprec) THEN
            ! To avoid problems of divisions by zero
            IF ((scalar(i-1,j,k)-scalar(i-2,j,k))==0._rprec) THEN
              scalar_x(i,j,k)=scalar(i-1,j,k)
            ELSE
              ! Gradient ratio
              rx=(scalar(i,j,k)-scalar(i-1,j,k))/(scalar(i-1,j,k)-scalar(i-2,j,k))
  	      ! Bound the gradient
	      rx=MAX(0._rprec,MIN(MIN(2._rprec*rx,0.75_rprec*rx+0.25_rprec),4._rprec))
	      ! Intrepolate the value
              scalar_x(i,j,k)=scalar(i-1,j,k)+rx*(scalar(i-1,j,k)-scalar(i-2,j,k))/2._rprec
            END IF
          ELSE
            ! To avoid problems of divisions by zero
            IF ((scalar(i,j,k)-scalar(i+1,j,k))==0._rprec) THEN
              scalar_x(i,j,k)=scalar(i,j,k)
            ELSE
              ! Gradient ratio
              rx=(scalar(i-1,j,k)-scalar(i,j,k))/(scalar(i,j,k)-scalar(i+1,j,k))
	      ! Bound the gradient
	      rx=MAX(0._rprec,MIN(MIN(2._rprec*rx,0.75_rprec*rx+0.25_rprec),4._rprec))
              ! Intrepolate the value
              scalar_x(i,j,k)=scalar(i,j,k)+rx*(scalar(i,j,k)-scalar(i+1,j,k))/2._rprec
            END IF
          END IF
        END DO
        
        ! nodes using BC
        
        ! Second node (i=2)
        ! Choose the appropriate upwind direction
        IF (u_int(2,j,k) >= 0._rprec) THEN  
          ! To avoid problems of divisions by zero
          IF ((scalar(1,j,k)-ghost_x0(j,k))==0._rprec) THEN
            scalar_x(2,j,k)=scalar(1,j,k)
          ELSE
            ! Gradient ratio
            rx=(scalar(2,j,k)-scalar(1,j,k))/(scalar(1,j,k)-ghost_x0(j,k))
            ! Bound the gradient
            rx=MAX(0._rprec,MIN(MIN(2._rprec*rx,0.75_rprec*rx+0.25_rprec),4._rprec))
            ! Intrepolate the value
            scalar_x(2,j,k)=scalar(1,j,k)+rx*(scalar(1,j,k)-ghost_x0(j,k))/2._rprec
          END IF
        ELSE    
          ! To avoid problems of divisions by zero
          IF ((scalar(2,j,k)-scalar(3,j,k))==0._rprec) THEN
            scalar_x(2,j,k)=scalar(2,j,k)
          ELSE
            ! Gradient ratio
            rx=(scalar(1,j,k)-scalar(2,j,k))/(scalar(2,j,k)-scalar(3,j,k))
            ! Bound the gradient
            rx=MAX(0._rprec,MIN(MIN(2._rprec*rx,0.75_rprec*rx+0.25_rprec),4._rprec))
            ! Intrepolate the value
            scalar_x(2,j,k)=scalar(2,j,k)+rx*(scalar(2,j,k)-scalar(3,j,k))/2._rprec
          END IF
        END IF

        ! Before last node (i=nx)
        ! Choose the appropriate upwind direction
        IF (u_int(nx,j,k) >= 0._rprec) THEN  
          ! To avoid problems of divisions by zero
          IF ((scalar(nx-1,j,k)-scalar(nx-2,j,k))==0._rprec) THEN
            scalar_x(nx,j,k)=scalar(nx-1,j,k)
          ELSE
            ! Gradient ratio
            rx=(scalar(nx,j,k)-scalar(nx-1,j,k))/(scalar(nx-1,j,k)-scalar(nx-2,j,k))
            ! Bound the gradient
            rx=MAX(0._rprec,MIN(MIN(2._rprec*rx,0.75_rprec*rx+0.25_rprec),4._rprec))
            ! Intrepolate the value
            scalar_x(nx,j,k)=scalar(nx-1,j,k)+rx*(scalar(nx-1,j,k)-scalar(nx-2,j,k))/2._rprec
          END IF
        ELSE    
          ! To avoid problems of divisions by zero
          IF ((scalar(nx,j,k)-ghost_xLx(j,k))==0._rprec) THEN
            scalar_x(nx,j,k)=scalar(nx,j,k)
          ELSE
            ! Gradient ratio
            rx=(scalar(nx-1,j,k)-scalar(nx,j,k))/(scalar(nx,j,k)-ghost_xLx(j,k))
            ! Bound the gradient
            rx=MAX(0._rprec,MIN(MIN(2._rprec*rx,0.75_rprec*rx+0.25_rprec),4._rprec))
            ! Intrepolate the value
            scalar_x(nx,j,k)=scalar(nx,j,k)+rx*(scalar(nx,j,k)-ghost_xLx(j,k))/2._rprec
          END IF
        END IF

      END DO        
    END DO
    
    
    
    ! Interpolated concentrations in y
    DO k=1,nz-1
      DO i=1,nx
              
        ! interior nodes
        DO j=3,ny-1       
          ! Choose the appropriate upwind direction
          IF (v_int(i,j,k) >= 0._rprec) THEN  
            ! To avoid problems of divisions by zero
            IF ((scalar(i,j-1,k)-scalar(i,j-2,k))==0._rprec) THEN
              scalar_y(i,j,k)=scalar(i,j-1,k)
            ELSE
              ! Gradient ratio
              ry=(scalar(i,j,k)-scalar(i,j-1,k))/(scalar(i,j-1,k)-scalar(i,j-2,k))
              ! Bound the gradient
              ry=MAX(0._rprec,MIN(MIN(2._rprec*ry,0.75_rprec*ry+0.25_rprec),4._rprec))
              ! Intrepolate the value
              scalar_y(i,j,k)=scalar(i,j-1,k)+ry*(scalar(i,j-1,k)-scalar(i,j-2,k))/2._rprec
            END IF
          ELSE      
            ! To avoid problems of divisions by zero
            IF ((scalar(i,j,k)-scalar(i,j+1,k))==0._rprec) THEN
              scalar_y(i,j,k)=scalar(i,j,k)
            ELSE
              ! Gradient ratio
              ry=(scalar(i,j-1,k)-scalar(i,j,k))/(scalar(i,j,k)-scalar(i,j+1,k))
              ! Bound the gradient
              ry=MAX(0._rprec,MIN(MIN(2._rprec*ry,0.75_rprec*ry+0.25_rprec),4._rprec))
              ! Intrepolate the value
              scalar_y(i,j,k)=scalar(i,j,k)+ry*(scalar(i,j,k)-scalar(i,j+1,k))/2._rprec
            END IF
          END IF
        END DO
        
        ! nodes using BC
        
        ! Second node (j=2)
        ! Choose the appropriate upwind direction
        IF (v_int(i,2,k) >= 0._rprec) THEN  
          ! To avoid problems of divisions by zero
          IF ((scalar(i,1,k)-ghost_y0(i,k))==0._rprec) THEN
            scalar_y(i,2,k)=scalar(i,1,k)
          ELSE
            ! Gradient ratio
            ry=(scalar(i,2,k)-scalar(i,1,k))/(scalar(i,1,k)-ghost_y0(i,k))
            ! Bound the gradient
            ry=MAX(0._rprec,MIN(MIN(2._rprec*ry,0.75_rprec*ry+0.25_rprec),4._rprec))
            ! Intrepolate the value
            scalar_y(i,2,k)=scalar(i,1,k)+ry*(scalar(i,1,k)-ghost_y0(i,k))/2._rprec
          END IF
        ELSE      
          ! To avoid problems of divisions by zero
          IF ((scalar(i,2,k)-scalar(i,3,k))==0._rprec) THEN
            scalar_y(i,2,k)=scalar(i,2,k)
          ELSE
            ! Gradient ratio
            ry=(scalar(i,1,k)-scalar(i,2,k))/(scalar(i,2,k)-scalar(i,3,k))
            ! Bound the gradient
            ry=MAX(0._rprec,MIN(MIN(2._rprec*ry,0.75_rprec*ry+0.25_rprec),4._rprec))
            ! Intrepolate the value
            scalar_y(i,2,k)=scalar(i,2,k)+ry*(scalar(i,2,k)-scalar(i,3,k))/2._rprec
          END IF
        END IF


        ! Before last node (j=ny)
        ! Choose the appropriate upwind direction
        IF (v_int(i,ny,k) >= 0._rprec) THEN  
          ! To avoid problems of divisions by zero
          IF ((scalar(i,ny-1,k)-scalar(i,ny-2,k))==0._rprec) THEN
            scalar_y(i,ny,k)=scalar(i,ny-1,k)
          ELSE
            ! Gradient ratio
            ry=(scalar(i,ny,k)-scalar(i,ny-1,k))/(scalar(i,ny-1,k)-scalar(i,ny-2,k))
            ! Bound the gradient
            ry=MAX(0._rprec,MIN(MIN(2._rprec*ry,0.75_rprec*ry+0.25_rprec),4._rprec))
            ! Intrepolate the value
            scalar_y(i,ny,k)=scalar(i,ny-1,k)+ry*(scalar(i,ny-1,k)-scalar(i,ny-2,k))/2._rprec
          END IF
        ELSE      
          ! To avoid problems of divisions by zero
          IF ((scalar(i,ny,k)-ghost_yLy(i,k))==0._rprec) THEN
            scalar_y(i,ny,k)=scalar(i,ny,k)
          ELSE
            ! Gradient ratio
            ry=(scalar(i,ny-1,k)-scalar(i,ny,k))/(scalar(i,ny,k)-ghost_yLy(i,k))
            ! Bound the gradient
            ry=MAX(0._rprec,MIN(MIN(2._rprec*ry,0.75_rprec*ry+0.25_rprec),4._rprec))
            ! Intrepolate the value
            scalar_y(i,ny,k)=scalar(i,ny,k)+ry*(scalar(i,ny,k)-ghost_yLy(i,k))/2._rprec
          END IF
        END IF

      END DO        
    END DO


    ! Interpolated concentrations in z
    ! Note: interpolated values at k=1 and k=nz are not needed, since the fluxes are
    ! imposed through BC and w'=0 anyway
    DO k=3,nz-2
      DO j=1,ny
        DO i=1,nx        
          ! Choose the appropriate upwind direction
          IF (w_int(i,j,k) >= 0._rprec) THEN  
            ! To avoid problems of divisions by zero
            IF ((scalar(i,j,k-1)-scalar(i,j,k-2))==0._rprec) THEN
              scalar_z(i,j,k)=scalar(i,j,k-1)
            ELSE
              ! Gradient ratio
              rz=(scalar(i,j,k)-scalar(i,j,k-1))/(scalar(i,j,k-1)-scalar(i,j,k-2))
              ! Bound the gradient
              rz=MAX(0._rprec,MIN(MIN(2._rprec*rz,0.75_rprec*rz+0.25_rprec),4._rprec))
              ! Intrepolate the value
              scalar_z(i,j,k)=scalar(i,j,k-1)+rz*(scalar(i,j,k-1)-scalar(i,j,k-2))/2._rprec
            END IF
          ELSE      
            ! To avoid problems of divisions by zero
            IF ((scalar(i,j,k)-scalar(i,j,k+1))==0._rprec) THEN
              scalar_z(i,j,k)=scalar(i,j,k)
            ELSE
              ! Gradient ratio
              rz=(scalar(i,j,k-1)-scalar(i,j,k))/(scalar(i,j,k)-scalar(i,j,k+1))
              ! Bound the gradient
              rz=MAX(0._rprec,MIN(MIN(2._rprec*rz,0.75_rprec*rz+0.25_rprec),4._rprec))
              ! Intrepolate the value
              scalar_z(i,j,k)=scalar(i,j,k)+rz*(scalar(i,j,k)-scalar(i,j,k+1))/2._rprec
            END IF
          END IF
        END DO
      END DO        
    END DO
    
    ! For the nodes k=2 and k=nz-1
    DO j=1,ny
      DO i=1,nx
        
        ! k=2
        ! Choose the appropriate upwind direction
        IF (w_int(i,j,2) >= 0._rprec) THEN  

          if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
!          ! For now use a simple interpolation in this case
            scalar_z(i,j,2)=(scalar(i,j,2)+scalar(i,j,1))/2._rprec
          else
           ! To avoid problems of divisions by zero
           IF ((scalar(i,j,1)-scalar(i,j,$lbz)==0._rprec)) THEN
             scalar_z(i,j,2)=scalar(i,j,1)
           ELSE
             ! Gradient ratio
             rz=(scalar(i,j,2)-scalar(i,j,1))/(scalar(i,j,1)-scalar(i,j,$lbz))
             ! Bound the gradient
             rz=MAX(0._rprec,MIN(MIN(2._rprec*rz,0.75_rprec*rz+0.25_rprec),4._rprec))
             ! Interpolate the value
             scalar_z(i,j,2)=scalar(i,j,1)+rz*(scalar(i,j,1)-scalar(i,j,$lbz))/2._rprec
           ENDIF
          endif
        ELSE      
          ! To avoid problems of divisions by zero
          IF ((scalar(i,j,2)-scalar(i,j,3))==0._rprec) THEN
            scalar_z(i,j,2)=scalar(i,j,2)
          ELSE
            ! Gradient ratio
            rz=(scalar(i,j,1)-scalar(i,j,2))/(scalar(i,j,2)-scalar(i,j,3))
            ! Bound the gradient
            rz=MAX(0._rprec,MIN(MIN(2._rprec*rz,0.75_rprec*rz+0.25_rprec),4._rprec))
            ! Intrepolate the value
            scalar_z(i,j,2)=scalar(i,j,2)+rz*(scalar(i,j,2)-scalar(i,j,3))/2._rprec
          END IF
        END IF

        ! k=nz-1
        ! Choose the appropriate upwind direction
        IF (w_int(i,j,nz-1) >= 0._rprec) THEN  
          ! To avoid problems of divisions by zero
          IF ((scalar(i,j,nz-2)-scalar(i,j,nz-3))==0._rprec) THEN
            scalar_z(i,j,nz-1)=scalar(i,j,nz-2)
          ELSE
            ! Gradient ratio
            rz=(scalar(i,j,nz-1)-scalar(i,j,nz-2))/(scalar(i,j,nz-2)-scalar(i,j,nz-3))
            ! Bound the gradient
            rz=MAX(0._rprec,MIN(MIN(2._rprec*rz,0.75_rprec*rz+0.25_rprec),4._rprec))
            ! Intrepolate the value
            scalar_z(i,j,nz-1)=scalar(i,j,nz-2)+rz*(scalar(i,j,nz-2)-scalar(i,j,nz-3))/2._rprec
          END IF
        ELSE      
          ! To avoid problems of divisions by zero
          IF ((scalar(i,j,nz-1)-scalar(i,j,nz))==0._rprec) THEN
            scalar_z(i,j,nz-1)=scalar(i,j,nz-1)
          ELSE
            ! Gradient ratio
            rz=(scalar(i,j,nz-2)-scalar(i,j,nz-1))/(scalar(i,j,nz-1)-scalar(i,j,nz))
            ! Bound the gradient
            rz=MAX(0._rprec,MIN(MIN(2._rprec*rz,0.75_rprec*rz+0.25_rprec),4._rprec))
            ! Intrepolate the value
            scalar_z(i,j,nz-1)=scalar(i,j,nz-1)+rz*(scalar(i,j,nz-1)-scalar(i,j,nz))/2._rprec
          END IF
        END IF

      END DO
    END DO      

$if ($MPI)
           call mpi_sendrecv(scalar(1,1,nz-2),ld*ny,MPI_RPREC,up,tag_counter+4, &
                            ghost_z0(1,1),ld*ny,MPI_RPREC,down,tag_counter+4, &
                            comm,status,ierr)
$endif

    if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
      scalar_z(1:nx,1:ny,1)=0._rprec
      !scalar_z(1:nx,1:ny,2)=(scalar(1:nx,1:ny,2)+scalar(1:nx,1:ny,1))/2._rprec
    else
      DO j=1,ny
       DO i=1,nx
        ! k=1
        ! Choose the appropriate upwind direction
        IF (w_int(i,j,1) >= 0._rprec) THEN
           IF ((scalar(i,j,$lbz)-ghost_z0(i,j)==0._rprec)) THEN
             scalar_z(i,j,1)=scalar(i,j,$lbz)
           ELSE
             ! Gradient ratio
             rz=(scalar(i,j,1)-scalar(i,j,$lbz))/(scalar(i,j,$lbz)-ghost_z0(i,j))
             ! Bound the gradient
             rz=MAX(0._rprec,MIN(MIN(2._rprec*rz,0.75_rprec*rz+0.25_rprec),4._rprec))
             ! Interpolate the value
             scalar_z(i,j,1)=scalar(i,j,$lbz)+rz*(scalar(i,j,$lbz)-ghost_z0(i,j))/2._rprec
           ENDIF

        ELSE
          ! To avoid problems of divisions by zero
          IF ((scalar(i,j,1)-scalar(i,j,2))==0._rprec) THEN
            scalar_z(i,j,1)=scalar(i,j,1)
          ELSE
            ! Gradient ratio
            rz=(scalar(i,j,$lbz)-scalar(i,j,1))/(scalar(i,j,1)-scalar(i,j,2))
            ! Bound the gradient
            rz=MAX(0._rprec,MIN(MIN(2._rprec*rz,0.75_rprec*rz+0.25_rprec),4._rprec))
            ! Intrepolate the value
            scalar_z(i,j,1)=scalar(i,j,1)+rz*(scalar(i,j,1)-scalar(i,j,2))/2._rprec
          END IF
        END IF
       END DO
      END DO
    endif
$if ($MPI)
    call mpi_sendrecv(scalar_z(1,1,1),(nx+1)*(ny+1),MPI_RPREC,down,tag_counter+1, &
                     scalar_z(1,1,nz),(nx+1)*(ny+1),MPI_RPREC,up,tag_counter+1, &
                     comm,status,ierr)
$endif
    if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
    scalar_z(1:nx,1:ny,nz)=scalar_z(1:nx,1:ny,nz-1)
    endif

    ! Compute face nodes for periodic boundary conditions
    IF (periodicbcx) THEN
    
      ! Interpolated concentrations in x
      DO k=1,nz-1
        DO j=1,ny
                
          ! First node (i=1)
          ! Choose the appropriate upwind direction
          IF (u_int(1,j,k) >= 0._rprec) THEN  
            ! To avoid problems of divisions by zero
            IF ((scalar(nx,j,k)-scalar(nx-1,j,k))==0._rprec) THEN
              scalar_x(1,j,k)=scalar(nx,j,k)
            ELSE
              ! Gradient ratio
              rx=(scalar(1,j,k)-scalar(nx,j,k))/(scalar(nx,j,k)-scalar(nx-1,j,k))
              ! Bound the gradient
              rx=MAX(0._rprec,MIN(MIN(2._rprec*rx,0.75_rprec*rx+0.25_rprec),4._rprec))
              ! Intrepolate the value
              scalar_x(1,j,k)=scalar(nx,j,k)+rx*(scalar(nx,j,k)-scalar(nx-1,j,k))/2._rprec
            END IF
          ELSE      
            ! To avoid problems of divisions by zero
            IF ((scalar(1,j,k)-scalar(2,j,k))==0._rprec) THEN
              scalar_x(1,j,k)=scalar(1,j,k)
            ELSE
              ! Gradient ratio
              rx=(scalar(nx,j,k)-scalar(1,j,k))/(scalar(1,j,k)-scalar(2,j,k))
              ! Bound the gradient
              rx=MAX(0._rprec,MIN(MIN(2._rprec*rx,0.75_rprec*rx+0.25_rprec),4._rprec))
              ! Intrepolate the value
              scalar_x(1,j,k)=scalar(1,j,k)+rx*(scalar(1,j,k)-scalar(2,j,k))/2._rprec
            END IF
          END IF

          ! Last node (i=nx+1) - this is the same as the first node!
          scalar_x(nx+1,j,k)=scalar_x(1,j,k)
          
        END DO        
      END DO

    END IF      


    IF (periodicbcy) THEN
      
      ! Interpolated concentrations in y
      DO k=1,nz-1
        DO i=1,nx
                
          ! First node (j=1)
          ! Choose the appropriate upwind direction
          IF (v_int(i,1,k) >= 0._rprec) THEN  
            ! To avoid problems of divisions by zero
            IF ((scalar(i,ny,k)-scalar(i,ny-1,k))==0._rprec) THEN
              scalar_y(i,1,k)=scalar(i,ny,k)
            ELSE
              ! Gradient ratio
              ry=(scalar(i,1,k)-scalar(i,ny,k))/(scalar(i,ny,k)-scalar(i,ny-1,k))
              ! Bound the gradient
              ry=MAX(0._rprec,MIN(MIN(2._rprec*ry,0.75_rprec*ry+0.25_rprec),4._rprec))
              ! Intrepolate the value
              scalar_y(i,1,k)=scalar(i,ny,k)+ry*(scalar(i,ny,k)-scalar(i,ny-1,k))/2._rprec
            END IF
          ELSE      
            ! To avoid problems of divisions by zero
            IF ((scalar(i,1,k)-scalar(i,2,k))==0._rprec) THEN
              scalar_y(i,1,k)=scalar(i,1,k)
            ELSE
              ! Gradient ratio
              ry=(scalar(i,ny,k)-scalar(i,1,k))/(scalar(i,1,k)-scalar(i,2,k))
              ! Bound the gradient
              ry=MAX(0._rprec,MIN(MIN(2._rprec*ry,0.75_rprec*ry+0.25_rprec),4._rprec))
              ! Intrepolate the value
              scalar_y(i,1,k)=scalar(i,1,k)+ry*(scalar(i,1,k)-scalar(i,2,k))/2._rprec
            END IF
          END IF

          ! Last node (j=ny+1) - this is the same as the first node!
          scalar_y(i,ny+1,k)=scalar_y(i,1,k)
          
        END DO        
      END DO


    END IF
    

  END IF
  
  
  ! In these 3 cases assign Kc_t=Nu_t/Sc (Sc from param.f90)
  ! Case 1 - initialization (jt<cs_count)
  ! Case 2 - No dynamic Sc (model_sc==1)
  ! Case 3 - Initial part of the run (jt<DYN_init)
  IF (jt<cs_count .OR. model_sc==1 .OR. jt<DYN_init) THEN


    ! This is the new interpolation for nu_t - center of finite volumes
    ! (interpolation in z)
    
    ! The first level does not need any interpolation
    if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    Kc_t(1:nx,1:ny,1)=Nu_t(1:nx,1:ny,1)
    else
    Kc_t(1:nx,1:ny,1)=(Nu_t(1:nx,1:ny,1)+Nu_t(1:nx,1:ny,2))/2._rprec
    end if
 
    ! From now on, interpolates in z
    Kc_t(1:nx,1:ny,2:nz-2)=(Nu_t(1:nx,1:ny,2:nz-2)+Nu_t(1:nx,1:ny,3:nz-1))/2._rprec

    ! The top one is only copied
    if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
    Kc_t(1:nx,1:ny,nz-1)=Nu_t(1:nx,1:ny,nz-1)
    else
    Kc_t(1:nx,1:ny,nz-1)=(Nu_t(1:nx,1:ny,nz-1)+Nu_t(1:nx,1:ny,nz))/2._rprec
    end if

    ! Divide by constant Schmidt number
    Kc_t=Kc_t/Sc
                

  ELSE
  
    ! Calculate dynamic coefficient Cs2/Sc
    IF (mod(jt,cs_count)==0) THEN
    
      ! Plane Averaged Scale Invariant (PASI);
      IF (model_sc==2) THEN
	CALL dynamic_sc_pasi(scalar,scalar_x,scalar_y,scalar_z,u_int,v_int,w_int)
      ! Plane Averaged Scale Dependent (PASD);
      ELSE IF (model_sc==3) THEN
	CALL dynamic_sc_pasd(scalar,scalar_x,scalar_y,scalar_z,u_int,v_int,w_int)
      ! Lagrangian Averaged Scale Invariant (LASI);
      ELSE IF (model_sc==4) THEN
	CALL dynamic_sc_lasi(scalar,scalar_x,scalar_y,scalar_z,u_int,v_int,w_int)
      ! Lagragian Averaged Scale Dependent (LASD);
      ELSE IF (model_sc==5) THEN
	CALL dynamic_sc_lasd(scalar,scalar_x,scalar_y,scalar_z,u_int,v_int,w_int)
      END IF

    END IF

    ! Filter size
    delta=(dx*dy*dz)**(1._rprec/3._rprec)

    ! Calculate eddy-diffusivity
    Kc_t=Cs2Sc*(delta**2)*magS_int

  END IF
  
$if ($MPI)
  call mpi_sendrecv (Kc_t(1, 1, 1), ld*ny, MPI_RPREC, down, tag_counter+9,  &
                     Kc_t(1, 1, nz), ld*ny, MPI_RPREC, up, tag_counter+9,   &
                     comm, status, ierr)
  call mpi_sendrecv (Kc_t(1, 1, nz-1), ld*ny, MPI_RPREC, up, tag_counter+10,  &
                     Kc_t(1, 1, 0), ld*ny, MPI_RPREC, down, tag_counter+10,   &
                     comm, status, ierr)
$endif

  !
  ! Step 4 - Assemble RHS
  !
  

  ! Now it is just standard centered finite differences
  ! Add the terms one by one...

  ! 1.1) Advection => -u*dC/dx
  DO k=1,nz-1
    DO j=1,ny
      DO i=1,nx
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Modified by Di Yang for ocean LES with Stokes drift
!DY Start here
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY  	RHS(i,j,k)=-(1._rprec/dx)*(u_int(i+1,j,k)*scalar_x(i+1,j,k)-u_int(i,j,k)*scalar_x(i,j,k))
         if(OCEAN_FLAG .and. STOKES_FLAG) then
            RHS(i,j,k)=-(1._rprec/dx)*(u_int(i+1,j,k)*scalar_x(i+1,j,k)-u_int(i,j,k)*scalar_x(i,j,k)&
                 +ust(k)*(scalar_x(i+1,j,k)-scalar_x(i,j,k)))
         else
            RHS(i,j,k)=-(1._rprec/dx)*(u_int(i+1,j,k)*scalar_x(i+1,j,k)-u_int(i,j,k)*scalar_x(i,j,k))
         endif
!+++++++++++++++
!DY End here
!+++++++++++++++
      END DO
    END DO
  END DO

!!$  ! 1.1) Advection => -u*dC/dx
!!$  DO k=1,nz-1
!!$    DO j=1,ny
!!$      DO i=1,nx
!!$  	RHS(i,j,k)=-(1._rprec/dx)*(u_int(i+1,j,k)*scalar_x(i+1,j,k)-u_int(i,j,k)*scalar_x(i,j,k))
!!$      END DO
!!$    END DO
!!$  END DO
  
  ! 1.2) Advection => -v*dC/dy
  DO k=1,nz-1
    DO j=1,ny
      DO i=1,nx
  	RHS(i,j,k)=RHS(i,j,k)-(1._rprec/dy)*(v_int(i,j+1,k)*scalar_y(i,j+1,k)-v_int(i,j,k)*scalar_y(i,j,k))
      END DO
    END DO
  END DO

  ! 1.3) Advection => -w*dC/dz
  DO k=1,nz-1
    DO j=1,ny
      DO i=1,nx
  	RHS(i,j,k)=RHS(i,j,k)-(1._rprec/dz)*(w_int(i,j,k+1)*scalar_z(i,j,k+1)-w_int(i,j,k)*scalar_z(i,j,k))
      END DO
    END DO
  END DO 
  ! Store SGS vertical flux
  res_vert(1:nx,1:ny,1:nz)=w_int(1:nx,1:ny,1:nz)*scalar_z(1:nx,1:ny,1:nz)
  
  ! 2.1) SGS diffusion => d((nu/Sc)*dc/dx)/dx
  DO k=1,nz-1
    DO j=1,ny
      DO i=2,nx-1
  	RHS(i,j,k)=RHS(i,j,k)+(1._rprec/dx**2)*(((Kc_t(i,j,k)+Kc_t(i+1,j,k))/2._rprec)*(scalar(i+1,j,k)-scalar(i,j,k)) &
	                                       -((Kc_t(i,j,k)+Kc_t(i-1,j,k))/2._rprec)*(scalar(i,j,k)-scalar(i-1,j,k)))
      END DO
      RHS(1,j,k)=RHS(1,j,k)+(1._rprec/dx**2)*(((Kc_t(1,j,k)+Kc_t(2,j,k))/2._rprec)*(scalar(2,j,k)-scalar(1,j,k)) &
	                                     -((Kc_t(1,j,k)+Kc_t(nx,j,k))/2._rprec)*(scalar(1,j,k)-ghost_x0(j,k)))
      RHS(nx,j,k)=RHS(nx,j,k)+(1._rprec/dx**2)*(((Kc_t(nx,j,k)+Kc_t(1,j,k))/2._rprec)*(ghost_xLx(j,k)-scalar(nx,j,k)) &
	                                       -((Kc_t(nx,j,k)+Kc_t(nx-1,j,k))/2._rprec)*(scalar(nx,j,k)-scalar(nx-1,j,k)))
    END DO
  END DO

  ! 2.2) SGS diffusion => d((nu/Sc)*dc/dy)/dy
  DO k=1,nz-1
    DO i=1,nx
      DO j=2,ny-1
  	RHS(i,j,k)=RHS(i,j,k)+(1._rprec/dy**2)*(((Kc_t(i,j,k)+Kc_t(i,j+1,k))/2._rprec)*(scalar(i,j+1,k)-scalar(i,j,k)) &
	                                       -((Kc_t(i,j,k)+Kc_t(i,j-1,k))/2._rprec)*(scalar(i,j,k)-scalar(i,j-1,k)))
      END DO
      RHS(i,1,k)=RHS(i,1,k)+(1._rprec/dy**2)*(((Kc_t(i,1,k)+Kc_t(i,2,k))/2._rprec)*(scalar(i,2,k)-scalar(i,1,k)) &
	                                     -((Kc_t(i,1,k)+Kc_t(i,ny,k))/2._rprec)*(scalar(i,1,k)-ghost_y0(i,k)))
      RHS(i,ny,k)=RHS(i,ny,k)+(1._rprec/dy**2)*(((Kc_t(i,ny,k)+Kc_t(i,1,k))/2._rprec)*(ghost_yLy(i,k)-scalar(i,ny,k)) &
	                                       -((Kc_t(i,ny,k)+Kc_t(i,ny-1,k))/2._rprec)*(scalar(i,ny,k)-scalar(i,ny-1,k)))
    END DO
  END DO

  ! 2.3) SGS diffusion => d((nu/Sc)*dc/dz)/dz
  DO j=1,ny
    DO i=1,nx
      DO k=2,nz-2

       RHS(i,j,k)=RHS(i,j,k)+(1._rprec/dz**2)*(((Kc_t(i,j,k+1)+Kc_t(i,j,k))/2._rprec)*(scalar(i,j,k+1)-scalar(i,j,k)) &
					      -((Kc_t(i,j,k)+Kc_t(i,j,k-1))/2._rprec)*(scalar(i,j,k)-scalar(i,j,k-1)))

       ! Store SGS vertical flux
       sgs_vert(i,j,k)=-(1._rprec/dz)*((Kc_t(i,j,k)+Kc_t(i,j,k-1))/2._rprec)*(scalar(i,j,k)-scalar(i,j,k-1))

      END DO

      if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
      ! Impose surface flux at bottom of element 1
      RHS(i,j,1)=RHS(i,j,1)+(1._rprec/dz**2)*(((Kc_t(i,j,2)+Kc_t(i,j,1))/2._rprec)*(scalar(i,j,2)-scalar(i,j,1))) &
                                            +(1._rprec/dz)*(P_surf_flux(i,j))
      ! Store surface flux
      sgs_vert(i,j,1)=P_surf_flux(i,j)
      else
      RHS(i,j,1)=RHS(i,j,1)+(1._rprec/dz**2)*(((Kc_t(i,j,2)+Kc_t(i,j,1))/2._rprec)*(scalar(i,j,2)-scalar(i,j,1)) &
                                              -((Kc_t(i,j,1)+Kc_t(i,j,$lbz))/2._rprec)*(scalar(i,j,1)-scalar(i,j,$lbz)))
      sgs_vert(i,j,1)=-(1._rprec/dz)*((Kc_t(i,j,1)+Kc_t(i,j,$lbz))/2._rprec)*(scalar(i,j,1)-scalar(i,j,$lbz))
      end if

      if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
      ! Impose zero flux at top of element nz-1
      RHS(i,j,nz-1)=RHS(i,j,nz-1)+(1._rprec/dz**2)*(0._rprec &
                                     -((Kc_t(i,j,nz-1)+Kc_t(i,j,nz-2))/2._rprec)*(scalar(i,j,nz-1)-scalar(i,j,nz-2)))
      else
      RHS(i,j,nz-1)=RHS(i,j,nz-1)+(1._rprec/dz**2)*(((Kc_t(i,j,nz)+Kc_t(i,j,nz-1))/2._rprec)*(scalar(i,j,nz)-scalar(i,j,nz-1)) &
                                              -((Kc_t(i,j,nz-1)+Kc_t(i,j,nz-2))/2._rprec)*(scalar(i,j,nz-1)-scalar(i,j,nz-2)))
      end if

      sgs_vert(i,j,nz-1)=-(1._rprec/dz)*((Kc_t(i,j,nz-1)+Kc_t(i,j,nz-2))/2._rprec)*(scalar(i,j,nz-1)-scalar(i,j,nz-2))

      DO k=1,nz
       ! Store SGS x-direction flux
       sgs_x(i,j,k)=-(1._rprec/dx)*((Kc_t(i+1,j,k)+Kc_t(i,j,k)/2._rprec)*(scalar(i+1,j,k)-scalar(i,j,k)))
      END DO

    END DO
  END DO
  
 
  ! Add source term if needed
  IF (pointsource .AND. (jt_total>=ini_src .AND. jt_total<end_src)) THEN

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Block 1 starts here
!DY Changed by Di Yang
!DY The original version add source term at jz=zps in the processor coord = 0.
!DY This is problematic when we use a large number of processors.
!DY Under certain conditions, the point source will not be in the first processor, 
!DY and zps will be larger than nz, which causes an "array bound" bug.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang, start here
     ncpu_source = int(zps/nz)
     zps_local = zps-ncpu_source*nz+1

     if ((.not. USE_MPI) .or. (USE_MPI .and. coord == ncpu_source)) then

        IF (single_point) THEN

           RHS(xps,yps,zps_local)=RHS(xps,yps,zps_local)+(0.75_rprec)*Q_src
           RHS(xps,yps,zps_local+1)=RHS(xps,yps,zps_local+1)+(0.25_rprec)*Q_src

        ELSE
           
           DO j=1,ny
              DO i = 1,nx
                 RHS(i,j,zps_local)=RHS(i,j,zps_local)+Q_src
              END DO
           END DO

        END IF

     end if
!DY Added by Di Yang, end here

!DY     if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then

!DY        IF (single_point) THEN      
    ! Source in 1 point
!    RHS(xps,yps,zps)=RHS(xps,yps,zps)+Q_src

!DY           RHS(xps,yps,zps)=RHS(xps,yps,zps)+(0.75_rprec)*Q_src
!DY           RHS(xps,yps,zps+1)=RHS(xps,yps,zps+1)+(0.25_rprec)*Q_src

!DY        ELSE

!DY           DO j=1,ny
!DY              DO i = 1,nx
!DY                 RHS(i,j,zps)=RHS(i,j,zps)+Q_src
!DY              END DO
!DY           END DO

!DY        END IF

!DY     end if

!++++++++++++++++++++++++
!DY Block 1 end here
!++++++++++++++++++++++++

  END IF
  

  
  ! Calculate total flux at x=Lx for mass balance
  ! Save previous value for time integration
  flux_out_prev=flux_out

  flux_out=(SUM(u_int(nx+1,1:ny,1:nz-1)*scalar_x(nx+1,1:ny,1:nz-1))*dt*dy*dz)
  flux_out=flux_out+(SUM(v_int(1:nx,1,1:nz-1)*scalar_x(1:nx,1,1:nz-1))*dt*dx*dz)

END SUBROUTINE pollen_RHS_calc2






! This is a simplified version of step_scalar()
! Chamecki 10/04/2006
subroutine step_pollen2(scalar,RHS_pre,RHS_post)
implicit none
integer:: i,j,k
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

  real(kind=rprec),dimension(ld,ny,$lbz:nz)::scalar, RHS_pre, RHS_post

  do k=1,nz-1
    do j=1,ny
      do i=1,nx
        scalar(i,j,k)= scalar(i,j,k)+dt*(1.5_rprec*RHS_pre(i,j,k)-0.5_rprec*RHS_post(i,j,k))
      end do
    end do
  end do     

  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
  scalar(:,:,nz)=scalar(:,:,nz-1)
  end if

end subroutine step_pollen2





end module scalars_module
