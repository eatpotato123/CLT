! This CLT program is upgraded form the old version clt of swang, 
! 4th order finite difference method is employed in the R and Z directions (upgraded by wzhang), 
! while in the phi(y) direction, either finite difference or pseudo-spectrum method is used,
! in the time-advance, 4th order Runge-Kutta scheme is chosen.
! Ref: Physics of Plasmas 22, 122504 (2015); doi: 10.1063/1.4936977

! recently added features (upgrader by hwzhang):
! 1. cut-cell method is used in the boundary,
! 2. can be used for any triangularity like East or the Circle case. Fixed boundary is used,
! 3. can be used for any nprx and nprz ( grids and nprx(z) should be divisible),
! 4. can be used for including SOL region from EFIT-gfile equilibriums,
! 5. subroutines for RMP-EAST coils is added.

! main parameters and notations
! mmode=2
! hall=flase
! fdiout=5.e-2
! mxt=256,myt=32,mzt=256
! cfl=1.2
! kap0=5.e-5
! x(1-8):total of rho,p,vx,vy(v_phi),vz,bx,by(b_phi),bz
! x1(1-8):pertubation of rho,p,vx,vy(v_phi),vz,bx,by(b_phi),bz
! xint(1-8):equilibrium of rho,p,vx,vy(v_phi),vz,bx,by(b_phi),bz
! cur(1-3):pertubation of current _(x,y(phi),z)
! cint(1-3):equilibrium of current _(x,y(phi),z)
! ef(1-3):pertubation of e-field _(x,y(phi),z)

! ****coordinates***
! xx(mx),yy(my),zz(mz): coordinates(r,phi,z) in each processes; 
! xxt(mxt),yyt(myt),zzt(mzt) for total; 
! xxst(n2th+5,npsi),zzst(n2th+5,npsi) in (theta,psi) grid; 
! xxs(n2th+5,mps4:mps),zzs(n2th+5,mps4:mps) in (theta,psi) bandary grid;

! thxz(mx,mz): theta coordinates in (r,z) grid; tht(mxt,mzt) for total; 
! tpxz(mx,mz): r.z.<->s(psi).p(pol). transit angle in (r,z); tpt(mxt,mzt) for total;
! tcxz(mx,mz): tc=ing(jcb/r^2)dth 

! thst(n2th+5): theta coordinates in (theta,psi) grid;
! tpst(n2th+5,npsi): r.z.<->s(psi).p(pol). transit angle in (theta,psi); tps(n2th+5,mps4:mps) for bndry;
! usually,
! th=arc/n2th
! tp=atan2(psdz,psdx)=atan2(bx,-bz) => cos(tp)=psdx/|grps|=-bz/bp;sin(tp)=psdz/|grps|=bx/bp;






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!   !!!!!!!!!!   !!!!!!!      !!!!!!    !!!!!!!!!!!
!!!!!!!!!!!    !!!!!!!  !  !!!!!!!  !!!!  !!!!!  !!!!!!!!!!!!
!!!!!!!!!!!  !  !!!!!  !!  !!!!!!!  !!!!! !!!!!  !!!!!!!!!!!!
!!!!!!!!!!!  !!!  !! !!!!  !!!!!!!  !!!! !!!!!!  !!!!!!!!!!!!
!!!!!!!!!!!  !!!!!  !!!!!  !!!!!!!      !!!!!!!  !!!!!!!!!!!!
!!!!!!!!!!!  !!!!!!!!!!!!  !!!!!!!  !!!!!!!!!!!  !!!!!!!!!!!!
!!!!!!!!!!!  !!!!!!!!!!!!  !!!!!!!  !!!!!!!!!!    !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module for MPI transfer subroutines.
! written by S. WANG
! sortted by H.W. Zhang during 2018 AUG.
! contact changhw@zju.edu.cn if there is any problem.

! including subroutines as follows
! mpi_transfersm(ws,mm)
! mpi_transfer8(w8)
! mpi_transfer3(w3)
! mpi_transfer1(w1)
! mpi_transfersy1(wy1,wyt)
! mpi_transfersm_2d(ws,mm)
! mpi_transfer8_2d(w8)
! mpi_transfer3_2d(w3)
! mpi_transfer1_2d(w1)
! mpi_transfersm_one_layer(ws,mm)



!ws*************************************************
      subroutine mpi_transfersm(ws,mm)
      use declare
      real*8, dimension(mx,mz,my,mm) :: ws
      real*8, dimension(mz,my,mm) :: wsx1,wsx2
      real*8, dimension(mx,my,mm) :: wsz1,wsz2
      real*8, dimension(mx,mz,mm) :: wsy1,wsy2
      include 'mpif.h'

!       
! send w8 up unless i'm at the top, then receive from below
     

	if (nrankxz.ge.nrkz(nrank)*nprx .and. nrankxz.lt.nrkz(nrank)*nprx+nprx-1) then
      wsx1(:,:,:)=ws(ix_last-2,:,:,:)
      wsx2(:,:,:)=ws(ix_last-3,:,:,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wsx1, myz*mm, mpi_double_precision, nrank + 1, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wsx2, myz*mm, mpi_double_precision, nrank + 1, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrankxz.gt.nrkz(nrank)*nprx .and. nrankxz.le.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wsx1, myz*mm, mpi_double_precision, nrank - 1, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wsx2, myz*mm, mpi_double_precision, nrank - 1, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      ws(ix_first+1,:,:,:)=wsx1(:,:,:)
      ws(ix_first,:,:,:)=wsx2(:,:,:)
	endif
	
      
! send w8 down unless i'm at the bottom

      
	if (nrankxz.gt.nrkz(nrank)*nprx .and. nrankxz.le.nrkz(nrank)*nprx+nprx-1) then
      wsx1(:,:,:)=ws(ix_first+2,:,:,:)
      wsx2(:,:,:)=ws(ix_first+3,:,:,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wsx1, myz*mm, mpi_double_precision, nrank - 1, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wsx2, myz*mm, mpi_double_precision, nrank - 1, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrankxz.ge.nrkz(nrank)*nprx .and. nrankxz.lt.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wsx1, myz*mm, mpi_double_precision, nrank + 1, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wsx2, myz*mm, mpi_double_precision, nrank + 1, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      ws(ix_last-1,:,:,:)=wsx1(:,:,:)
      ws(ix_last,:,:,:)=wsx2(:,:,:)
	endif
!
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	if (nrankxz.lt.nsizexz-nprx) then
      wsz1(:,:,:)=ws(:,iz_last-2,:,:)
      wsz2(:,:,:)=ws(:,iz_last-3,:,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wsz1, myx*mm, mpi_double_precision, nrank + nprx, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wsz2, myx*mm, mpi_double_precision, nrank + nprx, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrankxz.ge.nprx ) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wsz1, myx*mm, mpi_double_precision, nrank - nprx, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wsz2, myx*mm, mpi_double_precision, nrank - nprx, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      ws(:,iz_first+1,:,:)=wsz1(:,:,:)
      ws(:,iz_first,:,:)=wsz2(:,:,:)
	endif
	     
! send w8 down unless i'm at the bottom

      
	if (nrankxz.ge.nprx ) then
      wsz1(:,:,:)=ws(:,iz_first+2,:,:)
      wsz2(:,:,:)=ws(:,iz_first+3,:,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wsz1, myx*mm, mpi_double_precision, nrank - nprx, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wsz2, myx*mm, mpi_double_precision, nrank - nprx, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrankxz.lt.nsizexz-nprx) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wsz1, myx*mm, mpi_double_precision, nrank + nprx, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wsz2, myx*mm, mpi_double_precision, nrank + nprx, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      ws(:,iz_last-1,:,:)=wsz1(:,:,:)
      ws(:,iz_last,:,:)=wsz2(:,:,:)
	endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        if(npry .gt. 1) then
	if (nrky(nrank).lt. npry-1) then
      wsy1(:,:,:)=ws(:,:,iy_last-2,:)
      wsy2(:,:,:)=ws(:,:,iy_last-3,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wsy1, mxz*mm, mpi_double_precision, nrank + nprxz, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wsy2, mxz*mm, mpi_double_precision, nrank + nprxz, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrky(nrank) .ge. 1 )  then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wsy1, mxz*mm, mpi_double_precision, nrank - nprxz, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wsy2, mxz*mm, mpi_double_precision, nrank - nprxz, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      ws(:,:,iy_first+1,:)=wsy1(:,:,:)
      ws(:,:,iy_first,:)=wsy2(:,:,:)
	endif


  	if (nrky(nrank).eq. npry-1) then
      wsy1(:,:,:)=ws(:,:,iy_last-2,:)
      wsy2(:,:,:)=ws(:,:,iy_last-3,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wsy1, mxz*mm, mpi_double_precision, nrank -(npry-1)*nprxz, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wsy2, mxz*mm, mpi_double_precision, nrank -(npry-1)*nprxz, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrky(nrank) .eq. 0 )  then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wsy1, mxz*mm, mpi_double_precision, nrank +(npry-1)*nprxz, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wsy2, mxz*mm, mpi_double_precision, nrank +(npry-1)*nprxz, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      ws(:,:,iy_first+1,:)=wsy1(:,:,:)
      ws(:,:,iy_first,:)=wsy2(:,:,:)
	endif
	     
! send w8 down unless i'm at the bottom

      
	if (nrky(nrank) .ge. 1 ) then
      wsy1(:,:,:)=ws(:,:,iy_first+2,:)
      wsy2(:,:,:)=ws(:,:,iy_first+3,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wsy1, mxz*mm, mpi_double_precision, nrank - nprxz, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wsy2, mxz*mm, mpi_double_precision, nrank - nprxz, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrky(nrank).lt. npry-1) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wsy1, mxz*mm, mpi_double_precision, nrank + nprxz, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wsy2, mxz*mm, mpi_double_precision, nrank + nprxz, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      ws(:,:,iy_last-1,:)=wsy1(:,:,:)
      ws(:,:,iy_last,:)=wsy2(:,:,:)
	endif

    if (nrky(nrank) .eq. 0 ) then
      wsy1(:,:,:)=ws(:,:,iy_first+2,:)
      wsy2(:,:,:)=ws(:,:,iy_first+3,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wsy1, mxz*mm, mpi_double_precision, nrank +(npry-1)*nprxz, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wsy2, mxz*mm, mpi_double_precision, nrank +(npry-1)*nprxz, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrky(nrank).eq. npry-1) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wsy1, mxz*mm, mpi_double_precision, nrank -(npry-1)*nprxz, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wsy2, mxz*mm, mpi_double_precision, nrank -(npry-1)*nprxz, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      ws(:,:,iy_last-1,:)=wsy1(:,:,:)
      ws(:,:,iy_last,:)=wsy2(:,:,:)
	endif
      else
      ws(:,:,iy_first+1,:)=ws(:,:,iy_last-2,:)
      ws(:,:,iy_first,:)=ws(:,:,iy_last-3,:)
      ws(:,:,iy_last-1,:)=ws(:,:,iy_first+2,:)
      ws(:,:,iy_last,:)=ws(:,:,iy_first+3,:)
      endif

      return
      end 

!ws****************************************************************************************
      subroutine mpi_transfer8(w8)
      use declare
      real*8, dimension(mx,mz,my,8) :: w8
      real*8, dimension(mz,my,8) :: wfx1,wfx2
      real*8, dimension(mx,my,8) :: wfz1,wfz2
      real*8, dimension(mx,mz,8) :: wfy1,wfy2
      include 'mpif.h'

!       
! send w8 up unless i'm at the top, then receive from below
     
	if (nrankxz.ge.nrkz(nrank)*nprx .and. nrankxz.lt.nrkz(nrank)*nprx+nprx-1) then
      wfx1(:,:,:)=w8(ix_last-2,:,:,:)
      wfx2(:,:,:)=w8(ix_last-3,:,:,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wfx1, myz8, mpi_double_precision, nrank + 1, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wfx2, myz8, mpi_double_precision, nrank + 1, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrankxz.gt.nrkz(nrank)*nprx .and. nrankxz.le.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wfx1, myz8, mpi_double_precision, nrank - 1, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wfx2, myz8, mpi_double_precision, nrank - 1, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w8(ix_first+1,:,:,:)=wfx1(:,:,:)
      w8(ix_first,:,:,:)=wfx2(:,:,:)
	endif
	
      
! send w8 down unless i'm at the bottom

      
	if (nrankxz.gt.nrkz(nrank)*nprx .and. nrankxz.le.nrkz(nrank)*nprx+nprx-1) then
      wfx1(:,:,:)=w8(ix_first+2,:,:,:)
      wfx2(:,:,:)=w8(ix_first+3,:,:,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wfx1, myz8, mpi_double_precision, nrank - 1, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wfx2, myz8, mpi_double_precision, nrank - 1, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrankxz.ge.nrkz(nrank)*nprx .and. nrankxz.lt.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wfx1, myz8, mpi_double_precision, nrank + 1, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wfx2, myz8, mpi_double_precision, nrank + 1, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w8(ix_last-1,:,:,:)=wfx1(:,:,:)
      w8(ix_last,:,:,:)=wfx2(:,:,:)
	endif
!
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	if (nrankxz.lt.nsizexz-nprx) then
      wfz1(:,:,:)=w8(:,iz_last-2,:,:)
      wfz2(:,:,:)=w8(:,iz_last-3,:,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wfz1, myx8, mpi_double_precision, nrank + nprx, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wfz2, myx8, mpi_double_precision, nrank + nprx, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrankxz.ge.nprx ) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wfz1, myx8, mpi_double_precision, nrank - nprx, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wfz2, myx8, mpi_double_precision, nrank - nprx, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w8(:,iz_first+1,:,:)=wfz1(:,:,:)
      w8(:,iz_first,:,:)=wfz2(:,:,:)
	endif
	     
! send w8 down unless i'm at the bottom

      
	if (nrankxz.ge.nprx ) then
      wfz1(:,:,:)=w8(:,iz_first+2,:,:)
      wfz2(:,:,:)=w8(:,iz_first+3,:,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wfz1, myx8, mpi_double_precision, nrank - nprx, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wfz2, myx8, mpi_double_precision, nrank - nprx, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrankxz.lt.nsizexz-nprx) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wfz1, myx8, mpi_double_precision, nrank + nprx, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wfz2, myx8, mpi_double_precision, nrank + nprx, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w8(:,iz_last-1,:,:)=wfz1(:,:,:)
      w8(:,iz_last,:,:)=wfz2(:,:,:)
	endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      if (npry .gt. 1) then
	if (nrky(nrank).lt. npry-1) then
      wfy1(:,:,:)=w8(:,:,iy_last-2,:)
      wfy2(:,:,:)=w8(:,:,iy_last-3,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wfy1, mxz8, mpi_double_precision, nrank + nprxz, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wfy2, mxz8, mpi_double_precision, nrank + nprxz, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrky(nrank) .ge. 1 )  then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wfy1, mxz8, mpi_double_precision, nrank - nprxz, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wfy2, mxz8, mpi_double_precision, nrank - nprxz, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w8(:,:,iy_first+1,:)=wfy1(:,:,:)
      w8(:,:,iy_first,:)=wfy2(:,:,:)
	endif


  	if (nrky(nrank).eq. npry-1) then
      wfy1(:,:,:)=w8(:,:,iy_last-2,:)
      wfy2(:,:,:)=w8(:,:,iy_last-3,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wfy1, mxz8, mpi_double_precision, nrank -(npry-1)*nprxz, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wfy2, mxz8, mpi_double_precision, nrank -(npry-1)*nprxz, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrky(nrank) .eq. 0 )  then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wfy1, mxz8, mpi_double_precision, nrank +(npry-1)*nprxz, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wfy2, mxz8, mpi_double_precision, nrank +(npry-1)*nprxz, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w8(:,:,iy_first+1,:)=wfy1(:,:,:)
      w8(:,:,iy_first,:)=wfy2(:,:,:)
	endif
	     
! send w8 down unless i'm at the bottom

      
	if (nrky(nrank) .ge. 1 ) then
      wfy1(:,:,:)=w8(:,:,iy_first+2,:)
      wfy2(:,:,:)=w8(:,:,iy_first+3,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wfy1, mxz8, mpi_double_precision, nrank - nprxz, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wfy2, mxz8, mpi_double_precision, nrank - nprxz, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrky(nrank).lt. npry-1) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wfy1, mxz8, mpi_double_precision, nrank + nprxz, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wfy2, mxz8, mpi_double_precision, nrank + nprxz, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w8(:,:,iy_last-1,:)=wfy1(:,:,:)
      w8(:,:,iy_last,:)=wfy2(:,:,:)
	endif

    if (nrky(nrank) .eq. 0 ) then
      wfy1(:,:,:)=w8(:,:,iy_first+2,:)
      wfy2(:,:,:)=w8(:,:,iy_first+3,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wfy1, mxz8, mpi_double_precision, nrank +(npry-1)*nprxz, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wfy2, mxz8, mpi_double_precision, nrank +(npry-1)*nprxz, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrky(nrank).eq. npry-1) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wfy1, mxz8, mpi_double_precision, nrank -(npry-1)*nprxz, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wfy2, mxz8, mpi_double_precision, nrank -(npry-1)*nprxz, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w8(:,:,iy_last-1,:)=wfy1(:,:,:)
      w8(:,:,iy_last,:)=wfy2(:,:,:)
	endif
      else
      w8(:,:,iy_first+1,:)=w8(:,:,iy_last-2,:)
      w8(:,:,iy_first,:)=w8(:,:,iy_last-3,:)
      w8(:,:,iy_last-1,:)=w8(:,:,iy_first+2,:)
      w8(:,:,iy_last,:)=w8(:,:,iy_first+3,:)
      endif
      return
      end 

!ws****************************************************************************************
      subroutine mpi_transfer3(w3)
      use declare
      real*8, dimension(mx,mz,my,3) :: w3
      real*8, dimension(mz,my,3) :: wcx1,wcx2
      real*8, dimension(mx,my,3) :: wcz1,wcz2
      real*8, dimension(mx,mz,3) :: wcy1,wcy2
      include 'mpif.h'

       
! send w3 up unless i'm at the top, then receive from below
     
	if (nrankxz.ge.nrkz(nrank)*nprx .and. nrankxz.lt.nrkz(nrank)*nprx+nprx-1) then
      wcx1(:,:,:)=w3(ix_last-2,:,:,:)
      wcx2(:,:,:)=w3(ix_last-3,:,:,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wcx1, myz3, mpi_double_precision, nrank + 1, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wcx2, myz3, mpi_double_precision, nrank + 1, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrankxz.gt.nrkz(nrank)*nprx .and. nrankxz.le.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wcx1, myz3, mpi_double_precision, nrank - 1, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wcx2, myz3, mpi_double_precision, nrank - 1, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w3(ix_first+1,:,:,:)=wcx1(:,:,:)
      w3(ix_first,:,:,:)=wcx2(:,:,:)
	endif
	
      
! send w3 down unless i'm at the bottom

      
	if (nrankxz.gt.nrkz(nrank)*nprx .and. nrankxz.le.nrkz(nrank)*nprx+nprx-1) then
      wcx1(:,:,:)=w3(ix_first+2,:,:,:)
      wcx2(:,:,:)=w3(ix_first+3,:,:,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wcx1, myz3, mpi_double_precision, nrank - 1, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wcx2, myz3, mpi_double_precision, nrank - 1, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrankxz.ge.nrkz(nrank)*nprx .and. nrankxz.lt.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wcx1, myz3, mpi_double_precision, nrank + 1, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wcx2, myz3, mpi_double_precision, nrank + 1, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w3(ix_last-1,:,:,:)=wcx1(:,:,:)
      w3(ix_last,:,:,:)=wcx2(:,:,:)
	endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	if (nrankxz.lt.nsizexz-nprx) then
      wcz1(:,:,:)=w3(:,iz_last-2,:,:)
      wcz2(:,:,:)=w3(:,iz_last-3,:,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wcz1, myx3, mpi_double_precision, nrank + nprx, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wcz2, myx3, mpi_double_precision, nrank + nprx, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrankxz.ge.nprx ) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wcz1, myx3, mpi_double_precision, nrank - nprx, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wcz2, myx3, mpi_double_precision, nrank - nprx, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w3(:,iz_first+1,:,:)=wcz1(:,:,:)
      w3(:,iz_first,:,:)=wcz2(:,:,:)
	endif
	     
! send w3 down unless i'm at the bottom

      
	if (nrankxz.ge.nprx ) then
      wcz1(:,:,:)=w3(:,iz_first+2,:,:)
      wcz2(:,:,:)=w3(:,iz_first+3,:,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wcz1, myx3, mpi_double_precision, nrank - nprx, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wcz2, myx3, mpi_double_precision, nrank - nprx, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrankxz.lt.nsizexz-nprx) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wcz1, myx3, mpi_double_precision, nrank + nprx, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wcz2, myx3, mpi_double_precision, nrank + nprx, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w3(:,iz_last-1,:,:)=wcz1(:,:,:)
      w3(:,iz_last,:,:)=wcz2(:,:,:)
	endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      if(npry .gt. 1) then
	if (nrky(nrank).lt. npry-1) then
      wcy1(:,:,:)=w3(:,:,iy_last-2,:)
      wcy2(:,:,:)=w3(:,:,iy_last-3,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wcy1, mxz3, mpi_double_precision, nrank + nprxz, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wcy2, mxz3, mpi_double_precision, nrank + nprxz, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrky(nrank) .ge. 1 )  then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wcy1, mxz3, mpi_double_precision, nrank - nprxz, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wcy2, mxz3, mpi_double_precision, nrank - nprxz, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w3(:,:,iy_first+1,:)=wcy1(:,:,:)
      w3(:,:,iy_first,:)=wcy2(:,:,:)
	endif


  	if (nrky(nrank).eq. npry-1) then
      wcy1(:,:,:)=w3(:,:,iy_last-2,:)
      wcy2(:,:,:)=w3(:,:,iy_last-3,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wcy1, mxz3, mpi_double_precision, nrank -(npry-1)*nprxz, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wcy2, mxz3, mpi_double_precision, nrank -(npry-1)*nprxz, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrky(nrank) .eq. 0 )  then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wcy1, mxz3, mpi_double_precision, nrank +(npry-1)*nprxz, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wcy2, mxz3, mpi_double_precision, nrank +(npry-1)*nprxz, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w3(:,:,iy_first+1,:)=wcy1(:,:,:)
      w3(:,:,iy_first,:)=wcy2(:,:,:)
	endif
	     
! send w3 down unless i'm at the bottom

      
	if (nrky(nrank) .ge. 1 ) then
      wcy1(:,:,:)=w3(:,:,iy_first+2,:)
      wcy2(:,:,:)=w3(:,:,iy_first+3,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wcy1, mxz3, mpi_double_precision, nrank - nprxz, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wcy2, mxz3, mpi_double_precision, nrank - nprxz, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrky(nrank).lt. npry-1) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wcy1, mxz3, mpi_double_precision, nrank + nprxz, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wcy2, mxz3, mpi_double_precision, nrank + nprxz, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w3(:,:,iy_last-1,:)=wcy1(:,:,:)
      w3(:,:,iy_last,:)=wcy2(:,:,:)
	endif

    if (nrky(nrank) .eq. 0 ) then
      wcy1(:,:,:)=w3(:,:,iy_first+2,:)
      wcy2(:,:,:)=w3(:,:,iy_first+3,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wcy1, mxz3, mpi_double_precision, nrank +(npry-1)*nprxz, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wcy2, mxz3, mpi_double_precision, nrank +(npry-1)*nprxz, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrky(nrank).eq. npry-1) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wcy1, mxz3, mpi_double_precision, nrank -(npry-1)*nprxz, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wcy2, mxz3, mpi_double_precision, nrank -(npry-1)*nprxz, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w3(:,:,iy_last-1,:)=wcy1(:,:,:)
      w3(:,:,iy_last,:)=wcy2(:,:,:)
	endif
      else
      w3(:,:,iy_first+1,:)=w3(:,:,iy_last-2,:)
      w3(:,:,iy_first,:)=w3(:,:,iy_last-3,:)
      w3(:,:,iy_last-1,:)=w3(:,:,iy_first+2,:)
      w3(:,:,iy_last,:)=w3(:,:,iy_first+3,:)
      endif
      return
      end

!ws****************************************************************************************
      subroutine mpi_transfer1(w1)
      use declare
      real*8, dimension(mx,mz,my) :: w1
      real*8, dimension(mz,my) :: w1x1,w1x2
      real*8, dimension(mx,my) :: w1z1,w1z2
      real*8, dimension(mx,mz) :: w1y1,w1y2
      include 'mpif.h'

!       
! send w1 up unless i'm at the top, then receive from below
     
	if (nrankxz.ge.nrkz(nrank)*nprx .and. nrankxz.lt.nrkz(nrank)*nprx+nprx-1) then
      w1x1(:,:)=w1(ix_last-2,:,:)
      w1x2(:,:)=w1(ix_last-3,:,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( w1x1, myz, mpi_double_precision, nrank + 1, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( w1x2, myz, mpi_double_precision, nrank + 1, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrankxz.gt.nrkz(nrank)*nprx .and. nrankxz.le.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( w1x1, myz, mpi_double_precision, nrank - 1, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( w1x2, myz, mpi_double_precision, nrank - 1, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w1(ix_first+1,:,:)=w1x1(:,:)
      w1(ix_first,:,:)=w1x2(:,:)
	endif
	
      
! send w1 down unless i'm at the bottom

      
	if (nrankxz.gt.nrkz(nrank)*nprx .and. nrankxz.le.nrkz(nrank)*nprx+nprx-1) then
      w1x1(:,:)=w1(ix_first+2,:,:)
      w1x2(:,:)=w1(ix_first+3,:,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( w1x1, myz, mpi_double_precision, nrank - 1, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( w1x2, myz, mpi_double_precision, nrank - 1, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrankxz.ge.nrkz(nrank)*nprx .and. nrankxz.lt.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( w1x1, myz, mpi_double_precision, nrank + 1, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( w1x2, myz, mpi_double_precision, nrank + 1, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w1(ix_last-1,:,:)=w1x1(:,:)
      w1(ix_last,:,:)=w1x2(:,:)
	endif
!
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	if (nrankxz.lt.nsizexz-nprx) then
      w1z1(:,:)=w1(:,iz_last-2,:)
      w1z2(:,:)=w1(:,iz_last-3,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( w1z1, myx, mpi_double_precision, nrank + nprx, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( w1z2, myx, mpi_double_precision, nrank + nprx, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrankxz.ge.nprx ) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( w1z1, myx, mpi_double_precision, nrank - nprx, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( w1z2, myx, mpi_double_precision, nrank - nprx, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w1(:,iz_first+1,:)=w1z1(:,:)
      w1(:,iz_first,:)=w1z2(:,:)
	endif
	     
! send w1 down unless i'm at the bottom

      
	if (nrankxz.ge.nprx ) then
      w1z1(:,:)=w1(:,iz_first+2,:)
      w1z2(:,:)=w1(:,iz_first+3,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( w1z1, myx, mpi_double_precision, nrank - nprx, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( w1z2, myx, mpi_double_precision, nrank - nprx, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrankxz.lt.nsizews-nprx) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( w1z1, myx, mpi_double_precision, nrank + nprx, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( w1z2, myx, mpi_double_precision, nrank + nprx, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w1(:,iz_last-1,:)=w1z1(:,:)
      w1(:,iz_last,:)=w1z2(:,:)
	endif


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      if(npry .gt. 1) then
	if (nrky(nrank).lt. npry-1) then
      w1y1(:,:)=w1(:,:,iy_last-2)
      w1y2(:,:)=w1(:,:,iy_last-3)
!mpi   ----------------------------------------------------------------
	call mpi_send( w1y1, mxz, mpi_double_precision, nrank + nprxz, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( w1y2, mxz, mpi_double_precision, nrank + nprxz, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrky(nrank) .ge. 1 )  then
!mpi   ----------------------------------------------------------------
	call mpi_recv( w1y1, mxz, mpi_double_precision, nrank - nprxz, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( w1y2, mxz, mpi_double_precision, nrank - nprxz, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w1(:,:,iy_first+1)=w1y1(:,:)
      w1(:,:,iy_first)=w1y2(:,:)
	endif


  	if (nrky(nrank).eq. npry-1) then
      w1y1(:,:)=w1(:,:,iy_last-2)
      w1y2(:,:)=w1(:,:,iy_last-3)
!mpi   ----------------------------------------------------------------
	call mpi_send( w1y1, mxz, mpi_double_precision, nrank -(npry-1)*nprxz, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( w1y2, mxz, mpi_double_precision, nrank -(npry-1)*nprxz, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrky(nrank) .eq. 0 )  then
!mpi   ----------------------------------------------------------------
	call mpi_recv( w1y1, mxz, mpi_double_precision, nrank +(npry-1)*nprxz, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( w1y2, mxz, mpi_double_precision, nrank +(npry-1)*nprxz, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w1(:,:,iy_first+1)=w1y1(:,:)
      w1(:,:,iy_first)=w1y2(:,:)
	endif
	     
! send w1 down unless i'm at the bottom

      
	if (nrky(nrank) .ge. 1 ) then
      w1y1(:,:)=w1(:,:,iy_first+2)
      w1y2(:,:)=w1(:,:,iy_first+3)
!mpi   ----------------------------------------------------------------
	call mpi_send( w1y1, mxz, mpi_double_precision, nrank - nprxz, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( w1y2, mxz, mpi_double_precision, nrank - nprxz, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrky(nrank).lt. npry-1) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( w1y1, mxz, mpi_double_precision, nrank + nprxz, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( w1y2, mxz, mpi_double_precision, nrank + nprxz, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w1(:,:,iy_last-1)=w1y1(:,:)
      w1(:,:,iy_last)=w1y2(:,:)
	endif

    if (nrky(nrank) .eq. 0 ) then
      w1y1(:,:)=w1(:,:,iy_first+2)
      w1y2(:,:)=w1(:,:,iy_first+3)
!mpi   ----------------------------------------------------------------
	call mpi_send( w1y1, mxz, mpi_double_precision, nrank +(npry-1)*nprxz, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( w1y2, mxz, mpi_double_precision, nrank +(npry-1)*nprxz, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrky(nrank).eq. npry-1) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( w1y1, mxz, mpi_double_precision, nrank -(npry-1)*nprxz, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( w1y2, mxz, mpi_double_precision, nrank -(npry-1)*nprxz, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w1(:,:,iy_last-1)=w1y1(:,:)
      w1(:,:,iy_last)=w1y2(:,:)
	endif
      else
      w1(:,:,iy_first+1)=w1(:,:,iy_last-2)
      w1(:,:,iy_first)=w1(:,:,iy_last-3)
      w1(:,:,iy_last-1)=w1(:,:,iy_first+2)
      w1(:,:,iy_last)=w1(:,:,iy_first+3)
      endif
      return
      end

!ws****************************************************************************************
      subroutine mpi_transfersy1(wy1,wyt)
      use declare
      real*8, dimension(my) :: wy1
      real*8, dimension(myt) :: wyt
      include 'mpif.h'
      if(npry .gt. 1) then
      call mpi_allgather(wy1(3:my-2),mym,mpi_double_precision,wyt,mym,mpi_double_precision,mycomm_y,ierror)
      else
      wyt(1:mym)=wy1(3:my-2)
      endif
      return
      end

!ws****************************************************************************************
      subroutine mpi_transfersm_2d(ws,mm)
      use declare
      real*8, dimension(mx,mz,my,mm) :: ws
      real*8, dimension(mz,my,mm) :: wsx1,wsx2
      real*8, dimension(mx,my,mm) :: wsz1,wsz2
      include 'mpif.h'

!       
! send w8 up unless i'm at the top, then receive from below
     
	if (nrank.ge.nrkz(nrank)*nprx .and. nrank.lt.nrkz(nrank)*nprx+nprx-1) then
      wsx1(:,:,:)=ws(ix_last-2,:,:,:)
      wsx2(:,:,:)=ws(ix_last-3,:,:,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wsx1, myz*mm, mpi_double_precision, nrank + 1, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wsx2, myz*mm, mpi_double_precision, nrank + 1, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrank.gt.nrkz(nrank)*nprx .and. nrank.le.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wsx1, myz*mm, mpi_double_precision, nrank - 1, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wsx2, myz*mm, mpi_double_precision, nrank - 1, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      ws(ix_first+1,:,:,:)=wsx1(:,:,:)
      ws(ix_first,:,:,:)=wsx2(:,:,:)
	endif
	
      
! send w8 down unless i'm at the bottom

      
	if (nrank.gt.nrkz(nrank)*nprx .and. nrank.le.nrkz(nrank)*nprx+nprx-1) then
      wsx1(:,:,:)=ws(ix_first+2,:,:,:)
      wsx2(:,:,:)=ws(ix_first+3,:,:,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wsx1, myz*mm, mpi_double_precision, nrank - 1, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wsx2, myz*mm, mpi_double_precision, nrank - 1, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrank.ge.nrkz(nrank)*nprx .and. nrank.lt.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wsx1, myz*mm, mpi_double_precision, nrank + 1, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wsx2, myz*mm, mpi_double_precision, nrank + 1, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      ws(ix_last-1,:,:,:)=wsx1(:,:,:)
      ws(ix_last,:,:,:)=wsx2(:,:,:)
	endif
!
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	if (nrank.lt.nsize-nprx) then
      wsz1(:,:,:)=ws(:,iz_last-2,:,:)
      wsz2(:,:,:)=ws(:,iz_last-3,:,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wsz1, myx*mm, mpi_double_precision, nrank + nprx, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wsz2, myx*mm, mpi_double_precision, nrank + nprx, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrank.ge.nprx ) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wsz1, myx*mm, mpi_double_precision, nrank - nprx, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wsz2, myx*mm, mpi_double_precision, nrank - nprx, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      ws(:,iz_first+1,:,:)=wsz1(:,:,:)
      ws(:,iz_first,:,:)=wsz2(:,:,:)
	endif
	     
! send w8 down unless i'm at the bottom

      
	if (nrank.ge.nprx ) then
      wsz1(:,:,:)=ws(:,iz_first+2,:,:)
      wsz2(:,:,:)=ws(:,iz_first+3,:,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wsz1, myx*mm, mpi_double_precision, nrank - nprx, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wsz2, myx*mm, mpi_double_precision, nrank - nprx, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrank.lt.nsize-nprx) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wsz1, myx*mm, mpi_double_precision, nrank + nprx, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wsz2, myx*mm, mpi_double_precision, nrank + nprx, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      ws(:,iz_last-1,:,:)=wsz1(:,:,:)
      ws(:,iz_last,:,:)=wsz2(:,:,:)
	endif
      return
      end 

!ws****************************************************************************************
      subroutine mpi_transfer8_2d(w8)
      use declare
      real*8, dimension(mx,mz,my,8) :: w8
      real*8, dimension(mz,my,8) :: wfx1,wfx2
      real*8, dimension(mx,my,8) :: wfz1,wfz2
      include 'mpif.h'

!       
! send w8 up unless i'm at the top, then receive from below
     
	if (nrank.ge.nrkz(nrank)*nprx .and. nrank.lt.nrkz(nrank)*nprx+nprx-1) then
      wfx1(:,:,:)=w8(ix_last-2,:,:,:)
      wfx2(:,:,:)=w8(ix_last-3,:,:,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wfx1, myz8, mpi_double_precision, nrank + 1, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wfx2, myz8, mpi_double_precision, nrank + 1, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrank.gt.nrkz(nrank)*nprx .and. nrank.le.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wfx1, myz8, mpi_double_precision, nrank - 1, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wfx2, myz8, mpi_double_precision, nrank - 1, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w8(ix_first+1,:,:,:)=wfx1(:,:,:)
      w8(ix_first,:,:,:)=wfx2(:,:,:)
	endif
	
      
! send w8 down unless i'm at the bottom

      
	if (nrank.gt.nrkz(nrank)*nprx .and. nrank.le.nrkz(nrank)*nprx+nprx-1) then
      wfx1(:,:,:)=w8(ix_first+2,:,:,:)
      wfx2(:,:,:)=w8(ix_first+3,:,:,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wfx1, myz8, mpi_double_precision, nrank - 1, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wfx2, myz8, mpi_double_precision, nrank - 1, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrank.ge.nrkz(nrank)*nprx .and. nrank.lt.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wfx1, myz8, mpi_double_precision, nrank + 1, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wfx2, myz8, mpi_double_precision, nrank + 1, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w8(ix_last-1,:,:,:)=wfx1(:,:,:)
      w8(ix_last,:,:,:)=wfx2(:,:,:)
	endif
!
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	if (nrank.lt.nsize-nprx) then
      wfz1(:,:,:)=w8(:,iz_last-2,:,:)
      wfz2(:,:,:)=w8(:,iz_last-3,:,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wfz1, myx8, mpi_double_precision, nrank + nprx, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wfz2, myx8, mpi_double_precision, nrank + nprx, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrank.ge.nprx ) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wfz1, myx8, mpi_double_precision, nrank - nprx, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wfz2, myx8, mpi_double_precision, nrank - nprx, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w8(:,iz_first+1,:,:)=wfz1(:,:,:)
      w8(:,iz_first,:,:)=wfz2(:,:,:)
	endif
	     
! send w8 down unless i'm at the bottom

      
	if (nrank.ge.nprx ) then
      wfz1(:,:,:)=w8(:,iz_first+2,:,:)
      wfz2(:,:,:)=w8(:,iz_first+3,:,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wfz1, myx8, mpi_double_precision, nrank - nprx, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wfz2, myx8, mpi_double_precision, nrank - nprx, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrank.lt.nsize-nprx) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wfz1, myx8, mpi_double_precision, nrank + nprx, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wfz2, myx8, mpi_double_precision, nrank + nprx, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w8(:,iz_last-1,:,:)=wfz1(:,:,:)
      w8(:,iz_last,:,:)=wfz2(:,:,:)
	endif
      return
      end 

!ws****************************************************************************************
      subroutine mpi_transfer3_2d(w3)
      use declare
      real*8, dimension(mx,mz,my,3) :: w3
      real*8, dimension(mz,my,3) :: wcx1,wcx2
      real*8, dimension(mx,my,3) :: wcz1,wcz2
      include 'mpif.h'

       
! send w3 up unless i'm at the top, then receive from below
     
	if (nrank.ge.nrkz(nrank)*nprx .and. nrank.lt.nrkz(nrank)*nprx+nprx-1) then
      wcx1(:,:,:)=w3(ix_last-2,:,:,:)
      wcx2(:,:,:)=w3(ix_last-3,:,:,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wcx1, myz3, mpi_double_precision, nrank + 1, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wcx2, myz3, mpi_double_precision, nrank + 1, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrank.gt.nrkz(nrank)*nprx .and. nrank.le.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wcx1, myz3, mpi_double_precision, nrank - 1, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wcx2, myz3, mpi_double_precision, nrank - 1, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w3(ix_first+1,:,:,:)=wcx1(:,:,:)
      w3(ix_first,:,:,:)=wcx2(:,:,:)
	endif
	
      
! send w3 down unless i'm at the bottom

      
	if (nrank.gt.nrkz(nrank)*nprx .and. nrank.le.nrkz(nrank)*nprx+nprx-1) then
      wcx1(:,:,:)=w3(ix_first+2,:,:,:)
      wcx2(:,:,:)=w3(ix_first+3,:,:,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wcx1, myz3, mpi_double_precision, nrank - 1, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wcx2, myz3, mpi_double_precision, nrank - 1, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrank.ge.nrkz(nrank)*nprx .and. nrank.lt.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wcx1, myz3, mpi_double_precision, nrank + 1, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wcx2, myz3, mpi_double_precision, nrank + 1, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w3(ix_last-1,:,:,:)=wcx1(:,:,:)
      w3(ix_last,:,:,:)=wcx2(:,:,:)
	endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	if (nrank.lt.nsize-nprx) then
      wcz1(:,:,:)=w3(:,iz_last-2,:,:)
      wcz2(:,:,:)=w3(:,iz_last-3,:,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wcz1, myx3, mpi_double_precision, nrank + nprx, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wcz2, myx3, mpi_double_precision, nrank + nprx, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrank.ge.nprx ) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wcz1, myx3, mpi_double_precision, nrank - nprx, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wcz2, myx3, mpi_double_precision, nrank - nprx, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w3(:,iz_first+1,:,:)=wcz1(:,:,:)
      w3(:,iz_first,:,:)=wcz2(:,:,:)
	endif
	     
! send w3 down unless i'm at the bottom

      
	if (nrank.ge.nprx ) then
      wcz1(:,:,:)=w3(:,iz_first+2,:,:)
      wcz2(:,:,:)=w3(:,iz_first+3,:,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wcz1, myx3, mpi_double_precision, nrank - nprx, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wcz2, myx3, mpi_double_precision, nrank - nprx, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrank.lt.nsize-nprx) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wcz1, myx3, mpi_double_precision, nrank + nprx, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wcz2, myx3, mpi_double_precision, nrank + nprx, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w3(:,iz_last-1,:,:)=wcz1(:,:,:)
      w3(:,iz_last,:,:)=wcz2(:,:,:)
	endif
      return
      end

!ws****************************************************************************************
      subroutine mpi_transfer1_2d(w1)
      use declare
      real*8, dimension(mx,mz,my) :: w1
      real*8, dimension(mz,my) :: w1x1,w1x2
      real*8, dimension(mx,my) :: w1z1,w1z2
      include 'mpif.h'

!       
! send w1 up unless i'm at the top, then receive from below
     
	if (nrank.ge.nrkz(nrank)*nprx .and. nrank.lt.nrkz(nrank)*nprx+nprx-1) then
      w1x1(:,:)=w1(ix_last-2,:,:)
      w1x2(:,:)=w1(ix_last-3,:,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( w1x1, myz, mpi_double_precision, nrank + 1, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( w1x2, myz, mpi_double_precision, nrank + 1, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrank.gt.nrkz(nrank)*nprx .and. nrank.le.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( w1x1, myz, mpi_double_precision, nrank - 1, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( w1x2, myz, mpi_double_precision, nrank - 1, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w1(ix_first+1,:,:)=w1x1(:,:)
      w1(ix_first,:,:)=w1x2(:,:)
	endif
	
      
! send w1 down unless i'm at the bottom

      
	if (nrank.gt.nrkz(nrank)*nprx .and. nrank.le.nrkz(nrank)*nprx+nprx-1) then
      w1x1(:,:)=w1(ix_first+2,:,:)
      w1x2(:,:)=w1(ix_first+3,:,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( w1x1, myz, mpi_double_precision, nrank - 1, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( w1x2, myz, mpi_double_precision, nrank - 1, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrank.ge.nrkz(nrank)*nprx .and. nrank.lt.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( w1x1, myz, mpi_double_precision, nrank + 1, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( w1x2, myz, mpi_double_precision, nrank + 1, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w1(ix_last-1,:,:)=w1x1(:,:)
      w1(ix_last,:,:)=w1x2(:,:)
	endif
!
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	if (nrank.lt.nsize-nprx) then
      w1z1(:,:)=w1(:,iz_last-2,:)
      w1z2(:,:)=w1(:,iz_last-3,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( w1z1, myx, mpi_double_precision, nrank + nprx, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( w1z2, myx, mpi_double_precision, nrank + nprx, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrank.ge.nprx ) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( w1z1, myx, mpi_double_precision, nrank - nprx, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( w1z2, myx, mpi_double_precision, nrank - nprx, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w1(:,iz_first+1,:)=w1z1(:,:)
      w1(:,iz_first,:)=w1z2(:,:)
	endif
	     
! send w1 down unless i'm at the bottom

      
	if (nrank.ge.nprx ) then
      w1z1(:,:)=w1(:,iz_first+2,:)
      w1z2(:,:)=w1(:,iz_first+3,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( w1z1, myx, mpi_double_precision, nrank - nprx, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( w1z2, myx, mpi_double_precision, nrank - nprx, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrank.lt.nsize-nprx) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( w1z1, myx, mpi_double_precision, nrank + nprx, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( w1z2, myx, mpi_double_precision, nrank + nprx, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      w1(:,iz_last-1,:)=w1z1(:,:)
      w1(:,iz_last,:)=w1z2(:,:)
	endif
      return
      end

!hw*************************************************
!hw*************************************************
!only send the outermost layer, used for smooth case.(Haowei, 2018.01.12)
      subroutine mpi_transfersm_one_layer(ws,mm)
      use declare
      real*8, dimension(mx,mz,my,mm) :: ws
      real*8, dimension(mz,my,mm) :: wsx1,wsx2
      real*8, dimension(mx,my,mm) :: wsz1,wsz2
      real*8, dimension(mx,mz,mm) :: wsy1,wsy2
      include 'mpif.h'

!       
! send w8 up unless i'm at the top, then receive from below
     
	if (nrankxz.ge.nrkz(nrank)*nprx .and. nrankxz.lt.nrkz(nrank)*nprx+nprx-1) then
!      wsx1(:,:,:)=ws(ix_last-2,:,:,:)
      wsx2(:,:,:)=ws(ix_last-3,:,:,:)
!mpi   ----------------------------------------------------------------
!	call mpi_send( wsx1, myz*mm, mpi_double_precision, nrank + 1, 0,  &
!		      mpi_comm_world,ierror )
	call mpi_send( wsx2, myz*mm, mpi_double_precision, nrank + 1, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrankxz.gt.nrkz(nrank)*nprx .and. nrankxz.le.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
!	call mpi_recv( wsx1, myz*mm, mpi_double_precision, nrank - 1, 0,  &
!		      mpi_comm_world, status,ierror )
	call mpi_recv( wsx2, myz*mm, mpi_double_precision, nrank - 1, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
!      ws(ix_first+1,:,:,:)=wsx1(:,:,:)
      ws(ix_first,:,:,:)=wsx2(:,:,:)
	endif
	
      
! send w8 down unless i'm at the bottom

      
	if (nrankxz.gt.nrkz(nrank)*nprx .and. nrankxz.le.nrkz(nrank)*nprx+nprx-1) then
!      wsx1(:,:,:)=ws(ix_first+2,:,:,:)
      wsx2(:,:,:)=ws(ix_first+3,:,:,:)
!mpi   ----------------------------------------------------------------
!	call mpi_send( wsx1, myz*mm, mpi_double_precision, nrank - 1, 1,  &
!		      mpi_comm_world,ierror )
	call mpi_send( wsx2, myz*mm, mpi_double_precision, nrank - 1, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrankxz.ge.nrkz(nrank)*nprx .and. nrankxz.lt.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
!	call mpi_recv( wsx1, myz*mm, mpi_double_precision, nrank + 1, 1,  &
!		      mpi_comm_world, status,ierror )
	call mpi_recv( wsx2, myz*mm, mpi_double_precision, nrank + 1, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
!      ws(ix_last-1,:,:,:)=wsx1(:,:,:)
      ws(ix_last,:,:,:)=wsx2(:,:,:)
	endif
!
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	if (nrankxz.lt.nsizexz-nprx) then
!      wsz1(:,:,:)=ws(:,iz_last-2,:,:)
      wsz2(:,:,:)=ws(:,iz_last-3,:,:)
!mpi   ----------------------------------------------------------------
!	call mpi_send( wsz1, myx*mm, mpi_double_precision, nrank + nprx, 0,  &
!		      mpi_comm_world,ierror )
	call mpi_send( wsz2, myx*mm, mpi_double_precision, nrank + nprx, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrankxz.ge.nprx ) then
!mpi   ----------------------------------------------------------------
!	call mpi_recv( wsz1, myx*mm, mpi_double_precision, nrank - nprx, 0,  &
!		      mpi_comm_world, status,ierror )
	call mpi_recv( wsz2, myx*mm, mpi_double_precision, nrank - nprx, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
!      ws(:,iz_first+1,:,:)=wsz1(:,:,:)
      ws(:,iz_first,:,:)=wsz2(:,:,:)
	endif
	     
! send w8 down unless i'm at the bottom

      
	if (nrankxz.ge.nprx ) then
!      wsz1(:,:,:)=ws(:,iz_first+2,:,:)
      wsz2(:,:,:)=ws(:,iz_first+3,:,:)
!mpi   ----------------------------------------------------------------
!	call mpi_send( wsz1, myx*mm, mpi_double_precision, nrank - nprx, 1,  &
!		      mpi_comm_world,ierror )
	call mpi_send( wsz2, myx*mm, mpi_double_precision, nrank - nprx, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrankxz.lt.nsizexz-nprx) then
!mpi   ----------------------------------------------------------------
!	call mpi_recv( wsz1, myx*mm, mpi_double_precision, nrank + nprx, 1,  &
!		      mpi_comm_world, status,ierror )
	call mpi_recv( wsz2, myx*mm, mpi_double_precision, nrank + nprx, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
!      ws(:,iz_last-1,:,:)=wsz1(:,:,:)
      ws(:,iz_last,:,:)=wsz2(:,:,:)
	endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        if(npry .gt. 1) then
	if (nrky(nrank).lt. npry-1) then
!      wsy1(:,:,:)=ws(:,:,iy_last-2,:)
      wsy2(:,:,:)=ws(:,:,iy_last-3,:)
!mpi   ----------------------------------------------------------------
!	call mpi_send( wsy1, mxz*mm, mpi_double_precision, nrank + nprxz, 0,  &
!		      mpi_comm_world,ierror )
	call mpi_send( wsy2, mxz*mm, mpi_double_precision, nrank + nprxz, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrky(nrank) .ge. 1 )  then
!mpi   ----------------------------------------------------------------
!	call mpi_recv( wsy1, mxz*mm, mpi_double_precision, nrank - nprxz, 0,  &
!		      mpi_comm_world, status,ierror )
	call mpi_recv( wsy2, mxz*mm, mpi_double_precision, nrank - nprxz, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
!      ws(:,:,iy_first+1,:)=wsy1(:,:,:)
      ws(:,:,iy_first,:)=wsy2(:,:,:)
	endif


  	if (nrky(nrank).eq. npry-1) then
!      wsy1(:,:,:)=ws(:,:,iy_last-2,:)
      wsy2(:,:,:)=ws(:,:,iy_last-3,:)
!mpi   ----------------------------------------------------------------
!	call mpi_send( wsy1, mxz*mm, mpi_double_precision, nrank -(npry-1)*nprxz, 0,  &
!		      mpi_comm_world,ierror )
	call mpi_send( wsy2, mxz*mm, mpi_double_precision, nrank -(npry-1)*nprxz, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrky(nrank) .eq. 0 )  then
!mpi   ----------------------------------------------------------------
!	call mpi_recv( wsy1, mxz*mm, mpi_double_precision, nrank +(npry-1)*nprxz, 0,  &
!		      mpi_comm_world, status,ierror )
	call mpi_recv( wsy2, mxz*mm, mpi_double_precision, nrank +(npry-1)*nprxz, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
!      ws(:,:,iy_first+1,:)=wsy1(:,:,:)
      ws(:,:,iy_first,:)=wsy2(:,:,:)
	endif
	     
! send w8 down unless i'm at the bottom

      
	if (nrky(nrank) .ge. 1 ) then
!      wsy1(:,:,:)=ws(:,:,iy_first+2,:)
      wsy2(:,:,:)=ws(:,:,iy_first+3,:)
!mpi   ----------------------------------------------------------------
!	call mpi_send( wsy1, mxz*mm, mpi_double_precision, nrank - nprxz, 1,  &
!		      mpi_comm_world,ierror )
	call mpi_send( wsy2, mxz*mm, mpi_double_precision, nrank - nprxz, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrky(nrank).lt. npry-1) then
!mpi   ----------------------------------------------------------------
!	call mpi_recv( wsy1, mxz*mm, mpi_double_precision, nrank + nprxz, 1,  &
!		      mpi_comm_world, status,ierror )
	call mpi_recv( wsy2, mxz*mm, mpi_double_precision, nrank + nprxz, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
!      ws(:,:,iy_last-1,:)=wsy1(:,:,:)
      ws(:,:,iy_last,:)=wsy2(:,:,:)
	endif

    if (nrky(nrank) .eq. 0 ) then
!      wsy1(:,:,:)=ws(:,:,iy_first+2,:)
      wsy2(:,:,:)=ws(:,:,iy_first+3,:)
!mpi   ----------------------------------------------------------------
!	call mpi_send( wsy1, mxz*mm, mpi_double_precision, nrank +(npry-1)*nprxz, 1,  &
!		      mpi_comm_world,ierror )
	call mpi_send( wsy2, mxz*mm, mpi_double_precision, nrank +(npry-1)*nprxz, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrky(nrank).eq. npry-1) then
!mpi   ----------------------------------------------------------------
!	call mpi_recv( wsy1, mxz*mm, mpi_double_precision, nrank -(npry-1)*nprxz, 1,  &
!		      mpi_comm_world, status,ierror )
	call mpi_recv( wsy2, mxz*mm, mpi_double_precision, nrank -(npry-1)*nprxz, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
!      ws(:,:,iy_last-1,:)=wsy1(:,:,:)
      ws(:,:,iy_last,:)=wsy2(:,:,:)
	endif
      else
!      ws(:,:,iy_first+1,:)=ws(:,:,iy_last-2,:)
      ws(:,:,iy_first,:)=ws(:,:,iy_last-3,:)
!      ws(:,:,iy_last-1,:)=ws(:,:,iy_first+2,:)
      ws(:,:,iy_last,:)=ws(:,:,iy_first+3,:)
      endif

      return
      end 

