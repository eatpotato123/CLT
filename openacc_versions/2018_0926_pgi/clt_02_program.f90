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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the main program for the CLT code, time iterations
! call input, initia, setdt, stepon, right, recrd, diagn, etc.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

include 'clt_01_modules.f90'

program tokrz
    use declare
    use declare_for_openacc
    use declare_oxpoint
!	use omp_lib
    USE omp_lib
    include 'mpif.h'
!	include 'fftw3.f'
!#include <finclude/petscsys.h>
             
!      petscerrorcode ierr

!      integer status(mpi_status_size)
! mpi   -----------------------------------------------------------------
! initiate mpi
    call mpi_init(ierror)    
! clt-k using petsc
!   call petscinitialize(petsc_null_character,ierr)

    call mpi_comm_size(mpi_comm_world,nsize,ierror)
    call mpi_comm_rank(mpi_comm_world,nrank,ierror)
    !$acc set device_num(0)

    if (nsize.ne.npr) then
    print*,'the number of processors is not equal to', npr
    call mpi_abort(mpi_comm_world,1)
    stop
    endif
!mpi   ------------------------------------------------------------------
!    if(nrank==0) timestart=mpi_wtime()
! initializing
    nstep=0
    nstop=1000
    nrcd=1000
!    nrcd=100
    ndgn=10000
    neng=10000
    ncd=10
    nper=100
    nrank_dgn=nsize-1
    mtor_dgn=4
    
    call input
    nsizexz=nsize/npry

    do nrk=0,nsize-1
    nrky(nrk)=nrk/nprxz
    nrkxz(nrk)=nrk-nrky(nrk)*nprxz

    nrkz(nrk)=nrkxz(nrk)/nprx
    nrkx(nrk)=nrkxz(nrk)-nrkz(nrk)*nprx
    enddo
    nranky=nrky(nrank)
    nrankxz=nrkxz(nrank)
    nrankx=nrkx(nrank)
    nrankz=nrkz(nrank)
  
    do nrk=0,nsizexz-1
    ix_min(nrk)=1
    ix_max(nrk)=mx
    iz_min(nrk)=1
    iz_max(nrk)=mz
    if(nrkx(nrk)==0)  ix_min(nrk)=ix_min(nrk)+2
    if(nrkx(nrk)==nprx-1)  ix_max(nrk)=ix_max(nrk)-2
    if(nrkz(nrk)==0)  iz_min(nrk)=iz_min(nrk)+2
    if(nrkz(nrk)==nprz-1)  iz_max(nrk)=iz_max(nrk)-2
    enddo



!mpi_comm
!    call mpi_comm_group(mpi_comm_world,group_world,ierror)
!    do kxz=0,nsizexz-1
!    do ky=0,npry-1
!    nrank_in_gy(ky,kxz)=kxz+ky*nprxz
!    enddo
!
!    call mpi_group_incl(group_world,npry,nrank_in_gy(:,kxz),group_y(kxz),ierror)
!    call mpi_comm_create(mpi_comm_world,group_y(kxz),comm_y(kxz),ierror)
!    call mpi_comm_size(comm_y(kxz), nysize(kxz), ierror)
!    call mpi_comm_rank(comm_y(kxz), nyrank(kxz), ierror)
!
!    call mpi_group_rank(group_y(kxz),ngy_rank(kxz),ierror) 
!    enddo
    call mpi_comm_split(mpi_comm_world,nrankxz,0,mycomm_y,ierror)
    call mpi_comm_rank(mycomm_y,nrank_gy,ierror) 
    call mpi_comm_size(mycomm_y,nsize_gy,ierror) 
    call mpi_allgather(nrank,1,mpi_integer,nrankiny,1,mpi_integer,mycomm_y,ierror)
    write(*,*) nrank,'commy=',mycomm_y,nrankxz,nrank_gy,nrankiny(:)
!mpi_comm


    ix_first=1
    ix_last=mx
    iz_first=1
    iz_last=mz
    iy_first=1
    iy_last=my
    if(nrkx(nrank)==0)  then
    ix_first=ix_first+2
    xx(2)=xxt(1)-dxt(1)
    xx(1)=xx(2)-dxt(1)
    endif
    if(nrkx(nrank)==nprx-1) then
    ix_last=ix_last-2
    xx(mx-1)=xxt(mxt)+dxt(mxt)
    xx(mx)=xx(mx-1)+dxt(mxt)
    endif
    if(nrkz(nrank)==0) then
    iz_first=iz_first+2
    zz(2)=zzt(1)-dzt(1)
    zz(1)=zz(2)-dzt(1)
    endif
    if(nrkz(nrank)==nprz-1) then
    iz_last=iz_last-2
    zz(mz-1)=zzt(mzt)+dzt(mzt)
    zz(mz)=zz(mz-1)+dzt(mzt)
    endif

    do jx=ix_first,ix_last
    xx(jx)=xxt(nrkx(nrank)*mxm+jx-2)
    enddo
    do jz=iz_first,iz_last
    zz(jz)=zzt(nrkz(nrank)*mzm+jz-2)
    enddo
    do jy=iy_first,iy_last
    yy(jy)=yyt(nrky(nrank)*mym+jy-2)
    enddo
!
    call initia
    call init_dgn
 
    call init_ps1

! haowei add for cut cell initialize
    call find_bnd_grd
!	call find_bnd_grd_in_each_proc
!    call map_nova_to_bnd_grd ! select case is not allowed in SUNWAY
    call map_int_to_bnd_grd
!   call map_nova_to_bnd_grd_in_each_proc
	call find_ir_pnt_bndx
	call find_ir_pnt_bndz
	call decide_grd_type_in_each_proc
	call calculate_dnfm_coeff_for_ir_bndx_point
!	call decide_hyperbolic_ratio
	call decide_hyperbolic_ratio_v2
	call data_type_weight
	call jx_jz_list_bndry8_cut_cell_v2_fixed
  211 format(2(1x,e12.5)) 

	if(nrank.eq.0) then
	open(unit=193,file='grd_ir_bndx.dat',status='unknown',form='formatted')
	open(unit=194,file='grd_ir_bndz.dat',status='unknown',form='formatted')
!	write(193,511)(((xxt(jx),zzt(jz),grd_type(jx,jz,1)*1.d0,gdtp_bndx(jx,jz,1)*1.d0,gdtp_bndx(jx,jz,2)*1.d0),jx=1,mxt),jz=1,mzt)
!	write(194,511)(((xxt(jx),zzt(jz),grd_type(jx,jz,2)*1.d0,gdtp_bndz(jx,jz,1)*1.d0,gdtp_bndz(jx,jz,2)*1.d0),jx=1,mxt),jz=1,mzt)
	close(193)
	close(194)
  511 format(5(1x,e12.5)) 
	endif

!hw*************************************************************
!hw*************************************************************
! to find out the region for IR points
! integer ix_first_irpt,ix_last_irpt,iz_first_irpt,iz_last_irpt ! the rank region for  the IR points, gdtp_ep(jx,jz,1or4)==2
	ix_first_irpt=ix_last
	ix_last_irpt=ix_first
	iz_first_irpt=iz_last
	iz_last_irpt=iz_first
	
	do jz=iz_first,iz_last
	do jx=ix_first,ix_last

	if(((gdtp_ep(jx,jz,1).ne.1).and.(gdtp_ep(jx,jz,1).ne.4)).or.((gdtp_ep(jx,jz,4).ne.1).and.(gdtp_ep(jx,jz,4).ne.4))) then
		  ix_first_irpt=min(ix_first_irpt,jx)
		  ix_last_irpt=max(ix_last_irpt,jx)
		  iz_first_irpt=min(iz_first_irpt,jz)
		  iz_last_irpt=max(iz_last_irpt,jz)
	endif

	enddo
	enddo

	if(ix_first_irpt.gt.ix_last_irpt) ix_last_irpt=ix_first_irpt
	if(iz_first_irpt.gt.iz_last_irpt) iz_last_irpt=iz_first_irpt

	ix_first_irpt2=max(ix_first_irpt,ix_first+2)
	ix_last_irpt2=min(ix_last_irpt,ix_last-2)
	iz_first_irpt2=max(iz_first_irpt,iz_first+2)
	iz_last_irpt2=min(iz_last_irpt,iz_last-2)
	
!	print*,'ixiz_irpt=',ix_first_irpt,ix_last_irpt,iz_first_irpt,iz_last_irpt,nrank
!hw*************************************************************
!hw*************************************************************



!	if(rmp_east) call rmp_east_pert
	if(rmp_east) call rmp_east_pert_withA
	

  !   if(lpetsc .or. lpll .ge. 4) call init_a  
!    x1=x
!    call calculate_ps1
!    call recrd_ps1
!    goto 900
	if(initia_from_pgfile) then
	pssmw=5.*(psia-psmin)/200.d0
	endif
      pstransw=pssmw
      pstrans=psia-pssmw !*cos(time)
      
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
!	if(psi(jx,jz).lt.psia1) then ! keep the diffusion term cofficient only in closed flux, haowei, 20180402
      if(gdtp_ep(jx,jz,1).ne.4) then
        cfsmb(jx,jz)=cfsm0+0.5*cfsmout*(1+tanh((psi(jx,jz)-pstrans)/pstransw))
        etaint(jx,jz)=eta0 !*(1+(etacut-1)*tanh((tmint(jx,jz)/tm00)**(-1.5)/(etacut-1)))
        eta(jx,jz,:)=etaint(jx,jz)

!	  if(psi(jx,jz).gt.psia1) eta(jx,jz,:)=1.d-3

        fmu(jx,jz)=fmu0+0.5*fmuout*(1+tanh((psi(jx,jz)-pstrans)/pstransw))+fmuc*(1-tanh((psi(jx,jz)-psmin)/pstransw**0.5))
        pmu(jx,jz)=pmu0+0.5*pmuout*(1+tanh((psi(jx,jz)-pstrans)/pstransw))+pmuc*(1-tanh((psi(jx,jz)-psmin)/pstransw**0.5))
        kap(jx,jz)=kap0+0.5*kapout*(1+tanh((psi(jx,jz)-pstrans)/pstransw))+kapc*(1-tanh((psi(jx,jz)-psmin)/pstransw**0.5))
        kap_ll(jx,jz)=kapll0 !*(tanh((psi(jx,jz)-psmin)/pstransw**0.5))
        kap_pp(jx,jz)=kap(jx,jz)
        
        fdi(jx,jz)=0.5*fdiout*(1-tanh((psi(jx,jz)-pstrans+1.*pstransw)/pstransw))

!        etb(jx,jz)=0.5*etbout*(1+tanh((psi(jx,jz)-pstrans)/pstransw))+etbc*(1-tanh((psi(jx,jz)-psmin)/pstransw**0.5))
        etb(jx,jz)=0.5*etbout*(1-tanh((xx(jx)-xmin-0.1*(xzero-xmin))/(0.03*(xzero-xmin))))*exp(-zz(jz)**2/(0.2*(zmax-zmin))**2) !+etbc*(1-tanh((psi(jx,jz))/pstransw**0.5))
        cdb(jx,jz)=cdb0+0.5*cdbout*(1+tanh((psi(jx,jz)-pstrans)/pstransw))
        etax(jx,jz,:)=0.
        etaz(jx,jz,:)=0.
        etay(jx,jz,:)=0.
        cdbx(jx,jz)=0.5*cdbout/cosh((psi(jx,jz)-pstrans)/pstransw)**2*psi_dx(jx,jz)/pstransw
        cdbz(jx,jz)=0.5*cdbout/cosh((psi(jx,jz)-pstrans)/pstransw)**2*psi_dz(jx,jz)/pstransw
        fmux(jx,jz)=0.5*fmuout/cosh((psi(jx,jz)-pstrans)/pstransw)**2*psi_dx(jx,jz) &
                /pstransw-fmuc/cosh((psi(jx,jz)-psmin)/pstransw**0.5)**2*psi_dx(jx,jz)/pstransw**0.5
        fmuz(jx,jz)=0.5*fmuout/cosh((psi(jx,jz)-pstrans)/pstransw)**2*psi_dz(jx,jz) &
                /pstransw-fmuc/cosh((psi(jx,jz)-psmin)/pstransw**0.5)**2*psi_dz(jx,jz)/pstransw**0.5

        pmux(jx,jz)=0.5*pmuout/cosh((psi(jx,jz)-pstrans)/pstransw)**2*psi_dx(jx,jz) &
                /pstransw-pmuc/cosh((psi(jx,jz)-psmin)/pstransw**0.5)**2*psi_dx(jx,jz)/pstransw**0.5
        pmuz(jx,jz)=0.5*pmuout/cosh((psi(jx,jz)-pstrans)/pstransw)**2*psi_dz(jx,jz) &
                /pstransw-pmuc/cosh((psi(jx,jz)-psmin)/pstransw**0.5)**2*psi_dz(jx,jz)/pstransw**0.5

        kapx(jx,jz)=0.5*kapout/cosh((psi(jx,jz)-pstrans)/pstransw)**2*psi_dx(jx,jz) &
                /pstransw-kapc/cosh((psi(jx,jz)-psmin)/pstransw**0.5)**2*psi_dx(jx,jz)/pstransw**0.5
        kapz(jx,jz)=0.5*kapout/cosh((psi(jx,jz)-pstrans)/pstransw)**2*psi_dz(jx,jz) &
                /pstransw-kapc/cosh((psi(jx,jz)-psmin)/pstransw**0.5)**2*psi_dz(jx,jz)/pstransw**0.5
        kap_ll_dx(jx,jz)=0
        kap_ll_dz(jx,jz)=0
        kap_pp_dx(jx,jz)=kapx(jx,jz)
        kap_pp_dz(jx,jz)=kapz(jx,jz)
!        etbx(jx,jz)=0.5*etbout/cosh((psi(jx,jz)-pstrans)/pstransw)**2*psi_dx(jx,jz)/pstransw-etbc/cosh((psi(jx,jz)-psmin)/pstransw**0.5)**2*psi_dx(jx,jz)/pstransw**0.5
!        etbz(jx,jz)=0.5*etbout/cosh((psi(jx,jz)-pstrans)/pstransw)**2*psi_dz(jx,jz)/pstransw-etbc/cosh((psi(jx,jz)-psmin)/pstransw**0.5)**2*psi_dz(jx,jz)/pstransw**0.5
        etbx(jx,jz)=-0.5*etbout/cosh((xx(jx)-xmin-0.1*(xzero-xmin))/(0.02*(xzero-xmin)))**2/(0.02*(xzero-xmin))*exp(-zz(jz)**2/(0.2*(zmax-zmin))**2)
        etbz(jx,jz)=0.5*etbout*(1-tanh((xx(jx)-xmin-0.1*(xzero-xmin))/(0.02*(xzero-xmin)))) &
                *exp(-zz(jz)**2/(0.2*(zmax-zmin))**2)*(-2*zz(jz)/(0.2*(zmax-zmin))**2)

	else ! all set to zero
        cfsmb(jx,jz)=0.
        etaint(jx,jz)=0.
        eta(jx,jz,:)=0.
        fmu(jx,jz)=0.
        pmu(jx,jz)=0.
        kap(jx,jz)=0.
        kap_ll(jx,jz)=0.
        kap_pp(jx,jz)=0.
        fdi(jx,jz)=0.
        etb(jx,jz)=0.
        cdb(jx,jz)=0.
        etax(jx,jz,:)=0.
        etaz(jx,jz,:)=0.
        etay(jx,jz,:)=0.
        cdbx(jx,jz)=0.
        cdbz(jx,jz)=0.
        fmux(jx,jz)=0.
        fmuz(jx,jz)=0.

        pmux(jx,jz)=0.
        pmuz(jx,jz)=0.

        kapx(jx,jz)=0.
        kapz(jx,jz)=0.
        kap_ll_dx(jx,jz)=0
        kap_ll_dz(jx,jz)=0
        kap_pp_dx(jx,jz)=0.
        kap_pp_dz(jx,jz)=0.
!        etbx(jx,jz)=0.5*etbout/cosh((psi(jx,jz)-pstrans)/pstransw)**2*psi_dx(jx,jz)/pstransw-etbc/cosh((psi(jx,jz)-psmin)/pstransw**0.5)**2*psi_dx(jx,jz)/pstransw**0.5
!        etbz(jx,jz)=0.5*etbout/cosh((psi(jx,jz)-pstrans)/pstransw)**2*psi_dz(jx,jz)/pstransw-etbc/cosh((psi(jx,jz)-psmin)/pstransw**0.5)**2*psi_dz(jx,jz)/pstransw**0.5
        etbx(jx,jz)=0.
        etbz(jx,jz)=0.
	endif

      enddo
      enddo

      if(lrstrt) then
      call readin
      do jy=iy_first,iy_last
      x1(:,:,jy,:)=x(:,:,jy,:)-xint(:,:,:)
      enddo
      end if
      timein=time
      nstin=nst

!      call calculate_ps1
!      call recrd_ps1
!
    if(nrank==0) then
    if(nstep==0) then
!    open(unit=11,file='dgnmax1.dat',status='unknown',form='formatted')
    open(unit=111,file='dgnmax.dat',status='unknown',form='formatted')
    open(unit=112,file='dgnmin.dat',status='unknown',form='formatted')
!    open(unit=13,file='current.dat',status='unknown',form='formatted')
!    open(unit=15,file='bfield.dat',status='unknown',form='formatted')
    open(unit=12,file='energy.dat',status='unknown',form='formatted') 
    open(unit=14,file='divb.dat',status='unknown',form='formatted')
    open(unit=18,file='nstime.dat',status='unknown',form='formatted') 
    open(unit=211,file='eymaxmin.dat',status='unknown',form='formatted')    
    open(unit=311,file='oxpoint.dat',status='unknown',form='formatted')
    open(unit=413,file='debug.dat',status='unknown',form='formatted')
    else
    open(unit=11,file='dgnmax1.dat',status='unknown',form='formatted',position='append')
    open(unit=111,file='dgnmax.dat',status='unknown',form='formatted',position='append')
    open(unit=112,file='dgnmin.dat',status='unknown',form='formatted',position='append')
!    open(unit=13,file='current.dat',status='unknown',form='formatted',position='append')
!   open(unit=15,file='bfield.dat',status='unknown',form='formatted',position='append')
    open(unit=12,file='energy.dat',status='unknown',form='formatted',position='append') 
    open(unit=14,file='divb.dat',status='unknown',form='formatted',position='append') 
    open(unit=18,file='nstime.dat',status='unknown',form='formatted',position='append')
    open(unit=211,file='eymaxmin.dat',status='unknown',form='formatted',position='append')
    open(unit=311,file='oxpoint.dat',status='unknown',form='formatted',position='append')
    open(unit=413,file='debug.dat',status='unknown',form='formatted',position='append')
    endif
    endif

    if(nrank==nrank_mode) then
    if(nstep==0) then
    open(unit=16,file='xmode.dat',status='unknown',form='formatted') 
    else
    open(unit=16,file='xmode.dat',status='unknown',form='formatted',position='append') 
    endif
    endif
    
  !curdriven
    if(lcd .gt. 0) then
    tcds=time+100
    tcdd=500
    tcde=tcds+10000
    delcd=0.03*(psmax-psmin)
    psshift=0 !0.01*(psmax-psmin)
!    ps100=0.001*(psmax-psmin)
    ps1=0
    br_lim=4.02e-4
    lcdox=2
    if(.not. lrstrt_cd) then
      call find_oxpoint_1st
      call diagn_brmax0
    endif
    call distribution_cd
    if(lrstrt_cd) then 
      call readin_cud
!      call mpi_barrier(mpi_comm_world,ierror)
      call find_cud_ox
      call distribution_cd_oxpoint(cud(:,:,3,2))
    endif
   !_cos     
    endif

!curdriven
    if(eta_from_t) call calculate_eta
    call recrd_dssp  
    
    
!    	call find_bnd_grd
!!	call find_bnd_grd_in_each_proc
!     	call map_nova_to_bnd_grd
!!    	call map_nova_to_bnd_grd_in_each_proc
!	call find_ir_pnt_bndx
!	call find_ir_pnt_bndz
!	call decide_grd_type_in_each_proc
!	call calculate_dnfm_coeff_for_ir_bndx_point
!!	call decide_hyperbolic_ratio
!	call decide_hyperbolic_ratio_v2
!	call data_type_weight
!	if(nrank.eq.0) then
!	open(unit=193,file='grd_ir_bndx.dat',status='unknown',form='formatted')
!	open(unit=194,file='grd_ir_bndz.dat',status='unknown',form='formatted')
!	write(193,511)(((xxt(jx),zzt(jz),grd_type(jx,jz,1)*1.d0,gdtp_bndx(jx,jz,1)*1.d0,gdtp_bndx(jx,jz,2)*1.d0),jx=1,mxt),jz=1,mzt)
!	write(194,511)(((xxt(jx),zzt(jz),grd_type(jx,jz,2)*1.d0,gdtp_bndz(jx,jz,1)*1.d0,gdtp_bndz(jx,jz,2)*1.d0),jx=1,mxt),jz=1,mzt)
!	close(193)
!	close(194)
!  511 format(5(1x,e12.5)) 
!	endif


	if(rmp_east) then
	allocate(b_rmp_bndx_tmp1(nbndx,my,3))
	allocate(b_rmp_bndz_tmp1(nbndz,my,3))
	allocate(b_rmp_bndx_tmp2(nbndx,my,3))
	allocate(b_rmp_bndz_tmp2(nbndz,my,3))

! nrmp=1, shifted pi/2 between b_tmp1 and b_tmp2		  	
	b_rmp_tmp1=b_rmp_out
	b_rmp_tmp2(:,:,1:myt*3/4+1,:)=b_rmp_out(:,:,myt/4+1:myt+1,:)
	b_rmp_tmp2(:,:,myt*3/4+2:myt,:)=b_rmp_out(:,:,2:myt/4,:)
	b_rmp_tmp2(:,:,myt+1:myt+4,:)=b_rmp_tmp2(:,:,1:4,:)
		  	
	b_rmp_bndx_tmp1=b_rmp_bndx
	b_rmp_bndx_tmp2(:,1:myt*3/4+1,:)=b_rmp_bndx(:,myt/4+1:myt+1,:)
	b_rmp_bndx_tmp2(:,myt*3/4+2:myt,:)=b_rmp_bndx(:,2:myt/4,:)
	b_rmp_bndx_tmp2(:,myt+1:myt+4,:)=b_rmp_bndx_tmp2(:,1:4,:)

	b_rmp_bndz_tmp1=b_rmp_bndz
	b_rmp_bndz_tmp2(:,1:myt*3/4+1,:)=b_rmp_bndz(:,myt/4+1:myt+1,:)
	b_rmp_bndz_tmp2(:,myt*3/4+2:myt,:)=b_rmp_bndz(:,2:myt/4,:)
	b_rmp_bndz_tmp2(:,myt+1:myt+4,:)=b_rmp_bndz_tmp2(:,1:4,:)

	omega_rmp=0.d0
	start_rmp=0.d0
	finish_rmp=800.d0
!	time_old=start_rmp
	tau_rmp=10.d0
	endif

      if(nrank==0) timestart=mpi_wtime()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! start the loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!$acc update device(gdtp_ep)
	!$acc update device(hypb_ratio)
	!$acc update device(bndx_grd,bndz_grd)
	!$acc update device(x_8bndz,x_8bndx,x1_8bndz,x1_8bndx,ef_3bndx,ef_3bndz,xint8_bndz,xint_8bndx,xint_8bndz)
	!$acc update device(b_rmp_bndx,b_rmp_bndz)
	!$acc update device(xy_8bndx,xy2_8bndx)
	!$acc update device(xy_8bndz,xy2_8bndz)
      !$acc update device(az1_irx,bz1_irx,cz1_irx,dz1_irx)
      !$acc update device(ax1_irz,bx1_irz,cx1_irz,dx1_irz)
      !$acc update device(azbp_irx,bzbp_irx,czbp_irx,dzbp_irx)
      !$acc update device(azbm_irx,bzbm_irx,czbm_irx,dzbm_irx)
     	!$acc update device(axbp_irz,bxbp_irz,cxbp_irz,dxbp_irz)
      !$acc update device(axbm_irz,bxbm_irz,cxbm_irz,dxbm_irz)
      !$acc update device(a2zbm_irx,b2zbm_irx,c2zbm_irx)
      !$acc update device(a2zbp_irx,b2zbp_irx,c2zbp_irx)
      !$acc update device(a2xbm_irz,b2xbm_irz,c2xbm_irz)
      !$acc update device(a2xbp_irz,b2xbp_irz,c2xbp_irz)
      !$acc update device(az2_irx,bz2_irx,cz2_irx,dz2_irx)
      !$acc update device(ax2_irz,bx2_irz,cx2_irz,dx2_irz)
	!$acc update device(yy)
	!$acc update device(xx)
	!$acc update device(zz)
	!$acc update device(rvx,rey, drvx_dx_tmp, dex_dy_tmp, dez_dy_tmp)
	!$acc update device(rvx_bndx,rey_bndx)
	!$acc update device(rvx_bndz,rey_bndz)
	!$acc update device(rby)
	!$acc update device(wwy)
	!$acc update device(x1r_tmp,xr2_tmp,x1z_tmp,xz2_tmp,x1_tmp)
	!$acc update device(dt)
	!$acc update device(ft_rmp)
	!$acc update device(eta)
	!$acc update device(fmu,fmux,fmuz)
	!$acc update device(pmu,pmux,pmuz)
	!$acc update device(kap,kapx,kapz)
      !$acc update device(ax1,bx1,cx1,dx1)
      !$acc update device(ax2,bx2,cx2,dx2)
      !$acc update device(az1,bz1,cz1,dz1)
      !$acc update device(az2,bz2,cz2,dz2)
	!$acc update device(ay1,by1,cy1,dy1)
	!$acc update device(ay2,by2,cy2,dy2)
	!$acc update device(x,xm,xfold,xdif,x1)
	!$acc update device(xr,xy,xz,xr2,xy2,xz2,x1r,x1z)
	!$acc update device(cur,ef)
	!$acc update device(eta1)
	!$acc update device(xint,xint_dx,xint_dz)
	!$acc update device(xint3,xint4,xint5,xint6,xint7,xint8)
	!$acc update device(xint1_dx,xint2_dx,xint3_dx,xint4_dx,xint5_dx)
	!$acc update device(xint1_dz,xint2_dz,xint3_dz,xint4_dz,xint5_dz)
	!$acc update device(cint)
	!$acc update device(cint1,cint2,cint3)
	!$acc update device(b_rmp_out)
	!$acc update device(fdi)
	!$acc update device(jxjz_list1_bndry)
	!$acc update device(jxjz_list2_bndry)
	!$acc update device(jxjz_list3_bndry)
	!$acc update device(jxjz_list4_bndry)
	!$acc enter data copyin(wfx1,wfx2,wfz1,wfz2,wfy1,wfy2,wcx1,wcx2,wcz1,wcz2,wcy1,wcy2)
100 continue
	
! updated the rotation model 0601, haowei
	if(rmp_east) then ! open the rmp after 160tA

!	ft_rmp(1)=1.d0-exp(-max(0.d0,time-dt-start_rmp)/tau_rmp) ! the old time
!	ft_rmp(2)=1.d0-exp(-max(0.d0,time-start_rmp)/tau_rmp) ! the new time

	if(time.lt.finish_rmp) then
	ft_rmp(1)=1.d0*max(0.d0,min(time-dt-start_rmp,tau_rmp))/tau_rmp ! the old time
	ft_rmp(2)=1.d0*max(0.d0,min(time-start_rmp,tau_rmp))/tau_rmp ! the new time
	else
	ft_rmp(1)=1.d0-1.d0*max(0.d0,min(time-dt-finish_rmp,tau_rmp))/tau_rmp
	ft_rmp(2)=1.d0-1.d0*max(0.d0,min(time-finish_rmp,tau_rmp))/tau_rmp
	endif

	x(:,:,:,6:8)=x(:,:,:,6:8)-b_rmp_out(:,:,:,1:3)*ft_rmp(1)
   	x_8bndx(:,:,6:8)=x_8bndx(:,:,6:8)-b_rmp_bndx(:,:,1:3)*ft_rmp(1)
   	x_8bndz(:,:,6:8)=x_8bndz(:,:,6:8)-b_rmp_bndz(:,:,1:3)*ft_rmp(1)
	
	if(omega_rmp.ne.0.d0) then
	b_rmp_out=b_rmp_tmp1*cos(omega_rmp*(time-start_rmp))+b_rmp_tmp2*sin(omega_rmp*(time-start_rmp))
	b_rmp_bndx=b_rmp_bndx_tmp1*cos(omega_rmp*(time-start_rmp))+b_rmp_bndx_tmp2*sin(omega_rmp*(time-start_rmp))
	b_rmp_bndz=b_rmp_bndz_tmp1*cos(omega_rmp*(time-start_rmp))+b_rmp_bndz_tmp2*sin(omega_rmp*(time-start_rmp))
	endif

! B=B+dB_rmp
	x(:,:,:,6:8)=x(:,:,:,6:8)+b_rmp_out(:,:,:,1:3)*ft_rmp(2)
   	x_8bndx(:,:,6:8)=x_8bndx(:,:,6:8)+b_rmp_bndx(:,:,1:3)*ft_rmp(2)
	x_8bndz(:,:,6:8)=x_8bndz(:,:,6:8)+b_rmp_bndz(:,:,1:3)*ft_rmp(2)
	endif



    call setdt
!mpi   ----------------------------------------------------------------
! here collect all time step dtm and send send minimum dt to all preocess
    call mpi_allreduce(dt1,dt,1,mpi_double_precision,mpi_min, &
                        mpi_comm_world,ierror)
	!$acc update device(dt)
!mpi   ----------------------------------------------------------------
    if((nrank.eq.0).and.(mod(nstep,100).eq.0)) write(*,*) nstep,'t=',time,'dt=',dt
!    if(nrank.eq.0) open(unit=413,file='debug.dat',status='unknown',form='formatted',position='append')
!    if(nrank.eq.0) write(413,*) nstep,'t=',time,'dt=',dt,'ft_rmp=', ft_rmp(2)
!    if(nrank.eq.0) close(413)
    if(dt.lt.1.e-5) goto 900
!     if(nstep .eq. 1)  call recrd1
!     if(nstep .eq. 10)  call recrd10
!     if(nstep .eq. 100)  call recrd100
    if(mod(nstep,nrcd).eq.0) then
	!$acc update host(x)
    call recrd


    if(lcd .gt. 0) call recrd_cud
    if(nrank.eq.0) write(18,*) nst,ncase,nstep,time
	nst=nst+1
    endif

    if(mod(nstep,ndgn).eq.0) then
    do jdgn=1,mdgn
!    call diagn_nmmode(ef(:,:,:,2),jdgn)
    enddo

!    call diagn
!    call diagn_max 
    call diagn_maxmin 
    call diagnatxmode               
    endif

!energy-----------------------------------------------------------------

!      if(mod(nstep,neng).eq.0) then
!      if((mod(nstep,neng).eq.0).and.(nstep.gt.neng)) then
!      
!      call mpi_reduce(ft,ft1,1,mpi_double_precision,mpi_sum,0, &
!                       mpi_comm_world,ierror)
!      call mpi_reduce(gt,gt1,1,mpi_double_precision,mpi_sum,0, &
!                       mpi_comm_world,ierror)
!      call mpi_reduce(ht,ht1,1,mpi_double_precision,mpi_sum,0, &
!                       mpi_comm_world,ierror)
!!mpi   -----------------------------------------------------------------
!      if(nrank.eq.0) then
!      if(nstep.ne.0) then
!      g_rate=(gt1-gtold)/(gt1+gtold)/(time-timeold)
!      gtold=gt1
!      timeold=time
!!      open(unit=12,file='energy.dat',status='unknown',form='formatted',position='append') 
!      !write(12,400) time,ft1-ft0,gt1-gt0,ht1-ht0,g_rate
!400   format(5(1x,e12.5))
!      endif
!      if(nstep.eq.0) then
!      ft0=ft1
!      gt0=gt1
!      ht0=ht1
!      endif
!
!      endif          
!      endif
    
    if(soundwave) then
    call stepon_atfs
    else
!    call stepon
!	call stepon_cut_cell
	call stepon_openacc
!      print*, 'success 853',nstep
    endif
!cd_oxpoint3:ey
    if(cd_oxpoint .and. lcdox==3 .and. mod(nstep,ncd).eq.0 ) call distribution_cd_oxpoint(ef(:,:,3,2))
    
!      print*, 'success 859',nstep
    nstep=nstep+1    
    if(nstep.gt.nstop) goto 900
    goto 100

900 close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    close(18)
    close(111)
    close(112)
    close(211)
    close(311)
    close(413)
    if(nrank==0) then
    timeend=mpi_wtime()
    write(*,*) 'run time=',timeend-timestart
    endif
!clt-k using petsc
    !   call petsc_clean_var
     !  call petscfinalize(ierr)

	!$acc exit data delete(wfx1,wfx2,wfz1,wfz2,wfy1,wfy2,wcx1,wcx2,wcz1,wcz2,wcy1,wcy2)
!mpi   ------------------------------------------------------------------
    call mpi_comm_free(mycomm_y,ierror)
! end mpi
    call mpi_finalize(ierror)
!mpi   -----------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! finish the loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    stop
    end
