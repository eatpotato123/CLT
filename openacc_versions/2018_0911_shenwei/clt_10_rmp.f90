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
!!!!!      !!!!!!!   !!!!!!!!!!   !!!!!!!      !!!!!!!!!!!!!!
!!!!!  !!!  !!!!!!    !!!!!!!  !  !!!!!!!  !!!!  !!!!!!!!!!!!
!!!!!  !!!! !!!!!!  !  !!!!!  !!  !!!!!!!  !!!!! !!!!!!!!!!!!
!!!!!  !!  !!!!!!!  !!!  !! !!!!  !!!!!!!  !!!! !!!!!!!!!!!!!
!!!!!  !!! !!!!!!!  !!!!!  !!!!!  !!!!!!!      !!!!!!!!!!!!!!
!!!!!  !!!! !!!!!!  !!!!!!!!!!!!  !!!!!!!  !!!!!!!!!!!!!!!!!!
!!!!!  !!!!! !!!!!  !!!!!!!!!!!!  !!!!!!!  !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module for RMP coils in EAST tokamak with real parameters.
! the RMP field is calculated with the vector-potential (V-P) A.
! the RMP field is applied only at the SOL region outside the LCS.
! including SOL region, use fixed BCs close to the first wall.
! including the x-point.
! added by H.W. Zhang during 2018 Mar-Jul.
! contact changhw@zju.edu.cn if there is any problem.

! including subroutines as follows          -|
! rmp_east_pert_withA                        |
! -calculate the A for RMP coils           |> these three subroutines used with Div(B_rmp).eq.0
! convt_At_rmp                               |
! -convert V-P A to B_rmp                  |
! map_Bt_rmp_to_bnd_grd                     -|
! -map B_rmp to the boundary points

! rmp_east_pert          -|
! -calculate the B_rmp for RMP coils        |
! rmp_east_pert_bndx                          |> these three subroutines not used for Div(B_rmp).ne.0
! -calculate the B_rmp for boundary points  |
! rmp_east_pert_bndz                          |
! -calculate the B_rmp for boundary points -|


!hwzhang*****************************************************************************
!hwzhang*****************************************************************************
!subroutines for RMP-EAST
!these three subroutines used with Div(B_rmp).eq.0
!hwzhang*****************************************************************************
!hwzhang*****************************************************************************

!hw**************************************************************
	subroutine rmp_east_pert_withA
      use declare
	implicit none
	integer i, j, n_up, n_low, nphi, nab, k
	real*8, dimension(8,2) :: rmp_set
	real*8, dimension(mxt,mzt,my) :: Ax_rmp, Ay_rmp, Az_rmp
	real*8 i0, phi_m, phi_p, raa, zaa, rbb, zbb, rcc, zcc, rdd, zdd, ii0
	real*8 ra1,za1,rb1,zb1,rc1,zc1,rd1,zd1
	real*8 dl, dlx, dly, dlz, rx, ry, rz, phi, dphi, r0, phi0, z0
	real*8 u0, b0_norm, a0_norm, lab, l, rd
      include 'mpif.h'
	
	u0=4.d0*pi*1.d-7
!	b0_norm=2.2745d0 ! T
	b0_norm=2.246028872496384d0 ! T
	a0_norm=1.d0 ! m
!	i0=100 ! A
	i0=i0_rmp
	n_up=8
	n_low=8
	nphi=100
	nab=20

! M Jia et al 2016 Plasma Phys. Control. Fusion 58 055010 
	ra1=2.092
	za1=0.759
	rb1=2.278
	zb1=0.577
	rc1=2.278
	zc1=-0.577
	rd1=2.092
	zd1=-0.759

!	lab=sqrt((raa-rbb)**2+(zaa-zbb)**2)


	Ax_rmp(:,:,:)=0.
	Ay_rmp(:,:,:)=0.
	Az_rmp(:,:,:)=0.
	
	rmp_set(:,1)=(/1., 1., 1., 1., -1., -1., -1., -1./)
	rmp_set(:,2)=(/1., 1., 1., 1., -1., -1., -1., -1./)
!	rmp_set(:,1)=(/4.6, 11.3, 11.3, 4.6, -4.6, -11.3, -11.3, -4.6/)
!	rmp_set(:,2)=(/4.6, 11.3, 11.3, 4.6, -4.6, -11.3, -11.3, -4.6/) ! phi=0
!	rmp_set(:,2)=(/-4.6, 4.6, 11.3, 11.3, 4.6, -4.6, -11.3, -11.3/) ! phi=45
!	rmp_set(:,2)=(/-11.3, -4.6, 4.6, 11.3, 11.3, 4.6, -4.6, -11.3/) ! phi=90
!	rmp_set(:,2)=(/-11.3, -11.3, -4.6, 4.6, 11.3, 11.3, 4.6, -4.6/) ! phi=135
!	rmp_set(:,2)=(/-4.6, -11.3, -11.3, -4.6, 4.6, 11.3, 11.3, 4.6/) ! phi=180
!	rmp_set(:,2)=(/4.6, -4.6, -11.3, -11.3, -4.6, 4.6, 11.3, 11.3/) ! phi=225
!	rmp_set(:,2)=(/11.3, 4.6, -4.6, -11.3, -11.3, -4.6, 4.6, 11.3/) ! phi=270
!	rmp_set(:,2)=(/11.3, 11.3, 4.6, -4.6, -11.3, -11.3, -4.6, 4.6/) ! phi=315
!	rmp_set(:,2)=(/4.6, 11.3, 11.3, 4.6, -4.6, -11.3, -11.3, -4.6/) ! phi=360


! for the grid point ***************************
	do 11 k=1,2
	if(k.eq.1) then ! for the A-B upper
		  raa=ra1
		  zaa=za1
		  rbb=rb1
		  zbb=zb1
	endif
	if(k.eq.2) then ! for the C-D lower
		  raa=rc1
		  zaa=zc1
		  rbb=rd1
		  zbb=zd1
	endif

	lab=sqrt((raa-rbb)**2+(zaa-zbb)**2)

	do 1 jx=1,mxt
	do 1 jz=1,mzt
	do 1 jy=iy_first,iy_last

	r0=xxt(jx)*a0_norm
	phi0=yy(jy)
	z0=zzt(jz)*a0_norm


	do 2 i=1,n_up

	! for coil A-B upper(k=2 is for C-D lower)

	if(k.eq.1) then ! for the A-B upper
	phi_m=(rmp_phi_up+360.d0/n_up*(i*1.d0-1.d0)-18.5d0)/180.d0*pi
	phi_p=(rmp_phi_up+360.d0/n_up*(i*1.d0-1.d0)+18.5d0)/180.d0*pi
	dphi=(phi_p-phi_m)/nphi
	endif
	if(k.eq.2) then ! for the C-D lower
	phi_m=(rmp_phi_low+360.d0/n_up*(i*1.d0-1.d0)-18.5d0)/180.d0*pi
	phi_p=(rmp_phi_low+360.d0/n_up*(i*1.d0-1.d0)+18.5d0)/180.d0*pi
	dphi=(phi_p-phi_m)/nphi
	endif

	do 3 j=1,nphi
	phi=phi_m-0.5*dphi+j*dphi
	! for A toroidal part
	ii0=rmp_set(i,k)*i0
	rx=r0*cos(phi0)-raa*cos(phi)
	ry=r0*sin(phi0)-raa*sin(phi)
	rz=z0-zaa
	rd=sqrt(rx**2+ry**2+rz**2)

	dlx=-raa*sin(phi)*dphi
	dly=raa*cos(phi)*dphi
	dlz=0.

	Ax_rmp(jx,jz,jy)=Ax_rmp(jx,jz,jy)+u0*ii0/4./pi*1./rd*dlx
	Ay_rmp(jx,jz,jy)=Ay_rmp(jx,jz,jy)+u0*ii0/4./pi*1./rd*dly
	Az_rmp(jx,jz,jy)=Az_rmp(jx,jz,jy)+u0*ii0/4./pi*1./rd*dlz

	! for B toroidal part
	ii0=-rmp_set(i,k)*i0
	rx=r0*cos(phi0)-rbb*cos(phi)
	ry=r0*sin(phi0)-rbb*sin(phi)
	rz=z0-zbb
	rd=sqrt(rx**2+ry**2+rz**2)

	dlx=-rbb*sin(phi)*dphi
	dly=rbb*cos(phi)*dphi
	dlz=0.

	Ax_rmp(jx,jz,jy)=Ax_rmp(jx,jz,jy)+u0*ii0/4./pi*1./rd*dlx
	Ay_rmp(jx,jz,jy)=Ay_rmp(jx,jz,jy)+u0*ii0/4./pi*1./rd*dly
	Az_rmp(jx,jz,jy)=Az_rmp(jx,jz,jy)+u0*ii0/4./pi*1./rd*dlz
    3 continue

	! for left and right poroidal part, phi=phi_m or phi_p
	dl=lab/nab

	do 4 j=1,nab
	l=0-0.5*dl+j*dl
	! for left poroidal part
	ii0=-rmp_set(i,k)*i0
	phi=phi_m
	rx=r0*cos(phi0)-(raa+(rbb-raa)/lab*l)*cos(phi)
	ry=r0*sin(phi0)-(raa+(rbb-raa)/lab*l)*sin(phi)
	rz=z0-zaa-(zbb-zaa)/lab*l
	rd=sqrt(rx**2+ry**2+rz**2)

	dlx=(rbb-raa)*cos(phi)*dl/lab
	dly=(rbb-raa)*sin(phi)*dl/lab
	dlz=(zbb-zaa)*dl/lab

	Ax_rmp(jx,jz,jy)=Ax_rmp(jx,jz,jy)+u0*ii0/4./pi*1./rd*dlx
	Ay_rmp(jx,jz,jy)=Ay_rmp(jx,jz,jy)+u0*ii0/4./pi*1./rd*dly
	Az_rmp(jx,jz,jy)=Az_rmp(jx,jz,jy)+u0*ii0/4./pi*1./rd*dlz
	

	! for right poroidal part
	ii0=rmp_set(i,k)*i0
	phi=phi_p
	rx=r0*cos(phi0)-(raa+(rbb-raa)/lab*l)*cos(phi)
	ry=r0*sin(phi0)-(raa+(rbb-raa)/lab*l)*sin(phi)
	rz=z0-zaa-(zbb-zaa)/lab*l
	rd=sqrt(rx**2+ry**2+rz**2)

	dlx=(rbb-raa)*cos(phi)*dl/lab
	dly=(rbb-raa)*sin(phi)*dl/lab
	dlz=(zbb-zaa)*dl/lab

	Ax_rmp(jx,jz,jy)=Ax_rmp(jx,jz,jy)+u0*ii0/4./pi*1./rd*dlx
	Ay_rmp(jx,jz,jy)=Ay_rmp(jx,jz,jy)+u0*ii0/4./pi*1./rd*dly
	Az_rmp(jx,jz,jy)=Az_rmp(jx,jz,jy)+u0*ii0/4./pi*1./rd*dlz
    4 continue

    2 continue

    1 continue

   11 continue


	do 12 jx=1,mxt
	do 12 jz=1,mzt
	do 12 jy=iy_first,iy_last
	phi0=yy(jy)
	! normalization and transfer coord
	At_rmp(jx,jz,jy,1)=Ax_rmp(jx,jz,jy)*cos(phi0)+Ay_rmp(jx,jz,jy)*sin(phi0)
	At_rmp(jx,jz,jy,2)=Ay_rmp(jx,jz,jy)*cos(phi0)-Ax_rmp(jx,jz,jy)*sin(phi0)
	At_rmp(jx,jz,jy,3)=Az_rmp(jx,jz,jy)


!	At_rmp(jx,jz,jy,:)=At_rmp(jx,jz,jy,:)*(1.+tanh((pst(jx,jz)-1.07)*50.))/2.
	At_rmp(jx,jz,jy,:)=At_rmp(jx,jz,jy,:)/a0_norm/b0_norm
! this method cannot be used below, for At_rmp is for total Box and hypb_raito
! is for sub-box.
!	hypb_value=2.984926524705880E-002
!	At_rmp(jx,jz,jy,3)=At_rmp(jx,jz,jy,;)*(1.-hypb_ratio(jx,jz))
!	At_rmp(jx,jz,jy,3)=At_rmp(jx,jz,jy,;)*(1.-max(0.d0,dtanh((dist_to_bnd(jx,jz)-4.d-2)*30.d0)))

   12 continue

!haowei:  calculate Bt_rmp=\cross A_rmp
	call convt_At_rmp
	do jx=1,mxt
	do jz=1,mzt
	do jy=iy_first,iy_last
	Bt_rmp(jx,jz,jy,1)=dAy_rmp(jx,jz,jy,3)/xxt(jx)-dAz_rmp(jx,jz,jy,2)
	Bt_rmp(jx,jz,jy,2)=dAz_rmp(jx,jz,jy,1)-dAx_rmp(jx,jz,jy,3)
	Bt_rmp(jx,jz,jy,3)=At_rmp(jx,jz,jy,2)/xxt(jx)+dAx_rmp(jx,jz,jy,2)-dAy_rmp(jx,jz,jy,1)/xxt(jx)
!	Bt_rmp(jx,jz,jy,1)=At_rmp(jx,jz,jy,1)
!	Bt_rmp(jx,jz,jy,2)=At_rmp(jx,jz,jy,2)
!	Bt_rmp(jx,jz,jy,3)=At_rmp(jx,jz,jy,3)

	enddo
	enddo
	enddo

	do jy=iy_first,iy_last
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
      b_rmp(jx,jz,jy,1)=Bt_rmp(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2,jy,1)
      b_rmp(jx,jz,jy,2)=Bt_rmp(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2,jy,2)
      b_rmp(jx,jz,jy,3)=Bt_rmp(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2,jy,3)
	enddo
	enddo
	enddo

! calculate the b_rmp_outside, the fixed part
	do jx=1,mxt
	do jz=1,mzt
	do jy=iy_first,iy_last
	At_rmp(jx,jz,jy,:)=At_rmp(jx,jz,jy,:)*(1.+tanh((pst(jx,jz)-1.09)*50.))/2.
	enddo
	enddo
	enddo

!haowei:  calculate Bt_rmp=\cross A_rmp
	call convt_At_rmp
	do jx=1,mxt
	do jz=1,mzt
	do jy=iy_first,iy_last
	Bt_rmp(jx,jz,jy,1)=dAy_rmp(jx,jz,jy,3)/xxt(jx)-dAz_rmp(jx,jz,jy,2)
	Bt_rmp(jx,jz,jy,2)=dAz_rmp(jx,jz,jy,1)-dAx_rmp(jx,jz,jy,3)
	Bt_rmp(jx,jz,jy,3)=At_rmp(jx,jz,jy,2)/xxt(jx)+dAx_rmp(jx,jz,jy,2)-dAy_rmp(jx,jz,jy,1)/xxt(jx)

	enddo
	enddo
	enddo

	do jy=iy_first,iy_last
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
      b_rmp_out(jx,jz,jy,1)=Bt_rmp(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2,jy,1)
      b_rmp_out(jx,jz,jy,2)=Bt_rmp(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2,jy,2)
      b_rmp_out(jx,jz,jy,3)=Bt_rmp(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2,jy,3)
	enddo
	enddo
	enddo

	call recrd_rmp_vacuum ! nst=0 is used to recrd the vacuum rmp field
	call recrd_rmp_vacuum_2 ! nst=0 is used to recrd the vacuum rmp field

! keep only b_rmp_out part
!	b_rmp=b_rmp_out
	b_rmp_in(:,:,:,:)=b_rmp(:,:,:,:)-b_rmp_out(:,:,:,:)

!    	x(:,:,:,6:8)=x(:,:,:,6:8)+b_rmp(:,:,:,:)
!      do jy=iy_first,iy_last
!      x1(:,:,jy,:)=x(:,:,jy,:)-xint(:,:,:)
!      enddo


	call map_Bt_rmp_to_bnd_grd
!	call rmp_east_pert_bndx(u0, b0_norm, a0_norm, i0, n_up, n_low, nphi, nab,&
!	ra1, za1, rb1, zb1, rc1, zc1, rd1, zd1, rmp_set)
!	call rmp_east_pert_bndz(u0, b0_norm, a0_norm, i0, n_up, n_low, nphi, nab,&
!	ra1, za1, rb1, zb1, rc1, zc1, rd1, zd1, rmp_set)

!	if(nrank.eq.0) then
!	open(unit=1,file='bx_rmp_bndz.dat',status='unknown',form='formatted')
!	open(unit=2,file='bx_rmp_bndx.dat',status='unknown',form='formatted')
!	write(1,33) ((bndz_grd(jx,1),bndz_grd(jx,2),b_rmp_bndz(jx,1,1)),jx=1,nbndz)
!	write(2,33) ((bndx_grd(jx,1),bndx_grd(jx,2),b_rmp_bndx(jx,1,1)),jx=1,nbndx)
!   33 format(3(1x,e12.5))
!    	close(1)
!    	close(2)
!	endif
      
    	return
	end

!hw***************************************************************************************
	subroutine convt_At_rmp
	use declare
!	implicit none
	real*8, dimension(mxt,mzt,3) :: wsy1,wsy2
	integer, parameter :: mxzt=mxt*mzt
      include 'mpif.h'
!  d1fc= d f / dx  with fourth-order accuracy central difference
      d1fc(fm2,fm1,f0,fp1,fp2,a,b,c,d)= &
       a*(fp1-f0)+b*(f0-fm1)+c*(fp2-f0)+d*(f0-fm2)
!  d1fp= d f / dx  with  one-sided  difference involving 0  1 2 and 3
!  points
      d1fp(fp3,fp2,fp1,f0,a,b,c)= &
       a*(fp1-f0)+b*(fp2-f0)+c*(fp3-f0)
!  d1fm= d f / dx  with one-sided difference involving -3 -2 -1 and 0
!  points
      d1fm(fm3,fm2,fm1,f0,a,b,c)= &
       a*(f0-fm1)+b*(f0-fm2)+c*(f0-fm3)
!  d1fbp= d f / dx  with  one-sided-bias  difference involving -1 0  1 and 2
!  points
      d1fbp(fp2,fp1,f0,fm1,a,b,c)= &
       a*(fm1-f0)+b*(fp1-f0)+c*(fp2-f0)
!  d1fbm= d f / dx  with  one-sided-bias  difference involving -2 -1  0 and 1
!  points
      d1fbm(fm2,fm1,f0,fp1,a,b,c)= &
       a*(f0-fp1)+b*(f0-fm1)+c*(f0-fm2)

	do 1 m=1,3

	do 2 jy=iy_first,iy_last
	do 2 jz=1,mzt

	do jx=3,mxt-2
	dAx_rmp(jx,jz,jy,m)=d1fc(At_rmp(jx-2,jz,jy,m),At_rmp(jx-1,jz,jy,m),At_rmp(jx,jz,jy,m),&
		  At_rmp(jx+1,jz,jy,m),At_rmp(jx+2,jz,jy,m),axt1(jx),bxt1(jx),cxt1(jx),dxt1(jx))
	enddo

	dAx_rmp(1,jz,jy,m)=d1fp(At_rmp(4,jz,jy,m),At_rmp(3,jz,jy,m),At_rmp(2,jz,jy,m),&
		  At_rmp(1,jz,jy,m),axtp(1),bxtp(1),cxtp(1))
	dAx_rmp(mxt,jz,jy,m)=d1fm(At_rmp(mxt-3,jz,jy,m),At_rmp(mxt-2,jz,jy,m),At_rmp(mxt-1,jz,jy,m),&
		  At_rmp(mxt,jz,jy,m),axtm(mxt),bxtm(mxt),cxtm(mxt))
	dAx_rmp(2,jz,jy,m)=d1fbp(At_rmp(4,jz,jy,m),At_rmp(3,jz,jy,m),At_rmp(2,jz,jy,m),&
		  At_rmp(1,jz,jy,m),axtbp(2),bxtbp(2),cxtbp(2))
	dAx_rmp(mxt-1,jz,jy,m)=d1fbm(At_rmp(mxt-3,jz,jy,m),At_rmp(mxt-2,jz,jy,m),At_rmp(mxt-1,jz,jy,m),&
		  At_rmp(mxt,jz,jy,m),axtbm(mxt-1),bxtbm(mxt-1),cxtbm(mxt-1))

    2 continue

    	do 3 jy=iy_first,iy_last
	do 3 jx=1,mzt

	do jz=3,mzt-2
	dAz_rmp(jx,jz,jy,m)=d1fc(At_rmp(jx,jz-2,jy,m),At_rmp(jx,jz-1,jy,m),At_rmp(jx,jz,jy,m),&
		  At_rmp(jx,jz+1,jy,m),At_rmp(jx,jz+2,jy,m),azt1(jz),bzt1(jz),czt1(jz),dzt1(jz))
	enddo

	dAz_rmp(jx,1,jy,m)=d1fp(At_rmp(jx,4,jy,m),At_rmp(jx,3,jy,m),At_rmp(jx,2,jy,m),&
		  At_rmp(jx,1,jy,m),aztp(1),bztp(1),cztp(1))
	dAz_rmp(jx,mzt,jy,m)=d1fm(At_rmp(jx,mzt-3,jy,m),At_rmp(jx,mzt-2,jy,m),At_rmp(jx,mzt-1,jy,m),&
		  At_rmp(jx,mzt,jy,m),aztm(mzt),bztm(mzt),cztm(mzt))
	dAz_rmp(jx,2,jy,m)=d1fbp(At_rmp(jx,4,jy,m),At_rmp(jx,3,jy,m),At_rmp(jx,2,jy,m),&
		  At_rmp(jx,1,jy,m),aztbp(2),bztbp(2),cztbp(2))
	dAz_rmp(jx,mzt-1,jy,m)=d1fbm(At_rmp(jx,mzt-3,jy,m),At_rmp(jx,mzt-2,jy,m),At_rmp(jx,mzt-1,jy,m),&
		  At_rmp(jx,mzt,jy,m),aztbm(mzt-1),bztbm(mzt-1),cztbm(mzt-1))

    3 continue

	do 4 jz=1,mzt
	do 4 jx=1,mxt

	do jy=iy_first+2,iy_last-2
	dAy_rmp(jx,jz,jy,m)=d1fc(At_rmp(jx,jz,jy-2,m),At_rmp(jx,jz,jy-1,m),At_rmp(jx,jz,jy,m),&
		  At_rmp(jx,jz,jy+1,m),At_rmp(jx,jz,jy+2,m),ay1(jy),by1(jy),cy1(jy),dy1(jy))
	enddo
 
    4 continue

    1 continue



! send the dAy_rmp to nearest cpu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      if(npry .gt. 1) then
	if (nrky(nrank).lt. npry-1) then
      wsy1(:,:,:)=dAy_rmp(:,:,iy_last-2,:)
      wsy2(:,:,:)=dAy_rmp(:,:,iy_last-3,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wsy1, mxzt*3, mpi_double_precision, nrank + nprxz, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wsy2, mxzt*3, mpi_double_precision, nrank + nprxz, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrky(nrank) .ge. 1 )  then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wsy1, mxzt*3, mpi_double_precision, nrank - nprxz, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wsy2, mxzt*3, mpi_double_precision, nrank - nprxz, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      dAy_rmp(:,:,iy_first+1,:)=wsy1(:,:,:)
      dAy_rmp(:,:,iy_first,:)=wsy2(:,:,:)
	endif


  	if (nrky(nrank).eq. npry-1) then
      wsy1(:,:,:)=dAy_rmp(:,:,iy_last-2,:)
      wsy2(:,:,:)=dAy_rmp(:,:,iy_last-3,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wsy1, mxzt*3, mpi_double_precision, nrank -(npry-1)*nprxz, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wsy2, mxzt*3, mpi_double_precision, nrank -(npry-1)*nprxz, 0,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrky(nrank) .eq. 0 )  then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wsy1, mxzt*3, mpi_double_precision, nrank +(npry-1)*nprxz, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wsy2, mxzt*3, mpi_double_precision, nrank +(npry-1)*nprxz, 0,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      dAy_rmp(:,:,iy_first+1,:)=wsy1(:,:,:)
      dAy_rmp(:,:,iy_first,:)=wsy2(:,:,:)
	endif
	     
! send w8 down unless i'm at the bottom

      
	if (nrky(nrank) .ge. 1 ) then
      wsy1(:,:,:)=dAy_rmp(:,:,iy_first+2,:)
      wsy2(:,:,:)=dAy_rmp(:,:,iy_first+3,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wsy1, mxzt*3, mpi_double_precision, nrank - nprxz, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wsy2, mxzt*3, mpi_double_precision, nrank - nprxz, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrky(nrank).lt. npry-1) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wsy1, mxzt*3, mpi_double_precision, nrank + nprxz, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wsy2, mxzt*3, mpi_double_precision, nrank + nprxz, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      dAy_rmp(:,:,iy_last-1,:)=wsy1(:,:,:)
      dAy_rmp(:,:,iy_last,:)=wsy2(:,:,:)
	endif

    if (nrky(nrank) .eq. 0 ) then
      wsy1(:,:,:)=dAy_rmp(:,:,iy_first+2,:)
      wsy2(:,:,:)=dAy_rmp(:,:,iy_first+3,:)
!mpi   ----------------------------------------------------------------
	call mpi_send( wsy1, mxzt*3, mpi_double_precision, nrank +(npry-1)*nprxz, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wsy2, mxzt*3, mpi_double_precision, nrank +(npry-1)*nprxz, 1,  &
		      mpi_comm_world,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrky(nrank).eq. npry-1) then
!mpi   ----------------------------------------------------------------
	call mpi_recv( wsy1, mxzt*3, mpi_double_precision, nrank -(npry-1)*nprxz, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wsy2, mxzt*3, mpi_double_precision, nrank -(npry-1)*nprxz, 1,  &
		      mpi_comm_world, status,ierror )
!mpi   ----------------------------------------------------------------
      dAy_rmp(:,:,iy_last-1,:)=wsy1(:,:,:)
      dAy_rmp(:,:,iy_last,:)=wsy2(:,:,:)
	endif
      else
      dAy_rmp(:,:,iy_first+1,:)=dAy_rmp(:,:,iy_last-2,:)
      dAy_rmp(:,:,iy_first,:)=dAy_rmp(:,:,iy_last-3,:)
      dAy_rmp(:,:,iy_last-1,:)=dAy_rmp(:,:,iy_first+2,:)
      dAy_rmp(:,:,iy_last,:)=dAy_rmp(:,:,iy_first+3,:)
      endif

	return
	end

!hw***********************************************************************************
	subroutine map_Bt_rmp_to_bnd_grd
	use declare
      integer, parameter :: nq = 33
      integer, parameter :: nr = 150 !int( sqrt(ndat/3) )
!      integer, parameter :: nr = int( sqrt(ndata*1./3.) )
      integer, parameter :: nw = 39
    !
      integer lcell(nr,nr)
      integer lnext(ndata),lnext12(ndat12),lnext34(ndat34)
      real*8 risq(ndata),risq12(ndat12),risq34(ndat34)
      real*8 aw(5,ndata),aw12(5,ndat12),aw34(5,ndat34)
      real*8 rimax,ximin,zimin,dxi,dzi
!      real*8, dimension(mxt,mzt) :: bx_dx,bx_dz,bz_dx,bz_dz,by_dx,by_dz,uy_dx,uy_dz,tht_dx,tht_dz
!      real*8, dimension(mxt,mzt) :: pt_dx,pt_dz,rh_dx,rh_dz
!	real*8, dimension(nbndx) :: bx_bndx_dx, bx_bndx_dz
!      real*8, dimension(mxt,mzt) :: bx,bxdx,bxdz,bz,bzdx,bzdz
!      real*8, dimension(mxt,mzt) :: bxdx_dx,bxdx_dz,bxdz_dx,bxdz_dz,bzdx_dx,bzdx_dz,bzdz_dx,bzdz_dz
!      integer icell(1,1)
!      integer inext(9)
!      real*8 xin(9),zin(9),qin(9),rinsq(9),ain(5,9)
!       real*8 xout,zout,qout,qdxout,qdzout
       integer iam,iap, itag
	 integer ier

	 real*8, dimension(ndata) :: tmp_int, xx_int, zz_int
	 real*8, dimension(nbndx) :: tmp_bndx_out1, tmp_bndx_out2, tmp_bndx_out3
	 real*8, dimension(nbndz) :: tmp_bndz_out1, tmp_bndz_out2, tmp_bndz_out3
	 real*8, dimension(mxt,mzt) :: xxt2d, zzt2d
      include 'mpif.h'
!  d1fc= d f / dx  with fourth-order accuracy central difference
      d1fc(fm2,fm1,f0,fp1,fp2,a,b,c,d)= &
       a*(fp1-f0)+b*(f0-fm1)+c*(fp2-f0)+d*(f0-fm2)   
	
	allocate(b_rmp_bndx(nbndx,my,3))
	allocate(b_rmp_bndz(nbndz,my,3))
	b_rmp_bndx(:,:,:)=0.
	b_rmp_bndz(:,:,:)=0.

! interpolation, input ndata, xx_int, zz_int, f_int, 
! output f(jx,jz), and its x & z directional derivatives at xxt and zzt grids.
	if(firstmap) then

	do jx=1,mxt
	do jz=1,mzt
	xxt2d(jx,jz)=xxt(jx)
	zzt2d(jx,jz)=zzt(jz)
	enddo
	enddo

	xx_int=reshape(xxt2d, [ndata])
	zz_int=reshape(zzt2d, [ndata])

	do 1 jy=iy_first,iy_last
	do 1 m=1,3
	tmp_int(:)=reshape(Bt_rmp(:,:,jy,m), [ndata])

      call qshep2 ( ndata, xx_int, zz_int, tmp_int, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      do jx=1,nbndx
      call qs2grd ( bndx_grd(jx,1), bndx_grd(jx,2), ndata, xx_int, zz_int, tmp_int, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, tmp_bndx_out1(jx), tmp_bndx_out2(jx), tmp_bndx_out3(jx), ier )
      write(*,*) 'itag=', itag, 'jx=', jx, ier
      enddo
      do jz=1,nbndz
      call qs2grd ( bndz_grd(jz,1), bndz_grd(jz,2), ndata, xx_int, zz_int, tmp_int, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, tmp_bndz_out1(jz), tmp_bndz_out2(jz), tmp_bndz_out3(jz), ier )
      write(*,*) 'itag=', itag, 'jz=', jz, ier
      enddo

	b_rmp_bndx(:,jy,m)=tmp_bndx_out1(:)
	b_rmp_bndz(:,jy,m)=tmp_bndz_out1(:)

    1 continue
    	endif

!   	x_8bndx(:,:,6:8)=x_8bndx(:,:,6:8)+b_rmp_bndx(:,:,1:3)
!   	x_8bndz(:,:,6:8)=x_8bndz(:,:,6:8)+b_rmp_bndz(:,:,1:3)


	if(nrank.eq.0) then
	open(unit=1,file='bx_rmp_bndz.dat',status='unknown',form='formatted')
	open(unit=2,file='bx_rmp_bndx.dat',status='unknown',form='formatted')
	write(1,3) (bndz_grd(jx,1),bndz_grd(jx,2),b_rmp_bndz(jx,3,1),jx=1,nbndz)
	write(2,3) (bndx_grd(jx,1),bndx_grd(jx,2),b_rmp_bndx(jx,3,1),jx=1,nbndx)
    3 format(3(1x,e12.5))
    	close(1)
    	close(2)
	endif
    	
    	return
	end


!hwzhang*****************************************************************************
!hwzhang*****************************************************************************
!subroutines for RMP-EAST
!these three subroutines not used for Div(B_rmp).ne.0
!hwzhang*****************************************************************************
!hwzhang*****************************************************************************

!hw**************************************************************
	subroutine rmp_east_pert
      use declare
	implicit none
	integer i, j, n_up, n_low, nphi, nab, k
	real*8, dimension(8,2) :: rmp_set
	real*8 i0, phi_m, phi_p, raa, zaa, rbb, zbb, rcc, zcc, rdd, zdd, ii0
	real*8 ra1,za1,rb1,zb1,rc1,zc1,rd1,zd1
	real*8 dl, dlx, dly, dlz, rx, ry, rz, phi, dphi, r0, phi0, z0
	real*8 u0, b0_norm, a0_norm, lab, l, rd
      include 'mpif.h'
	
	u0=4.d0*pi*1.d-7
!	b0_norm=2.2745d0 ! T
	b0_norm=2.246028872496384d0 ! T
	a0_norm=1.d0 ! m
!	i0=100 ! A
	i0=i0_rmp ! A
	n_up=8
	n_low=8
	nphi=100
	nab=20

! M Jia et al 2016 Plasma Phys. Control. Fusion 58 055010 
	ra1=2.092
	za1=0.759
	rb1=2.278
	zb1=0.577
	rc1=2.278
	zc1=-0.577
	rd1=2.092
	zd1=-0.759

!	lab=sqrt((raa-rbb)**2+(zaa-zbb)**2)


	bx_rmp(:,:,:)=0
	by_rmp(:,:,:)=0
	bz_rmp(:,:,:)=0
	
	rmp_set(:,1)=(/1., 1., 1., 1., -1., -1., -1., -1./)
	rmp_set(:,2)=(/1., 1., 1., 1., -1., -1., -1., -1./)
!	rmp_set(:,1)=(/4.6, 11.3, 11.3, 4.6, -4.6, -11.3, -11.3, -4.6/)
!	rmp_set(:,2)=(/4.6, 11.3, 11.3, 4.6, -4.6, -11.3, -11.3, -4.6/) ! phi=0
!	rmp_set(:,2)=(/-4.6, 4.6, 11.3, 11.3, 4.6, -4.6, -11.3, -11.3/) ! phi=45
!	rmp_set(:,2)=(/-11.3, -4.6, 4.6, 11.3, 11.3, 4.6, -4.6, -11.3/) ! phi=90
!	rmp_set(:,2)=(/-11.3, -11.3, -4.6, 4.6, 11.3, 11.3, 4.6, -4.6/) ! phi=135
!	rmp_set(:,2)=(/-4.6, -11.3, -11.3, -4.6, 4.6, 11.3, 11.3, 4.6/) ! phi=180
!	rmp_set(:,2)=(/4.6, -4.6, -11.3, -11.3, -4.6, 4.6, 11.3, 11.3/) ! phi=225
!	rmp_set(:,2)=(/11.3, 4.6, -4.6, -11.3, -11.3, -4.6, 4.6, 11.3/) ! phi=270
!	rmp_set(:,2)=(/11.3, 11.3, 4.6, -4.6, -11.3, -11.3, -4.6, 4.6/) ! phi=315
!	rmp_set(:,2)=(/4.6, 11.3, 11.3, 4.6, -4.6, -11.3, -11.3, -4.6/) ! phi=360

! for the grid point ***************************
	do 11 k=1,2
	if(k.eq.1) then ! for the A-B upper
		  raa=ra1
		  zaa=za1
		  rbb=rb1
		  zbb=zb1
	endif
	if(k.eq.2) then ! for the C-D lower
		  raa=rc1
		  zaa=zc1
		  rbb=rd1
		  zbb=zd1
	endif

	lab=sqrt((raa-rbb)**2+(zaa-zbb)**2)

	do 1 jx=ix_first,ix_last
	do 1 jz=iz_first,iz_last
	do 1 jy=iy_first,iy_last

	r0=xx(jx)*a0_norm
	phi0=yy(jy)
	z0=zz(jz)*a0_norm


	do 2 i=1,n_up

	! for coil A-B upper(k=2 is for C-D lower)

	if(k.eq.1) then ! for the A-B upper
	phi_m=(rmp_phi_up+360.d0/n_up*(i*1.d0-1.d0)-18.5d0)/180.d0*pi
	phi_p=(rmp_phi_up+360.d0/n_up*(i*1.d0-1.d0)+18.5d0)/180.d0*pi
	dphi=(phi_p-phi_m)/nphi
	endif
	if(k.eq.2) then ! for the C-D lower
	phi_m=(rmp_phi_low+360.d0/n_up*(i*1.d0-1.d0)-18.5d0)/180.d0*pi
	phi_p=(rmp_phi_low+360.d0/n_up*(i*1.d0-1.d0)+18.5d0)/180.d0*pi
	dphi=(phi_p-phi_m)/nphi
	endif


	do 3 j=1,nphi
	phi=phi_m-0.5*dphi+j*dphi
	! for A toroidal part
	ii0=rmp_set(i,k)*i0
	rx=r0*cos(phi0)-raa*cos(phi)
	ry=r0*sin(phi0)-raa*sin(phi)
	rz=z0-zaa
	rd=sqrt(rx**2+ry**2+rz**2)

	dlx=-raa*sin(phi)*dphi
	dly=raa*cos(phi)*dphi
	dlz=0.

	bx_rmp(jx,jz,jy)=bx_rmp(jx,jz,jy)+u0*ii0/4./pi*1./rd**3*(dly*rz-dlz*ry)
	by_rmp(jx,jz,jy)=by_rmp(jx,jz,jy)+u0*ii0/4./pi*1./rd**3*(dlz*rx-dlx*rz)
	bz_rmp(jx,jz,jy)=bz_rmp(jx,jz,jy)+u0*ii0/4./pi*1./rd**3*(dlx*ry-dly*rx)

	! for B toroidal part
	ii0=-rmp_set(i,k)*i0
	rx=r0*cos(phi0)-rbb*cos(phi)
	ry=r0*sin(phi0)-rbb*sin(phi)
	rz=z0-zbb
	rd=sqrt(rx**2+ry**2+rz**2)

	dlx=-rbb*sin(phi)*dphi
	dly=rbb*cos(phi)*dphi
	dlz=0.

	bx_rmp(jx,jz,jy)=bx_rmp(jx,jz,jy)+u0*ii0/4./pi*1./rd**3*(dly*rz-dlz*ry)
	by_rmp(jx,jz,jy)=by_rmp(jx,jz,jy)+u0*ii0/4./pi*1./rd**3*(dlz*rx-dlx*rz)
	bz_rmp(jx,jz,jy)=bz_rmp(jx,jz,jy)+u0*ii0/4./pi*1./rd**3*(dlx*ry-dly*rx)
    3 continue

	! for left and right poroidal part, phi=phi_m or phi_p
	dl=lab/nab

	do 4 j=1,nab
	l=0-0.5*dl+j*dl
	! for left poroidal part
	ii0=-rmp_set(i,k)*i0
	phi=phi_m
	rx=r0*cos(phi0)-(raa+(rbb-raa)/lab*l)*cos(phi)
	ry=r0*sin(phi0)-(raa+(rbb-raa)/lab*l)*sin(phi)
	rz=z0-zaa-(zbb-zaa)/lab*l
	rd=sqrt(rx**2+ry**2+rz**2)

	dlx=(rbb-raa)*cos(phi)*dl/lab
	dly=(rbb-raa)*sin(phi)*dl/lab
	dlz=(zbb-zaa)*dl/lab

	bx_rmp(jx,jz,jy)=bx_rmp(jx,jz,jy)+u0*ii0/4./pi*1./rd**3*(dly*rz-dlz*ry)
	by_rmp(jx,jz,jy)=by_rmp(jx,jz,jy)+u0*ii0/4./pi*1./rd**3*(dlz*rx-dlx*rz)
	bz_rmp(jx,jz,jy)=bz_rmp(jx,jz,jy)+u0*ii0/4./pi*1./rd**3*(dlx*ry-dly*rx)
	

	! for right poroidal part
	ii0=rmp_set(i,k)*i0
	phi=phi_p
	rx=r0*cos(phi0)-(raa+(rbb-raa)/lab*l)*cos(phi)
	ry=r0*sin(phi0)-(raa+(rbb-raa)/lab*l)*sin(phi)
	rz=z0-zaa-(zbb-zaa)/lab*l
	rd=sqrt(rx**2+ry**2+rz**2)

	dlx=(rbb-raa)*cos(phi)*dl/lab
	dly=(rbb-raa)*sin(phi)*dl/lab
	dlz=(zbb-zaa)*dl/lab

	bx_rmp(jx,jz,jy)=bx_rmp(jx,jz,jy)+u0*ii0/4./pi*1./rd**3*(dly*rz-dlz*ry)
	by_rmp(jx,jz,jy)=by_rmp(jx,jz,jy)+u0*ii0/4./pi*1./rd**3*(dlz*rx-dlx*rz)
	bz_rmp(jx,jz,jy)=bz_rmp(jx,jz,jy)+u0*ii0/4./pi*1./rd**3*(dlx*ry-dly*rx)
    4 continue

    2 continue

    1 continue

   11 continue


	do 12 jx=ix_first,ix_last
	do 12 jz=iz_first,iz_last
	do 12 jy=iy_first,iy_last
	phi0=yy(jy)
	! normalization and transfer coord
	b_rmp(jx,jz,jy,1)=bx_rmp(jx,jz,jy)*cos(phi0)+by_rmp(jx,jz,jy)*sin(phi0)
	b_rmp(jx,jz,jy,2)=by_rmp(jx,jz,jy)*cos(phi0)-bx_rmp(jx,jz,jy)*sin(phi0)
	b_rmp(jx,jz,jy,3)=bz_rmp(jx,jz,jy)

	b_rmp(jx,jz,jy,:)=b_rmp(jx,jz,jy,:)/b0_norm
!	b_rmp(jx,jz,jy,:)=b_rmp(jx,jz,jy,:)*(1.+tanh((psi(jx,jz)-psia)*50.))/2.

   12 continue

    	x(:,:,:,6:8)=x(:,:,:,6:8)+b_rmp(:,:,:,:)
!      do jy=iy_first,iy_last
!      x1(:,:,jy,:)=x(:,:,jy,:)-xint(:,:,:)
!      enddo


	call rmp_east_pert_bndx(u0, b0_norm, a0_norm, i0, n_up, n_low, nphi, nab,&
	ra1, za1, rb1, zb1, rc1, zc1, rd1, zd1, rmp_set)
	call rmp_east_pert_bndz(u0, b0_norm, a0_norm, i0, n_up, n_low, nphi, nab,&
	ra1, za1, rb1, zb1, rc1, zc1, rd1, zd1, rmp_set)
      
    	return
	end

!hw**************************************************************
	subroutine rmp_east_pert_bndx(u0, b0_norm, a0_norm, i0, n_up, n_low, nphi, nab,&
	ra1, za1, rb1, zb1, rc1, zc1, rd1, zd1, rmp_set)
      use declare
	implicit none
	integer i, j, n_up, n_low, nphi, nab, k
	real*8, dimension(8,2) :: rmp_set
	real*8 i0, phi_m, phi_p, raa, zaa, rbb, zbb, rcc, zcc, rdd, zdd, ii0
	real*8 ra1,za1,rb1,zb1,rc1,zc1,rd1,zd1
	real*8 dl, dlx, dly, dlz, rx, ry, rz, phi, dphi, r0, phi0, z0
	real*8 u0, b0_norm, a0_norm, lab, l, rd
      include 'mpif.h'

	allocate(b_rmp_bndx(nbndx,my,3))
	
!	u0=4.d0*pi*1.d-7
!	b0_norm=2.2745d0 ! T
!	a0_norm=0.8 ! m
!	i0=4000 ! A
!	n_up=8
!	n_low=8
!	nphi=100
!	nab=20

!	ra1=2.095
!	za1=0.765
!	rb1=2.28
!	zb1=0.586
!	rc1=2.28
!	zc1=-0.586
!	rd1=2.095
!	zd1=-0.765

!	lab=sqrt((raa-rbb)**2+(zaa-zbb)**2)


	b_rmp_bndx(:,:,:)=0.
	
!	rmp_set(:,1)=(/1, -1, 1, -1, 1, -1, 1, -1/)
!	rmp_set(:,2)=(/1, -1, 1, -1, 1, -1, 1, -1/)

! for the grid point ***************************
	do 11 k=1,2
	if(k.eq.1) then ! for the A-B upper
		  raa=ra1
		  zaa=za1
		  rbb=rb1
		  zbb=zb1
	endif
	if(k.eq.2) then ! for the C-D lower
		  raa=rc1
		  zaa=zc1
		  rbb=rd1
		  zbb=zd1
	endif

	lab=sqrt((raa-rbb)**2+(zaa-zbb)**2)

	do 1 jx=1,nbndx
	do 1 jy=iy_first,iy_last

	r0=bndx_grd(jx,1)*a0_norm
	phi0=yy(jy)
	z0=bndx_grd(jx,2)*a0_norm


	do 2 i=1,n_up

	! for coil A-B upper(k=2 is for C-D lower)

	if(k.eq.1) then ! for the A-B upper
	phi_m=(rmp_phi_up+360.d0/n_up*(i*1.d0-1.d0)-18.5d0)/180.d0*pi
	phi_p=(rmp_phi_up+360.d0/n_up*(i*1.d0-1.d0)+18.5d0)/180.d0*pi
	dphi=(phi_p-phi_m)/nphi
	endif
	if(k.eq.2) then ! for the C-D lower
	phi_m=(rmp_phi_low+360.d0/n_up*(i*1.d0-1.d0)-18.5d0)/180.d0*pi
	phi_p=(rmp_phi_low+360.d0/n_up*(i*1.d0-1.d0)+18.5d0)/180.d0*pi
	dphi=(phi_p-phi_m)/nphi
	endif

	do 3 j=1,nphi
	phi=phi_m-0.5*dphi+j*dphi
	! for A toroidal part
	ii0=rmp_set(i,k)*i0
	rx=r0*cos(phi0)-raa*cos(phi)
	ry=r0*sin(phi0)-raa*sin(phi)
	rz=z0-zaa
	rd=sqrt(rx**2+ry**2+rz**2)

	dlx=-raa*sin(phi)*dphi
	dly=raa*cos(phi)*dphi
	dlz=0.

	b_rmp_bndx(jx,jy,1)=b_rmp_bndx(jx,jy,1)+u0*ii0/4./pi*1./rd**3*(dly*rz-dlz*ry)
	b_rmp_bndx(jx,jy,1)=b_rmp_bndx(jx,jy,1)+u0*ii0/4./pi*1./rd**3*(dlz*rx-dlx*rz)
	b_rmp_bndx(jx,jy,1)=b_rmp_bndx(jx,jy,1)+u0*ii0/4./pi*1./rd**3*(dlx*ry-dly*rx)

	! for B toroidal part
	ii0=-rmp_set(i,k)*i0
	rx=r0*cos(phi0)-rbb*cos(phi)
	ry=r0*sin(phi0)-rbb*sin(phi)
	rz=z0-zbb
	rd=sqrt(rx**2+ry**2+rz**2)

	dlx=-rbb*sin(phi)*dphi
	dly=rbb*cos(phi)*dphi
	dlz=0.

	b_rmp_bndx(jx,jy,1)=b_rmp_bndx(jx,jy,1)+u0*ii0/4./pi*1./rd**3*(dly*rz-dlz*ry)
	b_rmp_bndx(jx,jy,1)=b_rmp_bndx(jx,jy,1)+u0*ii0/4./pi*1./rd**3*(dlz*rx-dlx*rz)
	b_rmp_bndx(jx,jy,1)=b_rmp_bndx(jx,jy,1)+u0*ii0/4./pi*1./rd**3*(dlx*ry-dly*rx)
    3 continue

	! for left and right poroidal part, phi=phi_m or phi_p
	dl=lab/nab

	do 4 j=1,nab
	l=0-0.5*dl+j*dl
	! for left poroidal part
	ii0=-rmp_set(i,k)*i0
	phi=phi_m
	rx=r0*cos(phi0)-(raa+(rbb-raa)/lab*l)*cos(phi)
	ry=r0*sin(phi0)-(raa+(rbb-raa)/lab*l)*sin(phi)
	rz=z0-zaa-(zbb-zaa)/lab*l
	rd=sqrt(rx**2+ry**2+rz**2)

	dlx=(rbb-raa)*cos(phi)*dl/lab
	dly=(rbb-raa)*sin(phi)*dl/lab
	dlz=(zbb-zaa)*dl/lab

	b_rmp_bndx(jx,jy,1)=b_rmp_bndx(jx,jy,1)+u0*ii0/4./pi*1./rd**3*(dly*rz-dlz*ry)
	b_rmp_bndx(jx,jy,1)=b_rmp_bndx(jx,jy,1)+u0*ii0/4./pi*1./rd**3*(dlz*rx-dlx*rz)
	b_rmp_bndx(jx,jy,1)=b_rmp_bndx(jx,jy,1)+u0*ii0/4./pi*1./rd**3*(dlx*ry-dly*rx)
	

	! for right poroidal part
	ii0=rmp_set(i,k)*i0
	phi=phi_p
	rx=r0*cos(phi0)-(raa+(rbb-raa)/lab*l)*cos(phi)
	ry=r0*sin(phi0)-(raa+(rbb-raa)/lab*l)*sin(phi)
	rz=z0-zaa-(zbb-zaa)/lab*l
	rd=sqrt(rx**2+ry**2+rz**2)

	dlx=(rbb-raa)*cos(phi)*dl/lab
	dly=(rbb-raa)*sin(phi)*dl/lab
	dlz=(zbb-zaa)*dl/lab

	b_rmp_bndx(jx,jy,1)=b_rmp_bndx(jx,jy,1)+u0*ii0/4./pi*1./rd**3*(dly*rz-dlz*ry)
	b_rmp_bndx(jx,jy,1)=b_rmp_bndx(jx,jy,1)+u0*ii0/4./pi*1./rd**3*(dlz*rx-dlx*rz)
	b_rmp_bndx(jx,jy,1)=b_rmp_bndx(jx,jy,1)+u0*ii0/4./pi*1./rd**3*(dlx*ry-dly*rx)
    4 continue

    2 continue

    1 continue

   11 continue


	do 12 jx=1,nbndx
	do 12 jy=iy_first,iy_last
	phi0=yy(jy)
	! normalization and transfer coord
	b_rmp_bndx(jx,jy,1)=b_rmp_bndx(jx,jy,1)*cos(phi0)+b_rmp_bndx(jx,jy,2)*sin(phi0)
	b_rmp_bndx(jx,jy,2)=b_rmp_bndx(jx,jy,2)*cos(phi0)-b_rmp_bndx(jx,jy,1)*sin(phi0)
	b_rmp_bndx(jx,jy,3)=b_rmp_bndx(jx,jy,3)

	b_rmp_bndx(jx,jy,:)=b_rmp_bndx(jx,jy,:)/b0_norm

   12 continue

   	x_8bndx(:,:,6:8)=x_8bndx(:,:,6:8)+b_rmp_bndx(:,:,1:3)
!	do jy=1,my
!	x1_8bndx(:,jy,:)=x_8bndx(:,jy,:)-xint_8bndx(:,:)
!	enddo

    	return
	end

!hw**************************************************************
	subroutine rmp_east_pert_bndz(u0, b0_norm, a0_norm, i0, n_up, n_low, nphi, nab,&
	ra1, za1, rb1, zb1, rc1, zc1, rd1, zd1, rmp_set)
      use declare
	implicit none
	integer i, j, n_up, n_low, nphi, nab, k
	real*8, dimension(8,2) :: rmp_set
	real*8 i0, phi_m, phi_p, raa, zaa, rbb, zbb, rcc, zcc, rdd, zdd, ii0
	real*8 ra1,za1,rb1,zb1,rc1,zc1,rd1,zd1
	real*8 dl, dlx, dly, dlz, rx, ry, rz, phi, dphi, r0, phi0, z0
	real*8 u0, b0_norm, a0_norm, lab, l, rd
      include 'mpif.h'

	allocate(b_rmp_bndz(nbndz,my,3))
	
!	u0=4.d0*pi*1.d-7
!	b0_norm=2.2745d0 ! T
!	a0_norm=0.8 ! m
!	i0=4000 ! A
!	n_up=8
!	n_low=8
!	nphi=100
!	nab=20

!	ra1=2.095
!	za1=0.765
!	rb1=2.28
!	zb1=0.586
!	rc1=2.28
!	zc1=-0.586
!	rd1=2.095
!	zd1=-0.765

!	lab=sqrt((raa-rbb)**2+(zaa-zbb)**2)

	b_rmp_bndz(:,:,:)=0.
	
!	rmp_set(:,1)=(/1, -1, 1, -1, 1, -1, 1, -1/)
!	rmp_set(:,2)=(/1, -1, 1, -1, 1, -1, 1, -1/)

! for the grid point ***************************
	do 11 k=1,2
	if(k.eq.1) then ! for the A-B upper
		  raa=ra1
		  zaa=za1
		  rbb=rb1
		  zbb=zb1
	endif
	if(k.eq.2) then ! for the C-D lower
		  raa=rc1
		  zaa=zc1
		  rbb=rd1
		  zbb=zd1
	endif

	lab=sqrt((raa-rbb)**2+(zaa-zbb)**2)

	do 1 jx=1,nbndz
	do 1 jy=iy_first,iy_last

	r0=bndz_grd(jx,1)*a0_norm
	phi0=yy(jy)
	z0=bndz_grd(jx,2)*a0_norm


	do 2 i=1,n_up

	! for coil A-B upper(k=2 is for C-D lower)

	if(k.eq.1) then ! for the A-B upper
	phi_m=(rmp_phi_up+360.d0/n_up*(i*1.d0-1.d0)-18.5d0)/180.d0*pi
	phi_p=(rmp_phi_up+360.d0/n_up*(i*1.d0-1.d0)+18.5d0)/180.d0*pi
	dphi=(phi_p-phi_m)/nphi
	endif
	if(k.eq.2) then ! for the C-D lower
	phi_m=(rmp_phi_low+360.d0/n_up*(i*1.d0-1.d0)-18.5d0)/180.d0*pi
	phi_p=(rmp_phi_low+360.d0/n_up*(i*1.d0-1.d0)+18.5d0)/180.d0*pi
	dphi=(phi_p-phi_m)/nphi
	endif

	do 3 j=1,nphi
	phi=phi_m-0.5*dphi+j*dphi
	! for A toroidal part
	ii0=rmp_set(i,k)*i0
	rx=r0*cos(phi0)-raa*cos(phi)
	ry=r0*sin(phi0)-raa*sin(phi)
	rz=z0-zaa
	rd=sqrt(rx**2+ry**2+rz**2)

	dlx=-raa*sin(phi)*dphi
	dly=raa*cos(phi)*dphi
	dlz=0.

	b_rmp_bndz(jx,jy,1)=b_rmp_bndz(jx,jy,1)+u0*ii0/4./pi*1./rd**3*(dly*rz-dlz*ry)
	b_rmp_bndz(jx,jy,1)=b_rmp_bndz(jx,jy,1)+u0*ii0/4./pi*1./rd**3*(dlz*rx-dlx*rz)
	b_rmp_bndz(jx,jy,1)=b_rmp_bndz(jx,jy,1)+u0*ii0/4./pi*1./rd**3*(dlx*ry-dly*rx)

	! for B toroidal part
	ii0=-rmp_set(i,k)*i0
	rx=r0*cos(phi0)-rbb*cos(phi)
	ry=r0*sin(phi0)-rbb*sin(phi)
	rz=z0-zbb
	rd=sqrt(rx**2+ry**2+rz**2)

	dlx=-rbb*sin(phi)*dphi
	dly=rbb*cos(phi)*dphi
	dlz=0.

	b_rmp_bndz(jx,jy,1)=b_rmp_bndz(jx,jy,1)+u0*ii0/4./pi*1./rd**3*(dly*rz-dlz*ry)
	b_rmp_bndz(jx,jy,1)=b_rmp_bndz(jx,jy,1)+u0*ii0/4./pi*1./rd**3*(dlz*rx-dlx*rz)
	b_rmp_bndz(jx,jy,1)=b_rmp_bndz(jx,jy,1)+u0*ii0/4./pi*1./rd**3*(dlx*ry-dly*rx)
    3 continue

	! for left and right poroidal part, phi=phi_m or phi_p
	dl=lab/nab

	do 4 j=1,nab
	l=0-0.5*dl+j*dl
	! for left poroidal part
	ii0=-rmp_set(i,k)*i0
	phi=phi_m
	rx=r0*cos(phi0)-(raa+(rbb-raa)/lab*l)*cos(phi)
	ry=r0*sin(phi0)-(raa+(rbb-raa)/lab*l)*sin(phi)
	rz=z0-zaa-(zbb-zaa)/lab*l
	rd=sqrt(rx**2+ry**2+rz**2)

	dlx=(rbb-raa)*cos(phi)*dl/lab
	dly=(rbb-raa)*sin(phi)*dl/lab
	dlz=(zbb-zaa)*dl/lab

	b_rmp_bndz(jx,jy,1)=b_rmp_bndz(jx,jy,1)+u0*ii0/4./pi*1./rd**3*(dly*rz-dlz*ry)
	b_rmp_bndz(jx,jy,1)=b_rmp_bndz(jx,jy,1)+u0*ii0/4./pi*1./rd**3*(dlz*rx-dlx*rz)
	b_rmp_bndz(jx,jy,1)=b_rmp_bndz(jx,jy,1)+u0*ii0/4./pi*1./rd**3*(dlx*ry-dly*rx)
	

	! for right poroidal part
	ii0=rmp_set(i,k)*i0
	phi=phi_p
	rx=r0*cos(phi0)-(raa+(rbb-raa)/lab*l)*cos(phi)
	ry=r0*sin(phi0)-(raa+(rbb-raa)/lab*l)*sin(phi)
	rz=z0-zaa-(zbb-zaa)/lab*l
	rd=sqrt(rx**2+ry**2+rz**2)

	dlx=(rbb-raa)*cos(phi)*dl/lab
	dly=(rbb-raa)*sin(phi)*dl/lab
	dlz=(zbb-zaa)*dl/lab

	b_rmp_bndz(jx,jy,1)=b_rmp_bndz(jx,jy,1)+u0*ii0/4./pi*1./rd**3*(dly*rz-dlz*ry)
	b_rmp_bndz(jx,jy,1)=b_rmp_bndz(jx,jy,1)+u0*ii0/4./pi*1./rd**3*(dlz*rx-dlx*rz)
	b_rmp_bndz(jx,jy,1)=b_rmp_bndz(jx,jy,1)+u0*ii0/4./pi*1./rd**3*(dlx*ry-dly*rx)
    4 continue

    2 continue

    1 continue

   11 continue


	do 12 jx=1,nbndz
	do 12 jy=iy_first,iy_last
	phi0=yy(jy)
	! normalization and transfer coord
	b_rmp_bndz(jx,jy,1)=b_rmp_bndz(jx,jy,1)*cos(phi0)+b_rmp_bndz(jx,jy,2)*sin(phi0)
	b_rmp_bndz(jx,jy,2)=b_rmp_bndz(jx,jy,2)*cos(phi0)-b_rmp_bndz(jx,jy,1)*sin(phi0)
	b_rmp_bndz(jx,jy,3)=b_rmp_bndz(jx,jy,3)

	b_rmp_bndz(jx,jy,:)=b_rmp_bndz(jx,jy,:)/b0_norm

   12 continue

   	x_8bndz(:,:,6:8)=x_8bndz(:,:,6:8)+b_rmp_bndz(:,:,1:3)
!	do jy=1,my
!	x1_8bndz(:,jy,:)=x_8bndz(:,jy,:)-xint_8bndz(:,:)
!	enddo

    	return
	end

