	subroutine right_openacc
	integer itag
	real*8 drvx_dx,drey_dx,dez_dx,dex_dy,dez_dy,dex_dz,dey_dz
	real*8, dimension(mx,mz,my) :: rvx,rey, drvx_dx_tmp, dex_dy_tmp, dez_dy_tmp
	real*8, dimension(nbndx,my) :: rvx_bndx, rey_bndx
	real*8, dimension(nbndz,my) :: rvx_bndz, rey_bndz
	real*8 d1fc, d1f2, d1fm, d1fbp ,d1fbm, d1xf2

	!$acc data copyout(xdif) &
	!$acc create(drvx_dx_tmp,rey,rey_bndx,rey_bndz,rvx_bndz,rvx_bndx,rvx,dez_dy_tmp,dex_dy_tmp) &
	!$acc copyin(ef,bndx_grd,bndz_grd,x_8bndz,x_8bndx,ef_3bndx,ef_3bndz,xx,x,by1,ay1,cy1,dy1, &
	!$acc x1,fmu,fmux,gdtp_ep,kap,kapx,kapz,xint5_dz,xy,xy2,xz2,pmu,pmux, &
	!$acc cx1_irz,ax1_irz,bx1_irz,fmuz,ax1,bx1,cx1,dx1,pmuz,cxbm_irz, &
	!$acc axbm_irz,bxbm_irz,azbp_irx,cur,x,xr2,xx,x1r,xint1_dx,xint2_dz, &
	!$acc xint3,xint4_dz,xint5_dx,x1z,xint1_dz,xint2_dx,cint1,cint3,cint2,bzbp_irx, & 
	!$acc xint3_dz,xint3_dx,xint4,xint4_dx,xint5,dx1_irz,cxbp_irz,bxbp_irz,axbp_irz, &
	!$acc dz1_irx,czbp_irx,cz1_irx,az1_irx,bz1_irx, &
	!$acc az1,bz1,cz1,dz1,czbm_irx,azbm_irx,bzbm_irx) &
	!$acc copyin(viscous,implicitv,invaryrho,invaryp) 

	endsubroutine right_openacc



!hw*************************************************************************
      subroutine current_openacc
	integer itag
      real*8 drby_dx, rby_tmp
      real*8, dimension(mx,mz,my) :: rby
	real*8 d1f2, d1fc, d1fp, d1fbp, d1fbm, d1xf2

 	!$acc data copyin(bxbm_irz,axbm_irz,cxbm_irz,dx1,cx1,bx1,ax1,x1,az1,bx1_irz,ax1_irz,cx1_irz,xy,xx,dx1_irz,cxbp_irz,bxbp_irz,axbp_irz,bz1,cz1) &
 	!$acc copy(cur) &
 	!$acc copyin(bndz_grd,dz1,gdtp_ep,x1_8bndz) &
 	!$acc copyout(rby)

	endsubroutine current_openacc




!hw*************************************************************************
	subroutine convt_openacc
	integer itag, i
      real*8, dimension(my) :: wwy 
      real*8, dimension(mx,mz,my) :: x1r_tmp,xr2_tmp,x1z_tmp,xz2_tmp,x1_tmp
      real*8 d1f2, d1f2m, d1f2p, d1fc, d2fc, d1fp, d1fm, d1fbp, d1fbm
      real*8 d2f2, d2fbp, d2fbm, timestart1, timeend1

	!$acc data copyout(xy_8bndz) &
	!$acc copyin(xint_dx,xint_dz) &
	!$acc copy(x1z) &
	!$acc copyin(az2,bz2,cz2,dz2,az1,bz1,cz1,dz1,ax2,bx2,cx2,dx2,ax1,bx1,cx1,dx1,czbp_irx,c2zbp_irx,b2zbp_irx,a2zbp_irx,x1_8bndz,dx1_irz,cx1_irz,bx1_irz,ax1_irz) &
	!$acc copyout(xr) &
	!$acc copyin(dy2,cy2,by2,ay2) &
	!$acc copy(x1r) &
	!$acc copyin(dz1_irx,cz1_irx,bz1_irx,az1_irx) &
	!$acc copy(xz2) &
	!$acc copyout(xz) &
	!$acc copyin(dz2_irx,azbp_irx,bzbp_irx,bx2_irz,ax2_irz,cx2_irz,bxbm_irz,axbm_irz,cxbm_irz,b2xbm_irz,a2xbm_irz,c2xbm_irz,a2xbp_irz,b2xbp_irz,c2xbp_irz,bz2_irx,az2_irx,cz2_irx,bzbm_irx,azbm_irx,czbm_irx,b2zbm_irx,a2zbm_irx,c2zbm_irx) &
	!$acc copyout(xy2,xy2_8bndz) &
	!$acc copyin(x1,gdtp_ep,dx2_irz,cxbp_irz,bxbp_irz,axbp_irz) &
	!$acc copy(xr2) &
	!$acc copyin(dy1,cy1,by1,ay1) &
	!$acc copyout(xy) &
	!$acc copyin(x1_8bndx) &
	!$acc copyout(xy_8bndx,xy2_8bndx) 

	endsubroutine convt_openacc



!hw*************************************************************************
      subroutine efield_openacc
	integer itag

	!$acc data copyin(eta,cur,fdi,x,xint4,xint6,xint7,xint3) &
	!$acc copyin(ef) &
	!$acc copyout(ef) &
	!$acc copyin(x1z,xint5,cint,x1,eta1,xy,xx,xint8,x1r)


    	!$acc data copyin(hypb_ratio,gdtp_ep) &
    	!$acc copy(zz,xx) &
    	!$acc copy(bndx_grd,ef_3bndz,bndz_grd,ef,ef_3bndx)

	endsubroutine efield_openacc



!hw*************************************************************************
	subroutine bndry8_cut_cell_v2_fixed
      integer ibnd, itag

	!$acc data copyin(x_8bndx,xint_8bndz,xint_8bndx) &
	!$acc copy(x1_8bndx,x1_8bndz) &
	!$acc copyin(xint) &
	!$acc copy(x1) &
	!$acc copyin(ft_rmp,x_8bndz,x,b_rmp_bndz,b_rmp_out,b_rmp_bndx)


	!$acc data copyin(b_rmp_bndx) &
	!$acc copy(bndz_grd) &
	!$acc copyin(b_rmp_bndz,hypb_ratio) &
	!$acc copy(x_8bndz,x1_8bndz) &
	!$acc copyin(xint) &
	!$acc copy(x_8bndx) &
	!$acc copyin(xint_8bndz,xint_8bndx,gdtp_ep) &
	!$acc copy(x,x1_8bndx,bndx_grd) &
	!$acc copyin(ft_rmp,b_rmp_out) &
	!$acc copy(x1) &
	!$acc copy(zz,xx)
	endsubroutine bndry8_cut_cell_v2_fixed

!hw******************************************************************************
      subroutine stepon_openacc

	!$acc data copyin(xdif) &
	!$acc copy(x) &
	!$acc copyout(xfold,xm)

	endsubroutine stepon_openacc







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! all update device		  

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
	!$acc update device(xx)
	!$acc update device(zz)
	!$acc update device(rvx,rey, drvx_dx_tmp, dex_dy_tmp, dez_dy_tmp)
	!$acc update device(rvx_bndx,rey_bndx)
	!$acc update device(rvx_bndz,rey_bndz)
	!$acc update device(rby)
	!$acc update device(wwy)
	!$acc update device(x1r_tmp,xr2_tmp,x1z_tmp,xz2_tmp,x1_tmp)
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





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! all declare copyin		  

	!$acc declare copyin(gdtp_ep)
	!$acc declare copyin(hypb_ratio)
	!$acc declare copyin(bndx_grd,bndz_grd)
	!$acc declare copyin(x_8bndz,x_8bndx,x1_8bndz,x1_8bndx,ef_3bndx,ef_3bndz,xint8_bndz,xint_8bndx,xint_8bndz)
	!$acc declare copyin(b_rmp_bndx,b_rmp_bndz)
	!$acc declare copyin(xy_8bndx,xy2_8bndx)
	!$acc declare copyin(xy_8bndz,xy2_8bndz)
      !$acc declare copyin(az1_irx,bz1_irx,cz1_irx,dz1_irx)
      !$acc declare copyin(ax1_irz,bx1_irz,cx1_irz,dx1_irz)
      !$acc declare copyin(azbp_irx,bzbp_irx,czbp_irx,dzbp_irx)
      !$acc declare copyin(azbm_irx,bzbm_irx,czbm_irx,dzbm_irx)
     	!$acc declare copyin(axbp_irz,bxbp_irz,cxbp_irz,dxbp_irz)
      !$acc declare copyin(axbm_irz,bxbm_irz,cxbm_irz,dxbm_irz)
      !$acc declare copyin(a2zbm_irx,b2zbm_irx,c2zbm_irx)
      !$acc declare copyin(a2zbp_irx,b2zbp_irx,c2zbp_irx)
      !$acc declare copyin(a2xbm_irz,b2xbm_irz,c2xbm_irz)
      !$acc declare copyin(a2xbp_irz,b2xbp_irz,c2xbp_irz)
      !$acc declare copyin(az2_irx,bz2_irx,cz2_irx,dz2_irx)
      !$acc declare copyin(ax2_irz,bx2_irz,cx2_irz,dx2_irz)
	!$acc declare copyin(xx)
	!$acc declare copyin(zz)
	!$acc declare copyin(rvx,rey, drvx_dx_tmp, dex_dy_tmp, dez_dy_tmp)
	!$acc declare copyin(rvx_bndx,rey_bndx)
	!$acc declare copyin(rvx_bndz,rey_bndz)
	!$acc declare copyin(rby)
	!$acc declare copyin(wwy)
	!$acc declare copyin(x1r_tmp,xr2_tmp,x1z_tmp,xz2_tmp,x1_tmp)
	!$acc declare copyin(ft_rmp)
	!$acc declare copyin(eta)
	!$acc declare copyin(fmu,fmux,fmuz)
	!$acc declare copyin(pmu,pmux,pmuz)
	!$acc declare copyin(kap,kapx,kapz)
      !$acc declare copyin(ax1,bx1,cx1,dx1)
      !$acc declare copyin(ax2,bx2,cx2,dx2)
      !$acc declare copyin(az1,bz1,cz1,dz1)
      !$acc declare copyin(az2,bz2,cz2,dz2)
	!$acc declare copyin(ay1,by1,cy1,dy1)
	!$acc declare copyin(ay2,by2,cy2,dy2)
	!$acc declare copyin(x,xm,xfold,xdif,x1)
	!$acc declare copyin(xr,xy,xz,xr2,xy2,xz2,x1r,x1z)
	!$acc declare copyin(cur,ef)
	!$acc declare copyin(eta1)
	!$acc declare copyin(xint,xint_dx,xint_dz)
	!$acc declare copyin(xint3,xint4,xint5,xint6,xint7,xint8)
	!$acc declare copyin(xint1_dx,xint2_dx,xint3_dx,xint4_dx,xint5_dx)
	!$acc declare copyin(xint1_dz,xint2_dz,xint3_dz,xint4_dz,xint5_dz)
	!$acc declare copyin(cint)
	!$acc declare copyin(cint1,cint2,cint3)
	!$acc declare copyin(b_rmp_out)
	!$acc declare copyin(fdi)







!ws****************************************************************************************
      subroutine mpi_transfer8_x1_openacc
      use declare
!      real*8, dimension(mx,mz,my,8) :: w8
!      real*8, dimension(mz,my,8) :: wfx1,wfx2
!      real*8, dimension(mx,my,8) :: wfz1,wfz2
!      real*8, dimension(mx,mz,8) :: wfy1,wfy2
      include 'mpif.h'

!       
	!$acc kernels present(x1,w8,wfx1,wfx2,wfy1,wfy2,wfz1,wfz2)
	w8(:,:,:,:)=x1(:,:,:,:)
	!$acc end kernels

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


	x1(:,:,:,:)=w8(:,:,:,:)

      return
      end 
