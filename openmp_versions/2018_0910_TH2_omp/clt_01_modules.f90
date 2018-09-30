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
! mxt=256,myt=64,mzt=256
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
! the modules for the CLT code, define variables and parameters
! including: 
! declare_oxpoint, 
! declare_parameter, 
! bnd_grd_set(for cut-cell, by hw),
! declare_grid, 
! declare.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module declare_oxpoint
      integer jxto,jzto,jxtx,jztx,jyo,nrkyo
      real*8 xx_o,zz_o,ps_o,rr_o,tc_o,yy_o
      real*8 xx_x,zz_x,ps_x,rr_x,tc_x,yy_x
      real*8 br_max,br_max0,br_lim
!      logical br_rise
end module declare_oxpoint

module declare_parameter
      integer, parameter :: mxt=128,myt=32,mzt=128,npsi=111,nthe=101,mbm=4*mxt+4*mzt,mr=4,mdgn=1
	integer, parameter :: npsi_pgfile=129,nlim=86,nbs=87
      integer, parameter :: npsip=npsi+1,nthe2=nthe+2,nthe3=nthe+3,n2th=2*(nthe-1),ndat=(npsi-1)*(n2th-1)
      integer, parameter :: ndat12=(npsi-1)*(nthe+4),ndat34=(npsi-1)*(nthe+4)
	integer, parameter :: ndata=mxt*mzt
      integer, parameter :: mps=npsip,ids=1,mpsa=npsi,nda=4*ids,mps4=mpsa-nda
      integer, parameter :: nprx=4,nprz=4,npry=1,nprxz=nprx*nprz,npr=nprxz*npry !,nprbm=min(npr,2*nprx+2*nprz)
      integer, parameter :: mx=mxt/nprx+4, mxm=mxt/nprx,mxn=mx*nprx
      integer, parameter :: mz=mzt/nprz+4, mzm=mzt/nprz,mzn=mz*nprz
      integer, parameter :: my=myt/npry+4, mym=myt/npry,myn=my*npry
      integer, parameter :: myx=my*mx,myz=my*mz,mxz=mx*mz,myx8=myx*8,myz8=myz*8,mxz8=mxz*8,myx3=myx*3,myz3=myz*3,mxz3=mxz*3
      integer, parameter :: mbm_nrk=4*mx+4*mz,nvmax=min(mxt,mzt)/2
      real*8, parameter :: fac=1./myt,pi=3.1415926535898
!      complex*16, parameter :: c0=cmplx(0,0),c1=cmplx(0,1.d0)

      integer nrank, nsize,nrankxz,nsizexz,nranky,nrankz,nrankx      
      integer, dimension(0:npr-1) :: nrky,nrkx,nrkz,nrkxz

end module declare_parameter


! hw**************************************************************************
module bnd_grd_set
      use declare_parameter
! two case set for the nxzs and zzero as below	
!	if(initia_from_pgfile)
!	integer, parameter :: nxzs=100*(n2th+5), n7=7, mpsa5=5, m2xt=2*mxt, m2zt=2*mzt, m2x=2*mx, m2z=2*mz
!	else
	integer, parameter :: nxzs=n2th+5, n7=7, mpsa5=5, m2xt=2*mxt, m2zt=2*mzt, m2x=2*mx, m2z=2*mz
!	endif
!	real*8, parameter :: zzero=0.d0
      integer ix_first_irpt,ix_last_irpt,iz_first_irpt,iz_last_irpt ! the rank region for  the IR points, gdtp_ep(jx,jz,1or4)==2
	! to replace the ix_first ...
      integer ix_first_irpt2,ix_last_irpt2,iz_first_irpt2,iz_last_irpt2 ! the rank region for  the IR points, gdtp_ep(jx,jz,1or4)==2
	! to replace the ix_first+2 ...
	real*8, parameter :: grd_type_ratio=0.5d0
	real*8, dimension(m2xt,n7) :: bnd_x
	real*8, dimension(m2zt,n7) :: bnd_z
	real*8, dimension(m2x,n7) :: bnd_x_ep
	real*8, dimension(m2z,n7) :: bnd_z_ep
	real*8, dimension(m2xt) :: bnd_tmpx
	real*8, dimension(m2zt) :: bnd_tmpz
	real*8, dimension(m2x) :: bnd_tmpx_ep
	real*8, dimension(m2z) :: bnd_tmpz_ep
	real*8, dimension(nxzs,mpsa5) :: xxs5, zzs5
	real*8, dimension(nxzs) :: theta, theta_tmp, xb, zb
	integer, dimension(mxt,mzt,2) :: grd_type 
!1-5 means regualr, irregular, boundary, dropped(outside). dropped(inside) for bndx in z direction or for bndz in x direction
	integer, dimension(mxt,mzt,2) :: gdtp_bndx 
!for grd_type(:,:,1), gdtp_bndx(:,:,1) -> -2, -1, 0, +1, +2, gdtp(:,:,2) -> the nearest rank of bndry point in total bndx_grd
	integer, dimension(mxt,mzt,2) :: gdtp_bndz !for grd_type(:,:,2)
!for grd_type(:,:,1), gdtp_bndz(:,:,1) -> -2, -1, 0, +1, +2, gdtp(:,:,2) -> the nearest rank of bndry point in total bndz_grd
	integer, dimension(mx,mz,6) :: gdtp_ep ! the grid type for each point in every processor
!gdtp_ep(:,:,1-3): grd_type(:,:,1) and gdtp_bndx(:,:,1-2)	
!gdtp_ep(:,:,4-6): grd_type(:,:,2) and gdtp_bndz(:,:,1-2)	
	real*8, dimension(mx,mz) :: hypb_ratio
	integer ngrdb,nbndx,nbndz,nbndx_ep,nbndz_ep
	real*8, allocatable :: bndx_grd(:,:), bndz_grd(:,:), bx_bndx(:), bx_bndz(:), bz_bndx(:), bz_bndz(:), &
				  bxdx_bndx(:), bxdx_bndz(:), bxdz_bndx(:), bxdz_bndz(:), bzdx_bndx(:), bzdx_bndz(:), &
				  bzdz_bndx(:), bzdz_bndz(:), by_bndx(:), by_bndz(:), bydx_bndx(:), bydx_bndz(:), &
				  bydz_bndx(:), bydz_bndz(:), pdx_bndx(:), pdx_bndz(:), pdz_bndx(:), pdz_bndz(:), &
				  cy_bndx(:), cy_bndz(:), cx_bndx(:), cx_bndz(:), cz_bndx(:), cz_bndz(:), &
				  uy_bndx(:), uy_bndz(:), uydx_bndx(:), uydx_bndz(:), uydz_bndx(:), uydz_bndz(:), &
				  pt_bndx(:), pt_bndz(:), ptdx_bndx(:), ptdx_bndz(:), ptdz_bndx(:), ptdz_bndz(:), &
				  rh_bndx(:), rh_bndz(:), rhdx_bndx(:), rhdx_bndz(:), rhdz_bndx(:), rhdz_bndz(:)

	real*8, allocatable :: x_8bndx(:,:,:), x1_8bndx(:,:,:), xint_8bndx(:,:), cur_3bndx(:,:,:), cint_3bndx(:,:), ef_3bndx(:,:,:), &
				  x_8bndz(:,:,:), x1_8bndz(:,:,:), xint_8bndz(:,:), cur_3bndz(:,:,:), cint_3bndz(:,:), ef_3bndz(:,:,:), &
				  updated_bndx(:), updated_bndz(:)

	real*8, allocatable :: b_rmp_bndx(:,:,:), b_rmp_bndz(:,:,:), b_rmp_bndx_tmp1(:,:,:), b_rmp_bndz_tmp1(:,:,:), &
		  b_rmp_bndx_tmp2(:,:,:), b_rmp_bndz_tmp2(:,:,:)

      real*8, allocatable :: xy_8bndx(:,:,:),xy2_8bndx(:,:,:)
      real*8, allocatable :: xy_8bndz(:,:,:),xy2_8bndz(:,:,:)

	real*8, allocatable :: bndx_grd_ep(:,:), bndz_grd_ep(:,:), bx_bndx_ep(:), bx_bndz_ep(:), bz_bndx_ep(:), bz_bndz_ep(:), &
				  bxdx_bndx_ep(:), bxdx_bndz_ep(:), bxdz_bndx_ep(:), bxdz_bndz_ep(:), bzdx_bndx_ep(:), bzdx_bndz_ep(:), &
				  bzdz_bndx_ep(:), bzdz_bndz_ep(:), by_bndx_ep(:), by_bndz_ep(:), bydx_bndx_ep(:), bydx_bndz_ep(:), &
				  bydz_bndx_ep(:), bydz_bndz_ep(:), pdx_bndx_ep(:), pdx_bndz_ep(:), pdz_bndx_ep(:), pdz_bndz_ep(:), &
				  cy_bndx_ep(:), cy_bndz_ep(:), cx_bndx_ep(:), cx_bndz_ep(:), cz_bndx_ep(:), cz_bndz_ep(:), &
				  uy_bndx_ep(:), uy_bndz_ep(:), uydx_bndx_ep(:), uydx_bndz_ep(:), uydz_bndx_ep(:), uydz_bndz_ep(:), &
				  pt_bndx_ep(:), pt_bndz_ep(:), ptdx_bndx_ep(:), ptdx_bndz_ep(:), ptdz_bndx_ep(:), ptdz_bndz_ep(:), &
				  rh_bndx_ep(:), rh_bndz_ep(:), rhdx_bndx_ep(:), rhdx_bndz_ep(:), rhdz_bndx_ep(:), rhdz_bndz_ep(:)
      real*8, dimension(mx,mz) :: az1_irx,bz1_irx,cz1_irx,dz1_irx
      real*8, dimension(mx,mz) :: ax1_irz,bx1_irz,cx1_irz,dx1_irz
      real*8, dimension(mx,mz) :: azbp_irx,bzbp_irx,czbp_irx,dzbp_irx
      real*8, dimension(mx,mz) :: azbm_irx,bzbm_irx,czbm_irx,dzbm_irx
      real*8, dimension(mx,mz) :: axbp_irz,bxbp_irz,cxbp_irz,dxbp_irz
      real*8, dimension(mx,mz) :: axbm_irz,bxbm_irz,cxbm_irz,dxbm_irz
      real*8, dimension(mx,mz) :: a2zbm_irx,b2zbm_irx,c2zbm_irx !d2fbm
      real*8, dimension(mx,mz) :: a2zbp_irx,b2zbp_irx,c2zbp_irx !d2fbp
      real*8, dimension(mx,mz) :: a2xbm_irz,b2xbm_irz,c2xbm_irz !d2fbm
      real*8, dimension(mx,mz) :: a2xbp_irz,b2xbp_irz,c2xbp_irz !d2fbp
      real*8, dimension(mx,mz) :: az2_irx,bz2_irx,cz2_irx,dz2_irx
      real*8, dimension(mx,mz) :: ax2_irz,bx2_irz,cx2_irz,dx2_irz
      real*8, dimension(mx,mz) :: dist_to_bnd
      real*8, dimension(mx,mz) :: theta_ep
	real*8, dimension(mx,mz) :: wdifcx,wdifcz
	real*8, dimension(my,8) :: wdix,wdiz
    real*8, dimension(mx,mz,6) :: type_weight

	real*8, allocatable :: axm_bndz(:), bxm_bndz(:), cxm_bndz(:), axp_bndz(:), bxp_bndz(:), cxp_bndz(:)
	real*8, allocatable :: azm_bndx(:), bzm_bndx(:), czm_bndx(:), azp_bndx(:), bzp_bndx(:), czp_bndx(:)
end module bnd_grd_set

module declare_grid
	use bnd_grd_set
      real*8, dimension(mxt) :: xxt,dxt
      real*8, dimension(mzt) :: zzt,dzt
      real*8, dimension(-1:myt+2) :: yyt,dyt
      real*8, dimension(my) :: yy,dy
      real*8, dimension(mx) :: xx,dx
      real*8, dimension(mz) :: zz,dz

      real*8, dimension(mxt,mzt) :: pst,pst_dx,pst_dz,rr2t,rrt,tht,tpt,tcht
      real*8, dimension(mx,mz)   :: psi,psi_dx,psi_dz,rr2,rr,thxz,tpxz,tcxz
      real*8, dimension(n2th+5):: thst
      real*8, dimension(n2th+5,npsi):: psst,xxst,zzst,tpst,tst,rst 

      real*8 aa,aa2,ab,psia,psiam,psia1,psmin,psmax,time
      real*8 zzero,xzero,xmin,xmax,zmin,zmax,dxx,dyy,dzz,xmg,zmg,xdim,zdim
end module declare_grid

!module declare_eq
!
!
!end module declare_eq


module declare
      use declare_grid

      integer ix_first, ix_last,iz_first, iz_last,iy_first,iy_last, jxs, jxe, jzs
      integer jze,mxa,mza,jxamin,jxamax,jzamin,jzamax,nmm,nmp,npp,npm
      integer jx,jy,jz,m, jr,jt,ka,mb,mt,irk,ma1
      complex*16 temps,tempr

      logical lrstrt,smooth,correct,hall,uniformx,uniformz,pressure,resisitive,etaj_in_e,br_rise
      logical firstmap,smoothc,smoothx1,smoothbout,smoothef,smoothpll,smoothp1ll
      logical halfx,symmetryx,symmetryz,implicitb,implicitv,divb,viscous,analysis,spectral,filt
      logical soundwave,ef_mode,rshear,bootstrap,curdriven,conductpll
	logical initia_from_pgfile, rmp_east, rho_unif
      logical rotation,constp,constu,consto,rho_from_p,eta_from_t,lrstrt_cd,invaryrho,invaryp,ohm_heat,cd_oxpoint,lpetsc 
      logical linear_mhd !using linear mhd response
      logical,dimension(8) :: lbndxfix
      logical,dimension(3) :: lbndcfix
!      integer nmode,mmode
      real*8 nmode,mmode,qmode,psmode,rrmode,xxmode,ps1mode,xx1mode,q0,qmin,qmax !,asm,bsm,csm
      real*8 gamma,dt,dtm,cfl
	real*8, dimension(2) :: ft_rmp ! ft_rmp(1) the old time weight, ft_rmp(2) the new time weight 
	real*8 omega_rmp, start_rmp, tau_rmp, time_old, finish_rmp
	real*8 rmp_phi_up, rmp_phi_low, i0_rmp
      real*8 eta0,fmu0,fmuout,pmu0,pmuout,kap0,kapout,kapll0,cdb0,cdbout,cfsm0,cfsmout
      real*8 cs_atf,fmu_atf,fmuc,pmuc,kapc,etac,etaout,etacut,etbc,etbout,csmp0,csmp0all
      real*8 pssm,pssmw,pstrans,pstransw
      real*8 alamda, alpha, alpha1, caf, caf0, caf1, epsilon, di,cext1

      real*8 wt1,ft1,gt1,ht1
      real*8 dbmax0, dbmax1, dbmin0, dbmin1
      real*8 fxmax0, fxmax1, fxmin0, fxmin1
      real*8 fzmax0, fzmax1, fzmin0, fzmin1
      real*8 fymax0, fymax1, fymin0, fymin1

      real*8 by01,bz01,xsf1,vy00,vz00,ft0,gt0,ht0,bz00,cj00,cf0,rch,uy0,oy0,p00,tm00,tm01
      real*8 cxp1,cxp2,cxm1,cxm2,wt,ft,gt,ht,gtold,g_rate,timeold
      real*8 dtx,dty,dtz, dt1, dt2, dt3, dtime, timein, width,dts
      real*8 vx, vy, vz, va2, cs2, vpx, vpy, vpz
      real*8 fmin, ferrm, ferrm1, dbmax, dbmin, ratiob,x_dbmax,x_dbmin,z_dbmax,z_dbmin
      real*8 tt1, tt2, tt
      real*8 f1,f2,f3,h1,h2,h3,g1,g2,g3,ca
      real*8 fm1,fm2,fm3,f0,fp1,fp2,fp3,a,b,c,d,xm1,x0,xp1
      real*8 alphar,prho,arho
      real*8 cip,cd00,fcd,tcds,tcde,tcdd,delcd,cb00,fbs,tbss,tbsd,ps100,psshift
      integer lbs,lcd,lrot,lpll,lscheme,lbnd,iden,idivb,lcdox
      
      integer nend,nstep,nsp,ncase,nst,nintd,np,nst1,ndstep,iter,nrcd,nstop,nstin,ndgn,neng,nper,ncycl_atfs,nploop,ncd
      integer nsmthxi,nsmthxti,nsmthxe,nsmthxte,mxrs

      
      integer, dimension(mxt) :: jzap,jzam !,it_zm,it_zp
      integer, dimension(mzt) :: jxap,jxam !,it_xm,it_xp

      integer, dimension(myt) :: jym,jyp,jym2,jyp2
      integer, dimension(nvmax) :: jxmm,jzmm,jxmp,jzmp,jxpp,jzpp,jxpm,jzpm 

      real*8, dimension(mxt) :: amxt0,bmxt0,cmxt0,zzam,zzap 
      real*8, dimension(mzt) :: amzt0,bmzt0,cmzt0,xxam,xxap
      real*8, dimension(mxt,mzt) :: amt,bmt,cmt,etat,fmut,txzt,rxzt
      real*8, dimension(mxt,mzt) :: g,gp,p,pp,qsf,omrot,omprot

      real*8, dimension(mx,mz) :: gx,cj,cj_dx,cj_dz,xbxz,tmint,tmint_dx,tmint_dz,bp0,bb0,wx2r,wz2r,wx2p,wz2p
      real*8, dimension(mx,mz,my) :: tm,bf2,bf2dx,bf2dz,bf2dy
      real*8, dimension(mx,mz,my) ::eta,etax,etaz,etay
      real*8, dimension(mx,mz) ::fmu,fmux,fmuz,cfsmb,etaint
      real*8, dimension(mx,mz) ::pmu,pmux,pmuz
      real*8, dimension(mx,mz) ::kap,kapx,kapz
      real*8, dimension(mx,mz) ::kap_pp,kap_pp_dx,kap_pp_dz
      real*8, dimension(mx,mz) ::kap_ll,kap_ll_dx,kap_ll_dz
      real*8, dimension(mx,mz) ::cdb,cdbx,cdbz
      real*8, dimension(mx,mz) ::etb,etbx,etbz
      
      real*8, dimension(mps4:mps) :: ps,asm,bsm,csm
      integer, dimension(mps4:mps) :: jxsmin,jxsmax,jzsmin,jzsmax,ip_s

      real*8, dimension(mbm) :: xxb,zzb,psb,thb,tpb,wbrx,wbrz,wbpx,wbpz
      real*8, dimension(mbm) :: asm_b,bsm_b,csm_b
      integer, dimension(mbm) :: jbx,jbz,it_b,nrankb,jxi

!      real*8, dimension(4*mbm) :: xxa1,zza1,psa1,tha1,asm_a1,bsm_a1,csm_a1
!      integer, dimension(4*mbm) :: jxa1,jza1,it_a1

      real*8, dimension(mbm,mps4:mps) :: x1s,xhs

      real*8, dimension(mx) :: ax1,bx1,cx1,dx1
      real*8, dimension(mx) :: ax2,bx2,cx2,dx2
      real*8, dimension(mx) :: axp,bxp,cxp,axm,bxm,cxm
      real*8, dimension(mx) :: axbp,bxbp,cxbp,axbm,bxbm,cxbm
      
      real*8, dimension(mz) :: az1,bz1,cz1,dz1
      real*8, dimension(mz) :: az2,bz2,cz2,dz2
      real*8, dimension(mz) :: azp,bzp,czp,azm,bzm,czm
      real*8, dimension(mz) :: azbp,bzbp,czbp,azbm,bzbm,czbm

      real*8, dimension(my) :: ay1,by1,cy1,dy1
      real*8, dimension(my) :: ay2,by2,cy2,dy2


      real*8, dimension(mxt) :: axt1,bxt1,cxt1,dxt1
      real*8, dimension(mxt) :: axt2,bxt2,cxt2,dxt2
      real*8, dimension(mxt) :: axtp,bxtp,cxtp,axtm,bxtm,cxtm
      real*8, dimension(mxt) :: axtbp,bxtbp,cxtbp,axtbm,bxtbm,cxtbm
      
      real*8, dimension(mzt) :: azt1,bzt1,czt1,dzt1
      real*8, dimension(mzt) :: azt2,bzt2,czt2,dzt2
      real*8, dimension(mzt) :: aztp,bztp,cztp,aztm,bztm,cztm
      real*8, dimension(mzt) :: aztbp,bztbp,cztbp,aztbm,bztbm,cztbm

!      real*8, dimension(myt) :: ayt1,byt1,cyt1,dyt1
!      real*8, dimension(myt) :: ayt2,byt2,cyt2,dyt2


      real*8, dimension(mx,mz,my,8) :: x,xm,xfold,xdif,x1
      real*8, dimension(mx,mz,my,8) :: xr,xy,xz,xr2,xy2,xz2,x1r,x1z
      real*8, dimension(mx,mz,my,3) :: cur,cux,cuy,cuz,xh,ef,efx,efy,efz,eta1j,perb,cub,cud !,divb_clean
      real*8, dimension(mx,mz,my) :: vr,vp,br,bp,cr,cp,divb_x,divb_y,divb_z,ps1 ,eta1 !,eta1x,eta1y,eta1z
      real*8, dimension(mx,mz,my) :: bx_xy,bx_xz,by_yx,by_yz,bz_zx,bz_zy,dvb,fcx,fcz
      real*8, dimension(mx,mz,8) :: xint,xint_dx,xint_dz
      real*8, dimension(mx,mz) ::	xint1,xint2,xint3,xint4,xint5,xint6,xint7,xint8
      real*8, dimension(mx,mz) ::	xint1_dx,xint2_dx,xint3_dx,xint4_dx,xint5_dx,xint6_dx,xint7_dx,xint8_dx
      real*8, dimension(mx,mz) ::	xint1_dz,xint2_dz,xint3_dz,xint4_dz,xint5_dz,xint6_dz,xint7_dz,xint8_dz
      real*8, dimension(mx,mz,3) :: cint,cint_dx,cint_dz,fint
      real*8, dimension(mx,mz) :: cint1,cint2,cint3,cint1_dx,cint2_dx,cint3_dx,cint1_dz,cint2_dz,cint3_dz
      real*8, dimension(mx,mz) :: cbp0,cub0,bp0c
      real*8, dimension(mx,mz,my) :: fn_cdy

      real*8, dimension(mx,mz,myt/2+1,8) :: xkc,xks
      real*8, dimension(mx,mz) :: w0
      real*8, dimension(myt) :: data0
      complex*16, dimension(myt/2+1) :: spec,spec1      
      real*8, dimension(myt/2+1) :: fmode

      real*8, dimension(ndat) :: ps_nova,xx_nova,zz_nova,bx_nova,bz_nova,bxdx_nova,bzdx_nova,bxdz_nova,bzdz_nova,th_nova
      real*8, dimension(ndat) :: pt_nova,ptdx_nova,ptdz_nova,rh_nova,rhdx_nova,rhdz_nova
      real*8, dimension(ndat) :: by_nova,bydx_nova,bydz_nova,pdx_nova,pdz_nova,cx_nova,cz_nova,cy_nova,uy_nova,uydx_nova,uydz_nova
      real*8, dimension(npsi) :: psival_nova,q_nova,qp_nova,p_nova,pp_nova,g_nova,gp_nova,f_nova
      real*8, dimension(npsi) :: fp_nova,fb_nova,fbp_nova,omrot_nova,omprot_nova
      real*8, dimension(ndat12) :: th12_nova,xx12_nova,zz12_nova
      real*8, dimension(ndat34) :: th34_nova,xx34_nova,zz34_nova

      real*8, dimension(n2th+5,npsi):: bxst,bxdxst,bxdzst,bzst,bzdxst,bzdzst,byst,bydxst,bydzst,pdxst
      real*8, dimension(n2th+5,npsi):: pdzst,cxst,czst,cyst,uyst,uydxst,uydzst

      real*8, dimension(n2th+5,mps4:mps,my,8) :: x1st
      real*8, dimension(n2th+5,mps4:mps,my,3) :: xhst
      real*8, dimension(n2th+5,mps4:mps) :: xxs,zzs,tps,wbxr,wbxt,wbzr,wbzp

      real*8, dimension(mxt,mzt) :: cx_dx,cx_dz,cz_dx,cz_dz,cy_dx,cy_dz
      real*8, dimension(mxt,mzt) :: bx,bxdx,bxdz,bz,bzdx,bzdz,by,bydx,bydz,pdx,pdz,cx,cz,cy,ux,uy,uz, &
		  uxdx, uxdz, uydx,uydz, uzdx, uzdz, bpol
      real*8, dimension(mxt,mzt) :: pt,ptdx,ptdz,rh,rhdx,rhdz
	real*8, dimension(mx,mz,my) :: bx_rmp, by_rmp, bz_rmp
!	real*8, dimension(mxt,mzt,myt) :: Ax_rmp, Ay_rmp, Az_rmp
	real*8, dimension(mx,mz,my,3) :: b_rmp
	real*8, dimension(mx,mz,my,3) :: b_rmp_out
	real*8, dimension(mx,mz,my,3) :: b_rmp_in, b_rmp_tmp1, b_rmp_tmp2
	real*8, dimension(mxt,mzt,my,3) :: At_rmp, dAx_rmp, dAy_rmp, dAz_rmp
	real*8, dimension(mxt,mzt,my,3) :: Bt_rmp

!npr
!npr_boundary
     integer, dimension(n2th+5,mps4:mps) :: nrankts
     integer  irecv,isend,mrkb,nrk
     integer, dimension(nprxz) :: nrkb,nrkb1
     integer, dimension(0:nprxz-1) :: mb_nrk,itbmin,itbmax,ix_min,ix_max,iz_min,iz_max,inrkb
     integer, dimension(mbm_nrk,0:nprxz-1) :: ib_nrk,itb_nrk,jbx_nrk,jbz_nrk
     integer, dimension(mps4:mps) :: mrks
     integer, dimension(mps4:mps,0:nprxz-1) :: nrks,mts_nrk,itsmin,itsmax
     integer, dimension(n2th+5,mps4:mps,0:nprxz-1) :: its_nrk
     integer, dimension(nprxz,mps4:mps) :: nsend
     integer, dimension(nprxz,mps4:mps,5) :: nranksend,ittransmin,ittransmax

     integer nrank_mode, nrkx_mode,nrkz_mode,jxtmode,jztmode,jxmode,jzmode    
!npr_dgn      
     integer  nrank_dgn,jdgn,mdgn_rs,mtor_dgn
     integer, dimension(n2th+5,mdgn) :: nrkdgn
     integer, dimension(n2th+5,0:npr-1,mdgn) :: itdgn_nrk
     integer, dimension(0:npr-1,mdgn) :: mtdgn_nrk,itdgnmin,itdgnmax
     integer, dimension(npr,mdgn) :: nranksend_dgn
     integer, dimension(mdgn) :: mrkdgn
     real*8, dimension(mdgn) :: qdgn,psdgn
     real*8, dimension(n2th+5,mdgn) :: xxdgn,zzdgn

!npr_groupy   
     integer group_world,mycomm_y,nsize_gy,nrank_gy,nrankiny(npry)
!     integer group_y(0:nprxz-1),comm_y(0:nprxz-1),nrank_in_gy(0:npry-1,0:nprxz-1),ngy_rank(0:nprxz-1),nysize(0:nprxz-1),nyrank(0:nprxz-1),nrankiny(npry)
!     integer group_xz(0:npry-1),comm_xz(0:npry-1),nrank_in_gxz(0:nprxz-1,0:npry-1),ngxz_rank(0:npry-1)
     integer idmgx,idmgzp,idmgzm,nrank_mgp,nrank_mgm,jxmg,jzmgp,jzmgm

     integer nline,mcycline     
     
!hall effects
     real*8 fdiout,fdi0
    real*8, dimension(mx,mz) ::fdi,fdix,fdiz

    end module declare
