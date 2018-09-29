!该程序是由swang最新版程序修改而来，改动方面：1.差分四阶空间精度，2边界2层待处理点(修改后），3加入hall effect （测试p=03时hall effect 旋转）
    !mmode=2
    !hall=true
    !fdiout=5.e-2
    !mxt=256,myt=32,mzt=256
     !cfl=1.5
    ! kap0=5.e-5
    !x(1-8):total of rho,p,vx,vy(v_phi),vz,bx,by(b_phi),bz
!x1(1-8):pertubation of rho,p,vx,vy(v_phi),vz,bx,by(b_phi),bz
!xint(1-8):equilibrium of rho,p,vx,vy(v_phi),vz,bx,by(b_phi),bz
!cur(1-3):pertubation of current _(x,y(phi),z)
!cint(1-3):equilibrium of current _(x,y(phi),z)
!Ef(1-3):pertubation of E-field _(x,y(phi),z)

!****coordinates***
!xx(mx),yy(my),zz(mz): coordinates(R,phi,Z) in each processes; 
!xxt(mxt),yyt(myt),zzt(mzt) for total; 
!xxst(n2th+5,npsi),zzst(n2th+5,npsi) in (theta,psi) grid; 
!xxs(n2th+5,mps4:mps),zzs(n2th+5,mps4:mps) in (theta,psi) bandary grid;

!thxz(mx,mz): theta coordinates in (R,Z) grid; tht(mxt,mzt) for total; 
!tpxz(mx,mz): R.Z.<->S(psi).P(pol). transit angle in (R,Z); tpt(mxt,mzt) for total;
!tcxz(mx,mz): tc=Ing(Jcb/R^2)dth 

!thst(n2th+5): theta coordinates in (theta,psi) grid;
!tpst(n2th+5,npsi): R.Z.<->S(psi).P(pol). transit angle in (theta,psi); tps(n2th+5,mps4:mps) for bndry;
!usually,
   !th=arc/n2th
   !tp=atan2(psdz,psdx)=atan2(bx,-bz) => cos(tp)=psdx/|grps|=-bz/bp;sin(tp)=psdz/|grps|=bx/bp;
MODULE DECLARE_OXpoint
      integer jxto,jzto,jxtx,jztx,jyo,nrkyo
      real*8 xx_O,zz_O,ps_O,rr_O,tc_O,yy_O
      real*8 xx_X,zz_X,ps_X,rr_X,tc_X,yy_X
      real*8 br_max,br_max0,br_lim
      logical br_rise
END MODULE DECLARE_OXpoint

MODULE DECLARE_parameter
      integer, parameter :: mxt=256,myt=32,mzt=256,npsi=111,nthe=101,mbm=4*mxt+4*mzt,mr=4,mdgn=1
      integer, parameter :: npsip=npsi+1,nthe2=nthe+2,nthe3=nthe+3,n2th=2*(nthe-1),ndat=(npsi-1)*(n2th-1),ndat12=(npsi-1)*(nthe+4),ndat34=(npsi-1)*(nthe+4)
      integer, parameter :: mps=npsip,ids=1,mpsa=npsi,nda=4*ids,mps4=mpsa-nda
      integer, parameter :: nprx=8,nprz=8,npry=1,nprxz=nprx*nprz,npr=nprxz*npry !,nprbm=min(npr,2*nprx+2*nprz)
      integer, parameter :: mx=mxt/nprx+4, mxm=mxt/nprx,mxn=mx*nprx
      integer, parameter :: mz=mzt/nprz+4, mzm=mzt/nprz,mzn=mz*nprz
      integer, parameter :: my=myt/npry+4, mym=myt/npry,myn=my*npry
      integer, parameter :: myx=my*mx,myz=my*mz,mxz=mx*mz,myx8=myx*8,myz8=myz*8,mxz8=mxz*8,myx3=myx*3,myz3=myz*3,mxz3=mxz*3
      integer, parameter :: mbm_nrk=4*mx+4*mz,nvmax=min(mxt,mzt)/2
      real*8, parameter :: fac=1./myt,pi=dacos(-1.d0)
      complex*16, parameter :: c0=cmplx(0,0),c1=cmplx(0,1.d0)

      integer nrank, nsize,nrankxz,nsizexz,nranky,nrankz,nrankx      
      integer, dimension(0:npr-1) :: nrky,nrkx,nrkz,nrkxz

END MODULE DECLARE_parameter

MODULE DECLARE_grid
      use DECLARE_parameter
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
      real*8 xzero,xmin,xmax,zmin,zmax,dxx,dyy,dzz,xmg,zmg
END MODULE DECLARE_grid

!MODULE DECLARE_eq
!
!
!END MODULE DECLARE_eq


MODULE DECLARE
      use DECLARE_grid

      integer ix_first, ix_last,iz_first, iz_last,iy_first,iy_last, jxs, jxe, jzs, jze,mxa,mza,jxamin,jxamax,jzamin,jzamax,nmm,nmp,npp,npm
      integer jx,jy,jz,m, jr,jt,ka,mb,mt,irk,ma1
      complex*16 temps,tempr

      logical lrstrt,smooth,correct,hall,uniformx,uniformz,pressure,resisitive,etaJ_in_E,firstmap,smoothc,smoothx1,smoothbout,smoothEf,smoothpll,smoothp1ll
      logical halfx,symmetryx,symmetryz,implicitb,implicitv,divb,viscous,analysis,spectral,filt,soundwave,Ef_mode,rshear,bootstrap,curdriven,conductpll
      logical rotation,constp,constu,consto,rho_from_p,eta_from_t,lrstrt_cd,invaryrho,invaryp,ohm_heat,cd_OXpoint,lPETSC 
      logical,dimension(8) :: lbndxfix
      logical,dimension(3) :: lbndcfix
!      integer nmode,mmode
      real*8 nmode,mmode,qmode,psmode,rrmode,xxmode,ps1mode,xx1mode,q0,qmin,qmax !,asm,bsm,csm
      real*8 gamma,dt,dtm,cfl
      real*8 eta0,fmu0,fmuout,pmu0,pmuout,kap0,kapout,kapll0,cdb0,cdbout,cfsm0,cfsmout,cs_atf,fmu_atf,fmuc,pmuc,kapc,etac,etaout,etacut,etbc,etbout,csmp0,csmp0all
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
      real*8 cIp,cd00,fcd,tcds,tcde,tcdd,delcd,cb00,fbs,tbss,tbsd,ps100,psshift
      integer lbs,lcd,lrot,lpll,lscheme,lbnd,iden,idivb,lcdox
      
      integer nend,nstep,nsp,ncase,nst,nint,np,nst1,ndstep,iter,nrcd,nstop,nstin,ndgn,neng,nper,ncycl_atfs,nploop,ncd
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

      real*8, dimension(mx,mz,my,8) :: x,xm,xfold,xdif,x1
      real*8, dimension(mx,mz,my,8) :: xr,xy,xz,xr2,xy2,xz2,x1r,x1z
      real*8, dimension(mx,mz,my,3) :: cur,cux,cuy,cuz,xh,Ef,Efx,Efy,Efz,eta1J,perb,cub,cud !,divb_clean
      real*8, dimension(mx,mz,my) :: vr,vp,br,bp,cr,cp,divb_x,divb_y,divb_z,ps1 ,eta1 !,eta1x,eta1y,eta1z
      real*8, dimension(mx,mz,my) :: bx_xy,bx_xz,by_yx,by_yz,bz_zx,bz_zy,dvb,fcx,fcz
      real*8, dimension(mx,mz,8) :: xint,xint_dx,xint_dz
      real*8, dimension(mx,mz,3) :: cint,cint_dx,cint_dz,fint
      real*8, dimension(mx,mz) :: cbp0,cub0,bp0c
      real*8, dimension(mx,mz,my) :: fn_cdy

      real*8, dimension(mx,mz,myt/2+1,8) :: xkc,xks
      real*8, dimension(mx,mz) :: w0
      real*8, dimension(myt) :: data0
      complex*16, dimension(myt/2+1) :: spec,spec1      
      real*8, dimension(myt/2+1) :: fmode

      real*8, dimension(ndat) :: ps_NOVA,xx_NOVA,zz_NOVA,bx_NOVA,bz_NOVA,bxdx_NOVA,bzdx_NOVA,bxdz_NOVA,bzdz_NOVA,th_NOVA
      real*8, dimension(ndat) :: pt_NOVA,ptdx_NOVA,ptdz_NOVA,rh_NOVA,rhdx_NOVA,rhdz_NOVA
      real*8, dimension(ndat) :: by_NOVA,bydx_NOVA,bydz_NOVA,pdx_NOVA,pdz_NOVA,cx_NOVA,cz_NOVA,cy_NOVA,uy_NOVA,uydx_NOVA,uydz_NOVA
      real*8, dimension(npsi) :: psival_NOVA,q_NOVA,qp_NOVA,p_NOVA,pp_NOVA,g_NOVA,gp_NOVA,f_NOVA,fp_NOVA,fb_NOVA,fbp_NOVA,omrot_NOVA,omprot_NOVA
      real*8, dimension(ndat12) :: th12_NOVA,xx12_NOVA,zz12_NOVA
      real*8, dimension(ndat34) :: th34_NOVA,xx34_NOVA,zz34_NOVA

      real*8, dimension(n2th+5,npsi):: bxst,bxdxst,bxdzst,bzst,bzdxst,bzdzst,byst,bydxst,bydzst,pdxst,pdzst,cxst,czst,cyst,uyst,uydxst,uydzst

      real*8, dimension(n2th+5,mps4:mps,my,8) :: x1st
      real*8, dimension(n2th+5,mps4:mps,my,3) :: xhst
      real*8, dimension(n2th+5,mps4:mps) :: xxs,zzs,tps,wbxr,wbxt,wbzr,wbzp

      real*8, dimension(mxt,mzt) :: cx_dx,cx_dz,cz_dx,cz_dz,cy_dx,cy_dz
      real*8, dimension(mxt,mzt) :: bx,bxdx,bxdz,bz,bzdx,bzdz,by,bydx,bydz,pdx,pdz,cx,cz,cy,uy,uydx,uydz,bpol
      real*8, dimension(mxt,mzt) :: pt,ptdx,ptdz,rh,rhdx,rhdz

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

    END MODULE DECLARE

    
 program TokRZ
    USE DECLARE
    USE DECLARE_OXpoint
      include 'mpif.h'
      INCLUDE 'fftw3.f'
!#include <finclude/petscsys.h>
             
!      PetscErrorCode ierr

      integer status(mpi_status_size)
!mpi   -----------------------------------------------------------------
! Initiate MPI
       call mpi_init(ierror)    
!clt-k using PETSC
!       call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

       call mpi_comm_size(mpi_comm_world,nsize,ierror)
       call mpi_comm_rank(mpi_comm_world,nrank,ierror)
       if (nsize.ne.npr) then
       print*,'the number of processors is not equal to', npr
       call mpi_abort(mpi_comm_world,1)
       stop
       endif
!mpi   ------------------------------------------------------------------
    if(nrank==0) timestart=MPI_WTIME()
! initializing
    nstep=0
    nstop=400000
    nrcd=10000
    ndgn=10
    neng=10
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



!MPI_Comm
!       call MPI_Comm_group(MPI_COMM_WORLD,group_world,ierror)
!       do kxz=0,nsizexz-1
!       do ky=0,npry-1
!       nrank_in_gy(ky,kxz)=kxz+ky*nprxz
!       enddo
!
!       call MPI_Group_incl(group_world,npry,nrank_in_gy(:,kxz),group_y(kxz),ierror)
!       call MPI_Comm_create(MPI_COMM_WORLD,group_y(kxz),comm_y(kxz),ierror)
!       call MPI_Comm_size(comm_y(kxz), nysize(kxz), ierror)
!       call MPI_Comm_rank(comm_y(kxz), nyrank(kxz), ierror)
!
!       call MPI_Group_rank(group_y(kxz),ngy_rank(kxz),ierror) 
!       enddo
       call MPI_Comm_Split(MPI_Comm_World,nrankxz,0,mycomm_y,ierror)
       call MPI_Comm_Rank(mycomm_y,nrank_gy,ierror) 
       call MPI_Comm_Size(mycomm_y,nsize_gy,ierror) 
       call MPI_Allgather(nrank,1,MPI_INTEGER,nrankiny,1,MPI_INTEGER,mycomm_y,ierror)
       write(*,*) nrank,'commy=',mycomm_y,nrankxz,nrank_gy,nrankiny(:)
!MPI_Comm


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

  !   if(lPETSC .or. lpll .ge. 4) call init_A  
!    x1=x
!    call calculate_ps1
!    call recrd_ps1
!    goto 900
      pstransw=pssmw
      pstrans=psia-pssmw !*cos(time)
      
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
        cfsmb(jx,jz)=cfsm0+0.5*cfsmout*(1+tanh((psi(jx,jz)-pstrans)/pstransw))
        etaint(jx,jz)=eta0 !*(1+(etacut-1)*tanh((tmint(jx,jz)/tm00)**(-1.5)/(etacut-1)))
        eta(jx,jz,:)=etaint(jx,jz)
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
        fmux(jx,jz)=0.5*fmuout/cosh((psi(jx,jz)-pstrans)/pstransw)**2*psi_dx(jx,jz)/pstransw-fmuc/cosh((psi(jx,jz)-psmin)/pstransw**0.5)**2*psi_dx(jx,jz)/pstransw**0.5
        fmuz(jx,jz)=0.5*fmuout/cosh((psi(jx,jz)-pstrans)/pstransw)**2*psi_dz(jx,jz)/pstransw-fmuc/cosh((psi(jx,jz)-psmin)/pstransw**0.5)**2*psi_dz(jx,jz)/pstransw**0.5

        pmux(jx,jz)=0.5*pmuout/cosh((psi(jx,jz)-pstrans)/pstransw)**2*psi_dx(jx,jz)/pstransw-pmuc/cosh((psi(jx,jz)-psmin)/pstransw**0.5)**2*psi_dx(jx,jz)/pstransw**0.5
        pmuz(jx,jz)=0.5*pmuout/cosh((psi(jx,jz)-pstrans)/pstransw)**2*psi_dz(jx,jz)/pstransw-pmuc/cosh((psi(jx,jz)-psmin)/pstransw**0.5)**2*psi_dz(jx,jz)/pstransw**0.5

        kapx(jx,jz)=0.5*kapout/cosh((psi(jx,jz)-pstrans)/pstransw)**2*psi_dx(jx,jz)/pstransw-kapc/cosh((psi(jx,jz)-psmin)/pstransw**0.5)**2*psi_dx(jx,jz)/pstransw**0.5
        kapz(jx,jz)=0.5*kapout/cosh((psi(jx,jz)-pstrans)/pstransw)**2*psi_dz(jx,jz)/pstransw-kapc/cosh((psi(jx,jz)-psmin)/pstransw**0.5)**2*psi_dz(jx,jz)/pstransw**0.5
        kap_ll_dx(jx,jz)=0
        kap_ll_dz(jx,jz)=0
        kap_pp_dx(jx,jz)=kapx(jx,jz)
        kap_pp_dz(jx,jz)=kapz(jx,jz)
!        etbx(jx,jz)=0.5*etbout/cosh((psi(jx,jz)-pstrans)/pstransw)**2*psi_dx(jx,jz)/pstransw-etbc/cosh((psi(jx,jz)-psmin)/pstransw**0.5)**2*psi_dx(jx,jz)/pstransw**0.5
!        etbz(jx,jz)=0.5*etbout/cosh((psi(jx,jz)-pstrans)/pstransw)**2*psi_dz(jx,jz)/pstransw-etbc/cosh((psi(jx,jz)-psmin)/pstransw**0.5)**2*psi_dz(jx,jz)/pstransw**0.5
        etbx(jx,jz)=-0.5*etbout/cosh((xx(jx)-xmin-0.1*(xzero-xmin))/(0.02*(xzero-xmin)))**2/(0.02*(xzero-xmin))*exp(-zz(jz)**2/(0.2*(zmax-zmin))**2)
        etbz(jx,jz)=0.5*etbout*(1-tanh((xx(jx)-xmin-0.1*(xzero-xmin))/(0.02*(xzero-xmin))))*exp(-zz(jz)**2/(0.2*(zmax-zmin))**2)*(-2*zz(jz)/(0.2*(zmax-zmin))**2)
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
    open(unit=211,file='Eymaxmin.dat',status='unknown',form='formatted')    
    open(unit=311,file='OXpoint.dat',status='unknown',form='formatted')
    else
    open(unit=11,file='dgnmax1.dat',status='unknown',form='formatted',position='append')
    open(unit=111,file='dgnmax.dat',status='unknown',form='formatted',position='append')
    open(unit=112,file='dgnmin.dat',status='unknown',form='formatted',position='append')
!    open(unit=13,file='current.dat',status='unknown',form='formatted',position='append')
!   open(unit=15,file='bfield.dat',status='unknown',form='formatted',position='append')
    open(unit=12,file='energy.dat',status='unknown',form='formatted',position='append') 
    open(unit=14,file='divb.dat',status='unknown',form='formatted',position='append') 
    open(unit=18,file='nstime.dat',status='unknown',form='formatted',position='append')
    open(unit=211,file='Eymaxmin.dat',status='unknown',form='formatted',position='append')
    open(unit=311,file='OXpoint.dat',status='unknown',form='formatted',position='append')
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
!      ps100=0.001*(psmax-psmin)
      ps1=0
      br_lim=4.02e-4
      lcdox=2
        if(.not. lrstrt_cd) then
          call find_OXpoint_1st
          call diagn_brmax0
        endif
        call distribution_cd
        if(lrstrt_cd) then 
          call readin_cud
!          call MPI_Barrier(MPI_Comm_World,ierror)
          call find_cud_OX
          call distribution_cd_OXpoint(cud(:,:,3,2))
        endif
   !_cos     
      endif

!curdriven
      if(eta_from_t) call calculate_eta
      call recrd_dssp  
    
    
100 continue

      call setdt
!mpi   ----------------------------------------------------------------
! Here collect all time step dtm and send send minimum dt to all preocess
       CALL MPI_ALLREDUCE(dt1,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                        MPI_COMM_WORLD,IERROR)
!mpi   ----------------------------------------------------------------
    if(nrank.eq.0) write(*,*) nstep,'t=',time,'dt=',dt
    if(dt.lt.1.e-5) goto 900
!     if(nstep .eq. 1)  call recrd1
!     if(nstep .eq. 10)  call recrd10
!     if(nstep .eq. 100)  call recrd100
      if(mod(nstep,nrcd).eq.0) then
      call recrd
      if(lcd .gt. 0) call recrd_cud
      if(nrank.eq.0) write(18,*) nst,ncase,nstep,time
	  nst=nst+1
      endif

      if(mod(nstep,ndgn).eq.0) then
      do jdgn=1,mdgn
      call diagn_nmmode(Ef(:,:,:,2),jdgn)
      enddo

!      call diagn
!      call diagn_max 
      call diagn_maxmin 
      call diagnatxmode               
      endif

!energy-----------------------------------------------------------------

      if(mod(nstep,neng).eq.0) then
      call energy
      
      CALL MPI_REDUCE(ft,ft1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_REDUCE(gt,gt1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_REDUCE(ht,ht1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                       MPI_COMM_WORLD,IERROR)
!mpi   -----------------------------------------------------------------
      if(nrank.eq.0) then
      if(nstep.ne.0) then
      g_rate=(gt1-gtold)/(gt1+gtold)/(time-timeold)
      gtold=gt1
      timeold=time
!      open(unit=12,file='energy.dat',status='unknown',form='formatted',position='append') 
      write(12,400) time,ft1-ft0,gt1-gt0,ht1-ht0,g_rate
400   format(5(1x,e12.5))
      endif
      if(nstep.eq.0) then
      ft0=ft1
      gt0=gt1
      ht0=ht1
      endif

      endif          
      endif
    
    if(soundwave) then
    call stepon_atfs
    else
    call stepon
    endif
!cd_OXpoint3:Ey
    if(cd_OXpoint .and. lcdox==3 .and. mod(nstep,ncd).eq.0 ) call distribution_cd_OXpoint(Ef(:,:,3,2))
    
    nstep=nstep+1    
    if(nstep.gt.nstop) goto 900
    goto 100

900   close(11)
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
      if(nrank==0) then
      timeend=MPI_WTIME()
      write(*,*) 'run time=',timeend-timestart
      endif
!clt-k using PETSC
    !   call PETSC_clean_var
     !  call PetscFinalize(ierr)

!mpi   ------------------------------------------------------------------
       call MPI_Comm_free(mycomm_y,ierror)
! End MPI
       call mpi_finalize(ierror)
!mpi   -----------------------------------------------------------------
    stop
    end
!ws*************************************************************************
    subroutine input
! --------------------------
!  This routine inputs parameters and define basic variables
!   LRSTRT: =.f., starting from t=0; =.t., continueing from
!           steps as given bz NST.
!   NEND:   the final steps intended for the current run ,
!           including the steps from the previous run.
!   NSP:    time step interval for diagnostic plots.
!   NSMTHX:  the number of rows starting from incoming
!           boundary for the smoothing region in r-direction.
!   ETA:    exact inverse of magnetic Renolds number.
!   GAMMA:  adiabatic constant.
!   BETA:   ratio of kinetic to magnetic pressure
!   NCASE:  case number of the run.
!   NST:    beginning data file number for the current run.
!   NINT:   data file number increment.
! --------------------------
!
      USE DECLARE
      include 'mpif.h'
!
      cfl=1.5
      cext1=0.5
      caf=0.75d0
      lPETSC=.false. !.false. 
      constp=.false. !.true.
      p00=0.0001
      rotation=.true.
      lrot=2
      constu=.false. !.true.
      uy0=0.00
      consto=.false. !.true.
      oy0=0.00

      lrstrt=.false.
      ncase=24
      nend=47
      nst=0
      nint=1
!!mode      
      nmode=1
      mmode=2
      qmode=1.0*mmode/nmode
      qdgn(1)=qmode
!!mode

      firstmap=.true.
      rshear=.false. !.true.
      spectral=.false. 
      filt=.true.
      smooth=.false. !.true.
      smoothc=.false. !.true.
      smoothx1=.false. !.true.
      smoothbout=.false. !.false. !
      smoothEf=.false.
      smoothpll=.false. !.true.
      smoothp1ll=.false. !.true.
      invaryrho=.false. !.true.
      invaryp=.false. !.true.

!      correct=.true.
      uniformx=.true.
      uniformz=.true.
!      halfx=.false.
      symmetryx=.false.
      symmetryz=.true.
      analysis=.false. !.true.

      resisitive=.true.
      Ef_mode=.true. !.false. 
      etaJ_in_E=.true.  
      viscous=.true.
      hall=.true.
      pressure=.false.
      soundwave=.false. !.true.
      eta_from_t=.false. !.true.
      rho_from_p=.false. !.true.

      ohm_heat =.false. 
      implicitb=.false. !.true.
      implicitv=.false. !.true.  

!!parameter
      gamma=5./3.
      
      !fmu0=1.e-5
      !fmuout=1.e-4
      !fmuc=0.e-6
      !pmu0=1.e-5
      !pmuout=1.e-4
      !pmuc=0.e-6
      !kap0=1.e-5
      !kapout=1.e-4
      !kapc=0.e-6
      !eta0=1.1e-6
      
      fmu0=1.e-6
      fmuout=1.e-4
      fmuc=0.e-6
      pmu0=1.e-6
      pmuout=1.e-4
      pmuc=0.e-6
       kap0=5.e-5
      kapout=1.e-4
      kapc=0.e-6
      eta0=1.e-5
      
      !fmu0   =2.5e-5
      !fmuout =5.e-4
      !fmuc   =0.e-6
      !pmu0   =1.e-4
      !pmuout =0.e-4
      !pmuc   =0.e-6
      !kap0   =2.5e-5
      !kapout =0.e-4
      !kapc   =0.e-6
      !kapll0 =1.e8*kap0
      !eta0   =1.e-5
      etacut =1.e2
!      eta10=1.e-2
      etbout =1.e-4
      etbc   =0

      cfsm0=1./2
      cfsmout=0 !1./24

      fdiout=5.e-2

      cdb0=0 !1.e-6 !eta0 !
      cdbout=100.0*cdb0

!!divb
      divb=.false. !.true.
      idivb=1

!!boundary --->
      lbnd=30 !0:fix; 1:free; 3:free v,fix b; 10:fix vr,br; free others; 20:fix p,vr,br, free others

      lbndxfix(1)=.false. !rho
      lbndxfix(2)=.false. !p
      lbndxfix(3)=.true.  !vr
      lbndxfix(4)=.false. !vy
      lbndxfix(5)=.false. !vp
      lbndxfix(6)=.true.  !br
      lbndxfix(7)=.false. !by
      lbndxfix(8)=.false. !bp

      lbndcfix(1)=.true.
      lbndcfix(2)=.false. 
      lbndcfix(3)=.false.
!!boundary <---

!!density
      iden=1 !iden=1: rho=(1.00000-alphar*rsq**prho)**arho
      alphar=0
      prho=1.
      arho=1.
!!density  

!!bootstrap current --->
      bootstrap=.false. !.true.
      lbs=1 !1:pertb; 10:total; 20:total(t);
      fbs=0.7
      tbss=200
      tbsd=500

!!current drive --->
      curdriven =.false.    
      lrstrt_cd =.false. !.true.   
      cd_OXpoint=.false.  
      br_rise   =.false.
      lcd=0 !0:no curdrive; 1:cd in Efield; 2:cd in current
      fcd=0.01

!!conductpll   --->  
      conductpll=.false.
      lpll=0    !1=pll_subloop   !2=pll_smthpline        !3=pll_soundwave   !4=pll_PETSC        
      lscheme=1 !1=Euler         !1=smthp_traceline      !1=Euler
                !2=RK4           !2=smthp_traceline_v1   !2=RK4 
                !3=Lax           !3=smthp_traceline_v2   !3=Lax 
                !4=Implict       !4=smthp2_traceline     !4=Implict
                !5=              !5=smthp_traceline_5p   !5=pt
                !6=              !6=smthp2_traceline_tm  
                !7=              !7=smthp_traceline_spec
      nploop=1
      csmp0=1./4
      csmp0all=1./4

      nline=(mxt+mzt)/my*(int(qmax)+1)
      mcycline=1     
      cs_atf=10.
      fmu_atf=1.e-4
      ncycl_atfs=10
!!conductpll   
    
!!grid &/|read          

      aa=1
      aa2=aa*aa
      xzero=5.
      xmin=xzero-1.1
      xmax=xzero+1.1
      zmin=-1.1
      zmax=1.1

      nsmthxi=3*mxt/4-1
      nsmthxti=4*mxt/7-1
      nsmthxe=3*mxt/4-1
      nsmthxte=4*mxt/7-1

      epsilon=1./4.
      by01=1.e-3
      bz01=-by01/epsilon
      alpha=1.
      alpha1=0.0

      if(.not.analysis) call read_nova

      call gridpnt
!!!!

      cd00=fcd*cIp
      
!      if(nrank==0) then
!      open(unit=14,file='grid.dat',status='unknown',form='formatted')
!      write(14,99)(xxt(jx),dxt(jx),jx=1,mxt),(yy(jy),dy(jy),jy=1,my)&
!     ,(zz(jz),dz(jz),jz=1,mz)
!   99 format(4(1x,f13.7))
!      close(14)
!      endif

!!filter
      if(filt) then
      do jy=1,myt/2
      if(jy.gt.myt/3) fmode(jy)=0.
      if(jy.le.myt/3) fmode(jy)=1.
      enddo
      else
      do jy=1,myt/2
      fmode(jy)=1.
      enddo
      endif

      time=0.
      nstep=0

      if(nrank==0) call print_input
!
      return
      end
!ws*************************************************************************
    subroutine print_input

      USE DECLARE
      include 'mpif.h'
!
      write(*,*) 'ncase     =',ncase
      write(*,*) 'constp    =',constp,'p00=',p00
      write(*,*) 'constu    =',constu,'uy0=',uy0
      write(*,*) 'lrstrt    =',lrstrt,'nst=',nst 
      write(*,*) 'firstmap  =',firstmap
      write(*,*) 'spectral  =',spectral
      write(*,*) 'filt      =',filt
      write(*,*) 'smooth    =',smooth
      write(*,*) 'smoothc   =',smoothc
      write(*,*) 'smoothx1  =',smoothx1
      write(*,*) 'smoothp1  =',smoothp1ll,'csmp=',csmp0
      write(*,*) 'smoothp   =',smoothpll,'csmp=',csmp0all

      write(*,*) 'symmetryx =',symmetryx
      write(*,*) 'symmetryz =',symmetryz
      write(*,*) 'analysis  =',analysis
      write(*,*) 'resisitive=',resisitive
      write(*,*) 'Ef_mode   =',Ef_mode
      write(*,*) 'etaJ_in_E =',etaJ_in_E   
      write(*,*) 'viscous   =',viscous
      write(*,*) 'hall      =',hall
      write(*,*) 'pressure  =',pressure
      write(*,*) 'soundwave =',soundwave
      write(*,*) 'bootstrap =',bootstrap,'lbs=',lbs,'fbs=',fbs 
      write(*,*) 'curdriven =',curdriven,'lcd=',lcd,'fcd=',fcd
      write(*,*) 'eta_from_t=',eta_from_t
      write(*,*) 'rho_from_p=',rho_from_p

      write(*,*) 'implicitb =',implicitb 
      write(*,*) 'implicitv =',implicitv
      write(*,*) 'divb      =',divb
      write(*,*) 'idivb     =',idivb
      write(*,*) 'lbnd      =',lbnd, '!0:fix; 1:free; 3:free v,fix b; 10:fix vr,br, free others'
      write(*,*) 'iden=',iden,'alphar=',alphar,'prho=',prho,'arho=',arho
      write(*,*) '!iden=1: rho=(1.00000-alphar*rsq**prho)**arho'
      write(*,*) '---------------------------'
      write(*,*) 'epsilon=',epsilon
      write(*,*) 'xzero=',xzero
      write(*,*) 'xmin =',xmin
      write(*,*) 'xmax =',xmax
      write(*,*) 'zmin =',zmin
      write(*,*) 'zmax =',zmax
      write(*,*) 'psmin=',psmin
      write(*,*) 'psmax=',psmax
      write(*,*) '---------------------------'
      write(*,*) 'gamma=',gamma
      write(*,*) 'eta0 =',eta0
      write(*,*) 'fmu0 =',fmu0
      write(*,*) 'pmu0 =',pmu0
      write(*,*) 'kap0 =',kap0
      write(*,*) 'cdb0 =',cdb0 
      write(*,*) 'csma =',cfsmout
      write(*,*) 'di   =',di

      write(*,*) 'cs_atf =',cs_atf
      write(*,*) 'fmu_atf=',fmu_atf
      write(*,*) 'ncl_atf=',ncycl_atfs
      write(*,*) '---------------------------'
      write(*,*) 'q0=',q0,'qmin=',qmin,'qmax=',qmax
      write(*,*) 'qmode=',qmode
      write(*,*) 'qdgn=',(qdgn(jdgn),jdgn=1,mdgn)

      open(unit=19,file='input_para.dat',status='unknown',form='formatted') 
      write(19,*) 'ncase     =',ncase
      write(19,*) 'constp    =',constp,'p00=',p00
      write(19,*) 'constu    =',constu,'uy0=',uy0
      write(19,*) 'lrstrt    =',lrstrt,'nst=',nst 
      write(19,*) 'firstmap  =',firstmap
      write(19,*) 'spectral  =',spectral
      write(19,*) 'filt      =',filt
      write(19,*) 'smooth    =',smooth
      write(19,*) 'smoothc   =',smoothc
      write(19,*) 'smoothx1  =',smoothx1
      write(19,*) 'smoothp1  =',smoothp1ll,'csmp=',csmp0
      write(19,*) 'smoothp   =',smoothpll,'csmp=',csmp0all
      write(19,*) 'symmetryx =',symmetryx
      write(19,*) 'symmetryz =',symmetryz
      write(19,*) 'analysis  =',analysis
      write(19,*) 'resisitive=',resisitive
      write(19,*) 'Ef_mode   =',Ef_mode
      write(19,*) 'etaJ_in_E =',etaJ_in_E   
      write(19,*) 'viscous   =',viscous
      write(19,*) 'hall      =',hall
      write(19,*) 'pressure  =',pressure
      write(19,*) 'soundwave =',soundwave
      write(19,*) 'bootstrap =',bootstrap,'lbs=',lbs,'fbs=',fbs 
      write(19,*) 'curdriven =',curdriven,'lcd=',lcd,'fcd=',fcd
      write(19,*) 'eta_from_t=',eta_from_t
      write(19,*) 'rho_from_p=',rho_from_p

      write(19,*) 'implicitb =',implicitb 
      write(19,*) 'implicitv =',implicitv
      write(19,*) 'divb      =',divb
      write(19,*) 'idivb     =',idivb
      write(19,*) 'lbnd      =',lbnd, '!0:fix; 1:free; 3:free v,fix b; 10:fix vr,br, free others'
      write(19,*) 'iden=',iden,'alphar=',alphar,'prho=',prho,'arho=',arho
      write(19,*) '!iden=1: rho=(1.00000-alphar*rsq**prho)**arho'
      write(19,*) '---------------------------'
      write(19,*) 'epsilon=',epsilon
      write(19,*) 'xzero=',xzero
      write(19,*) 'xmin =',xmin
      write(19,*) 'xmax =',xmax
      write(19,*) 'zmin =',zmin
      write(19,*) 'zmax =',zmax
      write(19,*) 'psmin=',psmin
      write(19,*) 'psmax=',psmax
      write(19,*) '---------------------------'
      write(19,*) 'gamma=',gamma
      write(19,*) 'eta0 =',eta0,'etac=',etac,'etaout=',etaout
      write(19,*) 'fmu0 =',fmu0,'fmuc=',fmuc,'fmuout=',fmuout
      write(19,*) 'pmu0 =',pmu0,'pmuc=',pmuc,'pmuout=',pmuout
      write(19,*) 'kap0 =',kap0,'kapc=',pmuc,'kapout=',pmuout
      write(19,*) 'cdb0 =',cdb0 
      write(19,*) 'csma =',cfsmout
      write(19,*) 'di   =',di

      write(19,*) 'cs_atf =',cs_atf
      write(19,*) 'fmu_atf=',fmu_atf
      write(19,*) 'ncl_atf=',ncycl_atfs
      write(19,*) '---------------------------'
      write(19,*) 'q0=',q0,'qmin=',qmin,'qmax=',qmax
      write(19,*) 'qmode=',qmode
      write(19,*) 'qdgn=',(qdgn(jdgn),jdgn=1,mdgn)

      close(19)
      return
      end
!ws****************************************************************************
      subroutine initia
!
!----------------
! Defines coordinates system and specifies initial configuration.
! Normlization convention:
!   1. Density --- normalised to asymtotic value, i.e., rho=1
!   2. Magnetic field --- normalised to asymtotic value, i.e.,
!                         b0=1.
!   3. Velocity --- normalised to asymtotic Alfven speed, VA=1, a
!                   natural result of 1. and 2.
!   4. Length --- normalised to a.
!   5. Time --- normalised to a/VA.
!---------------
!
     USE DECLARE
     include 'mpif.h'
     real*8 rhom,rhomp
     real*8 cIy,cIy_nrk,bp0max,bp0max1,bpc
     real*8 tm0dx,tm0dz,oy0dx,oy0dz,expo2,expo2_dx,expo2_dz
     real*8, dimension(mx,mz) :: oyint,oypint
!  d1fc= d f / dx  with fourth-order accuracy central difference
      d1fc(fm2,fm1,f0,fp1,fp2,a,b,c,d)= &
       a*(fp1-f0)+b*(f0-fm1)+c*(fp2-f0)+d*(f0-fm2)
!  d1f2= d f / dx  with second-order accuracy central difference
      d1f2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(fp1-f0) &
         -(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)
!  d1fm= d f / dx  with  one-sided difference involving -2 -1 and 0
!  points
      d1f2m(fm2,fm1,f0,xm2,xm1,x0)= &
        ( (xm2-x0)/(xm1-x0)*(fm1-f0) &
         -(xm1-x0)/(xm2-x0)*(fm2-f0) ) / (xm2-xm1)

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
!      
! 2. coeficient ax1,az1,....dx1, dz1 related to the first 
! differential with a fourth order accuracy
!
! 5. coeficient axp,axm....cxp,cxm work with left- or right-side first
! differential or boundary conditions with a third order accuracy
!
!

      do 1 jx=ix_first+2,ix_last-2
      dxp1=xx(jx+1)-xx(jx)
      dxm1=xx(jx)-xx(jx-1)
      dxp2=xx(jx+2)-xx(jx)
      dxm2=xx(jx)-xx(jx-2)
      f1= dxp1**2+dxp1*dxm2
      f2= dxp1**3-dxp1*dxm2**2
      f3= dxp1**4+dxp1*dxm2**3
      g1=-dxm1**2+dxm1*dxm2
      g2= dxm1**3-dxm1*dxm2**2
      g3=-dxm1**4+dxm1*dxm2**3
      h1= dxp2**2+dxp2*dxm2
      h2= dxp2**3-dxp2*dxm2**2
      h3= dxp2**4+dxp2*dxm2**3
      ca1= (f1*h3-f3*h1)*(g2*h3-g3*h2)-(g1*h3-g3*h1)*(f2*h3-f3*h2)
      ax2(jx)=2.*h3*(g2*h3-g3*h2)/ca1
      bx2(jx)=2.*h3*(h2*f3-h3*f2)/ca1
      cx2(jx)=2.*h3*(f2*g3-f3*g2)/ca1
      dx2(jx)=-(dxp1*ax2(jx)+dxm1*bx2(jx) &
         +dxp2*cx2(jx))/dxm2
    1 continue
    
      do 2 jx=ix_first+2,ix_last-2
      dxp1=xx(jx+1)-xx(jx)
      dxm1=xx(jx)-xx(jx-1)
      dxp2=xx(jx+2)-xx(jx)
      dxm2=xx(jx)-xx(jx-2)
      f1= dxp1   +dxp1**4/dxm2**3
      f2= dxp1**2-dxp1**4/dxm2**2
      f3= dxp1**3+dxp1**4/dxm2
      g1= dxm1   -dxm1**4/dxm2**3
      g2=-dxm1**2+dxm1**4/dxm2**2
      g3= dxm1**3-dxm1**4/dxm2
      h1= dxp2   +dxp2**4/dxm2**3
      h2= dxp2**2-dxp2**4/dxm2**2
      h3= dxp2**3+dxp2**4/dxm2
      ca1= (f1*h3-f3*h1)*(g2*h3-g3*h2)-(g1*h3-g3*h1)*(f2*h3-f3*h2)
      ax1(jx)=h3*(g2*h3-g3*h2)/ca1
      bx1(jx)=h3*(h2*f3-h3*f2)/ca1
      cx1(jx)=h3*(f2*g3-f3*g2)/ca1
      dx1(jx)=(dxp1**2*ax1(jx)-dxm1**2*bx1(jx) &
         +dxp2**2*cx1(jx))/dxm2**2
    2 continue
    

      do 3 jx=ix_first,ix_last-3
      dxp1=xx(jx+1)-xx(jx)
      dxp2=xx(jx+2)-xx(jx)
      dxp3=xx(jx+3)-xx(jx)
      f1=dxp1-dxp1**3/dxp3**2
      f2=dxp1**2-dxp1**3/dxp3
      g1=dxp2-dxp2**3/dxp3**2
      g2=dxp2**2-dxp2**3/dxp3
      ca1=f1*g2-f2*g1
      axp(jx)=g2/ca1
      bxp(jx)=-f2/ca1
      cxp(jx)=(1-axp(jx)*dxp1-bxp(jx)*dxp2)/dxp3
    3 continue
      do 4 jx=ix_first+3,ix_last
      dxm1=xx(jx)-xx(jx-1)
      dxm2=xx(jx)-xx(jx-2)
      dxm3=xx(jx)-xx(jx-3)
      f1=dxm1-dxm1**3/dxm3**2
      f2=dxm1**2-dxm1**3/dxm3
      g1=dxm2-dxm2**3/dxm3**2
      g2=dxm2**2-dxm2**3/dxm3
      ca1=f1*g2-f2*g1
      axm(jx)=g2/ca1
      bxm(jx)=-f2/ca1
      cxm(jx)=(1-axm(jx)*dxm1-bxm(jx)*dxm2)/dxm3
    4 continue
      do 5 jx=ix_first+1,ix_last-2
      dxm1=xx(jx)-xx(jx-1)
      dxp1=xx(jx+1)-xx(jx)
      dxp2=xx(jx+2)-xx(jx)
      f1=-dxm1+dxm1**3/dxp2**2
      f2=dxm1**2+dxm1**3/dxp2
      g1=dxp1-dxp1**3/dxp2**2
      g2=dxp1**2-dxp1**3/dxp2
      ca1=f1*g2-f2*g1
      axbp(jx)=g2/ca1
      bxbp(jx)=-f2/ca1
      cxbp(jx)=(1+axbp(jx)*dxm1-bxbp(jx)*dxp1)/dxp2
    5 continue
      do 6 jx=ix_first+2,ix_last-1
      dxp1=xx(jx+1)-xx(jx)
      dxm1=xx(jx)-xx(jx-1)
      dxm2=xx(jx)-xx(jx-2)
      f1=-dxp1+dxp1**3/dxm2**2
      f2=dxp1**2+dxp1**3/dxm2
      g1=dxm1-dxm1**3/dxm2**2
      g2=dxm1**2-dxm1**3/dxm2
      ca1=f1*g2-f2*g1
      axbm(jx)=g2/ca1
      bxbm(jx)=-f2/ca1
      cxbm(jx)=(1+axbm(jx)*dxp1-bxbm(jx)*dxm1)/dxm2
    6 continue

      do 21 jz=iz_first+2,iz_last-2
      dzp1=zz(jz+1)-zz(jz)
      dzm1=zz(jz)-zz(jz-1)
      dzp2=zz(jz+2)-zz(jz)
      dzm2=zz(jz)-zz(jz-2)
      f1= dzp1**2+dzp1*dzm2
      f2= dzp1**3-dzp1*dzm2**2
      f3= dzp1**4+dzp1*dzm2**3
      g1=-dzm1**2+dzm1*dzm2
      g2= dzm1**3-dzm1*dzm2**2
      g3=-dzm1**4+dzm1*dzm2**3
      h1= dzp2**2+dzp2*dzm2
      h2= dzp2**3-dzp2*dzm2**2
      h3= dzp2**4+dzp2*dzm2**3
      ca1= (f1*h3-f3*h1)*(g2*h3-g3*h2)-(g1*h3-g3*h1)*(f2*h3-f3*h2)
      az2(jz)=2.*h3*(g2*h3-g3*h2)/ca1
      bz2(jz)=2.*h3*(h2*f3-h3*f2)/ca1
      cz2(jz)=2.*h3*(f2*g3-f3*g2)/ca1
      dz2(jz)=-(dzp1*az2(jz)+dzm1*bz2(jz) &
         +dzp2*cz2(jz))/dzm2
   21 continue
    
      do 22 jz=iz_first+2,iz_last-2
      dzp1=zz(jz+1)-zz(jz)
      dzm1=zz(jz)-zz(jz-1)
      dzp2=zz(jz+2)-zz(jz)
      dzm2=zz(jz)-zz(jz-2)
      f1= dzp1   +dzp1**4/dzm2**3
      f2= dzp1**2-dzp1**4/dzm2**2
      f3= dzp1**3+dzp1**4/dzm2
      g1= dzm1   -dzm1**4/dzm2**3
      g2=-dzm1**2+dzm1**4/dzm2**2
      g3= dzm1**3-dzm1**4/dzm2
      h1= dzp2   +dzp2**4/dzm2**3
      h2= dzp2**2-dzp2**4/dzm2**2
      h3= dzp2**3+dzp2**4/dzm2
      ca1= (f1*h3-f3*h1)*(g2*h3-g3*h2)-(g1*h3-g3*h1)*(f2*h3-f3*h2)
      az1(jz)=h3*(g2*h3-g3*h2)/ca1
      bz1(jz)=h3*(h2*f3-h3*f2)/ca1
      cz1(jz)=h3*(f2*g3-f3*g2)/ca1
      dz1(jz)=(dzp1**2*az1(jz)-dzm1**2*bz1(jz) &
         +dzp2**2*cz1(jz))/dzm2**2
   22 continue
    

      do 23 jz=iz_first,iz_last-3
      dzp1=zz(jz+1)-zz(jz)
      dzp2=zz(jz+2)-zz(jz)
      dzp3=zz(jz+3)-zz(jz)
      f1=dzp1-dzp1**3/dzp3**2
      f2=dzp1**2-dzp1**3/dzp3
      g1=dzp2-dzp2**3/dzp3**2
      g2=dzp2**2-dzp2**3/dzp3
      ca1=f1*g2-f2*g1
      azp(jz)=g2/ca1
      bzp(jz)=-f2/ca1
      czp(jz)=(1-azp(jz)*dzp1-bzp(jz)*dzp2)/dzp3
   23 continue
      do 24 jz=iz_first+3,iz_last
      dzm1=zz(jz)-zz(jz-1)
      dzm2=zz(jz)-zz(jz-2)
      dzm3=zz(jz)-zz(jz-3)
      f1=dzm1-dzm1**3/dzm3**2
      f2=dzm1**2-dzm1**3/dzm3
      g1=dzm2-dzm2**3/dzm3**2
      g2=dzm2**2-dzm2**3/dzm3
      ca1=f1*g2-f2*g1
      azm(jz)=g2/ca1
      bzm(jz)=-f2/ca1
      czm(jz)=(1-azm(jz)*dzm1-bzm(jz)*dzm2)/dzm3
   24 continue
      do 25 jz=iz_first+1,iz_last-2
      dzm1=zz(jz)-zz(jz-1)
      dzp1=zz(jz+1)-zz(jz)
      dzp2=zz(jz+2)-zz(jz)
      f1=-dzm1+dzm1**3/dzp2**2
      f2=dzm1**2+dzm1**3/dzp2
      g1=dzp1-dzp1**3/dzp2**2
      g2=dzp1**2-dzp1**3/dzp2
      ca1=f1*g2-f2*g1
      azbp(jz)=g2/ca1
      bzbp(jz)=-f2/ca1
      czbp(jz)=(1+azbp(jz)*dzm1-bzbp(jz)*dzp1)/dzp2
   25 continue
      do 26 jz=iz_first+2,iz_last-1
      dzp1=zz(jz+1)-zz(jz)
      dzm1=zz(jz)-zz(jz-1)
      dzm2=zz(jz)-zz(jz-2)
      f1=-dzp1+dzp1**3/dzm2**2
      f2=dzp1**2+dzp1**3/dzm2
      g1=dzm1-dzm1**3/dzm2**2
      g2=dzm1**2-dzm1**3/dzm2
      ca1=f1*g2-f2*g1
      azbm(jz)=g2/ca1
      bzbm(jz)=-f2/ca1
      czbm(jz)=(1+azbm(jz)*dzp1-bzbm(jz)*dzm1)/dzm2
   26 continue

      do 31 jy=iy_first+2,iy_last-2
!      do 31 jy=1,my
      dyp1=yy(jy+1)-yy(jy)
      dym1=yy(jy)-yy(jy-1)
      dyp2=yy(jy+2)-yy(jy)
      dym2=yy(jy)-yy(jy-2)
      f1= dyp1**2+dyp1*dym2
      f2= dyp1**3-dyp1*dym2**2
      f3= dyp1**4+dyp1*dym2**3
      g1=-dym1**2+dym1*dym2
      g2= dym1**3-dym1*dym2**2
      g3=-dym1**4+dym1*dym2**3
      h1= dyp2**2+dyp2*dym2
      h2= dyp2**3-dyp2*dym2**2
      h3= dyp2**4+dyp2*dym2**3
      ca1= (f1*h3-f3*h1)*(g2*h3-g3*h2)-(g1*h3-g3*h1)*(f2*h3-f3*h2)
      ay2(jy)=2.*h3*(g2*h3-g3*h2)/ca1
      by2(jy)=2.*h3*(h2*f3-h3*f2)/ca1
      cy2(jy)=2.*h3*(f2*g3-f3*g2)/ca1
      dy2(jy)=-(dyp1*ay2(jy)+dym1*by2(jy) &
         +dyp2*cy2(jy))/dym2
   31 continue
    
      do 32 jy=iy_first+2,iy_last-2
!      do 32 jy=1,my
      dyp1=yy(jy+1)-yy(jy)
      dym1=yy(jy)-yy(jy-1)
      dyp2=yy(jy+2)-yy(jy)
      dym2=yy(jy)-yy(jy-2)
      f1= dyp1   +dyp1**4/dym2**3
      f2= dyp1**2-dyp1**4/dym2**2
      f3= dyp1**3+dyp1**4/dym2
      g1= dym1   -dym1**4/dym2**3
      g2=-dym1**2+dym1**4/dym2**2
      g3= dym1**3-dym1**4/dym2
      h1= dyp2   +dyp2**4/dym2**3
      h2= dyp2**2-dyp2**4/dym2**2
      h3= dyp2**3+dyp2**4/dym2
      ca1= (f1*h3-f3*h1)*(g2*h3-g3*h2)-(g1*h3-g3*h1)*(f2*h3-f3*h2)
      ay1(jy)=h3*(g2*h3-g3*h2)/ca1
      by1(jy)=h3*(h2*f3-h3*f2)/ca1
      cy1(jy)=h3*(f2*g3-f3*g2)/ca1
      dy1(jy)=(dyp1**2*ay1(jy)-dym1**2*by1(jy) &
         +dyp2**2*cy1(jy))/dym2**2
   32 continue

      amxt0(1)=0.
      bmxt0(1)=-1.e15
      cmxt0(1)=0.
      amxt0(mxt)=0.
      bmxt0(mxt)=-1.e15
      cmxt0(mxt)=0.
      amzt0(1)=0.
      bmzt0(1)=-1.e15
      cmzt0(1)=0.
      amzt0(mzt)=0.
      bmzt0(mzt)=-1.e15
      cmzt0(mzt)=0.
      do jx=2,mxt-1
      amxt0(jx)=(2.-(xxt(jx+1)-xxt(jx))/xxt(jx))/(xxt(jx)-xxt(jx-1)) &
       /(xxt(jx+1)-xxt(jx-1))
      cmxt0(jx)=(2.+(xxt(jx)-xxt(jx-1))/xxt(jx))/(xxt(jx+1)-xxt(jx)) &
       /(xxt(jx+1)-xxt(jx-1))
!      amxt0(jx)=(xxt(jx)+xxt(jx-1))/xxt(jx)/(xxt(jx)-xxt(jx-1)) &
!       /(xxt(jx+1)-xxt(jx-1))
!      cmxt0(jx)=(xxt(jx+1)+xxt(jx))/xxt(jx)/(xxt(jx+1)-xxt(jx)) &
!       /(xxt(jx+1)-xxt(jx-1))
      bmxt0(jx)=-(amxt0(jx)+cmxt0(jx))
      enddo
      do jz=2,mzt-1
      amzt0(jz)=(2.-(zzt(jz+1)-zzt(jz))/zzt(jz))/(zzt(jz)-zzt(jz-1)) &
       /(zzt(jz+1)-zzt(jz-1))
      cmzt0(jz)=(2.+(zzt(jz)-zzt(jz-1))/zzt(jz))/(zzt(jz+1)-zzt(jz)) &
       /(zzt(jz+1)-zzt(jz-1))
!      amzt0(jz)=(zzt(jz)+zzt(jz-1))/zzt(jz)/(zzt(jz)-zzt(jz-1)) &
!       /(zzt(jz+1)-zzt(jz-1))
!      cmzt0(jz)=(zzt(jz+1)+zzt(jz))/zzt(jz)/(zzt(jz+1)-zzt(jz)) &
!       /(zzt(jz+1)-zzt(jz-1))
      bmzt0(jz)=-(amzt0(jz)+cmzt0(jz))
      enddo

      call map_nova
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
        xint(jx,jz,1)=rh(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        xint(jx,jz,2)=pt(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2) !-p_NOVA(mpsa)
        xint(jx,jz,3)=0
        xint(jx,jz,4)=uy(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        xint(jx,jz,5)=0
        xint(jx,jz,6)=bx(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        xint(jx,jz,8)=bz(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        xint(jx,jz,7)=by(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)

        cint(jx,jz,1)=cx(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        cint(jx,jz,3)=cz(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        cint(jx,jz,2)=cy(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)

        xint_dx(jx,jz,1)=rhdx(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        xint_dx(jx,jz,2)=ptdx(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        xint_dx(jx,jz,3)=0
        xint_dx(jx,jz,4)=uydx(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        xint_dx(jx,jz,5)=0
        xint_dx(jx,jz,6)=bxdx(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        xint_dx(jx,jz,8)=bzdx(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        xint_dx(jx,jz,7)=bydx(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)

        cint_dx(jx,jz,1)=cx_dx(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        cint_dx(jx,jz,3)=cz_dx(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        cint_dx(jx,jz,2)=cy_dx(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)

        xint_dz(jx,jz,1)=rhdz(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        xint_dz(jx,jz,2)=ptdz(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        xint_dz(jx,jz,3)=0
        xint_dz(jx,jz,4)=uydz(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        xint_dz(jx,jz,5)=0
        xint_dz(jx,jz,6)=bxdz(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        xint_dz(jx,jz,8)=bzdz(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        xint_dz(jx,jz,7)=bydz(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)

        cint_dz(jx,jz,1)=cx_dz(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        cint_dz(jx,jz,3)=cz_dz(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        cint_dz(jx,jz,2)=cy_dz(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)

        oyint(jx,jz)=omrot(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        oypint(jx,jz)=omprot(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
      enddo
      enddo
      if(constp)  xint(:,:,2)=p00
      if(constu)  xint(:,:,4)=uy0

      if(consto) then
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
      xint(jx,jz,4)=oy0*xx(jx)
      xint_dx(jx,jz,4)=oy0
      xint_dz(jx,jz,4)=0
      enddo
      enddo
      endif

!      select case(lrot)
!      case(1)
!      do jz=iz_first,iz_last
!      do jx=ix_first,ix_last
!        xint(jx,jz,2)=xint(jx,jz,2)+xint(jx,jz,1)*xint(jx,jz,4)**2/2.d0
!        xint_dx(jx,jz,2)=xint_dx(jx,jz,2)+xint_dx(jx,jz,1)*xint(jx,jz,4)**2/2.d0+xint(jx,jz,1)*xint(jx,jz,4)*xint_dx(jx,jz,4)
!        xint_dz(jx,jz,2)=xint_dz(jx,jz,2)+xint_dz(jx,jz,1)*xint(jx,jz,4)**2/2.d0+xint(jx,jz,1)*xint(jx,jz,4)*xint_dz(jx,jz,4)
!
!!        tmint(jx,jz)=xint(jx,jz,2)/xint(jx,jz,1)
!      enddo
!      enddo
!
!
!      case(2)
!      tmint(:,:)=xint(:,:,2)/xint(:,:,1)
!      do jz=iz_first,iz_last
!      do jx=ix_first,ix_last
!      oy0dx=-oypint(jx,jz)*xint(jx,jz,8)*xx(jx)
!      oy0dz=oypint(jx,jz)*xint(jx,jz,6)*xx(jx)
!      tm0dx=xint_dx(jx,jz,2)/xint(jx,jz,1)-xint(jx,jz,2)/xint(jx,jz,1)**2*xint_dx(jx,jz,1)
!      tm0dz=xint_dz(jx,jz,2)/xint(jx,jz,1)-xint(jx,jz,2)/xint(jx,jz,1)**2*xint_dz(jx,jz,1)
!
!      expo2=exp(0.5*(xx(jx)**2-xmg**2)*oyint(jx,jz)**2/tmint(jx,jz))
!      expo2_dx=(xx(jx)**2-xmg**2)/tmint(jx,jz)*oyint(jx,jz)*oy0dx &
!                -0.5*(xx(jx)**2-xmg**2)*oyint(jx,jz)**2/tmint(jx,jz)**2*tm0dx &
!                +oyint(jx,jz)**2/tmint(jx,jz)*xx(jx)
!      expo2_dz=(xx(jx)**2-xmg**2)/tmint(jx,jz)*oyint(jx,jz)*oy0dz &
!                -0.5*(xx(jx)**2-xmg**2)*oyint(jx,jz)**2/tmint(jx,jz)**2*tm0dz
!
!      xint(jx,jz,1)=xint(jx,jz,1)*expo2
!      xint_dx(jx,jz,1)=expo2*xint_dx(jx,jz,1)+xint(jx,jz,1)*expo2_dx
!      xint_dz(jx,jz,1)=expo2*xint_dz(jx,jz,1)+xint(jx,jz,1)*expo2_dz
!
!      xint(jx,jz,2)=xint(jx,jz,2)*expo2
!      xint_dx(jx,jz,2)=expo2*xint_dx(jx,jz,2)+xint(jx,jz,2)*expo2_dx
!      xint_dz(jx,jz,2)=expo2*xint_dz(jx,jz,2)+xint(jx,jz,2)*expo2_dz
!
!      enddo
!      enddo
!
!      end select

      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
        tmint(jx,jz)=xint(jx,jz,2)/xint(jx,jz,1)
        tmint_dx(jx,jz)=xint_dx(jx,jz,2)/xint(jx,jz,1)-xint(jx,jz,2)*xint_dx(jx,jz,1)/xint(jx,jz,1)**2
        tmint_dz(jx,jz)=xint_dz(jx,jz,2)/xint(jx,jz,1)-xint(jx,jz,2)*xint_dz(jx,jz,1)/xint(jx,jz,1)**2

        cj(jx,jz)=cint(jx,jz,2)   
        bp0(jx,jz)=sqrt(xint(jx,jz,6)**2+xint(jx,jz,8)**2)
        bb0(jx,jz)=sqrt(xint(jx,jz,6)**2+xint(jx,jz,7)**2+xint(jx,jz,8)**2)
        wx2r(jx,jz)=-xint(jx,jz,8)/bp0(jx,jz)
        wz2r(jx,jz)=xint(jx,jz,6)/bp0(jx,jz)
        wx2p(jx,jz)=-xint(jx,jz,6)/bp0(jx,jz)
        wz2p(jx,jz)=-xint(jx,jz,8)/bp0(jx,jz) 

!        wx2r(jx,jz)=dcos(tpxz(jx,jz))
!        wz2r(jx,jz)=dsin(tpxz(jx,jz))
!        wx2p(jx,jz)=-dsin(tpxz(jx,jz))
!        wz2p(jx,jz)=dcos(tpxz(jx,jz))
      enddo
      enddo

!ws:avoid Bp~0: bp0->bp0c
      bp0max1=maxval(bp0)
      CALL MPI_ALLREDUCE(bp0max1,bp0max,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                        MPI_COMM_WORLD,IERROR)

      bpc=0.05*bp0max
      bp0c=bp0
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
      if(rr(jx,jz).lt.0.3 .and. bp0(jx,jz).lt.bpc) then
      bp0c(jx,jz)=bpc/2+0.5*bp0(jx,jz)**2/bpc
        wx2r(jx,jz)=-xint(jx,jz,8)/bp0c(jx,jz)
        wz2r(jx,jz)=xint(jx,jz,6)/bp0c(jx,jz)
        wx2p(jx,jz)=-xint(jx,jz,6)/bp0c(jx,jz)
        wz2p(jx,jz)=-xint(jx,jz,8)/bp0c(jx,jz) 
      endif
      enddo
      enddo
!ws:avoid Bp~0

!compute I_phi++++++++++
      cIy_nrk=0
      do jz=iz_first+2,iz_last-2
      do jx=ix_first+2,ix_last-2
      if(psi(jx,jz) .lt. psia) then 
      cIy_nrk=cIy_nrk+cj(jx,jz)*(xx(jx+1)-xx(jx-1))*(zz(jz+1)-zz(jz-1))/4
      endif
      enddo
      enddo

      CALL MPI_ALLREDUCE(cIy_nrk,cIy,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                        MPI_COMM_WORLD,IERROR)
      if(nrank==0) write(*,*) 'Ip=',cIp,'Iy=',cIy
!compute I_phi-----------

      tm01=maxval(tmint)
      CALL MPI_ALLREDUCE(tm01,tm00,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                        MPI_COMM_WORLD,IERROR)  
                          
      pstrans=0.6*psia
      pstransw=0.1*(psmax-psmin)
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
        eta(jx,jz,:)=eta0*(1+(etacut-1)*tanh(((tmint(jx,jz)/tm00)**(-1.5)-1)/(etacut-1)))
        fmu(jx,jz)=fmu0 !+0.5*fmuout*(1+tanh((psi(jx,jz)-pstrans)/pstransw))
        pmu(jx,jz)=pmu0 !+0.5*pmuout*(1+tanh((psi(jx,jz)-pstrans)/pstransw))
        kap(jx,jz)=kap0
        cdb(jx,jz)=cdb0
      enddo
      enddo

      if(bootstrap) then

      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
       p0dr=wx2r(jx,jz)*xint_dx(jx,jz,2)+wz2r(jx,jz)*xint_dz(jx,jz,2)
!      p0dr=(1.-xint(jx,jz,6)/bp0(jx,jz))*xint_dx(jx,jz,2)+(1.-xint(jx,jz,8)/bp0(jx,jz))*xint_dz(jx,jz,2)  
       cbp0(jx,jz)=-(rr(jx,jz)/xx(jx))**0.5/bp0c(jx,jz)
       cub0(jx,jz)=cbp0(jx,jz)*p0dr
!       cub0(jx,jz)=-(rr(jx,jz)/xx(jx))**0.5*p0dr/bp0(jx,jz)
!       cub0(jx,jz)=fbs*(cint(jx,jz,1)*xint(jx,jz,6)+cint(jx,jz,2)*xint(jx,jz,7)+cint(jx,jz,3)*xint(jx,jz,8))/bb0(jx,jz)
!       cbp0(jx,jz)=cub0(jx,jz)/p0dr
      enddo
      enddo
      endif
!      if(bootstrap) then
!
!      do jz=iz_first,iz_last
!      do jx=ix_first,ix_last
!       p0dr=wx2r(jx,jz)*xint_dz(jx,jz,2)+wz2r(jx,jz)*xint_dz(jx,jz,2)
!!      p0dr=(1.-xint(jx,jz,6)/bp0(jx,jz))*xint_dx(jx,jz,2)+(1.-xint(jx,jz,8)/bp0(jx,jz))*xint_dz(jx,jz,2)       
!       cub0(jx,jz)=(rr(jx,jz)/xzero)**0.5*p0dr/bp0(jx,jz)
!!       cub0(jx,jz)=fbs*(cint(jx,jz,1)*xint(jx,jz,6)+cint(jx,jz,2)*xint(jx,jz,7)+cint(jx,jz,3)*xint(jx,jz,8))/bb0(jx,jz)
!       cbp0(jx,jz)=cub0(jx,jz)/p0dr
!      enddo
!      enddo
!      endif

!      cj01=maxval(abs(cj))*xzero
!      do jz=iz_first,iz_last
!      do jx=ix_first,ix_last
!      if(psi(jx,jz) .lt. psiam)  eta(jx,jz,jy)=eta0*abs(cj01/cj(jx,jz))/xx(jx)
!      enddo
!      enddo


      do jy=iy_first,iy_last
      x(:,:,jy,:)=xint(:,:,:)
      cur(:,:,jy,:)=0
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
        br(jx,jz,jy)=x(jx,jz,jy,6)*wx2r(jx,jz)+x(jx,jz,jy,8)*wz2r(jx,jz)
        bp(jx,jz,jy)=-x(jx,jz,jy,6)*wx2p(jx,jz)+x(jx,jz,jy,8)*wz2p(jx,jz)
      enddo
      enddo

      enddo
      
      fint(:,:,:)=0
      do 12 jz=iz_first+2,iz_last-2      
      do 12 jx=ix_first+2,ix_last-2
      if(psi(jx,jz) .lt. psival_NOVA(mpsa-2)) then
!      fint(jx,jz,1)=cint(jx,jz,2)*xint(jx,jz,8)-cint(jx,jz,3)*xint(jx,jz,7) &
!       -xint_dx(jx,jz,2)+xint(jx,jz,1)*xint(jx,jz,4)**2/xx(jx)
!      fint(jx,jz,2)=cint(jx,jz,3)*xint(jx,jz,6)-cint(jx,jz,1)*xint(jx,jz,8)
!      fint(jx,jz,3)=cint(jx,jz,1)*xint(jx,jz,7)-cint(jx,jz,2)*xint(jx,jz,6) &
!       -xint_dz(jx,jz,2)

      fint(jx,jz,1)=(cint(jx,jz,2)*xint(jx,jz,8)-cint(jx,jz,3)*xint(jx,jz,7) &
       -xint_dx(jx,jz,2)+xint(jx,jz,1)*xint(jx,jz,4)**2/xx(jx)) &
       /(cint(jx,jz,2)*xint(jx,jz,8)-cint(jx,jz,3)*xint(jx,jz,7))
      fint(jx,jz,2)=cint(jx,jz,3)*xint(jx,jz,6)-cint(jx,jz,1)*xint(jx,jz,8)
      fint(jx,jz,3)=cint(jx,jz,1)*xint(jx,jz,7)-cint(jx,jz,2)*xint(jx,jz,6) &
       -xint_dz(jx,jz,2) &
       /(cint(jx,jz,1)*xint(jx,jz,7)-cint(jx,jz,2)*xint(jx,jz,6))
      w0(jx,jz)=xint_dx(jx,jz,6)+xint(jx,jz,6)/xx(jx)+xint_dz(jx,jz,8)
      endif
   12 continue

      call funmax(fint(:,:,1),fxmax0,fxmin0,x_dbmax,x_dbmin,z_dbmax,z_dbmin,xx,zz,mx,mz)   
      call funmax(fint(:,:,2),fymax0,fymin0,x_dbmax,x_dbmin,z_dbmax,z_dbmin,xx,zz,mx,mz)
      call funmax(fint(:,:,3),fzmax0,fzmin0,x_dbmax,x_dbmin,z_dbmax,z_dbmin,xx,zz,mx,mz)  
      
      call funmax(w0,dbmax0,dbmin0,x_dbmax,x_dbmin,z_dbmax,z_dbmin,xx,zz,mx,mz)
 !mpi   ----------------------------------------------------------------
      CALL MPI_REDUCE(fxmax0,fxmax1,1,MPI_DOUBLE_PRECISION,MPI_MAX,0, &
            MPI_COMM_WORLD,IERROR)
      CALL MPI_REDUCE(fxmin0,fxmin1,1,MPI_DOUBLE_PRECISION,MPI_MIN,0, &
            MPI_COMM_WORLD,IERROR)

      CALL MPI_REDUCE(fzmax0,fzmax1,1,MPI_DOUBLE_PRECISION,MPI_MAX,0, &
            MPI_COMM_WORLD,IERROR)
      CALL MPI_REDUCE(fzmin0,fzmin1,1,MPI_DOUBLE_PRECISION,MPI_MIN,0, &
            MPI_COMM_WORLD,IERROR)

      CALL MPI_REDUCE(fymax0,fymax1,1,MPI_DOUBLE_PRECISION,MPI_MAX,0, &
            MPI_COMM_WORLD,IERROR)
      CALL MPI_REDUCE(fymin0,fymin1,1,MPI_DOUBLE_PRECISION,MPI_MIN,0, &
            MPI_COMM_WORLD,IERROR)

      CALL MPI_REDUCE(dbmax0,dbmax1,1,MPI_DOUBLE_PRECISION,MPI_MAX,0, &
            MPI_COMM_WORLD,IERROR)
      CALL MPI_REDUCE(dbmin0,dbmin1,1,MPI_DOUBLE_PRECISION,MPI_MIN,0, &
            MPI_COMM_WORLD,IERROR)
!mpi   ----------------------------------------------------------------
      if(nrank==0) then
      write(*,*) 'dbmax=',dbmax1,'dbmin=',dbmin1
      write(*,*) 'fxmax=',fxmax1,'fxmin=',fxmin1
      write(*,*) 'fzmax=',fzmax1,'fzmin=',fzmin1
      write(*,*) 'fymax=',fymax1,'fymin=',fymin1
      endif

      call recrd_init
!      call recrd_dssp

      ht0=0.
      gt0=0.
      ft0=0.
      call energy
!mpi   -----------------------------------------------------------------      
      CALL MPI_REDUCE(ft,ft1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_REDUCE(gt,gt1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_REDUCE(ht,ht1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                       MPI_COMM_WORLD,IERROR)
!mpi   -----------------------------------------------------------------

      if(nrank.eq.0) then
      open(unit=17,file='energy_init.dat',status='unknown',form='formatted')
      write(17,*) "Magnetic,Kinetic,Heat,Total" 
      write(17,*) "eqm:"
      write(17,400) ft1,gt1,ht1,ft1+gt1+ht1
400   format(4(1x,e12.5))
      endif
      ft0=ft1
      gt0=gt1
      ht0=ht1

      call perturbation

      return
      end
!ws**************************************************************************
      subroutine perturbation
      USE DECLARE
!      real*8, dimension(mx,mz,my,3) :: eta1J
      real*8 dREy_dx,dEy_dx,dEz_dx,dEx_dy,dEz_dy,dEx_dz,dEy_dz,per_vp !,eta1
      include 'mpif.h'
!  d1fc= d f / dx  with fourth-order accuracy central difference
      d1fc(fm2,fm1,f0,fp1,fp2,a,b,c,d)= &
       a*(fp1-f0)+b*(f0-fm1)+c*(fp2-f0)+d*(f0-fm2)
!  d1f2= d f / dx  with second-order accuracy central difference
      d1f2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(fp1-f0) &
         -(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)
!  d1xf2= d Rf / dx  with second-order accuracy central difference
      d1xf2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(xp1*fp1-x0*f0) &
         -(xp1-x0)/(xm1-x0)*(xm1*fm1-x0*f0))/(xm1-xp1)

      deltax=0.03*(psmax-psmin)
      psi0=1.e-6
      eta10=eta0

      do 15 jy=iy_first,iy_last
      do 15 jz=iz_first,iz_last
      do 15 jx=ix_first,ix_last
      if(psi(jx,jz) .lt. psiam) then
!!      vr(jx,jz,jy)=psi0*(1.-rr(jx,jz))*dtanh((rr(jx,jz)-rrmode)/deltax) &
!!       *dcos(yy(jy)-qmode*th(jx,jz))
!!      vp(jx,jz,jy)=-psi0*(rr(jx,jz)*(1.-rr(jx,jz))/deltax &
!!       /dcosh((rr(jx,jz)-rrmode)/deltax)/dcosh((rr(jx,jz)-rrmode)/deltax) &
!!       +(1.-2.*rr(jx,jz))*dtanh((rr(jx,jz)-rrmode)/deltax)) &
!!        *dsin(yy(jy)-qmode*th(jx,jz))
!!      br(jx,jz,jy)=psi0*(exp(-(rr(jx,jz)-rrmode)**2/deltax**2)-exp(-(aa-rrmode)**2/deltax**2)) &
!!       *dcos(yy(jy)-qmode*th(jx,jz))
!!      bp(jx,jz,jy)=0.
!!
!!      x(jx,jz,jy,6)=xint(jx,jz,6)+br(jx,jz,jy)*dcos(th(jx,jz))-bp(jx,jz,jy)*dsin(th(jx,jz))
!!      x(jx,jz,jy,8)=xint(jx,jz,8)+br(jx,jz,jy)*dsin(th(jx,jz))+bp(jx,jz,jy)*dcos(th(jx,jz))
!
!
!      x(jx,jz,jy,4)=xint(jx,jz,4)+psi0*(exp(-(psi(jx,jz)-psmode)**2/deltax**2)-exp(-(psiam-psmode)**2/deltax**2)) &
!!      *dsin(thxz(jx,jz))
!       *dcos(yy(jy)-qmode*thxz(jx,jz))
!!        *dcos(qmode*th(jx,jz))
!!      x(jx,jz,jy,3)=sin(pi*(xx(jx)-xx_xm(jz,mpsa))/(xx_xp(jz,mpsa)-xx_xm(jz,mpsa)))
!      endif
!      ps1(jx,jz,jy)=psi0*(exp(-(psi(jx,jz)-psmode)**2/deltax**2))*dcos(yy(jy)-qmode*thxz(jx,jz))
!       per_vp=-psi0*(exp(-(psi(jx,jz)-psmode)**2/deltax**2))*dsin(nmode*yy(jy)+mmode*thxz(jx,jz))
!       x(jx,jz,jy,3)=-per_vp*dsin(thxz(jx,jz))
!       x(jx,jz,jy,5)=per_vp*dcos(thxz(jx,jz))
      eta1(jx,jz,jy)=eta10*(exp(-(psi(jx,jz)-psmode)**2/deltax**2))*dcos(nmode*yy(jy)+mmode*thxz(jx,jz))
!      if(rshear) then
!      eta1(jx,jz,jy)=eta10*(exp(-(psi(jx,jz)-psmode)**2/deltax**2)+exp(-(psi(jx,jz)-ps1mode)**2/deltax**2))*((1+dcos(nmode*yy(jy)+mmode*thxz(jx,jz)))/2)**4
!      else
!      eta1(jx,jz,jy)=eta10*(exp(-(psi(jx,jz)-psmode)**2/deltax**2))*((1+dcos(nmode*yy(jy)+mmode*thxz(jx,jz)))/2)**4
!      endif
      do m=1,3
      eta1J(jx,jz,jy,m)=eta1(jx,jz,jy)*cint(jx,jz,m)
      enddo
      endif
   15 continue

      do 16 jy=iy_first+2,iy_last-2
      do 16 jz=iz_first+2,iz_last-2
      do 16 jx=ix_first+2,ix_last-2
      if(psi(jx,jz) .le. psiam) then
!      eta1y(jx,jz,jy)=-eta10*(exp(-(psi(jx,jz)-psmode)**2/deltax**2))*dsin(yy(jy)-qmode*thxz(jx,jz))
!      eta1x(jx,jz,jy)=d1fc(eta1(jx-2,jz,jy),eta1(jx-1,jz,jy),eta1(jx,jz,jy) &
!         ,eta1(jx+1,jz,jy),eta1(jx+2,jz,jy),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
!      eta1z(jx,jz,jy)=d1fc(eta1(jx,jz-2,jy),eta1(jx,jz-1,jy),eta1(jx,jz,jy) &
!         ,eta1(jx,jz+1,jy),eta1(jx,jz+2,jy),az1(jz),bz1(jz),cz1(jz),dz1(jz))

      dREy_dx=xx(jx)*d1fc(eta1J(jx-2,jz,jy,2),eta1J(jx-1,jz,jy,2),eta1J(jx,jz,jy,2),eta1J(jx+1,jz,jy,2),eta1J(jx+2,jz,jy,2),ax1(jx),bx1(jx),cx1(jx),dx1(jx))+eta1J(jx,jz,jy,2)
!      dEy_dx =d1f2(eta1J(jx-1,jz,jy,2),eta1J(jx,jz,jy,2),eta1J(jx+1,jz,jy,2),xx(jx-1),xx(jx),xx(jx+1))
      dEz_dx =d1fc(eta1J(jx-2,jz,jy,3),eta1J(jx-1,jz,jy,3),eta1J(jx,jz,jy,3),eta1J(jx+1,jz,jy,3),eta1J(jx+2,jz,jy,3),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
      dEx_dz =d1fc(eta1J(jx,jz-2,jy,1),eta1J(jx,jz-1,jy,1),eta1J(jx,jz,jy,1),eta1J(jx,jz+1,jy,1),eta1J(jx,jz+2,jy,1),az1(jz),bz1(jz),cz1(jz),dz1(jz))
      dEy_dz =d1fc(eta1J(jx,jz-2,jy,2),eta1J(jx,jz-1,jy,2),eta1J(jx,jz,jy,2),eta1J(jx,jz+1,jy,2),eta1J(jx,jz+2,jy,2),az1(jz),bz1(jz),cz1(jz),dz1(jz))
      dEx_dy =d1fc(eta1J(jx,jz,jy-2,1),eta1J(jx,jz,jy-1,1),eta1J(jx,jz,jy,1),eta1J(jx,jz,jy+1,1),eta1J(jx,jz,jy+2,1),ay1(jy),by1(jy),cy1(jy),dy1(jy))
      dEz_dy =d1fc(eta1J(jx,jz,jy-2,1),eta1J(jx,jz,jy-1,3),eta1J(jx,jz,jy,3),eta1J(jx,jz,jy+1,3),eta1J(jx,jz,jy+2,1),ay1(jy),by1(jy),cy1(jy),dy1(jy))

!      dREy_dx=d1xf2(eta1J(jx-1,jz,jy,2),eta1J(jx,jz,jy,2),eta1J(jx+1,jz,jy,2),xx(jx-1),xx(jx),xx(jx+1))
!!      dEy_dx =d1f2(eta1J(jx-1,jz,jy,2),eta1J(jx,jz,jy,2),eta1J(jx+1,jz,jy,2),xx(jx-1),xx(jx),xx(jx+1))
!      dEz_dx =d1f2(eta1J(jx-1,jz,jy,3),eta1J(jx,jz,jy,3),eta1J(jx+1,jz,jy,3),xx(jx-1),xx(jx),xx(jx+1))
!      dEx_dz =d1f2(eta1J(jx,jz-1,jy,1),eta1J(jx,jz,jy,1),eta1J(jx,jz+1,jy,1),zz(jz-1),zz(jz),zz(jz+1))
!      dEy_dz =d1f2(eta1J(jx,jz-1,jy,2),eta1J(jx,jz,jy,2),eta1J(jx,jz+1,jy,2),zz(jz-1),zz(jz),zz(jz+1))
!      dEx_dy =d1f2(eta1J(jx,jz,jy-1,1),eta1J(jx,jz,jy,1),eta1J(jx,jz,jy+1,1),yy(jy-1),yy(jy),yy(jy+1))
!      dEz_dy =d1f2(eta1J(jx,jz,jy-1,3),eta1J(jx,jz,jy,3),eta1J(jx,jz,jy+1,3),yy(jy-1),yy(jy),yy(jy+1))
          
      perb(jx,jz,jy,1)=-dEz_dy/xx(jx)+dEy_dz  
      perb(jx,jz,jy,2)=-dEx_dz+dEz_dx
!      perb(jx,jz,jy,3)=(dEx_dy-eta1J(jx,jz,jy,2))/xx(jx)-dEy_dx
      perb(jx,jz,jy,3)=(dEx_dy-dREy_dx)/xx(jx)

!      x1(jx,jz,jy,6)=d1fc(ps1(jx,jz-2,jy),ps1(jx,jz-1,jy),ps1(jx,jz,jy) &
!         ,ps1(jx,jz+1,jy),ps1(jx,jz+2,jy),az1(jz),bz1(jz),cz1(jz),dz1(jz))/xx(jx)
!
!      x1(jx,jz,jy,8)=-d1fc(ps1(jx-2,jz,jy),ps1(jx-1,jz,jy),ps1(jx,jz,jy) &
!         ,ps1(jx+1,jz,jy),ps1(jx+2,jz,jy),ax1(jx),bx1(jx),cx1(jx),dx1(jx))/xx(jx)
!      x(jx,jz,jy,6)=xint(jx,jz,6)+x1(jx,jz,jy,6)
!      x(jx,jz,jy,8)=xint(jx,jz,8)+x1(jx,jz,jy,8) 
      endif
   16 continue
!      open(unit=116,file='eta1'//cn1(nrank),status='unknown',form='formatted')
!      write(116,400)(((eta1(jx,jz,jy),eta1x(jx,jz,jy),eta1z(jx,jz,jy),eta1y(jx,jz,jy),jx=ix_first,ix_last),jz=iz_first,iz_last),jy=1,my)
!      close(116)

      call energy
!mpi   -----------------------------------------------------------------      
      CALL MPI_REDUCE(ft,ft1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_REDUCE(gt,gt1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_REDUCE(ht,ht1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                       MPI_COMM_WORLD,IERROR)
!mpi   -----------------------------------------------------------------
      if(nrank.eq.0) then
      timeold=time
      gtold=gt1
      open(unit=17,file='energy_init.dat',status='unknown',form='formatted') 
      write(17,*) "perb:"
      write(17,400) ft1-ft0,gt1-gt0,ht1-ht0,ft1+gt1+ht1-ft0-gt0-ht0
400   format(4(1x,e12.5))
      close(17)
      endif
      return
      end
!ws******************************************************************************
      subroutine stepon
!
!     This routine time-advances X's bz fourth order in time and second
!     order in space Runge-Kotta differential scheme.
!     note: X is alwazs the up-to-date value while Xm being the
!           intermediate value, and Xdif is increment
!
!
      USE DECLARE
      include 'mpif.h'

!      dts=dt/ncycl_atfs
      ms=1
      me=8
      ml=1
      tt=time
      tt1=time+dt/6.
      tt2=time
      irk=1
      call right
       xfold(:,:,:,:)=x(:,:,:,:)
       xm(:,:,:,:)=xfold(:,:,:,:)+xdif(:,:,:,:)*dt/6.
       x(:,:,:,:)=xfold(:,:,:,:)+xdif(:,:,:,:)*dt/2.
!
      tt=time+dt/2.
      tt1=time+dt/2.
      tt2=time+dt/6.
      irk=2
        call right
        xm(:,:,:,:)=xm(:,:,:,:)+xdif(:,:,:,:)*dt/3.
        x(:,:,:,:)=xfold(:,:,:,:)+xdif(:,:,:,:)*dt/2.
!
      tt1=time+5.*dt/6.
      tt2=time+dt/2.
      irk=3
        call right
        xm(:,:,:,:)=xm(:,:,:,:)+xdif(:,:,:,:)*dt/3.
        x(:,:,:,:)=xfold(:,:,:,:)+xdif(:,:,:,:)*dt
!
      time=time+dt
      tt1=time+dt
      tt2=time+5.*dt/6.
      irk=4
        call right
        x(:,:,:,:)=xm(:,:,:,:)+xdif(:,:,:,:)*dt/6.

      caf=0.75d0*(0.5+0.5*dtanh((time-40)/5.))
!      call bndry_x_ex(lbnd)
      call bndry8_x_ex(lbnd)
  !    if(conductpll) call pllconduct(lpll)
!      if(smoothpll) call smthp_traceline_5p(1)
      if(eta_from_t) call calculate_eta
      
      return
      end

!ws***************************************************************
     subroutine right
      USE DECLARE
      include 'mpif.h'
      real*8 dRvx_dx,dREy_dx,dEz_dx,dEx_dy,dEz_dy,dEx_dz,dEy_dz
      real*8, dimension(mx,mz,my) :: Rvx,REy
!
!
!  define statement functions
!  d2fc= d2 f / dx2 with third-order accuracy central difference
!      d2fc(fm2,fm1,f0,fp1,fp2,a,b,c,d)= &
!       a*(fp1-f0)+b*(f0-fm1)+c*(fp2-f0)+d*(f0-fm2)
!  d1fc= d f / dx  with fourth-order accuracy central difference
      d1fc(fm2,fm1,f0,fp1,fp2,a,b,c,d)= &
       a*(fp1-f0)+b*(f0-fm1)+c*(fp2-f0)+d*(f0-fm2)
!  d1f2= d f / dx  with second-order accuracy central difference
      d1f2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(fp1-f0) &
         -(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)
!  d1fm= d f / dx  with  one-sided difference involving -2 -1 and 0
!  points
      d1fm(fm2,fm1,f0,xm2,xm1,x0)= &
        ( (xm2-x0)/(xm1-x0)*(fm1-f0) &
         -(xm1-x0)/(xm2-x0)*(fm2-f0) ) / (xm2-xm1)
!  d1fbp= d f / dx  with  one-sided-bias  difference involving -1 0  1 and 2
!  points
      d1fbp(fp2,fp1,f0,fm1,a,b,c)= &
       a*(fm1-f0)+b*(fp1-f0)+c*(fp2-f0)
!  d1fbm= d f / dx  with  one-sided-bias  difference involving -2 -1  0 and 1
!  points
      d1fbm(fm2,fm1,f0,fp1,a,b,c)= &
       a*(f0-fp1)+b*(f0-fm1)+c*(f0-fm2)
!  d1xf2= d Rf / dx  with second-order accuracy central difference
      d1xf2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(xp1*fp1-x0*f0) &
         -(xp1-x0)/(xm1-x0)*(xm1*fm1-x0*f0))/(xm1-xp1)
!
!       integer status(mpi_status_size)

!      do 20 jy=1,my
!      x1(:,:,jy,:)=x(:,:,jy,:)-xint(:,:,:)
!   20 continue  
!      call bndry_x1_ex(lbnd)
!      call bndry8_ex(x1,lbnd)
      call bndry8_x_ex(lbnd)
!      if(smoothp1ll) call smthp1_traceline(3)
      if(rho_from_p) then
      do jy=iy_first,iy_last
      do jz=iz_first,iz_last      
      do jx=ix_first,ix_last
      if(psi(jx,jz).lt.psia1 .and. x(jx,jz,jy,2).gt.0) then
      x(jx,jz,jy,1)=xint(jx,jz,1)*(xint(jx,jz,2)/x(jx,jz,jy,2))**gamma
      else
      x(jx,jz,jy,1)=xint(jx,jz,1)
      endif

      enddo
      enddo
      enddo

      endif

      call convt
      call current(1)
!      call convtdivb
      if(Ef_mode) call Efield
      do jy=iy_first,iy_last
      do jz=iz_first,iz_last      
      do jx=ix_first,ix_last
      Rvx(jx,jz,jy)=xx(jx)*x(jx,jz,jy,3)
      REy(jx,jz,jy)=xx(jx)*Ef(jx,jz,jy,2)
      enddo
      enddo
      enddo

      do 1 jy=iy_first+2,iy_last-2
      do 1 jz=iz_first+2,iz_last-2      
      do 1 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.psia1) then
        dRvx_dx=d1fc(Rvx(jx-2,jz,jy),Rvx(jx-1,jz,jy),Rvx(jx,jz,jy),Rvx(jx+1,jz,jy),Rvx(jx+2,jz,jy),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
      xdif(jx,jz,jy,1)=-x1(jx,jz,jy,3)*xr(jx,jz,jy,1) &
       -x(jx,jz,jy,1)*dRvx_dx/xx(jx) &
!rplc       -x(jx,jz,jy,3)*x(jx,jz,jy,1)/xx(jx)-x(jx,jz,jy,1)*xr(jx,jz,jy,3) &
       -(xy(jx,jz,jy,1)*x(jx,jz,jy,4)+x(jx,jz,jy,1)*xy(jx,jz,jy,4))/xx(jx) &
       -x1(jx,jz,jy,5)*xz(jx,jz,jy,1)-x(jx,jz,jy,1)*x1z(jx,jz,jy,5) &
!ws:for poloidal flow
       -xint(jx,jz,3)*x1r(jx,jz,jy,1)-x1(jx,jz,jy,1)*(xint(jx,jz,3)/xx(jx)+xint_dx(jx,jz,3)) &
       -xint(jx,jz,5)*x1z(jx,jz,jy,1)-x1(jx,jz,jy,1)*xint_dz(jx,jz,5)

      xdif(jx,jz,jy,2)=-x1(jx,jz,jy,3)*xr(jx,jz,jy,2) &
       -gamma*x(jx,jz,jy,2)*(dRvx_dx/xx(jx)) &
!rplc       -gamma*x(jx,jz,jy,2)*(xr(jx,jz,jy,3)+x(jx,jz,jy,3)/xx(jx)) &
       -(xy(jx,jz,jy,2)*x(jx,jz,jy,4)+gamma*x(jx,jz,jy,2)*xy(jx,jz,jy,4))/xx(jx) &
       -x1(jx,jz,jy,5)*xz(jx,jz,jy,2)-gamma*x(jx,jz,jy,2)*x1z(jx,jz,jy,5) &    
!ws:for poloidal flow
       -xint(jx,jz,3)*x1r(jx,jz,jy,2)-gamma*x1(jx,jz,jy,2)*(xint(jx,jz,3)/xx(jx)+xint_dx(jx,jz,3)) &
       -xint(jx,jz,5)*x1z(jx,jz,jy,2)-gamma*x1(jx,jz,jy,2)*xint_dz(jx,jz,5)       
          
!       +(gamma-1)*eta(jx,jz,jy)*(cur(jx,jz,jy,1)**2+cur(jx,jz,jy,2)**2+cur(jx,jz,jy,3)**2)
     
!      xdif(jx,jz,jy,3)=-x(jx,jz,jy,3)*xr(jx,jz,jy,3)-x(jx,jz,jy,5)*xz(jx,jz,jy,3) &
!       -xy(jx,jz,jy,3)*x(jx,jz,jy,4)/xx(jx)+x(jx,jz,jy,4)**2/xx(jx) &
!       +(cur(jx,jz,jy,2)*x(jx,jz,jy,8)+cint(jx,jz,2)*(x(jx,jz,jy,8)-xint(jx,jz,8)) &
!       - cur(jx,jz,jy,3)*x(jx,jz,jy,7)-cint(jx,jz,3)*(x(jx,jz,jy,7)-xint(jx,jz,7)) &
!!       +((cur(jx,jz,jy,2)-cint(jx,jz,2))*x(jx,jz,jy,8)+cint(jx,jz,2)*(x(jx,jz,jy,8)-xint(jx,jz,8)) &
!!       - (cur(jx,jz,jy,3)-cint(jx,jz,3))*x(jx,jz,jy,7)-cint(jx,jz,3)*(x(jx,jz,jy,7)-xint(jx,jz,7)) &
!!!       -d1fc(xh(jx-2,jz,jy,2),xh(jx-1,jz,jy,2),xh(jx,jz,jy,2), &
!!!       xh(jx+1,jz,jy,2),xh(jx+2,jz,jy,2),ax1(jx),bx1(jx), &
!!!       cx1(jx),dx1(jx))
!!!       -d1f2(xh(jx-1,jz,jy,2),xh(jx,jz,jy,2),xh(jx+1,jz,jy,2),xx(jx-1),xx(jx),xx(jx+1)) &
!       -x1r(jx,jz,jy,2)-xint(jx,jz,1)*xint(jx,jz,4)**2/xx(jx))/x(jx,jz,jy,1)
             
      xdif(jx,jz,jy,3)=-x(jx,jz,jy,3)*x1r(jx,jz,jy,3)-x(jx,jz,jy,5)*x1z(jx,jz,jy,3) &
       -x(jx,jz,jy,4)*xy(jx,jz,jy,3)/xx(jx)+x(jx,jz,jy,4)*x1(jx,jz,jy,4)/xx(jx) &
       +(cur(jx,jz,jy,2)*x(jx,jz,jy,8)+cint(jx,jz,2)*x1(jx,jz,jy,8) &
       - cur(jx,jz,jy,3)*x(jx,jz,jy,7)-cint(jx,jz,3)*x1(jx,jz,jy,7) &
       -x1r(jx,jz,jy,2))/x(jx,jz,jy,1) &
       -x1(jx,jz,jy,3)*xint_dx(jx,jz,3)-x1(jx,jz,jy,5)*xint_dz(jx,jz,3) &
       +x1(jx,jz,jy,4)*xint(jx,jz,4)/xx(jx) &
       +(-xint(jx,jz,3)*xint_dx(jx,jz,3)-xint(jx,jz,5)*xint_dz(jx,jz,3) &
       +xint(jx,jz,4)*xint(jx,jz,4)/xx(jx))*x1(jx,jz,jy,1)/x(jx,jz,jy,1)

!      xdif(jx,jz,jy,4)=-x(jx,jz,jy,3)*xr(jx,jz,jy,4)-x(jx,jz,jy,5)*xz(jx,jz,jy,4) &
!       -xy(jx,jz,jy,4)*x(jx,jz,jy,4)/xx(jx)-x(jx,jz,jy,3)*x(jx,jz,jy,4)/xx(jx) &
!       +(cur(jx,jz,jy,3)*x(jx,jz,jy,6)+cint(jx,jz,3)*(x(jx,jz,jy,6)-xint(jx,jz,6)) &
!       - cur(jx,jz,jy,1)*x(jx,jz,jy,8)-cint(jx,jz,1)*(x(jx,jz,jy,8)-xint(jx,jz,8)) &   
!!       +((cur(jx,jz,jy,3)-cint(jx,jz,3))*x(jx,jz,jy,6)+cint(jx,jz,3)*(x(jx,jz,jy,6)-xint(jx,jz,6)) &
!!       - (cur(jx,jz,jy,1)-cint(jx,jz,1))*x(jx,jz,jy,8)-cint(jx,jz,1)*(x(jx,jz,jy,8)-xint(jx,jz,8)) &           
!       -xy(jx,jz,jy,2)/xx(jx))/x(jx,jz,jy,1)

      xdif(jx,jz,jy,4)=-x(jx,jz,jy,3)*x1r(jx,jz,jy,4)-x(jx,jz,jy,5)*x1z(jx,jz,jy,4) &
       -x(jx,jz,jy,4)*xy(jx,jz,jy,4)/xx(jx)-x(jx,jz,jy,4)*x1(jx,jz,jy,3)/xx(jx) &
       +(cur(jx,jz,jy,3)*x(jx,jz,jy,6)+cint(jx,jz,3)*x1(jx,jz,jy,6) &
       - cur(jx,jz,jy,1)*x(jx,jz,jy,8)-cint(jx,jz,1)*x1(jx,jz,jy,8) &   
       -xy(jx,jz,jy,2)/xx(jx))/x(jx,jz,jy,1) &
       -x1(jx,jz,jy,3)*xint_dx(jx,jz,4)-x1(jx,jz,jy,5)*xint_dz(jx,jz,4) &
       -x1(jx,jz,jy,4)*xint(jx,jz,3)/xx(jx) &
       +(-xint(jx,jz,3)*xint_dx(jx,jz,4)-xint(jx,jz,5)*xint_dz(jx,jz,4) &
       -xint(jx,jz,4)*xint(jx,jz,3)/xx(jx))*x1(jx,jz,jy,1)/x(jx,jz,jy,1)
     
!      xdif(jx,jz,jy,5)=-x(jx,jz,jy,3)*xr(jx,jz,jy,5)-x(jx,jz,jy,5)*xz(jx,jz,jy,5) &
!       -xy(jx,jz,jy,5)*x(jx,jz,jy,4)/xx(jx) &     
!       +(cur(jx,jz,jy,1)*x(jx,jz,jy,7)+cint(jx,jz,1)*(x(jx,jz,jy,7)-xint(jx,jz,7)) &
!       - cur(jx,jz,jy,2)*x(jx,jz,jy,6)-cint(jx,jz,2)*(x(jx,jz,jy,6)-xint(jx,jz,6)) & 
!!       +((cur(jx,jz,jy,1)-cint(jx,jz,1))*x(jx,jz,jy,7)-cint(jx,jz,1)*(x(jx,jz,jy,7)-xint(jx,jz,7)) &
!!       - (cur(jx,jz,jy,2)-cint(jx,jz,2))*x(jx,jz,jy,6)-cint(jx,jz,2)*(x(jx,jz,jy,6)-xint(jx,jz,6)) &     
!!!       -d1fc(xh(jx,jz-2,jy,2),xh(jx,jz-1,jy,2),xh(jx,jz,jy,2), &
!!!       xh(jx,jz+1,jy,2),xh(jx,jz+2,jy,2),az1(jz),bz1(jz), &
!!!       cz1(jz),dz1(jz))
!!!       -d1f2(xh(jx,jz-1,jy,2),xh(jx,jz,jy,2),xh(jx,jz+1,jy,2),zz(jz-1),zz(jz),zz(jz+1)) &
!       -x1z(jx,jz,jy,2))/x(jx,jz,jy,1)
       
      xdif(jx,jz,jy,5)=-x(jx,jz,jy,3)*x1r(jx,jz,jy,5)-x(jx,jz,jy,5)*x1z(jx,jz,jy,5) &
       -x(jx,jz,jy,4)*xy(jx,jz,jy,5)/xx(jx) &     
       +(cur(jx,jz,jy,1)*x(jx,jz,jy,7)+cint(jx,jz,1)*x1(jx,jz,jy,7) &
       - cur(jx,jz,jy,2)*x(jx,jz,jy,6)-cint(jx,jz,2)*x1(jx,jz,jy,6) & 
       -x1z(jx,jz,jy,2))/x(jx,jz,jy,1) &
       -x1(jx,jz,jy,3)*xint_dx(jx,jz,5)-x1(jx,jz,jy,5)*xint_dz(jx,jz,5) &
       +(-xint(jx,jz,3)*xint_dx(jx,jz,5)-xint(jx,jz,5)*xint_dz(jx,jz,5)) &
       *x1(jx,jz,jy,1)/x(jx,jz,jy,1)
          
     
!      xdif(jx,jz,jy,6)=-x(jx,jz,jy,3)*xr(jx,jz,jy,6)-x(jx,jz,jy,6)*x(jx,jz,jy,3)/xx(jx) &
!       -(x(jx,jz,jy,4)*xy(jx,jz,jy,6)-x(jx,jz,jy,7)*xy(jx,jz,jy,3)+x(jx,jz,jy,6)*xy(jx,jz,jy,4))/xx(jx) &
!       -x(jx,jz,jy,5)*xz(jx,jz,jy,6)+x(jx,jz,jy,8)*xz(jx,jz,jy,3)-x(jx,jz,jy,6)*xz(jx,jz,jy,5) !&
!!J X grd_eta
!!       -(cur(jx,jz,jy,2)-cj(jx,jz))*etaz(jx,jz,jy)
!     
!      xdif(jx,jz,jy,7)=-x(jx,jz,jy,3)*xr(jx,jz,jy,7)+x(jx,jz,jy,6)*xr(jx,jz,jy,4) &
!       -x(jx,jz,jy,6)*x(jx,jz,jy,4)/xx(jx)-x(jx,jz,jy,7)*xr(jx,jz,jy,3) &
!       -xy(jx,jz,jy,7)*x(jx,jz,jy,4)/xx(jx) &
!       -xz(jx,jz,jy,7)*x(jx,jz,jy,5)+x(jx,jz,jy,8)*xz(jx,jz,jy,4)-x(jx,jz,jy,7)*xz(jx,jz,jy,5) !&
!!J X grd_eta
!!       -cur(jx,jz,jy,3)*etax(jx,jz,jy)+cur(jx,jz,jy,1)*etaz(jx,jz,jy)
!      
!      xdif(jx,jz,jy,8)=-x(jx,jz,jy,3)*xr(jx,jz,jy,8)-x(jx,jz,jy,3)*x(jx,jz,jy,8)/xx(jx) &
!       -x(jx,jz,jy,8)*xr(jx,jz,jy,3)+x(jx,jz,jy,6)*xr(jx,jz,jy,5) &
!       -(xy(jx,jz,jy,8)*x(jx,jz,jy,4)-xy(jx,jz,jy,5)*x(jx,jz,jy,7)+xy(jx,jz,jy,4)*x(jx,jz,jy,8))/xx(jx) &
!       -xz(jx,jz,jy,8)*x(jx,jz,jy,5) !&
!!J X grd_eta
!!       +(cur(jx,jz,jy,2)-cj(jx,jz))*etax(jx,jz,jy)

!      xdif(jx,jz,jy,6)=-Efy(jx,jz,jy,3)/xx(jx)+Efz(jx,jz,jy,2) !&
!!                      +cdb0*divb_x(jx,jz,jy)
!    
!      xdif(jx,jz,jy,7)=-Efz(jx,jz,jy,1)+Efx(jx,jz,jy,3) !&
!!                      +cdb0*divb_y(jx,jz,jy)/xx(jx)
!
!      xdif(jx,jz,jy,8)=Efy(jx,jz,jy,1)/xx(jx)-Efx(jx,jz,jy,2)-Ef(jx,jz,jy,2)/xx(jx) !&
!!                      +cdb0*divb_z(jx,jz,jy)

      if(Ef_mode) then
      !dREy_dx=d1xf2(Ef(jx-1,jz,jy,2),Ef(jx,jz,jy,2),Ef(jx+1,jz,jy,2),xx(jx-1),xx(jx),xx(jx+1))
      !dEz_dx =d1f2(Ef(jx-1,jz,jy,3),Ef(jx,jz,jy,3),Ef(jx+1,jz,jy,3),xx(jx-1),xx(jx),xx(jx+1))
      !dEx_dz =d1f2(Ef(jx,jz-1,jy,1),Ef(jx,jz,jy,1),Ef(jx,jz+1,jy,1),zz(jz-1),zz(jz),zz(jz+1))
      !dEy_dz =d1f2(Ef(jx,jz-1,jy,2),Ef(jx,jz,jy,2),Ef(jx,jz+1,jy,2),zz(jz-1),zz(jz),zz(jz+1))
      !dEx_dy =d1f2(Ef(jx,jz,jy-1,1),Ef(jx,jz,jy,1),Ef(jx,jz,jy+1,1),yy(jy-1),yy(jy),yy(jy+1))
      !dEz_dy =d1f2(Ef(jx,jz,jy-1,3),Ef(jx,jz,jy,3),Ef(jx,jz,jy+1,3),yy(jy-1),yy(jy),yy(jy+1))

     dREy_dx=d1fc(REy(jx-2,jz,jy),REy(jx-1,jz,jy),REy(jx,jz,jy) &
         ,REy(jx+1,jz,jy),REy(jx+2,jz,jy),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
     dEz_dx =d1fc(Ef(jx-2,jz,jy,3),Ef(jx-1,jz,jy,3),Ef(jx,jz,jy,3) &
         ,Ef(jx+1,jz,jy,3),Ef(jx+2,jz,jy,3),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
     dEx_dz =d1fc(Ef(jx,jz-2,jy,1),Ef(jx,jz-1,jy,1),Ef(jx,jz,jy,1) &
         ,Ef(jx,jz+1,jy,1),Ef(jx,jz+2,jy,1),az1(jz),bz1(jz),cz1(jz),dz1(jz))
     dEy_dz =d1fc(Ef(jx,jz-2,jy,2),Ef(jx,jz-1,jy,2),Ef(jx,jz,jy,2) &
         ,Ef(jx,jz+1,jy,2),Ef(jx,jz+2,jy,2),az1(jz),bz1(jz),cz1(jz),dz1(jz))
     dEx_dy =d1fc(Ef(jx,jz,jy-2,1),Ef(jx,jz,jy-1,1),Ef(jx,jz,jy,1) &
         ,Ef(jx,jz,jy+1,1),Ef(jx,jz,jy+2,1),ay1(jy),by1(jy),cy1(jy),dy1(jy)) 
     dEz_dy =d1fc(Ef(jx,jz,jy-2,3),Ef(jx,jz,jy-1,3),Ef(jx,jz,jy,3) &
         ,Ef(jx,jz,jy+1,3),Ef(jx,jz,jy+2,3),ay1(jy),by1(jy),cy1(jy),dy1(jy)) 
          
      xdif(jx,jz,jy,6)=-dEz_dy/xx(jx)+dEy_dz !&
!                       +eta(jx,jz,jy)*cint_dz(jx,jz,2)  
      xdif(jx,jz,jy,7)=-dEx_dz+dEz_dx !&
!                       +eta(jx,jz,jy)*(-cint_dz(jx,jz,1)+cint_dx(jx,jz,3))
      xdif(jx,jz,jy,8)=(dEx_dy-dREy_dx)/xx(jx) !&
!                       +eta(jx,jz,jy)*(-cint_dx(jx,jz,2)-cint(jx,jz,2)/xx(jx))

      else
      xdif(jx,jz,jy,6)=(x(jx,jz,jy,3)*xy(jx,jz,jy,7)+x(jx,jz,jy,7)*xy(jx,jz,jy,3) &
                       -x(jx,jz,jy,4)*xy(jx,jz,jy,6)-x(jx,jz,jy,6)*xy(jx,jz,jy,4))/xx(jx) &
                      -(x(jx,jz,jy,5)*xz(jx,jz,jy,6)+x(jx,jz,jy,6)*xz(jx,jz,jy,5) &
                       -x(jx,jz,jy,3)*xz(jx,jz,jy,8)-x(jx,jz,jy,8)*xz(jx,jz,jy,3))


      xdif(jx,jz,jy,7)=(x(jx,jz,jy,4)*xz(jx,jz,jy,8)+x(jx,jz,jy,8)*xz(jx,jz,jy,4) &
                       -x(jx,jz,jy,5)*xz(jx,jz,jy,7)-x(jx,jz,jy,7)*xz(jx,jz,jy,5)) &
                      -(x(jx,jz,jy,3)*xr(jx,jz,jy,7)+x(jx,jz,jy,7)*xr(jx,jz,jy,3) &
                       -x(jx,jz,jy,4)*xr(jx,jz,jy,6)-x(jx,jz,jy,6)*xr(jx,jz,jy,4))
      
!      xdif(jx,jz,jy,8)=(x(jx,jz,jy,5)*xr(jx,jz,jy,6)+x(jx,jz,jy,6)*xr(jx,jz,jy,5) &
!                       -x(jx,jz,jy,3)*xr(jx,jz,jy,8)-x(jx,jz,jy,8)*xr(jx,jz,jy,3)) &    
!                      -(x(jx,jz,jy,3)*x(jx,jz,jy,8)-x(jx,jz,jy,5)*x(jx,jz,jy,6) &
!                       +x(jx,jz,jy,4)*xy(jx,jz,jy,8)+x(jx,jz,jy,8)*xy(jx,jz,jy,4) &
!                       -x(jx,jz,jy,5)*xy(jx,jz,jy,7)-x(jx,jz,jy,7)*xy(jx,jz,jy,5))/xx(jx)

      xdif(jx,jz,jy,8)=(x(jx,jz,jy,5)*dRBx_dx/xx(jx)+x(jx,jz,jy,6)*xr(jx,jz,jy,5) &
                       -x(jx,jz,jy,3)*xr(jx,jz,jy,8)-x(jx,jz,jy,8)*dRvx_dx/xx(jx)) &    
                      +(x(jx,jz,jy,4)*xy(jx,jz,jy,8)+x(jx,jz,jy,8)*xy(jx,jz,jy,4) &
                       -x(jx,jz,jy,5)*xy(jx,jz,jy,7)-x(jx,jz,jy,7)*xy(jx,jz,jy,5))/xx(jx)

       endif



       else    
       do m=1,8
       xdif(jx,jz,jy,m)=0.
       enddo
!
       endif

    1 continue
       
!       do m=1,8
!       do jy=1,my
!       do jx=jxamin,jxamax
!       xdif(jx,jzam(jx)-1,jy,m)=0
!       xdif(jx,jzap(jx)+1,jy,m)=0
!       enddo
!       do jz=jzamin,jzamax
!       xdif(jxam(jz)-1,jz,jy,m)=0
!       xdif(jxap(jz)+1,jz,jy,m)=0
!       enddo
!
!       enddo
!       enddo

!      write(*,*) 'ws'

!      if(resisitive) then
      if(Ef_mode .and. smoothbout) then
!      if(.not.implicitb) then
      do 41 jy=iy_first+2,iy_last-2
      do 41 jz=iz_first+2,iz_last-2
      do 41 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.psia1) then
!!ws:dB/dt=...+eta*grd2 B
      xdif(jx,jz,jy,6)=xdif(jx,jz,jy,6)+cur(jx,jz,jy,2)*etbz(jx,jz) &
       +etb(jx,jz)*(xr2(jx,jz,jy,6)+x1r(jx,jz,jy,6)/xx(jx)+xy2(jx,jz,jy,6)/xx(jx)**2+xz2(jx,jz,jy,6) &
        -x1(jx,jz,jy,6)/xx(jx)**2-2.0*xy(jx,jz,jy,7)/xx(jx)**2)
     
      xdif(jx,jz,jy,7)=xdif(jx,jz,jy,7)+cur(jx,jz,jy,3)*etbx(jx,jz)-cur(jx,jz,jy,1)*etaz(jx,jz,jy) &
       +etb(jx,jz)*(xr2(jx,jz,jy,7)+x1r(jx,jz,jy,7)/xx(jx)+xy2(jx,jz,jy,7)/xx(jx)**2+xz2(jx,jz,jy,7) &
        -x1(jx,jz,jy,7)/xx(jx)**2+2.0*xy(jx,jz,jy,6)/xx(jx)**2)
            
      xdif(jx,jz,jy,8)=xdif(jx,jz,jy,8)-cur(jx,jz,jy,2)*etbx(jx,jz) &
       +etb(jx,jz)*(xr2(jx,jz,jy,8)+x1r(jx,jz,jy,8)/xx(jx)+xy2(jx,jz,jy,8)/xx(jx)**2+xz2(jx,jz,jy,8))
      endif
   41 continue
      endif


      if(resisitive) then
      if(.not.etaJ_in_E) then
!      if(.not.implicitb) then
      do 4 jy=iy_first+2,iy_last-2
      do 4 jz=iz_first+2,iz_last-2
      do 4 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.psia1) then
!      xdif(jx,jz,jy,2)=xdif(jx,jz,jy,2) &       
!       +(gamma-1)*eta(jx,jz,jy)*(cur(jx,jz,jy,1)**2+cur(jx,jz,jy,2)**2+cur(jx,jz,jy,3)**2)

!      xdif(jx,jz,jy,6)=xdif(jx,jz,jy,6)-(cur(jx,jz,jy,2)-cint(jx,jz,2))*etaz(jx,jz,jy) &
!       +eta(jx,jz,jy)*(cuy(jx,jz,jy,3)/xx(jx)-cuz(jx,jz,jy,2))
!     
!      xdif(jx,jz,jy,7)=xdif(jx,jz,jy,7)-cur(jx,jz,jy,3)*etax(jx,jz,jy)+cur(jx,jz,jy,1)*etaz(jx,jz,jy) &
!       +eta(jx,jz,jy)*(cuz(jx,jz,jy,1)-cux(jx,jz,jy,3))
!            
!      xdif(jx,jz,jy,8)=xdif(jx,jz,jy,8)+(cur(jx,jz,jy,2)-cint(jx,jz,2))*etax(jx,jz,jy) &
!       +eta(jx,jz,jy)*(cux(jx,jz,jy,2)+(cur(jx,jz,jy,2)-cint(jx,jz,2))/xx(jx)-cuy(jx,jz,jy,1)/xx(jx))

!      xdif(jx,jz,jy,6)=xdif(jx,jz,jy,6)+(cur(jx,jz,jy,2)+cint(jx,jz,2))*etaz(jx,jz,jy) &
!       -eta(jx,jz,jy)*(cuy(jx,jz,jy,3)/xx(jx)-cuz(jx,jz,jy,2)-cint_dz(jx,jz,2))
!     
!      xdif(jx,jz,jy,7)=xdif(jx,jz,jy,7)+(cur(jx,jz,jy,3)+cint(jx,jz,3))*etax(jx,jz,jy)-(cur(jx,jz,jy,1)+cint(jx,jz,1))*etaz(jx,jz,jy) &
!       -eta(jx,jz,jy)*(cuz(jx,jz,jy,1)-cux(jx,jz,jy,3)+cint_dz(jx,jz,1)-cint_dx(jx,jz,3))
!            
!      xdif(jx,jz,jy,8)=xdif(jx,jz,jy,8)-(cur(jx,jz,jy,2)+cint(jx,jz,2))*etax(jx,jz,jy) &
!       -eta(jx,jz,jy)*(cux(jx,jz,jy,2)+cint_dx(jx,jz,2)+(cur(jx,jz,jy,2)+cint(jx,jz,2))/xx(jx)-cuy(jx,jz,jy,1)/xx(jx))

!!ws:dB/dt=...-eta*grd X J
!      xdif(jx,jz,jy,6)=xdif(jx,jz,jy,6)+cur(jx,jz,jy,2)*etaz(jx,jz,jy) &
!       -eta(jx,jz,jy)*(cuy(jx,jz,jy,3)/xx(jx)-cuz(jx,jz,jy,2))
!     
!      xdif(jx,jz,jy,7)=xdif(jx,jz,jy,7)+cur(jx,jz,jy,3)*etax(jx,jz,jy)-cur(jx,jz,jy,1)*etaz(jx,jz,jy) &
!       -eta(jx,jz,jy)*(cuz(jx,jz,jy,1)-cux(jx,jz,jy,3))
!            
!      xdif(jx,jz,jy,8)=xdif(jx,jz,jy,8)-cur(jx,jz,jy,2)*etax(jx,jz,jy) &
!       -eta(jx,jz,jy)*(cux(jx,jz,jy,2)+cur(jx,jz,jy,2)/xx(jx)-cuy(jx,jz,jy,1)/xx(jx))
!!ws:dB/dt=...+eta*grd2 B
      xdif(jx,jz,jy,6)=xdif(jx,jz,jy,6)+cur(jx,jz,jy,2)*etaz(jx,jz,jy) &
       +eta(jx,jz,jy)*(xr2(jx,jz,jy,6)+x1r(jx,jz,jy,6)/xx(jx)+xy2(jx,jz,jy,6)/xx(jx)**2+xz2(jx,jz,jy,6) &
        -x1(jx,jz,jy,6)/xx(jx)**2-2.0*xy(jx,jz,jy,7)/xx(jx)**2)
     
      xdif(jx,jz,jy,7)=xdif(jx,jz,jy,7)+cur(jx,jz,jy,3)*etax(jx,jz,jy)-cur(jx,jz,jy,1)*etaz(jx,jz,jy) &
       +eta(jx,jz,jy)*(xr2(jx,jz,jy,7)+x1r(jx,jz,jy,7)/xx(jx)+xy2(jx,jz,jy,7)/xx(jx)**2+xz2(jx,jz,jy,7) &
        -x1(jx,jz,jy,7)/xx(jx)**2+2.0*xy(jx,jz,jy,6)/xx(jx)**2)
            
      xdif(jx,jz,jy,8)=xdif(jx,jz,jy,8)-cur(jx,jz,jy,2)*etax(jx,jz,jy) &
       +eta(jx,jz,jy)*(xr2(jx,jz,jy,8)+x1r(jx,jz,jy,8)/xx(jx)+xy2(jx,jz,jy,8)/xx(jx)**2+xz2(jx,jz,jy,8))

!      xdif(jx,jz,jy,6)=xdif(jx,jz,jy,6)+cur(jx,jz,jy,2)*etaz(jx,jz,jy)+eta(jx,jz,jy) &
!       *((xy2(jx,jz,jy,6)-xy(jx,jz,jy,7))/xx(jx)**2+xz2(jx,jz,jy,6) &
!        -(by_yx(jx,jz,jy)/xx(jx)+bz_zx(jx,jz,jy)) )
!     
!      xdif(jx,jz,jy,7)=xdif(jx,jz,jy,7)+cur(jx,jz,jy,3)*etax(jx,jz,jy)-cur(jx,jz,jy,1)*etaz(jx,jz,jy)+eta(jx,jz,jy) &
!       *(xr2(jx,jz,jy,7)+x1r(jx,jz,jy,7)/xx(jx)+xz2(jx,jz,jy,7) &
!        -(x1(jx,jz,jy,7)-xy(jx,jz,jy,6))/xx(jx)**2 &
!        -(bx_xy(jx,jz,jy)+bz_zy(jx,jz,jy))/xx(jx) )
!            
!      xdif(jx,jz,jy,8)=xdif(jx,jz,jy,8)-cur(jx,jz,jy,2)*etax(jx,jz,jy)+eta(jx,jz,jy) &
!       *(xr2(jx,jz,jy,8)+x1r(jx,jz,jy,8)/xx(jx)+xy2(jx,jz,jy,8)/xx(jx)**2 &
!        -(bx_xz(jx,jz,jy)+(by_yz(jx,jz,jy)+x1z(jx,jz,jy,6))/xx(jx)) )

      if(nstep.lt.nper) then
      xdif(jx,jz,jy,6)=xdif(jx,jz,jy,6)+perb(jx,jz,jy,1)
      xdif(jx,jz,jy,7)=xdif(jx,jz,jy,7)+perb(jx,jz,jy,2)
      xdif(jx,jz,jy,8)=xdif(jx,jz,jy,8)+perb(jx,jz,jy,3)  
!      xdif(jx,jz,jy,6)=xdif(jx,jz,jy,6)+cint(jx,jz,2)*eta1z(jx,jz,jy)-cint(jx,jz,3)*eta1y(jx,jz,jy) &
!        -eta1(jx,jz,jy)*(-cint_dz(jx,jz,2))
!      xdif(jx,jz,jy,7)=xdif(jx,jz,jy,7)+cint(jx,jz,3)*eta1x(jx,jz,jy)-cint(jx,jz,1)*eta1z(jx,jz,jy) &
!        -eta1(jx,jz,jy)*(cint_dz(jx,jz,1)-cint_dx(jx,jz,3))
!      xdif(jx,jz,jy,8)=xdif(jx,jz,jy,8)+cint(jx,jz,1)*eta1y(jx,jz,jy)-cint(jx,jz,2)*eta1x(jx,jz,jy) &
!        -eta1(jx,jz,jy)*(cint_dx(jx,jz,2)+cint(jx,jz,2)/xx(jx))  
      endif
      
      endif
    4 continue
 

!      do 4 jy=1,my
!      do 4 jz=iz_first+1,iz_last-1
!      do 4 jx=ix_first+1,ix_last-1
!!      if(rr(jx,jz).lt.aa) then
!      xdif(jx,jz,jy,2)=xdif(jx,jz,jy,2) &       
!       +(gamma-1)*eta(jx,jz,jy)*(cur(jx,jz,jy,1)**2+cur(jx,jz,jy,2)**2+cur(jx,jz,jy,3)**2)
!
!      xdif(jx,jz,jy,6)=xdif(jx,jz,jy,6)-(cur(jx,jz,jy,2)-cj(jx,jz))*etaz(jx,jz,jy) &
!       +eta(jx,jz,jy)*(xr2(jx,jz,jy,6)+xr(jx,jz,jy,6)/xx(jx)+xz2(jx,jz,jy,6) &
!       +(xy2(jx,jz,jy,6)-x(jx,jz,jy,6)-2*xy(jx,jz,jy,7))/xx(jx)**2)
!     
!      xdif(jx,jz,jy,7)=xdif(jx,jz,jy,7)-cur(jx,jz,jy,3)*etax(jx,jz,jy)+cur(jx,jz,jy,1)*etaz(jx,jz,jy) &
!       +eta(jx,jz,jy)*(xr2(jx,jz,jy,7)+xr(jx,jz,jy,7)/xx(jx)+xz2(jx,jz,jy,7) &
!       +(xy2(jx,jz,jy,7)-x(jx,jz,jy,7)+2*xy(jx,jz,jy,6))/xx(jx)**2)
!            
!      xdif(jx,jz,jy,8)=xdif(jx,jz,jy,8)+(cur(jx,jz,jy,2)-cj(jx,jz))*etax(jx,jz,jy) &
!       +eta(jx,jz,jy)*(xr2(jx,jz,jy,8)+xr(jx,jz,jy,8)/xx(jx)+xz2(jx,jz,jy,8) &
!       +xy2(jx,jz,jy,8)/xx(jx)**2)
!!      endif
!    4 continue
      endif
      endif
!      call bndry_xdif      
            
!      if(resisitive) then   
!      do 42 jy=1,my
!      do 42 jz=iz_first,iz_last
!      do 42 jx=ix_first,ix_last
!      if(psi(jx,jz).lt.psia1) then
!      xdif(jx,jz,jy,2)=xdif(jx,jz,jy,2) &       
!       +(gamma-1)*eta(jx,jz,jy)*(cur(jx,jz,jy,1)**2+cur(jx,jz,jy,2)**2+cur(jx,jz,jy,3)**2)
!      endif
!   42 continue
!      call bndry_xdif_p
!      endif

      if(divb) then
      do 2 jy=iy_first+2,iy_last-2
      do 2 jz=iz_first+2,iz_last-2
      do 2 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.psia1) then
      xdif(jx,jz,jy,6)=xdif(jx,jz,jy,6)+cdb(jx,jz)*divb_x(jx,jz,jy)+cdbx(jx,jz)*dvb(jx,jz,jy)      
    
      xdif(jx,jz,jy,7)=xdif(jx,jz,jy,7)+cdb(jx,jz)*divb_y(jx,jz,jy)/xx(jx)

      xdif(jx,jz,jy,8)=xdif(jx,jz,jy,8)+cdb(jx,jz)*divb_z(jx,jz,jy)+cdbz(jx,jz)*dvb(jx,jz,jy)
      endif
    2 continue
      endif

      if(viscous) then     
      if(.not.implicitv)then
!      call vorticity
      do 5 jy=iy_first+2,iy_last-2
      do 5 jz=iz_first+2,iz_last-2      
      do 5 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.psia1) then
      xdif(jx,jz,jy,1)=xdif(jx,jz,jy,1)+pmu(jx,jz) &
       *(xr2(jx,jz,jy,1)+x1r(jx,jz,jy,1)/xx(jx)+xy2(jx,jz,jy,1)/xx(jx)**2+xz2(jx,jz,jy,1)) &
       +pmux(jx,jz)*x1r(jx,jz,jy,1)+pmuz(jx,jz)*x1z(jx,jz,jy,1) 

      xdif(jx,jz,jy,2)=xdif(jx,jz,jy,2)+kap(jx,jz) &
       *(xr2(jx,jz,jy,2)+x1r(jx,jz,jy,2)/xx(jx)+xy2(jx,jz,jy,2)/xx(jx)**2+xz2(jx,jz,jy,2)) &
       +kapx(jx,jz)*x1r(jx,jz,jy,2)+kapz(jx,jz)*x1z(jx,jz,jy,2)

      xdif(jx,jz,jy,3)=xdif(jx,jz,jy,3)+fmu(jx,jz) &
       *(xr2(jx,jz,jy,3)+x1r(jx,jz,jy,3)/xx(jx)+xy2(jx,jz,jy,3)/xx(jx)**2+xz2(jx,jz,jy,3) &
        -x1(jx,jz,jy,3)/xx(jx)**2-2.0*xy(jx,jz,jy,4)/xx(jx)**2) &
       +fmux(jx,jz)*x1r(jx,jz,jy,3)+fmuz(jx,jz)*x1z(jx,jz,jy,3)
     
      xdif(jx,jz,jy,4)=xdif(jx,jz,jy,4)+fmu(jx,jz) &
       *(xr2(jx,jz,jy,4)+x1r(jx,jz,jy,4)/xx(jx)+xy2(jx,jz,jy,4)/xx(jx)**2+xz2(jx,jz,jy,4) &
        -x1(jx,jz,jy,4)/xx(jx)**2+2.0*xy(jx,jz,jy,3)/xx(jx)**2) &
       +fmux(jx,jz)*x1r(jx,jz,jy,4)+fmuz(jx,jz)*x1z(jx,jz,jy,4)
     
      xdif(jx,jz,jy,5)=xdif(jx,jz,jy,5)+fmu(jx,jz) &
       *(xr2(jx,jz,jy,5)+x1r(jx,jz,jy,5)/xx(jx)+xy2(jx,jz,jy,5)/xx(jx)**2+xz2(jx,jz,jy,5)) &
       +fmux(jx,jz)*x1r(jx,jz,jy,5)+fmuz(jx,jz)*x1z(jx,jz,jy,5) 
      endif
    5 continue
      endif   
      endif

      if(ohm_heat) then
      do 6 jy=iy_first+2,iy_last-2
      do 6 jz=iz_first+2,iz_last-2      
      do 6 jx=ix_first+2,ix_last-2
!      xdif(jx,jz,jy,2)=xdif(jx,jz,jy,2)+eta(jx,jz,jy)*((cur(jx,jz,jy,1)+cint(jx,jz,1))**2 &
!       +(cur(jx,jz,jy,2)+cint(jx,jz,2))**2+(cur(jx,jz,jy,3)+cint(jx,jz,3))**2) &
!       -etaint(jx,jz)*(cint(jx,jz,1)**2+cint(jx,jz,2)**2+cint(jx,jz,3)**2)

      xdif(jx,jz,jy,2)=xdif(jx,jz,jy,2)+eta(jx,jz,jy)*((cur(jx,jz,jy,1)+cint(jx,jz,1))*cur(jx,jz,jy,1) &
       +(cur(jx,jz,jy,2)+cint(jx,jz,2))*cur(jx,jz,jy,2)+(cur(jx,jz,jy,3)+cint(jx,jz,3))*cur(jx,jz,jy,3)) 


    6 continue
      endif
      
!!      if(hall) then
!!      if(pressure) then
!!      do 33 jy=1,my
!!      do 33 jz=iz_first,iz_last      
!!      do 33 jx=ix_first+2,ix_last-2
!!      xdif(jx,jz,jy,6)=xdif(jx,jz,jy,6)-fdi(jx)/x(jx,jz,jy,1)**2 &
!!       /xx(jx)*(xy(jx,jz,jy,1)*xz(jx,jz,jy,2) &
!!       -xz(jx,jz,jy,1)*xy(jx,jz,jy,2))
!!     
!!       xdif(jx,jz,jy,7)=xdif(jx,jz,jy,7)-fdi(jx)/x(jx,jz,jy,1)**2 &
!!       *xz(jx,jz,jy,1)*xr(jx,jz,jy,2)-xz(jx,jz,jy,2) &
!!       *xr(jx,jz,jy,1)-d1fc(fdi(jx-2),fdi(jx-1),fdi(jx),fdi(jx+1), &
!!       fdi(jx+2),ax1(jx),bx1(jx),cx1(jx),dx1(jx)) &
!!       *xz(jx,jz,jy,2)/x(jx,jz,jy,1)
!!      
!!       xdif(jx,jz,jy,8)=xdif(jx,jz,jy,8)-fdi(jx)/x(jx,jz,jy,1)**2 &
!!       *(xy(jx,jz,jy,2)*xr(jx,jz,jy,1)-xy(jx,jz,jy,1)*xr(jx,jz,jy,2)) &
!!       /xx(jx)+d1fc(fdi(jx-2),fdi(jx-1), &
!!       fdi(jx),fdi(jx+1),fdi(jx+2),ax1(jx),bx1(jx),cx1(jx),dx1(jx)) &
!!       *xy(jx,jz,jy,2)/xx(jx)/x(jx,jz,jy,1)
!!   33 continue
!!      else
!!      do 2 jy=1,my
!!      do 2 jz=iz_first,iz_last      
!!      do 2 jx=ix_first+2,ix_last-2
!!       xh(jx,jz,jy,1)=cur(jx,jz,jy,1)*xr(jx,jz,jy,1) &
!!       +cur(jx,jz,jy,2)*xy(jx,jz,jy,1)/xx(jx)+xz(jx,jz,jy,1) &
!!       *cur(jx,jz,jy,3)
!!       xh(jx,jz,jy,2)=x(jx,jz,jy,6)*xr(jx,jz,jy,1) &
!!       +x(jx,jz,jy,7)*xy(jx,jz,jy,1)/xx(jx)+ &
!!       x(jx,jz,jy,8)*xz(jx,jz,jy,1)
!!    2 continue
!!      
!!      do 3 jy=1,my
!!      do 3 jz=iz_first,iz_last
!!      do 3 jx=ix_first+2,ix_last-2
!!      xdif(jx,jz,jy,6)=xdif(jx,jz,jy,6)-fdi(jx)/x(jx,jz,jy,1)**2* &
!!       (xh(jx,jz,jy,1)*x(jx,jz,jy,6)-xh(jx,jz,jy,2)*cur(jx,jz,jy,1) &
!!       +1/xx(jx)*(xy(jx,jz,jy,1)*xz(jx,jz,jy,2) &
!!       -xz(jx,jz,jy,1)*xy(jx,jz,jy,2)))-fdi(jx)/x(jx,jz,jy,1) &
!!       *(x(jx,jz,jy,6)*d1fc(cur(jx-2,jz,jy,1),cur(jx-1,jz,jy,1), &
!!       cur(jx,jz,jy,1),cur(jx+1,jz,jy,1),cur(jx+2,jz,jy,1), &
!!       ax1(jx),bx1(jx),cx1(jx),dx1(jx))-cur(jx,jz,jy,1) &
!!       *xr(jx,jz,jy,6)+(x(jx,jz,jy,7)*cuy(jx,jz,jy,1)-xy(jx,jz,jy,6) &
!!       *cur(jx,jz,jy,2))/xx(jx)+x(jx,jz,jy,8) &
!!       *cuz(jx,jz,jy,1)-xz(jx,jz,jy,6)*cur(jx,jz,jy,3))
!!     
!!       xdif(jx,jz,jy,7)=xdif(jx,jz,jy,7)-fdi(jx)/x(jx,jz,jy,1)**2 &
!!       *(xh(jx,jz,jy,1)*x(jx,jz,jy,7)-xh(jx,jz,jy,2)*cur(jx,jz,jy,2) &
!!       +xz(jx,jz,jy,1)*xr(jx,jz,jy,2)-xz(jx,jz,jy,2) &
!!       *xr(jx,jz,jy,1))-fdi(jx)/x(jx,jz,jy,1)*(x(jx,jz,jy,6)* &
!!       d1fc(cur(jx-2,jz,jy,2),cur(jx-1,jz,jy,2),cur(jx,jz,jy,2), &
!!       cur(jx+1,jz,jy,2),cur(jx+2,jz,jy,2),ax1(jx),bx1(jx),cx1(jx), &
!!       dx1(jx))-cur(jx,jz,jy,1)*xr(jx,jz,jy,7)+(x(jx,jz,jy,7)*cuy(jx,jz,jy,2) &
!!       -xy(jx,jz,jy,7)*cur(jx,jz,jy,2)-x(jx,jz,jy,6)*cur(jx,jz,jy,2) &
!!       +x(jx,jz,jy,7)*cur(jx,jz,jy,1))/xx(jx)+ &
!!       x(jx,jz,jy,8)*cuz(jx,jz,jy,2)-xz(jx,jz,jy,7) &
!!       *cur(jx,jz,jy,3))+d1fc(fdi(jx-2),fdi(jx-1),fdi(jx),fdi(jx+1), &
!!       fdi(jx+2),ax1(jx),bx1(jx),cx1(jx),dx1(jx))*(cur(jx,jz,jy,1) &
!!       *x(jx,jz,jy,7)-cur(jx,jz,jy,2)*x(jx,jz,jy,6)- &
!!       xz(jx,jz,jy,2))/x(jx,jz,jy,1)
!!      
!!       xdif(jx,jz,jy,8)=xdif(jx,jz,jy,8)-fdi(jx)/x(jx,jz,jy,1)**2 &
!!       *(xh(jx,jz,jy,1)*x(jx,jz,jy,8)-xh(jx,jz,jy,2)*cur(jx,jz,jy,3) &
!!       +(xy(jx,jz,jy,2)*xr(jx,jz,jy,1)-xy(jx,jz,jy,1)*xr(jx,jz,jy,2))/xx(jx)) &
!!       -fdi(jx)/x(jx,jz,jy,1)*(x(jx,jz,jy,6)*d1fc(cur(jx-2,jz,jy,3),cur(jx-1,jz,jy,3), &
!!       cur(jx,jz,jy,3),cur(jx+1,jz,jy,3),cur(jx+2,jz,jy,3),ax1(jx), &
!!       bx1(jx),cx1(jx),dx1(jx))-cur(jx,jz,jy,1)*xr(jx,jz,jy,8)+(x(jx,jz,jy,7) &
!!       *cuy(jx,jz,jy,3)-xy(jx,jz,jy,8)*cur(jx,jz,jy,2))/xx(jx) &
!!       +x(jx,jz,jy,8)*cuz(jx,jz,jy,3)-xz(jx,jz,jy,8) &
!!       *cur(jx,jz,jy,3)-d1fc(fdi(jx-2),fdi(jx-1),fdi(jx),fdi(jx+1), &
!!       fdi(jx+2),ax1(jx),bx1(jx),cx1(jx),dx1(jx))*(cur(jx,jz,jy,3) &
!!       *x(jx,jz,jy,6)-cur(jx,jz,jy,1)*x(jx,jz,jy,8) &
!!       -xy(jx,jz,jy,2)/xx(jx))/x(jx,jz,jy,1)
!!    3 continue
!!      endif
!!      endif
!!            
!!      if(.not.implicitb) then
!!      do 4 jy=1,my
!!      do 4 jz=iz_first,iz_last
!!      do 4 jx=ix_first+2,ix_last-2
!!      xdif(jx,jz,jy,6)=xdif(jx,jz,jy,6)-eta(jx) &
!!       *(cuy(jx,jz,jy,3)/xx(jx)-cuz(jx,jz,jy,2))
!!     
!!      xdif(jx,jz,jy,7)=xdif(jx,jz,jy,7)-eta(jx) &
!!       *(cuz(jx,jz,jy,1)-d1fc(cur(jx-2,jz,jy,3),cur(jx-1,jz,jy,3),cur(jx,jz,jy,3), &
!!       cur(jx+1,jz,jy,3),cur(jx+2,jz,jy,3),ax1(jx),bx1(jx),cx1(jx), &
!!       dx1(jx))
!!            
!!      xdif(jx,jz,jy,8)=xdif(jx,jz,jy,8)-eta(jx)/xx(jx)*(cur(jx,jz,jy,2) &
!!       -cuy(jx,jz,jy,1)+xx(jx)*d1fc(cur(jx-2,jz,jy,2),cur(jx-1,jz,jy,2) &
!!       ,cur(jx,jz,jy,2),cur(jx+1,jz,jy,2),cur(jx+2,jz,jy,2),ax1(jx), &
!!       bx1(jx),cx1(jx),dx1(jx)))
!!    4 continue
!!      endif
!!                  
!!      if(viscous) then     
!!      if(.not.implicitv)then
!!      call vorticity
!!      do 5 jy=1,my
!!      do 5 jz=iz_first,iz_last      
!!      do 5 jx=ix_first+2,ix_last-2
!!      xdif(jx,jz,jy,3)=xdif(jx,jz,jy,3)-fmu0 &
!!       *(voy(jx,jz,jy,3)/xx(jx)-voz(jx,jz,jy,2))
!!     
!!      xdif(jx,jz,jy,4)=xdif(jx,jz,jy,4)-fmu0 &
!!       *(voz(jx,jz,jy,1)-d1fc(vor(jx-2,jz,jy,3),vor(jx-1,jz,jy,3),vor(jx,jz,jy,3), &
!!       vor(jx+1,jz,jy,3),vor(jx+2,jz,jy,3),ax1(jx),bx1(jx),cx1(jx), &
!!       dx1(jx))
!!     
!!      xdif(jx,jz,jy,5)=xdif(jx,jz,jy,5)-fmu0/xx(jx)*(vor(jx,jz,jy,2) &
!!       -voy(jx,jz,jy,1)+xx(jx)*d1fc(vor(jx-2,jz,jy,2),vor(jx-1,jz,jy,2) &
!!       ,vor(jx,jz,jy,2),vor(jx+1,jz,jy,2),vor(jx+2,jz,jy,2),ax1(jx), &
!!       bx1(jx),cx1(jx),dx1(jx)))
!!    5 continue
!!      endif   
!!      endif   
      
!     call bndry
      if(invaryrho) xdif(:,:,:,1)=0
      if(invaryp)   xdif(:,:,:,2)=0
      call mpi_transfersm(xdif,8)


!      if(mod(nstep,nrcd).eq.0) then
!      call recrd_xdif
!      call recrd_cv
!      endif
!      if(nstep.eq.0) call recrd_xdif
      return
      end

!ws***************************************************************************************

     subroutine convt
      USE DECLARE
      include 'mpif.h'
      real*8, dimension(my) :: wwy 
!
!  define statement functions
!  d1f2= d f / dx  with second-order accuracy central difference
      d1f2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(fp1-f0) &
         -(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)
!  d1f2m= d f / dx  with second-order one-sided  difference involving -2 -1 and 0
!  points
      d1f2m(fm2,fm1,f0,xm2,xm1,x0)= &
        ((xm2-x0)/(xm1-x0)*(fm1-f0) &
         -(xm1-x0)/(xm2-x0)*(fm2-f0))/(xm2-xm1)
!  d1f2p= d f / dx  with second-order one-sided  difference involving 2 1 and 0
!  points
      d1f2p(fp2,fp1,f0,xp2,xp1,x0)= &
        ((xp2-x0)/(xp1-x0)*(fp1-f0) &
         -(xp1-x0)/(xp2-x0)*(fp2-f0))/(xp2-xp1)
!  d1fc= d f / dx  with fourth-order accuracy central difference
      d1fc(fm2,fm1,f0,fp1,fp2,a,b,c,d)= &
       a*(fp1-f0)+b*(f0-fm1)+c*(fp2-f0)+d*(f0-fm2)
!  d2fc= d2 f / dx2   with third-order accuracy central difference
      d2fc(fm2,fm1,f0,fp1,fp2,a,b,c,d)= &
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
      d1fbp(fm1,f0,fp1,fp2,a,b,c)= &
       a*(fm1-f0)+b*(fp1-f0)+c*(fp2-f0)
!  d1fbm= d f / dx  with  one-sided-bias  difference involving -2 -1  0 and 1
!  points
      d1fbm(fm2,fm1,f0,fp1,a,b,c)= &
       a*(f0-fp1)+b*(f0-fm1)+c*(f0-fm2)
!  d2f2= d2f / dx2  with second-order accuracy central difference
      d2f2(fm1,f0,fp1,xm1,x0,xp1)= &
       2*((fp1-f0)/(xp1-x0)-(f0-fm1)/(x0-xm1))/(xp1-xm1)
!  d2fbp= d2 f / dx2  with  one-sided-bias  difference involving -1 0  1 and 2
!  points
      d2fbp(fm1,f0,fp1,fp2,a,b,c)= &
       2*(a*(fm1-f0)+b*(fp1-f0)+c*(fp2-f0))
!  d2fbm= d2 f / dx2  with  one-sided-bias  difference involving -2 -1  0 and 1
!  points
      d2fbm(fm2,fm1,f0,fp1,a,b,c)= &
       2*(a*(fp1-f0)+b*(fm1-f0)+c*(fm2-f0))

       integer status(mpi_status_size)
!      init1=1
!      init2=1
!      jxs = ix_first
!      jxe = ix_last
!      jzs = iz_first
!      jze = iz_last
!      if (nrank == 0)        jxs=ix_first
!      if (nrank == nsize - 1) jxe =ix_last

 !      call mpi_transfer8(x1)
      if(spectral) then
      do 18 jx=ix_first,ix_last
      do 18 jz=iz_first,iz_last
      do 17 m=1,8
      wwy(:)=x1(jx,jz,:,m)
      call mpi_transfersy1(wwy,data0)

!      do 2 jy=1,my
!        data0(jy)=x1(jx,jz,jy,m)
!    2 continue
!!    
!!   76 CALL DRCFT2(init1,data,my+2,spec,my/2+1,my,mz,-1,scale,AUX11,naux1, &
!!                AUX12,naux2)
!!      if (init1.eq.1) then
!!       init1 = init1-1
!!       goto 76

!!      endif
!!
      call dfftw_plan_dft_r2c_1d(plan,myt,data0,spec,FFTW_ESTIMATE)
	  call dfftw_execute_dft_r2c(plan,data0,spec)
	  call dfftw_destroy_plan(plan)

      do 3 jy=1,myt/2+1
      spec1(jy)=spec(jy)*fac
    3 continue
      
      if(filt) then

      do 5 jy=1,myt/2+1
!      if(m.eq.2) then
!      spec1(jy)=spec1(jy) 
!      else
      spec1(jy)=spec1(jy)*fmode(jy)
!      endif
      spec(jy)=spec1(jy) 
    5 continue  
!      spec(myt/2+1)=c0

!    77 CALL DCRFT2(init2,spec,my/2+1,data,my+2,my,mz,1,scale,AUX21,naux1,  &
!               AUX22,naux2)
!      if (init2.eq.1) then
!       init2 = init2-1
!       goto 77
!      endif

      call dfftw_plan_dft_c2r_1d(plan,myt,spec,data0,FFTW_ESTIMATE)
	  call dfftw_execute_dft_c2r(plan,spec,data0)
	  call dfftw_destroy_plan(plan)

      do 7 jy=iy_first,iy_last
      if(nrky(nrank).eq.0 .and. jy.le.2) then
      x1(jx,jz,jy,m)=data0(myt+jy-2)
      else if(nrky(nrank).eq.npry-1 .and. jy.ge.my-1) then
      x1(jx,jz,jy,m)=data0(jy-2-mym)
      else
      x1(jx,jz,jy,m)=data0(nrky(nrank)*mym+jy-2)
      endif
    7 continue 

      endif

      do ky=1,myt/2+1
      xkc(jx,jz,ky,m)=real(spec1(ky))
      xks(jx,jz,ky,m)=-imag(spec1(ky))
      enddo


      do 10 i=1,myt/2+1
      spec(i)=c1*(i-1)*spec1(i)
   10 continue
!      spec(myt/2+1)=c0
!      
!   78 CALL DCRFT2(init2,spec,my/2+1,data,my+2,my,mz,1,scale,AUX21,naux1,  &
!               AUX22,naux2)
!      if (init2.eq.1) then
!       init2 = init2-1
!       goto 78
!      endif
!
      call dfftw_plan_dft_c2r_1d(plan,myt,spec,data0,FFTW_ESTIMATE)
	  call dfftw_execute_dft_c2r(plan,spec,data0)
	  call dfftw_destroy_plan(plan)

      do 12 jy=iy_first+2,iy_last-2
      xy(jx,jz,jy,m)=data0(nrky(nrank)*mym+jy-2)
   12 continue
   
   do 13 i=1,myt/2+1
      spec(i)=-(i-1)**2*spec1(i)
   13 continue
!      spec(myt/2+1)=c0
!      
!   78 CALL DCRFT2(init2,spec,my/2+1,data,my+2,my,mz,1,scale,AUX21,naux1,  &
!               AUX22,naux2)
!      if (init2.eq.1) then
!       init2 = init2-1
!       goto 78
!      endif
!
      call dfftw_plan_dft_c2r_1d(plan,myt,spec,data0,FFTW_ESTIMATE)
	  call dfftw_execute_dft_c2r(plan,spec,data0)
	  call dfftw_destroy_plan(plan)

      do 14 jy=iy_first+2,iy_last-2
      xy2(jx,jz,jy,m)=data0(nrky(nrank)*mym+jy-2)
   14 continue     
!      
   17 continue
     
   18 continue

      else !!not spectral
      do 19 m=1,8
      do 15 jz=iz_first,iz_last
      do 15 jx=ix_first,ix_last
!      do jy=2,my-1
!      xy(jx,jz,jy,m)=d1f2(x1(jx,jz,jy-1,m),x1(jx,jz,jy,m),x1(jx,jz,jy+1,m),yy(jy-1),yy(jy),yy(jy+1))      
!      xy2(jx,jz,jy,m)=d2f2(x1(jx,jz,jy-1,m),x1(jx,jz,jy,m),x1(jx,jz,jy+1,m),yy(jy-1),yy(jy),yy(jy+1))
!      enddo
!      xy(jx,jz,1,m)=d1f2(x1(jx,jz,my,m),x1(jx,jz,1,m),x1(jx,jz,2,m),yy(0),yy(1),yy(2))
!      xy(jx,jz,my,m)=d1f2(x1(jx,jz,my-1,m),x1(jx,jz,my,m),x1(jx,jz,1,m),yy(my-1),yy(my),yy(my+1))
!      xy2(jx,jz,1,m)=d2f2(x1(jx,jz,my,m),x1(jx,jz,1,m),x1(jx,jz,2,m),yy(0),yy(1),yy(2))
!      xy2(jx,jz,my,m)=d2f2(x1(jx,jz,my-1,m),x1(jx,jz,my,m),x1(jx,jz,1,m),yy(my-1),yy(my),yy(my+1))
      do jy=iy_first+2,iy_last-2
      xy(jx,jz,jy,m) =d1fc(x1(jx,jz,jy-2,m),x1(jx,jz,jy-1,m),x1(jx,jz,jy,m),x1(jx,jz,jy+1,m),x1(jx,jz,jy+2,m),ay1(jy),by1(jy),cy1(jy),dy1(jy))      
      xy2(jx,jz,jy,m)=d2fc(x1(jx,jz,jy-2,m),x1(jx,jz,jy-1,m),x1(jx,jz,jy,m),x1(jx,jz,jy+1,m),x1(jx,jz,jy+2,m),ay2(jy),by2(jy),cy2(jy),dy2(jy))
!      xy(jx,jz,jy,m) =d1f2(x1(jx,jz,jy-1,m),x1(jx,jz,jy,m),x1(jx,jz,jy+1,m),yy(jy-1),yy(jy),yy(jy+1))
!      xy2(jx,jz,jy,m)=d2f2(x1(jx,jz,jy-1,m),x1(jx,jz,jy,m),x1(jx,jz,jy+1,m),yy(jy-1),yy(jy),yy(jy+1))
      enddo

   15 continue
   19 continue 
      endif !spectral

      do 30 m=1,8
      do 30 jy=iy_first,iy_last
      do 21 jz=iz_first+2,iz_last-2
      do 21 jx=ix_first+2,ix_last-2
      !x1r(jx,jz,jy,m)=d1f2(x1(jx-1,jz,jy,m),x1(jx,jz,jy,m) &
      !   ,x1(jx+1,jz,jy,m),xx(jx-1),xx(jx),xx(jx+1))
      !
      !xr2(jx,jz,jy,m)=d2f2(x1(jx-1,jz,jy,m),x1(jx,jz,jy,m) &
      !   ,x1(jx+1,jz,jy,m),xx(jx-1),xx(jx),xx(jx+1))
      !x1z(jx,jz,jy,m)=d1f2(x1(jx,jz-1,jy,m),x1(jx,jz,jy,m) &
      !   ,x1(jx,jz+1,jy,m),zz(jz-1),zz(jz),zz(jz+1))
      !xz2(jx,jz,jy,m)=d2f2(x1(jx,jz-1,jy,m),x1(jx,jz,jy,m) &
      !   ,x1(jx,jz+1,jy,m),zz(jz-1),zz(jz),zz(jz+1))

     x1r(jx,jz,jy,m)=d1fc(x1(jx-2,jz,jy,m),x1(jx-1,jz,jy,m),x1(jx,jz,jy,m) &
         ,x1(jx+1,jz,jy,m),x1(jx+2,jz,jy,m),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
     xr2(jx,jz,jy,m)=d2fc(x1(jx-2,jz,jy,m),x1(jx-1,jz,jy,m),x1(jx,jz,jy,m) &
         ,x1(jx+1,jz,jy,m),x1(jx+2,jz,jy,m),ax2(jx),bx2(jx),cx2(jx),dx2(jx))
     x1z(jx,jz,jy,m)=d1fc(x1(jx,jz-2,jy,m),x1(jx,jz-1,jy,m),x1(jx,jz,jy,m) &
         ,x1(jx,jz+1,jy,m),x1(jx,jz+2,jy,m),az1(jz),bz1(jz),cz1(jz),dz1(jz))
     xz2(jx,jz,jy,m)=d2fc(x1(jx,jz-2,jy,m),x1(jx,jz-1,jy,m),x1(jx,jz,jy,m) &
         ,x1(jx,jz+1,jy,m),x1(jx,jz+2,jy,m),az2(jz),bz2(jz),cz2(jz),dz2(jz))

      xr(jx,jz,jy,m)=xint_dx(jx,jz,m)+x1r(jx,jz,jy,m)  
      xz(jx,jz,jy,m)=xint_dz(jx,jz,m)+x1z(jx,jz,jy,m)


   21 continue
   30 continue

    return
    end
!ws*************************************************************************
      subroutine current(kf)
      USE DECLARE
      include 'mpif.h'
      real*8 dRBy_dx
      real*8, dimension(mx,mz,my) :: RBy
!
!  define statement functions
!  d1f2= d f / dx  with second-order accuracy central difference
      d1f2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(fp1-f0) &
         -(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)
!  d1fc= d f / dx  with fourth-order accuracy central difference
      d1fc(fm2,fm1,f0,fp1,fp2,a,b,c,d)= &
       a*(fp1-f0)+b*(f0-fm1)+c*(fp2-f0)+d*(f0-fm2)
!  d2fc= d2 f / dx2   with third-order accuracy central difference
!      d2fc(fm2,fm1,f0,fp1,fp2,a,b,c,d)= &
!       a*(fp1-f0)+b*(f0-fm1)+c*(fp2-f0)+d*(f0-fm2)
!  d1fp= d f / dx  with  one-sided  difference involving 0  1 2 and 3
!  points
      d1fp(fp3,fp2,fp1,f0,a,b,c)= &
       a*(fp1-f0)+b*(fp2-f0)+c*(fp3-f0)
!  d1fbp= d f / dx  with  one-sided-bias  difference involving -1 0  1 and 2
!  points
      d1fbp(fp2,fp1,f0,fm1,a,b,c)= &
       a*(fm1-f0)+b*(fp1-f0)+c*(fp2-f0)
!  d1fbm= d f / dx  with  one-sided-bias  difference involving -2 -1  0 and 1
!  points
      d1fbm(fm2,fm1,f0,fp1,a,b,c)= &
       a*(f0-fp1)+b*(f0-fm1)+c*(f0-fm2)
!
!  d1xf2= d Rf / dx  with second-order accuracy central difference
      d1xf2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(xp1*fp1-x0*f0) &
         -(xp1-x0)/(xm1-x0)*(xm1*fm1-x0*f0))/(xm1-xp1)

       integer status(mpi_status_size)

!      do 10 m=1,3
      do 10 jy=iy_first,iy_last
      do 10 jz=iz_first,iz_last
      do 10 jx=ix_first,ix_last
      RBy(jx,jz,jy)=xx(jx)*x1(jx,jz,jy,7)
   10 continue

      do 1 jy=iy_first+2,iy_last-2
      do 1 jz=iz_first+2,iz_last-2
      do 1 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.psia1) then
      cur(jx,jz,jy,1)=xy(jx,jz,jy,8)/xx(jx)-x1z(jx,jz,jy,7) 
      cur(jx,jz,jy,2)=x1z(jx,jz,jy,6)-x1r(jx,jz,jy,8)
!      cur(jx,jz,jy,3)=x1r(jx,jz,jy,7)+(x1(jx,jz,jy,7)-xy(jx,jz,jy,6))/xx(jx)
      dRBy_dx=d1fc(RBy(jx-2,jz,jy),RBy(jx-1,jz,jy),RBy(jx,jz,jy) &
         ,RBy(jx+1,jz,jy),RBy(jx+2,jz,jy),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
  !    dRBy_dx=d1xf2(x1(jx-1,jz,jy,7),x1(jx,jz,jy,7),x1(jx+1,jz,jy,7),xx(jx-1),xx(jx),xx(jx+1))
      cur(jx,jz,jy,3)=(dRBy_dx-xy(jx,jz,jy,6))/xx(jx)

!      cur(jx,jz,jy,1)=d1f2(xh(jx,jz-1,jy,2),xh(jx,jz,jy,2),xh(jx,jz+1,jy,2),zz(jz-1),zz(jz),zz(jz+1)) &
!       -xy(jx,jz,jy,8)/xx(jx) 
!      cur(jx,jz,jy,2)=d1f2(xh(jx-1,jz,jy,3),xh(jx,jz,jy,3),xh(jx+1,jz,jy,3),xx(jx-1),xx(jx),xx(jx+1)) &
!       -d1f2(xh(jx,jz-1,jy,1),xh(jx,jz,jy,1),xh(jx,jz+1,jy,1),zz(jz-1),zz(jz),zz(jz+1))
!      cur(jx,jz,jy,3)=-d1f2(xh(jx-1,jz,jy,2),xh(jx,jz,jy,2),xh(jx+1,jz,jy,2),xx(jx-1),xx(jx),xx(jx+1)) &
!       -(xh(jx,jz,jy,2)-xy(jx,jz,jy,6))/xx(jx)
      endif
    1 continue
!      call bndry_cur_ex(lbnd)
      call bndry3_ex(cur,lbnd)
!      call mpi_transfer3(cur) 
!      call convtc
      if(lcd .eq. 2) then 
      call current_driven
      cur(:,:,:,:)=cur(:,:,:,:)+cud(:,:,:,:)
      endif
!      write(*,*)'cur'
      return
      end    
!ws********************************************************************
      subroutine convtc
      USE DECLARE
      real*8, dimension(my) :: wwy
      include 'mpif.h'
!
!  define statement functions
!  d1f2= d f / dx  with second-order accuracy central difference
      d1f2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(fp1-f0) &
         -(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)
!  d1f2m= d f / dx  with second-order one-sided  difference involving -2 -1 and 0
!  points
      d1f2m(fm2,fm1,f0,xm2,xm1,x0)= &
        ((xm2-x0)/(xm1-x0)*(fm1-f0) &
         -(xm1-x0)/(xm2-x0)*(fm2-f0))/(xm2-xm1)
!  d1f2p= d f / dx  with second-order one-sided  difference involving 2 1 and 0
!  points
      d1f2p(fp2,fp1,f0,xp2,xp1,x0)= &
        ((xp2-x0)/(xp1-x0)*(fp1-f0) &
         -(xp1-x0)/(xp2-x0)*(fp2-f0))/(xp2-xp1)

!  d1fc= d f / dx  with fourth-order accuracy central difference
      d1fc(fm2,fm1,f0,fp1,fp2,a,b,c,d)= &
       a*(fp1-f0)+b*(f0-fm1)+c*(fp2-f0)+d*(f0-fm2)
!  d2fc= d2 f / dx2   with third-order accuracy central difference
!      d2fc(fm2,fm1,f0,fp1,fp2,a,b,c,d)= &
!       a*(fp1-f0)+b*(f0-fm1)+c*(fp2-f0)+d*(f0-fm2)
!  d1fp= d f / dx  with  one-sided  difference involving 0  1 2 and 3
!  points
      d1fp(fp3,fp2,fp1,f0,a,b,c)= &
       a*(fp1-f0)+b*(fp2-f0)+c*(fp3-f0)
!  d1fbp= d f / dx  with  one-sided-bias  difference involving -1 0  1 and 2
!  points
      d1fbp(fm1,f0,fp1,fp2,a,b,c)= &
       a*(fm1-f0)+b*(fp1-f0)+c*(fp2-f0)
!  d1fbm= d f / dx  with  one-sided-bias  difference involving -2 -1  0 and 1
!  points
      d1fbm(fm2,fm1,f0,fp1,a,b,c)= &
       a*(f0-fp1)+b*(f0-fm1)+c*(f0-fm2)

       integer status(mpi_status_size)
!      init1=1
!      init2=1
!       jxs = ix_first
!      jxe = ix_last
!      jzs = iz_first
!      jze = iz_last
!      if (nrank == 0)        jxs=ix_first
!      if (nrank == nsize - 1) jxe =ix_last
      if(spectral) then   
      do 18 jz=ix_first,ix_last
      do 18 jx=iz_first,iz_last
      do 17 m=1,3
      wwy(:)=cur(jx,jz,:,m)
      call mpi_transfersy1(wwy,data0)
!      do 1 jy=1,my
!        data0(jy)=cur(jx,jz,jy,m)
!    1 continue
!    
!   76 CALL DRCFT2(init1,data,my+2,spec,my/2+1,my,mz,-1,scale,AUX11,naux1, &
!                AUX12,naux2)
!      if (init1.eq.1) then
!       init1 = init1-1
!       goto 76
!      endif
!
      call dfftw_plan_dft_r2c_1d(plan,myt,data0,spec,FFTW_ESTIMATE)
	  call dfftw_execute_dft_r2c(plan,data0,spec)
	  call dfftw_destroy_plan(plan)

      do 3 i=1,myt/2+1
      spec1(i)=spec(i)*fac !*fmode(i)
!      spec1(i,j)=spec(i,j)*fac*fmode1(jx,i)
    3 continue
!      do 2 i=1,my/2
!      spec1(i,1)=spec(i,1)*fac*fmode(jx,i)
!      spec1(i,2)=spec(i,2)*fac*fmode(jx,i)
!      spec1(i,mz)=spec(i,mz)*fac*fmode(jx,i)
!    2 continue
!
      if(filt) then
      do i=1,myt/2+1
      spec1(i)=spec1(i)*fmode(i)
      spec(i)=spec1(i)
      enddo
!      spec(myt/2+1)=c0

!   77 CALL DCRFT2(init2,spec,my/2+1,data,my+2,my,mz,1,scale,AUX21,naux1,  &
!               AUX22,naux2)
!      if (init2.eq.1) then
!       init2 = init2-1
!       goto 77
!      endif
      
	  call dfftw_plan_dft_c2r_1d(plan,myt,spec,data0,FFTW_ESTIMATE)
	  call dfftw_execute_dft_c2r(plan,spec,data0)
	  call dfftw_destroy_plan(plan)

      do 9 jy=iy_first,iy_last
      if(nrky(nrank).eq.0 .and. jy.le.2) then
      cur(jx,jz,jy,m)=data0(myt+jy-2)
      else if(nrky(nrank).eq.npry-1 .and. jy.ge.my-1) then
      cur(jx,jz,jy,m)=data0(jy-2-mym)
      else
      cur(jx,jz,jy,m)=data0(nrky(nrank)*mym+jy-2)
      endif
    9 continue
      endif
        
      do 10 i=1,myt/2+1
      spec(i)=c1*(i-1)*spec1(i)
   10 continue
!      spec(myt/2+1)=c0

!   78 CALL DCRFT2(init2,spec,my/2+1,data,my+2,my,mz,1,scale,AUX21,naux1,  &
!               AUX22,naux2)
!      if (init2.eq.1) then
!       init2 = init2-1
!       goto 78
!      endif
!
      call dfftw_plan_dft_c2r_1d(plan,myt,spec,data0,FFTW_ESTIMATE)
	  call dfftw_execute_dft_c2r(plan,spec,data0)
	  call dfftw_destroy_plan(plan)
      
      do 14 jy=iy_first+2,iy_last-2
      cuy(jx,jz,jy,m)=data0(nrky(nrank)*mym+jy-2)
   14 continue
!      
   17 continue
      
   18 continue
      
      else !!not spectral
      do 19 m=1,3
      do 15 jz=iz_first,iz_last
      do 15 jx=ix_first,ix_last
!      do jy=2,my-1
!      cuy(jx,jz,jy,m)=d1f2(cur(jx,jz,jy-1,m),cur(jx,jz,jy,m),cur(jx,jz,jy+1,m),yy(jy-1),yy(jy),yy(jy+1))      
!      enddo
!      cuy(jx,jz,1,m)=d1f2(cur(jx,jz,my,m),cur(jx,jz,1,m),cur(jx,jz,2,m),yy(0),yy(1),yy(2))
!      cuy(jx,jz,my,m)=d1f2(cur(jx,jz,my-1,m),cur(jx,jz,my,m),cur(jx,jz,1,m),yy(my-1),yy(my),yy(my+1))
      do jy=iy_first+2,iy_last-2
!      cuy(jx,jz,jy,m)=d1fc(cur(jx,jz,jy-2,m),cur(jx,jz,jy-1,m),cur(jx,jz,jy,m),cur(jx,jz,jy+1,m),cur(jx,jz,jy+2,m),ay1(jy),by1(jy),cy1(jy),dy1(jy)) 
      cuy(jx,jz,jy,m)=d1f2(cur(jx,jz,jy-1,m),cur(jx,jz,jy,m),cur(jx,jz,jy+1,m),yy(jy-1),yy(jy),yy(jy+1))     
      enddo

   15 continue
   19 continue 
      endif !!spectral

      do 30 m=1,3
      do 30 jy=iy_first,iy_last
      do 21 jz=iz_first+2,iz_last-2
      do 21 jx=ix_first+2,ix_last-2
      cux(jx,jz,jy,m)=d1f2(cur(jx-1,jz,jy,m),cur(jx,jz,jy,m),cur(jx+1,jz,jy,m),xx(jx-1),xx(jx),xx(jx+1))
      cuz(jx,jz,jy,m)=d1f2(cur(jx,jz-1,jy,m),cur(jx,jz,jy,m),cur(jx,jz+1,jy,m),zz(jz-1),zz(jz),zz(jz+1))
   21 continue
   30 continue

      return
      end
!ws*************************************
      subroutine current_boot(ibs)
      USE DECLARE
      include 'mpif.h'
      integer ibs
      real*8 bp2,bpp,bb2,bb,pdr,cbs,tbs
      select case(ibs)
!      case(1)
!        do 1 jy=iy_first+2,iy_last-2
!        do 1 jz=iz_first+2,iz_last-2
!        do 1 jx=ix_first+2,ix_last-2
!        if(psi(jx,jz).lt.psia1) then
!        bp2=x(jx,jz,jy,6)**2+x(jx,jz,jy,8)**2
!        bb2=bp2+x(jx,jz,jy,7)**2
!        bpp=sqrt(bp2)
!        bb=sqrt(bb2)
!!        pdr=(1-x(jx,jz,jy,6)/bpp)*xr(jx,jz,jy,2)+(1-x(jx,jz,jy,8)/bpp)*xz(jx,jz,jy,2)
!        pdr=-x(jx,jz,jy,8)/bpp*xr(jx,jz,jy,2)+x(jx,jz,jy,6)/bpp*xz(jx,jz,jy,2)
!        cbs=cbp0(jx,jz)*pdr/bb
!        cub(jx,jz,jy,1)=fbs*(cbs*x(jx,jz,jy,6)-cub0(jx,jz)*xint(jx,jz,6)/bb0(jx,jz))
!        cub(jx,jz,jy,2)=fbs*(cbs*x(jx,jz,jy,7)-cub0(jx,jz)*xint(jx,jz,7)/bb0(jx,jz))
!        cub(jx,jz,jy,3)=fbs*(cbs*x(jx,jz,jy,8)-cub0(jx,jz)*xint(jx,jz,8)/bb0(jx,jz))
!        endif
!   1    continue
      case(1)
        do 1 jy=iy_first+2,iy_last-2
        do 1 jz=iz_first+2,iz_last-2
        do 1 jx=ix_first+2,ix_last-2
        if(psi(jx,jz).lt.psia1) then
        bp2=x(jx,jz,jy,6)**2+x(jx,jz,jy,8)**2
        bb2=bp2+x(jx,jz,jy,7)**2
        bpp=sqrt(bp2)
        bb=sqrt(bb2)
!        pdr=(1-x(jx,jz,jy,6)/bpp)*xr(jx,jz,jy,2)+(1-x(jx,jz,jy,8)/bpp)*xz(jx,jz,jy,2)
!        if(rr(jx,jz).gt. 0.05) then
!        rrthc=(tanh((rr(jx,zj)-rrc)/rrw)-tanh(-rrc/rrw))/tanh((1.0-rrc)/rrw)-tanh(-rrc/rrw)
!        pdr=-x(jx,jz,jy,8)/bpp*xr(jx,jz,jy,2)+x(jx,jz,jy,6)/bpp*xz(jx,jz,jy,2)
        pdr=wx2r(jx,jz)*xr(jx,jz,jy,2)+wz2r(jx,jz)*xz(jx,jz,jy,2)
        cbs=cbp0(jx,jz)*pdr/bb
!        cbs=-(rr(jx,jz)/xx(jx))**0.5*pdr/bpp/bb
        cub(jx,jz,jy,1)=fbs*(cbs*x(jx,jz,jy,6)-cub0(jx,jz)*xint(jx,jz,6)/bb0(jx,jz))
        cub(jx,jz,jy,2)=fbs*(cbs*x(jx,jz,jy,7)-cub0(jx,jz)*xint(jx,jz,7)/bb0(jx,jz))
        cub(jx,jz,jy,3)=fbs*(cbs*x(jx,jz,jy,8)-cub0(jx,jz)*xint(jx,jz,8)/bb0(jx,jz))
!        endif

        endif
   1    continue

      case(10)
        do 10 jy=iy_first+2,iy_last-2
        do 10 jz=iz_first+2,iz_last-2
        do 10 jx=ix_first+2,ix_last-2
        if(psi(jx,jz).lt.psia1) then
        bp2=x(jx,jz,jy,6)**2+x(jx,jz,jy,8)**2
        bb2=bp2+x(jx,jz,jy,7)**2
        bpp=sqrt(bp2)
        bb=sqrt(bb2)
!        pdr=(1-x(jx,jz,jy,6)/bpp)*xr(jx,jz,jy,2)+(1-x(jx,jz,jy,8)/bpp)*xz(jx,jz,jy,2)
        pdr=-x(jx,jz,jy,8)/bpp*xr(jx,jz,jy,2)+x(jx,jz,jy,6)/bpp*xz(jx,jz,jy,2)
        cbs=cbp0(jx,jz)*pdr/bb
        cub(jx,jz,jy,1)=fbs*(cbs*x(jx,jz,jy,6))
        cub(jx,jz,jy,2)=fbs*(cbs*x(jx,jz,jy,7))
        cub(jx,jz,jy,3)=fbs*(cbs*x(jx,jz,jy,8))
        endif
  10    continue

      case(20)
        tbs=0.5*(tanh(20*(time-tbss)/tbsd)+1)
        do 20 jy=iy_first+2,iy_last-2
        do 20 jz=iz_first+2,iz_last-2
        do 20 jx=ix_first+2,ix_last-2
        if(psi(jx,jz).lt.psia1) then
        bp2=x(jx,jz,jy,6)**2+x(jx,jz,jy,8)**2
        bb2=bp2+x(jx,jz,jy,7)**2
        bpp=sqrt(bp2)
        bb=sqrt(bb2)
 !       pdr=(1-x(jx,jz,jy,6)/bpp)*xr(jx,jz,jy,2)+(1-x(jx,jz,jy,8)/bpp)*xz(jx,jz,jy,2)
         pdr=-x(jx,jz,jy,8)/bpp*xr(jx,jz,jy,2)+x(jx,jz,jy,6)/bpp*xz(jx,jz,jy,2)
        cbs=cbp0(jx,jz)*pdr/bb      
        cub(jx,jz,jy,1)=fbs*(cbs*x(jx,jz,jy,6))*tbs
        cub(jx,jz,jy,2)=fbs*(cbs*x(jx,jz,jy,7))*tbs
        cub(jx,jz,jy,3)=fbs*(cbs*x(jx,jz,jy,8))*tbs
        endif
  20    continue
      case(30)
        tbs=0.5*(tanh(20*(time-tbss)/tbsd)+1)
        do 30 jy=iy_first+2,iy_last-2
        do 30 jz=iz_first+2,iz_last-2
        do 30 jx=ix_first+2,ix_last-2
        if(psi(jx,jz).lt.psia1) then
        bp2=x(jx,jz,jy,6)**2+x(jx,jz,jy,8)**2
        bb2=bp2+x(jx,jz,jy,7)**2
        bpp=sqrt(bp2)
        bb=sqrt(bb2)
!        pdr=(1-x(jx,jz,jy,6)/bpp)*xr(jx,jz,jy,2)+(1-x(jx,jz,jy,8)/bpp)*xz(jx,jz,jy,2)
        pdr=-x(jx,jz,jy,8)/bpp*xr(jx,jz,jy,2)+x(jx,jz,jy,6)/bpp*xz(jx,jz,jy,2)
        cbs=cb00*(-(rr(jx,jz)/xx(jx))**0.5*pdr/bpp)
        cub(jx,jz,jy,1)=cbs*x(jx,jz,jy,6)/bb
        cub(jx,jz,jy,2)=cbs*x(jx,jz,jy,7)/bb
        cub(jx,jz,jy,3)=cbs*x(jx,jz,jy,8)/bb
        endif
 30     continue

       case(32)
        tbs=0.5*(tanh(20*(time-tbss)/tbsd)+1)
        do 32 jy=iy_first+2,iy_last-2
        do 32 jz=iz_first+2,iz_last-2
        do 32 jx=ix_first+2,ix_last-2
        if(psi(jx,jz).lt.psia1) then
 
        pdr=-xint(jx,jz,8)/bp0(jx,jz)*xr(jx,jz,jy,2)+xint(jx,jz,6)/bp0(jx,jz)*xz(jx,jz,jy,2)
        cub(jx,jz,jy,1)=0
        cub(jx,jz,jy,2)=cb00*(-(rr(jx,jz)/xx(jx))**0.5*pdr/bpp)
        cub(jx,jz,jy,3)=0
        endif
 32     continue
      end select

      call bndry3_ex(cub,lbnd)
      return
    end
!ws************************************************************************
    !wzhang************************************************************
      subroutine Efield
      USE DECLARE
      include 'mpif.h'


      do 1 jy=iy_first,iy_last
      do 1 jz=iz_first,iz_last
      do 1 jx=ix_first,ix_last
!      Ef(jx,jz,jy,1)=-x(jx,jz,jy,4)*x(jx,jz,jy,8)+x(jx,jz,jy,5)*x(jx,jz,jy,7)+eta(jx,jz,jy)*cur(jx,jz,jy,1)+xint(jx,jz,4)*xint(jx,jz,8)-xint(jx,jz,5)*xint(jx,jz,7)
!      Ef(jx,jz,jy,2)=-x(jx,jz,jy,5)*x(jx,jz,jy,6)+x(jx,jz,jy,3)*x(jx,jz,jy,8)+eta(jx,jz,jy)*cur(jx,jz,jy,2)+xint(jx,jz,5)*xint(jx,jz,6)-xint(jx,jz,3)*xint(jx,jz,8)
!      Ef(jx,jz,jy,3)=-x(jx,jz,jy,3)*x(jx,jz,jy,7)+x(jx,jz,jy,4)*x(jx,jz,jy,6)+eta(jx,jz,jy)*cur(jx,jz,jy,3)+xint(jx,jz,3)*xint(jx,jz,7)-xint(jx,jz,4)*xint(jx,jz,6)

      Ef(jx,jz,jy,1)=-x1(jx,jz,jy,4)*x(jx,jz,jy,8)-xint(jx,jz,4)*x1(jx,jz,jy,8)+x1(jx,jz,jy,5)*x(jx,jz,jy,7)+xint(jx,jz,5)*x1(jx,jz,jy,7)
      Ef(jx,jz,jy,2)=-x1(jx,jz,jy,5)*x(jx,jz,jy,6)-xint(jx,jz,5)*x1(jx,jz,jy,6)+x1(jx,jz,jy,3)*x(jx,jz,jy,8)+xint(jx,jz,3)*x1(jx,jz,jy,8)
      Ef(jx,jz,jy,3)=-x1(jx,jz,jy,3)*x(jx,jz,jy,7)-xint(jx,jz,3)*x1(jx,jz,jy,7)+x1(jx,jz,jy,4)*x(jx,jz,jy,6)+xint(jx,jz,4)*x1(jx,jz,jy,6)
!      Ef(jx,jz,jy,1)=Ef(jx,jz,jy,1)+eta1(jx,jz,jy)*cint(jx,jz,1)
!      Ef(jx,jz,jy,2)=Ef(jx,jz,jy,2)+eta1(jx,jz,jy)*cint(jx,jz,2)
!      Ef(jx,jz,jy,3)=Ef(jx,jz,jy,3)+eta1(jx,jz,jy)*cint(jx,jz,3)

1     continue
      
       if(hall) then
      do 12 jy=iy_first+2,iy_last-2
      do 12 jz=iz_first+2,iz_last-2
      do 12 jx=ix_first+2,ix_last-2

      Ef(jx,jz,jy,1)=Ef(jx,jz,jy,1)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,2)*x(jx,jz,jy,8)-cur(jx,jz,jy,3)*x(jx,jz,jy,7)+cint(jx,jz,2)*x1(jx,jz,jy,8)-cint(jx,jz,3)*x1(jx,jz,jy,7)-x1r(jx,jz,jy,2))
      Ef(jx,jz,jy,2)=Ef(jx,jz,jy,2)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,3)*x(jx,jz,jy,6)-cur(jx,jz,jy,1)*x(jx,jz,jy,8)+cint(jx,jz,3)*x1(jx,jz,jy,6)-cint(jx,jz,1)*x1(jx,jz,jy,8)-xy(jx,jz,jy,2)/xx(jx))    
      Ef(jx,jz,jy,3)=Ef(jx,jz,jy,3)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,1)*x(jx,jz,jy,7)-cur(jx,jz,jy,2)*x(jx,jz,jy,6)+cint(jx,jz,1)*x1(jx,jz,jy,7)-cint(jx,jz,2)*x1(jx,jz,jy,6)-x1z(jx,jz,jy,2))
      
12     continue
       endif
       
!   !turn off the gradp term in hall terms to see the importance of it   
!       if(hall) then
!      do 12 jy=iy_first+2,iy_last-2
!      do 12 jz=iz_first+2,iz_last-2
!      do 12 jx=ix_first+2,ix_last-2
!
!      Ef(jx,jz,jy,1)=Ef(jx,jz,jy,1)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,2)*x(jx,jz,jy,8)-cur(jx,jz,jy,3)*x(jx,jz,jy,7)+cint(jx,jz,2)*x1(jx,jz,jy,8)-cint(jx,jz,3)*x1(jx,jz,jy,7))
!      Ef(jx,jz,jy,2)=Ef(jx,jz,jy,2)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,3)*x(jx,jz,jy,6)-cur(jx,jz,jy,1)*x(jx,jz,jy,8)+cint(jx,jz,3)*x1(jx,jz,jy,6)-cint(jx,jz,1)*x1(jx,jz,jy,8))    
!      Ef(jx,jz,jy,3)=Ef(jx,jz,jy,3)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,1)*x(jx,jz,jy,7)-cur(jx,jz,jy,2)*x(jx,jz,jy,6)+cint(jx,jz,1)*x1(jx,jz,jy,7)-cint(jx,jz,2)*x1(jx,jz,jy,6))
!      
!12     continue
!       endif
   !cd_OXpoint1:vs
      if(cd_OXpoint .and. lcdox==1 .and. mod(nstep,ncd).eq.0 .and. irk.eq.1) call distribution_cd_OXpoint(Ef(:,:,3,2))

      if(resisitive .and. etaJ_in_E) then
      do m=1,3
      Ef(:,:,:,m)=Ef(:,:,:,m)+eta(:,:,:)*cur(:,:,:,m)
      enddo
      endif
   !cd_OXpoint2:Ey_nocd
      if(cd_OXpoint .and. lcdox==2 .and. mod(nstep,ncd).eq.0 .and. irk.eq.1) call distribution_cd_OXpoint(Ef(:,:,3,2))

      if(bootstrap) then
      call current_boot(lbs)      
      do m=1,3
      Ef(:,:,:,m)=Ef(:,:,:,m)-eta(:,:,:)*cub(:,:,:,m)
      enddo
      endif
      
      if(lcd .eq. 1) then
      call current_driven
      do m=1,3
      Ef(:,:,:,m)=Ef(:,:,:,m)-eta(:,:,:)*cud(:,:,:,m)
      enddo
      endif

      if(nstep.lt.nper) then
      do 11 m=1,3
      do 11 jy=iy_first,iy_last
      do 11 jz=iz_first,iz_last
      do 11 jx=ix_first,ix_last
      Ef(jx,jz,jy,m)=Ef(jx,jz,jy,m)+eta1(jx,jz,jy)*(cint(jx,jz,m)+cur(jx,jz,jy,m))
  11  continue
      endif 




!      call Ef_atlastgrid_r1p0_v1(3)

!      if(resisitive) then
!      if(etaJ_in_E) then
!      do 2 jy=1,my
!      do 2 jz=iz_first,iz_last
!      do 2 jx=ix_first,ix_last
!      Ef(jx,jz,jy,1)=Ef(jx,jz,jy,1)+eta(jx,jz,jy)*(cur(jx,jz,jy,1))
!      Ef(jx,jz,jy,2)=Ef(jx,jz,jy,2)+eta(jx,jz,jy)*(cur(jx,jz,jy,2))
!      Ef(jx,jz,jy,3)=Ef(jx,jz,jy,3)+eta(jx,jz,jy)*(cur(jx,jz,jy,3))
!   2  continue
!      endif
!      endif

!      write(*,*) 'Ef'
       
!      call convte
!      call recrd_Ef
       !***************revised**************************************
       call mpi_transfersm(Ef(:,:,:,:),3)
       !**********************************************************
      if(smoothEf) call smthEf_dis_v2(3)

      return
      end
!**************************************
      subroutine convte
      USE DECLARE
      real*8, dimension(my) :: wwy
      include 'mpif.h'
!
!  define statement functions
!  d1f2= d f / dx  with second-order accuracy central difference
      d1f2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(fp1-f0) &
         -(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)
!  d1f2m= d f / dx  with second-order one-sided  difference involving -2 -1 and 0
!  points
      d1f2m(fm2,fm1,f0,xm2,xm1,x0)= &
        ((xm2-x0)/(xm1-x0)*(fm1-f0) &
         -(xm1-x0)/(xm2-x0)*(fm2-f0))/(xm2-xm1)
!  d1f2p= d f / dx  with second-order one-sided  difference involving 2 1 and 0
!  points
      d1f2p(fp2,fp1,f0,xp2,xp1,x0)= &
        ((xp2-x0)/(xp1-x0)*(fp1-f0) &
         -(xp1-x0)/(xp2-x0)*(fp2-f0))/(xp2-xp1)
!  d1fc= d f / dx  with fourth-order accuracy central difference
      d1fc(fm2,fm1,f0,fp1,fp2,a,b,c,d)= &
       a*(fp1-f0)+b*(f0-fm1)+c*(fp2-f0)+d*(f0-fm2)
!  d2fc= d2 f / dx2   with third-order accuracy central difference
!      d2fc(fm2,fm1,f0,fp1,fp2,a,b,c,d)= &
!       a*(fp1-f0)+b*(f0-fm1)+c*(fp2-f0)+d*(f0-fm2)
!  d1fp= d f / dx  with  one-sided  difference involving 0  1 2 and 3
!  points
      d1fp(fp3,fp2,fp1,f0,a,b,c)= &
       a*(fp1-f0)+b*(fp2-f0)+c*(fp3-f0)
!  d1fbp= d f / dx  with  one-sided-bias  difference involving -1 0  1 and 2
!  points
      d1fbp(fm1,f0,fp1,fp2,a,b,c)= &
       a*(fm1-f0)+b*(fp1-f0)+c*(fp2-f0)
!  d1fbm= d f / dx  with  one-sided-bias  difference involving -2 -1  0 and 1
!  points
      d1fbm(fm2,fm1,f0,fp1,a,b,c)= &
       a*(f0-fp1)+b*(f0-fm1)+c*(f0-fm2)

       integer status(mpi_status_size)
!      init1=1
!      init2=1
!      jxs = ix_first
!      jxe = ix_last
!      jzs = iz_first
!      jze = iz_last

      if(spectral) then      
      do 18 jz=ix_first,ix_last
      do 18 jx=iz_first,iz_last
      do 17 m=1,3
      wwy(:)=Ef(jx,jz,:,m)
      call mpi_transfersy1(wwy,data0)      
!      do 2 jy=1,my
!        data0(jy)=Ef(jx,jz,jy,m)
!    2 continue
!    
!   76 CALL DRCFT2(init1,data,my+2,spec,my/2+1,my,mz,-1,scale,AUX11,naux1, &
!                AUX12,naux2)
!      if (init1.eq.1) then
!       init1 = init1-1
!       goto 76
!      endif
!
      call dfftw_plan_dft_r2c_1d(plan,myt,data0,spec,FFTW_ESTIMATE)
	  call dfftw_execute_dft_r2c(plan,data0,spec)
	  call dfftw_destroy_plan(plan)

      do 3 jy=1,myt/2+1
      spec1(jy)=spec(jy)*fac  !
    3 continue
      
      if(filt) then
      do 5 jy=1,myt/2+1
      spec1(jy)=spec1(jy)*fmode(jy) 
      spec(jy)=spec1(jy)
    5 continue  
!      spec(myt/2+1)=c0

!    77 CALL DCRFT2(init2,spec,my/2+1,data,my+2,my,mz,1,scale,AUX21,naux1,  &
!               AUX22,naux2)
!      if (init2.eq.1) then
!       init2 = init2-1
!       goto 77
!      endif

      call dfftw_plan_dft_c2r_1d(plan,myt,spec,data0,FFTW_ESTIMATE)
	  call dfftw_execute_dft_c2r(plan,spec,data0)
	  call dfftw_destroy_plan(plan)

      do 7 jy=iy_first,iy_last
      if(nrky(nrank).eq.0 .and. jy.le.2) then
      Ef(jx,jz,jy,m)=data0(myt+jy-2)
      else if(nrky(nrank).eq.npry-1 .and. jy.ge.my-1) then
      Ef(jx,jz,jy,m)=data0(jy-2-mym)
      else
      Ef(jx,jz,jy,m)=data0(nrky(nrank)*mym+jy-2)
      endif
    7 continue 
      endif
      
      do 10 i=1,myt/2+1
      spec(i)=c1*(i-1)*spec1(i)
   10 continue
!      spec(myt/2+1)=c0
!      
!   78 CALL DCRFT2(init2,spec,my/2+1,data,my+2,my,mz,1,scale,AUX21,naux1,  &
!               AUX22,naux2)
!      if (init2.eq.1) then
!       init2 = init2-1
!       goto 78
!      endif
!
      call dfftw_plan_dft_c2r_1d(plan,myt,spec,data0,FFTW_ESTIMATE)
	  call dfftw_execute_dft_c2r(plan,spec,data0)
	  call dfftw_destroy_plan(plan)

      do 12 jy=iy_first+2,iy_last-2
      Efy(jx,jz,jy,m)=data0(nrky(nrank)*mym+jy-2)
   12 continue   
!      
   17 continue
     
   18 continue
      else !!not spectral
      do 19 m=1,3
      do 15 jz=iz_first,iz_last
      do 15 jx=ix_first,ix_last
      do jy=iy_first+2,iy_last-2
      Efy(jx,jz,jy,m)=d1f2(Ef(jx,jz,jy-1,m),Ef(jx,jz,jy,m),Ef(jx,jz,jy+1,m),yy(jy-1),yy(jy),yy(jy+1))  
!      Efy(jx,jz,jy,m)=d1fc(Ef(jx,jz,jy-2,m),Ef(jx,jz,jy-1,m),Ef(jx,jz,jy,m) &
!          ,Ef(jx,jz,jy+1,m),Ef(jx,jz,jy+2,m),ay1(jy),by1(jy),cy1(jy),dy1(jy))              
      enddo
   15 continue
   19 continue 
      endif !!spectral

      do 30 m=1,3
      do 30 jy=iy_first,iy_last
      do 21 jz=iz_first+2,iz_last-2
      do 21 jx=ix_first+2,ix_last-2
!      if(rr(jx,jz).le.aa) then
      Efx(jx,jz,jy,m)=d1f2(Ef(jx-1,jz,jy,m),Ef(jx,jz,jy,m) &
         ,Ef(jx+1,jz,jy,m),xx(jx-1),xx(jx),xx(jx+1))
      Efz(jx,jz,jy,m)=d1f2(Ef(jx,jz-1,jy,m),Ef(jx,jz,jy,m) &
         ,Ef(jx,jz+1,jy,m),zz(jz-1),zz(jz),zz(jz+1))

!     Efx(jx,jz,jy,m)=d1fc(Ef(jx-2,jz,jy,m),Ef(jx-1,jz,jy,m),Ef(jx,jz,jy,m) &
!         ,Ef(jx+1,jz,jy,m),Ef(jx+2,jz,jy,m),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
!     Efz(jx,jz,jy,m)=d1fc(Ef(jx,jz-2,jy,m),Ef(jx,jz-1,jy,m),Ef(jx,jz,jy,m) &
!         ,Ef(jx,jz+1,jy,m),Ef(jx,jz+2,jy,m),az1(jz),bz1(jz),cz1(jz),dz1(jz))

   21 continue
   30 continue
!      call recrd_convte
      return
      end
!ws***********************************************************************
      subroutine convtdivb
      USE DECLARE
      real*8, dimension(mx,mz,my) :: d2dvb
!      real*8, dimension(n2th+5,mpsa-nda:mpsa,my) :: dbst
!      real*8 bx_xy,bx_xz,by_yx,by_yz,bz_zx,bz_zy
      character*12 output
      character*3 cn
      character*3 cn1
      include 'mpif.h'
!
!  d1f2= d f / dx  with second-order accuracy central difference
      d1f2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(fp1-f0) &
         -(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)
!  d1f2m= d f / dx  with second-order one-sided  difference involving -2 -1 and 0
!  points
      d1f2m(fm2,fm1,f0,xm2,xm1,x0)= &
        ((xm2-x0)/(xm1-x0)*(fm1-f0) &
         -(xm1-x0)/(xm2-x0)*(fm2-f0))/(xm2-xm1)
!  d1f2p= d f / dx  with second-order one-sided  difference involving 2 1 and 0
!  points
      d1f2p(fp2,fp1,f0,xp2,xp1,x0)= &
        ((xp2-x0)/(xp1-x0)*(fp1-f0) &
         -(xp1-x0)/(xp2-x0)*(fp2-f0))/(xp2-xp1)

!  d2fxz= d2 f / dxdz
      d2fxz(fm1m1,fp1m1,fm1p1,fp1p1,xm1,x0,xp1,zm1,z0,zp1)= &
       ( fp1p1-fm1p1-fp1m1+fm1m1)/((xp1-xm1)*(zp1-zm1))

!  d1fc= d f / dx  with fourth-order accuracy central difference
      d1fc(fm2,fm1,f0,fp1,fp2,a,b,c,d)= &
       a*(fp1-f0)+b*(f0-fm1)+c*(fp2-f0)+d*(f0-fm2)
!  d2fc= d2 f / dx2   with third-order accuracy central difference
!      d2fc(fm2,fm1,f0,fp1,fp2,a,b,c,d)= &
!       a*(fp1-f0)+b*(f0-fm1)+c*(fp2-f0)+d*(f0-fm2)
!  d1fp= d f / dx  with  one-sided  difference involving 0  1 2 and 3
!  points
      d1fp(fp3,fp2,fp1,f0,a,b,c)= &
       a*(fp1-f0)+b*(fp2-f0)+c*(fp3-f0)
!  d1fbp= d f / dx  with  one-sided-bias  difference involving -1 0  1 and 2
!  points
      d1fbp(fp2,fp1,f0,fm1,a,b,c)= &
       a*(fm1-f0)+b*(fp1-f0)+c*(fp2-f0)
!  d1fbm= d f / dx  with  one-sided-bias  difference involving -2 -1  0 and 1
!  points
      d1fbm(fm2,fm1,f0,fp1,a,b,c)= &
       a*(f0-fp1)+b*(f0-fm1)+c*(f0-fm2)
!  d1xf2= d Rf / dx  with second-order accuracy central difference
      d1xf2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(xp1*fp1-x0*f0) &
         -(xp1-x0)/(xm1-x0)*(xm1*fm1-x0*f0))/(xm1-xp1)

       integer status(mpi_status_size)
!      init1=1
!      init2=1
!      jxs = ix_first
!      jxe = ix_last
!      jzs = iz_first
!      jze = iz_last
!      if (nrank == 0)        jxs=ix_first
!      if (nrank == nsize - 1) jxe =ix_last

      
      do 1 jy=iy_first+2,iy_last-2
      do 1 jz=iz_first+2,iz_last-2
      do 1 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.psia1) then
      dvb(jx,jz,jy)=(d1xf2(x1(jx-1,jz,jy,6),x1(jx,jz,jy,6),x1(jx+1,jz,jy,6),xx(jx-1),xx(jx),xx(jx+1)) &
                   +xy(jx,jz,jy,7))/xx(jx)+x1z(jx,jz,jy,8)
!      dvb(jx,jz,jy)=x1r(jx,jz,jy,6)+(x1(jx,jz,jy,6)+xy(jx,jz,jy,7))/xx(jx)+x1z(jx,jz,jy,8)
      endif
   1  continue
         
      if(mod(nstep,ndgn).eq.1 .and. irk.eq.1) then
      w0(:,:)=dvb(:,:,1)
      call funmax(w0,dbmax0,dbmin0,x_dbmax,x_dbmin,z_dbmax,z_dbmin,xx,zz,mx,mz)
 !mpi   ----------------------------------------------------------------
      CALL MPI_REDUCE(dbmax0,dbmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0, &
            MPI_COMM_WORLD,IERROR)
      CALL MPI_REDUCE(dbmin0,dbmin,1,MPI_DOUBLE_PRECISION,MPI_MIN,0, &
            MPI_COMM_WORLD,IERROR)
!mpi   ----------------------------------------------------------------
      if(nrank==0) then
!      open(unit=14,file='divb.dat',status='unknown',form='formatted')
      write(14,41)time,dbmin,x_dbmin,z_dbmin,dbmax,x_dbmax,z_dbmax
   41 format(7(1x,e12.5))
      endif
      endif
      
      if(idivb==1) then
      do 81 jy=iy_first+2,iy_last-2
      do 81 jz=iz_first+2,iz_last-2
      do 81 jx=ix_first+2,ix_last-2
!      bx_xz=d2fxz(x1(jx-1,jz-1,jy,6),x1(jx+1,jz-1,jy,6),x1(jx-1,jz+1,jy,6),x1(jx+1,jz+1,jy,6),xx(jx-1),xx(jx),xx(jx+1),zz(jz-1),zz(jz),zz(jz+1))
      bx_xz(jx,jz,jy)=(x1(jx+1,jz+1,jy,6)-x1(jx-1,jz+1,jy,6)-x1(jx+1,jz-1,jy,6)+x1(jx-1,jz-1,jy,6))/((xx(jx+1)-xx(jx-1))*(zz(jz+1)-zz(jz-1)))
      bx_xy(jx,jz,jy)=(x1(jx+1,jz,jy+1,6)-x1(jx-1,jz,jy+1,6)-x1(jx+1,jz,jy-1,6)+x1(jx-1,jz,jy-1,6))/((xx(jx+1)-xx(jx-1))*(yy(jy+1)-yy(jy-1)))
      by_yx(jx,jz,jy)=(x1(jx+1,jz,jy+1,7)-x1(jx-1,jz,jy+1,7)-x1(jx+1,jz,jy-1,7)+x1(jx-1,jz,jy-1,7))/((xx(jx+1)-xx(jx-1))*(yy(jy+1)-yy(jy-1)))
      by_yz(jx,jz,jy)=(x1(jx,jz+1,jy+1,7)-x1(jx,jz-1,jy+1,7)-x1(jx,jz+1,jy-1,7)+x1(jx,jz-1,jy-1,7))/((zz(jz+1)-zz(jz-1))*(yy(jy+1)-yy(jy-1)))
      bz_zx(jx,jz,jy)=(x1(jx+1,jz+1,jy,8)-x1(jx-1,jz+1,jy,8)-x1(jx+1,jz-1,jy,8)+x1(jx-1,jz-1,jy,8))/((xx(jx+1)-xx(jx-1))*(zz(jz+1)-zz(jz-1)))
      bz_zy(jx,jz,jy)=(x1(jx,jz+1,jy+1,8)-x1(jx,jz+1,jy-1,8)-x1(jx,jz-1,jy+1,8)+x1(jx,jz-1,jy-1,8))/((yy(jy+1)-yy(jy-1))*(zz(jz+1)-zz(jz-1)))

      divb_x(jx,jz,jy)=xr2(jx,jz,jy,6)+x1r(jx,jz,jy,6)/xx(jx)-x1(jx,jz,jy,6)/xx(jx)**2 &
                      +by_yx(jx,jz,jy)/xx(jx)-xy(jx,jz,jy,7)/xx(jx)**2+bz_zx(jx,jz,jy)
      divb_y(jx,jz,jy)=bx_xy(jx,jz,jy)+xy(jx,jz,jy,6)/xx(jx)+xy2(jx,jz,jy,7)/xx(jx)+bz_zy(jx,jz,jy)                       
      divb_z(jx,jz,jy)=bx_xz(jx,jz,jy)+x1z(jx,jz,jy,6)/xx(jx)+by_yz(jx,jz,jy)/xx(jx)+xz2(jx,jz,jy,8)
  81  continue

      else !!ws:idivb!=1
      
      call valbm_atlastgrid_v1(dvb(:,:,:),1,1)  
      if(spectral) then  
      do 18 jz=ix_first,ix_last
      do 18 jx=iz_first,iz_last

     call mpi_transfersy1(dvb(jx,jz,:),data0)
!      do 2 jy=1,my
!        data0(jy)=dvb(jx,jz,jy)
!    2 continue

!    
!   76 CALL DRCFT2(init1,data,my+2,spec,my/2+1,my,mz,-1,scale,AUX11,naux1, &
!                AUX12,naux2)
!      if (init1.eq.1) then
!       init1 = init1-1
!       goto 76
!      endif
!
      call dfftw_plan_dft_r2c_1d(plan,myt,data0,spec,FFTW_ESTIMATE)
	  call dfftw_execute_dft_r2c(plan,data0,spec)
	  call dfftw_destroy_plan(plan)

      do 3 jy=1,myt/2
      spec1(jy)=spec(jy)*fac  !*fmode(jy)
    3 continue
!      if(filt) then
!      do 5 jy=1,my/2
!      spec(jy)=spec1(jy)*fmode(jy) 
!    5 continue  
!      spec(my/2+1)=c0
!
!!    77 CALL DCRFT2(init2,spec,my/2+1,data,my+2,my,mz,1,scale,AUX21,naux1,  &
!!               AUX22,naux2)
!!      if (init2.eq.1) then
!!       init2 = init2-1
!!       goto 77
!!      endif
!
!      call dfftw_plan_dft_c2r_1d(plan,my,spec,data0,FFTW_ESTIMATE)
!	  call dfftw_execute_dft_c2r(plan,spec,data0)
!	  call dfftw_destroy_plan(plan)
!
!      do 7 jy=1,my
!      dvb(jx,jz,jy)=data0(jy)
!    7 continue 
!      endif
      
      do 10 i=1,myt/2+1
      spec(i)=c1*(i-1)*spec1(i) !*fmode(i)
   10 continue
!      spec(myt/2+1)=c0
!      
!   78 CALL DCRFT2(init2,spec,my/2+1,data,my+2,my,mz,1,scale,AUX21,naux1,  &
!               AUX22,naux2)
!      if (init2.eq.1) then
!       init2 = init2-1
!       goto 78
!      endif
!
      call dfftw_plan_dft_c2r_1d(plan,myt,spec,data0,FFTW_ESTIMATE)
	  call dfftw_execute_dft_c2r(plan,spec,data0)
	  call dfftw_destroy_plan(plan)

      do 12 jy=iy_first+2,iy_last-2
      divb_y(jx,jz,jy)=data0(nrky(nrank)*mym+jy-2)
   12 continue   
!      
     
   18 continue
      else

      do 15 jz=iz_first,iz_last
      do 15 jx=ix_first,ix_last
      do jy=iy_first+2,iy_last-2
      divb_y(jx,jz,jy)=d1f2(dvb(jx,jz,jy-1),dvb(jx,jz,jy),dvb(jx,jz,jy+1),yy(jy-1),yy(jy),yy(jy+1))      
      enddo
 
   15 continue

      endif

      do 30 jy=iy_first,iy_last
      do 21 jz=iz_first+2,iz_last-2
      do 21 jx=ix_first+2,ix_last-2
!      if(rr(jx,jz).le.aa) then
      divb_x(jx,jz,jy)=d1f2(dvb(jx-1,jz,jy),dvb(jx,jz,jy) &
         ,dvb(jx+1,jz,jy),xx(jx-1),xx(jx),xx(jx+1)) 
      divb_z(jx,jz,jy)=d1f2(dvb(jx,jz-1,jy),dvb(jx,jz,jy) &
         ,dvb(jx,jz+1,jy),zz(jz-1),zz(jz),zz(jz+1))  

   21 continue

   30 continue

      endif !!ws:idivb==1
!
!      do 82 jy=1,my
!      do 82 jz=iz_first+1,iz_last-1
!      do 82 jx=ix_first+1,ix_last-1
!      divb_clean(jx,jz,jy,1)=cdb(jx,jz)*divb_x(jx,jz,jy)+cdbx(jx,jz)*dvb(jx,jz,jy)
!
!      divb_clean(jx,jz,jy,2)=cdb(jx,jz)*divb_y(jx,jz,jy)/xx(jx)
!                           
!      divb_clean(jx,jz,jy,3)=cdb(jx,jz)*divb_z(jx,jz,jy)+cdbz(jx,jz)*dvb(jx,jz,jy)
!  82  continue
!      call valb_atlastgrid(divb_clean(:,:,:,:),xhst(:,:,:,:),3,3,0) 

      if(mod(nstep,nrcd).eq.0) then
      do 83 jy=iy_first+2,iy_last-2
      do 83 jz=iz_first+2,iz_last-2
      do 83 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
       d2dvb(jx,jz,jy)=d1f2(divb_x(jx-1,jz,jy),divb_x(jx,jz,jy),divb_x(jx+1,jz,jy),xx(jx-1),xx(jx),xx(jx+1))+divb_x(jx,jz,jy)/xx(jx) &
                      +d1f2(divb_y(jx,jz,jy-1),divb_y(jx,jz,jy),divb_y(jx,jz,jy+1),yy(jy-1),yy(jy),yy(jy+1))/xx(jx)**2 &
                      +d1f2(divb_z(jx,jz-1,jy),divb_z(jx,jz,jy),divb_z(jx,jz+1,jy),zz(jz-1),zz(jz),zz(jz+1)) 
      endif
  83  continue   
      output='d2dvb'//cn1(nrank)//cn(nst)
      open(unit=92,file=output,status='unknown',form='unformatted')
      write(92) d2dvb(:,:,my/2)
      close(92)
      output='dvb'//cn1(nrank)//cn(nst)
      open(unit=91,file=output,status='unknown',form='unformatted')
      write(91) dvb(:,:,my/2)
      close(91)

      endif

!      call recrd_convte
      return
      end


!ws**********************************************************************
      subroutine diagn
      USE DECLARE
      include 'mpif.h'
      real*8 cxmax,x_cxmax,z_cxmax,cymax,x_cymax,z_cymax,czmax,x_czmax,z_czmax,crmax,x_crmax,z_crmax,cpmax,x_cpmax,z_cpmax
      real*8 cxmin,x_cxmin,z_cxmin,cymin,x_cymin,z_cymin,czmin,x_czmin,z_czmin,crmin,x_crmin,z_crmin,cpmin,x_cpmin,z_cpmin
      real*8 bxmax,x_bxmax,z_bxmax,bymax,x_bymax,z_bymax,bzmax,x_bzmax,z_bzmax,brmax,x_brmax,z_brmax,bpmax,x_bpmax,z_bpmax
      real*8 bxmin,x_bxmin,z_bxmin,bymin,x_bymin,z_bymin,bzmin,x_bzmin,z_bzmin,brmin,x_brmin,z_brmin,bpmin,x_bpmin,z_bpmin

      real*8 cxmax1,x_cxmax1,z_cxmax1,cymax1,x_cymax1,z_cymax1,czmax1,x_czmax1,z_czmax1,crmax1,x_crmax1,z_crmax1,cpmax1,x_cpmax1,z_cpmax1
      real*8 cxmin1,x_cxmin1,z_cxmin1,cymin1,x_cymin1,z_cymin1,czmin1,x_czmin1,z_czmin1,crmin1,x_crmin1,z_crmin1,cpmin1,x_cpmin1,z_cpmin1
      real*8 bxmax1,x_bxmax1,z_bxmax1,bymax1,x_bymax1,z_bymax1,bzmax1,x_bzmax1,z_bzmax1,brmax1,x_brmax1,z_brmax1,bpmax1,x_bpmax1,z_bpmax1
      real*8 bxmin1,x_bxmin1,z_bxmin1,bymin1,x_bymin1,z_bymin1,bzmin1,x_bzmin1,z_bzmin1,brmin1,x_brmin1,z_brmin1,bpmin1,x_bpmin1,z_bpmin1

!
!  d1fc= d f / dx  with fourth-order accuracy central difference
!      d1fc(fm2,fm1,f0,fp1,fp2,a,b,c,d)= &
!       a*(fp1-f0)+b*(f0-fm1)+c*(fp2-f0)+d*(f0-fm2)
!  define statement functions
!  d2fc= d2 f / dx2   with central difference
!      d2fc(fm1,f0,fp1,xm1,x0,xp1)= &
!       2.*((fp1-f0)/(xp1-x0)-(f0-fm1)/(x0-xm1))/(xp1-xm1)
!      call current(0)
!      call convt
      do 20 jz=iz_first,iz_last
      do 20 jx=ix_first,ix_last
        w0(jx,jz)=cur(jx,jz,1,2)
   20 continue
      call funmax(w0,cymax,cymin,x_cymax,x_cymin,z_cymax,z_cymin,xx,zz,mx,mz)

      do 10 jz=iz_first,iz_last
      do 10 jx=ix_first,ix_last
        w0(jx,jz)=cur(jx,jz,1,1)
   10 continue
      call funmax(w0,cxmax,cxmin,x_cxmax,x_cxmin,z_cxmax,z_cxmin,xx,zz,mx,mz)
      do 30 jz=iz_first,iz_last
      do 30 jx=ix_first,ix_last
        w0(jx,jz)=cur(jx,jz,1,3)
   30 continue
      call funmax(w0,czmax,czmin,x_czmax,x_czmin,z_czmax,z_czmin,xx,zz,mx,mz)

      do 40 jz=iz_first,iz_last
      do 40 jx=ix_first,ix_last
        w0(jx,jz)=cur(jx,jz,1,1)*wx2r(jx,jz)+cur(jx,jz,1,3)*wz2r(jx,jz)
   40 continue
      call funmax(w0,crmax,crmin,x_crmax,x_crmin,z_crmax,z_crmin,xx,zz,mx,mz)
      
      do 50 jz=iz_first,iz_last
      do 50 jx=ix_first,ix_last
        w0(jx,jz)=cur(jx,jz,1,1)*wx2p(jx,jz)+cur(jx,jz,1,3)*wz2p(jx,jz)
   50 continue
      call funmax(w0,cpmax,cpmin,x_cpmax,x_cpmin,z_cpmax,z_cpmin,xx,zz,mx,mz)


!      open(unit=13,file='current.dat',status='unknown',form='formatted',position='append')
!      write(13,1000)time,cymax,x_cymax,z_cymax,cxmax,x_cxmax,z_cxmax,czmax,x_czmax,z_czmax,crmax,x_crmax,z_crmax,cpmax,x_cpmax,z_cpmax
!
!1000  format(16(1x,e12.5))

      do 21 jz=iz_first,iz_last
      do 21 jx=ix_first,ix_last
        w0(jx,jz)=x1(jx,jz,1,7)
   21 continue
      call funmax(w0,bymax,bymin,x_bymax,x_bymin,z_bymax,z_bymin,xx,zz,mx,mz)

      do 11 jz=iz_first,iz_last
      do 11 jx=ix_first,ix_last
        w0(jx,jz)=x1(jx,jz,1,6)
   11 continue
      call funmax(w0,bxmax,bxmin,x_bxmax,x_bxmin,z_bxmax,z_bxmin,xx,zz,mx,mz)
      do 31 jz=iz_first,iz_last
      do 31 jx=ix_first,ix_last
        w0(jx,jz)=x1(jx,jz,1,8)
   31 continue
      call funmax(w0,bzmax,bzmin,x_bzmax,x_bzmin,z_bzmax,z_bzmin,xx,zz,mx,mz)

      do 41 jz=iz_first,iz_last
      do 41 jx=ix_first,ix_last
        w0(jx,jz)=x1(jx,jz,1,6)*wx2r(jx,jz)+x1(jx,jz,1,8)*wz2r(jx,jz)
   41 continue
      call funmax(w0,brmax,brmin,x_brmax,x_brmin,z_brmax,z_brmin,xx,zz,mx,mz)
      
      do 51 jz=iz_first,iz_last
      do 51 jx=ix_first,ix_last
        w0(jx,jz)=x1(jx,jz,1,6)*wx2p(jx,jz)+x1(jx,jz,1,8)*wz2p(jx,jz)
   51 continue
      call funmax(w0,bpmax,bpmin,x_bpmax,x_bpmin,z_bpmax,z_bpmin,xx,zz,mx,mz)


!      open(unit=15,file='bfield.dat',status='unknown',form='formatted',position='append')
!      write(15,1000)time,bymax,x_bymax,z_bymax,bxmax,x_bxmax,z_bxmax,bzmax,x_bzmax,z_bzmax,brmax,x_brmax,z_brmax,bpmax,x_bpmax,z_bpmax
!
!mpi   -----------------------------------------------------------------
      CALL MPI_ALLREDUCE(cymax,cymax1,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(x_cymax,x_cymax1,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(z_cymax,z_cymax1,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                       MPI_COMM_WORLD,IERROR)

      CALL MPI_ALLREDUCE(cxmax,cxmax1,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(x_cxmax,x_cxmax1,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(z_cxmax,z_cxmax1,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                       MPI_COMM_WORLD,IERROR)

      CALL MPI_ALLREDUCE(czmax,czmax1,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(x_czmax,x_czmax1,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(z_czmax,z_czmax1,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                       MPI_COMM_WORLD,IERROR)

      CALL MPI_ALLREDUCE(crmax,crmax1,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(x_crmax,x_crmax1,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(z_crmax,z_crmax1,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                       MPI_COMM_WORLD,IERROR)

      CALL MPI_ALLREDUCE(cpmax,cpmax1,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(x_cpmax,x_cpmax1,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(z_cpmax,z_cpmax1,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                       MPI_COMM_WORLD,IERROR)
!bfeild   -----------------------------------------------------------------
      CALL MPI_ALLREDUCE(bymax,bymax1,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(x_bymax,x_bymax1,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(z_bymax,z_bymax1,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                       MPI_COMM_WORLD,IERROR)

      CALL MPI_ALLREDUCE(bxmax,bxmax1,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(x_bxmax,x_bxmax1,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(z_bxmax,z_bxmax1,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                       MPI_COMM_WORLD,IERROR)

      CALL MPI_ALLREDUCE(bzmax,bzmax1,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(x_bzmax,x_bzmax1,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(z_bzmax,z_bzmax1,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                       MPI_COMM_WORLD,IERROR)

      CALL MPI_ALLREDUCE(brmax,brmax1,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(x_brmax,x_brmax1,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(z_brmax,z_brmax1,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                       MPI_COMM_WORLD,IERROR)

      CALL MPI_ALLREDUCE(bpmax,bpmax1,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(x_bpmax,x_bpmax1,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(z_bpmax,z_bpmax1,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                       MPI_COMM_WORLD,IERROR) 

      if(nrank.eq.0) then
!      open(unit=13,file='current.dat',status='unknown',form='formatted',position='append')
      write(13,1000)time,cymax1,x_cymax1,z_cymax1,cxmax1,x_cxmax1,z_cxmax1,czmax1,x_czmax1,z_czmax1,crmax1,x_crmax1,z_crmax1,cpmax1,x_cpmax1,z_cpmax1
!      open(unit=15,file='bfield.dat',status='unknown',form='formatted',position='append')
      write(15,1000)time,bymax1,x_bymax1,z_bymax1,bxmax1,x_bxmax1,z_bxmax1,bzmax1,x_bzmax1,z_bzmax1,brmax1,x_brmax1,z_brmax1,bpmax1,x_bpmax1,z_bpmax1
1000  format(16(1x,e12.5))

!      print*,cxmax1,cymax1,czmax1,cxmin1,cymin1,czmin1,time
      endif


      call diagnatxmode

      return
      end
!ws*************************************************      
      subroutine diagnatxmode
      USE DECLARE
      real*8 vsmode
      include 'mpif.h' 
    
      if(nrank==nrank_mode) then
      vsmode=x(jxmode,jzmode,1,3)*x(jxmode,jzmode,1,8)-x(jxmode,jzmode,1,5)*x(jxmode,jzmode,1,6)
      write(16,1300) xx(jxmode),time,vsmode,(x1(jxmode,jzmode,1,i),i=1,8),(cur(jxmode,jzmode,1,ic),ic=1,3),(Ef(jxmode,jzmode,1,ie),ie=1,3)
1300  format(17(1x,e12.5))
      endif
      return
      end
      
!ws*************************************************      
      subroutine energy
      USE DECLARE
      include 'mpif.h'
!
      real*8, dimension(mx,mz) :: rho_0
      real*8, dimension(mz,my) :: fyz,gyz,hyz !,wyz
      real*8, dimension(my) :: ffz,ggz,hhz,vc,wwy  !,wwz
      real*8, dimension(mx,mz,myt/2+1) :: gc,gs,gn
      real*8, dimension(mz,myt/2+1) :: gcyz,gsyz,gnyz
      real*8, dimension(myt/2+1) :: gcy,gsy,gcyt,gsyt,gny,gnyt
      real*8 a00, b00, c00, d00, gc00, gs00, gn00
!
!  define statement functions
!  d2fc= d2 f / dx2   with central difference
!      d2fc(fm1,f0,fp1,xm1,x0,xp1)= &
!       2.*((fp1-f0)/(xp1-x0)-(f0-fm1)/(x0-xm1))/(xp1-xm1)
!  d1fc= d f / dx  with  central difference
!      d1fc(fm1,f0,fp1,xm1,x0,xp1)= &
!        ((xm1-x0)/(xp1-x0)*(fp1-f0) &
!         -(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)
         
       integer status(mpi_status_size)
!       init1=1
!       init2=1
!
!      do 1 jx=i_first+1,i_last-1
!      xh(jx,:,:,1)=-d1fc(x(jx-1,:,:,2),x(jx,:,:,2) &
!        ,x(jx+1,:,:,2),xx(jx-1),xx(jx),xx(jx+1))+cur(jx,:,:,2) &
!        *x(jx,:,:,8)-cur(jx,:,:,3)*x(jx,:,:,7)
!      xh(jx,:,:,2)=-xy(jx,:,:,2)/xx(jx) &
!       +cur(jx,:,:,3)*x(jx,:,:,6)-cur(jx,:,:,1)*x(jx,:,:,8) 
!      xh(jx,:,:,3)=-tau(jx,:)*xz(jx,:,:,2) &
!       +cur(jx,:,:,1)*x(jx,:,:,7)-cur(jx,:,:,2)*x(jx,:,:,6)
!    1 continue
!
      
      do 2 jy=iy_first,iy_last
      do 2 jz=iz_first,iz_last
      do 2 jx=ix_first,ix_last
 !     w(jx,jz,jy)=x(jx,jz,jy,6)**2*xx(jx)/2.
      xh(jx,jz,jy,1)=0.5*(x(jx,jz,jy,6)**2+x(jx,jz,jy,7)**2 &
                   +x(jx,jz,jy,8)**2)*xx(jx)
      xh(jx,jz,jy,2)=0.5*(x(jx,jz,jy,3)**2+x(jx,jz,jy,4)**2 &
                   +x(jx,jz,jy,5)**2)*x(jx,jz,jy,1)*xx(jx)
      xh(jx,jz,jy,3)=x(jx,jz,jy,2)/(gamma-1.)*xx(jx)
    2 continue
    
      
      do 4 jy=iy_first,iy_last
      do 4 jz=iz_first,iz_last
 !     call integ(w(1,jz,jy),a00,xx,ix_first,ix_last,mx)
      call integ(xh(1,jz,jy,1),b00,xx,3,mx-2,mx)
      call integ(xh(1,jz,jy,2),c00,xx,3,mx-2,mx)
      call integ(xh(1,jz,jy,3),d00,xx,3,mx-2,mx)
 !     wyz(jz,jy)=a00
      fyz(jz,jy)=b00
      gyz(jz,jy)=c00
      hyz(jz,jy)=d00
   4  continue
      do 5 jy=iy_first,iy_last
 !      call integ(wyz(1,jy),a00,zz,iz_first,iz_last,mz)
      call integ(fyz(1,jy),b00,zz,3,mz-2,mz)
      call integ(gyz(1,jy),c00,zz,3,mz-2,mz)
      call integ(hyz(1,jy),d00,zz,3,mz-2,mz)
 !     wwz(jy)=a00
      ffz(jy)=b00
      ggz(jy)=c00
      hhz(jy)=d00
    5 continue
 !     call integ(wwz,wt,yy,1,my,my)
      call integ(ffz,ft,yy,3,my-2,my)
      call integ(ggz,gt,yy,3,my-2,my)
      call integ(hhz,ht,yy,3,my-2,my)

      gn=0
      do 18 jz=iz_first,iz_last
      do 18 jx=ix_first,ix_last

!      do jy=1,my
!      vc(jy)=sqrt(x(jx,jz,jy,3)**2+x(jx,jz,jy,4)**2+x(jx,jz,jy,5)**2)
!      enddo
      do m=3,5
          vc(:)=x(jx,jz,:,m)
          call mpi_transfersy1(vc,data0)
    !        data0(1:myt)=xh(jx,jz,3:my-2,2)
    !    
          call dfftw_plan_dft_r2c_1d(plan,myt,data0,spec,FFTW_ESTIMATE)
	      call dfftw_execute_dft_r2c(plan,data0,spec)
	      call dfftw_destroy_plan(plan)
      
          do ky=1,myt/2+1
          spec1(ky)=spec(ky)*fac
          gc(jx,jz,ky)=real(spec1(ky))*xx(jx)
          gs(jx,jz,ky)=-imag(spec1(ky))*xx(jx)
          gn(jx,jz,ky)=gn(jx,jz,ky)+abs(spec1(ky))**2*xx(jx)
          enddo
      enddo
! add by J.Zhu 2015.09.10 to calc n=1 component of mass density rho_0 
      wwy(:)=x(jx,jz,:,1) 
      call mpi_transfersy1(wwy,data0)
!    
      call dfftw_plan_dft_r2c_1d(plan,myt,data0,spec,FFTW_ESTIMATE)
	  call dfftw_execute_dft_r2c(plan,data0,spec)
	  call dfftw_destroy_plan(plan) 
      
      spec1(1)=spec(1)*fac
      rho_0(jx,jz) = real(spec1(1))*xx(jx)
      
      do ky=1,myt/2+1
          gn(jx,jz,ky)=rho_0(jx,jz)*gn(jx,jz,ky)
      enddo
   18 continue 
      
      do 41 ky=1,myt/2+1
      do 41 jz=iz_first,iz_last
 !     call integ(w(1,jz,jy),a00,xx,ix_first,ix_last,mx)
      call integ(gc(1,jz,ky),gc00,xx,3,mx-2,mx)
      call integ(gs(1,jz,ky),gs00,xx,3,mx-2,mx)
      call integ(gn(1,jz,ky),gn00,xx,3,mx-2,mx)
 !     wyz(jz,jy)=a00
      gcyz(jz,ky)=gc00
      gsyz(jz,ky)=gs00
      gnyz(jz,ky)=gn00
   41 continue
      do 51 ky=1,myt/2+1
      call integ(gcyz(1,ky),gc00,zz,3,mz-2,mz)
      call integ(gsyz(1,ky),gs00,zz,3,mz-2,mz)
      call integ(gnyz(1,ky),gn00,zz,3,mz-2,mz)
      gcy(ky)=gc00
      gsy(ky)=gs00
      gny(ky)=gn00
   51 continue
       
      CALL MPI_REDUCE(gcy,gcyt,myt/2,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_REDUCE(gsy,gsyt,myt/2,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_REDUCE(gny,gnyt,myt/2,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                       MPI_COMM_WORLD,IERROR)
!      do ky=1,myt/2
!      CALL MPI_REDUCE(gcy(ky),gcyt(ky),1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
!                       MPI_COMM_WORLD,IERROR)
!      CALL MPI_REDUCE(gsy(ky),gsyt(ky),1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
!                       MPI_COMM_WORLD,IERROR)
!      enddo

      if(nrank.eq.0) then
      if(nstep==0) then
      open(unit=21,file='energy_kyc.dat',status='unknown',form='formatted') 
      open(unit=22,file='energy_kys.dat',status='unknown',form='formatted')
      open(unit=23,file='energy_kyn.dat',status='unknown',form='formatted')
      else
      open(unit=21,file='energy_kyc.dat',status='unknown',form='formatted',position='append') 
      open(unit=22,file='energy_kys.dat',status='unknown',form='formatted',position='append') 
      open(unit=23,file='energy_kyn.dat',status='unknown',form='formatted',position='append') 
      endif

      write(21,400) time,(gcyt(ky),ky=1,myt/2+1)
      write(22,400) time,(gsyt(ky),ky=1,myt/2+1)
      write(23,400) time,(gnyt(ky),ky=1,myt/2+1)
400   format(<myt/2+2>(1x,e12.5))
      endif

      return
      end
!
      subroutine integ(fin,fout,xl,ms,me,nx)
      include 'mpif.h'
      real*8  fout
      real*8, dimension(nx) :: fin,xl
      fout=0
      do 1 jx=ms+1,me
      fout=fout+(fin(jx-1)+fin(jx))*(xl(jx)-xl(jx-1))/2.
    1 continue
      return
      end

      subroutine funmax(f,fmax,fmin,xlmax,xlmin,zlmax,zlmin,x,z,mx,mz)
      include 'mpif.h'
      real*8 fmax,fmin,xlmax,xlmin,zlmax,zlmin
      real*8, dimension(mx,mz) :: f
      real*8, dimension(mx) :: x
      real*8, dimension(mz) :: z
      fmax=-10000.
      fmin=10000.
      do 10 jz=1,mz
      do 10 jx=1,mx
      if(f(jx,jz).gt.fmax) then
      fmax=f(jx,jz)
      xlmax=x(jx)
      zlmax=z(jz)
      jxmax=jx
      jzmax=jz
      endif
      if(f(jx,jz).lt.fmin) then
      fmin=f(jx,jz)
      xlmin=x(jx)
      zlmin=z(jz)
      jxmin=jx
      jzmin=jz
      endif
   10 continue
      return
      end
!
      subroutine funm(f,fmax,fmin,my)
      real*8 fmax,fmin
      real*8, dimension(my) :: f
      fmax=-1.d5
      fmin=1.d5
      do 10 jy=1,my
      fmax=dmax1(f(jy),fmax)
      fmin=dmin1(f(jy),fmin)
   10 continue
!      if(dabs(fmax).le.1.e-18) fmax=0.
!      if(dabs(fmin).le.1.e-18) fmin=0.
      return
      end
!ws************************************************************
      subroutine gridpnt
      USE DECLARE
      include 'mpif.h'
      real*8, dimension(mxt) :: xxht
      real*8, dimension(mzt) :: zzht
      real*8  dwx,dwz,dxmin,dzmin


!
      dwx=3./4
      dwz=3./4
      width=3.
      xxt(1)=xmin
      xmax1=xmg
      
      mx1=mxt/2-mxt/16
      if(uniformx) then
      dxx=(xmax-xmin)/(mxt-1.)
      do 21 jx=1,mxt-1
      xxt(jx+1)=xxt(jx)+dxx
   21 continue
      dxt(:)=dxx

      else
      dxx=(xmax-xmin)/(mxt-1.)
      xxht(1)=xmin
      do jx=1,mxt-1
      xxht(jx+1)=xxht(jx)+dxx
      enddo

      xxt(1)=xmin
      do jx=2,mxt
      xxt(jx)=xmin+(xmax-xmin)*(sinh((xxht(jx)-xmg)/dwx)-sinh((xmin-xmg)/dwx))/(sinh((xmax-xmg)/dwx)-sinh((xmin-xmg)/dwx))
      dxt(jx-1)=xxt(jx)-xxt(jx-1)
      enddo
      dxt(mxt)=dxt(mxt-1)

!      dxmin=0.0025          
!      dxmax=(2.*xmax1-mx1*dxmin)/float(mx1)
!      dxp=0.5*(dxmax+dxmin)
!      dxm=0.5*(dxmax-dxmin)
!      do 22 jx=1,mx1-1
!      dxt(jx+1)=dxp+dxm*dtanh((width*jx-width*mx1/2.)/mx1)
!      xxt(jx+1)=xxt(jx)+dxt(jx)
!   22 continue
!      dxt(1)=dxt(2)
!      do 23 jx=1,mx1
!      help1(jx)=dxt(mx1-jx+1)
!      help(jx)=xxt(mx1)-xxt(mx1-jx+1)
!   23 continue
!      xmax2=xmax-xmax1
!      mx2=mxt-mx1
!      dxmax=(2.*xmax2-mx2*dxmin)/float(mx2)
!      dxp=0.5*(dxmax+dxmin)
!      dxm=0.5*(dxmax-dxmin)
!      do 24 jx=1,mx2-1
!      dxt(jx)=dxp+dxm*dtanh((width*jx-width*mx2/2.)/mx2)
!      xxt(jx+1)=xxt(jx)+dxt(jx)
!   24 continue
!      delta=dxt(1)-help1(mx1)
!      do 25 jx=1,mx2-1
!      dxt(jx)=dxt(jx)-delta
!      xxt(jx+1)=xxt(jx)+dxt(jx)
!   25 continue
!      do 26 jx=mx1+1,mxt
!      help1(jx)=dxt(jx-mx1)
!      help(jx)=help(mx1)+xxt(jx-mx1)
!   26 continue
!      delta=(xmax-help(mxt))/float(mxt-1)
!      do 27 jx=1,mxt-1
!      dxt(jx)=help1(jx)+delta
!      xxt(jx+1)=xxt(jx)+dxt(jx)
!   27 continue
!      delta=(xmax-xxt(mxt))/float(mxt-1)
!      do 28 jx=1,mxt-1
!      dxt(jx)=dxt(jx)+delta
!      xxt(jx+1)=xxt(jx)+dxt(jx)
!   28 continue
      endif
!
      if(symmetryz) zmin=-zmax
      zzt(1)=zmin
      zmax1=zmg
      mz1=mzt/2-mzt/16
      if(uniformz) then
      dzz=(zmax-zmin)/(mzt-1.)
      do 31 jz=1,mzt-1
      zzt(jz+1)=zzt(jz)+dzz
   31 continue
      dzt(:)=dzz
      else
      dzz=(zmax-zmin)/(mzt-1.)
      zzht(1)=zmin
      do jz=1,mzt-1
      zzht(jz+1)=zzht(jz)+dzz
      enddo

      zzt(1)=zmin
      do jz=2,mzt
      zzt(jz)=zmin+(zmax-zmin)*(sinh((zzht(jz)-zmg)/dwz)-sinh((zmin-zmg)/dwz))/(sinh((zmax-zmg)/dwz)-sinh((zmin-zmg)/dwz))
      dzt(jz-1)=zzt(jz)-zzt(jz-1)
      enddo
      dzt(mzt)=dzt(mzt-1)

!      dzmin=0.0025
!      dzmax=(2.*zmax1-mz1*dzmin)/float(mz1)
!      dzp=0.5*(dzmax+dzmin)
!      dzm=0.5*(dzmax-dzmin)
!      do 32 jz=1,mz1-1
!      dzt(jz+1)=dzp+dzm*dtanh((width*jz-width*mz1/2.)/mz1)
!      zzt(jz+1)=zzt(jz)+dzt(jz)
!   32 continue
!      dzt(1)=dzt(2)
!      do 33 jz=1,mz1
!      help1(jz)=dzt(mz1-jz+1)
!      help(jz)=zzt(mz1)-zzt(mz1-jz+1)
!   33 continue
!      zmax2=zmax-zmax1
!      mz2=mzt-mz1
!      dzmax=(2.*zmax2-mz2*dzmin)/float(mz2)
!      dzp=0.5*(dzmax+dzmin)
!      dzm=0.5*(dzmax-dzmin)
!      do 34 jz=1,mz2-1
!      dzt(jz)=dzp+dzm*dtanh((width*jz-width*mz2/2.)/mz2)
!      zzt(jz+1)=zzt(jz)+dzt(jz)
!   34 continue
!      delta=dzt(1)-help1(mz1)
!      do 35 jz=1,mz2-1
!      dzt(jz)=dzt(jz)-delta
!      zzt(jz+1)=zzt(jz)+dzt(jz)
!   35 continue
!      do 36 jz=mz1+1,mzt
!      help1(jz)=dzt(jz-mz1)
!      help(jz)=help(mz1)+zzt(jz-mz1)
!   36 continue
!      delta=(zmax-help(mzt))/float(mzt-1)
!      do 37 jz=1,mzt-1
!      dzt(jz)=help1(jz)+delta
!      zzt(jz+1)=zzt(jz)+dzt(jz)
!   37 continue
!      delta=(zmax-zzt(mzt))/float(mzt-1)
!      do 38 jz=1,mzt-1
!      dzt(jz)=dzt(jz)+delta
!      zzt(jz+1)=zzt(jz)+dzt(jz)
!   38 continue
      endif

      dyy=2.*pi/myt
      do 11 jy=-1,myt+2
      yyt(jy)=(jy-1)*dyy
      dyt(jy)=dyy
   11 continue
      
      do jz=1,mzt
      do jx=1,mxt
      txzt(jx,jz)=atan2(zzt(jz),xxt(jx)-xmg)
      if(zzt(jz) .lt. 0) txzt(jx,jz)=txzt(jx,jz)+2*pi
      rxzt(jx,jz)=sqrt((xxt(jx)-xmg)**2+zzt(jz)**2)
      enddo  
      enddo

!      do jy=1,my
!      jyp(jy)=jy+1
!      jym(jy)=jy-1
!      jyp2(jy)=jy+2
!      jym2(jy)=jy-2
!      enddo
!      jyp(my)=1
!      jym(1)=my
!      jyp2(my)=2
!      jym2(1)=my-1
!      jyp2(my-1)=1
!      jym2(2)=my

      if(nrank.eq.0) then
      open(unit=201,file='gridxx.dat',status='unknown',form='formatted')
      write(201,200)(xxt(jx),jx=1,mxt)
      close(201)
      open(unit=202,file='gridzz.dat',status='unknown',form='formatted')
      write(202,200)(zzt(jz),jz=1,mzt)
      close(202)
      open(unit=203,file='gridyy.dat',status='unknown',form='formatted')
      write(203,200)(yyt(jy),jy=1,myt)
      close(203)
 200  format(1x,e12.5)
      endif

      return
      end
!ws************************************************************
      subroutine setdt
      USE DECLARE
      include 'mpif.h'
!
      dt1=100.
      
      do 1 jy=iy_first+2,iy_last-2
      do 1 jz=iz_first+2,iz_last-2
      do 1 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.psia1) then
      vx=x(jx,jz,jy,3)+x(jx,jz,jy,6)**2/x(jx,jz,jy,1)
      vy=x(jx,jz,jy,4)+x(jx,jz,jy,7)**2/x(jx,jz,jy,1)
      vz=x(jx,jz,jy,5)+x(jx,jz,jy,8)**2/x(jx,jz,jy,1)
      va2=(x(jx,jz,jy,6)**2+x(jx,jz,jy,7)**2  &
          +x(jx,jz,jy,8)**2)/x(jx,jz,jy,1)
      cs2=gamma*x(jx,jz,jy,2)/x(jx,jz,jy,1)
       vpx=dabs(vx)+sqrt(dabs(cs2+va2))
       vpy=dabs(vy)+sqrt(dabs(cs2+va2))
       vpz=dabs(vz)+sqrt(dabs(cs2+va2))
       dtx=dabs(xx(jx)-xx(jx-1))/(vpx/cfl)
       dtz=dabs(zz(jz)-zz(jz-1))/(vpz/cfl)
       dty=dabs(xx(jx)*(yy(jy)-yy(jy-1)))/(vpy/cfl)
 !      dty=xmin/(0.5*my*vpy)

      dt2=dmin1(dtx,dtz)
      dt3=dmin1(dty,dt2)
      dt1=dmin1(dt1,dt3)
      endif
    1 continue
      return
      end
!ws**************************************************************
     subroutine map_nova
     USE DECLARE
      integer, parameter :: nq = 33
      integer, parameter :: nr = 85 !int( sqrt(ndat/3) )
      integer, parameter :: nw = 39
    !
      integer lcell(nr,nr)
      integer lnext(ndat),lnext12(ndat12),lnext34(ndat34)
      real*8 risq(ndat),risq12(ndat12),risq34(ndat34)
      real*8 aw(5,ndat),aw12(5,ndat12),aw34(5,ndat34)
      real*8 rimax,ximin,zimin,dxi,dzi
      real*8, dimension(mxt,mzt) :: bx_dx,bx_dz,bz_dx,bz_dz,by_dx,by_dz,uy_dx,uy_dz,tht_dx,tht_dz
      real*8, dimension(mxt,mzt) :: pt_dx,pt_dz,rh_dx,rh_dz
!      real*8, dimension(mxt,mzt) :: bx,bxdx,bxdz,bz,bzdx,bzdz
!      real*8, dimension(mxt,mzt) :: bxdx_dx,bxdx_dz,bxdz_dx,bxdz_dz,bzdx_dx,bzdx_dz,bzdz_dx,bzdz_dz
!      integer icell(1,1)
!      integer inext(9)
!      real*8 xin(9),zin(9),qin(9),rinsq(9),ain(5,9)
       real*8 xout,zout,qout,qdxout,qdzout
       integer iam,iap
      include 'mpif.h'
!  d1fc= d f / dx  with fourth-order accuracy central difference
      d1fc(fm2,fm1,f0,fp1,fp2,a,b,c,d)= &
       a*(fp1-f0)+b*(f0-fm1)+c*(fp2-f0)+d*(f0-fm2)   
!-----------------------------------------------------------

     write(*,*) ndat
     call qshep2 ( ndat, xx_NOVA, zz_NOVA, ps_NOVA, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
     write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_NOVA, zz_NOVA, ps_NOVA, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, pst(jx,jz), pst_dx(jx,jz), pst_dz(jx,jz), ier )
 !     write(*,*) ier
!       pst(jx,mzt-jz+1)=pst(jx,jz)
!       pst_dx(jx,mzt-jz+1)=pst_dx(jx,jz)
!       pst_dz(jx,mzt-jz+1)=-pst_dz(jx,jz)
      tpt(jx,jz)=atan2(pst_dz(jx,jz),pst_dx(jx,jz))
      if(zzt(jz).lt.0) tpt(jx,jz)=tpt(jx,jz)+2*pi
      enddo
      enddo

     call qshep2 ( ndat12, xx12_NOVA, zz12_NOVA, th12_NOVA, nq, nw, nr, lcell, lnext12, ximin, zimin, &
        dxi, dzi, rimax, risq12, aw12, ier )
!     write(*,*) ier
      do jz=mzt/2+1,mzt
      do jx=1,mxt
      call qs2grd ( xxt(jx), zzt(jz), ndat12, xx12_NOVA, zz12_NOVA, th12_NOVA, nr, lcell, lnext12, ximin, &
        zimin, dxi, dzi, rimax, risq12, aw12, tht(jx,jz), tht_dx(jx,jz), tht_dz(jx,jz), ier )
      enddo
      enddo  
      
      call qshep2 ( ndat34, xx34_NOVA, zz34_NOVA, th34_NOVA, nq, nw, nr, lcell, lnext34, ximin, zimin, &
        dxi, dzi, rimax, risq34, aw34, ier )
!     write(*,*) ier
      do jz=1,mzt/2
      do jx=1,mxt
      call qs2grd ( xxt(jx), zzt(jz), ndat34, xx34_NOVA, zz34_NOVA, th34_NOVA, nr, lcell, lnext34, ximin, &
        zimin, dxi, dzi, rimax, risq34, aw34, tht(jx,jz), tht_dx(jx,jz), tht_dz(jx,jz), ier )
      enddo
      enddo      
          
!      jz=mzt/2+1
!      do jx=1,mxt
!      pst_dz(jx,jz)=d1fc(pst(jx,jz-2),pst(jx,jz-1),pst(jx,jz) &
!         ,pst(jx,jz+1),pst(jx,jz+2),az1(jz),bz1(jz),cz1(jz),dz1(jz))
!      pst_dz(jx,mzt-jz+1)=-pst_dz(jx,jz)
!      enddo

      if(firstmap) then
!!ws:bx
      call qshep2 ( ndat, xx_NOVA, zz_NOVA, bx_NOVA, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_NOVA, zz_NOVA, bx_NOVA, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, bx(jx,jz), bx_dx(jx,jz), bx_dz(jx,jz), ier )
      enddo
      enddo

!!ws:bz
      call qshep2 ( ndat, xx_NOVA, zz_NOVA, bz_NOVA, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_NOVA, zz_NOVA, bz_NOVA, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, bz(jx,jz), bz_dx(jx,jz), bz_dz(jx,jz), ier )
      enddo
      enddo

!!ws:bxdx
      call qshep2 ( ndat, xx_NOVA, zz_NOVA, bxdx_NOVA, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_NOVA, zz_NOVA, bxdx_NOVA, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, bxdx(jx,jz), qdxout,qdzout, ier )
      enddo
      enddo

!!ws:bxdz
      call qshep2 ( ndat, xx_NOVA, zz_NOVA, bxdz_NOVA, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_NOVA, zz_NOVA, bxdz_NOVA, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, bxdz(jx,jz), qdxout,qdzout, ier )
 !     write(*,*) ier
      enddo
      enddo

!!ws:bzdx
      call qshep2 ( ndat, xx_NOVA, zz_NOVA, bzdx_NOVA, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_NOVA, zz_NOVA, bzdx_NOVA, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, bzdx(jx,jz), qdxout,qdzout, ier )
 !     write(*,*) ier
      enddo
      enddo

!!ws:bzdz
      call qshep2 ( ndat, xx_NOVA, zz_NOVA, bzdz_NOVA, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_NOVA, zz_NOVA, bzdz_NOVA, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, bzdz(jx,jz), qdxout,qdzout, ier )
 !     write(*,*) ier
      enddo
      enddo

!!ws:by
      call qshep2 ( ndat, xx_NOVA, zz_NOVA, by_NOVA, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_NOVA, zz_NOVA, by_NOVA, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, by(jx,jz), by_dx(jx,jz), by_dz(jx,jz), ier )
 !     write(*,*) ier
      enddo
      enddo
!!ws:bydx
      call qshep2 ( ndat, xx_NOVA, zz_NOVA, bydx_NOVA, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_NOVA, zz_NOVA, bydx_NOVA, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, bydx(jx,jz), qdxout,qdzout, ier )
 !     write(*,*) ier
      enddo
      enddo
!!ws:bydz
      call qshep2 ( ndat, xx_NOVA, zz_NOVA, bydz_NOVA, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_NOVA, zz_NOVA, bydz_NOVA, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, bydz(jx,jz), qdxout,qdzout, ier )
 !     write(*,*) ier
      enddo
      enddo

!!ws:pdx
      call qshep2 ( ndat, xx_NOVA, zz_NOVA, pdx_NOVA, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_NOVA, zz_NOVA, pdx_NOVA, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, pdx(jx,jz), qdxout,qdzout, ier )
 !     write(*,*) ier
      enddo
      enddo
!!ws:pdz
      call qshep2 ( ndat, xx_NOVA, zz_NOVA, pdz_NOVA, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_NOVA, zz_NOVA, pdz_NOVA, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, pdz(jx,jz), qdxout,qdzout, ier )
 !     write(*,*) ier
      enddo
      enddo
!!ws:cy
      call qshep2 ( ndat, xx_NOVA, zz_NOVA, cy_NOVA, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_NOVA, zz_NOVA, cy_NOVA, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, cy(jx,jz), cy_dx(jx,jz), cy_dz(jx,jz), ier )
 !     write(*,*) ier
      enddo
      enddo
!!ws:cx
      call qshep2 ( ndat, xx_NOVA, zz_NOVA, cx_NOVA, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_NOVA, zz_NOVA, cx_NOVA, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, cx(jx,jz), cx_dx(jx,jz), cx_dz(jx,jz), ier )
 !     write(*,*) ier
      enddo
      enddo
!!ws:cz
      call qshep2 ( ndat, xx_NOVA, zz_NOVA, cz_NOVA, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_NOVA, zz_NOVA, cz_NOVA, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, cz(jx,jz), cz_dx(jx,jz), cz_dz(jx,jz), ier )
 !     write(*,*) ier
      enddo
      enddo
  
 !!ws:uy
      call qshep2 ( ndat, xx_NOVA, zz_NOVA, uy_NOVA, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_NOVA, zz_NOVA, uy_NOVA, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, uy(jx,jz), uy_dx(jx,jz), uy_dz(jx,jz), ier )
 !     write(*,*) ier
      enddo
      enddo
!!ws:uydx
      call qshep2 ( ndat, xx_NOVA, zz_NOVA, uydx_NOVA, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_NOVA, zz_NOVA, uydx_NOVA, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, uydx(jx,jz), qdxout,qdzout, ier )
 !     write(*,*) ier
      enddo
      enddo
!!ws:uydz
      call qshep2 ( ndat, xx_NOVA, zz_NOVA, uydz_NOVA, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_NOVA, zz_NOVA, uydz_NOVA, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, uydz(jx,jz), qdxout,qdzout, ier )
 !     write(*,*) ier
      enddo
      enddo
 !!ws:pt
      call qshep2 ( ndat, xx_NOVA, zz_NOVA, pt_NOVA, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_NOVA, zz_NOVA, pt_NOVA, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, pt(jx,jz), pt_dx(jx,jz), pt_dz(jx,jz), ier )
 !     write(*,*) ier
      enddo
      enddo
!!ws:ptdx
      call qshep2 ( ndat, xx_NOVA, zz_NOVA, ptdx_NOVA, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_NOVA, zz_NOVA, ptdx_NOVA, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, ptdx(jx,jz), qdxout,qdzout, ier )
 !     write(*,*) ier
      enddo
      enddo
!!ws:ptdz
      call qshep2 ( ndat, xx_NOVA, zz_NOVA, ptdz_NOVA, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_NOVA, zz_NOVA, ptdz_NOVA, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, ptdz(jx,jz), qdxout,qdzout, ier )
 !     write(*,*) ier
      enddo
      enddo 
!!ws:rh
      call qshep2 ( ndat, xx_NOVA, zz_NOVA, rh_NOVA, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_NOVA, zz_NOVA, rh_NOVA, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, rh(jx,jz), rh_dx(jx,jz),rh_dz(jx,jz), ier )
 !     write(*,*) ier
      enddo
      enddo
!!ws:rhdx
      call qshep2 ( ndat, xx_NOVA, zz_NOVA, rhdx_NOVA, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_NOVA, zz_NOVA, rhdx_NOVA, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, rhdx(jx,jz), qdxout,qdzout, ier )
 !     write(*,*) ier
      enddo
      enddo
!!ws:rhdz
      call qshep2 ( ndat, xx_NOVA, zz_NOVA, rhdz_NOVA, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_NOVA, zz_NOVA, rhdz_NOVA, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, rhdz(jx,jz), qdxout,qdzout, ier )
 !     write(*,*) ier
      enddo
      enddo 
      
          
      if(nrank.eq.0) then
      open(unit=99,file='mapall',status='unknown',form='unformatted')
      write(99) bx,bxdx,bxdz,bz,bzdx,bzdz,by,bydx,bydz,cx,cy,cz,pt,ptdx,ptdz,rh,rhdx,rhdz,cx_dx,cx_dz,cz_dx,cz_dz,cy_dx,cy_dz,uy,uydx,uydz
      close(99)
      endif

      else

      open(unit=99,file='mapall',status='unknown',form='unformatted')
      read(99)  bx,bxdx,bxdz,bz,bzdx,bzdz,by,bydx,bydz,cx,cy,cz,pt,ptdx,ptdz,rh,rhdx,rhdz,cx_dx,cx_dz,cz_dx,cz_dz,cy_dx,cy_dz,uy,uydx,uydz
      close(99)
      endif
      psia=psival_NOVA(mpsa)
      weight=1.
      psiam=weight*psia+(1-weight)*psival_NOVA(mpsa-1)

      do jx=1,mxt
      do jz=1,mzt
      bpol(jx,jz)=sqrt(bx(jx,jz)**2+bz(jx,jz)**2)
      rr2t(jx,jz)=(pst(jx,jz)-psmin)/(psia-psmin)
      rrt(jx,jz)=sqrt(rr2t(jx,jz))
      enddo
      enddo
      
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
      psi(jx,jz)   =   pst(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
      psi_dx(jx,jz)=pst_dx(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
      psi_dz(jx,jz)=pst_dz(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
      rr2(jx,jz)   =  rr2t(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
      rr(jx,jz)    =   rrt(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
      thxz(jx,jz)  =   tht(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
      tpxz=atan2(psi_dz(jx,jz),psi_dx(jx,jz))
      if(zz(jz).lt.0) tpxz(jx,jz)=tpxz(jx,jz)+2*pi
      enddo
      enddo

          

    do js=mpsa,mps4,-1
    do jx=1+1,mxt-1 
    if(pst(jx,mzt/2).lt.psival_NOVA(js) .and. pst(jx-1,mzt/2).ge.psival_NOVA(js)) jxsmin(js)=jx
    if(pst(jx,mzt/2).lt.psival_NOVA(js) .and. pst(jx+1,mzt/2).ge.psival_NOVA(js)) jxsmax(js)=jx
    enddo
    do jz=1+1,mzt-1 
    if(pst(mxt/2+1,jz).lt.psival_NOVA(js) .and. pst(mxt/2+1,jz-1).ge.psival_NOVA(js)) jzsmin(js)=jz
    if(pst(mxt/2+1,jz).lt.psival_NOVA(js) .and. pst(mxt/2+1,jz+1).ge.psival_NOVA(js)) jzsmax(js)=jz
    enddo
    enddo

    do jx=1+1,mxt-1 
    if(pst(jx,mzt/2).lt.psiam .and. pst(jx-1,mzt/2).ge.psiam) jxamin=jx
    if(pst(jx,mzt/2).lt.psiam .and. pst(jx+1,mzt/2).ge.psiam) jxamax=jx
    enddo
    do jz=1+1,mzt-1 
    if(pst(mxt/2+1,jz).lt.psiam .and. pst(mxt/2+1,jz-1).ge.psiam) jzamin=jz
    if(pst(mxt/2+1,jz).lt.psiam .and. pst(mxt/2+1,jz+1).ge.psiam) jzamax=jz
    enddo

    do jx=jxamin,jxamax

    jzp=mzt
    do while(pst(jx,jzp).ge.psiam)
    jzp=jzp-1
    enddo
    jzap(jx)=jzp
    
    jzm=1
    do while(pst(jx,jzm).ge.psiam)
    jzm=jzm+1
    enddo
    jzam(jx)=jzm
    
    enddo

    do jz=jzamin,jzamax

    jxp=mxt
    do while(pst(jxp,jz).ge.psiam)
    jxp=jxp-1
    enddo
    jxap(jz)=jxp
    
    jxm=1
    do while(pst(jxm,jz).ge.psiam)
    jxm=jxm+1
    enddo
    jxam(jz)=jxm
    
    enddo


!ws150415:xxam(z),xxap(z),zzam(x),zzap(x)------------
    js=mpsa
    iam=nthe2
    iap=nthe2 
    do jx=jxamin,jxamax         
      do while(xxst(iap,js).lt.xxt(jx))
      iap=iap-1
      enddo
      call interp1d3l(zzst(iap-1,js),zzst(iap,js),zzst(iap+1,js),zzst(iap+2,js), &
                      xxst(iap-1,js),xxst(iap,js),xxst(iap+1,js),xxst(iap+2,js),xxt(jx),zzap(jx))


      do while(xxst(iam,js).lt.xxt(jx))
      iam=iam+1
      enddo
      call interp1d3l(zzst(iam-2,js),zzst(iam-1,js),zzst(iam,js),zzst(iam+1,js), &
                      xxst(iam-2,js),xxst(iam-1,js),xxst(iam,js),xxst(iam+1,js),xxt(jx),zzam(jx))

    enddo


    iam=nthe2
    iap=3
    do jz=mzt/2+1,jzamax

      do while(zzst(iam,js).lt.zzt(jz))
      iam=iam-1
      enddo
      call interp1d3l(xxst(iam-1,js),xxst(iam,js),xxst(iam+1,js),xxst(iam+2,js), &
                      zzst(iam-1,js),zzst(iam,js),zzst(iam+1,js),zzst(iam+2,js),zzt(jz),xxam(jz))
      
      do while(zzst(iap,js).lt.zzt(jz))
      iap=iap+1
      enddo
      call interp1d3l(xxst(iap-2,js),xxst(iap-1,js),xxst(iap,js),xxst(iap+1,js), &
                      zzst(iap-2,js),zzst(iap-1,js),zzst(iap,js),zzst(iap+1,js),zzt(jz),xxap(jz))

    enddo


     iam=nthe2
     iap=n2th+5
     do jz=mzt/2,jzamin,-1
 
      do while(zzst(iam,js).gt.zzt(jz))
      iam=iam+1
      enddo
      call interp1d3l(xxst(iam-2,js),xxst(iam-1,js),xxst(iam,js),xxst(iam+1,js), &
                      zzst(iam-2,js),zzst(iam-1,js),zzst(iam,js),zzst(iam+1,js),zzt(jz),xxam(jz))

      do while(zzst(iap,js).gt.zzt(jz))
      iap=iap-1
      enddo
      call interp1d3l(xxst(iap-1,js),xxst(iap,js),xxst(iap+1,js),xxst(iap+2,js), &
                      zzst(iap-1,js),zzst(iap,js),zzst(iap+1,js),zzst(iap+2,js),zzt(jz),xxap(jz))
    enddo
!ws150415:xxam(z),xxap(z),zzam(x),zzap(x)------------
    


      do jz=1,mzt
      do jx=1,mxt
      if( pst(jx,jz).lt. psival_NOVA(npsi-1)) then
      j=1
      do while(psival_NOVA(j) .lt. pst(jx,jz))
      j=j+1
      enddo
      if(j==1) j=j+2
      if(j==2) j=j+1
      call interp1d3l(q_NOVA(j-2),q_NOVA(j-1),q_NOVA(j),q_NOVA(j+1), &
                    psival_NOVA(j-2),psival_NOVA(j-1),psival_NOVA(j),psival_NOVA(j+1),pst(jx,jz),qsf(jx,jz))
      call interp1d3l(g_NOVA(j-2),g_NOVA(j-1),g_NOVA(j),g_NOVA(j+1), &
                    psival_NOVA(j-2),psival_NOVA(j-1),psival_NOVA(j),psival_NOVA(j+1),pst(jx,jz),g(jx,jz))
      call interp1d3l(p_NOVA(j-2),p_NOVA(j-1),p_NOVA(j),p_NOVA(j+1), &
                    psival_NOVA(j-2),psival_NOVA(j-1),psival_NOVA(j),psival_NOVA(j+1),pst(jx,jz),p(jx,jz))
      call interp1d3l(gp_NOVA(j-2),gp_NOVA(j-1),gp_NOVA(j),gp_NOVA(j+1), &
                    psival_NOVA(j-2),psival_NOVA(j-1),psival_NOVA(j),psival_NOVA(j+1),pst(jx,jz),gp(jx,jz))
      call interp1d3l(pp_NOVA(j-2),pp_NOVA(j-1),pp_NOVA(j),pp_NOVA(j+1), &
                    psival_NOVA(j-2),psival_NOVA(j-1),psival_NOVA(j),psival_NOVA(j+1),pst(jx,jz),pp(jx,jz))

     
      call interp1d3l(omrot_NOVA(j-2),omrot_NOVA(j-1),omrot_NOVA(j),omrot_NOVA(j+1), &
                    psival_NOVA(j-2),psival_NOVA(j-1),psival_NOVA(j),psival_NOVA(j+1),pst(jx,jz),omrot(jx,jz))
      call interp1d3l(omprot_NOVA(j-2),omprot_NOVA(j-1),omprot_NOVA(j),omprot_NOVA(j+1), &
                    psival_NOVA(j-2),psival_NOVA(j-1),psival_NOVA(j),psival_NOVA(j+1),pst(jx,jz),omprot(jx,jz))

      endif

      if( pst(jx,jz).ge. psival_NOVA(npsi-1) .and. pst(jx,jz).lt.psival_NOVA(npsi)) then
      call interp1d3l(q_NOVA(npsi-3),q_NOVA(npsi-2),q_NOVA(npsi-1),q_NOVA(npsi), &
                    psival_NOVA(npsi-3),psival_NOVA(npsi-2),psival_NOVA(npsi-1),psival_NOVA(npsi),pst(jx,jz),qsf(jx,jz))
      call interp1d3l(g_NOVA(npsi-3),g_NOVA(npsi-2),g_NOVA(npsi-1),g_NOVA(npsi), &
                    psival_NOVA(npsi-3),psival_NOVA(npsi-2),psival_NOVA(npsi-1),psival_NOVA(npsi),pst(jx,jz),g(jx,jz))
      call interp1d3l(p_NOVA(npsi-3),p_NOVA(npsi-2),p_NOVA(npsi-1),p_NOVA(npsi), &
                    psival_NOVA(npsi-3),psival_NOVA(npsi-2),psival_NOVA(npsi-1),psival_NOVA(npsi),pst(jx,jz),p(jx,jz))
      call interp1d3l(gp_NOVA(npsi-3),gp_NOVA(npsi-2),gp_NOVA(npsi-1),gp_NOVA(npsi), &
                    psival_NOVA(npsi-3),psival_NOVA(npsi-2),psival_NOVA(npsi-1),psival_NOVA(npsi),pst(jx,jz),gp(jx,jz))
      call interp1d3l(pp_NOVA(npsi-3),pp_NOVA(npsi-2),pp_NOVA(npsi-1),pp_NOVA(npsi), &
                    psival_NOVA(npsi-3),psival_NOVA(npsi-2),psival_NOVA(npsi-1),psival_NOVA(npsi),pst(jx,jz),pp(jx,jz))

      call interp1d3l(omrot_NOVA(npsi-3),omrot_NOVA(npsi-2),omrot_NOVA(npsi-1),omrot_NOVA(npsi), &
                    psival_NOVA(npsi-3),psival_NOVA(npsi-2),psival_NOVA(npsi-1),psival_NOVA(npsi),pst(jx,jz),omrot(jx,jz))
      call interp1d3l(omprot_NOVA(npsi-3),omprot_NOVA(npsi-2),omprot_NOVA(npsi-1),omprot_NOVA(npsi), &
                    psival_NOVA(npsi-3),psival_NOVA(npsi-2),psival_NOVA(npsi-1),psival_NOVA(npsi),pst(jx,jz),omprot(jx,jz))
      endif
            
!      uy(jx,jz)=xxt(jx)*omrot(jx,jz)
!      uydx(jx,jz)=-omprot(jx,jz)*bz(jx,jz)*xxt(jx)**2+omrot(jx,jz)
!      uydz(jx,jz)= omprot(jx,jz)*bx(jx,jz)*xxt(jx)**2
      enddo
      enddo


      
      if(rshear) then
      j=1
      do while(q_NOVA(j) .ge. qmode)
      j=j+1
      enddo
      if(j==1) j=j+2
      if(j==2) j=j+1
      call interp1d3l(psival_NOVA(j-2),psival_NOVA(j-1),psival_NOVA(j),psival_NOVA(j+1), &
                    q_NOVA(j-2),q_NOVA(j-1),q_NOVA(j),q_NOVA(j+1),qmode,ps1mode)
      call interp1d3l(xxst(3,j-2),xxst(3,j-1),xxst(3,j),xxst(3,j+1), &
                    q_NOVA(j-2),q_NOVA(j-1),q_NOVA(j),q_NOVA(j+1),qmode,xx1mode)
      else
      j=1
      endif

      do while(q_NOVA(j) .lt. qmode)
      j=j+1
      enddo
      if(j==1) j=j+2
      if(j==2) j=j+1
      call interp1d3l(psival_NOVA(j-2),psival_NOVA(j-1),psival_NOVA(j),psival_NOVA(j+1), &
                    q_NOVA(j-2),q_NOVA(j-1),q_NOVA(j),q_NOVA(j+1),qmode,psmode)

      call interp1d3l(xxst(3,j-2),xxst(3,j-1),xxst(3,j),xxst(3,j+1), &
                    q_NOVA(j-2),q_NOVA(j-1),q_NOVA(j),q_NOVA(j+1),qmode,xxmode)
      jxtmode=floor((xxmode-xmin)/(xmax-xmin)*(mxt-1))+1
      jztmode=mzt/2
      nrkx_mode=(jxtmode-1)/mxm
      nrkz_mode=(jztmode-1)/mzm
      jxmode=jxtmode-nrkx_mode*mxm+2
      jzmode=jztmode-nrkz_mode*mzm+2
      nrank_mode=nrkz_mode*nprx+nrkx_mode
      rrmode=sqrt((psmode-psmin)/(psia-psmin))
      
      if(nrank==nrank_mode) then
      write(*,*) 'mode:','q=',qmode,'x=',xxmode,'nrk=',nrank_mode,nrkx_mode,nrkz_mode
      write(*,*) 'jxt=',jxtmode,'jx=',jxmode,xxt(jxtmode),xx(jxmode)
      write(*,*) 'jzt=',jztmode,'jz=',jzmode,zzt(jztmode),zz(jzmode)
      endif

    if(nrank.eq.0) then

      open(unit=205,file='pst_qpg.dat',status='unknown',form='formatted')
      write(205,100)(((pst(jx,jz),pst_dx(jx,jz),pst_dz(jx,jz),qsf(jx,jz),p(jx,jz),pp(jx,jz),g(jx,jz),gp(jx,jz)),jx=1,mxt),jz=1,mzt)
 100  format(8(1x,e12.5))
      close(205)

      open(unit=204,file='gridts.dat',status='unknown',form='formatted')
      write(204,200)( (tpt(jx,jz),jx=1,mxt),jz=1,mzt)
      close(204)
      open(unit=206,file='gridth.dat',status='unknown',form='formatted')
      write(206,200)( (tht(jx,jz),jx=1,mxt),jz=1,mzt)
      close(206)
      open(unit=207,file='gridtc.dat',status='unknown',form='formatted')
      write(207,3000)( ((tht(jx,jz),tht_dx(jx,jz),tht_dz(jx,jz)),jx=1,mxt),jz=1,mzt)
      close(207)
  
200   format((1x,e12.5))
3000  format(3(1x,e12.5))      
      open(unit=301,file='bx0.dat',status='unknown',form='formatted')
      write(301,500)(((bx(jx,jz),bx_dx(jx,jz),bxdx(jx,jz),bx_dz(jx,jz),bxdz(jx,jz)),jx=1,mxt),jz=1,mzt)
 500  format(5(1x,e12.5))
      close(301)
      open(unit=302,file='bz0.dat',status='unknown',form='formatted')
      write(302,500)(((bz(jx,jz),bz_dx(jx,jz),bzdx(jx,jz),bz_dz(jx,jz),bzdz(jx,jz)),jx=1,mxt),jz=1,mzt)
      close(302)
      open(unit=303,file='c0.dat',status='unknown',form='formatted')
      write(303,900)(((cx(jx,jz),cx_dx(jx,jz),cx_dz(jx,jz),cy(jx,jz),cy_dx(jx,jz),cy_dz(jx,jz),cz(jx,jz),cz_dx(jx,jz),cz_dz(jx,jz)),jx=1,mxt),jz=1,mzt)
 900  format(9(1x,e12.5))
      close(303)
      open(unit=304,file='pt0.dat',status='unknown',form='formatted')
      write(304,500)(((pt(jx,jz),pt_dx(jx,jz),ptdx(jx,jz),pt_dz(jx,jz),ptdz(jx,jz)),jx=1,mxt),jz=1,mzt)
      close(304)
      open(unit=305,file='rh0.dat',status='unknown',form='formatted')
      write(305,500)(((rh(jx,jz),rh_dx(jx,jz),rhdx(jx,jz),rh_dz(jx,jz),rhdz(jx,jz)),jx=1,mxt),jz=1,mzt)
      close(305)

      open(unit=306,file='uy0.dat',status='unknown',form='formatted')
      write(306,500)(((uy(jx,jz),uy_dx(jx,jz),uydx(jx,jz),uy_dz(jx,jz),uydz(jx,jz)),jx=1,mxt),jz=1,mzt)
      close(306)
     endif
!ws150303
     call estimate_pst1
  !   call readmap_wch
!ws150303
     call last_grid
     return
     end
!ws****************************************************************************
     subroutine mapeq_st2xz(fst,fxz,jjx,jjz,itp,isp,rs)
     USE DECLARE
     integer jjx,jjz,itp,isp,is,ks
     real*8, dimension(n2th+5,npsi):: fst
     real*8, dimension(mxt,mzt):: fxz
     real*8, dimension(4):: rs,fsxz

     do ks=1,4,1
     is=isp+ks-3
     call interp1d3l(fst(itp-2,is),fst(itp-1,is),fst(itp,is),fst(itp+1,is), &
                    tst(itp-2,is),tst(itp-1,is),tst(itp,is),tst(itp+1,is),txzt(jjx,jjz),fsxz(ks) )
     enddo
     call interp1d3l(fsxz(1),fsxz(2),fsxz(3),fsxz(4), &
                    rs(1),rs(2),rs(3),rs(4),rxzt(jjx,jjz),fxz(jjx,jjz) )

     return
     end
!ws**************************************************************
     subroutine map_nova_v2
      USE DECLARE
      integer, dimension(mxt,mzt):: it_inp,is_inp
      integer, dimension(npsi):: itj
      real*8, dimension(npsi):: rsj
      real*8, dimension(4):: rsxz,tsxz
      include 'mpif.h'
!  d1fc= d f / dx  with fourth-order accuracy central difference
      d1fc(fm2,fm1,f0,fp1,fp2,a,b,c,d)= &
       a*(fp1-f0)+b*(f0-fm1)+c*(fp2-f0)+d*(f0-fm2)   
!-----------------------------------------------------------
     do jz=1,mzt
     do jx=1,mxt

     do j=1,npsi
     i=1
     do while(tst(i,j) .lt. txzt(jx,jz))
     i=i+1
     enddo
     call interp1d3l(rst(i-2,j),rst(i-1,j),rst(i,j),rst(i+1,j), &
                    tst(i-2,j),tst(i-1,j),tst(i,j),tst(i+1,j),txzt(jx,jz),rsj(j))
     itj(j)=i
     enddo

     j=1
     do while(rsj(j) .lt. rxzt(jx,jz) .and. j .lt. npsi)
     j=j+1
     enddo
     is_inp(jx,jz)=j
     if(j.le.npsi) it_inp(jx,jz)=itj(j)
     
     enddo
     enddo

     if(firstmap) then
     do jz=1,mzt
     do jx=1,mxt
     if(is_inp(jx,jz).le. npsi) then
     itp=it_inp(jx,jz)
     isp=is_inp(jx,jz)

     if(isp==2) isp=isp+1
     if(isp==npsi) isp=npsi-1
     do ks=1,4,1
     is=isp+ks-3
     call interp1d3l(rst(itp-2,is),rst(itp-1,is),rst(itp,is),rst(itp+1,is), &
                    tst(itp-2,is),tst(itp-1,is),tst(itp,is),tst(itp+1,is),txzt(jx,jz),rsxz(ks) )
     call interp1d3l(thst(itp-2),thst(itp-1),thst(itp),thst(itp+1), &
                    tst(itp-2,is),tst(itp-1,is),tst(itp,is),tst(itp+1,is),txzt(jx,jz),tsxz(ks) )
     enddo
     call interp1d3l(psival_NOVA(isp-2),psival_NOVA(isp-1),psival_NOVA(isp),psival_NOVA(isp+1), &
                    rsxz(1),rsxz(2),rsxz(3),rsxz(4),rxzt(jx,jz),pst(jx,jz) )
     call interp1d3l(tsxz(1),tsxz(2),tsxz(3),tsxz(4), &
                    rsxz(1),rsxz(2),rsxz(3),rsxz(4),rxzt(jx,jz),tht(jx,jz) )

     call mapeq_st2xz(bxst,bx,jx,jz,itp,isp,rsxz)
     call mapeq_st2xz(bxdxst,bxdx,jx,jz,itp,isp,rsxz)
     call mapeq_st2xz(bxdzst,bxdz,jx,jz,itp,isp,rsxz)

     call mapeq_st2xz(bzst,bz,jx,jz,itp,isp,rsxz)
     call mapeq_st2xz(bzdxst,bzdx,jx,jz,itp,isp,rsxz)
     call mapeq_st2xz(bzdzst,bzdz,jx,jz,itp,isp,rsxz)

     call mapeq_st2xz(byst,by,jx,jz,itp,isp,rsxz)
     call mapeq_st2xz(bydxst,bydx,jx,jz,itp,isp,rsxz)
     call mapeq_st2xz(bydzst,bydz,jx,jz,itp,isp,rsxz)

     call mapeq_st2xz(pdxst,pdx,jx,jz,itp,isp,rsxz)
     call mapeq_st2xz(pdzst,pdz,jx,jz,itp,isp,rsxz)

     call mapeq_st2xz(cxst,cx,jx,jz,itp,isp,rsxz)
     call mapeq_st2xz(czst,cz,jx,jz,itp,isp,rsxz)
     call mapeq_st2xz(cyst,cy,jx,jz,itp,isp,rsxz)

     call mapeq_st2xz(uyst,uy,jx,jz,itp,isp,rsxz)
     call mapeq_st2xz(uydxst,uydx,jx,jz,itp,isp,rsxz)
     call mapeq_st2xz(uydzst,uydz,jx,jz,itp,isp,rsxz)
     
     tpt(jx,jz)=atan2(bx(jx,jz),-bz(jx,jz))
     if(zzt(jz).lt.0) tpt(jx,jz)=tpt(jx,jz)+2*pi

     endif

     enddo
     enddo

      
      if(nrank.eq.0) then
      open(unit=98,file='mapst',status='unknown',form='unformatted')
      write(98) pst,tht,tpt
      close(98)

      open(unit=99,file='map',status='unknown',form='unformatted')
      write(99) bx,bxdx,bxdz,bz,bzdx,bzdz,by,bydx,bydz,cx,cy,cz,pdx,pdz,cx_dx,cx_dz,cz_dx,cz_dz,cy_dx,cy_dz,uy,uydx,uydz
      close(99)
      endif

      else
      open(unit=98,file='mapst',status='unknown',form='unformatted')
      read(98) pst,tht,tpt
      close(98)
      open(unit=99,file='map',status='unknown',form='unformatted')
      read(99)  bx,bxdx,bxdz,bz,bzdx,bzdz,by,bydx,bydz,cx,cy,cz,pdx,pdz,cx_dx,cx_dz,cz_dx,cz_dz,cy_dx,cy_dz,uy,uydx,uydz
      close(99)
      endif

      do jx=1,mxt
      do jz=1,mzt
      bpol(jx,jz)=sqrt(bx(jx,jz)**2+bz(jx,jz)**2)
      enddo
      enddo
      
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
      psi(jx,jz)=pst(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
      psi_dx(jx,jz)=-bz(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)*xx(jx)
      psi_dz(jx,jz)=bx(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)*xx(jx)
      rr2(jx,jz)=abs((psi(jx,jz)-psmin)/psmin)
      rr(jx,jz)=sqrt(rr2(jx,jz))
      thxz(jx,jz)=tht(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
      tpxz(jx,jz)=tpt(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
      if(zz(jz).lt.0) tpxz(jx,jz)=tpxz(jx,jz)+2*pi
      enddo
      enddo

      psia=psival_NOVA(mpsa)
      weight=1.
      psiam=weight*psia+(1-weight)*psival_NOVA(mpsa-1)
      

    do js=mpsa,mps4,-1
    do jx=1+1,mxt-1 
    if(pst(jx,mzt/2).lt.psival_NOVA(js) .and. pst(jx-1,mzt/2).ge.psival_NOVA(js)) jxsmin(js)=jx
    if(pst(jx,mzt/2).lt.psival_NOVA(js) .and. pst(jx+1,mzt/2).ge.psival_NOVA(js)) jxsmax(js)=jx
    enddo
    do jz=1+1,mzt-1 
    if(pst(mxt/2+1,jz).lt.psival_NOVA(js) .and. pst(mxt/2+1,jz-1).ge.psival_NOVA(js)) jzsmin(js)=jz
    if(pst(mxt/2+1,jz).lt.psival_NOVA(js) .and. pst(mxt/2+1,jz+1).ge.psival_NOVA(js)) jzsmax(js)=jz
    enddo
    enddo

    do jx=1+1,mxt-1 
    if(pst(jx,mzt/2).lt.psiam .and. pst(jx-1,mzt/2).ge.psiam) jxamin=jx
    if(pst(jx,mzt/2).lt.psiam .and. pst(jx+1,mzt/2).ge.psiam) jxamax=jx
    enddo
    do jz=1+1,mzt-1 
    if(pst(mxt/2+1,jz).lt.psiam .and. pst(mxt/2+1,jz-1).ge.psiam) jzamin=jz
    if(pst(mxt/2+1,jz).lt.psiam .and. pst(mxt/2+1,jz+1).ge.psiam) jzamax=jz
    enddo

    do jx=jxamin,jxamax

    jzp=mzt
    do while(pst(jx,jzp).ge.psiam)
    jzp=jzp-1
    enddo
    jzap(jx)=jzp
    
    jzm=1
    do while(pst(jx,jzm).ge.psiam)
    jzm=jzm+1
    enddo
    jzam(jx)=jzm
    
    enddo

    do jz=jzamin,jzamax

    jxp=mxt
    do while(pst(jxp,jz).ge.psiam)
    jxp=jxp-1
    enddo
    jxap(jz)=jxp
    
    jxm=1
    do while(pst(jxm,jz).ge.psiam)
    jxm=jxm+1
    enddo
    jxam(jz)=jxm
    
    enddo

      do jz=1,mzt
      do jx=1,mxt
      if( pst(jx,jz).lt. psival_NOVA(npsi-1)) then
      j=1
      do while(psival_NOVA(j) .lt. pst(jx,jz))
      j=j+1
      enddo
      if(j==1) j=j+2
      if(j==2) j=j+1
      call interp1d3l(q_NOVA(j-2),q_NOVA(j-1),q_NOVA(j),q_NOVA(j+1), &
                    psival_NOVA(j-2),psival_NOVA(j-1),psival_NOVA(j),psival_NOVA(j+1),pst(jx,jz),qsf(jx,jz))
      call interp1d3l(g_NOVA(j-2),g_NOVA(j-1),g_NOVA(j),g_NOVA(j+1), &
                    psival_NOVA(j-2),psival_NOVA(j-1),psival_NOVA(j),psival_NOVA(j+1),pst(jx,jz),g(jx,jz))
      call interp1d3l(p_NOVA(j-2),p_NOVA(j-1),p_NOVA(j),p_NOVA(j+1), &
                    psival_NOVA(j-2),psival_NOVA(j-1),psival_NOVA(j),psival_NOVA(j+1),pst(jx,jz),p(jx,jz))
      call interp1d3l(gp_NOVA(j-2),gp_NOVA(j-1),gp_NOVA(j),gp_NOVA(j+1), &
                    psival_NOVA(j-2),psival_NOVA(j-1),psival_NOVA(j),psival_NOVA(j+1),pst(jx,jz),gp(jx,jz))
      call interp1d3l(pp_NOVA(j-2),pp_NOVA(j-1),pp_NOVA(j),pp_NOVA(j+1), &
                    psival_NOVA(j-2),psival_NOVA(j-1),psival_NOVA(j),psival_NOVA(j+1),pst(jx,jz),pp(jx,jz))

     
      call interp1d3l(omrot_NOVA(j-2),omrot_NOVA(j-1),omrot_NOVA(j),omrot_NOVA(j+1), &
                    psival_NOVA(j-2),psival_NOVA(j-1),psival_NOVA(j),psival_NOVA(j+1),pst(jx,jz),omrot(jx,jz))
      call interp1d3l(omprot_NOVA(j-2),omprot_NOVA(j-1),omprot_NOVA(j),omprot_NOVA(j+1), &
                    psival_NOVA(j-2),psival_NOVA(j-1),psival_NOVA(j),psival_NOVA(j+1),pst(jx,jz),omprot(jx,jz))

      endif

      if( pst(jx,jz).ge. psival_NOVA(npsi-1) .and. pst(jx,jz).lt.psival_NOVA(npsi)) then
      call interp1d3l(q_NOVA(npsi-3),q_NOVA(npsi-2),q_NOVA(npsi-1),q_NOVA(npsi), &
                    psival_NOVA(npsi-3),psival_NOVA(npsi-2),psival_NOVA(npsi-1),psival_NOVA(npsi),pst(jx,jz),qsf(jx,jz))
      call interp1d3l(g_NOVA(npsi-3),g_NOVA(npsi-2),g_NOVA(npsi-1),g_NOVA(npsi), &
                    psival_NOVA(npsi-3),psival_NOVA(npsi-2),psival_NOVA(npsi-1),psival_NOVA(npsi),pst(jx,jz),g(jx,jz))
      call interp1d3l(p_NOVA(npsi-3),p_NOVA(npsi-2),p_NOVA(npsi-1),p_NOVA(npsi), &
                    psival_NOVA(npsi-3),psival_NOVA(npsi-2),psival_NOVA(npsi-1),psival_NOVA(npsi),pst(jx,jz),p(jx,jz))
      call interp1d3l(gp_NOVA(npsi-3),gp_NOVA(npsi-2),gp_NOVA(npsi-1),gp_NOVA(npsi), &
                    psival_NOVA(npsi-3),psival_NOVA(npsi-2),psival_NOVA(npsi-1),psival_NOVA(npsi),pst(jx,jz),gp(jx,jz))
      call interp1d3l(pp_NOVA(npsi-3),pp_NOVA(npsi-2),pp_NOVA(npsi-1),pp_NOVA(npsi), &
                    psival_NOVA(npsi-3),psival_NOVA(npsi-2),psival_NOVA(npsi-1),psival_NOVA(npsi),pst(jx,jz),pp(jx,jz))

      call interp1d3l(omrot_NOVA(npsi-3),omrot_NOVA(npsi-2),omrot_NOVA(npsi-1),omrot_NOVA(npsi), &
                    psival_NOVA(npsi-3),psival_NOVA(npsi-2),psival_NOVA(npsi-1),psival_NOVA(npsi),pst(jx,jz),omrot(jx,jz))
      call interp1d3l(omprot_NOVA(npsi-3),omprot_NOVA(npsi-2),omprot_NOVA(npsi-1),omprot_NOVA(npsi), &
                    psival_NOVA(npsi-3),psival_NOVA(npsi-2),psival_NOVA(npsi-1),psival_NOVA(npsi),pst(jx,jz),omprot(jx,jz))
      endif
            
!      uy(jx,jz)=xxt(jx)*omrot(jx,jz)
!      uydx(jx,jz)=-omprot(jx,jz)*bz(jx,jz)*xxt(jx)**2+omrot(jx,jz)
!      uydz(jx,jz)= omprot(jx,jz)*bx(jx,jz)*xxt(jx)**2
      enddo
      enddo


      
      if(rshear) then
      j=1
      do while(q_NOVA(j) .ge. qmode)
      j=j+1
      enddo
      if(j==1) j=j+2
      if(j==2) j=j+1
      call interp1d3l(psival_NOVA(j-2),psival_NOVA(j-1),psival_NOVA(j),psival_NOVA(j+1), &
                    q_NOVA(j-2),q_NOVA(j-1),q_NOVA(j),q_NOVA(j+1),qmode,ps1mode)
      call interp1d3l(xxst(3,j-2),xxst(3,j-1),xxst(3,j),xxst(3,j+1), &
                    q_NOVA(j-2),q_NOVA(j-1),q_NOVA(j),q_NOVA(j+1),qmode,xx1mode)
      else
      j=1
      endif

      do while(q_NOVA(j) .lt. qmode)
      j=j+1
      enddo
      if(j==1) j=j+2
      if(j==2) j=j+1
      call interp1d3l(psival_NOVA(j-2),psival_NOVA(j-1),psival_NOVA(j),psival_NOVA(j+1), &
                    q_NOVA(j-2),q_NOVA(j-1),q_NOVA(j),q_NOVA(j+1),qmode,psmode)

      call interp1d3l(xxst(3,j-2),xxst(3,j-1),xxst(3,j),xxst(3,j+1), &
                    q_NOVA(j-2),q_NOVA(j-1),q_NOVA(j),q_NOVA(j+1),qmode,xxmode)
      jxtmode=floor((xxmode-xmin)/(xmax-xmin)*(mxt-1))+1
      jztmode=mzt/2
      nrkx_mode=(jxtmode-1)/mxm
      nrkz_mode=(jztmode-1)/mzm
      jxmode=jxtmode-nrkx_mode*mxm+2
      jzmode=jztmode-nrkz_mode*mzm+2
      nrank_mode=nrkz_mode*nprx+nrkx_mode
      if(nrank==nrank_mode) then
      write(*,*) 'mode:','q=',qmode,'x=',xxmode,'nrk=',nrank_mode,nrkx_mode,nrkz_mode
      write(*,*) 'jxt=',jxtmode,'jx=',jxmode,xxt(jxtmode),xx(jxmode)
      write(*,*) 'jzt=',jztmode,'jz=',jzmode,zzt(jztmode),zz(jzmode)
      endif

    if(nrank.eq.0) then


      open(unit=207,file='istp.dat',status='unknown',form='formatted')
      write(207,2000)(((is_inp(jx,jz),it_inp(jx,jz)),jx=1,mxt),jz=1,mzt)
2000  format(2(1x,i5))
      close(207)

      open(unit=205,file='pst_qpg.dat',status='unknown',form='formatted')
      write(205,800)(((pst(jx,jz),txzt(jx,jz),rxzt(jx,jz),qsf(jx,jz),p(jx,jz),pp(jx,jz),g(jx,jz),gp(jx,jz)),jx=1,mxt),jz=1,mzt)
 800  format(8(1x,e12.5))
      close(205)

      open(unit=204,file='gridts.dat',status='unknown',form='formatted')
      write(204,200)( (tpt(jx,jz),jx=1,mxt),jz=1,mzt)
      close(204)
      open(unit=206,file='gridth.dat',status='unknown',form='formatted')
      write(206,200)( (tht(jx,jz),jx=1,mxt),jz=1,mzt)
      close(206)

 200  format((1x,e12.5))
      
      open(unit=301,file='bx0.dat',status='unknown',form='formatted')
      write(301,500)(((bx(jx,jz),bxdx(jx,jz),bxdz(jx,jz)),jx=1,mxt),jz=1,mzt)
 500  format(3(1x,e12.5))
      close(301)
      open(unit=302,file='bz0.dat',status='unknown',form='formatted')
      write(302,500)(((bz(jx,jz),bzdx(jx,jz),bzdz(jx,jz)),jx=1,mxt),jz=1,mzt)
      close(302)
      open(unit=303,file='c0.dat',status='unknown',form='formatted')
      write(303,500)(((cx(jx,jz),cy(jx,jz),cz(jx,jz)),jx=1,mxt),jz=1,mzt)
      close(303)
     endif

     call last_grid
     return
     end
!ws*************************************************************************************          
     subroutine last_grid
     USE DECLARE
     include 'mpif.h' 
     real*8, dimension(mbm) :: psb2
     integer, dimension(npr) :: itb1
     integer itmp,ii
     character*10 output
     character*3 cn1

      ib=0
      do jx=2,mxt-1
      do jz=2,mzt-1      
      if(pst(jx,jz).lt.psiam .and. (pst(jx-1,jz).ge.psiam .or. pst(jx+1,jz).ge.psiam &
                             .or. pst(jx,jz-1).ge.psiam .or. pst(jx,jz+1).ge.psiam)) then
        ib=ib+1
        jbx(ib)=jx
        jbz(ib)=jz
        xxb(ib)=xxt(jx)
        zzb(ib)=zzt(jz)
        psb(ib)=pst(jx,jz)
        thb(ib)=tht(jx,jz)
        tpb(ib)=tpt(jx,jz)

!        wbrx(ib)=dcos(tpb(ib))
!        wbpx(ib)=-dsin(tpb(ib))
!        wbrz(ib)=dsin(tpb(ib))
!        wbpz(ib)=dcos(tpb(ib))

        wbrx(ib)=-bz(jx,jz)/bpol(jx,jz) !cos(tp)=psdx=-bz
        wbpx(ib)=-bx(jx,jz)/bpol(jx,jz) !sin(tp)=psdz=bx
        wbrz(ib)=bx(jx,jz)/bpol(jx,jz)
        wbpz(ib)=-bz(jx,jz)/bpol(jx,jz)

      endif
      enddo
      enddo
      
      do jz=3,mzt-2    
      do jx=3,mxt-2
      if((pst(jx-1,jz) .lt.psiam .and. pst(jx+1,jz) .lt.psiam .and. pst(jx,jz-1) .lt.psiam .and. pst(jx,jz+1).lt.psiam)&
        .and. (pst(jx-2,jz).ge.psiam .or. pst(jx+2,jz).ge.psiam .or. pst(jx,jz-2).ge.psiam .or. pst(jx,jz+2).ge.psiam)) then
        ib=ib+1
        jbx(ib)=jx
        jbz(ib)=jz
        xxb(ib)=xxt(jx)
        zzb(ib)=zzt(jz)
        psb(ib)=pst(jx,jz)
        thb(ib)=tht(jx,jz)
        tpb(ib)=tpt(jx,jz)

!        wbrx(ib)=dcos(tpb(ib))
!        wbpx(ib)=-dsin(tpb(ib))
!        wbrz(ib)=dsin(tpb(ib))
!        wbpz(ib)=dcos(tpb(ib))

        wbrx(ib)=-bz(jx,jz)/bpol(jx,jz) !cos(tp)=psdx=-bz
        wbpx(ib)=-bx(jx,jz)/bpol(jx,jz) !sin(tp)=psdz=bx
        wbrz(ib)=bx(jx,jz)/bpol(jx,jz)
        wbpz(ib)=-bz(jx,jz)/bpol(jx,jz)

      endif
      enddo
      enddo
      mb=ib
      if(nrank.eq.0) then
     write(*,*) 'mb=',mb
      endif
      !if 256*256 mb ~=1400,if 512*512 mb~=2800,so mbm at last eq.4mxt+4mzt
      psbmin=minval(psb(1:mb))    
      psbmax=maxval(psb(1:mb))

!     ps(mpsa)=psia
!     ps(mpsa-1)=psbmin
!     dps=psia-psbmin
!     do js=mpsa-2,mpsa-nda,-1
!     ps(js)=ps(js+1)-dps
!     enddo

      ib2=0
      do jx=3,mxt-2
      do jz=3,mzt-2      
       if((pst(jx-2,jz) .lt.psiam .and. pst(jx+2,jz) .lt.psiam .and. pst(jx,jz-2) .lt.psiam .and. pst(jx,jz+2).lt.psiam) .and. &
         (pst(jx-3,jz) .ge.psiam  .or. pst(jx+3,jz) .ge.psiam  .or. pst(jx,jz-3) .ge.psiam  .or. pst(jx,jz+3).ge.psiam  .or.  &
          pst(jx-1,jz-3).ge.psiam .or. pst(jx+1,jz-3).ge.psiam .or. pst(jx-3,jz-1).ge.psiam .or. pst(jx-3,jz+1).ge.psiam .or.  &
          pst(jx-1,jz+3).ge.psiam .or. pst(jx+1,jz+3).ge.psiam .or. pst(jx+3,jz-1).ge.psiam .or. pst(jx+3,jz+1).ge.psiam))   then

        ib2=ib2+1
!        jb2x(ib2)=jx
!        jb2z(ib2)=jz
!        xxb2(ib2)=xxt(jx)
!        zzb2(ib2)=zzt(jz)
        psb2(ib2)=pst(jx,jz)
!        thb2(ib2)=thxz(jx,jz)
 !       nrankb(ib)=(jbx(ib)-1)/mx
      endif
      enddo
      enddo
      mb2=ib2

      psb2min=minval(psb2(1:mb2))    
      psb2max=maxval(psb2(1:mb2))
 
     ps(mpsa)=psia
     ps(mpsa-2)=psb2min
     ps(mpsa-1)=(ps(mpsa)+ps(mpsa-2))/2
     dps=(psia-psb2min)/2
     do js=mpsa-3,mpsa-nda,-1
     ps(js)=ps(js+1)-dps
     enddo

     pssm=psia-5.*dps
     pssmw=5.*dps
!     do jx=ix_first,ix_last
!     do jz=iz_first,iz_last
!     if(psi(jx,jz).lt.psia) then
!     cfsmb(jx,jz)=0.5*(1.+dtanh((psi(jx,jz)-pssm)/pssmw))/20.
!     endif
!     enddo
!     enddo
      output='cfsmb'//cn1(nrank)
      open(unit=151,file=output,status='unknown',form='formatted')
      write(151,100)((cfsmb(jx,jz),jx=ix_first,ix_last),jz=iz_first,iz_last)
100  format((1x,e12.5))
      close(151)

     do i=1,n2th+5
     xxs(i,mpsa)=xxst(i,mpsa)
     zzs(i,mpsa)=zzst(i,mpsa)
     tps(i,mpsa)=tpst(i,mpsa)
     enddo
     
     j=npsi-1
     do js=mpsa-1,mpsa-nda,-1
         do while(ps(js).lt. psival_NOVA(j))
         j=j-1
         enddo
         if(j==npsi-1) j=j-1
         ip_s(js)=j
         do i=1,n2th+5
         call interp1d3l(xxst(i,j+2),xxst(i,j+1),xxst(i,j),xxst(i,j-1), &
                        psival_NOVA(j+2),psival_NOVA(j+1),psival_NOVA(j),psival_NOVA(j-1),ps(js),xxs(i,js))
         call interp1d3l(zzst(i,j+2),zzst(i,j+1),zzst(i,j),zzst(i,j-1), &
                        psival_NOVA(j+2),psival_NOVA(j+1),psival_NOVA(j),psival_NOVA(j-1),ps(js),zzs(i,js))
         call interp1d3l(tpst(i,j+2),tpst(i,j+1),tpst(i,j),tpst(i,j-1), &
                        psival_NOVA(j+2),psival_NOVA(j+1),psival_NOVA(j),psival_NOVA(j-1),ps(js),tps(i,js))
         enddo
     enddo

      do js=mpsa,mpsa-nda,-1
      do jt=1,n2th+5
      wbxr(jt,js)=dcos(tps(jt,js))
      wbzr(jt,js)=dsin(tps(jt,js))
      wbxt(jt,js)=-dsin(tps(jt,js))
      wbzp(jt,js)=dcos(tps(jt,js))
      enddo
      enddo

     if(nrank.eq.0) then
      open(unit=161,file='xxs',status='unknown',form='formatted')
      write(161,500)((xxs(i,js),js=mpsa,mpsa-4,-1),i=1,n2th+5)
      close(161)
      open(unit=162,file='zzs',status='unknown',form='formatted')
      write(162,500)((zzs(i,js),js=mpsa,mpsa-4,-1),i=1,n2th+5)
      close(162)
      open(unit=163,file='tps',status='unknown',form='formatted')
      write(163,600)((tps(i,js),js=mpsa,mpsa-4,-1),thst(i),i=1,n2th+5)
      close(163)

500   format(5(1x,e12.5)) 
600   format(6(1x,e12.5)) 
     endif       
!     do js=mpsa,mpsa-nda,-1
!     ps(js)=psival_NOVA(js)
!     enddo

      do js=mpsa,mpsa-nda+3,-1
      dzm1=ps(js)-ps(js-1)
      dzm2=ps(js)-ps(js-2)
      dzm3=ps(js)-ps(js-3)
      f1=dzm1-dzm1**3/dzm3**2
      f2=dzm1**2-dzm1**3/dzm3
      g1=dzm2-dzm2**3/dzm3**2
      g2=dzm2**2-dzm2**3/dzm3
      ca1=f1*g2-f2*g1
      asm(js)=g2/ca1
      bsm(js)=-f2/ca1
      csm(js)=(1-asm(js)*dzm1-bsm(js)*dzm2)/dzm3
      enddo

     if(nrank.eq.0) then
     do js=mpsa,mpsa-nda,-1
     write(*,*) js,ps(js),psival_NOVA(js)
     write(*,*) js,asm(js),bsm(js),csm(js)
     enddo
     write(*,*) 'dps=',dps
     endif
!      ia1=0
!      do jx=2,mxt-1
!      do jz=2,mzt-1      
!      if(pst(jx,jz).lt.ps(mpsa) .and. pst(jx,jz).ge.ps(mpsa-1) ) then
!        ia1=ia1+1
!        jxa1(ia1)=jx
!        jza1(ia1)=jz
!        xxa1(ia1)=xxt(jx)
!        zza1(ia1)=zzt(jz)
!        psa1(ia1)=pst(jx,jz)
!        tha1(ia1)=tpt(jx,jz)
! !       nrankb(ib)=(jbx(ib)-1)/mx
!      endif
!      enddo
!      enddo
!      ma1=ia1
!
!      js=mpsa
!      do ja1=1,ma1,1
!      dzm1=psa1(ja1)-ps(js-1)
!      dzm2=psa1(ja1)-ps(js-2)
!      dzm3=psa1(ja1)-ps(js-3)
!      f1=dzm1-dzm1**3/dzm3**2
!      f2=dzm1**2-dzm1**3/dzm3
!      g1=dzm2-dzm2**3/dzm3**2
!      g2=dzm2**2-dzm2**3/dzm3
!      ca1=f1*g2-f2*g1
!      asm_a1(ja1)=g2/ca1
!      bsm_a1(ja1)=-f2/ca1
!      csm_a1(ja1)=(1-asm_a1(ja1)*dzm1-bsm_a1(ja1)*dzm2)/dzm3
!      
!      do jt=3,3+n2th
!      if(tha1(ia1).le.thst(jt) .and. tha1(ia1).gt.thst(jt-1)) it_a1(ia1)=jt 
!      enddo
!      enddo
      js=mpsa
      do jb=1,mb,1
      dzm1=psb(jb)-ps(js-1)
      dzm2=psb(jb)-ps(js-2)
      dzm3=psb(jb)-ps(js-3)
      f1=dzm1-dzm1**3/dzm3**2
      f2=dzm1**2-dzm1**3/dzm3
      g1=dzm2-dzm2**3/dzm3**2
      g2=dzm2**2-dzm2**3/dzm3
      ca1=f1*g2-f2*g1
      asm_b(jb)=g2/ca1
      bsm_b(jb)=-f2/ca1
      csm_b(jb)=(1-asm_b(jb)*dzm1-bsm_b(jb)*dzm2)/dzm3

      do jt=3,3+n2th
      if(thb(jb).le.thst(jt) .and. thb(jb).gt.thst(jt-1)) it_b(jb)=jt 
      if(it_b(jb).le.3 .and. zzb(jb).lt. 0 ) it_b(jb)=it_b(jb)+n2th
      enddo
      enddo

!      do ib=1,mb
!        i=it_b(ib)
!        do js=mpsa,mpsa-nda,-1     
!        call interp1d3l(xxs(i-2,js),xxs(i-1,js),xxs(i,js),xxs(i+1,js), &
!                        thst(i-2),thst(i-1),thst(i),thst(i+1),thb(ib),xxbs(ib,js))
!        call interp1d3l(zzs(i-2,js),zzs(i-1,js),zzs(i,js),zzs(i+1,js), &
!                        thst(i-2),thst(i-1),thst(i),thst(i+1),thb(ib),zzbs(ib,js))
!        enddo
!      enddo

        imm=0
        imp=0
        ipp=0
        ipm=0
        do jz=jzamin,jzamax
        do jx=jxamin,jxamax
        if(pst(jx,jz).lt.psiam .and. pst(jx-1,jz).ge.psiam .and. pst(jx,jz-1).ge.psiam) then
        imm=imm+1
        jxmm(imm)=jx
        jzmm(imm)=jz
        endif
        if(pst(jx,jz).lt.psiam .and. pst(jx-1,jz).ge.psiam .and. pst(jx,jz+1).ge.psiam) then
        imp=imp+1
        jxmp(imp)=jx
        jzmp(imp)=jz
        endif
        if(pst(jx,jz).lt.psiam .and. pst(jx+1,jz).ge.psiam .and. pst(jx,jz+1).ge.psiam) then
        ipp=ipp+1
        jxpp(ipp)=jx
        jzpp(ipp)=jz
        endif
        if(pst(jx,jz).lt.psiam .and. pst(jx+1,jz).ge.psiam .and. pst(jx,jz-1).ge.psiam) then
        ipm=ipm+1
        jxpm(ipm)=jx
        jzpm(ipm)=jz
        endif
        enddo
        enddo
        nmm=imm
        nmp=imp
        npp=ipp
        npm=ipm

      i=0
      do nrk=0,nsizexz-1
      ikb=0
      do jb=1,mb,1      
!      if((xxb(jb) .ge. xxt(nrkx(nrk)*mxm+ix_min(nrk)-2).and.xxb(jb).le. xxt(nrkx(nrk)*mxm+ix_max(nrk)-2)) .and. &
!         (zzb(jb) .ge. zzt(nrkz(nrk)*mzm+iz_min(nrk)-2).and.zzb(jb).le. zzt(nrkz(nrk)*mzm+iz_max(nrk)-2)) ) then
      if((xxb(jb) .ge. xxt(nrkx(nrk)*mxm+1).and.xxb(jb).lt. xxt(nrkx(nrk)*mxm+mxm)+dxt(nrkx(nrk)*mxm+mxm)) .and. &
         (zzb(jb) .ge. zzt(nrkz(nrk)*mzm+1).and.zzb(jb).lt. zzt(nrkz(nrk)*mzm+mzm)+dzt(nrkz(nrk)*mzm+mzm)) ) then
        ikb=ikb+1
        ib_nrk(ikb,nrk)=jb 
	    itb_nrk(ikb,nrk)=it_b(jb)    
        jbx_nrk(ikb,nrk)=jbx(jb)+2-nrkx(nrk)*mxm
        jbz_nrk(ikb,nrk)=jbz(jb)+2-nrkz(nrk)*mzm
        nrankb(jb)=nrk
        
      endif
      enddo   
      mb_nrk(nrk)=ikb
      if(mb_nrk(nrk).gt.0) then
      i=i+1
      nrkb(i)=nrk
      itbmin(nrk)=max(minval(itb_nrk(1:mb_nrk(nrk),nrk))-3,1)       
      itbmax(nrk)=min(maxval(itb_nrk(1:mb_nrk(nrk),nrk))+3,n2th+5)  

      endif
      enddo
      mrkb=i

      do i=1,mrkb
      itb1(i)=itbmin(nrkb(i))
      nrkb1(i)=nrkb(i)
      enddo

      do jj=1,mrkb-1
      do ii=1,mrkb-jj
      if(itb1(ii).gt.itb1(ii+1)) then
      itmp=itb1(ii)
      itb1(ii)=itb1(ii+1)
      itb1(ii+1)=itmp

      itmp=nrkb1(ii)
      nrkb1(ii)=nrkb1(ii+1)
      nrkb1(ii+1)=itmp
      endif
      enddo
      enddo

      do ii=1,mrkb
      inrkb(nrkb1(ii))=ii
      enddo

     if(nrank.eq.0) then
     do ii=1,mrkb
     write(*,*) inrkb(nrkb1(ii)),itbmin(nrkb1(ii)),itbmax(nrkb1(ii)),nrkb1(ii)
     enddo
     endif

      do js=mpsa,mpsa-nda,-1
      i=0
      do nrk=0,nsizexz-1
      iks=0
      do jt=1,n2th+5     
!      if((xxs(jt,js) .ge. xxt(nrkx(nrk)*mxm+ix_min(nrk)-2).and.xxs(jt,js).le. xxt(nrkx(nrk)*mxm+ix_max(nrk)-2)) .and. &
!         (zzs(jt,js) .ge. zzt(nrkz(nrk)*mzm+iz_min(nrk)-2).and.zzs(jt,js).le. zzt(nrkz(nrk)*mzm+iz_max(nrk)-2)) ) then
      if((xxs(jt,js) .ge. xxt(nrkx(nrk)*mxm+1).and.xxs(jt,js).lt. xxt(nrkx(nrk)*mxm+mxm)+dxt(nrkx(nrk)*mxm+mxm)) .and. &
         (zzs(jt,js) .ge. zzt(nrkz(nrk)*mzm+1).and.zzs(jt,js).lt. zzt(nrkz(nrk)*mzm+mzm)+dzt(nrkz(nrk)*mzm+mzm)) ) then
        nrankts(jt,js)=nrk
        if(jt.gt.3 .and. jt.le. 3+n2th) then
        iks=iks+1
        its_nrk(iks,js,nrk)=jt
        endif
      endif
      enddo   
      mts_nrk(js,nrk)=iks
      if(mts_nrk(js,nrk).gt.0) then
      i=i+1
      nrks(js,nrk)=nrk
      itsmin(js,nrk)=its_nrk(1,js,nrk)
      itsmax(js,nrk)=its_nrk(mts_nrk(js,nrk),js,nrk)       
      endif
      enddo
      mrks(js)=i
      enddo


     do irecv=1,mrkb    
     do js=mpsa-2,mpsa-nda,-1
     isend=1
!     nrkxsend(irecv,js)=int((xxs(itbmin(nrkb(irecv)),js)-xxt(2)+1.e-7)/(xxt(mxt)-xxt(1))*(nprx))
!     nrkzsend(irecv,js)=int((zzs(itbmin(nrkb(irecv)),js)-zzt(2)+1.e-7)/(zzt(mxt)-zzt(1))*(nprz))
!     nranksend(irecv,js,isend)=nrkzsend(irecv,js)*nprx+nrkxsend(irecv,js)
     nranksend(irecv,js,isend)=nrankts(itbmin(nrkb(irecv)),js)
     
     ittransmin(irecv,js,isend)=itbmin(nrkb(irecv))
     do while (itbmin(nrkb(irecv))+n2th .le. itsmax(js,nranksend(irecv,js,isend)) )
     ittransmax(irecv,js,isend)=itsmax(js,nranksend(irecv,js,isend))-n2th 
     ittransmin(irecv,js,isend+1)=ittransmax(irecv,js,isend)+1
     nranksend(irecv,js,isend+1)=nrankts(ittransmin(irecv,js,isend+1),js)    
     isend=isend+1
     enddo

     do while((itsmax(js,nranksend(irecv,js,isend))-itbmax(nrkb(irecv)).lt.0).and. (itsmax(js,nranksend(irecv,js,isend))-itbmax(nrkb(irecv)).gt.-nthe))
     ittransmax(irecv,js,isend)=itsmax(js,nranksend(irecv,js,isend))
!     nrkxsend(irecv,js)=int((xxs(ittransmax(irecv,js,isend)+1,js)-xxt(2))/(xxt(mxt)-xxt(1))*(nprx))
!     nrkzsend(irecv,js)=int((zzs(ittransmax(irecv,js,isend)+1,js)-zzt(2))/(zzt(mxt)-zzt(1))*(nprz))
!     nranksend(irecv,js,isend+1)=nrkzsend(irecv,js)*nprx+nrkxsend(irecv,js)
     ittransmin(irecv,js,isend+1)=ittransmax(irecv,js,isend)+1

     nranksend(irecv,js,isend+1)=nrankts(ittransmin(irecv,js,isend+1),js)
     isend=isend+1
     enddo
     ittransmax(irecv,js,isend)=itbmax(nrkb(irecv))     
     nsend(irecv,js)=isend  
     enddo
     enddo

      psia1=psiam
 !    psia1=ps(mpsa-1)
     return
     end
!ws*******************************************************************************
     subroutine read_nova
     USE DECLARE 
     real*8 rhom,rhomp
     real*8 xplmin,xplmax,zplmax,aguess,xmag,xmaj,xzmax,xatpi,xofset,aratio,bzero,curtotal,curnorm
     real*8, dimension(n2th+5,npsi):: bpst,wstx2r,wstz2r,wstx2p,wstz2p
     real*8, dimension(n2th+5,npsi,3):: fffst
     integer ipi
      include 'mpif.h' 

      do jt=1,n2th+5
      thst(jt)=pi*(jt-3)/(nthe-1)
      enddo
      ipi=nthe+2
!      open(888,file='psi_xz.dat')
!      do 20 jj=1,npsip
!      do 30 ij=1,nthe3
!      read(888,1000) j,i,psst(i,j),xxst(i,j),zzst(i,j) 
!   30 continue
!   20 continue


      open(888,file='psi_xz.dat')
      read(888,3000)
 3000 format(1h1,1x,'xplmin   xplmax   zplmax   arad   x-zero' &
      '   xmagax   xmaj  xzmax  xatpi  xofset  a-ratio b-zero' &
      '   Ip  Ipnorm' )

      read(888,4000) xplmin,xplmax,zplmax,aguess,xzero,xmag,xmaj,xzmax,xatpi,xofset,aratio,bzero,curtotal,curnorm
      
      aa=aguess    
      b0=1
      xzero=xzero/aa
      xmg=xmag/aa
      xmin=xplmin/aa
      xmax=xplmax/aa
      zmax=zplmax/aa
      zmin=-zmax/aa
      zmg=0.0
      cIp=curnorm*xzero*b0
      do j=2,npsi
      do i=3,nthe+1
      read(888,1000) js,jt,psst(i,j),xxst(i,j),zzst(i,j),bxst(i,j),bxdxst(i,j),bxdzst(i,j),bzst(i,j),bzdxst(i,j),bzdzst(i,j)

      xxst(i,j)=xxst(i,j)/aa
      zzst(i,j)=zzst(i,j)/aa
      psst(i,j)=psst(i,j)/(b0*aa**2)
      bxst(i,j)=bxst(i,j)/b0
      bzst(i,j)=bzst(i,j)/b0
      bxdxst(i,j)=bxdxst(i,j)/(b0/aa)
      bxdzst(i,j)=bxdzst(i,j)/(b0/aa)
      bzdxst(i,j)=bzdxst(i,j)/(b0/aa)
      bzdzst(i,j)=bzdzst(i,j)/(b0/aa)

      jd=(j-2)*(n2th-1)+i-2
      tst(i,j)=atan2(zzst(i,j),xxst(i,j)-xmg)          
      rst(i,j)=sqrt(zzst(i,j)**2+(xxst(i,j)-xmg)**2) 
      
      th_NOVA(jd)=thst(i)
      xx_NOVA(jd)=xxst(i,j)
      zz_NOVA(jd)=zzst(i,j)
      ps_NOVA(jd)=psst(i,j)            
      bx_NOVA(jd)=bxst(i,j)
      bz_NOVA(jd)=bzst(i,j)
      bxdx_NOVA(jd)=bxdxst(i,j)
      bxdz_NOVA(jd)=bxdzst(i,j)
      bzdx_NOVA(jd)=bzdxst(i,j)
      bzdz_NOVA(jd)=bzdzst(i,j)

      if(i.gt.3) then
      im=2*nthe+2-(i-2)
      xxst(im,j)=xxst(i,j)
      zzst(im,j)=-zzst(i,j)
      psst(im,j)=psst(i,j)
      bxst(im,j)=-bxst(i,j)
      bzst(im,j)=bzst(i,j)
      bxdxst(im,j)=-bxdxst(i,j)
      bxdzst(im,j)=bxdzst(i,j)
      bzdxst(im,j)=bzdxst(i,j)
      bzdzst(im,j)=-bzdzst(i,j)
      tst(im,j)=2*pi-tst(i,j)
      rst(im,j)=rst(i,j)

      jdm=(j-2)*(n2th-1)+im-3
      th_NOVA(jdm)=thst(im)
      xx_NOVA(jdm)=xxst(im,j)
      zz_NOVA(jdm)=zzst(im,j)
      ps_NOVA(jdm)=psst(im,j)
      bx_NOVA(jdm)=bxst(im,j)
      bz_NOVA(jdm)=bzst(im,j)
      bxdx_NOVA(jdm)=bxdxst(im,j)
      bxdz_NOVA(jdm)=bxdzst(im,j)
      bzdx_NOVA(jdm)=bzdxst(im,j)
      bzdz_NOVA(jdm)=bzdzst(im,j)
      endif       
      enddo

      xxst(1,j)=xxst(1+n2th,j)
      zzst(1,j)=zzst(1+n2th,j)
      psst(1,j)=psst(1+n2th,j)
      bxst(1,j)=bxst(1+n2th,j)
      bzst(1,j)=bzst(1+n2th,j)
      bxdxst(1,j)=bxdxst(1+n2th,j)
      bzdxst(1,j)=bzdxst(1+n2th,j)
      bxdzst(1,j)=bxdzst(1+n2th,j)
      bzdzst(1,j)=bzdzst(1+n2th,j)
      tst(1,j)=tst(1+n2th,j)-2*pi
      rst(1,j)=rst(1+n2th,j)

      xxst(2,j)=xxst(2+n2th,j)
      zzst(2,j)=zzst(2+n2th,j)
      psst(2,j)=psst(2+n2th,j)
      bxst(2,j)=bxst(2+n2th,j)
      bzst(2,j)=bzst(2+n2th,j)
      bxdxst(2,j)=bxdxst(2+n2th,j)
      bzdxst(2,j)=bzdxst(2+n2th,j)
      bxdzst(2,j)=bxdzst(2+n2th,j)
      bzdzst(2,j)=bzdzst(2+n2th,j)

      tst(2,j)=tst(2+n2th,j)-2*pi
      rst(2,j)=rst(2+n2th,j)

      xxst(3+n2th,j)=xxst(3,j)
      zzst(3+n2th,j)=zzst(3,j)
      psst(3+n2th,j)=psst(3,j)
      bxst(3+n2th,j)=bxst(3,j)
      bzst(3+n2th,j)=bzst(3,j)
      bxdxst(3+n2th,j)=bxdxst(3,j)
      bzdxst(3+n2th,j)=bzdxst(3,j)
      bxdzst(3+n2th,j)=bxdzst(3,j)
      bzdzst(3+n2th,j)=bzdzst(3,j)
      tst(3+n2th,j)=tst(3,j)+2*pi
      rst(3+n2th,j)=rst(3,j)

      xxst(4+n2th,j)=xxst(4,j)
      zzst(4+n2th,j)=zzst(4,j)
      psst(4+n2th,j)=psst(4,j)
      bxst(4+n2th,j)=bxst(4,j)
      bzst(4+n2th,j)=bzst(4,j)
      bxdxst(4+n2th,j)=bxdxst(4,j)
      bzdxst(4+n2th,j)=bzdxst(4,j)
      bxdzst(4+n2th,j)=bxdzst(4,j)
      bzdzst(4+n2th,j)=bzdzst(4,j)
      tst(4+n2th,j)=tst(4,j)+2*pi
      rst(4+n2th,j)=rst(4,j)

      xxst(5+n2th,j)=xxst(5,j)
      zzst(5+n2th,j)=zzst(5,j)
      psst(5+n2th,j)=psst(5,j)
      bxst(5+n2th,j)=bxst(5,j)
      bzst(5+n2th,j)=bzst(5,j)
      bxdxst(5+n2th,j)=bxdxst(5,j)
      bzdxst(5+n2th,j)=bzdxst(5,j)
      bxdzst(5+n2th,j)=bxdzst(5,j)
      bzdzst(5+n2th,j)=bzdzst(5,j)
      tst(5+n2th,j)=tst(5,j)+2*pi
      rst(5+n2th,j)=rst(5,j)

      zzst(ipi,j)=0
      call interp1d3l(xxst(ipi-2,j),xxst(ipi-1,j),xxst(ipi+1,j),xxst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),xxst(ipi,j))
      call interp1d3l(psst(ipi-2,j),psst(ipi-1,j),psst(ipi+1,j),psst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),psst(ipi,j))
      call interp1d3l(bxst(ipi-2,j),bxst(ipi-1,j),bxst(ipi+1,j),bxst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),bxst(ipi,j))
      call interp1d3l(bzst(ipi-2,j),bzst(ipi-1,j),bzst(ipi+1,j),bzst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),bzst(ipi,j))
      call interp1d3l(bxdxst(ipi-2,j),bxdxst(ipi-1,j),bxdxst(ipi+1,j),bxdxst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),bxdxst(ipi,j))
      call interp1d3l(bxdzst(ipi-2,j),bxdzst(ipi-1,j),bxdzst(ipi+1,j),bxdzst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),bxdzst(ipi,j))
      call interp1d3l(bzdxst(ipi-2,j),bzdxst(ipi-1,j),bzdxst(ipi+1,j),bzdxst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),bzdxst(ipi,j))
      call interp1d3l(bzdzst(ipi-2,j),bzdzst(ipi-1,j),bzdzst(ipi+1,j),bzdzst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),bzdzst(ipi,j))
      tst(ipi,j)=pi
      rst(ipi,j)=abs(xxst(ipi,j)-xmg)
      enddo
!      thst(ipi)=pi

      close(888)

      if(.not. rotation) then
      open(889,file='q_p_g.dat')
      do 41 j=1,npsi
      read(889,2100) jj,psival_NOVA(j),q_NOVA(j),qp_NOVA(j),p_NOVA(j),pp_NOVA(j),g_NOVA(j),gp_NOVA(j),f_NOVA(j),fp_NOVA(j),fb_NOVA(j),fbp_NOVA(j)
      psival_NOVA(j)=psival_NOVA(j)/(b0*aa**2)
      p_NOVA(j)=p_NOVA(j)/(b0**2)
      g_NOVA(j)=g_NOVA(j)/b0
      qp_NOVA(j)=qp_NOVA(j)*(b0*aa**2)
      pp_NOVA(j)=pp_NOVA(j)/(b0**2)*(b0*aa**2)
      gp_NOVA(j)=gp_NOVA(j)/b0*(b0*aa**2)
   41 continue
      close(889) 
      omrot_NOVA(:)=0
      omprot_NOVA(:)=0   
 
      else
      open(889,file='q_p_g.dat')
      do 40 j=1,npsi
      read(889,2000) jj,psival_NOVA(j),q_NOVA(j),qp_NOVA(j),p_NOVA(j),pp_NOVA(j),g_NOVA(j),gp_NOVA(j),f_NOVA(j),fp_NOVA(j),fb_NOVA(j),fbp_NOVA(j),omrot_NOVA(j),omprot_NOVA(j)     
      psival_NOVA(j)=psival_NOVA(j)/(b0*aa**2)
      p_NOVA(j)=p_NOVA(j)/(b0**2)
      g_NOVA(j)=g_NOVA(j)/b0
      qp_NOVA(j)=qp_NOVA(j)*(b0*aa**2)
      pp_NOVA(j)=pp_NOVA(j)/(b0**2)*(b0*aa**2)
      gp_NOVA(j)=gp_NOVA(j)/b0*(b0*aa**2)
      omrot_NOVA(j)=omrot_NOVA(j)/(b0*aa)
      omprot_NOVA(j)=omprot_NOVA(j)/(b0*aa)*(b0*aa**2)    
     
   40 continue
      close(889)      
      endif

      xxst(:,1)=xmg
      zzst(:,1)=0
      psst(:,1)=psival_NOVA(1)
      bxst(:,1)=0
      bzst(:,1)=0
      tst(:,1)=tst(:,2)
      rst(:,1)=0
 
       psia =psival_NOVA(mpsa)
       psmin=minval(ps_NOVA)
       psmax=maxval(ps_NOVA)

       qmin=minval(q_NOVA)
       qmax=maxval(q_NOVA)
       q0=q_NOVA(1)
 
      call interp1d3l(bxdxst(ipi,3),bxdxst(ipi,2),bxdxst(3,2),bxdxst(3,3), &
                      xxst(ipi,3),xxst(ipi,2),xxst(3,2),xxst(3,3),xmg,bxdxst(3,1))
      call interp1d3l(bxdzst(ipi,3),bxdzst(ipi,2),bxdzst(3,2),bxdzst(3,3), &
                      xxst(ipi,3),xxst(ipi,2),xxst(3,2),xxst(3,3),xmg,bxdzst(3,1))

      call interp1d3l(bzdxst(ipi,3),bzdxst(ipi,2),bzdxst(3,2),bzdxst(3,3), &
                      xxst(ipi,3),xxst(ipi,2),xxst(3,2),xxst(3,3),xmg,bzdxst(3,1))
      call interp1d3l(bzdzst(ipi,3),bzdzst(ipi,2),bzdzst(3,2),bzdzst(3,3), &
                      xxst(ipi,3),xxst(ipi,2),xxst(3,2),xxst(3,3),xmg,bzdzst(3,1))

!      call cubic(x1,x2,x3,x4,y1,y2,y3,y4,y,ans,dans,is,ierr)
      bxdxst(:,1)=bxdxst(3,1)
      bxdzst(:,1)=bxdzst(3,1)
      bzdxst(:,1)=bzdxst(3,1)
      bzdzst(:,1)=bzdzst(3,1)
      
      do i=1,ipi+2
      do j=2,npsi
      jd=(i-1)*(npsi-1)+j-1
      th12_NOVA(jd)=thst(i)
      xx12_NOVA(jd)=xxst(i,j)
      zz12_NOVA(jd)=zzst(i,j)

      th34_NOVA(jd)=thst(i+ipi-3)
      xx34_NOVA(jd)=xxst(i+ipi-3,j)
      zz34_NOVA(jd)=zzst(i+ipi-3,j)
      enddo
      enddo


      do j=1,npsi
      do i=1,n2th+5
      byst(i,j)=xzero*g_NOVA(j)/xxst(i,j)
      bydxst(i,j)=-xzero*(gp_NOVA(j)*bzst(i,j)+g_NOVA(j)/xxst(i,j)**2)
      bydzst(i,j)=xzero*gp_NOVA(j)*bxst(i,j)
      pdxst(i,j)=-pp_NOVA(j)*bzst(i,j)*xxst(i,j)
      pdzst(i,j)=pp_NOVA(j)*bxst(i,j)*xxst(i,j)

      cxst(i,j)=-xzero*gp_NOVA(j)*bxst(i,j)
      czst(i,j)=-xzero*gp_NOVA(j)*bzst(i,j)
      cyst(i,j)=bxdzst(i,j)-bzdxst(i,j)
!      cyst(i,j)=xxst(i,j)*pp_NOVA(j)+xzero**2*g_NOVA(j)*gp_NOVA(j)/xxst(i,j)

      uyst(i,j)=omrot_NOVA(j)*xxst(i,j)
      uydxst(i,j)=-omprot_NOVA(j)*bzst(i,j)*xxst(i,j)**2+omrot_NOVA(j)
      uydzst(i,j)= omprot_NOVA(j)*bxst(i,j)*xxst(i,j)**2

      if(j.ge.2 .and. i.ge.3 .and. i.le.n2th+2 .and. i.ne.ipi) then
      if(i.lt.ipi) jd=(j-2)*(n2th-1)+i-2
      if(i.gt.ipi) jd=(j-2)*(n2th-1)+i-3
      by_NOVA(jd)=byst(i,j)
      bydx_NOVA(jd)=bydxst(i,j)
      bydz_NOVA(jd)=bydzst(i,j)
      pdx_NOVA(jd)=pdxst(i,j)
      pdz_NOVA(jd)=pdzst(i,j)
      cy_NOVA(jd)=cyst(i,j)
      cx_NOVA(jd)=cxst(i,j)
      cz_NOVA(jd)=czst(i,j)

      uy_NOVA(jd)=uyst(i,j)
      uydx_NOVA(jd)=uydxst(i,j)
      uydz_NOVA(jd)=uydzst(i,j)

      rh_NOVA(jd)=rhom(ps_NOVA(jd))
      rhdx_NOVA(jd)=-rhomp(ps_NOVA(jd))*bzst(i,j)*xxst(i,j)
      rhdz_NOVA(jd)=rhomp(ps_NOVA(jd))*bxst(i,j)*xxst(i,j)
      
      pt_NOVA(jd)=p_NOVA(j)
      ptdx_NOVA(jd)=pdx_NOVA(jd)
      ptdz_NOVA(jd)=pdz_NOVA(jd)
      
      endif
      enddo
      enddo

      if(nrank.eq.0) then
      do j=1,npsi
      do i=3,nthe+1
      fffst(i,j,1)=cyst(i,j)*bzst(i,j)-czst(i,j)*byst(i,j)-pdxst(i,j)-uyst(i,j)*uydxst(i,j)+uyst(i,j)**2/xxst(i,j)
      fffst(i,j,2)=czst(i,j)*bxst(i,j)-cxst(i,j)*bzst(i,j)
      fffst(i,j,3)=cxst(i,j)*byst(i,j)-cyst(i,j)*bxst(i,j)-pdzst(i,j)-uyst(i,j)*uydzst(i,j)
      enddo
      enddo

      open(unit=101,file='fst_NOVA.dat',status='unknown',form='formatted')
      write(101,300)(((fffst(i,j,m),m=1,3),i=3,nthe+1),j=2,npsi)
 300  format(3(1x,e12.5))

      endif

!       psmin=minval(ps_NOVA)
!       psmax=maxval(ps_NOVA)

!       qmin=minval(q_NOVA)
!       qmax=maxval(q_NOVA)
!       q0=q_NOVA(1)
!      xmin=minval(xx_NOVA)
!      xmax=maxval(xx_NOVA)
!      zmin=minval(zz_NOVA)
!      zmax=maxval(zz_NOVA)
        epsilon=1./aratio
     
      do j=2,npsi
      do i=1,n2th+5
        tpst(i,j)=atan2(bxst(i,j),-bzst(i,j))
        if(i .gt. ipi) tpst(i,j)=tpst(i,j)+2*pi
        if(i .eq. ipi) tpst(i,j)=pi
        bpst(i,j)=sqrt(bxst(i,j)**2+bzst(i,j)**2)
        wstx2r(i,j)=-bzst(i,j)/bpst(i,j)
        wstz2r(i,j)=bxst(i,j)/bpst(i,j)
        wstx2p(i,j)=-bxst(i,j)/bpst(i,j)
        wstz2p(i,j)=-bzst(i,j)/bpst(i,j)
      enddo
      enddo




      if(nrank.eq.0) then
      open(890,file='psdat.dat')
      write(890,5000) (xx_NOVA(j),zz_NOVA(j),th_NOVA(j),ps_NOVA(j),bx_NOVA(j),bxdx_NOVA(j),bxdz_NOVA(j),bz_NOVA(j),bzdx_NOVA(j),bzdz_NOVA(j),j=1,ndat)
 5000 format(10(1x,e17.9))
      close(890)
      open(891,file='stdat.dat')
      write(891,5100) (((thst(i),psst(i,j),xxst(i,j),zzst(i,j),tpst(i,j),tst(i,j),rst(i,j)),i=1,n2th+5),j=1,npsi)
 5100 format(7(1x,e17.9))
      close(891)

      endif 
       
 1000 format(1x,i5,i5,9(1x,e17.9))
 2000 format(1x,i5,13(1x,e17.9))
 2100 format(1x,i5,11(1x,e17.9))
 4000 format(14(1x,e17.9))         
      return
      end

!ws*******************************************************************************
     subroutine read_nova_v2
     USE DECLARE 
     real*8 xplmin,xplmax,zplmax,aguess,xmag,xmaj,xzmax,xatpi,xofset,aratio,bzero,curtotal,curnorm
     real*8, dimension(nthe,npsi):: psst12,xxst12,zzst12,bxst12,bxdxst12,bxdzst12,bzst12,bzdxst12,bzdzst12
     real*8, dimension(n2th+5,npsi):: bpst,wstx2r,wstz2r,wstx2p,wstz2p
     real*8, dimension(n2th+5,npsi,3):: fffst
     integer ipi
      include 'mpif.h' 

      do jt=1,n2th+5
      thst(jt)=pi*(jt-3)/(nthe-1)
      enddo
      ipi=nthe+2
!      open(888,file='psi_xz.dat')
!      do 20 jj=1,npsip
!      do 30 ij=1,nthe3
!      read(888,1000) j,i,psst(i,j),xxst(i,j),zzst(i,j) 
!   30 continue
!   20 continue

      open(unit=889,file='eq_nova',status='unknown',form='unformatted')
      read(889) xplmin,xplmax,zplmax,aguess,xzero,xmg,xmaj,xzmax,xatpi,xofset,aratio
      read(889) psst12,xxst12,zzst12,bxst12,bxdxst12,bxdzst12,bzst12,bzdxst12,bzdzst12
      read(889) psival_NOVA,q_NOVA,qp_NOVA,p_NOVA,pp_NOVA,g_NOVA,gp_NOVA,f_NOVA,fp_NOVA,fb_NOVA,fbp_NOVA,omrot_NOVA,omprot_NOVA     
      close(889)      
    
      
      do j=2,npsi
      do i=3,nthe+1
      psst(i,j)=psst12(i,j)
      xxst(i,j)=xxst12(i,j)
      zzst(i,j)=zzst12(i,j)
      bxst(i,j)=bxst12(i,j)
      bxdxst(i,j)=bxdxst12(i,j)
      bxdzst(i,j)=bxdzst12(i,j)
      bzst(i,j)=bzst12(i,j)
      bzdxst(i,j)=bzdxst12(i,j)
      bzdzst(i,j)=bzdzst12(i,j)

      tst(i,j)=atan2(zzst(i,j),xxst(i,j)-xmg)          
      rst(i,j)=sqrt(zzst(i,j)**2+(xxst(i,j)-xmg)**2) 
      
      jd=(j-2)*(n2th-1)+i-2   
      th_NOVA(jd)=thst(i)
      xx_NOVA(jd)=xxst(i,j)
      zz_NOVA(jd)=zzst(i,j)
      ps_NOVA(jd)=psst(i,j)            
      bx_NOVA(jd)=bxst(i,j)
      bz_NOVA(jd)=bzst(i,j)
      bxdx_NOVA(jd)=bxdxst(i,j)
      bxdz_NOVA(jd)=bxdzst(i,j)
      bzdx_NOVA(jd)=bzdxst(i,j)
      bzdz_NOVA(jd)=bzdzst(i,j)

      if(i.gt.3) then
      im=2*nthe+2-(i-2)
      xxst(im,j)=xxst(i,j)
      zzst(im,j)=-zzst(i,j)
      psst(im,j)=psst(i,j)
      bxst(im,j)=-bxst(i,j)
      bzst(im,j)=bzst(i,j)
      bxdxst(im,j)=-bxdxst(i,j)
      bxdzst(im,j)=bxdzst(i,j)
      bzdxst(im,j)=bzdxst(i,j)
      bzdzst(im,j)=-bzdzst(i,j)
      tst(im,j)=2*pi-tst(i,j)
      rst(im,j)=rst(i,j)

      jdm=(j-2)*(n2th-1)+im-3
      th_NOVA(jdm)=thst(im)
      xx_NOVA(jdm)=xxst(im,j)
      zz_NOVA(jdm)=zzst(im,j)
      ps_NOVA(jdm)=psst(im,j)
      bx_NOVA(jdm)=bxst(im,j)
      bz_NOVA(jdm)=bzst(im,j)
      bxdx_NOVA(jdm)=bxdxst(im,j)
      bxdz_NOVA(jdm)=bxdzst(im,j)
      bzdx_NOVA(jdm)=bzdxst(im,j)
      bzdz_NOVA(jdm)=bzdzst(im,j)
      endif       
      enddo

      xxst(1,j)=xxst(1+n2th,j)
      zzst(1,j)=zzst(1+n2th,j)
      psst(1,j)=psst(1+n2th,j)
      bxst(1,j)=bxst(1+n2th,j)
      bzst(1,j)=bzst(1+n2th,j)
      tst(1,j)=tst(1+n2th,j)-2*pi
      rst(1,j)=rst(1+n2th,j)

      xxst(2,j)=xxst(2+n2th,j)
      zzst(2,j)=zzst(2+n2th,j)
      psst(2,j)=psst(2+n2th,j)
      bxst(2,j)=bxst(2+n2th,j)
      bzst(2,j)=bzst(2+n2th,j)
      tst(2,j)=tst(2+n2th,j)-2*pi
      rst(2,j)=rst(2+n2th,j)

      xxst(3+n2th,j)=xxst(3,j)
      zzst(3+n2th,j)=zzst(3,j)
      psst(3+n2th,j)=psst(3,j)
      bxst(3+n2th,j)=bxst(3,j)
      bzst(3+n2th,j)=bzst(3,j)
      tst(3+n2th,j)=tst(3,j)+2*pi
      rst(3+n2th,j)=rst(3,j)

      xxst(4+n2th,j)=xxst(4,j)
      zzst(4+n2th,j)=zzst(4,j)
      psst(4+n2th,j)=psst(4,j)
      bxst(4+n2th,j)=bxst(4,j)
      bzst(4+n2th,j)=bzst(4,j)
      tst(4+n2th,j)=tst(4,j)+2*pi
      rst(4+n2th,j)=rst(4,j)

      xxst(5+n2th,j)=xxst(5,j)
      zzst(5+n2th,j)=zzst(5,j)
      psst(5+n2th,j)=psst(5,j)
      bxst(5+n2th,j)=bxst(5,j)
      bzst(5+n2th,j)=bzst(5,j)
      tst(5+n2th,j)=tst(5,j)+2*pi
      rst(5+n2th,j)=rst(5,j)

      zzst(ipi,j)=0
      call interp1d3l(xxst(ipi-2,j),xxst(ipi-1,j),xxst(ipi+1,j),xxst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),xxst(ipi,j))
      call interp1d3l(psst(ipi-2,j),psst(ipi-1,j),psst(ipi+1,j),psst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),psst(ipi,j))
      call interp1d3l(bxst(ipi-2,j),bxst(ipi-1,j),bxst(ipi+1,j),bxst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),bxst(ipi,j))
      call interp1d3l(bzst(ipi-2,j),bzst(ipi-1,j),bzst(ipi+1,j),bzst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),bzst(ipi,j))
      call interp1d3l(bxdxst(ipi-2,j),bxdxst(ipi-1,j),bxdxst(ipi+1,j),bxdxst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),bxdxst(ipi,j))
      call interp1d3l(bxdzst(ipi-2,j),bxdzst(ipi-1,j),bxdzst(ipi+1,j),bxdzst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),bxdzst(ipi,j))
      call interp1d3l(bzdxst(ipi-2,j),bzdxst(ipi-1,j),bzdxst(ipi+1,j),bzdxst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),bzdxst(ipi,j))
      call interp1d3l(bzdzst(ipi-2,j),bzdzst(ipi-1,j),bzdzst(ipi+1,j),bzdzst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),bzdzst(ipi,j))
      tst(ipi,j)=pi
      rst(ipi,j)=abs(xxst(ipi,j)-xmg)
      enddo
!      thst(ipi)=pi

      close(888)

      xxst(:,1)=xmg
      zzst(:,1)=0
      psst(:,1)=psival_NOVA(1)
      bxst(:,1)=0
      bzst(:,1)=0
      tst(:,1)=tst(:,2)
      rst(:,1)=0

      call interp1d3l(bxdxst(ipi,3),bxdxst(ipi,2),bxdxst(3,2),bxdxst(3,3), &
                      xxst(ipi,3),xxst(ipi,2),xxst(3,2),xxst(3,3),xmg,bxdxst(3,1))
      call interp1d3l(bxdzst(ipi,3),bxdzst(ipi,2),bxdzst(3,2),bxdzst(3,3), &
                      xxst(ipi,3),xxst(ipi,2),xxst(3,2),xxst(3,3),xmg,bxdzst(3,1))

      call interp1d3l(bzdxst(ipi,3),bzdxst(ipi,2),bzdxst(3,2),bzdxst(3,3), &
                      xxst(ipi,3),xxst(ipi,2),xxst(3,2),xxst(3,3),xmg,bzdxst(3,1))
      call interp1d3l(bzdzst(ipi,3),bzdzst(ipi,2),bzdzst(3,2),bzdzst(3,3), &
                      xxst(ipi,3),xxst(ipi,2),xxst(3,2),xxst(3,3),xmg,bzdzst(3,1))

!      call cubic(x1,x2,x3,x4,y1,y2,y3,y4,y,ans,dans,is,ierr)
      bxdxst(:,1)=bxdxst(3,1)
      bxdzst(:,1)=bxdzst(3,1)
      bzdxst(:,1)=bzdxst(3,1)
      bzdzst(:,1)=bzdzst(3,1)
      
      do i=1,ipi+2
      do j=2,npsi
      jd=(i-1)*(npsi-1)+j-1
      th12_NOVA(jd)=thst(i)
      xx12_NOVA(jd)=xxst(i,j)
      zz12_NOVA(jd)=zzst(i,j)

      th34_NOVA(jd)=thst(i+ipi-3)
      xx34_NOVA(jd)=xxst(i+ipi-3,j)
      zz34_NOVA(jd)=zzst(i+ipi-3,j)
      enddo
      enddo


      do j=1,npsi
      do i=1,n2th+5
      byst(i,j)=xzero*g_NOVA(j)/xxst(i,j)
      bydxst(i,j)=-xzero*(gp_NOVA(j)*bzst(i,j)+g_NOVA(j)/xxst(i,j)**2)
      bydzst(i,j)=xzero*gp_NOVA(j)*bxst(i,j)
      pdxst(i,j)=-pp_NOVA(j)*bzst(i,j)*xxst(i,j)
      pdzst(i,j)=pp_NOVA(j)*bxst(i,j)*xxst(i,j)

      cxst(i,j)=-xzero*gp_NOVA(j)*bxst(i,j)
      czst(i,j)=-xzero*gp_NOVA(j)*bzst(i,j)
      cyst(i,j)=bxdzst(i,j)-bzdxst(i,j)
!      cyst(i,j)=xxst(i,j)*pp_NOVA(j)+xzero**2*g_NOVA(j)*gp_NOVA(j)/xxst(i,j)

      uyst(i,j)=omrot_NOVA(j)*xxst(i,j)
      uydxst(i,j)=-omprot_NOVA(j)*bzst(i,j)*xxst(i,j)**2+omrot_NOVA(j)
      uydzst(i,j)= omprot_NOVA(j)*bxst(i,j)*xxst(i,j)**2

      if(j.ge.2 .and. i.ge.3 .and. i.le.n2th+2 .and. i.ne.ipi) then
      if(i.lt.ipi) jd=(j-2)*(n2th-1)+i-2
      if(i.gt.ipi) jd=(j-2)*(n2th-1)+i-3
      by_NOVA(jd)=byst(i,j)
      bydx_NOVA(jd)=bydxst(i,j)
      bydz_NOVA(jd)=bydzst(i,j)
      pdx_NOVA(jd)=pdxst(i,j)
      pdz_NOVA(jd)=pdzst(i,j)
      cy_NOVA(jd)=cyst(i,j)
      cx_NOVA(jd)=cxst(i,j)
      cz_NOVA(jd)=czst(i,j)

      uy_NOVA(jd)=uyst(i,j)
      uydx_NOVA(jd)=uydxst(i,j)
      uydz_NOVA(jd)=uydzst(i,j)
      endif
      enddo
      enddo

      if(nrank.eq.0) then
      do j=1,npsi
      do i=3,nthe+1
      fffst(i,j,1)=cyst(i,j)*bzst(i,j)-czst(i,j)*byst(i,j)-pdxst(i,j)-uyst(i,j)*uydxst(i,j)+uyst(i,j)**2/xxst(i,j)
      fffst(i,j,2)=czst(i,j)*bxst(i,j)-cxst(i,j)*bzst(i,j)
      fffst(i,j,3)=cxst(i,j)*byst(i,j)-cyst(i,j)*bxst(i,j)-pdzst(i,j)-uyst(i,j)*uydzst(i,j)
      enddo
      enddo

      open(unit=101,file='fst_NOVA.dat',status='unknown',form='formatted')
      write(101,300)(((fffst(i,j,m),m=1,3),i=3,nthe+1),j=2,npsi)
 300  format(3(1x,e12.5))

      endif

       xmin=xplmin
       xmax=xplmax
       zmax=zplmax
       zmin=-zmax
       zmg=0.0
       psmin=minval(ps_NOVA)
       psmax=maxval(ps_NOVA)

       qmin=minval(q_NOVA)
       qmax=maxval(q_NOVA)
       q0=q_NOVA(1)
!      xmin=minval(xx_NOVA)
!      xmax=maxval(xx_NOVA)
!      zmin=minval(zz_NOVA)
!      zmax=maxval(zz_NOVA)
        epsilon=1./aratio
     
      do j=2,npsi
      do i=1,n2th+5
        tpst(i,j)=atan2(bxst(i,j),-bzst(i,j))
        if(i .gt. ipi) tpst(i,j)=tpst(i,j)+2*pi
        if(i .eq. ipi) tpst(i,j)=pi
        bpst(i,j)=sqrt(bxst(i,j)**2+bzst(i,j)**2)
        wstx2r(i,j)=-bzst(i,j)/bpst(i,j)
        wstz2r(i,j)=bxst(i,j)/bpst(i,j)
        wstx2p(i,j)=-bxst(i,j)/bpst(i,j)
        wstz2p(i,j)=-bzst(i,j)/bpst(i,j)
      enddo
      enddo




      if(nrank.eq.0) then
      open(890,file='psdat.dat')
      write(890,5000) (xx_NOVA(j),zz_NOVA(j),th_NOVA(j),ps_NOVA(j),bx_NOVA(j),bxdx_NOVA(j),bxdz_NOVA(j),bz_NOVA(j),bzdx_NOVA(j),bzdz_NOVA(j),j=1,ndat)
 5000 format(10(1x,e17.9))
      close(890)
      open(891,file='stdat.dat')
      write(891,5100) (((thst(i),psst(i,j),xxst(i,j),zzst(i,j),tpst(i,j),tst(i,j),rst(i,j)),i=1,n2th+5),j=1,npsi)
 5100 format(7(1x,e17.9))
      close(891)

      endif 
              
      return
      end
!ws*******************************************************************************
     subroutine read_nova_tm
     USE DECLARE 
     real*8 xplmin,xplmax,zplmax,aguess,xmag,xmaj,xzmax,xatpi,xofset,aratio,bzero,curtotal,curnorm
     real*8, dimension(n2th+5,npsi):: bpst,wstx2r,wstz2r,wstx2p,wstz2p
     real*8, dimension(n2th+5,npsi):: ptst,ptdxst,ptdzst,rhst,rhdxst,rhdzst
     real*8, dimension(n2th+5,npsi,3):: fffst
     integer ipi
      include 'mpif.h' 

      do jt=1,n2th+5
      thst(jt)=pi*(jt-3)/(nthe-1)
      enddo
      ipi=nthe+2
!      open(888,file='psi_xz.dat')
!      do 20 jj=1,npsip
!      do 30 ij=1,nthe3
!      read(888,1000) j,i,psst(i,j),xxst(i,j),zzst(i,j) 
!   30 continue
!   20 continue


      open(888,file='psi_xz_tm.dat')
      read(888,3000)
 3000 format(1h1,1x,'xplmin   xplmax   zplmax   arad   x-zero' &
      '   xmagax   xmaj  xzmax  xatpi  xofset  a-ratio b-zero' &
      '   Ip  Ipnorm' )

      read(888,4000) xplmin,xplmax,zplmax,aguess,xzero,xmag,xmaj,xzmax,xatpi,xofset,aratio,bzero,curtotal,curnorm
      
      aa=aguess    
      b0=1
      xzero=xzero/aa
      xmg=xmag/aa
      xmin=xplmin/aa
      xmax=xplmax/aa
      zmax=zplmax/aa
      zmin=-zmax/aa
      zmg=0.0
      cIp=curnorm*xzero*b0
      do j=2,npsi
      do i=3,nthe+1
      read(888,1500) js,jt,psst(i,j),xxst(i,j),zzst(i,j),bxst(i,j),bxdxst(i,j),bxdzst(i,j),bzst(i,j),bzdxst(i,j),bzdzst(i,j) &
                      ,ptst(i,j),ptdxst(i,j),ptdzst(i,j),rhst(i,j),rhdxst(i,j),rhdzst(i,j)

      xxst(i,j)=xxst(i,j)/aa
      zzst(i,j)=zzst(i,j)/aa
      psst(i,j)=psst(i,j)/(b0*aa**2)
      bxst(i,j)=bxst(i,j)/b0
      bzst(i,j)=bzst(i,j)/b0
      bxdxst(i,j)=bxdxst(i,j)/(b0/aa)
      bxdzst(i,j)=bxdzst(i,j)/(b0/aa)
      bzdxst(i,j)=bzdxst(i,j)/(b0/aa)
      bzdzst(i,j)=bzdzst(i,j)/(b0/aa)
      
      ptst(i,j)=ptst(i,j)/b0**2
      ptdxst(i,j)=ptdxst(i,j)/(b0**2/aa)
      ptdzst(i,j)=ptdzst(i,j)/(b0**2/aa)
      rhdxst(i,j)=rhdxst(i,j)*aa
      rhdzst(i,j)=rhdzst(i,j)*aa

      jd=(j-2)*(n2th-1)+i-2
      tst(i,j)=atan2(zzst(i,j),xxst(i,j)-xmg)          
      rst(i,j)=sqrt(zzst(i,j)**2+(xxst(i,j)-xmg)**2) 
      
      th_NOVA(jd)=thst(i)
      xx_NOVA(jd)=xxst(i,j)
      zz_NOVA(jd)=zzst(i,j)
      ps_NOVA(jd)=psst(i,j)            
      bx_NOVA(jd)=bxst(i,j)
      bz_NOVA(jd)=bzst(i,j)
      bxdx_NOVA(jd)=bxdxst(i,j)
      bxdz_NOVA(jd)=bxdzst(i,j)
      bzdx_NOVA(jd)=bzdxst(i,j)
      bzdz_NOVA(jd)=bzdzst(i,j)

      pt_NOVA(jd)=ptst(i,j)
      ptdx_NOVA(jd)=ptdxst(i,j)
      ptdz_NOVA(jd)=ptdzst(i,j)

      rh_NOVA(jd)=rhst(i,j)
      rhdx_NOVA(jd)=rhdxst(i,j)
      rhdz_NOVA(jd)=rhdzst(i,j)

      if(i.gt.3) then
      im=2*nthe+2-(i-2)
      xxst(im,j)=xxst(i,j)
      zzst(im,j)=-zzst(i,j)
      psst(im,j)=psst(i,j)
      bxst(im,j)=-bxst(i,j)
      bzst(im,j)=bzst(i,j)
      bxdxst(im,j)=-bxdxst(i,j)
      bxdzst(im,j)=bxdzst(i,j)
      bzdxst(im,j)=bzdxst(i,j)
      bzdzst(im,j)=-bzdzst(i,j)

      ptst(im,j)=ptst(i,j)
      ptdxst(im,j)=ptdxst(i,j)
      ptdzst(im,j)=-ptdzst(i,j)
      rhst(im,j)=rhst(i,j)
      rhdxst(im,j)=rhdxst(i,j)
      rhdzst(im,j)=-rhdzst(i,j)

      tst(im,j)=2*pi-tst(i,j)
      rst(im,j)=rst(i,j)

      jdm=(j-2)*(n2th-1)+im-3
      th_NOVA(jdm)=thst(im)
      xx_NOVA(jdm)=xxst(im,j)
      zz_NOVA(jdm)=zzst(im,j)
      ps_NOVA(jdm)=psst(im,j)
      bx_NOVA(jdm)=bxst(im,j)
      bz_NOVA(jdm)=bzst(im,j)
      bxdx_NOVA(jdm)=bxdxst(im,j)
      bxdz_NOVA(jdm)=bxdzst(im,j)
      bzdx_NOVA(jdm)=bzdxst(im,j)
      bzdz_NOVA(jdm)=bzdzst(im,j)

      pt_NOVA(jdm)=ptst(im,j)
      ptdx_NOVA(jdm)=ptdxst(im,j)
      ptdz_NOVA(jdm)=ptdzst(im,j)
      rh_NOVA(jdm)=rhst(im,j)
      rhdx_NOVA(jdm)=rhdxst(im,j)
      rhdz_NOVA(jdm)=rhdzst(im,j)
      endif       
      enddo

      xxst(1,j)=xxst(1+n2th,j)
      zzst(1,j)=zzst(1+n2th,j)
      psst(1,j)=psst(1+n2th,j)
      bxst(1,j)=bxst(1+n2th,j)
      bzst(1,j)=bzst(1+n2th,j)
      bxdxst(1,j)=bxdxst(1+n2th,j)
      bzdxst(1,j)=bzdxst(1+n2th,j)
      bxdzst(1,j)=bxdzst(1+n2th,j)
      bzdzst(1,j)=bzdzst(1+n2th,j)
      tst(1,j)=tst(1+n2th,j)-2*pi
      rst(1,j)=rst(1+n2th,j)

      ptst(1,j)=ptst(1+n2th,j)
      ptdxst(1,j)=ptdxst(1+n2th,j)
      ptdzst(1,j)=ptdzst(1+n2th,j)
      rhst(1,j)=rhst(1+n2th,j)
      rhdxst(1,j)=rhdxst(1+n2th,j)
      rhdzst(1,j)=rhdzst(1+n2th,j)


      xxst(2,j)=xxst(2+n2th,j)
      zzst(2,j)=zzst(2+n2th,j)
      psst(2,j)=psst(2+n2th,j)
      bxst(2,j)=bxst(2+n2th,j)
      bzst(2,j)=bzst(2+n2th,j)
      bxdxst(2,j)=bxdxst(2+n2th,j)
      bzdxst(2,j)=bzdxst(2+n2th,j)
      bxdzst(2,j)=bxdzst(2+n2th,j)
      bzdzst(2,j)=bzdzst(2+n2th,j)
      tst(2,j)=tst(2+n2th,j)-2*pi
      rst(2,j)=rst(2+n2th,j)

      ptst(2,j)=ptst(2+n2th,j)
      ptdxst(2,j)=ptdxst(2+n2th,j)
      ptdzst(2,j)=ptdzst(2+n2th,j)
      rhst(2,j)=rhst(2+n2th,j)
      rhdxst(2,j)=rhdxst(2+n2th,j)
      rhdzst(2,j)=rhdzst(2+n2th,j)

      xxst(3+n2th,j)=xxst(3,j)
      zzst(3+n2th,j)=zzst(3,j)
      psst(3+n2th,j)=psst(3,j)
      bxst(3+n2th,j)=bxst(3,j)
      bzst(3+n2th,j)=bzst(3,j)
      bxdxst(3+n2th,j)=bxdxst(3,j)
      bzdxst(3+n2th,j)=bzdxst(3,j)
      bxdzst(3+n2th,j)=bxdzst(3,j)
      bzdzst(3+n2th,j)=bzdzst(3,j)
      tst(3+n2th,j)=tst(3,j)+2*pi
      rst(3+n2th,j)=rst(3,j)

      ptst(3+n2th,j)=ptst(3,j)
      ptdxst(3+n2th,j)=ptdxst(3,j)
      ptdzst(3+n2th,j)=ptdzst(3,j)
      rhst(3+n2th,j)=rhst(3,j)
      rhdxst(3+n2th,j)=rhdxst(3,j)
      rhdzst(3+n2th,j)=rhdzst(3,j)

      xxst(4+n2th,j)=xxst(4,j)
      zzst(4+n2th,j)=zzst(4,j)
      psst(4+n2th,j)=psst(4,j)
      bxst(4+n2th,j)=bxst(4,j)
      bzst(4+n2th,j)=bzst(4,j)
      bxdxst(4+n2th,j)=bxdxst(4,j)
      bzdxst(4+n2th,j)=bzdxst(4,j)
      bxdzst(4+n2th,j)=bxdzst(4,j)
      bzdzst(4+n2th,j)=bzdzst(4,j)
      tst(4+n2th,j)=tst(4,j)+2*pi
      rst(4+n2th,j)=rst(4,j)

      ptst(4+n2th,j)=ptst(4,j)
      ptdxst(4+n2th,j)=ptdxst(4,j)
      ptdzst(4+n2th,j)=ptdzst(4,j)
      rhst(4+n2th,j)=rhst(4,j)
      rhdxst(4+n2th,j)=rhdxst(4,j)
      rhdzst(4+n2th,j)=rhdzst(4,j)

      xxst(5+n2th,j)=xxst(5,j)
      zzst(5+n2th,j)=zzst(5,j)
      psst(5+n2th,j)=psst(5,j)
      bxst(5+n2th,j)=bxst(5,j)
      bzst(5+n2th,j)=bzst(5,j)
      bxdxst(5+n2th,j)=bxdxst(5,j)
      bzdxst(5+n2th,j)=bzdxst(5,j)
      bxdzst(5+n2th,j)=bxdzst(5,j)
      bzdzst(5+n2th,j)=bzdzst(5,j)
      tst(5+n2th,j)=tst(5,j)+2*pi
      rst(5+n2th,j)=rst(5,j)

      ptst(5+n2th,j)=ptst(5,j)
      ptdxst(5+n2th,j)=ptdxst(5,j)
      ptdzst(5+n2th,j)=ptdzst(5,j)
      rhst(5+n2th,j)=rhst(5,j)
      rhdxst(5+n2th,j)=rhdxst(5,j)
      rhdzst(5+n2th,j)=rhdzst(5,j)

      zzst(ipi,j)=0
      call interp1d3l(xxst(ipi-2,j),xxst(ipi-1,j),xxst(ipi+1,j),xxst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),xxst(ipi,j))
      call interp1d3l(psst(ipi-2,j),psst(ipi-1,j),psst(ipi+1,j),psst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),psst(ipi,j))
      call interp1d3l(bxst(ipi-2,j),bxst(ipi-1,j),bxst(ipi+1,j),bxst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),bxst(ipi,j))
      call interp1d3l(bzst(ipi-2,j),bzst(ipi-1,j),bzst(ipi+1,j),bzst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),bzst(ipi,j))
      call interp1d3l(bxdxst(ipi-2,j),bxdxst(ipi-1,j),bxdxst(ipi+1,j),bxdxst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),bxdxst(ipi,j))
      call interp1d3l(bxdzst(ipi-2,j),bxdzst(ipi-1,j),bxdzst(ipi+1,j),bxdzst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),bxdzst(ipi,j))
      call interp1d3l(bzdxst(ipi-2,j),bzdxst(ipi-1,j),bzdxst(ipi+1,j),bzdxst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),bzdxst(ipi,j))
      call interp1d3l(bzdzst(ipi-2,j),bzdzst(ipi-1,j),bzdzst(ipi+1,j),bzdzst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),bzdzst(ipi,j))

      call interp1d3l(ptst(ipi-2,j),ptst(ipi-1,j),ptst(ipi+1,j),ptst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),ptst(ipi,j))
      call interp1d3l(ptdxst(ipi-2,j),ptdxst(ipi-1,j),ptdxst(ipi+1,j),ptdxst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),ptdxst(ipi,j))
      call interp1d3l(ptdzst(ipi-2,j),ptdzst(ipi-1,j),ptdzst(ipi+1,j),ptdzst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),ptdzst(ipi,j))
      call interp1d3l(rhst(ipi-2,j),rhst(ipi-1,j),rhst(ipi+1,j),rhst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),rhst(ipi,j))
      call interp1d3l(rhdxst(ipi-2,j),rhdxst(ipi-1,j),rhdxst(ipi+1,j),rhdxst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),rhdxst(ipi,j))
      call interp1d3l(rhdzst(ipi-2,j),rhdzst(ipi-1,j),rhdzst(ipi+1,j),rhdzst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),rhdzst(ipi,j))

      tst(ipi,j)=pi
      rst(ipi,j)=abs(xxst(ipi,j)-xmg)
      enddo
!      thst(ipi)=pi

      close(888)

      if(.not. rotation) then
      open(889,file='q_p_g.dat')
      do 41 j=1,npsi
      read(889,2100) jj,psival_NOVA(j),q_NOVA(j),qp_NOVA(j),p_NOVA(j),pp_NOVA(j),g_NOVA(j),gp_NOVA(j),f_NOVA(j),fp_NOVA(j),fb_NOVA(j),fbp_NOVA(j)
      psival_NOVA(j)=psival_NOVA(j)/(b0*aa**2)
      p_NOVA(j)=p_NOVA(j)/(b0**2)
      g_NOVA(j)=g_NOVA(j)/b0
      qp_NOVA(j)=qp_NOVA(j)*(b0*aa**2)
      pp_NOVA(j)=pp_NOVA(j)/(b0**2)*(b0*aa**2)
      gp_NOVA(j)=gp_NOVA(j)/b0*(b0*aa**2)
   41 continue
      close(889) 
      omrot_NOVA(:)=0
      omprot_NOVA(:)=0   
 
      else
      open(889,file='q_p_g.dat')
      do 40 j=1,npsi
      read(889,2000) jj,psival_NOVA(j),q_NOVA(j),qp_NOVA(j),p_NOVA(j),pp_NOVA(j),g_NOVA(j),gp_NOVA(j),f_NOVA(j),fp_NOVA(j),fb_NOVA(j),fbp_NOVA(j),omrot_NOVA(j),omprot_NOVA(j)     
      psival_NOVA(j)=psival_NOVA(j)/(b0*aa**2)
      p_NOVA(j)=p_NOVA(j)/(b0**2)
      g_NOVA(j)=g_NOVA(j)/b0
      qp_NOVA(j)=qp_NOVA(j)*(b0*aa**2)
      pp_NOVA(j)=pp_NOVA(j)/(b0**2)*(b0*aa**2)
      gp_NOVA(j)=gp_NOVA(j)/b0*(b0*aa**2)
      omrot_NOVA(j)=omrot_NOVA(j)/(b0*aa)
      omprot_NOVA(j)=omprot_NOVA(j)/(b0*aa)*(b0*aa**2)    
     
   40 continue
      close(889)      
      endif

      xxst(:,1)=xmg
      zzst(:,1)=0
      psst(:,1)=psival_NOVA(1)
      bxst(:,1)=0
      bzst(:,1)=0
      tst(:,1)=tst(:,2)
      rst(:,1)=0

      call interp1d3l(bxdxst(ipi,3),bxdxst(ipi,2),bxdxst(3,2),bxdxst(3,3), &
                      xxst(ipi,3),xxst(ipi,2),xxst(3,2),xxst(3,3),xmg,bxdxst(3,1))
      call interp1d3l(bxdzst(ipi,3),bxdzst(ipi,2),bxdzst(3,2),bxdzst(3,3), &
                      xxst(ipi,3),xxst(ipi,2),xxst(3,2),xxst(3,3),xmg,bxdzst(3,1))

      call interp1d3l(bzdxst(ipi,3),bzdxst(ipi,2),bzdxst(3,2),bzdxst(3,3), &
                      xxst(ipi,3),xxst(ipi,2),xxst(3,2),xxst(3,3),xmg,bzdxst(3,1))
      call interp1d3l(bzdzst(ipi,3),bzdzst(ipi,2),bzdzst(3,2),bzdzst(3,3), &
                      xxst(ipi,3),xxst(ipi,2),xxst(3,2),xxst(3,3),xmg,bzdzst(3,1))

!      call cubic(x1,x2,x3,x4,y1,y2,y3,y4,y,ans,dans,is,ierr)
      bxdxst(:,1)=bxdxst(3,1)
      bxdzst(:,1)=bxdzst(3,1)
      bzdxst(:,1)=bzdxst(3,1)
      bzdzst(:,1)=bzdzst(3,1)

      
      call interp1d3l(ptst(ipi,3),ptst(ipi,2),ptst(3,2),ptst(3,3), &
                      xxst(ipi,3),xxst(ipi,2),xxst(3,2),xxst(3,3),xmg,ptst(3,1))
      call interp1d3l(ptdxst(ipi,3),ptdxst(ipi,2),ptdxst(3,2),ptdxst(3,3), &
                      xxst(ipi,3),xxst(ipi,2),xxst(3,2),xxst(3,3),xmg,ptdxst(3,1))
      call interp1d3l(ptdzst(ipi,3),ptdzst(ipi,2),ptdzst(3,2),ptdzst(3,3), &
                      xxst(ipi,3),xxst(ipi,2),xxst(3,2),xxst(3,3),xmg,ptdzst(3,1))
      call interp1d3l(rhst(ipi,3),rhst(ipi,2),rhst(3,2),rhst(3,3), &
                      xxst(ipi,3),xxst(ipi,2),xxst(3,2),xxst(3,3),xmg,rhst(3,1))
      call interp1d3l(rhdxst(ipi,3),rhdxst(ipi,2),rhdxst(3,2),rhdxst(3,3), &
                      xxst(ipi,3),xxst(ipi,2),xxst(3,2),xxst(3,3),xmg,rhdxst(3,1))
      call interp1d3l(rhdzst(ipi,3),rhdzst(ipi,2),rhdzst(3,2),rhdzst(3,3), &
                      xxst(ipi,3),xxst(ipi,2),xxst(3,2),xxst(3,3),xmg,rhdzst(3,1))
      
      ptst(:,1)=ptst(3,1)
      ptdxst(:,1)=ptdxst(3,1)
      ptdzst(:,1)=ptdzst(3,1)
      rhst(:,1)=rhst(3,1)
      rhdxst(:,1)=rhdxst(3,1)
      rhdzst(:,1)=rhdzst(3,1)
                           
      do i=1,ipi+2
      do j=2,npsi
      jd=(i-1)*(npsi-1)+j-1
      th12_NOVA(jd)=thst(i)
      xx12_NOVA(jd)=xxst(i,j)
      zz12_NOVA(jd)=zzst(i,j)

      th34_NOVA(jd)=thst(i+ipi-3)
      xx34_NOVA(jd)=xxst(i+ipi-3,j)
      zz34_NOVA(jd)=zzst(i+ipi-3,j)
      enddo
      enddo


      do j=1,npsi
      do i=1,n2th+5
      byst(i,j)=xzero*g_NOVA(j)/xxst(i,j)
      bydxst(i,j)=-xzero*(gp_NOVA(j)*bzst(i,j)+g_NOVA(j)/xxst(i,j)**2)
      bydzst(i,j)=xzero*gp_NOVA(j)*bxst(i,j)

      cxst(i,j)=-xzero*gp_NOVA(j)*bxst(i,j)
      czst(i,j)=-xzero*gp_NOVA(j)*bzst(i,j)
      cyst(i,j)=bxdzst(i,j)-bzdxst(i,j)
!      cyst(i,j)=xxst(i,j)*pp_NOVA(j)+xzero**2*g_NOVA(j)*gp_NOVA(j)/xxst(i,j)

      uyst(i,j)=omrot_NOVA(j)*xxst(i,j)
      uydxst(i,j)=-omprot_NOVA(j)*bzst(i,j)*xxst(i,j)**2+omrot_NOVA(j)
      uydzst(i,j)= omprot_NOVA(j)*bxst(i,j)*xxst(i,j)**2

      if(j.ge.2 .and. i.ge.3 .and. i.le.n2th+2 .and. i.ne.ipi) then
      if(i.lt.ipi) jd=(j-2)*(n2th-1)+i-2
      if(i.gt.ipi) jd=(j-2)*(n2th-1)+i-3
      by_NOVA(jd)=byst(i,j)
      bydx_NOVA(jd)=bydxst(i,j)
      bydz_NOVA(jd)=bydzst(i,j)

      cy_NOVA(jd)=cyst(i,j)
      cx_NOVA(jd)=cxst(i,j)
      cz_NOVA(jd)=czst(i,j)

      uy_NOVA(jd)=uyst(i,j)
      uydx_NOVA(jd)=uydxst(i,j)
      uydz_NOVA(jd)=uydzst(i,j)
      endif
      enddo
      enddo

      if(nrank.eq.0) then
      do j=1,npsi
      do i=3,nthe+1
      fffst(i,j,1)=cyst(i,j)*bzst(i,j)-czst(i,j)*byst(i,j)-pdxst(i,j)-uyst(i,j)*uydxst(i,j)+uyst(i,j)**2/xxst(i,j)
      fffst(i,j,2)=czst(i,j)*bxst(i,j)-cxst(i,j)*bzst(i,j)
      fffst(i,j,3)=cxst(i,j)*byst(i,j)-cyst(i,j)*bxst(i,j)-pdzst(i,j)-uyst(i,j)*uydzst(i,j)
      enddo
      enddo

      open(unit=101,file='fst_NOVA.dat',status='unknown',form='formatted')
      write(101,300)(((fffst(i,j,m),m=1,3),i=3,nthe+1),j=2,npsi)
 300  format(3(1x,e12.5))

      endif

       psmin=minval(ps_NOVA)
       psmax=maxval(ps_NOVA)

       qmin=minval(q_NOVA)
       qmax=maxval(q_NOVA)
       q0=q_NOVA(1)
!      xmin=minval(xx_NOVA)
!      xmax=maxval(xx_NOVA)
!      zmin=minval(zz_NOVA)
!      zmax=maxval(zz_NOVA)
        epsilon=1./aratio
     
      do j=2,npsi
      do i=1,n2th+5
        tpst(i,j)=atan2(bxst(i,j),-bzst(i,j))
        if(i .gt. ipi) tpst(i,j)=tpst(i,j)+2*pi
        if(i .eq. ipi) tpst(i,j)=pi
        bpst(i,j)=sqrt(bxst(i,j)**2+bzst(i,j)**2)
        wstx2r(i,j)=-bzst(i,j)/bpst(i,j)
        wstz2r(i,j)=bxst(i,j)/bpst(i,j)
        wstx2p(i,j)=-bxst(i,j)/bpst(i,j)
        wstz2p(i,j)=-bzst(i,j)/bpst(i,j)
      enddo
      enddo




      if(nrank.eq.0) then
      open(890,file='psdat.dat')
      write(890,5000) (xx_NOVA(j),zz_NOVA(j),th_NOVA(j),ps_NOVA(j),bx_NOVA(j),bxdx_NOVA(j),bxdz_NOVA(j),bz_NOVA(j),bzdx_NOVA(j),bzdz_NOVA(j),j=1,ndat)
 5000 format(10(1x,e17.9))
      close(890)
      open(891,file='stdat.dat')
      write(891,5100) (((thst(i),psst(i,j),xxst(i,j),zzst(i,j),tpst(i,j),tst(i,j),rst(i,j)),i=1,n2th+5),j=1,npsi)
 5100 format(7(1x,e17.9))
      close(891)

      endif 
       
 1500 format(1x,i5,i5,15(1x,e17.9))
 2000 format(1x,i5,13(1x,e17.9))
 2100 format(1x,i5,11(1x,e17.9))
 4000 format(14(1x,e17.9))         
      return
      end

!ws*******************************************************************************
!ws************************************************************************
      subroutine recrd1
      USE DECLARE
      include 'mpif.h'
!
      character*12 output
      character*3 cn
      character*3 cn1
      
      output='per001_'//cn1(nrank)
      open(unit=7,file=output,status='unknown',form='unformatted')
      write(7)ncase,nstep,time
      write(7)x
      close(7)
    
      return
      end

      subroutine recrd10
      USE DECLARE
      include 'mpif.h'
!
      character*12 output
      character*3 cn
      character*3 cn1
      
      output='per010_'//cn1(nrank)
      open(unit=7,file=output,status='unknown',form='unformatted')
      write(7)ncase,nstep,time
      write(7)x
      close(7)
    
      return
      end
      
      subroutine recrd100
      USE DECLARE
      include 'mpif.h'
!
      character*12 output
      character*3 cn
      character*3 cn1
      
      output='per100_'//cn1(nrank)
      open(unit=7,file=output,status='unknown',form='unformatted')
      write(7)ncase,nstep,time
      write(7)x
      close(7)
    
      return
      end
!ws************************************************************************
      subroutine recrd
      USE DECLARE
      include 'mpif.h'
!
      character*9 output
      character*3 cn
      character*3 cn1
      
      output='tk'//cn1(nrank)//cn(nst)
      open(unit=7,file=output,status='unknown',form='unformatted')
      write(7)ncase,nstep,time
      write(7)x
      close(7)
    
      return
      end
!
     subroutine recrd_init
      USE DECLARE
      include 'mpif.h'
!
      character*8 output
      character*3 cn
      character*3 cn1

      output='int'//cn1(nrank)     
      open(unit=99,file=output,status='unknown',form='unformatted')
      write(99)ncase,nstep,time
      write(99)xint,cint
      close(99)
      return
      end  

      subroutine recrd_dssp
      USE DECLARE
      include 'mpif.h'
!
      character*8 output
      character*3 cn
      character*3 cn1

      output='dssp'//cn1(nrank)     
      open(unit=98,file=output,status='unknown',form='unformatted')
      write(98)ncase,nstep,time
      write(98)eta,fmu,pmu !,kap_pal,kap_perp
      close(98)
      return
      end 

      subroutine recrd_dbg
      USE DECLARE
      include 'mpif.h'
!
      character*12 output
      character*3 cn
      character*3 cn1
      
      do 10 jy=iy_first,iy_last
      do 10 jz=iz_first,iz_last
      do 10 jx=ix_first,ix_last 
      do m=1,8
      x1(jx,jz,jy,m)=x(jx,jz,jy,m)-xint(jx,jz,m)
      enddo
      do m=1,3
      xh(jx,jz,jy,m)=cur(jx,jz,jy,m)
      enddo
   10 continue
       
!      call vertex_x1_2
!      call vertex_x1_0_b1
!      call vertex_x1_0_1
      do 20 jy=iy_first,iy_last
      do 20 jz=iz_first,iz_last
      do 20 jx=ix_first,ix_last          
      vr(jx,jz,jy)=x1(jx,jz,jy,3)*wx2r(jx,jz)+x1(jx,jz,jy,5)*wz2r(jx,jz)
      vp(jx,jz,jy)=x1(jx,jz,jy,3)*wx2p(jx,jz)+x1(jx,jz,jy,5)*wz2p(jx,jz)
      br(jx,jz,jy)=x1(jx,jz,jy,6)*wx2r(jx,jz)+x1(jx,jz,jy,8)*wz2r(jx,jz)
      bp(jx,jz,jy)=x1(jx,jz,jy,6)*wx2p(jx,jz)+x1(jx,jz,jy,8)*wz2p(jx,jz)
!       
      cr(jx,jz,jy)=xh(jx,jz,jy,1)*wx2r(jx,jz)+xh(jx,jz,jy,3)*wz2r(jx,jz)
      cp(jx,jz,jy)=xh(jx,jz,jy,1)*wx2p(jx,jz)+xh(jx,jz,jy,3)*wz2p(jx,jz)
   20 continue

      output='x'//cn1(nrank)//cn(nst)
      open(unit=7,file=output,status='unknown',form='formatted') 
      jy=1   
      write(7,200)((((x1(jx,jz,jy,i),i=1,8),br(jx,jz,jy),bp(jx,jz,jy),vr(jx,jz,jy),vp(jx,jz,jy),vr(jx,jz,jy)*xx(jx)*bp0(jx,jz)),jx=ix_first,ix_last),jz=iz_first,iz_last)
 200  format(13(1x,e12.5)) 
      close(7)

      output='xy_x'//cn1(nrank)//cn(nst)
      open(unit=17,file=output,status='unknown',form='formatted') 
      jz=mz/2   
      write(17,200)((((x1(jx,jz,jy,i),i=1,8),br(jx,jz,jy),bp(jx,jz,jy),vr(jx,jz,jy),vp(jx,jz,jy),vr(jx,jz,jy)*xx(jx)*bp0(jx,jz)),jx=ix_first,ix_last),jy=1,my)
      close(17)
!      output='xrt'//cn1(nrank)//cn(nst)
!      open(unit=8,file=output,status='unknown',form='formatted')
!      write(8,300)(((xrt(jr,jt,1,i)-xrtint(jr,jt,i),i=1,8),jr=1,mr),jt=1,mt)
! 300  format(8(1x,e12.5)) 
!      close(8)
      output='cur'//cn1(nrank)//cn(nst)
      open(unit=9,file=output,status='unknown',form='formatted')
      write(9,400)((((cur(jx,jz,1,i),i=1,3),cr(jx,jz,1),cp(jx,jz,1)),jx=ix_first,ix_last),jz=iz_first,iz_last)
 400  format(5(1x,e12.5)) 
      close(9)
      output='xy_c'//cn1(nrank)//cn(nst)
      open(unit=19,file=output,status='unknown',form='formatted')
      jz=mz/2
      write(19,400)((((cur(jx,jz,jy,i),i=1,3),cr(jx,jz,jy),cp(jx,jz,jy)),jx=ix_first,ix_last),jy=1,my)
      close(19)
      output='Ef'//cn1(nrank)//cn(nst)
      open(unit=10,file=output,status='unknown',form='formatted')
      write(10,300)(((Ef(jx,jz,1,i),i=1,3),jx=ix_first,ix_last),jz=iz_first,iz_last)
 300  format(3(1x,e12.5)) 
      close(10)

!      call recrd_x1_a
!
!      output='x_dx'//cn1(nrank)//cn(nst)
!      open(unit=72,file=output,status='unknown',form='formatted')
!      write(72,300)(((xr(jx,jz,1,i),i=1,8),jx=ix_first,ix_last),jz=iz_first,iz_last) 
!      close(72)
!      output='x_dz'//cn1(nrank)//cn(nst)
!      open(unit=73,file=output,status='unknown',form='formatted')
!      write(73,300)(((xz(jx,jz,1,i),i=1,8),jx=ix_first,ix_last),jz=iz_first,iz_last) 
!      close(73)
!      output='x_dy'//cn1(nrank)//cn(nst)
!      open(unit=74,file=output,status='unknown',form='formatted')
!      write(74,300)(((xz(jx,jz,1,i),i=1,8),jx=ix_first,ix_last),jz=iz_first,iz_last) 
!      close(74)      
!      output='xrt_dr'//cn1(nrank)//cn(nst)
!      open(unit=82,file=output,status='unknown',form='formatted')
!      write(82,300)(((d1xrt1_dr(jr,jt,1,i),i=1,8),jr=1,mr),jt=1,mt)
!      close(82)
!      output='xrt_dt'//cn1(nrank)//cn(nst)
!      open(unit=83,file=output,status='unknown',form='formatted')
!      write(83,300)(((d1xrt1_dt(jr,jt,1,i),i=1,8),jr=1,mr),jt=1,mt)
!      close(83)
!       output='xrt_dy'//cn1(nrank)//cn(nst)
!      open(unit=84,file=output,status='unknown',form='formatted')
!      write(84,300)(((d1xrt_dy(jr,jt,1,i),i=1,8),jr=1,mr),jt=1,mt)
!      close(84)
      return
      end
!ws**************************************************

      subroutine recrd_per
      USE DECLARE
      include 'mpif.h'
!
      character*12 output
      character*3 cn
      character*3 cn1
      
      do 10 jy=iy_first,iy_last
      do 10 jz=iz_first,iz_last
      do 10 jx=ix_first,ix_last 
      do m=1,8
      x1(jx,jz,jy,m)=x(jx,jz,jy,m)-xint(jx,jz,m)
      enddo
      do m=1,3
      xh(jx,jz,jy,m)=cur(jx,jz,jy,m)
      enddo
   10 continue
       
!      call vertex_x1_2
!      call vertex_x1_0_b1
!      call vertex_x1_0_1
      do 20 jy=iy_first,iy_last
      do 20 jz=iz_first,iz_last
      do 20 jx=ix_first,ix_last          
      vr(jx,jz,jy)=x1(jx,jz,jy,3)*wx2r(jx,jz)+x1(jx,jz,jy,5)*wz2r(jx,jz)
      vp(jx,jz,jy)=x1(jx,jz,jy,3)*wx2p(jx,jz)+x1(jx,jz,jy,5)*wz2p(jx,jz)
      br(jx,jz,jy)=x1(jx,jz,jy,6)*wx2r(jx,jz)+x1(jx,jz,jy,8)*wz2r(jx,jz)
      bp(jx,jz,jy)=x1(jx,jz,jy,6)*wx2p(jx,jz)+x1(jx,jz,jy,8)*wz2p(jx,jz)
!       
      cr(jx,jz,jy)=xh(jx,jz,jy,1)*wx2r(jx,jz)+xh(jx,jz,jy,3)*wz2r(jx,jz)
      cp(jx,jz,jy)=xh(jx,jz,jy,1)*wx2p(jx,jz)+xh(jx,jz,jy,3)*wz2p(jx,jz)
   20 continue

      output='per_x'
      open(unit=7,file=output,status='unknown',form='formatted')   
      write(7,200)(((((x1(jx,jz,jy,i),i=1,8),br(jx,jz,jy),bp(jx,jz,jy),vr(jx,jz,jy),vp(jx,jz,jy),vr(jx,jz,jy)*xx(jx)*bp0(jx,jz)),jx=ix_first,ix_last),jz=iz_first,iz_last),jy=1,my)
 200  format(13(1x,e12.5)) 
      close(7)

      output='per_cur'
      open(unit=9,file=output,status='unknown',form='formatted')
      write(9,400)(((((cur(jx,jz,jy,i),i=1,3),cr(jx,jz,jy),cp(jx,jz,jy)),jx=ix_first,ix_last),jz=iz_first,iz_last),jy=1,my)
 400  format(5(1x,e12.5)) 
      close(9)
      return
      end

!ws***********************************************
     subroutine recrd_xdif
      USE DECLARE
      include 'mpif.h'
      character*12 output
      character*3 cn
      character*3 cn1
      output='xdif'//cn(nst)//cn1(irk)
      open(unit=71,file=output,status='unknown',form='formatted')
      write(71,200) (((( xdif(jx,jz,jy,i),i=1,8),jx=ix_first,ix_last),jz=iz_first,iz_last),jy=1,my)
 200  format(8(1x,e12.5)) 
      close(71)
!      output='xrtdif'//cn(nst)//cn1(irk)
!      open(unit=81,file=output,status='unknown',form='formatted')
!      write(81,300)(((xrtdif(jr,jt,1,i),i=1,8),jr=1,mr),jt=1,mt)
! 300  format(8(1x,e12.5))
!      close(81)
!      output='Efd'//cn(nst)//cn1(irk)
!      open(unit=72,file=output,status='unknown',form='formatted')
!      write(72,200)((((Ef(jx,jz,1,i),Efx(jx,jz,1,i),Efz(jx,jz,1,i)),i=1,3),jx=ix_first,ix_last),jz=iz_first,iz_last)
!      close(72)
      return
      end
!*********************************************************
     subroutine recrd_cv
      USE DECLARE
      include 'mpif.h'
      character*12 output
      character*3 cn
      character*3 cn1
      output='xdx'//cn(nst)//cn1(irk)
      open(unit=71,file=output,status='unknown',form='formatted')
      write(71,200) ((( x1r(jx,jz,1,i),i=1,8),(xr2(jx,jz,1,i),i=1,8),jx=ix_first,ix_last),jz=iz_first,iz_last)
 200  format(16(1x,e12.5)) 
      close(71)
      output='xdz'//cn(nst)//cn1(irk)
      open(unit=72,file=output,status='unknown',form='formatted')
      write(72,200) ((( x1z(jx,jz,1,i),i=1,8),(xz2(jx,jz,1,i),i=1,8),jx=ix_first,ix_last),jz=iz_first,iz_last)
      close(72)
      output='xdy'//cn(nst)//cn1(irk)
      open(unit=73,file=output,status='unknown',form='formatted')
      write(73,200) ((( xy(jx,jz,1,i),i=1,8),(xy2(jx,jz,1,i),i=1,8),jx=ix_first,ix_last),jz=iz_first,iz_last)
      close(73)

      jz=mz/2
      output='xy_xdx'//cn(nst)//cn1(irk)
      open(unit=74,file=output,status='unknown',form='formatted')
      write(74,200) ((( x1r(jx,jz,jy,i),i=1,8),(xr2(jx,jz,jy,i),i=1,8),jx=ix_first,ix_last),jy=1,my)
      close(74)
      output='xy_xdz'//cn(nst)//cn1(irk)
      open(unit=75,file=output,status='unknown',form='formatted')
      write(75,200) ((( x1z(jx,jz,jy,i),i=1,8),(xz2(jx,jz,jy,i),i=1,8),jx=ix_first,ix_last),jy=1,my)
      close(75)
      output='xy_xdy'//cn(nst)//cn1(irk)
      open(unit=76,file=output,status='unknown',form='formatted')
      write(76,200) ((( xy(jx,jz,jy,i),i=1,8),(xy2(jx,jz,jy,i),i=1,8),jx=ix_first,ix_last),jy=1,my)
      close(76)

      return
      end
 
!********************************************* 
     subroutine recrd_x1st
      USE DECLARE
      include 'mpif.h'
      character*14 output
      character*3 cn
      character*3 cn1
      output='x1st'//cn(nst)
      open(unit=71,file=output,status='unknown',form='formatted')
      write(71,200)(((x1st(jt,js,1,i),i=1,8),js=mpsa,mpsa-nda,-1),jt=1,n2th+5)
 200  format(8(1x,e12.5)) 
      close(71)
      output='byst'//cn(nst)
      open(unit=77,file=output,status='unknown',form='formatted')
      write(77,400)((x1st(jt,js,1,7),js=mpsa,mpsa-nda,-1),jt=1,n2th+5)
 400  format(5(1x,e12.5)) 
      close(77)
      
      output='cxst'//cn(nst)
      open(unit=81,file=output,status='unknown',form='formatted')
      write(81,400)((xhst(jt,js,1,1),js=mpsa,mpsa-nda,-1),jt=1,n2th+5)
      close(81)
      return
      end
!ws*************************************        
      subroutine readin
      USE DECLARE
      include 'mpif.h'
!
      character*9 output
      character*3 cn
      character*3 cn1
      
      output='tk'//cn1(nrank)//cn(nst)
      open(unit=7,file=output,status='unknown',form='unformatted')
      read(7)ncase,nstep,time
      read(7)x
      close(7)
      return
      end
!ws*************************************

      character*3 function cn(n)
!
!-----assume that n is no greater than 999
!
!
!-----separate the digits
!
      n1=n/100
      n2=(n-100*n1)/10
      n3=n-100*n1-10*n2
!
!-----stick together cn using char function
!
      n1=n1+48
      n2=n2+48
      n3=n3+48
      cn(1:1)=char(n1)
      cn(2:2)=char(n2)
      cn(3:3)=char(n3)
!
      return
      end
!
      character*3 function cn1(n)
!
!-----assume that n is no greater than 999
!
!
!-----separate the digits
!
      n1=n/100
      n2=(n-100*n1)/10
      n3=n-100*n1-10*n2
!
!-----stick together cn using char function
!
      n1=n1+48
      n2=n2+48
      n3=n3+48
      cn1(1:1)=char(n1)
      cn1(2:2)=char(n2)
      cn1(3:3)=char(n3)
!
      return
      end
!ws*************************************************
!ws*************************************************
!****************************************************************************************
      subroutine mpi_transfersm(ws,mm)
      USE DECLARE
      real*8, dimension(mx,mz,my,mm) :: ws
      real*8, dimension(mz,my,mm) :: wsx1,wsx2
      real*8, dimension(mx,my,mm) :: wsz1,wsz2
      real*8, dimension(mx,mz,mm) :: wsy1,wsy2
      include 'mpif.h'

!       
! Send w8 up unless I'm at the top, then receive from below
     
	if (nrankxz.ge.nrkz(nrank)*nprx .and. nrankxz.lt.nrkz(nrank)*nprx+nprx-1) then
      wsx1(:,:,:)=ws(ix_last-2,:,:,:)
      wsx2(:,:,:)=ws(ix_last-3,:,:,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( wsx1, myz*mm, MPI_DOUBLE_PRECISION, nrank + 1, 0,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wsx2, myz*mm, MPI_DOUBLE_PRECISION, nrank + 1, 0,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrankxz.gt.nrkz(nrank)*nprx .and. nrankxz.le.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wsx1, myz*mm, MPI_DOUBLE_PRECISION, nrank - 1, 0,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wsx2, myz*mm, MPI_DOUBLE_PRECISION, nrank - 1, 0,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      ws(ix_first+1,:,:,:)=wsx1(:,:,:)
      ws(ix_first,:,:,:)=wsx2(:,:,:)
	endif
	
      
! Send w8 down unless I'm at the bottom

      
	if (nrankxz.gt.nrkz(nrank)*nprx .and. nrankxz.le.nrkz(nrank)*nprx+nprx-1) then
      wsx1(:,:,:)=ws(ix_first+2,:,:,:)
      wsx2(:,:,:)=ws(ix_first+3,:,:,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( wsx1, myz*mm, MPI_DOUBLE_PRECISION, nrank - 1, 1,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wsx2, myz*mm, MPI_DOUBLE_PRECISION, nrank - 1, 1,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrankxz.ge.nrkz(nrank)*nprx .and. nrankxz.lt.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wsx1, myz*mm, MPI_DOUBLE_PRECISION, nrank + 1, 1,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wsx2, myz*mm, MPI_DOUBLE_PRECISION, nrank + 1, 1,  &
		      MPI_COMM_WORLD, status,ierror )
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
	CALL MPI_Send( wsz1, myx*mm, MPI_DOUBLE_PRECISION, nrank + nprx, 0,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wsz2, myx*mm, MPI_DOUBLE_PRECISION, nrank + nprx, 0,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrankxz.ge.nprx ) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wsz1, myx*mm, MPI_DOUBLE_PRECISION, nrank - nprx, 0,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wsz2, myx*mm, MPI_DOUBLE_PRECISION, nrank - nprx, 0,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      ws(:,iz_first+1,:,:)=wsz1(:,:,:)
      ws(:,iz_first,:,:)=wsz2(:,:,:)
	endif
	     
! Send w8 down unless I'm at the bottom

      
	if (nrankxz.ge.nprx ) then
      wsz1(:,:,:)=ws(:,iz_first+2,:,:)
      wsz2(:,:,:)=ws(:,iz_first+3,:,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( wsz1, myx*mm, MPI_DOUBLE_PRECISION, nrank - nprx, 1,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wsz2, myx*mm, MPI_DOUBLE_PRECISION, nrank - nprx, 1,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrankxz.lt.nsizexz-nprx) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wsz1, myx*mm, MPI_DOUBLE_PRECISION, nrank + nprx, 1,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wsz2, myx*mm, MPI_DOUBLE_PRECISION, nrank + nprx, 1,  &
		      MPI_COMM_WORLD, status,ierror )
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
	CALL MPI_Send( wsy1, mxz*mm, MPI_DOUBLE_PRECISION, nrank + nprxz, 0,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wsy2, mxz*mm, MPI_DOUBLE_PRECISION, nrank + nprxz, 0,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrky(nrank) .ge. 1 )  then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wsy1, mxz*mm, MPI_DOUBLE_PRECISION, nrank - nprxz, 0,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wsy2, mxz*mm, MPI_DOUBLE_PRECISION, nrank - nprxz, 0,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      ws(:,:,iy_first+1,:)=wsy1(:,:,:)
      ws(:,:,iy_first,:)=wsy2(:,:,:)
	endif


  	if (nrky(nrank).eq. npry-1) then
      wsy1(:,:,:)=ws(:,:,iy_last-2,:)
      wsy2(:,:,:)=ws(:,:,iy_last-3,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( wsy1, mxz*mm, MPI_DOUBLE_PRECISION, nrank -(npry-1)*nprxz, 0,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wsy2, mxz*mm, MPI_DOUBLE_PRECISION, nrank -(npry-1)*nprxz, 0,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrky(nrank) .eq. 0 )  then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wsy1, mxz*mm, MPI_DOUBLE_PRECISION, nrank +(npry-1)*nprxz, 0,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wsy2, mxz*mm, MPI_DOUBLE_PRECISION, nrank +(npry-1)*nprxz, 0,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      ws(:,:,iy_first+1,:)=wsy1(:,:,:)
      ws(:,:,iy_first,:)=wsy2(:,:,:)
	endif
	     
! Send w8 down unless I'm at the bottom

      
	if (nrky(nrank) .ge. 1 ) then
      wsy1(:,:,:)=ws(:,:,iy_first+2,:)
      wsy2(:,:,:)=ws(:,:,iy_first+3,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( wsy1, mxz*mm, MPI_DOUBLE_PRECISION, nrank - nprxz, 1,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wsy2, mxz*mm, MPI_DOUBLE_PRECISION, nrank - nprxz, 1,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrky(nrank).lt. npry-1) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wsy1, mxz*mm, MPI_DOUBLE_PRECISION, nrank + nprxz, 1,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wsy2, mxz*mm, MPI_DOUBLE_PRECISION, nrank + nprxz, 1,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      ws(:,:,iy_last-1,:)=wsy1(:,:,:)
      ws(:,:,iy_last,:)=wsy2(:,:,:)
	endif

    if (nrky(nrank) .eq. 0 ) then
      wsy1(:,:,:)=ws(:,:,iy_first+2,:)
      wsy2(:,:,:)=ws(:,:,iy_first+3,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( wsy1, mxz*mm, MPI_DOUBLE_PRECISION, nrank +(npry-1)*nprxz, 1,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wsy2, mxz*mm, MPI_DOUBLE_PRECISION, nrank +(npry-1)*nprxz, 1,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrky(nrank).eq. npry-1) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wsy1, mxz*mm, MPI_DOUBLE_PRECISION, nrank -(npry-1)*nprxz, 1,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wsy2, mxz*mm, MPI_DOUBLE_PRECISION, nrank -(npry-1)*nprxz, 1,  &
		      MPI_COMM_WORLD, status,ierror )
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


!****************************************************************************************
      subroutine mpi_transfer8(w8)
      USE DECLARE
      real*8, dimension(mx,mz,my,8) :: w8
      real*8, dimension(mz,my,8) :: wfx1,wfx2
      real*8, dimension(mx,my,8) :: wfz1,wfz2
      real*8, dimension(mx,mz,8) :: wfy1,wfy2
      include 'mpif.h'

!       
! Send w8 up unless I'm at the top, then receive from below
     
	if (nrankxz.ge.nrkz(nrank)*nprx .and. nrankxz.lt.nrkz(nrank)*nprx+nprx-1) then
      wfx1(:,:,:)=w8(ix_last-2,:,:,:)
      wfx2(:,:,:)=w8(ix_last-3,:,:,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( wfx1, myz8, MPI_DOUBLE_PRECISION, nrank + 1, 0,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wfx2, myz8, MPI_DOUBLE_PRECISION, nrank + 1, 0,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrankxz.gt.nrkz(nrank)*nprx .and. nrankxz.le.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wfx1, myz8, MPI_DOUBLE_PRECISION, nrank - 1, 0,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wfx2, myz8, MPI_DOUBLE_PRECISION, nrank - 1, 0,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      w8(ix_first+1,:,:,:)=wfx1(:,:,:)
      w8(ix_first,:,:,:)=wfx2(:,:,:)
	endif
	
      
! Send w8 down unless I'm at the bottom

      
	if (nrankxz.gt.nrkz(nrank)*nprx .and. nrankxz.le.nrkz(nrank)*nprx+nprx-1) then
      wfx1(:,:,:)=w8(ix_first+2,:,:,:)
      wfx2(:,:,:)=w8(ix_first+3,:,:,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( wfx1, myz8, MPI_DOUBLE_PRECISION, nrank - 1, 1,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wfx2, myz8, MPI_DOUBLE_PRECISION, nrank - 1, 1,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrankxz.ge.nrkz(nrank)*nprx .and. nrankxz.lt.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wfx1, myz8, MPI_DOUBLE_PRECISION, nrank + 1, 1,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wfx2, myz8, MPI_DOUBLE_PRECISION, nrank + 1, 1,  &
		      MPI_COMM_WORLD, status,ierror )
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
	CALL MPI_Send( wfz1, myx8, MPI_DOUBLE_PRECISION, nrank + nprx, 0,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wfz2, myx8, MPI_DOUBLE_PRECISION, nrank + nprx, 0,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrankxz.ge.nprx ) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wfz1, myx8, MPI_DOUBLE_PRECISION, nrank - nprx, 0,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wfz2, myx8, MPI_DOUBLE_PRECISION, nrank - nprx, 0,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      w8(:,iz_first+1,:,:)=wfz1(:,:,:)
      w8(:,iz_first,:,:)=wfz2(:,:,:)
	endif
	     
! Send w8 down unless I'm at the bottom

      
	if (nrankxz.ge.nprx ) then
      wfz1(:,:,:)=w8(:,iz_first+2,:,:)
      wfz2(:,:,:)=w8(:,iz_first+3,:,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( wfz1, myx8, MPI_DOUBLE_PRECISION, nrank - nprx, 1,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wfz2, myx8, MPI_DOUBLE_PRECISION, nrank - nprx, 1,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrankxz.lt.nsizexz-nprx) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wfz1, myx8, MPI_DOUBLE_PRECISION, nrank + nprx, 1,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wfz2, myx8, MPI_DOUBLE_PRECISION, nrank + nprx, 1,  &
		      MPI_COMM_WORLD, status,ierror )
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
	CALL MPI_Send( wfy1, mxz8, MPI_DOUBLE_PRECISION, nrank + nprxz, 0,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wfy2, mxz8, MPI_DOUBLE_PRECISION, nrank + nprxz, 0,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrky(nrank) .ge. 1 )  then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wfy1, mxz8, MPI_DOUBLE_PRECISION, nrank - nprxz, 0,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wfy2, mxz8, MPI_DOUBLE_PRECISION, nrank - nprxz, 0,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      w8(:,:,iy_first+1,:)=wfy1(:,:,:)
      w8(:,:,iy_first,:)=wfy2(:,:,:)
	endif


  	if (nrky(nrank).eq. npry-1) then
      wfy1(:,:,:)=w8(:,:,iy_last-2,:)
      wfy2(:,:,:)=w8(:,:,iy_last-3,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( wfy1, mxz8, MPI_DOUBLE_PRECISION, nrank -(npry-1)*nprxz, 0,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wfy2, mxz8, MPI_DOUBLE_PRECISION, nrank -(npry-1)*nprxz, 0,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrky(nrank) .eq. 0 )  then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wfy1, mxz8, MPI_DOUBLE_PRECISION, nrank +(npry-1)*nprxz, 0,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wfy2, mxz8, MPI_DOUBLE_PRECISION, nrank +(npry-1)*nprxz, 0,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      w8(:,:,iy_first+1,:)=wfy1(:,:,:)
      w8(:,:,iy_first,:)=wfy2(:,:,:)
	endif
	     
! Send w8 down unless I'm at the bottom

      
	if (nrky(nrank) .ge. 1 ) then
      wfy1(:,:,:)=w8(:,:,iy_first+2,:)
      wfy2(:,:,:)=w8(:,:,iy_first+3,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( wfy1, mxz8, MPI_DOUBLE_PRECISION, nrank - nprxz, 1,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wfy2, mxz8, MPI_DOUBLE_PRECISION, nrank - nprxz, 1,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrky(nrank).lt. npry-1) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wfy1, mxz8, MPI_DOUBLE_PRECISION, nrank + nprxz, 1,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wfy2, mxz8, MPI_DOUBLE_PRECISION, nrank + nprxz, 1,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      w8(:,:,iy_last-1,:)=wfy1(:,:,:)
      w8(:,:,iy_last,:)=wfy2(:,:,:)
	endif

    if (nrky(nrank) .eq. 0 ) then
      wfy1(:,:,:)=w8(:,:,iy_first+2,:)
      wfy2(:,:,:)=w8(:,:,iy_first+3,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( wfy1, mxz8, MPI_DOUBLE_PRECISION, nrank +(npry-1)*nprxz, 1,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wfy2, mxz8, MPI_DOUBLE_PRECISION, nrank +(npry-1)*nprxz, 1,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrky(nrank).eq. npry-1) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wfy1, mxz8, MPI_DOUBLE_PRECISION, nrank -(npry-1)*nprxz, 1,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wfy2, mxz8, MPI_DOUBLE_PRECISION, nrank -(npry-1)*nprxz, 1,  &
		      MPI_COMM_WORLD, status,ierror )
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
!****************************************************************************************
      subroutine mpi_transfer3(w3)
      USE DECLARE
      real*8, dimension(mx,mz,my,3) :: w3
      real*8, dimension(mz,my,3) :: wcx1,wcx2
      real*8, dimension(mx,my,3) :: wcz1,wcz2
      real*8, dimension(mx,mz,3) :: wcy1,wcy2
      include 'mpif.h'

       
! Send w3 up unless I'm at the top, then receive from below
     
	if (nrankxz.ge.nrkz(nrank)*nprx .and. nrankxz.lt.nrkz(nrank)*nprx+nprx-1) then
      wcx1(:,:,:)=w3(ix_last-2,:,:,:)
      wcx2(:,:,:)=w3(ix_last-3,:,:,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( wcx1, myz3, MPI_DOUBLE_PRECISION, nrank + 1, 0,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wcx2, myz3, MPI_DOUBLE_PRECISION, nrank + 1, 0,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrankxz.gt.nrkz(nrank)*nprx .and. nrankxz.le.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wcx1, myz3, MPI_DOUBLE_PRECISION, nrank - 1, 0,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wcx2, myz3, MPI_DOUBLE_PRECISION, nrank - 1, 0,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      w3(ix_first+1,:,:,:)=wcx1(:,:,:)
      w3(ix_first,:,:,:)=wcx2(:,:,:)
	endif
	
      
! Send w3 down unless I'm at the bottom

      
	if (nrankxz.gt.nrkz(nrank)*nprx .and. nrankxz.le.nrkz(nrank)*nprx+nprx-1) then
      wcx1(:,:,:)=w3(ix_first+2,:,:,:)
      wcx2(:,:,:)=w3(ix_first+3,:,:,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( wcx1, myz3, MPI_DOUBLE_PRECISION, nrank - 1, 1,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wcx2, myz3, MPI_DOUBLE_PRECISION, nrank - 1, 1,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrankxz.ge.nrkz(nrank)*nprx .and. nrankxz.lt.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wcx1, myz3, MPI_DOUBLE_PRECISION, nrank + 1, 1,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wcx2, myz3, MPI_DOUBLE_PRECISION, nrank + 1, 1,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      w3(ix_last-1,:,:,:)=wcx1(:,:,:)
      w3(ix_last,:,:,:)=wcx2(:,:,:)
	endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	if (nrankxz.lt.nsizexz-nprx) then
      wcz1(:,:,:)=w3(:,iz_last-2,:,:)
      wcz2(:,:,:)=w3(:,iz_last-3,:,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( wcz1, myx3, MPI_DOUBLE_PRECISION, nrank + nprx, 0,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wcz2, myx3, MPI_DOUBLE_PRECISION, nrank + nprx, 0,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrankxz.ge.nprx ) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wcz1, myx3, MPI_DOUBLE_PRECISION, nrank - nprx, 0,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wcz2, myx3, MPI_DOUBLE_PRECISION, nrank - nprx, 0,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      w3(:,iz_first+1,:,:)=wcz1(:,:,:)
      w3(:,iz_first,:,:)=wcz2(:,:,:)
	endif
	     
! Send w3 down unless I'm at the bottom

      
	if (nrankxz.ge.nprx ) then
      wcz1(:,:,:)=w3(:,iz_first+2,:,:)
      wcz2(:,:,:)=w3(:,iz_first+3,:,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( wcz1, myx3, MPI_DOUBLE_PRECISION, nrank - nprx, 1,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wcz2, myx3, MPI_DOUBLE_PRECISION, nrank - nprx, 1,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrankxz.lt.nsizexz-nprx) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wcz1, myx3, MPI_DOUBLE_PRECISION, nrank + nprx, 1,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wcz2, myx3, MPI_DOUBLE_PRECISION, nrank + nprx, 1,  &
		      MPI_COMM_WORLD, status,ierror )
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
	CALL MPI_Send( wcy1, mxz3, MPI_DOUBLE_PRECISION, nrank + nprxz, 0,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wcy2, mxz3, MPI_DOUBLE_PRECISION, nrank + nprxz, 0,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrky(nrank) .ge. 1 )  then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wcy1, mxz3, MPI_DOUBLE_PRECISION, nrank - nprxz, 0,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wcy2, mxz3, MPI_DOUBLE_PRECISION, nrank - nprxz, 0,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      w3(:,:,iy_first+1,:)=wcy1(:,:,:)
      w3(:,:,iy_first,:)=wcy2(:,:,:)
	endif


  	if (nrky(nrank).eq. npry-1) then
      wcy1(:,:,:)=w3(:,:,iy_last-2,:)
      wcy2(:,:,:)=w3(:,:,iy_last-3,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( wcy1, mxz3, MPI_DOUBLE_PRECISION, nrank -(npry-1)*nprxz, 0,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wcy2, mxz3, MPI_DOUBLE_PRECISION, nrank -(npry-1)*nprxz, 0,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrky(nrank) .eq. 0 )  then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wcy1, mxz3, MPI_DOUBLE_PRECISION, nrank +(npry-1)*nprxz, 0,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wcy2, mxz3, MPI_DOUBLE_PRECISION, nrank +(npry-1)*nprxz, 0,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      w3(:,:,iy_first+1,:)=wcy1(:,:,:)
      w3(:,:,iy_first,:)=wcy2(:,:,:)
	endif
	     
! Send w3 down unless I'm at the bottom

      
	if (nrky(nrank) .ge. 1 ) then
      wcy1(:,:,:)=w3(:,:,iy_first+2,:)
      wcy2(:,:,:)=w3(:,:,iy_first+3,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( wcy1, mxz3, MPI_DOUBLE_PRECISION, nrank - nprxz, 1,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wcy2, mxz3, MPI_DOUBLE_PRECISION, nrank - nprxz, 1,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrky(nrank).lt. npry-1) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wcy1, mxz3, MPI_DOUBLE_PRECISION, nrank + nprxz, 1,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wcy2, mxz3, MPI_DOUBLE_PRECISION, nrank + nprxz, 1,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      w3(:,:,iy_last-1,:)=wcy1(:,:,:)
      w3(:,:,iy_last,:)=wcy2(:,:,:)
	endif

    if (nrky(nrank) .eq. 0 ) then
      wcy1(:,:,:)=w3(:,:,iy_first+2,:)
      wcy2(:,:,:)=w3(:,:,iy_first+3,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( wcy1, mxz3, MPI_DOUBLE_PRECISION, nrank +(npry-1)*nprxz, 1,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wcy2, mxz3, MPI_DOUBLE_PRECISION, nrank +(npry-1)*nprxz, 1,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrky(nrank).eq. npry-1) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wcy1, mxz3, MPI_DOUBLE_PRECISION, nrank -(npry-1)*nprxz, 1,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wcy2, mxz3, MPI_DOUBLE_PRECISION, nrank -(npry-1)*nprxz, 1,  &
		      MPI_COMM_WORLD, status,ierror )
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

!****************************************************************************************
      subroutine mpi_transfer1(w1)
      USE DECLARE
      real*8, dimension(mx,mz,my) :: w1
      real*8, dimension(mz,my) :: w1x1,w1x2
      real*8, dimension(mx,my) :: w1z1,w1z2
      real*8, dimension(mx,mz) :: w1y1,w1y2
      include 'mpif.h'

!       
! Send w1 up unless I'm at the top, then receive from below
     
	if (nrankxz.ge.nrkz(nrank)*nprx .and. nrankxz.lt.nrkz(nrank)*nprx+nprx-1) then
      w1x1(:,:)=w1(ix_last-2,:,:)
      w1x2(:,:)=w1(ix_last-3,:,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( w1x1, myz, MPI_DOUBLE_PRECISION, nrank + 1, 0,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( w1x2, myz, MPI_DOUBLE_PRECISION, nrank + 1, 0,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrankxz.gt.nrkz(nrank)*nprx .and. nrankxz.le.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( w1x1, myz, MPI_DOUBLE_PRECISION, nrank - 1, 0,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( w1x2, myz, MPI_DOUBLE_PRECISION, nrank - 1, 0,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      w1(ix_first+1,:,:)=w1x1(:,:)
      w1(ix_first,:,:)=w1x2(:,:)
	endif
	
      
! Send w1 down unless I'm at the bottom

      
	if (nrankxz.gt.nrkz(nrank)*nprx .and. nrankxz.le.nrkz(nrank)*nprx+nprx-1) then
      w1x1(:,:)=w1(ix_first+2,:,:)
      w1x2(:,:)=w1(ix_first+3,:,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( w1x1, myz, MPI_DOUBLE_PRECISION, nrank - 1, 1,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( w1x2, myz, MPI_DOUBLE_PRECISION, nrank - 1, 1,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrankxz.ge.nrkz(nrank)*nprx .and. nrankxz.lt.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( w1x1, myz, MPI_DOUBLE_PRECISION, nrank + 1, 1,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( w1x2, myz, MPI_DOUBLE_PRECISION, nrank + 1, 1,  &
		      MPI_COMM_WORLD, status,ierror )
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
	CALL MPI_Send( w1z1, myx, MPI_DOUBLE_PRECISION, nrank + nprx, 0,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( w1z2, myx, MPI_DOUBLE_PRECISION, nrank + nprx, 0,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrankxz.ge.nprx ) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( w1z1, myx, MPI_DOUBLE_PRECISION, nrank - nprx, 0,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( w1z2, myx, MPI_DOUBLE_PRECISION, nrank - nprx, 0,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      w1(:,iz_first+1,:)=w1z1(:,:)
      w1(:,iz_first,:)=w1z2(:,:)
	endif
	     
! Send w1 down unless I'm at the bottom

      
	if (nrankxz.ge.nprx ) then
      w1z1(:,:)=w1(:,iz_first+2,:)
      w1z2(:,:)=w1(:,iz_first+3,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( w1z1, myx, MPI_DOUBLE_PRECISION, nrank - nprx, 1,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( w1z2, myx, MPI_DOUBLE_PRECISION, nrank - nprx, 1,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrankxz.lt.nsizews-nprx) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( w1z1, myx, MPI_DOUBLE_PRECISION, nrank + nprx, 1,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( w1z2, myx, MPI_DOUBLE_PRECISION, nrank + nprx, 1,  &
		      MPI_COMM_WORLD, status,ierror )
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
	CALL MPI_Send( w1y1, mxz, MPI_DOUBLE_PRECISION, nrank + nprxz, 0,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( w1y2, mxz, MPI_DOUBLE_PRECISION, nrank + nprxz, 0,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrky(nrank) .ge. 1 )  then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( w1y1, mxz, MPI_DOUBLE_PRECISION, nrank - nprxz, 0,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( w1y2, mxz, MPI_DOUBLE_PRECISION, nrank - nprxz, 0,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      w1(:,:,iy_first+1)=w1y1(:,:)
      w1(:,:,iy_first)=w1y2(:,:)
	endif


  	if (nrky(nrank).eq. npry-1) then
      w1y1(:,:)=w1(:,:,iy_last-2)
      w1y2(:,:)=w1(:,:,iy_last-3)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( w1y1, mxz, MPI_DOUBLE_PRECISION, nrank -(npry-1)*nprxz, 0,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( w1y2, mxz, MPI_DOUBLE_PRECISION, nrank -(npry-1)*nprxz, 0,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrky(nrank) .eq. 0 )  then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( w1y1, mxz, MPI_DOUBLE_PRECISION, nrank +(npry-1)*nprxz, 0,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( w1y2, mxz, MPI_DOUBLE_PRECISION, nrank +(npry-1)*nprxz, 0,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      w1(:,:,iy_first+1)=w1y1(:,:)
      w1(:,:,iy_first)=w1y2(:,:)
	endif
	     
! Send w1 down unless I'm at the bottom

      
	if (nrky(nrank) .ge. 1 ) then
      w1y1(:,:)=w1(:,:,iy_first+2)
      w1y2(:,:)=w1(:,:,iy_first+3)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( w1y1, mxz, MPI_DOUBLE_PRECISION, nrank - nprxz, 1,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( w1y2, mxz, MPI_DOUBLE_PRECISION, nrank - nprxz, 1,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrky(nrank).lt. npry-1) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( w1y1, mxz, MPI_DOUBLE_PRECISION, nrank + nprxz, 1,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( w1y2, mxz, MPI_DOUBLE_PRECISION, nrank + nprxz, 1,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      w1(:,:,iy_last-1)=w1y1(:,:)
      w1(:,:,iy_last)=w1y2(:,:)
	endif

    if (nrky(nrank) .eq. 0 ) then
      w1y1(:,:)=w1(:,:,iy_first+2)
      w1y2(:,:)=w1(:,:,iy_first+3)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( w1y1, mxz, MPI_DOUBLE_PRECISION, nrank +(npry-1)*nprxz, 1,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( w1y2, mxz, MPI_DOUBLE_PRECISION, nrank +(npry-1)*nprxz, 1,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrky(nrank).eq. npry-1) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( w1y1, mxz, MPI_DOUBLE_PRECISION, nrank -(npry-1)*nprxz, 1,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( w1y2, mxz, MPI_DOUBLE_PRECISION, nrank -(npry-1)*nprxz, 1,  &
		      MPI_COMM_WORLD, status,ierror )
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

!****************************************************************************************
      subroutine mpi_transfersy1(wy1,wyt)
      USE DECLARE
      real*8, dimension(my) :: wy1
      real*8, dimension(myt) :: wyt
      include 'mpif.h'
      if(npry .gt. 1) then
      call MPI_Allgather(wy1(3:my-2),mym,MPI_DOUBLE_PRECISION,wyt,mym,MPI_DOUBLE_PRECISION,mycomm_y,ierror)
      else
      wyt(1:mym)=wy1(3:my-2)
      endif
      return
      end


!****************************************************************************************
      subroutine mpi_transfersm_2d(ws,mm)
      USE DECLARE
      real*8, dimension(mx,mz,my,mm) :: ws
      real*8, dimension(mz,my,mm) :: wsx1,wsx2
      real*8, dimension(mx,my,mm) :: wsz1,wsz2
      include 'mpif.h'

!       
! Send w8 up unless I'm at the top, then receive from below
     
	if (nrank.ge.nrkz(nrank)*nprx .and. nrank.lt.nrkz(nrank)*nprx+nprx-1) then
      wsx1(:,:,:)=ws(ix_last-2,:,:,:)
      wsx2(:,:,:)=ws(ix_last-3,:,:,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( wsx1, myz*mm, MPI_DOUBLE_PRECISION, nrank + 1, 0,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wsx2, myz*mm, MPI_DOUBLE_PRECISION, nrank + 1, 0,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrank.gt.nrkz(nrank)*nprx .and. nrank.le.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wsx1, myz*mm, MPI_DOUBLE_PRECISION, nrank - 1, 0,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wsx2, myz*mm, MPI_DOUBLE_PRECISION, nrank - 1, 0,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      ws(ix_first+1,:,:,:)=wsx1(:,:,:)
      ws(ix_first,:,:,:)=wsx2(:,:,:)
	endif
	
      
! Send w8 down unless I'm at the bottom

      
	if (nrank.gt.nrkz(nrank)*nprx .and. nrank.le.nrkz(nrank)*nprx+nprx-1) then
      wsx1(:,:,:)=ws(ix_first+2,:,:,:)
      wsx2(:,:,:)=ws(ix_first+3,:,:,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( wsx1, myz*mm, MPI_DOUBLE_PRECISION, nrank - 1, 1,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wsx2, myz*mm, MPI_DOUBLE_PRECISION, nrank - 1, 1,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrank.ge.nrkz(nrank)*nprx .and. nrank.lt.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wsx1, myz*mm, MPI_DOUBLE_PRECISION, nrank + 1, 1,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wsx2, myz*mm, MPI_DOUBLE_PRECISION, nrank + 1, 1,  &
		      MPI_COMM_WORLD, status,ierror )
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
	CALL MPI_Send( wsz1, myx*mm, MPI_DOUBLE_PRECISION, nrank + nprx, 0,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wsz2, myx*mm, MPI_DOUBLE_PRECISION, nrank + nprx, 0,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrank.ge.nprx ) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wsz1, myx*mm, MPI_DOUBLE_PRECISION, nrank - nprx, 0,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wsz2, myx*mm, MPI_DOUBLE_PRECISION, nrank - nprx, 0,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      ws(:,iz_first+1,:,:)=wsz1(:,:,:)
      ws(:,iz_first,:,:)=wsz2(:,:,:)
	endif
	     
! Send w8 down unless I'm at the bottom

      
	if (nrank.ge.nprx ) then
      wsz1(:,:,:)=ws(:,iz_first+2,:,:)
      wsz2(:,:,:)=ws(:,iz_first+3,:,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( wsz1, myx*mm, MPI_DOUBLE_PRECISION, nrank - nprx, 1,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wsz2, myx*mm, MPI_DOUBLE_PRECISION, nrank - nprx, 1,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrank.lt.nsize-nprx) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wsz1, myx*mm, MPI_DOUBLE_PRECISION, nrank + nprx, 1,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wsz2, myx*mm, MPI_DOUBLE_PRECISION, nrank + nprx, 1,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      ws(:,iz_last-1,:,:)=wsz1(:,:,:)
      ws(:,iz_last,:,:)=wsz2(:,:,:)
	endif
      return
      end 


!****************************************************************************************
      subroutine mpi_transfer8_2d(w8)
      USE DECLARE
      real*8, dimension(mx,mz,my,8) :: w8
      real*8, dimension(mz,my,8) :: wfx1,wfx2
      real*8, dimension(mx,my,8) :: wfz1,wfz2
      include 'mpif.h'

!       
! Send w8 up unless I'm at the top, then receive from below
     
	if (nrank.ge.nrkz(nrank)*nprx .and. nrank.lt.nrkz(nrank)*nprx+nprx-1) then
      wfx1(:,:,:)=w8(ix_last-2,:,:,:)
      wfx2(:,:,:)=w8(ix_last-3,:,:,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( wfx1, myz8, MPI_DOUBLE_PRECISION, nrank + 1, 0,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wfx2, myz8, MPI_DOUBLE_PRECISION, nrank + 1, 0,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrank.gt.nrkz(nrank)*nprx .and. nrank.le.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wfx1, myz8, MPI_DOUBLE_PRECISION, nrank - 1, 0,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wfx2, myz8, MPI_DOUBLE_PRECISION, nrank - 1, 0,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      w8(ix_first+1,:,:,:)=wfx1(:,:,:)
      w8(ix_first,:,:,:)=wfx2(:,:,:)
	endif
	
      
! Send w8 down unless I'm at the bottom

      
	if (nrank.gt.nrkz(nrank)*nprx .and. nrank.le.nrkz(nrank)*nprx+nprx-1) then
      wfx1(:,:,:)=w8(ix_first+2,:,:,:)
      wfx2(:,:,:)=w8(ix_first+3,:,:,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( wfx1, myz8, MPI_DOUBLE_PRECISION, nrank - 1, 1,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wfx2, myz8, MPI_DOUBLE_PRECISION, nrank - 1, 1,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrank.ge.nrkz(nrank)*nprx .and. nrank.lt.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wfx1, myz8, MPI_DOUBLE_PRECISION, nrank + 1, 1,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wfx2, myz8, MPI_DOUBLE_PRECISION, nrank + 1, 1,  &
		      MPI_COMM_WORLD, status,ierror )
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
	CALL MPI_Send( wfz1, myx8, MPI_DOUBLE_PRECISION, nrank + nprx, 0,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wfz2, myx8, MPI_DOUBLE_PRECISION, nrank + nprx, 0,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrank.ge.nprx ) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wfz1, myx8, MPI_DOUBLE_PRECISION, nrank - nprx, 0,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wfz2, myx8, MPI_DOUBLE_PRECISION, nrank - nprx, 0,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      w8(:,iz_first+1,:,:)=wfz1(:,:,:)
      w8(:,iz_first,:,:)=wfz2(:,:,:)
	endif
	     
! Send w8 down unless I'm at the bottom

      
	if (nrank.ge.nprx ) then
      wfz1(:,:,:)=w8(:,iz_first+2,:,:)
      wfz2(:,:,:)=w8(:,iz_first+3,:,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( wfz1, myx8, MPI_DOUBLE_PRECISION, nrank - nprx, 1,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wfz2, myx8, MPI_DOUBLE_PRECISION, nrank - nprx, 1,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrank.lt.nsize-nprx) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wfz1, myx8, MPI_DOUBLE_PRECISION, nrank + nprx, 1,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wfz2, myx8, MPI_DOUBLE_PRECISION, nrank + nprx, 1,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      w8(:,iz_last-1,:,:)=wfz1(:,:,:)
      w8(:,iz_last,:,:)=wfz2(:,:,:)
	endif
      return
      end 
!****************************************************************************************
      subroutine mpi_transfer3_2d(w3)
      USE DECLARE
      real*8, dimension(mx,mz,my,3) :: w3
      real*8, dimension(mz,my,3) :: wcx1,wcx2
      real*8, dimension(mx,my,3) :: wcz1,wcz2
      include 'mpif.h'

       
! Send w3 up unless I'm at the top, then receive from below
     
	if (nrank.ge.nrkz(nrank)*nprx .and. nrank.lt.nrkz(nrank)*nprx+nprx-1) then
      wcx1(:,:,:)=w3(ix_last-2,:,:,:)
      wcx2(:,:,:)=w3(ix_last-3,:,:,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( wcx1, myz3, MPI_DOUBLE_PRECISION, nrank + 1, 0,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wcx2, myz3, MPI_DOUBLE_PRECISION, nrank + 1, 0,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrank.gt.nrkz(nrank)*nprx .and. nrank.le.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wcx1, myz3, MPI_DOUBLE_PRECISION, nrank - 1, 0,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wcx2, myz3, MPI_DOUBLE_PRECISION, nrank - 1, 0,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      w3(ix_first+1,:,:,:)=wcx1(:,:,:)
      w3(ix_first,:,:,:)=wcx2(:,:,:)
	endif
	
      
! Send w3 down unless I'm at the bottom

      
	if (nrank.gt.nrkz(nrank)*nprx .and. nrank.le.nrkz(nrank)*nprx+nprx-1) then
      wcx1(:,:,:)=w3(ix_first+2,:,:,:)
      wcx2(:,:,:)=w3(ix_first+3,:,:,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( wcx1, myz3, MPI_DOUBLE_PRECISION, nrank - 1, 1,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wcx2, myz3, MPI_DOUBLE_PRECISION, nrank - 1, 1,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrank.ge.nrkz(nrank)*nprx .and. nrank.lt.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wcx1, myz3, MPI_DOUBLE_PRECISION, nrank + 1, 1,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wcx2, myz3, MPI_DOUBLE_PRECISION, nrank + 1, 1,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      w3(ix_last-1,:,:,:)=wcx1(:,:,:)
      w3(ix_last,:,:,:)=wcx2(:,:,:)
	endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	if (nrank.lt.nsize-nprx) then
      wcz1(:,:,:)=w3(:,iz_last-2,:,:)
      wcz2(:,:,:)=w3(:,iz_last-3,:,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( wcz1, myx3, MPI_DOUBLE_PRECISION, nrank + nprx, 0,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wcz2, myx3, MPI_DOUBLE_PRECISION, nrank + nprx, 0,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrank.ge.nprx ) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wcz1, myx3, MPI_DOUBLE_PRECISION, nrank - nprx, 0,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wcz2, myx3, MPI_DOUBLE_PRECISION, nrank - nprx, 0,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      w3(:,iz_first+1,:,:)=wcz1(:,:,:)
      w3(:,iz_first,:,:)=wcz2(:,:,:)
	endif
	     
! Send w3 down unless I'm at the bottom

      
	if (nrank.ge.nprx ) then
      wcz1(:,:,:)=w3(:,iz_first+2,:,:)
      wcz2(:,:,:)=w3(:,iz_first+3,:,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( wcz1, myx3, MPI_DOUBLE_PRECISION, nrank - nprx, 1,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( wcz2, myx3, MPI_DOUBLE_PRECISION, nrank - nprx, 1,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrank.lt.nsize-nprx) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( wcz1, myx3, MPI_DOUBLE_PRECISION, nrank + nprx, 1,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( wcz2, myx3, MPI_DOUBLE_PRECISION, nrank + nprx, 1,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      w3(:,iz_last-1,:,:)=wcz1(:,:,:)
      w3(:,iz_last,:,:)=wcz2(:,:,:)
	endif
      return
      end

!****************************************************************************************
      subroutine mpi_transfer1_2d(w1)
      USE DECLARE
      real*8, dimension(mx,mz,my) :: w1
      real*8, dimension(mz,my) :: w1x1,w1x2
      real*8, dimension(mx,my) :: w1z1,w1z2
      include 'mpif.h'

!       
! Send w1 up unless I'm at the top, then receive from below
     
	if (nrank.ge.nrkz(nrank)*nprx .and. nrank.lt.nrkz(nrank)*nprx+nprx-1) then
      w1x1(:,:)=w1(ix_last-2,:,:)
      w1x2(:,:)=w1(ix_last-3,:,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( w1x1, myz, MPI_DOUBLE_PRECISION, nrank + 1, 0,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( w1x2, myz, MPI_DOUBLE_PRECISION, nrank + 1, 0,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrank.gt.nrkz(nrank)*nprx .and. nrank.le.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( w1x1, myz, MPI_DOUBLE_PRECISION, nrank - 1, 0,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( w1x2, myz, MPI_DOUBLE_PRECISION, nrank - 1, 0,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      w1(ix_first+1,:,:)=w1x1(:,:)
      w1(ix_first,:,:)=w1x2(:,:)
	endif
	
      
! Send w1 down unless I'm at the bottom

      
	if (nrank.gt.nrkz(nrank)*nprx .and. nrank.le.nrkz(nrank)*nprx+nprx-1) then
      w1x1(:,:)=w1(ix_first+2,:,:)
      w1x2(:,:)=w1(ix_first+3,:,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( w1x1, myz, MPI_DOUBLE_PRECISION, nrank - 1, 1,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( w1x2, myz, MPI_DOUBLE_PRECISION, nrank - 1, 1,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrank.ge.nrkz(nrank)*nprx .and. nrank.lt.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( w1x1, myz, MPI_DOUBLE_PRECISION, nrank + 1, 1,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( w1x2, myz, MPI_DOUBLE_PRECISION, nrank + 1, 1,  &
		      MPI_COMM_WORLD, status,ierror )
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
	CALL MPI_Send( w1z1, myx, MPI_DOUBLE_PRECISION, nrank + nprx, 0,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( w1z2, myx, MPI_DOUBLE_PRECISION, nrank + nprx, 0,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
	if (nrank.ge.nprx ) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( w1z1, myx, MPI_DOUBLE_PRECISION, nrank - nprx, 0,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( w1z2, myx, MPI_DOUBLE_PRECISION, nrank - nprx, 0,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      w1(:,iz_first+1,:)=w1z1(:,:)
      w1(:,iz_first,:)=w1z2(:,:)
	endif
	     
! Send w1 down unless I'm at the bottom

      
	if (nrank.ge.nprx ) then
      w1z1(:,:)=w1(:,iz_first+2,:)
      w1z2(:,:)=w1(:,iz_first+3,:)
!mpi   ----------------------------------------------------------------
	CALL MPI_Send( w1z1, myx, MPI_DOUBLE_PRECISION, nrank - nprx, 1,  &
		      MPI_COMM_WORLD,ierror )
	CALL MPI_Send( w1z2, myx, MPI_DOUBLE_PRECISION, nrank - nprx, 1,  &
		      MPI_COMM_WORLD,ierror )
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrank.lt.nsize-nprx) then
!mpi   ----------------------------------------------------------------
	CALL MPI_Recv( w1z1, myx, MPI_DOUBLE_PRECISION, nrank + nprx, 1,  &
		      MPI_COMM_WORLD, status,ierror )
	CALL MPI_Recv( w1z2, myx, MPI_DOUBLE_PRECISION, nrank + nprx, 1,  &
		      MPI_COMM_WORLD, status,ierror )
!mpi   ----------------------------------------------------------------
      w1(:,iz_last-1,:)=w1z1(:,:)
      w1(:,iz_last,:)=w1z2(:,:)
	endif
      return
      end

!ws******************************************************
!c     3h-quad - - quadratic interpolation routine
!c     3i-quad2 - - quadratic interpolation routine
!c     3j-cubic - - cubic interploation routine
!c     3k-newline - - forms arc length along a curve:linear version
!c     3l-newlin2 - - interpolates along a curve      
      subroutine quad(x1,x2,x3,y1,y2,y3,y,ans,ierr,deriv)
!cray  lcm(quad)
      d1 = (x1-x2)*(x1-x3)
      d2 = (x2-x1)*(x2-x3)
      d3 = (x3-x1)*(x3-x2)
      if(d1.eq.0 .or. d2.eq.0 .or. d3.eq.0) go to 350
      aco = y1/d1 + y2/d2 + y3/d3
      bco = -(y1*(x2+x3)/d1 + y2*(x1+x3)/d2 + y3*(x1+x2)/d3)
      cco = y1*x2*x3/d1 + y2*x1*x3/d2 + y3*x1*x2/d3 - y
      if(aco .eq. 0) go to 150
      disc = (bco**2 - 4.*aco*cco)
      if(disc.lt.0) go to 100
      disc = sqrt(disc)
      root1 = (-bco + disc)/(2.*aco)
      root2 = (-bco - disc)/(2.*aco)
      if(root1.lt.x1 .or. root1.gt.x3) go to 250
      if(root2.lt.x1 .or. root2.gt.x3) go to 260
  240 if(abs(root1-x2) .lt. abs(root2-x2)) go to 265
      go to 255
  250 continue
      if(root2.lt.x1 .or. root2.gt.x3) go to 240
  255 ans = root2
      go to 270
  260 continue
      if(root1.lt.x1 .or. root1.gt.x3) go to 240
  265 ans = root1
  270 continue
      denom = 2.*aco*ans + bco
      if(denom.eq.0) go to 100
      deriv = 1./denom
      ierr = 0
      return
  100 continue
      ierr = 1
      return
  150 continue
      ans = -cco/bco
      deriv = 1./bco
      return
  350 continue
      ans = 0.
      deriv = 0.
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     * * * number 3i * * *                                c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine interp1d2l(x1,x2,x3,y1,y2,y3,y,ans)
!cray  lcm(quad2)
      d1 = (y1-y2)*(y1-y3)
      d2 = (y2-y3)*(y2-y1)
      d3 = (y3-y1)*(y3-y2)
      ans = x1*(y-y2)*(y-y3)/d1 &
          + x2*(y-y3)*(y-y1)/d2 &
          + x3*(y-y1)*(y-y2)/d3

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     * * * number 3j * * *                                c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine interp1d3l(x1,x2,x3,x4,y1,y2,y3,y4,y,ans)
      real*8 x1,x2,x3,x4,y1,y2,y3,y4,y,ans
!cray  lcm(cubic)
!c
!c      this routine returns interpolated value 
!c
      d1 = (y1-y2)*(y1-y3)*(y1-y4)
      d2 = (y2-y1)*(y2-y3)*(y2-y4)
      d3 = (y3-y1)*(y3-y2)*(y3-y4)
      d4 = (y4-y1)*(y4-y2)*(y4-y3)
      ans = x1*(y-y2)*(y-y3)*(y-y4)/d1 &
          + x2*(y-y1)*(y-y3)*(y-y4)/d2 &
          + x3*(y-y1)*(y-y2)*(y-y4)/d3 &
          + x4*(y-y1)*(y-y2)*(y-y3)/d4 

      return
      end


      subroutine interp2d3l(x1,x2,z1,z2,w11,w12,w21,w22,xw,zw,ans)
      real*8 x1,x2,z1,z2,w11,w12,w21,w22,xw,zw,ans,beta1,beta2
!cray  !c
!c      this routine returns interpolated value
!
        beta1=(xw-x1)/(x2-x1)
        beta2=(zw-z1)/(z2-z1)
        ans=w11*(1-beta1)*(1-beta2) &
           +w21*beta1*(1-beta2) &
           +w12*(1-beta1)*beta2 &
           +w22*beta1*beta2

      return
      end

      subroutine interp2d4p(x11,x12,x21,x22,z11,z12,z21,z22,w11,w12,w21,w22,xw,zw,ans)
      real*8 x11,x12,x21,x22,z11,z12,z21,z22,w11,w12,w21,w22,xw,zw,ans,r11,r12,r21,r22,c11,c12,c21,c22,ccc
!cray  !c
!c      this routine returns interpolated value
!
        r11=sqrt((x11-xw)**2+(z11-zw)**2)
        r12=sqrt((x12-xw)**2+(z12-zw)**2)
        r21=sqrt((x21-xw)**2+(z21-zw)**2)
        r22=sqrt((x22-xw)**2+(z22-zw)**2)
        if(r11 .eq. 0) then
        ans=w11
        return
        else if(r12 .eq. 0) then 
        ans=w12
        return
        else if(r21 .eq. 0) then 
        ans=w21
        return
        else if(r22 .eq. 0) then 
        ans=w22
        return
        endif
        
        c11=1./r11
        c12=1./r12
        c21=1./r21
        c22=1./r22
        ccc=c11+c12+c21+c22

        ans=(w11*c11+w12*c12+w21*c21+w22*c22)/ccc

      return
      end
!******************************************************************
      subroutine interp_xz2ps_unif(bxz,xi,zi,nx,nz,bps,xs,zs,ns)
      implicit real*8 (b-h,o-z)
      dimension xi(nx),zi(nz),bxz(nx,nz),xs(ns),zs(ns),bps(ns)
      lx=1
      lz=1
    
      do 30 is=1,ns
      lx=int((xs(is)-xi(1))/(xi(nx)-xi(1))*(nx-1))+1
      lz=int((zs(is)-zi(1))/(zi(nz)-zi(1))*(nz-1))+1

      beta1=(xs(is)-xi(lx))/(xi(lx+1)-xi(lx))
      beta2=(zs(is)-zi(lz))/(zi(lz+1)-zi(lz))
    
      bps(is) =bxz(lx,lz)*(1-beta1)*(1-beta2) &
              +bxz(lx+1,lz)*beta1*(1-beta2) &
              +bxz(lx,lz+1)*(1-beta1)*beta2 &
              +bxz(lx+1,lz+1)*beta1*beta2
   30 continue

   
      return
      end
!******************************************************************
      subroutine interp_xz2ps(bxz,xi,zi,nx,nz,bps,xs,zs,ns)
      implicit real*8 (b-h,o-z)
      dimension xi(nx),zi(nz),bxz(nx,nz),xs(ns),zs(ns),bps(ns)
      lx=1
      lz=1
!mpi   ----------------------------------------------------------------
!      CALL MPI_ALLGather(ff,m3d,MPI_DOUBLE_COMPLEX,xfc,m3d, &
!            MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,IERROR)
!mpi   ----------------------------------------------------------------
 !     
      do 30 is=1,ns
      do 10 ix=1,nx-1
      if((xi(ix)-1.d-7).le.xs(is) .and. (xi(ix+1)+1.d-7).ge.xs(is)) then 
      beta1=(xs(is)-xi(ix))/(xi(ix+1)-xi(ix))
       lx=ix
!      stop'1'
      goto 15
      else
      endif
   10 continue
   15 do 20 iz=1,nz-1
      if((zi(iz)-1.d-7).le.zs(is).and.(zi(iz+1)+1.d-7).ge.zs(is)) then
      beta2=(zs(is)-zi(iz))/(zi(iz+1)-zi(iz))
      lz=iz
!	stop'2'
      goto 25
      else
      endif
   20 continue
    
   25 bps(is) =bxz(lx,lz)*(1-beta1)*(1-beta2) &
              +bxz(lx+1,lz)*beta1*(1-beta2) &
              +bxz(lx,lz+1)*(1-beta1)*beta2 &
              +bxz(lx+1,lz+1)*beta1*beta2
   30 continue

   
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     * * * number 3j * * *                                c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cubic(x1,x2,x3,x4,y1,y2,y3,y4,y,ans,dans,is,ierr)
!cray  lcm(cubic)
!c
!c      this routine returns interpolated value and first derivative
!c
      d1 = (y1-y2)*(y1-y3)*(y1-y4)
      d2 = (y2-y1)*(y2-y3)*(y2-y4)
      d3 = (y3-y1)*(y3-y2)*(y3-y4)
      d4 = (y4-y1)*(y4-y2)*(y4-y3)
      ans = x1*(y-y2)*(y-y3)*(y-y4)/d1 &
          + x2*(y-y1)*(y-y3)*(y-y4)/d2 &
          + x3*(y-y1)*(y-y2)*(y-y4)/d3 &
          + x4*(y-y1)*(y-y2)*(y-y3)/d4 

      if ( is .eq. 0 )  return
      if(is.eq.2) go to 5

      dans = x1/d1*((y-y3)*(y-y4)+(y-y2)*(y-y4)+(y-y2)*(y-y3)) &
           + x2/d2*((y-y3)*(y-y4)+(y-y1)*(y-y4)+(y-y1)*(y-y3)) &
           + x3/d3*((y-y2)*(y-y4)+(y-y1)*(y-y4)+(y-y1)*(y-y2)) &
           + x4/d4*((y-y2)*(y-y3)+(y-y1)*(y-y3)+(y-y1)*(y-y2))
      return
    5 continue
      dans = x1/d1*((y-y3)+(y-y4)+(y-y2)+(y-y4)+(y-y2)+(y-y3)) &
           + x2/d2*((y-y3)+(y-y4)+(y-y1)+(y-y4)+(y-y1)+(y-y3)) &
           + x3/d3*((y-y2)+(y-y4)+(y-y1)+(y-y4)+(y-y1)+(y-y2)) &
           + x4/d4*((y-y2)+(y-y3)+(y-y1)+(y-y3)+(y-y1)+(y-y2))
      return
   20 continue
      anum = (y2-y3)*(x1-x2)+(y1-y2)*(x3-x2)
      adem = (y1-y3)*(y1-y2)**2
      a = anum/adem
      b = (x3-x1)/(y3-y1) - 2.*y2*a
      dans = 2.*a*y + b
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     * * * number 3k * * *                                c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine newline(xp,zp,sp,np)
!cray  lcm(newline)
!c
!c     linear version
!c
      dimension xp(np),zp(np),sp(np)
      sp(1) = 0.
      do 10 j=2,np
      sp(j) = sp(j-1) + sqrt((xp(j)-xp(j-1))**2 + (zp(j)-zp(j-1))**2)
   10 continue
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     * * * number 3l * * *                                c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine newlin2(xp,zp,xpn,zpn,sp,spn,np,np2)
!cray  lcm(newlin2)
!c
!c     linear version
!c
      dimension xp(np2),zp(np2),xpn(np),zpn(np),sp(np2),spn(np)
      xpn(1) = xp(1)
      zpn(1) = zp(1)
      jj = 2
      do 10 j=1,np
    5 continue
      if(spn(j).lt.sp(jj)) go to 20
      if(jj.eq.np2) go to 20
      jj = jj + 1
      if(jj.eq.np2) go to 20
      go to 5
   20 continue
!c
!c     spn(j) lies between sp(jj-1) and sp(jj)
!c
      fac = (spn(j) - sp(jj-1))/(sp(jj) - sp(jj-1))
      xpn(j) = xp(jj-1) + fac*(xp(jj)-xp(jj-1))
      zpn(j) = zp(jj-1) + fac*(zp(jj)-zp(jj-1))
   10 continue
      return
      end

!****************************************************
      subroutine extap(x1,x2,x3,x4)
      real*8 x1,x2,x3,x4
      x4=3.*x3-3.*x2+x1
      ddx1=x3-x2
      ddx2=x2-x1
      ddx=x4-x3
      pm=ddx*ddx1
      if(pm.gt.0.) return
      if(ddx2.eq.0) x4=2.*x3-x2
      if(ddx2.eq.0) return
      ddx=(ddx1*ddx1)/ddx2
      x4=x3+ddx
      return
      end





!*************************************************************
      subroutine extrap_xz2rt(bxz,x0,xi,zi,nx,nz,brt,rh,th,nt)
      implicit real*8 (b-h,o-z)
      dimension xi(nx),zi(nz),bxz(nx,nz)
      dimension th(nt),brt(nt)
      real*8 rh,xl,zl
      lx=1
      lz=1
!mpi   ----------------------------------------------------------------
!      CALL MPI_ALLGather(ff,m3d,MPI_DOUBLE_COMPLEX,xfc,m3d, &
!            MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,IERROR)
!mpi   ----------------------------------------------------------------
 !     do 30 jr=1,nr
      do 30 jt=1,nt
      xl=x0+rh*dcos(th(jt))
      zl=rh*dsin(th(jt))
      do 10 ix=1,nx-1
      if((xi(ix)-1.d-7).le.xl .and. (xi(ix+1)+1.d-7).ge.xl) then 
      beta1=(xl-xi(ix))/(xi(ix+1)-xi(ix))
       lx=ix
!      stop'1'
      goto 15
      else
      endif
   10 continue
   15 do 20 iz=1,nz-1
      if((zi(iz)-1.d-7).le.zl.and.(zi(iz+1)+1.d-7).ge.zl) then
      beta2=(zl-zi(iz))/(zi(iz+1)-zi(iz))
      lz=iz
!	stop'2'
      goto 25
      else
      endif
   20 continue
    
   25 brt(jt)=bxz(lx,lz)*(1-beta1)*(1-beta2) &
              +bxz(lx+1,lz)*beta1*(1-beta2) &
              +bxz(lx,lz+1)*(1-beta1)*beta2 &
              +bxz(lx+1,lz+1)*beta1*beta2
   30 continue


!      brtt(jt)=

      return
      end

!*************************************************************
      subroutine extrap2_xz2rt(bxz,x0,xi,zi,nx,nz,brt,brt_dr,brt_dt,rh,th,nt)
      implicit real*8 (b-h,o-z)
      dimension xi(nx),zi(nz),bxz(nx,nz)
      dimension th(nt),brt(nt),brt_dr(nt),brt_dt(nt)
      real*8 rh,xl,zl,brt_dx,brt_dz

      lx=1
      lz=1
!mpi   ----------------------------------------------------------------
!      CALL MPI_ALLGather(ff,m3d,MPI_DOUBLE_COMPLEX,xfc,m3d, &
!            MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,IERROR)
!mpi   ----------------------------------------------------------------
 !     do 30 jr=1,nr
      do 30 jt=1,nt
      xl=x0+rh*dcos(th(jt))
      zl=rh*dsin(th(jt))
      do 10 ix=1,nx-1
      if((xi(ix)-1.d-7).le.xl .and. (xi(ix+1)+1.d-7).ge.xl) then 
      beta1=(xl-xi(ix))/(xi(ix+1)-xi(ix))
       lx=ix
!      stop'1'
      goto 15
      else
      endif
   10 continue
   15 do 20 iz=1,nz-1
      if((zi(iz)-1.d-7).le.zl.and.(zi(iz+1)+1.d-7).ge.zl) then
      beta2=(zl-zi(iz))/(zi(iz+1)-zi(iz))
      lz=iz
!	stop'2'
      goto 25
      else
      endif
   20 continue
    
   25 brt(jt)=bxz(lx,lz)*(1-beta1)*(1-beta2) &
              +bxz(lx+1,lz)*beta1*(1-beta2) &
              +bxz(lx,lz+1)*(1-beta1)*beta2 &
              +bxz(lx+1,lz+1)*beta1*beta2
      brt_dx=((bxz(lx+1,lz)-bxz(lx,lz))*(1-beta2) &
             +(bxz(lx+1,lz+1)-bxz(lx,lz+1))*beta2)/(xi(lx+1)-xi(lx))
      brt_dz=((bxz(lx,lz+1)-bxz(lx,lz))*(1-beta1) &
             +(bxz(lx+1,lz+1)-bxz(lx+1,lz))*beta1)/(zi(lz+1)-zi(lz))
      brt_dr(jt)=(xl*brt_dx+zl*brt_dz)/rh
      brt_dt(jt)=xl*brt_dz-zl*brt_dx
   30 continue


!      brtt(jt)=

      return
      end

      subroutine extrap_rt2xz(brt1,rh1,brt2,rh2,th,nt,bxz,ri,thi,nb)
      implicit real*8 (b-h,o-z)
      dimension ri(nb),thi(nb),bxz(nb)
      dimension th(nt),brt1(nt),brt2(nt)
!      lx=1
!      lz=1
      do 30 i=1,nb
      beta2=(ri(i)-rh1)/(rh2-rh1)
      do 10 it=1,nt-1
      if((th(it)-1.d-7).le.thi(i).and.(th(it+1)+1.d-7).ge.thi(i)) then 
      beta1=(thi(i)-th(it))/(th(it+1)-th(it))
      lt=it
      goto 25
      else
      endif
   10 continue
   
   25 bxz(i)=brt1(lt)*(1-beta1)*(1-beta2) &
             +brt1(lt+1)*beta1*(1-beta2) &
             +brt2(lt)*(1-beta1)*beta2 &
             +brt2(lt+1)*beta1*beta2 
   30 continue

      return
      end

      subroutine extrap1d0_r(brt2,rh2,brt3,rh3,nt)
      implicit real*8 (b-h,o-z)
      dimension brt2(nt),brt3(nt)
      do jt=1,nt
      brt3(jt)=brt2(jt)
      enddo
      return
      end

      subroutine extrap1d1_r(brt1,rh1,brt2,rh2,brt3,rh3,nt)
      implicit real*8 (b-h,o-z)
      dimension brt1(nt),brt2(nt),brt3(nt)
      beta1=(rh3-rh1)/(rh2-rh1)
      do 20 jt=1,nt
      brt3(jt)=brt1(jt)*(1-beta1)+brt2(jt)*beta1
   20 continue
      return
      end

!!ws:dmapb2+++++++++++++++++++++++++++++++++++++++++++++
!      subroutine intp1(px,ny,pxa,ns,nhf)
!      integer ny,ns,nhf
!      dimension px(ny),pxa(ns),ytemp(ny)
!      do 2 j=1,ny
!      ytemp(j)=(j-1.)/(ny-1)
!      if(nhf.eq.1) ytemp(j)=(j-0.5)/(ny-1)
!2     continue
!
!      j=4
!      do 3 jj=1,nosurf
!      pval=ps(jj)
!4     continue
!      if(j.eq.npsim) go to 5
!      if(pval.lt.ytemp(j)) go to 5
!      j=j+1
!      if(j.ge.npsim) go to 5
!      go to 4
!5     continue
!      call cubic(px(j-2),px(j-1),px(j),px(j+1) &
!      ytemp(j-2),ytemp(j-1),ytemp(j),ytemp(j+1) &
!      pval,pxa(jj),dum,0,ierr)
!3     continue
!      return
!      end
!
!
!      subroutine intp2(p2x,p2xa,nhf,nsignx,inv)
!      include 'clichpar.h'
!      dimension p2x(nq,ny),p2xb(nq,ns),p2xa(nt,ns)
!      common/var1/ npsi,nthe2,dth,dr,nptsi,jndexi
!     1,nosurf,npsim,npsip,nthe,nthe1,nthe3,nthe4,mth,pi
!      common/var3/fg(nsgrd),psig(nsgrd),f(ny),psival(ny),fp(ny)
!      common/work/ytemp(nsgrd),xtemp(nsgrd),ztemp(nsgrd)
!      do 2 j=1,npsi
!      ytemp(j)=(j-1.)/npsim
!      if(nhf.eq.1)
!     &ytemp(j)=(j-0.5)/npsim
!2     continue
!      imin=3
!      imax=nthe2
!      if(nsignx.eq.-1) imin=4
!      if(nsignx.eq.-1) imax=nthe2-1
!      do 1 i=imin,imax
!      do 11 kk=2,npsi
!      xtemp(kk)=p2x(i,kk)
!      if(inv.eq.1) xtemp(kk)=1./xtemp(kk)
!   11 continue
!      j=4
!      do 3 jj=1,nosurf
!      pval=psig(jj)
!4     continue
!      if(j.eq.npsim) go to 5
!      if(pval.lt.ytemp(j)) go to 5
!      j=j+1
!      if(j.ge.npsim) go to 5
!      go to 4
!5     continue
!6     call cubic
!     1(xtemp(j-2),xtemp(j-1),xtemp(j),xtemp(j+1)
!     2,ytemp(j-2),ytemp(j-1),ytemp(j),ytemp(j+1)
!     3,pval,xnew,dum,0,ierr)
!      if(inv.eq.1) xnew=1./xnew
!      p2xb(i,jj)=xnew
!3     continue
!1     continue
!      do 8 j=1,nosurf
!      p2xb(2,j)=p2xb(4,j)*nsignx
!      p2xb(1,j)=p2xb(5,j)*nsignx
!      p2xb(nthe3,j)=p2xb(nthe2-1,j)*nsignx
!      p2xb(nthe4,j)=p2xb(nthe2-2,j)*nsignx
!      if(nsignx.eq.-1) p2xb(3,j)=0.0
!      if(nsignx.eq.-1) p2xb(nthe2,j)=0.0
!8     continue
!7     continue
!      call intt(p2xb,p2xa,nsignx)
!      return
!      end
!      subroutine int2spl(p2x,p2xa,nhf,nsignx,inv)
!      include 'clichpar.h'
!      dimension p2x(nq,ny),p2xb(nq,ns),p2xa(nt,ns)
!      common/var1/ npsi,nthe2,dth,dr,nptsi,jndexi
!     1,nosurf,npsim,npsip,nthe,nthe1,nthe3,nthe4,mth,pi
!      common/var3/fg(nsgrd),psig(nsgrd),f(ny),psival(ny),fp(ny)
!      common/work/ytemp(nsgrd),xtemp(nsgrd),ztemp(nsgrd)
!      imin=3
!      imax=nthe2
!      if(nsignx.eq.-1) imin=4
!      if(nsignx.eq.-1) imax=nthe2-1
!      npsip2=npsi+2
!      do 1 i=imin,imax
!      do 11 kk=1,npsi
!      xtemp(kk)=p2x(i,kk)
!      if(inv.eq.1) xtemp(kk)=1./xtemp(kk)
!   11 continue
!c
!      call intspl(nosurf,npsip2,xtemp,ytemp,psig,ztemp)
!c
!      if(inv.ne.1) go to 2
!      do 3 jj=1,nosurf
!      p2xb(i,jj)=1./ytemp(jj)
!3     continue
!      go to 1
!    2 continue
!      do 4 jj=1,nosurf
!      p2xb(i,jj)=ytemp(jj)
!4     continue
!1     continue
!c
!      do 8 j=1,nosurf
!      p2xb(2,j)=p2xb(4,j)*nsignx
!      p2xb(1,j)=p2xb(5,j)*nsignx
!      p2xb(nthe3,j)=p2xb(nthe2-1,j)*nsignx
!      p2xb(nthe4,j)=p2xb(nthe2-2,j)*nsignx
!      if(nsignx.eq.-1) p2xb(3,j)=0.0
!      if(nsignx.eq.-1) p2xb(nthe2,j)=0.0
!8     continue
!c
!      call intt(p2xb,p2xa,nsignx)
!      return
!      end
!      subroutine intt(p2xb,p2xa,nsignx)
!      include 'clichpar.h'
!      dimension p2xb(nq,ns) ,p2xa(nt,ns)
!      common/work/ytemp(nsgrd),xtemp(nsgrd),ztemp(nsgrd)
!      common/var1/ npsi,nthe2,dth,dr,nptsi,jndexi
!     1,nosurf,npsim,npsip,nthe,nthe1,nthe3,nthe4,mth,pi
!      common/var4/tb(nq,ns)
!      mthx=mth/2+1
!      rmth=mth
!      mth1=mth+1
!      mth2=mth+2
!      do 1 j=1,nosurf
!      do 2 i=2,nthe3
!      ytemp(i)=tb(i,j)
!2     continue
!      i=4
!      do 3 ii=1,mthx
!      tval=(ii-1.)*2*pi/rmth
!4     continue
!      if(i.eq.nthe2) go to 5
!      if(tval.lt.ytemp(i)) go to 5
!      i=i+1
!      if(i.eq.nthe2) go to 5
!      go to 4
!5     continue
!      call cubic
!     1(p2xb(i-2,j),p2xb(i-1,j),p2xb(i,j),p2xb(i+1,j)
!     2,ytemp(i-2),ytemp(i-1),ytemp(i),ytemp(i+1)
!     3,tval,p2xa(ii,j),dub,0,ierr)
!3     continue
!1     continue
!      mthxm=mthx-1
!      do 7 j=1,nosurf
!      do 6 i=2,mthxm
!      ii=mth-i+2
!      p2xa(ii,j)=p2xa(i,j)*nsignx
!6     continue
!      if(nsignx.eq.-1) p2xa(1,j)=0.0
!      if(nsignx.eq.-1) p2xa(mthx,j)=0.0
!      p2xa(mth1,j)=p2xa(1,j)
!      p2xa(mth2,j)=p2xa(2,j)
!7     continue
!      return
!      end
!      subroutine extap(x1,x2,x3,x4)
!      include 'clichpar.h'
!      x4=3.*x3-3.*x2+x1
!      ddx1=x3-x2
!      ddx2=x2-x1
!      ddx=x4-x3
!      pm=ddx*ddx1
!      if(pm.gt.0.) return
!      if(ddx2.eq.0.) x4=2.*x3-x2
!      if(ddx2.eq.0.) return
!      ddx=(ddx1*ddx1)/ddx2
!      x4=x3+ddx
!      return
!      end
!      subroutine ex(d,n)
!      include 'clichpar.h'
!      common/var1/ npsi,nthe2,dth,dr,nptsi,jndexi
!     1,nosurf,npsim,npsip,nthe,nthe1,nthe3,nthe4,mth,pi
!      dimension d(nq,ny)
!      do 1 i=2,nthe3
!      if(n.eq.1)
!     1call extap(d(i,4),d(i,3),d(i,2),d(i,1))
!      if(n.eq.2)
!     1call extap(d(i,5),d(i,4),d(i,3),d(i,2))
!      if(n.eq.npsi)
!     1call extap(d(i,npsi-3),d(i,npsi-2),d(i,npsim),d(i,npsi))
!1     continue
!      return
!      end
!      subroutine sym(d,nsignx)
!      include 'clichpar.h'
!      common/var1/ npsi,nthe2,dth,dr,nptsi,jndexi
!     1,nosurf,npsim,npsip,nthe,nthe1,nthe3,nthe4,mth,pi
!      dimension d(nq,ny)
!      do 1 j=1,npsi
!      d(2,j)=d(4,j)*nsignx
!      d(nthe3,j)=d(nthe2-1,j)*nsignx
!      if(nsignx.ne.-1) go to 1
!      d(3,j)=0.0
!      d(nthe2,j)=0.0
!1     continue
!      return
!      end
!
!!ws:dmapb2----------------------------------------


!ws**********************************************************************

      subroutine smth_x1st(ms,me,jps,kk)
      USE DECLARE
      real*8, dimension(n2th+5) :: wst
      include 'mpif.h'
!
! second-order diffusion
      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)

!      difxz(fm1m1,fp1m1,fm1p1,fp1p1,xm1,x0,xp1,zm1,z0,zp1)= &
!       2.*( fp1p1-fm1p1-fp1m1-fm1m1)/((xp1-xm1)*(zp1-zm1)) &
!       *((xp1-x0)*(x0-xm1)*(zp1-z0)*(z0-zm1))**0.5
!
!      difxy(fm1,f0,fp1,xm1,x0,xp1,dyy)= &
!       (2.*((xm1-x0)/(xp1-x0)*(fp1-f0)-(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)/x0-f0/x0**2) &
!       *((xp1-x0)*(x0-xm1))**0.5*dyy
!
!      difzy(fm1,f0,fp1,zm1,z0,zp1,dyy,x0)= &
!       (2.*((zm1-z0)/(zp1-z0)*(fp1-f0)-(zp1-z0)/(zm1-z0)*(fm1-f0))/(zm1-zp1)/x0) &
!       *((zp1-z0)*(z0-zm1))**0.5*dyy


       integer status(mpi_status_size)
      do 10 k=1,kk
!      call convt
      do 11 m=ms,me
!      if(m.eq.1.or.m.eq.2.or.m.eq.5.or.m.eq.8) then
      do 12 jy=1,my

      do 21 jt=2,n2th+4
      wst(jt)=difc(x1st(jt-1,jps,jy,m),x1st(jt,jps,jy,m),x1st(jt+1,jps,jy,m),thst(jt-1),thst(jt),thst(jt+1))
   21 continue
   

      
      do 22 jt=2,n2th+4
      x1st(jt,jps,jy,m)=x1st(jt,jps,jy,m)+1./4.*wst(jt)
   22 continue

   12 continue 
!      endif
   11 continue
   10 continue
      return
      end

!ws*******************************************
      subroutine smth_ps1(bst,kk)
      USE DECLARE
      real*8, dimension(n2th+5) :: bst,wst
      include 'mpif.h'
!
! second-order diffusion
      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)

       integer status(mpi_status_size)
      do 10 k=1,kk
      do 21 jt=2,n2th+4
      wst(jt)=difc(bst(jt-1),bst(jt),bst(jt+1),thst(jt-1),thst(jt),thst(jt+1))
   21 continue
      
      do 22 jt=2,n2th+4
      bst(jt)=bst(jt)+1./4.*wst(jt)
   22 continue
      bst(1)=bst(n2th+1)
      bst(n2th+5)=bst(5)
   10 continue
      return
      end

!ws*******************************************
      subroutine smth_ps(bst,nt,kk)
      USE DECLARE
      integer nt
      real*8, dimension(nt) :: bst,wst
      include 'mpif.h'

       integer status(mpi_status_size)
      do 10 k=1,kk
      do 21 jt=3,nt-2
      wst(jt)=bst(jt)/2.+(bst(jt+1)+bst(jt-1))*3./16+(bst(jt+2)+bst(jt-2))/16.
   21 continue
      
      do 22 jt=3,nt-2
      bst(jt)=wst(jt)
   22 continue

   10 continue
      return
      end
!ws*******************************************
      subroutine smthEf_dis(kk)
      USE DECLARE
      real*8,dimension(mx,mz,my) :: fsm,wx2,wz2,wy2
      include 'mpif.h'

      dis2(fm1,f0,fp1,xm1,x0,xp1)= &
       ((fp1-f0)/(xp1-x0)-(f0-fm1)/(x0-xm1))

      do 10 k=1,kk
      do 11 m=1,3
      do 21 jy=1,my
      do 21 jz=iz_first+1,iz_last-1
      do 21 jx=ix_first+1,ix_last-1
      if(psi(jx-1,jz).lt.psia .and. psi(jx+1,jz).lt.psia .and. psi(jx,jz-1).lt.psia .and. psi(jx,jz+1).lt.psia) then
      wx2(jx,jz,jy)=dis2(Ef(jx-1,jz,jy,m),Ef(jx,jz,jy,m),Ef(jx+1,jz,jy,m),xx(jx-1),xx(jx),xx(jx+1))
      wz2(jx,jz,jy)=dis2(Ef(jx,jz-1,jy,m),Ef(jx,jz,jy,m),Ef(jx,jz+1,jy,m),zz(jz-1),zz(jz),zz(jz+1))
      endif
   21 continue
   
      do 15 jz=iz_first+1,iz_last-1
      do 15 jx=ix_first+1,ix_last-1
      do 15 jy=iy_first+1,iy_last-1
      wy2(jx,jz,jy)=dis2(Ef(jx,jz,jy-1,m),Ef(jx,jz,jy,m),Ef(jx,jz,jy+1,m),yy(jy-1),yy(jy),yy(jy+1))/xx(jx)

   15 continue

      do 13 jy=iy_first+1,iy_last-1
      do 13 jz=iz_first+1,iz_last-1
      do 13 jx=ix_first+1,ix_last-1
      Ef(jx,jz,jy,m)=Ef(jx,jz,jy,m)+cfsmb(jx,jz)*dxx*(wx2(jx,jz,jy)+wz2(jx,jz,jy)+wy2(jx,jz,jy))      
!           +(.5*(1.+dtanh(pi/2.-thetati(jx))))/20.*w(jx,jz,jy) &
!           +(.5*(1.+dtanh(pi/2.-thetate(jx))))/20.*w(jx,jz,jy)
   13 continue
   11 continue

      call valb3_atlastgrid_r1p0_v1(Ef)
      call mpi_transfersm(Ef(:,:,:,:),3)
   10 continue
      return
      end
!ws*******************************************
      subroutine smthEf_d2f(kk)
      USE DECLARE
      real*8,dimension(mx,mz,my) :: fsm,wx2,wz2,wy2
      include 'mpif.h'
! second-order diffusion
      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)

      dif2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)**2*(fp1-f0)-(xp1-x0)**2*(fm1-f0))/(xp1-xm1)/x0
!  d2f2= d2f / dx2  with second-order accuracy central difference
      d2f2(fm1,f0,fp1,xm1,x0,xp1)= &
       2*((fp1-f0)/(xp1-x0)-(f0-fm1)/(x0-xm1))/(xp1-xm1)

      d1f2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(fp1-f0) &
         -(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)
      do 10 k=1,kk
      do 11 m=1,3
      do 21 jy=1,my
      do 21 jz=iz_first+1,iz_last-1
      do 21 jx=ix_first+1,ix_last-1
      if(psi(jx-1,jz).lt.psia .and. psi(jx+1,jz).lt.psia .and. psi(jx,jz-1).lt.psia .and. psi(jx,jz+1).lt.psia) then
      wx2(jx,jz,jy)=d2f2(Ef(jx-1,jz,jy,m),Ef(jx,jz,jy,m),Ef(jx+1,jz,jy,m),xx(jx-1),xx(jx),xx(jx+1)) &
                   +d1f2(Ef(jx-1,jz,jy,m),Ef(jx,jz,jy,m),Ef(jx+1,jz,jy,m),xx(jx-1),xx(jx),xx(jx+1))/xx(jx)
      wz2(jx,jz,jy)=d2f2(Ef(jx,jz-1,jy,m),Ef(jx,jz,jy,m),Ef(jx,jz+1,jy,m),zz(jz-1),zz(jz),zz(jz+1))
      endif
   21 continue
   
      do 15 jz=iz_first+1,iz_last-1
      do 15 jx=ix_first+1,ix_last-1
      do 15 jy=iy_first+1,iy_last-1
      wy2(jx,jz,jy)=d2f2(Ef(jx,jz,jy-1,m),Ef(jx,jz,jy,m),Ef(jx,jz,jy+1,m),yy(jy-1),yy(jy),yy(jy+1))/xx(jx)**2

   15 continue

      do 13 jy=iy_first+1,iy_last-1
      do 13 jz=iz_first+1,iz_last-1
      do 13 jx=ix_first+1,ix_last-1
      Ef(jx,jz,jy,m)=Ef(jx,jz,jy,m)+cfsmb(jx,jz)*1.e-4*(wx2(jx,jz,jy)+wz2(jx,jz,jy)+wy2(jx,jz,jy))      
!           +(.5*(1.+dtanh(pi/2.-thetati(jx))))/20.*w(jx,jz,jy) &
!           +(.5*(1.+dtanh(pi/2.-thetate(jx))))/20.*w(jx,jz,jy)
   13 continue
   11 continue

      call valb3_atlastgrid_r1p0_v1(Ef)
      call mpi_transfersm(Ef(:,:,:,:),3)
   10 continue
      return
      end
!ws*******************************************
      subroutine smthEf(kk)
      USE DECLARE
      real*8,dimension(mx,mz,my) :: fsm,wx2,wz2,wy2
      include 'mpif.h'
! second-order diffusion
      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)

      dif2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)**2*(fp1-f0)-(xp1-x0)**2*(fm1-f0))/(xp1-xm1)/x0

      do 10 k=1,kk
      do 11 m=1,3
      do 21 jy=1,my
      do 21 jz=iz_first+1,iz_last-1
      do 21 jx=ix_first+1,ix_last-1
      if(psi(jx-1,jz).lt.psia .and. psi(jx+1,jz).lt.psia .and. psi(jx,jz-1).lt.psia .and. psi(jx,jz+1).lt.psia) then
      wx2(jx,jz,jy)=difc(Ef(jx-1,jz,jy,m),Ef(jx,jz,jy,m),Ef(jx+1,jz,jy,m),xx(jx-1),xx(jx),xx(jx+1)) &
                   +dif2(Ef(jx-1,jz,jy,m),Ef(jx,jz,jy,m),Ef(jx+1,jz,jy,m),xx(jx-1),xx(jx),xx(jx+1))
      wz2(jx,jz,jy)=difc(Ef(jx,jz-1,jy,m),Ef(jx,jz,jy,m),Ef(jx,jz+1,jy,m),zz(jz-1),zz(jz),zz(jz+1))
      endif
   21 continue
   
      do 15 jz=iz_first+1,iz_last-1
      do 15 jx=ix_first+1,ix_last-1
      do 15 jy=iy_first+1,iy_last-1
      wy2(jx,jz,jy)=difc(Ef(jx,jz,jy-1,m),Ef(jx,jz,jy,m),Ef(jx,jz,jy+1,m),yy(jy-1),yy(jy),yy(jy+1))

   15 continue

      do 13 jy=iy_first+1,iy_last-1
      do 13 jz=iz_first+1,iz_last-1
      do 13 jx=ix_first+1,ix_last-1
      Ef(jx,jz,jy,m)=Ef(jx,jz,jy,m)+cfsmb(jx,jz)*(wx2(jx,jz,jy)+wz2(jx,jz,jy)+wy2(jx,jz,jy))      
!           +(.5*(1.+dtanh(pi/2.-thetati(jx))))/20.*w(jx,jz,jy) &
!           +(.5*(1.+dtanh(pi/2.-thetate(jx))))/20.*w(jx,jz,jy)
   13 continue
   11 continue

      call valb3_atlastgrid_r1p0_v1(Ef)
      call mpi_transfersm(Ef(:,:,:,:),3)
   10 continue
      return
      end

!ws*******************************************
      subroutine smthEf_dis_v2(kk)
      USE DECLARE
        
        real*8 coeff_xm, coeff_xp, coeff_zm, coeff_zp, coeff_ym, coeff_yp
        real*8 coeff_total
        real*8, dimension(mx,mz,my,3) :: Ef_smooth
       include 'mpif.h' 
!        coeff_smooth = 0.6!large coefficient means heavily smooth

      do 10 k=1,kk
  
      do 21 jy=iy_first+2,iy_last-2
      do 21 jz=iz_first+2,iz_last-2
      do 21 jx=ix_first+2,ix_last-2
      if(psi(jx-1,jz).lt.psia .and. psi(jx+1,jz).lt.psia .and. psi(jx,jz-1).lt.psia .and. psi(jx,jz+1).lt.psia) then
            coeff_xm = 1.0/(xx(jx)-xx(jx-1))
            coeff_xp = 1.0/(xx(jx+1)-xx(jx))
            coeff_zm = 1.0/(zz(jz)-zz(jz-1))
            coeff_zp = 1.0/(zz(jz+1)-zz(jz))
            coeff_ym = 1.0/(xx(jx)*(yy(jy)-yy(jy-1)))
            coeff_yp = 1.0/(xx(jx)*(yy(jy+1)-yy(jy)))
            coeff_total = coeff_xm+coeff_xp+coeff_zm+coeff_zp+coeff_ym+coeff_yp
        Ef_smooth(jx,jz,jy,:) = (1.0 - cfsmb(jx,jz))*Ef(jx,jz,jy,:) &
                              + cfsmb(jx,jz)*(coeff_xm*Ef(jx-1,jz,jy,:) + coeff_xp*Ef(jx+1,jz,jy,:) &
                              +               coeff_zm*Ef(jx,jz-1,jy,:) + coeff_zp*Ef(jx,jz+1,jy,:) &
                              +               coeff_ym*Ef(jx,jz,jy-1,:) + coeff_yp*Ef(jx,jz,jy+1,:))/coeff_total
      endif
   21 continue
       
        Ef = Ef_smooth
      call valb3_atlastgrid_r1p0_v1(Ef)
      call mpi_transfersm(Ef(:,:,:,:),3)   
   10 continue
      return
      end
!ws*******************************************
      subroutine smthxzy_dis_v2(ws,mm,kk)
      USE DECLARE
        integer mm
        real*8 coeff_xm, coeff_xp, coeff_zm, coeff_zp, coeff_ym, coeff_yp
        real*8 coeff_total
        real*8, dimension(mx,mz,my,mm) :: ws, ws_smooth
       include 'mpif.h' 
!        coeff_smooth = 0.6!large coefficient means heavily smooth

      do 10 k=1,kk
  
      do 21 jy=iy_first+2,iy_last-2
      do 21 jz=iz_first+2,iz_last-2
      do 21 jx=ix_first+2,ix_last-2
      if(psi(jx-1,jz).lt.psia .and. psi(jx+1,jz).lt.psia .and. psi(jx,jz-1).lt.psia .and. psi(jx,jz+1).lt.psia) then
            coeff_xm = 1.0/(xx(jx)-xx(jx-1))
            coeff_xp = 1.0/(xx(jx+1)-xx(jx))
            coeff_zm = 1.0/(zz(jz)-zz(jz-1))
            coeff_zp = 1.0/(zz(jz+1)-zz(jz))
            coeff_ym = 1.0/(xx(jx)*(yy(jy)-yy(jy-1)))
            coeff_yp = 1.0/(xx(jx)*(yy(jy+1)-yy(jy)))
            coeff_total = coeff_xm+coeff_xp+coeff_zm+coeff_zp+coeff_ym+coeff_yp
        ws_smooth(jx,jz,jy,:) = (1.0 - cfsmb(jx,jz))*ws(jx,jz,jy,:) &
                                  + cfsmb(jx,jz)*(coeff_xm*ws(jx-1,jz,jy,:) + coeff_xp*ws(jx+1,jz,jy,:) &
                                  +               coeff_zm*ws(jx,jz-1,jy,:) + coeff_zp*ws(jx,jz+1,jy,:) &
                                  +               coeff_ym*ws(jx,jz,jy-1,:) + coeff_yp*ws(jx,jz,jy+1,:))/coeff_total
      endif
   21 continue
       
        ws = ws_smooth
!      call valb3_atlastgrid_r1p0_v1(Ef)
      call mpi_transfersm(ws,mm)   
   10 continue
      return
      end
!ws*******************************************
      subroutine smthxzy_dis(fsm,kk)
      USE DECLARE
      real*8,dimension(mx,mz,my) :: fsm,wx2,wz2,wy2
      include 'mpif.h'

      dis2(fm1,f0,fp1,xm1,x0,xp1)= &
       ((fp1-f0)/(xp1-x0)-(f0-fm1)/(x0-xm1))
      
      do 10 k=1,kk
      do 21 jy=1,my
      do 21 jz=iz_first+1,iz_last-1
      do 21 jx=ix_first+1,ix_last-1
      if(psi(jx-1,jz).lt.psia .and. psi(jx+1,jz).lt.psia .and. psi(jx,jz-1).lt.psia .and. psi(jx,jz+1).lt.psia) then
      wx2(jx,jz,jy)=dis2(fsm(jx-1,jz,jy),fsm(jx,jz,jy),fsm(jx+1,jz,jy),xx(jx-1),xx(jx),xx(jx+1))
      wz2(jx,jz,jy)=dis2(fsm(jx,jz-1,jy),fsm(jx,jz,jy),fsm(jx,jz+1,jy),zz(jz-1),zz(jz),zz(jz+1))
      endif
   21 continue
   
      do 15 jz=iz_first+1,iz_last-1
      do 15 jx=ix_first+1,ix_last-1
      do 15 jy=iy_first+1,iy_last-1
      wy2(jx,jz,jy)=dis2(fsm(jx,jz,jy-1),fsm(jx,jz,jy),fsm(jx,jz,jy+1),yy(jy-1),yy(jy),yy(jy+1))/xx(jx)

   15 continue

      do 13 jy=iy_first+1,iy_last-1
      do 13 jz=iz_first+1,iz_last-1
      do 13 jx=ix_first+1,ix_last-1
      fsm(jx,jz,jy)=fsm(jx,jz,jy)+cfsmb(jx,jz)*dxx*(wx2(jx,jz,jy)+wz2(jx,jz,jy)+wy2(jx,jz,jy))      
!           +(.5*(1.+dtanh(pi/2.-thetati(jx))))/20.*w(jx,jz,jy) &
!           +(.5*(1.+dtanh(pi/2.-thetate(jx))))/20.*w(jx,jz,jy)
   13 continue

   10 continue
      return
      end

!ws*******************************************
      subroutine smthxzy(fsm,kk)
      USE DECLARE
      real*8,dimension(mx,mz,my) :: fsm,wx2,wz2,wy2
      include 'mpif.h'
! second-order diffusion
      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)

!      difxz(fm1m1,fp1m1,fm1p1,fp1p1,xm1,x0,xp1,zm1,z0,zp1)= &
!       2.*( fp1p1-fm1p1-fp1m1-fm1m1)/((xp1-xm1)*(zp1-zm1)) &
!       *((xp1-x0)*(x0-xm1)*(zp1-z0)*(z0-zm1))**0.5
!
!      difxy(fm1,f0,fp1,xm1,x0,xp1,dyy)= &
!       (2.*((xm1-x0)/(xp1-x0)*(fp1-f0)-(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)/x0-f0/x0**2) &
!       *((xp1-x0)*(x0-xm1))**0.5*dyy
!
!      difzy(fm1,f0,fp1,zm1,z0,zp1,dyy,x0)= &
!       (2.*((zm1-z0)/(zp1-z0)*(fp1-f0)-(zp1-z0)/(zm1-z0)*(fm1-f0))/(zm1-zp1)/x0) &
!       *((zp1-z0)*(z0-zm1))**0.5*dyy


       integer status(mpi_status_size)
!      do 10 k=1,kk
!      do 12 jy=1,my
!      do 21 jz=jzamin,jzamax
!      do jx=jxam(jz)+1,jxap(jz)-1
!      wx2(jx,jz,jy)=difc(fsm(jx-1,jz,jy),fsm(jx,jz,jy),fsm(jx+1,jz,jy),xx(jx-1),xx(jx),xx(jx+1))
!      enddo
!   21 continue
!   
!      do 22 jx=jxamin,jxamax
!      do jz=jzam(jx)+1,jzap(jx)-1
!      wz2(jx,jz,jy)=difc(fsm(jx,jz-1,jy),fsm(jx,jz,jy),fsm(jx,jz+1,jy),zz(jz-1),zz(jz),zz(jz+1))
!      enddo
!   22 continue
!   12 continue
      do 10 k=1,kk
      do 21 jy=1,my
      do 21 jz=iz_first+1,iz_last-1
      do 21 jx=ix_first+1,ix_last-1
      if(psi(jx-1,jz).lt.psia .and. psi(jx+1,jz).lt.psia .and. psi(jx,jz-1).lt.psia .and. psi(jx,jz+1).lt.psia) then
      wx2(jx,jz,jy)=difc(fsm(jx-1,jz,jy),fsm(jx,jz,jy),fsm(jx+1,jz,jy),xx(jx-1),xx(jx),xx(jx+1))
      wz2(jx,jz,jy)=difc(fsm(jx,jz-1,jy),fsm(jx,jz,jy),fsm(jx,jz+1,jy),zz(jz-1),zz(jz),zz(jz+1))
      endif
   21 continue
   
      do 15 jz=iz_first+1,iz_last-1
      do 15 jx=ix_first+1,ix_last-1
      do 15 jy=iy_first+1,iy_last-1
      wy2(jx,jz,jy)=difc(fsm(jx,jz,jy-1),fsm(jx,jz,jy),fsm(jx,jz,jy+1),yy(jy-1),yy(jy),yy(jy+1))

   15 continue

      do 13 jy=iy_first+1,iy_last-1
      do 13 jz=iz_first+1,iz_last-1
      do 13 jx=ix_first+1,ix_last-1
      fsm(jx,jz,jy)=fsm(jx,jz,jy)+cfsmb(jx,jz)*(wx2(jx,jz,jy)+wz2(jx,jz,jy)+wy2(jx,jz,jy))      
!           +(.5*(1.+dtanh(pi/2.-thetati(jx))))/20.*w(jx,jz,jy) &
!           +(.5*(1.+dtanh(pi/2.-thetate(jx))))/20.*w(jx,jz,jy)
   13 continue

   10 continue
      return
      end

!ws*******************************************
      subroutine smthe4(fsm,kk)
      USE DECLARE
      real*8,dimension(mx,mz,my) :: fsm,wx2,wz2,wy2
      include 'mpif.h'
! second-order diffusion
      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)

!      difxz(fm1m1,fp1m1,fm1p1,fp1p1,xm1,x0,xp1,zm1,z0,zp1)= &
!       2.*( fp1p1-fm1p1-fp1m1-fm1m1)/((xp1-xm1)*(zp1-zm1)) &
!       *((xp1-x0)*(x0-xm1)*(zp1-z0)*(z0-zm1))**0.5
!
!      difxy(fm1,f0,fp1,xm1,x0,xp1,dyy)= &
!       (2.*((xm1-x0)/(xp1-x0)*(fp1-f0)-(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)/x0-f0/x0**2) &
!       *((xp1-x0)*(x0-xm1))**0.5*dyy
!
!      difzy(fm1,f0,fp1,zm1,z0,zp1,dyy,x0)= &
!       (2.*((zm1-z0)/(zp1-z0)*(fp1-f0)-(zp1-z0)/(zm1-z0)*(fm1-f0))/(zm1-zp1)/x0) &
!       *((zp1-z0)*(z0-zm1))**0.5*dyy


       integer status(mpi_status_size)
      do 10 k=1,kk
      do 12 jy=1,my
      do 21 jz=iz_first+1,iz_last-1
      do 21 jx=ix_first+1,ix_last-1
      if(psi(jx-1,jz).lt.psia .and. psi(jx+1,jz).lt.psia .and. psi(jx,jz-1).lt.psia .and. psi(jx,jz+1).lt.psia) then
      wx2(jx,jz,jy)=difc(fsm(jx-1,jz,jy),fsm(jx,jz,jy),fsm(jx+1,jz,jy),xx(jx-1),xx(jx),xx(jx+1))
      wz2(jx,jz,jy)=difc(fsm(jx,jz-1,jy),fsm(jx,jz,jy),fsm(jx,jz+1,jy),zz(jz-1),zz(jz),zz(jz+1))
      endif
   21 continue
   12 continue

      do 15 jz=iz_first+1,iz_last-1
      do 15 jx=ix_first+1,ix_last-1
      do 15 jy=iy_first+1,iy_last-1
      wy2(jx,jz,jy)=difc(fsm(jx,jz,jy-1),fsm(jx,jz,jy),fsm(jx,jz,jy+1),yy(jy-1),yy(jy),yy(jy+1))
   15 continue

      do 13 jy=iy_first+1,iy_last-1
      do 13 jz=iz_first+1,iz_last-1
      do 13 jx=ix_first+1,ix_last-1
      fsm(jx,jz,jy)=fsm(jx,jz,jy)+1./12*(wx2(jx,jz,jy)+wz2(jx,jz,jy)+wy2(jx,jz,jy))      
!           +(.5*(1.+dtanh(pi/2.-thetati(jx))))/20.*w(jx,jz,jy) &
!           +(.5*(1.+dtanh(pi/2.-thetate(jx))))/20.*w(jx,jz,jy)
   13 continue

   10 continue
      return
      end

 !ws*******************************************
      subroutine avrgt(fsm,kk)
      USE DECLARE
      real*8,dimension(mx,mz,my) :: fsm,wx2,wz2,wy2
      include 'mpif.h'
! second-order diffusion
      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)

!      difxz(fm1m1,fp1m1,fm1p1,fp1p1,xm1,x0,xp1,zm1,z0,zp1)= &
!       2.*( fp1p1-fm1p1-fp1m1-fm1m1)/((xp1-xm1)*(zp1-zm1)) &
!       *((xp1-x0)*(x0-xm1)*(zp1-z0)*(z0-zm1))**0.5
!
!      difxy(fm1,f0,fp1,xm1,x0,xp1,dyy)= &
!       (2.*((xm1-x0)/(xp1-x0)*(fp1-f0)-(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)/x0-f0/x0**2) &
!       *((xp1-x0)*(x0-xm1))**0.5*dyy
!
!      difzy(fm1,f0,fp1,zm1,z0,zp1,dyy,x0)= &
!       (2.*((zm1-z0)/(zp1-z0)*(fp1-f0)-(zp1-z0)/(zm1-z0)*(fm1-f0))/(zm1-zp1)/x0) &
!       *((zp1-z0)*(z0-zm1))**0.5*dyy


       integer status(mpi_status_size)
      do 10 k=1,kk
      do 12 jy=1,my
      do 21 jz=iz_first+1,iz_last-1
      do 21 jx=ix_first+1,ix_last-1
      if(psi(jx-1,jz).lt.psia .and. psi(jx+1,jz).lt.psia .and. psi(jx,jz-1).lt.psia .and. psi(jx,jz+1).lt.psia) then
      wx2(jx,jz,jy)=difc(fsm(jx-1,jz,jy),fsm(jx,jz,jy),fsm(jx+1,jz,jy),xx(jx-1),xx(jx),xx(jx+1))
      wz2(jx,jz,jy)=difc(fsm(jx,jz-1,jy),fsm(jx,jz,jy),fsm(jx,jz+1,jy),zz(jz-1),zz(jz),zz(jz+1))
      endif
   21 continue
   12 continue

      do 15 jz=iz_first+1,iz_last-1
      do 15 jx=ix_first+1,ix_last-1
      do 15 jy=iy_first+1,iy_last-1
      wy2(jx,jz,jy)=difc(fsm(jx,jz,jy-1),fsm(jx,jz,jy),fsm(jx,jz,jy+1),yy(jy-1),yy(jy),yy(jy+1))
   15 continue

      do 13 jy=iy_first+1,iy_last-1
      do 13 jz=iz_first+1,iz_last-1
      do 13 jx=ix_first+1,ix_last-1
      fsm(jx,jz,jy)=fsm(jx,jz,jy)+(1-caf)/1024.0*(wx2(jx,jz,jy)+wz2(jx,jz,jy)+wy2(jx,jz,jy))      
!           +(.5*(1.+dtanh(pi/2.-thetati(jx))))/20.*w(jx,jz,jy) &
!           +(.5*(1.+dtanh(pi/2.-thetate(jx))))/20.*w(jx,jz,jy)
   13 continue

   10 continue
      return
      end
!ws********************************************************************
!ws********************************************************************
subroutine getnp2 ( px, py, x, y, nr, lcell, lnext, xmin, ymin, &
  dx, dy, np, dsq )
!
!***********************************************************************
!
!! GETNP2 seeks the closest unmarked node to a point.
!
!
!  Discussion:
!
!    GETNP2 uses the cell method to find the closest unmarked node NP
!    to a specified point P, given a set of N nodes and the data structure 
!    defined by STORE2.
!
!    NP is then marked by negating LNEXT(NP).  Thus, the closest M nodes to
!    P may be determined by a sequence of M calls to this routine.  
!
!    If the point P is itself actually a node K, and you want to find the
!    nearest point to P that is not node K, then you must be sure to mark
!    node K before calling.
!
!    The search is begun in the cell containing or closest to P and proceeds 
!    outward in rectangular layers until all cells which contain points 
!    within distance R of P have been searched.  R is the distance from P to 
!    the first unmarked node encountered, or infinite if no unmarked nodes
!    are present.
!
!    Input parameters other than LNEXT are not altered by this routine.  
!    With the exception of ( PX, PY ) and the signs of LNEXT elements, 
!    these parameters should be unaltered from their values on output 
!    from subroutine STORE2.
!
!  Modified:
!
!    10 July 1999
!
!  Author:
!
!    Robert Renka,
!    University of North Texas
!
!  Parameters:
!
!    Input, real PX, PY, the (X,Y) coordinates of the point P whose
!    nearest unmarked neighbor is to be found.
!
!    Input, real X(N), Y(N), the coordinates of the nodes at which
!    data has been supplied.
!
!    Input, integer NR, the number of rows and columns in the cell grid.
!    NR must be at least 1.
!
!    Input, integer LCELL(NR,NR), array of nodal indices associated
!    with cells.
!
!    Input/output, integer LNEXT(N), contains next-node indices ( or their 
!    negatives ).  On return, if the output value of NP is nonzero, then
!    LNEXT(NP) will be negative.
!
!    Input, real XMIN, YMIN, DX, DY, the minimum nodal X, Y coordinates,
!    and the X, Y dimensions of a cell.  DX and DY must be positive.
!
!    Output, integer NP, the index into the vectors X and Y of the nearest
!    unmarked node to the point P.  NP will be 0 if all nodes are marked 
!    or if the values of NR, DX, DY are illegal.  LNEXT(NP) will be less
!    than 0 if NP is nonzero (this marks node NP as being used now).
!
!    Output, real DSQ, if NP is nonzero, then DSQ is the squared distance
!    between P and node NP.
!
!  Local Parameters:
!
!    first = true iff the first unmarked node has yet to be encountered,
!
!    imin,imax,jmin,jmax = cell indices defining the range of the search,
!
!    delx,dely = px-xmin and py-ymin,
!
!    i0,j0 = cell containing or closest to P,
!
!    i1,i2,j1,j2 = cell indices of the layer whose intersection with 
!    the range defined by imin,...,jmax is currently being searched.
!
  implicit none
!
  integer nr
!
  real*8 delx
  real*8 dely
  real*8 dsq
  real*8 dx
  real*8 dy
  logical first
  integer i
  integer i0
  integer i1
  integer i2
  integer imax
  integer imin
  integer j
  integer j0
  integer j1
  integer j2
  integer jmax
  integer jmin
  integer l
  integer lcell(nr,nr)
  integer lmin
  integer ln
  integer lnext(*)
  integer np
  real*8 px
  real*8 py
  real*8 r
  real*8 rsmin
  real*8 rsq
  real*8 x(*)
  real*8 xmin
  real*8 xp
  real*8 y(*)
  real*8 ymin
  real*8 yp
!
  xp = px
  yp = py
!
!  Test for invalid input parameters.
!
  if ( nr < 1 .or. dx <= 0.0E+00 .or. dy <= 0.0E+00 ) then
    np = 0
    dsq = 0.0E+00
  end if
!
!  Initialize parameters:
!
  first = .true.
  imin = 1
  imax = nr
  jmin = 1
  jmax = nr
  delx = xp - xmin
  dely = yp - ymin

  i0 = int ( delx / dx ) + 1
  i0 = max ( i0, 1 )
  i0 = min ( i0, nr )

  j0 = int ( dely / dy ) + 1
  j0 = max ( j0, 1 )
  j0 = min ( j0, nr )

  i1 = i0
  i2 = i0
  j1 = j0
  j2 = j0
!
!  Outer loop on layers, inner loop on layer cells, excluding
!  those outside the range (imin,imax) x (jmin,jmax).
!
1 continue

  do j = j1, j2

    if ( j > jmax ) go to 7
    if ( j < jmin ) go to 6

    do i = i1, i2

      if ( i > imax ) go to 6
      if ( i < imin ) go to 5

      if ( j /= j1 .and. j /= j2 .and. i /= i1 .and. i /= i2 ) then
        go to 5
      end if
!
!  Search cell (i,j) for unmarked nodes l.
!
      l = lcell(i,j)

      if ( l > 0 ) then
!
!  Loop on nodes in cell (i,j).
!
2       continue

        ln = lnext(l)
!
!  Node L is the first unmarked neighbor of P encountered.
!
!  Initialize lmin to the current candidate for np, and
!  rsmin to the squared distance from p to lmin.  imin,
!  imax, jmin, and jmax are updated to define the smal-
!  lest rectangle containing a circle of radius r =
!  sqrt(rsmin) centered at p, and contained in (1,nr) x
!  (1,nr) (except that, if p is outside the rectangle
!  defined by the nodes, it is possible that imin .gt.
!  nr, imax < 1, jmin > nr, or jmax < 1).
!
        if ( ln >= 0 ) then

          rsq = ( x(l) - xp )**2 + ( y(l) - yp )**2

          if ( first ) then

            lmin = l
            rsmin = rsq
            r = sqrt ( rsmin )

            imin = int ( ( delx - r ) / dx ) + 1
            imin = max ( imin, 1 )

            imax = int ( ( delx + r ) / dx ) + 1
            imax = min ( imax, nr )

            jmin = int ( ( dely - r ) / dy ) + 1
            jmin = max ( jmin, 1 )

            jmax = int ( ( dely + r ) / dy ) + 1
            jmax = min ( jmax, nr )

            first = .false.

          else

            if ( rsq < rsmin ) then
              lmin = l
              rsmin = rsq
            end if
 
          end if

        end if

        if ( abs ( ln ) /= l ) then
          l = abs ( ln )
          go to 2
        end if

      end if

5     continue

    end do

6   continue

  end do
!
!  Test for termination of loop on cell layers.
!
7 continue

  if ( i1 > imin .or. i2 < imax .or. j1 > jmin .or. j2 < jmax ) then

    i1 = i1 - 1
    i2 = i2 + 1
    j1 = j1 - 1
    j2 = j2 + 1
    go to 1

  end if

  if ( first ) then
    np = 0
    dsq = 0.0E+00
  else
    np = lmin
    dsq = rsmin
    lnext(lmin) = -lnext(lmin)
  end if

  return
end
subroutine givens ( a, b, c, s )
!
!***********************************************************************
!
!! GIVENS constructs a Givens plane rotation.
!
!
!  Discussion:
!
!    The transformation has the form of a 2 by 2 matrix G(C,S):
!
!      (   C  S )
!      ( - S  C )
!
!    where C*C + S*S = 1, which zeroes the second entry of the
!    the column vector ( A, B ) when C and S are properly chosen.
!    A call to GIVENS is normally followed by a call to ROTATE
!    which computes the product of G(C,S) with a 2 by N matrix.
!
!  Modified:
!
!    10 July 1999
!
!  Parameters:
!
!    Input/output, real*8 A, B.
!
!    On input, A and B define the 2-vector whose second entry (B) is
!    to be annihilated by a Givens rotation.
!
!    On output, A has been overwritten by a value
!      R = +/- SQRT ( A*A + B*B )
!    and B has been overwritten by a value Z which allows C
!    and S to be recovered as:
!
!      if | Z | <= 1, then
!        C = SQRT ( 1 - Z*Z ), 
!        S = Z
!      else if | Z | > 1 then
!        C = 1 / Z, 
!        S = SQRT ( 1 - C*C ).
!
!    Output, real*8 C, S, the components of the Givens transformation, 
!    which may be computed by:
!
!      C = +/- A / SQRT ( A*A + B*B )
!      S = +/- B / SQRT ( A*A + B*B )
!
!  Local parameters:
!
!  r =        c*a + s*b = +/-sqrt(a*a+b*b)
!  u,v =   variables used to scale a and b for computing r
!
!  abs(a) > abs(b)
!
!  Note that r has the sign of a, c > 0, and s has
!  sign(a)*sign(b).
!
  implicit none
!
  real*8 a
  real*8 b
  real*8 c
  real*8 r
  real*8 s
  real*8 u
  real*8 v
!
  if ( abs ( a ) > abs ( b ) ) then

    u = 2.0E+00 * a
    v = b / u
    r = sqrt ( 0.25E+00 + v * v ) * u
    c = a / r
    s = 2.0E+00 * v * c
    b = s
    a = r
!
!  abs(a) <= abs(b)
!
!  Store r in a.
!  Note that r has the sign of b, s > 0, and c has sign(a)*sign(b).
!
  else if ( b /= 0.0E+00 ) then

    u = 2.0E+00 * b
    v = a / u
    a = sqrt ( 0.25E+00 + v * v ) * u
    s = b / a
    c = 2.0E+00 * v * s

    if ( c /= 0.0E+00 ) then
      b = 1.0E+00 / c
    else
      b = 1.0E+00
    end if
!
!  a = b = 0.
!
  else

    c = 1.0E+00
    s = 0.0E+00

  end if

  return
end
subroutine qs2grd ( px, py, n, x, y, f, nr, lcell, lnext, xmin, &
  ymin, dx, dy, rmax, rsq, a, q, qx, qy, ier )
!
!***********************************************************************
!
!! QS2GRD evaluates the interpolant and its first spatial derivatives.
!
!
!  Discussion:
!
!    QS2GRD computes the value and the gradient at the point (PX,PY) 
!    of the interpolatory function Q, defined by QSHEP2 for a given set
!    of scattered data.  Q(X,Y) is a weighted sum of quadratic
!    nodal functions.
!
!    Input parameters are not altered by this subroutine.  The parameters 
!    other than PX and PY should be input unaltered from their values 
!    on output from QSHEP2.  This subroutine should not be called if a 
!    nonzero error flag was returned by QSHEP2.
!
!  Modified:
!
!    10 July 1999
!
!  Author:
!
!    Robert Renka,
!    University of North Texas
!
!  Parameters:
!
!    Input, real*8 PX, PY, the coordinates of the point at which the
!    interpolant and its derivatives are to be evaluated.
!
!    Input, integer N, the number of nodes and data values which
!    are to be interpolated.  N must be at least 6. 
!
!    Input, real*8 X(N), Y(N), the coordinates of the nodes at which
!    data has been supplied.
!
!    Input, real*8 F(N), the data values at the nodes.
!
!    Input, integer NR, the number of rows and columns in the cell 
!    grid.  Refer to subroutine STORE2 for details.  NR must be at least 1.
!
!    Input, integer LCELL(NR,NR), array of nodal indices associated
!    with cells.
!
!    Input, integer LNEXT(N), contains next-node indices.
!
!    Input, real*8 XMIN, YMIN, DX, DY, the minimum nodal X, Y coordinates,
!    and the X, Y dimensions of a cell.  Computed by QSHEP2.
!
!    Input, real*8 RMAX, the square root of the largest element in RSQ,
!    the maximum radius of influence.  Computed by QSHEP2.
!
!    Input, real*8 RSQ(N), the squared radii which enter into the weights 
!    defining the interpolant Q.  Computed by QSHEP2.
!
!    Input, real*8 A(5,N), the coefficients for the nodal functions 
!    defining the interpolant Q.  Computed by QSHEP2.
!
!    Output, real*8 Q, QX, QY, the value of the interpolant, and its derivatives
!    with respect to X and Y, at (PX,PY).
!
!    Output, integer IER, error indicator.
!    0, if no errors were encountered.
!    1, if N, NR, DX, DY or RMAX is invalid.
!    2, if no errors were encountered but (PX,PY) is not within the 
!       radius R(K) for any node K and thus Q = QX = QY = 0.
!
  implicit none
!
  integer n
  integer nr
!
  real*8 a(5,n)
  real*8 delx
  real*8 dely
  real*8 ds
  real*8 dx
  real*8 dy
  real*8 f(n)
  integer i
  integer ier
  integer imax
  integer imin
  integer j
  integer jmax
  integer jmin
  integer k
  integer kp
  integer lcell(nr,nr)
  integer lnext(n)
  real*8 px
  real*8 py
  real*8 q
  real*8 qk
  real*8 qkx
  real*8 qky
  real*8 qx
  real*8 qy
  real*8 rd
  real*8 rds
  real*8 rmax
  real*8 rs
  real*8 rsq(n)
  real*8 sw
  real*8 swq
  real*8 swqx
  real*8 swqy
  real*8 sws
  real*8 swx
  real*8 swy
  real*8 t
  real*8 w
  real*8 wx
  real*8 wy
  real*8 x(n)
  real*8 xmin
  real*8 xp
  real*8 y(n)
  real*8 ymin
  real*8 yp
!
  xp = px
  yp = py

  if ( n < 6 ) then
    ier = 1
    return
  else if ( nr < 1 ) then
    ier = 1
    return
  else if ( dx <= 0.0E+00 ) then
    ier = 1
    return
  else if ( dy <= 0.0E+00 ) then
    ier = 1
    return
  else if ( rmax < 0.0E+00 ) then
    ier = 1
    return
  end if
!
!  Set imin, imax, jmin, and jmax to cell indices defining
!  the range of the search for nodes whose radii include P.
!  The cells which must be searched are those inter-
!  sected by (or contained in) a circle of radius rmax
!  centered at p.
!
  imin = int ( ( xp - xmin - rmax ) / dx ) + 1
  imin = max ( imin, 1 )

  imax = int ( ( xp - xmin + rmax ) / dx ) + 1
  imax = min ( imax, nr )

  jmin = int ( ( yp - ymin - rmax ) / dy ) + 1
  jmin = max ( jmin, 1 )

  jmax = int ( ( yp - ymin + rmax ) / dy ) + 1
  jmax = min ( jmax, nr )
!
!  Test for no cells within the circle of radius RMAX.
!
  if ( imin > imax .or. jmin > jmax ) then
    q = 0.0E+00
    qx = 0.0E+00
    qy = 0.0E+00
    ier = 2
    return
  end if
!
!  Q = swq/sw = sum(w(k)*q(k))/sum(w(k)) where the sum is
!  from k = 1 to n, q(k) is the quadratic nodal function,
!  and w(k) = ((r-d)+/(r*d))**2 for radius r(k) and distance d(k).  Thus
!
!    qx = (swqx*sw - swq*swx)/sw**2  and
!    qy = (swqy*sw - swq*swy)/sw**2
!
!  where swqx and swx are partial derivatives with respect
!  to x of swq and sw, respectively.  swqy and swy are 
!  defined similarly.
!
  sw = 0.0E+00
  swx = 0.0E+00
  swy = 0.0E+00
  swq = 0.0E+00
  swqx = 0.0E+00
  swqy = 0.0E+00
!
!  Outer loop on cells (I,J).
!
  do j = jmin, jmax

    do i = imin, imax

      k = lcell(i,j)
!
!  Inner loop on nodes K.
!
      if ( k /= 0 ) then

1       continue

        delx = xp - x(k)
        dely = yp - y(k)
        ds = delx * delx + dely * dely
        rs = rsq(k)

        if ( ds == 0.0E+00 ) then
          q = f(k)
          qx = a(4,k)
          qy = a(5,k)
          ier = 0
          return
        end if

        if ( ds < rs ) then

          rds = rs * ds
          rd = sqrt ( rds )
          w = ( rs + ds - rd - rd ) / rds
          t = 2.0E+00 * ( rd - rs ) / ( ds * rds )
          wx = delx * t
          wy = dely * t
          qkx = 2.0E+00 * a(1,k) * delx + a(2,k) * dely
          qky = a(2,k) * delx + 2.0E+00 * a(3,k) * dely
          qk = ( qkx * delx + qky * dely ) / 2.0E+00
          qkx = qkx + a(4,k)
          qky = qky + a(5,k)
          qk = qk + a(4,k) * delx + a(5,k) * dely + f(k)
          sw = sw + w
          swx = swx + wx
          swy = swy + wy
          swq = swq + w * qk
          swqx = swqx + wx * qk + w * qkx
          swqy = swqy + wy * qk + w * qky

        end if

        kp = k
        k = lnext(kp)

        if ( k /= kp ) then
          go to 1
        end if

      end if

    end do

  end do
!
!  SW = 0 if and only if P is not within the radius R(K) for any node K.
!
  if ( sw /= 0.0E+00 ) then

    q = swq / sw
    sws = sw * sw
    qx = ( swqx * sw - swq * swx ) / sws
    qy = ( swqy * sw - swq * swy ) / sws
    ier = 0

  else

    q = 0.0E+00
    qx = 0.0E+00
    qy = 0.0E+00
    ier = 2

  end if

  return
end
subroutine qshep2 ( n, x, y, f, nq, nw, nr, lcell, lnext, xmin, &
  ymin, dx, dy, rmax, rsq, a, ier )
!
!***********************************************************************
!
!! QSHEP2 computes an interpolant to scattered data in the plane.
!
!
!  Discussion:
!
!    QSHEP2 computes a set of parameters A and RSQ defining a smooth, 
!    once continuously differentiable, bi-variate function Q(X,Y) which 
!    interpolates given data values F at scattered nodes (X,Y).  
!
!    The interpolant function Q(X,Y) may be evaluated at an arbitrary point 
!    by passing the parameters A and RSQ to the function QS2VAL.  The
!    first derivatives dQdX(X,Y) and dQdY(X,Y) may be evaluated by 
!    subroutine QS2GRD.
!
!    The interpolation scheme is a modified quadratic Shepard method:
!
!      Q = ( W(1) * Q(1) + W(2) * Q(2) + .. + W(N) * Q(N) ) 
!        / ( W(1)        + W(2)        + .. + W(N) )
!
!    for bivariate functions W(K) and Q(K).  The nodal functions are given by
!
!      Q(K)(X,Y) = 
!          F(K)
!        + A(1,K) * ( X - X(K) )**2 
!        + A(2,K) * ( X - X(K) ) * ( Y - Y(K) )
!        + A(3,K) * ( Y - Y(K) )**2 
!        + A(4,K) * ( X - X(K) )
!        + A(5,K) * ( Y - Y(K) ).
!
!    Thus, Q(K) is a quadratic function which interpolates the
!    data value at node K.  Its coefficients A(*,K) are obtained
!    by a weighted least squares fit to the closest NQ data
!    points with weights similar to W(K).  Note that the radius
!    of influence for the least squares fit is fixed for each
!    K, but varies with K.
!
!    The weights are taken to be
!
!      W(K)(X,Y) = ( (R(K)-D(K))+ / R(K) * D(K) )**2
!
!    where (R(K)-D(K))+ = 0 if R(K) <= D(K) and D(K)(X,Y) is
!    the euclidean distance between (X,Y) and (X(K),Y(K)).  The
!    radius of influence R(K) varies with K and is chosen so
!    that NW nodes are within the radius.  Note that W(K) is
!    not defined at node (X(K),Y(K)), but Q(X,Y) has limit F(K)
!    as (X,Y) approaches (X(K),Y(K)).
!
!  Author:
!
!    Robert Renka,
!    University of North Texas
!
!  Parameters:
!
!    Input, integer N, the number of nodes (X,Y) at which data values
!    are given.  N must be at least 6.
!
!    Input, real X(N), Y(N), the coordinates of the nodes at which
!    data has been supplied.
!
!    Input, real F(N), the data values.
!
!    Input, integer NQ, the number of data points to be used in the least
!    squares fit for coefficients defining the nodal functions Q(K).  
!    A highly recommended value is NQ = 13.  
!    NQ must be at least 5, and no greater than the minimum of 40 and N-1.
!
!    Input, integer NW, the number of nodes within (and defining) the radii
!    of influence R(K) which enter into the weights W(K).  For N 
!    sufficiently large, a recommended value is NW = 19.   NW must be
!    at least 1, and no greater than the minimum of 40 and N-1.
!
!    Input, integer NR, the number of rows and columns in the cell grid 
!    defined in subroutine STORE2.  A rectangle containing the nodes 
!    is partitioned into cells in order to increase search efficiency.  
!    NR = SQRT(N/3) is recommended.  NR must be at least 1.
!
!    Output, integer LCELL(NR,NR), array of nodal indices associated
!    with cells.
!
!    Output, integer LNEXT(N), contains next-node indices ( or their 
!    negatives ).
!
!    Output, real XMIN, YMIN, DX, DY, the minimum nodal X, Y coordinates,
!    and the X, Y dimensions of a cell.
!
!    Output, real RMAX, the square root of the largest element in RSQ,
!    the maximum radius of influence.
!
!    Output, real RSQ(N), the squared radii which enter into the weights 
!    defining the interpolant Q.
!
!    Output, real A(5,N), the coefficients for the nodal functions 
!    defining the interpolant Q.
!
!    Output, integer IER, error indicator.
!    0, if no errors were encountered.
!    1, if N, NQ, NW, or NR is out of range.
!    2, if duplicate nodes were encountered.
!    3, if all nodes are collinear.
!
!  Local parameters:
!
! av =        root-mean-square distance between k and the
!             nodes in the least squares system (unless
!             additional nodes are introduced for stabil-
!             ity).      the first 3 columns of the matrix
!             are scaled by 1/avsq, the last 2 by 1/av
! avsq =      av*av
! b =         transpose of the augmented regression matrix
! c =         first component of the plane rotation used to
!             zero the lower triangle of b**t -- computed
!             by subroutine givens
! ddx,ddy =   local variables for dx and dy
! dmin =      minimum of the magnitudes of the diagonal
!             elements of the regression matrix after
!             zeros are introduced below the diagonal
! DTOL =      tolerance for detecting an ill-conditioned
!             system.  the system is accepted when dmin
!             >= DTOL
! fk =        data value at node k -- f(k)
! i =         index for a, b, and npts
! ib =        do-loop index for back solve
! ierr =      error flag for the call to store2
! irow =      row index for b
! j =         index for a and b
! jp1 =       j+1
! k =         nodal function index and column index for a
! lmax =      maximum number of npts elements (must be con-
!             sistent with the dimension statement above)
! lnp =       current length of npts
! neq =       number of equations in the least squares fit
! nn,nnr =    local copies of n and nr
! np =        npts element
! npts =      array containing the indices of a sequence of
!             nodes to be used in the least squares fit
!             or to compute rsq.  the nodes are ordered
!             by distance from k and the last element
!             (usually indexed by lnp) is used only to
!             determine rq, or rsq(k) if nw > nq
! nqwmax =    max(nq,nw)
! rq =        radius of influence which enters into the
!             weights for q(k) (see subroutine setup2)
! rs =        squared distance between k and npts(lnp) --
!             used to compute rq and rsq(k)
! rsmx =      maximum rsq element encountered
! rsold =     squared distance between k and npts(lnp-1) --
!             used to compute a relative change in rs
!             between succeeding npts elements
! RTOL =      tolerance for detecting a sufficiently large
!             relative change in rs.  if the change is
!             not greater than RTOL, the nodes are
!             treated as being the same distance from k
! rws =       current value of rsq(k)
! s =         second component of the plane givens rotation
! SF =        marquardt stabilization factor used to damp
!             out the first 3 solution components (second
!             partials of the quadratic) when the system
!             is ill-conditioned.  as SF increases, the
!             fitting function approaches a linear
! sum2 =      sum of squared euclidean distances between
!             node k and the nodes used in the least
!             squares fit (unless additional nodes are
!             added for stability)
! t =         temporary variable for accumulating a scalar
!             product in the back solve
! xk,yk =     coordinates of node k -- x(k), y(k)
! xmn,ymn =   local variables for xmin and ymin
!
  implicit none
!
  integer n
  integer nr
!
  real*8 a(5,n)
  real*8 av
  real*8 avsq
  real*8 b(6,6)
  real*8 c
  real*8 ddx
  real*8 ddy
  real*8 dmin
  real*8, parameter :: dtol = 0.01E+00 !0.01E+00
  real*8 dx
  real*8 dy
  real*8 f(n)
  real*8 fk
  integer i
  integer ier
  integer ierr
  integer irow
  integer j
  integer jp1
  integer k
  integer lcell(nr,nr)
  integer lmax
  integer lnext(n)
  integer lnp
  integer neq
  integer nn
  integer nnr
  integer np
  integer npts(40)
  integer nq
  integer nqwmax
  integer nw
  real*8 rmax
  real*8 rq
  real*8 rs
  real*8 rsmx
  real*8 rsold
  real*8 rsq(n)
  real*8, parameter :: rtol = 1.0E-05
  real*8 rws
  real*8 s
  real*8, parameter :: SF = 1.0E+00
  real*8 sum2
  real*8 t
  real*8 x(n)
  real*8 xk
  real*8 xmin
  real*8 xmn
  real*8 y(n)
  real*8 yk
  real*8 ymin
  real*8 ymn
!
  nn = n
  nnr = nr
  nqwmax = max ( nq, nw )
  lmax = min ( 40, n-1 )

  if ( 5 > nq ) then
    ier = 1
    return
  else if ( 1 > nw ) then
    ier = 1
    return
  else if ( nqwmax > lmax ) then
    ier = 1
    return
  else if ( nr < 1 ) then
    ier = 1
    return
  end if
!
!  Create the cell data structure, and initialize RSMX.
!
  call store2 ( nn, x, y, nnr, lcell, lnext, xmn, ymn, ddx, ddy, ierr )

  if ( ierr /= 0 ) then
    xmin = xmn
    ymin = ymn
    dx = ddx
    dy = ddy
    ier = 3
    return
  end if

  rsmx = 0.0E+00
!
!  Outer loop on node K.
!
  do k = 1, nn

    xk = x(k)
    yk = y(k)
    fk = f(k)
!
!  Mark node K to exclude it from the search for nearest neighbors.
!
    lnext(k) = - lnext(k)
!
!  Initialize for loop on NPTS.
!
    rs = 0.0E+00
    sum2 = 0.0E+00
    rws = 0.0E+00
    rq = 0.0E+00
    lnp = 0
!
!  Compute NPTS, LNP, RWS, NEQ, RQ, and AVSQ.
!
1   continue

    sum2 = sum2 + rs

    if ( lnp == lmax ) then
      go to 3
    end if

    lnp = lnp + 1
    rsold = rs

    call getnp2 ( xk, yk, x, y, nnr, lcell, lnext, xmn, ymn, ddx, ddy, np, rs )

    if ( rs == 0.0E+00 ) then
      ier = 2
      return
    end if

    npts(lnp) = np

    if ( ( rs - rsold ) / rs < RTOL ) then
      go to 1
    end if

    if ( rws == 0.0E+00 .and. lnp > nw ) then
      rws = rs
    end if
!
!  RQ = 0 (not yet computed) and lnp > nq.     
!
!  RQ = sqrt(rs) is sufficiently large to (strictly) include nq nodes.  
!
!  The least squares fit will include NEQ = LNP - 1 equations for 
!  5 <= NQ <= NEQ < LMAX <= N-1.
!
    if ( rq == 0.0E+00 .and. lnp > nq ) then
      neq = lnp - 1
      rq = sqrt ( rs )
      avsq = sum2 / real ( neq )
    end if

    if ( lnp > nqwmax ) then
      go to 4
    else
      go to 1
    end if
!
!  All LMAX nodes are included in npts.   RWS and/or rq**2 is
!  (arbitrarily) taken to be 10 percent larger than the
!  distance rs to the last node included.
!
3   continue

    if ( rws == 0.0E+00 ) then
      rws = 1.1E+00 * rs
    end if

    if ( rq == 0.0E+00 ) then
      neq = lmax
      rq = sqrt ( 1.1E+00 * rs )
      avsq = sum2 / real ( neq )
    end if

4   continue
!
!  Store rsq(k), update rsmx if necessary, and compute av.
!
    rsq(k) = rws
    rsmx = max ( rsmx, rws )
    av = sqrt ( avsq )
!
!  Set up the augmented regression matrix (transposed) as the
!  columns of B, and zero out the lower triangle (upper
!  triangle of B) with givens rotations -- qr decomposition
!  with orthogonal matrix Q not stored.
!
    i = 0

5   continue

    i = i + 1
    np = npts(i)
    irow = min ( i, 6 )

    call setup2 ( xk, yk, fk, x(np), y(np), f(np), av, avsq, rq, b(1,irow) )

    if ( i == 1 ) then
      go to 5
    end if

    do j = 1, irow-1
      jp1 = j + 1
      call givens ( b(j,j), b(j,irow), c, s )
      call rotate ( 6-j, c, s, b(jp1,j), b(jp1,irow) )
    end do

    if ( i < neq ) then
      go to 5
    end if
!
!  Test the system for ill-conditioning.
!
    dmin =  min ( abs ( b(1,1) ), abs ( b(2,2) ), abs ( b(3,3) ), &
      abs ( b(4,4) ), abs ( b(5,5) ) )

    if ( dmin * rq >= DTOL ) then
      go to 13
    end if

    if ( neq == lmax ) then
      go to 10
    end if
!
!  Increase RQ and add another equation to the system to improve conditioning.  
!  The number of NPTS elements is also increased if necessary.
!
7   continue

    rsold = rs
    neq = neq + 1

    if ( neq == lmax ) then
      go to 9
    end if
!
!   NEQ < LNP.
!
    if ( neq /= lnp ) then
      np = npts(neq+1)
      rs = ( x(np) - xk )**2 + ( y(np) - yk )**2
      if ( ( rs - rsold ) / rs < rtol ) then
        go to 7
      end if
      rq = sqrt(rs)
      go to 5
    end if
!
!  Add an element to NPTS.
!
    lnp = lnp + 1
    call getnp2 ( xk, yk, x, y, nnr, lcell, lnext, xmn, ymn, ddx, ddy, np, rs )

    if ( np == 0 ) then
      ier = 2
      return
    end if

    npts(lnp) = np

    if ( ( rs - rsold ) / rs < rtol ) then
      go to 7
    end if

    rq = sqrt ( rs )
    go to 5

9   continue

    rq = sqrt ( 1.1E+00 * rs )
    go to 5
!
!  Stabilize the system by damping second partials.  Add multiples of the 
!  first three unit vectors to the first three equations.
!
10  continue

    do i = 1, 3

      b(i,6) = sf

      do j = i+1, 6
        b(j,6) = 0.0E+00
      end do

      do j = i, 5
        jp1 = j + 1
        call givens ( b(j,j), b(j,6), c, s )
        call rotate ( 6-j, c, s, b(jp1,j), b(jp1,6) )
      end do

    end do
!
!  Test the stabilized system for ill-conditioning.
!
    dmin = min ( abs ( b(1,1) ), abs ( b(2,2) ), abs ( b(3,3) ), &
      abs ( b(4,4) ), abs ( b(5,5) ) )

!    if ( dmin * rq < dtol ) then
!      xmin = xmn
!      ymin = ymn
!      dx = ddx
!      dy = ddy
!      ier = 3
!      return
!    end if
!
!  Solve the 5 by 5 triangular system for the coefficients.
!
13  continue

    do i = 5, 1, -1

      t = 0.0E+00

      do j = i+1, 5
        t = t + b(j,i) * a(j,k)
      end do

      a(i,k) = ( b(6,i) - t ) / b(i,i)

    end do
!
!  Scale the coefficients to adjust for the column scaling.
!
    do i = 1, 3
      a(i,k) = a(i,k) / avsq
    end do

    a(4,k) = a(4,k) / av
    a(5,k) = a(5,k) / av
!
!  Unmark K and the elements of NPTS.
!
    lnext(k) = - lnext(k)

    do i = 1, lnp
      np = npts(i)
      lnext(np) = - lnext(np)
    end do

  end do
!
!  No errors encountered.
!
  xmin = xmn
  ymin = ymn
  dx = ddx
  dy = ddy
  rmax = sqrt ( rsmx )
  ier = 0

  return
end
function qs2val ( px, py, n, x, y, f, nr, lcell, lnext, xmin, &
  ymin, dx, dy, rmax, rsq, a )
!
!***********************************************************************
!
!! QS2VAL evaluates the interpolant function at a point.
!
!
!  Discussion:
!
!    QS2VAL returns the value Q(PX,PY) where Q is the weighted sum of 
!    quadratic nodal functions defined by QSHEP2.  If the spatial 
!    derivatives of Q are also desired, call QS2GRD instead.
!
!    Input parameters are not altered by this function.  The
!    parameters other than PX and PY should be input unaltered
!    from their values on output from QSHEP2.  This function
!    should not be called if a nonzero error flag was returned
!    by QSHEP2.
!
!  Modified:
!
!    10 July 1999
!
!  Author:
!
!    Robert Renka,
!    University of North Texas
!
!  Parameters:
!
!    Input, real PX, PY, the (X,Y) coordinates of the point P at
!    which Q is to be evaluated.
!
!    Input, integer N, the number of nodes and data values to be 
!    interpolated.  N must be at least 6.
!
!    Input, real X(N), Y(N), the coordinates of the nodes at which
!    data has been supplied.
!
!    Input, real F(N), the data values at the nodes.
!
!    Input, integer NR, the number of rows and columns in the cell grid.
!    Refer to subroutine STORE2.  NR must be at least 1.
!
!    Input, integer LCELL(NR,NR), the array of nodal indices associated
!    with cells.  Refer to STORE2.
!
!    Input, integer LNEXT(N), the next-node indices.  Refer to STORE2.
!
!    Input, real XMIN, YMIN, DX, DY, the minimum nodal X, Y coordinates,
!    and the X, Y dimensions of a cell.  Computed by QSHEP2.
!
!    Input, real RMAX, the square root of the largest element in RSQ,
!    the maximum radius of influence.  Computed by QSHEP2.
!
!    Input, real RSQ(N), the squared radii which enter into the weights 
!    defining the interpolant Q.  Computed by QSHEP2.
!
!    Input, real A(5,N), the coefficients for the nodal functions 
!    defining the interpolant Q.  Computed by QSHEP2.
!
!    Output, real QS2VAL, the interpolated function value at (PX,PY).
!
  implicit none
!
  integer n
  integer nr
!
  real*8 a(5,n)
  real*8 delx
  real*8 dely
  real*8 dx
  real*8 dy
  real*8 f(n)
  integer i
  integer imax
  integer imin
  integer j
  integer jmax
  integer jmin
  real*8 ds
  integer k
  integer kp
  integer lcell(nr,nr)
  integer lnext(n)
  real*8 px
  real*8 py
  real*8 qs2val
  real*8 rd
  real*8 rds
  real*8 rmax
  real*8 rs
  real*8 rsq(n)
  real*8 sw
  real*8 swq
  real*8 w
  real*8 x(n)
  real*8 xmin
  real*8 xp
  real*8 y(n)
  real*8 ymin
  real*8 yp
!
  xp = px
  yp = py
  qs2val = 0.0E+00

  if ( n < 6  ) then
    return  
  else if ( nr < 1  ) then
    return
  else if ( dx <= 0.0E+00 ) then
    return
  else if ( dy <= 0.0E+00 ) then
    return
  else if ( rmax < 0.0E+00 ) then
    return
  end if
!
!  Set imin, imax, jmin, and jmax to cell indices defining
!  the range of the search for nodes whose radii include
!  p.  The cells which must be searched are those inter-
!  sected by (or contained in) a circle of radius rmax
!  centered at p.
!
  imin = int ( ( xp - xmin - rmax ) / dx ) + 1
  imin = max ( imin, 1 )

  imax = int ( ( xp - xmin + rmax ) / dx ) + 1
  imax = min ( imax, nr )

  jmin = int ( ( yp - ymin - rmax ) / dy ) + 1
  jmin = max ( jmin, 1 )

  jmax = int ( ( yp - ymin + rmax ) / dy ) + 1
  jmax = min ( jmax, nr )
!
!  Test for no cells within the circle of radius RMAX.
!
  if ( imin > imax .or. jmin > jmax ) then
    qs2val = 0.0E+00
    return
  end if
!
!  Accumulate weight values in SW and weighted nodal function
!  values in swq.  the weights are w(k) = ((r-d)+/(r*d))**2
!  for r**2 = rsq(k) and d = distance between p and node k.
!
  sw = 0.0E+00
  swq = 0.0E+00

  do j = jmin, jmax

    do i = imin, imax

      k = lcell(i,j)

      if ( k /= 0 ) then

1       continue

        delx = xp - x(k)
        dely = yp - y(k)
        ds = delx * delx + dely * dely
        rs = rsq(k)

        if ( ds < rs ) then

          if ( ds == 0.0E+00 ) then
            qs2val = f(k)
            return
          end if

          rds = rs * ds
          rd = sqrt ( rds )
          w = ( rs + ds - rd - rd ) / rds
          sw = sw + w

          swq = swq + w * ( f(k) + a(1,k) * delx * delx &
            + a(2,k) * delx * dely + a(3,k) * dely * dely &
            + a(4,k) * delx + a(5,k) * dely )

        end if

        kp = k
        k = lnext(kp)

        if ( k /= kp ) then
          go to 1
        end if

      end if

    end do

  end do
!
!  SW = 0 if and only if P is not within the radius R(K) for any node K.
!
  if ( sw == 0.0E+00 ) then
    qs2val = 0.0E+00
  else
    qs2val = swq / sw
  end if

  return
end
subroutine rotate ( n, c, s, x, y )
!
!***********************************************************************
!
!! ROTATE applies a Givens rotation.
!
!
!  Discussion:
!
!    The rotation has the form:
!
!      (   C  S )
!      ( - S  C )
!
!    and is essentially applied to a 2 by N matrix:
!
!      ( X(1) X(2) ... X(N) )
!      ( Y(1) Y(2) ... Y(N) )
!
!  Modified:
!
!    28 June 1999
!
!  Parameters:
!
!    Input, integer N, the dimension of the vectors.
!
!    Input, real C, S, the cosine and sine entries of the Givens
!    rotation matrix.  These may be determined by subroutine GIVENS.
!
!    Input/output, real X(N), Y(N), the rotated vectors. 
!
  implicit none
!
  integer n
!
  real*8 c
  integer i
  real*8 s
  real*8 x(n)
  real*8 xi
  real*8 y(n)
  real*8 yi
!
  if ( n <= 0 ) then
    return
  else if ( c == 1.0E+00 .and. s == 0.0E+00 ) then
    return
  end if

  do i = 1, n
    xi = x(i)
    yi = y(i)
    x(i) =   c * xi + s * yi
    y(i) = - s * xi + c * yi
  end do

  return
end
subroutine setup2 ( xk, yk, fk, xi, yi, fi, s1, s2, r, row )
!
!***********************************************************************
!
!! SETUP2 sets up a row of the least squares regression matrix.
!
!
!  Discussion:
!
!    SETUP2 sets up the I-th row of an augmented regression matrix for 
!    a weighted least-squares fit of a quadratic function Q(X,Y) to a set 
!    of data values F, where Q(XK,YK) = FK.  
!
!    The first 3 columns are quadratic terms, and are scaled by 1/S2.
!    The fourth and fifth columns represent linear terms, and are scaled 
!    by 1/S1.  
!
!    If D = 0, or D >= R, the weight is
!      0,
!    else if D < R, the weight is 
!      (R-D)/(R*D), 
!    where D is the distance between nodes I and K, and R is a maximum
!    influence distance.
!
!  Modified:
!
!    05 July 1999
!
!  Author:
!
!    Robert Renka,
!    University of North Texas
!
!  Parameters:
!
!    Input, real XK, YK, FK, the coordinates and value of the data
!    at data node K.
!
!    Input, real XI, YI, FI, the coorindates and value of the data
!    at data node I.
!
!    Input, real S1, S2, reciprocals of the scale factors.
!
!    Input, real R, the maximum radius of influence about node K.
!
!    Output, real ROW(6), a row of the augmented regression matrix.
!
  implicit none
!
  real*8 d
  real*8 dx
  real*8 dy
  real*8 fi
  real*8 fk
  integer i
  real*8 r
  real*8 row(6)
  real*8 s1
  real*8 s2
  real*8 w
  real*8 xi
  real*8 xk
  real*8 yi
  real*8 yk
!
  dx = xi - xk
  dy = yi - yk

  d = sqrt ( dx * dx + dy * dy )

  if ( d <= 0.0E+00 .or. d >= r ) then

    row(1:6) = 0.0E+00

  else

    w = ( r - d ) / r / d

    row(1) = dx * dx * w / s2
    row(2) = dx * dy * w / s2
    row(3) = dy * dy * w / s2
    row(4) = dx * w / s1
    row(5) = dy * w / s1
    row(6) = ( fi - fk ) * w

  end if

  return
end
subroutine store2 ( n, x, y, nr, lcell, lnext, xmin, ymin, dx, dy, ier )
!
!***********************************************************************
!
!! STORE2 creates a cell data structure for the scattered data.
!
!
!  Discussion:
!
!    STORE2 is given a set of N arbitrarily distributed nodes in the 
!    plane and creates a data structure for a cell-based method of 
!    solving closest-point problems.  The smallest rectangle containing 
!    all the nodes is partitioned into an NR by NR uniform grid of cells, 
!    and nodes are associated with cells.      
!
!    In particular, the data structure stores the indices of the nodes 
!    contained in each cell.  For a uniform random distribution of nodes, 
!    the nearest node to an arbitrary point can be determined in constant
!    expected time.
!
!  Modified:
!
!    05 July 1999
!
!  Author:
!
!    Robert Renka
!    University of North Texas
!
!  Parameters:
!
!    Input, integer N, the number of data nodes.  N must be at least 2.
!
!    Input, real X(N), Y(N), the coordinates of the data nodes.
!
!    Input, integer NR, the number of rows and columns in the grid.  The
!    cell density, or average number of data nodes per cell, is
!      D = N / ( NR * NR ).
!    A recommended value, based on empirical evidence, is 
!      D = 3. 
!    Hence, the corresponding value of NR is recommended to be about
!      NR = SQRT ( N / 3 ).  
!    NR must be at least 1.
!
!    Output, integer LCELL(NR,NR), an array set up so that LCELL(I,J)
!    contains the index (for X and Y) of the first data node (that is, the
!    data node with smallest index) in the (I,J) cell.  LCELL(I,J) will be 0 if 
!    no data nodes are contained in the (I,J) cell.  The upper right corner of 
!    the (I,J) cell has coordinates 
!      ( XMIN + I * DX, YMIN + J * DY ).
!
!    Output, integer LNEXT(N), an array of next-node indices.  LNEXT(K)
!    contains the index of the next node in the cell which contains node K, 
!    or LNEXT(K) = K if K is the last node in the cell.
!    The data nodes contained in a cell are ordered by their indices.
!    If, for example, cell (I,J) contains nodes 2, 3, and 5 and no others, 
!    then:
!
!      LCELL(I,J) = 2, (index of the first data node)
!
!      LNEXT(2) = 3, 
!      LNEXT(3) = 5,
!      LNEXT(5) = 5.
!
!    Output, real XMIN, YMIN, the X, Y coordinates of the lower left
!    corner of the rectangle defined by the data nodes.  The upper right 
!    corner is ( XMAX, YMAX ), where
!      XMAX = XMIN + NR * DX,
!      YMAX = YMIN + NR * DY.
!
!    Output, real DX, DY, the X and Y dimensions of the individual cells.
!      DX = ( XMAX - XMIN ) / NR
!      DY = ( YMAX - YMIN ) / NR,
!    where XMIN, XMAX, YMIN and YMAX are the extrema of X and Y.
!
!    Output, integer IER, an error indicator.
!    0, if no errors were encountered.
!    1, if N < 2 or NR < 1.
!    2, if DX = 0 or DY = 0.
!
  implicit none
!
  integer n
  integer nr
!
  real*8 dx
  real*8 dy
  integer i
  integer ier
  integer j
  integer k
  integer l
  integer lcell(nr,nr)
  integer lnext(n)
  real*8 x(n)
  real*8 xmax
  real*8 xmin
  real*8 y(n)
  real*8 ymax
  real*8 ymin
!
  ier = 0

  if ( n < 2 ) then
    ier = 1
    return
  end if

  if ( nr < 1 ) then
    ier = 1
    return
  end if
!
!  Compute the dimensions of the (X,Y) rectangle containing all the data nodes.
!
  xmin = minval ( x(1:n) )
  xmax = maxval ( x(1:n) )
  ymin = minval ( y(1:n) )
  ymax = maxval ( y(1:n) )
!
!  Compute the dimensions of a single cell.
!
  dx = ( xmax - xmin ) / real ( nr )
  dy = ( ymax - ymin ) / real ( nr )
!
!  Test for zero area.
!
  if ( dx == 0.0E+00 .or. dy == 0.0E+00 ) then
    ier = 2
    return
  end if
!
!  Initialize LCELL.
!
  do j = 1, nr
    do i = 1, nr
      lcell(i,j) = 0
    end do
  end do
!
!  Loop on nodes, storing indices in LCELL and LNEXT.
!
  do k = n, 1, -1

    i = int ( ( x(k) - xmin ) / dx ) + 1
    i = min ( i, nr )

    j = int ( ( y(k) - ymin ) / dy ) + 1
    j = min ( j, nr )

    l = lcell(i,j)

    if ( l /= 0 ) then
      lnext(k) = l
    else
      lnext(k) = k
    end if

    lcell(i,j) = k

  end do

  return
end
subroutine timestamp ( )
!
!*******************************************************************************
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none
!
  character ( len = 8 ) ampm
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
  character ( len = 5 ) zone
!
  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
!ws********************************************
!ws********************************************
      real*8 function rhom(psival)
      USE DECLARE
      real*8 psival,rsq
!
!   rsq is the normalized poloidal flux
!
      rsq=(psival-psmin)/(psia-psmin)
      if(iden.eq.1) rhom=(1.00000-alphar*rsq**prho)**arho
!    Gaussian density profile
      if(iden.eq.2) rhom=exp(-alphar*rsq)+prho*rsq*(arho-rsq)
! TFTR DT shot 66887 electron density profile
      if(iden.eq.3)rhom=(1.3+5.3*(1.-rsq)*(1.-0.95*rsq*(1.-rsq)))/6.6

      if(iden.eq.5)rhom=(1.00000-alphar*rsq-prho*rsq**3)
      return
      end
!ws********************************************
      real*8 function rhomp(psival)
      USE DECLARE
      real*8 psival,rsq
!
!   rsq is the normalized poloidal flux
!
      rsq=(psival-psmin)/(psia-psmin)
      if(iden.eq.1) rhomp=-arho*(1.00000-alphar*rsq**prho)**(arho-1.)*alphar*prho*rsq**(prho-1)/(psia-psmin)
!    Gaussian density profile
      if(iden.eq.2) rhomp=-alphar*exp(-alphar*rsq)+prho*(arho-2.*rsq)/(psia-psmin)
! TFTR DT shot 66887 electron density profile
      if(iden.eq.3) rhomp=(5.3/6.6)*(-(1.-0.95*rsq*(1.-rsq))-0.95*(1.-rsq)*(1.-2.*rsq))/(psia-psmin)

      if(iden.eq.5) rhomp=-alphar-prho*3*rsq**2/(psia-psmin)
      return
      end

!ws**********************************************
      subroutine mpi_test
      implicit real*8 (b-h,o-z)
      real*8 www
      integer mm
      dimension fxz(3,3)
      include 'mpif.h'
       integer status(mpi_status_size)
      if(nrank.eq.0) then
      do i=1,3
      do j=1,3
      fxz(i,j)=i+j
      enddo
      enddo
      www=1.
      CALL MPI_Send(www, 1, MPI_DOUBLE_PRECISION, 1, 1,  &
		               MPI_COMM_WORLD,ierror )
      CALL MPI_Send(fxz(1:2,:), 2*3, MPI_DOUBLE_PRECISION, 1, 2,  &
		               MPI_COMM_WORLD,ierror )
      endif
      if(nrank.eq.1) then    
      CALL MPI_Recv(www, 1, MPI_DOUBLE_PRECISION, 0, 1,  &
		               MPI_COMM_WORLD,status,ierror )
      CALL MPI_Recv(fxz(1:2,:), 2*3, MPI_DOUBLE_PRECISION, 0, 2,  &
		               MPI_COMM_WORLD,status,ierror )
      endif
      write(*,*) nrank,www,fxz(1,1),fxz(2,2)

      return
      end



!ws:bndry8&3
!****************************************************************
      subroutine map_xz2st(fxz,fst,mm)
      USE DECLARE
      implicit real*8 (b-h,o-z)
      integer mm,im
      dimension fxz(mx,mz,my,mm),fst(n2th+5,mps4:mps,my,mm)
      include 'mpif.h'

       integer status(mpi_status_size)

!
      do 2 js=mpsa-2,mpsa-nda,-1
      if(mts_nrk(js,nrankxz).gt.0) then
      do 21 im=1,mm
      do 21 jy=1,my      
      call interp_xz2ps(fxz(ix_first:ix_last,iz_first:iz_last,jy,im), &
                        xx(ix_first:ix_last),zz(iz_first:iz_last),ix_last-ix_first+1,iz_last-iz_first+1, &
                       fst(itsmin(js,nrankxz):itsmax(js,nrankxz),js,jy,im), &
                        xxs(itsmin(js,nrankxz):itsmax(js,nrankxz),js),zzs(itsmin(js,nrankxz):itsmax(js,nrankxz),js),mts_nrk(js,nrankxz))
      
   21 continue    
!!mpi
!	write(*,*) nrankxz,js,"map1 done" 
         do irecv=1,mrkb
         do isend=1,nsend(irecv,js)
         if(nrankxz.eq.nranksend(irecv,js,isend) .and. nrankxz.ne.nrkb(irecv)) then
         ltmin=ittransmin(irecv,js,isend)
         ltmax=ittransmax(irecv,js,isend)
         if (itbmin(nrkb(irecv))+n2th .le. itsmax(js,nranksend(irecv,js,isend))) then
         ltmin=ltmin+n2th
         ltmax=ltmax+n2th
         endif
         if (itbmax(nrkb(irecv))-n2th .ge. itsmin(js,nranksend(irecv,js,isend))) then
         ltmin=ltmin-n2th
         ltmax=ltmax-n2th
         endif
         do lt=ltmin,ltmax
         CALL MPI_Send(fst(lt,js,1:my,1:mm),my*mm, MPI_DOUBLE_PRECISION, nrkb(irecv)+nrky(nrank)*nprxz, isend,  &
		               MPI_COMM_WORLD,ierror )
         enddo
         endif
         enddo
         enddo

!	write(*,*) nrank,js,"mapsend done" 
      endif

      if(mb_nrk(nrankxz).gt.0) then  
         do irecv=1,mrkb
         do isend=1,nsend(irecv,js)
         if(nrankxz.eq.nrkb(irecv) .and. nrankxz.ne.nranksend(irecv,js,isend)) then
         ltmin=ittransmin(irecv,js,isend)
         ltmax=ittransmax(irecv,js,isend)
         do lt=ltmin,ltmax
         CALL MPI_Recv(fst(lt,js,1:my,1:mm),my*mm, MPI_DOUBLE_PRECISION, nranksend(irecv,js,isend)+nrky(nrank)*nprxz, isend,  &
		               MPI_COMM_WORLD,status,ierror )
         enddo
         endif
         enddo
         enddo
!	write(*,*) nrank,js,"maprecv done" 
       endif


!!mpi
    2 continue 
      return
      end
!****************************************************************
      subroutine smth_st_nrk(fstsm,js,mm,kk)
      USE DECLARE
      implicit real*8 (b-h,o-z)
      integer mm,js,kk,ltmin,ltmax,im
      dimension fstsm(n2th+5,mps4:mps,my,mm),wst(n2th+5)
      include 'mpif.h'
       integer status(mpi_status_size)

      do 11 k=1,kk
      do im=1,mm
      do jy=iy_first+2,iy_last-2
      do jt=itbmin(nrankxz)+2,itbmax(nrankxz)-2
      wst(jt)=fstsm(jt,js,jy,im)/2.+((fstsm(jt+1,js,jy,im)+fstsm(jt-1,js,jy,im))*3./16+(fstsm(jt+2,js,jy,im)+fstsm(jt-2,js,jy,im))/16.&
             +(fstsm(jt,js,jy+1,im)+fstsm(jt,js,jy-1,im))*3./16+(fstsm(jt,js,jy+2,im)+fstsm(jt,js,jy-2,im))/16.)/2.
      enddo
      
      do jt=itbmin(nrankxz)+2,itbmax(nrankxz)-2
      fstsm(jt,js,jy,im)=wst(jt)
      enddo

      enddo
      enddo

         if(inrkb(nrankxz) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrankxz)+1))
         ltmax=itbmin(nrkb1(inrkb(nrankxz)+1))+1
         CALL MPI_Send(fstsm(ltmin:ltmax,js,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrankxz)+1)+nrky(nrank)*nprxz, 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
   
         if(inrkb(nrankxz) .gt. 1) then
         ltmin=itbmin(nrankxz)
         ltmax=itbmin(nrankxz)+1
         CALL MPI_Recv(fstsm(ltmin:ltmax,js,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrankxz)-1)+nrky(nrank)*nprxz, 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrankxz) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         CALL MPI_Send(fstsm(ltmin:ltmax,js,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(1)+nrky(nrank)*nprxz, 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
         
         if(inrkb(nrankxz) .eq. 1) then
         ltmin=itbmin(nrankxz)
         ltmax=itbmin(nrankxz)+1
         CALL MPI_Recv(fstsm(ltmin:ltmax,js,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(mrkb)+nrky(nrank)*nprxz, 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

  !!ws  
         if(inrkb(nrankxz) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrankxz)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrankxz)-1))
         CALL MPI_Send(fstsm(ltmin:ltmax,js,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrankxz)-1)+nrky(nrank)*nprxz, 2,  &
		               MPI_COMM_WORLD,ierror )
         endif

         if(inrkb(nrankxz) .lt. mrkb) then
         ltmin=itbmax(nrankxz)-1
         ltmax=itbmax(nrankxz)
         CALL MPI_Recv(fstsm(ltmin:ltmax,js,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrankxz)+1)+nrky(nrank)*nprxz, 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   
         if(inrkb(nrankxz) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         CALL MPI_Send(fstsm(ltmin:ltmax,js,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(mrkb)+nrky(nrank)*nprxz, 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrankxz) .eq. mrkb) then
         ltmin=itbmax(nrankxz)-1
         ltmax=itbmax(nrankxz)
         CALL MPI_Recv(fstsm(ltmin:ltmax,js,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(1)+nrky(nrank)*nprxz, 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   11 continue
      return
      end

!****************************************************************
      subroutine valbm_atlastgrid_v1(fxz,mm,ibnd)
      USE DECLARE
      implicit real*8 (b-h,o-z)
      integer mm,ibnd
      dimension fxz(mx,mz,my,mm),fst(n2th+5,mps4:mps,my,mm),f1s(mbm_nrk,mps4:mps)
      include 'mpif.h'

       integer status(mpi_status_size)
!
      call map_xz2st(fxz,fst,mm)
      
      if(mb_nrk(nrankxz).gt.0) then 

      do 10 js=mpsa-2,mpsa-4,-1
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),js,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
!ws_smps--------------------------
      call smth_st_nrk(fst,js,mm,3)
!ws_smps--------------------
   10 continue

      if(ibnd==1) then
      is=mpsa-1
      do 31 m=1,mm    
      do 31 jy=1,my      
      do jt=itbmin(nrankxz),itbmax(nrankxz)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
!       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))  
       fst(jt,is,jy,m)=(4*fst(jt,is-1,jy,m)-fst(jt,is-2,jy,m))/3.0

       if(fst(jt,is,jy,m)*fst(jt,is-1,jy,m).le.0) fst(jt,is,jy,m)=0
!       fst(jt,is)=fst(jt,is-1)                   
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   31 continue 
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),is,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
      call smth_st_nrk(fst,is,mm,3)
!ws_smps--------------------
      is=mpsa
      do 41 m=1,mm    
      do 41 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
!       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))     
       fst(jt,is,jy,m)=(4*fst(jt,is-1,jy,m)-fst(jt,is-2,jy,m))/3.0

       if(fst(jt,is,jy,m)*fst(jt,is-1,jy,m).le.0) fst(jt,is,jy,m)=0
!       fst(jt,is)=fst(jt,is-1)                   
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   41 continue
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),is,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
!ws_smps--------------------
       call smth_st_nrk(fst,is,mm,3)
!ws_smps--------------------
      endif !!ibnd==1
      
      if(ibnd==0) then
      is=mpsa
      do 51 m=1,mm    
      do 51 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
       fst(jt,is,jy,m)=0.           
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   51 continue

      is=mpsa-1
      do 52 m=1,mm    
      do 52 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
      call interp1d3l(fst(jt,is-3,jy,m),fst(jt,is-2,jy,m),fst(jt,is-1,jy,m),fst(jt,is+1,jy,m), &
                        ps(is-3),ps(is-2),ps(is-1),ps(is+1),ps(is),fst(jt,is,jy,m))                     
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   52 continue

      call smth_st_nrk(fst,is,mm,3)
      endif !!ibnd==0 

      do 3 m=1,mm    
      do 3 jy=1,my           
      do ikb=1,mb_nrk(nrankxz)
        i=itb_nrk(ikb,nrankxz)
        ib=ib_nrk(ikb,nrankxz)
        do js=mpsa,mpsa-3,-1     
        call interp1d3l(fst(i-2,js,jy,m),fst(i-1,js,jy,m),fst(i,js,jy,m),fst(i+1,js,jy,m), &
                        thst(i-2),thst(i-1),thst(i),thst(i+1),thb(ib),f1s(ikb,js))
        enddo
        js=mpsa
        call interp1d3l(f1s(ikb,js),f1s(ikb,js-1),f1s(ikb,js-2),f1s(ikb,js-3), &
                        ps(js),ps(js-1),ps(js-2),ps(js-3),psb(ib),fxz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,m))
!        fxz(jbx(ib),jbz(ib))=(asm_b(ib)*x1s(ib,js-1)+bsm_b(ib)*x1s(ib,js-2)+csm_b(ib)*x1s(ib,js-3)) &
!                 /(asm_b(ib)+bsm_b(ib)+csm_b(ib))
      enddo
    3 continue 
      endif !!mb_nrk(nrankxz).gt.0
           
!      if(ibnd==2) then
!      
!      do ika1=1,ma1_nrk(nrankxz)
!        i=ita1_nrk(ika1,nrankxz)
!        ia1=ia1_nrk(ika1,nrankxz)
!        do js=mpsa,mpsa-3,1
!        call interp1d3l(fst(i-2,js),fst(i-1,js),fst(i,js),fst(i+1,js), &
!                        thst(i-2),thst(i-1),thst(i),thst(i+1),tha1(ia1),x1s(ika1,js,jy))
!        enddo
!        js=mpsa
!        call interp1d3l(x1s(ika1,js,jy),x1s(ika1,js-1,jy),x1s(ika1,js-2,jy),x1s(ika1,js-3,jy), &
!                        ps(js),ps(js-1),ps(js-2),ps(js-3),psa1(ia1),fxz(jxa1_nrk(ika1,nrankxz),jza1_nrk(ika1,nrankxz)))
!!        fxz(jxa1(ia1),jza1(ia1))=(asm_a1(ia1)*x1s(ia1,js-1)+bsm_a1(ia1)*x1s(ia1,js-2)+csm_a1(ia1)*x1s(ia1,js-3)) &
!!                  /(asm_a1(ia1)+bsm_a1(ia1)+csm_a1(ia1)) 
!      enddo
!      endif
      return
      end
!****************************************************************
      subroutine valb8_atlastgrid_r0p1_v2(f8xz)
      USE DECLARE
      implicit real*8 (b-h,o-z)
      real*8 vx1st,vz1st,bx1st,bz1st
      integer is
      dimension f8xz(mx,mz,my,8),fst(n2th+5,mps4:mps,my,8),f1s(mbm_nrk,mps4:mps),fsxz(mbm_nrk,8) !
      include 'mpif.h'

       integer status(mpi_status_size)

      call map_xz2st(f8xz,fst,8)
         
      if(mb_nrk(nrankxz).gt.0) then      

      do jy=1,my
      do js=mpsa-2,mpsa-4,-1
      do jt=itbmin(nrankxz),itbmax(nrankxz)
      vx1st=fst(jt,js,jy,3)
      vz1st=fst(jt,js,jy,5)
      bx1st=fst(jt,js,jy,6)
      bz1st=fst(jt,js,jy,8)
      fst(jt,js,jy,3)=vx1st*wbxr(jt,js)+vz1st*wbzr(jt,js)
      fst(jt,js,jy,5)=vx1st*wbxt(jt,js)+vz1st*wbzp(jt,js)
      fst(jt,js,jy,6)=bx1st*wbxr(jt,js)+bz1st*wbzr(jt,js)
      fst(jt,js,jy,8)=bx1st*wbxt(jt,js)+bz1st*wbzp(jt,js)
       
!      vrst=fst(jt,js,jy,3)*dcos(thst(jt))+fst(jt,js,jy,5)*dsin(thst(jt))
!      vpst=-fst(jt,js,jy,3)*dsin(thst(jt))+fst(jt,js,jy,5)*dcos(thst(jt))
!      brst=fst(jt,js,jy,6)*dcos(thst(jt))+fst(jt,js,jy,8)*dsin(thst(jt))
!      bpst=-fst(jt,js,jy,6)*dsin(thst(jt))+fst(jt,js,jy,8)*dcos(thst(jt))
!      fst(jt,js,jy,3)=vrst
!      fst(jt,js,jy,5)=vpst
!      fst(jt,js,jy,6)=brst
!      fst(jt,js,jy,8)=bpst
      enddo
      enddo
      enddo    

      do 10 js=mpsa-2,mpsa-4,-1
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),js,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
!ws_smps--------------------------
      call smth_st_nrk(fst,js,8,3)
!ws_smps--------------------
   10 continue

      do 61 m=1,8

      if(m.eq.2 .or. m.eq.3 .or. m.eq.6) then !!(m=2,3,6) ibnd==0
      is=mpsa 
      do 51 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
       fst(jt,is,jy,m)=0.           
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   51 continue

      is=mpsa-1 
      do 52 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
        call interp1d3l(fst(jt,is-3,jy,m),fst(jt,is-2,jy,m),fst(jt,is-1,jy,m),fst(jt,is+1,jy,m), &
                        ps(is-3),ps(is-2),ps(is-1),ps(is+1),ps(is),fst(jt,is,jy,m))                   
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   52 continue

      call smth_st_nrk(fst(:,:,:,m),is,1,3)

      else  !!(m=1,4,5,7,8,) bnd==1  
      is=mpsa-1    
      do 31 jy=1,my      
      do jt=itbmin(nrankxz),itbmax(nrankxz)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
!       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is)) 
       fst(jt,is,jy,m)=(4*fst(jt,is-1,jy,m)-fst(jt,is-2,jy,m))/3.0
!       fst(jt,is)=fst(jt,is-1) 
       if(fst(jt,is,jy,m)*fst(jt,is-1,jy,m).le.0) fst(jt,is,jy,m)=0.                  
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   31 continue 
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),is,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
      call smth_st_nrk(fst(:,:,:,m),is,1,3)
!ws_smps--------------------
      is=mpsa
 
      do 41 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
!       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))
       fst(jt,is,jy,m)=(4*fst(jt,is-1,jy,m)-fst(jt,is-2,jy,m))/3.0
!       fst(jt,is)=fst(jt,is-1) 
       if(fst(jt,is,jy,m)*fst(jt,is-1,jy,m).le.0) fst(jt,is,jy,m)=0.                  
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   41 continue
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),is,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
!ws_smps--------------------
      call smth_st_nrk(fst(:,:,:,m),is,1,3)
!	write(*,*) nrankxz,m,"smth done"
!ws_smps--------------------
      endif !!(m=3,6)ibnd==0
   61 continue
      
      do 3 jy=1,my           
      do ikb=1,mb_nrk(nrankxz)
        i=itb_nrk(ikb,nrankxz)
        ib=ib_nrk(ikb,nrankxz) 
        do m=1,8 
        do js=mpsa,mpsa-3,-1     
        call interp1d3l(fst(i-2,js,jy,m),fst(i-1,js,jy,m),fst(i,js,jy,m),fst(i+1,js,jy,m), &
                        thst(i-2),thst(i-1),thst(i),thst(i+1),thb(ib),f1s(ikb,js))
        enddo
        js=mpsa
        call interp1d3l(f1s(ikb,js),f1s(ikb,js-1),f1s(ikb,js-2),f1s(ikb,js-3), &
                        ps(js),ps(js-1),ps(js-2),ps(js-3),psb(ib),fsxz(ikb,m))
        enddo
!        fxz(jbx(ib),jbz(ib))=(asm_b(ib)*x1s(ib,js-1)+bsm_b(ib)*x1s(ib,js-2)+csm_b(ib)*x1s(ib,js-3)) &
!                 /(asm_b(ib)+bsm_b(ib)+csm_b(ib))
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,1)=fsxz(ikb,1)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,2)=fsxz(ikb,2)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,3)=fsxz(ikb,3)*wbrx(ib)+fsxz(ikb,5)*wbpx(ib)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,4)=fsxz(ikb,4)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,5)=fsxz(ikb,3)*wbrz(ib)+fsxz(ikb,5)*wbpz(ib)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,6)=fsxz(ikb,6)*wbrx(ib)+fsxz(ikb,8)*wbpx(ib)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,7)=fsxz(ikb,7)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,8)=fsxz(ikb,6)*wbrz(ib)+fsxz(ikb,8)*wbpz(ib)
      enddo
    3 continue 
      endif !!mb_nrk(nrankxz).gt.0
           
!      if(ibnd==2) then
!      
!      do ika1=1,ma1_nrk(nrankxz)
!        i=ita1_nrk(ika1,nrankxz)
!        ia1=ia1_nrk(ika1,nrankxz)
!        do js=mpsa,mpsa-3,1
!        call interp1d3l(fst(i-2,js),fst(i-1,js),fst(i,js),fst(i+1,js), &
!                        thst(i-2),thst(i-1),thst(i),thst(i+1),tha1(ia1),x1s(ika1,js,jy))
!        enddo
!        js=mpsa
!        call interp1d3l(x1s(ika1,js,jy),x1s(ika1,js-1,jy),x1s(ika1,js-2,jy),x1s(ika1,js-3,jy), &
!                        ps(js),ps(js-1),ps(js-2),ps(js-3),psa1(ia1),fxz(jxa1_nrk(ika1,nrankxz),jza1_nrk(ika1,nrankxz)))
!!        fxz(jxa1(ia1),jza1(ia1))=(asm_a1(ia1)*x1s(ia1,js-1)+bsm_a1(ia1)*x1s(ia1,js-2)+csm_a1(ia1)*x1s(ia1,js-3)) &
!!                  /(asm_a1(ia1)+bsm_a1(ia1)+csm_a1(ia1)) 
!      enddo
!      endif
      return
      end
!****************************************************************
!****************************************************************
      subroutine valb8_atlastgrid(f8xz)
      USE DECLARE
      implicit real*8 (b-h,o-z)
      real*8 vx1st,vz1st,bx1st,bz1st
      integer is
      dimension f8xz(mx,mz,my,8),fst(n2th+5,mps4:mps,my,8),f1s(mbm_nrk,mps4:mps),fsxz(mbm_nrk,8) !
      include 'mpif.h'

       integer status(mpi_status_size)

      call map_xz2st(f8xz,fst,8)
         
      if(mb_nrk(nrankxz).gt.0) then      

      do jy=1,my
      do js=mpsa-2,mpsa-4,-1
      do jt=itbmin(nrankxz),itbmax(nrankxz)
      vx1st=fst(jt,js,jy,3)
      vz1st=fst(jt,js,jy,5)
      bx1st=fst(jt,js,jy,6)
      bz1st=fst(jt,js,jy,8)
      fst(jt,js,jy,3)=vx1st*wbxr(jt,js)+vz1st*wbzr(jt,js)
      fst(jt,js,jy,5)=vx1st*wbxt(jt,js)+vz1st*wbzp(jt,js)
      fst(jt,js,jy,6)=bx1st*wbxr(jt,js)+bz1st*wbzr(jt,js)
      fst(jt,js,jy,8)=bx1st*wbxt(jt,js)+bz1st*wbzp(jt,js)
       
!      vrst=fst(jt,js,jy,3)*dcos(thst(jt))+fst(jt,js,jy,5)*dsin(thst(jt))
!      vpst=-fst(jt,js,jy,3)*dsin(thst(jt))+fst(jt,js,jy,5)*dcos(thst(jt))
!      brst=fst(jt,js,jy,6)*dcos(thst(jt))+fst(jt,js,jy,8)*dsin(thst(jt))
!      bpst=-fst(jt,js,jy,6)*dsin(thst(jt))+fst(jt,js,jy,8)*dcos(thst(jt))
!      fst(jt,js,jy,3)=vrst
!      fst(jt,js,jy,5)=vpst
!      fst(jt,js,jy,6)=brst
!      fst(jt,js,jy,8)=bpst
      enddo
      enddo
      enddo    

      do 10 js=mpsa-2,mpsa-4,-1
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),js,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
!ws_smps--------------------------
      call smth_st_nrk(fst,js,8,3)
!ws_smps--------------------
   10 continue

      do 61 m=1,8

      if(lbndxfix(m)) then !!(m=3,6) ibnd==0
      is=mpsa 
      do 51 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
       fst(jt,is,jy,m)=0.           
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   51 continue

      is=mpsa-1 
      do 52 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
        !call interp1d3l(fst(jt,is-3,jy,m),fst(jt,is-2,jy,m),fst(jt,is-1,jy,m),fst(jt,is+1,jy,m), &
        !                ps(is-3),ps(is-2),ps(is-1),ps(is+1),ps(is),fst(jt,is,jy,m))                   
           call interp1d2l(fst(jt,is-2,jy,m),fst(jt,is-1,jy,m),fst(jt,is+1,jy,m), &
                        ps(is-2),ps(is-1),ps(is+1),ps(is),fst(jt,is,jy,m))    
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   52 continue

      call smth_st_nrk(fst(:,:,:,m),is,1,3)

      else  !!(m=1,2,4,5,7,8,) bnd==1  
      is=mpsa-1    
      do 31 jy=1,my      
      do jt=itbmin(nrankxz),itbmax(nrankxz)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
!       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is)) 
       fst(jt,is,jy,m)=(4*fst(jt,is-1,jy,m)-fst(jt,is-2,jy,m))/3.0
!       fst(jt,is)=fst(jt,is-1) 
       if(fst(jt,is,jy,m)*fst(jt,is-1,jy,m).le.0) fst(jt,is,jy,m)=0.                  
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   31 continue 
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),is,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
      call smth_st_nrk(fst(:,:,:,m),is,1,3)
!ws_smps--------------------
      is=mpsa
 
      do 41 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
!       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))
       fst(jt,is,jy,m)=(4*fst(jt,is-1,jy,m)-fst(jt,is-2,jy,m))/3.0
!       fst(jt,is)=fst(jt,is-1) 
       if(fst(jt,is,jy,m)*fst(jt,is-1,jy,m).le.0) fst(jt,is,jy,m)=0.                  
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   41 continue
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),is,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
!ws_smps--------------------
      call smth_st_nrk(fst(:,:,:,m),is,1,3)
!	write(*,*) nrankxz,m,"smth done"
!ws_smps--------------------
      endif !!(m=3,6)ibnd==0
   61 continue
      
      do 3 jy=1,my           
      do ikb=1,mb_nrk(nrankxz)
        i=itb_nrk(ikb,nrankxz)
        ib=ib_nrk(ikb,nrankxz) 
        do m=1,8 
        do js=mpsa,mpsa-3,-1     
        call interp1d3l(fst(i-2,js,jy,m),fst(i-1,js,jy,m),fst(i,js,jy,m),fst(i+1,js,jy,m), &
                        thst(i-2),thst(i-1),thst(i),thst(i+1),thb(ib),f1s(ikb,js))
        enddo
        js=mpsa
        call interp1d3l(f1s(ikb,js),f1s(ikb,js-1),f1s(ikb,js-2),f1s(ikb,js-3), &
                        ps(js),ps(js-1),ps(js-2),ps(js-3),psb(ib),fsxz(ikb,m))
        enddo
!        fxz(jbx(ib),jbz(ib))=(asm_b(ib)*x1s(ib,js-1)+bsm_b(ib)*x1s(ib,js-2)+csm_b(ib)*x1s(ib,js-3)) &
!                 /(asm_b(ib)+bsm_b(ib)+csm_b(ib))
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,1)=fsxz(ikb,1)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,2)=fsxz(ikb,2)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,3)=fsxz(ikb,3)*wbrx(ib)+fsxz(ikb,5)*wbpx(ib)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,4)=fsxz(ikb,4)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,5)=fsxz(ikb,3)*wbrz(ib)+fsxz(ikb,5)*wbpz(ib)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,6)=fsxz(ikb,6)*wbrx(ib)+fsxz(ikb,8)*wbpx(ib)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,7)=fsxz(ikb,7)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,8)=fsxz(ikb,6)*wbrz(ib)+fsxz(ikb,8)*wbpz(ib)
      enddo
    3 continue 
      endif !!mb_nrk(nrankxz).gt.0
           
!      if(ibnd==2) then
!      
!      do ika1=1,ma1_nrk(nrankxz)
!        i=ita1_nrk(ika1,nrankxz)
!        ia1=ia1_nrk(ika1,nrankxz)
!        do js=mpsa,mpsa-3,1
!        call interp1d3l(fst(i-2,js),fst(i-1,js),fst(i,js),fst(i+1,js), &
!                        thst(i-2),thst(i-1),thst(i),thst(i+1),tha1(ia1),x1s(ika1,js,jy))
!        enddo
!        js=mpsa
!        call interp1d3l(x1s(ika1,js,jy),x1s(ika1,js-1,jy),x1s(ika1,js-2,jy),x1s(ika1,js-3,jy), &
!                        ps(js),ps(js-1),ps(js-2),ps(js-3),psa1(ia1),fxz(jxa1_nrk(ika1,nrankxz),jza1_nrk(ika1,nrankxz)))
!!        fxz(jxa1(ia1),jza1(ia1))=(asm_a1(ia1)*x1s(ia1,js-1)+bsm_a1(ia1)*x1s(ia1,js-2)+csm_a1(ia1)*x1s(ia1,js-3)) &
!!                  /(asm_a1(ia1)+bsm_a1(ia1)+csm_a1(ia1)) 
!      enddo
!      endif
      return
      end
!****************************************************************
      subroutine valb8_atlastgrid_r0p1_v1(f8xz)
      USE DECLARE
      implicit real*8 (b-h,o-z)
      real*8 vx1st,vz1st,bx1st,bz1st
      integer is
      dimension f8xz(mx,mz,my,8),fst(n2th+5,mps4:mps,my,8),f1s(mbm_nrk,mps4:mps),fsxz(mbm_nrk,8) !
      include 'mpif.h'

       integer status(mpi_status_size)

      call map_xz2st(f8xz,fst,8)
         
      if(mb_nrk(nrankxz).gt.0) then      

      do jy=1,my
      do js=mpsa-2,mpsa-4,-1
      do jt=itbmin(nrankxz),itbmax(nrankxz)
      vx1st=fst(jt,js,jy,3)
      vz1st=fst(jt,js,jy,5)
      bx1st=fst(jt,js,jy,6)
      bz1st=fst(jt,js,jy,8)
      fst(jt,js,jy,3)=vx1st*wbxr(jt,js)+vz1st*wbzr(jt,js)
      fst(jt,js,jy,5)=vx1st*wbxt(jt,js)+vz1st*wbzp(jt,js)
      fst(jt,js,jy,6)=bx1st*wbxr(jt,js)+bz1st*wbzr(jt,js)
      fst(jt,js,jy,8)=bx1st*wbxt(jt,js)+bz1st*wbzp(jt,js)
       
!      vrst=fst(jt,js,jy,3)*dcos(thst(jt))+fst(jt,js,jy,5)*dsin(thst(jt))
!      vpst=-fst(jt,js,jy,3)*dsin(thst(jt))+fst(jt,js,jy,5)*dcos(thst(jt))
!      brst=fst(jt,js,jy,6)*dcos(thst(jt))+fst(jt,js,jy,8)*dsin(thst(jt))
!      bpst=-fst(jt,js,jy,6)*dsin(thst(jt))+fst(jt,js,jy,8)*dcos(thst(jt))
!      fst(jt,js,jy,3)=vrst
!      fst(jt,js,jy,5)=vpst
!      fst(jt,js,jy,6)=brst
!      fst(jt,js,jy,8)=bpst
      enddo
      enddo
      enddo    

      do 10 js=mpsa-2,mpsa-4,-1
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),js,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
!ws_smps--------------------------
      call smth_st_nrk(fst,js,8,3)
!ws_smps--------------------
   10 continue

      do 61 m=1,8

      if(m.eq.3 .or. m.eq.6) then !!(m=3,6) ibnd==0
      is=mpsa 
      do 51 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
       fst(jt,is,jy,m)=0.           
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   51 continue

      is=mpsa-1 
      do 52 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
        call interp1d3l(fst(jt,is-3,jy,m),fst(jt,is-2,jy,m),fst(jt,is-1,jy,m),fst(jt,is+1,jy,m), &
                        ps(is-3),ps(is-2),ps(is-1),ps(is+1),ps(is),fst(jt,is,jy,m))                   
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   52 continue

      call smth_st_nrk(fst(:,:,:,m),is,1,3)

      else  !!(m=1,2,4,5,7,8,) bnd==1  
      is=mpsa-1    
      do 31 jy=1,my      
      do jt=itbmin(nrankxz),itbmax(nrankxz)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
!       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is)) 
       fst(jt,is,jy,m)=(4*fst(jt,is-1,jy,m)-fst(jt,is-2,jy,m))/3.0
!       fst(jt,is)=fst(jt,is-1) 
       if(fst(jt,is,jy,m)*fst(jt,is-1,jy,m).le.0) fst(jt,is,jy,m)=0.                  
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   31 continue 
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),is,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
      call smth_st_nrk(fst(:,:,:,m),is,1,3)
!ws_smps--------------------
      is=mpsa
 
      do 41 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
!       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))
       fst(jt,is,jy,m)=(4*fst(jt,is-1,jy,m)-fst(jt,is-2,jy,m))/3.0
!       fst(jt,is)=fst(jt,is-1) 
       if(fst(jt,is,jy,m)*fst(jt,is-1,jy,m).le.0) fst(jt,is,jy,m)=0.                  
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   41 continue
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),is,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
!ws_smps--------------------
      call smth_st_nrk(fst(:,:,:,m),is,1,3)
!	write(*,*) nrankxz,m,"smth done"
!ws_smps--------------------
      endif !!(m=3,6)ibnd==0
   61 continue
      
      do 3 jy=1,my           
      do ikb=1,mb_nrk(nrankxz)
        i=itb_nrk(ikb,nrankxz)
        ib=ib_nrk(ikb,nrankxz) 
        do m=1,8 
        do js=mpsa,mpsa-3,-1     
        call interp1d3l(fst(i-2,js,jy,m),fst(i-1,js,jy,m),fst(i,js,jy,m),fst(i+1,js,jy,m), &
                        thst(i-2),thst(i-1),thst(i),thst(i+1),thb(ib),f1s(ikb,js))
        enddo
        js=mpsa
        call interp1d3l(f1s(ikb,js),f1s(ikb,js-1),f1s(ikb,js-2),f1s(ikb,js-3), &
                        ps(js),ps(js-1),ps(js-2),ps(js-3),psb(ib),fsxz(ikb,m))
        enddo
!        fxz(jbx(ib),jbz(ib))=(asm_b(ib)*x1s(ib,js-1)+bsm_b(ib)*x1s(ib,js-2)+csm_b(ib)*x1s(ib,js-3)) &
!                 /(asm_b(ib)+bsm_b(ib)+csm_b(ib))
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,1)=fsxz(ikb,1)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,2)=fsxz(ikb,2)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,3)=fsxz(ikb,3)*wbrx(ib)+fsxz(ikb,5)*wbpx(ib)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,4)=fsxz(ikb,4)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,5)=fsxz(ikb,3)*wbrz(ib)+fsxz(ikb,5)*wbpz(ib)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,6)=fsxz(ikb,6)*wbrx(ib)+fsxz(ikb,8)*wbpx(ib)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,7)=fsxz(ikb,7)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,8)=fsxz(ikb,6)*wbrz(ib)+fsxz(ikb,8)*wbpz(ib)
      enddo
    3 continue 
      endif !!mb_nrk(nrankxz).gt.0
           
!      if(ibnd==2) then
!      
!      do ika1=1,ma1_nrk(nrankxz)
!        i=ita1_nrk(ika1,nrankxz)
!        ia1=ia1_nrk(ika1,nrankxz)
!        do js=mpsa,mpsa-3,1
!        call interp1d3l(fst(i-2,js),fst(i-1,js),fst(i,js),fst(i+1,js), &
!                        thst(i-2),thst(i-1),thst(i),thst(i+1),tha1(ia1),x1s(ika1,js,jy))
!        enddo
!        js=mpsa
!        call interp1d3l(x1s(ika1,js,jy),x1s(ika1,js-1,jy),x1s(ika1,js-2,jy),x1s(ika1,js-3,jy), &
!                        ps(js),ps(js-1),ps(js-2),ps(js-3),psa1(ia1),fxz(jxa1_nrk(ika1,nrankxz),jza1_nrk(ika1,nrankxz)))
!!        fxz(jxa1(ia1),jza1(ia1))=(asm_a1(ia1)*x1s(ia1,js-1)+bsm_a1(ia1)*x1s(ia1,js-2)+csm_a1(ia1)*x1s(ia1,js-3)) &
!!                  /(asm_a1(ia1)+bsm_a1(ia1)+csm_a1(ia1)) 
!      enddo
!      endif
      return
      end
!****************************************************************
      subroutine valb3_atlastgrid(f3xz)
      USE DECLARE
      implicit real*8 (b-h,o-z)
      real*8 cx1st,cz1st
      integer is
      dimension f3xz(mx,mz,my,3),fst(n2th+5,mps4:mps,my,3),f1s(mbm_nrk,mps4:mps),fsxz(mbm_nrk,3)
      include 'mpif.h'

       integer status(mpi_status_size)

      call map_xz2st(f3xz,fst,3)
         
      if(mb_nrk(nrankxz).gt.0) then      

      do jy=1,my
      do js=mpsa-2,mpsa-4,-1
      do jt=itbmin(nrankxz),itbmax(nrankxz)
      cx1st=fst(jt,js,jy,1)
      cz1st=fst(jt,js,jy,3)
      fst(jt,js,jy,1)=cx1st*wbxr(jt,js)+cz1st*wbzr(jt,js)
      fst(jt,js,jy,3)=cx1st*wbxt(jt,js)+cz1st*wbzp(jt,js)
      enddo
      enddo
      enddo    

      do 10 js=mpsa-2,mpsa-4,-1
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),js,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
!ws_smps--------------------------
      call smth_st_nrk(fst,js,3,3)
!ws_smps--------------------
   10 continue

      do 61 m=1,3

      if(lbndcfix(m)) then !!(m=1) ibnd==0
      is=mpsa 
      do 51 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
       fst(jt,is,jy,m)=0.           
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   51 continue

      is=mpsa-1 
      do 52 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
        call interp1d3l(fst(jt,is-3,jy,m),fst(jt,is-2,jy,m),fst(jt,is-1,jy,m),fst(jt,is+1,jy,m), &
                        ps(is-3),ps(is-2),ps(is-1),ps(is+1),ps(is),fst(jt,is,jy,m))                   
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   52 continue

      call smth_st_nrk(fst(:,:,:,m),is,1,3)

      else  !!(m=2,3) ibnd==1  
      is=mpsa-1    
      do 31 jy=1,my      
      do jt=itbmin(nrankxz),itbmax(nrankxz)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
!       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is)) 
       fst(jt,is,jy,m)=(4*fst(jt,is-1,jy,m)-fst(jt,is-2,jy,m))/3.0
!       fst(jt,is)=fst(jt,is-1)
       if(fst(jt,is,jy,m)*fst(jt,is-1,jy,m).le.0) fst(jt,is,jy,m)=0                   
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   31 continue 
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),is,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
      call smth_st_nrk(fst(:,:,:,m),is,1,3)
!ws_smps--------------------
      is=mpsa
 
      do 41 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
!       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))
       fst(jt,is,jy,m)=(4*fst(jt,is-1,jy,m)-fst(jt,is-2,jy,m))/3.0
!       fst(jt,is)=fst(jt,is-1) 
       if(fst(jt,is,jy,m)*fst(jt,is-1,jy,m).le.0) fst(jt,is,jy,m)=0                  
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   41 continue
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),is,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
!ws_smps--------------------
      call smth_st_nrk(fst(:,:,:,m),is,1,3)
!	write(*,*) nrankxz,m,"smth done"
!ws_smps--------------------
      endif !!(m=1)ibnd==0
   61 continue
      
      do 3 jy=1,my           
      do ikb=1,mb_nrk(nrankxz)
        i=itb_nrk(ikb,nrankxz)
        ib=ib_nrk(ikb,nrankxz) 
        do m=1,3
        do js=mpsa,mpsa-3,-1     
        call interp1d3l(fst(i-2,js,jy,m),fst(i-1,js,jy,m),fst(i,js,jy,m),fst(i+1,js,jy,m), &
                        thst(i-2),thst(i-1),thst(i),thst(i+1),thb(ib),f1s(ikb,js))
        enddo
        js=mpsa
        call interp1d3l(f1s(ikb,js),f1s(ikb,js-1),f1s(ikb,js-2),f1s(ikb,js-3), &
                        ps(js),ps(js-1),ps(js-2),ps(js-3),psb(ib),fsxz(ikb,m))
        enddo
!        fxz(jbx(ib),jbz(ib))=(asm_b(ib)*x1s(ib,js-1)+bsm_b(ib)*x1s(ib,js-2)+csm_b(ib)*x1s(ib,js-3)) &
!                 /(asm_b(ib)+bsm_b(ib)+csm_b(ib))
        f3xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,1)=fsxz(ikb,1)*wbrx(ib)+fsxz(ikb,3)*wbpx(ib)
        f3xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,2)=fsxz(ikb,2)
        f3xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,3)=fsxz(ikb,1)*wbrz(ib)+fsxz(ikb,3)*wbpz(ib)

      enddo
    3 continue 
      endif !!mb_nrk(nrankxz).gt.0
           
!      if(ibnd==2) then
!      
!      do ika1=1,ma1_nrk(nrankxz)
!        i=ita1_nrk(ika1,nrankxz)
!        ia1=ia1_nrk(ika1,nrankxz)
!        do js=mpsa,mpsa-3,1
!        call interp1d3l(fst(i-2,js),fst(i-1,js),fst(i,js),fst(i+1,js), &
!                        thst(i-2),thst(i-1),thst(i),thst(i+1),tha1(ia1),x1s(ika1,js,jy))
!        enddo
!        js=mpsa
!        call interp1d3l(x1s(ika1,js,jy),x1s(ika1,js-1,jy),x1s(ika1,js-2,jy),x1s(ika1,js-3,jy), &
!                        ps(js),ps(js-1),ps(js-2),ps(js-3),psa1(ia1),fxz(jxa1_nrk(ika1,nrankxz),jza1_nrk(ika1,nrankxz)))
!!        fxz(jxa1(ia1),jza1(ia1))=(asm_a1(ia1)*x1s(ia1,js-1)+bsm_a1(ia1)*x1s(ia1,js-2)+csm_a1(ia1)*x1s(ia1,js-3)) &
!!                  /(asm_a1(ia1)+bsm_a1(ia1)+csm_a1(ia1)) 
!      enddo
!      endif
      return
      end
!****************************************************************
      subroutine valb3_atlastgrid_r0p1_v1(f3xz)
      USE DECLARE
      implicit real*8 (b-h,o-z)
      real*8 cx1st,cz1st
      integer is
      dimension f3xz(mx,mz,my,3),fst(n2th+5,mps4:mps,my,3),f1s(mbm_nrk,mps4:mps),fsxz(mbm_nrk,3)
      include 'mpif.h'

       integer status(mpi_status_size)

      call map_xz2st(f3xz,fst,3)
         
      if(mb_nrk(nrankxz).gt.0) then      

      do jy=1,my
      do js=mpsa-2,mpsa-4,-1
      do jt=itbmin(nrankxz),itbmax(nrankxz)
      cx1st=fst(jt,js,jy,1)
      cz1st=fst(jt,js,jy,3)
      fst(jt,js,jy,1)=cx1st*wbxr(jt,js)+cz1st*wbzr(jt,js)
      fst(jt,js,jy,3)=cx1st*wbxt(jt,js)+cz1st*wbzp(jt,js)
      enddo
      enddo
      enddo    

      do 10 js=mpsa-2,mpsa-4,-1
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),js,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
!ws_smps--------------------------
      call smth_st_nrk(fst,js,3,3)
!ws_smps--------------------
   10 continue

      do 61 m=1,3

      if(m.eq.1) then !!(m=1) ibnd==0
      is=mpsa 
      do 51 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
       fst(jt,is,jy,m)=0.           
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   51 continue

      is=mpsa-1 
      do 52 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
        call interp1d3l(fst(jt,is-3,jy,m),fst(jt,is-2,jy,m),fst(jt,is-1,jy,m),fst(jt,is+1,jy,m), &
                        ps(is-3),ps(is-2),ps(is-1),ps(is+1),ps(is),fst(jt,is,jy,m))                   
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   52 continue

      call smth_st_nrk(fst(:,:,:,m),is,1,3)

      else  !!(m=2,3) ibnd==1  
      is=mpsa-1    
      do 31 jy=1,my      
      do jt=itbmin(nrankxz),itbmax(nrankxz)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
!       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is)) 
       fst(jt,is,jy,m)=(4*fst(jt,is-1,jy,m)-fst(jt,is-2,jy,m))/3.0
!       fst(jt,is)=fst(jt,is-1)
       if(fst(jt,is,jy,m)*fst(jt,is-1,jy,m).le.0) fst(jt,is,jy,m)=0                   
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   31 continue 
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),is,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
      call smth_st_nrk(fst(:,:,:,m),is,1,3)
!ws_smps--------------------
      is=mpsa
 
      do 41 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
!       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))
       fst(jt,is,jy,m)=(4*fst(jt,is-1,jy,m)-fst(jt,is-2,jy,m))/3.0
!       fst(jt,is)=fst(jt,is-1) 
       if(fst(jt,is,jy,m)*fst(jt,is-1,jy,m).le.0) fst(jt,is,jy,m)=0                  
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   41 continue
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),is,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
!ws_smps--------------------
      call smth_st_nrk(fst(:,:,:,m),is,1,3)
!	write(*,*) nrankxz,m,"smth done"
!ws_smps--------------------
      endif !!(m=1)ibnd==0
   61 continue
      
      do 3 jy=1,my           
      do ikb=1,mb_nrk(nrankxz)
        i=itb_nrk(ikb,nrankxz)
        ib=ib_nrk(ikb,nrankxz) 
        do m=1,3
        do js=mpsa,mpsa-3,-1     
        call interp1d3l(fst(i-2,js,jy,m),fst(i-1,js,jy,m),fst(i,js,jy,m),fst(i+1,js,jy,m), &
                        thst(i-2),thst(i-1),thst(i),thst(i+1),thb(ib),f1s(ikb,js))
        enddo
        js=mpsa
        call interp1d3l(f1s(ikb,js),f1s(ikb,js-1),f1s(ikb,js-2),f1s(ikb,js-3), &
                        ps(js),ps(js-1),ps(js-2),ps(js-3),psb(ib),fsxz(ikb,m))
        enddo
!        fxz(jbx(ib),jbz(ib))=(asm_b(ib)*x1s(ib,js-1)+bsm_b(ib)*x1s(ib,js-2)+csm_b(ib)*x1s(ib,js-3)) &
!                 /(asm_b(ib)+bsm_b(ib)+csm_b(ib))
        f3xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,1)=fsxz(ikb,1)*wbrx(ib)+fsxz(ikb,3)*wbpx(ib)
        f3xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,2)=fsxz(ikb,2)
        f3xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,3)=fsxz(ikb,1)*wbrz(ib)+fsxz(ikb,3)*wbpz(ib)

      enddo
    3 continue 
      endif !!mb_nrk(nrankxz).gt.0
           
!      if(ibnd==2) then
!      
!      do ika1=1,ma1_nrk(nrankxz)
!        i=ita1_nrk(ika1,nrankxz)
!        ia1=ia1_nrk(ika1,nrankxz)
!        do js=mpsa,mpsa-3,1
!        call interp1d3l(fst(i-2,js),fst(i-1,js),fst(i,js),fst(i+1,js), &
!                        thst(i-2),thst(i-1),thst(i),thst(i+1),tha1(ia1),x1s(ika1,js,jy))
!        enddo
!        js=mpsa
!        call interp1d3l(x1s(ika1,js,jy),x1s(ika1,js-1,jy),x1s(ika1,js-2,jy),x1s(ika1,js-3,jy), &
!                        ps(js),ps(js-1),ps(js-2),ps(js-3),psa1(ia1),fxz(jxa1_nrk(ika1,nrankxz),jza1_nrk(ika1,nrankxz)))
!!        fxz(jxa1(ia1),jza1(ia1))=(asm_a1(ia1)*x1s(ia1,js-1)+bsm_a1(ia1)*x1s(ia1,js-2)+csm_a1(ia1)*x1s(ia1,js-3)) &
!!                  /(asm_a1(ia1)+bsm_a1(ia1)+csm_a1(ia1)) 
!      enddo
!      endif
      return
      end

!****************************************************************
      subroutine valb3_atlastgrid_r1p0_v1(f3xz)
      USE DECLARE
      implicit real*8 (b-h,o-z)
      real*8 cx1st,cz1st
      integer is
      dimension f3xz(mx,mz,my,3),fst(n2th+5,mps4:mps,my,3),f1s(mbm_nrk,mps4:mps),fsxz(mbm_nrk,3)
      include 'mpif.h'

       integer status(mpi_status_size)

      call map_xz2st(f3xz,fst,3)
         
      if(mb_nrk(nrankxz).gt.0) then      

      do jy=1,my
      do js=mpsa-2,mpsa-4,-1
      do jt=itbmin(nrankxz),itbmax(nrankxz)
      cx1st=fst(jt,js,jy,1)
      cz1st=fst(jt,js,jy,3)
      fst(jt,js,jy,1)=cx1st*wbxr(jt,js)+cz1st*wbzr(jt,js)
      fst(jt,js,jy,3)=cx1st*wbxt(jt,js)+cz1st*wbzp(jt,js)
      enddo
      enddo
      enddo    

      do 10 js=mpsa-2,mpsa-4,-1
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),js,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
!ws_smps--------------------------
      call smth_st_nrk(fst,js,3,3)
!ws_smps--------------------
   10 continue

      do 61 m=1,3

      if(m.eq.2 .or. m.eq.3) then !!(m=2,3) ibnd==0
      is=mpsa 
      do 51 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
       fst(jt,is,jy,m)=0.           
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   51 continue

      is=mpsa-1 
      do 52 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
        call interp1d3l(fst(jt,is-3,jy,m),fst(jt,is-2,jy,m),fst(jt,is-1,jy,m),fst(jt,is+1,jy,m), &
                        ps(is-3),ps(is-2),ps(is-1),ps(is+1),ps(is),fst(jt,is,jy,m))                   
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   52 continue

      call smth_st_nrk(fst(:,:,:,m),is,1,3)

      else  !!(m=1) ibnd==1  
      is=mpsa-1    
      do 31 jy=1,my      
      do jt=itbmin(nrankxz),itbmax(nrankxz)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
!       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is)) 
       fst(jt,is,jy,m)=(4*fst(jt,is-1,jy,m)-fst(jt,is-2,jy,m))/3.0 
!       fst(jt,is)=fst(jt,is-1) 
       if(fst(jt,is,jy,m)*fst(jt,is-1,jy,m).le.0) fst(jt,is,jy,m)=0          
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   31 continue 
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),is,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
      call smth_st_nrk(fst(:,:,:,m),is,1,3)
!ws_smps--------------------
      is=mpsa
 
      do 41 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
!       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))
       fst(jt,is,jy,m)=(4*fst(jt,is-1,jy,m)-fst(jt,is-2,jy,m))/3.0
!       fst(jt,is)=fst(jt,is-1)   
       if(fst(jt,is,jy,m)*fst(jt,is-1,jy,m).le.0) fst(jt,is,jy,m)=0     
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   41 continue
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),is,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
!ws_smps--------------------
      call smth_st_nrk(fst(:,:,:,m),is,1,3)
!	write(*,*) nrankxz,m,"smth done"
!ws_smps--------------------
      endif !!(m=2,3)ibnd==0
   61 continue
      
      do 3 jy=1,my           
      do ikb=1,mb_nrk(nrankxz)
        i=itb_nrk(ikb,nrankxz)
        ib=ib_nrk(ikb,nrankxz) 
        do m=1,3
        do js=mpsa,mpsa-3,-1     
        call interp1d3l(fst(i-2,js,jy,m),fst(i-1,js,jy,m),fst(i,js,jy,m),fst(i+1,js,jy,m), &
                        thst(i-2),thst(i-1),thst(i),thst(i+1),thb(ib),f1s(ikb,js))
        enddo
        js=mpsa
        call interp1d3l(f1s(ikb,js),f1s(ikb,js-1),f1s(ikb,js-2),f1s(ikb,js-3), &
                        ps(js),ps(js-1),ps(js-2),ps(js-3),psb(ib),fsxz(ikb,m))
        enddo
!        fxz(jbx(ib),jbz(ib))=(asm_b(ib)*x1s(ib,js-1)+bsm_b(ib)*x1s(ib,js-2)+csm_b(ib)*x1s(ib,js-3)) &
!                 /(asm_b(ib)+bsm_b(ib)+csm_b(ib))
        f3xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,1)=fsxz(ikb,1)*wbrx(ib)+fsxz(ikb,3)*wbpx(ib)
        f3xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,2)=fsxz(ikb,2)
        f3xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,3)=fsxz(ikb,1)*wbrz(ib)+fsxz(ikb,3)*wbpz(ib)

      enddo
    3 continue 
      endif !!mb_nrk(nrankxz).gt.0
           
!      if(ibnd==2) then
!      
!      do ika1=1,ma1_nrk(nrankxz)
!        i=ita1_nrk(ika1,nrankxz)
!        ia1=ia1_nrk(ika1,nrankxz)
!        do js=mpsa,mpsa-3,1
!        call interp1d3l(fst(i-2,js),fst(i-1,js),fst(i,js),fst(i+1,js), &
!                        thst(i-2),thst(i-1),thst(i),thst(i+1),tha1(ia1),x1s(ika1,js,jy))
!        enddo
!        js=mpsa
!        call interp1d3l(x1s(ika1,js,jy),x1s(ika1,js-1,jy),x1s(ika1,js-2,jy),x1s(ika1,js-3,jy), &
!                        ps(js),ps(js-1),ps(js-2),ps(js-3),psa1(ia1),fxz(jxa1_nrk(ika1,nrankxz),jza1_nrk(ika1,nrankxz)))
!!        fxz(jxa1(ia1),jza1(ia1))=(asm_a1(ia1)*x1s(ia1,js-1)+bsm_a1(ia1)*x1s(ia1,js-2)+csm_a1(ia1)*x1s(ia1,js-3)) &
!!                  /(asm_a1(ia1)+bsm_a1(ia1)+csm_a1(ia1)) 
!      enddo
!      endif
      return
      end
!ws************************************************************
     subroutine bndry8_x_ex(ibnd)
      USE DECLARE
      integer ibnd
      include 'mpif.h'

      do 2 jy=1,my
      x1(:,:,jy,:)=x(:,:,jy,:)-xint(:,:,:)
   2  continue
      
      call bndry8_ex(x1,ibnd)        
      if(smoothp1ll) call smthp1_traceline(3)  
      do 11 jy=1,my
      x(:,:,jy,:)=x1(:,:,jy,:)+xint(:,:,:)
   11 continue
  
      return
      end
!ws***************************************
      subroutine bndry8_ex(f8xz,ibnd)
      USE DECLARE
      integer ibnd
      real*8,dimension(mx,mz,my,8) :: f8xz
      include 'mpif.h'

      select case(ibnd)
      case(0)
      call valbm_atlastgrid_v1(f8xz(:,:,:,:),8,0)
      case(1)
      call valbm_atlastgrid_v1(f8xz(:,:,:,:),8,1)         
      case(3)
      call valbm_atlastgrid_v1(f8xz(:,:,:,1:5),5,1)
      call valbm_atlastgrid_v1(f8xz(:,:,:,6:8),3,0)
      case(10)
      call valb8_atlastgrid_r0p1_v1(f8xz)
      case(20)
      call valb8_atlastgrid_r0p1_v2(f8xz)
!      call x1_atlastgrid_r0p1_v1(3)
      case default
      call valb8_atlastgrid(f8xz)
      end select  

      call mpi_transfersm(f8xz(:,:,:,:),8) 
      if(smoothx1) then 
      do m=1,8    
       call smthxzy(f8xz(:,:,:,m),1)
      enddo
      endif

      return
      end       

!ws***************************************
      subroutine bndry3_ex(f3xz,ibnd)
      USE DECLARE
      integer ibnd
      real*8,dimension(mx,mz,my,3) :: f3xz
      include 'mpif.h'

      select case(ibnd)
      case(0)
      call valbm_atlastgrid_v1(f3xz(:,:,:,:),3,0)
      case(1)
      call valbm_atlastgrid_v1(f3xz(:,:,:,:),3,1)  
             
      case(10)
      call valb3_atlastgrid_r0p1_v1(f3xz)
!      call cur_atlastgrid_r0p1_v1(3) 
      case(11)
      call valb3_atlastgrid_r1p0_v1(f3xz)
      case default
      call valb3_atlastgrid(f3xz)
      end select  
        
      call mpi_transfersm(f3xz(:,:,:,:),3) 

      if(smoothc) then 
      do m=1,3    
       call smthxzy(f3xz(:,:,:,m),1)
      enddo
      endif

      return
      end 
      
!ws******************************************************************************
!ws:artificial soundwave 
!ws******************************************************************************
      subroutine stepon_atfs
!
!     This routine time-advances X's bz fourth order in time and second
!     order in space Runge-Kotta differential scheme.
!     note: X is alwazs the up-to-date value while Xm being the
!           intermediate value, and Xdif is increment
!
!
      USE DECLARE
      include 'mpif.h'
!
!      nst=0
!      call recrd_dbg
      dts=dt/ncycl_atfs
      ms=1
      me=8
      ml=1
      tt=time
      tt1=time+dt/6.
      tt2=time
      irk=1
      call right
!      call artif_sound(dt/2)
       xfold(:,:,:,:)=x(:,:,:,:)
       xm(:,:,:,:)=xfold(:,:,:,:)+xdif(:,:,:,:)*dt/6.
       x(:,:,:,:)=xfold(:,:,:,:)+xdif(:,:,:,:)*dt/2.
!       call artif_sound_replace(2)
!       xm(:,:,:,2)=xm(:,:,:,2)+p_atfs(:,:,:)/3.
!       x(:,:,:,2)=x(:,:,:,2)+p_atfs(:,:,:)
!      nst=1
!      call recrd_dbg
!
      tt=time+dt/2.
      tt1=time+dt/2.
      tt2=time+dt/6.
      irk=2
        call right
!        call artif_sound(dt/2)
        xm(:,:,:,:)=xm(:,:,:,:)+xdif(:,:,:,:)*dt/3.
        x(:,:,:,:)=xfold(:,:,:,:)+xdif(:,:,:,:)*dt/2.
!        call artif_sound_replace(2)
!        xm(:,:,:,2)=xm(:,:,:,2)+p_atfs(:,:,:)*2./3.
!        x(:,:,:,2)=x(:,:,:,2)+p_atfs(:,:,:)

!      nst=2
!      call recrd_dbg
!
      tt1=time+5.*dt/6.
      tt2=time+dt/2.
      irk=3
        call right
!       call artif_sound(dt)
        xm(:,:,:,:)=xm(:,:,:,:)+xdif(:,:,:,:)*dt/3.
        x(:,:,:,:)=xfold(:,:,:,:)+xdif(:,:,:,:)*dt
!       call artif_sound_replace(1)
!        xm(:,:,:,2)=xm(:,:,:,2)+p_atfs(:,:,:)/3.
!        x(:,:,:,2)=x(:,:,:,2)+p_atfs(:,:,:)

!      nst=3
!      call recrd_dbg
!
      time=time+dt
      tt1=time+dt
      tt2=time+5.*dt/6.
      irk=4
        call right
!        call artif_sound(dt/6)
        x(:,:,:,:)=xm(:,:,:,:)+xdif(:,:,:,:)*dt/6.
!        call artif_sound_replace(6)
!        x(:,:,:,2)=x(:,:,:,2)+p_atfs(:,:,:)
!        call artif_sound(1)
!      call artif_sound_replace(1)
!      call artif_sound_replace_Lax(1)
      call artif_sound_replace_RK(1)
!      call artif_sound_implicity(1)

      caf=0.75d0*(0.5+0.5*dtanh((time-40)/5.))
!      call bndry8_x_ex(lbnd)

      return
      end
!ws******************************************************************************
      subroutine stepon_atfs_RK
!
!     This routine time-advances X's bz fourth order in time and second
!     order in space Runge-Kotta differential scheme.
!     note: X is alwazs the up-to-date value while Xm being the
!           intermediate value, and Xdif is increment
!
!
      USE DECLARE
      include 'mpif.h'
!
!      nst=0
!      call recrd_dbg
      dts=dt/ncycl_atfs
      ms=1
      me=8
      ml=1
      tt=time
      tt1=time+dt/6.
      tt2=time
      irk=1
      call right
!      call artif_sound(dt/2)
       xfold(:,:,:,:)=x(:,:,:,:)
       xm(:,:,:,:)=xfold(:,:,:,:)+xdif(:,:,:,:)*dt/6.
       x(:,:,:,:)=xfold(:,:,:,:)+xdif(:,:,:,:)*dt/2.
       call artif_sound_replace(2)
!       xm(:,:,:,2)=xm(:,:,:,2)+p_atfs(:,:,:)/3.
!       x(:,:,:,2)=x(:,:,:,2)+p_atfs(:,:,:)
!      nst=1
!      call recrd_dbg
!
      tt=time+dt/2.
      tt1=time+dt/2.
      tt2=time+dt/6.
      irk=2
        call right
!        call artif_sound(dt/2)
        xm(:,:,:,:)=xm(:,:,:,:)+xdif(:,:,:,:)*dt/3.
        x(:,:,:,:)=xfold(:,:,:,:)+xdif(:,:,:,:)*dt/2.
        call artif_sound_replace(2)
!        xm(:,:,:,2)=xm(:,:,:,2)+p_atfs(:,:,:)*2./3.
!        x(:,:,:,2)=x(:,:,:,2)+p_atfs(:,:,:)

!      nst=2
!      call recrd_dbg
!
      tt1=time+5.*dt/6.
      tt2=time+dt/2.
      irk=3
        call right
!       call artif_sound(dt)
        xm(:,:,:,:)=xm(:,:,:,:)+xdif(:,:,:,:)*dt/3.
        x(:,:,:,:)=xfold(:,:,:,:)+xdif(:,:,:,:)*dt
        call artif_sound_replace(1)
!        xm(:,:,:,2)=xm(:,:,:,2)+p_atfs(:,:,:)/3.
!        x(:,:,:,2)=x(:,:,:,2)+p_atfs(:,:,:)

!      nst=3
!      call recrd_dbg
!
      time=time+dt
      tt1=time+dt
      tt2=time+5.*dt/6.
      irk=4
        call right
!        call artif_sound(dt/6)
        x(:,:,:,:)=xm(:,:,:,:)+xdif(:,:,:,:)*dt/6.
        call artif_sound_replace(6)
!        x(:,:,:,2)=x(:,:,:,2)+p_atfs(:,:,:)
!        call artif_sound(1)

      caf=0.75d0*(0.5+0.5*dtanh((time-40)/5.))
      call bndry8_x_ex(lbnd)

      return
      end
!ws**************************************************************************
     subroutine artif_sound(n)
      USE DECLARE
      include 'mpif.h'
      real*8, dimension(mx,mz,my,2) :: u,um !, p_atfs
!      real*8 dts
      integer n,ncyc,nc
      d1f2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(fp1-f0) &
         -(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)

      d2f2(fm1,f0,fp1,xm1,x0,xp1)= &
       2*((fp1-f0)/(xp1-x0)-(f0-fm1)/(x0-xm1))/(xp1-xm1)

      dts=dt/ncycl_atfs
      ncyc=ncycl_atfs/n
      do jy=1,my
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
      u(jx,jz,jy,1)=(x(jx,jz,jy,2)-xint(jx,jz,2))/x(jx,jz,jy,1)
      u(jx,jz,jy,2)=0
      enddo
      enddo
      enddo

      udx=0
      udz=0
      udy=0
      udx2=0
      udz2=0
      udy2=0

      do nc=1,ncyc
      um(:,:,:,:)=u(:,:,:,:)
      do jy=iy_first+1,iy_last-1
      do jz=iz_first+1,iz_last-1
      do jx=ix_first+1,ix_last-1
      tmdx=d1f2(um(jx-1,jz,jy,1),um(jx,jz,jy,1),um(jx+1,jz,jy,1),xx(jx-1),xx(jx),xx(jx+1))
      tmdz=d1f2(um(jx,jz-1,jy,1),um(jx,jz,jy,1),um(jx,jz+1,jy,1),zz(jz-1),zz(jz),zz(jz+1))
      tmdy=d1f2(um(jx,jz,jy-1,1),um(jx,jz,jy,1),um(jx,jz,jy+1,1),yy(jy-1),yy(jy),yy(jy+1))
       

      u(jx,jz,jy,2)=um(jx,jz,jy,2)+dts*cs_atf*( &
                    x(jx,jz,jy,6)*tmdx &
                   +x(jx,jz,jy,8)*tmdz &
                   +x(jx,jz,jy,7)*tmdy/xx(jx)) &
                   +fmu_atf*(udx2+udx/xx(jx)+udy2/xx(jx)**2+udz2)

      udx=d1f2(um(jx-1,jz,jy,2),um(jx,jz,jy,2),um(jx+1,jz,jy,2),xx(jx-1),xx(jx),xx(jx+1))
      udz=d1f2(um(jx,jz-1,jy,2),um(jx,jz,jy,2),um(jx,jz+1,jy,2),zz(jz-1),zz(jz),zz(jz+1))
      udy=d1f2(um(jx,jz,jy-1,2),um(jx,jz,jy,2),um(jx,jz,jy+1,2),yy(jy-1),yy(jy),yy(jy+1))
      udx2=d2f2(um(jx-1,jz,jy,2),um(jx,jz,jy,2),um(jx+1,jz,jy,2),xx(jx-1),xx(jx),xx(jx+1))
      udz2=d2f2(um(jx,jz-1,jy,2),um(jx,jz,jy,2),um(jx,jz+1,jy,2),zz(jz-1),zz(jz),zz(jz+1))
      udy2=d2f2(um(jx,jz,jy-1,2),um(jx,jz,jy,2),um(jx,jz,jy+1,2),yy(jy-1),yy(jy),yy(jy+1))

      u(jx,jz,jy,1)=(um(jx+1,jz,jy,1)+um(jx-1,jz,jy,1)+um(jx,jz+1,jy,1)+um(jx,jz-1,jy,1)+um(jx,jz,jy+1,1)+um(jx,jz,jy-1,1))/6. &
                   +dts*cs_atf/x(jx,jz,jy,1)*( &
                    x(jx,jz,jy,6)*udx &
                   +x(jx,jz,jy,8)*udz &
                   +x(jx,jz,jy,7)*udy/xx(jx) )
      enddo
      enddo
      enddo
      call mpi_transfersm(u,2)
      enddo



      do jy=1,my
      do jx=ix_first,ix_last
      do jz=iz_first,iz_last
      x(jx,jz,jy,2)=u(jx,jz,jy,1)*x(jx,jz,jy,1)+xint(jx,jz,2)
      enddo
      enddo
      enddo
!      call mpi_transfer1(p_atfs)
      return
      end
!ws***********************************************************
     subroutine artif_sound_replace(n)
      USE DECLARE
      include 'mpif.h'
      real*8, dimension(mx,mz,my,2) :: u,um,udif
!      real*8 dts
      real*8, dimension(mx,mz,my,8) :: x_asw,dx_asw
!      real*8 dts
      integer n,ncyc,nc

!      dts=dt/ncycl_atfs
      ncyc=ncycl_atfs/n

      x_asw(:,:,:,:)=xfold(:,:,:,:)
      dx_asw(:,:,:,:)=(x(:,:,:,:)-xfold(:,:,:,:))/ncyc
      do jy=1,my
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
      u(jx,jz,jy,1)=(xfold(jx,jz,jy,2)-xint(jx,jz,2))/xfold(jx,jz,jy,1)
      dx_asw(jx,jz,jy,2)=(x(jx,jz,jy,2)/x(jx,jz,jy,1)-xfold(jx,jz,jy,2)/xfold(jx,jz,jy,1))/ncyc
      u(jx,jz,jy,2)=0

      enddo
      enddo
      enddo
            
      do nc=1,ncyc
      x_asw(:,:,:,:)=x_asw(:,:,:,:)+dx_asw(:,:,:,:)        
      u(:,:,:,1)=u(:,:,:,1)+dx_asw(:,:,:,2)
      call right_atfs(x_asw,u,udif)
      u(:,:,:,:)=u(:,:,:,:)+udif(:,:,:,:)*dts 
      call mpi_transfersm(u,2)
      enddo

      do jy=1,my
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
      x(jx,jz,jy,2)=u(jx,jz,jy,1)*x(jx,jz,jy,1)+xint(jx,jz,2)
      enddo
      enddo
      enddo

!      if(nrankxz==nrank_mode) then
!      open(unit=501,file='u1.dat',status='unknown',form='formatted')
!      write(501,100)((u(jx,jz,1,1),jx=1,mx),jz=1,mz)
!      close(501)
!      open(unit=502,file='up.dat',status='unknown',form='formatted')
!      write(502,100)((x(jx,jz,1,2),jx=1,mx),jz=1,mz)
!      close(502)
!100   format(1(1x,e12.5))  
!      endif
!      call mpi_transfer1(p_atfs)
      return
      end
!ws***********************************************************
     subroutine artif_sound_replace_Lax(n)
      USE DECLARE
      include 'mpif.h'
      real*8, dimension(mx,mz,my,2) :: u,um,udif
      real*8, dimension(mx,mz,my,8) :: x_asw,dx_asw
!      real*8 dts
      integer n,ncyc,nc

!      dts=dt/ncycl_atfs
      ncyc=ncycl_atfs/n

      x_asw(:,:,:,:)=xfold(:,:,:,:)
      dx_asw(:,:,:,:)=(x(:,:,:,:)-xfold(:,:,:,:))/ncyc
      do jy=1,my
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
      u(jx,jz,jy,1)=(xfold(jx,jz,jy,2)-xint(jx,jz,2))/xfold(jx,jz,jy,1)
      dx_asw(jx,jz,jy,2)=(x(jx,jz,jy,2)/x(jx,jz,jy,1)-xfold(jx,jz,jy,2)/xfold(jx,jz,jy,1))/ncyc
      u(jx,jz,jy,2)=0

      enddo
      enddo
      enddo
            
      do nc=1,ncyc 
      x_asw(:,:,:,:)=x_asw(:,:,:,:)+dx_asw(:,:,:,:)      
      do jy=iy_first+1,iy_last-1
      do jz=iz_first+1,iz_last-1
      do jx=ix_first+1,ix_last-1
      um(jx,jz,jy,1)=(u(jx+1,jz,jy,1)+u(jx-1,jz,jy,1)+u(jx,jz+1,jy,1)+u(jx,jz-1,jy,1)+u(jx,jz,jy+1,1)+u(jx,jz,jy-1,1))/6 &
                    +dx_asw(jx,jz,jy,2)
      um(jx,jz,jy,2)=u(jx,jz,jy,2)
      enddo 
      enddo
      enddo
    
      call right_atfs(x_asw,u,udif)
!ws:wrong      u(:,:,:,:)=u(:,:,:,:)+udif(:,:,:,:)*dts 
      u(:,:,:,:)=um(:,:,:,:)+udif(:,:,:,:)*dts
      call mpi_transfersm(u,2)
      enddo

      do jy=1,my
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
      x(jx,jz,jy,2)=u(jx,jz,jy,1)*x(jx,jz,jy,1)+xint(jx,jz,2)
      enddo
      enddo
      enddo

!      if(nrank==nrank_mode) then
!      open(unit=501,file='u1.dat',status='unknown',form='formatted')
!      write(501,100)((u(jx,jz,1,1),jx=1,mx),jz=1,mz)
!      close(501)
!      open(unit=502,file='up.dat',status='unknown',form='formatted')
!      write(502,100)((x(jx,jz,1,2),jx=1,mx),jz=1,mz)
!      close(502)
!100   format(1(1x,e12.5))  
!      endif
!      call mpi_transfer1(p_atfs)
      return
      end
!ws***********************************************************
     subroutine artif_sound_replace_RK(n)
      USE DECLARE
      include 'mpif.h'
      real*8, dimension(mx,mz,my,2) :: u,udif,ufold,um
      real*8, dimension(mx,mz,my,8) :: x_asw,dx_asw
!      real*8 dts
      integer n,ncyc,nc

!      dts=dt/ncycl_atfs
      ncyc=ncycl_atfs/n

      x_asw(:,:,:,:)=xfold(:,:,:,:)
      dx_asw(:,:,:,:)=(x(:,:,:,:)-xfold(:,:,:,:))/ncyc
      do jy=1,my
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
      u(jx,jz,jy,1)=(xfold(jx,jz,jy,2)-xint(jx,jz,2))/xfold(jx,jz,jy,1)
      dx_asw(jx,jz,jy,2)=(x(jx,jz,jy,2)/x(jx,jz,jy,1)-xfold(jx,jz,jy,2)/xfold(jx,jz,jy,1))/ncyc
      u(jx,jz,jy,2)=0

      enddo
      enddo
      enddo
            
      do nc=1,ncyc      
        x_asw(:,:,:,:)=x_asw(:,:,:,:)+dx_asw(:,:,:,:)
        ufold(:,:,:,1)=u(:,:,:,1)+dx_asw(:,:,:,2)
        ufold(:,:,:,2)=u(:,:,:,2)
      call right_atfs(x_asw,u,udif)
        um(:,:,:,:)=ufold(:,:,:,:)+udif(:,:,:,:)*dts/6.
        u(:,:,:,:)=ufold(:,:,:,:)+udif(:,:,:,:)*dts/2.
      call mpi_transfersm(u,2)
      call right_atfs(x_asw,u,udif)
        um(:,:,:,:)=um(:,:,:,:)+udif(:,:,:,:)*dts/3.
        u(:,:,:,:)=ufold(:,:,:,:)+udif(:,:,:,:)*dts/2.
      call mpi_transfersm(u,2)
      call right_atfs(x_asw,u,udif)
        um(:,:,:,:)=um(:,:,:,:)+udif(:,:,:,:)*dts/3.
        u(:,:,:,:)=ufold(:,:,:,:)+udif(:,:,:,:)*dts
      call mpi_transfersm(u,2)      
      call right_atfs(x_asw,u,udif)
        u(:,:,:,:)=um(:,:,:,:)+udif(:,:,:,:)*dts/6.
      call mpi_transfersm(u,2)
      enddo 


      do jy=1,my
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
      x(jx,jz,jy,2)=u(jx,jz,jy,1)*x(jx,jz,jy,1)+xint(jx,jz,2)
      enddo
      enddo
      enddo

!      if(nrank==nrank_mode) then
!      open(unit=501,file='u1.dat',status='unknown',form='formatted')
!      write(501,100)((u(jx,jz,1,1),jx=1,mx),jz=1,mz)
!      close(501)
!      open(unit=502,file='up.dat',status='unknown',form='formatted')
!      write(502,100)((x(jx,jz,1,2),jx=1,mx),jz=1,mz)
!      close(502)
!100   format(1(1x,e12.5))  
!      endif
!      call mpi_transfer1(p_atfs)
      return
      end

!ws***********************************************************
     subroutine right_atfs(x_as,u,udif)
      USE DECLARE
      include 'mpif.h'
      real*8, dimension(mx,mz,my,2) :: u,udif
      real*8, dimension(mx,mz,my,8) :: x_as
      d1f2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(fp1-f0) &
         -(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)

      d2f2(fm1,f0,fp1,xm1,x0,xp1)= &
       2*((fp1-f0)/(xp1-x0)-(f0-fm1)/(x0-xm1))/(xp1-xm1)
      do jy=iy_first+1,iy_last-1
      do jz=iz_first+1,iz_last-1
      do jx=ix_first+1,ix_last-1
      if(psi(jx,jz).lt.psia1) then
      tmdx=d1f2(u(jx-1,jz,jy,1),u(jx,jz,jy,1),u(jx+1,jz,jy,1),xx(jx-1),xx(jx),xx(jx+1))
      tmdz=d1f2(u(jx,jz-1,jy,1),u(jx,jz,jy,1),u(jx,jz+1,jy,1),zz(jz-1),zz(jz),zz(jz+1))
      tmdy=d1f2(u(jx,jz,jy-1,1),u(jx,jz,jy,1),u(jx,jz,jy+1,1),yy(jy-1),yy(jy),yy(jy+1))
      udx= d1f2(u(jx-1,jz,jy,2),u(jx,jz,jy,2),u(jx+1,jz,jy,2),xx(jx-1),xx(jx),xx(jx+1))
      udz= d1f2(u(jx,jz-1,jy,2),u(jx,jz,jy,2),u(jx,jz+1,jy,2),zz(jz-1),zz(jz),zz(jz+1))
      udy= d1f2(u(jx,jz,jy-1,2),u(jx,jz,jy,2),u(jx,jz,jy+1,2),yy(jy-1),yy(jy),yy(jy+1))
      udx2=d2f2(u(jx-1,jz,jy,2),u(jx,jz,jy,2),u(jx+1,jz,jy,2),xx(jx-1),xx(jx),xx(jx+1))
      udz2=d2f2(u(jx,jz-1,jy,2),u(jx,jz,jy,2),u(jx,jz+1,jy,2),zz(jz-1),zz(jz),zz(jz+1))
      udy2=d2f2(u(jx,jz,jy-1,2),u(jx,jz,jy,2),u(jx,jz,jy+1,2),yy(jy-1),yy(jy),yy(jy+1)) 

      udif(jx,jz,jy,1)=cs_atf/x_as(jx,jz,jy,1)*( &
                    x_as(jx,jz,jy,6)*udx &
                   +x_as(jx,jz,jy,8)*udz &
                   +x_as(jx,jz,jy,7)*udy/xx(jx) )
      udif(jx,jz,jy,2)=cs_atf*( &
                    x_as(jx,jz,jy,6)*tmdx &
                   +x_as(jx,jz,jy,8)*tmdz &
                   +x_as(jx,jz,jy,7)*tmdy/xx(jx)) &
                   +fmu_atf*(udx2+udx/xx(jx)+udy2/xx(jx)**2+udz2)
      endif

      enddo
      enddo
      enddo

!      call mpi_transfersm(udif,2)

      return
      end

!ws***********************************************************
     subroutine artif_sound_implicity(n)
      USE DECLARE
      include 'mpif.h'
      real*8, dimension(mx,mz,my,2) :: u,udif,ufold,um
      real*8, dimension(mx,mz,my) :: aay,bby,ccy,css,caa
      real*8, dimension(my) :: ssy,uuy
!      real*8 dts
      integer n,ncyc,nc

!      dts=dt/ncycl_atfs
      ncyc=ncycl_atfs/n
      
      cdts=dts*cs_atf/dyy
      
      do jy=iy_first+1,iy_last-1
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
      u(jx,jz,jy,1)=(x(jx,jz,jy,2)-xint(jx,jz,2))/x(jx,jz,jy,1)
!      u(jx,jz,jy,1)=x(jx,jz,jy,2)/x(jx,jz,jy,1)
      u(jx,jz,jy,2)=0

      cdx=cdts/xx(jx)
      css(jx,jz,jy)=cdx*x(jx,jz,jy,7)/x(jx,jz,jy,1)
      caa(jx,jz,jy)=cdx*x(jx,jz,jy,7)
      cbb=cdx*(x(jx,jz,jy,7)+x(jx,jz,jy+1,7))/2.
      ccc=cdx*(x(jx,jz,jy,7)+x(jx,jz,jy-1,7))/2.
      
      aay(jx,jz,jy)=1+css(jx,jz,jy)*caa(jx,jz,jy)*2.
      bby(jx,jz,jy)=-css(jx,jz,jy)*cbb
      ccy(jx,jz,jy)=-css(jx,jz,jy)*ccc
      enddo
      enddo
      enddo
          
            
      do nc=1,ncyc      
      call right_xz_atfs(u,udif)
      um(:,:,:,:)=u+udif(:,:,:,:)*dts

      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
      if(psi(jx,jz).lt.psia1) then
      do jy=iy_first+1,iy_last-1
      ssy(jy)=um(jx,jz,jy,1)+css(jx,jz,jy)*(um(jx,jz,jy+1,2)-um(jx,jz,jy-1,2))/2.
      enddo
      call tridag_real_period(aay(jx,jz,:),bby(jx,jz,:),ccy(jx,jz,:),ssy,uuy,my) 
      do jy=iy_first+1,iy_last-1     
      u(jx,jz,jy,1)=uuy(jy)
      u(jx,jz,jy,2)=um(jx,jz,jy,2)+caa(jx,jz,jy)*(uuy(jy+1)-uuy(jy-1))/2
      enddo
      endif
      enddo
      enddo
!WS131117
!      call valbm_atlastgrid_v1(u(:,:,:,:),2,1)
      call mpi_transfersm(u,2)
      enddo 

       call smthe4(u(:,:,:,1),1)

      do jy=1,my
      do jx=ix_first,ix_last
      do jz=iz_first,iz_last
      if(psi(jx,jz).lt.psia1) then
      x(jx,jz,jy,2)=u(jx,jz,jy,1)*x(jx,jz,jy,1)+xint(jx,jz,2)
!      x(jx,jz,jy,2)=u(jx,jz,jy,1)*x(jx,jz,jy,1) !+xint(jx,jz,2)
      endif
      enddo
      enddo
      enddo
!      call mpi_transfer1(p_atfs)
      return
      end
 !ws***********************************************************
     subroutine right_xz_atfs(u,udif)
      USE DECLARE
      include 'mpif.h'
      real*8, dimension(mx,mz,my,2) :: u,udif
      d1f2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(fp1-f0) &
         -(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)

      d2f2(fm1,f0,fp1,xm1,x0,xp1)= &
       2*((fp1-f0)/(xp1-x0)-(f0-fm1)/(x0-xm1))/(xp1-xm1)
      do jy=iy_first+1,iy_last-1
      do jz=iz_first+1,iz_last-1
      do jx=ix_first+1,ix_last-1
      tmdx=d1f2(u(jx-1,jz,jy,1),u(jx,jz,jy,1),u(jx+1,jz,jy,1),xx(jx-1),xx(jx),xx(jx+1))
      tmdz=d1f2(u(jx,jz-1,jy,1),u(jx,jz,jy,1),u(jx,jz+1,jy,1),zz(jz-1),zz(jz),zz(jz+1))
      tmdy=d1f2(u(jx,jz,jy-1,1),u(jx,jz,jy,1),u(jx,jz,jy+1,1),yy(jy-1),yy(jy),yy(jy+1))
      udx= d1f2(u(jx-1,jz,jy,2),u(jx,jz,jy,2),u(jx+1,jz,jy,2),xx(jx-1),xx(jx),xx(jx+1))
      udz= d1f2(u(jx,jz-1,jy,2),u(jx,jz,jy,2),u(jx,jz+1,jy,2),zz(jz-1),zz(jz),zz(jz+1))
      udy= d1f2(u(jx,jz,jy-1,2),u(jx,jz,jy,2),u(jx,jz,jy+1,2),yy(jy-1),yy(jy),yy(jy+1))
      udx2=d2f2(u(jx-1,jz,jy,2),u(jx,jz,jy,2),u(jx+1,jz,jy,2),xx(jx-1),xx(jx),xx(jx+1))
      udz2=d2f2(u(jx,jz-1,jy,2),u(jx,jz,jy,2),u(jx,jz+1,jy,2),zz(jz-1),zz(jz),zz(jz+1))
      udy2=d2f2(u(jx,jz,jy-1,2),u(jx,jz,jy,2),u(jx,jz,jy+1,2),yy(jy-1),yy(jy),yy(jy+1)) 

      udif(jx,jz,jy,1)=cs_atf/x(jx,jz,jy,1)*( &
                    x(jx,jz,jy,6)*udx &
!                   +x(jx,jz,jy,7)*udy/xx(jx) &  
                   +x(jx,jz,jy,8)*udz  )
      udif(jx,jz,jy,2)=cs_atf*( &
                    x(jx,jz,jy,6)*tmdx &
!                   +x(jx,jz,jy,7)*tmdy/xx(jx) &
                   +x(jx,jz,jy,8)*tmdz ) &
                   +fmu_atf*(udx2+udx/xx(jx)+udy2/xx(jx)**2+udz2)

      enddo
      enddo
      enddo

!      call mpi_transfersm(udif,2)

      return
      end
!ws:******************************
      SUBROUTINE tridag_real_period(a,b,c,s,u,n)
      include 'mpif.h'
!      complex*16 r(n),u(n)
      INTEGER n,NMAX
      REAL*8 a(n),b(n),c(n),s(n),u(n)
      REAL*8 a1(n),b1(n),c1(n),s1(n),vu(n),qu(n),pu(n)
!      PARAMETER (NMAX=5001)
      INTEGER j
      real*8 bet  !,gam(NMAX)

      if(a(1).eq.0.)pause 'tridag: rewrite equations'
      qu(1)=b(1)/a(1)
      pu(1)=c(1)/a(1)
      vu(1)=s(1)/a(1)
      do j=2,n-1
      bet=a(j)-c(j)*qu(j-1)
      if(bet.eq.0.)pause 'tridag failed'
      qu(j)=b(j)/bet
      pu(j)=-c(j)*pu(j-1)/bet
      vu(j)=(s(j)-c(j)*vu(j-1))/bet
      enddo

      b1(1)=b(n)
      a1(1)=a(n)
      s1(1)=s(n)
      do j=2,n-1
      b1(j)=-b1(j-1)*qu(j-1)
      a1(j)=a1(j-1)-b1(j-1)*pu(j-1)
      s1(j)=s1(j-1)-b1(j-1)*vu(j-1)
      enddo
      
      u(n)=(s1(n-1)-(b1(n-1)+c(n))*vu(n-1))/(a1(n-1)-(b1(n-1)+c(n))*(qu(n-1)+pu(n-1)))
      do j=n-1,1,-1
      u(j)=vu(j)-qu(j)*u(j+1)-pu(j)*u(n)
      enddo

      return
      END  

!ws*******************************************************
    subroutine init_dgn
     USE DECLARE
     include 'mpif.h' 
     integer itmp,ii,ikdgn
     character*10 output
     character*3 cn1
      
!      qdgn(jdgn)=jdgn

      j=1 
      do jdgn=1,mdgn_rs         
      do while(q_NOVA(j) .ge. qdgn(jdgn))
      j=j+1
      enddo
      if(j==1) j=j+2
      if(j==2) j=j+1
      call interp1d3l(psival_NOVA(j-2),psival_NOVA(j-1),psival_NOVA(j),psival_NOVA(j+1), &
                    q_NOVA(j-2),q_NOVA(j-1),q_NOVA(j),q_NOVA(j+1),qdgn(jdgn),psdgn(jdgn))
       do i=1,n2th+5
         call interp1d3l(xxst(i,j-2),xxst(i,j-1),xxst(i,j),xxst(i,j+1), &
                    q_NOVA(j-2),q_NOVA(j-1),q_NOVA(j),q_NOVA(j+1), qdgn(jdgn),xxdgn(i,jdgn))
         call interp1d3l(zzst(i,j-2),zzst(i,j-1),zzst(i,j),zzst(i,j+1), &
                    q_NOVA(j-2),q_NOVA(j-1),q_NOVA(j),q_NOVA(j+1), qdgn(jdgn),zzdgn(i,jdgn))
       enddo
      enddo
      
      do jdgn=mdgn_rs+1,mdgn
      do while(q_NOVA(j) .lt. qdgn(jdgn))
      j=j+1
      enddo
      if(j==1) j=j+2
      if(j==2) j=j+1
      call interp1d3l(psival_NOVA(j-2),psival_NOVA(j-1),psival_NOVA(j),psival_NOVA(j+1), &
                    q_NOVA(j-2),q_NOVA(j-1),q_NOVA(j),q_NOVA(j+1), qdgn(jdgn),psdgn(jdgn))
       do i=1,n2th+5
         call interp1d3l(xxst(i,j-2),xxst(i,j-1),xxst(i,j),xxst(i,j+1), &
                    q_NOVA(j-2),q_NOVA(j-1),q_NOVA(j),q_NOVA(j+1), qdgn(jdgn),xxdgn(i,jdgn))
         call interp1d3l(zzst(i,j-2),zzst(i,j-1),zzst(i,j),zzst(i,j+1), &
                    q_NOVA(j-2),q_NOVA(j-1),q_NOVA(j),q_NOVA(j+1), qdgn(jdgn),zzdgn(i,jdgn))
       enddo

      enddo  


      do jdgn=1,mdgn
      i=0
      do nrk=0,nsize-1
      ikdgn=0
      do jt=1,n2th+5     
!      if((xxs(jt,js) .ge. xxt(nrkx(nrk)*mxm+ix_min(nrk)-2).and.xxs(jt,js).le. xxt(nrkx(nrk)*mxm+ix_max(nrk)-2)) .and. &
!         (zzs(jt,js) .ge. zzt(nrkz(nrk)*mzm+iz_min(nrk)-2).and.zzs(jt,js).le. zzt(nrkz(nrk)*mzm+iz_max(nrk)-2)) ) then
      if((xxdgn(jt,jdgn) .ge. xxt(nrkx(nrk)*mxm+1).and.xxdgn(jt,jdgn).lt. xxt(nrkx(nrk)*mxm+mxm)+dxt(nrkx(nrk)*mxm+mxm)) .and. &
         (zzdgn(jt,jdgn) .ge. zzt(nrkz(nrk)*mzm+1).and.zzdgn(jt,jdgn).lt. zzt(nrkz(nrk)*mzm+mzm)+dzt(nrkz(nrk)*mzm+mzm)) ) then
        nrkdgn(jt,jdgn)=nrk
        if(jt.gt.3 .and. jt.le. 3+n2th) then
        ikdgn=ikdgn+1
        itdgn_nrk(ikdgn,nrk,jdgn)=jt
        endif
      endif
      enddo   
      mtdgn_nrk(nrk,jdgn)=ikdgn

      if(mtdgn_nrk(nrk,jdgn).gt.0) then
      i=i+1
      nranksend_dgn(i,jdgn)=nrk
      itdgnmin(nrk,jdgn)=itdgn_nrk(1,nrk,jdgn)
      itdgnmax(nrk,jdgn)=itdgn_nrk(mtdgn_nrk(nrk,jdgn),nrk,jdgn)       
      endif
      enddo
      mrkdgn(jdgn)=i
      enddo
     
     if(nrank.eq.nrank_dgn) then
     do jdgn=1,mdgn
     write(*,*) jdgn,'q=',qdgn(jdgn),'ps=',psdgn(jdgn),'mrk=',mrkdgn(jdgn)
     do ii=1,mrkdgn(jdgn)
     write(*,*) 'nrk=',nranksend_dgn(ii,jdgn),'it=',itdgnmin(nranksend_dgn(ii,jdgn),jdgn),itdgnmax(nranksend_dgn(ii,jdgn),jdgn),mtdgn_nrk(nranksend_dgn(ii,jdgn),jdgn)
     enddo
     enddo
     endif


      return
      end


!****************************************************************
      subroutine diagn_nmmode(fxz,jdg)
      USE DECLARE
      implicit real*8 (b-h,o-z)
      integer ms,me,jdg,ii,nfile,nfiler,nfilei,mtsend,ltmin,ltmax
       character*3 cn_dgn
       character*2 cn_ntor
       character*12 outputm,outputr,outputi
      dimension fxz(mx,mz,my),fst(n2th+5,my),fst_recv(n2th+5,myt)
!      real*8, dimension(n2th,my) :: fdgn
      real*8, dimension(n2th,myt) :: data_dgn      
      complex*16, dimension(n2th/2+1,myt) :: spec_dgn
!      complex*16, dimension(n2th/2+1,my) :: fdgn_spec
      include 'mpif.h'

       integer status(mpi_status_size)

!
!      do 2 jdg=1,mdgn
      if(mtdgn_nrk(nrank,jdg).gt.0) then
!      do 21 im=1,mm
      do 21 jy=1,my      
      call interp_xz2ps(fxz(ix_first:ix_last,iz_first:iz_last,jy), &
                        xx(ix_first:ix_last),zz(iz_first:iz_last),ix_last-ix_first+1,iz_last-iz_first+1, &
                       fst(itdgnmin(nrank,jdg):itdgnmax(nrank,jdg),jy), &
                        xxdgn(itdgnmin(nrank,jdg):itdgnmax(nrank,jdg),jdg),zzdgn(itdgnmin(nrank,jdg):itdgnmax(nrank,jdg),jdg),mtdgn_nrk(nrank,jdg))
      
   21 continue    
      endif       
!!mpi
!	write(*,*) nrank,jdg,"map1 done" 

         do ii=1,mrkdgn(jdg)
         ltmin=itdgnmin(nranksend_dgn(ii,jdg),jdg)
         ltmax=itdgnmax(nranksend_dgn(ii,jdg),jdg)
         mtsend=mtdgn_nrk(nranksend_dgn(ii,jdg),jdg)
         
         if(nrank .eq. nranksend_dgn(ii,jdg)) then
         CALL MPI_Send(fst(ltmin:ltmax,3:my-2),mtsend*mym, MPI_DOUBLE_PRECISION, nrank_dgn, ii,  &
		               MPI_COMM_WORLD,ierror )
         endif

         if(nrank .eq. nrank_dgn) then
         CALL MPI_Recv(fst_recv(ltmin:ltmax,1+nrky(nrank)*mym:mym+nrky(nrank)*mym),mtsend*mym, MPI_DOUBLE_PRECISION, nranksend_dgn(ii,jdg),ii,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         enddo

!	write(*,*) nrank,js,"mapsend done" 

!!mpi
!    2 continue 

      if(nrank .eq. nrank_dgn) then

      data_dgn(2:n2th,:)=fst_recv(4:n2th+2,:)
      data_dgn(1,:)=fst_recv(n2th+3,:)

!      do jdg=1,mdgn
!      do im=ms,me

      call dfftw_plan_dft_r2c_2d(plan,n2th,myt,data_dgn,spec_dgn,FFTW_ESTIMATE)
	  call dfftw_execute_dft_r2c(plan,data_dgn,spec_dgn)
	  call dfftw_destroy_plan(plan)

!      fdgn_spec(:,:,im)=spec_dgn(:,:)
!      enddo
!      enddo
      spec_dgn=spec_dgn/myt/n2th

      write(cn_dgn,'(i3.3)') jdg 
      do ntor=0,mtor_dgn
      write(cn_ntor,'(i2.2)') ntor
     
      minpol=ntor
      maxpol=min(int(qmax*ntor),16)
      mpol=maxpol-minpol+1
      nfile=610+ntor
      nfiler=640+ntor
      nfilei=670+ntor
      outputm='dgnm'//cn_dgn//'n'//cn_ntor
      outputr='dgnr'//cn_dgn//'n'//cn_ntor
      outputi='dgni'//cn_dgn//'n'//cn_ntor
      if(nstep==0) then
        open(unit=nfile,file=outputm,status='unknown',form='formatted')
        open(unit=nfiler,file=outputr,status='unknown',form='formatted')
        open(unit=nfilei,file=outputi,status='unknown',form='formatted') 
      else
        open(unit=nfile,file=outputm,status='unknown',form='formatted',position='append')
        open(unit=nfiler,file=outputr,status='unknown',form='formatted',position='append')
        open(unit=nfilei,file=outputi,status='unknown',form='formatted',position='append')
      endif
      write(UNIT=nfile,FMT=500) time,(sqrt(real(spec_dgn(npol+1,ntor+1))**2+aimag(spec_dgn(npol+1,ntor+1))**2),npol=minpol,maxpol)      
      write(UNIT=nfiler,FMT=500) time,(real(spec_dgn(npol+1,ntor+1)),npol=minpol,maxpol)      
      write(UNIT=nfilei,FMT=500) time,(aimag(spec_dgn(npol+1,ntor+1)),npol=minpol,maxpol)
      enddo
500   format(<mpol+1>(1x,e14.5E4))
      endif

      return
      end

!ws**********************************************************************
      subroutine diagn_max
      USE DECLARE
      include 'mpif.h'
      real*8 vsmax,Eymax,rhomax,pmax,cymax,crmax,cpmax,bymax,brmax,bpmax,vymax,vrmax,vpmax
      real*8 vsmax1,Eymax1,rhomax1,pmax1,cymax1,crmax1,cpmax1,bymax1,brmax1,bpmax1,vymax1,vrmax1,vpmax1
      real*8, dimension(mx,mz,my) :: vs
!  
!  d1fc= d f / dx  with fourth-order accuracy central difference
!      d1fc(fm2,fm1,f0,fp1,fp2,a,b,c,d)= &
!       a*(fp1-f0)+b*(f0-fm1)+c*(fp2-f0)+d*(f0-fm2)
!  define statement functions
!  d2fc= d2 f / dx2   with central difference
!      d2fc(fm1,f0,fp1,xm1,x0,xp1)= &
!       2.*((fp1-f0)/(xp1-x0)-(f0-fm1)/(x0-xm1))/(xp1-xm1)
!      call current(0)
!      call convt

      rhomax=maxval(x1(:,:,:,1))
        pmax=maxval(x1(:,:,:,2))
       vymax=maxval(x1(:,:,:,4))
       bymax=maxval(x1(:,:,:,7))
       cymax=maxval(cur(:,:,:,2))
       Eymax=maxval(Ef(:,:,:,2))
      

      do 40 jz=iz_first,iz_last
      do 40 jx=ix_first,ix_last
        cr(jx,jz,:)=cur(jx,jz,:,1)*wx2r(jx,jz)+cur(jx,jz,:,3)*wz2r(jx,jz)
   40 continue
      crmax=maxval(cr(:,:,:))

      
      do 50 jz=iz_first,iz_last
      do 50 jx=ix_first,ix_last
        cp(jx,jz,:)=cur(jx,jz,:,1)*wx2p(jx,jz)+cur(jx,jz,:,3)*wz2p(jx,jz)
   50 continue
      cpmax=maxval(cp(:,:,:))

      do 41 jz=iz_first,iz_last
      do 41 jx=ix_first,ix_last
        br(jx,jz,:)=x1(jx,jz,:,6)*wx2r(jx,jz)+x1(jx,jz,:,8)*wz2r(jx,jz)
   41 continue
      brmax=maxval(br(:,:,:))
     
      do 51 jz=iz_first,iz_last
      do 51 jx=ix_first,ix_last
        bp(jx,jz,:)=x1(jx,jz,:,6)*wx2p(jx,jz)+x1(jx,jz,:,8)*wz2p(jx,jz)
   51 continue
      bpmax=maxval(bp(:,:,:))
      
      do 42 jz=iz_first,iz_last
      do 42 jx=ix_first,ix_last
        vr(jx,jz,:)=x1(jx,jz,:,3)*wx2r(jx,jz)+x1(jx,jz,:,5)*wz2r(jx,jz)
   42 continue
      vrmax=maxval(vr(:,:,:))
     
      do 52 jz=iz_first,iz_last
      do 52 jx=ix_first,ix_last
        vp(jx,jz,:)=x1(jx,jz,:,3)*wx2p(jx,jz)+x1(jx,jz,:,5)*wz2p(jx,jz)
   52 continue
      vpmax=maxval(vp(:,:,:))

      do 53 jy=1,my
      do 53 jz=iz_first,iz_last
      do 53 jx=ix_first,ix_last
        vs(jx,jz,jy)=x1(jx,jz,jy,3)*x(jx,jz,jy,8)-x1(jx,jz,jy,5)*x(jx,jz,jy,6)
   53 continue
      vsmax=maxval(vs(:,:,:))
       

!      open(unit=15,file='bfield.dat',status='unknown',form='formatted',position='append')
!      write(15,1000)time,bymax,x_bymax,z_bymax,bxmax,x_bxmax,z_bxmax,bzmax,x_bzmax,z_bzmax,brmax,x_brmax,z_brmax,bpmax,x_bpmax,z_bpmax
!
!mpi   -----------------------------------------------------------------
      CALL MPI_ALLREDUCE(cymax,cymax1,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                       MPI_COMM_WORLD,IERROR)

      CALL MPI_ALLREDUCE(crmax,crmax1,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                       MPI_COMM_WORLD,IERROR)

      CALL MPI_ALLREDUCE(cpmax,cpmax1,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                       MPI_COMM_WORLD,IERROR)

!bfeild   -----------------------------------------------------------------
      CALL MPI_ALLREDUCE(bymax,bymax1,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                       MPI_COMM_WORLD,IERROR)

      CALL MPI_ALLREDUCE(brmax,brmax1,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                       MPI_COMM_WORLD,IERROR)

      CALL MPI_ALLREDUCE(bpmax,bpmax1,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                       MPI_COMM_WORLD,IERROR)

!v   -----------------------------------------------------------------
      CALL MPI_ALLREDUCE(vymax,vymax1,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                       MPI_COMM_WORLD,IERROR)

      CALL MPI_ALLREDUCE(vrmax,vrmax1,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                       MPI_COMM_WORLD,IERROR)

      CALL MPI_ALLREDUCE(vpmax,vpmax1,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                       MPI_COMM_WORLD,IERROR)

      CALL MPI_ALLREDUCE(vsmax,vsmax1,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                       MPI_COMM_WORLD,IERROR)

      CALL MPI_ALLREDUCE(Eymax,Eymax1,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(rhomax,rhomax1,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                       MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(pmax,pmax1,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                       MPI_COMM_WORLD,IERROR)


      if(nrank.eq.0) then
      write(11,1000)time,vsmax1,Eymax1,rhomax1,pmax1,cymax1,crmax1,cpmax1,bymax1,brmax1,bpmax1,vymax1,vrmax1,vpmax1
1000  format(14(1x,e12.5))

!      print*,cxmax1,cymax1,czmax1,cxmin1,cymin1,czmin1,time
      endif


!      call diagnatxmode

      return
      end
!ws**********************************************************************
      subroutine diagn_maxmin
      USE DECLARE
      USE DECLARE_OXpoint
      include 'mpif.h'
      real*8,dimension(13) :: wsmax,wsmax1,wsmin,wsmin1
      real*8, dimension(mx,mz,my) :: vs
!  
!  d1fc= d f / dx  with fourth-order accuracy central difference
!      d1fc(fm2,fm1,f0,fp1,fp2,a,b,c,d)= &
!       a*(fp1-f0)+b*(f0-fm1)+c*(fp2-f0)+d*(f0-fm2)
!  define statement functions
!  d2fc= d2 f / dx2   with central difference
!      d2fc(fm1,f0,fp1,xm1,x0,xp1)= &
!       2.*((fp1-f0)/(xp1-x0)-(f0-fm1)/(x0-xm1))/(xp1-xm1)
!      call current(0)
!      call convt


       vs(:,:,:)=x1(:,:,:,3)*x(:,:,:,8)-x1(:,:,:,5)*x(:,:,:,6)
      do jy=iy_first,iy_last
       cr(:,:,jy)=cur(:,:,jy,1)*wx2r(:,:)+cur(:,:,jy,3)*wz2r(:,:)
       cp(:,:,jy)=cur(:,:,jy,1)*wx2p(:,:)+cur(:,:,jy,3)*wz2p(:,:)
       br(:,:,jy)= x1(:,:,jy,6)*wx2r(:,:)+ x1(:,:,jy,8)*wz2r(:,:)
       bp(:,:,jy)= x1(:,:,jy,6)*wx2p(:,:)+ x1(:,:,jy,8)*wz2p(:,:)
       vr(:,:,jy)= x1(:,:,jy,3)*wx2r(:,:)+ x1(:,:,jy,5)*wz2r(:,:)
       vp(:,:,jy)= x1(:,:,jy,3)*wx2p(:,:)+ x1(:,:,jy,5)*wz2p(:,:)
      enddo

       wsmax(1) =maxval( vs(:,:,:))       
       wsmax(2) =maxval( Ef(:,:,:,2))
       wsmax(3) =maxval( x1(:,:,:,1))
       wsmax(4) =maxval( x1(:,:,:,2))       
       wsmax(5) =maxval(cur(:,:,:,2))
       wsmax(6) =maxval( cr(:,:,:))
       wsmax(7) =maxval( cp(:,:,:))
       wsmax(8) =maxval( x1(:,:,:,7))
       wsmax(9) =maxval( br(:,:,:))
       wsmax(10)=maxval( bp(:,:,:))
       wsmax(11)=maxval( x1(:,:,:,4))
       wsmax(12)=maxval( vr(:,:,:))
       wsmax(13)=maxval( vp(:,:,:))

       wsmin(1) =minval( vs(:,:,:))
       wsmin(2) =minval( Ef(:,:,:,2))
       wsmin(3) =minval( x1(:,:,:,1))
       wsmin(4) =minval( x1(:,:,:,2))
       wsmin(5) =minval(cur(:,:,:,2))
       wsmin(6) =minval( cr(:,:,:))
       wsmin(7) =minval( cp(:,:,:))
       wsmin(8) =minval( x1(:,:,:,7))
       wsmin(9) =minval( br(:,:,:))
       wsmin(10)=minval( bp(:,:,:))
       wsmin(11)=minval( x1(:,:,:,4))
       wsmin(12)=minval( vr(:,:,:))        
       wsmin(13)=minval( vp(:,:,:))

!      open(unit=15,file='bfield.dat',status='unknown',form='formatted',position='append')
!      write(15,1000)time,bymax,x_bymax,z_bymax,bxmax,x_bxmax,z_bxmax,bzmax,x_bzmax,z_bzmax,brmax,x_brmax,z_brmax,bpmax,x_bpmax,z_bpmax
!
!mpi   -----------------------------------------------------------------
!      CALL MPI_ALLREDUCE(wsmax,wsmax1,13,MPI_DOUBLE_PRECISION,MPI_MAX, &
!                       MPI_COMM_WORLD,IERROR)
!      CALL MPI_ALLREDUCE(wsmin,wsmin1,13,MPI_DOUBLE_PRECISION,MPI_MIN, &
!                       MPI_COMM_WORLD,IERROR)

      CALL MPI_ALLREDUCE(wsmax,wsmax1,13,MPI_DOUBLE_PRECISION,MPI_MAX, &
            MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(wsmin,wsmin1,13,MPI_DOUBLE_PRECISION,MPI_MIN, &
            MPI_COMM_WORLD,IERROR)


      if(nrank.eq.0) then
      write(111,1000)time,wsmax1(:)
      write(112,1000)time,wsmax1(:)

1000  format(14(1x,e12.5))

!      print*,cxmax1,cymax1,czmax1,cxmin1,cymin1,czmin1,time
      endif

      br_max=wsmax1(9)
!      call diagnatxmode

      return
      end
!ws**********************************************************************
      subroutine diagn_brmax0
      USE DECLARE
      USE DECLARE_OXpoint
      include 'mpif.h'
      real*8  brmax

      do jy=iy_first,iy_last
       br(:,:,jy)= x1(:,:,jy,6)*wx2r(:,:)+ x1(:,:,jy,8)*wz2r(:,:)
      enddo
       brmax =maxval( br(:,:,:))
      CALL MPI_ALLREDUCE(brmax,br_max0,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
            MPI_COMM_WORLD,IERROR)
      return
      end

!****************************************************************
      subroutine valb_atlastgrid(fxz,fst,mm,kk,ibnd)
      USE DECLARE
      implicit real*8 (b-h,o-z)
      integer mm,kk,ibnd
      dimension fxz(mx,mz,my,mm),fst(n2th+5,mps4:mps,my,mm),f1s(mbm_nrk,mps4:mps),wst(n2th+5)
      include 'mpif.h'

       integer status(mpi_status_size)

!
      do 2 js=mpsa-2,mpsa-nda,-1
      if(mts_nrk(js,nrank).gt.0) then
      do 21 m=1,mm
      do 21 jy=1,my      
      call interp_xz2ps(fxz(ix_first:ix_last,iz_first:iz_last,jy,m), &
                        xx(ix_first:ix_last),zz(iz_first:iz_last),ix_last-ix_first+1,iz_last-iz_first+1, &
                       fst(itsmin(js,nrank):itsmax(js,nrank),js,jy,m), &
                        xxs(itsmin(js,nrank):itsmax(js,nrank),js),zzs(itsmin(js,nrank):itsmax(js,nrank),js),mts_nrk(js,nrank))
      
   21 continue    
!!mpi
         do irecv=1,mrkb
         do isend=1,nsend(irecv,js)
         if(nrank.eq.nranksend(irecv,js,isend) .and. nrank.ne.nrkb(irecv)) then
         ltmin=ittransmin(irecv,js,isend)
         ltmax=ittransmax(irecv,js,isend)
         if (itbmin(nrkb(irecv))+n2th .le. itsmax(js,nranksend(irecv,js,isend))) then
         ltmin=ltmin+n2th
         ltmax=ltmax+n2th
         endif
         if (itbmax(nrkb(irecv))-n2th .ge. itsmin(js,nranksend(irecv,js,isend))) then
         ltmin=ltmin-n2th
         ltmax=ltmax-n2th
         endif
         do lt=ltmin,ltmax
         CALL MPI_Send(fst(lt,js,1:my,1:mm),my*mm, MPI_DOUBLE_PRECISION, nrkb(irecv), isend,  &
		               MPI_COMM_WORLD,ierror )
         enddo
         endif
         enddo
         enddo
      endif

      if(mb_nrk(nrank).gt.0) then  
         do irecv=1,mrkb
         do isend=1,nsend(irecv,js)
         if(nrank.eq.nrkb(irecv) .and. nrank.ne.nranksend(irecv,js,isend)) then
         ltmin=ittransmin(irecv,js,isend)
         ltmax=ittransmax(irecv,js,isend)
         do lt=ltmin,ltmax
         CALL MPI_Recv(fst(lt,js,1:my,1:mm),my*mm, MPI_DOUBLE_PRECISION, nranksend(irecv,js,isend), isend,  &
		               MPI_COMM_WORLD,status,ierror )
         enddo
         endif
         enddo
         enddo
       endif
!!mpi
    2 continue 
      
      if(mb_nrk(nrank).gt.0) then     

      do 10 js=mpsa-2,mpsa-4,-1
!      call smth_ps_nrk(fst(itbmin(nrank):itbmax(nrank),js,:,:),itbmax(nrank)-itbmin(nrank)+1,3)
!ws_smps--------------------------
      do 11 k=1,kk
      do m=1,mm
      do jy=iy_first+2,iy_last-2
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      wst(jt)=fst(jt,js,jy,m)/2.+(fst(jt+1,js,jy,m)+fst(jt-1,js,jy,m))*3./16+(fst(jt+2,js,jy,m)+fst(jt-2,js,jy,m))/16.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,js,jy,m)=wst(jt)
      enddo

      enddo
      enddo

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         CALL MPI_Send(fst(ltmin:ltmax,js,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,js,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         CALL MPI_Send(fst(ltmin:ltmax,js,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,js,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         CALL MPI_Send(fst(ltmin:ltmax,js,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 2,  &
		               MPI_COMM_WORLD,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,js,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         CALL MPI_Send(fst(ltmin:ltmax,js,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,js,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   11 continue
!ws_smps--------------------
   10 continue

      if(ibnd==1) then
      is=mpsa-1
      do 31 m=1,mm    
      do 31 jy=1,my      
      do jt=itbmin(nrank),itbmax(nrank)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))  
!       fst(jt,is)=fst(jt,is-1)                   
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   31 continue 
!      call smth_ps_nrk(fst(itbmin(nrank):itbmax(nrank),is,:,:),itbmax(nrank)-itbmin(nrank)+1,3)
     do 32 k=1,kk
      do m=1,mm
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      wst(jt)=fst(jt,is,jy,m)/2.+(fst(jt+1,is,jy,m)+fst(jt-1,is,jy,m))*3./16+(fst(jt+2,is,jy,m)+fst(jt-2,is,jy,m))/16.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,is,jy,m)=wst(jt)
      enddo

      enddo
      enddo

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 2,  &
		               MPI_COMM_WORLD,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   32 continue
!ws_smps--------------------
      is=mpsa
      do 41 m=1,mm    
      do 41 jy=1,my  
      do jt=itbmin(nrank),itbmax(nrank)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))
!       fst(jt,is)=fst(jt,is-1)                   
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   41 continue
!      call smth_ps_nrk(fst(itbmin(nrank):itbmax(nrank),is,:,:),itbmax(nrank)-itbmin(nrank)+1,3)
!ws_smps--------------------
     do 42 k=1,kk
      do m=1,mm
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      wst(jt)=fst(jt,is,jy,m)/2.+(fst(jt+1,is,jy,m)+fst(jt-1,is,jy,m))*3./16+(fst(jt+2,is,jy,m)+fst(jt-2,is,jy,m))/16.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,is,jy,m)=wst(jt)
      enddo

      enddo
      enddo

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 2,  &
		               MPI_COMM_WORLD,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   42 continue
!ws_smps--------------------
      endif !!ibnd==1
      
      if(ibnd==0) then
      is=mpsa
      do 51 m=1,mm    
      do 51 jy=1,my  
      do jt=itbmin(nrank),itbmax(nrank)
       fst(jt,is,jy,m)=0.           
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   51 continue

      is=mpsa-1
      do 52 m=1,mm    
      do 52 jy=1,my  
      do jt=itbmin(nrank),itbmax(nrank)
      call interp1d3l(fst(jt,is-3,jy,m),fst(jt,is-2,jy,m),fst(jt,is-1,jy,m),fst(jt,is+1,jy,m), &
                        ps(is-3),ps(is-2),ps(is-1),ps(is+1),ps(is),fst(jt,is,jy,m))                     
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   52 continue

      do 53 k=1,kk
      do m=1,mm
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      wst(jt)=fst(jt,is,jy,m)/2.+(fst(jt+1,is,jy,m)+fst(jt-1,is,jy,m))*3./16+(fst(jt+2,is,jy,m)+fst(jt-2,is,jy,m))/16.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,is,jy,m)=wst(jt)
      enddo

      enddo
      enddo

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 2,  &
		               MPI_COMM_WORLD,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, MPI_DOUBLE_PRECISION, nrkb1(1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   53 continue
      endif !!ibnd==1 

      do 3 m=1,mm    
      do 3 jy=1,my           
      do ikb=1,mb_nrk(nrank)
        i=itb_nrk(ikb,nrank)
        ib=ib_nrk(ikb,nrank)
        do js=mpsa,mpsa-3,-1     
        call interp1d3l(fst(i-2,js,jy,m),fst(i-1,js,jy,m),fst(i,js,jy,m),fst(i+1,js,jy,m), &
                        thst(i-2),thst(i-1),thst(i),thst(i+1),thb(ib),f1s(ikb,js))
        enddo
        js=mpsa
        call interp1d3l(f1s(ikb,js),f1s(ikb,js-1),f1s(ikb,js-2),f1s(ikb,js-3), &
                        ps(js),ps(js-1),ps(js-2),ps(js-3),psb(ib),fxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,m))
!        fxz(jbx(ib),jbz(ib))=(asm_b(ib)*x1s(ib,js-1)+bsm_b(ib)*x1s(ib,js-2)+csm_b(ib)*x1s(ib,js-3)) &
!                 /(asm_b(ib)+bsm_b(ib)+csm_b(ib))
      enddo
    3 continue 
      endif !!mb_nrk(nrank).gt.0
           
!      if(ibnd==2) then
!      
!      do ika1=1,ma1_nrk(nrank)
!        i=ita1_nrk(ika1,nrank)
!        ia1=ia1_nrk(ika1,nrank)
!        do js=mpsa,mpsa-3,1
!        call interp1d3l(fst(i-2,js),fst(i-1,js),fst(i,js),fst(i+1,js), &
!                        thst(i-2),thst(i-1),thst(i),thst(i+1),tha1(ia1),x1s(ika1,js,jy))
!        enddo
!        js=mpsa
!        call interp1d3l(x1s(ika1,js,jy),x1s(ika1,js-1,jy),x1s(ika1,js-2,jy),x1s(ika1,js-3,jy), &
!                        ps(js),ps(js-1),ps(js-2),ps(js-3),psa1(ia1),fxz(jxa1_nrk(ika1,nrank),jza1_nrk(ika1,nrank)))
!!        fxz(jxa1(ia1),jza1(ia1))=(asm_a1(ia1)*x1s(ia1,js-1)+bsm_a1(ia1)*x1s(ia1,js-2)+csm_a1(ia1)*x1s(ia1,js-3)) &
!!                  /(asm_a1(ia1)+bsm_a1(ia1)+csm_a1(ia1)) 
!      enddo
!      endif
      return
      end
!****************************************************************
!****************************************************************
      subroutine x1_atlastgrid_r0p1_v1(kk)
      USE DECLARE
      implicit real*8 (b-h,o-z)
      integer kk
      dimension fst(n2th+5,mps4:mps,my,8),f1s(mbm_nrk,mps4:mps),fsxz(mbm_nrk,8),wst(n2th+5) !fxz(mx,mz,my,mm),
      include 'mpif.h'

       integer status(mpi_status_size)

!
      do 2 js=mpsa-2,mpsa-nda,-1
      if(mts_nrk(js,nrank).gt.0) then
      do 21 m=1,8
      do 21 jy=1,my      
      call interp_xz2ps(x1(ix_first:ix_last,iz_first:iz_last,jy,m), &
                        xx(ix_first:ix_last),zz(iz_first:iz_last),ix_last-ix_first+1,iz_last-iz_first+1, &
                       fst(itsmin(js,nrank):itsmax(js,nrank),js,jy,m), &
                        xxs(itsmin(js,nrank):itsmax(js,nrank),js),zzs(itsmin(js,nrank):itsmax(js,nrank),js),mts_nrk(js,nrank))
      
   21 continue    
!!mpi
         do irecv=1,mrkb
         do isend=1,nsend(irecv,js)
         if(nrank.eq.nranksend(irecv,js,isend) .and. nrank.ne.nrkb(irecv)) then
         ltmin=ittransmin(irecv,js,isend)
         ltmax=ittransmax(irecv,js,isend)
         if (itbmin(nrkb(irecv))+n2th .le. itsmax(js,nranksend(irecv,js,isend))) then
         ltmin=ltmin+n2th
         ltmax=ltmax+n2th
         endif
         if (itbmax(nrkb(irecv))-n2th .ge. itsmin(js,nranksend(irecv,js,isend))) then
         ltmin=ltmin-n2th
         ltmax=ltmax-n2th
         endif
         do lt=ltmin,ltmax
         CALL MPI_Send(fst(lt,js,1:my,:),my*8, MPI_DOUBLE_PRECISION, nrkb(irecv), isend,  &
		               MPI_COMM_WORLD,ierror )
         enddo
         endif
         enddo
         enddo
      endif

      if(mb_nrk(nrank).gt.0) then  
         do irecv=1,mrkb
         do isend=1,nsend(irecv,js)
         if(nrank.eq.nrkb(irecv) .and. nrank.ne.nranksend(irecv,js,isend)) then
         ltmin=ittransmin(irecv,js,isend)
         ltmax=ittransmax(irecv,js,isend)
         do lt=ltmin,ltmax
         CALL MPI_Recv(fst(lt,js,1:my,:),my*8, MPI_DOUBLE_PRECISION, nranksend(irecv,js,isend), isend,  &
		               MPI_COMM_WORLD,status,ierror )
         enddo
         endif
         enddo
         enddo
       endif
!!mpi
    2 continue 
           
      if(mb_nrk(nrank).gt.0) then      

      do jy=1,my
      do js=mpsa-2,mpsa-4,-1
      do jt=itbmin(nrank),itbmax(nrank)
      vx1st=fst(jt,js,jy,3)
      vz1st=fst(jt,js,jy,5)
      bx1st=fst(jt,js,jy,6)
      bz1st=fst(jt,js,jy,8)
      fst(jt,js,jy,3)=vx1st*dcos(thst(jt))+vz1st*dsin(thst(jt))
      fst(jt,js,jy,5)=-vx1st*dsin(thst(jt))+vz1st*dcos(thst(jt))
      fst(jt,js,jy,6)=bx1st*dcos(thst(jt))+bz1st*dsin(thst(jt))
      fst(jt,js,jy,8)=-bx1st*dsin(thst(jt))+bz1st*dcos(thst(jt))
       
!      vrst=fst(jt,js,jy,3)*dcos(thst(jt))+fst(jt,js,jy,5)*dsin(thst(jt))
!      vpst=-fst(jt,js,jy,3)*dsin(thst(jt))+fst(jt,js,jy,5)*dcos(thst(jt))
!      brst=fst(jt,js,jy,6)*dcos(thst(jt))+fst(jt,js,jy,8)*dsin(thst(jt))
!      bpst=-fst(jt,js,jy,6)*dsin(thst(jt))+fst(jt,js,jy,8)*dcos(thst(jt))
!      fst(jt,js,jy,3)=vrst
!      fst(jt,js,jy,5)=vpst
!      fst(jt,js,jy,6)=brst
!      fst(jt,js,jy,8)=bpst
      enddo
      enddo
      enddo    

      do 10 js=mpsa-2,mpsa-4,-1
!      call smth_ps_nrk(fst(itbmin(nrank):itbmax(nrank),js,:,:),itbmax(nrank)-itbmin(nrank)+1,3)
!ws_smps--------------------------
      do 11 k=1,kk
      do m=1,8
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      wst(jt)=fst(jt,js,jy,m)/2.+((fst(jt+1,js,jy,m)+fst(jt-1,js,jy,m))*3./16+(fst(jt+2,js,jy,m)+fst(jt-2,js,jy,m))/16.&
             +(fst(jt,js,jy+1,m)+fst(jt,js,jy-1,m))*3./16+(fst(jt,js,jy+2,m)+fst(jt,js,jy-2,m))/16.)/2.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,js,jy,m)=wst(jt)
      enddo

      enddo
      enddo

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         CALL MPI_Send(fst(ltmin:ltmax,js,1:my,:), 2*my*8, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,js,1:my,:), 2*my*8, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         CALL MPI_Send(fst(ltmin:ltmax,js,1:my,:), 2*my*8, MPI_DOUBLE_PRECISION, nrkb1(1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,js,1:my,:), 2*my*8, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         CALL MPI_Send(fst(ltmin:ltmax,js,1:my,:), 2*my*8, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 2,  &
		               MPI_COMM_WORLD,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,js,1:my,:), 2*my*8, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         CALL MPI_Send(fst(ltmin:ltmax,js,1:my,:), 2*my*8, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,js,1:my,:), 2*my*8, MPI_DOUBLE_PRECISION, nrkb1(1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   11 continue
!ws_smps--------------------
   10 continue

      do 61 m=1,8

      if(m.eq.3 .or. m.eq.6) then !!(m=3,6) ibnd==0
      is=mpsa 
      do 51 jy=1,my  
      do jt=itbmin(nrank),itbmax(nrank)
       fst(jt,is,jy,m)=0.           
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   51 continue

      is=mpsa-1 
      do 52 jy=1,my  
      do jt=itbmin(nrank),itbmax(nrank)
        call interp1d3l(fst(jt,is-3,jy,m),fst(jt,is-2,jy,m),fst(jt,is-1,jy,m),fst(jt,is+1,jy,m), &
                        ps(is-3),ps(is-2),ps(is-1),ps(is+1),ps(is),fst(jt,is,jy,m))                   
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   52 continue

      do 53 k=1,kk
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
        wst(jt)=fst(jt,is,jy,m)/2.+((fst(jt+1,is,jy,m)+fst(jt-1,is,jy,m))*3./16+(fst(jt+2,is,jy,m)+fst(jt-2,is,jy,m))/16.&
             +(fst(jt,is,jy+1,m)+fst(jt,is,jy-1,m))*3./16+(fst(jt,is,jy+2,m)+fst(jt,is,jy-2,m))/16.)/2.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,is,jy,m)=wst(jt)
      enddo

      enddo

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 2,  &
		               MPI_COMM_WORLD,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   53 continue

      else  !!(m=1,2,4,5,7,8,) bnd==1  
      is=mpsa-1    
      do 31 jy=1,my      
      do jt=itbmin(nrank),itbmax(nrank)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))  
!       fst(jt,is)=fst(jt,is-1)                   
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   31 continue 
!      call smth_ps_nrk(fst(itbmin(nrank):itbmax(nrank),is,:,:),itbmax(nrank)-itbmin(nrank)+1,3)
      do 32 k=1,kk
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
       wst(jt)=fst(jt,is,jy,m)/2.+((fst(jt+1,is,jy,m)+fst(jt-1,is,jy,m))*3./16+(fst(jt+2,is,jy,m)+fst(jt-2,is,jy,m))/16.&
             +(fst(jt,is,jy+1,m)+fst(jt,is,jy-1,m))*3./16+(fst(jt,is,jy+2,m)+fst(jt,is,jy-2,m))/16.)/2.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,is,jy,m)=wst(jt)
      enddo

      enddo
      
         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 2,  &
		               MPI_COMM_WORLD,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   32 continue
!ws_smps--------------------
      is=mpsa
 
      do 41 jy=1,my  
      do jt=itbmin(nrank),itbmax(nrank)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))
!       fst(jt,is)=fst(jt,is-1)                   
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   41 continue
!      call smth_ps_nrk(fst(itbmin(nrank):itbmax(nrank),is,:,:),itbmax(nrank)-itbmin(nrank)+1,3)
!ws_smps--------------------
      do 42 k=1,kk
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      wst(jt)=fst(jt,is,jy,m)/2.+((fst(jt+1,is,jy,m)+fst(jt-1,is,jy,m))*3./16+(fst(jt+2,is,jy,m)+fst(jt-2,is,jy,m))/16.&
             +(fst(jt,is,jy+1,m)+fst(jt,is,jy-1,m))*3./16+(fst(jt,is,jy+2,m)+fst(jt,is,jy-2,m))/16.)/2.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,is,jy,m)=wst(jt)
      enddo

      enddo

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 2,  &
		               MPI_COMM_WORLD,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   42 continue
!ws_smps--------------------
      endif !!(m=3,6)ibnd==0
   61 continue

      do 3 jy=1,my           
      do ikb=1,mb_nrk(nrank)
        i=itb_nrk(ikb,nrank)
        ib=ib_nrk(ikb,nrank) 
        do m=1,8 
        do js=mpsa,mpsa-3,-1     
        call interp1d3l(fst(i-2,js,jy,m),fst(i-1,js,jy,m),fst(i,js,jy,m),fst(i+1,js,jy,m), &
                        thst(i-2),thst(i-1),thst(i),thst(i+1),thb(ib),f1s(ikb,js))
        enddo
        js=mpsa
        call interp1d3l(f1s(ikb,js),f1s(ikb,js-1),f1s(ikb,js-2),f1s(ikb,js-3), &
                        ps(js),ps(js-1),ps(js-2),ps(js-3),psb(ib),fsxz(ikb,m))
        enddo
!        fxz(jbx(ib),jbz(ib))=(asm_b(ib)*x1s(ib,js-1)+bsm_b(ib)*x1s(ib,js-2)+csm_b(ib)*x1s(ib,js-3)) &
!                 /(asm_b(ib)+bsm_b(ib)+csm_b(ib))
        x1(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,1)=fsxz(ikb,1)
        x1(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,2)=fsxz(ikb,2)
        x1(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,3)=fsxz(ikb,3)*dcos(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank))) &
                                                      -fsxz(ikb,5)*dsin(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank)))
        x1(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,4)=fsxz(ikb,4)
        x1(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,5)=fsxz(ikb,3)*dsin(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank))) &
                                                      +fsxz(ikb,5)*dcos(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank)))
        x1(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,6)=fsxz(ikb,6)*dcos(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank))) &
                                                      -fsxz(ikb,8)*dsin(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank)))
        x1(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,7)=fsxz(ikb,7)
        x1(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,8)=fsxz(ikb,6)*dsin(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank))) &
                                                      +fsxz(ikb,8)*dcos(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank)))
      enddo
    3 continue 
      endif !!mb_nrk(nrank).gt.0
           
!      if(ibnd==2) then
!      
!      do ika1=1,ma1_nrk(nrank)
!        i=ita1_nrk(ika1,nrank)
!        ia1=ia1_nrk(ika1,nrank)
!        do js=mpsa,mpsa-3,1
!        call interp1d3l(fst(i-2,js),fst(i-1,js),fst(i,js),fst(i+1,js), &
!                        thst(i-2),thst(i-1),thst(i),thst(i+1),tha1(ia1),x1s(ika1,js,jy))
!        enddo
!        js=mpsa
!        call interp1d3l(x1s(ika1,js,jy),x1s(ika1,js-1,jy),x1s(ika1,js-2,jy),x1s(ika1,js-3,jy), &
!                        ps(js),ps(js-1),ps(js-2),ps(js-3),psa1(ia1),fxz(jxa1_nrk(ika1,nrank),jza1_nrk(ika1,nrank)))
!!        fxz(jxa1(ia1),jza1(ia1))=(asm_a1(ia1)*x1s(ia1,js-1)+bsm_a1(ia1)*x1s(ia1,js-2)+csm_a1(ia1)*x1s(ia1,js-3)) &
!!                  /(asm_a1(ia1)+bsm_a1(ia1)+csm_a1(ia1)) 
!      enddo
!      endif
      return
      end

!ws************************************************************
!****************************************************************
      subroutine Ef_atlastgrid_r1p0_v1(kk)
      USE DECLARE
      implicit real*8 (b-h,o-z)
      integer kk
      dimension fst(n2th+5,mps4:mps,my,3),f1s(mbm_nrk,mps4:mps),fsxz(mbm_nrk,3),wst(n2th+5) !fxz(mx,mz,my,mm),
      include 'mpif.h'

       integer status(mpi_status_size)

!
      do 2 js=mpsa-2,mpsa-nda,-1
      if(mts_nrk(js,nrank).gt.0) then
      do 21 m=1,3
      do 21 jy=1,my      
      call interp_xz2ps(Ef(ix_first:ix_last,iz_first:iz_last,jy,m), &
                        xx(ix_first:ix_last),zz(iz_first:iz_last),ix_last-ix_first+1,iz_last-iz_first+1, &
                       fst(itsmin(js,nrank):itsmax(js,nrank),js,jy,m), &
                        xxs(itsmin(js,nrank):itsmax(js,nrank),js),zzs(itsmin(js,nrank):itsmax(js,nrank),js),mts_nrk(js,nrank))
      
   21 continue    
!!mpi
         do irecv=1,mrkb
         do isend=1,nsend(irecv,js)
         if(nrank.eq.nranksend(irecv,js,isend) .and. nrank.ne.nrkb(irecv)) then
         ltmin=ittransmin(irecv,js,isend)
         ltmax=ittransmax(irecv,js,isend)
         if (itbmin(nrkb(irecv))+n2th .le. itsmax(js,nranksend(irecv,js,isend))) then
         ltmin=ltmin+n2th
         ltmax=ltmax+n2th
         endif
         if (itbmax(nrkb(irecv))-n2th .ge. itsmin(js,nranksend(irecv,js,isend))) then
         ltmin=ltmin-n2th
         ltmax=ltmax-n2th
         endif
         do lt=ltmin,ltmax
         CALL MPI_Send(fst(lt,js,1:my,:),my*3, MPI_DOUBLE_PRECISION, nrkb(irecv), isend,  &
		               MPI_COMM_WORLD,ierror )
         enddo
         endif
         enddo
         enddo
      endif

      if(mb_nrk(nrank).gt.0) then  
         do irecv=1,mrkb
         do isend=1,nsend(irecv,js)
         if(nrank.eq.nrkb(irecv) .and. nrank.ne.nranksend(irecv,js,isend)) then
         ltmin=ittransmin(irecv,js,isend)
         ltmax=ittransmax(irecv,js,isend)
         do lt=ltmin,ltmax
         CALL MPI_Recv(fst(lt,js,1:my,:),my*3, MPI_DOUBLE_PRECISION, nranksend(irecv,js,isend), isend,  &
		               MPI_COMM_WORLD,status,ierror )
         enddo
         endif
         enddo
         enddo
       endif
!!mpi
    2 continue 
           
      if(mb_nrk(nrank).gt.0) then      

      do jy=1,my
      do js=mpsa-2,mpsa-4,-1
      do jt=itbmin(nrank),itbmax(nrank)
      Exst=fst(jt,js,jy,1)
      Ezst=fst(jt,js,jy,3)

      fst(jt,js,jy,1)=Exst*dcos(thst(jt))+Ezst*dsin(thst(jt))
      fst(jt,js,jy,3)=-Exst*dsin(thst(jt))+Ezst*dcos(thst(jt))
       
!      vrst=fst(jt,js,jy,3)*dcos(thst(jt))+fst(jt,js,jy,5)*dsin(thst(jt))
!      vpst=-fst(jt,js,jy,3)*dsin(thst(jt))+fst(jt,js,jy,5)*dcos(thst(jt))
!      brst=fst(jt,js,jy,6)*dcos(thst(jt))+fst(jt,js,jy,8)*dsin(thst(jt))
!      bpst=-fst(jt,js,jy,6)*dsin(thst(jt))+fst(jt,js,jy,8)*dcos(thst(jt))
!      fst(jt,js,jy,3)=vrst
!      fst(jt,js,jy,5)=vpst
!      fst(jt,js,jy,6)=brst
!      fst(jt,js,jy,8)=bpst
      enddo
      enddo
      enddo    

      do 10 js=mpsa-2,mpsa-4,-1
!      call smth_ps_nrk(fst(itbmin(nrank):itbmax(nrank),js,:,:),itbmax(nrank)-itbmin(nrank)+1,3)
!ws_smps--------------------------
      do 11 k=1,kk
      do m=1,3
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      wst(jt)=fst(jt,js,jy,m)/2.+((fst(jt+1,js,jy,m)+fst(jt-1,js,jy,m))*3./16+(fst(jt+2,js,jy,m)+fst(jt-2,js,jy,m))/16.&
             +(fst(jt,js,jy+1,m)+fst(jt,js,jy-1,m))*3./16+(fst(jt,js,jy+2,m)+fst(jt,is,jy-2,m))/16.)/2.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,js,jy,m)=wst(jt)
      enddo

      enddo
      enddo

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         CALL MPI_Send(fst(ltmin:ltmax,js,1:my,:), 2*my*3, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,js,1:my,:), 2*my*3, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         CALL MPI_Send(fst(ltmin:ltmax,js,1:my,:), 2*my*3, MPI_DOUBLE_PRECISION, nrkb1(1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,js,1:my,:), 2*my*3, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         CALL MPI_Send(fst(ltmin:ltmax,js,1:my,:), 2*my*3, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 2,  &
		               MPI_COMM_WORLD,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,js,1:my,:), 2*my*3, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         CALL MPI_Send(fst(ltmin:ltmax,js,1:my,:), 2*my*3, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,js,1:my,:), 2*my*3, MPI_DOUBLE_PRECISION, nrkb1(1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   11 continue
!ws_smps--------------------
   10 continue

      do 61 m=1,3

      if(m.eq.2 .or. m.eq.3) then !!(m=2,3) ibnd==0
      is=mpsa 
      do 51 jy=1,my  
      do jt=itbmin(nrank),itbmax(nrank)
       fst(jt,is,jy,m)=0.           
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   51 continue

      is=mpsa-1 
      do 52 jy=1,my  
      do jt=itbmin(nrank),itbmax(nrank)
        call interp1d3l(fst(jt,is-3,jy,m),fst(jt,is-2,jy,m),fst(jt,is-1,jy,m),fst(jt,is+1,jy,m), &
                        ps(is-3),ps(is-2),ps(is-1),ps(is+1),ps(is),fst(jt,is,jy,m))                   
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   52 continue

      do 53 k=1,kk
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      wst(jt)=fst(jt,is,jy,m)/2.+((fst(jt+1,is,jy,m)+fst(jt-1,is,jy,m))*3./16+(fst(jt+2,is,jy,m)+fst(jt-2,is,jy,m))/16.&
             +(fst(jt,is,jy+1,m)+fst(jt,is,jy-1,m))*3./16+(fst(jt,is,jy+2,m)+fst(jt,is,jy-2,m))/16.)/2.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,is,jy,m)=wst(jt)
      enddo

      enddo

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 2,  &
		               MPI_COMM_WORLD,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   53 continue

      else  !!(m=1) bnd==1  
      is=mpsa-1    
      do 31 jy=1,my      
      do jt=itbmin(nrank),itbmax(nrank)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))  
!       fst(jt,is)=fst(jt,is-1)                   
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   31 continue 
!      call smth_ps_nrk(fst(itbmin(nrank):itbmax(nrank),is,:,:),itbmax(nrank)-itbmin(nrank)+1,3)
      do 32 k=1,kk
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      wst(jt)=fst(jt,is,jy,m)/2.+((fst(jt+1,is,jy,m)+fst(jt-1,is,jy,m))*3./16+(fst(jt+2,is,jy,m)+fst(jt-2,is,jy,m))/16.&
             +(fst(jt,is,jy+1,m)+fst(jt,is,jy-1,m))*3./16+(fst(jt,is,jy+2,m)+fst(jt,is,jy-2,m))/16.)/2.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,is,jy,m)=wst(jt)
      enddo

      enddo
      
         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 2,  &
		               MPI_COMM_WORLD,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   32 continue
!ws_smps--------------------
      is=mpsa
 
      do 41 jy=1,my  
      do jt=itbmin(nrank),itbmax(nrank)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))
!       fst(jt,is)=fst(jt,is-1)                   
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   41 continue
!      call smth_ps_nrk(fst(itbmin(nrank):itbmax(nrank),is,:,:),itbmax(nrank)-itbmin(nrank)+1,3)
!ws_smps--------------------
      do 42 k=1,kk
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      wst(jt)=fst(jt,is,jy,m)/2.+((fst(jt+1,is,jy,m)+fst(jt-1,is,jy,m))*3./16+(fst(jt+2,is,jy,m)+fst(jt-2,is,jy,m))/16.&
             +(fst(jt,is,jy+1,m)+fst(jt,is,jy-1,m))*3./16+(fst(jt,is,jy+2,m)+fst(jt,is,jy-2,m))/16.)/2.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,is,jy,m)=wst(jt)
      enddo

      enddo

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 2,  &
		               MPI_COMM_WORLD,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   42 continue
!ws_smps--------------------
      endif !!(m=2,3)ibnd==0
   61 continue

      do 3 jy=1,my           
      do ikb=1,mb_nrk(nrank)
        i=itb_nrk(ikb,nrank)
        ib=ib_nrk(ikb,nrank) 
        do m=1,3 
        do js=mpsa,mpsa-3,-1     
        call interp1d3l(fst(i-2,js,jy,m),fst(i-1,js,jy,m),fst(i,js,jy,m),fst(i+1,js,jy,m), &
                        thst(i-2),thst(i-1),thst(i),thst(i+1),thb(ib),f1s(ikb,js))
        enddo
        js=mpsa
        call interp1d3l(f1s(ikb,js),f1s(ikb,js-1),f1s(ikb,js-2),f1s(ikb,js-3), &
                        ps(js),ps(js-1),ps(js-2),ps(js-3),psb(ib),fsxz(ikb,m))
        enddo
!        fxz(jbx(ib),jbz(ib))=(asm_b(ib)*x1s(ib,js-1)+bsm_b(ib)*x1s(ib,js-2)+csm_b(ib)*x1s(ib,js-3)) &
!                 /(asm_b(ib)+bsm_b(ib)+csm_b(ib))

        Ef(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,1)=fsxz(ikb,1)*dcos(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank))) &
                                                      -fsxz(ikb,3)*dsin(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank)))
        Ef(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,2)=fsxz(ikb,2)
        Ef(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,3)=fsxz(ikb,1)*dsin(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank))) &
                                                      +fsxz(ikb,3)*dcos(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank)))

      enddo
    3 continue 
      endif !!mb_nrk(nrank).gt.0
           
!      if(ibnd==2) then
!      
!      do ika1=1,ma1_nrk(nrank)
!        i=ita1_nrk(ika1,nrank)
!        ia1=ia1_nrk(ika1,nrank)
!        do js=mpsa,mpsa-3,1
!        call interp1d3l(fst(i-2,js),fst(i-1,js),fst(i,js),fst(i+1,js), &
!                        thst(i-2),thst(i-1),thst(i),thst(i+1),tha1(ia1),x1s(ika1,js,jy))
!        enddo
!        js=mpsa
!        call interp1d3l(x1s(ika1,js,jy),x1s(ika1,js-1,jy),x1s(ika1,js-2,jy),x1s(ika1,js-3,jy), &
!                        ps(js),ps(js-1),ps(js-2),ps(js-3),psa1(ia1),fxz(jxa1_nrk(ika1,nrank),jza1_nrk(ika1,nrank)))
!!        fxz(jxa1(ia1),jza1(ia1))=(asm_a1(ia1)*x1s(ia1,js-1)+bsm_a1(ia1)*x1s(ia1,js-2)+csm_a1(ia1)*x1s(ia1,js-3)) &
!!                  /(asm_a1(ia1)+bsm_a1(ia1)+csm_a1(ia1)) 
!      enddo
!      endif
      return
      end
!****************************************************************
      subroutine cur_atlastgrid_r0p1_v1(kk)
      USE DECLARE
      implicit real*8 (b-h,o-z)
      integer kk
      dimension fst(n2th+5,mps4:mps,my,3),f1s(mbm_nrk,mps4:mps),fsxz(mbm_nrk,3),wst(n2th+5) !fxz(mx,mz,my,mm),
      include 'mpif.h'

       integer status(mpi_status_size)

!
      do 2 js=mpsa-2,mpsa-nda,-1
      if(mts_nrk(js,nrank).gt.0) then
      do 21 m=1,3
      do 21 jy=1,my      
      call interp_xz2ps(cur(ix_first:ix_last,iz_first:iz_last,jy,m), &
                        xx(ix_first:ix_last),zz(iz_first:iz_last),ix_last-ix_first+1,iz_last-iz_first+1, &
                       fst(itsmin(js,nrank):itsmax(js,nrank),js,jy,m), &
                        xxs(itsmin(js,nrank):itsmax(js,nrank),js),zzs(itsmin(js,nrank):itsmax(js,nrank),js),mts_nrk(js,nrank))
      
   21 continue    
!!mpi
         do irecv=1,mrkb
         do isend=1,nsend(irecv,js)
         if(nrank.eq.nranksend(irecv,js,isend) .and. nrank.ne.nrkb(irecv)) then
         ltmin=ittransmin(irecv,js,isend)
         ltmax=ittransmax(irecv,js,isend)
         if (itbmin(nrkb(irecv))+n2th .le. itsmax(js,nranksend(irecv,js,isend))) then
         ltmin=ltmin+n2th
         ltmax=ltmax+n2th
         endif
         if (itbmax(nrkb(irecv))-n2th .ge. itsmin(js,nranksend(irecv,js,isend))) then
         ltmin=ltmin-n2th
         ltmax=ltmax-n2th
         endif
         do lt=ltmin,ltmax
         CALL MPI_Send(fst(lt,js,1:my,:),my*3, MPI_DOUBLE_PRECISION, nrkb(irecv), isend,  &
		               MPI_COMM_WORLD,ierror )
         enddo
         endif
         enddo
         enddo
      endif

      if(mb_nrk(nrank).gt.0) then  
         do irecv=1,mrkb
         do isend=1,nsend(irecv,js)
         if(nrank.eq.nrkb(irecv) .and. nrank.ne.nranksend(irecv,js,isend)) then
         ltmin=ittransmin(irecv,js,isend)
         ltmax=ittransmax(irecv,js,isend)
         do lt=ltmin,ltmax
         CALL MPI_Recv(fst(lt,js,1:my,:),my*3, MPI_DOUBLE_PRECISION, nranksend(irecv,js,isend), isend,  &
		               MPI_COMM_WORLD,status,ierror )
         enddo
         endif
         enddo
         enddo
       endif
!!mpi
    2 continue 
           
      if(mb_nrk(nrank).gt.0) then      

      do jy=1,my
      do js=mpsa-2,mpsa-4,-1
      do jt=itbmin(nrank),itbmax(nrank)
      cx1st=fst(jt,js,jy,1)
      cz1st=fst(jt,js,jy,3)

      fst(jt,js,jy,1)=cx1st*dcos(thst(jt))+cz1st*dsin(thst(jt))
      fst(jt,js,jy,3)=-cx1st*dsin(thst(jt))+cz1st*dcos(thst(jt))
       
!      vrst=fst(jt,js,jy,3)*dcos(thst(jt))+fst(jt,js,jy,5)*dsin(thst(jt))
!      vpst=-fst(jt,js,jy,3)*dsin(thst(jt))+fst(jt,js,jy,5)*dcos(thst(jt))
!      brst=fst(jt,js,jy,6)*dcos(thst(jt))+fst(jt,js,jy,8)*dsin(thst(jt))
!      bpst=-fst(jt,js,jy,6)*dsin(thst(jt))+fst(jt,js,jy,8)*dcos(thst(jt))
!      fst(jt,js,jy,3)=vrst
!      fst(jt,js,jy,5)=vpst
!      fst(jt,js,jy,6)=brst
!      fst(jt,js,jy,8)=bpst
      enddo
      enddo
      enddo    

      do 10 js=mpsa-2,mpsa-4,-1
!      call smth_ps_nrk(fst(itbmin(nrank):itbmax(nrank),js,:,:),itbmax(nrank)-itbmin(nrank)+1,3)
!ws_smps--------------------------
      do 11 k=1,kk
      do m=1,3
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      wst(jt)=fst(jt,js,jy,m)/2.+((fst(jt+1,js,jy,m)+fst(jt-1,js,jy,m))*3./16+(fst(jt+2,js,jy,m)+fst(jt-2,js,jy,m))/16.&
             +(fst(jt,js,jy+1,m)+fst(jt,js,jy-1,m))*3./16+(fst(jt,js,jy+2,m)+fst(jt,js,jy-2,m))/16.)/2.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,js,jy,m)=wst(jt)
      enddo

      enddo
      enddo

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         CALL MPI_Send(fst(ltmin:ltmax,js,1:my,:), 2*my*3, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,js,1:my,:), 2*my*3, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         CALL MPI_Send(fst(ltmin:ltmax,js,1:my,:), 2*my*3, MPI_DOUBLE_PRECISION, nrkb1(1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,js,1:my,:), 2*my*3, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         CALL MPI_Send(fst(ltmin:ltmax,js,1:my,:), 2*my*3, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 2,  &
		               MPI_COMM_WORLD,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,js,1:my,:), 2*my*3, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         CALL MPI_Send(fst(ltmin:ltmax,js,1:my,:), 2*my*3, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,js,1:my,:), 2*my*3, MPI_DOUBLE_PRECISION, nrkb1(1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   11 continue
!ws_smps--------------------
   10 continue

      do 61 m=1,3

      if(m.eq.1) then !!(m=1) ibnd==0
      is=mpsa 
      do 51 jy=1,my  
      do jt=itbmin(nrank),itbmax(nrank)
       fst(jt,is,jy,m)=0.           
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   51 continue

      is=mpsa-1 
      do 52 jy=1,my  
      do jt=itbmin(nrank),itbmax(nrank)
        call interp1d3l(fst(jt,is-3,jy,m),fst(jt,is-2,jy,m),fst(jt,is-1,jy,m),fst(jt,is+1,jy,m), &
                        ps(is-3),ps(is-2),ps(is-1),ps(is+1),ps(is),fst(jt,is,jy,m))                   
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   52 continue

      do 53 k=1,kk
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      wst(jt)=fst(jt,is,jy,m)/2.+((fst(jt+1,is,jy,m)+fst(jt-1,is,jy,m))*3./16+(fst(jt+2,is,jy,m)+fst(jt-2,is,jy,m))/16.&
             +(fst(jt,is,jy+1,m)+fst(jt,is,jy-1,m))*3./16+(fst(jt,is,jy+2,m)+fst(jt,is,jy-2,m))/16.)/2.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,is,jy,m)=wst(jt)
      enddo

      enddo

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 2,  &
		               MPI_COMM_WORLD,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   53 continue

      else  !!(m=2,3) bnd==1  
      is=mpsa-1    
      do 31 jy=1,my      
      do jt=itbmin(nrank),itbmax(nrank)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))  
!       fst(jt,is)=fst(jt,is-1)                   
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   31 continue 
!      call smth_ps_nrk(fst(itbmin(nrank):itbmax(nrank),is,:,:),itbmax(nrank)-itbmin(nrank)+1,3)
      do 32 k=1,kk
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      wst(jt)=fst(jt,is,jy,m)/2.+((fst(jt+1,is,jy,m)+fst(jt-1,is,jy,m))*3./16+(fst(jt+2,is,jy,m)+fst(jt-2,is,jy,m))/16.&
             +(fst(jt,is,jy+1,m)+fst(jt,is,jy-1,m))*3./16+(fst(jt,is,jy+2,m)+fst(jt,is,jy-2,m))/16.)/2.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,is,jy,m)=wst(jt)
      enddo

      enddo
      
         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 2,  &
		               MPI_COMM_WORLD,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   32 continue
!ws_smps--------------------
      is=mpsa
 
      do 41 jy=1,my  
      do jt=itbmin(nrank),itbmax(nrank)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))
!       fst(jt,is)=fst(jt,is-1)                   
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   41 continue
!      call smth_ps_nrk(fst(itbmin(nrank):itbmax(nrank),is,:,:),itbmax(nrank)-itbmin(nrank)+1,3)
!ws_smps--------------------
      do 42 k=1,kk
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      wst(jt)=fst(jt,is,jy,m)/2.+((fst(jt+1,is,jy,m)+fst(jt-1,is,jy,m))*3./16+(fst(jt+2,is,jy,m)+fst(jt-2,is,jy,m))/16.&
             +(fst(jt,is,jy+1,m)+fst(jt,is,jy-1,m))*3./16+(fst(jt,is,jy+2,m)+fst(jt,is,jy-2,m))/16.)/2.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,is,jy,m)=wst(jt)
      enddo

      enddo

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(1), 1,  &
		               MPI_COMM_WORLD,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 1,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)-1), 2,  &
		               MPI_COMM_WORLD,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(inrkb(nrank)+1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         CALL MPI_Send(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(mrkb), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         CALL MPI_Recv(fst(ltmin:ltmax,is,1:my,m), 2*my, MPI_DOUBLE_PRECISION, nrkb1(1), 2,  &
		               MPI_COMM_WORLD,status,ierror )
         endif
   42 continue
!ws_smps--------------------
      endif !!(m=2,3)ibnd==0
   61 continue

      do 3 jy=1,my           
      do ikb=1,mb_nrk(nrank)
        i=itb_nrk(ikb,nrank)
        ib=ib_nrk(ikb,nrank) 
        do m=1,3 
        do js=mpsa,mpsa-3,-1     
        call interp1d3l(fst(i-2,js,jy,m),fst(i-1,js,jy,m),fst(i,js,jy,m),fst(i+1,js,jy,m), &
                        thst(i-2),thst(i-1),thst(i),thst(i+1),thb(ib),f1s(ikb,js))
        enddo
        js=mpsa
        call interp1d3l(f1s(ikb,js),f1s(ikb,js-1),f1s(ikb,js-2),f1s(ikb,js-3), &
                        ps(js),ps(js-1),ps(js-2),ps(js-3),psb(ib),fsxz(ikb,m))
        enddo
!        fxz(jbx(ib),jbz(ib))=(asm_b(ib)*x1s(ib,js-1)+bsm_b(ib)*x1s(ib,js-2)+csm_b(ib)*x1s(ib,js-3)) &
!                 /(asm_b(ib)+bsm_b(ib)+csm_b(ib))

        cur(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,1)=fsxz(ikb,1)*dcos(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank))) &
                                                       -fsxz(ikb,3)*dsin(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank)))
        cur(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,2)=fsxz(ikb,2)
        cur(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,3)=fsxz(ikb,1)*dsin(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank))) &
                                                       +fsxz(ikb,3)*dcos(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank)))

      enddo
    3 continue 
      endif !!mb_nrk(nrank).gt.0
           
!      if(ibnd==2) then
!      
!      do ika1=1,ma1_nrk(nrank)
!        i=ita1_nrk(ika1,nrank)
!        ia1=ia1_nrk(ika1,nrank)
!        do js=mpsa,mpsa-3,1
!        call interp1d3l(fst(i-2,js),fst(i-1,js),fst(i,js),fst(i+1,js), &
!                        thst(i-2),thst(i-1),thst(i),thst(i+1),tha1(ia1),x1s(ika1,js,jy))
!        enddo
!        js=mpsa
!        call interp1d3l(x1s(ika1,js,jy),x1s(ika1,js-1,jy),x1s(ika1,js-2,jy),x1s(ika1,js-3,jy), &
!                        ps(js),ps(js-1),ps(js-2),ps(js-3),psa1(ia1),fxz(jxa1_nrk(ika1,nrank),jza1_nrk(ika1,nrank)))
!!        fxz(jxa1(ia1),jza1(ia1))=(asm_a1(ia1)*x1s(ia1,js-1)+bsm_a1(ia1)*x1s(ia1,js-2)+csm_a1(ia1)*x1s(ia1,js-3)) &
!!                  /(asm_a1(ia1)+bsm_a1(ia1)+csm_a1(ia1)) 
!      enddo
!      endif
      return
      end
!ws************************************************************
!ws*****************************************************************

     subroutine smthp_traceline_v1(kk)
      USE DECLARE
      include 'mpif.h'
      real*8,dimension(3) :: bbb
      real*8,dimension(2) :: ppp,sl
      real*8,dimension(mx,mz,my) :: wsm
      real*8 sl0,xxp,zzp,yyp,dyyp,dy1p,dx1p,dz1p,rxp,rxm,rzp,rzm,ryp,rym
      real*8 volume,volppp,volmpp,volppm,volmpm,volpmp,volmmp,volpmm,volmmm
      integer kk,k,isgn,jjx,jjz,jjy,im,j_add
      
      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)
!      volume=dyy*dxx*dzz
!
!      dyyp=yy(jy+1)-yy(jy)
!      dyym=yy(jy-1)-yy(jy)

!      dy1p=dyyp/nn
!      dy1m=dyym/nn
      sl0=0
      do 3 k=1,kk

      do 2 jy=1,my
      do 2 jz=iz_first+2,iz_last-2
      do 2 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
      sl(:)=0

!      jjy=jy
!      jjx=jx
!      jjz=jz

      do 1 isgn=1,2

      dyyp=yy(jy+(-1)**isgn)-yy(jy)

      dy1p=dyyp/nline    
      
100   continue
      xxp=xx(jx)
      zzp=zz(jz)
      yyp=yy(jy)
      bbb(:)=x(jx,jz,jy,6:8)

      do nc=1,mcycline
      dx1p=xxp*dy1p*bbb(1)/bbb(2)
      dz1p=xxp*dy1p*bbb(3)/bbb(2)    			
	!计算前进一步后的新位置
      xxp=xxp+dx1p
	  zzp=zzp+dz1p
      yyp=yyp+dy1p 
      if(yyp .gt. 2*pi) yyp=yyp-2*pi
      if(yyp .lt. 0) yyp=yyp+2*pi

      if( (xxp.gt.xx(jx+1)) .or. (xxp.lt.xx(jx-1)) &
      .or.(zzp.gt.zz(jz+1)) .or. (zzp.lt.zz(jz-1))) then
      if(nc.eq.1) then
      dy1p=dy1p*0.5
      goto 100
      else
      goto 200
      endif
      endif

      sl(isgn)=sl(isgn)+(-1)**isgn*sqrt(dx1p**2+dz1p**2+xxp**2*dy1p**2)
      
      jjx=jx+j_add(xxp,xx,mx,jx)
	  jjz=jz+j_add(zzp,zz,mz,jz)
	  jjy=jy+j_add(yyp,yy,my,jy)
      if(jjy.lt.1) jjy=jjy+my
      if(jjy.gt.my)jjy=jjy-my

      volume=(yy(jjy+1)-yy(jjy))*(xx(jjx+1)-xx(jjx))*(zz(jjz+1)-zz(jjz))

      rxp=xx(jjx+1)-xxp
      rxm=xxp-xx(jjx)
      rzp=zz(jjz+1)-zzp
      rzm=zzp-zz(jjz)
      ryp=yy(jjy+1)-yyp
      rym=yyp-yy(jjy)

      volppp=rxp*rzp*ryp
      volppm=rxp*rzp*rym
      volpmp=rxp*rzm*ryp
      volpmm=rxp*rzm*rym

      volmpp=rxm*rzp*ryp
      volmpm=rxm*rzp*rym
      volmmp=rxm*rzm*ryp
      volmmm=rxm*rzm*rym

      do im=6,8
      bbb(im-5)=(volppp*x(jjx,jjz,jjy,im)       +volmpp*x(jjx+1,jjz,jjy,im) &
	           +volppm*x(jjx,jjz,jjy+1,im)  +volmpm*x(jjx+1,jjz,jjy+1,im) &
	           +volpmp*x(jjx,jjz+1,jjy,im)     +volmmp*x(jjx+1,jjz+1,jjy,im)&
	           +volpmm*x(jjx,jjz+1,jjy+1,im)+volmmm*x(jjx+1,jjz+1,jjy+1,im))/volume
      enddo

      enddo

200   continue
      ppp(isgn)=(volppp*x1(jjx,jjz,jjy,2)       +volmpp*x1(jjx+1,jjz,jjy,2) &
	            +volppm*x1(jjx,jjz,jjy+1,2)  +volmpm*x1(jjx+1,jjz,jjy+1,2) &
	            +volpmp*x1(jjx,jjz+1,jjy,2)     +volmmp*x1(jjx+1,jjz+1,jjy,2)&
	            +volpmm*x1(jjx,jjz+1,jjy+1,2)+volmmm*x1(jjx+1,jjz+1,jjy+1,2))/volume         
                
1     continue

      wsm(jx,jz,jy)=csmp0*difc(ppp(1),x1(jx,jz,jy,2),ppp(2),sl(1),sl0,sl(2))
      endif
2     continue  


      x1(:,:,:,2)=x1(:,:,:,2)+wsm(:,:,:)

      call mpi_transfersm(x1(:,:,:,2),1)

3     continue 
      
      return
      end


!ws*****************************************************************

     subroutine smthp_traceline_v2(kk)
      USE DECLARE
      include 'mpif.h'
      real*8,dimension(2,mx,mz,my) :: sl,ppp,volppp,volmpp,volppm,volmpm,volpmp,volmmp,volpmm,volmmm
      real*8,dimension(mx,mz,my) :: wsm
      real*8,dimension(3) :: bbb
      integer,dimension(2,mx) :: jjx
      integer,dimension(2,mz) :: jjz
      integer,dimension(2,my) :: jjy
      real*8 sl0,xxp,zzp,yyp,dyyp,dy1p,dx1p,dz1p,rxp,rxm,rzp,rzm,ryp,rym,volume
      integer kk,isgn,sgn

      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)
!      volume=dyy*dxx*dzz
!
!      dyyp=yy(jy+1)-yy(jy)
!      dyym=yy(jy-1)-yy(jy)

!      dy1p=dyyp/nn
!      dy1m=dyym/nn
      sl0=0
!      volume=dyy*dxx*dzz

      do 2 jy=1,my
      do 2 jz=iz_first+2,iz_last-2
      do 2 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
      sl(:,jx,jz,jy)=0

      do 1 isgn=1,2

      dyyp=yy(jy+(-1)**isgn)-yy(jy)
!      
      dy1p=dyyp/nline 

100   continue

      xxp=xx(jx)
      zzp=zz(jz)
      yyp=yy(jy)
      
      bbb(:)=x(jx,jz,jy,6:8)
      
      do nc=1,mcycline
      dx1p=xxp*dy1p*bbb(1)/bbb(2)
      dz1p=xxp*dy1p*bbb(3)/bbb(2)

    			
	!计算前进一步后的新位置
      xxp=xxp+dx1p
	  zzp=zzp+dz1p
      yyp=yyp+dy1p 
      if(yyp.gt. 2*pi) yyp=yyp-2*pi
      if(yyp.lt. 0) yyp=yyp+2*pi

      if( (xxp.gt.xx(jx+1)) .or. (xxp.lt.xx(jx-1)) &
      .or.(zzp.gt.zz(jz+1)) .or. (zzp.lt.zz(jz-1))) then
      if(nc.eq.1) then
      dy1p=dy1p*0.5
      goto 100
      else
      goto 200
      endif
      endif

      sl(isgn,jx,jz,jy)=sl(isgn,jx,jz,jy)+(-1)**isgn*sqrt(dx1p**2+dz1p**2+xxp**2*dy1p**2)
      
      jjx(isgn,jx)=jx+floor((xxp-xx(jx))/abs(xx(jx)-xx(jx+sgn(xxp-xx(jx))) ))
	  jjz(isgn,jz)=jz+floor((zzp-zz(jz))/abs(zz(jz)-zz(jz+sgn(zzp-zz(jz))) ))   
	  jjy(isgn,jy)=jy+floor((yyp-yy(jy))/abs(yy(jy)-yy(jy+sgn(yyp-zz(jy))) ))
      if(jjy(isgn,jy).lt.1) jjy(isgn,jy)=jjy(isgn,jy)+my
      if(jjy(isgn,jy).gt.my)jjy(isgn,jy)=jjy(isgn,jy)-my

      rxp=xx(jjx(isgn,jx)+1)-xxp
      rxm=xxp-xx(jjx(isgn,jx))
      rzp=zz(jjz(isgn,jz)+1)-zzp
      rzm=zzp-zz(jjz(isgn,jz))
      ryp=yy(jjy(isgn,jy)+1)-yyp
      rym=yyp-yy(jjy(isgn,jy))
      volume=(xx(jjx(isgn,jx)+1)-xx(jjx(isgn,jx)))*(zz(jjz(isgn,jz)+1)-zz(jjz(isgn,jz)))*(yy(jjy(isgn,jy)+1)-yy(jjy(isgn,jy)))

      volppp(isgn,jx,jz,jy)=rxp*rzp*ryp/volume
      volppm(isgn,jx,jz,jy)=rxp*rzp*rym/volume
      volpmp(isgn,jx,jz,jy)=rxp*rzm*ryp/volume
      volpmm(isgn,jx,jz,jy)=rxp*rzm*rym/volume

      volmpp(isgn,jx,jz,jy)=rxm*rzp*ryp/volume
      volmpm(isgn,jx,jz,jy)=rxm*rzp*rym/volume
      volmmp(isgn,jx,jz,jy)=rxm*rzm*ryp/volume
      volmmm(isgn,jx,jz,jy)=rxm*rzm*rym/volume

      do m=6,8
      bbb(m-5)=(volppp(isgn,jx,jz,jy)*x(jjx(isgn,jx),jjz(isgn,jz),jjy(isgn,jy),m) &
               +volmpp(isgn,jx,jz,jy)*x(jjx(isgn,jx)+1,jjz(isgn,jz),jjy(isgn,jy),m) &
	           +volppm(isgn,jx,jz,jy)*x(jjx(isgn,jx),jjz(isgn,jz),jyp(jjy(isgn,jy)),m) &
               +volmpm(isgn,jx,jz,jy)*x(jjx(isgn,jx)+1,jjz(isgn,jz),jyp(jjy(isgn,jy)),m) &
	           +volpmp(isgn,jx,jz,jy)*x(jjx(isgn,jx),jjz(isgn,jz)+1,jjy(isgn,jy),m) &
               +volmmp(isgn,jx,jz,jy)*x(jjx(isgn,jx)+1,jjz(isgn,jz)+1,jjy(isgn,jy),m) &
	           +volpmm(isgn,jx,jz,jy)*x(jjx(isgn,jx),jjz(isgn,jz)+1,jyp(jjy(isgn,jy)),m) &
               +volmmm(isgn,jx,jz,jy)*x(jjx(isgn,jx)+1,jjz(isgn,jz)+1,jyp(jjy(isgn,jy)),m))

      enddo

      enddo

200   continue
                    
                
1     continue


      endif
2     continue  

      do 3 k=1,kk
      
      do 21 jy=1,my
      do 21 jz=iz_first+2,iz_last-2
      do 21 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
      ppp(isgn,jx,jz,jy)=(volppp(isgn,jx,jz,jy)*x(jjx(isgn,jx),jjz(isgn,jz),jjy(isgn,jy),2) &
               +volmpp(isgn,jx,jz,jy)*x(jjx(isgn,jx)+1,jjz(isgn,jz),jjy(isgn,jy),2) &
	           +volppm(isgn,jx,jz,jy)*x(jjx(isgn,jx),jjz(isgn,jz),jyp(jjy(isgn,jy)),2) &
               +volmpm(isgn,jx,jz,jy)*x(jjx(isgn,jx)+1,jjz(isgn,jz),jyp(jjy(isgn,jy)),2) &
	           +volpmp(isgn,jx,jz,jy)*x(jjx(isgn,jx),jjz(isgn,jz)+1,jjy(isgn,jy),2) &
               +volmmp(isgn,jx,jz,jy)*x(jjx(isgn,jx)+1,jjz(isgn,jz)+1,jjy(isgn,jy),2) &
	           +volpmm(isgn,jx,jz,jy)*x(jjx(isgn,jx),jjz(isgn,jz)+1,jyp(jjy(isgn,jy)),2) &
               +volmmm(isgn,jx,jz,jy)*x(jjx(isgn,jx)+1,jjz(isgn,jz)+1,jyp(jjy(isgn,jy)),2))

      wsm(jx,jz,jy)=csmp0*difc(ppp(1,jx,jz,jy),x1(jx,jz,jy,2),ppp(2,jx,jz,jy),sl(1,jx,jz,jy),sl0,sl(2,jx,jz,jy))
      endif
21    continue 

      x1(:,:,:,2)=x1(:,:,:,2)+wsm(:,:,:)

      call mpi_transfersm(x1(:,:,:,2),1)

3     continue 
      
      return
      end

!ws*****************************************************************

     subroutine tracelineat(jxl,jzl,jyl,lx,lz,ly,sl,vol)
      USE DECLARE
      include 'mpif.h'
      integer,dimension(2) :: lx,lz,ly
      real*8,dimension(2) :: sl
      real*8,dimension(8,2) :: vol
      real*8,dimension(3) :: bbb
      real*8 sl0,xxp,zzp,yyp,dyyp,dy1p,dx1p,dz1p,rxp,rxm,rzp,rzm,ryp,rym
      real*8 volume,volppp,volmpp,volppm,volmpm,volpmp,volmmp,volpmm,volmmm
      integer isgn,jjx,jjz,jjy,jxl,jzl,jyl,j_add
       
      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)
!      volume=dyy*dxx*dzz
!
!      dyyp=yy(jy+1)-yy(jy)
!      dyym=yy(jy-1)-yy(jy)

!      dy1p=dyyp/nn
!      dy1m=dyym/nn
      sl0=0
      sl(:)=0.

      do 1 isgn=1,2


      dyyp=yy(jyl+(-1)**isgn)-yy(jyl)
      dy1p=dyyp/nline    

100   continue   
      xxp=xx(jxl)
      zzp=zz(jzl)
      yyp=yy(jyl)
      
      bbb(:)=x(jxl,jzl,jyl,6:8)

      do nc=1,mcycline
      dx1p=xxp*dy1p*bbb(1)/bbb(2)
      dz1p=xxp*dy1p*bbb(3)/bbb(2)
    			
	!计算前进一步后的新位置
      xxp=xxp+dx1p
	  zzp=zzp+dz1p
      yyp=yyp+dy1p 
      if(yyp.gt. 2*pi) yyp=yyp-2*pi
      if(yyp.lt. 0) yyp=yyp+2*pi

!      if( (xxp.gt.xx(jxl+1)) .or. (xxp.lt.xx(jxl-1)) &
!      .or.(zzp.gt.zz(jzl+1)) .or. (zzp.lt.zz(jzl-1))) then
!      if(nc.eq.1) then
!      dy1p=dy1p*0.5
!      goto 100
!      else
!      goto 200
!      endif
!      endif

      sl(isgn)=sl(isgn)+(-1)**isgn*sqrt(dx1p**2+dz1p**2+xxp**2*dy1p**2)

      jjx=floor((xxp-xx(ix_first))/dxx)+ix_first
	  jjz=floor((zzp-zz(iz_first))/dzz)+iz_first   
	  jjy=floor((yyp-yy(iy_first))/dyy)+iy_first      
!      jjx=jxl+j_add(xxp,xx,mx,jxl)
!	  jjz=jzl+j_add(zzp,zz,mz,jzl)
!	  jjy=jyl+j_add(yyp,yy,my,jyl)
      if(jjy.lt.1) jjy=jjy+my
      if(jjy.gt.my)jjy=jjy-my
!      write(*,*) nrank,'nc',nc,jjy,xxp,jjx,zzp,jjz

      rxp=xx(jjx+1)-xxp
      rxm=xxp-xx(jjx)
      rzp=zz(jjz+1)-zzp
      rzm=zzp-zz(jjz)
      ryp=yy(jjy+1)-yyp
      rym=yyp-yy(jjy)

      volume=(xx(jjx+1)-xx(jjx))*(zz(jjz+1)-zz(jjz))*(yy(jjy+1)-yy(jjy))

      volppp=rxp*rzp*ryp/volume
      volppm=rxp*rzp*rym/volume
      volpmp=rxp*rzm*ryp/volume
      volpmm=rxp*rzm*rym/volume

      volmpp=rxm*rzp*ryp/volume
      volmpm=rxm*rzp*rym/volume
      volmmp=rxm*rzm*ryp/volume
      volmmm=rxm*rzm*rym/volume

      do m=6,8
      bbb(m-5)=(volppp*x(jjx,jjz,jjy,m)       +volmpp*x(jjx+1,jjz,jjy,m) &
	           +volppm*x(jjx,jjz,jjy+1,m)  +volmpm*x(jjx+1,jjz,jjy+1,m) &
	           +volpmp*x(jjx,jjz+1,jjy,m)     +volmmp*x(jjx+1,jjz+1,jjy,m)&
	           +volpmm*x(jjx,jjz+1,jjy+1,m)+volmmm*x(jjx+1,jjz+1,jjy+1,m))/volume
      enddo

      enddo

200   continue

      vol(1,isgn)=volppp
      vol(2,isgn)=volmpp
      vol(3,isgn)=volppm
      vol(4,isgn)=volmpm
      vol(5,isgn)=volpmp
      vol(6,isgn)=volmmp
      vol(7,isgn)=volpmm
      vol(8,isgn)=volmmm
      
      lx(isgn)=jjx
      lz(isgn)=jjz
      ly(isgn)=jjy     
                
1     continue

      return
      end
!ws********************************************
      subroutine smthp1_traceline(kk)
      USE DECLARE
      include 'mpif.h'
      integer,dimension(2) :: lx,lz,ly
      real*8,dimension(2) :: sl,ppp
      real*8,dimension(8,2) :: vol
      integer,dimension(2,mx,mz,my) :: llx,llz,lly
      real*8,dimension(2,mx,mz,my) :: ssl
      real*8,dimension(8,2,mx,mz,my) :: voll
      real*8,dimension(mx,mz,my) :: wsm
      integer kk,k,isgn,jjx,jjz,jjy
      real*8 ssl0

      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)

      ssl0=0.

      do 2 jy=1,my
      do 2 jz=iz_first+2,iz_last-2
      do 2 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
      call tracelineat(jx,jz,jy,lx,lz,ly,sl,vol)

      llx(:,jx,jz,jy)=lx(:)
      llz(:,jx,jz,jy)=lz(:)
      lly(:,jx,jz,jy)=ly(:)
      ssl(:,jx,jz,jy)=sl(:)
      voll(:,:,jx,jz,jy)=vol(:,:)
      endif
2     continue  

      do 3 k=1,kk
      
      do 21 jy=1,my
      do 21 jz=iz_first+2,iz_last-2
      do 21 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
      do isgn=1,2
      jjx=llx(isgn,jx,jz,jy)
      jjz=llz(isgn,jx,jz,jy)
      jjy=lly(isgn,jx,jz,jy)

      ppp(isgn)=voll(1,isgn,jx,jz,jy)*x1(jjx,jjz,jjy,2)  &
               +voll(2,isgn,jx,jz,jy)*x1(jjx+1,jjz,jjy,2) &
	           +voll(3,isgn,jx,jz,jy)*x1(jjx,jjz,jjy+1,2) &
               +voll(4,isgn,jx,jz,jy)*x1(jjx+1,jjz,jjy+1,2) &
	           +voll(5,isgn,jx,jz,jy)*x1(jjx,jjz+1,jjy,2) &
               +voll(6,isgn,jx,jz,jy)*x1(jjx+1,jjz+1,jjy,2) &
	           +voll(7,isgn,jx,jz,jy)*x1(jjx,jjz+1,jjy+1,2) &
               +voll(8,isgn,jx,jz,jy)*x1(jjx+1,jjz+1,jjy+1,2)
      enddo

      wsm(jx,jz,jy)=csmp0*difc(ppp(1),x1(jx,jz,jy,2),ppp(2),ssl(1,jx,jz,jy),ssl0,ssl(2,jx,jz,jy))
      endif
21    continue 
      x1(:,:,:,2)=x1(:,:,:,2)+wsm(:,:,:)

      call mpi_transfersm(x1(:,:,:,2),1)
!      do jy=1,my
!      x(:,:,jy,2)=x1(:,:,jy,2)+xint(:,:,2)
!      enddo
3     continue 
      
      return
      end

!ws********************************************
      subroutine smthp_traceline(kk)
      USE DECLARE
      include 'mpif.h'
      integer,dimension(2) :: lx,lz,ly
      real*8,dimension(2) :: sl,ppp
      real*8,dimension(8,2) :: vol
      integer,dimension(2,mx,mz,my) :: llx,llz,lly
      real*8,dimension(2,mx,mz,my) :: ssl
      real*8,dimension(8,2,mx,mz,my) :: voll
      real*8,dimension(mx,mz,my) :: wsm
      integer kk,k,isgn,jjx,jjz,jjy
      real*8 ssl0,csmp

      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)

      ssl0=0.

      do 2 jy=iy_first+2,iy_last-2
      do 2 jz=iz_first+2,iz_last-2
      do 2 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
!      if(psi(jx,jz).lt.psia1) then
      call tracelineat2(jx,jz,jy,lx,lz,ly,sl,vol)

      llx(:,jx,jz,jy)=lx(:)
      llz(:,jx,jz,jy)=lz(:)
      lly(:,jx,jz,jy)=ly(:)
      ssl(:,jx,jz,jy)=sl(:)
      voll(:,:,jx,jz,jy)=vol(:,:)
      endif
2     continue  

      do 3 k=1,kk
      
      do 21 jy=iy_first+2,iy_last-2
      do 21 jz=iz_first+2,iz_last-2
      do 21 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
!      if(psi(jx,jz).lt.psia1) then
      do isgn=1,2
      jjx=llx(isgn,jx,jz,jy)
      jjz=llz(isgn,jx,jz,jy)
      jjy=lly(isgn,jx,jz,jy)

      ppp(isgn)=voll(1,isgn,jx,jz,jy)*x(jjx,jjz,jjy,2)  &
               +voll(2,isgn,jx,jz,jy)*x(jjx+1,jjz,jjy,2) &
	       +voll(3,isgn,jx,jz,jy)*x(jjx,jjz,jjy+1,2) &
               +voll(4,isgn,jx,jz,jy)*x(jjx+1,jjz,jjy+1,2) &
	       +voll(5,isgn,jx,jz,jy)*x(jjx,jjz+1,jjy,2) &
               +voll(6,isgn,jx,jz,jy)*x(jjx+1,jjz+1,jjy,2) &
	       +voll(7,isgn,jx,jz,jy)*x(jjx,jjz+1,jjy+1,2) &
               +voll(8,isgn,jx,jz,jy)*x(jjx+1,jjz+1,jjy+1,2)
      enddo
      csmp=csmp0all*(1-tanh((rr(jx,jz)-0.9)*40))/2.
!      csmp=csmp0all*(1-tanh((psi(jx,jz)-pstrans)/pstransw))/2.
      wsm(jx,jz,jy)=csmp*difc(ppp(1),x(jx,jz,jy,2),ppp(2),ssl(1,jx,jz,jy),ssl0,ssl(2,jx,jz,jy))
      endif
21    continue 
      x(:,:,:,2)=x(:,:,:,2)+wsm(:,:,:)

      call mpi_transfersm(x(:,:,:,2),1)

3     continue 
      
      return
      end

!ws*****************************************************************

     subroutine tracelineat2(jxl,jzl,jyl,lx,lz,ly,sl,vol)
      USE DECLARE
      include 'mpif.h'
      integer,dimension(2) :: lx,lz,ly
      real*8,dimension(2) :: sl
      real*8,dimension(8,2) :: vol
      real*8,dimension(3) :: bbb
      real*8 sl0,xxp,zzp,yyp,dyyp,dy1p,dx1p,dz1p,rxp,rxm,rzp,rzm,ryp,rym
      real*8 xxp0,zzp0,yyp0,xxp12,zzp12,yyp12
      real*8 volume,volppp,volmpp,volppm,volmpm,volpmp,volmmp,volpmm,volmmm
      integer isgn,jjx,jjz,jjy,jxl,jzl,jyl,j_add
       
      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)
!      volume=dyy*dxx*dzz
!
!      dyyp=yy(jy+1)-yy(jy)
!      dyym=yy(jy-1)-yy(jy)

!      dy1p=dyyp/nn
!      dy1m=dyym/nn
      sl0=0
      sl(:)=0.

      do 1 isgn=1,2


      dyyp=yy(jyl+(-1)**isgn)-yy(jyl)
      dy1p=dyyp/nline    

100   continue   
      xxp0=xx(jxl)
      zzp0=zz(jzl)
      yyp0=yy(jyl)

      bbb(:)=x(jxl,jzl,jyl,6:8)

      do nc=1,mcycline
!      dy1p=dyyp/nline/xxp0*xzero
      dx1p=xxp0*dy1p*bbb(1)/bbb(2)
      dz1p=xxp0*dy1p*bbb(3)/bbb(2)
    			
	!计算前进一步后的新位置
      xxp=xxp0+dx1p
	  zzp=zzp0+dz1p
      yyp=yyp0+dy1p 
      if(yyp.gt. 2*pi) yyp=yyp-2*pi
      if(yyp.lt. 0) yyp=yyp+2*pi

!      if( (xxp.gt.xx(jxl+1)) .or. (xxp.lt.xx(jxl-1)) &
!      .or.(zzp.gt.zz(jzl+1)) .or. (zzp.lt.zz(jzl-1))) then
!      if(nc.eq.1) then
!      dy1p=dy1p*0.5
!      goto 100
!      else
!      goto 200
!      endif
!      endif

!      sl(isgn)=sl(isgn)+(-1)**isgn*sqrt(dx1p**2+dz1p**2+xxp**2*dy1p**2)
      jjx=floor((xxp-xx(ix_first))/dxx)+ix_first
	  jjz=floor((zzp-zz(iz_first))/dzz)+iz_first   
	  jjy=floor((yyp-yy(iy_first))/dyy)+iy_first      
!      jjx=jxl+j_add(xxp,xx,mx,jxl)
!	  jjz=jzl+j_add(zzp,zz,mz,jzl)
!	  jjy=jyl+j_add(yyp,yy,my,jyl)
      if(jjy.lt.1) jjy=jjy+myt
      if(jjy.gt.myt)jjy=jjy-myt
!      write(*,*) nrank,'nc',nc,jjy,xxp,jjx,zzp,jjz

      rxp=xx(jjx+1)-xxp
      rxm=xxp-xx(jjx)
      rzp=zz(jjz+1)-zzp
      rzm=zzp-zz(jjz)
      ryp=yy(jjy+1)-yyp
      rym=yyp-yy(jjy)

      volume=(xx(jjx+1)-xx(jjx))*(zz(jjz+1)-zz(jjz))*(yy(jjy+1)-yy(jjy))

      volppp=rxp*rzp*ryp/volume
      volppm=rxp*rzp*rym/volume
      volpmp=rxp*rzm*ryp/volume
      volpmm=rxp*rzm*rym/volume

      volmpp=rxm*rzp*ryp/volume
      volmpm=rxm*rzp*rym/volume
      volmmp=rxm*rzm*ryp/volume
      volmmm=rxm*rzm*rym/volume

      do m=6,8
      bbb(m-5)=(volppp*x(jjx,jjz,jjy,m)       +volmpp*x(jjx+1,jjz,jjy,m) &
	           +volppm*x(jjx,jjz,jjy+1,m)  +volmpm*x(jjx+1,jjz,jjy+1,m) &
	           +volpmp*x(jjx,jjz+1,jjy,m)     +volmmp*x(jjx+1,jjz+1,jjy,m)&
	           +volpmm*x(jjx,jjz+1,jjy+1,m)+volmmm*x(jjx+1,jjz+1,jjy+1,m))/volume
      enddo


      xxp12=xxp0+0.5*dx1p
	  zzp12=zzp0+0.5*dz1p
      yyp12=yyp0+0.5*dy1p

!      dy1p=dyyp/nline/xxp*xzero
      dx1p=xxp*dy1p*bbb(1)/bbb(2)
      dz1p=xxp*dy1p*bbb(3)/bbb(2)

      xxp=xxp12+0.5*dx1p
      zzp=zzp12+0.5*dz1p
      yyp=yyp12+0.5*dy1p
      if(yyp.gt. 2*pi) yyp=yyp-2*pi
      if(yyp.lt. 0) yyp=yyp+2*pi

      sl(isgn)=sl(isgn)+(-1)**isgn*sqrt(dx1p**2+dz1p**2+xxp**2*dy1p**2)
      
!      jjx=jxl+j_add(xxp,xx,mx,jxl)
!	  jjz=jzl+j_add(zzp,zz,mz,jzl)
!	  jjy=jyl+j_add(yyp,yy,my,jyl)

      jjx=floor((xxp-xx(ix_first))/dxx)+ix_first
	  jjz=floor((zzp-zz(iz_first))/dzz)+iz_first   
	  jjy=floor((yyp-yy(iy_first))/dyy)+iy_first
      if(jjy.lt.1) jjy=jjy+myt
      if(jjy.gt.myt)jjy=jjy-myt
!      write(*,*) nrank,'nc',nc,jjy,xxp,jjx,zzp,jjz

      rxp=xx(jjx+1)-xxp
      rxm=xxp-xx(jjx)
      rzp=zz(jjz+1)-zzp
      rzm=zzp-zz(jjz)
      ryp=yy(jjy+1)-yyp
      rym=yyp-yy(jjy)

      volume=(xx(jjx+1)-xx(jjx))*(zz(jjz+1)-zz(jjz))*(yy(jjy+1)-yy(jjy))

      volppp=rxp*rzp*ryp/volume
      volppm=rxp*rzp*rym/volume
      volpmp=rxp*rzm*ryp/volume
      volpmm=rxp*rzm*rym/volume

      volmpp=rxm*rzp*ryp/volume
      volmpm=rxm*rzp*rym/volume
      volmmp=rxm*rzm*ryp/volume
      volmmm=rxm*rzm*rym/volume

      do m=6,8
      bbb(m-5)=(volppp*x(jjx,jjz,jjy,m)       +volmpp*x(jjx+1,jjz,jjy,m) &
	           +volppm*x(jjx,jjz,jjy+1,m)  +volmpm*x(jjx+1,jjz,jjy+1,m) &
	           +volpmp*x(jjx,jjz+1,jjy,m)     +volmmp*x(jjx+1,jjz+1,jjy,m)&
	           +volpmm*x(jjx,jjz+1,jjy+1,m)+volmmm*x(jjx+1,jjz+1,jjy+1,m))/volume
      enddo

      xxp0=xxp
      zzp0=zzp
      yyp0=yyp

      enddo

200   continue

      vol(1,isgn)=volppp
      vol(2,isgn)=volmpp
      vol(3,isgn)=volppm
      vol(4,isgn)=volmpm
      vol(5,isgn)=volpmp
      vol(6,isgn)=volmmp
      vol(7,isgn)=volpmm
      vol(8,isgn)=volmmm
      
      lx(isgn)=jjx
      lz(isgn)=jjz
      ly(isgn)=jjy     
                
1     continue

      return
      end

!ws*****************************************************************

!     subroutine tracelineat2_dx(jxl,jzl,jyl,lx,lz,ly,sl,vol)
!      USE DECLARE
!      include 'mpif.h'
!      integer,dimension(2) :: lx,lz,ly
!      real*8,dimension(2) :: sl
!      real*8,dimension(8,2) :: vol
!      real*8,dimension(3) :: bbb
!      real*8 sl0,xxp,zzp,yyp,dyyp,dy1p,dx1p,dz1p,rxp,rxm,rzp,rzm,ryp,rym
!      real*8 xxp0,zzp0,yyp0,xxp12,zzp12,yyp12
!      real*8 volume,volppp,volmpp,volppm,volmpm,volpmp,volmmp,volpmm,volmmm
!      integer isgn,jjx,jjz,jjy,jxl,jzl,jyl
!       
!      difc(fm1,f0,fp1,xm1,x0,xp1)= &
!       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)
!!      volume=dyy*dxx*dzz
!!
!!      dyyp=yy(jy+1)-yy(jy)
!!      dyym=yy(jy-1)-yy(jy)
!
!!      dy1p=dyyp/nn
!!      dy1m=dyym/nn
!      sl0=0
!      sl(:)=0.
!
!      do 1 isgn=1,2
!
!
!!      dyyp=yy(jyl+(-1)**isgn)-yy(jyl)
!!      dy1p=dyyp/nline
!      
!      dx1p=xx(jxl+(-1)**isgn)-xx(jxl)    
!
!100   continue   
!      xxp0=xx(jxl)
!      zzp0=zz(jzl)
!      yyp0=yy(jyl)
!
!      bbb(:)=x(jxl,jzl,jyl,6:8)
!
!      do nc=1,mcycline
!!      dy1p=dyyp/nline/xxp0*xzero
!!      dx1p=xxp0*dy1p*bbb(1)/bbb(2)
!!      dz1p=xxp0*dy1p*bbb(3)/bbb(2)
!       if(abs(bbb(1)) .lt. 1.e-7)
!       dz1p=dx1p*bbb(3)/bbb(1)
!       dy1p=dx1p*bbb(2)/bbb(1)/xxp0		
!	!计算前进一步后的新位置
!      xxp=xxp0+dx1p
!	  zzp=zzp0+dz1p
!      yyp=yyp0+dy1p 
!      if(yyp.gt. 2*pi) yyp=yyp-2*pi
!      if(yyp.lt. 0) yyp=yyp+2*pi
!
!!      if( (xxp.gt.xx(jxl+1)) .or. (xxp.lt.xx(jxl-1)) &
!!      .or.(zzp.gt.zz(jzl+1)) .or. (zzp.lt.zz(jzl-1))) then
!!      if(nc.eq.1) then
!!      dy1p=dy1p*0.5
!!      goto 100
!!      else
!!      goto 200
!!      endif
!!      endif
!
!!      sl(isgn)=sl(isgn)+(-1)**isgn*sqrt(dx1p**2+dz1p**2+xxp**2*dy1p**2)
!      
!      jjx=floor((xxp-xx(ix_first))/dxx)+ix_first
!	  jjz=floor((zzp-zz(iz_first))/dzz)+iz_first   
!	  jjy=floor((yyp-yy(1))/dyy)+1
!      if(jjy.lt.1) jjy=jjy+my
!      if(jjy.gt.my)jjy=jjy-my
!!      write(*,*) nrank,'nc',nc,jjy,xxp,jjx,zzp,jjz
!
!      rxp=xx(jjx+1)-xxp
!      rxm=xxp-xx(jjx)
!      rzp=zz(jjz+1)-zzp
!      rzm=zzp-zz(jjz)
!      ryp=yy(jjy+1)-yyp
!      rym=yyp-yy(jjy)
!
!      volume=(xx(jjx+1)-xx(jjx))*(zz(jjz+1)-zz(jjz))*(yy(jjy+1)-yy(jjy))
!
!      volppp=rxp*rzp*ryp/volume
!      volppm=rxp*rzp*rym/volume
!      volpmp=rxp*rzm*ryp/volume
!      volpmm=rxp*rzm*rym/volume
!
!      volmpp=rxm*rzp*ryp/volume
!      volmpm=rxm*rzp*rym/volume
!      volmmp=rxm*rzm*ryp/volume
!      volmmm=rxm*rzm*rym/volume
!
!      do m=6,8
!      bbb(m-5)=(volppp*x(jjx,jjz,jjy,m)       +volmpp*x(jjx+1,jjz,jjy,m) &
!	           +volppm*x(jjx,jjz,jjy+1,m)  +volmpm*x(jjx+1,jjz,jjy+1,m) &
!	           +volpmp*x(jjx,jjz+1,jjy,m)     +volmmp*x(jjx+1,jjz+1,jjy,m)&
!	           +volpmm*x(jjx,jjz+1,jjy+1,m)+volmmm*x(jjx+1,jjz+1,jjy+1,m))/volume
!      enddo
!
!
!      xxp12=xxp0+0.5*dx1p
!	  zzp12=zzp0+0.5*dz1p
!      yyp12=yyp0+0.5*dy1p
!
!!      dy1p=dyyp/nline/xxp*xzero
!!      dx1p=xxp*dy1p*bbb(1)/bbb(2)
!!      dz1p=xxp*dy1p*bbb(3)/bbb(2)
!      dz1p=dx1p*bbb(3)/bbb(1)
!      dy1p=dx1p*bbb(2)/bbb(1)/xxp	
!      
!      xxp=xxp12+0.5*dx1p
!      zzp=zzp12+0.5*dz1p
!      yyp=yyp12+0.5*dy1p
!      if(yyp.gt. 2*pi) yyp=yyp-2*pi
!      if(yyp.lt. 0) yyp=yyp+2*pi
!
!      sl(isgn)=sl(isgn)+(-1)**isgn*sqrt(dx1p**2+dz1p**2+xxp**2*dy1p**2)
!      
!      jjx=floor((xxp-xx(ix_first))/dxx)+ix_first
!	  jjz=floor((zzp-zz(iz_first))/dzz)+iz_first   
!	  jjy=floor((yyp-yy(1))/dyy)+1
!      if(jjy.lt.1) jjy=jjy+my
!      if(jjy.gt.my)jjy=jjy-my
!!      write(*,*) nrank,'nc',nc,jjy,xxp,jjx,zzp,jjz
!
!      rxp=xx(jjx+1)-xxp
!      rxm=xxp-xx(jjx)
!      rzp=zz(jjz+1)-zzp
!      rzm=zzp-zz(jjz)
!      ryp=yy(jjy+1)-yyp
!      rym=yyp-yy(jjy)
!
!      volume=(xx(jjx+1)-xx(jjx))*(zz(jjz+1)-zz(jjz))*(yy(jjy+1)-yy(jjy))
!
!      volppp=rxp*rzp*ryp/volume
!      volppm=rxp*rzp*rym/volume
!      volpmp=rxp*rzm*ryp/volume
!      volpmm=rxp*rzm*rym/volume
!
!      volmpp=rxm*rzp*ryp/volume
!      volmpm=rxm*rzp*rym/volume
!      volmmp=rxm*rzm*ryp/volume
!      volmmm=rxm*rzm*rym/volume
!
!      do m=6,8
!      bbb(m-5)=(volppp*x(jjx,jjz,jjy,m)       +volmpp*x(jjx+1,jjz,jjy,m) &
!	           +volppm*x(jjx,jjz,jjy+1,m)  +volmpm*x(jjx+1,jjz,jjy+1,m) &
!	           +volpmp*x(jjx,jjz+1,jjy,m)     +volmmp*x(jjx+1,jjz+1,jjy,m)&
!	           +volpmm*x(jjx,jjz+1,jjy+1,m)+volmmm*x(jjx+1,jjz+1,jjy+1,m))/volume
!      enddo
!
!      xxp0=xxp
!      zzp0=zzp
!      yyp0=yyp
!
!      enddo
!
!200   continue
!
!      vol(1,isgn)=volppp
!      vol(2,isgn)=volmpp
!      vol(3,isgn)=volppm
!      vol(4,isgn)=volmpm
!      vol(5,isgn)=volpmp
!      vol(6,isgn)=volmmp
!      vol(7,isgn)=volpmm
!      vol(8,isgn)=volmmm
!      
!      lx(isgn)=jjx
!      lz(isgn)=jjz
!      ly(isgn)=jjy     
!                
!1     continue
!
!      return
!      end
!ws********************************************
      subroutine smthp2_traceline(kk)
      USE DECLARE
      include 'mpif.h'
      integer,dimension(2) :: lx,lz,ly
      real*8,dimension(2) :: sl,ppp,ppb
      real*8,dimension(8,2) :: vol
      integer,dimension(2,mx,mz,my) :: llx,llz,lly
      real*8,dimension(2,mx,mz,my) :: ssl
      real*8,dimension(8,2,mx,mz,my) :: voll
      real*8,dimension(mx,mz,my) :: wsm,pb !,wsm1,wsm2
      integer kk,k,isgn,jjx,jjz,jjy
      real*8 ssl0,bb

      !  dif1= (d f / dx)*dx2 = [((x0-xm1)/(xp1-x0)*(fp1-f0)-(xp1-x0)/(x0-xm1)*(f0-fm1))/(xp1-xm1)] *(xp1-x0)*(x0-xm1) 
      dif1(fm1,f0,fp1,xm1,x0,xp1)= &
       ((x0-xm1)**2*(fp1-f0)-(xp1-x0)**2*(f0-fm1))/(xp1-xm1)


      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)

      ssl0=0.

      do 2 jy=1,my
      do 2 jz=iz_first,iz_last
      do 2 jx=ix_first,ix_last
      if(psi(jx,jz).lt.ps(mpsa-2)) then
      call tracelineat2(jx,jz,jy,lx,lz,ly,sl,vol)

      llx(:,jx,jz,jy)=lx(:)
      llz(:,jx,jz,jy)=lz(:)
      lly(:,jx,jz,jy)=ly(:)
      ssl(:,jx,jz,jy)=sl(:)
      voll(:,:,jx,jz,jy)=vol(:,:)
      endif      
      bb=sqrt(x(jx,jz,jy,6)**2+x(jx,jz,jy,7)**2+x(jx,jz,jy,8)**2)
      pb(jx,jz,jy)=(x1(jx,jz,jy,6)*xint_dx(jx,jz,2)+x1(jx,jz,jy,8)*xint_dz(jx,jz,2))/bb
2     continue  

      do 3 k=1,kk
      
      do 21 jy=1,my
      do 21 jz=iz_first+2,iz_last-2
      do 21 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
      do isgn=1,2
      jjx=llx(isgn,jx,jz,jy)
      jjz=llz(isgn,jx,jz,jy)
      jjy=lly(isgn,jx,jz,jy)

      ppp(isgn)=voll(1,isgn,jx,jz,jy)*x1(jjx,jjz,jjy,2)  &
               +voll(2,isgn,jx,jz,jy)*x1(jjx+1,jjz,jjy,2) &
	           +voll(3,isgn,jx,jz,jy)*x1(jjx,jjz,jjy+1,2) &
               +voll(4,isgn,jx,jz,jy)*x1(jjx+1,jjz,jjy+1,2) &
	           +voll(5,isgn,jx,jz,jy)*x1(jjx,jjz+1,jjy,2) &
               +voll(6,isgn,jx,jz,jy)*x1(jjx+1,jjz+1,jjy,2) &
	           +voll(7,isgn,jx,jz,jy)*x1(jjx,jjz+1,jjy+1,2) &
               +voll(8,isgn,jx,jz,jy)*x1(jjx+1,jjz+1,jjy+1,2)

      ppb(isgn)=voll(1,isgn,jx,jz,jy)*pb(jjx,jjz,jjy)  &
               +voll(2,isgn,jx,jz,jy)*pb(jjx+1,jjz,jjy) &
	           +voll(3,isgn,jx,jz,jy)*pb(jjx,jjz,jjy+1) &
               +voll(4,isgn,jx,jz,jy)*pb(jjx+1,jjz,jjy+1) &
	           +voll(5,isgn,jx,jz,jy)*pb(jjx,jjz+1,jjy) &
               +voll(6,isgn,jx,jz,jy)*pb(jjx+1,jjz+1,jjy) &
	           +voll(7,isgn,jx,jz,jy)*pb(jjx,jjz+1,jjy+1) &
               +voll(8,isgn,jx,jz,jy)*pb(jjx+1,jjz+1,jjy+1)

      enddo
!       wsm1(jx,jz,jy)=difc(ppp(1),x1(jx,jz,jy,2),ppp(2),ssl(1,jx,jz,jy),ssl0,ssl(2,jx,jz,jy))
!       wsm2(jx,jz,jy)=dif1(ppb(1),pb(jx,jz,jy),ppb(2),ssl(1,jx,jz,jy),ssl0,ssl(2,jx,jz,jy))
!	wsm(jx,jz,jy)=csmp0all*(wsm1(jx,jz,jy)+wsm2(jx,jz,jy))
      wsm(jx,jz,jy)=csmp0all*(difc(ppp(1),x1(jx,jz,jy,2),ppp(2),ssl(1,jx,jz,jy),ssl0,ssl(2,jx,jz,jy)) &
                             +dif1(ppb(1),pb(jx,jz,jy),ppb(2),ssl(1,jx,jz,jy),ssl0,ssl(2,jx,jz,jy)))
      endif
21    continue 
      x(:,:,:,2)=x(:,:,:,2)+wsm(:,:,:)

      call mpi_transfersm(x(:,:,:,2),1)

3     continue 
      
      return
      end

!ws********************************************
      subroutine smthp_traceline_5p(kk)
      USE DECLARE
      include 'mpif.h'
      integer,dimension(4) :: lx,lz,ly
      real*8,dimension(4) :: sl,ppp
      real*8,dimension(8,4) :: vol
      integer,dimension(4,mx,mz,my) :: llx,llz,lly
      real*8,dimension(4,mx,mz,my) :: ssl
      real*8,dimension(8,4,mx,mz,my) :: voll
      real*8,dimension(mx,mz,my) :: wsm
      integer kk,k,isgn,jjx,jjz,jjy
      real*8 ssl0,csmp

      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)

      ssl0=0.

      do 2 jy=1,my
      do 2 jz=iz_first+2,iz_last-2
      do 2 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
!      if(psi(jx,jz).lt.psia1) then
      call tracelineat_5p(jx,jz,jy,lx,lz,ly,sl,vol)

      llx(:,jx,jz,jy)=lx(:)
      llz(:,jx,jz,jy)=lz(:)
      lly(:,jx,jz,jy)=ly(:)
      ssl(:,jx,jz,jy)=sl(:)
      voll(:,:,jx,jz,jy)=vol(:,:)
      endif
2     continue  

      do 3 k=1,kk
      
      do 21 jy=1,my
      do 21 jz=iz_first+2,iz_last-2
      do 21 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
!      if(psi(jx,jz).lt.psia1) then
      do isgn=1,4
      jjx=llx(isgn,jx,jz,jy)
      jjz=llz(isgn,jx,jz,jy)
      jjy=lly(isgn,jx,jz,jy)

      ppp(isgn)=voll(1,isgn,jx,jz,jy)*x(jjx,jjz,jjy,2)  &
               +voll(2,isgn,jx,jz,jy)*x(jjx+1,jjz,jjy,2) &
	       +voll(3,isgn,jx,jz,jy)*x(jjx,jjz,jjy+1,2) &
               +voll(4,isgn,jx,jz,jy)*x(jjx+1,jjz,jjy+1,2) &
	       +voll(5,isgn,jx,jz,jy)*x(jjx,jjz+1,jjy,2) &
               +voll(6,isgn,jx,jz,jy)*x(jjx+1,jjz+1,jjy,2) &
	       +voll(7,isgn,jx,jz,jy)*x(jjx,jjz+1,jjy+1,2) &
               +voll(8,isgn,jx,jz,jy)*x(jjx+1,jjz+1,jjy+1,2)
      enddo
      csmp=csmp0all*(1-tanh((rr(jx,jz)-0.9)*40))/2.
!      csmp=csmp0all*(1-tanh((psi(jx,jz)-pstrans)/pstransw))/2.
      wsm(jx,jz,jy)=csmp*(0.75*difc(ppp(1),x(jx,jz,jy,2),ppp(2),ssl(1,jx,jz,jy),ssl0,ssl(2,jx,jz,jy)) &
                         +0.25*difc(ppp(3),x(jx,jz,jy,2),ppp(4),ssl(3,jx,jz,jy),ssl0,ssl(4,jx,jz,jy)))
      endif
21    continue 
      x(:,:,:,2)=x(:,:,:,2)+wsm(:,:,:)

      call mpi_transfersm(x(:,:,:,2),1)

3     continue 
      
      return
      end

!ws*****************************************************************

     subroutine tracelineat_5p(jxl,jzl,jyl,lx,lz,ly,sl,vol)
      USE DECLARE
      include 'mpif.h'
      integer,dimension(4) :: lx,lz,ly
      real*8,dimension(4) :: sl
      real*8,dimension(8,4) :: vol
      real*8,dimension(3) :: bbb
      real*8 sl0,xxp,zzp,yyp,dyyp,dy1p,dx1p,dz1p,rxp,rxm,rzp,rzm,ryp,rym
      real*8 xxp0,zzp0,yyp0,xxp12,zzp12,yyp12
      real*8 volume,volppp,volmpp,volppm,volmpm,volpmp,volmmp,volpmm,volmmm
      integer isgn,jjx,jjz,jjy,jxl,jzl,jyl,j_add
       
      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)
!      volume=dyy*dxx*dzz
!
!      dyyp=yy(jy+1)-yy(jy)
!      dyym=yy(jy-1)-yy(jy)

!      dy1p=dyyp/nn
!      dy1m=dyym/nn
      sl0=0
      sl(:)=0.

      do 1 isgn=1,2


      dyyp=yy(jyl+(-1)**isgn)-yy(jyl)
      dy1p=dyyp/nline    

100   continue   
      xxp0=xx(jxl)
      zzp0=zz(jzl)
      yyp0=yy(jyl)

      bbb(:)=x(jxl,jzl,jyl,6:8)
      do 2 kp=1,2
      do nc=1,mcycline
!      dy1p=dyyp/nline/xxp0*xzero
      dx1p=xxp0*dy1p*bbb(1)/bbb(2)
      dz1p=xxp0*dy1p*bbb(3)/bbb(2)
    			
	!计算前进一步后的新位置
      xxp=xxp0+dx1p
	  zzp=zzp0+dz1p
      yyp=yyp0+dy1p 
      if(yyp.gt. 2*pi) yyp=yyp-2*pi
      if(yyp.lt. 0) yyp=yyp+2*pi

!      if( (xxp.gt.xx(jxl+1)) .or. (xxp.lt.xx(jxl-1)) &
!      .or.(zzp.gt.zz(jzl+1)) .or. (zzp.lt.zz(jzl-1))) then
!      if(nc.eq.1) then
!      dy1p=dy1p*0.5
!      goto 100
!      else
!      goto 200
!      endif
!      endif

!      sl(isgn)=sl(isgn)+(-1)**isgn*sqrt(dx1p**2+dz1p**2+xxp**2*dy1p**2)
      
      jjx=jxl+j_add(xxp,xx,mx,jxl)
	  jjz=jzl+j_add(zzp,zz,mz,jzl)
	  jjy=jyl+j_add(yyp,yy,my,jyl)
      if(jjy.lt.1) jjy=jjy+my
      if(jjy.gt.my)jjy=jjy-my
!      write(*,*) nrank,'nc',nc,jjy,xxp,jjx,zzp,jjz

      rxp=xx(jjx+1)-xxp
      rxm=xxp-xx(jjx)
      rzp=zz(jjz+1)-zzp
      rzm=zzp-zz(jjz)
      ryp=yy(jjy+1)-yyp
      rym=yyp-yy(jjy)

      volume=(xx(jjx+1)-xx(jjx))*(zz(jjz+1)-zz(jjz))*(yy(jjy+1)-yy(jjy))

      volppp=rxp*rzp*ryp/volume
      volppm=rxp*rzp*rym/volume
      volpmp=rxp*rzm*ryp/volume
      volpmm=rxp*rzm*rym/volume

      volmpp=rxm*rzp*ryp/volume
      volmpm=rxm*rzp*rym/volume
      volmmp=rxm*rzm*ryp/volume
      volmmm=rxm*rzm*rym/volume

      do m=6,8
      bbb(m-5)=(volppp*x(jjx,jjz,jjy,m)       +volmpp*x(jjx+1,jjz,jjy,m) &
	           +volppm*x(jjx,jjz,jjy+1,m)  +volmpm*x(jjx+1,jjz,jjy+1,m) &
	           +volpmp*x(jjx,jjz+1,jjy,m)     +volmmp*x(jjx+1,jjz+1,jjy,m)&
	           +volpmm*x(jjx,jjz+1,jjy+1,m)+volmmm*x(jjx+1,jjz+1,jjy+1,m))/volume
      enddo


      xxp12=xxp0+0.5*dx1p
	  zzp12=zzp0+0.5*dz1p
      yyp12=yyp0+0.5*dy1p

!      dy1p=dyyp/nline/xxp*xzero
      dx1p=xxp*dy1p*bbb(1)/bbb(2)
      dz1p=xxp*dy1p*bbb(3)/bbb(2)

      xxp=xxp12+0.5*dx1p
      zzp=zzp12+0.5*dz1p
      yyp=yyp12+0.5*dy1p
      if(yyp.gt. 2*pi) yyp=yyp-2*pi
      if(yyp.lt. 0) yyp=yyp+2*pi

      sl(isgn+(kp-1)*2)=sl(isgn+(kp-1)*2)+(-1)**isgn*sqrt(dx1p**2+dz1p**2+xxp**2*dy1p**2)
      
      jjx=jxl+j_add(xxp,xx,mx,jxl)
	  jjz=jzl+j_add(zzp,zz,mz,jzl)
	  jjy=jyl+j_add(yyp,yy,my,jyl)
      if(jjy.lt.1) jjy=jjy+my
      if(jjy.gt.my)jjy=jjy-my
!      write(*,*) nrank,'nc',nc,jjy,xxp,jjx,zzp,jjz

      rxp=xx(jjx+1)-xxp
      rxm=xxp-xx(jjx)
      rzp=zz(jjz+1)-zzp
      rzm=zzp-zz(jjz)
      ryp=yy(jjy+1)-yyp
      rym=yyp-yy(jjy)

      volume=(xx(jjx+1)-xx(jjx))*(zz(jjz+1)-zz(jjz))*(yy(jjy+1)-yy(jjy))

      volppp=rxp*rzp*ryp/volume
      volppm=rxp*rzp*rym/volume
      volpmp=rxp*rzm*ryp/volume
      volpmm=rxp*rzm*rym/volume

      volmpp=rxm*rzp*ryp/volume
      volmpm=rxm*rzp*rym/volume
      volmmp=rxm*rzm*ryp/volume
      volmmm=rxm*rzm*rym/volume

      do m=6,8
      bbb(m-5)=(volppp*x(jjx,jjz,jjy,m)       +volmpp*x(jjx+1,jjz,jjy,m) &
	           +volppm*x(jjx,jjz,jjy+1,m)  +volmpm*x(jjx+1,jjz,jjy+1,m) &
	           +volpmp*x(jjx,jjz+1,jjy,m)     +volmmp*x(jjx+1,jjz+1,jjy,m)&
	           +volpmm*x(jjx,jjz+1,jjy+1,m)+volmmm*x(jjx+1,jjz+1,jjy+1,m))/volume
      enddo

      xxp0=xxp
      zzp0=zzp
      yyp0=yyp

      enddo

200   continue

      vol(1,isgn+(kp-1)*2)=volppp
      vol(2,isgn+(kp-1)*2)=volmpp
      vol(3,isgn+(kp-1)*2)=volppm
      vol(4,isgn+(kp-1)*2)=volmpm
      vol(5,isgn+(kp-1)*2)=volpmp
      vol(6,isgn+(kp-1)*2)=volmmp
      vol(7,isgn+(kp-1)*2)=volpmm
      vol(8,isgn+(kp-1)*2)=volmmm
      
      lx(isgn+(kp-1)*2)=jjx
      lz(isgn+(kp-1)*2)=jjz
      ly(isgn+(kp-1)*2)=jjy 
      
2     continue  
                
1     continue

      return
      end

      integer function sgn(val)
      real*8 val
      if(val .ge. 0) sgn=1
      if(val .lt. 0) sgn=-1
      return
      end


      integer function j_add(val,ary,nn,jj)
      integer nn,jj,jl
      real*8 val,ary(nn)
     
      jl=0
      if(val.gt.ary(jj)) then
      do while(val .gt. ary(jj+jl))
      jl=jl+1
!      if(jl .gt. 2) then      
!      write(*,*) nrank,jj,'error:j_add gt 2'
!      return
!      endif
      enddo
      j_add=jl-1
      endif

      if(val.lt.ary(jj)) then
      do while(val .lt. ary(jj+jl))
      jl=jl-1
!      if(jl .lt. -2) then      
!      write(*,*) nrank,jj,'error:j_add lt -2'
!      return
!      endif
      enddo
      j_add=jl
      endif

      return
      end


!ws********************************************
      subroutine smthp2_traceline_tm(kk)
      USE DECLARE
      include 'mpif.h'
      integer,dimension(2) :: lx,lz,ly
      real*8,dimension(2) :: sl,ppp,ppb
      real*8,dimension(8,2) :: vol
      integer,dimension(2,mx,mz,my) :: llx,llz,lly
      real*8,dimension(2,mx,mz,my) :: ssl
      real*8,dimension(8,2,mx,mz,my) :: voll
      real*8,dimension(mx,mz,my) :: wsm,pb,tm1 !,wsm1,wsm2
      integer kk,k,isgn,jjx,jjz,jjy
      real*8 ssl0,bb

      !  dif1= (d f / dx)*dx2 = [((x0-xm1)/(xp1-x0)*(fp1-f0)-(xp1-x0)/(x0-xm1)*(f0-fm1))/(xp1-xm1)] *(xp1-x0)*(x0-xm1) 
      dif1(fm1,f0,fp1,xm1,x0,xp1)= &
       ((x0-xm1)**2*(fp1-f0)-(xp1-x0)**2*(f0-fm1))/(xp1-xm1)


      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)

      ssl0=0.

      do 2 jy=1,my
      do 2 jz=iz_first,iz_last
      do 2 jx=ix_first,ix_last
      if(psi(jx,jz).lt.ps(mpsa-2)) then
      call tracelineat2(jx,jz,jy,lx,lz,ly,sl,vol)

      llx(:,jx,jz,jy)=lx(:)
      llz(:,jx,jz,jy)=lz(:)
      lly(:,jx,jz,jy)=ly(:)
      ssl(:,jx,jz,jy)=sl(:)
      voll(:,:,jx,jz,jy)=vol(:,:)
      endif      
      bb=sqrt(x(jx,jz,jy,6)**2+x(jx,jz,jy,7)**2+x(jx,jz,jy,8)**2)
      pb(jx,jz,jy)=(x1(jx,jz,jy,6)*xint_dx(jx,jz,2)+x1(jx,jz,jy,8)*xint_dz(jx,jz,2))/bb/x(jx,jz,jy,1)
      tm1(jx,jz,jy)=x1(jx,jz,jy,2)/x(jx,jz,jy,1)
2     continue  

      do 3 k=1,kk
      
      do 21 jy=1,my
      do 21 jz=iz_first+2,iz_last-2
      do 21 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
      do isgn=1,2
      jjx=llx(isgn,jx,jz,jy)
      jjz=llz(isgn,jx,jz,jy)
      jjy=lly(isgn,jx,jz,jy)

      ppp(isgn)=voll(1,isgn,jx,jz,jy)*tm1(jjx,jjz,jjy)  &
               +voll(2,isgn,jx,jz,jy)*tm1(jjx+1,jjz,jjy) &
	           +voll(3,isgn,jx,jz,jy)*tm1(jjx,jjz,jjy+1) &
               +voll(4,isgn,jx,jz,jy)*tm1(jjx+1,jjz,jjy+1) &
	           +voll(5,isgn,jx,jz,jy)*tm1(jjx,jjz+1,jjy) &
               +voll(6,isgn,jx,jz,jy)*tm1(jjx+1,jjz+1,jjy) &
	           +voll(7,isgn,jx,jz,jy)*tm1(jjx,jjz+1,jjy+1) &
               +voll(8,isgn,jx,jz,jy)*tm1(jjx+1,jjz+1,jjy+1)

      ppb(isgn)=voll(1,isgn,jx,jz,jy)*pb(jjx,jjz,jjy)  &
               +voll(2,isgn,jx,jz,jy)*pb(jjx+1,jjz,jjy) &
	           +voll(3,isgn,jx,jz,jy)*pb(jjx,jjz,jjy+1) &
               +voll(4,isgn,jx,jz,jy)*pb(jjx+1,jjz,jjy+1) &
	           +voll(5,isgn,jx,jz,jy)*pb(jjx,jjz+1,jjy) &
               +voll(6,isgn,jx,jz,jy)*pb(jjx+1,jjz+1,jjy) &
	           +voll(7,isgn,jx,jz,jy)*pb(jjx,jjz+1,jjy+1) &
               +voll(8,isgn,jx,jz,jy)*pb(jjx+1,jjz+1,jjy+1)

      enddo
!       wsm1(jx,jz,jy)=difc(ppp(1),x1(jx,jz,jy,2),ppp(2),ssl(1,jx,jz,jy),ssl0,ssl(2,jx,jz,jy))
!       wsm2(jx,jz,jy)=dif1(ppb(1),pb(jx,jz,jy),ppb(2),ssl(1,jx,jz,jy),ssl0,ssl(2,jx,jz,jy))
!	wsm(jx,jz,jy)=csmp0all*(wsm1(jx,jz,jy)+wsm2(jx,jz,jy))
      wsm(jx,jz,jy)=csmp0all*(difc(ppp(1),x1(jx,jz,jy,2),ppp(2),ssl(1,jx,jz,jy),ssl0,ssl(2,jx,jz,jy)) &
                             +dif1(ppb(1),pb(jx,jz,jy),ppb(2),ssl(1,jx,jz,jy),ssl0,ssl(2,jx,jz,jy)))*x(jx,jz,jy,1)
      endif
21    continue 
      x(:,:,:,2)=x(:,:,:,2)+wsm(:,:,:)

      call mpi_transfersm(x(:,:,:,2),1)

3     continue 
      
      return
      end

!ws*****************************************************************

     subroutine tracelineat_spec1(jxl,jzl,jyl,lx,lz,ly,sl,yypl,area)
      USE DECLARE
      include 'mpif.h'
      integer,dimension(2) :: lx,lz,ly
      real*8,dimension(2) :: sl,yypl
      real*8,dimension(4,2) :: area
!      real*8,dimension(8,2) :: vol
      real*8,dimension(3) :: bbb
      real*8 bkc,bks
      real*8 sl0,xxp,zzp,yyp,dyyp,dy1p,dx1p,dz1p,rxp,rxm,rzp,rzm,ryp,rym,app,apm,amp,amm,atotal
!      real*8 volume,volppp,volmpp,volppm,volmpm,volpmp,volmmp,volpmm,volmmm
      integer isgn,jjx,jjz,jjy,jxl,jzl,jyl
       
      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)
!      volume=dyy*dxx*dzz
!
!      dyyp=yy(jy+1)-yy(jy)
!      dyym=yy(jy-1)-yy(jy)

!      dy1p=dyyp/nn
!      dy1m=dyym/nn
      sl0=0
      sl(:)=0.

      do 1 isgn=1,2


      dyyp=yy(jyl+(-1)**isgn)-yy(jyl)
      dy1p=dyyp/nline    

100   continue   
      xxp=xx(jxl)
      zzp=zz(jzl)
      yyp=yy(jyl)
      
      bbb(:)=x(jxl,jzl,jyl,6:8)

      do nc=1,mcycline
      dx1p=xxp*dy1p*bbb(1)/bbb(2)
      dz1p=xxp*dy1p*bbb(3)/bbb(2)
    			
	!计算前进一步后的新位置
      xxp=xxp+dx1p
	  zzp=zzp+dz1p
      yyp=yyp+dy1p 
      if(yyp.gt. 2*pi) yyp=yyp-2*pi
      if(yyp.lt. 0) yyp=yyp+2*pi

!      if( (xxp.gt.xx(jxl+1)) .or. (xxp.lt.xx(jxl-1)) &
!      .or.(zzp.gt.zz(jzl+1)) .or. (zzp.lt.zz(jzl-1))) then
!      if(nc.eq.1) then
!      dy1p=dy1p*0.5
!      goto 100
!      else
!      goto 200
!      endif
!      endif

      sl(isgn)=sl(isgn)+(-1)**isgn*sqrt(dx1p**2+dz1p**2+xxp**2*dy1p**2)
      
      jjx=floor((xxp-xx(ix_first))/dxx)+ix_first
	  jjz=floor((zzp-zz(iz_first))/dzz)+iz_first   
	  jjy=floor((yyp-yy(1))/dyy)+1
!      if(jjy.lt.1) jjy=jjy+my
!      if(jjy.gt.my)jjy=jjy-my
!      write(*,*) nrank,'nc',nc,jjy,xxp,jjx,zzp,jjz

      rxp=xx(jjx+1)-xxp
      rxm=xxp-xx(jjx)
      rzp=zz(jjz+1)-zzp
      rzm=zzp-zz(jjz)
      atotal=(xx(jjx+1)-xx(jjx))*(zz(jjz+1)-zz(jjz))

      app=rxp*rzp/atotal
      apm=rxp*rzm/atotal
      amp=rxm*rzp/atotal
      amm=rxm*rzm/atotal

      
      do m=1,3
      bbb(m)=app*xint(jjx,jjz,m+5)+apm*xint(jjx,jjz+1,m+5)+amp*xint(jjx+1,jjz,m+5)+amm*xint(jjx+1,jjz+1,m+5)
      do lk=1,myt/2
      bkc=app*xkc(jjx,jjz,lk,m+5)+apm*xkc(jjx,jjz+1,lk,m+5)+amp*xkc(jjx+1,jjz,lk,m+5)+amm*xkc(jjx+1,jjz+1,lk,m+5)
      bks=app*xks(jjx,jjz,lk,m+5)+apm*xks(jjx,jjz+1,lk,m+5)+amp*xks(jjx+1,jjz,lk,m+5)+amm*xks(jjx+1,jjz+1,lk,m+5)
      bbb(m)=bbb(m)+bkc*cos((lk-1)*yyp)+bks*sin((lk-1)*yyp)
      enddo
      enddo
!      ryp=yy(jjy+1)-yyp
!      rym=yyp-yy(jjy)
!
!      volume=(xx(jjx+1)-xx(jjx))*(zz(jjz+1)-zz(jjz))*(yy(jjy+1)-yy(jjy))
!
!      volppp=rxp*rzp*ryp/volume
!      volppm=rxp*rzp*rym/volume
!      volpmp=rxp*rzm*ryp/volume
!      volpmm=rxp*rzm*rym/volume
!
!      volmpp=rxm*rzp*ryp/volume
!      volmpm=rxm*rzp*rym/volume
!      volmmp=rxm*rzm*ryp/volume
!      volmmm=rxm*rzm*rym/volume
!
!      do m=6,8
!      bbb(m-5)=(volppp*x(jjx,jjz,jjy,m)       +volmpp*x(jjx+1,jjz,jjy,m) &
!	           +volppm*x(jjx,jjz,jyp(jjy),m)  +volmpm*x(jjx+1,jjz,jyp(jjy),m) &
!	           +volpmp*x(jjx,jjz+1,jjy,m)     +volmmp*x(jjx+1,jjz+1,jjy,m)&
!	           +volpmm*x(jjx,jjz+1,jyp(jjy),m)+volmmm*x(jjx+1,jjz+1,jyp(jjy),m))/volume
!      enddo

      enddo

200   continue

!      vol(1,isgn)=volppp
!      vol(2,isgn)=volmpp
!      vol(3,isgn)=volppm
!      vol(4,isgn)=volmpm
!      vol(5,isgn)=volpmp
!      vol(6,isgn)=volmmp
!      vol(7,isgn)=volpmm
!      vol(8,isgn)=volmmm
      area(1,isgn)=app
      area(2,isgn)=apm
      area(3,isgn)=amp
      area(4,isgn)=amm
      yypl(isgn)=yyp
      lx(isgn)=jjx
      lz(isgn)=jjz
      ly(isgn)=jjy     
                
1     continue

      return
      end

!ws********************************************
      subroutine smthp_traceline_spec(kk)
      USE DECLARE
      include 'mpif.h'
      integer,dimension(2) :: lx,lz,ly
      real*8,dimension(2) :: sl,yypl,ppp
      real*8,dimension(4,2) :: area
      integer,dimension(2,mx,mz,my) :: llx,llz,lly
      real*8,dimension(2,mx,mz,my) :: ssl,yyl
!      real*8,dimension(8,2,mx,mz,my) :: voll
      real*8,dimension(4,2,mx,mz,my) :: arl
      real*8,dimension(mx,mz,my) :: wsm
      integer kk,k,isgn,jjx,jjz,jjy
      real*8 ssl0,csmp,pkc,pks

      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)

      ssl0=0.

      do 2 jy=1,my
      do 2 jz=iz_first+2,iz_last-2
      do 2 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
!      if(psi(jx,jz).lt.psia1) then
      call tracelineat_spec1(jx,jz,jy,lx,lz,ly,sl,yypl,area)

      llx(:,jx,jz,jy)=lx(:)
      llz(:,jx,jz,jy)=lz(:)
      lly(:,jx,jz,jy)=ly(:)
      ssl(:,jx,jz,jy)=sl(:)
      yyl(:,jx,jz,jy)=yypl(:)
      arl(:,:,jx,jz,jy)=area(:,:)
      endif
2     continue  

      do 3 k=1,kk
      
      do 21 jy=1,my
      do 21 jz=iz_first+2,iz_last-2
      do 21 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
!      if(psi(jx,jz).lt.psia1) then
      do isgn=1,2
      jjx=llx(isgn,jx,jz,jy)
      jjz=llz(isgn,jx,jz,jy)
!      jjy=lly(isgn,jx,jz,jy)
      ppp(isgn)=arl(1,isgn,jx,jz,jy)*xint(jjx,jjz,2) &
         +arl(2,isgn,jx,jz,jy)*xint(jjx,jjz+1,2) &
         +arl(3,isgn,jx,jz,jy)*xint(jjx+1,jjz,2) &
         +arl(4,isgn,jx,jz,jy)*xint(jjx+1,jjz+1,2)
      do ky=1,myt/2
      pkc=arl(1,isgn,jx,jz,jy)*xkc(jjx,jjz,ky,2) &
         +arl(2,isgn,jx,jz,jy)*xkc(jjx,jjz+1,ky,2) &
         +arl(3,isgn,jx,jz,jy)*xkc(jjx+1,jjz,ky,2) &
         +arl(4,isgn,jx,jz,jy)*xkc(jjx+1,jjz+1,ky,2)

      pks=arl(1,isgn,jx,jz,jy)*xks(jjx,jjz,ky,2) &
         +arl(2,isgn,jx,jz,jy)*xks(jjx,jjz+1,ky,2) &
         +arl(3,isgn,jx,jz,jy)*xks(jjx+1,jjz,ky,2) &
         +arl(4,isgn,jx,jz,jy)*xks(jjx+1,jjz+1,ky,2)

      ppp(isgn)=ppp(isgn)+pkc*cos((ky-1)*yyl(isgn,jx,jz,jy))+pks*sin((ky-1)*yyl(isgn,jx,jz,jy))
      enddo
!      ppp(isgn)=voll(1,isgn,jx,jz,jy)*x(jjx,jjz,jjy,2)  &
!               +voll(2,isgn,jx,jz,jy)*x(jjx+1,jjz,jjy,2) &
!	       +voll(3,isgn,jx,jz,jy)*x(jjx,jjz,jyp(jjy),2) &
!               +voll(4,isgn,jx,jz,jy)*x(jjx+1,jjz,jyp(jjy),2) &
!	       +voll(5,isgn,jx,jz,jy)*x(jjx,jjz+1,jjy,2) &
!               +voll(6,isgn,jx,jz,jy)*x(jjx+1,jjz+1,jjy,2) &
!	       +voll(7,isgn,jx,jz,jy)*x(jjx,jjz+1,jyp(jjy),2) &
!               +voll(8,isgn,jx,jz,jy)*x(jjx+1,jjz+1,jyp(jjy),2)
      enddo
      csmp=csmp0all*(1-tanh((rr(jx,jz)-0.9)*40))/2.
!      csmp=csmp0all*(1-tanh((psi(jx,jz)-pstrans)/pstransw))/2.
      wsm(jx,jz,jy)=csmp*difc(ppp(1),x(jx,jz,jy,2),ppp(2),ssl(1,jx,jz,jy),ssl0,ssl(2,jx,jz,jy))
      endif
21    continue 
      x(:,:,:,2)=x(:,:,:,2)+wsm(:,:,:)

      call mpi_transfersm(x(:,:,:,2),1)

3     continue 
      
      return
      end
!ws********************************************
      subroutine smthp1_traceline_spec(kk)
      USE DECLARE
      include 'mpif.h'
      integer,dimension(2) :: lx,lz,ly
      real*8,dimension(2) :: sl,yypl,ppp
      real*8,dimension(4,2) :: area
      integer,dimension(2,mx,mz,my) :: llx,llz,lly
      real*8,dimension(2,mx,mz,my) :: ssl,yyl
!      real*8,dimension(8,2,mx,mz,my) :: voll
      real*8,dimension(4,2,mx,mz,my) :: arl
      real*8,dimension(mx,mz,my) :: wsm
      integer kk,k,isgn,jjx,jjz,jjy
      real*8 ssl0,csmp,pkc,pks

      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)

      ssl0=0.

      do 2 jy=1,my
      do 2 jz=iz_first+2,iz_last-2
      do 2 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
!      if(psi(jx,jz).lt.psia1) then
      call tracelineat_spec1(jx,jz,jy,lx,lz,ly,sl,yypl,area)

      llx(:,jx,jz,jy)=lx(:)
      llz(:,jx,jz,jy)=lz(:)
      lly(:,jx,jz,jy)=ly(:)
      ssl(:,jx,jz,jy)=sl(:)
      yyl(:,jx,jz,jy)=yypl(:)
      arl(:,:,jx,jz,jy)=area(:,:)
      endif
2     continue  

      do 3 k=1,kk
      
      do 21 jy=1,my
      do 21 jz=iz_first+2,iz_last-2
      do 21 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
!      if(psi(jx,jz).lt.psia1) then
      do isgn=1,2
      jjx=llx(isgn,jx,jz,jy)
      jjz=llz(isgn,jx,jz,jy)
!      jjy=lly(isgn,jx,jz,jy)
      ppp(isgn)=0
      do ky=1,myt/2
      pkc=arl(1,isgn,jx,jz,jy)*xkc(jjx,jjz,ky,2) &
         +arl(2,isgn,jx,jz,jy)*xkc(jjx,jjz+1,ky,2) &
         +arl(3,isgn,jx,jz,jy)*xkc(jjx+1,jjz,ky,2) &
         +arl(4,isgn,jx,jz,jy)*xkc(jjx+1,jjz+1,ky,2)

      pks=arl(1,isgn,jx,jz,jy)*xks(jjx,jjz,ky,2) &
         +arl(2,isgn,jx,jz,jy)*xks(jjx,jjz+1,ky,2) &
         +arl(3,isgn,jx,jz,jy)*xks(jjx+1,jjz,ky,2) &
         +arl(4,isgn,jx,jz,jy)*xks(jjx+1,jjz+1,ky,2)

      ppp(isgn)=ppp(isgn)+pkc*cos((ky-1)*yyl(isgn,jx,jz,jy))+pks*sin((ky-1)*yyl(isgn,jx,jz,jy))
      enddo
!      ppp(isgn)=voll(1,isgn,jx,jz,jy)*x(jjx,jjz,jjy,2)  &
!               +voll(2,isgn,jx,jz,jy)*x(jjx+1,jjz,jjy,2) &
!	       +voll(3,isgn,jx,jz,jy)*x(jjx,jjz,jyp(jjy),2) &
!               +voll(4,isgn,jx,jz,jy)*x(jjx+1,jjz,jyp(jjy),2) &
!	       +voll(5,isgn,jx,jz,jy)*x(jjx,jjz+1,jjy,2) &
!               +voll(6,isgn,jx,jz,jy)*x(jjx+1,jjz+1,jjy,2) &
!	       +voll(7,isgn,jx,jz,jy)*x(jjx,jjz+1,jyp(jjy),2) &
!               +voll(8,isgn,jx,jz,jy)*x(jjx+1,jjz+1,jyp(jjy),2)
      enddo
!      csmp=csmp0all*(1-tanh((rr(jx,jz)-0.9)*40))/2.
!      csmp=csmp0all*(1-tanh((psi(jx,jz)-pstrans)/pstransw))/2.
      wsm(jx,jz,jy)=csmp0*difc(ppp(1),x1(jx,jz,jy,2),ppp(2),ssl(1,jx,jz,jy),ssl0,ssl(2,jx,jz,jy))
      endif
21    continue 
      x1(:,:,:,2)=x1(:,:,:,2)+wsm(:,:,:)

      call mpi_transfersm(x1(:,:,:,2),1)

3     continue 
      
      return
      end




!ws********************************************
      subroutine calculate_eta
      USE DECLARE
      include 'mpif.h'
!  d1fc= d f / dx  with fourth-order accuracy central difference
      d1fc(fm2,fm1,f0,fp1,fp2,a,b,c,d)= &
       a*(fp1-f0)+b*(f0-fm1)+c*(fp2-f0)+d*(f0-fm2)

      do jy=iy_first,iy_last
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
      if(psi(jx,jz) .lt. psia) then
      tm(jx,jz,jy)=x(jx,jz,jy,2)/x(jx,jz,jy,1)
      eta(jx,jz,jy)=eta0*(1.0+(etacut-1.0)*tanh(((tm(jx,jz,jy)/tm00)**(-1.5)-1.0)/(etacut-1.0)))
      endif
      enddo
      enddo
      enddo
      
      do jy=iy_first+2,iy_last-2  
      do jz=iz_first+2,iz_last-2
      do jx=ix_first+2,ix_last-2 
      if(psi(jx,jz) .lt. psia) then         
      etax(jx,jz,jy)=d1fc(eta(jx-2,jz,jy),eta(jx-1,jz,jy),eta(jx,jz,jy) &
         ,eta(jx+1,jz,jy),eta(jx+2,jz,jy),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
      etaz(jx,jz,jy)=d1fc(eta(jx,jz-2,jy),eta(jx,jz-1,jy),eta(jx,jz,jy) &
         ,eta(jx,jz+1,jy),eta(jx,jz+2,jy),az1(jz),bz1(jz),cz1(jz),dz1(jz))
      etay(jx,jz,jy)=d1fc(eta(jx,jz,jy-2),eta(jx,jz,jy-1),eta(jx,jz,jy) &
         ,eta(jx,jz,jy+1),eta(jx,jz,jy+2),ay1(jy),by1(jy),cy1(jy),dy1(jy))
      endif
      enddo
      enddo
      enddo

      return
      end


      subroutine current_driven
      USE DECLARE
      include 'mpif.h'
      real*8 tcd,bb2,bb,cud0
!      call calculate_ps1
!      call distribution_cd
!      if(cd_Opoint) call distribution_cd_Opoint
      tcd=0.5*(tanh(pi*(time-tcds)/tcdd)+tanh(pi*(tcde-time)/tcdd))

      do 1 jy=iy_first,iy_last
      do 1 jz=iz_first,iz_last
      do 1 jx=ix_first,ix_last
      if(psi(jx,jz) .lt. psia) then 
!      cud(jx,jz,jy,2)=cd00*fn_cdy(jx,jz)*tcd
!      cud(jx,jz,jy,1)=cud(jx,jz,jy,2)*x(jx,jz,jy,6)/x(jx,jz,jy,7)
!      cud(jx,jz,jy,3)=cud(jx,jz,jy,2)*x(jx,jz,jy,8)/x(jx,jz,jy,7)
        bb2=x(jx,jz,jy,6)**2+x(jx,jz,jy,8)**2+x(jx,jz,jy,7)**2
        bb=sqrt(bb2)
      cud0=cd00*fn_cdy(jx,jz,jy)*tcd
      cud(jx,jz,jy,2)=cud0*x(jx,jz,jy,7)/bb
      cud(jx,jz,jy,1)=cud0*x(jx,jz,jy,6)/bb
      cud(jx,jz,jy,3)=cud0*x(jx,jz,jy,8)/bb
      endif
   1  continue
      
      return
      end

      subroutine distribution_cd
      USE DECLARE
      USE DECLARE_OXpoint

      real*8 fnV,fnV1
      include 'mpif.h'      
      
      if(lrstrt_cd) then
      open(unit=777,file='cd.dat',status='unknown',form='formatted')
      read(777,3001) 
      read(777,*) cIp,fcd,tcds,tcde,tcdd,psmode,delcd,fnV,psshift,br_max0
      close(777)
3001  format('cIp,fcd,tcds,tcde,tcdd,psmode,delcd,fnV')
      do jy=iy_first,iy_last
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
!      phase=nmode*yy(jy)+mmode*thxz(jx,jz)+0
!      ps1(jx,jz,jy)
      fn_cdy(jx,jz,jy)=exp(-(psi(jx,jz)+ps1(jx,jz,jy)-(psmode+psshift))**2/delcd**2)/fnV
      enddo
      enddo
      enddo

      write(*,*) 'brmax',br_max0,br_max,nrank
      return
      endif


      do jy=iy_first,iy_last
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
!      phase=nmode*yy(jy)+mmode*thxz(jx,jz)+0
      fn_cdy(jx,jz,jy)=exp(-(psi(jx,jz)+ps1(jx,jz,jy)-(psmode+psshift))**2/delcd**2)
      enddo
      enddo
      enddo

      fnV1=0
!      do jy=iy_first+2,iy_last-2  
      do jz=iz_first+2,iz_last-2
      do jx=ix_first+2,ix_last-2 
      if(psi(jx,jz) .lt. psia) then 
!      fnV1=fnV1+fn_cdy(jx,jz)*(xx(jx+1)-xx(jx-1))*(zz(jz+1)-zz(jz-1))/4
        bb2=x(jx,jz,3,6)**2+x(jx,jz,3,8)**2+x(jx,jz,3,7)**2
        bb=sqrt(bb2)
      fnV1=fnV1+fn_cdy(jx,jz,3)*x(jx,jz,3,7)/bb*(xx(jx+1)-xx(jx-1))*(zz(jz+1)-zz(jz-1))/4
      endif
      enddo
      enddo
!      enddo

      CALL MPI_ALLREDUCE(fnV1,fnV,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                        MPI_COMM_WORLD,IERROR)

      fn_cdy(:,:,:)=fn_cdy(:,:,:)/fnV
      if(nrank==0 .and. time .le. timein+1.e-6) then
      
      write(*,*) '***current drive:cd=fcd*Ip*fn,Ing(fn*dS)=1***'
      write(*,*) 'Ip    =',cIp
      write(*,*) 'fcd   =',fcd
      write(*,*) 'tcds  =',tcds
      write(*,*) 'tcde  =',tcde
      write(*,*) 'tcdd  =',tcdd
      write(*,*) 'psmode=',psmode
      write(*,*) 'delcd =',delcd
      write(*,*) 'fnV   =',fnV
      write(*,*) 'pssft =',psshift
      write(*,*) '***current drive***'
      open(unit=777,file='cd.dat',status='unknown',form='formatted')
      write(777,3000) 
      write(777,*) cIp,fcd,tcds,tcde,tcdd,psmode,delcd,fnV,psshift,br_max0
3000  format('cIp,fcd,tcds,tcde,tcdd,psmode,delcd,fnV')
      close(777)
      endif

      return
      end
!ws*****************************************************
      subroutine distribution_cd_cos
      USE DECLARE
      real*8 fnV,fnV1
      include 'mpif.h'      
      
      if(lrstrt_cd) then
      open(unit=777,file='cd.dat',status='unknown',form='formatted')
      read(777,3001) 
      read(777,*) cIp,fcd,tcds,tcde,tcdd,psmode,delcd,fnV,psshift
      close(777)
3001  format('cIp,fcd,tcds,tcde,tcdd,psmode,delcd,fnV')
      do jy=iy_first,iy_last
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
!      phase=nmode*yy(jy)+mmode*thxz(jx,jz)+0
!      ps1(jx,jz,jy)
      fn_cdy(jx,jz,jy)=exp(-(psi(jx,jz)+ps1(jx,jz,jy)-(psmode+psshift))**2/delcd**2)*(1+dcos(nmode*yy(jy)+mmode*thxz(jx,jz)))/fnV
      enddo
      enddo
      enddo
      return
      endif


      do jy=iy_first,iy_last
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
!      phase=nmode*yy(jy)+mmode*thxz(jx,jz)+0
      fn_cdy(jx,jz,jy)=exp(-(psi(jx,jz)+ps1(jx,jz,jy)-(psmode+psshift))**2/delcd**2)*(1+dcos(nmode*yy(jy)+mmode*thxz(jx,jz)))
      enddo
      enddo
      enddo

      fnV1=0
!      do jy=iy_first+2,iy_last-2  
      do jz=iz_first+2,iz_last-2
      do jx=ix_first+2,ix_last-2 
      if(psi(jx,jz) .lt. psia) then 
!      fnV1=fnV1+fn_cdy(jx,jz)*(xx(jx+1)-xx(jx-1))*(zz(jz+1)-zz(jz-1))/4
        bb2=x(jx,jz,3,6)**2+x(jx,jz,3,8)**2+x(jx,jz,3,7)**2
        bb=sqrt(bb2)
      fnV1=fnV1+fn_cdy(jx,jz,3)*x(jx,jz,3,7)/bb*(xx(jx+1)-xx(jx-1))*(zz(jz+1)-zz(jz-1))/4
      endif
      enddo
      enddo
!      enddo

      CALL MPI_ALLREDUCE(fnV1,fnV,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                        MPI_COMM_WORLD,IERROR)

      fn_cdy(:,:,:)=fn_cdy(:,:,:)/fnV
      if(nrank==0 .and. time .le. timein+1.e-6) then
      write(*,*) '***current drive:cd=fcd*Ip*fn,Ing(fn*dS)=1***'
      write(*,*) 'Ip    =',cIp
      write(*,*) 'fcd   =',fcd
      write(*,*) 'tcds  =',tcds
      write(*,*) 'tcde  =',tcde
      write(*,*) 'tcdd  =',tcdd
      write(*,*) 'psmode=',psmode
      write(*,*) 'delcd =',delcd
      write(*,*) 'fnV   =',fnV
      write(*,*) 'pssft =',psshift
      write(*,*) '***current drive***'
      open(unit=777,file='cd.dat',status='unknown',form='formatted')
      write(777,3000) 
      write(777,*) cIp,fcd,tcds,tcde,tcdd,psmode,delcd,fnV,psshift
3000  format('cIp,fcd,tcds,tcde,tcdd,psmode,delcd,fnV')
      close(777)
      endif

      return
      end
!ws*****************************************************
      subroutine calculate_cb00
      USE DECLARE 
      include 'mpif.h'
      real*8 fnV,fnV1
      real*8,dimension(mx,mz) :: fn_bs 
!      do jy=iy_first,iy_last
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
       if(bp0(jx,jz)==0) then
       p0dr=0
       else
       p0dr=-xint(jx,jz,8)/bp0(jx,jz)*xint_dx(jx,jz,2)+xint(jx,jz,6)/bp0(jx,jz)*xint_dz(jx,jz,2)
       endif       
       fn_bs(jx,jz)=-(rr(jx,jz)/xx(jx))**0.5*p0dr/bp0(jx,jz)
      enddo
      enddo
!      enddo

      fnV1=0
!      do jy=iy_first+2,iy_last-2  
      do jz=iz_first+2,iz_last-2
      do jx=ix_first+2,ix_last-2 
      if(psi(jx,jz) .lt. psia) then 
      fnV1=fnV1+fn_bs(jx,jz)*(xx(jx+1)-xx(jx-1))*(zz(jz+1)-zz(jz-1))/4
      endif
      enddo
      enddo
!      enddo

      CALL MPI_ALLREDUCE(fnV1,fnV,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                        MPI_COMM_WORLD,IERROR)

      cb00=fbs*cIp/fnV
      return
      end
!ws********************************************
      subroutine calculate_ps1
      USE DECLARE
      include 'mpif.h'

      real*8,dimension(mz,my) :: wsx1s,wsx1r
      integer status(mpi_status_size)
      integer ii

      ps1=0
!      do jy=iy_first,iy_last
      do jz=iz_first,iz_last
      do jx=ix_first+2,ix_last-2
      if(psi(jx-1,jz).ge.psia .and. psi(jx,jz).lt.psia ) then
      ps1(jx,jz,:)=-(xx(jx)*x1(jx,jz,:,8))/2*(xx(jx)-xxam(nrkz(nrank)*mzm+jz-2))
      endif
      if(psi(jx,jz).lt.psia) then
      ps1(jx+1,jz,:)=ps1(jx,jz,:)-(xx(jx)*x1(jx,jz,:,8)+xx(jx+1)*x1(jx+1,jz,:,8))/2*(xx(jx+1)-xx(jx))
      endif
      enddo
      enddo
!     enddo

      if(nrkx(nrank) .eq. 0) then
      wsx1s=ps1(ix_last-1,:,:)
      CALL MPI_Send( wsx1s,mz*my, MPI_DOUBLE_PRECISION,nrank+1, 0,  &
		      MPI_COMM_WORLD,ierror )
      endif

      do ii=1,nprx-1
      if(nrkx(nrank).eq. ii) then
      CALL MPI_Recv( wsx1r,mz*my, MPI_DOUBLE_PRECISION,nrank-1, ii-1,  &
		      MPI_COMM_WORLD, status,ierror )
      do jx=ix_first+2,ix_last-2
      ps1(jx,:,:)=ps1(jx,:,:)+wsx1r(:,:)
      enddo

      wsx1s(:,:)=ps1(ix_last-1,:,:)
      if(nrkx(nrank).lt. nprx-1) then
      CALL MPI_Send( wsx1s,mz*my, MPI_DOUBLE_PRECISION,nrank+1, ii,  &
		      MPI_COMM_WORLD,ierror )
      endif
      endif

      enddo


      return
      end



!ws********************************************
      subroutine init_ps1
      USE DECLARE
      include 'mpif.h'
      integer jjxt,jjztp,jjztm

      jjxt=1
      do while(xxt(jjxt) .lt. xmg)
      jjxt=jjxt+1
      enddo  
      idmgx=(jjxt-1)/mxm
      jxmg=jjxt+2-idmgx*mxm

      jjztp=mzt/2+1
      idmgzp=(jjztp-1)/mzm
      jzmgp=jjztp+2-idmgzp*mzm

      jjztm=mzt/2
      idmgzm=(jjztm-1)/mzm
      jzmgm=jjztm+2-idmgzm*mzm

      nrank_mgp=nranky*nprxz+idmgzp*nprx+idmgx
      nrank_mgm=nranky*nprxz+idmgzm*nprx+idmgx

      return
      end

!ws********************************************
      subroutine calculate_ps1_mgax
      USE DECLARE
      include 'mpif.h'
      real*8,dimension(mx,2,my) :: wsz2,wsz2s,wsz2r
      real*8,dimension(2,mz,my) :: wsx2,wsx2s,wsx2r
      integer status(mpi_status_size)
      integer kx,kz,nranksend_zp,nrankrecv_zp,nranksend_zm,nrankrecv_zm,nranksendx
!   ps1(xmg,0)=0
!      write(*,*) nrank,nrank_mgp,nrank_mgm,'jmg=',jxmg,jzmgp,jzmgm,idmgx,idmgzp,idmgzm

      if(nrank .eq. nrank_mgp) then
      do jy=iy_first,iy_last
      ps1(jxmg,jzmgp,jy)=-xx(jxmg)*(x1(jxmg,jzmgp-1,jy,8)+x1(jxmg,jzmgp,jy,8))/4*(xx(jxmg)-xmg) &
                         +xx(jxmg)*(x1(jxmg,jzmgp-1,jy,6)+3*x1(jxmg,jzmgp,jy,6))/4*(zz(jzmgp)-zmg)
      do jz=jzmgp,iz_last-2
          ps1(jxmg,jz+1,jy)=ps1(jxmg,jz,jy)+xx(jxmg)*(x1(jxmg,jz,jy,6)+x1(jxmg,jz+1,jy,6))/2*(zz(jz+1)-zz(jz))
          do jx=jxmg,ix_last-2
          ps1(jx+1,jz,jy)=ps1(jx,jz,jy)-(xx(jx)*x1(jx,jz,jy,8)+xx(jx+1)*x1(jx+1,jz,jy,8))/2*(xx(jx+1)-xx(jx))
          enddo
          do jx=jxmg,ix_first+2,-1
          ps1(jx-1,jz,jy)=ps1(jx,jz,jy)-(xx(jx)*x1(jx,jz,jy,8)+xx(jx-1)*x1(jx-1,jz,jy,8))/2*(xx(jx-1)-xx(jx))
          enddo

      enddo
      enddo
!      write(*,*) 'ws1p',nrank

      endif

      if(nrank .eq. nrank_mgm) then
      do jy=iy_first,iy_last
      ps1(jxmg,jzmgm,jy)=-xx(jxmg)*(x1(jxmg,jzmgm+1,jy,8)+x1(jxmg,jzmgm,jy,8))/4*(xx(jxmg)-xmg) &
                         +xx(jxmg)*(x1(jxmg,jzmgm+1,jy,6)+3*x1(jxmg,jzmgm,jy,6))/4*(zz(jzmgm)-zmg)
      do jz=jzmgm,iz_first+2,-1
      ps1(jxmg,jz-1,jy)=ps1(jxmg,jz,jy)+xx(jxmg)*(x1(jxmg,jz,jy,6)+x1(jxmg,jz-1,jy,6))/2*(zz(jz-1)-zz(jz))
      do jx=jxmg,ix_last-2
      ps1(jx+1,jz,jy)=ps1(jx,jz,jy)-(xx(jx)*x1(jx,jz,jy,8)+xx(jx+1)*x1(jx+1,jz,jy,8))/2*(xx(jx+1)-xx(jx))
      enddo
      do jx=jxmg,ix_first+2,-1
      ps1(jx-1,jz,jy)=ps1(jx,jz,jy)-(xx(jx)*x1(jx,jz,jy,8)+xx(jx-1)*x1(jx-1,jz,jy,8))/2*(xx(jx-1)-xx(jx))
      enddo

      enddo

      enddo
!      write(*,*) 'ws1m',nrank

      endif

!      call MPI_Barrier(MPI_COMM_WORLD,ierror)

      do kz=0,nprz/2-2
!      nranksend_zp=nrank_mgp+kz*nprx
!      nrankrecv_zp=nrank_mgp+(kz+1)*nprx
!      nranksend_zm=nrank_mgm-kz*nprx
!      nrankrecv_zm=nrank_mgm-(kz+1)*nprx
      
      if(nrank==nrank_mgp+kz*nprx) then
!     write(*,*) 'ws11p',nrank,kz
      wsz2(:,:,:)=ps1(:,iz_last-3:iz_last-2,1:my)
      CALL MPI_Send( wsz2,2*mx*my, MPI_DOUBLE_PRECISION,nrank+nprx, kz,  &
		      MPI_COMM_WORLD,ierror )
!        write(*,*) 'ws11ps',nrank,kz
      endif
!      CALL MPI_SendRecv(wsz2s,mx*2*my, MPI_DOUBLE_PRECISION, nrankrecv_zp, kz,  &
!                        wsz2r,mx*2*my, MPI_DOUBLE_PRECISION, nranksend_zp, kz,  &
!		               MPI_COMM_WORLD,status,ierror )

      if(nrank==nrank_mgm-kz*nprx) then
!      write(*,*) 'ws11m',nrank,kz
      wsz2(:,:,:)=ps1(:,iz_first+2:iz_first+3,1:my)
      CALL MPI_Send( wsz2,2*mx*my, MPI_DOUBLE_PRECISION,nrank-nprx, kz,  &
		      MPI_COMM_WORLD,ierror )
!      write(*,*) 'ws11ms',nrank,kz
      endif
!      CALL MPI_SendRecv(wsz2s,mx*2*my, MPI_DOUBLE_PRECISION, nrankrecv_zm, kz,  &
!                        wsz2r,mx*2*my, MPI_DOUBLE_PRECISION, nranksend_zm, kz,  &
!		               MPI_COMM_WORLD,status,ierror )
      if(nrank==nrank_mgp+(kz+1)*nprx) then
!      write(*,*) 'ws42',nrank,kz
      CALL MPI_Recv( wsz2 ,2*mx*my, MPI_DOUBLE_PRECISION,nrank-nprx, kz,  &
		      MPI_COMM_WORLD, status,ierror )
        ps1(:,iz_first:iz_first+1,1:my)=wsz2(:,:,:)
!      write(*,*) 'ws43',nrank,kz
      do jz=iz_first+2,iz_last-2
          ps1(jxmg,jz,:)=ps1(jxmg,jz-1,:)+xx(jxmg)*(x1(jxmg,jz-1,:,6)+x1(jxmg,jz,:,6))/2*(zz(jz)-zz(jz-1))
          do jx=jxmg,ix_last-2
          ps1(jx+1,jz,:)=ps1(jx,jz,:)-(xx(jx)*x1(jx,jz,:,8)+xx(jx+1)*x1(jx+1,jz,:,8))/2*(xx(jx+1)-xx(jx))
          enddo
          do jx=jxmg,ix_first+2,-1
          ps1(jx-1,jz,:)=ps1(jx,jz,:)-(xx(jx)*x1(jx,jz,:,8)+xx(jx-1)*x1(jx-1,jz,:,8))/2*(xx(jx-1)-xx(jx))
          enddo

      enddo
!      write(*,*) 'ws4',nrank,kz
      
      endif


      if(nrank==nrank_mgm-(kz+1)*nprx) then
      CALL MPI_Recv(wsz2,mx*2*my, MPI_DOUBLE_PRECISION,nrank+nprx, kz,  &
		      MPI_COMM_WORLD, status,ierror )
      ps1(:,iz_last-1:iz_last,1:my)=wsz2(:,:,:)
      do jz=iz_last-2,iz_first+2,-1
          ps1(jxmg,jz,:)=ps1(jxmg,jz+1,:)+xx(jxmg)*(x1(jxmg,jz+1,:,6)+x1(jxmg,jz,:,6))/2*(zz(jz)-zz(jz+1))
          do jx=jxmg,ix_last-2
          ps1(jx+1,jz,:)=ps1(jx,jz,:)-(xx(jx)*x1(jx,jz,:,8)+xx(jx+1)*x1(jx+1,jz,:,8))/2*(xx(jx+1)-xx(jx))
          enddo
          do jx=jxmg,ix_first+2,-1
          ps1(jx-1,jz,:)=ps1(jx,jz,:)-(xx(jx)*x1(jx,jz,:,8)+xx(jx-1)*x1(jx-1,jz,:,8))/2*(xx(jx-1)-xx(jx))
          enddo

      enddo
      
      endif
 !     call MPI_Barrier(MPI_COMM_WORLD,ierror)
      enddo

 !     call MPI_Barrier(MPI_COMM_WORLD,ierror)

      if(nrkx(nrank)==idmgx) then
!      wsx2s(:,:,:)=ps1(ix_last-3:ix_last-2,:,1:my)
      CALL MPI_Send( ps1(ix_last-3:ix_last-2,:,1:my),2*mz*my, MPI_DOUBLE_PRECISION,nrank+1, 1,  &
		      MPI_COMM_WORLD,ierror )
      CALL MPI_Send( ps1(ix_first+2:ix_first+3,:,1:my),2*mz*my, MPI_DOUBLE_PRECISION,nrank-1, 1,  &
		      MPI_COMM_WORLD,ierror )
      endif

      do kx=1,idmgx
!      if(nrank==nrank_mgp+kx) wsx2s(:,:,:)=ps1(ix_last-3:ix_last-2,:,1:my)
!      CALL MPI_SendRecv(wsx2s,2*mz*my, MPI_DOUBLE_PRECISION, nrank_mgp+kx+1, kx,  &
!                        wsx2r,2*mz*my, MPI_DOUBLE_PRECISION, nrank_mgp+kx, kx,  &
!		               MPI_COMM_WORLD,status,ierror )

      if(nrkx(nrank)==idmgx+kx) then
          CALL MPI_Recv( wsx2r,2*mz*my, MPI_DOUBLE_PRECISION,nrank-1, kx,  &
		          MPI_COMM_WORLD, status,ierror )
          ps1(ix_first:ix_first+1,:,1:my)=wsx2r(:,:,:)
          do jx=ix_first,ix_last-2
          ps1(jx+1,:,:)=ps1(jx,:,:)-(xx(jx)*x1(jx,:,:,8)+xx(jx+1)*x1(jx+1,:,:,8))/2*(xx(jx+1)-xx(jx))
          enddo
 !         write(*,*) 'ws2',nrank,kx
          if(nrkx(nrank).lt. nprx-1) then 
          CALL MPI_Send( ps1(ix_last-3:ix_last-2,:,1:my),2*mz*my, MPI_DOUBLE_PRECISION,nrank+1, kx+1,  &
		          MPI_COMM_WORLD,ierror )
          endif
      endif

      if(nrkx(nrank)==idmgx-kx) then
          CALL MPI_Recv( wsx2r,2*mz*my, MPI_DOUBLE_PRECISION,nrank+1, kx,  &
		          MPI_COMM_WORLD, status,ierror )
         ps1(ix_last-1:ix_last,:,1:my)=wsx2r(:,:,:)
          do jx=ix_last,ix_first+2,-1
          ps1(jx-1,:,:)=ps1(jx,:,:)-(xx(jx)*x1(jx,:,:,8)+xx(jx-1)*x1(jx-1,:,:,8))/2*(xx(jx-1)-xx(jx))
          enddo
 !         write(*,*) 'ws3',nrank,kx
          if(nrkx(nrank).gt. 0) then 
          CALL MPI_Send( ps1(ix_first+2:ix_first+3,:,1:my),2*mz*my, MPI_DOUBLE_PRECISION,nrank-1, kx+1,  &
		          MPI_COMM_WORLD,ierror )
          endif
      endif

      enddo

      call smthxzy_dis_v2(ps1,1,3)

      return
      end
!ws********************************************
      subroutine estimate_pst1
      USE DECLARE
      include 'mpif.h'
      real*8, dimension(mxt,mzt) :: pst1
      real*8 psmg,pstgx,bxgx,bzgx
      integer jjxt,jjztp,jjztm
!   ps1(xmg,0)=0
      jjxt=1
      do while(xxt(jjxt) .lt. xmg)
      jjxt=jjxt+1
      enddo  
      psmg=0
      bzgx=(bz(jjxt,mzt/2)+bz(jjxt,mzt/2+1))/2
      bxgx=(bx(jjxt,mzt/2)+bx(jjxt,mzt/2+1))/2
      pstgx=psmg-(xxt(jjxt)*bzgx)/2*(xxt(jjxt)-xmg)

      jjztp=mzt/2+1
      pst1(jjxt,jjztp)=pstgx+xxt(jjxt)*(bxgx+bx(jjxt,jjztp))/2*(zzt(jjztp)-0)

      do jz=jjztp,mzt-1

      pst1(jjxt,jz+1)=pst1(jjxt,jz)+xxt(jjxt)*(bx(jjxt,jz)+bx(jjxt,jz+1))/2*(zzt(jz+1)-zzt(jz))
      do jx=jjxt,mxt-1
      pst1(jx+1,jz)=pst1(jx,jz)-(xxt(jx)*bz(jx,jz)+xxt(jx+1)*bz(jx+1,jz))/2*(xxt(jx+1)-xxt(jx))
      enddo
      do jx=jjxt,2,-1
      pst1(jx-1,jz)=pst1(jx,jz)-(xxt(jx)*bz(jx,jz)+xxt(jx-1)*bz(jx-1,jz))/2*(xxt(jx-1)-xxt(jx))
      enddo

      enddo

      jjztm=mzt/2
      pst1(jjxt,jjztm)=pstgx+xxt(jjxt)*(bxgx+bx(jjxt,jjztm))/2*(zzt(jjztm)-0)

      do jz=jjztm,2,-1

      pst1(jjxt,jz-1)=pst1(jjxt,jz)+xxt(jjxt)*(bx(jjxt,jz)+bx(jjxt,jz-1))/2*(zzt(jz-1)-zzt(jz))
      do jx=jjxt,mxt-1
      pst1(jx+1,jz)=pst1(jx,jz)-(xxt(jx)*bz(jx,jz)+xxt(jx+1)*bz(jx+1,jz))/2*(xxt(jx+1)-xxt(jx))
      enddo
      do jx=jjxt,2,-1
      pst1(jx-1,jz)=pst1(jx,jz)-(xxt(jx)*bz(jx,jz)+xxt(jx-1)*bz(jx-1,jz))/2*(xxt(jx-1)-xxt(jx))
      enddo

      enddo
      if(nrank .eq. 0) then
      open(unit=905,file='pst1.dat',status='unknown',form='formatted')
      write(905,100)(((pst(jx,jz),pst1(jx,jz)),jx=1,mxt),jz=1,mzt)
 100  format(2(1x,e12.5))
      close(905)
      endif
      
      return
      end
!ws************************************************************************
      subroutine recrd_ps1
      USE DECLARE
      include 'mpif.h'
!
      character*9 output
      character*3 cn
      character*3 cn1
      
      output='ps'//cn1(nrank)//cn(nst)
      open(unit=72,file=output,status='unknown',form='unformatted')
      write(72)ncase,nstep,time
      write(72)ps1
      close(72)
    
      return
      end

!ws************************************************************************
      subroutine recrd_cud
      USE DECLARE
      include 'mpif.h'
!
      character*9 output
      character*3 cn
      character*3 cn1
      
      output='cd'//cn1(nrank)//cn(nst)
      open(unit=73,file=output,status='unknown',form='unformatted')
      write(73)ncase,nstep,time
      write(73)cud
      close(73)
    
      return
      end
      subroutine readin_cud
      USE DECLARE
      include 'mpif.h'
!
      character*9 output
      character*3 cn
      character*3 cn1

      output='cd'//cn1(nrank)//cn(nst)
      open(unit=73,file=output,status='unknown',form='unformatted')
      read(73)ncase,nstep,time
      read(73)cud
      
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
      if(psi(jx,jz) .ge. psia) then
      cud(jx,jz,:,:)=0    
      endif
      enddo
      enddo
      close(73)     
     
      return
      end
    
      subroutine find_cud_OX
      USE DECLARE
      USE DECLARE_OXpoint      
          call find_maxmin(cud(:,:,3,2),0,wmax,wmin,jxto,jzto,jxtx,jztx)
      yy_O= 0 
      ps_O= pst(jxto,jzto)
      rr_O= rrt(jxto,jzto) 
      tc_O= 0 

      yy_X= 0 !yy_O
      ps_X= psmode
      rr_X= rrmode
      tc_X= pi/2 
      return
      end
          
!ws*******************************************************************************
     subroutine readmap_wch
     USE DECLARE 
      integer, parameter :: nq = 33
      integer, parameter :: nr = 85 !int( sqrt(ndat/3) )
      integer, parameter :: nw = 39
    !
      integer lcell(nr,nr)
      integer lnext(ndat),lnext12(ndat12),lnext34(ndat34)
      real*8 risq(ndat),risq12(ndat12),risq34(ndat34)
      real*8 tch_NOVA(ndat),tch12_NOVA(ndat12),tch34_NOVA(ndat34),tcdx_NOVA(ndat),tcdz_NOVA(ndat),aj_NOVA(ndat)
      real*8 aw(5,ndat),aw12(5,ndat12),aw34(5,ndat34)
      real*8 rimax,ximin,zimin,dxi,dzi,wwdx,wwdz
      real*8, dimension(n2th+5,npsi):: tcst,tcdxst,tcdzst,ajst
       real*8, dimension(mxt,mzt) :: tchdx,tchdz,tcht_dx,tcht_dz,ajt,ajt_dx,ajt_dz
      integer ipi
      include 'mpif.h' 

      do jt=1,n2th+5
      thst(jt)=pi*(jt-3)/(nthe-1)
      enddo
      ipi=nthe+2
!      open(888,file='psi_xz.dat')
!      do 20 jj=1,npsip
!      do 30 ij=1,nthe3
!      read(888,1000) j,i,psst(i,j),xxst(i,j),zzst(i,j) 
!   30 continue
!   20 continue


      open(777,file='wch.dat')
      do j=2,npsi
      do i=3,nthe+1
      read(777,"(4(1x,e17.9))") tcst(i,j),tcdxst(i,j),tcdzst(i,j),ajst(i,j)
      
      jd=(j-2)*(n2th-1)+i-2   
      tch_NOVA(jd)=tcst(i,j)
      tcdx_NOVA(jd)=tcdxst(i,j)
      tcdz_NOVA(jd)=tcdzst(i,j)
      aj_NOVA(jd)=ajst(i,j)

      if(i.gt.3) then
      im=2*nthe+2-(i-2)
      tcst(im,j)=2*pi-tcst(i,j)
      tcdxst(im,j)=-tcdxst(i,j)
      tcdzst(im,j)=tcdzst(i,j)
      ajst(im,j)=ajst(i,j)

      jdm=(j-2)*(n2th-1)+im-3
      tch_NOVA(jdm)=tcst(im,j)
      tcdx_NOVA(jdm)=tcdxst(im,j)
      tcdz_NOVA(jdm)=tcdzst(im,j)
      aj_NOVA(jdm)=ajst(im,j)

      endif       
      enddo

      tcst(1,j)=tcst(1+n2th,j)-2*pi
      tcdxst(1,j)=tcdxst(1+n2th,j)
      tcdzst(1,j)=tcdzst(1+n2th,j)
      ajst(1,j)=ajst(1+n2th,j)

      tcst(2,j)=tcst(2+n2th,j)-2*pi
      tcdxst(2,j)=tcdxst(2+n2th,j)
      tcdzst(2,j)=tcdzst(2+n2th,j)
      ajst(2,j)=ajst(2+n2th,j)

      tcst(3+n2th,j)=tcst(3,j)+2*pi
      tcdxst(3+n2th,j)=tcdxst(3,j)
      tcdzst(3+n2th,j)=tcdzst(3,j)
      ajst(3+n2th,j)=ajst(3,j)

      tcst(4+n2th,j)=tcst(4,j)+2*pi
      tcdxst(4+n2th,j)=tcdxst(4,j)
      tcdzst(4+n2th,j)=tcdzst(4,j)
      ajst(4+n2th,j)=ajst(4,j)

      tcst(5+n2th,j)=tcst(5,j)+2*pi
      tcdxst(5+n2th,j)=tcdxst(5,j)
      tcdzst(5+n2th,j)=tcdzst(5,j)
      ajst(5+n2th,j)=ajst(5,j)

      tcst(ipi,j)=pi
      call interp1d3l(ajst(ipi-2,j),ajst(ipi-1,j),ajst(ipi+1,j),ajst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),ajst(ipi,j))
      call interp1d3l(tcdxst(ipi-2,j),tcdxst(ipi-1,j),tcdxst(ipi+1,j),tcdxst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),tcdxst(ipi,j))
      call interp1d3l(tcdzst(ipi-2,j),tcdzst(ipi-1,j),tcdzst(ipi+1,j),tcdzst(ipi+2,j), &
                      zzst(ipi-2,j),zzst(ipi-1,j),zzst(ipi+1,j),zzst(ipi+2,j),zzst(ipi,j),tcdzst(ipi,j))
      enddo
!      thst(ipi)=pi

      close(777)

      tcst(:,1)=tcst(:,2)

      call interp1d3l(ajst(ipi,3),ajst(ipi,2),ajst(3,2),ajst(3,3), &
                      xxst(ipi,3),xxst(ipi,2),xxst(3,2),xxst(3,3),xmg,ajst(3,1))
      call interp1d3l(tcdxst(ipi,3),tcdxst(ipi,2),tcdxst(3,2),tcdxst(3,3), &
                      xxst(ipi,3),xxst(ipi,2),xxst(3,2),xxst(3,3),xmg,tcdxst(3,1))
      call interp1d3l(tcdzst(ipi,3),tcdzst(ipi,2),tcdzst(3,2),tcdzst(3,3), &
                      xxst(ipi,3),xxst(ipi,2),xxst(3,2),xxst(3,3),xmg,tcdzst(3,1))
!      call cubic(x1,x2,x3,x4,y1,y2,y3,y4,y,ans,dans,is,ierr)
      ajst(:,1)=ajst(3,1)
      tcdxst(:,1)=tcdxst(3,1)
      tcdzst(:,1)=tcdzst(3,1)
      
      do i=1,ipi+2
      do j=2,npsi
      jd=(i-1)*(npsi-1)+j-1
      tch12_NOVA(jd)=tcst(i,j)
      tch34_NOVA(jd)=tcst(i+ipi-3,j)
      enddo
      enddo

      call qshep2 ( ndat12, xx12_NOVA, zz12_NOVA, tch12_NOVA, nq, nw, nr, lcell, lnext12, ximin, zimin, &
        dxi, dzi, rimax, risq12, aw12, ier )
!     write(*,*) ier
      do jz=mzt/2+1,mzt
      do jx=1,mxt
      call qs2grd ( xxt(jx), zzt(jz), ndat12, xx12_NOVA, zz12_NOVA, tch12_NOVA, nr, lcell, lnext12, ximin, &
        zimin, dxi, dzi, rimax, risq12, aw12, tcht(jx,jz), tcht_dx(jx,jz), tcht_dz(jx,jz), ier )
      enddo
      enddo  
!      
      call qshep2 ( ndat34, xx34_NOVA, zz34_NOVA, tch34_NOVA, nq, nw, nr, lcell, lnext34, ximin, zimin, &
        dxi, dzi, rimax, risq34, aw34, ier )
!     write(*,*) ier
      do jz=1,mzt/2
      do jx=1,mxt
      call qs2grd ( xxt(jx), zzt(jz), ndat34, xx34_NOVA, zz34_NOVA, tch34_NOVA, nr, lcell, lnext34, ximin, &
        zimin, dxi, dzi, rimax, risq34, aw34, tcht(jx,jz), tcht_dx(jx,jz), tcht_dz(jx,jz), ier )
      enddo
      enddo

      call qshep2 ( ndat, xx_NOVA, zz_NOVA, tcdx_NOVA, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_NOVA, zz_NOVA, tcdx_NOVA, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, tchdx(jx,jz), wwdx, wwdz, ier )
      enddo
      enddo
      call qshep2 ( ndat, xx_NOVA, zz_NOVA, tcdz_NOVA, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_NOVA, zz_NOVA, tcdz_NOVA, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, tchdz(jx,jz), wwdx, wwdz, ier )
      enddo
      enddo

      call qshep2 ( ndat, xx_NOVA, zz_NOVA, aj_NOVA, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_NOVA, zz_NOVA, aj_NOVA, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, ajt(jx,jz), ajt_dx(jx,jz), ajt_dz(jx,jz), ier )
      enddo
      enddo
      
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
      tcxz(jx,jz)=tcht(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
      enddo
      enddo
      
      if(nrank.eq.0) then

      open(791,file='wtcst.dat')
      write(791,2100) (((tcst(i,j),tcdxst(i,j),tcdzst(i,j),ajst(i,j)),i=1,n2th+5),j=1,npsi)
 2100 format(4(1x,e17.9))
      close(791)

      open(unit=792,file='wtcxz.dat',status='unknown',form='formatted')
      write(792,300)(((tcht(jx,jz),tchdx(jx,jz),tchdz(jx,jz)),jx=1,mxt),jz=1,mzt)
 300  format(3(1x,e12.5))
      close(792)

      open(unit=793,file='ajxz.dat',status='unknown',form='formatted')
      write(793,300)(((ajt(jx,jz),ajt_dx(jx,jz),ajt_dz(jx,jz)),jx=1,mxt),jz=1,mzt)
      close(793)
      endif 
       
       
      return
      end

!ws*******************************************************************************

!ws++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!ws pllconduct
!ws++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine pllconduct(ipll)
      integer ipll

      select case(ipll)
      case(1)
      call pll_subloop(nploop)
      case(2)
      call pll_smthpline(nploop)
      case(3)
      call pll_soundwave(nploop)
      case(4)
      !call pll_PETSC
      case(5)
     ! call pll_PETSC_T1
!      case default
      end select

      return
      end

      !ws*****************************************************
      subroutine pll_subloop(npll)
      USE DECLARE
      integer npll
      select case(lscheme)
      case(1) !Euler
      call pll_subloop_EL(npll)
      case(2) !RK4
      call pll_subloop_RK(npll)
      case(4)
      call pllconduct_implicity
!      case default
      end select
      return
      end

!ws*****************************************************
      subroutine pll_subloop_EL(npll)
      USE DECLARE
      integer npll
      real*8, dimension(mx,mz,my) :: p0dif,p1dif
      real*8 betx,betz,betcut,dbet

      betcut=0.2
      dbet=0.02
      do jy=iy_first,iy_last
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last  
      betx=abs(x(jx,jz,jy,6)/x(jx,jz,jy,7)*xx(jx)*dy(jy)/dx(jx))
      betz=abs(x(jx,jz,jy,8)/x(jx,jz,jy,7)*xx(jx)*dy(jy)/dz(jz))
      fcx(jx,jz,jy)=(1+tanh((betx-betcut)/dbet))/2.0
      fcz(jx,jz,jy)=(1+tanh((betz-betcut)/dbet))/2.0
      enddo
      enddo
      enddo  

      call right_pll_p0(p0dif)
      
      do k=1,npll
      call right_pll_p1(p1dif) 
!      call right_pll(pdif) 
      x(:,:,:,2)=x(:,:,:,2)+p1dif(:,:,:)*dt+p0dif(:,:,:)*dt
      enddo

      return
      end

!ws*****************************************************
      subroutine pll_subloop_RK(npll)
      USE DECLARE
      integer npll
      real*8, dimension(mx,mz,my) :: pdif,pfold,pm
      do k=1,npll
      call right_pll_p1(pdif) 
        pfold(:,:,:)=x(:,:,:,2)
        pm(:,:,:)=pfold(:,:,:)+pdif(:,:,:)*dt/6.
        x(:,:,:,2)=pfold(:,:,:)+pdif(:,:,:)*dt/2.
!
      call right_pll_p1(pdif) 
        pm(:,:,:)=pm(:,:,:)+pdif(:,:,:)*dt/3.
        x(:,:,:,2)=pfold(:,:,:)+pdif(:,:,:)*dt/2.
!
      call right_pll_p1(pdif)
        pm(:,:,:)=pm(:,:,:)+pdif(:,:,:)*dt/3.
        x(:,:,:,2)=pfold(:,:,:)+pdif(:,:,:)*dt
!
      call right_pll_p1(pdif)
        x(:,:,:,2)=pm(:,:,:)+pdif(:,:,:)*dt/6.

!      call bndry_p
      enddo

      return
      end
!ws*****************************************************
      subroutine pll_subloop_y(npll)
      USE DECLARE
      integer npll
      real*8, dimension(mx,mz,my) :: p1dif
      
      do k=1,npll
      call right_py_p1(p1dif) 
!      call right_pll(pdif) 
      x(:,:,:,2)=x(:,:,:,2)+p1dif(:,:,:)*dt
      enddo

      return
      end


!ws*****************************************************
      subroutine bndry_p
      USE DECLARE
      do jy=1,my
      x1(:,:,jy,2)=x(:,:,jy,2)-xint(:,:,2)
      enddo
      call valbm_atlastgrid_v1(x1(:,:,:,2),1,1)
      call mpi_transfersm(x1(:,:,:,2),1) 
      do jy=1,my
      x(:,:,jy,2)=x1(:,:,jy,2)+xint(:,:,2)
      enddo
      return
      end
!ws*****************************************************
      subroutine right_pll(udif)
      USE DECLARE
      real*8, dimension(mx,mz,my,3) :: hf
      real*8, dimension(mx,mz,my) :: udif
      real*8 bfx,bfz,bfy,bb2,udx,udz,udy,bgrdu,hfb

!  d1f2= d f / dx  with second-order accuracy central difference
      d1f2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(fp1-f0) &
         -(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)
!  d2f2= d2f / dx2  with second-order accuracy central difference
      d2f2(fm1,f0,fp1,xm1,x0,xp1)= &
       2*((fp1-f0)/(xp1-x0)-(f0-fm1)/(x0-xm1))/(xp1-xm1)
!  d1xf2= d Rf / dx  with second-order accuracy central difference
      d1xf2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(xp1*fp1-x0*f0) &
         -(xp1-x0)/(xm1-x0)*(xm1*fm1-x0*f0))/(xm1-xp1)

      call bndry_p
 !     u(:,:,:)=x(:,:,:,2)

      do jy=iy_first+1,iy_last-1
      do jz=iz_first+1,iz_last-1
      do jx=ix_first+1,ix_last-1
      bfx=x(jx,jz,jy,6)
      bfz=x(jx,jz,jy,8)
      bfy=x(jx,jz,jy,7)
      bb2=bfx**2+bfz**2+bfx**2

      udx = d1f2(x(jx-1,jz,jy,2),x(jx,jz,jy,2),x(jx+1,jz,jy,2),xx(jx-1),xx(jx),xx(jx+1))
      udz = d1f2(x(jx,jz-1,jy,2),x(jx,jz,jy,2),x(jx,jz+1,jy,2),zz(jz-1),zz(jz),zz(jz+1))
      udy = d1f2(x(jx,jz,jy-1,2),x(jx,jz,jy,2),x(jx,jz,jy+1,2),yy(jy-1),yy(jy),yy(jy+1))
      bgrdu=bfx*udx+bfz*udz+bfy*udy/xx(jx)

      hfb=-kap_ll(jx,jz)*bgrdu/bb2
      hf(jx,jz,jy,1)=hfb*bfx
      hf(jx,jz,jy,2)=hfb*bfy
      hf(jx,jz,jy,3)=hfb*bfz
      enddo
      enddo
      enddo

      do jy=iy_first+2,iy_last-2
      do jz=iz_first+2,iz_last-2
      do jx=ix_first+2,ix_last-2  
      if(psi(jx,jz).lt.psia1) then    
      udif(jx,jz,jy)=-d1xf2(hf(jx-1,jz,jy,1),hf(jx,jz,jy,1),hf(jx+1,jz,jy,1),xx(jx-1),xx(jx),xx(jx+1))/xx(jx) &
           -d1f2(hf(jx,jz-1,jy,3),hf(jx,jz,jy,3),hf(jx,jz+1,jy,3),zz(jz-1),zz(jz),zz(jz+1)) &
           -d1f2(hf(jx,jz,jy-1,2),hf(jx,jz,jy,2),hf(jx,jz,jy+1,2),yy(jy-1),yy(jy),yy(jy+1))/xx(jx)
      endif
      enddo
      enddo
      enddo
      return
      end
!ws*****************************************************
      subroutine right_pll_p0(udif)
      USE DECLARE
      real*8, dimension(mx,mz,my,3) :: hf
      real*8, dimension(mx,mz,my) :: udif
      real*8 bfx,bfz,bfy,bb2,kpl,udx,udz,udy,bgrdu,hfb

!  d1f2= d f / dx  with second-order accuracy central difference
      d1f2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(fp1-f0) &
         -(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)
!  d2f2= d2f / dx2  with second-order accuracy central difference
      d2f2(fm1,f0,fp1,xm1,x0,xp1)= &
       2*((fp1-f0)/(xp1-x0)-(f0-fm1)/(x0-xm1))/(xp1-xm1)
!  d1xf2= d Rf / dx  with second-order accuracy central difference
      d1xf2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(xp1*fp1-x0*f0) &
         -(xp1-x0)/(xm1-x0)*(xm1*fm1-x0*f0))/(xm1-xp1)

!      call bndry_p
!      u(:,:,:)=x1(:,:,:,2)
      do jy=iy_first+1,iy_last-1
      do jz=iz_first+1,iz_last-1
      do jx=ix_first+1,ix_last-1
      bfx=x(jx,jz,jy,6)
      bfz=x(jx,jz,jy,8)
      bfy=x(jx,jz,jy,7)
      bb2=bfx**2+bfz**2+bfx**2

      bgrdu=x1(jx,jz,jy,6)*xint_dx(jx,jz,2)+x1(jx,jz,jy,8)*xint_dz(jx,jz,2)

!      hfb=-kap_ll(jx,jz)*bgrdu/bb2
      hfb=-kap_ll(jx,jz)*bgrdu/bb2
      hf(jx,jz,jy,1)=hfb*bfx
      hf(jx,jz,jy,2)=hfb*bfy
      hf(jx,jz,jy,3)=hfb*bfz
      enddo
      enddo
      enddo

      do jy=iy_first+2,iy_last-2
      do jz=iz_first+2,iz_last-2
      do jx=ix_first+2,ix_last-2 
      if(psi(jx,jz).lt.psia1) then     
      udif(jx,jz,jy)=(-fcx(jx,jz,jy)*d1xf2(hf(jx-1,jz,jy,1),hf(jx,jz,jy,1),hf(jx+1,jz,jy,1),xx(jx-1),xx(jx),xx(jx+1))/xx(jx) &
           -fcz(jx,jz,jy)*d1f2(hf(jx,jz-1,jy,3),hf(jx,jz,jy,3),hf(jx,jz+1,jy,3),zz(jz-1),zz(jz),zz(jz+1)) &
           -d1f2(hf(jx,jz,jy-1,2),hf(jx,jz,jy,2),hf(jx,jz,jy+1,2),yy(jy-1),yy(jy),yy(jy+1))/xx(jx))
      endif
      enddo
      enddo
      enddo
      return
      end
!ws*****************************************************
      subroutine right_pll_p1(udif)
      USE DECLARE
      real*8, dimension(mx,mz,my,3) :: hf
      real*8, dimension(mx,mz,my) :: udif
      real*8 bfx,bfz,bfy,bb2,kpl,udx,udz,udy,bgrdu,hfb

!  d1f2= d f / dx  with second-order accuracy central difference
      d1f2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(fp1-f0) &
         -(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)
!  d2f2= d2f / dx2  with second-order accuracy central difference
      d2f2(fm1,f0,fp1,xm1,x0,xp1)= &
       2*((fp1-f0)/(xp1-x0)-(f0-fm1)/(x0-xm1))/(xp1-xm1)
!  d1xf2= d Rf / dx  with second-order accuracy central difference
      d1xf2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(xp1*fp1-x0*f0) &
         -(xp1-x0)/(xm1-x0)*(xm1*fm1-x0*f0))/(xm1-xp1)

      call bndry_p
!      u(:,:,:)=x1(:,:,:,2)
      do jy=iy_first+1,iy_last-1
      do jz=iz_first+1,iz_last-1
      do jx=ix_first+1,ix_last-1
      bfx=x(jx,jz,jy,6)
      bfz=x(jx,jz,jy,8)
      bfy=x(jx,jz,jy,7)
      bb2=bfx**2+bfz**2+bfx**2

      udx = d1f2(x1(jx-1,jz,jy,2),x1(jx,jz,jy,2),x1(jx+1,jz,jy,2),xx(jx-1),xx(jx),xx(jx+1))
      udz = d1f2(x1(jx,jz-1,jy,2),x1(jx,jz,jy,2),x1(jx,jz+1,jy,2),zz(jz-1),zz(jz),zz(jz+1))
      udy = d1f2(x1(jx,jz,jy-1,2),x1(jx,jz,jy,2),x1(jx,jz,jy+1,2),yy(jy-1),yy(jy),yy(jy+1))
      bgrdu=fcx(jx,jz,jy)*bfx*udx+fcz(jx,jz,jy)*bfz*udz+bfy*udy/xx(jx) ! &
!           +x1(jx,jz,jy,6)*xint_dx(jx,jz,2)+x1(jx,jz,jy,8)*xint_dz(jx,jz,2)

!      hfb=-kap_ll(jx,jz)*bgrdu/bb2
      hfb=-kap_ll(jx,jz)*bgrdu/bb2
      hf(jx,jz,jy,1)=hfb*bfx
      hf(jx,jz,jy,2)=hfb*bfy
      hf(jx,jz,jy,3)=hfb*bfz
      enddo
      enddo
      enddo

      do jy=iy_first+2,iy_last-2
      do jz=iz_first+2,iz_last-2
      do jx=ix_first+2,ix_last-2 
      if(psi(jx,jz).lt.psia1) then     
      udif(jx,jz,jy)=(-fcx(jx,jz,jy)*d1xf2(hf(jx-1,jz,jy,1),hf(jx,jz,jy,1),hf(jx+1,jz,jy,1),xx(jx-1),xx(jx),xx(jx+1))/xx(jx) &
           -fcz(jx,jz,jy)*d1f2(hf(jx,jz-1,jy,3),hf(jx,jz,jy,3),hf(jx,jz+1,jy,3),zz(jz-1),zz(jz),zz(jz+1)) &
           -d1f2(hf(jx,jz,jy-1,2),hf(jx,jz,jy,2),hf(jx,jz,jy+1,2),yy(jy-1),yy(jy),yy(jy+1))/xx(jx))
      endif
      enddo
      enddo
      enddo
      return
      end
!ws*****************************************************
      subroutine right_py_p1(udif)
      USE DECLARE
      real*8, dimension(mx,mz,my,3) :: hf
      real*8, dimension(mx,mz,my) :: udif
      real*8 bfx,bfz,bfy,bb2,kpl,udx,udz,udy,bgrdu,hfb

!  d1f2= d f / dx  with second-order accuracy central difference
      d1f2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(fp1-f0) &
         -(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)
!  d2f2= d2f / dx2  with second-order accuracy central difference
      d2f2(fm1,f0,fp1,xm1,x0,xp1)= &
       2*((fp1-f0)/(xp1-x0)-(f0-fm1)/(x0-xm1))/(xp1-xm1)
!  d1xf2= d Rf / dx  with second-order accuracy central difference
      d1xf2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(xp1*fp1-x0*f0) &
         -(xp1-x0)/(xm1-x0)*(xm1*fm1-x0*f0))/(xm1-xp1)

      call bndry_p
!      u(:,:,:)=x1(:,:,:,2)
      do jy=iy_first+1,iy_last-1
      do jz=iz_first+1,iz_last-1
      do jx=ix_first+1,ix_last-1
      bfx=x(jx,jz,jy,6)
      bfz=x(jx,jz,jy,8)
      bfy=x(jx,jz,jy,7)
      bb2=bfx**2+bfz**2+bfx**2

 !     udx = d1f2(x1(jx-1,jz,jy,2),x1(jx,jz,jy,2),x1(jx+1,jz,jy,2),xx(jx-1),xx(jx),xx(jx+1))
 !     udz = d1f2(x1(jx,jz-1,jy,2),x1(jx,jz,jy,2),x1(jx,jz+1,jy,2),zz(jz-1),zz(jz),zz(jz+1))
      udy = d1f2(x1(jx,jz,jy-1,2),x1(jx,jz,jy,2),x1(jx,jz,jy+1,2),yy(jy-1),yy(jy),yy(jy+1))
      bgrdu=bfy*udy/xx(jx) 

      hfb=-kap_ll(jx,jz)*bgrdu/bb2
!      hf(jx,jz,jy,1)=hfb*bfx
      hf(jx,jz,jy,2)=hfb*bfy
!      hf(jx,jz,jy,3)=hfb*bfz
      enddo
      enddo
      enddo

      do jy=iy_first+2,iy_last-2
      do jz=iz_first+2,iz_last-2
      do jx=ix_first+2,ix_last-2 
      if(psi(jx,jz).lt.psia1) then     
      udif(jx,jz,jy)=( &!-d1xf2(hf(jx-1,jz,jy,1),hf(jx,jz,jy,1),hf(jx+1,jz,jy,1),xx(jx-1),xx(jx),xx(jx+1))/xx(jx) &
!           -d1f2(hf(jx,jz-1,jy,3),hf(jx,jz,jy,3),hf(jx,jz+1,jy,3),zz(jz-1),zz(jz),zz(jz+1)) &
           -d1f2(hf(jx,jz,jy-1,2),hf(jx,jz,jy,2),hf(jx,jz,jy+1,2),yy(jy-1),yy(jy),yy(jy+1))/xx(jx))
      endif
      enddo
      enddo
      enddo
      return
      end
!!ws*****************************************************
!      subroutine right_pll_p1_v2(udif)
!      USE DECLARE
!      real*8, dimension(mx,mz,my,3) :: bf
!      real*8, dimension(mx,mz,my) :: u,udif,klrb2
!!  d1f2= d f / dx  with second-order accuracy central difference
!      d1f2(fm1,f0,fp1,xm1,x0,xp1)= &
!        ((xm1-x0)/(xp1-x0)*(fp1-f0) &
!         -(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)
!!  d2f2= d2f / dx2  with second-order accuracy central difference
!      d2f2(fm1,f0,fp1,xm1,x0,xp1)= &
!       2*((fp1-f0)/(xp1-x0)-(f0-fm1)/(x0-xm1))/(xp1-xm1)
!
!      u(:,:,:)=x(:,:,:,2)
!      bf(:,:,:,:)=x(:,:,:,6:8)
!
!      do jy=iy_first,iy_last
!      do jz=iz_first,iz_last
!      do jx=ix_first,ix_last
!      klrb2(jx,jz,jy)=kap_ll(jx,jz)/(bf(jx,jz,jy,1)**2+bf(jx,jz,jy,2)**2+bf(jx,jz,jy,3)**2)
!      enddo
!      enddo
!      enddo
!
!      do jy=iy_first+2,iy_last-2
!      do jz=iz_first+2,iz_last-2
!      do jx=ix_first+2,ix_last-2
!
!      if(psi(jx,jz).lt.psia1) then
!      bfx=x(jx,jz,jy,6)
!      bfy=x(jx,jz,jy,7)
!      bfz=x(jx,jz,jy,8)
!      bb2=bf2(jx,jz,jy)
!      bfx_dx=d1f2(x(jx-1,jz,jy,6),x(jx,jz,jy,6),x(jx+1,jz,jy,6),xx(jx-1),xx(jx),xx(jx+1))
!      bfx_dz=d1f2(x(jx,jz-1,jy,6),x(jx,jz,jy,6),x(jx,jz+1,jy,6),zz(jz-1),zz(jz),zz(jz+1))
!      bfx_dy=d1f2(x(jx,jz,jy-1,6),x(jx,jz,jy,6),x(jx,jz,jy+1,6),yy(jy-1),yy(jy),yy(jy+1))
!      bfz_dx=d1f2(x(jx-1,jz,jy,8),x(jx,jz,jy,8),x(jx+1,jz,jy,8),xx(jx-1),xx(jx),xx(jx+1))
!      bfz_dz=d1f2(x(jx,jz-1,jy,8),x(jx,jz,jy,8),x(jx,jz+1,jy,8),zz(jz-1),zz(jz),zz(jz+1))
!      bfz_dy=d1f2(x(jx,jz,jy-1,8),x(jx,jz,jy,8),x(jx,jz,jy+1,8),yy(jy-1),yy(jy),yy(jy+1))
!      bfy_dx=d1f2(x(jx-1,jz,jy,7),x(jx,jz,jy,7),x(jx+1,jz,jy,7),xx(jx-1),xx(jx),xx(jx+1))
!      bfy_dz=d1f2(x(jx,jz-1,jy,7),x(jx,jz,jy,7),x(jx,jz+1,jy,7),zz(jz-1),zz(jz),zz(jz+1))
!      bfy_dy=d1f2(x(jx,jz,jy-1,7),x(jx,jz,jy,7),x(jx,jz,jy+1,7),yy(jy-1),yy(jy),yy(jy+1))
!
!      bgbx=bfx*xr(jx,jz,jy,6)+bfz*xz(jx,jz,jy,6)+bfy*xy(jx,jz,jy,6)/xx(jx)
!      bgbz=bfx*xr(jx,jz,jy,8)+bfz*xz(jx,jz,jy,8)+bfy*xy(jx,jz,jy,8)/xx(jx)
!      bgby=bfx*(xr(jx,jz,jy,7)-bfy/xx(jx))+bfz*xz(jx,jz,jy,7)+bfy*xy(jx,jz,jy,7)/xx(jx)
!      
!      bgb1x=bfx*x1r(jx,jz,jy,6)+bfz*x1z(jx,jz,jy,6)+bfy*xy(jx,jz,jy,6)/xx(jx)
!      bgb1z=bfx*x1r(jx,jz,jy,8)+bfz*x1z(jx,jz,jy,8)+bfy*xy(jx,jz,jy,8)/xx(jx)
!      bgb1y=bfx*(x1r(jx,jz,jy,7)-bfy/xx(jx))+bfz*x1z(jx,jz,jy,7)+bfy*xy(jx,jz,jy,7)/xx(jx)
!
!!      bgb0x=bfx*xint_dx(jx,jz,6)+bfz*xint_dz(jx,jz,6)
!!      bgb0z=bfx*xint_dx(jx,jz,8)+bfz*xint_dz(jx,jz,8)
!!      bgb0y=bfx*(xint_dx(jx,jz,7)-xint(jx,jz,7)/xx(jx))+bfz*xint_dz(jx,jz,7)
!
!      bdklrb2=1.0/bb2*(bfx*kap_ll_dx(jx,jz)+bfz*kap_ll_dz(jx,jz) &
!             -kap_ll(jx,jz)/bb2*(bfx*bf2dx(jx,jz,jy)+bfz*bf2dx(jx,jz,jy)+bfy*bf2dy(jx,jz,jy)/xx(jx)) )
!
!      udx = d1f2(x1(jx-1,jz,jy,2),x1(jx,jz,jy,2),x1(jx+1,jz,jy,2),xx(jx-1),xx(jx),xx(jx+1))
!      udz = d1f2(x1(jx,jz-1,jy,2),x1(jx,jz,jy,2),x1(jx,jz+1,jy,2),zz(jz-1),zz(jz),zz(jz+1))
!      udy = d1f2(x1(jx,jz,jy-1,2),x1(jx,jz,jy,2),x1(jx,jz,jy+1,2),yy(jy-1),yy(jy),yy(jy+1))
!      ud2xx=d2f2(x1(jx-1,jz,jy,2),x1(jx,jz,jy,2),x1(jx+1,jz,jy,2),xx(jx-1),xx(jx),xx(jx+1))
!      ud2zz=d2f2(x1(jx,jz-1,jy,2),x1(jx,jz,jy,2),x1(jx,jz+1,jy,2),zz(jz-1),zz(jz),zz(jz+1))
!      ud2yy=d2f2(x1(jx,jz,jy-1,2),x1(jx,jz,jy,2),x1(jx,jz,jy+1,2),yy(jy-1),yy(jy),yy(jy+1)) 
!      ud2xz=(x1(jx+1,jz+1,jy,2)+x1(jx-1,jz-1,jy,2)-x1(jx+1,jz-1,jy,2)-x1(jx-1,jz+1,jy,2))/dxz/4
!      ud2zy=(x1(jx,jz+1,jy+1,2)+x1(jx,jz-1,jy-1,2)-x1(jx,jz-1,jy+1,2)-x1(jx,jz+1,jy-1,2))/dzy/4
!      ud2yx=(x1(jx+1,jz,jy+1,2)+x1(jx,jz-1,jy-1,2)-x1(jx+1,jz,jy-1,2)-x1(jx-1,jz,jy+1,2))/dxy/4
!
!      udif(jx,jz,jy)=-kap_ll(jx,jz)/bf2(jx,jz,jy)*( &
!          bfx**2*ud2xx+bfz**2*ud2zz+bfy**2*ud2yy/xx(jx)**2 &
!         +2*bfx*bfz*ud2xz+2*bfz*bfy*ud2zy/xx(jx)+2*bfy*bfx*ud2yx/xx(jx) &
!         +(bfx*bfx_dx+bfz*bfx_dz+bfy*bfx_dy/xx(jx))*udx &
!         +(bfx*bfz_dx+bfz*bfz_dz+bfy*bfz_dy/xx(jx))*udz &
!         +(bfx*(bfy_dx-bfy/xx(jx))+bfz*bfy_dz+bfy*bfy_dy/xx(jx))/xx(jx)*udy) &
!         -bdklrb2*(bfx*udx+bfz*udz+bfy*udy/xx(jx)) 
!      endif
!
!      enddo
!      enddo
!      enddo
!
!      return
!      end
!ws*****************************************************
      subroutine pllconduct_implicity
      USE DECLARE
      real*8, dimension(mx,mz,my,3) :: bf,bfdel
      real*8, dimension(mx,mz,my) :: u,u0,udel,www,klrb2
      real*8, dimension(mx) :: ccx0,ccxp,ccxm,ssx,ufx
      real*8, dimension(mz) :: ccz0,cczp,cczm,ssz,ufz
      real*8, dimension(my) :: ccy0,ccyp,ccym,ssy,ufy
!  d1f2= d f / dx  with second-order accuracy central difference
      d1f2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(fp1-f0) &
         -(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)
!  d2f2= d2f / dx2  with second-order accuracy central difference
      d2f2(fm1,f0,fp1,xm1,x0,xp1)= &
       2*((fp1-f0)/(xp1-x0)-(f0-fm1)/(x0-xm1))/(xp1-xm1)

      bf(:,:,:,:)=xfold(:,:,:,6:8)
      bfdel(:,:,:,:)=x(:,:,:,6:8)-xfold(:,:,:,6:8)
      u(:,:,:)=xfold(:,:,:,2)
      udel(:,:,:)=x(:,:,:,2)-xfold(:,:,:,2)
      u0=u
      do jy=iy_first,iy_last
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
      klrb2(jx,jz,jy)=kap_ll(jx,jz)/(bf(jx,jz,jy,1)**2+bf(jx,jz,jy,2)**2+bf(jx,jz,jy,3)**2)
      enddo
      enddo
      enddo
!      call heat(u,www)

      do jy=iy_first+2,iy_last-2
      do jz=iz_first+2,iz_last-2
      do jx=jxam(jz-1),jxap(jz-1)

      bfx=bf(jx,jz,jy,1)
      bfy=bf(jx,jz,jy,2)
      bfz=bf(jx,jz,jy,3)
      bfx_dx=d1f2(bf(jx-1,jz,jy,1),bf(jx,jz,jy,1),bf(jx+1,jz,jy,1),xx(jx-1),xx(jx),xx(jx+1))
      bfx_dz=d1f2(bf(jx,jz-1,jy,1),bf(jx,jz,jy,1),bf(jx,jz+1,jy,1),zz(jz-1),zz(jz),zz(jz+1))
      bfx_dy=d1f2(bf(jx,jz,jy-1,1),bf(jx,jz,jy,1),bf(jx,jz,jy+1,1),yy(jy-1),yy(jy),yy(jy+1))
      bfz_dx=d1f2(bf(jx-1,jz,jy,3),bf(jx,jz,jy,3),bf(jx+1,jz,jy,3),xx(jx-1),xx(jx),xx(jx+1))
      bfz_dz=d1f2(bf(jx,jz-1,jy,3),bf(jx,jz,jy,3),bf(jx,jz+1,jy,3),zz(jz-1),zz(jz),zz(jz+1))
      bfz_dy=d1f2(bf(jx,jz,jy-1,3),bf(jx,jz,jy,3),bf(jx,jz,jy+1,3),yy(jy-1),yy(jy),yy(jy+1))
      bfy_dx=d1f2(bf(jx-1,jz,jy,2),bf(jx,jz,jy,2),bf(jx+1,jz,jy,2),xx(jx-1),xx(jx),xx(jx+1))
      bfy_dz=d1f2(bf(jx,jz-1,jy,2),bf(jx,jz,jy,2),bf(jx,jz+1,jy,2),zz(jz-1),zz(jz),zz(jz+1))
      bfy_dy=d1f2(bf(jx,jz,jy-1,2),bf(jx,jz,jy,2),bf(jx,jz,jy+1,2),yy(jy-1),yy(jy),yy(jy+1))
  
      klrb2_dx=d1f2(klrb2(jx-1,jz,jy),klrb2(jx,jz,jy),klrb2(jx+1,jz,jy),xx(jx-1),xx(jx),xx(jx+1))
      klrb2_dz=d1f2(klrb2(jx,jz-1,jy),klrb2(jx,jz,jy),klrb2(jx,jz+1,jy),zz(jz-1),zz(jz),zz(jz+1))
      klrb2_dy=d1f2(klrb2(jx,jz,jy-1),klrb2(jx,jz,jy),klrb2(jx,jz,jy+1),yy(jy-1),yy(jy),yy(jy+1))
      bdklrb2=bfx*klrb2_dx+bfz*klrb2_dz+bfy*klrb2_dy/xx(jx)
      kp=kap_pp(jx,jz)
      kp_dx=kap_pp_dx(jx,jz)

      udx = d1f2(u(jx-1,jz,jy),u(jx,jz,jy),u(jx+1,jz,jy),xxt(jx-1),xxt(jx),xxt(jx+1))
      udz = d1f2(u(jx,jz-1,jy),u(jx,jz,jy),u(jx,jz+1,jy),zzt(jz-1),zzt(jz),zzt(jz+1))
      udy = d1f2(u(jx,jz,jym(jy)),u(jx,jz,jy),u(jx,jz,jyp(jy)),yy(jy-1),yy(jy),yy(jy+1))
      ud2xx=d2f2(u(jx-1,jz,jy),u(jx,jz,jy),u(jx+1,jz,jy),xxt(jx-1),xxt(jx),xxt(jx+1))
      ud2zz=d2f2(u(jx,jz-1,jy),u(jx,jz,jy),u(jx,jz+1,jy),zzt(jz-1),zzt(jz),zzt(jz+1))
      ud2yy=d2f2(u(jx,jz,jym(jy)),u(jx,jz,jy),u(jx,jz,jyp(jy)),yyt(jy-1),yyt(jy),yyt(jy+1)) 
      ud2xz=(u(jx+1,jz+1,jy)+u(jx-1,jz-1,jy)-u(jx+1,jz-1,jy)-u(jx-1,jz+1,jy))/dxz/4
      ud2zy=(u(jx,jz+1,jy+1)+u(jx,jz-1,jy-1)-u(jx,jz-1,jy+1)-u(jx,jz+1,jy-1))/dzy/4
      ud2yx=(u(jx+1,jz,jy+1)+u(jx,jz-1,jy-1)-u(jx+1,jz,jy-1)-u(jx-1,jz,jy+1))/dxy/4

      www(jx,jz,jy)=-klrb2(jx,jz,jy)*( &
          bfx**2*ud2xx+bfz**2*ud2zz+bfy**2*ud2yy/xx(jx)**2 &
         +2*bfx*bfz*ud2xz+2*bfz*bfy*ud2zy/xx(jx)+2*bfy*bfx*ud2yx/xx(jx) &
         +(bfx*bfx_dx+bfz*bfx_dz+bfy*bfx_dy/xx(jx))*udx &
         +(bfx*bfz_dx+bfz*bfz_dz+bfy*bfz_dy/xx(jx))*udz &
         +(bfx*(bfy_dx-bfy/xx(jx))+bfz*bfy_dz+bfy*bfy_dy/xx(jx))/xx(jx)*udy) &
         -bdklrb2*(bfx*udx+bfz*udz+bfy*udy/xx(jx)) &
         -kap_pp(jx,jz)*(ud2xx+ud2zz+ud2yy/xx(jx)**2+udx/xx(jx)) &
         -kap_pp_dx(jx,jz)*udx-kap_pp_dz(jx,jz)*udz

      w1x =-klrb2(jx,jz,jy)*(bfx*bfx_dx+bfz*bfx_dz+bfy*bfx_dy/xx(jx))-bfx*bdklrb2-kp_dx-kp/xx(jx)
      w2xx=-klrb2(jx,jz,jy)*bfx*bfx-kp
      w2xz=-klrb2(jx,jz,jy)*bfx*bfz*2
      w2xy=-klrb2(jx,jz,jy)*bfx*bfy*2/xx(jx) 
 
      ccx0(jx)=3/dt-w2xx*2/dxx2
      ccxp(jx)=w2xx/dxx2+w2xz/dxz/4+w2xy/dxy/4+w1x/dxx/2
      ccxm(jx)=w2xx/dxx2-w2xz/dxz/4-w2xy/dxy/4-w1x/dxx/2
      ssx(jx) =(udel(jx,jz,jy)+u(jx,jz,jy))*3/dt-www(jx,jz,jy)+w2xx*ud2xx &
          -w2xz/dxz/4*(u0(jx-1,jz,jy)-u0(jx+1,jz,jy)-u0(jx-1,jz-1,jy)+u0(jx+1,jz-1,jy)+u(jx-1,jz-1,jy)-u(jx+1,jz-1,jy)) &
          -w2xy/dxy/4*(u0(jx-1,jz,jy)-u0(jx+1,jz,jy)-u0(jx-1,jz,jy-1)+u0(jx+1,jz,jy-1)+u(jx-1,jz,jy-1)-u(jx+1,jz,jy-1))
      enddo
      call tridag_real_period(ccxm,ccx0,ccxp,ssx,ufx,nx)
      do jx=jxam(jz-1),jxap(jz-1)
      u(jx,jz,jy)=ufx(jx)
      enddo
     
      enddo
      enddo

      bf(:,:,:,:)=bf(:,:,:,:)+bfdel(:,:,:,:)/3
      do jy=iy_first,iy_last
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
      klrb2(jx,jz,jy)=kap_ll(jx,jz)/(bf(jx,jz,jy,1)**2+bf(jx,jz,jy,2)**2+bf(jx,jz,jy,3)**2)
      enddo
      enddo
      enddo

 !     call heat(u,www)

      do jy=iy_first+2,iy_last-2
      do jx=ix_first+2,ix_last-2

      do jz=jzam(jx-1),jzap(jx-1)
      bfx=bf(jx,jz,jy,1)
      bfy=bf(jx,jz,jy,2)
      bfz=bf(jx,jz,jy,3)
      bfx_dx=d1f2(bf(jx-1,jz,jy,1),bf(jx,jz,jy,1),bf(jx+1,jz,jy,1),xx(jx-1),xx(jx),xx(jx+1))
      bfx_dz=d1f2(bf(jx,jz-1,jy,1),bf(jx,jz,jy,1),bf(jx,jz+1,jy,1),zz(jz-1),zz(jz),zz(jz+1))
      bfx_dy=d1f2(bf(jx,jz,jy-1,1),bf(jx,jz,jy,1),bf(jx,jz,jy+1,1),yy(jy-1),yy(jy),yy(jy+1))
      bfz_dx=d1f2(bf(jx-1,jz,jy,3),bf(jx,jz,jy,3),bf(jx+1,jz,jy,3),xx(jx-1),xx(jx),xx(jx+1))
      bfz_dz=d1f2(bf(jx,jz-1,jy,3),bf(jx,jz,jy,3),bf(jx,jz+1,jy,3),zz(jz-1),zz(jz),zz(jz+1))
      bfz_dy=d1f2(bf(jx,jz,jy-1,3),bf(jx,jz,jy,3),bf(jx,jz,jy+1,3),yy(jy-1),yy(jy),yy(jy+1))
      bfy_dx=d1f2(bf(jx-1,jz,jy,2),bf(jx,jz,jy,2),bf(jx+1,jz,jy,2),xx(jx-1),xx(jx),xx(jx+1))
      bfy_dz=d1f2(bf(jx,jz-1,jy,2),bf(jx,jz,jy,2),bf(jx,jz+1,jy,2),zz(jz-1),zz(jz),zz(jz+1))
      bfy_dy=d1f2(bf(jx,jz,jy-1,2),bf(jx,jz,jy,2),bf(jx,jz,jy+1,2),yy(jy-1),yy(jy),yy(jy+1))
  
      klrb2_dx=d1f2(klrb2(jx-1,jz,jy),klrb2(jx,jz,jy),klrb2(jx+1,jz,jy),xx(jx-1),xx(jx),xx(jx+1))
      klrb2_dz=d1f2(klrb2(jx,jz-1,jy),klrb2(jx,jz,jy),klrb2(jx,jz+1,jy),zz(jz-1),zz(jz),zz(jz+1))
      klrb2_dy=d1f2(klrb2(jx,jz,jy-1),klrb2(jx,jz,jy),klrb2(jx,jz,jy+1),yy(jy-1),yy(jy),yy(jy+1))
      bdklrb2=bfx*klrb2_dx+bfz*klrb2_dz+bfy*klrb2_dy/xx(jx)
      kp=kap_pp(jx,jz)
      kp_dx=kap_pp_dx(jx,jz)

      udx = d1f2(u(jx-1,jz,jy),u(jx,jz,jy),u(jx+1,jz,jy),xxt(jx-1),xxt(jx),xxt(jx+1))
      udz = d1f2(u(jx,jz-1,jy),u(jx,jz,jy),u(jx,jz+1,jy),zzt(jz-1),zzt(jz),zzt(jz+1))
      udy = d1f2(u(jx,jz,jym(jy)),u(jx,jz,jy),u(jx,jz,jyp(jy)),yy(jy-1),yy(jy),yy(jy+1))
      ud2xx=d2f2(u(jx-1,jz,jy),u(jx,jz,jy),u(jx+1,jz,jy),xxt(jx-1),xxt(jx),xxt(jx+1))
      ud2zz=d2f2(u(jx,jz-1,jy),u(jx,jz,jy),u(jx,jz+1,jy),zzt(jz-1),zzt(jz),zzt(jz+1))
      ud2yy=d2f2(u(jx,jz,jym(jy)),u(jx,jz,jy),u(jx,jz,jyp(jy)),yyt(jy-1),yyt(jy),yyt(jy+1)) 
      ud2xz=(u(jx+1,jz+1,jy)+u(jx-1,jz-1,jy)-u(jx+1,jz-1,jy)-u(jx-1,jz+1,jy))/dxz/4
      ud2zy=(u(jx,jz+1,jy+1)+u(jx,jz-1,jy-1)-u(jx,jz-1,jy+1)-u(jx,jz+1,jy-1))/dzy/4
      ud2yx=(u(jx+1,jz,jy+1)+u(jx,jz-1,jy-1)-u(jx+1,jz,jy-1)-u(jx-1,jz,jy+1))/dxy/4

      www(jx,jz,jy)=-klrb2(jx,jz,jy)*( &
          bfx**2*ud2xx+bfz**2*ud2zz+bfy**2*ud2yy/xx(jx)**2 &
         +2*bfx*bfz*ud2xz+2*bfz*bfy*ud2zy/xx(jx)+2*bfy*bfx*ud2yx/xx(jx) &
         +(bfx*bfx_dx+bfz*bfx_dz+bfy*bfx_dy/xx(jx))*udx &
         +(bfx*bfz_dx+bfz*bfz_dz+bfy*bfz_dy/xx(jx))*udz &
         +(bfx*(bfy_dx-bfy/xx(jx))+bfz*bfy_dz+bfy*bfy_dy/xx(jx))/xx(jx)*udy) &
         -bdklrb2*(bfx*udx+bfz*udz+bfy*udy/xx(jx)) &
         -kap_pp(jx,jz)*(ud2xx+ud2zz+ud2yy/xx(jx)**2+udx/xx(jx)) &
         -kap_pp_dx(jx,jz)*udx-kap_pp_dz(jx,jz)*udz
 
      w1z =-klrb2(jx,jz,jy)*(bfx*bfz_dx+bfz*bfz_dz+bfy*bfz_dy/xx(jx))-bfz*bdklrb2-kp_dz
      w2zz=-klrb2(jx,jz,jy)*bfz*bfz-kp
      w2zx=-klrb2(jx,jz,jy)*bfz*bfx*2
      w2zy=-klrb2(jx,jz,jy)*bfz*bfy*2/xx(jx)

      ccz0(jz)=3/dt-w2zz*2/dzz2
      cczp(jz)=w2zz/dzz2+w2zx/dxz/4+w2zy/dzy/4+w1z/dzz/2
      cczm(jz)=w2zz/dzz2-w2zx/dxz/4-w2zy/dzy/4-w1z/dzz/2
      ssz(jz) =(udel(jx,jz,jy)+u(jx,jz,jy))*3/dt-www(jx,jz,jy)+w2zz*ud2zz &
          -w2zx/dxz/4*(u0(jx,jz-1,jy)-u0(jx,jz+1,jy)-u0(jx-1,jz-1,jy)+u0(jx-1,jz+1,jy)+u(jx-1,jz-1,jy)-u(jx-1,jz+1,jy)) &
          -w2zy/dzy/4*(u0(jx,jz-1,jy)-u0(jx,jz+1,jy)-u0(jx,jz-1,jy-1)+u0(jx,jz+1,jy-1)+u(jx,jz-1,jy-1)-u(jx,jz+1,jy-1))
      enddo
      call tridag_real_period(cczm,ccz0,cczp,ssz,ufz,nz)
      do jz=jzam(jx-1),jzap(jx-1)
      u(jx,jz,jy)=ufz(jz)
      enddo
       
      enddo
      enddo
      
      bf(:,:,:,:)=bf(:,:,:,:)+bfdel(:,:,:,:)/3
      do jy=iy_first,iy_last
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
      klrb2(jx,jz,jy)=kap_ll(jx,jz)/(bf(jx,jz,jy,1)**2+bf(jx,jz,jy,2)**2+bf(jx,jz,jy,3)**2)
      enddo
      enddo
      enddo
 !     call heat(u,www)
      
      do jz=iz_first+2,iz_last-2
      do jx=ix_first+2,ix_last-2
      if(pst(jx,jz).lt.psia1) then
      do jy=iy_first,iy_last
      bfx=bf(jx,jz,jy,1)
      bfy=bf(jx,jz,jy,2)
      bfz=bf(jx,jz,jy,3)
      bfx_dx=d1f2(bf(jx-1,jz,jy,1),bf(jx,jz,jy,1),bf(jx+1,jz,jy,1),xx(jx-1),xx(jx),xx(jx+1))
      bfx_dz=d1f2(bf(jx,jz-1,jy,1),bf(jx,jz,jy,1),bf(jx,jz+1,jy,1),zz(jz-1),zz(jz),zz(jz+1))
      bfx_dy=d1f2(bf(jx,jz,jy-1,1),bf(jx,jz,jy,1),bf(jx,jz,jy+1,1),yy(jy-1),yy(jy),yy(jy+1))
      bfz_dx=d1f2(bf(jx-1,jz,jy,3),bf(jx,jz,jy,3),bf(jx+1,jz,jy,3),xx(jx-1),xx(jx),xx(jx+1))
      bfz_dz=d1f2(bf(jx,jz-1,jy,3),bf(jx,jz,jy,3),bf(jx,jz+1,jy,3),zz(jz-1),zz(jz),zz(jz+1))
      bfz_dy=d1f2(bf(jx,jz,jy-1,3),bf(jx,jz,jy,3),bf(jx,jz,jy+1,3),yy(jy-1),yy(jy),yy(jy+1))
      bfy_dx=d1f2(bf(jx-1,jz,jy,2),bf(jx,jz,jy,2),bf(jx+1,jz,jy,2),xx(jx-1),xx(jx),xx(jx+1))
      bfy_dz=d1f2(bf(jx,jz-1,jy,2),bf(jx,jz,jy,2),bf(jx,jz+1,jy,2),zz(jz-1),zz(jz),zz(jz+1))
      bfy_dy=d1f2(bf(jx,jz,jy-1,2),bf(jx,jz,jy,2),bf(jx,jz,jy+1,2),yy(jy-1),yy(jy),yy(jy+1))
  
      klrb2_dx=d1f2(klrb2(jx-1,jz,jy),klrb2(jx,jz,jy),klrb2(jx+1,jz,jy),xx(jx-1),xx(jx),xx(jx+1))
      klrb2_dz=d1f2(klrb2(jx,jz-1,jy),klrb2(jx,jz,jy),klrb2(jx,jz+1,jy),zz(jz-1),zz(jz),zz(jz+1))
      klrb2_dy=d1f2(klrb2(jx,jz,jy-1),klrb2(jx,jz,jy),klrb2(jx,jz,jy+1),yy(jy-1),yy(jy),yy(jy+1))
      bdklrb2=bfx*klrb2_dx+bfz*klrb2_dz+bfy*klrb2_dy/xx(jx)
      kp=kap_pp(jx,jz)
      kp_dx=kap_pp_dx(jx,jz)

      udx = d1f2(u(jx-1,jz,jy),u(jx,jz,jy),u(jx+1,jz,jy),xxt(jx-1),xxt(jx),xxt(jx+1))
      udz = d1f2(u(jx,jz-1,jy),u(jx,jz,jy),u(jx,jz+1,jy),zzt(jz-1),zzt(jz),zzt(jz+1))
      udy = d1f2(u(jx,jz,jym(jy)),u(jx,jz,jy),u(jx,jz,jyp(jy)),yy(jy-1),yy(jy),yy(jy+1))
      ud2xx=d2f2(u(jx-1,jz,jy),u(jx,jz,jy),u(jx+1,jz,jy),xxt(jx-1),xxt(jx),xxt(jx+1))
      ud2zz=d2f2(u(jx,jz-1,jy),u(jx,jz,jy),u(jx,jz+1,jy),zzt(jz-1),zzt(jz),zzt(jz+1))
      ud2yy=d2f2(u(jx,jz,jym(jy)),u(jx,jz,jy),u(jx,jz,jyp(jy)),yyt(jy-1),yyt(jy),yyt(jy+1)) 
      ud2xz=(u(jx+1,jz+1,jy)+u(jx-1,jz-1,jy)-u(jx+1,jz-1,jy)-u(jx-1,jz+1,jy))/dxz/4
      ud2zy=(u(jx,jz+1,jy+1)+u(jx,jz-1,jy-1)-u(jx,jz-1,jy+1)-u(jx,jz+1,jy-1))/dzy/4
      ud2yx=(u(jx+1,jz,jy+1)+u(jx,jz-1,jy-1)-u(jx+1,jz,jy-1)-u(jx-1,jz,jy+1))/dxy/4

      www(jx,jz,jy)=-klrb2(jx,jz,jy)*( &
          bfx**2*ud2xx+bfz**2*ud2zz+bfy**2*ud2yy/xx(jx)**2 &
         +2*bfx*bfz*ud2xz+2*bfz*bfy*ud2zy/xx(jx)+2*bfy*bfx*ud2yx/xx(jx) &
         +(bfx*bfx_dx+bfz*bfx_dz+bfy*bfx_dy/xx(jx))*udx &
         +(bfx*bfz_dx+bfz*bfz_dz+bfy*bfz_dy/xx(jx))*udz &
         +(bfx*(bfy_dx-bfy/xx(jx))+bfz*bfy_dz+bfy*bfy_dy/xx(jx))/xx(jx)*udy) &
         -bdklrb2*(bfx*udx+bfz*udz+bfy*udy/xx(jx)) &
         -kap_pp(jx,jz)*(ud2xx+ud2zz+ud2yy/xx(jx)**2+udx/xx(jx)) &
         -kap_pp_dx(jx,jz)*udx-kap_pp_dz(jx,jz)*udz

      w1y =(-klrb2(jx,jz,jy)*(bfx*(bfy_dx-bfy/xx(jx))+bfz*bfy_dz+bfy*bfy_dy/xx(jx))-bfy*bdklrb2-kp_dy/xx(jx))/xx(jx)
      w2yy=(-klrb2(jx,jz,jy)*bfy*bfy-kp)/xx(jx)**2
      w2yx=-klrb2(jx,jz,jy)*bfy*bfx*2/xx(jx)
      w2yz=-klrb2(jx,jz,jy)*bfy*bfz*2/xx(jx)
      
      ccy0(jy)=3/dt-w2yy*2/dyy2
      ccyp(jy)=w2yy/dyy2+w2yx/dxy/4+w2yz/dzy/4+w1y/dyy/2
      ccym(jy)=w2yy/dyy2-w2yx/dxy/4-w2yz/dzy/4-w1y/dyy/2
      ssy(jy) =(udel(jx,jz,jy)+u(jx,jz,jy))*3/dt-www(jx,jz,jy)+w2yy*ud2yy &
          -w2yx/dxy/4*(u0(jx,jz,jy-1)-u0(jx,jz,jy+1)-u0(jx-1,jz,jy-1)+u0(jx-1,jz,jy+1)+u(jx-1,jz,jy-1)-u(jx-1,jz,jy+1)) &
          -w2yz/dzy/4*(u0(jx,jz,jy-1)-u0(jx,jz,jy+1)-u0(jx,jz-1,jy-1)+u0(jx,jz-1,jy+1)+u(jx,jz-1,jy-1)-u(jx,jz-1,jy+1))


      enddo
      call tridag_real_period(ccym,ccy0,ccyp,ssy,ufy,ny)
      do jy=1,my
      u(jx,jz,jy)=ufy(jy)
      enddo
      endif

      enddo
      enddo

      return
      end
!ws*****************************************************
      !subroutine pll_PETSC
      !USE DECLARE
      !call convtb
      !call set_A_S(dt)
      !call possion_solver_3D(x(:,:,:,2))
      !call mpi_transfersm(x(:,:,:,2),1)
      !return
      !end
!ws*****************************************************
      !subroutine pll_PETSC_T1
      !USE DECLARE
      !call convtb
      !call set_A_S(dt)
      !call possion_solver_3D_T1(x(:,:,:,2))
      !call mpi_transfersm(x(:,:,:,2),1)
      !return
      !end
!ws************************************************************************
      subroutine pll_smthpline(npll)
      USE DECLARE
      integer npll

      select case(lscheme)
      case(1)
      call smthp_traceline(npll)
      case(2)
      call smthp_traceline_v1(npll)
      case(3)
      call smthp_traceline_v2(npll)
      case(4) 
      call smthp2_traceline(npll)
      case(5)
      call smthp_traceline_5p(npll)
      case(6)
      call smthp2_traceline_tm(npll)
      case(7)
      call smthp_traceline_spec(npll)
!      case default
      end select
      return
      end

!ws************************************************************************
      subroutine pll_soundwave(npll)
      USE DECLARE
      integer npll

      select case(lscheme)
      case(1)
      call artif_sound_replace(npll)
      case(2)
      call artif_sound_replace_RK(npll)
      case(3)
      call artif_sound_replace_Lax(npll)
      case(4)
      call artif_sound_implicity(npll)
!      case(5)
!      call artif_sound_pt(npll)
!      case default
      end select
      return
      end
!ws*****************************************************

      subroutine convtb
      USE DECLARE
      real*8, dimension(my) :: wwy
      include 'mpif.h'
!
!  define statement functions
!  d1f2= d f / dx  with second-order accuracy central difference
      d1f2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(fp1-f0) &
         -(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)
!  d1f2m= d f / dx  with second-order one-sided  difference involving -2 -1 and 0
!  points
      d1f2m(fm2,fm1,f0,xm2,xm1,x0)= &
        ((xm2-x0)/(xm1-x0)*(fm1-f0) &
         -(xm1-x0)/(xm2-x0)*(fm2-f0))/(xm2-xm1)
!  d1f2p= d f / dx  with second-order one-sided  difference involving 2 1 and 0
!  points
      d1f2p(fp2,fp1,f0,xp2,xp1,x0)= &
        ((xp2-x0)/(xp1-x0)*(fp1-f0) &
         -(xp1-x0)/(xp2-x0)*(fp2-f0))/(xp2-xp1)
!  d1fc= d f / dx  with fourth-order accuracy central difference
      d1fc(fm2,fm1,f0,fp1,fp2,a,b,c,d)= &
       a*(fp1-f0)+b*(f0-fm1)+c*(fp2-f0)+d*(f0-fm2)
!  d2fc= d2 f / dx2   with third-order accuracy central difference
      d2fc(fm2,fm1,f0,fp1,fp2,a,b,c,d)= &
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
      d1fbp(fm1,f0,fp1,fp2,a,b,c)= &
       a*(fm1-f0)+b*(fp1-f0)+c*(fp2-f0)
!  d1fbm= d f / dx  with  one-sided-bias  difference involving -2 -1  0 and 1
!  points
      d1fbm(fm2,fm1,f0,fp1,a,b,c)= &
       a*(f0-fp1)+b*(f0-fm1)+c*(f0-fm2)
!  d2f2= d2f / dx2  with second-order accuracy central difference
      d2f2(fm1,f0,fp1,xm1,x0,xp1)= &
       2*((fp1-f0)/(xp1-x0)-(f0-fm1)/(x0-xm1))/(xp1-xm1)
!  d2fbp= d2 f / dx2  with  one-sided-bias  difference involving -1 0  1 and 2
!  points
      d2fbp(fm1,f0,fp1,fp2,a,b,c)= &
       2*(a*(fm1-f0)+b*(fp1-f0)+c*(fp2-f0))
!  d2fbm= d2 f / dx2  with  one-sided-bias  difference involving -2 -1  0 and 1
!  points
      d2fbm(fm2,fm1,f0,fp1,a,b,c)= &
       2*(a*(fp1-f0)+b*(fm1-f0)+c*(fm2-f0))

       integer status(mpi_status_size)
!      init1=1
!      init2=1
!      jxs = ix_first
!      jxe = ix_last
!      jzs = iz_first
!      jze = iz_last
!      if (nrank == 0)        jxs=ix_first
!      if (nrank == nsize - 1) jxe =ix_last

 !      call mpi_transfer8(x1)
      if(spectral) then
      do 18 jx=ix_first,ix_last
      do 18 jz=iz_first,iz_last
      do 17 m=6,8
      wwy(:)=x1(jx,jz,:,m)
      call mpi_transfersy1(wwy,data0)

!      do 2 jy=1,my
!        data0(jy)=x1(jx,jz,jy,m)
!    2 continue
!!    
!!   76 CALL DRCFT2(init1,data,my+2,spec,my/2+1,my,mz,-1,scale,AUX11,naux1, &
!!                AUX12,naux2)
!!      if (init1.eq.1) then
!!       init1 = init1-1
!!       goto 76
!!      endif
!!
      call dfftw_plan_dft_r2c_1d(plan,myt,data0,spec,FFTW_ESTIMATE)
	  call dfftw_execute_dft_r2c(plan,data0,spec)
	  call dfftw_destroy_plan(plan)

      do 3 jy=1,myt/2
      spec1(jy)=spec(jy)*fac
    3 continue
      
      if(filt) then
      do 5 jy=1,myt/2
      if(m.eq.2) then
      spec(jy)=spec1(jy) 
      else
      spec(jy)=spec1(jy)*fmode(jy)
      endif
    5 continue  
      spec(myt/2+1)=c0

!    77 CALL DCRFT2(init2,spec,my/2+1,data,my+2,my,mz,1,scale,AUX21,naux1,  &
!               AUX22,naux2)
!      if (init2.eq.1) then
!       init2 = init2-1
!       goto 77
!      endif

      call dfftw_plan_dft_c2r_1d(plan,myt,spec,data0,FFTW_ESTIMATE)
	  call dfftw_execute_dft_c2r(plan,spec,data0)
	  call dfftw_destroy_plan(plan)

      do 7 jy=iy_first,iy_last
      x1(jx,jz,jy,m)=data0(nrky(nrank)*mym+jy-2)
    7 continue 
      endif

      do ky=1,myt/2
      xkc(jx,jz,ky,m)=real(spec1(ky))
      xks(jx,jz,ky,m)=-imag(spec1(ky))
      enddo


      do 10 i=1,myt/2
      if(m.eq.2) then
      spec(i)=c1*(i-1)*spec1(i)
      else
      spec(i)=c1*(i-1)*spec1(i)*fmode(i)
      endif
   10 continue
      spec(myt/2+1)=c0
!      
!   78 CALL DCRFT2(init2,spec,my/2+1,data,my+2,my,mz,1,scale,AUX21,naux1,  &
!               AUX22,naux2)
!      if (init2.eq.1) then
!       init2 = init2-1
!       goto 78
!      endif
!
      call dfftw_plan_dft_c2r_1d(plan,myt,spec,data0,FFTW_ESTIMATE)
	  call dfftw_execute_dft_c2r(plan,spec,data0)
	  call dfftw_destroy_plan(plan)

      do 12 jy=iy_first,iy_last
      xy(jx,jz,jy,m)=data0(nrky(nrank)*mym+jy-2)
   12 continue
   
   do 13 i=1,myt/2
      spec(i)=-(i-1)**2*spec1(i)
   13 continue
      spec(myt/2+1)=c0
!      
!   78 CALL DCRFT2(init2,spec,my/2+1,data,my+2,my,mz,1,scale,AUX21,naux1,  &
!               AUX22,naux2)
!      if (init2.eq.1) then
!       init2 = init2-1
!       goto 78
!      endif
!
      call dfftw_plan_dft_c2r_1d(plan,myt,spec,data0,FFTW_ESTIMATE)
	  call dfftw_execute_dft_c2r(plan,spec,data0)
	  call dfftw_destroy_plan(plan)

      do 14 jy=iy_first,iy_last
      xy2(jx,jz,jy,m)=data0(nrky(nrank)*mym+jy-2)
   14 continue     
!      
   17 continue
     
   18 continue

      else !!not spectral
      do 19 m=6,8
      do 15 jz=iz_first,iz_last
      do 15 jx=ix_first,ix_last
!      do jy=2,my-1
!      xy(jx,jz,jy,m)=d1f2(x1(jx,jz,jy-1,m),x1(jx,jz,jy,m),x1(jx,jz,jy+1,m),yy(jy-1),yy(jy),yy(jy+1))      
!      xy2(jx,jz,jy,m)=d2f2(x1(jx,jz,jy-1,m),x1(jx,jz,jy,m),x1(jx,jz,jy+1,m),yy(jy-1),yy(jy),yy(jy+1))
!      enddo
!      xy(jx,jz,1,m)=d1f2(x1(jx,jz,my,m),x1(jx,jz,1,m),x1(jx,jz,2,m),yy(0),yy(1),yy(2))
!      xy(jx,jz,my,m)=d1f2(x1(jx,jz,my-1,m),x1(jx,jz,my,m),x1(jx,jz,1,m),yy(my-1),yy(my),yy(my+1))
!      xy2(jx,jz,1,m)=d2f2(x1(jx,jz,my,m),x1(jx,jz,1,m),x1(jx,jz,2,m),yy(0),yy(1),yy(2))
!      xy2(jx,jz,my,m)=d2f2(x1(jx,jz,my-1,m),x1(jx,jz,my,m),x1(jx,jz,1,m),yy(my-1),yy(my),yy(my+1))
      do jy=iy_first+2,iy_last-2
!      xy(jx,jz,jy,m) =d1fc(x1(jx,jz,jy-2,m),x1(jx,jz,jy-1,m),x1(jx,jz,jy,m),x1(jx,jz,jy+1,m),x1(jx,jz,jy+2,m),ay1(jy),by1(jy),cy1(jy),dy1(jy))      
!      xy2(jx,jz,jy,m)=d2fc(x1(jx,jz,jy-2,m),x1(jx,jz,jy-1,m),x1(jx,jz,jy,m),x1(jx,jz,jy+1,m),x1(jx,jz,jy+2,m),ay2(jy),by2(jy),cy2(jy),dy2(jy))
      xy(jx,jz,jy,m) =d1f2(x1(jx,jz,jy-1,m),x1(jx,jz,jy,m),x1(jx,jz,jy+1,m),yy(jy-1),yy(jy),yy(jy+1))
      xy2(jx,jz,jy,m)=d2f2(x1(jx,jz,jy-1,m),x1(jx,jz,jy,m),x1(jx,jz,jy+1,m),yy(jy-1),yy(jy),yy(jy+1))
      enddo

   15 continue
   19 continue 
      endif !spectral

      do 30 m=6,8
      do 30 jy=iy_first,iy_last
      do 21 jz=iz_first+2,iz_last-2
      do 21 jx=ix_first+2,ix_last-2
      x1r(jx,jz,jy,m)=d1f2(x1(jx-1,jz,jy,m),x1(jx,jz,jy,m) &
         ,x1(jx+1,jz,jy,m),xx(jx-1),xx(jx),xx(jx+1))

      xr2(jx,jz,jy,m)=d2f2(x1(jx-1,jz,jy,m),x1(jx,jz,jy,m) &
         ,x1(jx+1,jz,jy,m),xx(jx-1),xx(jx),xx(jx+1))
      x1z(jx,jz,jy,m)=d1f2(x1(jx,jz-1,jy,m),x1(jx,jz,jy,m) &
         ,x1(jx,jz+1,jy,m),zz(jz-1),zz(jz),zz(jz+1))
      xz2(jx,jz,jy,m)=d2f2(x1(jx,jz-1,jy,m),x1(jx,jz,jy,m) &
         ,x1(jx,jz+1,jy,m),zz(jz-1),zz(jz),zz(jz+1))

!     x1r(jx,jz,jy,m)=d1fc(x1(jx-2,jz,jy,m),x1(jx-1,jz,jy,m),x1(jx,jz,jy,m) &
!         ,x1(jx+1,jz,jy,m),x1(jx+2,jz,jy,m),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
!     xr2(jx,jz,jy,m)=d2fc(x1(jx-2,jz,jy,m),x1(jx-1,jz,jy,m),x1(jx,jz,jy,m) &
!         ,x1(jx+1,jz,jy,m),x1(jx+2,jz,jy,m),ax2(jx),bx2(jx),cx2(jx),dx2(jx))
!     x1z(jx,jz,jy,m)=d1fc(x1(jx,jz-2,jy,m),x1(jx,jz-1,jy,m),x1(jx,jz,jy,m) &
!         ,x1(jx,jz+1,jy,m),x1(jx,jz+2,jy,m),az1(jz),bz1(jz),cz1(jz),dz1(jz))
!     xz2(jx,jz,jy,m)=d1fc(x1(jx,jz-2,jy,m),x1(jx,jz-1,jy,m),x1(jx,jz,jy,m) &
!         ,x1(jx,jz+1,jy,m),x1(jx,jz+2,jy,m),az2(jz),bz2(jz),cz2(jz),dz2(jz))

      xr(jx,jz,jy,m)=xint_dx(jx,jz,m)+x1r(jx,jz,jy,m)  
      xz(jx,jz,jy,m)=xint_dz(jx,jz,m)+x1z(jx,jz,jy,m)


   21 continue
   30 continue
   
      do jy=iy_first,iy_last
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
      bf2(jx,jz,jy)=x(jx,jz,jy,6)**2+x(jx,jz,jy,7)**2+x(jx,jz,jy,8)**2
      bf2dx(jx,jz,jy)=2*(x(jx,jz,jy,6)*xr(jx,jz,jy,6)+x(jx,jz,jy,7)*xr(jx,jz,jy,7)+x(jx,jz,jy,8)*xr(jx,jz,jy,8))
      bf2dz(jx,jz,jy)=2*(x(jx,jz,jy,6)*xz(jx,jz,jy,6)+x(jx,jz,jy,7)*xz(jx,jz,jy,7)+x(jx,jz,jy,8)*xz(jx,jz,jy,8))
      bf2dy(jx,jz,jy)=2*(x(jx,jz,jy,6)*xy(jx,jz,jy,6)+x(jx,jz,jy,7)*xy(jx,jz,jy,7)+x(jx,jz,jy,8)*xy(jx,jz,jy,8))
      enddo
      enddo
      enddo 


      return
      end

!ws**********************************************************************
      subroutine find_Opoint_Z0
      USE DECLARE
      USE DECLARE_OXpoint
      integer jxtom,jxtop,jxom(1),jxop(1),jxo(1)
      real*8  bzmaxm,bzmaxp,bzmaxo
      real*8, dimension(mx) :: bo
      include 'mpif.h'  
       
      bo(:)=Ef(:,jzmode,3,2)
  

      if(nrank==nrank_mode-1) then
      bzmaxm=maxval(bo)
      jxom  =maxloc(bo)
      jxtom =jxom(1)+(nrkx_mode-1)*mxm-2

      CALL MPI_Send( bzmaxm, 1, MPI_DOUBLE_PRECISION, nrank_mode, 1,  &
		      MPI_COMM_WORLD,ierror )
      CALL MPI_Send( jxtom,   1, MPI_INTEGER, nrank_mode, 2,  &
		      MPI_COMM_WORLD,ierror )
      endif

      if(nrank==nrank_mode+1) then
      bzmaxp=maxval(bo)
      jxop  =maxloc(bo)
      jxtop =jxop(1)+(nrkx_mode+1)*mxm-2

      CALL MPI_Send( bzmaxp, 1, MPI_DOUBLE_PRECISION, nrank_mode, 3,  &
		      MPI_COMM_WORLD,ierror )
      CALL MPI_Send( jxtop,   1, MPI_INTEGER, nrank_mode, 4,  &
		      MPI_COMM_WORLD,ierror )
      endif

      if(nrank==nrank_mode) then
      bzmaxo=maxval(bo)
      jxo   =maxloc(bo)
      jxto  =jxo(1)+(nrkx_mode)*mxm-2

      CALL MPI_Recv( bzmaxm, 1, MPI_DOUBLE_PRECISION, nrank_mode-1, 1,  &
		      MPI_COMM_WORLD, status,ierror  )
      CALL MPI_Recv( jxtom,  1, MPI_INTEGER, nrank_mode-1, 2,  &
		      MPI_COMM_WORLD, status,ierror  )
      CALL MPI_Recv( bzmaxp, 1, MPI_DOUBLE_PRECISION, nrank_mode+1, 3,  &
		      MPI_COMM_WORLD, status,ierror  )
      CALL MPI_Recv( jxtop,  1, MPI_INTEGER, nrank_mode+1, 4,  &
		      MPI_COMM_WORLD, status,ierror  )
      
      if(bzmaxo.lt.bzmaxm) then
      bzmaxo=bzmaxm
      jxto=jxtom
      endif
      if(bzmaxo.lt.bzmaxp) then
      bzmaxo=bzmaxp
      jxto=jxtop
      endif

      endif
      
      call MPI_Bcast(jxto, 1, MPI_INTEGER, nrank_mode, MPI_COMM_WORLD,ierror )

      jzto=jztmode

      xx_O=xxt(jxto+1)
      zz_O=zzt(jzto)
      ps_O=pst(jxto+1,jzto)
      rr_O=sqrt((ps_O-psmin)/(psmax-psmin))
      ! th_O=tcht(jxto,jzto)
      if(nrank.eq.0) then
      write(311,700) time,jxto,jzto,xx_O,zz_O,ps_O,rr_O
700   format((1x,e12.5,2i,4(1x,e12.5)))
      endif

      return
      end

!ws*****************************************************
      subroutine distribution_cd_OXpoint(ww)
      USE DECLARE
      USE DECLARE_OXpoint
      real*8, dimension(mx,mz) :: ww
      real*8 fnV,fnV1,pscd,phase,cdc,br_old

      include 'mpif.h'      
      
      nrkyo=0    
      jyo=3
      call find_OXpoint(ww)
      if(br_max .ge. br_old) br_rise=.true.
      
      if(br_rise .or. br_max .lt. br_lim) then
      cdc=0. 
      else
      cdc=(br_max-br_lim)/(br_max0-br_lim)
      cdc=1.0*tanh(pi*cdc)/tanh(pi)
!      cdc=cdc**0.5
      endif
      
      br_old=br_max
      
      do jy=iy_first,iy_last
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
      phase=nmode*(yy(jy)-yy_X)+mmode*(tcxz(jx,jz)-tc_X)
!      pscd=ps_O+(ps_X-ps_O)*(1-dcos(nmode*(yy(jy)-yy_O)+mmode*(tcxz(jx,jz)-tc_O)))
!      fn_cdy(jx,jz,jy)=exp(-(psi(jx,jz)-pscd)**2/delcd**2)*(1.0+dcos(nmode*(yy(jy)-yy_O)+mmode*(tcxz(jx,jz)-tc_O)))
      pscd=ps_X+(ps_O-ps_X)*(1-dcos(phase))/2+psshift
      fn_cdy(jx,jz,jy)=exp(-(psi(jx,jz)-pscd)**2/delcd**2)*(1.0-cdc*dcos(phase))

      enddo
      enddo
      enddo

      fnV1=0
!      do jy=iy_first+2,iy_last-2  
      do jz=iz_first+2,iz_last-2
      do jx=ix_first+2,ix_last-2 
      if(psi(jx,jz) .lt. psia) then 
!      fnV1=fnV1+fn_cdy(jx,jz)*(xx(jx+1)-xx(jx-1))*(zz(jz+1)-zz(jz-1))/4
        bb2=x(jx,jz,3,6)**2+x(jx,jz,3,8)**2+x(jx,jz,3,7)**2
        bb=sqrt(bb2)
      fnV1=fnV1+fn_cdy(jx,jz,3)*x(jx,jz,3,7)/bb*(xx(jx+1)-xx(jx-1))*(zz(jz+1)-zz(jz-1))/4
      endif
      enddo
      enddo
!      enddo

      CALL MPI_ALLREDUCE(fnV1,fnV,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                        MPI_COMM_WORLD,IERROR)

      fn_cdy(:,:,:)=fn_cdy(:,:,:)/fnV

      return
      end
!ws**********************************************************************
      subroutine find_maxmin_1(ww,nky,wmax,wmin,jxlmax,jzlmax,jxlmin,jzlmin)
      use DECLARE_grid
      integer nky,jxlmax,jzlmax,jxlmin,jzlmin
      real*8  wmax,wmin
      real*8, dimension(mx,mz) :: ww
      real*8, dimension(6,0:npr-1) :: wm
      include 'mpif.h'  

      call funmax(ww,wm(1,nrank),wm(4,nrank),wm(2,nrank),wm(5,nrank),wm(3,nrank),wm(6,nrank),xx,zz,mx,mz)

!      write(*,600) nrank,wm(:,nrank)
!
      call MPI_Allgather(wm(:,nrank),6,MPI_DOUBLE_PRECISION,wm,6,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERROR)
!
!      if(nrank==11) then
!      write(*,*) '*********************after allgather*********************'
!      do ik=0,npr-1
!      write(*,600) ik,(wm(ii,ik),ii=1,6)
!      enddo
!600   format(i3,6(1x,e12.5)) 
!      endif
      wmax=-10000.
      wmin=10000.

      do 10 nrk=nky*nprxz,(nky+1)*nprxz-1
      if(wm(1,nrk).gt.wmax) then
      wmax=wm(1,nrk)
      xlmax=wm(2,nrk)
      zlmax=wm(3,nrk)
      endif
      if(wm(4,nrk).lt.wmin) then
      wmin=wm(4,nrk)
      xlmin=wm(5,nrk)
      zlmin=wm(6,nrk)
      endif
   10 continue

      jxlmax=int((xlmax-xxt(1))/dxx+0.01)+1
      jzlmax=int((zlmax-zzt(1))/dzz+0.01)+1

      jxlmin=int((xlmin-xxt(1))/dxx+0.01)+1
      jzlmin=int((zlmin-zzt(1))/dzz+0.01)+1

      return
      end

!ws**********************************************************************
      subroutine find_maxmin(ww,nky,wmax,wmin,jxlmax,jzlmax,jxlmin,jzlmin)
      use DECLARE_grid
      integer nky,jxlmax,jzlmax,jxlmin,jzlmin
      real*8  wmax,wmin,xlmax,zlmax,xlmin,zlmin
      real*8, dimension(mx,mz) :: ww
      real*8, dimension(6,0:npr-1) :: wm
      include 'mpif.h'  

      wm(1,nrank)=-10000.
      wm(4,nrank)=10000.
      do 10 jz=3,mz-2
      do 10 jx=3,mx-2
      if(ww(jx,jz).gt.wm(1,nrank)) then
      wm(1,nrank)=ww(jx,jz)
      jxlmax=jx
      jzlmax=jz
      endif
      if(ww(jx,jz).lt.wm(4,nrank)) then
      wm(4,nrank)=ww(jx,jz)
      jxlmin=jx
      jzlmin=jz
      endif
   10 continue
      wm(2,nrank)=real(nrankx*mxm+jxlmax-2)
      wm(3,nrank)=real(nrankz*mzm+jzlmax-2)

      wm(5,nrank)=real(nrankx*mxm+jxlmin-2)
      wm(6,nrank)=real(nrankz*mzm+jzlmin-2)

!      write(*,600) nrank,wm(:,nrank)
!
      call MPI_Allgather(wm(:,nrank),6,MPI_DOUBLE_PRECISION,wm,6,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERROR)
!
!      if(nrank==11) then
!      write(*,*) '*********************after allgather*********************'
!      open(unit=211,file='maxmin.dat',status='unknown',form='formatted')
!      do ik=0,npr-1
!      write(*,600) ik,(wm(ii,ik),ii=1,6)
!      write(211,600) ik,(wm(ii,ik),ii=1,6)
!      enddo
! 
!      close(211)
!      endif

      wmax=-10000.
      wmin=10000.

      do 11 nrk=nky*nprxz,(nky+1)*nprxz-1
      if(wm(1,nrk).gt.wmax) then
      wmax=wm(1,nrk)
      xlmax=wm(2,nrk)
      zlmax=wm(3,nrk)
      endif
      if(wm(4,nrk).lt.wmin) then
      wmin=wm(4,nrk)
      xlmin=wm(5,nrk)
      zlmin=wm(6,nrk)
      endif
   11 continue

      jxlmax=int(xlmax+0.01)
      jzlmax=int(zlmax+0.01)

      jxlmin=int(xlmin+0.01)
      jzlmin=int(zlmin+0.01)

!      write(*,*) '*********************after maxmin*********************'
      if(nrank==0) then
      write(211,700) time,wmax,jxlmax,jzlmax,wmin,jxlmin,jzlmin
      endif

!600   format(i3,6(1x,e12.5))
700   format(e12.5,2(1x,e12.5,2(1x,i3)))
      return
      end

!ws**********************************************************************
      subroutine find_OXpoint_backup(ww)
      USE DECLARE
      USE DECLARE_OXpoint

      real*8  wmax,wmin
      real*8, dimension(mx,mz) :: ww

      call find_maxmin(ww,nrkyo,wmax,wmin,jxto,jzto,jxtx,jztx)
!????
!      if(jxto .lt. mxt/2) then
!      jxto=jxto-1
!      else
!      jxto=jxto+1
!      endif
!
!      if(jztx .lt. mzt/2) then
!      jztx=jztx-1
!      else
!      jztx=jztx+1
!      endif
!????

      yy_O= yyt(nrkyo*mym+jyo-2)
      ps_O= pst(jxto,jzto)
      rr_O= rrt(jxto,jzto)
      tc_O=tcht(jxto,jzto)      
      tc_O=ANINT(tc_O/pi*myt)*pi/myt

      yy_X= yy_O
      ps_X= pst(jxtx,jztx)
      rr_X= rrt(jxtx,jztx)
      tc_X=tcht(jxtx,jztx)
      tc_X=ANINT(tc_X/pi*myt)*pi/myt
!      write(*,*) '*********************Opint*********************'
!      write(*,*) nrank,jxto,jzto,ps_O,rr_O,tc_O,yy_O
!
!      write(*,*) '*********************Xpint*********************'
!      write(*,*) nrank,jxtx,jztx,ps_X,rr_X,tc_X,yy_X

      if(nrank.eq.0) then
      write(311,1100) time,jxto,jzto,ps_O,rr_O,tc_O,jxtx,jztx,ps_X,rr_X,tc_X
1100  format((1x,e12.5,2(2i,3(1x,e12.5))))
      endif

      return
    end
    
!ws**********************************************************************
      subroutine find_OXpoint_1st
      USE DECLARE
      USE DECLARE_OXpoint

      real*8  wmax,wmin
      real*8, dimension(mx,mz) :: ww
      
       ww(:,:)=x1(:,:,3,3)*x(:,:,3,8)-x1(:,:,3,5)*x(:,:,3,6)
      call find_maxmin(ww,nrkyo,wmax,wmin,jxto,jzto,jxtx,jztx)
!????
!      if(jxto .lt. mxt/2) then
!      jxto=jxto-1
!      else
!      jxto=jxto+1
!      endif
!
!      if(jztx .lt. mzt/2) then
!      jztx=jztx-1
!      else
!      jztx=jztx+1
!      endif
!????

      yy_O= 0 
      ps_O= pst(jxto,jzto)
      rr_O= rrt(jxto,jzto) 
      tc_O= 0 

      yy_X= 0 !yy_O
      ps_X= psmode
      rr_X= rrmode
      tc_X= pi/2 !ANINT(tc_X/pi*myt)*pi/myt
!      write(*,*) '*********************Opint*********************'
!      write(*,*) nrank,jxto,jzto,ps_O,rr_O,tc_O,yy_O
!
!      write(*,*) '*********************Xpint*********************'
!      write(*,*) nrank,jxtx,jztx,ps_X,rr_X,tc_X,yy_X

      if(nrank.eq.0) then
      write(311,1100) time,jxto,jzto,ps_O,rr_O,tc_O,jxtx,jztx,ps_X,rr_X,tc_X
      
1100  format((1x,e12.5,2(2i,3(1x,e12.5))))
      endif

      return
    end
    
    !ws**********************************************************************
      subroutine find_OXpoint(ww)
      USE DECLARE
      USE DECLARE_OXpoint

      real*8  wmax,wmin
      real*8, dimension(mx,mz) :: ww

      call find_maxmin(ww,nrkyo,wmax,wmin,jxto,jzto,jxtx,jztx)

      if(pst(jxto,jzto) .gt. ps_O) then
      ps_O= pst(jxto,jzto)
      rr_O= rrt(jxto,jzto)
      endif

 
      if(nrank.eq.0) then
      write(311,1100) time,jxto,jzto,ps_O,rr_O,tc_O,jxtx,jztx,ps_X,rr_X,tc_X
1100  format((1x,e12.5,2(2i,3(1x,e12.5))))
      endif

      return
    end

