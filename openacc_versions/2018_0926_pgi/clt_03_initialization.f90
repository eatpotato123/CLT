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





!!! major parts for initialization in clt is  as follows

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!    !!!!!    !!!!!  !!!!!!    !!!!!!          !!!!
!!!!!  !!!!!!  !  !!!!  !!!!!!!  !!!!!!!!!!!  !!!!!!!!
!!!!!  !!!!!!  !!  !!!  !!!!!!!  !!!!!!!!!!!  !!!!!!!!
!!!!!  !!!!!!  !!!  !!  !!!!!!!  !!!!!!!!!!!  !!!!!!!!
!!!!!  !!!!!!  !!!!  !  !!!!!!!  !!!!!!!!!!!  !!!!!!!!
!!!!!  !!!!!!  !!!!!    !!!!!!!  !!!!!!!!!!!  !!!!!!!!
!!!!    !!!!!  !!!!!!   !!!!!!    !!!!!!!!!!  !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module for initial equilibrium and initialization.
! including two kinds of equilibriums (from NOVA or gfile of k-EFIT)
! the read_pgfile is added by H.W. Zhang during 2018 Mar-Jul.
! contact changhw@zju.edu.cn if there is any problem.

! including subroutines as follows
! input
! -input the options and parameter setting up
! print_input
! initia
! init_ps1
! perturbation
! read_nova
! read_nova_v2 (not used)
! read_nova_tm (not used)
! gridpnt
! map_nova
! mapeq_st2xz
! map_nova_v2 (not used)
! readmap_wch
! estimate_pst1
! last_grid
! read_pgfile                        
! -read equilibrium from k-EFIT gfile 
!   -(modified by MATLAB) 
! map_nova_to_bnd_grd (not used, added by H.W. Zhang during 2018 Mar-Jul.)
! map_nova_to_bnd_grd_in_each_proc (not used, added by H.W. Zhang during 2018 Mar-Jul.) 
! map_int_to_bnd_grd (added by H.W. Zhang during 2018 Mar-Jul.) 



!ws*************************************************************************
     subroutine input
! --------------------------
!  this routine inputs parameters and define basic variables
!   lrstrt: =.f., starting from t=0; =.t., continueing from
!           steps as given bz nst.
!   nend:   the final steps intended for the current run ,
!           including the steps from the previous run.
!   nsp:    time step interval for diagnostic plots.
!   nsmthx:  the number of rows starting from incoming
!           boundary for the smoothing region in r-direction.
!   eta:    exact inverse of magnetic renolds number.
!   gamma:  adiabatic constant.
!   beta:   ratio of kinetic to magnetic pressure
!   ncase:  case number of the run.
!   nst:    beginning data file number for the current run.
!   nintd:   data file number increment.
! --------------------------
!
      use declare
      include 'mpif.h'
!
      cfl=1.2
      cext1=0.5
      caf=0.75d0
      lpetsc=.false. !.false. 
      constp=.false. !.true.
      p00=0.0001
      rotation=.false.
      lrot=2
      constu=.false. !.true.
      uy0=0.00
      consto=.false. !.true.
      oy0=0.00

      lrstrt=.false.
      ncase=24
      nend=47
      nst=0
      nintd=1
!!mode      
      nmode=1
      mmode=2
      qmode=1.0*mmode/nmode
      qdgn(1)=qmode
!!mode

	initia_from_pgfile=.false.
      firstmap=.true.
      rshear=.false. !.true.
      spectral=.false. 
      filt=.true.
      smooth=.false. !.true.
      smoothc=.false. !.true.
      smoothx1=.false. !.true.
      smoothbout=.false. !.false. !
      smoothef=.false.
      smoothpll=.false. !.true.
      smoothp1ll=.false. !.true.
      invaryrho=.false. !.true.
      invaryp=.false. !.true.

	rho_unif=.false.
!      correct=.true.
      uniformx=.true.
      uniformz=.true.
!      halfx=.false.
      symmetryx=.false.

	if(initia_from_pgfile) then
      symmetryz=.false.
	else
      symmetryz=.true.
	zzero=0.d0
	endif

      analysis=.false. !.true.

      resisitive=.true.
      ef_mode=.true. !.false. 
      etaj_in_e=.true.  
      viscous=.true.
      hall=.false.
      pressure=.false.
      soundwave=.false. !.true.
      eta_from_t=.false. !.true.
      rho_from_p=.false. !.true.
	linear_mhd=.false.

      ohm_heat =.false. 
      implicitb=.false. !.true.
      implicitv=.false. !.true.  
	rmp_east=.false.
	rmp_phi_up=0.d0
	rmp_phi_low=0.d0
	i0_rmp=100 ! A ~ usually = kA order


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
      cd_oxpoint=.false.  
      br_rise=.false.
      lcd=0 !0:no curdrive; 1:cd in efield; 2:cd in current
      fcd=0.01

!!conductpll   --->  
      conductpll=.false.
      lpll=0    !1=pll_subloop   !2=pll_smthpline        !3=pll_soundwave   !4=pll_petsc        
      lscheme=1 !1=euler         !1=smthp_traceline      !1=euler
                !2=rk4           !2=smthp_traceline_v1   !2=rk4 
                !3=lax           !3=smthp_traceline_v2   !3=lax 
                !4=implict       !4=smthp2_traceline     !4=implict
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

	if(initia_from_pgfile) then
	call read_pgfile
	else
      if(.not.analysis) call read_nova
	endif

      call gridpnt
!!!!

      cd00=fcd*cip
      
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
      dt=0.
      nstep=0

      if(nrank==0) call print_input
!
      return
      end

!ws*************************************************************************
     subroutine print_input

      use declare
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
      write(*,*) 'ef_mode   =',ef_mode
      write(*,*) 'etaj_in_e =',etaj_in_e   
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
      write(*,*) 'zzero=',zzero
      write(*,*) 'xmg=',xmg
      write(*,*) 'zmg=',zmg
      write(*,*) 'xdim=',xdim
      write(*,*) 'zdim=',zdim
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
      write(19,*) 'ef_mode   =',ef_mode
      write(19,*) 'etaj_in_e =',etaj_in_e   
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
! defines coordinates system and specifies initial configuration.
! normlization convention:
!   1. density --- normalised to asymtotic value, i.e., rho=1
!   2. magnetic field --- normalised to asymtotic value, i.e.,
!                         b0=1.
!   3. velocity --- normalised to asymtotic alfven speed, va=1, a
!                   natural result of 1. and 2.
!   4. length --- normalised to a.
!   5. time --- normalised to a/va.
!---------------
!
     use declare
     include 'mpif.h'
     real*8 rhom,rhomp
     real*8 ciy,ciy_nrk,bp0max,bp0max1,bpc
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


      do jx=3,mxt-2
      dxp1=xxt(jx+1)-xxt(jx)
      dxm1=xxt(jx)-xxt(jx-1)
      dxp2=xxt(jx+2)-xxt(jx)
      dxm2=xxt(jx)-xxt(jx-2)
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
      axt2(jx)=2.*h3*(g2*h3-g3*h2)/ca1
      bxt2(jx)=2.*h3*(h2*f3-h3*f2)/ca1
      cxt2(jx)=2.*h3*(f2*g3-f3*g2)/ca1
      dxt2(jx)=-(dxp1*axt2(jx)+dxm1*bxt2(jx) &
         +dxp2*cxt2(jx))/dxm2
 	enddo
    
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
    
      do jx=3,mxt-2
      dxp1=xxt(jx+1)-xxt(jx)
      dxm1=xxt(jx)-xxt(jx-1)
      dxp2=xxt(jx+2)-xxt(jx)
      dxm2=xxt(jx)-xxt(jx-2)
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
      axt1(jx)=h3*(g2*h3-g3*h2)/ca1
      bxt1(jx)=h3*(h2*f3-h3*f2)/ca1
      cxt1(jx)=h3*(f2*g3-f3*g2)/ca1
      dxt1(jx)=(dxp1**2*axt1(jx)-dxm1**2*bxt1(jx) &
         +dxp2**2*cxt1(jx))/dxm2**2
 	enddo

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

      do jx=1,mxt-3
      dxp1=xxt(jx+1)-xxt(jx)
      dxp2=xxt(jx+2)-xxt(jx)
      dxp3=xxt(jx+3)-xxt(jx)
      f1=dxp1-dxp1**3/dxp3**2
      f2=dxp1**2-dxp1**3/dxp3
      g1=dxp2-dxp2**3/dxp3**2
      g2=dxp2**2-dxp2**3/dxp3
      ca1=f1*g2-f2*g1
      axtp(jx)=g2/ca1
      bxtp(jx)=-f2/ca1
      cxtp(jx)=(1-axtp(jx)*dxp1-bxtp(jx)*dxp2)/dxp3
	enddo

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

      do jx=4,mxt
      dxm1=xxt(jx)-xxt(jx-1)
      dxm2=xxt(jx)-xxt(jx-2)
      dxm3=xxt(jx)-xxt(jx-3)
      f1=dxm1-dxm1**3/dxm3**2
      f2=dxm1**2-dxm1**3/dxm3
      g1=dxm2-dxm2**3/dxm3**2
      g2=dxm2**2-dxm2**3/dxm3
      ca1=f1*g2-f2*g1
      axtm(jx)=g2/ca1
      bxtm(jx)=-f2/ca1
      cxtm(jx)=(1-axtm(jx)*dxm1-bxtm(jx)*dxm2)/dxm3
      enddo

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

      do jx=2,mxt-2
      dxm1=xxt(jx)-xxt(jx-1)
      dxp1=xxt(jx+1)-xxt(jx)
      dxp2=xxt(jx+2)-xxt(jx)
      f1=-dxm1+dxm1**3/dxp2**2
      f2=dxm1**2+dxm1**3/dxp2
      g1=dxp1-dxp1**3/dxp2**2
      g2=dxp1**2-dxp1**3/dxp2
      ca1=f1*g2-f2*g1
      axtbp(jx)=g2/ca1
      bxtbp(jx)=-f2/ca1
      cxtbp(jx)=(1+axtbp(jx)*dxm1-bxtbp(jx)*dxp1)/dxp2
      enddo

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

      do jx=3,mxt-1
      dxp1=xxt(jx+1)-xxt(jx)
      dxm1=xxt(jx)-xxt(jx-1)
      dxm2=xxt(jx)-xxt(jx-2)
      f1=-dxp1+dxp1**3/dxm2**2
      f2=dxp1**2+dxp1**3/dxm2
      g1=dxm1-dxm1**3/dxm2**2
      g2=dxm1**2-dxm1**3/dxm2
      ca1=f1*g2-f2*g1
      axtbm(jx)=g2/ca1
      bxtbm(jx)=-f2/ca1
      cxtbm(jx)=(1+axtbm(jx)*dxp1-bxtbm(jx)*dxm1)/dxm2
      enddo

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


      do jz=3,mzt-2
      dzp1=zzt(jz+1)-zzt(jz)
      dzm1=zzt(jz)-zzt(jz-1)
      dzp2=zzt(jz+2)-zzt(jz)
      dzm2=zzt(jz)-zzt(jz-2)
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
      azt2(jz)=2.*h3*(g2*h3-g3*h2)/ca1
      bzt2(jz)=2.*h3*(h2*f3-h3*f2)/ca1
      czt2(jz)=2.*h3*(f2*g3-f3*g2)/ca1
      dzt2(jz)=-(dzp1*azt2(jz)+dzm1*bzt2(jz) &
         +dzp2*czt2(jz))/dzm2
 	enddo
    
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
    
      do jz=3,mzt-2
      dzp1=zzt(jz+1)-zzt(jz)
      dzm1=zzt(jz)-zzt(jz-1)
      dzp2=zzt(jz+2)-zzt(jz)
      dzm2=zzt(jz)-zzt(jz-2)
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
      azt1(jz)=h3*(g2*h3-g3*h2)/ca1
      bzt1(jz)=h3*(h2*f3-h3*f2)/ca1
      czt1(jz)=h3*(f2*g3-f3*g2)/ca1
      dzt1(jz)=(dzp1**2*azt1(jz)-dzm1**2*bzt1(jz) &
         +dzp2**2*czt1(jz))/dzm2**2
 	enddo

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

      do jz=1,mzt-3
      dzp1=zzt(jz+1)-zzt(jz)
      dzp2=zzt(jz+2)-zzt(jz)
      dzp3=zzt(jz+3)-zzt(jz)
      f1=dzp1-dzp1**3/dzp3**2
      f2=dzp1**2-dzp1**3/dzp3
      g1=dzp2-dzp2**3/dzp3**2
      g2=dzp2**2-dzp2**3/dzp3
      ca1=f1*g2-f2*g1
      aztp(jz)=g2/ca1
      bztp(jz)=-f2/ca1
      cztp(jz)=(1-aztp(jz)*dzp1-bztp(jz)*dzp2)/dzp3
	enddo

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

      do jz=4,mzt
      dzm1=zzt(jz)-zzt(jz-1)
      dzm2=zzt(jz)-zzt(jz-2)
      dzm3=zzt(jz)-zzt(jz-3)
      f1=dzm1-dzm1**3/dzm3**2
      f2=dzm1**2-dzm1**3/dzm3
      g1=dzm2-dzm2**3/dzm3**2
      g2=dzm2**2-dzm2**3/dzm3
      ca1=f1*g2-f2*g1
      aztm(jz)=g2/ca1
      bztm(jz)=-f2/ca1
      cztm(jz)=(1-aztm(jz)*dzm1-bztm(jz)*dzm2)/dzm3
	enddo


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

      do jz=2,mzt-2
      dzm1=zzt(jz)-zzt(jz-1)
      dzp1=zzt(jz+1)-zzt(jz)
      dzp2=zzt(jz+2)-zzt(jz)
      f1=-dzm1+dzm1**3/dzp2**2
      f2=dzm1**2+dzm1**3/dzp2
      g1=dzp1-dzp1**3/dzp2**2
      g2=dzp1**2-dzp1**3/dzp2
      ca1=f1*g2-f2*g1
      aztbp(jz)=g2/ca1
      bztbp(jz)=-f2/ca1
      cztbp(jz)=(1+aztbp(jz)*dzm1-bztbp(jz)*dzp1)/dzp2
	enddo

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

      do jz=3,mzt-1
      dzp1=zzt(jz+1)-zzt(jz)
      dzm1=zzt(jz)-zzt(jz-1)
      dzm2=zzt(jz)-zzt(jz-2)
      f1=-dzp1+dzp1**3/dzm2**2
      f2=dzp1**2+dzp1**3/dzm2
      g1=dzm1-dzm1**3/dzm2**2
      g2=dzm1**2-dzm1**3/dzm2
      ca1=f1*g2-f2*g1
      aztbm(jz)=g2/ca1
      bztbm(jz)=-f2/ca1
      cztbm(jz)=(1+aztbm(jz)*dzp1-bztbm(jz)*dzm1)/dzm2
	enddo

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

!      do jy=3,myt-2
!!      do 31 jy=1,my
!      dyp1=yyt(jy+1)-yyt(jy)
!      dym1=yyt(jy)-yyt(jy-1)
!      dyp2=yyt(jy+2)-yyt(jy)
!      dym2=yyt(jy)-yyt(jy-2)
!      f1= dyp1**2+dyp1*dym2
!      f2= dyp1**3-dyp1*dym2**2
!      f3= dyp1**4+dyp1*dym2**3
!      g1=-dym1**2+dym1*dym2
!      g2= dym1**3-dym1*dym2**2
!      g3=-dym1**4+dym1*dym2**3
!      h1= dyp2**2+dyp2*dym2
!      h2= dyp2**3-dyp2*dym2**2
!      h3= dyp2**4+dyp2*dym2**3
!      ca1= (f1*h3-f3*h1)*(g2*h3-g3*h2)-(g1*h3-g3*h1)*(f2*h3-f3*h2)
!      ayt2(jy)=2.*h3*(g2*h3-g3*h2)/ca1
!      byt2(jy)=2.*h3*(h2*f3-h3*f2)/ca1
!      cyt2(jy)=2.*h3*(f2*g3-f3*g2)/ca1
!      dyt2(jy)=-(dyp1*ayt2(jy)+dym1*byt2(jy) &
!         +dyp2*cyt2(jy))/dym2
! 	enddo
    
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

!      do jy=3,myt-2
!!      do 32 jy=1,my
!      dyp1=yyt(jy+1)-yyt(jy)
!      dym1=yyt(jy)-yyt(jy-1)
!      dyp2=yyt(jy+2)-yyt(jy)
!      dym2=yyt(jy)-yyt(jy-2)
!      f1= dyp1   +dyp1**4/dym2**3
!      f2= dyp1**2-dyp1**4/dym2**2
!      f3= dyp1**3+dyp1**4/dym2
!      g1= dym1   -dym1**4/dym2**3
!      g2=-dym1**2+dym1**4/dym2**2
!      g3= dym1**3-dym1**4/dym2
!      h1= dyp2   +dyp2**4/dym2**3
!      h2= dyp2**2-dyp2**4/dym2**2
!      h3= dyp2**3+dyp2**4/dym2
!      ca1= (f1*h3-f3*h1)*(g2*h3-g3*h2)-(g1*h3-g3*h1)*(f2*h3-f3*h2)
!      ayt1(jy)=h3*(g2*h3-g3*h2)/ca1
!      byt1(jy)=h3*(h2*f3-h3*f2)/ca1
!      cyt1(jy)=h3*(f2*g3-f3*g2)/ca1
!      dyt1(jy)=(dyp1**2*ayt1(jy)-dym1**2*byt1(jy) &
!         +dyp2**2*cyt1(jy))/dym2**2
! 	enddo

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
	
	if(initia_from_pgfile) then
	ux(:,:)=0.d0
	uz(:,:)=0.d0
	uxdx(:,:)=0.d0
	uzdx(:,:)=0.d0
	uxdz(:,:)=0.d0
	uzdz(:,:)=0.d0
	endif

	if(.not.initia_from_pgfile) call map_nova

      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
        xint(jx,jz,1)=rh(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        xint(jx,jz,2)=pt(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2) !-p_nova(mpsa)
	if(.not.initia_from_pgfile) then
        xint(jx,jz,3)=0
	else
        xint(jx,jz,3)=ux(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
	endif
        xint(jx,jz,4)=uy(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
	if(.not.initia_from_pgfile) then
        xint(jx,jz,5)=0
	else
        xint(jx,jz,5)=uz(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
	endif
        xint(jx,jz,6)=bx(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        xint(jx,jz,8)=bz(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        xint(jx,jz,7)=by(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)

        cint(jx,jz,1)=cx(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        cint(jx,jz,3)=cz(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        cint(jx,jz,2)=cy(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)

        xint_dx(jx,jz,1)=rhdx(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        xint_dx(jx,jz,2)=ptdx(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
	if(.not.initia_from_pgfile) then
        xint_dx(jx,jz,3)=0
	else	  
        xint_dx(jx,jz,3)=uxdx(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
	endif	  
        xint_dx(jx,jz,4)=uydx(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
	if(.not.initia_from_pgfile) then
        xint_dx(jx,jz,5)=0
	else
        xint_dx(jx,jz,5)=uzdx(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
	endif
        xint_dx(jx,jz,6)=bxdx(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        xint_dx(jx,jz,8)=bzdx(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        xint_dx(jx,jz,7)=bydx(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)

        cint_dx(jx,jz,1)=cx_dx(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        cint_dx(jx,jz,3)=cz_dx(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        cint_dx(jx,jz,2)=cy_dx(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)

        xint_dz(jx,jz,1)=rhdz(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        xint_dz(jx,jz,2)=ptdz(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
	if(.not.initia_from_pgfile) then
        xint_dz(jx,jz,3)=0
	else	  
        xint_dz(jx,jz,3)=uxdz(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
	endif
        xint_dz(jx,jz,4)=uydz(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
	if(.not.initia_from_pgfile) then
        xint_dz(jx,jz,5)=0
	else	  
        xint_dz(jx,jz,5)=uzdz(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
	endif
        xint_dz(jx,jz,6)=bxdz(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        xint_dz(jx,jz,8)=bzdz(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        xint_dz(jx,jz,7)=bydz(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)

        cint_dz(jx,jz,1)=cx_dz(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        cint_dz(jx,jz,3)=cz_dz(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        cint_dz(jx,jz,2)=cy_dz(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)

        oyint(jx,jz)=omrot(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        oypint(jx,jz)=omprot(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)


	  
	  if(initia_from_pgfile) then
        psi(jx,jz)   =   pst(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        psi_dx(jx,jz)=pst_dx(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        psi_dz(jx,jz)=pst_dz(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        rr2(jx,jz)   =  rr2t(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        rr(jx,jz)    =   rrt(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        thxz(jx,jz)  =   tht(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2)
        tpxz(jx,jz)  =atan2(psi_dz(jx,jz),psi_dx(jx,jz))
        if(zz(jz).lt.0) tpxz(jx,jz)=tpxz(jx,jz)+2*pi
        endif

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

! xint1-8, cint1-3 for openacc use
	xint1(:,:)=xint(:,:,1)
	xint2(:,:)=xint(:,:,2)
	xint3(:,:)=xint(:,:,3)
	xint4(:,:)=xint(:,:,4)
	xint5(:,:)=xint(:,:,5)
	xint6(:,:)=xint(:,:,6)
	xint7(:,:)=xint(:,:,7)
	xint8(:,:)=xint(:,:,8)

	xint1_dx(:,:)=xint_dx(:,:,1)
	xint2_dx(:,:)=xint_dx(:,:,2)
	xint3_dx(:,:)=xint_dx(:,:,3)
	xint4_dx(:,:)=xint_dx(:,:,4)
	xint5_dx(:,:)=xint_dx(:,:,5)
	xint6_dx(:,:)=xint_dx(:,:,6)
	xint7_dx(:,:)=xint_dx(:,:,7)
	xint8_dx(:,:)=xint_dx(:,:,8)

	xint1_dz(:,:)=xint_dz(:,:,1)
	xint2_dz(:,:)=xint_dz(:,:,2)
	xint3_dz(:,:)=xint_dz(:,:,3)
	xint4_dz(:,:)=xint_dz(:,:,4)
	xint5_dz(:,:)=xint_dz(:,:,5)
	xint6_dz(:,:)=xint_dz(:,:,6)
	xint7_dz(:,:)=xint_dz(:,:,7)
	xint8_dz(:,:)=xint_dz(:,:,8)


	cint1(:,:)=cint(:,:,1)
	cint2(:,:)=cint(:,:,2)
	cint3(:,:)=cint(:,:,3)

	cint1_dx(:,:)=cint_dx(:,:,1)
	cint2_dx(:,:)=cint_dx(:,:,2)
	cint3_dx(:,:)=cint_dx(:,:,3)

	cint1_dz(:,:)=cint_dz(:,:,1)
	cint2_dz(:,:)=cint_dz(:,:,2)
	cint3_dz(:,:)=cint_dz(:,:,3)
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
! hw add if to avoid divided by zero
      
      if(xint(jx,jz,1).ne.0.d0) then
        tmint(jx,jz)=xint(jx,jz,2)/xint(jx,jz,1)
        tmint_dx(jx,jz)=xint_dx(jx,jz,2)/xint(jx,jz,1)-xint(jx,jz,2)*xint_dx(jx,jz,1)/xint(jx,jz,1)**2
        tmint_dz(jx,jz)=xint_dz(jx,jz,2)/xint(jx,jz,1)-xint(jx,jz,2)*xint_dz(jx,jz,1)/xint(jx,jz,1)**2

        cj(jx,jz)=cint(jx,jz,2)   
        bp0(jx,jz)=sqrt(xint(jx,jz,6)**2+xint(jx,jz,8)**2)
        bb0(jx,jz)=sqrt(xint(jx,jz,6)**2+xint(jx,jz,7)**2+xint(jx,jz,8)**2)
      endif

      if(bp0(jx,jz).ne.0.d0) then
        wx2r(jx,jz)=-xint(jx,jz,8)/bp0(jx,jz)
        wz2r(jx,jz)=xint(jx,jz,6)/bp0(jx,jz)
        wx2p(jx,jz)=-xint(jx,jz,6)/bp0(jx,jz)
        wz2p(jx,jz)=-xint(jx,jz,8)/bp0(jx,jz) 

!        wx2r(jx,jz)=dcos(tpxz(jx,jz))
!        wz2r(jx,jz)=dsin(tpxz(jx,jz))
!        wx2p(jx,jz)=-dsin(tpxz(jx,jz))
!        wz2p(jx,jz)=dcos(tpxz(jx,jz))
      endif
      enddo
      enddo

!ws:avoid bp~0: bp0->bp0c
      bp0max1=maxval(bp0)
      call mpi_allreduce(bp0max1,bp0max,1,mpi_double_precision,mpi_max, &
                        mpi_comm_world,ierror)

      bpc=0.05*bp0max
      bp0c=bp0
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
      if(rr(jx,jz).lt.0.3 .and. bp0(jx,jz).lt.bpc) then
! hw add if to avoid divided by zero
      if(bpc*bp0c(jx,jz).ne.0.d0) then
      bp0c(jx,jz)=bpc/2+0.5*bp0(jx,jz)**2/bpc
        wx2r(jx,jz)=-xint(jx,jz,8)/bp0c(jx,jz)
        wz2r(jx,jz)=xint(jx,jz,6)/bp0c(jx,jz)
        wx2p(jx,jz)=-xint(jx,jz,6)/bp0c(jx,jz)
        wz2p(jx,jz)=-xint(jx,jz,8)/bp0c(jx,jz) 
      endif
      endif
      enddo
      enddo
!ws:avoid bp~0

!compute i_phi++++++++++
      ciy_nrk=0
      do jz=iz_first+2,iz_last-2
      do jx=ix_first+2,ix_last-2
      if(psi(jx,jz) .lt. psia) then 
      ciy_nrk=ciy_nrk+cj(jx,jz)*(xx(jx+1)-xx(jx-1))*(zz(jz+1)-zz(jz-1))/4
      endif
      enddo
      enddo

      call mpi_allreduce(ciy_nrk,ciy,1,mpi_double_precision,mpi_sum, &
                        mpi_comm_world,ierror)
      if(nrank==0) write(*,*) 'ip=',cip,'iy=',ciy
!compute i_phi-----------

      tm01=maxval(tmint)
      call mpi_allreduce(tm01,tm00,1,mpi_double_precision,mpi_max, &
                        mpi_comm_world,ierror)  
                          
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


      x1(:,:,:,:)=0.d0
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
      if(psi(jx,jz) .lt. psival_nova(mpsa-2)) then
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
      call mpi_reduce(fxmax0,fxmax1,1,mpi_double_precision,mpi_max,0, &
            mpi_comm_world,ierror)
      call mpi_reduce(fxmin0,fxmin1,1,mpi_double_precision,mpi_min,0, &
            mpi_comm_world,ierror)

      call mpi_reduce(fzmax0,fzmax1,1,mpi_double_precision,mpi_max,0, &
            mpi_comm_world,ierror)
      call mpi_reduce(fzmin0,fzmin1,1,mpi_double_precision,mpi_min,0, &
            mpi_comm_world,ierror)

      call mpi_reduce(fymax0,fymax1,1,mpi_double_precision,mpi_max,0, &
            mpi_comm_world,ierror)
      call mpi_reduce(fymin0,fymin1,1,mpi_double_precision,mpi_min,0, &
            mpi_comm_world,ierror)

      call mpi_reduce(dbmax0,dbmax1,1,mpi_double_precision,mpi_max,0, &
            mpi_comm_world,ierror)
      call mpi_reduce(dbmin0,dbmin1,1,mpi_double_precision,mpi_min,0, &
            mpi_comm_world,ierror)
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
!mpi   -----------------------------------------------------------------      
!      call mpi_reduce(ft,ft1,1,mpi_double_precision,mpi_sum,0, &
!                       mpi_comm_world,ierror)
!      call mpi_reduce(gt,gt1,1,mpi_double_precision,mpi_sum,0, &
!                       mpi_comm_world,ierror)
!      call mpi_reduce(ht,ht1,1,mpi_double_precision,mpi_sum,0, &
!                       mpi_comm_world,ierror)
!mpi   -----------------------------------------------------------------

!      if(nrank.eq.0) then
!      open(unit=17,file='energy_init.dat',status='unknown',form='formatted')
!      write(17,*) "magnetic,kinetic,heat,total" 
!      write(17,*) "eqm:"
!      write(17,400) ft1,gt1,ht1,ft1+gt1+ht1
!400   format(4(1x,e12.5))
!      endif
!      ft0=ft1
!      gt0=gt1
!      ht0=ht1

      call perturbation

!	if(rmp_east) call rmp_east_pert

      return
      end

!ws********************************************
     subroutine init_ps1
      use declare
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
	  
!ws**************************************************************************
     subroutine perturbation
      use declare
!      real*8, dimension(mx,mz,my,3) :: eta1j
      real*8 drey_dx,dey_dx,dez_dx,dex_dy,dez_dy,dex_dz,dey_dz,per_vp !,eta1
      include 'mpif.h'
!  d1fc= d f / dx  with fourth-order accuracy central difference
      d1fc(fm2,fm1,f0,fp1,fp2,a,b,c,d)= &
       a*(fp1-f0)+b*(f0-fm1)+c*(fp2-f0)+d*(f0-fm2)
!  d1f2= d f / dx  with second-order accuracy central difference
      d1f2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(fp1-f0) &
         -(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)
!  d1xf2= d rf / dx  with second-order accuracy central difference
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
!	psmode=0.050299150612909
!	psmode=0.070146296723848
!      eta1(jx,jz,jy)=eta10*(exp(-(psi(jx,jz)-psmode)**2/deltax**2))*dcos(nmode*yy(jy)+mmode*thxz(jx,jz))
!	psmode=0.262538435913944
!	psmode=0.206509429951520
!      eta1(jx,jz,jy)=eta1(jx,jz,jy)+eta10*(exp(-(psi(jx,jz)-psmode)**2/deltax**2))*dcos(nmode*(yy(jy)-pi)+mmode*thxz(jx,jz))
	if(initia_from_pgfile) then
!	psmode=0.494853880167460
	psmode=0.527707001845421
      eta1(jx,jz,jy)=eta1(jx,jz,jy)+eta10*(exp(-(psi(jx,jz)-psmode)**2/deltax**2))*dcos(nmode*(yy(jy)-pi)+mmode*thxz(jx,jz))
	else
      eta1(jx,jz,jy)=eta10*(exp(-(psi(jx,jz)-psmode)**2/deltax**2))*dcos(nmode*(yy(jy)-pi)+mmode*thxz(jx,jz))
	endif
!      if(rshear) then
!      eta1(jx,jz,jy)=eta10*(exp(-(psi(jx,jz)-psmode)**2/deltax**2)+exp(-(psi(jx,jz)-ps1mode)**2/deltax**2))*((1+dcos(nmode*yy(jy)+mmode*thxz(jx,jz)))/2)**4
!      else
!      eta1(jx,jz,jy)=eta10*(exp(-(psi(jx,jz)-psmode)**2/deltax**2))*((1+dcos(nmode*yy(jy)+mmode*thxz(jx,jz)))/2)**4
!      endif
      do m=1,3
      eta1j(jx,jz,jy,m)=eta1(jx,jz,jy)*cint(jx,jz,m)
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

      drey_dx=xx(jx)*d1fc(eta1j(jx-2,jz,jy,2),eta1j(jx-1,jz,jy,2),eta1j(jx,jz,jy,2),eta1j(jx+1,jz,jy,2), &
                eta1j(jx+2,jz,jy,2),ax1(jx),bx1(jx),cx1(jx),dx1(jx))+eta1j(jx,jz,jy,2)
!      dey_dx =d1f2(eta1j(jx-1,jz,jy,2),eta1j(jx,jz,jy,2),eta1j(jx+1,jz,jy,2),xx(jx-1),xx(jx),xx(jx+1))
      dez_dx =d1fc(eta1j(jx-2,jz,jy,3),eta1j(jx-1,jz,jy,3),eta1j(jx,jz,jy,3),eta1j(jx+1,jz,jy,3),eta1j(jx+2,jz,jy,3) &
                ,ax1(jx),bx1(jx),cx1(jx),dx1(jx))
      dex_dz =d1fc(eta1j(jx,jz-2,jy,1),eta1j(jx,jz-1,jy,1),eta1j(jx,jz,jy,1),eta1j(jx,jz+1,jy,1),eta1j(jx,jz+2,jy,1) &
                ,az1(jz),bz1(jz),cz1(jz),dz1(jz))
      dey_dz =d1fc(eta1j(jx,jz-2,jy,2),eta1j(jx,jz-1,jy,2),eta1j(jx,jz,jy,2),eta1j(jx,jz+1,jy,2),eta1j(jx,jz+2,jy,2) &
                ,az1(jz),bz1(jz),cz1(jz),dz1(jz))
      dex_dy =d1fc(eta1j(jx,jz,jy-2,1),eta1j(jx,jz,jy-1,1),eta1j(jx,jz,jy,1),eta1j(jx,jz,jy+1,1),eta1j(jx,jz,jy+2,1) &
                ,ay1(jy),by1(jy),cy1(jy),dy1(jy))
      dez_dy =d1fc(eta1j(jx,jz,jy-2,1),eta1j(jx,jz,jy-1,3),eta1j(jx,jz,jy,3),eta1j(jx,jz,jy+1,3),eta1j(jx,jz,jy+2,1) &
                ,ay1(jy),by1(jy),cy1(jy),dy1(jy))

!      drey_dx=d1xf2(eta1j(jx-1,jz,jy,2),eta1j(jx,jz,jy,2),eta1j(jx+1,jz,jy,2),xx(jx-1),xx(jx),xx(jx+1))
!!      dey_dx =d1f2(eta1j(jx-1,jz,jy,2),eta1j(jx,jz,jy,2),eta1j(jx+1,jz,jy,2),xx(jx-1),xx(jx),xx(jx+1))
!      dez_dx =d1f2(eta1j(jx-1,jz,jy,3),eta1j(jx,jz,jy,3),eta1j(jx+1,jz,jy,3),xx(jx-1),xx(jx),xx(jx+1))
!      dex_dz =d1f2(eta1j(jx,jz-1,jy,1),eta1j(jx,jz,jy,1),eta1j(jx,jz+1,jy,1),zz(jz-1),zz(jz),zz(jz+1))
!      dey_dz =d1f2(eta1j(jx,jz-1,jy,2),eta1j(jx,jz,jy,2),eta1j(jx,jz+1,jy,2),zz(jz-1),zz(jz),zz(jz+1))
!      dex_dy =d1f2(eta1j(jx,jz,jy-1,1),eta1j(jx,jz,jy,1),eta1j(jx,jz,jy+1,1),yy(jy-1),yy(jy),yy(jy+1))
!      dez_dy =d1f2(eta1j(jx,jz,jy-1,3),eta1j(jx,jz,jy,3),eta1j(jx,jz,jy+1,3),yy(jy-1),yy(jy),yy(jy+1))
          
      perb(jx,jz,jy,1)=-dez_dy/xx(jx)+dey_dz  
      perb(jx,jz,jy,2)=-dex_dz+dez_dx
!      perb(jx,jz,jy,3)=(dex_dy-eta1j(jx,jz,jy,2))/xx(jx)-dey_dx
      perb(jx,jz,jy,3)=(dex_dy-drey_dx)/xx(jx)

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

!mpi   -----------------------------------------------------------------      
!      call mpi_reduce(ft,ft1,1,mpi_double_precision,mpi_sum,0, &
!                       mpi_comm_world,ierror)
!      call mpi_reduce(gt,gt1,1,mpi_double_precision,mpi_sum,0, &
!                       mpi_comm_world,ierror)
!      call mpi_reduce(ht,ht1,1,mpi_double_precision,mpi_sum,0, &
!                       mpi_comm_world,ierror)
!mpi   -----------------------------------------------------------------
!      if(nrank.eq.0) then
!      timeold=time
!      gtold=gt1
!      open(unit=17,file='energy_init.dat',status='unknown',form='formatted') 
!      write(17,*) "perb:"
!      write(17,400) ft1-ft0,gt1-gt0,ht1-ht0,ft1+gt1+ht1-ft0-gt0-ht0
!400   format(4(1x,e12.5))
!      close(17)
!      endif
      return
      end

!ws*******************************************************************************
     subroutine read_nova
     use declare 
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
      '   ip  ipnorm' )

      read(888,4000) xplmin,xplmax,zplmax,aguess,xzero,xmag,xmaj,xzmax,xatpi,xofset,aratio,bzero,curtotal,curnorm
      
      aa=aguess    
      aa=1.d0
      b0=1
      xzero=xzero/aa
      xmg=xmag/aa
      xmin=xplmin/aa
      xmax=xplmax/aa
      zmax=zplmax/aa
      zmin=-zmax

	xdim=(xmax-xmin)/2.d0
	zdim=(zmax-zmin)/2.d0
	

	xmin=xmin-(xmax-xmin)/mxt
	xmax=xmax+(xmax-xmin)/mxt
	zmin=zmin-(zmax-zmin)/mzt
	zmax=zmax+(zmax-zmin)/mzt

      zmg=0.0
      cip=curnorm*xzero*b0
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
      
      th_nova(jd)=thst(i)
      xx_nova(jd)=xxst(i,j)
      zz_nova(jd)=zzst(i,j)
      ps_nova(jd)=psst(i,j)            
      bx_nova(jd)=bxst(i,j)
      bz_nova(jd)=bzst(i,j)
      bxdx_nova(jd)=bxdxst(i,j)
      bxdz_nova(jd)=bxdzst(i,j)
      bzdx_nova(jd)=bzdxst(i,j)
      bzdz_nova(jd)=bzdzst(i,j)

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
      th_nova(jdm)=thst(im)
      xx_nova(jdm)=xxst(im,j)
      zz_nova(jdm)=zzst(im,j)
      ps_nova(jdm)=psst(im,j)
      bx_nova(jdm)=bxst(im,j)
      bz_nova(jdm)=bzst(im,j)
      bxdx_nova(jdm)=bxdxst(im,j)
      bxdz_nova(jdm)=bxdzst(im,j)
      bzdx_nova(jdm)=bzdxst(im,j)
      bzdz_nova(jdm)=bzdzst(im,j)
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
! to obtain the value at ipi=103, pi, which is missed in nova data	
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
      read(889,2100) jj,psival_nova(j),q_nova(j),qp_nova(j),p_nova(j),pp_nova(j),g_nova(j),gp_nova(j),f_nova(j),fp_nova(j),fb_nova(j),fbp_nova(j)
      psival_nova(j)=psival_nova(j)/(b0*aa**2)
      p_nova(j)=p_nova(j)/(b0**2)
      g_nova(j)=g_nova(j)/b0
      qp_nova(j)=qp_nova(j)*(b0*aa**2)
      pp_nova(j)=pp_nova(j)/(b0**2)*(b0*aa**2)
      gp_nova(j)=gp_nova(j)/b0*(b0*aa**2)
   41 continue
      close(889) 
      omrot_nova(:)=0
      omprot_nova(:)=0   
 
      else
      open(889,file='q_p_g.dat')
      do 40 j=1,npsi
      read(889,2000) jj,psival_nova(j),q_nova(j),qp_nova(j),p_nova(j),pp_nova(j),g_nova(j),gp_nova(j),f_nova(j),fp_nova(j),fb_nova(j),fbp_nova(j),omrot_nova(j),omprot_nova(j)     
      psival_nova(j)=psival_nova(j)/(b0*aa**2)
      p_nova(j)=p_nova(j)/(b0**2)
      g_nova(j)=g_nova(j)/b0
      qp_nova(j)=qp_nova(j)*(b0*aa**2)
      pp_nova(j)=pp_nova(j)/(b0**2)*(b0*aa**2)
      gp_nova(j)=gp_nova(j)/b0*(b0*aa**2)
      omrot_nova(j)=omrot_nova(j)/(b0*aa)
      omprot_nova(j)=omprot_nova(j)/(b0*aa)*(b0*aa**2)    
     
   40 continue
      close(889)      
      endif

      xxst(:,1)=xmg
      zzst(:,1)=0
      psst(:,1)=psival_nova(1)
      bxst(:,1)=0
      bzst(:,1)=0
      tst(:,1)=tst(:,2)
      rst(:,1)=0
 
       psia =psival_nova(mpsa)
       psmin=minval(ps_nova)
       psmax=maxval(ps_nova)

       qmin=minval(q_nova)
       qmax=maxval(q_nova)
       q0=q_nova(1)
 
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
      th12_nova(jd)=thst(i)
      xx12_nova(jd)=xxst(i,j)
      zz12_nova(jd)=zzst(i,j)

      th34_nova(jd)=thst(i+ipi-3)
      xx34_nova(jd)=xxst(i+ipi-3,j)
      zz34_nova(jd)=zzst(i+ipi-3,j)
      enddo
      enddo


      do j=1,npsi
      do i=1,n2th+5
      byst(i,j)=xzero*g_nova(j)/xxst(i,j)
      bydxst(i,j)=-xzero*(gp_nova(j)*bzst(i,j)+g_nova(j)/xxst(i,j)**2)
      bydzst(i,j)=xzero*gp_nova(j)*bxst(i,j)
      pdxst(i,j)=-pp_nova(j)*bzst(i,j)*xxst(i,j)
      pdzst(i,j)=pp_nova(j)*bxst(i,j)*xxst(i,j)

      cxst(i,j)=-xzero*gp_nova(j)*bxst(i,j)
      czst(i,j)=-xzero*gp_nova(j)*bzst(i,j)
      cyst(i,j)=bxdzst(i,j)-bzdxst(i,j)
!      cyst(i,j)=xxst(i,j)*pp_nova(j)+xzero**2*g_nova(j)*gp_nova(j)/xxst(i,j)

      uyst(i,j)=omrot_nova(j)*xxst(i,j)
      uydxst(i,j)=-omprot_nova(j)*bzst(i,j)*xxst(i,j)**2+omrot_nova(j)
      uydzst(i,j)= omprot_nova(j)*bxst(i,j)*xxst(i,j)**2

      if(j.ge.2 .and. i.ge.3 .and. i.le.n2th+2 .and. i.ne.ipi) then
      if(i.lt.ipi) jd=(j-2)*(n2th-1)+i-2
      if(i.gt.ipi) jd=(j-2)*(n2th-1)+i-3
      by_nova(jd)=byst(i,j)
      bydx_nova(jd)=bydxst(i,j)
      bydz_nova(jd)=bydzst(i,j)
      pdx_nova(jd)=pdxst(i,j)
      pdz_nova(jd)=pdzst(i,j)
      cy_nova(jd)=cyst(i,j)
      cx_nova(jd)=cxst(i,j)
      cz_nova(jd)=czst(i,j)

      uy_nova(jd)=uyst(i,j)
      uydx_nova(jd)=uydxst(i,j)
      uydz_nova(jd)=uydzst(i,j)

      rh_nova(jd)=rhom(ps_nova(jd))
      rhdx_nova(jd)=-rhomp(ps_nova(jd))*bzst(i,j)*xxst(i,j)
      rhdz_nova(jd)=rhomp(ps_nova(jd))*bxst(i,j)*xxst(i,j)
      
      pt_nova(jd)=p_nova(j)
      ptdx_nova(jd)=pdx_nova(jd)
      ptdz_nova(jd)=pdz_nova(jd)
      
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

      open(unit=101,file='fst_nova.dat',status='unknown',form='formatted')
      write(101,300)(((fffst(i,j,m),m=1,3),i=3,nthe+1),j=2,npsi)
 300  format(3(1x,e12.5))

      endif

!       psmin=minval(ps_nova)
!       psmax=maxval(ps_nova)

!       qmin=minval(q_nova)
!       qmax=maxval(q_nova)
!       q0=q_nova(1)
!      xmin=minval(xx_nova)
!      xmax=maxval(xx_nova)
!      zmin=minval(zz_nova)
!      zmax=maxval(zz_nova)
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
      write(890,5000) (xx_nova(j),zz_nova(j),th_nova(j),ps_nova(j),bx_nova(j),bxdx_nova(j),bxdz_nova(j),bz_nova(j),bzdx_nova(j),bzdz_nova(j),j=1,ndat)
 5000 format(10(1x,e17.9))
      close(890)
      open(891,file='stdat.dat')
      write(891,5100) ((thst(i),psst(i,j),xxst(i,j),zzst(i,j),tpst(i,j),tst(i,j),rst(i,j),i=1,n2th+5),j=1,npsi)
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
     use declare 
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
      read(889) psival_nova,q_nova,qp_nova,p_nova,pp_nova,g_nova,gp_nova,f_nova,fp_nova,fb_nova,fbp_nova,omrot_nova,omprot_nova     
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
      th_nova(jd)=thst(i)
      xx_nova(jd)=xxst(i,j)
      zz_nova(jd)=zzst(i,j)
      ps_nova(jd)=psst(i,j)            
      bx_nova(jd)=bxst(i,j)
      bz_nova(jd)=bzst(i,j)
      bxdx_nova(jd)=bxdxst(i,j)
      bxdz_nova(jd)=bxdzst(i,j)
      bzdx_nova(jd)=bzdxst(i,j)
      bzdz_nova(jd)=bzdzst(i,j)

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
      th_nova(jdm)=thst(im)
      xx_nova(jdm)=xxst(im,j)
      zz_nova(jdm)=zzst(im,j)
      ps_nova(jdm)=psst(im,j)
      bx_nova(jdm)=bxst(im,j)
      bz_nova(jdm)=bzst(im,j)
      bxdx_nova(jdm)=bxdxst(im,j)
      bxdz_nova(jdm)=bxdzst(im,j)
      bzdx_nova(jdm)=bzdxst(im,j)
      bzdz_nova(jdm)=bzdzst(im,j)
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
      psst(:,1)=psival_nova(1)
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
      th12_nova(jd)=thst(i)
      xx12_nova(jd)=xxst(i,j)
      zz12_nova(jd)=zzst(i,j)

      th34_nova(jd)=thst(i+ipi-3)
      xx34_nova(jd)=xxst(i+ipi-3,j)
      zz34_nova(jd)=zzst(i+ipi-3,j)
      enddo
      enddo


      do j=1,npsi
      do i=1,n2th+5
      byst(i,j)=xzero*g_nova(j)/xxst(i,j)
      bydxst(i,j)=-xzero*(gp_nova(j)*bzst(i,j)+g_nova(j)/xxst(i,j)**2)
      bydzst(i,j)=xzero*gp_nova(j)*bxst(i,j)
      pdxst(i,j)=-pp_nova(j)*bzst(i,j)*xxst(i,j)
      pdzst(i,j)=pp_nova(j)*bxst(i,j)*xxst(i,j)

      cxst(i,j)=-xzero*gp_nova(j)*bxst(i,j)
      czst(i,j)=-xzero*gp_nova(j)*bzst(i,j)
      cyst(i,j)=bxdzst(i,j)-bzdxst(i,j)
!      cyst(i,j)=xxst(i,j)*pp_nova(j)+xzero**2*g_nova(j)*gp_nova(j)/xxst(i,j)

      uyst(i,j)=omrot_nova(j)*xxst(i,j)
      uydxst(i,j)=-omprot_nova(j)*bzst(i,j)*xxst(i,j)**2+omrot_nova(j)
      uydzst(i,j)= omprot_nova(j)*bxst(i,j)*xxst(i,j)**2

      if(j.ge.2 .and. i.ge.3 .and. i.le.n2th+2 .and. i.ne.ipi) then
      if(i.lt.ipi) jd=(j-2)*(n2th-1)+i-2
      if(i.gt.ipi) jd=(j-2)*(n2th-1)+i-3
      by_nova(jd)=byst(i,j)
      bydx_nova(jd)=bydxst(i,j)
      bydz_nova(jd)=bydzst(i,j)
      pdx_nova(jd)=pdxst(i,j)
      pdz_nova(jd)=pdzst(i,j)
      cy_nova(jd)=cyst(i,j)
      cx_nova(jd)=cxst(i,j)
      cz_nova(jd)=czst(i,j)

      uy_nova(jd)=uyst(i,j)
      uydx_nova(jd)=uydxst(i,j)
      uydz_nova(jd)=uydzst(i,j)
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

      open(unit=101,file='fst_nova.dat',status='unknown',form='formatted')
      write(101,300)(((fffst(i,j,m),m=1,3),i=3,nthe+1),j=2,npsi)
 300  format(3(1x,e12.5))

      endif

       xmin=xplmin
       xmax=xplmax
       zmax=zplmax
       zmin=-zmax
       zmg=0.0
       psmin=minval(ps_nova)
       psmax=maxval(ps_nova)

       qmin=minval(q_nova)
       qmax=maxval(q_nova)
       q0=q_nova(1)
!      xmin=minval(xx_nova)
!      xmax=maxval(xx_nova)
!      zmin=minval(zz_nova)
!      zmax=maxval(zz_nova)
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
      write(890,5000) (xx_nova(j),zz_nova(j),th_nova(j),ps_nova(j),bx_nova(j),bxdx_nova(j),bxdz_nova(j),bz_nova(j),bzdx_nova(j),bzdz_nova(j),j=1,ndat)
 5000 format(10(1x,e17.9))
      close(890)
      open(891,file='stdat.dat')
      write(891,5100) ((thst(i),psst(i,j),xxst(i,j),zzst(i,j),tpst(i,j),tst(i,j),rst(i,j),i=1,n2th+5),j=1,npsi)
 5100 format(7(1x,e17.9))
      close(891)

      endif 
              
      return
      end

!ws*******************************************************************************
     subroutine read_nova_tm
     use declare 
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
      '   ip  ipnorm' )

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
      cip=curnorm*xzero*b0
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
      
      th_nova(jd)=thst(i)
      xx_nova(jd)=xxst(i,j)
      zz_nova(jd)=zzst(i,j)
      ps_nova(jd)=psst(i,j)            
      bx_nova(jd)=bxst(i,j)
      bz_nova(jd)=bzst(i,j)
      bxdx_nova(jd)=bxdxst(i,j)
      bxdz_nova(jd)=bxdzst(i,j)
      bzdx_nova(jd)=bzdxst(i,j)
      bzdz_nova(jd)=bzdzst(i,j)

      pt_nova(jd)=ptst(i,j)
      ptdx_nova(jd)=ptdxst(i,j)
      ptdz_nova(jd)=ptdzst(i,j)

      rh_nova(jd)=rhst(i,j)
      rhdx_nova(jd)=rhdxst(i,j)
      rhdz_nova(jd)=rhdzst(i,j)

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
      th_nova(jdm)=thst(im)
      xx_nova(jdm)=xxst(im,j)
      zz_nova(jdm)=zzst(im,j)
      ps_nova(jdm)=psst(im,j)
      bx_nova(jdm)=bxst(im,j)
      bz_nova(jdm)=bzst(im,j)
      bxdx_nova(jdm)=bxdxst(im,j)
      bxdz_nova(jdm)=bxdzst(im,j)
      bzdx_nova(jdm)=bzdxst(im,j)
      bzdz_nova(jdm)=bzdzst(im,j)

      pt_nova(jdm)=ptst(im,j)
      ptdx_nova(jdm)=ptdxst(im,j)
      ptdz_nova(jdm)=ptdzst(im,j)
      rh_nova(jdm)=rhst(im,j)
      rhdx_nova(jdm)=rhdxst(im,j)
      rhdz_nova(jdm)=rhdzst(im,j)
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
      read(889,2100) jj,psival_nova(j),q_nova(j),qp_nova(j),p_nova(j),pp_nova(j),g_nova(j),gp_nova(j),f_nova(j),fp_nova(j),fb_nova(j),fbp_nova(j)
      psival_nova(j)=psival_nova(j)/(b0*aa**2)
      p_nova(j)=p_nova(j)/(b0**2)
      g_nova(j)=g_nova(j)/b0
      qp_nova(j)=qp_nova(j)*(b0*aa**2)
      pp_nova(j)=pp_nova(j)/(b0**2)*(b0*aa**2)
      gp_nova(j)=gp_nova(j)/b0*(b0*aa**2)
   41 continue
      close(889) 
      omrot_nova(:)=0
      omprot_nova(:)=0   
 
      else
      open(889,file='q_p_g.dat')
      do 40 j=1,npsi
      read(889,2000) jj,psival_nova(j),q_nova(j),qp_nova(j),p_nova(j),pp_nova(j),g_nova(j),gp_nova(j),f_nova(j),fp_nova(j),fb_nova(j),fbp_nova(j),omrot_nova(j),omprot_nova(j)     
      psival_nova(j)=psival_nova(j)/(b0*aa**2)
      p_nova(j)=p_nova(j)/(b0**2)
      g_nova(j)=g_nova(j)/b0
      qp_nova(j)=qp_nova(j)*(b0*aa**2)
      pp_nova(j)=pp_nova(j)/(b0**2)*(b0*aa**2)
      gp_nova(j)=gp_nova(j)/b0*(b0*aa**2)
      omrot_nova(j)=omrot_nova(j)/(b0*aa)
      omprot_nova(j)=omprot_nova(j)/(b0*aa)*(b0*aa**2)    
     
   40 continue
      close(889)      
      endif

      xxst(:,1)=xmg
      zzst(:,1)=0
      psst(:,1)=psival_nova(1)
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
      th12_nova(jd)=thst(i)
      xx12_nova(jd)=xxst(i,j)
      zz12_nova(jd)=zzst(i,j)

      th34_nova(jd)=thst(i+ipi-3)
      xx34_nova(jd)=xxst(i+ipi-3,j)
      zz34_nova(jd)=zzst(i+ipi-3,j)
      enddo
      enddo


      do j=1,npsi
      do i=1,n2th+5
      byst(i,j)=xzero*g_nova(j)/xxst(i,j)
      bydxst(i,j)=-xzero*(gp_nova(j)*bzst(i,j)+g_nova(j)/xxst(i,j)**2)
      bydzst(i,j)=xzero*gp_nova(j)*bxst(i,j)

      cxst(i,j)=-xzero*gp_nova(j)*bxst(i,j)
      czst(i,j)=-xzero*gp_nova(j)*bzst(i,j)
      cyst(i,j)=bxdzst(i,j)-bzdxst(i,j)
!      cyst(i,j)=xxst(i,j)*pp_nova(j)+xzero**2*g_nova(j)*gp_nova(j)/xxst(i,j)

      uyst(i,j)=omrot_nova(j)*xxst(i,j)
      uydxst(i,j)=-omprot_nova(j)*bzst(i,j)*xxst(i,j)**2+omrot_nova(j)
      uydzst(i,j)= omprot_nova(j)*bxst(i,j)*xxst(i,j)**2

      if(j.ge.2 .and. i.ge.3 .and. i.le.n2th+2 .and. i.ne.ipi) then
      if(i.lt.ipi) jd=(j-2)*(n2th-1)+i-2
      if(i.gt.ipi) jd=(j-2)*(n2th-1)+i-3
      by_nova(jd)=byst(i,j)
      bydx_nova(jd)=bydxst(i,j)
      bydz_nova(jd)=bydzst(i,j)

      cy_nova(jd)=cyst(i,j)
      cx_nova(jd)=cxst(i,j)
      cz_nova(jd)=czst(i,j)

      uy_nova(jd)=uyst(i,j)
      uydx_nova(jd)=uydxst(i,j)
      uydz_nova(jd)=uydzst(i,j)
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

      open(unit=101,file='fst_nova.dat',status='unknown',form='formatted')
      write(101,300)(((fffst(i,j,m),m=1,3),i=3,nthe+1),j=2,npsi)
 300  format(3(1x,e12.5))

      endif

       psmin=minval(ps_nova)
       psmax=maxval(ps_nova)

       qmin=minval(q_nova)
       qmax=maxval(q_nova)
       q0=q_nova(1)
!      xmin=minval(xx_nova)
!      xmax=maxval(xx_nova)
!      zmin=minval(zz_nova)
!      zmax=maxval(zz_nova)
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
      write(890,5000) (xx_nova(j),zz_nova(j),th_nova(j),ps_nova(j),bx_nova(j),bxdx_nova(j),bxdz_nova(j),bz_nova(j),bzdx_nova(j),bzdz_nova(j),j=1,ndat)
 5000 format(10(1x,e17.9))
      close(890)
      open(891,file='stdat.dat')
      write(891,5100) ((thst(i),psst(i,j),xxst(i,j),zzst(i,j),tpst(i,j),tst(i,j),rst(i,j),i=1,n2th+5),j=1,npsi)
 5100 format(7(1x,e17.9))
      close(891)

      endif 
       
 1500 format(1x,i5,i5,15(1x,e17.9))
 2000 format(1x,i5,13(1x,e17.9))
 2100 format(1x,i5,11(1x,e17.9))
 4000 format(14(1x,e17.9))         
      return
      end

!ws************************************************************
     subroutine gridpnt
      use declare
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
	
!ws**************************************************************
     subroutine map_nova
     use declare
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

! interpolation, input ndat, xx_nova, zz_nova, f_nova, 
! output f(jx,jz), and its x & z directional derivatives at xxt and zzt grids.
     write(*,*) ndat
     call qshep2 ( ndat, xx_nova, zz_nova, ps_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
     write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_nova, zz_nova, ps_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, pst(jx,jz), pst_dx(jx,jz), pst_dz(jx,jz), ier )
 !     write(*,*) ier
!       pst(jx,mzt-jz+1)=pst(jx,jz)
!       pst_dx(jx,mzt-jz+1)=pst_dx(jx,jz)
!       pst_dz(jx,mzt-jz+1)=-pst_dz(jx,jz)
      tpt(jx,jz)=atan2(pst_dz(jx,jz),pst_dx(jx,jz))
      if(zzt(jz).lt.0) tpt(jx,jz)=tpt(jx,jz)+2*pi
      enddo
      enddo

     call qshep2 ( ndat12, xx12_nova, zz12_nova, th12_nova, nq, nw, nr, lcell, lnext12, ximin, zimin, &
        dxi, dzi, rimax, risq12, aw12, ier )
!     write(*,*) ier
      do jz=mzt/2+1,mzt
      do jx=1,mxt
      call qs2grd ( xxt(jx), zzt(jz), ndat12, xx12_nova, zz12_nova, th12_nova, nr, lcell, lnext12, ximin, &
        zimin, dxi, dzi, rimax, risq12, aw12, tht(jx,jz), tht_dx(jx,jz), tht_dz(jx,jz), ier )
      enddo
      enddo  
      
      call qshep2 ( ndat34, xx34_nova, zz34_nova, th34_nova, nq, nw, nr, lcell, lnext34, ximin, zimin, &
        dxi, dzi, rimax, risq34, aw34, ier )
!     write(*,*) ier
      do jz=1,mzt/2
      do jx=1,mxt
      call qs2grd ( xxt(jx), zzt(jz), ndat34, xx34_nova, zz34_nova, th34_nova, nr, lcell, lnext34, ximin, &
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
      call qshep2 ( ndat, xx_nova, zz_nova, bx_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_nova, zz_nova, bx_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, bx(jx,jz), bx_dx(jx,jz), bx_dz(jx,jz), ier )
      enddo
      enddo

!!ws:bz
      call qshep2 ( ndat, xx_nova, zz_nova, bz_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_nova, zz_nova, bz_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, bz(jx,jz), bz_dx(jx,jz), bz_dz(jx,jz), ier )
      enddo
      enddo

!!ws:bxdx
      call qshep2 ( ndat, xx_nova, zz_nova, bxdx_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_nova, zz_nova, bxdx_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, bxdx(jx,jz), qdxout,qdzout, ier )
      enddo
      enddo

!!ws:bxdz
      call qshep2 ( ndat, xx_nova, zz_nova, bxdz_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_nova, zz_nova, bxdz_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, bxdz(jx,jz), qdxout,qdzout, ier )
 !     write(*,*) ier
      enddo
      enddo

!!ws:bzdx
      call qshep2 ( ndat, xx_nova, zz_nova, bzdx_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_nova, zz_nova, bzdx_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, bzdx(jx,jz), qdxout,qdzout, ier )
 !     write(*,*) ier
      enddo
      enddo

!!ws:bzdz
      call qshep2 ( ndat, xx_nova, zz_nova, bzdz_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_nova, zz_nova, bzdz_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, bzdz(jx,jz), qdxout,qdzout, ier )
 !     write(*,*) ier
      enddo
      enddo

!!ws:by
      call qshep2 ( ndat, xx_nova, zz_nova, by_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_nova, zz_nova, by_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, by(jx,jz), by_dx(jx,jz), by_dz(jx,jz), ier )
 !     write(*,*) ier
      enddo
      enddo
!!ws:bydx
      call qshep2 ( ndat, xx_nova, zz_nova, bydx_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_nova, zz_nova, bydx_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, bydx(jx,jz), qdxout,qdzout, ier )
 !     write(*,*) ier
      enddo
      enddo
!!ws:bydz
      call qshep2 ( ndat, xx_nova, zz_nova, bydz_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_nova, zz_nova, bydz_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, bydz(jx,jz), qdxout,qdzout, ier )
 !     write(*,*) ier
      enddo
      enddo

!!ws:pdx
      call qshep2 ( ndat, xx_nova, zz_nova, pdx_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_nova, zz_nova, pdx_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, pdx(jx,jz), qdxout,qdzout, ier )
 !     write(*,*) ier
      enddo
      enddo
!!ws:pdz
      call qshep2 ( ndat, xx_nova, zz_nova, pdz_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_nova, zz_nova, pdz_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, pdz(jx,jz), qdxout,qdzout, ier )
 !     write(*,*) ier
      enddo
      enddo
!!ws:cy
      call qshep2 ( ndat, xx_nova, zz_nova, cy_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_nova, zz_nova, cy_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, cy(jx,jz), cy_dx(jx,jz), cy_dz(jx,jz), ier )
 !     write(*,*) ier
      enddo
      enddo
!!ws:cx
      call qshep2 ( ndat, xx_nova, zz_nova, cx_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_nova, zz_nova, cx_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, cx(jx,jz), cx_dx(jx,jz), cx_dz(jx,jz), ier )
 !     write(*,*) ier
      enddo
      enddo
!!ws:cz
      call qshep2 ( ndat, xx_nova, zz_nova, cz_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_nova, zz_nova, cz_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, cz(jx,jz), cz_dx(jx,jz), cz_dz(jx,jz), ier )
 !     write(*,*) ier
      enddo
      enddo
  
 !!ws:uy
      call qshep2 ( ndat, xx_nova, zz_nova, uy_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_nova, zz_nova, uy_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, uy(jx,jz), uy_dx(jx,jz), uy_dz(jx,jz), ier )
 !     write(*,*) ier
      enddo
      enddo
!!ws:uydx
      call qshep2 ( ndat, xx_nova, zz_nova, uydx_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_nova, zz_nova, uydx_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, uydx(jx,jz), qdxout,qdzout, ier )
 !     write(*,*) ier
      enddo
      enddo
!!ws:uydz
      call qshep2 ( ndat, xx_nova, zz_nova, uydz_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_nova, zz_nova, uydz_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, uydz(jx,jz), qdxout,qdzout, ier )
 !     write(*,*) ier
      enddo
      enddo
 !!ws:pt
      call qshep2 ( ndat, xx_nova, zz_nova, pt_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_nova, zz_nova, pt_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, pt(jx,jz), pt_dx(jx,jz), pt_dz(jx,jz), ier )
 !     write(*,*) ier
      enddo
      enddo
!!ws:ptdx
      call qshep2 ( ndat, xx_nova, zz_nova, ptdx_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_nova, zz_nova, ptdx_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, ptdx(jx,jz), qdxout,qdzout, ier )
 !     write(*,*) ier
      enddo
      enddo
!!ws:ptdz
      call qshep2 ( ndat, xx_nova, zz_nova, ptdz_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_nova, zz_nova, ptdz_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, ptdz(jx,jz), qdxout,qdzout, ier )
 !     write(*,*) ier
      enddo
      enddo 
!!ws:rh
      call qshep2 ( ndat, xx_nova, zz_nova, rh_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_nova, zz_nova, rh_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, rh(jx,jz), rh_dx(jx,jz),rh_dz(jx,jz), ier )
 !     write(*,*) ier
      enddo
      enddo
!!ws:rhdx
      call qshep2 ( ndat, xx_nova, zz_nova, rhdx_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_nova, zz_nova, rhdx_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, rhdx(jx,jz), qdxout,qdzout, ier )
 !     write(*,*) ier
      enddo
      enddo
!!ws:rhdz
      call qshep2 ( ndat, xx_nova, zz_nova, rhdz_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_nova, zz_nova, rhdz_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, rhdz(jx,jz), qdxout,qdzout, ier )
 !     write(*,*) ier
      enddo
      enddo 
      
          
      if(nrank.eq.0) then
      open(unit=99,file='mapall',status='unknown',form='unformatted')
      write(99) bx,bxdx,bxdz,bz,bzdx,bzdz,by,bydx,bydz,cx,cy,cz,pt,ptdx,ptdz,rh,rhdx,rhdz,cx_dx, &
                cx_dz,cz_dx,cz_dz,cy_dx,cy_dz,uy,uydx,uydz
      close(99)
      endif

      else

      open(unit=99,file='mapall',status='unknown',form='unformatted')
      read(99)  bx,bxdx,bxdz,bz,bzdx,bzdz,by,bydx,bydz,cx,cy,cz,pt,ptdx,ptdz,rh,rhdx,rhdz,cx_dx, &
                cx_dz,cz_dx,cz_dz,cy_dx,cy_dz,uy,uydx,uydz
      close(99)
      endif

      psia=psival_nova(mpsa)
      weight=1.
      psiam=weight*psia+(1-weight)*psival_nova(mpsa-1)

      do jx=1,mxt
      do jz=1,mzt
      bpol(jx,jz)=sqrt(max(0.d0,bx(jx,jz)**2+bz(jx,jz)**2))
      rr2t(jx,jz)=(pst(jx,jz)-psmin)/(psia-psmin)
      rrt(jx,jz)=sqrt(max(0.d0,rr2t(jx,jz)))
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
    if(pst(jx,mzt/2).lt.psival_nova(js) .and. pst(jx-1,mzt/2).ge.psival_nova(js)) jxsmin(js)=jx
    if(pst(jx,mzt/2).lt.psival_nova(js) .and. pst(jx+1,mzt/2).ge.psival_nova(js)) jxsmax(js)=jx
    enddo
    do jz=1+1,mzt-1 
    if(pst(mxt/2+1,jz).lt.psival_nova(js) .and. pst(mxt/2+1,jz-1).ge.psival_nova(js)) jzsmin(js)=jz
    if(pst(mxt/2+1,jz).lt.psival_nova(js) .and. pst(mxt/2+1,jz+1).ge.psival_nova(js)) jzsmax(js)=jz
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
      if( pst(jx,jz).lt. psival_nova(npsi-1)) then
      j=1
      do while(psival_nova(j) .lt. pst(jx,jz))
      j=j+1
      enddo
      if(j==1) j=j+2
      if(j==2) j=j+1
      call interp1d3l(q_nova(j-2),q_nova(j-1),q_nova(j),q_nova(j+1), &
                    psival_nova(j-2),psival_nova(j-1),psival_nova(j),psival_nova(j+1),pst(jx,jz),qsf(jx,jz))
      call interp1d3l(g_nova(j-2),g_nova(j-1),g_nova(j),g_nova(j+1), &
                    psival_nova(j-2),psival_nova(j-1),psival_nova(j),psival_nova(j+1),pst(jx,jz),g(jx,jz))
      call interp1d3l(p_nova(j-2),p_nova(j-1),p_nova(j),p_nova(j+1), &
                    psival_nova(j-2),psival_nova(j-1),psival_nova(j),psival_nova(j+1),pst(jx,jz),p(jx,jz))
      call interp1d3l(gp_nova(j-2),gp_nova(j-1),gp_nova(j),gp_nova(j+1), &
                    psival_nova(j-2),psival_nova(j-1),psival_nova(j),psival_nova(j+1),pst(jx,jz),gp(jx,jz))
      call interp1d3l(pp_nova(j-2),pp_nova(j-1),pp_nova(j),pp_nova(j+1), &
                    psival_nova(j-2),psival_nova(j-1),psival_nova(j),psival_nova(j+1),pst(jx,jz),pp(jx,jz))

     
      call interp1d3l(omrot_nova(j-2),omrot_nova(j-1),omrot_nova(j),omrot_nova(j+1), &
                    psival_nova(j-2),psival_nova(j-1),psival_nova(j),psival_nova(j+1),pst(jx,jz),omrot(jx,jz))
      call interp1d3l(omprot_nova(j-2),omprot_nova(j-1),omprot_nova(j),omprot_nova(j+1), &
                    psival_nova(j-2),psival_nova(j-1),psival_nova(j),psival_nova(j+1),pst(jx,jz),omprot(jx,jz))

      endif

      if( pst(jx,jz).ge. psival_nova(npsi-1) .and. pst(jx,jz).lt.psival_nova(npsi)) then
      call interp1d3l(q_nova(npsi-3),q_nova(npsi-2),q_nova(npsi-1),q_nova(npsi), &
                    psival_nova(npsi-3),psival_nova(npsi-2),psival_nova(npsi-1),psival_nova(npsi),pst(jx,jz),qsf(jx,jz))
      call interp1d3l(g_nova(npsi-3),g_nova(npsi-2),g_nova(npsi-1),g_nova(npsi), &
                    psival_nova(npsi-3),psival_nova(npsi-2),psival_nova(npsi-1),psival_nova(npsi),pst(jx,jz),g(jx,jz))
      call interp1d3l(p_nova(npsi-3),p_nova(npsi-2),p_nova(npsi-1),p_nova(npsi), &
                    psival_nova(npsi-3),psival_nova(npsi-2),psival_nova(npsi-1),psival_nova(npsi),pst(jx,jz),p(jx,jz))
      call interp1d3l(gp_nova(npsi-3),gp_nova(npsi-2),gp_nova(npsi-1),gp_nova(npsi), &
                    psival_nova(npsi-3),psival_nova(npsi-2),psival_nova(npsi-1),psival_nova(npsi),pst(jx,jz),gp(jx,jz))
      call interp1d3l(pp_nova(npsi-3),pp_nova(npsi-2),pp_nova(npsi-1),pp_nova(npsi), &
                    psival_nova(npsi-3),psival_nova(npsi-2),psival_nova(npsi-1),psival_nova(npsi),pst(jx,jz),pp(jx,jz))

      call interp1d3l(omrot_nova(npsi-3),omrot_nova(npsi-2),omrot_nova(npsi-1),omrot_nova(npsi), &
                    psival_nova(npsi-3),psival_nova(npsi-2),psival_nova(npsi-1),psival_nova(npsi),pst(jx,jz),omrot(jx,jz))
      call interp1d3l(omprot_nova(npsi-3),omprot_nova(npsi-2),omprot_nova(npsi-1),omprot_nova(npsi), &
                    psival_nova(npsi-3),psival_nova(npsi-2),psival_nova(npsi-1),psival_nova(npsi),pst(jx,jz),omprot(jx,jz))
      endif
            
!      uy(jx,jz)=xxt(jx)*omrot(jx,jz)
!      uydx(jx,jz)=-omprot(jx,jz)*bz(jx,jz)*xxt(jx)**2+omrot(jx,jz)
!      uydz(jx,jz)= omprot(jx,jz)*bx(jx,jz)*xxt(jx)**2
      enddo
      enddo


      
      if(rshear) then
      j=1
      do while(q_nova(j) .ge. qmode)
      j=j+1
      enddo
      if(j==1) j=j+2
      if(j==2) j=j+1
      call interp1d3l(psival_nova(j-2),psival_nova(j-1),psival_nova(j),psival_nova(j+1), &
                    q_nova(j-2),q_nova(j-1),q_nova(j),q_nova(j+1),qmode,ps1mode)
      call interp1d3l(xxst(3,j-2),xxst(3,j-1),xxst(3,j),xxst(3,j+1), &
                    q_nova(j-2),q_nova(j-1),q_nova(j),q_nova(j+1),qmode,xx1mode)
      else
      j=1
      endif

      do while(q_nova(j) .lt. qmode)
      j=j+1
      enddo
      if(j==1) j=j+2
      if(j==2) j=j+1
      call interp1d3l(psival_nova(j-2),psival_nova(j-1),psival_nova(j),psival_nova(j+1), &
                    q_nova(j-2),q_nova(j-1),q_nova(j),q_nova(j+1),qmode,psmode)

      call interp1d3l(xxst(3,j-2),xxst(3,j-1),xxst(3,j),xxst(3,j+1), &
                    q_nova(j-2),q_nova(j-1),q_nova(j),q_nova(j+1),qmode,xxmode)
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
      write(205,100)((pst(jx,jz),pst_dx(jx,jz),pst_dz(jx,jz),qsf(jx,jz),p(jx,jz),pp(jx,jz),g(jx,jz),gp(jx,jz),jx=1,mxt),jz=1,mzt)
 100  format(8(1x,e12.5))
      close(205)

      open(unit=204,file='gridts.dat',status='unknown',form='formatted')
      write(204,200)( (tpt(jx,jz),jx=1,mxt),jz=1,mzt)
      close(204)
      open(unit=206,file='gridth.dat',status='unknown',form='formatted')
      write(206,200)( (tht(jx,jz),jx=1,mxt),jz=1,mzt)
      close(206)
      open(unit=207,file='gridtc.dat',status='unknown',form='formatted')
      write(207,3000)( (tht(jx,jz),tht_dx(jx,jz),tht_dz(jx,jz),jx=1,mxt),jz=1,mzt)
      close(207)
  
200   format((1x,e12.5))
3000  format(3(1x,e12.5))      
      open(unit=301,file='bx0.dat',status='unknown',form='formatted')
      write(301,500)((bx(jx,jz),bx_dx(jx,jz),bxdx(jx,jz),bx_dz(jx,jz),bxdz(jx,jz),jx=1,mxt),jz=1,mzt)
 500  format(5(1x,e12.5))
      close(301)
      open(unit=302,file='bz0.dat',status='unknown',form='formatted')
      write(302,500)((bz(jx,jz),bz_dx(jx,jz),bzdx(jx,jz),bz_dz(jx,jz),bzdz(jx,jz),jx=1,mxt),jz=1,mzt)
      close(302)
      open(unit=303,file='c0.dat',status='unknown',form='formatted')
      write(303,900)((cx(jx,jz),cx_dx(jx,jz),cx_dz(jx,jz),cy(jx,jz),cy_dx(jx,jz),cy_dz(jx,jz),cz(jx,jz),cz_dx(jx,jz),cz_dz(jx,jz),jx=1,mxt),jz=1,mzt)
 900  format(9(1x,e12.5))
      close(303)
      open(unit=304,file='pt0.dat',status='unknown',form='formatted')
      write(304,500)((pt(jx,jz),pt_dx(jx,jz),ptdx(jx,jz),pt_dz(jx,jz),ptdz(jx,jz),jx=1,mxt),jz=1,mzt)
      close(304)
      open(unit=305,file='rh0.dat',status='unknown',form='formatted')
      write(305,500)((rh(jx,jz),rh_dx(jx,jz),rhdx(jx,jz),rh_dz(jx,jz),rhdz(jx,jz),jx=1,mxt),jz=1,mzt)
      close(305)

      open(unit=306,file='uy0.dat',status='unknown',form='formatted')
      write(306,500)((uy(jx,jz),uy_dx(jx,jz),uydx(jx,jz),uy_dz(jx,jz),uydz(jx,jz),jx=1,mxt),jz=1,mzt)
      close(306)
     endif
!ws150303
     call estimate_pst1
     call readmap_wch
!ws150303
     call last_grid
     return
     end

!ws****************************************************************************
     subroutine mapeq_st2xz(fst,fxz,jjx,jjz,itp,isp,rs)
     use declare
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
      use declare
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
     call interp1d3l(psival_nova(isp-2),psival_nova(isp-1),psival_nova(isp),psival_nova(isp+1), &
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

      psia=psival_nova(mpsa)
      weight=1.
      psiam=weight*psia+(1-weight)*psival_nova(mpsa-1)
      

    do js=mpsa,mps4,-1
    do jx=1+1,mxt-1 
    if(pst(jx,mzt/2).lt.psival_nova(js) .and. pst(jx-1,mzt/2).ge.psival_nova(js)) jxsmin(js)=jx
    if(pst(jx,mzt/2).lt.psival_nova(js) .and. pst(jx+1,mzt/2).ge.psival_nova(js)) jxsmax(js)=jx
    enddo
    do jz=1+1,mzt-1 
    if(pst(mxt/2+1,jz).lt.psival_nova(js) .and. pst(mxt/2+1,jz-1).ge.psival_nova(js)) jzsmin(js)=jz
    if(pst(mxt/2+1,jz).lt.psival_nova(js) .and. pst(mxt/2+1,jz+1).ge.psival_nova(js)) jzsmax(js)=jz
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
      if( pst(jx,jz).lt. psival_nova(npsi-1)) then
      j=1
      do while(psival_nova(j) .lt. pst(jx,jz))
      j=j+1
      enddo
      if(j==1) j=j+2
      if(j==2) j=j+1
      call interp1d3l(q_nova(j-2),q_nova(j-1),q_nova(j),q_nova(j+1), &
                    psival_nova(j-2),psival_nova(j-1),psival_nova(j),psival_nova(j+1),pst(jx,jz),qsf(jx,jz))
      call interp1d3l(g_nova(j-2),g_nova(j-1),g_nova(j),g_nova(j+1), &
                    psival_nova(j-2),psival_nova(j-1),psival_nova(j),psival_nova(j+1),pst(jx,jz),g(jx,jz))
      call interp1d3l(p_nova(j-2),p_nova(j-1),p_nova(j),p_nova(j+1), &
                    psival_nova(j-2),psival_nova(j-1),psival_nova(j),psival_nova(j+1),pst(jx,jz),p(jx,jz))
      call interp1d3l(gp_nova(j-2),gp_nova(j-1),gp_nova(j),gp_nova(j+1), &
                    psival_nova(j-2),psival_nova(j-1),psival_nova(j),psival_nova(j+1),pst(jx,jz),gp(jx,jz))
      call interp1d3l(pp_nova(j-2),pp_nova(j-1),pp_nova(j),pp_nova(j+1), &
                    psival_nova(j-2),psival_nova(j-1),psival_nova(j),psival_nova(j+1),pst(jx,jz),pp(jx,jz))

     
      call interp1d3l(omrot_nova(j-2),omrot_nova(j-1),omrot_nova(j),omrot_nova(j+1), &
                    psival_nova(j-2),psival_nova(j-1),psival_nova(j),psival_nova(j+1),pst(jx,jz),omrot(jx,jz))
      call interp1d3l(omprot_nova(j-2),omprot_nova(j-1),omprot_nova(j),omprot_nova(j+1), &
                    psival_nova(j-2),psival_nova(j-1),psival_nova(j),psival_nova(j+1),pst(jx,jz),omprot(jx,jz))

      endif

      if( pst(jx,jz).ge. psival_nova(npsi-1) .and. pst(jx,jz).lt.psival_nova(npsi)) then
      call interp1d3l(q_nova(npsi-3),q_nova(npsi-2),q_nova(npsi-1),q_nova(npsi), &
                    psival_nova(npsi-3),psival_nova(npsi-2),psival_nova(npsi-1),psival_nova(npsi),pst(jx,jz),qsf(jx,jz))
      call interp1d3l(g_nova(npsi-3),g_nova(npsi-2),g_nova(npsi-1),g_nova(npsi), &
                    psival_nova(npsi-3),psival_nova(npsi-2),psival_nova(npsi-1),psival_nova(npsi),pst(jx,jz),g(jx,jz))
      call interp1d3l(p_nova(npsi-3),p_nova(npsi-2),p_nova(npsi-1),p_nova(npsi), &
                    psival_nova(npsi-3),psival_nova(npsi-2),psival_nova(npsi-1),psival_nova(npsi),pst(jx,jz),p(jx,jz))
      call interp1d3l(gp_nova(npsi-3),gp_nova(npsi-2),gp_nova(npsi-1),gp_nova(npsi), &
                    psival_nova(npsi-3),psival_nova(npsi-2),psival_nova(npsi-1),psival_nova(npsi),pst(jx,jz),gp(jx,jz))
      call interp1d3l(pp_nova(npsi-3),pp_nova(npsi-2),pp_nova(npsi-1),pp_nova(npsi), &
                    psival_nova(npsi-3),psival_nova(npsi-2),psival_nova(npsi-1),psival_nova(npsi),pst(jx,jz),pp(jx,jz))

      call interp1d3l(omrot_nova(npsi-3),omrot_nova(npsi-2),omrot_nova(npsi-1),omrot_nova(npsi), &
                    psival_nova(npsi-3),psival_nova(npsi-2),psival_nova(npsi-1),psival_nova(npsi),pst(jx,jz),omrot(jx,jz))
      call interp1d3l(omprot_nova(npsi-3),omprot_nova(npsi-2),omprot_nova(npsi-1),omprot_nova(npsi), &
                    psival_nova(npsi-3),psival_nova(npsi-2),psival_nova(npsi-1),psival_nova(npsi),pst(jx,jz),omprot(jx,jz))
      endif
            
!      uy(jx,jz)=xxt(jx)*omrot(jx,jz)
!      uydx(jx,jz)=-omprot(jx,jz)*bz(jx,jz)*xxt(jx)**2+omrot(jx,jz)
!      uydz(jx,jz)= omprot(jx,jz)*bx(jx,jz)*xxt(jx)**2
      enddo
      enddo


      
      if(rshear) then
      j=1
      do while(q_nova(j) .ge. qmode)
      j=j+1
      enddo
      if(j==1) j=j+2
      if(j==2) j=j+1
      call interp1d3l(psival_nova(j-2),psival_nova(j-1),psival_nova(j),psival_nova(j+1), &
                    q_nova(j-2),q_nova(j-1),q_nova(j),q_nova(j+1),qmode,ps1mode)
      call interp1d3l(xxst(3,j-2),xxst(3,j-1),xxst(3,j),xxst(3,j+1), &
                    q_nova(j-2),q_nova(j-1),q_nova(j),q_nova(j+1),qmode,xx1mode)
      else
      j=1
      endif

      do while(q_nova(j) .lt. qmode)
      j=j+1
      enddo
      if(j==1) j=j+2
      if(j==2) j=j+1
      call interp1d3l(psival_nova(j-2),psival_nova(j-1),psival_nova(j),psival_nova(j+1), &
                    q_nova(j-2),q_nova(j-1),q_nova(j),q_nova(j+1),qmode,psmode)

      call interp1d3l(xxst(3,j-2),xxst(3,j-1),xxst(3,j),xxst(3,j+1), &
                    q_nova(j-2),q_nova(j-1),q_nova(j),q_nova(j+1),qmode,xxmode)
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
      write(207,2000)((is_inp(jx,jz),it_inp(jx,jz),jx=1,mxt),jz=1,mzt)
2000  format(2(1x,i5))
      close(207)

      open(unit=205,file='pst_qpg.dat',status='unknown',form='formatted')
      write(205,800)((pst(jx,jz),txzt(jx,jz),rxzt(jx,jz),qsf(jx,jz),p(jx,jz),pp(jx,jz),g(jx,jz),gp(jx,jz),jx=1,mxt),jz=1,mzt)
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
      write(301,500)((bx(jx,jz),bxdx(jx,jz),bxdz(jx,jz),jx=1,mxt),jz=1,mzt)
 500  format(3(1x,e12.5))
      close(301)
      open(unit=302,file='bz0.dat',status='unknown',form='formatted')
      write(302,500)((bz(jx,jz),bzdx(jx,jz),bzdz(jx,jz),jx=1,mxt),jz=1,mzt)
      close(302)
      open(unit=303,file='c0.dat',status='unknown',form='formatted')
      write(303,500)((cx(jx,jz),cy(jx,jz),cz(jx,jz),jx=1,mxt),jz=1,mzt)
      close(303)
     endif

     call last_grid
     return
     end

!ws*******************************************************************************
     subroutine readmap_wch
     use declare 
      integer, parameter :: nq = 33
      integer, parameter :: nr = 85 !int( sqrt(ndat/3) )
      integer, parameter :: nw = 39
    !
      integer lcell(nr,nr)
      integer lnext(ndat),lnext12(ndat12),lnext34(ndat34)
      real*8 risq(ndat),risq12(ndat12),risq34(ndat34)
      real*8 tch_nova(ndat),tch12_nova(ndat12),tch34_nova(ndat34),tcdx_nova(ndat),tcdz_nova(ndat),aj_nova(ndat)
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
      tch_nova(jd)=tcst(i,j)
      tcdx_nova(jd)=tcdxst(i,j)
      tcdz_nova(jd)=tcdzst(i,j)
      aj_nova(jd)=ajst(i,j)

      if(i.gt.3) then
      im=2*nthe+2-(i-2)
      tcst(im,j)=2*pi-tcst(i,j)
      tcdxst(im,j)=-tcdxst(i,j)
      tcdzst(im,j)=tcdzst(i,j)
      ajst(im,j)=ajst(i,j)

      jdm=(j-2)*(n2th-1)+im-3
      tch_nova(jdm)=tcst(im,j)
      tcdx_nova(jdm)=tcdxst(im,j)
      tcdz_nova(jdm)=tcdzst(im,j)
      aj_nova(jdm)=ajst(im,j)

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
      tch12_nova(jd)=tcst(i,j)
      tch34_nova(jd)=tcst(i+ipi-3,j)
      enddo
      enddo

      call qshep2 ( ndat12, xx12_nova, zz12_nova, tch12_nova, nq, nw, nr, lcell, lnext12, ximin, zimin, &
        dxi, dzi, rimax, risq12, aw12, ier )
!     write(*,*) ier
      do jz=mzt/2+1,mzt
      do jx=1,mxt
      call qs2grd ( xxt(jx), zzt(jz), ndat12, xx12_nova, zz12_nova, tch12_nova, nr, lcell, lnext12, ximin, &
        zimin, dxi, dzi, rimax, risq12, aw12, tcht(jx,jz), tcht_dx(jx,jz), tcht_dz(jx,jz), ier )
      enddo
      enddo  
!      
      call qshep2 ( ndat34, xx34_nova, zz34_nova, tch34_nova, nq, nw, nr, lcell, lnext34, ximin, zimin, &
        dxi, dzi, rimax, risq34, aw34, ier )
!     write(*,*) ier
      do jz=1,mzt/2
      do jx=1,mxt
      call qs2grd ( xxt(jx), zzt(jz), ndat34, xx34_nova, zz34_nova, tch34_nova, nr, lcell, lnext34, ximin, &
        zimin, dxi, dzi, rimax, risq34, aw34, tcht(jx,jz), tcht_dx(jx,jz), tcht_dz(jx,jz), ier )
      enddo
      enddo

      call qshep2 ( ndat, xx_nova, zz_nova, tcdx_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_nova, zz_nova, tcdx_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, tchdx(jx,jz), wwdx, wwdz, ier )
      enddo
      enddo
      call qshep2 ( ndat, xx_nova, zz_nova, tcdz_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_nova, zz_nova, tcdz_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, tchdz(jx,jz), wwdx, wwdz, ier )
      enddo
      enddo

      call qshep2 ( ndat, xx_nova, zz_nova, aj_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      write(*,*) ier
      do jx=1,mxt
      do jz=1,mzt
      call qs2grd ( xxt(jx), zzt(jz), ndat, xx_nova, zz_nova, aj_nova, nr, lcell, lnext, ximin, &
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
      write(791,2100) ((tcst(i,j),tcdxst(i,j),tcdzst(i,j),ajst(i,j),i=1,n2th+5),j=1,npsi)
 2100 format(4(1x,e17.9))
      close(791)

      open(unit=792,file='wtcxz.dat',status='unknown',form='formatted')
      write(792,300)((tcht(jx,jz),tchdx(jx,jz),tchdz(jx,jz),jx=1,mxt),jz=1,mzt)
 300  format(3(1x,e12.5))
      close(792)

      open(unit=793,file='ajxz.dat',status='unknown',form='formatted')
      write(793,300)((ajt(jx,jz),ajt_dx(jx,jz),ajt_dz(jx,jz),jx=1,mxt),jz=1,mzt)
      close(793)
      endif 
       
       
      return
      end

!ws********************************************
     subroutine estimate_pst1
      use declare
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
      write(905,100)((pst(jx,jz),pst1(jx,jz),jx=1,mxt),jz=1,mzt)
 100  format(2(1x,e12.5))
      close(905)
      endif
      
      return
      end

!ws*************************************************************************************          
     subroutine last_grid
     use declare
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
         do while(ps(js).lt. psival_nova(j))
         j=j-1
         enddo
         if(j==npsi-1) j=j-1
         ip_s(js)=j
         do i=1,n2th+5
         call interp1d3l(xxst(i,j+2),xxst(i,j+1),xxst(i,j),xxst(i,j-1), &
                        psival_nova(j+2),psival_nova(j+1),psival_nova(j),psival_nova(j-1),ps(js),xxs(i,js))
         call interp1d3l(zzst(i,j+2),zzst(i,j+1),zzst(i,j),zzst(i,j-1), &
                        psival_nova(j+2),psival_nova(j+1),psival_nova(j),psival_nova(j-1),ps(js),zzs(i,js))
         call interp1d3l(tpst(i,j+2),tpst(i,j+1),tpst(i,j),tpst(i,j-1), &
                        psival_nova(j+2),psival_nova(j+1),psival_nova(j),psival_nova(j-1),ps(js),tps(i,js))
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
!     ps(js)=psival_nova(js)
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
     write(*,*) js,ps(js),psival_nova(js)
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
	return
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

!hw**************************************************************
     subroutine read_pgfile
	use declare
	implicit none
	integer i, j
	real*8 weight
	real*8 PSIMAG, PSIBRY, RRMIN, RRMAX, ZZMIN, ZZMAX, RMAXIS, ZMAXIS
	real*8, dimension(npsi_pgfile) :: PSI_1D, Q_1D, NE_1D, NI_1D, PRES_1D, PE_1D, PI_1D, TE_1D, TI_1D
	real*8, dimension(nlim) :: RLIM, ZLIM
	real*8, dimension(nbs) :: RBS, ZBS
	real*8, dimension(mxt,mzt) :: RR_2D, ZZ_2D, THT_2D, PSI_2D, QRZ_2D, NE_2D, NI_2D, &
		  PRES_2D, PE_2D, PI_2D, VR_2D, VT_2D, VZ_2D, &
		  BR_2D, BT_2D, BZ_2D, JR_2D, JT_2D, JZ_2D, &
		  NE_DRRZ, NE_DZRZ, NI_DRRZ, NI_DZRZ, P_DRRZ, P_DZRZ, &
		  PE_DRRZ, PE_DZRZ, PI_DRRZ, PI_DZRZ, VR_DRRZ, VR_DZRZ, &
		  VT_DRRZ, VT_DZRZ, VZ_DRRZ, VZ_DZRZ, BR_DRRZ, BR_DZRZ, &
		  BT_DRRZ, BT_DZRZ, BZ_DRRZ, BZ_DZRZ, JR_DRRZ, JR_DZRZ, &
		  JT_DRRZ, JT_DZRZ, JZ_DRRZ, JZ_DZRZ, PSI_DRRZ, PSI_DZRZ, &
		  ERR_R_2D, ERR_T_2D, ERR_Z_2D, DIVB_2D
     real*8, dimension(mbm) :: psb2
     integer, dimension(npr) :: itb1
     integer ib2, mb2
     real*8 psb2min, psb2max, dps
     character*10 output
     character*3 cn1

	include 'mpif.h'
! read eq_pgfile_1d.dat, functions of psi
	open(1,file='eq_pgfile_1d.dat')
	read(1,*)
!    2 format(1h1,1x,'PSIMAG, PSIBRY, RMIN, RMAX, ZMIN, ZMAX, RMAXIS, ZMAXIS')
    	read(1,3) PSIMAG, PSIBRY, RRMIN, RRMAX, ZZMIN, ZZMAX, RMAXIS, ZMAXIS
    3 format(8(1x,e17.9))         
    	read(1,*)
!    4 format(1h1,1x,'PSI_1D, Q_1D, NE_1D, NI_1D, PRES_1D, PE_1D, PI_1D, TE_1D, TI_1D'	
	do i=1,npsi_pgfile
	read(1,5) PSI_1D(i), Q_1D(i), NE_1D(i), NI_1D(i), PRES_1D(i), PE_1D(i), PI_1D(i), TE_1D(i), TI_1D(i)
	enddo
    5 format(9(1x,e17.9))
   	close(1)

! read eq_pgfile_rr.dat, grids for xxt, i.e., rr
	open(6,file='eq_pgfile_rr.dat')
	do i=1,mxt
	read(6,7) xxt(i)
	enddo
    7 format(1(1x,e17.9))
    	close(6)

! read eq_pgfile_zz.dat, grids for zzt, i.e., zz
	open(8,file='eq_pgfile_zz.dat')
	do i=1,mzt
	read(8,9) zzt(i)
	enddo
    9 format(1(1x,e17.9))
    	close(8)


! read eq_pgfile_rzlim.dat, limter grids	
	open(10,file='eq_pgfile_rzlim.dat')
	do i=1,nlim
	read(10,11) RLIM(i), ZLIM(i)
	enddo
   11 format(2(1x,e17.9))
   	close(10)

! read eq_pgfile_rzbs.dat, closed flux surface grids
	open(12,file='eq_pgfile_rzbs.dat')
	do i=1,nbs
	read(12,13) RBS(i), ZBS(i)
	enddo
   13 format(2(1x,e17.9))
   	close(12)

! read eq_pgfile_2d.dat, data in mxt*mzt grids for 
! RR, ZZ, PSI, Q, NE, NI, PRES, PE, PI, VR, VT, VZ, BR, BT, BZ, JR, JT, JZ, 
! [d/dR and d/dZ of above variables(NE-JZ)], ERR_R, ERR_T, ERR_Z, DIVB
	open(14,file='eq_pgfile_2d.dat')
	read(14,*)
	do j=1,mzt
	do i=1,mxt
	read(14,15) RR_2D(i,j), ZZ_2D(i,j),THT_2D(i,j), PSI_2D(i,j), QRZ_2D(i,j), NE_2D(i,j), NI_2D(i,j), &
		  PRES_2D(i,j), PE_2D(i,j), PI_2D(i,j), VR_2D(i,j), VT_2D(i,j), VZ_2D(i,j), &
		  BR_2D(i,j), BT_2D(i,j), BZ_2D(i,j), JR_2D(i,j), JT_2D(i,j), JZ_2D(i,j), &
		  NE_DRRZ(i,j), NE_DZRZ(i,j), NI_DRRZ(i,j), NI_DZRZ(i,j), P_DRRZ(i,j), P_DZRZ(i,j), &
		  PE_DRRZ(i,j), PE_DZRZ(i,j), PI_DRRZ(i,j), PI_DZRZ(i,j), VR_DRRZ(i,j), VR_DZRZ(i,j), &
		  VT_DRRZ(i,j), VT_DZRZ(i,j), VZ_DRRZ(i,j), VZ_DZRZ(i,j), BR_DRRZ(i,j), BR_DZRZ(i,j), &
		  BT_DRRZ(i,j), BT_DZRZ(i,j), BZ_DRRZ(i,j), BZ_DZRZ(i,j), JR_DRRZ(i,j), JR_DZRZ(i,j), &
		  JT_DRRZ(i,j), JT_DZRZ(i,j), JZ_DRRZ(i,j), JZ_DZRZ(i,j), PSI_DRRZ(i,j), PSI_DZRZ(i,j), &
		  ERR_R_2D(i,j), ERR_T_2D(i,j), ERR_Z_2D(i,j), DIVB_2D(i,j)
	enddo
	enddo
   15 format(53(1x,e17.9))
   	close(14)

! read eq_pgfile_rzlim_convex.dat, convex limter grids, number is nxzs
	open(16,file='eq_pgfile_rzlim_convex.dat')
	do i=1,nxzs
	read(16,17) xxs5(i,5), zzs5(i,5)
	xxs5(i,1:4)=xxs5(i,5)
	zzs5(i,1:4)=zzs5(i,5)
	enddo
   17 format(2(1x,e17.9))

   	close(16)
! after read the files, transform the necessary data into clt form
	aa=1 ! normalized minor radius
	aa2=aa*aa ! normalized minor radius
!	xzero=(maxval(RBS)+minval(RBS))/2.d0
	xmax=RRMAX
	xmin=RRMIN
	zmax=ZZMAX
	zmin=ZZMIN
	xmg=RMAXIS
	zmg=ZMAXIS
	xzero=xmg
	if(initia_from_pgfile) zzero=zmg
	xdim=(xmax-xmin)/2.d0
	zdim=(zmax-zmin)/2.d0
	psia=PSIBRY
	psmin=PSIMAG
	psmax=PSIBRY
	qmin=minval(Q_1D)
	qmax=maxval(Q_1D)
	q0=Q_1D(1)
!	pssmw=5.*(psia-psmin)/10.d0


	
	rh(:,:)=NI_2D(:,:)

	if(rho_unif) rh(:,:)=1.d0

	pt(:,:)=PRES_2D(:,:)
	ux(:,:)=VR_2D(:,:)
	uy(:,:)=VT_2D(:,:)
	uz(:,:)=VZ_2D(:,:)
	bx(:,:)=BR_2D(:,:)
	by(:,:)=BT_2D(:,:)
	bz(:,:)=BZ_2D(:,:)

	cx(:,:)=JR_2D(:,:)
	cy(:,:)=JT_2D(:,:)
	cz(:,:)=JZ_2D(:,:)

	rhdx(:,:)=NI_DRRZ(:,:)
	rhdz(:,:)=NI_DZRZ(:,:)
	ptdx(:,:)=P_DRRZ(:,:)
	ptdz(:,:)=P_DZRZ(:,:)

	uxdx(:,:)=VR_DRRZ(:,:)
	uxdz(:,:)=VR_DZRZ(:,:)
	uydx(:,:)=VT_DRRZ(:,:)
	uydz(:,:)=VT_DZRZ(:,:)
	uzdx(:,:)=VZ_DRRZ(:,:)
	uzdz(:,:)=VZ_DZRZ(:,:)

	bxdx(:,:)=BR_DRRZ(:,:)
	bxdz(:,:)=BR_DZRZ(:,:)
	bydx(:,:)=BT_DRRZ(:,:)
	bydz(:,:)=BT_DZRZ(:,:)
	bzdx(:,:)=BZ_DRRZ(:,:)
	bzdz(:,:)=BZ_DZRZ(:,:)

	cx_dx(:,:)=JR_DRRZ(:,:)
	cx_dz(:,:)=JR_DZRZ(:,:)
	cy_dx(:,:)=JT_DRRZ(:,:)
	cy_dz(:,:)=JT_DZRZ(:,:)
	cz_dx(:,:)=JZ_DRRZ(:,:)
	cz_dz(:,:)=JZ_DZRZ(:,:)


	omrot(:,:)=0.d0
	omprot(:,:)=0.d0
	
	weight=1.
	psiam=weight*psia+(1-weight)*PSI_1D(npsi_pgfile-1)
	psia1=psiam

	pst(:,:)=PSI_2D(:,:)
	pst_dx(:,:)=PSI_DRRZ(:,:)
	pst_dz(:,:)=PSI_DZRZ(:,:)
	tht(:,:)=THT_2D(:,:)

      do jx=1,mxt
      do jz=1,mzt
      bpol(jx,jz)=sqrt(bx(jx,jz)**2+bz(jx,jz)**2)
      rr2t(jx,jz)=(pst(jx,jz)-psmin)/(psia-psmin)
      rrt(jx,jz)=sqrt(rr2t(jx,jz))

      tpt(jx,jz)=atan2(pst_dz(jx,jz),pst_dx(jx,jz))
      if(tpt(jx,jz).lt.0) tpt(jx,jz)=tpt(jx,jz)+2*pi

      enddo
      enddo


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
 
!     ps(mpsa)=psia
!     ps(mpsa-2)=psb2min
!     ps(mpsa-1)=(ps(mpsa)+ps(mpsa-2))/2
     dps=(psia-psb2min)/2
!     do js=mpsa-3,mpsa-nda,-1
!     ps(js)=ps(js+1)-dps
!     enddo

     pssm=psia-5.*dps
     pssmw=5.*dps



      open(unit=204,file='gridts.dat',status='unknown',form='formatted')
      write(204,200)( (tpt(jx,jz),jx=1,mxt),jz=1,mzt)
 200   format((1x,e12.5))
      close(204)

      open(unit=205,file='pst_qpg.dat',status='unknown',form='formatted')
      write(205,100)((pst(jx,jz),pst_dx(jx,jz),pst_dz(jx,jz),QRZ_2D(jx,jz),pt(jx,jz),0,0,0,jx=1,mxt),jz=1,mzt)
 100  format(8(1x,e12.5))
      close(205)



	do j=2,npsi_pgfile-1
	if(qmode.ge.Q_1D(j-1) .and. qmode.le.Q_1D(j)) then

	if(j==2) then
      call interp1d3l(PSI_1D(1),PSI_1D(2),PSI_1D(3),PSI_1D(4), &
                    Q_1D(1),Q_1D(2),Q_1D(3),Q_1D(4),qmode,psmode)
      else

      call interp1d3l(PSI_1D(j-2),PSI_1D(j-1),PSI_1D(j),PSI_1D(j+1), &
                    Q_1D(j-2),Q_1D(j-1),Q_1D(j),Q_1D(j+1),qmode,psmode)
	endif

	if(nrank.eq.0) print*,'24321,psmode=',psmode,qmode,cx(101,200),cy(101,200),cz(101,200)
		
	endif

	enddo


	return
	end

!hw**********************************************************************
     subroutine map_nova_to_bnd_grd
	use declare
      integer, parameter :: nq = 33
      integer, parameter :: nr = 85 !int( sqrt(ndat/3) )
      integer, parameter :: nw = 39
    !
      integer lcell(nr,nr)
      integer lnext(ndat),lnext12(ndat12),lnext34(ndat34)
      real*8 risq(ndat),risq12(ndat12),risq34(ndat34)
      real*8 aw(5,ndat),aw12(5,ndat12),aw34(5,ndat34)
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

	 real*8, dimension(ndat) :: tmp_nova
	 real*8, dimension(nbndx) :: tmp_bndx_out1, tmp_bndx_out2, tmp_bndx_out3
	 real*8, dimension(nbndz) :: tmp_bndz_out1, tmp_bndz_out2, tmp_bndz_out3
      include 'mpif.h'
!  d1fc= d f / dx  with fourth-order accuracy central difference
      d1fc(fm2,fm1,f0,fp1,fp2,a,b,c,d)= &
       a*(fp1-f0)+b*(f0-fm1)+c*(fp2-f0)+d*(f0-fm2)   
!-----------------------------------------------------------

! allocate the bndgrds
	allocate(bndx_grd(nbndx,n7))
	allocate(bndz_grd(nbndz,n7))
	bndx_grd(:,:)=bnd_x(1:nbndx,:)
	bndz_grd(:,:)=bnd_z(1:nbndz,:)

	allocate(bx_bndx(nbndx))
	allocate(bx_bndz(nbndz))	
	allocate(bz_bndx(nbndx))
	allocate(bz_bndz(nbndz))
	allocate(bxdx_bndx(nbndx))
	allocate(bxdx_bndz(nbndz))
	allocate(bxdz_bndx(nbndx))
	allocate(bxdz_bndz(nbndz))
	allocate(bzdx_bndx(nbndx))
	allocate(bzdx_bndz(nbndz))
	allocate(bzdz_bndx(nbndx))
	allocate(bzdz_bndz(nbndz))
	allocate(by_bndx(nbndx))
	allocate(by_bndz(nbndz))
	allocate(bydx_bndx(nbndx))
	allocate(bydx_bndz(nbndz))
	allocate(bydz_bndx(nbndx))
	allocate(bydz_bndz(nbndz))
	allocate(pdx_bndx(nbndx))
	allocate(pdx_bndz(nbndz))
	allocate(pdz_bndx(nbndx))
	allocate(pdz_bndz(nbndz))
	allocate(cy_bndx(nbndx))
	allocate(cy_bndz(nbndz))
	allocate(cx_bndx(nbndx))
	allocate(cx_bndz(nbndz))
	allocate(cz_bndx(nbndx))
	allocate(cz_bndz(nbndz))
	allocate(uy_bndx(nbndx))
	allocate(uy_bndz(nbndz))
	allocate(uydx_bndx(nbndx))
	allocate(uydx_bndz(nbndz))
	allocate(uydz_bndx(nbndx))
	allocate(uydz_bndz(nbndz))
	allocate(pt_bndx(nbndx))
	allocate(pt_bndz(nbndz))
	allocate(ptdx_bndx(nbndx))
	allocate(ptdx_bndz(nbndz))
	allocate(ptdz_bndx(nbndx))
	allocate(ptdz_bndz(nbndz))
	allocate(rh_bndx(nbndx))
	allocate(rh_bndz(nbndz))
	allocate(rhdx_bndx(nbndx))
	allocate(rhdx_bndz(nbndz))
	allocate(rhdz_bndx(nbndx))
	allocate(rhdz_bndz(nbndz))

	allocate(x_8bndx(nbndx,my,8))
	allocate(x_8bndz(nbndz,my,8))
	allocate(x1_8bndx(nbndx,my,8))
	allocate(x1_8bndz(nbndz,my,8))
	allocate(xint_8bndx(nbndx,8))
	allocate(xint_8bndz(nbndz,8))
	allocate(cur_3bndx(nbndx,my,3))
	allocate(cur_3bndz(nbndz,my,3))
	allocate(cint_3bndx(nbndx,3))
	allocate(cint_3bndz(nbndz,3))
	allocate(ef_3bndx(nbndx,my,3))
	allocate(ef_3bndz(nbndz,my,3))

	allocate(updated_bndx(nbndx))
	allocate(updated_bndz(nbndz))

	allocate(xy_8bndx(nbndx,my,8))
	allocate(xy2_8bndx(nbndx,my,8))
	allocate(xy_8bndz(nbndz,my,8))
	allocate(xy2_8bndz(nbndz,my,8))


! interpolation, input ndat, xx_nova, zz_nova, f_nova, 
! output f(jx,jz), and its x & z directional derivatives at xxt and zzt grids.
	if(firstmap) then
!!!hw: bx_bndx & bx_bndz  use the grid point of bndx_grd(:,1:2) bndz_grd(:,1:2)
!      call qshep2 ( ndat, xx_nova, zz_nova, bx_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
!        dxi, dzi, rimax, risq, aw, ier )
!      write(*,*) ier
!      do jx=1,nbndx
!      call qs2grd ( bndx_grd(jx,1), bndx_grd(jx,2), ndat, xx_nova, zz_nova, bx_nova, nr, lcell, lnext, ximin, &
!        zimin, dxi, dzi, rimax, risq, aw, bx_bndx(jx), bx_bndx_dx(jx), bx_bndx_dz(jx), ier )
!      enddo
!      do jz=1,nbndz
!      call qs2grd ( bndz_grd(jz,1), bndz_grd(jz,2), ndat, xx_nova, zz_nova, bx_nova, nr, lcell, lnext, ximin, &
!        zimin, dxi, dzi, rimax, risq, aw, bx_bndx(jz), bx_bndz_dx(jz), bx_bndz_dz(jz), ier )
!      enddo


        do itag=1,23

        if(itag.eq.1)  tmp_nova(:)=bx_nova(:)
        if(itag.eq.2)  tmp_nova(:)=bz_nova(:)
        if(itag.eq.3)  tmp_nova(:)=bxdx_nova(:)
        if(itag.eq.4)  tmp_nova(:)=bxdz_nova(:)
        if(itag.eq.5)  tmp_nova(:)=bzdx_nova(:)
        if(itag.eq.6)  tmp_nova(:)=bzdz_nova(:)
        if(itag.eq.7)  tmp_nova(:)=by_nova(:)
        if(itag.eq.8)  tmp_nova(:)=bydx_nova(:)
        if(itag.eq.9)  tmp_nova(:)=bydz_nova(:)
        if(itag.eq.10)  tmp_nova(:)=pdx_nova(:)
        if(itag.eq.11)  tmp_nova(:)=pdz_nova(:)
        if(itag.eq.12)  tmp_nova(:)=cy_nova(:)
        if(itag.eq.13)  tmp_nova(:)=cx_nova(:)
        if(itag.eq.14)  tmp_nova(:)=cz_nova(:)
        if(itag.eq.15)  tmp_nova(:)=uy_nova(:)
        if(itag.eq.16)  tmp_nova(:)=uydx_nova(:)
        if(itag.eq.17)  tmp_nova(:)=uydz_nova(:)
        if(itag.eq.18)  tmp_nova(:)=pt_nova(:)
        if(itag.eq.19)  tmp_nova(:)=ptdx_nova(:)
        if(itag.eq.20)  tmp_nova(:)=ptdz_nova(:)
        if(itag.eq.21)  tmp_nova(:)=rh_nova(:)
        if(itag.eq.22)  tmp_nova(:)=rhdx_nova(:)
        if(itag.eq.23)  tmp_nova(:)=rhdz_nova(:)

!	select case(itag)
!	case(1)
!		  tmp_nova(:)=bx_nova(:)
!	case(2)
!		  tmp_nova(:)=bz_nova(:)
!	case(3)
!		  tmp_nova(:)=bxdx_nova(:)
!	case(4)
!		  tmp_nova(:)=bxdz_nova(:)
!	case(5)
!		  tmp_nova(:)=bzdx_nova(:)
!	case(6)
!		  tmp_nova(:)=bzdz_nova(:)
!	case(7)
!		  tmp_nova(:)=by_nova(:)
!	case(8)
!		  tmp_nova(:)=bydx_nova(:)
!	case(9)
!		  tmp_nova(:)=bydz_nova(:)
!	case(10)
!		  tmp_nova(:)=pdx_nova(:)
!	case(11)
!		  tmp_nova(:)=pdz_nova(:)
!	case(12)
!		  tmp_nova(:)=cy_nova(:)
!	case(13)
!		  tmp_nova(:)=cx_nova(:)
!	case(14)
!		  tmp_nova(:)=cz_nova(:)
!	case(15)
!		  tmp_nova(:)=uy_nova(:)
!	case(16)
!		  tmp_nova(:)=uydx_nova(:)
!	case(17)
!		  tmp_nova(:)=uydz_nova(:)
!	case(18)
!		  tmp_nova(:)=pt_nova(:)
!	case(19)
!		  tmp_nova(:)=ptdx_nova(:)
!	case(20)
!		  tmp_nova(:)=ptdz_nova(:)
!	case(21)
!		  tmp_nova(:)=rh_nova(:)
!	case(22)
!		  tmp_nova(:)=rhdx_nova(:)
!	case(23)
!		  tmp_nova(:)=rhdz_nova(:)
!	end select

      call qshep2 ( ndat, xx_nova, zz_nova, tmp_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )


      do jx=1,nbndx
      call qs2grd ( bndx_grd(jx,1), bndx_grd(jx,2), ndat, xx_nova, zz_nova, tmp_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, tmp_bndx_out1(jx), tmp_bndx_out2(jx), tmp_bndx_out3(jx), ier )
      enddo
      write(*,*) 'itag=', itag, 'jx=', jx, ier
      do jz=1,nbndz
      call qs2grd ( bndz_grd(jz,1), bndz_grd(jz,2), ndat, xx_nova, zz_nova, tmp_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, tmp_bndz_out1(jz), tmp_bndz_out2(jz), tmp_bndz_out3(jz), ier )
      enddo
      write(*,*) 'itag=', itag, 'jz=', jz, ier



        if(itag.eq.1)  bx_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.1)  bx_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.2)  bz_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.2)  bz_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.3)  bxdx_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.3)  bxdx_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.4)  bxdz_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.4)  bxdz_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.5)  bzdx_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.5)  bzdx_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.6)  bzdz_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.6)  bzdz_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.7)  by_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.7)  by_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.8)  bydx_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.8)  bydx_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.9)  bydz_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.9)  bydz_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.10)  pdx_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.10)  pdx_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.11)  pdz_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.11)  pdz_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.12)  cy_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.12)  cy_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.13)  cx_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.13)  cx_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.14)  cz_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.14)  cz_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.15)  uy_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.15)  uy_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.16)  uydx_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.16)  uydx_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.17)  uydz_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.17)  uydz_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.18)  pt_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.18)  pt_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.19)  ptdx_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.19)  ptdx_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.20)  ptdz_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.20)  ptdz_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.21)  rh_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.21)  rh_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.22)  rhdx_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.22)  rhdx_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.23)  rhdz_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.23)  rhdz_bndz(:)=tmp_bndz_out1(:)


!	select case(itag)
!	case(1)
!		  bx_bndx(:)=tmp_bndx_out1(:)
!		  bx_bndz(:)=tmp_bndz_out1(:)
!	case(2)
!		  bz_bndx(:)=tmp_bndx_out1(:)
!		  bz_bndz(:)=tmp_bndz_out1(:)
!	case(3)
!		  bxdx_bndx(:)=tmp_bndx_out1(:)
!		  bxdx_bndz(:)=tmp_bndz_out1(:)
!	case(4)
!		  bxdz_bndx(:)=tmp_bndx_out1(:)
!		  bxdz_bndz(:)=tmp_bndz_out1(:)
!	case(5)
!		  bzdx_bndx(:)=tmp_bndx_out1(:)
!		  bzdx_bndz(:)=tmp_bndz_out1(:)
!	case(6)
!		  bzdz_bndx(:)=tmp_bndx_out1(:)
!		  bzdz_bndz(:)=tmp_bndz_out1(:)
!	case(7)
!		  by_bndx(:)=tmp_bndx_out1(:)
!		  by_bndz(:)=tmp_bndz_out1(:)
!	case(8)
!		  bydx_bndx(:)=tmp_bndx_out1(:)
!		  bydx_bndz(:)=tmp_bndz_out1(:)
!	case(9)
!		  bydz_bndx(:)=tmp_bndx_out1(:)
!		  bydz_bndz(:)=tmp_bndz_out1(:)
!	case(10)
!		  pdx_bndx(:)=tmp_bndx_out1(:)
!		  pdx_bndz(:)=tmp_bndz_out1(:)
!	case(11)
!		  pdz_bndx(:)=tmp_bndx_out1(:)
!		  pdz_bndz(:)=tmp_bndz_out1(:)
!	case(12)
!		  cy_bndx(:)=tmp_bndx_out1(:)
!		  cy_bndz(:)=tmp_bndz_out1(:)
!	case(13)
!		  cx_bndx(:)=tmp_bndx_out1(:)
!		  cx_bndz(:)=tmp_bndz_out1(:)
!	case(14)
!		  cz_bndx(:)=tmp_bndx_out1(:)
!		  cz_bndz(:)=tmp_bndz_out1(:)
!	case(15)
!		  uy_bndx(:)=tmp_bndx_out1(:)
!		  uy_bndz(:)=tmp_bndz_out1(:)
!	case(16)
!		  uydx_bndx(:)=tmp_bndx_out1(:)
!		  uydx_bndz(:)=tmp_bndz_out1(:)
!	case(17)
!		  uydz_bndx(:)=tmp_bndx_out1(:)
!		  uydz_bndz(:)=tmp_bndz_out1(:)
!	case(18)
!		  pt_bndx(:)=tmp_bndx_out1(:)
!		  pt_bndz(:)=tmp_bndz_out1(:)
!	case(19)
!		  ptdx_bndx(:)=tmp_bndx_out1(:)
!		  ptdx_bndz(:)=tmp_bndz_out1(:)
!	case(20)
!		  ptdz_bndx(:)=tmp_bndx_out1(:)
!		  ptdz_bndz(:)=tmp_bndz_out1(:)
!	case(21)
!		  rh_bndx(:)=tmp_bndx_out1(:)
!		  rh_bndz(:)=tmp_bndz_out1(:)
!	case(22)
!		  rhdx_bndx(:)=tmp_bndx_out1(:)
!		  rhdx_bndz(:)=tmp_bndz_out1(:)
!	case(23)
!		  rhdz_bndx(:)=tmp_bndx_out1(:)
!		  rhdz_bndz(:)=tmp_bndz_out1(:)
!	end select

	enddo

	xint_8bndx(:,1)=rh_bndx(:)
	xint_8bndx(:,2)=pt_bndx(:)
	xint_8bndx(:,3)=0.d0
	xint_8bndx(:,4)=uy_bndx(:)
	xint_8bndx(:,5)=0.d0
	xint_8bndx(:,6)=bx_bndx(:)
	xint_8bndx(:,7)=by_bndx(:)
	xint_8bndx(:,8)=bz_bndx(:)
	cint_3bndx(:,1)=cx_bndx(:)
	cint_3bndx(:,2)=cy_bndx(:)
	cint_3bndx(:,3)=cz_bndx(:)


	xint_8bndz(:,1)=rh_bndz(:)
	xint_8bndz(:,2)=pt_bndz(:)
	xint_8bndz(:,3)=0.d0
	xint_8bndz(:,4)=uy_bndz(:)
	xint_8bndz(:,5)=0.d0
	xint_8bndz(:,6)=bx_bndz(:)
	xint_8bndz(:,7)=by_bndz(:)
	xint_8bndz(:,8)=bz_bndz(:)
	cint_3bndz(:,1)=cx_bndz(:)
	cint_3bndz(:,2)=cy_bndz(:)
	cint_3bndz(:,3)=cz_bndz(:)


	do jy=iy_first,iy_last
	x_8bndx(:,jy,:)=xint_8bndx(:,:)
	cur_3bndx(:,jy,:)=0.d0
	x1_8bndx(:,jy,:)=0.d0
	ef_3bndx(:,jy,:)=0.d0

	x_8bndz(:,jy,:)=xint_8bndz(:,:)
	cur_3bndz(:,jy,:)=0.d0
	x1_8bndz(:,jy,:)=0.d0
	ef_3bndz(:,jy,:)=0.d0

	enddo


	if(nrank.eq.0) then
	open(unit=1,file='bx_bndz.dat',status='unknown',form='formatted')
	write(1,2) (bndz_grd(jx,1),bndz_grd(jx,2),bx_bndz(jx),jx=1,nbndz)
    2 format(3(1x,e12.5))
    	close(1)
	endif


	else
	endif

	return
	end

!hw**************************************************************
     subroutine map_nova_to_bnd_grd_in_each_proc
	use declare
      integer, parameter :: nq = 33
      integer, parameter :: nr = 85 !int( sqrt(ndat/3) )
      integer, parameter :: nw = 39
    !
      integer lcell(nr,nr)
      integer lnext(ndat),lnext12(ndat12),lnext34(ndat34)
      real*8 risq(ndat),risq12(ndat12),risq34(ndat34)
      real*8 aw(5,ndat),aw12(5,ndat12),aw34(5,ndat34)
      real*8 rimax,ximin,zimin,dxi,dzi
!      real*8, dimension(mxt,mzt) :: bx_dx,bx_dz,bz_dx,bz_dz,by_dx,by_dz,uy_dx,uy_dz,tht_dx,tht_dz
 !     real*8, dimension(mxt,mzt) :: pt_dx,pt_dz,rh_dx,rh_dz
!	real*8, dimension(nbndx) :: bx_bndx_dx, bx_bndx_dz
!      real*8, dimension(mxt,mzt) :: bx,bxdx,bxdz,bz,bzdx,bzdz
!      real*8, dimension(mxt,mzt) :: bxdx_dx,bxdx_dz,bxdz_dx,bxdz_dz,bzdx_dx,bzdx_dz,bzdz_dx,bzdz_dz
!      integer icell(1,1)
!      integer inext(9)
!      real*8 xin(9),zin(9),qin(9),rinsq(9),ain(5,9)
!       real*8 xout,zout,qout,qdxout,qdzout
       integer iam,iap, itag
	 integer ier

	 real*8, dimension(ndat) :: tmp_nova
	 real*8, dimension(nbndx_ep) :: tmp_bndx_out1, tmp_bndx_out2, tmp_bndx_out3
	 real*8, dimension(nbndz_ep) :: tmp_bndz_out1, tmp_bndz_out2, tmp_bndz_out3
      include 'mpif.h'
!  d1fc= d f / dx  with fourth-order accuracy central difference
      d1fc(fm2,fm1,f0,fp1,fp2,a,b,c,d)= &
       a*(fp1-f0)+b*(f0-fm1)+c*(fp2-f0)+d*(f0-fm2)   
!-----------------------------------------------------------

! allocate the bndgrds
	allocate(bndx_grd_ep(nbndx_ep,n7))
	allocate(bndz_grd_ep(nbndz_ep,n7))
	bndx_grd_ep(:,:)=bnd_x_ep(1:nbndx,:)
	bndz_grd_ep(:,:)=bnd_z_ep(1:nbndz,:)

	allocate(bx_bndx_ep(nbndx_ep))
	allocate(bx_bndz_ep(nbndz_ep))	
	allocate(bz_bndx_ep(nbndx_ep))
	allocate(bz_bndz_ep(nbndz_ep))
	allocate(bxdx_bndx_ep(nbndx_ep))
	allocate(bxdx_bndz_ep(nbndz_ep))
	allocate(bxdz_bndx_ep(nbndx_ep))
	allocate(bxdz_bndz_ep(nbndz_ep))
	allocate(bzdx_bndx_ep(nbndx_ep))
	allocate(bzdx_bndz_ep(nbndz_ep))
	allocate(bzdz_bndx_ep(nbndx_ep))
	allocate(bzdz_bndz_ep(nbndz_ep))
	allocate(by_bndx_ep(nbndx_ep))
	allocate(by_bndz_ep(nbndz_ep))
	allocate(bydx_bndx_ep(nbndx_ep))
	allocate(bydx_bndz_ep(nbndz_ep))
	allocate(bydz_bndx_ep(nbndx_ep))
	allocate(bydz_bndz_ep(nbndz_ep))
	allocate(pdx_bndx_ep(nbndx_ep))
	allocate(pdx_bndz_ep(nbndz_ep))
	allocate(pdz_bndx_ep(nbndx_ep))
	allocate(pdz_bndz_ep(nbndz_ep))
	allocate(cy_bndx_ep(nbndx_ep))
	allocate(cy_bndz_ep(nbndz_ep))
	allocate(cx_bndx_ep(nbndx_ep))
	allocate(cx_bndz_ep(nbndz_ep))
	allocate(cz_bndx_ep(nbndx_ep))
	allocate(cz_bndz_ep(nbndz_ep))
	allocate(uy_bndx_ep(nbndx_ep))
	allocate(uy_bndz_ep(nbndz_ep))
	allocate(uydx_bndx_ep(nbndx_ep))
	allocate(uydx_bndz_ep(nbndz_ep))
	allocate(uydz_bndx_ep(nbndx_ep))
	allocate(uydz_bndz_ep(nbndz_ep))
	allocate(pt_bndx_ep(nbndx_ep))
	allocate(pt_bndz_ep(nbndz_ep))
	allocate(ptdx_bndx_ep(nbndx_ep))
	allocate(ptdx_bndz_ep(nbndz_ep))
	allocate(ptdz_bndx_ep(nbndx_ep))
	allocate(ptdz_bndz_ep(nbndz_ep))
	allocate(rh_bndx_ep(nbndx_ep))
	allocate(rh_bndz_ep(nbndz_ep))
	allocate(rhdx_bndx_ep(nbndx_ep))
	allocate(rhdx_bndz_ep(nbndz_ep))
	allocate(rhdz_bndx_ep(nbndx_ep))
	allocate(rhdz_bndz_ep(nbndz_ep))


! interpolation, input ndat, xx_nova, zz_nova, f_nova, 
! output f(jx,jz), and its x & z directional derivatives at xxt and zzt grids.
	if(firstmap) then
!!!hw: bx_bndx & bx_bndz  use the grid point of bndx_grd(:,1:2) bndz_grd(:,1:2)
!      call qshep2 ( ndat, xx_nova, zz_nova, bx_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
!        dxi, dzi, rimax, risq, aw, ier )
!      write(*,*) ier
!      do jx=1,nbndx
!      call qs2grd ( bndx_grd(jx,1), bndx_grd(jx,2), ndat, xx_nova, zz_nova, bx_nova, nr, lcell, lnext, ximin, &
!        zimin, dxi, dzi, rimax, risq, aw, bx_bndx(jx), bx_bndx_dx(jx), bx_bndx_dz(jx), ier )
!      enddo
!      do jz=1,nbndz
!      call qs2grd ( bndz_grd(jz,1), bndz_grd(jz,2), ndat, xx_nova, zz_nova, bx_nova, nr, lcell, lnext, ximin, &
!        zimin, dxi, dzi, rimax, risq, aw, bx_bndx(jz), bx_bndz_dx(jz), bx_bndz_dz(jz), ier )
!      enddo


	do itag=1,23

!	select case(itag)
!	case(1)
!		  tmp_nova(:)=bx_nova(:)
!	case(2)
!		  tmp_nova(:)=bz_nova(:)
!	case(3)
!		  tmp_nova(:)=bxdx_nova(:)
!	case(4)
!		  tmp_nova(:)=bxdz_nova(:)
!	case(5)
!		  tmp_nova(:)=bzdx_nova(:)
!	case(6)
!		  tmp_nova(:)=bzdz_nova(:)
!	case(7)
!		  tmp_nova(:)=by_nova(:)
!	case(8)
!		  tmp_nova(:)=bydx_nova(:)
!	case(9)
!		  tmp_nova(:)=bydz_nova(:)
!	case(10)
!		  tmp_nova(:)=pdx_nova(:)
!	case(11)
!		  tmp_nova(:)=pdz_nova(:)
!	case(12)
!		  tmp_nova(:)=cy_nova(:)
!	case(13)
!		  tmp_nova(:)=cx_nova(:)
!	case(14)
!		  tmp_nova(:)=cz_nova(:)
!	case(15)
!		  tmp_nova(:)=uy_nova(:)
!	case(16)
!		  tmp_nova(:)=uydx_nova(:)
!	case(17)
!		  tmp_nova(:)=uydz_nova(:)
!	case(18)
!		  tmp_nova(:)=pt_nova(:)
!	case(19)
!		  tmp_nova(:)=ptdx_nova(:)
!	case(20)
!		  tmp_nova(:)=ptdz_nova(:)
!	case(21)
!		  tmp_nova(:)=rh_nova(:)
!	case(22)
!		  tmp_nova(:)=rhdx_nova(:)
!	case(23)
!		  tmp_nova(:)=rhdz_nova(:)
!	end select

	if(itag.eq.1)    tmp_nova(:)=bx_nova(:)
	if(itag.eq.2)    tmp_nova(:)=bz_nova(:)
	if(itag.eq.3)    tmp_nova(:)=bxdx_nova(:)
	if(itag.eq.4)    tmp_nova(:)=bxdz_nova(:)
	if(itag.eq.5)    tmp_nova(:)=bzdx_nova(:)
	if(itag.eq.6)    tmp_nova(:)=bzdz_nova(:)
	if(itag.eq.7)    tmp_nova(:)=by_nova(:)
	if(itag.eq.8)    tmp_nova(:)=bydx_nova(:)
	if(itag.eq.9)    tmp_nova(:)=bydz_nova(:)
	if(itag.eq.10)   tmp_nova(:)=pdx_nova(:)
	if(itag.eq.11)   tmp_nova(:)=pdz_nova(:)
	if(itag.eq.12)   tmp_nova(:)=cy_nova(:)
	if(itag.eq.13)   tmp_nova(:)=cx_nova(:)
	if(itag.eq.14)   tmp_nova(:)=cz_nova(:)
	if(itag.eq.15)   tmp_nova(:)=uy_nova(:)
	if(itag.eq.16)   tmp_nova(:)=uydx_nova(:)
	if(itag.eq.17)   tmp_nova(:)=uydz_nova(:)
	if(itag.eq.18)   tmp_nova(:)=pt_nova(:)
	if(itag.eq.19)   tmp_nova(:)=ptdx_nova(:)
	if(itag.eq.20)   tmp_nova(:)=ptdz_nova(:)
	if(itag.eq.21)   tmp_nova(:)=rh_nova(:)
	if(itag.eq.22)   tmp_nova(:)=rhdx_nova(:)
	if(itag.eq.23)   tmp_nova(:)=rhdz_nova(:)


      call qshep2 ( ndat, xx_nova, zz_nova, tmp_nova, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )
      do jx=1,nbndx_ep
      call qs2grd ( bndx_grd_ep(jx,1), bndx_grd_ep(jx,2), ndat, xx_nova, zz_nova, tmp_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, tmp_bndx_out1(jx), tmp_bndx_out2(jx), tmp_bndx_out3(jx), ier )
      write(*,*) 'itag=', itag, 'jx=', jx, ier
      enddo
      do jz=1,nbndz_ep
      call qs2grd ( bndz_grd_ep(jz,1), bndz_grd_ep(jz,2), ndat, xx_nova, zz_nova, tmp_nova, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, tmp_bndz_out1(jz), tmp_bndz_out2(jz), tmp_bndz_out3(jz), ier )
      write(*,*) 'itag=', itag, 'jz=', jz, ier
      enddo


!	select case(itag)
!	case(1)
!		  bx_bndx_ep(:)=tmp_bndx_out1(:)
!		  bx_bndz_ep(:)=tmp_bndz_out1(:)
!	case(2)
!		  bz_bndx_ep(:)=tmp_bndx_out1(:)
!		  bz_bndz_ep(:)=tmp_bndz_out1(:)
!	case(3)
!		  bxdx_bndx_ep(:)=tmp_bndx_out1(:)
!		  bxdx_bndz_ep(:)=tmp_bndz_out1(:)
!	case(4)
!		  bxdz_bndx_ep(:)=tmp_bndx_out1(:)
!		  bxdz_bndz_ep(:)=tmp_bndz_out1(:)
!	case(5)
!		  bzdx_bndx_ep(:)=tmp_bndx_out1(:)
!		  bzdx_bndz_ep(:)=tmp_bndz_out1(:)
!	case(6)
!		  bzdz_bndx_ep(:)=tmp_bndx_out1(:)
!		  bzdz_bndz_ep(:)=tmp_bndz_out1(:)
!	case(7)
!		  by_bndx_ep(:)=tmp_bndx_out1(:)
!		  by_bndz_ep(:)=tmp_bndz_out1(:)
!	case(8)
!		  bydx_bndx_ep(:)=tmp_bndx_out1(:)
!		  bydx_bndz_ep(:)=tmp_bndz_out1(:)
!	case(9)
!		  bydz_bndx_ep(:)=tmp_bndx_out1(:)
!		  bydz_bndz_ep(:)=tmp_bndz_out1(:)
!	case(10)
!		  pdx_bndx_ep(:)=tmp_bndx_out1(:)
!		  pdx_bndz_ep(:)=tmp_bndz_out1(:)
!	case(11)
!		  pdz_bndx_ep(:)=tmp_bndx_out1(:)
!		  pdz_bndz_ep(:)=tmp_bndz_out1(:)
!	case(12)
!		  cy_bndx_ep(:)=tmp_bndx_out1(:)
!		  cy_bndz_ep(:)=tmp_bndz_out1(:)
!	case(13)
!		  cx_bndx_ep(:)=tmp_bndx_out1(:)
!		  cx_bndz_ep(:)=tmp_bndz_out1(:)
!	case(14)
!		  cz_bndx_ep(:)=tmp_bndx_out1(:)
!		  cz_bndz_ep(:)=tmp_bndz_out1(:)
!	case(15)
!		  uy_bndx_ep(:)=tmp_bndx_out1(:)
!		  uy_bndz_ep(:)=tmp_bndz_out1(:)
!	case(16)
!		  uydx_bndx_ep(:)=tmp_bndx_out1(:)
!		  uydx_bndz_ep(:)=tmp_bndz_out1(:)
!	case(17)
!		  uydz_bndx_ep(:)=tmp_bndx_out1(:)
!		  uydz_bndz_ep(:)=tmp_bndz_out1(:)
!	case(18)
!		  pt_bndx_ep(:)=tmp_bndx_out1(:)
!		  pt_bndz_ep(:)=tmp_bndz_out1(:)
!	case(19)
!		  ptdx_bndx_ep(:)=tmp_bndx_out1(:)
!		  ptdx_bndz_ep(:)=tmp_bndz_out1(:)
!	case(20)
!		  ptdz_bndx_ep(:)=tmp_bndx_out1(:)
!		  ptdz_bndz_ep(:)=tmp_bndz_out1(:)
!	case(21)
!		  rh_bndx_ep(:)=tmp_bndx_out1(:)
!		  rh_bndz_ep(:)=tmp_bndz_out1(:)
!	case(22)
!		  rhdx_bndx_ep(:)=tmp_bndx_out1(:)
!		  rhdx_bndz_ep(:)=tmp_bndz_out1(:)
!	case(23)
!		  rhdz_bndx_ep(:)=tmp_bndx_out1(:)
!		  rhdz_bndz_ep(:)=tmp_bndz_out1(:)
!	end select

        if(itag.eq.1)  bx_bndx_ep(:)=tmp_bndx_out1(:)
        if(itag.eq.1)  bx_bndz_ep(:)=tmp_bndz_out1(:)
        if(itag.eq.2)  bz_bndx_ep(:)=tmp_bndx_out1(:)
        if(itag.eq.2)  bz_bndz_ep(:)=tmp_bndz_out1(:)
        if(itag.eq.3)  bxdx_bndx_ep(:)=tmp_bndx_out1(:)
        if(itag.eq.3)  bxdx_bndz_ep(:)=tmp_bndz_out1(:)
        if(itag.eq.4)  bxdz_bndx_ep(:)=tmp_bndx_out1(:)
        if(itag.eq.4)  bxdz_bndz_ep(:)=tmp_bndz_out1(:)
        if(itag.eq.5)  bzdx_bndx_ep(:)=tmp_bndx_out1(:)
        if(itag.eq.5)  bzdx_bndz_ep(:)=tmp_bndz_out1(:)
        if(itag.eq.6)  bzdz_bndx_ep(:)=tmp_bndx_out1(:)
        if(itag.eq.6)  bzdz_bndz_ep(:)=tmp_bndz_out1(:)
        if(itag.eq.7)  by_bndx_ep(:)=tmp_bndx_out1(:)
        if(itag.eq.7)  by_bndz_ep(:)=tmp_bndz_out1(:)
        if(itag.eq.8)  bydx_bndx_ep(:)=tmp_bndx_out1(:)
        if(itag.eq.8)  bydx_bndz_ep(:)=tmp_bndz_out1(:)
        if(itag.eq.9)  bydz_bndx_ep(:)=tmp_bndx_out1(:)
        if(itag.eq.9)  bydz_bndz_ep(:)=tmp_bndz_out1(:)
        if(itag.eq.10) pdx_bndx_ep(:)=tmp_bndx_out1(:)
        if(itag.eq.10) pdx_bndz_ep(:)=tmp_bndz_out1(:)
        if(itag.eq.11) pdz_bndx_ep(:)=tmp_bndx_out1(:)
        if(itag.eq.11) pdz_bndz_ep(:)=tmp_bndz_out1(:)
        if(itag.eq.12) cy_bndx_ep(:)=tmp_bndx_out1(:)
        if(itag.eq.12) cy_bndz_ep(:)=tmp_bndz_out1(:)
        if(itag.eq.13) cx_bndx_ep(:)=tmp_bndx_out1(:)
        if(itag.eq.13) cx_bndz_ep(:)=tmp_bndz_out1(:)
        if(itag.eq.14) cz_bndx_ep(:)=tmp_bndx_out1(:)
        if(itag.eq.14) cz_bndz_ep(:)=tmp_bndz_out1(:)
        if(itag.eq.15) uy_bndx_ep(:)=tmp_bndx_out1(:)
        if(itag.eq.15) uy_bndz_ep(:)=tmp_bndz_out1(:)
        if(itag.eq.16) uydx_bndx_ep(:)=tmp_bndx_out1(:)	
        if(itag.eq.16) uydx_bndz_ep(:)=tmp_bndz_out1(:)	
        if(itag.eq.17) uydz_bndx_ep(:)=tmp_bndx_out1(:)	
        if(itag.eq.17) uydz_bndz_ep(:)=tmp_bndz_out1(:)	
        if(itag.eq.18) pt_bndx_ep(:)=tmp_bndx_out1(:)
        if(itag.eq.18) pt_bndz_ep(:)=tmp_bndz_out1(:)
        if(itag.eq.19) ptdx_bndx_ep(:)=tmp_bndx_out1(:)	
        if(itag.eq.19) ptdx_bndz_ep(:)=tmp_bndz_out1(:)	
        if(itag.eq.20) ptdz_bndx_ep(:)=tmp_bndx_out1(:)	
        if(itag.eq.20) ptdz_bndz_ep(:)=tmp_bndz_out1(:)	
        if(itag.eq.21) rh_bndx_ep(:)=tmp_bndx_out1(:)
        if(itag.eq.21) rh_bndz_ep(:)=tmp_bndz_out1(:)
        if(itag.eq.22) rhdx_bndx_ep(:)=tmp_bndx_out1(:)	
        if(itag.eq.22) rhdx_bndz_ep(:)=tmp_bndz_out1(:)	
        if(itag.eq.23) rhdz_bndx_ep(:)=tmp_bndx_out1(:)	
        if(itag.eq.23) rhdz_bndz_ep(:)=tmp_bndz_out1(:)	
	

	enddo

!	if(nrank.eq.0) then
!	open(unit=1,file='bx_bndz.dat',status='unknown',form='formatted')
!	write(1,2) ((bndz_grd(jx,1),bndz_grd(jx,2),bx_bndz(jx)),jx=1,nbndz)
!    2 format(3(1x,e12.5))
!    	close(1)
!	endif


	else
	endif

	return
	end

!hw**************************************************************
     subroutine map_int_to_bnd_grd
	use declare
      integer, parameter :: nq = 33
!      integer, parameter :: nr = 150 !int( sqrt(ndat/3) )
      integer, parameter :: nr = int( sqrt(ndata*1./3.) )
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
!-----------------------------------------------------------
! allocate the bndgrds
	allocate(bndx_grd(nbndx,n7))
	allocate(bndz_grd(nbndz,n7))
	bndx_grd(:,:)=bnd_x(1:nbndx,:)
	bndz_grd(:,:)=bnd_z(1:nbndz,:)

	allocate(bx_bndx(nbndx))
	allocate(bx_bndz(nbndz))	
	allocate(bz_bndx(nbndx))
	allocate(bz_bndz(nbndz))
	allocate(bxdx_bndx(nbndx))
	allocate(bxdx_bndz(nbndz))
	allocate(bxdz_bndx(nbndx))
	allocate(bxdz_bndz(nbndz))
	allocate(bzdx_bndx(nbndx))
	allocate(bzdx_bndz(nbndz))
	allocate(bzdz_bndx(nbndx))
	allocate(bzdz_bndz(nbndz))
	allocate(by_bndx(nbndx))
	allocate(by_bndz(nbndz))
	allocate(bydx_bndx(nbndx))
	allocate(bydx_bndz(nbndz))
	allocate(bydz_bndx(nbndx))
	allocate(bydz_bndz(nbndz))
	allocate(pdx_bndx(nbndx))
	allocate(pdx_bndz(nbndz))
	allocate(pdz_bndx(nbndx))
	allocate(pdz_bndz(nbndz))
	allocate(cy_bndx(nbndx))
	allocate(cy_bndz(nbndz))
	allocate(cx_bndx(nbndx))
	allocate(cx_bndz(nbndz))
	allocate(cz_bndx(nbndx))
	allocate(cz_bndz(nbndz))
	allocate(uy_bndx(nbndx))
	allocate(uy_bndz(nbndz))
	allocate(uydx_bndx(nbndx))
	allocate(uydx_bndz(nbndz))
	allocate(uydz_bndx(nbndx))
	allocate(uydz_bndz(nbndz))
	allocate(pt_bndx(nbndx))
	allocate(pt_bndz(nbndz))
	allocate(ptdx_bndx(nbndx))
	allocate(ptdx_bndz(nbndz))
	allocate(ptdz_bndx(nbndx))
	allocate(ptdz_bndz(nbndz))
	allocate(rh_bndx(nbndx))
	allocate(rh_bndz(nbndz))
	allocate(rhdx_bndx(nbndx))
	allocate(rhdx_bndz(nbndz))
	allocate(rhdz_bndx(nbndx))
	allocate(rhdz_bndz(nbndz))

	allocate(x_8bndx(nbndx,my,8))
	allocate(x_8bndz(nbndz,my,8))
	allocate(x1_8bndx(nbndx,my,8))
	allocate(x1_8bndz(nbndz,my,8))
	allocate(xint_8bndx(nbndx,8))
	allocate(xint_8bndz(nbndz,8))
	allocate(cur_3bndx(nbndx,my,3))
	allocate(cur_3bndz(nbndz,my,3))
	allocate(cint_3bndx(nbndx,3))
	allocate(cint_3bndz(nbndz,3))
	allocate(ef_3bndx(nbndx,my,3))
	allocate(ef_3bndz(nbndz,my,3))

	allocate(updated_bndx(nbndx))
	allocate(updated_bndz(nbndz))

	allocate(xy_8bndx(nbndx,my,8))
	allocate(xy2_8bndx(nbndx,my,8))
	allocate(xy_8bndz(nbndz,my,8))
	allocate(xy2_8bndz(nbndz,my,8))
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

	do itag=1,23
        print*, 'itag=',itag

!	select case(itag)
!	case(1)
!		  tmp_int(:)=reshape(bx, [ndata])
!	case(2)
!		  tmp_int(:)=reshape(bz, [ndata])
!	case(3)
!		  tmp_int(:)=reshape(bxdx, [ndata])
!	case(4)
!		  tmp_int(:)=reshape(bxdz, [ndata])
!	case(5)
!		  tmp_int(:)=reshape(bzdx, [ndata])
!	case(6)
!		  tmp_int(:)=reshape(bzdz, [ndata])
!	case(7)
!		  tmp_int(:)=reshape(by, [ndata])
!	case(8)
!		  tmp_int(:)=reshape(bydx, [ndata])
!	case(9)
!		  tmp_int(:)=reshape(bydz, [ndata])
!	case(10)
!		  tmp_int(:)=reshape(ptdx, [ndata])
!	case(11)
!		  tmp_int(:)=reshape(ptdz, [ndata])
!	case(12)
!		  tmp_int(:)=reshape(cy, [ndata])
!	case(13)
!		  tmp_int(:)=reshape(cx, [ndata])
!	case(14)
!		  tmp_int(:)=reshape(cz, [ndata])
!	case(15)
!		  tmp_int(:)=reshape(uy, [ndata])
!	case(16)
!		  tmp_int(:)=reshape(uydx, [ndata])
!	case(17)
!		  tmp_int(:)=reshape(uydz, [ndata])
!	case(18)
!		  tmp_int(:)=reshape(pt, [ndata])
!	case(19)
!		  tmp_int(:)=reshape(ptdx, [ndata])
!	case(20)
!		  tmp_int(:)=reshape(ptdz, [ndata])
!	case(21)
!		  tmp_int(:)=reshape(rh, [ndata])
!	case(22)
!		  tmp_int(:)=reshape(rhdx, [ndata])
!	case(23)
!		  tmp_int(:)=reshape(rhdz, [ndata])
!	end select
	


	if(itag.eq.1)    tmp_int(:)=reshape(bx, [ndata])
	if(itag.eq.2)    tmp_int(:)=reshape(bz, [ndata])
	if(itag.eq.3)    tmp_int(:)=reshape(bxdx, [ndata])
	if(itag.eq.4)    tmp_int(:)=reshape(bxdz, [ndata])
	if(itag.eq.5)    tmp_int(:)=reshape(bzdx, [ndata])
	if(itag.eq.6)    tmp_int(:)=reshape(bzdz, [ndata])
	if(itag.eq.7)    tmp_int(:)=reshape(by, [ndata])
	if(itag.eq.8)    tmp_int(:)=reshape(bydx, [ndata])
	if(itag.eq.9)    tmp_int(:)=reshape(bydz, [ndata])
	if(itag.eq.10)   tmp_int(:)=reshape(ptdx, [ndata])
	if(itag.eq.11)   tmp_int(:)=reshape(ptdz, [ndata])
	if(itag.eq.12)   tmp_int(:)=reshape(cy, [ndata])
	if(itag.eq.13)   tmp_int(:)=reshape(cx, [ndata])
	if(itag.eq.14)   tmp_int(:)=reshape(cz, [ndata])
	if(itag.eq.15)   tmp_int(:)=reshape(uy, [ndata])
	if(itag.eq.16)   tmp_int(:)=reshape(uydx, [ndata])
	if(itag.eq.17)   tmp_int(:)=reshape(uydz, [ndata])
	if(itag.eq.18)   tmp_int(:)=reshape(pt, [ndata])
	if(itag.eq.19)   tmp_int(:)=reshape(ptdx, [ndata])
	if(itag.eq.20)   tmp_int(:)=reshape(ptdz, [ndata])
	if(itag.eq.21)   tmp_int(:)=reshape(rh, [ndata])
	if(itag.eq.22)   tmp_int(:)=reshape(rhdx, [ndata])
	if(itag.eq.23)   tmp_int(:)=reshape(rhdz, [ndata])



      call qshep2 ( ndata, xx_int, zz_int, tmp_int, nq, nw, nr, lcell, lnext, ximin, zimin, &
        dxi, dzi, rimax, risq, aw, ier )

      do jx=1,nbndx
      call qs2grd ( bndx_grd(jx,1), bndx_grd(jx,2), ndata, xx_int, zz_int, tmp_int, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, tmp_bndx_out1(jx), tmp_bndx_out2(jx), tmp_bndx_out3(jx), ier )
      enddo
      write(*,*) 'itag=', itag, 'jx=', jx, ier

      do jz=1,nbndz
      call qs2grd ( bndz_grd(jz,1), bndz_grd(jz,2), ndata, xx_int, zz_int, tmp_int, nr, lcell, lnext, ximin, &
        zimin, dxi, dzi, rimax, risq, aw, tmp_bndz_out1(jz), tmp_bndz_out2(jz), tmp_bndz_out3(jz), ier )
      enddo
      write(*,*) 'itag=', itag, 'jz=', jz, ier

!	select case(itag)
!	case(1)
!		  bx_bndx(:)=tmp_bndx_out1(:)
!		  bx_bndz(:)=tmp_bndz_out1(:)
!	case(2)
!		  bz_bndx(:)=tmp_bndx_out1(:)
!		  bz_bndz(:)=tmp_bndz_out1(:)
!	case(3)
!		  bxdx_bndx(:)=tmp_bndx_out1(:)
!		  bxdx_bndz(:)=tmp_bndz_out1(:)
!	case(4)
!		  bxdz_bndx(:)=tmp_bndx_out1(:)
!		  bxdz_bndz(:)=tmp_bndz_out1(:)
!	case(5)
!		  bzdx_bndx(:)=tmp_bndx_out1(:)
!		  bzdx_bndz(:)=tmp_bndz_out1(:)
!	case(6)
!		  bzdz_bndx(:)=tmp_bndx_out1(:)
!		  bzdz_bndz(:)=tmp_bndz_out1(:)
!	case(7)
!		  by_bndx(:)=tmp_bndx_out1(:)
!		  by_bndz(:)=tmp_bndz_out1(:)
!	case(8)
!		  bydx_bndx(:)=tmp_bndx_out1(:)
!		  bydx_bndz(:)=tmp_bndz_out1(:)
!	case(9)
!		  bydz_bndx(:)=tmp_bndx_out1(:)
!		  bydz_bndz(:)=tmp_bndz_out1(:)
!	case(10)
!		  pdx_bndx(:)=tmp_bndx_out1(:)
!		  pdx_bndz(:)=tmp_bndz_out1(:)
!	case(11)
!		  pdz_bndx(:)=tmp_bndx_out1(:)
!		  pdz_bndz(:)=tmp_bndz_out1(:)
!	case(12)
!		  cy_bndx(:)=tmp_bndx_out1(:)
!		  cy_bndz(:)=tmp_bndz_out1(:)
!	case(13)
!		  cx_bndx(:)=tmp_bndx_out1(:)
!		  cx_bndz(:)=tmp_bndz_out1(:)
!	case(14)
!		  cz_bndx(:)=tmp_bndx_out1(:)
!		  cz_bndz(:)=tmp_bndz_out1(:)
!	case(15)
!		  uy_bndx(:)=tmp_bndx_out1(:)
!		  uy_bndz(:)=tmp_bndz_out1(:)
!	case(16)
!		  uydx_bndx(:)=tmp_bndx_out1(:)
!		  uydx_bndz(:)=tmp_bndz_out1(:)
!	case(17)
!		  uydz_bndx(:)=tmp_bndx_out1(:)
!		  uydz_bndz(:)=tmp_bndz_out1(:)
!	case(18)
!		  pt_bndx(:)=tmp_bndx_out1(:)
!		  pt_bndz(:)=tmp_bndz_out1(:)
!	case(19)
!		  ptdx_bndx(:)=tmp_bndx_out1(:)
!		  ptdx_bndz(:)=tmp_bndz_out1(:)
!	case(20)
!		  ptdz_bndx(:)=tmp_bndx_out1(:)
!		  ptdz_bndz(:)=tmp_bndz_out1(:)
!	case(21)
!		  rh_bndx(:)=tmp_bndx_out1(:)
!		  rh_bndz(:)=tmp_bndz_out1(:)
!	case(22)
!		  rhdx_bndx(:)=tmp_bndx_out1(:)
!		  rhdx_bndz(:)=tmp_bndz_out1(:)
!	case(23)
!		  rhdz_bndx(:)=tmp_bndx_out1(:)
!		  rhdz_bndz(:)=tmp_bndz_out1(:)
!	end select

	
        if(itag.eq.1) 	bx_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.1) 	bx_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.2) 	bz_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.2) 	bz_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.3) 	bxdx_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.3) 	bxdx_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.4) 	bxdz_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.4) 	bxdz_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.5) 	bzdx_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.5) 	bzdx_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.6) 	bzdz_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.6) 	bzdz_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.7) 	by_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.7) 	by_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.8) 	bydx_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.8) 	bydx_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.9) 	bydz_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.9) 	bydz_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.10)	pdx_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.10)	pdx_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.11)	pdz_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.11)	pdz_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.12)	cy_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.12)	cy_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.13)	cx_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.13)	cx_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.14)	cz_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.14)	cz_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.15)	uy_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.15)	uy_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.16)	uydx_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.16)	uydx_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.17)	uydz_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.17)	uydz_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.18)	pt_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.18)	pt_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.19)	ptdx_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.19)	ptdx_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.20)	ptdz_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.20)	ptdz_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.21)	rh_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.21)	rh_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.22)	rhdx_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.22)	rhdx_bndz(:)=tmp_bndz_out1(:)
        if(itag.eq.23)	rhdz_bndx(:)=tmp_bndx_out1(:)
        if(itag.eq.23)	rhdz_bndz(:)=tmp_bndz_out1(:)
	
	enddo

	xint_8bndx(:,1)=rh_bndx(:)
	xint_8bndx(:,2)=pt_bndx(:)
	xint_8bndx(:,3)=0.d0
	xint_8bndx(:,4)=uy_bndx(:)
	xint_8bndx(:,5)=0.d0
	xint_8bndx(:,6)=bx_bndx(:)
	xint_8bndx(:,7)=by_bndx(:)
	xint_8bndx(:,8)=bz_bndx(:)
	cint_3bndx(:,1)=cx_bndx(:)
	cint_3bndx(:,2)=cy_bndx(:)
	cint_3bndx(:,3)=cz_bndx(:)


	xint_8bndz(:,1)=rh_bndz(:)
	xint_8bndz(:,2)=pt_bndz(:)
	xint_8bndz(:,3)=0.d0
	xint_8bndz(:,4)=uy_bndz(:)
	xint_8bndz(:,5)=0.d0
	xint_8bndz(:,6)=bx_bndz(:)
	xint_8bndz(:,7)=by_bndz(:)
	xint_8bndz(:,8)=bz_bndz(:)
	cint_3bndz(:,1)=cx_bndz(:)
	cint_3bndz(:,2)=cy_bndz(:)
	cint_3bndz(:,3)=cz_bndz(:)


	do jy=iy_first,iy_last
	x_8bndx(:,jy,:)=xint_8bndx(:,:)
	cur_3bndx(:,jy,:)=0.d0
	x1_8bndx(:,jy,:)=0.d0
	ef_3bndx(:,jy,:)=0.d0

	x_8bndz(:,jy,:)=xint_8bndz(:,:)
	cur_3bndz(:,jy,:)=0.d0
	x1_8bndz(:,jy,:)=0.d0
	ef_3bndz(:,jy,:)=0.d0

	enddo


	if(nrank.eq.0) then
	open(unit=1,file='bx_bndz.dat',status='unknown',form='formatted')
	write(1,2) (bndz_grd(jx,1),bndz_grd(jx,2),bx_bndz(jx),jx=1,nbndz)
    2 format(3(1x,e12.5))
    	close(1)
	endif


	else
	endif

	return
	end
	
