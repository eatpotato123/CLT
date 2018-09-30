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





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!          !!!!  !!!!!!!!!!!!!!!!      !!!!!!!!          !!
!!  !!!!!!!!!!!!  !!!!!!!!!!!!!!   !!!!  !!!!!!!  !!!!!!!!!!
!!  !!!!!!!!!!!!  !!!!!!!!!!!!!!!!  !!!!!!!!!!!!  !!!!!!!!!!
!!        !!!!!!  !!!!!!!!!!!!!!!!!!   !!!!!!!!!        !!!!
!!  !!!!!!!!!!!!  !!!!!!!!!!!!!!!!!!!!   !!!!!!!  !!!!!!!!!!
!!  !!!!!!!!!!!!  !!!!!!!!!!!!!!   !!!   !!!!!!!  !!!!!!!!!!
!!          !!!!           !!!!!!!     !!!!!!!!!          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module for else left subroutines. (not used)
! written by S. WANG
! sortted by H.W. Zhang during 2018 AUG.
! contact changhw@zju.edu.cn if there is any problem.

! intp1(px,ny,pxa,ns,nhf)
! intp2(p2x,p2xa,nhf,nsignx,inv)
! int2spl(p2x,p2xa,nhf,nsignx,inv)
! intt(p2xb,p2xa,nsignx)
! extap(x1,x2,x3,x4)
! ex(d,n)
! sym(d,nsignx)
! stepon_atfs
! stepon_atfs_rk
! artif_sound(n)
! artif_sound_replace(n)
! artif_sound_replace_lax(n)
! artif_sound_replace_rk(n)
! right_atfs(x_as,u,udif)
! artif_sound_implicity(n)
! right_xz_atfs(u,udif)
! tridag_real_period(a,b,c,s,u,n)
! current_driven
! distribution_cd
! distribution_cd_cos
! calculate_cb00
! calculate_ps1
! calculate_ps1_mgax
! recrd_ps1
! recrd_cud
! readin_cud
! find_cud_ox
! pllconduct(ipll)
! pll_subloop(npll)
! pll_subloop_el(npll)
! pll_subloop_rk(npll)
! pll_subloop_y(npll)
! bndry_p
! right_pll(udif)
! right_pll_p0(udif)
! right_pll_p1(udif)
! right_py_p1(udif)
! right_pll_p1_v2(udif)
! pllconduct_implicity
! pll_petsc
! pll_petsc_t1
! pll_smthpline(npll)
! pll_soundwave(npll)
! convtb
! find_opoint_z0
! distribution_cd_oxpoint(ww)
! find_maxmin_1(ww,nky,wmax,wmin,jxlmax,jzlmax,jxlmin,jzlmin)
! find_maxmin(ww,nky,wmax,wmin,jxlmax,jzlmax,jxlmin,jzlmin)
! find_oxpoint_backup(ww)
! find_oxpoint_1st
! find_oxpoint(ww)
! character*3 function cn(n)
! character*3 function cn1(n)
! real*8 function rhom(psival)
! real*8 function rhomp(psival)
! real*8 function rhom(psival)
! real*8 function rhomp(psival)
! integer function sgn(val)
! integer function j_add(val,ary,nn,jj)



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
!ws*****************************************************
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

!ws:dmapb2+++++++++++++++++++++++++++++++++++++++++++++
      subroutine intp1(px,ny,pxa,ns,nhf)
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
      end

!ws****************************************************
      subroutine intp2(p2x,p2xa,nhf,nsignx,inv)
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
      end
      subroutine int2spl(p2x,p2xa,nhf,nsignx,inv)
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
      end
      subroutine intt(p2xb,p2xa,nsignx)
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
      end
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
      subroutine ex(d,n)
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
      end
      subroutine sym(d,nsignx)
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
      end


!!ws:dmapb2----------------------------------------
!ws********************************************
!ws********************************************
	real*8 function rhom(psival)
      use declare
      real*8 psival,rsq
!
!   rsq is the normalized poloidal flux
!
      rsq=(psival-psmin)/(psia-psmin)
      if(iden.eq.1) rhom=(1.00000-alphar*rsq**prho)**arho
!    gaussian density profile
      if(iden.eq.2) rhom=exp(-alphar*rsq)+prho*rsq*(arho-rsq)
! tftr dt shot 66887 electron density profile
      if(iden.eq.3)rhom=(1.3+5.3*(1.-rsq)*(1.-0.95*rsq*(1.-rsq)))/6.6

      if(iden.eq.5)rhom=(1.00000-alphar*rsq-prho*rsq**3)
      return
      end

!ws********************************************
      real*8 function rhomp(psival)
      use declare
      real*8 psival,rsq
!
!   rsq is the normalized poloidal flux
!
      rsq=(psival-psmin)/(psia-psmin)
      if(iden.eq.1) rhomp=-arho*(1.00000-alphar*rsq**prho)**(arho-1.)*alphar*prho*rsq**(prho-1)/(psia-psmin)
!    gaussian density profile
      if(iden.eq.2) rhomp=-alphar*exp(-alphar*rsq)+prho*(arho-2.*rsq)/(psia-psmin)
! tftr dt shot 66887 electron density profile
      if(iden.eq.3) rhomp=(5.3/6.6)*(-(1.-0.95*rsq*(1.-rsq))-0.95*(1.-rsq)*(1.-2.*rsq))/(psia-psmin)

      if(iden.eq.5) rhomp=-alphar-prho*3*rsq**2/(psia-psmin)
      return
      end
   
!ws******************************************************************************
!ws:artificial soundwave 
!ws******************************************************************************
      subroutine stepon_atfs
!
!     this routine time-advances x's bz fourth order in time and second
!     order in space runge-kotta differential scheme.
!     note: x is alwazs the up-to-date value while xm being the
!           intermediate value, and xdif is increment
!
!
      use declare
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
!      call artif_sound_replace_lax(1)
      call artif_sound_replace_rk(1)
!      call artif_sound_implicity(1)

      caf=0.75d0*(0.5+0.5*dtanh((time-40)/5.))
!      call bndry8_x_ex(lbnd)

      return
      end

!ws******************************************************************************
      subroutine stepon_atfs_rk
!
!     this routine time-advances x's bz fourth order in time and second
!     order in space runge-kotta differential scheme.
!     note: x is alwazs the up-to-date value while xm being the
!           intermediate value, and xdif is increment
!
!
      use declare
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
      use declare
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
      use declare
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
	subroutine artif_sound_replace_lax(n)
      use declare
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
	subroutine artif_sound_replace_rk(n)
      use declare
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
      use declare
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
      use declare
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
!ws131117
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
      use declare
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

!ws:***************************************************
	subroutine tridag_real_period(a,b,c,s,u,n)
      include 'mpif.h'
!      complex*16 r(n),u(n)
      integer n,nmax
      real*8 a(n),b(n),c(n),s(n),u(n)
      real*8 a1(n),b1(n),c1(n),s1(n),vu(n),qu(n),pu(n)
!      parameter (nmax=5001)
      integer j
      real*8 bet  !,gam(nmax)

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
      end  

!ws*****************************************************
      integer function sgn(val)
      real*8 val
      if(val .ge. 0) sgn=1
      if(val .lt. 0) sgn=-1
      return
      end

!ws*****************************************************
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

!ws*****************************************************
      subroutine current_driven
      use declare
      include 'mpif.h'
      real*8 tcd,bb2,bb,cud0
!      call calculate_ps1
!      call distribution_cd
!      if(cd_opoint) call distribution_cd_opoint
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

!ws****************************************************
      subroutine distribution_cd
      use declare
      use declare_oxpoint

      real*8 fnv,fnv1
      include 'mpif.h'      
      
      if(lrstrt_cd) then
      open(unit=777,file='cd.dat',status='unknown',form='formatted')
      read(777,3001) 
      read(777,*) cip,fcd,tcds,tcde,tcdd,psmode,delcd,fnv,psshift,br_max0
      close(777)
3001  format('cip,fcd,tcds,tcde,tcdd,psmode,delcd,fnv')
      do jy=iy_first,iy_last
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
!      phase=nmode*yy(jy)+mmode*thxz(jx,jz)+0
!      ps1(jx,jz,jy)
      fn_cdy(jx,jz,jy)=exp(-(psi(jx,jz)+ps1(jx,jz,jy)-(psmode+psshift))**2/delcd**2)/fnv
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

      fnv1=0
!      do jy=iy_first+2,iy_last-2  
      do jz=iz_first+2,iz_last-2
      do jx=ix_first+2,ix_last-2 
      if(psi(jx,jz) .lt. psia) then 
!      fnv1=fnv1+fn_cdy(jx,jz)*(xx(jx+1)-xx(jx-1))*(zz(jz+1)-zz(jz-1))/4
        bb2=x(jx,jz,3,6)**2+x(jx,jz,3,8)**2+x(jx,jz,3,7)**2
        bb=sqrt(bb2)
      fnv1=fnv1+fn_cdy(jx,jz,3)*x(jx,jz,3,7)/bb*(xx(jx+1)-xx(jx-1))*(zz(jz+1)-zz(jz-1))/4
      endif
      enddo
      enddo
!      enddo

      call mpi_allreduce(fnv1,fnv,1,mpi_double_precision,mpi_sum, &
                        mpi_comm_world,ierror)

      fn_cdy(:,:,:)=fn_cdy(:,:,:)/fnv
      if(nrank==0 .and. time .le. timein+1.e-6) then
      
      write(*,*) '***current drive:cd=fcd*ip*fn,ing(fn*ds)=1***'
      write(*,*) 'ip    =',cip
      write(*,*) 'fcd   =',fcd
      write(*,*) 'tcds  =',tcds
      write(*,*) 'tcde  =',tcde
      write(*,*) 'tcdd  =',tcdd
      write(*,*) 'psmode=',psmode
      write(*,*) 'delcd =',delcd
      write(*,*) 'fnv   =',fnv
      write(*,*) 'pssft =',psshift
      write(*,*) '***current drive***'
      open(unit=777,file='cd.dat',status='unknown',form='formatted')
      write(777,3000) 
      write(777,*) cip,fcd,tcds,tcde,tcdd,psmode,delcd,fnv,psshift,br_max0
3000  format('cip,fcd,tcds,tcde,tcdd,psmode,delcd,fnv')
      close(777)
      endif

      return
      end

!ws*****************************************************
      subroutine distribution_cd_cos
      use declare
      real*8 fnv,fnv1
      include 'mpif.h'      
      
      if(lrstrt_cd) then
      open(unit=777,file='cd.dat',status='unknown',form='formatted')
      read(777,3001) 
      read(777,*) cip,fcd,tcds,tcde,tcdd,psmode,delcd,fnv,psshift
      close(777)
3001  format('cip,fcd,tcds,tcde,tcdd,psmode,delcd,fnv')
      do jy=iy_first,iy_last
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
!      phase=nmode*yy(jy)+mmode*thxz(jx,jz)+0
!      ps1(jx,jz,jy)
      fn_cdy(jx,jz,jy)=exp(-(psi(jx,jz)+ps1(jx,jz,jy)-(psmode+psshift))**2/delcd**2)*(1+dcos(nmode*yy(jy)+mmode*thxz(jx,jz)))/fnv
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

      fnv1=0
!      do jy=iy_first+2,iy_last-2  
      do jz=iz_first+2,iz_last-2
      do jx=ix_first+2,ix_last-2 
      if(psi(jx,jz) .lt. psia) then 
!      fnv1=fnv1+fn_cdy(jx,jz)*(xx(jx+1)-xx(jx-1))*(zz(jz+1)-zz(jz-1))/4
        bb2=x(jx,jz,3,6)**2+x(jx,jz,3,8)**2+x(jx,jz,3,7)**2
        bb=sqrt(bb2)
      fnv1=fnv1+fn_cdy(jx,jz,3)*x(jx,jz,3,7)/bb*(xx(jx+1)-xx(jx-1))*(zz(jz+1)-zz(jz-1))/4
      endif
      enddo
      enddo
!      enddo

      call mpi_allreduce(fnv1,fnv,1,mpi_double_precision,mpi_sum, &
                        mpi_comm_world,ierror)

      fn_cdy(:,:,:)=fn_cdy(:,:,:)/fnv
      if(nrank==0 .and. time .le. timein+1.e-6) then
      write(*,*) '***current drive:cd=fcd*ip*fn,ing(fn*ds)=1***'
      write(*,*) 'ip    =',cip
      write(*,*) 'fcd   =',fcd
      write(*,*) 'tcds  =',tcds
      write(*,*) 'tcde  =',tcde
      write(*,*) 'tcdd  =',tcdd
      write(*,*) 'psmode=',psmode
      write(*,*) 'delcd =',delcd
      write(*,*) 'fnv   =',fnv
      write(*,*) 'pssft =',psshift
      write(*,*) '***current drive***'
      open(unit=777,file='cd.dat',status='unknown',form='formatted')
      write(777,3000) 
      write(777,*) cip,fcd,tcds,tcde,tcdd,psmode,delcd,fnv,psshift
3000  format('cip,fcd,tcds,tcde,tcdd,psmode,delcd,fnv')
      close(777)
      endif

      return
      end

!ws*****************************************************
      subroutine calculate_cb00
      use declare 
      include 'mpif.h'
      real*8 fnv,fnv1
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

      fnv1=0
!      do jy=iy_first+2,iy_last-2  
      do jz=iz_first+2,iz_last-2
      do jx=ix_first+2,ix_last-2 
      if(psi(jx,jz) .lt. psia) then 
      fnv1=fnv1+fn_bs(jx,jz)*(xx(jx+1)-xx(jx-1))*(zz(jz+1)-zz(jz-1))/4
      endif
      enddo
      enddo
!      enddo

      call mpi_allreduce(fnv1,fnv,1,mpi_double_precision,mpi_sum, &
                        mpi_comm_world,ierror)

      cb00=fbs*cip/fnv
      return
      end

!ws********************************************
      subroutine calculate_ps1
      use declare
      include 'mpif.h'

      real*8,dimension(mz,my) :: wsx1s,wsx1r
!      integer status(mpi_status_size)
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
      call mpi_send( wsx1s,mz*my, mpi_double_precision,nrank+1, 0,  &
		      mpi_comm_world,ierror )
      endif

      do ii=1,nprx-1
      if(nrkx(nrank).eq. ii) then
      call mpi_recv( wsx1r,mz*my, mpi_double_precision,nrank-1, ii-1,  &
		      mpi_comm_world, status,ierror )
      do jx=ix_first+2,ix_last-2
      ps1(jx,:,:)=ps1(jx,:,:)+wsx1r(:,:)
      enddo

      wsx1s(:,:)=ps1(ix_last-1,:,:)
      if(nrkx(nrank).lt. nprx-1) then
      call mpi_send( wsx1s,mz*my, mpi_double_precision,nrank+1, ii,  &
		      mpi_comm_world,ierror )
      endif
      endif

      enddo


      return
      end

!ws********************************************
      subroutine calculate_ps1_mgax
      use declare
      include 'mpif.h'
      real*8,dimension(mx,2,my) :: wsz2,wsz2s,wsz2r
      real*8,dimension(2,mz,my) :: wsx2,wsx2s,wsx2r
!      integer status(mpi_status_size)
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

!      call mpi_barrier(mpi_comm_world,ierror)

      do kz=0,nprz/2-2
!      nranksend_zp=nrank_mgp+kz*nprx
!      nrankrecv_zp=nrank_mgp+(kz+1)*nprx
!      nranksend_zm=nrank_mgm-kz*nprx
!      nrankrecv_zm=nrank_mgm-(kz+1)*nprx
      
      if(nrank==nrank_mgp+kz*nprx) then
!     write(*,*) 'ws11p',nrank,kz
      wsz2(:,:,:)=ps1(:,iz_last-3:iz_last-2,1:my)
      call mpi_send( wsz2,2*mx*my, mpi_double_precision,nrank+nprx, kz,  &
		      mpi_comm_world,ierror )
!        write(*,*) 'ws11ps',nrank,kz
      endif
!      call mpi_sendrecv(wsz2s,mx*2*my, mpi_double_precision, nrankrecv_zp, kz,  &
!                        wsz2r,mx*2*my, mpi_double_precision, nranksend_zp, kz,  &
!		               mpi_comm_world,status,ierror )

      if(nrank==nrank_mgm-kz*nprx) then
!      write(*,*) 'ws11m',nrank,kz
      wsz2(:,:,:)=ps1(:,iz_first+2:iz_first+3,1:my)
      call mpi_send( wsz2,2*mx*my, mpi_double_precision,nrank-nprx, kz,  &
		      mpi_comm_world,ierror )
!      write(*,*) 'ws11ms',nrank,kz
      endif
!      call mpi_sendrecv(wsz2s,mx*2*my, mpi_double_precision, nrankrecv_zm, kz,  &
!                        wsz2r,mx*2*my, mpi_double_precision, nranksend_zm, kz,  &
!		               mpi_comm_world,status,ierror )
      if(nrank==nrank_mgp+(kz+1)*nprx) then
!      write(*,*) 'ws42',nrank,kz
      call mpi_recv( wsz2 ,2*mx*my, mpi_double_precision,nrank-nprx, kz,  &
		      mpi_comm_world, status,ierror )
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
      call mpi_recv(wsz2,mx*2*my, mpi_double_precision,nrank+nprx, kz,  &
		      mpi_comm_world, status,ierror )
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
 !     call mpi_barrier(mpi_comm_world,ierror)
      enddo

 !     call mpi_barrier(mpi_comm_world,ierror)

      if(nrkx(nrank)==idmgx) then
!      wsx2s(:,:,:)=ps1(ix_last-3:ix_last-2,:,1:my)
      call mpi_send( ps1(ix_last-3:ix_last-2,:,1:my),2*mz*my, mpi_double_precision,nrank+1, 1,  &
		      mpi_comm_world,ierror )
      call mpi_send( ps1(ix_first+2:ix_first+3,:,1:my),2*mz*my, mpi_double_precision,nrank-1, 1,  &
		      mpi_comm_world,ierror )
      endif

      do kx=1,idmgx
!      if(nrank==nrank_mgp+kx) wsx2s(:,:,:)=ps1(ix_last-3:ix_last-2,:,1:my)
!      call mpi_sendrecv(wsx2s,2*mz*my, mpi_double_precision, nrank_mgp+kx+1, kx,  &
!                        wsx2r,2*mz*my, mpi_double_precision, nrank_mgp+kx, kx,  &
!		               mpi_comm_world,status,ierror )

      if(nrkx(nrank)==idmgx+kx) then
          call mpi_recv( wsx2r,2*mz*my, mpi_double_precision,nrank-1, kx,  &
		          mpi_comm_world, status,ierror )
          ps1(ix_first:ix_first+1,:,1:my)=wsx2r(:,:,:)
          do jx=ix_first,ix_last-2
          ps1(jx+1,:,:)=ps1(jx,:,:)-(xx(jx)*x1(jx,:,:,8)+xx(jx+1)*x1(jx+1,:,:,8))/2*(xx(jx+1)-xx(jx))
          enddo
 !         write(*,*) 'ws2',nrank,kx
          if(nrkx(nrank).lt. nprx-1) then 
          call mpi_send( ps1(ix_last-3:ix_last-2,:,1:my),2*mz*my, mpi_double_precision,nrank+1, kx+1,  &
		          mpi_comm_world,ierror )
          endif
      endif

      if(nrkx(nrank)==idmgx-kx) then
          call mpi_recv( wsx2r,2*mz*my, mpi_double_precision,nrank+1, kx,  &
		          mpi_comm_world, status,ierror )
         ps1(ix_last-1:ix_last,:,1:my)=wsx2r(:,:,:)
          do jx=ix_last,ix_first+2,-1
          ps1(jx-1,:,:)=ps1(jx,:,:)-(xx(jx)*x1(jx,:,:,8)+xx(jx-1)*x1(jx-1,:,:,8))/2*(xx(jx-1)-xx(jx))
          enddo
 !         write(*,*) 'ws3',nrank,kx
          if(nrkx(nrank).gt. 0) then 
          call mpi_send( ps1(ix_first+2:ix_first+3,:,1:my),2*mz*my, mpi_double_precision,nrank-1, kx+1,  &
		          mpi_comm_world,ierror )
          endif
      endif

      enddo

      call smthxzy_dis_v2(ps1,1,3)

      return
      end

!ws************************************************************************
      subroutine recrd_ps1
      use declare
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
      use declare
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

!ws*****************************************************	
      subroutine readin_cud
      use declare
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

!ws****************************************************	
      subroutine find_cud_ox
      use declare
      use declare_oxpoint      
          call find_maxmin(cud(:,:,3,2),0,wmax,wmin,jxto,jzto,jxtx,jztx)
      yy_o= 0 
      ps_o= pst(jxto,jzto)
      rr_o= rrt(jxto,jzto) 
      tc_o= 0 

      yy_x= 0 !yy_o
      ps_x= psmode
      rr_x= rrmode
      tc_x= pi/2 
      return
      end
      
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
      !call pll_petsc
      case(5)
     ! call pll_petsc_t1
!      case default
      end select

      return
      end

!ws*****************************************************
      subroutine pll_subloop(npll)
      use declare
      integer npll
      select case(lscheme)
      case(1) !euler
      call pll_subloop_el(npll)
      case(2) !rk4
      call pll_subloop_rk(npll)
      case(4)
      call pllconduct_implicity
!      case default
      end select
      return
      end

!ws*****************************************************
      subroutine pll_subloop_el(npll)
      use declare
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
      subroutine pll_subloop_rk(npll)
      use declare
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
      use declare
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
      use declare
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
      use declare
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
!  d1xf2= d rf / dx  with second-order accuracy central difference
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
      use declare
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
!  d1xf2= d rf / dx  with second-order accuracy central difference
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
      use declare
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
!  d1xf2= d rf / dx  with second-order accuracy central difference
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
      use declare
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
!  d1xf2= d rf / dx  with second-order accuracy central difference
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

!ws*****************************************************
      subroutine right_pll_p1_v2(udif)
!      use declare
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
      end

!ws*****************************************************
      subroutine pllconduct_implicity
      use declare
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
      subroutine pll_petsc
      !use declare
      !call convtb
      !call set_a_s(dt)
      !call possion_solver_3d(x(:,:,:,2))
      !call mpi_transfersm(x(:,:,:,2),1)
      !return
      !end
!ws*****************************************************
      !subroutine pll_petsc_t1
      !use declare
      !call convtb
      !call set_a_s(dt)
      !call possion_solver_3d_t1(x(:,:,:,2))
      !call mpi_transfersm(x(:,:,:,2),1)
      !return
      end

!ws************************************************************************
      subroutine pll_smthpline(npll)
      use declare
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
      use declare
      integer npll

      select case(lscheme)
      case(1)
      call artif_sound_replace(npll)
      case(2)
      call artif_sound_replace_rk(npll)
      case(3)
      call artif_sound_replace_lax(npll)
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
      use declare
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
      d1fbp(fp2,fp1,f0,fm1,a,b,c)= &
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

!       integer status(mpi_status_size)
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
      subroutine find_opoint_z0
      use declare
      use declare_oxpoint
      integer jxtom,jxtop,jxom(1),jxop(1),jxo(1)
      real*8  bzmaxm,bzmaxp,bzmaxo
      real*8, dimension(mx) :: bo
      include 'mpif.h'  
       
      bo(:)=ef(:,jzmode,3,2)
  

      if(nrank==nrank_mode-1) then
      bzmaxm=maxval(bo)
      jxom  =maxloc(bo)
      jxtom =jxom(1)+(nrkx_mode-1)*mxm-2

      call mpi_send( bzmaxm, 1, mpi_double_precision, nrank_mode, 1,  &
		      mpi_comm_world,ierror )
      call mpi_send( jxtom,   1, mpi_integer, nrank_mode, 2,  &
		      mpi_comm_world,ierror )
      endif

      if(nrank==nrank_mode+1) then
      bzmaxp=maxval(bo)
      jxop  =maxloc(bo)
      jxtop =jxop(1)+(nrkx_mode+1)*mxm-2

      call mpi_send( bzmaxp, 1, mpi_double_precision, nrank_mode, 3,  &
		      mpi_comm_world,ierror )
      call mpi_send( jxtop,   1, mpi_integer, nrank_mode, 4,  &
		      mpi_comm_world,ierror )
      endif

      if(nrank==nrank_mode) then
      bzmaxo=maxval(bo)
      jxo   =maxloc(bo)
      jxto  =jxo(1)+(nrkx_mode)*mxm-2

      call mpi_recv( bzmaxm, 1, mpi_double_precision, nrank_mode-1, 1,  &
		      mpi_comm_world, status,ierror  )
      call mpi_recv( jxtom,  1, mpi_integer, nrank_mode-1, 2,  &
		      mpi_comm_world, status,ierror  )
      call mpi_recv( bzmaxp, 1, mpi_double_precision, nrank_mode+1, 3,  &
		      mpi_comm_world, status,ierror  )
      call mpi_recv( jxtop,  1, mpi_integer, nrank_mode+1, 4,  &
		      mpi_comm_world, status,ierror  )
      
      if(bzmaxo.lt.bzmaxm) then
      bzmaxo=bzmaxm
      jxto=jxtom
      endif
      if(bzmaxo.lt.bzmaxp) then
      bzmaxo=bzmaxp
      jxto=jxtop
      endif

      endif
      
      call mpi_bcast(jxto, 1, mpi_integer, nrank_mode, mpi_comm_world,ierror )

      jzto=jztmode

      xx_o=xxt(jxto+1)
      zz_o=zzt(jzto)
      ps_o=pst(jxto+1,jzto)
      rr_o=sqrt((ps_o-psmin)/(psmax-psmin))
      ! th_o=tcht(jxto,jzto)
      if(nrank.eq.0) then
      write(311,700) time,jxto,jzto,xx_o,zz_o,ps_o,rr_o
!700   format((1x,e12.5,2i,4(1x,e12.5)))
700   format(7(1x,e12.5))
      endif

      return
      end

!ws*****************************************************
      subroutine distribution_cd_oxpoint(ww)
      use declare
      use declare_oxpoint
      real*8, dimension(mx,mz) :: ww
      real*8 fnv,fnv1,pscd,phase,cdc,br_old

      include 'mpif.h'      
      
      nrkyo=0    
      jyo=3
      call find_oxpoint(ww)
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
      phase=nmode*(yy(jy)-yy_x)+mmode*(tcxz(jx,jz)-tc_x)
!      pscd=ps_o+(ps_x-ps_o)*(1-dcos(nmode*(yy(jy)-yy_o)+mmode*(tcxz(jx,jz)-tc_o)))
!      fn_cdy(jx,jz,jy)=exp(-(psi(jx,jz)-pscd)**2/delcd**2)*(1.0+dcos(nmode*(yy(jy)-yy_o)+mmode*(tcxz(jx,jz)-tc_o)))
      pscd=ps_x+(ps_o-ps_x)*(1-dcos(phase))/2+psshift
      fn_cdy(jx,jz,jy)=exp(-(psi(jx,jz)-pscd)**2/delcd**2)*(1.0-cdc*dcos(phase))

      enddo
      enddo
      enddo

      fnv1=0
!      do jy=iy_first+2,iy_last-2  
      do jz=iz_first+2,iz_last-2
      do jx=ix_first+2,ix_last-2 
      if(psi(jx,jz) .lt. psia) then 
!      fnv1=fnv1+fn_cdy(jx,jz)*(xx(jx+1)-xx(jx-1))*(zz(jz+1)-zz(jz-1))/4
        bb2=x(jx,jz,3,6)**2+x(jx,jz,3,8)**2+x(jx,jz,3,7)**2
        bb=sqrt(bb2)
      fnv1=fnv1+fn_cdy(jx,jz,3)*x(jx,jz,3,7)/bb*(xx(jx+1)-xx(jx-1))*(zz(jz+1)-zz(jz-1))/4
      endif
      enddo
      enddo
!      enddo

      call mpi_allreduce(fnv1,fnv,1,mpi_double_precision,mpi_sum, &
                        mpi_comm_world,ierror)

      fn_cdy(:,:,:)=fn_cdy(:,:,:)/fnv

      return
      end

!ws**********************************************************************
      subroutine find_maxmin_1(ww,nky,wmax,wmin,jxlmax,jzlmax,jxlmin,jzlmin)
      use declare_grid
      integer nky,jxlmax,jzlmax,jxlmin,jzlmin
      real*8  wmax,wmin
      real*8, dimension(mx,mz) :: ww
      real*8, dimension(6,0:npr-1) :: wm
      include 'mpif.h'  

      call funmax(ww,wm(1,nrank),wm(4,nrank),wm(2,nrank),wm(5,nrank),wm(3,nrank),wm(6,nrank),xx,zz,mx,mz)

!      write(*,600) nrank,wm(:,nrank)
!
      call mpi_allgather(wm(:,nrank),6,mpi_double_precision,wm,6,mpi_double_precision,mpi_comm_world,ierror)
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
      use declare_grid
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
      call mpi_allgather(wm(:,nrank),6,mpi_double_precision,wm,6,mpi_double_precision,mpi_comm_world,ierror)
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
      subroutine find_oxpoint_backup(ww)
      use declare
      use declare_oxpoint

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

      yy_o= yyt(nrkyo*mym+jyo-2)
      ps_o= pst(jxto,jzto)
      rr_o= rrt(jxto,jzto)
      tc_o=tcht(jxto,jzto)      
      tc_o=anint(tc_o/pi*myt)*pi/myt

      yy_x= yy_o
      ps_x= pst(jxtx,jztx)
      rr_x= rrt(jxtx,jztx)
      tc_x=tcht(jxtx,jztx)
      tc_x=anint(tc_x/pi*myt)*pi/myt
!      write(*,*) '*********************opint*********************'
!      write(*,*) nrank,jxto,jzto,ps_o,rr_o,tc_o,yy_o
!
!      write(*,*) '*********************xpint*********************'
!      write(*,*) nrank,jxtx,jztx,ps_x,rr_x,tc_x,yy_x

      if(nrank.eq.0) then
      write(311,1100) time,jxto,jzto,ps_o,rr_o,tc_o,jxtx,jztx,ps_x,rr_x,tc_x
!1100  format((1x,e12.5,2(2i,3(1x,e12.5))))
1100  format(11(1x,e12.5))
      endif

      return
    end
    
!ws**********************************************************************
      subroutine find_oxpoint_1st
      use declare
      use declare_oxpoint

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

      yy_o= 0 
      ps_o= pst(jxto,jzto)
      rr_o= rrt(jxto,jzto) 
      tc_o= 0 

      yy_x= 0 !yy_o
      ps_x= psmode
      rr_x= rrmode
      tc_x= pi/2 !anint(tc_x/pi*myt)*pi/myt
!      write(*,*) '*********************opint*********************'
!      write(*,*) nrank,jxto,jzto,ps_o,rr_o,tc_o,yy_o
!
!      write(*,*) '*********************xpint*********************'
!      write(*,*) nrank,jxtx,jztx,ps_x,rr_x,tc_x,yy_x

      if(nrank.eq.0) then
      write(311,1100) time,jxto,jzto,ps_o,rr_o,tc_o,jxtx,jztx,ps_x,rr_x,tc_x
      
!1100  format((1x,e12.5,2(2i,3(1x,e12.5))))
1100  format(11(1x,e12.5))
      endif

      return
    end
    
!ws**********************************************************************
      subroutine find_oxpoint(ww)
      use declare
      use declare_oxpoint

      real*8  wmax,wmin
      real*8, dimension(mx,mz) :: ww

      call find_maxmin(ww,nrkyo,wmax,wmin,jxto,jzto,jxtx,jztx)

      if(pst(jxto,jzto) .gt. ps_o) then
      ps_o= pst(jxto,jzto)
      rr_o= rrt(jxto,jzto)
      endif

 
      if(nrank.eq.0) then
      write(311,1100) time,jxto,jzto,ps_o,rr_o,tc_o,jxtx,jztx,ps_x,rr_x,tc_x
!1100  format((1x,e12.5,2(2i,3(1x,e12.5))))
1100  format(11(1x,e12.5))
      endif

      return
    end

