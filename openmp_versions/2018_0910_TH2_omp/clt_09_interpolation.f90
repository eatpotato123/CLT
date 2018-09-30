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




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!    !!!!!!!!!          !!!!!!!!      !!!!!!!!!!!!!
!!!!!!!!  !!!!!!!!!!!!!!  !!!!!!!!!!!!  !!!!  !!!!!!!!!!!
!!!!!!!!  !!!!!!!!!!!!!!  !!!!!!!!!!!!  !!!!! !!!!!!!!!!!
!!!!!!!!  !!!!!!!!!!!!!!  !!!!!!!!!!!!  !!!! !!!!!!!!!!!!
!!!!!!!!  !!!!!!!!!!!!!!  !!!!!!!!!!!!      !!!!!!!!!!!!!
!!!!!!!!  !!!!!!!!!!!!!!  !!!!!!!!!!!!  !!!!!!!!!!!!!!!!!
!!!!!!!    !!!!!!!!!!!!!  !!!!!!!!!!!!  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module for MPI interpolation subroutines.
! written by S. WANG
! sortted by H.W. Zhang during 2018 AUG.
! contact changhw@zju.edu.cn if there is any problem.

! including subroutines as follows
! quad(x1,x2,x3,y1,y2,y3,y,ans,ierr,deriv)
! interp1d2l(x1,x2,x3,y1,y2,y3,y,ans)
! interp1d3l(x1,x2,x3,x4,y1,y2,y3,y4,y,ans)
! interp2d3l(x1,x2,z1,z2,w11,w12,w21,w22,xw,zw,ans)
! interp2d4p(x11,x12,x21,x22,z11,z12,z21,z22,w11,w12,w21,w22,xw,zw,ans)
! interp_xz2ps_unif(bxz,xi,zi,nx,nz,bps,xs,zs,ns)
! interp_xz2ps(bxz,xi,zi,nx,nz,bps,xs,zs,ns)
! cubic(x1,x2,x3,x4,y1,y2,y3,y4,y,ans,dans,is,ierr)
! newline(xp,zp,sp,np)
! newlin2(xp,zp,xpn,zpn,sp,spn,np,np2)
! extap(x1,x2,x3,x4)
! extrap_xz2rt(bxz,x0,xi,zi,nx,nz,brt,rh,th,nt)
! extrap2_xz2rt(bxz,x0,xi,zi,nx,nz,brt,brt_dr,brt_dt,rh,th,nt)
! extrap_rt2xz(brt1,rh1,brt2,rh2,th,nt,bxz,ri,thi,nb)
! extrap1d0_r(brt2,rh2,brt3,rh3,nt)
! extrap1d1_r(brt1,rh1,brt2,rh2,brt3,rh3,nt)
! gauss_solve(A,b,x,N)
! uptri(A,b,x,N)
! lag_intp1d4p(x0,f0,x1,f1,x2,f2,x3,f3,x_lag,f_lag)
! shell_sort(n,a,m,b)
! piksrt(n,arr,m,b)



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
      real*8 x1,x2,x3,y1,y2,y3,y,ans
      real*8 d1,d2,d3,tmp_add
!cray  lcm(quad2)
      ans=0.d0
!      tmp_add=x1+x2+x3+y1+y2+y3+y
!      if(isnan(tmp_add)) then 
!        if(nstep.eq.0)  print*,x1,x2,x3,y1,y2,y3,y
!      else
      d1 = (y1-y2)*(y1-y3)
      d2 = (y2-y3)*(y2-y1)
      d3 = (y3-y1)*(y3-y2)
!! hw add if to avoid divide zeor situations      
!      if((x1.gt.-1.d10).and.(x1.lt.1.d10).and. &
!         (x2.gt.-1.d10).and.(x2.lt.1.d10).and. &
!         (x3.gt.-1.d10).and.(x3.lt.1.d10).and. &
!         (y1.gt.-1.d10).and.(y1.lt.1.d10).and. &
!         (y2.gt.-1.d10).and.(y2.lt.1.d10).and. &
!         (y3.gt.-1.d10).and.(y3.lt.1.d10).and. &
!         (y.gt.-1.d10).and.(y.lt.1.d10)) then
      if(dabs(d1*d2*d3).ge.0.d0) then
      ans = x1*(y-y2)*(y-y3)/d1 &
          + x2*(y-y3)*(y-y1)/d2 &
          + x3*(y-y1)*(y-y2)/d3
      endif
!      endif

!      endif

      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     * * * number 3j * * *                                c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine interp1d3l(x1,x2,x3,x4,y1,y2,y3,y4,y,ans)
      real*8 x1,x2,x3,x4,y1,y2,y3,y4,y,ans
      real*8 d1,d2,d3,d4
!cray  lcm(cubic)
!c
!c      this routine returns interpolated value 
!c
      d1 = (y1-y2)*(y1-y3)*(y1-y4)
      d2 = (y2-y1)*(y2-y3)*(y2-y4)
      d3 = (y3-y1)*(y3-y2)*(y3-y4)
      d4 = (y4-y1)*(y4-y2)*(y4-y3)
! hw add if to avoid divide zeor situations      
      if(d1*d2*d3*d4.ne.0.d0) then
      ans = x1*(y-y2)*(y-y3)*(y-y4)/d1 &
          + x2*(y-y1)*(y-y3)*(y-y4)/d2 &
          + x3*(y-y1)*(y-y2)*(y-y4)/d3 &
          + x4*(y-y1)*(y-y2)*(y-y3)/d4 
      endif

      return
      end

!ws*************************************************************
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

!ws*************************************************************
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

!ws******************************************************************
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

!ws******************************************************************
      subroutine interp_xz2ps(bxz,xi,zi,nx,nz,bps,xs,zs,ns)
      implicit real*8 (b-h,o-z)
      dimension xi(nx),zi(nz),bxz(nx,nz),xs(ns),zs(ns),bps(ns)
      lx=1
      lz=1
!mpi   ----------------------------------------------------------------
!      call mpi_allgather(ff,m3d,mpi_double_complex,xfc,m3d, &
!            mpi_double_complex,mpi_comm_world,ierror)
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

!ws****************************************************
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

!ws*************************************************************
      subroutine extrap_xz2rt(bxz,x0,xi,zi,nx,nz,brt,rh,th,nt)
      implicit real*8 (b-h,o-z)
      dimension xi(nx),zi(nz),bxz(nx,nz)
      dimension th(nt),brt(nt)
      real*8 rh,xl,zl
      lx=1
      lz=1
!mpi   ----------------------------------------------------------------
!      call mpi_allgather(ff,m3d,mpi_double_complex,xfc,m3d, &
!            mpi_double_complex,mpi_comm_world,ierror)
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

!ws*************************************************************
      subroutine extrap2_xz2rt(bxz,x0,xi,zi,nx,nz,brt,brt_dr,brt_dt,rh,th,nt)
      implicit real*8 (b-h,o-z)
      dimension xi(nx),zi(nz),bxz(nx,nz)
      dimension th(nt),brt(nt),brt_dr(nt),brt_dt(nt)
      real*8 rh,xl,zl,brt_dx,brt_dz

      lx=1
      lz=1
!mpi   ----------------------------------------------------------------
!      call mpi_allgather(ff,m3d,mpi_double_complex,xfc,m3d, &
!            mpi_double_complex,mpi_comm_world,ierror)
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

!ws*************************************************************	
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

!ws*************************************************************	
      subroutine extrap1d0_r(brt2,rh2,brt3,rh3,nt)
      implicit real*8 (b-h,o-z)
      dimension brt2(nt),brt3(nt)
      do jt=1,nt
      brt3(jt)=brt2(jt)
      enddo
      return
      end

!ws*************************************************************	
      subroutine extrap1d1_r(brt1,rh1,brt2,rh2,brt3,rh3,nt)
      implicit real*8 (b-h,o-z)
      dimension brt1(nt),brt2(nt),brt3(nt)
      beta1=(rh3-rh1)/(rh2-rh1)
      do 20 jt=1,nt
      brt3(jt)=brt1(jt)*(1-beta1)+brt2(jt)*beta1
   20 continue
      return
      end

!hw******************************************************************
subroutine gauss_solve(A,b,x,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-8
!-----------------------------------------------------
!  Purpose   :  楂樻柉娑堝幓娉?!                 Ax=b
!-----------------------------------------------------
!  Input  parameters  :
!       1.   A(N,N)绯绘暟鐭╅樀
!       2.   b(N)鍙冲悜閲?!       3.   N鏂圭▼缁存暟
!  Output parameters  :
!       1.  x  鏂圭▼鐨勬牴
!       2.
!  Common parameters  :
!
!----------------------------------------------------

implicit real*8(a-z)

integer::i,k,N

real*8::A(N,N),b(N),x(N)

real*8::Aup(N,N),bup(N)

!Ab涓哄骞跨煩闃? [Ab]
real*8::Ab(N,N+1)

Ab(1:N,1:N)=A

Ab(:,N+1)=b


!-------------------------------
!  杩欐鏄?楂樻柉娑堝幓娉曠殑鏍稿績閮ㄥ垎
do k=1,N-1

   do i=k+1,N
   
     temp=Ab(i,k)/Ab(k,k)
     
     Ab(i,:)=Ab(i,:)-temp*Ab(k,:)
   
   end do

end do

!-----------------------------
! 缁忚繃涓婁竴姝ワ紝Ab宸茬粡鍖栦负濡備笅褰㈠紡鐨勭煩闃?!            | *  *  *  *  # |
!     [A b]= | 0  *  *  *  # |
!            | 0  0  *  *  # |
!            | 0  0  0  *  # |
!
Aup(:,:)=Ab(1:N,1:N)

bup(:)=Ab(:,N+1)

!璋冪敤鐢ㄤ笂涓夎鏂圭▼缁勭殑鍥炲甫鏂规硶
call uptri(Aup,bup,x,n)

end subroutine gauss_solve

!hw**********************************************************************
subroutine uptri(A,b,x,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-8
!-----------------------------------------------------
!  Purpose   :  涓婁笁瑙掓柟绋嬬粍鐨勫洖甯︽柟娉?!                 Ax=b
!-----------------------------------------------------
!  Input  parameters  :
!       1.   A(N,N)绯绘暟鐭╅樀
!       2.   b(N)鍙冲悜閲?!       3.   N鏂圭▼缁存暟
!  Output parameters  :
!       1.  x  鏂圭▼鐨勬牴
!       2.
!  Common parameters  :
!
!----------------------------------------------------

implicit real*8(a-z)

integer::i,j,N

real*8::A(N,N),b(N),x(N)

x(N)=b(N)/A(N,N)

!鍥炲甫閮ㄥ垎
do i=n-1,1,-1
   
    x(i)=b(i)
   do j=i+1,N
    x(i)=x(i)-a(i,j)*x(j)
   end do
    x(i)=x(i)/A(i,i)

end do

end subroutine uptri

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
!hw**********************************************************************
! Lagrange's interpolation in 1 dimensional with 4 points
	subroutine lag_intp1d4p(x0,f0,x1,f1,x2,f2,x3,f3,x_lag,f_lag)
	implicit none
	real*8 x0,f0,x1,f1,x2,f2,x3,f3,x_lag,f_lag
	real*8 l0,l1,l2,l3,dtmp

! hw add if to avoid divided by zero
        if(x2.eq.x1) stop'lag_intp1d4p error'
        dtmp=(x0-x1)*(x0-x2)*(x0-x3)*(x1-x2)*(x1-x3)*(x2-x3)
        if(dtmp.ne.0.d0) then
	l0=(x_lag-x1)*(x_lag-x2)*(x_lag-x3)/(x0-x1)/(x0-x2)/(x0-x3)
	l1=(x_lag-x0)*(x_lag-x2)*(x_lag-x3)/(x1-x0)/(x1-x2)/(x1-x3)
	l2=(x_lag-x0)*(x_lag-x1)*(x_lag-x3)/(x2-x0)/(x2-x1)/(x2-x3)
	l3=(x_lag-x0)*(x_lag-x1)*(x_lag-x2)/(x3-x0)/(x3-x1)/(x3-x2)

	f_lag=f0*l0+f1*l1+f2*l2+f3*l3
        else

!	if((f_lag.gt.-1.d3).and.(f_lag.lt.1.d3)) then
!	else
		  f_lag=f1+(f2-f1)/(x2-x1)*(x_lag-x1)
		  print*,'hahalag2p',x_lag,x1,x2,f1,f2,f_lag
!		  stop
!        endif
	endif
	return
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
!hw**********************************************************************
	subroutine shell_sort(n,a,m,b)
	implicit none
	integer n, m
	real*8 a(n)
	real*8 b(n,m)
!	copy from Numerical recipes in FORTRAN 77, vol. 1 Page323
!	Sorts an array a(1:n) into ascending numerical order by Shell's method (diminishing in-
!	crement sort). n is input; a is replaced on output by its sorted rearrangement.
!	b is the additaion array that need to be sorted in the rank of a with m
!	elements in each row.
	integer i,j,inc
	real*8 v
	real*8, dimension(1,m):: v2
	inc=1
    1 inc=3*inc+1
    	if(inc.le.n) goto 1
    2 continue
    	inc=inc/3
	do i=inc+1,n
	v=a(i)
	v2(1,:)=b(i,:)
	j=i
    3	if(a(j-inc).gt.v) then
		  a(j)=a(j-inc)
		  b(j,:)=b(j-inc,:)
		  j=j-inc
		  if(j.le.inc) goto 4
		  goto 3
	endif
    4 a(j)=v
  	b(j,:)=v2(1,:)
      enddo
	if(inc.gt.1) goto 2
	return
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
!hw**********************************************************************
	subroutine piksrt(n,arr,m,b)
	implicit none
	integer n,m
	real*8 arr(n)
	real*8 b(n,m)
	integer i,j
	real*8 a
	real*8, dimension(1,m):: v
	do j=2,n
	a=arr(j)
	v(1,:)=b(j,:)
	do i=j-1,1,-1
	if(arr(i).le.a) goto 10
	arr(i+1)=arr(i)
	b(i+1,:)=b(i,:)
	enddo
	i=0
   10 arr(i+1)=a
   	b(i+1,:)=v(1,:)
      enddo
	return
	end
	
