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
!!!!!!!!!    !!!!!!!!!!!    !!!!!!!!!!!!  !!!!!!!!!!!!!!!
!!!!!!  !!!!!!!!!!!!!!!  !!  !!!!!!!!!!!  !!!!!!!!!!!!!!!
!!!!   !!!!!!!!!!!!!!!  !!!!  !!!!!!!!!!  !!!!!!!!!!!!!!!
!!!   !!!!!!!!!!!!!!!  !!!!!!  !!!!!!!!!  !!!!!!!!!!!!!!!
!!!!   !!!!!!!!!!!!!            !!!!!!!!  !!!!!!!!!!!!!!!
!!!!!    !!!!!!!!!!  !!!!!!!!!!  !!!!!!!  !!!!!!!!!!!!!!!
!!!!!!!!    !!!!!!  !!!!!!!!!!!!  !!!!!!           !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module for calculate subroutines (iteration).
! old method written by S. WANG.
! new method (cut cell) written by H.W. Zhang during 2018 Mar-Jul.
! sortted by H.W. Zhang during 2018 AUG.
! contact changhw@zju.edu.cn if there is any problem.

! including subroutines as follows
! setdt
! calculate_eta

! stepon                         -|
! right                           |
! convt                           |
! current(kf)                     |
! convtc (not used)               |> old method by S. WANG
! current_boot (not used)         |
! efield                          |  
! convte                          |
! convtdivb (not used)           -|

! right_cut_cell                 -|
! stepon_cut_cell                 |
! convt_cut_cell                  |> cut cell method by H.W. Zhang
! current_cut_cell                |
! efield_cut_cell (not used)      |
! efield_cut_cell_v2_fixed       -|

! stepon_openacc                -|
! right_openacc                  |
! current_openacc                |> openacc method by H.W. Zhang
! convt_openacc                  |
! efield_openacc                -|

!ws************************************************************
      subroutine setdt
      use declare
      include 'mpif.h'
!
      dt1=100.
      
      do 1 jy=iy_first+2,iy_last-2
      do 1 jz=iz_first+2,iz_last-2
      do 1 jx=ix_first+2,ix_last-2
!      if(psi(jx,jz).lt.psia1) then
      if(gdtp_ep(jx,jz,1).ne.4) then
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

!ws********************************************
      subroutine calculate_eta
      use declare
      include 'mpif.h'
!  d1fc= d f / dx  with fourth-order accuracy central difference
      d1fc(fm2,fm1,f0,fp1,fp2,a,b,c,d)= &
       a*(fp1-f0)+b*(f0-fm1)+c*(fp2-f0)+d*(f0-fm2)

      do jy=iy_first,iy_last
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
!      if(psi(jx,jz) .lt. psia) then
      if(gdtp_ep(jx,jz,1).ne.4) then
	if(lrstrt) then
      tm(jx,jz,jy)=xint(jx,jz,2)/xint(jx,jz,1)
	else
      tm(jx,jz,jy)=x(jx,jz,jy,2)/x(jx,jz,jy,1)
	endif
      eta(jx,jz,jy)=eta0*(1.0+(etacut-1.0)*tanh(((tm(jx,jz,jy)/tm00)**(-1.5)-1.0)/(etacut-1.0)))
      endif
      enddo
      enddo
      enddo
      
      do jy=iy_first+2,iy_last-2  
      do jz=iz_first+2,iz_last-2
      do jx=ix_first+2,ix_last-2 
!      if(psi(jx,jz) .lt. psia) then         
      if(gdtp_ep(jx,jz,1).ne.4) then
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

!swang*****************************************************************************
!swang*****************************************************************************
!subroutines for old method
!swang*****************************************************************************
!swang*****************************************************************************
	
!ws******************************************************************************
      subroutine stepon
!
!     this routine time-advances x's bz fourth order in time and second
!     order in space runge-kotta differential scheme.
!     note: x is alwazs the up-to-date value while xm being the
!           intermediate value, and xdif is increment
!
!
      use declare
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
!      if(eta_from_t) call calculate_eta
      
      return
      end

!ws***************************************************************
	subroutine right
      use declare
      include 'mpif.h'
      real*8 drvx_dx,drey_dx,dez_dx,dex_dy,dez_dy,dex_dz,dey_dz
      real*8, dimension(mx,mz,my) :: rvx,rey
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
!  d1xf2= d rf / dx  with second-order accuracy central difference
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
      if(ef_mode) call efield
      do jy=iy_first,iy_last
      do jz=iz_first,iz_last      
      do jx=ix_first,ix_last
      rvx(jx,jz,jy)=xx(jx)*x(jx,jz,jy,3)
      rey(jx,jz,jy)=xx(jx)*ef(jx,jz,jy,2)
      enddo
      enddo
      enddo

      do 1 jy=iy_first+2,iy_last-2
      do 1 jz=iz_first+2,iz_last-2      
      do 1 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.psia1) then

        drvx_dx=d1fc(rvx(jx-2,jz,jy),rvx(jx-1,jz,jy),rvx(jx,jz,jy),rvx(jx+1,jz,jy),rvx(jx+2,jz,jy),ax1(jx),bx1(jx),cx1(jx),dx1(jx))

!clt option for linear mhd simulation!!!!
      if(linear_mhd) then      

        xdif(jx,jz,jy,1)=-x1(jx,jz,jy,3)*xint_dx(jx,jz,1) &
            -xint(jx,jz,1)*drvx_dx/xx(jx) &
            -(xy(jx,jz,jy,1)*xint(jx,jz,4)+xint(jx,jz,1)*xy(jx,jz,jy,4))/xx(jx) &
            -x1(jx,jz,jy,5)*xint_dz(jx,jz,1)-xint(jx,jz,1)*x1z(jx,jz,jy,5) &
    !ws:for poloidal flow
            -xint(jx,jz,3)*x1r(jx,jz,jy,1)-x1(jx,jz,jy,1)*(xint(jx,jz,3)/xx(jx)+xint_dx(jx,jz,3)) &
            -xint(jx,jz,5)*x1z(jx,jz,jy,1)-x1(jx,jz,jy,1)*xint_dz(jx,jz,5)

            xdif(jx,jz,jy,2)=-x1(jx,jz,jy,3)*xint_dx(jx,jz,2) &
            -gamma*xint(jx,jz,2)*(drvx_dx/xx(jx)) &
            -(xy(jx,jz,jy,2)*xint(jx,jz,4)+gamma*xint(jx,jz,2)*xy(jx,jz,jy,4))/xx(jx) &
            -x1(jx,jz,jy,5)*xint_dz(jx,jz,2)-gamma*xint(jx,jz,2)*x1z(jx,jz,jy,5) &    
    !ws:for poloidal flow
            -xint(jx,jz,3)*x1r(jx,jz,jy,2)-gamma*x1(jx,jz,jy,2)*(xint(jx,jz,3)/xx(jx)+xint_dx(jx,jz,3)) &
            -xint(jx,jz,5)*x1z(jx,jz,jy,2)-gamma*x1(jx,jz,jy,2)*xint_dz(jx,jz,5)       
     
             
            xdif(jx,jz,jy,3)=-xint(jx,jz,3)*x1r(jx,jz,jy,3)-xint(jx,jz,5)*x1z(jx,jz,jy,3) &
            -xint(jx,jz,4)*xy(jx,jz,jy,3)/xx(jx)+xint(jx,jz,4)*x1(jx,jz,jy,4)/xx(jx) &
            +((cur(jx,jz,jy,2))*xint(jx,jz,8)+(cint(jx,jz,2))*x1(jx,jz,jy,8) &
            - (cur(jx,jz,jy,3))*xint(jx,jz,7)-(cint(jx,jz,3))*x1(jx,jz,jy,7) &
            -x1r(jx,jz,jy,2))/xint(jx,jz,1) &
            -x1(jx,jz,jy,3)*xint_dx(jx,jz,3)-x1(jx,jz,jy,5)*xint_dz(jx,jz,3) &
            +x1(jx,jz,jy,4)*xint(jx,jz,4)/xx(jx) &
            +(-xint(jx,jz,3)*xint_dx(jx,jz,3)-xint(jx,jz,5)*xint_dz(jx,jz,3) &
            +xint(jx,jz,4)*xint(jx,jz,4)/xx(jx))*x1(jx,jz,jy,1)/xint(jx,jz,1)


            xdif(jx,jz,jy,4)=-xint(jx,jz,3)*x1r(jx,jz,jy,4)-xint(jx,jz,5)*x1z(jx,jz,jy,4) &
            -xint(jx,jz,4)*xy(jx,jz,jy,4)/xx(jx)-xint(jx,jz,4)*x1(jx,jz,jy,3)/xx(jx) &
            +((cur(jx,jz,jy,3))*xint(jx,jz,6)+(cint(jx,jz,3))*x1(jx,jz,jy,6) &
            - (cur(jx,jz,jy,1))*xint(jx,jz,8)-(cint(jx,jz,1))*x1(jx,jz,jy,8) &   
            -xy(jx,jz,jy,2)/xx(jx))/xint(jx,jz,1) &
            -x1(jx,jz,jy,3)*xint_dx(jx,jz,4)-x1(jx,jz,jy,5)*xint_dz(jx,jz,4) &
            -x1(jx,jz,jy,4)*xint(jx,jz,3)/xx(jx) &
            +(-xint(jx,jz,3)*xint_dx(jx,jz,4)-xint(jx,jz,5)*xint_dz(jx,jz,4) &
            -xint(jx,jz,4)*xint(jx,jz,3)/xx(jx))*x1(jx,jz,jy,1)/xint(jx,jz,1)
     
       
            xdif(jx,jz,jy,5)=-xint(jx,jz,3)*x1r(jx,jz,jy,5)-xint(jx,jz,5)*x1z(jx,jz,jy,5) &
            -xint(jx,jz,4)*xy(jx,jz,jy,5)/xx(jx) &     
            +((cur(jx,jz,jy,1))*xint(jx,jz,7)+(cint(jx,jz,1))*x1(jx,jz,jy,7) &
            - (cur(jx,jz,jy,2))*xint(jx,jz,6)-(cint(jx,jz,2))*x1(jx,jz,jy,6) & 
            -x1z(jx,jz,jy,2))/xint(jx,jz,1) &
            -x1(jx,jz,jy,3)*xint_dx(jx,jz,5)-x1(jx,jz,jy,5)*xint_dz(jx,jz,5) &
            +(-xint(jx,jz,3)*xint_dx(jx,jz,5)-xint(jx,jz,5)*xint_dz(jx,jz,5)) &
            *x1(jx,jz,jy,1)/xint(jx,jz,1)


	else

      xdif(jx,jz,jy,1)=-x1(jx,jz,jy,3)*xr(jx,jz,jy,1) &
       -x(jx,jz,jy,1)*drvx_dx/xx(jx) &
!rplc       -x(jx,jz,jy,3)*x(jx,jz,jy,1)/xx(jx)-x(jx,jz,jy,1)*xr(jx,jz,jy,3) &
       -(xy(jx,jz,jy,1)*x(jx,jz,jy,4)+x(jx,jz,jy,1)*xy(jx,jz,jy,4))/xx(jx) &
       -x1(jx,jz,jy,5)*xz(jx,jz,jy,1)-x(jx,jz,jy,1)*x1z(jx,jz,jy,5) &
!ws:for poloidal flow
       -xint(jx,jz,3)*x1r(jx,jz,jy,1)-x1(jx,jz,jy,1)*(xint(jx,jz,3)/xx(jx)+xint_dx(jx,jz,3)) &
       -xint(jx,jz,5)*x1z(jx,jz,jy,1)-x1(jx,jz,jy,1)*xint_dz(jx,jz,5)

      xdif(jx,jz,jy,2)=-x1(jx,jz,jy,3)*xr(jx,jz,jy,2) &
       -gamma*x(jx,jz,jy,2)*(drvx_dx/xx(jx)) &
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

	endif
          
     
!      xdif(jx,jz,jy,6)=-x(jx,jz,jy,3)*xr(jx,jz,jy,6)-x(jx,jz,jy,6)*x(jx,jz,jy,3)/xx(jx) &
!       -(x(jx,jz,jy,4)*xy(jx,jz,jy,6)-x(jx,jz,jy,7)*xy(jx,jz,jy,3)+x(jx,jz,jy,6)*xy(jx,jz,jy,4))/xx(jx) &
!       -x(jx,jz,jy,5)*xz(jx,jz,jy,6)+x(jx,jz,jy,8)*xz(jx,jz,jy,3)-x(jx,jz,jy,6)*xz(jx,jz,jy,5) !&
!!j x grd_eta
!!       -(cur(jx,jz,jy,2)-cj(jx,jz))*etaz(jx,jz,jy)
!     
!      xdif(jx,jz,jy,7)=-x(jx,jz,jy,3)*xr(jx,jz,jy,7)+x(jx,jz,jy,6)*xr(jx,jz,jy,4) &
!       -x(jx,jz,jy,6)*x(jx,jz,jy,4)/xx(jx)-x(jx,jz,jy,7)*xr(jx,jz,jy,3) &
!       -xy(jx,jz,jy,7)*x(jx,jz,jy,4)/xx(jx) &
!       -xz(jx,jz,jy,7)*x(jx,jz,jy,5)+x(jx,jz,jy,8)*xz(jx,jz,jy,4)-x(jx,jz,jy,7)*xz(jx,jz,jy,5) !&
!!j x grd_eta
!!       -cur(jx,jz,jy,3)*etax(jx,jz,jy)+cur(jx,jz,jy,1)*etaz(jx,jz,jy)
!      
!      xdif(jx,jz,jy,8)=-x(jx,jz,jy,3)*xr(jx,jz,jy,8)-x(jx,jz,jy,3)*x(jx,jz,jy,8)/xx(jx) &
!       -x(jx,jz,jy,8)*xr(jx,jz,jy,3)+x(jx,jz,jy,6)*xr(jx,jz,jy,5) &
!       -(xy(jx,jz,jy,8)*x(jx,jz,jy,4)-xy(jx,jz,jy,5)*x(jx,jz,jy,7)+xy(jx,jz,jy,4)*x(jx,jz,jy,8))/xx(jx) &
!       -xz(jx,jz,jy,8)*x(jx,jz,jy,5) !&
!!j x grd_eta
!!       +(cur(jx,jz,jy,2)-cj(jx,jz))*etax(jx,jz,jy)

!      xdif(jx,jz,jy,6)=-efy(jx,jz,jy,3)/xx(jx)+efz(jx,jz,jy,2) !&
!!                      +cdb0*divb_x(jx,jz,jy)
!    
!      xdif(jx,jz,jy,7)=-efz(jx,jz,jy,1)+efx(jx,jz,jy,3) !&
!!                      +cdb0*divb_y(jx,jz,jy)/xx(jx)
!
!      xdif(jx,jz,jy,8)=efy(jx,jz,jy,1)/xx(jx)-efx(jx,jz,jy,2)-ef(jx,jz,jy,2)/xx(jx) !&
!!                      +cdb0*divb_z(jx,jz,jy)

      if(ef_mode) then
      !drey_dx=d1xf2(ef(jx-1,jz,jy,2),ef(jx,jz,jy,2),ef(jx+1,jz,jy,2),xx(jx-1),xx(jx),xx(jx+1))
      !dez_dx =d1f2(ef(jx-1,jz,jy,3),ef(jx,jz,jy,3),ef(jx+1,jz,jy,3),xx(jx-1),xx(jx),xx(jx+1))
      !dex_dz =d1f2(ef(jx,jz-1,jy,1),ef(jx,jz,jy,1),ef(jx,jz+1,jy,1),zz(jz-1),zz(jz),zz(jz+1))
      !dey_dz =d1f2(ef(jx,jz-1,jy,2),ef(jx,jz,jy,2),ef(jx,jz+1,jy,2),zz(jz-1),zz(jz),zz(jz+1))
      !dex_dy =d1f2(ef(jx,jz,jy-1,1),ef(jx,jz,jy,1),ef(jx,jz,jy+1,1),yy(jy-1),yy(jy),yy(jy+1))
      !dez_dy =d1f2(ef(jx,jz,jy-1,3),ef(jx,jz,jy,3),ef(jx,jz,jy+1,3),yy(jy-1),yy(jy),yy(jy+1))

     drey_dx=d1fc(rey(jx-2,jz,jy),rey(jx-1,jz,jy),rey(jx,jz,jy) &
         ,rey(jx+1,jz,jy),rey(jx+2,jz,jy),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
     dez_dx =d1fc(ef(jx-2,jz,jy,3),ef(jx-1,jz,jy,3),ef(jx,jz,jy,3) &
         ,ef(jx+1,jz,jy,3),ef(jx+2,jz,jy,3),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
     dex_dz =d1fc(ef(jx,jz-2,jy,1),ef(jx,jz-1,jy,1),ef(jx,jz,jy,1) &
         ,ef(jx,jz+1,jy,1),ef(jx,jz+2,jy,1),az1(jz),bz1(jz),cz1(jz),dz1(jz))
     dey_dz =d1fc(ef(jx,jz-2,jy,2),ef(jx,jz-1,jy,2),ef(jx,jz,jy,2) &
         ,ef(jx,jz+1,jy,2),ef(jx,jz+2,jy,2),az1(jz),bz1(jz),cz1(jz),dz1(jz))
     dex_dy =d1fc(ef(jx,jz,jy-2,1),ef(jx,jz,jy-1,1),ef(jx,jz,jy,1) &
         ,ef(jx,jz,jy+1,1),ef(jx,jz,jy+2,1),ay1(jy),by1(jy),cy1(jy),dy1(jy)) 
     dez_dy =d1fc(ef(jx,jz,jy-2,3),ef(jx,jz,jy-1,3),ef(jx,jz,jy,3) &
         ,ef(jx,jz,jy+1,3),ef(jx,jz,jy+2,3),ay1(jy),by1(jy),cy1(jy),dy1(jy)) 
          
      xdif(jx,jz,jy,6)=-dez_dy/xx(jx)+dey_dz !&
!                       +eta(jx,jz,jy)*cint_dz(jx,jz,2)  
      xdif(jx,jz,jy,7)=-dex_dz+dez_dx !&
!                       +eta(jx,jz,jy)*(-cint_dz(jx,jz,1)+cint_dx(jx,jz,3))
      xdif(jx,jz,jy,8)=(dex_dy-drey_dx)/xx(jx) !&
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

      xdif(jx,jz,jy,8)=(x(jx,jz,jy,5)*drvx_dx/xx(jx)+x(jx,jz,jy,6)*xr(jx,jz,jy,5) &
                       -x(jx,jz,jy,3)*xr(jx,jz,jy,8)-x(jx,jz,jy,8)*drvx_dx/xx(jx)) &    
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
      if(ef_mode .and. smoothbout) then
!      if(.not.implicitb) then
      do 41 jy=iy_first+2,iy_last-2
      do 41 jz=iz_first+2,iz_last-2
      do 41 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.psia1) then
!!ws:db/dt=...+eta*grd2 b
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
      if(.not.etaj_in_e) then
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

!!ws:db/dt=...-eta*grd x j
!      xdif(jx,jz,jy,6)=xdif(jx,jz,jy,6)+cur(jx,jz,jy,2)*etaz(jx,jz,jy) &
!       -eta(jx,jz,jy)*(cuy(jx,jz,jy,3)/xx(jx)-cuz(jx,jz,jy,2))
!     
!      xdif(jx,jz,jy,7)=xdif(jx,jz,jy,7)+cur(jx,jz,jy,3)*etax(jx,jz,jy)-cur(jx,jz,jy,1)*etaz(jx,jz,jy) &
!       -eta(jx,jz,jy)*(cuz(jx,jz,jy,1)-cux(jx,jz,jy,3))
!            
!      xdif(jx,jz,jy,8)=xdif(jx,jz,jy,8)-cur(jx,jz,jy,2)*etax(jx,jz,jy) &
!       -eta(jx,jz,jy)*(cux(jx,jz,jy,2)+cur(jx,jz,jy,2)/xx(jx)-cuy(jx,jz,jy,1)/xx(jx))
!!ws:db/dt=...+eta*grd2 b
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
      use declare
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
      xy(jx,jz,jy,m) =d1fc(x1(jx,jz,jy-2,m),x1(jx,jz,jy-1,m),x1(jx,jz,jy,m),x1(jx,jz,jy+1,m),x1(jx,jz,jy+2,m) &
                ,ay1(jy),by1(jy),cy1(jy),dy1(jy))      
      xy2(jx,jz,jy,m)=d2fc(x1(jx,jz,jy-2,m),x1(jx,jz,jy-1,m),x1(jx,jz,jy,m),x1(jx,jz,jy+1,m),x1(jx,jz,jy+2,m) &
                ,ay2(jy),by2(jy),cy2(jy),dy2(jy))
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
      use declare
      include 'mpif.h'
      real*8 drby_dx
      real*8, dimension(mx,mz,my) :: rby
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
!  d1xf2= d rf / dx  with second-order accuracy central difference
      d1xf2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(xp1*fp1-x0*f0) &
         -(xp1-x0)/(xm1-x0)*(xm1*fm1-x0*f0))/(xm1-xp1)

!       integer status(mpi_status_size)

!      do 10 m=1,3
      do 10 jy=iy_first,iy_last
      do 10 jz=iz_first,iz_last
      do 10 jx=ix_first,ix_last
      rby(jx,jz,jy)=xx(jx)*x1(jx,jz,jy,7)
   10 continue

      do 1 jy=iy_first+2,iy_last-2
      do 1 jz=iz_first+2,iz_last-2
      do 1 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.psia1) then
      cur(jx,jz,jy,1)=xy(jx,jz,jy,8)/xx(jx)-x1z(jx,jz,jy,7) 
      cur(jx,jz,jy,2)=x1z(jx,jz,jy,6)-x1r(jx,jz,jy,8)
!      cur(jx,jz,jy,3)=x1r(jx,jz,jy,7)+(x1(jx,jz,jy,7)-xy(jx,jz,jy,6))/xx(jx)
      drby_dx=d1fc(rby(jx-2,jz,jy),rby(jx-1,jz,jy),rby(jx,jz,jy) &
         ,rby(jx+1,jz,jy),rby(jx+2,jz,jy),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
  !    drby_dx=d1xf2(x1(jx-1,jz,jy,7),x1(jx,jz,jy,7),x1(jx+1,jz,jy,7),xx(jx-1),xx(jx),xx(jx+1))
      cur(jx,jz,jy,3)=(drby_dx-xy(jx,jz,jy,6))/xx(jx)

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

!       integer status(mpi_status_size)
!      init1=1
!      init2=1
!       jxs = ix_first
!      jxe = ix_last
!      jzs = iz_first
!      jze = iz_last
!      if (nrank == 0)        jxs=ix_first
!      if (nrank == nsize - 1) jxe =ix_last
      if(spectral) then   

      
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
      use declare
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
      subroutine efield
      use declare
      include 'mpif.h'


      do 1 jy=iy_first,iy_last
      do 1 jz=iz_first,iz_last
      do 1 jx=ix_first,ix_last
!      ef(jx,jz,jy,1)=-x(jx,jz,jy,4)*x(jx,jz,jy,8)+x(jx,jz,jy,5)*x(jx,jz,jy,7)+eta(jx,jz,jy)*cur(jx,jz,jy,1)+xint(jx,jz,4)*xint(jx,jz,8)-xint(jx,jz,5)*xint(jx,jz,7)
!      ef(jx,jz,jy,2)=-x(jx,jz,jy,5)*x(jx,jz,jy,6)+x(jx,jz,jy,3)*x(jx,jz,jy,8)+eta(jx,jz,jy)*cur(jx,jz,jy,2)+xint(jx,jz,5)*xint(jx,jz,6)-xint(jx,jz,3)*xint(jx,jz,8)
!      ef(jx,jz,jy,3)=-x(jx,jz,jy,3)*x(jx,jz,jy,7)+x(jx,jz,jy,4)*x(jx,jz,jy,6)+eta(jx,jz,jy)*cur(jx,jz,jy,3)+xint(jx,jz,3)*xint(jx,jz,7)-xint(jx,jz,4)*xint(jx,jz,6)

      if(linear_mhd) then

          ef(jx,jz,jy,1)=-x1(jx,jz,jy,4)*xint(jx,jz,8)-xint(jx,jz,4)*x1(jx,jz,jy,8)+x1(jx,jz,jy,5)*xint(jx,jz,7) &
                +xint(jx,jz,5)*x1(jx,jz,jy,7)
          ef(jx,jz,jy,2)=-x1(jx,jz,jy,5)*xint(jx,jz,6)-xint(jx,jz,5)*x1(jx,jz,jy,6)+x1(jx,jz,jy,3)*xint(jx,jz,8) &
                +xint(jx,jz,3)*x1(jx,jz,jy,8)
          ef(jx,jz,jy,3)=-x1(jx,jz,jy,3)*xint(jx,jz,7)-xint(jx,jz,3)*x1(jx,jz,jy,7)+x1(jx,jz,jy,4)*xint(jx,jz,6) &
                +xint(jx,jz,4)*x1(jx,jz,jy,6)

	else

      ef(jx,jz,jy,1)=-x1(jx,jz,jy,4)*x(jx,jz,jy,8)-xint(jx,jz,4)*x1(jx,jz,jy,8)+x1(jx,jz,jy,5)*x(jx,jz,jy,7) &
                +xint(jx,jz,5)*x1(jx,jz,jy,7)
      ef(jx,jz,jy,2)=-x1(jx,jz,jy,5)*x(jx,jz,jy,6)-xint(jx,jz,5)*x1(jx,jz,jy,6)+x1(jx,jz,jy,3)*x(jx,jz,jy,8) &
                +xint(jx,jz,3)*x1(jx,jz,jy,8)
      ef(jx,jz,jy,3)=-x1(jx,jz,jy,3)*x(jx,jz,jy,7)-xint(jx,jz,3)*x1(jx,jz,jy,7)+x1(jx,jz,jy,4)*x(jx,jz,jy,6) &
                +xint(jx,jz,4)*x1(jx,jz,jy,6)

	endif
!      ef(jx,jz,jy,1)=ef(jx,jz,jy,1)+eta1(jx,jz,jy)*cint(jx,jz,1)
!      ef(jx,jz,jy,2)=ef(jx,jz,jy,2)+eta1(jx,jz,jy)*cint(jx,jz,2)
!      ef(jx,jz,jy,3)=ef(jx,jz,jy,3)+eta1(jx,jz,jy)*cint(jx,jz,3)

1     continue
      
       if(hall) then
      do 12 jy=iy_first+2,iy_last-2
      do 12 jz=iz_first+2,iz_last-2
      do 12 jx=ix_first+2,ix_last-2

      if(linear_mhd) then

      ef(jx,jz,jy,1)=ef(jx,jz,jy,1)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,2)*x1(jx,jz,jy,8)-cur(jx,jz,jy,3)*x1(jx,jz,jy,7) &
                +cint(jx,jz,2)*x1(jx,jz,jy,8)-cint(jx,jz,3)*x1(jx,jz,jy,7)-x1r(jx,jz,jy,2))
      ef(jx,jz,jy,2)=ef(jx,jz,jy,2)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,3)*x1(jx,jz,jy,6)-cur(jx,jz,jy,1)*x1(jx,jz,jy,8) &
                +cint(jx,jz,3)*x1(jx,jz,jy,6)-cint(jx,jz,1)*x1(jx,jz,jy,8)-xy(jx,jz,jy,2)/xx(jx))    
      ef(jx,jz,jy,3)=ef(jx,jz,jy,3)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,1)*x1(jx,jz,jy,7)-cur(jx,jz,jy,2)*x1(jx,jz,jy,6) &
                +cint(jx,jz,1)*x1(jx,jz,jy,7)-cint(jx,jz,2)*x1(jx,jz,jy,6)-x1z(jx,jz,jy,2))

	else

      ef(jx,jz,jy,1)=ef(jx,jz,jy,1)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,2)*x(jx,jz,jy,8)-cur(jx,jz,jy,3)*x(jx,jz,jy,7) &
                +cint(jx,jz,2)*x1(jx,jz,jy,8)-cint(jx,jz,3)*x1(jx,jz,jy,7)-x1r(jx,jz,jy,2))
      ef(jx,jz,jy,2)=ef(jx,jz,jy,2)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,3)*x(jx,jz,jy,6)-cur(jx,jz,jy,1)*x(jx,jz,jy,8) &
                +cint(jx,jz,3)*x1(jx,jz,jy,6)-cint(jx,jz,1)*x1(jx,jz,jy,8)-xy(jx,jz,jy,2)/xx(jx))    
      ef(jx,jz,jy,3)=ef(jx,jz,jy,3)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,1)*x(jx,jz,jy,7)-cur(jx,jz,jy,2)*x(jx,jz,jy,6) &
                +cint(jx,jz,1)*x1(jx,jz,jy,7)-cint(jx,jz,2)*x1(jx,jz,jy,6)-x1z(jx,jz,jy,2))

	endif
      
12     continue
       endif
       
!   !turn off the gradp term in hall terms to see the importance of it   
!       if(hall) then
!      do 12 jy=iy_first+2,iy_last-2
!      do 12 jz=iz_first+2,iz_last-2
!      do 12 jx=ix_first+2,ix_last-2
!
!      ef(jx,jz,jy,1)=ef(jx,jz,jy,1)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,2)*x(jx,jz,jy,8)-cur(jx,jz,jy,3)*x(jx,jz,jy,7)+cint(jx,jz,2)*x1(jx,jz,jy,8)-cint(jx,jz,3)*x1(jx,jz,jy,7))
!      ef(jx,jz,jy,2)=ef(jx,jz,jy,2)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,3)*x(jx,jz,jy,6)-cur(jx,jz,jy,1)*x(jx,jz,jy,8)+cint(jx,jz,3)*x1(jx,jz,jy,6)-cint(jx,jz,1)*x1(jx,jz,jy,8))    
!      ef(jx,jz,jy,3)=ef(jx,jz,jy,3)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,1)*x(jx,jz,jy,7)-cur(jx,jz,jy,2)*x(jx,jz,jy,6)+cint(jx,jz,1)*x1(jx,jz,jy,7)-cint(jx,jz,2)*x1(jx,jz,jy,6))
!      
!12     continue
!       endif
   !cd_oxpoint1:vs
      if(cd_oxpoint .and. lcdox==1 .and. mod(nstep,ncd).eq.0 .and. irk.eq.1) call distribution_cd_oxpoint(ef(:,:,3,2))

      if(resisitive .and. etaj_in_e) then
      do m=1,3
      ef(:,:,:,m)=ef(:,:,:,m)+eta(:,:,:)*cur(:,:,:,m)
      enddo
      endif
   !cd_oxpoint2:ey_nocd
      if(cd_oxpoint .and. lcdox==2 .and. mod(nstep,ncd).eq.0 .and. irk.eq.1) call distribution_cd_oxpoint(ef(:,:,3,2))

      if(bootstrap) then
      call current_boot(lbs)      
      do m=1,3
      ef(:,:,:,m)=ef(:,:,:,m)-eta(:,:,:)*cub(:,:,:,m)
      enddo
      endif
      
      if(lcd .eq. 1) then
      call current_driven
      do m=1,3
      ef(:,:,:,m)=ef(:,:,:,m)-eta(:,:,:)*cud(:,:,:,m)
      enddo
      endif

      if(nstep.lt.nper) then
      do 11 m=1,3
      do 11 jy=iy_first,iy_last
      do 11 jz=iz_first,iz_last
      do 11 jx=ix_first,ix_last
      ef(jx,jz,jy,m)=ef(jx,jz,jy,m)+eta1(jx,jz,jy)*(cint(jx,jz,m)+cur(jx,jz,jy,m))
  11  continue
      endif 




!      call ef_atlastgrid_r1p0_v1(3)

!      if(resisitive) then
!      if(etaj_in_e) then
!      do 2 jy=1,my
!      do 2 jz=iz_first,iz_last
!      do 2 jx=ix_first,ix_last
!      ef(jx,jz,jy,1)=ef(jx,jz,jy,1)+eta(jx,jz,jy)*(cur(jx,jz,jy,1))
!      ef(jx,jz,jy,2)=ef(jx,jz,jy,2)+eta(jx,jz,jy)*(cur(jx,jz,jy,2))
!      ef(jx,jz,jy,3)=ef(jx,jz,jy,3)+eta(jx,jz,jy)*(cur(jx,jz,jy,3))
!   2  continue
!      endif
!      endif

!      write(*,*) 'ef'
       
!      call convte
!      call recrd_ef
       !***************revised**************************************
       call mpi_transfersm(ef(:,:,:,:),3)
       !**********************************************************
      if(smoothef) call smthef_dis_v2(3)

      return
      end

!ws***************************************************************************
      subroutine convte
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

!       integer status(mpi_status_size)
!      init1=1
!      init2=1
!      jxs = ix_first
!      jxe = ix_last
!      jzs = iz_first
!      jze = iz_last

      if(spectral) then      

      else !!not spectral
      do 19 m=1,3
      do 15 jz=iz_first,iz_last
      do 15 jx=ix_first,ix_last
      do jy=iy_first+2,iy_last-2
      efy(jx,jz,jy,m)=d1f2(ef(jx,jz,jy-1,m),ef(jx,jz,jy,m),ef(jx,jz,jy+1,m),yy(jy-1),yy(jy),yy(jy+1))  
!      efy(jx,jz,jy,m)=d1fc(ef(jx,jz,jy-2,m),ef(jx,jz,jy-1,m),ef(jx,jz,jy,m) &
!          ,ef(jx,jz,jy+1,m),ef(jx,jz,jy+2,m),ay1(jy),by1(jy),cy1(jy),dy1(jy))              
      enddo
   15 continue
   19 continue 
      endif !!spectral

      do 30 m=1,3
      do 30 jy=iy_first,iy_last
      do 21 jz=iz_first+2,iz_last-2
      do 21 jx=ix_first+2,ix_last-2
!      if(rr(jx,jz).le.aa) then
      efx(jx,jz,jy,m)=d1f2(ef(jx-1,jz,jy,m),ef(jx,jz,jy,m) &
         ,ef(jx+1,jz,jy,m),xx(jx-1),xx(jx),xx(jx+1))
      efz(jx,jz,jy,m)=d1f2(ef(jx,jz-1,jy,m),ef(jx,jz,jy,m) &
         ,ef(jx,jz+1,jy,m),zz(jz-1),zz(jz),zz(jz+1))

!     efx(jx,jz,jy,m)=d1fc(ef(jx-2,jz,jy,m),ef(jx-1,jz,jy,m),ef(jx,jz,jy,m) &
!         ,ef(jx+1,jz,jy,m),ef(jx+2,jz,jy,m),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
!     efz(jx,jz,jy,m)=d1fc(ef(jx,jz-2,jy,m),ef(jx,jz-1,jy,m),ef(jx,jz,jy,m) &
!         ,ef(jx,jz+1,jy,m),ef(jx,jz+2,jy,m),az1(jz),bz1(jz),cz1(jz),dz1(jz))

   21 continue
   30 continue
!      call recrd_convte
      return
      end
	
!ws************************************************************
	subroutine convtdivb ! deleted now for nofftw in sunway(0810)
	end
	
	
!hwzhang*****************************************************************************
!hwzhang*****************************************************************************
!subroutines for cut cell method
!hwzhang*****************************************************************************
!hwzhang*****************************************************************************	
	
!hw******************************************************************
	subroutine right_cut_cell
      use declare
	integer itag
      real*8 drvx_dx,drey_dx,dez_dx,dex_dy,dez_dy,dex_dz,dey_dz
      real*8, dimension(mx,mz,my) :: rvx,rey
	real*8, dimension(nbndx,my) :: rvx_bndx, rey_bndx
	real*8, dimension(nbndz,my) :: rvx_bndz, rey_bndz
        real*8 d1fc, d1f2, d1fm, d1fbp ,d1fbm, d1xf2
      include 'mpif.h'
! R R: calculated as subroutine right
! R IR2: the IR2 part need to be calculated with the boundary point by d1fc, order 4th
! R IR1: the IR1 part need to be calculated with the boundary point by d1fbm, d1fbp, order 3th
! IR2 IR2: both two IR2 parts need to be calculated with the boundary point by d1fc, order 4th
! IR1 IR1: both tow IR1 parts need to be calculated with the boundary point by d1fbm, d1fbp, order 3th
! IR1 IR2: the IR1 part need to be calculated with the boundary point by d1fbm, d1fbp, order 3th
! 	     the IR2 parts need to be calculated with the boundary point by d1fc, order 4th
! R or IR1 or IR2 with D: interpolation in the dropped direction
! D D: do not calculate

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
        ((xm2-x0)/(xm1-x0)*(fm1-f0) &
         -(xm1-x0)/(xm2-x0)*(fm2-f0))/(xm2-xm1)
!  d1fbp= d f / dx  with  one-sided-bias  difference involving -1 0  1 and 2
!  points
      d1fbp(fp2,fp1,f0,fm1,a,b,c)= &
       a*(fm1-f0)+b*(fp1-f0)+c*(fp2-f0)
!  d1fbm= d f / dx  with  one-sided-bias  difference involving -2 -1  0 and 1
!  points
      d1fbm(fm2,fm1,f0,fp1,a,b,c)= &
       a*(f0-fp1)+b*(f0-fm1)+c*(f0-fm2)
!  d1xf2= d rf / dx  with second-order accuracy central difference
      d1xf2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(xp1*fp1-x0*f0) &
         -(xp1-x0)/(xm1-x0)*(xm1*fm1-x0*f0))/(xm1-xp1)

      call bndry8_cut_cell_v2_fixed
 	

      if(rho_from_p) then
      do jy=iy_first,iy_last
      do jz=iz_first,iz_last      
      do jx=ix_first,ix_last
!      if(psi(jx,jz).lt.psia1 .and. x(jx,jz,jy,2).gt.0) then
      if( gdtp_ep(jx,jz,1).ne.4 .and. x(jx,jz,jy,2).gt.0) then ! gdtp_ep=4 means the outside point
      x(jx,jz,jy,1)=xint(jx,jz,1)*(xint(jx,jz,2)/x(jx,jz,jy,2))**gamma
      else
      x(jx,jz,jy,1)=xint(jx,jz,1)
      endif

	enddo
	enddo
	enddo

	endif

      call convt_cut_cell
      call current_cut_cell
!      call convtdivb
	! haowei need to obtain the current for bnd grids, i.e., cur_3bndz,
	! cur_3bndx
      if(ef_mode) call efield_cut_cell_v2_fixed
	! haowei need to obatin the efield for bnd grids, i.e., ef_3bndz, ef_3bndx

      do jy=iy_first,iy_last
      do jz=iz_first,iz_last      
      do jx=ix_first,ix_last
      rvx(jx,jz,jy)=xx(jx)*x(jx,jz,jy,3)
      rey(jx,jz,jy)=xx(jx)*ef(jx,jz,jy,2)
      enddo
      enddo
      enddo

	do jy=iy_first,iy_last
	!calculate the rvx_bndz and rey_bndz
	do itag=1,nbndz
	rvx_bndz(itag,jy)=bndz_grd(itag,1)*x_8bndz(itag,jy,3)
	rey_bndz(itag,jy)=bndz_grd(itag,1)*ef_3bndz(itag,jy,2)
	enddo
	!calculate the rvx_bndx and rey_bndx
	do itag=1,nbndx
	rvx_bndx(itag,jy)=bndx_grd(itag,1)*x_8bndx(itag,jy,3)
	rey_bndx(itag,jy)=bndx_grd(itag,1)*ef_3bndx(itag,jy,2)
	enddo
	enddo


      do 1 jy=iy_first+2,iy_last-2
      do 1 jz=iz_first+2,iz_last-2      
      do 1 jx=ix_first+2,ix_last-2
!      if((psi(jx,jz).lt.psia1).and.(max(gdtp_ep(jx,jz,1),gdtp_ep(jx,jz,4)).ne.5)) then
      if((gdtp_ep(jx,jz,1).ne.4).and.(max(gdtp_ep(jx,jz,1),gdtp_ep(jx,jz,4)).ne.5)) then
        drvx_dx=d1fc(rvx(jx-2,jz,jy),rvx(jx-1,jz,jy),rvx(jx,jz,jy),rvx(jx+1,jz,jy),rvx(jx+2,jz,jy),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
	! hw for IR in x direction, drvx_dx need to be recalculated
	if(gdtp_ep(jx,jz,4).eq.2) then 
		  itag=gdtp_ep(jx,jz,6)
		  if(gdtp_ep(jx,jz,5).eq.-2) then ! x1(jx-2,jz,jy,m)=bndz
			    drvx_dx=d1fc(rvx_bndz(itag,jy),rvx(jx-1,jz,jy),rvx(jx,jz,jy)&
					,rvx(jx+1,jz,jy),rvx(jx+2,jz,jy),ax1_irz(jx,jz),bx1_irz(jx,jz)&
					,cx1_irz(jx,jz),dx1_irz(jx,jz))
		  endif
		  if(gdtp_ep(jx,jz,5).eq.2) then ! x1(jx+2,jz,jy,m)=bndz
			    drvx_dx=d1fc(rvx(jx-2,jz,jy),rvx(jx-1,jz,jy),rvx(jx,jz,jy)&
					,rvx(jx+1,jz,jy),rvx_bndz(itag,jy),ax1_irz(jx,jz),bx1_irz(jx,jz)&
					,cx1_irz(jx,jz),dx1_irz(jx,jz))
		  endif
		  if(gdtp_ep(jx,jz,5).eq.-1) then !x1(jx-1,jz,jy,m)=bndz d1fbp d2fbp
			    drvx_dx=d1fbp(rvx(jx+2,jz,jy),rvx(jx+1,jz,jy),rvx(jx,jz,jy)&
					,rvx_bndz(itag,jy),axbp_irz(jx,jz),bxbp_irz(jx,jz),cxbp_irz(jx,jz))
		  endif
		  if(gdtp_ep(jx,jz,5).eq.1) then !x1(jx+1,jz,jy,m)=bndz d1fbm d2fbm
			    drvx_dx=d1fbm(rvx(jx-2,jz,jy),rvx(jx-1,jz,jy),rvx(jx,jz,jy)&
					,rvx_bndz(itag,jy),axbm_irz(jx,jz),bxbm_irz(jx,jz),cxbm_irz(jx,jz))
		  endif
	endif


!clt option for linear mhd simulation!!!!
      if(linear_mhd) then      

        xdif(jx,jz,jy,1)=-x1(jx,jz,jy,3)*xint_dx(jx,jz,1) &
            -xint(jx,jz,1)*drvx_dx/xx(jx) &
            -(xy(jx,jz,jy,1)*xint(jx,jz,4)+xint(jx,jz,1)*xy(jx,jz,jy,4))/xx(jx) &
            -x1(jx,jz,jy,5)*xint_dz(jx,jz,1)-xint(jx,jz,1)*x1z(jx,jz,jy,5) &
    !ws:for poloidal flow
            -xint(jx,jz,3)*x1r(jx,jz,jy,1)-x1(jx,jz,jy,1)*(xint(jx,jz,3)/xx(jx)+xint_dx(jx,jz,3)) &
            -xint(jx,jz,5)*x1z(jx,jz,jy,1)-x1(jx,jz,jy,1)*xint_dz(jx,jz,5)

            xdif(jx,jz,jy,2)=-x1(jx,jz,jy,3)*xint_dx(jx,jz,2) &
            -gamma*xint(jx,jz,2)*(drvx_dx/xx(jx)) &
            -(xy(jx,jz,jy,2)*xint(jx,jz,4)+gamma*xint(jx,jz,2)*xy(jx,jz,jy,4))/xx(jx) &
            -x1(jx,jz,jy,5)*xint_dz(jx,jz,2)-gamma*xint(jx,jz,2)*x1z(jx,jz,jy,5) &    
    !ws:for poloidal flow
            -xint(jx,jz,3)*x1r(jx,jz,jy,2)-gamma*x1(jx,jz,jy,2)*(xint(jx,jz,3)/xx(jx)+xint_dx(jx,jz,3)) &
            -xint(jx,jz,5)*x1z(jx,jz,jy,2)-gamma*x1(jx,jz,jy,2)*xint_dz(jx,jz,5)       
     
             
            xdif(jx,jz,jy,3)=-xint(jx,jz,3)*x1r(jx,jz,jy,3)-xint(jx,jz,5)*x1z(jx,jz,jy,3) &
            -xint(jx,jz,4)*xy(jx,jz,jy,3)/xx(jx)+xint(jx,jz,4)*x1(jx,jz,jy,4)/xx(jx) &
            +((cur(jx,jz,jy,2))*xint(jx,jz,8)+(cint(jx,jz,2))*x1(jx,jz,jy,8) &
            - (cur(jx,jz,jy,3))*xint(jx,jz,7)-(cint(jx,jz,3))*x1(jx,jz,jy,7) &
            -x1r(jx,jz,jy,2))/xint(jx,jz,1) &
            -x1(jx,jz,jy,3)*xint_dx(jx,jz,3)-x1(jx,jz,jy,5)*xint_dz(jx,jz,3) &
            +x1(jx,jz,jy,4)*xint(jx,jz,4)/xx(jx) &
            +(-xint(jx,jz,3)*xint_dx(jx,jz,3)-xint(jx,jz,5)*xint_dz(jx,jz,3) &
            +xint(jx,jz,4)*xint(jx,jz,4)/xx(jx))*x1(jx,jz,jy,1)/xint(jx,jz,1)


            xdif(jx,jz,jy,4)=-xint(jx,jz,3)*x1r(jx,jz,jy,4)-xint(jx,jz,5)*x1z(jx,jz,jy,4) &
            -xint(jx,jz,4)*xy(jx,jz,jy,4)/xx(jx)-xint(jx,jz,4)*x1(jx,jz,jy,3)/xx(jx) &
            +((cur(jx,jz,jy,3))*xint(jx,jz,6)+(cint(jx,jz,3))*x1(jx,jz,jy,6) &
            - (cur(jx,jz,jy,1))*xint(jx,jz,8)-(cint(jx,jz,1))*x1(jx,jz,jy,8) &   
            -xy(jx,jz,jy,2)/xx(jx))/xint(jx,jz,1) &
            -x1(jx,jz,jy,3)*xint_dx(jx,jz,4)-x1(jx,jz,jy,5)*xint_dz(jx,jz,4) &
            -x1(jx,jz,jy,4)*xint(jx,jz,3)/xx(jx) &
            +(-xint(jx,jz,3)*xint_dx(jx,jz,4)-xint(jx,jz,5)*xint_dz(jx,jz,4) &
            -xint(jx,jz,4)*xint(jx,jz,3)/xx(jx))*x1(jx,jz,jy,1)/xint(jx,jz,1)
     
       
            xdif(jx,jz,jy,5)=-xint(jx,jz,3)*x1r(jx,jz,jy,5)-xint(jx,jz,5)*x1z(jx,jz,jy,5) &
            -xint(jx,jz,4)*xy(jx,jz,jy,5)/xx(jx) &     
            +((cur(jx,jz,jy,1))*xint(jx,jz,7)+(cint(jx,jz,1))*x1(jx,jz,jy,7) &
            - (cur(jx,jz,jy,2))*xint(jx,jz,6)-(cint(jx,jz,2))*x1(jx,jz,jy,6) & 
            -x1z(jx,jz,jy,2))/xint(jx,jz,1) &
            -x1(jx,jz,jy,3)*xint_dx(jx,jz,5)-x1(jx,jz,jy,5)*xint_dz(jx,jz,5) &
            +(-xint(jx,jz,3)*xint_dx(jx,jz,5)-xint(jx,jz,5)*xint_dz(jx,jz,5)) &
            *x1(jx,jz,jy,1)/xint(jx,jz,1)


	else


      xdif(jx,jz,jy,1)=-x1(jx,jz,jy,3)*xr(jx,jz,jy,1) &
       -x(jx,jz,jy,1)*drvx_dx/xx(jx) &
!rplc       -x(jx,jz,jy,3)*x(jx,jz,jy,1)/xx(jx)-x(jx,jz,jy,1)*xr(jx,jz,jy,3) &
       -(xy(jx,jz,jy,1)*x(jx,jz,jy,4)+x(jx,jz,jy,1)*xy(jx,jz,jy,4))/xx(jx) &
       -x1(jx,jz,jy,5)*xz(jx,jz,jy,1)-x(jx,jz,jy,1)*x1z(jx,jz,jy,5) &
!ws:for poloidal flow
       -xint(jx,jz,3)*x1r(jx,jz,jy,1)-x1(jx,jz,jy,1)*(xint(jx,jz,3)/xx(jx)+xint_dx(jx,jz,3)) &
       -xint(jx,jz,5)*x1z(jx,jz,jy,1)-x1(jx,jz,jy,1)*xint_dz(jx,jz,5)

      xdif(jx,jz,jy,2)=-x1(jx,jz,jy,3)*xr(jx,jz,jy,2) &
       -gamma*x(jx,jz,jy,2)*(drvx_dx/xx(jx)) &
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
          
	endif
     
!      xdif(jx,jz,jy,6)=-x(jx,jz,jy,3)*xr(jx,jz,jy,6)-x(jx,jz,jy,6)*x(jx,jz,jy,3)/xx(jx) &
!       -(x(jx,jz,jy,4)*xy(jx,jz,jy,6)-x(jx,jz,jy,7)*xy(jx,jz,jy,3)+x(jx,jz,jy,6)*xy(jx,jz,jy,4))/xx(jx) &
!       -x(jx,jz,jy,5)*xz(jx,jz,jy,6)+x(jx,jz,jy,8)*xz(jx,jz,jy,3)-x(jx,jz,jy,6)*xz(jx,jz,jy,5) !&
!!j x grd_eta
!!       -(cur(jx,jz,jy,2)-cj(jx,jz))*etaz(jx,jz,jy)
!     
!      xdif(jx,jz,jy,7)=-x(jx,jz,jy,3)*xr(jx,jz,jy,7)+x(jx,jz,jy,6)*xr(jx,jz,jy,4) &
!       -x(jx,jz,jy,6)*x(jx,jz,jy,4)/xx(jx)-x(jx,jz,jy,7)*xr(jx,jz,jy,3) &
!       -xy(jx,jz,jy,7)*x(jx,jz,jy,4)/xx(jx) &
!       -xz(jx,jz,jy,7)*x(jx,jz,jy,5)+x(jx,jz,jy,8)*xz(jx,jz,jy,4)-x(jx,jz,jy,7)*xz(jx,jz,jy,5) !&
!!j x grd_eta
!!       -cur(jx,jz,jy,3)*etax(jx,jz,jy)+cur(jx,jz,jy,1)*etaz(jx,jz,jy)
!      
!      xdif(jx,jz,jy,8)=-x(jx,jz,jy,3)*xr(jx,jz,jy,8)-x(jx,jz,jy,3)*x(jx,jz,jy,8)/xx(jx) &
!       -x(jx,jz,jy,8)*xr(jx,jz,jy,3)+x(jx,jz,jy,6)*xr(jx,jz,jy,5) &
!       -(xy(jx,jz,jy,8)*x(jx,jz,jy,4)-xy(jx,jz,jy,5)*x(jx,jz,jy,7)+xy(jx,jz,jy,4)*x(jx,jz,jy,8))/xx(jx) &
!       -xz(jx,jz,jy,8)*x(jx,jz,jy,5) !&
!!j x grd_eta
!!       +(cur(jx,jz,jy,2)-cj(jx,jz))*etax(jx,jz,jy)

!      xdif(jx,jz,jy,6)=-efy(jx,jz,jy,3)/xx(jx)+efz(jx,jz,jy,2) !&
!!                      +cdb0*divb_x(jx,jz,jy)
!    
!      xdif(jx,jz,jy,7)=-efz(jx,jz,jy,1)+efx(jx,jz,jy,3) !&
!!                      +cdb0*divb_y(jx,jz,jy)/xx(jx)
!
!      xdif(jx,jz,jy,8)=efy(jx,jz,jy,1)/xx(jx)-efx(jx,jz,jy,2)-ef(jx,jz,jy,2)/xx(jx) !&
!!                      +cdb0*divb_z(jx,jz,jy)

      if(ef_mode) then
      !drey_dx=d1xf2(ef(jx-1,jz,jy,2),ef(jx,jz,jy,2),ef(jx+1,jz,jy,2),xx(jx-1),xx(jx),xx(jx+1))
      !dez_dx =d1f2(ef(jx-1,jz,jy,3),ef(jx,jz,jy,3),ef(jx+1,jz,jy,3),xx(jx-1),xx(jx),xx(jx+1))
      !dex_dz =d1f2(ef(jx,jz-1,jy,1),ef(jx,jz,jy,1),ef(jx,jz+1,jy,1),zz(jz-1),zz(jz),zz(jz+1))
      !dey_dz =d1f2(ef(jx,jz-1,jy,2),ef(jx,jz,jy,2),ef(jx,jz+1,jy,2),zz(jz-1),zz(jz),zz(jz+1))
      !dex_dy =d1f2(ef(jx,jz,jy-1,1),ef(jx,jz,jy,1),ef(jx,jz,jy+1,1),yy(jy-1),yy(jy),yy(jy+1))
      !dez_dy =d1f2(ef(jx,jz,jy-1,3),ef(jx,jz,jy,3),ef(jx,jz,jy+1,3),yy(jy-1),yy(jy),yy(jy+1))

     drey_dx=d1fc(rey(jx-2,jz,jy),rey(jx-1,jz,jy),rey(jx,jz,jy) &
         ,rey(jx+1,jz,jy),rey(jx+2,jz,jy),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
     dez_dx =d1fc(ef(jx-2,jz,jy,3),ef(jx-1,jz,jy,3),ef(jx,jz,jy,3) &
         ,ef(jx+1,jz,jy,3),ef(jx+2,jz,jy,3),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
     dex_dz =d1fc(ef(jx,jz-2,jy,1),ef(jx,jz-1,jy,1),ef(jx,jz,jy,1) &
         ,ef(jx,jz+1,jy,1),ef(jx,jz+2,jy,1),az1(jz),bz1(jz),cz1(jz),dz1(jz))
     dey_dz =d1fc(ef(jx,jz-2,jy,2),ef(jx,jz-1,jy,2),ef(jx,jz,jy,2) &
         ,ef(jx,jz+1,jy,2),ef(jx,jz+2,jy,2),az1(jz),bz1(jz),cz1(jz),dz1(jz))
     dex_dy =d1fc(ef(jx,jz,jy-2,1),ef(jx,jz,jy-1,1),ef(jx,jz,jy,1) &
         ,ef(jx,jz,jy+1,1),ef(jx,jz,jy+2,1),ay1(jy),by1(jy),cy1(jy),dy1(jy)) 
     dez_dy =d1fc(ef(jx,jz,jy-2,3),ef(jx,jz,jy-1,3),ef(jx,jz,jy,3) &
         ,ef(jx,jz,jy+1,3),ef(jx,jz,jy+2,3),ay1(jy),by1(jy),cy1(jy),dy1(jy)) 

	! hw for IR in x direction, drey_dx dez_dx need to be recalculated
	if(gdtp_ep(jx,jz,4).eq.2) then 
		  itag=gdtp_ep(jx,jz,6)
		  if(gdtp_ep(jx,jz,5).eq.-2) then ! x1(jx-2,jz,jy,m)=bndz
			    drey_dx=d1fc(rey_bndz(itag,jy),rey(jx-1,jz,jy),rey(jx,jz,jy)&
					,rey(jx+1,jz,jy),rey(jx+2,jz,jy),ax1_irz(jx,jz),bx1_irz(jx,jz)&
					,cx1_irz(jx,jz),dx1_irz(jx,jz))
			    dez_dx=d1fc(ef_3bndz(itag,jy,3),ef(jx-1,jz,jy,3),ef(jx,jz,jy,3)&
					,ef(jx+1,jz,jy,3),ef(jx+2,jz,jy,3),ax1_irz(jx,jz),bx1_irz(jx,jz)&
					,cx1_irz(jx,jz),dx1_irz(jx,jz))
		  endif
		  if(gdtp_ep(jx,jz,5).eq.2) then ! x1(jx+2,jz,jy,m)=bndz
			    drey_dx=d1fc(rey(jx-2,jz,jy),rey(jx-1,jz,jy),rey(jx,jz,jy)&
					,rey(jx+1,jz,jy),rey_bndz(itag,jy),ax1_irz(jx,jz),bx1_irz(jx,jz)&
					,cx1_irz(jx,jz),dx1_irz(jx,jz))
			    dez_dx=d1fc(ef(jx-2,jz,jy,3),ef(jx-1,jz,jy,3),ef(jx,jz,jy,3)&
					,ef(jx+1,jz,jy,3),ef_3bndz(itag,jy,3),ax1_irz(jx,jz),bx1_irz(jx,jz)&
					,cx1_irz(jx,jz),dx1_irz(jx,jz))
		  endif
		  if(gdtp_ep(jx,jz,5).eq.-1) then !x1(jx-1,jz,jy,m)=bndz d1fbp d2fbp
			    drey_dx=d1fbp(rey(jx+2,jz,jy),rey(jx+1,jz,jy),rey(jx,jz,jy)&
					,rey_bndz(itag,jy),axbp_irz(jx,jz),bxbp_irz(jx,jz),cxbp_irz(jx,jz))
			    dez_dx=d1fbp(ef(jx+2,jz,jy,3),ef(jx+1,jz,jy,3),ef(jx,jz,jy,3)&
					,ef_3bndz(itag,jy,3),axbp_irz(jx,jz),bxbp_irz(jx,jz),cxbp_irz(jx,jz))
		  endif
		  if(gdtp_ep(jx,jz,5).eq.1) then !x1(jx+1,jz,jy,m)=bndz d1fbm d2fbm
			    drey_dx=d1fbm(rey(jx-2,jz,jy),rey(jx-1,jz,jy),rey(jx,jz,jy)&
					,rey_bndz(itag,jy),axbm_irz(jx,jz),bxbm_irz(jx,jz),cxbm_irz(jx,jz))
			    dez_dx=d1fbm(ef(jx-2,jz,jy,3),ef(jx-1,jz,jy,3),ef(jx,jz,jy,3)&
					,ef_3bndz(itag,jy,3),axbm_irz(jx,jz),bxbm_irz(jx,jz),cxbm_irz(jx,jz))
		  endif
	endif

	! for IR in z direction, dex_dz dey_dz need to be recalculated
	if(gdtp_ep(jx,jz,1).eq.2) then 
		  itag=gdtp_ep(jx,jz,3)
		  if(gdtp_ep(jx,jz,2).eq.-2) then ! x1(jx,jz-2,jy,m)=bndx
			    dex_dz=d1fc(ef_3bndx(itag,jy,1),ef(jx,jz-1,jy,1),ef(jx,jz,jy,1)&
					,ef(jx,jz+1,jy,1),ef(jx,jz+2,jy,1),az1_irx(jx,jz),bz1_irx(jx,jz)&
					,cz1_irx(jx,jz),dz1_irx(jx,jz))
			    dey_dz=d1fc(ef_3bndx(itag,jy,2),ef(jx,jz-1,jy,2),ef(jx,jz,jy,2)&
					,ef(jx,jz+1,jy,2),ef(jx,jz+2,jy,2),az1_irx(jx,jz),bz1_irx(jx,jz)&
					,cz1_irx(jx,jz),dz1_irx(jx,jz))
		  endif
		  if(gdtp_ep(jx,jz,2).eq.2) then ! x1(jx,jz+2,jy,m)=bndx
			    dex_dz=d1fc(ef(jx,jz-2,jy,1),ef(jx,jz-1,jy,1),ef(jx,jz,jy,1)&
					,ef(jx,jz+1,jy,1),ef_3bndx(itag,jy,1),az1_irx(jx,jz),bz1_irx(jx,jz)&
					,cz1_irx(jx,jz),dz1_irx(jx,jz))
			    dey_dz=d1fc(ef(jx,jz-2,jy,2),ef(jx,jz-1,jy,2),ef(jx,jz,jy,2)&
					,ef(jx,jz+1,jy,2),ef_3bndx(itag,jy,2),az1_irx(jx,jz),bz1_irx(jx,jz)&
					,cz1_irx(jx,jz),dz1_irx(jx,jz))
		  endif
		  if(gdtp_ep(jx,jz,2).eq.-1) then ! x1(jx,jz-1,jy,m)=bndx d1fbp d2fbp
			    dex_dz=d1fbp(ef(jx,jz+2,jy,1),ef(jx,jz+1,jy,1),ef(jx,jz,jy,1)&
					,ef_3bndx(itag,jy,1),azbp_irx(jx,jz),bzbp_irx(jx,jz),czbp_irx(jx,jz))
			    dey_dz=d1fbp(ef(jx,jz+2,jy,2),ef(jx,jz+1,jy,2),ef(jx,jz,jy,2)&
					,ef_3bndx(itag,jy,2),azbp_irx(jx,jz),bzbp_irx(jx,jz),czbp_irx(jx,jz))
		  endif
		  if(gdtp_ep(jx,jz,2).eq.1) then ! x1(jx,jz+1,jy,m)=bndx d1fbm d2fbm
			    dex_dz=d1fbm(ef(jx,jz-2,jy,1),ef(jx,jz-1,jy,1),ef(jx,jz,jy,1)&
					,ef_3bndx(itag,jy,1),azbm_irx(jx,jz),bzbm_irx(jx,jz),czbm_irx(jx,jz))
			    dey_dz=d1fbm(ef(jx,jz-2,jy,2),ef(jx,jz-1,jy,2),ef(jx,jz,jy,2)&
					,ef_3bndx(itag,jy,2),azbm_irx(jx,jz),bzbm_irx(jx,jz),czbm_irx(jx,jz))
		  endif
	endif
          
      xdif(jx,jz,jy,6)=-dez_dy/xx(jx)+dey_dz !&
!                       +eta(jx,jz,jy)*cint_dz(jx,jz,2)  
      xdif(jx,jz,jy,7)=-dex_dz+dez_dx !&
!                       +eta(jx,jz,jy)*(-cint_dz(jx,jz,1)+cint_dx(jx,jz,3))
      xdif(jx,jz,jy,8)=(dex_dy-drey_dx)/xx(jx) !&
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

      xdif(jx,jz,jy,8)=(x(jx,jz,jy,5)*drvx_dx/xx(jx)+x(jx,jz,jy,6)*xr(jx,jz,jy,5) &
                       -x(jx,jz,jy,3)*xr(jx,jz,jy,8)-x(jx,jz,jy,8)*drvx_dx/xx(jx)) &    
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
      if(ef_mode .and. smoothbout) then
!      if(.not.implicitb) then
      do 41 jy=iy_first+2,iy_last-2
      do 41 jz=iz_first+2,iz_last-2
      do 41 jx=ix_first+2,ix_last-2
!      if(psi(jx,jz).lt.psia1) then ! keep the etbx and etbz only in closed flux
      if(gdtp_ep(jx,jz,1).ne.4) then
!!ws:db/dt=...+eta*grd2 b
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
      if(.not.etaj_in_e) then
!      if(.not.implicitb) then
      do 4 jy=iy_first+2,iy_last-2
      do 4 jz=iz_first+2,iz_last-2
      do 4 jx=ix_first+2,ix_last-2
!      if(psi(jx,jz).lt.psia1) then  ! keep the etaj only in closed flux
      if(gdtp_ep(jx,jz,1).ne.4) then
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

!!ws:db/dt=...-eta*grd x j
!      xdif(jx,jz,jy,6)=xdif(jx,jz,jy,6)+cur(jx,jz,jy,2)*etaz(jx,jz,jy) &
!       -eta(jx,jz,jy)*(cuy(jx,jz,jy,3)/xx(jx)-cuz(jx,jz,jy,2))
!     
!      xdif(jx,jz,jy,7)=xdif(jx,jz,jy,7)+cur(jx,jz,jy,3)*etax(jx,jz,jy)-cur(jx,jz,jy,1)*etaz(jx,jz,jy) &
!       -eta(jx,jz,jy)*(cuz(jx,jz,jy,1)-cux(jx,jz,jy,3))
!            
!      xdif(jx,jz,jy,8)=xdif(jx,jz,jy,8)-cur(jx,jz,jy,2)*etax(jx,jz,jy) &
!       -eta(jx,jz,jy)*(cux(jx,jz,jy,2)+cur(jx,jz,jy,2)/xx(jx)-cuy(jx,jz,jy,1)/xx(jx))
!!ws:db/dt=...+eta*grd2 b
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
!      if(psi(jx,jz).lt.psia1) then
      if(gdtp_ep(jx,jz,1).ne.4) then
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
!      if(psi(jx,jz).lt.psia1) then ! keep the viscous term only in closed flux
      if(gdtp_ep(jx,jz,1).ne.4) then
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

!hw******************************************************************************
      subroutine stepon_cut_cell
!
!     this routine time-advances x's bz fourth order in time and second
!     order in space runge-kotta differential scheme.
!     note: x is alwazs the up-to-date value while xm being the
!           intermediate value, and xdif is increment
!
!
      use declare
	implicit none
      include 'mpif.h'

!      dts=dt/ncycl_atfs

      tt=time
      tt1=time+dt/6.
      tt2=time
      irk=1
      call right_cut_cell
       xfold(:,:,:,:)=x(:,:,:,:)
       xm(:,:,:,:)=xfold(:,:,:,:)+xdif(:,:,:,:)*dt/6.
       x(:,:,:,:)=xfold(:,:,:,:)+xdif(:,:,:,:)*dt/2.
!
      tt=time+dt/2.
      tt1=time+dt/2.
      tt2=time+dt/6.
      irk=2
        call right_cut_cell
        xm(:,:,:,:)=xm(:,:,:,:)+xdif(:,:,:,:)*dt/3.
        x(:,:,:,:)=xfold(:,:,:,:)+xdif(:,:,:,:)*dt/2.
!
      tt1=time+5.*dt/6.
      tt2=time+dt/2.
      irk=3
        call right_cut_cell
        xm(:,:,:,:)=xm(:,:,:,:)+xdif(:,:,:,:)*dt/3.
        x(:,:,:,:)=xfold(:,:,:,:)+xdif(:,:,:,:)*dt
!
      time=time+dt
      tt1=time+dt
      tt2=time+5.*dt/6.
      irk=4
        call right_cut_cell
        x(:,:,:,:)=xm(:,:,:,:)+xdif(:,:,:,:)*dt/6.

      caf=0.75d0*(0.5+0.5*dtanh((time-40)/5.))
!      call bndry_x_ex(lbnd)
!      call bndry8_x_ex(lbnd)
	call bndry8_cut_cell_v2_fixed
  !    if(conductpll) call pllconduct(lpll)
!      if(smoothpll) call smthp_traceline_5p(1)
!      if(eta_from_t) call calculate_eta
      
      return
      end
	
!hw***************************************************************************************
	subroutine convt_cut_cell
      use declare
	integer itag, i
      real*8, dimension(my) :: wwy 
      real*8, dimension(mx,mz,my) :: x1r_tmp,xr2_tmp,x1z_tmp,xz2_tmp,x1_tmp
      real*8 d1f2, d1f2m, d1f2p, d1fc, d2fc, d1fp, d1fm, d1fbp, d1fbm
      real*8 d2f2, d2fbp, d2fbm, timestart1, timeend1
      include 'mpif.h'
! R R: calculated as subroutine right
! R IR2: the IR2 part need to be calculated with the boundary point by d1fc, order 4th
! R IR1: the IR1 part need to be calculated with the boundary point by d1fbm, d1fbp, order 3th
! IR2 IR2: both two IR2 parts need to be calculated with the boundary point by d1fc, order 4th
! IR1 IR1: both tow IR1 parts need to be calculated with the boundary point by d1fbm, d1fbp, order 3th
! IR1 IR2: the IR1 part need to be calculated with the boundary point by d1fbm, d1fbp, order 3th
! 	     the IR2 parts need to be calculated with the boundary point by d1fc, order 4th
! R or IR1 or IR2 with D: interpolation in the dropped direction
! D D: do not calculate

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
	 !!not spectral, need to add spectral if necessary, haowei
      do 19 m=1,8
      do 15 jz=iz_first,iz_last
      do 15 jx=ix_first,ix_last
      do jy=iy_first+2,iy_last-2
      xy(jx,jz,jy,m) =d1fc(x1(jx,jz,jy-2,m),x1(jx,jz,jy-1,m),x1(jx,jz,jy,m),x1(jx,jz,jy+1,m),x1(jx,jz,jy+2,m),ay1(jy),by1(jy),cy1(jy),dy1(jy))      
      xy2(jx,jz,jy,m)=d2fc(x1(jx,jz,jy-2,m),x1(jx,jz,jy-1,m),x1(jx,jz,jy,m),x1(jx,jz,jy+1,m),x1(jx,jz,jy+2,m),ay2(jy),by2(jy),cy2(jy),dy2(jy))
      enddo
   15 continue
! for the xy_8bndx and xy2_8bndx and xy_8bndz and xy2_8bndz
	do 16 jx=1,nbndx
	do jy=iy_first+2,iy_last-2
	xy_8bndx(jx,jy,m)=d1fc(x1_8bndx(jx,jy-2,m),x1_8bndx(jx,jy-1,m),x1_8bndx(jx,jy,m),x1_8bndx(jx,jy+1,m),x1_8bndx(jx,jy+2,m),&
		  ay1(jy),by1(jy),cy1(jy),dy1(jy))
	xy2_8bndx(jx,jy,m)=d2fc(x1_8bndx(jx,jy-2,m),x1_8bndx(jx,jy-1,m),x1_8bndx(jx,jy,m),x1_8bndx(jx,jy+1,m),x1_8bndx(jx,jy+2,m),&
		  ay2(jy),by2(jy),cy2(jy),dy2(jy))
	xy_8bndz(jx,jy,m)=d1fc(x1_8bndz(jx,jy-2,m),x1_8bndz(jx,jy-1,m),x1_8bndz(jx,jy,m),x1_8bndz(jx,jy+1,m),x1_8bndz(jx,jy+2,m),&
		  ay1(jy),by1(jy),cy1(jy),dy1(jy))
	xy2_8bndz(jx,jy,m)=d2fc(x1_8bndz(jx,jy-2,m),x1_8bndz(jx,jy-1,m),x1_8bndz(jx,jy,m),x1_8bndz(jx,jy+1,m),x1_8bndz(jx,jy+2,m),&
		  ay2(jy),by2(jy),cy2(jy),dy2(jy))
	enddo
   16 continue
   19 continue 

      !not spectral

      do m=1,8
      do jy=iy_first,iy_last
      do jz=iz_first+2,iz_last-2
      do jx=ix_first+2,ix_last-2

     x1r(jx,jz,jy,m)=d1fc(x1(jx-2,jz,jy,m),x1(jx-1,jz,jy,m),x1(jx,jz,jy,m) &
         ,x1(jx+1,jz,jy,m),x1(jx+2,jz,jy,m),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
     xr2(jx,jz,jy,m)=d2fc(x1(jx-2,jz,jy,m),x1(jx-1,jz,jy,m),x1(jx,jz,jy,m) &
         ,x1(jx+1,jz,jy,m),x1(jx+2,jz,jy,m),ax2(jx),bx2(jx),cx2(jx),dx2(jx))
     x1z(jx,jz,jy,m)=d1fc(x1(jx,jz-2,jy,m),x1(jx,jz-1,jy,m),x1(jx,jz,jy,m) &
         ,x1(jx,jz+1,jy,m),x1(jx,jz+2,jy,m),az1(jz),bz1(jz),cz1(jz),dz1(jz))
     xz2(jx,jz,jy,m)=d2fc(x1(jx,jz-2,jy,m),x1(jx,jz-1,jy,m),x1(jx,jz,jy,m) &
         ,x1(jx,jz+1,jy,m),x1(jx,jz+2,jy,m),az2(jz),bz2(jz),cz2(jz),dz2(jz))


      enddo
      enddo
      enddo
      enddo


      do 30 m=1,8
      do 30 jy=iy_first,iy_last
      do 21 jz=iz_first+2,iz_last-2
      do 21 jx=ix_first+2,ix_last-2

	! for IR in x direction, x1r xr2 need to be recalculated
	if(gdtp_ep(jx,jz,4).eq.2) then 
		  itag=gdtp_ep(jx,jz,6)
		  if(gdtp_ep(jx,jz,5).eq.-2) then ! x1(jx-2,jz,jy,m)=bndz
			    x1r(jx,jz,jy,m)=d1fc(x1_8bndz(itag,jy,m),x1(jx-1,jz,jy,m),x1(jx,jz,jy,m)&
					,x1(jx+1,jz,jy,m),x1(jx+2,jz,jy,m),ax1_irz(jx,jz),bx1_irz(jx,jz)&
					,cx1_irz(jx,jz),dx1_irz(jx,jz))
			    xr2(jx,jz,jy,m)=d2fc(x1_8bndz(itag,jy,m),x1(jx-1,jz,jy,m),x1(jx,jz,jy,m)&
					,x1(jx+1,jz,jy,m),x1(jx+2,jz,jy,m),ax2_irz(jx,jz),bx2_irz(jx,jz)&
					,cx2_irz(jx,jz),dx2_irz(jx,jz))
		  endif
		  if(gdtp_ep(jx,jz,5).eq.2) then ! x1(jx+2,jz,jy,m)=bndz
			    x1r(jx,jz,jy,m)=d1fc(x1(jx-2,jz,jy,m),x1(jx-1,jz,jy,m),x1(jx,jz,jy,m)&
					,x1(jx+1,jz,jy,m),x1_8bndz(itag,jy,m),ax1_irz(jx,jz),bx1_irz(jx,jz)&
					,cx1_irz(jx,jz),dx1_irz(jx,jz))
			    xr2(jx,jz,jy,m)=d2fc(x1(jx-2,jz,jy,m),x1(jx-1,jz,jy,m),x1(jx,jz,jy,m)&
					,x1(jx+1,jz,jy,m),x1_8bndz(itag,jy,m),ax2_irz(jx,jz),bx2_irz(jx,jz)&
					,cx2_irz(jx,jz),dx2_irz(jx,jz))
		  endif
		  if(gdtp_ep(jx,jz,5).eq.-1) then !x1(jx-1,jz,jy,m)=bndz d1fbp d2fbp
			    x1r(jx,jz,jy,m)=d1fbp(x1(jx+2,jz,jy,m),x1(jx+1,jz,jy,m),x1(jx,jz,jy,m)&
					,x1_8bndz(itag,jy,m),axbp_irz(jx,jz),bxbp_irz(jx,jz),cxbp_irz(jx,jz))
			    xr2(jx,jz,jy,m)=d2fbp(x1_8bndz(itag,jy,m),x1(jx,jz,jy,m),x1(jx+1,jz,jy,m)&
					,x1(jx+2,jz,jy,m),a2xbp_irz(jx,jz),b2xbp_irz(jx,jz),c2xbp_irz(jx,jz))
		  endif
		  if(gdtp_ep(jx,jz,5).eq.1) then !x1(jx+1,jz,jy,m)=bndz d1fbm d2fbm
			    x1r(jx,jz,jy,m)=d1fbm(x1(jx-2,jz,jy,m),x1(jx-1,jz,jy,m),x1(jx,jz,jy,m)&
					,x1_8bndz(itag,jy,m),axbm_irz(jx,jz),bxbm_irz(jx,jz),cxbm_irz(jx,jz))
			    xr2(jx,jz,jy,m)=d2fbm(x1(jx-2,jz,jy,m),x1(jx-1,jz,jy,m),x1(jx,jz,jy,m)&
					,x1_8bndz(itag,jy,m),a2xbm_irz(jx,jz),b2xbm_irz(jx,jz),c2xbm_irz(jx,jz))
		  endif
	endif

	! for IR in z direction, x1z xz2 need to be recalculated
	if(gdtp_ep(jx,jz,1).eq.2) then 
		  itag=gdtp_ep(jx,jz,3)
		  if(gdtp_ep(jx,jz,2).eq.-2) then ! x1(jx,jz-2,jy,m)=bndx
			    x1z(jx,jz,jy,m)=d1fc(x1_8bndx(itag,jy,m),x1(jx,jz-1,jy,m),x1(jx,jz,jy,m)&
					,x1(jx,jz+1,jy,m),x1(jx,jz+2,jy,m),az1_irx(jx,jz),bz1_irx(jx,jz)&
					,cz1_irx(jx,jz),dz1_irx(jx,jz))
			    xz2(jx,jz,jy,m)=d2fc(x1_8bndx(itag,jy,m),x1(jx,jz-1,jy,m),x1(jx,jz,jy,m)&
					,x1(jx,jz+1,jy,m),x1(jx,jz+2,jy,m),az2_irx(jx,jz),bz2_irx(jx,jz)&
					,cz2_irx(jx,jz),dz2_irx(jx,jz))
		  endif
		  if(gdtp_ep(jx,jz,2).eq.2) then ! x1(jx,jz+2,jy,m)=bndx
			    x1z(jx,jz,jy,m)=d1fc(x1(jx,jz-2,jy,m),x1(jx,jz-1,jy,m),x1(jx,jz,jy,m)&
					,x1(jx,jz+1,jy,m),x1_8bndx(itag,jy,m),az1_irx(jx,jz),bz1_irx(jx,jz)&
					,cz1_irx(jx,jz),dz1_irx(jx,jz))
			    xz2(jx,jz,jy,m)=d2fc(x1(jx,jz-2,jy,m),x1(jx,jz-1,jy,m),x1(jx,jz,jy,m)&
					,x1(jx,jz+1,jy,m),x1_8bndx(itag,jy,m),az2_irx(jx,jz),bz2_irx(jx,jz)&
					,cz2_irx(jx,jz),dz2_irx(jx,jz))
		  endif
		  if(gdtp_ep(jx,jz,2).eq.-1) then ! x1(jx,jz-1,jy,m)=bndx d1fbp d2fbp
			    x1z(jx,jz,jy,m)=d1fbp(x1(jx,jz+2,jy,m),x1(jx,jz+1,jy,m),x1(jx,jz,jy,m)&
					,x1_8bndx(itag,jy,m),azbp_irx(jx,jz),bzbp_irx(jx,jz),czbp_irx(jx,jz))
			    xz2(jx,jz,jy,m)=d2fbp(x1_8bndx(itag,jy,m),x1(jx,jz,jy,m),x1(jx,jz+1,jy,m)&
					,x1(jx,jz+2,jy,m),a2zbp_irx(jx,jz),b2zbp_irx(jx,jz),c2zbp_irx(jx,jz))
		  endif
		  if(gdtp_ep(jx,jz,2).eq.1) then ! x1(jx,jz+1,jy,m)=bndx d1fbm d2fbm
			    x1z(jx,jz,jy,m)=d1fbm(x1(jx,jz-2,jy,m),x1(jx,jz-1,jy,m),x1(jx,jz,jy,m)&
					,x1_8bndx(itag,jy,m),azbm_irx(jx,jz),bzbm_irx(jx,jz),czbm_irx(jx,jz))
			    xz2(jx,jz,jy,m)=d2fbm(x1(jx,jz-2,jy,m),x1(jx,jz-1,jy,m),x1(jx,jz,jy,m)&
					,x1_8bndx(itag,jy,m),a2zbm_irx(jx,jz),b2zbm_irx(jx,jz),c2zbm_irx(jx,jz))
		  endif
	endif
	xr(jx,jz,jy,m)=xint_dx(jx,jz,m)+x1r(jx,jz,jy,m)  
      xz(jx,jz,jy,m)=xint_dz(jx,jz,m)+x1z(jx,jz,jy,m)
	
!	if((gdtp_ep(jx,jz,1).eq.5).or.(gdtp_ep(jx,jz,4).eq.5)) then
!		  ! need to do for the dropped but inside point by interpolation,
!		  ! but no derivation for x1r xr2 x1z xz2, haha, need not to
!		  ! calculate the interpolation
!		  ! haowei 12/18
!		  if(gdtp_ep(jx-1,jz,1).eq.2) then ! in x direction, it is R R IR IR DI B DO, the DI need to be interpolated by jx-2, jx-1, jx, B
!		  endif
!
!		  if(gdtp_ep(jx+1,jz,1).eq.2) then ! in x direction, it is DO B DI IR IR R R, the DI need to be interpolated by B, jx, jx+1, jx+2
!		  endif
!
!		  if(gdtp_ep(jx,jz-1,4).eq.2) then ! in z direction, it is R R IR IR DI B DO, the DI need to be interpolated by jz-2, jz-1, jz, B
!		  endif
!
!		  if(gdtp_ep(jx,jz+1,4).eq.2) then ! in z direction, it is DO B DI IR IR R R, the DI need to be interpolated by B, jz, jz+1, jz+2
!		  endif
!	endif

   21 continue
   30 continue

    return
    end

!hw*************************************************************************
      subroutine current_cut_cell
      use declare
	integer itag
      real*8 drby_dx, rby_tmp
      real*8, dimension(mx,mz,my) :: rby
      include 'mpif.h'
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
!  d1xf2= d rf / dx  with second-order accuracy central difference
      d1xf2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(xp1*fp1-x0*f0) &
         -(xp1-x0)/(xm1-x0)*(xm1*fm1-x0*f0))/(xm1-xp1)

!       integer status(mpi_status_size)

!      do 10 m=1,3
      do 10 jy=iy_first,iy_last
      do 10 jz=iz_first,iz_last
      do 10 jx=ix_first,ix_last
      rby(jx,jz,jy)=xx(jx)*x1(jx,jz,jy,7)
   10 continue

      do 1 jy=iy_first+2,iy_last-2
      do 1 jz=iz_first+2,iz_last-2
      do 1 jx=ix_first+2,ix_last-2
!      if(psi(jx,jz).lt.psia1) then
      if(gdtp_ep(jx,jz,1).ne.4) then
      cur(jx,jz,jy,1)=xy(jx,jz,jy,8)/xx(jx)-x1z(jx,jz,jy,7) 
      cur(jx,jz,jy,2)=x1z(jx,jz,jy,6)-x1r(jx,jz,jy,8)
!      cur(jx,jz,jy,3)=x1r(jx,jz,jy,7)+(x1(jx,jz,jy,7)-xy(jx,jz,jy,6))/xx(jx)
      drby_dx=d1fc(rby(jx-2,jz,jy),rby(jx-1,jz,jy),rby(jx,jz,jy) &
         ,rby(jx+1,jz,jy),rby(jx+2,jz,jy),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
  !    drby_dx=d1xf2(x1(jx-1,jz,jy,7),x1(jx,jz,jy,7),x1(jx+1,jz,jy,7),xx(jx-1),xx(jx),xx(jx+1))
  	
	! for IR in x direction, drby_dx need to be recalculated
  	if(gdtp_ep(jx,jz,4).eq.2) then
		  itag=gdtp_ep(jx,jz,6)
		  rby_tmp=bndz_grd(itag,1)*x1_8bndz(itag,jy,7)
		  if(gdtp_ep(jx,jz,5).eq.-2) then !jx-2 is the bndz
			    drby_dx=d1fc(rby_tmp,rby(jx-1,jz,jy),rby(jx,jz,jy) &
					,rby(jx+1,jz,jy),rby(jx+2,jz,jy) &
					,ax1_irz(jx,jz),bx1_irz(jx,jz),cx1_irz(jx,jz),dx1_irz(jx,jz))
		  endif
		  if(gdtp_ep(jx,jz,5).eq.2) then !jx+2 is the bndz
			    drby_dx=d1fc(rby(jx-2,jz,jy),rby(jx-1,jz,jy),rby(jx,jz,jy) &
					,rby(jx+1,jz,jy),rby_tmp &
					,ax1_irz(jx,jz),bx1_irz(jx,jz),cx1_irz(jx,jz),dx1_irz(jx,jz))
		  endif
		  if(gdtp_ep(jx,jz,5).eq.-1) then !jx-1 is the bndz
			    drby_dx=d1fbp(rby(jx+2,jz,jy),rby(jx+1,jz,jy),rby(jx,jz,jy),rby_tmp &
					,axbp_irz(jx,jz),bxbp_irz(jx,jz),cxbp_irz(jx,jz))
		  endif
		  if(gdtp_ep(jx,jz,5).eq.1) then !jx+1 is the bndz
			    drby_dx=d1fbm(rby(jx-2,jz,jy),rby(jx-1,jz,jy),rby(jx,jz,jy),rby_tmp &
					,axbm_irz(jx,jz),bxbm_irz(jx,jz),cxbm_irz(jx,jz))
		  endif
	endif

      cur(jx,jz,jy,3)=(drby_dx-xy(jx,jz,jy,6))/xx(jx)

!      cur(jx,jz,jy,1)=d1f2(xh(jx,jz-1,jy,2),xh(jx,jz,jy,2),xh(jx,jz+1,jy,2),zz(jz-1),zz(jz),zz(jz+1)) &
!       -xy(jx,jz,jy,8)/xx(jx) 
!      cur(jx,jz,jy,2)=d1f2(xh(jx-1,jz,jy,3),xh(jx,jz,jy,3),xh(jx+1,jz,jy,3),xx(jx-1),xx(jx),xx(jx+1)) &
!       -d1f2(xh(jx,jz-1,jy,1),xh(jx,jz,jy,1),xh(jx,jz+1,jy,1),zz(jz-1),zz(jz),zz(jz+1))
!      cur(jx,jz,jy,3)=-d1f2(xh(jx-1,jz,jy,2),xh(jx,jz,jy,2),xh(jx+1,jz,jy,2),xx(jx-1),xx(jx),xx(jx+1)) &
!       -(xh(jx,jz,jy,2)-xy(jx,jz,jy,6))/xx(jx)
      endif
    1 continue
    	call mpi_transfersm(cur(:,:,:,:),3)
	if(smoothc) then
		  do m=1,3
		  call smthxzy(cur(:,:,:,m),1)
		  enddo
	endif
!      call bndry_cur_ex(lbnd)
!      call bndry3_ex(cur,lbnd)
	! need to be replaced by bndry3_ex_cut_cell   ! haowei 1219
!      call mpi_transfer3(cur) 
!      call convtc
      if(lcd .eq. 2) then 
      call current_driven
      cur(:,:,:,:)=cur(:,:,:,:)+cud(:,:,:,:)
      endif
!      write(*,*)'cur'
      return
      end    
	
!hw************************************************************************
!ws************************************************************************
!wzhang************************************************************
      subroutine efield_cut_cell
      use declare
	implicit none
	integer itag
      include 'mpif.h'


      do 1 jy=iy_first,iy_last
      do 1 jz=iz_first,iz_last
      do 1 jx=ix_first,ix_last
!      ef(jx,jz,jy,1)=-x(jx,jz,jy,4)*x(jx,jz,jy,8)+x(jx,jz,jy,5)*x(jx,jz,jy,7)+eta(jx,jz,jy)*cur(jx,jz,jy,1)+xint(jx,jz,4)*xint(jx,jz,8)-xint(jx,jz,5)*xint(jx,jz,7)
!      ef(jx,jz,jy,2)=-x(jx,jz,jy,5)*x(jx,jz,jy,6)+x(jx,jz,jy,3)*x(jx,jz,jy,8)+eta(jx,jz,jy)*cur(jx,jz,jy,2)+xint(jx,jz,5)*xint(jx,jz,6)-xint(jx,jz,3)*xint(jx,jz,8)
!      ef(jx,jz,jy,3)=-x(jx,jz,jy,3)*x(jx,jz,jy,7)+x(jx,jz,jy,4)*x(jx,jz,jy,6)+eta(jx,jz,jy)*cur(jx,jz,jy,3)+xint(jx,jz,3)*xint(jx,jz,7)-xint(jx,jz,4)*xint(jx,jz,6)


      if(linear_mhd) then

          ef(jx,jz,jy,1)=-x1(jx,jz,jy,4)*xint(jx,jz,8)-xint(jx,jz,4)*x1(jx,jz,jy,8)+x1(jx,jz,jy,5)*xint(jx,jz,7)+xint(jx,jz,5)*x1(jx,jz,jy,7)
          ef(jx,jz,jy,2)=-x1(jx,jz,jy,5)*xint(jx,jz,6)-xint(jx,jz,5)*x1(jx,jz,jy,6)+x1(jx,jz,jy,3)*xint(jx,jz,8)+xint(jx,jz,3)*x1(jx,jz,jy,8)
          ef(jx,jz,jy,3)=-x1(jx,jz,jy,3)*xint(jx,jz,7)-xint(jx,jz,3)*x1(jx,jz,jy,7)+x1(jx,jz,jy,4)*xint(jx,jz,6)+xint(jx,jz,4)*x1(jx,jz,jy,6)

	else

      ef(jx,jz,jy,1)=-x1(jx,jz,jy,4)*x(jx,jz,jy,8)-xint(jx,jz,4)*x1(jx,jz,jy,8)+x1(jx,jz,jy,5)*x(jx,jz,jy,7)+xint(jx,jz,5)*x1(jx,jz,jy,7)
      ef(jx,jz,jy,2)=-x1(jx,jz,jy,5)*x(jx,jz,jy,6)-xint(jx,jz,5)*x1(jx,jz,jy,6)+x1(jx,jz,jy,3)*x(jx,jz,jy,8)+xint(jx,jz,3)*x1(jx,jz,jy,8)
      ef(jx,jz,jy,3)=-x1(jx,jz,jy,3)*x(jx,jz,jy,7)-xint(jx,jz,3)*x1(jx,jz,jy,7)+x1(jx,jz,jy,4)*x(jx,jz,jy,6)+xint(jx,jz,4)*x1(jx,jz,jy,6)
	endif
!      ef(jx,jz,jy,1)=ef(jx,jz,jy,1)+eta1(jx,jz,jy)*cint(jx,jz,1)
!      ef(jx,jz,jy,2)=ef(jx,jz,jy,2)+eta1(jx,jz,jy)*cint(jx,jz,2)
!      ef(jx,jz,jy,3)=ef(jx,jz,jy,3)+eta1(jx,jz,jy)*cint(jx,jz,3)

1     continue
      
       if(hall) then
      do 12 jy=iy_first+2,iy_last-2
      do 12 jz=iz_first+2,iz_last-2
      do 12 jx=ix_first+2,ix_last-2

      if(linear_mhd) then

      ef(jx,jz,jy,1)=ef(jx,jz,jy,1)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,2)*x1(jx,jz,jy,8)-cur(jx,jz,jy,3)*x1(jx,jz,jy,7)+cint(jx,jz,2)*x1(jx,jz,jy,8)-cint(jx,jz,3)*x1(jx,jz,jy,7)-x1r(jx,jz,jy,2))
      ef(jx,jz,jy,2)=ef(jx,jz,jy,2)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,3)*x1(jx,jz,jy,6)-cur(jx,jz,jy,1)*x1(jx,jz,jy,8)+cint(jx,jz,3)*x1(jx,jz,jy,6)-cint(jx,jz,1)*x1(jx,jz,jy,8)-xy(jx,jz,jy,2)/xx(jx))    
      ef(jx,jz,jy,3)=ef(jx,jz,jy,3)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,1)*x1(jx,jz,jy,7)-cur(jx,jz,jy,2)*x1(jx,jz,jy,6)+cint(jx,jz,1)*x1(jx,jz,jy,7)-cint(jx,jz,2)*x1(jx,jz,jy,6)-x1z(jx,jz,jy,2))

	else

      ef(jx,jz,jy,1)=ef(jx,jz,jy,1)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,2)*x(jx,jz,jy,8)-cur(jx,jz,jy,3)*x(jx,jz,jy,7)+cint(jx,jz,2)*x1(jx,jz,jy,8)-cint(jx,jz,3)*x1(jx,jz,jy,7)-x1r(jx,jz,jy,2))
      ef(jx,jz,jy,2)=ef(jx,jz,jy,2)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,3)*x(jx,jz,jy,6)-cur(jx,jz,jy,1)*x(jx,jz,jy,8)+cint(jx,jz,3)*x1(jx,jz,jy,6)-cint(jx,jz,1)*x1(jx,jz,jy,8)-xy(jx,jz,jy,2)/xx(jx))    
      ef(jx,jz,jy,3)=ef(jx,jz,jy,3)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,1)*x(jx,jz,jy,7)-cur(jx,jz,jy,2)*x(jx,jz,jy,6)+cint(jx,jz,1)*x1(jx,jz,jy,7)-cint(jx,jz,2)*x1(jx,jz,jy,6)-x1z(jx,jz,jy,2))
      
	endif
12     continue
       endif
       
!   !turn off the gradp term in hall terms to see the importance of it   
!       if(hall) then
!      do 12 jy=iy_first+2,iy_last-2
!      do 12 jz=iz_first+2,iz_last-2
!      do 12 jx=ix_first+2,ix_last-2
!
!      ef(jx,jz,jy,1)=ef(jx,jz,jy,1)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,2)*x(jx,jz,jy,8)-cur(jx,jz,jy,3)*x(jx,jz,jy,7)+cint(jx,jz,2)*x1(jx,jz,jy,8)-cint(jx,jz,3)*x1(jx,jz,jy,7))
!      ef(jx,jz,jy,2)=ef(jx,jz,jy,2)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,3)*x(jx,jz,jy,6)-cur(jx,jz,jy,1)*x(jx,jz,jy,8)+cint(jx,jz,3)*x1(jx,jz,jy,6)-cint(jx,jz,1)*x1(jx,jz,jy,8))    
!      ef(jx,jz,jy,3)=ef(jx,jz,jy,3)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,1)*x(jx,jz,jy,7)-cur(jx,jz,jy,2)*x(jx,jz,jy,6)+cint(jx,jz,1)*x1(jx,jz,jy,7)-cint(jx,jz,2)*x1(jx,jz,jy,6))
!      
!12     continue
!       endif
   !cd_oxpoint1:vs
      if(cd_oxpoint .and. lcdox==1 .and. mod(nstep,ncd).eq.0 .and. irk.eq.1) call distribution_cd_oxpoint(ef(:,:,3,2))

      if(resisitive .and. etaj_in_e) then
      do m=1,3
      ef(:,:,:,m)=ef(:,:,:,m)+eta(:,:,:)*cur(:,:,:,m)
      enddo
      endif
   !cd_oxpoint2:ey_nocd
      if(cd_oxpoint .and. lcdox==2 .and. mod(nstep,ncd).eq.0 .and. irk.eq.1) call distribution_cd_oxpoint(ef(:,:,3,2))

      if(bootstrap) then
      call current_boot(lbs)      
      do m=1,3
      ef(:,:,:,m)=ef(:,:,:,m)-eta(:,:,:)*cub(:,:,:,m)
      enddo
      endif
      
      if(lcd .eq. 1) then
      call current_driven
      do m=1,3
      ef(:,:,:,m)=ef(:,:,:,m)-eta(:,:,:)*cud(:,:,:,m)
      enddo
      endif

      if(nstep.lt.nper) then
      do 11 m=1,3
      do 11 jy=iy_first,iy_last
      do 11 jz=iz_first,iz_last
      do 11 jx=ix_first,ix_last
      ef(jx,jz,jy,m)=ef(jx,jz,jy,m)+eta1(jx,jz,jy)*(cint(jx,jz,m)+cur(jx,jz,jy,m))
  11  continue
      endif 




!      call ef_atlastgrid_r1p0_v1(3)

!      if(resisitive) then
!      if(etaj_in_e) then
!      do 2 jy=1,my
!      do 2 jz=iz_first,iz_last
!      do 2 jx=ix_first,ix_last
!      ef(jx,jz,jy,1)=ef(jx,jz,jy,1)+eta(jx,jz,jy)*(cur(jx,jz,jy,1))
!      ef(jx,jz,jy,2)=ef(jx,jz,jy,2)+eta(jx,jz,jy)*(cur(jx,jz,jy,2))
!      ef(jx,jz,jy,3)=ef(jx,jz,jy,3)+eta(jx,jz,jy)*(cur(jx,jz,jy,3))
!   2  continue
!      endif
!      endif

!      write(*,*) 'ef'
       
!      call convte
!      call recrd_ef
       !***************revised**************************************
       call mpi_transfersm(ef(:,:,:,:),3)
    ! haowei need to obatin the efield for bnd grids, i.e., ef_3bndz, ef_3bndx

    do 3 jz=iz_first,iz_last
!    
    do jx=ix_first+2,ix_last

    itag=gdtp_ep(jx,jz,6)
    if((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz,4).eq.5)) ef(jx,jz,:,:)=0.d0

    if(gdtp_ep(jx,jz,5).eq.1) then ! jx+1 is the boundary point, use the jx, jx-1, jx-2 to calculate the bndz point

    do 21 m=1,3
    do 21 jy=1,my
!    ef_3bndz(itag,jy,m)=(axm_bndz(itag)*ef(jx,jz,jy,m)+bxm_bndz(itag)*ef(jx-1,jz,jy,m)+ &
!            cxm_bndz(itag)*ef(jx-2,jz,jy,m))/(axm_bndz(itag)+bxm_bndz(itag)+cxm_bndz(itag))
    ! fixed
    ef_3bndz(itag,jy,m)=0.d0
  21 continue

    else if((gdtp_ep(jx,jz,4).eq.5).and.(gdtp_ep(jx,jz,1).ne.5).and.(gdtp_ep(jx-1,jz,5).eq.1)) then ! jx is the inside dropped point
    
    do 22 m=1,3
    do 22 jy=1,my
    call interp1d2l(ef(jx-2,jz,jy,m),ef(jx-1,jz,jy,m),ef_3bndz(itag,jy,m), &
            xx(jx-2),xx(jx-1),bndz_grd(itag,1),xx(jx),ef(jx,jz,jy,m))
  22 continue
    endif
    enddo
!
    do jx=ix_first,ix_last-2

    itag=gdtp_ep(jx,jz,6)
    if((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz,4).eq.5)) ef(jx,jz,:,:)=0.d0

    if(gdtp_ep(jx,jz,5).eq.-1) then ! jx-1 is the boundary point, use the jx, jx+1, jx+2 to calculate the bndz point

    do 23 m=1,3
    do 23 jy=1,my
!    ef_3bndz(itag,jy,m)=(axp_bndz(itag)*ef(jx,jz,jy,m)+bxp_bndz(itag)*ef(jx+1,jz,jy,m)+ &
!            cxp_bndz(itag)*ef(jx+2,jz,jy,m))/(axp_bndz(itag)+bxp_bndz(itag)+cxp_bndz(itag))
    ! fixed
    ef_3bndz(itag,jy,m)=0.d0
  23 continue

    else if((gdtp_ep(jx,jz,4).eq.5).and.(gdtp_ep(jx,jz,1).ne.5).and.(gdtp_ep(jx+1,jz,5).eq.-1)) then ! jx is the inside dropped point

    do 24 m=1,3
    do 24 jy=1,my
    call interp1d2l(ef(jx+2,jz,jy,m),ef(jx+1,jz,jy,m),ef_3bndz(itag,jy,m), &
            xx(jx+2),xx(jx+1),bndz_grd(itag,1),xx(jx),ef(jx,jz,jy,m))
  24 continue
    endif
    enddo

    3 continue
! 

    do 4 jx=ix_first,ix_last

    do jz=iz_first+2,iz_last

    itag=gdtp_ep(jx,jz,3)
    if((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz,4).eq.5)) ef(jx,jz,:,:)=0.d0

    if(gdtp_ep(jx,jz,2).eq.1) then ! jz+1 is the boundary point, use the jz, jz-1, jz-2 to calculate the bndx point

    do 25 m=1,3
    do 25 jy=1,my
!    ef_3bndx(itag,jy,m)=(azm_bndx(itag)*ef(jx,jz,jy,m)+bzm_bndx(itag)*ef(jx,jz-1,jy,m)+ &
!            czm_bndx(itag)*ef(jx,jz-2,jy,m))/(azm_bndx(itag)+bzm_bndx(itag)+czm_bndx(itag))
    ! fixed
    ef_3bndx(itag,jy,m)=0.d0
  25 continue
    else if((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz,4).ne.5).and.(gdtp_ep(jx,jz-1,2).eq.1)) then ! jz is the inside dropped point
    do 26 m=1,3
    do 26 jy=1,my
    call interp1d2l(ef(jx,jz-2,jy,m),ef(jx,jz-1,jy,m),ef_3bndx(itag,jy,m), &
            zz(jz-2),zz(jz-1),bndx_grd(itag,2),zz(jz),ef(jx,jz,jy,m))
  26 continue
    endif
    enddo

!
    do jz=iz_first,iz_last-2

    itag=gdtp_ep(jx,jz,3)
    if((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz,4).eq.5)) ef(jx,jz,:,:)=0.d0

    if(gdtp_ep(jx,jz,2).eq.-1) then ! jz-1 is the boundary point, use the jz, jz+1, jz+2 to calculate the bndx point

    do 27 m=1,3
    do 27 jy=1,my
!    ef_3bndx(itag,jy,m)=(azp_bndx(itag)*ef(jx,jz,jy,m)+bzp_bndx(itag)*ef(jx,jz+1,jy,m)+ &
!            czp_bndx(itag)*ef(jx,jz+2,jy,m))/(azp_bndx(itag)+bzp_bndx(itag)+czp_bndx(itag))
    ! fixed
    ef_3bndx(itag,jy,m)=0.d0
  27 continue

    else if((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz,4).ne.5).and.(gdtp_ep(jx,jz+1,2).eq.-1)) then ! jz is the inside dropped point
    do 28 m=1,3
    do 28 jy=1,my
    call interp1d2l(ef(jx,jz+2,jy,m),ef(jx,jz+1,jy,m),ef_3bndx(itag,jy,m), &
            zz(jz+2),zz(jz+1),bndx_grd(itag,2),zz(jz),ef(jx,jz,jy,m))
  28 continue
    endif
    enddo

    4 continue

!	do 10 m=1,3
!      do 10 jy=1,my
!      do 10 jz=iz_first,iz_last
!      do 10 jx=ix_first,ix_last
!      if(psi(jx,jz).lt.psia1) then
!	if((hypb_ratio(jx,jz).ge.0.d0).and.(hypb_ratio(jx,jz).le.1.d0)) then
!	ef(jx,jz,jy,m)=ef(jx,jz,jy,m)*hypb_ratio(jx,jz)
!	endif
!	endif
!   10 continue
!       call mpi_transfersm(ef(:,:,:,:),3)
       !**********************************************************
      if(smoothef) call smthef_dis_v2(3)

!	call smth_irpt_with_difc(4,2,1,3,0.9d0)
!	call smth_irpt_with_difc_v2(4,2,1,3,0.9d0)

      return
      end


!hw************************************************************************
!ws************************************************************************
!wzhang************************************************************
      subroutine efield_cut_cell_v2_fixed
      use declare
	implicit none
	integer itag
      include 'mpif.h'


      do 1 jy=iy_first,iy_last
      do 1 jz=iz_first,iz_last
      do 1 jx=ix_first,ix_last
!      ef(jx,jz,jy,1)=-x(jx,jz,jy,4)*x(jx,jz,jy,8)+x(jx,jz,jy,5)*x(jx,jz,jy,7)+eta(jx,jz,jy)*cur(jx,jz,jy,1)+xint(jx,jz,4)*xint(jx,jz,8)-xint(jx,jz,5)*xint(jx,jz,7)
!      ef(jx,jz,jy,2)=-x(jx,jz,jy,5)*x(jx,jz,jy,6)+x(jx,jz,jy,3)*x(jx,jz,jy,8)+eta(jx,jz,jy)*cur(jx,jz,jy,2)+xint(jx,jz,5)*xint(jx,jz,6)-xint(jx,jz,3)*xint(jx,jz,8)
!      ef(jx,jz,jy,3)=-x(jx,jz,jy,3)*x(jx,jz,jy,7)+x(jx,jz,jy,4)*x(jx,jz,jy,6)+eta(jx,jz,jy)*cur(jx,jz,jy,3)+xint(jx,jz,3)*xint(jx,jz,7)-xint(jx,jz,4)*xint(jx,jz,6)

      if(linear_mhd) then

          ef(jx,jz,jy,1)=-x1(jx,jz,jy,4)*xint(jx,jz,8)-xint(jx,jz,4)*x1(jx,jz,jy,8)+x1(jx,jz,jy,5)*xint(jx,jz,7)+xint(jx,jz,5)*x1(jx,jz,jy,7)
          ef(jx,jz,jy,2)=-x1(jx,jz,jy,5)*xint(jx,jz,6)-xint(jx,jz,5)*x1(jx,jz,jy,6)+x1(jx,jz,jy,3)*xint(jx,jz,8)+xint(jx,jz,3)*x1(jx,jz,jy,8)
          ef(jx,jz,jy,3)=-x1(jx,jz,jy,3)*xint(jx,jz,7)-xint(jx,jz,3)*x1(jx,jz,jy,7)+x1(jx,jz,jy,4)*xint(jx,jz,6)+xint(jx,jz,4)*x1(jx,jz,jy,6)

	else

      ef(jx,jz,jy,1)=-x1(jx,jz,jy,4)*x(jx,jz,jy,8)-xint(jx,jz,4)*x1(jx,jz,jy,8)+x1(jx,jz,jy,5)*x(jx,jz,jy,7)+xint(jx,jz,5)*x1(jx,jz,jy,7)
      ef(jx,jz,jy,2)=-x1(jx,jz,jy,5)*x(jx,jz,jy,6)-xint(jx,jz,5)*x1(jx,jz,jy,6)+x1(jx,jz,jy,3)*x(jx,jz,jy,8)+xint(jx,jz,3)*x1(jx,jz,jy,8)
      ef(jx,jz,jy,3)=-x1(jx,jz,jy,3)*x(jx,jz,jy,7)-xint(jx,jz,3)*x1(jx,jz,jy,7)+x1(jx,jz,jy,4)*x(jx,jz,jy,6)+xint(jx,jz,4)*x1(jx,jz,jy,6)

	endif
!      ef(jx,jz,jy,1)=ef(jx,jz,jy,1)+eta1(jx,jz,jy)*cint(jx,jz,1)
!      ef(jx,jz,jy,2)=ef(jx,jz,jy,2)+eta1(jx,jz,jy)*cint(jx,jz,2)
!      ef(jx,jz,jy,3)=ef(jx,jz,jy,3)+eta1(jx,jz,jy)*cint(jx,jz,3)

1     continue
      
       if(hall) then
      do 12 jy=iy_first+2,iy_last-2
      do 12 jz=iz_first+2,iz_last-2
      do 12 jx=ix_first+2,ix_last-2

      if(linear_mhd) then

      ef(jx,jz,jy,1)=ef(jx,jz,jy,1)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,2)*x1(jx,jz,jy,8)-cur(jx,jz,jy,3)*x1(jx,jz,jy,7)+cint(jx,jz,2)*x1(jx,jz,jy,8)-cint(jx,jz,3)*x1(jx,jz,jy,7)-x1r(jx,jz,jy,2))
      ef(jx,jz,jy,2)=ef(jx,jz,jy,2)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,3)*x1(jx,jz,jy,6)-cur(jx,jz,jy,1)*x1(jx,jz,jy,8)+cint(jx,jz,3)*x1(jx,jz,jy,6)-cint(jx,jz,1)*x1(jx,jz,jy,8)-xy(jx,jz,jy,2)/xx(jx))    
      ef(jx,jz,jy,3)=ef(jx,jz,jy,3)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,1)*x1(jx,jz,jy,7)-cur(jx,jz,jy,2)*x1(jx,jz,jy,6)+cint(jx,jz,1)*x1(jx,jz,jy,7)-cint(jx,jz,2)*x1(jx,jz,jy,6)-x1z(jx,jz,jy,2))

	else

      ef(jx,jz,jy,1)=ef(jx,jz,jy,1)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,2)*x(jx,jz,jy,8)-cur(jx,jz,jy,3)*x(jx,jz,jy,7)+cint(jx,jz,2)*x1(jx,jz,jy,8)-cint(jx,jz,3)*x1(jx,jz,jy,7)-x1r(jx,jz,jy,2))
      ef(jx,jz,jy,2)=ef(jx,jz,jy,2)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,3)*x(jx,jz,jy,6)-cur(jx,jz,jy,1)*x(jx,jz,jy,8)+cint(jx,jz,3)*x1(jx,jz,jy,6)-cint(jx,jz,1)*x1(jx,jz,jy,8)-xy(jx,jz,jy,2)/xx(jx))    
      ef(jx,jz,jy,3)=ef(jx,jz,jy,3)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,1)*x(jx,jz,jy,7)-cur(jx,jz,jy,2)*x(jx,jz,jy,6)+cint(jx,jz,1)*x1(jx,jz,jy,7)-cint(jx,jz,2)*x1(jx,jz,jy,6)-x1z(jx,jz,jy,2))
      
	endif
12     continue
       endif
       
!   !turn off the gradp term in hall terms to see the importance of it   
!       if(hall) then
!      do 12 jy=iy_first+2,iy_last-2
!      do 12 jz=iz_first+2,iz_last-2
!      do 12 jx=ix_first+2,ix_last-2
!
!      ef(jx,jz,jy,1)=ef(jx,jz,jy,1)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,2)*x(jx,jz,jy,8)-cur(jx,jz,jy,3)*x(jx,jz,jy,7)+cint(jx,jz,2)*x1(jx,jz,jy,8)-cint(jx,jz,3)*x1(jx,jz,jy,7))
!      ef(jx,jz,jy,2)=ef(jx,jz,jy,2)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,3)*x(jx,jz,jy,6)-cur(jx,jz,jy,1)*x(jx,jz,jy,8)+cint(jx,jz,3)*x1(jx,jz,jy,6)-cint(jx,jz,1)*x1(jx,jz,jy,8))    
!      ef(jx,jz,jy,3)=ef(jx,jz,jy,3)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,1)*x(jx,jz,jy,7)-cur(jx,jz,jy,2)*x(jx,jz,jy,6)+cint(jx,jz,1)*x1(jx,jz,jy,7)-cint(jx,jz,2)*x1(jx,jz,jy,6))
!      
!12     continue
!       endif
   !cd_oxpoint1:vs
      if(cd_oxpoint .and. lcdox==1 .and. mod(nstep,ncd).eq.0 .and. irk.eq.1) call distribution_cd_oxpoint(ef(:,:,3,2))

      if(resisitive .and. etaj_in_e) then
      do m=1,3
      ef(:,:,:,m)=ef(:,:,:,m)+eta(:,:,:)*cur(:,:,:,m)
      enddo
      endif
   !cd_oxpoint2:ey_nocd
      if(cd_oxpoint .and. lcdox==2 .and. mod(nstep,ncd).eq.0 .and. irk.eq.1) call distribution_cd_oxpoint(ef(:,:,3,2))

      if(bootstrap) then
      call current_boot(lbs)      
      do m=1,3
      ef(:,:,:,m)=ef(:,:,:,m)-eta(:,:,:)*cub(:,:,:,m)
      enddo
      endif
      
      if(lcd .eq. 1) then
      call current_driven
      do m=1,3
      ef(:,:,:,m)=ef(:,:,:,m)-eta(:,:,:)*cud(:,:,:,m)
      enddo
      endif

      if(nstep.lt.nper) then
      do 11 m=1,3
      do 11 jy=iy_first,iy_last
      do 11 jz=iz_first,iz_last
      do 11 jx=ix_first,ix_last
      ef(jx,jz,jy,m)=ef(jx,jz,jy,m)+eta1(jx,jz,jy)*(cint(jx,jz,m)+cur(jx,jz,jy,m))
  11  continue
      endif 




!      call ef_atlastgrid_r1p0_v1(3)

!      if(resisitive) then
!      if(etaj_in_e) then
!      do 2 jy=1,my
!      do 2 jz=iz_first,iz_last
!      do 2 jx=ix_first,ix_last
!      ef(jx,jz,jy,1)=ef(jx,jz,jy,1)+eta(jx,jz,jy)*(cur(jx,jz,jy,1))
!      ef(jx,jz,jy,2)=ef(jx,jz,jy,2)+eta(jx,jz,jy)*(cur(jx,jz,jy,2))
!      ef(jx,jz,jy,3)=ef(jx,jz,jy,3)+eta(jx,jz,jy)*(cur(jx,jz,jy,3))
!   2  continue
!      endif
!      endif

!      write(*,*) 'ef'
       
!      call convte
!      call recrd_ef
       !***************revised**************************************
       call mpi_transfersm(ef(:,:,:,:),3)
    ! haowei need to obatin the efield for bnd grids, i.e., ef_3bndz, ef_3bndx

    do 3 jz=iz_first,iz_last
    do 3 jx=ix_first,ix_last
    itag=gdtp_ep(jx,jz,6)
    if((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz,4).eq.5)) ef(jx,jz,:,:)=0.d0

    if((gdtp_ep(jx,jz,5).eq.1).or.(gdtp_ep(jx,jz,5).eq.-1)) then ! jx+-1 is the boundary point, use the jx, jx-1, jx-2 to calculate the bndz point

    ef_3bndz(itag,:,:)=0.d0

    else if((gdtp_ep(jx,jz,4).eq.5).and.(gdtp_ep(jx,jz,1).ne.5).and.(gdtp_ep(jx-1,jz,5).eq.1).and.(jx.gt.ix_first+1)) then ! jx is the inside dropped point
    
    do 22 m=1,3
    do 22 jy=1,my
    call interp1d2l(ef(jx-2,jz,jy,m),ef(jx-1,jz,jy,m),ef_3bndz(itag,jy,m), &
!    call interp1d2l(ef(jx-2,jz,jy,m),ef(jx-1,jz,jy,m),0.d0, &
            xx(jx-2),xx(jx-1),bndz_grd(itag,1),xx(jx),ef(jx,jz,jy,m))
  22 continue

    else if((gdtp_ep(jx,jz,4).eq.5).and.(gdtp_ep(jx,jz,1).ne.5).and.(gdtp_ep(jx+1,jz,5).eq.-1).and.(jx.lt.ix_last-1)) then ! jx is the inside dropped point

    do 24 m=1,3
    do 24 jy=1,my
    call interp1d2l(ef(jx+2,jz,jy,m),ef(jx+1,jz,jy,m),ef_3bndz(itag,jy,m), &
!    call interp1d2l(ef(jx+2,jz,jy,m),ef(jx+1,jz,jy,m),0.d0, &
            xx(jx+2),xx(jx+1),bndz_grd(itag,1),xx(jx),ef(jx,jz,jy,m))
  24 continue
    endif

!
    itag=gdtp_ep(jx,jz,3)
    if((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz,4).eq.5)) ef(jx,jz,:,:)=0.d0

    if((gdtp_ep(jx,jz,2).eq.1).or.(gdtp_ep(jx,jz,2).eq.-1)) then ! jz+1 is the boundary point, use the jz, jz-1, jz-2 to calculate the bndx point
    ef_3bndx(itag,:,:)=0.d0

    else if((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz,4).ne.5).and.(gdtp_ep(jx,jz-1,2).eq.1).and.(jz.gt.iz_first+1)) then ! jz is the inside dropped point
    do 26 m=1,3
    do 26 jy=1,my
    call interp1d2l(ef(jx,jz-2,jy,m),ef(jx,jz-1,jy,m),ef_3bndx(itag,jy,m), &
!    call interp1d2l(ef(jx,jz-2,jy,m),ef(jx,jz-1,jy,m),0.d0, &
            zz(jz-2),zz(jz-1),bndx_grd(itag,2),zz(jz),ef(jx,jz,jy,m))
  26 continue

    else if((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz,4).ne.5).and.(gdtp_ep(jx,jz+1,2).eq.-1).and.(jz.lt.iz_last-1)) then ! jz is the inside dropped point
    do 28 m=1,3
    do 28 jy=1,my
    call interp1d2l(ef(jx,jz+2,jy,m),ef(jx,jz+1,jy,m),ef_3bndx(itag,jy,m), &
!    call interp1d2l(ef(jx,jz+2,jy,m),ef(jx,jz+1,jy,m),0.d0, &
            zz(jz+2),zz(jz+1),bndx_grd(itag,2),zz(jz),ef(jx,jz,jy,m))
  28 continue
    endif

    3 continue

        if(nstep.eq.0) ef(:,:,:,:)=0.d0
	do 10 m=1,3
      do 10 jy=1,my
      do 10 jz=iz_first,iz_last
      do 10 jx=ix_first,ix_last
!      if(psi(jx,jz).lt.psia1) then
!	if((hypb_ratio(jx,jz).ge.0.d0).and.(hypb_ratio(jx,jz).le.1.d0)) then
	ef(jx,jz,jy,m)=ef(jx,jz,jy,m)*hypb_ratio(jx,jz)
!	endif
!	endif
   10 continue
!       call mpi_transfersm(ef(:,:,:,:),3)
       !**********************************************************
      if(smoothef) call smthef_dis_v2(3)

!	call smth_irpt_with_difc(4,2,1,3,0.9d0)
!	call smth_irpt_with_difc_v2(4,2,1,3,0.9d0)

      return
      end
!hw*****************************************************************
	
	
	

!hwzhang*****************************************************************************
!hwzhang*****************************************************************************
!subroutines for openacc method
!hwzhang*****************************************************************************
!hwzhang*****************************************************************************	
	subroutine right_openacc
	use declare
	integer itag
	real*8 drvx_dx,drey_dx,dez_dx,dex_dy,dez_dy,dex_dz,dey_dz
	real*8, dimension(mx,mz,my) :: rvx,rey, drvx_dx_tmp, dex_dy_tmp, dez_dy_tmp
	real*8, dimension(nbndx,my) :: rvx_bndx, rey_bndx
	real*8, dimension(nbndz,my) :: rvx_bndz, rey_bndz
	real*8 d1fc, d1f2, d1fm, d1fbp ,d1fbm, d1xf2
	real*8 t_start, t_end
	include 'mpif.h'
! R R: calculated as subroutine right
! R IR2: the IR2 part need to be calculated with the boundary point by d1fc, order 4th
! R IR1: the IR1 part need to be calculated with the boundary point by d1fbm, d1fbp, order 3th
! IR2 IR2: both two IR2 parts need to be calculated with the boundary point by d1fc, order 4th
! IR1 IR1: both tow IR1 parts need to be calculated with the boundary point by d1fbm, d1fbp, order 3th
! IR1 IR2: the IR1 part need to be calculated with the boundary point by d1fbm, d1fbp, order 3th
! 	     the IR2 parts need to be calculated with the boundary point by d1fc, order 4th
! R or IR1 or IR2 with D: interpolation in the dropped direction
! D D: do not calculate
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
		((xm2-x0)/(xm1-x0)*(fm1-f0) &
		-(xm1-x0)/(xm2-x0)*(fm2-f0))/(xm2-xm1)
!  d1fbp= d f / dx  with  one-sided-bias  difference involving -1 0  1 and 2
!  points
	d1fbp(fp2,fp1,f0,fm1,a,b,c)= &
		a*(fm1-f0)+b*(fp1-f0)+c*(fp2-f0)
!  d1fbm= d f / dx  with  one-sided-bias  difference involving -2 -1  0 and 1
!  points
	d1fbm(fm2,fm1,f0,fp1,a,b,c)= &
		a*(f0-fp1)+b*(f0-fm1)+c*(f0-fm2)
!  d1xf2= d rf / dx  with second-order accuracy central difference
	d1xf2(fm1,f0,fp1,xm1,x0,xp1)= &
		((xm1-x0)/(xp1-x0)*(xp1*fp1-x0*f0) &
		-(xp1-x0)/(xm1-x0)*(xm1*fm1-x0*f0))/(xm1-xp1)

	call bndry8_cut_cell_v2_fixed
!	print*,'success 0', nrank
	
	if(rho_from_p) then
	do 11 jy=iy_first,iy_last
	do 11 jz=iz_first,iz_last      
	do 11 jx=ix_first,ix_last
!      if(psi(jx,jz).lt.psia1 .and. x(jx,jz,jy,2).gt.0) then
      if( gdtp_ep(jx,jz,1).ne.4 .and. x(jx,jz,jy,2).gt.0) then ! gdtp_ep=4 means the outside point
      x(jx,jz,jy,1)=xint(jx,jz,1)*(xint(jx,jz,2)/x(jx,jz,jy,2))**gamma
      else
      x(jx,jz,jy,1)=xint(jx,jz,1)
      endif
   11 continue
	endif	   

!	print*,'success 1', nrank
	call convt_openacc
!	print*,'success 2', nrank
      call current_openacc
!	print*,'success 3', nrank
      call efield_openacc ! default if(ef_mode) 
!	print*,'success 4', nrank
 	
      !$acc parallel loop &
      !$acc local(jx,jz,jy,itag) &
      !$acc copyout(rvx, rey, rvx_bndz, rvx_bndx, rey_bndz, rey_bndx) &
      !$acc copyin(xx,x)

      !acc copyin(xx,x(*,*,*,3),ef(*,*,*,2),bndz_grd(*,1),bndx_grd(*,1),x_8bndz(*,*,3),ef_3bndz(*,*,2),x_8bndx(*,*,3),ef_3bndx(*,*,2),iy_first,iy_last,iz_first,iz_last,ix_first,ix_last,nbndx,nbndz)
      !acc annotate(entire(xx)) & 
      !acc annotate(dimension(bndx_grd(nbndx,7),bndz_grd(nbndz,7),ef_3bndz(nbndz,my,3),ef_3bndx(nbndx,my,3),x_8bndz(nbndz,my,8),x_8bndx(nbndx,my,8)))
      do jy=iy_first+2,iy_last-2
      do jz=iz_first,iz_last      
      do jx=ix_first,ix_last
      rvx(jx,jz,jy)=xx(jx)*x(jx,jz,jy,3)
	rey(jx,jz,jy)=xx(jx)*ef(jx,jz,jy,2)
      enddo
      enddo
	!calculate the rvx_bndz and rey_bndz
	do itag=1,nbndz
	rvx_bndz(itag,jy)=bndz_grd(itag,1)*x_8bndz(itag,jy,3)
	rey_bndz(itag,jy)=bndz_grd(itag,1)*ef_3bndz(itag,jy,2)	
	enddo
	!calculate the rvx_bndx and rey_bndx
	do itag=1,nbndx
	rvx_bndx(itag,jy)=bndx_grd(itag,1)*x_8bndx(itag,jy,3)
	rey_bndx(itag,jy)=bndx_grd(itag,1)*ef_3bndx(itag,jy,2)
	enddo
      enddo
      !$acc end parallel loop

!	print*,'success 5', nrank

! total(mx*mz*15+500+mx*5+mz*4+3)=6683
      !$acc parallel loop &
      !$acc local(jx,jz,jy,itag) &
      !$acc copyout(drvx_dx_tmp) &
      !$acc copyin(rvx,gdtp_ep(*,*,4:6),rvx_bndz,ax1_irz,bx1_irz,cx1_irz,dx1_irz,axbp_irz,bxbp_irz,cxbp_irz,axbm_irz,bxbm_irz,cxbm_irz,xx,ax1,bx1,cx1,dx1,iy_first,iy_last,ix_first,ix_last,iz_first,iz_last) &
      !$acc annotate(entire(ax1,bx1,cx1,dx1,cxbm_irz, bxbm_irz, axbm_irz, gdtp_ep, cxbp_irz, bxbp_irz, axbp_irz, dx1_irz, cx1_irz, bx1_irz, ax1_irz))
      do jy=iy_first+2,iy_last-2
      do jz=iz_first+2,iz_last-2      
      do jx=ix_first+2,ix_last-2
	if(gdtp_ep(jx,jz,4).eq.2) then 
		  itag=gdtp_ep(jx,jz,6)
		  if(gdtp_ep(jx,jz,5).eq.-2) then ! x1(jx-2,jz,jy,m)=bndz
			    drvx_dx_tmp(jx,jz,jy)=d1fc(rvx_bndz(itag,jy),rvx(jx-1,jz,jy),rvx(jx,jz,jy)&
					,rvx(jx+1,jz,jy),rvx(jx+2,jz,jy),ax1_irz(jx,jz),bx1_irz(jx,jz)&
					,cx1_irz(jx,jz),dx1_irz(jx,jz))
		  endif
		  if(gdtp_ep(jx,jz,5).eq.2) then ! x1(jx+2,jz,jy,m)=bndz
			    drvx_dx_tmp(jx,jz,jy)=d1fc(rvx(jx-2,jz,jy),rvx(jx-1,jz,jy),rvx(jx,jz,jy)&
					,rvx(jx+1,jz,jy),rvx_bndz(itag,jy),ax1_irz(jx,jz),bx1_irz(jx,jz)&
					,cx1_irz(jx,jz),dx1_irz(jx,jz))
		  endif
		  if(gdtp_ep(jx,jz,5).eq.-1) then !x1(jx-1,jz,jy,m)=bndz d1fbp d2fbp
			    drvx_dx_tmp(jx,jz,jy)=d1fbp(rvx(jx+2,jz,jy),rvx(jx+1,jz,jy),rvx(jx,jz,jy)&
					,rvx_bndz(itag,jy),axbp_irz(jx,jz),bxbp_irz(jx,jz),cxbp_irz(jx,jz))
		  endif
		  if(gdtp_ep(jx,jz,5).eq.1) then !x1(jx+1,jz,jy,m)=bndz d1fbm d2fbm
			    drvx_dx_tmp(jx,jz,jy)=d1fbm(rvx(jx-2,jz,jy),rvx(jx-1,jz,jy),rvx(jx,jz,jy)&
					,rvx_bndz(itag,jy),axbm_irz(jx,jz),bxbm_irz(jx,jz),cxbm_irz(jx,jz))
		  endif
	else ! for no IR points, use regular d1fc
	drvx_dx_tmp(jx,jz,jy)=d1fc(rvx(jx-2,jz,jy),rvx(jx-1,jz,jy),rvx(jx,jz,jy),rvx(jx+1,jz,jy),rvx(jx+2,jz,jy),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
	endif
	enddo
	enddo
	enddo
      !$acc end parallel loop
!	print*,'success 6', nrank

! we choose mx=mz=20, (mx=mz=256/16+4)
! each core has 64K LDM, requires less than numbers of 8192 real*8 data is copyed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!														!
!---------------------------- Eq1, ddt\rho=R.H.S -----------------------------------!
!														!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! first calculate the ddt(\rho)
! total (mx*mz*18+mx*5+mz*4+3)=7383
	!$acc parallel loop &
	!$acc local(jx,jz,jy) &
	!$acc copyout(xdif) &
	!$acc copyin(drvx_dx_tmp,x1r,x1z,x,x1,xint3,xint5,xint1_dx,xint3_dx,xint1_dz,xint5_dz,xy,xx,ax1,bx1,cx1,dx1,az1,bz1,cz1,dz1,ix_first,ix_last,iz_first,iz_last,iy_first,iy_last) &
	!$acc annotate(entire(xx,ax1,bx1,cx1,dx1,az1,bz1,cz1,dz1,xint3,xint5,xint1_dx,xint3_dx,xint1_dz,xint5_dz)) &
	!$acc annotate(slice(x1r(*,*,jy,1),x1z(*,*,jy,1),x1z(*,*,jy,5),x(*,*,jy,1),x(*,*,jy,4),x1(*,*,jy,1),x1(*,*,jy,3),x1(*,*,jy,5),xy(*,*,jy,1),xy(*,*,jy,4),xdif(*,*,jy,1)))
	do 1 jy=iy_first+2,iy_last-2
	do 1 jz=iz_first+2,iz_last-2      
	do 1 jx=ix_first+2,ix_last-2
	
      xdif(jx,jz,jy,1)=-x1(jx,jz,jy,3)*(x1r(jx,jz,jy,1)+xint1_dx(jx,jz)) &
       -x(jx,jz,jy,1)*drvx_dx_tmp(jx,jz,jy)/xx(jx) &
!rplc       -x(jx,jz,jy,3)*x(jx,jz,jy,1)/xx(jx)-x(jx,jz,jy,1)*xr(jx,jz,jy,3) &
       -(xy(jx,jz,jy,1)*x(jx,jz,jy,4)+x(jx,jz,jy,1)*xy(jx,jz,jy,4))/xx(jx) &
       -x1(jx,jz,jy,5)*(x1z(jx,jz,jy,1)+xint1_dz(jx,jz))-x(jx,jz,jy,1)*x1z(jx,jz,jy,5) &
!ws:for poloidal flow
       -xint3(jx,jz)*x1r(jx,jz,jy,1)-x1(jx,jz,jy,1)*(xint3(jx,jz)/xx(jx)+xint3_dx(jx,jz)) &
       -xint5(jx,jz)*x1z(jx,jz,jy,1)-x1(jx,jz,jy,1)*xint5_dz(jx,jz)	
		
    1 continue	
      !$acc end parallel loop
!	print*,'success 7', nrank
    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!														!
!------------------------------- Eq2, ddtP=R.H.S -----------------------------------!
!														!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    	
! second calculate the ddt(p)
! total(mx*mz*(14+1)+mx*5+mz*4+3)=7383
	!$acc parallel loop &
	!$acc local(jx,jz,jy) &
	!$acc copyout(xdif) &
	!$acc copyin(drvx_dx_tmp,x1r,x1z,x,x1,xint3,xint5,xint2_dx,xint3_dx,xint2_dz,xint5_dz,xy,xx,ix_first,ix_last,iz_first,iz_last,iy_first,iy_last,gamma) &
	!$acc annotate(entire(xx,xint3,xint5,xint2_dx,xint3_dx,xint2_dz,xint5_dz)) &
	!$acc annotate(slice(x1r(*,*,jy,2),x1z(*,*,jy,2),x1z(*,*,jy,5),x(*,*,jy,2),x(*,*,jy,4),x1(*,*,jy,2),x1(*,*,jy,3),x1(*,*,jy,5),xy(*,*,jy,2),xy(*,*,jy,4),xdif(*,*,jy,2)))
	do 2 jy=iy_first+2,iy_last-2
	do 2 jz=iz_first+2,iz_last-2      
	do 2 jx=ix_first+2,ix_last-2

      xdif(jx,jz,jy,2)=-x1(jx,jz,jy,3)*(xint2_dx(jx,jz)+x1r(jx,jz,jy,2)) &
		-gamma*x(jx,jz,jy,2)*(drvx_dx_tmp(jx,jz,jy)/xx(jx)) &
!rplc       -gamma*x(jx,jz,jy,2)*(xr(jx,jz,jy,3)+x(jx,jz,jy,3)/xx(jx)) &
		-(xy(jx,jz,jy,2)*x(jx,jz,jy,4)+gamma*x(jx,jz,jy,2)*xy(jx,jz,jy,4))/xx(jx) &
		-x1(jx,jz,jy,5)*(xint2_dz(jx,jz)+x1z(jx,jz,jy,2))-gamma*x(jx,jz,jy,2)*x1z(jx,jz,jy,5) &    
!ws:for poloidal flow
		-xint3(jx,jz)*x1r(jx,jz,jy,2)-gamma*x1(jx,jz,jy,2)*(xint3(jx,jz)/xx(jx)+xint3_dx(jx,jz)) &
		-xint5(jx,jz)*x1z(jx,jz,jy,2)-gamma*x1(jx,jz,jy,2)*xint5_dz(jx,jz)	
	
    2 continue
      !$acc end parallel loop
    
!	print*,'success 8', nrank
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!														!
!------------------------------ Eq3, ddtVx=R.H.S, part1-----------------------------!
!														!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! calculate ddtvx, part1
! total(mx*mz*17+mx*5+mz*4+3)=7383
	!$acc parallel loop &
	!$acc local(jx,jz,jy) &
	!$acc copyout(xdif) &
	!$acc copyin(x1r,x1z,x1,x,xy,cur,cint2,cint3,xx,ix_first,ix_last,iz_first,iz_last,iy_first,iy_last) &
	!$acc annotate(entire(xx,cint2,cint3)) &
	!$acc annotate(slice(x1r(*,*,jy,2),x1r(*,*,jy,3),x1z(*,*,jy,3),x(*,*,jy,3),x(*,*,jy,4),x(*,*,jy,5),x(*,*,jy,8),x(*,*,jy,8),x(*,*,jy,1),x1(*,*,jy,4),x1(*,*,jy,7),x1(*,*,jy,8),xy(*,*,jy,3),cur(*,*,jy,2),cur(*,*,jy,3),xdif(*,*,jy,3)))
	do 3 jy=iy_first+2,iy_last-2
	do 3 jz=iz_first+2,iz_last-2      
	do 3 jx=ix_first+2,ix_last-2
	
      xdif(jx,jz,jy,3)=-x(jx,jz,jy,3)*x1r(jx,jz,jy,3)-x(jx,jz,jy,5)*x1z(jx,jz,jy,3) &
		-x(jx,jz,jy,4)*xy(jx,jz,jy,3)/xx(jx)+x(jx,jz,jy,4)*x1(jx,jz,jy,4)/xx(jx) &
		+(cur(jx,jz,jy,2)*x(jx,jz,jy,8)+cint2(jx,jz)*x1(jx,jz,jy,8) &
		- cur(jx,jz,jy,3)*x(jx,jz,jy,7)-cint3(jx,jz)*x1(jx,jz,jy,7) &
		-x1r(jx,jz,jy,2))/x(jx,jz,jy,1) 
    3 continue
      !$acc end parallel loop
!	print*,'success 9', nrank

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!														!
!------------------------------ Eq4, ddtVy=R.H.S, part1-----------------------------!
!														!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
! calculate ddtvy, part1
! total(mx*mz*17+mx*5+mz*4+3)=7383    
	!$acc parallel loop &
	!$acc local(jx,jz,jy) &
	!$acc copyout(xdif) &
	!$acc copyin(x1r,x1z,x1,x,xy,cur,cint1,cint3,xx,ix_first,ix_last,iz_first,iz_last,iy_first,iy_last) &
	!$acc annotate(entire(xx,cint1,cint3)) &
	!$acc annotate(slice(x1r(*,*,jy,4),x1z(*,*,jy,4),x(*,*,jy,1),x(*,*,jy,3),x(*,*,jy,4),x(*,*,jy,5),x(*,*,jy,6),x(*,*,jy,8),x(*,*,jy,1),x1(*,*,jy,3),x1(*,*,jy,6),x1(*,*,jy,8),xy(*,*,jy,2),xy(*,*,jy,4),cur(*,*,jy,1),cur(*,*,jy,3),xdif(*,*,jy,4)))
	do 4 jy=iy_first+2,iy_last-2 
	do 4 jz=iz_first+2,iz_last-2      
	do 4 jx=ix_first+2,ix_last-2
	
	xdif(jx,jz,jy,4)=-x(jx,jz,jy,3)*x1r(jx,jz,jy,4)-x(jx,jz,jy,5)*x1z(jx,jz,jy,4) &
		-x(jx,jz,jy,4)*xy(jx,jz,jy,4)/xx(jx)-x(jx,jz,jy,4)*x1(jx,jz,jy,3)/xx(jx) &
		+(cur(jx,jz,jy,3)*x(jx,jz,jy,6)+cint3(jx,jz)*x1(jx,jz,jy,6) &
		- cur(jx,jz,jy,1)*x(jx,jz,jy,8)-cint1(jx,jz)*x1(jx,jz,jy,8) &   
		-xy(jx,jz,jy,2)/xx(jx))/x(jx,jz,jy,1)
    4 continue		
      !$acc end parallel loop

!	print*,'success 10', nrank

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!														!
!------------------------------ Eq5, ddtVz=R.H.S, part1-----------------------------!
!														!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate ddtvz, part1
! total(mx*mz*16+mx*5+mz*4+3)=6983
	!$acc parallel loop &
	!$acc local(jx,jz,jy) &
	!$acc copyout(xdif) &
	!$acc copyin(x1r,x1z,x1,x,xy,cur,cint1,cint2,xx,ix_first,ix_last,iz_first,iz_last,iy_first,iy_last) &
	!$acc annotate(entire(xx,cint1,cint2)) &
	!$acc annotate(slice(x1r(*,*,jy,5),x1z(*,*,jy,2),x1z(*,*,jy,5),x(*,*,jy,3),x(*,*,jy,4),x(*,*,jy,5),x(*,*,jy,6),x(*,*,jy,7),x(*,*,jy,1),x1(*,*,jy,6),x1(*,*,jy,7),xy(*,*,jy,5),cur(*,*,jy,1),cur(*,*,jy,2),xdif(*,*,jy,5)))
	do 5 jy=iy_first+2,iy_last-2
	do 5 jz=iz_first+2,iz_last-2      
	do 5 jx=ix_first+2,ix_last-2    
	
	xdif(jx,jz,jy,5)=-x(jx,jz,jy,3)*x1r(jx,jz,jy,5)-x(jx,jz,jy,5)*x1z(jx,jz,jy,5) &
		-x(jx,jz,jy,4)*xy(jx,jz,jy,5)/xx(jx) &     
		+(cur(jx,jz,jy,1)*x(jx,jz,jy,7)+cint1(jx,jz)*x1(jx,jz,jy,7) &
		- cur(jx,jz,jy,2)*x(jx,jz,jy,6)-cint2(jx,jz)*x1(jx,jz,jy,6) & 
		-x1z(jx,jz,jy,2))/x(jx,jz,jy,1)
		
    5 continue		
      !$acc end parallel loop
!	print*,'success 11', nrank
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!														!
!------------------------------ Eq345, ddtV=R.H.S, part2~---------------------------!
!														!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 	
! calculate ddtvx, part2
! total(mx*mz*17+mx*1)=6820
	!$acc parallel loop &
	!$acc local(jx,jz,jy) &
	!$acc copy(xdif) &
	!$acc copyin(x1,x,xint3,xint4,xint5,xint3_dx,xint4_dx,xint5_dx,xint3_dz,xint4_dz,xint5_dz,xx,ix_first,ix_last,iz_first,iz_last,iy_first,iy_last) &
	!$acc annotate(entire(xx,xint3,xint4,xint5,xint3_dx,xint4_dx,xint5_dx,xint3_dz,xint4_dz,xint5_dz)) &
	!$acc annotate(slice(x(*,*,jy,1),x1(*,*,jy,1),x1(*,*,jy,3),x1(*,*,jy,4),x1(*,*,jy,5),xdif(*,*,jy,3),xdif(*,*,jy,4),xdif(*,*,jy,5)))

	do 345 jy=iy_first+2,iy_last-2
	do 345 jz=iz_first+2,iz_last-2      
	do 345 jx=ix_first+2,ix_last-2
	
       xdif(jx,jz,jy,3)=xdif(jx,jz,jy,3)-x1(jx,jz,jy,3)*xint3_dx(jx,jz)-x1(jx,jz,jy,5)*xint3_dz(jx,jz) &
       +x1(jx,jz,jy,4)*xint4(jx,jz)/xx(jx) &
       +(-xint3(jx,jz)*xint3_dx(jx,jz)-xint5(jx,jz)*xint3_dz(jx,jz) &
       +xint4(jx,jz)*xint4(jx,jz)/xx(jx))*x1(jx,jz,jy,1)/x(jx,jz,jy,1)	
	 
	 xdif(jx,jz,jy,4)=xdif(jx,jz,jy,4)-x1(jx,jz,jy,3)*xint4_dx(jx,jz)-x1(jx,jz,jy,5)*xint4_dz(jx,jz) &
		-x1(jx,jz,jy,4)*xint3(jx,jz)/xx(jx) &
		+(-xint3(jx,jz)*xint4_dx(jx,jz)-xint5(jx,jz)*xint4_dz(jx,jz) &
		-xint4(jx,jz)*xint3(jx,jz)/xx(jx))*x1(jx,jz,jy,1)/x(jx,jz,jy,1)
		
	xdif(jx,jz,jy,5)=xdif(jx,jz,jy,5)-x1(jx,jz,jy,3)*xint5_dx(jx,jz)-x1(jx,jz,jy,5)*xint5_dz(jx,jz) &
		+(-xint3(jx,jz)*xint5_dx(jx,jz)-xint5(jx,jz)*xint5_dz(jx,jz)) &
		*x1(jx,jz,jy,1)/x(jx,jz,jy,1)		
	 
  345 continue	
      !$acc end parallel loop
!	print*,'success 12', nrank
    
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!														!
!------------------------------ Eq6-8, ddtB=R.H.S, part1, for R points--------------!
!														!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate ddtB
! use ef_mode as default
	!acc parallel loop &
	!acc local(jx,jz,jy) &
	!acc copyin(ay1,by1,cy1,dy1,iy_first,iy_last,iz_first,iz_last,ix_first,ix_last) &
	!acc annotate(entire(ay1,by1,cy1,dy1))
	do jy=iy_first+2,iy_last-2
	!acc data copyin(ef(*,*,jy-2:jy+2,*)) &
	!acc copyout(dex_dy_tmp(*,*,jy),dez_dy_tmp(*,*,jy)) &
	!acc present(ay1,by1,cy1,dy1,iy_first,iy_last,iz_first,iz_last,ix_first,ix_last,jx,jz,jy)
	do jz=iz_first+2,iz_last-2      
	do jx=ix_first+2,ix_last-2
	dex_dy_tmp(jx,jz,jy) =d1fc(ef(jx,jz,jy-2,1),ef(jx,jz,jy-1,1),ef(jx,jz,jy,1) &
		,ef(jx,jz,jy+1,1),ef(jx,jz,jy+2,1),ay1(jy),by1(jy),cy1(jy),dy1(jy)) 
	dez_dy_tmp(jx,jz,jy) =d1fc(ef(jx,jz,jy-2,3),ef(jx,jz,jy-1,3),ef(jx,jz,jy,3) &
		,ef(jx,jz,jy+1,3),ef(jx,jz,jy+2,3),ay1(jy),by1(jy),cy1(jy),dy1(jy))
	enddo
	enddo
      !acc end data
	enddo
      !acc end parallel loop

! total(mx*mz*8+mx*5+mz*4)=3380
	!$acc parallel loop &
	!$acc local(jx,jz,jy,drey_dx,dez_dx,dex_dz,dey_dz) &
	!$acc copyout(xdif) &
	!$acc copyin(dex_dy_tmp,dez_dy_tmp,ef,xx,ax1,bx1,cx1,dx1,az1,bz1,cz1,dz1,ix_first,ix_last,iz_first,iz_last,iy_first,iy_last) &
	!$acc annotate(entire(xx,ax1,bx1,cx1,dx1,az1,bz1,cz1,dz1)) &
	!$acc annotate(slice(xdif(*,*,jy,6),xdif(*,*,jy,7),xdif(*,*,jy,8)))
	do 7 jy=iy_first+2,iy_last-2
	do 7 jz=iz_first+2,iz_last-2
	do 7 jx=ix_first+2,ix_last-2

	drey_dx=d1fc(xx(jx)*ef(jx-2,jz,jy,2),xx(jx)*ef(jx-1,jz,jy,2),xx(jx)*ef(jx,jz,jy,2) &
		,xx(jx)*ef(jx+1,jz,jy,2),xx(jx)*ef(jx+2,jz,jy,2),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
	dez_dx =d1fc(ef(jx-2,jz,jy,3),ef(jx-1,jz,jy,3),ef(jx,jz,jy,3) &
		,ef(jx+1,jz,jy,3),ef(jx+2,jz,jy,3),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
	dex_dz =d1fc(ef(jx,jz-2,jy,1),ef(jx,jz-1,jy,1),ef(jx,jz,jy,1) &
		,ef(jx,jz+1,jy,1),ef(jx,jz+2,jy,1),az1(jz),bz1(jz),cz1(jz),dz1(jz))
	dey_dz =d1fc(ef(jx,jz-2,jy,2),ef(jx,jz-1,jy,2),ef(jx,jz,jy,2) &
		,ef(jx,jz+1,jy,2),ef(jx,jz+2,jy,2),az1(jz),bz1(jz),cz1(jz),dz1(jz))
!	dex_dy =d1fc(ef(jx,jz,jy-2,1),ef(jx,jz,jy-1,1),ef(jx,jz,jy,1) &
!		,ef(jx,jz,jy+1,1),ef(jx,jz,jy+2,1),ay1(jy),by1(jy),cy1(jy),dy1(jy)) 
!	dez_dy =d1fc(ef(jx,jz,jy-2,3),ef(jx,jz,jy-1,3),ef(jx,jz,jy,3) &
!		,ef(jx,jz,jy+1,3),ef(jx,jz,jy+2,3),ay1(jy),by1(jy),cy1(jy),dy1(jy)) 	
!	dex_dy =dex_dy_tmp(jx,jz,jy)
!	dez_dy =dez_dy_tmp(jx,jz,jy)
	
      xdif(jx,jz,jy,6)=-dez_dy_tmp(jx,jz,jy)/xx(jx)+dey_dz 
      xdif(jx,jz,jy,7)=-dex_dz+dez_dx 
      xdif(jx,jz,jy,8)=(dex_dy_tmp(jx,jz,jy)-drey_dx)/xx(jx) 
	
    7 continue
      !$acc end parallel loop
    
!	print*,'success 13', nrank
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!														!
!------------------------------ Eq6-8, ddtB=R.H.S, part2, for IR points-------------!
!														!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
! hard to use openacc
! add a subroutine to decide the jx and jz region for IR points
	!$acc parallel loop &
	!$acc local(jx,jz,jy,drey_dx,dez_dx,dex_dz,itag) &
	!$acc copy(xdif) &
	!$acc copyin(dex_dy_tmp,ef,gdtp_ep,rey,ax1_irz,bx1_irz,cx1_irz,dx1_irz,axbp_irz,bxbp_irz,cxbp_irz,axbm_irz,bxbm_irz,cxbm_irz,xx,az1,bz1,cz1,dz1,ix_first_irpt2,ix_last_irpt2,iz_first_irpt2,iz_last_irpt2,iy_first,iy_last) &
	!$acc annotate(entire(xx,az1,bz1,cz1,dz1,gdtp_ep,cxbm_irz,bxbm_irz,axbm_irz,cxbp_irz,bxbp_irz,axbp_irz,dx1_irz,cx1_irz,bx1_irz,ax1_irz)) &
	!$acc annotate(slice(ef(*,*,jy,1),ef(*,*,jy,3),xdif(*,*,jy,7),xdif(*,*,jy,8)))
      !acc annotate(dimension(rey_bndz(nbndz,my),ef_3bndz(nbndz,my,3))) &
	do 71 jy=iy_first+2,iy_last-2
	do 71 jz=iz_first_irpt2,iz_last_irpt2
	do 71 jx=ix_first_irpt2,ix_last_irpt2
!	do 71 jz=iz_first+2,iz_last-2
!	do 71 jx=ix_first+2,ix_last-2

	! hw for IR in x direction, drey_dx dez_dx need to be recalculated
	if(gdtp_ep(jx,jz,4).eq.2) then 
	
		dex_dz =d1fc(ef(jx,jz-2,jy,1),ef(jx,jz-1,jy,1),ef(jx,jz,jy,1) &
			,ef(jx,jz+1,jy,1),ef(jx,jz+2,jy,1),az1(jz),bz1(jz),cz1(jz),dz1(jz))
!		dey_dz =d1fc(ef(jx,jz-2,jy,2),ef(jx,jz-1,jy,2),ef(jx,jz,jy,2) &
!			,ef(jx,jz+1,jy,2),ef(jx,jz+2,jy,2),az1(jz),bz1(jz),cz1(jz),dz1(jz))
	
		  itag=gdtp_ep(jx,jz,6)
		  if(gdtp_ep(jx,jz,5).eq.-2) then ! x1(jx-2,jz,jy,m)=bndz
!			    drey_dx=d1fc(rey_bndz(itag,jy),rey(jx-1,jz,jy),rey(jx,jz,jy)&
			    drey_dx=d1fc(0.d0,rey(jx-1,jz,jy),rey(jx,jz,jy)&
					,rey(jx+1,jz,jy),rey(jx+2,jz,jy),ax1_irz(jx,jz),bx1_irz(jx,jz)&
					,cx1_irz(jx,jz),dx1_irz(jx,jz))
!			    dez_dx=d1fc(ef_3bndz(itag,jy,3),ef(jx-1,jz,jy,3),ef(jx,jz,jy,3)&
			    dez_dx=d1fc(0.d0,ef(jx-1,jz,jy,3),ef(jx,jz,jy,3)&
					,ef(jx+1,jz,jy,3),ef(jx+2,jz,jy,3),ax1_irz(jx,jz),bx1_irz(jx,jz)&
					,cx1_irz(jx,jz),dx1_irz(jx,jz))
		  endif
		  if(gdtp_ep(jx,jz,5).eq.2) then ! x1(jx+2,jz,jy,m)=bndz
			    drey_dx=d1fc(rey(jx-2,jz,jy),rey(jx-1,jz,jy),rey(jx,jz,jy)&
!					,rey(jx+1,jz,jy),rey_bndz(itag,jy),ax1_irz(jx,jz),bx1_irz(jx,jz)&
					,rey(jx+1,jz,jy),0.d0,ax1_irz(jx,jz),bx1_irz(jx,jz)&
					,cx1_irz(jx,jz),dx1_irz(jx,jz))
			    dez_dx=d1fc(ef(jx-2,jz,jy,3),ef(jx-1,jz,jy,3),ef(jx,jz,jy,3)&
!					,ef(jx+1,jz,jy,3),ef_3bndz(itag,jy,3),ax1_irz(jx,jz),bx1_irz(jx,jz)&
					,ef(jx+1,jz,jy,3),0.d0,ax1_irz(jx,jz),bx1_irz(jx,jz)&
					,cx1_irz(jx,jz),dx1_irz(jx,jz))
		  endif
		  if(gdtp_ep(jx,jz,5).eq.-1) then !x1(jx-1,jz,jy,m)=bndz d1fbp d2fbp
			    drey_dx=d1fbp(rey(jx+2,jz,jy),rey(jx+1,jz,jy),rey(jx,jz,jy)&
!					,rey_bndz(itag,jy),axbp_irz(jx,jz),bxbp_irz(jx,jz),cxbp_irz(jx,jz))
					,0.d0,axbp_irz(jx,jz),bxbp_irz(jx,jz),cxbp_irz(jx,jz))
			    dez_dx=d1fbp(ef(jx+2,jz,jy,3),ef(jx+1,jz,jy,3),ef(jx,jz,jy,3)&
!					,ef_3bndz(itag,jy,3),axbp_irz(jx,jz),bxbp_irz(jx,jz),cxbp_irz(jx,jz))
					,0.d0,axbp_irz(jx,jz),bxbp_irz(jx,jz),cxbp_irz(jx,jz))
		  endif
		  if(gdtp_ep(jx,jz,5).eq.1) then !x1(jx+1,jz,jy,m)=bndz d1fbm d2fbm
			    drey_dx=d1fbm(rey(jx-2,jz,jy),rey(jx-1,jz,jy),rey(jx,jz,jy)&
!					,rey_bndz(itag,jy),axbm_irz(jx,jz),bxbm_irz(jx,jz),cxbm_irz(jx,jz))
					,0.d0,axbm_irz(jx,jz),bxbm_irz(jx,jz),cxbm_irz(jx,jz))
			    dez_dx=d1fbm(ef(jx-2,jz,jy,3),ef(jx-1,jz,jy,3),ef(jx,jz,jy,3)&
!					,ef_3bndz(itag,jy,3),axbm_irz(jx,jz),bxbm_irz(jx,jz),cxbm_irz(jx,jz))
					,0.d0,axbm_irz(jx,jz),bxbm_irz(jx,jz),cxbm_irz(jx,jz))
		  endif
!      xdif(jx,jz,jy,6)=-dez_dy_tmp(jx,jz,jy)/xx(jx)+dey_dz
      xdif(jx,jz,jy,7)=-dex_dz+dez_dx
      xdif(jx,jz,jy,8)=(dex_dy_tmp(jx,jz,jy)-drey_dx)/xx(jx)		  
	endif
   71 continue
      !$acc end parallel loop

!	print*,'success 14', nrank

! hard to use openacc
! add a subroutine to decide the jx and jz region for IR points
	!$acc parallel loop &
	!$acc local(jx,jz,jy,dez_dx,dex_dz,dey_dz,itag) &
	!$acc copy(xdif) &
	!$acc copyin(dez_dy_tmp,ef,gdtp_ep,az1_irx,bz1_irx,cz1_irx,dz1_irx,azbp_irx,bzbp_irx,czbp_irx,azbm_irx,bzbm_irx,czbm_irx,xx,ax1,bx1,cx1,dx1,ix_first_irpt2,ix_last_irpt2,iz_first_irpt2,iz_last_irpt2,iy_first,iy_last) &
	!$acc annotate(entire(xx,ax1,bx1,cx1,dx1,gdtp_ep,czbm_irx,bzbm_irx,azbm_irx,czbp_irx,bzbp_irx,azbp_irx,dz1_irx,cz1_irx,bz1_irx,az1_irx)) &
	!$acc annotate(slice(ef(*,*,jy,1),ef(*,*,jy,2),ef(*,*,jy,3),xdif(*,*,jy,6),xdif(*,*,jy,7)))
      !acc annotate(dimension(ef_3bndx(nbndx,my,3))) $ ef is set to zero at the boundary
	do 72 jy=iy_first+2,iy_last-2
	do 72 jz=iz_first_irpt2,iz_last_irpt2
	do 72 jx=ix_first_irpt2,ix_last_irpt2

	! for IR in z direction, dex_dz dey_dz need to be recalculated
	if(gdtp_ep(jx,jz,1).eq.2) then 
	
!		drey_dx=d1fc(rey(jx-2,jz,jy),rey(jx-1,jz,jy),rey(jx,jz,jy) &
!			,rey(jx+1,jz,jy),rey(jx+2,jz,jy),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
		dez_dx =d1fc(ef(jx-2,jz,jy,3),ef(jx-1,jz,jy,3),ef(jx,jz,jy,3) &
			,ef(jx+1,jz,jy,3),ef(jx+2,jz,jy,3),ax1(jx),bx1(jx),cx1(jx),dx1(jx))	
	
		  itag=gdtp_ep(jx,jz,3)
		  if(gdtp_ep(jx,jz,2).eq.-2) then ! x1(jx,jz-2,jy,m)=bndx
!			    dex_dz=d1fc(ef_3bndx(itag,jy,1),ef(jx,jz-1,jy,1),ef(jx,jz,jy,1)&
			    dex_dz=d1fc(0.d0,ef(jx,jz-1,jy,1),ef(jx,jz,jy,1)&
					,ef(jx,jz+1,jy,1),ef(jx,jz+2,jy,1),az1_irx(jx,jz),bz1_irx(jx,jz)&
					,cz1_irx(jx,jz),dz1_irx(jx,jz))
!			    dey_dz=d1fc(ef_3bndx(itag,jy,2),ef(jx,jz-1,jy,2),ef(jx,jz,jy,2)&
			    dey_dz=d1fc(0.d0,ef(jx,jz-1,jy,2),ef(jx,jz,jy,2)&
					,ef(jx,jz+1,jy,2),ef(jx,jz+2,jy,2),az1_irx(jx,jz),bz1_irx(jx,jz)&
					,cz1_irx(jx,jz),dz1_irx(jx,jz))
		  endif
		  if(gdtp_ep(jx,jz,2).eq.2) then ! x1(jx,jz+2,jy,m)=bndx
			    dex_dz=d1fc(ef(jx,jz-2,jy,1),ef(jx,jz-1,jy,1),ef(jx,jz,jy,1)&
!					,ef(jx,jz+1,jy,1),ef_3bndx(itag,jy,1),az1_irx(jx,jz),bz1_irx(jx,jz)&
					,ef(jx,jz+1,jy,1),0.d0,az1_irx(jx,jz),bz1_irx(jx,jz)&
					,cz1_irx(jx,jz),dz1_irx(jx,jz))
			    dey_dz=d1fc(ef(jx,jz-2,jy,2),ef(jx,jz-1,jy,2),ef(jx,jz,jy,2)&
!					,ef(jx,jz+1,jy,2),ef_3bndx(itag,jy,2),az1_irx(jx,jz),bz1_irx(jx,jz)&
					,ef(jx,jz+1,jy,2),0.d0,az1_irx(jx,jz),bz1_irx(jx,jz)&
					,cz1_irx(jx,jz),dz1_irx(jx,jz))
		  endif
		  if(gdtp_ep(jx,jz,2).eq.-1) then ! x1(jx,jz-1,jy,m)=bndx d1fbp d2fbp
			    dex_dz=d1fbp(ef(jx,jz+2,jy,1),ef(jx,jz+1,jy,1),ef(jx,jz,jy,1)&
!					,ef_3bndx(itag,jy,1),azbp_irx(jx,jz),bzbp_irx(jx,jz),czbp_irx(jx,jz))
					,0.d0,azbp_irx(jx,jz),bzbp_irx(jx,jz),czbp_irx(jx,jz))
			    dey_dz=d1fbp(ef(jx,jz+2,jy,2),ef(jx,jz+1,jy,2),ef(jx,jz,jy,2)&
!					,ef_3bndx(itag,jy,2),azbp_irx(jx,jz),bzbp_irx(jx,jz),czbp_irx(jx,jz))
					,0.d0,azbp_irx(jx,jz),bzbp_irx(jx,jz),czbp_irx(jx,jz))
		  endif
		  if(gdtp_ep(jx,jz,2).eq.1) then ! x1(jx,jz+1,jy,m)=bndx d1fbm d2fbm
			    dex_dz=d1fbm(ef(jx,jz-2,jy,1),ef(jx,jz-1,jy,1),ef(jx,jz,jy,1)&
!					,ef_3bndx(itag,jy,1),azbm_irx(jx,jz),bzbm_irx(jx,jz),czbm_irx(jx,jz))
					,0.d0,azbm_irx(jx,jz),bzbm_irx(jx,jz),czbm_irx(jx,jz))
			    dey_dz=d1fbm(ef(jx,jz-2,jy,2),ef(jx,jz-1,jy,2),ef(jx,jz,jy,2)&
!					,ef_3bndx(itag,jy,2),azbm_irx(jx,jz),bzbm_irx(jx,jz),czbm_irx(jx,jz))
					,0.d0,azbm_irx(jx,jz),bzbm_irx(jx,jz),czbm_irx(jx,jz))
		  endif
	xdif(jx,jz,jy,6)=-dez_dy_tmp(jx,jz,jy)/xx(jx)+dey_dz
      xdif(jx,jz,jy,7)=-dex_dz+dez_dx
!      xdif(jx,jz,jy,8)=(dex_dy_tmp(jx,jz,jy)-drey_dx)/xx(jx)
	endif
   72 continue    
      !$acc end parallel loop
    
!	print*,'success 15', nrank
    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!														!
!------------------------------ Eq1-5, viscous terms--------------------------------!
!														!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    
    
	if(viscous) then     
      if(.not.implicitv)then
!      call vorticity

! total(mx*mz*18+mx+3)=7223
	!$acc parallel loop &
	!$acc local(jx,jz,jy) &
	!$acc copy(xdif) &
	!$acc copyin(xr2,x1r,xz2,x1z,xy2,pmu,pmux,pmuz,kap,kapx,kapz,xx,ix_first,ix_last,iz_first,iz_last,iy_first,iy_last) &
	!$acc annotate(entire(xx,pmu,pmux,pmuz,kap,kapx,kapz)) &
	!$acc annotate(slice(xr2(*,*,jy,1),xr2(*,*,jy,2),x1r(*,*,jy,1),x1r(*,*,jy,2),xz2(*,*,jy,1),xz2(*,*,jy,2),x1z(*,*,jy,1),x1z(*,*,jy,2),xy2(*,*,jy,1),xy2(*,*,jy,2),xdif(*,*,jy,1),xdif(*,*,jy,2)))
      do 8 jy=iy_first+2,iy_last-2
      do 8 jz=iz_first+2,iz_last-2      
      do 8 jx=ix_first+2,ix_last-2
!      if(psi(jx,jz).lt.psia1) then ! keep the viscous term only in closed flux
      !if(gdtp_ep(jx,jz,1).ne.4) then	
      xdif(jx,jz,jy,1)=xdif(jx,jz,jy,1)+pmu(jx,jz) &
       *(xr2(jx,jz,jy,1)+x1r(jx,jz,jy,1)/xx(jx)+xy2(jx,jz,jy,1)/xx(jx)**2+xz2(jx,jz,jy,1)) &
       +pmux(jx,jz)*x1r(jx,jz,jy,1)+pmuz(jx,jz)*x1z(jx,jz,jy,1) 

      xdif(jx,jz,jy,2)=xdif(jx,jz,jy,2)+kap(jx,jz) &
       *(xr2(jx,jz,jy,2)+x1r(jx,jz,jy,2)/xx(jx)+xy2(jx,jz,jy,2)/xx(jx)**2+xz2(jx,jz,jy,2)) &
       +kapx(jx,jz)*x1r(jx,jz,jy,2)+kapz(jx,jz)*x1z(jx,jz,jy,2)
    8 continue
      !$acc end parallel loop
	 
!	print*,'success 16', nrank
	 
! total(mx*mz*19+mx+3)=7603
	!$acc parallel loop &
	!$acc local(jx,jz,jy) &
	!$acc copy(xdif) &
	!$acc copyin(x1,xy,xr2,x1r,xz2,xy2,fmu,xx,ix_first,ix_last,iz_first,iz_last,iy_first,iy_last) &
	!$acc annotate(entire(xx,fmu)) &
	!$acc annotate(slice(x1(*,*,jy,3),x1(*,*,jy,4),xy(*,*,jy,3),xy(*,*,jy,4),xr2(*,*,jy,3),xr2(*,*,jy,4)))
	!$acc annotate(xr2(*,*,jy,5),x1r(*,*,jy,3),x1r(*,*,jy,4),x1r(*,*,jy,5),xz2(*,*,jy,3),xz2(*,*,jy,4),xz2(*,*,jy,5),xy2(*,*,jy,3),xy2(*,*,jy,4),xy2(*,*,jy,5),xdif(*,*,jy,3),xdif(*,*,jy,4),xdif(*,*,jy,5)))
      do 9 jy=iy_first+2,iy_last-2
      do 9 jz=iz_first+2,iz_last-2      
      do 9 jx=ix_first+2,ix_last-2
	
      xdif(jx,jz,jy,3)=xdif(jx,jz,jy,3)+fmu(jx,jz) &
       *(xr2(jx,jz,jy,3)+x1r(jx,jz,jy,3)/xx(jx)+xy2(jx,jz,jy,3)/xx(jx)**2+xz2(jx,jz,jy,3) &
        -x1(jx,jz,jy,3)/xx(jx)**2-2.0*xy(jx,jz,jy,4)/xx(jx)**2)
!       +fmux(jx,jz)*x1r(jx,jz,jy,3)+fmuz(jx,jz)*x1z(jx,jz,jy,3)
     
      xdif(jx,jz,jy,4)=xdif(jx,jz,jy,4)+fmu(jx,jz) &
       *(xr2(jx,jz,jy,4)+x1r(jx,jz,jy,4)/xx(jx)+xy2(jx,jz,jy,4)/xx(jx)**2+xz2(jx,jz,jy,4) &
        -x1(jx,jz,jy,4)/xx(jx)**2+2.0*xy(jx,jz,jy,3)/xx(jx)**2)
!       +fmux(jx,jz)*x1r(jx,jz,jy,4)+fmuz(jx,jz)*x1z(jx,jz,jy,4)
     
      xdif(jx,jz,jy,5)=xdif(jx,jz,jy,5)+fmu(jx,jz) &
       *(xr2(jx,jz,jy,5)+x1r(jx,jz,jy,5)/xx(jx)+xy2(jx,jz,jy,5)/xx(jx)**2+xz2(jx,jz,jy,5))
!       +fmux(jx,jz)*x1r(jx,jz,jy,5)+fmuz(jx,jz)*x1z(jx,jz,jy,5) 
!      endif
    9 continue
      !$acc end parallel loop
!	print*,'success 17', nrank
    
! LDM for openacc is not enough, so seperate it out    
	!$acc parallel loop &
	!$acc local(jx,jz,jy) &
	!$acc copy(xdif) &
	!$acc copyin(fmux,fmuz,x1r,x1z,ix_first,ix_last,iz_first,iz_last,iy_first,iy_last) &
	!$acc annotate(entire(xx,fmux,fmuz)) &
	!$acc annotate(slice(x1r(*,*,jy,3),x1r(*,*,jy,4),x1r(*,*,jy,5),x1z(*,*,jy,3),x1z(*,*,jy,4),x1z(*,*,jy,5),xdif(*,*,jy,3),xdif(*,*,jy,4),xdif(*,*,jy,5)))
	do 10 jy=iy_first+2,iy_last-2
      do 10 jz=iz_first+2,iz_last-2      
      do 10 jx=ix_first+2,ix_last-2
	xdif(jx,jz,jy,3)=xdif(jx,jz,jy,3)+fmux(jx,jz)*x1r(jx,jz,jy,3)+fmuz(jx,jz)*x1z(jx,jz,jy,3) 
	xdif(jx,jz,jy,4)=xdif(jx,jz,jy,4)+fmux(jx,jz)*x1r(jx,jz,jy,4)+fmuz(jx,jz)*x1z(jx,jz,jy,4) 
	xdif(jx,jz,jy,5)=xdif(jx,jz,jy,5)+fmux(jx,jz)*x1r(jx,jz,jy,5)+fmuz(jx,jz)*x1z(jx,jz,jy,5) 
   10 continue
      !$acc end parallel loop
      
!	print*,'success 18', nrank
      endif   
      endif
	
! set the outside point into zero	
!    	do jz=iz_first_irpt2,iz_last_irpt2
!    	do jx=ix_first_irpt2,ix_last_irpt2
	do jz=iz_first+2,iz_last-2
	do jx=ix_first+2,ix_last-2
	if(max(gdtp_ep(jx,jz,1),gdtp_ep(jx,jz,4)).ge.4) xdif(jx,jz,:,:)=0.d0
	enddo
	enddo

      if(invaryrho) xdif(:,:,:,1)=0
      if(invaryp)   xdif(:,:,:,2)=0
      call mpi_transfersm(xdif,8)	
	
!	print*,'success 19', nrank
      return
      endsubroutine right_openacc
	
!

!hw*************************************************************************
      subroutine current_openacc
      use declare
	integer itag
      real*8 drby_dx, rby_tmp
      real*8, dimension(mx,mz,my) :: rby
	real*8 d1f2, d1fc, d1fp, d1fbp, d1fbm, d1xf2
      include 'mpif.h'
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
!  d1xf2= d rf / dx  with second-order accuracy central difference
      d1xf2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(xp1*fp1-x0*f0) &
         -(xp1-x0)/(xm1-x0)*(xm1*fm1-x0*f0))/(xm1-xp1)

!       integer status(mpi_status_size)

      do 10 jy=iy_first,iy_last
      do 10 jz=iz_first,iz_last
      do 10 jx=ix_first,ix_last
      rby(jx,jz,jy)=xx(jx)*x1(jx,jz,jy,7)
   10 continue

! total(mx*mz*9+mx*5+mz*4)=3780
	!$acc parallel loop &
	!$acc local(jx,jz,jy,drby_dx) &
	!$acc copyout(cur) &
	!$acc copyin(xy,x1,rby,xx,ax1,bx1,cx1,dx1,az1,bz1,cz1,dz1,ix_first,ix_last,iz_first,iz_last,iy_first,iy_last) &
	!$acc annotate(entire(xx,ax1,bx1,cx1,dx1,az1,bz1,cz1,dz1)) &
	!$acc annotate(slice(xy(*,*,jy,6),xy(*,*,jy,8),x1(*,*,jy,6),x1(*,*,jy,7),x1(*,*,jy,8)))
      do 1 jy=iy_first+2,iy_last-2
      do 1 jz=iz_first+2,iz_last-2
      do 1 jx=ix_first+2,ix_last-2
	
      cur(jx,jz,jy,1)=xy(jx,jz,jy,8)/xx(jx) &
		-d1fc(x1(jx,jz-2,jy,7),x1(jx,jz-1,jy,7),x1(jx,jz,jy,7) &
		,x1(jx,jz+1,jy,7),x1(jx,jz+2,jy,7),az1(jz),bz1(jz),cz1(jz),dz1(jz))
	   
      cur(jx,jz,jy,2)=d1fc(x1(jx,jz-2,jy,6),x1(jx,jz-1,jy,6),x1(jx,jz,jy,6) &
		,x1(jx,jz+1,jy,6),x1(jx,jz+2,jy,6),az1(jz),bz1(jz),cz1(jz),dz1(jz)) &
		-d1fc(x1(jx-2,jz,jy,8),x1(jx-1,jz,jy,8),x1(jx,jz,jy,8) &
		,x1(jx+1,jz,jy,8),x1(jx+2,jz,jy,8),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
		
	drby_dx=d1fc(rby(jx-2,jz,jy),rby(jx-1,jz,jy),rby(jx,jz,jy) &
		,rby(jx+1,jz,jy),rby(jx+2,jz,jy),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
		
      cur(jx,jz,jy,3)=(drby_dx-xy(jx,jz,jy,6))/xx(jx)

    1 continue
      !$acc end parallel loop
    
    
! current for IR points, calculate without openACC
	!$acc parallel loop &
	!$acc local(jx,jz,jy,itag,rby_tmp,drby_dx) &
	!$acc copy(cur) &
	!$acc copyin(rby,xy,gdtp_ep,ax1_irz,bx1_irz,cx1_irz,dx1_irz,axbp_irz,bxbp_irz,cxbp_irz,axbm_irz,bxbm_irz,cxbm_irz,xx,iy_first,iy_last,iz_first_irpt2,iz_last_irpt2,ix_first_irpt2,ix_last_irpt2,x1_8bndz(*,*,7),bndz_grd(*,1)) &
	!$acc annotate(dimension(x1_8bndz(nbndz,my,8),bndz_grd(nbndz,7))) &
	!$acc annotate(entire(gdtp_ep,ax1_irz,bx1_irz,cx1_irz,dx1_irz,axbp_irz,bxbp_irz,cxbp_irz,axbm_irz,bxbm_irz,cxbm_irz,xx)) &
	!$acc annotate(slice(cur(*,*,jy,3),xy(*,*,jy,6)))
      do 2 jy=iy_first+2,iy_last-2
	do 2 jz=iz_first_irpt2,iz_last_irpt2
	do 2 jx=ix_first_irpt2,ix_last_irpt2
!      do 2 jz=iz_first+2,iz_last-2
!      do 2 jx=ix_first+2,ix_last-2
	
!	if(gdtp_ep(jx,jz,1).eq.2) then
!		drby_dx=d1fc(rby(jx-2,jz,jy),rby(jx-1,jz,jy),rby(jx,jz,jy) &
!			,rby(jx+1,jz,jy),rby(jx+2,jz,jy),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
!		cur(jx,jz,jy,1)=xy(jx,jz,jy,8)/xx(jx)-x1z(jx,jz,jy,7) 
!		cur(jx,jz,jy,2)=x1z(jx,jz,jy,6)-x1r(jx,jz,jy,8)
!		cur(jx,jz,jy,3)=(drby_dx-xy(jx,jz,jy,6))/xx(jx)					
!	endif
	
	if(gdtp_ep(jx,jz,4).eq.2) then
		itag=gdtp_ep(jx,jz,6)
		rby_tmp=bndz_grd(itag,1)*x1_8bndz(itag,jy,7)
		if(gdtp_ep(jx,jz,5).eq.-2) then !jx-2 is the bndz
			drby_dx=d1fc(rby_tmp,rby(jx-1,jz,jy),rby(jx,jz,jy) &
					,rby(jx+1,jz,jy),rby(jx+2,jz,jy) &
					,ax1_irz(jx,jz),bx1_irz(jx,jz),cx1_irz(jx,jz),dx1_irz(jx,jz))
		endif
		if(gdtp_ep(jx,jz,5).eq.2) then !jx+2 is the bndz
			drby_dx=d1fc(rby(jx-2,jz,jy),rby(jx-1,jz,jy),rby(jx,jz,jy) &
					,rby(jx+1,jz,jy),rby_tmp &
					,ax1_irz(jx,jz),bx1_irz(jx,jz),cx1_irz(jx,jz),dx1_irz(jx,jz))
		endif
		if(gdtp_ep(jx,jz,5).eq.-1) then !jx-1 is the bndz
			drby_dx=d1fbp(rby(jx+2,jz,jy),rby(jx+1,jz,jy),rby(jx,jz,jy),rby_tmp &
					,axbp_irz(jx,jz),bxbp_irz(jx,jz),cxbp_irz(jx,jz))
		endif
		if(gdtp_ep(jx,jz,5).eq.1) then !jx+1 is the bndz
			drby_dx=d1fbm(rby(jx-2,jz,jy),rby(jx-1,jz,jy),rby(jx,jz,jy),rby_tmp &
					,axbm_irz(jx,jz),bxbm_irz(jx,jz),cxbm_irz(jx,jz))
		endif
		
!		cur(jx,jz,jy,1)=xy(jx,jz,jy,8)/xx(jx)-x1z(jx,jz,jy,7) 
!		cur(jx,jz,jy,2)=x1z(jx,jz,jy,6)-x1r(jx,jz,jy,8)
		cur(jx,jz,jy,3)=(drby_dx-xy(jx,jz,jy,6))/xx(jx)
        endif

    2 continue    
      !$acc end parallel loop
        
    	call mpi_transfersm(cur(:,:,:,:),3)
	if(smoothc) then
		  do m=1,3
		  call smthxzy(cur(:,:,:,m),1)
		  enddo
	endif
!      call bndry_cur_ex(lbnd)
!      call bndry3_ex(cur,lbnd)
	! need to be replaced by bndry3_ex_cut_cell   ! haowei 1219
!      call mpi_transfer3(cur) 
!      call convtc
      if(lcd .eq. 2) then 
      call current_driven
      cur(:,:,:,:)=cur(:,:,:,:)+cud(:,:,:,:)
      endif
!      write(*,*)'cur'
      return
      endsubroutine current_openacc    	
!


!hw***************************************************************************************
!calculate the xy, xy2 for all points, and x1z, x2z, x1r, x2r for IR points
	subroutine convt_openacc
      use declare
	integer itag, i
      real*8, dimension(my) :: wwy 
      real*8, dimension(mx,mz,my) :: x1r_tmp,xr2_tmp,x1z_tmp,xz2_tmp,x1_tmp
      real*8 d1f2, d1f2m, d1f2p, d1fc, d2fc, d1fp, d1fm, d1fbp, d1fbm
      real*8 d2f2, d2fbp, d2fbm, timestart1, timeend1
      include 'mpif.h'
! R R: calculated as subroutine right
! R IR2: the IR2 part need to be calculated with the boundary point by d1fc, order 4th
! R IR1: the IR1 part need to be calculated with the boundary point by d1fbm, d1fbp, order 3th
! IR2 IR2: both two IR2 parts need to be calculated with the boundary point by d1fc, order 4th
! IR1 IR1: both tow IR1 parts need to be calculated with the boundary point by d1fbm, d1fbp, order 3th
! IR1 IR2: the IR1 part need to be calculated with the boundary point by d1fbm, d1fbp, order 3th
! 	     the IR2 parts need to be calculated with the boundary point by d1fc, order 4th
! R or IR1 or IR2 with D: interpolation in the dropped direction
! D D: do not calculate

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
	 !!not spectral, need to add spectral if necessary, haowei

	!$acc parallel loop &
	!$acc local(jx,jz,jy,m) &
	!$acc copyin(ay1,by1,cy1,dy1,ay2,by2,cy2,dy2,iy_first,iy_last,iz_first,iz_last,ix_first,ix_last) &
	!$acc annotate(entire(ay1,by1,cy1,dy1,ay2,by2,cy2,dy2))
	do 15 jy=iy_first+2,iy_last-2
      do 15 m=1,8
	!$acc data copyin(x1(*,*,jy-2:jy+2,m)) &
	!$acc copyout(xy(*,*,jy,m),xy2(*,*,jy,m)) &
	!$acc present(ay1,by1,cy1,dy1,ay2,by2,cy2,dy2,iy_first,iy_last,iz_first,iz_last,ix_first,ix_last,jx,jz,jy)
      do 16 jz=iz_first,iz_last
      do 16 jx=ix_first,ix_last
      xy(jx,jz,jy,m) =d1fc(x1(jx,jz,jy-2,m),x1(jx,jz,jy-1,m),x1(jx,jz,jy,m),x1(jx,jz,jy+1,m),x1(jx,jz,jy+2,m),ay1(jy),by1(jy),cy1(jy),dy1(jy))      
      xy2(jx,jz,jy,m)=d2fc(x1(jx,jz,jy-2,m),x1(jx,jz,jy-1,m),x1(jx,jz,jy,m),x1(jx,jz,jy+1,m),x1(jx,jz,jy+2,m),ay2(jy),by2(jy),cy2(jy),dy2(jy))
   16 continue
      !$acc end data

! for the xy_8bndx and xy2_8bndx and xy_8bndz and xy2_8bndz
	!$acc data copyin(x1_8bndx(*,jy-2:jy+2,m),nbndx) &
	!$acc copyout(xy_8bndx(*,jy,m),xy2_8bndx(*,jy,m)) &
	!$acc annotate(dimension(x1_8bndx(nbndx,my,8),xy_8bndx(nbndx,my,8),xy2_8bndx(nbndx,my,8))) &
	!$acc present(ay1,by1,cy1,dy1,ay2,by2,cy2,dy2,iy_first,iy_last,iz_first,iz_last,ix_first,ix_last,jx,jz,jy)
	do 17 jx=1,nbndx
	xy_8bndx(jx,jy,m)=d1fc(x1_8bndx(jx,jy-2,m),x1_8bndx(jx,jy-1,m),x1_8bndx(jx,jy,m),x1_8bndx(jx,jy+1,m),x1_8bndx(jx,jy+2,m),&
		  ay1(jy),by1(jy),cy1(jy),dy1(jy))
	xy2_8bndx(jx,jy,m)=d2fc(x1_8bndx(jx,jy-2,m),x1_8bndx(jx,jy-1,m),x1_8bndx(jx,jy,m),x1_8bndx(jx,jy+1,m),x1_8bndx(jx,jy+2,m),&
		  ay2(jy),by2(jy),cy2(jy),dy2(jy))
   17 continue
      !$acc end data

	!$acc data copyin(x1_8bndz(*,jy-2:jy+2,m),nbndz) &
	!$acc copyout(xy_8bndz(*,jy,m),xy2_8bndz(*,jy,m)) &
	!$acc annotate(dimension(x1_8bndz(nbndz,my,8),xy_8bndz(nbndz,my,8),xy2_8bndz(nbndz,my,8))) &
	!$acc present(ay1,by1,cy1,dy1,ay2,by2,cy2,dy2,iy_first,iy_last,iz_first,iz_last,ix_first,ix_last,jx,jz,jy)
	do 18 jx=1,nbndz
	xy_8bndz(jx,jy,m)=d1fc(x1_8bndz(jx,jy-2,m),x1_8bndz(jx,jy-1,m),x1_8bndz(jx,jy,m),x1_8bndz(jx,jy+1,m),x1_8bndz(jx,jy+2,m),&
		  ay1(jy),by1(jy),cy1(jy),dy1(jy))
	xy2_8bndz(jx,jy,m)=d2fc(x1_8bndz(jx,jy-2,m),x1_8bndz(jx,jy-1,m),x1_8bndz(jx,jy,m),x1_8bndz(jx,jy+1,m),x1_8bndz(jx,jy+2,m),&
		  ay2(jy),by2(jy),cy2(jy),dy2(jy))
   18 continue
      !$acc end data
   15 continue
      !$acc end parallel loop

      !not spectral

	!$acc parallel loop &
	!$acc local(jx,jz,jy,m) &
	!$acc swapout(x1r(dimension order:1,2,4,3)) &
	!$acc swapin(x1(dimension order:1,2,4,3)) &
	!$acc copyin(ax1,bx1,cx1,dx1,iy_first,iy_last,iz_first,iz_last,ix_first,ix_last) &
	!$acc annotate(entire(ax1,bx1,cx1,dx1))
      do 1 jy=iy_first,iy_last
      do 1 m=1,8
      do 1 jz=iz_first+2,iz_last-2
      do 1 jx=ix_first+2,ix_last-2
      x1r(jx,jz,jy,m)=d1fc(x1(jx-2,jz,jy,m),x1(jx-1,jz,jy,m),x1(jx,jz,jy,m) &
          ,x1(jx+1,jz,jy,m),x1(jx+2,jz,jy,m),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
    1 continue
      !$acc end parallel loop

	!$acc parallel loop &
	!$acc local(jx,jz,jy,m) &
	!$acc swapout(xr2(dimension order:1,2,4,3)) &
	!$acc swapin(x1(dimension order:1,2,4,3)) &
	!$acc copyin(ax2,bx2,cx2,dx2,iy_first,iy_last,iz_first,iz_last,ix_first,ix_last) &
	!$acc annotate(entire(ax2,bx2,cx2,dx2))
      do 2 jy=iy_first,iy_last
      do 2 m=1,8
      do 2 jz=iz_first+2,iz_last-2
      do 2 jx=ix_first+2,ix_last-2
      xr2(jx,jz,jy,m)=d2fc(x1(jx-2,jz,jy,m),x1(jx-1,jz,jy,m),x1(jx,jz,jy,m) &
          ,x1(jx+1,jz,jy,m),x1(jx+2,jz,jy,m),ax2(jx),bx2(jx),cx2(jx),dx2(jx))
    2 continue
      !$acc end parallel loop

	!$acc parallel loop &
	!$acc local(jx,jz,jy,m) &
	!$acc swapout(x1z(dimension order:1,2,4,3)) &
	!$acc swapin(x1(dimension order:1,2,4,3)) &
	!$acc copyin(az1,bz1,cz1,dz1,iy_first,iy_last,iz_first,iz_last,ix_first,ix_last) &
	!$acc annotate(entire(az1,bz1,cz1,dz1))
      do 3 jy=iy_first,iy_last
      do 3 m=1,8
      do 3 jz=iz_first+2,iz_last-2
      do 3 jx=ix_first+2,ix_last-2
      x1z(jx,jz,jy,m)=d1fc(x1(jx,jz-2,jy,m),x1(jx,jz-1,jy,m),x1(jx,jz,jy,m) &
          ,x1(jx,jz+1,jy,m),x1(jx,jz+2,jy,m),az1(jz),bz1(jz),cz1(jz),dz1(jz))
    3 continue
      !$acc end parallel loop
	
	!$acc parallel loop &
	!$acc local(jx,jz,jy,m) &
	!$acc swapout(xz2(dimension order:1,2,4,3)) &
	!$acc swapin(x1(dimension order:1,2,4,3)) &
	!$acc copyin(az2,bz2,cz2,dz2,iy_first,iy_last,iz_first,iz_last,ix_first,ix_last) &
	!$acc annotate(entire(az2,bz2,cz2,dz2))
      do 4 jy=iy_first,iy_last
      do 4 m=1,8
      do 4 jz=iz_first+2,iz_last-2
      do 4 jx=ix_first+2,ix_last-2
      xz2(jx,jz,jy,m)=d2fc(x1(jx,jz-2,jy,m),x1(jx,jz-1,jy,m),x1(jx,jz,jy,m) &
          ,x1(jx,jz+1,jy,m),x1(jx,jz+2,jy,m),az2(jz),bz2(jz),cz2(jz),dz2(jz))
    4 continue
      !$acc end parallel loop

!      do m=1,8
!      do jy=iy_first,iy_last
!      do jz=iz_first+2,iz_last-2
!      do jx=ix_first+2,ix_last-2
!
!      x1r(jx,jz,jy,m)=d1fc(x1(jx-2,jz,jy,m),x1(jx-1,jz,jy,m),x1(jx,jz,jy,m) &
!          ,x1(jx+1,jz,jy,m),x1(jx+2,jz,jy,m),ax1(jx),bx1(jx),cx1(jx),dx1(jx))
!      xr2(jx,jz,jy,m)=d2fc(x1(jx-2,jz,jy,m),x1(jx-1,jz,jy,m),x1(jx,jz,jy,m) &
!          ,x1(jx+1,jz,jy,m),x1(jx+2,jz,jy,m),ax2(jx),bx2(jx),cx2(jx),dx2(jx))
!      x1z(jx,jz,jy,m)=d1fc(x1(jx,jz-2,jy,m),x1(jx,jz-1,jy,m),x1(jx,jz,jy,m) &
!          ,x1(jx,jz+1,jy,m),x1(jx,jz+2,jy,m),az1(jz),bz1(jz),cz1(jz),dz1(jz))
!      xz2(jx,jz,jy,m)=d2fc(x1(jx,jz-2,jy,m),x1(jx,jz-1,jy,m),x1(jx,jz,jy,m) &
!          ,x1(jx,jz+1,jy,m),x1(jx,jz+2,jy,m),az2(jz),bz2(jz),cz2(jz),dz2(jz))
!  
!      enddo
!      enddo
!      enddo
!      enddo

! can it be vectorization ????????????????????????????????????????????????????????
! d1fc(x1(jx,jz,:,:)) ??
! cost about 25s for total 165s

	!$acc parallel loop &
	!$acc local(jx,jz,jy,m,itag) &
	!$acc copyin(iy_first,iy_last,iz_first_irpt2,iz_last_irpt2,ix_first_irpt2,ix_last_irpt2,gdtp_ep,ax1_irz,bx1_irz,cx1_irz,dx1_irz,ax2_irz,bx2_irz,cx2_irz,dx2_irz) &
	!$acc annotate(entire(gdtp_ep,ax1_irz,bx1_irz,cx1_irz,dx1_irz,ax2_irz,bx2_irz,cx2_irz,dx2_irz))
	do 5 jy=iy_first,iy_last
      do 51 m=1,8
	!$acc data copyin(x1(*,*,jy,m),x1_8bndz(*,jy,m)) &
	!$acc copy(x1r(*,*,jy,m),xr2(*,*,jy,m)) &
	!$acc annotate(dimension(x1_8bndz(nbndz,my,8))) &
	!$acc present(jx,jz,jy,m,itag,iy_first,iy_last,iz_first_irpt2,iz_last_irpt2,ix_first_irpt2,ix_last_irpt2,gdtp_ep,ax1_irz,bx1_irz,cx1_irz,dx1_irz,ax2_irz,bx2_irz,cx2_irz,dx2_irz)
	do 52 jz=iz_first_irpt2,iz_last_irpt2
	do 52 jx=ix_first_irpt2,ix_last_irpt2
	! for IR in x direction, x1r xr2 need to be recalculated
	if(gdtp_ep(jx,jz,4).eq.2) then 
	
		  itag=gdtp_ep(jx,jz,6)
		  if(gdtp_ep(jx,jz,5).eq.-2) then ! x1(jx-2,jz,jy,m)=bndz
			    x1r(jx,jz,jy,m)=d1fc(x1_8bndz(itag,jy,m),x1(jx-1,jz,jy,m),x1(jx,jz,jy,m)&
					,x1(jx+1,jz,jy,m),x1(jx+2,jz,jy,m),ax1_irz(jx,jz),bx1_irz(jx,jz)&
					,cx1_irz(jx,jz),dx1_irz(jx,jz))
			    xr2(jx,jz,jy,m)=d2fc(x1_8bndz(itag,jy,m),x1(jx-1,jz,jy,m),x1(jx,jz,jy,m)&
					,x1(jx+1,jz,jy,m),x1(jx+2,jz,jy,m),ax2_irz(jx,jz),bx2_irz(jx,jz)&
					,cx2_irz(jx,jz),dx2_irz(jx,jz))
		  else if(gdtp_ep(jx,jz,5).eq.2) then ! x1(jx+2,jz,jy,m)=bndz
			    x1r(jx,jz,jy,m)=d1fc(x1(jx-2,jz,jy,m),x1(jx-1,jz,jy,m),x1(jx,jz,jy,m)&
					,x1(jx+1,jz,jy,m),x1_8bndz(itag,jy,m),ax1_irz(jx,jz),bx1_irz(jx,jz)&
					,cx1_irz(jx,jz),dx1_irz(jx,jz))
			    xr2(jx,jz,jy,m)=d2fc(x1(jx-2,jz,jy,m),x1(jx-1,jz,jy,m),x1(jx,jz,jy,m)&
					,x1(jx+1,jz,jy,m),x1_8bndz(itag,jy,m),ax2_irz(jx,jz),bx2_irz(jx,jz)&
					,cx2_irz(jx,jz),dx2_irz(jx,jz))
		  endif
	endif

   52 continue
   	!$acc end data
   51 continue
   5  continue
   	!$acc end parallel loop


	!$acc parallel loop &
	!$acc local(jx,jz,jy,m,itag) &
	!$acc copyin(iy_first,iy_last,iz_first_irpt2,iz_last_irpt2,ix_first_irpt2,ix_last_irpt2,gdtp_ep,axbp_irz,bxbp_irz,cxbp_irz,a2xbp_irz,b2xbp_irz,c2xbp_irz,axbm_irz,bxbm_irz,cxbm_irz,a2xbm_irz,b2xbm_irz,c2xbm_irz) &
	!$acc annotate(entire(gdtp_ep,axbp_irz,bxbp_irz,cxbp_irz,a2xbp_irz,b2xbp_irz,c2xbp_irz,axbm_irz,bxbm_irz,cxbm_irz,a2xbm_irz,b2xbm_irz,c2xbm_irz))
	do 6 jy=iy_first,iy_last
      do 61 m=1,8
	!$acc data copyin(x1(*,*,jy,m),x1_8bndz(*,jy,m)) &
	!$acc copy(x1r(*,*,jy,m),xr2(*,*,jy,m)) &
	!$acc annotate(dimension(x1_8bndz(nbndz,my,8))) &
	!$acc present(jx,jz,jy,m,itag,iy_first,iy_last,iz_first_irpt2,iz_last_irpt2,ix_first_irpt2,ix_last_irpt2,gdtp_ep,axbp_irz,bxbp_irz,cxbp_irz,a2xbp_irz,b2xbp_irz,c2xbp_irz,axbm_irz,bxbm_irz,cxbm_irz,a2xbm_irz,b2xbm_irz,c2xbm_irz)
	do 62 jz=iz_first_irpt2,iz_last_irpt2
	do 62 jx=ix_first_irpt2,ix_last_irpt2
	! for IR in x direction, x1r xr2 need to be recalculated
	if(gdtp_ep(jx,jz,4).eq.2) then 
	
		  itag=gdtp_ep(jx,jz,6)
              if(gdtp_ep(jx,jz,5).eq.-1) then !x1(jx-1,jz,jy,m)=bndz d1fbp d2fbp
			    x1r(jx,jz,jy,m)=d1fbp(x1(jx+2,jz,jy,m),x1(jx+1,jz,jy,m),x1(jx,jz,jy,m)&
					,x1_8bndz(itag,jy,m),axbp_irz(jx,jz),bxbp_irz(jx,jz),cxbp_irz(jx,jz))
			    xr2(jx,jz,jy,m)=d2fbp(x1_8bndz(itag,jy,m),x1(jx,jz,jy,m),x1(jx+1,jz,jy,m)&
					,x1(jx+2,jz,jy,m),a2xbp_irz(jx,jz),b2xbp_irz(jx,jz),c2xbp_irz(jx,jz))
		  else if(gdtp_ep(jx,jz,5).eq.1) then !x1(jx+1,jz,jy,m)=bndz d1fbm d2fbm
			    x1r(jx,jz,jy,m)=d1fbm(x1(jx-2,jz,jy,m),x1(jx-1,jz,jy,m),x1(jx,jz,jy,m)&
					,x1_8bndz(itag,jy,m),axbm_irz(jx,jz),bxbm_irz(jx,jz),cxbm_irz(jx,jz))
			    xr2(jx,jz,jy,m)=d2fbm(x1(jx-2,jz,jy,m),x1(jx-1,jz,jy,m),x1(jx,jz,jy,m)&
					,x1_8bndz(itag,jy,m),a2xbm_irz(jx,jz),b2xbm_irz(jx,jz),c2xbm_irz(jx,jz))
		  endif
	endif

   62 continue
   	!$acc end data
   61 continue
   6  continue
   	!$acc end parallel loop


	!$acc parallel loop &
	!$acc local(jx,jz,jy,m,itag) &
	!$acc copyin(iy_first,iy_last,iz_first_irpt2,iz_last_irpt2,ix_first_irpt2,ix_last_irpt2,gdtp_ep,az1_irx,bz1_irx,cz1_irx,dz1_irx,az2_irx,bz2_irx,cz2_irx,dz2_irx) &
	!$acc annotate(entire(gdtp_ep,az1_irx,bz1_irx,cz1_irx,dz1_irx,az2_irx,bz2_irx,cz2_irx,dz2_irx))
	do 7 jy=iy_first,iy_last
      do 71 m=1,8
	!$acc data copyin(x1(*,*,jy,m),x1_8bndx(*,jy,m)) &
	!$acc copy(x1z(*,*,jy,m),xz2(*,*,jy,m)) &
	!$acc annotate(dimension(x1_8bndx(nbndx,my,8))) &
	!$acc present(jx,jz,jy,m,itag,iy_first,iy_last,iz_first_irpt2,iz_last_irpt2,ix_first_irpt2,ix_last_irpt2,gdtp_ep,az1_irx,bz1_irx,cz1_irx,dz1_irx,az2_irx,bz2_irx,cz2_irx,dz2_irx)
	do 72 jz=iz_first_irpt2,iz_last_irpt2
	do 72 jx=ix_first_irpt2,ix_last_irpt2

	! for IR in z direction, x1z xz2 need to be recalculated
	if(gdtp_ep(jx,jz,1).eq.2) then 
	
		  itag=gdtp_ep(jx,jz,3)
		  if(gdtp_ep(jx,jz,2).eq.-2) then ! x1(jx,jz-2,jy,m)=bndx
			    x1z(jx,jz,jy,m)=d1fc(x1_8bndx(itag,jy,m),x1(jx,jz-1,jy,m),x1(jx,jz,jy,m)&
					,x1(jx,jz+1,jy,m),x1(jx,jz+2,jy,m),az1_irx(jx,jz),bz1_irx(jx,jz)&
					,cz1_irx(jx,jz),dz1_irx(jx,jz))
			    xz2(jx,jz,jy,m)=d2fc(x1_8bndx(itag,jy,m),x1(jx,jz-1,jy,m),x1(jx,jz,jy,m)&
					,x1(jx,jz+1,jy,m),x1(jx,jz+2,jy,m),az2_irx(jx,jz),bz2_irx(jx,jz)&
					,cz2_irx(jx,jz),dz2_irx(jx,jz))
		  else if(gdtp_ep(jx,jz,2).eq.2) then ! x1(jx,jz+2,jy,m)=bndx
			    x1z(jx,jz,jy,m)=d1fc(x1(jx,jz-2,jy,m),x1(jx,jz-1,jy,m),x1(jx,jz,jy,m)&
					,x1(jx,jz+1,jy,m),x1_8bndx(itag,jy,m),az1_irx(jx,jz),bz1_irx(jx,jz)&
					,cz1_irx(jx,jz),dz1_irx(jx,jz))
			    xz2(jx,jz,jy,m)=d2fc(x1(jx,jz-2,jy,m),x1(jx,jz-1,jy,m),x1(jx,jz,jy,m)&
					,x1(jx,jz+1,jy,m),x1_8bndx(itag,jy,m),az2_irx(jx,jz),bz2_irx(jx,jz)&
					,cz2_irx(jx,jz),dz2_irx(jx,jz))
		  endif
	endif

   72 continue
   	!$acc end data
   71 continue
   7  continue
   	!$acc end parallel loop




	!$acc parallel loop &
	!$acc local(jx,jz,jy,m,itag) &
	!$acc copyin(iy_first,iy_last,iz_first_irpt2,iz_last_irpt2,ix_first_irpt2,ix_last_irpt2,gdtp_ep,azbp_irx,bzbp_irx,czbp_irx,a2zbp_irx,b2zbp_irx,c2zbp_irx,azbm_irx,bzbm_irx,czbm_irx,a2zbm_irx,b2zbm_irx,c2zbm_irx) &
	!$acc annotate(entire(gdtp_ep,azbp_irx,bzbp_irx,czbp_irx,a2zbp_irx,b2zbp_irx,c2zbp_irx,azbm_irx,bzbm_irx,czbm_irx,a2zbm_irx,b2zbm_irx,c2zbm_irx))
	do 8 jy=iy_first,iy_last
      do 81 m=1,8
	!$acc data copyin(x1(*,*,jy,m),x1_8bndx(*,jy,m)) &
	!$acc copy(x1z(*,*,jy,m),xz2(*,*,jy,m)) &
	!$acc annotate(dimension(x1_8bndx(nbndx,my,8))) &
	!$acc present(jx,jz,jy,m,itag,iy_first,iy_last,iz_first_irpt2,iz_last_irpt2,ix_first_irpt2,ix_last_irpt2,gdtp_ep,azbp_irx,bzbp_irx,czbp_irx,a2zbp_irx,b2zbp_irx,c2zbp_irx,azbm_irx,bzbm_irx,czbm_irx,a2zbm_irx,b2zbm_irx,c2zbm_irx)
	do 82 jz=iz_first_irpt2,iz_last_irpt2
	do 82 jx=ix_first_irpt2,ix_last_irpt2

	! for IR in z direction, x1z xz2 need to be recalculated
	if(gdtp_ep(jx,jz,1).eq.2) then 
	
		  itag=gdtp_ep(jx,jz,3)
	        if(gdtp_ep(jx,jz,2).eq.-1) then ! x1(jx,jz-1,jy,m)=bndx d1fbp d2fbp
			    x1z(jx,jz,jy,m)=d1fbp(x1(jx,jz+2,jy,m),x1(jx,jz+1,jy,m),x1(jx,jz,jy,m)&
					,x1_8bndx(itag,jy,m),azbp_irx(jx,jz),bzbp_irx(jx,jz),czbp_irx(jx,jz))
			    xz2(jx,jz,jy,m)=d2fbp(x1_8bndx(itag,jy,m),x1(jx,jz,jy,m),x1(jx,jz+1,jy,m)&
					,x1(jx,jz+2,jy,m),a2zbp_irx(jx,jz),b2zbp_irx(jx,jz),c2zbp_irx(jx,jz))
		  else if(gdtp_ep(jx,jz,2).eq.1) then ! x1(jx,jz+1,jy,m)=bndx d1fbm d2fbm
			    x1z(jx,jz,jy,m)=d1fbm(x1(jx,jz-2,jy,m),x1(jx,jz-1,jy,m),x1(jx,jz,jy,m)&
					,x1_8bndx(itag,jy,m),azbm_irx(jx,jz),bzbm_irx(jx,jz),czbm_irx(jx,jz))
			    xz2(jx,jz,jy,m)=d2fbm(x1(jx,jz-2,jy,m),x1(jx,jz-1,jy,m),x1(jx,jz,jy,m)&
					,x1_8bndx(itag,jy,m),a2zbm_irx(jx,jz),b2zbm_irx(jx,jz),c2zbm_irx(jx,jz))
		  endif
	endif

   82 continue
   	!$acc end data
   81 continue
   8  continue
   	!$acc end parallel loop


   	do 30 m=1,8
	do 30 jy=iy_first,iy_last
	xr(:,:,jy,m)=xint_dx(:,:,m)+x1r(:,:,jy,m)  
	xz(:,:,jy,m)=xint_dz(:,:,m)+x1z(:,:,jy,m)
   30 continue

    return
    endsubroutine convt_openacc
!

!hw************************************************************************
!ws************************************************************************
!wzhang************************************************************
      subroutine efield_openacc
      use declare
	implicit none
	integer itag
      include 'mpif.h'

	!$acc parallel loop &
	!$acc local(jx,jz,jy) &
	!$acc swapout(ef(dimension order:1,2,4,3)) &
	!$acc copyin(x,x1,xint3,xint4,xint5,xint6,xint7,xint8,linear_mhd,iy_first,iy_last,iz_first,iz_last,ix_first,ix_last) &
	!$acc annotate(entire(xint3,xint4,xint5,xint6,xint7,xint8)) &
	!$acc annotate(slice(x1(*,*,jy,3),x1(*,*,jy,4),x1(*,*,jy,5),x1(*,*,jy,6),x1(*,*,jy,7),x1(*,*,jy,8),x(*,*,jy,6),x(*,*,jy,7),x(*,*,jy,8)))
      do 1 jy=iy_first+2,iy_last-2
!      do 1 jy=iy_first,iy_last
      do 1 jz=iz_first,iz_last
      do 1 jx=ix_first,ix_last

      if(linear_mhd) then

!          ef(jx,jz,jy,1)=-x1(jx,jz,jy,4)*xint(jx,jz,8)-xint(jx,jz,4)*x1(jx,jz,jy,8)+x1(jx,jz,jy,5)*xint(jx,jz,7)+xint(jx,jz,5)*x1(jx,jz,jy,7)
!          ef(jx,jz,jy,2)=-x1(jx,jz,jy,5)*xint(jx,jz,6)-xint(jx,jz,5)*x1(jx,jz,jy,6)+x1(jx,jz,jy,3)*xint(jx,jz,8)+xint(jx,jz,3)*x1(jx,jz,jy,8)
!          ef(jx,jz,jy,3)=-x1(jx,jz,jy,3)*xint(jx,jz,7)-xint(jx,jz,3)*x1(jx,jz,jy,7)+x1(jx,jz,jy,4)*xint(jx,jz,6)+xint(jx,jz,4)*x1(jx,jz,jy,6)
          ef(jx,jz,jy,1)=-x1(jx,jz,jy,4)*xint8(jx,jz)-xint4(jx,jz)*x1(jx,jz,jy,8)+x1(jx,jz,jy,5)*xint7(jx,jz)+xint5(jx,jz)*x1(jx,jz,jy,7)
          ef(jx,jz,jy,2)=-x1(jx,jz,jy,5)*xint6(jx,jz)-xint5(jx,jz)*x1(jx,jz,jy,6)+x1(jx,jz,jy,3)*xint8(jx,jz)+xint3(jx,jz)*x1(jx,jz,jy,8)
          ef(jx,jz,jy,3)=-x1(jx,jz,jy,3)*xint7(jx,jz)-xint3(jx,jz)*x1(jx,jz,jy,7)+x1(jx,jz,jy,4)*xint6(jx,jz)+xint4(jx,jz)*x1(jx,jz,jy,6)

	else

!      ef(jx,jz,jy,1)=-x1(jx,jz,jy,4)*x(jx,jz,jy,8)-xint(jx,jz,4)*x1(jx,jz,jy,8)+x1(jx,jz,jy,5)*x(jx,jz,jy,7)+xint(jx,jz,5)*x1(jx,jz,jy,7)
!      ef(jx,jz,jy,2)=-x1(jx,jz,jy,5)*x(jx,jz,jy,6)-xint(jx,jz,5)*x1(jx,jz,jy,6)+x1(jx,jz,jy,3)*x(jx,jz,jy,8)+xint(jx,jz,3)*x1(jx,jz,jy,8)
!      ef(jx,jz,jy,3)=-x1(jx,jz,jy,3)*x(jx,jz,jy,7)-xint(jx,jz,3)*x1(jx,jz,jy,7)+x1(jx,jz,jy,4)*x(jx,jz,jy,6)+xint(jx,jz,4)*x1(jx,jz,jy,6)
      ef(jx,jz,jy,1)=-x1(jx,jz,jy,4)*x(jx,jz,jy,8)-xint4(jx,jz)*x1(jx,jz,jy,8)+x1(jx,jz,jy,5)*x(jx,jz,jy,7)+xint5(jx,jz)*x1(jx,jz,jy,7)
      ef(jx,jz,jy,2)=-x1(jx,jz,jy,5)*x(jx,jz,jy,6)-xint5(jx,jz)*x1(jx,jz,jy,6)+x1(jx,jz,jy,3)*x(jx,jz,jy,8)+xint3(jx,jz)*x1(jx,jz,jy,8)
      ef(jx,jz,jy,3)=-x1(jx,jz,jy,3)*x(jx,jz,jy,7)-xint3(jx,jz)*x1(jx,jz,jy,7)+x1(jx,jz,jy,4)*x(jx,jz,jy,6)+xint4(jx,jz)*x1(jx,jz,jy,6)

	endif

    1 continue
      !$acc end parallel loop
      
       if(hall) then
	!$acc parallel loop &
	!$acc local(jx,jz,jy) &
	!$acc swap(ef(dimension order:1,2,4,3)) &
	!$acc swapin(cur(dimension order:1,2,4,3)) &
	!$acc copyin(x,x1,fdi,cint,xy,x1r,x1z,xx,linear_mhd,iy_first,iy_last,iz_first,iz_last,ix_first,ix_last) &
	!$acc annotate(entire(fdi,cint,xx)) &
	!$acc annotate(slice(x1(*,*,jy,6),x1(*,*,jy,7),x1(*,*,jy,8),x(*,*,jy,1),x(*,*,jy,6),x(*,*,jy,7),x(*,*,jy,8),xy(*,*,jy,2),x1z(*,*,jy,2),x1r(*,*,jy,2)))
      do 12 jy=iy_first+2,iy_last-2
      do 12 jz=iz_first+2,iz_last-2
      do 12 jx=ix_first+2,ix_last-2

      if(linear_mhd) then

      ef(jx,jz,jy,1)=ef(jx,jz,jy,1)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,2)*x1(jx,jz,jy,8)-cur(jx,jz,jy,3)*x1(jx,jz,jy,7)+cint(jx,jz,2)*x1(jx,jz,jy,8)-cint(jx,jz,3)*x1(jx,jz,jy,7)-x1r(jx,jz,jy,2))
      ef(jx,jz,jy,2)=ef(jx,jz,jy,2)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,3)*x1(jx,jz,jy,6)-cur(jx,jz,jy,1)*x1(jx,jz,jy,8)+cint(jx,jz,3)*x1(jx,jz,jy,6)-cint(jx,jz,1)*x1(jx,jz,jy,8)-xy(jx,jz,jy,2)/xx(jx))    
      ef(jx,jz,jy,3)=ef(jx,jz,jy,3)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,1)*x1(jx,jz,jy,7)-cur(jx,jz,jy,2)*x1(jx,jz,jy,6)+cint(jx,jz,1)*x1(jx,jz,jy,7)-cint(jx,jz,2)*x1(jx,jz,jy,6)-x1z(jx,jz,jy,2))

	else

      ef(jx,jz,jy,1)=ef(jx,jz,jy,1)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,2)*x(jx,jz,jy,8)-cur(jx,jz,jy,3)*x(jx,jz,jy,7)+cint(jx,jz,2)*x1(jx,jz,jy,8)-cint(jx,jz,3)*x1(jx,jz,jy,7)-x1r(jx,jz,jy,2))
      ef(jx,jz,jy,2)=ef(jx,jz,jy,2)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,3)*x(jx,jz,jy,6)-cur(jx,jz,jy,1)*x(jx,jz,jy,8)+cint(jx,jz,3)*x1(jx,jz,jy,6)-cint(jx,jz,1)*x1(jx,jz,jy,8)-xy(jx,jz,jy,2)/xx(jx))    
      ef(jx,jz,jy,3)=ef(jx,jz,jy,3)+fdi(jx,jz)/x(jx,jz,jy,1)*(cur(jx,jz,jy,1)*x(jx,jz,jy,7)-cur(jx,jz,jy,2)*x(jx,jz,jy,6)+cint(jx,jz,1)*x1(jx,jz,jy,7)-cint(jx,jz,2)*x1(jx,jz,jy,6)-x1z(jx,jz,jy,2))
      
	endif
   12 continue
      !$acc end parallel loop
	endif
       

      if(resisitive .and. etaj_in_e) then
	!$acc parallel loop &
	!$acc local(jy,m) &
	!$acc swap(ef(dimension order:1,2,4,3)) &
	!$acc swapin(cur(dimension order:1,2,4,3)) &
	!$acc copyin(eta,iy_first,iy_last)
      do jy=iy_first+2,iy_last-2
!      do jy=iy_first,iy_last
      do m=1,3
      ef(:,:,jy,m)=ef(:,:,jy,m)+eta(:,:,jy)*cur(:,:,jy,m)
      enddo
	enddo
      !$acc end parallel loop
      endif
	

      if(bootstrap) then
      call current_boot(lbs)      
      do m=1,3
      ef(:,:,:,m)=ef(:,:,:,m)-eta(:,:,:)*cub(:,:,:,m)
      enddo
      endif


      if(nstep.lt.nper) then
      do 11 m=1,3
      do 11 jy=iy_first,iy_last
      do 11 jz=iz_first,iz_last
      do 11 jx=ix_first,ix_last
      ef(jx,jz,jy,m)=ef(jx,jz,jy,m)+eta1(jx,jz,jy)*(cint(jx,jz,m)+cur(jx,jz,jy,m))
  11  continue
      endif 

       !***************revised**************************************
	call mpi_transfersm(ef(:,:,:,:),3)
    ! haowei need to obatin the efield for bnd grids, i.e., ef_3bndz, ef_3bndx

	ef_3bndx(:,:,:)=0.d0 ! fixed the ef_3bndx
	ef_3bndz(:,:,:)=0.d0 ! fixed the ef_3bndz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Boundary for efield
	if(hall) then ! Efield is complete when hall=.false.
	
!	do 3 jz=iz_first,iz_last
!	do 3 jx=ix_first,ix_last
	do 3 jz=iz_first_irpt,iz_last_irpt
	do 3 jx=ix_first_irpt,ix_last_irpt
	itag=gdtp_ep(jx,jz,6)
	if((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz,4).eq.5)) ef(jx,jz,:,:)=0.d0
	
	if((gdtp_ep(jx,jz,5).eq.1).or.(gdtp_ep(jx,jz,5).eq.-1)) then ! jx+-1 is the boundary point, use the jx, jx-1, jx-2 to calculate the bndz point
	
	ef_3bndz(itag,:,:)=0.d0
	
	else if((gdtp_ep(jx,jz,4).eq.5).and.(gdtp_ep(jx,jz,1).ne.5).and.(gdtp_ep(jx-1,jz,5).eq.1).and.(jx.gt.ix_first+1)) then ! jx is the inside dropped point
	
	do 22 m=1,3
	do 22 jy=1,my
	call interp1d2l(ef(jx-2,jz,jy,m),ef(jx-1,jz,jy,m),ef_3bndz(itag,jy,m), &
		xx(jx-2),xx(jx-1),bndz_grd(itag,1),xx(jx),ef(jx,jz,jy,m))
   22 continue

	else if((gdtp_ep(jx,jz,4).eq.5).and.(gdtp_ep(jx,jz,1).ne.5).and.(gdtp_ep(jx+1,jz,5).eq.-1).and.(jx.lt.ix_last-1)) then ! jx is the inside dropped point

	do 24 m=1,3
	do 24 jy=1,my
	call interp1d2l(ef(jx+2,jz,jy,m),ef(jx+1,jz,jy,m),ef_3bndz(itag,jy,m), &
		xx(jx+2),xx(jx+1),bndz_grd(itag,1),xx(jx),ef(jx,jz,jy,m))
   24 continue
	endif

!
	itag=gdtp_ep(jx,jz,3)
	if((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz,4).eq.5)) ef(jx,jz,:,:)=0.d0
	
	if((gdtp_ep(jx,jz,2).eq.1).or.(gdtp_ep(jx,jz,2).eq.-1)) then ! jz+1 is the boundary point, use the jz, jz-1, jz-2 to calculate the bndx point
	ef_3bndx(itag,:,:)=0.d0
	
	else if((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz,4).ne.5).and.(gdtp_ep(jx,jz-1,2).eq.1).and.(jz.gt.iz_first+1)) then ! jz is the inside dropped point
	do 26 m=1,3
	do 26 jy=1,my
	call interp1d2l(ef(jx,jz-2,jy,m),ef(jx,jz-1,jy,m),ef_3bndx(itag,jy,m), &
		zz(jz-2),zz(jz-1),bndx_grd(itag,2),zz(jz),ef(jx,jz,jy,m))
   26 continue

	else if((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz,4).ne.5).and.(gdtp_ep(jx,jz+1,2).eq.-1).and.(jz.lt.iz_last-1)) then ! jz is the inside dropped point
	do 28 m=1,3
	do 28 jy=1,my
	call interp1d2l(ef(jx,jz+2,jy,m),ef(jx,jz+1,jy,m),ef_3bndx(itag,jy,m), &
		zz(jz+2),zz(jz+1),bndx_grd(itag,2),zz(jz),ef(jx,jz,jy,m))
   28 continue
	endif

    3 continue

    	endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Boundary for efield

	if(nstep.eq.0) ef(:,:,:,:)=0.d0
	!$acc parallel loop tile(jy:2) &
	!$acc local(jx,jz,jy,m) &
	!$acc swap(ef(dimension order:1,2,4,3)) &
	!$acc copyin(hypb_ratio,iy_first,iy_last,iz_first,iz_last,ix_first,ix_last) &
	!$acc annotate(entire(hypb_ratio))
      do 10 jy=iy_first,iy_last
	do 10 m=1,3
      do 10 jz=iz_first,iz_last
      do 10 jx=ix_first,ix_last
!	if(psi(jx,jz).lt.psia1) then
!	if((hypb_ratio(jx,jz).ge.0.d0).and.(hypb_ratio(jx,jz).le.1.d0)) then
	ef(jx,jz,jy,m)=ef(jx,jz,jy,m)*hypb_ratio(jx,jz)
!	endif
!	endif
   10 continue
      !$acc end parallel loop
!       call mpi_transfersm(ef(:,:,:,:),3)
       !**********************************************************
      if(smoothef) call smthef_dis_v2(3)

!	call smth_irpt_with_difc(4,2,1,3,0.9d0)
!	call smth_irpt_with_difc_v2(4,2,1,3,0.9d0)

      return
      endsubroutine efield_openacc
!hw*****************************************************************    


!hw******************************************************************************
      subroutine stepon_openacc
!
!     this routine time-advances x's bz fourth order in time and second
!     order in space runge-kotta differential scheme.
!     note: x is alwazs the up-to-date value while xm being the
!           intermediate value, and xdif is increment
!
!
      use declare
	implicit none
      include 'mpif.h'

!      dts=dt/ncycl_atfs

      tt=time
      tt1=time+dt/6.
      tt2=time
      irk=1
      call right_openacc
       xfold(:,:,:,:)=x(:,:,:,:)
       xm(:,:,:,:)=xfold(:,:,:,:)+xdif(:,:,:,:)*dt/6.
       x(:,:,:,:)=xfold(:,:,:,:)+xdif(:,:,:,:)*dt/2.
!
      tt=time+dt/2.
      tt1=time+dt/2.
      tt2=time+dt/6.
      irk=2
        call right_openacc
        xm(:,:,:,:)=xm(:,:,:,:)+xdif(:,:,:,:)*dt/3.
        x(:,:,:,:)=xfold(:,:,:,:)+xdif(:,:,:,:)*dt/2.
!
      tt1=time+5.*dt/6.
      tt2=time+dt/2.
      irk=3
        call right_openacc
        xm(:,:,:,:)=xm(:,:,:,:)+xdif(:,:,:,:)*dt/3.
        x(:,:,:,:)=xfold(:,:,:,:)+xdif(:,:,:,:)*dt
!
      time=time+dt
      tt1=time+dt
      tt2=time+5.*dt/6.
      irk=4
        call right_openacc
        x(:,:,:,:)=xm(:,:,:,:)+xdif(:,:,:,:)*dt/6.

      caf=0.75d0*(0.5+0.5*dtanh((time-40)/5.))
!      call bndry_x_ex(lbnd)
!      call bndry8_x_ex(lbnd)
	call bndry8_cut_cell_v2_fixed
  !    if(conductpll) call pllconduct(lpll)
!      if(smoothpll) call smthp_traceline_5p(1)
!      if(eta_from_t) call calculate_eta
      
      return
      endsubroutine stepon_openacc
	

