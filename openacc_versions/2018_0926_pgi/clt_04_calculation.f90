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
	!$acc parallel present(x,xx,yy,zz,gdtp_ep)
	!$acc loop reduction(min:dt1) independent collapse(3)
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
	!$acc end parallel
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
	use declare_for_openacc
	use omp_lib
!	integer itag
!	real*8 drvx_dx,drey_dx,dez_dx,dex_dy,dez_dy,dex_dz,dey_dz
!	real*8, dimension(mx,mz,my) :: rvx,rey, drvx_dx_tmp, dex_dy_tmp, dez_dy_tmp
!	real*8, dimension(nbndx,my) :: rvx_bndx, rey_bndx
!	real*8, dimension(nbndz,my) :: rvx_bndz, rey_bndz
	real*8 d1fc, d1f2, d1fm, d1fbp ,d1fbm, d1xf2
	include 'mpif.h'
	!omp threadprivate(gdtp_ep,ax1_irz,bx1_irz,cx1_irz,dx1_irz,axbp_irz,bxbp_irz,cxbp_irz,axbm_irz,bxbm_irz,cxbm_irz,xx,ax1,bx1,cx1,dx1)
	!omp threadprivate(xint1_dx,xint1_dz,xint3_dx,xint5_dz,xint3,xint5,iz_first,iz_last,ix_first,ix_last)
	!omp threadprivate(xint2_dx,xint2_dz,gamma)
	!omp threadprivate(cint1,cint2,cint3)
	!omp threadprivate(xint3_dz,xint4,xint4_dx,xint4_dz,xint5_dx)
	!omp threadprivate(az1,bz1,cz1,dz1)
	!omp threadprivate(az1_irx,bz1_irx,cz1_irx,dz1_irx,azbp_irx,bzbp_irx,czbp_irx,azbm_irx,bzbm_irx,czbm_irx)
	!omp threadprivate(viscous,implicitv,pmu,pmux,pmuz,fmu,fmux,fmuz)
	
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
	
	!$acc kernels
	if(rho_from_p) then
	!$acc loop independent collapse(3)
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
	!$acc end kernels

	call convt_openacc
      call current_openacc
      call efield_openacc ! default if(ef_mode) 
 	

	!$omp parallel default(shared) &
	!$ private(jx,jz,jy)
	!$omp do

      !$acc kernels
	!$acc loop independent collapse(1)
      do jy=iy_first+2,iy_last-2
	!$acc loop independent collapse(2)
      do jz=iz_first,iz_last      
      do jx=ix_first,ix_last
      rvx(jx,jz,jy)=xx(jx)*x(jx,jz,jy,3)
	rey(jx,jz,jy)=xx(jx)*ef(jx,jz,jy,2)
      enddo
      enddo
	!$acc loop independent collapse(1)
	do itag=1,nbndz
	rvx_bndz(itag,jy)=bndz_grd(itag,1)*x_8bndz(itag,jy,3)
	rey_bndz(itag,jy)=bndz_grd(itag,1)*ef_3bndz(itag,jy,2)	
	enddo
	!calculate the rvx_bndx and rey_bndx
	!$acc loop independent collapse(1)
	do itag=1,nbndx
	rvx_bndx(itag,jy)=bndx_grd(itag,1)*x_8bndx(itag,jy,3)
	rey_bndx(itag,jy)=bndx_grd(itag,1)*ef_3bndx(itag,jy,2)
	enddo
      enddo
	!$omp end do
    	!$omp end parallel


	!$acc loop independent collapse(3)
      do jy=iy_first+2,iy_last-2
      do jz=iz_first+2,iz_last-2      
      do jx=ix_first+2,ix_last-2
	dex_dy_tmp(jx,jz,jy) =d1fc(ef(jx,jz,jy-2,1),ef(jx,jz,jy-1,1),ef(jx,jz,jy,1) &
		,ef(jx,jz,jy+1,1),ef(jx,jz,jy+2,1),ay1(jy),by1(jy),cy1(jy),dy1(jy)) 
	dez_dy_tmp(jx,jz,jy) =d1fc(ef(jx,jz,jy-2,3),ef(jx,jz,jy-1,3),ef(jx,jz,jy,3) &
		,ef(jx,jz,jy+1,3),ef(jx,jz,jy+2,3),ay1(jy),by1(jy),cy1(jy),dy1(jy))
	enddo
	enddo
	enddo

	!$omp parallel default(shared) &
	!$ private(jx,jz,jy,drey_dx,dez_dx,dex_dz,dey_dz)
	!$omp do
	!$acc loop independent collapse(3)
      do 1 jy=iy_first+2,iy_last-2
      do 1 jz=iz_first+2,iz_last-2      
      do 1 jx=ix_first+2,ix_last-2
! total(mx*mz*15+500+mx*5+mz*4+3)=6683

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

	
      xdif(jx,jz,jy,1)=-x1(jx,jz,jy,3)*(x1r(jx,jz,jy,1)+xint1_dx(jx,jz)) &
       -x(jx,jz,jy,1)*drvx_dx_tmp(jx,jz,jy)/xx(jx) &
!rplc       -x(jx,jz,jy,3)*x(jx,jz,jy,1)/xx(jx)-x(jx,jz,jy,1)*xr(jx,jz,jy,3) &
       -(xy(jx,jz,jy,1)*x(jx,jz,jy,4)+x(jx,jz,jy,1)*xy(jx,jz,jy,4))/xx(jx) &
       -x1(jx,jz,jy,5)*(x1z(jx,jz,jy,1)+xint1_dz(jx,jz))-x(jx,jz,jy,1)*x1z(jx,jz,jy,5) &
!ws:for poloidal flow
       -xint3(jx,jz)*x1r(jx,jz,jy,1)-x1(jx,jz,jy,1)*(xint3(jx,jz)/xx(jx)+xint3_dx(jx,jz)) &
       -xint5(jx,jz)*x1z(jx,jz,jy,1)-x1(jx,jz,jy,1)*xint5_dz(jx,jz)	
		
    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!														!
!------------------------------- Eq2, ddtP=R.H.S -----------------------------------!
!														!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    	
! second calculate the ddt(p)
! total(mx*mz*(14+1)+mx*5+mz*4+3)=7383
      xdif(jx,jz,jy,2)=-x1(jx,jz,jy,3)*(xint2_dx(jx,jz)+x1r(jx,jz,jy,2)) &
		-gamma*x(jx,jz,jy,2)*(drvx_dx_tmp(jx,jz,jy)/xx(jx)) &
!rplc       -gamma*x(jx,jz,jy,2)*(xr(jx,jz,jy,3)+x(jx,jz,jy,3)/xx(jx)) &
		-(xy(jx,jz,jy,2)*x(jx,jz,jy,4)+gamma*x(jx,jz,jy,2)*xy(jx,jz,jy,4))/xx(jx) &
		-x1(jx,jz,jy,5)*(xint2_dz(jx,jz)+x1z(jx,jz,jy,2))-gamma*x(jx,jz,jy,2)*x1z(jx,jz,jy,5) &    
!ws:for poloidal flow
		-xint3(jx,jz)*x1r(jx,jz,jy,2)-gamma*x1(jx,jz,jy,2)*(xint3(jx,jz)/xx(jx)+xint3_dx(jx,jz)) &
		-xint5(jx,jz)*x1z(jx,jz,jy,2)-gamma*x1(jx,jz,jy,2)*xint5_dz(jx,jz)	
	
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!														!
!------------------------------ Eq3, ddtVx=R.H.S, part1-----------------------------!
!														!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! calculate ddtvx, part1
! total(mx*mz*17+mx*5+mz*4+3)=7383
	
      xdif(jx,jz,jy,3)=-x(jx,jz,jy,3)*x1r(jx,jz,jy,3)-x(jx,jz,jy,5)*x1z(jx,jz,jy,3) &
		-x(jx,jz,jy,4)*xy(jx,jz,jy,3)/xx(jx)+x(jx,jz,jy,4)*x1(jx,jz,jy,4)/xx(jx) &
		+(cur(jx,jz,jy,2)*x(jx,jz,jy,8)+cint2(jx,jz)*x1(jx,jz,jy,8) &
		- cur(jx,jz,jy,3)*x(jx,jz,jy,7)-cint3(jx,jz)*x1(jx,jz,jy,7) &
		-x1r(jx,jz,jy,2))/x(jx,jz,jy,1) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!														!
!------------------------------ Eq4, ddtVy=R.H.S, part1-----------------------------!
!														!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
! calculate ddtvy, part1
! total(mx*mz*17+mx*5+mz*4+3)=7383    
	
	xdif(jx,jz,jy,4)=-x(jx,jz,jy,3)*x1r(jx,jz,jy,4)-x(jx,jz,jy,5)*x1z(jx,jz,jy,4) &
		-x(jx,jz,jy,4)*xy(jx,jz,jy,4)/xx(jx)-x(jx,jz,jy,4)*x1(jx,jz,jy,3)/xx(jx) &
		+(cur(jx,jz,jy,3)*x(jx,jz,jy,6)+cint3(jx,jz)*x1(jx,jz,jy,6) &
		- cur(jx,jz,jy,1)*x(jx,jz,jy,8)-cint1(jx,jz)*x1(jx,jz,jy,8) &   
		-xy(jx,jz,jy,2)/xx(jx))/x(jx,jz,jy,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!														!
!------------------------------ Eq5, ddtVz=R.H.S, part1-----------------------------!
!														!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate ddtvz, part1
! total(mx*mz*16+mx*5+mz*4+3)=6983
	
	xdif(jx,jz,jy,5)=-x(jx,jz,jy,3)*x1r(jx,jz,jy,5)-x(jx,jz,jy,5)*x1z(jx,jz,jy,5) &
		-x(jx,jz,jy,4)*xy(jx,jz,jy,5)/xx(jx) &     
		+(cur(jx,jz,jy,1)*x(jx,jz,jy,7)+cint1(jx,jz)*x1(jx,jz,jy,7) &
		- cur(jx,jz,jy,2)*x(jx,jz,jy,6)-cint2(jx,jz)*x1(jx,jz,jy,6) & 
		-x1z(jx,jz,jy,2))/x(jx,jz,jy,1)
		
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!														!
!------------------------------ Eq345, ddtV=R.H.S, part2~---------------------------!
!														!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 	
! calculate ddtvx, part2
! total(mx*mz*17+mx*1)=6820
	
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
	 
    
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!														!
!------------------------------ Eq6-8, ddtB=R.H.S, part1, for R points--------------!
!														!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate ddtB
! use ef_mode as default

! total(mx*mz*8+mx*5+mz*4)=3380

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
	
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!														!
!------------------------------ Eq6-8, ddtB=R.H.S, part2, for IR points-------------!
!														!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
! hard to use openacc
! add a subroutine to decide the jx and jz region for IR points
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
	endif


! hard to use openacc
! add a subroutine to decide the jx and jz region for IR points

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
	endif
    
    
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
!      if(psi(jx,jz).lt.psia1) then ! keep the viscous term only in closed flux
      xdif(jx,jz,jy,1)=xdif(jx,jz,jy,1)+pmu(jx,jz) &
       *(xr2(jx,jz,jy,1)+x1r(jx,jz,jy,1)/xx(jx)+xy2(jx,jz,jy,1)/xx(jx)**2+xz2(jx,jz,jy,1)) &
       +pmux(jx,jz)*x1r(jx,jz,jy,1)+pmuz(jx,jz)*x1z(jx,jz,jy,1) 

      xdif(jx,jz,jy,2)=xdif(jx,jz,jy,2)+kap(jx,jz) &
       *(xr2(jx,jz,jy,2)+x1r(jx,jz,jy,2)/xx(jx)+xy2(jx,jz,jy,2)/xx(jx)**2+xz2(jx,jz,jy,2)) &
       +kapx(jx,jz)*x1r(jx,jz,jy,2)+kapz(jx,jz)*x1z(jx,jz,jy,2)
	 
! total(mx*mz*19+mx+3)=7603
	
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
    
! LDM for openacc is not enough, so seperate it out    
	xdif(jx,jz,jy,3)=xdif(jx,jz,jy,3)+fmux(jx,jz)*x1r(jx,jz,jy,3)+fmuz(jx,jz)*x1z(jx,jz,jy,3) 
	xdif(jx,jz,jy,4)=xdif(jx,jz,jy,4)+fmux(jx,jz)*x1r(jx,jz,jy,4)+fmuz(jx,jz)*x1z(jx,jz,jy,4) 
	xdif(jx,jz,jy,5)=xdif(jx,jz,jy,5)+fmux(jx,jz)*x1r(jx,jz,jy,5)+fmuz(jx,jz)*x1z(jx,jz,jy,5) 
      endif   
      endif

!    1 continue
!	!$omp end do
!    	!$omp end parallel
	
! set the outside point into zero	
	! loop independent collapse(4)
!	do m=1,8
!	do jy=iy_first+2,iy_last-2
!	do jz=iz_first+2,iz_last-2
!	do jx=ix_first+2,ix_last-2
	if(max(gdtp_ep(jx,jz,1),gdtp_ep(jx,jz,4)).ge.4) xdif(jx,jz,jy,:)=0.d0
!	enddo
!	enddo
!	enddo
!	enddo

    1 continue
	!$omp end do
    	!$omp end parallel

      if(invaryrho) xdif(:,:,:,1)=0
      if(invaryp)   xdif(:,:,:,2)=0
	!$acc end kernels

	! update host(xdif(ix_last-3:ix_last-2,:,:,:),xdif(ix_first+2:ix_first+3,:,:,:))
	! update host(xdif(:,iz_last-3:iz_last-2,:,:),xdif(:,iz_first+2:iz_first+3,:,:))
	! update host(xdif(:,:,iy_last-3:iy_last-2,:),xdif(:,:,iy_first+2:iy_first+3,:))
!      call mpi_transfersm(xdif,8)	
	! update device(xdif(ix_last-1:ix_last,:,:,:),xdif(ix_first:ix_first+1,:,:,:))
	! update device(xdif(:,iz_last-1:iz_last,:,:),xdif(:,iz_first:iz_first+1,:,:))
	! update device(xdif(:,:,iy_last-1:iy_last,:),xdif(:,:,iy_first:iy_first+1,:))

!      call mpi_transfer8_xdif_openacc
	
      return
      endsubroutine right_openacc
	
!

!hw*************************************************************************
      subroutine current_openacc
      use declare
	use declare_for_openacc
!	integer itag
!     real*8 drby_dx, rby_tmp
!     real*8, dimension(mx,mz,my) :: rby
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

 	
 	!$acc kernels
	!$acc loop independent collapse(3)
      do 10 jy=iy_first,iy_last
      do 10 jz=iz_first,iz_last
      do 10 jx=ix_first,ix_last
      rby(jx,jz,jy)=xx(jx)*x1(jx,jz,jy,7)
   10 continue

! total(mx*mz*9+mx*5+mz*4)=3780

	!$omp parallel default(shared) &
	!$ private(jx,jz,jy,drby_dx,rby_tmp,itag)
	!$omp do
	!$acc loop independent collapse(3)
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

    
    
! current for IR points, calculate without openACC
	
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

    1 continue
	!$acc end kernels
	!$omp end do
    	!$omp end parallel
        
	! update host(cur(ix_last-3:ix_last-2,:,:,:),cur(ix_first+2:ix_first+3,:,:,:))
	! update host(cur(:,iz_last-3:iz_last-2,:,:),cur(:,iz_first+2:iz_first+3,:,:))
	! update host(cur(:,:,iy_last-3:iy_last-2,:),cur(:,:,iy_first+2:iy_first+3,:))
!    	call mpi_transfersm(cur(:,:,:,:),3)
	! update device(cur(ix_last-1:ix_last,:,:,:),cur(ix_first:ix_first+1,:,:,:))
	! update device(cur(:,iz_last-1:iz_last,:,:),cur(:,iz_first:iz_first+1,:,:))
	! update device(cur(:,:,iy_last-1:iy_last,:),cur(:,:,iy_first:iy_first+1,:))

!      call mpi_transfer3_cur_openacc

      return
      endsubroutine current_openacc    	
!


!hw***************************************************************************************
!calculate the xy, xy2 for all points, and x1z, x2z, x1r, x2r for IR points
	subroutine convt_openacc
      use declare
	use declare_for_openacc
!	integer itag, i
!     real*8, dimension(my) :: wwy 
!     real*8, dimension(mx,mz,my) :: x1r_tmp,xr2_tmp,x1z_tmp,xz2_tmp,x1_tmp
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


	!$omp parallel default(shared) &
	!$ private(jx,jz,jy,m)
	!$omp do
	!$acc kernels
	!$acc loop independent collapse(2)
      do 15 m=1,8
	do 15 jy=iy_first+2,iy_last-2
	!$acc loop independent collapse(2)
      do 16 jz=iz_first,iz_last
      do 16 jx=ix_first,ix_last
      xy(jx,jz,jy,m) =d1fc(x1(jx,jz,jy-2,m),x1(jx,jz,jy-1,m),x1(jx,jz,jy,m),x1(jx,jz,jy+1,m),x1(jx,jz,jy+2,m),ay1(jy),by1(jy),cy1(jy),dy1(jy))      
      xy2(jx,jz,jy,m)=d2fc(x1(jx,jz,jy-2,m),x1(jx,jz,jy-1,m),x1(jx,jz,jy,m),x1(jx,jz,jy+1,m),x1(jx,jz,jy+2,m),ay2(jy),by2(jy),cy2(jy),dy2(jy))
   16 continue

! for the xy_8bndx and xy2_8bndx and xy_8bndz and xy2_8bndz
	!$acc loop independent collapse(1)
	do 17 jx=1,nbndx
	xy_8bndx(jx,jy,m)=d1fc(x1_8bndx(jx,jy-2,m),x1_8bndx(jx,jy-1,m),x1_8bndx(jx,jy,m),x1_8bndx(jx,jy+1,m),x1_8bndx(jx,jy+2,m),&
		  ay1(jy),by1(jy),cy1(jy),dy1(jy))
	xy2_8bndx(jx,jy,m)=d2fc(x1_8bndx(jx,jy-2,m),x1_8bndx(jx,jy-1,m),x1_8bndx(jx,jy,m),x1_8bndx(jx,jy+1,m),x1_8bndx(jx,jy+2,m),&
		  ay2(jy),by2(jy),cy2(jy),dy2(jy))
   17 continue

	!$acc loop independent collapse(1)
	do 18 jx=1,nbndz
	xy_8bndz(jx,jy,m)=d1fc(x1_8bndz(jx,jy-2,m),x1_8bndz(jx,jy-1,m),x1_8bndz(jx,jy,m),x1_8bndz(jx,jy+1,m),x1_8bndz(jx,jy+2,m),&
		  ay1(jy),by1(jy),cy1(jy),dy1(jy))
	xy2_8bndz(jx,jy,m)=d2fc(x1_8bndz(jx,jy-2,m),x1_8bndz(jx,jy-1,m),x1_8bndz(jx,jy,m),x1_8bndz(jx,jy+1,m),x1_8bndz(jx,jy+2,m),&
		  ay2(jy),by2(jy),cy2(jy),dy2(jy))
   18 continue
   15 continue
	!$omp end do
    	!$omp end parallel


      !not spectral


	!$omp parallel default(shared) &
	!$ private(jx,jz,jy,m,itag)
	!$omp do
	!$acc loop independent collapse(4)
      do 1 m=1,8
      do 1 jy=iy_first,iy_last
      do 1 jz=iz_first+2,iz_last-2
      do 1 jx=ix_first+2,ix_last-2
      x1r(jx,jz,jy,m)=d1fc(x1(jx-2,jz,jy,m),x1(jx-1,jz,jy,m),x1(jx,jz,jy,m) &
          ,x1(jx+1,jz,jy,m),x1(jx+2,jz,jy,m),ax1(jx),bx1(jx),cx1(jx),dx1(jx))

      xr2(jx,jz,jy,m)=d2fc(x1(jx-2,jz,jy,m),x1(jx-1,jz,jy,m),x1(jx,jz,jy,m) &
          ,x1(jx+1,jz,jy,m),x1(jx+2,jz,jy,m),ax2(jx),bx2(jx),cx2(jx),dx2(jx))

      x1z(jx,jz,jy,m)=d1fc(x1(jx,jz-2,jy,m),x1(jx,jz-1,jy,m),x1(jx,jz,jy,m) &
          ,x1(jx,jz+1,jy,m),x1(jx,jz+2,jy,m),az1(jz),bz1(jz),cz1(jz),dz1(jz))
	
      xz2(jx,jz,jy,m)=d2fc(x1(jx,jz-2,jy,m),x1(jx,jz-1,jy,m),x1(jx,jz,jy,m) &
          ,x1(jx,jz+1,jy,m),x1(jx,jz+2,jy,m),az2(jz),bz2(jz),cz2(jz),dz2(jz))

! can it be vectorization ????????????????????????????????????????????????????????
! d1fc(x1(jx,jz,:,:)) ??
! cost about 25s for total 165s

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

!    1 continue
!	!$omp end do
!    	!$omp end parallel

	! loop independent collapse(4)
!   	do 30 m=1,8
!	do 30 jy=iy_first,iy_last
!	do 30 jz=iz_first,iz_last
!	do 30 jx=ix_first,ix_last
	xr(jx,jz,jy,m)=xint_dx(jx,jz,m)+x1r(jx,jz,jy,m)  
	xz(jx,jz,jy,m)=xint_dz(jx,jz,m)+x1z(jx,jz,jy,m)
 !  30 continue
    1 continue
	!$omp end do
    	!$omp end parallel
	!$acc end kernels

    return
    endsubroutine convt_openacc
!

!hw************************************************************************
!ws************************************************************************
!wzhang************************************************************
      subroutine efield_openacc
      use declare
	use declare_for_openacc
	implicit none
!	integer itag
      include 'mpif.h'

	interface
      subroutine interp1d2l(x1,x2,x3,y1,y2,y3,y,ans)
	!$acc routine seq
      real*8 x1,x2,x3,y1,y2,y3,y,ans
      real*8 d1,d2,d3,tmp_add
	end subroutine
	end interface


	!$omp parallel default(shared) &
	!$ private(jx,jz,jy,m)
	!$omp do
	!$acc kernels
	!$acc loop independent collapse(1)
      do 1 jy=iy_first+2,iy_last-2
	!$acc loop independent collapse(2)
      do 2 jz=iz_first,iz_last
      do 2 jx=ix_first,ix_last

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

    2 continue
      
       if(hall) then
	!$acc loop independent collapse(2)
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
	endif
       

      if(resisitive .and. etaj_in_e) then
	!$acc loop independent collapse(3)
      do m=1,3
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
      ef(jx,jz,jy,m)=ef(jx,jz,jy,m)+eta(jx,jz,jy)*cur(jx,jz,jy,m)
      enddo
      enddo
      enddo
      endif
	
    1 continue

	!$omp end do
    	!$omp end parallel

!      if(bootstrap) then
!      call current_boot(lbs)      
!      do m=1,3
!      ef(:,:,:,m)=ef(:,:,:,m)-eta(:,:,:)*cub(:,:,:,m)
!      enddo
!      endif

      if(nstep.lt.nper) then
	!$acc loop independent collapse(3)
      do 11 m=1,3
      do 11 jy=iy_first,iy_last
      do 11 jz=iz_first,iz_last
      do 11 jx=ix_first,ix_last
      ef(jx,jz,jy,m)=ef(jx,jz,jy,m)+eta1(jx,jz,jy)*(cint(jx,jz,m)+cur(jx,jz,jy,m))
  11  continue
      endif 
	!$acc end kernels

       !***************revised**************************************
	! update host(ef(ix_last-3:ix_last-2,:,:,:),ef(ix_first+2:ix_first+3,:,:,:))
	! update host(ef(:,iz_last-3:iz_last-2,:,:),ef(:,iz_first+2:iz_first+3,:,:))
	! update host(ef(:,:,iy_last-3:iy_last-2,:),ef(:,:,iy_first+2:iy_first+3,:))
!	call mpi_transfersm(ef(:,:,:,:),3)
	! update device(ef(ix_last-1:ix_last,:,:,:),ef(ix_first:ix_first+1,:,:,:))
	! update device(ef(:,iz_last-1:iz_last,:,:),ef(:,iz_first:iz_first+1,:,:))
	! update device(ef(:,:,iy_last-1:iy_last,:),ef(:,:,iy_first:iy_first+1,:))

      call mpi_transfer3_ef_openacc

    ! haowei need to obatin the efield for bnd grids, i.e., ef_3bndz, ef_3bndx
	
    	!$acc kernels
	!$acc loop independent collapse(2)
    	do m=1,3
    	do jy=iy_first,iy_last
	!$acc loop independent collapse(1)
	do itag=1,nbndx
	ef_3bndx(itag,jy,m)=0.d0 ! fixed the ef_3bndx
    	enddo
	!$acc loop independent collapse(1)
	do itag=1,nbndz
	ef_3bndz(itag,jy,m)=0.d0 ! fixed the ef_3bndz
    	enddo
    	enddo
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Boundary for efield
	if(hall) then ! Efield is complete when hall=.false.
	
!	do 3 jz=iz_first,iz_last
!	do 3 jx=ix_first,ix_last
	!$acc loop independent collapse(2)
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
	!$omp parallel default(shared) &
	!$ private(jx,jz,jy,m)
	!$omp do
	!$acc loop independent collapse(4)
	do 10 m=1,3
      do 10 jy=iy_first,iy_last
      do 10 jz=iz_first,iz_last
      do 10 jx=ix_first,ix_last
!	if(psi(jx,jz).lt.psia1) then
!	if((hypb_ratio(jx,jz).ge.0.d0).and.(hypb_ratio(jx,jz).le.1.d0)) then
	ef(jx,jz,jy,m)=ef(jx,jz,jy,m)*hypb_ratio(jx,jz)
!	endif
!	endif
   10 continue
	!$acc end kernels
	!$omp end do
    	!$omp end parallel
!       call mpi_transfersm(ef(:,:,:,:),3)
       !**********************************************************

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
	!$acc parallel default(present)
	!$acc loop independent collapse(4) 
	do 1 m=1,8
	do 1 jy=iy_first,iy_last
	do 1 jz=iz_first,iz_last
	do 1 jx=ix_first,ix_last
       xfold(jx,jz,jy,m)=x(jx,jz,jy,m)
       xm(jx,jz,jy,m)=xfold(jx,jz,jy,m)+xdif(jx,jz,jy,m)*dt/6.
       x(jx,jz,jy,m)=xfold(jx,jz,jy,m)+xdif(jx,jz,jy,m)*dt/2.
!       xfold(:,:,:,:)=x(:,:,:,:)
!       xm(:,:,:,:)=xfold(:,:,:,:)+xdif(:,:,:,:)*dt/6.
!       x(:,:,:,:)=xfold(:,:,:,:)+xdif(:,:,:,:)*dt/2.
    1 continue
	!$acc end parallel
!
      tt=time+dt/2.
      tt1=time+dt/2.
      tt2=time+dt/6.
      irk=2
        call right_openacc
	!$acc parallel default(present)
	!$acc loop independent collapse(4) 
	do 2 m=1,8
	do 2 jy=iy_first,iy_last
	do 2 jz=iz_first,iz_last
	do 2 jx=ix_first,ix_last
       xm(jx,jz,jy,m)=xm(jx,jz,jy,m)+xdif(jx,jz,jy,m)*dt/3.
       x(jx,jz,jy,m)=xfold(jx,jz,jy,m)+xdif(jx,jz,jy,m)*dt/2.
!        xm(:,:,:,:)=xm(:,:,:,:)+xdif(:,:,:,:)*dt/3.
!        x(:,:,:,:)=xfold(:,:,:,:)+xdif(:,:,:,:)*dt/2.
    2 continue
	!$acc end parallel
!
      tt1=time+5.*dt/6.
      tt2=time+dt/2.
      irk=3
        call right_openacc
	!$acc parallel default(present)
	!$acc loop independent collapse(4) 
	do 3 m=1,8
	do 3 jy=iy_first,iy_last
	do 3 jz=iz_first,iz_last
	do 3 jx=ix_first,ix_last
       xm(jx,jz,jy,m)=xm(jx,jz,jy,m)+xdif(jx,jz,jy,m)*dt/3.
       x(jx,jz,jy,m)=xfold(jx,jz,jy,m)+xdif(jx,jz,jy,m)*dt
!        xm(:,:,:,:)=xm(:,:,:,:)+xdif(:,:,:,:)*dt/3.
!        x(:,:,:,:)=xfold(:,:,:,:)+xdif(:,:,:,:)*dt
    3 continue
	!$acc end parallel
!
      time=time+dt
      tt1=time+dt
      tt2=time+5.*dt/6.
      irk=4
        call right_openacc
	!$acc parallel default(present)
	!$acc loop independent collapse(4) 
	do 4 m=1,8
	do 4 jy=iy_first,iy_last
	do 4 jz=iz_first,iz_last
	do 4 jx=ix_first,ix_last
        x(jx,jz,jy,m)=xm(jx,jz,jy,m)+xdif(jx,jz,jy,m)*dt/6.
!        x(:,:,:,:)=xm(:,:,:,:)+xdif(:,:,:,:)*dt/6.
    4 continue
	!$acc end parallel

      caf=0.75d0*(0.5+0.5*dtanh((time-40)/5.))
!      call bndry_x_ex(lbnd)
!      call bndry8_x_ex(lbnd)
	call bndry8_cut_cell_v2_fixed
  !    if(conductpll) call pllconduct(lpll)
!      if(smoothpll) call smthp_traceline_5p(1)
!      if(eta_from_t) call calculate_eta
      
      return
      endsubroutine stepon_openacc
	



!ws****************************************************************************************
      subroutine mpi_transfer8_xdif_openacc
      use declare
	use declare_for_openacc
!      real*8, dimension(mx,mz,my,8) :: w8
!      real*8, dimension(mz,my,8) :: wfx1,wfx2
!      real*8, dimension(mx,my,8) :: wfz1,wfz2
!      real*8, dimension(mx,mz,8) :: wfy1,wfy2
      include 'mpif.h'

! send xdif up unless i'm at the top, then receive from below
	! data create(wfx1,wfx2,wfz1,wfz2,wfy1,wfy2)
     
	if (nrankxz.ge.nrkz(nrank)*nprx .and. nrankxz.lt.nrkz(nrank)*nprx+nprx-1) then
	!$acc kernels present(xdif,wfx1,wfx2,wfy1,wfy2,wfz1,wfz2)
      wfx1(:,:,:)=xdif(ix_last-2,:,:,:)
      wfx2(:,:,:)=xdif(ix_last-3,:,:,:)
	!$acc end kernels
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wfx1,wfx2) if_present
	call mpi_send( wfx1, myz8, mpi_double_precision, nrank + 1, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wfx2, myz8, mpi_double_precision, nrank + 1, 0,  &
		      mpi_comm_world,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	endif
	if (nrankxz.gt.nrkz(nrank)*nprx .and. nrankxz.le.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wfx1,wfx2) if_present
	call mpi_recv( wfx1, myz8, mpi_double_precision, nrank - 1, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wfx2, myz8, mpi_double_precision, nrank - 1, 0,  &
		      mpi_comm_world, status,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	!$acc kernels present(xdif,wfx1,wfx2,wfy1,wfy2,wfz1,wfz2)
      xdif(ix_first+1,:,:,:)=wfx1(:,:,:)
      xdif(ix_first,:,:,:)=wfx2(:,:,:)
	!$acc end kernels
	endif
	
      
! send xdif down unless i'm at the bottom

      
	if (nrankxz.gt.nrkz(nrank)*nprx .and. nrankxz.le.nrkz(nrank)*nprx+nprx-1) then
	!$acc kernels present(xdif,wfx1,wfx2,wfy1,wfy2,wfz1,wfz2)
      wfx1(:,:,:)=xdif(ix_first+2,:,:,:)
      wfx2(:,:,:)=xdif(ix_first+3,:,:,:)
	!$acc end kernels
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wfx1,wfx2) if_present
	call mpi_send( wfx1, myz8, mpi_double_precision, nrank - 1, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wfx2, myz8, mpi_double_precision, nrank - 1, 1,  &
		      mpi_comm_world,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrankxz.ge.nrkz(nrank)*nprx .and. nrankxz.lt.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wfx1,wfx2) if_present
	call mpi_recv( wfx1, myz8, mpi_double_precision, nrank + 1, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wfx2, myz8, mpi_double_precision, nrank + 1, 1,  &
		      mpi_comm_world, status,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	!$acc kernels present(xdif,wfx1,wfx2,wfy1,wfy2,wfz1,wfz2)
      xdif(ix_last-1,:,:,:)=wfx1(:,:,:)
      xdif(ix_last,:,:,:)=wfx2(:,:,:)
	!$acc end kernels
	endif
!
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	if (nrankxz.lt.nsizexz-nprx) then
	!$acc kernels present(xdif,wfx1,wfx2,wfy1,wfy2,wfz1,wfz2)
      wfz1(:,:,:)=xdif(:,iz_last-2,:,:)
      wfz2(:,:,:)=xdif(:,iz_last-3,:,:)
	!$acc end kernels
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wfz1,wfz2) if_present
	call mpi_send( wfz1, myx8, mpi_double_precision, nrank + nprx, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wfz2, myx8, mpi_double_precision, nrank + nprx, 0,  &
		      mpi_comm_world,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	endif
	if (nrankxz.ge.nprx ) then
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wfz1,wfz2) if_present
	call mpi_recv( wfz1, myx8, mpi_double_precision, nrank - nprx, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wfz2, myx8, mpi_double_precision, nrank - nprx, 0,  &
		      mpi_comm_world, status,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	!$acc kernels present(xdif,wfx1,wfx2,wfy1,wfy2,wfz1,wfz2)
      xdif(:,iz_first+1,:,:)=wfz1(:,:,:)
      xdif(:,iz_first,:,:)=wfz2(:,:,:)
	!$acc end kernels
	endif
	     
! send xdif down unless i'm at the bottom

      
	if (nrankxz.ge.nprx ) then
	!$acc kernels present(xdif,wfx1,wfx2,wfy1,wfy2,wfz1,wfz2)
      wfz1(:,:,:)=xdif(:,iz_first+2,:,:)
      wfz2(:,:,:)=xdif(:,iz_first+3,:,:)
	!$acc end kernels
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wfz1,wfz2) if_present
	call mpi_send( wfz1, myx8, mpi_double_precision, nrank - nprx, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wfz2, myx8, mpi_double_precision, nrank - nprx, 1,  &
		      mpi_comm_world,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrankxz.lt.nsizexz-nprx) then
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wfz1,wfz2) if_present
	call mpi_recv( wfz1, myx8, mpi_double_precision, nrank + nprx, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wfz2, myx8, mpi_double_precision, nrank + nprx, 1,  &
		      mpi_comm_world, status,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	!$acc kernels present(xdif,wfx1,wfx2,wfy1,wfy2,wfz1,wfz2)
      xdif(:,iz_last-1,:,:)=wfz1(:,:,:)
      xdif(:,iz_last,:,:)=wfz2(:,:,:)
	!$acc end kernels
	endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      if (npry .gt. 1) then
	if (nrky(nrank).lt. npry-1) then
	!$acc kernels present(xdif,wfx1,wfx2,wfy1,wfy2,wfz1,wfz2)
      wfy1(:,:,:)=xdif(:,:,iy_last-2,:)
      wfy2(:,:,:)=xdif(:,:,iy_last-3,:)
	!$acc end kernels
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wfy1,wfy2) if_present
	call mpi_send( wfy1, mxz8, mpi_double_precision, nrank + nprxz, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wfy2, mxz8, mpi_double_precision, nrank + nprxz, 0,  &
		      mpi_comm_world,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	endif
	if (nrky(nrank) .ge. 1 )  then
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wfy1,wfy2) if_present
	call mpi_recv( wfy1, mxz8, mpi_double_precision, nrank - nprxz, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wfy2, mxz8, mpi_double_precision, nrank - nprxz, 0,  &
		      mpi_comm_world, status,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	!$acc kernels present(xdif,wfx1,wfx2,wfy1,wfy2,wfz1,wfz2)
      xdif(:,:,iy_first+1,:)=wfy1(:,:,:)
      xdif(:,:,iy_first,:)=wfy2(:,:,:)
	!$acc end kernels
	endif


  	if (nrky(nrank).eq. npry-1) then
	!$acc kernels present(xdif,wfx1,wfx2,wfy1,wfy2,wfz1,wfz2)
      wfy1(:,:,:)=xdif(:,:,iy_last-2,:)
      wfy2(:,:,:)=xdif(:,:,iy_last-3,:)
	!$acc end kernels
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wfy1,wfy2) if_present
	call mpi_send( wfy1, mxz8, mpi_double_precision, nrank -(npry-1)*nprxz, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wfy2, mxz8, mpi_double_precision, nrank -(npry-1)*nprxz, 0,  &
		      mpi_comm_world,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	endif
	if (nrky(nrank) .eq. 0 )  then
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wfy1,wfy2) if_present
	call mpi_recv( wfy1, mxz8, mpi_double_precision, nrank +(npry-1)*nprxz, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wfy2, mxz8, mpi_double_precision, nrank +(npry-1)*nprxz, 0,  &
		      mpi_comm_world, status,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	!$acc kernels present(xdif,wfx1,wfx2,wfy1,wfy2,wfz1,wfz2)
      xdif(:,:,iy_first+1,:)=wfy1(:,:,:)
      xdif(:,:,iy_first,:)=wfy2(:,:,:)
	!$acc end kernels
	endif
	     
! send xdif down unless i'm at the bottom

      
	if (nrky(nrank) .ge. 1 ) then
	!$acc kernels present(xdif,wfx1,wfx2,wfy1,wfy2,wfz1,wfz2)
      wfy1(:,:,:)=xdif(:,:,iy_first+2,:)
      wfy2(:,:,:)=xdif(:,:,iy_first+3,:)
	!$acc end kernels
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wfy1,wfy2) if_present
	call mpi_send( wfy1, mxz8, mpi_double_precision, nrank - nprxz, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wfy2, mxz8, mpi_double_precision, nrank - nprxz, 1,  &
		      mpi_comm_world,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrky(nrank).lt. npry-1) then
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wfy1,wfy2) if_present
	call mpi_recv( wfy1, mxz8, mpi_double_precision, nrank + nprxz, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wfy2, mxz8, mpi_double_precision, nrank + nprxz, 1,  &
		      mpi_comm_world, status,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	!$acc kernels present(xdif,wfx1,wfx2,wfy1,wfy2,wfz1,wfz2)
      xdif(:,:,iy_last-1,:)=wfy1(:,:,:)
      xdif(:,:,iy_last,:)=wfy2(:,:,:)
	!$acc end kernels
	endif

    if (nrky(nrank) .eq. 0 ) then
	!$acc kernels present(xdif,wfx1,wfx2,wfy1,wfy2,wfz1,wfz2)
      wfy1(:,:,:)=xdif(:,:,iy_first+2,:)
      wfy2(:,:,:)=xdif(:,:,iy_first+3,:)
	!$acc end kernels
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wfy1,wfy2) if_present
	call mpi_send( wfy1, mxz8, mpi_double_precision, nrank +(npry-1)*nprxz, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wfy2, mxz8, mpi_double_precision, nrank +(npry-1)*nprxz, 1,  &
		      mpi_comm_world,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrky(nrank).eq. npry-1) then
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wfy1,wfy2) if_present
	call mpi_recv( wfy1, mxz8, mpi_double_precision, nrank -(npry-1)*nprxz, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wfy2, mxz8, mpi_double_precision, nrank -(npry-1)*nprxz, 1,  &
		      mpi_comm_world, status,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	!$acc kernels present(xdif,wfx1,wfx2,wfy1,wfy2,wfz1,wfz2)
      xdif(:,:,iy_last-1,:)=wfy1(:,:,:)
      xdif(:,:,iy_last,:)=wfy2(:,:,:)
	!$acc end kernels
	endif
      else
	!$acc kernels present(xdif,wfx1,wfx2,wfy1,wfy2,wfz1,wfz2)
      xdif(:,:,iy_first+1,:)=xdif(:,:,iy_last-2,:)
      xdif(:,:,iy_first,:)=xdif(:,:,iy_last-3,:)
      xdif(:,:,iy_last-1,:)=xdif(:,:,iy_first+2,:)
      xdif(:,:,iy_last,:)=xdif(:,:,iy_first+3,:)
	!$acc end kernels
      endif
	! end data

      return
      end 


!hw**************************************************************



!ws****************************************************************************************
      subroutine mpi_transfer3_ef_openacc
      use declare
	use declare_for_openacc
!      real*8, dimension(mx,mz,my,3) :: ef
!      real*8, dimension(mz,my,3) :: wcx1,wcx2
!      real*8, dimension(mx,my,3) :: wcz1,wcz2
!      real*8, dimension(mx,mz,3) :: wcy1,wcy2
      include 'mpif.h'

       
! send ef up unless i'm at the top, then receive from below
     
	if (nrankxz.ge.nrkz(nrank)*nprx .and. nrankxz.lt.nrkz(nrank)*nprx+nprx-1) then
	!$acc kernels present(ef,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      wcx1(:,:,:)=ef(ix_last-2,:,:,:)
      wcx2(:,:,:)=ef(ix_last-3,:,:,:)
	!$acc end kernels
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcx1,wcx2) if_present
	call mpi_send( wcx1, myz3, mpi_double_precision, nrank + 1, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wcx2, myz3, mpi_double_precision, nrank + 1, 0,  &
		      mpi_comm_world,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	endif
	if (nrankxz.gt.nrkz(nrank)*nprx .and. nrankxz.le.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcx1,wcx2) if_present
	call mpi_recv( wcx1, myz3, mpi_double_precision, nrank - 1, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wcx2, myz3, mpi_double_precision, nrank - 1, 0,  &
		      mpi_comm_world, status,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	!$acc kernels present(ef,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      ef(ix_first+1,:,:,:)=wcx1(:,:,:)
      ef(ix_first,:,:,:)=wcx2(:,:,:)
	!$acc end kernels
	endif
	
      
! send ef down unless i'm at the bottom

      
	if (nrankxz.gt.nrkz(nrank)*nprx .and. nrankxz.le.nrkz(nrank)*nprx+nprx-1) then
	!$acc kernels present(ef,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      wcx1(:,:,:)=ef(ix_first+2,:,:,:)
      wcx2(:,:,:)=ef(ix_first+3,:,:,:)
	!$acc end kernels
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcx1,wcx2) if_present
	call mpi_send( wcx1, myz3, mpi_double_precision, nrank - 1, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wcx2, myz3, mpi_double_precision, nrank - 1, 1,  &
		      mpi_comm_world,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrankxz.ge.nrkz(nrank)*nprx .and. nrankxz.lt.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcx1,wcx2) if_present
	call mpi_recv( wcx1, myz3, mpi_double_precision, nrank + 1, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wcx2, myz3, mpi_double_precision, nrank + 1, 1,  &
		      mpi_comm_world, status,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	!$acc kernels present(ef,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      ef(ix_last-1,:,:,:)=wcx1(:,:,:)
      ef(ix_last,:,:,:)=wcx2(:,:,:)
	!$acc end kernels
	endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	if (nrankxz.lt.nsizexz-nprx) then
	!$acc kernels present(ef,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      wcz1(:,:,:)=ef(:,iz_last-2,:,:)
      wcz2(:,:,:)=ef(:,iz_last-3,:,:)
	!$acc end kernels
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcz1,wcz2) if_present
	call mpi_send( wcz1, myx3, mpi_double_precision, nrank + nprx, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wcz2, myx3, mpi_double_precision, nrank + nprx, 0,  &
		      mpi_comm_world,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	endif
	if (nrankxz.ge.nprx ) then
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcz1,wcz2) if_present
	call mpi_recv( wcz1, myx3, mpi_double_precision, nrank - nprx, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wcz2, myx3, mpi_double_precision, nrank - nprx, 0,  &
		      mpi_comm_world, status,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	!$acc kernels present(ef,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      ef(:,iz_first+1,:,:)=wcz1(:,:,:)
      ef(:,iz_first,:,:)=wcz2(:,:,:)
	!$acc end kernels
	endif
	     
! send ef down unless i'm at the bottom

      
	if (nrankxz.ge.nprx ) then
	!$acc kernels present(ef,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      wcz1(:,:,:)=ef(:,iz_first+2,:,:)
      wcz2(:,:,:)=ef(:,iz_first+3,:,:)
	!$acc end kernels
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcz1,wcz2) if_present
	call mpi_send( wcz1, myx3, mpi_double_precision, nrank - nprx, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wcz2, myx3, mpi_double_precision, nrank - nprx, 1,  &
		      mpi_comm_world,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrankxz.lt.nsizexz-nprx) then
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcz1,wcz2) if_present
	call mpi_recv( wcz1, myx3, mpi_double_precision, nrank + nprx, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wcz2, myx3, mpi_double_precision, nrank + nprx, 1,  &
		      mpi_comm_world, status,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	!$acc kernels present(ef,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      ef(:,iz_last-1,:,:)=wcz1(:,:,:)
      ef(:,iz_last,:,:)=wcz2(:,:,:)
	!$acc end kernels
	endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      if(npry .gt. 1) then
	if (nrky(nrank).lt. npry-1) then
	!$acc kernels present(ef,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      wcy1(:,:,:)=ef(:,:,iy_last-2,:)
      wcy2(:,:,:)=ef(:,:,iy_last-3,:)
	!$acc end kernels
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcy1,wcy2) if_present
	call mpi_send( wcy1, mxz3, mpi_double_precision, nrank + nprxz, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wcy2, mxz3, mpi_double_precision, nrank + nprxz, 0,  &
		      mpi_comm_world,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	endif
	if (nrky(nrank) .ge. 1 )  then
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcy1,wcy2) if_present
	call mpi_recv( wcy1, mxz3, mpi_double_precision, nrank - nprxz, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wcy2, mxz3, mpi_double_precision, nrank - nprxz, 0,  &
		      mpi_comm_world, status,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	!$acc kernels present(ef,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      ef(:,:,iy_first+1,:)=wcy1(:,:,:)
      ef(:,:,iy_first,:)=wcy2(:,:,:)
	!$acc end kernels
	endif


  	if (nrky(nrank).eq. npry-1) then
	!$acc kernels present(ef,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      wcy1(:,:,:)=ef(:,:,iy_last-2,:)
      wcy2(:,:,:)=ef(:,:,iy_last-3,:)
	!$acc end kernels
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcy1,wcy2) if_present
	call mpi_send( wcy1, mxz3, mpi_double_precision, nrank -(npry-1)*nprxz, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wcy2, mxz3, mpi_double_precision, nrank -(npry-1)*nprxz, 0,  &
		      mpi_comm_world,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	endif
	if (nrky(nrank) .eq. 0 )  then
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcy1,wcy2) if_present
	call mpi_recv( wcy1, mxz3, mpi_double_precision, nrank +(npry-1)*nprxz, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wcy2, mxz3, mpi_double_precision, nrank +(npry-1)*nprxz, 0,  &
		      mpi_comm_world, status,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	!$acc kernels present(ef,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      ef(:,:,iy_first+1,:)=wcy1(:,:,:)
      ef(:,:,iy_first,:)=wcy2(:,:,:)
	!$acc end kernels
	endif
	     
! send ef down unless i'm at the bottom

      
	if (nrky(nrank) .ge. 1 ) then
	!$acc kernels present(ef,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      wcy1(:,:,:)=ef(:,:,iy_first+2,:)
      wcy2(:,:,:)=ef(:,:,iy_first+3,:)
	!$acc end kernels
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcy1,wcy2) if_present
	call mpi_send( wcy1, mxz3, mpi_double_precision, nrank - nprxz, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wcy2, mxz3, mpi_double_precision, nrank - nprxz, 1,  &
		      mpi_comm_world,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrky(nrank).lt. npry-1) then
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcy1,wcy2) if_present
	call mpi_recv( wcy1, mxz3, mpi_double_precision, nrank + nprxz, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wcy2, mxz3, mpi_double_precision, nrank + nprxz, 1,  &
		      mpi_comm_world, status,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	!$acc kernels present(ef,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      ef(:,:,iy_last-1,:)=wcy1(:,:,:)
      ef(:,:,iy_last,:)=wcy2(:,:,:)
	!$acc end kernels
	endif

    if (nrky(nrank) .eq. 0 ) then
	!$acc kernels present(ef,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      wcy1(:,:,:)=ef(:,:,iy_first+2,:)
      wcy2(:,:,:)=ef(:,:,iy_first+3,:)
	!$acc end kernels
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcy1,wcy2) if_present
	call mpi_send( wcy1, mxz3, mpi_double_precision, nrank +(npry-1)*nprxz, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wcy2, mxz3, mpi_double_precision, nrank +(npry-1)*nprxz, 1,  &
		      mpi_comm_world,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrky(nrank).eq. npry-1) then
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcy1,wcy2) if_present
	call mpi_recv( wcy1, mxz3, mpi_double_precision, nrank -(npry-1)*nprxz, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wcy2, mxz3, mpi_double_precision, nrank -(npry-1)*nprxz, 1,  &
		      mpi_comm_world, status,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	!$acc kernels present(ef,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      ef(:,:,iy_last-1,:)=wcy1(:,:,:)
      ef(:,:,iy_last,:)=wcy2(:,:,:)
	!$acc end kernels
	endif
      else
	!$acc kernels present(ef,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      ef(:,:,iy_first+1,:)=ef(:,:,iy_last-2,:)
      ef(:,:,iy_first,:)=ef(:,:,iy_last-3,:)
      ef(:,:,iy_last-1,:)=ef(:,:,iy_first+2,:)
      ef(:,:,iy_last,:)=ef(:,:,iy_first+3,:)
	!$acc end kernels
      endif
      return
      end


!ws****************************************************************************************
      subroutine mpi_transfer3_cur_openacc
      use declare
	use declare_for_openacc
!      real*8, dimension(mx,mz,my,3) :: cur
!      real*8, dimension(mz,my,3) :: wcx1,wcx2
!      real*8, dimension(mx,my,3) :: wcz1,wcz2
!      real*8, dimension(mx,mz,3) :: wcy1,wcy2
      include 'mpif.h'

       
! send cur up unless i'm at the top, then receive from below
     
	if (nrankxz.ge.nrkz(nrank)*nprx .and. nrankxz.lt.nrkz(nrank)*nprx+nprx-1) then
	!$acc kernels present(cur,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      wcx1(:,:,:)=cur(ix_last-2,:,:,:)
      wcx2(:,:,:)=cur(ix_last-3,:,:,:)
	!$acc end kernels
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcx1,wcx2) if_present
	call mpi_send( wcx1, myz3, mpi_double_precision, nrank + 1, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wcx2, myz3, mpi_double_precision, nrank + 1, 0,  &
		      mpi_comm_world,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	endif
	if (nrankxz.gt.nrkz(nrank)*nprx .and. nrankxz.le.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcx1,wcx2) if_present
	call mpi_recv( wcx1, myz3, mpi_double_precision, nrank - 1, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wcx2, myz3, mpi_double_precision, nrank - 1, 0,  &
		      mpi_comm_world, status,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	!$acc kernels present(cur,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      cur(ix_first+1,:,:,:)=wcx1(:,:,:)
      cur(ix_first,:,:,:)=wcx2(:,:,:)
	!$acc end kernels
	endif
	
      
! send cur down unless i'm at the bottom

      
	if (nrankxz.gt.nrkz(nrank)*nprx .and. nrankxz.le.nrkz(nrank)*nprx+nprx-1) then
	!$acc kernels present(cur,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      wcx1(:,:,:)=cur(ix_first+2,:,:,:)
      wcx2(:,:,:)=cur(ix_first+3,:,:,:)
	!$acc end kernels
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcx1,wcx2) if_present
	call mpi_send( wcx1, myz3, mpi_double_precision, nrank - 1, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wcx2, myz3, mpi_double_precision, nrank - 1, 1,  &
		      mpi_comm_world,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrankxz.ge.nrkz(nrank)*nprx .and. nrankxz.lt.nrkz(nrank)*nprx+nprx-1) then
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcx1,wcx2) if_present
	call mpi_recv( wcx1, myz3, mpi_double_precision, nrank + 1, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wcx2, myz3, mpi_double_precision, nrank + 1, 1,  &
		      mpi_comm_world, status,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	!$acc kernels present(cur,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      cur(ix_last-1,:,:,:)=wcx1(:,:,:)
      cur(ix_last,:,:,:)=wcx2(:,:,:)
	!$acc end kernels
	endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	if (nrankxz.lt.nsizexz-nprx) then
	!$acc kernels present(cur,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      wcz1(:,:,:)=cur(:,iz_last-2,:,:)
      wcz2(:,:,:)=cur(:,iz_last-3,:,:)
	!$acc end kernels
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcz1,wcz2) if_present
	call mpi_send( wcz1, myx3, mpi_double_precision, nrank + nprx, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wcz2, myx3, mpi_double_precision, nrank + nprx, 0,  &
		      mpi_comm_world,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	endif
	if (nrankxz.ge.nprx ) then
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcz1,wcz2) if_present
	call mpi_recv( wcz1, myx3, mpi_double_precision, nrank - nprx, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wcz2, myx3, mpi_double_precision, nrank - nprx, 0,  &
		      mpi_comm_world, status,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	!$acc kernels present(cur,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      cur(:,iz_first+1,:,:)=wcz1(:,:,:)
      cur(:,iz_first,:,:)=wcz2(:,:,:)
	!$acc end kernels
	endif
	     
! send cur down unless i'm at the bottom

      
	if (nrankxz.ge.nprx ) then
	!$acc kernels present(cur,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      wcz1(:,:,:)=cur(:,iz_first+2,:,:)
      wcz2(:,:,:)=cur(:,iz_first+3,:,:)
	!$acc end kernels
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcz1,wcz2) if_present
	call mpi_send( wcz1, myx3, mpi_double_precision, nrank - nprx, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wcz2, myx3, mpi_double_precision, nrank - nprx, 1,  &
		      mpi_comm_world,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrankxz.lt.nsizexz-nprx) then
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcz1,wcz2) if_present
	call mpi_recv( wcz1, myx3, mpi_double_precision, nrank + nprx, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wcz2, myx3, mpi_double_precision, nrank + nprx, 1,  &
		      mpi_comm_world, status,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	!$acc kernels present(cur,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      cur(:,iz_last-1,:,:)=wcz1(:,:,:)
      cur(:,iz_last,:,:)=wcz2(:,:,:)
	!$acc end kernels
	endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      if(npry .gt. 1) then
	if (nrky(nrank).lt. npry-1) then
	!$acc kernels present(cur,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      wcy1(:,:,:)=cur(:,:,iy_last-2,:)
      wcy2(:,:,:)=cur(:,:,iy_last-3,:)
	!$acc end kernels
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcy1,wcy2) if_present
	call mpi_send( wcy1, mxz3, mpi_double_precision, nrank + nprxz, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wcy2, mxz3, mpi_double_precision, nrank + nprxz, 0,  &
		      mpi_comm_world,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	endif
	if (nrky(nrank) .ge. 1 )  then
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcy1,wcy2) if_present
	call mpi_recv( wcy1, mxz3, mpi_double_precision, nrank - nprxz, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wcy2, mxz3, mpi_double_precision, nrank - nprxz, 0,  &
		      mpi_comm_world, status,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	!$acc kernels present(cur,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      cur(:,:,iy_first+1,:)=wcy1(:,:,:)
      cur(:,:,iy_first,:)=wcy2(:,:,:)
	!$acc end kernels
	endif


  	if (nrky(nrank).eq. npry-1) then
	!$acc kernels present(cur,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      wcy1(:,:,:)=cur(:,:,iy_last-2,:)
      wcy2(:,:,:)=cur(:,:,iy_last-3,:)
	!$acc end kernels
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcy1,wcy2) if_present
	call mpi_send( wcy1, mxz3, mpi_double_precision, nrank -(npry-1)*nprxz, 0,  &
		      mpi_comm_world,ierror )
	call mpi_send( wcy2, mxz3, mpi_double_precision, nrank -(npry-1)*nprxz, 0,  &
		      mpi_comm_world,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	endif
	if (nrky(nrank) .eq. 0 )  then
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcy1,wcy2) if_present
	call mpi_recv( wcy1, mxz3, mpi_double_precision, nrank +(npry-1)*nprxz, 0,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wcy2, mxz3, mpi_double_precision, nrank +(npry-1)*nprxz, 0,  &
		      mpi_comm_world, status,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	!$acc kernels present(cur,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      cur(:,:,iy_first+1,:)=wcy1(:,:,:)
      cur(:,:,iy_first,:)=wcy2(:,:,:)
	!$acc end kernels
	endif
	     
! send cur down unless i'm at the bottom

      
	if (nrky(nrank) .ge. 1 ) then
	!$acc kernels present(cur,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      wcy1(:,:,:)=cur(:,:,iy_first+2,:)
      wcy2(:,:,:)=cur(:,:,iy_first+3,:)
	!$acc end kernels
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcy1,wcy2) if_present
	call mpi_send( wcy1, mxz3, mpi_double_precision, nrank - nprxz, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wcy2, mxz3, mpi_double_precision, nrank - nprxz, 1,  &
		      mpi_comm_world,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrky(nrank).lt. npry-1) then
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcy1,wcy2) if_present
	call mpi_recv( wcy1, mxz3, mpi_double_precision, nrank + nprxz, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wcy2, mxz3, mpi_double_precision, nrank + nprxz, 1,  &
		      mpi_comm_world, status,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	!$acc kernels present(cur,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      cur(:,:,iy_last-1,:)=wcy1(:,:,:)
      cur(:,:,iy_last,:)=wcy2(:,:,:)
	!$acc end kernels
	endif

    if (nrky(nrank) .eq. 0 ) then
	!$acc kernels present(cur,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      wcy1(:,:,:)=cur(:,:,iy_first+2,:)
      wcy2(:,:,:)=cur(:,:,iy_first+3,:)
	!$acc end kernels
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcy1,wcy2) if_present
	call mpi_send( wcy1, mxz3, mpi_double_precision, nrank +(npry-1)*nprxz, 1,  &
		      mpi_comm_world,ierror )
	call mpi_send( wcy2, mxz3, mpi_double_precision, nrank +(npry-1)*nprxz, 1,  &
		      mpi_comm_world,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	endif
		      
	if (nrky(nrank).eq. npry-1) then
!mpi   ----------------------------------------------------------------
	!$acc host_data use_device(wcy1,wcy2) if_present
	call mpi_recv( wcy1, mxz3, mpi_double_precision, nrank -(npry-1)*nprxz, 1,  &
		      mpi_comm_world, status,ierror )
	call mpi_recv( wcy2, mxz3, mpi_double_precision, nrank -(npry-1)*nprxz, 1,  &
		      mpi_comm_world, status,ierror )
	!$acc end host_data
!mpi   ----------------------------------------------------------------
	!$acc kernels present(cur,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      cur(:,:,iy_last-1,:)=wcy1(:,:,:)
      cur(:,:,iy_last,:)=wcy2(:,:,:)
	!$acc end kernels
	endif
      else
	!$acc kernels present(cur,wcx1,wcx2,wcy1,wcy2,wcz1,wcz2)
      cur(:,:,iy_first+1,:)=cur(:,:,iy_last-2,:)
      cur(:,:,iy_first,:)=cur(:,:,iy_last-3,:)
      cur(:,:,iy_last-1,:)=cur(:,:,iy_first+2,:)
      cur(:,:,iy_last,:)=cur(:,:,iy_first+3,:)
	!$acc end kernels
      endif
      return
      end
