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




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!      !!!!!!!          !!!!!!!!!!!!!    !!!!
!!!!!  !!!  !!!!!!  !!!!!!!!!!!!!!!!!!  !!!!!!!!!
!!!!!  !!!! !!!!!!  !!!!!!!!!!!!!!!!   !!!!!!!!!!
!!!!!  !!  !!!!!!!        !!!!!!!!!   !!!!!!!!!!!
!!!!!  !!! !!!!!!!  !!!!!!!!!!!!!!!!   !!!!!!!!!!
!!!!!  !!!! !!!!!!  !!!!!!!!!!!!!!!!!    !!!!!!!!
!!!!!  !!!!! !!!!!          !!!!!!!!!!!!    !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module for record subroutines.
! written by S. WANG
! sortted by H.W. Zhang during 2018 AUG.
! contact changhw@zju.edu.cn if there is any problem.

! including subroutines as follows
! recrd1
! recrd10
! recrd100
! recrd
! recrd_rmp_vacuum (record the total vacuum RMP field)
! recrd_rmp_vacuum_2 (record the added outside RMP field)
! recrd_init
! recrd_dssp
! recrd_dbg
! recrd_per
! recrd_xdif
! recrd_cv
! recrd_x1st
! readin



!ws************************************************************************
      subroutine recrd1
      use declare
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

!ws*********************************************
      subroutine recrd10
      use declare
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

!ws*********************************************      
      subroutine recrd100
      use declare
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
      use declare
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
!hw************************************************************************
      subroutine recrd_rmp_vacuum
      use declare
      include 'mpif.h'
!
      character*9 output
      character*3 cn
      character*3 cn1

	x(:,:,:,6:8)=x(:,:,:,6:8)+b_rmp(:,:,:,1:3)
      
      output='tk'//cn1(nrank)//cn(0)
      open(unit=7,file=output,status='unknown',form='unformatted')
      write(7)ncase,nstep,time
      write(7)x
      close(7)

	x(:,:,:,6:8)=x(:,:,:,6:8)-b_rmp(:,:,:,1:3)
    
      return
      end

!hw************************************************************************
      subroutine recrd_rmp_vacuum_2
      use declare
      include 'mpif.h'
!
      character*9 output
      character*3 cn
      character*3 cn1

	x(:,:,:,6:8)=x(:,:,:,6:8)+b_rmp_out(:,:,:,1:3)
      
      output='tk'//cn1(nrank)//cn(1)
      open(unit=7,file=output,status='unknown',form='unformatted')
      write(7)ncase,nstep,time
      write(7)x
      close(7)

	x(:,:,:,6:8)=x(:,:,:,6:8)-b_rmp_out(:,:,:,1:3)
    
      return
      end

!ws*********************************************	
     subroutine recrd_init
      use declare
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

!ws*********************************************	
      subroutine recrd_dssp
      use declare
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

!ws*********************************************
      subroutine recrd_dbg
      use declare
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
      write(7,200)(((x1(jx,jz,jy,i),i=1,8),br(jx,jz,jy),bp(jx,jz,jy),vr(jx,jz,jy),vp(jx,jz,jy),vr(jx,jz,jy)*xx(jx)*bp0(jx,jz),jx=ix_first,ix_last),jz=iz_first,iz_last)
 200  format(13(1x,e12.5)) 
      close(7)

      output='xy_x'//cn1(nrank)//cn(nst)
      open(unit=17,file=output,status='unknown',form='formatted') 
      jz=mz/2   
      write(17,200)(((x1(jx,jz,jy,i),i=1,8),br(jx,jz,jy),bp(jx,jz,jy),vr(jx,jz,jy),vp(jx,jz,jy),vr(jx,jz,jy)*xx(jx)*bp0(jx,jz),jx=ix_first,ix_last),jy=1,my)
      close(17)
!      output='xrt'//cn1(nrank)//cn(nst)
!      open(unit=8,file=output,status='unknown',form='formatted')
!      write(8,300)(((xrt(jr,jt,1,i)-xrtint(jr,jt,i),i=1,8),jr=1,mr),jt=1,mt)
! 300  format(8(1x,e12.5)) 
!      close(8)
      output='cur'//cn1(nrank)//cn(nst)
      open(unit=9,file=output,status='unknown',form='formatted')
      write(9,400)(((cur(jx,jz,1,i),i=1,3),cr(jx,jz,1),cp(jx,jz,1),jx=ix_first,ix_last),jz=iz_first,iz_last)
 400  format(5(1x,e12.5)) 
      close(9)
      output='xy_c'//cn1(nrank)//cn(nst)
      open(unit=19,file=output,status='unknown',form='formatted')
      jz=mz/2
      write(19,400)(((cur(jx,jz,jy,i),i=1,3),cr(jx,jz,jy),cp(jx,jz,jy),jx=ix_first,ix_last),jy=1,my)
      close(19)
      output='ef'//cn1(nrank)//cn(nst)
      open(unit=10,file=output,status='unknown',form='formatted')
      write(10,300)(((ef(jx,jz,1,i),i=1,3),jx=ix_first,ix_last),jz=iz_first,iz_last)
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
      use declare
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
      write(7,200)((((x1(jx,jz,jy,i),i=1,8),br(jx,jz,jy),bp(jx,jz,jy),vr(jx,jz,jy),vp(jx,jz,jy),vr(jx,jz,jy)*xx(jx)*bp0(jx,jz),jx=ix_first,ix_last),jz=iz_first,iz_last),jy=1,my)
 200  format(13(1x,e12.5)) 
      close(7)

      output='per_cur'
      open(unit=9,file=output,status='unknown',form='formatted')
      write(9,400)((((cur(jx,jz,jy,i),i=1,3),cr(jx,jz,jy),cp(jx,jz,jy),jx=ix_first,ix_last),jz=iz_first,iz_last),jy=1,my)
 400  format(5(1x,e12.5)) 
      close(9)
      return
      end

!ws***********************************************
     subroutine recrd_xdif
      use declare
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
!      output='efd'//cn(nst)//cn1(irk)
!      open(unit=72,file=output,status='unknown',form='formatted')
!      write(72,200)((((ef(jx,jz,1,i),efx(jx,jz,1,i),efz(jx,jz,1,i)),i=1,3),jx=ix_first,ix_last),jz=iz_first,iz_last)
!      close(72)
      return
      end

!ws*********************************************************
     subroutine recrd_cv
      use declare
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
 
!ws********************************************* 
     subroutine recrd_x1st
      use declare
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
      use declare
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



