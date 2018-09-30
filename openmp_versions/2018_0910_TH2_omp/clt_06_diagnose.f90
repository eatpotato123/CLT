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



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!      !!!!!!!    !!!!!!!!!    !!!!!!!!!!!!!!!        !!!!!!!!
!!!!!  !!!!  !!!!!!  !!!!!!!!!  !!  !!!!!!!!!!!!   !!!!!!!!!!!!!!!
!!!!!  !!!!!  !!!!!  !!!!!!!!  !!!!  !!!!!!!!!!  !!!!!!!!!!!!!!!!!
!!!!!  !!!!!!  !!!!  !!!!!!!  !!!!!!  !!!!!!!!  !!!!!!!!    !!!!!!
!!!!!  !!!!!  !!!!!  !!!!!!            !!!!!!!  !!!!!!!!!!  !!!!!!
!!!!!  !!!!  !!!!!!  !!!!!  !!!!!!!!!!  !!!!!!!  !!!!!!!!  !!!!!!!
!!!!!      !!!!!!!    !!!  !!!!!!!!!!!!  !!!!!!!!         !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module for diagnose subroutines.
! written by S. WANG
! sortted by H.W. Zhang during 2018 AUG.
! contact changhw@zju.edu.cn if there is any problem.

! including subroutines as follows
! diagn
! diagnatxmode
! energy
! integ
! funmax
! funm
! init_dgn
! diagn_nmmode(fxz,jdg)
! diagn_max
! diagn_maxmin
! diagn_brmax0


!ws**********************************************************************
      subroutine diagn
      use declare
      include 'mpif.h'
      real*8 cxmax,x_cxmax,z_cxmax,cymax,x_cymax,z_cymax,czmax,x_czmax,z_czmax,crmax,x_crmax,z_crmax,cpmax,x_cpmax,z_cpmax
      real*8 cxmin,x_cxmin,z_cxmin,cymin,x_cymin,z_cymin,czmin,x_czmin,z_czmin,crmin,x_crmin,z_crmin,cpmin,x_cpmin,z_cpmin
      real*8 bxmax,x_bxmax,z_bxmax,bymax,x_bymax,z_bymax,bzmax,x_bzmax,z_bzmax,brmax,x_brmax,z_brmax,bpmax,x_bpmax,z_bpmax
      real*8 bxmin,x_bxmin,z_bxmin,bymin,x_bymin,z_bymin,bzmin,x_bzmin,z_bzmin,brmin,x_brmin,z_brmin,bpmin,x_bpmin,z_bpmin

      real*8 cxmax1,x_cxmax1,z_cxmax1,cymax1,x_cymax1,z_cymax1,czmax1,x_czmax1,z_czmax1,crmax1,x_crmax1,z_crmax1,cpmax1, &
                x_cpmax1,z_cpmax1
      real*8 cxmin1,x_cxmin1,z_cxmin1,cymin1,x_cymin1,z_cymin1,czmin1,x_czmin1,z_czmin1,crmin1,x_crmin1,z_crmin1,cpmin1, &
                x_cpmin1,z_cpmin1
      real*8 bxmax1,x_bxmax1,z_bxmax1,bymax1,x_bymax1,z_bymax1,bzmax1,x_bzmax1,z_bzmax1,brmax1,x_brmax1,z_brmax1,bpmax1, &
                x_bpmax1,z_bpmax1
      real*8 bxmin1,x_bxmin1,z_bxmin1,bymin1,x_bymin1,z_bymin1,bzmin1,x_bzmin1,z_bzmin1,brmin1,x_brmin1,z_brmin1,bpmin1, &
                x_bpmin1,z_bpmin1

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
      call mpi_allreduce(cymax,cymax1,1,mpi_double_precision,mpi_max, &
                       mpi_comm_world,ierror)
      call mpi_allreduce(x_cymax,x_cymax1,1,mpi_double_precision,mpi_min, &
                       mpi_comm_world,ierror)
      call mpi_allreduce(z_cymax,z_cymax1,1,mpi_double_precision,mpi_min, &
                       mpi_comm_world,ierror)

      call mpi_allreduce(cxmax,cxmax1,1,mpi_double_precision,mpi_max, &
                       mpi_comm_world,ierror)
      call mpi_allreduce(x_cxmax,x_cxmax1,1,mpi_double_precision,mpi_min, &
                       mpi_comm_world,ierror)
      call mpi_allreduce(z_cxmax,z_cxmax1,1,mpi_double_precision,mpi_min, &
                       mpi_comm_world,ierror)

      call mpi_allreduce(czmax,czmax1,1,mpi_double_precision,mpi_max, &
                       mpi_comm_world,ierror)
      call mpi_allreduce(x_czmax,x_czmax1,1,mpi_double_precision,mpi_min, &
                       mpi_comm_world,ierror)
      call mpi_allreduce(z_czmax,z_czmax1,1,mpi_double_precision,mpi_min, &
                       mpi_comm_world,ierror)

      call mpi_allreduce(crmax,crmax1,1,mpi_double_precision,mpi_max, &
                       mpi_comm_world,ierror)
      call mpi_allreduce(x_crmax,x_crmax1,1,mpi_double_precision,mpi_min, &
                       mpi_comm_world,ierror)
      call mpi_allreduce(z_crmax,z_crmax1,1,mpi_double_precision,mpi_min, &
                       mpi_comm_world,ierror)

      call mpi_allreduce(cpmax,cpmax1,1,mpi_double_precision,mpi_max, &
                       mpi_comm_world,ierror)
      call mpi_allreduce(x_cpmax,x_cpmax1,1,mpi_double_precision,mpi_min, &
                       mpi_comm_world,ierror)
      call mpi_allreduce(z_cpmax,z_cpmax1,1,mpi_double_precision,mpi_min, &
                       mpi_comm_world,ierror)
!bfeild   -----------------------------------------------------------------
      call mpi_allreduce(bymax,bymax1,1,mpi_double_precision,mpi_max, &
                       mpi_comm_world,ierror)
      call mpi_allreduce(x_bymax,x_bymax1,1,mpi_double_precision,mpi_min, &
                       mpi_comm_world,ierror)
      call mpi_allreduce(z_bymax,z_bymax1,1,mpi_double_precision,mpi_min, &
                       mpi_comm_world,ierror)

      call mpi_allreduce(bxmax,bxmax1,1,mpi_double_precision,mpi_max, &
                       mpi_comm_world,ierror)
      call mpi_allreduce(x_bxmax,x_bxmax1,1,mpi_double_precision,mpi_min, &
                       mpi_comm_world,ierror)
      call mpi_allreduce(z_bxmax,z_bxmax1,1,mpi_double_precision,mpi_min, &
                       mpi_comm_world,ierror)

      call mpi_allreduce(bzmax,bzmax1,1,mpi_double_precision,mpi_max, &
                       mpi_comm_world,ierror)
      call mpi_allreduce(x_bzmax,x_bzmax1,1,mpi_double_precision,mpi_min, &
                       mpi_comm_world,ierror)
      call mpi_allreduce(z_bzmax,z_bzmax1,1,mpi_double_precision,mpi_min, &
                       mpi_comm_world,ierror)

      call mpi_allreduce(brmax,brmax1,1,mpi_double_precision,mpi_max, &
                       mpi_comm_world,ierror)
      call mpi_allreduce(x_brmax,x_brmax1,1,mpi_double_precision,mpi_min, &
                       mpi_comm_world,ierror)
      call mpi_allreduce(z_brmax,z_brmax1,1,mpi_double_precision,mpi_min, &
                       mpi_comm_world,ierror)

      call mpi_allreduce(bpmax,bpmax1,1,mpi_double_precision,mpi_max, &
                       mpi_comm_world,ierror)
      call mpi_allreduce(x_bpmax,x_bpmax1,1,mpi_double_precision,mpi_min, &
                       mpi_comm_world,ierror)
      call mpi_allreduce(z_bpmax,z_bpmax1,1,mpi_double_precision,mpi_min, &
                       mpi_comm_world,ierror) 

      if(nrank.eq.0) then
!      open(unit=13,file='current.dat',status='unknown',form='formatted',position='append')
      write(13,1000)time,cymax1,x_cymax1,z_cymax1,cxmax1,x_cxmax1,z_cxmax1,czmax1,x_czmax1,z_czmax1,crmax1, &
                x_crmax1,z_crmax1,cpmax1,x_cpmax1,z_cpmax1
!      open(unit=15,file='bfield.dat',status='unknown',form='formatted',position='append')
      write(15,1000)time,bymax1,x_bymax1,z_bymax1,bxmax1,x_bxmax1,z_bxmax1,bzmax1,x_bzmax1,z_bzmax1,brmax1, &
                x_brmax1,z_brmax1,bpmax1,x_bpmax1,z_bpmax1
1000  format(16(1x,e12.5))

!      print*,cxmax1,cymax1,czmax1,cxmin1,cymin1,czmin1,time
      endif


      call diagnatxmode

      return
      end

!ws*************************************************      
      subroutine diagnatxmode
      use declare
      real*8 vsmode
      include 'mpif.h' 
    
      if(nrank==nrank_mode) then
      vsmode=x(jxmode,jzmode,1,3)*x(jxmode,jzmode,1,8)-x(jxmode,jzmode,1,5)*x(jxmode,jzmode,1,6)
      write(16,1300) xx(jxmode),time,vsmode,(x1(jxmode,jzmode,1,i),i=1,8),(cur(jxmode,jzmode,1,ic),ic=1,3), &
                (ef(jxmode,jzmode,1,ie),ie=1,3)
1300  format(17(1x,e12.5))
      endif
      return
      end

!ws*******************************************************	
	subroutine energy ! deleted now for nofftw in sunway(0810)
	end

!ws*******************************************************	
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

!ws*******************************************************	
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

!ws*******************************************************
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

!ws*******************************************************
	subroutine init_dgn
     use declare
     include 'mpif.h' 
     integer itmp,ii,ikdgn
     character*10 output
     character*3 cn1
      
!      qdgn(jdgn)=jdgn

      j=1 
      do jdgn=1,mdgn_rs         
      do while(q_nova(j) .ge. qdgn(jdgn))
      j=j+1
      enddo
      if(j==1) j=j+2
      if(j==2) j=j+1
      call interp1d3l(psival_nova(j-2),psival_nova(j-1),psival_nova(j),psival_nova(j+1), &
                    q_nova(j-2),q_nova(j-1),q_nova(j),q_nova(j+1),qdgn(jdgn),psdgn(jdgn))
       do i=1,n2th+5
         call interp1d3l(xxst(i,j-2),xxst(i,j-1),xxst(i,j),xxst(i,j+1), &
                    q_nova(j-2),q_nova(j-1),q_nova(j),q_nova(j+1), qdgn(jdgn),xxdgn(i,jdgn))
         call interp1d3l(zzst(i,j-2),zzst(i,j-1),zzst(i,j),zzst(i,j+1), &
                    q_nova(j-2),q_nova(j-1),q_nova(j),q_nova(j+1), qdgn(jdgn),zzdgn(i,jdgn))
       enddo
      enddo
      
      do jdgn=mdgn_rs+1,mdgn
      do while(q_nova(j) .lt. qdgn(jdgn))
      j=j+1
      enddo
      if(j==1) j=j+2
      if(j==2) j=j+1
      call interp1d3l(psival_nova(j-2),psival_nova(j-1),psival_nova(j),psival_nova(j+1), &
                    q_nova(j-2),q_nova(j-1),q_nova(j),q_nova(j+1), qdgn(jdgn),psdgn(jdgn))
       do i=1,n2th+5
         call interp1d3l(xxst(i,j-2),xxst(i,j-1),xxst(i,j),xxst(i,j+1), &
                    q_nova(j-2),q_nova(j-1),q_nova(j),q_nova(j+1), qdgn(jdgn),xxdgn(i,jdgn))
         call interp1d3l(zzst(i,j-2),zzst(i,j-1),zzst(i,j),zzst(i,j+1), &
                    q_nova(j-2),q_nova(j-1),q_nova(j),q_nova(j+1), qdgn(jdgn),zzdgn(i,jdgn))
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

!ws*******************************************************
	subroutine diagn_nmmode ! deleted now for nofftw in sunway(0810)
	end

!ws**********************************************************************
      subroutine diagn_max
      use declare
      include 'mpif.h'
      real*8 vsmax,eymax,rhomax,pmax,cymax,crmax,cpmax,bymax,brmax,bpmax,vymax,vrmax,vpmax
      real*8 vsmax1,eymax1,rhomax1,pmax1,cymax1,crmax1,cpmax1,bymax1,brmax1,bpmax1,vymax1,vrmax1,vpmax1
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
       eymax=maxval(ef(:,:,:,2))
      

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
      call mpi_allreduce(cymax,cymax1,1,mpi_double_precision,mpi_max, &
                       mpi_comm_world,ierror)

      call mpi_allreduce(crmax,crmax1,1,mpi_double_precision,mpi_max, &
                       mpi_comm_world,ierror)

      call mpi_allreduce(cpmax,cpmax1,1,mpi_double_precision,mpi_max, &
                       mpi_comm_world,ierror)

!bfeild   -----------------------------------------------------------------
      call mpi_allreduce(bymax,bymax1,1,mpi_double_precision,mpi_max, &
                       mpi_comm_world,ierror)

      call mpi_allreduce(brmax,brmax1,1,mpi_double_precision,mpi_max, &
                       mpi_comm_world,ierror)

      call mpi_allreduce(bpmax,bpmax1,1,mpi_double_precision,mpi_max, &
                       mpi_comm_world,ierror)

!v   -----------------------------------------------------------------
      call mpi_allreduce(vymax,vymax1,1,mpi_double_precision,mpi_max, &
                       mpi_comm_world,ierror)

      call mpi_allreduce(vrmax,vrmax1,1,mpi_double_precision,mpi_max, &
                       mpi_comm_world,ierror)

      call mpi_allreduce(vpmax,vpmax1,1,mpi_double_precision,mpi_max, &
                       mpi_comm_world,ierror)

      call mpi_allreduce(vsmax,vsmax1,1,mpi_double_precision,mpi_max, &
                       mpi_comm_world,ierror)

      call mpi_allreduce(eymax,eymax1,1,mpi_double_precision,mpi_max, &
                       mpi_comm_world,ierror)
      call mpi_allreduce(rhomax,rhomax1,1,mpi_double_precision,mpi_max, &
                       mpi_comm_world,ierror)
      call mpi_allreduce(pmax,pmax1,1,mpi_double_precision,mpi_max, &
                       mpi_comm_world,ierror)


      if(nrank.eq.0) then
      write(11,1000)time,vsmax1,eymax1,rhomax1,pmax1,cymax1,crmax1,cpmax1,bymax1,brmax1,bpmax1,vymax1,vrmax1,vpmax1
1000  format(14(1x,e12.5))

!      print*,cxmax1,cymax1,czmax1,cxmin1,cymin1,czmin1,time
      endif


!      call diagnatxmode

      return
      end

!ws**********************************************************************
      subroutine diagn_maxmin
      use declare
      use declare_oxpoint
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
       wsmax(2) =maxval( ef(:,:,:,2))
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
       wsmin(2) =minval( ef(:,:,:,2))
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
!      call mpi_allreduce(wsmax,wsmax1,13,mpi_double_precision,mpi_max, &
!                       mpi_comm_world,ierror)
!      call mpi_allreduce(wsmin,wsmin1,13,mpi_double_precision,mpi_min, &
!                       mpi_comm_world,ierror)

      call mpi_allreduce(wsmax,wsmax1,13,mpi_double_precision,mpi_max, &
            mpi_comm_world,ierror)
      call mpi_allreduce(wsmin,wsmin1,13,mpi_double_precision,mpi_min, &
            mpi_comm_world,ierror)


      if(nrank.eq.0) then
      write(111,1000)time,wsmax1(:)
      write(112,1000)time,wsmin1(:)

1000  format(14(1x,e12.5))

!      print*,cxmax1,cymax1,czmax1,cxmin1,cymin1,czmin1,time
      endif

      br_max=wsmax1(9)
!      call diagnatxmode

      return
      end

!ws**********************************************************************
      subroutine diagn_brmax0
      use declare
      use declare_oxpoint
      include 'mpif.h'
      real*8  brmax

      do jy=iy_first,iy_last
       br(:,:,jy)= x1(:,:,jy,6)*wx2r(:,:)+ x1(:,:,jy,8)*wz2r(:,:)
      enddo
       brmax =maxval( br(:,:,:))
      call mpi_allreduce(brmax,br_max0,1,mpi_double_precision,mpi_max, &
            mpi_comm_world,ierror)
      return
      end
	
	
	
	
	