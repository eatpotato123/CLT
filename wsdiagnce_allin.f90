MODULE DECLARE
    integer,parameter::mxt=256,mzt=256,myt=16,n=200*2,np=mxt,nnn=1000
    integer,parameter::nprx=4,nprz=8,npry=1,mnst=73
    integer,parameter::mxm=mxt/nprx,mzm=mzt/nprz,mym=myt/npry,mx=mxm+4,mz=mzm+4,my=mym+4
!    real*8,parameter::pi=3.141592654 !,rmin=3.0,rmax=5.0,zmin=-1.5,zmax=1.5
    real*8, parameter :: pi=dacos(-1.d0)
    !**************************************************************************************************************
    !the variables below are used for trace3d
    logical,parameter:: itrace3d=.false.!.true.
    integer,parameter::gap=10
    real*8,parameter::maxlnum=myt*nnn*n/gap/4
    real*8 lcr(maxlnum),lcfy(maxlnum),lcz(maxlnum)
    integer lnum
    !********************************************************************************************************************
    !the variables below are used for q_compute
    integer np00(np),nt00(np)
    real*8 q00(np)
    !*******************************************************************************************************************
    real*8 rmin,rmax,zmin,zmax
    
    integer i,j,k,numr,numfy,numz,zmid,times,halfs,nn,loc
    real*8 detr,detfy,detz,dr,dfy,dz,volume,dr2,dz2
    real*8 r0,r1,fy0,fy1,z0,z1,time,r12,fy12,z12
    real*8 br0,bfy0,bz0,br1,bfy1,bz1
    real*8 r(mxt),fy(myt),z(mzt)
    real*8 xx(mxt),zz(mzt),yy(0:myt+1)
    real*8 br(mxt,myt,mzt),bfy(mxt,myt,mzt),bz(mxt,myt,mzt)
    real*8 pointr(np,n),pointz(np,n)    
    real*8 pointr1(np,n,4),pointz1(np,n,4)
    integer nrkx,nrkz,nrky,nrank      
    integer m,jx,jy,jz
    character*12 output
    character*3 cn1,cn
    integer ierr
    integer numprocs,nstrk
    
END MODULE DECLARE

!mxt,mzt,myt分别是R,Z,fy方向的格点数
!n为磁力线所绕的圈数。即在fy=0面上与fy=pi面上的点数。
!np表示在R轴上所选的起始点个数，也即是画的磁面的个数。（如果改格点的话需要重新计算np)
!nn表示运算进行的步数
!nnn表示每步所走的大小是detfy/nnn
!nprx,nprz分别为R,z方向分块的数
!nst表示时间


program    tracefieldline
    use declare
    implicit none
    include 'mpif.h'
    
    integer nst,nst0
    

    
    call mpi_init(ierr)
    call mpi_comm_size(mpi_comm_world,numprocs,ierr)
    call mpi_comm_rank(mpi_comm_world,nstrk,ierr)
    
    !do nst=0,30,1
        nst0=0
        nst=nst0+nstrk
        call getgrid()
        call getparalleldata(nst)
        call trace()
    !    call outputdata(nst)
        call outputdata1(nst)
        if(itrace3d)    then
            call output3d(nst)
        endif
        call outq(nst)
!        if(nstrk==0) call getper
    !enddo
        
    call mpi_finalize(ierr)
end




subroutine getgrid()
    use declare
    implicit none

     
    open(unit=201,file='gridxx.dat',status='unknown',form='formatted')
    read(201,200)(xx(jx),jx=1,mxt)
    close(201)
    open(unit=202,file='gridzz.dat',status='unknown',form='formatted')
    read(202,200)(zz(jz),jz=1,mzt)
    close(202)
    open(unit=203,file='gridyy.dat',status='unknown',form='formatted')
    read(203,200)(yy(jy),jy=1,myt)
    close(203)
      yy(0)=yy(1)-(yy(2)-yy(1))
      yy(myt+1)=yy(myt)+(yy(myt)-yy(myt-1))
     rmin=minval(xx)
     rmax=maxval(xx)
     zmin=minval(zz)
     zmax=maxval(zz)
200  format(1x,e12.5)

end subroutine


!读入数据子程序
subroutine getdata(nst)
    use declare
    implicit none
    integer nst,ncase,nstep
    real*8,dimension(mxt,mzt,myt,8) :: xt
    real*8,dimension(mx,mz,my,8) :: x   
   
    
   !读取数据文件中的数据
    write(cn,'(i3.3)')nst
    output='tk'//cn
    open(unit=7,file=output,status='unknown',form='unformatted')
    read(7)ncase,nstep,time
    read(7)x
    close(7)
    
    !由于输出数据的顺序不一致，转换一下
    do k=1,myt,1
        do j=1,mzt,1
            do i=1,mxt,1
                br(i,k,j)=x(i,j,k,6)
                bfy(i,k,j)=x(i,j,k,7)
                bz(i,k,j)=x(i,j,k,8)
            enddo
        enddo
    enddo
end

subroutine getparalleldata(nst)
    use declare
    implicit none
    include 'mpif.h'
    integer nst,ncase,nstep,ip,ic,ie

    real*8,dimension(mxt,mzt,myt,13) :: x1
    real*8,dimension(mxt,mzt,myt,8) :: xt
    real*8,dimension(mx,mz,my,8) :: x
    real*8,dimension(mxt,mzt,8) :: xtint
    real*8,dimension(mx,mz,8) :: xint
    real*8,dimension(mx,mz,3) :: cint
    real*8,dimension(mxt,mzt,myt,3) :: cur,Ef
    real*8,dimension(mxt,mzt) :: tsxz,dbt,d2dbt,etat,psi
    real*8,dimension(mxt,mzt) :: pmut,fmut
    real*8,dimension(mx,mz,1) :: dvb,d2dvb
    real*8,dimension(mx,mz,my) :: eta
    real*8,dimension(mx,mz) :: fmu,pmu
    character*9 output1
    character*3 cny,cn11
    integer, dimension(myt) :: jym,jyp
    integer,dimension(numprocs) :: nnst
    real*8,dimension(numprocs) :: tnst
    real*8 bzdy,bydz,bxdz,bzdx,Rbydx,bxdy,psia1,pst_dx,pst_dz,qsf,p,pp,g,gp,kap_pal,kap_perp
    real*8 fm1,fm2,fm3,f0,fp1,fp2,fp3,a,b,c,d,xm1,x0,xp1,d1f2,d1xf2
!  d1f2= d f / dx  with second-order accuracy central difference
      d1f2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(fp1-f0) &
         -(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)
!  d1xf2= d Rf / dx  with second-order accuracy central difference
      d1xf2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(xp1*fp1-x0*f0) &
         -(xp1-x0)/(xm1-x0)*(xm1*fm1-x0*f0))/(xm1-xp1)

     write(cn,'(i3.3)')nst 
     write(cn11,'(i3.3)')nst+1 


      do jy=1,myt
      jyp(jy)=jy+1
      jym(jy)=jy-1
      enddo
      jyp(myt)=1
      jym(1)=myt

      open(unit=204,file='gridts.dat',status='unknown',form='formatted')
      read(204,200)((tsxz(jx,jz),jx=1,mxt),jz=1,mzt)
      close(204)
 200  format(1x,e12.5)

      open(unit=205,file='pst_qpg.dat',status='unknown',form='formatted')
      read(205,100)(((psi(jx,jz),pst_dx,pst_dz,qsf,p,pp,g,gp),jx=1,mxt),jz=1,mzt)
 100  format(8(1x,e12.5))
      close(205)


     do nrky=0,npry-1   
     do nrkx=0,nprx-1
        do nrkz=0,nprz-1
            nrank=nrkz*nprx+nrkx+nrky*nprz*nprx
            
            write(cn1,'(i3.3)')nrank
            write(cn,'(i3.3)')nst

!      output='dvb'//cn1//cn11
!      open(unit=91,file=output,status='unknown',form='unformatted')
!      read(91) dvb(:,:,1)
!      close(91)
!      output='d2dvb'//cn1//cn11
!      open(unit=92,file=output,status='unknown',form='unformatted')
!      read(92) d2dvb(:,:,1)
!      close(92)
!         do jz=1,mzm
!           do jx=1,mxm
!             dbt(nrkx*mxm+jx,nrkz*mzm+jz)=dvb(jx+2,jz+2,1)
!             d2dbt(nrkx*mxm+jx,nrkz*mzm+jz)=d2dvb(jx+2,jz+2,1)
!           end do
!         end do
      output='dssp'//cn1   
      open(unit=98,file=output,status='unknown',form='unformatted')
      read(98)ncase,nstep,time
      read(98)eta,fmu,pmu !,kap_pal,kap_perp
      close(98)

      do jz=1,mzm
      do jx=1,mxm
         etat(nrkx*mxm+jx,nrkz*mzm+jz)=eta(jx+2,jz+2,3)
!         pmut(nrkx*mxm+jx,nrkz*mzm+jz)=pmu(jx+2,jz+2)
!         fmut(nrkx*mxm+jx,nrkz*mzm+jz)=fmu(jx+2,jz+2)
      enddo
      enddo

!	output='dssp_tot'
!	open(unit=22,file=output,status='unknown',form='formatted')
!	write(22,2201)((etat(jx,jz),pmut(jx,jz),fmut(jx,jz),jx=1,mxt),jz=1,mzt)
! 2201 format(3(1x,e12.4E4))
! 	close(22)


        output='int'//cn1   
          open(unit=99,file=output,status='unknown',form='unformatted')
          read(99)ncase,nstep,time
          read(99)xint,cint
          close(99)
        do m=1,8
            do jz=1,mzm
               do jx=1,mxm
                 xtint(nrkx*mxm+jx,nrkz*mzm+jz,m)=xint(jx+2,jz+2,m)
               enddo
            enddo
        enddo
    !

            
           !output='tk'//cn1(nrank)//cn(nst)
           output='tk'//cn1//cn
            open(unit=7,file=output,status='unknown',form='unformatted')
            read(7)ncase,nstep,time
            read(7)x
            close(7)

            do m=1,8
                do jy=1,mym
                    do jz=1,mzm
                        do jx=1,mxm
                           ! xt(nrkx*mxm+jx,nrkz*mzm+jz,jy,m)=x(jx+2,jz+2,jy+2,m)
                             xt(nrkx*mxm+jx,nrkz*mzm+jz,nrky*mym+jy,m)=x(jx+2,jz+2,jy+2,m)
                        enddo
                    enddo
                enddo
            enddo
        enddo
     enddo
     enddo



    !由于输出数据的顺序不一致，转换一下
    do k=1,myt,1
        do j=1,mzt,1
            do i=1,mxt,1
                br(i,k,j)=xt(i,j,k,6)
                bfy(i,k,j)=xt(i,j,k,7)
                bz(i,k,j)=xt(i,j,k,8)
            enddo
        enddo
    enddo

      do jy=1,myt
      do jz=1,mzt
      do jx=1,mxt
      do m=1,8
      x1(jx,jz,jy,m)=xt(jx,jz,jy,m)-xtint(jx,jz,m)
      enddo
      x1(jx,jz,jy,9)=x1(jx,jz,jy,6)*dcos(tsxz(jx,jz))+x1(jx,jz,jy,8)*dsin(tsxz(jx,jz))
      x1(jx,jz,jy,10)=-x1(jx,jz,jy,6)*dsin(tsxz(jx,jz))+x1(jx,jz,jy,8)*dcos(tsxz(jx,jz))          
      x1(jx,jz,jy,11)=x1(jx,jz,jy,3)*dcos(tsxz(jx,jz))+x1(jx,jz,jy,5)*dsin(tsxz(jx,jz))
      x1(jx,jz,jy,12)=-x1(jx,jz,jy,3)*dsin(tsxz(jx,jz))+x1(jx,jz,jy,5)*dcos(tsxz(jx,jz))
      x1(jx,jz,jy,13)=x1(jx,jz,jy,3)*xt(jx,jz,jy,8)-x1(jx,jz,jy,5)*xt(jx,jz,jy,6)
      enddo
      enddo
      enddo

      psia1=-1.e-4
      psia1=1.
      psia1=100.
      do jz=2,mzt-1
      do jx=2,mxt-1
      if(psi(jx,jz).lt.psia1) then

      do jy=1,myt
      bzdy=d1f2(x1(jx,jz,jym(jy),8),x1(jx,jz,jy,8),x1(jx,jz,jyp(jy),8),yy(jy-1),yy(jy),yy(jy+1)) 
      bydz=d1f2(x1(jx,jz-1,jy,7),x1(jx,jz,jy,7),x1(jx,jz+1,jy,7),zz(jz-1),zz(jz),zz(jz+1))
      bxdz=d1f2(x1(jx,jz-1,jy,6),x1(jx,jz,jy,6),x1(jx,jz+1,jy,6),zz(jz-1),zz(jz),zz(jz+1))
      bzdx=d1f2(x1(jx-1,jz,jy,8),x1(jx,jz,jy,8),x1(jx+1,jz,jy,8),xx(jx-1),xx(jx),xx(jx+1))
      bxdy=d1f2(x1(jx,jz,jym(jy),6),x1(jx,jz,jy,6),x1(jx,jz,jyp(jy),6),yy(jy-1),yy(jy),yy(jy+1))
      Rbydx=d1xf2(x1(jx-1,jz,jy,7),x1(jx,jz,jy,7),x1(jx+1,jz,jy,7),xx(jx-1),xx(jx),xx(jx+1))      
      cur(jx,jz,jy,1)=bzdy/xx(jx)-bydz
      cur(jx,jz,jy,2)=bxdz-bzdx
      cur(jx,jz,jy,3)=(Rbydx-bxdy)/xx(jx)
      enddo

      endif
      enddo
      enddo

      do jy=1,myt
      do jz=1,mzt
      do jx=1,mxt
      Ef(jx,jz,jy,1)=-x1(jx,jz,jy,4)*xt(jx,jz,jy,8)-xtint(jx,jz,4)*x1(jx,jz,jy,8)+x1(jx,jz,jy,5)*xt(jx,jz,jy,7)+xtint(jx,jz,5)*x1(jx,jz,jy,7)+etat(jx,jz)*cur(jx,jz,jy,1)
      Ef(jx,jz,jy,2)=-x1(jx,jz,jy,5)*xt(jx,jz,jy,6)-xtint(jx,jz,5)*x1(jx,jz,jy,6)+x1(jx,jz,jy,3)*xt(jx,jz,jy,8)+xtint(jx,jz,3)*x1(jx,jz,jy,8)+etat(jx,jz)*cur(jx,jz,jy,2)
      Ef(jx,jz,jy,3)=-x1(jx,jz,jy,3)*xt(jx,jz,jy,7)-xtint(jx,jz,3)*x1(jx,jz,jy,7)+x1(jx,jz,jy,4)*xt(jx,jz,jy,6)+xtint(jx,jz,4)*x1(jx,jz,jy,6)+etat(jx,jz)*cur(jx,jz,jy,3)
      enddo
      enddo
      enddo 


!      if(mod(nst,10).eq.0) then
!      output='x3d'//cn
!      open(unit=9,file=output,status='unknown',form='formatted')
!      write(9,51)((((xt(jx,jz,jy,i),i=1,8),jx=1,mxt),jz=1,mzt),jy=1,my)
! 51   format(8(1x,e12.5))
!      close(9)
!      output='x13d'//cn
!      open(unit=96,file=output,status='unknown',form='formatted')
!      write(96,901)((((x1(jx,jz,jy,i),i=1,8),x1(jx,jz,jy,13),jx=1,mxt),jz=1,mzt),jy=1,my)
!901   format(9(1x,e12.5))
!      close(96)
!      endif


	output='x3d'//cn
	open(unit=9,file=output,status='unknown',form='binary')
	write(9) ((((xt(jx,jz,jy,i),i=1,8),jx=1,mxt),jz=1,mzt),jy=1,myt)
	close(9)

	output='x13d'//cn
	open(unit=96,file=output,status='unknown',form='binary')
	write(96) ((((x1(jx,jz,jy,i),i=1,8),jx=1,mxt),jz=1,mzt),jy=1,myt)
	close(96)

      if(mod(nst,10).eq.0) then
      k=0
      do jy=1,my,my/4
      write(cny,'(i3.3)') k
      output='x12dk'//cn//cny
      open(unit=k+11,file=output,status='unknown',form='formatted')
      write(k+11,900) (((x1(jx,jz,jy,i),i=1,8),x1(jx,jz,jy,13),jx=1,mxt),jz=1,mzt)
      close(k+11)
      k=k+1
      enddo
900   format(9(1x,e12.4E4))
      endif

!      output='db'//cn11
!      open(unit=93,file=output,status='unknown',form='formatted')
!      write(93,300)(((dbt(jx,jz),d2dbt(jx,jz)),jx=1,mxt),jz=1,mzt)
!300   format(2(1x,e12.5))
!      close(93)

      output='x12d'//cn
      open(unit=19,file=output,status='unknown',form='formatted')
      write(19,1300)(((x1(jx,jz,1,i),i=1,13),jx=1,mxt),jz=1,mzt)
1300  format(13(1x,e12.4E4))
      close(19)
      
      output='pt'//cn
      open(unit=44,file=output,status='unknown',form='formatted')
      write(44,1800)((xt(jx,jz,1,2),jx=1,mxt),jz=1,mzt)
1800  format(1x,e12.4E4)
      close(44)
      
      output='xt'//cn
      open(unit=45,file=output,status='unknown',form='formatted')
      write(45,1801)(((xt(jx,jz,1,i),i=1,8),jx=1,mxt),jz=1,mzt)
1801  format(8(1x,e12.4E4))
      close(45)

      output='x12dyz'//cn
      jx=mxt/2
      open(unit=9,file=output,status='unknown',form='formatted')
      write(9,50)(((x1(jx,jz,jy,i),i=1,8),jz=1,mzt),jy=1,myt)
 50   format(8(1x,e12.4E4))
      close(9)
      output='x12dyx'//cn
      jz=mzt/2
      open(unit=91,file=output,status='unknown',form='formatted')
      write(91,50)(((x1(jx,jz,jy,i),i=1,8),jx=1,mxt),jy=1,myt)
      close(91)

      output='xce12d'//cn
      open(unit=81,file=output,status='unknown',form='formatted')
      write(81,600)((((cur(jx,jz,1,ic),ic=1,3),(Ef(jx,jz,1,ie),ie=1,3)),jx=1,mxt),jz=1,mzt)
 600  format(6(1x,e12.4E4))
      close(81)
      output='xce12dy'//cn
      jx=mxt/2
      open(unit=82,file=output,status='unknown',form='formatted')
      write(82,600)((((cur(jx,jz,jy,ic),ic=1,3),(Ef(jx,jz,jy,ie),ie=1,3)),jz=1,mzt),jy=1,myt)
      close(82)


      CALL MPI_Gather(nst,1,MPI_INTEGER,nnst,1, &
            MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_Gather(time,1,MPI_DOUBLE_PRECISION,tnst,1, &
            MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      if(nstrk==0) then
      output='xtint'
      open(unit=29,file=output,status='unknown',form='formatted')
      write(29,50)(((xtint(jx,jz,i),i=1,8),jx=1,mxt),jz=1,mzt)
      close(29)
      open(unit=39,file='nstt.dat',status='unknown',form='formatted')
      write(39,"(1x,i5,1x,e12.5)") (nnst(ip),tnst(ip),ip=1,numprocs)
      close(39)
      endif  
    
end










!读入并行程序输出的数据
subroutine getparalleldata_y(nst)
    use declare
    implicit none
    
    integer nst,ncase,nstep
    real*8,dimension(mxt,mzt,myt,8) :: xt
    real*8,dimension(mx,mz,my,8) :: x  

   do nrky=0,npry-1   
     do nrkx=0,nprx-1
        do nrkz=0,nprz-1
            nrank=nrkz*nprx+nrkx+nrky*nprz*nprx
            
            write(cn1,'(i3.3)')nrank
            write(cn,'(i3.3)')nst 
            
           !output='tk'//cn1(nrank)//cn(nst)
           output='tk'//cn1//cn
            open(unit=7,file=output,status='unknown',form='unformatted')
            read(7)ncase,nstep,time
            read(7)x
            close(7)

            do m=1,8
                do jy=1,mym
                    do jz=1,mzm
                        do jx=1,mxm
                            xt(nrkx*mxm+jx,nrkz*mzm+jz,nrky*mym+jy,m)=x(jx+2,jz+2,jy+2,m)
                        enddo
                    enddo
                enddo
            enddo
        enddo

    enddo
  enddo

    !由于输出数据的顺序不一致，转换一下
    do k=1,myt,1
        do j=1,mzt,1
            do i=1,mxt,1
                br(i,k,j)=xt(i,j,k,6)
                bfy(i,k,j)=xt(i,j,k,7)
                bz(i,k,j)=xt(i,j,k,8)
            enddo
        enddo
    enddo
    
end



















    
!输出fy=0面的磁面结构
subroutine outputdata(nst)
    use declare
    implicit none
    
    integer nst
    character*15 outputdatar,outputdataz
    
    write(cn,'(i3.3)')nst
    outputdatar='outputdatar'//cn 
    open(unit=20,file=outputdatar,status='replace',form='formatted')
    
    do j=1,n*2,1
        do  i=1,np,1
            write(20,*) pointr(i,j)
       enddo
    enddo
    close(20)
   
    outputdataz='outputdataz'//cn 
    open(unit=30,file=outputdataz,status='replace',form='formatted')
    
    
    do j=1,n*2,1
       do  i=1,np,1
            write(30,*) pointz(i,j)
        enddo
    enddo
    
    close(30)
end







!输出fy=pi面的磁面结构
subroutine outputdata1(nst)
    use declare
    implicit none
    
    integer nst
    character*15 outputdatar1,outputdataz1
    
    write(cn,'(i3.3)')nst
    outputdatar1='outputdatar1'//cn 
    open(unit=25,file=outputdatar1,status='replace',form='formatted')
    
    
    do j=1,n,1
        do  i=1,np,1
            write(25,400) (pointr1(i,j,k),k=1,4)
       enddo
    enddo
    close(25)
    
    outputdataz1='outputdataz1'//cn 
    open(unit=35,file=outputdataz1,status='replace',form='formatted')
    
    
    
    do j=1,n,1
       do  i=1,np,1
            write(35,400) (pointz1(i,j,k),k=1,4)
        enddo
    enddo     
    close(35)
400  format(4(1x,e14.7))
end
    


subroutine trace()
    use declare
    implicit none
    integer i0,j0
    integer numr0,numz0
    
	!detr=(rmax-rmin)/(mxt-1) !计算网格的小大
	!detfy=2*pi/myt
	!detz=(zmax-zmin)/(mzt-1)
	
	!计算每个网格的坐标
	!do i=1,mxt,1
	    !r(i)=rmin+(i-1)*detr
	!enddo
	
	!do j=1,myt,1
	 !   fy(j)=(j-1)*detfy
	!enddo
	
	!do k=1,mzt,1
	 !   z(k)=zmin+(k-1)*detz
	!enddo
	
	!固定dfy
	!dfy=detfy/nnn
	
	!计算体积元
	!volume=detr*detfy*detz
	
	!zmid=mzt/2

    detfy=2*pi/myt
	!对网格点赋值
    r=xx
    z=zz
    do j=1,myt,1
	    fy(j)=(j-1)*detfy
	enddo
	
    dfy=detfy/nnn
    zmid=mzt/2
!    zmid=128

	 loc=0
!	    do i=5,mxt-5,3
	    do i=1,mxt,1
		    !选择初始点
		    r0=r(i)
		    fy0=0
		    z0=z(zmid+loc)

		    !初始点的磁场
		    br0=br(i,1,zmid+loc)
		    bfy0=bfy(i,1,zmid+loc)
		    bz0=bz(i,1,zmid+loc)
			
            !used in q_comput
            nt00(i)=0
            np00(i)=0
            !*****************************
            
		    do nn=1,myt*nnn*n,1
		        !计算每个步长r方向和z方向前进的距离
			    dr=br0*r0*dfy/bfy0
			    dz=bz0*r0*dfy/bfy0
    			
			    !计算前进一步后的新位置
			    r1=r0+dr
			    fy1=fy0+dfy
			    if(fy1>=2*pi) fy1=fy1-2*pi   
			    z1=z0+dz
    			
			     !判断此时在numr->numr+1,numfy->numfy+1,numz->numz+1格点之间
			    numfy=nn/nnn+1
                if(nn==1)    then
                    do i0=1,mxt-1,1
                        if(r1>=r(i0).and.r1<r(i0+1)) then
                            numr=i0
                            exit
                        endif
                    enddo
                   

                    do j0=1,mzt-1,1
                        if(z1>=z(j0).and.z1<z(j0+1)) then
                            numz=j0
                            exit
                        endif
                    enddo
                    numr0=numr
                    numz0=numz

                else
            
                   
                    do i0=numr0-1,mxt-1,1
                        if(r1>=r(i0).and.r1<r(i0+1)) then
                            numr=i0
                            exit
                        endif
                    enddo
                    
                        
                   
                    do j0=numz0-1,mzt-1,1
                        if(z1>=z(j0).and.z1<z(j0+1)) then
                            numz=j0
                            exit
                        endif
                    enddo
                   
                    
                    

                endif

              
                
		
			    do while(numfy>=myt+1)
			        numfy=numfy-myt
			    enddo
    			
                
			    !插值法计算此处的磁场
			    if(numfy+1<=myt)    then !如果没有超过2pi
			        !下面考虑由于计算精度问题可能出现的跑到边界的问题
			        !case1
			        if((0<numz.and.numz<mzt).and.(0<numr.and.numr<mxt))   then
                        volume=0.5*detfy*(r(numr+1)-r(numr))*(r(numr+1)+r(numr))*(z(numz+1)-z(numz))
			            br1=0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(Fy(numfy+1)-fy1)*(Z(numz+1)-z1)/volume*br(numr,numfy,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(Fy(numfy+1)-fy1)*(Z(numz+1)-z1)/volume*br(numr+1,numfy,numz) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*br(numr,numfy+1,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*br(numr+1,numfy+1,numz) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(Fy(numfy+1)-fy1)*(z1-Z(numz))/volume*br(numr,numfy,numz+1)+0.5*(R(numr)+r1)*(r1-R(numr))*(Fy(numfy+1)-fy1)*(z1-Z(numz))/volume*br(numr+1,numfy,numz+1) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))*(z1-Z(numz))/volume*br(numr,numfy+1,numz+1)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))*(z1-Z(numz))/volume*br(numr+1,numfy+1,numz+1)

			            bfy1=0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(Fy(numfy+1)-fy1)*(Z(numz+1)-z1)/volume*bfy(numr,numfy,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(Fy(numfy+1)-fy1)*(Z(numz+1)-z1)/volume*bfy(numr+1,numfy,numz) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*bfy(numr,numfy+1,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*bfy(numr+1,numfy+1,numz) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(Fy(numfy+1)-fy1)*(z1-Z(numz))/volume*bfy(numr,numfy,numz+1)+0.5*(R(numr)+r1)*(r1-R(numr))*(Fy(numfy+1)-fy1)*(z1-Z(numz))/volume*bfy(numr+1,numfy,numz+1) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))*(z1-Z(numz))/volume*bfy(numr,numfy+1,numz+1)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))*(z1-Z(numz))/volume*bfy(numr+1,numfy+1,numz+1)

                        bz1=0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(Fy(numfy+1)-fy1)*(Z(numz+1)-z1)/volume*bz(numr,numfy,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(Fy(numfy+1)-fy1)*(Z(numz+1)-z1)/volume*bz(numr+1,numfy,numz) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*bz(numr,numfy+1,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*bz(numr+1,numfy+1,numz) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(Fy(numfy+1)-fy1)*(z1-Z(numz))/volume*bz(numr,numfy,numz+1)+0.5*(R(numr)+r1)*(r1-R(numr))*(Fy(numfy+1)-fy1)*(z1-Z(numz))/volume*bz(numr+1,numfy,numz+1) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))*(z1-Z(numz))/volume*bz(numr,numfy+1,numz+1)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))*(z1-Z(numz))/volume*bz(numr+1,numfy+1,numz+1)
				    !case2       
				    elseif((0<numz.and.numz<mzt).and.numr==mxt)   then
                        volume=detfy*(z(numz+1)-z(numz))
				        br1=(Fy(numfy+1)-fy1)*(Z(numz+1)-z1)/volume*br(numr,numfy,numz) &
				            +(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*br(numr,numfy+1,numz) &
				            +(Fy(numfy+1)-fy1)*(z1-Z(numz))/volume*br(numr,numfy,numz+1) &
				            +(fy1-Fy(numfy))*(z1-Z(numz))/volume*br(numr,numfy+1,numz+1)

			            bfy1=(Fy(numfy+1)-fy1)*(Z(numz+1)-z1)/volume*bfy(numr,numfy,numz) &
				            +(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*bfy(numr,numfy+1,numz) &
				            +(Fy(numfy+1)-fy1)*(z1-Z(numz))/volume*bfy(numr,numfy,numz+1) &
				            +(fy1-Fy(numfy))*(z1-Z(numz))/volume*bfy(numr,numfy+1,numz+1)

			            bz1=(Fy(numfy+1)-fy1)*(Z(numz+1)-z1)/volume*bz(numr,numfy,numz) &
				            +(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*bz(numr,numfy+1,numz) &
				            +(Fy(numfy+1)-fy1)*(z1-Z(numz))/volume*bz(numr,numfy,numz+1) &
				            +(fy1-Fy(numfy))*(z1-Z(numz))/volume*bz(numr,numfy+1,numz+1)
				    
			        !case3    
			        elseif(numz==mzt.and.(0<numr.and.numr<mxt))  then
                        volume=0.5*detfy*(r(numr+1)-r(numr))*(r(numr+1)+r(numr))
			            br1=0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(Fy(numfy+1)-fy1)/volume*br(numr,numfy,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(Fy(numfy+1)-fy1)/volume*br(numr+1,numfy,numz) &
				           +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))/volume*br(numr,numfy+1,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))/volume*br(numr+1,numfy+1,numz) 
    				       

			            bfy1=0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(Fy(numfy+1)-fy1)/volume*bfy(numr,numfy,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(Fy(numfy+1)-fy1)/volume*bfy(numr+1,numfy,numz) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))/volume*bfy(numr,numfy+1,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))/volume*bfy(numr+1,numfy+1,numz) 
    				       

			            bz1=0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(Fy(numfy+1)-fy1)/volume*bz(numr,numfy,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(Fy(numfy+1)-fy1)/volume*bz(numr+1,numfy,numz) &
				           +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))/volume*bz(numr,numfy+1,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))/volume*bz(numr+1,numfy+1,numz)
                    endif
    	
			    elseif(numfy==myt) then   !如果刚好在2pi的两边(numfy==myt)
			        !case1
			        if((0<numz.and.numz<mzt).and.(0<numr.and.numr<mxt))   then
                        volume=0.5*detfy*(r(numr+1)-r(numr))*(r(numr+1)+r(numr))*(z(numz+1)-z(numz))
			            br1=0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(2*pi-fy1)*(Z(numz+1)-z1)/volume*br(numr,numfy,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(2*pi-fy1)*(Z(numz+1)-z1)/volume*br(numr+1,numfy,numz) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*br(numr,1,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*br(numr+1,1,numz) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(2*pi-fy1)*(z1-Z(numz))/volume*br(numr,numfy,numz+1)+0.5*(R(numr)+r1)*(r1-R(numr))*(2*pi-fy1)*(z1-Z(numz))/volume*br(numr+1,numfy,numz+1) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))*(z1-Z(numz))/volume*br(numr,1,numz+1)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))*(z1-Z(numz))/volume*br(numr+1,1,numz+1)

                        bfy1=0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(2*pi-fy1)*(Z(numz+1)-z1)/volume*bfy(numr,numfy,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(2*pi-fy1)*(Z(numz+1)-z1)/volume*bfy(numr+1,numfy,numz) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*bfy(numr,1,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*bfy(numr+1,1,numz) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(2*pi-fy1)*(z1-Z(numz))/volume*bfy(numr,numfy,numz+1)+0.5*(R(numr)+r1)*(r1-R(numr))*(2*pi-fy1)*(z1-Z(numz))/volume*bfy(numr+1,numfy,numz+1) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))*(z1-Z(numz))/volume*bfy(numr,1,numz+1)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))*(z1-Z(numz))/volume*bfy(numr+1,1,numz+1)
                        
                        bz1=0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(2*pi-fy1)*(Z(numz+1)-z1)/volume*bz(numr,numfy,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(2*pi-fy1)*(Z(numz+1)-z1)/volume*bz(numr+1,numfy,numz) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*bz(numr,1,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*bz(numr+1,1,numz) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(2*pi-fy1)*(z1-Z(numz))/volume*bz(numr,numfy,numz+1)+0.5*(R(numr)+r1)*(r1-R(numr))*(2*pi-fy1)*(z1-Z(numz))/volume*bz(numr+1,numfy,numz+1) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))*(z1-Z(numz))/volume*bz(numr,1,numz+1)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))*(z1-Z(numz))/volume*bz(numr+1,1,numz+1)

			            
				     !case2   
			        elseif((0<numz.and.numz<mzt).and.numr==mxt)   then
                            volume=detfy*(z(numz+1)-z(numz))
				            br1=(2*pi-fy1)*(Z(numz+1)-z1)/volume*br(numr,numfy,numz) &
				                +(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*br(numr,1,numz) &
				                +(2*pi-fy1)*(z1-Z(numz))/volume*br(numr,numfy,numz+1) &
				                +(fy1-Fy(numfy))*(z1-Z(numz))/volume*br(numr,1,numz+1)

			                bfy1=(2*pi-fy1)*(Z(numz+1)-z1)/volume*bfy(numr,numfy,numz) &
				                +(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*bfy(numr,1,numz) &
				                +(2*pi-fy1)*(z1-Z(numz))/volume*bfy(numr,numfy,numz+1) &
				                +(fy1-Fy(numfy))*(z1-Z(numz))/volume*bfy(numr,1,numz+1)

			                bz1=(2*pi-fy1)*(Z(numz+1)-z1)/volume*bz(numr,numfy,numz) &
				                +(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*bz(numr,1,numz) &
				                +(2*pi-fy1)*(z1-Z(numz))/volume*bz(numr,numfy,numz+1) &
				                +(fy1-Fy(numfy))*(z1-Z(numz))/volume*bz(numr,1,numz+1)
				    
			        !case3      
				    elseif(numz==mzt.and.(0<numr.and.numr<mxt))  then 
                            volume=0.5*detfy*(r(numr+1)-r(numr))*(r(numr+1)+r(numr))  
			                br1=(R(numr+1)-r1)*(2*pi-fy1)/volume*br(numr,numfy,numz)+(r1-R(numr))*(2*pi-fy1)/volume*br(numr+1,numfy,numz) &
				                +(R(numr+1)-r1)*(fy1-Fy(numfy))/volume*br(numr,1,numz)+(r1-R(numr))*(fy1-Fy(numfy))/volume*br(numr+1,1,numz) 
    				       
			                bfy1=(R(numr+1)-r1)*(2*pi-fy1)/volume*bfy(numr,numfy,numz)+(r1-R(numr))*(2*pi-fy1)/volume*bfy(numr+1,numfy,numz) &
				                +(R(numr+1)-r1)*(fy1-Fy(numfy))/volume*bfy(numr,1,numz)+(r1-R(numr))*(fy1-Fy(numfy))/volume*bfy(numr+1,1,numz) 
    				       
			                bz1=(R(numr+1)-r1)*(2*pi-fy1)/volume*bz(numr,numfy,numz)+(r1-R(numr))*(2*pi-fy1)/volume*bz(numr+1,numfy,numz) &
				                +(R(numr+1)-r1)*(fy1-Fy(numfy))/volume*bz(numr,1,numz)+(r1-R(numr))*(fy1-Fy(numfy))/volume*bz(numr+1,1,numz) 
				          
				    
    				 
                   endif
    				       
			    endif
    			
			    !先按照初始点的磁场走半步
			    r12=r0+0.5*dr
			    fy12=fy0+0.5*dfy
			    if(fy12>=2*pi)  fy12=fy12-2*pi
			    z12=z0+0.5*dz
    			
			    !按照预估的终点估算后半步
			    dr2=br1*r1*0.5*dfy/bfy1
			    dz2=bz1*r1*0.5*dfy/bfy1
    			
			    !从中间步出发走完后半部
			    r1=r12+dr2
			    fy1=fy12+0.5*dfy
			    if(fy1>=2*pi)   fy1=fy1-2*pi
			    z1=z12+dz2
    			
			    !判断此时在numr->numr+1,numfy->numfy+1,numz->numz+1格点之间
			     !判断此时在numr->numr+1,numfy->numfy+1,numz->numz+1格点之间
			    numfy=nn/nnn+1
                
                do i0=numr0-1,mxt-1,1
                    if(r1>=r(i0).and.r1<r(i0+1)) then
                        numr=i0
                        exit
                    endif
                enddo
              
               
                do j0=numz0-1,mzt-1,1
                    if(z1>=z(j0).and.z1<z(j0+1)) then
                        numz=j0
                        exit
                    endif
                enddo
                

                numr0=numr
                numz0=numz
		
			    do while(numfy>=myt+1)
			        numfy=numfy-myt
			    enddo

			    !插值法计算此处的磁场
			    if(numfy+1<=myt)    then !如果没有超过2pi
			        !下面考虑由于计算精度问题可能出现的跑到边界的问题
			        !case1
			        if((0<numz.and.numz<mzt).and.(0<numr.and.numr<mxt))   then
                        volume=0.5*detfy*(r(numr+1)-r(numr))*(r(numr+1)+r(numr))*(z(numz+1)-z(numz))
			            br1=0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(Fy(numfy+1)-fy1)*(Z(numz+1)-z1)/volume*br(numr,numfy,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(Fy(numfy+1)-fy1)*(Z(numz+1)-z1)/volume*br(numr+1,numfy,numz) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*br(numr,numfy+1,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*br(numr+1,numfy+1,numz) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(Fy(numfy+1)-fy1)*(z1-Z(numz))/volume*br(numr,numfy,numz+1)+0.5*(R(numr)+r1)*(r1-R(numr))*(Fy(numfy+1)-fy1)*(z1-Z(numz))/volume*br(numr+1,numfy,numz+1) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))*(z1-Z(numz))/volume*br(numr,numfy+1,numz+1)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))*(z1-Z(numz))/volume*br(numr+1,numfy+1,numz+1)

			            bfy1=0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(Fy(numfy+1)-fy1)*(Z(numz+1)-z1)/volume*bfy(numr,numfy,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(Fy(numfy+1)-fy1)*(Z(numz+1)-z1)/volume*bfy(numr+1,numfy,numz) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*bfy(numr,numfy+1,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*bfy(numr+1,numfy+1,numz) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(Fy(numfy+1)-fy1)*(z1-Z(numz))/volume*bfy(numr,numfy,numz+1)+0.5*(R(numr)+r1)*(r1-R(numr))*(Fy(numfy+1)-fy1)*(z1-Z(numz))/volume*bfy(numr+1,numfy,numz+1) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))*(z1-Z(numz))/volume*bfy(numr,numfy+1,numz+1)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))*(z1-Z(numz))/volume*bfy(numr+1,numfy+1,numz+1)

                        bz1=0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(Fy(numfy+1)-fy1)*(Z(numz+1)-z1)/volume*bz(numr,numfy,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(Fy(numfy+1)-fy1)*(Z(numz+1)-z1)/volume*bz(numr+1,numfy,numz) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*bz(numr,numfy+1,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*bz(numr+1,numfy+1,numz) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(Fy(numfy+1)-fy1)*(z1-Z(numz))/volume*bz(numr,numfy,numz+1)+0.5*(R(numr)+r1)*(r1-R(numr))*(Fy(numfy+1)-fy1)*(z1-Z(numz))/volume*bz(numr+1,numfy,numz+1) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))*(z1-Z(numz))/volume*bz(numr,numfy+1,numz+1)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))*(z1-Z(numz))/volume*bz(numr+1,numfy+1,numz+1)
				    !case2       
				    elseif((0<numz.and.numz<mzt).and.numr==mxt)   then
                        volume=detfy*(z(numz+1)-z(numz))
				        br1=(Fy(numfy+1)-fy1)*(Z(numz+1)-z1)/volume*br(numr,numfy,numz) &
				            +(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*br(numr,numfy+1,numz) &
				            +(Fy(numfy+1)-fy1)*(z1-Z(numz))/volume*br(numr,numfy,numz+1) &
				            +(fy1-Fy(numfy))*(z1-Z(numz))/volume*br(numr,numfy+1,numz+1)

			            bfy1=(Fy(numfy+1)-fy1)*(Z(numz+1)-z1)/volume*bfy(numr,numfy,numz) &
				            +(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*bfy(numr,numfy+1,numz) &
				            +(Fy(numfy+1)-fy1)*(z1-Z(numz))/volume*bfy(numr,numfy,numz+1) &
				            +(fy1-Fy(numfy))*(z1-Z(numz))/volume*bfy(numr,numfy+1,numz+1)

			            bz1=(Fy(numfy+1)-fy1)*(Z(numz+1)-z1)/volume*bz(numr,numfy,numz) &
				            +(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*bz(numr,numfy+1,numz) &
				            +(Fy(numfy+1)-fy1)*(z1-Z(numz))/volume*bz(numr,numfy,numz+1) &
				            +(fy1-Fy(numfy))*(z1-Z(numz))/volume*bz(numr,numfy+1,numz+1)
				    
			        !case3    
			        elseif(numz==mzt.and.(0<numr.and.numr<mxt))  then
                        volume=0.5*detfy*(r(numr+1)-r(numr))*(r(numr+1)+r(numr))
			            br1=0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(Fy(numfy+1)-fy1)/volume*br(numr,numfy,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(Fy(numfy+1)-fy1)/volume*br(numr+1,numfy,numz) &
				           +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))/volume*br(numr,numfy+1,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))/volume*br(numr+1,numfy+1,numz) 
    				       

			            bfy1=0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(Fy(numfy+1)-fy1)/volume*bfy(numr,numfy,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(Fy(numfy+1)-fy1)/volume*bfy(numr+1,numfy,numz) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))/volume*bfy(numr,numfy+1,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))/volume*bfy(numr+1,numfy+1,numz) 
    				       

			            bz1=0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(Fy(numfy+1)-fy1)/volume*bz(numr,numfy,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(Fy(numfy+1)-fy1)/volume*bz(numr+1,numfy,numz) &
				           +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))/volume*bz(numr,numfy+1,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))/volume*bz(numr+1,numfy+1,numz)
                    endif
    	
			    elseif(numfy==myt) then   !如果刚好在2pi的两边(numfy==myt)
			        !case1
			        if((0<numz.and.numz<mzt).and.(0<numr.and.numr<mxt))   then
                        volume=0.5*detfy*(r(numr+1)-r(numr))*(r(numr+1)+r(numr))*(z(numz+1)-z(numz))
			            br1=0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(2*pi-fy1)*(Z(numz+1)-z1)/volume*br(numr,numfy,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(2*pi-fy1)*(Z(numz+1)-z1)/volume*br(numr+1,numfy,numz) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*br(numr,1,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*br(numr+1,1,numz) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(2*pi-fy1)*(z1-Z(numz))/volume*br(numr,numfy,numz+1)+0.5*(R(numr)+r1)*(r1-R(numr))*(2*pi-fy1)*(z1-Z(numz))/volume*br(numr+1,numfy,numz+1) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))*(z1-Z(numz))/volume*br(numr,1,numz+1)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))*(z1-Z(numz))/volume*br(numr+1,1,numz+1)

                        bfy1=0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(2*pi-fy1)*(Z(numz+1)-z1)/volume*bfy(numr,numfy,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(2*pi-fy1)*(Z(numz+1)-z1)/volume*bfy(numr+1,numfy,numz) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*bfy(numr,1,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*bfy(numr+1,1,numz) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(2*pi-fy1)*(z1-Z(numz))/volume*bfy(numr,numfy,numz+1)+0.5*(R(numr)+r1)*(r1-R(numr))*(2*pi-fy1)*(z1-Z(numz))/volume*bfy(numr+1,numfy,numz+1) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))*(z1-Z(numz))/volume*bfy(numr,1,numz+1)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))*(z1-Z(numz))/volume*bfy(numr+1,1,numz+1)
                        
                        bz1=0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(2*pi-fy1)*(Z(numz+1)-z1)/volume*bz(numr,numfy,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(2*pi-fy1)*(Z(numz+1)-z1)/volume*bz(numr+1,numfy,numz) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*bz(numr,1,numz)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*bz(numr+1,1,numz) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(2*pi-fy1)*(z1-Z(numz))/volume*bz(numr,numfy,numz+1)+0.5*(R(numr)+r1)*(r1-R(numr))*(2*pi-fy1)*(z1-Z(numz))/volume*bz(numr+1,numfy,numz+1) &
				            +0.5*(R(numr+1)+r1)*(R(numr+1)-r1)*(fy1-Fy(numfy))*(z1-Z(numz))/volume*bz(numr,1,numz+1)+0.5*(R(numr)+r1)*(r1-R(numr))*(fy1-Fy(numfy))*(z1-Z(numz))/volume*bz(numr+1,1,numz+1)

			            
				     !case2   
			        elseif((0<numz.and.numz<mzt).and.numr==mxt)   then
                            volume=detfy*(z(numz+1)-z(numz))
				            br1=(2*pi-fy1)*(Z(numz+1)-z1)/volume*br(numr,numfy,numz) &
				                +(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*br(numr,1,numz) &
				                +(2*pi-fy1)*(z1-Z(numz))/volume*br(numr,numfy,numz+1) &
				                +(fy1-Fy(numfy))*(z1-Z(numz))/volume*br(numr,1,numz+1)

			                bfy1=(2*pi-fy1)*(Z(numz+1)-z1)/volume*bfy(numr,numfy,numz) &
				                +(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*bfy(numr,1,numz) &
				                +(2*pi-fy1)*(z1-Z(numz))/volume*bfy(numr,numfy,numz+1) &
				                +(fy1-Fy(numfy))*(z1-Z(numz))/volume*bfy(numr,1,numz+1)

			                bz1=(2*pi-fy1)*(Z(numz+1)-z1)/volume*bz(numr,numfy,numz) &
				                +(fy1-Fy(numfy))*(Z(numz+1)-z1)/volume*bz(numr,1,numz) &
				                +(2*pi-fy1)*(z1-Z(numz))/volume*bz(numr,numfy,numz+1) &
				                +(fy1-Fy(numfy))*(z1-Z(numz))/volume*bz(numr,1,numz+1)
				    
			        !case3      
				    elseif(numz==mzt.and.(0<numr.and.numr<mxt))  then 
                            volume=0.5*detfy*(r(numr+1)-r(numr))*(r(numr+1)+r(numr))  
			                br1=(R(numr+1)-r1)*(2*pi-fy1)/volume*br(numr,numfy,numz)+(r1-R(numr))*(2*pi-fy1)/volume*br(numr+1,numfy,numz) &
				                +(R(numr+1)-r1)*(fy1-Fy(numfy))/volume*br(numr,1,numz)+(r1-R(numr))*(fy1-Fy(numfy))/volume*br(numr+1,1,numz) 
    				       
			                bfy1=(R(numr+1)-r1)*(2*pi-fy1)/volume*bfy(numr,numfy,numz)+(r1-R(numr))*(2*pi-fy1)/volume*bfy(numr+1,numfy,numz) &
				                +(R(numr+1)-r1)*(fy1-Fy(numfy))/volume*bfy(numr,1,numz)+(r1-R(numr))*(fy1-Fy(numfy))/volume*bfy(numr+1,1,numz) 
    				       
			                bz1=(R(numr+1)-r1)*(2*pi-fy1)/volume*bz(numr,numfy,numz)+(r1-R(numr))*(2*pi-fy1)/volume*bz(numr+1,numfy,numz) &
				                +(R(numr+1)-r1)*(fy1-Fy(numfy))/volume*bz(numr,1,numz)+(r1-R(numr))*(fy1-Fy(numfy))/volume*bz(numr+1,1,numz) 
				          
				    
    				 
                   endif
    				       
			    endif

    			
			    !判断是否已经绕行一周
			    if(mod(nn,myt*nnn)==1)  then
				    times=nn/(myt*nnn)
				    pointr((i-1)/1+1,loc*n+times+1)=r0
				    pointz((i-1)/1+1,loc*n+times+1)=z0
			    endif

			    !判断是否已经绕行半周
			    if(mod(nn,myt*nnn)==1)  then
			        halfs=nn/(myt*nnn)
			        pointr1((i-1)/1+1,loc*n+halfs+1,1)=r0
			        pointz1((i-1)/1+1,loc*n+halfs+1,1)=z0
			    endif


			    if(mod(nn,myt*nnn)==myt*nnn/8+1)  then
			        halfs=nn/(myt*nnn)
			        pointr1((i-1)/1+1,loc*n+halfs+1,2)=r0
			        pointz1((i-1)/1+1,loc*n+halfs+1,2)=z0
			    endif

			    if(mod(nn,myt*nnn)==myt*nnn/4+1)  then
			        halfs=nn/(myt*nnn)
			        pointr1((i-1)/1+1,loc*n+halfs+1,3)=r0
			        pointz1((i-1)/1+1,loc*n+halfs+1,3)=z0
			    endif

			    if(mod(nn,myt*nnn)==myt*nnn/2+1)  then
			        halfs=nn/(myt*nnn)
			        pointr1((i-1)/1+1,loc*n+halfs+1,4)=r0
			        pointz1((i-1)/1+1,loc*n+halfs+1,4)=z0
			    endif
    			
    			
    			!*****************************************************************************************************************
    		    !this secession is used for 3dtrace
    		if(itrace3d)  then
    		if(i==167)  then
    		    if(nn<=maxlnum) then
                    if(mod(nn,gap)==1)  then
                    lnum=nn/gap+1
                    lcr(lnum)=r0
			        lcfy(lnum)=fy0
			        lcz(lnum)=z0
                    lnum=lnum+1
                    endif
                endif
            endif
            endif
            !*************************************************************************************************
    	    !this secession is used for q_comput
    	    if(z0==0)   then
    	        np00((i-1)/1+1)=np00((i-1)/1+1)+1
    	        nt00((i-1)/1+1)=nn/(myt*nnn)
    	    endif
    	    if(z0*z1<0) then
    	        np00((i-1)/1+1)=np00((i-1)/1+1)+1
    	        nt00((i-1)/1+1)=nn/(myt*nnn)
    	    endif
    	    
    	    
    	    
    	    
    	    !**********************************************************************************************
			    !将本次得到的坐标和磁场作为下一次的初始点
			    r0=r1;
			    fy0=fy1;
			    z0=z1;

			    br0=br1;
			    bfy0=bfy1;
			    bz0=bz1;
    			
		    enddo
        enddo
	
end


subroutine getper
    use declare
    implicit none
    include 'mpif.h'
    integer nstp,ncase,nstep,ip,iw
    real*8,dimension(mxt,mzt,myt,13) :: x1
    real*8,dimension(mxt,mzt,myt,8) :: xt
    real*8,dimension(mx,mz,my,8) :: x
    real*8,dimension(mxt,mzt,8) :: xtint
    real*8,dimension(mx,mz,8) :: xint
    real*8,dimension(mx,mz,3) :: cint
    real*8,dimension(mxt,mzt) :: tsxz,dbt,d2dbt
    real*8,dimension(mx,mz,1) :: dvb,d2dvb
    character*9 output1
    character*3 cny,cn11
    integer,dimension(numprocs) :: nnst
    real*8,dimension(numprocs) :: tnst

     do iw=0,2
            nstp=10**iw  
     write(cn,'(i3.3)')nstp 
     do nrkx=0,nprx-1
        do nrkz=0,nprz-1
            nrank=nrkz*nprx+nrkx
            
            write(cn1,'(i3.3)')nrank

       output='int'//cn1   
          open(unit=99,file=output,status='unknown',form='unformatted')
          read(99)ncase,nstep,time
          read(99)xint,cint
          close(99)
        do m=1,8
            do jz=1,mzm
               do jx=1,mxm
                 xtint(nrkx*mxm+jx,nrkz*mzm+jz,m)=xint(jx+2,jz+2,m)
               end do
            end do
        enddo
    !

            
           !output='tk'//cn1(nrank)//cn(nst)
           output='per'//cn//'_'//cn1
            open(unit=7,file=output,status='unknown',form='unformatted')
            read(7)ncase,nstep,time
            read(7)x
            close(7)

            do m=1,8
                do jy=1,mym
                    do jz=1,mzm
                        do jx=1,mxm
                            xt(nrkx*mxm+jx,nrkz*mzm+jz,jy,m)=x(jx+2,jz+2,jy,m)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo


     do jy=1,myt
      do jz=1,mzt
      do jx=1,mxt
      do m=1,8
      x1(jx,jz,jy,m)=xt(jx,jz,jy,m)-xtint(jx,jz,m)
      enddo
      x1(jx,jz,jy,9)=x1(jx,jz,jy,6)*dcos(tsxz(jx,jz))+x1(jx,jz,jy,8)*dsin(tsxz(jx,jz))
      x1(jx,jz,jy,10)=-x1(jx,jz,jy,6)*dsin(tsxz(jx,jz))+x1(jx,jz,jy,8)*dcos(tsxz(jx,jz))          
      x1(jx,jz,jy,11)=x1(jx,jz,jy,3)*dcos(tsxz(jx,jz))+x1(jx,jz,jy,5)*dsin(tsxz(jx,jz))
      x1(jx,jz,jy,12)=-x1(jx,jz,jy,3)*dsin(tsxz(jx,jz))+x1(jx,jz,jy,5)*dcos(tsxz(jx,jz))
      x1(jx,jz,jy,13)=x1(jx,jz,jy,3)*xt(jx,jz,jy,8)-x1(jx,jz,jy,5)*xt(jx,jz,jy,6)
      enddo
      enddo
      enddo
!      if(mod(nst,10).eq.0) then
!!      output='x3d'//cn(nst)
!!      open(unit=9,file=output,status='unknown',form='formatted')
!!      write(9,50)((((xt(jx,jz,jy,i),i=1,8),jx=1,mxt),jz=1,mzt),jy=1,my)
!! 50   format(8(1x,e12.5))
!!      close(9)
!      output='x13d'//cn
!      open(unit=96,file=output,status='unknown',form='formatted')
!      write(96,900)((((x1(jx,jz,jy,i),i=1,8),x1(jx,jz,jy,13),jx=1,mxt),jz=1,mzt),jy=1,my)
!900   format(9(1x,e12.5))
!      close(96)
!      endif


      k=0
      do jy=1,my,my/4
      write(cny,'(i3.3)') k
      output='xp12dk'//cn//cny
      open(unit=k+11,file=output,status='unknown',form='formatted')
      write(k+11,900) (((x1(jx,jz,jy,i),i=1,8),x1(jx,jz,jy,13),jx=1,mxt),jz=1,mzt)
      close(k+11)
      k=k+1
      enddo
900   format(9(1x,e12.5))

      output='xp12d'//cn
      open(unit=19,file=output,status='unknown',form='formatted')
      write(19,1300)(((x1(jx,jz,1,i),i=1,13),jx=1,mxt),jz=1,mzt)
1300  format(13(1x,e12.5))
      close(19)
      output='xp12dy'//cn
      jx=mxt/2
      open(unit=9,file=output,status='unknown',form='formatted')
      write(9,50)(((x1(jx,jz,jy,i),i=1,8),jz=1,mzt),jy=1,myt)
 50   format(8(1x,e12.5))
      close(9)

      enddo

    
end

!ws*******************************************
      subroutine smthxzyt(fsm,kk)
      USE DECLARE
      real*8,dimension(mxt,mzt,myt) :: fsm,wx2,wz2,wy2
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
      do 21 jy=2,myt-1
      do 21 jz=2,mzt-1
      do 21 jx=2,mzt-1
!      if(psi(jx-1,jz).lt.psia .and. psi(jx+1,jz).lt.psia .and. psi(jx,jz-1).lt.psia .and. psi(jx,jz+1).lt.psia) then
      wx2(jx,jz,jy)=difc(fsm(jx-1,jz,jy),fsm(jx,jz,jy),fsm(jx+1,jz,jy),xx(jx-1),xx(jx),xx(jx+1))
      wz2(jx,jz,jy)=difc(fsm(jx,jz-1,jy),fsm(jx,jz,jy),fsm(jx,jz+1,jy),zz(jz-1),zz(jz),zz(jz+1))
      wy2(jx,jz,jy)=difc(fsm(jx,jz,jy-1),fsm(jx,jz,jy),fsm(jx,jz,jy+1),yy(jy-1),yy(jy),yy(jy+1))
!      endif
   21 continue
   
!      do 15 jz=iz_first+1,iz_last-1
!      do 15 jx=ix_first+1,ix_last-1
!      do 15 jy=2,myt-1
!      wy2(jx,jz,jy)=difc(fsm(jx,jz,jy-1),fsm(jx,jz,jy),fsm(jx,jz,jy+1),yy(jy-1),yy(jy),yy(jy+1))
!
!   15 continue

      do 13 jy=2,myt-1
      do 13 jz=2,mzt-1
      do 13 jx=2,mxt-1
      fsm(jx,jz,jy)=fsm(jx,jz,jy)+1./24*(wx2(jx,jz,jy)+wz2(jx,jz,jy)+wy2(jx,jz,jy))      
!           +(.5*(1.+dtanh(pi/2.-thetati(jx))))/20.*w(jx,jz,jy) &
!           +(.5*(1.+dtanh(pi/2.-thetate(jx))))/20.*w(jx,jz,jy)
   13 continue

   10 continue
      return
      end

subroutine output3d(nst)
    use declare
    implicit none

    integer nst
    character*15 output3dr,output3dfy,output3dz
    
    write(cn,'(i3.3)')nst

!输出磁力线上点的坐标分别存入3个data
    output3dr='output3dr'//cn 
    open(unit=40,file=output3dr,status='replace',form='formatted')
    
    do i=1,maxlnum,1
        write(40,*) lcr(i)
    enddo
    close(40)
    
    output3dfy='output3dfy'//cn 
    open(unit=50,file=output3dfy,status='replace',form='formatted')
    
    do i=1,maxlnum,1
        write(50,*) lcfy(i)
    enddo
    close(50)
    

    output3dz='output3dz'//cn 
    open(unit=60,file=output3dz,status='replace',form='formatted')
    
    do i=1,maxlnum,1
        write(60,*) lcz(i)
    enddo
    close(60)
end

subroutine outq(nst)
    use declare
    implicit none
    
    integer nst
    character*15 outputq

    q00=2.0*nt00/np00
    write(cn,'(i3.3)')nst
    
    
    outputq='outq'//cn 
    open(unit=45,file=outputq,status='replace',form='formatted')
    
    do i=1,np,1
        !write(45,*) q00(i),nt00(i),np00(i)
        write(45,*) q00(i),xx((i-1)*1+1)
    enddo
    close(45)
end





