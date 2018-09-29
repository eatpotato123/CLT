!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! openacc demo for clt struct, only the dt is reduced into cpu during
! calculation. data for boundary calculation is updated in host.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	module declare
	integer, parameter :: mx=1280, mz=1280
	integer jx,jz
	real*8, dimension(mx,mz) :: x, A, B, xdif
	real*8 time
	real*8 dt
	!$acc declare create(x,A,B,xdif) copyin(mx,mz)
	end module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	program main
	use declare
	integer irank

	call initialize

	do 1 irank=1,10000
	call setdt
	time=time+dt
!	print*, dt, time, irank
	call stepon
    1 continue

	!$acc update host(x)
	print*,x(100,99), time, dt

	endprogram main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine initialize
	use declare
	!$acc parallel default(present)
	!$acc loop independent collapse(2)
	do 1 jz=1,mz
	do 1 jx=1,mx
	A(jx,jz)=jx*1.d0
	B(jx,jz)=jz*1.d0
	x(jx,jz)=A(jx,jz)+B(jx,jz)
	xdif(jx,jz)=0.d0
    1 continue
	!$acc end parallel
	dt=0.d0
	time=0.d0
	!$acc update host(A,B,x,xdif)
	return
	endsubroutine initialize


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine setdt
	use declare
	dt=1.d0
	!$acc parallel present(x)
	!$acc loop reduction(min:dt) independent collapse(2)
	do 1 jz=1,mz
	do 1 jx=1,mx
	dt=min(dt,1.d0/dabs(x(jx,jz)+1.d-6))
    1 continue
	!$acc end parallel

	return
	endsubroutine setdt


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine stepon
	use declare
	call right
	!$acc parallel copyin(dt) default(present)
	!$acc loop independent collapse(2)
	do 1 jz=1,mz
	do 1 jx=1,mx
	x(jx,jz)=x(jx,jz)+xdif(jx,jz)*dt
    1 continue
	!$acc end parallel
	! update host(x)
	call boundary
	! update device(x)
	return
	endsubroutine stepon


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine right
	use declare
	! update host(x)
	call boundary
	! update device(x)
	call convt
	!$acc parallel default(present)
	!$acc loop independent collapse(2)
	do 1 jz=2,mz-1
	do 1 jx=2,mx-1
	xdif(jx,jz)=(A(jx+1,jz)-A(jx-1,jz))/2.d0 &
	           +(B(jx,jz+1)-B(jx,jz-1))/2.d0
    1 continue
	!$acc end parallel
	return
	endsubroutine right

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine convt
	use declare

	!$acc parallel default(present)
	!$acc loop independent collapse(2)
	do 1 jz=2,mz-1
	do 1 jx=2,mx-1
	A(jx,jz)=x(jx,jz)
	B(jx,jz)=x(jx,jz)
    1 continue
	!$acc end parallel
	return
	endsubroutine convt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine boundary
	use declare
	! kernels default(present)
	!$acc update host(x(1:2,:),x(mx-1:mx,:),x(:,1:2),x(:,mz-1:mz))
	x(1,:)=x(2,:)
	x(mx,:)=x(mx-1,:)
	x(:,1)=x(:,2)
	x(:,mz)=x(:,mz-1)
	!$acc update device(x(1:2,:),x(mx-1:mx,:),x(:,1:2),x(:,mz-1:mz))
	! end kernels
	return
	endsubroutine boundary






