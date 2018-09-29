	include 'demo_part1.f90'
!	include 'demo_part3.f90'
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

