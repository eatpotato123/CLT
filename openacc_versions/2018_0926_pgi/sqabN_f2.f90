program main
	use tseqm
	implicit none
	integer, parameter :: N=100
	integer :: idx
	real x(N)

	!$acc parallel loop copyout(x)
	do idx=1,N
	x(idx)=sqab(idx*1.0)
	enddo
	print*, "x(3)=", x(3)
endprogram	
