
	module declare
	integer, parameter :: mx=1280, mz=1280
	integer jx,jz
	real*8, dimension(mx,mz) :: x, A, B, xdif
	real*8 time
	real*8 dt
	!$acc declare create(x,A,B,xdif) copyin(mx,mz)
	end module
