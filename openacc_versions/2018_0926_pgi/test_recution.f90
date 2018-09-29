module declare
	integer, parameter:: N=100
	integer i, j, a(N,N)
	real*8 ired
end module

subroutine main
	use declare
	implicit none

	do i=1,N
	do j=1,N
	a(i,j)=i
	enddo
	enddo
	ired = 0
	!$acc parallel default(present)
	!$acc loop reduction(min:ired) independent collapse(2)
	do 1 i=1,N
	do 1 j=1,N
	ired = min(ired,a(i,j)*1.d0)
    1 continue
	!$acc end parallel
	print*, "ired =", ired
end
