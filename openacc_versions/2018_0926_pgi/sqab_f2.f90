module tseqm
	contains
	real function sqab(a)
		  !$acc routine seq
		  real,intent(in):: a
		  sqab=sqrt(abs(a))
	end function
end module
