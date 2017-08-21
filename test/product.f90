program product
  implicit none
  integer, parameter :: n = 100000
  real(8) :: a(0:n-1), b(0:n-1)
  real(8) :: result
  integer :: i
  do i = 0, n-1
     a(i) = i*1.0d0
     b(i) = ((n-1)-i)*1.0d0
  end do
  do i = 0,n-1
     result = result + a(i)*b(i)
  end do
end program
