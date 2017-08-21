program omp_product
  implicit none
  integer, parameter :: n = 10000
  real(8) :: a(0:n-1), b(0:n-1)
  real(8) :: result
  integer :: i
  integer, parameter :: chunk = 2500
  do i = 0, n-1
     a(i) = i*1.0d0
     b(i) = ((n-1)-i)*1.0d0
  end do
  !$omp parallel do &
  !$omp& default(shared) private(i) &
  !$omp& reduction(+:result)
  do i = 0,n-1
     result = result + a(i)*b(i)
  end do
  !$omp end parallel do
end program omp_product
