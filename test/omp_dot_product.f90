subroutine omp_dot_product(a, b, prod, n)
  implicit none
  integer, intent(in) :: n
  real(8), intent(in) :: a(n), b(n)
  real(8), intent(out) :: prod
  integer :: i, chunk
  integer, parameter :: chunksize = 10
  prod = 0.0d0
  chunk = chunksize
  !$omp parallel do
  !$omp& default(shared) private(i)
  !$omp& schedule(static, chunk)
  !$omp& reduction(+:prod)
  do i = 1,n
     prod = prod + (a(i)*b(i))
  end do
  !$omp end parallel do
end subroutine omp_dot_product
