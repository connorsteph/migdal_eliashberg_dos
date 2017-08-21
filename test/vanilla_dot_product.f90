subroutine vanilla_dot_product(a, b, prod, n)
  implicit none
  integer, intent(in) :: n
  real(8), intent(in) :: a(n), b(n)
  real(8), intent(out) :: prod
  integer :: i
  prod = 0.0d0
  do i = 1,n
     prod = prod + a(i)*b(i)
  end do
end subroutine
