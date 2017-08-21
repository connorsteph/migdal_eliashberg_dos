module phi_mat
  implicit none
  real(8), parameter :: pi = 4.d0*datan(1.d0)
contains 
  subroutine phi_matrix(t, lambda, w_e, mu, dos_mu, dee, emin, &
       zeta, chi, dos, q_list, p_list, matrix, nc, nee)
    real(8), intent(in) :: t, lambda, w_e, dee, emin, mu, dos_mu
    integer, intent(in) :: nc, nee
    real(8), intent(in) :: dos(0:nee-1), q_list(nee-1), p_list(nee-1)
    real(8), intent(in) :: zeta(nc), chi(nc)
    real(8), intent(out) :: matrix(Nc, Nc)
    integer :: i, j
    real(8) :: w_m(nc), A0(nc), dos_integral, w_e2
    real(8) :: lam_even
    w_e2 = w_e*w_e
    !$omp parallel default(shared) private(j, i, lam_even, dos_integral)
    !$omp do
    do j = 1,nc
       call quad_interp(q_list, p_list, zeta(j), chi(j), &
            dee, emin, mu, dos_integral, nee)
       !call quad(dos, emin, mu, dee, zeta(j), chi(j), dos_integral, nee)
       A0(j) = dos_integral
       w_m(j) = pi*t*(2.d0*dble(j) - 1.d0)
    end do
    !$omp end do

    !$omp do 
    do j = 1, nc
       do i = 1, nc
          matrix(i,j) = 0.0d0
          lam_even = w_e2 * ( 1.d0/( w_e2 + ( w_m(i) - w_m(j) )**2 ) &
               + 1.d0/( w_e2 + ( w_m(i) + w_m(j) )**2 ) )
          !matrix(i, j) = t*lambda/dos_mu*lam_even*A0(i)*A0(j)
          !matrix(j, i) = matrix(i, j)
          matrix(i, j) = t*lambda/dos_mu*lam_even*A0(i)
          ! if (i.eq.j) then
          !    matrix(i,j) = matrix(i,j) - A0(j)
          ! endif
       end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine phi_matrix

  subroutine quad(dos, emin, mu, dee, zeta_n, chi_n, &
       integral, nee)
    integer, intent(in) :: nee
    real(8), intent(in) :: emin, dee, zeta_n, chi_n, mu
    real(8), intent(in) :: dos(0:nee-1)
    real(8), intent(out) :: integral
    real(8) :: ebar_j, ebar_j1, ebar_j2, ebar_j12
    integer :: j
    integral = 0.0d0
    do j = 0, nee-2
       ebar_j = (emin + dee*dble(j) - mu + chi_n)
       ebar_j2 = ebar_j*ebar_j
       ebar_j1 = (emin + dee*(dble(j)+1.d0) - mu + chi_n)
       ebar_j12 = ebar_j1*ebar_j1
       integral = integral + 0.5d0*dee/pi*(dos(j)/(ebar_j2 &
            + zeta_n**2) + dos(j+1)/(ebar_j12 + zeta_n**2)) 
    end do
  end subroutine quad


  subroutine quad_interp(q_list, p_list, zeta_n, chi_n, &
       dee, emin, mu, integral, nee)
    integer, intent(in) :: nee
    real(8), intent(in) :: q_list(nee-1), p_list(nee-1)
    real(8), intent(in) :: zeta_n, chi_n, dee, emin, mu
    real(8), intent(out) :: integral
    integer :: i
    real(8) :: F_b, F_a
    integral = 0.0d0
!!$omp parallel do default(shared) private(i, F_b, F_a) reduction(+:integral) 
    do i = 1, nee-1
       call anti_derivative(emin+dee*i, p_list(i), q_list(i), mu - chi_n, zeta_n, F_b)
       call anti_derivative(emin+dee*(i-1), p_list(i), q_list(i), mu - chi_n, zeta_n, F_a)
       integral = integral + (F_b - F_a)
    end do
!!$omp end parallel do

   ! call anti_derivative(emin+dee*(nee-1), p_list(1), q_list(1), mu - chi_n, zeta_n, F_b)
   ! call anti_derivative(emin+dee*(1-1), p_list(1), q_list(1), mu - chi_n, zeta_n, F_a)
   ! integral = integral + (F_b - F_a)
  end subroutine quad_interp

  subroutine anti_derivative(x, p_i, q_i, a, b, val)
  ! arguments refer to the following integrand:
  ! ( p_i*x + q_i ) / ( ( x - a )^2 + b^2 )
  ! the numerator is the linear interpolation of the elect. DOS
  ! p_i = ( g_(i+1) - g_i ) / dee, where g_i is the dos evaluated at e_i  
  ! q_i =  e_i * ( 1 - p_i )
  ! a = (mu - chi_n) for integrals in zeta and phi
    real(8) :: val, x, p_i, q_i, a, b
    val = 1.d0/pi*0.5d0*p_i*dlog( (x-a) * (x-a) + b*b) &
         - (a*p_i + q_i)*datan((a-x)/abs(b))/b/pi 
  end subroutine anti_derivative

end module phi_mat
