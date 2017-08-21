module mu_equations
  implicit none
  real(8), parameter :: pi = 3.141592653589793d0

contains

  !************************************************************************
  !     Mu root equation - used to find mu for given tc
  subroutine mu_root_eqn(t, lambda, w_e, mu, dos_mu, n, emin, &
       dee, damp, maxiter, tol, dos, q_list, p_list, output, iter, &
       nc, nee)

    real(8), intent(in) :: t, lambda, w_e, dee, emin, mu, dos_mu, &
         n, damp, tol
    integer, intent(in) :: nc, nee, maxiter
    real(8), intent(in) :: dos(0:nee-1), q_list(nee-1), p_list(nee-1)
    real(8), intent(out) :: output
    integer, intent(out) :: iter
    real(8) :: ssum, integral
    real(8) :: zeta(nc), chi(nc)
    integer :: i
    ssum = 0.0d0
    call simult_mu_func(t, lambda, w_e, mu, dos_mu, emin, &
         dee, damp, maxiter, tol, dos, q_list, p_list, zeta, chi, iter, nc, nee)
    !$omp parallel do default(shared) private(i, integral) reduction(+:ssum)
    do i = 1, nc
       call quad_interp_c(q_list, p_list, zeta(i), chi(i), &
            dee, emin, mu, integral, nee)
       ssum = ssum + 2.d0*integral
    end do
    !$omp end parallel do
    output = n - (1.d0 - 2.d0*t*ssum)
  end subroutine mu_root_eqn
  !************************************************************************
  !     Section to simultaneously converge zeta and chi for a given
  !     mu and tc
  subroutine simult_mu_func(t, lambda, w_e, mu, dos_mu, emin, &
       dee, damp, maxiter, tol, dos, q_list, p_list, &
       new_zeta, new_chi, iter, nc, nee)

    real(8), intent(in) :: t, lambda, w_e, dee, emin, mu, &
         damp, tol, dos_mu
    integer, intent(in) :: nc, nee, maxiter
    real(8), intent(in) :: dos(nee), q_list(nee-1), p_list(nee-1)
    real(8), intent(out) :: new_zeta(nc), new_chi(nc)
    integer, intent(out) :: iter
    real(8) :: old_chi(nc), old_zeta(nc), init_zeta(nc), init_chi(nc)
    real(8) :: chi_diff, zeta_diff
    integer :: i
    logical :: converged = .false.
    iter = -999
    do i = 1,nc 
       init_chi(i) = 0.d0
       init_zeta = pi*t*(2.d0*i - 1.d0)
    end do
    call zeta(t, lambda, w_e, mu, dos_mu, dee, emin, &
         dos, q_list, p_list, init_zeta, init_chi, new_zeta, nc, nee)
    call chi(t, lambda, w_e, mu, dos_mu, dee, emin, &
         damp, dos, q_list, p_list, new_zeta, init_chi, new_chi, nc, nee)
    call f_compare(init_zeta, new_zeta, zeta_diff, nc)
    call f_compare(init_chi, new_chi, chi_diff, nc)
    if ((max(zeta_diff, chi_diff)).le.tol) then
       iter = 0
       converged = .true.
    end if
    if (.not.converged) then 
       do i = 1, maxiter
          old_chi = new_chi
          old_zeta = new_zeta
          call zeta(t, lambda, w_e, mu, dos_mu, dee, emin,&
               dos, q_list, p_list, old_zeta, old_chi, new_zeta, nc, nee)
          call chi(t, lambda, w_e, mu, dos_mu, dee, emin, &
               damp, dos, q_list, p_list, new_zeta, old_chi, new_chi, nc, nee)
          call f_compare(old_zeta, new_zeta, zeta_diff, nc)
          call f_compare(old_chi, new_chi, chi_diff, nc)
          if ((max(zeta_diff, chi_diff)).le.tol) then
             iter = i
             exit
          end if
          if (i.eq.maxiter) then
             iter = -i
          end if
       end do
    end if
  end subroutine simult_mu_func
  !************************************************************************
  !    Section to calculate chi
  subroutine chi(t, lambda, w_e, mu, dos_mu, dee, emin, &
       damp, dos, q_list, p_list, zeta, old_chi, new_chi, nc, nee)

    real(8), intent(in) :: t, lambda, w_e, dee, emin, damp, &
         mu, dos_mu
    integer, intent(in) :: nc, nee
    real(8), intent(in) :: dos(0:nee-1), q_list(nee-1), p_list(nee-1)
    real(8), intent(in) :: zeta(nc), old_chi(nc)
    real(8), intent(out) :: new_chi(nc)
    integer :: j
    real(8) :: matsu_sum, w_m
    do j = 1,nc
       w_m = pi*t*(2.d0*j-1.d0)
       call chi_sum(t, w_e, w_m, mu, dos, q_list, p_list, dee, &
            emin, zeta, old_chi, matsu_sum, nc, nee)
       new_chi(j) = (1-damp)*(-lambda/dos_mu*t*matsu_sum) &
            + damp*old_chi(j)
    end do
  end subroutine chi

  subroutine chi_sum(t, w_e, w_m, mu, dos, q_list, p_list, dee, &
       emin, zeta, old_chi, ssum, nc, nee)

    real(8), intent(in) :: t, w_e, w_m, dee, emin, mu
    integer, intent(in) :: nc, nee
    real(8), intent(in) :: zeta(nc), old_chi(nc)
    real(8), intent(in) :: dos(0:nee-1), q_list(nee-1), p_list(nee-1)       
    real(8), intent(out) :: ssum
    real(8) :: w_n, integral, w_e2, lam_even
    integer :: n
    ssum = 0.0d0
    w_e2 = w_e*w_e
    !$omp parallel do default(shared) private(n, w_n, lam_even, integral) reduction(+:ssum)
    do n = 1,nc
       integral = 0.0d0
       w_n = pi*t*(2*n-1.d0)
       lam_even = (w_e2)*(1.d0/(w_e2+(w_m-w_n)**2)&
            +1.d0/(w_e2+(w_m+w_n)**2))
       call quad_interp_c(q_list, p_list, zeta(n), old_chi(n), &
            dee, emin, mu, integral, nee)
       !       call quad_c(dos, emin, emax, mu, dee, zeta(n), old_chi(n), &
       !            w_e, integral, nee)
       ssum = ssum + (lam_even*integral)
    end do
    !$omp end parallel do
  end subroutine chi_sum

  subroutine quad_c(dos, emin, mu, dee, zeta_n, chi_n, &
       integral, nee)

    integer, intent(in) :: nee
    real(8), intent(in) :: emin, dee, zeta_n, chi_n, mu
    real(8), intent(in) :: dos(0:nee-1)
    real(8), intent(out) :: integral
    real(8) :: ebar_j, ebar_j1, ebar_j2, ebar_j12
    integer :: j
    integral = 0.0d0
    !!$omp parallel do default(shared) private(j, ebar_j, ebar_j2, ebar_j1, ebar_j12) &
    !!$omp & reduction(+:integral)
    do j = 0, nee-2
       ebar_j = (emin + dee*j - mu + chi_n)
       ebar_j2 = ebar_j*ebar_j
       ebar_j1 = (emin + dee*(j+1) - mu + chi_n)
       ebar_j12 = ebar_j1*ebar_j1
       integral = integral + 0.5d0*dee/pi*(dos(j)*ebar_j/(ebar_j2 &
            + zeta_n**2) + dos(j+1)*ebar_j1/(ebar_j12 + zeta_n**2))
    end do
    !!$omp end parallel do
  end subroutine quad_c
  !****************************************************************************
  !    Section to calculate Zeta
  subroutine zeta(t, lambda, w_e, mu, dos_mu, dee, emin, &
       dos, q_list, p_list, old_zeta, chi, new_zeta, nc, nee)

    real(8), intent(in) :: t, lambda, w_e, dee, emin, &
         mu, dos_mu
    integer, intent(in) :: nc, nee
    real(8), intent(in) :: dos(0:nee-1), q_list(nee-1), p_list(nee-1)
    real(8), intent(in) :: old_zeta(nc), chi(nc)
    real(8), intent(out) :: new_zeta(nc)
    integer :: m
    real(8) :: matsu_sum, w_m
    do m = 1,nc
       w_m = pi*t*(2.d0*m-1.d0)
       call zeta_sum(t, w_e, w_m, mu, dos, q_list, p_list, dee, &
            emin, old_zeta, chi, matsu_sum, nc, nee)
       new_zeta(m) = (w_m + lambda/dos_mu*t*matsu_sum)
    end do
  end subroutine zeta

  subroutine zeta_sum(t, w_e, w_m, mu, dos, q_list, p_list, dee, &
       emin, old_zeta, chi, ssum, nc, nee)
    real(8), intent(in) :: t, w_e, w_m, dee, emin, mu
    integer, intent(in) :: nc, nee
    real(8), intent(in) :: old_zeta(nc), chi(nc)
    real(8), intent(in) :: dos(0:nee-1), q_list(nee-1), p_list(nee-1)       
    real(8), intent(out) :: ssum
    real(8) :: w_n, integral, w_e2, lam_odd
    integer :: n
    ssum = 0.0d0
    w_e2 = w_e*w_e
    integral = 0.0d0
    !$omp parallel do default(shared) private(n, w_n, lam_odd, integral) reduction(+:ssum)
    do n = 1,nc
       w_n = pi*t*(2.d0*n-1.d0)
       lam_odd = (w_e2)*(1.d0/(w_e2+(w_m-w_n)**2) &
            -1.d0/(w_e2+(w_m+w_n)**2))
       call quad_interp_z(q_list, p_list, old_zeta(n), chi(n), &
            dee, emin, mu, integral, nee)
       ! call quad_z(dos, emin, emax, mu, dee, old_zeta(n), chi(n), &
       !      w_e, integral, nee)
       ssum = ssum + lam_odd*old_zeta(n)*integral
    end do
    !$omp end parallel do
  end subroutine zeta_sum

  subroutine quad_z(dos, emin, mu, dee, zeta_n, chi_n, &
       integral, nee)

    integer, intent(in) :: nee
    real(8), intent(in) :: emin, dee, zeta_n, chi_n, mu
    real(8), intent(in) :: dos(0:nee-1)
    real(8), intent(out) :: integral
    real(8) :: ebar_j, ebar_j1, ebar_j2, ebar_j12
    integer :: j
    integral = 0.0d0
    !$omp parallel do default(shared) private(j, ebar_j, ebar_j2, ebar_j1, ebar_j12)  reduction(+:integral)
    do j = 0, nee-2
       ebar_j = (emin + dee*j - mu + chi_n)
       ebar_j2 = ebar_j*ebar_j
       ebar_j1 = (emin + dee*(j+1) - mu + chi_n)
       ebar_j12 = ebar_j1*ebar_j1
       integral = integral + 0.5d0*dee/pi*(dos(j)/(ebar_j2 &
            + zeta_n**2) + dos(j+1)/(ebar_j12 + zeta_n**2))
    end do
    !$omp end parallel do
  end subroutine quad_z
  
  subroutine quad_interp_z(q_list, p_list, zeta_n, chi_n, &
       dee, emin, mu, integral, nee)
    integer, intent(in) :: nee
    real(8), intent(in) :: q_list(1: nee-1), p_list(1: nee-1)
    real(8), intent(in) :: zeta_n, chi_n, dee, emin, mu
    real(8), intent(out) :: integral
    integer :: i
    real(8) :: F_b, F_a
    integral = 0.0d0
    do i = 1, nee-1
       call anti_derivative_z(emin+dee*i, p_list(i), q_list(i), mu - chi_n, zeta_n, F_b)
       call anti_derivative_z(emin+dee*(i-1), p_list(i), q_list(i), mu - chi_n, zeta_n, F_a)       
       integral = integral + (F_b - F_a)
    end do
    ! call anti_derivative_z(emin+dee*(nee-1), p_list(1), q_list(1), mu - chi_n, zeta_n, F_b)
    ! call anti_derivative_z(emin+dee, p_list(1), q_list(1), mu - chi_n, zeta_n, F_a)
    ! integral = integral + (F_b - F_a)
  end subroutine quad_interp_z

  subroutine anti_derivative_z(x, p_i, q_i, a, b, val)
    ! arguments refer to the following integrand:
    ! ( p_i*x + q_i ) / ( ( x - a )^2 + b^2 )
    ! the numerator is the linear interpolation of the elect. DOS
    ! p_i = ( g_(i+1) - g_i ) / dee, where g_i is the dos evaluated at e_i  
    ! q_i =  e_i * ( 1 - p_i )
    ! a = (mu - chi_n) for integrals in zeta and phi
    real(8) :: val, x, p_i, q_i, a, b
    val = 1.d0*0.5d0*p_i*dlog( (x-a) * (x-a) + b * b ) &
         - (a*p_i + q_i)*datan((a-x)/b)/b 
  end subroutine anti_derivative_z

  subroutine quad_interp_c(q_list, p_list, zeta_n, chi_n,&
       dee, emin, mu, integral, nee)
    integer, intent(in) :: nee
    real(8), intent(in) :: q_list(1: nee-1), p_list(1: nee-1)
    real(8), intent(in) :: zeta_n, chi_n, dee, emin, mu
    real(8), intent(out) :: integral
    integer :: i
    real(8) :: F_b, F_a
    integral = 0.0d0
    do i = 1, nee-1
       call anti_derivative_c(emin+dee*i, p_list(i), q_list(i), mu - chi_n, zeta_n, F_b)
       call anti_derivative_c(emin+dee*(i-1), p_list(i), q_list(i), mu - chi_n, zeta_n, F_a)       
       integral = integral + (F_b - F_a)
    end do  
    ! call anti_derivative_c(emin+dee*(nee-1), p_list(1), q_list(1), mu - chi_n, zeta_n, F_b)
    ! call anti_derivative_c(emin+dee*(1-1), p_list(1), q_list(1), mu - chi_n, zeta_n, F_a)
    ! integral = integral + (F_b - F_a)
  end subroutine quad_interp_c

  subroutine anti_derivative_c(x, p_i, q_i, a, b, val)
    ! arguments refer to the following integrand:
    ! ( p_i*x + q_i ) * ( x - a ) / ( ( x - a )^2 + b^2 )
    ! p_i = ( g_(i+1) - g_i ) / dee, where g_i is the dos evaluated at e_i  
    ! q_i =  e_i * ( 1 - p_i )
    ! a = (mu - chi_n) for integrals in zeta and phi
    real(8) :: val, x, p_i, q_i, a, b
    val = 0.5d0*(a*p_i + q_i)*dlog( (x-a)*(x-a) + b*b ) &
         + b*p_i*datan( ( a - x ) / b ) + p_i*( x - a )
    val = val
  end subroutine anti_derivative_c

  subroutine anti_derivative_num_z(x, p, q, a, b, val)
    ! arguments refer to the following integrand:
    ! ( p*x + q ) * ( x - a )^2 / ( ( x - a )^2 + b^2 )
    ! the numerator is the linear interpolation of the elect. DOS
    ! p = ( g_(i+1) - g_i ) / dee, where g_i is the dos evaluated at e_i  
    ! q =  e_i * ( 1 - p_i )
    ! a = (mu - chi_n) for integrals in zeta and phi
    ! b = zeta_n
    real(8) :: val, x, p, q, a, b
    val = 0.5d0*( b*b - p*dlog( (x-a)*(x-a) + b*b ) &
         + 2*b*( a*p + q )*datan( ( a-x ) / b ) + ( x - a )*( a*p + p*x + 2*q ) ) 
  end subroutine anti_derivative_num_z


  !************************************************************************
  !    Function comparator
  subroutine f_compare(v1, v2, result, length)

    integer, intent(in) :: length
    real(8), intent(in) :: v1(0:length-1), v2(0:length-1)
    real(8), intent(out) :: result
    integer :: i
    result = 0.0d0
    do i = 0, length-1
       result = result + abs((v1(i)/v2(i))**2-1.d0)
       ! result = result + ((v1(i) - v2(i))**2)
    end do
    result = result/(length*1.0d0)
  end subroutine f_compare


end module mu_equations
