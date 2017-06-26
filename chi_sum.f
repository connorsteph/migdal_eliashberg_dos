      subroutine chi(t, g, w_e, dee, emin, emax, mu, dos_mu,
     -      damp, dos, zeta, new_chi, old_chi, nc, nee)
      implicit none
      real*8, intent(in) :: t, g, w_e, dee, emin, emax, damp, mu, dos_mu
      integer, intent(in) :: nc, nee
      real*8, intent(in) :: dos(0:nee-1)
      real*8, intent(in) :: zeta(0:nc-1), old_chi(0:nc-1)
      real*8, intent(out) :: new_chi(0:nc-1)
      integer :: j
      real*8, parameter :: pi = 3.1415926535897
      real*8 :: matsu_sum, lambda, w_m
      external :: chi_sum
      lambda = 2*dos_mu*g**2/w_e
      do j = 0,nc-1
         w_m = pi*t*(2*(j+1)-1)
         call chi_sum(t, g, w_e, w_m, dos, dee,
     -        emin, emax, mu, old_zeta, old_chi, matsu_sum, nc, nee)
         new_chi(j) = (1-damp)*(-lambda/dos_mu*t*pi*matsu_sum)
     -        + damp*old_chi(j)
      end do
      end subroutine
      
      subroutine chi_init(t, g, w_e, dee, emin,
     -     emax,  mu, dos_mu, dos, zeta, new_chi, nc, nee)
      implicit none
      real*8, intent(in) :: t, g, w_e, dee, emin, emax, mu, dos_mu
      integer, intent(in) :: nc, nee
      real*8, intent(in) :: dos(0:nee-1)
      real*8, intent(in) :: zeta(0:nc-1)
      real*8, intent(out) :: new_chi(0:nc-1)
      integer :: j
      real*8, parameter :: pi = 3.1415926535897
      real*8 :: matsu_sum, lambda, w_m
      external :: chi_sum_init
      lambda = 2*dos_mu*g**2/w_e
      do j = 0,nc-1
         w_m = pi*t*(2*(j+1)-1)
         call chi_sum_init(t, g, w_e, w_m, mu, dee,
     -        emin, emax, dos, zeta, init_chi, matsu_sum, nc, nee)
         new_chi(j) = -lambda/dos_mu*t*pi*matsu_sum
      end do
      end subroutine
      
      subroutine chi_sum(t, g, w_e, w_m, dos, dee,
     - emin, emax, mu, zeta, old_chi, ssum, nc, nee)
      implicit none
      real*8, intent(in) :: t, g, w_e, w_m, dee, emin, emax
      integer, intent(in) :: nc, nee
      real*8, intent(in) :: zeta(0:nc-1), old_chi(0:nc-1)
      real*8, intent(in) :: dos(0:nee-1)       
      real*8, intent(out) :: ssum
      real*8 :: w_n, integral, w_e2, lam_odd
      integer :: n
      external :: quad
      real*8, parameter :: pi = 3.1415926535897
      ssum = 0.0d0
      w_e2 = w_e*w_e
      do n = 0,nc-1
         integral = 0.0d0
         w_n = pi*t*(2*(n+1)-1)
         lam_odd = (w_e2)*(1/(w_e2+(w_m-w_n)**2)
     -        -1/(w_e2+(w_m+w_n)**2))
         call quad(dos, emin, emax, dee, zeta(n), old_chi(n),
     -        w_e, integral, nee)
         ssum = ssum + lam_odd*integral
      end do
      end subroutine

      subroutine quad(dos, emin, emax, dee, zeta_n, chi_n,
     -     w_e, integral, nee)
      implicit none
      integer, intent(in) :: nee
      real*8, intent(in) :: emin, emax, dee, zeta_n, chi_n, w_e
      real*8, intent(in) :: dos(0:nee-1)
      real*8, intent(out) :: integral
      real*8 :: ebar_j, ebar_j_1
      integer :: j
      real*8, parameter :: pi = 3.1415926535897
      integral = 0.0d0
      do j = 0, nee-2
         ebar_j = emin + dee*j - mu
         ebar_j_1 = emin + dee*(j+1) - mu
         integral = integral + 0.5d0*dee/pi*(dos(j)/((ebar_j + chi_n)**2
     -        + zeta_n**2) + dos(j+1)/((ebar_j_1 + chi_n)**2
     -        + zeta_n**2))
      end do
      end subroutine

      subroutine chi_sum_init(t, g, w_e, w_m, mu, dee,
     - emin, emax, dos, zeta, init_chi, ssum, nc, nee)
      implicit none
      real*8, intent(in) :: t, g, w_e, w_m, dee, emin, emax
      integer, intent(in) :: nc, nee
      real*8, intent(in) :: dos(0:nee-1)       
      real*8, intent(out) :: ssum
      real*8 :: w_n, integral, w_e2, lam_odd
      integer :: n
      external :: quad_init
      real*8, parameter :: pi = 3.141592653589793
      ssum = 0.0d0
      w_e2 = w_e*w_e
      do n = 0,nc-1
         integral = 0.0d0
         w_n = pi*t*(2*(n+1)-1)
         lam_odd = (w_e2)*(1/(w_e2+(w_m-w_n)**2)
     -        -1/(w_e2+(w_m+w_n)**2))
         call quad_init(dos, emin, emax, mu, dee, zeta_n,
     -        init_chi_n, w_e, integral, nee)
         ssum = ssum + lam_odd*integral
      end do
      end subroutine

      subroutine quad_init(dos, emin, emax, mu, dee, zeta_n, chi_n,
     -     w_e, integral, nee)
      implicit none
      integer, intent(in) :: nee
      real*8, intent(in) :: emin, emax, dee, zeta_n, chi_n, w_e
      real*8, intent(in) :: dos(0:nee-1)
      real*8, intent(out) :: integral
      real*8 :: ebar_j, ebar_j_1
      integer :: j
      real*8, parameter :: pi = 3.141592653589793
      integral = 0.0d0
      do j = 0, nee-2
         ebar_j = emin + dee*j - mu
         ebar_j_1 = emin + dee*(j+1) - mu
         integral = integral + 0.5d0*dee/pi*(dos(j)/((ebar_j + chi_n)**2
     -        + zeta_n**2) + dos(j+1)/((ebar_j_1 + chi_n)**2
     -        + zeta_n**2))
      end do
      end subroutine
