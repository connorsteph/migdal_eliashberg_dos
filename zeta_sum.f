      subroutine zeta(t, g, w_e, mu, dee, emin, emax, dos_mu,
     -      damp, dos, old_zeta, chi, new_zeta, nc, nee)
      implicit none
      real*8, intent(in) :: t, g, w_e, dee, emin, emax, damp, mu, dos_mu
      integer, intent(in) :: nc, nee
      real*8, intent(in) :: dos(0:nee-1)
      real*8, intent(in) :: old_zeta(0:nc-1), chi(0:nc-1)
      real*8, intent(out) :: new_zeta(0:nc-1)
      integer :: j
      real*8, parameter :: pi = 3.1415926535897
      real*8 :: matsu_sum, lambda, w_m
      external :: zeta_sum
      lambda = 2*dos_mu*g**2/w_e
      do j = 0,nc-1
         w_m = pi*t*(2*(j+1)-1)
         call zeta_sum(t, g, w_e, w_m, mu, dos, dee,
     -        emin, emax, old_zeta, chi, matsu_sum, nc, nee)
         new_zeta(j) = (1-damp)*(w_m + lambda/dos_mu*t*pi*matsu_sum)
     -        + damp*old_zeta(j)
      end do
      end subroutine
      
      subroutine zeta_init(t, g, w_e, mu, dee, emin,
     -     emax,  dos_mu, dos, chi, new_zeta, nc, nee)
      implicit none
      real*8, intent(in) :: t, g, w_e, dee, emin, emax, mu, dos_mu
      integer, intent(in) :: nc, nee
      real*8, intent(in) :: dos(0:nee-1)
      real*8, intent(in) ::  chi(0:nc-1)
      real*8, intent(out) :: new_zeta(0:nc-1)
      integer :: j
      real*8, parameter :: pi = 3.1415926535897
      real*8 :: matsu_sum, lambda, w_m
      external :: zeta_sum_init
      lambda = 2*dos_mu*g**2/w_e
      do j = 0,nc-1
         w_m = pi*t*(2*(j+1)-1)
         call zeta_sum_init(t, g, w_e, w_m, mu, dos, dee,
     -        emin, emax, chi, matsu_sum, nc, nee)
         new_zeta(j) = w_m + lambda/dos_mu*t*pi*matsu_sum
      end do
      end subroutine
      
      subroutine zeta_sum(t, g, w_e, w_m, mu, dos, dee,
     - emin, emax, old_zeta, chi, ssum, nc, nee)
      implicit none
      real*8, intent(in) :: t, g, w_e, w_m, dee, emin, emax, mu
      integer, intent(in) :: nc, nee
      real*8, intent(in) :: old_zeta(0:nc-1), chi(0:nc-1)
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
         call quad(dos, emin, emax, dee, mu, old_zeta(n), chi(n),
     -        w_e, integral, nee)
         ssum = ssum + lam_odd*old_zeta(n)*integral
      end do
      end subroutine

      subroutine quad(dos, emin, emax, mu, dee, zeta_n, chi_n,
     -     w_e, integral, nee)
      implicit none
      integer, intent(in) :: nee
      real*8, intent(in) :: emin, emax, dee, zeta_n, chi_n, w_e, mu
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

      subroutine zeta_sum_init(t, g, w_e, w_m, mu, dos, dee,
     - emin, emax, chi, ssum, nc, nee)
      implicit none
      real*8, intent(in) :: t, g, w_e, w_m, dee, emin, emax, mu
      integer, intent(in) :: nc, nee
      real*8, intent(in) :: dos(0:nee-1)
      real*8, intent(in) :: chi(0:nc-1)
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
         call quad_init(dos, emin, emax, mu, dee, w_n,
     -        w_e, chi(n), integral, nee)
         ssum = ssum + lam_odd*w_n*integral
      end do
      end subroutine

      subroutine quad_init(dos, emin, emax, mu, dee, zeta_n,
     -     w_e, chi_n, integral, nee)
      implicit none
      integer, intent(in) :: nee
      real*8, intent(in) :: emin, emax, dee, zeta_n, chi_n, w_e, mu
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

