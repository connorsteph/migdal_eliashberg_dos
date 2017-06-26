      subroutine n(t, g, w_e, dee, emin, emax, mu,
     -      dos, zeta, chi, n_occ, nc, nee)
      implicit none
      real*8, intent(in) :: t, g, w_e, dee, emin, emax, mu
      integer, intent(in) :: nc, nee
      real*8, intent(in) :: dos(0:nee-1)
      real*8, intent(in) :: zeta(0:nc-1), chi(0:nc-1)
      real*8, intent(out) :: n_occ
      real*8, parameter :: pi = 3.1415926535897
      real*8 :: ssum, w_m, integral
      integer :: j
      external :: quad
      ssum = 0.0d0
      do j = 0, nc-1
         call quad(dos, emin, emax, mu, dee, zeta(j), chi(j), w_e,
     -        integral, nee)
         ssum = ssum + integral
      end do
      n_occ = 1 - 2*pi*t*ssum
      end subroutine

      
      subroutine quad(dos, emin, emax, mu, dee, zeta_n, chi_n,
     -     w_e, integral, nee)
      implicit none
      integer, intent(in) :: nee
      real*8, intent(in) :: emin, emax, dee, zeta_n, w_e, chi_n, mu
      real*8, intent(in) :: dos(0:nee-1)
      real*8, intent(out) :: integral
      real*8 :: ebar_j2, ebar_j12
      integer :: j
      real*8, parameter :: pi = 3.141592653589793
      integral = 0.0d0
      do j = 0, nee-2
         ebar_j2 = (emin + dee*j - mu + chi_n)*(emin + dee*j - mu
     -        + chi_n)
         ebar_j12 = (emin + dee*(j+1) - mu
     -        + chi_n)*(emin + dee*(j+1) - mu + chi_n)
         integral = integral + 0.5d0*dee/pi*(dos(j)*ebar_j2/(ebar_j2
     -    + zeta_n**2) + dos(j+1)*ebar_j12/(ebar_j12 + zeta_n**2))
      end do
      end subroutine
