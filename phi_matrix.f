      subroutine phi_matrix(t, g, w_e, mu, dee, emin, emax,
     -     dos_mu, zeta, chi, dos, matrix, nc, nee)
      implicit none
      real*8, intent(in) :: t, g, w_e, dee, emin, emax, dos_mu, mu
      integer, intent(in) :: nc, nee
      real*8, intent(in) :: dos(0:nee-1)
      real*8, intent(in) :: zeta(0:nc-1), chi(0:nc-1)
      real*8, intent(out) :: matrix(0 : nc-1, 0 : nc-1)
      integer :: i, j
      real*8, parameter :: pi = 3.1415926535897
      real*8 :: lambda, w_i, w_j, dos_integral, w_e2
      real*8 :: lam_even
      external :: quad
      lambda = 2*dos_mu*g**2/w_e
      w_e2 = w_e*w_e
      do j = 0,nc-1
         w_j = pi*t*(2*(j+1) - 1)
         do i = 0,nc-1
            w_i = pi*t*(2*(i+1) - 1)
            lam_even = (w_e2)*(1/(w_e2+(w_i-w_j)**2)
     -        +1/(w_e2+(w_i+w_j)**2))
            call quad(dos, emin, emax, mu, dee, zeta(j), chi(j),
     -        w_e, dos_integral, nee)
            matrix(i, j) = pi*t*lambda/dos_mu
     - *lam_even*dos_integral
         end do
      end do
      end subroutine


      subroutine quad(dos, emin, emax, mu, dee, zeta_n, chi_n,
     -     w_e, integral, nee)
      implicit none
      integer, intent(in) :: nee
      real*8, intent(in) :: emin, emax, dee, zeta_n, w_e, chi_n, mu
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
