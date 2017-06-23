      subroutine phi_matrix(t, g, w_e, dee, emin, emax,
     -     dos_mu, zeta, dos, matrix, nc, nee)
      implicit none
      real*8, intent(in) :: t, g, w_e, dee, emin, emax, dos_mu
      integer, intent(in) :: nc, nee
      real*8, intent(in) :: dos(0:nee-1)
      real*8, intent(in) :: zeta(0:nc-1)
      real*8, intent(out) :: matrix(0 : nc-1, 0 : nc-1)
      integer :: i, j
      real*8, parameter :: pi = 3.1415926535897
      real*8 :: lambda, w_i, w_j, dos_integral, w_e2
      external :: quad
      lambda = 2*dos_mu*g**2/w_e
      w_e2 = w_e*w_e
      do j = 0,nc-1
         w_j = pi*t*(2*j - 1)
         do i = 0,nc-1
            w_i = pi*t*(2*i - 1)
            call quad(dos, emin, emax, dee, zeta(j),
     -        w_e, dos_integral, nee)
            matrix(i, j) = pi*t*lambda/dos_mu*w_e2/(w_e2 +
     -           (w_i - w_j)**2)*dos_integral
         end do
      end do
      end subroutine


      subroutine quad(dos, emin, emax, dee, zeta_n,
     -     w_e, integral, nee)
      implicit none
      integer, intent(in) :: nee
      real*8, intent(in) :: emin, emax, dee, zeta_n, w_e
      real*8, intent(in) :: dos(0:nee-1)
      real*8, intent(out) :: integral
      integer :: j
      real*8, parameter :: pi = 3.141592653589793
      integral = 0.0d0
      do j = 0, nee-2
         integral = integral + 0.5d0*dee/pi*(dos(j)/((emin
     -        +dee*(j))**2+zeta_n**2)+dos(j+1)/((emin
     -        +dee*(j+1))**2+zeta_n**2))
      end do
      end subroutine
