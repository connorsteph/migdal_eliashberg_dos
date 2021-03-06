      subroutine phi(t, g, w_e, dee, emin, emax, dos_mu, p_damp,
     -     new_phi, old_phi, zeta, dos, nc, nee)
      implicit none
      real*8, intent(in) :: t, g, w_e, dee, emin, emax, p_damp, dos_mu
      integer, intent(in) :: nc, nee
      real*8, intent(in) :: dos(0:nee-1)
      real*8, intent(in) :: old_phi(0:nc-1), zeta(0:nc-1)
      real*8, intent(out) :: new_phi(0:nc-1)
      integer :: j
      real*8, parameter :: pi = 3.1415926535897
      real*8 :: matsu_sum, lambda, w_m
      external :: phi_sum
      lambda = 2*dos_mu*g**2/w_e
      do j = 0,nc-1
         w_m = pi*t*(2*(j+1)-1)
         call phi_sum(t, g, w_e, w_m, dee,
     - emin, emax, matsu_sum, old_phi, zeta, dos, nc, nee)
         new_phi(j) = (1-p_damp)*(lambda/dos_mu*t*pi*matsu_sum)
     -        + p_damp*old_phi(j)
      end do
      do j = 0,nc-1
         new_phi(j) = new_phi(j)/new_phi(0)
      end do
      end subroutine
      
      subroutine phi_init(t, g, w_e, dee, emin, emax, dos_mu,
     -     init_phi, zeta, dos, new_phi, nc, nee)
      implicit none
      real*8, intent(in) :: t, g, w_e, dee, emin, emax, dos_mu
      integer, intent(in) :: nc, nee
      real*8, intent(in) :: dos(0:nee-1)
      real*8, intent(in) :: init_phi(0:nc-1), zeta(0:nc-1)
      real*8, intent(out) :: new_phi(0:nc-1)
      integer :: j
      real*8, parameter :: pi = 3.1415926535897
      real*8 :: matsu_sum, lambda, w_m
      external :: phi_sum_init
      lambda = 2*dos_mu*g**2/w_e
      do j = 0,nc-1
         w_m = pi*t*(2*(j+1)-1)
         call phi_sum_init(t, g, w_e, w_m, dee,
     - emin, emax, matsu_sum, init_phi, zeta, dos, nc, nee)
         new_phi(j) = lambda/dos_mu*t*pi*matsu_sum
      end do
      do j = 0,nc-1
         new_phi(j) = new_phi(j)/new_phi(0)
      end do
      end subroutine
      
      subroutine phi_sum(t, g, w_e, w_m, dee,
     - emin, emax, matsu_sum, old_phi, zeta, dos, nc, nee)
      implicit none
      real*8, intent(in) :: t, g, w_e, w_m, dee, emin, emax
      integer, intent(in) :: nc, nee
      real*8, intent(in) :: old_phi(0:nc-1), zeta(0:nc-1)
      real*8, intent(in) :: dos(0:nee-1)       
      real*8, intent(out) :: matsu_sum
      real*8 :: w_n, integral, w_e2, lam_even
      integer :: n
      external :: quad
      real*8, parameter :: pi = 3.1415926535897
      matsu_sum = 0.0d0
      w_e2 = w_e*w_e
      do n = 0,nc-1
         integral = 0.0d0
         w_n = pi*t*(2*(n+1)-1)
         lam_even = (w_e2)*(1/(w_e2+(w_m-w_n)**2)
     -        +1/(w_e2+(w_m+w_n)**2))
         call quad(dos, emin, emax, dee, zeta(n),
     -        w_e, integral, nee)
         matsu_sum = matsu_sum + lam_even*old_phi(n)*integral
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

      subroutine phi_sum_init(t, g, w_e, w_m, dee,
     - emin, emax, matsu_sum, init_phi, zeta, dos, nc, nee)
      implicit none
      real*8, intent(in) :: t, g, w_e, w_m, dee, emin, emax
      integer, intent(in) :: nc, nee
      real*8, intent(in) :: dos(0:nee-1)
      real*8, intent(in) :: init_phi(0:nc-1), zeta(0:nc-1)
      real*8, intent(out) :: matsu_sum
      real*8 :: w_n, integral, w_e2, lam_even
      integer :: n
      external :: quad_init
      real*8, parameter :: pi = 3.1415926535897
      matsu_sum = 0.0d0
      w_e2 = w_e*w_e
      do n = 0,nc-1
         integral = 0.0d0
         w_n = pi*t*(2*(n+1)-1)
         lam_even = (w_e2)*(1/(w_e2+(w_m-w_n)**2)
     -        +1/(w_e2+(w_m+w_n)**2))
         call quad_init(dos, emin, emax, dee, zeta(n),
     -        w_e, integral, nee)
         matsu_sum = matsu_sum + lam_even*init_phi(n)*integral
      end do
      end subroutine

      subroutine quad_init(dos, emin, emax, dee, zeta_n,
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
