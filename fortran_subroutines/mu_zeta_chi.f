c************************************************************************
c     Mu root equation - used to find mu for given tc
      subroutine mu_root_eqn(t, lambda, w_e, dee, emin, emax, mu,
     -     dos_mu, n, maxiter, damp, tol, dos, init_zeta, init_chi,
     -     result, nc, nee)
      implicit none
      real*8, intent(in) :: t, lambda, w_e, dee, emin, emax, mu, dos_mu,
     -     n, damp, tol
      integer, intent(in) :: nc, nee, maxiter
      real*8, intent(in) :: dos(0:nee-1)
      real*8, intent(in) :: init_zeta(0:nc-1), init_chi(0:nc-1)
      real*8, intent(out) :: result
      real*8, parameter :: pi = 3.1415926535897
      real*8 :: ssum, integral
      real*8 :: zeta(0:nc-1), chi(0:nc-1)
      integer :: i, iter
      external :: quad_c, simult_mu_func
      ssum = 0.0d0
      call simult_mu_func(t, lambda, w_e, mu, dos_mu, emin, emax,
     -     dee, damp, maxiter, tol, dos, init_zeta, init_chi, zeta,
     -     chi, iter, nc, nee)
      do i = 0, nc-1
         call quad_c(dos, emin, emax, mu, dee, zeta(i), chi(i),
     -     w_e, integral, nee)
         ssum = ssum + 2*integral
      end do
      result = n - (1 - 2*pi*t*ssum)
      end subroutine
c************************************************************************
c     Section to simultaneously converge zeta and chi for a given
c     mu and tc
      subroutine simult_mu_func(t, lambda, w_e, mu, dos_mu, emin, emax,
     -     dee, damp, maxiter, tol, dos, init_zeta, init_chi,
     -     new_zeta, new_chi, iter, nc, nee)
      implicit none
      real*8, intent(in) :: t, lambda, w_e, dee, emin, emax, mu,
     -     damp, tol, dos_mu
      integer, intent(in) :: nc, nee, maxiter
      real*8, intent(in) :: dos(0:nee-1)
      real*8, intent(in) :: init_zeta(0:nc-1), init_chi(0:nc-1)
      real*8, intent(out) :: new_zeta(0:nc-1), new_chi(0:nc-1)
      integer, intent(out) :: iter
      real*8 :: old_chi(0:nc-1), old_zeta(0:nc-1)
      real*8, parameter :: pi = 3.1415926535897
      real*8 :: chi_diff, zeta_diff
      integer :: i
      external :: chi, zeta, f_compare
      iter = 0
      call zeta(t, lambda, w_e, mu, dos_mu, dee, emin, emax,
     -     0.3, dos, init_zeta, init_chi, new_zeta, nc, nee)
      call chi(t, lambda, w_e, mu, dos_mu, dee, emin, emax,
     -      damp, dos, new_zeta, init_chi, new_chi, nc, nee)
      call f_compare(init_zeta, new_zeta, zeta_diff, nc)
      call f_compare(init_chi, new_chi, chi_diff, nc)
      if ((max(zeta_diff, chi_diff)).le.tol) then
         goto 20
      end if
      do i = 1, maxiter
         old_chi = new_chi
         old_zeta = new_zeta
         call zeta(t, lambda, w_e, mu, dos_mu, dee, emin, emax,
     -        0.3, dos, old_zeta, old_chi, new_zeta, nc, nee)
         call chi(t, lambda, w_e, mu, dos_mu, dee, emin, emax,
     -        damp, dos, new_zeta, old_chi, new_chi, nc, nee)
         call f_compare(old_zeta, new_zeta, zeta_diff, nc)
         call f_compare(old_chi, new_chi, chi_diff, nc)
         if ((max(zeta_diff, chi_diff)).le.tol) then
            iter = i
            goto 20 
         end if
         if (i.eq.maxiter) then
            iter = -i
         end if
      end do
 20   continue
      end subroutine 
c************************************************************************
c     Section to calculate chi
      subroutine chi(t, lambda, w_e, mu, dos_mu, dee, emin, emax,
     -      damp, dos, zeta, old_chi, new_chi, nc, nee)
      implicit none
      real*8, intent(in) :: t, lambda, w_e, dee, emin, emax, damp,
     -     mu, dos_mu
      integer, intent(in) :: nc, nee
      real*8, intent(in) :: dos(0:nee-1)
      real*8, intent(in) :: zeta(0:nc-1), old_chi(0:nc-1)
      real*8, intent(out) :: new_chi(0:nc-1)
      integer :: j
      real*8, parameter :: pi = 3.1415926535897
      real*8 :: matsu_sum, w_m
      external :: chi_sum
      do j = 0,nc-1
         w_m = pi*t*(2*(j+1)-1)
         call chi_sum(t, w_e, w_m, mu, dos, dee,
     -        emin, emax, zeta, old_chi, matsu_sum, nc, nee)
         new_chi(j) = (1-damp)*(-lambda/dos_mu*t*pi*matsu_sum)
     -        + damp*old_chi(j)
      end do
      end subroutine

      subroutine chi_sum(t, w_e, w_m, mu, dos, dee,
     - emin, emax, zeta, old_chi, ssum, nc, nee)
      implicit none
      real*8, intent(in) :: t, w_e, w_m, dee, emin, emax, mu
      integer, intent(in) :: nc, nee
      real*8, intent(in) :: zeta(0:nc-1), old_chi(0:nc-1)
      real*8, intent(in) :: dos(0:nee-1)       
      real*8, intent(out) :: ssum
      real*8 :: w_n, integral, w_e2, lam_even
      integer :: n
      external :: quad
      real*8, parameter :: pi = 3.1415926535897
      ssum = 0.0d0
      w_e2 = w_e*w_e
      do n = 0,nc-1
         integral = 0.0d0
         w_n = pi*t*(2*(n+1)-1)
         lam_even = (w_e2)*(1/(w_e2+(w_m-w_n)**2)
     -        +1/(w_e2+(w_m+w_n)**2))
         call quad_c(dos, emin, emax, mu, dee, zeta(n), old_chi(n),
     -        w_e, integral, nee)
         ssum = ssum + lam_even*integral
      end do
      end subroutine

      subroutine quad_c(dos, emin, emax, mu, dee, zeta_n, chi_n,
     -     w_e, integral, nee)
      implicit none
      integer, intent(in) :: nee
      real*8, intent(in) :: emin, emax, dee, zeta_n, chi_n, w_e, mu
      real*8, intent(in) :: dos(0:nee-1)
      real*8, intent(out) :: integral
      real*8 :: ebar_j, ebar_j1, ebar_j2, ebar_j12
      integer :: j
      real*8, parameter :: pi = 3.1415926535897
      integral = 0.0d0
      do j = 0, nee-2
         ebar_j = (emin + dee*j - mu + chi_n)
         ebar_j2 = ebar_j*ebar_j
         ebar_j1 = (emin + dee*(j+1) - mu + chi_n)
         ebar_j12 = ebar_j1*ebar_j1
         integral = integral + 0.5d0*dee/pi*(dos(j)*ebar_j/(ebar_j2
     -    + zeta_n**2) + dos(j+1)*ebar_j1/(ebar_j12 + zeta_n**2))
      end do
      end subroutine
c****************************************************************************
c     Section to calculate Zeta
      subroutine zeta(t, lambda, w_e, mu, dos_mu, dee, emin, emax,
     -      damp, dos, old_zeta, chi, new_zeta, nc, nee)
      implicit none
      real*8, intent(in) :: t, lambda, w_e, dee, emin, emax, damp,
     -     mu, dos_mu
      integer, intent(in) :: nc, nee
      real*8, intent(in) :: dos(0:nee-1)
      real*8, intent(in) :: old_zeta(0:nc-1), chi(0:nc-1)
      real*8, intent(out) :: new_zeta(0:nc-1)
      integer :: j
      real*8, parameter :: pi = 3.1415926535897
      real*8 :: matsu_sum, w_m
      external :: zeta_sum
      do j = 0,nc-1
         w_m = pi*t*(2*(j+1)-1)
         call zeta_sum(t, w_e, w_m, mu, dos, dee,
     -        emin, emax, old_zeta, chi, matsu_sum, nc, nee)
         new_zeta(j) = (1-damp)*(w_m + lambda/dos_mu*t*pi*matsu_sum)
     -        + damp*old_zeta(j)
      end do
      end subroutine

      subroutine zeta_sum(t, w_e, w_m, mu, dos, dee,
     - emin, emax, old_zeta, chi, ssum, nc, nee)
      implicit none
      real*8, intent(in) :: t, w_e, w_m, dee, emin, emax, mu
      integer, intent(in) :: nc, nee
      real*8, intent(in) :: old_zeta(0:nc-1), chi(0:nc-1)
      real*8, intent(in) :: dos(0:nee-1)       
      real*8, intent(out) :: ssum
      real*8 :: w_n, integral, w_e2, lam_odd
      integer :: n
      external :: quad_z
      real*8, parameter :: pi = 3.1415926535897
      ssum = 0.0d0
      w_e2 = w_e*w_e
      do n = 0,nc-1
         integral = 0.0d0
         w_n = pi*t*(2*(n+1)-1)
         lam_odd = (w_e2)*(1/(w_e2+(w_m-w_n)**2)
     -        -1/(w_e2+(w_m+w_n)**2))
         call quad_z(dos, emin, emax, mu, dee, old_zeta(n), chi(n),
     -        w_e, integral, nee)
         ssum = ssum + lam_odd*old_zeta(n)*integral
      end do
      end subroutine
      
      subroutine quad_z(dos, emin, emax, mu, dee, zeta_n, chi_n,
     -     w_e, integral, nee)
      implicit none
      integer, intent(in) :: nee
      real*8, intent(in) :: emin, emax, dee, zeta_n, chi_n, w_e, mu
      real*8, intent(in) :: dos(0:nee-1)
      real*8, intent(out) :: integral
      real*8 :: ebar_j, ebar_j1, ebar_j2, ebar_j12
      integer :: j
      real*8, parameter :: pi = 3.1415926535897
      integral = 0.0d0
      do j = 0, nee-2
         ebar_j = (emin + dee*j - mu + chi_n)
         ebar_j2 = ebar_j*ebar_j
         ebar_j1 = (emin + dee*(j+1) - mu + chi_n)
         ebar_j12 = ebar_j1*ebar_j1
         integral = integral + 0.5d0*dee/pi*(dos(j)/(ebar_j2
     -    + zeta_n**2) + dos(j+1)/(ebar_j12 + zeta_n**2))
      end do
      end subroutine
c************************************************************************
c     Function comparator
      subroutine f_compare(v1, v2, result, length)
      implicit none
      integer, intent(in) :: length
      real*8, intent(in) :: v1(0:length-1), v2(0:length-1)
      real*8, intent(out) :: result
      integer :: i
      result = 0.0d0
      do i = 0, length-1
         result = result + abs((v1(i)/v2(i)) - 1)/(length*1.0d0)
      end do
      end subroutine
