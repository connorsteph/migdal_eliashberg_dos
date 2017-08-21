      program denbccsum
      use prmts
      implicit none
      real(dp) eemin,eemax,fact,del2,pih
      real(dp) summp,summm,summt,ee,eeb,den3,den3p,den3m,wt,prod,arg
      real(dp) argp,argm,dx
      integer nee,i,ix,iy,iz,ie
      real(dp),dimension(1:4512) :: xx,coss,cos5
      character (len=30) :: in_file = 'bcc_config.init'
      call init(in_file)
      open(unit=9,file=out_file)
      write(*,*)'ttp,ttx,del,nx = ',ttp,ttx,del,nx
      write(9,*)'ttp,ttx,del,nx '
      write(9,*) ttp,ttx,del,nx
      dx = 1.0d0/dfloat(nx) !delta x in integral
      pih=pi/2
c     populates domain array xx, and the cos(xx), as well as cos(xx/2)
c     xx is positive
      do i = 1,nx
        xx(i) = (i - 0.5d0)*dx ! using midpoint trapezoidal ! xx goes from dx/2 to (1000-1/2)dx
        coss(i) = dcos(pi*xx(i))
        cos5(i) = dcos(pih*xx(i))
      end do
c     here the parameters of the problem are set, as well as the sampling of the region in energy
      eemin = -8*ttp - 6*ttx     !lower bound for energy  ! everything in terms of the bcc hop
      eemax = 8*ttp - 6*ttx !upper bound
      nee = nint((eemax - eemin)/dee) + 1 !number of points in energy
      del2 = del*del
      fact = 0.125d0*dx*dx*dx/(ttp*dsqrt(pi*del2)) !coefficient outside of integral
      summp = 0.0d0
      summm = 0.0d0
      do ie = 1,nee             !energy values go from emin to emax, in steps of dee
         ee = eemin + (ie - 1)*dee 
         eeb = 0.125d0*ee/ttp   !ebar (e/(8t))
         den3p = 0.0d0
         den3m = 0.0d0   
         !Begin triple integral using midpoint trapezoidal, along a wedge-like region in space
         do ix = 1,nx           !integrate along all points in x
            do iy = 1,ix        !integrate along  y, up to ix
               do iz = 1,iy     !integrate along z, up to iy
                  wt = 3.0d0    !default weighting of points
                  prod = cos5(ix)*cos5(iy)*cos5(iz)
                  arg = 0.25d0*ttx*(coss(ix) + coss(iy) + coss(iz))/ttp + eeb
                  argp = arg + prod
                  argm = arg - prod
c     assigning weighting to points depending on index
                  if ((ix.ne.iy).and.(iy.ne.iz).and.(ix.ne.iz)) then
                     wt = 6.0d0
                  else if ((ix.eq.iy).and.(iy.eq.iz)) then
                     wt = 1.0d0
                  endif
                  den3p = den3p + wt*fact*dexp(-argp*argp/del2) !delta function with plus arg
                  den3m = den3m + wt*fact*dexp(-argm*argm/del2) !delta function with minus arg
               end do
            end do
         end do
         if (ie.ne.1) then      !checks if we are beyond the first loop in energy vals, if so then add
            summp = summp + den3p 
            summm = summm + den3m
         else
            summp = summp + 0.5d0*den3p !otherwise, weighting of 1/2 for for loop
            summm = summm + 0.5d0*den3m
         endif
         den3 = 0.5d0*(den3p + den3m)
c     write(6,99)ee,den3p,den3m,den3
c     write(9,99)ee,den3p,den3m,den3
         write(9,*) ee,den3
         if (verbose.eq.1) then 
            write(*,*) ee,den3
         end if
      end do
      summp = summp*dee
      summm = summm*dee
      summt = summp + summm
c      write(9,*)'summp, summm,summt = ',summp,summm,summt
      write(*,*)'#may30 summp, summm,summt = ',summp,summm,summt
      write(*,*) 'nx,ttx,ttp,dee,del',nx,ttx,ttp,dee,del
 99   format(1x,f10.5,1x,2(f13.7,1x,2(f13.7,1x),3x))
      stop
      end
