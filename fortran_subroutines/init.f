      subroutine init(in_file)
      use prmts
      implicit none
      character (len = 30) :: in_file
      character c*11,str*54
      write(*,*) 'Reading file: ', in_file
      write(*,*) '**********************************************'
      open(unit=9,file=in_file,form='formatted')
      read(9,10) c,ttp
      write(*,*) c,ttp
      read(9,10) c,ttx
      write(*,*) c,ttx
      read(9,10) c,dee
      write(*,*) c,dee
      read(9,10) c,del
      write(*,*) c,del
      read(9,20) c,nx
      write(*,*) c,nx
      read(9,20) c,verbose
      write(*,*) c,verbose
      read(9,30) c,out_file
      write(*,*) c,out_file
      write(*,*) '**********************************************'

      close(9)
      
 10   format(A11,F8.3) !floating point
 20   format(A11,I10) !integer
 30   format(A11,A40) !character
      
      end subroutine  init
