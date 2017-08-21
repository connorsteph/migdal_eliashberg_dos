program test
  implicit none
  integer tid, omp_get_thread_num
  !$omp parallel private(tid)
  tid = omp_get_thread_num()
  print *, 'hello world from thread = ', tid
  !$omp end parallel
end program test

