ext_suffix := $(shell python3-config --extension-suffix)

all: ./phi_matrix$(ext_suffix) ./omp_phi_matrix$(ext_suffix) ./omp_mu_zeta_chi$(ext_suffix) ./mu_zeta_chi$(ext_suffix)

./phi_matrix$(ext_suffix): ./fortran_subroutines/phi_matrix.f90
	f2py3 -c -m phi_matrix ./fortran_subroutines/phi_matrix.f90

./omp_phi_matrix$(ext_suffix): ./fortran_subroutines/phi_matrix.f90
	f2py3 -c -m --f90flags='-fopenmp' omp_phi_matrix ./fortran_subroutines/phi_matrix.f90 -lgomp

./mu_zeta_chi$(ext_suffix): ./fortran_subroutines/mu_zeta_chi.f90
	f2py3 -c -m mu_zeta_chi ./fortran_subroutines/mu_zeta_chi.f90

./omp_mu_zeta_chi$(ext_suffix): ./fortran_subroutines/mu_zeta_chi.f90
	f2py3 -c -m --f90flags='-fopenmp' omp_mu_zeta_chi ./fortran_subroutines/mu_zeta_chi.f90 -lgomp

.PHONY: clean

clean: 
	rm -f *$(ext_suffix)