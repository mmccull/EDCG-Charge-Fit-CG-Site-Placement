
prefix = /usr/local/bin
mpif90 = /opt/local/lib/openmpi/bin/mpif90
f90 = gfortran
#flags = -ftree-vectorize -O3 -fopenmp -ffree-line-length-none
flags = -O3 -ffree-line-length-none -fopenmp
lapacklib = /usr/local/lib

fit : edcg_charge_site_placement.f90 dcd.f90 pca.f90
	$(f90) -c edcg_charge_site_placement.f90 dcd.f90 pca.f90  $(flags) -L$(lapacklib) -llapack -lblas
	$(f90)  edcg_charge_site_placement.o dcd.o pca.o -o edcg_charge_site_placement.x  $(flags) -L$(lapacklib) -llapack -lblas
	cp edcg_charge_site_placement.x ../test/actin_test

clean :
	rm -f *.o *.x *.mod

