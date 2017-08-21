opts = gfortran -o
optsc = gfortran -c -O4 -ffixed-line-length-128

all: o1 main

prmts.o: prmts.f
	$(optsc) prmts.f

init.o: init.f
	$(optsc) init.f

denbccsum.o: denbccsum.f
	$(optsc) denbccsum.f

o1: prmts.o init.o denbccsum.o

main: o1
	$(opts) main prmts.o init.o denbccsum.o
clean:
	rm -f *.o
	rm -f *.mod
	rm -f main
