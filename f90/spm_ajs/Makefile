FLAGS = #-fast
COMPILER = ifort
TARGET= spm_main

spm_main: spm_main.o global_data.mod compute_rates.o particle_event.o replacement.o collision.o rand240.o
	${COMPILER} -o ${TARGET} ${FLAGS} spm_main.o global_data.o compute_rates.o particle_event.o replacement.o collision.o rand240.o

spm_main.o: spm_main.f90 global_data.mod
	${COMPILER} -c ${FLAGS} spm_main.f90

global_data.mod: global_data.f90 rand240.o
	${COMPILER} -c ${FLAGS} global_data.f90

compute_rates.o: compute_rates.f90
	${COMPILER} -c ${FLAGS} compute_rates.f90

particle_event.o: particle_event.f90
	${COMPILER} -c ${FLAGS} particle_event.f90

replacement.o: replacement.f90
	${COMPILER} -c ${FLAGS} replacement.f90 

collision.o: collision.f90
	${COMPILER} -c ${FLAGS} collision.f90 

rand240.o: rand240.f90
	${COMPILER} -c ${FLAGS} rand240.f90

clean:
	rm *.o *.mod

backup:
	tar -czvf spm_src.tar.gz * --exclude *.gz #*.f *.h Makefile start
