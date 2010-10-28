FLAGS = #-fast
COMPILER = gfortran
TARGET= spm_new_main

spm_new_main: spm_new_main.o global_data.mod rand240.o kernel.o theta.o #collision.o
	${COMPILER} -o ${TARGET} ${FLAGS} spm_new_main.o global_data.o rand240.o kernel.o theta.o #collision.o

spm_new_main.o: nr_mods global_data.mod spm_new_main.f90 
	${COMPILER} -c ${FLAGS} spm_new_main.f90

global_data.mod: global_data.f90 rand240.o
	${COMPILER} -c ${FLAGS} global_data.f90

#compute_rates.o: compute_rates.f90
#	${COMPILER} -c ${FLAGS} compute_rates.f90

#particle_event.o: particle_event.f90
#	${COMPILER} -c ${FLAGS} particle_event.f90

#replacement.o: replacement.f90
#	${COMPILER} -c ${FLAGS} replacement.f90 

#collision.o: collision.f90
#	${COMPILER} -c ${FLAGS} collision.f90 

rand240.o: rand240.f90
	${COMPILER} -c ${FLAGS} rand240.f90

kernel.o: kernel.f90
	${COMPILER} -c ${FLAGS} kernel.f90

theta.o: theta.f90
	${COMPILER} -c ${FLAGS} theta.f90

nr_mods: NR_Modules.f90
	${COMPILER} -c ${FLAGS} NR_Modules.f90

clean:
	rm *.o *.mod

backup:
	tar -czvf spm_src.tar.gz * --exclude *.gz #*.f *.h Makefile start
