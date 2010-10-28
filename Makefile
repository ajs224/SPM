# SPM Makefile
# A. J. Smith (ajs224@cam.ac.uk)
#
FLAGS =-I./include -O3 #-g -fast
COMPILER = g++
TARGET= spm-c++
SRC=./source

spm-c++: spm.o kernel.o theta.o moments.o Particle.o m_in.o mfa_params.o
	${COMPILER} -o ${TARGET} ${FLAGS} mfa.o kernel.o theta.o moments.o Particle.o m_in.o mfa_params.o

spm.o: ${SRC}/spm.c++
	${COMPILER} -c ${FLAGS} ${SRC}/spm.c++

kernel.o: ${SRC}/kernel.c++
	${COMPILER} -c ${FLAGS} ${SRC}/kernel.c++

moments.o: ${SRC}/moments.c++
	${COMPILER} -c ${FLAGS} ${SRC}/moments.c++

theta.o: ${SRC}/theta.c++
	${COMPILER} -c ${FLAGS} ${SRC}/theta.c++

m_in.o: ${SRC}/m_in.c++
	${COMPILER} -c ${FLAGS} ${SRC}/m_in.c++


Particle.o: ${SRC}/Particle.c++
	${COMPILER} -c ${FLAGS} ${SRC}/Particle.c++


mfa_params.o: ${SRC}/mfa_params.c++
	${COMPILER} -c ${FLAGS} ${SRC}/mfa_params.c++

clean:
	rm ${TARGET} *.o

clean-data:
	rm data/*.txt

backup:
	tar -czvf spm_src.tar.gz * --exclude *.gz #*.f *.h Makefile start
