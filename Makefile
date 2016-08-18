# MFA Makefile
# A. J. Smith (ajs224@cam.ac.uk)
#
FLAGS =-I./include -O3 -ggdb #-g -fast
COMPILER = g++
TARGET= spm-c++
SRC=./source

spm-c++: spm.o kernel.o theta.o moments.o Particle.o m_in.o mfa_params.o random.o parse_args.o
	${COMPILER} -o ${TARGET} ${FLAGS} spm.o kernel.o theta.o moments.o \
	Particle.o m_in.o mfa_params.o random.o parse_args.o

spm.o: ${SRC}/spm.c++
	${COMPILER} -c ${FLAGS} ${SRC}/spm.c++

driver.o: ${SRC}/driver.c++
	${COMPILER} -c ${FLAGS} ${SRC}/driver.c++

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

random.o: ${SRC}/random.c++
	${COMPILER} -c ${FLAGS} ${SRC}/random.c++

parse_args.o: ${SRC}/parse_args.c++
	${COMPILER} -c ${FLAGS} ${SRC}/parse_args.c++

clean:
	rm ${TARGET} *.o

clean-data:
	rm data/*.txt

backup:
	tar -czvf mfa_src.tar.gz * --exclude *.gz #*.f *.h Makefile start
