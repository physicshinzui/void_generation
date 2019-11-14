FC = gfortran
#FC = ifort
#OPT = -fpp -g -CB  -traceback -fpe0  
#OPT = -fpp -g -CB  -traceback -fpe0 -D DEBUG 
OPT = -cpp -O2 -static 
LOAD = pdb_void_disp
SRC =	Enum_Type.f90		\
		ASA.f                   \
		pdb_void_disp.f 

ALL	:	${SRC}
	${FC} ${OPT} -o ${LOAD} ${SRC} $(LIBS)

clean	:
	rm -f *.o *.mod ${LOAD}
