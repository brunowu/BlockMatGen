ALL: blib exec 

#compilation and various flags
DIRS    = 
EXEC    = generateur
CFLAGS	= 
FFLAGS	= 
CPPFLAGS	= 
FPPFLAGS	=
CLEANFILES  = ${EXEC}
OFILES= ${wildcard ./*.o}

#Execution flags
MAT = utm300.mtx_300x300_3155nnz1
MPI_NODES = 2
NB = 2
MATTYPE = matblock

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

blib :
	-@echo "BEGINNING TO COMPILE LIBRARIES "
	-@echo "========================================="
	-@${OMAKE}  PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} ACTION=libfast tree
	-@echo "Completed building libraries"
	-@echo "========================================="

distclean :
	-@echo "Cleaning application and libraries"
	-@echo "========================================="
	-@${OMAKE} PETSC_ARCH=${PETSC_ARCH}  PETSC_DIR=${PETSC_DIR} clean
	-${RM} ${OFILES}
	-@echo "Finised cleaning application and libraries"
	-@echo "========================================="	

exec: gen.o
	-@echo "BEGINNING TO COMPILE APPLICATION "
	-@echo "========================================="
	-@${CLINKER} -o ${EXEC} gen.o ${PETSC_LIB}
	-@echo "Completed building application"
	-@echo "========================================="

runa:
	-@${MPIEXEC} -np ${MPI_NODES} ./generateur -mfile ${MAT} -nb ${NB} -type ${MATTYPE}

