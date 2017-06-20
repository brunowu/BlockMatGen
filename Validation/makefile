ALL: blib exec 

#compilation and various flags
DIRS    = 
EXEC    = compMat
CFLAGS	= 
FFLAGS	= 
CPPFLAGS	= 
FPPFLAGS	=
CLEANFILES  = ${EXEC}
OFILES= ${wildcard ./*.o}


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

exec: read_matrix.o
	-@echo "BEGINNING TO COMPILE APPLICATION "
	-@echo "========================================="
	-@${CLINKER} -o ${EXEC} read_matrix.o ${PETSC_LIB}
	-@echo "Completed building application"
	-@echo "========================================="

#-ksp_monitor_true_residual -eps_monitor



