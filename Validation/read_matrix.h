/*CREATED BY Xinzhe WU 2017*/
#ifndef _READ_MATRIX_H
#define _READ_MATRIX_H

#include "petscmat.h"
#include "petscerror.h"
#include <stdlib.h>
#include <stdio.h>

#endif
PetscErrorCode read_matrix(Mat * A, char *filea);
PetscErrorCode compMat(Mat A, Mat B);
