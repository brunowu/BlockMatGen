#ifndef _GEN_H
#define _GEN_H

#include "petscmat.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct _MatrixInfo{
	int n;
	int m;
	int nnz;
} MatrixInfo;

PetscReal uniform_distribution(PetscReal rangeLow, PetscReal rangeHigh);
PetscErrorCode MatGenbyOther();
PetscErrorCode read_matrix_vector(Mat * A);

#endif
