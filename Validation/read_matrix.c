/*CREATED BY Xinzhe WU 2017*/
#include "read_matrix.h"

static char help[] = "Matrix Market File Validation";

PetscErrorCode read_matrix(Mat * A, char *filea){
  //  char filea[PETSC_MAX_PATH_LEN];
  PetscErrorCode ierr;
  PetscViewer fd;
  PetscInt size,sizea;

  /* read matrix file */
  PetscPrintf(PETSC_COMM_WORLD,"Loading Matrix : %s\n",filea);
  ierr=PetscViewerBinaryOpen(PETSC_COMM_WORLD,filea,FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD,A);CHKERRQ(ierr);
  ierr=MatLoad(*A,fd);CHKERRQ(ierr);
  ierr=PetscViewerDestroy(&fd);CHKERRQ(ierr);
  ierr=MatGetSize(*A,&size,&sizea);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"Loaded Matrix of size : %d %d\n",size,sizea);

  return 0;
}

PetscErrorCode compMat(Mat A, Mat B){
  
  PetscErrorCode ierr;
  PetscInt ma,na,mb,nb;
  PetscBool flag;
  Vec Id,res;
  Mat minus;

  ierr = MatGetSize(A, &ma,&na);
  ierr = MatGetSize(B,&mb, &nb);
/*
  VecCreate(PETSC_COMM_WORLD, &Id);
  VecSetSizes(Id, PETSC_DECIDE, ma);
  VecSet(Id, 1.0);
  VecDuplicate(Id, &res);

  MatCreate(PETSC_COMM_WORLD, &minus);
  MatSetSizes(minus, PETSC_DECIDE, PETSC_DECIDE, ma, na);
*/
  if( (ma == mb) && (na == nb)){

    ierr = MatEqual(A, B, &flag);
    
    if(flag)
      {
	PetscPrintf(PETSC_COMM_WORLD, ">>>>The two matrix are equal\n\n");
      }
    else
      {
	PetscPrintf(PETSC_COMM_WORLD, ">>>>The two matrix are NOT equal\n\n");
	ierr = MatAXPY(B, -1.0, A, DIFFERENT_NONZERO_PATTERN);
	ierr = MatView(B, PETSC_VIEWER_STDOUT_WORLD);
      }
  }
  else{
    PetscPrintf(PETSC_COMM_WORLD, "****Error: Not the same dimension matrix\n\n");
    return 0;
  }
  
  
  return 0;
}

int main(int argc, char ** argv){

  PetscInt sizen;
  Mat A,B;
  PetscErrorCode ierr;

  PetscInitialize(&argc,&argv,(char *)0,help);
  MPI_Comm_size(PETSC_COMM_WORLD,&sizen);
  if (sizen > 1) SETERRQ(PETSC_COMM_WORLD,1,"Uniprocessor only\n");

  ierr = read_matrix(&A,"Block_matrix_nb_1_300x300_3.21615e-314_nnz");
  ierr = read_matrix(&B, "utm300.mtx_300x300_3155nnz1");

  //ierr = MatView(A, PETSC_VIEWER_STDOUT_WORLD);	
  
  ierr = compMat(A, B);
  
  PetscFinalize();

  return 0;
}
